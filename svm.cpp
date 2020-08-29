#include "shamir.hpp"

#include "ncomm.hpp"
#include "ztk.hpp"
#include "crypto.hpp"
#include "timing.hpp"

#include <iostream>
#include <fstream>
#include <string>

using namespace ncomm;
using namespace std;
using namespace goat;

const size_t K = 64;
const size_t D = 4;

// number of features in the model
const size_t nfeatures = 3073;
// number of images we perform a prediction on
const size_t nimages = 100;

typedef GR<K, D> share_t;

vector<string> split(string line) {
    vector<string> r;
    stringstream s (line);
    string v;
    while (getline(s, v, ' ')) {
	r.push_back(v);
    }
    return r;
}

vector<vector<share_t>> load_model(ShamirProtocol<K,D>& prot) {

    // read model
    string line;
    ifstream f ("svm/model_2_16.bin");
    if (!f.is_open())
	throw -1;

    auto n = prot.number_of_parties();
    auto t = prot.threshold();

    vector<Z2k<K>> vals;

    while (getline(f, line)) {
	auto inputs = split(line);
	for (auto i: inputs)
	    vals.emplace_back(Z2k<K>{stoul(i)});
    }

    cout << "model size:" << vals.size() << "\n";

    vector<vector<share_t>> x (n);

    // encode the model twice in each value we share
    for (auto &v: vals) {
	GR<K,D> model {{v, v, Z2k<K>::zero, Z2k<K>::zero}};
	auto shares = prot.share(model, t);

	for (size_t i = 0; i < n; i++) {
	    x[i].emplace_back(shares[i]);
	}
    }

    cout << "done generating shares for model\n";

    return x;
}

vector<vector<share_t>> load_images(ShamirProtocol<K,D>& prot) {

    auto n = prot.number_of_parties();
    auto t = prot.threshold();
    (void)n; (void)t;

    string line;
    ifstream f ("svm/images_2_8.bin");
    if (!f.is_open())
	throw -1;

    vector<Z2k<K>> vals;

    while(getline(f, line)) {
	auto inputs = split(line);
	for (auto i: inputs)
	    vals.emplace_back(Z2k<K>{stoul(i)});
    }

    cout << "images size:" << vals.size() << "\n";

    vector<vector<share_t>> x (n);

    auto nimages_ = vals.size() / nfeatures;

    assert(nimages == nimages_);

    cout << "number of features: " << nfeatures << "\n"
	 << "number of images: " << nimages << "\n";

    // each GR element encodes two images
    for (size_t i = 0; i < vals.size(); i += 2) {
	auto i0 = vals[i];
	auto i1 = vals[i+1];
	GR<K,D> img {{i0, i1, Z2k<K>::zero, Z2k<K>::zero}};
	auto shares = prot.share(img, t);
	for (size_t j = 0; j < n; j++) {
	    x[j].emplace_back(shares[j]);
	}
    }

    cout << "done generating shares for images\n";
    cout << "each party gets " << x[0].size() << " shares\n";
    return x;
}

int main(int argc, char** argv) {

    if (argc < 3) {
	cout << "usage: " << argv[0] << " <party ID> <network info>\n";
	return -1;
    }

    partyid_t id = stoul(argv[1]);
    Network nw (id, argv[2]);
    nw.connect();
    auto prot = ShamirProtocol<K,D>(nw);

    const size_t n = prot.number_of_parties();
    const size_t t = prot.threshold();

    vector<share_t> model (nfeatures * 10);
    vector<share_t> images ((nimages / 2) * nfeatures);

    if (id == 0) {
	auto model_shares = load_model(prot);
	assert(model.size() == model_shares[0].size());

	model = model_shares[id];

	// send model
	vector<vector<unsigned char>> sbuf (n);
	for (size_t i = 1; i < n; i++) {
	    const auto m = model_shares[i].size();
	    sbuf[i].resize(m * GR<K,D>::byte_size());
	    for (size_t j = 0; j < m; j++)
		model_shares[i][j].pack(sbuf[i].data() + (j * GR<K,D>::byte_size()));
	}

	for (size_t i = 1; i < n; i++)
	    nw.send_to(i, sbuf[i]);

	// receive images
	vector<unsigned char> imagesbytes (images.size() * GR<K,D>::byte_size());
	nw.recv_from((partyid_t)1, imagesbytes);
	for (size_t i = 0; i < images.size(); i++)
	    images[i] = GR<K,D>{imagesbytes.data() + (i*GR<K,D>::byte_size())};

    } else if (id == 1) {

	auto images_shares = load_images(prot);
	assert(images.size() == images_shares[0].size());

	images = images_shares[id];

	// obtain model
	vector<unsigned char> modelbytes (nfeatures * 10 * GR<K,D>::byte_size());
	nw.recv_from((partyid_t)0, modelbytes);

	for (size_t i = 0; i < model.size(); i++)
	    model[i] = GR<K,D>{modelbytes.data() + (i*GR<K,D>::byte_size())};

	// send images
	vector<vector<unsigned char>> sbuf(n);
	for (size_t i = 0; i < n; i++) {
	    if (i == id)
		continue;

	    const auto m = images_shares[i].size();
	    sbuf[i].resize(m * GR<K,D>::byte_size());
	    for (size_t j = 0; j < m; j++)
		images_shares[i][j].pack(sbuf[i].data() + (j * GR<K,D>::byte_size()));
	}

	for (size_t i = 0; i < n; i++) {
	    if (i == id)
		continue;
	    nw.send_to(i, sbuf[i]);
	}

    } else {

	// receive model
	vector<unsigned char> modelbytes (nfeatures * 10 * GR<K,D>::byte_size());
	nw.recv_from((partyid_t)0, modelbytes);
	for (size_t i = 0; i < model.size(); i++)
	    model[i] = GR<K,D>{modelbytes.data() + (i*GR<K,D>::byte_size())};

	// receive images
	vector<unsigned char> imagesbytes (images.size() * GR<K,D>::byte_size());
	nw.recv_from((partyid_t)1, imagesbytes);
	for (size_t i = 0; i < images.size(); i++)
	    images[i] = GR<K,D>{imagesbytes.data() + (i*GR<K,D>::byte_size())};
    }

    cout << "input done\n";

    // reshape model into the proper form
    vector<vector<share_t>> mdl (10);
    (void)t;
    for (size_t i = 0; i < 10; i++) {
	mdl[i].resize(nfeatures);
	for (size_t j = 0; j < nfeatures; j++)
	    mdl[i][j] = model[i*nfeatures + j];
    }

    // ditto for the images
    vector<vector<share_t>> img (nfeatures);
    for (size_t i = 0; i < nfeatures; i++) {
	img[i].resize(50);
	for (size_t j = 0; j < 50; j++)
	    img[i][j] = images[i*50 + j];
    }

    auto result = prot.matmul(mdl, img);
    cout << "result: " << result.size() << " x " << result[0].size() << "\n";

    // convert to bit shares over F_2^d.



    // auto m = prot.open(result[0], 2*t);

    // if (id == 0) {
    // 	for (size_t i = 0; i < 10; i++) {
    // 	    cout << m[i][0] << " " << m[i][2] << "\n";
    // 	}
    // }


}
