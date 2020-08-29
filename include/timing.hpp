#ifndef _TIMING_HPP
#define _TIMING_HPP

#ifdef WITH_TIMINGS

#include <chrono>
#include <iostream>

#define TIMER_BEGIN(id)						\
    auto __tc_beg##id = std::chrono::system_clock::now()	\

#define TIMER_END(id)							\
    do {								\
	auto __tc_end##id = std::chrono::system_clock::now();		\
	std::chrono::duration<double> __tc_d {__tc_end##id - __tc_beg##id}; \
	std::cerr << #id << ": " << __tc_d.count() << " s.\n";		\
    } while(0)

#else

#define TIMER_BEGIN(id)
#define TIMER_END(id)

#endif

#endif // _TIMING_HPP
