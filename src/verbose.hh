#ifndef VERBOSE_HH
#define VERBOSE_HH

#include <sys/resource.h> // getrusage
#include "stdlib.h"
#include "string.h"
#include <iostream>
#include <chrono>
#include <thread>
#include <atomic>

namespace verb {

    long long int mem_usage_kb_(){ // in kB
        FILE* file = fopen("/proc/self/status", "r");
        if (file) {
            long long int result = -1;
            char line[128];

            while (fgets(line, 128, file) != NULL){
                if (strncmp(line, "VmRSS:", 6) == 0){
                    // parse line :
                    int i = strlen(line);
                    const char* p = line;
                    while (*p <'0' || *p > '9') p++;
                    line[i-3] = '\0'; // assumes line ends in " kB"
                    result = atoll(p);
                    break;
                }
            }
            fclose(file);
            return result;
        } else {
            struct rusage usage;
            getrusage(RUSAGE_SELF, &usage);
            return usage.ru_maxrss / 1000;
        }
    }

    double time_() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    }

    std::atomic<double> t_now(0.);
    std::atomic<long long int> mem_now(0);
    double t_init = 0., t_wall = 0., t_last = 0., t_lap = 0., t_prog = 1.;
    int verbosity = 0;
    
    std::thread update_now;
    
    void begin(int verb = 1, int wall_sec = 7 * 24 * 3600) {
        verbosity = verb;
        t_init = time_();
        t_now.store(t_init, std::memory_order_release);
        mem_now.store(mem_usage_kb_(), std::memory_order_release);
        t_lap = t_init;
        t_last = t_init - 1; 
        t_wall = t_init + wall_sec;
        auto upd = []() {
            std::chrono::milliseconds ms10(10);
            while (t_now.load(std::memory_order_relaxed) < t_wall) {
                std::this_thread::sleep_for(ms10);
                t_now.store(time_(), std::memory_order_release);
                mem_now.store(mem_usage_kb_(), std::memory_order_release);
            }
        };
        update_now = std::thread(upd);
        std::cerr << "-- start\n";
    }

    double time_now_s() { return t_now.load(std::memory_order_acquire); }
    long long int mem_now_kb() {
        return mem_now.load(std::memory_order_acquire);
    }

    double top(double t1, std::string msg) {
        double t2 = time_();
        t_now.store(t2, std::memory_order_release);
        if (verbosity >= 1) {
            std::cerr << "-- " << msg << " " << (t2 - t1) << "s";
            long long int mem = mem_usage_kb_();
            std::cerr << " " << mem / 1000 << "m\n";
            mem_now.store(mem, std::memory_order_release);
            std::cerr.flush();
        }
        return t2;
    }

    void lap(std::string msg) {
        t_lap = top(t_lap, msg);
        t_prog = 1.;
    }

    double lap_time() {
        return time_() - t_lap;
    }

    bool delay(int millis=1000) {
        if (verbosity >= 1 && (time_now_s() - t_last) * 1000. >= millis) {
            t_last = time_now_s();
            return true;
        }
        return false;
    }

    bool progress(float fact = 1.0) {
        if (verbosity >= 1 && time_now_s() >= t_last + fact*t_prog) {
            t_last = time_now_s();
            t_prog *= 1.2;
            return true;
        }
        return false;
    }

    std::ostream nowhere(NULL); // NULL buffer -> bad bit set
    
    std::ostream & cerr(std::string msg, int verb = 2) {
        if (verb > verbosity) return nowhere;
        std::cerr << "  "
                  << (time_now_s() - t_init) <<"s "
                  << mem_now_kb() / 1000 << "m "
                  << msg << (msg.size() > 0 ? " " : "") <<": ";
        return std::cerr;
    }

    std::ostream & cerr(int verb = 2) {
        return cerr("", verb);
    }

    void end() {
        t_wall = t_init;
        top(t_init, "end");
        update_now.join();
    }

}


#endif // VERBOSE_HH
