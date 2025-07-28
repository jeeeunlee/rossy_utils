#pragma once

#include <chrono>


// Uncomment the following line to enable the Timer functionality
#define ROSSY_TIME

#ifdef ROSSY_TIME
class Clock {
    public:
        [[nodiscard]] Clock(){start();}
        ~Clock(){}

        void start() { ini_time_ = std::chrono::high_resolution_clock::now(); }
        // return in milliseconds
        double stop() {
            end_time_ = std::chrono::high_resolution_clock::now();
            duration_ = std::chrono::duration_cast<std::chrono::microseconds>(end_time_-ini_time_);
            return double(duration_.count());
        }
        void printElapsedSec(const std::string_view message){
            std::cout << message << stop()*1e-6 << " sec " << std::endl;
            start();
        }

        void printElapsedMiliSec(const std::string_view message){
            std::cout << message << stop()*1e-3 << " msec " << std::endl;
            start();
        }

        void printElapsedMicroSec(const std::string_view message){
            std::cout << message << stop() << " µsec " << std::endl;
            start();
        }



    private:
        std::chrono::microseconds duration_;
        std::chrono::high_resolution_clock::time_point ini_time_, end_time_;
};
#else
// Empty Clock class when ENABLE_TIMER is not defined
class Clock {
    public:
        [[nodiscard]] Clock(){}
        ~Clock(){}
        void start() {}
        // return in milliseconds
        double stop() const { return 0.0; }
        void printElapsedSec(const std::string_view message){}
        void printElapsedMiliSec(const std::string_view message){}
        void printElapsedMicroSec(const std::string_view message){}

};
#endif  // ENABLE_TIMER

