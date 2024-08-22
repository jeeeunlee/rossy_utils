#pragma once

#include <chrono>

class Clock {
    public:
        Clock(){}
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
