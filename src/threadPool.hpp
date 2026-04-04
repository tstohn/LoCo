// thread_pool.hpp
#pragma once
#include <vector>
#include <thread>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <atomic>
#include "loco_io.h"

class ThreadPool {
public:
    explicit ThreadPool(size_t threads) : stop(false) {
        for (size_t i = 0; i < threads; ++i) {
            workers.emplace_back([this]() {
                while (true) {
                    std::function<void()> task;

                    {
                        std::unique_lock<std::mutex> lock(this->queue_mutex);
                        this->condition.wait(lock, [this]() {
                            return this->stop || !this->tasks.empty();
                        });

                        if (this->stop && this->tasks.empty())
                            return;

                        task = std::move(this->tasks.front());
                        this->tasks.pop();
                    }

                    task();
                }
            });
        }
    }

    template<class F>
    void enqueue(F&& f) {
        active_tasks++;
        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            tasks.emplace([this, func = std::forward<F>(f)]() {
                func();
                active_tasks--;
            });
        }
        condition.notify_one();
    }

    void wait_for_tasks() {
        while (true) {
            if (active_tasks == 0) {
                std::lock_guard<std::mutex> lock(queue_mutex);
                if (tasks.empty())
                    return;
            }
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
    }

    ~ThreadPool() {
        {
            std::lock_guard<std::mutex> lock(queue_mutex);
            stop = true;
        }
        condition.notify_all();
        for (std::thread &worker : workers)
            worker.join();
    }

private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex queue_mutex;
    std::condition_variable condition;
    std::atomic<int> active_tasks{0};
    bool stop;
};