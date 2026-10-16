// minimal stand-ins for the Windows API used by the DGG/SVG code (timers only) -- macOS port
#pragma once
#include <chrono>
#include <cstdint>
typedef int64_t __int64;
typedef union { int64_t QuadPart; } LARGE_INTEGER;
inline bool QueryPerformanceFrequency(LARGE_INTEGER* f) { f->QuadPart = 1000000000LL; return true; }
inline bool QueryPerformanceCounter(LARGE_INTEGER* c) {
    c->QuadPart = std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::steady_clock::now().time_since_epoch()).count();
    return true;
}
inline unsigned long GetTickCount() {
    return (unsigned long)std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now().time_since_epoch()).count();
}
#include <cstdio>
#include <string>
inline bool DeleteFile(const char* f) { return std::remove(f) == 0; }
inline bool DeleteFile(const std::string& f) { return std::remove(f.c_str()) == 0; }
