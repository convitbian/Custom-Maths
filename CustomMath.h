#ifndef CUSTOM_MATH_H
#define CUSTOM_MATH_H

#include <iostream>

using namespace std;

// for trig functions, input in degree instead of radian
const float pi = 3.1415926535897932385;
const float ln10 = 2.30258509299;

// most of the functions are based on Taylor's Series for higher speed but lower accuracy
int c_factor(int x);
float c_sqrt(float number);  // based on Quake III's inverse square root code
float c_arctan(float x, float y);  // return value -90 <= x <= 90
float c_arcsin(float ratio);  // using Lagrange Interpolation
float c_arccos(float ratio);
float c_sin(float x);
float c_cos(float x);
float c_abs(float x);
float c_toRad(float degree);
float c_toDeg(float rad);
float c_pow(float x, int exp);
float c_round(float number);
float c_floor(float number);
float c_ceil(float number);
float c_ln(float x);  // ln(x * 10^n) = ln(x) + n * ln(10)
float c_log(float base, float x);
float c_root(int nth_root, float x);  // https://en.wikipedia.org/wiki/Nth_root

// custom math functions
int c_factor(int x) {
    int num = 1;
    while(x > 1) {
        num *= x;
        x--;
    }
    return num;
}

float c_sqrt(float number) {
    const float xhalf = 0.5f * number;

    union {  // get bits for floating value
        float x;
        int i;
    } u;

    u.x = number;
    
    // more accurate but slower (I guess)
    u.i = 0x5f3759df - (u.i >> 1);  // gives initial guess y0
    return number * u.x * (1.5f - xhalf * u.x * u.x);  // Newton step, repeating increases accuracy

    // faster but not too accurate
    // u.i = (1 << 29) + (u.i >> 1) - (1 << 22); 
    // return u.x;
}

float c_arctan(float x, float y) {
    if(x == 0) {
        if(y >= 0) return 90;
        else if(y < 0) return -90;
    }
    float ratio = y / x;

    float degree_rad;
    if(c_abs(ratio) < 1) {
        for(int n = 0; n < 6; n++) degree_rad += c_pow(-1, n) * (c_pow(ratio, 2 * n + 1) / (2 * n + 1));
        return (180 / pi) * degree_rad;
    }
    else if(c_abs(ratio) > 1) {
        for(int n = 0; n < 6; n++) degree_rad += c_pow(-1, n) * (1 / (c_pow(ratio, 2 * n + 1) * (2 * n + 1)));
        return (180 / pi) * (((pi * c_abs(ratio)) / (2 * (ratio))) - degree_rad);
    }
    else if(c_abs(ratio) == 1) return 45;
    return 0;
}

float LagrangePolynomial_0_8_0_9 (float x) {
    float asin_0_8 = 53.1301;
    float asin_0_85 = 58.2117;
    float asin_0_9 = 64.158;

    return ((asin_0_8 * (x - 0.85) * (x - 0.9)) / ((0.8 - 0.85) * (0.8 - 0.9)) + 
            (asin_0_85 * (x - 0.8) * (x - 0.9)) / ((0.85 - 0.8) * (0.85 - 0.9)) +
            (asin_0_9 * (x - 0.8) * (x - 0.85)) / ((0.9 - 0.8) * (0.9 - 0.85)));
}

float LagrangePolynomial_0_9_1(float x) {
    float asin_0_9 = 64.158;
    float asin_0_925 = 67.6684;
    float asin_0_95 = 71.8051;
    float asin_0_975 = 77.1614;
    float asin_0_995 = 84.2680;
    float asin_1 = 90;

    return (asin_0_9*((x-0.925)*(x-0.95)*(x-0.975)*(x-0.995)*(x-1))/((0.9-0.925)*(0.9-0.95)*(0.9-0.975)*(0.9-0.995)*(0.9-1))+
            asin_0_925*((x-0.9)*(x-0.95)*(x-0.975)*(x-0.995)*(x-1))/((0.925-0.9)*(0.925-0.95)*(0.925-0.975)*(0.925-0.995)*(0.925-1))+
            asin_0_95*((x-0.9)*(x-0.925)*(x-0.975)*(x-0.995)*(x-1))/((0.95-0.9)*(0.95-0.925)*(0.95-0.975)*(0.95-0.995)*(0.95-1))+
            asin_0_975*((x-0.9)*(x-0.925)*(x-0.95)*(x-0.995)*(x-1))/((0.975-0.9)*(0.975-0.925)*(0.975-0.95)*(0.975-0.995)*(0.975-1))+
            asin_0_995*((x-0.9)*(x-0.925)*(x-0.95)*(x-0.975)*(x-1))/((0.995-0.9)*(0.995-0.925)*(0.995-0.95)*(0.995-0.975)*(0.995-1))+
            asin_1*((x-0.9)*(x-0.925)*(x-0.95)*(x-0.975)*(x-0.995))/((1-0.9)*(1-0.925)*(1-0.95)*(1-0.975)*(1-0.995)));
}

float c_arcsin(float ratio) {
    if(ratio >= 1) return 90;
    else if(ratio <= -1) return -90;

    if(c_abs(ratio) < 0.8) {
        float rad = 0;
        for(int n = 0; n < 7; n++) {
            rad += (c_factor(2 * n) * c_pow(ratio, 2 * n + 1)) / (c_pow(c_factor(n), 2) * c_pow(2, 2 * n) * (2 * n + 1));
        }

        float degree = c_toDeg(rad);
        return degree;
    }
    else if(c_abs(ratio) < 0.9) {
        float degree = LagrangePolynomial_0_8_0_9(c_abs(ratio));
        if(ratio < 0) return -degree;
        else return degree;
    }
    else if(c_abs(ratio) < 1) {
        float degree = LagrangePolynomial_0_9_1(c_abs(ratio));
        if(ratio < 0) return -degree;
        else return degree;
    }

    return 0;
}

float c_arccos(float ratio) {
    if(c_abs(ratio) >= 1) return 0;

    float degree = c_arcsin(ratio);

    return 90 - degree;
}

float c_sin(float x) {
    float deg = x;
    while (deg > 180 || deg < -180) {
        if(deg > 180) deg -= 360;
        else if(deg < -180) deg += 360;
    }

    float ratio = 0;
    float rad = deg * (pi / 180);

    for(int n = 0; n < 6; n++) {
        ratio += ((c_pow(-1, n) * c_pow(rad, 2 * n + 1)) / c_factor(2 * n + 1));
    }
    return ratio;
}

float c_cos(float x) {
    float deg = x;
    while (deg > 180 || deg < -180) {
        if(deg > 180) deg -= 360;
        else if(deg < -180) deg += 360;
    }

    float ratio = 0;
    float rad = deg * (pi / 180);

    for(int n = 0; n < 6; n++) {
        ratio += ((c_pow(-1, n) * c_pow(rad, 2 * n)) / c_factor(2 * n));
    }
    return ratio;
}

float c_abs(float x) {
    if(x < 0) return -x;
    return x;
}

float c_toRad(float degree) {
    return (degree / 180) * pi;
}

float c_toDeg(float rad) {
    return (rad / pi) * 180;
}

float c_pow(float x, int exp) {
    if(exp == 0) return 1;
    float num = x;
    for(int i = 1; i < c_abs(exp); i++) num *= x;

    if(exp < 0) return 1 / num;
    return num;
}

float c_round(float number) {t
    float decimal = number - (int)number;
    if(decimal * 10 >= 5) return (int)number + 1;
    
    return (int)number;
}

float c_floor(float number) {
    return (int)number;
}

float c_ceil(float number) {
    float decimal = number - (int)number;
    if(decimal > 0.00001) return (int)number + 1;

    return number;
}

float c_ln(float x) {
    float power = 0;
    if(x >= 10) {
        while(x >= 10) {
            x /= 10;
            power++;
        }
    }
    else if(x < 1) {
        while(x < 1) {
            x *= 10;
            power--;
        }
    }

    return (500 * c_root(500, x) - 500) + (power * ln10);
}

float c_log(float base, float x) {
    return c_ln(x) / c_ln(base);
}

float find_root(float nth_root, float x, float guess, int count) {
    float answer = 0;
    if(count == 10) return guess;
    else {
        answer = (((nth_root - 1) / nth_root) * guess) + ((x / nth_root) * (1 / c_pow(guess, nth_root - 1)));
    }
    return find_root(nth_root, x, answer, count + 1);
}

float c_root(int nth_root, float x) {
    float guess = 1;
    float lower_value, higher_value;
    while(true) {
        if(c_pow(guess, nth_root) < x) guess++;
        else {
            higher_value = guess;
            lower_value = guess - 1;
            break;
        }
    }

    if(c_abs(c_pow(lower_value, nth_root) - x) < c_abs(c_pow(higher_value, nth_root) - guess)) guess = lower_value;
    else guess = higher_value;

    return find_root((float)nth_root, x, guess, 0);
}

#endif
