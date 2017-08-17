#!/bin/bash
g++-7 -Ofast -shared -fPIC -Wall -lpthread -std=c++11 sampleUnitCell.cpp -o sampleUnitCell.so
