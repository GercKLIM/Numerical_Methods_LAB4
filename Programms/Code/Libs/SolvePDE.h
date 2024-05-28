//
// Created by Иван on 5/24/2024.
//

#ifndef CODE_SOLVEPDE_H
#define CODE_SOLVEPDE_H

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cmath>
#include <functional>
#include "algebra.cpp"
#include "FileIO.h"
#include "PDEProblem.h"

bool LongTransScheme(const PDEProblem& problem, const std::string& filename);
#endif //CODE_SOLVEPDE_H
