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
#include "FileIO.h"
#include "PDEProblem.h"

bool LongTransScheme(const PDEProblem& problem, const std::string& filename="UntitledTest");

bool LongTransScheme_for_tables(const PDEProblem &problem, const string &filename, const double & EPS = 1e-2);
#endif //CODE_SOLVEPDE_H
