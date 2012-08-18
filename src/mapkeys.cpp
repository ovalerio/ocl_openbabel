/*
 * mapkeys.cpp
 *
 *  Created on: 2 Aug 2012
 *      Author: ovalerio
 */
#include <string>

#include <openbabel/mapkeys.h>


std::string MapKeys::TIME_BOND_CALCULATIONS = "Time Bond Calculation";
std::string MapKeys::TIME_ANGLE_CALCULATIONS = "Time Angle Calculation";
std::string MapKeys::TIME_STRBND_CALCULATIONS = "Time StrBnd Calculation";
std::string MapKeys::TIME_TORSION_CALCULATIONS = "Time Torsion Calculation";
std::string MapKeys::TIME_OOP_CALCULATIONS = "Time OOP Calculation";
std::string MapKeys::TIME_VDW_CALCULATIONS = "Time VDW Calculation";
std::string MapKeys::TIME_ELECTROSTATIC_CALCULATIONS = "Time Electrostatic Calculation";

std::string MapKeys::TOTAL_BOND_CALCULATIONS = "Num Bond Calculations";
std::string MapKeys::TOTAL_ANGLE_CALCULATIONS = "Num Angle Calculations";
std::string MapKeys::TOTAL_STRBND_CALCULATIONS = "Num StrBnd Calculations";
std::string MapKeys::TOTAL_TORSION_CALCULATIONS = "Num Torsion Calculations";
std::string MapKeys::TOTAL_OOP_CALCULATIONS = "Num OOP Calculations";
std::string MapKeys::TOTAL_VDW_CALCULATIONS = "Num VDW Calculations";
std::string MapKeys::TOTAL_ELECTROSTATIC_CALCULATIONS = "Num Electrostatic Calculations";

std::string MapKeys::MEM_BOND_CALCULATIONS = "Mem Bond Calculations";
std::string MapKeys::MEM_ANGLE_CALCULATIONS = "Mem Angle Calculations";
std::string MapKeys::MEM_STRBND_CALCULATIONS = "Mem StrBnd Calculations";
std::string MapKeys::MEM_TORSION_CALCULATIONS = "Mem Torsion Calculations";
std::string MapKeys::MEM_OOP_CALCULATIONS = "Mem OOP Calculations";
std::string MapKeys::MEM_VDW_CALCULATIONS = "Mem VDW Calculations";
std::string MapKeys::MEM_ELECTROSTATIC_CALCULATIONS = "Mem Electrostatic Calculations";

/**
 *
 * Definition of static variables is given in a separate class to avoid  multiple definition errors.
 * See : http://ubuntuforums.org/showthread.php?t=836043
 *
 * You didn't both define and declare the static member in the header file for the class, did you? Like:

 * file foo.h:
 * Code:

 * class Foo
 * {
 *  static int sInt;
 * }
 *
 * static int Foo::sInt;

 * That would cause every source file that includes foo.h to try to allocate its own static member.
 * Move the definition line to the cpp file instead.
 */
