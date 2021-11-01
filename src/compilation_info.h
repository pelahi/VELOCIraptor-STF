//
// ICRAR - International Centre for Radio Astronomy Research
// (c) UWA - The University of Western Australia, 2021
// Copyright by UWA (in the framework of the ICRAR)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//

/**
 * @file
 *
 * Stores VELOCIraptor compilation information
 */
#include <string>

namespace vr
{

/**
 * Returns the C++ compilation flags used for building VELOCIraptor
 * @return The C++ compilation flags used for building VELOCIraptor
 */
std::string get_cxx_flags();

/**
 * Returns the macros defined for building VELOCIraptor
 * @return The macros defined for building VELOCIraptor
 */
std::string get_macros();

/**
 * Returns the compilation date of VELOCIraptor
 * @return The compilation date of VELOCIraptor
 */
std::string get_compilation_date();

} // namespace vr