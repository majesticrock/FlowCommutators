#include "coulomb_transform.hpp"
#include "phonon_higher_order.hpp"
#include "phonon_first_order.hpp"
#include "commute_eta_eta.hpp"

using namespace mrock::symbolic_operators;

int main(int argc, char** argv) {
	// When running, comment out everything but the function you need.
	// By default everything is listed so that compilation of all variants can be tested.
	phonon_first_order();
	phonon_higher_order();
	coulomb_transform();
	commute_eta_eta();
	return 0;
}