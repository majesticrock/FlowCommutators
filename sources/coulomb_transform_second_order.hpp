#pragma once
#include <mrock/symbolic_operators/Term.hpp>
#include <vector>
#include <array>
#include <algorithm>


inline void coulomb_transform_second_order() {
	using namespace mrock::symbolic_operators;

	const Term H_C(IntFractional(1, 2),
		Coefficient("V", MomentumList({ 'k', 'l', 'q' }), IndexWrapper{}, false, false),
		SumContainer{ MomentumSum({ 'k', 'l', 'q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
		std::vector<Operator>({
			Operator(Momentum('k'), Index::Sigma, true),
			Operator(Momentum('l'), Index::SigmaPrime, true),
			Operator(Momentum("l-q"), Index::SigmaPrime, false),
			Operator(Momentum("k+q"), Index::Sigma, false),
			}));
	const Term H_C_mixed(IntFractional(1, 2),
		Coefficient("V", MomentumList({ Momentum('l'), Momentum('k'), Momentum('q', 1) }), IndexWrapper{}, false, false),
		SumContainer{ MomentumSum({ 'k', 'l', 'q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
		std::vector<Operator>({
			Operator(Momentum("k"), Index::Sigma, true),
			Operator(Momentum("l"), Index::SigmaPrime, true),
			Operator(Momentum("l-q"), Index::SigmaPrime, false),
			Operator(Momentum("k+q"), Index::Sigma, false),
			}));
	const std::vector<Term> H_sym{ H_C, H_C_mixed };

	const std::vector<Term> CUT_eta({ 
		Term(1,
			Coefficient("A", MomentumList({ 'P', 'Q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'P', 'Q' }), IndexSum(Index::GeneralSpin_S) },
			std::vector<Operator>({
				Operator::Boson(Momentum('Q', -1), true),
				Operator(momentum_symbols({ MomentumSymbol(1, 'P'), MomentumSymbol(1, 'Q') }), Index::GeneralSpin_S, true),
				Operator('P', 1, false, Index::GeneralSpin_S, false)
			})),
		Term(-1,
			Coefficient("B", MomentumList({ 'P', 'Q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'P', 'Q' }), IndexSum(Index::GeneralSpin_S) },
			std::vector<Operator>({
				 Operator::Boson(Momentum('Q'), false),
				 Operator(momentum_symbols({ MomentumSymbol(1, 'P'), MomentumSymbol(1, 'Q') }), Index::GeneralSpin_S, true),
				 Operator('P', 1, false, Index::GeneralSpin_S, false)
			})),
		Term(1, Coefficient("A_2", MomentumList({ 'K', 'L', 'P', 'Q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'K', 'L', 'P', 'Q' }), IndexSum({ Index::GeneralSpin_S, Index::GeneralSpin_SPrime }) },
			std::vector<Operator>({
				Operator::Boson(Momentum('P', -1), true),
				Operator(Momentum("K+P"), Index::GeneralSpin_S, true),
				Operator('L', 1, false, Index::GeneralSpin_SPrime, true),
				Operator(Momentum("L-Q"), Index::GeneralSpin_SPrime, false),
				Operator(Momentum("K+Q"), Index::GeneralSpin_S, false),
			})),
		Term(-1, Coefficient("B_2", MomentumList({ 'K', 'L', 'P', 'Q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'K', 'L', 'P', 'Q' }), IndexSum({ Index::GeneralSpin_S, Index::GeneralSpin_SPrime }) },
			std::vector<Operator>({
				Operator::Boson(Momentum('P', 1), false),
				Operator(Momentum("K+P"), Index::GeneralSpin_S, true),
				Operator('L', 1, false, Index::GeneralSpin_SPrime, true),
				Operator(Momentum("L-Q"), Index::GeneralSpin_SPrime, false),
				Operator(Momentum("K+Q"), Index::GeneralSpin_S, false),
			}))
	});
	std::cout << "The Coulomb interaction term is" << std::endl;
	std::cout << "\\begin{align*}\n\tH_\\mathrm{C} = "
		<< H_sym
		<< ",\n\\end{align*}" << std::endl;
	std::cout << "Here, we used the symmetry of the potential $V(k, l, q) = V(l, k, q)$ and $V(k, l, q) = V(k, l, -q)$."
		<< "We can even define\n\\begin{equation}\n\tV^\\mathrm{sym} (k, l, q) \\equiv (1/2) \\cdot (V(k, l, q) + V(l, k, -q))\n\\end{equation}\n"
		<< "because this quantity will appear throughout the flow." << std::endl;
	std::cout << "Our generator for the CUT is " << std::endl;
	std::cout << "\\begin{align*}\n\t\\eta = "
		<< CUT_eta
		<< ", \\end{align*}" << std::endl;
	std::cout << "where\n\\begin{align}\n\tA(P,Q) &= M_{P,Q} \\sgn ( \\epsilon_{P+Q} - \\epsilon_P + \\omega_{-Q} ) \\\\\n"
		<< "\tB(P,Q) &= M_{P+Q,-Q}^* \\sgn ( \\epsilon_{P+Q} - \\epsilon_P - \\omega_Q )\n\\end{align}\n" << std::endl;

	std::vector<Term> commutation_result;
	commutator(commutation_result, CUT_eta, H_sym);
	clean_up(commutation_result);

	std::erase_if(commutation_result, [](const Term& term) { return term.operators.size() > 6U; });

	for (auto& term : commutation_result) {
		assert(term.operators.size() == 5U);
		if (term.operators[1].first_index() != Index::Sigma) {
			if (term.operators[2].first_index() == Index::Sigma) {
				term.rename_indizes(term.operators[1].first_index(), Index::UndefinedIndex);
				term.rename_indizes(term.operators[2].first_index(), Index::SigmaPrime);
			}
			term.rename_indizes(term.operators[1].first_index(), Index::Sigma);
		}
		if (term.operators[2].first_index() != Index::SigmaPrime) {
			term.rename_indizes(term.operators[2].first_index(), Index::SigmaPrime);
		}

		// Boson
		term.rename_momenta('p', 'z');
		Momentum* current_momentum = &(term.operators[0].momentum);
		if(current_momentum->size() > 1U) {
			const MomentumSymbol select = current_momentum->front();
			if(select.factor < 0) {
				term.invert_momentum_sum(select.name);
			}
			const Momentum replacement = Momentum('p') - ((*current_momentum) - Momentum(select.name));
			term.transform_momentum_sum(select.name, replacement, 'p');
		}
		else {
			const MomentumSymbol::name_type current = term.operators.front().momentum.front().name;
			if (current != 'p') {
				term.swap_momenta('p', current);
			}
		}
		if ((term.operators.front().momentum.front().factor > 0) == term.operators.front().is_daggered) {
			term.invert_momentum_sum('p');
		}
		
		term.rename_momenta('q', 'x');
		// first fermion c_{k+p}^+
		current_momentum = &(term.operators[1].momentum);
		if (current_momentum->size() != 2U) {
			const MomentumSymbol select = (current_momentum->front().name != 'p' 
				? current_momentum->front() : (*current_momentum)[1]);
			if(select.factor < 0) {
				term.invert_momentum_sum(select.name);
			}
			const Momentum replacement = Momentum('k') - ((*current_momentum) - Momentum(select.name)) + Momentum('p');
			term.transform_momentum_sum(select.name, replacement, 'k');
		}
		else {
			if ((*current_momentum)[0].name == 'p' || (*current_momentum)[1].name == 'k') {
				std::swap((*current_momentum)[0], (*current_momentum)[1]);
			}
			if ((*current_momentum)[1].name != 'p') {
				term.rename_momenta((*current_momentum)[0].name, '*');
				if ((*current_momentum)[0].factor < 0) {
					term.invert_momentum_sum('*');
				}
				const MomentumSymbol select = (*current_momentum)[0];
				const Momentum replacement = Momentum('k') + Momentum('p') - Momentum((*current_momentum)[1]);
				term.transform_momentum_sum(select.name, replacement, 'k');
			}
			else {
				term.rename_momenta((*current_momentum)[0].name, 'k');
				if((*current_momentum)[0].factor < 0) {
					term.invert_momentum_sum('k');
				}
			}
		}
		
		// second fermion c_l^+
		current_momentum = &(term.operators[2].momentum);
		if (current_momentum->size() != 1U) {
			const MomentumSymbol select = *std::find_if_not(current_momentum->begin(), current_momentum->end(), 
				[](const MomentumSymbol& _pair) {
					return (_pair.name == 'p' || _pair.name == 'k');
				});
			if (select.factor < 0) {
				term.invert_momentum_sum(select.name);
			}
			const Momentum replacement = Momentum('l') - ((*current_momentum) - Momentum(select.name));
			term.transform_momentum_sum(select.name, replacement, 'l');
		}
		else {
			term.rename_momenta(current_momentum->front().name, 'l');
		}
		if((*current_momentum)[0].factor < 0) {
			term.invert_momentum_sum('l');
		}

		// third fermion c_{l-q}
		current_momentum = &(term.operators[3].momentum);
		const MomentumSymbol select = *std::find_if_not(current_momentum->begin(), current_momentum->end(), 
			[](const MomentumSymbol& _pair) {
				return (_pair.name == 'p' || _pair.name == 'k' || _pair.name == 'l');
			});
		if(select.factor < 0) {
			term.invert_momentum_sum(select.name);
		}
		const Momentum replacement = Momentum('l') - Momentum('q') - ((*current_momentum) - Momentum(select.name));
		term.transform_momentum_sum(select.name, replacement, 'q');
		// fourth should be fixed automatically

		// There should not be an x anywhere
		term.rename_momenta('x', '$');
		term.rename_momenta('z', 'x');
		/* if (term.sums.momenta.size() == 5U && term.coefficients[0].momenta[0].size() == 3U) {
			//term.invert_momentum_sum('x');
			Momentum replacement = Momentum('y') - Momentum(term.coefficients[0].momenta[0][1]);
			term.transform_momentum_sum('x', replacement, 'y');
			term.rename_momenta('y', 'x');
		} */
		

		// Sort the sum indizes to be identical
		std::ranges::sort(term.sums.momenta);
		std::ranges::sort(term.sums.spins);
	}
	clear_duplicates(commutation_result);
	std::ranges::sort(commutation_result, [](const Term& A, const Term& B) {
		if (A.operators[0].is_daggered && !B.operators[0].is_daggered) return true;
		if (A.coefficients[0].name < B.coefficients[0].name) return true;
		if (A.coefficients[1].name < B.coefficients[1].name) return true;
		return false;
	});

	std::cout << "The first commutation, while omitting terms with 6 fermionic operators, yields" << std::endl;
	std::cout << "\\begin{align*}\n\t[\\eta, H_\\mathrm{C}] = "
		<< commutation_result
		<< ". \\end{align*}" << std::endl;
	return;
	std::cout << "Here, we recognize our $V^\\mathrm{sym}$! Symplifying the expression to" << std::endl;
	for (auto& term : commutation_result) {
		Coefficient& coeff = term.coefficients[1];
		if (coeff.momenta[0].uses('l')) {
			std::swap(coeff.momenta[0], coeff.momenta[1]);
		}
		if (coeff.momenta[2] == Momentum('q', -1)) {
			coeff.momenta[2].flip_momentum();
		}
		coeff.name = "V^\\mathrm{sym}";
	}
	clear_duplicates(commutation_result);
	std::cout << "\\begin{align*}\n\t[\\eta, H_\\mathrm{C}] = "
		<< commutation_result
		<< ", \\end{align*}" << std::endl;
	std::cout << "which has no two-electron interaction terms."
		<< " Thus, we include the resulting terms into our Hamiltonian as " << std::endl;

	const std::vector<Term> H_prime(
		{
			Term(1, Coefficient("C_1", MomentumList({ 'k', 'l', 'p', 'q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'k', 'l', 'p', 'q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
			std::vector<Operator>({
				Operator::Boson(Momentum('p', -1), true),
				Operator(momentum_symbols({ MomentumSymbol(1, 'k'), MomentumSymbol(1, 'p')}), Index::Sigma, true),
				Operator('l', 1, false, Index::SigmaPrime, true),
				Operator(momentum_symbols({ MomentumSymbol(1, 'l'), MomentumSymbol(-1, 'q') }), Index::SigmaPrime, false),
				Operator(momentum_symbols({ MomentumSymbol(1, 'k'), MomentumSymbol(1, 'q') }), Index::Sigma, false),
			})),
			Term(1, Coefficient("C_2", MomentumList({ 'k', 'l', 'p', 'q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'k', 'l', 'p', 'q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
			std::vector<Operator>({
				Operator::Boson(Momentum('p', 1), false),
				Operator(momentum_symbols({ MomentumSymbol(1, 'k'), MomentumSymbol(1, 'p')}), Index::Sigma, true),
				Operator('l', 1, false, Index::SigmaPrime, true),
				Operator(momentum_symbols({ MomentumSymbol(1, 'l'), MomentumSymbol(-1, 'q') }), Index::SigmaPrime, false),
				Operator(momentum_symbols({ MomentumSymbol(1, 'k'), MomentumSymbol(1, 'q') }), Index::Sigma, false),
			})),
		});


	std::cout << "\\begin{align*}\n\tH'= "
		<< H_prime
		<< ", \\end{align*}" << std::endl;

	std::vector<Term> second_commutation_result;
	commutator(second_commutation_result, CUT_eta, H_prime);
	clean_up(second_commutation_result);
	for (auto& term : second_commutation_result) {
		term.rename_momenta('p', 'l');
		term.rename_momenta('q', 'p');
		term.rename_momenta('r', 'k');
		term.rename_momenta('s', 'q');
		term.rename_momenta('t', 'Q');
	}

	std::cout << "We commute again with $H'$ and obtain" << std::endl;
	std::cout << "\\begin{align*}\n\t[\\eta, H'] = "
		<< second_commutation_result
		<< "\\end{align*}" << std::endl;
	std::cout << "Here, we want to omit all terms with bosonic creation/annihilation operators, claiming them to be small. "
		<< "We proceed in the same manner with hexatic terms." << std::endl;
	std::cout << "It is also nice to directly sort the summations so that the operators appear as they do in the Coulomb interaction. This yields" << std::endl;

	std::erase_if(second_commutation_result, [](const Term& term) { return (term.contains_boson() || term.operators.size() == 6U); });

	for (auto& term : second_commutation_result) {
		{
			Momentum& first_momentum = term.operators[0].momentum;
			if (first_momentum != Momentum('k')) {
				if (first_momentum.momentum_list.size() == 1U) {
					const MomentumSymbol::name_type original = first_momentum.momentum_list[0].name;
					term.rename_momenta(original, 'x');
					term.rename_momenta('k', original);
					term.rename_momenta('x', 'k');
				}
				else {
					int k_pos = first_momentum.isUsed('k');
					assert(k_pos > -1);
					Momentum replacement = Momentum('x') - first_momentum + Momentum(first_momentum.momentum_list[k_pos]);
					term.transform_momentum_sum(first_momentum.momentum_list[k_pos].name, replacement, 'x');
					term.rename_momenta('x', 'k');
				}
			}
		}
		{
			Momentum& second_momentum = term.operators[1].momentum;
			if (second_momentum != Momentum('l')) {
				int l_pos = second_momentum.isUsed('l');
				assert(l_pos > -1);
				Momentum replacement = Momentum('x') - second_momentum + Momentum(second_momentum.momentum_list[l_pos]);
				term.transform_momentum_sum(second_momentum.momentum_list[l_pos].name, replacement, 'x');
				term.rename_momenta('x', 'l');
			}
		}
		{
			Momentum& third_momentum = term.operators[2].momentum;
			int q_pos = third_momentum.isUsed('q');
			assert(q_pos > -1);
			int l_pos = third_momentum.isUsed('l');
			assert(l_pos > -1);
			Momentum replacement = Momentum('x') - third_momentum + Momentum(third_momentum.momentum_list[q_pos]) + Momentum(third_momentum.momentum_list[l_pos]);
			if (third_momentum.momentum_list[q_pos].factor < 0) {
				replacement.flip_momentum();
			}
			term.transform_momentum_sum(third_momentum.momentum_list[q_pos].name, replacement, 'x');
			term.rename_momenta('x', 'q');

			q_pos = third_momentum.isUsed('q');
			assert(q_pos == 1);
			if (third_momentum.momentum_list[q_pos].factor > 0) {
				term.invert_momentum_sum('q');
			}
		}
	}
	{
		// We are allowed to change the p in any manner we like - to obtain the same notation across all terms we do the following:
		Term* term_ptr = &(second_commutation_result[0]);
		// \sum_{ p l k q } \sum_{ \sigma \sigma' } A( k+q, -k+p )  C_2( p, l, k-p, q )   c_{ k, \sigma }^\dagger  c_{ l, \sigma' }^\dagger  c_{ l-q, \sigma' } c_{ k+q, \sigma }
		Momentum replacement = Momentum('x') + Momentum('k');
		term_ptr->transform_momentum_sum('p', replacement, 'x');
		term_ptr->rename_momenta('x', 'p');

		term_ptr = &(second_commutation_result[1]);
		// \sum_{ p l k q } \sum_{ \sigma \sigma' } A( l-q, -k+p )  C_2( p, l, k-p, k-p+q )   c_{ k, \sigma }^\dagger  c_{ l, \sigma' }^\dagger  c_{ l-q, \sigma' } c_{ k+q, \sigma }
		term_ptr->transform_momentum_sum('p', replacement, 'x');
		term_ptr->rename_momenta('x', 'p');

		term_ptr = &(second_commutation_result[3]);
		// \sum_{ p l k q } \sum_{ \sigma \sigma' } B( k+l-p, -k+p )  C_1( p, k+l-p, k-p, k-p+q )   c_{ k, \sigma }^\dagger  c_{ l, \sigma' }^\dagger  c_{ l-q, \sigma' } c_{ k+q, \sigma } 
		term_ptr->transform_momentum_sum('p', replacement, 'x');
		term_ptr->rename_momenta('x', 'p');
	}

	std::cout << "\\begin{align*}\n\t[\\eta, H'] \\approx "
		<< second_commutation_result
		<< "\\end{align*}" << std::endl;

	std::cout << "Inserting the definitions for $A$ and $B$ yields" << std::endl;
	std::cout << "\\begin{align*}\n\t\\partial_{\\lambda} V(\\lambda; k, l, q) \\approx - \\sum_{p} \\bigg\\{ ";

	std::array<Coefficient, 4> M_lambdas;
	M_lambdas.fill(Coefficient("M_{\\lambda}"));
	std::array<std::array<Coefficient, 3>, 4> in_signums;
	in_signums.fill({ Coefficient("\\epsilon"), Coefficient("\\epsilon"), Coefficient("\\omega") });
	for (size_t i = 0U; i < second_commutation_result.size(); ++i) {
		// The first coefficient is always A or B
		const Coefficient& term_coeff = second_commutation_result[i].coefficients[0];
		assert(term_coeff.name == "A" || term_coeff.name == "B");
		M_lambdas[i].inversion_symmetry = false;
		M_lambdas[i].momenta.resize(2);
		if (term_coeff.name == "A") {
			M_lambdas[i].momenta[0] = term_coeff.momenta[0];
			M_lambdas[i].momenta[1] = term_coeff.momenta[1];
		}
		else {
			M_lambdas[i].momenta[0] = term_coeff.momenta[0] + term_coeff.momenta[1];
			M_lambdas[i].momenta[1] = -term_coeff.momenta[1];
			M_lambdas[i].is_daggered = true;
		}

		in_signums[i][0].momenta.push_back(term_coeff.momenta[0] + term_coeff.momenta[1]);
		in_signums[i][1].momenta.push_back(term_coeff.momenta[0]);
		in_signums[i][2].momenta.push_back(term_coeff.name == "A" ? -term_coeff.momenta[1] : term_coeff.momenta[1]);

		if (i > 0) {
			std::cout << "\t+ ";
		}
		std::cout << "&" << M_lambdas[i] << " \\sgn \\left[ ";
		for (int j = 0; j < 3; ++j) {
			if (j == 1) {
				std::cout << " - ";
			}
			else if (j == 2) {
				std::cout << (term_coeff.name == "A" ? "+" : "-");
			}
			std::cout << in_signums[i][j];
		}
		std::cout << " \\right]" << second_commutation_result[i].coefficients[1];
		if (i < second_commutation_result.size() - 1) std::cout << "\\\\\n";
	}
	std::cout << " \\bigg\\}\n\\end{align*}" << std::endl;
}