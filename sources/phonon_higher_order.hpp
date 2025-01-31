#pragma once
#include <mrock/symbolic_operators/Term.hpp>
#include <vector>
#include <utility>

inline void phonon_higher_order() {
	using namespace mrock::symbolic_operators;
	const std::string PHONON_V_NAME = "V_\\mathrm{ph}";

	const std::vector<Term> eta_prime({
		Term(1, Coefficient("A_2", MomentumList({ 'K', 'L', 'P', 'Q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'K', 'L', 'P', 'Q' }), IndexSum({ Index::GeneralSpin_S, Index::GeneralSpin_SPrime }) },
			std::vector<Operator>({
				Operator::Boson(Momentum('P', -1), true),
				Operator(momentum_symbols({ MomentumSymbol(1, 'K'), MomentumSymbol(1, 'P')}), Index::GeneralSpin_S, true),
				Operator('L', 1, false, Index::GeneralSpin_SPrime, true),
				Operator(momentum_symbols({ MomentumSymbol(1, 'L'), MomentumSymbol(-1, 'Q') }), Index::GeneralSpin_SPrime, false),
				Operator(momentum_symbols({ MomentumSymbol(1, 'K'), MomentumSymbol(1, 'Q') }), Index::GeneralSpin_S, false),
			})),
		Term(-1, Coefficient("B_2", MomentumList({ 'K', 'L', 'P', 'Q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'K', 'L', 'P', 'Q' }), IndexSum({ Index::GeneralSpin_S, Index::GeneralSpin_SPrime }) },
			std::vector<Operator>({
				Operator::Boson(Momentum('P', 1), false),
				Operator(momentum_symbols({ MomentumSymbol(1, 'K'), MomentumSymbol(1, 'P')}), Index::GeneralSpin_S, true),
				Operator('L', 1, false, Index::GeneralSpin_SPrime, true),
				Operator(momentum_symbols({ MomentumSymbol(1, 'L'), MomentumSymbol(-1, 'Q') }), Index::GeneralSpin_SPrime, false),
				Operator(momentum_symbols({ MomentumSymbol(1, 'K'), MomentumSymbol(1, 'Q') }), Index::GeneralSpin_S, false),
			}))
		});
	const std::vector<Term> H_original({
		// c^+ c
		Term(1, Coefficient("\\epsilon", Momentum('k')),
			SumContainer{ MomentumSum({'k'}), IndexSum({ Index::Sigma }) },
			std::vector<Operator>({
				Operator(Momentum('k'), Index::Sigma, true),
				Operator(Momentum('k'), Index::Sigma, false),
			})),
		// b^+ b
		Term(1, Coefficient("\\omega", Momentum('q')),
			SumContainer{ MomentumSum({'q'}), {} },
			std::vector<Operator>({
				Operator::Boson(Momentum('q'), true),
				Operator::Boson(Momentum('q'), false)
			})),
		// el-ph interaction
		Term(1, Coefficient("M", MomentumList({ 'k', 'q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'k', 'q' }), IndexSum(Index::Sigma) },
			std::vector<Operator>({
				Operator::Boson(Momentum('q', -1), true),
				Operator(momentum_symbols({ MomentumSymbol(1, 'k'), MomentumSymbol(1, 'q') }), Index::Sigma, true),
				Operator('k', 1, false, Index::Sigma, false)
			})),
		Term(-1,Coefficient("M", MomentumList({ 'k', 'q' }), IndexWrapper{}, false, false, true),
			SumContainer{ MomentumSum({ 'k', 'q' }), IndexSum(Index::Sigma) },
			std::vector<Operator>({
				Operator::Boson(Momentum('q'), false),
				Operator(momentum_symbols({ MomentumSymbol(1, 'k'), MomentumSymbol(1, 'q') }), Index::Sigma, true),
				Operator('k', 1, false, Index::Sigma, false)
			})),
		// el-el interaction
		Term(1, Coefficient(PHONON_V_NAME, MomentumList({ 'k', 'l', 'q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'k', 'l', 'q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
			std::vector<Operator>({
				Operator('k', 1, false, Index::Sigma, true),
				Operator('l', 1, false, Index::SigmaPrime, true),
				Operator(momentum_symbols({ MomentumSymbol(1, 'l'), MomentumSymbol(-1, 'q') }), Index::SigmaPrime, false),
				Operator(momentum_symbols({ MomentumSymbol(1, 'k'), MomentumSymbol(1, 'q') }), Index::Sigma, false),
			}))
		});

	std::vector<Term> commutation_result;
	commutator(commutation_result, eta_prime, H_original);
	clean_up(commutation_result);

	std::erase_if(commutation_result, [](const Term& term) -> bool {
		if (term.operators.size() > 5U) return true;
		return false;
		});
	std::ranges::sort(commutation_result, [](const Term& l, const Term& r) {
		return l.operators.size() < r.operators.size();
		});
	
	for (auto& term : commutation_result) {
		if(!term.operators.front().is_fermion) {
			// Boson
			assert(term.operators.size() == 5U);
			if(term.operators.front().momentum.size() > 1U) {
				throw std::runtime_error("Expected the phonon to have only one momentum!");
			}
			const MomentumSymbol::name_type current = term.operators.front().momentum.front().name;
			if (current != 'p') {
				term.swap_momenta('p', current);
			}
			if ((term.operators.front().momentum.front().factor > 0) == term.operators.front().is_daggered) {
				term.invert_momentum_sum('p');
			}

			term.rename_momenta('q', 'x');
			// first fermion c_{k+p}^+
			Momentum* current_momentum = &(term.operators[1].momentum);
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
			// fourth should be fixed automatically */

			term.rename_indizes(term.operators[3].first_index(), Index::UndefinedIndex);
			term.rename_indizes(term.operators[1].first_index(), Index::Sigma);
			term.rename_indizes(term.operators[3].first_index(), Index::SigmaPrime);
		}
		else {
			term.rename_momenta('q', 'x');
			// first fermion c_k^+
			Momentum* current_momentum = &(term.operators[0].momentum);
			if (current_momentum->size() != 1U) {
				const MomentumSymbol select = current_momentum->front();
				if(select.factor < 0) {
					term.invert_momentum_sum(select.name);
				}
				const Momentum replacement = Momentum('k') - ((*current_momentum) - Momentum(select.name));
				term.transform_momentum_sum(select.name, replacement, 'k');
			}
			else {
				term.rename_momenta((*current_momentum)[0].name, 'k');
				if((*current_momentum)[0].factor < 0) {
					term.invert_momentum_sum('k');
				}
			}
			
			// second fermion c_l^+
			current_momentum = &(term.operators[1].momentum);
			if (current_momentum->size() != 1U) {
				const MomentumSymbol select = *std::find_if_not(current_momentum->begin(), current_momentum->end(), 
					[](const MomentumSymbol& _pair) {
						return (_pair.name == 'k');
					});
				if (select.factor < 0) {
					term.invert_momentum_sum(select.name);
				}
				const Momentum replacement = Momentum('l') - ((*current_momentum) - Momentum(select.name));
				term.transform_momentum_sum(select.name, replacement, 'l');
			}
			else {
				term.rename_momenta(current_momentum->front().name, 'l');
				if((*current_momentum)[0].factor < 0) {
					term.invert_momentum_sum('l');
				}
			}
			

			// third fermion c_{l-q}
			current_momentum = &(term.operators[2].momentum);
			const MomentumSymbol select = *std::find_if_not(current_momentum->begin(), current_momentum->end(), 
				[](const MomentumSymbol& _pair) {
					return (_pair.name == 'k' || _pair.name == 'l');
				});
			if(select.factor < 0) {
				term.invert_momentum_sum(select.name);
			}
			const Momentum replacement = Momentum('l') - Momentum('q') - ((*current_momentum) - Momentum(select.name));
			term.transform_momentum_sum(select.name, replacement, 'q');
			// fourth should be fixed automatically */

			term.rename_indizes(term.operators[2].first_index(), Index::UndefinedIndex);
			term.rename_indizes(term.operators[0].first_index(), Index::Sigma);
			term.rename_indizes(term.operators[2].first_index(), Index::SigmaPrime);
		}

		std::ranges::sort(term.sums.momenta);
		std::ranges::sort(term.sums.spins);

		if(!term.operators.front().is_fermion) {
			if (term.sums.momenta.size() > 4U) {
				term.rename_momenta(term.sums.momenta[4], 'x');
			}
		}
		else {
			if (term.sums.momenta.size() > 3U) {
				term.rename_momenta(term.sums.momenta[3], 'p');
			}
		}
	}
	clear_duplicates(commutation_result);
	std::ranges::sort(commutation_result, [](const Term& A, const Term& B) {
		if (A.operators.size() < B.operators.size()) return true;
		if (A.operators.front().is_daggered && !B.operators.front().is_daggered) return true;
		return (A.coefficients[1].name < B.coefficients[1].name);
	});

	std::cout << "\\begin{align*}\n\t[\\eta', H] = "
		<< commutation_result
		<< ". \\end{align*}" << std::endl;
	std::cout << "Similar to Lenz and Wegner, we ignore all right sides that include $V_\\mathrm{Ph}$ because those are of order $\\mathcal{O}(M^3)$."
		<< "We are left with" << std::endl;

	std::erase_if(commutation_result, [&PHONON_V_NAME](const Term& term) {
		return term.coefficients[1].name == PHONON_V_NAME;
	});

	std::cout << "\\begin{align*}\n\t[\\eta', H] \\approx "
		<< commutation_result
		<< ". \\end{align*}" << std::endl;
}