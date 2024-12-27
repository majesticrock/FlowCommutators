#pragma once
#include <mrock/SymbolicOperators/Term.hpp>
#include <vector>
#include <utility>

inline void phonon_first_order() {
    using namespace mrock::SymbolicOperators;
    const std::string PHONON_V_NAME = "V_\\mathrm{ph}";

    const std::vector<Term> CUT_eta({ 
		Term(1,
			std::vector<Coefficient>({ 
				Coefficient("A", MomentumList({ 'P', 'Q' }), IndexWrapper{}, false, false),
				Coefficient("M", MomentumList({ 'P', 'Q' }), IndexWrapper{}, false, false),
			}),
			SumContainer{ MomentumSum({ 'P', 'Q' }), IndexSum(Index::GeneralSpin_S) },
			std::vector<Operator>({
				Operator::Boson(Momentum('Q', -1), true),
				Operator(momentum_pairs({ std::make_pair(1, 'P'), std::make_pair(1, 'Q') }), Index::GeneralSpin_S, true),
				Operator('P', 1, false, Index::GeneralSpin_S, false)
			})),
		Term(1,
			std::vector<Coefficient>({ 
				Coefficient("B", MomentumList({ 'P', 'Q' }), IndexWrapper{}, false, false),
				Coefficient("M", MomentumList({ Momentum("P+Q"), Momentum('Q', -1) }), IndexWrapper{}, false, false, true)
			}),
			SumContainer{ MomentumSum({ 'P', 'Q' }), IndexSum(Index::GeneralSpin_S) },
			std::vector<Operator>({
				 Operator::Boson(Momentum('Q'), false),
				 Operator(momentum_pairs({ std::make_pair(1, 'P'), std::make_pair(1, 'Q') }), Index::GeneralSpin_S, true),
				 Operator('P', 1, false, Index::GeneralSpin_S, false)
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
				Operator(momentum_pairs({ std::make_pair(1, 'k'), std::make_pair(1, 'q') }), Index::Sigma, true),
				Operator('k', 1, false, Index::Sigma, false)
			})),
		Term(1, Coefficient("M", MomentumList({ Momentum("k+q"), Momentum('q', -1) }), IndexWrapper{}, false, false, true),
			SumContainer{ MomentumSum({ 'k', 'q' }), IndexSum(Index::Sigma) },
			std::vector<Operator>({
				Operator::Boson(Momentum('q'), false),
				Operator(momentum_pairs({ std::make_pair(1, 'k'), std::make_pair(1, 'q') }), Index::Sigma, true),
				Operator('k', 1, false, Index::Sigma, false)
			})),
		// el-el interaction
		//Term(1, Coefficient(PHONON_V_NAME, MomentumList({ 'k', 'l', 'q' }), IndexWrapper{}, false, false),
		//	SumContainer{ MomentumSum({ 'k', 'l', 'q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
		//	std::vector<Operator>({
		//		Operator('k', 1, false, Index::Sigma, true),
		//		Operator('l', 1, false, Index::SigmaPrime, true),
		//		Operator(momentum_pairs({ std::make_pair(1, 'l'), std::make_pair(-1, 'q') }), Index::SigmaPrime, false),
		//		Operator(momentum_pairs({ std::make_pair(1, 'k'), std::make_pair(1, 'q') }), Index::Sigma, false),
		//	}))
		});

	std::cout << "\\begin{align*}\n\tH = "
		<< H_original
		<< "\\end{align*}" << std::endl;
	std::cout << "\\begin{align*}\n\t\\eta = "
		<< CUT_eta
		<< "\\end{align*}" << std::endl;

    std::vector<Term> commutation_result;
	commutator(commutation_result, CUT_eta, H_original);
	cleanUp(commutation_result);

	for(auto& term : commutation_result) {
		if (std::find(term.sums.spins.begin(), term.sums.spins.end(), Index::GeneralSpin_S) != term.sums.spins.end()) {
			term.rename_indizes(Index::GeneralSpin_S, Index::SigmaPrime);
		}

		if (term.operators.size() == 2U) {
			Operator const * target_op = &term.operators.front();
			if (target_op->momentum.size() > 1U) {
				const char target = target_op->momentum.back().second;
				if (target_op->momentum.back().first < 0) {
					term.invert_momentum_sum(target);
				}
				term.transform_momentum_sum(target,
					Momentum('x') - target_op->momentum + Momentum(target_op->momentum.back()), 'x');
				term.rename_momenta('x', target);
			}

			if (term.operators.front().momentum.front().second == 'q') {
				term.swap_momenta('q', 'p');
			}
		}

		if (term.operators.size() == 3U) {
			Operator const * target_op = &term.operators.front();
			assert(!target_op->is_fermion);
			if (target_op->is_daggered && target_op->momentum.front().first > 0) {
				term.invert_momentum_sum(target_op->momentum.front().second);
			}
			if (!target_op->is_daggered && target_op->momentum.front().first < 0) {
				term.invert_momentum_sum(target_op->momentum.front().second);
			}

			target_op = &term.operators.back();
			if (target_op->momentum.size() > 1U) {
				char target;
				if (target_op->momentum.back().second == term.operators.front().momentum.front().second) {
					target = target_op->momentum.front().second;
					if (target_op->momentum.front().first < 0) {
						term.invert_momentum_sum(target);
					}
				} 
				else {
					target = target_op->momentum.back().second;
					if (target_op->momentum.back().first < 0) {
						term.invert_momentum_sum(target);
					}
				}
				
				term.transform_momentum_sum(target,
					Momentum('x') - target_op->momentum + Momentum(target), 'x');
				term.rename_momenta('x', target);
			}

			if (term.operators.front().momentum.front().second != 'q') {
				term.swap_momenta('p', 'q');
			}
		}

		if (term.operators.size() == 4U) {
			if (term.count_fermions() == 2) {
				char good = '0';
				if (term.operators[2].momentum.size() == 1U) {
					good = term.operators[2].momentum.front().second;
					if (term.operators[3].momentum.size() == 1U) {
						continue;
					}				
				}
				if (term.operators[3].momentum.size() == 1U) {
					good = term.operators[3].momentum.front().second;
				}

				for (int i = 2; i < 4; ++i) {
					Operator const * target_op = &term.operators[i];
					if (target_op->momentum.size() == 1U) continue;
					momentum_pair const& target_mom = (target_op->momentum[0].second == good ? target_op->momentum[1] : target_op->momentum[0]);
					const char target = target_mom.second;
					if (target_mom.first < 0) {
						term.invert_momentum_sum(target);
					}

					term.transform_momentum_sum(target,
						Momentum('x') - target_op->momentum + Momentum(target), 'x');
					term.rename_momenta('x', target);
					good = target;
				}
			}
			else {
				term.swap_momenta('r', 'q');
			}
		}

		term.rename_momenta('p', 'k');
		term.rename_momenta('r', 'l');
	}

	std::ranges::sort(commutation_result, [](const Term& l, const Term& r) {
		if (l.count_fermions() < r.count_fermions()) return true;
		if (l.count_fermions() > r.count_fermions()) return false;
		if (l.operators.size() < r.operators.size()) return true;
		if (l.operators.size() > r.operators.size()) return false;
		if (!l.operators.front().is_daggered && r.operators.front().is_daggered) return true;
		return false;
		});

    std::cout << "\\begin{align*}\n\t[\\eta, H] = "
		<< commutation_result
		<< ". \\end{align*}" << std::endl;

	std::cout << "Following the arguments of Lenz and Wegner, we omit any terms of the form $b^{(\\dagger)} b^{(\\dagger)} c^\\dagger c$."
		<< " Thereby, we arrive at " << std::endl;
	std::erase_if(commutation_result, [](const Term& term) { return (term.operators.size() == 4U && term.count_bosons() == 2); });
	std::cout << "\\begin{align*}\n\t[\\eta, H] = "
		<< commutation_result
		<< ". \\end{align*}" << std::endl;
}