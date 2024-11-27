#pragma once
#include <SymbolicOperators/Term.hpp>
#include <vector>
#include <utility>

inline void phonon_first_order() {
    using namespace SymbolicOperators;
    const std::string PHONON_V_NAME = "V_\\mathrm{ph}";

    const std::vector<Term> CUT_eta({ 
		Term(1,
			Coefficient("A", MomentumList({ 'P', 'Q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'P', 'Q' }), IndexSum(Index::GeneralSpin_S) },
			std::vector<Operator>({
				Operator::Boson(Momentum('Q', -1), true),
				Operator(momentum_pairs({ std::make_pair(1, 'P'), std::make_pair(1, 'Q') }), Index::GeneralSpin_S, true),
				Operator('P', 1, false, Index::GeneralSpin_S, false)
			})),
		Term(-1,
			Coefficient("B", MomentumList({ 'P', 'Q' }), IndexWrapper{}, false, false),
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
		Term(-1,Coefficient("M", MomentumList({ 'k', 'q' }), IndexWrapper{}, false, false, true),
			SumContainer{ MomentumSum({ 'k', 'q' }), IndexSum(Index::Sigma) },
			std::vector<Operator>({
				Operator::Boson(Momentum('q'), false),
				Operator(momentum_pairs({ std::make_pair(1, 'k'), std::make_pair(1, 'q') }), Index::Sigma, true),
				Operator('k', 1, false, Index::Sigma, false)
			})),
		// el-el interaction
		Term(1, Coefficient(PHONON_V_NAME, MomentumList({ 'k', 'l', 'q' }), IndexWrapper{}, false, false),
			SumContainer{ MomentumSum({ 'k', 'l', 'q' }), IndexSum({ Index::Sigma, Index::SigmaPrime }) },
			std::vector<Operator>({
				Operator('k', 1, false, Index::Sigma, true),
				Operator('l', 1, false, Index::SigmaPrime, true),
				Operator(momentum_pairs({ std::make_pair(1, 'l'), std::make_pair(-1, 'q') }), Index::SigmaPrime, false),
				Operator(momentum_pairs({ std::make_pair(1, 'k'), std::make_pair(1, 'q') }), Index::Sigma, false),
			}))
		});

    std::vector<Term> commutation_result;
	commutator(commutation_result, CUT_eta, H_original);
	cleanUp(commutation_result);

    std::ranges::sort(commutation_result, [](const Term& l, const Term& r) {
		return l.operators.size() < r.operators.size();
		});
    std::cout << "\\begin{align*}\n\t[\\eta', H] = "
		<< commutation_result
		<< ". \\end{align*}" << std::endl;
}