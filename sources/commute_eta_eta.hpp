#pragma once
#include <mrock/symbolic_operators/Term.hpp>
#include <vector>
#include <utility>

inline void commute_eta_eta() {
    using namespace mrock::symbolic_operators;

    const std::vector<Term> CUT_eta({ 
	    Term(1,
	    	std::vector<Coefficient>({ 
	    		Coefficient("\\alpha", MomentumList({ 'k', 'l' }), IndexWrapper{}, false, false),
	    		Coefficient("M_{\\ell}", MomentumList({ 'k', 'l' }), IndexWrapper{}, false, false),
	    	}),
	    	SumContainer{ MomentumSum({ 'k', 'l' }), IndexSum(Index::Sigma) },
	    	std::vector<Operator>({
	    		Operator::Boson(Momentum('l', -1), true),
	    		Operator(momentum_symbols({ MomentumSymbol(1, 'k'), MomentumSymbol(1, 'l') }), Index::Sigma, true),
	    		Operator('k', 1, false, Index::Sigma, false)
	    	})),
	    Term(1,
	    	std::vector<Coefficient>({ 
	    		Coefficient("\\beta", MomentumList({ 'k', 'l' }), IndexWrapper{}, false, false),
	    		Coefficient("M_{\\ell}", MomentumList({ Momentum("k+l"), Momentum('l', -1) }), IndexWrapper{}, false, false, true)
	    	}),
	    	SumContainer{ MomentumSum({ 'k', 'l' }), IndexSum(Index::Sigma) },
	    	std::vector<Operator>({
	    		 Operator::Boson(Momentum('l'), false),
	    		 Operator(momentum_symbols({ MomentumSymbol(1, 'k'), MomentumSymbol(1, 'l') }), Index::Sigma, true),
	    		 Operator('k', 1, false, Index::Sigma, false)
	    	}))
	    });

    const std::vector<Term> CUT_eta_prime({ 
        Term(1,
            std::vector<Coefficient>({ 
                Coefficient("\\alpha", MomentumList({ 'p', 'q' }), IndexWrapper{}, false, false),
                Coefficient("M_{\\ell'}", MomentumList({ 'p', 'q' }), IndexWrapper{}, false, false),
            }),
            SumContainer{ MomentumSum({ 'p', 'q' }), IndexSum(Index::SigmaPrime) },
            std::vector<Operator>({
                Operator::Boson(Momentum('q', -1), true),
                Operator(momentum_symbols({ MomentumSymbol(1, 'p'), MomentumSymbol(1, 'q') }), Index::SigmaPrime, true),
                Operator('p', 1, false, Index::SigmaPrime, false)
            })),
        Term(1,
            std::vector<Coefficient>({ 
                Coefficient("\\beta", MomentumList({ 'p', 'q' }), IndexWrapper{}, false, false),
                Coefficient("M_{\\ell'}", MomentumList({ Momentum("p+q"), Momentum('q', -1) }), IndexWrapper{}, false, false, true)
            }),
            SumContainer{ MomentumSum({ 'p', 'q' }), IndexSum(Index::SigmaPrime) },
            std::vector<Operator>({
                 Operator::Boson(Momentum('q'), false),
                 Operator(momentum_symbols({ MomentumSymbol(1, 'p'), MomentumSymbol(1, 'q') }), Index::SigmaPrime, true),
                 Operator('p', 1, false, Index::SigmaPrime, false)
            }))
        });

    std::cout << "\\begin{align*}\n\t\\eta(\\ell) = "
		<< CUT_eta
		<< "\\end{align*}" << std::endl;

    std::cout << "\\begin{align*}\n\t\\eta(\\ell') = "
		<< CUT_eta_prime
		<< "\\end{align*}" << std::endl;

    std::vector<Term> commutation_result = commutator(CUT_eta, CUT_eta_prime);
    clean_up(commutation_result);

    for (auto& term : commutation_result) {
        if (term.operators.size() == 2U) {
			Operator const * target_op = &term.operators.front();
			if (target_op->momentum.size() > 1U) {
				const MomentumSymbol::name_type target = target_op->momentum.back().name;
				if (target_op->momentum.back().factor < 0) {
					term.invert_momentum_sum(target);
				}
				term.transform_momentum_sum(target,
					Momentum('x') - target_op->momentum + Momentum(target_op->momentum.back()), 'x');
				term.rename_momenta('x', target);
			}

			if (term.operators.front().momentum.front().name == 'q') {
				term.swap_momenta('q', 'p');
			}
		}

        if (term.operators.size() == 4U) {
			if (term.count_fermions() == 2) {
				MomentumSymbol::name_type good = '0';
				if (term.operators[2].momentum.size() == 1U) {
					good = term.operators[2].momentum.front().name;
					if (term.operators[3].momentum.size() == 1U) {
						continue;
					}				
				}
				if (term.operators[3].momentum.size() == 1U) {
					good = term.operators[3].momentum.front().name;
				}

				for (int i = 2; i < 4; ++i) {
					Operator const * target_op = &term.operators[i];
					if (target_op->momentum.size() == 1U) continue;
					MomentumSymbol const& target_mom = (target_op->momentum[0].name == good ? target_op->momentum[1] : target_op->momentum[0]);
					const MomentumSymbol::name_type target = target_mom.name;
					if (target_mom.factor < 0) {
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
    }
    clean_up(commutation_result);
    for (auto& term : commutation_result) {
        if (term.operators.size() == 2U) {
            auto const * coeff_momentum = &(term.coefficients.front().momenta.back().front());
            if (coeff_momentum->factor < 0) {
                if (term.operators.front().momentum.front().name != coeff_momentum->name
                    && term.operators.back().momentum.front().name != coeff_momentum->name) {
                    term.invert_momentum_sum(coeff_momentum->name);
                }
            }
        }
        else if (term.operators.size() == 4U) {
            if (term.count_fermions() == 2) {
                if (term.operators[2].momentum.front().name != 'p') {
                    term.swap_momenta(term.operators[2].momentum.front().name, term.operators[3].momentum.front().name);
                }
            }
            else {
                
            }
        }
        term.rename_momenta('p', 'k');
		term.rename_momenta('r', 'l');

        std::ranges::sort(term.sums.momenta);
    }


    std::ranges::sort(commutation_result, [](const Term& l, const Term& r) {
		if (l.count_fermions() < r.count_fermions()) return true;
		if (l.count_fermions() > r.count_fermions()) return false;
		if (l.operators.size() < r.operators.size()) return true;
		if (l.operators.size() > r.operators.size()) return false;
		if (!l.operators.front().is_daggered && r.operators.front().is_daggered) return true;
		return false;
		});

    std::cout << "\\begin{align*}\n\t[\\eta (\\ell), \\eta (\\ell')] = "
		<< commutation_result
		<< ". \\end{align*}" << std::endl;
}