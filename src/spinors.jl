
mutable struct Spinor
    spins::Array{Array{Int}}
    coefficients::Array{Complex}
end

function +(spinorA::Spinor, spinorB::Spinor)
    spinsC = deepcopy(spinorA.spins)
    coefficientsC = deepcopy(spinorA.coefficients)
    for (coeffB,spinB) in zip(spinorB.coefficients,spinorB.spins)
        n = findfirst(isequal(spinB),spinorA.spins)
        if n !== nothing
            coefficientsC[n] += coeffB
        else
            push!(spinsC,spinB)
            push!(coefficientsC,coeffB)
        end
    end
    Spinor(spinsC,coefficientsC)         
end
function *(spinorA::Spinor, spinorB::Spinor)
    result = 0
    conj!(spinorA.coefficients)
    for (coeffA,spinA) in zip(spinorA.coefficients,spinorA.spins)
        for (coeffB,spinB) in zip(spinorB.coefficients,spinorB.spins)
            if spinA == spinB
                result += coeffA*coeffB
            end
        end
    end
    result
end
*(α::Any, spinorB::Spinor) = Spinor(spinorB.spins,α.*spinorB.coefficients)
*(spinorB::Spinor, α::Any) = α*spinorB