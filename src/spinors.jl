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
    for (coeffA,spinA) in zip(spinorA.coefficients,spinorA.spins)
        for (coeffB,spinB) in zip(spinorB.coefficients,spinorB.spins)
            if spinA == spinB
                result += conj(coeffA)*coeffB
            end
        end
    end
    return result
end
function normalize(spinor::Spinor)
    newSpinor = deepcopy(spinor)
    newSpinor.coefficients = spinor.coefficients ./sqrt(spinor*spinor) #spinor on left is already daggered. See *(::Spinor,::Spinor)
    return newSpinor
end
function isequal(spinorA::Spinor, spinorB::Spinor)
    spinorA = normalize(spinorA)
    spinorB = normalize(spinorB)
    return spinorA*spinorB ≈ 1.0
end
==(spinorA::Spinor, spinorB::Spinor) = isequal(spinorA, spinorB)

*(α::Any, spinorB::Spinor) = Spinor(spinorB.spins,α.*spinorB.coefficients)
*(spinorB::Spinor, α::Any) = α*spinorB

function getInitKet(L,singlet)
    if singlet
        initialSpinor = Spinor([[[1,0]; ones(Int64,L-2)]],[1/sqrt(2)]) + Spinor([[[0,1]; ones(Int64,L-2)]],[-1/sqrt(2)])
    else
        initialSpinor = Spinor([[[0]; ones(Int64,L-1)]],[1.0])
    end

    basis = [Spinor([(reverse(digits(i, base=2, pad=L)))],[1.0]) for i in 0:2^L-1]

    initKet = [initialSpinor] .* basis

    return initKet
end