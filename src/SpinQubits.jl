module SpinQubits

    using LinearAlgebra, Statistics, Distributions, Test, Plots

    import Base: +, *, isequal, ==
    
    export calculateFidelities, saveFidelities, plotter!, readmathematica

    include("spinors.jl")
    include("operators.jl")
    include("utils.jl")
    include("tensors.jl")
    include("IO.jl")

    function calculateFidelities(L::Int64, β::Float64, γ0::Float64, disGam, sigmas, nReals::Int64, spacing::Float64; singlet=false)
    
        # construct exponents; Scale στ to our units
        j0 = 1.0
        exponents = collect(range(0.0,3.0,step=spacing))
        δt = sigmas[3]*1e-9
        tShortestPetta = 23e-9
        δrel = δt/tShortestPetta
        tShortestMine = pi/*(4.0 * j0 * 10.0^exponents[end])
        sigmas[3] = δrel*tShortestMine

        nIterations = (maximum(sigmas) == 0) ? 1 : nReals
        
        # Initialize collections for building hamiltonian
        jtensor = getjtensor(L,β)
        γm = getγtensor(L)
        numJs = Int(L*(L-1)/2) # n + (n-1) + (n-2) + ... = n(n+1)/2
        js = zeros(numJs) 
        
        sequence = singlet ? [collect(2:L-1);collect((L-1):-1:2)] : [collect(1:L-1);collect((L-1):-1:1)]

        # Find initial ket; pre allocate resulting kets after evolution
        initKet = getInitKet(L, singlet)
        finalKet::Vector{ComplexF64} = deepcopy(initKet)
        currentKet::Vector{ComplexF64} = deepcopy(initKet)

        # Pre allocation
        singleExpFidelities = zeros(nIterations)
        fidelities = zeros(length(exponents))
        D = zeros(2^L,2^L)
        R = zeros(2^L,2^L) 
        ham::Matrix{ComplexF64} = zeros(2^L,2^L)
        diss::Matrix{ComplexF64} = Diagonal(zeros(2^L))
        UD::Matrix{ComplexF64} = zeros(2^L,2^L)
        U::Matrix{ComplexF64} = zeros(2^L,2^L)
        trueU::Matrix{ComplexF64} = zeros(2^L,2^L)
        j0s = fill(j0,(numJs, nIterations))
        γs = fill(disGam,nIterations)
        τs = zeros(nIterations)
        for expIndex in 1:length(exponents)
            
            exponent = exponents[expIndex]
            expIndex % 25 == 0 ? println("Calculating ",expIndex,"th exponent out of ",length(exponents)) : nothing
            jSWAP = 10.0^exponent
            
            incorporateNoise!(j0s, γs, τs, sigmas, disGam, jSWAP, j0)

            for i in 1:nIterations
                currentKet .= initKet
                finalKet .= initKet
            
                for j in 1:(2^L)
                    diss[j,j] = (-im * (-im*γs[i]) * τs[i])
                end
                diagexp!(diss)

                for swapIndex in sequence
                    zeros!(U)
                    zeros!(ham)
                    currentKet .= finalKet
                    js[:] .= @view j0s[:,i]
                    baseIndex = 0
                    for k in 1:L-1 # k is the interspin distance
                        if swapIndex <= L-k
                            js[baseIndex + swapIndex] *= jSWAP/j0
                        end
                        if swapIndex > (k-1) && k > 1
                            js[baseIndex + (swapIndex - (k-1))] *= jSWAP/j0
                        end
                        baseIndex += L-k
                    end
                    Ham!(ham,L,jtensor,γm,js,γ0)

                    eigenObject = eigen!(ham) 
                    D .= Diagonal(eigenObject.values)
                    
                    R .= eigenObject.vectors
                    UD .= (-im .* D .* τs[i])
                    diagexp!(UD) 
                    mul!(ham,UD,transpose(R)) # ham here is just used as dummy memory space, not the hamiltonian
                    mul!(U,R,ham) # ham here is just used as dummy memory space, not the hamiltonian
                    
                    mul!(trueU,U,diss) # ham here is just used as dummy memory space, not the hamiltonian
                    mul!(finalKet,trueU,currentKet)
                end # sequence
                singleExpFidelities[i] = abs2(initKet'*finalKet)           
            end # average
            fidelities[expIndex] = mean(singleExpFidelities)
        end # exponents
        return 10 .^exponents,1 .-fidelities
    end
end
