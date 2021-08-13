module SpinQubits

    using LinearAlgebra, Statistics, Distributions, Test, Plots

    import Base: +, *
    
    export calculateFidelities, saveFidelities, plotter, readmathematica

    include("spinors.jl")
    include("operators.jl")
    include("utils.jl")
    include("tensors.jl")

    function calculateFidelities(L::Int64, β::Float64, γ0::Float64, disGam, sigmas, nReals::Int64, spacing::Float64)
    
        σJ = sigmas[1]
        σγ = sigmas[2]
        στ = sigmas[3] #deviation in nanoseconds
        j0 = 1.0
        
        jtensor = getjtensor(L,β)
        γm = getγtensor(L)
        D = zeros(2^L,2^L)
        R = zeros(2^L,2^L) 

        exponents = collect(range(3.0,3.0,step=spacing))
        nIterations = (maximum(sigmas) == 0) ? 1 : nReals
        singleExpFidelities = zeros(nIterations)
        fidelities = zeros(length(exponents))
        basis = [Spinor([(reverse(digits(i, base=2, pad=L)))],[1]) for i in 0:2^L-1]
        kets::Array{Array{Int,1},1} = map(x->x.spins[1],basis)
        sequence = [collect(1:L-1);collect((L-1):-1:1)]
        
        δt = στ*1e-9
        tShortestPetta = 23e-9
        δrel = δt/tShortestPetta
        tShortestMine = pi/*(4.0 * j0 * 10.0^exponents[end])
        scaledστ = δrel*tShortestMine

        initialSpinKet = [[0]; ones(Int64,L-1)]
        numJs = Int(L*(L-1)/2) # n + (n-1) + (n-2) + ... = n(n+1)/2
        js = zeros(numJs) 
        ham::Matrix{ComplexF64} = zeros(2^L,2^L)
        diss::Matrix{ComplexF64} = Diagonal(zeros(2^L))
        udiss::Matrix{ComplexF64} = Diagonal(zeros(2^L))
        UD::Matrix{ComplexF64} = zeros(2^L,2^L)
        U::Matrix{ComplexF64} = zeros(2^L,2^L)
        trueU::Matrix{ComplexF64} = zeros(2^L,2^L)
        initKet = zeros(2^L)
        initKet[findfirst(isequal(initialSpinKet),kets)] = 1.0
        finalKet::Vector{ComplexF64} = deepcopy(initKet)
        currentKet::Vector{ComplexF64} = deepcopy(initKet)

        
        for expIndex in 1:length(exponents)
            
            exponent = exponents[expIndex]
            expIndex % 25 == 0 ? println("Calculating ",expIndex,"th exponent out of ",length(exponents)) : nothing
            jSWAP = 10.0^exponent
        
            δjs = randn(numJs,nIterations).*σJ
            if σγ == 0
                γs = fill(disGam,nIterations)
            else
                γdist = truncated(Normal(disGam,σγ*j0),0.0,Inf) 
                γs = rand(γdist,nIterations)
            end
            δτs = randn(nIterations).*scaledστ

            j0s = j0.*(1.0 .+ δjs)
            τs = pi/(4.0 * jSWAP) .+ δτs

            for i in 1:nIterations
                currentKet .= initKet
                finalKet .= initKet
            
                for j in 1:(2^L)
                    diss[j,j] = (-im * (-im*γs[i]) * τs[i])
                end
                for j in 1:(2^L)
                    udiss[j,j] = exp(diss[j,j])
                end

                for swapIndex in sequence
                    for m in 1:2^L
                        for n in 1:2^L
                        U[m,n] = 0.0
                        ham[m,n] = 0.0 
                        end
                    end
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
                    diagexp!(UD, U) #UD is the input, U is the output
                    mul!(ham,U,transpose(R)) # ham here is just used as dummy memory space, not the hamiltonian
                    mul!(U,R,ham) # ham here is just used as dummy memory space, not the hamiltonian
                    
                    mul!(trueU,U,udiss) # ham here is just used as dummy memory space, not the hamiltonian
                    mul!(finalKet,trueU,currentKet)
                end # sequence
                singleExpFidelities[i] = abs2(initKet'*finalKet)           
            end # average
            fidelities[expIndex] = mean(singleExpFidelities)
        end # exponents
        return 10 .^exponents,1 .-fidelities
    end

    calculateFidelities(3, 0.01, 0.0, 0.0, zeros(3), 10, 0.01)
    #include("IO.jl")
end
