module SpinQubits

    using LinearAlgebra
    using Statistics
    using Distributions

    const GAMMA = 6.796381511118120480641144794583e-5;

    mutable struct Spinor
        spins::Array{Array{Int}}
        coefficients::Array{Complex}
    end


    import Base.+
    import Base.*
    spinor1 = Spinor([[1,0,0,0]],[1.0])
    spinor2 = Spinor([[1,1,0,0]],[1.0])
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
    function σx(n::Int,spinor::Spinor)
        newspinor = deepcopy(spinor)
        for state in newspinor.spins
            state[n] = (state[n] + 1) % 2
        end
        return newspinor
    end
    function σy(n::Int,spinor::Spinor)
        newspinor = deepcopy(spinor)
        for (i,state) in enumerate(newspinor.spins)
            state[n] = (state[n] + 1) % 2
            newspinor.coefficients[i] *= spinor.spins[i][n] == 1 ? im : -im
        end
        return newspinor
    end
    function σz(n::Int,spinor::Spinor)
        newspinor = deepcopy(spinor)
        for i in 1:length(newspinor.coefficients)
            newspinor.coefficients[i] *= ((spinor.spins[i][n] == 1) ? 1 : -1)
        end
        return newspinor
    end
    σiσj(i,j,spinor::Spinor) = σx(i,σx(j, spinor)) + σy(i,σy(j, spinor)) + σz(i,σz(j, spinor))
    dipoledipole(i,j,spinor::Spinor) = 1/abs(i-j)^3 * ((-1)*σx(i,σx(j, spinor)) + (-1)*σy(i,σy(j, spinor)) + 2*σz(i,σz(j, spinor)))


    function getγm(L)
        h = zeros(2^L,2^L)
        gamma = GAMMA
        basis = [Spinor([(reverse(digits(i, base=2, pad=L)))],[1]) for i in 0:2^L-1]
        for (i,spinorI) in enumerate(basis)
            for (j,spinorJ) in enumerate(basis)
                h[i, j] = sum( [k == l ? 0 : (-1/2*(spinorI*dipoledipole(k,l,spinorJ))) for k in 1:L, l in 1:L] )
            end
        end
        return h
    end

    function getjtensor(L,β)
        jTensor = zeros(2^L, 2^L, Int(L*(L-1)/2))
        gamma = GAMMA
        basis = [Spinor([(reverse(digits(i, base=2, pad=L)))],[1]) for i in 0:2^L-1]
        for (i,spinorI) in enumerate(basis)
            for (j,spinorJ) in enumerate(basis)
                baseIndex = 0
                for k in 1:L-1
                    for n in 1:L-k
                        jTensor[i, j, baseIndex + n] = β^(k-1)*spinorI*σiσj(n,n+k,spinorJ)         
                    end
                    baseIndex += L-k
                end
            end
        end
        return jTensor
    end

    function Ham!(ham,L,jt,γm,jcouplings,γ)
        @views begin
            for i in 1:Int(L*(L-1)/2)
                ham .= ham .+ jcouplings[i].*jt[:,:,i]
            end
            ham .= ham .+ γ.*γm
            nothing
        end
    end

    function diagexp!(array, result)
        for i in 1:length(@view result[:,1])
            result[i,i] = exp(array[i,i])
        end
        nothing
    end
        
    function calculateFidelities(L::Int64, β::Float64, γ0::Float64, disGam, sigmas, nReals::Int64, spacing::Float64)
    
        σJ = sigmas[1]
        σγ = sigmas[2]
        στ = sigmas[3] #deviation in nanoseconds
        j0 = 1.0
        
        jtensor = getjtensor(L,β)
        γm = getγm(L)
        D = zeros(2^L,2^L)
        R = zeros(2^L,2^L) 

        exponents = collect(range(0.0,3.0,step=spacing))
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
                        if swapIndex > k && k > 1
                            js[baseIndex + (swapIndex - k)] *= jSWAP/j0
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


    function saveFidelities(L,BETA,DISGAM,sigmas,nREALS,SPACING)

        gamString = rpad(DISGAM,4,"0")
        nRealsPrime = maximum(sigmas) == 0 ? 0 : nREALS
        sigString = string("σJ",rpad(sigmas[1],4,"0"),"_","σγ",rpad(sigmas[2],4,"0"),"_","στ",rpad(sigmas[3],4,"0"),"_",lpad(nRealsPrime,5,"0"))
        filename = string("jdata/",L,"_up_β",rpad(BETA,4,"0"),"_γ",gamString,"_",sigString,"_",SPACING,".m")
        if !isfile(filename) || parse(Int,split(read(`wc -c $filename`, String)," ")[1]) == 0
            println("No file found. Calculating...")
            f = open(filename,"w")
            theseExponents = collect(range(0.0,3.0,step=SPACING))
            data = calculateFidelities(L,BETA,0.0,DISGAM,sigmas,nREALS,SPACING)
            println("Done.")
            print(f, "{")
            for i in 1:length(theseExponents)
                print(f,"{",data[1][i],",",replace(string(data[2][i]),"e"=>"*^"),"}")
                i == length(theseExponents) ? print(f,"}") : print(f,",")
            end
            close(f)
        else
            println("File already found. Skipping.")
        end
        nothing
    end

end
