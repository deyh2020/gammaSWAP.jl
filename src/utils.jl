const GAMMA = 6.796381511118120480641144794583e-5;

function Ham!(ham,L,jt,γm,jcouplings,γ)
    @views begin
        for i in 1:Int(L*(L-1)/2)
            ham .= ham .+ jcouplings[i].*jt[:,:,i]
        end
        ham .= ham .+ γ.*γm
        nothing
    end
end

function incorporateNoise!(j0s, γs, τs, sigmas, disGam, jSWAP, j0)
    if sigmas[2] > 0.0
        γdist = truncated(Normal(disGam,sigmas[2]*j0s[1]),0.0,Inf) # Prob. distribution of γ. σγ ∝ J^1
        for i in eachindex(γs)
            γs[i] = rand(γdist)
        end
    end
    for i in eachindex(j0s)
        j0s[i] = j0 * (1.0 + randn()*sigmas[1]) # This isn't truncated because σJ is so small ~0.01
    end
    for i in eachindex(τs)
        τs[i] = pi/(4.0 * jSWAP) + randn().*sigmas[3]
    end
end

function diagexp!(matrix)
    for i in 1:length(@view matrix[:,1])
        matrix[i,i] = exp(matrix[i,i])
    end
    nothing
end

function zeros!(A)
    for i in eachindex(A)
        A[i] = 0.0
    end
    nothing
end