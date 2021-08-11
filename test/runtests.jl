using SpinQubits
using Test

# Testing spinors

# Testing operators

# Testing utils

# Testing tensors

# Testing integration
@testset "Noiseless fidelities" begin
	calculateFidelities(4,0.01,0.0,0.0,[0.0, 0.0, 0.0],10, 0.01)[2][end] ≈ 0.0018048111455146731
	calculateFidelities(4,0.01,0.0,0.0,[0.0, 0.0, 0.0],10, 0.01)[2][end] ≈ 1.4824046652317513e-5
end