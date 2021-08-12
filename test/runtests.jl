using SpinQubits, Test


@testset "Noiseless fidelities" begin
	@test calculateFidelities(4, 0.00, 0.0, 0.0, zeros(3), 10, 1.0)[2][end] ≈ 0.0000148240466523175 rtol=1e-8
	@test calculateFidelities(4, 0.01, 0.0, 0.0, zeros(3), 10, 1.0)[2][end] ≈ 0.0018048111455146731 rtol=1e-8
	@test calculateFidelities(4, 0.05, 0.0, 0.0, zeros(3), 10, 1.0)[2][end] ≈ 0.0396766971366405 rtol=1e-8
end

# Testing spinors

# Testing operators

# Testing utils

# Testing tensors

# Testing integration
#=@testset "Noiseless fidelities" begin
	#@test calculateFidelities(4,0.01,0.0,0.0,[0.0, 0.0, 0.0],10, 0.01)[2][end] ≈ 0.0018048111455146731
	#@test calculateFidelities(4,0.01,0.0,0.0,[0.0, 0.0, 0.0],10, 0.01)[2][end] ≈ 1.4824046652317513e-5
#end=#