using SpinQubits, Test

@testset verbose=true "All tests" begin
	@testset "Comparison w/ Mathematica tests" begin
		@testset "Noiseless fidelities" begin
			@test calculateFidelities(4, 0.00, 0.0, 0.0, zeros(3), 10, 3.0)[2][end] ≈ 0.0000148240466523175 rtol=1e-8
			@test calculateFidelities(4, 0.01, 0.0, 0.0, zeros(3), 10, 3.0)[2][end] ≈ 0.0018048111455146731 rtol=1e-8
			@test calculateFidelities(4, 0.05, 0.0, 0.0, zeros(3), 10, 3.0)[2][end] ≈ 0.0396766971366405 rtol=1e-8
		end
		@testset "Noisy fidelities" begin
			@test calculateFidelities(4, 0.01, 0.0, 0.0, [0.01, 0.00, 0.00], 5000, 3.0)[2][end] ≈ 0.004784666103416213 rtol=5e-2
			@test calculateFidelities(4, 0.02, 0.0, 0.0, [0.01, 0.00, 0.00], 5000, 3.0)[2][end] ≈ 0.00951159418717562 rtol=5e-2
		end
	end



	@testset "Tensors" begin
		@testset "J tensors" begin
			@test all(sum(SpinQubits.getjtensor(2, 0.0),dims=3) .≈ [1.0  0.0  0.0 0.0;
											0.0 -1.0  2.0 0.0;
											0.0  2.0 -1.0 0.0;
											0.0  0.0  0.0 1.0])
			@test all(sum(SpinQubits.getjtensor(2, 0.1),dims=3) .≈ [1.0  0.0  0.0 0.0;
											0.0 -1.0  2.0 0.0;
											0.0  2.0 -1.0 0.0;
											0.0  0.0  0.0 1.0])
			@test all(sum(SpinQubits.getjtensor(3, 0.01),dims=3) .≈ [2.01 0.0	0.0	0.0	0.0	0.0	0.0	0.0;
											0.0 -0.01 2.0	0.0	0.02 0.0 0.0 0.0;
											0.0 2 -1.99 0.0 2.0	0.0	0.0	0.0;
											0.0 0.0 0.0 -0.01	0.0 2 0.02	0.0;
											0.0 0.02	2.0 0.0 -0.01 0.0 0.0 0.0;
											0.0 0.0 0.0 2.0 0.0 -1.99 2.0 0.0;
											0.0 0.0 0.0 0.02 0.0 2.0 -0.01	0.0;
											0.0 0.0 0.0 0.0 0.0 0.0 0.0 2.01])			  							  
		end
	end
end