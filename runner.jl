#!/usr/bin/env julia

using SpinQubits, ArgParse	

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--length", "-L"
            arg_type = Int
	        default = 4
        "--singlet", "-s"
            arg_type = Bool
            default = false
        "--sigmas", "-σ"
		    arg_type = String 
	        default = "0.0,0.0,0.0"
        "--dirspec"
	        arg_type = String
	        default = ""
	    "--nReals"
	        arg_type = Int
	        default = 10000
    end

    return parse_args(s)
end

args = parse_commandline()
sigmas = args["sigmas"]
sigmas = replace(sigmas, "[" => "")
sigmas = replace(sigmas, "]" => "")
sigArray = parse.(Float64,split(sigmas,","))
disGam = sigArray[2]
spacing = maximum(sigArray) == 0.0 ? 0.01 : 0.03
realSigArray = [sigArray[1], 0.0, sigArray[3]]

calculateFidelities(4, 0.01, 0.0, 0.0, zeros(3), 10, 0.01; verbose=false) # Compilation run

#saveFidelities( 4, 0.01, 0.1, zeros(3), 10, 0.01; singlet=false)

saveFidelities(args["length"], 0.01, disGam, realSigArray, args["nReals"], spacing ;index=args["dirspec"], singlet=args["singlet"])
