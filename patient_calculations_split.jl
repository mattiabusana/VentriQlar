# Inclusions and loading modules

using DelimitedFiles
using TimerOutputs
using StatsBase

clearconsole()

println("===========================================")
println("VENTILATION-PERFUSION MISMATCH CALCULATOR")
println("Mattia Busana M.D.")
println("Copyright 2020")
println("Version 1.0")
println("===========================================")

to = TimerOutput()

include("physiology_functions.jl")
include("fndvaq.jl")
include("gvnvaq.jl")
include("split_vaq_distribution_model.jl")
include("find_r_minimum.jl")
include("distribution_blood_compartments.jl")
include("fndvaq_modified.jl")
include("atrium_composition.jl")



# General settings
iteration_number = 100000  #Default =
iterations_block_2 = 100000  #Default =
n_comparts = 500  #Default =
tolerance_1 = 0.1  #Tolerance at STEP 1, default = 10%
tolerance_2 = 0.05  #Tolerance at STEP 2, default = 5%
dummy_gradient = false   #True to escape from debugging and start calculations

# Patient asset input: ###################### CHANGE PATIENT NAME ######################
patient_file_name = "prova2"
qt = 5
shunt = 0.02
fio2 = 0.21
hb_global = 14
be_global = 0
dead_space = 0.3
po2_ven_input = 40
pco2_ven_input = 50

# Patient asset TARGET
po2_art_patient = 90
pco2_art_patient = 40


# Header writing, STEP 1
path = "Results//" * patient_file_name * "_step_1.csv"
headers =
    ["vo2" "vco2" "ve" "va" "vq_global" "po2_art" "so2_art" "ph_art" "pco2_art" "mean_qt_1" "mean_qt_2" "log_sd_1" "log_sd_2" "qt_ratio" "solution"]
open(path, "w") do io
    writedlm(io, headers, ";")
end


# Find minimum VA/Q compartment and related R
min_vaq, min_r =
    find_r_minimum(fio2, hb_global, be_global, pco2_ven_input, pco2_ven_input)


# Extremes of distribution and distribution generation
dist_min = min_vaq
dist_max = 100
q_rest = qt - (qt * shunt)


distribution = compartments_distribution(n_comparts, dist_min, dist_max)

@timeit to "VA/Q compartments" begin

    # Generate arrays for each variable
    array_vq, array_o2, array_co2, array_ph = vaq_distribution_model_split(
        distribution,
        min_vaq,
        min_r,
        fio2,
        hb_global,
        be_global,
        dead_space,
        po2_ven_input,
        pco2_ven_input,
    )

end

# Check each compartment was found
@assert length(array_vq) == length(distribution)

# Calculate mix point
global sat_mix = severinghaus(po2_ven_input) * 100

for i = 1:3
    global ph_mix =
        ph_from_be_zander(be_global, pco2_ven_input, hb_global, sat_mix)
    global sat_mix =
        zander_saturation(po2_ven_input, hb_global, pco2_ven_input, ph_mix)
    global cao2_mix = oxygen_content(sat_mix, hb_global, po2_ven_input)
    global caco2_mix =
        douglas_to_content(ph_mix, pco2_ven_input, hb_global, sat_mix)

end

# Space of random choices
space_vq_1 = collect(log_distribution(min_vaq, 100, 100))
space_vq_2 = collect(log_distribution(min_vaq, 100, 100))
space_sd_1 = collect(range(0.3, 2, step = 0.1))
space_sd_2 = collect(range(0.3, 2, step = 0.1))
space_qt_ratio = collect(range(0, 1, step = 0.05))

# Array for solutions and counter solutions below tolerance
global solutions = Array{Float64}(undef, 0, 14)
counter_step_1 = 0

# Start STEP 1
println("")
println("")
println("###### START ITERATION STEP 1 ######")

@timeit to "Random compartments choice, STEP 1" begin

    # Random sampling in the defined space
    for index = 1:iteration_number
        choice_vq_1 = sample(space_vq_1)
        choice_vq_2 = sample(space_vq_2)
        choice_sd_1 = sample(space_sd_1)
        choice_sd_2 = sample(space_sd_2)
        choice_space_qt_ratio = sample(space_qt_ratio)


        if index % 10000 == 0
            println("------------------------------")
            println("STEP 1, iteration n° $index/$iteration_number")
        end

        # Generate blood compartments
        blood_compartments = blood_volume_compartments(
            distribution,
            qt,
            shunt,
            choice_vq_1,
            choice_vq_2,
            choice_sd_1,
            choice_sd_2,
            choice_space_qt_ratio,
        )
        output = atrium_composition(
            blood_compartments,
            qt,
            shunt,
            array_o2,
            array_co2,
            array_ph,
            cao2_mix,
            caco2_mix,
            ph_mix,
            dead_space,
            choice_vq_1,
            choice_vq_2,
            choice_sd_1,
            choice_sd_2,
            choice_space_qt_ratio,
        )


        # Distance from the target
        if dummy_gradient == false   # It may be useful to debug file saving and the second block
            gradient_po2 = abs((po2_art_patient - output[6])) / po2_art_patient
            gradient_pco2 =
                abs((pco2_art_patient - output[9])) / pco2_art_patient
        else
            gradient_po2 = 0.09
            gradient_pco2 = 0.09
        end


        array_solution = Array{Float64}(undef, 1, 15)
        array_solution[1:14] = output

        vo2 = output[1]
        vco2 = output[2]

        r = vco2 / vo2

        #Adding solution if below tolerance with its index
        if gradient_po2 <= tolerance_1 && gradient_pco2 <= tolerance_1 && r <= 1

            global counter_step_1 += 1
            array_solution[15] = counter_step_1
            global solutions = vcat(solutions, output)
            open(path, "a") do io
                writedlm(io, array_solution, ";")
            end
        end

    end

end


# End STEP 1

println("")
println("")
println("###### START ITERATION STEP 2 ######")

# Start STEP 2: everything is repeated essentially in the same way
path_2 = "Results//" * patient_file_name * "_step_2.csv"
headers =
    ["vo2" "vco2" "ve" "va" "vq_global" "po2_art" "so2_art" "ph_art" "pco2_art" "mean_qt_1" "mean_qt_2" "log_sd_1" "log_sd_2" "qt_ratio" "solution"]
open(path_2, "w") do io
    writedlm(io, headers, ";")
end

# Parameters for range selection over the given solution
percent_bound_vq = 0.1
scans = 500
sd_bound = 0.1
step_sd = 0.01
qt_bound = 0.10
step_ratio = 0.01

# Indexes to keep track among solutions
row_index = 0
length_sol = length(solutions[:, 1])

#= For EACH SOLUTION below tolerance in STEP 1, the paramters of the VA/Q distribution are varied
within the limits specified above. In STEP 2 the tolerance is lowered to fine tune the
results.
=#

@timeit to "Solutions, STEP 2" begin

    for row_solution in eachrow(solutions)
        global row_index += 1

        # Identification of the parameters

        sol_vq_1 = row_solution[10]
        sol_vq_2 = row_solution[11]
        sol_sd_1 = row_solution[12]
        sol_sd_2 = row_solution[13]
        sol_qt_ratio = row_solution[14]

        #= Calculation of the upper and lower boundaries of the new ranges among
        which the random choice takes place.
        =#
        lim_up_vq_1 = sol_vq_1 + (sol_vq_1 * percent_bound_vq)
        lim_low_vq_1 = sol_vq_1 - (sol_vq_1 * percent_bound_vq)
        dist_vq_1 = log_distribution(lim_low_vq_1, lim_up_vq_1, scans)

        lim_up_vq_2 = sol_vq_2 + (sol_vq_2 * percent_bound_vq)
        lim_low_vq_2 = sol_vq_2 - (sol_vq_2 * percent_bound_vq)
        dist_vq_2 = log_distribution(lim_low_vq_2, lim_up_vq_2, scans)

        lim_up_sd_1 = sol_sd_1 + sd_bound
        lim_low_sd_1 = sol_sd_1 - sd_bound
        dist_sd_1 = collect(range(lim_low_sd_1, lim_up_sd_1, step = step_sd))

        lim_up_sd_2 = sol_sd_2 + sd_bound
        lim_low_sd_2 = sol_sd_2 - sd_bound
        dist_sd_2 = collect(range(lim_low_sd_2, lim_up_sd_2, step = step_sd))

        lim_up_ratio = sol_qt_ratio + qt_bound
        lim_low_ratio = sol_qt_ratio - qt_bound
        if lim_low_ratio <= 0
            lim_low_ratio = 0
        end
        dist_ratio =
            collect(range(lim_low_ratio, lim_up_ratio, step = step_ratio))

        #= Within each solution, each parameter is chosen randomly in the ranges
        defined above =#
        @timeit to "Random compartments choice, STEP 2" begin

            for index = 1:iterations_block_2
                choice_2_vq_1 = sample(dist_vq_1)
                choice_2_vq_2 = sample(dist_vq_2)
                choice_2_sd_1 = sample(dist_sd_1)
                choice_2_sd_2 = sample(dist_sd_2)
                choice_2_space_qt_ratio = sample(dist_ratio)

                if index % 10000 == 0
                    println("------------------------------")
                    println("STEP 2, solution n° $row_index/$length_sol, iteration n° $index/$iterations_block_2")
                end

                blood_compartments = blood_volume_compartments(
                    distribution,
                    qt,
                    shunt,
                    choice_2_vq_1,
                    choice_2_vq_2,
                    choice_2_sd_1,
                    choice_2_sd_2,
                    choice_2_space_qt_ratio,
                )


                output = atrium_composition(
                    blood_compartments,
                    qt,
                    shunt,
                    array_o2,
                    array_co2,
                    array_ph,
                    cao2_mix,
                    caco2_mix,
                    ph_mix,
                    dead_space,
                    choice_2_vq_1,
                    choice_2_vq_2,
                    choice_2_sd_1,
                    choice_2_sd_2,
                    choice_2_space_qt_ratio,
                )


                if dummy_gradient == false
                    gradient_po2 =
                        abs((po2_art_patient - output[6])) / po2_art_patient
                    gradient_pco2 =
                        abs((pco2_art_patient - output[9])) / pco2_art_patient
                else
                    gradient_po2 = 0.04
                    gradient_pco2 = 0.04
                end

                array_solution = Array{Float64}(undef, 1, 15)
                array_solution[1:14] = output
                array_solution[15] = row_index

                vo2 = output[1]
                vco2 = output[2]

                r = vco2 / vo2

                if gradient_po2 <= tolerance_2 &&
                   gradient_pco2 <= tolerance_2 &&
                   r <= 1
                    open(path_2, "a") do io
                        writedlm(io, array_solution, ";")
                    end

                end


            end

        end

    end
end

println("")
show(to, linechars = :ascii)
