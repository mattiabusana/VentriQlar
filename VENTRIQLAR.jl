# Inclusions and loading modules

using DelimitedFiles
using TimerOutputs
using StatsBase

clearconsole()

println("============================================")
println("-----------------VentriQlar-----------------")
println("")
println("Ventilation-perfusion mismatch calculator")
println("------------------------------------------")
println("Mattia Busana M.D.")
println("Copyright 2020")
println("Version 1.0")
println("============================================")

to = TimerOutput()

include("physiology_functions.jl")
include("gvnvaq.jl")
include("split_vaq_distribution_model.jl")
include("find_vaq_range.jl")
include("distribution_blood_compartments.jl")
include("fndvaq_modified.jl")
include("atrium_composition.jl")



# General settings
iteration_number = 1_000_000  #Default =
n_comparts = 100  #Default =
tolerance_1 = 0.1  #Tolerance, default = 10%
dummy_gradient = false   #True to escape from debugging and start calculations

# Patient asset input
patient_file_name = "pt1"
qt = 8.064
shunt = 0.4
fio2 = 0.9
hb_global = 9.5
be_global = 15.2
dead_space = 0.3
po2_ven_input = 38
pco2_ven_input = 73.8

# Patient asset TARGET
po2_art_patient = 105
pco2_art_patient = 71


# Header writing, STEP 1
path = "Results//" * patient_file_name * ".csv"
headers =
    ["vo2" "vco2" "ve" "va" "vq_global" "po2_art" "so2_art" "ph_art" "pco2_art" "mean_qt_1" "mean_qt_2" "log_sd_1" "log_sd_2" "qt_ratio" "distance" "solution" "qt" "shunt" "fio2" "hb" "be" "dead_space" "po2_ven" "pco2_ven" "min_vaq" "po2_patient" "pco2_patient" "gradient_po2" "gradient_pco2" "max_vaq"]
open(path, "w") do io
    writedlm(io, headers, ",")
end


# Find minimum and maximum VA/Q compartment and related R
global min_vaq, min_r =
    find_vq_min(fio2, hb_global, be_global, po2_ven_input, pco2_ven_input)

global max_vaq, max_r =
    find_vq_max(fio2, hb_global, be_global, po2_ven_input, pco2_ven_input)


# Extremes of distribution and distribution generation
dist_min = min_vaq
dist_max = max_vaq
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

space_vq_1 = collect(log_distribution(min_vaq, max_vaq, 5000))
space_vq_2 = collect(log_distribution(min_vaq, max_vaq, 5000))
space_sd_1 = collect(range(0.3, 2, step = 0.005))
space_sd_2 = collect(range(0.3, 2, step = 0.005))
space_qt_ratio = collect(range(0, 1, step = 0.005))

# Array for solutions and counter solutions below tolerance
global solutions = Array{Float64}(undef, 0, 14)
counter_step_1 = 0

# Start STEP 1
println("")
println("")
println("###### START ITERATION STEP 1 ######")

@timeit to "Random compartments choice" begin

    # Random sampling in the defined space
    for index = 1:iteration_number
        @timeit to "Single iteration" begin
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
                gradient_po2 =
                    abs((po2_art_patient - output[6])) / po2_art_patient
                gradient_pco2 =
                    abs((pco2_art_patient - output[9])) / pco2_art_patient
                euclidean_distance = sqrt(
                    (po2_art_patient - output[6])^2 +
                    (pco2_art_patient - output[9])^2,
                )
            else
                gradient_po2 = 0.09
                gradient_pco2 = 0.09
            end

            array_solution = Array{Float64}(undef, 1, 30)
            array_solution[1:14] = output

            vo2 = output[1]
            vco2 = output[2]

            r = vco2 / vo2

            #Adding solution if below tolerance with its index
            if gradient_po2 <= tolerance_1 &&
               gradient_pco2 <= tolerance_1 &&
               r <= 1

                global counter_step_1 += 1
                array_solution[15] = euclidean_distance
                array_solution[16] = counter_step_1
                array_solution[17] = qt
                array_solution[18] = shunt
                array_solution[19] = fio2
                array_solution[20] = hb_global
                array_solution[21] = be_global
                array_solution[22] = dead_space
                array_solution[23] = po2_ven_input
                array_solution[24] = pco2_ven_input
                array_solution[25] = min_vaq
                array_solution[26] = po2_art_patient
                array_solution[27] = pco2_art_patient
                array_solution[28] = gradient_po2
                array_solution[29] = gradient_pco2
                array_solution[30] = max_vaq

                global solutions = vcat(solutions, output)
                open(path, "a") do io
                    writedlm(io, array_solution, ",")
                end
            end

        end

    end

end

println("")
show(to, linechars = :ascii)
