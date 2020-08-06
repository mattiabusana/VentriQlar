using DelimitedFiles
using TimerOutputs
using StatsBase

clearconsole()
to = TimerOutput()
@timeit to "Total time" begin

    include("vaq_distribution_model.jl")

    # General settings
    iteration_number = 500   #Default = 500
    iterations_block_2 = 50  #Default = 50
    comparts = 50  #Default = 50
    tolerance_1 = 0.15  #Tolerance at STEP 1, default = 15%
    tolerance_2 = 0.1  #Tolerance at STEP 2, default = 10%
    dummy_gradient = false   #True to escape from debugging and start calculations

    # Patient asset input: ###################### CHANGE PATIENT NAME ######################
    patient_file_name = "prova"
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


    # Space of random choices
    space_vq_1 = collect(log_distribution(0.1, 100, 100))
    space_vq_2 = collect(log_distribution(0.1, 100, 100))
    space_sd_1 = collect(range(0.3, 2, step = 0.1))
    space_sd_2 = collect(range(0.3, 2, step = 0.1))
    space_qt_ratio = collect(range(0, 1, step = 0.05))

    # Array for solutions and counter solutions below tolerance
    global solutions = Array{Float64}(undef, 0, 14)
    counter_step_1 = 0

    # Start STEP 1
    println("")
    println("")
    println("###### START STEP 1 ######")

    # Random sampling in the defined space
    for index = 1:iteration_number
        choice_vq_1 = sample(space_vq_1)
        choice_vq_2 = sample(space_vq_2)
        choice_sd_1 = sample(space_sd_1)
        choice_sd_2 = sample(space_sd_2)
        choice_space_qt_ratio = sample(space_qt_ratio)

        println("")
        println("------------------------------")
        println("STEP 1, iteration n° $index/$iteration_number")

        @timeit to "Single random choice, block 1" begin

            output = vaq_distribution_model(
                comparts,
                qt,
                shunt,
                choice_vq_1,
                choice_vq_2,
                choice_sd_1,
                choice_sd_2,
                choice_space_qt_ratio,
                fio2,
                hb_global,
                be_global,
                dead_space,
                po2_ven_input,
                pco2_ven_input,
                print_results = false,
            )
        end

        # Distance from the target
        if dummy_gradient == false   # It may be useful to debug file saving and the second block
            gradient_po2 = abs((po2_art_patient - output[6])) / po2_art_patient
            gradient_pco2 =
                abs((pco2_art_patient - output[9])) / pco2_art_patient
        else
            gradient_po2 = 0.14
            gradient_pco2 = 0.14
        end

        println("Gradient PO2 = ", gradient_po2)
        println("Gradient PCO2 = ", gradient_pco2)

        array_solution = Array{Float64}(undef, 1, 15)
        array_solution[1:14] = output

        vo2 = output[1]
        vco2 = output[2]

        r = vco2 / vo2

        #Adding solution if below tolerance with its index
        if gradient_po2 <= tolerance_1 &&
           gradient_pco2 <= tolerance_1 & r <= 1.1

            global counter_step_1 += 1
            array_solution[15] = counter_step_1
            global solutions = vcat(solutions, output)
            open(path, "a") do io
                writedlm(io, array_solution, ";")
            end

        end

    end


    # End STEP 1

    println("")
    println("")
    println("###### START STEP 2 ######")

    # Start STEP 2: everything is repeated essentially in the same way
    path_2 = "Results//" * patient_file_name * "_step_2.csv"
    headers =
        ["vo2" "vco2" "ve" "va" "vq_global" "po2_art" "so2_art" "ph_art" "pco2_art" "mean_qt_1" "mean_qt_2" "log_sd_1" "log_sd_2" "qt_ratio" "solution"]
    open(path_2, "w") do io
        writedlm(io, headers, ";")
    end

    # Parameters for range selection over the given solution
    percent_bound_vq = 0.2
    scans = 20
    sd_bound = 0.2
    step_sd = 0.05
    qt_bound = 0.20
    step_ratio = 0.02

    # Indexes to keep track among solutions

    row_index = 0
    length_sol = length(solutions[:, 1])

    #= For EACH SOLUTION below tolerance in STEP 1, the paramters of the VA/Q distribution are varied
    within the limits specified above. In STEP 2 the tolerance is lowered to fine tune the
    results.
    =#
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
        for index = 1:iterations_block_2
            choice_2_vq_1 = sample(dist_vq_1)
            choice_2_vq_2 = sample(dist_vq_2)
            choice_2_sd_1 = sample(dist_sd_1)
            choice_2_sd_2 = sample(dist_sd_2)
            choice_2_space_qt_ratio = sample(dist_ratio)

            println("")
            println("------------------------------")
            println("STEP 2, solution n° $row_index/$length_sol, iteration n° $index/$iterations_block_2")

            @timeit to "Single random choice, block 2" begin

                output = vaq_distribution_model(
                    comparts,
                    qt,
                    shunt,
                    choice_2_vq_1,
                    choice_2_vq_2,
                    choice_2_sd_1,
                    choice_2_sd_2,
                    choice_2_space_qt_ratio,
                    fio2,
                    hb_global,
                    be_global,
                    dead_space,
                    po2_ven_input,
                    pco2_ven_input,
                    print_results = false,
                )
            end

            if dummy_gradient == false
                gradient_po2 =
                    abs((po2_art_patient - output[6])) / po2_art_patient
                gradient_pco2 =
                    abs((pco2_art_patient - output[9])) / pco2_art_patient
            else
                gradient_po2 = 0.07
                gradient_pco2 = 0.07
            end

            println("Gradient PO2 = ", gradient_po2)
            println("Gradient PCO2 = ", gradient_pco2)

            array_solution = Array{Float64}(undef, 1, 15)
            array_solution[1:14] = output
            array_solution[15] = row_index

            if gradient_po2 <= tolerance_2 && gradient_pco2 <= tolerance_2
                open(path_2, "a") do io
                    writedlm(io, array_solution, ";")
                end

            end


        end

    end

end
println("")
println("######### END SIMULATION #########")
println("")

# Final timing
println(to)
