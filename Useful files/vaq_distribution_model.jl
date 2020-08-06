function vaq_distribution_model(
    n_comparts,
    qt,
    shunt,
    mean_qt_1,
    mean_qt_2,
    log_sd_1,
    log_sd_2,
    qt_ratio,
    fio2,
    hb_global,
    be_global,
    dead_space,
    po2_ven_input,
    pco2_ven_input;
    print_results = true,
)

    # Timing ann code inclusion

    include("distribution_gen.jl")
    include("physiology_functions.jl")
    include("fndvaq.jl")
    include("gvnvaq.jl")

    # Limits of distribution

    min_vaq, min_r = find_r_minimum(
        fio2,
        hb_global,
        be_global,
        pco2_ven_input,
        pco2_ven_input,
    )

    dist_min = min_vaq
    dist_max = 100
    q_rest = qt - (qt * shunt)


    global last_r = 1


    # Generate distributions

    distribution, blood_compartments = distribution_gen(
        n_comparts,
        qt,
        shunt,
        mean_qt_1,
        mean_qt_2,
        log_sd_1,
        log_sd_2,
        qt_ratio,
        dist_min,
        dist_max,
    )


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

    # Single compartment calculation

    global array_vq = []
    global array_o2 = []
    global array_co2 = []
    global array_ph = []

    println("")
    println("--- START CALCULATIONS ---")

    for (index, vq_compart) in enumerate(distribution)
        println("")
        println("$index/", length(distribution), " compartments")
        vq_calculated, o2_cont, co2_cont, ph, r_index = gvnvaq(
            vq_compart,
            fio2,
            hb_global,
            be_global,
            po2_ven_input,
            pco2_ven_input,
            last_r,
            minimum_r = min_r,
        )

        global last_r = r_index

        push!(array_vq, vq_calculated)
        push!(array_o2, o2_cont)
        push!(array_co2, co2_cont)
        push!(array_ph, ph)

        #= Print single compartment composition
        println("VQ = ", vq_calculated)
        println("O2 content = ", o2_cont)
        println("CO2 content = ", co2_cont)
        println("pH = ", ph)
        println("----")
        =#
    end

    # Check each compartment was found
    @assert length(array_vq) == length(distribution)


    # Calculate the blood composition in the left atrium
    cont_o2_atrium =
        (sum(array_o2 .* blood_compartments) + cao2_mix * (shunt * qt)) / qt
    cont_co2_atrium =
        (sum(array_co2 .* blood_compartments) + caco2_mix * (shunt * qt)) / qt
    ph_atrium =
        (sum(array_ph .* blood_compartments) + ph_mix * (shunt * qt)) / qt

    global content_i = 0
    global po2_atrium = 1

    while content_i <= cont_o2_atrium
        global sat_atrium =
            zander_saturation(po2_atrium, hb_global, po2_ven_input, ph_atrium)
        global content_i = oxygen_content(sat_atrium, hb_global, po2_atrium)
        global po2_atrium += 1
    end

    pco2_atrium =
        douglas_to_pressure(ph_atrium, cont_co2_atrium, hb_global, sat_atrium)

    # Calculate global variables

    vo2 = qt * 10 * (cont_o2_atrium - cao2_mix)
    vco2 = qt * 10 * (caco2_mix - cont_co2_atrium)
    r_total = vco2 / vo2
    va_global = sum(blood_compartments .* array_vq)
    ve_global = va_global / (1 - dead_space)
    qt_global = sum(blood_compartments) + (qt * shunt)
    vq_global = va_global / qt_global

    array_results =
        [vo2 vco2 ve_global va_global vq_global po2_atrium sat_atrium ph_atrium pco2_atrium mean_qt_1 mean_qt_2 log_sd_1 log_sd_2 qt_ratio]

    if print_results == true
        println("")
        println("")
        println("VO2 = ", vo2)
        println("VCO2  = ", vco2)
        println("R total= ", r_total)
        println("---")
        println("Qt total = ", qt_global)
        println("VE = ", ve_global)
        println("VA = ", va_global)
        println("VA/Q global = ", vq_global)
        println("---")
        println("PaO2 = ", po2_atrium)
        println("SatO2 = ", sat_atrium)
        println("pH = ", ph_atrium)
        println("PaCO2 = ", pco2_atrium)
        println("")
        println("--- END CALCULATIONS ---")
        println("")

        return array_results
    else
        println("")
        println("")
        println("--- END CALCULATIONS ---")
        println("")
        return array_results
    end
end
