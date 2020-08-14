function vaq_distribution_model_split(
    distribution,
    min_vaq,
    min_r,
    max_r,
    fio2,
    hb_global,
    be_global,
    dead_space,
    po2_ven_input,
    pco2_ven_input,
)

    global last_r = 1


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

    global array_vq = []
    global array_o2 = []
    global array_co2 = []
    global array_ph = []

    println("")
    println("--- Calculating compartments ---")

    for (index, vq_compart) in enumerate(distribution)
        println("")
        println("$index/", length(distribution), " compartments")
        #println("VAQ = ", vq_compart)
        #println("R start = ", last_r)
        vq_calculated, o2_cont, co2_cont, ph, r_index = gvnvaq(
            vq_compart,
            fio2,
            hb_global,
            be_global,
            po2_ven_input,
            pco2_ven_input,
            last_r,
            minimum_r = min_r,
            maximum_r = max_r,
        )

        global last_r = r_index

        push!(array_vq, vq_calculated)
        push!(array_o2, o2_cont)
        push!(array_co2, co2_cont)
        push!(array_ph, ph)

    end

    return array_vq, array_o2, array_co2, array_ph

end
