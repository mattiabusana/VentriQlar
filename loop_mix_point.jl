function loop_mix_point(po2_ven_input, pco2_ven_input, hb_global, be_global)

    include("physiology_functions.jl")


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

    return ph_mix, sat_mix, cao2_mix, caco2_mix

end
