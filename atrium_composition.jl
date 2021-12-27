function atrium_composition(
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
    mean_qt_1,
    mean_qt_2,
    log_sd_1,
    log_sd_2,
    qt_ratio,
)

################################
array_ph = 10 .^(-array_ph)
ph_mix = 10 ^(-ph_mix)

    cont_o2_atrium =
        (sum(array_o2 .* blood_compartments) + cao2_mix * (shunt * qt)) / qt
    cont_co2_atrium =
        (sum(array_co2 .* blood_compartments) + caco2_mix * (shunt * qt)) / qt
    ph_atrium =
        (sum(array_ph .* blood_compartments) + ph_mix * (shunt * qt)) / qt


    ph_atrium = -log10(ph_atrium)


#############################

    global content_i = 0
    global po2_atrium = 1

    while content_i <= cont_o2_atrium
        global sat_atrium =
            zander_saturation(po2_atrium, hb_global, po2_ven_input, ph_atrium)
        global content_i = oxygen_content(sat_atrium, hb_global, po2_atrium)
        global po2_atrium += 0.1
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

    #println("PO2 = ", po2_atrium)  #Debug only
    #println("PCO2 = ", pco2_atrium)
    #println("")

    array_results =
        [vo2 vco2 ve_global va_global vq_global po2_atrium sat_atrium ph_atrium pco2_atrium mean_qt_1 mean_qt_2 log_sd_1 log_sd_2 qt_ratio]

    return array_results

end
