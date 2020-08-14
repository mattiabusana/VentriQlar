clearconsole()

po2_mix = 40
pco2_mix = 50
global sat_mix = severinghaus(po2_mix) * 100

for i = 1:3

    global ph_mix = ph_from_be_zander(be_global, pco2_mix, hb_global, sat_mix)
    global sat_mix = zander_saturation(po2_mix, hb_global, pco2_mix, ph_mix)
    global o2_cont_mix = oxygen_content(sat_mix, hb_global, po2_mix)
    global co2_cont_mix = douglas_to_content(ph_mix, pco2_mix, hb_global, sat_mix)

end


println(sat_mix)
