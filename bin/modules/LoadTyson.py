



def LoadTyson(sim_config,wd):
    k1norm = 0.015
    k2 = 0
    k3 = 200
    k4prime = 0.018
    k4 = 180

    k5 = 0
    k6 = 1
    k7 = 0.6

    k9 = 100 * k6
    k8 = 100 * k9

    params = [k1norm, k2, k3, k4prime, k4, k5, k6, k7, k8, k9]



    Y = 0.1
    YP = 0.1
    C2 = 0.5
    CP = 0.5
    M = 0.1
    pM = 0.1

    y0 = [Y, YP, C2, CP, M, pM]  # Initial conditions for x1, x2, ..., x6


    # Solve the system of ODEs
    species_all = ['cyclin', 'cyclin-P', 'cdc2', 'cdc2-P', 'cyclin-P/cdc2', 'cyclin-P/cdc2-P']
    cc_marker = str(sim_config["cc_marker"])

    th_preinc = int(sim_config["preinc_time"])

    kwargs_default = {'th': th_preinc, 'spdata': y0, 'params': params}
    model_specs = {'species_all': species_all, 'cc_marker': cc_marker}

    return(model_specs,kwargs_default)