#!/Library/Frameworks/Python.framework/Versions/3.6/bin/python3
import pandas as pd 

# FMO2 cutoffs
def fmo2_cutoff(fmo2_csv_file):
    # Read the FMO2 datas into a dataframe
    # Rename some columns for working with datas in a more easier way
    fmo2_file = fmo2_csv_file
    fmo2_data = pd.read_csv(fmo2_file)
    fmo2_data = pd.DataFrame(fmo2_data)
    fmo2_data = fmo2_data.rename(columns = {"corr/uncorr": "corr_uncorr", "dDIJ*VIJ,unc": "coul_bath"})
    
    # Calculate the two-body Hartree-Fock energy of each dimer
    # Then save them to a new column called E-HF-delFMO2 in the data frame
    E_HF_delFMO2 = [ row.corr_uncorr + row.coul_bath for index, row in fmo2_data.iterrows() ]
    fmo2_data["E_HF_delFMO2"] = E_HF_delFMO2
    
    # Determine the benchmark value (Ebm) for the HF energy of the two-body interaction (Delta FMO2)
    total_corr_uncorr = fmo2_data["corr_uncorr"].sum()
    total_coul_bath = fmo2_data["coul_bath"].sum()
    Ebm_HF_delFMO2 = total_corr_uncorr + total_coul_bath

    # Find the R threshold values out of the R column in the data frame
    # For each R value in the threshold list, sum the HF two-body energies if this energy is in the range of the R value 
    # Then, find the difference of this HF energy to the benchmark value (in kJ/mol)
    R_threshold = fmo2_data["R"].drop_duplicates().sort_values()
    cut_off_data = pd.DataFrame(R_threshold)

    E_HF_delFMO2_in_range = []
    for i in R_threshold:
        E_HF_for_this_R = 0.0
        for index, row in fmo2_data.iterrows():
            if row.R <= i:
                E_HF_for_this_R += row.E_HF_delFMO2
        E_HF_delFMO2_in_range.append(E_HF_for_this_R)
    cut_off_data["E_HF_delFMO2_in_range"] = E_HF_delFMO2_in_range

    diff_EHFdelFMO2_to_Ebm = []
    for j in cut_off_data["E_HF_delFMO2_in_range"]:
        diff_EHF = (j - Ebm_HF_delFMO2)*2625.5
        diff_EHF = '{:.8f}'.format(diff_EHF)
        diff_EHFdelFMO2_to_Ebm.append(diff_EHF)
    cut_off_data["diff_EHF_delFMO2_to_Ebm"] = diff_EHFdelFMO2_to_Ebm
    return(cut_off_data)

# FMO3 cutoff with 3R model
def fmo3_3r_cutoff(fmo3_csv_file):
    # Read the FMO3 datas into a dataframe
    # Rename some columns for working with datas in a more easier way
    fmo3_file = fmo3_csv_file
    fmo3_data = pd.read_csv(fmo3_file)
    fmo3_data = pd.DataFrame(fmo3_data)
    fmo3_data = fmo3_data.rename(columns = {"corr/uncorr": "corr_uncorr", "dDIJK*VIJK": "coul_bath"})
    
    # Determine the benchmark value (Ebm) for the HF energy of the two-body interaction (Delta FMO3)
    total_corr_uncorr = bz_p1_uc1_fmo3["corr_uncorr"].sum()
    total_coul_bath = bz_p1_uc1_fmo3["coul_bath"].sum()
    Ebm_HF_delFMO3 = total_corr_uncorr + total_coul_bath
    
    # Calculate the three-body Hartree-Fock energy of each trimer
    # Then save them to a new column called E-HF-delFMO3 in the data frame
    E_HF_delFMO3 = [ row.corr_uncorr + row.coul_bath for index, row in bz_p1_uc1_fmo3.iterrows() ]
    bz_p1_uc1_fmo3["E_HF_delFMO3"] = E_HF_delFMO3
    
    # Make a new data frame with the IJ index and R values from the FMO2 data frame
    IJ = [ str(row.I) + "-" + str(row.J) for index, row in fmo2_data.iterrows() ]
    IJ_R_FMO2_data["IJ"] = IJ
    IJ_R_FMO2_data["R"] = fmo2_data["R"]
    
    # Create new columns for indexing the dimers IJ, IK and JK
    IJ = [ str(row.I) + "-" + str(row.J) for index, row in fmo3_data.iterrows() ]
    IK = [ str(row.I) + "-" + str(row.K) for index, row in fmo3_data.iterrows() ]
    JK = [ str(row.J) + "-" + str(row.K) for index, row in fmo3_data.iterrows() ]
    fmo3_data["IJ"] = IJ
    fmo3_data["IK"] = IK 
    fmo3_data["JK"] = JK
    
    # Three new columns for three ranges: R_IJ, R_IK, R_JK
    # They correspond to three columns IJ, IK and JK. 
    # For each value of each column, look up this value in the IJ_R_FMO2_data, and take out the R value.
    R_IJ = []
    R_IK = []
    R_JK = [] 

    for q in fmo3_data["IJ"]: 
        for index, row in IJ_R_FMO2_data.iterrows():
            if row.IJ == q:
                R_IJ.append(row.R)

    for r in fmo3_data["IK"]:
        for index, row in IJ_R_FMO2_data.iterrows():
            if row.IJ == r:
                R_IK.append(row.R)

    for t in fmo3_data["JK"]:
        for index, row in IJ_R_FMO2_data.iterrows():
            if row.IJ == t:
                R_JK.append(row.R)

    fmo3_data["R_IJ"] = R_IJ 
    fmo3_data["R_IK"] = R_IK 
    fmo3_data["R_JK"] = R_JK
    
    # Compare R_IJ, R_IK and R_JK 
    # Sort them from lowest to highest values (R1, R2 and R3)
    R1 = []    # The lowest R value
    R2 = []    # The medium R value
    R3 = []    # The largest R value 

    for index, row in fmo3_data.iterrows():
        compare_R = []
        compare_R.append(row.R_IJ)
        compare_R.append(row.R_IK)
        compare_R.append(row.R_JK)
        compare_R.sort()
        R1.append(compare_R[0])
        R2.append(compare_R[1])
        R3.append(compare_R[2])

    fmo3_data["R1"] = R1
    fmo3_data["R2"] = R2    
    fmo3_data["R3"] = R3
    
    # 3R Model
    threeR_threshold = fmo3_data["R3"].drop_duplicates().sort_values()
    threeR_cut_off_data = pd.DataFrame(threeR_threshold)

    E_HF_delFMO3_in_range = []

    for i in threeR_threshold:
        E_HF_for_this_R = 0.0
        for index, row in fmo3_data.iterrows():
            if row.R3 <= i:
                E_HF_for_this_R += row.E_HF_delFMO3
        E_HF_delFMO3_in_range.append(E_HF_for_this_R)
    threeR_cut_off_data["E_HF_delFMO3_in_range"] = E_HF_delFMO3_in_range

    diff_EHFdelFMO3_to_Ebm = []

    for j in threeR_cut_off_data["E_HF_delFMO3_in_range"]:
        diff_EHF = (j - Ebm_HF_delFMO3)*2625.5
        diff_EHF = '{:.8f}'.format(diff_EHF)
        diff_EHFdelFMO3_to_Ebm.append(diff_EHF)
    threeR_cut_off_data["diff_EHF_delFMO3_to_Ebm"] = diff_EHFdelFMO3_to_Ebm
    
    return(threeR_cut_off_data)

# FMO3 cutoff with 2R model
def fmo3_2r_cutoff(fmo3_csv_file):
    # Read the FMO3 datas into a dataframe
    # Rename some columns for working with datas in a more easier way
    fmo3_file = fmo3_csv_file
    fmo3_data = pd.read_csv(fmo3_file)
    fmo3_data = pd.DataFrame(fmo3_data)
    fmo3_data = fmo3_data.rename(columns = {"corr/uncorr": "corr_uncorr", "dDIJK*VIJK": "coul_bath"})
    
    # Determine the benchmark value (Ebm) for the HF energy of the two-body interaction (Delta FMO3)
    total_corr_uncorr = bz_p1_uc1_fmo3["corr_uncorr"].sum()
    total_coul_bath = bz_p1_uc1_fmo3["coul_bath"].sum()
    Ebm_HF_delFMO3 = total_corr_uncorr + total_coul_bath
    
    # Calculate the three-body Hartree-Fock energy of each trimer
    # Then save them to a new column called E-HF-delFMO3 in the data frame
    E_HF_delFMO3 = [ row.corr_uncorr + row.coul_bath for index, row in bz_p1_uc1_fmo3.iterrows() ]
    bz_p1_uc1_fmo3["E_HF_delFMO3"] = E_HF_delFMO3
    
    # Make a new data frame with the IJ index and R values from the FMO2 data frame
    IJ = [ str(row.I) + "-" + str(row.J) for index, row in fmo2_data.iterrows() ]
    IJ_R_FMO2_data["IJ"] = IJ
    IJ_R_FMO2_data["R"] = fmo2_data["R"]
    
    # Create new columns for indexing the dimers IJ, IK and JK
    IJ = [ str(row.I) + "-" + str(row.J) for index, row in fmo3_data.iterrows() ]
    IK = [ str(row.I) + "-" + str(row.K) for index, row in fmo3_data.iterrows() ]
    JK = [ str(row.J) + "-" + str(row.K) for index, row in fmo3_data.iterrows() ]
    fmo3_data["IJ"] = IJ
    fmo3_data["IK"] = IK 
    fmo3_data["JK"] = JK
    
    # Three new columns for three ranges: R_IJ, R_IK, R_JK
    # They correspond to three columns IJ, IK and JK. 
    # For each value of each column, look up this value in the IJ_R_FMO2_data, and take out the R value.
    R_IJ = []
    R_IK = []
    R_JK = [] 

    for q in fmo3_data["IJ"]: 
        for index, row in IJ_R_FMO2_data.iterrows():
            if row.IJ == q:
                R_IJ.append(row.R)

    for r in fmo3_data["IK"]:
        for index, row in IJ_R_FMO2_data.iterrows():
            if row.IJ == r:
                R_IK.append(row.R)

    for t in fmo3_data["JK"]:
        for index, row in IJ_R_FMO2_data.iterrows():
            if row.IJ == t:
                R_JK.append(row.R)

    fmo3_data["R_IJ"] = R_IJ 
    fmo3_data["R_IK"] = R_IK 
    fmo3_data["R_JK"] = R_JK
    
    # Compare R_IJ, R_IK and R_JK 
    # Sort them from lowest to highest values (R1, R2 and R3)
    R1 = []    # The lowest R value
    R2 = []    # The medium R value
    R3 = []    # The largest R value 

    for index, row in fmo3_data.iterrows():
        compare_R = []
        compare_R.append(row.R_IJ)
        compare_R.append(row.R_IK)
        compare_R.append(row.R_JK)
        compare_R.sort()
        R1.append(compare_R[0])
        R2.append(compare_R[1])
        R3.append(compare_R[2])

    fmo3_data["R1"] = R1
    fmo3_data["R2"] = R2    
    fmo3_data["R3"] = R3
    
    # 2R Model
    twoR_threshold = fmo3_data["R2"].drop_duplicates().sort_values()
    twoR_cut_off_data = pd.DataFrame(twoR_threshold)

    E_HF_delFMO3_in_range = []

    for i in twoR_threshold:
        E_HF_for_this_R = 0.0
        for index, row in fmo3_data.iterrows():
            if row.R2 <= i:
                E_HF_for_this_R += row.E_HF_delFMO3
        E_HF_delFMO3_in_range.append(E_HF_for_this_R)
    twoR_cut_off_data["E_HF_delFMO3_in_range"] = E_HF_delFMO3_in_range

    diff_EHFdelFMO3_to_Ebm = []

    for j in twoR_cut_off_data["E_HF_delFMO3_in_range"]:
        diff_EHF = (j - Ebm_HF_delFMO3)*2625.5
        diff_EHF = '{:.8f}'.format(diff_EHF)
        diff_EHFdelFMO3_to_Ebm.append(diff_EHF)
    twoR_cut_off_data["diff_EHF_delFMO3_to_Ebm"] = diff_EHFdelFMO3_to_Ebm

    return(twoR_cut_off_data)

# FMO3 cutoff with 1R model
def fmo3_1r_cutoff(fmo3_csv_file):    
    # Read the FMO3 datas into a dataframe
    # Rename some columns for working with datas in a more easier way
    fmo3_file = fmo3_csv_file
    fmo3_data = pd.read_csv(fmo3_file)
    fmo3_data = pd.DataFrame(fmo3_data)
    fmo3_data = fmo3_data.rename(columns = {"corr/uncorr": "corr_uncorr", "dDIJK*VIJK": "coul_bath"})
    
    # Determine the benchmark value (Ebm) for the HF energy of the two-body interaction (Delta FMO3)
    total_corr_uncorr = bz_p1_uc1_fmo3["corr_uncorr"].sum()
    total_coul_bath = bz_p1_uc1_fmo3["coul_bath"].sum()
    Ebm_HF_delFMO3 = total_corr_uncorr + total_coul_bath
    
    # Calculate the three-body Hartree-Fock energy of each trimer
    # Then save them to a new column called E-HF-delFMO3 in the data frame
    E_HF_delFMO3 = [ row.corr_uncorr + row.coul_bath for index, row in bz_p1_uc1_fmo3.iterrows() ]
    bz_p1_uc1_fmo3["E_HF_delFMO3"] = E_HF_delFMO3
    
    # Make a new data frame with the IJ index and R values from the FMO2 data frame
    IJ = [ str(row.I) + "-" + str(row.J) for index, row in fmo2_data.iterrows() ]
    IJ_R_FMO2_data["IJ"] = IJ
    IJ_R_FMO2_data["R"] = fmo2_data["R"]
    
    # Create new columns for indexing the dimers IJ, IK and JK
    IJ = [ str(row.I) + "-" + str(row.J) for index, row in fmo3_data.iterrows() ]
    IK = [ str(row.I) + "-" + str(row.K) for index, row in fmo3_data.iterrows() ]
    JK = [ str(row.J) + "-" + str(row.K) for index, row in fmo3_data.iterrows() ]
    fmo3_data["IJ"] = IJ
    fmo3_data["IK"] = IK 
    fmo3_data["JK"] = JK
    
    # Three new columns for three ranges: R_IJ, R_IK, R_JK
    # They correspond to three columns IJ, IK and JK. 
    # For each value of each column, look up this value in the IJ_R_FMO2_data, and take out the R value.
    R_IJ = []
    R_IK = []
    R_JK = [] 

    for q in fmo3_data["IJ"]: 
        for index, row in IJ_R_FMO2_data.iterrows():
            if row.IJ == q:
                R_IJ.append(row.R)

    for r in fmo3_data["IK"]:
        for index, row in IJ_R_FMO2_data.iterrows():
            if row.IJ == r:
                R_IK.append(row.R)

    for t in fmo3_data["JK"]:
        for index, row in IJ_R_FMO2_data.iterrows():
            if row.IJ == t:
                R_JK.append(row.R)

    fmo3_data["R_IJ"] = R_IJ 
    fmo3_data["R_IK"] = R_IK 
    fmo3_data["R_JK"] = R_JK
    
    # Compare R_IJ, R_IK and R_JK 
    # Sort them from lowest to highest values (R1, R2 and R3)
    R1 = []    # The lowest R value
    R2 = []    # The medium R value
    R3 = []    # The largest R value 

    for index, row in fmo3_data.iterrows():
        compare_R = []
        compare_R.append(row.R_IJ)
        compare_R.append(row.R_IK)
        compare_R.append(row.R_JK)
        compare_R.sort()
        R1.append(compare_R[0])
        R2.append(compare_R[1])
        R3.append(compare_R[2])

    fmo3_data["R1"] = R1
    fmo3_data["R2"] = R2    
    fmo3_data["R3"] = R3
    
    # 1R Model
    oneR_threshold = fmo3_data["R1"].drop_duplicates().sort_values()
    oneR_cut_off_data = pd.DataFrame(oneR_threshold)

    E_HF_delFMO3_in_range = []

    for i in oneR_threshold:
        E_HF_for_this_R = 0.0
        for index, row in fmo3_data.iterrows():
            if row.R1 <= i:
                E_HF_for_this_R += row.E_HF_delFMO3
        E_HF_delFMO3_in_range.append(E_HF_for_this_R)
    oneR_cut_off_data["E_HF_delFMO3_in_range"] = E_HF_delFMO3_in_range

    diff_EHFdelFMO3_to_Ebm = []

    for j in oneR_cut_off_data["E_HF_delFMO3_in_range"]:
        diff_EHF = (j - Ebm_HF_delFMO3)*2625.5
        diff_EHF = '{:.8f}'.format(diff_EHF)
        diff_EHFdelFMO3_to_Ebm.append(diff_EHF)
    oneR_cut_off_data["diff_EHF_delFMO3_to_Ebm"] = diff_EHFdelFMO3_to_Ebm

    return(oneR_cut_off_data)
    