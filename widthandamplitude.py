import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

if __name__ == '__main__':

    dlineList = [1,2]
    rabiList = [0.5]#, 1.]#, 5.]
    # Bdirection = 'Bx'
    BdirectionList = ['Bx', 'By']


    x_range = 2.0

    transition_long_list = []
    rabi_long_list = []
    Bdirection_long_list = []
    max_abs_long_list = []
    min_abs_long_list = []
    max_B_long_list = []
    min_B_long_list = []
    slope_long_list = []
    symmetry_long_list = []
    dline_long_list = []
    amplitude_long_list = []
    width_long_list = []

    for dline in dlineList:
        if dline == 1:
            transitionList = ['3-3', '3-4', '4-3', '4-4']
        elif dline == 2:
            transitionList = ['3-2', '3-3', '3-4', '4-3', '4-4', '4-5']
        for Bdirection in BdirectionList:
            for rabi in rabiList:
                for transition in transitionList:
                    theorPath = f'./D{dline}pump{Bdirection}/Pump_D{dline}-Cs133_Probe_D{dline}-Cs133/{transition}/Absorption_Cs133-D{dline}-PumpRabi={rabi:.2f}-shift=0.0-Ge=0'

                    try:
                        theorData = pd.read_csv(theorPath, header=None, delimiter=' ',
                                            names=['B', 'Total', 'comp1', 'comp2', 'diff'])
                    except:
                        continue

                    max_x_range = 1e20
                    if x_range:
                        max_x_range = x_range
                    condition = np.abs(theorData['B']) <= max_x_range

                    theordata_crop = theorData[condition]
                    theordata_crop.reset_index(drop=True)

                    # print(theordata_crop)

                    max_abs = max(theordata_crop['comp1'])
                    min_abs = min(theordata_crop['comp1'])

                    max_B = theordata_crop[theordata_crop['comp1'] == max_abs]['B'].values[0]
                    min_B = theordata_crop[theordata_crop['comp1'] == min_abs]['B'].values[0]

                    width = max_B - min_B
                    amplitude = max_abs - min_abs

                    slope = amplitude / width

                    symmetry = (1-max([max_B, min_abs]))/(1-min([max_B, min_abs]))

                    print(f'rabi = {rabi}, transition = {transition}')
                    print(f'max = {max_B} G :  {max_abs}')
                    print(f'min = {min_B} G :  {min_abs}')

                    x = theordata_crop['B']
                    y = theordata_crop['comp1']

                    # plt.plot(x,y, label = f'rabi={rabi} MHz, {transition}')

                    transition_long_list.append(transition)
                    rabi_long_list.append(rabi)
                    Bdirection_long_list.append(Bdirection)
                    max_abs_long_list.append(max_abs)
                    min_abs_long_list.append(min_abs)
                    max_B_long_list.append(max_B)
                    min_B_long_list.append(min_B)
                    slope_long_list.append(slope)
                    symmetry_long_list.append(symmetry)
                    dline_long_list.append(dline)
                    amplitude_long_list.append(amplitude)
                    width_long_list.append(width)

    data = {
        'transition' : transition_long_list,
        'Rabi': rabi_long_list,
        'Bdirection' : Bdirection_long_list,
        'max Abs': max_abs_long_list,
        'min Abs': min_abs_long_list,
        'min B': min_B_long_list,
        'max B': max_B_long_list,
        'slope': slope_long_list,
        'symmetry': symmetry_long_list,
        'D line': dline_long_list,
        'amplitude': amplitude_long_list,
        'width': width_long_list
    }

    parameters_dataframe = pd.DataFrame(data)
    print(parameters_dataframe)

    max_amplitude = parameters_dataframe[parameters_dataframe['amplitude'] == max(parameters_dataframe['amplitude'])]
    print(max_amplitude)

    plt.plot(parameters_dataframe['amplitude'], label=f"{parameters_dataframe[['D line','Bdirection','transition'] ]}")

    plt.legend()
    plt.show()