import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# theorPath = '../D1pumpBx/Pump_D1-Cs133_Probe_D1-Cs133/gamma=0.0190-DSteps=150-DScan=2.00sigma-lWidth=2.0/Pump_3-3-theta=1.57-phi=1.57-pol=0-det=0_Probe_3-3-det=0/Absorption_theta=1.57-0.79_phi=0.79-1.57_pol=0-0_gammaprobe=0.0190_lWidthpr=2.0/Absorption_Cs133-D1-PumpRabi=0.50-shift=0.0-Ge=0'
# theorPath = 'C:\\Users\\laima\\OneDrive - University of Latvia\\Pētījumi\\Atomi\\Cs magnetometrs\Calculations-v1\\D1pumpBx\\Pump_D1-Cs133_Probe_D1-Cs133\gamma=0.0190-DSteps=150-DScan=2.00sigma-lWidth=2.0\\Pump_3-3-theta=1.57-phi=1.57-pol=0-det=0_Probe_3-3-det=0\\Absorption_theta=1.57-0.79_phi=0.79-1.57_pol=0-0_gammaprobe=0.0190_lWidthpr=2.0\\Absorption_Cs133-D1-PumpRabi=0.50-shift=0.0-Ge=0'
import scipy.signal


def plotTransitionDependence(Bdirection, rabiList, transitionList, probeList=None, dline=1, x_range=None):
    rabiList = np.array(rabiList)
    transitionList = np.array(transitionList)
    if probeList == None:
        probeList = transitionList
    fig = plt.figure()
    fig.canvas.set_window_title(Bdirection)
    plt.suptitle(f'D{dline}-{Bdirection}')
    nrows = 2
    ncols = len(rabiList) // 2

    for idx, rabi in enumerate(rabiList):
        plt.subplot(nrows, ncols, idx + 1)
        plt.title(f'Rabi = {rabi} MHz')
        for transition in transitionList:
            # theorPath = '../D1pumpBx/Pump_D1-Cs133_Probe_D1-Cs133/Absorption_Cs133-D1-PumpRabi=0.50-shift=0.0-Ge=0'
            theorPath = f'./D{dline}pump{Bdirection}/Pump_D{dline}-Cs133_Probe_D{dline}-Cs133/{transition}/Absorption_Cs133-D{dline}-PumpRabi={rabi:.2f}-shift=0.0-Ge=0'

            theorData = pd.read_csv(theorPath, header=None, delimiter=' ',
                                    names=['B', 'Total', 'comp1', 'comp2', 'diff'])

            print(theorData)

            max_x_range = 1e20
            if x_range:
                max_x_range = x_range
            condition = np.abs(theorData['B']) <= max_x_range

            # Normalize data

            # theorData['comp1_smooth'] = scipy.signal.savgol_filter(theorData['comp1'], 101, 2)

            min_all = min(theorData['comp1'][condition])
            min_all_index = theorData[theorData['comp1'] == min_all].index[0]
            B_min_all = theorData.loc[min_all_index, 'B']

            if B_min_all < 0.0:
                max_norm = max(theorData[condition].loc[:min_all_index, 'comp1'])
            else:
                max_norm = max(theorData[condition].loc[min_all_index:, 'comp1'])

            # max_norm = min(theorData[np.abs(theorData['B']) < 0.01]['comp1'])
            zero_index = theorData[np.round(theorData['B'], 2) == 0.0].index[0]
            max_norm = theorData.loc[zero_index, 'comp1']

            print(theorData[np.abs(theorData['B']) == 0.0], zero_index, max_norm)

            theorData['comp1_norm'] = theorData.loc[:, 'comp1'] / max_norm
            #
            x = theorData['B'][condition].values
            # x = np.abs(theorData['B'][condition])
            y = theorData['comp1'][condition].values
            # y = theorData['comp1_norm'][condition]

            norm_y = np.mean([y[0], y[-1]])
            y /= norm_y

            # Normalize data
            # max_neg = theorData['comp1']

            if transition == '3-3':
                y = 2. - y

            plt.plot(x, y, label=f'{transition}')
            # plt.plot(x, theorData[condition]['comp1_smooth']/max_norm, label=rabi)
            plt.legend()


def plotRabiDependence(Bdirection, rabiList, transitionList, probeList=None, dline=1, x_range=None, comp='comp1',
                       excitation_params=None, observation1=None, observation2=None):
    if not excitation_params:
        excitation_params = {'theta': 1.57, 'phi': 0.79, 'pol': 0}
    theta_exc = excitation_params['theta']
    phi_exc = excitation_params['phi']
    pol_exc = excitation_params['pol']
    print(f'excitation_params: {excitation_params}')
    if not observation1:
        observation1 = {'theta': 1.57, 'phi': 0.79, 'pol': 0}
    theta_obs1 = observation1['theta']
    phi_obs1 = observation1['phi']
    pol_obs1 = observation1['pol']
    if not observation2:
        observation2 = {'theta': 0.79, 'phi': 1.57, 'pol': 0}
    theta_obs2 = observation2['theta']
    phi_obs2 = observation2['phi']
    pol_obs2 = observation2['pol']


    rabiList = np.array(rabiList)
    transitionList = np.array(transitionList)
    if probeList == None:
        probeList = transitionList

    fig = plt.figure()
    fig.canvas.set_window_title(Bdirection)
    plt.suptitle(f'D{dline}-{Bdirection}-{comp}')
    nrows = 2 if len(transitionList) > 1 else 1
    ncols = len(transitionList) // 2 if len(transitionList) > 1 else 1
    print(nrows, ncols)

    dir0 = f'/home/laima/Documents/Pētījumi/LZP_Cs_3D/Calculations/'

    for idx, transition in enumerate(transitionList):
        plt.subplot(nrows, ncols, idx + 1)
        plt.title(f'pump {transition}, probe {probeList[idx]}')
        for rabi in rabiList:

            t1 = f'D1pumpBx/Pump_D1-Cs133_Probe_D1-Cs133/{transition}/' \
                 f'Absorption_Cs133-D1-PumpRabi={rabi:.2f}-shift=0.0-Ge=0'
            # theorPath = '../D1pumpBx/Pump_D1-Cs133_Probe_D1-Cs133/Absorption_Cs133-D1-PumpRabi=0.50-shift=0.0-Ge=0'
            # theorPath = f'D{dline}pump{Bdirection}/Pump_D{dline}-Cs133_Probe_D{dline}-Cs133/{transition}/' \
            #             f'Absorption_Cs133-D{dline}-PumpRabi={rabi:.2f}-shift=0.0-Ge=0 '
            # print(t1==theorPath)
            # print('D1pumpBx/Pump_D1-Cs133_Probe_D1-Cs133/3-3/Absorption_theta=1.57-0.79_phi=0.79-1.57_pol=0-0_gammaprobe=0.0190_lWidthpr=2.0/Absorption_Cs133-D1-PumpRabi=0.50-shift=0.0-Ge=0')
            # theorPath = dir0 + f'D{dline}pump{Bdirection}_probe/Pump_D{dline}-Cs133_Probe_D{dline}-Cs133/' \
            #                    f'gamma=0.0190-DSteps=150-DScan=2.00sigma-lWidth=2.0/' \
            #                    f'Pump_{transition}-theta=0.79-phi=1.57-pol=0-det=0_Probe_{transition}-det=0/' \
            #                    f'Absorption_theta=1.57-0.79_phi=0.79-1.57_pol=0-0_gammaprobe=0.0190_lWidthpr=2.0/' \
            #                    f'Absorption_Cs133-D{dline}-PumpRabi={rabi:.2f}-shift=0.0-Ge=0'
            theorPath = dir0 + f'D{dline}pump{Bdirection}_compensation/Pump_D{dline}-Cs133_Probe_D{dline}-Cs133/' \
                               f'gamma=0.0190-DSteps=150-DScan=2.00sigma-lWidth=2.0/' \
                               f'Pump_{transition}-theta={theta_exc:.2f}-phi={phi_exc:.2f}-pol={pol_exc:1d}-det=0_Probe_{transition}-det=0/' \
                               f'Absorption_theta={theta_obs1:.2f}-{theta_obs2:.2f}_phi={phi_obs1:.2f}-{phi_obs2:.2f}_pol={pol_obs1:1d}-{pol_obs2:1d}_gammaprobe=0.0190_lWidthpr=2.0/' \
                               f'Absorption_Cs133-D{dline}-PumpRabi={rabi:.2f}-shift=0.0-Ge=0'
            # theorPath = dir0 + f'D{dline}pump{Bdirection}_probe/Pump_D{dline}-Cs133_Probe_D{dline}-Cs133/' \
            #                    f'gamma=0.0190-DSteps=150-DScan=2.00sigma-lWidth=2.0/' \
            #                    f'Pump_{transition}-theta=1.57-phi=0.79-pol=0-det=0_Probe_{transition}-det=0/' \
            #                    f'Absorption_theta=1.57-0.79_phi=0.79-1.57_pol=0-0_gammaprobe=0.0190_lWidthpr=2.0/' \
            #                    f'Absorption_Cs133-D{dline}-PumpRabi={rabi:.2f}-shift=0.0-Ge=0'

            theorData = pd.read_csv(theorPath, header=None, delimiter=' ',
                                    names=['B', 'Total', 'comp1', 'comp2', 'diff'])
            # theorData = pd.read_csv(t1, header=None, delimiter=' ',
            #                         names=['B', 'Total', 'comp1', 'comp2', 'diff'])
            # theorData = pd.read_csv(theorPath, header=None, delimiter=' ',
            #                         names=['B', 'Total', 'comp2', 'comp1', 'diff'])

            print(theorData)

            max_x_range = 1e20
            if x_range:
                max_x_range = x_range
            condition = np.abs(theorData['B']) <= max_x_range

            # Normalize data

            # theorData['comp1_smooth'] = scipy.signal.savgol_filter(theorData['comp1'], 101, 2)

            min_all = min(theorData[comp][condition])
            min_all_index = theorData[theorData[comp] == min_all].index[0]
            B_min_all = theorData.loc[min_all_index, 'B']

            if B_min_all < 0.0:
                max_norm = max(theorData[condition].loc[:min_all_index, comp])
            else:
                max_norm = max(theorData[condition].loc[min_all_index:, comp])

            # max_norm = min(theorData[np.abs(theorData['B']) < 0.01]['comp1'])
            zero_index = theorData[np.round(theorData['B'], 2) == 0.0].index[0]
            max_norm = theorData.loc[zero_index, comp]

            print(theorData[np.abs(theorData['B']) == 0.0], zero_index, max_norm)

            theorData[f'{comp}_norm'] = theorData.loc[:, comp] / max_norm
            #
            x = theorData['B'][condition].values
            # x = np.abs(theorData['B'][condition])
            # y = theorData['comp1'][condition].values
            # y = theorData[f'{comp}_norm'][condition]
            y = theorData[f'{comp}'][condition]

            # norm_y = np.mean([y[0], y[-1]])
            # y /= norm_y

            # Normalize data
            # max_neg = theorData['comp1']

            plt.plot(x, -y, marker='o', label=rabi)
            # plt.plot(x, theorData[condition]['comp1_smooth']/max_norm, label=rabi)
            plt.legend()


if __name__ == '__main__':
    rabiList = [0.0, 0.1, 0.5, 1, 2, 5]  # , 10, 20, 50, 100]
    # rabiList = [0.5, 1, 5]
    transitionList = ['4-3']  # ['3-3', '3-4', '4-4', '4-3']
    probeList = ['4-3']  # ['3-3', '3-4', '4-4', '4-3']
    x_range = 10

    excitation_params = {'theta': 0, 'phi': 0., 'pol': 0}
    observation1 = {'theta': 0, 'phi': 0, 'pol': 0}
    observation2 = {'theta': 1.57, 'phi': 1.57, 'pol': 0}

    # plotRabiDependence('Bx', rabiList, transitionList, probeList, x_range=2, comp='comp1')
    plotRabiDependence('Bz', rabiList, transitionList, probeList, x_range=x_range, comp='comp1',
                       excitation_params=excitation_params, observation1=observation1, observation2=observation2)

    # plotRabiDependence('Bz', rabiList, transitionList, probeList, x_range=x_range, comp='comp2',
    #                    excitation_params=excitation_params, observation1=observation1, observation2=observation2)
    #
    # plotRabiDependence('Bz', rabiList, transitionList, probeList, x_range=x_range, comp='Total',
    #                    excitation_params=excitation_params, observation1=observation1, observation2=observation2)
    #
    # excitation_params = {'theta': 1.57, 'phi': 0.79, 'pol': 0}
    # observation1 = {'theta': 0.9553166181245093, 'phi': 1.3089969389957472, 'pol': 0}  # 45 deg
    # observation2 = {'theta': 0.6154797086703874, 'phi': 2.356194490192345, 'pol': 0}  # 60 deg
    #
    # # plotRabiDependence('Bx', rabiList, transitionList, probeList, x_range=2, comp='comp1')
    # plotRabiDependence('Bx', rabiList, transitionList, probeList, x_range=x_range, comp='comp1',
    #                    excitation_params=excitation_params, observation1=observation1, observation2=observation2)
    #
    # plotRabiDependence('Bx', rabiList, transitionList, probeList, x_range=x_range, comp='comp2',
    #                    excitation_params=excitation_params, observation1=observation1, observation2=observation2)
    #
    # plotRabiDependence('Bx', rabiList, transitionList, probeList, x_range=x_range, comp='Total',
    #                    excitation_params=excitation_params, observation1=observation1, observation2=observation2)

    plt.show()
