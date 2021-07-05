import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# theorPath = '../D1pumpBx/Pump_D1-Cs133_Probe_D1-Cs133/gamma=0.0190-DSteps=150-DScan=2.00sigma-lWidth=2.0/Pump_3-3-theta=1.57-phi=1.57-pol=0-det=0_Probe_3-3-det=0/Absorption_theta=1.57-0.79_phi=0.79-1.57_pol=0-0_gammaprobe=0.0190_lWidthpr=2.0/Absorption_Cs133-D1-PumpRabi=0.50-shift=0.0-Ge=0'
# theorPath = 'C:\\Users\\laima\\OneDrive - University of Latvia\\Pētījumi\\Atomi\\Cs magnetometrs\Calculations-v1\\D1pumpBx\\Pump_D1-Cs133_Probe_D1-Cs133\gamma=0.0190-DSteps=150-DScan=2.00sigma-lWidth=2.0\\Pump_3-3-theta=1.57-phi=1.57-pol=0-det=0_Probe_3-3-det=0\\Absorption_theta=1.57-0.79_phi=0.79-1.57_pol=0-0_gammaprobe=0.0190_lWidthpr=2.0\\Absorption_Cs133-D1-PumpRabi=0.50-shift=0.0-Ge=0'
import scipy.signal

def plotTransitionDependence(Bdirection, rabiList, transitionList, probeList=None, dline=1, x_range = None):
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




def plotRabiDependence(Bdirection, rabiList, transitionList, probeList=None, dline=1, x_range = None):
    rabiList = np.array(rabiList)
    transitionList = np.array(transitionList)
    if probeList == None:
        probeList = transitionList

    fig = plt.figure()
    fig.canvas.set_window_title(Bdirection)
    plt.suptitle(f'D{dline}-{Bdirection}')
    nrows = 2
    ncols = len(transitionList)//2

    for idx, transition in enumerate(transitionList):
        plt.subplot(nrows,ncols,idx+1)
        plt.title(f'pump {transition}, probe {probeList[idx]}')
        for rabi in rabiList:

            # theorPath = '../D1pumpBx/Pump_D1-Cs133_Probe_D1-Cs133/Absorption_Cs133-D1-PumpRabi=0.50-shift=0.0-Ge=0'
            theorPath = f'./D{dline}pump{Bdirection}/Pump_D{dline}-Cs133_Probe_D{dline}-Cs133/{transition}/Absorption_Cs133-D{dline}-PumpRabi={rabi:.2f}-shift=0.0-Ge=0'


            theorData = pd.read_csv(theorPath, header=None, delimiter=' ', names=['B', 'Total', 'comp1', 'comp2', 'diff'])

            print(theorData)



            max_x_range = 1e20
            if x_range:
                max_x_range = x_range
            condition = np.abs(theorData['B'])<=max_x_range

            # Normalize data

            # theorData['comp1_smooth'] = scipy.signal.savgol_filter(theorData['comp1'], 101, 2)
            
            min_all = min(theorData['comp1'][condition])
            min_all_index = theorData[theorData['comp1'] == min_all].index[0]
            B_min_all = theorData.loc[min_all_index,'B']

            if B_min_all < 0.0:
                max_norm = max(theorData[condition].loc[:min_all_index, 'comp1'])
            else:
                max_norm = max(theorData[condition].loc[min_all_index:, 'comp1'])

            # max_norm = min(theorData[np.abs(theorData['B']) < 0.01]['comp1'])
            zero_index = theorData[np.round(theorData['B'],2) == 0.0].index[0]
            max_norm = theorData.loc[zero_index, 'comp1']

            print(theorData[np.abs(theorData['B']) == 0.0], zero_index, max_norm)

            theorData['comp1_norm'] = theorData.loc[:, 'comp1'] / max_norm
            #
            x = theorData['B'][condition].values
            # x = np.abs(theorData['B'][condition])
            # y = theorData['comp1'][condition].values
            y = theorData['comp1_norm'][condition]

            # norm_y = np.mean([y[0], y[-1]])
            # y /= norm_y


            #Normalize data
            # max_neg = theorData['comp1']

            plt.plot(x, y, label=rabi)
            # plt.plot(x, theorData[condition]['comp1_smooth']/max_norm, label=rabi)
            plt.legend()

if __name__ == '__main__':
    rabiList = [0.5, 1, 5]#, 10, 20, 50, 100]
    transitionList = ['3-3', '3-4', '4-4', '4-3']
    probeList = ['3-3', '3-4', '4-4', '4-3']
    plotRabiDependence('Bx', rabiList, transitionList, probeList, x_range = 2)

    plt.show()