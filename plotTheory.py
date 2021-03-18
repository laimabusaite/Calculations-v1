import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# theorPath = '../D1pumpBx/Pump_D1-Cs133_Probe_D1-Cs133/gamma=0.0190-DSteps=150-DScan=2.00sigma-lWidth=2.0/Pump_3-3-theta=1.57-phi=1.57-pol=0-det=0_Probe_3-3-det=0/Absorption_theta=1.57-0.79_phi=0.79-1.57_pol=0-0_gammaprobe=0.0190_lWidthpr=2.0/Absorption_Cs133-D1-PumpRabi=0.50-shift=0.0-Ge=0'
# theorPath = 'C:\\Users\\laima\\OneDrive - University of Latvia\\Pētījumi\\Atomi\\Cs magnetometrs\Calculations-v1\\D1pumpBx\\Pump_D1-Cs133_Probe_D1-Cs133\gamma=0.0190-DSteps=150-DScan=2.00sigma-lWidth=2.0\\Pump_3-3-theta=1.57-phi=1.57-pol=0-det=0_Probe_3-3-det=0\\Absorption_theta=1.57-0.79_phi=0.79-1.57_pol=0-0_gammaprobe=0.0190_lWidthpr=2.0\\Absorption_Cs133-D1-PumpRabi=0.50-shift=0.0-Ge=0'


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

            plt.plot(theorData['B'][condition], theorData['comp1'][condition], label=rabi)
            plt.legend()

if __name__ == '__main__':
    rabiList = [0.5, 1, 5]#, 10, 20, 50, 100]
    transitionList = ['3-3', '3-4', '4-4', '4-3']
    probeList = ['3-3', '3-4', '4-4', '4-3']
    plotRabiDependence('Bx', rabiList, transitionList, probeList)

    plt.show()