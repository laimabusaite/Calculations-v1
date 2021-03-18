import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# theorPath = '../D1pumpBx/Pump_D1-Cs133_Probe_D1-Cs133/gamma=0.0190-DSteps=150-DScan=2.00sigma-lWidth=2.0/Pump_3-3-theta=1.57-phi=1.57-pol=0-det=0_Probe_3-3-det=0/Absorption_theta=1.57-0.79_phi=0.79-1.57_pol=0-0_gammaprobe=0.0190_lWidthpr=2.0/Absorption_Cs133-D1-PumpRabi=0.50-shift=0.0-Ge=0'
# theorPath = 'C:\\Users\\laima\\OneDrive - University of Latvia\\Pētījumi\\Atomi\\Cs magnetometrs\Calculations-v1\\D1pumpBx\\Pump_D1-Cs133_Probe_D1-Cs133\gamma=0.0190-DSteps=150-DScan=2.00sigma-lWidth=2.0\\Pump_3-3-theta=1.57-phi=1.57-pol=0-det=0_Probe_3-3-det=0\\Absorption_theta=1.57-0.79_phi=0.79-1.57_pol=0-0_gammaprobe=0.0190_lWidthpr=2.0\\Absorption_Cs133-D1-PumpRabi=0.50-shift=0.0-Ge=0'


def plotBy(rabiList, transitionList, probeList=None):
    rabiList = np.array(rabiList)
    transitionList = np.array(transitionList)
    if probeList == None:
        probeList = transitionList

    fig = plt.figure()
    fig.canvas.set_window_title('By')
    plt.suptitle('By')
    for idx, transition in enumerate(transitionList):
        plt.subplot(2,2,idx+1)
        plt.title(f'pump {transition}, probe {probeList[idx]}')
        for rabi in rabiList:

            # theorPath = '../D1pumpBx/Pump_D1-Cs133_Probe_D1-Cs133/Absorption_Cs133-D1-PumpRabi=0.50-shift=0.0-Ge=0'
            theorPath = f'./D1pumpBy/Pump_D1-Cs133_Probe_D1-Cs133/{transition}/Absorption_Cs133-D1-PumpRabi={rabi:.2f}-shift=0.0-Ge=0'#.format(rabi)


            theorData = pd.read_csv(theorPath, header=None, delimiter=' ', names=['B', 'Total', 'comp1', 'comp2', 'diff'])

            print(theorData)

            plt.plot(theorData['B'], theorData['comp1'], label=rabi)
            plt.legend()

if __name__ == '__main__':
    rabiList = [0.5, 1, 5, 10, 20, 50, 100]
    transitionList = ['3-3', '3-4', '4-4', '4-3']
    probeList = ['3-3', '3-3', '4-4', '4-4']
    plotBy(rabiList, transitionList)

    plt.show()