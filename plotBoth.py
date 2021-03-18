import matplotlib.pyplot as plt
from plotBx import plotBx
from plotBy import plotBy
from plotTheory import plotRabiDependence

if __name__ == "__main__":

    rabiList = [0.5, 1, 5]#, 10, 20, 50, 100]
    # transitionList = ['3-2', '3-3', '3-4', '4-3', '4-4', '4-5']
    # probeList = transitionList #['3-3', '3-4', '4-4', '4-3']
    # # plotBx(rabiList, transitionList,probeList)
    # # plotBy(rabiList, transitionList, probeList)
    # plotRabiDependence('Bx', rabiList, transitionList, probeList, dline=2, x_range = 10)
    # plotRabiDependence('By', rabiList, transitionList, probeList, dline=2)

    transitionList = ['3-3', '3-4', '4-3', '4-4']
    probeList = transitionList  # ['3-3', '3-4', '4-4', '4-3']
    # plotBx(rabiList, transitionList,probeList)
    # plotBy(rabiList, transitionList, probeList)
    plotRabiDependence('Bx', rabiList, transitionList, probeList, dline=1, x_range = 10)
    # plotRabiDependence('By', rabiList, transitionList, probeList, dline=2)

    plt.show()