from matplotlib import pyplot as plt
import argparse

def text2sig(sigFile):

    file = open(sigFile)
    sigList = []
    for sig in file:
        sig = sig.strip('\n')
        if sig != '':
            sigList.append(float(sig))
    file.close()

    return sigList


def drawSigPlot(sigFile, outFig):
    sig = text2sig(sigFile)
    plt.plot(sig)
    plt.savefig(outFig)

def initializationParameters():
    parser = argparse.ArgumentParser()

    parser.add_argument('--sig-txt', type = str, required = True,
                        help = 'Indicates a input txt file that containing a nanopore signal.')
    parser.add_argument('--out-fig', type = str, required = True,
                        help = 'Indicates a output figure file corresponding to the input nanopore signal.')

    args = parser.parse_args()
    return args

def main():
    args = initializationParameters()
    drawSigPlot(sigFile = args.sig_txt, outFig = args.out_fig)

if __name__ == "__main__":
    main()