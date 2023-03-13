import numpy as np
import sys
if 'module' not in sys.path:
    sys.path.append('module')

from dtw_semi_global import semi_global_dtw
from dtw_semi_global import semi_global_dtw_with_rescaling


def drawFigureForTest(position1, position2):

    signal1 = get_signal_file(filetxt_path = 'testData/signal_0.txt')[position1:position2]

    signal2 = get_signal_file(filetxt_path = 'testData/top_adapter_signal.txt')

    plt.plot(normalise(signal1))
    plt.plot(normalise(signal2))

    plt.savefig('testData/test.png')
    plt.show()

def drawFigureForCheckbarcodeSig(position1):

    signal1 = get_signal_file(filetxt_path = 'testData/signal_0.txt')[position1:position1 + 250]

    signal2 = get_signal_file(filetxt_path = '/data/qijunhai/barcodeDesignData/test_new15bp_noise/BarcodeSignal15bp831/timeSeries_0.txt')

    plt.plot(normalise(signal1))
    plt.plot(normalise(signal2))

    plt.savefig('testData/test.png')
    plt.show()


def normalise(signal):
    if len(signal) == 0:
        return signal
    mean = np.mean(signal)
    stdev = np.std(signal)
    if stdev > 0.0:
        return (signal - mean) / stdev
    else:
        return signal - mean

def get_signal_file(filetxt_path):
    signal_list = list()
    with open(filetxt_path, 'r') as f:
        signal_file = f.readlines()
        for line in signal_file:
            signal_Value = float(line.rstrip())
            signal_list.append(signal_Value)

    signal = np.array(signal_list)
    return normalise(signal)
    # return signal

def fromLongRefFindShortQuery(ref_signal, query_signal):

    # query_signal = get_signal_file(query_signal_path)
    # ref_signal = get_signal_file(ref_siganl_path)

    normalise_query_signal = normalise(query_signal)
    normalise_ref_signal = normalise(ref_signal[0:2000])

    # topRef_signal = normalise_ref_signal[0:1000]

    position_start, position_end = \
        semi_global_dtw(normalise_ref_signal, normalise_query_signal)

    return (position_start, position_end + 10)

def fromLongRefFindShortQuery_test(ref_signal_path, query_signal_path):

    query_signal = get_signal_file(query_signal_path)
    ref_signal = get_signal_file(ref_signal_path)

    normalise_query_signal = normalise(query_signal)
    normalise_ref_signal = normalise(ref_signal[0:2000])

    # topRef_signal = normalise_ref_signal[0:1000]

    position_start, position_end = \
        semi_global_dtw(normalise_ref_signal, normalise_query_signal)

    return (position_start, position_end + 10)


def RefromLongRefFindShortQuery(ref_signal, query_signal):

    # query_signal = get_signal_file(query_signal_path)
    # ref_signal = get_signal_file(ref_siganl_path)

    normalise_query_signal = normalise(query_signal)
    normalise_ref_signal = normalise(ref_signal)

    # topRef_signal = normalise_ref_signal[0:1000]

    position_start, position_end = \
        semi_global_dtw_with_rescaling(normalise_ref_signal, normalise_query_signal)

    return (position_start, position_end + 10)

def RefromLongRefFindShortQuery_twoEnd(ref_signal_path, bottom_signal_path, topAndAdapter_signal_path):
    query_signal1 = get_signal_file(topAndAdapter_signal_path)
    # print(len(query_signal1))
    query_signal2 = get_signal_file(bottom_signal_path)

    ref_signal = get_signal_file(ref_signal_path)

    normalise_query_signal1 = normalise(query_signal1)
    normalise_query_signal2 = normalise(query_signal2)
    normalise_ref_signal = normalise(ref_signal[0:2000])

    # topRef_signal = normalise_ref_signal[0:1000]

    position_start1, position_end1 = \
        semi_global_dtw(normalise_ref_signal, normalise_query_signal1)

    position_start2, position_end2 = \
        semi_global_dtw(normalise_ref_signal, normalise_query_signal2)

    # print(position_start1, position_end1, position_start2, position_end2)
    return (position_end1, position_start2, ref_signal[position_end1 + 60:position_start2-40])

if __name__ == "__main__":
    pass


