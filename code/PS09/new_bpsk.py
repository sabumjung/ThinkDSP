import scipy
from scipy.io import wavfile
import numpy as np
import matplotlib.pyplot as plt
import thinkdsp
from array import array
import bpsk as original_bpsk


# This function converts a string into a numpy array of bits
# note that it is assumed that each character is 7 bits long here
def string2NPArray(s):
    bits = np.array([])
    for a in bytearray(s, 'ascii'):
        for b in range(0,7):
            bits = np.append(bits,float((a>>(7-b-1))&1))
    return bits

# This function converts a numpy array of bits to a string
# note that it is assumed that each character is 7 bits long here
def NPbits2String(bits):
    S = ""
    for a in np.arange(0, len(bits)/7):
        tmp = 0
        for k in np.arange(0,7):
            b = bits[a*7+k]
            tmp = tmp + (2**(6-k))*b
        S = S + chr(int(tmp))
    return S

# this function is used to help convert numpy array data into a format
# suitable for writing into a wave file
def convert_to_int16(sig):
    # convert into int16  to write as wave
    sig = (sig/np.max(sig))*(2**14)
    sig = sig.astype('int16')
    return sig
    

# this is a utility function that  finds the start and  end 
# of transmission in the numpy array of samples xf
# The function looks for the first instance where the entries of xf
# go above threshold and returns the index into xf where this happens
# in start_idx
# The function looks for the last instance where the entries of xf
# are above threshold and returns the index into xf where this happens
# in end_idx
# 
# You will probably have to do some trial and error to get the threshold right
# one possibility is to se the threshold equal to some factor of the maximum value
# in the input signal,  e.g. 0.3 * maximum value in xf
#
def find_start_and_end(xf, threshold = 2000): 
    import numpy as np    
    start_idx = -1
 
    for k in range(0, len(xf)):
        if(np.abs(xf[k])) > threshold:
            start_idx = k
            break

    if(start_idx  < 0):
        print "Unable to detect start of transmission"
        return -1
    
    end_idx = -1
    
    for k in range(0, len(xf)):
        if(np.abs(xf[len(xf)-k-1])) > threshold:
            end_idx = len(xf)-k-1
            break

    if(end_idx < 0):
        print "Unable to detect end of transmission"
        return -1

    return start_idx, end_idx

def my_gen_bpsk_1(bits, sample_rate=8000, symbol_period=.01, freq=1000):
    symbol_len = symbol_period * sample_rate
    square_signal = bits_to_signal(bits, symbol_len)
    transmission_duration = len(square_signal)*symbol_period
    cos_signal = thinkdsp.CosSignal(freq=freq, offset=0)
    print transmission_duration
    print sample_rate
    cos_wave = cos_signal.make_wave(duration = transmission_duration, framerate = sample_rate)
    print cos_wave
    print square_signal
    print len(cos_wave.ys)
    print len(square_signal)
    return np.multiply(cos_wave.ys, square_signal)

def plot_bpsk():
    my_bits = string2NPArray('ab')
    sig,cos,ts,square_signal = my_gen_bpsk(my_bits)
    plt.subplot(3,1,1)
    plt.plot(ts,square_signal)
    plt.subplot(3,1,2)
    plt.plot(ts,cos)
    plt.subplot(3,1,3)
    plt.plot(ts,sig)
    plt.show()

def bits_to_signal(bits, symbol_len):
    """converts array of bits to a signal array, where -1 and 1 represent 0 and 1, and
    control and transmission start/end bits are put in"""
    symbol_len = int(symbol_len)
    symbol_zeros = np.array([0]*symbol_len)
    bits = [-1 if x==0 else x for x in bits]
    square_signal = np.array([1]*250) #first bit transmitted is always one to help with synchronization
    for bit in bits:
        square_signal = np.hstack((square_signal, np.array([bit]*symbol_len)))
    #append begin/end transmission bits
    square_signal = np.hstack((symbol_zeros, square_signal, symbol_zeros))
    return square_signal

def my_gen_bpsk(bits, rate=8000, symbol_len=250, freq=1000):
    square_signal = bits_to_signal(bits, symbol_len)
    transmission_end_t = (float(symbol_len)/rate)*len(square_signal)
    times = np.linspace(0,transmission_end_t,num=len(square_signal))
    cos_signal = np.cos(freq*times)
    bpsk_signal = np.multiply(cos_signal, square_signal)
    return bpsk_signal, cos_signal, times, square_signal

def generate_bpsk_signal(bits, rate = 8820, symbol_len = 250, freq = 1000):
    square_signal = bits_to_signal(bits, symbol_len)
    transmission_end_t = len(square_signal)/(1.0*rate)*2*np.pi
    times = np.linspace(0,transmission_end_t,num=len(square_signal))
    cos_signal = 7500*np.cos(freq*times)
    bpsk_signal = np.multiply(cos_signal, square_signal)
    return bpsk_signal, cos_signal

def low_pass_filter(signal, rate, cutoff):
    w = thinkdsp.Wave(signal, rate)
    s = w.make_spectrum()
    s.low_pass(cutoff)
    return s.make_wave().ys

def find_phase(fft):
    max_component = max(fft, key=lambda f: abs(f))
    return np.angle(max_component)

def demodulate(signal, samples_per_fft):
    #note: it is known that this is not actually the demodulate step, demodulation is the multiplication by carrier frequency cosine, 
    #but we don't really want to mess with the namespace now
    phases = np.array([])
    for i in np.arange(0,len(signal),samples_per_fft):
        try:
            chunk = signal[i:i+samples_per_fft]
        except:
            chunk = signal[i:]
        fft = scipy.fftpack.fft(chunk)
        phases = np.hstack((phases, find_phase(fft)))
    return phases

def plot_demodulate(signal, samples_per_fft):
    phases = np.array([])
    for i in np.arange(0,len(signal),samples_per_fft):
        try:
            chunk = signal[i:i+samples_per_fft]
        except:
            chunk = signal[i:]
        fft = scipy.fftpack.fft(chunk)
        phases = np.hstack((phases, find_phase(fft)))
    plt.plot(phases)
    plt.show()
    
def decode_phases(phases):
    normalized = phases - np.mean(phases)
    decoded_phases = [1 if p > 0 else 0 for p in normalized]
    if decoded_phases[0] == 1:
        #flip the bits if the 0th bit is 1 instead of zero.
        decoded_phases = [1 if p == 0 else 0 for p in decoded_phases]
    return decoded_phases

def cut_signal(x,detection_threshold_factor):
    """snips out relevant part of signal"""
    max_val = np.amax(np.abs(x))
    print 'max val signal'
    print max_val
    beginning, end = find_start_and_end(x, detection_threshold_factor*max_val)
    return x[beginning:end],beginning,end

def decode_bits(bits):
    """chops off transmission start/end bits and flips bits as nesseccesary"""
    #bits = bits[1:-1]
    if bits[1] == 0:
        bits = [1 if b == 0 else 0 for b in bits]
    return bits[1:]

def decode_bpsk_signal(x, freq=1000, rate = 8000, symbol_len = 250, detection_threshold_factor = 0.3, LPFbw = 320):
    x,beginning,end = cut_signal(x,detection_threshold_factor)
    transmission_end_t = len(x)/(1.0*rate)*2*np.pi
    times = np.linspace(0,transmission_end_t,num=len(x))
    phase_offsets = np.linspace(0,2*np.pi,16)
    filtered_signals = []
    for phase_offset in phase_offsets:
        cos_signal = np.cos(freq*times + phase_offset)
        recovered_signal = np.multiply(x, cos_signal)
        filtered_signal = low_pass_filter(recovered_signal, rate, LPFbw)
        filtered_signals.append(filtered_signal)
    best_signal = max(filtered_signals, key=lambda signal: max(signal))
    best_signal,beginning,end = cut_signal(best_signal,detection_threshold_factor)
    phases = demodulate(best_signal, samples_per_fft=symbol_len)
    coded_bits = decode_phases(phases)
    bits = decode_bits(coded_bits)
    return bits

def plot_decode_bpsk_signal(x, freq=1000, rate = 8000, symbol_len = 250, detection_threshold_factor = 0.3, LPFbw = 320):
    x,beginning,end = cut_signal(x,detection_threshold_factor)
    transmission_end_t = len(x)/(1.0*rate)*2*np.pi
    times = np.linspace(0,transmission_end_t,num=len(x))
    phase_offsets = np.linspace(0,2*np.pi,16)
    filtered_signals = []
    for phase_offset in phase_offsets:
        cos_signal = np.cos(freq*times + phase_offset)
        recovered_signal = np.multiply(x, cos_signal)
        filtered_signal = low_pass_filter(recovered_signal, rate, LPFbw)
        filtered_signals.append(filtered_signal)
    best_signal = max(filtered_signals, key=lambda signal: max(signal))
    best_signal,beginning,end = cut_signal(best_signal,detection_threshold_factor)
    plt.plot(best_signal)
    phases = plot_demodulate(best_signal, samples_per_fft=symbol_len)
    return phases
    coded_bits = decode_phases(phases)
    bits = decode_bits(coded_bits)
    return bits

def test_decode():
    bits = string2NPArray('a')
    bpsk_signal, cos_signal, times, square_signal = my_gen_bpsk(bits)
    phases, decoded = decode_bpsk_signal(bpsk_signal)
    plt.subplot(3,1,1)
    plt.plot(square_signal)
    plt.axis()
    plt.subplot(3,1,2)
    plt.plot(phases)
    plt.subplot(3,1,3)
    plt.plot(decoded)
    plt.show()

def plot_wav_file(file_name):
    fs, x = wavfile.read(file_name)
    plt.plot(x)
    plt.show()

if __name__ == '__main__':
    fs, x = wavfile.read('AcousticModemTx.wav')
    signal,beginning,end = cut_signal(x, detection_threshold_factor=1.0)
    print beginning,end
    plt.subplot(2,1,1)
    plt.plot(x)
    plt.subplot(2,1,2)
    plt.plot(x[beginning:end])
    plt.show()
    #bits,beginning,end = decode_bpsk_signal(x, freq= 1000, rate = fs, symbol_len = 250, detection_threshold_factor = 0.4, LPFbw = 320)
    #print bits
    #message_string = NPbits2String(bits)
    #print message_string
    #plot_wav_file('AcousticModemTx.wav')
    #bits = string2NPArray('autist')
    #signal,cos = generate_bpsk_signal(bits, rate = 8820, symbol_len = 250, freq = 1000)
    #wavfile.write('AcousticModemRx.wav', 8820, convert_to_int16(signal))
    #original_signal = original_bpsk.generate_bpsk_signal(bits, rate = 8820, symbol_len = 250, freq = 1000)
    #wavfile.write('OriginalAcousticModemRx.wav', 8820, convert_to_int16(signal))

    #plt.subplot(3,1,1)
    #plt.plot(signal[:2000])
    #plt.subplot(3,1,2)
    #plt.plot(original_signal[:2000])
    #plt.subplot(3,1,3)
    #plt.plot(cos)
    #plt.show()
