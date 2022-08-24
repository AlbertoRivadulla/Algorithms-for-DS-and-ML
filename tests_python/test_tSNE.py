import _ctypes
import ctypes
import numpy as np
import matplotlib.pyplot as plt

def tSNE_in_Cpp( dataIn, dimReduced, perpl, gdRate, gdMomentum ):
    # Load the library
    # lib = ctypes.cdll.LoadLibrary( "../DSAlgorithms/build/libDSAlgorithms.so" )
    lib = ctypes.cdll.LoadLibrary( "../DSAlgorithms/build/libDSAlgorithms.dylib" )
    # Array for the output of the computation in C++
    reducedDataOut = np.zeros( (len(dataIn), dimReduced), dtype=np.double )
    # Make sure that both arrays are contiguous, and that their data type is correct
    if not dataIn.flags["C_CONTIGUOUS"]:
        dataIn = np.ascontiguousarray( dataIn )
    if not reducedDataOut.flags["C_CONTIGUOUS"]:
        reducedDataOut = np.ascontiguousarray( reducedDataOut )
    # if not dataIn.dtype == 'float':
    #     dataIn = dataIn.astype( 'float' )
    # if not reducedDataOut.dtype == 'float':
    #     reducedDataOut = reducedDataOut.astype( 'float' )

    # Setup the types of the arguments of the function
    lib.computetSNE.argtypes = [ np.ctypeslib.ndpointer( dtype=ctypes.c_double, flags="C_CONTIGUOUS" ),
                                 ctypes.c_int, ctypes.c_int,
                                 np.ctypeslib.ndpointer( dtype=ctypes.c_double, flags="C_CONTIGUOUS" ),
                                 ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double ]
    # Setup the type of the return value of the function
    lib.computetSNE.restype = None # Void return type

    # Call the function
    lib.computetSNE( dataIn, len(dataIn), len(dataIn[0]), reducedDataOut,
                     dimReduced, perpl, gdRate, gdMomentum )

    # Return the result of the function
    return reducedDataOut

# Load the dataset, with the labels
data = []
# with open( 'Datasets/mnist500_X.txt', 'r' ) as file:
with open( 'Datasets/mnist2500_X.txt', 'r' ) as file:
    for line in file:
        line = line.strip().split( "   " )
        data.append( [ float(nr) for nr in line ] )
labels = []
# with open( 'Datasets/mnist500_labels.txt', 'r' ) as file:
with open( 'Datasets/mnist2500_labels.txt', 'r' ) as file:
    for line in file:
        line = line.strip().split( "   " )[0]
        labels.append( int(float(line) ) )
data = np.array( data, dtype=np.double )
labels = np.array( labels )

# Take only a few of the data points
nPoints = 2500
data = data[ : nPoints ]
labels = labels[ : nPoints ]

# # Load the dataset
# file = open('Datasets/test_temperature.dat')
# data = []
# for line in file:
#     data.append( [] )
#     line_aux = line.split( ' ' )
#     for i in range( len(line_aux) ):
#         # Remove the first (city name) and last (direction) elements of each line
#         # if i != 0 and i < len(line_aux) - 1:
#         if i != 0 and i < len(line_aux) - 5:
#             data[-1].append( float(line_aux[i]) )
# data = np.array( data )
# file.close()

# Parameters for t-SNE
nr_dimensions    = 2
perplexity       = 30.
gradientRate     = 500
gradientMomentum = 0.8

# Perform t-SNE on the data
reducedData = tSNE_in_Cpp( data, nr_dimensions, perplexity, gradientRate, gradientMomentum )

# Plot the results
fig1 = plt.figure( figsize=(10, 10), dpi=100 )
# fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.set_title( 'Events projected onto a 2-dimensional space using t-SNE' )
# Plot the points corresponding to each class with a different color
for cl in range( max(labels) + 1 ):
    mask = labels == cl
    ax1.scatter( reducedData[mask, 0], reducedData[mask, 1], label='Class {}'.format(cl) )

ax1.set_xlabel(r'$U_1$')
ax1.set_ylabel(r'$U_2$')
ax1.legend(loc='best')

plt.show()
