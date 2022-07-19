import _ctypes
import ctypes
import numpy as np
import matplotlib.pyplot as plt

# Load the library
lib = ctypes.cdll.LoadLibrary( "../DSAlgorithms/build/libDSAlgorithms.so" )

def tSNE_in_Cpp( data, dimReduced, perplexity, gdRate, gdMomentum ):
    # Array for the output of the computation in C++
    reducedData = np.zeros( (len(data), dimReduced), dtype=np.double )
    # Make sure that both arrays are contiguous, and that their data type is correct
    if not data.flags["C_CONTIGUOUS"]:
        data = np.ascontiguousarray( data )
    if not reducedData.flags["C_CONTIGUOUS"]:
        reducedData = np.ascontiguousarray( reducedData )
    if not data.dtype == 'float':
        data = data.astype( 'float' )
    if not reducedData.dtype == 'float':
        reducedData = reducedData.astype( 'float' )

    # Setup the types of the arguments of the function
    lib.computetSNE.argtypes = [ np.ctypeslib.ndpointer( dtype=ctypes.c_double, flags="C_CONTIGUOUS" ),
                                 ctypes.c_int, ctypes.c_int,
                                 np.ctypeslib.ndpointer( dtype=ctypes.c_double, flags="C_CONTIGUOUS" ),
                                 ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double ]
    # Setup the type of the return value of the function
    lib.computetSNE.restype = None # Void return type

    # Call the function
    # lib.computePCA( data, len(data), len(data[0]), reducedData, dimReduced )
    lib.computetSNE( data, len(data), len(data[0]), reducedData,
                     dimReduced, perplexity, gdRate, gdMomentum )

    # Return the result of the function
    return reducedData

# Load the dataset, with the labels

data = []
with open( 'Datasets/mnist500_X.txt', 'r' ) as file:
    for line in file:
        line = line.strip().split( "   " )
        data.append( [ float(nr) for nr in line ] )

labels = []
with open( 'Datasets/mnist500_labels.txt', 'r' ) as file:
    for line in file:
        line = line.strip().split( "   " )[0]
        labels.append( int(float(line) ) )

data = np.array( data )
labels = np.array( labels )

# Take only a few of the data points
nPoints = 100
data = data[ : nPoints ]
labels = labels[ : nPoints ]

# Parameters for t-SNE
nr_dimensions = 2
perplexity = 30.
gradientRate = 0.5
gradientMomentum = 0.5

# Perform t-SNE on the data
reducedData = tSNE_in_Cpp( data, nr_dimensions, perplexity, gradientRate, gradientMomentum )

# print(reducedData)

# Plot the results
x = reducedData[:, 0]
y = reducedData[:, 1]
print(x)
print(y)

print(len(x))
print(len(y))

# plt.scatter( x, y )
# plt.show()

# # Plot the results
# # fig1 = plt.figure( figsize=(10, 10), dpi=100 )
# fig1 = plt.figure()
# ax1 = fig1.add_subplot(111)
# ax1.set_title( 'Events projected onto a 2-dimensional space using t-SNE' )
# # Plot the points corresponding to each class with a different color
# for cl in range( max(labels) + 1 ):
#     mask = labels == cl
#     # print(cl)
#     # print(reducedData[mask,0])
#     # print(reducedData[mask,1])
#     # print()
#     # print()
#     # print([reducedData[mask, 0], reducedData[mask,1]])
#     # ax1.scatter( reducedData[mask, 0], reducedData[mask, 1], label="Class {}".format(cl) )
#     ax1.scatter( reducedData[mask, 0], reducedData[mask, 1], label='Class {}'.format(cl) )
#     # ax1.scatter( reducedData[mask, 0], reducedData[mask, 1] )
#
# # ax1.scatter(reducedData[:, 0], reducedData[:, 1])
# ax1.set_xlabel(r'$U_1$')
# ax1.set_ylabel(r'$U_2$')
# ax1.legend(loc='best')
#
# plt.show()
