#This README explains exactly how to **compile and run** the code. 

# Requirements You only need: - **GCC compiler** and Two single-header image libraries: - stb_image.h - stb_image_write.h Download them directly in the same folder where your main.c file is located:

#commands to download:
wget https://raw.githubusercontent.com/nothings/stb/master/stb_image.h
wget https://raw.githubusercontent.com/nothings/stb/master/stb_image_write.h


Compilation
Go to the folder that contains your code and in terminal, run:

gcc main.c -o svd -lm
-lm links the math library required by functions.

This will create an executable named svd.

Running the Program
Start the program:

./svd
You will be required to enter:

Enter input file name:
Enter output file name:
Enter rank k:

The program will:

Read the input image,
Compute the top-k SVD,
Reconstruct the compressed image,
Save the output as a JPG file.

Choosing the Rank k
Small k leads to more compression but image becomes blurry

Large k leads to better quality but less compression

Output Information
The program prints:

Minimum, maximum, mean pixel value of reconstructed image
Absolute Frobenius error ||A − Ak||
Relative Frobenius error ||A − Ak|| / ||A||
Compression ratio estimate

The reconstructed image will be found in the figs folder.


