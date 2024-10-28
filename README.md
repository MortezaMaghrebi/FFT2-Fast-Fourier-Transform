**Title:** C# Image Processing Fourier Transform Library

**Description:**

This library provides efficient implementations of the Fast Fourier Transform (FFT) for image processing applications in C#. It enables you to analyze and manipulate images in the frequency domain, unlocking various image processing techniques like filtering, noise reduction, edge detection, and more.

## Fourier Class

This document describes the `Fourier` class, which provides functionalities for performing Fast Fourier Transforms (FFTs) on one and two dimensional data.

###  Features

* Complex number representation (`number` struct) for real and imaginary parts.
* Various Discrete Fourier Transform (DFT) functions:
    * `DFT1`: 1D DFT for a complex number array.
    * `iDFT1`: 1D Inverse DFT for a complex number array.
    * `DFT2`: 2D DFT for a complex number matrix.
    * `iDFT2`: 2D Inverse DFT for a complex number matrix.
* Fast Fourier Transform (FFT) functions (more efficient than DFT):
    * `FFT1`: 1D FFT for a complex number array.
    * `iFFT1`: 1D Inverse FFT for a complex number array.
    * `FFT2`: 2D FFT for a complex number matrix.
    * `iFFT2`: 2D Inverse FFT for a complex number matrix.
* Utility functions for complex number manipulation:
    * `Summation`: Adds two complex numbers.
    * `Subtraction`: Subtracts two complex numbers.
    * `Multiplication`: Multiplies two complex numbers.
    * `Division`: Divides two complex numbers.
    * `Conjucate`: Returns the complex conjugate of a number.
* Conversion functions for creating complex number representations from real number arrays and matrices.
    * `complex(double[])`: Creates a complex number array from a real number array.
    * `complex(double[,])`: Creates a complex number matrix from a real number matrix.
    * `complex(double[,], double[,])`: Creates a complex number matrix from separate real and imaginary parts (magnitude and phase).
    * `complex(Bitmap)`: Converts a grayscale image (Bitmap) to a complex number matrix.
* Get real, imaginary, abs and angle of the 1D and 2D signals.
    * `real(number[] array)` and `real(number[,] array)`:** Extracts the real part of each complex number.
    * `image(number[] array)` and `image(number[,] array)`:** Extracts the imaginary part of each complex number.
    * `abs(number[] array)` and `abs(number[,] array)`:** Calculates the magnitude (absolute value) of each complex number.
    * `angle(number[] array)` and `angle(number[,] array)`:** Calculates the phase angle of each complex number.


### Usage

The class provides various functions for performing different types of Fourier Transforms. Here's a general breakdown of their usage:

* **DFT/iDFT**:
    * Input: A complex number array/matrix (`number[]` or `number[,]`).
    * Output: The transformed complex number array/matrix.
    * Example: `number[] transformed = DFT1(data);`
* **FFT/iFFT**:
    * Input: A complex number array/matrix (`number[]` or `number[,]`).
    * Output: The transformed complex number array/matrix.
    * Example: `number[] transformed = FFT1(data);` (generally faster than DFT)
* **Complex number conversion**:
    * These functions convert real number representations to complex number representations for use with the transform functions.
    * Example: `number[] complex_data = complex(real_data);`

## ImageProcessingFourier Class

**Features:**

* **2D Fast Fourier Transform (FFT2):** Efficiently transforms an image into its frequency domain representation.
* **Inverse Fast Fourier Transform (IFFT2):** Reconstructs an image from its frequency domain representation.
* **Frequency Domain Filtering:** Applies filters to modify image properties based on frequency components.
* **Magnitude and Phase Extraction:** Separates the magnitude (intensity) and phase (orientation) information in the frequency domain.
* **Normalization:** Normalizes magnitude and phase images for visualization.

**Usage:**

1. **Reference the libraries:**

   ```csharp
   using ImageProcessingFourier;
   using Fourier;
   ```

2. **Create an instance of `ImageProcessingFourier.FFT2`:**

   ```csharp
   ImageProcessingFourier.FFT2 IPF = new ImageProcessingFourier.FFT2();
   ```

3. **Perform FFT on your image:**

   ```csharp
   Bitmap image = // Load your image using System.Drawing
   IPF.fft2(image);
   ```

4. **Apply filtering (optional):**

   ```csharp
   double[,] filter = // Define your filter (e.g., low-pass, high-pass)
   IPF.Filter(filter);
   ```

5. **Perform inverse FFT to reconstruct the image:**

   ```csharp
   Bitmap processedImage = IPF.ifft2();
   // Display or save the processed image
   ```
