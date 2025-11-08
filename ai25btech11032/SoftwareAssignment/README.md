# Image Compression using Truncated SVD

This project demonstrates how Singular Value Decomposition (SVD) can be used to
compress grayscale images by keeping only the top `k` dominant singular values.
The idea is based on the fact that an image can be treated as a matrix, and most
of the important structure in the image is captured by the first few singular
components.

# Problem Overview
A grayscale image of size `m × n` can be represented as a matrix `A`.  
The SVD of `A` is:

\[
A = U \Sigma V^T
\]

To compress the image, we keep only the top `k` singular values:

\[
A_k = \sum_{i=1}^{k} \sigma_i u_i v_i^T
\]

This reduces storage and can still reconstruct an image that visually looks close
to the original for suitable `k`.

# Method Used
We do not use built-in SVD libraries.  
Instead, we compute the top `k` singular values and vectors using:

- **Power Iteration** to extract the dominant singular pair `(σ, u, v)`
- **Deflation** to remove this rank-1 component from the matrix before the next iteration

This allows us to get the top `k` modes without full SVD.

# Results
We tested on three images:
- `globe.jpg`
- `einstein.jpg`
- `greyscale.png`

# Observations
- Small `k` → blurry image, high error  
- Increasing `k` → sharper image, decreasing error  
- Images with **sharp edges** (e.g., Einstein) compress well at small `k`
- Images with **smooth gradients** (e.g., Globe) require larger `k`
