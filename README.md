# image-matrix-completion
Open-sourced project for CS 5220

This project takes a "noisy" image, and attempts to recover the original using a low-rank matrix approxipation.
This low-rank approximation is generated via ALS matrix completion. See the report in the file report.pdf for a more technical explanation on the topic.


This project can be compiled with the command `gcc -O3 matrix.c png.c -lpng -llapacke -lcblas -lomp -fopenmp -o image_completion`

![Masked/Noisy version of Bella the dog](bella_masked.png)

![Recovered version of Bella the dog](bella_recovered.png)
