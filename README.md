# Fourier2D
Model-based image enhancement, inspired by an exoplanet imaging instrument called a coronagraph, that involves (inverse) 2D Fourier transforms.

# MATLAB program
After cloning the repository and making it the working folder, run the main script, 'coronagraph.m', a program that includes subfunctions. It requires the CC0-licensed image, 'Airy_disk_D65.png' [1] (and the helper function, 'dotdot.m'). Additional documentation about the program is available in Section III of a manuscript that will be available via arXiv.

# Example output
After time and space measurements complete for key steps in an image enhancement, the program outputs two graphs, 'secsvspels.pdf' and 'bytesvspels.pdf'. It then sets up and executes a model-based image enhancement, i.e., a phase aberration correction. Outputs include two images, 'groundtruth.jpg' and 'coronagraph.jpg', and a video, 'coronagraph.avi'. Examples were produced with a Dell Latitude E7450.

# See also
[1] Wikipedia contributors, “Airy disk,” Wikipedia, The Free Encyclopedia, Nov. 2023 (https://en.wikipedia.org/w/index.php?title=Airy_disk&oldid=1187106015)
