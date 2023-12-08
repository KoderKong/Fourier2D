# Fourier2D
Model-based image enhancement, inspired by an exoplanet imaging instrument called a coronagraph, that involves 2D Fourier transforms.

# MATLAB program
After cloning the repository and making it the working folder, run the main script, 'coronagraph.m', a program that includes subfunctions. It requires the CC0-licensed image, 'Airy_disk_D65.png' [1] (and the helper function, 'dotdot.m'). Additional documentation about the program is available in Section III of a manuscript that is available via arXiv [2].

# Example output
After time and space measurements complete for key steps in an image enhancement, the program outputs two graphs, 'secsvspels.pdf' and 'bytesvspels.pdf'. It then sets up and executes a model-based image enhancement. Outputs include two images, 'groundtruth.jpg' and 'coronagraph.jpg', and a video, 'coronagraph.avi'. Examples were produced with a Dell Latitude E7450.

# See also
[1] Wikipedia contributors, "Airy disk," Wikipedia, The Free Encyclopedia, ID:1187106015, pp. 1-8, Nov. 2023 (https://en.wikipedia.org/w/index.php?title=Airy_disk&oldid=1187106015)
[2] Dileepan Joseph, "Ricci-Notation Tensor Framework for Model-Based Approaches to Imaging," arXiv:2312.04018v1 [cs.MS], pp. 1-37, Dec. 2023 (https://arxiv.org/abs/2312.04018)
