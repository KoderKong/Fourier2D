# Fourier2D
Model-based image enhancement, inspired by an exoplanet imaging instrument called a coronagraph, that involves 2D Fourier transforms.

# MATLAB program
After cloning the repository and making it the working folder, run the main script, 'coronagraph.m'. It requires an included CC0-licensed image, 'Airy_disk_D65.png' [1], as well as MATLAB's Computer Vision and Optimization Toolboxes. Additional documentation about the program is available in Section III of an upcoming publication [2].

# Example output
After time and space measurements complete for key steps in an image enhancement, the program outputs two graphs, 'secsvspels.pdf' and 'bytesvspels.pdf'. It then sets up and executes a model-based image enhancement. Outputs include two images, 'groundtruth.png' and 'coronagraph.png', and a video, 'coronagraph.avi'. Examples were produced with a Dell Latitude E7450.

# See also
[1] Wikipedia contributors, "Airy disk," Wikipedia, The Free Encyclopedia, pp. 1-8, Nov. 2023 (https://en.wikipedia.org/w/index.php?title=Airy_disk&oldid=1187106015)

[2] Dileepan Joseph, "Ricci-Notation Tensor Framework for Model-based Approaches to Imaging," Journal of Imaging Science and Technology, vol. 68, no. 4, pp. 040504 1-17, Jul. 2024 (https://doi.org/10.2352/J.ImagingSci.Technol.2024.68.4.040504, https://arxiv.org/abs/2312.04018)
