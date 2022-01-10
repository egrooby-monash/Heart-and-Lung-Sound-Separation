# Heart-and-Lung-Sound-Separation

Primary aim of the code is to denoise, and separate noisy neonatal chest sounds into heart and lung sounds. The MATLAB script “example_code.m” shows how to implement existing denoising and sound separation methods presented in the folder “Past Work Code” and the proposed non-negative matrix factorisation (NMF) and non-negative matrix co-factorisation (NMCF) methods presented in the folder “NMF and NMCF”. The folder “Signal Quality Evaluation” includes the code for the bss_eval toolbox which was used in the paper to assess the effectiveness of the sound separation methods. Additionally, in our previous work we have developed our own heart and lung sound signal quality evaluation, which can be found here https://github.com/egrooby-monash/Heart-and-Lung-Signal-Quality-Estimation. 

Note: This code is dependent on the code previously developed here https://github.com/egrooby-monash/Heart-and-Lung-Signal-Quality-Estimation

Please cite “Grooby, Ethan, et al. "A New Non-Negative Matrix Co-Factorisation Approach for Noisy Neonatal Chest Sound Separation." 2021 43rd Annual International Conference of the IEEE Engineering in Medicine & Biology Society (EMBC). IEEE, 2021.” when using this code. 

Additionally, you may also want to cite "Grooby, Ethan, et al. "Neonatal heart and lung sound quality assessment for robust heart and breathing rate estimation for telehealth applications." IEEE Journal of Biomedical and Health Informatics (2020)." and "Grooby, Ethan, et al. "Real-Time Multi-Level Neonatal Heart and Lung Sound Quality Assessment for Telehealth Applications." arXiv preprint arXiv:2109.15127 (2021).". 
