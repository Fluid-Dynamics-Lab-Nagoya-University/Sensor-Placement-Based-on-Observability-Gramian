# Sensor-Placement-Based-on-Observability-Gramian
This repository contains Matlab (R2022a) codes to reproduce results for "Efficient sensor node selection for observability Gramian optimization."

## License
[MIT-License](https://github.com/Aerodynamics-Lab/Sensor-Placement-Based-on-Observability-Gramian/blob/main/LICENSE)

## Code structures
---
### Main program  
- P1_main.mlx  

### Exporting graphics
- P2_plots.mlx

### Function  
#### Preprocessing  
- F_PIV_svd.m  

#### Sensor selection  
- F_sensor_random.m
- F_sensor_dg_gram.m  
- F_sensor_dg_gramgrad.m  
- F_sensor_CRSNC_gram.m  
- F_sensor_sdp_gram.m  

#### Calculation
- F_calc_obsgram.m  

      
## How to cite
If you use the codes in your work, please cite the software itself and relevent paper.
### General software reference:
```bibtex
@misc{yamada2023github,
      author = {Keigo Yamada},
      title = {Sensor Placement Based on Observability Gramian},
      howpublished = {Available online},
      year = {2023},
      url = {https://github.com/Aerodynamics-Lab/Sensor-Placement-Based-on-Observability-Gramian}
}
```
<!--
### Relevent paper reference:
```bibtex
@ARTICLE{yamada2021fast,
  author = {Keigo Yamada and Yuji Saito and Koki Nankai and Taku Nonomura and
	Keisuke Asai and Daisuke Tsubakino},
  title = {Fast greedy optimization of sensor selection in measurement with
	correlated noise},
  journal = {Mechanical Systems and Signal Processing},
  year = {2021},
  volume = {158},
  pages = {107619},
  abstract = {A greedy algorithm is proposed for sparse-sensor selection in reduced-order
	sensing that contains correlated noise in measurement. The sensor
	selection is carried out by maximizing the determinant of the Fisher
	information matrix in a Bayesian estimation operator. The Bayesian
	estimation with a covariance matrix of the measurement noise and
	a prior probability distribution of estimating parameters, which
	are given by the modal decomposition of high dimensional data, robustly
	works even in the presence of the correlated noise. After computational
	efficiency of the algorithm is improved by a low-rank approximation
	of the noise covariance matrix, the proposed algorithms are applied
	to various problems. The proposed method yields more accurate reconstruction
	than the previously presented method with the determinant-based greedy
	algorithm, with reasonable increase in computational time.},
  doi = {https://doi.org/10.1016/j.ymssp.2021.107619},
  issn = {0888-3270},
  keywords = {Data processing, Sensor placement optimization, Greedy algorithm,
	Bayesian state estimation},
  url = {https://www.sciencedirect.com/science/article/pii/S0888327021000145}
}

```
-->
## Author
YAMADA Keigo

[Experimental Aerodynamics Laboratory](http://www.aero.mech.tohoku.ac.jp/en/)
Department of Aerospace Engineering
Graduate School of Engineering

Tohoku University, Sendai, JAPAN

E-mail: keigo.yamada.t5@dc.tohoku.ac.jp

Github: [k5-Yamada](https://github.com/k5-Yamada)
