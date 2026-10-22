# Instructions for downloading necessary packages

Available at Medical Physics: [https://doi.org/10.1002/mp.70670](https://doi.org/10.1002/mp.70670)
Also available as a pre-print on medRxiv: [https://doi.org/10.64898/2025.12.03.25341501](https://doi.org/10.64898/2025.12.03.25341501)

Note, this code has been developed on MATLAB version R2023a. It is untested how this code will behave in other versions of MATLAB.

1. Download [fiNUFFT](https://finufft.readthedocs.io/en/latest/matlab.html) and place the contents into `./finufft/matlab/` ... at minimum you need `finufft_plan.m` for reference
2. Download the contents of [Mark Bydder's PRESCO package](https://github.com/marcsous/pdff.git) and place them into `./contrib/presco/pdff-master/`
3. Download the contents of [Mark Bydder's unwrapping package](https://github.com/marcsous/unwrap.git) and place them into `./contrib/presco/unwrap-master/`
4. Download the `lipari.mat` and `navia.mat` files from [Fabio Crameri](https://www.fabiocrameri.ch/colourmaps/) and place them into `./colormaps/`
5. If desired, download the sample data from [FigShare](https://doi.org/10.6084/m9.figshare.30779273) and place the contents into `./data_in/`.
