# Honor_Thesis
Includes some core files of the research project

[Abstract] The Heston stochastic volatility model is one of the most widely studied stochastic volatility models, in which the variance follows a Cox–Ingersoll–Ross process. Estimating this model under the physical measure is challenging, as the likelihood function involves high-dimensional integral. While an approximate analytical solution for the likelihood function exists, the task of maximizing the function remains difficult in practice. Furthermore, these approximate solutions are invalid if any modifications or extensions of the Heston model are considered, such as extending the model to higher dimensions. Being full-information, plug-and-play, and frequentist, iterated filtering algorithms are adopted to estimate the volatility process of the Heston model. We use the S&P500 index as an example, estimating model parameters and their confidence intervals. The results indicate that the estimated volatility of the S&P500 index matches the pattern of the VIX index. An application in options pricing is also given. We then demonstrate the benefit of iterated-filtering methods by extending to a multi-dimensional panel of Heston models, estimating the volatility processes of four emerging market indices. The results illustrate that the volatility processes of these emerging market indices share the same sensitivity to their corresponding price processes with a 95% confidence level.

[Files] 
DATA:
(i) SPX.csv: Data we used for 1-d model (SPX from yahoo finance)
(ii) brazil.csv: Data we used in the Panel POMP model (BVSP from yahoo finance)
(iii) india.csv: Data we used in the Panel POMP model (BSESN from Bombay Exchange)
(iv) indonesia.cvs: Data we used in the Panel POMP model (JSKE from yahoo finance)
(v) mexico.csv: Data we used in the Panel POMP model (MXX from yahoo finance)

CODE:
(vi) 1d_global_search.R: 1d POMP model (profile not included)
(vii) panel_global_search.R: Panel POMP model (profile not included)
(viii) Estimate_vol_process.R: Estimate volatility process
(ix) shared1_figs.R: generate profile figure
(x) getBestMods.R: get parameters estimations based on our modified profile dataset
(xi) sde_sim.ipynb: jupiter notebook (python) for heston model simulation using Euler scheme