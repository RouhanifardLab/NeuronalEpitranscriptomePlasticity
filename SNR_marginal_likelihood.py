import pandas as pd
from scipy.special import loggamma
import sys
import numpy as np

def log_beta_prior():
    return loggamma(0.5) + loggamma(0.5) - loggamma(1)

def log_beta_posterior(n, h):
    return loggamma(h + 0.5) + loggamma(n - h + 0.5) - loggamma(n + 0.5)

def log_marginal_likelihood(n, h):
    return log_beta_posterior(n, h) - log_beta_prior()

def log_marginal_likelihood_ratio(n1, h1, n2, h2):
    return log_marginal_likelihood(n1, h1) + log_marginal_likelihood(n2, h2) - log_marginal_likelihood(n1 + n2, h1 + h2)

def main(in_csv, out_csv):
    
    df = pd.read_csv(in_csv)
    
    # Extracting column indices from the DataFrame
    drs_u = df.columns.get_loc("DRS_U")
    drs_c = df.columns.get_loc("DRS_C")
    ivt_u = df.columns.get_loc("IVT_U")
    ivt_c = df.columns.get_loc("IVT_C")
    

    snr_array = [log_marginal_likelihood_ratio(drs_u + drs_c, 
                                               drs_c, 
                                               ivt_u + ivt_c, 
                                               ivt_c) for drs_u, drs_c, ivt_u, ivt_c in zip(df.iloc[:, drs_u].to_list(),
                                                                                            df.iloc[:, drs_c].to_list(),
                                                                                            df.iloc[:, ivt_u].to_list(),
                                                                                            df.iloc[:, ivt_c].to_list())]
    print(snr_array)
    df["SNR"] = snr_array
    df.to_csv(out_csv, index=False)

    

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
    