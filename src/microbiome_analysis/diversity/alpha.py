# diversity/alpha.py
from skbio.diversity import alpha_diversity

def alpha_table(counts_df):
    return pd.DataFrame({
        "richness": alpha_diversity("observed_otus", counts_df.values, ids=counts_df.index),
        "shannon":  alpha_diversity("shannon",        counts_df.values, ids=counts_df.index),
        "simpson":  alpha_diversity("simpson",        counts_df.values, ids=counts_df.index),
    })

# diversity/beta.py
from skbio.diversity import beta_diversity
from skbio.stats.distance import permanova, DistanceMatrix

def beta_and_test(counts_df, metadata, group_col="body_site"):
    dm = DistanceMatrix(beta_diversity("braycurtis", counts_df.values,
                                       ids=counts_df.index), ids=counts_df.index)
    res = permanova(dm, metadata.loc[dm.ids, group_col].values, permutations=999)
    return dm, res  # res has pseudo-F, p-value, R²

