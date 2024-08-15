### BayesSpace Clustering


# This wont work for a SPATAData dataset so turn this chunk off if needed
# maybe set the wd to where the dataset is to ensure it is doing this for the right dataset?
# This takes around 11 mins to run.

# Ensure that the number of clusters is similar in both
# clustering methods used for an accurate comparison
# This works well but commenting it out now for fast execution

main_spata2_obj_hi_res <-
  runBayesSpaceClustering(
    object = main_spata2_obj_hi_res,
    name = "bayes_space" # the name of the output grouping variable
  )

# This will tell us if it worked or not
# a factor called "bayes_space" will have been created

getGroupingOptions(main_spata2_obj_hi_res)

plotSurface(
  object = main_spata2_obj_hi_res,
  color_by = "bayes_space",
  pt_size = 2,
  pt_clrp = "uc"
) +
  labs(subtitle = "Bayes Space Clusters")



### HW k-Means Clustering

num_clusters_bayes <- length(getGroupNames(main_spata2_obj_hi_res, grouping_variable = "bayes_space"))

print(paste0("There were ", num_clusters_bayes, " clusters in bayes_space clustering so the same number for k is assigned for k-means clustering so that the clusters are comparable"))

# The number of centers needs to be changed as required
# No other changes need to be made
# Code taken from the SPATA2 clustering tutorial
# https://themilolab.github.io/SPATA2/articles/spata-v2-clustering.html
# Double check it is doing the right thing

# Centers is what sets the value of k 

kmeans_res <- 
  stats::kmeans(
    x = getPcaMtr(main_spata2_obj_hi_res), 
    centers = num_clusters_bayes, 
    algorithm = "Hartigan-Wong"
  )

head(kmeans_res[["cluster"]])

cluster_df <- 
  as.data.frame(kmeans_res[["cluster"]]) %>% 
  tibble::rownames_to_column(var = "barcodes") %>% 
  magrittr::set_colnames(value = c("barcodes", "kmeans_4_HW")) %>% 
  tibble::as_tibble()

cluster_df[["kmeans_4_HW"]] <- as.factor(cluster_df[["kmeans_4_HW"]])

cluster_df

# This will tell us if it worked or not
# a factor called "kmeans_4_HW" will have been created

getGroupingOptions(main_spata2_obj_hi_res)

# overwrite = TRUE ensures that centers can be changed as needed 
# and code will still work when this chunk is executed again

main_spata2_obj_hi_res <- 
  addFeatures(
    object = main_spata2_obj_hi_res,
    feature_df = cluster_df,
    overwrite = TRUE
  )


# This will tell us if it worked or not
# a factor called "bayes_space will have been created

plotSurface(
  object = main_spata2_obj_hi_res,
  color_by = "kmeans_4_HW",
  pt_clrp = "uc",
  pt_size = 2,
)+
  labs(subtitle = "HW Kmeans Clusters")


