## The optimal neural network set up for denoising

## Denoising (Assessment)

# # Run assessment to determine parameters
#
main_spata2_obj_hi_res <-
  runAutoencoderAssessment(
    object = main_spata2_obj_hi_res,
    activations = c("relu", "selu", "sigmoid"),
    bottlenecks = c(32, 40, 48, 56, 64),
    epochs = 20,
    layers = c(128, 64, 32),
    dropout = 0.1
  )


### Continued in the main script