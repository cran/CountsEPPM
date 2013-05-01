Model.Counts <-
function(parameter,model.type,model,covariates.matrix.mean,
                   covariates.matrix.variance,offset.mean,offset.variance,
                   scale.factor.model,vnmax) {
   if (model.type=="mean only") { 
                   output <- Model.Faddy(parameter,model.type,model,covariates.matrix.mean,
                                       offset.mean,vnmax)
                                } # end of over-dispersed models
   if (model.type=="mean and variance") { 
          if (model=="general") { 
                   output <- Model.FaddyJMV.general(parameter,covariates.matrix.mean,
                                         covariates.matrix.variance,offset.mean,
                                         offset.variance,scale.factor.model,vnmax) }
          if (model=="limiting") { 
                   output <- Model.FaddyJMV.limiting(parameter,covariates.matrix.mean,
                                         covariates.matrix.variance,offset.mean,
                                         offset.variance,scale.factor.model,vnmax) }
                                } # end of under-dispersed model
   return(output)                          }
