Model.Counts <-
function(parameter,model.type,model.name,link,covariates.matrix.mean,
                   covariates.matrix.scalef,offset.mean,offset.scalef,
                   fixed.b,vnmax) {
   if (model.type=="mean only") { 
                   output <- Model.Faddy(parameter,model.name,link,covariates.matrix.mean,
                                       offset.mean,fixed.b,vnmax)
                                } # end of mean only model
   if (model.type=="mean and scale-factor") { 
          if ((model.name=="general") | (model.name=="general fixed b")) {
                   output <- Model.FaddyJMV.general(parameter,link,covariates.matrix.mean,
                                         covariates.matrix.scalef,offset.mean,
                                         offset.scalef,fixed.b,vnmax) }
          if (model.name=="limiting") { 
                   output <- Model.FaddyJMV.limiting(parameter,link,covariates.matrix.mean,
                                         covariates.matrix.scalef,offset.mean,
                                         offset.scalef,vnmax) }
                                } # end of mean and scale-factor model
   return(output)                          }
