LL.gradient <-
function(parameter,model.type,model.name,link=link,list.data,
                     covariates.matrix.mean,covariates.matrix.scalef,
                     offset.mean,offset.scalef,ltvalue,utvalue,
                     fixed.b,weights,grad.method) {
   if (grad.method=="simple") { method.args=list(eps=1.e-4)
                       } else { method.args=list(r=6,eps=1.e-4) }
   gradient <- grad(LL.Regression.Counts,x=parameter,
                model.type=model.type,model.name=model.name,link=link,list.data=list.data,
                covariates.matrix.mean=covariates.matrix.mean,
                covariates.matrix.scalef=covariates.matrix.scalef,
                offset.mean=offset.mean,offset.scalef=offset.scalef,
                ltvalue=ltvalue,utvalue=utvalue,
                fixed.b=fixed.b,weights=weights,
                method=grad.method,method.args=method.args)
   return(gradient) }
