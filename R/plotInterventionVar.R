plotInterventionVar <- function(estIntVars, trueIntVars = NULL){
  if(!is.null(trueIntVars)){
    gg.df <- melt(as.data.frame(estIntVars))
    gg.df <- cbind(gg.df, rep(1:nrow(estIntVars), times = ncol(estIntVars)))
    
    
    true.vars <- melt(trueIntVars)[,3]
    gg.df <- cbind(gg.df, true.vars)
    
    gg.df <- as.data.frame(gg.df)
    
    colnames(gg.df) <- c("variable", "est.variance", "env", "true.variance")
    
    gg.df <- melt(gg.df, c("variable","env"), variable.name = "Variance")
    
    p <- ggplot(gg.df, aes(x = env, y = value, aes(group = variable, color = Variance)))
    p <- p + theme_classic(base_size = 8)
    p <- p + geom_point(aes(shape = Variance, color = Variance), size = 0.5)
    p <- p + geom_line(aes(type = Variance, color = Variance), size = 0.5)
    p <- p + facet_wrap(~variable, scales = "free") 
    p <- p + scale_x_discrete(labels = 1:nrow(estIntVars))
    p <- p + theme(axis.text.x=element_text(angle=90, hjust=1))
  }else{
    gg.df <- melt(as.data.frame(estIntVars))
    gg.df <- cbind(gg.df, rep(1:nrow(estIntVars), times = ncol(estIntVars)))
    gg.df <- as.data.frame(gg.df)
    colnames(gg.df) <- c("variable", "est.variance", "env")
    
    p <- ggplot(gg.df, aes(x = env, y = est.variance, aes(group = variable))) + theme_classic(base_size = 8)
    p <- p + geom_line(aes(group = variable, color = variable), size = 0.5)
    p <- p + facet_wrap(~variable, scales = "free") 
    p <- p + scale_x_discrete(labels = 1:nrow(estIntVars))
    p <- p + theme(axis.text.x=element_text(angle=90, hjust=1))
  }
  p
}