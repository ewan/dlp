experiments <- list()
experiments$emmog <- mclustexpt("emmog",
  function(data, obs) {
    list(
      G=1:25,
      modelNames="VVV"
    )
  }
)
experiments$emmog$obsreq <- F
experiments$flgfa <- jqgexpt("flgfa",
  function(data, obs) {
    phi <- 200*diag(ncol(obs)+1)
    for (i in 2:(ncol(obs)+1)) {
      phi[i,i] <- 1
    }
    list(
      model.fn=flgfa,
      jqgmodel="org/jqgibbs/models/FLGFAModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(ncol(obs)+1,ncol(data))),
        S=cov(data),
        Psi=cov(data),
        kappa=as.double(ncol(data)+7),
        Phi=phi,
        lambda=as.double(ncol(obs)+1+20),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=2,
        alb=.5
      ),
      init=list(
        A=array(0, c(ncol(obs)+1,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=array(diag(ncol(obs)+1), c(ncol(obs)+1,ncol(obs)+1,1)),
        M=array(0, c(ncol(obs)+1,ncol(data))),
        B=matrix(as.integer(as.matrix(cbind(1,obs))),nrow(obs),ncol(obs)+1),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=5000,
      finalsave=F
    )
  }
)
experiments$flgfa$obsreq <- T
experiments$flgfaa <- jqgexpt("flgfaa",
  function(data, obs) {
    phi <- diag(ncol(obs)+1)
    for (i in 2:(ncol(obs)+1)) {
      phi[i,i] <- 0.5
    }
    list(
      model.fn=flgfa,
      jqgmodel="org/jqgibbs/models/FLGFAModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(ncol(obs)+1,ncol(data))),
        S=cov(data),
        Psi=cov(data),
        kappa=as.double(ncol(data)+7),
        Phi=phi,
        lambda=as.double(ncol(obs)+1+20),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=2,
        alb=.5
      ),
      init=list(
        A=array(0, c(ncol(obs)+1,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=array(diag(ncol(obs)+1), c(ncol(obs)+1,ncol(obs)+1,1)),
        M=array(0, c(ncol(obs)+1,ncol(data))),
        B=matrix(as.integer(as.matrix(cbind(1,obs))),nrow(obs),ncol(obs)+1),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=5000,
      finalsave=F
    )
  }
)
experiments$flgfaa$obsreq <- T
experiments$flgfa1a <- jqgexpt("flgfa1a",
  function(data, obs) {
    phi <- diag(2)
    phi[2,2] <- 0.5
    list(
      model.fn=flgfa,
      jqgmodel="org/jqgibbs/models/FLGFAModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(2,ncol(data))),
        S=cov(data),
        Psi=cov(data),
        kappa=as.double(ncol(data)+7),
        Phi=phi,
        lambda=(2+20),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=2,
        alb=.5
      ),
      init=list(
        A=array(0, c(2,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=array(diag(2), c(2,2,1)),
        M=array(0, c(2,ncol(data))),
        B=matrix(as.integer(as.matrix(cbind(1,obs[,1,drop=F]))),nrow(obs),2),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=5000,
      finalsave=T
    )
  }
)
experiments$flgfa1a$obsreq <- T
experiments$flgfa1 <- jqgexpt("flgfa1",
  function(data, obs) {
    phi <- 200*diag(2)
    phi[2,2] <- 1
    list(
      model.fn=flgfa,
      jqgmodel="org/jqgibbs/models/FLGFAModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(2,ncol(data))),
        S=cov(data),
        Psi=cov(data),
        kappa=as.double(ncol(data)+7),
        Phi=phi,
        lambda=as.double(2+20),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=2,
        alb=.5
      ),
      init=list(
        A=array(0, c(2,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=array(diag(2), c(2,2,1)),
        M=array(0, c(2,ncol(data))),
        B=matrix(as.integer(as.matrix(cbind(1,obs[,1,drop=F]))),nrow(obs),2),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=5000,
      finalsave=T
    )
  }
)
experiments$flgfa1$obsreq <- T
experiments$dpmod <- jqgexpt("dpmod",
  function(data, obs) {
    phi <- array(200, c(1,1))
    list(
      model.fn=flgfa,
      jqgmodel="org/jqgibbs/models/FLGFDModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(1,ncol(data))),
        S=cov(data),
        Psi=0.1*cov(data),
        kappa=as.double(ncol(data)),
        Phi=phi,
        lambda=as.double(22),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=2,
        alb=.5
      ),
      init=list(
        A=array(0, c(1,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=array(200, c(1,1)),
        M=array(0, c(1,ncol(data))),
        B=array(as.integer(1), c(nrow(data),1)),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=2000,
      finalsave=F
    )
  }
)
experiments$dpmod$obsreq <- F
experiments$dpmoda <- jqgexpt("dpmoda",
  function(data, obs) {
    phi <- array(1, c(1,1))
    list(
      model.fn=flgfa,
      jqgmodel="org/jqgibbs/models/FLGFDModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(1,ncol(data))),
        S=cov(data),
        Psi=0.1*cov(data),
        kappa=as.double(ncol(data)),
        Phi=phi,
        lambda=as.double(2+20),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=2,
        alb=.5
      ),
      init=list(
        A=array(0, c(1,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=array(1, c(1,1)),
        M=array(0, c(1,ncol(data))),
        B=array(as.integer(1), c(nrow(data),1)),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=5000,
      finalsave=F
    )
  }
)
experiments$dpmoda$obsreq <- F
experiments$dpmodb <- jqgexpt("dpmodb",
  function(data, obs) {
    phi <- array(1, c(1,1))
    list(
      model.fn=flgfa,
      jqgmodel="org/jqgibbs/models/FLGFDModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(1,ncol(data))),
        S=cov(data),
        Psi=cov(data),
        kappa=as.double(ncol(data)+1+10),
        Phi=phi,
        lambda=as.double(2+10),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=2,
        alb=.5
      ),
      init=list(
        A=array(0, c(1,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=array(1, c(1,1)),
        M=array(0, c(1,ncol(data))),
        B=array(1, c(nrow(data),1)),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=5000,
      finalsave=F
    )
  }
)
experiments$dpmodb$obsreq <- F
experiments$flgfd <- jqgexpt("flgfd",
  function(data, obs) {
    phi <- 200*diag(ncol(obs)+1)
    for (i in 2:(ncol(obs)+1)) {
      phi[i,i] <- 1
    }
    list(
      model.fn=flgfa,
      jqgmodel="org/jqgibbs/models/FLGFDModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(ncol(obs)+1,ncol(data))),
        S=cov(data),
        Psi=0.1*cov(data),
        kappa=as.double(ncol(data)),
        Phi=phi,
        lambda=as.double(ncol(obs)+1+20),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=2,
        alb=.5
      ),
      init=list(
        A=array(0, c(ncol(obs)+1,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=diag(ncol(obs)+1),
        M=array(0, c(ncol(obs)+1,ncol(data))),
        B=matrix(as.integer(as.matrix(cbind(1,obs))),nrow(obs),ncol(obs)+1),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=5000,
      finalsave=F
    )
  }
)
experiments$flgfd$obsreq <- T
experiments$flgfda <- jqgexpt("flgfda",
  function(data, obs) {
    phi <- diag(ncol(obs)+1)
    for (i in 2:(ncol(obs)+1)) {
      phi[i,i] <- 0.5
    }
    list(
      model.fn=flgfa,
      jqgmodel="org/jqgibbs/models/FLGFDModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(ncol(obs)+1,ncol(data))),
        S=cov(data),
        Psi=0.1*cov(data),
        kappa=as.double(ncol(data)),
        Phi=phi,
        lambda=as.double(ncol(obs)+1+20),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=2,
        alb=.5
      ),
      init=list(
        A=array(0, c(ncol(obs)+1,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=diag(ncol(obs)+1),
        M=array(0, c(ncol(obs)+1,ncol(data))),
        B=matrix(as.integer(as.matrix(cbind(1,obs))),nrow(obs),ncol(obs)+1),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=5000,
      finalsave=F
    )
  }
)
experiments$flgfda$obsreq <- T
experiments$flgfdb <- jqgexpt("flgfdb",
  function(data, obs) {
    phi <- diag(ncol(obs)+1)
    for (i in 2:(ncol(obs)+1)) {
      phi[i,i] <- 0.5
    }
    list(
      model.fn=flgfa,
      jqgmodel="org/jqgibbs/models/FLGFDModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(ncol(obs)+1,ncol(data))),
        S=cov(data),
        Psi=cov(data),
        kappa=as.double(ncol(data)+1+10),
        Phi=phi,
        lambda=as.double(ncol(obs)+1+10),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=1000,
        alb=1000
      ),
      init=list(
        A=array(0, c(ncol(obs)+1,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=diag(ncol(obs)+1),
        M=array(0, c(ncol(obs)+1,ncol(data))),
        B=as.matrix(cbind(1,obs)),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=5000,
      finalsave=F
    )
  }
)
experiments$flgfdb$obsreq <- T
experiments$flgfea <- jqgexpt("flgfea",
  function(data, obs) {
    phi <- diag(2)
    phi[2,2] <- 0.5
    list(
      model.fn=flgfe,
      jqgmodel="org/jqgibbs/models/FLGFEModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(2,ncol(data))),
        S=cov(data),
        Psi=0.1*cov(data),
        kappa=as.double(ncol(data)),
        Phi=phi,
        lambda=as.double(22),
        bshape1=360,
        bshape2=640,
        xa=1,
        xb=1,
        ala=2,
        alb=.5,
        be=4,
        ga=6
      ),
      init=list(
        A=array(0, c(2,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=diag(2),
        M=array(0, c(2,ncol(data))),
        B=as.integer(rep(0,nrow(data))),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=5000,
      finalsave=F
    )
  }
)
experiments$flgfea$obsreq <- F

experiments$flgfd1 <- jqgexpt("flgfd1",
  function(data, obs) {
    phi <- 200*diag(2)
    phi[2,2] <- 1
    list(
      model.fn=flgfd,
      jqgmodel="org/jqgibbs/models/FLGFDModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(2,ncol(data))),
        S=cov(data),
        Psi=0.1*cov(data),
        kappa=as.double(ncol(data)),
        Phi=phi,
        lambda=as.double(2+20),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=2,
        alb=.5
      ),
      init=list(
        A=array(0, c(2,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=phi,
        M=array(0, c(2,ncol(data))),
        B=matrix(as.integer(as.matrix(cbind(1,obs[,1,drop=F]))),nrow(obs),2),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=2000,
      finalsave=T
    )
  }
)
experiments$flgfd1$obsreq <- F

experiments$flgfd2a <- jqgexpt("flgfd2a",
  function(data, obs) {
    phi <- diag(3)
    phi[2,2] <- 0.5
    phi[3,3] <- 0.5
    list(
      model.fn=flgfa,
      jqgmodel="org/jqgibbs/models/FLGFDModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(3,ncol(data))),
        S=cov(data),
        Psi=0.1*cov(data),
        kappa=as.double(ncol(data)),
        Phi=phi,
        lambda=as.double(3+20),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=2,
        alb=.5
      ),
      init=list(
        A=array(0, c(3,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=phi,
        M=array(0, c(3,ncol(data))),
        B=matrix(as.integer(as.matrix(cbind(1,obs[,1:2,drop=F]))),nrow(obs),3),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=5000,
      finalsave=T
    )
  }
)
experiments$flgfd2a$obsreq <- F

experiments$flgfd1a <- jqgexpt("flgfd1a",
  function(data, obs) {
    phi <- diag(2)
    phi[2,2] <- 0.5
    list(
      model.fn=flgfa,
      jqgmodel="org/jqgibbs/models/FLGFDModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(2,ncol(data))),
        S=cov(data),
        Psi=0.1*cov(data),
        kappa=as.double(ncol(data)),
        Phi=phi,
        lambda=as.double(2+20),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=2,
        alb=.5
      ),
      init=list(
        A=array(0, c(2,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=phi,
        M=array(0, c(2,ncol(data))),
        B=matrix(as.integer(as.matrix(cbind(1,obs[,1,drop=F]))),nrow(obs),2),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=5000,
      finalsave=T
    )
  }
)
experiments$flgfd1a$obsreq <- F
experiments$flgfc1 <- jqgexpt("flgfc1",
  function(data, obs) {
    phi <- 200*diag(2)
    phi[2,2] <- 1
    list(
      model.fn=flgfa,
      jqgmodel="org/jqgibbs/models/FLGFCModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(2,ncol(data))),
        S=cov(data),
        Psi=cov(data),
        kappa=as.double(ncol(data)+7),
        Phi=phi,
        lambda=as.double(2+20),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=2,
        alb=.5
      ),
      init=list(
        A=array(0, c(2,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=array(diag(2), c(2,2,1)),
        M=array(0, c(2,ncol(data))),
        B=matrix(as.integer(as.matrix(cbind(1,sample(0:1,nrow(data),replace=T)))),nrow(data),2),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=5000,
      finalsave=T
    )
  }
)
experiments$flgfc1$obsreq <- F
experiments$flgfc1a <- jqgexpt("flgfc1a",
  function(data, obs) {
    phi <- diag(2)
    phi[2,2] <- 0.5
    list(
      model.fn=flgfa,
      jqgmodel="org/jqgibbs/models/FLGFCModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(2,ncol(data))),
        S=cov(data),
        Psi=cov(data),
        kappa=as.double(ncol(data)+7),
        Phi=phi,
        lambda=as.double(2+20),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=2,
        alb=.5
      ),
      init=list(
        A=array(0, c(2,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=array(diag(2), c(2,2,1)),
        M=array(0, c(2,ncol(data))),
        B=matrix(as.integer(as.matrix(cbind(1,sample(0:1,nrow(data),replace=T)))),nrow(data),2),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=5000,
      finalsave=T
    )
  }
)
experiments$flgfc1a$obsreq <- F
experiments$flgfc3 <- jqgexpt("flgfc3",
  function(data, obs) {
    phi <- 200*diag(4)
    phi[2,2] <- 1
    phi[3,3] <- 1
    phi[4,4] <- 1
    list(
      model.fn=flgfa,
      jqgmodel="org/jqgibbs/models/FLGFCModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(4,ncol(data))),
        S=cov(data),
        Psi=cov(data),
        kappa=as.double(ncol(data)+7),
        Phi=phi,
        lambda=(4+20),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=2,
        alb=.5
      ),
      init=list(
        A=array(0, c(4,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=array(diag(4), c(4,4,1)),
        M=array(0, c(4,ncol(data))),
        B=matrix(as.integer(as.matrix(cbind(1,sample(0:1,nrow(data),replace=T),
                                              sample(0:1,nrow(data),replace=T),
                                              sample(0:1,nrow(data),replace=T)))),nrow(data),4),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=5000,
      finalsave=T
    )
  }
)
experiments$flgfc3$obsreq <- F
experiments$dpmog <- jqgexpt("dpmog",
  function(data, obs) {
    list(
      model.fn=flgfa,
      jqgmodel="org/jqgibbs/models/FLGFAModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(1,ncol(data))),
        S=cov(data),
        Psi=cov(data),
        kappa=as.double(ncol(data)+7),
        Phi=array(200, c(1,1)),
        lambda=(22),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=2,
        alb=.5
      ),
      init=list(
        A=array(0, c(1,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=array(1, c(1,1,1)),
        M=array(0, c(1,ncol(data))),
        B=array(as.integer(1), c(nrow(data),1)),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=600,
      finalsave=F
    )
  }
)
experiments$dpmog$obsreq <- F
experiments$dpmoga <- jqgexpt("dpmoga",
  function(data, obs) {
    list(
      model.fn=flgfa,
      jqgmodel="org/jqgibbs/models/FLGFAModel",
      samplerClass="org/jqgibbs/GenericSampler",
      hypers=list(
        W=array(0, c(1,ncol(data))),
        S=cov(data),
        Psi=cov(data),
        kappa=as.double(ncol(data)+7),
        Phi=array(1, c(1,1)),
        lambda=(20),
        bshape1=1,
        bshape2=1,
        xa=1,
        xb=1,
        ala=2,
        alb=.5
      ),
      init=list(
        A=array(0, c(1,ncol(data),1)),
        Sg=array(diag(ncol(data)), c(ncol(data),ncol(data),1)),
        Omega=array(1, c(1,1,1)),
        M=array(0, c(1,ncol(data))),
        B=array(as.integer(1), c(nrow(data),1)),
        Z=as.integer(rep(0,nrow(data))),
        x=0.5,
        al=1
      ),
      length=500,
      burnin=1200,
      lag=10,
      savetime=600,
      finalsave=F
    )
  }
)
experiments$dpmoga$obsreq <- F

