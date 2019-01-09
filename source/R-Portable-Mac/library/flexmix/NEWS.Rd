\name{NEWS}
\title{News for Package 'flexmix'}
\section{Changes in flexmix version 2.3-14}{
  \itemize{
    \item A bug fixed for \code{FLXMRcondlogit()} to ensure that the
    formula for the strata is also stored and can be used for
    predicting, etc. new data. Thanks to Peter Fraser-Mackenziefor
    for pointing the problem out.
    \item A bug fixed for \code{FLXMRglmfix()} which occurred if several
    components were simultaneously omitted. Thanks to Moritz Berger for
    pointing the problem out.
    \item For mixtures of mixed effects models and with censored data
    only weighted ML estimation is implemented. This is now checked and
    an error is given if a different estimation method is requested.
    \item A generic function \code{sigma()} is available for R >= 3.3.0
    and thus different versions of \code{sigma()} need to be called
    depending on the R version. Thanks also to Stephen Martin for
    pointing the issue out.
    \item Components are now generated using functions instead of
    expressions.
    \item Functions from the base packages \pkg{stats}, \pkg{graphics}
    and \pkg{grDevices} are now correctly imported before being used.
    \item Function \code{FLXMCdist1} implements model drivers for
    univariate mixtures of different distributions.
  }
}
\section{Changes in flexmix version 2.3-13}{
  \itemize{
    \item A model driver for mixtures of regression models with
      (adaptive) lasso and elastic net penalizations for the coefficients
      building on \pkg{glmnet} was provided by Frederic Mortier and
      Nicolas Picard.
    \item A bug in a coerce method to class \code{"FLXnested"} fixed.
  }
}
\section{Changes in flexmix version 2.3-12}{
  \itemize{
    \item Construction of model matrices changed to re-use levels of
      factors while fitting for prediction. Thanks to Robert Poos for
      pointing the problem out.
    \item Package mgcv of version at least 1.8-0 is required in order
      to avoid loading of the package for formula evaluation.
   }
}
\section{Changes in flexmix version 2.3-11}{
  \itemize{
    \item Examples changed to be less time consuming.
    \item Bug fixed in ungroupPriors() and getPriors() to work with a
      grouping where only one unique value is contained. Thanks to
      Christine Koehler for pointing the problem out.
    \item The \code{logLik()} method for \code{"flexmix"} objects now
      also has a \code{newdata} argument.
    \item In the M-step only the parameters of the previously fitted
      components are passed over to avoid nesting of
      environments. Thanks to Dominik Ernst for pointing the problem
      out.
    \item Function \code{boot()} was fixed to work correctly with
      grouped data if the groups are kept and if fixed effects are
      fitted.
    \item Authors@R used in DESCRIPTION. Deepayan Sarkar listed as
      contributor due to internal code copied from package lattice used
      for the plots in flexmix.
    \item Model driver \code{FLXMRlmer()} adapted for \pkg{lme4} (>=
      1.0).
  }
}
\section{Changes in flexmix version 2.3-10}{
  \itemize{
    \item Parallel processing disabled in \code{stepFlexmix()}.
  }
}
\section{Changes in flexmix version 2.3-9}{
  \itemize{
    \item Package dependencies, imports and suggests modified suitably
      for R >= 2.15.0.
    \item NEWS file adapted to a NEWS.Rd file.
    \item Parallel processing is enabled in \code{stepFlexmix()}. See
      \pkg{flexclust} for details.
    \item New model drivers \code{FLXMRmultinom()} and \code{FLXMRcondlogit()} are
      included which allow to fit finite mixtures of multinomial logit
      and conditional logit models. Identifiability problems might arise
      for this model class. Details on sufficient identifiability
      conditions are given in Gruen and Leisch (2008).
    \item A bug in \code{FLXMRlmm()} was fixed which did not allow to specify
      restrictions on the variances of the random effects and / or the
      residuals. Thanks to Gregory Matthews for pointing the problem out.
  }
}
\section{Changes in flexmix version 2.3-8}{
  \itemize{
    \item The fit function in the M-step by default now is called with an
      argument containing the fitted component. This allows to use the
      parameter estimates from the previous step for initialization. Fit
      functions which do not require this now need a \code{...} argument.
      Thanks to Hannah Frick and Achim Zeileis for requesting this
      feature.
    \item Function \code{initFlexmix()} was added which is an alternative to
      \code{stepFlexmix()} where first several short runs of EM, SEM or CEM are
      performed followed by a long run of EM. 
  }
}
\section{Changes in flexmix version 2.3-7}{
  \itemize{
    \item A bug fixed in \code{predict()} and \code{fitted()} if a concomitant variable
      model is specified and \code{aggregate = TRUE}. Thanks to Julia Schiffner
      for pointing this out.
    \item A bug fixed in \code{FLXMRmgcv()} if observations were removed in the
      M-step because their a-posterior probabilities were smaller than
      eps. Thanks to Ghislain Geniaux for pointing the problem out.
  }
}
\section{Changes in flexmix version 2.3-6}{
  \itemize{
    \item Vignettes moved from /inst/doc to /vignettes.
    \item \code{stepFlexmix()} can now be called with a concomitant variable model
      \code{FLXPmultinom()} for \code{k} starting with 1 without getting an error. The
      concomitant variable model is internally replaced by
      \code{FLXPconstant()}.
    \item The \code{boot()} method for \code{"flexmix"} objects is extended to mixture
      models with concomitant variables and mixtures of linear mixed
      models.
    \item A bug fixed in the \code{summary()} method for \code{"flexmix"} objects. The
      column post > 0 did not give the correct results if weights were
      used for fitting the mixture.
    \item A bug fixed in the \code{unique()} method for \code{"stepFlexmix"} objects. This
      only occurred if components were dropped as well as if the EM
      algorithm did not converge for all repetitions. Thanks to
      Sebastian Meyer for pointing this out.
  }
}
\section{Changes in flexmix version 2.3-5}{
  \itemize{
    \item A bug fixed in \code{posterior()}. Fixed priors were always used, also if
      a concomitant variable model was present.
    \item A method added for \code{prior()} such that if newdata is supplied and
      the object is of class \code{"flexmix"} the prior class probabilities for
      each observation are returned.
  }
}
\section{Changes in flexmix version 2.3-4}{
  \itemize{
    \item A generic method for \code{nobs()} is introduced in \pkg{stats4} for \R
      2.13.0. \pkg{flexmix} now does not define this generic function and
      \code{logLik()}, \code{AIC()} and \code{BIC()} methods were modified to better exploit
      already available methods. Thanks to Prof. Brian D. Ripley for
      suggesting the modification.
  }
}
\section{Changes in flexmix version 2.3-3}{
  \itemize{
    \item A bug for \code{boot()} fixed for \code{"flexmix"} objects with an unbalanced 
      grouping variable. Thanks to Laszlo Sandor for pointing this out.
  }
}
\section{Changes in flexmix version 2.3-2}{
  \itemize{
    \item A bug for \code{rflexmix()} fixed for \code{"flexmix"} objects with a
      concomitant variable model. Thanks to Greg Petroski for pointing
      this out.
  }
}
\section{Changes in flexmix version 2.3-1}{
  \itemize{
    \item Functionality for bootstrapping finite mixture models added.
  }
}
\section{Changes in flexmix version 2.2-11}{
  \itemize{
    \item More generics and methods exported to use the \code{refit()} method when
      extending \pkg{flexmix} in other packages.
  }
}
\section{Changes in flexmix version 2.2-10}{
  \itemize{
    \item For long formulas \code{FLXMRglmfix()} did not work properly due to the
      splitting of the formula into several parts by \code{deparse()}. This is
      fixed by pasting them together again. Thanks to Dustin Tingley for
      the bug report.
  }
}
\section{Changes in flexmix version 2.2-9}{
  \itemize{
    \item A new model driver \code{FLXMRmgcv()} is added which allows to fit
      regularized linear models in the components.
    \item More generics and methods exported to allow extending \pkg{flexmix} in
      other packages.
  }
}
\section{Changes in flexmix version 2.2-8}{
  \itemize{
    \item The a-posteriori probabilities are now also determined as changed
      for \code{FLXfit()} for version 2.2-6 for \code{refit()}.
    \item Bug fixed for FLXfillconcomitant and refit when weights and
      grouping are present. A check was added that weights are identical
      within groups.
    \item Function \code{group()} is now exported.
  }
}
\section{Changes in flexmix version 2.2-7}{
  \itemize{
    \item Bug in the \code{FLXgetModelmatrix()} method for the \code{"FLXMRlmm"} class
      fixed when determining identical random effects covariates for the
      grouping.
    \item A new model driver for finite mixtures of linear mixed effects
      models with left-censored observations is added.
  }
}
\section{Changes in flexmix version 2.2-6}{
  \itemize{
    \item Determination of the a-posteriori probabilities made numerically
      more stable for small likelihoods. Thanks to Nicolas Picard for
      the code patch.
    \item \code{summary()} for \code{"FLXMRmstep"} objects now returns similar output for
      \code{which = "concomitant"} as for \pkg{flexmix} version 2.0-1.
    \item New demo driver \code{FLXMCnorm1()} for univariate Gaussian clustering.
    \item Non-postive values for the maximum number of iterations for the
      \code{"FLXcontrol"} object are not valid. A validity check for this is
      now included.
  }
}
\section{Changes in flexmix version 2.2-5}{
  \itemize{
    \item Model class \code{"FLXMRfix"} introduced which is a subclass
      of \code{"FLXMR"} and a superclass for \code{"FLXMRglmfix"} which also
      extends \code{"FLXMRglm"}.
    \item Model driver \code{FLXMCfactanal()} added which allows to
      fit finite mixtures of Gaussian distributions where the
      variance-covariance matrix is estimated using factor analyzers.
    \item Comparison of formulas now done using \code{identical()}.
  }
}
\section{Changes in flexmix version 2.2-4}{
  \itemize{
    \item Model drivers \code{FLXMRlmer()} and \code{FLXMRlmm()} added for fitting finite
      mixtures of linear mixed effects models.
    \item \code{EIC()} added as additional information criterion for assessing model fit.
    \item Bug fixed in \code{plot()} method for \code{"flexmix"} objects introduced in version 2.2-3.
  }
}
\section{Changes in flexmix version 2.2-3}{
  \itemize{
    \item New model driver \code{FLXMCmvcombi()} which is a combination of 
      Gaussian and binary. 
    \item \code{parameters()} now also has a which argument in order to allow to
      access the parameters of the concomitant variable model.
    \item Bug fixed in \code{refit()}.
    \item \code{nobs()} now returns the number of rows in the data.frame and not the
      number of individuals (similar as for example by lme).
  }
}
\section{Changes in flexmix version 2.2-0}{
  \itemize{
    \item vignette describing Version 2 added
    \item isTRUE(\code{all.equal()}) replaced with \code{identical()}.
    \item Bug fixed for prior in \code{flexmix()}.
    \item New function \code{relabel()} to sort components (generic is in modeltools).
    \item New example data generator \code{ExLinear()}.
    \item Fixed a bug in handling groups (gave an error for empty design matrices).
    \item Added new model \code{FLXMRrobglm()} for robust estimation of GLMs.
  }
}
\section{Changes in flexmix version 2.1-0}{
  \itemize{
    \item Renamed \code{cluster()} to \code{clusters()} to avoid conflict with \code{cluster()} 
      from package survival
    \item Bug fixed in internal functions using S4 generics and methods.
  }
}
\section{Changes in flexmix version 2.0-2}{
  \itemize{
    \item \code{refit()} now has a method argument. For method \code{"optim"} the
      variance-covariance matrix is determined using \code{optim()} to maximize the
      likelihood initialized in the solution found by the EM algorithm. 
      Method \code{"mstep"} refits the component specific and concomitant models 
      treating the posterior probabilities as given, i.e. performs an M-step 
      of the EM algorithm.
  }
}
\section{Changes in flexmix version 2.0-1}{
  \itemize{
    \item \code{Lapply()} added which allows to apply a function to each component
      of a finite mixture
    \item \code{KLdiv()} modified to allow for determination with a discrete and a
      continuous version of the KL divergence
  }
}
\section{Changes in flexmix version 2.0-0}{
  \itemize{
    \item Model driver for zero-inflated component specific models.
    \item Latent class analysis for binary multivariate data is now
      possible to estimate for truncated data where the number of 
      observations with pattern only zeros is missing.
    \item new argument newdata for \code{cluster()}
    \item new \code{unique()} method for \code{"stepFlexmix"} objects
  }
}
\section{Changes in flexmix version 1.9-0}{
  \itemize{
    \item New class definitions for component specific models and
      concomitant variable models.
    \item \code{fitted()} and \code{predict()} now have an aggregate argument in order to
      be able to determine the aggregated values over all components.
    \item The package has now a vignette presenting several applications of
      finite mixtures of regression models with varying and fixed effects 
      on artificial and real data which can be a accessed using the
      command vignette("regression-examples").
    \item The vignette "flexmix-intro" was adapted to reflect the changes
      made in the package.
    \item \code{stepFlexmix()} now returns an object of class \code{"stepFlexmix"} which has
      a \code{print()} and \code{plot()} method. In addition \code{getModel()} can be used to
      select an appropriate model.
    \item \code{flexmix()} now has a weights argument for multiple identical
      observations.
    \item New model drivers for latent class analysis with Bernoulli and
      Poisson distributed multivariate observations.
    \item Variants of the EM algorithm have been modified to correspond to
      CEM and SEM. These names can now also be used for specifying the
      classify slot of the \code{"FLXcontrol"} object.
  }
}
\section{Changes in flexmix version 1.8-1}{
  \itemize{
    \item The package can now fit concomitant variable models.
    \item New M-step driver for regression models with varying and fixed
      effects.
    \item ICL information criterion added.
  }
}
\section{Changes in flexmix version 1.1-2}{
  \itemize{
    \item Fixed a bug that made the log-likelihood infinity for observations
      where all posteriors are numerically zero.
    \item Fixed a bug for formulae with dots.
    \item \code{posterior()} now has a newdata argument.
    \item New demo driver for model-based clustering of binary data.
    \item Adapted to changes in \code{summary.glm()} of \R version 2.3.0.
  }
}
\section{Changes in flexmix version 1.1-1}{
  \itemize{
    \item The \code{cluster} argument of \code{flexmix()} may now also be a matrix of 
      posterior probabilities.
    \item Fixed a bug to make size table work in case of empty clusters.
    \item Fixed a bug in likelihood computation for grouped observations.
    \item The artificial NPreg data now also have a binomial response, 
      added example to help(\code{"flexmix"}).
  }
}
\section{Changes in flexmix version 1.1-0}{
  \itemize{
    \item The \code{FLXglm()} driver now has an offset argument.
    \item New data set seizure as example for a Poisson GLM
      with an offset.
    \item \code{fitted()} can be used to extract fitted values from \code{"flexmix"}
      and \code{"FLXrefit"} objects.
    \item New accessor methods \code{cluster()} and \code{posterior()}.
    \item The package now uses lazy loading and has a namespace.
  }
}
\section{Changes in flexmix version 1.0-0}{
  \itemize{
    \item The package has now an introductionary vignette which can be
      accessed using the command vignette("flexmix-intro"). The vignette
      has been published in the Journal of Statistical Software, Volume
      11, Issue 8 (\url{http://www.jstatsoft.org/v11/i08}), and the paper should be used as
      citation for \pkg{flexmix}, run \code{citation("flexmix")} in \R 2.0.0 or newer
      for details.
    \item Several typos in help pages have been fixed.
  }
}
\section{Changes in flexmix version 0.9-1}{
  \itemize{
    \item Adjusted for \R 2.0.0.
    \item Fixed a bug in the \code{summary()} and \code{plot()} methods
      of \code{"flexmix"} objects in case of empty clusters.
    \item \code{stepFlexmix()} takes two new arguments: \code{compare} allows fo find
      minimum AIC/BIC solutions in addition to maximum likelihood,
      \code{verbose} gives some information about progress.
    \item Use a default of \code{verbose = 0} in \code{FLXcontrol()} (better in combination 
      with \code{stepFlexmix()}). 
  }
}
\section{Changes in flexmix version 0.9-0}{
  \itemize{
    \item New \code{summary()} and \code{plot()} methods for \code{"flexmix"} objects.
    \item \code{"FLXglm"} objects can now be automatically \code{refit()}ted to get table of
      significance tests for coefficients.
    \item New function \code{stepFlexmix()} for more automated model search.
    \item The artificial example data now have functions to
      create them and a pre-stored data sets, new function \code{plotEll()} to
      plot 2d Gaussians with confidence ellipses.
    \item New function \code{KLdiv()} to compute Kullback-Leibler divergence of
      components.
    \item The calculation of the degrees of freedom for \code{FLXmclust()} was wrong
  }
}
\section{Changes in flexmix version 0.7-1}{
  \itemize{
    \item Fixed some codoc problems (missing aliases).
    \item First version released on CRAN: 0.7-0
  }
}
