Author: D Yoan L Mekontchou Yombat 

## Viscous Fluid Registration Naive Implementation (Navier Stockes Based PDE)

This repository contains a naive matlab implementation of the Viscous Fluid
model for Image Registration.

## Goal And Aim of Image Registration (Non-Rigid/Non-Linear)

Nonlinear registration is the warping of one object onto another, using measures
of image similarity or corresponding features extracted from the images to guide
the deformation. It is widely used in brain imaging for computational anatomy
studies. Tensor-based morphometry (TBM), for example, identifies systematic
differences in brain structure by statistically analyzing deformations that align
images from many subjects to a common template.

Registration methods commonly combine two terms: a similarity measure
(distance or measure of agreement between two images) that drives the transformation
and a regularizer that ensures its smoothness. This extra term is
added to the registration function being optimized, to enforce desirable transformation
properties such as smoothness, invertibility and inverse-consistency.
For instance, the similarity criterion is regarded as a body force introduced into
mechanical equations that govern linear elastic motion (Hooke’s Law) in or
viscous fluid equations (Navier-Stokes equation) in. Other algorithms rely on
Gaussian filtering or enforce particular properties of the deformation such
as diffeomorphic trajectories.

When registering structural magnetic resonance brain images, the information
available (voxel intensity, pre-defined landmarks) is rather limited and correspondence
mappings are not unique. Consequently, a realistic model is needed
to achieve deformations that are closer to an independently defined ground truth.
This can be done for instance, as we chose to do here, by incorporating statistical
information on the data set into the deformation.

## Viscous Image Registratin

The goal of fluid registration is to determine a mapping from one image (the target) to another image (the source) in the form of a displacement field, u, defined throughout the target volume. The way this problem is set up has been described in some detail in previous publications (Christensen et al 1996, Freeborough and Fox 1998, Crum et al 2001) to which we refer the interested reader for more detail. Briefly though, the displacement field is modelled as a time-dependent viscous fluid flow on the target image, driven by image-derived forces that act to improve a measure of image similarity. The response of the fluid to the applied forces at an instant is obtained by solving the Navier–Stokes equation for a compressible viscous fluid. Such a fluid does not exist in nature but is a convenient transformation model for registration.

## This Repository

Contained in this repository is a sample solution of the PDE driving this model
and serves as a reference for any one currently performing research in the
field.


# TODO 
1. Implement Linear Interpolation For Reference frame mapping
2. Work On Parameter Optimization



