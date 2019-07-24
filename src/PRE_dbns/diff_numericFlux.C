diff --git a/numericFlux/numericFlux.C b/home/thomas/foam/foam-extend-4.0/src/dbns/numericFlux/numericFlux.C
index a86259e..5c53823 100644
--- a/numericFlux/numericFlux.C
+++ b/home/thomas/foam/foam-extend-4.0/src/dbns/numericFlux/numericFlux.C
@@ -25,7 +25,6 @@ License
 
 #include "numericFlux.H"
 #include "MDLimiter.H"
-#include "arbMesh.H"
 
 // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
 
@@ -86,7 +85,7 @@ Foam::numericFlux<Flux, Limiter>::numericFlux
 // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //
 
 template<class Flux, class Limiter>
-void Foam::numericFlux<Flux, Limiter>::computeFlux(arbMesh& amsh)
+void Foam::numericFlux<Flux, Limiter>::computeFlux()
 {
     // Get face-to-cell addressing: face area point from owner to neighbour
     const unallocLabelList& owner = this->mesh().owner();
@@ -140,36 +139,17 @@ void Foam::numericFlux<Flux, Limiter>::computeFlux(arbMesh& amsh)
     const volVectorField& ULimiter = vectorULimiter.phiLimiter();
     const volScalarField& TLimiter = scalarTLimiter.phiLimiter();
 
-    //Get mesh information
-    const fvMesh& mesh = this->mesh();
-    const surfaceVectorField& Cf = mesh.Cf();
-    vector xyztmp = vector::zero, xyztmp2 = vector::zero;
-    vector& xyzOwn = xyztmp;
-    vector& xyzNei = xyztmp2;
-    
     // Calculate fluxes at internal faces
     forAll (owner, faceI)
-    {		   
+    {
         const label own = owner[faceI];
         const label nei = neighbour[faceI];
-	
+
         const vector deltaRLeft = faceCentre[faceI] - cellCentre[own];
         const vector deltaRRight = faceCentre[faceI] - cellCentre[nei];
 
-	// MODIFICATIONS HERE
-	// get coordinates of cell center for own and nei
-	xyzOwn[0] = Cf[own].x();
-	xyzOwn[1] = Cf[own].y();
-	xyzOwn[2] = Cf[own].z();
-	xyzNei[0] = Cf[nei].x();
-	xyzNei[1] = Cf[nei].y();
-	xyzNei[2] = Cf[nei].z();
-
-	vector dotXtemp = amsh.vw(xyzOwn)/2 + amsh.vw(xyzNei)/2;
-	const vector& dotX = dotXtemp;
-	
         // calculate fluxes with reconstructed primitive variables at faces
-	Flux::evaluateFlux
+        Flux::evaluateFlux
         (
             rhoFlux_[faceI],
             rhoUFlux_[faceI],
@@ -185,15 +165,10 @@ void Foam::numericFlux<Flux, Limiter>::computeFlux(arbMesh& amsh)
             Cv[own],
             Cv[nei],
             Sf[faceI],
-            magSf[faceI],
-	    xyzOwn,
-	    xyzNei,
-	    amsh,
-	    faceI,
-	    dotX
+            magSf[faceI]
         );
     }
-    
+
     // Update boundary field and values
     forAll (rhoFlux_.boundaryField(), patchi)
     {
@@ -226,10 +201,6 @@ void Foam::numericFlux<Flux, Limiter>::computeFlux(arbMesh& amsh)
         const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
         const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];
 
-	
-	vector dotXtemp = amsh.vw(xyzOwn)/2 + amsh.vw(xyzNei)/2;
-	const vector& dotX = dotXtemp;
-	
         if (pp.coupled())
         {
             // Coupled patch
@@ -322,12 +293,7 @@ void Foam::numericFlux<Flux, Limiter>::computeFlux(arbMesh& amsh)
                     pCv[facei],
                     pCv[facei],
                     pSf[facei],
-                    pMagSf[facei],
-		    xyzOwn,
-		    xyzNei,
-		    amsh,
-		    facei,
-		    dotX
+                    pMagSf[facei]
                 );
             }
         }
@@ -352,13 +318,8 @@ void Foam::numericFlux<Flux, Limiter>::computeFlux(arbMesh& amsh)
                     pCv[facei],
                     pCv[facei],
                     pSf[facei],
-                    pMagSf[facei],
-		    xyzOwn,
-		    xyzNei,
-		    amsh,
-		    facei,
-		    dotX
-		 );
+                    pMagSf[facei]
+                );
             }
         }
     }
