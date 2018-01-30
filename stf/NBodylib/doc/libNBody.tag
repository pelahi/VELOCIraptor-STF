<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile>
  <compound kind="file">
    <name>Analysis.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Analysis/</path>
    <filename>Analysis_8h</filename>
    <includes id="Density_8h" name="Density.h" local="no" imported="no">Density.h</includes>
    <includes id="Energy_8h" name="Energy.h" local="no" imported="no">Energy.h</includes>
    <includes id="Morphology_8h" name="Morphology.h" local="no" imported="no">Morphology.h</includes>
    <includes id="Power_8h" name="Power.h" local="no" imported="no">Power.h</includes>
  </compound>
  <compound kind="file">
    <name>Density.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Analysis/</path>
    <filename>Density_8cxx</filename>
    <includes id="NBody_8h" name="NBody.h" local="no" imported="no">NBody.h</includes>
    <includes id="Density_8h" name="Density.h" local="no" imported="no">Density.h</includes>
    <includes id="KDTree_8h" name="KDTree.h" local="no" imported="no">KDTree.h</includes>
    <namespace>NBody</namespace>
    <member kind="function">
      <type>void</type>
      <name>DensitySmooth</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a50be98286eb2be7bcf0c5ecee4014290</anchor>
      <arglist>(System &amp;S, int BucketSize, int NumNearest)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Density.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Analysis/</path>
    <filename>Density_8h</filename>
    <includes id="NBody_8h" name="NBody.h" local="no" imported="no">NBody.h</includes>
    <includes id="NBodyMath_8h" name="NBodyMath.h" local="no" imported="no">NBodyMath.h</includes>
    <namespace>NBody</namespace>
    <member kind="function">
      <type>void</type>
      <name>DensityBinEqualMass</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a33ec5f91802b46dbafee873789ae1c14</anchor>
      <arglist>(System &amp;S, int NumBins, Double_t *RBin, Double_t *DBin)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DensityBinLog</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ac9c4e8f661670737da3d02d478d209b2</anchor>
      <arglist>(const System &amp;S, int NumBins, Double_t *RBin, Double_t *DBin)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DensityBinLogEllip</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>af98d9605f8e1837b03a38780d9eb939f</anchor>
      <arglist>(const System &amp;S, int NumBins, Double_t *RBin, Double_t *DBin, Double_t q, Double_t s)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DensityBinNorm</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a9f216fbe3fb75d346ade3b7c5c0051cc</anchor>
      <arglist>(const System &amp;S, int NumBins, Double_t *RBin, Double_t *DBin)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DensitySmooth</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a50be98286eb2be7bcf0c5ecee4014290</anchor>
      <arglist>(System &amp;S, int BucketSize, int NumNearest)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DensitySmoothBall</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a288b85c5e87864a1b58565517e8631d3</anchor>
      <arglist>(System &amp;S, Double_t R)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Energy.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Analysis/</path>
    <filename>Energy_8h</filename>
    <includes id="NBody_8h" name="NBody.h" local="no" imported="no">NBody.h</includes>
    <namespace>NBody</namespace>
    <member kind="enumeration">
      <type></type>
      <name>Orbit</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a392718bf49bf3a990af266aba7b82079</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>Radial</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a392718bf49bf3a990af266aba7b82079a8865a16d55e6b8111c5a511abf31d2eb</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>Isotropic</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a392718bf49bf3a990af266aba7b82079a208f65742765cbe9f1a8791aa40124d2</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>CalcEnergy</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ab06ccd1e9ac542d1220cc9f31bf8cd09</anchor>
      <arglist>(const System &amp;S, Double_t *E, Double_t &amp;Etot)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>CalcPotDirect</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a993e7fc42b40e66b1724424b6a905da0</anchor>
      <arglist>(System &amp;S, Double_t eps=0.0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>CalcPotShell</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a80628dea731cb2661ec3e8c576802dd9</anchor>
      <arglist>(System &amp;S, Double_t eps=0.0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DenStates</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a6505f11a44a67cad8b78721c70c0afe1</anchor>
      <arglist>(const System &amp;S, Double_t *g, const Double_t *Ebin, int NumBins, Orbit orbit)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DiffEnergyDistEqual</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>af5520cd37fea91613efa03b151f20445</anchor>
      <arglist>(System &amp;S, const Double_t *E, Double_t *Ebin, Double_t *Mbin, int NumBins)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DiffEnergyDistNorm</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>afc959d8071bb31d168af94437dd0bff1</anchor>
      <arglist>(const System &amp;S, const Double_t *E, Double_t *Ebin, Double_t *Mbin, int NumBins)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>NumberDensity</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a99346c99c7e4ce25703996edc297c277</anchor>
      <arglist>(const System &amp;S, const Double_t *E, Double_t *Ebin, Double_t *Jbin, Double_t *Nbin, int NumBins)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Morphology.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Analysis/</path>
    <filename>Morphology_8cxx</filename>
    <includes id="Morphology_8h" name="Morphology.h" local="no" imported="no">Morphology.h</includes>
    <namespace>NBody</namespace>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphology</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a68397ab51b06718e9cf3ad77239b8a0e</anchor>
      <arglist>(const Int_t nbodies, Particle *p, Double_t &amp;q, Double_t &amp;s, Double_t Error, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphology</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a6c01dbea286d85b56ebcfdbe3ec5b845</anchor>
      <arglist>(const Int_t nbodies, Particle *p, Double_t &amp;q, Double_t &amp;s, Double_t Error, Matrix &amp;eigenvec, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphology</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aaf436aa1d0e2e8876cd3f14c4a169be3</anchor>
      <arglist>(System &amp;S, Double_t &amp;q, Double_t &amp;s, Double_t Error, int verbose=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphology</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aea80e105c9f44f5ae05535fefe9eef50</anchor>
      <arglist>(System &amp;S, Double_t &amp;q, Double_t &amp;s, Double_t Error, Matrix &amp;eigenvec, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphologyWithMass</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aee812d1a04cbad75feed79462dada9c6</anchor>
      <arglist>(const Int_t nbodies, Particle *p, Double_t &amp;q, Double_t &amp;s, Double_t Error, int mdenom=0, int verbose=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphologyWithMass</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a0bd395e9e22d1a7fb1958f91cd527adf</anchor>
      <arglist>(const Int_t nbodies, Particle *p, Double_t &amp;q, Double_t &amp;s, Double_t Error, Matrix &amp;eigenvec, int mdenom, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphologyWithMass</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aa7e450ff359b11a2a4be90915541603f</anchor>
      <arglist>(System &amp;S, Double_t &amp;q, Double_t &amp;s, Double_t Error, int mdenom, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphologyWithMass</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a64b03d8b923b719741f49f89cb7b8513</anchor>
      <arglist>(System &amp;S, Double_t &amp;q, Double_t &amp;s, Double_t Error, Matrix &amp;eigenvec, int mdenom, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetInertiaTensor</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a7ad395a51dc78c9a00653bb054145ea5</anchor>
      <arglist>(const Int_t nbodies, Particle *p, Double_t &amp;a, Double_t &amp;b, Double_t &amp;c, Matrix &amp;eigvec)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetInertiaTensor</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a87f960a3886d621231e641ed13b98a66</anchor>
      <arglist>(const Int_t n, Particle *p, Double_t &amp;a, Double_t &amp;b, Double_t &amp;c, Matrix &amp;eigenvec, Matrix &amp;I)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetInertiaTensor</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ab68b7ccda248a59fcaadda35c962569b</anchor>
      <arglist>(System &amp;S, Double_t &amp;a, Double_t &amp;b, Double_t &amp;c, Matrix &amp;eigvec)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetInertiaTensor</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a8b1d672d1bf4b868f1087809f54a882e</anchor>
      <arglist>(System &amp;S, Double_t &amp;a, Double_t &amp;b, Double_t &amp;c, Matrix &amp;eigenvec, Matrix &amp;I)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>Rotate</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a2ce569f54674485a595c76cff99110c1</anchor>
      <arglist>(const Matrix R, Coordinate x)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Morphology.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Analysis/</path>
    <filename>Morphology_8h</filename>
    <includes id="NBodyMath_8h" name="NBodyMath.h" local="no" imported="no">NBodyMath.h</includes>
    <includes id="NBody_8h" name="NBody.h" local="no" imported="no">NBody.h</includes>
    <namespace>NBody</namespace>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphology</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a68397ab51b06718e9cf3ad77239b8a0e</anchor>
      <arglist>(const Int_t nbodies, Particle *p, Double_t &amp;q, Double_t &amp;s, Double_t Error, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphology</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a6c01dbea286d85b56ebcfdbe3ec5b845</anchor>
      <arglist>(const Int_t nbodies, Particle *p, Double_t &amp;q, Double_t &amp;s, Double_t Error, Matrix &amp;eigenvec, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphology</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aaf436aa1d0e2e8876cd3f14c4a169be3</anchor>
      <arglist>(System &amp;S, Double_t &amp;q, Double_t &amp;s, Double_t Error, int verbose=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphology</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aea80e105c9f44f5ae05535fefe9eef50</anchor>
      <arglist>(System &amp;S, Double_t &amp;q, Double_t &amp;s, Double_t Error, Matrix &amp;eigenvec, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphologyWithMass</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aee812d1a04cbad75feed79462dada9c6</anchor>
      <arglist>(const Int_t nbodies, Particle *p, Double_t &amp;q, Double_t &amp;s, Double_t Error, int mdenom=0, int verbose=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphologyWithMass</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a0bd395e9e22d1a7fb1958f91cd527adf</anchor>
      <arglist>(const Int_t nbodies, Particle *p, Double_t &amp;q, Double_t &amp;s, Double_t Error, Matrix &amp;eigenvec, int mdenom, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphologyWithMass</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aa7e450ff359b11a2a4be90915541603f</anchor>
      <arglist>(System &amp;S, Double_t &amp;q, Double_t &amp;s, Double_t Error, int mdenom, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphologyWithMass</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a64b03d8b923b719741f49f89cb7b8513</anchor>
      <arglist>(System &amp;S, Double_t &amp;q, Double_t &amp;s, Double_t Error, Matrix &amp;eigenvec, int mdenom, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetInertiaTensor</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a7ad395a51dc78c9a00653bb054145ea5</anchor>
      <arglist>(const Int_t nbodies, Particle *p, Double_t &amp;a, Double_t &amp;b, Double_t &amp;c, Matrix &amp;eigvec)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetInertiaTensor</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a87f960a3886d621231e641ed13b98a66</anchor>
      <arglist>(const Int_t n, Particle *p, Double_t &amp;a, Double_t &amp;b, Double_t &amp;c, Matrix &amp;eigenvec, Matrix &amp;I)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetInertiaTensor</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ab68b7ccda248a59fcaadda35c962569b</anchor>
      <arglist>(System &amp;S, Double_t &amp;a, Double_t &amp;b, Double_t &amp;c, Matrix &amp;eigvec)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetInertiaTensor</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a8b1d672d1bf4b868f1087809f54a882e</anchor>
      <arglist>(System &amp;S, Double_t &amp;a, Double_t &amp;b, Double_t &amp;c, Matrix &amp;eigenvec, Matrix &amp;I)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Power.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Analysis/</path>
    <filename>Power_8cxx</filename>
    <includes id="NBody_8h" name="NBody.h" local="no" imported="no">NBody.h</includes>
    <namespace>NBody</namespace>
    <member kind="function">
      <type>void</type>
      <name>CalcPowSpectrum</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a8f7dff7ee7fc024bdf41ebe5ba9729d9</anchor>
      <arglist>(Double_t *rho, Double_t *p, int Ng, Double_t BoxLength)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DensityGridNGP</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a9e3bd604a9b342a16fae9e633c87a2a3</anchor>
      <arglist>(System &amp;S, Double_t *rho, int Ng, Double_t &amp;L)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Power.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Analysis/</path>
    <filename>Power_8h</filename>
    <includes id="NBody_8h" name="NBody.h" local="no" imported="no">NBody.h</includes>
    <namespace>NBody</namespace>
    <member kind="function">
      <type>void</type>
      <name>CalcPowSpectrum</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a8f7dff7ee7fc024bdf41ebe5ba9729d9</anchor>
      <arglist>(Double_t *rho, Double_t *p, int Ng, Double_t BoxLength)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DensityGridNGP</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a9e3bd604a9b342a16fae9e633c87a2a3</anchor>
      <arglist>(System &amp;S, Double_t *rho, int Ng, Double_t &amp;L)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Cosmology.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Cosmology/</path>
    <filename>Cosmology_8cxx</filename>
    <includes id="Cosmology_8h" name="Cosmology.h" local="no" imported="no">Cosmology.h</includes>
    <namespace>Cosmology</namespace>
    <member kind="function">
      <type>Double_t</type>
      <name>aHIntFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ad54a65fe4587f2e0a6b59bbce0fdc95e</anchor>
      <arglist>(Double_t a, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>aHIntFuncGSL</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a2f1cb4bb3d240599ed6633aa55858069</anchor>
      <arglist>(double a, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>aHIntFuncGSLMonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a8a3d0e95fcf0b95e41ac04b17d633acd</anchor>
      <arglist>(double a, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>aintegral</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>aca839619c098dab3d2d7e25e66f7a41d</anchor>
      <arglist>(const Double_t a1, const Double_t a2, const Double_t omega, const Double_t xlambda)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a269b55d49c4b8c5f81e22bb44798439c</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CEH98diff</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>adbcb8747ff665b7f0bacdfecce9adc63</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GrowthFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a701cbd7c6a984ed69d0256b9f4c00bd4</anchor>
      <arglist>(const Double_t a, const Double_t om, const Double_t ola)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GrowthFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a55eac79c3a605b928ec60eadd742d5e6</anchor>
      <arglist>(const Double_t a, const Double_t om, const Double_t ode, const Double_t ok, Double_t alpha, Double_t wde)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>HubbleFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>aad1be2f6641090f6757308631fb89e14</anchor>
      <arglist>(const Double_t a, const Double_t Ho, const Double_t om, const Double_t ola)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>HubbleFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a225aeef989bca704ce6b29bb20350640</anchor>
      <arglist>(const Double_t a, const Double_t Ho, const Double_t om, const Double_t ode, const Double_t ok, Double_t alpha, Double_t wde)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>HubbleIntFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a9fe760bfcd7471d2e11747ec4dfd6f68</anchor>
      <arglist>(Double_t a, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>HubbleIntFuncGSL</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a2bdb4d7670764893b4939cc8939fb356</anchor>
      <arglist>(double a, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>HubbleIntFuncGSLMonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>adec37f0116bfce0e36e634f8c4c7e135</anchor>
      <arglist>(double a, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegralkkPower</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a1db6a3f36e27864f82e0f77c38024393</anchor>
      <arglist>(int powerspectype, void *params, int integrationtype, bool logint)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegralSigma</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>afdcdec92d70adf7517fa8e192da28ef7</anchor>
      <arglist>(int powerspectype, Double_t R, void *params, int integrationtype, bool logint)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>LEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a04e6ad2e407977cc24ea195e9608344f</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>LEH98diff</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a112c703ccbb28974bc7177075ba2a0b2</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>neffEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a21f2aa6ccc0c33681417aed69d7f2a2c</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>neffWDMEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a2a2540ccb558a0585d833ec65cefb754</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Theta27, const Double_t Rd)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PBBKS</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>aad648da6daf2f7edb66289973eeedbbf</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Amp, const Double_t Omegam, const Double_t h, const Double_t Omegab)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PBBKSint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a15f2e2aafca0943851390548b82d3181</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PBBKSintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>af5e737552805567ad00e2f643476cbce</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PBBKSlogint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a59908f21fd15b932f1887f52cffa4e7b</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PBBKSlogintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a7c8948083a12484e110d9aa0342ac50b</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ab57708819b9a94327a84f9a8533447af</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t deltaH, const Double_t clight, const Double_t Ho, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a4dc841298fd503e20c695595c20df522</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Amp, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PEH98int</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a91da985b0dce318c620c8752da4f0867</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PEH98intgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a362daec612b3644d387bb0790bf216b7</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PEH98logint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a1fd1d0312797e78d0c81b412f711b5a7</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PEH98logintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ad48a21b0f47e026dff4fdec31807bf9f</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PGreen04</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a23ecc8f326d6040c773a6b5b2f2c1c14</anchor>
      <arglist>(const Double_t k, const Double_t Amp, const Double_t wm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PGreen04int</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a0a769ad4f78be1caabfddb824d69b634</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PGreen04intgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a749e1c3f250dc2e041012416f4acc3a7</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PGreen04logint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a96ee6732f3299f1bcc5a29c76122efee</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PGreen04logintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a619e6d7b3b9781419cc421a02c61d69d</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PK</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ad19289dc0b293f5610752d0eb217d9d8</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Amp)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PKint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a447874d1be6d0d5ad7c8b5250a041a58</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PKintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a7442fb1ee5be26ba69dbd8ac3b3d8ae5</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PKlogint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a576f1c4fb9280de51327a79e7e699b6a</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PKlogintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a7b1762a0283dee484cd5fa27c0e36442</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PWDMEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ac2dfabffb61c67897592a80cb8392c6d</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Amp, const Double_t Omegam, const Double_t h, const Double_t Theta27, const Double_t Rd)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PWDMEH98int</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a3b11af9a8a327f3276ddb85c36c8dc0a</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PWDMEH98intgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>adb48d2dabf29e33254083dc72d712392</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PWDMEH98logint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a62718d753953db960a604cdbd44a7612</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PWDMEH98logintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a5712d5ea7742344da26917ffd865a998</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>qBBSK</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a0eba524990c82b0326851b21db328a4d</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Omegab)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>qEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>aec718bb6cc3bb0bd15cfbe4e2a680fac</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>ScaleFactorFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a7bfec341af541d5b5108f74fc1087d33</anchor>
      <arglist>(const Double_t a, const Double_t om, const Double_t xla)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPBBKSint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a57479f00b1098f16bc9fa65683afa8e1</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPBBKSintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a34ae2e8704ba5160fe07648374e615a1</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPBBKSlogint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a3e4aff04bf9b71811838976106e03f36</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPBBKSlogintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a42163bd76913b101b18fe2fa8bba998d</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPEH98int</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a036565e1f1933766443b574c835caa6e</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPEH98intgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a1e32599583cbf58e2d6939edb1cef3b7</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPEH98logint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a587cb6c3c6e3b7666ef423a9c02fde33</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPEH98logintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a603facf330094e9ddc8bbddee9c16c71</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPGreen04int</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a27e7d4d22fead09731ccc9a49ba1091b</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPGreen04intgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ae150998c5eb14d2f3ae7ff6f2ebfc55c</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPGreen04logint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ab4b312eb3132c03d73ac85234a73d59d</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPGreen04logintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a91a6a037350785c0c0baf6d192a9d52e</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPKint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a47a816ef67cf310979aafb1fdf713cbd</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPKintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a69f5bc18907a19c212dede36182e5160</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPKlogint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a09dee238033d6b83fc8734e334530d51</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPKlogintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a093c1a07b23163c1e6809a1c40c9530e</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPWDMEH98int</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a384971b260d3ebb5cae2b0948fdbfe90</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPWDMEH98intgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a6a032be1299ad331a3a2b876abae086b</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPWDMEH98logint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a7c9c34451e8ea7fbac6eaa090d5fd6d3</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPWDMEH98logintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a5bd7e310476241c5e1100c6294cdc853</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>TBBKS</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a5368ed2e2aed7d2eb0d74c077f74dbb6</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Omegab)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>TEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a8d3f6583b853ba18cb080d36b62e133c</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Timeh</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a6f7793fa35350ec7c841b6b8a68da25a</anchor>
      <arglist>(const Double_t a, const Double_t om, const Double_t ola)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Timet</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ace1a8005ff6a3e92567df819694831b5</anchor>
      <arglist>(Double_t a, const Double_t om, const Double_t ola)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>WKR2</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>af72eeaf43aa0d79d371c0fb6492e3a45</anchor>
      <arglist>(const Double_t k, const Double_t R)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Cosmology.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Cosmology/</path>
    <filename>Cosmology_8h</filename>
    <includes id="NBody_8h" name="NBody.h" local="no" imported="no">NBody.h</includes>
    <includes id="NBodyMath_8h" name="NBodyMath.h" local="no" imported="no">NBodyMath.h</includes>
    <namespace>Cosmology</namespace>
    <member kind="function">
      <type>Double_t</type>
      <name>aintegral</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>aca839619c098dab3d2d7e25e66f7a41d</anchor>
      <arglist>(const Double_t a1, const Double_t a2, const Double_t omega, const Double_t xlambda)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a269b55d49c4b8c5f81e22bb44798439c</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CEH98diff</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>adbcb8747ff665b7f0bacdfecce9adc63</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GrowthFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>aea982866216250805ede1d530a632b5b</anchor>
      <arglist>(const Double_t a, const Double_t Ho, const Double_t om, const Double_t ode, const Double_t ok, Double_t alpha=0.0, Double_t wde=-1.0)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GrowthFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a701cbd7c6a984ed69d0256b9f4c00bd4</anchor>
      <arglist>(const Double_t a, const Double_t om, const Double_t ola)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>HubbleFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>aad1be2f6641090f6757308631fb89e14</anchor>
      <arglist>(const Double_t a, const Double_t Ho, const Double_t om, const Double_t ola)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>HubbleFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a225aeef989bca704ce6b29bb20350640</anchor>
      <arglist>(const Double_t a, const Double_t Ho, const Double_t om, const Double_t ode, const Double_t ok, Double_t alpha, Double_t wde)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegralkkPower</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a1db6a3f36e27864f82e0f77c38024393</anchor>
      <arglist>(int powerspectype, void *params, int integrationtype, bool logint)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegralSigma</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>afdcdec92d70adf7517fa8e192da28ef7</anchor>
      <arglist>(int powerspectype, Double_t R, void *params, int integrationtype, bool logint)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>LEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a04e6ad2e407977cc24ea195e9608344f</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>LEH98diff</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a112c703ccbb28974bc7177075ba2a0b2</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>neffEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a21f2aa6ccc0c33681417aed69d7f2a2c</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>neffWDMEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a2a2540ccb558a0585d833ec65cefb754</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Theta27, const Double_t Rd)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PBBKS</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>aad648da6daf2f7edb66289973eeedbbf</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Amp, const Double_t Omegam, const Double_t h, const Double_t Omegab)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ab57708819b9a94327a84f9a8533447af</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t deltaH, const Double_t clight, const Double_t Ho, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a4dc841298fd503e20c695595c20df522</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Amp, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PFromFile</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a449ab31b3b09ea16112fcdb1b33e7474</anchor>
      <arglist>(const Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PGreen04</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a23ecc8f326d6040c773a6b5b2f2c1c14</anchor>
      <arglist>(const Double_t k, const Double_t Amp, const Double_t wm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PK</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ad19289dc0b293f5610752d0eb217d9d8</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Amp)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PWDMEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ac2dfabffb61c67897592a80cb8392c6d</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Amp, const Double_t Omegam, const Double_t h, const Double_t Theta27, const Double_t Rd)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>qBBSK</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a0eba524990c82b0326851b21db328a4d</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Omegab)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>qEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>aec718bb6cc3bb0bd15cfbe4e2a680fac</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>ScaleFactorFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a7bfec341af541d5b5108f74fc1087d33</anchor>
      <arglist>(const Double_t a, const Double_t om, const Double_t xla)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>TBBKS</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a5368ed2e2aed7d2eb0d74c077f74dbb6</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Omegab)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>TEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a8d3f6583b853ba18cb080d36b62e133c</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Timeh</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a6f7793fa35350ec7c841b6b8a68da25a</anchor>
      <arglist>(const Double_t a, const Double_t om, const Double_t ola)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Timet</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ace1a8005ff6a3e92567df819694831b5</anchor>
      <arglist>(Double_t a, const Double_t om, const Double_t ola)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>WKR2</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>af72eeaf43aa0d79d371c0fb6492e3a45</anchor>
      <arglist>(const Double_t k, const Double_t R)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>InitCond.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/InitCond/</path>
    <filename>InitCond_8cxx</filename>
    <includes id="NBody_8h" name="NBody.h" local="no" imported="no">NBody.h</includes>
    <includes id="NBodyMath_8h" name="NBodyMath.h" local="no" imported="no">NBodyMath.h</includes>
    <includes id="InitCond_8h" name="InitCond.h" local="no" imported="no">InitCond.h</includes>
    <namespace>NBody</namespace>
    <member kind="function">
      <type>System *</type>
      <name>CubicGrid</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a503b86d87a69b1aa2e88056bcbd204b9</anchor>
      <arglist>(int N, Double_t BoxLength, Double_t Mtotal, Double_t time)</arglist>
    </member>
    <member kind="function">
      <type>System *</type>
      <name>UniformSphere</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a126d0532db237e67a2d2b881b8c3765d</anchor>
      <arglist>(int N, Double_t Radius, Double_t Mtotal, Double_t time)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>InitCond.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/InitCond/</path>
    <filename>InitCond_8h</filename>
    <includes id="NBody_8h" name="NBody.h" local="no" imported="no">NBody.h</includes>
    <includes id="NBodyMath_8h" name="NBodyMath.h" local="no" imported="no">NBodyMath.h</includes>
    <namespace>NBody</namespace>
    <member kind="function">
      <type>System *</type>
      <name>CubicGrid</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a503b86d87a69b1aa2e88056bcbd204b9</anchor>
      <arglist>(int N, Double_t BoxLength, Double_t Mtotal, Double_t time)</arglist>
    </member>
    <member kind="function">
      <type>System *</type>
      <name>UniformSphere</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a126d0532db237e67a2d2b881b8c3765d</anchor>
      <arglist>(int N, Double_t Radius, Double_t Mtotal, Double_t time)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>DistFunc.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/KDTree/</path>
    <filename>DistFunc_8h</filename>
    <includes id="NBodyMath_8h" name="NBodyMath.h" local="no" imported="no">NBodyMath.h</includes>
    <namespace>NBody</namespace>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a6d94b57122fcedb536c6b0c3885eca1e</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aec050d8b2be9712f3aac8386b639e4a5</anchor>
      <arglist>(const Double_t *v1, const Double_t *v2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>acfdf50f94d221815619e97831984154d</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a038a5a120237444305f0050dd3d408d4</anchor>
      <arglist>(const Double_t *v1, const Double_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PhaseDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ae693d0baeb7a12462799531c165418d0</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceProjSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a4536ce36d79eb5a3da958e19103d0fa1</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>afbe92863ed7180298c0e31faa8328198</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2, const Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a1e08c5a585ca16acdf190a1a05d08941</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricwithTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aa67cba5b7f8b27631fd824fd159a959c</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a0de40474be2bee286bd8a799f1fea32f</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>abe4c2dd5784b67ec2a7b0365bdae32bd</anchor>
      <arglist>(const Real_t *v1, const Real_t *v2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a87ca1c9704eba746948b1bdf92526c04</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a4c984ded3aff2bbf070daeb7dad067ea</anchor>
      <arglist>(const Real_t *v1, const Real_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PhaseDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a636f53a810a026a82738ecb198690ef9</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceProjSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a4c2b9f1bb732c80246da765970590ede</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a767e2056d87ce3e5b5b7d995cacc36d5</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2, const Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a99134f8221cd0922cf37c0ad07cb6304</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricwithTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>af2cf025d5cb8a4fa8f1fc935fd9f97c2</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a38e6b832ee8edf498f941e23f0b8e184</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a283d307a17468ba310500eb9cacb257c</anchor>
      <arglist>(const Double_t *v1, const Real_t *v2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a295f0d65cda2b9bcfa1fcfa5809210e9</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aceed47806f85296d2c8133fbaab9649e</anchor>
      <arglist>(const Double_t *v1, const Real_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PhaseDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a19e9ad32388807681ab4fb0d46ea2840</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceProjSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aa7b3973b7923f8086f3e3a0135395141</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a936161e0981e5aa3e6d1e26b8c2e012d</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2, const Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ae58dfacab58a84130ab85ba5a490b1f5</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricwithTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a8b8d97c6f5ba79d1d5b1016df95b9ec2</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a37c40a63a326927d2ed6689c4b9f643e</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a1a36b17b45103ac1420e06fba4d75624</anchor>
      <arglist>(const Real_t *v1, const Double_t *v2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a9f1247c5c83dbb5399095d44e4d0f1e1</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a438c6f4d2d62774d1411a844564925d3</anchor>
      <arglist>(const Real_t *v1, const Double_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PhaseDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a3ca04510e1122e8475d23fde1a96ead1</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceProjSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a3c54b28ff5ddbf31aa0e473bd34db24e</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aa719edd1c96a433cff1d5389bdc9dd1c</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2, const Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a50d4baedaab96f71c4d8e44170db0571</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricwithTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ad223d7bdc6917b8dc13f95331b99432d</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflection1D</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>acfcca7d0b1e4c399e4e01cfa6062d658</anchor>
      <arglist>(const Coordinate &amp;x0, Coordinate &amp;xp, Coordinate p, int k)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflection2D</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a53422b8fc7e733ec5c610281baaae7e2</anchor>
      <arglist>(const Coordinate &amp;x0, Coordinate &amp;xp, Coordinate p, int k1, int k2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflectionND</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a963d7d6b33ec1c4a582474e41e21bbcb</anchor>
      <arglist>(const Coordinate &amp;x0, Coordinate &amp;xp, Coordinate p, int ndim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflection1D</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a7424acea8a5356b27372a9578e9a155c</anchor>
      <arglist>(const Particle &amp;x0, Particle &amp;xp, Coordinate p, int k)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflection2D</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a98bd2d2266c8b4a93dcb130488884979</anchor>
      <arglist>(const Particle &amp;x0, Particle &amp;xp, Coordinate p, int k1, int k2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflectionND</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>af5007ee2e8aa728826863caa82d70f53</anchor>
      <arglist>(const Particle &amp;x0, Particle &amp;xp, Coordinate p, int ndim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a6d94b57122fcedb536c6b0c3885eca1e</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aec050d8b2be9712f3aac8386b639e4a5</anchor>
      <arglist>(const Double_t *v1, const Double_t *v2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>acfdf50f94d221815619e97831984154d</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a038a5a120237444305f0050dd3d408d4</anchor>
      <arglist>(const Double_t *v1, const Double_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PhaseDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ae693d0baeb7a12462799531c165418d0</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceProjSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a4536ce36d79eb5a3da958e19103d0fa1</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>afbe92863ed7180298c0e31faa8328198</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2, const Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a1e08c5a585ca16acdf190a1a05d08941</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricwithTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aa67cba5b7f8b27631fd824fd159a959c</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a0de40474be2bee286bd8a799f1fea32f</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>abe4c2dd5784b67ec2a7b0365bdae32bd</anchor>
      <arglist>(const Real_t *v1, const Real_t *v2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a87ca1c9704eba746948b1bdf92526c04</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a4c984ded3aff2bbf070daeb7dad067ea</anchor>
      <arglist>(const Real_t *v1, const Real_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PhaseDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a636f53a810a026a82738ecb198690ef9</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceProjSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a4c2b9f1bb732c80246da765970590ede</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a767e2056d87ce3e5b5b7d995cacc36d5</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2, const Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a99134f8221cd0922cf37c0ad07cb6304</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricwithTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>af2cf025d5cb8a4fa8f1fc935fd9f97c2</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a38e6b832ee8edf498f941e23f0b8e184</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a283d307a17468ba310500eb9cacb257c</anchor>
      <arglist>(const Double_t *v1, const Real_t *v2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a295f0d65cda2b9bcfa1fcfa5809210e9</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aceed47806f85296d2c8133fbaab9649e</anchor>
      <arglist>(const Double_t *v1, const Real_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PhaseDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a19e9ad32388807681ab4fb0d46ea2840</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceProjSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aa7b3973b7923f8086f3e3a0135395141</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a936161e0981e5aa3e6d1e26b8c2e012d</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2, const Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ae58dfacab58a84130ab85ba5a490b1f5</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricwithTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a8b8d97c6f5ba79d1d5b1016df95b9ec2</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a37c40a63a326927d2ed6689c4b9f643e</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a1a36b17b45103ac1420e06fba4d75624</anchor>
      <arglist>(const Real_t *v1, const Double_t *v2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a9f1247c5c83dbb5399095d44e4d0f1e1</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a438c6f4d2d62774d1411a844564925d3</anchor>
      <arglist>(const Real_t *v1, const Double_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PhaseDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a3ca04510e1122e8475d23fde1a96ead1</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceProjSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a3c54b28ff5ddbf31aa0e473bd34db24e</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aa719edd1c96a433cff1d5389bdc9dd1c</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2, const Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a50d4baedaab96f71c4d8e44170db0571</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricwithTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ad223d7bdc6917b8dc13f95331b99432d</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflection1D</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>acfcca7d0b1e4c399e4e01cfa6062d658</anchor>
      <arglist>(const Coordinate &amp;x0, Coordinate &amp;xp, Coordinate p, int k)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflection2D</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a53422b8fc7e733ec5c610281baaae7e2</anchor>
      <arglist>(const Coordinate &amp;x0, Coordinate &amp;xp, Coordinate p, int k1, int k2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflectionND</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a963d7d6b33ec1c4a582474e41e21bbcb</anchor>
      <arglist>(const Coordinate &amp;x0, Coordinate &amp;xp, Coordinate p, int ndim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflection1D</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a7424acea8a5356b27372a9578e9a155c</anchor>
      <arglist>(const Particle &amp;x0, Particle &amp;xp, Coordinate p, int k)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflection2D</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a98bd2d2266c8b4a93dcb130488884979</anchor>
      <arglist>(const Particle &amp;x0, Particle &amp;xp, Coordinate p, int k1, int k2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflectionND</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>af5007ee2e8aa728826863caa82d70f53</anchor>
      <arglist>(const Particle &amp;x0, Particle &amp;xp, Coordinate p, int ndim)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>FOFFunc.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/KDTree/</path>
    <filename>FOFFunc_8h</filename>
    <namespace>NBody</namespace>
    <member kind="typedef">
      <type>int(*</type>
      <name>FOFcheckfunc</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ae68c57f56f1373303562475da97d3763</anchor>
      <arglist>)(Particle &amp;, Double_t *)</arglist>
    </member>
    <member kind="typedef">
      <type>int(*</type>
      <name>FOFcompfunc</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a4b80079c12e359efc6b9d5d65b469aa9</anchor>
      <arglist>)(Particle &amp;, Particle &amp;, Double_t *)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>FOF3d</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a793b2e06c4f5156974cab32d93177501</anchor>
      <arglist>(Particle &amp;a, Particle &amp;b, Double_t *params)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>FOF6d</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a60efce1bf35d4f2c094927fcfdcd09f3</anchor>
      <arglist>(Particle &amp;a, Particle &amp;b, Double_t *params)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>FOFVel</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a6f1486955ca76c8b4df2a3b8b9196b2e</anchor>
      <arglist>(Particle &amp;a, Particle &amp;b, Double_t *params)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>Pnocheck</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aa88761eca5c73784bcb9298e54bd4fd6</anchor>
      <arglist>(Particle &amp;a, Double_t *params)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>KDCalcSmoothQuantities.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/KDTree/</path>
    <filename>KDCalcSmoothQuantities_8cxx</filename>
    <includes id="KDTree_8h" name="KDTree.h" local="no" imported="no">KDTree.h</includes>
    <namespace>NBody</namespace>
  </compound>
  <compound kind="file">
    <name>KDFindNearest.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/KDTree/</path>
    <filename>KDFindNearest_8cxx</filename>
    <includes id="KDTree_8h" name="KDTree.h" local="no" imported="no">KDTree.h</includes>
    <namespace>NBody</namespace>
  </compound>
  <compound kind="file">
    <name>KDFOF.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/KDTree/</path>
    <filename>KDFOF_8cxx</filename>
    <includes id="KDTree_8h" name="KDTree.h" local="no" imported="no">KDTree.h</includes>
    <namespace>NBody</namespace>
  </compound>
  <compound kind="file">
    <name>KDLeafNode.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/KDTree/</path>
    <filename>KDLeafNode_8cxx</filename>
    <includes id="KDNode_8h" name="KDNode.h" local="no" imported="no">KDNode.h</includes>
    <namespace>NBody</namespace>
  </compound>
  <compound kind="file">
    <name>KDNode.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/KDTree/</path>
    <filename>KDNode_8h</filename>
    <includes id="NBody_8h" name="NBody.h" local="no" imported="no">NBody.h</includes>
    <includes id="NBodyMath_8h" name="NBodyMath.h" local="no" imported="no">NBodyMath.h</includes>
    <includes id="PriorityQueue_8h" name="PriorityQueue.h" local="no" imported="no">PriorityQueue.h</includes>
    <includes id="DistFunc_8h" name="DistFunc.h" local="no" imported="no">DistFunc.h</includes>
    <includes id="FOFFunc_8h" name="FOFFunc.h" local="no" imported="no">FOFFunc.h</includes>
    <class kind="class">NBody::LeafNode</class>
    <class kind="class">NBody::Node</class>
    <class kind="class">NBody::SplitNode</class>
    <namespace>NBody</namespace>
    <member kind="define">
      <type>#define</type>
      <name>NPHASEDIM</name>
      <anchorfile>KDNode_8h.html</anchorfile>
      <anchor>aa0f9d08ebd731cd3e7f50393285f9382</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>NSPACEDIM</name>
      <anchorfile>KDNode_8h.html</anchorfile>
      <anchor>ab0927bea92a0c45ff803d07ea0881f53</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>int</type>
      <name>Int_tree_t</name>
      <anchorfile>KDNode_8h.html</anchorfile>
      <anchor>a7e4452867d3aa2953e75247d91150bf9</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>unsigned int</type>
      <name>UInt_tree_t</name>
      <anchorfile>KDNode_8h.html</anchorfile>
      <anchor>aa069d5cc66dc8c6160dd5c0cc66d766c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>KDSplitNode.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/KDTree/</path>
    <filename>KDSplitNode_8cxx</filename>
    <includes id="KDNode_8h" name="KDNode.h" local="no" imported="no">KDNode.h</includes>
    <namespace>NBody</namespace>
  </compound>
  <compound kind="file">
    <name>KDTree.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/KDTree/</path>
    <filename>KDTree_8cxx</filename>
    <includes id="KDTree_8h" name="KDTree.h" local="no" imported="no">KDTree.h</includes>
    <namespace>NBody</namespace>
  </compound>
  <compound kind="file">
    <name>KDTree.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/KDTree/</path>
    <filename>KDTree_8h</filename>
    <includes id="NBody_8h" name="NBody.h" local="no" imported="no">NBody.h</includes>
    <includes id="PriorityQueue_8h" name="PriorityQueue.h" local="no" imported="no">PriorityQueue.h</includes>
    <includes id="KDNode_8h" name="KDNode.h" local="no" imported="no">KDNode.h</includes>
    <includes id="NBodyMath_8h" name="NBodyMath.h" local="no" imported="no">NBodyMath.h</includes>
    <includes id="SmoothingKernels_8h" name="SmoothingKernels.h" local="no" imported="no">SmoothingKernels.h</includes>
    <includes id="FOFFunc_8h" name="FOFFunc.h" local="no" imported="no">FOFFunc.h</includes>
    <class kind="class">NBody::KDTree</class>
    <namespace>NBody</namespace>
  </compound>
  <compound kind="file">
    <name>PriorityQueue.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/KDTree/</path>
    <filename>PriorityQueue_8h</filename>
    <includes id="NBodyMath_8h" name="NBodyMath.h" local="no" imported="no">NBodyMath.h</includes>
    <class kind="class">NBody::NPriorityQueue</class>
    <class kind="struct">NBody::NPriorityQueue::nqueue</class>
    <class kind="class">NBody::PriorityQueue</class>
    <class kind="struct">NBody::PriorityQueue::queue_member</class>
    <namespace>NBody</namespace>
  </compound>
  <compound kind="file">
    <name>SmoothingKernels.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/KDTree/</path>
    <filename>SmoothingKernels_8h</filename>
    <class kind="struct">NBody::smoothfunc</class>
    <namespace>NBody</namespace>
    <member kind="function">
      <type>Double_t</type>
      <name>WEpan</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a1f1b0b18be2cc3c6ec759e8dfc92f806</anchor>
      <arglist>(Double_t r, Double_t h)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>WGauss</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ac68030ab7c00989fe9e5666dfb61033e</anchor>
      <arglist>(Double_t r, Double_t h)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>WSPH</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ac9bcae63f04efc56464010976d84b580</anchor>
      <arglist>(Double_t r, Double_t h)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>WTH</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a073d954241c808cc19effdd73d2d1ebd</anchor>
      <arglist>(Double_t r, Double_t h)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Coordinate.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Coordinate_8cxx</filename>
    <includes id="Coordinate_8h" name="Coordinate.h" local="no" imported="no">Coordinate.h</includes>
    <namespace>Math</namespace>
    <member kind="function">
      <type>Coordinate</type>
      <name>operator*</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a899ae08543e7360b976e44a7b69d9a1d</anchor>
      <arglist>(Double_t a, const Coordinate &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>ac3d40a22a9698f7c68a2018c0971eb02</anchor>
      <arglist>(std::ostream &amp;outs, Coordinate c)</arglist>
    </member>
    <member kind="function">
      <type>std::istream &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>acc7c91fc21d8fd273fbf783e8ea29727</anchor>
      <arglist>(std::istream &amp;ins, Coordinate &amp;c)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Coordinate.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Coordinate_8h</filename>
    <includes id="Precision_8h" name="Precision.h" local="no" imported="no">Precision.h</includes>
    <class kind="struct">Math::Coord</class>
    <class kind="class">Math::Coordinate</class>
    <namespace>Math</namespace>
  </compound>
  <compound kind="file">
    <name>Coordinate2D.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Coordinate2D_8cxx</filename>
    <includes id="Coordinate2D_8h" name="Coordinate2D.h" local="no" imported="no">Coordinate2D.h</includes>
    <namespace>Math</namespace>
    <member kind="function">
      <type>Coordinate2D</type>
      <name>operator*</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>aa998c693107b8b5538e0ef2bcab79498</anchor>
      <arglist>(Double_t a, const Coordinate2D &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a9ca86e51f9d04aa5d8c0e23b81c6141d</anchor>
      <arglist>(std::ostream &amp;outs, Coordinate2D c)</arglist>
    </member>
    <member kind="function">
      <type>std::istream &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a03d140240a639aa68092a221a2f9aa4c</anchor>
      <arglist>(std::istream &amp;ins, Coordinate2D &amp;c)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Coordinate2D.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Coordinate2D_8h</filename>
    <includes id="Precision_8h" name="Precision.h" local="no" imported="no">Precision.h</includes>
    <class kind="class">Math::Coordinate2D</class>
    <namespace>Math</namespace>
  </compound>
  <compound kind="file">
    <name>ExtremaRootFinding.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>ExtremaRootFinding_8h</filename>
    <includes id="Precision_8h" name="Precision.h" local="no" imported="no">Precision.h</includes>
    <includes id="Function_8h" name="Function.h" local="no" imported="no">Function.h</includes>
    <includes id="Matrix_8h" name="Matrix.h" local="no" imported="no">Matrix.h</includes>
    <includes id="GMatrix_8h" name="GMatrix.h" local="no" imported="no">GMatrix.h</includes>
    <includes id="Interpolate_8h" name="Interpolate.h" local="no" imported="no">Interpolate.h</includes>
    <namespace>Math</namespace>
  </compound>
  <compound kind="file">
    <name>Fitting.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Fitting_8cxx</filename>
    <includes id="Fitting_8h" name="Fitting.h" local="no" imported="no">Fitting.h</includes>
    <namespace>Math</namespace>
    <member kind="function">
      <type>Double_t</type>
      <name>FitNonLinLS</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a7d7cb8f24f0499e5a3ec12286485aba6</anchor>
      <arglist>(const math_function fitfunc, const math_function *difffuncs, const int nparams, Double_t *params, GMatrix &amp;covar, const int npoints, const Double_t x[], const Double_t y[], GMatrix *Weight=NULL, Double_t error=1e-3, Double_t cl=0.95, int *fixparam=NULL, int binned=1, int maxiterations=1000, int iestimateerror=0)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>OptimalBins</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>af1c80002b9862fc16637882937b62664</anchor>
      <arglist>(const int npoints, Double_t *points, Double_t xmin, Double_t xmax, Double_t *weights=NULL)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Fitting.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Fitting_8h</filename>
    <includes id="Precision_8h" name="Precision.h" local="no" imported="no">Precision.h</includes>
    <includes id="Function_8h" name="Function.h" local="no" imported="no">Function.h</includes>
    <includes id="Matrix_8h" name="Matrix.h" local="no" imported="no">Matrix.h</includes>
    <includes id="GMatrix_8h" name="GMatrix.h" local="no" imported="no">GMatrix.h</includes>
    <includes id="Interpolate_8h" name="Interpolate.h" local="no" imported="no">Interpolate.h</includes>
    <includes id="ExtremaRootFinding_8h" name="ExtremaRootFinding.h" local="no" imported="no">ExtremaRootFinding.h</includes>
    <namespace>Math</namespace>
    <member kind="function">
      <type>Double_t</type>
      <name>FitNonLinLS</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a7d7cb8f24f0499e5a3ec12286485aba6</anchor>
      <arglist>(const math_function fitfunc, const math_function *difffuncs, const int nparams, Double_t *params, GMatrix &amp;covar, const int npoints, const Double_t x[], const Double_t y[], GMatrix *Weight=NULL, Double_t error=1e-3, Double_t cl=0.95, int *fixparam=NULL, int binned=1, int maxiterations=1000, int iestimateerror=0)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>OptimalBins</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>af1c80002b9862fc16637882937b62664</anchor>
      <arglist>(const int npoints, Double_t *points, Double_t xmin, Double_t xmax, Double_t *weights=NULL)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Function.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Function_8h</filename>
    <class kind="struct">Math::math_function</class>
    <class kind="struct">Math::math_multidim_function</class>
    <namespace>Math</namespace>
  </compound>
  <compound kind="file">
    <name>GMatrix.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>GMatrix_8cxx</filename>
    <includes id="GMatrix_8h" name="GMatrix.h" local="no" imported="no">GMatrix.h</includes>
    <namespace>Math</namespace>
  </compound>
  <compound kind="file">
    <name>GMatrix.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>GMatrix_8h</filename>
    <includes id="Precision_8h" name="Precision.h" local="no" imported="no">Precision.h</includes>
    <includes id="Coordinate_8h" name="Coordinate.h" local="no" imported="no">Coordinate.h</includes>
    <includes id="Matrix_8h" name="Matrix.h" local="no" imported="no">Matrix.h</includes>
    <class kind="class">Math::GMatrix</class>
    <namespace>Math</namespace>
  </compound>
  <compound kind="file">
    <name>Integrate.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Integrate_8cxx</filename>
    <includes id="Integrate_8h" name="Integrate.h" local="no" imported="no">Integrate.h</includes>
    <namespace>Math</namespace>
    <member kind="define">
      <type>#define</type>
      <name>ALPH</name>
      <anchorfile>Integrate_8cxx.html</anchorfile>
      <anchor>aef36de74f149c1c490a11183dec82801</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MXDIM</name>
      <anchorfile>Integrate_8cxx.html</anchorfile>
      <anchor>a6540625776b644130c080d90967c7863</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>NDMX</name>
      <anchorfile>Integrate_8cxx.html</anchorfile>
      <anchor>a01449e9b435a7ae91a05f1ff9b7709fc</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>TINY</name>
      <anchorfile>Integrate_8cxx.html</anchorfile>
      <anchor>acf1c38f71f39386356edb151a131ad11</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateClosed</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a9426b041fe580caa17e4dae6b22f18ab</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b, const int n)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateData</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a7b662bfcf4a1b7da754e282578226742</anchor>
      <arglist>(const Double_t x[], const Double_t f[], const int lower, const int upper)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateQTrap</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>ae585139c1b1184047125bc4109072549</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b, const Double_t epsrel, int IMAX)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateRomberg</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a8654394ccfc022538a69c7af8e2dbd7f</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b, const Double_t epsrel, const int n, const int K)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateRombergMonte</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>ad462b2a23bfaa735abcb41fd97de868c</anchor>
      <arglist>(math_multidim_function *f, Double_t *a, Double_t *b, const Double_t epsrel, const int sampling, const int n, const int K, int interations=6)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateSimpleTrapezoidal</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a4ce1ced0a7bd11555915dc52bac4e87c</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b, const int n)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateSimpson</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>aa8d1091f877e020438189780fabc1448</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateTrap</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a5e5e5c79213f576436b082d589169dd9</anchor>
      <arglist>(const Double_t x[], const Double_t f[], const int a, const int b)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateTrapezoidal</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>af28d827090a237e07c9cf82bd2a68b11</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b, const int n)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateVegasMonte</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>aa9ba450da02048058df4c70fee238405</anchor>
      <arglist>(math_multidim_function *fm, Double_t *a, Double_t *b, const int numintervals, Double_t *error, Double_t *chisq, int interations)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>IntegrateVegasMonte</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>aae99a70a37aa9c1848ff7f972d00d9e3</anchor>
      <arglist>(gsl_monte_function *gslfm, double *a, double *b, const int numintervals, double *error, double *chisq)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rebin</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>aca549ace0630f4eee107187d3e5ccb41</anchor>
      <arglist>(Double_t rc, int nd, Double_t r[], Double_t xin[], Double_t xi[])</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>vegas</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a7f438f7e967b6c3d9cfe88d8b370a5ee</anchor>
      <arglist>(math_multidim_function *fxn, Double_t regn[], int ndim, int init, unsigned long ncall, int itmx, int nprn, long int idum, Double_t *tgral, Double_t *sd, Double_t *chi2a)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Integrate.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Integrate_8h</filename>
    <includes id="Precision_8h" name="Precision.h" local="no" imported="no">Precision.h</includes>
    <includes id="Function_8h" name="Function.h" local="no" imported="no">Function.h</includes>
    <includes id="Interpolate_8h" name="Interpolate.h" local="no" imported="no">Interpolate.h</includes>
    <includes id="Random_8h" name="Random.h" local="no" imported="no">Random.h</includes>
    <namespace>Math</namespace>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateClosed</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a9426b041fe580caa17e4dae6b22f18ab</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b, const int n)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateData</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a7b662bfcf4a1b7da754e282578226742</anchor>
      <arglist>(const Double_t x[], const Double_t f[], const int lower, const int upper)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateQTrap</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>ae585139c1b1184047125bc4109072549</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b, const Double_t epsrel, int IMAX)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateRomberg</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a8654394ccfc022538a69c7af8e2dbd7f</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b, const Double_t epsrel, const int n, const int K)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateRombergMonte</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>ad462b2a23bfaa735abcb41fd97de868c</anchor>
      <arglist>(math_multidim_function *f, Double_t *a, Double_t *b, const Double_t epsrel, const int sampling, const int n, const int K, int interations=6)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateSimpleTrapezoidal</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a4ce1ced0a7bd11555915dc52bac4e87c</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b, const int n)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateSimpson</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>aa8d1091f877e020438189780fabc1448</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateTrap</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a5e5e5c79213f576436b082d589169dd9</anchor>
      <arglist>(const Double_t x[], const Double_t f[], const int a, const int b)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateTrapezoidal</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>af28d827090a237e07c9cf82bd2a68b11</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b, const int n)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateVegasMonte</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>aa9ba450da02048058df4c70fee238405</anchor>
      <arglist>(math_multidim_function *fm, Double_t *a, Double_t *b, const int numintervals, Double_t *error, Double_t *chisq, int interations)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>IntegrateVegasMonte</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>aae99a70a37aa9c1848ff7f972d00d9e3</anchor>
      <arglist>(gsl_monte_function *gslfm, double *a, double *b, const int numintervals, double *error, double *chisq)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rebin</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>aca549ace0630f4eee107187d3e5ccb41</anchor>
      <arglist>(Double_t rc, int nd, Double_t r[], Double_t xin[], Double_t xi[])</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>vegas</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a7f438f7e967b6c3d9cfe88d8b370a5ee</anchor>
      <arglist>(math_multidim_function *fxn, Double_t regn[], int ndim, int init, unsigned long ncall, int itmx, int nprn, long int idum, Double_t *tgral, Double_t *sd, Double_t *chi2a)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Interpolate.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Interpolate_8cxx</filename>
    <includes id="Interpolate_8h" name="Interpolate.h" local="no" imported="no">Interpolate.h</includes>
    <namespace>Math</namespace>
    <member kind="function">
      <type>void</type>
      <name>PolyInt</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a261cfcce47b95f62b4cb4a77e4cedeea</anchor>
      <arglist>(Double_t *xa, Double_t *ya, int n, Double_t x, Double_t &amp;y, Double_t &amp;dy)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Interpolate.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Interpolate_8h</filename>
    <includes id="NBodyMath_8h" name="NBodyMath.h" local="no" imported="no">NBodyMath.h</includes>
    <includes id="Precision_8h" name="Precision.h" local="no" imported="no">Precision.h</includes>
    <namespace>Math</namespace>
    <member kind="function">
      <type>void</type>
      <name>PolyInt</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a261cfcce47b95f62b4cb4a77e4cedeea</anchor>
      <arglist>(Double_t *xa, Double_t *ya, int n, Double_t x, Double_t &amp;y, Double_t &amp;dy)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Matrix.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Matrix_8cxx</filename>
    <includes id="Coordinate_8h" name="Coordinate.h" local="no" imported="no">Coordinate.h</includes>
    <includes id="Matrix_8h" name="Matrix.h" local="no" imported="no">Matrix.h</includes>
    <namespace>Math</namespace>
    <member kind="function">
      <type>Matrix</type>
      <name>operator*</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>af26e1bc6214f1a7e3a1164709c764cb8</anchor>
      <arglist>(Double_t a, const Matrix &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a7d33cb24d8a9027be45c2b933736118d</anchor>
      <arglist>(std::ostream &amp;stream, const Matrix &amp;m)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Matrix.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Matrix_8h</filename>
    <includes id="Precision_8h" name="Precision.h" local="no" imported="no">Precision.h</includes>
    <includes id="Coordinate_8h" name="Coordinate.h" local="no" imported="no">Coordinate.h</includes>
    <class kind="class">Math::Matrix</class>
    <namespace>Math</namespace>
  </compound>
  <compound kind="file">
    <name>Matrix2D.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Matrix2D_8cxx</filename>
    <includes id="Coordinate2D_8h" name="Coordinate2D.h" local="no" imported="no">Coordinate2D.h</includes>
    <includes id="Matrix2D_8h" name="Matrix2D.h" local="no" imported="no">Matrix2D.h</includes>
    <namespace>Math</namespace>
    <member kind="function">
      <type>Matrix2D</type>
      <name>operator*</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>abf0582063e3317432a5c627cbc0b5c3e</anchor>
      <arglist>(Double_t a, const Matrix2D &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a52310bb38b8cd3769957842c350fd94f</anchor>
      <arglist>(std::ostream &amp;stream, const Matrix2D &amp;m)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Matrix2D.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Matrix2D_8h</filename>
    <includes id="Precision_8h" name="Precision.h" local="no" imported="no">Precision.h</includes>
    <includes id="Coordinate2D_8h" name="Coordinate2D.h" local="no" imported="no">Coordinate2D.h</includes>
    <class kind="class">Math::Matrix2D</class>
    <namespace>Math</namespace>
  </compound>
  <compound kind="file">
    <name>NBodyMath.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>NBodyMath_8h</filename>
    <includes id="Precision_8h" name="Precision.h" local="no" imported="no">Precision.h</includes>
    <includes id="Coordinate_8h" name="Coordinate.h" local="no" imported="no">Coordinate.h</includes>
    <includes id="Coordinate2D_8h" name="Coordinate2D.h" local="no" imported="no">Coordinate2D.h</includes>
    <includes id="Matrix2D_8h" name="Matrix2D.h" local="no" imported="no">Matrix2D.h</includes>
    <includes id="Matrix_8h" name="Matrix.h" local="no" imported="no">Matrix.h</includes>
    <includes id="GMatrix_8h" name="GMatrix.h" local="no" imported="no">GMatrix.h</includes>
    <includes id="Function_8h" name="Function.h" local="no" imported="no">Function.h</includes>
    <includes id="Integrate_8h" name="Integrate.h" local="no" imported="no">Integrate.h</includes>
    <includes id="Interpolate_8h" name="Interpolate.h" local="no" imported="no">Interpolate.h</includes>
    <includes id="Fitting_8h" name="Fitting.h" local="no" imported="no">Fitting.h</includes>
    <includes id="SpecialFunctions_8h" name="SpecialFunctions.h" local="no" imported="no">SpecialFunctions.h</includes>
    <includes id="Random_8h" name="Random.h" local="no" imported="no">Random.h</includes>
    <includes id="Statistics_8h" name="Statistics.h" local="no" imported="no">Statistics.h</includes>
    <namespace>Math</namespace>
  </compound>
  <compound kind="file">
    <name>Precision.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Precision_8h</filename>
    <member kind="define">
      <type>#define</type>
      <name>MAXVALUE</name>
      <anchorfile>Precision_8h.html</anchorfile>
      <anchor>a28c7f4454d27d9a294ebe7eb12a54e8c</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MINVALUE</name>
      <anchorfile>Precision_8h.html</anchorfile>
      <anchor>ae6be43e86a075f1cec470ab1ec1fe3c3</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>double</type>
      <name>Double_t</name>
      <anchorfile>Precision_8h.html</anchorfile>
      <anchor>a72d9b9fcbdedda8192c902d5244cff12</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>int</type>
      <name>Int_t</name>
      <anchorfile>Precision_8h.html</anchorfile>
      <anchor>a5dcf592183fe647745a2741671e58b52</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>float</type>
      <name>Real_t</name>
      <anchorfile>Precision_8h.html</anchorfile>
      <anchor>a5657d17fda4eb7b338c46648b104a805</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>unsigned int</type>
      <name>UInt_t</name>
      <anchorfile>Precision_8h.html</anchorfile>
      <anchor>a7c1bc4939263cb6c0f48e434d77ac258</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Random.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Random_8cxx</filename>
    <includes id="Random_8h" name="Random.h" local="no" imported="no">Random.h</includes>
    <namespace>Math</namespace>
    <member kind="define">
      <type>#define</type>
      <name>AM</name>
      <anchorfile>Random_8cxx.html</anchorfile>
      <anchor>ad301e6a88b1c01108f4867f2ea6f683c</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>EPS</name>
      <anchorfile>Random_8cxx.html</anchorfile>
      <anchor>a6ebf6899d6c1c8b7b9d09be872c05aae</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IA1</name>
      <anchorfile>Random_8cxx.html</anchorfile>
      <anchor>a6ef2749dca39c605c3d033f788afe6e3</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IA2</name>
      <anchorfile>Random_8cxx.html</anchorfile>
      <anchor>a372a58d7e9e25912fd79e7afaa06cc7a</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IM1</name>
      <anchorfile>Random_8cxx.html</anchorfile>
      <anchor>a78325bdf48423acef7c012567628b391</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IM2</name>
      <anchorfile>Random_8cxx.html</anchorfile>
      <anchor>ad8c519de7e5de4ae35344ddcf21fd062</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IMM1</name>
      <anchorfile>Random_8cxx.html</anchorfile>
      <anchor>a87a6e0054f9d827c979c43aa0d5e621a</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IQ1</name>
      <anchorfile>Random_8cxx.html</anchorfile>
      <anchor>a9afa86ff22da69bda72fe271ae71ad46</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IQ2</name>
      <anchorfile>Random_8cxx.html</anchorfile>
      <anchor>abae040385946a6acdf5e10d2efd86f3d</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IR1</name>
      <anchorfile>Random_8cxx.html</anchorfile>
      <anchor>a7b2b32709f9770a283701ffcf3723497</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>IR2</name>
      <anchorfile>Random_8cxx.html</anchorfile>
      <anchor>a5a08e4f5cb3582e623cc14a6c92d48de</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>NDIV</name>
      <anchorfile>Random_8cxx.html</anchorfile>
      <anchor>a62339d74dd5d9d00480e1a288cf88fe8</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>NTAB</name>
      <anchorfile>Random_8cxx.html</anchorfile>
      <anchor>a0e93cfb2d62849853fd34957ba6e6fdc</anchor>
      <arglist></arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>RNMX</name>
      <anchorfile>Random_8cxx.html</anchorfile>
      <anchor>aa7436c9270ffb06f8c1eae8d2e605cec</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>nran2</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a431cdd9a33c8aa71fb9062287782c9e3</anchor>
      <arglist>(long *idum)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>ran2</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a9e6edeee1c9bf8200968ff55ce476a71</anchor>
      <arglist>(long *idum)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Random.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Random_8h</filename>
    <includes id="Precision_8h" name="Precision.h" local="no" imported="no">Precision.h</includes>
    <namespace>Math</namespace>
    <member kind="function">
      <type>Double_t</type>
      <name>nran2</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a431cdd9a33c8aa71fb9062287782c9e3</anchor>
      <arglist>(long *idum)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>ran2</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a9e6edeee1c9bf8200968ff55ce476a71</anchor>
      <arglist>(long *idum)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>SpecialFunctions.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>SpecialFunctions_8cxx</filename>
    <includes id="SpecialFunctions_8h" name="SpecialFunctions.h" local="no" imported="no">SpecialFunctions.h</includes>
    <namespace>Math</namespace>
    <member kind="function">
      <type>int</type>
      <name>Factorial</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a2d0f9a58a07490b72ada749b89751631</anchor>
      <arglist>(int k)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>SpecialFunctions.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>SpecialFunctions_8h</filename>
    <includes id="Precision_8h" name="Precision.h" local="no" imported="no">Precision.h</includes>
    <namespace>Math</namespace>
    <member kind="function">
      <type>int</type>
      <name>Factorial</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a2d0f9a58a07490b72ada749b89751631</anchor>
      <arglist>(int k)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Statistics.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Statistics_8cxx</filename>
    <includes id="Statistics_8h" name="Statistics.h" local="no" imported="no">Statistics.h</includes>
    <namespace>Math</namespace>
  </compound>
  <compound kind="file">
    <name>Statistics.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/Math/</path>
    <filename>Statistics_8h</filename>
    <includes id="Precision_8h" name="Precision.h" local="no" imported="no">Precision.h</includes>
    <includes id="Coordinate_8h" name="Coordinate.h" local="no" imported="no">Coordinate.h</includes>
    <includes id="Matrix_8h" name="Matrix.h" local="no" imported="no">Matrix.h</includes>
    <includes id="GMatrix_8h" name="GMatrix.h" local="no" imported="no">GMatrix.h</includes>
    <namespace>Math</namespace>
    <member kind="function">
      <type>Double_t</type>
      <name>average</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a16bfec5be9cdda9d8b961158b64c22cf</anchor>
      <arglist>(Int_t n, Double_t *x)</arglist>
    </member>
    <member kind="function">
      <type>Double_t *</type>
      <name>average</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a30303faecc8bb40e09a853e032aef6a3</anchor>
      <arglist>(Int_t parnum, Int_t datnum, Double_t **x)</arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>covariance</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a8187a177d6158b91c01027c46da62a84</anchor>
      <arglist>(Int_t parnum, Int_t datnum, Double_t *xm, Double_t **x)</arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>covariance</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a03718a69fdf5fdc8818fe91288fdffcd</anchor>
      <arglist>(Int_t parnum, Int_t datnum, Double_t **x)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>variance</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>ab65424600034618edf519814ca108eaa</anchor>
      <arglist>(Int_t n, Double_t xm, Double_t *x)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>variance</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a16b603565b32814503cbab2cd3686ef7</anchor>
      <arglist>(Int_t n, Double_t *x)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>NBody.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/NBody/</path>
    <filename>NBody_8h</filename>
    <includes id="NBodyMath_8h" name="NBodyMath.h" local="no" imported="no">NBodyMath.h</includes>
    <includes id="Particle_8h" name="Particle.h" local="no" imported="no">Particle.h</includes>
    <includes id="System_8h" name="System.h" local="no" imported="no">System.h</includes>
    <namespace>NBody</namespace>
  </compound>
  <compound kind="file">
    <name>Particle.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/NBody/</path>
    <filename>Particle_8cxx</filename>
    <includes id="Particle_8h" name="Particle.h" local="no" imported="no">Particle.h</includes>
    <namespace>NBody</namespace>
    <member kind="function">
      <type>int</type>
      <name>DenCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a8573bb8cb7c6e930be22cc67752c56d3</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a6b5b0e13bcc9826748a451ecafabd1ba</anchor>
      <arglist>(ostream &amp;outs, const Particle &amp;p)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PIDCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a92da7f4fdf5e1a352cf7cb0ca181df0d</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>IDCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a26671a95adb3af98801f7623670951a2</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>RadCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a47ac977b109b91226271f35bc274bc7d</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>TypeCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ab9c42487f7eb685d00f481080ba6ffc7</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PotCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a72c2f2c6ad30d74216ad07683cda33b8</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PIDCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a92da7f4fdf5e1a352cf7cb0ca181df0d</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>IDCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a26671a95adb3af98801f7623670951a2</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>RadCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a47ac977b109b91226271f35bc274bc7d</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>TypeCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ab9c42487f7eb685d00f481080ba6ffc7</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PotCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a72c2f2c6ad30d74216ad07683cda33b8</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>Particle.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/NBody/</path>
    <filename>Particle_8h</filename>
    <includes id="NBodyMath_8h" name="NBodyMath.h" local="no" imported="no">NBodyMath.h</includes>
    <class kind="class">NBody::GasParticle</class>
    <class kind="class">NBody::Particle</class>
    <class kind="class">NBody::StarParticle</class>
    <namespace>NBody</namespace>
    <member kind="function">
      <type>int</type>
      <name>PIDCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a92da7f4fdf5e1a352cf7cb0ca181df0d</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>IDCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a26671a95adb3af98801f7623670951a2</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>RadCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a47ac977b109b91226271f35bc274bc7d</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>TypeCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ab9c42487f7eb685d00f481080ba6ffc7</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PotCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a72c2f2c6ad30d74216ad07683cda33b8</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="define">
      <type>#define</type>
      <name>MASSVAL</name>
      <anchorfile>Particle_8h.html</anchorfile>
      <anchor>a62fc6b0c1e122d201a60a8fa8fb2f901</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>Double_t</type>
      <name>DoublePos_t</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a234ec5660f1ec8cad4c118bf35f0ceca</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PIDCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a92da7f4fdf5e1a352cf7cb0ca181df0d</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>IDCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a26671a95adb3af98801f7623670951a2</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>RadCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a47ac977b109b91226271f35bc274bc7d</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>TypeCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ab9c42487f7eb685d00f481080ba6ffc7</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PotCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a72c2f2c6ad30d74216ad07683cda33b8</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>System.cxx</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/NBody/</path>
    <filename>System_8cxx</filename>
    <includes id="System_8h" name="System.h" local="no" imported="no">System.h</includes>
    <namespace>NBody</namespace>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ac19a01089d346407d2d5e58bf1a4b749</anchor>
      <arglist>(ostream &amp;outs, const System &amp;S)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>System.h</name>
    <path>/home/pelahi/myresearch/streamfinder/code/repo/VELOCIraptor-STF/stf/NBodylib/src/NBody/</path>
    <filename>System_8h</filename>
    <includes id="NBodyMath_8h" name="NBodyMath.h" local="no" imported="no">NBodyMath.h</includes>
    <includes id="Particle_8h" name="Particle.h" local="no" imported="no">Particle.h</includes>
    <class kind="class">NBody::System</class>
    <namespace>NBody</namespace>
  </compound>
  <compound kind="struct">
    <name>Math::Coord</name>
    <filename>structMath_1_1Coord.html</filename>
    <member kind="function">
      <type>Coord &amp;</type>
      <name>operator=</name>
      <anchorfile>structMath_1_1Coord.html</anchorfile>
      <anchor>a466328db2c7171cbcfa14b4dc4e56c6b</anchor>
      <arglist>(const Coord &amp;q)</arglist>
    </member>
    <member kind="function">
      <type>Coord &amp;</type>
      <name>operator=</name>
      <anchorfile>structMath_1_1Coord.html</anchorfile>
      <anchor>ae0d2b1e1301e8a5fe64f1d95c4c663b9</anchor>
      <arglist>(const Double_t q[3])</arglist>
    </member>
    <member kind="variable">
      <type>Double_t</type>
      <name>pos</name>
      <anchorfile>structMath_1_1Coord.html</anchorfile>
      <anchor>a13d405774fcaf61667f2ee7f584f448a</anchor>
      <arglist>[3]</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Math::Coordinate</name>
    <filename>classMath_1_1Coordinate.html</filename>
    <member kind="function">
      <type></type>
      <name>Coordinate</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a704807d0b749ba5b04dcb0e1e5591898</anchor>
      <arglist>(const Double_t x=0, const Double_t y=0, const Double_t z=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Coordinate</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a381c8935e37ad8b94e100575fa3dc99f</anchor>
      <arglist>(const Double_t *array)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Coordinate</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a787637d3debc7561fe3eee651042b322</anchor>
      <arglist>(const Real_t *array)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Coordinate</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a3287d790b4b261541ae774037df3414f</anchor>
      <arglist>(const Coord q)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Coordinate</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a4ed04a87c224621fea74c2a46dd338f7</anchor>
      <arglist>(const Coordinate &amp;q)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>X</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a3a3c6403a9f72d94d3298a89c76fc1e5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Y</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a2355804d9f9e46c56aa14f9960e17710</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Z</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>ab7960403767c9857969e7054f76258e8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t *</type>
      <name>GetCoord</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>ac104f94da86e1bbbed16d2f437d8fcae</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Length</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>ae18d0e7bb69443bbe4f261adc9697202</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>Normal</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a93ae92b0fc72055f5a611235d3ebefa5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>Cross</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>adf30dda6a9f84da641d00b840ff97898</anchor>
      <arglist>(Coordinate c2)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>coord</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a1a4958363a8a5041351bb067c45b8473</anchor>
      <arglist>[3]</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Coordinate</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a704807d0b749ba5b04dcb0e1e5591898</anchor>
      <arglist>(const Double_t x=0, const Double_t y=0, const Double_t z=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Coordinate</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a381c8935e37ad8b94e100575fa3dc99f</anchor>
      <arglist>(const Double_t *array)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Coordinate</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a787637d3debc7561fe3eee651042b322</anchor>
      <arglist>(const Real_t *array)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Coordinate</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a3287d790b4b261541ae774037df3414f</anchor>
      <arglist>(const Coord q)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Coordinate</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a4ed04a87c224621fea74c2a46dd338f7</anchor>
      <arglist>(const Coordinate &amp;q)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>X</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a3a3c6403a9f72d94d3298a89c76fc1e5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Y</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a2355804d9f9e46c56aa14f9960e17710</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Z</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>ab7960403767c9857969e7054f76258e8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t *</type>
      <name>GetCoord</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>ac104f94da86e1bbbed16d2f437d8fcae</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Coordinate &amp;</type>
      <name>operator=</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a0c13c3cae4b37e6670c7ccd5fa493249</anchor>
      <arglist>(const Double_t *array)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate &amp;</type>
      <name>operator=</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a351eaf6d341db370cf00050e8891cd9b</anchor>
      <arglist>(const Real_t *array)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate &amp;</type>
      <name>operator=</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a92eb62a14dce461e516cd299bce7dcb3</anchor>
      <arglist>(const Coord &amp;q)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate &amp;</type>
      <name>operator=</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a3cac4ce2563547145bac17632018bdc1</anchor>
      <arglist>(const Coordinate q)</arglist>
    </member>
    <member kind="function">
      <type>Double_t &amp;</type>
      <name>operator[]</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a03bc2739d319000cde8d002e13651f9f</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>const Double_t &amp;</type>
      <name>operator[]</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>ab06abbd67dd2fefc6108be95902f14f7</anchor>
      <arglist>(int i) const </arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>operator*</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>acc1b2448915caaa34c24dc75eaab7243</anchor>
      <arglist>(Double_t a)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>operator/</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>aad3a0587f02ab906ad08c13945b70689</anchor>
      <arglist>(Double_t a)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>operator*</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a91a21fbd4fea969674b9e8b16bf38ded</anchor>
      <arglist>(const Coordinate &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>operator+</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a8077ecee1186dabd5c0d4577b3d8313f</anchor>
      <arglist>(const Coordinate &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>operator-</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a312c712414bf22268a073ff9869c6592</anchor>
      <arglist>(const Coordinate &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator+=</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a75d9a414794114275d4c7cd477f22cf6</anchor>
      <arglist>(const Coordinate &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator*=</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a7e8692223f0f5042d92d479323a288f4</anchor>
      <arglist>(const Coordinate &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>operator*=</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a49faf9f56766277bcb4703995ec3cd81</anchor>
      <arglist>(const Double_t &amp;a)</arglist>
    </member>
    <member kind="friend">
      <type>friend Coordinate</type>
      <name>operator*</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a1d8b05251f0b93185de7fac0b8033666</anchor>
      <arglist>(Double_t a, const Coordinate &amp;c)</arglist>
    </member>
    <member kind="friend">
      <type>friend std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a0b7db74177811858ed94ee221f1eec8d</anchor>
      <arglist>(std::ostream &amp;outs, Coordinate c)</arglist>
    </member>
    <member kind="friend">
      <type>friend std::istream &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a20b6c043969a51b146cc89b3567dd81c</anchor>
      <arglist>(std::istream &amp;ins, Coordinate &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Length</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>ae18d0e7bb69443bbe4f261adc9697202</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>Normal</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>a93ae92b0fc72055f5a611235d3ebefa5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>Cross</name>
      <anchorfile>classMath_1_1Coordinate.html</anchorfile>
      <anchor>adf30dda6a9f84da641d00b840ff97898</anchor>
      <arglist>(Coordinate c2)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Math::Coordinate2D</name>
    <filename>classMath_1_1Coordinate2D.html</filename>
    <member kind="function">
      <type></type>
      <name>Coordinate2D</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>abdda75c08bd00bc03f7ac304dd1c2fcc</anchor>
      <arglist>(const Double_t x=0, const Double_t y=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Coordinate2D</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>a33b14467507adb902923aee4da1e512b</anchor>
      <arglist>(const Double_t *array)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Coordinate2D</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>a37841e7a2a869a7b34cb5c2b57ed0440</anchor>
      <arglist>(const Real_t *array)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>X</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>ad60e75117d4abcea8b4ed6e52632b1b1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Y</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>a573c286e11048761785a433ccb1c1a44</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t *</type>
      <name>Coord</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>adc9e4724cadd8f6103453988851b57a1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Length</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>addc90dac48f581f0d480756e5c8de0ee</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Coordinate2D</type>
      <name>Normal</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>a9bc1391beea6a04d4ce05d198ea73313</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>coord</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>a519fba42d31f620cc105e144d2e6cbd2</anchor>
      <arglist>[2]</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Coordinate2D</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>abdda75c08bd00bc03f7ac304dd1c2fcc</anchor>
      <arglist>(const Double_t x=0, const Double_t y=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Coordinate2D</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>a33b14467507adb902923aee4da1e512b</anchor>
      <arglist>(const Double_t *array)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Coordinate2D</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>a37841e7a2a869a7b34cb5c2b57ed0440</anchor>
      <arglist>(const Real_t *array)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>X</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>ad60e75117d4abcea8b4ed6e52632b1b1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Y</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>a573c286e11048761785a433ccb1c1a44</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t *</type>
      <name>Coord</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>adc9e4724cadd8f6103453988851b57a1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Coordinate2D &amp;</type>
      <name>operator=</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>a7c0df8443c2b18eb0ddfd25d0b8b0b07</anchor>
      <arglist>(const Double_t *array)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate2D &amp;</type>
      <name>operator=</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>a596849468e5115353477fc27901544b6</anchor>
      <arglist>(const Real_t *array)</arglist>
    </member>
    <member kind="function">
      <type>Double_t &amp;</type>
      <name>operator[]</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>af7ff91e0e2017b0099879d42a95d41bc</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>const Double_t &amp;</type>
      <name>operator[]</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>ae95c220cc122504606b4410076999629</anchor>
      <arglist>(int i) const </arglist>
    </member>
    <member kind="function">
      <type>Coordinate2D</type>
      <name>operator*</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>ad72fd1cace8c0920516b9d3f7ea93a2c</anchor>
      <arglist>(Double_t a)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate2D</type>
      <name>operator/</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>a66e056a2f7cfb3bd4e3bfe6b6c5a3c17</anchor>
      <arglist>(Double_t a)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>operator*</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>af1dc3cb7e62d295bc852508526444e3f</anchor>
      <arglist>(const Coordinate2D &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate2D</type>
      <name>operator+</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>ad29bc7e5a726bd7c6de7219d293721c2</anchor>
      <arglist>(const Coordinate2D &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate2D</type>
      <name>operator-</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>a6f2779b53ea64ef98b2ca5898cb20eef</anchor>
      <arglist>(const Coordinate2D &amp;c)</arglist>
    </member>
    <member kind="friend">
      <type>friend Coordinate2D</type>
      <name>operator*</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>ac9fdba206fad81be46c98136679ec13a</anchor>
      <arglist>(Double_t a, const Coordinate2D &amp;c)</arglist>
    </member>
    <member kind="friend">
      <type>friend std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>a6aee2d3c935fc48107f6a8b4722b78ef</anchor>
      <arglist>(std::ostream &amp;outs, Coordinate2D c)</arglist>
    </member>
    <member kind="friend">
      <type>friend std::istream &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>adaee0fdaa279b836835e4a679a5df710</anchor>
      <arglist>(std::istream &amp;ins, Coordinate2D &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Length</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>addc90dac48f581f0d480756e5c8de0ee</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Coordinate2D</type>
      <name>Normal</name>
      <anchorfile>classMath_1_1Coordinate2D.html</anchorfile>
      <anchor>a9bc1391beea6a04d4ce05d198ea73313</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>NBody::GasParticle</name>
    <filename>classNBody_1_1GasParticle.html</filename>
    <base>NBody::Particle</base>
    <member kind="function">
      <type></type>
      <name>GasParticle</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a1a712f9df20a6d85e4f343ba17816bd0</anchor>
      <arglist>(Double_t Mass=0, Double_t x=0, Double_t y=0, Double_t z=0, Double_t vx=0, Double_t vy=0, Double_t vz=0, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Double_t Temp=0, Double_t Ui=0, Double_t Pi=0, Double_t NE=0, Double_t NH0=0, Double_t Zi=0, Double_t SFR=0, Double_t LGS=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GasParticle</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a7bd2256911964d2f21b04dfdd9cd01d2</anchor>
      <arglist>(Double_t Mass, Double_t *NewPos, Double_t *NewVel, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Double_t Temp=0, Double_t Ui=0, Double_t Pi=0, Double_t NE=0, Double_t NH0=0, Double_t Zi=0, Double_t SFR=0, Double_t LGS=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GasParticle</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a6657bde6cfa4d0cdf72cc14b7033b13b</anchor>
      <arglist>(const GasParticle &amp;p)</arglist>
    </member>
    <member kind="function">
      <type>GasParticle &amp;</type>
      <name>operator=</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>ab5bbd24304ad7a0b7604dde5e046456e</anchor>
      <arglist>(const GasParticle &amp;part)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator==</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a457dd75a5b4f4fad689cc947a7d8f5c0</anchor>
      <arglist>(const GasParticle &amp;p)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator!=</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>ad4a0026b8fc3df45c6d4127f49837bbb</anchor>
      <arglist>(const GasParticle &amp;p) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetTemp</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a4c2159dcf14222a318cae8cbbe23aa3e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetTemp</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a510e7325574374b5b6cf185faf3076a9</anchor>
      <arglist>(Double_t Temp)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetPressure</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>af968d83366dfeb31d502882d953b3427</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetPressure</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>aed6563b2276d086178e2a5d410c4b0c4</anchor>
      <arglist>(Double_t Pi)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetU</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a654b88f60ea4a75b0784f0ea67e4492d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetU</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a1ba0ab0cc1ddd4ef52c1bd6738b123eb</anchor>
      <arglist>(Double_t Ui)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetNe</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>ac1ecf73928bc666740fbc4960e4eaf1b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetNe</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>ac8cfa335e805a7d38aa87f143df632ed</anchor>
      <arglist>(Double_t NE)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetNh0</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a62a9954b2a6dcf44d2e8dfb58caacf91</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetNh0</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a8e2679d0bef1260b1356924421fb8373</anchor>
      <arglist>(Double_t NH0)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetZ</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a56e96a7fb0b29359d0556159dc6b609e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetZ</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>affbe493486c3b220ee0ec4fc1a68bb78</anchor>
      <arglist>(Double_t Zi)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Getsfr</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a5f95da60eb86b46ee9acfa1dd92100c8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Setsfr</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a09c0ce251d1bed9f80ee36702ccbf1c8</anchor>
      <arglist>(Double_t SFR)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetEntropy</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a5f9082a71e75801991b40a3fd08abf6f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetEntropy</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>ad0cedf346084220ca6c2621ad5aee566</anchor>
      <arglist>(Double_t lgs)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>temp</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a404980d3a2550165d198c3a76fde3bc0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>lgS</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a7ee70f1aa654d803d49415fd0832d89b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>P</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a58670787746667094e9f938897c8e1ae</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>Ne</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a0b4a6f5a3169f0a3b981675f04b8f31c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>Nh0</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a4e79e8b81c959026d6669bdbfbdb3719</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>metal</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>afa7e6dcd5a41555194715300c7934023</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>sfr</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>aca147af3d1bc6d6483e07e3bd6e193b7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>U</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a746892525a1272de9f7b44cc5a5748af</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>temp</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a404980d3a2550165d198c3a76fde3bc0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>lgS</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a7ee70f1aa654d803d49415fd0832d89b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>P</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a58670787746667094e9f938897c8e1ae</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>Ne</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a0b4a6f5a3169f0a3b981675f04b8f31c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>Nh0</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a4e79e8b81c959026d6669bdbfbdb3719</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>metal</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>afa7e6dcd5a41555194715300c7934023</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>sfr</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>aca147af3d1bc6d6483e07e3bd6e193b7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>U</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a746892525a1272de9f7b44cc5a5748af</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GasParticle</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a1a712f9df20a6d85e4f343ba17816bd0</anchor>
      <arglist>(Double_t Mass=0, Double_t x=0, Double_t y=0, Double_t z=0, Double_t vx=0, Double_t vy=0, Double_t vz=0, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Double_t Temp=0, Double_t Ui=0, Double_t Pi=0, Double_t NE=0, Double_t NH0=0, Double_t Zi=0, Double_t SFR=0, Double_t LGS=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GasParticle</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a7bd2256911964d2f21b04dfdd9cd01d2</anchor>
      <arglist>(Double_t Mass, Double_t *NewPos, Double_t *NewVel, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Double_t Temp=0, Double_t Ui=0, Double_t Pi=0, Double_t NE=0, Double_t NH0=0, Double_t Zi=0, Double_t SFR=0, Double_t LGS=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GasParticle</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a6657bde6cfa4d0cdf72cc14b7033b13b</anchor>
      <arglist>(const GasParticle &amp;p)</arglist>
    </member>
    <member kind="function">
      <type>GasParticle &amp;</type>
      <name>operator=</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>ab5bbd24304ad7a0b7604dde5e046456e</anchor>
      <arglist>(const GasParticle &amp;part)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator==</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a457dd75a5b4f4fad689cc947a7d8f5c0</anchor>
      <arglist>(const GasParticle &amp;p)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator!=</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>ad4a0026b8fc3df45c6d4127f49837bbb</anchor>
      <arglist>(const GasParticle &amp;p) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetTemp</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a4c2159dcf14222a318cae8cbbe23aa3e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetTemp</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a510e7325574374b5b6cf185faf3076a9</anchor>
      <arglist>(Double_t Temp)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetPressure</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>af968d83366dfeb31d502882d953b3427</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetPressure</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>aed6563b2276d086178e2a5d410c4b0c4</anchor>
      <arglist>(Double_t Pi)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetU</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a654b88f60ea4a75b0784f0ea67e4492d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetU</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a1ba0ab0cc1ddd4ef52c1bd6738b123eb</anchor>
      <arglist>(Double_t Ui)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetNe</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>ac1ecf73928bc666740fbc4960e4eaf1b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetNe</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>ac8cfa335e805a7d38aa87f143df632ed</anchor>
      <arglist>(Double_t NE)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetNh0</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a62a9954b2a6dcf44d2e8dfb58caacf91</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetNh0</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a8e2679d0bef1260b1356924421fb8373</anchor>
      <arglist>(Double_t NH0)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetZ</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a56e96a7fb0b29359d0556159dc6b609e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetZ</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>affbe493486c3b220ee0ec4fc1a68bb78</anchor>
      <arglist>(Double_t Zi)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Getsfr</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a5f95da60eb86b46ee9acfa1dd92100c8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Setsfr</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a09c0ce251d1bed9f80ee36702ccbf1c8</anchor>
      <arglist>(Double_t SFR)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetEntropy</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>a5f9082a71e75801991b40a3fd08abf6f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetEntropy</name>
      <anchorfile>classNBody_1_1GasParticle.html</anchorfile>
      <anchor>ad0cedf346084220ca6c2621ad5aee566</anchor>
      <arglist>(Double_t lgs)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Math::GMatrix</name>
    <filename>classMath_1_1GMatrix.html</filename>
    <member kind="function">
      <type></type>
      <name>GMatrix</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>adbad0708a02d06fb1443cef98993a4e4</anchor>
      <arglist>(int r, int c, Double_t *data=NULL)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GMatrix</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a06cdbdf589374d0e2614119e618db1c1</anchor>
      <arglist>(const Coordinate &amp;c)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GMatrix</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a3e8077f199f916e077b46d17e0fe0918</anchor>
      <arglist>(const Matrix &amp;m)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GMatrix</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a10b86fbedc99382688002a8eec556203</anchor>
      <arglist>(const GMatrix &amp;m)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GMatrix</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>aa29be6dee5eb02558a1e88c88f50a174</anchor>
      <arglist>(int N)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~GMatrix</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a6ff081194d0c57b6238874c4787b3363</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>Row</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>aed816931e895ce8d77fe8a5c8b549021</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>Col</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a71112b3d5ab5951ae58aad91ab4f6eca</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>Rank</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>ae0cbee7a43e363a6ceb8680750f3e69a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Trace</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a5561d923d0cb0824811c07be9493b2c2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Transpose</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>aabcc99abbacb785388f910963b9d247a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>TransposeInPlace</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a4750210a08c0db83979f1d2d99905148</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Diag</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>adaeb63d0b2ce41163c891ed74e86fbcf</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Pivot</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a284671ddfb10ce9d0775c30678c2c6c7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>PivotInPlace</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a8ca24741d48e99c43ffb178fc2a19b18</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>SubMatrix</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>ae673400d01ba60945851d5d9ee47f518</anchor>
      <arglist>(int ra, int rb, int ca, int cb) const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>ExtractMatrix</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a967ad7d90205797738d933f9d093564f</anchor>
      <arglist>(int r, int c) const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>isSymmetric</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a41243e775e6d3b22949d3df12165d705</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>isZero</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>aa4944521be0ba0626e8afe3d7e751284</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Det</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a58da56a89b05edd3525c269dc2b69428</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>LUDecomp</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a851c93aaec2724cbd30625b9bb8ec1c3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Lower_Triangular_Inverse</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>ac7d38c92b777013dc8a1b59374a11882</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Unit_Upper_Triangular_Inverse</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a55e13eb6ccb1c0e281cfb6347102a128</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Adjugate</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a415600fca3ed1802c9a989e89bd3cc89</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Inverse</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a65a3101c47343fb53d474125ea892fd3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>InversewithPivot</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a63207bfe7029fa438d498dd7eb7cf632</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ROTATE</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>ab9407d80d7d07a83adf01f99f4370cb9</anchor>
      <arglist>(GMatrix &amp;m, int i, int j, int k, int l, Double_t sinrot, Double_t tanrot) const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Jacobi</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a161b2f98900c7e98a0d69131633c91b9</anchor>
      <arglist>(Double_t tol=1e-2) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Jacobi</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a3c41ae8d33ae865680859412992ec082</anchor>
      <arglist>(GMatrix &amp;eigenval, GMatrix &amp;Prot, Double_t tol=1e-2) const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Eigenvalues</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a0fe3e119840c2bb802e11f4d9f06dbe3</anchor>
      <arglist>(Double_t tol=1e-2) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Eigenvalvec</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a474b6a4fb50a405430691e0c85893dcc</anchor>
      <arglist>(GMatrix &amp;eigenval, GMatrix &amp;eigenvector, Double_t tol=1e-2) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>int</type>
      <name>compare_Double</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a6d6263146aa3adbf04d4e3f4db58cc90</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>col</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a65be03caa33cd483e47a6291f406f39e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t *</type>
      <name>matrix</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>aa53f90cdd0dbc75d65b96fcb31698ee9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>row</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a5f2efa17859a5972be93237b7eab403d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GMatrix</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>adbad0708a02d06fb1443cef98993a4e4</anchor>
      <arglist>(int r, int c, Double_t *data=NULL)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GMatrix</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a06cdbdf589374d0e2614119e618db1c1</anchor>
      <arglist>(const Coordinate &amp;c)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GMatrix</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a3e8077f199f916e077b46d17e0fe0918</anchor>
      <arglist>(const Matrix &amp;m)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GMatrix</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a10b86fbedc99382688002a8eec556203</anchor>
      <arglist>(const GMatrix &amp;m)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>GMatrix</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>aa29be6dee5eb02558a1e88c88f50a174</anchor>
      <arglist>(int N)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~GMatrix</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a6ff081194d0c57b6238874c4787b3363</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>Row</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>aed816931e895ce8d77fe8a5c8b549021</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>Col</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a71112b3d5ab5951ae58aad91ab4f6eca</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>Rank</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>ae0cbee7a43e363a6ceb8680750f3e69a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t &amp;</type>
      <name>operator()</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a5a3ca134e485e53d8c2e9de65fc505ac</anchor>
      <arglist>(int i, int j)</arglist>
    </member>
    <member kind="function">
      <type>const Double_t &amp;</type>
      <name>operator()</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>aceab2f1c4cf77aef229ed9e375be429a</anchor>
      <arglist>(int i, int j) const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>operator+</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>aadfe2e4423d77e108da84540ad6852fd</anchor>
      <arglist>(const GMatrix &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>operator-</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a3ade6d14ae6c85513bcdcd6b26807bc5</anchor>
      <arglist>(const GMatrix &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>operator*</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a5a194e3f119f7dc715d0e66fd6469272</anchor>
      <arglist>(Double_t a)</arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>operator*</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a60141e82b341325f00342c3dd7fd3605</anchor>
      <arglist>(const GMatrix &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>GMatrix &amp;</type>
      <name>operator=</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a73175180c6dd6a60b7ee26cd79d2a873</anchor>
      <arglist>(const GMatrix m)</arglist>
    </member>
    <member kind="friend">
      <type>friend GMatrix</type>
      <name>operator*</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a07ef05e12bfadf82c8dca0751b9b2b5b</anchor>
      <arglist>(Double_t a, const GMatrix &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Trace</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a5561d923d0cb0824811c07be9493b2c2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Transpose</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>aabcc99abbacb785388f910963b9d247a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>TransposeInPlace</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a4750210a08c0db83979f1d2d99905148</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Diag</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>adaeb63d0b2ce41163c891ed74e86fbcf</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Pivot</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a284671ddfb10ce9d0775c30678c2c6c7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>PivotInPlace</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a8ca24741d48e99c43ffb178fc2a19b18</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>SubMatrix</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>ae673400d01ba60945851d5d9ee47f518</anchor>
      <arglist>(int ra, int rb, int ca, int cb) const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>ExtractMatrix</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a967ad7d90205797738d933f9d093564f</anchor>
      <arglist>(int r, int c) const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>isSymmetric</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a41243e775e6d3b22949d3df12165d705</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>isZero</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>aa4944521be0ba0626e8afe3d7e751284</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Det</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a58da56a89b05edd3525c269dc2b69428</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>LUDecomp</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a851c93aaec2724cbd30625b9bb8ec1c3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Lower_Triangular_Inverse</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>ac7d38c92b777013dc8a1b59374a11882</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Unit_Upper_Triangular_Inverse</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a55e13eb6ccb1c0e281cfb6347102a128</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Adjugate</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a415600fca3ed1802c9a989e89bd3cc89</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Inverse</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a65a3101c47343fb53d474125ea892fd3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>InversewithPivot</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a63207bfe7029fa438d498dd7eb7cf632</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ROTATE</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>ab9407d80d7d07a83adf01f99f4370cb9</anchor>
      <arglist>(GMatrix &amp;m, int i, int j, int k, int l, Double_t sinrot, Double_t tanrot) const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Jacobi</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a161b2f98900c7e98a0d69131633c91b9</anchor>
      <arglist>(Double_t tol=1e-2) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Jacobi</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a3c41ae8d33ae865680859412992ec082</anchor>
      <arglist>(GMatrix &amp;eigenval, GMatrix &amp;Prot, Double_t tol=1e-2) const </arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>Eigenvalues</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a0fe3e119840c2bb802e11f4d9f06dbe3</anchor>
      <arglist>(Double_t tol=1e-2) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Eigenvalvec</name>
      <anchorfile>classMath_1_1GMatrix.html</anchorfile>
      <anchor>a474b6a4fb50a405430691e0c85893dcc</anchor>
      <arglist>(GMatrix &amp;eigenval, GMatrix &amp;eigenvector, Double_t tol=1e-2) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>NBody::KDTree</name>
    <filename>classNBody_1_1KDTree.html</filename>
    <member kind="function">
      <type></type>
      <name>KDTree</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad68f262da92174fcc384a9d5a8774a40</anchor>
      <arglist>(Particle *p, Int_t numparts, Int_t bucket_size=16, int TreeType=TPHYS, int KernType=KEPAN, int KernRes=1000, int SplittingCriterion=0, int Aniso=0, int ScaleSpace=0, Double_t *Period=NULL, Double_t **metric=NULL)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>KDTree</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad34cea8c3ba5f18de5c424072f969b29</anchor>
      <arglist>(System &amp;s, Int_t bucket_size=16, int TreeType=TPHYS, int KernType=KEPAN, int KernRes=1000, int SplittingCriterion=0, int Aniso=0, int ScaleSpace=0, Double_t **metric=NULL)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~KDTree</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aa2abf2018780d1be3d52f44c4899e52c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetNumNodes</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a95891d021490d3f805e6a9e97ecb6c6d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetNumLeafNodes</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a32bff63406186f876844783e6fd07835</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetBucketSize</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a94a2c83edc5de03b3804e1bd5d82818f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetTreeType</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ab23fa720d57157f80d861c3e238ed27e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetKernType</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a6c7ac9d27697119a0526c2cd921fbf32</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetKernNorm</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a0bb7ab867cd4335ad859c5d995ccdf10</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Node *</type>
      <name>GetRoot</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad60c3768328d7e2153b5ac3cdaf0fc55</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetPeriod</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aa3b489bb80090cde772713f02d7bb07d</anchor>
      <arglist>(int j)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearest</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ae3c79e5fb889d4a4d14ecc0ae8374f7f</anchor>
      <arglist>(Int_t **nn, Double_t **dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>add9c46bd3e5a4bfa8b0402c876e903dd</anchor>
      <arglist>(Int_t **nn, Double_t **dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ab3cf73fd48ae33e6d9530e8c054386a5</anchor>
      <arglist>(Int_t **nn, Double_t **dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>afeb98226c436f301b68dcfe5a9311c4e</anchor>
      <arglist>(Int_t **nn, Double_t **dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a845d0a6df00755a297ea34eb3b07f747</anchor>
      <arglist>(FOFcompfunc cmp, Double_t *params, Int_t **nn, Double_t **dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearest</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a1d249268e725c375ae50eb2a6dcec7c0</anchor>
      <arglist>(Int_t tt, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a52a3be0cfc421225695c34f330f6172c</anchor>
      <arglist>(Int_t tt, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a23098f6033aff6729bc92c04edf2a8f3</anchor>
      <arglist>(Int_t tt, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a35ed1f74def57d237ff518244144090a</anchor>
      <arglist>(Int_t tt, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearest</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a5b3ef1437766578990c3074cd1536f85</anchor>
      <arglist>(Double_t *x, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad457bf15291b9a388113f13303e98a2b</anchor>
      <arglist>(Double_t *x, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ae7a94c4106653116bf53ca7a7619ea25</anchor>
      <arglist>(Double_t *v, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a633aa4e77cf2ff05eb47364bf2cf430b</anchor>
      <arglist>(Double_t *x, Double_t *v, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a410b409750426b307aef5b24ecc88958</anchor>
      <arglist>(Coordinate x, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ae8ac7e4c970c59c1ba9d2d3804de0026</anchor>
      <arglist>(Coordinate v, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a73b1b578760cb174e9617fb5efb1362a</anchor>
      <arglist>(Coordinate x, Coordinate v, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a7e0a8bd0362b3c478c72034755cf6aea</anchor>
      <arglist>(Int_t tt, FOFcompfunc cmp, Double_t *params, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>af98daa0ca7fcc6fb12c409b6761ea20a</anchor>
      <arglist>(Particle p, FOFcompfunc cmp, Double_t *params, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBall</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a86888b7be77da9ad42691cd87db929a5</anchor>
      <arglist>(Int_t tt, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a57f312f21caa5aa1316561d1db1ec217</anchor>
      <arglist>(Int_t tt, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a57b980a39e4c0741abc70abf74635a01</anchor>
      <arglist>(Int_t tt, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPhase</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ac69e103e8b8da8bc02591caa7eccb135</anchor>
      <arglist>(Int_t tt, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBall</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a3a9113a3650d488f4176b3f942db1d7c</anchor>
      <arglist>(Double_t *x, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a452e2431da28b8753bcf7af0e41a2bae</anchor>
      <arglist>(Double_t *x, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a988f3e538c2337d4a49254824cc7dd25</anchor>
      <arglist>(Double_t *v, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPhase</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad0d4bb038fd1cde03c2bee06867fa9c9</anchor>
      <arglist>(Double_t *x, Double_t *v, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBall</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a126899aaba1fa6a5ddac9d20ed2b4e2d</anchor>
      <arglist>(Coordinate x, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a372afa16b3f280268d7025e8e6f360b0</anchor>
      <arglist>(Coordinate x, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ae0b83f4d4bcb3f251b93d68c3ece3ff8</anchor>
      <arglist>(Coordinate v, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPhase</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a3bb62c923f6a8c2c265d012cc53b6a5b</anchor>
      <arglist>(Coordinate x, Coordinate v, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterion</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a3d4560dd709924c2b326fa375271c8fa</anchor>
      <arglist>(Int_t tt, FOFcompfunc p, Double_t *params, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterion</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a63d820687992c468f0890a77ea2da3b1</anchor>
      <arglist>(Int_t tt, FOFcompfunc p, Double_t *params, Int_t imark, Int_t *nn)</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a0a70772eb802c0fc39566d8da80a23af</anchor>
      <arglist>(Int_t tt, Double_t fdist2, Int_t *tagged)</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ab9cde2cb2f1b101a8073891a7a0c9a9b</anchor>
      <arglist>(Coordinate x, Double_t fdist2, Int_t *tagged)</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a4480b407eb873fdb0ae977b911e27251</anchor>
      <arglist>(Double_t *x, Double_t fdist2, Int_t *tagged)</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>SearchCriterionTagged</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a870da6c7c1bae8e5babc3a507d2eb4a8</anchor>
      <arglist>(Int_t tt, FOFcompfunc cmp, Double_t *params, Int_t *tagged)</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>SearchCriterionTagged</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a91d66c4b62c2ba71e80791c310f9b795</anchor>
      <arglist>(Particle &amp;p, FOFcompfunc cmp, Double_t *params, Int_t *tagged)</arglist>
    </member>
    <member kind="function">
      <type>Node *</type>
      <name>FindLeafNode</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a0cc5ead98ebb862fdf506ccb9d00277a</anchor>
      <arglist>(Int_t tt)</arglist>
    </member>
    <member kind="function">
      <type>Node *</type>
      <name>FindLeafNode</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a5586bfa61c5d45993a95cb475558e547</anchor>
      <arglist>(Double_t *x)</arglist>
    </member>
    <member kind="function">
      <type>Node *</type>
      <name>FindLeafNode</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a2a2a13a32a18b93a2db01b93d24d3fe4</anchor>
      <arglist>(Double_t search[6][2])</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>CalcDensity</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad4e03bf55723a599850fe3e15af5add2</anchor>
      <arglist>(Int_t Nsmooth=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>CalcVelDensity</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a3e42b28de44aad77516d091fe72c053b</anchor>
      <arglist>(Int_t Nsmooth=64, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>Double_t *</type>
      <name>CalcVelDensityWithPhysDensity</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a2e7ea3a186de7baf943bead4e1968a67</anchor>
      <arglist>(Int_t Nsmooth=64, Int_t Nsearch=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate *</type>
      <name>CalcSmoothVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>af25297a80fe72a8396c74266e19c6d5f</anchor>
      <arglist>(Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Matrix *</type>
      <name>CalcSmoothVelDisp</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a97d4ba62011bf08948ac272217d1313d</anchor>
      <arglist>(Coordinate *smvel, Int_t Nsmooth=64, int densityset=1, int meanvelset=1)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate *</type>
      <name>CalcSmoothVelSkew</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>af0cd9c8fdadc3ac92b941aea719a8b9d</anchor>
      <arglist>(Coordinate *smvel, Matrix *smveldisp, Int_t Nsmooth=64, int densityset=1, int meanvelset=1, int veldispset=1)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate *</type>
      <name>CalcSmoothVelKurtosis</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ae058f09de1d4a983288c5c5f02375c1e</anchor>
      <arglist>(Coordinate *smvel, Matrix *smveldisp, Int_t Nsmooth=64, int densityset=1, int meanvelset=1, int veldispset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcDensityParticle</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ae07ac97de0b123ecfe2e2e81fc707524</anchor>
      <arglist>(Int_t target, Int_t Nsmooth=64)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcVelDensityParticle</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a44a92a218c9601df2f9cd61e27771000</anchor>
      <arglist>(Int_t target, Int_t Nsmooth=64, Int_t Nsearch=64, int iflag=0, PriorityQueue *pq=NULL, PriorityQueue *pq2=NULL, Int_t *nnIDs=NULL, Double_t *vdist=NULL)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcVelDensityWithPhysDensityParticle</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a52a014e7b489f45461b64218abc2c205</anchor>
      <arglist>(Int_t target, Int_t Nsmooth=64, Int_t Nsearch=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>CalcSmoothVelParticle</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aa2a21e48da4f5283c109d2c852fd91d2</anchor>
      <arglist>(Int_t target, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>CalcSmoothVelDispParticle</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a226fe168e0c666e33b725e97346dd2f5</anchor>
      <arglist>(Int_t target, Coordinate smvel, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcDensityPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a7beae71edc364f8fc4d9d4d3748340cc</anchor>
      <arglist>(Double_t *x, Int_t Nsmooth=64, Double_t *v=NULL)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcVelDensityPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a9db41f8594c17279ae6f3313e7006ec2</anchor>
      <arglist>(Double_t *x, Double_t *v, Int_t Nsmooth=64, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>CalcSmoothVelPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a7380400e1847af4a03db178219bec168</anchor>
      <arglist>(Double_t *x, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>CalcSmoothVelDispPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a72e6fcdc1e7febb1580980e658be7cb4</anchor>
      <arglist>(Double_t *x, Coordinate smvel, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcDensityPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad522decde9a620c9611a6784487de3ea</anchor>
      <arglist>(Coordinate x, Int_t Nsmooth=64, Coordinate v=Coordinate(0.))</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcVelDensityPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aa65db49c3b6fbb0cefc50c22764c9c7b</anchor>
      <arglist>(Coordinate x, Int_t Nsmooth=64, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>CalcSmoothVelPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aac3058c4b52c702c4646b7acb4994048</anchor>
      <arglist>(Coordinate x, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>CalcSmoothVelDispPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a6545172972f2e9e0d3e000b19faba207</anchor>
      <arglist>(Coordinate x, Coordinate smvel, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcDensityPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aeae8308d507973f11f7e4f486a91b8a8</anchor>
      <arglist>(Real_t *x, Int_t Nsmooth=64, Real_t *v=NULL)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcVelDensityPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>afa3ece24a460a6bfe6476aa5879825b7</anchor>
      <arglist>(Real_t *x, Real_t *v, Int_t Nsmooth=64, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>CalcSmoothVelPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a559378eedd2b016e287ab4c6966ce9a1</anchor>
      <arglist>(Real_t *x, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>CalcSmoothVelDispPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad382ad3084e42eecf474adc5f1e300ab</anchor>
      <arglist>(Real_t *x, Coordinate smvel, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t *</type>
      <name>CalcSmoothLocalMean</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a7f6c6ac46b636d85c7b5cc0bb65297fc</anchor>
      <arglist>(Double_t *weight, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t *</type>
      <name>CalcSmoothLocalDisp</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a6dd9b6a26aa7c7208fc2f7a82d510bec</anchor>
      <arglist>(Double_t *weight, Double_t *localmean, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcSmoothLocalMeanParticle</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a9631ff5badc7d042c9dbaee09013f8b8</anchor>
      <arglist>(Int_t target, Double_t *weight, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcSmoothLocalDispParticle</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>adef2f3a2bb53461f94bbecb16512c5ff</anchor>
      <arglist>(Int_t target, Double_t *weight, Double_t localmean, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcSmoothLocalMeanPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a669096967ed3e8386c1b9b886eff3180</anchor>
      <arglist>(Double_t *x, Double_t *weight, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcSmoothLocalDispPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a6b007b5d7df17540b24b3b69e6eb6b31</anchor>
      <arglist>(Double_t *x, Double_t *weight, Double_t localmean, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcSmoothLocalMeanPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a516c1ea1abd8060c8ab0fc39cce9e79e</anchor>
      <arglist>(Coordinate x, Double_t *weight, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcSmoothLocalDispPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>abb91b843b890016b4223fab569314143</anchor>
      <arglist>(Coordinate x, Double_t *weight, Double_t localmean, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcSmoothLocalValue</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a512192c8bd3d50bfd76056860b3015e9</anchor>
      <arglist>(Int_t Nsmooth, PriorityQueue *pq, Double_t *weight)</arglist>
    </member>
    <member kind="function">
      <type>Int_t *</type>
      <name>FOF</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a178903f3a706d9d36085d044396a35c7</anchor>
      <arglist>(Double_t fdist, Int_t &amp;numgroup, Int_t minnum=8, int order=0, Int_tree_t *pHead=NULL, Int_tree_t *pNext=NULL, Int_tree_t *pTail=NULL, Int_tree_t *pLen=NULL)</arglist>
    </member>
    <member kind="function">
      <type>Int_t *</type>
      <name>FOFCriterion</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a54433b1636f0e5c75bfa9b27fa1bc66c</anchor>
      <arglist>(FOFcompfunc p, Double_t *params, Int_t &amp;numgroups, Int_t minnum=8, int order=0, int ipcheckflag=0, FOFcheckfunc check=Pnocheck, Int_tree_t *pHead=NULL, Int_tree_t *pNext=NULL, Int_tree_t *pTail=NULL, Int_tree_t *pLen=NULL)</arglist>
    </member>
    <member kind="function">
      <type>Int_t *</type>
      <name>FOFCriterionSetBasisForLinks</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a68af3ba6a95b7c3b6a2d9d743d4ebc00</anchor>
      <arglist>(FOFcompfunc cmp, Double_t *params, Int_t &amp;numgroup, Int_t minnum=8, int order=0, int ipcheckflag=0, FOFcheckfunc check=Pnocheck, Int_tree_t *pHead=NULL, Int_tree_t *pNext=NULL, Int_tree_t *pTail=NULL, Int_tree_t *pLen=NULL)</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>FOFCriterionParticle</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a80277ea2f0337cba28616243685fb279</anchor>
      <arglist>(FOFcompfunc p, Int_t *pfof, Int_t target, Int_t iGroup, Double_t *params, Int_tree_t *pGroupHead, Int_tree_t *Fifo, Int_tree_t *pHead, Int_tree_t *pTail, Int_tree_t *pNext, Int_tree_t *pLen)</arglist>
    </member>
    <member kind="function">
      <type>Int_t *</type>
      <name>FOFNNCriterion</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a2f1a25ef3f1f8eca6c51c660703e26da</anchor>
      <arglist>(FOFcompfunc p, Double_t *params, Int_t numNN, Int_t **nnIDs, Int_t &amp;numgroups, Int_t minnum=8)</arglist>
    </member>
    <member kind="function">
      <type>Int_t *</type>
      <name>FOFNNDistCriterion</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a8172fbe3a807e627e3410027a90b0d99</anchor>
      <arglist>(FOFcompfunc p, Double_t *params, Int_t numNN, Int_t **nnIDs, Double_t **dist2, Double_t disfunc(Int_t, Double_t *), Int_t npc, Int_t *npca, Int_t &amp;numgroups, Int_t minnum=8)</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>TPHYS</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a3a37ad5dc668deb4dcb30a73ddb4ffaf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>TPROJ</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a74df3e8aa3f76cc6612f17e7dcf3a187</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>TVEL</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a2514a82816fafb1a61112647ecb8edd5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>TPHS</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a2bd32605634282709d8f3a6e0474b9bc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>TMETRIC</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a053b0fd24c3261b441c4776ac84a33c2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>KSPH</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a583a85ca3b644b4243d0243ff43f5b10</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>KGAUSS</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a6a2e3bc67ca903e47589cc74cf061207</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>KEPAN</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a14680ea01d84c82a930ae484ce983fe3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>KTH</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aa73ae200687c8a638e5d241bb5524d9b</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>Node *</type>
      <name>BuildNodes</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a8247f166032c6d40a1034b4adcdde96d</anchor>
      <arglist>(Int_t start, Int_t end)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>ScaleSpace</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>afb34afbe30171832d3316c20315671ed</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>int</type>
      <name>TreeTypeCheck</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a75d9ae25e0169f9878c1a10bb5151057</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>KernelConstruction</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a885b4b998c0d3f70e000ca34f09b9169</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>SpreadestPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ae05333d352d876ce064a9367cc75e765</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t *bnd)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>SpreadestVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a7d1f24e2981cc26899098296bd140111</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t *bnd)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>SpreadestPhs</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad3ab962356ec6265931e8b0e13118692</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t *bnd)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>BoundaryandMeanPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a5cf924f1a6bd70d8879564ff59670377</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t *bnd)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>BoundaryandMeanVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aa64d34d16856c8cc13344ef009378e56</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t *bnd)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>BoundaryandMeanPhs</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a28fec631a862b83783aeff72af8b685f</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t *bnd)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>DispersionPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aa6e08bf249a5aa09acf928c37035f8d0</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t mean)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>DispersionVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a132182b036ea58425562684ccc217a2b</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t mean)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>DispersionPhs</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a5978085e40e951636edf0a708e46e62c</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t mean)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>EntropyPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a1dbab902a9cc099f8f09e64c9fc583f4</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t low, Double_t up, Double_t nbins, Double_t *ni)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>EntropyVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a88c01aa506632afbd9dc059c8ad08906</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t low, Double_t up, Double_t nbins, Double_t *ni)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>EntropyPhs</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a19e7466f4ccefd6f5e1cd62e80f6217d</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t low, Double_t up, Double_t nbins, Double_t *ni)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>MedianPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a75c54d1a61a27ce5fe052294c24e3487</anchor>
      <arglist>(int d, Int_t k, Int_t start, Int_t end, bool balanced=true)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>MedianVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a518800579c8aeff8a136d84da504ac4b</anchor>
      <arglist>(int d, Int_t k, Int_t start, Int_t end, bool balanced=true)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>MedianPhs</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad625f0d10fdb7c06aac10b99131c085d</anchor>
      <arglist>(int d, Int_t k, Int_t start, Int_t end, bool balanced=true)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>Wsm</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>af702721d847b96b922092986ae970ac8</anchor>
      <arglist>(Double_t r, Int_t i, Int_t size, Double_t delta, Double_t *x)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>CalculateMetricSpacing</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a3db4480b0e0ad6215cd33507c00064ab</anchor>
      <arglist>(Int_t target, int treetype, Double_t *smetric)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>CalculateMetricSpacing</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a804d5e2cf646615ec4cb673f1e757049</anchor>
      <arglist>(const Double_t *x, const Double_t *v, int treetype, Double_t *smetric)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>CalculateMetricSpacing</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>abe82debeb779bf53f4817f72ebd4a0a2</anchor>
      <arglist>(const Real_t *x, const Real_t *v, int treetype, Double_t *smetric)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>CalculateMetricTensor</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a4713a6decf855de60a6e50e326216e7f</anchor>
      <arglist>(Int_t target, int treetype, Double_t *smetric, Double_t *metric, GMatrix &amp;gmetric)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>CalculateMetricTensor</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a8df252b1017292b815a085e2628721a0</anchor>
      <arglist>(const Double_t *x, const Double_t *v, int treetype, Double_t *smetric, Double_t *metric, GMatrix &amp;gmetric)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>CalculateMetricTensor</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ace41c3032fb804dc2417d278cfb10713</anchor>
      <arglist>(const Real_t *x, const Real_t *v, int treetype, Double_t *smetric, Double_t *metric, GMatrix &amp;gmetric)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>LoadNN</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a033a72ca0920c26739451fabf8512eb0</anchor>
      <arglist>(const Int_t ns, PriorityQueue *pq, Int_t *nn, Double_t *dist)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>anisotropic</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a460d37a06f894b7f62bcc0689265562d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Int_t</type>
      <name>b</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a565157282e9f8b9111d76c24656e5837</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Particle *</type>
      <name>bucket</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a17532bc511393bf212c431dcf491dbf3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t *</type>
      <name>derKernel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a53d7ef719f246d1556cda1818c3cec5a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>ivol</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ab914d19757850daca7a06435e07e46ee</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>ixvar</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>acfeb1f0d1486a08d32066e644e8b6ae7</anchor>
      <arglist>[MAXND]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t *</type>
      <name>Kernel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a868140107c03b058a19fec4ef6bf4450</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>smoothfunc</type>
      <name>kernfunc</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a10f13cf300d95ecb1ae8a8b344ea56c9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>kernfunctype</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a90fbcb17d3edfdf2d2eed9b241857026</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>kernnorm</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>af00442a6b6f40551799908c4f2ea2c69</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>kernres</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad6f3307176568a40e4df5e3f78d9a9fb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t **</type>
      <name>metric</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>af0089d8b0a6125d938723a2af7f5f7c9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>ND</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a2efecb5969da59830b59b894c2251f41</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t *</type>
      <name>nientropy</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a153dfcbde713c31f8825edaf7938fe6e</anchor>
      <arglist>[MAXND]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Int_t</type>
      <name>numleafnodes</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ae8e00754520360b4a8aa9cd7122c8c6a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Int_t</type>
      <name>numnodes</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a20d2f9caf81d7d6a133eed551bd67aab</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Int_t</type>
      <name>numparts</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a16423d8cf7072d3d0fd22c807317e266</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t *</type>
      <name>period</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a509c6647ed8e15f4336d8398818ad52b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Node *</type>
      <name>root</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a1360a87e53536870e2d1c2364370ffe2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>scalespace</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a5a26649d58e19a855018f760a5995f9a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>splittingcriterion</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad1214f63d6eaa1598a3234cf781cc596</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>treetype</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a4d27393ca72f405d855ae11ad9141aa7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>vol</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a5fbb8d39093c6475259ba745efbd64ac</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>xmean</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>af46b2fe996c0f70de0aede109e659e06</anchor>
      <arglist>[MAXND]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>xvar</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a82b27d51d7e98cdac58f466e14f6b3e5</anchor>
      <arglist>[MAXND]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t(NBody::KDTree::*</type>
      <name>bmfunc</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a5e6fdd1be0fb12c2648ef9e844763e2a</anchor>
      <arglist>)(int, Int_t, Int_t, Double_t *)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t(NBody::KDTree::*</type>
      <name>dispfunc</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>af73cb9c00df90d290e7cdb32456ed59a</anchor>
      <arglist>)(int, Int_t, Int_t, Double_t)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t(NBody::KDTree::*</type>
      <name>spreadfunc</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a5f4c05b95c0bf377dc7b9beb721cae0c</anchor>
      <arglist>)(int, Int_t, Int_t, Double_t *)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t(NBody::KDTree::*</type>
      <name>entropyfunc</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aa47732e9da2599bb9f0124f485a77cd9</anchor>
      <arglist>)(int, Int_t, Int_t, Double_t, Double_t, Double_t, Double_t *)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t(NBody::KDTree::*</type>
      <name>medianfunc</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a8f098b1e73e30dc428d9c95389d3164b</anchor>
      <arglist>)(int, Int_t, Int_t, Int_t, bool)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>spreada</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>af1b5328053cedb390a4a43389816f335</anchor>
      <arglist>[MAXND]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>meana</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a67ababedc69ce2809e4ee462fae4f662</anchor>
      <arglist>[MAXND]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>vara</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a109a9aaae0f68b6ef32179c46d3b9521</anchor>
      <arglist>[MAXND]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>entropya</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ae2317c1ff27a2983a8009aaa29f83743</anchor>
      <arglist>[MAXND]</arglist>
    </member>
    <member kind="variable" protection="private" static="yes">
      <type>static const int</type>
      <name>MAXND</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ab1ca2d1ca718bc77603f7c57abf0b8f1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>TPHYS</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a3a37ad5dc668deb4dcb30a73ddb4ffaf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>TPROJ</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a74df3e8aa3f76cc6612f17e7dcf3a187</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>TVEL</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a2514a82816fafb1a61112647ecb8edd5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>TPHS</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a2bd32605634282709d8f3a6e0474b9bc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>TMETRIC</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a053b0fd24c3261b441c4776ac84a33c2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>KSPH</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a583a85ca3b644b4243d0243ff43f5b10</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>KGAUSS</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a6a2e3bc67ca903e47589cc74cf061207</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>KEPAN</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a14680ea01d84c82a930ae484ce983fe3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>KTH</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aa73ae200687c8a638e5d241bb5524d9b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t(NBody::KDTree::*</type>
      <name>bmfunc</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a5e6fdd1be0fb12c2648ef9e844763e2a</anchor>
      <arglist>)(int, Int_t, Int_t, Double_t *)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t(NBody::KDTree::*</type>
      <name>dispfunc</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>af73cb9c00df90d290e7cdb32456ed59a</anchor>
      <arglist>)(int, Int_t, Int_t, Double_t)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t(NBody::KDTree::*</type>
      <name>spreadfunc</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a5f4c05b95c0bf377dc7b9beb721cae0c</anchor>
      <arglist>)(int, Int_t, Int_t, Double_t *)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t(NBody::KDTree::*</type>
      <name>entropyfunc</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aa47732e9da2599bb9f0124f485a77cd9</anchor>
      <arglist>)(int, Int_t, Int_t, Double_t, Double_t, Double_t, Double_t *)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t(NBody::KDTree::*</type>
      <name>medianfunc</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a8f098b1e73e30dc428d9c95389d3164b</anchor>
      <arglist>)(int, Int_t, Int_t, Int_t, bool)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>spreada</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>af1b5328053cedb390a4a43389816f335</anchor>
      <arglist>[MAXND]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>meana</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a67ababedc69ce2809e4ee462fae4f662</anchor>
      <arglist>[MAXND]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>vara</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a109a9aaae0f68b6ef32179c46d3b9521</anchor>
      <arglist>[MAXND]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>entropya</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ae2317c1ff27a2983a8009aaa29f83743</anchor>
      <arglist>[MAXND]</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Node *</type>
      <name>BuildNodes</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a8247f166032c6d40a1034b4adcdde96d</anchor>
      <arglist>(Int_t start, Int_t end)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>ScaleSpace</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>afb34afbe30171832d3316c20315671ed</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>int</type>
      <name>TreeTypeCheck</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a75d9ae25e0169f9878c1a10bb5151057</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>KernelConstruction</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a885b4b998c0d3f70e000ca34f09b9169</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>KDTree</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad68f262da92174fcc384a9d5a8774a40</anchor>
      <arglist>(Particle *p, Int_t numparts, Int_t bucket_size=16, int TreeType=TPHYS, int KernType=KEPAN, int KernRes=1000, int SplittingCriterion=0, int Aniso=0, int ScaleSpace=0, Double_t *Period=NULL, Double_t **metric=NULL)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>KDTree</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad34cea8c3ba5f18de5c424072f969b29</anchor>
      <arglist>(System &amp;s, Int_t bucket_size=16, int TreeType=TPHYS, int KernType=KEPAN, int KernRes=1000, int SplittingCriterion=0, int Aniso=0, int ScaleSpace=0, Double_t **metric=NULL)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~KDTree</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aa2abf2018780d1be3d52f44c4899e52c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetNumNodes</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a95891d021490d3f805e6a9e97ecb6c6d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetNumLeafNodes</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a32bff63406186f876844783e6fd07835</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetBucketSize</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a94a2c83edc5de03b3804e1bd5d82818f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetTreeType</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ab23fa720d57157f80d861c3e238ed27e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetKernType</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a6c7ac9d27697119a0526c2cd921fbf32</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetKernNorm</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a0bb7ab867cd4335ad859c5d995ccdf10</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Node *</type>
      <name>GetRoot</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad60c3768328d7e2153b5ac3cdaf0fc55</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetPeriod</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aa3b489bb80090cde772713f02d7bb07d</anchor>
      <arglist>(int j)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearest</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ae3c79e5fb889d4a4d14ecc0ae8374f7f</anchor>
      <arglist>(Int_t **nn, Double_t **dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>add9c46bd3e5a4bfa8b0402c876e903dd</anchor>
      <arglist>(Int_t **nn, Double_t **dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ab3cf73fd48ae33e6d9530e8c054386a5</anchor>
      <arglist>(Int_t **nn, Double_t **dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>afeb98226c436f301b68dcfe5a9311c4e</anchor>
      <arglist>(Int_t **nn, Double_t **dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a845d0a6df00755a297ea34eb3b07f747</anchor>
      <arglist>(FOFcompfunc cmp, Double_t *params, Int_t **nn, Double_t **dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearest</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a1d249268e725c375ae50eb2a6dcec7c0</anchor>
      <arglist>(Int_t tt, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a52a3be0cfc421225695c34f330f6172c</anchor>
      <arglist>(Int_t tt, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a23098f6033aff6729bc92c04edf2a8f3</anchor>
      <arglist>(Int_t tt, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a35ed1f74def57d237ff518244144090a</anchor>
      <arglist>(Int_t tt, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearest</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a5b3ef1437766578990c3074cd1536f85</anchor>
      <arglist>(Double_t *x, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad457bf15291b9a388113f13303e98a2b</anchor>
      <arglist>(Double_t *x, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ae7a94c4106653116bf53ca7a7619ea25</anchor>
      <arglist>(Double_t *v, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a633aa4e77cf2ff05eb47364bf2cf430b</anchor>
      <arglist>(Double_t *x, Double_t *v, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a410b409750426b307aef5b24ecc88958</anchor>
      <arglist>(Coordinate x, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ae8ac7e4c970c59c1ba9d2d3804de0026</anchor>
      <arglist>(Coordinate v, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a73b1b578760cb174e9617fb5efb1362a</anchor>
      <arglist>(Coordinate x, Coordinate v, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a7e0a8bd0362b3c478c72034755cf6aea</anchor>
      <arglist>(Int_t tt, FOFcompfunc cmp, Double_t *params, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>af98daa0ca7fcc6fb12c409b6761ea20a</anchor>
      <arglist>(Particle p, FOFcompfunc cmp, Double_t *params, Int_t *nn, Double_t *dist2, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBall</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a86888b7be77da9ad42691cd87db929a5</anchor>
      <arglist>(Int_t tt, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a57f312f21caa5aa1316561d1db1ec217</anchor>
      <arglist>(Int_t tt, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a57b980a39e4c0741abc70abf74635a01</anchor>
      <arglist>(Int_t tt, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPhase</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ac69e103e8b8da8bc02591caa7eccb135</anchor>
      <arglist>(Int_t tt, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBall</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a3a9113a3650d488f4176b3f942db1d7c</anchor>
      <arglist>(Double_t *x, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a452e2431da28b8753bcf7af0e41a2bae</anchor>
      <arglist>(Double_t *x, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a988f3e538c2337d4a49254824cc7dd25</anchor>
      <arglist>(Double_t *v, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPhase</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad0d4bb038fd1cde03c2bee06867fa9c9</anchor>
      <arglist>(Double_t *x, Double_t *v, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBall</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a126899aaba1fa6a5ddac9d20ed2b4e2d</anchor>
      <arglist>(Coordinate x, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a372afa16b3f280268d7025e8e6f360b0</anchor>
      <arglist>(Coordinate x, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ae0b83f4d4bcb3f251b93d68c3ece3ff8</anchor>
      <arglist>(Coordinate v, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPhase</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a3bb62c923f6a8c2c265d012cc53b6a5b</anchor>
      <arglist>(Coordinate x, Coordinate v, Double_t fdist2, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterion</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a3d4560dd709924c2b326fa375271c8fa</anchor>
      <arglist>(Int_t tt, FOFcompfunc p, Double_t *params, Int_t imark, Int_t *nn, Double_t *dist2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterion</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a63d820687992c468f0890a77ea2da3b1</anchor>
      <arglist>(Int_t tt, FOFcompfunc p, Double_t *params, Int_t imark, Int_t *nn)</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a0a70772eb802c0fc39566d8da80a23af</anchor>
      <arglist>(Int_t tt, Double_t fdist2, Int_t *tagged)</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ab9cde2cb2f1b101a8073891a7a0c9a9b</anchor>
      <arglist>(Coordinate x, Double_t fdist2, Int_t *tagged)</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a4480b407eb873fdb0ae977b911e27251</anchor>
      <arglist>(Double_t *x, Double_t fdist2, Int_t *tagged)</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>SearchCriterionTagged</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a870da6c7c1bae8e5babc3a507d2eb4a8</anchor>
      <arglist>(Int_t tt, FOFcompfunc cmp, Double_t *params, Int_t *tagged)</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>SearchCriterionTagged</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a91d66c4b62c2ba71e80791c310f9b795</anchor>
      <arglist>(Particle &amp;p, FOFcompfunc cmp, Double_t *params, Int_t *tagged)</arglist>
    </member>
    <member kind="function">
      <type>Node *</type>
      <name>FindLeafNode</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a0cc5ead98ebb862fdf506ccb9d00277a</anchor>
      <arglist>(Int_t tt)</arglist>
    </member>
    <member kind="function">
      <type>Node *</type>
      <name>FindLeafNode</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a5586bfa61c5d45993a95cb475558e547</anchor>
      <arglist>(Double_t *x)</arglist>
    </member>
    <member kind="function">
      <type>Node *</type>
      <name>FindLeafNode</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a2a2a13a32a18b93a2db01b93d24d3fe4</anchor>
      <arglist>(Double_t search[6][2])</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>CalcDensity</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad4e03bf55723a599850fe3e15af5add2</anchor>
      <arglist>(Int_t Nsmooth=64)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>CalcVelDensity</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a3e42b28de44aad77516d091fe72c053b</anchor>
      <arglist>(Int_t Nsmooth=64, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>Double_t *</type>
      <name>CalcVelDensityWithPhysDensity</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a2e7ea3a186de7baf943bead4e1968a67</anchor>
      <arglist>(Int_t Nsmooth=64, Int_t Nsearch=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate *</type>
      <name>CalcSmoothVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>af25297a80fe72a8396c74266e19c6d5f</anchor>
      <arglist>(Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Matrix *</type>
      <name>CalcSmoothVelDisp</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a97d4ba62011bf08948ac272217d1313d</anchor>
      <arglist>(Coordinate *smvel, Int_t Nsmooth=64, int densityset=1, int meanvelset=1)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate *</type>
      <name>CalcSmoothVelSkew</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>af0cd9c8fdadc3ac92b941aea719a8b9d</anchor>
      <arglist>(Coordinate *smvel, Matrix *smveldisp, Int_t Nsmooth=64, int densityset=1, int meanvelset=1, int veldispset=1)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate *</type>
      <name>CalcSmoothVelKurtosis</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ae058f09de1d4a983288c5c5f02375c1e</anchor>
      <arglist>(Coordinate *smvel, Matrix *smveldisp, Int_t Nsmooth=64, int densityset=1, int meanvelset=1, int veldispset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcDensityParticle</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ae07ac97de0b123ecfe2e2e81fc707524</anchor>
      <arglist>(Int_t target, Int_t Nsmooth=64)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcVelDensityParticle</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a44a92a218c9601df2f9cd61e27771000</anchor>
      <arglist>(Int_t target, Int_t Nsmooth=64, Int_t Nsearch=64, int iflag=0, PriorityQueue *pq=NULL, PriorityQueue *pq2=NULL, Int_t *nnIDs=NULL, Double_t *vdist=NULL)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcVelDensityWithPhysDensityParticle</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a52a014e7b489f45461b64218abc2c205</anchor>
      <arglist>(Int_t target, Int_t Nsmooth=64, Int_t Nsearch=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>CalcSmoothVelParticle</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aa2a21e48da4f5283c109d2c852fd91d2</anchor>
      <arglist>(Int_t target, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>CalcSmoothVelDispParticle</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a226fe168e0c666e33b725e97346dd2f5</anchor>
      <arglist>(Int_t target, Coordinate smvel, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcDensityPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a7beae71edc364f8fc4d9d4d3748340cc</anchor>
      <arglist>(Double_t *x, Int_t Nsmooth=64, Double_t *v=NULL)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcVelDensityPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a9db41f8594c17279ae6f3313e7006ec2</anchor>
      <arglist>(Double_t *x, Double_t *v, Int_t Nsmooth=64, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>CalcSmoothVelPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a7380400e1847af4a03db178219bec168</anchor>
      <arglist>(Double_t *x, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>CalcSmoothVelDispPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a72e6fcdc1e7febb1580980e658be7cb4</anchor>
      <arglist>(Double_t *x, Coordinate smvel, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcDensityPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad522decde9a620c9611a6784487de3ea</anchor>
      <arglist>(Coordinate x, Int_t Nsmooth=64, Coordinate v=Coordinate(0.))</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcVelDensityPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aa65db49c3b6fbb0cefc50c22764c9c7b</anchor>
      <arglist>(Coordinate x, Int_t Nsmooth=64, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>CalcSmoothVelPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aac3058c4b52c702c4646b7acb4994048</anchor>
      <arglist>(Coordinate x, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>CalcSmoothVelDispPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a6545172972f2e9e0d3e000b19faba207</anchor>
      <arglist>(Coordinate x, Coordinate smvel, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcDensityPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aeae8308d507973f11f7e4f486a91b8a8</anchor>
      <arglist>(Real_t *x, Int_t Nsmooth=64, Real_t *v=NULL)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcVelDensityPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>afa3ece24a460a6bfe6476aa5879825b7</anchor>
      <arglist>(Real_t *x, Real_t *v, Int_t Nsmooth=64, Int_t Nsearch=64)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>CalcSmoothVelPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a559378eedd2b016e287ab4c6966ce9a1</anchor>
      <arglist>(Real_t *x, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>CalcSmoothVelDispPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad382ad3084e42eecf474adc5f1e300ab</anchor>
      <arglist>(Real_t *x, Coordinate smvel, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t *</type>
      <name>CalcSmoothLocalMean</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a7f6c6ac46b636d85c7b5cc0bb65297fc</anchor>
      <arglist>(Double_t *weight, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t *</type>
      <name>CalcSmoothLocalDisp</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a6dd9b6a26aa7c7208fc2f7a82d510bec</anchor>
      <arglist>(Double_t *weight, Double_t *localmean, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcSmoothLocalMeanParticle</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a9631ff5badc7d042c9dbaee09013f8b8</anchor>
      <arglist>(Int_t target, Double_t *weight, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcSmoothLocalDispParticle</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>adef2f3a2bb53461f94bbecb16512c5ff</anchor>
      <arglist>(Int_t target, Double_t *weight, Double_t localmean, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcSmoothLocalMeanPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a669096967ed3e8386c1b9b886eff3180</anchor>
      <arglist>(Double_t *x, Double_t *weight, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcSmoothLocalDispPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a6b007b5d7df17540b24b3b69e6eb6b31</anchor>
      <arglist>(Double_t *x, Double_t *weight, Double_t localmean, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcSmoothLocalMeanPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a516c1ea1abd8060c8ab0fc39cce9e79e</anchor>
      <arglist>(Coordinate x, Double_t *weight, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcSmoothLocalDispPosition</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>abb91b843b890016b4223fab569314143</anchor>
      <arglist>(Coordinate x, Double_t *weight, Double_t localmean, Int_t Nsmooth=64, int densityset=1)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CalcSmoothLocalValue</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a512192c8bd3d50bfd76056860b3015e9</anchor>
      <arglist>(Int_t Nsmooth, PriorityQueue *pq, Double_t *weight)</arglist>
    </member>
    <member kind="function">
      <type>Int_t *</type>
      <name>FOF</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a178903f3a706d9d36085d044396a35c7</anchor>
      <arglist>(Double_t fdist, Int_t &amp;numgroup, Int_t minnum=8, int order=0, Int_tree_t *pHead=NULL, Int_tree_t *pNext=NULL, Int_tree_t *pTail=NULL, Int_tree_t *pLen=NULL)</arglist>
    </member>
    <member kind="function">
      <type>Int_t *</type>
      <name>FOFCriterion</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a54433b1636f0e5c75bfa9b27fa1bc66c</anchor>
      <arglist>(FOFcompfunc p, Double_t *params, Int_t &amp;numgroups, Int_t minnum=8, int order=0, int ipcheckflag=0, FOFcheckfunc check=Pnocheck, Int_tree_t *pHead=NULL, Int_tree_t *pNext=NULL, Int_tree_t *pTail=NULL, Int_tree_t *pLen=NULL)</arglist>
    </member>
    <member kind="function">
      <type>Int_t *</type>
      <name>FOFCriterionSetBasisForLinks</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a68af3ba6a95b7c3b6a2d9d743d4ebc00</anchor>
      <arglist>(FOFcompfunc cmp, Double_t *params, Int_t &amp;numgroup, Int_t minnum=8, int order=0, int ipcheckflag=0, FOFcheckfunc check=Pnocheck, Int_tree_t *pHead=NULL, Int_tree_t *pNext=NULL, Int_tree_t *pTail=NULL, Int_tree_t *pLen=NULL)</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>FOFCriterionParticle</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a80277ea2f0337cba28616243685fb279</anchor>
      <arglist>(FOFcompfunc p, Int_t *pfof, Int_t target, Int_t iGroup, Double_t *params, Int_tree_t *pGroupHead, Int_tree_t *Fifo, Int_tree_t *pHead, Int_tree_t *pTail, Int_tree_t *pNext, Int_tree_t *pLen)</arglist>
    </member>
    <member kind="function">
      <type>Int_t *</type>
      <name>FOFNNCriterion</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a2f1a25ef3f1f8eca6c51c660703e26da</anchor>
      <arglist>(FOFcompfunc p, Double_t *params, Int_t numNN, Int_t **nnIDs, Int_t &amp;numgroups, Int_t minnum=8)</arglist>
    </member>
    <member kind="function">
      <type>Int_t *</type>
      <name>FOFNNDistCriterion</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a8172fbe3a807e627e3410027a90b0d99</anchor>
      <arglist>(FOFcompfunc p, Double_t *params, Int_t numNN, Int_t **nnIDs, Double_t **dist2, Double_t disfunc(Int_t, Double_t *), Int_t npc, Int_t *npca, Int_t &amp;numgroups, Int_t minnum=8)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>SpreadestPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ae05333d352d876ce064a9367cc75e765</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t *bnd)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>SpreadestVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a7d1f24e2981cc26899098296bd140111</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t *bnd)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>SpreadestPhs</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad3ab962356ec6265931e8b0e13118692</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t *bnd)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>BoundaryandMeanPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a5cf924f1a6bd70d8879564ff59670377</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t *bnd)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>BoundaryandMeanVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aa64d34d16856c8cc13344ef009378e56</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t *bnd)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>BoundaryandMeanPhs</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a28fec631a862b83783aeff72af8b685f</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t *bnd)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>DispersionPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>aa6e08bf249a5aa09acf928c37035f8d0</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t mean)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>DispersionVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a132182b036ea58425562684ccc217a2b</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t mean)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>DispersionPhs</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a5978085e40e951636edf0a708e46e62c</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t mean)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>EntropyPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a1dbab902a9cc099f8f09e64c9fc583f4</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t low, Double_t up, Double_t nbins, Double_t *ni)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>EntropyVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a88c01aa506632afbd9dc059c8ad08906</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t low, Double_t up, Double_t nbins, Double_t *ni)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>EntropyPhs</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a19e7466f4ccefd6f5e1cd62e80f6217d</anchor>
      <arglist>(int j, Int_t start, Int_t end, Double_t low, Double_t up, Double_t nbins, Double_t *ni)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>MedianPos</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a75c54d1a61a27ce5fe052294c24e3487</anchor>
      <arglist>(int d, Int_t k, Int_t start, Int_t end, bool balanced=true)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>MedianVel</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a518800579c8aeff8a136d84da504ac4b</anchor>
      <arglist>(int d, Int_t k, Int_t start, Int_t end, bool balanced=true)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>MedianPhs</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ad625f0d10fdb7c06aac10b99131c085d</anchor>
      <arglist>(int d, Int_t k, Int_t start, Int_t end, bool balanced=true)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>Double_t</type>
      <name>Wsm</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>af702721d847b96b922092986ae970ac8</anchor>
      <arglist>(Double_t r, Int_t i, Int_t size, Double_t delta, Double_t *x)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>CalculateMetricSpacing</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a3db4480b0e0ad6215cd33507c00064ab</anchor>
      <arglist>(Int_t target, int treetype, Double_t *smetric)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>CalculateMetricSpacing</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a804d5e2cf646615ec4cb673f1e757049</anchor>
      <arglist>(const Double_t *x, const Double_t *v, int treetype, Double_t *smetric)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>CalculateMetricSpacing</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>abe82debeb779bf53f4817f72ebd4a0a2</anchor>
      <arglist>(const Real_t *x, const Real_t *v, int treetype, Double_t *smetric)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>CalculateMetricTensor</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a4713a6decf855de60a6e50e326216e7f</anchor>
      <arglist>(Int_t target, int treetype, Double_t *smetric, Double_t *metric, GMatrix &amp;gmetric)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>CalculateMetricTensor</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a8df252b1017292b815a085e2628721a0</anchor>
      <arglist>(const Double_t *x, const Double_t *v, int treetype, Double_t *smetric, Double_t *metric, GMatrix &amp;gmetric)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>CalculateMetricTensor</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>ace41c3032fb804dc2417d278cfb10713</anchor>
      <arglist>(const Real_t *x, const Real_t *v, int treetype, Double_t *smetric, Double_t *metric, GMatrix &amp;gmetric)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>LoadNN</name>
      <anchorfile>classNBody_1_1KDTree.html</anchorfile>
      <anchor>a033a72ca0920c26739451fabf8512eb0</anchor>
      <arglist>(const Int_t ns, PriorityQueue *pq, Int_t *nn, Double_t *dist)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>NBody::LeafNode</name>
    <filename>classNBody_1_1LeafNode.html</filename>
    <base>NBody::Node</base>
    <member kind="function">
      <type></type>
      <name>LeafNode</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>aa01ec2361f1dc2f4b6125c2b7378e3f6</anchor>
      <arglist>(Int_t id, Int_t new_bucket_start, Int_t new_bucket_end, Double_t bnd[6][2], unsigned short ndim)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~LeafNode</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a165aab7f298bfc6217af4e0dedc79746</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a81d995cdca80025ce90eb82b804c60a8</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a13ed52efa6d980028d20f50d64d57084</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ab7036e0fc27ee7fcfd3d341d3b8c3527</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ab3cade1bc2fe791b371e2b0f44fec971</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>add29a894675fdf371148cc2ef37b1a85</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *v, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ab9646159f5e1d396103797c7b4344cc2</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a3642b342933b319cbe046e0b028e3dc3</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a78acd03eccbf459658ffa2d8fb4f55e5</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate v, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ae74b48150fde6bc1bb8e018a74553047</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a4812017e4d78bd112c249df8f6e48b80</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ae3453759822d2473c73071efbbc4476d</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ad9522b3daa0f2b1b585d836482739b47</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>acd2c38b7186f0847eac4dfa2cfab8d69</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ab400ad1f179a38406865685d923aba41</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ace62b1ef6f826b80778c26ab71bb29e4</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ae2c1e95b14e9ee68c991b05c4655bb7e</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a735d31b34bf84be06b7a0bf4de4276e2</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>aa4c20639dc116b5013cb50442f2b333f</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a916f22ac3b1c82a6e8e910e5ef1d3df9</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a17dcd482b4d5be861112bd360549d7e5</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Coordinate x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ad6f971b760418ea6ab14617b42dc9e24</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Int_t t, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a66b2469cce2a107bc0596339622b843d</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *x, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>aa73fe7db4d12135d86fc4ce773800631</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Coordinate x, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterion</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a6e414a1405fdaa9dc6d530f562b0820e</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterion</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a087e12cda66c24bf1e8deae27a1b43f4</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>abf127aff0478ff9e5dbf5d62b36af1a6</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a5ef9fce165b5f79aa70145422d627c19</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionNoDist</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a5ccaed27c5f0961fe396b216490bf548</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionNoDist</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>adfe8039f4f981bf9998eab7ae0ae590e</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchBall</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ac21ea4f55e5dfc6661fd5816c163d48d</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchCriterion</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a93b2531a65d43ff833129cc520b68004</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchCriterionSetBasisForLinks</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a894080d638e68cebd846d60cc05dfdfa</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ad546bb9a16996653517e2ec631645e9a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>aa9525299d9445b1706fffe97a2991d7a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a01d818deef9ff1e438051d5fc424536e</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>af6b2aeecd62667d780bc33d1110ed523</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a195642e50da39a2d6d2759e40aa44bf6</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *v, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ad675315c973e8c956052010430c3b255</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a8263166af4e01f84b35daeb8c50c6d9d</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ac14b20c7f989f9f6d0a8101f6f5c6900</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate v, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>aab5a452bfa386510411c16f0a16a749a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a72f689f1d28dd5e8a88b15da30043936</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ae889d53c06b0857304b824570aca72da</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>aaf01dca60ea2914f7e15733f7bddefbb</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ab3a2f4291a0e02866d74ec357a577c04</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>af3e3f6ccffef06cb33b653656f6978aa</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ac67f328427b43474862aa9281f9a659e</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterionPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a2cb15d47d2dd9d9f0c0c2bb8a86d3ea3</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterionPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a5d9de7f4cec82891cd83f982b8eaec92</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>aed873611259322032115dd3fae332f1f</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a6584ae9ff17e7ff9c651d435a44ae9f1</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Double_t *x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a1113eb1ddef3e1403efc29b6205035c8</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Coordinate x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>af9b10d964d021cd40f20a3611aa4a48a</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Int_t t, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ae10eaf63495e908ef5ac3aa65663512d</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Double_t *x, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>aca58eb05f6656385f3ebccf438a159a7</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Coordinate x, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a9a587a1c355ab72b6ae7f580612470c1</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a78e040de8107035e30b0a6c2898ad08d</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionPeriodicTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>af702a0204b27db907c5b869ceb41ca48</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionPeriodicTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a810c466853c4d27f4c461c3569efc3bb</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionNoDistPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a38cbaba2cdbd2bce2da8ee11110411ea</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionNoDistPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ab9bab8b9767f47747ee020fc192be982</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchBallPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ac10f5bf8e27eed9a42ac5f3055cf2172</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a537244b9935c02b28f72acaad3ac547b</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchCriterionSetBasisForLinksPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>aa33b8136ac8c4f76f48f98fcd3cd647f</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a81d995cdca80025ce90eb82b804c60a8</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a13ed52efa6d980028d20f50d64d57084</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ab7036e0fc27ee7fcfd3d341d3b8c3527</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ab3cade1bc2fe791b371e2b0f44fec971</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>add29a894675fdf371148cc2ef37b1a85</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *v, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ab9646159f5e1d396103797c7b4344cc2</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a3642b342933b319cbe046e0b028e3dc3</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a78acd03eccbf459658ffa2d8fb4f55e5</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate v, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ae74b48150fde6bc1bb8e018a74553047</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a4812017e4d78bd112c249df8f6e48b80</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ae3453759822d2473c73071efbbc4476d</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ad9522b3daa0f2b1b585d836482739b47</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>acd2c38b7186f0847eac4dfa2cfab8d69</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ab400ad1f179a38406865685d923aba41</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ace62b1ef6f826b80778c26ab71bb29e4</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ae2c1e95b14e9ee68c991b05c4655bb7e</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a735d31b34bf84be06b7a0bf4de4276e2</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>aa4c20639dc116b5013cb50442f2b333f</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a916f22ac3b1c82a6e8e910e5ef1d3df9</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a17dcd482b4d5be861112bd360549d7e5</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Coordinate x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ad6f971b760418ea6ab14617b42dc9e24</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Int_t t, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a66b2469cce2a107bc0596339622b843d</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *x, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>aa73fe7db4d12135d86fc4ce773800631</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Coordinate x, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterion</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a6e414a1405fdaa9dc6d530f562b0820e</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterion</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a087e12cda66c24bf1e8deae27a1b43f4</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>abf127aff0478ff9e5dbf5d62b36af1a6</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a5ef9fce165b5f79aa70145422d627c19</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionNoDist</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a5ccaed27c5f0961fe396b216490bf548</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionNoDist</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>adfe8039f4f981bf9998eab7ae0ae590e</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchBall</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ac21ea4f55e5dfc6661fd5816c163d48d</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchCriterion</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a93b2531a65d43ff833129cc520b68004</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchCriterionSetBasisForLinks</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a894080d638e68cebd846d60cc05dfdfa</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ad546bb9a16996653517e2ec631645e9a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>aa9525299d9445b1706fffe97a2991d7a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a01d818deef9ff1e438051d5fc424536e</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>af6b2aeecd62667d780bc33d1110ed523</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a195642e50da39a2d6d2759e40aa44bf6</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *v, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ad675315c973e8c956052010430c3b255</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a8263166af4e01f84b35daeb8c50c6d9d</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ac14b20c7f989f9f6d0a8101f6f5c6900</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate v, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>aab5a452bfa386510411c16f0a16a749a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a72f689f1d28dd5e8a88b15da30043936</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ae889d53c06b0857304b824570aca72da</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>aaf01dca60ea2914f7e15733f7bddefbb</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ab3a2f4291a0e02866d74ec357a577c04</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>af3e3f6ccffef06cb33b653656f6978aa</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ac67f328427b43474862aa9281f9a659e</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterionPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a2cb15d47d2dd9d9f0c0c2bb8a86d3ea3</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterionPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a5d9de7f4cec82891cd83f982b8eaec92</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>aed873611259322032115dd3fae332f1f</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a6584ae9ff17e7ff9c651d435a44ae9f1</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Double_t *x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a1113eb1ddef3e1403efc29b6205035c8</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Coordinate x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>af9b10d964d021cd40f20a3611aa4a48a</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Int_t t, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ae10eaf63495e908ef5ac3aa65663512d</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Double_t *x, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>aca58eb05f6656385f3ebccf438a159a7</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Coordinate x, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a9a587a1c355ab72b6ae7f580612470c1</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a78e040de8107035e30b0a6c2898ad08d</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionPeriodicTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>af702a0204b27db907c5b869ceb41ca48</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionPeriodicTagged</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a810c466853c4d27f4c461c3569efc3bb</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionNoDistPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a38cbaba2cdbd2bce2da8ee11110411ea</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionNoDistPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ab9bab8b9767f47747ee020fc192be982</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchBallPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>ac10f5bf8e27eed9a42ac5f3055cf2172</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>a537244b9935c02b28f72acaad3ac547b</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchCriterionSetBasisForLinksPeriodic</name>
      <anchorfile>classNBody_1_1LeafNode.html</anchorfile>
      <anchor>aa33b8136ac8c4f76f48f98fcd3cd647f</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>Math::math_function</name>
    <filename>structMath_1_1math__function.html</filename>
    <member kind="variable">
      <type>Double_t(*</type>
      <name>function</name>
      <anchorfile>structMath_1_1math__function.html</anchorfile>
      <anchor>a5f37bfcb60d52aebc0f7b173063f788d</anchor>
      <arglist>)(Double_t x, void *params)</arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>params</name>
      <anchorfile>structMath_1_1math__function.html</anchorfile>
      <anchor>accf92edda842a1f6fc91757381992c71</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>Math::math_multidim_function</name>
    <filename>structMath_1_1math__multidim__function.html</filename>
    <member kind="variable">
      <type>Double_t(*</type>
      <name>function</name>
      <anchorfile>structMath_1_1math__multidim__function.html</anchorfile>
      <anchor>a721b8b820cc5e42adc510edb49ddcb68</anchor>
      <arglist>)(Double_t *x, int ndim, void *params)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>ndim</name>
      <anchorfile>structMath_1_1math__multidim__function.html</anchorfile>
      <anchor>adee5590a123b7cc39ca94c61747bb79c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>void *</type>
      <name>params</name>
      <anchorfile>structMath_1_1math__multidim__function.html</anchorfile>
      <anchor>af523b7ef989e4de3e6738499a4bc981b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Math::Matrix</name>
    <filename>classMath_1_1Matrix.html</filename>
    <member kind="function">
      <type></type>
      <name>Matrix</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a131174c997e4afa01160d3c74fb60285</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Matrix</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a5917235e695f303c76991499f1b00ade</anchor>
      <arglist>(Double_t a)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Matrix</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a40a49df40afb5ea03233160fbd9ff820</anchor>
      <arglist>(Double_t a, Double_t b, Double_t c, Double_t d, Double_t e, Double_t f, Double_t g, Double_t h, Double_t i)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Matrix</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>ae5dd1dcd6dbb98bf0b061574b738d038</anchor>
      <arglist>(const Double_t m[9])</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Matrix</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a4765978886227259ef7d2d4691ba9d9b</anchor>
      <arglist>(const Matrix &amp;m)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>matrix</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a8a9e0b87049adb51ce709562462d4e5f</anchor>
      <arglist>[9]</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Matrix</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a131174c997e4afa01160d3c74fb60285</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Matrix</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a5917235e695f303c76991499f1b00ade</anchor>
      <arglist>(Double_t a)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Matrix</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a40a49df40afb5ea03233160fbd9ff820</anchor>
      <arglist>(Double_t a, Double_t b, Double_t c, Double_t d, Double_t e, Double_t f, Double_t g, Double_t h, Double_t i)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Matrix</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>ae5dd1dcd6dbb98bf0b061574b738d038</anchor>
      <arglist>(const Double_t m[9])</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Matrix</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a4765978886227259ef7d2d4691ba9d9b</anchor>
      <arglist>(const Matrix &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>Double_t &amp;</type>
      <name>operator()</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a3af4c8761c8774f1e7a43d5d3f131f17</anchor>
      <arglist>(int i, int j)</arglist>
    </member>
    <member kind="function">
      <type>const Double_t &amp;</type>
      <name>operator()</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a6d7387d92b1ea3128d332b587b694878</anchor>
      <arglist>(int i, int j) const </arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>operator+</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>ab53cb55fc29a3a4234e8d4a373684c5a</anchor>
      <arglist>(const Matrix &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>operator-</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a197f6bd1196d39b9bb36c34d728251cc</anchor>
      <arglist>(const Matrix &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>operator*</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a046ca1f5774405ff1680f03a1a556ee4</anchor>
      <arglist>(Double_t a)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>operator*</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a4009014cea7b94bdb7e4bcc1fd9c31f1</anchor>
      <arglist>(const Coordinate &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>operator*</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a7c4d33826b1250d44ea76d7d65294843</anchor>
      <arglist>(const Matrix &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>Matrix &amp;</type>
      <name>operator=</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a85791e27dd275cb7f51b65b4eb353ec4</anchor>
      <arglist>(const Matrix m)</arglist>
    </member>
    <member kind="friend">
      <type>friend Matrix</type>
      <name>operator*</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a628bd0363ac2e4a63da36953ef048d03</anchor>
      <arglist>(Double_t a, const Matrix &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Trace</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a75af96132f172a1fb50a68183a557ede</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>Transpose</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a8b0179a7f462dbda7af679a944c8f7d2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>TransposeInPlace</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>aa3793f2447c338935c082477101ae43f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Det</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>abaaaa4a0eddb178254b70204f9363860</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>Adjugate</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a4eb2f17e5ddc5de92e8fc5d133950081</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>LUDecomp</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a42e42a2eff11c1bc0cf1b2f8ff998517</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>Inverse</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a6c84f3fece37fb5ef571c04c6f688828</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>Eigenvalues</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a92a6d9a4d702d83bee33fdb2cd3a3c43</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>Eigenvectors</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>a081983835830f3b025fd7d8c9a05895c</anchor>
      <arglist>(const Coordinate &amp;e) const </arglist>
    </member>
    <member kind="friend">
      <type>friend std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classMath_1_1Matrix.html</anchorfile>
      <anchor>aec0a552024111f8cfa5233c531431e66</anchor>
      <arglist>(std::ostream &amp;stream, const Matrix &amp;m)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>Math::Matrix2D</name>
    <filename>classMath_1_1Matrix2D.html</filename>
    <member kind="function">
      <type></type>
      <name>Matrix2D</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a70d5c96135271649eb5b6ffb71a23967</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Matrix2D</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a599f008a46d18c8adeaa822540f1cdd1</anchor>
      <arglist>(Double_t a)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Matrix2D</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>aeae19c0e29aff5334a5aa5358c51f4e8</anchor>
      <arglist>(Double_t a, Double_t b, Double_t c, Double_t d)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Matrix2D</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a97a76813389cd9783bf85d2be34b8fd4</anchor>
      <arglist>(Double_t m[4])</arglist>
    </member>
    <member kind="function">
      <type>Matrix2D</type>
      <name>Transpose</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>aca866eeebcea903d982dbeebe69ac461</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>TransposeInPlace</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a3dc16a8e65878c113482340486a3fd77</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Det</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a6c8aeb365642710d259e232edee224c1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Matrix2D</type>
      <name>Adjugate</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>ac91c76081f8da74ca83cc9481b711961</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Matrix2D</type>
      <name>Inverse</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a86c5d6cd4d53890f0642404797209581</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Coordinate2D</type>
      <name>Eigenvalues</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>aaa664bc1ad317b6aee9201f5fc39012b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Matrix2D</type>
      <name>Eigenvectors</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>aa172365cbf934d7446c9ec8b04e3a7d9</anchor>
      <arglist>(const Coordinate2D &amp;e) const </arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>matrix2d</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a4464a37be49dcc3f8bee45f6e946f672</anchor>
      <arglist>[4]</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Matrix2D</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a70d5c96135271649eb5b6ffb71a23967</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Matrix2D</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a599f008a46d18c8adeaa822540f1cdd1</anchor>
      <arglist>(Double_t a)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Matrix2D</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>aeae19c0e29aff5334a5aa5358c51f4e8</anchor>
      <arglist>(Double_t a, Double_t b, Double_t c, Double_t d)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Matrix2D</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a97a76813389cd9783bf85d2be34b8fd4</anchor>
      <arglist>(Double_t m[4])</arglist>
    </member>
    <member kind="function">
      <type>Double_t &amp;</type>
      <name>operator()</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>aad0ff64d193cc6a8ac471724ce040770</anchor>
      <arglist>(int i, int j)</arglist>
    </member>
    <member kind="function">
      <type>const Double_t &amp;</type>
      <name>operator()</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a07ce871bec31fe044dae5031222e316b</anchor>
      <arglist>(int i, int j) const </arglist>
    </member>
    <member kind="function">
      <type>Matrix2D</type>
      <name>operator+</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>aa7be4718c254184d10cee9f2c0be3456</anchor>
      <arglist>(const Matrix2D &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>Matrix2D</type>
      <name>operator-</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a03664d0a16e15a3cf14e6cbd0b1ef4ad</anchor>
      <arglist>(const Matrix2D &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>Matrix2D</type>
      <name>operator*</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a8d62e9fb2abc0405d8a3faba7d8f4a9a</anchor>
      <arglist>(Double_t a)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate2D</type>
      <name>operator*</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>addb31bd02729d0c4e816f2ba58fadead</anchor>
      <arglist>(const Coordinate2D &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>Matrix2D</type>
      <name>operator*</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a43efb6241566ff22802b8eba1359c304</anchor>
      <arglist>(const Matrix2D &amp;m)</arglist>
    </member>
    <member kind="friend">
      <type>friend Matrix2D</type>
      <name>operator*</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a73aa034cbe317c5ef757a4ec331b5814</anchor>
      <arglist>(Double_t a, const Matrix2D &amp;m)</arglist>
    </member>
    <member kind="friend">
      <type>friend std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>ac542e2c21aba305752338fa03a8985be</anchor>
      <arglist>(std::ostream &amp;stream, const Matrix2D &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>Matrix2D</type>
      <name>Transpose</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>aca866eeebcea903d982dbeebe69ac461</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>TransposeInPlace</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a3dc16a8e65878c113482340486a3fd77</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Det</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a6c8aeb365642710d259e232edee224c1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Matrix2D</type>
      <name>Adjugate</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>ac91c76081f8da74ca83cc9481b711961</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Matrix2D</type>
      <name>Inverse</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>a86c5d6cd4d53890f0642404797209581</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Coordinate2D</type>
      <name>Eigenvalues</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>aaa664bc1ad317b6aee9201f5fc39012b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Matrix2D</type>
      <name>Eigenvectors</name>
      <anchorfile>classMath_1_1Matrix2D.html</anchorfile>
      <anchor>aa172365cbf934d7446c9ec8b04e3a7d9</anchor>
      <arglist>(const Coordinate2D &amp;e) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>NBody::Node</name>
    <filename>classNBody_1_1Node.html</filename>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~Node</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>aac74fd996823e8faaf2c313721052e7b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Int_t</type>
      <name>GetID</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a93ba32779e74a1367bfb4174701470d2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Double_t</type>
      <name>GetBoundary</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a53a734804a729b4a1c51ddf843073680</anchor>
      <arglist>(int i, int j)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Int_t</type>
      <name>GetCount</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a87d5456efd138e49a2f91589cb032891</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Int_t</type>
      <name>GetStart</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a1a6605aa40e2e63aaa270d46365483fe</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Int_t</type>
      <name>GetEnd</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a7a36b9bc394f5edec7bb896c96eba905</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>SetID</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>ab934f75663ea0ae072d03b12deae10fa</anchor>
      <arglist>(Int_tree_t id)</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a3b3f8a41185425ea6cbbee5b9ea10972</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>af0953143e3d6d149acd89352afc798e7</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a68ebea3d82814c5c152508317a634f25</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>addf2b54c8e6a43e66cf2b37aa5197a3e</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>af684e0a70f7e36ad3db44f4aacb6d3de</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *v, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>ac149cbed37d85658be1e9ff1d87ffaba</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>af86f0ea0d18a0db3985df39999af7446</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>ad88e57f8b223e19d1eac927168a20584</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate v, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a50229944e9aa75caac8d0d569d06123b</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a9d1bd6dcc1375754a64df366c065495c</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a610096b92a96254b1f17d2fc23e5b2b6</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Particle &amp;p, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a8cb9a4ebac2799a514d0bda44f57ed70</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, Double_t *metric)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>ab44a2ca8fa32e1e92a1ad67dcff44456</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>abaefa85fc8a4bb858d3649321a04ea42</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v, Double_t *metric)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a82fc800b4055ce2d64641aa4cb75d65a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v, Double_t *metric)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a9c3588ec2550b41c09ebbb8248232596</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>ad1a0c9db66fd66f5c598e366b6ec7ed1</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a4d36091f5e3cf867df5ea1e3c8fdb4ea</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>afb4f2a21b24bf827049164446cfd8bc9</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a2c40cb88bc559d3128db8a88f8044ef4</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a8d050cec3f72007323a07f199c6c6f33</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a83c79ca820d385721b164bc90d17914b</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *v, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a6f9de6b00d8fd050e180a502eaadeab6</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a2229703148cef9b786daa99e3f1499d4</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a4ff68d0f47703273885a6d4bba37b00a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate v, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a15906f3e7a3ec08183cbbe1562acfaea</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a6e0dcc5aada1f7455d52bc13407dd617</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t *metric)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a1e148d287f6e282e6ba809e56e3b3964</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a4e847278c1a96c5be13f56d5a697cfbf</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t *metric)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a32fad4563643a0323ce30fcd6fc03bff</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t *metric)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a2a7bf63206140f1eee44d3ab44c3974f</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a700424fa526d9f10afe9eb1bdc9cd879</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestCriterionPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a238dd70e5199f5f7fa14af24b7f24b1d</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestCriterionPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a8373862881b003e13baf22d06afbd162</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a178abe5ddb955dfb157d5bcc4d3a5ce1</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>ad894d91b7b27f182ed7ea50e09ac2ffd</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *x, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>af30f19c1c7aef54a897461fd8546d303</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Coordinate x, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a2ef10cf0d530d2a4e6cec5bb512c06c1</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a4275e206231420c4bef0bece60b053a7</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Double_t *x, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>aabeccb7e8d9d2230b1ddfab7b17b33ee</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Coordinate x, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a41db2f1e82854ae15b71fd18d0d04724</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Int_t t, Int_t &amp;nt, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>acb02e52aef60689fb360de6831c8282e</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *x, Int_t &amp;nt, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a3f35abc29cb9be478c6793a53007dbdb</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Coordinate x, Int_t &amp;nt, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a3244409aa6d9f4f37c56432d98008c19</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Int_t t, Int_t &amp;nt, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a628b32b9328082d0c118401dcdd6e88e</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Double_t *x, Int_t &amp;nt, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>afde6a3ceb4612b3c5fed9e1ee1180a94</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Coordinate x, Int_t &amp;nt, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterion</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>abd6d1ea6021678e7f916077d61b98e39</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterion</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a9d16ec24ac7cc6035e140f78ae600ace</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Particle &amp;p, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>aeebf074510d533ef0ed658bef29cc764</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a231cfb7da0b692c5966504962e82bf2c</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>ad1ed855daa78356eacba7a678d72eed2</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a7a3fcffb5bbf25c3eb77e61a66d7262a</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Particle &amp;p, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionPeriodicTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a0d02d5e5cf82a727ba25a224908661bc</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Double_t *period, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionPeriodicTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a897d2d4a1f1ebb2d25d2675eb0796ad4</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionNoDist</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a190d7b10e61c19bd498e2ee2e81b7299</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionNoDist</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a49388b4aaecad98be12a0d1bf946156a</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Particle &amp;p, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionNoDistPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a76f6d03cd21d4b46fad029603c6ae21c</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Double_t *period, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionNoDistPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a58158e2855e9283bfcf23340bbd5e2f3</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FOFSearchBall</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>afdd286e8544503823b765245408f9dcb</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FOFSearchCriterion</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a19695e3b37a632c5340b4d90117fa699</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FOFSearchCriterionSetBasisForLinks</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a4a8ab5182dc5347ad4b5a57ed8faaa45</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FOFSearchBallPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>ae6007f657c60a02cc0f32b407f0881a0</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FOFSearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>addd673b25447f7e0dea2defbefddaea6</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FOFSearchCriterionSetBasisForLinksPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>af5832ef5c9b04aed4b52502faf550b82</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)=0</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>UInt_tree_t</type>
      <name>bucket_end</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a120b64926b4c5429d37801b10b9bf56c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>UInt_tree_t</type>
      <name>bucket_start</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a9893f82cb070c4b5418f191908ecd49b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>UInt_tree_t</type>
      <name>count</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a70025a991557f1e81ef0ba04c8350dbf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>UInt_tree_t</type>
      <name>nid</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>afc0ae669efdc9623a7396c18f0bef2f5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>unsigned short</type>
      <name>numdim</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a7314b3efd16e0ce777d65faa24d0811e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>DoublePos_t</type>
      <name>xbnd</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>abbd6f376fa1a9c7548efec9c44736c63</anchor>
      <arglist>[6][2]</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Int_t</type>
      <name>GetID</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a93ba32779e74a1367bfb4174701470d2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Double_t</type>
      <name>GetBoundary</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a53a734804a729b4a1c51ddf843073680</anchor>
      <arglist>(int i, int j)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Int_t</type>
      <name>GetCount</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a87d5456efd138e49a2f91589cb032891</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Int_t</type>
      <name>GetStart</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a1a6605aa40e2e63aaa270d46365483fe</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Int_t</type>
      <name>GetEnd</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a7a36b9bc394f5edec7bb896c96eba905</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>SetID</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>ab934f75663ea0ae072d03b12deae10fa</anchor>
      <arglist>(Int_tree_t id)</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a3b3f8a41185425ea6cbbee5b9ea10972</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>af0953143e3d6d149acd89352afc798e7</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a68ebea3d82814c5c152508317a634f25</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>addf2b54c8e6a43e66cf2b37aa5197a3e</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>af684e0a70f7e36ad3db44f4aacb6d3de</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *v, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>ac149cbed37d85658be1e9ff1d87ffaba</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>af86f0ea0d18a0db3985df39999af7446</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>ad88e57f8b223e19d1eac927168a20584</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate v, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a50229944e9aa75caac8d0d569d06123b</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a9d1bd6dcc1375754a64df366c065495c</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a610096b92a96254b1f17d2fc23e5b2b6</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Particle &amp;p, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a8cb9a4ebac2799a514d0bda44f57ed70</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, Double_t *metric)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>ab44a2ca8fa32e1e92a1ad67dcff44456</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>abaefa85fc8a4bb858d3649321a04ea42</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v, Double_t *metric)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a82fc800b4055ce2d64641aa4cb75d65a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v, Double_t *metric)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a9c3588ec2550b41c09ebbb8248232596</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>ad1a0c9db66fd66f5c598e366b6ec7ed1</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a4d36091f5e3cf867df5ea1e3c8fdb4ea</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>afb4f2a21b24bf827049164446cfd8bc9</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a2c40cb88bc559d3128db8a88f8044ef4</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a8d050cec3f72007323a07f199c6c6f33</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a83c79ca820d385721b164bc90d17914b</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *v, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a6f9de6b00d8fd050e180a502eaadeab6</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a2229703148cef9b786daa99e3f1499d4</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a4ff68d0f47703273885a6d4bba37b00a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate v, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a15906f3e7a3ec08183cbbe1562acfaea</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a6e0dcc5aada1f7455d52bc13407dd617</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t *metric)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a1e148d287f6e282e6ba809e56e3b3964</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a4e847278c1a96c5be13f56d5a697cfbf</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t *metric)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a32fad4563643a0323ce30fcd6fc03bff</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t *metric)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a2a7bf63206140f1eee44d3ab44c3974f</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a700424fa526d9f10afe9eb1bdc9cd879</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestCriterionPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a238dd70e5199f5f7fa14af24b7f24b1d</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FindNearestCriterionPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a8373862881b003e13baf22d06afbd162</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a178abe5ddb955dfb157d5bcc4d3a5ce1</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>ad894d91b7b27f182ed7ea50e09ac2ffd</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *x, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>af30f19c1c7aef54a897461fd8546d303</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Coordinate x, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a2ef10cf0d530d2a4e6cec5bb512c06c1</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a4275e206231420c4bef0bece60b053a7</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Double_t *x, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>aabeccb7e8d9d2230b1ddfab7b17b33ee</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Coordinate x, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a41db2f1e82854ae15b71fd18d0d04724</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Int_t t, Int_t &amp;nt, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>acb02e52aef60689fb360de6831c8282e</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *x, Int_t &amp;nt, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a3f35abc29cb9be478c6793a53007dbdb</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Coordinate x, Int_t &amp;nt, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a3244409aa6d9f4f37c56432d98008c19</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Int_t t, Int_t &amp;nt, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a628b32b9328082d0c118401dcdd6e88e</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Double_t *x, Int_t &amp;nt, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>afde6a3ceb4612b3c5fed9e1ee1180a94</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Coordinate x, Int_t &amp;nt, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterion</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>abd6d1ea6021678e7f916077d61b98e39</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterion</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a9d16ec24ac7cc6035e140f78ae600ace</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Particle &amp;p, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>aeebf074510d533ef0ed658bef29cc764</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a231cfb7da0b692c5966504962e82bf2c</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>ad1ed855daa78356eacba7a678d72eed2</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a7a3fcffb5bbf25c3eb77e61a66d7262a</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Particle &amp;p, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionPeriodicTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a0d02d5e5cf82a727ba25a224908661bc</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Double_t *period, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionPeriodicTagged</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a897d2d4a1f1ebb2d25d2675eb0796ad4</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionNoDist</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a190d7b10e61c19bd498e2ee2e81b7299</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionNoDist</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a49388b4aaecad98be12a0d1bf946156a</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Particle &amp;p, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionNoDistPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a76f6d03cd21d4b46fad029603c6ae21c</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Double_t *period, Int_t t, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>SearchCriterionNoDistPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a58158e2855e9283bfcf23340bbd5e2f3</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FOFSearchBall</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>afdd286e8544503823b765245408f9dcb</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FOFSearchCriterion</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a19695e3b37a632c5340b4d90117fa699</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FOFSearchCriterionSetBasisForLinks</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>a4a8ab5182dc5347ad4b5a57ed8faaa45</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FOFSearchBallPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>ae6007f657c60a02cc0f32b407f0881a0</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FOFSearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>addd673b25447f7e0dea2defbefddaea6</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>FOFSearchCriterionSetBasisForLinksPeriodic</name>
      <anchorfile>classNBody_1_1Node.html</anchorfile>
      <anchor>af5832ef5c9b04aed4b52502faf550b82</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>NBody::NPriorityQueue</name>
    <filename>classNBody_1_1NPriorityQueue.html</filename>
    <base>NBody::PriorityQueue</base>
    <class kind="struct">NBody::NPriorityQueue::nqueue</class>
    <member kind="function">
      <type></type>
      <name>NPriorityQueue</name>
      <anchorfile>classNBody_1_1NPriorityQueue.html</anchorfile>
      <anchor>aef6f282f004503c734c190e67717f116</anchor>
      <arglist>(Int_t max, Int_t Nval)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~NPriorityQueue</name>
      <anchorfile>classNBody_1_1NPriorityQueue.html</anchorfile>
      <anchor>a49a202d6673de9a81f624c53079cb983</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Pop</name>
      <anchorfile>classNBody_1_1NPriorityQueue.html</anchorfile>
      <anchor>a67634b6247887a277d36aa4f4b3295f6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Push</name>
      <anchorfile>classNBody_1_1NPriorityQueue.html</anchorfile>
      <anchor>ace38940e9845b63ff35f4c33e1ac5de7</anchor>
      <arglist>(Int_t p, Double_t dist, Double_t *pvalue)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>TopValue</name>
      <anchorfile>classNBody_1_1NPriorityQueue.html</anchorfile>
      <anchor>a688465abf95ddb1f1828656f133183c5</anchor>
      <arglist>(Int_t N)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>nqueue *</type>
      <name>npq</name>
      <anchorfile>classNBody_1_1NPriorityQueue.html</anchorfile>
      <anchor>a81977ffbe593369cfa764bfe919576d9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Int_t</type>
      <name>nvalues</name>
      <anchorfile>classNBody_1_1NPriorityQueue.html</anchorfile>
      <anchor>aaa1dc4dc26fe0b4d9450d844a8f6e0f4</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>NBody::NPriorityQueue::nqueue</name>
    <filename>structNBody_1_1NPriorityQueue_1_1nqueue.html</filename>
    <member kind="variable">
      <type>Double_t *</type>
      <name>values</name>
      <anchorfile>structNBody_1_1NPriorityQueue_1_1nqueue.html</anchorfile>
      <anchor>a4a102ac2e2ae8b55db1d0eb8a8b7be31</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>NBody::Particle</name>
    <filename>classNBody_1_1Particle.html</filename>
    <member kind="function">
      <type></type>
      <name>Particle</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a2e2241dee17e54647ed361b7e2f7a244</anchor>
      <arglist>(Double_t Mass=0, Double_t x=0, Double_t y=0, Double_t z=0, Double_t vx=0, Double_t vy=0, Double_t vz=0, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Int_t PID=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Particle</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a2c4a59ce622d2162e72a0c221153562b</anchor>
      <arglist>(Double_t Mass, Double_t *NewPos, Double_t *NewVel, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Int_t PID=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Particle</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a651fe60e65dc4a85ca53dc7a4fb67228</anchor>
      <arglist>(const Particle &amp;p)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Particle</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>aeebcd8320b77acf861e5eb36859839f0</anchor>
      <arglist>(std::istream &amp;F)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Particle</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>aa8d71f6b30785ed774a6e95673a1b269</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Write</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ae36a1ceec72d21ad960bc762ef48a539</anchor>
      <arglist>(std::ostream &amp;F) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetMass</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a657584618fe8e026192c371ae05dc68b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetMass</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>accef7a06c08dc808f2fbbc72a8e4e8e6</anchor>
      <arglist>(const Double_t &amp;Mass)</arglist>
    </member>
    <member kind="function">
      <type>const Double_t *</type>
      <name>GetPosition</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a8768c4f8727f73b66c8f441be728c017</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetPosition</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a6f9b21fdf5c2415a9174d46af9c09190</anchor>
      <arglist>(const int &amp;i) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>X</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>aaa516a974da9e23ac785489be26f54ee</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Y</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>acc1af1e687213aadf80a625da696e428</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Z</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a5a2efc5a2cb8bda85fe545d90bdabee2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetPosition</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a48169e163f6870a060582eb132679ccc</anchor>
      <arglist>(const int &amp;i, const Double_t &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetPosition</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a12b1f35526902001da14bb0268ecdf4d</anchor>
      <arglist>(const Double_t *p)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetPosition</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a212d76497af397f89b900b0334912331</anchor>
      <arglist>(const Double_t &amp;x, const Double_t &amp;y, const Double_t &amp;z)</arglist>
    </member>
    <member kind="function">
      <type>const Double_t *</type>
      <name>GetVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>af76b5400448f8cd0edc710c861cd2cbd</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a1f6983e8da7b0802c54b3da419a19aa5</anchor>
      <arglist>(const int &amp;i) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Vx</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>acbb926efe0a80f169ae035f5e9c5bdc0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Vy</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>aeab7a9f0729f84789abdbfb88ad20092</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Vz</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a116697ab59ba9bb30ebb095c624084ce</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a3caf8889214e89a70f98f00606ff3d4b</anchor>
      <arglist>(const int &amp;i, const Double_t &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ab2ad34cb3516802329372982697450b3</anchor>
      <arglist>(const Double_t *v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a8eb97f76ebf90ffb3395dd59dcc9187b</anchor>
      <arglist>(const Double_t &amp;vx, const Double_t &amp;vy, const Double_t &amp;vz)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetPhase</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>abe29ffbda2722393f2ea894d15fae341</anchor>
      <arglist>(const int &amp;i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetPhase</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>aa846d2f6d07d2f4b118aad6c7c305ac1</anchor>
      <arglist>(const int &amp;i, const Double_t &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetPID</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ab3df27c8f52d4c8e71542bceccda3d7f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetPID</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a413ecb4cf16d6fc075a60443b73987b2</anchor>
      <arglist>(const Int_t &amp;i)</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetID</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ac335d0985ee8ae1752a57d42ef820668</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetID</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a2661d341921f867d9a5a9d5d63c35bf7</anchor>
      <arglist>(Int_t i)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>GetType</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a6ac14552c23e880c846e6c83bf24f101</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetType</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>afe8c3063ef31f331c6c81a6ea69fdb33</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetDensity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ae95e0dd4f389b0ef6803aae9065778e4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetDensity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a98d88bd1818a5cfffe9ae1c246ac9413</anchor>
      <arglist>(const Double_t &amp;Rho)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetPotential</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a825dc7a3a482dd9eb5e0e7fcaa8f1315</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetPotential</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a75b7fb29664e10eadffe86cf1990ff12</anchor>
      <arglist>(const Double_t &amp;Phi)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ScalePositions</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a1d8a243499cb35b76158f3248cc02fbf</anchor>
      <arglist>(Double_t &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ScaleVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ad4b31026e2f6901baff915138a4e0391</anchor>
      <arglist>(Double_t &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ScalePhase</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a5bd7b27208c71402b07174214067c77c</anchor>
      <arglist>(Double_t &amp;x, Double_t &amp;v)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Radius</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a75c1fcdb3012f6da8d580ff3b98a0a12</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Radius2</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a433800f4ca81bd78e2a7cc78182bc612</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Theta</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a2344b6643ae3779bbd3f3d4c7eaedd9a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Phi</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>af9390bab9664e61d76c0885ba1b18c28</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CylRadius</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ae8dc2e028c07ae6376a2af2798610b3e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>AbsoluteVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a4c6444297c61ed5a263285ec5c37b1d2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>RadialVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a408d244678e51f02777ae48f9de8e1b7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CircularVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a610365ecd8ad5617f15eb98ef67e68d5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CylRadialVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ac7789140d42335d6760126d32ccaacaf</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CylTangentialVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a413cc99bf734ab02a3bfd913bad4e25e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>AngularMomentum</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a93e6df0b31d8e3e50d2533505d0889ce</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Double_t</type>
      <name>mass</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ae03bef61a91f31fd4ffceae86591429d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Double_t</type>
      <name>phi</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a54f3e17c3aae989646133f61aa1f6d73</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Double_t</type>
      <name>position</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ad28aa6a467675d8eea52a768dd5f184e</anchor>
      <arglist>[3]</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Double_t</type>
      <name>rho</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ae7b7ae9f068717e1c3f42035ae9ce67f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Double_t</type>
      <name>velocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>aa00c0856ac5b1832f2c0a591f8d8d9f6</anchor>
      <arglist>[3]</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Int_t</type>
      <name>pid</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a5dc0b36c678984bcd55269c17e5f4841</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Int_t</type>
      <name>id</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ac2faa7e0e9c37a21d9eec22925939650</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>type</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ad2a4619c08a9afe25b9d4a2a8d204a90</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Int_t</type>
      <name>pid</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a5dc0b36c678984bcd55269c17e5f4841</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Int_t</type>
      <name>id</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ac2faa7e0e9c37a21d9eec22925939650</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>type</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ad2a4619c08a9afe25b9d4a2a8d204a90</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Particle</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a2e2241dee17e54647ed361b7e2f7a244</anchor>
      <arglist>(Double_t Mass=0, Double_t x=0, Double_t y=0, Double_t z=0, Double_t vx=0, Double_t vy=0, Double_t vz=0, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Int_t PID=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Particle</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a2c4a59ce622d2162e72a0c221153562b</anchor>
      <arglist>(Double_t Mass, Double_t *NewPos, Double_t *NewVel, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Int_t PID=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Particle</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a651fe60e65dc4a85ca53dc7a4fb67228</anchor>
      <arglist>(const Particle &amp;p)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Particle</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>aeebcd8320b77acf861e5eb36859839f0</anchor>
      <arglist>(std::istream &amp;F)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Particle</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>aa8d71f6b30785ed774a6e95673a1b269</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Particle &amp;</type>
      <name>operator=</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a83e78a4ceb106644f8f32578eb55ceba</anchor>
      <arglist>(const Particle &amp;part)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator==</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a23708753246a0d3097a0767d9ce59837</anchor>
      <arglist>(const Particle &amp;p) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator!=</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ac5d5bc15d893b4352efdf609bbff8b61</anchor>
      <arglist>(const Particle &amp;p) const </arglist>
    </member>
    <member kind="friend">
      <type>friend ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a0c15391e5c766a616b4a39ca0c41421d</anchor>
      <arglist>(ostream &amp;outs, const Particle &amp;p)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Write</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ae36a1ceec72d21ad960bc762ef48a539</anchor>
      <arglist>(std::ostream &amp;F) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetMass</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a657584618fe8e026192c371ae05dc68b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetMass</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>accef7a06c08dc808f2fbbc72a8e4e8e6</anchor>
      <arglist>(const Double_t &amp;Mass)</arglist>
    </member>
    <member kind="function">
      <type>const Double_t *</type>
      <name>GetPosition</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a8768c4f8727f73b66c8f441be728c017</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetPosition</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a6f9b21fdf5c2415a9174d46af9c09190</anchor>
      <arglist>(const int &amp;i) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>X</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>aaa516a974da9e23ac785489be26f54ee</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Y</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>acc1af1e687213aadf80a625da696e428</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Z</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a5a2efc5a2cb8bda85fe545d90bdabee2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetPosition</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a48169e163f6870a060582eb132679ccc</anchor>
      <arglist>(const int &amp;i, const Double_t &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetPosition</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a12b1f35526902001da14bb0268ecdf4d</anchor>
      <arglist>(const Double_t *p)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetPosition</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a212d76497af397f89b900b0334912331</anchor>
      <arglist>(const Double_t &amp;x, const Double_t &amp;y, const Double_t &amp;z)</arglist>
    </member>
    <member kind="function">
      <type>const Double_t *</type>
      <name>GetVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>af76b5400448f8cd0edc710c861cd2cbd</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a1f6983e8da7b0802c54b3da419a19aa5</anchor>
      <arglist>(const int &amp;i) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Vx</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>acbb926efe0a80f169ae035f5e9c5bdc0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Vy</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>aeab7a9f0729f84789abdbfb88ad20092</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Vz</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a116697ab59ba9bb30ebb095c624084ce</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a3caf8889214e89a70f98f00606ff3d4b</anchor>
      <arglist>(const int &amp;i, const Double_t &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ab2ad34cb3516802329372982697450b3</anchor>
      <arglist>(const Double_t *v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a8eb97f76ebf90ffb3395dd59dcc9187b</anchor>
      <arglist>(const Double_t &amp;vx, const Double_t &amp;vy, const Double_t &amp;vz)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetPhase</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>abe29ffbda2722393f2ea894d15fae341</anchor>
      <arglist>(const int &amp;i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetPhase</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>aa846d2f6d07d2f4b118aad6c7c305ac1</anchor>
      <arglist>(const int &amp;i, const Double_t &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetPID</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ab3df27c8f52d4c8e71542bceccda3d7f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetPID</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a413ecb4cf16d6fc075a60443b73987b2</anchor>
      <arglist>(const Int_t &amp;i)</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetID</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ac335d0985ee8ae1752a57d42ef820668</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetID</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a2661d341921f867d9a5a9d5d63c35bf7</anchor>
      <arglist>(Int_t i)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>GetType</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a6ac14552c23e880c846e6c83bf24f101</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetType</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>afe8c3063ef31f331c6c81a6ea69fdb33</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetDensity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ae95e0dd4f389b0ef6803aae9065778e4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetDensity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a98d88bd1818a5cfffe9ae1c246ac9413</anchor>
      <arglist>(const Double_t &amp;Rho)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetPotential</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a825dc7a3a482dd9eb5e0e7fcaa8f1315</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetPotential</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a75b7fb29664e10eadffe86cf1990ff12</anchor>
      <arglist>(const Double_t &amp;Phi)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ScalePositions</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a1d8a243499cb35b76158f3248cc02fbf</anchor>
      <arglist>(Double_t &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ScaleVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ad4b31026e2f6901baff915138a4e0391</anchor>
      <arglist>(Double_t &amp;x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ScalePhase</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a5bd7b27208c71402b07174214067c77c</anchor>
      <arglist>(Double_t &amp;x, Double_t &amp;v)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Radius</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a75c1fcdb3012f6da8d580ff3b98a0a12</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Radius2</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a433800f4ca81bd78e2a7cc78182bc612</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Theta</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a2344b6643ae3779bbd3f3d4c7eaedd9a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Phi</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>af9390bab9664e61d76c0885ba1b18c28</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CylRadius</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ae8dc2e028c07ae6376a2af2798610b3e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>AbsoluteVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a4c6444297c61ed5a263285ec5c37b1d2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>RadialVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a408d244678e51f02777ae48f9de8e1b7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CircularVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a610365ecd8ad5617f15eb98ef67e68d5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CylRadialVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>ac7789140d42335d6760126d32ccaacaf</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CylTangentialVelocity</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a413cc99bf734ab02a3bfd913bad4e25e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>AngularMomentum</name>
      <anchorfile>classNBody_1_1Particle.html</anchorfile>
      <anchor>a93e6df0b31d8e3e50d2533505d0889ce</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>NBody::PriorityQueue</name>
    <filename>classNBody_1_1PriorityQueue.html</filename>
    <class kind="struct">NBody::PriorityQueue::queue_member</class>
    <member kind="function">
      <type></type>
      <name>PriorityQueue</name>
      <anchorfile>classNBody_1_1PriorityQueue.html</anchorfile>
      <anchor>a8856eaaa82abf504121f1d217c84a0d9</anchor>
      <arglist>(Int_t max)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~PriorityQueue</name>
      <anchorfile>classNBody_1_1PriorityQueue.html</anchorfile>
      <anchor>aaeb5cd9028e94508ed093817698931d1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>Empty</name>
      <anchorfile>classNBody_1_1PriorityQueue.html</anchorfile>
      <anchor>a573564f0b2192444887684aee4827194</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>MaxSize</name>
      <anchorfile>classNBody_1_1PriorityQueue.html</anchorfile>
      <anchor>a0830b58422db5d7eddffaa32cf9751a6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Pop</name>
      <anchorfile>classNBody_1_1PriorityQueue.html</anchorfile>
      <anchor>a01515732034596da0ad9a257ad78654e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Push</name>
      <anchorfile>classNBody_1_1PriorityQueue.html</anchorfile>
      <anchor>a8390a0706a909093bcd698af435d0220</anchor>
      <arglist>(Int_t p, Double_t dist)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Reset</name>
      <anchorfile>classNBody_1_1PriorityQueue.html</anchorfile>
      <anchor>a8ecd812d7df547da1c0857b1aa959257</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>Size</name>
      <anchorfile>classNBody_1_1PriorityQueue.html</anchorfile>
      <anchor>a519a2d35c3efa8ea14e7d36337fcc507</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>TopPriority</name>
      <anchorfile>classNBody_1_1PriorityQueue.html</anchorfile>
      <anchor>ae0c95361a5cc073fafff1229f1fa666f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>TopQueue</name>
      <anchorfile>classNBody_1_1PriorityQueue.html</anchorfile>
      <anchor>aa32533669227c5fe7e899ed4812614ec</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Int_t</type>
      <name>max_size</name>
      <anchorfile>classNBody_1_1PriorityQueue.html</anchorfile>
      <anchor>abfae47eb1b70ae69f90d2654d5a083b2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Int_t</type>
      <name>n</name>
      <anchorfile>classNBody_1_1PriorityQueue.html</anchorfile>
      <anchor>af8d16add6021faf5088aea57b9e5d0c3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>queue_member *</type>
      <name>pq</name>
      <anchorfile>classNBody_1_1PriorityQueue.html</anchorfile>
      <anchor>ac2d56450fe2b2f8a3e7a57fdaa2a4f22</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>NBody::PriorityQueue::queue_member</name>
    <filename>structNBody_1_1PriorityQueue_1_1queue__member.html</filename>
    <member kind="variable">
      <type>Double_t</type>
      <name>priority</name>
      <anchorfile>structNBody_1_1PriorityQueue_1_1queue__member.html</anchorfile>
      <anchor>a8dd19c7385520baeb3488484c14ae550</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Int_t</type>
      <name>queue</name>
      <anchorfile>structNBody_1_1PriorityQueue_1_1queue__member.html</anchorfile>
      <anchor>ac0b6550c39e490b38d92a8ca142385ae</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>NBody::smoothfunc</name>
    <filename>structNBody_1_1smoothfunc.html</filename>
    <member kind="variable">
      <type>Double_t(*</type>
      <name>function</name>
      <anchorfile>structNBody_1_1smoothfunc.html</anchorfile>
      <anchor>a251ea72c5dcd7a7a94b2e4976d900ed2</anchor>
      <arglist>)(Double_t r, Double_t h)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>NBody::SplitNode</name>
    <filename>classNBody_1_1SplitNode.html</filename>
    <base>NBody::Node</base>
    <member kind="function">
      <type></type>
      <name>SplitNode</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a6d7a574c6280b495ab4fe4c543e2121e</anchor>
      <arglist>(Int_t id, int d, Double_t p, Int_t Count, Double_t bnd[6][2], Int_t new_bucket_start, Int_t new_bucket_end, unsigned short ndim, Node *initial_left=NULL, Node *initial_right=NULL)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~SplitNode</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ac914b6a234a3c6dc72b3aaa58aa14e9b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetCutValue</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a37c7e975f949666e1c08a55f17c9263a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>GetCutDim</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a4a3241a630c97e0f57efcd3f67f661d0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Node *</type>
      <name>GetLeft</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>aa007e2b7b0f77dd51314e1761fa09ef2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Node *</type>
      <name>GetRight</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a307d0b68eed27acd1438a6e80ddff897</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ab6f90612d701225b456f79afefe6121b</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a11275405f27f0067bf66d99ee94fd52d</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>adf164b244a596a0d3c3668feb59b0252</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>abb2ccaadd88467f734dae873c548145a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a81f08c658f871ed7e240b36734afd715</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *v, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a2967790c1aba7be0f4041a118d2009f1</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ae297ba85f233ab0b44da06e125992f74</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ad2f771a2ee71f20bb58a33631f512d45</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate v, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a9b4dbdcb09c201e01332a07897a5a279</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a20e794f26ba35446820750bf78af92e9</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ad2d1e95f1b8a6a397c1a72d701bef8bf</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a1c9341b3830046eb78da8d08b680568a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>abefbc97d52c5ca1ca094958dceb57dff</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ad80c1e4e80c4e87be8c2ed7e342589a3</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a6e09856e08c12d5109c40b649952da60</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a32fa13b48dc0b9c6e606e3de1d9498e0</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a65d1b82b54298690b1a4e33bbe318a73</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a53c1c36c0365a58b58d36d5024012f3c</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a7db0da5352da7d18b2a7c8ae7c2575da</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ad57658f349f6b634a07bc5c57941f212</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Coordinate x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a6ea7aed967225d435e7ac375925a08e0</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Int_t t, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a020469782c0e662d0564e8feb4010865</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *x, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a09f1c094e97a6c226be46b18a9b81c74</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Coordinate x, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterion</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ad3b059581af34323bc6ebfce4766f793</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterion</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>acce66a323005756c843a79ef025dd42a</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a5120ff0a3d306f591a98d9172fb450c0</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a2c74db34a53837c8b85325540eb5d8d2</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionNoDist</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a89384fc2352a9a72fb1057be4f0c51da</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionNoDist</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>af83e71b78a101377bc8d3afa1600c264</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchBall</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a8b1ec3920e4fd3b7db8ff464d50d464d</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchCriterion</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ae3f459b7581a67588302c344c86f650b</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchCriterionSetBasisForLinks</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>aa6c623f205e1267f39f1326a3f574e7e</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a53512c6a8bbf5d96cfb70419f51c3ed0</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a71bf42965d0e7d24419eeedaa160c868</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a21b0e2faae6290016a53be2af6ad404d</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a2a290853d6aa01a5c3e2e5e15d8b42bf</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>af50803122d2dbfb2c2e65e9eaa3f161e</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *v, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>aa66d4ad60af9155eb4ac0137d7bbef1b</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a73c05c76e58a093dcfb00dc7f3fd60c1</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a5a10e935a80ee9d50c503e5d26fe528a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate v, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a8fe3b3b4d2bc00dc6504339afbf9e743</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a6ae89d5497f0265bbc29fd5081964266</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a8aa3497acf26832b4620006f1d4e41af</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a5790ecdbd770be49c9d0c0ba6f592e4d</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a70bdaf06733b25f3a10aff36ef5ffa2f</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a1a08a5952a7a4a82d5906650090955dc</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>aafecb5d7f97e72c1cd18059015fef39a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterionPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>aad29939db5d0fc83bb4fdfa539382959</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterionPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>af02ce26ba68429ca81cc03bd68d17b14</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a1530f356826dadec85edafb894681dfe</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a5469c6d447a534793ffa56056f4d60ed</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Double_t *x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>aa3352b16f2a651c02cd505679dbf44ff</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Coordinate x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a150eef303f87f5185ddae338f3f4aad9</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Int_t t, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a3fcf57378621912dc7d85a9f53035643</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Double_t *x, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a4ca5b81d669b9f43696e01aebcfb7f75</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Coordinate x, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>acf71c612b57dd9f035289eba13e8c608</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a7fd729a134235babbaa0288553341c7f</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionPeriodicTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>abca9e34f04a91c74811918b487cb0098</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionPeriodicTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a5229cf2013b5bced2bfe83694ac08f40</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionNoDistPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>af087f7bb261aeb4d1328f3ddfe05df57</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionNoDistPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>adf12057eb50621e87aad1cbd41e5cf09</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchBallPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a28bc4b14c07286386f0c4c1e42f3fb50</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ac833da83ddca0ab4febf30eca267b1a2</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchCriterionSetBasisForLinksPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a72924c518bc16c8089af676e9057c0c8</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>cut_dim</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a2c3871ea77fd583920a92d88dd9b9486</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>cut_val</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ae264c61eb2a46bb746900d30693ae360</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Node *</type>
      <name>left</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a0b757d9a7f3b637c8d986cb48325d257</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Node *</type>
      <name>right</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a4723c5f4c903fbb04c1d44f9e361d74e</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetCutValue</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a37c7e975f949666e1c08a55f17c9263a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>GetCutDim</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a4a3241a630c97e0f57efcd3f67f661d0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Node *</type>
      <name>GetLeft</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>aa007e2b7b0f77dd51314e1761fa09ef2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Node *</type>
      <name>GetRight</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a307d0b68eed27acd1438a6e80ddff897</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ab6f90612d701225b456f79afefe6121b</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a11275405f27f0067bf66d99ee94fd52d</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>adf164b244a596a0d3c3668feb59b0252</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>abb2ccaadd88467f734dae873c548145a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a81f08c658f871ed7e240b36734afd715</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *v, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a2967790c1aba7be0f4041a118d2009f1</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPos</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ae297ba85f233ab0b44da06e125992f74</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVel</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ad2f771a2ee71f20bb58a33631f512d45</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate v, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhase</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a9b4dbdcb09c201e01332a07897a5a279</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a20e794f26ba35446820750bf78af92e9</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ad2d1e95f1b8a6a397c1a72d701bef8bf</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a1c9341b3830046eb78da8d08b680568a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetric</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>abefbc97d52c5ca1ca094958dceb57dff</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ad80c1e4e80c4e87be8c2ed7e342589a3</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensor</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a6e09856e08c12d5109c40b649952da60</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a32fa13b48dc0b9c6e606e3de1d9498e0</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterion</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a65d1b82b54298690b1a4e33bbe318a73</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a53c1c36c0365a58b58d36d5024012f3c</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a7db0da5352da7d18b2a7c8ae7c2575da</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPos</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ad57658f349f6b634a07bc5c57941f212</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Coordinate x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a6ea7aed967225d435e7ac375925a08e0</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Int_t t, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a020469782c0e662d0564e8feb4010865</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *x, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a09f1c094e97a6c226be46b18a9b81c74</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Coordinate x, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterion</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ad3b059581af34323bc6ebfce4766f793</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterion</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>acce66a323005756c843a79ef025dd42a</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a5120ff0a3d306f591a98d9172fb450c0</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a2c74db34a53837c8b85325540eb5d8d2</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionNoDist</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a89384fc2352a9a72fb1057be4f0c51da</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionNoDist</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>af83e71b78a101377bc8d3afa1600c264</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchBall</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a8b1ec3920e4fd3b7db8ff464d50d464d</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchCriterion</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ae3f459b7581a67588302c344c86f650b</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchCriterionSetBasisForLinks</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>aa6c623f205e1267f39f1326a3f574e7e</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a53512c6a8bbf5d96cfb70419f51c3ed0</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a71bf42965d0e7d24419eeedaa160c868</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a21b0e2faae6290016a53be2af6ad404d</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a2a290853d6aa01a5c3e2e5e15d8b42bf</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>af50803122d2dbfb2c2e65e9eaa3f161e</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *v, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>aa66d4ad60af9155eb4ac0137d7bbef1b</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPosPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a73c05c76e58a093dcfb00dc7f3fd60c1</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestVelPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a5a10e935a80ee9d50c503e5d26fe528a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate v, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestPhasePeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a8fe3b3b4d2bc00dc6504339afbf9e743</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a6ae89d5497f0265bbc29fd5081964266</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a8aa3497acf26832b4620006f1d4e41af</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a5790ecdbd770be49c9d0c0ba6f592e4d</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a70bdaf06733b25f3a10aff36ef5ffa2f</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a1a08a5952a7a4a82d5906650090955dc</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Double_t *x, Double_t *v, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestMetricwithTensorPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>aafecb5d7f97e72c1cd18059015fef39a</anchor>
      <arglist>(Double_t rd, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Coordinate x, Coordinate v, Double_t *m0, Double_t *m1, GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterionPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>aad29939db5d0fc83bb4fdfa539382959</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FindNearestCriterionPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>af02ce26ba68429ca81cc03bd68d17b14</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, PriorityQueue *pq, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a1530f356826dadec85edafb894681dfe</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a5469c6d447a534793ffa56056f4d60ed</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Double_t *x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>aa3352b16f2a651c02cd505679dbf44ff</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Coordinate x, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a150eef303f87f5185ddae338f3f4aad9</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Int_t t, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a3fcf57378621912dc7d85a9f53035643</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Double_t *x, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchBallPosPeriodicTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a4ca5b81d669b9f43696e01aebcfb7f75</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Particle *bucket, Int_t *tagged, Double_t *off, Double_t *period, Coordinate x, Int_t &amp;nt, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>acf71c612b57dd9f035289eba13e8c608</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a7fd729a134235babbaa0288553341c7f</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *pdist2, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionPeriodicTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>abca9e34f04a91c74811918b487cb0098</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionPeriodicTagged</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a5229cf2013b5bced2bfe83694ac08f40</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Particle *bucket, Int_t &amp;nt, Int_t *tagged, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionNoDistPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>af087f7bb261aeb4d1328f3ddfe05df57</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Double_t *period, Int_t t, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SearchCriterionNoDistPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>adf12057eb50621e87aad1cbd41e5cf09</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Particle *bucket, Int_t *Group, Double_t *off, Double_t *period, Particle &amp;p, int dim=3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchBallPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a28bc4b14c07286386f0c4c1e42f3fb50</anchor>
      <arglist>(Double_t rd, Double_t fdist2, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchCriterionPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>ac833da83ddca0ab4febf30eca267b1a2</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>FOFSearchCriterionSetBasisForLinksPeriodic</name>
      <anchorfile>classNBody_1_1SplitNode.html</anchorfile>
      <anchor>a72924c518bc16c8089af676e9057c0c8</anchor>
      <arglist>(Double_t rd, FOFcompfunc cmp, FOFcheckfunc check, Double_t *params, Int_t iGroup, Int_t nActive, Particle *bucket, Int_t *Group, Int_tree_t *Len, Int_tree_t *Head, Int_tree_t *Tail, Int_tree_t *Next, short *BucketFlag, Int_tree_t *Fifo, Int_t &amp;iTail, Double_t *off, Double_t *period, Int_t target)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>NBody::StarParticle</name>
    <filename>classNBody_1_1StarParticle.html</filename>
    <base>NBody::Particle</base>
    <member kind="function">
      <type></type>
      <name>StarParticle</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>a32795333ac4ab28efd00bb5b0d0f3b89</anchor>
      <arglist>(Double_t Mass=0, Double_t x=0, Double_t y=0, Double_t z=0, Double_t vx=0, Double_t vy=0, Double_t vz=0, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Double_t TF=0, Double_t Zi=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>StarParticle</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>a93bf28c7f7277879d21b267dc7e41233</anchor>
      <arglist>(Double_t Mass, Double_t *NewPos, Double_t *NewVel, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Double_t TF=0, Double_t Zi=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>StarParticle</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>af65f5eea33a2f1bb1cc00d2b6863349d</anchor>
      <arglist>(const StarParticle &amp;p)</arglist>
    </member>
    <member kind="function">
      <type>StarParticle &amp;</type>
      <name>operator=</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>af9d5c8ce19a810ba85e3223d2fda994a</anchor>
      <arglist>(const StarParticle &amp;part)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator==</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>a2ea657855db5b5cd8bc46edb95e65ce2</anchor>
      <arglist>(const StarParticle &amp;p)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator!=</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>ad65fe3f5b32ae6d69a6ea25128c5aa91</anchor>
      <arglist>(const StarParticle &amp;p) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetFormationTime</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>a5507047bad91660957fba6a475597999</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetFormationTime</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>a1fe47a28820438549d2e1210da4c029b</anchor>
      <arglist>(Double_t Tform)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetZ</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>aebc025bb46bd6b5c51f949dcd9a56f20</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetZ</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>a0a0292dfae87e6682acddd9ebfb9651b</anchor>
      <arglist>(Double_t zz)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>metal</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>a435078f7fc0c732602856d6e380e8aa2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>Double_t</type>
      <name>tform</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>ad651422055c499f1e0dda3308ffd4829</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>StarParticle</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>a32795333ac4ab28efd00bb5b0d0f3b89</anchor>
      <arglist>(Double_t Mass=0, Double_t x=0, Double_t y=0, Double_t z=0, Double_t vx=0, Double_t vy=0, Double_t vz=0, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Double_t TF=0, Double_t Zi=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>StarParticle</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>a93bf28c7f7277879d21b267dc7e41233</anchor>
      <arglist>(Double_t Mass, Double_t *NewPos, Double_t *NewVel, Int_t ID=0, int type=0, Double_t Rho=0, Double_t Phi=0, Double_t TF=0, Double_t Zi=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>StarParticle</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>af65f5eea33a2f1bb1cc00d2b6863349d</anchor>
      <arglist>(const StarParticle &amp;p)</arglist>
    </member>
    <member kind="function">
      <type>StarParticle &amp;</type>
      <name>operator=</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>af9d5c8ce19a810ba85e3223d2fda994a</anchor>
      <arglist>(const StarParticle &amp;part)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator==</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>a2ea657855db5b5cd8bc46edb95e65ce2</anchor>
      <arglist>(const StarParticle &amp;p)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator!=</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>ad65fe3f5b32ae6d69a6ea25128c5aa91</anchor>
      <arglist>(const StarParticle &amp;p) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetFormationTime</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>a5507047bad91660957fba6a475597999</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetFormationTime</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>a1fe47a28820438549d2e1210da4c029b</anchor>
      <arglist>(Double_t Tform)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetZ</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>aebc025bb46bd6b5c51f949dcd9a56f20</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetZ</name>
      <anchorfile>classNBody_1_1StarParticle.html</anchorfile>
      <anchor>a0a0292dfae87e6682acddd9ebfb9651b</anchor>
      <arglist>(Double_t zz)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>NBody::System</name>
    <filename>classNBody_1_1System.html</filename>
    <member kind="function">
      <type></type>
      <name>System</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a830fbb3709ffc4c12f2e7451ce35270e</anchor>
      <arglist>(Int_t num=1, Double_t t=0.0, Double_t *period=NULL, Int_t ng=0, Int_t ns=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>System</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ad9b58a37d8bf18ff6c616fb4fbc457c9</anchor>
      <arglist>(Int_t num, Particle *p, Double_t t=0.0, Double_t *period=NULL, Int_t ng=0, Int_t ns=0, GasParticle *gp=NULL, StarParticle *sp=NULL)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>System</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a49da90019817c63a129fb1d2a9001712</anchor>
      <arglist>(const System &amp;s)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~System</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ac0b59d7375bf7b8cec99945b681feaf7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>TotalMass</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ad1d8e856a8bab272c52582ef2bcda3dd</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MaxLength</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aff7f3ff197367b808c678ffab23daadc</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MaxRadius</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a2514c5516b4b1ec84a0593d2022566e9</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MinRadius</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>af62bccbd62ee66c1b3145cd291593f84</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t *</type>
      <name>RadialLimits</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aad5a851bb2efeb92e9d0e38a4dbd1a54</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>AverageRadius</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a3cd8e8c921ef68b36deb16ca06226c89</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>AverageSpeed</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a6554e41dd8e585bc53895fb60df27211</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>AverageVelocity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a02c9a0c987fed87a5e29e905038cc33c</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>AverageRadialVelocity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ab693fb1cde2c7f235a13f4aee1124661</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MaxRadialVelocity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>acc5977cdca9cc6f6322ffafefd9f4806</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MaxCircularVelocity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a2c2cbb2049d868519a90a7fede858153</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>AverageDensity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a8c552c851ae7c783a6846886b6396dd3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>RegionDensity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a76bfb306c6026f2a304139cfa58d93c2</anchor>
      <arglist>(Double_t R, Double_t x, Double_t y, Double_t z) const </arglist>
    </member>
    <member kind="function">
      <type>Coord</type>
      <name>AngularMomentum</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>af4746eda2310a53b8ab79dce3b733c62</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Coord</type>
      <name>AngularMomentum</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aa776e4f1b922627bb78a3d89306c4522</anchor>
      <arglist>(Coord) const </arglist>
    </member>
    <member kind="function">
      <type>Coord</type>
      <name>AngularMomentum</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aa1ac1ee5a90ee817354748437b7ba0d9</anchor>
      <arglist>(Coordinate) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SortByID</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a805a3da6c40d1f660b5f833222013165</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SortByRadius</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aed5cb393db9538e850dd7db9c2cdeeaa</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SortByDistance</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a8cd5cfdd64bc6c3a375c0a069732333c</anchor>
      <arglist>(Coordinate c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SortByDistance</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a514888d05e2edca0f5f936e68a52ec1f</anchor>
      <arglist>(Coordinate c, int start, int n)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SortByEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ae26e69b8ea7491e07ec8be2987016ac5</anchor>
      <arglist>(Double_t G, Double_t eps)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SortByEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a4bfad3b86f150c491d9279afccb8a5b9</anchor>
      <arglist>(Double_t G, Double_t eps, int start, int n)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SortByDensity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a9b74a089ecd134f9846f232884e8d9fa</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SortByDensity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a39441652159ce43a4902045144c2148b</anchor>
      <arglist>(int start, int n)</arglist>
    </member>
    <member kind="function">
      <type>Coord</type>
      <name>CM</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a45573b5f5529d4b5d816349499d4c5c0</anchor>
      <arglist>(Double_t tolerance=0.0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AdjustForCM</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a09377e523306f330b95c778ae63e48e1</anchor>
      <arglist>(Double_t tolerance=0.0)</arglist>
    </member>
    <member kind="function">
      <type>Coord</type>
      <name>CMVel</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>adabe1682b68df10112ea07c9a136dd1a</anchor>
      <arglist>(Double_t tolerance=0.0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AdjustForCMVel</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aa00818e55f49911c81d6dec0958a243d</anchor>
      <arglist>(Double_t tolerance=0.0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AdjustPosition</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a6d82a74ac4067b7ea40a8861e10ee4d6</anchor>
      <arglist>(Coord)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AdjustPosition</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a3f0677eca5a6a9e17b58ba55fca31fcc</anchor>
      <arglist>(Coordinate)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AdjustVelocity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ab955207d46aa650ae2db1bac6bbbcb66</anchor>
      <arglist>(Coord)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AdjustVelocity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aaa5c92b7e95afe8b70203960cb4278d1</anchor>
      <arglist>(Coordinate)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>KineticEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a3c3ca4fe980687aaf850af7059458c2f</anchor>
      <arglist>(Int_t i) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>KineticEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a34e0418fca052b017ce1ba783774d835</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>KineticEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ae1c893bd9170ea7e2ed2da5fe3d74595</anchor>
      <arglist>(Coord, Int_t) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>KineticEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a6a464c1fa0b538898249cbed595441a4</anchor>
      <arglist>(Coord) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PotentialEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ab3763e2839f276c8d6ee7c2f53e23649</anchor>
      <arglist>(Int_t i, Int_t j, Double_t eps=0.0) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PotentialEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a0cb4e8f58c8c5f8d8db21df78ca2b37e</anchor>
      <arglist>(Int_t i, Double_t eps=0.0) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PotentialEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a7d73932eec5e6f5a642cbc1ccb7e3a04</anchor>
      <arglist>(Double_t eps=0.0, bool sphericalapprox=false) const </arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>FindParticle</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ace301b05392e7ccab56320374c38f6d8</anchor>
      <arglist>(const Particle &amp;p) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AddParticle</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ad2af64c8008adbdbc6e64d59e5547116</anchor>
      <arglist>(Particle part)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>RemoveParticle</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>addb46f1aa30205fe69d8cc4b0a5da4a7</anchor>
      <arglist>(Particle part)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>RemoveParticle</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aca3fa8a69ddb7bf624ba28faaad4acd5</anchor>
      <arglist>(Int_t i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AddSystem</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a4961f58090fa4cb88b5be5015f283aa5</anchor>
      <arglist>(const System &amp;s, const Double_t *o)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AddSystem</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a202c2ddf1bdd5f125098ba0f89847caa</anchor>
      <arglist>(const System &amp;s, const Double_t *o, const Double_t *vo)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ExtractSphere</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>af4756e63b898bb1ec9404e2dcdff0c31</anchor>
      <arglist>(Double_t radius)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Write</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a3dc72719dfc07b206b34c652d61b8ec5</anchor>
      <arglist>(ofstream &amp;F) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Write</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ab80d43d340bb5c5d017bf2dd5d7ab6c8</anchor>
      <arglist>(char *Filename) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Write</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a3ae136ff83cf04dd1b950a9c81f1b926</anchor>
      <arglist>(FILE *stream=stdout) const </arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Coordinate</type>
      <name>period</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a32b1b6f5e7be5dd9735f9c0a19186758</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Double_t</type>
      <name>time</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a379dd3983f99b07a3ba72bdf75e6bd89</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Int_t</type>
      <name>numparts</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ae54672e38c5e6ddaea89e863f65f1819</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Int_t</type>
      <name>ndark</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a26541081cffc28070f63f4c5fd62c2e4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Int_t</type>
      <name>ngas</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a05b5dbb672212d6a7349de6a3485087e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Int_t</type>
      <name>nstar</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a87d5dccb4e8047597c0b162a7f0dc7d0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Particle *</type>
      <name>particle</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a6e17976dad6cba33e13887fa67c3c899</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>GasParticle *</type>
      <name>gparticle</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a71d48bd39df2213b546693aef8921f3f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>StarParticle *</type>
      <name>sparticle</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a81aec42113611b830002aa193d90def1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Int_t</type>
      <name>numparts</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ae54672e38c5e6ddaea89e863f65f1819</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Int_t</type>
      <name>ndark</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a26541081cffc28070f63f4c5fd62c2e4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Int_t</type>
      <name>ngas</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a05b5dbb672212d6a7349de6a3485087e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Int_t</type>
      <name>nstar</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a87d5dccb4e8047597c0b162a7f0dc7d0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>Particle *</type>
      <name>particle</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a6e17976dad6cba33e13887fa67c3c899</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>GasParticle *</type>
      <name>gparticle</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a71d48bd39df2213b546693aef8921f3f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>StarParticle *</type>
      <name>sparticle</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a81aec42113611b830002aa193d90def1</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>System</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a830fbb3709ffc4c12f2e7451ce35270e</anchor>
      <arglist>(Int_t num=1, Double_t t=0.0, Double_t *period=NULL, Int_t ng=0, Int_t ns=0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>System</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ad9b58a37d8bf18ff6c616fb4fbc457c9</anchor>
      <arglist>(Int_t num, Particle *p, Double_t t=0.0, Double_t *period=NULL, Int_t ng=0, Int_t ns=0, GasParticle *gp=NULL, StarParticle *sp=NULL)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>System</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a49da90019817c63a129fb1d2a9001712</anchor>
      <arglist>(const System &amp;s)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~System</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ac0b59d7375bf7b8cec99945b681feaf7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetNumParts</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a3c1123dc2b53b41a5e6abacaab5aafa1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetNumDark</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>af426636c7358285cb8af011a4e49c1de</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetNumGas</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a77987c5c9f2c3f0e274d7fe1441cec47</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>GetNumStar</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aa7f9926558757c3bf32105b816d0af61</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GetTime</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>add247263b798f103f2a9aa80113f2d8c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetTime</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a67226178e17105331bf003338833832a</anchor>
      <arglist>(Double_t t)</arglist>
    </member>
    <member kind="function">
      <type>Particle *</type>
      <name>Parts</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ae499b811f9c351521938baa795806561</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Particle &amp;</type>
      <name>Part</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a922fc38fb6f2b212ba0fcc54cccdfe7d</anchor>
      <arglist>(Int_t i) const </arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>GetPeriod</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a2eda6bd46e088801b7cb822492c0851b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetPeriod</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a220aa7e59f8a2341961507954b91ee6d</anchor>
      <arglist>(Coordinate p)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetPeriod</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aeab17d146d5dfdb94ec1f5e105e0da53</anchor>
      <arglist>(Coord p)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SetPeriod</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a636c365769d5128668c25d939f2034cb</anchor>
      <arglist>(Double_t *p)</arglist>
    </member>
    <member kind="function">
      <type>System &amp;</type>
      <name>operator=</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ad19f7bf478af04fa5b7654956bd3b0d1</anchor>
      <arglist>(const System &amp;)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator==</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a82c9e957707e728e89f435e56e1c6f5b</anchor>
      <arglist>(const System &amp;) const </arglist>
    </member>
    <member kind="function">
      <type>Particle &amp;</type>
      <name>operator[]</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a7f7bd6d34a5c6207bdb2f71378db860b</anchor>
      <arglist>(Int_t i)</arglist>
    </member>
    <member kind="function">
      <type>Particle</type>
      <name>operator[]</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a7e8d1c8d59a2450c2b2596ca3cafd6f9</anchor>
      <arglist>(Int_t i) const </arglist>
    </member>
    <member kind="friend">
      <type>friend ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aecf8342d3ff0247ae51ab9d5668ada48</anchor>
      <arglist>(ostream &amp;outs, const System &amp;s)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>TotalMass</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ad1d8e856a8bab272c52582ef2bcda3dd</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MaxLength</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aff7f3ff197367b808c678ffab23daadc</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MaxRadius</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a2514c5516b4b1ec84a0593d2022566e9</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MinRadius</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>af62bccbd62ee66c1b3145cd291593f84</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t *</type>
      <name>RadialLimits</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aad5a851bb2efeb92e9d0e38a4dbd1a54</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>AverageRadius</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a3cd8e8c921ef68b36deb16ca06226c89</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>AverageSpeed</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a6554e41dd8e585bc53895fb60df27211</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>AverageVelocity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a02c9a0c987fed87a5e29e905038cc33c</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>AverageRadialVelocity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ab693fb1cde2c7f235a13f4aee1124661</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MaxRadialVelocity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>acc5977cdca9cc6f6322ffafefd9f4806</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MaxCircularVelocity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a2c2cbb2049d868519a90a7fede858153</anchor>
      <arglist>(bool cmframe=true) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>AverageDensity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a8c552c851ae7c783a6846886b6396dd3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>RegionDensity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a76bfb306c6026f2a304139cfa58d93c2</anchor>
      <arglist>(Double_t R, Double_t x, Double_t y, Double_t z) const </arglist>
    </member>
    <member kind="function">
      <type>Coord</type>
      <name>AngularMomentum</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>af4746eda2310a53b8ab79dce3b733c62</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Coord</type>
      <name>AngularMomentum</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aa776e4f1b922627bb78a3d89306c4522</anchor>
      <arglist>(Coord) const </arglist>
    </member>
    <member kind="function">
      <type>Coord</type>
      <name>AngularMomentum</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aa1ac1ee5a90ee817354748437b7ba0d9</anchor>
      <arglist>(Coordinate) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SortByID</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a805a3da6c40d1f660b5f833222013165</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SortByRadius</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aed5cb393db9538e850dd7db9c2cdeeaa</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SortByDistance</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a8cd5cfdd64bc6c3a375c0a069732333c</anchor>
      <arglist>(Coordinate c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SortByDistance</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a514888d05e2edca0f5f936e68a52ec1f</anchor>
      <arglist>(Coordinate c, int start, int n)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SortByEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ae26e69b8ea7491e07ec8be2987016ac5</anchor>
      <arglist>(Double_t G, Double_t eps)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SortByEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a4bfad3b86f150c491d9279afccb8a5b9</anchor>
      <arglist>(Double_t G, Double_t eps, int start, int n)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SortByDensity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a9b74a089ecd134f9846f232884e8d9fa</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>SortByDensity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a39441652159ce43a4902045144c2148b</anchor>
      <arglist>(int start, int n)</arglist>
    </member>
    <member kind="function">
      <type>Coord</type>
      <name>CM</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a45573b5f5529d4b5d816349499d4c5c0</anchor>
      <arglist>(Double_t tolerance=0.0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AdjustForCM</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a09377e523306f330b95c778ae63e48e1</anchor>
      <arglist>(Double_t tolerance=0.0)</arglist>
    </member>
    <member kind="function">
      <type>Coord</type>
      <name>CMVel</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>adabe1682b68df10112ea07c9a136dd1a</anchor>
      <arglist>(Double_t tolerance=0.0) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AdjustForCMVel</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aa00818e55f49911c81d6dec0958a243d</anchor>
      <arglist>(Double_t tolerance=0.0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AdjustPosition</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a6d82a74ac4067b7ea40a8861e10ee4d6</anchor>
      <arglist>(Coord)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AdjustPosition</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a3f0677eca5a6a9e17b58ba55fca31fcc</anchor>
      <arglist>(Coordinate)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AdjustVelocity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ab955207d46aa650ae2db1bac6bbbcb66</anchor>
      <arglist>(Coord)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AdjustVelocity</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aaa5c92b7e95afe8b70203960cb4278d1</anchor>
      <arglist>(Coordinate)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>KineticEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a3c3ca4fe980687aaf850af7059458c2f</anchor>
      <arglist>(Int_t i) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>KineticEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a34e0418fca052b017ce1ba783774d835</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>KineticEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ae1c893bd9170ea7e2ed2da5fe3d74595</anchor>
      <arglist>(Coord, Int_t) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>KineticEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a6a464c1fa0b538898249cbed595441a4</anchor>
      <arglist>(Coord) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PotentialEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ab3763e2839f276c8d6ee7c2f53e23649</anchor>
      <arglist>(Int_t i, Int_t j, Double_t eps=0.0) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PotentialEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a0cb4e8f58c8c5f8d8db21df78ca2b37e</anchor>
      <arglist>(Int_t i, Double_t eps=0.0) const </arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PotentialEnergy</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a7d73932eec5e6f5a642cbc1ccb7e3a04</anchor>
      <arglist>(Double_t eps=0.0, bool sphericalapprox=false) const </arglist>
    </member>
    <member kind="function">
      <type>Int_t</type>
      <name>FindParticle</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ace301b05392e7ccab56320374c38f6d8</anchor>
      <arglist>(const Particle &amp;p) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AddParticle</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ad2af64c8008adbdbc6e64d59e5547116</anchor>
      <arglist>(Particle part)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>RemoveParticle</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>addb46f1aa30205fe69d8cc4b0a5da4a7</anchor>
      <arglist>(Particle part)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>RemoveParticle</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>aca3fa8a69ddb7bf624ba28faaad4acd5</anchor>
      <arglist>(Int_t i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AddSystem</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a4961f58090fa4cb88b5be5015f283aa5</anchor>
      <arglist>(const System &amp;s, const Double_t *o)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>AddSystem</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a202c2ddf1bdd5f125098ba0f89847caa</anchor>
      <arglist>(const System &amp;s, const Double_t *o, const Double_t *vo)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ExtractSphere</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>af4756e63b898bb1ec9404e2dcdff0c31</anchor>
      <arglist>(Double_t radius)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Write</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a3dc72719dfc07b206b34c652d61b8ec5</anchor>
      <arglist>(ofstream &amp;F) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Write</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>ab80d43d340bb5c5d017bf2dd5d7ab6c8</anchor>
      <arglist>(char *Filename) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>Write</name>
      <anchorfile>classNBody_1_1System.html</anchorfile>
      <anchor>a3ae136ff83cf04dd1b950a9c81f1b926</anchor>
      <arglist>(FILE *stream=stdout) const </arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>Cosmology</name>
    <filename>namespaceCosmology.html</filename>
    <member kind="function">
      <type>Double_t</type>
      <name>aHIntFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ad54a65fe4587f2e0a6b59bbce0fdc95e</anchor>
      <arglist>(Double_t a, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>aHIntFuncGSL</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a2f1cb4bb3d240599ed6633aa55858069</anchor>
      <arglist>(double a, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>aHIntFuncGSLMonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a8a3d0e95fcf0b95e41ac04b17d633acd</anchor>
      <arglist>(double a, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>aintegral</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>aca839619c098dab3d2d7e25e66f7a41d</anchor>
      <arglist>(const Double_t a1, const Double_t a2, const Double_t omega, const Double_t xlambda)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a269b55d49c4b8c5f81e22bb44798439c</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>CEH98diff</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>adbcb8747ff665b7f0bacdfecce9adc63</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GrowthFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>aea982866216250805ede1d530a632b5b</anchor>
      <arglist>(const Double_t a, const Double_t Ho, const Double_t om, const Double_t ode, const Double_t ok, Double_t alpha=0.0, Double_t wde=-1.0)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GrowthFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a701cbd7c6a984ed69d0256b9f4c00bd4</anchor>
      <arglist>(const Double_t a, const Double_t om, const Double_t ola)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>GrowthFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a55eac79c3a605b928ec60eadd742d5e6</anchor>
      <arglist>(const Double_t a, const Double_t om, const Double_t ode, const Double_t ok, Double_t alpha, Double_t wde)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>HubbleFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>aad1be2f6641090f6757308631fb89e14</anchor>
      <arglist>(const Double_t a, const Double_t Ho, const Double_t om, const Double_t ola)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>HubbleFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a225aeef989bca704ce6b29bb20350640</anchor>
      <arglist>(const Double_t a, const Double_t Ho, const Double_t om, const Double_t ode, const Double_t ok, Double_t alpha, Double_t wde)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>HubbleIntFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a9fe760bfcd7471d2e11747ec4dfd6f68</anchor>
      <arglist>(Double_t a, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>HubbleIntFuncGSL</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a2bdb4d7670764893b4939cc8939fb356</anchor>
      <arglist>(double a, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>HubbleIntFuncGSLMonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>adec37f0116bfce0e36e634f8c4c7e135</anchor>
      <arglist>(double a, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegralkkPower</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a1db6a3f36e27864f82e0f77c38024393</anchor>
      <arglist>(int powerspectype, void *params, int integrationtype, bool logint)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegralSigma</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>afdcdec92d70adf7517fa8e192da28ef7</anchor>
      <arglist>(int powerspectype, Double_t R, void *params, int integrationtype, bool logint)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>LEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a04e6ad2e407977cc24ea195e9608344f</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>LEH98diff</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a112c703ccbb28974bc7177075ba2a0b2</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>neffEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a21f2aa6ccc0c33681417aed69d7f2a2c</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>neffWDMEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a2a2540ccb558a0585d833ec65cefb754</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Theta27, const Double_t Rd)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PBBKS</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>aad648da6daf2f7edb66289973eeedbbf</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Amp, const Double_t Omegam, const Double_t h, const Double_t Omegab)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PBBKSint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a15f2e2aafca0943851390548b82d3181</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PBBKSintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>af5e737552805567ad00e2f643476cbce</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PBBKSlogint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a59908f21fd15b932f1887f52cffa4e7b</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PBBKSlogintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a7c8948083a12484e110d9aa0342ac50b</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ab57708819b9a94327a84f9a8533447af</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t deltaH, const Double_t clight, const Double_t Ho, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a4dc841298fd503e20c695595c20df522</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Amp, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PEH98int</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a91da985b0dce318c620c8752da4f0867</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PEH98intgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a362daec612b3644d387bb0790bf216b7</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PEH98logint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a1fd1d0312797e78d0c81b412f711b5a7</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PEH98logintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ad48a21b0f47e026dff4fdec31807bf9f</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PFromFile</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a449ab31b3b09ea16112fcdb1b33e7474</anchor>
      <arglist>(const Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PGreen04</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a23ecc8f326d6040c773a6b5b2f2c1c14</anchor>
      <arglist>(const Double_t k, const Double_t Amp, const Double_t wm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PGreen04int</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a0a769ad4f78be1caabfddb824d69b634</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PGreen04intgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a749e1c3f250dc2e041012416f4acc3a7</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PGreen04logint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a96ee6732f3299f1bcc5a29c76122efee</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PGreen04logintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a619e6d7b3b9781419cc421a02c61d69d</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PK</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ad19289dc0b293f5610752d0eb217d9d8</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Amp)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PKint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a447874d1be6d0d5ad7c8b5250a041a58</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PKintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a7442fb1ee5be26ba69dbd8ac3b3d8ae5</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PKlogint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a576f1c4fb9280de51327a79e7e699b6a</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PKlogintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a7b1762a0283dee484cd5fa27c0e36442</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PWDMEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ac2dfabffb61c67897592a80cb8392c6d</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Amp, const Double_t Omegam, const Double_t h, const Double_t Theta27, const Double_t Rd)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PWDMEH98int</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a3b11af9a8a327f3276ddb85c36c8dc0a</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PWDMEH98intgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>adb48d2dabf29e33254083dc72d712392</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PWDMEH98logint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a62718d753953db960a604cdbd44a7612</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>PWDMEH98logintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a5712d5ea7742344da26917ffd865a998</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>qBBSK</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a0eba524990c82b0326851b21db328a4d</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Omegab)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>qEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>aec718bb6cc3bb0bd15cfbe4e2a680fac</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>ScaleFactorFunc</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a7bfec341af541d5b5108f74fc1087d33</anchor>
      <arglist>(const Double_t a, const Double_t om, const Double_t xla)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPBBKSint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a57479f00b1098f16bc9fa65683afa8e1</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPBBKSintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a34ae2e8704ba5160fe07648374e615a1</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPBBKSlogint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a3e4aff04bf9b71811838976106e03f36</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPBBKSlogintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a42163bd76913b101b18fe2fa8bba998d</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPEH98int</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a036565e1f1933766443b574c835caa6e</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPEH98intgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a1e32599583cbf58e2d6939edb1cef3b7</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPEH98logint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a587cb6c3c6e3b7666ef423a9c02fde33</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPEH98logintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a603facf330094e9ddc8bbddee9c16c71</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPGreen04int</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a27e7d4d22fead09731ccc9a49ba1091b</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPGreen04intgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ae150998c5eb14d2f3ae7ff6f2ebfc55c</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPGreen04logint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ab4b312eb3132c03d73ac85234a73d59d</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPGreen04logintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a91a6a037350785c0c0baf6d192a9d52e</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPKint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a47a816ef67cf310979aafb1fdf713cbd</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPKintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a69f5bc18907a19c212dede36182e5160</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPKlogint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a09dee238033d6b83fc8734e334530d51</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPKlogintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a093c1a07b23163c1e6809a1c40c9530e</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPWDMEH98int</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a384971b260d3ebb5cae2b0948fdbfe90</anchor>
      <arglist>(Double_t k, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPWDMEH98intgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a6a032be1299ad331a3a2b876abae086b</anchor>
      <arglist>(double *ka, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>sigmaPWDMEH98logint</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a7c9c34451e8ea7fbac6eaa090d5fd6d3</anchor>
      <arglist>(Double_t logk, void *params)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sigmaPWDMEH98logintgslmonte</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a5bd7e310476241c5e1100c6294cdc853</anchor>
      <arglist>(double *logk, size_t dim, void *params)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>TBBKS</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a5368ed2e2aed7d2eb0d74c077f74dbb6</anchor>
      <arglist>(const Double_t k, const Double_t ns, const Double_t Omegam, const Double_t h, const Double_t Omegab)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>TEH98</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a8d3f6583b853ba18cb080d36b62e133c</anchor>
      <arglist>(const Double_t k, const Double_t Omegam, const Double_t h, const Double_t Theta27)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Timeh</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>a6f7793fa35350ec7c841b6b8a68da25a</anchor>
      <arglist>(const Double_t a, const Double_t om, const Double_t ola)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>Timet</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>ace1a8005ff6a3e92567df819694831b5</anchor>
      <arglist>(Double_t a, const Double_t om, const Double_t ola)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>WKR2</name>
      <anchorfile>namespaceCosmology.html</anchorfile>
      <anchor>af72eeaf43aa0d79d371c0fb6492e3a45</anchor>
      <arglist>(const Double_t k, const Double_t R)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>Math</name>
    <filename>namespaceMath.html</filename>
    <class kind="struct">Math::Coord</class>
    <class kind="class">Math::Coordinate</class>
    <class kind="class">Math::Coordinate2D</class>
    <class kind="class">Math::GMatrix</class>
    <class kind="struct">Math::math_function</class>
    <class kind="struct">Math::math_multidim_function</class>
    <class kind="class">Math::Matrix</class>
    <class kind="class">Math::Matrix2D</class>
    <member kind="function">
      <type>Double_t</type>
      <name>average</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a16bfec5be9cdda9d8b961158b64c22cf</anchor>
      <arglist>(Int_t n, Double_t *x)</arglist>
    </member>
    <member kind="function">
      <type>Double_t *</type>
      <name>average</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a30303faecc8bb40e09a853e032aef6a3</anchor>
      <arglist>(Int_t parnum, Int_t datnum, Double_t **x)</arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>covariance</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a8187a177d6158b91c01027c46da62a84</anchor>
      <arglist>(Int_t parnum, Int_t datnum, Double_t *xm, Double_t **x)</arglist>
    </member>
    <member kind="function">
      <type>GMatrix</type>
      <name>covariance</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a03718a69fdf5fdc8818fe91288fdffcd</anchor>
      <arglist>(Int_t parnum, Int_t datnum, Double_t **x)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>Factorial</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a2d0f9a58a07490b72ada749b89751631</anchor>
      <arglist>(int k)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>FitNonLinLS</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a7d7cb8f24f0499e5a3ec12286485aba6</anchor>
      <arglist>(const math_function fitfunc, const math_function *difffuncs, const int nparams, Double_t *params, GMatrix &amp;covar, const int npoints, const Double_t x[], const Double_t y[], GMatrix *Weight=NULL, Double_t error=1e-3, Double_t cl=0.95, int *fixparam=NULL, int binned=1, int maxiterations=1000, int iestimateerror=0)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateClosed</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a9426b041fe580caa17e4dae6b22f18ab</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b, const int n)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateData</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a7b662bfcf4a1b7da754e282578226742</anchor>
      <arglist>(const Double_t x[], const Double_t f[], const int lower, const int upper)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateQTrap</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>ae585139c1b1184047125bc4109072549</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b, const Double_t epsrel, int IMAX)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateRomberg</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a8654394ccfc022538a69c7af8e2dbd7f</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b, const Double_t epsrel, const int n, const int K)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateRombergMonte</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>ad462b2a23bfaa735abcb41fd97de868c</anchor>
      <arglist>(math_multidim_function *f, Double_t *a, Double_t *b, const Double_t epsrel, const int sampling, const int n, const int K, int interations=6)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateSimpleTrapezoidal</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a4ce1ced0a7bd11555915dc52bac4e87c</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b, const int n)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateSimpson</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>aa8d1091f877e020438189780fabc1448</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateTrap</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a5e5e5c79213f576436b082d589169dd9</anchor>
      <arglist>(const Double_t x[], const Double_t f[], const int a, const int b)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateTrapezoidal</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>af28d827090a237e07c9cf82bd2a68b11</anchor>
      <arglist>(const math_function *f, const Double_t a, const Double_t b, const int n)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>IntegrateVegasMonte</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>aa9ba450da02048058df4c70fee238405</anchor>
      <arglist>(math_multidim_function *fm, Double_t *a, Double_t *b, const int numintervals, Double_t *error, Double_t *chisq, int interations)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>IntegrateVegasMonte</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>aae99a70a37aa9c1848ff7f972d00d9e3</anchor>
      <arglist>(gsl_monte_function *gslfm, double *a, double *b, const int numintervals, double *error, double *chisq)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>nran2</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a431cdd9a33c8aa71fb9062287782c9e3</anchor>
      <arglist>(long *idum)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate2D</type>
      <name>operator*</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>aa998c693107b8b5538e0ef2bcab79498</anchor>
      <arglist>(Double_t a, const Coordinate2D &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>operator*</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a899ae08543e7360b976e44a7b69d9a1d</anchor>
      <arglist>(Double_t a, const Coordinate &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>Matrix</type>
      <name>operator*</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>af26e1bc6214f1a7e3a1164709c764cb8</anchor>
      <arglist>(Double_t a, const Matrix &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>Matrix2D</type>
      <name>operator*</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>abf0582063e3317432a5c627cbc0b5c3e</anchor>
      <arglist>(Double_t a, const Matrix2D &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>ac3d40a22a9698f7c68a2018c0971eb02</anchor>
      <arglist>(std::ostream &amp;outs, Coordinate c)</arglist>
    </member>
    <member kind="function">
      <type>std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a9ca86e51f9d04aa5d8c0e23b81c6141d</anchor>
      <arglist>(std::ostream &amp;outs, Coordinate2D c)</arglist>
    </member>
    <member kind="function">
      <type>std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a52310bb38b8cd3769957842c350fd94f</anchor>
      <arglist>(std::ostream &amp;stream, const Matrix2D &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a7d33cb24d8a9027be45c2b933736118d</anchor>
      <arglist>(std::ostream &amp;stream, const Matrix &amp;m)</arglist>
    </member>
    <member kind="function">
      <type>std::istream &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a03d140240a639aa68092a221a2f9aa4c</anchor>
      <arglist>(std::istream &amp;ins, Coordinate2D &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>std::istream &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>acc7c91fc21d8fd273fbf783e8ea29727</anchor>
      <arglist>(std::istream &amp;ins, Coordinate &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>OptimalBins</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>af1c80002b9862fc16637882937b62664</anchor>
      <arglist>(const int npoints, Double_t *points, Double_t xmin, Double_t xmax, Double_t *weights=NULL)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>PolyInt</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a261cfcce47b95f62b4cb4a77e4cedeea</anchor>
      <arglist>(Double_t *xa, Double_t *ya, int n, Double_t x, Double_t &amp;y, Double_t &amp;dy)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>ran2</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a9e6edeee1c9bf8200968ff55ce476a71</anchor>
      <arglist>(long *idum)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>rebin</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>aca549ace0630f4eee107187d3e5ccb41</anchor>
      <arglist>(Double_t rc, int nd, Double_t r[], Double_t xin[], Double_t xi[])</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>variance</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>ab65424600034618edf519814ca108eaa</anchor>
      <arglist>(Int_t n, Double_t xm, Double_t *x)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>variance</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a16b603565b32814503cbab2cd3686ef7</anchor>
      <arglist>(Int_t n, Double_t *x)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>vegas</name>
      <anchorfile>namespaceMath.html</anchorfile>
      <anchor>a7f438f7e967b6c3d9cfe88d8b370a5ee</anchor>
      <arglist>(math_multidim_function *fxn, Double_t regn[], int ndim, int init, unsigned long ncall, int itmx, int nprn, long int idum, Double_t *tgral, Double_t *sd, Double_t *chi2a)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>NBody</name>
    <filename>namespaceNBody.html</filename>
    <class kind="class">NBody::GasParticle</class>
    <class kind="class">NBody::KDTree</class>
    <class kind="class">NBody::LeafNode</class>
    <class kind="class">NBody::Node</class>
    <class kind="class">NBody::NPriorityQueue</class>
    <class kind="class">NBody::Particle</class>
    <class kind="class">NBody::PriorityQueue</class>
    <class kind="struct">NBody::smoothfunc</class>
    <class kind="class">NBody::SplitNode</class>
    <class kind="class">NBody::StarParticle</class>
    <class kind="class">NBody::System</class>
    <member kind="typedef">
      <type>int(*</type>
      <name>FOFcheckfunc</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ae68c57f56f1373303562475da97d3763</anchor>
      <arglist>)(Particle &amp;, Double_t *)</arglist>
    </member>
    <member kind="typedef">
      <type>int(*</type>
      <name>FOFcompfunc</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a4b80079c12e359efc6b9d5d65b469aa9</anchor>
      <arglist>)(Particle &amp;, Particle &amp;, Double_t *)</arglist>
    </member>
    <member kind="typedef">
      <type>Double_t</type>
      <name>DoublePos_t</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a234ec5660f1ec8cad4c118bf35f0ceca</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumeration">
      <type></type>
      <name>Orbit</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a392718bf49bf3a990af266aba7b82079</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>Radial</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a392718bf49bf3a990af266aba7b82079a8865a16d55e6b8111c5a511abf31d2eb</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>Isotropic</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a392718bf49bf3a990af266aba7b82079a208f65742765cbe9f1a8791aa40124d2</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>CalcEnergy</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ab06ccd1e9ac542d1220cc9f31bf8cd09</anchor>
      <arglist>(const System &amp;S, Double_t *E, Double_t &amp;Etot)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>CalcPotDirect</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a993e7fc42b40e66b1724424b6a905da0</anchor>
      <arglist>(System &amp;S, Double_t eps=0.0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>CalcPotShell</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a80628dea731cb2661ec3e8c576802dd9</anchor>
      <arglist>(System &amp;S, Double_t eps=0.0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>CalcPowSpectrum</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a8f7dff7ee7fc024bdf41ebe5ba9729d9</anchor>
      <arglist>(Double_t *rho, Double_t *p, int Ng, Double_t BoxLength)</arglist>
    </member>
    <member kind="function">
      <type>System *</type>
      <name>CubicGrid</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a503b86d87a69b1aa2e88056bcbd204b9</anchor>
      <arglist>(int N, Double_t BoxLength, Double_t Mtotal, Double_t time)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>DenCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a8573bb8cb7c6e930be22cc67752c56d3</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DensityBinEqualMass</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a33ec5f91802b46dbafee873789ae1c14</anchor>
      <arglist>(System &amp;S, int NumBins, Double_t *RBin, Double_t *DBin)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DensityBinLog</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ac9c4e8f661670737da3d02d478d209b2</anchor>
      <arglist>(const System &amp;S, int NumBins, Double_t *RBin, Double_t *DBin)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DensityBinLogEllip</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>af98d9605f8e1837b03a38780d9eb939f</anchor>
      <arglist>(const System &amp;S, int NumBins, Double_t *RBin, Double_t *DBin, Double_t q, Double_t s)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DensityBinNorm</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a9f216fbe3fb75d346ade3b7c5c0051cc</anchor>
      <arglist>(const System &amp;S, int NumBins, Double_t *RBin, Double_t *DBin)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DensityGridNGP</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a9e3bd604a9b342a16fae9e633c87a2a3</anchor>
      <arglist>(System &amp;S, Double_t *rho, int Ng, Double_t &amp;L)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DensitySmooth</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a50be98286eb2be7bcf0c5ecee4014290</anchor>
      <arglist>(System &amp;S, int BucketSize, int NumNearest)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DensitySmoothBall</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a288b85c5e87864a1b58565517e8631d3</anchor>
      <arglist>(System &amp;S, Double_t R)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DenStates</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a6505f11a44a67cad8b78721c70c0afe1</anchor>
      <arglist>(const System &amp;S, Double_t *g, const Double_t *Ebin, int NumBins, Orbit orbit)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DiffEnergyDistEqual</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>af5520cd37fea91613efa03b151f20445</anchor>
      <arglist>(System &amp;S, const Double_t *E, Double_t *Ebin, Double_t *Mbin, int NumBins)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DiffEnergyDistNorm</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>afc959d8071bb31d168af94437dd0bff1</anchor>
      <arglist>(const System &amp;S, const Double_t *E, Double_t *Ebin, Double_t *Mbin, int NumBins)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>FOF3d</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a793b2e06c4f5156974cab32d93177501</anchor>
      <arglist>(Particle &amp;a, Particle &amp;b, Double_t *params)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>FOF6d</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a60efce1bf35d4f2c094927fcfdcd09f3</anchor>
      <arglist>(Particle &amp;a, Particle &amp;b, Double_t *params)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>FOFVel</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a6f1486955ca76c8b4df2a3b8b9196b2e</anchor>
      <arglist>(Particle &amp;a, Particle &amp;b, Double_t *params)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphology</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a68397ab51b06718e9cf3ad77239b8a0e</anchor>
      <arglist>(const Int_t nbodies, Particle *p, Double_t &amp;q, Double_t &amp;s, Double_t Error, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphology</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a6c01dbea286d85b56ebcfdbe3ec5b845</anchor>
      <arglist>(const Int_t nbodies, Particle *p, Double_t &amp;q, Double_t &amp;s, Double_t Error, Matrix &amp;eigenvec, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphology</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aaf436aa1d0e2e8876cd3f14c4a169be3</anchor>
      <arglist>(System &amp;S, Double_t &amp;q, Double_t &amp;s, Double_t Error, int verbose=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphology</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aea80e105c9f44f5ae05535fefe9eef50</anchor>
      <arglist>(System &amp;S, Double_t &amp;q, Double_t &amp;s, Double_t Error, Matrix &amp;eigenvec, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphologyWithMass</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aee812d1a04cbad75feed79462dada9c6</anchor>
      <arglist>(const Int_t nbodies, Particle *p, Double_t &amp;q, Double_t &amp;s, Double_t Error, int mdenom=0, int verbose=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphologyWithMass</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a0bd395e9e22d1a7fb1958f91cd527adf</anchor>
      <arglist>(const Int_t nbodies, Particle *p, Double_t &amp;q, Double_t &amp;s, Double_t Error, Matrix &amp;eigenvec, int mdenom, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphologyWithMass</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aa7e450ff359b11a2a4be90915541603f</anchor>
      <arglist>(System &amp;S, Double_t &amp;q, Double_t &amp;s, Double_t Error, int mdenom, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetGlobalMorphologyWithMass</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a64b03d8b923b719741f49f89cb7b8513</anchor>
      <arglist>(System &amp;S, Double_t &amp;q, Double_t &amp;s, Double_t Error, Matrix &amp;eigenvec, int mdenom, int verbose)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetInertiaTensor</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a7ad395a51dc78c9a00653bb054145ea5</anchor>
      <arglist>(const Int_t nbodies, Particle *p, Double_t &amp;a, Double_t &amp;b, Double_t &amp;c, Matrix &amp;eigvec)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetInertiaTensor</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a87f960a3886d621231e641ed13b98a66</anchor>
      <arglist>(const Int_t n, Particle *p, Double_t &amp;a, Double_t &amp;b, Double_t &amp;c, Matrix &amp;eigenvec, Matrix &amp;I)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetInertiaTensor</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ab68b7ccda248a59fcaadda35c962569b</anchor>
      <arglist>(System &amp;S, Double_t &amp;a, Double_t &amp;b, Double_t &amp;c, Matrix &amp;eigvec)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>GetInertiaTensor</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a8b1d672d1bf4b868f1087809f54a882e</anchor>
      <arglist>(System &amp;S, Double_t &amp;a, Double_t &amp;b, Double_t &amp;c, Matrix &amp;eigenvec, Matrix &amp;I)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>NumberDensity</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a99346c99c7e4ce25703996edc297c277</anchor>
      <arglist>(const System &amp;S, const Double_t *E, Double_t *Ebin, Double_t *Jbin, Double_t *Nbin, int NumBins)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a6b5b0e13bcc9826748a451ecafabd1ba</anchor>
      <arglist>(ostream &amp;outs, const Particle &amp;p)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ac19a01089d346407d2d5e58bf1a4b749</anchor>
      <arglist>(ostream &amp;outs, const System &amp;S)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>Pnocheck</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aa88761eca5c73784bcb9298e54bd4fd6</anchor>
      <arglist>(Particle &amp;a, Double_t *params)</arglist>
    </member>
    <member kind="function">
      <type>Coordinate</type>
      <name>Rotate</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a2ce569f54674485a595c76cff99110c1</anchor>
      <arglist>(const Matrix R, Coordinate x)</arglist>
    </member>
    <member kind="function">
      <type>System *</type>
      <name>UniformSphere</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a126d0532db237e67a2d2b881b8c3765d</anchor>
      <arglist>(int N, Double_t Radius, Double_t Mtotal, Double_t time)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>WEpan</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a1f1b0b18be2cc3c6ec759e8dfc92f806</anchor>
      <arglist>(Double_t r, Double_t h)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>WGauss</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ac68030ab7c00989fe9e5666dfb61033e</anchor>
      <arglist>(Double_t r, Double_t h)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>WSPH</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ac9bcae63f04efc56464010976d84b580</anchor>
      <arglist>(Double_t r, Double_t h)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>WTH</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a073d954241c808cc19effdd73d2d1ebd</anchor>
      <arglist>(Double_t r, Double_t h)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a6d94b57122fcedb536c6b0c3885eca1e</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aec050d8b2be9712f3aac8386b639e4a5</anchor>
      <arglist>(const Double_t *v1, const Double_t *v2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>acfdf50f94d221815619e97831984154d</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a038a5a120237444305f0050dd3d408d4</anchor>
      <arglist>(const Double_t *v1, const Double_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PhaseDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ae693d0baeb7a12462799531c165418d0</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceProjSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a4536ce36d79eb5a3da958e19103d0fa1</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>afbe92863ed7180298c0e31faa8328198</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2, const Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a1e08c5a585ca16acdf190a1a05d08941</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricwithTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aa67cba5b7f8b27631fd824fd159a959c</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a0de40474be2bee286bd8a799f1fea32f</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>abe4c2dd5784b67ec2a7b0365bdae32bd</anchor>
      <arglist>(const Real_t *v1, const Real_t *v2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a87ca1c9704eba746948b1bdf92526c04</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a4c984ded3aff2bbf070daeb7dad067ea</anchor>
      <arglist>(const Real_t *v1, const Real_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PhaseDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a636f53a810a026a82738ecb198690ef9</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceProjSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a4c2b9f1bb732c80246da765970590ede</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a767e2056d87ce3e5b5b7d995cacc36d5</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2, const Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a99134f8221cd0922cf37c0ad07cb6304</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricwithTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>af2cf025d5cb8a4fa8f1fc935fd9f97c2</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a38e6b832ee8edf498f941e23f0b8e184</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a283d307a17468ba310500eb9cacb257c</anchor>
      <arglist>(const Double_t *v1, const Real_t *v2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a295f0d65cda2b9bcfa1fcfa5809210e9</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aceed47806f85296d2c8133fbaab9649e</anchor>
      <arglist>(const Double_t *v1, const Real_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PhaseDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a19e9ad32388807681ab4fb0d46ea2840</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceProjSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aa7b3973b7923f8086f3e3a0135395141</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a936161e0981e5aa3e6d1e26b8c2e012d</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2, const Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ae58dfacab58a84130ab85ba5a490b1f5</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricwithTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a8b8d97c6f5ba79d1d5b1016df95b9ec2</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a37c40a63a326927d2ed6689c4b9f643e</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a1a36b17b45103ac1420e06fba4d75624</anchor>
      <arglist>(const Real_t *v1, const Double_t *v2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a9f1247c5c83dbb5399095d44e4d0f1e1</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a438c6f4d2d62774d1411a844564925d3</anchor>
      <arglist>(const Real_t *v1, const Double_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PhaseDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a3ca04510e1122e8475d23fde1a96ead1</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceProjSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a3c54b28ff5ddbf31aa0e473bd34db24e</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aa719edd1c96a433cff1d5389bdc9dd1c</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2, const Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a50d4baedaab96f71c4d8e44170db0571</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricwithTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ad223d7bdc6917b8dc13f95331b99432d</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflection1D</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>acfcca7d0b1e4c399e4e01cfa6062d658</anchor>
      <arglist>(const Coordinate &amp;x0, Coordinate &amp;xp, Coordinate p, int k)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflection2D</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a53422b8fc7e733ec5c610281baaae7e2</anchor>
      <arglist>(const Coordinate &amp;x0, Coordinate &amp;xp, Coordinate p, int k1, int k2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflectionND</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a963d7d6b33ec1c4a582474e41e21bbcb</anchor>
      <arglist>(const Coordinate &amp;x0, Coordinate &amp;xp, Coordinate p, int ndim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflection1D</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a7424acea8a5356b27372a9578e9a155c</anchor>
      <arglist>(const Particle &amp;x0, Particle &amp;xp, Coordinate p, int k)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflection2D</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a98bd2d2266c8b4a93dcb130488884979</anchor>
      <arglist>(const Particle &amp;x0, Particle &amp;xp, Coordinate p, int k1, int k2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflectionND</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>af5007ee2e8aa728826863caa82d70f53</anchor>
      <arglist>(const Particle &amp;x0, Particle &amp;xp, Coordinate p, int ndim)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PIDCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a92da7f4fdf5e1a352cf7cb0ca181df0d</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>IDCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a26671a95adb3af98801f7623670951a2</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>RadCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a47ac977b109b91226271f35bc274bc7d</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>TypeCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ab9c42487f7eb685d00f481080ba6ffc7</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PotCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a72c2f2c6ad30d74216ad07683cda33b8</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a6d94b57122fcedb536c6b0c3885eca1e</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aec050d8b2be9712f3aac8386b639e4a5</anchor>
      <arglist>(const Double_t *v1, const Double_t *v2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>acfdf50f94d221815619e97831984154d</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a038a5a120237444305f0050dd3d408d4</anchor>
      <arglist>(const Double_t *v1, const Double_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PhaseDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ae693d0baeb7a12462799531c165418d0</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceProjSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a4536ce36d79eb5a3da958e19103d0fa1</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>afbe92863ed7180298c0e31faa8328198</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2, const Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a1e08c5a585ca16acdf190a1a05d08941</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricwithTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aa67cba5b7f8b27631fd824fd159a959c</anchor>
      <arglist>(const Double_t *p1, const Double_t *p2, const Double_t *v1, const Double_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a0de40474be2bee286bd8a799f1fea32f</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>abe4c2dd5784b67ec2a7b0365bdae32bd</anchor>
      <arglist>(const Real_t *v1, const Real_t *v2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a87ca1c9704eba746948b1bdf92526c04</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a4c984ded3aff2bbf070daeb7dad067ea</anchor>
      <arglist>(const Real_t *v1, const Real_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PhaseDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a636f53a810a026a82738ecb198690ef9</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceProjSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a4c2b9f1bb732c80246da765970590ede</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a767e2056d87ce3e5b5b7d995cacc36d5</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2, const Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a99134f8221cd0922cf37c0ad07cb6304</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricwithTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>af2cf025d5cb8a4fa8f1fc935fd9f97c2</anchor>
      <arglist>(const Real_t *p1, const Real_t *p2, const Real_t *v1, const Real_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a38e6b832ee8edf498f941e23f0b8e184</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a283d307a17468ba310500eb9cacb257c</anchor>
      <arglist>(const Double_t *v1, const Real_t *v2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a295f0d65cda2b9bcfa1fcfa5809210e9</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aceed47806f85296d2c8133fbaab9649e</anchor>
      <arglist>(const Double_t *v1, const Real_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PhaseDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a19e9ad32388807681ab4fb0d46ea2840</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceProjSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aa7b3973b7923f8086f3e3a0135395141</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a936161e0981e5aa3e6d1e26b8c2e012d</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2, const Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ae58dfacab58a84130ab85ba5a490b1f5</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricwithTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a8b8d97c6f5ba79d1d5b1016df95b9ec2</anchor>
      <arglist>(const Double_t *p1, const Real_t *p2, const Double_t *v1, const Real_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a37c40a63a326927d2ed6689c4b9f643e</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a1a36b17b45103ac1420e06fba4d75624</anchor>
      <arglist>(const Real_t *v1, const Double_t *v2, int dim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a9f1247c5c83dbb5399095d44e4d0f1e1</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>VelDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a438c6f4d2d62774d1411a844564925d3</anchor>
      <arglist>(const Real_t *v1, const Double_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PhaseDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a3ca04510e1122e8475d23fde1a96ead1</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>DistanceProjSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a3c54b28ff5ddbf31aa0e473bd34db24e</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>aa719edd1c96a433cff1d5389bdc9dd1c</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2, const Double_t *metric)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a50d4baedaab96f71c4d8e44170db0571</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>MetricwithTensorDistSqd</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ad223d7bdc6917b8dc13f95331b99432d</anchor>
      <arglist>(const Real_t *p1, const Double_t *p2, const Real_t *v1, const Double_t *v2, const Double_t *m0, const Double_t *m1, const GMatrix gm)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflection1D</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>acfcca7d0b1e4c399e4e01cfa6062d658</anchor>
      <arglist>(const Coordinate &amp;x0, Coordinate &amp;xp, Coordinate p, int k)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflection2D</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a53422b8fc7e733ec5c610281baaae7e2</anchor>
      <arglist>(const Coordinate &amp;x0, Coordinate &amp;xp, Coordinate p, int k1, int k2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflectionND</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a963d7d6b33ec1c4a582474e41e21bbcb</anchor>
      <arglist>(const Coordinate &amp;x0, Coordinate &amp;xp, Coordinate p, int ndim)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflection1D</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a7424acea8a5356b27372a9578e9a155c</anchor>
      <arglist>(const Particle &amp;x0, Particle &amp;xp, Coordinate p, int k)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflection2D</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a98bd2d2266c8b4a93dcb130488884979</anchor>
      <arglist>(const Particle &amp;x0, Particle &amp;xp, Coordinate p, int k1, int k2)</arglist>
    </member>
    <member kind="function">
      <type>Double_t</type>
      <name>PeriodicReflectionND</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>af5007ee2e8aa728826863caa82d70f53</anchor>
      <arglist>(const Particle &amp;x0, Particle &amp;xp, Coordinate p, int ndim)</arglist>
    </member>
    <member kind="typedef">
      <type>Double_t</type>
      <name>DoublePos_t</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a234ec5660f1ec8cad4c118bf35f0ceca</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PIDCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a92da7f4fdf5e1a352cf7cb0ca181df0d</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>IDCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a26671a95adb3af98801f7623670951a2</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>RadCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a47ac977b109b91226271f35bc274bc7d</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>TypeCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>ab9c42487f7eb685d00f481080ba6ffc7</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>PotCompare</name>
      <anchorfile>namespaceNBody.html</anchorfile>
      <anchor>a72c2f2c6ad30d74216ad07683cda33b8</anchor>
      <arglist>(const void *a, const void *b)</arglist>
    </member>
  </compound>
  <compound kind="group">
    <name>libNBody</name>
    <title>NBody library</title>
    <filename>group__libNBody.html</filename>
    <docanchor file="group__libNBody" title="License">license</docanchor>
    <docanchor file="group__libNBody" title="Description">Description</docanchor>
    <docanchor file="group__libNBody" title="Getting started">prelim</docanchor>
    <docanchor file="group__libNBody" title="Compilation">install</docanchor>
    <docanchor file="group__libNBody" title="Including the library">howtorun</docanchor>
  </compound>
  <compound kind="page">
    <name>libNBody-makeflags</name>
    <title>Makefile.config of libNBody</title>
    <filename>libNBody-makeflags</filename>
    <docanchor file="libNBody-makeflags" title="Parallel and computational flags">secmake1</docanchor>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>Reference documentation for libNBody</title>
    <filename>index</filename>
  </compound>
</tagfile>
