/*! \file allvars.cxx
 *  \brief provides instances of global variables
 */

#include "allvars.h"

//-- Structures and external variables
///external pointer to keep track of structure levels and parent
StrucLevelData *psldata;


///\name define write routines for the property data structure
//{@
///write (append) the properties data to an already open binary file
void PropData::WriteBinary(fstream &Fout, Options&opt){
    long long lval;
    long unsigned idval;
    unsigned int ival;
    double val, val3[3],val9[9];
    float fval;
    idval=haloid;
    Fout.write((char*)&idval,sizeof(idval));
    lval=ibound;
    Fout.write((char*)&lval,sizeof(idval));
    lval=iminpot;
    Fout.write((char*)&lval,sizeof(idval));
    lval=hostid;
    Fout.write((char*)&lval,sizeof(idval));
    idval=numsubs;
    Fout.write((char*)&idval,sizeof(idval));
    idval=num;
    Fout.write((char*)&idval,sizeof(idval));
    ival=stype;
    Fout.write((char*)&ival,sizeof(ival));
    if (opt.iKeepFOF==1) {
        idval=directhostid;
        Fout.write((char*)&idval,sizeof(idval));
        idval=hostfofid;
        Fout.write((char*)&idval,sizeof(idval));
    }

    val=gMvir;
    Fout.write((char*)&val,sizeof(val));

    for (int k=0;k<3;k++) val3[k]=gcm[k];
    Fout.write((char*)val3,sizeof(val)*3);
    for (int k=0;k<3;k++) val3[k]=gposmbp[k];
    Fout.write((char*)val3,sizeof(val)*3);
    for (int k=0;k<3;k++) val3[k]=gposminpot[k];
    Fout.write((char*)val3,sizeof(val)*3);
    for (int k=0;k<3;k++) val3[k]=gcmvel[k];
    Fout.write((char*)val3,sizeof(val)*3);
    for (int k=0;k<3;k++) val3[k]=gvelmbp[k];
    Fout.write((char*)val3,sizeof(val)*3);
    for (int k=0;k<3;k++) val3[k]=gvelminpot[k];
    Fout.write((char*)val3,sizeof(val)*3);


    val=gmass;
    Fout.write((char*)&val,sizeof(val));
    val=gMFOF;
    Fout.write((char*)&val,sizeof(val));
    val=gM200m;
    Fout.write((char*)&val,sizeof(val));
    val=gM200c;
    Fout.write((char*)&val,sizeof(val));
    val=gMBN98;
    Fout.write((char*)&val,sizeof(val));

    val=Efrac;
    Fout.write((char*)&val,sizeof(val));

    val=gRvir;
    Fout.write((char*)&val,sizeof(val));
    val=gsize;
    Fout.write((char*)&val,sizeof(val));
    val=gR200m;
    Fout.write((char*)&val,sizeof(val));
    val=gR200c;
    Fout.write((char*)&val,sizeof(val));
    val=gRBN98;
    Fout.write((char*)&val,sizeof(val));
    val=gRhalfmass;
    Fout.write((char*)&val,sizeof(val));
    val=gRmaxvel;
    Fout.write((char*)&val,sizeof(val));

    val=gmaxvel;
    Fout.write((char*)&val,sizeof(val));
    val=gsigma_v;
    Fout.write((char*)&val,sizeof(val));
    for (int k=0;k<3;k++) for (int n=0;n<3;n++) val9[k*3+n]=gveldisp(k,n);
    Fout.write((char*)val9,sizeof(val)*9);

    val=glambda_B;
    Fout.write((char*)&val,sizeof(val));
    for (int k=0;k<3;k++) val3[k]=gJ[k];
    Fout.write((char*)val3,sizeof(val)*3);

    val=gq;
    Fout.write((char*)&val,sizeof(val));
    val=gs;
    Fout.write((char*)&val,sizeof(val));
    for (int k=0;k<3;k++) for (int n=0;n<3;n++) val9[k*3+n]=geigvec(k,n);
    Fout.write((char*)val9,sizeof(val)*9);

    val=cNFW;
    Fout.write((char*)&val,sizeof(val));
    val=Krot;
    Fout.write((char*)&val,sizeof(val));
    val=T;
    Fout.write((char*)&val,sizeof(val));
    val=Pot;
    Fout.write((char*)&val,sizeof(val));

    val=RV_sigma_v;
    Fout.write((char*)&val,sizeof(val));
    for (int k=0;k<3;k++) for (int n=0;n<3;n++) val9[k*3+n]=RV_veldisp(k,n);
    Fout.write((char*)val9,sizeof(val)*9);

    val=RV_lambda_B;
    Fout.write((char*)&val,sizeof(val));
    for (int k=0;k<3;k++) val3[k]=RV_J[k];
    Fout.write((char*)val3,sizeof(val)*3);

    val=RV_q;
    Fout.write((char*)&val,sizeof(val));
    val=RV_s;
    Fout.write((char*)&val,sizeof(val));
    for (int k=0;k<3;k++) for (int n=0;n<3;n++) val9[k*3+n]=RV_eigvec(k,n);
    Fout.write((char*)val9,sizeof(val)*9);

    if (opt.iextrahalooutput) {
        for (int k=0;k<3;k++) val3[k]=gJ200m[k];
        Fout.write((char*)val3,sizeof(val)*3);
        for (int k=0;k<3;k++) val3[k]=gJ200c[k];
        Fout.write((char*)val3,sizeof(val)*3);
        for (int k=0;k<3;k++) val3[k]=gJBN98[k];
        Fout.write((char*)val3,sizeof(val)*3);
        if (opt.iInclusiveHalo>0) {
            val=gM200m_excl;
            Fout.write((char*)&val,sizeof(val));
            val=gM200c_excl;
            Fout.write((char*)&val,sizeof(val));
            val=gMBN98_excl;
            Fout.write((char*)&val,sizeof(val));

            val=gR200m_excl;
            Fout.write((char*)&val,sizeof(val));
            val=gR200c_excl;
            Fout.write((char*)&val,sizeof(val));
            val=gRBN98_excl;
            Fout.write((char*)&val,sizeof(val));

            for (int k=0;k<3;k++) val3[k]=gJ200m_excl[k];
            Fout.write((char*)val3,sizeof(val)*3);
            for (int k=0;k<3;k++) val3[k]=gJ200c_excl[k];
            Fout.write((char*)val3,sizeof(val)*3);
            for (int k=0;k<3;k++) val3[k]=gJBN98_excl[k];
            Fout.write((char*)val3,sizeof(val)*3);
        }
    }
#ifdef GASON
    idval=n_gas;
    Fout.write((char*)&idval,sizeof(idval));
    val=M_gas;
    Fout.write((char*)&val,sizeof(val));
    val=M_gas_rvmax;
    Fout.write((char*)&val,sizeof(val));
    val=M_gas_30kpc;
    Fout.write((char*)&val,sizeof(val));
    val=M_gas_500c;
    Fout.write((char*)&val,sizeof(val));

    for (int k=0;k<3;k++) val3[k]=cm_gas[k];
    Fout.write((char*)val3,sizeof(val)*3);
    for (int k=0;k<3;k++) val3[k]=cmvel_gas[k];
    Fout.write((char*)val3,sizeof(val)*3);

    val=Efrac_gas;
    Fout.write((char*)&val,sizeof(val));

    val=Rhalfmass_gas;
    Fout.write((char*)&val,sizeof(val));

    for (int k=0;k<3;k++) for (int n=0;n<3;n++) val9[k*3+n]=veldisp_gas(k,n);
    Fout.write((char*)val9,sizeof(val)*9);

    for (int k=0;k<3;k++) val3[k]=L_gas[k];
    Fout.write((char*)val3,sizeof(val)*3);

    val=q_gas;
    Fout.write((char*)&val,sizeof(val));
    val=s_gas;
    Fout.write((char*)&val,sizeof(val));
    for (int k=0;k<3;k++) for (int n=0;n<3;n++) val9[k*3+n]=eigvec_gas(k,n);
    Fout.write((char*)val9,sizeof(val)*9);

    val=Krot_gas;
    Fout.write((char*)&val,sizeof(val));
    val=Temp_mean_gas;
    Fout.write((char*)&val,sizeof(val));

#ifdef STARON
    val=Z_mean_gas;
    Fout.write((char*)&val,sizeof(val));
    val=SFR_gas;
    Fout.write((char*)&val,sizeof(val));
#endif

if (opt.iextragasoutput) {
    val=M_200mean_gas;
    Fout.write((char*)&val,sizeof(val));
    val=M_200crit_gas;
    Fout.write((char*)&val,sizeof(val));
    val=M_BN98_gas;
    Fout.write((char*)&val,sizeof(val));
    for (int k=0;k<3;k++) val3[k]=L_200mean_gas[k];
    Fout.write((char*)val3,sizeof(val)*3);
    for (int k=0;k<3;k++) val3[k]=L_200crit_gas[k];
    Fout.write((char*)val3,sizeof(val)*3);
    for (int k=0;k<3;k++) val3[k]=L_BN98_gas[k];
    Fout.write((char*)val3,sizeof(val)*3);
    if (opt.iInclusiveHalo>0) {
        val=M_200mean_excl_gas;
        Fout.write((char*)&val,sizeof(val));
        val=M_200crit_excl_gas;
        Fout.write((char*)&val,sizeof(val));
        val=M_BN98_excl_gas;
        Fout.write((char*)&val,sizeof(val));
        for (int k=0;k<3;k++) val3[k]=L_200mean_excl_gas[k];
        Fout.write((char*)val3,sizeof(val)*3);
        for (int k=0;k<3;k++) val3[k]=L_200crit_excl_gas[k];
        Fout.write((char*)val3,sizeof(val)*3);
        for (int k=0;k<3;k++) val3[k]=L_BN98_excl_gas[k];
        Fout.write((char*)val3,sizeof(val)*3);
    }
}
#endif

#ifdef STARON
    idval=n_star;
    Fout.write((char*)&idval,sizeof(idval));
    val=M_star;
    Fout.write((char*)&val,sizeof(val));
    val=M_star_rvmax;
    Fout.write((char*)&val,sizeof(val));
    val=M_star_30kpc;
    Fout.write((char*)&val,sizeof(val));
    val=M_star_500c;
    Fout.write((char*)&val,sizeof(val));

    for (int k=0;k<3;k++) val3[k]=cm_star[k];
    Fout.write((char*)val3,sizeof(val)*3);
    for (int k=0;k<3;k++) val3[k]=cmvel_star[k];
    Fout.write((char*)val3,sizeof(val)*3);

    val=Efrac_star;
    Fout.write((char*)&val,sizeof(val));

    val=Rhalfmass_star;
    Fout.write((char*)&val,sizeof(val));

    for (int k=0;k<3;k++) for (int n=0;n<3;n++) val9[k*3+n]=veldisp_star(k,n);
    Fout.write((char*)val9,sizeof(val)*9);

    for (int k=0;k<3;k++) val3[k]=L_star[k];
    Fout.write((char*)val3,sizeof(val)*3);

    val=q_star;
    Fout.write((char*)&val,sizeof(val));
    val=s_star;
    Fout.write((char*)&val,sizeof(val));
    for (int k=0;k<3;k++) for (int n=0;n<3;n++) val9[k*3+n]=eigvec_star(k,n);
    Fout.write((char*)val9,sizeof(val)*9);

    val=Krot_star;
    Fout.write((char*)&val,sizeof(val));
    val=t_mean_star;
    Fout.write((char*)&val,sizeof(val));
    val=Z_mean_star;
    Fout.write((char*)&val,sizeof(val));

    if (opt.iextrastaroutput) {
        val=M_200mean_star;
        Fout.write((char*)&val,sizeof(val));
        val=M_200crit_star;
        Fout.write((char*)&val,sizeof(val));
        val=M_BN98_star;
        Fout.write((char*)&val,sizeof(val));
        for (int k=0;k<3;k++) val3[k]=L_200mean_star[k];
        Fout.write((char*)val3,sizeof(val)*3);
        for (int k=0;k<3;k++) val3[k]=L_200crit_star[k];
        Fout.write((char*)val3,sizeof(val)*3);
        for (int k=0;k<3;k++) val3[k]=L_BN98_star[k];
        Fout.write((char*)val3,sizeof(val)*3);
        if (opt.iInclusiveHalo>0) {
            val=M_200mean_excl_star;
            Fout.write((char*)&val,sizeof(val));
            val=M_200crit_excl_star;
            Fout.write((char*)&val,sizeof(val));
            val=M_BN98_excl_star;
            Fout.write((char*)&val,sizeof(val));
            for (int k=0;k<3;k++) val3[k]=L_200mean_excl_star[k];
            Fout.write((char*)val3,sizeof(val)*3);
            for (int k=0;k<3;k++) val3[k]=L_200crit_excl_star[k];
            Fout.write((char*)val3,sizeof(val)*3);
            for (int k=0;k<3;k++) val3[k]=L_BN98_excl_star[k];
            Fout.write((char*)val3,sizeof(val)*3);
        }
    }
#endif

#ifdef BHON
    idval=n_bh;
    Fout.write((char*)&idval,sizeof(idval));
    val=M_bh;
    Fout.write((char*)&val,sizeof(val));
#endif
#ifdef HIGHRES
    idval=n_interloper;
    Fout.write((char*)&idval,sizeof(idval));
    val=M_interloper;
    Fout.write((char*)&val,sizeof(val));
#endif

#if defined(GASON) && defined(STARON)
    val=M_gas_sf;
    Fout.write((char*)&val,sizeof(val));
    val=Rhalfmass_gas_sf;
    Fout.write((char*)&val,sizeof(val));
    val=sigV_gas_sf;
    Fout.write((char*)&val,sizeof(val));
    for (int k=0;k<3;k++) val3[k]=L_gas_sf[k];
    Fout.write((char*)val3,sizeof(val)*3);
    val=Krot_gas_sf;
    Fout.write((char*)&val,sizeof(val));
    val=Temp_mean_gas_sf;
    Fout.write((char*)&val,sizeof(val));
    val=Z_mean_gas_sf;
    Fout.write((char*)&val,sizeof(val));

    if (opt.iextragasoutput) {
        val=M_200mean_gas_sf;
        Fout.write((char*)&val,sizeof(val));
        val=M_200crit_gas_sf;
        Fout.write((char*)&val,sizeof(val));
        val=M_BN98_gas_sf;
        Fout.write((char*)&val,sizeof(val));
        for (int k=0;k<3;k++) val3[k]=L_200mean_gas_sf[k];
        Fout.write((char*)val3,sizeof(val)*3);
        for (int k=0;k<3;k++) val3[k]=L_200crit_gas_sf[k];
        Fout.write((char*)val3,sizeof(val)*3);
        for (int k=0;k<3;k++) val3[k]=L_BN98_gas_sf[k];
        Fout.write((char*)val3,sizeof(val)*3);
        if (opt.iInclusiveHalo>0) {
            val=M_200mean_excl_gas_sf;
            Fout.write((char*)&val,sizeof(val));
            val=M_200crit_excl_gas_sf;
            Fout.write((char*)&val,sizeof(val));
            val=M_BN98_excl_gas_sf;
            Fout.write((char*)&val,sizeof(val));
            for (int k=0;k<3;k++) val3[k]=L_200mean_excl_gas_sf[k];
            Fout.write((char*)val3,sizeof(val)*3);
            for (int k=0;k<3;k++) val3[k]=L_200crit_excl_gas_sf[k];
            Fout.write((char*)val3,sizeof(val)*3);
            for (int k=0;k<3;k++) val3[k]=L_BN98_excl_gas_sf[k];
            Fout.write((char*)val3,sizeof(val)*3);
        }
    }

    val=M_gas_nsf;
    Fout.write((char*)&val,sizeof(val));
    val=Rhalfmass_gas_nsf;
    Fout.write((char*)&val,sizeof(val));
    val=sigV_gas_nsf;
    Fout.write((char*)&val,sizeof(val));
    for (int k=0;k<3;k++) val3[k]=L_gas_nsf[k];
    Fout.write((char*)val3,sizeof(val)*3);
    val=Krot_gas_nsf;
    Fout.write((char*)&val,sizeof(val));
    val=Temp_mean_gas_nsf;
    Fout.write((char*)&val,sizeof(val));
    val=Z_mean_gas_nsf;
    Fout.write((char*)&val,sizeof(val));

    if (opt.iextragasoutput) {
        val=M_200mean_gas_nsf;
        Fout.write((char*)&val,sizeof(val));
        val=M_200crit_gas_nsf;
        Fout.write((char*)&val,sizeof(val));
        val=M_BN98_gas_nsf;
        Fout.write((char*)&val,sizeof(val));
        for (int k=0;k<3;k++) val3[k]=L_200mean_gas_nsf[k];
        Fout.write((char*)val3,sizeof(val)*3);
        for (int k=0;k<3;k++) val3[k]=L_200crit_gas_nsf[k];
        Fout.write((char*)val3,sizeof(val)*3);
        for (int k=0;k<3;k++) val3[k]=L_BN98_gas_nsf[k];
        Fout.write((char*)val3,sizeof(val)*3);
        if (opt.iInclusiveHalo>0) {
            val=M_200mean_excl_gas_nsf;
            Fout.write((char*)&val,sizeof(val));
            val=M_200crit_excl_gas_nsf;
            Fout.write((char*)&val,sizeof(val));
            val=M_BN98_excl_gas_nsf;
            Fout.write((char*)&val,sizeof(val));
            for (int k=0;k<3;k++) val3[k]=L_200mean_excl_gas_nsf[k];
            Fout.write((char*)val3,sizeof(val)*3);
            for (int k=0;k<3;k++) val3[k]=L_200crit_excl_gas_nsf[k];
            Fout.write((char*)val3,sizeof(val)*3);
            for (int k=0;k<3;k++) val3[k]=L_BN98_excl_gas_nsf[k];
            Fout.write((char*)val3,sizeof(val)*3);
        }
    }
#endif

#ifdef GASON
    if (opt.gas_internalprop_names.size()+ opt.gas_chem_names.size()+opt.gas_chemproduction_names.size()>0) {
        for (auto &extrafield:opt.gas_internalprop_names)
        {
            val = hydroprop.GetInternalProperties(extrafield);
            Fout.write((char*)&val,sizeof(val));
        }
        for (auto &extrafield:opt.gas_chem_names)
        {
            val = hydroprop.GetChemistry(extrafield);
            Fout.write((char*)&val,sizeof(val));
        }
        for (auto &extrafield:opt.gas_chemproduction_names)
        {
            val = hydroprop.GetChemistryProduction(extrafield);
            Fout.write((char*)&val,sizeof(val));
        }
    }
#endif
#ifdef STARON
    if (opt.star_internalprop_names.size()+opt.star_chem_names.size()+opt.star_chemproduction_names.size()>0) {
        for (auto &extrafield:opt.star_internalprop_names)
        {
            val = starprop.GetInternalProperties(extrafield);
            Fout.write((char*)&val,sizeof(val));
        }
        for (auto &extrafield:opt.star_chem_names)
        {
            val = starprop.GetChemistry(extrafield);
            Fout.write((char*)&val,sizeof(val));
        }
        for (auto &extrafield:opt.star_chemproduction_names)
        {
            val = starprop.GetChemistryProduction(extrafield);
            Fout.write((char*)&val,sizeof(val));
        }
    }
#endif
#ifdef BHON
    if (opt.bh_internalprop_names.size()+opt.bh_chem_names.size()+opt.bh_chemproduction_names.size()>0) {
        for (auto &extrafield:opt.bh_internalprop_names)
        {
            val = bhprop.GetInternalProperties(extrafield);
            Fout.write((char*)&val,sizeof(val));
        }
        for (auto &extrafield:opt.bh_chem_names)
        {
            val = bhprop.GetChemistry(extrafield);
            Fout.write((char*)&val,sizeof(val));
        }
        for (auto &extrafield:opt.bh_chemproduction_names)
        {
            val = bhprop.GetChemistryProduction(extrafield);
            Fout.write((char*)&val,sizeof(val));
        }
    }
#endif
#ifdef EXTRADMON
    if (opt.extra_dm_internalprop_names.size()>0) {
        for (auto &extrafield:opt.extra_dm_internalprop_names)
        {
            val = extradmprop.GetExtraProperties(extrafield);
            Fout.write((char*)&val,sizeof(val));
        }
    }
#endif

    if (opt.iaperturecalc && opt.aperturenum>0){
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_npart[j],sizeof(int));
        }
#ifdef GASON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_npart_gas[j],sizeof(int));
        }
#ifdef STARON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_npart_gas_sf[j],sizeof(int));
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_npart_gas_nsf[j],sizeof(int));
        }
#endif
#endif
#ifdef STARON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_npart_star[j],sizeof(int));
        }
#endif
#ifdef BHON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_npart_bh[j],sizeof(int));
        }
#endif

        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_mass[j],sizeof(val));
        }
#ifdef GASON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_mass_gas[j],sizeof(val));
        }
#ifdef STARON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_mass_gas_sf[j],sizeof(val));
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_mass_gas_nsf[j],sizeof(val));
        }
#endif
#endif
#ifdef STARON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_mass_star[j],sizeof(val));
        }
#endif
#ifdef BHON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_mass_bh[j],sizeof(val));
        }
#endif
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_rhalfmass[j],sizeof(val));
        }
#ifdef GASON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_rhalfmass_gas[j],sizeof(val));
        }
#ifdef STARON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_rhalfmass_gas_sf[j],sizeof(val));
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_rhalfmass_gas_nsf[j],sizeof(val));
        }
#endif
#endif
#ifdef STARON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_rhalfmass_star[j],sizeof(val));
        }
#endif
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_veldisp[j],sizeof(val));
        }
#ifdef GASON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_veldisp_gas[j],sizeof(val));
        }
#ifdef STARON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_veldisp_gas_sf[j],sizeof(val));
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_veldisp_gas_nsf[j],sizeof(val));
        }
#endif
#endif
#ifdef STARON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_veldisp_star[j],sizeof(val));
        }
#endif
#if defined(GASON) && defined(STARON)
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_SFR_gas[j],sizeof(val));
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_Z_gas[j],sizeof(val));
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_Z_gas_sf[j],sizeof(val));
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_Z_gas_nsf[j],sizeof(val));
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_Z_star[j],sizeof(val));
        }
#endif
#ifdef GASON
        if (opt.gas_extraprop_aperture_calc) {
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.gas_internalprop_output_names_aperture){
                    val = aperture_properties_gas[j].GetInternalProperties(x);
                    Fout.write((char*)&val,sizeof(val));
                }
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.gas_chem_output_names_aperture){
                    val = aperture_properties_gas[j].GetChemistry(x);
                    Fout.write((char*)&val,sizeof(val));
                }
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.gas_chemproduction_output_names_aperture){
                    val = aperture_properties_gas[j].GetChemistryProduction(x);
                    Fout.write((char*)&val,sizeof(val));
                }
            }
        }
#endif
#ifdef STARON
        if (opt.star_extraprop_aperture_calc) {
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.star_internalprop_output_names_aperture){
                    val = aperture_properties_star[j].GetInternalProperties(x);
                    Fout.write((char*)&val,sizeof(val));
                }
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.star_chem_output_names_aperture){
                    val = aperture_properties_star[j].GetChemistry(x);
                    Fout.write((char*)&val,sizeof(val));
                }
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.star_chemproduction_output_names_aperture){
                    val = aperture_properties_star[j].GetChemistryProduction(x);
                    Fout.write((char*)&val,sizeof(val));
                }
            }
        }
#endif
#ifdef BHON
        if (opt.bh_extraprop_aperture_calc) {
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.bh_internalprop_output_names_aperture){
                    val = aperture_properties_bh[j].GetInternalProperties(x);
                    Fout.write((char*)&val,sizeof(val));
                }
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.bh_chem_output_names_aperture){
                    val = aperture_properties_bh[j].GetChemistry(x);
                    Fout.write((char*)&val,sizeof(val));
                }
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.bh_chemproduction_output_names_aperture){
                    val = aperture_properties_bh[j].GetChemistryProduction(x);
                    Fout.write((char*)&val,sizeof(val));
                }
            }
        }
#endif
#ifdef EXTRADMON
        if (opt.extra_dm_extraprop_aperture_calc) {
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.extra_dm_internalprop_output_names_aperture){
                    val = aperture_properties_extra_dm[j].GetExtraProperties(x);
                    Fout.write((char*)&val,sizeof(val));
                }
            }
        }
#endif
    }

    if (opt.iaperturecalc && opt.apertureprojnum>0){
        for (auto k=0;k<3;k++) {
            for (auto j=0;j<opt.apertureprojnum;j++) {
                Fout.write((char*)&aperture_mass_proj[j][k],sizeof(val));
            }
#ifdef GASON
            for (auto j=0;j<opt.apertureprojnum;j++) {
                Fout.write((char*)&aperture_mass_proj_gas[j][k],sizeof(val));
            }
#ifdef STARON
            for (auto j=0;j<opt.apertureprojnum;j++) {
                Fout.write((char*)&aperture_mass_proj_gas_sf[j][k],sizeof(val));
            }
            for (auto j=0;j<opt.apertureprojnum;j++) {
                Fout.write((char*)&aperture_mass_proj_gas_nsf[j][k],sizeof(val));
            }
#endif
#endif
#ifdef STARON
            for (auto j=0;j<opt.apertureprojnum;j++) {
                Fout.write((char*)&aperture_mass_proj_star[j][k],sizeof(val));
            }
#endif
            for (auto j=0;j<opt.apertureprojnum;j++) {
                Fout.write((char*)&aperture_rhalfmass_proj[j][k],sizeof(val));
            }
#ifdef GASON
            for (auto j=0;j<opt.apertureprojnum;j++) {
                Fout.write((char*)&aperture_rhalfmass_proj_gas[j][k],sizeof(val));
            }
#ifdef STARON
            for (auto j=0;j<opt.apertureprojnum;j++) {
                Fout.write((char*)&aperture_rhalfmass_proj_gas_sf[j][k],sizeof(val));
            }
            for (auto j=0;j<opt.apertureprojnum;j++) {
                Fout.write((char*)&aperture_rhalfmass_proj_gas_nsf[j][k],sizeof(val));
            }
#endif
#endif
#ifdef STARON
            for (auto j=0;j<opt.apertureprojnum;j++) {
                Fout.write((char*)&aperture_rhalfmass_proj_star[j][k],sizeof(val));
            }
#endif
#if defined(GASON) && defined(STARON)
            for (auto j=0;j<opt.apertureprojnum;j++) {
                Fout.write((char*)&aperture_SFR_proj_gas[j][k],sizeof(val));
            }
            for (auto j=0;j<opt.apertureprojnum;j++) {
                Fout.write((char*)&aperture_Z_proj_gas[j][k],sizeof(val));
            }
            for (auto j=0;j<opt.apertureprojnum;j++) {
                Fout.write((char*)&aperture_Z_proj_gas_sf[j][k],sizeof(val));
            }
            for (auto j=0;j<opt.apertureprojnum;j++) {
                Fout.write((char*)&aperture_Z_proj_gas_nsf[j][k],sizeof(val));
            }
            for (auto j=0;j<opt.apertureprojnum;j++) {
                Fout.write((char*)&aperture_Z_proj_star[j][k],sizeof(val));
            }
#endif
        }
    }
    if (opt.SOnum>0){
        for (auto j=0;j<opt.SOnum;j++) {
            Fout.write((char*)&SO_mass[j],sizeof(int));
        }
        for (auto j=0;j<opt.SOnum;j++) {
            Fout.write((char*)&SO_radius[j],sizeof(int));
        }
#ifdef GASON
        if (opt.iextragasoutput && opt.iextrahalooutput)
        for (auto j=0;j<opt.SOnum;j++) {
            Fout.write((char*)&SO_mass_gas[j],sizeof(int));
        }
#ifdef STARON
#endif
#endif
#ifdef STARON
        if (opt.iextrastaroutput && opt.iextrahalooutput)
        for (auto j=0;j<opt.SOnum;j++) {
            Fout.write((char*)&SO_mass_star[j],sizeof(int));
        }
#endif
    }
    if (opt.SOnum>0 && opt.iextrahalooutput){
        for (auto j=0;j<opt.SOnum;j++) {
            for (auto k=0;k<3;k++) Fout.write((char*)&SO_angularmomentum[j][k],sizeof(int));

        }
#ifdef GASON
        if (opt.iextragasoutput)
        for (auto j=0;j<opt.SOnum;j++) {
            for (auto k=0;k<3;k++) Fout.write((char*)&SO_angularmomentum_gas[j][k],sizeof(int));
        }
#ifdef STARON
#endif
#endif
#ifdef STARON
        if (opt.iextrastaroutput)
        for (auto j=0;j<opt.SOnum;j++) {
            for (auto k=0;k<3;k++) Fout.write((char*)&SO_angularmomentum_star[j][k],sizeof(int));
        }
#endif
    }
}

void PropData::WriteAscii(fstream &Fout, Options&opt){
    double val;
    Fout<<haloid<<" ";
    Fout<<ibound<<" ";
    Fout<<iminpot<<" ";
    Fout<<hostid<<" ";
    Fout<<numsubs<<" ";
    Fout<<num<<" ";
    Fout<<stype<<" ";
    if (opt.iKeepFOF==1) {
        Fout<<directhostid<<" ";
        Fout<<hostfofid<<" ";
    }
    Fout<<gMvir<<" ";
    for (int k=0;k<3;k++) Fout<<gcm[k]<<" ";
    for (int k=0;k<3;k++) Fout<<gposmbp[k]<<" ";
    for (int k=0;k<3;k++) Fout<<gposminpot[k]<<" ";
    for (int k=0;k<3;k++) Fout<<gcmvel[k]<<" ";
    for (int k=0;k<3;k++) Fout<<gvelmbp[k]<<" ";
    for (int k=0;k<3;k++) Fout<<gvelminpot[k]<<" ";
    Fout<<gmass<<" ";
    Fout<<gMFOF<<" ";
    Fout<<gM200m<<" ";
    Fout<<gM200c<<" ";
    Fout<<gMBN98<<" ";
    Fout<<Efrac<<" ";
    Fout<<gRvir<<" ";
    Fout<<gsize<<" ";
    Fout<<gR200m<<" ";
    Fout<<gR200c<<" ";
    Fout<<gRBN98<<" ";
    Fout<<gRhalfmass<<" ";
    Fout<<gRmaxvel<<" ";
    Fout<<gmaxvel<<" ";
    Fout<<gsigma_v<<" ";
    for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<gveldisp(k,n)<<" ";
    Fout<<glambda_B<<" ";
    for (int k=0;k<3;k++) Fout<<gJ[k]<<" ";
    Fout<<gq<<" ";
    Fout<<gs<<" ";
    for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<geigvec(k,n)<<" ";
    Fout<<cNFW<<" ";
    Fout<<Krot<<" ";
    Fout<<T<<" ";
    Fout<<Pot<<" ";

    Fout<<RV_sigma_v<<" ";
    for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<RV_veldisp(k,n)<<" ";
    Fout<<RV_lambda_B<<" ";
    for (int k=0;k<3;k++) Fout<<RV_J[k]<<" ";
    Fout<<RV_q<<" ";
    Fout<<RV_s<<" ";
    for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<RV_eigvec(k,n)<<" ";

    if (opt.iextrahalooutput) {
        for (int k=0;k<3;k++) Fout<<gJ200m[k]<<" ";
        for (int k=0;k<3;k++) Fout<<gJ200c[k]<<" ";
        for (int k=0;k<3;k++) Fout<<gJBN98[k]<<" ";
        if (opt.iInclusiveHalo>0) {
            Fout<<gM200m_excl<<" ";
            Fout<<gM200c_excl<<" ";
            Fout<<gMBN98_excl<<" ";
            Fout<<gR200m_excl<<" ";
            Fout<<gR200c_excl<<" ";
            Fout<<gRBN98_excl<<" ";
            for (int k=0;k<3;k++) Fout<<gJ200m_excl[k]<<" ";
            for (int k=0;k<3;k++) Fout<<gJ200c_excl[k]<<" ";
            for (int k=0;k<3;k++) Fout<<gJBN98_excl[k]<<" ";
        }
    }

#ifdef GASON
    Fout<<n_gas<<" ";
    Fout<<M_gas<<" ";
    Fout<<M_gas_rvmax<<" ";
    Fout<<M_gas_30kpc<<" ";
    //Fout<<M_gas_50kpc<<" ";
    Fout<<M_gas_500c<<" ";
    for (int k=0;k<3;k++) Fout<<cm_gas[k]<<" ";
    for (int k=0;k<3;k++) Fout<<cmvel_gas[k]<<" ";
    Fout<<Efrac_gas<<" ";
    Fout<<Rhalfmass_gas<<" ";
    for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<veldisp_gas(k,n)<<" ";
    for (int k=0;k<3;k++) Fout<<L_gas[k]<<" ";
    Fout<<q_gas<<" ";
    Fout<<s_gas<<" ";
    for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<eigvec_gas(k,n)<<" ";
    Fout<<Krot_gas<<" ";
    Fout<<Temp_mean_gas<<" ";
#ifdef STARON
    Fout<<Z_mean_gas<<" ";
    Fout<<SFR_gas<<" ";
#endif
if (opt.iextragasoutput) {
    Fout<<M_200mean_gas<<" ";
    Fout<<M_200crit_gas<<" ";
    Fout<<M_BN98_gas<<" ";
    for (int k=0;k<3;k++) Fout<<L_200mean_gas[k]<<" ";
    for (int k=0;k<3;k++) Fout<<L_200crit_gas[k]<<" ";
    for (int k=0;k<3;k++) Fout<<L_BN98_gas[k]<<" ";
    if (opt.iInclusiveHalo>0) {
        Fout<<M_200mean_excl_gas<<" ";
        Fout<<M_200crit_excl_gas<<" ";
        Fout<<M_BN98_excl_gas<<" ";
        for (int k=0;k<3;k++) Fout<<L_200mean_excl_gas[k]<<" ";
        for (int k=0;k<3;k++) Fout<<L_200crit_excl_gas[k]<<" ";
        for (int k=0;k<3;k++) Fout<<L_BN98_excl_gas[k]<<" ";
    }
}
#endif

#ifdef STARON
    Fout<<n_star<<" ";
    Fout<<M_star<<" ";
    Fout<<M_star_rvmax<<" ";
    Fout<<M_star_30kpc<<" ";
    //Fout<<M_star_50kpc<<" ";
    Fout<<M_star_500c<<" ";
    for (int k=0;k<3;k++) Fout<<cm_star[k]<<" ";
    for (int k=0;k<3;k++) Fout<<cmvel_star[k]<<" ";
    Fout<<Efrac_star<<" ";
    Fout<<Rhalfmass_star<<" ";
    for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<veldisp_star(k,n)<<" ";
    for (int k=0;k<3;k++) Fout<<L_star[k]<<" ";
    Fout<<q_star<<" ";
    Fout<<s_star<<" ";
    for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<eigvec_star(k,n)<<" ";
    Fout<<Krot_star<<" ";
    Fout<<t_mean_star<<" ";
    Fout<<Z_mean_star<<" ";
    if (opt.iextrastaroutput) {
        Fout<<M_200mean_star<<" ";
        Fout<<M_200crit_star<<" ";
        Fout<<M_BN98_star<<" ";
        for (int k=0;k<3;k++) Fout<<L_200mean_star[k]<<" ";
        for (int k=0;k<3;k++) Fout<<L_200crit_star[k]<<" ";
        for (int k=0;k<3;k++) Fout<<L_BN98_star[k]<<" ";
        if (opt.iInclusiveHalo>0) {
            Fout<<M_200mean_excl_star<<" ";
            Fout<<M_200crit_excl_star<<" ";
            Fout<<M_BN98_excl_star<<" ";
            for (int k=0;k<3;k++) Fout<<L_200mean_excl_star[k]<<" ";
            for (int k=0;k<3;k++) Fout<<L_200crit_excl_star[k]<<" ";
            for (int k=0;k<3;k++) Fout<<L_BN98_excl_star[k]<<" ";
        }
    }
#endif

#ifdef BHON
    Fout<<n_bh<<" ";
    Fout<<M_bh<<" ";
#endif
#ifdef HIGHRES
    Fout<<n_interloper<<" ";
    Fout<<M_interloper<<" ";
    if (opt.iextrainterloperoutput) {
        Fout<<M_200mean_interloper<<" ";
        Fout<<M_200crit_interloper<<" ";
        Fout<<M_BN98_interloper<<" ";
        if (opt.iInclusiveHalo>0) {
            Fout<<M_200mean_excl_interloper<<" ";
            Fout<<M_200crit_excl_interloper<<" ";
            Fout<<M_BN98_excl_interloper<<" ";
        }
    }
#endif

#if defined(GASON) && defined(STARON)
    Fout<<M_gas_sf<<" ";
    Fout<<Rhalfmass_gas_sf<<" ";
    Fout<<sigV_gas_sf<<" ";
    for (int k=0;k<3;k++) Fout<<L_gas_sf[k]<<" ";
    Fout<<Krot_gas_sf<<" ";
    Fout<<Temp_mean_gas_sf<<" ";
    Fout<<Z_mean_gas_sf<<" ";
    if (opt.iextragasoutput) {
        Fout<<M_200mean_gas_sf<<" ";
        Fout<<M_200crit_gas_sf<<" ";
        Fout<<M_BN98_gas_sf<<" ";
        for (int k=0;k<3;k++) Fout<<L_200mean_gas_sf[k]<<" ";
        for (int k=0;k<3;k++) Fout<<L_200crit_gas_sf[k]<<" ";
        for (int k=0;k<3;k++) Fout<<L_BN98_gas_sf[k]<<" ";
        if (opt.iInclusiveHalo>0) {
            Fout<<M_200mean_excl_gas_sf<<" ";
            Fout<<M_200crit_excl_gas_sf<<" ";
            Fout<<M_BN98_excl_gas_sf<<" ";
            for (int k=0;k<3;k++) Fout<<L_200mean_excl_gas_sf[k]<<" ";
            for (int k=0;k<3;k++) Fout<<L_200crit_excl_gas_sf[k]<<" ";
            for (int k=0;k<3;k++) Fout<<L_BN98_excl_gas_sf[k]<<" ";
        }
    }
    Fout<<M_gas_nsf<<" ";
    Fout<<Rhalfmass_gas_nsf<<" ";
    Fout<<sigV_gas_nsf<<" ";
    for (int k=0;k<3;k++) Fout<<L_gas_nsf[k]<<" ";
    Fout<<Krot_gas_nsf<<" ";
    Fout<<Temp_mean_gas_nsf<<" ";
    Fout<<Z_mean_gas_nsf<<" ";
    if (opt.iextragasoutput) {
        Fout<<M_200mean_gas_nsf<<" ";
        Fout<<M_200crit_gas_nsf<<" ";
        Fout<<M_BN98_gas_nsf<<" ";
        for (int k=0;k<3;k++) Fout<<L_200mean_gas_nsf[k]<<" ";
        for (int k=0;k<3;k++) Fout<<L_200crit_gas_nsf[k]<<" ";
        for (int k=0;k<3;k++) Fout<<L_BN98_gas_nsf[k]<<" ";
        if (opt.iInclusiveHalo>0) {
            Fout<<M_200mean_excl_gas_nsf<<" ";
            Fout<<M_200crit_excl_gas_nsf<<" ";
            Fout<<M_BN98_excl_gas_nsf<<" ";
            for (int k=0;k<3;k++) Fout<<L_200mean_excl_gas_nsf[k]<<" ";
            for (int k=0;k<3;k++) Fout<<L_200crit_excl_gas_nsf[k]<<" ";
            for (int k=0;k<3;k++) Fout<<L_BN98_excl_gas_nsf[k]<<" ";
        }
    }
#endif
#ifdef GASON
    if (opt.gas_internalprop_names.size()+opt.gas_chem_names.size()+opt.gas_chemproduction_names.size()>0) {
        for (auto &extrafield:opt.gas_internalprop_names)
            Fout<<hydroprop.GetInternalProperties(extrafield)<<" ";
        for (auto &extrafield:opt.gas_chem_names)
            Fout<<hydroprop.GetChemistry(extrafield)<<" ";
        for (auto &extrafield:opt.gas_chemproduction_names)
            Fout<<hydroprop.GetChemistryProduction(extrafield)<<" ";
    }
#endif
#ifdef STARON
    if (opt.star_internalprop_names.size()+opt.star_chem_names.size()+opt.star_chemproduction_names.size()>0) {
        for (auto &extrafield:opt.star_internalprop_names)
            Fout<<starprop.GetInternalProperties(extrafield)<<" ";
        for (auto &extrafield:opt.star_chem_names)
            Fout<<starprop.GetChemistry(extrafield)<<" ";
        for (auto &extrafield:opt.star_chemproduction_names)
            Fout<<starprop.GetChemistryProduction(extrafield)<<" ";
    }
#endif
#ifdef BHON
    if (opt.bh_internalprop_names.size()+opt.bh_chem_names.size()+opt.bh_chemproduction_names.size()>0) {
        for (auto &extrafield:opt.bh_internalprop_names)
            Fout<<bhprop.GetInternalProperties(extrafield)<<" ";
        for (auto &extrafield:opt.bh_chem_names)
            Fout<<bhprop.GetChemistry(extrafield)<<" ";
        for (auto &extrafield:opt.bh_chemproduction_names)
            Fout<<bhprop.GetChemistryProduction(extrafield)<<" ";
    }
#endif
#ifdef EXTRADMON
    if (opt.extra_dm_internalprop_names.size()>0) {
        for (auto &extrafield:opt.extra_dm_internalprop_names)
            Fout<<extradmprop.GetExtraProperties(extrafield)<<" ";
    }
#endif

    if (opt.iaperturecalc && opt.aperturenum>0){
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_npart[j]<<" ";
        }
#ifdef GASON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_npart_gas[j]<<" ";
        }
#ifdef STARON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_npart_gas_sf[j]<<" ";
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_npart_gas_nsf[j]<<" ";
        }
#endif
#endif
#ifdef STARON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_npart_star[j]<<" ";
        }
#endif
#ifdef BHON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_npart_bh[j]<<" ";
        }
#endif
#ifdef HIGHRES
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_npart_interloper[j]<<" ";
        }
#endif
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_mass[j]<<" ";
        }
#ifdef GASON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_mass_gas[j]<<" ";
        }
#ifdef STARON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_mass_gas_sf[j]<<" ";
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_mass_gas_nsf[j]<<" ";
        }
#endif
#endif
#ifdef STARON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_mass_star[j]<<" ";
        }
#endif
#ifdef BHON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_mass_bh[j]<<" ";
        }
#endif
#ifdef HIGHRES
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_mass_interloper[j]<<" ";
        }
#endif
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_rhalfmass[j]<<" ";
        }
#ifdef GASON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_rhalfmass_gas[j]<<" ";
        }
#ifdef STARON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_rhalfmass_gas_sf[j]<<" ";
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_rhalfmass_gas_nsf[j]<<" ";
        }
#endif
#endif
#ifdef STARON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_rhalfmass_star[j]<<" ";
        }
#endif
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_veldisp[j]<<" ";
        }
#ifdef GASON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_veldisp_gas[j]<<" ";
        }
#ifdef STARON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_veldisp_gas_sf[j]<<" ";
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_veldisp_gas_nsf[j]<<" ";
        }
#endif
#endif
#ifdef STARON
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_veldisp_star[j]<<" ";
        }
#endif
#if defined(GASON) && defined(STARON)
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_SFR_gas[j]<<" ";
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_Z_gas[j]<<" ";
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_Z_gas_sf[j]<<" ";
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_Z_gas_nsf[j]<<" ";
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_Z_star[j]<<" ";
        }
#endif

#ifdef GASON
        if (opt.gas_extraprop_aperture_calc) {
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.gas_internalprop_output_names_aperture){
                    val = aperture_properties_gas[j].GetInternalProperties(x);
                    Fout<<val<<" ";
                }
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.gas_chem_output_names_aperture){
                    val = aperture_properties_gas[j].GetChemistry(x);
                    Fout<<val<<" ";
                }
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.gas_chemproduction_output_names_aperture){
                    val = aperture_properties_gas[j].GetChemistryProduction(x);
                    Fout<<val<<" ";
                }
            }
        }
#endif
#ifdef STARON
        if (opt.star_extraprop_aperture_calc) {
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.star_internalprop_output_names_aperture){
                    val = aperture_properties_star[j].GetInternalProperties(x);
                    Fout<<val<<" ";
                }
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.star_chem_output_names_aperture){
                    val = aperture_properties_star[j].GetChemistry(x);
                    Fout<<val<<" ";
                }
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.star_chemproduction_output_names_aperture){
                    val = aperture_properties_star[j].GetChemistryProduction(x);
                    Fout<<val<<" ";
                }
            }
        }
#endif
#ifdef BHON
        if (opt.bh_extraprop_aperture_calc) {
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.bh_internalprop_output_names_aperture){
                    val = aperture_properties_bh[j].GetInternalProperties(x);
                    Fout<<val<<" ";
                }
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.bh_chem_output_names_aperture){
                    val = aperture_properties_bh[j].GetChemistry(x);
                    Fout<<val<<" ";
                }
            }
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.bh_chemproduction_output_names_aperture){
                    val = aperture_properties_bh[j].GetChemistryProduction(x);
                    Fout<<val<<" ";
                }
            }
        }
#endif
#ifdef EXTRADMON
        if (opt.extra_dm_extraprop_aperture_calc) {
            for (auto j=0;j<opt.aperturenum;j++) {
                for (auto &x:opt.extra_dm_internalprop_output_names_aperture){
                    val = aperture_properties_extra_dm[j].GetExtraProperties(x);
                    Fout<<val<<" ";
                }
            }
        }
#endif
    }
    if (opt.iaperturecalc && opt.apertureprojnum>0) {
        for (auto k=0;k<3;k++) {
        for (auto j=0;j<opt.apertureprojnum;j++){
            Fout<<aperture_mass_proj[j][k]<<" ";
        }
        #ifdef GASON
        for (auto j=0;j<opt.apertureprojnum;j++){
            Fout<<aperture_mass_proj_gas[j][k]<<" ";
        }
        #ifdef STARON
        for (auto j=0;j<opt.apertureprojnum;j++){
            Fout<<aperture_mass_proj_gas_sf[j][k]<<" ";
        }
        for (auto j=0;j<opt.apertureprojnum;j++){
            Fout<<aperture_mass_proj_gas_nsf[j][k]<<" ";
        }
        #endif
        #endif
        #ifdef STARON
        for (auto j=0;j<opt.apertureprojnum;j++){
            Fout<<aperture_mass_proj_star[j][k]<<" ";
        }
        #endif
        for (auto j=0;j<opt.apertureprojnum;j++){
            Fout<<aperture_rhalfmass_proj[j][k]<<" ";
        }
        #ifdef GASON
        for (auto j=0;j<opt.apertureprojnum;j++){
            Fout<<aperture_rhalfmass_proj_gas[j][k]<<" ";
        }
        #ifdef STARON
        for (auto j=0;j<opt.apertureprojnum;j++){
            Fout<<aperture_rhalfmass_proj_gas_sf[j][k]<<" ";
        }
        for (auto j=0;j<opt.apertureprojnum;j++){
            Fout<<aperture_rhalfmass_proj_gas_nsf[j][k]<<" ";
        }
        #endif
        #endif
        #ifdef STARON
        for (auto j=0;j<opt.apertureprojnum;j++){
            Fout<<aperture_rhalfmass_proj_star[j][k]<<" ";
        }
        #endif
        #if defined(GASON) && defined(STARON)
        for (auto j=0;j<opt.apertureprojnum;j++) {
            Fout<<aperture_SFR_proj_gas[j][k]<<" ";
        }
        for (auto j=0;j<opt.apertureprojnum;j++) {
            Fout<<aperture_Z_proj_gas[j][k]<<" ";
        }
        for (auto j=0;j<opt.apertureprojnum;j++) {
            Fout<<aperture_Z_proj_gas_sf[j][k]<<" ";
        }
        for (auto j=0;j<opt.apertureprojnum;j++) {
            Fout<<aperture_Z_proj_gas_nsf[j][k]<<" ";
        }
        for (auto j=0;j<opt.apertureprojnum;j++) {
            Fout<<aperture_Z_proj_star[j][k]<<" ";
        }
        #endif
        }
    }
    if (opt.SOnum>0){
        for (auto j=0;j<opt.SOnum;j++) {
            Fout<<SO_mass[j]<<" ";
        }
        for (auto j=0;j<opt.SOnum;j++) {
            Fout<<SO_radius[j]<<" ";
        }
#ifdef GASON
        if (opt.iextragasoutput && opt.iextrahalooutput)
        for (auto j=0;j<opt.SOnum;j++) {
            Fout<<SO_mass_gas[j]<<" ";
        }
#ifdef STARON
#endif
#endif
#ifdef STARON
        if (opt.iextrastaroutput && opt.iextrahalooutput)
        for (auto j=0;j<opt.SOnum;j++) {
            Fout<<SO_mass_star[j]<<" ";
        }
#endif
#ifdef HIGHRES
        if (opt.iextrainterloperoutput && opt.iextrahalooutput)
        for (auto j=0;j<opt.SOnum;j++) {
            Fout<<SO_mass_interloper[j]<<" ";
        }
#endif
    }
    if (opt.SOnum>0 && opt.iextrahalooutput){
        for (auto j=0;j<opt.SOnum;j++) {
            for (auto k=0;k<3;k++) Fout<<SO_angularmomentum[j][k]<<" ";
        }
#ifdef GASON
        if (opt.iextragasoutput)
        for (auto j=0;j<opt.SOnum;j++) {
            for (auto k=0;k<3;k++) Fout<<SO_angularmomentum_gas[j][k]<<" ";
        }
#ifdef STARON
#endif
#endif
#ifdef STARON
        if (opt.iextrastaroutput)
        for (auto j=0;j<opt.SOnum;j++) {
            for (auto k=0;k<3;k++) Fout<<SO_angularmomentum_star[j][k]<<" ";
        }
#endif
    }
    Fout<<endl;
}

PropDataHeader::PropDataHeader(Options&opt){
    int sizeval;
#ifdef USEHDF
    // vector<PredType> desiredproprealtype;
    vector<hid_t> hdfdesiredproprealtype;
    if (sizeof(Double_t)==sizeof(double)) hdfdesiredproprealtype.push_back(H5T_NATIVE_DOUBLE);
    else hdfdesiredproprealtype.push_back(H5T_NATIVE_FLOAT);

#endif
#ifdef USEADIOS
    vector<ADIOS_DATATYPES> desiredadiosproprealtype;
    if (sizeof(Double_t)==sizeof(double)) desiredadiosproprealtype.push_back(ADIOS_DATATYPES::adios_double);
    else desiredadiosproprealtype.push_back(ADIOS_DATATYPES::adios_real);
#endif

    headerdatainfo.push_back("ID");
    headerdatainfo.push_back("ID_mbp");
    headerdatainfo.push_back("ID_minpot");
    headerdatainfo.push_back("hostHaloID");
    headerdatainfo.push_back("numSubStruct");
    headerdatainfo.push_back("npart");
    headerdatainfo.push_back("Structuretype");
    if (opt.iKeepFOF==1){
        headerdatainfo.push_back("hostDirectHaloID");
        headerdatainfo.push_back("hostFOFID");
    }

    //if using hdf, store the type
#ifdef USEHDF
    hdfpredtypeinfo.push_back(H5T_NATIVE_ULONG);
    hdfpredtypeinfo.push_back(H5T_NATIVE_LONG);
    hdfpredtypeinfo.push_back(H5T_NATIVE_LONG);
    hdfpredtypeinfo.push_back(H5T_NATIVE_LONG);
    hdfpredtypeinfo.push_back(H5T_NATIVE_ULONG);
    hdfpredtypeinfo.push_back(H5T_NATIVE_ULONG);
    hdfpredtypeinfo.push_back(H5T_NATIVE_INT);
    if (opt.iKeepFOF==1){
        hdfpredtypeinfo.push_back(H5T_NATIVE_LONG);
        hdfpredtypeinfo.push_back(H5T_NATIVE_LONG);
    }
#endif
#ifdef USEADIOS
    adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_long);
    adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_long);
    adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_long);
    adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_long);
    adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_long);
    adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_long);
    adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_integer);
    if (opt.iKeepFOF==1){
        adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_long);
        adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_long);
    }
#endif

    headerdatainfo.push_back("Mvir");
    headerdatainfo.push_back("Xc");
    headerdatainfo.push_back("Yc");
    headerdatainfo.push_back("Zc");
    headerdatainfo.push_back("Xcmbp");
    headerdatainfo.push_back("Ycmbp");
    headerdatainfo.push_back("Zcmbp");
    headerdatainfo.push_back("Xcminpot");
    headerdatainfo.push_back("Ycminpot");
    headerdatainfo.push_back("Zcminpot");
    headerdatainfo.push_back("VXc");
    headerdatainfo.push_back("VYc");
    headerdatainfo.push_back("VZc");
    headerdatainfo.push_back("VXcmbp");
    headerdatainfo.push_back("VYcmbp");
    headerdatainfo.push_back("VZcmbp");
    headerdatainfo.push_back("VXcminpot");
    headerdatainfo.push_back("VYcminpot");
    headerdatainfo.push_back("VZcminpot");
    headerdatainfo.push_back("Mass_tot");
    headerdatainfo.push_back("Mass_FOF");
    headerdatainfo.push_back("Mass_200mean");
    headerdatainfo.push_back("Mass_200crit");
    headerdatainfo.push_back("Mass_BN98");
    headerdatainfo.push_back("Efrac");
    headerdatainfo.push_back("Rvir");
    headerdatainfo.push_back("R_size");
    headerdatainfo.push_back("R_200mean");
    headerdatainfo.push_back("R_200crit");
    headerdatainfo.push_back("R_BN98");
    headerdatainfo.push_back("R_HalfMass");
    headerdatainfo.push_back("Rmax");
    headerdatainfo.push_back("Vmax");
    headerdatainfo.push_back("sigV");
    headerdatainfo.push_back("veldisp_xx");
    headerdatainfo.push_back("veldisp_xy");
    headerdatainfo.push_back("veldisp_xz");
    headerdatainfo.push_back("veldisp_yx");
    headerdatainfo.push_back("veldisp_yy");
    headerdatainfo.push_back("veldisp_yz");
    headerdatainfo.push_back("veldisp_zx");
    headerdatainfo.push_back("veldisp_zy");
    headerdatainfo.push_back("veldisp_zz");
    headerdatainfo.push_back("lambda_B");
    headerdatainfo.push_back("Lx");
    headerdatainfo.push_back("Ly");
    headerdatainfo.push_back("Lz");
    headerdatainfo.push_back("q");
    headerdatainfo.push_back("s");
    headerdatainfo.push_back("eig_xx");
    headerdatainfo.push_back("eig_xy");
    headerdatainfo.push_back("eig_xz");
    headerdatainfo.push_back("eig_yx");
    headerdatainfo.push_back("eig_yy");
    headerdatainfo.push_back("eig_yz");
    headerdatainfo.push_back("eig_zx");
    headerdatainfo.push_back("eig_zy");
    headerdatainfo.push_back("eig_zz");
    headerdatainfo.push_back("cNFW");
    headerdatainfo.push_back("Krot");
    headerdatainfo.push_back("Ekin");
    headerdatainfo.push_back("Epot");

    //some properties within RVmax
    headerdatainfo.push_back("RVmax_sigV");
    headerdatainfo.push_back("RVmax_veldisp_xx");
    headerdatainfo.push_back("RVmax_veldisp_xy");
    headerdatainfo.push_back("RVmax_veldisp_xz");
    headerdatainfo.push_back("RVmax_veldisp_yx");
    headerdatainfo.push_back("RVmax_veldisp_yy");
    headerdatainfo.push_back("RVmax_veldisp_yz");
    headerdatainfo.push_back("RVmax_veldisp_zx");
    headerdatainfo.push_back("RVmax_veldisp_zy");
    headerdatainfo.push_back("RVmax_veldisp_zz");
    headerdatainfo.push_back("RVmax_lambda_B");
    headerdatainfo.push_back("RVmax_Lx");
    headerdatainfo.push_back("RVmax_Ly");
    headerdatainfo.push_back("RVmax_Lz");
    headerdatainfo.push_back("RVmax_q");
    headerdatainfo.push_back("RVmax_s");
    headerdatainfo.push_back("RVmax_eig_xx");
    headerdatainfo.push_back("RVmax_eig_xy");
    headerdatainfo.push_back("RVmax_eig_xz");
    headerdatainfo.push_back("RVmax_eig_yx");
    headerdatainfo.push_back("RVmax_eig_yy");
    headerdatainfo.push_back("RVmax_eig_yz");
    headerdatainfo.push_back("RVmax_eig_zx");
    headerdatainfo.push_back("RVmax_eig_zy");
    headerdatainfo.push_back("RVmax_eig_zz");

    if (opt.iextrahalooutput) {
        headerdatainfo.push_back("Lx_200mean");
        headerdatainfo.push_back("Ly_200mean");
        headerdatainfo.push_back("Lz_200mean");
        headerdatainfo.push_back("Lx_200crit");
        headerdatainfo.push_back("Ly_200crit");
        headerdatainfo.push_back("Lz_200crit");
        headerdatainfo.push_back("Lx_BN98");
        headerdatainfo.push_back("Ly_BN98");
        headerdatainfo.push_back("Lz_BN98");
        if (opt.iInclusiveHalo>0) {
            headerdatainfo.push_back("Mass_200mean_excl");
            headerdatainfo.push_back("Mass_200crit_excl");
            headerdatainfo.push_back("Mass_BN98_excl");
            headerdatainfo.push_back("R_200mean_excl");
            headerdatainfo.push_back("R_200crit_excl");
            headerdatainfo.push_back("R_BN98_excl");
            headerdatainfo.push_back("Lx_200mean_excl");
            headerdatainfo.push_back("Ly_200mean_excl");
            headerdatainfo.push_back("Lz_200mean_excl");
            headerdatainfo.push_back("Lx_200crit_excl");
            headerdatainfo.push_back("Ly_200crit_excl");
            headerdatainfo.push_back("Lz_200crit_excl");
            headerdatainfo.push_back("Lx_BN98_excl");
            headerdatainfo.push_back("Ly_BN98_excl");
            headerdatainfo.push_back("Lz_BN98_excl");
        }
    }

#ifdef USEHDF
    sizeval=hdfpredtypeinfo.size();
    // for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(desiredproprealtype[0]);
    for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
    sizeval=adiospredtypeinfo.size();
    for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif

#ifdef GASON
    headerdatainfo.push_back("n_gas");
#ifdef USEHDF
    // predtypeinfo.push_back(PredType::STD_U64LE);
    hdfpredtypeinfo.push_back(H5T_NATIVE_ULONG);
#endif
#ifdef USEADIOS
    adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_long);
#endif
    headerdatainfo.push_back("M_gas");
    headerdatainfo.push_back("M_gas_Rvmax");
    headerdatainfo.push_back("M_gas_30kpc");
    //headerdatainfo.push_back("M_gas_50kpc");
    headerdatainfo.push_back("M_gas_500c");
    headerdatainfo.push_back("Xc_gas");
    headerdatainfo.push_back("Yc_gas");
    headerdatainfo.push_back("Zc_gas");
    headerdatainfo.push_back("VXc_gas");
    headerdatainfo.push_back("VYc_gas");
    headerdatainfo.push_back("VZc_gas");
    headerdatainfo.push_back("Efrac_gas");
    headerdatainfo.push_back("R_HalfMass_gas");
    headerdatainfo.push_back("veldisp_xx_gas");
    headerdatainfo.push_back("veldisp_xy_gas");
    headerdatainfo.push_back("veldisp_xz_gas");
    headerdatainfo.push_back("veldisp_yx_gas");
    headerdatainfo.push_back("veldisp_yy_gas");
    headerdatainfo.push_back("veldisp_yz_gas");
    headerdatainfo.push_back("veldisp_zx_gas");
    headerdatainfo.push_back("veldisp_zy_gas");
    headerdatainfo.push_back("veldisp_zz_gas");
    headerdatainfo.push_back("Lx_gas");
    headerdatainfo.push_back("Ly_gas");
    headerdatainfo.push_back("Lz_gas");
    headerdatainfo.push_back("q_gas");
    headerdatainfo.push_back("s_gas");
    headerdatainfo.push_back("eig_xx_gas");
    headerdatainfo.push_back("eig_xy_gas");
    headerdatainfo.push_back("eig_xz_gas");
    headerdatainfo.push_back("eig_yx_gas");
    headerdatainfo.push_back("eig_yy_gas");
    headerdatainfo.push_back("eig_yz_gas");
    headerdatainfo.push_back("eig_zx_gas");
    headerdatainfo.push_back("eig_zy_gas");
    headerdatainfo.push_back("eig_zz_gas");
    headerdatainfo.push_back("Krot_gas");
    headerdatainfo.push_back("T_gas");
#ifdef STARON
    headerdatainfo.push_back("Zmet_gas");
    headerdatainfo.push_back("SFR_gas");
#endif
if (opt.iextragasoutput) {
    headerdatainfo.push_back("Mass_200mean_gas");
    headerdatainfo.push_back("Mass_200crit_gas");
    headerdatainfo.push_back("Mass_BN98_gas");
    headerdatainfo.push_back("Lx_200c_gas");
    headerdatainfo.push_back("Ly_200c_gas");
    headerdatainfo.push_back("Lz_200c_gas");
    headerdatainfo.push_back("Lx_200m_gas");
    headerdatainfo.push_back("Ly_200m_gas");
    headerdatainfo.push_back("Lz_200m_gas");
    headerdatainfo.push_back("Lx_BN98_gas");
    headerdatainfo.push_back("Ly_BN98_gas");
    headerdatainfo.push_back("Lz_BN98_gas");
    if (opt.iInclusiveHalo>0) {
        headerdatainfo.push_back("Mass_200mean_excl_gas");
        headerdatainfo.push_back("Mass_200crit_excl_gas");
        headerdatainfo.push_back("Mass_BN98_excl_gas");
        headerdatainfo.push_back("Lx_200c_excl_gas");
        headerdatainfo.push_back("Ly_200c_excl_gas");
        headerdatainfo.push_back("Lz_200c_excl_gas");
        headerdatainfo.push_back("Lx_200m_excl_gas");
        headerdatainfo.push_back("Ly_200m_excl_gas");
        headerdatainfo.push_back("Lz_200m_excl_gas");
        headerdatainfo.push_back("Lx_BN98_excl_gas");
        headerdatainfo.push_back("Ly_BN98_excl_gas");
        headerdatainfo.push_back("Lz_BN98_excl_gas");
    }
}
#ifdef USEHDF
    sizeval=hdfpredtypeinfo.size();
    // for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(desiredproprealtype[0]);
    for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
    sizeval=adiospredtypeinfo.size();
    for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
#endif

#ifdef STARON
    headerdatainfo.push_back("n_star");
#ifdef USEHDF
    // predtypeinfo.push_back(PredType::STD_U64LE);
    hdfpredtypeinfo.push_back(H5T_NATIVE_ULONG);
#endif
#ifdef USEADIOS
    adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_long);
#endif
    headerdatainfo.push_back("M_star");
    headerdatainfo.push_back("M_star_Rvmax");
    headerdatainfo.push_back("M_star_30kpc");
    //headerdatainfo.push_back("M_star_50kpc");
    headerdatainfo.push_back("M_star_500c");
    headerdatainfo.push_back("Xc_star");
    headerdatainfo.push_back("Yc_star");
    headerdatainfo.push_back("Zc_star");
    headerdatainfo.push_back("VXc_star");
    headerdatainfo.push_back("VYc_star");
    headerdatainfo.push_back("VZc_star");
    headerdatainfo.push_back("Efrac_star");
    headerdatainfo.push_back("R_HalfMass_star");
    headerdatainfo.push_back("veldisp_xx_star");
    headerdatainfo.push_back("veldisp_xy_star");
    headerdatainfo.push_back("veldisp_xz_star");
    headerdatainfo.push_back("veldisp_yx_star");
    headerdatainfo.push_back("veldisp_yy_star");
    headerdatainfo.push_back("veldisp_yz_star");
    headerdatainfo.push_back("veldisp_zx_star");
    headerdatainfo.push_back("veldisp_zy_star");
    headerdatainfo.push_back("veldisp_zz_star");
    headerdatainfo.push_back("Lx_star");
    headerdatainfo.push_back("Ly_star");
    headerdatainfo.push_back("Lz_star");
    headerdatainfo.push_back("q_star");
    headerdatainfo.push_back("s_star");
    headerdatainfo.push_back("eig_xx_star");
    headerdatainfo.push_back("eig_xy_star");
    headerdatainfo.push_back("eig_xz_star");
    headerdatainfo.push_back("eig_yx_star");
    headerdatainfo.push_back("eig_yy_star");
    headerdatainfo.push_back("eig_yz_star");
    headerdatainfo.push_back("eig_zx_star");
    headerdatainfo.push_back("eig_zy_star");
    headerdatainfo.push_back("eig_zz_star");
    headerdatainfo.push_back("Krot_star");
    headerdatainfo.push_back("tage_star");
    headerdatainfo.push_back("Zmet_star");
    if (opt.iextrastaroutput) {
        headerdatainfo.push_back("Mass_200mean_star");
        headerdatainfo.push_back("Mass_200crit_star");
        headerdatainfo.push_back("Mass_BN98_star");
        headerdatainfo.push_back("Lx_200c_star");
        headerdatainfo.push_back("Ly_200c_star");
        headerdatainfo.push_back("Lz_200c_star");
        headerdatainfo.push_back("Lx_200m_star");
        headerdatainfo.push_back("Ly_200m_star");
        headerdatainfo.push_back("Lz_200m_star");
        headerdatainfo.push_back("Lx_BN98_star");
        headerdatainfo.push_back("Ly_BN98_star");
        headerdatainfo.push_back("Lz_BN98_star");
        if (opt.iInclusiveHalo>0) {
            headerdatainfo.push_back("Mass_200mean_excl_star");
            headerdatainfo.push_back("Mass_200crit_excl_star");
            headerdatainfo.push_back("Mass_BN98_excl_star");
            headerdatainfo.push_back("Lx_200c_excl_star");
            headerdatainfo.push_back("Ly_200c_excl_star");
            headerdatainfo.push_back("Lz_200c_excl_star");
            headerdatainfo.push_back("Lx_200m_excl_star");
            headerdatainfo.push_back("Ly_200m_excl_star");
            headerdatainfo.push_back("Lz_200m_excl_star");
            headerdatainfo.push_back("Lx_BN98_excl_star");
            headerdatainfo.push_back("Ly_BN98_excl_star");
            headerdatainfo.push_back("Lz_BN98_excl_star");
        }
    }
#ifdef USEHDF
    sizeval=hdfpredtypeinfo.size();
    // for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(desiredproprealtype[0]);
    for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
    sizeval=adiospredtypeinfo.size();
    for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
#endif

#ifdef BHON
    headerdatainfo.push_back("n_bh");
#ifdef USEHDF
    // predtypeinfo.push_back(PredType::STD_U64LE);
    hdfpredtypeinfo.push_back(H5T_NATIVE_ULONG);
#endif
#ifdef USEADIOS
    adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_long);
#endif
    headerdatainfo.push_back("M_bh");
#ifdef USEHDF
    sizeval=hdfpredtypeinfo.size();
    // for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(desiredproprealtype[0]);
    for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
    sizeval=adiospredtypeinfo.size();
    for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
#endif


#ifdef HIGHRES
    headerdatainfo.push_back("n_interloper");
#ifdef USEHDF
    // predtypeinfo.push_back(PredType::STD_U64LE);
    hdfpredtypeinfo.push_back(H5T_NATIVE_ULONG);
#endif
#ifdef USEADIOS
    adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_long);
#endif
    headerdatainfo.push_back("M_interloper");
    if (opt.iextrainterloperoutput) {
        headerdatainfo.push_back("Mass_200mean_interloper");
        headerdatainfo.push_back("Mass_200crit_interloper");
        headerdatainfo.push_back("Mass_BN98_interloper");
        if (opt.iInclusiveHalo>0) {
            headerdatainfo.push_back("Mass_200mean_excl_interloper");
            headerdatainfo.push_back("Mass_200crit_excl_interloper");
            headerdatainfo.push_back("Mass_BN98_excl_interloper");
        }
    }
#ifdef USEHDF
    sizeval=hdfpredtypeinfo.size();
    // for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(desiredproprealtype[0]);
    for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
    sizeval=adiospredtypeinfo.size();
    for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
#endif

#if defined(GASON) && defined(STARON)
    headerdatainfo.push_back("M_gas_sf");
    headerdatainfo.push_back("R_HalfMass_gas_sf");
    headerdatainfo.push_back("sigV_gas_sf");
    headerdatainfo.push_back("Lx_gas_sf");
    headerdatainfo.push_back("Ly_gas_sf");
    headerdatainfo.push_back("Lz_gas_sf");
    headerdatainfo.push_back("Krot_gas_sf");
    headerdatainfo.push_back("T_gas_sf");
    headerdatainfo.push_back("Zmet_gas_sf");
    if (opt.iextragasoutput) {
        headerdatainfo.push_back("Mass_200mean_gas_sf");
        headerdatainfo.push_back("Mass_200crit_gas_sf");
        headerdatainfo.push_back("Mass_BN98_gas_sf");
        headerdatainfo.push_back("Lx_200c_gas_sf");
        headerdatainfo.push_back("Ly_200c_gas_sf");
        headerdatainfo.push_back("Lz_200c_gas_sf");
        headerdatainfo.push_back("Lx_200m_gas_sf");
        headerdatainfo.push_back("Ly_200m_gas_sf");
        headerdatainfo.push_back("Lz_200m_gas_sf");
        headerdatainfo.push_back("Lx_BN98_gas_sf");
        headerdatainfo.push_back("Ly_BN98_gas_sf");
        headerdatainfo.push_back("Lz_BN98_gas_sf");
        if (opt.iInclusiveHalo>0) {
            headerdatainfo.push_back("Mass_200mean_excl_gas_sf");
            headerdatainfo.push_back("Mass_200crit_excl_gas_sf");
            headerdatainfo.push_back("Mass_BN98_excl_gas_sf");
            headerdatainfo.push_back("Lx_200c_excl_gas_sf");
            headerdatainfo.push_back("Ly_200c_excl_gas_sf");
            headerdatainfo.push_back("Lz_200c_excl_gas_sf");
            headerdatainfo.push_back("Lx_200m_excl_gas_sf");
            headerdatainfo.push_back("Ly_200m_excl_gas_sf");
            headerdatainfo.push_back("Lz_200m_excl_gas_sf");
            headerdatainfo.push_back("Lx_BN98_excl_gas_sf");
            headerdatainfo.push_back("Ly_BN98_excl_gas_sf");
            headerdatainfo.push_back("Lz_BN98_excl_gas_sf");
        }
    }
    headerdatainfo.push_back("M_gas_nsf");
    headerdatainfo.push_back("R_HalfMass_gas_nsf");
    headerdatainfo.push_back("sigV_gas_nsf");
    headerdatainfo.push_back("Lx_gas_nsf");
    headerdatainfo.push_back("Ly_gas_nsf");
    headerdatainfo.push_back("Lz_gas_nsf");
    headerdatainfo.push_back("Krot_gas_nsf");
    headerdatainfo.push_back("T_gas_nsf");
    headerdatainfo.push_back("Zmet_gas_nsf");
    if (opt.iextragasoutput) {
        headerdatainfo.push_back("Mass_200mean_gas_nsf");
        headerdatainfo.push_back("Mass_200crit_gas_nsf");
        headerdatainfo.push_back("Mass_BN98_gas_nsf");
        headerdatainfo.push_back("Lx_200c_gas_nsf");
        headerdatainfo.push_back("Ly_200c_gas_nsf");
        headerdatainfo.push_back("Lz_200c_gas_nsf");
        headerdatainfo.push_back("Lx_200m_gas_nsf");
        headerdatainfo.push_back("Ly_200m_gas_nsf");
        headerdatainfo.push_back("Lz_200m_gas_nsf");
        headerdatainfo.push_back("Lx_BN98_gas_nsf");
        headerdatainfo.push_back("Ly_BN98_gas_nsf");
        headerdatainfo.push_back("Lz_BN98_gas_nsf");
        if (opt.iInclusiveHalo>0) {
            headerdatainfo.push_back("Mass_200mean_excl_gas_nsf");
            headerdatainfo.push_back("Mass_200crit_excl_gas_nsf");
            headerdatainfo.push_back("Mass_BN98_excl_gas_nsf");
            headerdatainfo.push_back("Lx_200c_excl_gas_nsf");
            headerdatainfo.push_back("Ly_200c_excl_gas_nsf");
            headerdatainfo.push_back("Lz_200c_excl_gas_nsf");
            headerdatainfo.push_back("Lx_200m_excl_gas_nsf");
            headerdatainfo.push_back("Ly_200m_excl_gas_nsf");
            headerdatainfo.push_back("Lz_200m_excl_gas_nsf");
            headerdatainfo.push_back("Lx_BN98_excl_gas_nsf");
            headerdatainfo.push_back("Ly_BN98_excl_gas_nsf");
            headerdatainfo.push_back("Lz_BN98_excl_gas_nsf");
        }
    }
#ifdef USEHDF
    sizeval=hdfpredtypeinfo.size();
    // for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(desiredproprealtype[0]);
    for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
    sizeval=adiospredtypeinfo.size();
    for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
#endif

    //if extra hydro properties are calculated
#ifdef GASON
    if (opt.gas_internalprop_names.size()+opt.gas_chem_names.size()+opt.gas_chemproduction_names.size() > 0)
    {
        for (auto x:opt.gas_internalprop_output_names) headerdatainfo.push_back(x+string("_gas"));
        for (auto x:opt.gas_chem_output_names) headerdatainfo.push_back(x+string("_gas"));
        for (auto x:opt.gas_chemproduction_output_names) headerdatainfo.push_back(x+string("_gas"));
#ifdef USEHDF
        sizeval=hdfpredtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
    }
#endif
#ifdef STARON
    if (opt.star_internalprop_names.size()+opt.star_chem_names.size()+opt.star_chemproduction_names.size() > 0)
    {
        for (auto x:opt.star_internalprop_output_names) headerdatainfo.push_back(x+string("_star"));
        for (auto x:opt.star_chem_output_names) headerdatainfo.push_back(x+string("_star"));
        for (auto x:opt.star_chemproduction_output_names) headerdatainfo.push_back(x+string("_star"));
#ifdef USEHDF
        sizeval=hdfpredtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
    }
#endif
#ifdef BHON
    if (opt.bh_internalprop_names.size()+opt.bh_chem_names.size()+opt.bh_chemproduction_names.size() > 0)
    {
        for (auto x:opt.bh_internalprop_output_names) headerdatainfo.push_back(x+string("_bh"));
        for (auto x:opt.bh_chem_output_names) headerdatainfo.push_back(x+string("_bh"));
        for (auto x:opt.bh_chemproduction_output_names) headerdatainfo.push_back(x+string("_bh"));
#ifdef USEHDF
        sizeval=hdfpredtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
    }
#endif
#ifdef EXTRADMON
    if (opt.extra_dm_internalprop_names.size() > 0)
    {
        for (auto x:opt.extra_dm_internalprop_output_names) headerdatainfo.push_back(x+string("_extra_dm"));
#ifdef USEHDF
        sizeval=hdfpredtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
    }
#endif

    //if aperture information calculated also include
    if (opt.iaperturecalc>0 && opt.aperturenum>0) {
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_npart_")+opt.aperture_names_kpc[i]+string("_kpc")));
#ifdef GASON
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_npart_gas_")+opt.aperture_names_kpc[i]+string("_kpc")));
#ifdef STARON
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_npart_gas_sf_")+opt.aperture_names_kpc[i]+string("_kpc")));
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_npart_gas_nsf_")+opt.aperture_names_kpc[i]+string("_kpc")));
#endif
#endif
#ifdef STARON
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_npart_star_")+opt.aperture_names_kpc[i]+string("_kpc")));
#endif
#ifdef BHON
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_npart_bh_")+opt.aperture_names_kpc[i]+string("_kpc")));
#endif
#ifdef HIGHRES
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_npart_interloper_")+opt.aperture_names_kpc[i]+string("_kpc")));
#endif
#ifdef USEHDF
        sizeval=hdfpredtypeinfo.size();
        // for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(PredType::STD_U32LE);
        for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(H5T_NATIVE_UINT);
#endif
#ifdef USEADIOS
        sizeval=adiospredtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(ADIOS_DATATYPES::adios_unsigned_int);
#endif
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_mass_")+opt.aperture_names_kpc[i]+string("_kpc")));
#ifdef GASON
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_mass_gas_")+opt.aperture_names_kpc[i]+string("_kpc")));
#ifdef STARON
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_mass_gas_sf_")+opt.aperture_names_kpc[i]+string("_kpc")));
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_mass_gas_nsf_")+opt.aperture_names_kpc[i]+string("_kpc")));
#endif
#endif
#ifdef STARON
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_mass_star_")+opt.aperture_names_kpc[i]+string("_kpc")));
#endif
#ifdef BHON
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_mass_bh_")+opt.aperture_names_kpc[i]+string("_kpc")));
#endif
#ifdef HIGHRES
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_mass_interloper_")+opt.aperture_names_kpc[i]+string("_kpc")));
#endif

        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_rhalfmass_")+opt.aperture_names_kpc[i]+string("_kpc")));
#ifdef GASON
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_rhalfmass_gas_")+opt.aperture_names_kpc[i]+string("_kpc")));
#ifdef STARON
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_rhalfmass_gas_sf_")+opt.aperture_names_kpc[i]+string("_kpc")));
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_rhalfmass_gas_nsf_")+opt.aperture_names_kpc[i]+string("_kpc")));
#endif
#endif
#ifdef STARON
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_rhalfmass_star_")+opt.aperture_names_kpc[i]+string("_kpc")));
#endif

        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_veldisp_")+opt.aperture_names_kpc[i]+string("_kpc")));
#ifdef GASON
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_veldisp_gas_")+opt.aperture_names_kpc[i]+string("_kpc")));
#ifdef STARON
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_veldips_gas_sf_")+opt.aperture_names_kpc[i]+string("_kpc")));
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_veldisp_gas_nsf_")+opt.aperture_names_kpc[i]+string("_kpc")));
#endif
#endif
#ifdef STARON
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_veldisp_star_")+opt.aperture_names_kpc[i]+string("_kpc")));
#endif
#if defined(GASON) && defined(STARON)
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_SFR_gas_")+opt.aperture_names_kpc[i]+string("_kpc")));
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_Zmet_gas_")+opt.aperture_names_kpc[i]+string("_kpc")));
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_Zmet_gas_sf_")+opt.aperture_names_kpc[i]+string("_kpc")));
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_Zmet_gas_nsf_")+opt.aperture_names_kpc[i]+string("_kpc")));
        for (auto i=0; i<opt.aperturenum;i++)
            headerdatainfo.push_back((string("Aperture_Zmet_star_")+opt.aperture_names_kpc[i]+string("_kpc")));
#endif

#ifdef USEHDF
        sizeval=hdfpredtypeinfo.size();
        // for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(desiredproprealtype[0]);
        for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
        sizeval=adiospredtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif

        //add extra property data in apertures
#ifdef GASON
        if (opt.gas_extraprop_aperture_calc) {
            for (auto i=0; i<opt.aperturenum;i++) {
                for (auto &x:opt.gas_internalprop_output_names_aperture) {
                    headerdatainfo.push_back((string("Aperture_")+x+string("_gas_")+opt.aperture_names_kpc[i]+string("_kpc")));
                }
            }
            for (auto i=0; i<opt.aperturenum;i++) {
                for (auto &x:opt.gas_chem_output_names_aperture) {
                    headerdatainfo.push_back((string("Aperture_")+x+string("_gas_")+opt.aperture_names_kpc[i]+string("_kpc")));
                }
            }
            for (auto i=0; i<opt.aperturenum;i++) {
                for (auto &x:opt.gas_chemproduction_output_names_aperture) {
                    headerdatainfo.push_back((string("Aperture_")+x+string("_gas_")+opt.aperture_names_kpc[i]+string("_kpc")));
                }
            }
#ifdef USEHDF
            sizeval=hdfpredtypeinfo.size();
            for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
            sizeval=adiospredtypeinfo.size();
            for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
        }
#endif
#ifdef STARON
        if (opt.star_extraprop_aperture_calc) {
            for (auto i=0; i<opt.aperturenum;i++) {
                for (auto &x:opt.star_internalprop_output_names_aperture) {
                    headerdatainfo.push_back((string("Aperture_")+x+string("_star_")+opt.aperture_names_kpc[i]+string("_kpc")));
                }
            }
            for (auto i=0; i<opt.aperturenum;i++) {
                for (auto &x:opt.star_chem_output_names_aperture) {
                    headerdatainfo.push_back((string("Aperture_")+x+string("_star_")+opt.aperture_names_kpc[i]+string("_kpc")));
                }
            }
            for (auto i=0; i<opt.aperturenum;i++) {
                for (auto &x:opt.star_chemproduction_output_names_aperture) {
                    headerdatainfo.push_back((string("Aperture_")+x+string("_star_")+opt.aperture_names_kpc[i]+string("_kpc")));
                }
            }
#ifdef USEHDF
            sizeval=hdfpredtypeinfo.size();
            for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
            sizeval=adiospredtypeinfo.size();
            for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
        }
#endif
#ifdef BHON
        if (opt.bh_extraprop_aperture_calc) {
            for (auto i=0; i<opt.aperturenum;i++) {
                for (auto &x:opt.bh_internalprop_output_names_aperture) {
                    headerdatainfo.push_back((string("Aperture_")+x+string("_bh_")+opt.aperture_names_kpc[i]+string("_kpc")));
                }
            }
            for (auto i=0; i<opt.aperturenum;i++) {
                for (auto &x:opt.bh_chem_output_names_aperture) {
                    headerdatainfo.push_back((string("Aperture_")+x+string("_bh_")+opt.aperture_names_kpc[i]+string("_kpc")));
                }
            }
            for (auto i=0; i<opt.aperturenum;i++) {
                for (auto &x:opt.bh_chemproduction_output_names_aperture) {
                    headerdatainfo.push_back((string("Aperture_")+x+string("_bh_")+opt.aperture_names_kpc[i]+string("_kpc")));
                }
            }
#ifdef USEHDF
            sizeval=hdfpredtypeinfo.size();
            for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
            sizeval=adiospredtypeinfo.size();
            for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
        }
#endif
#ifdef EXTRADMON
        if (opt.extra_dm_extraprop_aperture_calc) {
            for (auto i=0; i<opt.aperturenum;i++) {
                for (auto &x:opt.extra_dm_internalprop_output_names_aperture) {
                    headerdatainfo.push_back((string("Aperture_")+x+string("_extra_dm_")+opt.aperture_names_kpc[i]+string("_kpc")));
                }
            }
#ifdef USEHDF
            sizeval=hdfpredtypeinfo.size();
            for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
            sizeval=adiospredtypeinfo.size();
            for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
        }
#endif

    }
    if (opt.iaperturecalc>0 && opt.apertureprojnum>0) {
        for (auto k=0;k<3;k++) {
        string projname = "Projected_aperture_"+to_string(k+1)+"_";
        for (auto i=0; i<opt.apertureprojnum;i++)
            headerdatainfo.push_back(projname+string("mass_")+opt.aperture_proj_names_kpc[i]+string("_kpc"));
#ifdef GASON
        for (auto i=0; i<opt.apertureprojnum;i++)
            headerdatainfo.push_back(projname+string("mass_gas_")+opt.aperture_proj_names_kpc[i]+string("_kpc"));
#ifdef STARON
        for (auto i=0; i<opt.apertureprojnum;i++)
            headerdatainfo.push_back(projname+string("mass_gas_sf_")+opt.aperture_proj_names_kpc[i]+string("_kpc"));
        for (auto i=0; i<opt.apertureprojnum;i++)
            headerdatainfo.push_back(projname+string("mass_gas_nsf_")+opt.aperture_proj_names_kpc[i]+string("_kpc"));
#endif
#endif
#ifdef STARON
        for (auto i=0; i<opt.apertureprojnum;i++)
            headerdatainfo.push_back(projname+string("mass_star_")+opt.aperture_proj_names_kpc[i]+string("_kpc"));
#endif

        for (auto i=0; i<opt.apertureprojnum;i++)
            headerdatainfo.push_back(projname+string("rhalfmass_")+opt.aperture_proj_names_kpc[i]+string("_kpc"));
#ifdef GASON
        for (auto i=0; i<opt.apertureprojnum;i++)
            headerdatainfo.push_back(projname+string("rhalfmass_gas_")+opt.aperture_proj_names_kpc[i]+string("_kpc"));
#ifdef STARON
        for (auto i=0; i<opt.apertureprojnum;i++)
            headerdatainfo.push_back(projname+string("rhalfmass_gas_sf_")+opt.aperture_proj_names_kpc[i]+string("_kpc"));
        for (auto i=0; i<opt.apertureprojnum;i++)
            headerdatainfo.push_back(projname+string("rhalfmass_gas_nsf_")+opt.aperture_proj_names_kpc[i]+string("_kpc"));
#endif
#endif
#ifdef STARON
        for (auto i=0; i<opt.apertureprojnum;i++)
            headerdatainfo.push_back(projname+string("rhalfmass_star_")+opt.aperture_proj_names_kpc[i]+string("_kpc"));
#endif
#if defined(GASON) && defined(STARON)
        for (auto i=0; i<opt.apertureprojnum;i++)
            headerdatainfo.push_back(projname+string("SFR_gas_")+opt.aperture_proj_names_kpc[i]+string("_kpc"));
        for (auto i=0; i<opt.apertureprojnum;i++)
            headerdatainfo.push_back(projname+string("Zmet_gas_")+opt.aperture_proj_names_kpc[i]+string("_kpc"));
        for (auto i=0; i<opt.apertureprojnum;i++)
            headerdatainfo.push_back(projname+string("Zmet_gas_sf_")+opt.aperture_proj_names_kpc[i]+string("_kpc"));
        for (auto i=0; i<opt.apertureprojnum;i++)
            headerdatainfo.push_back(projname+string("Zmet_gas_nsf_")+opt.aperture_proj_names_kpc[i]+string("_kpc"));
        for (auto i=0; i<opt.apertureprojnum;i++)
            headerdatainfo.push_back(projname+string("Zmet_star_")+opt.aperture_proj_names_kpc[i]+string("_kpc"));
#endif
        }
#ifdef USEHDF
        sizeval=hdfpredtypeinfo.size();
        // for (int i=sizeval;i<headerdatainfo.size();i++) predtypeinfo.push_back(desiredproprealtype[0]);
        for (int i=sizeval;i<headerdatainfo.size();i++) hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
        sizeval=adiospredtypeinfo.size();
        for (int i=sizeval;i<headerdatainfo.size();i++) adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
    }

    //if aperture information calculated also include
    if (opt.SOnum>0) {
        for (auto i=0; i<opt.SOnum;i++) {
            headerdatainfo.push_back((string("SO_Mass_")+opt.SOthresholds_names_crit[i]+string("_rhocrit")));
#ifdef USEHDF
            // predtypeinfo.push_back(desiredproprealtype[0]);
            hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
            adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
        }
        for (auto i=0; i<opt.SOnum;i++) {
            headerdatainfo.push_back((string("SO_R_")+opt.SOthresholds_names_crit[i]+string("_rhocrit")));
#ifdef USEHDF
            // predtypeinfo.push_back(desiredproprealtype[0]);
            hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
            adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
        }
#ifdef GASON
        if (opt.iextragasoutput && opt.iextrahalooutput) {
            for (auto i=0; i<opt.SOnum;i++) {
                headerdatainfo.push_back((string("SO_Mass_gas_")+opt.SOthresholds_names_crit[i]+string("_rhocrit")));
#ifdef USEHDF
                // predtypeinfo.push_back(desiredproprealtype[0]);
                hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
                adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
            }
#ifdef STARON
#endif
        }
#endif

#ifdef STARON
        if (opt.iextrastaroutput && opt.iextrahalooutput) {
            for (auto i=0; i<opt.SOnum;i++) {
                headerdatainfo.push_back((string("SO_Mass_star_")+opt.SOthresholds_names_crit[i]+string("_rhocrit")));
#ifdef USEHDF
                // predtypeinfo.push_back(desiredproprealtype[0]);
                hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
                adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
            }
        }
#endif
#ifdef HIGHRES
        if (opt.iextrainterloperoutput && opt.iextrahalooutput) {
            for (auto i=0; i<opt.SOnum;i++) {
                headerdatainfo.push_back((string("SO_Mass_interloper_")+opt.SOthresholds_names_crit[i]+string("_rhocrit")));
#ifdef USEHDF
                // predtypeinfo.push_back(desiredproprealtype[0]);
                hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
                adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
            }
        }
#endif
    }
    if (opt.SOnum>0 && opt.iextrahalooutput) {
        for (auto i=0; i<opt.SOnum;i++) {
            headerdatainfo.push_back((string("SO_Lx_")+opt.SOthresholds_names_crit[i]+string("_rhocrit")));
            headerdatainfo.push_back((string("SO_Ly_")+opt.SOthresholds_names_crit[i]+string("_rhocrit")));
            headerdatainfo.push_back((string("SO_Lz_")+opt.SOthresholds_names_crit[i]+string("_rhocrit")));
            for (auto k=0;k<3;k++) {
#ifdef USEHDF
                // predtypeinfo.push_back(desiredproprealtype[0]);
                hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
                adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
            }
        }
#ifdef GASON
        if (opt.iextragasoutput) {
            for (auto i=0; i<opt.SOnum;i++) {
                headerdatainfo.push_back((string("SO_Lx_gas_")+opt.SOthresholds_names_crit[i]+string("_rhocrit")));
                headerdatainfo.push_back((string("SO_Ly_gas_")+opt.SOthresholds_names_crit[i]+string("_rhocrit")));
                headerdatainfo.push_back((string("SO_Lz_gas_")+opt.SOthresholds_names_crit[i]+string("_rhocrit")));
                for (auto k=0;k<3;k++) {
#ifdef USEHDF
                    // predtypeinfo.push_back(desiredproprealtype[0]);
                    hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
                    adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
                }
            }
#ifdef STARON
#endif
        }
#endif

#ifdef STARON
        if (opt.iextrastaroutput) {
            for (auto i=0; i<opt.SOnum;i++) {
                headerdatainfo.push_back((string("SO_Lx_star_")+opt.SOthresholds_names_crit[i]+string("_rhocrit")));
                headerdatainfo.push_back((string("SO_Ly_star_")+opt.SOthresholds_names_crit[i]+string("_rhocrit")));
                headerdatainfo.push_back((string("SO_Lz_star_")+opt.SOthresholds_names_crit[i]+string("_rhocrit")));
                for (auto k=0;k<3;k++) {
#ifdef USEHDF
                    // predtypeinfo.push_back(desiredproprealtype[0]);
                    hdfpredtypeinfo.push_back(hdfdesiredproprealtype[0]);
#endif
#ifdef USEADIOS
                    adiospredtypeinfo.push_back(desiredadiosproprealtype[0]);
#endif
                }
            }
        }
#endif
    }
}

//@}
