/*! \file allvars.cxx
 *  \brief provides instances of global variables
 */

#include <cassert>
#include <set>

#include "allvars.h"

//-- Structures and external variables
///external pointer to keep track of structure levels and parent
StrucLevelData *psldata;


///\name define routines for the HeaderUnitInfo
HeaderUnitInfo::HeaderUnitInfo(string s)
{
    string delimiter(":"), token;
    vector<string> units;
    size_t pos;
    while ((pos = s.find(delimiter)) != string::npos) {
        token = s.substr(0, pos);
        units.push_back(token);
        s.erase(0, pos + delimiter.length());
    }
    if (units.size()==5) {
        massdim = stof(units[0]);
        lengthdim = stof(units[1]);
        velocitydim = stof(units[2]);
        timedim = stof(units[3]);
	tempdim = stof(units[4]);
        extrainfo = string("");
    }
    else {
        massdim = 0;
        lengthdim = 0;
        velocitydim = 0;
        timedim = 0;
	tempdim = 0;
        extrainfo = s;
    }
}


///\name define write routines for the property data structure
//{@
///write (append) the properties data to an already open binary file
void PropData::WriteBinary(fstream &Fout, Options&opt){
    long long lval;
    long unsigned idval;
    unsigned int ival;
    double val, val3[3],val9[9];
    float fval;
    int sonum_hotgas;

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
    val=gRhalf200m;
    Fout.write((char*)&val,sizeof(val));
    val=gRhalf200c;
    Fout.write((char*)&val,sizeof(val));
    val=gRhalfBN98;
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
    val=cNFW200c;
    Fout.write((char*)&val,sizeof(val));
    val=cNFW200m;
    Fout.write((char*)&val,sizeof(val));
    val=cNFWBN98;
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
#if (defined(GASON)) || (defined(GASON) && defined(SWIFTINTERFACE))
    /*writing M_gas_highT and other related quantities*/
    val=M_gas_highT;
    Fout.write((char*)&val,sizeof(val));
    val=Temp_mean_gas_highT;
    Fout.write((char*)&val,sizeof(val));
    val=Z_mean_gas_highT;
    Fout.write((char*)&val,sizeof(val));
    val=M_gas_highT_incl;
    Fout.write((char*)&val,sizeof(val));
    val=Temp_mean_gas_highT_incl;
    Fout.write((char*)&val,sizeof(val));
    val=Z_mean_gas_highT_incl;
    Fout.write((char*)&val,sizeof(val));

    sonum_hotgas = opt.aperture_hotgas_normalised_to_overdensity.size();
    if (sonum_hotgas>0){
        for (auto j=0;j<sonum_hotgas;j++) {
            Fout.write((char*)&SO_totalmass_highT[j],sizeof(int));
        }
        for (auto j=0;j<sonum_hotgas;j++) {
            Fout.write((char*)&SO_mass_highT[j],sizeof(int));
        }
        for (auto j=0;j<sonum_hotgas;j++) {
            Fout.write((char*)&SO_Temp_mean_gas_highT[j],sizeof(int));
        }
        for (auto j=0;j<sonum_hotgas;j++) {
            Fout.write((char*)&SO_Z_mean_gas_highT[j],sizeof(int));
        }
    }
#endif
    val=M_tot_incl;
    Fout.write((char*)&val,sizeof(val));
#ifdef GASON
    val=M_gas_incl;
    Fout.write((char*)&val,sizeof(val));
#ifdef STARON
    val=M_gas_nsf_incl;
    Fout.write((char*)&val,sizeof(val));
    val=M_gas_sf_incl;
    Fout.write((char*)&val,sizeof(val));
#endif
#endif
#if STARON
    val=M_star_incl;
    Fout.write((char*)&val,sizeof(val));
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
	#if (defined(GASON)) || (defined(GASON) && defined(SWIFTINTERFACE))
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_M_gas_highT[j],sizeof(val));
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_Temp_mean_gas_highT[j],sizeof(val));
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout.write((char*)&aperture_Z_mean_gas_highT[j],sizeof(val));
        }

        if (opt.SOnum>0){
            for (auto j=0;j<opt.SOnum;j++) {
                Fout<<SO_mass[j]<<" ";
            }
            for (auto j=0;j<opt.SOnum;j++) {
                Fout<<SO_radius[j]<<" ";
            }
            for (auto j=0;j<opt.SOnum;j++) {
                Fout<<SO_mass_gas[j]<<" ";
            }
            for (auto j=0;j<opt.SOnum;j++) {
                Fout<<SO_mass_star[j]<<" ";
            }
	    #ifdef HIGHRES
            for (auto j=0;j<opt.SOnum;j++) {
                Fout<<SO_mass_interloper[j]<<" ";
            }
	    #endif
        }
	#endif
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
    int sonum_hotgas;
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
    Fout<<gRhalf200m<<" ";
    Fout<<gRhalf200c<<" ";
    Fout<<gRhalfBN98<<" ";
    Fout<<gmaxvel<<" ";
    Fout<<gsigma_v<<" ";
    for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<gveldisp(k,n)<<" ";
    Fout<<glambda_B<<" ";
    for (int k=0;k<3;k++) Fout<<gJ[k]<<" ";
    Fout<<gq<<" ";
    Fout<<gs<<" ";
    for (int k=0;k<3;k++) for (int n=0;n<3;n++) Fout<<geigvec(k,n)<<" ";
    Fout<<cNFW<<" ";
    Fout<<cNFW200c<<" ";
    Fout<<cNFW200m<<" ";
    Fout<<cNFWBN98<<" ";
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
#if (defined(GASON)) || (defined(GASON) && defined(SWIFTINTERFACE))
    Fout<<M_gas_highT<<" ";
    Fout<<Temp_mean_gas_highT<<" ";
    Fout<<Z_mean_gas_highT<<" ";
    Fout<<M_gas_highT_incl<<" ";
    Fout<<Temp_mean_gas_highT_incl<<" ";
    Fout<<Z_mean_gas_highT_incl<<" ";

    sonum_hotgas = opt.aperture_hotgas_normalised_to_overdensity.size();
    if (sonum_hotgas>0){
        for (auto j=0;j<sonum_hotgas;j++) {
            Fout.write((char*)&SO_totalmass_highT[j],sizeof(int));
        }
        for (auto j=0;j<sonum_hotgas;j++) {
            Fout.write((char*)&SO_mass_highT[j],sizeof(int));
        }
        for (auto j=0;j<sonum_hotgas;j++) {
            Fout.write((char*)&SO_Temp_mean_gas_highT[j],sizeof(int));
        }
        for (auto j=0;j<sonum_hotgas;j++) {
            Fout.write((char*)&SO_Z_mean_gas_highT[j],sizeof(int));
        }
    }

#endif
    Fout<<M_tot_incl<<" ";
#ifdef GASON
    Fout<<M_gas_incl<<" ";
#ifdef STARON
    Fout<<M_gas_nsf_incl<<" ";
    Fout<<M_gas_sf_incl<<" ";
#endif
#endif
#ifdef STARON
    Fout<<M_star_incl<<" ";
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
        #if (defined(GASON)) || (defined(GASON) && defined(SWIFTINTERFACE))
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_M_gas_highT[j]<<" ";
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_Temp_mean_gas_highT[j]<<" ";
        }
        for (auto j=0;j<opt.aperturenum;j++) {
            Fout<<aperture_Z_mean_gas_highT[j]<<" ";
        }
        #endif
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

template <typename T>
struct output_traits;

#ifdef USEHDF
#ifdef USEADIOS
#define OUTPUT_TRAITS(CTYPE, HDF5_TYPE, ADIOS_TYPE)                     \
    template <>                                                         \
    struct output_traits<CTYPE>                                         \
    {                                                                   \
        const static hid_t hdf5_type;                                   \
        const static ADIOS_DATATYPES adios_type;                        \
    };                                                                  \
    const hid_t output_traits<CTYPE>::hdf5_type = HDF5_TYPE;            \
    const ADIOS_DATATYPES output_traits<CTYPE>::adios_type = ADIOS_TYPE
#else
#define OUTPUT_TRAITS(CTYPE, HDF5_TYPE, ADIOS_TYPE)         \
    template <>                                             \
    struct output_traits<CTYPE>                             \
    {                                                       \
        const static hid_t hdf5_type;                       \
    };                                                      \
    const hid_t output_traits<CTYPE>::hdf5_type = HDF5_TYPE
#endif // USEADIOS
#endif // USEHDF

OUTPUT_TRAITS(int, H5T_NATIVE_INT, adios_integer);
OUTPUT_TRAITS(unsigned int, H5T_NATIVE_UINT, adios_unsigned_integer);
OUTPUT_TRAITS(long, H5T_NATIVE_LONG, adios_long);
OUTPUT_TRAITS(unsigned long, H5T_NATIVE_ULONG, adios_unsigned_long);
OUTPUT_TRAITS(float, H5T_NATIVE_FLOAT, adios_real);
OUTPUT_TRAITS(double, H5T_NATIVE_DOUBLE, adios_double);


const HeaderUnitInfo PropDataHeader::UNITLESS;
const HeaderUnitInfo PropDataHeader::MASS {1};
const HeaderUnitInfo PropDataHeader::TEMPERATURE {4};
const HeaderUnitInfo PropDataHeader::LENGTH {0, 1};
const HeaderUnitInfo PropDataHeader::VELOCITY {0, 0, 1};
const HeaderUnitInfo PropDataHeader::VELOCITY_2D {0, 0, 2};
const HeaderUnitInfo PropDataHeader::ANGULAR_MOMENTUM {1, 1, 1};
const HeaderUnitInfo PropDataHeader::ENERGY {1, 0, 2, 0};
const HeaderUnitInfo PropDataHeader::MASS_OVER_TIME {1, 0, 0, -1};
const HeaderUnitInfo PropDataHeader::TIME {0, 0, 0, 1};

template <typename T>
void PropDataHeader::declare_dataset(std::string name, HeaderUnitInfo unit_info)
{
    headerdatainfo.emplace_back(std::move(name));
    unitdatainfo.emplace_back(std::move(unit_info));
#ifdef USEHDF
    hdfpredtypeinfo.emplace_back(output_traits<T>::hdf5_type);
#endif // USEHDF
#ifdef USEADIOS
    adiospredtypeinfo.push_back(output_traits<T>::adios_type);
#endif // USEADIOS
}

void PropDataHeader::declare_dataset(std::string name, HeaderUnitInfo unit_info)
{
	declare_dataset<Double_t>(std::move(name), std::move(unit_info));
}

void PropDataHeader::declare_xyz_datasets(std::string prefix, HeaderUnitInfo unit_info)
{
	declare_dataset(prefix + "x", unit_info);
	declare_dataset(prefix + "y", unit_info);
	declare_dataset(prefix + "z", unit_info);
}

void PropDataHeader::declare_xyz_datasets(std::string prefix, std::string suffix, HeaderUnitInfo unit_info)
{
	declare_dataset(prefix + "x" + suffix, unit_info);
	declare_dataset(prefix + "y" + suffix, unit_info);
	declare_dataset(prefix + "z" + suffix, unit_info);
}

void PropDataHeader::declare_XYZ_datasets(std::string suffix, HeaderUnitInfo unit_info)
{
	declare_dataset("X" + suffix, unit_info);
	declare_dataset("Y" + suffix, unit_info);
	declare_dataset("Z" + suffix, unit_info);
}

void PropDataHeader::declare_XYZ_datasets(std::string prefix, std::string suffix, HeaderUnitInfo unit_info)
{
	declare_dataset(prefix + "X" + suffix, unit_info);
	declare_dataset(prefix + "Y" + suffix, unit_info);
	declare_dataset(prefix + "Z" + suffix, unit_info);
}

void PropDataHeader::declare_xyz2_datasets(std::string prefix, HeaderUnitInfo unit_info)
{
	declare_xyz_datasets(prefix + "x", unit_info);
	declare_xyz_datasets(prefix + "y", unit_info);
	declare_xyz_datasets(prefix + "z", unit_info);
}

void PropDataHeader::declare_xyz2_datasets(std::string prefix, std::string suffix, HeaderUnitInfo unit_info)
{
	declare_xyz_datasets(prefix + "x", suffix, unit_info);
	declare_xyz_datasets(prefix + "y", suffix, unit_info);
	declare_xyz_datasets(prefix + "z", suffix, unit_info);
}

void PropDataHeader::declare_datasets(const std::vector<string> &names, const std::vector<string> &units, std::string suffix)
{
	assert(names.size() == units.size());
	if (!suffix.empty()) {
		suffix = "_" + suffix;
	}
	for (std::size_t i = 0; i != names.size(); i++) {
		declare_dataset(names[i] + suffix, HeaderUnitInfo{units[i]});
	}
}

PropDataHeader::PropDataHeader(const Options &opt)
{
	declare_all_datasets(opt);
	// double-check there are no duplicate dataset names and that all fields
	// are the same size
	assert(std::set<std::string>(headerdatainfo.begin(), headerdatainfo.end()).size() == headerdatainfo.size());
	assert(headerdatainfo.size() == unitdatainfo.size());
#ifdef USEHDF
    assert(headerdatainfo.size() == hdfpredtypeinfo.size());
#endif // USEHDF
#ifdef USEADIOS
    assert(headerdatainfo.size() == adiospredtypeinfo.size());
#endif // USEADIOS
}

void PropDataHeader::declare_all_datasets(const Options &opt)
{
    declare_dataset<unsigned long>("ID");
    declare_dataset<long>("ID_mbp");
    declare_dataset<long>("ID_minpot");
    declare_dataset<long>("hostHaloID");
    declare_dataset<unsigned long>("numSubStruct");
    declare_dataset<unsigned long>("npart");
    declare_dataset<int>("Structuretype");
    if (opt.iKeepFOF==1){
        declare_dataset<long>("hostDirectHaloID");
        declare_dataset<long>("hostFOFID");
    }

    declare_dataset("Mvir", MASS);
    declare_XYZ_datasets("c", LENGTH);
    declare_XYZ_datasets("cmbp", LENGTH);
    declare_XYZ_datasets("cminpot", LENGTH);
    declare_XYZ_datasets("V", "c", VELOCITY);
    declare_XYZ_datasets("V", "cmbp", VELOCITY);
    declare_XYZ_datasets("V", "cminpot", VELOCITY);
    declare_dataset("Mass_tot", MASS);
    declare_dataset("Mass_FOF", MASS);
    declare_dataset("Mass_200mean", MASS);
    declare_dataset("Mass_200crit", MASS);
    declare_dataset("Mass_BN98", MASS);
    declare_dataset("Efrac");
    declare_dataset("Rvir", LENGTH);
    declare_dataset("R_size", LENGTH);
    declare_dataset("R_200mean", LENGTH);
    declare_dataset("R_200crit", LENGTH);
    declare_dataset("R_BN98", LENGTH);
    declare_dataset("R_HalfMass", LENGTH);
    declare_dataset("Rmax", LENGTH);
    declare_dataset("R_HalfMass_200mean", LENGTH);
    declare_dataset("R_HalfMass_200crit", LENGTH);
    declare_dataset("R_HalfMass_BN98", LENGTH);
    declare_dataset("Vmax", VELOCITY);
    declare_dataset("sigV", VELOCITY);
    declare_xyz2_datasets("veldisp_", VELOCITY_2D);
    declare_dataset("lambda_B");
    declare_xyz_datasets("L", ANGULAR_MOMENTUM);
    declare_dataset("q");
    declare_dataset("s");
    declare_xyz2_datasets("eig_");
    declare_dataset("cNFW");
    declare_dataset("cNFW_200crit");
    declare_dataset("cNFW_200mean");
    declare_dataset("cNFW_BN98");
    declare_dataset("Krot");
    declare_dataset("Ekin", ENERGY);
    declare_dataset("Epot", ENERGY);

    //some properties within RVmax
    declare_dataset("RVmax_sigV", VELOCITY);
    declare_xyz2_datasets("RVmax_veldisp_", VELOCITY_2D);
    declare_dataset("RVmax_lambda_B");
    declare_xyz_datasets("RVmax_L", ANGULAR_MOMENTUM);
    declare_dataset("RVmax_q");
    declare_dataset("RVmax_s");
    declare_xyz2_datasets("RVmax_eig_");

    if (opt.iextrahalooutput) {
        declare_xyz_datasets("L", "_200mean", ANGULAR_MOMENTUM);
        declare_xyz_datasets("L", "_200crit", ANGULAR_MOMENTUM);
        declare_xyz_datasets("L", "_BN98", ANGULAR_MOMENTUM);
        if (opt.iInclusiveHalo > 0) {
            declare_dataset("Mass_200mean_excl", MASS);
            declare_dataset("Mass_200crit_excl", MASS);
            declare_dataset("Mass_BN98_excl", MASS);
            declare_dataset("R_200mean_excl", LENGTH);
            declare_dataset("R_200crit_excl", LENGTH);
            declare_dataset("R_BN98_excl", LENGTH);
            declare_xyz_datasets("L", "_200mean_excl", ANGULAR_MOMENTUM);
            declare_xyz_datasets("L", "_200crit_excl", ANGULAR_MOMENTUM);
            declare_xyz_datasets("L", "_BN98_excl", ANGULAR_MOMENTUM);
        }
    }

#ifdef GASON
    declare_dataset<unsigned long>("n_gas");
    declare_dataset("Mass_gas", MASS);
    declare_dataset("Mass_gas_Rvmax", MASS);
    declare_dataset("Mass_gas_30kpc", MASS);
    declare_dataset("Mass_gas_500c", MASS);
    declare_XYZ_datasets("c_gas", LENGTH);
    declare_XYZ_datasets("V", "c_gas", VELOCITY);
    declare_dataset("Efrac_gas");
    declare_dataset("R_HalfMass_gas", LENGTH);
    declare_xyz2_datasets("veldisp_", "_gas", VELOCITY_2D);
    declare_xyz_datasets("L", "_gas", ANGULAR_MOMENTUM);
    declare_dataset("q_gas");
    declare_dataset("s_gas");
    declare_xyz2_datasets("eig_", "_gas");
    declare_dataset("Krot_gas");
    declare_dataset("T_gas", TEMPERATURE);
#ifdef STARON
    declare_dataset("Zmet_gas");
    declare_dataset("SFR_gas", MASS_OVER_TIME);
#endif
    if (opt.iextragasoutput) {
        declare_dataset("Mass_200mean_gas", MASS);
        declare_dataset("Mass_200crit_gas", MASS);
        declare_dataset("Mass_BN98_gas", MASS);
        declare_xyz_datasets("L", "_200c_gas", ANGULAR_MOMENTUM);
        declare_xyz_datasets("L", "_200m_gas", ANGULAR_MOMENTUM);
        declare_xyz_datasets("L", "_BN98_gas", ANGULAR_MOMENTUM);
        if (opt.iInclusiveHalo>0) {
            declare_dataset("Mass_200mean_excl_gas", MASS);
            declare_dataset("Mass_200crit_excl_gas", MASS);
            declare_dataset("Mass_BN98_excl_gas", MASS);
            declare_xyz_datasets("L", "_200c_excl_gas", ANGULAR_MOMENTUM);
            declare_xyz_datasets("L", "_200m_excl_gas", ANGULAR_MOMENTUM);
            declare_xyz_datasets("L", "_BN98_excl_gas", ANGULAR_MOMENTUM);
        }
    }
#endif // GASON

#ifdef STARON
    declare_dataset<unsigned long>("n_star");
    declare_dataset("Mass_star", MASS);
    declare_dataset("Mass_star_Rvmax", MASS);
    declare_dataset("Mass_star_30kpc", MASS);
    //declare_dataset("M_star_50kpc", MASS);
    declare_dataset("Mass_star_500c", MASS);
    declare_XYZ_datasets("c_star", LENGTH);
    declare_XYZ_datasets("V", "c_star", VELOCITY);
    declare_dataset("Efrac_star");
    declare_dataset("R_HalfMass_star", LENGTH);
    declare_xyz2_datasets("veldisp_", "_star", VELOCITY_2D);
    declare_xyz_datasets("L", "_star", ANGULAR_MOMENTUM);
    declare_dataset("q_star");
    declare_dataset("s_star");
    declare_xyz2_datasets("eig_", "_star");
    declare_dataset("Krot_star");
    declare_dataset("tage_star", TIME);
    declare_dataset("Zmet_star");
    if (opt.iextrastaroutput) {
        declare_dataset("Mass_200mean_star", MASS);
        declare_dataset("Mass_200crit_star", MASS);
        declare_dataset("Mass_BN98_star", MASS);
        declare_xyz_datasets("L", "_200c_star", ANGULAR_MOMENTUM);
        declare_xyz_datasets("L", "_200m_star", ANGULAR_MOMENTUM);
        declare_xyz_datasets("L", "_BN98_star", ANGULAR_MOMENTUM);
        if (opt.iInclusiveHalo>0) {
            declare_dataset("Mass_200mean_excl_star", MASS);
            declare_dataset("Mass_200crit_excl_star", MASS);
            declare_dataset("Mass_BN98_excl_star", MASS);
            declare_xyz_datasets("L", "_200c_excl_star", ANGULAR_MOMENTUM);
            declare_xyz_datasets("L", "_200m_excl_star", ANGULAR_MOMENTUM);
            declare_xyz_datasets("L", "_BN98_excl_star", ANGULAR_MOMENTUM);
        }
    }
#endif // STARON

#ifdef BHON
    declare_dataset<unsigned long>("n_bh");
    declare_dataset("Mass_bh", MASS);
#endif // BHON

#ifdef HIGHRES
    declare_dataset<unsigned long>("n_interloper");
    declare_dataset("Mass_interloper", MASS);
    if (opt.iextrainterloperoutput) {
        declare_dataset("Mass_200mean_interloper", MASS);
        declare_dataset("Mass_200crit_interloper", MASS);
        declare_dataset("Mass_BN98_interloper", MASS);
        if (opt.iInclusiveHalo>0) {
            declare_dataset("Mass_200mean_excl_interloper", MASS);
            declare_dataset("Mass_200crit_excl_interloper", MASS);
            declare_dataset("Mass_BN98_excl_interloper", MASS);
        }
    }
#endif // HIGHRES

#if defined(GASON) && defined(STARON)
    declare_dataset("Mass_gas_sf", MASS);
    declare_dataset("R_HalfMass_gas_sf", LENGTH);
    declare_dataset("sigV_gas_sf", VELOCITY);
    declare_xyz_datasets("L", "_gas_sf", ANGULAR_MOMENTUM);
    declare_dataset("Krot_gas_sf");
    declare_dataset("T_gas_sf", TEMPERATURE);
    declare_dataset("Zmet_gas_sf");
    if (opt.iextragasoutput) {
        declare_dataset("Mass_200mean_gas_sf", MASS);
        declare_dataset("Mass_200crit_gas_sf", MASS);
        declare_dataset("Mass_BN98_gas_sf", MASS);
        declare_xyz_datasets("L", "_200c_gas_sf", ANGULAR_MOMENTUM);
        declare_xyz_datasets("L", "_200m_gas_sf", ANGULAR_MOMENTUM);
        declare_xyz_datasets("L", "_BN98_gas_sf", ANGULAR_MOMENTUM);
        if (opt.iInclusiveHalo>0) {
            declare_dataset("Mass_200mean_excl_gas_sf", MASS);
            declare_dataset("Mass_200crit_excl_gas_sf", MASS);
            declare_dataset("Mass_BN98_excl_gas_sf", MASS);
            declare_xyz_datasets("L", "_200c_excl_gas_sf", ANGULAR_MOMENTUM);
            declare_xyz_datasets("L", "_200m_excl_gas_sf", ANGULAR_MOMENTUM);
            declare_xyz_datasets("L", "_BN98_excl_gas_sf", ANGULAR_MOMENTUM);
        }
    }
    declare_dataset("Mass_gas_nsf", MASS);
    declare_dataset("R_HalfMass_gas_nsf", LENGTH);
    declare_dataset("sigV_gas_nsf", VELOCITY);
    declare_xyz_datasets("L", "_gas_nsf", ANGULAR_MOMENTUM);
    declare_dataset("Krot_gas_nsf");
    declare_dataset("T_gas_nsf", TEMPERATURE);
    declare_dataset("Zmet_gas_nsf");

    if (opt.iextragasoutput) {
        declare_dataset("Mass_200mean_gas_nsf", MASS);
        declare_dataset("Mass_200crit_gas_nsf", MASS);
        declare_dataset("Mass_BN98_gas_nsf", MASS);
        declare_xyz_datasets("L", "_200c_gas_nsf", ANGULAR_MOMENTUM);
        declare_xyz_datasets("L", "_200m_gas_nsf", ANGULAR_MOMENTUM);
        declare_xyz_datasets("L", "_BN98_gas_nsf", ANGULAR_MOMENTUM);
        if (opt.iInclusiveHalo>0) {
            declare_dataset("Mass_200mean_excl_gas_nsf", MASS);
            declare_dataset("Mass_200crit_excl_gas_nsf", MASS);
            declare_dataset("Mass_BN98_excl_gas_nsf", MASS);
            declare_xyz_datasets("L", "_200c_excl_gas_nsf", ANGULAR_MOMENTUM);
            declare_xyz_datasets("L", "_200m_excl_gas_nsf", ANGULAR_MOMENTUM);
            declare_xyz_datasets("L", "_BN98_excl_gas_nsf", ANGULAR_MOMENTUM);
        }
    }
#endif // defined(GASON) && defined(STARON)

#if (defined(GASON)) || (defined(GASON) && defined(SWIFTINTERFACE))
    declare_dataset("Mass_gas_highT_excl", MASS);
    declare_dataset("T_gas_highT_excl", TEMPERATURE);
    declare_dataset("Zmet_gas_highT_excl");
    declare_dataset("Mass_gas_highT_incl", MASS);
    declare_dataset("T_gas_highT_incl", TEMPERATURE);
    declare_dataset("Zmet_gas_highT_incl");

    int sonum_hotgas = opt.aperture_hotgas_normalised_to_overdensity.size();
    if (sonum_hotgas > 0) {
        for (auto aperture: opt.aperture_hotgas_normalised_to_overdensity) {
            std::ostringstream name;
            name << "SO_Mass_highT_" << aperture << "_times_"
                 << opt.hot_gas_overdensity_normalisation << "_rhocrit";
            declare_dataset(name.str(), MASS);
        }
        for (auto aperture: opt.aperture_hotgas_normalised_to_overdensity) {
            std::ostringstream name;
            name << "SO_Mass_gas_highT_" << aperture << "_times_"
                 << opt.hot_gas_overdensity_normalisation << "_rhocrit";
            declare_dataset(name.str(), MASS);
        }
        for (auto aperture: opt.aperture_hotgas_normalised_to_overdensity) {
            std::ostringstream name;
            name << "SO_T_gas_highT_" << aperture << "_times_"
                 << opt.hot_gas_overdensity_normalisation << "_rhocrit";
            declare_dataset(name.str(), TEMPERATURE);
        }
        for (auto aperture: opt.aperture_hotgas_normalised_to_overdensity) {
            std::ostringstream name;
            name << "SO_Zmet_gas_highT_" << aperture << "_times_"
                 << opt.hot_gas_overdensity_normalisation << "_rhocrit";
            declare_dataset(name.str());
        }
    }
#endif // (defined(GASON)) || (defined(GASON) && defined(SWIFTINTERFACE))

    declare_dataset("Mass_tot_incl", MASS);

#ifdef GASON
    declare_dataset("Mass_gas_incl", MASS);
#ifdef STARON
    declare_dataset("Mass_gas_nsf_incl", MASS);
    declare_dataset("Mass_gas_sf_incl", MASS);
#endif // STARON
#endif // GASON
#ifdef STARON
    declare_dataset("Mass_star_incl", MASS);
#endif // STARON

    //if extra hydro properties are calculated
#ifdef GASON
    declare_datasets(opt.gas_internalprop_output_names, opt.gas_internalprop_output_units, "gas");
    declare_datasets(opt.gas_chem_output_names, opt.gas_chem_output_units, "gas");
    declare_datasets(opt.gas_chemproduction_output_names, opt.gas_chemproduction_output_units, "gas");
#endif // GASON
#ifdef STARON
    declare_datasets(opt.star_internalprop_output_names, opt.star_internalprop_output_units, "star");
    declare_datasets(opt.star_chem_output_names, opt.star_chem_output_units, "star");
    declare_datasets(opt.star_chemproduction_output_names, opt.star_chemproduction_output_units, "star");
#endif // STARON
#ifdef BHON
    declare_datasets(opt.bh_internalprop_output_names, opt.bh_internalprop_output_units, "bh");
    declare_datasets(opt.bh_chem_output_names, opt.bh_chem_output_units, "bh");
    declare_datasets(opt.bh_chemproduction_output_names, opt.bh_chemproduction_output_units, "bh");
#endif // BHON
#ifdef EXTRADMON
    declare_datasets(opt.extra_dm_internalprop_output_names, opt.extra_dm_internalprop_output_units, "extra_dm");
#endif // EXTRADMON

    //if aperture information calculated also include
    if (opt.iaperturecalc>0 && opt.aperturenum>0) {
        // declare are the unsigned integer aperture fields 
        auto declare_aperture_counts_datasets = [&](const std::string &name)
        {
            for (auto i = 0; i < opt.aperturenum; i++)
                declare_dataset<unsigned int>("Aperture_" + name + "_" + opt.aperture_names_kpc[i] + "_kpc");
        };

        declare_aperture_counts_datasets("npart");
#ifdef GASON
        declare_aperture_counts_datasets("npart_gas");
#ifdef STARON
        declare_aperture_counts_datasets("npart_gas_sf");
        declare_aperture_counts_datasets("npart_gas_nsf");
#endif // STARON
#endif // GASON
#ifdef STARON
        declare_aperture_counts_datasets("npart_star");
#endif // STARON
#ifdef BHON
        declare_aperture_counts_datasets("npart_bh");
#endif // BHON
#ifdef HIGHRES
        declare_aperture_counts_datasets("npart_interloper");
#endif // HIGHRES

        // declare double aperture fields 
        auto declare_aperture_datasets = [&](const std::string &name, HeaderUnitInfo unit_info=UNITLESS)
        {
            for (auto i = 0; i < opt.aperturenum; i++)
                declare_dataset("Aperture_" + name + "_" + opt.aperture_names_kpc[i] + "_kpc", unit_info);
        };

        declare_aperture_datasets("mass", MASS);
#ifdef GASON
        declare_aperture_datasets("mass_gas", MASS);
#ifdef STARON
        declare_aperture_datasets("mass_gas_sf", MASS);
        declare_aperture_datasets("mass_gas_nsf", MASS);
#endif // STARON
#endif // GASON
#ifdef STARON
        declare_aperture_datasets("mass_star", MASS);
#endif // STARON
#ifdef BHON
        declare_aperture_datasets("mass_bh", MASS);
#endif // BHON
#ifdef HIGHRES
        declare_aperture_datasets("mass_interloper", MASS);
#endif

        declare_aperture_datasets("rhalfmass", LENGTH);
#ifdef GASON
        declare_aperture_datasets("rhalfmass_gas", LENGTH);
#ifdef STARON
        declare_aperture_datasets("rhalfmass_gas_sf", LENGTH);
        declare_aperture_datasets("rhalfmass_gas_nsf", LENGTH);
#endif // STARON
#endif // GASON
#ifdef STARON
        declare_aperture_datasets("rhalfmass_star", LENGTH);
#endif // STARON

        declare_aperture_datasets("veldisp", VELOCITY_2D);
#ifdef GASON
        declare_aperture_datasets("veldisp_gas", VELOCITY_2D);
#ifdef STARON
        declare_aperture_datasets("veldisp_gas_sf", VELOCITY_2D);
        declare_aperture_datasets("veldisp_gas_nsf", VELOCITY_2D);
#endif // STARON
#endif // GASON
#ifdef STARON
        declare_aperture_datasets("veldisp_star", VELOCITY_2D);
#endif // STARON

#if defined(GASON) && defined(STARON)
        declare_aperture_datasets("SFR_gas", MASS_OVER_TIME);
        declare_aperture_datasets("Zmet_gas");
        declare_aperture_datasets("Zmet_gas_sf");
        declare_aperture_datasets("Zmet_gas_nsf");
        declare_aperture_datasets("Zmet_star");
#if (defined(GASON)) || (defined(GASON) && defined(SWIFTINTERFACE))
        declare_aperture_datasets("mass_highT", MASS);
        declare_aperture_datasets("T_highT", TEMPERATURE);
        declare_aperture_datasets("Zmet_highT");
#endif // (defined(GASON)) || (defined(GASON) && defined(SWIFTINTERFACE))
#endif // defined(GASON) && defined(STARON)

        //add extra property data in apertures
        auto declare_extra_properties_aperture_datasets = [&](const std::vector<std::string> &names, const std::vector<std::string> &units, const std::string &suffix)
        {
            assert(names.size() == units.size());
            std::vector<std::string> ds_names(names.size());
            for (auto i = 0; i < opt.aperturenum; i++) {
                std::transform(names.begin(), names.end(), ds_names.begin(), [&](const std::string &name)
                {
                    return "Aperture_" + name + "_" + suffix + "_" + opt.aperture_names_kpc[i] + "_kpc";
                });
                declare_datasets(ds_names, units);
            }
        };
#ifdef GASON
        if (opt.gas_extraprop_aperture_calc) {
            declare_extra_properties_aperture_datasets(opt.gas_internalprop_output_names_aperture, opt.gas_internalprop_output_units_aperture, "gas");
            declare_extra_properties_aperture_datasets(opt.gas_chem_output_names_aperture, opt.gas_chem_output_units_aperture, "gas");
            declare_extra_properties_aperture_datasets(opt.gas_chemproduction_output_names_aperture, opt.gas_chemproduction_output_units_aperture, "gas");
        }
#endif
#ifdef STARON
        if (opt.star_extraprop_aperture_calc) {
            declare_extra_properties_aperture_datasets(opt.star_internalprop_output_names_aperture, opt.star_internalprop_output_units_aperture, "star");
            declare_extra_properties_aperture_datasets(opt.star_chem_output_names_aperture, opt.star_chem_output_units_aperture, "star");
            declare_extra_properties_aperture_datasets(opt.star_chemproduction_output_names_aperture, opt.star_chemproduction_output_units_aperture, "star");
        }
#endif
#ifdef BHON
        if (opt.bh_extraprop_aperture_calc) {
            declare_extra_properties_aperture_datasets(opt.bh_internalprop_output_names_aperture, opt.bh_internalprop_output_units_aperture, "bh");
            declare_extra_properties_aperture_datasets(opt.bh_chem_output_names_aperture, opt.bh_chem_output_units_aperture, "bh");
            declare_extra_properties_aperture_datasets(opt.bh_chemproduction_output_names_aperture, opt.bh_chemproduction_output_units_aperture, "bh");
        }
#endif
#ifdef EXTRADMON
        if (opt.extra_dm_extraprop_aperture_calc) {
            declare_extra_properties_aperture_datasets(opt.extra_dm_internalprop_output_names_aperture, opt.extra_dm_internalprop_output_units_aperture, "extra_dm");
        }
#endif
    }

    if (opt.iaperturecalc>0 && opt.apertureprojnum>0) {
        for (auto k=0;k<3;k++) {
            auto declare_projected_aperture_datasets = [&](const std::string &name, HeaderUnitInfo unit_info=UNITLESS)
            {
                std::string prefix = "Projected_aperture_" + std::to_string(k + 1) + "_" + name + "_";
                for (std::size_t i = 0; i < opt.apertureprojnum; i++)
                    declare_dataset(prefix + opt.aperture_proj_names_kpc[i] + "_kpc", unit_info);
            };

            declare_projected_aperture_datasets("mass", MASS);
#ifdef GASON
            declare_projected_aperture_datasets("mass_gas", MASS);
#ifdef STARON
            declare_projected_aperture_datasets("mass_gas_sf", MASS);
            declare_projected_aperture_datasets("mass_gas_nsf", MASS);
#endif // STARON
#endif // GASON
#ifdef STARON
            declare_projected_aperture_datasets("mass_star", MASS);
#endif // STARON

            declare_projected_aperture_datasets("rhalfmass", LENGTH);
#ifdef GASON
            declare_projected_aperture_datasets("rhalfmass_gas", LENGTH);
#ifdef STARON
            declare_projected_aperture_datasets("rhalfmass_gas_sf", LENGTH);
            declare_projected_aperture_datasets("rhalfmass_gas_nsf", LENGTH);
#endif // STARON
#endif // GASON
#ifdef STARON
            declare_projected_aperture_datasets("rhalfmass_star", LENGTH);
#endif // STARON

#if defined(GASON) && defined(STARON)
            declare_projected_aperture_datasets("SFR_gas", MASS_OVER_TIME);
            declare_projected_aperture_datasets("Zmet_gas");
            declare_projected_aperture_datasets("Zmet_gas_sf");
            declare_projected_aperture_datasets("Zmet_gas_nsf");
            declare_projected_aperture_datasets("Zmet_star");
#endif // defined(GASON) && defined(STARON)
        }
    }

    // Spherical Overdensities
    auto declare_so_datasets = [&](const std::string &name, HeaderUnitInfo unit_info=UNITLESS)
    {
        for (auto i = 0; i < opt.SOnum; i++) {
            declare_dataset("SO_" + name + "_" + opt.SOthresholds_names_crit[i] + "_rhocrit", unit_info);
        }
    };
    auto declare_so_angular_momentum_datasets = [&](std::string name="")
    {
        if (!name.empty()) {
            name = "_" + name;
        }
        for (auto i = 0; i < opt.SOnum; i++) {
            declare_xyz_datasets("SO_L", name + "_" + opt.SOthresholds_names_crit[i] + "_rhocrit", ANGULAR_MOMENTUM);
        }
    };
    if (opt.SOnum>0) {
        declare_so_datasets("Mass", MASS);
        declare_so_datasets("R", LENGTH);
#ifdef GASON
        if (opt.iextragasoutput && opt.iextrahalooutput) {
            declare_so_datasets("Mass_gas", MASS);
        }
#endif // GASON
#ifdef STARON
        if (opt.iextrastaroutput && opt.iextrahalooutput) {
            declare_so_datasets("Mass_star", MASS);
        }
#endif // STARON
#ifdef HIGHRES
        if (opt.iextrainterloperoutput && opt.iextrahalooutput) {
            declare_so_datasets("Mass_interloper", MASS);
        }
#endif // HIGHRES
    }
    if (opt.SOnum>0 && opt.iextrahalooutput) {
        declare_so_angular_momentum_datasets();
#ifdef GASON
        if (opt.iextragasoutput) {
            declare_so_angular_momentum_datasets("gas");
        }
#endif // GASON
#ifdef STARON
        if (opt.iextrastaroutput) {
            declare_so_angular_momentum_datasets("star");
        }
#endif // STARON
    }
}

//@}
