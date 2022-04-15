import numpy as np
from copy import deepcopy as dc
from scipy.interpolate import griddata
from scipy.optimize import fsolve
from UnitConversions import UnitConverter as UC

class saturatedData():
    """
    I made this class for storing the saturated water table data, so I only have to read it once.
    """
    def __init__(self):
        self.readTable()
        self.Pc=max(self.PCol)
        self.Tc = float(griddata(self.PCol, self.TCol, self.Pc,method='cubic'))
        self.hc = float(griddata(self.PCol, self.hfCol, self.Pc,method='cubic'))
        self.sc = float(griddata(self.PCol, self.sfCol, self.Pc,method='cubic'))
        self.vc = float(griddata(self.PCol, self.vfCol, self.Pc,method='cubic'))

    def readTable(self):
        """
        I expect a tuple of columns of data in the order:
        ts, ps, hfs, hgs, sfs, sgs, vfs, vgs
        :param vals:
        :return:
        """
        self.TCol, self.PCol, self.hfCol, self.hgCol, self.sfCol, self.sgCol, self.vfCol, self.vgCol = np.loadtxt('sat_water_table.txt', skiprows=1, unpack=True)

    def getSatProps(self, P_Bar=None, P_kPa=None, T=None, vg=None, hg=None, sg=None, sf=None, method='cubic'):
        """
        Retrieve an object of satProps()
        :param P_Bar: pressure in bars
        :param P_kPA: pressure in kPa
        :param T: should be in C
        :return: satProps() object.  Note: P gets returned in kPa
        """
        # given PSat, calc props
        if P_Bar is not None:
            P=P_Bar
        elif P_kPa is not None:
            P=P_kPa/100.0

        sats = satProps()
        if vg is not None:
            sats.Psat = float(griddata(self.vgCol, self.PCol, vg, method=method))
            P=sats.Psat
        elif hg is not None:
            sats.Psat = float(griddata(self.hgCol, self.PCol, hg, method=method))
            P=sats.Psat
        elif sg is not None:
            sats.Psat = float(griddata(self.sgCol, self.PCol, sg, method=method))
            P=sats.Psat
        elif sf is not None:
            sats.Psat = float(griddata(self.sfCol, self.PCol, sf, method=method))
            P=sats.Psat

        if P is not None:
            sats.Psat = P
            sats.Tsat = float(griddata((self.PCol), self.TCol, (P),method=method))
            sats.hf = float(griddata((self.PCol), self.hfCol, (P),method=method))
            sats.hg = float(griddata((self.PCol), self.hgCol, (P),method=method))
            sats.sf = float(griddata((self.PCol), self.sfCol, (P),method=method))
            sats.sg = float(griddata((self.PCol), self.sgCol, (P),method=method))
            sats.vf = float(griddata((self.PCol), self.vfCol, (P),method=method))
            sats.vg = float(griddata((self.PCol), self.vgCol, (P),method=method))
            sats.hgf = sats.hg - sats.hf
            sats.sgf = sats.sg - sats.sf
            sats.vgf = sats.vg - sats.vf
            sats.Psat *= 100.0  # convert to kPa
            return sats
        # given TSat, calc props
        if T is not None:
            sats.Tsat = T
            sats.Psat = float(griddata((self.TCol), self.PCol, (T),method=method))
            sats.hf = float(griddata((self.TCol), self.hfCol, (T),method=method))
            sats.hg = float(griddata((self.TCol), self.hgCol, (T),method=method))
            sats.sf = float(griddata((self.TCol), self.sfCol, (T),method=method))
            sats.sg = float(griddata((self.TCol), self.sgCol, (T),method=method))
            sats.vf = float(griddata((self.TCol), self.vfCol, (T),method=method))
            sats.vg = float(griddata((self.TCol), self.vgCol, (T),method=method))
            sats.hgf = sats.hg - sats.hf
            sats.sgf = sats.sg - sats.sf
            sats.vgf = sats.vg - sats.vf
            sats.Psat*=100.0  # convert to kPa
            return sats
    
    def getState(self, P_Bar=None, P_kPa=None, T=None, h=None, s=None, v=None, x=None):
        """
        general function to retrieve the thermodynamic state object (TSat, PSat, h, s, v)
        For this to work properly, need to give a (pressure or a temperature) + one other property
        :param P: pressure in bar (but output in kPa)
        :param T: temperature in C
        :param h: specific enthalpy in kJ/kg
        :param s: specific entropy in kJ/kg*K
        :param v: specific volume in m^3/kg
        :param x: quality
        :return: a thermodynamic state object
        """
        state = stateProps()
        state.region = 'saturated'
        sats=satProps()

        if P_Bar is not None:
            sats = self.getSatProps(P_Bar=P_Bar)
        if P_kPa is not None:
            sats=self.getSatProps(P_kPa=P_kPa)
        if T is not None:
            sats = self.getSatProps(T=T)
        state.P = sats.Psat
        state.T = sats.Tsat
        if h is not None:
            state.h = h
            state.x = (state.h - sats.hf) / sats.hgf
            state.s = sats.sf + state.x * sats.sgf
            state.v = sats.vf + state.x * sats.vgf
        elif s is not None:
            state.s = s
            state.x = (state.s - sats.sf) / sats.sgf
            state.h = sats.hf + state.x * sats.hgf
            state.v = sats.vf + state.x * sats.vgf
        elif v is not None:
            state.v = v
            state.x = (state.v - sats.vf) / sats.vgf
            state.h = sats.hf + state.x * sats.hgf
            state.s = sats.sf + state.x * sats.sgf
        elif x is not None:
            state.x = x
            state.h = sats.hf + state.x * sats.hgf
            state.s = sats.sf + state.x * sats.sgf
            state.v = sats.vf + state.x * sats.vgf
        return dc(state)

    def getCols(self):
        return self.TCol, self.PCol, self.hfCol, self.hgCol, self.sfCol, self.sgCol, self.vfCol, self.vgCol

class superheatedData():
    """
    I made this class for storing the superheated water table data for easy retrieval.  The
    superheated data table from superheated_water_table.txt is loaded once at the
    beginning.  A state object is returned from calling getState.
    """
    def __init__(self):
        self.readTable()

    def readTable(self):
        """
        Read and store the data from superheated_water_table.txt
        :return:
        """
        self.TCol, self.hCol, self.sCol, self.PCol = np.loadtxt('superheated_water_table.txt', skiprows=1, unpack=True)

        # no specific volume data in this table, so assume ideal gas behavior and calculate specific volume in SI
        self.vCol = np.ndarray(np.shape(self.TCol))
        MW = UC.MW_Water  # kg/kmol
        R = UC.R / MW  # kJ/kg*K -> kN*m/kg*K
        # P in kPa->kN/m^2
        # T in K
        # v =RT/p ->m^3/kg
        for i in range(len(self.TCol)):
            self.vCol[i] = R * (self.TCol[i] + 273.15) / self.PCol[i]

    def getState(self, P=None, T=None, h=None, s=None, v=None, method='cubic'):
        """
        general function to retrieve the thermodynamic state object (T, P, h, s, v).  Since there are
        5 variables, there are 5!/3!2! = 20/2= 10 combinations to worry about.
        :param P: pressure in kPa
        :param T: temperature in C
        :param h: specific enthalpy in kJ/kg
        :param s: specific entropy in kJ/kg*K
        :param v: specific volume in m^3/kg
        :param x: quality
        :return: a thermodynamic state object
        """
        state = stateProps()
        state.x=1.0 #if superheated, mass fraction vapor is 1.0
        state.region = 'superheated'
        #combo 1
        if P is not None and T is not None:
            state.T=T
            state.P=P
            state.h=float(griddata((self.TCol,self.PCol), self.hCol, (state.T, state.P), method=method))
            state.s=float(griddata((self.TCol,self.PCol), self.sCol, (state.T, state.P), method=method))
            state.v=float(griddata((self.TCol,self.PCol), self.vCol, (state.T, state.P), method=method))
        #combo 2
        elif P is not None and h is not None:
            state.h=h
            state.P=P
            state.T=float(griddata((self.hCol,self.PCol), self.TCol, (state.h, state.P), method=method))
            state.s=float(griddata((self.hCol,self.PCol), self.sCol, (state.h, state.P), method=method))
            state.v=float(griddata((self.hCol,self.PCol), self.vCol, (state.h, state.P), method=method))
        #combo 3
        elif P is not None and s is not None:
            state.s=s
            state.P=P
            state.T=float(griddata((self.sCol,self.PCol), self.TCol, (state.s, state.P), method=method))
            state.h=float(griddata((self.sCol,self.PCol), self.hCol, (state.s, state.P), method=method))
            state.v=float(griddata((self.sCol,self.PCol), self.vCol, (state.s, state.P), method=method))
        #combo 4
        elif P is not None and v is not None:
            state.v=v
            state.P=P
            state.T=float(griddata((self.vCol,self.PCol), self.TCol, (state.v, state.P), method=method))
            state.s=float(griddata((self.vCol,self.PCol), self.sCol, (state.v, state.P), method=method))
            state.h=float(griddata((self.vCol,self.PCol), self.hCol, (state.v, state.P), method=method))
        #combo 5
        elif T is not None and h is not None:
            state.h=h
            state.T=T
            state.p=float(griddata((self.hCol,self.TCol), self.PCol, (state.h, state.T), method=method))
            state.s=float(griddata((self.hCol,self.TCol), self.sCol, (state.h, state.T), method=method))
            state.v=float(griddata((self.hCol,self.TCol), self.vCol, (state.h, state.T), method=method))
        # combo 6
        elif T is not None and s is not None:
            state.s = s
            state.T = T
            state.P = float(griddata((self.sCol, self.TCol), self.PCol, (state.s, state.T), method=method))
            state.h = float(griddata((self.sCol, self.TCol), self.hCol, (state.s, state.T), method=method))
            state.v = float(griddata((self.sCol, self.TCol), self.vCol, (state.s, state.T), method=method))
        # combo 7
        elif T is not None and v is not None:
            state.v = v
            state.T = T
            state.P = float(griddata((self.vCol, self.TCol), self.PCol, (state.v, state.T), method=method))
            state.s = float(griddata((self.vCol, self.TCol), self.sCol, (state.v, state.T), method=method))
            state.h = float(griddata((self.vCol, self.TCol), self.hCol, (state.v, state.T), method=method))
        # combo 8
        elif h is not None and s is not None:
            state.s = s
            state.h = h
            state.P = float(griddata((self.sCol, self.hCol), self.PCol, (state.s, state.h), method=method))
            state.T = float(griddata((self.sCol, self.hCol), self.TCol, (state.s, state.h), method=method))
            state.v = float(griddata((self.sCol, self.hCol), self.vCol, (state.s, state.h), method=method))
        # combo 9
        elif h is not None and v is not None:
            state.v = v
            state.h = h
            state.P = float(griddata((self.vCol, self.hCol), self.PCol, (state.v, state.h), method=method))
            state.T = float(griddata((self.vCol, self.hCol), self.TCol, (state.v, state.h), method=method))
            state.s = float(griddata((self.vCol, self.hCol), self.sCol, (state.v, state.h), method=method))
        # combo 10
        elif s is not None and v is not None:
            state.v = v
            state.s = s
            state.P = float(griddata((self.vCol, self.sCol), self.PCol, (state.v, state.s), method=method))
            state.T = float(griddata((self.vCol, self.sCol), self.TCol, (state.v, state.s), method=method))
            state.h = float(griddata((self.vCol, self.sCol), self.sCol, (state.v, state.s), method=method))
        return dc(state)

    def getCols(self):
        self.TCol, self.hCol, self.sCol, self.PCol

class subcooled():
    """
    A subcooled liquid can be modeled as an incompressible substance:
    u=uf(T), v=vf(T), s=sf(T), h=hf(T)+(P-PSat(T))*vf(T)
    """
    def __init__(self):
        self.satData=None

    def getState(self, PLowSat, PHighSat, P, T):
        """
        For Rankine, we exit pump at P=PHigh, T4=T3=Tsat for PLow.
        #case 1: PLow<=P<=PHigh
        For the states between P=PLow to PHigh, T=T3: h=h3+(P-PLow)*v3, v=v3 (i.e., incompressible), s=s3 (i.e., isentropic efficiency of pump = 1.0)
        #case 2: P=PHigh, T3<=T<=T5
        Between states (P=PHigh, T4) to (P=PHigh, T5=Tsat,PHigh), I will assume P and T vary linearly, so:
            z=(T-T4)/(T5-T4)
            h4=h3+(PHigh-PLow)v3
            h5, s5, v5=hf, sf, vf for PHigh
            h=h4 +(h5-h4)*z
            s=s3+(s5-s3)*z
            v=v3+(v5-v3)*z
        general function to retrieve the thermodynamic state object (TSat, P, h, s, v)
        :param PLow: in kPa
        :param PHigh: in kPa
        :param P: in kPa
        :param T: in C
        :return: a deep copy of a thermodynamic state object
        """
        state=stateProps()
        case = 1 if (T<=PLowSat.Tsat) else 2
        if case ==1:
            z=(P-PLowSat.Psat)/(PHighSat.Psat-PLowSat.Psat)
            state.h=PLowSat.hf+z*(PHighSat.Psat-PLowSat.Psat)*PLowSat.vf
            state.s=PLowSat.sf
            state.v=PLowSat.vf
            state.T=PLowSat.Tsat
            state.P=P
        else:
            z=(T-PLowSat.Tsat)/(PHighSat.Tsat-PLowSat.Tsat)
            h4=PLowSat.hf+(PHighSat.Psat-PLowSat.Psat)*PLowSat.vf
            s4=PLowSat.sf
            v4=PLowSat.vf
            h5=PHighSat.hf
            s5=PHighSat.sf
            v5=PHighSat.vf
            state.h=h4+z*(h5-h4)
            state.s=s4+z*(s5-s4)
            state.v=v4+z*(v5-v4)
            state.P=P
            state.T=T
        name='subcooled'
        state.x=-0.1
        return dc(state)

class satProps():
    """
    For storage and retrieval of saturated properties at a given isobar or isotherm
    """
    def __init__(self):
        self.Tsat = None
        self.Psat = None
        self.hf = None
        self.hg = None
        self.hgf = None
        self.sf = None
        self.sg = None
        self.sgf = None
        self.vf = None
        self.vg = None
        self.vgf = None

    def set(self, vals):
        self.Tsat, self.Psat, self.hf, self.hg, self.sf, self.sg, self.vf, self.vg = vals
        self.hgf = self.hg - self.hf
        self.sgf = self.sg - self.sf
        self.vgf = self.vg - self.vf

    def get(self):
        return (self.Tsat, self.Psat, self.hf, self.hg, self.hgf, self.sf, self.sg, self.sgf, self.vf, self.vg, self.vgf)

    def getTextOutput(self, SI=True):
        """
        Sets the self.txtOut string for display.
        :param SI:
        :return: self.txtOut
        """
        if SI is False:
            P = self.PSat*self.UC.bar_to_psi
            PUnits="psi"
            T = self.UC.C_to_F(self.TSat)
            TUnits="F"
            hf=self.hf*self.UC.kJperkg_to_BTUperlb
            hg=self.hg*self.UC.kJperkg_to_BTUperlb
            HUnits="BTU/lb"
            sf=self.sf*self.UC.kJperkgK_to_BTUperlbR
            sg=self.sg*self.UC.kJperkgK_to_BTUperlbR
            SUnits="BTU/lb*R"
            vf=self.vf*self.UC.m3perkg_to_ft3perlb
            vg=self.vg*self.UC.m3perkg_to_ft3perlb
            VUnits="ft^3/lb"
        else:
            P = self.PSat
            PUnits = "bar"
            T = self.TSat
            TUnits = "C"
            hf = self.hf
            hg = self.hg
            HUnits = "kJ/kg"
            sf = self.sf
            sg = self.sg
            SUnits = "kJ/kg*K"
            vf = self.vf
            vg = self.vg
            VUnits = "m^3/kg"

        self.txtOut = "PSat = {:0.2f} {}, TSat = {:0.2f} {}".format(P, PUnits, T, TUnits)
        self.txtOut += "\nhf = {:0.2f} {}, hg = {:0.2f} {}".format(hf, HUnits,hg,HUnits)
        self.txtOut += "\nsf = {:0.2f} {}, sg = {:0.2f} {}".format(sf,SUnits,sg,SUnits)
        self.txtOut += "\nvf = {:0.4f} {}, vg = {:0.2f} {}".format(vf,VUnits,vg,VUnits)
        return self.txtOut
    
class stateProps():
    """
    for storage and retrieval of a thermodynamic state
    T (C), P (kPa), h (kJ/kg), s (kJ/kg*K), v (m^3/kg), x (dimensionless)
    """
    def __init__(self):
        self.name = None
        self.T = None
        self.P = None
        self.h = None
        self.s = None
        self.v = None
        self.x = None
        self.region = None

    def getVal(self, name='T', SI=True ):
        if SI:
            hCF=1
            sCF=1
            vCF=1
            pCF=1
        else:
            hCF=UC.kJperkg_to_BTUperlb
            sCF=UC.kJperkgK_to_BTUperlbR
            vCF=UC.m3perkg_to_ft3perlb
            pCF=UC.kpa_to_psi

        n=name.lower()
        if n == 't':
            return self.T if SI else UC.C_to_F(self.T)
        if n == 'h':
            return self.h*hCF
        if n == 's':
            return self.s*sCF
        if n == 'v':
            return self.v*vCF
        if n == 'p':
            return self.P*pCF

    def print(self):
        if self.name is not None:
            print(self.name)
        if self.x is None or self.x < 0.0:
            print('Region: compressed liquid')
            print('p = {:0.2f} kPa'.format(self.P))
            print('h = {:0.2f} kJ/kg'.format(self.h))
        else:
            print('Region: ', self.region)
            print('p = {:0.2f} kPa'.format(self.P))
            print('T = {:0.1f} degrees C'.format(self.T))
            print('h = {:0.2f} kJ/kg'.format(self.h))
            print('s = {:0.4f} kJ/(kg K)'.format(self.s))
            print('v = {:0.6f} m^3/kg'.format(self.v))
            print('x = {:0.4f}'.format(self.x))
        print()

class SatPropsIsobar():
    def __init__(self, P, SI=True, interpMethod='cubic'):
        """
        Sets saturation properties for a given isobar
        :param P:  in kPa
        """
        tscol, pscol, hfcol, hgcol, sfcol, sgcol, vfcol, vgcol = np.loadtxt('sat_water_table.txt', skiprows=1, unpack=True)
        self.PSat=P/100
        self.TSat=float(griddata(pscol, tscol, self.PSat, method=interpMethod))
        self.hf=float(griddata(pscol, hfcol, self.PSat, method=interpMethod))
        self.hg=float(griddata(pscol, hgcol, self.PSat, method=interpMethod))
        self.sf=float(griddata(pscol, sfcol, self.PSat, method=interpMethod))
        self.sg=float(griddata(pscol, sgcol, self.PSat, method=interpMethod))
        self.vf=float(griddata(pscol, vfcol, self.PSat, method=interpMethod))
        self.vg=float(griddata(pscol, vgcol, self.PSat, method=interpMethod))
        # for doing unit conversions
        self.getTextOutput(SI=SI)

    def getTextOutput(self, SI=True):
        """
        Sets the self.txtOut string for display.
        :param SI:
        :return: self.txtOut
        """
        if SI is False:
            P = self.PSat*UC.bar_to_psi
            PUnits="psi"
            T = UC.C_to_F(self.TSat)
            TUnits="F"
            hf=self.hf*UC.kJperkg_to_BTUperlb
            hg=self.hg*UC.kJperkg_to_BTUperlb
            HUnits="BTU/lb"
            sf=self.sf*UC.kJperkgK_to_BTUperlbR
            sg=self.sg*UC.kJperkgK_to_BTUperlbR
            SUnits="BTU/lb*R"
            vf=self.vf*UC.m3perkg_to_ft3perlb
            vg=self.vg*UC.m3perkg_to_ft3perlb
            VUnits="ft^3/lb"
        else:
            P = self.PSat
            PUnits = "bar"
            T = self.TSat
            TUnits = "C"
            hf = self.hf
            hg = self.hg
            HUnits = "kJ/kg"
            sf = self.sf
            sg = self.sg
            SUnits = "kJ/kg*K"
            vf = self.vf
            vg = self.vg
            VUnits = "m^3/kg"

        self.txtOut = "PSat = {:0.2f} {}, TSat = {:0.2f} {}".format(P, PUnits, T, TUnits)
        self.txtOut += "\nhf = {:0.2f} {}, hg = {:0.2f} {}".format(hf, HUnits,hg,HUnits)
        self.txtOut += "\nsf = {:0.2f} {}, sg = {:0.2f} {}".format(sf,SUnits,sg,SUnits)
        self.txtOut += "\nvf = {:0.4f} {}, vg = {:0.2f} {}".format(vf,VUnits,vg,VUnits)
        return self.txtOut

class StateDataForPlotting:
    """
    I'm making this class for easy storage of data for plotting.
    """
    def __init__(self):
        self.T = []
        self.P = []
        self.h = []
        self.s = []
        self.v = []

    def clear(self):
        self.T.clear()
        self.P.clear()
        self.h.clear()
        self.s.clear()
        self.v.clear()

    def addPt(self, vals):
        """
        adds a thermodynamic state point to the list
        :param vals: a list or tuple with T, P, h, s, v in that order
        :return:
        """
        T, P, h, s, v = vals
        self.T.append(T)
        self.P.append(P)
        self.h.append(h)
        self.s.append(s)
        self.v.append(v)

    def getAxisLabel(self, W='T', SI=True):
        w=W.lower()
        if w == 't':
            return r'T $\left(^oC\right)$' if SI else r'T $\left(^oF\right)$'
        if w == 'h':
            return r'h $\left(\frac{kJ}{kg}\right)$' if SI else r'h $\left(\frac{BTU}{lb}\right)$'
        if w == 's':
            return r's $\left(\frac{kJ}{kg\cdot K}\right)$' if SI else r's $\left(\frac{BTU}{lb\cdot ^oR}\right)$'
        if w == 'v':
            return r'v $\left(\frac{m^3}{kg}\right)$' if SI else r'v $\left(\frac{ft^3}{lb}\right)$'
        if w == 'p':
            return r'P $\left(kPa\right)$' if SI else r'P (psi)'

    def getDataCol(self, W='T', SI=True):

        if SI:
            hCF=1
            sCF=1
            vCF=1
            pCF=1
        else:
            hCF=UC.kJperkg_to_BTUperlb
            sCF=UC.kJperkgK_to_BTUperlbR
            vCF=UC.m3perkg_to_ft3perlb
            pCF=UC.kpa_to_psi

        w=W.lower()
        if w=='t':
            return self.T if SI else [UC.C_to_F(t) for t in self.T]
        if w=='h':
            return np.array(self.h)*hCF
        if w=='s':
            return np.array(self.s)*sCF
        if w=='v':
            return np.array(self.v)*vCF
        if w=='p':
            return np.array(self.P)*pCF

class Steam_SI:
    def __init__(self, P=None, T=None, x=None, v=None, h=None, s=None, name=None):
        """
        This is a general steam class for sub-critical (i.e., superheated and saturated) properties of steam.
        The user may specify any two properties to calculate all other properties of the steam.
        Note: we have 6 properties, but only can specify two of them.  Combinations=6!/(2!4!)=15

        I handle all cases in self.calc

        :param P: Pressure (kPa)
        :param T: Temperature (C)
        :param x: Quality
        :param v: Specific Volume (kg/m^3)
        :param h: Enthalpy (kJ/kg)
        :param s: Entropy (kJ/(kg*K))
        :param name:
        """
        self.state=stateProps()
        self.state.P = P  # pressure - kPa
        self.state.T = T  # Temperature - degrees C
        self.state.x = x  # quality (a value between 0 and 1)
        self.state.v = v  # specific volume - m^3/kg
        self.state.h = h  # enthalpy - kJ/kg
        self.state.s = s  # entropy - kJ/(kg K)
        self.state.name = name # a useful identifier
        self.state.region = None # 'superheated' or 'saturated'
        self.RW=UC.R/UC.MW_Water #water gas constant kJ/(kg*K)
        self.SH_Data = superheatedData()
        self.Sat_Data = saturatedData()
        self.SC_Data = subcooled()
        self.getState(P=P,T=T,x=x,v=v,h=h,s=s,name=name)
        pass

    def getState(self, P=None, T=None, x=None, v=None, h=None, s=None, name=None):
        """
        In principle, there are 15 cases to handle any pair of variables.
        I depend on user to specify proper set of variables
        :param P: Pressure (kPa)
        :param T: Temperature (C)
        :param x: Quality
        :param v: Specific Volume (kg/m^3)
        :param h: Enthalpy (kJ/kg)
        :param s: Entropy (kJ/(kg*K))
        :param name:
        :return: True if properties calculated
        """
        #region select case
        case=None
        if P is not None:  # pressure is specified
            if T is not None:
                case="PT"
            elif x is not None:
                case="Px"
            elif v is not None:
                case="Pv"
            elif h is not None:
                case="Ph"
            elif s is not None:
                case="Ps"
        elif case is None and T is not None:   #temperature is specified
            if x is not None:
                case="Tx"
            elif v is not None:
                case="Tv"
            elif h is not None:
                case="Th"
            elif s is not None:
                case="Ts"
        elif case is None and x is not None:  # quality is specified
            if v is not None:
                case ="xv"
            elif h is not None:
                case="xh"
            elif s is not None:
                case="xs"
        elif case is None and v is not None:  # quality is specified
            if h is not None:
                case="vh"
            elif s is not None:
                case="vs"
        elif case is None and h is not None:  # enthalpy is specified
            if s is not None:
                case="hs"
        if case is None:
            return False
        # my 15 cases are:  PT, Px, Pv, Ph, Ps, Tx, Tv, Th, Ts, xv, xh, xs, vh, vs, hs
        #endregion

        #region calculate properties based on case
        # if P given (easy to check saturated vs. superheated)
        if case.__contains__("P"):
            satProps=self.Sat_Data.getSatProps(P_kPa=P)
            if case.__contains__("x"):  #quality given
                self.state=self.Sat_Data.getState(P_kPa=P,x=x)
                return dc(self.state)
            if case.__contains__("v"):  #find if saturated or superheated
                xval=(v-satProps.vf)/(satProps.vgf)
                if xval<=1:  #saturated
                    self.state= self.Sat_Data.getState(P_kPa=P,x=xval)
                    return dc(self.state)
                else:  #superheated
                    self.region = "superheated"
                    self.x=1.0
                    self.state=self.SH_Data.getState(P=P, v=v)
                    return dc(self.state)
            if case.__contains__("h"):  #find if saturated or superheated
                xval=(h-satProps.hf)/(satProps.hgf)
                if xval<=1:  #saturated
                    self.state=self.Sat_Data.getState(P_kPa=P, h=h)
                    return dc(self.state)
                else:  #superheated
                    self.state=self.SH_Data.getState(P=P, h=h)
                    return dc(self.state)
            if case.__contains__("s"):  # find if saturated or superheated
                xval = (s - satProps.sf) / (satProps.sgf)
                if xval <= 1 and xval >=0:  # saturated
                    self.state= self.Sat_Data.getState(P_kPa=P, s=s)
                    return dc(self.state)
                elif xval<0:  # subcooled
                    PLowSat=self.Sat_Data.getSatProps(sf=s)
                    PHighSat=self.Sat_Data.getSatProps(P_kPa=P)
                    self.state = self.SC_Data.getState(PLowSat=PLowSat,PHighSat=PHighSat,P=P,T=PLowSat.Tsat)
                    return dc(self.state)
                else:  # superheated
                    self.state= self.SH_Data.getState(P=P, s=s)
                    return dc(self.state)
            if case.__contains__("T"):  # find if satruated or superheated steam
                if T==satProps.Tsat:  # generally, this is an indeterminate case, but assume saturated vapor
                    self.state= self.Sat_Data.getState(P_kPa=P, x=1.0)
                    return dc(self.state)
                elif T>satProps.Tsat:  #superheated
                    self.state= self.SH_Data.getState(P=P, T=T)
                    return dc(self.state)
                else: # sub-cooled, so estimate properties
                    satProps= self.Sat_Data.getSatProps(T=T)
                    psat = satProps.Psat
                    self.state.region = "subcooled"
                    self.state.x=0
                    self.state.h=satProps.hf+(self.P-psat)*satProps.vf
                    self.s=satProps.sf
                    self.v=satProps.vf
                    return dc(self.state)
        # if T given (easy to check saturated vs. superheated)
        if case.__contains__("T"):
            # Using the known Temperature, interpolate on the saturation tables columns
            # at the known pressure
            satProps = self.Sat_Data.getSatProps(T=T)
            if case.__contains__("x"):  # quality given
                self.state= self.Sat_Data.getState(T=T,x=x)
                return dc(self.state)
            if case.__contains__("v"):  # find if saturated or superheated
                xval = (self.v - satProps.vf) / (satProps.vgf)
                if xval <= 1:  # saturated
                    self.state =  self.Sat_Data.getState(T=T, x=xval)
                    return dc(self.state)
                else:  # superheated
                    self.state =  self.SH_Data.getState(T=T, v=v)
                    return dc(self.state)
            if case.__contains__("h"):  # find if saturated or superheated
                xval = (self.h - satProps.hf) / (satProps.hgf)
                if xval <= 1:  # saturated
                    self.state =  self.Sat_Data.getState(T=T, h=h)
                    return dc(self.state)
                else:  # superheated
                    self.state =  self.SH_Data.getState(T=T, h=h)
                    return dc(self.state)
            if case.__contains__("s"):  # find if saturated or superheated
                xval = (self.s - satProps.sf) / (satProps.sgf)
                if xval <= 1:  # saturated
                    self.state =  self.Sat_Data.getState(T=T, s=s)
                    return dc(self.state)
                else:  # superheated
                    self.state =  self.SH_Data.getState(T=T, s=s)
                    return dc(self.state)
            return
        # if quality given (easy case) for saturated
        if case.__contains__("x"):
            if case.__contains__("v"):
                self.state =  self.Sat_Data.getState(x=x, v=v)
                return dc(self.state)
            if case.__contains__("h"):  #find if saturated or superheated
                self.state =  self.Sat_Data.getState(x=x, h=h)
                return dc(self.state)
            if case.__contains__("s"):  # find if saturated or superheated
                self.state =  self.Sat_Data.getState(x=x, s=s)
                return dc(self.state)
            return
        # if vh, vs, or hs (searching required to determine if saturated)
        if case.__contains__("v"):
            twophase=False
            if v < self.Sat_Data.vc:
                twophase=True
            else:  #might be two-phase
                SatProps=self.Sat_Data.getSatProps(vg=v)
                if case.__contains__("h"):
                    if self.h<=SatProps.hg:  # means two phase
                        twophase=True
                else:
                    if self.s <= SatProps.sg:  # means two phase
                        twophase=True
            if twophase:
                if case.__contains__("h"):
                    self.state = self.Sat_Data.getState(v=v, h=h)
                    return dc(self.state)
                if case.__contains__("s"):
                    self.state = self.Sat_Data.getState(v=v, s=s)
                    return dc(self.state)
            else:  # means superheated
                self.state = self.SH_Data.getState(v=v, h=h)
                return dc(self.state)
        if case.__contains__("h"):  # this case has h & s given
            #see Mollier diagram.
            SatProps=self.Sat_Data.getSatProps(sg=s)
            if self.h<=SatProps.hg:
                self.state = self.Sat_Data.getState(s=s, h=h)
                return dc(self.state)
            else:  # means superheated
                self.state = self.SH_Data.getState(s=s, h=h)
                return dc(self.state)
        #endregion

    def igl_v(self):
        # ideal gas law V=RT/P-> v=(R/MW)*(T+273)/(P)
        # calculate a column of specific volume for the superheated water table using ideal gas
        return self.RW * (self.T + 273) / self.P

    def print(self):
        if self.name is not None:
            print('Name: {}'.format(self.name))
        if self.region is not None:
            print('Region: {}'.format(self.region))   
        if self.P is not None:
            print('p = {:.2f} kPa'.format(self.P))
        if self.T is not None:
            print('T = {:.1f} degrees C'.format(self.T))
        if self.h is not None:
            print('h = {:.2f} kJ/kg'.format(self.h))
        if self.s is not None:
            print('s = {:.4f} kJ/(kg K)'.format(self.s))
        if self.v is not None:
            print('v = {:.6f} m^3/kg'.format(self.v))
        if self.x is not None:
            print('x = {:.4f}'.format(self.x))
        print('')

def main():
    inlet=Steam_SI(P=7350, name='Turbine Inlet')
    inlet.x =0.9  # 90 percent quality  - accessed DIRECTLY  - not through a method
    inlet.calc()
    inlet.print()
    h1 = inlet.h; s1 = inlet.s
    print(h1,s1,'\n')

    outlet=Steam_SI(P=100, s=inlet.s, name='Turbine Exit')
           #notice -  s=inlet.s
    outlet.print()


    another=Steam_SI(P=8575, h=2050, name='State 3')
    another.print()

    yetanother=Steam_SI(P=8900, h=41250, name='State 4')
    yetanother.print()

    g1=Steam_SI(P=800, x=1.0, name='Gap1')
    g1.print()
    g2=Steam_SI(P=g1.P, s=g1.s * 1.0012, name='Gap2')
    g2.print()
    g2=Steam_SI(P=g1.P, s=6.6699, name='Gap3')
    g2.print()




    #uncommenting the next two lines will cause a ValueError error
    #final1 = State(7000, T=250, name='State 5')
    #final1.print()

if __name__ == "__main__":
   main()

