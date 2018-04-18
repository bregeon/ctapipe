#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple Lidar data file reader (HESS txt file format)
       Basic Klett reconstruction

@author: Johan Bregeon
"""
import os
import sys
import numpy
import math
from datetime import datetime

import logging
logger = logging.getLogger(__name__)


class pLidarRun(object):
    """ Class to manage Lidar run

    Parameters
    ----------
    None
    """
    def __init__(self):
        self.FileName = None
        self.NBins = None
        self.RunNumber = None
        self.DateTime = None

    def readFile(self, file_path):
        """
        Wrapper to choose which method to use to read the input file
        from either an ascii or a FITS file, based on the file name extension

        fits reader not implemented yet as format is not defined yet

        Parameters
        ----------
        file_path : string
            a file path
        """
        self.FileName=file_path
        # choose which method to call according to file name extension
        ext=os.path.basename(self.FileName).split('.')[-1]
        if ext == 'fits':
            self.readFITSFile()
        elif ext in ['txt', 'text']:
            self.readTxtFile()
        else:
            logger.error('Unknown file extension: %s'%ext)
            logger.error('Can read only .txt or .fits files. Aborting')
            sys.exit(2)

    def readTxtFile(self):
        """
        Read Lidar data from an ascii file

        Parameters
        ----------
        None
        """
        cont=open(self.FileName,'r').readlines()
        self.RunNumber=100 #int(self.FileName.split('_')[1])
        self.DateTime=datetime.strptime(cont[0].strip(),'%a %b %d %H:%M:%S %Y')
        logger.info('-'*80)
        logger.info('Reading file %s, Date %s, Run %s'%\
              (self.FileName,self.DateTime,self.RunNumber))
        self.NPoints=len(cont[1:])
        self.RawData=numpy.ndarray((self.NPoints,3),'d')
        self.RawAltitude=numpy.ndarray((self.NPoints),'d')
        self.RawWL1=numpy.ndarray((self.NPoints),'d')
        self.RawWL2=numpy.ndarray((self.NPoints),'d')
        i=0
        for l in cont[1:]:
            data=l.split()
            alt=float(data[0])
            wl1=float(data[1])
            wl2=float(data[2])
            self.RawData[i]=(alt,wl1,wl2)
            self.RawAltitude[i]=alt
            self.RawWL1[i]=wl1
            self.RawWL2[i]=wl2
            i+=1
        #Get pointers to slices for convenience - but TGraph unhappy afterward
        #self.RawAltitude=self.RawData[:,0]
        #self.RawWL1=self.RawData[:,1]
        #self.RawWL2=self.RawData[:,2]

    def readFITSFile(self):
        """
        Read Lidar data from a FITS file,
        not implemented yet

        Parameters
        ----------
        None
        """
        logger.info('Reading FITS file')
        logger.error('Not implemented')
        sys.exit(2)

    def process(self, nBins=100):
        """
        Full processing including Klett inversion

        Parameters
        ----------
        nBins : int
            number of bins to use for calculations
        """
        self.NBins=nBins
        self.reduce()
        self.lnPower()
        self.binData(nBins)
        alpha_wl1, alpha_wl2=self.run_Klett()
        return (alpha_wl1, alpha_wl2)

    def reduce(self,altmin=0.800, altmax=10., bkgmin=20., bkgmax=25.):
        """
        Subtract background and calculate power

        Parameters
        ----------
        altmin : float
            minimum altitude to consider for signal
        altmax : float
            maximum altitute to consider for signal
        bkgmin : float
            minimum altitude to consider for background estimation
        bkgmax : float
            maximum altitute to consider for background estimation
        """
        logger.info('Subtract background and produce reduced data and power arrays')
        self.AltMin=altmin
        self.AltMax=altmax
        self.BkgMin=bkgmin
        self.BkgMax=bkgmax
        self.BkgMinIndex=self.NPoints-sum(self.RawAltitude>self.BkgMin)
        self.BkgMaxIndex=self.NPoints-sum(self.RawAltitude>self.BkgMax)
        self.AltMinIndex=self.NPoints-sum(self.RawAltitude>self.AltMin)
        self.AltMaxIndex=self.NPoints-sum(self.RawAltitude>self.AltMax)
        self.BkgWL1=self.RawWL1[self.BkgMinIndex:self.BkgMaxIndex].mean()
        self.BkgWL2=self.RawWL2[self.BkgMinIndex:self.BkgMaxIndex].mean()
        self.WL1=abs(self.RawWL1[self.AltMinIndex:self.AltMaxIndex]-self.BkgWL1)
        self.WL2=abs(self.RawWL2[self.AltMinIndex:self.AltMaxIndex]-self.BkgWL2)
        self.Alt=self.RawAltitude[self.AltMinIndex:self.AltMaxIndex]
        self.PW1=self.WL1*self.Alt**2
        self.PW2=self.WL2*self.Alt**2
        self.NData=self.Alt.size

    def lnPower(self):
        """
        Compute the natural logarithm of the signal power

        Parameters
        ----------
        None
        """
        logger.info('Calculate Ln of Power')
        self.LnPW1=numpy.ndarray((self.NData),'d')
        self.LnPW2=numpy.ndarray((self.NData),'d')
        for i in range(self.NData):
            self.LnPW1[i]=math.log(self.PW1[i])
            self.LnPW2[i]=math.log(self.PW2[i])

    def binData(self, nBins=50):
        """
        Bin data in interval of altitudes

        Parameters
        ----------
        nBins : int
            number of bins
        """
        logger.info('Bin data')
        self.NBins=nBins
        # get bin width in Log scale
        self.BinLnAltWidth=(math.log(self.AltMax)-math.log(self.AltMin))/self.NBins
        # bins edges array
        self.BinsAlt=numpy.ndarray((self.NBins+1),'d')
        for i in range(self.NBins+1):
            self.BinsAlt[i]=math.exp(math.log(self.AltMin)+i*self.BinLnAltWidth)
        # bins min, max, center
        self.BinsAltMin=self.BinsAlt[:-1]
        self.BinsAltMax=self.BinsAlt[1:]
        self.BinsAltCenter=numpy.ndarray((self.NBins),'d')
        for i in range(self.NBins):
            self.BinsAltCenter[i]=(self.BinsAlt[i]+self.BinsAlt[i+1])/2.

        # array to put average PW and LnPw for each bin in altitude
        self.BinnedPW1=numpy.ndarray((self.NBins),'d')
        self.BinnedPW2=numpy.ndarray((self.NBins),'d')
        self.BinnedLnPW1=numpy.ndarray((self.NBins),'d')
        self.BinnedLnPW2=numpy.ndarray((self.NBins),'d')
        self.BinnedNpoints=numpy.ndarray((self.NBins),'d')
        sumpw1=0
        sumpw2=0
        sumlnpw1=0
        sumlnpw2=0
        nPoints=0
        kAlt=0
        #print 'Bin \t altMin  altMean altMax  SumPW \t\t NPoints'
        for alt,pw1,pw2,lnpw1,lnpw2 in zip(self.Alt,self.PW1, self.PW2, self.LnPW1,self.LnPW2):
            if alt>self.BinsAlt[kAlt]:
                sumpw1+=pw1
                sumpw2+=pw2
                sumlnpw1+=lnpw1
                sumlnpw2+=lnpw2
                nPoints+=1
            if alt>self.BinsAlt[kAlt+1] or alt==self.Alt[-1]:
                #print '%d \t %.03f \t %.03f \t %.03f \t %.05f \t %d'%\
                #(kAlt,self.BinsAlt[kAlt],self.BinsAltCenter[kAlt],self.BinsAlt[kAlt+1], sumpw1, nPoints)
                self.BinnedNpoints[kAlt]=nPoints
                self.BinnedPW1[kAlt]=sumpw1/nPoints
                self.BinnedPW2[kAlt]=sumpw2/nPoints
                self.BinnedLnPW1[kAlt]=sumlnpw1/nPoints
                self.BinnedLnPW2[kAlt]=sumlnpw2/nPoints
                sumpw1=pw1
                sumpw2=pw2
                sumlnpw1=lnpw1
                sumlnpw2=lnpw2
                nPoints=1
                kAlt+=1
        logger.info('Binned data ready.')

    def run_Klett(self, r0=10, alpha0wl1=0.0038, alpha0wl2=0.018, k=1, l=1):
        """
        Klett implementation:
        P(r)=V(r)*r2
        alpha(r)=P(r)/( P(r0)/alpha0 - 2*sum{r0,r}{P(h).dh} )

        Parameters
        ----------
        r0 : float
            reference altitude
        alpha0wl1 : float
            alpha0 for wavelength 1
        alpha0wl2 : float
            alpha0 for wavelength 2
        k : float
            as in beta=l*alpha^k
        l : float
            as in beta=l*alpha^k
        """
        logger.info('Klett: starting inversion')
        # save parameters
        self.Klett_r0=r0
        self.Klett_alpha0wl1=alpha0wl1
        self.Klett_alpha0wl2=alpha0wl2
        self.Klett_k=k
        self.Klett_l=l
        # get altitude bin closest to r0, convert to int for TGraph
        self.AlphaNBins=int(self.NBins-sum(self.BinsAltCenter>r0)-1)
        r0=self.BinsAltCenter[self.AlphaNBins]
        # array for binned alpha(r)
        self.BinnedAlphaPW1=numpy.ndarray((self.AlphaNBins),'d')
        self.BinnedAlphaPW2=numpy.ndarray((self.AlphaNBins),'d')
        # first bin
        self.BinnedAlphaPW1[self.AlphaNBins-1]=alpha0wl1
        self.BinnedAlphaPW2[self.AlphaNBins-1]=alpha0wl2
        # iterate
        #print 'n   pw \t\t alpha_n \t int_pw \t int_alpha \t binAlt'
        for i in range(self.AlphaNBins-2, -1, -1):
            # WL1
            pw_m1=self.BinnedPW1[i+1]
            alpha_m1=self.BinnedAlphaPW1[i+1]
            int_pw_m1=(self.BinnedPW1[i+1]+self.BinnedPW1[i])/2.*(self.BinsAltCenter[i+1]-self.BinsAltCenter[i])
            self.BinnedAlphaPW1[i]=self.BinnedPW1[i]/(pw_m1/alpha_m1-2*int_pw_m1)
            #print '%d %0.5f \t %.05f \t %.05f \t %.05f \t %.03f'%\
            #(i, pw_m1, alpha_m1, int_pw_m1, self.BinnedAlphaPW1[i], self.BinsAltCenter[i])
            # WL2
            pw_m2=self.BinnedPW2[i+1]
            alpha_m2=self.BinnedAlphaPW2[i+1]
            int_pw_m2=(self.BinnedPW2[i+1]+self.BinnedPW2[i])/2.*(self.BinsAltCenter[i+1]-self.BinsAltCenter[i])
            self.BinnedAlphaPW2[i]=self.BinnedPW2[i]/(pw_m2/alpha_m2-2*int_pw_m2)

        # Truncate BinsAltCenter
        #self.BinsAltCenter= self.BinsAltCenter[:-1]
        logger.info('Klett: inversion done.')
        return (self.BinnedAlphaPW1, self.BinnedAlphaPW2)

    def TauAndQuality(self, h=4):
        """
        Integrate Alpha to get absorption for the first h km

        Parameters
        ----------
        h : float
            maximum altitude for integration in km
        """
        nBinToSum=int(self.AlphaNBins-sum(self.BinsAltCenter>h))
        self.Tau4WL1=sum(self.BinnedAlphaPW1[:nBinToSum]*(self.BinsAltMax[:nBinToSum]-self.BinsAltMin[:nBinToSum]))
        self.Tau4WL2=sum(self.BinnedAlphaPW2[:nBinToSum]*(self.BinsAltMax[:nBinToSum]-self.BinsAltMin[:nBinToSum]))
        self.IsGood=True
        if self.Tau4WL1<0.002 or self.Tau4WL2<0.01:
            self.IsGood=False
        return (self.Tau4WL1,self.Tau4WL2)

    def calcTau4(self):
        """
        Wrapper to Integrate Alpha to get absorption for the first 4 km

        Parameters
        ----------
        None
        """
        return self.TauAndQuality(h=4)

    def calcTransmission(self, h=4):
        """
        Calculate transmission from optical depth

        Parameters
        ----------
        h :
            maximum altitude to integrate to in km
        """
        self.TauAndQuality(h)
        return (math.exp(-2*self.Tau4WL1),math.exp(-2*self.Tau4WL2))

    def dumpToDict(self):
        """
        Writes summary results to a string dictionary

        Parameters
        ----------
        None
        """
        stringDict='%d'%self.RunNumber
        stringDict+=':{"FileName":"%s","DateTime":"%s", "NBins":%d,'%\
                     (os.path.basename(self.FileName), self.DateTime, self.NBins)
        stringDict+='"Bkg":(%.5f,%.5f), "Tau4":(%.5f,%.5f), "IsGood":%s},\n'%\
                     (self.BkgWL1,self.BkgWL2, self.Tau4WL1, self.Tau4WL2, self.IsGood)

        return stringDict

if __name__ == '__main__':

    # set CTAPIPE_SVC_PATH to a path containing the file below
    from traitlets import Unicode
    from ctapipe.utils import get_dataset
    lidar_file_path = get_dataset('hess_elastic_lidar_data.txt')

    # JB - can't make this work...
    #input_path = Unicode(get_dataset('atmprof26.dat'), allow_none=True,
    #                     help='Path to the atmospheric profile file, e.g. '
    #                          'atmprof26.dat').tag(config=True)

    r = pLidarRun()
    r.readFile(lidar_file_path)
    alpha_wl1, alpha_wl2 = r.process()
    print('Date Time: %s'%r.DateTime)
    print('Background: %.5f %.5f'%(r.BkgWL1, r.BkgWL2))
    print('Tau4: %.2f %.2f'%r.calcTau4())
    print('Prob: %.2f %.2f'%r.calcTransmission(h=4))
