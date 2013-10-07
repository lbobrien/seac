       program seac

c***********************************************************************
c                              SEAC
c         Simulation of Environmentally-Assisted Cracking
c***********************************************************************

c  Version : SEAC 1.1           950202 

c  Author  : John H. Chun

c  seac calculates crevice chemstry using:

c     dC/dt = -   Divergence
c             +/- Generation & Consumption by Chemcial Rx
c             +/- Generation & Consumption by Wall Electrochemcial Rx 
c             +/- Generation & Consumption by Wall Chemical Dissolution
c             +   Generation by Radiolysis

c     with Initial Conditions at t = TInit
c       Concentrations at all mesh points are set to Bulk Concentrations

c     and Boundaray Conditions at x = 0.0 (crevice mouth)
c       Fixed Concentration at Bulk Concentrations for all species.

c     and Boundary Conditions at x = L (crevice tip)
c       Zero flux

c     but these appear through the source term
c       Flux due to Electrochemical Rx
c       Flux due to Chemical Dissolution (or Adsorption)

c     Method of Lines is used to convert the original form of pde
c     dC/dt = d2C/dx2 to an ode dC/dt = f(C) by differencing d2C/dx2.
c     This increases the number of effective species to NSp*NMesh
c     since concentrations at each mesh point are linked to the
c     concentrations in the rest of mesh points.


c  MODIFICATION NOTES:


c  DEVELOPMENT STAGE

c  950202  Passivity and limiting current routine modified for numerical stab.
c  950114  Precipitation reactions re-added. jhc
c  950107  Function Steady added to test for steady state. jhc
c  950102  Tip flux is moved from the flux equation to the source term. jhc
c  941227  Chemical dissolutions now depend on reactant concentrations.  jhc
c  940827  Gap added to simulate Alkire's circular speciemen.  jhc
c  11/8/93 Major update commenced w/ SEAC 1.1.
c  8/22/92 Conception and creation of SEAC from RADICAL 1.2.

c  Author : JOHN H. CHUN
c           MASSACHUSETTS INSTITUTE of TECHNOLOGY
c           185 ALBANY ST. Rm NW22-197
c           (617) 253-3247/0807 fax
c           johnchun@mit.edu

c  Adapted from RADICAL 1.2, another fine product of John Chun.
c  Thanks to Dr. Allan Hindmarsh of LLNL for help with lsode.
c  Thanks to Dr. Steve Boerigter for many helpful discussions.
c  Thanks to Dr. Yutaka Watanabe for his expert opinions.
c  Thanks to Prof. Tetsuo Shoji for his guidance.
c  Thanks to Hiu Au for many friendly helps.

c***********************************************************************

      include 'seac1.1.fi'
      logical*4 IsItBatch,OpenInputFile

      Batch=IsItBatch()      !is it Batch Mode?
10    if (.not.OpenInputFile()) goto 999  !no input file 

      call Initialize        !reads input parameters and prepare them
      call CalcCrevice       !calculate crevice chemstry
      call Terminate         !writes run statistics and finish job

      if (Batch) goto 10     !continue batch process

999   print *
      print *,'**************** END OF A SEAC RUN *****************'
      print *,'**************** HAVE A GREAT DAY! *****************'
      print *
      print *,'Please hit return to close this window......'
      print *,Bel,Bel,Bel
      close(IBF)             !close batch file

      pause
      stop
      end  ! of SEAC


      logical*4 function IsItBatch()

c***********************************************************************
c     Version:        seac 1.1          3 Feb 1994
c     Author :        John H. Chun
c***********************************************************************
c     Called by seac

c     see if the batch file 'seac.bat' is present in the current folder.
c     If it is, abandon interactive mode and jump into batch mode.
c***********************************************************************

      include 'seac1.1.fi'
      character Confirm*10

      print *
      print 3
3     format (' ***** Welcome to SEAC - ',
     +'Simulation of Environmentally-Assisted Cracking *****')
      print *

      open (IBF,file=BatchFile,status='old',err=10)  !open input data file
      print 1
1     format('Are you sure you want to run in the '
     +       'batch mode? [yes/no] '$)
      read 2,Confirm
2     format(a10)

10    print *
      if ((Confirm(1:1).eq.'y').or.(Confirm(1:1).eq.'Y')) then
        IsItBatch=.true.
        print *,'Beginning Batch Mode Calculations...'
      else
        IsItBatch=.false.
        print *,'Running User Interactive Mode...'
      endif
      print *

      return
      end  ! of IsItBatch


      logical*4 function OpenInputFile()

c***********************************************************************
c     Version:        seac 1.1          3 Feb 1994
c     Author :        John H. Chun
c***********************************************************************
c     Called by seac

c     opens input file
c***********************************************************************

      include 'seac1.1.fi'

2     if (Batch) then
10      read (IBF,'(a80)',end=20) InFile
        if ((InFile(1:1).eq.' ').or.(InFile(1:1).eq.'#')) goto 10  !skip
      else
        type 1
1       format('Please type the input file name: '$)
        accept '(a80)',InFile
      endif

      open (IIF,file=InFile,status='old',err=30)  !open input data file
      print *
      print *,'Reading from the input file:   ',InFile
      print *
      OpenInputFile=.true.
      return

20    OpenInputFile=.false.
      return

30    print *
      print '(x,a29,a49)','Error opening the input file ',InFile
      if (Batch) then
        print *,"Let's skip to the next one.",Bel
      else
        print *,"Please try again.",Bel
      endif
      print *
      goto 2
      
      end  ! of OpenInputFile


      subroutine Initialize

c***********************************************************************
c     Version:        SEAC 1.1                    11/5/93
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by seac
c     CALLS ReadIn, SetUpRx, WriteOut

c     opens input and output files and reads input data.
c     prepares input data for calculation and prints them in output file.
c***********************************************************************

      include 'seac1.1.fi'

      Time1=secnds(0.0)      !start clock to measure total execution time

      call InitializeData    !initialize global parameters to default values
      call ReadIn            !read crevice data
      call InputProcessor    !process input reqiring immediate attention
      call CheckInput        !validate input data
      call SetUpRx           !prepare reactions for DifEq
      open(IOF,file=OutFile,status='unknown') !output file
      if (PlotOut) open(IPF,file=PlotFile,status='unknown') !plot file
      call X2MeshPoints      !nodalize x-axis with mesh points
      call WriteOut          !writes input parameters to output file
      call Adjust            !adjust data for cal

      return
      end  ! of Initialize


      subroutine InitializeData

c***********************************************************************
c     Version:        SEAC 1.1                    11/5/93
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Set default values.
c***********************************************************************

      include 'seac1.1.fi'

c     ASSIGN DEFAULT VALUES

      OutFile='SEAC.OUT'
      PlotFile='SEAC.PLOT'

      NMesh=0
      NMeshIn=0
      NMeshEx=0
      WidthMouth=0.d0
      WidthTip=0.d0
      Width(1)=0.d0         !test for tapered width
      Gap=0.d0
      WallSide=2.d0
      DifLayer=0.d0
      Diameter=0.d0
      RefDist=0.d0
      RefMesh=-1
      EndFlux=1.d0

      TInit=0.d0
      TFirst=0.d0
      TStep=10.d0

      TempRef=298.15d0
      Temp=298.15d0
      Density=1.d-99
      Viscosity=1.d-99
      PotMet=1.d-99
      PotInit=1.d-99
      PotFinal=1.d-99
      ChargeDens=0.d0

      NewtMax=500
      NewtXTol=1.d-10
      NewtFTol=1.d-20

      do j=1,IRX
        PotStd(j)=1.d-99
        PotPas(j)=1.d-99
        IPas(j)=1.d-99
        PotPit(j)=1.d-99
        ILimit(j)=1.d-99
        PotPasCoef(j)=1.d-99
        PotPitCoef(j)=1.d-99
      enddo

      MeshSpace=0
      PlotOut=.true.
      do i=1,IDE
        Debug(i)=.false.
      enddo
      Ppb=.false.
      WritePara=.true.
      CalcLimit=.true.
      SteadyMin=0.d0

      GammaAvg=0.d0
      NeutAvg=0.d0
      Velocity=0.d0
      do i=0,IPO
        GammaCoef(i)=1.d-99
        NeutCoef(i)=1.d-99
        VelCoef(i)=1.d-99
      enddo

      ATol=1.d-15
      RTol=1.d-6
      MF=25
      ITol=1
      IState=1
      ITask=4
      IOpt=0

      do m=1,IMP
        IExtNet(m)=0.d0
        dIEdP(m)=0.d0
        do i=1,ISP
          FluxE(i,m)=0.d0
        enddo  !through ISP
      enddo  !through IMP

      NARx=0

      end  ! of InitializeData


      subroutine ReadIn

c***********************************************************************
c     Version:        SEAC 1.1                    1/14/95
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by Initialize

c     reads the reaction matrix and reaction rate constants.
c***********************************************************************

      include 'seac1.1.fi'
      integer*4 i,j
      logical*4 Error
      character*16 TheName

      namelist /FileName/ OutFile,PlotFile
      namelist /Time/ TInit,TFinal,TFirst,TStep
      namelist /Dimension/ Length,Width,WidthMouth,WidthTip,WallSide,
     +                     Height,DifLayer,NMesh,NMeshIn,NMeshEx,
     +                     EndFlux,RefDist,RefMesh,Diameter,Gap
      namelist /State/ TempRef,Temp,Density,PotMet,PotInit,PotFinal,
     +                 Viscosity,ChargeDens
      namelist /Control/ PlotOut,MeshSpace,Ppb,WritePara,Debug,NewtMax,
     +                   NewtXTol,NewtFTol,CalcLimit,SteadyMin,
     +                   ITask,IOpt,RTol,ITol,ATol,IWork,MF
      namelist /Radiation/ GammaCoef,NeutCoef,GammaAvg,NeutAvg
      namelist /Flow/ Velocity,VelCoef

      rewind(IIF)            !always rewind just to make sure     
      read(IIF,110) TitleLine
110   format(a80)
      rewind(IIF)     
      read(IIF,nml=FileName) !read OutFile, PlotFile Names
      rewind(IIF)
      read(IIF,nml=Time)
      rewind(IIF)
      read(IIF,nml=Dimension)
      rewind(IIF)
      read(IIF,nml=State,end=202)
202   rewind(IIF)
      read(IIF,nml=Control,end=109)
109   rewind(IIF)
      read(IIF,nml=Radiation,end=113)
113   rewind(IIF)
      read(IIF,nml=Flow,end=111)

c  read init conc, charge, diffusion coefficients, molecular weights,
c  gamma, neutron g-values

111   TheName='$Species'
      call FindLine(TheName,IIF,Error)
      if (Error) then
        print *,'Species information missing in the input file.'
        print *,'Program terminated.',Bel,Bel
        pause
        stop
      endif
      NSp=NofLine(IIF)       !determine NSp
      call FindLine(TheName,IIF,Error)  !reset the input pointer
      do i=1,NSp
        read(IIF,220) SpName(i),ConcBulk(i),Charge(i),Diff(i),
     +                MolWt(i),GGamma(i),GNeut(i)
      enddo
220   format(x,a8,x,e10.3,x,i3,5(x,e10.3))

c  read chemical reaction set

      TheName='$Reaction'
      call ReadRx(TheName,1,i,NRx)

c  read precipitation reaction set

      NPRx1=NRx+1            !beginning index of NPRx
      TheName='$PptReaction'
      call ReadRx(TheName,NPRx1,NPRx2,NPRx)

c  read tip electrochemical reaction set

      NTRx1=NPRx2+1          !beginning index of NTRx
      TheName='$TipReaction'
      call ReadRx(TheName,NTRx1,NTRx2,NTRx)

c  read wall electrochemical reaction set

      NWRx1=NTRx2+1          !beginning index of NWRx
      TheName='$WallReaction'
      call ReadRx(TheName,NWRx1,NWRx2,NWRx)

c  read external surface electrochemical reaction set

      NSRx1=NWRx2+1          !beginning index of NSRx
      if ((NMeshEx.gt.0).and.(DifLayer.gt.0.d0)) then
        TheName='$SurfReaction'
        call ReadRx(TheName,NSRx1,NSRx2,NSRx)
      else
        NSRx=0
        NSRx2=NWRx2
      endif

c  read tip dissolution reaction set

      NTDis1=NSRx2+1         !beginning index of NTDis
      TheName='$TipDissolution'
      call ReadRx(TheName,NTDis1,NTDis2,NTDis)

c  read wall dissolution reaction set

      NWDis1=NTDis2+1        !beginning index of NWDis
      TheName='$WallDissolution'
      call ReadRx(TheName,NWDis1,NWDis2,NWDis)

c  read external surface dissolution reaction set

      NSDis1=NWDis2+1        !beginning index of NSDis
      if ((NMeshEx.gt.0).and.(DifLayer.gt.0.d0)) then
        TheName='$SurfDissolution'
        call ReadRx(TheName,NSDis1,NSDis2,NSDis)
      else
        NSDis=0
        NSDis2=NWDis2
      endif

c  read standard potential, passivation potential/current, pitting pot

      if ((NTRx+NWRx+NSRx).eq.0) goto 999
      TheName='$CorrosionData'
      call FindLine(TheName,IIF,Error)
      if (Error) goto 999
      NCorData=NofLine(IIF)      !determine number of corrosion data set
      call FindLine(TheName,IIF,Error)  !reset the input pointer
      do i=1,NCorData
        read(IIF,230) RxName(0),PotStd(0),PotPas(0),IPas(0),PotPit(0),
     +                PotPasCoef(0),PotPitCoef(0),ILimit(0)
        do j=NTRx1,max(NTRx2,NWRx2,NSRx2)
          if (RxName(0).eq.RxName(j)) then
            PotStd(j)=PotStd(0)
            PotPas(j)=PotPas(0)
            IPas(j)=IPas(0)
            PotPit(j)=PotPit(0)
            PotPasCoef(j)=PotPasCoef(0)
            PotPitCoef(j)=PotPitCoef(0)
            ILimit(j)=ILimit(0)
          endif
        enddo
      enddo
230   format(x,a3,x,7d9.2)

999   return
      end  !of ReadIn


      subroutine FindLine(TheName,IIF,Error)

c***********************************************************************
c     Version:        SEAC 1.1                    11/5/93
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by ReadIn

c     Scans input file and looks for 'TheName'.  Once it finds the text,
c     the input pointer points to the next line.     

c     Error is used to notify error condition and to transfer control.
c     Error=-1 if TheName not found.
c***********************************************************************

      integer*4 IIF
      logical*4 Error
      character*16 AName,TheName

      Error=.false.
      rewind(IIF)            !rewind input file pointer

100   read(IIF,110,end=120) AName !look for TheName
      if (AName.eq.TheName) return
      goto 100
110   format(x,a16)

120   Error=.true.           !notify the caller
      
      return
      end  !of FindLine


      subroutine ReadRx(TheName,Init,Final,NofRx)

c***********************************************************************
c     Version:        SEAC 1.1                    3/9/93
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by ReadIn
c***********************************************************************

      include 'seac1.1.fi'
      integer*4 Init,Final,NofRx,NofLine
      logical*4 Error
      character*16 TheName
      character*8 Sp(INR+INP)

      call FindLine(TheName,IIF,Error)
      if (Error) then
        NofRx=0
        Final=Init+NofRx-1
        return
      endif
      NofRx=NofLine(IIF)       !determine number of rx
      Final=Init+NofRx-1
      call FindLine(TheName,IIF,Error)  !reset the input pointer
      do j=Init,Final
        read(IIF,100) RxName(j),(Sp(k),k=1,INR+INP),RateConst(j),EA(j)
        if (RxName(j)//Sp(1).eq.'Sam as Tip') then  !'e' of 'Same' is not read
          if (TheName(6:).eq.'Reaction') then
            Init=NTRx1
            Final=NTRx2
            NofRx=NTRx
          else
            Init=NTDis1
            Final=NTDis2
            NofRx=NTDis
          endif
          return
        endif
        call MapSpecies(NSp,INR,INP,IR,IP,Sp,SpName,j,IRX)
      enddo  !for reactions
100   format(x,a3,x,3a8,x,4a8,2d9.2)
      
      NARx=NARx+NofRx          !number of all reactions

      return
      end  !of ReadRx


      integer*4 function NofLine(IIF)

c***********************************************************************
c     Version:        SEAC 1.1                    11/5/93
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by ReadIn
c     returns the number of lines just before '$' is reached.
c***********************************************************************

      integer*4 IIF
      character*80 Line

      NofLine=0
      read(IIF,110) Line
      do while (index(Line,'$').ne.2)  !increment until line is '$END'
        NofLine=NofLine+1
        read(IIF,110) Line
      enddo
      
110   format(a80)

      return
      end  !of NofLine


      subroutine MapSpecies(NSp,INR,INP,IR,IP,Sp,SpName,j,IRX)

c***********************************************************************
c     Version:        SEAC 1.1                    11/5/93
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by ReadIn
c     maps species names to numbers for use in DifEq.
c***********************************************************************

      integer*4 NSp,INR,INP,IR(IRX,INR),IP(IRX,INP),j,IRX,i,k
      character*8 Sp(INR+INP),SpName(0:NSp)

      do k=1,INR             !reactants
        IR(j,k)=0
        do i=1,NSp
          if (Sp(k).eq.SpName(i)) IR(j,k)=i
        enddo  !for species
      enddo  !for reactants

      do k=1,INP             !products
        IP(j,k)=0
        do i=1,NSp
          if (Sp(k+INR).eq.SpName(i)) IP(j,k)=i
        enddo  !for species
      enddo  !for products

      return
      end  !of MapSpecies


      subroutine InputProcessor

c***********************************************************************
c     Version:        SEAC 1.1                    3/13/94
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by Initialize
c     Process input data for self consistency.
c***********************************************************************

      include 'seac1.1.fi'
      real*8 Tm,Tc

c  set smallest meaningful concentrations

      Small=ATol*10.d0

c  determine the maximum order

      MaxOrdGamma=-1
      MaxOrdNeut=-1
      MaxOrdVel=-1
      do i=0,IPO
        if (GammaCoef(i).ne.1.d-99) MaxOrdGamma=i
        if (NeutCoef(i).ne.1.d-99) MaxOrdNeut=i
        if (VelCoef(i).ne.1.d-99) MaxOrdVel=i
      enddo

      if (WidthMouth.eq.0.d0) WidthMouth=Width(0)

c  check for external surface data consistency

      if (NMesh.gt.0) NMeshIn=NMesh
      NMesh=NMeshIn+NMeshEx

      if ((NMeshEx.le.0).or.(DifLayer.le.0.d0).or.(Height.le.0.d0)) then
        NMeshEx=0
        NMesh=NMeshIn
        Height=0.d0
        DifLayer=0.d0
        NSRx=0
        NSRx2=NSRx1-1
        NSDis=0
        NSDis2=NSDis1-1
      endif

      if ((NMeshIn.le.0).or.(Length.le.0.d0).or.(WidthMouth.le.0.d0))
     +  then
        NMeshIn=0
        NMesh=NMeshEx
        Length=0.d0
        WidthMouth=0.d0
        NTRx=0
        NTRx2=NTRx1-1
        NWRx=0
        NWRx2=NWRx1-1
        NTDis=0
        NTDis2=NTDis1-1
        NWDis=0
        NWDis2=NWDis1-1
        print '(a50,a)',
     +    '*** No Internal Calculations - Not Enough Data',Bel
      endif

c  water density
c  obtained from Masanori Takahashi,1991, valid 0-360 C within 0.1%

      if (Density.eq.1.d-99) then
        Tm=Temp-273.15d0-3.986d0
        Tc=647.3d0
        Density=(0.999973d0-Tm*Tm*1.d-6*(1.74224d0+482.502d0
     +         /(Tm+77.861d0)))*(1.d0-0.253d0*(Temp/Tc)**17.4122d0)
      endif

c  water viscosity
c  from D.D. Macdonald, 1992

      if (Viscosity.eq.1.d-99) Viscosity=exp(-1103.164d0/Temp
     +  +457155.3d0/Temp/Temp-6.140834d0)

c  Reynolds number

      Re=Velocity*Diameter/Viscosity

c  misc fixes

      if (RefMesh.lt.0) RefMesh=NMeshEx
      if (TFirst.le.0.d0) TFirst=TStep
      if (PotInit.eq.PotFinal) PotFinal=1.d-99
      if ((PotMet.eq.1.d-99).and.(PotInit.ne.1.d-99)) PotMet=PotInit

      return
      end  !of InputProcessor


      subroutine CheckInput

c***********************************************************************
c     Version:        SEAC 1.1                    11/5/93
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by Initialize

c     Validates input data to make sure everything is within limits.
c***********************************************************************

      include 'seac1.1.fi'
      real*8 Mass
      logical*4 Error

      Error=.false.
      print *

      if (NSp.gt.ISP) then
        print *,Bel,
     +  '### INPUT DATA ERROR - Too many species!  NSp must be <=',ISP
        Error=.true.
      endif

      if (NARx.gt.IRX) then
        print *,'### INPUT DATA ERROR - Too many reactions! ',
     +          'NRx+NTRx+NWRx must be <=',IRx,Bel
        Error=.true.
      endif        

      if (NMesh.lt.2) then
        print *,'### INPUT DATA ERROR - Too few mesh points! ',
     +          'NMeshIn + NMeshEx must be >= 2',Bel
        Error=.true.
      endif

      if (NMesh.gt.IMP) then
        print *,'### INPUT DATA ERROR - Too many mesh points! ',
     +          'NMeshIn + NMeshEx must be <=',IMP,Bel
        Error=.true.
      endif

      if (TInit.lt.0.d0) then
        print *,'### INPUT DATA ERROR - TInit must be >= 0.0 sec!'
        Error=.true.
      endif

      if (TFinal.lt.TInit) then
        print *,'### INPUT DATA ERROR - TFinal must be > TInit!'
        Error=.true.
      endif

      if ((TStep.ge.0.d0).and.(TStep.le.1.d0)) then
        print *,'### INPUT DATA ERROR - TStep cannot be between',
     +          ' 0 and 1!'
        Error=.true.
      endif

      if (Length.lt.0.d0) then
        print *,'### INPUT DATA ERROR - Length must be > 0.0 cm!'
        Error=.true.
      endif

      if (WidthMouth.lt.0.d0) then
        print *,'### INPUT DATA ERROR - Width must be > 0.0 cm!'
        Error=.true.
      endif

      if ((RefDist.gt.0.d0).and.(RefDist.le.DifLayer)) then
        print *,'### INPUT DATA ERROR - RefDist must be > DifLayer!'
        Error=.true.
      endif

      if (TempRef.lt.0.d0) then
        print *,'### INPUT DATA ERROR - TempRef must be > 0.0 k!'
        Error=.true.
      endif

      if (Temp.lt.0.d0) then
        print *,'### INPUT DATA ERROR - Temp must be > 0.0 k!'
        Error=.true.
      endif

      if (Density.lt.0.d0) then
        print *,'### INPUT DATA ERROR - Density must be > 0.0 g/cc!'
        Error=.true.
      endif

      if (GammaAvg.lt.0.d0) then
        print *,'### INPUT DATA ERROR - GammaAvg must be >= 0.0 rad/s'
        Error=.true.
      endif

      if (NeutAvg.lt.0.d0) then
        print *,'### INPUT DATA ERROR - NeutAvg must be >= 0.0 rad/s'
        Error=.true.
      endif

      if ((MF.ne.24).and.(MF.ne.25)) then
        print *,'### INPUT DATA ERROR - MF must be 24 or 25'
        Error=.true.
      endif

      do i=1,NSp             !check for all species data
        if (ConcBulk(i).lt.0.d0) then
          print *,'### INPUT DATA ERROR - ConcBulk for ',SpName(i),
     +            ' must be >= 0.0 mol/L'
          Error=.true.
        endif
        if (Diff(i).lt.0.d0) then
          print *,'### INPUT DATA ERROR - Diff. coef. for ',SpName(i),
     +            ' must be >= 0.0 cm2/s'
          Error=.true.
        endif
        if (MolWt(i).lt.0.d0) then
          print *,'### INPUT DATA ERROR - Molecular weight for ',
     +            SpName(i),' must be >= 0.0'
          Error=.true.
        endif
        if (GGamma(i).lt.0.d0) then
          print *,'### INPUT DATA ERROR - Gamma g-value for ',SpName(i),
     +            ' must be >= 0.0 /100eV'
          Error=.true.
        endif
        if (GNeut(i).lt.0.d0) then
          print *,'### INPUT DATA ERROR - Neutron g-value for ',
     +             SpName(i),' must be >= 0.0 /100eV'
          Error=.true.
        endif
      enddo  !through all species

        
      do j=1,NARx            !check for all reactions
        if (RateConst(j).lt.0.d0) print 110,RxName(j)
        if (EA(j).lt.0.d0) print 150,RxName(j)
      enddo  !through all reactions

110   format('*** INPUT DATA WARNING - Rate constant ',
     + 'in reaction ',a3,/25x'must be >= 0.0')
150   format('*** INPUT DATA WARNING - Activation energy (or partition
     + coef. or'/25x'equilibrium const.) in reaction ',a3,
     + ' must be >= 0.0')

c     Mass and charge balance for CHEMICAL REACTIONS

      Charge(0)=0
      MolWt(0)=0.d0
      j=min(NTDis,NWDis,NSDis)
      if (NRx.gt.0) j=1
      do while (j.le.NARx)         !do not check electrochemcial reactions
        NetCharge(j)=0
        Mass=0.d0
        do i=1,INP
          Mass=Mass+MolWt(IP(j,i))
          NetCharge(j)=NetCharge(j)+Charge(IP(j,i))
        enddo
        do i=1,INR
          Mass=Mass-MolWt(IR(j,i))
          NetCharge(j)=NetCharge(j)-Charge(IR(j,i))
        enddo
        if (mod(abs(dnint(Mass)),18.d0).ne.0) print 174,RxName(j),Mass
        if (NetCharge(j).ne.0) print 173,RxName(j)
        j=j+1
        if (j-1.eq.NPRx2) j=NTDis1   !skip electrochemcial reactions
      enddo

174   format('*** INPUT DATA WARNING - mass is not balanced in reaction'
     +       ,x,a3,x,f7.3)
173   format('*** INPUT DATA WARNING - Net charge for reaction'
     +       ,x,a3,' is not zero.')

c     Net Charge Transfer for Wall and Tip Electrochemical Rx

      j1=NTRx1
      j2=NSRx2
      if (NMeshIn.eq.0) j1=NSRx1
      if (NMeshEx.eq.0) j2=max(NTRx2,NWRx2)
      do j=j1,j2       !check for tip, wall, surface rx
        NetCharge(j)=0
        do i=1,INP
          NetCharge(j)=NetCharge(j)+Charge(IP(j,i))
        enddo
        do i=1,INR
          NetCharge(j)=NetCharge(j)-Charge(IR(j,i))
        enddo
        if (NetCharge(j).eq.0) then
          print 163,RxName(j)
          Error=.true.
        endif
      enddo

163   format('### INPUT DATA ERROR - Net charge for tip and wall '
     + 'electrochemical reaction ',a3,/23x'must be non-zero.')

999   if (Error) then
        print *
        print *
        print *,'### EXECUTION TERMINATED DUE TO INPUT ERROR!'
        print *
        print *,'--- Please correct input data and try again.'
        print *
        print *,'Please press return to close this window.',Bel
        pause
        stop
      endif

      return
      end  !of CheckInput


      subroutine SetUpRx

c***********************************************************************
c     Version:        SEAC 1.1                    11/5/93
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by Initialize

c     Sets up chemical reaction equations.
c     Defines dc/dt with all reactions considering whether
c     the species is produced or consumed, first or second order.

c     Koef :reaction coefficient; + for product, - for reactant.
c                                 1 for first order, 2 for second order.
c***********************************************************************

      include 'seac1.1.fi'
      dimension Ko(IRX,0:ISP)

c  initialize coefficients to zero

      do i=0,NSp
        do j=1,NARx
          Koef(j,i)=0
          Ko(j,i)=0
        enddo
      enddo

c  set up coefficient matricies(Koef) for DifEq

      do j=1,NARx

c       check for SECOND ORDER REACTANTS ******************************

        if (((IR(j,1).eq.IR(j,2)).or.(IR(j,2).eq.IR(j,3)))
     +     .and.(IR(j,2).ne.0)) then
          Koef(j,IR(j,2))=-2
          Ko(j,IR(j,2))=-2
        endif

c       check for FIRST ORDER REACTANTS *******************************

        do k=1,INR
          if ((IR(j,k).ne.0).and.(Koef(j,IR(j,k)).ne.-2)) then
            Koef(j,IR(j,k))=-1
            Ko(j,IR(j,k))=-1
          endif
        enddo

c       check for SECOND ORDER PRODUCTS *******************************

        if (((IP(j,1).eq.IP(j,2)).or.(IP(j,2).eq.IP(j,3)))
     +    .and.(IP(j,2).ne.0)) Koef(j,IP(j,2))=2
        if (((IP(j,2).eq.IP(j,3)).or.(IP(j,3).eq.IP(j,4)))
     +    .and.(IP(j,3).ne.0)) Koef(j,IP(j,3))=2

c       fill up the product matrix for FIRST ORDER PRODUCTS **********

        do k=1,INP
          if ((IP(j,k).ne.0).and.(Koef(j,IP(j,k)).ne.2))
     +       Koef(j,IP(j,k))=1
        enddo
      enddo

c     remove catalytic reactants and products ************************

      do i=1,NSp
        do j=1,NARx
          if ((Koef(j,i).ne.Ko(j,i)).and.(Ko(j,i).ne.0))
     +       Koef(j,i)=Koef(j,i)+Ko(j,i)
        enddo
      enddo

      do m=0,NMesh       !assign zero-order rx concentrations
        Conc(0,m)=1.d0
      enddo
      SpName(0)='        '
            
c     adjust chemical rx rate constants using Arrhenius law.

      if (Temp.ne.TempRef) then
        do j=1,NRx
          RateConst(j)=RateConst(j)*dexp(EA(j)*1.d3/GasConst*
     +                 (1.d0/TempRef-1.d0/Temp))
        enddo
      endif

      return
      end  !of SetUpRx


      subroutine WriteOut

c***********************************************************************
c     Version:        SEAC 1.1                    1/14/95
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by Initialize

c     Writes title, some input info and rx matrix.
c***********************************************************************

      include 'seac1.1.fi'
      logical*4 Duplicate
      character Today*9,Now*8,ConcUnit*5

      call date(Today)
      call time(Now)
      write(IOF,120) Today,Now
120   format(
     + 12x,   '_______________________________________________________'
     +,/ ,12x,'|                                                     |'
     +,/ ,12x,'|                      SEAC 1.1                       |'
     +,/ ,12x,'|   Simulation of Environmentally-Assisted Cracking   |'
     +,/ ,12x,'|                                                     |'
     +,/ ,12x,'|',17x,a9,2x,a8,17x,'|'
     +,/ ,12x,'_______________________________________________________')
      write(IOF,140)
      write(IOF,'(/34x,a)') 'INPUT DATA'
      write(IOF,140)
      write(IOF,'(a80)') TitleLine
      call DifEqStamp(IOF)
      write(IOF,134) InFile,OutFile
      if (PlotOut) write(IOF,136) PlotFile
140   format(/79('_')/)
134   format(/7x,'Input File Name                 = 'a35
     +       /7x,'Output File Name                = 'a35)
136   format( 7x,'Plot File Name                  = 'a35)

      write(IOF,110) TInit,TFinal,TStep
      if (TFirst.ne.TStep) write(IOF,116) TFirst
110   format(/7x,'Initial Time                    = '1pe14.6' sec'
     +       /7x,'Final Time                      = '  e14.6' sec'
     +       /7x,'Time Increment                  = '  e14.6' sec')
116   format(/7x,'First Time Step                 = '1pe14.6' sec')

      if (NMeshIn.gt.0) then
        write(IOF,105) Length
        if (WidthTip.eq.Width(1)) then
          write(IOF,107) WidthMouth
        else
          write(IOF,108) WidthMouth,WidthTip
        endif
        if (Gap.gt.0.d0) write(IOF,122) Gap
        write(IOF,111) NMeshIn
      endif
      if (NMeshEx.gt.0) write(IOF,109) Height,DifLayer,NMeshEx
      if (WallSide.eq.1.d0) write(IOF,117)
      if (EndFlux.eq.0.d0) write(IOF,114)
      if (Diameter.gt.0.d0) write(IOF,119) Diameter

105   format(/7x,'Crevice Length                  = '1pe14.6' cm')
107   format( 7x,'Crevice Width                   = '1pe14.6' cm')
108   format( 7x,'Crevice Width at Mouth          = '1pe14.6' cm'
     +       /7x,'Crevice Width at Tip            = '  e14.6' cm')
122   format( 7x,'Crevice Gap                     = '1pe14.6' cm')
111   format( 7x,'Mesh Points for Crevice         = 'i8)
109   format( 7x,'External Surface Height         = '1pe14.6' cm'
     +       /7x,'Diffusion Layer Thickness       = '  e14.6' cm'
     +       /7x,'Mesh Points for Extrnal Surface = 'i8)
119   format(/7x,'Pipe Diameter                   = '  f14.6' cm')
117   format( 7x,'Crevice is Single Sided.')
114   format( 7x,'End Flux is Set to Zero.')

      write(IOF,106) Temp,TempRef,Density,Viscosity,ChargeDens
106   format(/7x,'Temperature                     = 'f14.6' Kelvin'
     +       /7x,'Reference Temperature           = 'f14.6' Kelvin'
     +       /7x,'Water Density                   = 'f14.6' g/cc'
     +       /7x,'Water Viscosity                 = 'f14.6' cm2/s'
     +       /7x,'Charge Density of Metal         = 'f14.6' Coul/cm3'/)
      if (PotMet+PotInit+PotFinal.eq.3.d-99) then
        write(IOF,112)
      elseif (PotFinal.eq.1.d-99) then
        write(IOF,113) PotMet
      else
        write(IOF,118) PotInit,PotFinal
      endif
112   format( 7x,'Metal at Free Corrosion Potential.')
113   format( 7x,'Metal Potential                 = 'f14.6' Volts')
118   format( 7x,'Potentiostat Polariz Scan from'f10.6' V to'f10.6' V')

      if (.not.CalcLimit) write(IOF,121)
      if (MeshSpace.gt.0) write(IOF,115) MeshSpace
      if (.not.PlotOut) write(IOF,256)
      if (SteadyMin.gt.0.d0) write(IOF,257),SteadyMin
121   format( 7x,'Limiting Currents are Not Calculated.')
115   format( 7x,'Nonuniform Mesh Spacing Profile = 'i8)
256   format( 7x,'No Plot File is Generated.')
257   format( 7x,'Steady state threshold          = '1pe14.6)

      if (GammaAvg.gt.0.d0) then
        write(IOF,300) GammaAvg
        write(IOF,310) (i,GammaCoef(i),i=0,MaxOrdGamma)
      endif
300   format(/7x,'Average Gamma Dose Rate         = '1pe14.6' Rad/s')
310   format( 7X,'Gamma Coefficient',i3,12x,'= ',1pe14.6)

      if (NeutAvg.gt.0.d0) then
        write(IOF,320) NeutAvg
        write(IOF,330) (i,NeutCoef(i),i=0,MaxOrdNeut)
      endif
320   format(/7x,'Average Neutron Dose Rate       = '1pe14.6' Rad/s')
330   format( 7X,'Neutron Coefficient',i3,10x,'= ',1pe14.6)

      if (Velocity.ne.0.d0) then
        write(IOF,500) Velocity
        write(IOF,520) (i,VelCoef(i),i=0,MaxOrdVel)
      endif
500   format(/7x,'Flow Velocity at Crevice Mouth  = 'f14.6' cm/s')
520   format( 7X,'Velocity Coefficient',i3,9x,'= ',1pe14.6)

      write(IOF,250) NewtXTol,NewtFTol,NewtMax,ATol,RTol,MF
250   format(/7x,"Newton's X     Tolerance        = "1pe14.6
     +       /7x,'         Func  Tolerance        = 'e14.6
     +       /7x,'         Max   Iteration        = 'i8
     +       /7x,'lsode Absolute Tolerance        = 'e14.6
     +       /7x,'      Relative Tolerance        = 'e14.6
     +       /7x,'      Method Flag               = 'i8)

      if (ppb) then
        write(IOF,190) ' [ppb] '
      else
        write(IOF,190) '[mol/L]'
      endif
190   format(//x,'Species   Initial  Charge  Diffu'
     +           "       Mol      Gamma     Neutron"
     +        /x,' Name      Conc            Coef'
     +           '        Wt      G-Values   G-Values'
     +        /x,'          ',a7,'         [cm2/s]'
     +           '    [g/mol]  [#/100eV]  [#/100eV]'
     +        /x,'-------- ---------- --- ----------'
     +           ' ---------- ---------- ----------'/)
      do i=1,NSp
        write(IOF,210) SpName(i),ConcBulk(i),Charge(i),Diff(i),
     +                 MolWt(i),GGamma(i),GNeut(i)
      enddo
210   format(x,a8,1pe11.4,x,i3,e11.4,4(0pf11.6))

c     write chemical and electrochemical reactions ******************

      if (NARx.gt.0) then
        write(IOF,292)
        write(IOF,9)
        write(IOF,140)
      endif

292   format(////79('_'),/)
9     format(/20X,'Chemical & Electrochemical Reactions')

c     write CHEMICAL REACTIONS **************************************

      if (NRx.gt.0) then
        write(IOF,10)
        write(IOF,11) (RxName(j),(SpName(IR(j,k)),k=1,INR),
     +                (SpName(IP(j,k)),k=1,INP),RateConst(j),EA(j),
     +                 j=1,NRx)
      endif

10    format(//64x,'Rate',2x,'Activation'
     +        /21x,'Chemical Reactions',23x,'Constant',2x,'Energy'
     +        /70x,'(kJ/mol-k)'/)
11    format(x,a3,x,3a8,'>',4a8,1pe9.2,e9.2)

c     write PRECIPITATION REACTIONS ********************************

      if (NPRx.gt.0) then
        write(IOF,410)
        write(IOF,11) (RxName(j),(SpName(IR(j,k)),k=1,INR),
     +                (SpName(IP(j,k)),k=1,INP),RateConst(j),Ksp(j),
     +                 j=NPRx1,NPRx2)
      endif

410   format(//18x,'Precipitation Reactions',23x,'Rate',6x,'Ksp'/)

c     write crevice TIP ELECTROCHEMICAL REACTIONS *****************

      if (NTRx.gt.0) then
        write(IOF,12)
        write(IOF,11) (RxName(j),(SpName(IR(j,k)),k=1,INR),
     +                (SpName(IP(j,k)),k=1,INP),IRate(j),EACoef(j),
     +                 j=NTRx1,NTRx2)
      endif

12    format(//64x,'Rate   Partition'
     +        /11x,'Crevice Tip Electrochemical Reactions'
     +         15x,'Const    Coef'/)

c     write WALL ELECTROCHEMICAL REACTIONS ***************************

      if (NWRx.gt.0) then
        if (NWRx2.ne.NTRx2) then
          write(IOF,13)
          write(IOF,11) (RxName(j),(SpName(IR(j,k)),k=1,INR),
     +                  (SpName(IP(j,k)),k=1,INP),IRate(j),EACoef(j),
     +                   j=NWRx1,NWRx2)
        else
          write(IOF,14)
        endif
      endif

13    format(//64x,'Rate   Partition'
     +        /12x,'Crevice Wall Electrochemical Reactions'
     +         13x,'Const    Coef'/)
14    format(//12x,'Crevice Wall Electrochemical Reactions'
     +        /12x,'Same as Tip Electrochemical Reactions'/)

c     write SURFACE ELECTROCHEMICAL REACTIONS ************************

      if (NSRx.gt.0) then
        if (NSRx2.ne.NTRx2) then
          write(IOF,15)
          write(IOF,11) (RxName(j),(SpName(IR(j,k)),k=1,INR),
     +                  (SpName(IP(j,k)),k=1,INP),IRate(j),EACoef(j),
     +                   j=NSRx1,NSRx2)
        else
          write(IOF,16)
        endif
      endif

15    format(//64x,'Rate   Partition'
     +        /10x,'External Surface Electrochemical Reactions'
     +         11x,'Const    Coef'/)
16    format(//10x,'External Surface Electrochemical Reactions'
     +        /12x,'Same as Tip Electrochemical Reactions'/)

c     write CORROSION DATA *******************************************

      if (NCorData.gt.0) then
        write(IOF,211)
        j1=NTRx1
        j2=NSRx2
        if (NMeshIn.eq.0) j1=NSRx1
        if (NMeshEx.eq.0) j2=max(NTRx2,NWRx2)
        do j=j1,j2
          P=PotStd(j)+PotPas(j)+IPas(j)+PotPit(j)
     +     +PotPasCoef(j)+PotPitCoef(j)+ILimit(j)
          if (P.ne.7.d-99) then
            Duplicate=.false.
            do i=j1,j-1
              if (RxName(i).eq.RxName(j)) Duplicate=.true.
            enddo
            if (.not.Duplicate) write(IOF,212)
     +        RxName(j),PotStd(j),PotPas(j),IPas(j),PotPit(j),
     +        PotPasCoef(j),PotPitCoef(j),ILimit(j)
          endif
        enddo
      endif
211   format(//27x,'Corrosion Data'
     +       // 6x,'PotStd   PotPas     IPas     PotPit   '
     +             'PotPas   PotPit    ILimit'
     +        / 6x,'(Volts)  (Volts)   (A/cm2)   (Volts)  '
     +             ' Coef     Coef     (A/cm2)'
     +        / 6x,'-------  -------  ---------  -------  '
     +             '-------  -------  ---------')
212   format(x,a3,2(x,f8.5),x,1pe10.3,3(x,0pf8.5),x,1pe10.3)

c     write TIP CHEMICAL DISSOLUTION REACTIONS ***********************

      if (NTDis.gt.0) then
        write(IOF,17)
        write(IOF,411) (RxName(j),(SpName(IR(j,k)),k=1,INR),
     +                (SpName(IP(j,k)),k=1,INP),JDis(j),
     +                 j=NTDis1,NTDis2)
      endif

17    format(//12x,'Crevice Tip Chemical Dissolutions',19x,'Flux'/)
411   format(x,a3,x,3a8,'>',4a8,1pe9.2)

c     write WALL CHEMICAL DISSOLUTION REACTIONS **********************

      if (NWDis.gt.0) then
        if (NWDis2.ne.NTDis2) then
          write(IOF,18)
          write(IOF,411) (RxName(j),(SpName(IR(j,k)),k=1,INR),
     +                  (SpName(IP(j,k)),k=1,INP),JDis(j),
     +                   j=NWDis1,NWDis2)
        else
          write(IOF,19)
        endif
      endif

18    format(//12x,'Crevice Wall Chemical Dissolutions',18x,'Flux'/)
19    format(//12x,'Crevice Wall Chemical Dissolutions'
     +        /12x,'Same as Tip Chemical Dissolutions'/)

c     write EXTERNAL SURFACE DISSOLUTION REACTIONS ********************

      if (NSDis.gt.0) then
        if (NSDis2.ne.NTDis2) then
          write(IOF,20)
          write(IOF,411) (RxName(j),(SpName(IR(j,k)),k=1,INR),
     +                  (SpName(IP(j,k)),k=1,INP),JDis(j),
     +                   j=NSDis1,NSDis2)
        else
          write(IOF,21)
        endif
      endif

20    format(//12x,'External Surface Chemical Dissolutions',15x,'Flux'/)
21    format(//12x,'External Surface Chemical Dissolutions'
     +        /12x,'Same as Tip Chemical Dissolutions'/)

      write(IOF,292)
      write(IOF,400)
      write(IOF,140)

400   format(/25x,'Concentrations and Potential')

c     write species names to plot file *******************************

      if (PlotOut) then      !write CONCENTRATION RESULTS
        write(IPF,'(a80)') TitleLine
        call DifEqStamp(IPF)
        ConcUnit='mol/L'
        if (Ppb) ConcUnit='ppb'
        write(IPF,'(2(a11,a),99(a19,a))') 'Time [sec]',Tab,'x [cm]',Tab,
     +               (SpName(i)(:index(SpName(i)//' ',' '))//
     +               '['//ConcUnit//']',Tab,i=1,NSp),
     +               'Potential [Volts]',Tab,'Ext Pot [Volts]',Tab,
     +               'Conductivity [S/cm]',Tab,'Tot Current [A/cm2]',
     +               Tab,'Total Charge',Tab,'Wal Current [A]',Tab,
     +               'Anod Curr [A/cm2]',Tab,'Ext Curr [A]',Tab,
     +               'Gamma Dose [Rad/s]',Tab,'Neut Dose [Rad/s]',Tab,
     +               'Velocity [cm/s]'
      endif

      return
      end  !of WriteOut


      subroutine Adjust

c***********************************************************************
c     Version:        SEAC 1.1                    11/5/93
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by Initialize

c     Adjusts input data for calculation.
c***********************************************************************

      include 'seac1.1.fi'
      integer*4 ISign,i,j,j1,j2,m

c  calc Mobility using Nernst-Einstein equation

      RT=GasConst*Temp
      do i=1,NSp
        Mobil(i)=Diff(i)/RT
      enddo

c  modify IRate & EACoef to reduce computation

      FRT=Faraday/RT
      j1=NTRx1
      j2=NSRx2
      if (NMeshIn.eq.0) j1=NSRx1
      if (NMeshEx.eq.0) j2=max(NTRx2,NWRx2)
      do j=j1,j2
        ISign=NetCharge(j)/abs(NetCharge(j))
        EACoef(j)=ISign*FRT*EACoef(j)
        if (EACoef(j).eq.0.d0) EACoef(j)=1.d-100  !patched 941010 jhc
        if (PotStd(j).eq.1.d-99) then
          IRate(j)=ISign*IRate(j)*exp(-EACoef(j)*PotMet)
        else  !mixed potential model
          IRate(j)=ISign*IRate(j)*exp(-EACoef(j)*PotStd(j))
        endif
      enddo

c  calculate limiting current density

      do j=j1,j2                   !if not input, set to a huge number
        ISign=NetCharge(j)/abs(NetCharge(j))
        if (abs(PotPas(j)).le.1.d-99) PotPas(j)=1.d99
        if (abs(PotPit(j)).le.1.d-99) PotPit(j)=1.d99
        IPas(j)=ISign*IPas(j)
        if (ILimit(j).le.1.d-99) then
          ILimit(j)=1.d99
        else
          ILimit(j)=ISign*ILimit(j)
        endif
      enddo
      if (CalcLimit) call LimitingCurrent
      ILimit(0)=1.d99              !ILimit(0) is always large

c  convert From #Species/100eV to moles/L-rad. The formula is:
c  molecule/(100*eV)*mol/(6.02204e23*molecule)*eV/(1.602189e-19*joule)*
c  (0.01*joule/kg)/rad*(1*kg)/(1000*gram)*Density*gram/cc*(1000*cc)/liter
c  =1.036436307514464021e-9*(mol*Density)/(liter*rad)

      GConvert=1.0364363d-9*Density
      if ((GammaAvg.gt.0d0).or.(NeutAvg.gt.0.d0)) then
        call DoseShape       !get doseshape profile for crevice       
        do i=1,NSp
          GNeut(i)=GNeut(i)*GConvert
          GGamma(i)=GGamma(i)*GConvert
        enddo
      endif
      
c  conversion factors between mol/L and ppb for input/output

      do i=1,NSp
        ConvFact(i)=1.0d0
        if (Ppb) then        !convert ppb to mol/L for calc
          ConvFact(i)=MolWt(i)*1.0d6/Density
          ConcBulk(i)=ConcBulk(i)/ConvFact(i)
        endif
      enddo

      if (Velocity.ne.0.d0) call VelShape

c  initial conditions

      Pot(0)=PotMet
      do m=1,NMesh
        do i=1,NSp
          Co(i+NSp*(m-1))=ConcBulk(i)
        enddo
        Pot(m)=PotMet
      enddo

c  Flux at tip is ZERO -- Tip reactions are part of source term

      do i=1,NSp
        FluxH(i,NMesh)=0.d0
      enddo  !through all species

c  convert Henry's constant from [atm/mol frac] to [atm/mol/L] - not used

      do i=1,NSp
        Henry(i)=Henry(i)*Density*1.d3/MolWtH2O
      enddo

c  determine iHydro,iChlor

      iHydro=0
      iChlor=0
      do i=1,NSp
        if (SpName(i).eq.'H+') iHydro=i
        if (SpName(i).eq.'Cl-') iChlor=i
      enddo

c  set output file number for lsode errors

      call xsetun(IOF)

      return
      end  !of Adjust


      subroutine LimitingCurrent

c***********************************************************************
c     Version:        SEAC 1.1                    3/28/94 
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by Adjust
c***********************************************************************

      include 'seac1.1.fi'
      integer*4 k,l
      real*8 JLimit,Flux,Sc

      do j=NSRx1,NSRx2
        JLimit=1.d99              !set to some huge number
        do k=1,INR
          l=IR(j,k)
          if ((l.gt.0).and.(Velocity.eq.0.d0)) then  !flux across diff layer
            Flux=Diff(l)*ConcBulk(l)/DifLayer*1.d-3
          elseif (l.gt.0) then
            Sc=Viscosity/Diff(l)
            Flux=0.0165d0*Diff(l)*ConcBulk(l)*1.d-3*Re**.86d0*Sc**.33d0
     +          /Diameter
          endif
          if (IR(j,1)+IR(j,2)+IR(j,3).gt.0) JLimit=min(JLimit,Flux)
        enddo
        if (ILimit(j).ge.1.d98) ILimit(j)=JLimit*NetCharge(j)*Faraday
      enddo
         
      return
      end  !of LimitingCurrent


      subroutine DoseShape

c***********************************************************************
c     Version:        SEAC 1.1                    11/5/93
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by Adjust

c     Evaluate dose rates in the crevice using the input dose polynomial
c***********************************************************************

      include 'seac1.1.fi'

      do m=NMeshEx,NMesh
      
        if (MaxOrdGamma.ge.0) then !evaluate MaxOrd-th polynomial
          Gamma(m)=GammaCoef(MaxOrdGamma)
        else                 !trap error
          GammaAvg=0.d0
        endif

        if (MaxOrdNeut.ge.0) then
          Neutron(m)=NeutCoef(MaxOrdNeut)
        else
          NeutAvg=0.d0
        endif

        do i=MaxOrdGamma-1,0,-1
          Gamma(m)=x(m)*Gamma(m)+GammaCoef(i)
        enddo
        do i=MaxOrdNeut-1,0,-1
          Neutron(m)=x(m)*Neutron(m)+NeutCoef(i)
        enddo

        Gamma(m)=GammaAvg*Gamma(m)
        Neutron(m)=NeutAvg*Neutron(m)

        if (Gamma(m).lt.0.d0) Gamma(m)=0.d0
        if (Neutron(m).lt.0.d0) Neutron(m)=0.d0

      enddo

      do m=0,NMeshEx-1       !external surface = crack mouth
        Gamma(m)=Gamma(NMeshEx)
        Neutron(m)=Neutron(NMeshEx)
      enddo

      return
      end  !of DoseShape


      subroutine VelShape

c***********************************************************************
c     Version:        SEAC 1.1                    11/5/93
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by Adjust

c     Evaluate dose rates in the crevice using the input dose polynomial
c***********************************************************************

      include 'seac1.1.fi'

      do m=NMeshEx,NMesh
      
        if (MaxOrdVel.ge.0) then !evaluate MaxOrdVel-th polynomial
          VelH(m)=VelCoef(MaxOrdVel)
        else                 !trap error
          Velocity=0.d0
        endif

        do i=MaxOrdVel-1,0,-1
          if (m.lt.NMesh) then
            VelH(m)=(x(m)+x(m+1))/2.d0*VelH(m)+VelCoef(i)
          else
            VelH(m)=x(m)*VelH(m)+VelCoef(i)
          endif
        enddo

        VelH(m)=Velocity*VelH(m)

      enddo

      do m=0,NMeshEx-1     !external surface = diffusion layer = no flow
        VelH(m)=0.d0
      enddo

      return
      end  !of VelShape


      subroutine CalcCrevice

c***********************************************************************
c     Version:        SEAC 1.1                    11/5/93
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by seac
c     CALLS lsode, Y2Conc, CalcFlux, WriteCon

c     Evaluates radiolysis of a component using lsode.
c***********************************************************************

      include 'seac1.1.fi'
      logical*4 Steady       !check if steady state reached
      external DifEq,Jacob
      parameter (MLMax=2*ISP-1,MUMax=MLMax)
      parameter (LRWMax=22+10*IEQ+(2*MLMax+MUMax)*IEQ)
      real*8 RWork(LRWMax)

c      real*8 RWork(*)       !these 4 lines are for Macintosh computer only
c      pointer (pSpare,Spare),(pRWork,RWork)    !define pointer to Array
c      volatile RWork
c      integer*4 GetMemory,Size

      print *
      print *,'Calculating..........'
      print *,'To terminate press Command-period.'
      print *
      print *

c  initialize for lsode

      NStepTot=0
      NFuncTot=0
      NJacobTot=0
      do i=11,18             !reset step counters
        IWork(i)=0
      enddo
      RWork(1)=TFinal
      IState=1

      t     = TInit
      TOut  = TInit

      ML=2*NSp-1
      MU=ML
      NEq=NSp*NMesh
      LRW=22+10*NEq+(2*ML+MU)*NEq
      LIW=20+NEq

      IWork(1)=ML            !ML for banded jacobian
      IWork(2)=MU            !MU for banded jacobian

c  dynamic memory allocation for RWork - for Macintosh computer only

      Size=LRW
c     pRWork=GetMemory(Size) !get that memory space
c     pSpare=GetMemory(1024) !get 8kb spare memory space
c     if ((pRWork.eq.0).or.(pSpare.eq.0)) then  !error
c       print *,'Memory Allocation Error!'
c       print *,'Please increase application memory size by',
c    +           (LRW*8-Size)/1024+10,' kb and try again.'
c       print *,'SEAC Terminated.'
c       pause
c       stop
c     endif
c     call DisposeMemory(pSpare) !release spare memory


c**** MAIN LOOP of CalcCrevice begins here ****************************

280   call lsode(DifEq,NEq,Co,t,TOut,ITol,RTol,ATol,ITask,
     +           IState,IOpt,RWork,LRW,IWork,LIW,Jacob,MF)

      if (IState.eq.-1) then !lsode error - excessive work done
        IState=2
        print *,'*** WORKING HARD at CalcCrevice!  t =',t
        write(IOF,*) '*** WORKING HARD AT CalcCrevice!'
        goto 280             !reset and try again
      else if (IState.lt.0) then
        goto 999
      endif

      call C2Conc            !restore concentrations
      call CalcFluxData      !calc mesh-dep, but pot-indep parameters
      call CalcPotential     !calc pot using Newton's method
      call CalcFlux(1)       !recalculate at t since LSODE overshoots

      if ((t.eq.TInit).and.(PotInit.eq.1.d-99).and.
     +   (PotFinal.ne.1.d-99)) then  !corrosion pot
        PotCor=Pot(RefMesh)
        if (DifLayer.gt.0.d0) PotCor=PotE(RefMesh)
        PotMet=PotCor
        PotInit=PotCor
        PotFinal=PotInit+PotFinal
      endif

      NStep=IWork(11)-NStepTot !# steps taken in this time step
      NFunc=IWork(12)-NFuncTot !# f eval in this time step
      NJacob=IWork(13)-NJacobTot !# jacobian eval in this time step    

      call WriteConc(t)      !write conc to output & plot files

      NStepTot=IWork(11)     !# steps taken so far
      NFuncTot=IWork(12)     !# f eval so far
      NJacobTot=IWork(13)    !# jacobian eval      

c     exit loop upon TFinal

      if (t.ge.TFinal) goto 999

c     exit if steady state reached

      if (Steady(t)) goto 999

c     increment TOut and continue

      if (TStep.lt.0.d0) then !additive time step
        TOut = TOut-TStep
      elseif (TOut.eq.TInit) then !multiplicative time step
        TOut = TFirst
      else
        TOut = TOut*TStep
      endif
      if (TOut.gt.TFinal) TOut=TFinal
      goto 280

c**** continue through time steps until TFinal ************************

999   return
      end  !of CalcCrevice


      logical*4 function Steady(t)

c***********************************************************************
c     Version:        SEAC 1.1                    1/7/95
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by CalcCrevice

c     Check if steady state reached to within SteadyMin.
c***********************************************************************

      include 'seac1.1.fi'
      real*8 t,ConcChange,ConcLast(ISP,IMP),PotChange,PotLast(IMP)
      logical*4 OnceMore

      Steady=.false.
      if (SteadyMin.le.0.d0) return
      if (t.eq.TInit) then
        OnceMore=.false.
        goto 30   !skip first time step
      endif

      do m=1,NMesh
        do i=1,NSp
          ConcChange=(Conc(i,m)-ConcLast(i,m))/ConcLast(i,m)
          if ((dabs(ConcChange).gt.SteadyMin).and.
     +       (Conc(i,m).gt.Small)) goto 20
        enddo
        PotChange=(Pot(m)-PotLast(m))/PotLast(m)
        if (dabs(PotChange).gt.SteadyMin) goto 40
      enddo

      if (.not.OnceMore) then
        OnceMore=.true.
        print *,'Going one more time step to confirm STEADY STATE.'
        print 10,PotLast(NMesh),Pot(NMesh),PotChange,SteadyMin
      else
        Steady=.true.
        print *
        print *,'********* Concentrations Reached STEADY STATE!',Bel,Bel
        write(IOF,*)
        write(IOF,*) '********* Concentrations Reached STEADY STATE!'
      endif
      return

20    print *,SpName(i),' has not reached steady state at m =',m 
      print 10,ConcLast(i,m),Conc(i,m),ConcChange,SteadyMin
      goto 30

40    print *,'Potential has not reached steady state at m =',m 
      print 10,PotLast(m),Pot(m),PotChange,SteadyMin
10    format(' Old =',1pe11.4,'  New =',e11.4,'  Change =',e11.4,
     +       '  Threshold =',e11.4)

30    do m=1,NMesh
        do i=1,NSp
          ConcLast(i,m)=Conc(i,m)
        enddo
        PotLast(m)=Pot(m)
      enddo

      if (OnceMore) then
        print *,'Calculations are UNSTABLE!',Bel
        OnceMore=.false.
      endif

      return
      end  !of Steady


      subroutine WriteConc(t)

c***********************************************************************
c     Version:        SEAC 1.1                    11/5/93
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by CalcCrevice

c     Writes concentrations at each step to output file and plot file.
c     Output format for plotfile has been modified for macintosh-specifi
c     plotting method. Kaleidagraph, a macintosh plotting application is
c     Recommended for plotting concentration results in plotfile.
c***********************************************************************

      include 'seac1.1.fi'
      integer*4 m
      real*8 IWallTot,ISurfTot

c  write to output file

      if (Ppb) then
        print 120,t
        write(IOF,120) t
      else
        print 310,t
        write(IOF,310) t
      endif

310   format(10x,'Concentrations [moles/liter] at Time =',
     +        1pe11.4,' sec')
120   format(10x,'        Concentrations [ppb] at Time =',
     +        1pe11.4,' sec')

      m=NMesh
      Cond(m)=0.d0
      ZC(m)=0.d0
      do i=1,NSp
        Cond(m)=Cond(m)+Charge(i)*Charge(i)*Mobil(i)*Conc(i,m)
        ZC(m)=ZC(m)+Charge(i)*Conc(i,m)
      enddo
      Cond(m)=Cond(m)*Faraday*Faraday*1.d-3  !liter to cm3

      if (NMeshEx.gt.0) then
        PotE(0)=Pot(0)
        m=1
        if (PotMet.eq.1.d-99) m=0
        m1=max((NMeshEx-m)/2,1)
        do while (m.le.NMeshEx)
          print 331, x(m),Pot(m),t
          print 320, (SpName(i),ConvFact(i)*Conc(i,m),i=1,NSp)
          write(IOF,331) x(m),Pot(m),t
          write(IOF,320) (SpName(i),ConvFact(i)*Conc(i,m),i=1,NSp)
          if (WritePara) then
            print 302, Cond(m),ITotalH(m),ZC(m),IWallNet(m),IAnod(m)
            print 305, IExtNet(m)
            print 100, PotE(m)
            write(IOF,302) Cond(m),ITotalH(m),ZC(m),IWallNet(m),IAnod(m)
            write(IOF,305) IExtNet(m)
            write(IOF,100) PotE(m)
          endif
          m=m+m1
          if (m.eq.NMeshEx-1) m=NMeshEx
        enddo
      endif

      if (NMeshIn.gt.0) then
        m=NMeshEx+1
        if ((PotMet.eq.1.d-99).and.(NMeshEx.eq.0)) m=0
        m1=max((NMesh-m)/2,1)
        do while (m.le.NMesh)
          print 331, x(m),Pot(m),t
          print 320, (SpName(i),ConvFact(i)*Conc(i,m),i=1,NSp)
          write(IOF,331) x(m),Pot(m),t
          write(IOF,320) (SpName(i),ConvFact(i)*Conc(i,m),i=1,NSp)
          if (WritePara) then
            print 302, Cond(m),ITotalH(m),ZC(m),IWallNet(m),IAnod(m)
            write(IOF,302) Cond(m),ITotalH(m),ZC(m),IWallNet(m),IAnod(m)
            if ((GammaAvg.gt.0.d0).or.(NeutAvg.gt.0.d0)) then
              print 303, Gamma(m),Neutron(m)
              write(IOF,303) Gamma(m),Neutron(m)
            endif
            if ((Velocity.ne.0.d0).and.(MaxOrdVel.ge.0)) then
              print 304, VelH(m)
              write(IOF,304) VelH(m)
            endif
          endif
          m=m+m1
          if (m.eq.NMesh-1) m=NMesh
        enddo
      endif

      IWallTot=0.d0
      do m=0,NMesh
        if (m.lt.NMeshEx) then
          IWallTot=IWallTot+IWallNet(m)*WallSide
        elseif (m.eq.NMeshEx) then
          ISurfTot=IWallTot+ISurfMouth*WallSide
          IWallTot=ISurfTot+IWallNet(m)-ISurfMouth
        else
          IWallTot=IWallTot+IWallNet(m)
        endif
      enddo
      if (WritePara) then
        print 332, IWallTot,ISurfTot,IWallTot-ISurfTot
        print 330, PotMet,Pot(RefMesh)-PotMet,PotMet-PotCor
        if (ChargeDens.ne.0.d0) print 334, ITip/ChargeDens
        write(IOF,332) IWallTot,ISurfTot,IWallTot-ISurfTot
        write(IOF,330) PotMet,Pot(RefMesh)-PotMet,PotMet-PotCor
        if (ChargeDens.ne.0.d0) write(IOF,334) ITip/ChargeDens
      endif

      print 300, NStep,NFunc,NJacob
      print 290
      write(IOF,300) NStep,NFunc,NJacob
      write(IOF,290)

      if (Batch) then
        print *,'Batch job in progress ... reading from ',InFile
      else
        print *,'Reading from ',InFile
      endif

331   format(/x,'At x =',1pe12.5,' cm,  ECP = ',0pf10.7' V,  ',
     +          'Time =',1pe12.5,' sec'/)
320   format(2(5x,A8,' = ',1pe15.6,' **'))

302   format(/5x,'Conductivity (Siemens/cm)  = ',1pe14.6,
     +       /5x,'Total Current (A/cm2)      = ',e14.6
     +       /5x,'Total Charge (equiv/liter) = ',e14.6,' *'
     +       /5x,'Net Wall Current (A)       = ',e14.6
     +       /5x,'Anodic Curr Dens (A/cm2)   = ',e14.6)
305   format( 5x,'External Current (A)       = ',1pe14.6)
100   format( 5x,'External Potential (V)     = ',f10.6)
303   format( 5x,'Gamma dose rate (Rad/s)    = ',1pe14.6,
     +       /5x,'Neutron dose rate (Rad/s)  = ',e14.6)
304   format( 5x,'Flow velocity (cm/s)       = ',1pe14.6)
332   format(/5x,'Tot Wall/Ext/Int Curnt (A) = ',1pe14.6,2(' /',e13.6))
330   format( 5x,'RefE Pot/IR Drop/Polar (V) = ',f10.6,2(4x,' /',f9.6))
334   format( 5x,'Crack Growth Rate (cm/s)   = ',1pe14.6)
300   format(/5X,'Number of lsode steps/DifEq calls/jacobian evals  = ',
     +       i5,' /',i5,' /',i5)
290   format(79('_'),/)

c  write to plot file ********************************************

      if (PlotOut) then
        do m=0,NMesh
          write(IPF,'(99(1pe10.3,a))') t,Tab,x(m),Tab,
     +                   (ConvFact(i)*Conc(i,m),Tab,i=1,NSp),
     +                   Pot(m),Tab,PotE(m),Tab,Cond(m),Tab,
     +                   ITotalH(m),Tab,ZC(m),Tab,IWallNet(m),Tab,
     +                   IAnod(m),Tab,IExtNet(m),Tab,
     +                   Gamma(m),Tab,Neutron(m),Tab,VelH(m)
        enddo
        write(IPF,*)
      endif
      
      return
      end  !of WriteConc


      subroutine Terminate

c***********************************************************************
c     Version:        SEAC 1.1                    11/5/93
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by seac

c     Terminate writes final run statistics to the output file and
c     closes all files.
c***********************************************************************

      include 'seac1.1.fi'

c     print run statistics

      ET=secnds(Time1)       !elapsed time
      write(IOF,390) IWork(17),IWork(18),IWork(11),IWork(12),IWork(13),
     +               ET
390   format(/7x,'Required RWork Size   = 'i9
     +       /7x,'IWork Size            = 'i9
     +       /7x,'Number of Steps       = 'i9
     +       /7x,'# of Func Evals       = 'i9
     +       /7x,'# of Jacob Evals      = 'i9
     +       /7x,'Job Time              = 'f10.0,' seconds')

      if (IState.gt.0) then  !success
        print 395
        write(IOF,395)
      else
        print 400, IState
        write(IOF,400) IState
      endif
395   format(/' CONCENTRATION PROFILE HAS BEEN EVALUATED SUCCESSFULLY!')
400   format(/' ### LSODE ERROR ... IState =',i3)

      close(IIF)              !close input file
      close(IOF)              !close output file
      if (PlotOut) close(IPF) !close plot file
      print *,Bel,Bel

      return
      end  !of Terminate


      subroutine C2Conc

c***********************************************************************
c     Version:        SEAC 1.1                    11/5/93
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by DifEq, CalcCrevice
c
c     Converts C array to Conc array - part of the method of lines
c
c     Although it is easier for the human to visualize crevice concentrations
c     in terms of species and mesh points, they are just a series of
c     numbers for lsode.
c     Therefore the 2D array Conc(i,m), i.e., the concentration of
c     species i at mesh point m must be remapped to a 1D array for
c     lsode.
c***********************************************************************

      include 'seac1.1.fi'

      do m=1,NMesh
        do i=1,NSp
          Conc(i,m)=Co(i+NSp*(m-1))
        enddo
      enddo

      return
      end  !of C2Conc


      subroutine Jacob(Dum,t,ConcVec,ML,MU,pd,nrowpd)

c***********************************************************************
c     Version:        SEAC 1.1                    11/5/93
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by LSODE
c     This is a dummy subroutine
c***********************************************************************

      return
      end  !of Jacob


      subroutine DifEqStamp(IOF)

c***********************************************************************
c     Version:        SEAC 1.1                    1/14/95
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by WriteOut
c     Prints the latest version number of DifEq on output for future record
c***********************************************************************

      integer*4 IOF

      write(IOF,*) 'DifEq v1.1 on HP, 950202'

      return
      end  !of DifEqStamp


      subroutine DifEq(NEq,t,C,dCdt)

c***********************************************************************
c     Version:        SEAC 1.1                    1/14/95
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by lsode
c     CALLS C2Conc, CalcFluxData, CalcPotential, CalcFlux, CalcEChem
c
c     This is the heart of seac.  The transport equation, dC/dt, is
c     formulated here explicitly.

c     DifEq calculates:
c     dC/dt = -   Divergence
c             +/- Generation & Consumption by Chemcial Rx
c             +/- Generation & Consumption by Wall Electrochemcial Rx 
c             +/- Generation & Consumption by Wall Chemical Dissolution
c             +   Generation by Radiolysis

c     with Initial Conditions at t = TInit (in subroutine Adjust)
c       Concentrations at all mesh points are set to Bulk Concentrations

c     and Boundaray Conditions at x = 0.0 (crevice mouth)
c       Fixed Concentration at Bulk Concentrations for all species.
c       Fixed Potential at 0.0.

c     and Boundary Conditions at x = L (crevice tip)
c       Zero flux

c     but the following appear in the source term:
c       Flux due to Electrochemical Rx
c       Flux due to Chemical Dissolution (or Adsorption)

c     Method of Lines is used to convert the original form of pde
c     dC/dt = d2C/dx2 to an ode dC/dt = f(C) by differencing d2C/dx2.
c     This increases the number of effective species to NSp*NMesh
c     since concentrations at each mesh point are linked to the
c     concentrations in the rest of mesh points.
c***********************************************************************

      include 'seac1.1.fi'
      integer*4 i,j,m
      real*8 dCdt(NSp,NMesh),C(1) ! C ~ Co
      real*8 dx,Con
      
      if (Debug(9)) write(IOF,*) 'Entering DifEq AT t =',t      

      Ti=t
      call C2Conc            !convert C to Conc
      call CalcFluxData      !calc mesh-dep, but pot-indep parameters
      call CalcPotential     !calc potential using Newton's method
      call CalcFlux(0)       !calc flux

c     Evaluate dc/dt.
c     Outer loop iterates through species, and the inner
c     loop iterates over mesh points.

      do m=1,NMesh

        dx=(x(m+1)-x(m-1))/2.d0
        do i=1,NSp

c         Calc DIVERGENCE upto the TIP ********************************

          dJdx=(FluxH(i,m)-FluxH(i,m-1))/dx
          if ((WidthMouth.ne.WidthTip).and.(m.gt.NMeshEx))
     +      dJdx=(FluxH(i,m)*WidthH(m)-FluxH(i,m-1)*WidthH(m-1))
     +           /dx/Width(m)
          dCdt(i,m)=-dJdx*1.d3 !cc to liter

c         Calc External Flux ******************************************

          if (m.lt.NMeshEx) dCdt(i,m)=dCdt(i,m)+FluxE(i,m)/DifLayer*1.d3
          if (m.eq.NMeshEx) then
            dCdt(i,m)=(FluxH(i,m-1)*DifLayer*WallSide
     +                -FluxH(i,m)*WidthH(m)
     +                +FluxE(i,m)*(WidthMouth-x(m-1)/2.d0*WallSide))
     +                /VolMouth
            if (Gap.gt.0.d0)
     +        dCdt(i,m)=(FluxH(i,m-1)*DifLayer*WallSide
     +                  -FluxH(i,m)*WidthH(m)*Gap/WidthMouth
     +                  +FluxE(i,m)*(Gap-x(m-1)/2.d0*WallSide))
     +                  /VolMouth
          endif

c         Calc CHEMICAL REACTIONS *************************************

          do j=1,NRx
            if (Koef(j,i).ne.0) 
     +        dCdt(i,m)=dCdt(i,m)+RateConst(j)*Koef(j,i)*
     +                  Conc(IR(j,1),m)*Conc(IR(j,2),m)*Conc(IR(j,3),m)
          enddo

c         Calc PRECIPITATION REACTIONS ********************************

          do j=NPRx1,NPRx2
            Con=Conc(IR(j,1),m)*Conc(IR(j,2),m)*Conc(IR(j,3),m)
            if ((Con.gt.Ksp(j)).and.(Koef(j,i).ne.0))
     +        dCdt(i,m)=dCdt(i,m)+RateConst(j)*Koef(j,i)*(Con-Ksp(j))
          enddo

c         Calc WALL & TIP ELECTROCHEMICAL REACTIONS *******************

          call CalcEChem(i,m,dCdt(i,m))

c         Calc WALL & TIP CHEMICAL DISSOLUTION ************************

          if (m.lt.NMeshEx) then      !external surface
            do j=NSDis1,NSDis2
              if (Koef(j,i).ne.0)
     +          dCdt(i,m)=dCdt(i,m)+JDis(j)*Koef(j,i)/
     +                    Width(m)*1.d3*  !cc to liter
     +                  Conc(IR(j,1),m)*Conc(IR(j,2),m)*Conc(IR(j,3),m)
            enddo
          elseif (m.eq.NMeshEx) then  !crevice mouth
            do j=NSDis1,NSDis2
              if (Koef(j,i).ne.0)
     +          dCdt(i,m)=dCdt(i,m)-JDis(j)*Koef(j,i)*
     +                    x(m-1)/2.d0*WallSide/VolMouth*
     +                  Conc(IR(j,1),m)*Conc(IR(j,2),m)*Conc(IR(j,3),m)
            enddo
            do j=NWDis1,NWDis2
              if (Koef(j,i).ne.0)
     +          dCdt(i,m)=dCdt(i,m)+JDis(j)*Koef(j,i)*
     +                    x(m+1)/2.d0*WallSide/VolMouth*
     +                  Conc(IR(j,1),m)*Conc(IR(j,2),m)*Conc(IR(j,3),m)
            enddo
          else                        !inside crevice
            Wid=Width(m)
            if (Gap.gt.0.d0) Wid=Gap  !circular instead rectangular
            do j=NWDis1,NWDis2
              if (Koef(j,i).ne.0)
     +          dCdt(i,m)=dCdt(i,m)+JDis(j)*Koef(j,i)/
     +                    Wid*WallSide*1.d3*  !cc to liter
     +                  Conc(IR(j,1),m)*Conc(IR(j,2),m)*Conc(IR(j,3),m)
            enddo
            if (m.eq.NMesh) then  !calc tip rx
              do j=NTDis1,NTDis2
                if (Koef(j,i).ne.0)
     +            dCdt(i,m)=dCdt(i,m)+JDis(j)*Koef(j,i)/dx*1.d3*  !cc to liter
     +                  Conc(IR(j,1),m)*Conc(IR(j,2),m)*Conc(IR(j,3),m)
              enddo  !through all tip dissolution rx
            endif
          endif

c         Calc RADIOLYSIS *********************************************

          if ((GammaAvg.gt.0d0).or.(NeutAvg.gt.0.d0))
     +      dCdt(i,m)=dCdt(i,m)+GGamma(i)*Gamma(m)+GNeut(i)*Neutron(m)

          if (Debug(10)) write(IOF,10) i,m,dCdt(i,m),i,m,Conc(i,m),
     +      i,m,FluxH(i,m)
         
        enddo  !through all species
      enddo  !through all mesh points

10    format (' dCdt('I3','I3') = '1pe14.6'    Conc('I3','I3') ='e14.6
     +        '    FluxH('I3','I3') ='e14.6)

      if (Debug(9)) write(IOF,*) 'Exiting DifEq AT t =',t      

      return
      end  !of DifEq


      subroutine CalcFlux(PrintFlag)

c***********************************************************************
c     Version:        SEAC 1.1                    1/2/95
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by DifEq, CalcCrevice

c     Calculates flux using the Infinitely Dilute Model.
c     Values at midpoint between two mesh points end with an 'H'.
c     Ex. FluxH(m) = Flux(m+1/2)

c     All length and area units are in cm and cm2.  Since concentrations
c     are in moles/liter, watch out for conversions between
c     liter and 1000 cm3.

c     Flux at crevice tip is NOT part of the flux equation, but appears
c     in the source term along with wall reactions.
c***********************************************************************

      include 'seac1.1.fi'
      real*8 dx,CH,dCdx,dPdx
      integer*4 PrintFlag

      do m=0,NMesh-1

        dx=x(m+1)-x(m)
        dPdx=(Pot(m)-Pot(m+1))/dx  !remember Pot is ECP which is -electrostatic
        ITotalH(m)=0.d0
        do i=1,NSp
          CH=(Conc(i,m)+Conc(i,m+1))/2.d0
          dCdx=(Conc(i,m+1)-Conc(i,m))/dx
          FluxH(i,m)=(-Charge(i)*Mobil(i)*Faraday*CH*dPdx
     +                -Diff(i)*dCdx+CH*VelH(m))*1.d-3
          if ((m.eq.0).and.(EndFlux.eq.0.d0)) FluxH(i,m)=0.d0 
          if (Debug(1).and.(PrintFlag.eq.1)) then
            print 10,m,SpName(i),
     +        -Diff(i)*dCdx*1.d-3,
     +        -Charge(i)*Mobil(i)*Faraday*CH*dPdx*1.d-3,CH*VelH(m)*1.d-3
            write(IOF,10) m,SpName(i),
     +        -Diff(i)*dCdx*1.d-3,
     +        -Charge(i)*Mobil(i)*Faraday*CH*dPdx*1.d-3,CH*VelH(m)*1.d-3
          endif
          ITotalH(m)=ITotalH(m)+FluxH(i,m)*Charge(i)*Faraday
        enddo  !through NSp

      enddo  !through NMesh-1

10    format('Flux D P V ',i4,x,a8,3(1pe14.6))

      return
      end  !of CalcFlux


      subroutine CalcFluxData

c***********************************************************************
c     Version:        SEAC 1.1                    3/13/94
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by DifEq
c***********************************************************************

      include 'seac1.1.fi'
      real*8 dx,CH,dCdx,ZDdCdx

c  calc all mesh-dependent, but potential-independent parameters

      do m=0,NMesh-1

        dx=x(m+1)-x(m)
        Cond(m)=0.d0         !calc conductivity
        CondH=0.d0
        ZDdCdx=0.d0          !calc sum of ZDdCdx
        ZC(m)=0.d0           !calc sum of charge

        do i=1,NSp

          CH=(Conc(i,m)+Conc(i,m+1))/2.d0
          Cond(m)=Cond(m)+Charge(i)*Charge(i)*Mobil(i)*Conc(i,m)
          CondH=CondH+Charge(i)*Charge(i)*Mobil(i)*CH
          dCdx=(Conc(i,m+1)-Conc(i,m))/dx
          ZDdCdx=ZDdCdx+Charge(i)*Diff(i)*dCdx
          ZC(m)=ZC(m)+Charge(i)*Conc(i,m)

          if (Debug(14)) write(IOF,*) m,i,Conc(i,m+1),Conc(i,m),
     +                   Conc(i,m+1)-Conc(i,m),dx,dCdx,
     +                   Charge(i)*Diff(i)*dCdx,ZDdCdx

        enddo  !through NSp

        Cond(m)=Cond(m)*Faraday*Faraday*1.d-3  !liter to cm3
        CondH=CondH*Faraday*Faraday*1.d-3      !liter to cm3
        Wid=WidthH(m)
        Dif=DifLayer
        if (Gap.gt.0.d0) Wid=WidthH(m)*Gap
        if (Gap.gt.0.d0) Dif=WidthMouth*DifLayer
        CondxH(m)=CondH/dx*Wid
        if (m.lt.NMeshEx) CondxH(m)=CondH/dx*Dif
        FZDdCdx(m)=Faraday*ZDdCdx*Wid*1.d-3  !diffusion current
        if (m.lt.NMeshEx) FZDdCdx(m)=Faraday*ZDdCdx*Dif*1.d-3

      enddo  !through NMesh-1

      if (EndFlux.eq.0.d0) then
        CondxH(0)=0.d0
        FZDdCdx(0)=0.d0
      endif

      return
      end  !of CalcFluxData


      subroutine CalcExtFlux

c***********************************************************************
c     Version:        SEAC 1.1                    1/27/94
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by UsrFun
c***********************************************************************

      include 'seac1.1.fi'
      real*8 dx,dy,a,b
      real*8 C,dCdx,dPdx,CondH(IMP)

      dx=DifLayer
      if (PotFinal.ne.1.d-99) PotMet=(PotFinal-PotInit)
     +                              /(TFinal-TInit)*(Ti-TInit)+PotInit
      CondB=Cond(0)/(RefDist-DifLayer)
      do m=1,NMeshEx
        PotE(m)=PotMet
        if (RefDist.gt.0.d0) then
          CondH(m)=0.d0
          ZDdCdx=0.d0
          do i=1,NSp
            C=(Conc(i,m)+ConcBulk(i))/2.d0
            dCdx=(Conc(i,m)-ConcBulk(i))/dx
            CondH(m)=CondH(m)+Charge(i)*Charge(i)*Mobil(i)*C
            ZDdCdx=ZDdCdx+Charge(i)*Diff(i)*dCdx
          enddo  !through NSp
          CondH(m)=CondH(m)*Faraday*Faraday/dx*1.d-3
          DifCur=-Faraday*ZDdCdx*1.d-3
          PotE(m)=(DifCur+CondH(m)*Pot(m)+CondB*PotMet)/(CondH(m)+CondB)
        endif
      enddo

      do m=1,NMeshEx

        dy=(x(m+1)-x(m-1))/2.d0
        if (m.eq.NMeshEx) then
          dy=WidthMouth-x(m-1)/2.d0*WallSide
          if (Gap.gt.0.d0) dy=Gap-x(m-1)/2.d0*WallSide
        endif
        if (Gap.gt.0.d0) dy=dy*WidthMouth

        dPdx=(PotE(m)-Pot(m))/dx
1       IExtNet(m)=0.d0
        dIEdP(m)=0.d0
        a=0.d0
        b=0.d0

        do i=1,NSp
          C=(Conc(i,m)+ConcBulk(i))/2.d0
          dCdx=(Conc(i,m)-ConcBulk(i))/dx
          FluxE(i,m)=(-Charge(i)*Mobil(i)*Faraday*C*dPdx
     +                -Diff(i)*dCdx)*1.d-3
          IExtNet(m)=IExtNet(m)+FluxE(i,m)*Charge(i)*Faraday
          dIEdP(m)=dIEdP(m)+Charge(i)*Charge(i)*Mobil(i)
     +            *Faraday*Faraday*C/dx*1.d-3
          if (RefDist.gt.0.d0) dIEdP(m)=dIEdP(m)-CondH(m)
     +      /(CondH(m)+CondB)*Charge(i)*Charge(i)*Mobil(i)
     +      *Faraday*Faraday*C/dx*1.d-3
          a=a-Charge(i)*Charge(i)*Mobil(i)*Faraday*Faraday*C
          b=b+Diff(i)*dCdx*Charge(i)*Faraday
        enddo  !through NSp

c       fix dPdx so IExtNet=0

        if ((PotMet.eq.1.d-99).and.(dPdx.ne.b/a)) then
          dPdx=b/a
          goto 1
        endif

        IExtNet(m)=IExtNet(m)*dy
        dIEdP(m)=dIEdP(m)*dy
        if (PotMet.eq.1.d-99) dIEdP(m)=0.d0
        PotE(m)=dPdx*dx+Pot(m)

      enddo  !through NMeshEx

      if ((RefDist.gt.0.d0).and.(PotMet.ne.1.d-99)) Pot(0)=PotE(1)

      return
      end  !of CalcExtFlux


      subroutine CalcEChem(i,m,dCdt)

c***********************************************************************
c     Version:        SEAC 1.1                    1/3/95
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by DifEq, UsrFun
c     Limiting currents are calculated only for external surface.
c***********************************************************************

      include 'seac1.1.fi'
      integer*4 i,m
      real*8 dCdt,dx,a,d,Wid,area

      dx=(x(1)-x(0))/2.d0     !for m=0
      if (m.gt.0) dx=(x(m+1)-x(m-1))/2.d0
      area=dx
      Wid=Width(m)
      if (Gap.gt.0.d0) then  !circular instead of rectangular for Alkire
        area=dx*Wid
        Wid=Gap
        if (m.lt.NMeshEx) Wid=DifLayer
      endif

      IWallNet(m)=0.d0
      dIWdP(m)=0.d0
      IAnod(m)=0.d0

      if (m.lt.NMeshEx) then  !external surface

        call CalcEChemRx(i,m,NSRx1,NSRx2,dCdt,area,Wid*1.d-3,1)

      elseif (m.eq.NMeshEx) then  !crack mouth

        dx=-x(m-1)/2.d0
        area=dx
        if (Gap.gt.0.d0) area=dx*WidthMouth
        call CalcEChemRx(i,m,NSRx1,NSRx2,dCdt,area,
     +    VolMouth/dx/WallSide,1)
        ISurfMouth=IWallNet(m)
        dISMdP=dIWdP(m)
        a=IAnod(m)
        d=dx

        dx=x(m+1)/2.d0*WallSide
        area=dx
        if (Gap.gt.0.d0) area=dx*WidthMouth
        call CalcEChemRx(i,m,NWRx1,NWRx2,dCdt,area,VolMouth/dx,1)
        IAnod(m)=(a*d+(IAnod(m)-a)*dx)/(d+dx)

      else  !inside crack

        call CalcEChemRx(i,m,NWRx1,NWRx2,dCdt,
     +    area*WallSide,Wid/WallSide*1.d-3,1) !up/low wall
        if (m.eq.NMesh) then  !calc tip rx
          ITip=IAnod(m)       !for crack growth rate
          ITotalH(m)=IWallNet(m)  !current density at tip
          area=WidthTip
          if (Gap.gt.0.d0) area=WidthTip*Gap
          call CalcEChemRx(i,m,NTRx1,NTRx2,dCdt,area,dx*1.d-3,1)
          ITip=IAnod(m)-ITip
          IAnod(m)=ITip       !at m=NMesh, IAnod is only for tip
          ITotalH(m)=(ITotalH(m)-IWallNet(m))/area    !I & J opposite direction
        endif

      endif

      return
      end  !of CalcEChem


      subroutine CalcEChemRx(i,m,j1,j2,dCdt,area,dx,lim)

c***********************************************************************
c     Version:        SEAC 1.1                    2/2/95
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by CalcEChem

c     Calculates wall/tip electrochemical reactions for a given species
c     Includes Limiting Current, Passivation, Trans-Passivation
c     lim = 0 -> don't calc lim current;  lim = 1 -> calc lim current
c***********************************************************************

      include 'seac1.1.fi'
      integer*4 i,j,j1,j2,m,lim
      real*8 dCdt,area,dx
      real*8 IWall,ILim,JWall,JWallTot,PotPass,PotPitt
      real*8 W,dWdP,Anod,wa,dW

      W=0.d0
      dWdP=0.d0
      Anod=0.d0
      JWallTot=0.d0

      do j=j1,j2

        IWall=IRate(j)*exp(EACoef(j)*Pot(m))
     +       *Conc(IR(j,1),m)*Conc(IR(j,2),m)*Conc(IR(j,3),m)

c       Limiting current

        wa=IWall
        ILim=ILimit(j*lim)*exp(abs(EACoef(j))/EACoef(j)*Pot(m))
        IWall=1.d0/(1.d0/IWall+1.d0/ILim)
        dW=IWall*IWall*EACoef(j)*(1.d0/wa+1.d0/abs(EACoef(j))/ILim)

c       Passivity and Pitting

        PotPass=PotPas(j)+PotPasCoef(j)/FRT*log(Conc(iHydro,m))
        PotPitt=PotPit(j)-PotPitCoef(j)/FRT*log(Conc(iChlor,m))
        if (Pot(m).gt.PotPass) then
          wa=IWall
          if (IPas(j).lt.IWall) then
            IWall=IPas(j)
            dW=0.d0
          endif
          if (Pot(m).gt.PotPitt) then
            IWall=IWall*exp(EACoef(j)*(Pot(m)-PotPitt))
            dW=IWall*EACoef(j)
          endif
        endif

        if (wa.eq.0.d0) IWall=0.d0
        W=W+IWall
        if (wa.eq.0.d0) dW=0.d0
        dWdP=dWdP+dW
        if (NetCharge(j).gt.0) Anod=Anod+IWall

        if (Koef(j,i).ne.0) then
          JWall=Koef(j,i)*IWall/NetCharge(j)/Faraday
          JWallTot=JWallTot+JWall
        endif

      enddo

      IWallNet(m)=IWallNet(m)+W*area
      dIWdP(m)=dIWdP(m)+dWdP*area
      IAnod(m)=IAnod(m)+Anod
      if (i.gt.0) dCdt=dCdt+JWallTot/dx  !already in mol/L-sec

      return
      end  !of CalcEChemRx


      subroutine CalcPotential

c***********************************************************************
c     Version:        SEAC 1.1                    3/14/94
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by DifEq
c     Calculates potential using Newton's method based on nonlinear
c     circuit elements.
c***********************************************************************

      include 'seac1.1.fi'

      if (PotMet.eq.1.d-99) then  !metal at free corrosion potential
        call MNewt(NewtMax,Pot(0),NMesh+1,NewtXTol,NewtFTol)
      else                        !metal polarized to PotMet
        call MNewt(NewtMax,Pot(1),NMesh,NewtXTol,NewtFTol)
      endif

      return
      end  !of CalcPotential


      real*8 function PotEq(j,m)

c***********************************************************************
c     Version:        SEAC 1.1                    3/8/94
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by DifEq, CalcFlux
c     Calculates equilibrium potential using the Nernst equation.
c     For dissolved gas, the partial pressure is obtained from
c     concentration using Henry's law:

c       Part Pres [atm] = Henry's Const [atm/mol/L] * Conc [mol/L]

c     Ref: Himmelblau, DM;"Solubility of Inert Gases in Water", J Chemical
c     and Eng Data, 5 (1960) 10.
c     950102 not used in current implementation. jhc
c***********************************************************************

      include 'seac1.1.fi'
      real*8 Ratio
      integer*4 j,k,m

      if (PotStd(j).eq.1.d-99) then
        PotEq=PotMet
        return
      endif

      Ratio=1.d0
      do k=1,INP
        Ratio=Ratio*Conc(IP(j,k),m)
        if (Henry(IP(j,k)).ne.0.d0) Ratio=Ratio*Henry(IP(j,k))
      enddo
      do k=1,INR
        Ratio=Ratio/Conc(IR(j,k),m)
        if (Henry(IR(j,k)).ne.0.d0) Ratio=Ratio/Henry(IR(j,k))
      enddo
      PotEq=PotStd(j)-log(Ratio)/NetCharge(j)/FRT

      return
      end  !of PotEq


      subroutine MNewt(ntrial,xx,n,tolx,tolf)

c***********************************************************************
c     Version:        SEAC 1.1                    1/28/94
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by CalcFlux
c     Calls UsrFun, LUDcmp, LUBksb
c
c     Solve for potential in the crevice using Newton's method.
c     Press, William H., Brian P. Flannery, Saul A. Teukolsky,
c          William T. Vetterling; "Numerical Recipes - Fortran version",
c          1989, Cambridge University Press, p 269-272.
c***********************************************************************

      include 'seac1.1.fi'
      INTEGER*4 n,ntrial
      REAL*8 tolf,tolx,xx(n)
      INTEGER*4 i,k
      REAL*8 errf,errx,fvec(IMP+1),a(IMP+1),b(IMP+1),c(IMP+1),p(IMP+1)

      do k=1,ntrial
        call UsrFun(fvec,a,b,c)
        errf=0.d0
        do i=1,n
          errf=errf+abs(fvec(i))
        enddo
        if(errf.le.tolf) return
        do i=1,n
          p(i)=-fvec(i)
        enddo
        call TriDag(a,b,c,p,n)
        errx=0.d0
        do i=1,n
          errx=errx+abs(p(i))
          xx(i)=xx(i)+p(i)
        enddo
        if(errx.le.tolx) return
      enddo
      
      print *,'#### Potential did not converge after',ntrial,' trials!'
      print *,'Pot(NMesh) =',xx(n),'  at t =',Ti
      print *,'errf =',errf,' /',tolf
      print *,'errx =',errx,' /',tolx
      pause

      end  !of MNewt


      subroutine UsrFun(f,a,b,c)

c***********************************************************************
c     Version:        SEAC 1.1                    1/28/94
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by MNewt
c
c     Sets up f and its Jacobian to be solve for potential in the
c     crevice using Newton's method
c***********************************************************************

      include 'seac1.1.fi'
      integer*4 m,m1,n
      real*8 f(*),a(*),b(*),c(*),dummy

      m1=1                        !starting node for potential calculation
      if (PotMet.eq.1.d-99) m1=0  !metal at free corrosion potential

c  calc CURRENT from WALL & TIP ELECTROCHEMICAL REACTIONS *************

      call CalcExtFlux            !external surface flux
      do m=m1,NMesh
        call CalcEChem(0,m,dummy)
      enddo  !through all mesh points

c  calc function vector & Jacobian ************************************

      do m=m1,NMesh
        n=m+1-m1
        if (m.eq.0) then
          f(n)=FZDdCdx(m)+CondxH(m)*(Pot(m)-Pot(m+1))+IWallNet(m)
          b(n)=CondxH(m)+dIWdP(m)
          c(n)=-CondxH(m)
        elseif (m.lt.NMeshEx) then
          f(n)=-FZDdCdx(m-1)-CondxH(m-1)*(Pot(m-1)-Pot(m))
     +         +FZDdCdx(m)+CondxH(m)*(Pot(m)-Pot(m+1))
     +         +IWallNet(m)+IExtNet(m)
          b(n)=CondxH(m-1)+CondxH(m)+dIWdP(m)+dIEdP(m)
          c(n)=-CondxH(m)
        elseif (m.eq.NMeshEx) then
          f(n)=(-FZDdCdx(m-1)-CondxH(m-1)*(Pot(m-1)-Pot(m))+ISurfMouth)
     +          *WallSide+FZDdCdx(m)+CondxH(m)*(Pot(m)-Pot(m+1))
     +          +IWallNet(m)-ISurfMouth+IExtNet(m)
          b(n)=(CondxH(m-1)+dISMdP)*WallSide+CondxH(m)+dIWdP(m)-dISMdP
     +         +dIEdP(m)
          c(n)=-CondxH(m)
        elseif (m.lt.NMesh) then
          f(n)=-FZDdCdx(m-1)-CondxH(m-1)*(Pot(m-1)-Pot(m))
     +         +FZDdCdx(m)+CondxH(m)*(Pot(m)-Pot(m+1))
     +         +IWallNet(m)+IExtNet(m)
          b(n)=CondxH(m-1)+CondxH(m)+dIWdP(m)+dIEdP(m)
          c(n)=-CondxH(m)
        else
          f(n)=-FZDdCdx(m-1)-CondxH(m-1)*(Pot(m-1)-Pot(m))
     +         +IWallNet(m)+IExtNet(m)
          b(n)=CondxH(m-1)+dIWdP(m)+dIEdP(m)
        endif
        a(n)=-CondxH(m-1)
        if (m.eq.NMeshEx) a(n)=-CondxH(m-1)*WallSide

      enddo
      
      return
      end  !of UsrFun


      subroutine TriDag(a,b,c,u,n)

c***********************************************************************
c     Version:        SEAC 1.1                    1/29/94
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by MNewt
c
c     Solves tridiagonal matrix.
c     Press, William H., Brian P. Flannery, Saul A. Teukolsky,
c          William T. Vetterling; "Numerical Recipes - Fortran version",
c          1989, Cambridge University Press, p 40-41.
c***********************************************************************

      include 'seac1.1.fi'
      integer*4 n
      real*8 a(n),b(n),c(n),u(n)
      integer*4 j
      real*8 bet,gam(IMP)

      if(b(1).eq.0.d0)pause 'tridag: rewrite equations'
      bet=b(1)
      u(1)=u(1)/bet
      do 11 j=2,n
        gam(j)=c(j-1)/bet
        bet=b(j)-a(j)*gam(j)
        if(bet.eq.0.d0)pause 'tridag failed'
        u(j)=(u(j)-a(j)*u(j-1))/bet
11    continue
      do 12 j=n-1,1,-1
        u(j)=u(j)-gam(j+1)*u(j+1)
12    continue

      return
      end  !of TriDag


      subroutine X2MeshPoints

c***********************************************************************
c     Version:        SEAC 1.1                    11/5/93
c     Author :        JOHN H. CHUN
c***********************************************************************
c     Called by Initialize
c
c     Nodalize x-axis with mesh points
c***********************************************************************

      include 'seac1.1.fi'
      real*8 dx

      dx=0.d0

      if (NMeshEx.gt.0) then
        dx=Height/NMeshEx
        x(0)=-Height
      endif

      do m=1,NMeshEx
        x(m)=(m-NMeshEx)*dx
      enddo

      if (NMeshIn.gt.0) then
        dx=Length/NMeshIn
        x(NMesh)=Length
      endif

      if (MeshSpace.eq.0) then      !uniform mesh spacing
        do m=NMeshEx,NMesh-1
          x(m)=(m-NMeshEx)*dx
        enddo
      elseif (MeshSpace.eq.1) then  !non-uniform mesh spacing
        mm=0
        do m=1,NMesh
          mm=mm+m
        enddo
        x1=Length/mm
        x(0)=0.d0
        do m=1,NMesh-1
          x(m)=x(m-1)+x1
        enddo
      elseif (MeshSpace.eq.2) then
        mm=0
        do m=1,NMesh
          mm=mm+m
        enddo
        x1=Length/mm
        x(0)=0.d0
        do m=1,NMesh-1
          x(m)=x(m-1)+x1*m
        enddo
      else
        a=1.3d0              !a and b are shape factors
        b=5.d0
        do m=1,NMesh/2
          y1=datand(((m-1)*dx-Length/2.d0)*b/Length)*Pi/180.d0
          x(m)=a*Length/Pi*y1+Length/2.d0
          y2=datand(((NMesh-m+1)*dx-Length/2.d0)*b/Length)*Pi/180.d0
          x(NMesh-m)=a*Length/Pi*y2+Length/2.d0
        enddo
        if ((NMesh/2)*2.eq.NMesh) x(NMesh/2)=Length/2.d0 !corect for even NMesh
      endif
      
      x(NMesh+1)=x(NMesh)    !used for TIP flux

c  straight or tapered width

      if (Width(1)+WidthTip.le.0.d0) then
        do m=NMeshEx,NMesh
          Width(m)=WidthMouth
          WidthH(m)=WidthMouth
        enddo
      elseif ((Width(1).le.0.d0).and.(WidthTip.gt.0.d0)) then
        slope=(WidthMouth-WidthTip)/Length
        do m=NMeshEx,NMesh
          Width(m)=WidthMouth-slope*x(m)
          WidthH(m)=WidthMouth-slope*(x(m)+x(m+1))/2.d0
        enddo
      endif
      WidthTip=Width(NMesh)
      Width(NMesh+1)=WidthTip

      do m=0,NMeshEx-1
        Width(m)=DifLayer
        if (Gap.gt.0.d0) Width(m)=WidthMouth
      enddo

      if (NMeshIn.eq.0) then
        Width(NMesh)=0.d0
        WidthMouth=0.d0
        WidthTip=0.d0
      endif

c  volume at mouth

      m=NMeshEx
      Wid=WidthMouth
      if (Gap.gt.0.d0) Wid=Gap
      if (m.eq.0) then
        VolMouth=x(m+1)/2.d0*Wid*1.d-3
      else
        VolMouth=((Wid-x(m-1)/2.d0*WallSide)*DifLayer
     +            +x(m+1)/2.d0*Wid)*1.d-3
      endif

c  print for debug

      if (Debug(13)) then
        write(IOF,*)
        write(IOF,*) 'Mesh point spacing:'
        write(IOF,*)
        write(IOF,*) 'y1 =',y1,'  y2 =',y2
        do m=0,NMesh
          write(IOF,130) m,x(m),Width(m)
        enddo
        if (NMeshEx.gt.0) write(IOF,*) 'Mouth Volume ',VolMouth
        write(IOF,*)
      endif
130   format ('x('i3') = '1pe15.7,'  Width = ',e15.7)

      return
      end  !of X2MeshPoints
