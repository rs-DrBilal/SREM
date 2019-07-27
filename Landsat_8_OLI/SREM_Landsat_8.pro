;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

;   Code:
;     SREM_Landsat_OLI.pro
;   Purpose:
;     Atmospheric Correction/
;     Surface Reflectance Estimation
;     from Landsat 8 OLI Satellite
;
;   Input Data:
;     a. TOA reflectance
;     b. Sun and Satellite (view) Zenith Angles
;     C. Sun and Satellite (view) Azimuth Angles
;
;   Developed By:
;     Dr. Muhammad Bilal
;     Distinguished Professor at NUIST (Nanjing, China)
;     Email: muhammad.bilal@connect.polyu.hk
;
;   Reference:
;     https://doi.org/10.3390/rs11111344
;
;   Citation:
;     Bilal, M.; Nazeer, M.; Nichol, J.E.; Bleiweiss, M.P.; Qiu, Z.;
;     Jäkel, E.; Campbell, J.R.; Atique, L.;Huang, X.; Lolli, S.
;     A Simplified and Robust Surface
;     Reflectance Estimation Method (SREM) for Use
;     over Diverse Land Surfaces Using Multi-Sensor Data
;     Remote Sens. 2019, 11, 1344.
;
;   NOTE:
;     This code can be modified for DN/
;     TOA radiance. For modification, please
;     contact at my email (above mentioned)
;
;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

Pro SREM_Landsat_8

  COMPILE_OPT idl2
  ENVI, /RESTORE_BASE_SAVE_FILES
  ENVI_BATCH_INIT, /NO_STATUS_WINDOW
  print, '     ' + 'Batch Processing Start'
  print, '     ' + systime()

  ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

  ;  Define Input and Output Directories
  InputTOA=Dialog_Pickfile(Title='Select Landsat 8 Band3 TOA Reflectance Images',/DIRECTORY)
  OutputSR=Dialog_Pickfile(Title='Select SREM SR OutPut Folder',/DIRECTORY)

  ;  Search for inpiut data
  ;  Please change '.tif' format according to your data
  FTOA = File_Search(InputTOA + '*.tif',COUNT=count_in)
  nfiles = n_elements(FTOA)

  ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$ $$$$$$$$$$$$$$$$$
  ;
  ;  For Loop for batch processing
  ;  FOR k=0, 0 DO BEGIN
  FOR k=0, nfiles-1 DO BEGIN
    string_count = strtrim(nfiles,2)

    ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ;  Open Rtoa = TOA reflectance
    ENVI_OPEN_FILE, FTOA[k], r_FiD=FiD_toa
    IF (FiD_toa EQ -1) then return
    ENVI_FILE_QUERY, FiD_toa, dims=dims, nb=nb, ns=toa_ncols, nl=toa_nrows
    pos = [0]
    map_info = envi_get_map_info(FiD=FiD_toa)
    Rtoa = Double(ENVI_GET_DATA(FiD=FiD_toa, dims=dims, pos=[0]))

    ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    ;  Define szen = sun zenith angle
    szen =
    ;  Define sazm = sun azimuth angle
    sazm =
    ;  Define vzen = satellite zenith angle
    vzen =
    ;  Define vazm = satellite azimuth angle
    vazm =

    ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    ;  !DTOR = Degree to Radian
    ;  cos(vzen) = cvz
    ;  cos(szen) = csz
    ;  sin(vzen) = svz
    ;  sin(szen) = ssz
    CVZ = cos(vzen*!DTOR)
    CSZ = cos(szen*!DTOR)
    SVZ = sin(vzen*!DTOR)
    SSZ = sin(szen*!DTOR)

    ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    ;
    ;  calculate relative azimuth angle (RelAzm)
    RelAA = ABS((sazm) - (vazm))
    index = where(RelAA gt 180.d, count)
    if count gt 0 then begin
      RelAA[index] = 360.- RelAA[index]
    endif
    index = where(RelAA le 180.d, count)
    if count gt 0 then begin
      RelAA[index] = 180.- RelAA[index]
    endif

    ;  calculate cosinse of RelAA
    ;  cos (RelAzm) = CRA
    CRA = cos(RelAA*!DTOR)

    ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    ;  calculate scattering angle (SCA)
    ;  Reference: Eq. 6 in Levy et al. (2007)
    ;  DOI: https://doi.org/10.1029/2006JD007811
    SCA =  acos((-CSZ * CVZ) + ((SSZ * SVZ) * CRA))

    ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    ;  calculate Rayleigh Optical Depth (ROD)
    ;  CW = Band central wavelenght (CW) in micro meter
    ;       and shoule be changed accoring to input band
    ;  Reference: Eq. 9 in SREM paper
    CW = 0.561d ; band 3 (green)
    ROD = (0.008569D * ((CW)^(-4)))*(1.D + (0.0113D *(CW)^(-2)) + (0.00013D * (CW)^(-4)))

    ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    ;  calculate Rayleigh Phase Function (Pr)
    A = 0.9587256d
    B = 1.d - A
    ;  Reference: Eq. 10 in SREM paper
    Pr =  ((3.D * A)/(4.D + B)) * (1 + (cos(SCA) * cos(SCA)))

    ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    ;  calulate Rayleigh Reflectance (Rr)
    ;  Reference: Eq. 8 in SREM paper
    M = (1/CSZ) + (1/CVZ)
    ;  Reference: Eq. 7 in SREM paper
    Rr = Pr * (1 - exp(-M * ROD))/(4.D * (CSZ + CVZ))

    ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    ;  calcualte Total transmission (T)
    ;  Transmission on sun-surface path
    ;  Reference: Eq.12 in SREM paper
    Es = exp(-(ROD)/CSZ) * (exp((0.52D * ROD)/CSZ)-1)
    Ts = exp(-(ROD)/CSZ) + Es

    ;  Transmission on surface-satellite path
    ;  Reference: Eq. 13 in SREM paper
    Ev = exp(-(ROD)/CVZ) * (exp((0.52D * ROD)/CVZ)-1)
    Tv = exp(-(ROD)/CVZ) + Ev
    ;  Total transmission
    T = Ts * Tv

    ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    ;  calculate atmospheric backscattering ratio
    ;  Reference: Eq. 11 in SREM paper
    Satm = (0.92D * ROD) * exp(-ROD)

    ;  calculate SREM SR
    ;  Reference: Eq. 5 in SREM paper
    SREM = ((Rtoa - Rr) / ((((Rtoa - Rr) * Satm) + T)))

    ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    ;  Define output FileName
    FileName = STRJOIN(STRSPLIT(STRMID(FILE_BASENAME(FTOA[k], '.tif'), 0, 34), '.', /EXTRACT), '_')

    ;  write output File (SREM)
    out_name_out = OutputSR + FileName +'_SREM'
    map_info = envi_get_map_info(FiD=FiD_toa)
    ENVI_WRITE_ENVI_FILE, SREM, out_name=out_name_out, map_info=map_info

    envi_open_file, out_name_out, r_FiD=FiD_out
    envi_file_query, FiD_out, dims=dims, nb=nb, ns=ns, nl=nl

    out_name=OutputSR + FileName +'_SREM_Landsat_8_OLI_B3.tif'
    ENVI_OUTPUT_TO_EXTERNAL_FORMAT, FiD=FiD_out, pos=[0], dims=dims, out_name=out_name,/TIFF

    file_delete, out_name_out
    file_delete, out_name_out+'.hdr'

    ;$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    Print, FileName + ' ' + strtrim(k+1,2) + ' of ' + string_count + ' Processed'
  Endfor
  Print, '     ' + systime()
  Print, '     ' +  'Batch Processing End'
End