# -*- coding: utf-8 -*-
"""
Created on Mon Sep 17 15:17:58 2018

@author: bec
"""



#def FlipDict(dictionary):
#    dictb = dict((v,k) for k, v in dictionary.items())
#    return dictb

#def calc_basinmean(perc_fh, lu_fh):
#    """
#    Calculate the mean of a map after masking out the areas outside an basin defined by
#    its landusemap.
#
#    Parameters
#    ----------
#    perc_fh : str
#        Filehandle pointing to the map for which the mean needs to be determined.
#    lu_fh : str
#        Filehandle pointing to landusemap.
#
#    Returns
#    -------
#    percentage : float
#        The mean of the map within the border of the lu_fh.
#    """
#    output_folder = tempfile.mkdtemp()
#    perc_fh = MatchProjResNDV(lu_fh, np.array([perc_fh]), output_folder)
#    EWR = OpenAsArray(perc_fh[0])
#    LULC = OpenAsArray(lu_fh)
#    EWR[np.isnan(LULC)] = np.nan
#    percentage = np.nanmean(EWR)
#    shutil.rmtree(output_folder)
#    return percentage
#
#def GapFil(input_tif, footprint, output_folder, method = 'max'):
#    """
#    Gapfil a raster by filling with the minimum, maximum or median of nearby pixels.
#
#    Parameters
#    ----------
#    input_tif : str
#        Raster to be filled.
#    footprint : ndarray
#        Boolean array describing the area in which to look for the filling value.
#    output_folder : str
#        Folder to store gapfilled map.
#    method : str, optional
#        Method to use for gapfilling, options are 'max', 'min' or 'median'. Default is 'max'.
#        
#    Returns
#    -------
#    fh : str
#        Location of the gapfilled map.
#    """
#    driver, NDV, xsize, ysize, GeoT, Projection = GetGeoInfo(input_tif)
#    population = OpenAsArray(input_tif)
#    
#    if method == 'median':
#        population_gaps = ndimage.median_filter(population, footprint = footprint)
#    if method == 'max':
#        population_gaps = ndimage.maximum_filter(population, footprint = footprint)
#    if method == 'min':
#        population_gaps = ndimage.minimum_filter(population, footprint = footprint)
#        
#    population[np.isnan(population)] = population_gaps[np.isnan(population)]
#    
#    fn = os.path.split(input_tif)[1].replace('.tif','_gapfilled.tif')
#    fh = os.path.join(output_folder, fn)
#    
#    CreateGeoTiff(fh, population, driver, NDV, xsize, ysize, GeoT, Projection)
#    
#    return fh
#        
#def Multiply(fh1, fh2, fh3):
#    """
#    Multiply two maps with eachother and store the results in a new map.
#    
#    Parameters
#    ----------
#    fh1 : str
#        Filehandle pointing to map to be multiplied with fh2.
#    fh2 : str
#        Filehandle pointing to map to be multiplied with fh1.
#    fh3 : str
#        Filehandle indicating where to store the results.
#    """
#    FH1 = OpenAsArray(fh1)
#    FH2 = OpenAsArray(fh2)
#    
#    FH3 = FH1 * FH2
#    
#    if not os.path.exists(os.path.split(fh3)[0]):
#        os.makedirs(os.path.split(fh3)[0])
#    
#    driver, NDV, xsize, ysize, GeoT, Projection = GetGeoInfo(fh1)
#    CreateGeoTiff(fh3, FH3, driver, NDV, xsize, ysize, GeoT, Projection)
#        
#def Unzip(list_of_tuples):
#    """
#    Creates seperate lists from values inside tuples in a list.
#    
#    Parameters
#    ----------
#    list_of_tuples : list
#        List containing tuples.
#        
#    Returns
#    -------
#    out : list
#        List with arrays with the values of the tuples.
#        
#    Examples
#    --------
#    >>> Unzip([(2,3,4),(5,6,7),(1,2,3)])
#    [np.array([2, 5, 1]), np.array([3, 6, 2]), np.array([4, 7, 3])]
#    """
#    out = [np.array(list(t)) for t in zip(*list_of_tuples)]
#    return out
#        
#def GetMonthLabels():
#    """
#    Function to create a dictionary with two digit month labels, alternative to 
#    applying zfill(2) to a string.
#    
#    Returns
#    -------
#    month_labels : dict
#        Dictionary with two digit months labels.
#    """
#    month_labels = {1:'01',2:'02',3:'03',4:'04',5:'05',6:'06',7:'07',8:'08',9:'09',10:'10',11:'11',12:'12'}
#    return month_labels
#
#def MergeDictionaries(list_of_dictionaries):
#    """
#    Merge multiple dictionaries into one, gives a warning if keys are 
#    overwritten.
#    
#    Parameters
#    ----------
#    list_of_dictionaries : list
#        List containing the dictionaries to merge.
#        
#    Returns
#    -------
#    merged_dict : dict
#        The combined dictionary.
#    """
#    merged_dict = dict()
#    expected_length = 0
#    for dic in list_of_dictionaries:
#        expected_length += len(dic.keys())
#        merged_dict = dict(merged_dict.items() + dic.items())
#    if expected_length is not len(merged_dict):
#        print "WARNING: It seems some station(s) with similar keys have been overwritten ({0} != {1}), keys: {2}".format(expected_length, len(merged_dict))
#    return merged_dict    
#

#def AssertPresentKeys(dictionary1, dictionary2):
#    """
#    Check if the keys of dictionary 1 are present in dictionary 2.
#    
#    Parameters
#    ----------
#    dictionary1 : dict
#        Dictionary with keys.
#    dictionary2 : dict
#        Another dictionary with keys.
#    """
#    for key in dictionary1.keys():
#        assert key in dictionary2.keys(), "{0} key is missing in dictionary"
#



#def Ysum(fhs, fh3):
#    """
#    sum maps with each other and store the results in a new map.
#   
#    Parameters
#    ----------
#    fhs : list of maps to sum
#    fh3 : str
#        Filehandle indicating where to store the results.
#    """
#    FH3 = OpenAsArray(fhs[0]) * 0
#    for fh in fhs:
#        FH3 = np.nansum((FH3, OpenAsArray(fh)),0)
#   
#    if not os.path.exists(os.path.split(fh3)[0]):
#        os.makedirs(os.path.split(fh3)[0])
#   
#    driver, NDV, xsize, ysize, GeoT, Projection = GetGeoInfo(fhs[0])
#    CreateGeoTiff(fh3, FH3, driver, NDV, xsize, ysize, GeoT, Projection)
#    
#
#def Aggregate(fhs, fh = None):
#    """
#    Calculate the sum of multiple geotiffs.
#    
#    Parameters
#    ----------
#    fhs : list
#        List of filehandles to be aggregated.
#    fh  : str, optional
#        String specifying where to store output, default is None.
#        
#    Returns
#    -------
#    fh : str
#        Filehandle specifying where the aggregated map has been stored.
#    """
#    AssertProjResNDV([fhs])
#    if fh is None:
#        temp_folder = tempfile.mkdtemp()
#        fh = os.path.join(temp_folder,'temp.tif')
#    DATA = OpenAsArray(fhs[0])
#    for i in range(1,len(fhs)):
#        DATA += OpenAsArray(fhs[i])
#    driver, NDV, xsize, ysize, GeoT, Projection = GetGeoInfo(fhs[0])
#    CreateGeoTiff(fh, DATA, driver, NDV, xsize, ysize, GeoT, Projection)
#    return fh    
#
#def ZeroesDictionary(dictionary):
#    """
#    Creates a dictionary with the same keys as the input dictionary, but
#    all values are zero.
#    
#    Parameters
#    ----------
#    dictionary : dict
#        dictionary to be copied.
#    
#    Returns
#    -------
#    null_dictionary : dict
#        dictionary with zero values.
#    """
#    null_dictionary = dict()
#    for key in dictionary.keys():
#        null_dictionary[key] = 0.0
#    return null_dictionary
#            
#def ZonalStats(fhs, dates, output_dir, quantity, unit, location, color = '#6bb8cc'):
#    """
#    Calculate and plot some statictics of a timeseries of maps.
#    
#    Parameters
#    ----------
#    fhs : ndarray
#        Filehandles pointing to maps.
#    dates : ndarray
#        Datetime.date objects corresponding to fhs.
#    output_dir : str
#        Folder to save the graphs.
#    quantity : str
#        Quantity of the maps.
#    unit : str
#        Unit of the maps.
#    location : str
#        Location name of the maps.
#    color : str, optional
#        Color in which the graphs will be plotted, default is '#6bb8cc'.
#        
#    Returns
#    -------
#    monthly_average : float
#        Monthly spatial average.
#    yearly_average : float
#        Yearly spatial average.
#        
#    Examples
#    --------
#    >>> ZonalStats(p_fhs, p_dates, output_dir, 'Precipitation', 'mm/month', 'North-Vietnam')
#    
#    >>> ZonalStats(et_fhs, et_dates, output_dir, 'Evapotranspiration', 'mm/month', 'South-Vietnam')
#    
#    """
#    ts = np.array([])
#    
#    data_monthly_ts = dict()
#    data_monthly_counter = dict()
#    months = np.unique([date.month for date in dates])
#    
#    for month in months:
#        data_monthly_ts[month] = 0
#        data_monthly_counter[month] = 0
#    
#    data_yearly_ts = dict()
#    data_yearly_counter = dict()
#    years = np.unique([date.year for date in dates])
#    
#    for year in years:
#        data_yearly_ts[year] = 0
#        data_yearly_counter[year] = 0
#    
#    for date in dates:
#        
#        DATA = OpenAsArray(fhs[dates == date][0])
#        data = np.nanmean(DATA)
#        ts = np.append(ts, data)
#        data_monthly_ts[date.month] += data
#        data_monthly_counter[date.month] += 1
#        data_yearly_ts[date.year] += data
#        data_yearly_counter[date.year] += 1
#    
#    monthly_ts = np.array(data_monthly_ts.values()) / np.array(data_monthly_counter.values())
#    months = np.array(data_monthly_ts.keys())
#    
#    yearly_mask = np.array(data_yearly_counter.values()) == 12
#    yearly_ts = np.array(data_yearly_ts.values())[yearly_mask] / np.array(data_yearly_counter.values())[yearly_mask]
#    years = np.array(data_yearly_ts.keys())[yearly_mask]
#      
#    idx = np.argsort(dates)
#      
#    fig = plt.figure(figsize = (10,5))
#    plt.clf()
#    plt.grid(b=True, which='Major', color='0.65',linestyle='--', zorder = 0)
#    ax = plt.subplot(111)
#    ax.plot(dates[idx], ts[idx], '-k')
#    ax.fill_between(dates[idx], ts[idx], color = color)
#    ax.set_xlabel('Time')
#    ax.set_ylabel(quantity + ' ' + unit)
#    ax.set_title(quantity + ', ' + location)
#    fig.autofmt_xdate()
#    [i.set_zorder(10) for i in ax.spines.itervalues()]
#    plt.savefig(os.path.join(output_dir, quantity + '_' + location + '_ts.png'))
#    plt.close(fig)
#        
#    fig = plt.figure(figsize = (10,5))
#    plt.clf()
#    plt.grid(b=True, which='Major', color='0.65',linestyle='--', zorder = 0)
#    ax = plt.subplot(111)
#    ax.bar(months - 0.4, monthly_ts, 0.8, color = color)
#    ax.set_xlabel('Time [month]')
#    ax.set_xlim([0, max(months)+1])
#    ax.set_xticks(months)
#    ax.set_ylabel(quantity + ' ' + unit)
#    ax.set_title('Monthly average ' + quantity + ', ' + location)
#    [i.set_zorder(10) for i in ax.spines.itervalues()]
#    plt.savefig(os.path.join(output_dir, quantity + '_' + location + '_monthly.png'))
#    plt.close(fig)
#        
#    fig = plt.figure(figsize = (10,5))
#    plt.clf()
#    plt.grid(b=True, which='Major', color='0.65',linestyle='--', zorder = 0)
#    ax = plt.subplot(111)
#    ax.bar(years - 0.4, yearly_ts, 0.8, color = color)
#    ax.set_xlabel('Time [year]')
#    ax.set_xlim([min(years) - 1, max(years)+1])
#    ax.set_ylabel(quantity + ' ' + unit)
#    ax.set_title('Yearly average ' + quantity + ', ' + location)
#    [i.set_zorder(10) for i in ax.spines.itervalues()]
#    plt.savefig(os.path.join(output_dir, quantity + '_' + location + '_yearly.png'))
#    plt.close(fig)
#    
#    monthly_max = np.nanmax(monthly_ts)
#    monthly_average = np.nanmean(monthly_ts)
#    yearly_average = np.nanmean(yearly_ts)
#    
#    return monthly_max, monthly_average, yearly_average