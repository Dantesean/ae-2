def mp_worker(args):
    import os, sys, time, json, urllib3, requests
    (date, longitude, latitude) = args
    base_url = r"https://power.larc.nasa.gov/cgi-bin/v1/DataAccess.py?request=execute&identifier=SinglePoint&tempAverage=DAILY&parameters=T2M&startDate={date}&endDate={date}&lat={latitude}&lon={longitude}&outputList={output}&user=CFK0e7I4o9tX0PHPVoyPBfkyUNjSNuogV0frLnKl"
    output = "JSON"
    api_request_url = base_url.format(longitude=longitude, latitude=latitude, date=date, output=output.upper())

    # Python Memory Object
    json_response = json.loads(requests.get(api_request_url).content.decode('utf-8'))

    # Select the variable from the response
    response = json_response['features'][0]['properties']['parameter']['T2M']
    
    return response

def closest_point(point, points):
    from scipy.spatial.distance import cdist
    """ Find closest point from a list of points. """
    return points[cdist([point], points).argmin()]

def match_value(df, col1, x, col2, y, col3):
    """ Match value x from col1 row to value in col2. """
    return df[(df[col1] == x) & (df[col3] == y)][col2].values[0]

# def wrapper(point):
#     import pandas as pd
#     global df1 = pd.read_csv('data/large/frame.csv')
#     closest_point(point, list(df1['point']))

def workaround(df):
    df1 = config.df1
    df = config.df2
    
    return match_value(df1, 'point', df[0], 'T2M', df[1], 'DATE')
#       ns.ty = type(df[0])