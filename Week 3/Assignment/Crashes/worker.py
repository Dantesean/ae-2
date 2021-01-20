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