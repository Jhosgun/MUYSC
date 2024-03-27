import os
import requests


class TopographyData:

    def __init__(self, region_points):
        self.region_points = region_points

    def format_coordinate(self, value, positive, negative):
        prefix = positive if value >= 0 else negative
        if value >= 0:
            abs_value = abs(int(value))
        else:
            abs_value = abs(int(value))+1
        width = 2 if prefix in ('N', 'S') else 3
        return f"{prefix}{abs_value:0{width}}"

    def download_data(self):
        if len(self.region_points) != 4:
            raise ValueError("Expected 4 region points")

        lat1, lat2, lon1, lon2 = self.region_points

        lat1_str = self.format_coordinate(lat1, "N", "S")
        lat2_str = self.format_coordinate(lat2, "N", "S")
        lon1_str = self.format_coordinate(lon1, "E", "W")
        lon2_str = self.format_coordinate(lon2, "E", "W")

        files = []

        if lat1_str == lat2_str and lon1_str == lon2_str:
            files.append(f"{lat1_str}{lon1_str}.SRTMGL1.hgt.zip")
        else:
            for lat in (lat1_str, lat2_str):
                for lon in (lon1_str, lon2_str):
                    files.append(f"{lat}{lon}.SRTMGL1.hgt.zip")

        return files

    def download_file(self, download_url):
        req = requests.get(download_url)
        filename = req.url[download_url.rfind('/')+1:]

        with open(filename, 'wb') as f:
            for chunk in req.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

    def download_to_path(self, path):
        os.chdir(path)

        for filename in self.download_data():
            download_link = f"http://step.esa.int/auxdata/dem/SRTMGL1/{filename}"
            self.download_file(download_link)



