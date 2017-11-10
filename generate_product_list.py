#!/usr/bin/env python

# TODO improve docstrings
# TODO improve format
# TODO improve help text

from __future__ import absolute_import, division

import re
import urllib
import cachetools
from datetime import datetime
from pathlib import Path
import csv as lcsv
import itertools

import click

import geojson as gjson

from fastkml import kml as fkml
from fastkml.geometry import Geometry, Polygon
from pyproj import Proj


THREDDS_PRODUCT_LIST = [
    'ls5_nbar_albers',
    'ls7_nbar_albers',
    'ls8_nbar_albers',
    'ls5_nbart_albers',
    'ls7_nbart_albers',
    'ls8_nbart_albers',
    'ls5_pq_albers',
    'ls7_pq_albers',
    'ls8_pq_albers'
]

# 1986 is the start of our landsat 5 nbar collection, run until the present
NOW = datetime.now()
YEAR_RANGE = list(range(1986, NOW.year))


def get_grid_perimeter(coord_y, coord_x, coordinate_system='epsg:3577'):
    scaling_factor = 100000
    proj = Proj(init=coordinate_system)

    top_left = tuple([(coord_y + 1) * scaling_factor, coord_x * scaling_factor])
    bottom_right = tuple([coord_y * scaling_factor, (coord_x + 1) * scaling_factor])

    return list(itertools.chain(
        (proj(lon, top_left[1], inverse=True) for lon in range(top_left[0], bottom_right[0], -1000)),
        (proj(bottom_right[0], lat, inverse=True) for lat in range(top_left[1], bottom_right[1], 1000)),
        (proj(lon, bottom_right[1], inverse=True) for lon in range(bottom_right[0], top_left[0], 1000)),
        (proj(top_left[0], lat, inverse=True) for lat in range(bottom_right[1], top_left[1], -1000)),
        [proj(*top_left, inverse=True)]
    ))

def convert_grid_coords(polygon_bounding_points, coordinate_system='epsg:3577'):
    proj = Proj(init=coordinate_system)

    corrected_sets = []
    for coordinate_set in polygon_bounding_points:
        corrected_sets.append([proj(p[0], p[1], inverse=True) for p in coordinate_set])

    return corrected_sets

class LandSatPathRowDCField(object):
    CACHE_SIZE = 6100  # Set arbitrarily; based on observations in landsat 5
    NAME = 'wrs2_path_row'
    
    def __init__(self, datacube):
        self.cache = cachetools.LRUCache(self.CACHE_SIZE)
        self.datacube = datacube

    def get_value(self, dss_record):
        path_row = self.cache.get(dss_record.metadata.time.begin)

        if not path_row:
            pr_coords = (
                self.datacube.index.datasets.get(dss_record.id, include_sources=True)
                .metadata.sources['0']['image']['satellite_ref_point_start']
            )
            path_row = '_'.join((str(pr_coords['x']).zfill(3), str(pr_coords['y']).zfill(3)))
            self.cache.update({
                dss_record.metadata.time.begin: path_row
            })

        return path_row


class ObservationDateDCField(object):
    NAME = 'observation_date'

    @classmethod
    def get_value(cls, dss_record):
        """get_value: Returns the observation date for the record

        :param dss_record: A dataset record
        """
        try:
            return dss_record.time
        except AttributeError:
            return None


class CreationDateDCField(object):
    NAME = 'creation_date'

    @classmethod
    def get_value(cls, dss_record):
        """get_value: Returns the creation date of the record

        :param dss_record: A dataset record
        """
        try:
            return dss_record.metadata.creation_dt
        except AttributeError:
            return None


class CoordinateSetDCField(object):
    NAME = 'coordinate_set'
    REGEX_PATTERN = re.compile('file://(?:/[^/]*){6}/([-0-9]+_[-0-9]+)')

    @classmethod
    def get_value(cls, dss_record):
        """get_value: Returns the coordinate for the record

        :param dss_record: A dataset record
        """
        try:
            return cls.REGEX_PATTERN.match(dss_record.local_uri).groups()[0]
        except AttributeError:
            return None


class PolygonPointBounds(object):
     NAME = 'polygon_point_bounds'

     @classmethod
     def get_value(cls, dss_record):
         try:
             return dss_record.metadata.grid_spatial['valid_data']['coordinates']
         except AttributeError:
             return None


class NetCDFSliceDCField(object):
    NAME = 'netcdf_slice'

    THREDDS_SERVER = 'http://dapds00.nci.org.au/thredds/ncss'
    THREDDS_TS_FORMAT = '%Y-%m-%dT%H:%M:%SZ'
    LOCAL_URI_PREFIX = 'file:///g/data/'

    def __init__(self, bands):
        self.bands = bands

    def get_value(self, dss_record):
        """get_value: Returns the netcdf_slice url for the record

        :param dss_record: A dataset record
        """
        params = {
            'timeStride': 1,
            'horizStride': 1,
            'var': self.bands,
            'time_start': dss_record.time.begin.strftime(self.THREDDS_TS_FORMAT),
            'time_end': dss_record.time.end.strftime(self.THREDDS_TS_FORMAT)
        }

        return '{server}/{file_path}?{params}'.format(
            server=self.THREDDS_SERVER,
            file_path=dss_record.local_uri.replace(self.LOCAL_URI_PREFIX, ''),
            params=urllib.parse.urlencode(params, doseq=True)  # doseq handles arrays
        )


class GeoTiFFDCField(object):

    THREDDS_SERVER = 'http://dapds00.nci.org.au/thredds/wcs'
    THREDDS_TS_FORMAT = '%Y-%m-%dT%H:%M:%SZ'
    LOCAL_URI_PREFIX = 'file:///g/data/'

    def __init__(self, band):
        self.NAME = 'band_{}'.format(band)
        self.band = band

    def get_value(self, dss_record):
        """get_value: Returns the geotiff url for the record for a given band

        :param dss_record: A dataset record
        """
        params = {
            'service': 'WCS',
            'version': '1.0.0',
            'request': 'GetCoverage',
            'format': 'GeoTIFF',
            'coverage': self.band,
            'time': dss_record.time.begin.strftime(self.THREDDS_TS_FORMAT),
            self.band: '100.0'
        }

        return '{server}/{file_path}?{params}'.format(
            server=self.THREDDS_SERVER,
            file_path=dss_record.local_uri.replace(self.LOCAL_URI_PREFIX, ''),
            params=urllib.parse.urlencode(params, doseq=True)
        )


class SpatialReferenceDCField(object):
    NAME = 'spatial_reference'

    @classmethod
    def get_value(cls, dss_record):
        """get_value: Returns the spatial reference for the dataset record

        :param dss_record: A dataset record
        """
        try:
            return dss_record.crs
        except AttributeError:
            return None


class UUIDDCField(object):
    NAME = 'uuid'

    @classmethod
    def get_value(cls, dss_record):
        """get_value: Returns the uuid for the dataset record

        :param dss_record: A dataset record
        """
        try:
            return dss_record.metadata.id
        except AttributeError:
            return None


class CSVWriter(object):
    TS_OUTPUT_FMT = '%y-%m-%dT:%H:%M:%SZ'

    def __init__(self, outfile, fields):
        self.outfile = Path(outfile)
        self.fields = fields

        self._written_data = False
        self._fdesc = None
        self._writer = None

    def _write_headers(self):
        self._writer.writerow((map(lambda f: f.NAME, self.fields)))

    def __enter__(self):
        self.outfile.absolute().parent.mkdir(parents=True, exist_ok=True)
        self._fdesc = self.outfile.open(mode='w')
        self._writer = lcsv.writer(self._fdesc)
        self._write_headers()

        return self

    @classmethod
    def convert_to_str(cls, field_name, field_value):
        if field_name == 'observation_date':
            return field_value[0].strftime(cls.TS_OUTPUT_FMT)

        if isinstance(field_value, datetime):
            return field_value.strftime(cls.TS_OUTPUT_FMT)

        return str(field_value)

    def write(self, record_set):
        self._written_data = True
        self._writer.writerow(map(lambda f: self.convert_to_str(f.NAME, record_set[f.NAME]), self.fields))

    # TODO clean this up
    def __exit__(self, *args, **kwargs):
        if self._fdesc:
            self._fdesc.close()
        self._fdesc = None

        if not self._written_data:
            self.outfile.unlink()


class GeoJSONWriter(object):

    TS_OUTPUT_FMT = '%Y-%m-%dT:%H:%M:%SZ'
    
    def __init__(self, bands, projection, outfile):
        self.features = {}
        self.bands = bands
        self.projection = projection
        self.outfile = Path(outfile)

    def get_feature(self, coord):
        if coord in self.features:
            return self.features[coord]

        geo_polygon = gjson.Polygon(
            coordinates=[get_grid_perimeter(*[int(c) for c in coord.split('_')])],
            crs=gjson.crs.Named(properties={'crs':'EPSG:4236'}),
            validate=True
        )

        geo_feature = gjson.Feature(
            id=coord,
            geometry=geo_polygon,
            properties={}
        )

        geo_feature.properties.update({
            'coord': coord,
            'coord_spatial_reference': self.projection,
            'dataset_info': {}
        })

        self.features[coord] = geo_feature

        return geo_feature

    def write(self, record_set):
        coord_feature = self.get_feature(record_set.get('coordinate_set'))

        dataset_dict = {
            'uuid': record_set.get('uuid'),
            'wrs2_path_row': record_set.get('wrs2_path_row'),
            'netcdf_slice': record_set.get('netcdf_slice'),
            'observation_date': record_set.get('observation_date')[0].strftime(self.TS_OUTPUT_FMT),
            'creation_date': record_set.get('creation_date')
        }

        for band in self.bands:
            dataset_dict['band_{}'.format(band)] = record_set.get('band_{}'.format(band))

        coord_feature['properties']['dataset_info'][dataset_dict['observation_date']] = dataset_dict

    def save(self):
        feature_collection = gjson.FeatureCollection(
            features=list(self.features.values())
        )
        self.outfile.parent.mkdir(parents=True, exist_ok=True)
        self.outfile.write_text(gjson.dumps(feature_collection))


class KMLWriter(object):
    TS_OUTPUT_FMT = '%Y-%m-%dT:%H:%M:%SZ'
    NS = '{http://www.opengis.net/kml/2.2}'  # pylint: disable=invalid-name

    def __init__(self, outfile, fields, projection, bands):
        self.outfile = Path(outfile)
        self.fields = fields
        self.projection = projection
        self.bands = bands
        self.schema_set = {}

        self.kml_root = None
        self.coord_container = None
        self.coord_cache = {}

    def __enter__(self):
        self.outfile.absolute().parent.mkdir(parents=True, exist_ok=True)
        self._write_headers()

    def register_schemas(self, root_document):
        # Attributes at the document level
        document_metadata_schema = fkml.Schema(
            ns=self.NS,
            id='doc-metadata',
            fields=[
                {'type': 'string', 'name': 'gen_date', 'displayName': 'Generated date'},
                {'type': 'string', 'name': 'spatialref', 'displayName': 'Spatial Reference'}
            ]
        )

        root_document.append_schema(document_metadata_schema)
        self.schema_set['document_metadata'] = '#doc-metadata'

        # Attributes on the coordinate container
        coordinate_metadata_schema = fkml.Schema(
            ns=self.NS,
            id='coord-metadata',
            fields=[
                {'type': 'string', 'name': 'coord', 'displayName': 'coordinate'}
            ]
        )

        root_document.append_schema(coordinate_metadata_schema)
        self.schema_set['coordinate_metadata'] = '#coord-metadata'

        # Attributes on the placemark
        # Note that for QGIS all the placemarks in a container should share the same schema
        placemark_fields = [
            {'type': 'string', 'name': 'obs_date', 'displayName': 'Observation date'},
            {'type': 'string', 'name': 'uuid', 'displayName': 'UUID'},
            {'type': 'string', 'name': 'wrs2_path_row', 'displayName': 'WRS2 Path Row'},
            {'type': 'string', 'name': 'netcdf_slice', 'displayName': 'NetCDF3 Slice'}
        ]

        for band in self.bands:
            placemark_fields.append(
                {'type': 'string', 'name': 'band_{}'.format(band), 'displayName': band.capitalize()}
            )

        placemark_metadata_schema = fkml.Schema(
            ns=self.NS,
            id='placemark-metadata',
            fields=placemark_fields
        )
        self.schema_set['placemark_metadata'] = '#placemark-metadata'
        root_document.append_schema(placemark_metadata_schema)
 

    def configure_document_root(self):
        root_document = fkml.Document(ns=self.NS, id='root-doc')
        self.register_schemas(root_document)

        self.kml_root = fkml.KML(ns=self.NS)
        self.kml_root.append(root_document)

        self.coord_container = fkml.Folder(
            ns=self.NS,
            id='coordinates',
            name='coordinates'
        )

        self.coord_container.extended_data = fkml.SchemaData(
            ns=self.NS,
            schema_url=self.schema_set['document_metadata'],
            data=[
               {'name': 'gen_date', 'value': NOW.strftime(self.TS_OUTPUT_FMT)},
               {'name': 'spatial_reference', 'value': self.projection}
            ]
        )

        root_document.append(self.coord_container)

    def get_coord_folder(self, coord):
        folder = self.coord_cache.get(coord)
        if folder:
            return folder

        folder = fkml.Folder(ns=self.NS, name=coord, id='f-{}'.format(coord))

        folder.extended_data = fkml.SchemaData(
            ns=self.NS,
            schema_url=self.schema_set['coordinate_metadata'],
            data=[{'name': 'coord', 'value': coord}]
        )

        visible_marker = fkml.Placemark(ns=self.NS, name=coord, id='pl-{}'.format(coord))

        # Dummy data required by QGIS to display all columns
        dummy_data = [
            {'name': 'uuid', 'value': None},
            {'name': 'wrs2_path_row', 'value': None},
            {'name': 'netcdf_slice', 'value': None},
        ]

        for band in self.bands:
            dummy_data.append({'name': 'band_{}'.format(band), 'value': None})

        visible_marker.extended_data = fkml.SchemaData(
            ns=self.NS,
            schema_url=self.schema_set['placemark_metadata'],
            data=dummy_data
        )

        visible_marker.geometry = Geometry(
            ns=self.NS,
            geometry=Polygon(
                get_grid_perimeter(*[int(i) for i in coord.split('_')])
            )
        )

        visible_marker.description = ''
        folder.append(visible_marker)

        self.coord_container.append(folder)
        self.coord_cache[coord] = folder

        return folder

    def write(self, record_set):
        coord_folder = self.get_coord_folder(record_set.get('coordinate_set'))
        coord_mark = fkml.Placemark(ns=self.NS, name=record_set.get('observation_date')[0].strftime(self.TS_OUTPUT_FMT))
        coord_mark.begin, coord_mark.end = record_set.get('observation_date')
        coord_mark.visibility = 0

        coord_mark.geometry = Geometry(
            ns=self.NS,
            geometry=Polygon(
                *convert_grid_coords(record_set.get('polygon_point_bounds'))
            )
        )
         
        coord_data = [
            {'name': 'uuid', 'value':  record_set.get('uuid')},
            {'name': 'wrs2_path_row', 'value': record_set.get('wrs2_path_row')},
            {'name': 'netcdf_slice', 'value': record_set.get('netcdf_slice')}
        ]

        for band in self.bands:
            coord_data.append({'name': 'band_{}'.format(band), 'value': record_set.get('band_{}'.format(band))})

        coord_mark.extended_data = fkml.SchemaData(
            ns=self.NS,
            schema_url=self.schema_set['placemark_metadata'],
            data=coord_data
        )
        coord_folder.append(coord_mark)

    def save(self):
        self.outfile.parent.mkdir(parents=True, exist_ok=True)
        self.outfile.write_text(self.kml_root.to_string(prettyprint=True))


class DCReader(object):

    def __init__(self, datacube, query_params=None, fields=None):
        self.dc_result_set = None
        self.datacube = datacube
        self.query_params = query_params
        self.fields = fields

    def __iter__(self):
        self.dc_result_set = self.datacube.index.datasets.search(
            **self.query_params
        )
        return self

    def __next__(self):
        record = next(self.dc_result_set)
        return dict(zip(
            map(lambda f: f.NAME, self.fields),
            map(lambda f: f.get_value(record), self.fields)
        ))


@click.group()
@click.option('--year', '-y', multiple=True, type=int, default=YEAR_RANGE)
@click.option('--product', '-p', multiple=True, type=str,
              default=THREDDS_PRODUCT_LIST)
@click.option('--outdir', '-o', type=str, default='product_list')
@click.pass_context
def cli(context, year, product, outdir):
    context.obj['YEARS'] = year
    context.obj['PRODUCTS'] = product
    context.obj['OUTDIR'] = outdir

@cli.command()
@click.pass_context
def geojson(context):
    import datacube

    outdir = Path(context.obj['OUTDIR']).absolute()
    dc = datacube.Datacube()

    for product in context.obj['PRODUCTS']:
        product_info = dc.index.products.get_by_name(product)
        bands = product_info.measurements.keys()
        projection = str(product_info.grid_spec.crs)

        fields = [
            ObservationDateDCField,
            CreationDateDCField,
            CoordinateSetDCField,
            LandSatPathRowDCField(dc),
            NetCDFSliceDCField(bands),
            UUIDDCField,
        ]

        # Add a field for each band
        for band in bands:
            fields.append(GeoTiFFDCField(band))

        for year in context.obj['YEARS']:
            outfile = outdir / product / '{}__{}.json'.format(
                year, NOW.strftime('%Y-%m-%d')
            )
            writer = GeoJSONWriter(bands, projection, outfile)
            query_params = {
                'product': product,
                'time': datacube.model.Range(
                    datetime(year, 1, 1), datetime(year + 1, 1, 1)
                )
            }

            for record in DCReader(dc, query_params, fields):
                writer.write(record)
                break

            writer.save()


@cli.command()
@click.pass_context
def kml(context):
    import datacube

    outdir = Path(context.obj['OUTDIR']).absolute()
    dc = datacube.Datacube()

    for product in context.obj['PRODUCTS']:
        product_info = dc.index.products.get_by_name(product)
        bands = product_info.measurements.keys()
        projection = str(product_info.grid_spec.crs)

        fields = [
            ObservationDateDCField,
            CreationDateDCField,
            CoordinateSetDCField,
            LandSatPathRowDCField(dc),
            NetCDFSliceDCField(bands),
            PolygonPointBounds,
            UUIDDCField,
        ]

        # Add a field for each band
        for band in bands:
            fields.append(GeoTiFFDCField(band))

        for year in context.obj['YEARS']:
            writer = KMLWriter(
                outfile=outdir / product / '{}__{}.kml'.format(
                    year, NOW.strftime('%Y-%m-%d')
                ),
                fields=fields,
                projection=projection,
                bands=bands
            )
            writer.configure_document_root()
            query_params = {
                'product': product,
                'time': datacube.model.Range(
                    datetime(year, 1, 1), datetime(year + 1, 1, 1)
                )
            }

            for record in DCReader(dc, query_params, fields):
                writer.write(record)

            writer.save()


@cli.command()
@click.pass_context
def csv(context):
    import datacube

    dc = datacube.Datacube()
    outdir = Path(context.obj['OUTDIR']).absolute()

    for product in context.obj['PRODUCTS']:
        product_info = dc.index.products.get_by_name(product)
        bands = product_info.measurements.keys()
        projection = str(product_info.grid_spec.crs)

        fields = [
            ObservationDateDCField,
            UUIDDCField,
            CreationDateDCField,
            SpatialReferenceDCField,
            CoordinateSetDCField,
            LandSatPathRowDCField(dc),
            NetCDFSliceDCField(bands),
        ]

        # Add a field for each band
        for band in bands:
            fields.append(GeoTiFFDCField(band))

        for year in context.obj['YEARS']:
            outfile = Path(outdir / product / '{}__{}.csv'.format(
                year, NOW.strftime('%Y-%m-%d'))
            )

            query_params = {
                'product': product,
                'time': datacube.model.Range(
                    datetime(year, 1, 1), datetime(year + 1, 1, 1)
                )
            }
            with CSVWriter(outfile, fields) as writer:
                for record in DCReader(dc, query_params, fields):
                    writer.write(record)


if __name__ == '__main__':
    cli(obj={})
