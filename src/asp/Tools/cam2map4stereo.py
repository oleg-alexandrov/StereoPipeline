#!/usr/bin/env python
# __BEGIN_LICENSE__
#  Copyright (c) 2009-2013, United States Government as represented by the
#  Administrator of the National Aeronautics and Space Administration. All
#  rights reserved.
#
#  The NGT platform is licensed under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance with the
#  License. You may obtain a copy of the License at
#  http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing, software
#  distributed under the License is distributed on an "AS IS" BASIS,
#  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#  See the License for the specific language governing permissions and
#  limitations under the License.
# __END_LICENSE__

from __future__ import print_function
import os, optparse, subprocess, sys, tempfile

# The path to the ASP python files.
basepath    = os.path.abspath(sys.path[0])
pythonpath  = os.path.abspath(basepath + '/../Python')  # for dev ASP
libexecpath = os.path.abspath(basepath + '/../libexec') # for packaged ASP
sys.path.insert(0, basepath) # prepend to Python path
sys.path.insert(0, pythonpath)
sys.path.insert(0, libexecpath)

from asp_stereo_utils import get_asp_version
from asp_system_utils import mkdir_p

import asp_system_utils
asp_system_utils.verify_python_version_is_supported()

def man(option, opt, value, parser):
    print (parser.usage, file=sys.stderr)
    print ('''\
This program takes similar arguments as the ISIS3 cam2map program,
but takes two input images.  With no arguments, the program determines
the minimum overlap of the two images, and the worst common resolution,
and then map-projects the two images to this identical area and resolution.

Using the --lat and --lon options allows you to select portions of
the image, ideal for small initial test runs or cropping out just
that portion of the image you need.
''', file=sys.stderr)
    
    sys.exit()

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

class MapExists(Exception):
    def __init__(self, msg):
        self.msg = msg

class ImgInfo:
    def __init__( self ):
        self.resolution = None
        self.minlat = None
        self.maxlat = None
        self.minlon = None
        self.maxlon = None

def mapfile( cube, prefix, suffix ):
    '''???'''
    items = cube.split('.')
    mapname = '.'.join([items[0], suffix, 'cub'])

    # If the file exists, or it starts with a dot (say ../file.cub)
    if( os.path.exists(mapname) or mapname[0] == "."):
        items.insert(len(items)-1, suffix)
        mapname = '.'.join(items)

    # Make somedir/file.cub into myprefix-file.map.cub
    if prefix != "":
        mapname = prefix + "-" + os.path.basename(mapname)

    # Ensure that any output directory exists
    dirname = os.path.dirname(mapname)
    mkdir_p(dirname)

    if( os.path.exists(mapname) ):
        raise MapExists( mapname )

    return mapname

def camrange( cube, lonGroupName, options ):
    '''run camrange on cube'''
    info = ImgInfo();
    try:
        # Run the camrange function to write to a temp file
        tmpfile = tempfile.NamedTemporaryFile(dir='.', prefix="camrange",
                                              suffix=".txt", mode="w")
        cmd = 'camrange from= '+ cube +' to= '+ tmpfile.name

        print('Running ' + cmd)
        result = subprocess.call(cmd, shell=True)
        if not result == 0:
            raise Exception('ProcessError', 'Non zero return code ('+str(result)+')')
        os.system("cat %s" % tmpfile.name )

        # Path to the ISIS getkey tool
        getkey_path = os.environ['ISISROOT']+'/bin/getkey'

        # By default we get the longitude from the Universal field,
        #   but if the user specified a longitude region we need to match it.
        lonDomain    = '360'
        lonDirection = 'PositiveEast'
        if options.map:
            # Read the pertinent values from the map file
            lonDomainIn    = subprocess.Popen([getkey_path, 'from= '+options.map, "grpname= Mapping", "keyword= LongitudeDomain"   ], stdout=subprocess.PIPE,stderr=subprocess.PIPE, universal_newlines=True).communicate()[0].strip()
            lonDirectionIn = subprocess.Popen([getkey_path, 'from= '+options.map, "grpname= Mapping", "keyword= LongitudeDirection"], stdout=subprocess.PIPE,stderr=subprocess.PIPE, universal_newlines=True).communicate()[0].strip()
            if lonDomainIn:
                lonDomain = lonDomainIn
            if lonDirectionIn:
                lonDirection = lonDirectionIn
            # Determine the correct lon group
            if (lonDomain == '360') and (lonDirection == 'PositiveWest'):
                lonGroupName = 'PositiveWest360'
            if (lonDomain == '180') and (lonDirection == 'PositiveEast'):
                lonGroupName = 'PositiveEast180'
            if (lonDomain == '180') and (lonDirection == 'PositiveWest'):
                lonGroupName = 'PositiveWest180'
        lonGroupString = 'grpname= ' + lonGroupName
        fromString     = "from= "+ tmpfile.name

        # extract information
        info.resolution = subprocess.Popen([getkey_path, fromString, "grpname= PixelResolution",      "keyword= Lowest"          ], stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].strip()
        info.minlat     = subprocess.Popen([getkey_path, fromString, "grpname= UniversalGroundRange", "keyword= MinimumLatitude" ], stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].strip()
        info.maxlat     = subprocess.Popen([getkey_path, fromString, "grpname= UniversalGroundRange", "keyword= MaximumLatitude" ], stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].strip()
        info.minlon     = subprocess.Popen([getkey_path, fromString, lonGroupString,                  "keyword= MinimumLongitude"], stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].strip()
        info.maxlon     = subprocess.Popen([getkey_path, fromString, lonGroupString,                  "keyword= MaximumLongitude"], stdout=subprocess.PIPE, universal_newlines=True).communicate()[0].strip()

    except OSError as e:
        print ("Execution failed:", e, file=sys.stderr)
        raise e

    return info

def main():
    try:
        try:
            usage = "usage: cam2map4stereo.py [--help][--manual][--map mapfile][--pixres CAMERA|MAP|MPP|PPD][--resolution float][--interp NN|BI|CC][--lat min:max][--lon min:max][--prefix string] [--suffix string] image1.cub image2.cub\n" + get_asp_version()
            parser = optparse.OptionParser(usage=usage)
            parser.set_defaults(dryrun=False)
            parser.set_defaults(pixres='MPP')
            parser.set_defaults(prefix='')
            parser.set_defaults(suffix='map')
            parser.add_option("--manual", action="callback", callback=man,
                              help="Read the manual.")
            parser.add_option("-m", "--map", dest="map",
                              help="The mapfile to use for cam2map.")
            parser.add_option("-p", "--pixres", dest="pixres",
                              help="The pixel resolution mode to use for cam2map.")
            parser.add_option("-r", "--resolution", dest="resolution",
                              help="Resolution of final map for cam2map.")
            parser.add_option("-i", "--interp", dest="interp",
                              help="Pixel interpolation scheme for cam2map.")
            parser.add_option("-a", "--lat", dest="lat",
                              help="Latitude range for cam2map where LAT is minlat:maxlat.")
            parser.add_option("-o", "--lon", dest="lon",
                              help="Longitude range for cam2map where LON is minlon:maxlon.")
            parser.add_option("--prefix", dest="prefix",
                              help="Make all output files use this prefix. Default: no prefix.")
            parser.add_option("-s", "--suffix", dest="suffix",
                              help="Suffix that gets inserted in the output file names.")
            parser.add_option("-n", "--dry-run", dest="dryrun",
                              action="store_true",
                              help="Make calculations, but print cam2map command, but don't actually run it.")

            (options, args) = parser.parse_args()

            if not args: parser.error("need .cub files")


        except optparse.OptionError as msg:
            raise Usage(msg)

        # Call camrange to get bounds. If the lon bounds come as [0, 360] then try a different group.
        success = False
        lonGroupNames = ['UniversalGroundRange',  'PositiveEast180',
                         'PositiveWest180', 'PositiveWest360']
        for lonGroupName in lonGroupNames:
            image1 = camrange( args[0], lonGroupName, options)
            image2 = camrange( args[1], lonGroupName, options)
            lonRange = float(max(image1.maxlon, image2.maxlon)) - \
                       float(min(image1.minlon, image2.minlon))
            if lonRange >= 360.0:
                print("With lonGroupName = " + lonGroupName + " found a longitude range of " + \
                      str(lonRange) + ". Will try the next group from: " + (" ".join(lonGroupNames)))
            else:
                success = True
                break

        if not success:
            raise Exception('ProcessError', 'The longitude range keeps on coming as encompasing the entire planet. Something is wrong.')
            
        mapout = ImgInfo()

        mapout.minlat     = max( image1.minlat,     image2.minlat )
        mapout.maxlat     = min( image1.maxlat,     image2.maxlat )
        mapout.minlon     = max( image1.minlon,     image2.minlon )
        mapout.maxlon     = min( image1.maxlon,     image2.maxlon )
        mapout.resolution = max( image1.resolution, image2.resolution )

        if( options.resolution ): mapout.resolution = options.resolution
        if( options.lat ):
            latrange = options.lat.split(':')
            if latrange[0]: mapout.minlat = latrange[0]
            if latrange[1]: mapout.maxlat = latrange[1]
        if( options.lon ):
            lonrange = options.lon.split(':')
            if lonrange[0]: mapout.minlon = lonrange[0]
            if lonrange[1]: mapout.maxlon = lonrange[1]

        # call cam2map with the arguments

        cam2map = ['cam2map', 'from=' + args[0], 'to='
                   + mapfile( args[0], options.prefix, options.suffix )]

        if( options.map ):
            cam2map.append( 'map=' + options.map )
        if( options.interp):
            cam2map.append( 'interp=' + options.interp)

        cam2map.append( 'pixres=' + options.pixres )
        cam2map.append( 'resolution=' + mapout.resolution )

        cam2map.append( 'defaultrange=MAP' )
        cam2map.append( 'minlat=' + mapout.minlat )
        cam2map.append( 'maxlat=' + mapout.maxlat )
        cam2map.append( 'minlon=' + mapout.minlon )
        cam2map.append( 'maxlon=' + mapout.maxlon )

        # Run for first image

        # Need to put these together to keep ISIS from calling the GUI
        cam2map_cmd = ' '.join(cam2map)
        print(cam2map_cmd)
        if( not options.dryrun ):
            code = subprocess.call(cam2map_cmd, shell=True)
            if not code == 0:
                raise Exception('ProcessError', 'Non zero return code ('+str(code)+')')

        # Run for second image
        cam2map[1] = 'from=' + args[1]
        cam2map[2] = 'to='+ mapfile( args[1], options.prefix, options.suffix )
        cam2map_cmd = ' '.join(cam2map)
        print(cam2map_cmd)
        if( not options.dryrun ):
            code = subprocess.call(cam2map_cmd, shell=True)
            if not code == 0:
                raise Exception('ProcessError', 'Non zero return code ('+str(code)+')')

        print("Finished")
        return 0

    except Usage as err:
        print (err.msg, file=sys.stderr)
        return 2

    except MapExists as e:
        print ('The file '+ e.msg +' already exists, delete first.', file=sys.stderr)
        return 3

    # # To more easily debug this program, comment out this catch block.
    # except Exception as err:
    #     sys.stderr.write( str(err) + '\n' )
    #     return 1

if __name__ == "__main__":
    sys.exit(main())
