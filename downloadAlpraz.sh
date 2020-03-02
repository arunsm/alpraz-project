# A python script to download all the DICOMs associated with project "ALPRAZ_805556" from XNAT

import sys
import re
import os
sys.path.append("/data/joy/BBL/applications/python/lib/python2.6/site-packages/pyxnat-1.0.0.0-py2.6.egg/")
from pyxnat import Interface
xnat_config = "/home/amahad/.xnat.cfg"
central = Interface(config=xnat_config)
import urllib3
urllib3.disable_warnings()

def DownloadLoop():
		
	write_dir = "/data/joy/BBL/studies/alpraz/rawData2/"; # directory to download DICOMs to

	# creating constraints for querying data - REPLACE WITH XnatHelper function
	constraints = [('bbl:Sequence/PROJECT','LIKE','%'+'ALPRAZ'+'%'),'AND'];		
	seqs = central.select('bbl:Sequence',['bbl:Sequence/QLUX_QLUXNAME','bbl:Sequence/IMAGESCAN_ID','bbl:Sequence/SUBJECT_ID', 'bbl:Sequence/imageSession_ID', 'bbl:Sequence/date', 'bbl:Sequence/PROTOCOL', 'bbl:Sequence/PROJECT','bbl:Sequence/MR_SERIESDESCRIPTION','bbl:Sequence/MR_IMAGEORIENTATIONPATIENT','bbl:Sequence/QLUX_MATCHED']).where(constraints);

	for line in seqs:
		BBL_ID = line.get('subject_id');
	        seqname = line.get('qlux_qluxname');
	        scan_ID = line.get('session_id');
	        proj_name = line.get('project');
	        scandate = line.get('date');
	        seq_id = line.get('imagescan_id');
	        imgorient = line.get('mr_imageorientationpatient');
	        formname = line.get('mr_seriesdescription');
	        match = line.get('qlux_matched');
		download_rest='/projects/' + proj_name + '/subjects/' + BBL_ID + '/experiments/' + scan_ID + '/scans/' + seq_id + '/resources/DICOM/files';
		f = central.select(download_rest)
	
		for i in f:
			cwd = write_dir+"%s/%s/" % (BBL_ID, scan_ID)
			if not os.path.exists(cwd):
				os.makedirs(cwd)
            		print "Downloading file ", (str(cwd+i._urn))
			try:
            			if not os.path.isfile(str(cwd+i._urn)) and not os.path.isfile(str(cwd+i._urn).split(".gz")[0]):
       	       				i.get(cwd+i._urn)
			except:
				print "Could not process participant %s, scan %s" % (BBL_ID, scan_ID)
				with open('/data/joy/BBL/studies/alpraz/XNATDownloadErrors.txt', 'a') as f:
					f.write("%s, %s/" % (BBL_ID, scan_ID))

def main():
	DownloadLoop()

if __name__ == "__main__":
	main()

