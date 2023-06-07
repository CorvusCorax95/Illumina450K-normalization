import ui
import os
import shutil

# MAIN CLASS
## We just call the UI functions to build out streamlit app and prepare
# everything
download_dir = "download"
old_zip = "compressed_download.zip"

if os.path.exists(old_zip):
	os.remove(old_zip)

if not os.path.exists(download_dir):
	os.makedirs(download_dir)
else:
	try:
		shutil.rmtree(download_dir)
		os.makedirs(download_dir)
	except OSError as e:
		print("Error: %s - %s." % (e.filename, e.strerror))


ui.make_header()
ui.make_plots()
