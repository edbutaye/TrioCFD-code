ffmpeg -framerate 24 -i img%04d.png -vcodec libxvid -q:v 1 movie.mp4
visit -movie  -format png -sessionfile visit0000.session -framestep 1 -start 1 -end 5000 -output "plot/fig"
visit -movie  -format png -sessionfile visit0000.session -framestep 1 -start 1 -end 5000 -output "plot/fig"

ffmpeg -r 60  -i figure%04d.png -vcodec libx264 -crf 25 notcl.mp4
ffmpeg -r 25 -start_number 1 -i figure%04d.png  -vf "format=yuv444p, scale=1433:922" out.avi
