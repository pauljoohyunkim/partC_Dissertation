import cv2
import sys
import os

if __name__ == "__main__":

    fourcc = cv2.VideoWriter_fourcc("M","J","P","G")
    fps = 60
    directory = sys.argv[1]
    
    image = cv2.imread(os.path.join(directory, "0.png"))
    rows, columns, channels = image.shape

    writer = cv2.VideoWriter("output.mp4", fourcc, fps, (columns, rows))

    i = 0

    try:
        while (os.path.isfile(os.path.join(directory, f"{i}.png"))):
            image = cv2.imread(os.path.join(directory, f"{i}.png"))
            writer.write(image)
            i += 1
    except KeyboardInterrupt:
        writer.release()
        sys.exit(1)

    writer.release()
