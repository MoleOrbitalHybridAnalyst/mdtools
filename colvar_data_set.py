from read_colvar import *
from glob import glob

class DataSet:
    def __init__(self, image_dirs, label_files, 
                       image_ext = '.dat', label_ext = '', 
                       image_cvname = 'density', label_cvname = 'label',
                       verbose = True):
        self.images = []
        self.labels = []

        for image_dir, label_file in zip(image_dirs, label_files):
            df_label = read_colvar(label_file)
            try:
                self.labels.extend(df_label[label_cvname].values)
            except KeyError as kerr:
                raise KeyError("cannot find " + kerr.__str__() + " in " +\
                        label_file)
            if verbose:
                print("read in label file " + label_file)

            nimage = len(glob(image_dir+"/*"+image_ext))
            # check whether this image_dir contains consistent number 
            # of images with the number of lines in this label file:
            if nimage != len(df_label):
                raise Exception("inconsistent " + image_dir + " and " +\
                        label_file)
            for idx in range(nimage):
                df_image = read_colvar(image_dir+"/"+str(idx)+image_ext)
                try:
                    self.images.append(df_image[image_cvname].values)
                except KeyError as kerr:
                    raise KeyError("cannot find " + kerr.__str__() + " in " + \
                            image_dir+"/"+str(idx)+image_ext)

                # check whether each image has the same size:
                if len(self.images[0]) != len(df_image):
                    raise Exception("inconsistent image size of " +\
                            image_dir+"/"+str(idx)+image_ext)
            if verbose:
                print("read in", nimage, "images from " + image_dir)

        self.images = np.array(self.images)
        self.labels = np.array(self.labels)

        self.ndata = len(self.labels)
        if verbose:
            print("total number of data points read is", self.ndata)
