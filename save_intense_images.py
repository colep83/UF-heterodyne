# Will it be quicker to save images or save arrays

#what needs to be imported?
import os
import PySpin
import numpy
import argparse
from numpy import array

parser = argparse.ArgumentParser(description="Initialize HPC settings")

parser.add_argument('-exp', help='The exposure time of the camera in micro seconds. Bounds: 5 - 25000.',
                    type = float, default=4000)
parser.add_argument('-gain', help='The gain of the camera in dB. Bounds: 0 - 47.',
                    type = float, default = 0)
parser.add_argument('-NI', help='The number of images to be collected by the camera.',
                    type = int, default = 400)
parser.add_argument('-PSA', help='Phase shift algorithm to use. Currently only Novak', default = 'Novak')

args = parser.parse_args()

exp = args.exp
gain = args.gain
NUM_IMAGES = args.NI
phase = args.PSA #this variable currently serves no purpuse

Folder_Name = input("Enter the name of the folder you wish to save the images too\nthis will also be the beginning of the file name\n(Be careful to not overwrite another folder within\nthe Data directory by using the same name): ")

#TEMPORARY CHANGE COME BACK TO THIS!!!!!!!!!!
Data_dir='C:/Users/physics-svc-fulda/Desktop/HetPhaseCam/Data'
data_path = os.path.join(Data_dir,Folder_Name)
os.mkdir(data_path)
print('Directory created:', Folder_Name)

def configure_exp_gain(cam,exp,gain):
    #this function sets the exposure and gain as well as explicitly setting the ADC to 12 bit
    try:
        result = True
        nodemap = cam.GetNodeMap()

        # setting exposure
        if cam.ExposureAuto.GetAccessMode() != PySpin.RW:
            print('Unable to disable automatic exposure. Aborting...')
            return False

        cam.ExposureAuto.SetValue(PySpin.ExposureAuto_Off)
        #print('Automatic exposure disabled...')

        if cam.ExposureTime.GetAccessMode() != PySpin.RW:
            print('Unable to set exposure time. Aborting...')
            return False

        exposure_time_to_set = exp #measured in us
        exposure_time_to_set = min(cam.ExposureTime.GetMax(), exposure_time_to_set)
        cam.ExposureTime.SetValue(exposure_time_to_set)
        print('Shutter time set to %s us...\n' % exposure_time_to_set)


        # setting gain
        if cam.GainAuto.GetAccessMode() != PySpin.RW:
            print('Unable to disable automatic Gain. Aborting...')
            return False

        cam.GainAuto.SetValue(PySpin.GainAuto_Off)

        if cam.Gain.GetAccessMode() != PySpin.RW:
            print('Unable to set gain. Aborting...')
            return False

        gain_to_set = gain #measured in dB
        gain_to_set = min(cam.Gain.GetMax(), gain_to_set)
        cam.Gain.SetValue(gain_to_set)
        print('Gain set to %s dB...\n' % gain_to_set)

        #set ADC to 12 bit
        node_ADC = PySpin.CEnumerationPtr(nodemap.GetNode('AdcBitDepth'))
        if not PySpin.IsAvailable(node_ADC) or not PySpin.IsWritable(node_ADC):
            print('Unable to get ADC (node retrieval). Aborting...')
            return False

        node_ADC_setting = node_ADC.GetEntryByName('Bit12')
        if not PySpin.IsAvailable(node_ADC_setting) or not PySpin.IsReadable(
                node_ADC_setting):
            print('Unable to set ADC (enum entry retrieval). Aborting...')
            return False
        node_ADC.SetIntValue(node_ADC_setting.GetValue())


        #setting the pixel format
        node_pixel_format = PySpin.CEnumerationPtr(nodemap.GetNode('PixelFormat'))
        if PySpin.IsAvailable(node_pixel_format) and PySpin.IsWritable(node_pixel_format):

            # Retrieve the desired entry node from the enumeration node
            node_pixel_format_mono8 = PySpin.CEnumEntryPtr(node_pixel_format.GetEntryByName('Mono12p'))
            if PySpin.IsAvailable(node_pixel_format_mono8) and PySpin.IsReadable(node_pixel_format_mono8):

                # Retrieve the integer value from the entry node
                pixel_format_mono8 = node_pixel_format_mono8.GetValue()

                # Set integer as new value for enumeration node
                node_pixel_format.SetIntValue(pixel_format_mono8)

                print('Pixel format set to %s...' % node_pixel_format.GetCurrentEntry().GetSymbolic())

                #check bits per pixel
                ##if PySpin.IsAvailable(node_pixel_size) and PySpin.IsReadable(node_pixel_size):
                    #print(node_pixel_size)

            else:
                print('Pixel format mono 12 not available...')

        else:
            print('Pixel format not available...')


    except PySpin.SpinnakerException as ex:
        print('Error: %s' % ex)
        result = False

    return result


class TriggerType:
    HARDWARE = 2

CHOSEN_TRIGGER = TriggerType.HARDWARE

def configure_trigger(cam):
    """
    This function configures the camera to use a trigger. First, trigger mode is
    set to off in order to select the trigger source. Once the trigger source
    has been selected, trigger mode is then enabled, which has the camera
    capture only a single image upon the execution of the chosen trigger.

     :param cam: Camera to configure trigger for.
     :type cam: CameraPtr
     :return: True if successful, False otherwise.
     :rtype: bool
    """
    result = True

    print('*** CONFIGURING TRIGGER ***\n')



    try:
        # Ensure trigger mode off
        # The trigger must be disabled in order to configure whether the source
        # is software or hardware.
        nodemap = cam.GetNodeMap()
        node_trigger_mode = PySpin.CEnumerationPtr(nodemap.GetNode('TriggerMode'))
        if not PySpin.IsAvailable(node_trigger_mode) or not PySpin.IsReadable(node_trigger_mode):
            print('Unable to disable trigger mode (node retrieval). Aborting...')
            return False

        node_trigger_mode_off = node_trigger_mode.GetEntryByName('Off')
        if not PySpin.IsAvailable(node_trigger_mode_off) or not PySpin.IsReadable(node_trigger_mode_off):
            print('Unable to disable trigger mode (enum entry retrieval). Aborting...')
            return False

        node_trigger_mode.SetIntValue(node_trigger_mode_off.GetValue())

        print('Trigger mode disabled...')

        # Select trigger source
        # The trigger source must be set to hardware or software while trigger
        # mode is off.
        node_trigger_source = PySpin.CEnumerationPtr(nodemap.GetNode('TriggerSource'))
        if not PySpin.IsAvailable(node_trigger_source) or not PySpin.IsWritable(node_trigger_source):
            print('Unable to get trigger source (node retrieval). Aborting...')
            return False

        if CHOSEN_TRIGGER == TriggerType.HARDWARE:
            node_trigger_source_hardware = node_trigger_source.GetEntryByName('Line3') #are we at 40Hz???
			#THERE IS NO METHOD FOR TRIGGER ACTIVATION, NOT SURE WHAT IS ACTUALLY CAUSING THE TRIGGER...
			#rising edge, level high???
            if not PySpin.IsAvailable(node_trigger_source_hardware) or not PySpin.IsReadable(
                    node_trigger_source_hardware):
                print('Unable to set trigger source (enum entry retrieval). Aborting...')
                return False
            node_trigger_source.SetIntValue(node_trigger_source_hardware.GetValue())

        # Turn trigger mode on
        # Once the appropriate trigger source has been set, turn trigger mode
        # on in order to retrieve images using the trigger.
        node_trigger_mode_on = node_trigger_mode.GetEntryByName('On')
        if not PySpin.IsAvailable(node_trigger_mode_on) or not PySpin.IsReadable(node_trigger_mode_on):
            print('Unable to enable trigger mode (enum entry retrieval). Aborting...')
            return False

        node_trigger_mode.SetIntValue(node_trigger_mode_on.GetValue())
        print('Trigger mode turned back on...')


    except PySpin.SpinnakerException as ex:
        print('Error: %s' % ex)
        return False

    return result



def grab_next_image_by_trigger(nodemap, cam):  # not sure what this does, but its working so it will be left
    """
		This function acquires an image by executing the trigger node.

		:param cam: Camera to acquire images from.
		:param nodemap: Device nodemap.
		:type cam: CameraPtr
		:type nodemap: INodeMap
		:return: True if successful, False otherwise.
		:rtype: bool
    """

    try:
        result = True
			# Use trigger to capture image
			# The software trigger only feigns being executed by the Enter key;
			# what might not be immediately apparent is that there is not a
			# continuous stream of images being captured; in other examples that
			# acquire images, the camera captures a continuous stream of images.
			# When an image is retrieved, it is plucked from the stream.


        if CHOSEN_TRIGGER == TriggerType.HARDWARE:
			 	# don't need to see this every image
				## print('Use the hardware to trigger image acquisition.')
            a = 2 #will this make it happy? (CP)

    except PySpin.SpinnakerException as ex:
        print('Error: %s' % ex)
        return False

    return result


def acquire_images(cam, nodemap, nodemap_tldevice,NUM_IMAGES,data_path,Folder_Name):
    """
		This function acquires and saves 10 images from a device.
		Please see Acquisition example for more in-depth comments on acquiring images.

		:param cam: Camera to acquire images from.
		:param nodemap: Device nodemap.
		:param nodemap_tldevice: Transport layer device nodemap.
		:type cam: CameraPtr
		:type nodemap: INodeMap
		:type nodemap_tldevice: INodeMap
		:return: True if successful, False otherwise.
		:rtype: bool
    """

    print('*** IMAGE ACQUISITION ***\n')
    try:
        result = True

        # Set acquisition mode to continuous
        # In order to access the node entries, they have to be casted to a pointer type (CEnumerationPtr here)
        node_acquisition_mode = PySpin.CEnumerationPtr(nodemap.GetNode('AcquisitionMode'))
        if not PySpin.IsAvailable(node_acquisition_mode) or not PySpin.IsWritable(node_acquisition_mode):
            print('Unable to set acquisition mode to continuous (enum retrieval). Aborting...')
            return False

        # Retrieve entry node from enumeration node
        node_acquisition_mode_continuous = node_acquisition_mode.GetEntryByName('Continuous')
        if not PySpin.IsAvailable(node_acquisition_mode_continuous) or not PySpin.IsReadable(
                node_acquisition_mode_continuous):
            print('Unable to set acquisition mode to continuous (entry retrieval). Aborting...')
            return False

        # Retrieve integer value from entry node
        acquisition_mode_continuous = node_acquisition_mode_continuous.GetValue()

        # Set integer value from entry node as new value of enumeration node
        node_acquisition_mode.SetIntValue(acquisition_mode_continuous)

        print('Acquisition mode set to continuous...')

        #  Begin acquiring images
        cam.BeginAcquisition()

        print('Acquiring images...')


        #  Retrieve device serial number for filename
        #
        #  *** NOTES ***
        #  The device serial number is retrieved in order to keep cameras from
        #  overwriting one another. Grabbing image IDs could also accomplish
        #  this.
        device_serial_number = ''
        node_device_serial_number = PySpin.CStringPtr(nodemap_tldevice.GetNode('DeviceSerialNumber'))
        if PySpin.IsAvailable(node_device_serial_number) and PySpin.IsReadable(node_device_serial_number):
            device_serial_number = node_device_serial_number.GetValue()
            print('Device serial number retrieved as %s...' % device_serial_number)


        print('')
        picList = []
        timelist = []

        # Retrieve, convert, and save images
        for i in range(NUM_IMAGES):
            try:

                #  Retrieve the next image from the trigger
                result &= grab_next_image_by_trigger(nodemap, cam)

                #  Retrieve next received image
                image_result = cam.GetNextImage()

                #  Ensure image completion
                if image_result.IsIncomplete():
                    print('Image incomplete with image status %d ...' % image_result.GetImageStatus())

                else:

                    #  Print image information; height and width recorded in pixels
                    #
                    #  *** NOTES ***
                    #  Images have quite a bit of available metadata including
                    #  things such as CRC, image status, and offset values, to
                    #  name a few.
					#
					#Don't need this right now (CP)
                    ##width = image_result.GetWidth()
                    ##height = image_result.GetHeight()
                    ##print('Grabbed Image %d, width = %d, height = %d' % (i, width, height))

                    #  Convert image to mono 8
                    #
                    #  *** NOTES ***
                    #  Images can be converted between pixel formats by using
                    #  the appropriate enumeration value. Unlike the original
                    #  image, the converted one does not need to be released as
                    #  it does not affect the camera buffer.
                    #
                    #  When converting images, color processing algorithm is an
                    #  optional parameter.
                    image_converted = image_result.Convert(PySpin.PixelFormat_Mono16, PySpin.HQ_LINEAR)

					#add to piclist
                    imgarray = image_converted.GetNDArray()
                    #picList.append(imgarray) #save piclist



					#saving imgarray
                    name = Folder_Name + '_%d' % i
                    numpy.save('{}/{}'.format(data_path, name), imgarray)

#Removed JPG saving process from original save_acquire on pjflab as of 12/16/2019 (kjl)

            except PySpin.SpinnakerException as ex:
                print('Error: %s' % ex)
                return False




        # End acquisition
        #
        #  *** NOTES ***
        #  Ending acquisition appropriately helps ensure that devices clean up
        #  properly and do not need to be power-cycled to maintain integrity.
        cam.EndAcquisition()



    except PySpin.SpinnakerException as ex:
        print('Error: %s' % ex)
        return False

    return result


def reset_trigger(nodemap):
    """
    This function returns the camera to a normal state by turning off trigger mode.

    :param nodemap: Transport layer device nodemap.
    :type nodemap: INodeMap
    :returns: True if successful, False otherwise.
    :rtype: bool
    """
    try:
        result = True
        node_trigger_mode = PySpin.CEnumerationPtr(nodemap.GetNode('TriggerMode'))
        if not PySpin.IsAvailable(node_trigger_mode) or not PySpin.IsReadable(node_trigger_mode):
            print('Unable to disable trigger mode (node retrieval). Aborting...')
            return False

        node_trigger_mode_off = node_trigger_mode.GetEntryByName('Off')
        if not PySpin.IsAvailable(node_trigger_mode_off) or not PySpin.IsReadable(node_trigger_mode_off):
            print('Unable to disable trigger mode (enum entry retrieval). Aborting...')
            return False

        node_trigger_mode.SetIntValue(node_trigger_mode_off.GetValue())

        print('Trigger mode disabled...')

    except PySpin.SpinnakerException as ex:
        print('Error: %s' % ex)
        result = False

    return result


def print_device_info(nodemap):
    """
    This function prints the device information of the camera from the transport
    layer; please see NodeMapInfo example for more in-depth comments on printing
    device information from the nodemap.

    :param nodemap: Transport layer device nodemap.
    :type nodemap: INodeMap
    :returns: True if successful, False otherwise.
    :rtype: bool
    """

    print('*** DEVICE INFORMATION ***\n')

    try:
        result = True
        node_device_information = PySpin.CCategoryPtr(nodemap.GetNode('DeviceInformation'))

        if PySpin.IsAvailable(node_device_information) and PySpin.IsReadable(node_device_information):
            features = node_device_information.GetFeatures()
            for feature in features:
                node_feature = PySpin.CValuePtr(feature)
                print('%s: %s' % (node_feature.GetName(),
                                  node_feature.ToString() if PySpin.IsReadable(node_feature) else 'Node not readable'))

        else:
            print('Device control information not available.')

    except PySpin.SpinnakerException as ex:
        print('Error: %s' % ex)
        return False

    return result


def run_single_camera(cam):
    """
    This function acts as the body of the example; please see NodeMapInfo example
    for more in-depth comments on setting up cameras.

    :param cam: Camera to run on.
    :type cam: CameraPtr
    :return: True if successful, False otherwise.
    :rtype: bool
    """
    try:
        result = True
        err = False

        # Retrieve TL device nodemap and print device information
        nodemap_tldevice = cam.GetTLDeviceNodeMap()

        result &= print_device_info(nodemap_tldevice)

        # Initialize camera
        cam.Init()

        # Retrieve GenICam nodemap
        nodemap = cam.GetNodeMap()

        #configure settings
        if configure_exp_gain(cam,exp,gain) is False:
            return False

        # Configure trigger
        if configure_trigger(cam) is False:
            return False

        # Acquire images
        result &= acquire_images(cam, nodemap, nodemap_tldevice,NUM_IMAGES,data_path,Folder_Name)

        # Reset trigger
        result &= reset_trigger(nodemap)

        # Deinitialize camera
        cam.DeInit()

    except PySpin.SpinnakerException as ex:
        print('Error: %s' % ex)
        result = False

    return result


def main():
    """
    Example entry point; please see Enumeration example for more in-depth
    comments on preparing and cleaning up the system.

    :return: True if successful, False otherwise.
    :rtype: bool
    """

    # Since this application saves images in the current folder
    # we must ensure that we have permission to write to this folder.
    # If we do not have permission, fail right away.
    try:
        test_file = open('test.txt', 'w+')
    except IOError:
        print('Unable to write to current directory. Please check permissions.')
        input('Press Enter to exit...')
        return False

    test_file.close()
    os.remove(test_file.name)

    result = True

    # Retrieve singleton reference to system object
    system = PySpin.System.GetInstance()

    # Get current library version
    version = system.GetLibraryVersion()
    print('Library version: %d.%d.%d.%d' % (version.major, version.minor, version.type, version.build))

    # Retrieve list of cameras from the system
    cam_list = system.GetCameras()

    num_cameras = cam_list.GetSize()

    print('Number of cameras detected: %d' % num_cameras)

    # Finish if there are no cameras
    if num_cameras == 0:
        # Clear camera list before releasing system
        cam_list.Clear()

        # Release system instance
        system.ReleaseInstance()

        print('Not enough cameras!')
        input('Done! Press Enter to exit...')
        return False

    # Run example on each camera
    for i, cam in enumerate(cam_list):

        print('Running example for camera %d...' % i)

        result &= run_single_camera(cam)
        print('Camera %d example complete... \n' % i)

    # Release reference to camera
    # NOTE: Unlike the C++ examples, we cannot rely on pointer objects being automatically
    # cleaned up when going out of scope.
    # The usage of del is preferred to assigning the variable to None.
    del cam

    # Clear camera list before releasing system
    cam_list.Clear()

    # Release system instance
    system.ReleaseInstance()

    input('Done! Press Enter to exit...')
    return result


if __name__ == '__main__':
    main()
