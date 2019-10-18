#  Trigger.py shows how to trigger the camera. It relies on information
#  provided in the Enumeration, Acquisition, and NodeMapInfo examples.
#
#  It can also be helpful to familiarize yourself with the ImageFormatControl
#  and Exposure examples. As they are somewhat shorter and simpler, either
#  provides a strong introduction to camera customization.
#
#  This example shows the process of configuring, using, and cleaning up a
#  camera for use with both a software and a hardware trigger.


# ADDING BITS SLOWLY TO ORIGINAL tRIGGER.PY SPINNAKER EXAMPLE
# DOUBLE POUND (##) IS COMMENTED OUT ORIGINAL CODE

import os
import PySpin
import numpy
import matplotlib.pyplot as plt
import time
from PIL import Image
from numpy import array, empty, ravel, where, ones, reshape, arctan2
from matplotlib.pyplot import plot, draw, show, ion

NUM_IMAGES = 9000  # number of images to grab
#eventually change this to an input value 


class TriggerType:
    #SOFTWARE = 1
    HARDWARE = 2


CHOSEN_TRIGGER = TriggerType.HARDWARE
#we only want HARDWARE so permanently set this later (CP), or get rid of software code bits

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

        ##if CHOSEN_TRIGGER == TriggerType.SOFTWARE:
        ##    node_trigger_source_software = node_trigger_source.GetEntryByName('Software')
        ##    if not PySpin.IsAvailable(node_trigger_source_software) or not PySpin.IsReadable(
        ##            node_trigger_source_software):
        ##        print('Unable to set trigger source (enum entry retrieval). Aborting...')
        ##        return False
        ##    node_trigger_source.SetIntValue(node_trigger_source_software.GetValue())

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


def grab_next_image_by_trigger(nodemap, cam):  #is this necessary if we don't have a software trigger???
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

        #if CHOSEN_TRIGGER == TriggerType.SOFTWARE:
            # Get user input
        #    input('Press the Enter key to initiate software trigger.')

            # Execute software trigger
        #    node_softwaretrigger_cmd = PySpin.CCommandPtr(nodemap.GetNode('TriggerSoftware'))
        #    if not PySpin.IsAvailable(node_softwaretrigger_cmd) or not PySpin.IsWritable(node_softwaretrigger_cmd):
        #        print('Unable to execute trigger. Aborting...')
        #        return False

        #    node_softwaretrigger_cmd.Execute()

            # TODO: Blackfly and Flea3 GEV cameras need 2 second delay after software trigger

        if CHOSEN_TRIGGER == TriggerType.HARDWARE:
			# don't need to see this every image
            ## print('Use the hardware to trigger image acquisition.')
            a = 2 #will this make it happy? (CP)
			
#I AM REaLLY SCARED ABOUT THIS I DO NOT KNOW WHAT IS ACTIVATING THE TRIGGER. IS IT REALLY THE CLOCK GENERATOR?
    except PySpin.SpinnakerException as ex:
        print('Error: %s' % ex)
        return False

    return result

def Novak_phase(listo):
    #phase of each pixel, assuming: list of five images read in, equally centered and sized
    arr1 = numpy.ravel(numpy.array(listo[0],dtype='int'))
    arr2 = numpy.ravel(numpy.array(listo[1],dtype='int')) #converts to numpy arrays for faster operation
    arr3 = numpy.ravel(numpy.array(listo[2],dtype='int'))
    arr4 = numpy.ravel(numpy.array(listo[3],dtype='int'))
    arr5 = numpy.ravel(numpy.array(listo[4],dtype='int'))
    phase = numpy.empty(388800)

    mask = numpy.ones(388800,dtype=bool)

    cuts = numpy.where(arr1 < 15)

    mask[cuts] = False

    p1 = arr1[mask]
    p2 = arr2[mask]
    p3 = arr3[mask]
    p4 = arr4[mask]
    p5 = arr5[mask]	

    den = 2*p3-p1-p5

    A = p2-p4

    B = p1-p5

    num = numpy.sqrt(abs(4*A**2-B**2))

    pm = numpy.sign(A)

    pha = numpy.arctan2(pm*num,den)
	
    phase[~mask] = 0
    phase[mask] = pha

    phase = numpy.reshape(phase,(540,720))

    return phase
	
def fourpointphase(listo):
    #phase of each pixel, assuming: list of four images read in, equally centered and sized
    arr1 = numpy.ravel(numpy.array(listo[0],dtype='int'))
    arr2 = numpy.ravel(numpy.array(listo[1],dtype='int')) #converts to numpy arrays for faster operation
    arr3 = numpy.ravel(numpy.array(listo[2],dtype='int'))
    arr4 = numpy.ravel(numpy.array(listo[3],dtype='int'))

    phase = numpy.empty(388800)

    mask = numpy.ones(388800,dtype=bool)

    cuts = numpy.where(arr1 < 15)

    mask[cuts] = False

    p1 = arr1[mask]
    p2 = arr2[mask]
    p3 = arr3[mask]
    p4 = arr4[mask]

    num = p4 - p2
    den = p1 - p3
    pha = numpy.arctan2(num,den)

    phase[~mask] = 0
    phase[mask] = pha

    phase = numpy.reshape(phase,(540,720))

    return phase	
	
def acquire_images(cam, nodemap, nodemap_tldevice):
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
			
			
		#From Novak (CP)
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
                    image_converted = image_result.Convert(PySpin.PixelFormat_Mono8, PySpin.HQ_LINEAR)
					
					#add to piclist
                    imgarray = image_converted.GetNDArray()
                    picList.append(imgarray)

                    # Create a unique filename
					#maybe later? (CP)
                    ##if device_serial_number:
                    ##    filename = 'Trigger-%s-%d.jpg' % (device_serial_number, i)
                    ##else:  # if serial number is empty
                    ##    filename = 'Trigger-%d.jpg' % i

                    # Save image (no thanks) 
                    #
                    #  *** NOTES ***
                    #  The standard practice of the examples is to use device
                    #  serial numbers to keep images of one device from
                    #  overwriting those of another.
                    ##image_converted.Save(filename)
                    ##print('Image saved at %s\n' % filename)
					
					#  Release image (from Novak, LivePhase_V3)
                    image_result.Release()
                    if i%4 == 0 and i > 0:
                        if i%8 == 0:
                            faze = fourpointphase(picList) #Novak #4 point
                        
                            #print(faze.dtype)
                            #print(numpy.shape(faze))                     
                            plt.ion()						
                            plt.imshow(faze, cmap = 'jet')
                            #cbar = plt.colorbar()
                            #plt.clim(vmin=-numpy.pi,vmax=numpy.pi)
                            #cbar.set_label("Phase Shift (rad)")
                            plt.xlabel("Pixels(x)")
                            plt.ylabel("Pixels(y)")
                            #plt.xlim([300,500])
                            #plt.ylim([200,400])
                            plt.pause(0.00001)
                            plt.show()
                            plt.clf()
                        
                        del picList[0:4] 
						
					

            except PySpin.SpinnakerException as ex:
                print('Error: %s' % ex)
                return False

		#some more stuff from V3	
        timelist = numpy.asarray(timelist)
        timediffs = timelist[1:NUM_IMAGES]-timelist[0:NUM_IMAGES-1]
        print("Average time between image captures: ", numpy.mean(timediffs)," seconds")
        print("Standard deviation: ", numpy.std(timediffs)," seconds")
        print(timediffs[0:100])
        plt.ion()
        plt.plot(timediffs)
        plt.show()
		
		
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

        # Configure trigger
        if configure_trigger(cam) is False:
            return False

        # Acquire images
        result &= acquire_images(cam, nodemap, nodemap_tldevice)

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
