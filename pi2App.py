#pi2App 01 February, 2021
#Author: Matthias Schmid, LHEP, University of Bern
#Released under MIT licence

from tkinter import *
from tkinter import filedialog
from PIL import Image, ImageTk, ImageStat
import io
import picamera
import numpy
import time
import os
import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
from matplotlib.figure import Figure
from matplotlib import cm
import matplotlib.pyplot as plt
import cv2
import numpy as np
from scipy.optimize import curve_fit

class CaptureApp():
    
    def getIntensityCrosshair( self, event ):
        #lists that will contain intensity values
        xIntensity = []
        yIntensity = []
        #lists to plot x, y axis
        xPoints = numpy.linspace( 0, self.img_width-1, num=self.img_width )
        yPoints = numpy.linspace( 0, self.img_height-1, num=self.img_height )
        #event coordinates
        ix, iy = event.x, event.y
        px = self.imgCorr.load()
        
        #store x values
        for x in range( 0, self.img_width ):
            xIntensity.append( px[x, iy] )
        #store y values
        for y in range( 0, self.img_height ):
            yIntensity.append( px[ix, y] )

        #take unedited picture and load it again
        self.jmg = self.imgCorr.copy()
        px = self.jmg.load()
        #draw a crosshair
        for x in range( 0, self.img_width ):
            px[ x, iy ] = 255
        for y in range( 0, self.img_height ):
            px[ ix, y ] = 255

        #save picture with drawn crosshair in j
        #self.jmg.save( self.fpath+"EDIT"+".png" )
        #update canvas with crosshair picture
        self.jmgTk = ImageTk.PhotoImage(self.jmg)
        self.can.create_image( self.img_width/2, self.img_height/2, image = self.jmgTk )

        #clear existing axes
        self.ax1.clear()
        self.ax2.clear()
        
        #update axes of default plot
        #subplot for x values
        self.ax1.plot( xPoints, xIntensity )
        self.ax1.set_ylabel( "Intensity" )
        self.ax1.set_title( "X-direction", y = 1.08, fontsize=10 )
        self.ax1.set_xlabel( "Pixel No." )
        self.ax1.set_xlim( 0, self.img_width-1 )
        self.ax1.set_autoscalex_on( False )
        self.ax1b.set_autoscalex_on( False )
        #to scale the cm axis: use scaled limits
        self.ax1b.set_xlim( 0, (self.img_width-1)*self.pxFactorx )
        self.ax1.set_autoscaley_on( False )
        self.ax1.set_ylim( [0,255] )
        #subplot for y values
        self.ax2.plot( yPoints, yIntensity )
        self.ax2.set_title( "Y-direction", y = 1.08, fontsize=10 )
        self.ax2.set_xlabel( "Pixel No." )
        self.ax2.set_xlim( 0, self.img_height - 1 )
        self.ax2.set_autoscalex_on( False )
        self.ax2b.set_autoscalex_on( False )
        self.ax2b.set_xlim( 0, (self.img_height - 1 )*self.pxFactory )
        self.ax2.set_autoscaley_on( False )
        self.ax2.set_ylim( [0,255] )

        #update canvas
        self.plot_canvas.draw()

    def getIntensityPixel( self, event ):
        px=self.imgCorr.load() #use corrected image
        print( "Intensity at (%d, %d): %d" %( event.x, event.y, px[event.x, event.y] ) )

    def getTotalIntensity( self, image):
        #takes an image and returns its total intensity
        stat = ImageStat.Stat( image )
        return stat.sum[0]
    
    def getReducedImage( self, image ):
        img = image
        px = img.load()

        reducedImg = Image.new("L", (640, 480))
        redpx = reducedImg.load()
        for x in range(640):
            for y in range(480):
                redpx[x,y] = 0

        #Set all entries on foil to 1
        p1 = (105, 28)
        p2 = (12, 137)
        p3 = (5, 219)
        p4 = (20, 322)
        p5 = (118, 438)
        p6 = (280, 471)
        p7 = (450, 438)
        p8 = (567, 322)
        p9 = (582, 219)
        p10 = (567, 137)
        p11 = (487, 28)


        for y in range(0, p2[1]-p1[1]):#y from 0 to dist(p1[1], p2[1])
            for x in range(int(p1[0] + y*(p1[0]-p2[0])/(p1[1]-p2[1])), int(p11[0] + y*(p10[0]-p11[0])/(p10[1]-p11[1])) ):
                redpx[x,p1[1]+y] = px[x,p1[1]+y]
        
        for y in range(0, p3[1]-p2[1]):
            for x in range(int(p2[0] + y*(p2[0]-p3[0])/(p2[1]-p3[1])), int(p10[0] + y*(p9[0]-p10[0])/(p9[1]-p10[1])) ):
                redpx[x,p2[1]+y] = px[x,p2[1]+y]

        for y in range(0, p4[1]-p3[1]):
            for x in range(int(p3[0] + y*(p3[0]-p4[0])/(p3[1]-p4[1])), int(p9[0] + y*(p8[0]-p9[0])/(p8[1]-p9[1])) ):
                redpx[x,p3[1]+y] = px[x,p3[1]+y]
        
        for y in range(0, p5[1]-p4[1]):
            for x in range(int(p4[0] + y*(p4[0]-p5[0])/(p4[1]-p5[1])), int(p8[0] + y*(p7[0]-p8[0])/(p7[1]-p8[1])) ):
                redpx[x,p4[1]+y] = px[x,p4[1]+y]
        
        for y in range(0, p6[1]-p5[1]):
            for x in range(int(p5[0] + y*(p5[0]-p6[0])/(p5[1]-p6[1])), int(p7[0] + y*(p6[0]-p7[0])/(p6[1]-p7[1])) ):
                redpx[x,p5[1]+y] = px[x,p5[1]+y]

        return reducedImg
    
    def gauss(self, x, a, x0, sigma):    
        return a*np.exp( -(x-x0)**2/(2*sigma**2) )
    
    def getImageStats( self ):
        img = self.reducedImg
        px = img.load()
        xPoints = np.linspace(0, 639, num=640)
        yPoints = np.linspace(0, 479, num=480)
        xIntensity = []#intensity profile along x-axis
        yIntensity = []
        xSum = 0
        ySum = 0
        for x in range(640):
            for y in range(480):
                xSum += px[x,y]
            xIntensity.append(xSum)
            xSum = 0
            
        for y in range(480):
            for x in range(640):
                ySum += px[x,y]
            yIntensity.append(ySum)
            ySum = 0
            
        #Get weighted mean
        nx = sum(xIntensity)
        ny = sum(yIntensity)

        meanx = sum(xPoints * xIntensity )/ nx
        stdx = sum(xIntensity*(xPoints-meanx)**2)/ nx
        ax = max(xIntensity)
        meany = sum(yPoints*yIntensity) / ny
        stdy = sum(yIntensity*(yPoints-meany)**2)/ ny
        ay = max(yIntensity)        

        #Try fitting a Gaussian
        #Use statistical values above as starting parameters
        #popt values: a, mu, sigma
        try:
            poptx, pcovx = curve_fit( self.gauss, xPoints, xIntensity, p0=[ax, meanx, stdx])
            popty, pcovy = curve_fit( self.gauss, yPoints, yIntensity, p0=[ay, meany, stdy])

            #exclude fit if sigma x > 100 mm = 1155 px
            if abs(poptx[2]) > 1155:
                poptx = [0,0,0]
                print("Sigma x too large.")
            #exclude fit if sigma y > 100 mm = 734 px
            if abs(popty[2]) > 734:
                popty = [0,0,0]
                print("Sigma y too large.")
            self.meanx_txtvar.set("Mean x: %.1f mm" %( poptx[1]*self.pxFactorx )) #convert to mm
            self.meany_txtvar.set("Mean y: %.1f mm" %( popty[1]*self.pxFactory ))
            print("X: mu = %.1f, sigma = %.1f, a=%.1f" %(poptx[1], abs(poptx[2]), poptx[0]))
            print("Y: mu = %.1f, sigma = %.1f, a=%.1f" %(popty[1], abs(popty[2]), popty[0]))

            fwhmx = 2.3548*abs(poptx[2])*self.pxFactorx #convert to mm
            fwhmy = 2.3548*abs(popty[2])*self.pxFactory
            self.fwhmx_txtvar.set("FWHM x: %.1f mm" %( fwhmx ))
            self.fwhmy_txtvar.set("FWHM y: %.1f mm" %( fwhmy ))

            #clear existing axes
            self.ax1.clear()
            self.ax2.clear()
        
            #update axes of default plot
            #subplot for x values
            self.ax1.plot( xPoints, xIntensity )
            if poptx[ 2 ] != 0: #to prevent division by 0
                self.ax1.plot(xPoints, self.gauss(xPoints, *poptx), "r")
            self.ax1.set_ylabel( "Intensity" )
            self.ax1.set_title( "X-direction", y = 1.08, fontsize=10 )
            self.ax1.set_xlabel( "Pixel No." )
            self.ax1.set_xlim( 0, self.img_width-1 )
            self.ax1.set_autoscalex_on( False )
            self.ax1b.set_autoscalex_on( False )
            #self.ax1.ticklabel_format( self, style="sci" )
            #to scale the mm axis: use scaled limits
            self.ax1b.set_xlim( 0, (self.img_width-1)*self.pxFactorx )
            self.ax1.set_autoscaley_on( False )
            #subplot for y values
            self.ax2.plot( yPoints, yIntensity )
            if popty[ 2 ] != 0: #to prevent division by 0
                self.ax2.plot(yPoints, self.gauss(yPoints, *popty), "r")
            self.ax2.set_title( "Y-direction", y = 1.08, fontsize=10 )
            self.ax2.set_xlabel( "Pixel No." )
            self.ax2.set_xlim( 0, self.img_height - 1 )
            self.ax2.set_autoscalex_on( False )
            self.ax2b.set_autoscalex_on( False )
            self.ax2b.set_xlim( 0, (self.img_height - 1 )*self.pxFactory )
            self.ax2.set_autoscaley_on( False )

            self.plot_canvas.draw()
        except RuntimeError:
            print("No optimal fit parameters found.")

        
    def setExpTime( self ):
        self.expTime = self.expTime_entry.get()
        self.expTime_txtvar.set("Exposure time: " + str( self.expTime ) + " ms")

    def captureImage( self ):
        #Stop live view if running
        self.stopLive()
        camera = picamera.PiCamera()
        camera.resolution = ( 640, 480 )

        expTime_seconds = int( self.expTime ) * 10**( -3 )

        #Proceed with automatic or manual exposure time
        if self.expIsAuto.get() == True:
            camera.exposure_mode = "auto"
        else:
            try:
                frate = 1.0 / float( expTime_seconds )
                camera.framerate = frate
            except ValueError:
                print("Invalid frame rate (exposure time too low)!")
            except ZeroDivisionError:
                print("Invalid exposure time (0)!")
                
        time.sleep(1) #let gain values settle
        camera.awb_mode = "off" #disable automatic white balance
        camera.awb_gains = (2.0 , 2.0) #fix awb gains at arbitrary value

        #save to stream, file has to be manually saved by user
        imgStream = io.BytesIO()
        ##flip image vertically, so beam axes correspond to image axes
        #camera.vflip=True
        camera.capture( imgStream, format="png" )

        #update ET label with actual ET
        self.expTime_txtvar.set("Exposure Time: " + str( camera.exposure_speed/1000 ) + " ms" )

        #Get gain values to correct intensity
        analog_gain = camera.analog_gain
        digital_gain = camera.digital_gain
        total_gain = analog_gain*digital_gain
        awb = camera.awb_gains
        print( "Capturing (ET = %d millisecs, FPS: %d, gain: %d)..." %( float( camera.exposure_speed/1000 ),
              int( camera.framerate ) , float(total_gain)))#actual expTime and fps from picamera
        #update image to current
        self.img = Image.open( imgStream ).convert( mode = 'L')
        #Get normalized intensity
        #Important: not to be calculated from perspective corrected image
        totInt = self.getTotalIntensity( self.img )
        normInt = float( totInt ) / ( float( camera.exposure_speed ) * float(digital_gain) * float(analog_gain) )
        self.norm_int_txtvar = "Norm. Intensity: " + "{:.2e}".format( normInt )
        self.norm_int_label.configure( text = self.norm_int_txtvar )
        #update total intensity
        self.tot_int_txtvar = "Tot. Intensity: " + "{:.2e}".format( totInt )
        self.tot_int_label.configure( text = self.tot_int_txtvar )
        #update beam current
        #self.current_txtvar = "Beam Current: " + str(normInt*30.3-90) + " nA"
        
        current = normInt*30.3-90 #in nA
        currentError = np.sqrt((0.9*current) + 900)
        self.current_txtvar = u"Beam Current: %.0f +- %.0f nA" %(current, currentError)
        self.current_label.configure( text = self.current_txtvar )

        #perspective correction
        self.imgCv = np.array( self.img )
        self.imgCvCorr = cv2.warpPerspective( self.imgCv, self.corrMat, (640,  480) )
        self.imgCorr = Image.fromarray( self.imgCvCorr )
        self.imgTk = ImageTk.PhotoImage( self.imgCorr )        
        #coordinates are center
        self.can.create_image( self.img_width/2, self.img_height/2, image = self.imgTk )

        #Get reduced image form perspective corrected image
        self.reducedImg = self.getReducedImage( self.imgCorr )

        self.getImageStats()

        #close camera before next initiation
        camera.close()

    def saveImage( self ):
        fileDest = filedialog.asksaveasfilename( initialdir="/home/pi/Desktop", title="Select file" )
        if fileDest == "":#if user cancels, empty string
            pass
        elif fileDest == ():#cancel after selecting returns empty tuple!
            pass
        else:
            self.img.save( fileDest+".png", format="png" )
        

    def stopLive( self ):
        self.isLive.set(False)

        
    def getTime( self ):
        return time.strftime( "%H-%M-%S" )
        
    def expTimeSetMan( self ):
        self.expIsAuto.set(False)
        self.expTime_txtvar.set("Exposure time: -")
        self.expTime_entry.config(state=NORMAL)
        self.expTime_button.config(state=NORMAL)
        
    def expTimeSetAuto( self ):
        self.expIsAuto.set(True)
        self.expTime = 0
        self.expTime_txtvar.set("Exposure time: Auto")
        self.expTime_entry.config(state=DISABLED)
        self.expTime_button.config(state=DISABLED)

    def startLive( self ):
        self.isLive.set(True)
        print("Starting live view.")
        #Update the live view button
        self.live_button.config( text = "Stop Live View", command = lambda:self.stopLive() )
        #master.update()

        #go to live view loop
        self.waitLive()
        
    def waitLive( self ):
        #call getLiveImage if isLive flag is True
        if self.isLive.get():
            self.getLiveImage()
        else:
            print("Stopping live view.")
            #update live button
            self.live_button.config( text = "Start Live View", command = lambda:self.startLive() )
            master.update()
            return
        #waits so isLive value can be updated
        master.after(10, self.waitLive())
    
    def testFunc( self ):
        print("Test.")
        
    def getLiveImage( self ):
        stream = io.BytesIO()
        with picamera.PiCamera() as camera:
            camera.resolution = (640, 480)
            camera.capture(stream, format="jpeg")
        stream.seek(0)
        self.vidImage = Image.open(stream)
        self.vidImageTk = ImageTk.PhotoImage(self.vidImage)
        self.can.delete( "all" )
        self.can.create_image( 0, 0, anchor=NW, image = self.vidImageTk )
        master.update() 
    
    def closeApp( self ):
        master.quit()

    def __init__( self, main ):

        #Initiate image variables with black image
        self.img = Image.new("L", (640,480))
        self.imgCorr = self.img
        self.jmg = self.img

        #PhotoImage for display
        self.imgTk = ImageTk.PhotoImage( self.img )
        self.jmgTk = ImageTk.PhotoImage( self.jmg )
        
        #conversion factor from pixel to mm, has to be calibrated
        self.pxFactorx = 0.0866 #1 px = 0.0866 mm
        self.pxFactory = 0.136 #1 px = 0.136 mm

        #Matrix for perspective correction (obtained on 23 March 2020)
        self.corrMat = np.array([   [9.19601973e-01, -4.19156526e-01, 2.81428327e+01],
                                    [7.99019461e-03,  5.18085721e-01, 2.65480402e+01],
                                    [3.11816220e-05, -1.34525758e-03, 1.0           ]])
        
        #get image dimensions
        self.img_width = self.imgTk.width()
        self.img_height = self.imgTk.height()
        
        #Exposure time variables
        self.expTime = 0 #auto by default. User enters milliseconds, camera needs micro.
        self.expTime_txtvar = StringVar()
        self.expTime_txtvar.set("Exposure time: Auto")

        #----------
        #TKINTER
        #----------
        
        self.isLive = BooleanVar()
        self.isLive.set(True)
        self.expIsAuto = BooleanVar()
        self.expIsAuto.set(True)
        
        #build master window, place in upper left corner
        master.geometry( "+0+0" )
        master.title("pi2App")
        
        #Four frames contain the elements to display
        self.imageFrame = Frame( master, borderwidth=1, relief=RAISED )
        self.imageFrame.grid( row = 0, column = 0, sticky=W+S+E+N )
        self.controlFrame = Frame( master, borderwidth=1, relief=RAISED, width=640 )
        self.controlFrame.grid( row = 1, column = 0, sticky=W+S+E+N )
        self.graphFrame = Frame( master, borderwidth=1, relief=RAISED, height=480, width=640 )
        self.graphFrame.grid( row = 0, column = 1, sticky=W+S+E+N )
        self.analysisFrame = Frame( master, borderwidth=1, relief=RAISED, width=640 )
        self.analysisFrame.grid( row = 1, column = 1, sticky=W+S+E+N )
        
        #experimental toolbar frame
        self.toolbarFrame = Frame( master, borderwidth=1 )
        self.toolbarFrame.grid( row = 1, column = 2 )

        #IMAGE FRAME
        #canvas for image display
        self.can = Canvas( self.imageFrame, width = self.img_width, height = self.img_height )
        self.can.create_image( self.img_width/2, self.img_height/2, image = self.imgTk )
        self.can.grid( row = 0, column = 0 )
        self.can.bind( "<Double-Button-1>", self.getIntensityCrosshair )
        self.can.bind( "<Button-1>", self.getIntensityPixel )
        #filepath to store pictures
        self.fpath = "Images/"#+timestr+".png" will be added!
        
        #CONTROL FRAME
        #Entry, Label and button for shutter time
        self.expTime_entry = Entry( self.controlFrame, width = 10, state=DISABLED )
        self.expTime_entry.grid( row = 2, column = 3, sticky=W+E )
        self.expTime_button = Button( self.controlFrame, text = "Set (ms)", command = lambda:self.setExpTime(), state=DISABLED )
        self.expTime_button.grid( row = 2, column = 4, sticky=W )

        #button for new capture
        self.capture_button = Button( self.controlFrame, text = "Capture Image", command = lambda:self.captureImage() )
        self.capture_button.grid( row = 0, column = 1, columnspan = 2, sticky=W )
        #button to store the taken image
        self.save_button = Button( self.controlFrame, text = "Save Image", command=lambda:self.saveImage() )
        self.save_button.grid( row = 0, column = 3, sticky=W )
        #button to start live view
        self.live_button = Button( self.controlFrame, text = "Start Live View", command = lambda:self.startLive() )
        self.live_button.grid( row = 0, column = 0, rowspan=3, sticky=N+E+W+S )
        #button to show fit
        self.fit_button = Button( self.controlFrame, text = "Show Fit", command = lambda:self.getImageStats() )
        self.fit_button.grid( row = 0, column = 4, sticky=W )
        #radiobutton: exposure time manual/auto
        self.expRadioButton1 = Radiobutton( self.controlFrame, text="Auto", variable=self.expIsAuto, value=True, command=lambda:self.expTimeSetAuto() )
        self.expRadioButton1.grid( row = 1, column = 1, sticky=W)
        self.expRadioButton2 = Radiobutton( self.controlFrame, text="Manual:", variable=self.expIsAuto, value=False, command=lambda:self.expTimeSetMan() )
        self.expRadioButton2.grid( row = 2, column = 1, sticky=W)        

        #label to display exposure time
        self.expTime_label = Label( self.controlFrame, textvariable=self.expTime_txtvar)
        self.expTime_label.grid( row = 0, column = 5 )
        
        #ANALYSIS FRAME
        #label for total intensity
        self.tot_int_txtvar = "Tot. Intensity: -"
        self.tot_int_label = Label( self.analysisFrame, text = self.tot_int_txtvar )
        self.tot_int_label.grid( row = 0, column = 0 )

        #Label for averaged intensity
        self.norm_int_txtvar = "Norm. Intensity: -"
        self.norm_int_label = Label( self.analysisFrame, text = self.norm_int_txtvar )
        self.norm_int_label.grid( row = 1, column = 0 )

        #Label for beam current
        self.current_txtvar = "Beam Current: -"
        self.current_label = Label( self.analysisFrame, text = self.current_txtvar )
        self.current_label.grid( row = 0, column = 1 )

        #Labels for beam size
        self.meanx_txtvar = StringVar()
        self.meany_txtvar = StringVar()
        self.fwhmx_txtvar = StringVar()
        self.fwhmy_txtvar = StringVar()
        self.meanx_txtvar.set("Mean x: - mm")
        self.meany_txtvar.set("Mean y: - mm")
        self.fwhmx_txtvar.set("FWHM x: - mm")
        self.fwhmy_txtvar.set("FWHM y: - mm")
        self.meanx_label = Label( self.analysisFrame, textvariable = self.meanx_txtvar )
        self.meanx_label.grid( row = 0, column = 2)
        self.meany_label = Label( self.analysisFrame, textvariable = self.meany_txtvar )
        self.meany_label.grid( row = 1, column = 2 )
        self.fwhmx_label = Label( self.analysisFrame, textvariable = self.fwhmx_txtvar )
        self.fwhmx_label.grid( row = 0, column = 3 )
        self.fwhmy_label = Label( self.analysisFrame, textvariable = self.fwhmy_txtvar )
        self.fwhmy_label.grid( row = 1, column = 3 )
 

        #close button
        self.close_button = Button( self.analysisFrame, text = "Close", command = lambda:self.closeApp() )
        self.close_button.grid( row = 2, column = 2, sticky = "SE" )

        #GRAPH FRAME
        #default graph
        #may change size according to screen resolution
        self.fig = Figure()
        self.ax1 = self.fig.add_subplot( 121 )
        self.ax1.set_title( "X-direction", y = 1.1, fontsize = 10 ) #y raises title
        self.ax1.set_ylabel( "Intensity" )
        self.ax1.set_xlabel( "Pixel No." )
        self.ax1b = self.ax1.twiny() #a twin of ax1, y invisible
        self.ax1b.set_xlabel( "mm" )
        self.ax1b.set_xlim( 0, ( self.img_width-1 )*self.pxFactorx )
        plt.plot( [1,2], [1,2] )
        
        self.ax2 = self.fig.add_subplot( 122 )
        self.ax2.set_title( "Y-direction", y = 1.1, fontsize=10 )
        self.ax2.set_xlabel( "Pixel No." )
        self.ax2b = self.ax2.twiny() #a twin of ax1, y invisible
        self.ax2b.set_xlabel( "mm" )
        self.ax1b.set_xlim( 0, ( self.img_width-1 )*self.pxFactory )
        plt.plot( [1,2],[1,2] )

        #canvas for graph display
        self.plot_canvas = FigureCanvasTkAgg( self.fig, self.graphFrame )
        self.plot_canvas.show()
        self.plot_canvas.get_tk_widget().pack( expand=TRUE )
        """
        self.toolbar = NavigationToolbar2TkAgg( self.fig, self.toolbarFrame )
        self.toolbar.update()
        self.plot_canvas._tkcanvas.pack()
        """

master = Tk()
CaptureApp(master)
master.mainloop()

