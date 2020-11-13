package immunity;

import java.util.HashMap;
import java.util.Set;

import org.COPASI.CCompartment;
import org.COPASI.CCopasiDataModel;
import org.COPASI.CCopasiMessage;
import org.COPASI.CCopasiMethod;
import org.COPASI.CCopasiObjectName;
import org.COPASI.CCopasiParameter;
import org.COPASI.CCopasiReportSeparator;
import org.COPASI.CCopasiRootContainer;
import org.COPASI.CCopasiTask;
import org.COPASI.CMetab;
import org.COPASI.CModel;
import org.COPASI.CReaction;
import org.COPASI.CRegisteredObjectName;
import org.COPASI.CReportDefinition;
import org.COPASI.CReportDefinitionVector;
import org.COPASI.CTrajectoryMethod;
import org.COPASI.CTrajectoryProblem;
import org.COPASI.CTrajectoryTask;
import org.COPASI.ReportItemVector;

public class RabConversion {
	private static RabConversion instance = null;
	private CCopasiDataModel dataModel;
	private CModel model;
    private CReportDefinition report;
    private CTrajectoryTask trajectoryTask;
	private HashMap<String, CMetab> nameMetabs = new HashMap<String, CMetab>();
	
	public static RabConversion getInstance () {
		if (instance == null) {
			instance = new RabConversion();
		}
		
		return instance;
	}
	
	protected RabConversion(){
		
		// to defeat instantiation
		assert CCopasiRootContainer.getRoot() != null;
        // create a new datamodel
		dataModel = CCopasiRootContainer.addDatamodel();
        assert CCopasiRootContainer.getDatamodelList().size() == 1;
        
        String modelFileName = CellProperties.getInstance().getCopasiFiles().get("rabCopasi");//"rabs_conversion.cps";
        
        try
        {
            // load the model without progress report
            dataModel.loadModel(modelFileName);
        }
        catch (java.lang.Exception ex)
        {
            System.err.println("Error while loading the model from file named \"" + modelFileName + "\".");
            System.exit(1);
        }
        
        model = dataModel.getModel();
        assert model != null;
        
     // output number and names of all compartments
        int i, iMax = (int)model.getCompartments().size();
        for (i = 0;i < iMax;++i)
        {
            CCompartment compartment = model.getCompartment(i);
            assert compartment != null;
        }

        // output number and names of all metabolites
        iMax = (int)model.getMetabolites().size();
        for (i = 0;i < iMax;++i)
        {
            CMetab metab = model.getMetabolite(i);
            assert metab != null;
            nameMetabs.put(metab.getObjectName(), metab);
        }
        // SET INITIAL CONCENTRATIONS
        
        for (String s : nameMetabs.keySet()) {
        	CMetab metab = nameMetabs.get(s);
        }

        // output number and names of all reactions
        iMax = (int)model.getReactions().size();
        for (i = 0;i < iMax;++i)
        {
            CReaction reaction = model.getReaction(i);
            assert reaction != null;
        }
        
       setUpReport();
        setUpTask();
	}
	
	private void setUpReport() {
		// create a report with the correct filename and all the species against
        // time.
        CReportDefinitionVector reports = dataModel.getReportDefinitionList();
        // create a new report definition object
        report = reports.createReportDefinition("Report", "Output for timecourse");
        // set the task type for the report definition to timecourse
        report.setTaskType(CCopasiTask.timeCourse);
        // we don't want a table
        report.setIsTable(false);
        // the entries in the output should be seperated by a ", "
        report.setSeparator(new CCopasiReportSeparator(", "));
        
        // we need a handle to the header and the body
        // the header will display the ids of the metabolites and "time" for
        // the first column
        // the body will contain the actual timecourse data
        ReportItemVector header = report.getHeaderAddr();
        ReportItemVector body = report.getBodyAddr();
        
        body.add(new CRegisteredObjectName(model.getObject(new CCopasiObjectName("Reference=Time")).getCN().getString()));
 
     
	}
	
	private void setUpTask() {
		// get the trajectory task object
        trajectoryTask = (CTrajectoryTask)dataModel.getTask("Time-Course");
        assert trajectoryTask != null;
        // if there isn't one
        if (trajectoryTask == null)
        {
            // create a new one
            trajectoryTask = new CTrajectoryTask();

            // add the new time course task to the task list
            // this method makes sure that the object is now owned 
            // by the list and that it does not get deleted by SWIG
            dataModel.getTaskList().addAndOwn(trajectoryTask);
        }

        // run a deterministic time course
        trajectoryTask.setMethodType(CCopasiMethod.deterministic);

        // pass a pointer of the model to the problem
        trajectoryTask.getProblem().setModel(dataModel.getModel());

        // actiavate the task so that it will be run when the model is saved
        // and passed to CopasiSE
        trajectoryTask.setScheduled(true);

        // set the report for the task
        trajectoryTask.getReport().setReportDefinition(report);
        // set the output filename.  The file is in workspace/immunity
        trajectoryTask.getReport().setTarget("intracellularTransportTimeCourse.txt");
        // don't append output if the file exists, but overwrite the file
        trajectoryTask.getReport().setAppend(false);
        
        // get the problem for the task to set some parameters
        CTrajectoryProblem problem = (CTrajectoryProblem)trajectoryTask.getProblem();

        // simulate 600 steps
        problem.setStepNumber(50);
        // start at time 0
        dataModel.getModel().setInitialTime(0.0);
        // simulate a duration of 60 time units
        problem.setDuration(50);
        // tell the problem to actually generate time series data
        problem.setTimeSeriesRequested(true);

        // set some parameters for the LSODA method through the method
        CTrajectoryMethod method = (CTrajectoryMethod)trajectoryTask.getMethod();

        CCopasiParameter parameter = method.getParameter("Absolute Tolerance");
        assert parameter != null;
        assert parameter.getType() == CCopasiParameter.DOUBLE;
        parameter.setDblValue(1.0e-12);
	}
	
	public void setInitialConcentration(String name, double value) {
		if (!nameMetabs.containsKey(name)) {
			System.out.println(name + "\t does not exist as a metab");
		} else {
			CMetab m = nameMetabs.get(name);
			m.setInitialConcentration(value);
			m.refreshInitialValue();
		}
	}
		
	public void runTimeCourse() {
		// reapply the initial values
		model.applyInitialValues();
		
		boolean result=true;
        String processError = "";
		String processWarning = "";
		try
        {
            // now we run the actual trajectory
        	System.out.println("trajectoryTask.process RABS_CONVERSION"+ dataModel.getObjectDisplayName());
            result=trajectoryTask.process(true);
            processError = trajectoryTask.getProcessError();
            processWarning = trajectoryTask.getProcessWarning();
        }
        catch (java.lang.Exception ex)
        {
            System.err.println( "Error. Running the time course simulation failed." );
            System.out.println(processError);
            System.out.println(processWarning);
            // check if there are additional error messages
            if (CCopasiMessage.size() > 0)
            {
                // print the messages in chronological order
                System.err.println(CCopasiMessage.getAllMessageText(true));
            }
            System.exit(1);
        }
        if(result==false)
        {
            System.err.println( "An error occured while running the time course simulation." );
            System.out.println(processError);
            System.out.println(processWarning);
            // check if there are additional error messages
            if (CCopasiMessage.size() > 0)
            {
                // print the messages in chronological order
                System.err.println(CCopasiMessage.getAllMessageText(true));
            }
            System.exit(1);
        }

	}
	
	public double getConcentration(String name) {
		double d = 0.0;
		
		if (!nameMetabs.containsKey(name)) {
		} else {
			CMetab m = nameMetabs.get(name);
			d = m.getConcentration();
		}
		return d;
	}
	public Set<String> getMetabolites(){
		Set<String> metabolites =  nameMetabs.keySet();
		return metabolites;
	}
	
	public CTrajectoryTask getTrajectoryTask() {
		return trajectoryTask;
	}
	
}
