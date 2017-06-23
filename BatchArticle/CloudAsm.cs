using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.IO.Compression;
using System.Globalization;
using Microsoft.Azure.Batch.Apps.Cloud;

namespace BatchArticle
{
    public class ApplicationDefinition
    {
        public static readonly CloudApplication Application = new ParallelCloudApplication
        {
            ApplicationName = "ArticleTest",
            JobType = "ArticleTestJob",
            JobSplitterType = typeof(Splitter),
            TaskProcessorType = typeof(Processor)
        };
    }

    public class Splitter : JobSplitter
    {
        protected override IEnumerable<TaskSpecifier> Split(IJob job, JobSplitSettings settings)
        {
            var tasks = new List<TaskSpecifier>();
            var pars = job.Parameters;

            double K_max = Double.Parse(pars["K_max"], CultureInfo.InvariantCulture);
            double sigma = Double.Parse(pars["sigma"], CultureInfo.InvariantCulture);
            double Omega = Double.Parse(pars["Omega"], CultureInfo.InvariantCulture);

            double l0_init = Double.Parse(pars["l0_init"], CultureInfo.InvariantCulture);
            double l0_end = Double.Parse(pars["l0_end"], CultureInfo.InvariantCulture);
            double l0_steps = Double.Parse(pars["l0_steps"], CultureInfo.InvariantCulture);

            double l0_step = (l0_end-l0_init)/l0_steps;

            for (double l0 = l0_init; l0 <= l0_end; l0 += l0_step)
            {
                tasks.Add(makeTask(K_max, sigma, Omega, l0, String.Format("{0}_{1}_{2}_{3}", K_max, sigma, Omega, l0)));
            }

            return tasks;
        }

        private double[] parseParamVals(string line)
        {
            string[] valStrings = line.Split(' ');
            double[] vals = new double[valStrings.Length];

            for (var i = 0; i < valStrings.Length; ++i)
            {
                vals[i] = Double.Parse(valStrings[i], CultureInfo.InvariantCulture);
            }

            return vals;
        }
        private TaskSpecifier makeTask(double K_max, double sigma, double Omega, double l0, string index)
        {
            var pars = new Dictionary<string, string>() 
            {
                {"K_max", K_max.ToString(CultureInfo.InvariantCulture) },
                {"sigma", sigma.ToString(CultureInfo.InvariantCulture) },
                {"Omega", Omega.ToString(CultureInfo.InvariantCulture) },
                {"l0", l0.ToString(CultureInfo.InvariantCulture) },
                {"index", index}
            };
            return new TaskSpecifier
            {
                Parameters = pars
            };
        }
    }
    public class Processor : ParallelTaskProcessor
    {
        protected override TaskProcessResult RunExternalTaskProcess(ITask task, TaskExecutionSettings settings)
        {
            var parameters = task.Parameters;
            var outFile = String.Format("{0}.txt", parameters["index"]);

            var process = new ExternalProcess
            {
                CommandPath = ExecutablePath(@"Model.exe"),
                Arguments = String.Format("{0} {1} {2} {3} {4}",
                parameters["l0"],
                parameters["sigma"],
                parameters["K_max"],
                parameters["Omega"],
                outFile),
                WorkingDirectory = LocalStoragePath
            };

            try
            {
                ExternalProcessResult processOutput = process.Run();
                return TaskProcessResult.FromExternalProcessResult(processOutput, outFile);
            }
            catch (ExternalProcessException ex)
            {
                string outputInfo = "No program output";
                if (!string.IsNullOrEmpty(ex.StandardError) || !string.IsNullOrEmpty(ex.StandardOutput))
                {
                    outputInfo = Environment.NewLine + "stderr: " + ex.StandardError + Environment.NewLine + "stdout: " + ex.StandardOutput;
                }

                Log.Error("Failed to invoke command {0} {1}: exit code was {2}.  {3}", ex.CommandPath, ex.Arguments, ex.ExitCode, outputInfo);
            }
            catch (Exception ex)
            {
                Log.Error("Error in task processor: {0}", ex.ToString());
            }
            return new TaskProcessResult { Success = TaskProcessSuccess.RetryableFailure };
        }

        protected override JobResult RunExternalMergeProcess(ITask mergeTask, TaskExecutionSettings settings)
        {
            var inputFilter = string.Format(CultureInfo.InvariantCulture, "*.txt");
            var inputFiles = CollectFiles(LocalStoragePath, inputFilter);

            var completionFile = LocalPath("results.zip");

            var result = ZipOutputs(inputFiles, completionFile);

            return result;
        }

        private JobResult ZipOutputs(HashSet<string> inputs, string output)
        {
            try
            {
                using (System.IO.Compression.ZipArchive outputs = ZipFile.Open(output, ZipArchiveMode.Create))
                {
                    foreach (var input in inputs)
                    {
                        outputs.CreateEntryFromFile(input, Path.GetFileName(input), CompressionLevel.Optimal);
                    }
                }

                return new JobResult { OutputFile = output };
            }
            catch (Exception ex)
            {
                var error = string.Format("Failed to zip outputs: {0}", ex.ToString());
                return null;
            }
        }

        private HashSet<string> CollectFiles(string location, string inputFilter)
        {
            return new HashSet<string>(Directory.GetFiles(location, inputFilter));
        }
    }
}
