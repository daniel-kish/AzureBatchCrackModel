using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.IO;
using System.Globalization;
using Microsoft.Azure.Batch.Apps.Cloud;
using System.IO.Compression;

namespace ModelBatchTest
{
    public class ApplicationDefinition
    {
        public static readonly CloudApplication Application = new ParallelCloudApplication
        {
            ApplicationName = "ModelBatchTest",
            JobType = "ModelBatchTest",
            JobSplitterType = typeof(MySplitter),
            TaskProcessorType = typeof(MyProcessor)
        };
    }

    public class MySplitter : JobSplitter
    {
        protected override IEnumerable<TaskSpecifier> Split(IJob job, JobSplitSettings settings)
        {
            var tasks = new List<TaskSpecifier>();
            double low = Double.Parse(job.Parameters["low"]);
            double high = Double.Parse(job.Parameters["high"]);
            int num = Int32.Parse(job.Parameters["steps"]);
            double omega = 1.5;

            if (job.Parameters.Select(par => par.Key).Contains("omega"))
                omega = Double.Parse(job.Parameters["omega"]);

            double step = (high - low) / num;

            for (double x = low; x <= high; x += step)
            {
                tasks.Add(new TaskSpecifier{
                    Parameters = new Dictionary<string,string>{
                            {"l0", x.ToString()},
                            {"omega", omega.ToString()}
                        }
                });
            }
            return tasks;
        }
    }

    public class MyProcessor : ParallelTaskProcessor
    {
        protected override TaskProcessResult RunExternalTaskProcess(ITask task, TaskExecutionSettings settings)
        {
            var parameters = task.Parameters;
            var outFile = String.Format("{0}_{1}.txt", "out", task.TaskId);

            var process = new ExternalProcess{
                CommandPath = ExecutablePath(@"ModelBatchTest/Model.exe"),
                Arguments = String.Format("{0} {1} {2}", parameters["l0"], parameters["omega"], outFile),
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
            var inputFilter = string.Format(CultureInfo.InvariantCulture, "{0}*", "out");
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
