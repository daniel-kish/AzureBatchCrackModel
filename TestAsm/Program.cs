using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Globalization;
using System.IO;
using System.IO.Compression;

namespace TestAsm
{
    class Program
    {
        static string LocalStoragePath = @"C:\Users\Daniel\Documents\Visual Studio 2013\Projects\StochasticBatch\rawResults";
        static void Main(string[] args)
        {
            using (System.IO.Compression.ZipArchive zip = ZipFile.Open("results.zip", ZipArchiveMode.Create))
            {
                for (int i = 0; i < 2; ++i)
                {
                    var inputFilter = string.Format(CultureInfo.InvariantCulture, String.Format("{0}_*", i));
                    HashSet<string> inputFiles = CollectFiles(LocalStoragePath, inputFilter);
                    
                    zip.CreateEntry(String.Format("{0}/", i));
                    foreach (var file in inputFiles)
                    {
                        zip.CreateEntryFromFile(file, String.Format("{0}/{1}", i, Path.GetFileName(file)));
                    }

                }
            }
        }
        private static void ZipOutputs(HashSet<string> inputs, string output)
        {
            System.IO.Compression.ZipArchive outputs = ZipFile.Open(output, ZipArchiveMode.Create);
            foreach (var input in inputs)
            {
                outputs.CreateEntryFromFile(input, Path.GetFileName(input), CompressionLevel.Optimal);
                
            }
        }

        private static HashSet<string> CollectFiles(string location, string inputFilter)
        {
            return new HashSet<string>(Directory.GetFiles(location, inputFilter));
        }
    }
}
