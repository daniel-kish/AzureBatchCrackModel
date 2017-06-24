# AzureBatchCrackModel

## About
This project is roughly 60% of my bachelor's thesis.
Here I use Microsoft Azure Batch service for parallel processing of a certain rupture mechanics simulations (elaborated [here](https://github.com/daniel-kish/CrackModel)). Main thing here is cloud apps, so-called cloud assemblies - ZIP archives that consist of certain dlls that implement the behaviour of a batch service that is to be deployed in the cloud at a great scale. You can read more about batch processing with Azure Batch [here](https://azure.microsoft.com/en-us/services/batch/).

## Work relevance
The model stored in the [CrackModel repo](https://github.com/daniel-kish/CrackModel) is used for sampling, i.e. performing repetitive calculation with small variations in modelling parameters. Obviously all simulations are independent and the task becomes an "embarrassingly parallel" one. The Azure Batch service provides you with a way of "acquiring" multiple virtual machines and loading your application to some cloud storage. Then a huge queue of tasks is generated and the VMs start processing them in so-called batch mode.
Here I use (slighlty outdated) "Azure Batch Apps" technique. Using that technique one provides so to speak parts of an existing service: 1. Job splitter. Generates a queue of tasks based on general problem statement, for instance what parameters of the model are to be changed and in what ranges.
2. TaskProcessor - they are executed on each VM and specify how exactly any given task is handled.
3. MergeProcess - defines how do we gather and maybe reduce the results.
All of this parts are specified in a dll that is uploaded to some Azure Storage instance. Then the service can be run, managed and monitored through a web interface.
In essence this project has helped me to conduct a massive research of the modeling technique and explore the parameter space with great efficiency and performance.


## Results
Most of the results are covered in the thesis itself. Based on that work a couple of papers was published:

The gist of the analysis is presented in proceedings here:

D. L. Kishlakov, P. V. Tarakanov, G. V. Shashurin, Y. V. Berchun, "Cloud Applications Performance in Crack Growth Simulations of Pre-Hydrogenated Structure Components", Materials Science Forum, Vol. 844, pp. 97-102, 2016. DOI:10.4028/www.scientific.net/MSF.844.97

And in more detail in a russian journal article here:

Д. Л. Кишлаков [D. L. Kishlakov], П. В. Тараканов, Г. В. Шашурин, и Ю. В. Берчун. 2017. «Эффективность облачных вычислений в моделировании кинетики трещин в наводороженных элементах конструкций.» Информационные Технологии 23 (2): 113-120. 
URL: http://novtex.ru/IT/it2017/number_2_annot.html#6

And as always, a plot for dessert:

![Plot](/коэффициент_ускорения_вязкость.png)
