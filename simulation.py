from sim_with_cells import simulation_with_cells as sm
from cell_class import cell
from helix_class import chromosome_matrix as cm
import matplotlib.pyplot as plt
import statistics
import random
plt.style.use('seaborn-whitegrid')
base_pairs = 5500
# # root array is an array of 46 chromosome matrices
array = []
for i in range(100):
    root_array = [cm(base_pairs, base_pairs, base_pairs, base_pairs) for x in range(46)]
    # root cell will be used to initialize the simulation class
    root_cell = cell(root_array)
    full_sim = sm(root_cell, 7, 0)
    x = full_sim.start()
    array.append(x)
    print(i+1, x)
print(array)
# a_80 = []
# for i in range(50):
#     root_array = [cm(base_pairs, base_pairs, base_pairs, base_pairs) for x in range(46)]
#     root_cell = cell(root_array)
#     full_sim = sm(root_cell,6,0)
#     x = full_sim.start()
#     print(x)
#     a_80.append(x)
#     print(i)
# print(a_80)
# x = full_sim.population_doublings_array
# y = full_sim.percent_array
# print(y)
# plt.plot(x, y, 'o', color='black')
#
#
# # array = inst.start()
# # mother_cell = random.choice(array)
# # mother_sim = sm(mother_cell, 8, 16)
# # cell_culture = mother_sim.start()
# # iteration_array = []
# # for seed_cell in cell_culture:
# #     seed_sim = sm(seed_cell, 7, 0)
# #     # because of the upper bound, the start method will return a 200 cell array
# #     iteration_array.append(seed_sim.start())
# print (array)
# a_80 = [70.63419487174218, 69.58545967926322, 70.47801822215771, 69.68448202325138, 70.36964640376215, 70.6332808613458, 69.77950083300152, 69.78095810633718, 69.84856767017862, 69.33009801201172, 68.92525784731208, 69.37401700050928, 70.99786715097176, 71.53489462553851, 69.55794000897443, 69.54199137112721, 70.85775670213044, 69.63631346198822, 70.44934482713312, 70.44366772367026, 70.5809468965085, 70.50732985032255, 70.27809638137555, 70.09820101602192, 70.34653725837144, 70.78264348341263, 70.12060980385579, 70.8072330942108, 69.55815993704203, 69.34837797, 70.1194956115276, 69.96809802476272, 69.63829161658121, 70.46215068694724, 70.454950382137, 69.18744533597571, 68.89437306880058, 70.23433595662178, 70.08560863848808, 69.71634814309944, 69.60635267254776, 69.93729161260005, 70.57777032576186, 70.51598058396841, 70.7546268271515, 69.80841448411518, 68.95955190770353, 70.54858343968586, 70.63563779417119, 69.21219754700878]
# a_86 = [69.53144632517461, 68.30438401200165, 70.18818775012454, 69.39341146895187, 70.25508415372131, 70.64990859437984, 69.56788080405812, 69.00726728571337, 70.14671064566382, 70.12555114162156, 69.22057162511923, 69.92243655219569, 68.89524708422066, 70.52240061913551, 70.86137932158101, 69.55407686830824, 69.30760927552865, 68.8504319350982, 70.21072801220417, 69.33372268128235, 69.9341119400557, 69.58998764504007, 69.67770652152308, 70.10808280964741, 70.22300353919427, 68.84002592430662, 69.11620869705014, 67.69434221595696, 69.9626488941123, 68.54232712002475, 69.97332878741065, 70.22685475262816, 69.91876016505745, 68.92765618097796, 69.42276155464305, 69.66377774538117, 69.00389450328038, 69.66929294653436, 70.97710636514356, 69.22057934910981, 69.78335831343293, 69.37057629913794, 69.13518863453216, 70.08279251368239, 70.26658549693641, 70.28418671635045, 69.07565438529264, 70.45007321380706, 70.5471084602833, 69.95301062541908]
# a_90 = [69.2356470926902, 69.90807562950805, 70.17059103273235, 69.18283617928022, 69.23079691589871, 68.47577076449124, 69.3006703458134, 69.11666817000938, 69.52388356424278, 69.5608146364147, 69.63543651215353, 69.23929028438168, 69.51908653177604, 70.92492404163463, 70.21682319314027, 68.48270117449404, 70.26287618697587, 70.19553773519681, 69.30053065092898, 68.92275911743401, 70.07803831243689, 69.68174544293744, 69.73807761708011, 70.23719187621927, 69.46036550851946, 69.618612522305, 69.02932760743955, 68.86007888234712, 70.70737503014706, 70.14540695188127, 70.95640828074522, 68.76273764510205, 69.49287467652613, 71.58981332813396, 69.47518719256679, 69.41622507653302, 69.83580900146161, 70.36680877062145, 70.73650005430872, 70.17656390176437, 69.79435620152297, 69.5066661873792, 68.66379883472489, 69.4611601135736, 69.51457951222442, 69.67652798921141, 68.99915678565925, 70.47929115177214, 68.9760562174267, 69.98236150022903]
# a_95 = [69.42999475679466, 69.49134315916808, 69.07815400847066, 68.14339676880387, 69.4433683452424, 69.59607137237221, 69.92629518206924, 69.18679198584391, 69.1438304987934, 69.45075345819572, 69.05115317195255, 69.12609422933366, 68.6149412002804, 68.75279663322671, 68.93221550319691, 68.93130411140227, 69.3870292900244, 69.02310079799511, 69.05572103033951, 68.76813578164484, 69.21434962955409, 69.28245323759279, 69.75514371268426, 70.03692286764314, 69.56827280643898, 68.76278431712196, 69.30008054872336, 69.21783662234611, 69.96804128201686, 69.36200916778564, 69.00669112895802, 69.50173611459995, 68.47151281973981, 71.06604881252301, 69.06949686151263, 69.23393120367592, 70.61735588624234, 69.61395562560848, 68.93840497929382, 69.47308828850407, 69.59484768133323, 70.01431069084936, 68.4859082786629, 69.88764689838887, 69.19958841968322, 68.44758317346957, 69.26687600384366, 68.0597849911771, 69.2988153691759, 70.70099301240083]

# a_80 = [68.56275263847103, 69.22421953107377, 68.30930576719713, 68.86855035090798, 69.36956906569013, 69.27691950481548, 69.18596998674502, 67.92628854191958, 68.04987991069689, 70.08584950265148, 70.1232391128664, 70.23906229389704, 67.87665335021441, 69.04814942677535, 69.35017404643962, 69.59983734238298, 68.59684889143881, 69.45029651047501, 68.6755862187948, 68.98440298544311, 68.64131467998595, 68.97809670878269, 67.97324483662345, 69.30363702573365, 68.70504426487535, 68.11391931539373, 68.61400748611325, 69.26442817352866, 67.49782981614806, 68.91993170377707, 66.89803306650768, 68.48756304767201, 69.27534125282223, 68.30337445274306, 69.66302797085983, 68.88092651547755, 69.05607916986284, 68.68279800244785, 69.2802912101859, 69.05524378806403, 69.59584633908187, 68.18428218124119, 68.65806192763898, 68.71783121438827, 69.91345523953164, 69.72385164448677, 68.47323761119986, 69.28466790414292, 69.98096154267725, 69.12640270862336]

plt.hist(array, bins=20, alpha=0.80, label='alpha = 0.80')
# plt.hist(a_90, bins=10, alpha=0.5, label='alpha = 0.90')
# plt.hist(a_86, bins=10, alpha=0.5, label='alpha = 0.86')
# plt.hist(a_80, bins=10, alpha=0.5, label='alpha = 0.80')
# plt.hist(a_85, bins= 100,range= [0,120], alpha=0.5, label='alpha = 0.85')
# plt.hist(a_90, bins=100, range= [0,120], alpha=0.5, label='alpha = 0.9')
# plt.hist(a_95, bins=100,range= [0,120], alpha=0.5, label='alpha = 0.95')
# x = [0.8,0.85,0.9,0.95]
# y = [statistics.mean(a_80),statistics.mean(a_85),statistics.mean(a_90),statistics.mean(a_95)]
#
# plt.plot(x, y, 'o', color='blue');
#
plt.ylabel('count')
plt.xlabel("PDs")
plt.legend(loc='upper left')
plt.show()




#####################  Things to do #######################

# run the 200 cell culture by taking the chromosome that has the shortest telomere and  copying to make a new cell...


##New changes, copy the cell to the next level, but if it is senescent then don't unwind it -> this will help in the calculations
## Use PDs with the new jaggi formula that multiplies by

#things to do: implement bimodal distribution and compare it with the resampling procedure and thejust letting a cell run.
# When implementing the abrupt shortening if the cell divides (I believe considering prob of replication and senescence) then replicate all matrices without shortening
# and decide for abrupt shortening with nick's formula. In the future we want to vary the parameters and test for the varying peaks -> In the decision making process
# if l = 300 then the probability is 0. If the matrix undergoes abrupt shortening then
# Let the program run for alpha values incrementing from 0.8 to 0.94