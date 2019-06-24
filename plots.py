from sim_with_cells import simulation_with_cells as sm
from cell_class import cell
from helix_class import chromosome_matrix as cm
import matplotlib.pyplot as plt
import math
import random
import math

plt.style.use('seaborn-whitegrid')
# base_pairs = 8000
# probability = []
# dep_var = list(range(base_pairs, 0, -1))
# for i in dep_var:
#     B = 0.008561 #max height
#     C = 1776.6    # marks the point of the peak
#     D = 0.56846
#     probability.append(B * math.exp(-((math.log(i / C)) / D) ** 2))
#
#
# x = dep_var
# y = probability
# print(x)
# print(y)
# plt.title("probability distribution B = 8e-3")
# plt.plot(x, y, '.', color='blue')
# plt.ylabel('probability')
# plt.xlabel("shortest telomere length")
# plt.show()
bp_4000_length = [3937.388586956522, 3874.6277173913045, 3813.1374464307146, 3750.3684246915846, 3690.9354778661855, 3630.3148299829895, 3568.0101578312724, 3506.5992310611296, 3444.7526155104833, 3384.4061981240543, 3326.9186758853803, 3273.643433014885, 3222.6097062527883, 3177.402986953179, 3135.6130523483102, 3102.4463816122916, 3071.5088759634755, 3046.6734768690244, 3025.9406209972108, 3010.6516254311746, 2997.2001608306655, 2987.2186409416104, 2979.3544577706757, 2974.092824405338, 2969.5044954360956, 2967.6233166967754, 2964.0141658612147, 2962.138171154788, 2961.358156624023, 2959.889383942288, 2959.252354560907, 2959.0246996445508, 2958.9164074994455, 2958.659491277757, 2958.3380674139967, 2958.6565475352636, 2958.7642416205063, 2958.623199804097, 2958.7573712418284, 2958.4072006357487, 2958.0374486218093, 2958.2068296394136, 2957.35380425841, 2957.267198549823, 2957.2284058567384, 2956.5986709945514, 2957.317880668688, 2957.1243886105085, 2957.604800345771, 2957.7511951661477, 2957.9892408282794, 2957.8456278457434, 2957.7787622400256, 2958.179295404844, 2958.411354898692, 2958.3826632086952, 2957.9794211280805, 2957.8087749944734, 2957.8379547589775, 2957.8377991340785, 2957.8197157655663, 2957.586241106451, 2957.971519727592, 2957.299524562705, 2957.3477098947274, 2957.698839838288, 2957.8624175062055, 2957.507593489075, 2957.6195885761094, 2957.2536837892617, 2957.610723347953, 2957.615782157909, 2957.158227617553, 2957.9659615689015, 2958.0793051884584, 2958.8980046808956, 2958.5388857963076, 2958.4087703071455, 2958.5229934471654, 2958.224100523285, 2957.6730335467764, 2957.4537822710968, 2957.761218516913, 2957.646040537859, 2957.6461135218456, 2957.472328311533, 2957.4963750105603, 2957.059923295012, 2957.081726094965, 2957.073611083739, 2957.5397111445836, 2957.8580334710464, 2957.8430744690686, 2957.9916807829113, 2957.447474532139, 2957.6289890747307, 2958.0883529707276, 2958.1090707372823, 2958.9417097774585, 2959.4385624705224, 2959.1201076253305, 2958.562759856769]
# bp_4000_doublings = [1.0, 2.0, 3.0, 4.0, 4.954196310386876, 5.906890595608519, 6.906890595608519, 7.894817763307945, 9.374903629112072, 10.617604109603175, 11.711497654656995, 12.741620039262795, 13.758060395600001, 14.750506810804353, 15.715407616781166, 16.685030777975914, 17.654186009589537, 18.605938168122538, 19.573176272829702, 20.543999502723665, 21.518057560797413, 22.465142741470277, 23.432862787087107, 24.40822775107181, 25.386003757857157, 26.3706872220216, 27.348851712273124, 28.33085247209122, 29.29485126238463, 30.26112696517908, 31.24997459622127, 32.218835341417446, 33.171887223347326, 34.141178206313654, 35.07581553185454, 36.067307345745995, 37.0188643924927, 38.011763087075416, 38.974069202446785, 39.92520039722153, 40.90070512371671, 41.869975462385334, 42.84323871730197, 43.83278771713424, 44.801996734169535, 45.75522447151918, 46.7460414076744, 47.69447967458064, 48.66856452507628, 49.64997216138071, 50.62506613242026, 51.57413831015484, 52.52278896283669, 53.50423312710345, 54.44787228752209, 55.4425115691304, 56.41424064370227, 57.39016067085036, 58.339644972647385, 59.284209390931515, 60.2621507209563, 61.24123759395478, 62.21373978958114, 63.19433685369274, 64.14333529741828, 65.13058283602484, 66.07148677865302, 67.03458180275344, 67.99757989924666, 68.96482692336976, 69.96139815314106, 70.9398100006655, 71.89628508082193, 72.89014520163092, 73.83301192192272, 74.81566554955307, 75.75710346963794, 76.73767351099615, 77.67228814985339, 78.64655751920033, 79.64383492624457, 80.60548378703372, 81.58200104332529, 82.5574846019138, 83.51531986278374, 84.44984661961877, 85.40683828068015, 86.36388346299472, 87.3596972626908, 88.3405725123091, 89.29399319702037, 90.26776665856696, 91.21041635272981, 92.1873371731678, 93.13584906424245, 94.09140878803204, 95.07354065730966, 96.03470246143105, 97.01817374000429, 97.97716949658043, 98.95093333756799, 99.92360783523712]
# plt.scatter(bp_4000_doublings, bp_4000_length, color='red')
# bp_2000_length =[2000.0, 2020.8151688627804, 2039.383520127239, 2055.8429160169317, 2075.1961188363975, 2092.07320619887, 2111.0405449704003, 2128.7866745826936, 2146.93102453362, 2164.839761351793, 2183.1418485943736, 2201.9051894817144, 2219.4483420484607, 2237.9986849350084, 2255.069608866545, 2271.385358910261, 2289.437626447393, 2308.23463747403, 2324.8729983718963, 2342.497330922911, 2361.1437020188596, 2377.4497894206597, 2394.9234478434482, 2413.975009641264, 2432.5098035660976, 2450.450359961755, 2467.756572646173, 2486.5802550875733, 2507.589062640007, 2525.566472892985, 2543.319040334578, 2558.6108601313103, 2578.3286133205665, 2597.0259089963783, 2612.5360224099677, 2630.359042899876, 2646.79970295175, 2661.0216217459415, 2673.6468066237794, 2689.686719056105, 2708.1178219583935, 2722.356683084549, 2736.9902496138407, 2751.7642209503456, 2763.156254800365, 2777.5821310516867, 2788.07610044135, 2800.558346874763, 2814.3341735515623, 2823.386882372908, 2836.2494856369144, 2844.4773682982973, 2854.5047930546552, 2863.748495234845, 2872.4500772056163, 2879.407903506289, 2889.18568304435, 2896.325910971659, 2905.254662268029, 2912.3112185880145, 2918.3278708372427, 2924.394838007284, 2927.992444827928, 2930.8484916784987, 2934.644670376086, 2936.7700453101515, 2940.088038553412, 2941.8461826111634, 2943.0123873905095, 2944.3468059689903, 2946.183598774351, 2946.4448319556614, 2947.0235450182035, 2948.461565794327, 2950.062749211929, 2951.2986203824094, 2952.2485016158207, 2953.04066001917, 2953.625022434865, 2954.265960327613, 2955.140239966049, 2954.9053126009617, 2954.8111228656494, 2955.614991762071, 2956.016586431978, 2955.9835776389978, 2955.8241271715656, 2956.3722567201526, 2956.3970441990805, 2956.4544977215724, 2957.126826494896, 2957.4445499576495, 2957.3490334682288, 2957.7111927723786, 2958.07325071802, 2957.971929883947, 2957.6364811703415, 2957.9902038453565, 2957.6449007471642, 2957.5376288390084, 2958.325183085975, 2958.4257725278217, 2957.9260898097377, 2957.8114408893753, 2958.1181206920633, 2957.971770912887]
#
# bp_2000_doublings =[0.0, 1.0, 2.0, 2.807354922057604, 3.8073549220576037, 4.700439718141093, 5.643856189774724, 6.554588851677638, 7.475733430966398, 8.502956083262315, 9.78218079392331, 10.847946885011506, 11.855268842976466, 12.834130299581885, 13.750851809033675, 14.646785519237017, 15.594074765641665, 16.52319899137956, 17.44008318885459, 18.343312771781402, 19.221823188656703, 20.166844961653798, 21.069892811357406, 22.00111833583068, 22.93806467374107, 23.832653739381787, 24.71395171738445, 25.648585202186187, 26.52138395444535, 27.4309976083824, 28.365799700389246, 29.313401357825068, 30.27320097797515, 31.205710559167365, 32.11264012807812, 33.04998142876618, 33.98994140980251, 34.88906598923199, 35.83696818631361, 36.75469565858924, 37.683102092053005, 38.64176346112438, 39.59109141495807, 40.50814346777884, 41.45498922722341, 42.39105761233194, 43.32739789903616, 44.294300350318395, 45.2477695360463, 46.15687466658955, 47.09806611600493, 48.0602072507268, 49.009826470175724, 49.967498789143995, 50.90618943326701, 51.8206972802546, 52.76616277261972, 53.73625187727368, 54.714369053444216, 55.64475338274205, 56.5708555841448, 57.50924107373655, 58.4958920699224, 59.453534341087725, 60.415225254554414, 61.38040151872424, 62.31862813520782, 63.27758457686651, 64.23564087310118, 65.1846079001364, 66.16919631935957, 67.13158147698941, 68.1028094071151, 69.08136435415766, 70.02786996114722, 71.012468589581, 71.97631220611339, 72.91533866012513, 73.89329746184605, 74.86380866648474, 75.79475257275834, 76.7614153832651, 77.72193428186885, 78.70371408972548, 79.69472098577982, 80.64753450747556, 81.60382720690109, 82.55761936439197, 83.5474133405894, 84.48526462801114, 85.48407550062865, 86.43361877089609, 87.42115637641506, 88.37506319595737, 89.31903223051985, 90.27082786120181, 91.21947085107634, 92.20388698889381, 93.17477843254012, 94.1574489363971, 95.14173793074985, 96.09832279658922, 97.0793992841025, 98.02434917305466, 98.99523882860314, 99.9478649225106]
#
# plt.scatter(bp_2000_doublings, bp_2000_length, color = 'blue')
plt.plot(bp_4000_length, list(range(len(bp_4000_length))) )
plt.show()

# a_80 = [53.86826579445115, 57.70854034338113, 43.523561956057016, 58.608151038316926, 54.42867680168123, 56.15269300325011, 55.069699060250755, 53.64260113462012, 53.60832751848808, 52.28484463146517, 55.28806806916803, 52.706074858490375, 57.420600446638566, 58.76335826711987, 60.417784104295606, 59.99595872178685, 46.27530134938641, 56.03623926861076, 52.86851085836423, 58.00506019001829, 54.14492894357066, 54.17938298356164, 57.06320059990486, 52.2926437795485, 54.68495454176393, 54.834946610534686, 59.49070159431186, 59.67490080611354, 43.39231742277876, 61.23104352304039, 42.906890595608516, 53.926575354138144, 57.14177991603211, 40.0, 43.32192809488736, 52.87224855579928, 58.28202136598374, 59.33359119397369, 52.642702313957486, 50.74411884418107, 54.655853037603194, 39.0, 53.706339284992126, 53.29263327708281, 52.6514036257547, 56.18414839359967, 40.0, 57.793589128188955, 55.11285623692198, 53.15441492819169, 52.328773904074374, 58.817657982771536, 43.523561956057016, 56.957397632784335, 54.29702769985087, 51.670878895627396, 54.73749738808554, 42.0, 54.44046531493142, 56.20198741999741, 57.715174801278636, 54.758453668930855, 57.137075437862336, 57.53527162356123, 53.60298391742373, 57.31346825478068, 57.46137885531109, 58.226448445627675, 40.0, 58.24286959851434, 53.932221655731595, 53.19230461978367, 58.76196139159847, 52.63629705564843, 42.169925001442316, 54.9519359034663, 51.705571369135754, 53.94665402972507, 54.515098653169325, 55.5299278573127, 52.38762672149084, 53.888533777878905, 54.98461969386739, 55.03431340468384, 55.710388581098414, 59.17859254819085, 56.96441020562552, 56.47838277699788, 58.365327193890714, 57.840290745489064, 59.71623499828289, 54.65091055942931, 42.0, 56.1950496325678, 40.0, 52.515103084566896, 40.0, 59.08945392798377, 55.23977471637346, 56.1674998955752]x