% Calculate inertial parameters regressor of coriolis matrix for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% cmat_reg [(6*6)x(6*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RRPPRR8_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:26:16
% EndTime: 2019-03-09 09:27:03
% DurationCPUTime: 36.87s
% Computational Cost: add. (24373->983), mult. (49552->1292), div. (0->0), fcn. (53325->8), ass. (0->706)
t903 = qJD(5) + qJD(6);
t706 = sin(qJ(5));
t708 = cos(qJ(5));
t709 = cos(qJ(2));
t1053 = pkin(8) - qJ(3);
t707 = sin(qJ(2));
t1115 = t1053 * t707 - pkin(1);
t704 = cos(pkin(10));
t779 = t704 * t1115;
t703 = sin(pkin(10));
t780 = t703 * t1115;
t1062 = pkin(7) * t703;
t891 = -pkin(3) - t1062;
t802 = pkin(2) * t704 - t891;
t781 = pkin(4) + t802;
t783 = -pkin(2) * t703 + pkin(7) * t704 - qJ(4);
t294 = -t706 * t779 + t708 * t780 + (t706 * t781 + t708 * t783) * t709;
t992 = t708 * t707;
t995 = t706 * t707;
t573 = -t703 * t992 + t704 * t995;
t215 = -t573 * pkin(9) + t294;
t705 = sin(qJ(6));
t1008 = t705 * t215;
t1064 = cos(qJ(6));
t733 = t708 * t781;
t745 = t706 * t783;
t714 = t745 - t733;
t712 = (pkin(5) - t714) * t709;
t628 = t706 * t703 + t708 * t704;
t575 = t628 * t707;
t1060 = t575 * pkin(9);
t743 = t708 * t779;
t744 = t706 * t780;
t715 = t744 + t743;
t713 = -t715 - t1060;
t710 = t712 + t713;
t195 = t1064 * t710;
t111 = -t195 + t1008;
t1127 = t709 * t714;
t214 = t713 - t1127;
t888 = t1064 * t214;
t116 = t888 - t1008;
t1166 = t111 + t116;
t996 = t706 * t704;
t631 = t708 * t703 - t996;
t1000 = t705 * t631;
t606 = t1064 * t628;
t1117 = t606 + t1000;
t1097 = t1117 / 0.2e1;
t1167 = t709 * t1097;
t991 = t708 * t709;
t663 = t703 * t991;
t994 = t706 * t709;
t576 = t704 * t994 - t663;
t1003 = t705 * t576;
t577 = t628 * t709;
t884 = t1064 * t577;
t984 = t1003 / 0.2e1 - t884 / 0.2e1;
t255 = t1167 + t984;
t1174 = t255 * qJD(2);
t1159 = t1167 - t984;
t1004 = t705 * t575;
t553 = t1064 * t573;
t1119 = t553 + t1004;
t952 = qJD(1) * t709;
t879 = t1119 * t952;
t1169 = qJD(2) * t1159 + t879;
t1012 = t704 * t707;
t1016 = t703 * t707;
t993 = t707 * qJ(3);
t845 = -pkin(1) - t993;
t825 = t704 * t845;
t826 = t703 * t845;
t887 = t1064 * t215;
t112 = t887 + (-t706 * (pkin(8) * t1016 + t826) + t708 * (-pkin(8) * t1012 - t825) - t1060 + t712) * t705;
t1009 = t705 * t214;
t115 = -t887 - t1009;
t1165 = t112 + t115;
t1001 = t705 * t628;
t607 = t1064 * t631;
t1116 = t607 - t1001;
t1147 = t1116 * t709;
t1002 = t705 * t577;
t555 = t1064 * t576;
t828 = -t555 / 0.2e1 - t1002 / 0.2e1;
t1157 = t1147 / 0.2e1 + t828;
t1173 = qJD(1) * t1157;
t1172 = qJD(1) * t1159;
t1005 = t705 * t573;
t554 = t1064 * t575;
t1118 = t554 - t1005;
t1129 = t1118 ^ 2;
t1152 = t1119 ^ 2;
t1161 = -t1129 + t1152;
t1171 = qJD(1) * t1161;
t1158 = -t1147 / 0.2e1 + t828;
t1163 = qJD(2) * t1158;
t1130 = t1116 ^ 2;
t1153 = t1117 ^ 2;
t1160 = -t1130 + t1153;
t1170 = qJD(2) * t1160;
t881 = t1118 * t952;
t1156 = qJD(2) * t1157 + t881;
t848 = t1053 * t704;
t849 = t1053 * t703;
t510 = -t706 * t848 + t708 * t849;
t421 = -t631 * pkin(9) - t510;
t722 = t705 * t421;
t511 = -t706 * t849 - t708 * t848;
t422 = -t628 * pkin(9) + t511;
t885 = t1064 * t422;
t237 = t885 + t722;
t1168 = -t237 / 0.2e1;
t1067 = -t709 / 0.2e1;
t1164 = t1067 * t237;
t1150 = -t1117 / 0.2e1;
t1162 = t1118 * t1150;
t1137 = t1150 * t1119;
t1151 = -t1116 / 0.2e1;
t1138 = t1151 * t1118;
t1066 = t709 / 0.2e1;
t638 = -t704 * pkin(3) - t703 * qJ(4) - pkin(2);
t610 = t704 * pkin(4) - t638;
t508 = pkin(5) * t628 + t610;
t1093 = -t508 / 0.2e1;
t236 = t1064 * t421 - t705 * t422;
t1155 = t1066 * t236 - t1093 * t1119;
t1154 = t903 * t236;
t883 = t1064 * t706;
t998 = t705 * t708;
t634 = t883 + t998;
t1080 = t634 / 0.2e1;
t882 = t1064 * t708;
t999 = t705 * t706;
t768 = -t882 + t999;
t1082 = -t768 / 0.2e1;
t1149 = -t1118 / 0.2e1;
t1101 = t1118 / 0.2e1;
t1148 = t1119 / 0.2e1;
t950 = qJD(2) * t1116;
t874 = t1117 * t950;
t694 = t707 * qJ(4);
t677 = t704 * t694;
t830 = -pkin(7) + (-pkin(3) - pkin(4)) * t703;
t514 = t707 * t830 + t677;
t399 = pkin(5) * t573 + t514;
t1145 = t1118 * t399;
t1144 = t1119 * t399;
t957 = qJD(1) * t1118;
t880 = t1119 * t957;
t146 = t1119 * t1151 + t1162;
t1143 = t146 * qJD(6);
t1085 = -t631 / 0.2e1;
t1086 = -t628 / 0.2e1;
t320 = t1085 * t573 + t1086 * t575;
t1142 = t320 * qJD(5);
t1092 = -t554 / 0.2e1;
t540 = t554 / 0.2e1;
t1120 = t1092 + t540;
t1141 = qJD(5) * t1120;
t1088 = -t607 / 0.2e1;
t589 = t607 / 0.2e1;
t1121 = t1088 + t589;
t1140 = qJD(5) * t1121;
t1134 = -qJD(1) * t146 + t874;
t946 = qJD(2) * t631;
t1133 = -qJD(1) * t320 + t628 * t946;
t1132 = -qJD(2) * t146 + t880;
t953 = qJD(1) * t575;
t1131 = -qJD(2) * t320 + t573 * t953;
t853 = -t722 / 0.2e1;
t1094 = t1116 / 0.2e1;
t844 = qJD(5) + t952;
t1128 = t573 * t844;
t701 = t707 ^ 2;
t702 = t709 ^ 2;
t679 = t702 - t701;
t1075 = t704 / 0.2e1;
t1077 = -t703 / 0.2e1;
t1015 = t703 * t709;
t824 = -pkin(2) * t709 - t993;
t784 = -pkin(1) + t824;
t551 = -pkin(7) * t1015 + t704 * t784;
t1011 = t704 * t709;
t902 = pkin(7) * t1011;
t552 = t703 * t784 + t902;
t775 = t1075 * t552 + t1077 * t551;
t1055 = t709 * pkin(7);
t892 = t1055 / 0.2e1;
t341 = t892 - t775;
t699 = t703 ^ 2;
t700 = t704 ^ 2;
t671 = t700 + t699;
t640 = t671 * qJ(3);
t907 = t640 * qJD(2);
t1125 = -qJD(1) * t341 + t907;
t516 = t709 * t783 + t826;
t1023 = t516 * t704;
t678 = qJ(4) * t1011;
t1078 = -t678 / 0.2e1;
t517 = t709 * t802 - t825;
t264 = t892 + t1078 - t1023 / 0.2e1 + (pkin(3) * t1066 - t517 / 0.2e1) * t703;
t1124 = -qJD(1) * t264 + t907;
t944 = qJD(2) * t707;
t871 = t704 * t944;
t627 = t679 * t703;
t909 = t627 * qJD(1);
t1123 = t909 + t871;
t793 = t950 + t957;
t943 = qJD(2) * t709;
t933 = qJD(6) * t1116;
t1122 = qJD(6) * t1118;
t1114 = t573 ^ 2;
t1113 = t575 ^ 2;
t1112 = t628 ^ 2;
t1111 = t631 ^ 2;
t1110 = pkin(7) / 0.2e1;
t1109 = -t112 / 0.2e1;
t1107 = -t236 / 0.2e1;
t1105 = t237 / 0.2e1;
t1103 = -t1119 / 0.2e1;
t538 = -t553 / 0.2e1;
t539 = t553 / 0.2e1;
t1091 = -t575 / 0.2e1;
t1090 = t575 / 0.2e1;
t862 = t577 / 0.2e1;
t600 = t634 * t709;
t1089 = -t600 / 0.2e1;
t587 = -t606 / 0.2e1;
t588 = t606 / 0.2e1;
t1087 = t610 / 0.2e1;
t1084 = t631 / 0.2e1;
t1083 = t768 / 0.2e1;
t1081 = -t634 / 0.2e1;
t1079 = t663 / 0.2e1;
t1076 = t703 / 0.2e1;
t1074 = -t705 / 0.2e1;
t1073 = t705 / 0.2e1;
t1072 = -t706 / 0.2e1;
t1071 = -t707 / 0.2e1;
t1070 = t707 / 0.2e1;
t1069 = -t708 / 0.2e1;
t1068 = t708 / 0.2e1;
t60 = t1116 * t1149 + 0.2e1 * t1137 + t1138;
t1065 = t60 * qJD(3);
t1063 = pkin(5) * t707;
t1061 = t575 * pkin(5);
t1059 = t576 * pkin(5);
t1058 = t631 * pkin(5);
t1057 = t703 * pkin(3);
t1056 = t705 * pkin(5);
t1054 = qJD(5) / 0.2e1;
t46 = (t236 / 0.2e1 + t1107) * t634 - (t1105 + t1168) * t768;
t1052 = t46 * qJD(5);
t61 = t1097 * t1119 - t1137 + 0.2e1 * t1138;
t91 = -t1116 * t1118 + t1117 * t1119;
t1051 = t61 * qJD(5) + t91 * qJD(6);
t94 = (t1149 + t1101) * t634 - (t1103 + t1148) * t768;
t1050 = t94 * qJD(5);
t1049 = pkin(5) * qJD(5);
t1018 = t631 * t575;
t1019 = t628 * t573;
t162 = -t1018 - t1019;
t1047 = t162 * qJD(3);
t1046 = qJ(3) * t709;
t34 = -t111 * t1118 + t1119 * t112;
t1044 = qJD(1) * t34;
t56 = t1061 * t1119 + t115 * t709 + t1145;
t1043 = qJD(1) * t56;
t57 = -t1061 * t1118 + t116 * t709 + t1144;
t1042 = qJD(1) * t57;
t65 = -t111 * t709 + t1144;
t1041 = qJD(1) * t65;
t66 = -t112 * t709 + t1145;
t1040 = qJD(1) * t66;
t414 = t555 + t1002;
t418 = t884 - t1003;
t771 = t1081 * t414 + t1083 * t418;
t601 = t768 * t709;
t773 = t1089 * t1116 + t1150 * t601;
t97 = -t771 + t773;
t1039 = qJD(1) * t97;
t659 = t707 * pkin(2) - t1046;
t463 = (-pkin(8) * t709 - t659) * t704 + (-pkin(4) + t891) * t707;
t430 = t706 * t463;
t561 = -pkin(7) * t1012 + t703 * t659;
t518 = t561 + t694;
t495 = pkin(8) * t1015 + t518;
t466 = t708 * t495;
t302 = t466 + t430;
t222 = -pkin(9) * t576 + t302;
t1007 = t705 * t222;
t431 = t708 * t463;
t997 = t706 * t495;
t301 = t431 - t997;
t209 = -pkin(9) * t577 - t1063 + t301;
t889 = t1064 * t209;
t113 = t889 - t1007;
t1010 = t705 * t209;
t886 = t1064 * t222;
t114 = t886 + t1010;
t10 = t111 * t418 - t1118 * t113 - t1119 * t114 - t112 * t414;
t1037 = t10 * qJD(1);
t515 = t709 * t830 + t678;
t400 = t515 + t1059;
t14 = -t111 * t113 + t112 * t114 + t399 * t400;
t1034 = t14 * qJD(1);
t15 = t1061 * t399 - t111 * t115 + t112 * t116;
t1033 = t15 * qJD(1);
t23 = t111 * t707 + t1119 * t400 + t113 * t709 + t399 * t414;
t1032 = t23 * qJD(1);
t24 = -t1118 * t400 - t112 * t707 + t114 * t709 - t399 * t418;
t1031 = t24 * qJD(1);
t33 = t1012 * t399 - t111 * t600 + t112 * t601;
t1030 = t33 * qJD(1);
t1022 = t516 * t709;
t293 = t715 + t1127;
t53 = t293 * t577 - t294 * t576 - t301 * t575 - t302 * t573;
t1021 = t53 * qJD(1);
t55 = -t293 * t301 + t294 * t302 + t514 * t515;
t1020 = t55 * qJD(1);
t1017 = t631 * t709;
t1014 = t704 * t659;
t1013 = t704 * t701;
t78 = t293 * t707 + t301 * t709 + t514 * t576 + t515 * t573;
t990 = t78 * qJD(1);
t79 = -t294 * t707 + t302 * t709 - t514 * t577 - t515 * t575;
t989 = t79 * qJD(1);
t123 = (t1151 + t1094) * t634 - (t1150 + t1097) * t768;
t988 = t123 * qJD(5);
t143 = t1101 * t1117 + t1116 * t1148;
t987 = t143 * qJD(5) - t1143;
t144 = -t1094 * t1119 + t1162;
t986 = t144 * qJD(5) + t1143;
t406 = t539 + t538;
t474 = t588 + t587;
t983 = t991 * t1075 + t994 * t1076;
t851 = -t994 / 0.2e1;
t854 = -t1011 / 0.2e1;
t982 = t703 * t851 + t708 * t854;
t117 = -t1012 * t514 + t293 * t994 + t294 * t991;
t981 = qJD(1) * t117;
t119 = -t293 * t575 + t294 * t573;
t980 = qJD(1) * t119;
t864 = t511 * t1067;
t723 = t1084 * t514 + t1087 * t575 + t864;
t829 = t431 / 0.2e1 - t997 / 0.2e1;
t129 = -t723 + t829;
t979 = qJD(1) * t129;
t865 = t510 * t1067;
t724 = t865 + t514 * t628 / 0.2e1 + t573 * t1087;
t850 = -t430 / 0.2e1 - t466 / 0.2e1;
t130 = t724 + t850;
t978 = qJD(1) * t130;
t166 = -t1118 * t600 - t1119 * t601;
t977 = qJD(1) * t166;
t223 = t1119 * t707 - t414 * t709;
t976 = qJD(1) * t223;
t224 = -t1118 * t707 + t418 * t709;
t975 = qJD(1) * t224;
t861 = -t1017 / 0.2e1;
t247 = 0.2e1 * t862 * t708 + (t861 + t576 / 0.2e1) * t706;
t974 = qJD(1) * t247;
t973 = qJD(1) * t1158;
t968 = qJD(1) * t255;
t890 = pkin(7) + t1057;
t567 = t707 * t890 - t677;
t295 = pkin(3) * t1071 + (pkin(7) * t1071 + t567 / 0.2e1) * t703 + (-t659 / 0.2e1 + t1046 / 0.2e1 + t638 * t1070) * t704;
t966 = qJD(1) * t295;
t305 = t1012 * t1119 + t600 * t709;
t965 = qJD(1) * t305;
t306 = -t1012 * t1118 + t601 * t709;
t964 = qJD(1) * t306;
t322 = -t1012 * t517 + t1016 * t516;
t963 = qJD(1) * t322;
t344 = (t551 * t704 + t552 * t703) * t707;
t961 = qJD(1) * t344;
t345 = t1012 * t567 + t1022;
t960 = qJD(1) * t345;
t373 = -t573 * t991 + t575 * t994;
t959 = qJD(1) * t373;
t958 = qJD(1) * t1119;
t419 = t573 * t707 - t576 * t709;
t956 = qJD(1) * t419;
t420 = -t575 * t707 + t577 * t709;
t955 = qJD(1) * t420;
t954 = qJD(1) * t573;
t951 = qJD(2) * t1117;
t949 = qJD(2) * t508;
t948 = (t699 - t700) * t943;
t947 = qJD(2) * t628;
t945 = qJD(2) * t703;
t942 = qJD(3) * t709;
t941 = qJD(4) * t703;
t940 = qJD(4) * t704;
t939 = qJD(4) * t709;
t938 = qJD(5) * t1116;
t937 = qJD(5) * t1117;
t936 = qJD(5) * t628;
t935 = qJD(5) * t709;
t934 = qJD(6) * t1117;
t932 = qJD(6) * t508;
t135 = -t1118 * t414 - t1119 * t418;
t931 = t135 * qJD(1);
t519 = t707 * t891 - t1014;
t568 = t709 * t890 - t678;
t147 = t516 * t518 + t517 * t519 + t567 * t568;
t930 = t147 * qJD(1);
t163 = -t517 * t1011 - t519 * t1012 + (t518 * t707 + t1022) * t703;
t929 = t163 * qJD(1);
t164 = -t293 * t709 + t514 * t573;
t928 = t164 * qJD(1);
t165 = -t294 * t709 + t514 * t575;
t927 = t165 * qJD(1);
t901 = pkin(7) * t1016;
t560 = t901 + t1014;
t176 = (t551 * t709 + t560 * t707) * t704 + (t552 * t709 + t561 * t707) * t703;
t926 = t176 * qJD(1);
t815 = t567 * t709 + t568 * t707;
t177 = -t516 * t707 + t518 * t709 + t704 * t815;
t925 = t177 * qJD(1);
t178 = -t517 * t707 + t519 * t709 + t703 * t815;
t924 = t178 * qJD(1);
t225 = pkin(7) ^ 2 * t707 * t709 + t551 * t560 + t552 * t561;
t923 = t225 * qJD(1);
t304 = -t573 * t577 - t575 * t576;
t922 = t304 * qJD(1);
t329 = -t551 * t707 + (t560 - 0.2e1 * t901) * t709;
t921 = t329 * qJD(1);
t330 = t561 * t709 + (-t552 + 0.2e1 * t902) * t707;
t920 = t330 * qJD(1);
t852 = -t996 / 0.2e1;
t455 = t1079 + (t852 + t1085) * t709;
t919 = t455 * qJD(1);
t456 = t1079 + (t852 + t1084) * t709;
t918 = t456 * qJD(1);
t457 = t862 + t983;
t917 = t457 * qJD(1);
t458 = t862 + t982;
t916 = t458 * qJD(1);
t736 = t998 / 0.2e1 + t883 / 0.2e1;
t467 = (t1080 + t736) * t709;
t446 = t467 * qJD(1);
t735 = t882 / 0.2e1 - t999 / 0.2e1;
t468 = (t1082 + t735) * t709;
t448 = t468 * qJD(1);
t503 = t1012 * t573 + t702 * t706;
t915 = t503 * qJD(1);
t504 = t1012 * t575 + t702 * t708;
t914 = t504 * qJD(1);
t913 = t573 * qJD(3);
t912 = t573 * qJD(5);
t911 = t575 * qJD(3);
t566 = t575 * qJD(5);
t625 = t671 * t701;
t910 = t625 * qJD(1);
t908 = t627 * qJD(2);
t630 = t702 * t704 - t1013;
t616 = t630 * qJD(1);
t618 = t631 * qJD(5);
t906 = t671 * qJD(2);
t905 = t679 * qJD(1);
t904 = t707 * qJD(1);
t900 = -t1064 / 0.2e1;
t899 = t1064 / 0.2e1;
t898 = pkin(1) * t904;
t897 = pkin(1) * t952;
t896 = pkin(7) * t943;
t895 = t1061 / 0.2e1;
t894 = t1058 / 0.2e1;
t893 = pkin(5) * t1067;
t877 = t573 * t952;
t876 = t575 * t952;
t875 = t704 * t904;
t872 = t704 * t945;
t870 = qJD(3) * t1012;
t869 = t628 * t618;
t675 = t703 * t942;
t868 = t703 * t940;
t867 = t704 * t952;
t866 = t707 * t943;
t682 = t709 * t904;
t863 = -t577 / 0.2e1;
t860 = t1119 * t1076;
t859 = t1118 * t1076;
t858 = t573 * t1076;
t857 = t575 * t1076;
t856 = -t1012 / 0.2e1;
t855 = t1012 / 0.2e1;
t532 = -t1004 / 0.2e1;
t585 = -t1000 / 0.2e1;
t847 = t1064 * qJD(5);
t846 = t1064 * qJD(6);
t843 = -qJD(6) - t952;
t842 = qJD(2) * t610 - qJD(3);
t841 = qJD(1) * t703 * t1013;
t840 = t703 * t871;
t839 = t707 * t675;
t838 = t699 * t682;
t837 = t700 * t682;
t836 = t703 * t682;
t835 = t704 * t682;
t834 = 0.2e1 * t532;
t833 = 0.2e1 * t585;
t405 = t532 + t1004 / 0.2e1;
t473 = t585 + t1000 / 0.2e1;
t832 = -t886 / 0.2e1;
t831 = t707 * t900;
t827 = -t566 - t876;
t823 = t709 * t840;
t822 = t703 * t835;
t821 = -t638 * t709 + t993;
t730 = (t1073 * t601 + t600 * t899) * pkin(5);
t12 = (-t116 / 0.2e1 - t111 / 0.2e1) * t634 - (t1109 - t115 / 0.2e1) * t768 + t730;
t820 = -t12 * qJD(1) + t46 * qJD(2);
t718 = t1094 * t111 + t1107 * t1118 + t1119 * t1168 + t112 * t1150;
t721 = (t1110 + (pkin(3) / 0.2e1 + pkin(4) / 0.2e1) * t703) * t709 + t1078;
t720 = -t1059 / 0.2e1 + t721;
t19 = t718 + t720;
t77 = t1116 * t236 + t1117 * t237;
t819 = -qJD(1) * t19 + qJD(2) * t77;
t11 = -t1118 * t1165 - t1119 * t1166;
t818 = t11 * qJD(1) + t94 * qJD(4);
t817 = t518 * t704 + t519 * t703;
t816 = -t560 * t703 + t561 * t704;
t814 = qJD(1) * t94 + qJD(2) * t123;
t125 = t853 + t722 / 0.2e1;
t47 = (t1060 / 0.2e1 + t214 / 0.2e1 + t743 / 0.2e1 + t744 / 0.2e1 + (-pkin(5) + t745 / 0.2e1 - t733 / 0.2e1) * t709) * t705;
t813 = qJD(1) * t47 - qJD(2) * t125;
t746 = -t195 / 0.2e1 + t1064 * t893;
t49 = t888 / 0.2e1 + t746;
t812 = t49 * qJD(1);
t133 = -t1129 - t1152;
t811 = qJD(1) * t133 + qJD(2) * t60;
t810 = qJD(2) * t61 + t1171;
t158 = -t1130 - t1153;
t809 = qJD(1) * t60 + qJD(2) * t158;
t808 = qJD(1) * t61 + t1170;
t807 = qJD(2) * t91 + t1171;
t171 = -t1058 * t1117 - t1116 * t508;
t728 = t1093 * t1118 + t1151 * t399 - t1164;
t737 = -t1007 / 0.2e1 + t889 / 0.2e1;
t26 = (t1085 * t1119 + t1091 * t1117 + t831) * pkin(5) + t728 + t737;
t806 = qJD(1) * t26 + qJD(2) * t171;
t172 = -t1058 * t1116 + t1117 * t508;
t726 = t1097 * t399 + t1155;
t738 = -t1010 / 0.2e1 + t832;
t25 = (t1070 * t705 + t1085 * t1118 + t1091 * t1116) * pkin(5) + t726 + t738;
t805 = qJD(1) * t25 + qJD(2) * t172;
t804 = qJD(1) * t91 + t1170;
t249 = -t510 * t631 + t511 * t628;
t717 = t293 * t1084 + t294 * t1086 - t573 * t511 / 0.2e1 + t510 * t1090;
t63 = t717 + t721;
t803 = -qJD(1) * t63 + qJD(2) * t249;
t303 = -t1113 - t1114;
t801 = qJD(1) * t303 + qJD(2) * t162;
t348 = -t1111 - t1112;
t800 = qJD(1) * t162 + qJD(2) * t348;
t777 = t1081 * t1119 + t1083 * t1118;
t168 = t856 + t777;
t776 = t1081 * t1117 + t1083 * t1116;
t211 = t1077 + t776;
t799 = qJD(1) * t168 + qJD(2) * t211;
t740 = t1073 * t1119 + t1118 * t899;
t169 = (t1090 + t740) * pkin(5);
t739 = t1073 * t1117 + t1116 * t899;
t226 = (t1084 + t739) * pkin(5);
t798 = qJD(1) * t169 + qJD(2) * t226;
t212 = -t1018 + t1019;
t343 = -t1113 + t1114;
t797 = qJD(1) * t343 + qJD(2) * t212;
t425 = -t1111 + t1112;
t796 = qJD(1) * t212 + qJD(2) * t425;
t242 = -t1005 + 0.2e1 * t540;
t307 = -t1001 + 0.2e1 * t589;
t795 = qJD(1) * t242 + qJD(2) * t307;
t244 = 0.2e1 * t538 + t834;
t309 = 0.2e1 * t587 + t833;
t794 = qJD(1) * t244 + qJD(2) * t309;
t317 = t1004 + 0.2e1 * t539;
t358 = t1000 + 0.2e1 * t588;
t792 = qJD(1) * t317 + qJD(2) * t358;
t774 = t1069 * t575 + t1072 * t573;
t354 = t856 + t774;
t772 = t1069 * t631 + t1072 * t628;
t437 = t1077 + t772;
t791 = qJD(1) * t354 + qJD(2) * t437;
t790 = qJD(1) * t405 + qJD(2) * t473;
t789 = qJD(1) * t1120 + qJD(2) * t1121;
t788 = t947 + t954;
t787 = t946 + t953;
t786 = -qJD(5) * t1118 - t1122;
t785 = t903 * t1119;
t782 = t704 * t851 + t1079;
t778 = t1081 * t114 + t1083 * t113;
t716 = t1076 * t399 - t1089 * t236 + t1105 * t601 + t508 * t855;
t22 = t716 + t778;
t767 = qJD(1) * t22 + t508 * t945;
t725 = t1094 * t399 + t1101 * t508 + t1164;
t35 = -t725 + t737;
t766 = qJD(1) * t35 - t1116 * t949;
t727 = t1150 * t399 - t1155;
t36 = -t727 + t738;
t765 = qJD(1) * t36 + t1117 * t949;
t742 = t1076 * t514 + t610 * t855;
t76 = (t864 - t301 / 0.2e1) * t708 + (t865 - t302 / 0.2e1) * t706 + t742;
t764 = qJD(1) * t76 + t610 * t945;
t763 = qJD(2) * t143 + t1118 * t958;
t762 = qJD(1) * t143 + t1116 * t951;
t761 = qJD(2) * t144 - t880;
t760 = qJD(1) * t144 - t874;
t185 = t860 + (t1075 * t1117 + t1082) * t707;
t755 = qJD(1) * t185 + t1117 * t945;
t188 = t859 + (t1075 * t1116 + t1081) * t707;
t754 = qJD(1) * t188 + t1116 * t945;
t363 = t858 + (t1075 * t628 + t1068) * t707;
t749 = qJD(1) * t363 + t628 * t945;
t366 = t857 + (t1075 * t631 + t1072) * t707;
t748 = qJD(1) * t366 + t631 * t945;
t580 = (-0.1e1 / 0.2e1 + t699 / 0.2e1 - t700 / 0.2e1) * t707;
t747 = -qJD(1) * t580 + t872;
t623 = t875 + t945;
t741 = t1073 * t114 + t113 * t899;
t642 = t700 * t701 + t702;
t734 = qJD(1) * t642 + t840;
t719 = t1107 * t115 + t1109 * t236 + t1166 * t1168;
t1 = (t1085 * t399 + t1091 * t508 + t741) * pkin(5) + t719;
t51 = t1058 * t508;
t732 = -t1 * qJD(1) + t51 * qJD(2) + t46 * qJD(4);
t711 = t1103 * t236 - t1107 * t1119 + t1150 * t1166 + t1151 * t1165;
t729 = (t1074 * t414 + t418 * t900) * pkin(5);
t3 = t729 - t711;
t731 = -t3 * qJD(1) + t123 * qJD(4);
t689 = -t944 / 0.2e1;
t688 = -t904 / 0.2e1;
t687 = t904 / 0.2e1;
t652 = t671 * qJD(3);
t651 = t709 * t870;
t650 = t700 * t866;
t649 = t699 * t866;
t639 = t707 * t868;
t635 = t640 * qJD(3);
t622 = t1054 * t707 + t682;
t621 = -0.2e1 * t822;
t620 = 0.2e1 * t822;
t617 = t630 * qJD(2);
t611 = t625 * qJD(3);
t608 = qJD(2) * t699 + t703 * t875;
t581 = t700 * t1070 + t1071 * t699 + t1071;
t571 = t623 * t709;
t570 = -t704 * t943 + t836;
t569 = t682 + (t1054 + qJD(6) / 0.2e1) * t707;
t557 = (t700 * t904 + t872) * t709;
t556 = (t699 * t904 - t872) * t709;
t523 = t703 * t944 - t616;
t470 = t1067 * t634 + t709 * t736;
t469 = -t1067 * t768 + t709 * t735;
t462 = t861 + t782;
t461 = t1017 / 0.2e1 + t782;
t460 = t863 + t983;
t459 = t863 + t982;
t450 = t470 * qJD(4);
t449 = t469 * qJD(4);
t447 = t468 * qJD(4);
t445 = t467 * qJD(4);
t436 = t1077 - t772;
t365 = t857 + t631 * t855 + t995 / 0.2e1;
t364 = t858 + t628 * t855 - t992 / 0.2e1;
t359 = 0.2e1 * t1088 + t1001;
t357 = -t606 + t833;
t353 = t856 - t774;
t342 = t892 + t775;
t318 = 0.2e1 * t1092 + t1005;
t316 = -t553 + t834;
t314 = -t634 * t903 - t446;
t313 = t768 * t903 - t448;
t310 = t473 + t474;
t296 = qJ(3) * t854 + t638 * t856 + t567 * t1077 - t1014 / 0.2e1 + (-t1062 / 0.2e1 - pkin(3) / 0.2e1) * t707;
t265 = t1023 / 0.2e1 + t517 * t1076 + t1078 + (t1110 + t1057 / 0.2e1) * t709;
t246 = t576 * t1072 + t577 * t1069 + (t1068 * t628 + t1072 * t631) * t709;
t245 = t405 + t406;
t227 = -t1058 / 0.2e1 + t739 * pkin(5);
t210 = t1077 - t776;
t205 = t212 * qJD(5);
t187 = t1070 * t634 + t1116 * t855 + t859;
t186 = -t1071 * t768 + t1117 * t855 + t860;
t170 = -t1061 / 0.2e1 + t740 * pkin(5);
t167 = t856 - t777;
t132 = t723 + t829;
t131 = -t724 + t850;
t128 = -t885 + 0.2e1 * t853;
t96 = t771 + t773;
t75 = t302 * t706 / 0.2e1 + t301 * t1068 + (t1069 * t511 + t1072 * t510) * t709 + t742;
t64 = -t717 + t721;
t50 = t1008 - t888 / 0.2e1 + t746;
t48 = t705 * t893 - t887 + t710 * t1074 - t1009 / 0.2e1;
t38 = t725 + t737;
t37 = t727 + t738;
t28 = t1118 * t894 + t1116 * t895 + t832 + (t1063 / 0.2e1 - t209 / 0.2e1) * t705 - t726;
t27 = pkin(5) * t831 + t1117 * t895 + t1119 * t894 - t728 + t737;
t21 = t716 - t778;
t20 = -t718 + t720;
t13 = t1080 * t1166 + t1082 * t1165 + t730;
t4 = t729 + t711;
t2 = pkin(5) * t741 + t399 * t894 + t508 * t895 - t719;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t866, t679 * qJD(2), 0, -t866, 0, 0, -pkin(1) * t944, -pkin(1) * t943, 0, 0, t650, -0.2e1 * t823, -t617, t649, t908, -t866, -qJD(2) * t329 + t651, t330 * qJD(2) - t839, -qJD(2) * t176 + t611, qJD(2) * t225 - qJD(3) * t344, t650, -t617, 0.2e1 * t823, -t866, -t908, t649, qJD(2) * t178 - t701 * t868 + t651, -t163 * qJD(2) + t1016 * t939 + t611, -t177 * qJD(2) + t642 * qJD(4) + t839, qJD(2) * t147 - qJD(3) * t322 - qJD(4) * t345 (qJD(2) * t577 - t912) * t575, qJD(2) * t304 + qJD(5) * t343, t420 * qJD(2) - t709 * t912 (qJD(2) * t576 + t566) * t573, t419 * qJD(2) - t566 * t709, -t866, t78 * qJD(2) + t503 * qJD(4) + t165 * qJD(5) + t709 * t911, -t79 * qJD(2) + t504 * qJD(4) - t164 * qJD(5) - t709 * t913, qJD(2) * t53 + qJD(3) * t303 - qJD(4) * t373, qJD(2) * t55 + qJD(3) * t119 - qJD(4) * t117 (qJD(2) * t418 - t785) * t1118, qJD(2) * t135 + t1161 * t903, t224 * qJD(2) - t709 * t785 (qJD(2) * t414 - t786) * t1119, t223 * qJD(2) + t709 * t786, -t866, t23 * qJD(2) + t305 * qJD(4) + t56 * qJD(5) + t66 * qJD(6) + t1118 * t942, -t24 * qJD(2) - t306 * qJD(4) - t57 * qJD(5) - t65 * qJD(6) - t1119 * t942, qJD(2) * t10 + qJD(3) * t133 + qJD(4) * t166 + qJD(5) * t11, qJD(2) * t14 + qJD(3) * t34 + qJD(4) * t33 + qJD(5) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t682, t905, t943, -t682, -t944, 0, -t896 - t898, pkin(7) * t944 - t897, 0, 0, t557, t621 - t948, t523, t556, t1123, -t682, -t921 + (t703 * t824 - t902) * qJD(2) + t675, t703 * t896 + t920 + (qJD(2) * t824 + t942) * t704, qJD(2) * t816 - t926, t923 + (-pkin(2) * t1055 + qJ(3) * t816) * qJD(2) + t342 * qJD(3), t557, t523, t620 + t948, -t682, -t1123, t556, t924 + (-t568 * t704 - t703 * t821) * qJD(2) + t675 + t581 * qJD(4), qJD(2) * t817 - t929, -t568 * t945 - t925 + t639 + (qJD(2) * t821 - t942) * t704, t930 + (qJ(3) * t817 + t568 * t638) * qJD(2) + t265 * qJD(3) + t296 * qJD(4), t577 * t787 + t1142, t922 + (-t576 * t631 - t577 * t628) * qJD(2) + t205, qJD(5) * t460 - t631 * t944 + t955, t576 * t788 - t1142, qJD(5) * t462 + t628 * t944 + t956, -t622, t990 + (t510 * t707 + t515 * t628 + t576 * t610) * qJD(2) + t461 * qJD(3) + t364 * qJD(4) + t132 * qJD(5), -t989 + (t511 * t707 + t515 * t631 + t577 * t610) * qJD(2) + t459 * qJD(3) + t365 * qJD(4) + t131 * qJD(5), t1021 + (-t301 * t631 - t302 * t628 + t510 * t577 - t511 * t576) * qJD(2) + t246 * qJD(4) + t1047, t1020 + (-t301 * t510 + t302 * t511 + t515 * t610) * qJD(2) + t64 * qJD(3) + t75 * qJD(4), t418 * t793 + t986, t931 + (-t1116 * t414 - t1117 * t418) * qJD(2) + t1051, -t1116 * t944 - t255 * t903 + t975 (t951 + t958) * t414 + t987, t1117 * t944 + t1158 * t903 + t976, -t569, t1032 + (t1117 * t400 - t236 * t707 + t414 * t508) * qJD(2) + t1157 * qJD(3) + t186 * qJD(4) + t27 * qJD(5) + t38 * qJD(6), -t1031 + (t1116 * t400 + t237 * t707 + t418 * t508) * qJD(2) - t1159 * qJD(3) + t187 * qJD(4) + t28 * qJD(5) + t37 * qJD(6), t1037 + (-t1116 * t113 - t1117 * t114 - t236 * t418 - t237 * t414) * qJD(2) + t96 * qJD(4) + t4 * qJD(5) + t1065, t1034 + (t113 * t236 + t114 * t237 + t400 * t508) * qJD(2) + t20 * qJD(3) + t21 * qJD(4) + t2 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t571, -t570, t910, qJD(2) * t342 - t961, 0, 0, 0, 0, 0, 0, t571, t910, t570, qJD(2) * t265 - t963, 0, 0, 0, 0, 0, 0, t461 * qJD(2) + t876, t459 * qJD(2) - t877, t801, qJD(2) * t64 + qJD(4) * t353 + t980, 0, 0, 0, 0, 0, 0, t1141 + t1156, t245 * qJD(5) + t406 * qJD(6) - t1169, t811, qJD(2) * t20 + qJD(4) * t167 + qJD(5) * t170 + t1044; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t581 - t841, t836, t734, qJD(2) * t296 - t960, 0, 0, 0, 0, 0, 0, qJD(2) * t364 + t915, qJD(2) * t365 + t914, qJD(2) * t246 - t959, qJD(2) * t75 + qJD(3) * t353 - t981, 0, 0, 0, 0, 0, 0, qJD(2) * t186 + t470 * t903 + t965, qJD(2) * t187 + t469 * t903 - t964, qJD(2) * t96 + t1050 + t977, t1030 + t21 * qJD(2) + t167 * qJD(3) + (-t600 * t768 + t601 * t634) * qJD(4) + t13 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1131, t797, t460 * qJD(2) - t1128, t1131, t462 * qJD(2) + t827, t689, qJD(2) * t132 - qJD(5) * t294 + t927, qJD(2) * t131 + qJD(5) * t293 - t928, 0, 0, t761, t810, t316 * qJD(6) - t1119 * t844 - t1174, t763, t318 * qJD(6) - t1118 * t844 + t1163, t689, qJD(2) * t27 + qJD(3) * t1120 + qJD(5) * t115 + qJD(6) * t48 + t1043 + t450, qJD(2) * t28 + qJD(3) * t245 - qJD(5) * t116 + qJD(6) * t50 - t1042 + t449, t4 * qJD(2) + (t1064 * t1119 - t1118 * t705) * t1049 + t818, t1033 + t2 * qJD(2) + t170 * qJD(3) + t13 * qJD(4) + (t1064 * t115 + t116 * t705) * t1049; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1132, t807, t316 * qJD(5) + t1119 * t843 - t1174, t1132, t318 * qJD(5) + t1118 * t843 + t1163, t689, qJD(2) * t38 + qJD(5) * t48 - qJD(6) * t112 + t1040 + t450, qJD(2) * t37 + qJD(3) * t406 + qJD(5) * t50 + qJD(6) * t111 - t1041 + t449, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t682, -t905, 0, t682, 0, 0, t898, t897, 0, 0, -t837, t620, t616, -t838, -t909, t682, t921, -t920, t926, -qJD(3) * t341 - t923, -t837, t616, t621, t682, t909, -t838, -qJD(4) * t580 - t924, -t704 * t939 + t929, t639 + t925, -qJD(3) * t264 - qJD(4) * t295 - t930, -t577 * t953 + t1142, t205 - t922, -qJD(5) * t457 - t955, -t576 * t954 - t1142, -qJD(5) * t456 - t956, t622, -qJD(3) * t455 + qJD(4) * t363 - qJD(5) * t129 - t990, -qJD(3) * t458 + qJD(4) * t366 - qJD(5) * t130 + t989, qJD(4) * t247 - t1021 + t1047, -qJD(3) * t63 + qJD(4) * t76 - t1020, -t418 * t957 + t986, -t931 + t1051, -t1159 * t903 - t975, -t414 * t958 + t987, -t1157 * t903 - t976, t569, -qJD(3) * t1158 + qJD(4) * t185 - qJD(5) * t26 - qJD(6) * t35 - t1032, -qJD(3) * t255 + qJD(4) * t188 - qJD(5) * t25 - qJD(6) * t36 + t1031, qJD(4) * t97 - qJD(5) * t3 - t1037 + t1065, -qJD(3) * t19 + qJD(4) * t22 - qJD(5) * t1 - t1034; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t652, t635, 0, 0, 0, 0, 0, 0, t868, t652, t699 * qJD(4), -t638 * t941 + t635, -t869, t425 * qJD(5), 0, t869, 0, 0, t610 * t618 + t628 * t941, -t610 * t936 + t631 * t941, qJD(3) * t348, qJD(3) * t249 + t610 * t941 (-t934 - t937) * t1116, t903 * t1160, 0 (t933 + t938) * t1117, 0, 0, -qJD(5) * t171 + t1116 * t932 + t1117 * t941, -qJD(5) * t172 + t1116 * t941 - t1117 * t932, qJD(3) * t158, qJD(3) * t77 + qJD(5) * t51 + t508 * t941; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t906, t1125, 0, 0, 0, 0, 0, 0, 0, t906, 0, t1124, 0, 0, 0, 0, 0, 0, -t919, -t916, t800, qJD(4) * t436 + t803, 0, 0, 0, 0, 0, 0, -t973 + t1140, qJD(5) * t310 + qJD(6) * t474 - t968, t809, qJD(4) * t210 + qJD(5) * t227 + t819; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t747, -t867, t608, -t638 * t945 - t966, 0, 0, 0, 0, 0, 0, t749, t748, t974, qJD(3) * t436 + t764, 0, 0, 0, 0, 0, 0, t755, t754, t988 + t1039, qJD(3) * t210 + t1052 + t767; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1133, t796, -t917 - t936, t1133, -t618 - t918, t687, -qJD(5) * t511 + t610 * t946 - t979, qJD(5) * t510 - t610 * t947 - t978, 0, 0, t760, t808, qJD(6) * t357 - t1172 - t937, t762, qJD(6) * t359 - t1173 - t938, t687, qJD(3) * t1121 - qJD(5) * t237 + qJD(6) * t128 - t806, qJD(3) * t310 - t1154 - t805 (t1064 * t1117 - t1116 * t705) * t1049 + t731, t227 * qJD(3) + (-t1064 * t237 + t236 * t705) * t1049 + t732; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1134, t804, qJD(5) * t357 - t1172 - t934, t1134, qJD(5) * t359 - t1173 - t933, t687, qJD(5) * t128 - qJD(6) * t237 - t766, qJD(3) * t474 - t1154 - t765, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t835, t836, -t910, qJD(2) * t341 + t961, 0, 0, 0, 0, 0, 0, -t835, -t910, -t836, qJD(2) * t264 - t707 * t940 + t963, 0, 0, 0, 0, 0, 0, t455 * qJD(2) + t827, t458 * qJD(2) + t1128, -t801, qJD(2) * t63 + qJD(4) * t354 - t980, 0, 0, 0, 0, 0, 0, -t242 * qJD(5) - t1122 + t1163 - t881, -t244 * qJD(5) + t317 * qJD(6) + t1174 + t879, -t811, qJD(2) * t19 + qJD(4) * t168 - qJD(5) * t169 - t1044; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t906, -t1125, 0, 0, 0, 0, 0, 0, 0, -t906, 0, -t941 - t1124, 0, 0, 0, 0, 0, 0, -t618 + t919, t916 + t936, -t800, qJD(4) * t437 - t803, 0, 0, 0, 0, 0, 0, -qJD(5) * t307 - t933 + t973, -qJD(5) * t309 + qJD(6) * t358 + t968, -t809, qJD(4) * t211 - qJD(5) * t226 - t819; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t623, 0, 0, 0, 0, 0, 0, 0, 0, 0, t791, 0, 0, 0, 0, 0, 0, 0, 0, 0, t799; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t787, t788, 0, 0, 0, 0, 0, 0, 0, 0, -t795, -t794, 0, -t798; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t793, t792, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2) * t580 + t841, -t570, -t734, qJD(2) * t295 + t870 + t960, 0, 0, 0, 0, 0, 0, -t363 * qJD(2) - t706 * t935 - t915, -t366 * qJD(2) - t708 * t935 - t914, -qJD(2) * t247 + t959, -qJD(2) * t76 - qJD(3) * t354 + t981, 0, 0, 0, 0, 0, 0, -qJD(2) * t185 - t467 * t903 - t965, -qJD(2) * t188 - t468 * t903 + t964, -qJD(2) * t97 + t1050 - t977, -qJD(2) * t22 - qJD(3) * t168 - qJD(5) * t12 - t1030; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t747, t867, -t608, t966 + (qJD(2) * t638 + qJD(3)) * t703, 0, 0, 0, 0, 0, 0, -t749, -t748, -t974, -qJD(3) * t437 - t764, 0, 0, 0, 0, 0, 0, -t755, -t754, t988 - t1039, -qJD(3) * t211 + t1052 - t767; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t623, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t791, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t799; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t844 * t706, -t844 * t708, 0, 0, 0, 0, 0, 0, 0, 0, t314, t313, t814 (-t1064 * t634 - t705 * t768) * t1049 + t820; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t314, t313, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1131, -t797, t457 * qJD(2) + t877, -t1131, t456 * qJD(2) + t876, t689, t129 * qJD(2) + t706 * t939 + t911 - t927, t130 * qJD(2) + t708 * t939 - t913 + t928, 0, 0, -t761, -t810, t405 * qJD(6) + t1169, -t763, qJD(6) * t1120 + t1156, t689, qJD(2) * t26 + qJD(3) * t242 + qJD(6) * t47 - t1043 + t445, qJD(2) * t25 + qJD(3) * t244 + qJD(6) * t49 + t1042 + t447, qJD(2) * t3 - t818, qJD(2) * t1 + qJD(3) * t169 + qJD(4) * t12 - t1033; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1133, -t796, t917, -t1133, t918, t688, -t631 * t842 + t979, t628 * t842 + t978, 0, 0, -t760, -t808, qJD(6) * t473 + t1172, -t762, qJD(6) * t1121 + t1173, t688, qJD(3) * t307 - qJD(6) * t125 + t806, qJD(3) * t309 + t805, -t731, qJD(3) * t226 - t732; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t787, -t788, 0, 0, 0, 0, 0, 0, 0, 0, t795, t794, 0, t798; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t706 * t952, t708 * t952, 0, 0, 0, 0, 0, 0, 0, 0, t446, t448, -t814, -t820; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(6) * t1056, -pkin(5) * t846, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t790, 0, t789, 0, -t1056 * t903 + t813 (-t847 - t846) * pkin(5) + t812, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1132, -t807, -t405 * qJD(5) + t1169, -t1132, -t1141 + t1156, t689, qJD(2) * t35 + qJD(3) * t1118 - qJD(5) * t47 - t1040 + t445, qJD(2) * t36 - qJD(3) * t317 - qJD(5) * t49 + t1041 + t447, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1134, -t804, -qJD(5) * t473 + t1172, -t1134, t1173 - t1140, t688, qJD(3) * t1116 + qJD(5) * t125 + t766, -qJD(3) * t358 + t765, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t793, -t792, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t446, t448, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t790, 0, -t789, 0, t1049 * t705 - t813, pkin(5) * t847 - t812, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t5;
