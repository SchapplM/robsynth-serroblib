% Calculate inertial parameters regressor of inverse dynamics base forces vector with Newton-Euler for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% tauB_reg [6x(7*10)]
%   inertial parameter regressor of inverse dynamics base forces vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 11:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauB_reg = S6PRRRRR3_invdynB_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_invdynB_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_invdynB_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRR3_invdynB_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR3_invdynB_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_invdynB_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauB_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 11:05:23
% EndTime: 2019-05-05 11:06:27
% DurationCPUTime: 50.72s
% Computational Cost: add. (447403->1018), mult. (887647->1656), div. (0->0), fcn. (669314->14), ass. (0->728)
t1008 = sin(qJ(4));
t1013 = cos(qJ(4));
t1009 = sin(qJ(3));
t1109 = qJD(2) * t1009;
t959 = -t1013 * qJD(3) + t1008 * t1109;
t1014 = cos(qJ(3));
t1108 = qJD(2) * t1014;
t1095 = qJD(3) * t1108;
t1106 = qJDD(2) * t1009;
t963 = t1095 + t1106;
t1102 = t959 * qJD(4) - t1008 * qJDD(3) - t1013 * t963;
t987 = -qJD(4) + t1108;
t941 = t959 * t987;
t873 = t941 + t1102;
t1003 = cos(pkin(6));
t1001 = sin(pkin(6));
t1144 = g(3) - qJDD(1);
t1093 = t1001 * t1144;
t1000 = sin(pkin(12));
t1002 = cos(pkin(12));
t1094 = g(1) * t1000 - t1002 * g(2);
t1165 = t1003 * t1094 - t1093;
t1007 = sin(qJ(5));
t1012 = cos(qJ(5));
t1045 = qJDD(3) * t1013 - t1008 * t963;
t960 = qJD(3) * t1008 + t1013 * t1109;
t1020 = qJD(4) * t960 - t1045;
t912 = t1007 * t960 + t1012 * t959;
t799 = -t912 * qJD(5) - t1007 * t1020 - t1012 * t1102;
t979 = -qJD(5) + t987;
t887 = t912 * t979;
t763 = t887 - t799;
t1164 = t887 + t799;
t1010 = sin(qJ(2));
t1015 = cos(qJ(2));
t972 = g(1) * t1002 + g(2) * t1000;
t893 = -t1010 * t972 - t1015 * t1165;
t894 = t1010 * t1165 - t1015 * t972;
t823 = t1010 * t893 + t1015 * t894;
t1163 = t823 * t1001;
t1006 = sin(qJ(6));
t1011 = cos(qJ(6));
t1091 = -t1007 * t1102 + t1012 * t1020;
t914 = -t1007 * t959 + t1012 * t960;
t798 = -qJD(5) * t914 - t1091;
t847 = t1006 * t914 + t1011 * t912;
t1021 = qJD(6) * t847 - t1006 * t798 - t1011 * t799;
t973 = -qJD(6) + t979;
t826 = t847 * t973;
t1162 = -t1021 + t826;
t871 = t941 - t1102;
t1161 = t1000 * t1144;
t1160 = t1002 * t1144;
t849 = -t1006 * t912 + t1011 * t914;
t1147 = t847 * t849;
t1104 = qJDD(2) * t1014;
t990 = qJD(3) * t1109;
t964 = -t990 + t1104;
t956 = -qJDD(4) + t964;
t950 = -qJDD(5) + t956;
t942 = -qJDD(6) + t950;
t1017 = -t942 - t1147;
t1159 = t1006 * t1017;
t1146 = t912 * t914;
t1019 = -t950 - t1146;
t1158 = t1007 * t1019;
t1145 = t959 * t960;
t1022 = -t956 - t1145;
t1157 = t1008 * t1022;
t953 = t1001 * t1094;
t1083 = t1003 * t1144 + t953;
t1040 = t1010 * t1083;
t1156 = t1011 * t1017;
t1155 = t1012 * t1019;
t1154 = t1013 * t1022;
t1039 = t1015 * t1083;
t921 = -t1000 * t1094 - t1002 * t972;
t874 = (qJD(4) + t987) * t960 - t1045;
t764 = (qJD(5) + t979) * t914 + t1091;
t1092 = t1006 * t799 - t1011 * t798;
t664 = (qJD(6) + t973) * t849 + t1092;
t920 = -t1000 * t972 + t1002 * t1094;
t845 = t847 ^ 2;
t846 = t849 ^ 2;
t1153 = t912 ^ 2;
t911 = t914 ^ 2;
t1152 = t959 ^ 2;
t955 = t960 ^ 2;
t969 = t973 ^ 2;
t978 = t979 ^ 2;
t985 = t987 ^ 2;
t1151 = qJD(3) ^ 2;
t1150 = pkin(3) * t1009;
t1149 = pkin(3) * t1014;
t1148 = pkin(7) * t1003;
t1016 = qJD(2) ^ 2;
t881 = -t1016 * pkin(2) + qJDD(2) * pkin(8) + t894;
t844 = -t1009 * t1083 + t1014 * t881;
t1088 = -pkin(9) * t1009 - t1149;
t961 = t1088 * qJD(2);
t807 = -pkin(3) * t1151 + qJDD(3) * pkin(9) + t1108 * t961 + t844;
t1085 = t963 + t1095;
t1086 = -t964 + t990;
t880 = -qJDD(2) * pkin(2) - t1016 * pkin(8) + t893;
t821 = pkin(3) * t1086 - pkin(9) * t1085 + t880;
t735 = t1008 * t807 - t1013 * t821;
t700 = t1022 * pkin(4) + pkin(10) * t873 - t735;
t736 = t1008 * t821 + t1013 * t807;
t935 = -pkin(4) * t987 - pkin(10) * t960;
t706 = -pkin(4) * t1152 - pkin(10) * t1020 + t987 * t935 + t736;
t622 = t1007 * t706 - t1012 * t700;
t574 = t1019 * pkin(5) + pkin(11) * t763 - t622;
t623 = t1007 * t700 + t1012 * t706;
t882 = -pkin(5) * t979 - pkin(11) * t914;
t585 = -pkin(5) * t1153 + pkin(11) * t798 + t882 * t979 + t623;
t525 = t1006 * t574 + t1011 * t585;
t996 = t1009 ^ 2;
t997 = t1014 ^ 2;
t1143 = t996 + t997;
t934 = t1014 * t1083;
t806 = -t1151 * pkin(9) + t934 - qJDD(3) * pkin(3) + (qJD(2) * t961 + t881) * t1009;
t739 = t1020 * pkin(4) - t1152 * pkin(10) + t935 * t960 + t806;
t646 = -t798 * pkin(5) - pkin(11) * t1153 + t882 * t914 + t739;
t1141 = t1006 * t646;
t751 = t942 - t1147;
t1140 = t1006 * t751;
t1139 = t1006 * t973;
t524 = t1006 * t585 - t1011 * t574;
t473 = t1006 * t525 - t1011 * t524;
t1138 = t1007 * t473;
t1137 = t1007 * t739;
t828 = t950 - t1146;
t1136 = t1007 * t828;
t1135 = t1007 * t979;
t552 = t1007 * t623 - t1012 * t622;
t1134 = t1008 * t552;
t1133 = t1008 * t806;
t890 = t956 - t1145;
t1132 = t1008 * t890;
t1131 = t1008 * t987;
t1130 = t1009 * t880;
t986 = t1009 * t1016 * t1014;
t975 = -t986 + qJDD(3);
t1129 = t1009 * t975;
t976 = qJDD(3) + t986;
t1128 = t1009 * t976;
t1127 = t1011 * t646;
t1126 = t1011 * t751;
t1125 = t1011 * t973;
t1124 = t1012 * t473;
t1123 = t1012 * t739;
t1122 = t1012 * t828;
t1121 = t1012 * t979;
t1120 = t1013 * t552;
t1119 = t1013 * t806;
t1118 = t1013 * t890;
t1117 = t1013 * t987;
t1116 = t1014 * t880;
t1115 = t1014 * t975;
t1114 = t1016 * t996;
t1110 = qJD(2) * qJD(3);
t1107 = qJDD(2) * t1001;
t1105 = qJDD(2) * t1010;
t1103 = qJDD(2) * t1015;
t1101 = t1009 * t1147;
t1100 = t1009 * t1146;
t1099 = t1009 * t1145;
t1098 = t1014 * t1147;
t1097 = t1014 * t1146;
t1096 = t1014 * t1145;
t474 = t1006 * t524 + t1011 * t525;
t553 = t1007 * t622 + t1012 * t623;
t843 = t1009 * t881 + t934;
t767 = t1009 * t843 + t1014 * t844;
t1090 = t1010 * t986;
t1089 = t1015 * t986;
t765 = t1009 * t844 - t1014 * t843;
t966 = t1143 * qJDD(2);
t994 = t997 * t1016;
t970 = t994 + t1114;
t918 = -t1010 * t970 + t1015 * t966;
t1087 = pkin(7) * t918 - t1010 * t765;
t1043 = t1015 * t1016 + t1105;
t945 = t1043 * t1003;
t967 = -t1010 * t1016 + t1103;
t1082 = t1000 * t967 + t1002 * t945;
t905 = t1000 * t945 - t1002 * t967;
t667 = t1008 * t736 - t1013 * t735;
t668 = t1008 * t735 + t1013 * t736;
t434 = t1007 * t474 + t1124;
t435 = t1012 * t474 - t1138;
t406 = -t1008 * t434 + t1013 * t435;
t402 = t1009 * t646 + t1014 * t406;
t405 = t1008 * t435 + t1013 * t434;
t1081 = t1010 * t402 - t1015 * t405;
t663 = t1021 + t826;
t580 = -t1006 * t664 + t1011 * t663;
t582 = -t1006 * t663 - t1011 * t664;
t520 = t1007 * t582 + t1012 * t580;
t522 = -t1007 * t580 + t1012 * t582;
t471 = -t1008 * t520 + t1013 * t522;
t723 = -t845 - t846;
t463 = t1009 * t723 + t1014 * t471;
t469 = t1008 * t522 + t1013 * t520;
t1080 = t1010 * t463 - t1015 * t469;
t659 = (qJD(6) - t973) * t849 + t1092;
t581 = -t1006 * t659 + t1011 * t1162;
t583 = -t1006 * t1162 - t1011 * t659;
t521 = t1007 * t583 + t1012 * t581;
t523 = -t1007 * t581 + t1012 * t583;
t472 = -t1008 * t521 + t1013 * t523;
t770 = -t846 + t845;
t467 = -t1009 * t770 + t1014 * t472;
t470 = -t1008 * t523 - t1013 * t521;
t1079 = t1010 * t467 + t1015 * t470;
t492 = t1013 * t553 - t1134;
t487 = t1009 * t739 + t1014 * t492;
t491 = t1008 * t553 + t1120;
t1078 = t1010 * t487 - t1015 * t491;
t695 = -qJD(6) * t849 - t1092;
t655 = t1011 * t695 - t1139 * t847;
t656 = -t1006 * t695 - t1125 * t847;
t575 = t1007 * t656 + t1012 * t655;
t577 = -t1007 * t655 + t1012 * t656;
t517 = -t1008 * t575 + t1013 * t577;
t508 = t1014 * t517 - t1101;
t515 = -t1008 * t577 - t1013 * t575;
t1077 = t1010 * t508 + t1015 * t515;
t657 = -t1006 * t1021 - t1125 * t849;
t658 = -t1011 * t1021 + t1139 * t849;
t576 = t1007 * t658 + t1012 * t657;
t578 = -t1007 * t657 + t1012 * t658;
t518 = -t1008 * t576 + t1013 * t578;
t509 = t1014 * t518 + t1101;
t516 = -t1008 * t578 - t1013 * t576;
t1076 = t1010 * t509 + t1015 * t516;
t758 = -t969 - t845;
t689 = t1006 * t758 + t1156;
t690 = t1011 * t758 - t1159;
t606 = t1007 * t690 + t1012 * t689;
t607 = -t1007 * t689 + t1012 * t690;
t542 = -t1008 * t606 + t1013 * t607;
t511 = t1009 * t659 + t1014 * t542;
t541 = t1008 * t607 + t1013 * t606;
t1075 = t1010 * t511 - t1015 * t541;
t812 = -t846 - t969;
t714 = t1011 * t812 + t1140;
t715 = -t1006 * t812 + t1126;
t624 = t1007 * t715 + t1012 * t714;
t625 = -t1007 * t714 + t1012 * t715;
t558 = -t1008 * t624 + t1013 * t625;
t528 = t1009 * t1162 + t1014 * t558;
t557 = t1008 * t625 + t1013 * t624;
t1074 = t1010 * t528 - t1015 * t557;
t825 = -t846 + t969;
t716 = t1011 * t825 + t1159;
t718 = -t1006 * t825 + t1156;
t626 = t1007 * t718 + t1012 * t716;
t628 = -t1007 * t716 + t1012 * t718;
t564 = -t1008 * t626 + t1013 * t628;
t533 = -t1009 * t663 + t1014 * t564;
t562 = -t1008 * t628 - t1013 * t626;
t1073 = t1010 * t533 + t1015 * t562;
t824 = t845 - t969;
t717 = t1006 * t824 - t1126;
t719 = t1011 * t824 + t1140;
t627 = t1007 * t719 + t1012 * t717;
t629 = -t1007 * t717 + t1012 * t719;
t565 = -t1008 * t627 + t1013 * t629;
t534 = -t1009 * t664 + t1014 * t565;
t563 = -t1008 * t629 - t1013 * t627;
t1072 = t1010 * t534 + t1015 * t563;
t691 = -t1007 * t764 + t1012 * t763;
t693 = -t1007 * t763 - t1012 * t764;
t610 = -t1008 * t691 + t1013 * t693;
t802 = -t911 - t1153;
t589 = t1009 * t802 + t1014 * t610;
t608 = t1008 * t693 + t1013 * t691;
t1071 = t1010 * t589 - t1015 * t608;
t743 = (t1006 * t847 + t1011 * t849) * t973;
t744 = (-t1006 * t849 + t1011 * t847) * t973;
t676 = t1007 * t744 + t1012 * t743;
t677 = -t1007 * t743 + t1012 * t744;
t600 = -t1008 * t676 + t1013 * t677;
t595 = -t1009 * t942 + t1014 * t600;
t599 = -t1008 * t677 - t1013 * t676;
t1070 = t1010 * t595 + t1015 * t599;
t759 = (qJD(5) - t979) * t914 + t1091;
t692 = -t1007 * t759 + t1012 * t1164;
t694 = -t1007 * t1164 - t1012 * t759;
t611 = -t1008 * t692 + t1013 * t694;
t860 = -t911 + t1153;
t598 = -t1009 * t860 + t1014 * t611;
t609 = -t1008 * t694 - t1013 * t692;
t1069 = t1010 * t598 + t1015 * t609;
t838 = -t978 - t1153;
t749 = t1007 * t838 + t1155;
t750 = t1012 * t838 - t1158;
t683 = -t1008 * t749 + t1013 * t750;
t633 = t1009 * t759 + t1014 * t683;
t682 = t1008 * t750 + t1013 * t749;
t1068 = t1010 * t633 - t1015 * t682;
t635 = t1009 * t806 + t1014 * t668;
t1067 = t1010 * t635 - t1015 * t667;
t875 = -t911 - t978;
t773 = t1012 * t875 + t1136;
t774 = -t1007 * t875 + t1122;
t702 = -t1008 * t773 + t1013 * t774;
t640 = t1009 * t1164 + t1014 * t702;
t701 = t1008 * t774 + t1013 * t773;
t1066 = t1010 * t640 - t1015 * t701;
t886 = -t911 + t978;
t777 = t1012 * t886 + t1158;
t779 = -t1007 * t886 + t1155;
t709 = -t1008 * t777 + t1013 * t779;
t644 = -t1009 * t763 + t1014 * t709;
t707 = -t1008 * t779 - t1013 * t777;
t1065 = t1010 * t644 + t1015 * t707;
t885 = -t978 + t1153;
t778 = t1007 * t885 - t1122;
t780 = t1012 * t885 + t1136;
t710 = -t1008 * t778 + t1013 * t780;
t645 = -t1009 * t764 + t1014 * t710;
t708 = -t1008 * t780 - t1013 * t778;
t1064 = t1010 * t645 + t1015 * t708;
t754 = t1012 * t798 - t1135 * t912;
t755 = -t1007 * t798 - t1121 * t912;
t686 = -t1008 * t754 + t1013 * t755;
t652 = t1014 * t686 - t1100;
t684 = -t1008 * t755 - t1013 * t754;
t1063 = t1010 * t652 + t1015 * t684;
t756 = t1007 * t799 - t1121 * t914;
t757 = t1012 * t799 + t1135 * t914;
t687 = -t1008 * t756 + t1013 * t757;
t653 = t1014 * t687 + t1100;
t685 = -t1008 * t757 - t1013 * t756;
t1062 = t1010 * t653 + t1015 * t685;
t814 = (t1007 * t912 + t1012 * t914) * t979;
t815 = (-t1007 * t914 + t1012 * t912) * t979;
t738 = -t1008 * t814 + t1013 * t815;
t729 = -t1009 * t950 + t1014 * t738;
t737 = -t1008 * t815 - t1013 * t814;
t1061 = t1010 * t729 + t1015 * t737;
t793 = -t1008 * t873 - t1013 * t874;
t889 = t955 + t1152;
t746 = -t1009 * t889 + t1014 * t793;
t791 = -t1008 * t874 + t1013 * t873;
t1060 = t1010 * t746 - t1015 * t791;
t1059 = t1010 * t767 - t1015 * t880;
t869 = (-qJD(4) + t987) * t960 + t1045;
t792 = -t1008 * t871 + t1013 * t869;
t922 = -t955 + t1152;
t769 = -t1009 * t922 + t1014 * t792;
t790 = -t1008 * t869 - t1013 * t871;
t1058 = t1010 * t769 + t1015 * t790;
t907 = -t985 - t1152;
t832 = t1013 * t907 - t1157;
t772 = -t1009 * t869 + t1014 * t832;
t831 = t1008 * t907 + t1154;
t1057 = t1010 * t772 - t1015 * t831;
t917 = -t955 - t985;
t842 = -t1008 * t917 + t1118;
t776 = t1009 * t871 + t1014 * t842;
t841 = t1013 * t917 + t1132;
t1056 = t1010 * t776 - t1015 * t841;
t940 = -t955 + t985;
t858 = -t1008 * t940 + t1154;
t783 = -t1009 * t873 + t1014 * t858;
t856 = -t1013 * t940 - t1157;
t1055 = t1010 * t783 + t1015 * t856;
t939 = -t985 + t1152;
t859 = t1013 * t939 + t1132;
t784 = -t1009 * t874 + t1014 * t859;
t857 = -t1008 * t939 + t1118;
t1054 = t1010 * t784 + t1015 * t857;
t866 = t1008 * t1020 - t1117 * t959;
t810 = t1014 * t866 - t1099;
t865 = t1013 * t1020 + t1131 * t959;
t1053 = t1010 * t810 + t1015 * t865;
t868 = -t1013 * t1102 + t1131 * t960;
t811 = t1014 * t868 + t1099;
t867 = t1008 * t1102 + t1117 * t960;
t1052 = t1010 * t811 + t1015 * t867;
t879 = (-t1008 * t960 + t1013 * t959) * t987;
t862 = -t1009 * t956 + t1014 * t879;
t878 = (-t1008 * t959 - t1013 * t960) * t987;
t1051 = t1010 * t862 + t1015 * t878;
t1050 = t1010 * t894 - t1015 * t893;
t962 = 0.2e1 * t1095 + t1106;
t965 = -0.2e1 * t990 + t1104;
t916 = -t1009 * t962 + t1014 * t965;
t971 = t994 - t1114;
t1049 = t1010 * t916 + t1015 * t971;
t984 = -t994 - t1151;
t930 = t1014 * t984 - t1128;
t1048 = t1010 * t930 + t1015 * t965;
t982 = -t1114 - t1151;
t932 = -t1009 * t982 - t1115;
t1047 = t1010 * t932 - t1015 * t962;
t1046 = t1010 * t966 + t1015 * t970;
t957 = t1143 * t1110;
t1044 = -qJDD(3) * t1015 + t1010 * t957;
t958 = t1014 * t976;
t981 = -t1114 + t1151;
t931 = -t1009 * t981 + t958;
t1042 = -t1009 * t1103 + t1010 * t931;
t983 = t994 - t1151;
t929 = t1014 * t983 - t1129;
t1041 = t1010 * t929 - t1014 * t1103;
t457 = -pkin(5) * t646 + pkin(11) * t474;
t400 = -pkin(4) * t646 + pkin(10) * t435 - pkin(11) * t1138 + t1012 * t457;
t403 = -pkin(10) * t434 - pkin(11) * t1124 - t1007 * t457;
t382 = -pkin(9) * t405 - t1008 * t400 + t1013 * t403;
t394 = -pkin(3) * t405 - pkin(4) * t434 - pkin(5) * t473;
t401 = t1009 * t406 - t1014 * t646;
t376 = -pkin(8) * t401 - t1009 * t394 + t1014 * t382;
t379 = -pkin(2) * t401 + pkin(3) * t646 - pkin(9) * t406 - t1008 * t403 - t1013 * t400;
t392 = t1010 * t405 + t1015 * t402;
t1038 = pkin(7) * t392 + t1010 * t376 + t1015 * t379;
t449 = -pkin(5) * t723 + pkin(11) * t582 + t474;
t452 = -pkin(11) * t580 - t473;
t416 = -pkin(4) * t723 + pkin(10) * t522 + t1007 * t452 + t1012 * t449;
t418 = -pkin(10) * t520 - t1007 * t449 + t1012 * t452;
t393 = -pkin(9) * t469 - t1008 * t416 + t1013 * t418;
t439 = -pkin(3) * t469 - pkin(4) * t520 - pkin(5) * t580;
t462 = t1009 * t471 - t1014 * t723;
t388 = -pkin(8) * t462 - t1009 * t439 + t1014 * t393;
t391 = -pkin(2) * t462 + pkin(3) * t723 - pkin(9) * t471 - t1008 * t418 - t1013 * t416;
t432 = t1010 * t469 + t1015 * t463;
t1037 = pkin(7) * t432 + t1010 * t388 + t1015 * t391;
t554 = -pkin(5) * t659 + pkin(11) * t690 - t1127;
t590 = -pkin(11) * t689 + t1141;
t477 = -pkin(4) * t659 + pkin(10) * t607 + t1007 * t590 + t1012 * t554;
t483 = -pkin(10) * t606 - t1007 * t554 + t1012 * t590;
t428 = -pkin(9) * t541 - t1008 * t477 + t1013 * t483;
t455 = -pkin(3) * t541 - pkin(4) * t606 - pkin(5) * t689 + t524;
t510 = t1009 * t542 - t1014 * t659;
t407 = -pkin(8) * t510 - t1009 * t455 + t1014 * t428;
t420 = -pkin(2) * t510 + pkin(3) * t659 - pkin(9) * t542 - t1008 * t483 - t1013 * t477;
t476 = t1010 * t541 + t1015 * t511;
t1036 = pkin(7) * t476 + t1010 * t407 + t1015 * t420;
t559 = -pkin(5) * t1162 + pkin(11) * t715 + t1141;
t604 = -pkin(11) * t714 + t1127;
t480 = -pkin(4) * t1162 + pkin(10) * t625 + t1007 * t604 + t1012 * t559;
t490 = -pkin(10) * t624 - t1007 * t559 + t1012 * t604;
t437 = -pkin(9) * t557 - t1008 * t480 + t1013 * t490;
t460 = -pkin(3) * t557 - pkin(4) * t624 - pkin(5) * t714 + t525;
t527 = t1009 * t558 - t1014 * t1162;
t411 = -pkin(8) * t527 - t1009 * t460 + t1014 * t437;
t421 = -pkin(2) * t527 + pkin(3) * t1162 - pkin(9) * t558 - t1008 * t490 - t1013 * t480;
t482 = t1010 * t557 + t1015 * t528;
t1035 = pkin(7) * t482 + t1010 * t411 + t1015 * t421;
t543 = -pkin(4) * t739 + pkin(10) * t553;
t447 = -pkin(9) * t491 - pkin(10) * t1120 - t1008 * t543;
t461 = -pkin(3) * t491 - pkin(4) * t552;
t486 = t1009 * t492 - t1014 * t739;
t412 = -pkin(8) * t486 - t1009 * t461 + t1014 * t447;
t424 = -pkin(2) * t486 + pkin(3) * t739 - pkin(9) * t492 + pkin(10) * t1134 - t1013 * t543;
t446 = t1010 * t491 + t1015 * t487;
t1034 = pkin(7) * t446 + t1010 * t412 + t1015 * t424;
t526 = -pkin(4) * t802 + pkin(10) * t693 + t553;
t529 = -pkin(10) * t691 - t552;
t456 = -pkin(9) * t608 - t1008 * t526 + t1013 * t529;
t568 = -pkin(3) * t608 - pkin(4) * t691;
t588 = t1009 * t610 - t1014 * t802;
t440 = -pkin(8) * t588 - t1009 * t568 + t1014 * t456;
t448 = -pkin(2) * t588 + pkin(3) * t802 - pkin(9) * t610 - t1008 * t529 - t1013 * t526;
t535 = t1010 * t608 + t1015 * t589;
t1033 = pkin(7) * t535 + t1010 * t440 + t1015 * t448;
t630 = -pkin(4) * t759 + pkin(10) * t750 - t1123;
t675 = -pkin(10) * t749 + t1137;
t545 = -pkin(9) * t682 - t1008 * t630 + t1013 * t675;
t560 = -pkin(3) * t682 - pkin(4) * t749 + t622;
t632 = t1009 * t683 - t1014 * t759;
t479 = -pkin(8) * t632 - t1009 * t560 + t1014 * t545;
t504 = -pkin(2) * t632 + pkin(3) * t759 - pkin(9) * t683 - t1008 * t675 - t1013 * t630;
t571 = t1010 * t682 + t1015 * t633;
t1032 = pkin(7) * t571 + t1010 * t479 + t1015 * t504;
t636 = -pkin(4) * t1164 + pkin(10) * t774 + t1137;
t688 = -pkin(10) * t773 + t1123;
t556 = -pkin(9) * t701 - t1008 * t636 + t1013 * t688;
t566 = -pkin(3) * t701 - pkin(4) * t773 + t623;
t639 = t1009 * t702 - t1014 * t1164;
t481 = -pkin(8) * t639 - t1009 * t566 + t1014 * t556;
t512 = -pkin(2) * t639 + pkin(3) * t1164 - pkin(9) * t702 - t1008 * t688 - t1013 * t636;
t591 = t1010 * t701 + t1015 * t640;
t1031 = pkin(7) * t591 + t1010 * t481 + t1015 * t512;
t634 = t1009 * t668 - t1014 * t806;
t544 = -pkin(8) * t634 + (-pkin(9) * t1014 + t1150) * t667;
t567 = -pkin(2) * t634 + pkin(3) * t806 - pkin(9) * t668;
t570 = t1010 * t667 + t1015 * t635;
t1030 = pkin(7) * t570 + t1010 * t544 + t1015 * t567;
t631 = -pkin(9) * t791 - t667;
t745 = t1009 * t793 + t1014 * t889;
t579 = -pkin(8) * t745 + t1014 * t631 + t1150 * t791;
t596 = -pkin(2) * t745 - pkin(3) * t889 - pkin(9) * t793 - t668;
t699 = t1010 * t791 + t1015 * t746;
t1029 = pkin(7) * t699 + t1010 * t579 + t1015 * t596;
t705 = -pkin(3) * t831 + t735;
t740 = -pkin(9) * t831 + t1133;
t771 = t1009 * t832 + t1014 * t869;
t615 = -pkin(8) * t771 - t1009 * t705 + t1014 * t740;
t669 = -pkin(2) * t771 - pkin(3) * t869 - pkin(9) * t832 + t1119;
t722 = t1010 * t831 + t1015 * t772;
t1028 = pkin(7) * t722 + t1010 * t615 + t1015 * t669;
t713 = -pkin(3) * t841 + t736;
t742 = -pkin(9) * t841 + t1119;
t775 = t1009 * t842 - t1014 * t871;
t618 = -pkin(8) * t775 - t1009 * t713 + t1014 * t742;
t670 = -pkin(2) * t775 + pkin(3) * t871 - pkin(9) * t842 - t1133;
t724 = t1010 * t841 + t1015 * t776;
t1027 = pkin(7) * t724 + t1010 * t618 + t1015 * t670;
t926 = t1009 * t984 + t958;
t803 = -pkin(2) * t926 + t843;
t839 = -pkin(8) * t926 + t1130;
t883 = -t1010 * t965 + t1015 * t930;
t1026 = pkin(7) * t883 + t1010 * t839 + t1015 * t803;
t928 = t1014 * t982 - t1129;
t804 = -pkin(2) * t928 + t844;
t840 = -pkin(8) * t928 + t1116;
t884 = t1010 * t962 + t1015 * t932;
t1025 = pkin(7) * t884 + t1010 * t840 + t1015 * t804;
t936 = -t1009 * t964 - t1110 * t997;
t1024 = t1010 * t936 - t1089;
t937 = t1014 * t963 - t1110 * t996;
t1023 = t1010 * t937 + t1089;
t730 = t1010 * t880 + t1015 * t767;
t1018 = pkin(7) * t730 + (-pkin(2) * t1015 - pkin(8) * t1010) * t765;
t946 = t967 * t1003;
t944 = t967 * t1001;
t943 = t1043 * t1001;
t938 = qJDD(3) * t1010 + t1015 * t957;
t927 = t1014 * t981 + t1128;
t925 = t1009 * t983 + t1115;
t924 = t1085 * t1009;
t923 = t1086 * t1014;
t919 = t1044 * t1003;
t915 = t1009 * t965 + t1014 * t962;
t910 = t1046 * t1003;
t909 = t1046 * t1001;
t906 = -t1000 * t946 - t1002 * t1043;
t904 = -t1000 * t1043 + t1002 * t946;
t901 = t1015 * t937 - t1090;
t900 = t1015 * t936 + t1090;
t899 = t1009 * t1105 + t1015 * t931;
t898 = t1010 * t1104 + t1015 * t929;
t877 = -t1010 * t971 + t1015 * t916;
t864 = -t1039 + (t1001 * t943 + t1003 * t945) * pkin(7);
t863 = -t1040 + (-t1001 * t944 - t1003 * t946) * pkin(7);
t861 = t1009 * t879 + t1014 * t956;
t855 = -t1000 * t910 + t1002 * t918;
t854 = t1000 * t918 + t1002 * t910;
t853 = -t1001 * t924 + t1003 * t1023;
t852 = t1001 * t923 + t1003 * t1024;
t851 = -t1001 * t927 + t1003 * t1042;
t850 = -t1001 * t925 + t1003 * t1041;
t836 = -t1001 * t928 + t1003 * t1047;
t835 = -t1001 * t926 + t1003 * t1048;
t834 = t1001 * t1047 + t1003 * t928;
t833 = t1001 * t1048 + t1003 * t926;
t820 = -t1001 * t915 + t1003 * t1049;
t819 = pkin(2) * t965 + pkin(8) * t930 - t1116;
t818 = -pkin(2) * t962 + pkin(8) * t932 + t1130;
t813 = t823 * t1003;
t809 = t1009 * t868 - t1096;
t808 = t1009 * t866 + t1096;
t801 = -pkin(1) * t944 + t1001 * t893 + t1003 * t1039 - t1043 * t1148;
t800 = pkin(1) * t943 + t1001 * t894 - t1003 * t1040 - t1148 * t967;
t795 = t1001 * t953 + (t1093 + t1050) * t1003;
t794 = t1001 * t1050 - t1003 * t1083;
t789 = -t1010 * t878 + t1015 * t862;
t788 = -t1000 * t836 + t1002 * t884;
t787 = -t1000 * t835 + t1002 * t883;
t786 = t1000 * t884 + t1002 * t836;
t785 = t1000 * t883 + t1002 * t835;
t782 = t1009 * t859 + t1014 * t874;
t781 = t1009 * t858 + t1014 * t873;
t768 = t1009 * t792 + t1014 * t922;
t748 = -t1010 * t867 + t1015 * t811;
t747 = -t1010 * t865 + t1015 * t810;
t741 = pkin(2) * t970 + pkin(8) * t966 + t767;
t734 = -t1001 * t861 + t1003 * t1051;
t733 = -pkin(2) * t880 + pkin(8) * t767;
t732 = -t1010 * t857 + t1015 * t784;
t731 = -t1010 * t856 + t1015 * t783;
t728 = t1009 * t738 + t1014 * t950;
t727 = -pkin(1) * t794 + t1148 * t823;
t726 = -t1000 * t795 + t1002 * t823;
t725 = t1000 * t823 + t1002 * t795;
t721 = -t1015 * t765 + (-t1001 * t909 - t1003 * t910) * pkin(7);
t720 = (-t1001 * t794 - t1003 * t795) * pkin(7);
t712 = -t1001 * t809 + t1003 * t1052;
t711 = -t1001 * t808 + t1003 * t1053;
t703 = -t1010 * t790 + t1015 * t769;
t681 = -t1010 * t804 + t1015 * t840 + (-t1001 * t834 - t1003 * t836) * pkin(7);
t680 = -t1010 * t803 + t1015 * t839 + (-t1001 * t833 - t1003 * t835) * pkin(7);
t679 = -t1001 * t782 + t1003 * t1054;
t678 = -t1001 * t781 + t1003 * t1055;
t674 = -t1001 * t775 + t1003 * t1056;
t673 = t1001 * t1056 + t1003 * t775;
t672 = -t1001 * t765 + t1003 * t1059;
t671 = t1001 * t1059 + t1003 * t765;
t666 = -t1001 * t771 + t1003 * t1057;
t665 = t1001 * t1057 + t1003 * t771;
t654 = -t1010 * t737 + t1015 * t729;
t651 = t1009 * t687 - t1097;
t650 = t1009 * t686 + t1097;
t649 = -pkin(1) * t834 - t1001 * t818 + t1003 * t1025;
t648 = -pkin(1) * t833 - t1001 * t819 + t1003 * t1026;
t647 = -pkin(1) * t909 - t1001 * t741 + t1003 * t1087;
t643 = t1009 * t710 + t1014 * t764;
t642 = t1009 * t709 + t1014 * t763;
t641 = -t1001 * t768 + t1003 * t1058;
t638 = -t1001 * t745 + t1003 * t1060;
t637 = t1001 * t1060 + t1003 * t745;
t620 = -t1000 * t672 + t1002 * t730;
t619 = t1000 * t730 + t1002 * t672;
t617 = -t1000 * t674 + t1002 * t724;
t616 = t1000 * t724 + t1002 * t674;
t614 = -t1001 * t728 + t1003 * t1061;
t613 = -t1000 * t666 + t1002 * t722;
t612 = t1000 * t722 + t1002 * t666;
t605 = -pkin(2) * t841 + pkin(8) * t776 + t1009 * t742 + t1014 * t713;
t603 = -pkin(2) * t831 + pkin(8) * t772 + t1009 * t740 + t1014 * t705;
t602 = -t1010 * t708 + t1015 * t645;
t601 = -t1010 * t707 + t1015 * t644;
t597 = t1009 * t611 + t1014 * t860;
t594 = t1009 * t600 + t1014 * t942;
t593 = -t1010 * t685 + t1015 * t653;
t592 = -t1010 * t684 + t1015 * t652;
t587 = -t1000 * t638 + t1002 * t699;
t586 = t1000 * t699 + t1002 * t638;
t569 = pkin(8) * t746 + t1009 * t631 + (-pkin(2) - t1149) * t791;
t561 = (pkin(2) * t1010 - pkin(8) * t1015) * t765 + (-t1001 * t671 - t1003 * t672) * pkin(7);
t555 = -pkin(1) * t671 - t1001 * t733 + t1003 * t1018;
t551 = -t1001 * t643 + t1003 * t1064;
t550 = -t1001 * t642 + t1003 * t1065;
t549 = -t1001 * t651 + t1003 * t1062;
t548 = -t1001 * t650 + t1003 * t1063;
t547 = -t1001 * t639 + t1003 * t1066;
t546 = t1001 * t1066 + t1003 * t639;
t540 = -t1001 * t632 + t1003 * t1068;
t539 = t1001 * t1068 + t1003 * t632;
t538 = -t1010 * t609 + t1015 * t598;
t537 = -t1001 * t634 + t1003 * t1067;
t536 = t1001 * t1067 + t1003 * t634;
t532 = t1009 * t565 + t1014 * t664;
t531 = t1009 * t564 + t1014 * t663;
t530 = -t1010 * t599 + t1015 * t595;
t514 = pkin(8) * t635 + (-pkin(2) + t1088) * t667;
t513 = -t1010 * t670 + t1015 * t618 + (-t1001 * t673 - t1003 * t674) * pkin(7);
t507 = t1009 * t518 - t1098;
t506 = t1009 * t517 + t1098;
t505 = -t1010 * t669 + t1015 * t615 + (-t1001 * t665 - t1003 * t666) * pkin(7);
t503 = -t1000 * t547 + t1002 * t591;
t502 = t1000 * t591 + t1002 * t547;
t501 = -t1001 * t597 + t1003 * t1069;
t500 = -t1001 * t588 + t1003 * t1071;
t499 = t1001 * t1071 + t1003 * t588;
t498 = -t1000 * t540 + t1002 * t571;
t497 = t1000 * t571 + t1002 * t540;
t496 = -t1001 * t594 + t1003 * t1070;
t495 = -pkin(1) * t673 - t1001 * t605 + t1003 * t1027;
t494 = -t1000 * t537 + t1002 * t570;
t493 = t1000 * t570 + t1002 * t537;
t489 = -pkin(1) * t665 - t1001 * t603 + t1003 * t1028;
t488 = -t1010 * t596 + t1015 * t579 + (-t1001 * t637 - t1003 * t638) * pkin(7);
t485 = -t1010 * t563 + t1015 * t534;
t484 = -t1010 * t562 + t1015 * t533;
t478 = -pkin(2) * t701 + pkin(8) * t640 + t1009 * t556 + t1014 * t566;
t475 = -pkin(2) * t682 + pkin(8) * t633 + t1009 * t545 + t1014 * t560;
t468 = -pkin(1) * t637 - t1001 * t569 + t1003 * t1029;
t466 = t1009 * t472 + t1014 * t770;
t465 = -t1010 * t516 + t1015 * t509;
t464 = -t1010 * t515 + t1015 * t508;
t459 = -t1000 * t500 + t1002 * t535;
t458 = t1000 * t535 + t1002 * t500;
t454 = -t1001 * t532 + t1003 * t1072;
t453 = -t1001 * t531 + t1003 * t1073;
t451 = -t1001 * t527 + t1003 * t1074;
t450 = t1001 * t1074 + t1003 * t527;
t445 = -t1001 * t510 + t1003 * t1075;
t444 = t1001 * t1075 + t1003 * t510;
t443 = -t1010 * t567 + t1015 * t544 + (-t1001 * t536 - t1003 * t537) * pkin(7);
t442 = -t1001 * t507 + t1003 * t1076;
t441 = -t1001 * t506 + t1003 * t1077;
t438 = -pkin(2) * t608 + pkin(8) * t589 + t1009 * t456 + t1014 * t568;
t436 = -pkin(1) * t536 - t1001 * t514 + t1003 * t1030;
t433 = -t1010 * t470 + t1015 * t467;
t431 = -t1000 * t451 + t1002 * t482;
t430 = t1000 * t482 + t1002 * t451;
t429 = -t1010 * t512 + t1015 * t481 + (-t1001 * t546 - t1003 * t547) * pkin(7);
t427 = -t1001 * t486 + t1003 * t1078;
t426 = t1001 * t1078 + t1003 * t486;
t425 = -t1010 * t504 + t1015 * t479 + (-t1001 * t539 - t1003 * t540) * pkin(7);
t423 = -t1000 * t445 + t1002 * t476;
t422 = t1000 * t476 + t1002 * t445;
t419 = -pkin(1) * t546 - t1001 * t478 + t1003 * t1031;
t417 = -t1001 * t466 + t1003 * t1079;
t415 = -pkin(1) * t539 - t1001 * t475 + t1003 * t1032;
t414 = -t1001 * t462 + t1003 * t1080;
t413 = t1001 * t1080 + t1003 * t462;
t410 = -t1000 * t427 + t1002 * t446;
t409 = t1000 * t446 + t1002 * t427;
t408 = -pkin(2) * t557 + pkin(8) * t528 + t1009 * t437 + t1014 * t460;
t404 = -pkin(2) * t491 + pkin(8) * t487 + t1009 * t447 + t1014 * t461;
t399 = -pkin(2) * t541 + pkin(8) * t511 + t1009 * t428 + t1014 * t455;
t398 = -t1010 * t448 + t1015 * t440 + (-t1001 * t499 - t1003 * t500) * pkin(7);
t397 = -t1000 * t414 + t1002 * t432;
t396 = t1000 * t432 + t1002 * t414;
t395 = -pkin(1) * t499 - t1001 * t438 + t1003 * t1033;
t390 = -t1010 * t421 + t1015 * t411 + (-t1001 * t450 - t1003 * t451) * pkin(7);
t389 = -t1010 * t420 + t1015 * t407 + (-t1001 * t444 - t1003 * t445) * pkin(7);
t387 = -t1010 * t424 + t1015 * t412 + (-t1001 * t426 - t1003 * t427) * pkin(7);
t386 = -pkin(2) * t469 + pkin(8) * t463 + t1009 * t393 + t1014 * t439;
t385 = -pkin(1) * t450 - t1001 * t408 + t1003 * t1035;
t384 = -t1001 * t401 + t1003 * t1081;
t383 = t1001 * t1081 + t1003 * t401;
t381 = -pkin(1) * t444 - t1001 * t399 + t1003 * t1036;
t380 = -pkin(1) * t426 - t1001 * t404 + t1003 * t1034;
t378 = -t1000 * t384 + t1002 * t392;
t377 = t1000 * t392 + t1002 * t384;
t375 = -t1010 * t391 + t1015 * t388 + (-t1001 * t413 - t1003 * t414) * pkin(7);
t374 = -pkin(2) * t405 + pkin(8) * t402 + t1009 * t382 + t1014 * t394;
t373 = -pkin(1) * t413 - t1001 * t386 + t1003 * t1037;
t372 = -t1010 * t379 + t1015 * t376 + (-t1001 * t383 - t1003 * t384) * pkin(7);
t371 = -pkin(1) * t383 - t1001 * t374 + t1003 * t1038;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, t921, 0, 0, 0, 0, 0, 0, t906, t905, 0, t726, 0, 0, 0, 0, 0, 0, t787, t788, t855, t620, 0, 0, 0, 0, 0, 0, t613, t617, t587, t494, 0, 0, 0, 0, 0, 0, t498, t503, t459, t410, 0, 0, 0, 0, 0, 0, t423, t431, t397, t378; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, t920, 0, 0, 0, 0, 0, 0, t904, -t1082, 0, t725, 0, 0, 0, 0, 0, 0, t785, t786, t854, t619, 0, 0, 0, 0, 0, 0, t612, t616, t586, t493, 0, 0, 0, 0, 0, 0, t497, t502, t458, t409, 0, 0, 0, 0, 0, 0, t422, t430, t396, t377; 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1144, 0, 0, 0, 0, 0, 0, t944, -t943, 0, t794, 0, 0, 0, 0, 0, 0, t833, t834, t909, t671, 0, 0, 0, 0, 0, 0, t665, t673, t637, t536, 0, 0, 0, 0, 0, 0, t539, t546, t499, t426, 0, 0, 0, 0, 0, 0, t444, t450, t413, t383; 0, 0, 0, 0, 0, 0, 0, -g(3), g(2), 0, 0, 0, 0, 0, 0, 0, -t1161, -t1160, -t920, -qJ(1) * t920, 0, 0, -t905, 0, t906, t1000 * t1107, -qJ(1) * t904 - t1000 * t801 + t1002 * t863, qJ(1) * t1082 - t1000 * t800 + t1002 * t864, -t1000 * t813 - t1002 * t1050, -qJ(1) * t725 - t1000 * t727 + t1002 * t720, -t1000 * t853 + t1002 * t901, -t1000 * t820 + t1002 * t877, -t1000 * t851 + t1002 * t899, -t1000 * t852 + t1002 * t900, -t1000 * t850 + t1002 * t898, -t1000 * t919 + t1002 * t938, -qJ(1) * t785 - t1000 * t648 + t1002 * t680, -qJ(1) * t786 - t1000 * t649 + t1002 * t681, -qJ(1) * t854 - t1000 * t647 + t1002 * t721, -qJ(1) * t619 - t1000 * t555 + t1002 * t561, -t1000 * t712 + t1002 * t748, -t1000 * t641 + t1002 * t703, -t1000 * t678 + t1002 * t731, -t1000 * t711 + t1002 * t747, -t1000 * t679 + t1002 * t732, -t1000 * t734 + t1002 * t789, -qJ(1) * t612 - t1000 * t489 + t1002 * t505, -qJ(1) * t616 - t1000 * t495 + t1002 * t513, -qJ(1) * t586 - t1000 * t468 + t1002 * t488, -qJ(1) * t493 - t1000 * t436 + t1002 * t443, -t1000 * t549 + t1002 * t593, -t1000 * t501 + t1002 * t538, -t1000 * t550 + t1002 * t601, -t1000 * t548 + t1002 * t592, -t1000 * t551 + t1002 * t602, -t1000 * t614 + t1002 * t654, -qJ(1) * t497 - t1000 * t415 + t1002 * t425, -qJ(1) * t502 - t1000 * t419 + t1002 * t429, -qJ(1) * t458 - t1000 * t395 + t1002 * t398, -qJ(1) * t409 - t1000 * t380 + t1002 * t387, -t1000 * t442 + t1002 * t465, -t1000 * t417 + t1002 * t433, -t1000 * t453 + t1002 * t484, -t1000 * t441 + t1002 * t464, -t1000 * t454 + t1002 * t485, -t1000 * t496 + t1002 * t530, -qJ(1) * t422 - t1000 * t381 + t1002 * t389, -qJ(1) * t430 - t1000 * t385 + t1002 * t390, -qJ(1) * t396 - t1000 * t373 + t1002 * t375, -qJ(1) * t377 - t1000 * t371 + t1002 * t372; 0, 0, 0, 0, 0, 0, g(3), 0, -g(1), 0, 0, 0, 0, 0, 0, 0, t1160, -t1161, t921, qJ(1) * t921, 0, 0, t1082, 0, t904, -t1002 * t1107, qJ(1) * t906 + t1000 * t863 + t1002 * t801, qJ(1) * t905 + t1000 * t864 + t1002 * t800, -t1000 * t1050 + t1002 * t813, qJ(1) * t726 + t1000 * t720 + t1002 * t727, t1000 * t901 + t1002 * t853, t1000 * t877 + t1002 * t820, t1000 * t899 + t1002 * t851, t1000 * t900 + t1002 * t852, t1000 * t898 + t1002 * t850, t1000 * t938 + t1002 * t919, qJ(1) * t787 + t1000 * t680 + t1002 * t648, qJ(1) * t788 + t1000 * t681 + t1002 * t649, qJ(1) * t855 + t1000 * t721 + t1002 * t647, qJ(1) * t620 + t1000 * t561 + t1002 * t555, t1000 * t748 + t1002 * t712, t1000 * t703 + t1002 * t641, t1000 * t731 + t1002 * t678, t1000 * t747 + t1002 * t711, t1000 * t732 + t1002 * t679, t1000 * t789 + t1002 * t734, qJ(1) * t613 + t1000 * t505 + t1002 * t489, qJ(1) * t617 + t1000 * t513 + t1002 * t495, qJ(1) * t587 + t1000 * t488 + t1002 * t468, qJ(1) * t494 + t1000 * t443 + t1002 * t436, t1000 * t593 + t1002 * t549, t1000 * t538 + t1002 * t501, t1000 * t601 + t1002 * t550, t1000 * t592 + t1002 * t548, t1000 * t602 + t1002 * t551, t1000 * t654 + t1002 * t614, qJ(1) * t498 + t1000 * t425 + t1002 * t415, qJ(1) * t503 + t1000 * t429 + t1002 * t419, qJ(1) * t459 + t1000 * t398 + t1002 * t395, qJ(1) * t410 + t1000 * t387 + t1002 * t380, t1000 * t465 + t1002 * t442, t1000 * t433 + t1002 * t417, t1000 * t484 + t1002 * t453, t1000 * t464 + t1002 * t441, t1000 * t485 + t1002 * t454, t1000 * t530 + t1002 * t496, qJ(1) * t423 + t1000 * t389 + t1002 * t381, qJ(1) * t431 + t1000 * t390 + t1002 * t385, qJ(1) * t397 + t1000 * t375 + t1002 * t373, qJ(1) * t378 + t1000 * t372 + t1002 * t371; 0, 0, 0, 0, 0, 0, -g(2), g(1), 0, 0, 0, 0, 0, 0, 0, 0, t1094, t972, 0, 0, 0, 0, t943, 0, t944, t1003 * qJDD(2), pkin(1) * t946 - t1003 * t893 + (-pkin(7) * t1043 + t1039) * t1001, -pkin(1) * t945 - t1003 * t894 + (-pkin(7) * t967 - t1040) * t1001, t1163, pkin(1) * t795 + pkin(7) * t1163, t1001 * t1023 + t1003 * t924, t1001 * t1049 + t1003 * t915, t1001 * t1042 + t1003 * t927, t1001 * t1024 - t1003 * t923, t1001 * t1041 + t1003 * t925, t1044 * t1001, pkin(1) * t835 + t1001 * t1026 + t1003 * t819, pkin(1) * t836 + t1001 * t1025 + t1003 * t818, pkin(1) * t910 + t1001 * t1087 + t1003 * t741, pkin(1) * t672 + t1001 * t1018 + t1003 * t733, t1001 * t1052 + t1003 * t809, t1001 * t1058 + t1003 * t768, t1001 * t1055 + t1003 * t781, t1001 * t1053 + t1003 * t808, t1001 * t1054 + t1003 * t782, t1001 * t1051 + t1003 * t861, pkin(1) * t666 + t1001 * t1028 + t1003 * t603, pkin(1) * t674 + t1001 * t1027 + t1003 * t605, pkin(1) * t638 + t1001 * t1029 + t1003 * t569, pkin(1) * t537 + t1001 * t1030 + t1003 * t514, t1001 * t1062 + t1003 * t651, t1001 * t1069 + t1003 * t597, t1001 * t1065 + t1003 * t642, t1001 * t1063 + t1003 * t650, t1001 * t1064 + t1003 * t643, t1001 * t1061 + t1003 * t728, pkin(1) * t540 + t1001 * t1032 + t1003 * t475, pkin(1) * t547 + t1001 * t1031 + t1003 * t478, pkin(1) * t500 + t1001 * t1033 + t1003 * t438, pkin(1) * t427 + t1001 * t1034 + t1003 * t404, t1001 * t1076 + t1003 * t507, t1001 * t1079 + t1003 * t466, t1001 * t1073 + t1003 * t531, t1001 * t1077 + t1003 * t506, t1001 * t1072 + t1003 * t532, t1001 * t1070 + t1003 * t594, pkin(1) * t445 + t1001 * t1036 + t1003 * t399, pkin(1) * t451 + t1001 * t1035 + t1003 * t408, pkin(1) * t414 + t1001 * t1037 + t1003 * t386, pkin(1) * t384 + t1001 * t1038 + t1003 * t374;];
tauB_reg  = t1;