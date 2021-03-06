% Calculate minimal parameter regressor of Coriolis joint torque vector for
% S6RRRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
% MDP [38x1]
%   Minimal dynamic parameter vector (fixed base model)
%   see S6RRRRRR9_convert_par2_MPV_fixb.m
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRRR9_coriolisvecJ_fixb_mdp_slag_vp(qJ, qJD, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(38,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_coriolisvecJ_fixb_mdp_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_coriolisvecJ_fixb_mdp_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_coriolisvecJ_fixb_mdp_slag_vp: pkin has to be [13x1] (double)');
assert(isreal(MDP) && all(size(MDP) == [38 1]), ...
  'S6RRRRRR9_coriolisvecJ_fixb_mdp_slag_vp: MDP has to be [38x1] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:40:01
% EndTime: 2019-03-10 05:40:53
% DurationCPUTime: 30.29s
% Computational Cost: add. (31854->894), mult. (93270->1218), div. (0->0), fcn. (79137->14), ass. (0->356)
t965 = cos(qJ(3));
t1083 = qJD(3) * t965;
t954 = sin(pkin(7));
t1058 = t954 * t1083;
t955 = sin(pkin(6));
t1086 = qJD(1) * t955;
t966 = cos(qJ(2));
t1116 = t965 * t966;
t960 = sin(qJ(3));
t961 = sin(qJ(2));
t1120 = t960 * t961;
t956 = cos(pkin(7));
t992 = -t1120 * t956 + t1116;
t880 = t992 * t1086;
t1181 = -t880 + t1058;
t1124 = t956 * t960;
t1127 = t954 * t960;
t1123 = t956 * t965;
t941 = pkin(10) * t1127;
t1153 = pkin(2) * t1123 - t941;
t1143 = cos(pkin(6));
t1064 = pkin(1) * t1143;
t1028 = t966 * t1064;
t938 = qJD(1) * t1028;
t1020 = t955 * (-pkin(10) * t956 - pkin(9));
t996 = t961 * t1020;
t859 = qJD(1) * t996 + t938;
t1029 = t961 * t1064;
t972 = t1020 * t966 - t1029;
t860 = t972 * qJD(1);
t1144 = pkin(10) * t954;
t982 = (pkin(2) * t961 - t1144 * t966) * t955;
t893 = qJD(1) * t982;
t1180 = t1153 * qJD(3) - t860 * t1124 - t893 * t1127 - t965 * t859;
t797 = -t860 * t954 + t956 * t893;
t1118 = t961 * t965;
t1119 = t960 * t966;
t994 = t1118 * t956 + t1119;
t879 = t994 * t1086;
t1179 = pkin(3) * t879 - pkin(11) * t880 + t797 - (pkin(3) * t960 - pkin(11) * t965) * t954 * qJD(3);
t1062 = t961 * t1086;
t1027 = t954 * t1062;
t1178 = pkin(11) * t1027 - t1180;
t959 = sin(qJ(4));
t964 = cos(qJ(4));
t905 = t1127 * t959 - t964 * t956;
t1097 = -qJD(4) * t905 - t1027 * t959 + t1181 * t964;
t906 = t1127 * t964 + t956 * t959;
t1096 = qJD(4) * t906 + t1027 * t964 + t1181 * t959;
t1126 = t954 * t965;
t1090 = pkin(2) * t1124 + pkin(10) * t1126;
t1177 = t1090 * qJD(3) - t960 * t859;
t1026 = t960 * t1062;
t1061 = t966 * t1086;
t1024 = t956 * t1061;
t917 = t965 * t1024;
t1047 = t1143 * qJD(1);
t999 = t1047 + qJD(2);
t983 = t999 * t954;
t823 = -t965 * t983 + t1026 - t917;
t821 = qJD(4) + t823;
t1084 = qJD(3) * t960;
t1059 = t954 * t1084;
t1014 = -t879 + t1059;
t1093 = t860 * t1123 - (-pkin(3) * t1062 - t893 * t965) * t954 + t1177;
t1080 = qJD(4) * t964;
t1082 = qJD(4) * t959;
t897 = pkin(11) * t956 + t1090;
t898 = (-pkin(3) * t965 - pkin(11) * t960 - pkin(2)) * t954;
t1176 = -t898 * t1080 + t1082 * t897 + t1178 * t964 + t1179 * t959;
t1175 = -t897 * t1080 - t898 * t1082 + t1178 * t959 - t1179 * t964;
t1174 = -pkin(12) * t1014 + t1176;
t1125 = t955 * t966;
t978 = pkin(9) * t1125 + t1029;
t817 = t978 * qJD(1) + (t983 + t1024) * pkin(10);
t973 = pkin(2) * t1143 + t996;
t822 = qJD(2) * pkin(2) + qJD(1) * t973 + t938;
t889 = (-pkin(2) * t966 - t1144 * t961 - pkin(1)) * t955;
t875 = qJD(1) * t889;
t738 = t822 * t1124 + t875 * t1127 + t965 * t817;
t1173 = t738 - t821 * (pkin(4) * t959 - pkin(12) * t964);
t1172 = -t1096 * pkin(4) + pkin(12) * t1097 - t1093;
t993 = t1119 * t956 + t1118;
t980 = t993 * t955;
t825 = qJD(1) * t980 + t960 * t983;
t958 = sin(qJ(5));
t963 = cos(qJ(5));
t760 = -t823 * t958 * t964 - t963 * t825;
t1161 = -t958 * t1080 + t760;
t1117 = t963 * t964;
t761 = -t1117 * t823 + t825 * t958;
t1013 = t963 * t1080 - t761;
t1079 = qJD(5) * t958;
t1170 = -t959 * t1079 + t1013;
t774 = -t822 * t954 + t956 * t875;
t715 = pkin(3) * t823 - pkin(11) * t825 + t774;
t928 = t954 * t1061;
t1075 = t928 - qJD(3);
t974 = -t956 * t999 + t1075;
t721 = -pkin(11) * t974 + t738;
t648 = t715 * t964 - t959 * t721;
t643 = -pkin(4) * t821 - t648;
t783 = t964 * t825 - t959 * t974;
t744 = t783 * t958 - t963 * t821;
t633 = pkin(5) * t744 + t643;
t746 = t783 * t963 + t821 * t958;
t957 = sin(qJ(6));
t1132 = t746 * t957;
t962 = cos(qJ(6));
t670 = t962 * t744 + t1132;
t1169 = t633 * t670;
t874 = t964 * t974;
t781 = t825 * t959 + t874;
t779 = qJD(5) + t781;
t778 = qJD(6) + t779;
t1168 = t670 * t778;
t1078 = qJD(5) * t963;
t1166 = t1078 * t959 - t1161;
t865 = t1126 * t963 + t906 * t958;
t1104 = -qJD(5) * t865 + t1014 * t958 + t1097 * t963;
t1070 = t958 * t1126;
t1103 = -qJD(5) * t1070 - t1014 * t963 + t1078 * t906 + t1097 * t958;
t1004 = t744 * t957 - t962 * t746;
t1165 = t1004 * t778;
t1164 = t633 * t1004;
t1151 = qJD(5) + qJD(6);
t1162 = t1151 + t781;
t1074 = qJD(1) * qJD(2);
t1050 = t955 * t1074;
t1023 = t961 * t1050;
t1000 = t954 * t1023;
t1022 = t966 * t1050;
t1032 = -t956 * qJD(2) - qJD(3);
t979 = qJD(3) * t983;
t784 = qJD(3) * t917 + t1026 * t1032 + (t1022 + t979) * t965;
t1040 = -t964 * t1000 + t959 * t784;
t708 = qJD(4) * t783 + t1040;
t1160 = t708 * MDP(36) + (t1004 ^ 2 - t670 ^ 2) * MDP(33) - t670 * MDP(32) * t1004;
t951 = t955 ^ 2;
t1159 = -0.2e1 * t951 * t1074;
t1158 = MDP(6) * t966;
t1157 = t821 * t959;
t921 = t957 * t963 + t958 * t962;
t902 = t921 * t959;
t858 = t1028 + t973;
t793 = -t858 * t954 + t956 * t889;
t1049 = t1143 * t954;
t1150 = t956 * t1116 - t1120;
t849 = -t1049 * t965 - t1150 * t955;
t850 = t1049 * t960 + t980;
t733 = pkin(3) * t849 - pkin(11) * t850 + t793;
t842 = (t1125 * t956 + t1049) * pkin(10) + t978;
t1067 = t858 * t1124 + t889 * t1127 + t965 * t842;
t1048 = t1143 * t956;
t904 = t1125 * t954 - t1048;
t741 = -pkin(11) * t904 + t1067;
t1107 = t959 * t733 + t964 * t741;
t658 = pkin(12) * t849 + t1107;
t1148 = (t858 * t956 + t889 * t954) * t965 - t960 * t842;
t740 = pkin(3) * t904 - t1148;
t802 = t850 * t959 + t904 * t964;
t803 = t850 * t964 - t904 * t959;
t679 = pkin(4) * t802 - pkin(12) * t803 + t740;
t1114 = t963 * t658 + t958 * t679;
t1105 = -pkin(4) * t1014 - t1175;
t1156 = t1172 * t963;
t896 = t941 + (-pkin(2) * t965 - pkin(3)) * t956;
t804 = pkin(4) * t905 - pkin(12) * t906 + t896;
t1095 = t964 * t897 + t959 * t898;
t806 = -pkin(12) * t1126 + t1095;
t1100 = t958 * t804 + t963 * t806;
t1155 = MDP(5) * (t961 ^ 2 - t966 ^ 2);
t1071 = pkin(11) * t1082;
t1154 = t958 * t1071 - t1173 * t963;
t1152 = -t804 * t1078 + t1079 * t806 + t1172 * t958 + t1174 * t963;
t737 = (t822 * t956 + t875 * t954) * t965 - t960 * t817;
t766 = pkin(3) * t825 + pkin(11) * t823;
t1108 = t964 * t737 + t959 * t766;
t666 = pkin(12) * t825 + t1108;
t931 = -pkin(4) * t964 - pkin(12) * t959 - pkin(3);
t1149 = -t931 * t1078 + t1173 * t958 + t963 * t666;
t707 = -qJD(4) * t874 + t959 * t1000 - t1082 * t825 + t964 * t784;
t1147 = qJD(2) * t994 + qJD(3) * t993;
t968 = t1147 * t955;
t975 = t960 * t979;
t785 = qJD(1) * t968 + t975;
t642 = qJD(5) * t746 + t707 * t958 - t963 * t785;
t641 = t821 * t1078 - t1079 * t783 + t963 * t707 + t958 * t785;
t1044 = t641 * t957 + t962 * t642;
t602 = -qJD(6) * t1004 + t1044;
t967 = qJD(1) ^ 2;
t1146 = pkin(12) + pkin(13);
t1145 = pkin(5) * t959;
t1056 = t956 * t1083;
t1015 = qJD(2) * t1047;
t997 = pkin(1) * t1015;
t932 = t966 * t997;
t988 = qJD(2) * t996;
t839 = qJD(1) * t988 + t932;
t862 = t972 * qJD(2);
t840 = qJD(1) * t862;
t894 = qJD(2) * t982;
t884 = qJD(1) * t894;
t977 = -t822 * t1056 - t875 * t1058 + t1084 * t817 - t840 * t1124 - t884 * t1127 - t965 * t839;
t667 = pkin(11) * t1000 - t977;
t786 = -t840 * t954 + t956 * t884;
t692 = pkin(3) * t785 - pkin(11) * t784 + t786;
t610 = -t721 * t1080 - t715 * t1082 - t959 * t667 + t964 * t692;
t605 = -pkin(4) * t785 - t610;
t1142 = t605 * t958;
t1141 = t605 * t963;
t649 = t959 * t715 + t964 * t721;
t644 = pkin(12) * t821 + t649;
t720 = pkin(3) * t974 - t737;
t655 = t781 * pkin(4) - t783 * pkin(12) + t720;
t618 = -t644 * t958 + t963 * t655;
t612 = -pkin(13) * t746 + t618;
t606 = pkin(5) * t779 + t612;
t1140 = t606 * t962;
t619 = t644 * t963 + t655 * t958;
t613 = -pkin(13) * t744 + t619;
t1139 = t613 * t962;
t1138 = t641 * t958;
t1137 = t708 * t958;
t1136 = t708 * t963;
t1135 = t708 * t964;
t1134 = t744 * t779;
t1133 = t746 * t779;
t1131 = t781 * t821;
t1130 = t781 * t958;
t1129 = t783 * t821;
t1128 = t951 * t967;
t1122 = t958 * t959;
t1121 = t959 * t963;
t725 = pkin(4) * t783 + pkin(12) * t781;
t1115 = t963 * t648 + t958 * t725;
t866 = t906 * t963 - t1070;
t787 = t962 * t865 + t866 * t957;
t1111 = -qJD(6) * t787 - t1103 * t957 + t1104 * t962;
t788 = -t865 * t957 + t866 * t962;
t1110 = qJD(6) * t788 + t1103 * t962 + t1104 * t957;
t1109 = pkin(5) * t1103 + t1105;
t920 = t957 * t958 - t962 * t963;
t1102 = -t920 * t1080 - t1151 * t902 + t760 * t957 - t761 * t962;
t1077 = qJD(6) * t957;
t1101 = -t1077 * t1122 + (t1121 * t1151 - t1161) * t962 + t1170 * t957;
t1099 = t1162 * t920;
t1098 = t1162 * t921;
t731 = t959 * t737;
t665 = -pkin(4) * t825 - t766 * t964 + t731;
t1094 = pkin(5) * t1166 + pkin(11) * t1080 - t665;
t945 = pkin(11) * t1117;
t1089 = t958 * t931 + t945;
t1085 = qJD(2) * t955;
t1081 = qJD(4) * t963;
t1076 = qJD(6) * t962;
t1073 = pkin(1) * t1128;
t1068 = -t744 * t1076 + t962 * t641 - t957 * t642;
t1063 = qJD(5) * t1146;
t1060 = t961 * t1085;
t1057 = t956 * t1084;
t1053 = t779 * t1079;
t609 = t715 * t1080 - t1082 * t721 + t964 * t667 + t959 * t692;
t604 = pkin(12) * t785 + t609;
t678 = -t822 * t1057 - t875 * t1059 - t817 * t1083 + t840 * t1123 + t884 * t1126 - t960 * t839;
t668 = -pkin(3) * t1000 - t678;
t626 = pkin(4) * t708 - pkin(12) * t707 + t668;
t595 = -qJD(5) * t619 - t604 * t958 + t963 * t626;
t591 = pkin(5) * t708 - pkin(13) * t641 + t595;
t594 = t655 * t1078 - t1079 * t644 + t963 * t604 + t958 * t626;
t592 = -pkin(13) * t642 + t594;
t1046 = t962 * t591 - t957 * t592;
t611 = t613 * t1077;
t1045 = t957 * t591 - t611;
t1043 = -t658 * t958 + t963 * t679;
t1041 = t733 * t964 - t959 * t741;
t1039 = t963 * t804 - t806 * t958;
t798 = -t862 * t954 + t956 * t894;
t1038 = -t959 * t897 + t898 * t964;
t1037 = t779 * t963;
t1034 = t821 * t964;
t1033 = qJD(6) * t606 + t592;
t1031 = MDP(4) * t951 * t961 * t966;
t939 = qJD(2) * t1028;
t861 = t939 + t988;
t1030 = -t858 * t1057 - t889 * t1059 - t842 * t1083 - t960 * t861;
t1025 = t954 * t1060;
t1021 = -t649 + (t1079 + t1130) * pkin(5);
t1019 = qJD(3) * t1049;
t1017 = pkin(1) * t1159;
t805 = pkin(4) * t1126 - t1038;
t724 = t963 * t725;
t937 = t1146 * t963;
t1012 = pkin(5) * t783 + qJD(6) * t937 - t648 * t958 + t724 + (pkin(13) * t781 + t1063) * t963;
t869 = -pkin(13) * t1122 + t1089;
t1011 = -pkin(13) * t761 + qJD(6) * t869 - t1145 * t823 - t666 * t958 - (-pkin(13) * t1117 + t1145) * qJD(4) - (-t945 + (pkin(13) * t959 - t931) * t958) * qJD(5) - t1154;
t722 = -pkin(13) * t865 + t1100;
t1010 = -pkin(5) * t1096 + pkin(13) * t1104 + qJD(5) * t1100 + qJD(6) * t722 - t1174 * t958 + t1156;
t919 = t963 * t931;
t848 = -pkin(13) * t1121 + t919 + (-pkin(11) * t958 - pkin(5)) * t964;
t1009 = -qJD(6) * t848 - (-t1079 * t964 - t1081 * t959) * pkin(11) + t1149 + t1166 * pkin(13);
t936 = t1146 * t958;
t1008 = pkin(13) * t1130 + qJD(6) * t936 + t958 * t1063 + t1115;
t704 = pkin(5) * t905 - pkin(13) * t866 + t1039;
t1007 = pkin(13) * t1103 - qJD(6) * t704 + t1152;
t598 = t606 * t957 + t1139;
t758 = t803 * t963 + t849 * t958;
t617 = pkin(5) * t802 - pkin(13) * t758 + t1043;
t757 = t803 * t958 - t849 * t963;
t620 = -pkin(13) * t757 + t1114;
t1006 = t617 * t962 - t620 * t957;
t1005 = t617 * t957 + t620 * t962;
t696 = t962 * t757 + t758 * t957;
t697 = -t757 * t957 + t758 * t962;
t1001 = t954 ^ 2 * t1023;
t976 = t858 * t1056 + t889 * t1058 - t1084 * t842 + t862 * t1124 + t894 * t1127 + t965 * t861;
t685 = pkin(11) * t1025 + t976;
t795 = t1019 * t960 + t968;
t796 = t965 * t1019 + (t992 * qJD(2) + qJD(3) * t1150) * t955;
t699 = pkin(3) * t795 - pkin(11) * t796 + t798;
t995 = -t741 * t1080 - t733 * t1082 - t959 * t685 + t699 * t964;
t657 = -pkin(4) * t849 - t1041;
t991 = -t1078 * t779 - t1137;
t990 = -pkin(11) * t785 + t720 * t821;
t989 = -pkin(12) * t708 + t643 * t779;
t987 = t733 * t1080 - t1082 * t741 + t964 * t685 + t959 * t699;
t615 = pkin(12) * t795 + t987;
t686 = -t862 * t1123 + (-pkin(3) * t1060 - t894 * t965) * t954 - t1030;
t726 = qJD(4) * t803 - t1025 * t964 + t796 * t959;
t727 = -qJD(4) * t802 + t1025 * t959 + t796 * t964;
t632 = pkin(4) * t726 - pkin(12) * t727 + t686;
t985 = t679 * t1078 - t1079 * t658 + t963 * t615 + t958 * t632;
t601 = -t1077 * t746 + t1068;
t616 = -pkin(4) * t795 - t995;
t971 = -qJD(5) * t1114 - t615 * t958 + t963 * t632;
t589 = -qJD(6) * t598 + t1046;
t970 = qJD(3) * t974;
t969 = t999 * t978;
t949 = -pkin(5) * t963 - pkin(4);
t926 = (pkin(5) * t958 + pkin(11)) * t959;
t903 = t920 * t959;
t762 = pkin(5) * t865 + t805;
t703 = t708 * t905;
t684 = t708 * t802;
t652 = -qJD(5) * t757 + t727 * t963 + t795 * t958;
t651 = qJD(5) * t758 + t727 * t958 - t795 * t963;
t637 = pkin(5) * t757 + t657;
t608 = qJD(6) * t697 + t962 * t651 + t652 * t957;
t607 = -qJD(6) * t696 - t651 * t957 + t652 * t962;
t600 = pkin(5) * t651 + t616;
t599 = pkin(5) * t642 + t605;
t597 = -t613 * t957 + t1140;
t596 = -pkin(13) * t651 + t985;
t593 = pkin(5) * t726 - pkin(13) * t652 + t971;
t588 = t1033 * t962 + t1045;
t1 = [(-MDP(7) * t1060 + t1085 * t1158) * (0.2e1 * t1047 + qJD(2)) + (-((t862 * t956 + t894 * t954) * t965 + t1030) * t974 - t678 * t904 + t798 * t823 + t793 * t785 + t786 * t849 + t774 * t795 + (qJD(1) * t1148 + t737) * t1025) * MDP(16) + (-qJD(2) * t969 - t1015 * t978 + t1017 * t961) * MDP(9) + (t1041 * t785 + t610 * t849 + t648 * t795 + t668 * t802 + t686 * t781 + t740 * t708 + t720 * t726 + t821 * t995) * MDP(23) + (-t1107 * t785 - t609 * t849 - t649 * t795 + t668 * t803 + t686 * t783 + t740 * t707 + t720 * t727 - t821 * t987) * MDP(24) + (-(-pkin(9) * t1060 + t939) * t999 - (-pkin(9) * t1023 + t932) * t1143 + t966 * t1017) * MDP(10) + (-t1114 * t708 - t594 * t802 + t605 * t758 + t616 * t746 - t619 * t726 + t657 * t641 + t643 * t652 - t779 * t985) * MDP(31) + (-t1004 * t726 + t601 * t802 + t607 * t778 + t697 * t708) * MDP(34) + (-(qJD(6) * t1006 + t593 * t957 + t596 * t962) * t778 - t1005 * t708 - t588 * t802 - t598 * t726 - t600 * t1004 + t637 * t601 + t599 * t697 + t633 * t607) * MDP(38) + (-t1004 * t607 + t601 * t697) * MDP(32) + (t1004 * t608 - t601 * t696 - t602 * t697 - t607 * t670) * MDP(33) + (t976 * t974 - t977 * t904 + t798 * t825 + t793 * t784 + t786 * t850 + t774 * t796 + (-qJD(1) * t1067 - t738) * t1025) * MDP(17) + 0.2e1 * t1031 * t1074 + (t1043 * t708 + t595 * t802 + t605 * t757 + t616 * t744 + t618 * t726 + t657 * t642 + t643 * t651 + t779 * t971) * MDP(30) + (-t796 * t974 - t784 * t904 + (qJD(1) * t850 + t825) * t1025) * MDP(13) + (t795 * t974 + t785 * t904 + (-qJD(1) * t849 - t823) * t1025) * MDP(14) + (t707 * t803 + t727 * t783) * MDP(18) + (-t707 * t802 - t708 * t803 - t726 * t783 - t727 * t781) * MDP(19) + (t726 * t779 + t684) * MDP(29) + (-t642 * t802 - t651 * t779 - t708 * t757 - t726 * t744) * MDP(28) + (t641 * t802 + t652 * t779 + t708 * t758 + t726 * t746) * MDP(27) + (-t602 * t802 - t608 * t778 - t670 * t726 - t696 * t708) * MDP(35) + ((-qJD(6) * t1005 + t593 * t962 - t596 * t957) * t778 + t1006 * t708 + t589 * t802 + t597 * t726 + t600 * t670 + t637 * t602 + t599 * t696 + t633 * t608) * MDP(37) + (t726 * t778 + t684) * MDP(36) + (t641 * t758 + t652 * t746) * MDP(25) + (-t641 * t757 - t642 * t758 - t651 * t746 - t652 * t744) * MDP(26) + t1155 * t1159 + (-t928 + (t1048 - t904) * qJD(1) - t1032) * MDP(15) * t1025 + (-t708 * t849 - t726 * t821 - t781 * t795 - t785 * t802) * MDP(21) + (t707 * t849 + t727 * t821 + t783 * t795 + t785 * t803) * MDP(20) + (t785 * t849 + t795 * t821) * MDP(22) + (t784 * t850 + t796 * t825) * MDP(11) + (-t784 * t849 - t785 * t850 - t795 * t825 - t796 * t823) * MDP(12); (-t1096 * t619 - t1100 * t708 + t1104 * t643 + t1105 * t746 + t1152 * t779 - t594 * t905 + t605 * t866 + t805 * t641) * MDP(31) + (t1014 * t648 + t1038 * t785 + t1093 * t781 + t1096 * t720 - t610 * t1126 + t1175 * t821 + t668 * t905 + t896 * t708) * MDP(23) + (-pkin(9) * t1022 + qJD(1) * t969 + (t1073 - t997) * t961) * MDP(9) + (t960 * t1001 + t784 * t956 + t880 * t974 + (-t1062 * t825 - t965 * t970) * t954) * MDP(13) + (t965 * t1001 - t785 * t956 - t879 * t974 + (t1062 * t823 + t960 * t970) * t954) * MDP(14) + (-t1031 + (MDP(7) * t961 - t1158) * t955 * t1143) * t967 + (-t1014 * t649 + t1093 * t783 - t1095 * t785 + t1097 * t720 + t609 * t1126 + t1176 * t821 + t668 * t906 + t896 * t707) * MDP(24) + (-t1096 * t744 - t1103 * t779 - t642 * t905 - t708 * t865) * MDP(28) + (t1096 * t746 + t1104 * t779 + t641 * t905 + t708 * t866) * MDP(27) + (t1104 * t746 + t641 * t866) * MDP(25) + (-t1103 * t746 - t1104 * t744 - t641 * t865 - t642 * t866) * MDP(26) + ((t704 * t962 - t722 * t957) * t708 + t589 * t905 + t762 * t602 + t599 * t787 + (t1007 * t957 - t1010 * t962) * t778 + t1109 * t670 + t1110 * t633 + t1096 * t597) * MDP(37) + (-t1096 * t670 - t1110 * t778 - t602 * t905 - t708 * t787) * MDP(35) + (t1014 * t821 - t1126 * t785) * MDP(22) + (t1014 * t783 + t1097 * t821 - t1126 * t707 + t785 * t906) * MDP(20) + (-t1014 * t781 - t1096 * t821 + t1126 * t708 - t785 * t905) * MDP(21) + (t1004 * t1110 - t1111 * t670 - t601 * t787 - t602 * t788) * MDP(33) + (-(t704 * t957 + t722 * t962) * t708 - t588 * t905 + t762 * t601 + t599 * t788 + (t1007 * t962 + t1010 * t957) * t778 - t1109 * t1004 + t1111 * t633 - t1096 * t598) * MDP(38) + (-t1004 * t1096 + t1111 * t778 + t601 * t905 + t708 * t788) * MDP(34) + (-t1004 * t1111 + t601 * t788) * MDP(32) + (-t932 + (-pkin(9) * t1062 + t938) * t1047 + t966 * t1073 + t938 * qJD(2)) * MDP(10) + (t823 * t880 + t825 * t879 + (t784 * t965 - t785 * t960 + (-t823 * t965 - t825 * t960) * qJD(3)) * t954) * MDP(12) + (-t1096 * t783 - t1097 * t781 - t707 * t905 - t708 * t906) * MDP(19) + (t1096 * t778 + t703) * MDP(36) + (t1096 * t779 + t703) * MDP(29) + (t1097 * t783 + t707 * t906) * MDP(18) - (t1047 * t956 - t1075) * MDP(15) * t1027 + t1128 * t1155 + (t977 * t956 - t797 * t825 - t774 * t880 + (t774 * t1083 - pkin(2) * t784 + t786 * t960 + (-qJD(2) * t1090 + t738) * t1062) * t954 + t1180 * t974) * MDP(17) + (t784 * t1127 + t1181 * t825) * MDP(11) + (t678 * t956 - t954 * pkin(2) * t785 - t786 * t1126 - t797 * t823 + (qJD(2) * t1153 - t737) * t1027 + t1014 * t774 + ((t860 * t956 + t893 * t954) * t965 + t1177) * t974) * MDP(16) + (t1039 * t708 + t595 * t905 + t805 * t642 + t605 * t865 + (-t806 * t1078 + (-qJD(5) * t804 + t1174) * t958 - t1156) * t779 + t1105 * t744 + t1103 * t643 + t1096 * t618) * MDP(30); (t641 * t1121 + t1170 * t746) * MDP(25) + (-t643 * t760 - t665 * t744 + t919 * t708 + ((-qJD(5) * t931 + t666) * t958 + t1154) * t779 + (t643 * t958 * qJD(4) - t595 + (qJD(4) * t744 + t991) * pkin(11)) * t964 + (pkin(11) * t642 + t1078 * t643 + t618 * t821 + t1142) * t959) * MDP(30) + (-t1089 * t708 - t665 * t746 - t643 * t761 + t1149 * t779 + (t643 * t1081 + t594 + (qJD(4) * t746 + t1053) * pkin(11)) * t964 + (-t643 * t1079 + t1141 - t821 * t619 + (t1081 * t779 + t641) * pkin(11)) * t959) * MDP(31) + (-t1086 * t1147 - t975) * MDP(14) + (t1034 * t783 + t707 * t959) * MDP(18) + (-pkin(3) * t707 + t668 * t959 - t738 * t783 + (t1071 + t1108) * t821 + t990 * t964) * MDP(24) + (-pkin(3) * t708 - t668 * t964 - t738 * t781 + (t731 + (-pkin(11) * qJD(4) - t766) * t964) * t821 + t990 * t959) * MDP(23) + (MDP(11) * t823 + MDP(12) * t825 - MDP(14) * t974 - MDP(16) * t774 - MDP(20) * t783 + MDP(21) * t781 - MDP(22) * t821 - MDP(23) * t648 + MDP(24) * t649) * t825 + (-t823 * t974 + t784) * MDP(13) + MDP(15) * t1000 + (-t1004 * t1102 - t601 * t903) * MDP(32) + (t1004 * t1101 - t1102 * t670 - t601 * t902 + t602 * t903) * MDP(33) + (t1157 * t779 - t1135) * MDP(29) + (-t1004 * t1157 + t1102 * t778 - t601 * t964 - t708 * t903) * MDP(34) + ((t848 * t962 - t869 * t957) * t708 - t589 * t964 + t926 * t602 + t599 * t902 + (t1009 * t957 - t1011 * t962) * t778 + t1094 * t670 + t1101 * t633 + t597 * t1157) * MDP(37) + (t1157 * t778 - t1135) * MDP(36) + (-t1157 * t821 + t785 * t964) * MDP(21) + (-t1101 * t778 - t1157 * t670 + t602 * t964 - t708 * t902) * MDP(35) + (-(t848 * t957 + t869 * t962) * t708 + t588 * t964 + t926 * t601 - t599 * t903 + (t1009 * t962 + t1011 * t957) * t778 - t1094 * t1004 + t1102 * t633 - t598 * t1157) * MDP(38) + ((t707 - t1131) * t964 + (-t708 - t1129) * t959) * MDP(19) + (t1034 * t821 + t785 * t959) * MDP(20) + (-t738 * t974 + t678) * MDP(16) + (-t737 * t974 + t774 * t823 + t977) * MDP(17) + (t744 * t761 + t746 * t760 + (-t744 * t963 - t746 * t958) * t1080 + (-t1138 - t642 * t963 + (t744 * t958 - t746 * t963) * qJD(5)) * t959) * MDP(26) + (-t641 * t964 + t1013 * t779 + (t746 * t821 - t1053 + t1136) * t959) * MDP(27) - t823 ^ 2 * MDP(12) + (t642 * t964 + t1161 * t779 + (-t744 * t821 + t991) * t959) * MDP(28); -t781 ^ 2 * MDP(19) + (t707 + t1131) * MDP(20) + (-t1040 + t1129) * MDP(21) + t785 * MDP(22) + (t649 * t821 + t610) * MDP(23) + (t648 * t821 + t720 * t781 - t609) * MDP(24) + (t1037 * t746 + t1138) * MDP(25) + ((t641 - t1134) * t963 + (-t642 - t1133) * t958) * MDP(26) + (t1037 * t779 + t1137) * MDP(27) + (-t779 ^ 2 * t958 + t1136) * MDP(28) + (-pkin(4) * t642 - t1141 - t649 * t744 + (-pkin(12) * t1078 - t724) * t779 + (t648 * t779 + t989) * t958) * MDP(30) + (-pkin(4) * t641 + t1142 - t649 * t746 + (pkin(12) * t1079 + t1115) * t779 + t989 * t963) * MDP(31) + (t1004 * t1099 + t601 * t921) * MDP(32) + (t1004 * t1098 + t1099 * t670 - t601 * t920 - t602 * t921) * MDP(33) + (-t1099 * t778 + t708 * t921) * MDP(34) + (-t1098 * t778 - t708 * t920) * MDP(35) + ((-t936 * t962 - t937 * t957) * t708 + t949 * t602 + t599 * t920 + (t1008 * t957 - t1012 * t962) * t778 + t1021 * t670 + t1098 * t633) * MDP(37) + (-(-t936 * t957 + t937 * t962) * t708 + t949 * t601 + t599 * t921 + (t1008 * t962 + t1012 * t957) * t778 - t1021 * t1004 - t1099 * t633) * MDP(38) + (MDP(18) * t781 + MDP(19) * t783 - MDP(21) * qJD(4) - MDP(23) * t720 - MDP(27) * t746 + MDP(28) * t744 - MDP(29) * t779 - MDP(30) * t618 + MDP(31) * t619 + MDP(34) * t1004 + MDP(35) * t670 - MDP(36) * t778 - MDP(37) * t597 + MDP(38) * t598) * t783; t746 * t744 * MDP(25) + (-t744 ^ 2 + t746 ^ 2) * MDP(26) + (t641 + t1134) * MDP(27) + (t1133 - t642) * MDP(28) + t708 * MDP(29) + (t619 * t779 - t643 * t746 + t595) * MDP(30) + (t618 * t779 + t643 * t744 - t594) * MDP(31) + (t601 + t1168) * MDP(34) + (-t602 - t1165) * MDP(35) + (-(-t612 * t957 - t1139) * t778 + t1164 + (-t1077 * t778 - t670 * t746 + t708 * t962) * pkin(5) + t589) * MDP(37) + (t1169 + t611 + (-t613 * t778 - t591) * t957 + (t612 * t778 - t1033) * t962 + (t1004 * t746 - t1076 * t778 - t708 * t957) * pkin(5)) * MDP(38) + t1160; (t1068 + t1168) * MDP(34) + (-t1044 - t1165) * MDP(35) + (t598 * t778 + t1046 + t1164) * MDP(37) + (-t962 * t592 + t597 * t778 - t1045 + t1169) * MDP(38) + (-MDP(34) * t1132 + MDP(35) * t1004 - MDP(37) * t598 - MDP(38) * t1140) * qJD(6) + t1160;];
tauc  = t1;
