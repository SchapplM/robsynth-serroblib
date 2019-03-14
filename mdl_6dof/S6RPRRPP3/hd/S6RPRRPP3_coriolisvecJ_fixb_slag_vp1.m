% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [7x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:35:09
% EndTime: 2019-03-09 04:36:49
% DurationCPUTime: 92.24s
% Computational Cost: add. (41362->1185), mult. (55974->1488), div. (0->0), fcn. (53685->8), ass. (0->571)
t1066 = Icges(6,4) - Icges(5,5) - Icges(7,5);
t1065 = Icges(7,4) + Icges(6,5) - Icges(5,6);
t558 = cos(qJ(4));
t1105 = (-Icges(5,4) - Icges(6,6)) * t558;
t555 = sin(qJ(4));
t1104 = (Icges(5,4) - Icges(7,6)) * t555;
t1103 = -Icges(5,1) - Icges(6,2);
t1100 = Icges(5,2) + Icges(7,2);
t1099 = -t1104 + (Icges(5,1) + Icges(7,3)) * t558;
t1098 = -t1105 + (-Icges(5,2) - Icges(6,3)) * t555;
t556 = sin(qJ(3));
t559 = cos(qJ(3));
t1097 = t1066 * t559 + t1099 * t556;
t1096 = t1065 * t559 + t1098 * t556;
t1095 = Icges(6,1) + Icges(7,1) + Icges(5,3);
t1094 = t1065 * t555 - t1066 * t558;
t865 = t556 * t558;
t537 = Icges(7,6) * t865;
t867 = t555 * t556;
t1093 = Icges(7,2) * t867 - t1096 + t537;
t538 = Icges(6,6) * t867;
t1085 = -Icges(6,2) * t865 - t1097 + t538;
t1092 = (t1100 * t558 + t1104) * t556;
t1091 = (t1103 * t555 + t1105) * t556;
t554 = qJ(1) + pkin(9);
t550 = sin(t554);
t551 = cos(t554);
t866 = t555 * t559;
t411 = t550 * t866 + t551 * t558;
t392 = Icges(7,6) * t411;
t864 = t558 * t559;
t412 = t550 * t864 - t551 * t555;
t871 = t550 * t556;
t201 = Icges(7,5) * t871 + Icges(7,3) * t412 + t392;
t395 = Icges(6,6) * t411;
t211 = -Icges(6,4) * t871 + Icges(6,2) * t412 - t395;
t398 = Icges(5,4) * t411;
t226 = -Icges(5,1) * t412 - Icges(5,5) * t871 + t398;
t1078 = t201 + t211 - t226;
t414 = t550 * t555 + t551 * t864;
t869 = t551 * t556;
t413 = -t550 * t558 + t551 * t866;
t889 = Icges(7,6) * t413;
t203 = Icges(7,5) * t869 + Icges(7,3) * t414 + t889;
t397 = Icges(6,6) * t413;
t212 = Icges(6,4) * t869 - Icges(6,2) * t414 + t397;
t400 = Icges(5,4) * t413;
t227 = Icges(5,1) * t414 + Icges(5,5) * t869 - t400;
t1077 = t203 - t212 + t227;
t396 = Icges(6,6) * t412;
t204 = Icges(6,5) * t871 + Icges(6,3) * t411 - t396;
t393 = Icges(7,6) * t412;
t208 = -Icges(7,4) * t871 - Icges(7,2) * t411 - t393;
t399 = Icges(5,4) * t412;
t222 = -Icges(5,2) * t411 + Icges(5,6) * t871 + t399;
t1076 = t204 - t208 - t222;
t892 = Icges(6,6) * t414;
t206 = Icges(6,5) * t869 + Icges(6,3) * t413 - t892;
t394 = Icges(7,6) * t414;
t209 = Icges(7,4) * t869 + Icges(7,2) * t413 + t394;
t904 = Icges(5,4) * t414;
t224 = -Icges(5,2) * t413 + Icges(5,6) * t869 + t904;
t1075 = t206 + t209 - t224;
t213 = Icges(7,1) * t871 + Icges(7,4) * t411 + Icges(7,5) * t412;
t216 = Icges(6,1) * t871 - Icges(6,4) * t412 + Icges(6,5) * t411;
t219 = Icges(5,5) * t412 - Icges(5,6) * t411 + Icges(5,3) * t871;
t1074 = t213 + t216 + t219;
t215 = Icges(7,1) * t869 + Icges(7,4) * t413 + Icges(7,5) * t414;
t218 = Icges(6,1) * t869 - Icges(6,4) * t414 + Icges(6,5) * t413;
t221 = Icges(5,5) * t414 - Icges(5,6) * t413 + Icges(5,3) * t869;
t1049 = t215 + t218 + t221;
t1048 = t1094 * t556 - t1095 * t559;
t1090 = t1094 * t559 + t1095 * t556;
t887 = Icges(7,6) * t558;
t679 = Icges(7,2) * t555 + t887;
t1006 = (t679 - t1098) * t559 + t1065 * t556;
t891 = Icges(6,6) * t555;
t680 = Icges(6,2) * t558 - t891;
t1005 = (t680 + t1099) * t559 - t1066 * t556;
t1087 = (t1065 * t558 + t1066 * t555) * t556;
t1054 = t1074 * t869 + t1076 * t413 + t1078 * t414;
t1053 = t1049 * t869 + t1075 * t413 + t1077 * t414;
t811 = qJD(4) * t558;
t767 = t551 * t811;
t813 = qJD(4) * t555;
t768 = t550 * t813;
t815 = qJD(3) * t556;
t773 = t555 * t815;
t196 = qJD(1) * t411 + t551 * t773 - t559 * t767 - t768;
t819 = qJD(1) * t551;
t620 = t550 * t811 + t555 * t819;
t772 = t558 * t815;
t820 = qJD(1) * t550;
t782 = t559 * t820;
t810 = qJD(4) * t559;
t197 = t558 * t782 + (t555 * t810 + t772) * t551 - t620;
t814 = qJD(3) * t559;
t774 = t551 * t814;
t818 = qJD(1) * t556;
t783 = t550 * t818;
t621 = t774 - t783;
t105 = Icges(6,4) * t621 + Icges(6,2) * t197 - Icges(6,6) * t196;
t115 = -Icges(5,1) * t197 + Icges(5,4) * t196 + Icges(5,5) * t621;
t99 = Icges(7,5) * t621 - Icges(7,6) * t196 - Icges(7,3) * t197;
t1084 = -t105 + t115 + t99;
t198 = -t550 * t773 - t551 * t813 - t558 * t820 + t559 * t620;
t199 = qJD(1) * t414 - t550 * t772 - t559 * t768 - t767;
t1062 = t550 * t814;
t623 = t551 * t818 + t1062;
t100 = Icges(7,5) * t623 + Icges(7,6) * t198 + Icges(7,3) * t199;
t106 = Icges(6,4) * t623 - Icges(6,2) * t199 + Icges(6,6) * t198;
t116 = Icges(5,1) * t199 - Icges(5,4) * t198 + Icges(5,5) * t623;
t1083 = t100 - t106 + t116;
t101 = Icges(6,5) * t621 + Icges(6,6) * t197 - Icges(6,3) * t196;
t103 = Icges(7,4) * t621 - Icges(7,2) * t196 - Icges(7,6) * t197;
t113 = -Icges(5,4) * t197 + Icges(5,2) * t196 + Icges(5,6) * t621;
t1082 = t101 + t103 - t113;
t102 = Icges(6,5) * t623 - Icges(6,6) * t199 + Icges(6,3) * t198;
t104 = Icges(7,4) * t623 + Icges(7,2) * t198 + Icges(7,6) * t199;
t114 = Icges(5,4) * t199 - Icges(5,2) * t198 + Icges(5,6) * t623;
t1081 = t102 + t104 - t114;
t107 = Icges(7,1) * t621 - Icges(7,4) * t196 - Icges(7,5) * t197;
t109 = Icges(6,1) * t621 + Icges(6,4) * t197 - Icges(6,5) * t196;
t111 = -Icges(5,5) * t197 + Icges(5,6) * t196 + Icges(5,3) * t621;
t1080 = t107 + t109 + t111;
t108 = Icges(7,1) * t623 + Icges(7,4) * t198 + Icges(7,5) * t199;
t110 = Icges(6,1) * t623 - Icges(6,4) * t199 + Icges(6,5) * t198;
t112 = Icges(5,5) * t199 - Icges(5,6) * t198 + Icges(5,3) * t623;
t1079 = t108 + t110 + t112;
t979 = t1048 * t869 - t1085 * t414 + t1093 * t413;
t1073 = qJD(3) * t1090 + qJD(4) * t1087;
t812 = qJD(4) * t556;
t1072 = (Icges(6,3) * t558 + t891) * t812 + t1092 * qJD(4) + t1006 * qJD(3);
t1071 = (-Icges(7,3) * t555 + t887) * t812 + t1091 * qJD(4) + t1005 * qJD(3);
t1068 = -t1085 * t558 + t1093 * t555;
t1036 = rSges(7,1) + pkin(5);
t547 = t551 * pkin(7);
t476 = pkin(2) * t550 - t547;
t458 = qJD(1) * t476;
t495 = pkin(8) * t774;
t532 = pkin(7) * t819;
t557 = sin(qJ(1));
t775 = t551 * t815;
t522 = pkin(3) * t556 - pkin(8) * t559;
t816 = qJD(3) * t551;
t779 = t522 * t816;
t936 = pkin(1) * qJD(1);
t947 = pkin(3) * t559;
t523 = pkin(8) * t556 + t947;
t448 = t523 * t550;
t983 = qJD(1) * t448;
t1067 = -pkin(3) * t775 + t557 * t936 + t458 + t495 + t532 + t779 + t983;
t1056 = t1074 * t871 + t1076 * t411 + t1078 * t412;
t1055 = t1049 * t871 + t1075 * t411 + t1077 * t412;
t980 = t1048 * t871 - t1085 * t412 + t1093 * t411;
t1038 = t222 * t555 + t226 * t558;
t81 = -t1038 * t556 - t219 * t559;
t1064 = t550 * pkin(7);
t948 = pkin(1) * t557;
t748 = -t476 - t948;
t585 = (-t448 + t748) * qJD(1) - t779;
t715 = rSges(5,1) * t412 - rSges(5,2) * t411;
t236 = rSges(5,3) * t871 + t715;
t714 = rSges(5,1) * t558 - rSges(5,2) * t555;
t441 = -rSges(5,3) * t559 + t556 * t714;
t454 = -t550 * t812 + t816;
t539 = qJD(1) - t810;
t999 = -t236 * t539 - t441 * t454;
t90 = t585 + t999;
t1063 = t550 * t90;
t388 = qJD(5) * t413;
t117 = -t197 * pkin(4) - qJ(5) * t196 + t388;
t391 = t411 * qJ(5);
t946 = pkin(4) * t412;
t290 = t391 + t946;
t1059 = t290 * t539 + t1067 + t117 - t388;
t1017 = t1074 * t621 - t1076 * t196 - t1078 * t197 + t1079 * t869 + t1081 * t413 + t1083 * t414;
t1016 = t1049 * t621 - t1075 * t196 - t1077 * t197 + t1080 * t869 + t1082 * t413 + t1084 * t414;
t1015 = t1074 * t623 + t1076 * t198 + t1078 * t199 + t1079 * t871 + t1081 * t411 + t1083 * t412;
t1014 = t1049 * t623 + t1075 * t198 + t1077 * t199 + t1080 * t871 + t1082 * t411 + t1084 * t412;
t817 = qJD(3) * t550;
t453 = t551 * t812 + t817;
t981 = t1053 * t453 - t1054 * t454 + t539 * t979;
t1058 = t1048 * t623 + t1071 * t412 + t1072 * t411 + t1073 * t871 - t1085 * t199 + t1093 * t198;
t1057 = t1048 * t621 + t1071 * t414 + t1072 * t413 + t1073 * t869 + t1085 * t197 - t1093 * t196;
t675 = t201 * t558 - t208 * t555;
t83 = -t213 * t559 + t556 * t675;
t673 = t204 * t555 + t211 * t558;
t85 = -t216 * t559 + t556 * t673;
t1052 = t81 + t83 + t85;
t670 = -t224 * t555 + t227 * t558;
t82 = -t221 * t559 + t556 * t670;
t674 = t203 * t558 + t209 * t555;
t84 = -t215 * t559 + t556 * t674;
t672 = t206 * t555 - t212 * t558;
t86 = -t218 * t559 + t556 * t672;
t1051 = t82 + t84 + t86;
t978 = -t1048 * t559 + t1068 * t556;
t1047 = -t556 * t679 + t1096;
t1046 = t556 * t680 + t1097;
t1045 = (-t1068 + t1090) * t539 + (t1048 * t550 - t1038 + t673 + t675) * t454 + (-t1048 * t551 - t670 - t672 - t674) * t453;
t296 = t414 * pkin(4) + qJ(5) * t413;
t868 = t551 * t559;
t526 = pkin(3) * t868;
t450 = pkin(8) * t869 + t526;
t824 = t551 * pkin(2) + t1064;
t455 = qJD(1) * t824;
t560 = cos(qJ(1));
t549 = t560 * t936;
t830 = t549 + t455;
t735 = qJD(1) * t450 - t522 * t817 + t830;
t809 = qJD(5) * t411;
t709 = pkin(4) * t558 + qJ(5) * t555;
t996 = t709 * t556;
t1022 = -t539 * t296 + t453 * t996 - t735 - t809;
t710 = rSges(7,2) * t555 + rSges(7,3) * t558;
t1023 = -qJ(6) * t865 + t1036 * t559 - t556 * t710;
t384 = qJD(6) * t412;
t1035 = rSges(7,3) + qJ(6);
t855 = t413 * rSges(7,2) + t1035 * t414 + t1036 * t869;
t55 = t1023 * t453 + t855 * t539 - t1022 + t384;
t1044 = t550 * t551;
t982 = t1055 * t453 - t1056 * t454 + t539 * t980;
t1043 = t550 * t982 + t551 * t981;
t1042 = (qJD(3) * t1068 - t1073) * t559 + (t1071 * t558 + t1072 * t555 + (t1085 * t555 + t1093 * t558) * qJD(4) + t1048 * qJD(3)) * t556;
t1041 = Icges(7,3) - t1103;
t1040 = Icges(6,3) + t1100;
t1039 = -t1087 * t539 + (t1065 * t412 + t1066 * t411) * t454 + (-t1065 * t414 - t1066 * t413) * t453;
t561 = qJD(1) ^ 2;
t1037 = 0.2e1 * qJD(3);
t370 = t550 * t996;
t457 = t522 * t820;
t474 = t709 * t559;
t1021 = -t370 * t810 + t454 * t474 + t820 * t996 + t457;
t803 = qJD(3) * qJD(4);
t759 = t559 * t803;
t351 = qJD(1) * t453 + t550 * t759;
t352 = qJD(1) * t454 + t551 * t759;
t760 = t556 * t803;
t1019 = t1016 * t453 - t1017 * t454 + t1053 * t352 + t1054 * t351 + t1057 * t539 + t760 * t979;
t1018 = t1014 * t453 - t1015 * t454 + t1055 * t352 + t1056 * t351 + t1058 * t539 + t760 * t980;
t25 = (-qJD(3) * t1038 - t112) * t559 + (qJD(3) * t219 - t114 * t555 + t116 * t558 + (-t222 * t558 + t226 * t555) * qJD(4)) * t556;
t27 = (qJD(3) * t675 - t108) * t559 + (qJD(3) * t213 + t100 * t558 + t104 * t555 + (-t201 * t555 - t208 * t558) * qJD(4)) * t556;
t29 = (qJD(3) * t673 - t110) * t559 + (qJD(3) * t216 + t102 * t555 - t106 * t558 + (t204 * t558 - t211 * t555) * qJD(4)) * t556;
t1013 = t25 + t27 + t29;
t26 = (qJD(3) * t670 - t111) * t559 + (qJD(3) * t221 - t113 * t555 + t115 * t558 + (-t224 * t558 - t227 * t555) * qJD(4)) * t556;
t28 = (qJD(3) * t674 - t107) * t559 + (qJD(3) * t215 + t103 * t555 + t558 * t99 + (-t203 * t555 + t209 * t558) * qJD(4)) * t556;
t30 = (qJD(3) * t672 - t109) * t559 + (qJD(3) * t218 + t101 * t555 - t105 * t558 + (t206 * t558 + t212 * t555) * qJD(4)) * t556;
t1012 = t26 + t28 + t30;
t1011 = t1051 * t453 - t1052 * t454 + t539 * t978;
t1010 = t1047 * t550;
t1009 = t1047 * t551;
t1008 = t1046 * t550;
t1007 = t1046 * t551;
t1004 = t1045 * t556;
t1003 = t1048 * t539 + t1049 * t453 - t1074 * t454;
t973 = t1051 * t551 + t1052 * t550;
t1002 = t1051 * t550 - t1052 * t551;
t972 = t1053 * t551 + t1054 * t550;
t971 = t1055 * t551 + t1056 * t550;
t1001 = t1053 * t550 - t1054 * t551;
t1000 = t1055 * t550 - t1056 * t551;
t784 = -pkin(2) - t947;
t802 = -pkin(8) - t1036;
t624 = t802 * t556 + t784;
t997 = t550 * t624;
t651 = t411 * t539 + t454 * t867;
t992 = -t196 + t651;
t650 = -t413 * t539 + t453 * t867;
t991 = t198 + t650;
t939 = rSges(7,2) * t411;
t858 = t1035 * t412 + t1036 * t871 + t939;
t371 = t551 * t996;
t990 = t296 * t812 - t539 * t371;
t988 = -rSges(7,2) * t198 - t384;
t742 = t551 * rSges(3,1) - rSges(3,2) * t550;
t552 = Icges(4,4) * t559;
t685 = -Icges(4,2) * t556 + t552;
t501 = Icges(4,1) * t556 + t552;
t385 = qJD(6) * t414;
t987 = -t196 * rSges(7,2) - t1035 * t197 + t1036 * t774 + t385;
t235 = rSges(6,1) * t869 - t414 * rSges(6,2) + t413 * rSges(6,3);
t711 = rSges(6,2) * t558 - rSges(6,3) * t555;
t599 = rSges(6,1) * t559 + t556 * t711;
t66 = t235 * t539 + t453 * t599 - t1022;
t976 = (Icges(6,3) * t865 + t1085 + t1092 + t538) * t539 + (-t1040 * t412 + t1078 + t392 - t395 - t398) * t454 + (t1040 * t414 - t1077 + t397 + t400 - t889) * t453;
t975 = (-Icges(7,3) * t867 + t1091 + t1093 + t537) * t539 + (t1041 * t411 - t1076 - t393 + t396 + t399) * t454 + (-t1041 * t413 + t1075 + t394 - t892 - t904) * t453;
t974 = t1039 * t556;
t970 = t1042 * t539 + t760 * t978;
t870 = t550 * t559;
t886 = Icges(4,3) * t551;
t355 = Icges(4,5) * t870 - Icges(4,6) * t871 - t886;
t516 = Icges(4,4) * t871;
t899 = Icges(4,5) * t551;
t359 = Icges(4,1) * t870 - t516 - t899;
t894 = Icges(4,6) * t551;
t357 = Icges(4,4) * t870 - Icges(4,2) * t871 - t894;
t877 = t357 * t556;
t668 = -t359 * t559 + t877;
t138 = -t355 * t551 - t550 * t668;
t618 = t555 * t814 + t556 * t811;
t498 = Icges(4,5) * t559 - Icges(4,6) * t556;
t497 = Icges(4,5) * t556 + Icges(4,6) * t559;
t628 = qJD(3) * t497;
t905 = Icges(4,4) * t556;
t502 = Icges(4,1) * t559 - t905;
t360 = Icges(4,5) * t550 + t502 * t551;
t358 = Icges(4,6) * t550 + t551 * t685;
t876 = t358 * t556;
t667 = -t360 * t559 + t876;
t969 = -t551 * t628 + (-t498 * t550 + t667 + t886) * qJD(1);
t356 = Icges(4,3) * t550 + t498 * t551;
t968 = -t550 * t628 + (t356 + t668) * qJD(1);
t499 = Icges(4,2) * t559 + t905;
t662 = t499 * t556 - t501 * t559;
t967 = qJD(1) * t662 + t498 * qJD(3);
t859 = t198 * qJ(5) + t809;
t118 = pkin(4) * t199 + t859;
t807 = qJD(5) * t556;
t535 = t558 * t807;
t808 = qJD(5) * t555;
t536 = t559 * t808;
t966 = qJD(3) * t536 + qJD(4) * t535 + t453 * t118 + t352 * t290;
t841 = -Icges(4,2) * t870 + t359 - t516;
t843 = t501 * t550 + t357;
t965 = -t556 * t841 - t559 * t843;
t961 = t351 / 0.2e1;
t960 = t352 / 0.2e1;
t959 = -t453 / 0.2e1;
t958 = t453 / 0.2e1;
t957 = -t454 / 0.2e1;
t956 = t454 / 0.2e1;
t955 = -t539 / 0.2e1;
t954 = t539 / 0.2e1;
t950 = -rSges(6,1) - pkin(8);
t949 = -rSges(5,3) - pkin(8);
t553 = t560 * pkin(1);
t945 = rSges(4,1) * t559;
t944 = rSges(6,1) * t556;
t938 = rSges(5,3) * t556;
t145 = t411 * t453 + t413 * t454;
t806 = qJD(6) * t558;
t533 = t556 * t806;
t534 = t555 * t807;
t763 = t448 * t817 + t450 * t816 + qJD(2);
t688 = t453 * t290 + t534 + t763;
t793 = t296 + t855;
t41 = t453 * t858 + t454 * t793 + t533 + t688;
t935 = t145 * t41;
t934 = t25 * t454;
t933 = t26 * t453;
t932 = t27 * t454;
t931 = t28 * t453;
t930 = t29 * t454;
t929 = t30 * t453;
t581 = -t454 * t996 + t388 + t585;
t794 = -t290 - t858;
t54 = t1023 * t454 + t539 * t794 + t385 + t581;
t922 = t54 * t551;
t543 = t550 * rSges(4,3);
t921 = t551 * t90;
t712 = rSges(6,2) * t412 - rSges(6,3) * t411;
t232 = rSges(6,1) * t871 - t712;
t854 = t235 + t296;
t56 = t232 * t453 + t454 * t854 + t688;
t920 = t56 * t232;
t856 = -t232 - t290;
t65 = t454 * t599 + t539 * t856 + t581;
t919 = t65 * t599;
t918 = t81 * t351;
t917 = t82 * t352;
t916 = t83 * t351;
t915 = t84 * t352;
t914 = t85 * t351;
t913 = t86 * t352;
t909 = t118 * t869 + t290 * t774;
t825 = rSges(4,2) * t871 + t551 * rSges(4,3);
t361 = rSges(4,1) * t870 - t825;
t510 = rSges(4,1) * t556 + rSges(4,2) * t559;
t776 = t510 * t816;
t169 = -t776 + (-t361 + t748) * qJD(1);
t881 = t169 * t550;
t880 = t169 * t551;
t362 = rSges(4,1) * t868 - rSges(4,2) * t869 + t543;
t780 = t510 * t817;
t170 = qJD(1) * t362 - t780 + t830;
t440 = t510 * t551;
t879 = t170 * t440;
t489 = qJD(3) * t523;
t874 = t489 * t551;
t873 = t497 * t550;
t872 = t497 * t551;
t795 = rSges(6,1) * t774 + t197 * rSges(6,2) - t196 * rSges(6,3);
t120 = -rSges(6,1) * t783 + t795;
t863 = t117 + t120;
t862 = -t1036 * t783 + t987;
t861 = t1035 * t199 + t1036 * t623 - t988;
t288 = -pkin(4) * t411 + qJ(5) * t412;
t860 = t453 * t288 + t535;
t294 = -pkin(4) * t413 + qJ(5) * t414;
t853 = qJD(5) * t412 + t539 * t294;
t852 = t296 * t815 + t783 * t996;
t851 = t559 * t290 + t871 * t996;
t619 = -t555 * t812 + t558 * t814;
t301 = pkin(4) * t619 + qJ(5) * t618 + t534;
t444 = -t559 * t711 + t944;
t470 = (rSges(6,2) * t555 + rSges(6,3) * t558) * t556;
t313 = qJD(3) * t444 + qJD(4) * t470;
t850 = -t301 - t313;
t442 = t559 * t714 + t938;
t471 = (-rSges(5,1) * t555 - rSges(5,2) * t558) * t556;
t311 = qJD(3) * t442 + qJD(4) * t471;
t849 = -t311 - t489;
t443 = rSges(7,1) * t556 + t559 * t710;
t473 = (rSges(7,2) * t558 - rSges(7,3) * t555) * t556;
t848 = pkin(5) * t815 + qJ(6) * t619 + qJD(3) * t443 + qJD(4) * t473 + t533;
t847 = -t550 * t355 - t359 * t868;
t846 = t550 * t356 + t360 * t868;
t845 = t1023 * t550;
t844 = t1023 * t551;
t842 = -t501 * t551 - t358;
t840 = -t499 * t551 + t360;
t447 = t522 * t550;
t449 = t522 * t551;
t839 = -t447 * t817 - t449 * t816;
t837 = t550 * t448 + t551 * t450;
t469 = (-pkin(4) * t555 + qJ(5) * t558) * t556;
t836 = qJD(5) * t414 - t454 * t469;
t834 = -t441 - t522;
t833 = -pkin(5) * t556 - qJ(6) * t864 - t443;
t831 = t599 - t996;
t829 = rSges(4,2) * t783 + rSges(4,3) * t819;
t494 = t550 * pkin(3) * t815;
t828 = t494 - t549;
t827 = -t499 + t502;
t826 = t501 + t685;
t822 = qJD(1) * t447;
t821 = qJD(1) * t498;
t245 = -t550 * t662 - t872;
t805 = t245 * qJD(1);
t804 = qJD(1) * qJD(3);
t801 = -pkin(4) - t1035;
t800 = t561 * t948;
t799 = t561 * t553;
t798 = t117 + t862;
t797 = t454 * t290;
t796 = -t197 * rSges(5,1) + t196 * rSges(5,2) + rSges(5,3) * t774;
t622 = -t775 - t782;
t298 = pkin(3) * t622 - pkin(8) * t783 + t495;
t299 = pkin(8) * t623 + qJD(1) * t526 - t494;
t792 = t299 * t817 + (t298 + t983) * t816;
t791 = t551 * t298 + t550 * t299 + t448 * t819;
t790 = -t301 - t848;
t789 = -t489 + t850;
t238 = t414 * rSges(5,1) - t413 * rSges(5,2) + rSges(5,3) * t869;
t788 = -t996 + t1023;
t787 = -t522 + t831;
t785 = t553 + t824;
t770 = qJD(6) * t867;
t769 = t559 * t806;
t762 = -pkin(2) - t945;
t761 = t550 * t804;
t757 = t819 / 0.2e1;
t756 = -t817 / 0.2e1;
t753 = t816 / 0.2e1;
t752 = t815 / 0.2e1;
t747 = t547 - t948;
t746 = t41 * t858;
t745 = t54 * t1023;
t744 = t66 * t831;
t314 = t360 * t870;
t741 = t356 * t551 - t314;
t740 = -t355 + t876;
t738 = t559 * t118 + t301 * t871 + t623 * t996;
t737 = t550 * t290 + t551 * t296 + t837;
t736 = -t489 + t790;
t734 = -t522 + t788;
t727 = qJD(1) * (-pkin(2) * t820 + t532) - t800;
t726 = -t391 + t747;
t725 = t55 * t788;
t722 = qJD(4) * t752;
t409 = qJD(1) * t449;
t721 = -t523 * t817 - t409;
t713 = -rSges(6,2) * t199 + rSges(6,3) * t198;
t122 = rSges(6,1) * t623 + t713;
t610 = -t450 * t761 + t792;
t8 = t122 * t453 + t232 * t352 - t351 * t854 + t454 * t863 + t610 + t966;
t720 = t56 * t122 + t8 * t232;
t475 = rSges(3,1) * t550 + rSges(3,2) * t551;
t717 = -rSges(4,2) * t556 + t945;
t716 = rSges(5,1) * t199 - rSges(5,2) * t198;
t690 = t551 * t290 * t810 - t453 * t370 + t536 + t839;
t689 = t785 + t450;
t676 = -t170 * t550 - t880;
t669 = t236 * t551 - t238 * t550;
t183 = t357 * t559 + t359 * t556;
t184 = t358 * t559 + t360 * t556;
t666 = t361 * t550 + t362 * t551;
t661 = qJD(1) * t298 + t727;
t660 = -pkin(2) - t523;
t659 = t551 * t117 + t550 * t118 + t290 * t819 + t791;
t640 = -t489 * t550 - t522 * t819;
t439 = t510 * t550;
t637 = -t56 * t854 - t919;
t636 = t66 * t235 + t65 * t856;
t635 = -t523 * t816 + t822;
t632 = t556 * t950 + t784;
t631 = t556 * t949 + t784;
t630 = qJD(3) * t501;
t629 = qJD(3) * t499;
t139 = -t358 * t871 - t741;
t627 = (-t138 * t551 + t139 * t550) * qJD(3);
t140 = -t357 * t869 - t847;
t141 = -t358 * t869 + t846;
t626 = (-t140 * t551 + t141 * t550) * qJD(3);
t7 = -qJD(4) * t770 + t861 * t453 + t858 * t352 + (-t450 * t820 + t769) * qJD(3) + t798 * t454 - t793 * t351 + t792 + t966;
t611 = t41 * t861 + t7 * t858;
t609 = qJD(5) * t198 + t539 * t117 + t296 * t760 + t661;
t608 = -t556 * t840 + t559 * t842;
t607 = t522 * t761 + (-t299 - t455) * qJD(1) - t799;
t606 = t296 + t689;
t605 = -t41 * t793 - t745;
t604 = t54 * t794 + t55 * t855;
t587 = (-t556 * t826 + t559 * t827) * qJD(1);
t586 = -t489 + (-t806 - t808) * t556;
t582 = -qJD(5) * t196 - t454 * t301 + t351 * t996 + t607;
t248 = rSges(4,1) * t622 - rSges(4,2) * t774 + t829;
t249 = -qJD(3) * t439 + (t551 * t717 + t543) * qJD(1);
t574 = t248 * t551 + t249 * t550 + (t361 * t551 - t362 * t550) * qJD(1);
t79 = t236 * t453 + t238 * t454 + t763;
t91 = t238 * t539 - t441 * t453 + t735;
t570 = t79 * t669 + (-t551 * t91 + t1063) * t441;
t479 = t685 * qJD(3);
t480 = t502 * qJD(3);
t567 = qJD(1) * t497 - t479 * t556 + t480 * t559 + (-t499 * t559 - t501 * t556) * qJD(3);
t566 = (t744 + t920) * t551 + t637 * t550;
t565 = (t725 + t746) * t551 + t605 * t550;
t481 = t717 * qJD(3);
t340 = t599 * t551;
t338 = t599 * t550;
t336 = t441 * t551;
t335 = t441 * t550;
t297 = rSges(7,2) * t414 - rSges(7,3) * t413;
t295 = -rSges(5,1) * t413 - rSges(5,2) * t414;
t293 = rSges(6,2) * t413 + rSges(6,3) * t414;
t291 = rSges(7,2) * t412 - rSges(7,3) * t411;
t289 = -rSges(5,1) * t411 - rSges(5,2) * t412;
t287 = rSges(6,2) * t411 + rSges(6,3) * t412;
t247 = t290 * t869;
t246 = -t551 * t662 + t873;
t200 = t246 * qJD(1);
t168 = qJD(3) * t666 + qJD(2);
t129 = -t799 - t481 * t816 + (-t249 - t455 + t780) * qJD(1);
t128 = -t481 * t817 + (t248 - t776) * qJD(1) + t727;
t124 = rSges(5,3) * t623 + t716;
t123 = -rSges(5,3) * t783 + t796;
t93 = t567 * t550 - t551 * t967;
t92 = t550 * t967 + t567 * t551;
t88 = -qJD(3) * t667 + (-t551 * t629 + (-t550 * t685 + t894) * qJD(1)) * t559 + (-t551 * t630 + (-t502 * t550 + t899) * qJD(1)) * t556;
t87 = -qJD(3) * t668 + (qJD(1) * t358 - t550 * t629) * t559 + (qJD(1) * t360 - t550 * t630) * t556;
t80 = t574 * qJD(3);
t61 = t200 + t626;
t60 = t627 + t805;
t43 = -t124 * t539 - t311 * t454 + t351 * t441 + (-t236 * t812 - t874) * qJD(3) + t607;
t42 = t123 * t539 - t311 * t453 - t352 * t441 + (t238 * t812 + t640) * qJD(3) + t661;
t31 = t123 * t454 + t124 * t453 + t236 * t352 - t238 * t351 + t610;
t24 = -t313 * t454 - t351 * t599 + (-t118 - t122) * t539 + (t812 * t856 - t874) * qJD(3) + t582;
t23 = t120 * t539 + t850 * t453 + t831 * t352 + (t235 * t812 + t640) * qJD(3) + t609;
t10 = -qJD(6) * t197 - t848 * t454 - t1023 * t351 + (-t118 - t861) * t539 + (t794 * t812 - t874) * qJD(3) + t582;
t9 = qJD(6) * t199 + t862 * t539 + t790 * t453 + t788 * t352 + (t812 * t855 + t640) * qJD(3) + t609;
t1 = [(t200 + ((t139 - t314 + (t356 + t877) * t551 + t847) * t551 + t846 * t550) * qJD(3)) * t753 + (-t56 * t797 - (-t290 * t56 + t744) * t454 + t23 * (t606 + t235) + (t550 * t632 + t712 + t726 - t946) * t24 + (-t118 - t713 + t828 + t950 * t1062 + (t551 * t632 - t1064) * qJD(1)) * t65 + (t232 * t539 + t65 + t795 + (-t948 + (t660 - t944) * t550) * qJD(1) + t1059) * t66) * m(6) - (t87 + t93 + t61) * t816 / 0.2e1 + (t88 + t92) * t817 / 0.2e1 + (-t805 + ((t551 * t740 + t141 - t846) * t551 + (t550 * t740 + t140 + t741) * t550) * qJD(3) + t60) * t756 + (-qJD(3) * t662 + t479 * t559 + t480 * t556) * qJD(1) + t1057 * t958 + (t1058 + t981) * t957 + t981 * t956 + t979 * t960 + t980 * t961 + ((t184 + t246) * t551 + (t183 + t245) * t550) * t804 / 0.2e1 + t916 / 0.2e1 + t917 / 0.2e1 + t914 / 0.2e1 + t915 / 0.2e1 + t913 / 0.2e1 + t931 / 0.2e1 - t932 / 0.2e1 + t933 / 0.2e1 + t929 / 0.2e1 - t930 / 0.2e1 - t934 / 0.2e1 + (t129 * (t550 * t762 + t747 + t825) - t169 * t549 + t128 * (t362 + t785) + t170 * (t532 + t829) + (t510 * t881 - t879) * qJD(3) + ((-pkin(2) - t717) * t880 + (t169 * (-rSges(4,3) - pkin(7)) + t170 * t762) * t550) * qJD(1) - (-qJD(1) * t361 - t169 - t458 - t776) * t170) * m(4) + (t43 * (-t715 + t747) + t90 * (-t716 + t828) + t42 * (t689 + t238) + (t814 * t90 * t949 + t43 * t631) * t550 + (-pkin(7) * t1063 + t631 * t921) * qJD(1) + (t90 - t999 + t796 + (-t948 + (t660 - t938) * t550) * qJD(1) + t1067) * t91) * m(5) + t918 / 0.2e1 + (t9 * (t606 + t855) + t624 * t922 * qJD(1) + (t412 * t801 + t726 - t939 + t997) * t10 - t41 * t797 - (-t290 * t41 + t725) * t454 + ((-t948 + t997) * qJD(1) + t987 - t385 + t858 * t539 + t1059) * t55 + (t828 - t859 + t801 * t199 + (-pkin(7) * qJD(1) + t802 * t814) * t550 + t988 + t55) * t54) * m(7) + t970 + m(3) * ((-t475 * t561 - t800) * (t553 + t742) + (-t561 * t742 - t799) * (-t475 - t948)); m(4) * t80 + m(5) * t31 + m(6) * t8 + m(7) * t7; -qJD(1) * ((t556 * t827 + t559 * t826) * qJD(1) + ((t550 * t840 - t551 * t841) * t559 + (t550 * t842 + t551 * t843) * t556) * qJD(3)) / 0.2e1 + qJD(1) * (t550 * t88 - t551 * t87 + (t183 * t550 + t184 * t551) * qJD(1)) / 0.2e1 + ((-t816 * t873 - t821) * t551 + (t587 + (t608 * t550 + (t872 - t965) * t551) * qJD(3)) * t550) * t753 + ((-t817 * t872 + t821) * t550 + (t587 + (-t965 * t551 + (t873 + t608) * t550) * qJD(3)) * t551) * t756 + (-t586 * t922 - (t556 * t604 + t559 * t565) * qJD(4) - (t690 + t769 + (-t371 + t844) * t454 + t845 * t453) * t41 - (-t409 + t586 * t550 + t844 * t539 + (-t474 + t833) * t453 + t990) * t55 + t7 * t737 + t41 * t659 + (t10 * t734 + t7 * t855 + t41 * t862 + (t55 * t734 + t746) * qJD(1)) * t551 + (t9 * t734 + t55 * t736 + (-t745 + t41 * (-t450 - t793)) * qJD(1) + t611) * t550 + (-t822 - (t370 - t845) * t539 - t833 * t454 + t736 * t551 + t1021) * t54) * m(7) + (t8 * t737 + t56 * t659 + (t24 * t787 + t8 * t235 + t56 * t120 + (t66 * t787 + t920) * qJD(1)) * t551 + (t23 * t787 + t66 * t789 + (-t919 + t56 * (-t450 - t854)) * qJD(1) + t720) * t550 - t66 * (t340 * t539 - t534 * t550 + t721 + (-t444 - t474) * t453 + t990) - t56 * (t338 * t453 + t690 + (t340 - t371) * t454) - (t556 * t636 + t559 * t566) * qJD(4) + (t444 * t454 - t635 - (-t338 + t370) * t539 + (t789 + t534) * t551 + t1021) * t65) * m(6) + (t90 * t457 + t31 * t837 + (t31 * t238 + t90 * t849 + (qJD(1) * t91 + t43) * t834) * t551 + (qJD(1) * t441 * t90 + t31 * t236 + t42 * t834 + t91 * t849) * t550 - t90 * (t335 * t539 - t442 * t454 + t635) - t91 * (-t336 * t539 - t442 * t453 + t721) - ((-t236 * t90 + t238 * t91) * t556 + t570 * t559) * qJD(4) + (t791 + (qJD(1) * t236 + t123) * t551 + (t124 + (-t238 - t450) * qJD(1)) * t550 + t335 * t453 + t336 * t454 - t839) * t79) * m(5) + (t80 * t666 + t168 * t574 + t676 * t481 + (-t128 * t550 - t129 * t551 + (-t170 * t551 + t881) * qJD(1)) * t510 - (t169 * t439 - t879) * qJD(1) - (t168 * (-t439 * t550 - t440 * t551) + t676 * t717) * qJD(3)) * m(4) + t1000 * t961 + t1001 * t960 + ((t1054 * t870 + t979 * t556) * qJD(4) + ((qJD(4) * t1053 + t1003) * t559 + t1004) * t551 + (t1005 * t414 + t1006 * t413) * t539 + (t1008 * t414 - t1010 * t413) * t454 + (-t1007 * t414 + t1009 * t413) * t453) * t959 + (qJD(1) * t972 + t1016 * t550 - t1017 * t551) * t958 + (qJD(1) * t971 + t1014 * t550 - t1015 * t551) * t957 + ((t1055 * t868 + t980 * t556) * qJD(4) + ((qJD(4) * t1056 + t1003) * t559 + t1004) * t550 + (t1005 * t412 + t1006 * t411) * t539 + (t1008 * t412 - t1010 * t411) * t454 + (-t1007 * t412 + t1009 * t411) * t453) * t956 + ((qJD(4) * t973 - t1045) * t559 + ((t1005 * t558 + t1006 * t555 + t1048) * t539 + (t1008 * t558 - t1010 * t555 - t1074) * t454 + (-t1007 * t558 + t1009 * t555 + t1049) * t453 + t978 * qJD(4)) * t556) * t955 + (qJD(1) * t973 + t1012 * t550 - t1013 * t551) * t954 - t1011 * t812 / 0.2e1 + t1002 * t722 + (qJD(1) * t92 + (-t968 * t1044 + t550 ^ 2 * t969 + (t140 * t550 + t141 * t551) * qJD(1)) * t1037 + t1019) * t550 / 0.2e1 - (qJD(1) * t93 + (t551 ^ 2 * t968 - t969 * t1044 + (t138 * t550 + t139 * t551) * qJD(1)) * t1037 + t1018) * t551 / 0.2e1 + (t627 + t60 + t982) * t820 / 0.2e1 + (t626 + t61 + t981) * t757 - t1043 * t810 / 0.2e1; (-t54 * (-qJD(6) * t413 - t454 * t473 + t836 + (-t288 - t291) * t539) - t55 * (-qJD(6) * t411 + t297 * t539 + t853 + (-t469 - t473) * t453) - t41 * (t291 * t453 - t770 + t860 + (t294 + t297) * t454) - (t54 * t651 + t55 * t650 - t935) * qJ(6) + t10 * t851 + t54 * t738 + t55 * t852 + t7 * t247 + t41 * t909 + (qJD(3) * t565 + t10 * t858 + t54 * t861 - t55 * t798 - t793 * t9) * t559 + (t604 * qJD(3) + (qJD(1) * t605 + t55 * t790 + t788 * t9 + t611) * t551 + (-t10 * t1023 + t54 * t848 - t7 * t793 - t41 * t798 + (-t1023 * t55 + t41 * t794) * qJD(1)) * t550) * t556) * m(7) + (-t65 * (-t454 * t470 + (-t287 - t288) * t539 + t836) - t66 * (t293 * t539 + (-t469 - t470) * t453 + t853) - t56 * (t287 * t453 + (t293 + t294) * t454 + t860) + t24 * t851 + t65 * t738 + t66 * t852 + t8 * t247 + t56 * t909 + (qJD(3) * t566 + t65 * t122 - t23 * t854 + t24 * t232 - t66 * t863) * t559 + (t636 * qJD(3) + (qJD(1) * t637 + t23 * t831 + t66 * t850 + t720) * t551 + (-t24 * t599 + t65 * t313 - t8 * t854 - t56 * t863 + (t56 * t856 - t599 * t66) * qJD(1)) * t550) * t556) * m(6) + ((qJD(3) * t570 - t91 * t123 + t90 * t124 + t43 * t236 - t42 * t238) * t559 + (t90 * (-qJD(3) * t236 + t311 * t550) + t91 * (qJD(3) * t238 - t311 * t551) + t31 * t669 + t79 * (-t123 * t550 + t124 * t551 - t236 * t820 - t238 * t819) + (-t42 * t551 + t43 * t550 + (t550 * t91 + t921) * qJD(1)) * t441) * t556 - t90 * (-t289 * t539 - t454 * t471) - t91 * (t295 * t539 - t453 * t471) - t79 * (t289 * t453 + t295 * t454)) * m(5) + (t556 * t971 - t559 * t980) * t961 + (t556 * t972 - t559 * t979) * t960 + (t413 * t976 + t414 * t975 - t551 * t974) * t959 + ((qJD(3) * t972 - t1057) * t559 + (-qJD(1) * t1001 + t979 * qJD(3) + t1016 * t551 + t1017 * t550) * t556) * t958 + ((qJD(3) * t971 - t1058) * t559 + (-qJD(1) * t1000 + t980 * qJD(3) + t1014 * t551 + t1015 * t550) * t556) * t957 + (t411 * t976 + t412 * t975 - t550 * t974) * t956 + (t1039 * t559 + (t555 * t976 + t558 * t975) * t556) * t955 + ((qJD(3) * t973 - t1042) * t559 + (-qJD(1) * t1002 + t978 * qJD(3) + t1012 * t551 + t1013 * t550) * t556) * t954 - (t917 + t918 + t933 - t934 + t915 + t916 + t931 - t932 + t913 + t914 + t929 - t930 + t970) * t559 / 0.2e1 + t1018 * t871 / 0.2e1 + t1019 * t869 / 0.2e1 + t1011 * t752 - t981 * t783 / 0.2e1 + t982 * t556 * t757 + (t556 * t973 - t559 * t978) * t722 + t1043 * t814 / 0.2e1; (t10 * t413 + t41 * t618 + t411 * t9 + t54 * t992 + t55 * t991 + t7 * t867 - t935) * m(7) + (t8 * t867 + t23 * t411 + t24 * t413 + t991 * t66 + t992 * t65 + (-t145 + t618) * t56) * m(6); (t10 * t414 + t412 * t9 + t7 * t865 + (-t414 * t539 + t453 * t865 + t199) * t55 + (t412 * t539 + t454 * t865 - t197) * t54 + (-t412 * t453 - t414 * t454 + t619) * t41) * m(7);];
tauc  = t1(:);