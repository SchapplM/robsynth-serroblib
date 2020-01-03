% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPR7_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR7_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR7_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:09
% EndTime: 2019-12-31 21:17:00
% DurationCPUTime: 39.28s
% Computational Cost: add. (130283->946), mult. (130006->1250), div. (0->0), fcn. (139290->10), ass. (0->565)
t780 = qJ(2) + qJ(3);
t765 = sin(t780);
t1067 = qJ(4) * t765;
t766 = cos(t780);
t1103 = pkin(3) * t766;
t783 = -pkin(8) - qJ(4);
t995 = qJ(4) + t783;
t1205 = t995 * t765;
t785 = sin(qJ(1));
t1037 = t765 * t785;
t777 = pkin(9) + qJ(5);
t763 = sin(t777);
t1014 = t785 * t763;
t764 = cos(t777);
t787 = cos(qJ(1));
t670 = t766 * t1014 + t764 * t787;
t1013 = t785 * t764;
t671 = t766 * t1013 - t763 * t787;
t1201 = -t671 * rSges(6,1) + t670 * rSges(6,2);
t500 = -rSges(6,3) * t1037 + t1201;
t1036 = t765 * t787;
t1027 = t766 * t787;
t672 = -t763 * t1027 + t1013;
t673 = t764 * t1027 + t1014;
t861 = t673 * rSges(6,1) + t672 * rSges(6,2);
t501 = rSges(6,3) * t1036 + t861;
t781 = sin(pkin(9));
t1012 = t785 * t781;
t874 = pkin(4) * t1012 - t783 * t1036;
t1026 = t781 * t787;
t1028 = t766 * t785;
t782 = cos(pkin(9));
t759 = pkin(4) * t782 + pkin(3);
t875 = pkin(4) * t1026 - t759 * t1028;
t1102 = pkin(3) - t759;
t890 = t1102 * t766;
t724 = t1067 + t1103;
t778 = t785 ^ 2;
t779 = t787 ^ 2;
t934 = t778 + t779;
t948 = t934 * t724;
t285 = t948 + (t501 + t874 + (-t890 - t1067) * t787) * t787 + (-t500 - (t1103 + t1205) * t785 - t875) * t785;
t1029 = t766 * t783;
t1199 = t1102 * t765;
t1094 = rSges(6,3) * t766;
t1046 = t763 * t765;
t927 = rSges(6,2) * t1046;
t831 = t927 + t1094;
t1209 = t785 * t831;
t1045 = t764 * t765;
t929 = rSges(6,1) * t1045;
t715 = t785 * t929;
t571 = -t715 + t1209;
t940 = rSges(6,3) * t1027 + t787 * t927;
t572 = -t787 * t929 + t940;
t741 = qJ(4) * t1027;
t750 = pkin(3) * t1037;
t869 = qJ(4) * t1028 - t750;
t941 = t783 * t1028 + t759 * t1037;
t1194 = (-t741 + t572 + (-t1029 + t1199) * t787) * t787 + (-t869 - t941 + t571) * t785;
t949 = t785 * t869 + t787 * (-pkin(3) * t1036 + t741);
t317 = t949 + t1194;
t1097 = rSges(6,1) * t764;
t860 = -rSges(6,2) * t763 + t1097;
t612 = t765 * t860 - t1094;
t1203 = t995 * t766 - t1199 + t612;
t722 = pkin(3) * t765 - qJ(4) * t766;
t907 = t722 + t1203;
t438 = t907 * t785;
t1095 = rSges(6,3) * t765;
t613 = t766 * t860 + t1095;
t1202 = t890 + t1205 - t613;
t906 = -t724 + t1202;
t439 = t906 * t785;
t440 = t907 * t787;
t441 = t906 * t787;
t110 = t285 * t317 - t438 * t439 - t440 * t441;
t1025 = t782 * t787;
t700 = t766 * t1012 + t1025;
t1011 = t785 * t782;
t701 = t766 * t1011 - t1026;
t1200 = -t701 * rSges(5,1) + t700 * rSges(5,2);
t702 = -t766 * t1026 + t1011;
t703 = t766 * t1025 + t1012;
t864 = t703 * rSges(5,1) + t702 * rSges(5,2);
t371 = t785 * (rSges(5,3) * t1037 - t1200) + t787 * (rSges(5,3) * t1036 + t864) + t948;
t1096 = rSges(5,3) * t766;
t1038 = t765 * t782;
t930 = rSges(5,1) * t1038;
t728 = t785 * t930;
t1039 = t765 * t781;
t928 = rSges(5,2) * t1039;
t937 = rSges(5,3) * t1027 + t787 * t928;
t1204 = t785 * (-t728 + (t928 + t1096) * t785) + t787 * (-t787 * t930 + t937);
t388 = t949 + t1204;
t1098 = rSges(5,1) * t782;
t863 = -rSges(5,2) * t781 + t1098;
t632 = t765 * t863 - t1096;
t958 = t632 + t722;
t545 = t958 * t785;
t633 = rSges(5,3) * t765 + t766 * t863;
t957 = -t633 - t724;
t546 = t957 * t785;
t547 = t958 * t787;
t548 = t957 * t787;
t196 = t371 * t388 - t545 * t546 - t547 * t548;
t659 = rSges(4,1) * t1028 - rSges(4,2) * t1037 - t787 * rSges(4,3);
t886 = -rSges(4,2) * t1036 + t785 * rSges(4,3);
t511 = t785 * t659 + t787 * (rSges(4,1) * t1027 + t886);
t723 = rSges(4,1) * t765 + rSges(4,2) * t766;
t697 = t723 * t785;
t699 = t723 * t787;
t552 = -t785 * t697 - t787 * t699;
t1099 = rSges(4,1) * t766;
t725 = -rSges(4,2) * t765 + t1099;
t363 = t934 * t723 * t725 + t511 * t552;
t1262 = m(4) * t363 + m(5) * t196 + m(6) * t110;
t786 = cos(qJ(2));
t1104 = pkin(2) * t786;
t760 = pkin(1) + t1104;
t1175 = -pkin(7) - pkin(6);
t761 = t785 * t1175;
t776 = t787 * pkin(6);
t936 = -t787 * t1175 - t785 * t760;
t959 = -t785 * (pkin(1) * t785 - t776 + t936) + t787 * (-t785 * pkin(6) - t761 + (-pkin(1) + t760) * t787);
t251 = t285 + t959;
t784 = sin(qJ(2));
t1105 = pkin(2) * t784;
t889 = t722 + t1105;
t835 = t889 + t1203;
t426 = t835 * t785;
t428 = t835 * t787;
t103 = t251 * t317 - t426 * t439 - t428 * t441;
t324 = t371 + t959;
t873 = t632 + t889;
t513 = t873 * t785;
t515 = t873 * t787;
t169 = t324 * t388 - t513 * t546 - t515 * t548;
t818 = t723 + t1105;
t1207 = t818 * t787;
t1208 = t818 * t785;
t401 = t511 + t959;
t299 = t401 * t552 + (t1207 * t787 + t1208 * t785) * t725;
t1261 = m(4) * t299 + m(5) * t169 + m(6) * t103;
t1162 = -t766 / 0.2e1;
t1243 = t765 / 0.2e1;
t492 = Icges(6,5) * t673 + Icges(6,6) * t672 + Icges(6,3) * t1036;
t1076 = Icges(6,4) * t764;
t854 = -Icges(6,2) * t763 + t1076;
t606 = -Icges(6,6) * t766 + t765 * t854;
t568 = t606 * t787;
t1077 = Icges(6,4) * t763;
t856 = Icges(6,1) * t764 - t1077;
t608 = -Icges(6,5) * t766 + t765 * t856;
t570 = t608 * t787;
t849 = Icges(6,5) * t764 - Icges(6,6) * t763;
t604 = -Icges(6,3) * t766 + t765 * t849;
t1078 = Icges(6,4) * t673;
t495 = Icges(6,2) * t672 + Icges(6,6) * t1036 + t1078;
t646 = Icges(6,4) * t672;
t498 = Icges(6,1) * t673 + Icges(6,5) * t1036 + t646;
t844 = -t495 * t763 + t498 * t764;
t826 = -t604 * t787 - t844;
t246 = -t826 * t766 + (t568 * t763 - t570 * t764 + t492) * t765;
t607 = Icges(6,6) * t765 + t766 * t854;
t609 = Icges(6,5) * t765 + t766 * t856;
t1055 = t604 * t766;
t605 = Icges(6,3) * t765 + t766 * t849;
t1053 = t608 * t764;
t1054 = t606 * t763;
t839 = t1053 - t1054;
t823 = t605 - t839;
t799 = t765 * t823 + t1055;
t293 = t672 * t607 + t673 * t609 + t787 * t799;
t306 = -t823 * t766 + (-t607 * t763 + t609 * t764 + t604) * t765;
t490 = Icges(6,5) * t671 - Icges(6,6) * t670 + Icges(6,3) * t1037;
t1057 = t490 * t766;
t645 = Icges(6,4) * t671;
t493 = -Icges(6,2) * t670 + Icges(6,6) * t1037 + t645;
t644 = Icges(6,4) * t670;
t497 = -Icges(6,1) * t671 - Icges(6,5) * t1037 + t644;
t1245 = t493 * t763 + t497 * t764;
t331 = t1245 * t765 + t1057;
t1056 = t492 * t766;
t332 = t765 * t844 - t1056;
t375 = t604 * t1037 - t606 * t670 + t608 * t671;
t377 = t604 * t1036 + t672 * t606 + t673 * t608;
t390 = t765 * t839 - t1055;
t1260 = t306 * t1162 + t390 * t1243 + (-t331 + t375) * t1028 / 0.4e1 + (t332 + t377) * t1027 / 0.4e1 + t1036 * (t246 + t293) / 0.4e1;
t1177 = m(6) / 0.2e1;
t1179 = m(5) / 0.2e1;
t905 = t715 + t941;
t459 = (-t831 + t1105) * t785 + t905;
t815 = -t1029 + (-t759 - t1097) * t765;
t460 = (t815 - t1105) * t787 + t940;
t1087 = rSges(5,3) + qJ(4);
t817 = -t1087 * t766 - t928;
t938 = t728 + t750;
t484 = (t817 + t1105) * t785 + t938;
t866 = (-pkin(3) - t1098) * t765;
t904 = t741 + t937;
t485 = (t866 - t1105) * t787 + t904;
t1195 = -t1087 * t765 - t1103;
t446 = t1195 * t785 + t1200 + t936;
t447 = -t761 + (-t1195 + t760) * t787 + t864;
t972 = t446 * t1027 + t447 * t1028;
t407 = (-rSges(6,3) + t783) * t1037 + t875 + t936 + t1201;
t408 = -t761 + (t759 * t766 + t1095 + t760) * t787 + t861 + t874;
t978 = t407 * t1027 + t408 * t1028;
t984 = ((t459 * t787 + t460 * t785) * t765 + t978) * t1177 + ((t484 * t787 + t485 * t785) * t765 + t972) * t1179;
t986 = (-t426 * t1036 + t428 * t1037 + t978) * t1177 + (-t513 * t1036 + t515 * t1037 + t972) * t1179;
t26 = t986 - t984;
t1259 = qJD(1) * t26;
t474 = t905 - t1209;
t475 = t787 * t815 + t940;
t521 = t785 * t817 + t938;
t522 = t787 * t866 + t904;
t983 = ((t474 * t787 + t475 * t785) * t765 + t978) * t1177 + ((t521 * t787 + t522 * t785) * t765 + t972) * t1179;
t985 = (-t438 * t1036 + t440 * t1037 + t978) * t1177 + (-t545 * t1036 + t547 * t1037 + t972) * t1179;
t30 = t985 - t983;
t1258 = qJD(1) * t30;
t653 = (-Icges(6,5) * t763 - Icges(6,6) * t764) * t765;
t1032 = t766 * t653;
t540 = -rSges(6,1) * t670 - rSges(6,2) * t671;
t541 = rSges(6,1) * t672 - rSges(6,2) * t673;
t654 = (-Icges(6,2) * t764 - t1077) * t765;
t655 = (-Icges(6,1) * t763 - t1076) * t765;
t173 = (-(t608 / 0.2e1 + t654 / 0.2e1) * t763 + (t655 / 0.2e1 - t606 / 0.2e1) * t764) * t765 + m(6) * (-t407 * t540 + t408 * t541) - t1032 / 0.2e1;
t1257 = t173 * qJD(1);
t309 = t490 * t1037 - t493 * t670 - t497 * t671;
t310 = t492 * t1037 - t670 * t495 + t671 * t498;
t848 = t309 * t785 + t310 * t787;
t145 = -t375 * t766 + t765 * t848;
t194 = -t309 * t787 + t310 * t785;
t311 = t490 * t1036 + t672 * t493 - t673 * t497;
t312 = t492 * t1036 + t672 * t495 + t673 * t498;
t1248 = -t311 * t787 + t312 * t785;
t901 = t1037 / 0.4e1;
t903 = -t1037 / 0.4e1;
t1256 = (t901 + t903) * t1248;
t1160 = t785 / 0.2e1;
t1157 = t787 / 0.2e1;
t524 = Icges(5,5) * t701 - Icges(5,6) * t700 + Icges(5,3) * t1037;
t649 = Icges(4,4) * t1028 - Icges(4,2) * t1037 - Icges(4,6) * t787;
t757 = Icges(4,4) * t766;
t720 = Icges(4,1) * t765 + t757;
t1255 = t720 * t785 - t524 + t649;
t719 = -Icges(4,2) * t765 + t757;
t650 = Icges(4,6) * t785 + t719 * t787;
t1254 = Icges(5,5) * t703 + Icges(5,6) * t702 + Icges(5,3) * t1036 - t720 * t787 - t650;
t850 = Icges(5,5) * t782 - Icges(5,6) * t781;
t626 = -Icges(5,3) * t766 + t765 * t850;
t1079 = Icges(4,4) * t765;
t721 = Icges(4,1) * t766 - t1079;
t1253 = t626 + t721;
t652 = Icges(4,5) * t785 + t721 * t787;
t718 = Icges(4,2) * t766 + t1079;
t1252 = (Icges(5,4) * t703 + Icges(5,2) * t702 + Icges(5,6) * t1036) * t781 - (Icges(5,1) * t703 + Icges(5,4) * t702 + Icges(5,5) * t1036) * t782 - t652 + (-t626 + t718) * t787;
t744 = Icges(4,4) * t1037;
t651 = Icges(4,1) * t1028 - Icges(4,5) * t787 - t744;
t1251 = -(Icges(5,4) * t701 - Icges(5,2) * t700 + Icges(5,6) * t1037) * t781 + (Icges(5,1) * t701 - Icges(5,4) * t700 + Icges(5,5) * t1037) * t782 + t626 * t785 - Icges(4,2) * t1028 + t651 - t744;
t857 = Icges(5,1) * t782 - Icges(5,4) * t781;
t630 = -Icges(5,5) * t766 + t765 * t857;
t1250 = -t630 * t782 - t719 - t720;
t415 = t612 * t1037 - t500 * t766;
t847 = t785 * t311 + t312 * t787;
t1249 = -t377 * t766 + t765 * t847;
t855 = Icges(5,4) * t782 - Icges(5,2) * t781;
t628 = -Icges(5,6) * t766 + t765 * t855;
t1246 = Icges(5,3) * t765 + t628 * t781 + t766 * t850;
t1050 = t649 * t765;
t1158 = -t787 / 0.2e1;
t1230 = t1157 * t1248 + t194 * t1160;
t601 = t652 * t1028;
t717 = Icges(4,5) * t766 - Icges(4,6) * t765;
t1048 = t717 * t787;
t648 = Icges(4,3) * t785 + t1048;
t881 = t787 * t648 - t601;
t431 = -t650 * t1037 - t881;
t647 = Icges(4,5) * t1028 - Icges(4,6) * t1037 - Icges(4,3) * t787;
t964 = -t651 * t1027 - t785 * t647;
t432 = -t649 * t1036 - t964;
t963 = t652 * t1027 + t785 * t648;
t433 = -t650 * t1036 + t963;
t879 = t650 * t765 - t647;
t795 = (-t432 * t787 + t433 * t785) * t1157 + ((t431 - t601 + (t648 + t1050) * t787 + t964) * t787 + t963 * t785 + t1248) * t1158 + (-t194 + (-t524 * t1037 + t433 - t963 + (t647 + t879) * t787) * t787 + (-(t651 * t766 - t1050) * t787 + t431 + t432 + t881 + t524 * t1036 + t785 * t879) * t785) * t1160 + t1230;
t1180 = m(4) / 0.2e1;
t1244 = -t765 / 0.2e1;
t1242 = t766 / 0.2e1;
t1241 = t785 / 0.4e1;
t1240 = -t787 / 0.4e1;
t567 = t606 * t785;
t569 = t608 * t785;
t827 = -t604 * t785 + t1245;
t245 = -t827 * t766 + (t567 * t763 - t569 * t764 + t490) * t765;
t292 = -t607 * t670 + t609 * t671 + t785 * t799;
t1229 = t245 + t292;
t1227 = (-t718 + t1253) * t766 + (t1246 + t1250) * t765;
t1226 = t626 * t1037 - t628 * t700 + t630 * t701 + t649 * t766 + t651 * t765;
t563 = -t659 + t936;
t564 = -t761 + (t760 + t1099) * t787 + t886;
t820 = (-t563 * t787 - t564 * t785) * t725;
t979 = t548 * t446 + t546 * t447;
t980 = t441 * t407 + t439 * t408;
t912 = (-t438 * t460 - t440 * t459 + t980) * t1177 + (-t484 * t547 - t485 * t545 + t979) * t1179 + (t820 + (t1207 * t785 - t1208 * t787) * t723) * t1180;
t913 = (-t426 * t475 - t428 * t474 + t980) * t1177 + (-t513 * t522 - t515 * t521 + t979) * t1179 + (-t1207 * t697 + t1208 * t699 + t820) * t1180;
t1225 = t912 - t913;
t1080 = Icges(3,4) * t784;
t733 = Icges(3,2) * t786 + t1080;
t736 = Icges(3,1) * t786 - t1080;
t1222 = (t736 / 0.2e1 - t733 / 0.2e1) * t784;
t867 = t548 * t1036 + t546 * t1037 - t388 * t766;
t868 = t441 * t1036 + t439 * t1037 - t317 * t766;
t970 = -t515 * t1027 - t513 * t1028;
t908 = t765 * t324 + t970;
t975 = -t428 * t1027 - t426 * t1028;
t909 = t765 * t251 + t975;
t1085 = (t867 + t908) * t1179 + (t868 + t909) * t1177;
t704 = t934 * t765;
t160 = t251 * t704 + t975;
t973 = -t440 * t1027 - t438 * t1028;
t175 = t285 * t704 + t973;
t242 = t324 * t704 + t970;
t965 = -t547 * t1027 - t545 * t1028;
t282 = t371 * t704 + t965;
t1086 = (t282 + t242) * t1179 + (t175 + t160) * t1177;
t1219 = t1085 - t1086;
t1218 = (t1254 * t785 + t1255 * t787) * t766 + (t1251 * t787 + t1252 * t785) * t765;
t392 = (t540 * t787 - t541 * t785) * t765;
t1040 = t765 * t766;
t939 = t934 * t1040;
t1206 = (m(5) / 0.4e1 + m(6) / 0.4e1) * (t939 - t1040);
t773 = Icges(3,4) * t786;
t734 = -Icges(3,2) * t784 + t773;
t735 = Icges(3,1) * t784 + t773;
t534 = -Icges(6,5) * t670 - Icges(6,6) * t671;
t967 = -Icges(6,2) * t671 - t497 - t644;
t969 = -Icges(6,1) * t670 - t493 - t645;
t231 = t534 * t1037 - t967 * t670 + t969 * t671;
t535 = Icges(6,5) * t672 - Icges(6,6) * t673;
t966 = -Icges(6,2) * t673 + t498 + t646;
t968 = Icges(6,1) * t672 - t1078 - t495;
t232 = t535 * t1037 - t966 * t670 + t968 * t671;
t124 = -t231 * t787 + t232 * t785;
t233 = t534 * t1036 + t967 * t672 + t969 * t673;
t234 = t535 * t1036 + t966 * t672 + t968 * t673;
t125 = -t233 * t787 + t234 * t785;
t993 = t124 * t1158 + t125 * t1160;
t1125 = m(6) * (t285 * t765 + t868 + t973);
t1143 = m(5) * (t371 * t765 + t867 + t965);
t994 = t1125 / 0.2e1 + t1143 / 0.2e1;
t263 = -t535 * t766 + (-t966 * t763 + t968 * t764) * t765;
t1063 = t263 * t785;
t262 = -t534 * t766 + (-t967 * t763 + t969 * t764) * t765;
t1064 = t262 * t787;
t961 = t608 + t654;
t962 = -t606 + t655;
t307 = t653 * t1037 - t961 * t670 + t962 * t671;
t308 = t653 * t1036 + t961 * t672 + t962 * t673;
t877 = -t1064 / 0.4e1 + t1063 / 0.4e1 + t308 * t1241 + t307 * t1240;
t1007 = t787 * t1249;
t1021 = t785 * t145;
t1190 = t1007 / 0.4e1 + t1021 / 0.4e1 + t1249 * t1240 - t145 * t1241;
t687 = Icges(3,5) * t785 + t736 * t787;
t943 = -t733 * t787 + t687;
t1010 = t785 * t786;
t1024 = t784 * t785;
t754 = Icges(3,4) * t1024;
t686 = Icges(3,1) * t1010 - Icges(3,5) * t787 - t754;
t944 = -Icges(3,2) * t1010 + t686 - t754;
t685 = Icges(3,6) * t785 + t734 * t787;
t945 = -t735 * t787 - t685;
t684 = Icges(3,4) * t1010 - Icges(3,2) * t1024 - Icges(3,6) * t787;
t946 = t735 * t785 + t684;
t1189 = (-t943 * t785 + t787 * t944) * t784 + (t945 * t785 + t787 * t946) * t786;
t629 = Icges(5,6) * t765 + t766 * t855;
t631 = Icges(5,5) * t765 + t766 * t857;
t793 = -t607 * t1046 / 0.2e1 + t609 * t1045 / 0.2e1 - t629 * t1039 / 0.2e1 + t631 * t1038 / 0.2e1 + t718 * t1244 + (t604 + t1253) * t1243 + (t1053 - t1250) * t1242 + (t1054 + t605 + t1246) * t1162;
t1185 = 0.4e1 * qJD(1);
t1184 = 2 * qJD(2);
t1182 = 2 * qJD(3);
t843 = -t787 * t500 - t501 * t785;
t380 = t843 * t765;
t417 = t612 * t1036 + t766 * t501;
t910 = t380 * t317 + t415 * t441 - t417 * t439;
t316 = t843 * t766 + (t571 * t787 - t572 * t785) * t765;
t352 = (t612 * t785 + t571) * t766 + (t613 * t785 + t500) * t765;
t353 = (-t612 * t787 - t572) * t766 + (-t613 * t787 + t501) * t765;
t911 = t316 * t251 - t352 * t428 - t353 * t426;
t1172 = m(6) * (t910 + t911);
t37 = t285 * t316 - t352 * t440 - t353 * t438 + t910;
t1171 = m(6) * t37;
t409 = t785 * t540 + t541 * t787;
t202 = t251 * t409;
t222 = t285 * t409;
t658 = (-rSges(6,1) * t763 - rSges(6,2) * t764) * t765;
t1168 = m(6) * (t202 + t222 + ((t428 + t440) * t787 + (t426 + t438) * t785) * t658);
t977 = t415 * t1027 - t417 * t1028;
t1165 = m(6) * (-t316 * t766 + (t352 * t787 + t353 * t785 + t380) * t765 + t977);
t1163 = m(6) * (t316 * t380 + t352 * t415 - t353 * t417);
t1161 = -t785 / 0.2e1;
t1100 = rSges(3,1) * t786;
t891 = pkin(1) + t1100;
t935 = rSges(3,2) * t1024 + t787 * rSges(3,3);
t614 = -t785 * t891 + t776 + t935;
t1023 = t784 * t787;
t756 = rSges(3,2) * t1023;
t615 = -t756 + t891 * t787 + (rSges(3,3) + pkin(6)) * t785;
t737 = rSges(3,1) * t784 + rSges(3,2) * t786;
t711 = t737 * t785;
t712 = t737 * t787;
t1156 = m(3) * (t614 * t711 - t615 * t712);
t1149 = m(4) * (-t1207 * t564 + t1208 * t563);
t1148 = m(4) * (t563 * t697 - t564 * t699);
t1133 = m(5) * (t446 * t484 + t447 * t485);
t1132 = m(5) * (t446 * t521 + t447 * t522);
t436 = t447 * t1036;
t1130 = m(5) * (-t446 * t1037 + t436);
t982 = t352 * t407 + t353 * t408;
t1129 = m(6) * (t415 * t459 - t417 * t460 + t982);
t1127 = m(6) * (t415 * t474 - t417 * t475 + t982);
t822 = (-t407 * t787 - t408 * t785) * t658;
t1120 = m(6) * (-t426 * t541 + t428 * t540 + t822);
t1118 = m(6) * (-t438 * t541 + t440 * t540 + t822);
t1112 = m(6) * (t380 * t704 + t977);
t1111 = m(6) * (t407 * t459 + t408 * t460);
t1110 = m(6) * (t407 * t474 + t408 * t475);
t398 = t408 * t1036;
t1109 = m(6) * (-t407 * t1037 + t398);
t1108 = m(6) * (-t417 * t1036 - t415 * t1037);
t1107 = m(6) * (-t409 * t766 - t658 * t704);
t1106 = m(6) * t392;
t1061 = t331 * t785;
t846 = t332 * t787 - t1061;
t167 = -t390 * t766 + t765 * t846;
t801 = t765 * t827 + t1057;
t223 = t567 * t670 - t569 * t671 + t785 * t801;
t800 = t765 * t826 + t1056;
t224 = t568 * t670 - t570 * t671 + t785 * t800;
t44 = (-t292 + t848) * t766 + (t223 * t785 + t224 * t787 + t375) * t765;
t225 = -t672 * t567 - t673 * t569 + t787 * t801;
t226 = -t672 * t568 - t673 * t570 + t787 * t800;
t45 = (-t293 + t847) * t766 + (t225 * t785 + t226 * t787 + t377) * t765;
t62 = (-t306 + t846) * t766 + (t245 * t785 + t246 * t787 + t390) * t765;
t14 = t1163 + (t1007 / 0.2e1 + t1021 / 0.2e1 - t62 / 0.2e1) * t766 + (t45 * t1157 + t44 * t1160 + t167 / 0.2e1) * t765;
t228 = t1112 / 0.2e1;
t92 = t1165 / 0.2e1;
t60 = t228 + t92 - t1107 / 0.2e1;
t1101 = t60 * qJD(4) + t14 * qJD(5);
t356 = t1107 / 0.2e1;
t59 = t356 + t228 - t1165 / 0.2e1;
t1084 = t59 * qJD(5) + (-0.4e1 * t1206 + 0.2e1 * (t1177 + t1179) * (-t704 * t766 + t939)) * qJD(4);
t448 = 0.4e1 * t1206;
t58 = t356 + t92 - t1112 / 0.2e1;
t1083 = t448 * qJD(4) + t58 * qJD(5);
t1066 = t245 * t787;
t1065 = t246 * t785;
t1049 = t684 * t784;
t1009 = t786 * t787;
t682 = Icges(3,5) * t1010 - Icges(3,6) * t1024 - Icges(3,3) * t787;
t956 = -t686 * t1009 - t785 * t682;
t853 = Icges(3,5) * t786 - Icges(3,6) * t784;
t683 = Icges(3,3) * t785 + t787 * t853;
t955 = t687 * t1009 + t785 * t683;
t933 = qJD(2) + qJD(3);
t924 = t1168 / 0.2e1 + t993;
t902 = t1037 / 0.2e1;
t899 = t1036 / 0.2e1;
t888 = -t724 - t1104;
t887 = -t725 - t1104;
t634 = t687 * t1010;
t880 = t787 * t683 - t634;
t878 = t685 * t784 - t682;
t876 = t934 * t1105;
t872 = -t633 + t888;
t104 = m(5) * t242 + m(6) * t160;
t111 = m(5) * t282 + m(6) * t175;
t852 = -Icges(3,5) * t784 - Icges(3,6) * t786;
t851 = Icges(4,5) * t765 + Icges(4,6) * t766;
t834 = t888 + t1202;
t116 = -t223 * t787 + t224 * t785;
t117 = -t225 * t787 + t226 * t785;
t591 = t628 * t785;
t592 = t628 * t787;
t593 = t630 * t785;
t594 = t630 * t787;
t691 = t851 * t785;
t692 = t787 * t851;
t830 = (t117 + (-t702 * t592 - t703 * t594) * t785 - t778 * t692 + (t702 * t591 + t703 * t593 + t785 * t691 + t1218) * t787) * t1160 + (t116 - (t591 * t700 - t593 * t701) * t787 - t779 * t691 + (t700 * t592 - t701 * t594 + t787 * t692 + t1218) * t785) * t1158;
t819 = -t876 + t949;
t814 = t380 * t409 + (-t415 * t787 + t417 * t785) * t658;
t812 = t1190 + t1256;
t810 = t44 * t1158 + t116 * t902 + t45 * t1160 + t117 * t899 + (t1065 - t1066) * t1162 + (t331 * t787 + t332 * t785) * t1243 - t993 + t1230 * t766;
t803 = t1229 * t901 + t1260;
t90 = -t307 * t766 + (t231 * t785 + t232 * t787) * t765;
t91 = -t308 * t766 + (t233 * t785 + t234 * t787) * t765;
t797 = t124 * t902 + t125 * t899 + t167 * t1244 - t44 * t1037 / 0.2e1 - t45 * t1036 / 0.2e1 + t62 * t1242 + t91 * t1160 + t90 * t1158 - t1163 + (t1021 + t1007 + t1063 - t1064) * t1162;
t792 = t803 + t812 - t877;
t791 = t1229 * t903 - t1260 + t812 + t877;
t790 = -t1190 + t1256 + t803 + t877;
t789 = t1065 / 0.2e1 - t1066 / 0.2e1 + (t702 * t629 + t703 * t631 + t785 * t717 + t293 + t1227 * t787 - t1252 * t766 + (t592 * t781 - t594 * t782 + t1254) * t765) * t1160 + (-t629 * t700 + t631 * t701 - t1048 + t292 + t1227 * t785 + t1251 * t766 + (t591 * t781 - t593 * t782 - t1255) * t765) * t1158 - t795;
t788 = -t1061 / 0.2e1 - t793 + t1226 * t1161 + (t331 + t1226) * t1160;
t739 = -rSges(3,2) * t784 + t1100;
t706 = t852 * t787;
t705 = t852 * t785;
t642 = t887 * t787;
t640 = t887 * t785;
t516 = t872 * t787;
t514 = t872 * t785;
t486 = -t876 + t552;
t458 = -t685 * t1023 + t955;
t457 = -t684 * t1023 - t956;
t456 = -t685 * t1024 - t880;
t443 = -t658 * t1036 - t766 * t541;
t442 = t658 * t1037 + t540 * t766;
t429 = t834 * t787;
t427 = t834 * t785;
t389 = -t1106 / 0.2e1;
t378 = t819 + t1204;
t365 = -t457 * t787 + t458 * t785;
t364 = -(-t785 * (-t686 * t786 + t1049) - t787 * t682) * t787 + t456 * t785;
t355 = -t1032 + (-t961 * t763 + t962 * t764) * t765;
t314 = t1108 / 0.2e1;
t305 = t819 + t1194;
t205 = t1109 + t1130;
t193 = (t456 - t634 + (t683 + t1049) * t787 + t956) * t787 + t955 * t785;
t192 = (t787 * t878 + t458 - t955) * t787 + (t785 * t878 + t457 + t880) * t785;
t166 = t1118 / 0.2e1;
t158 = t1120 / 0.2e1;
t147 = t222 + (t438 * t785 + t440 * t787) * t658;
t131 = t202 + (t426 * t785 + t428 * t787) * t658;
t105 = t1127 / 0.2e1;
t101 = t314 + t1106 / 0.2e1;
t100 = t389 + t314;
t99 = t389 - t1108 / 0.2e1;
t97 = t1129 / 0.2e1;
t76 = t1110 + t1132 + t793 + t1148;
t63 = t1222 + (t735 / 0.2e1 + t734 / 0.2e1) * t786 + t793 + t1156 + t1149 + t1133 + t1111;
t56 = t59 * qJD(4);
t36 = t1171 / 0.2e1;
t34 = t1172 / 0.2e1;
t31 = m(6) * t147 + t993;
t29 = t983 + t985;
t27 = m(6) * t131 + t993;
t25 = t984 + t986;
t23 = t994 - t1219;
t22 = t1085 + t1086 - t994;
t21 = t994 + t1219;
t20 = t830 + t1262;
t19 = t20 * qJD(3);
t18 = t830 + t1261;
t16 = t36 - t1172 / 0.2e1 + t924;
t15 = t34 - t1171 / 0.2e1 + t924;
t11 = (t364 / 0.2e1 + t192 / 0.2e1) * t785 + (-t193 / 0.2e1 + t365 / 0.2e1) * t787 + t795;
t10 = t34 + t36 - t1168 / 0.2e1 + t810;
t9 = t795 + t1225;
t8 = t795 - t1225;
t7 = t105 + t792 - t1118 / 0.2e1;
t6 = t791 + t166 - t1127 / 0.2e1;
t5 = t105 + t166 + t790;
t4 = t97 + t792 - t1120 / 0.2e1;
t3 = t158 + t791 - t1129 / 0.2e1;
t2 = t158 + t97 + t790;
t1 = t789 + t912 + t913;
t12 = [t63 * qJD(2) + t76 * qJD(3) + t205 * qJD(4) + t173 * qJD(5), t63 * qJD(1) + t1 * qJD(3) + t25 * qJD(4) + t2 * qJD(5) + ((t563 * t642 + t564 * t640) * t1180 + (t446 * t516 + t447 * t514 - t484 * t515 - t485 * t513) * t1179 + (t407 * t429 + t408 * t427 - t426 * t460 - t428 * t459) * t1177 + m(3) * ((-t614 * t787 - t615 * t785) * t739 + (-t711 * t787 + t712 * t785) * t737) / 0.2e1) * t1184 + ((t784 * t945 + t786 * t943) * t1160 + t789 + t193 * t1157 + (t364 + t192) * t1161 + (-t784 * t946 + t786 * t944 + t365) * t1158 + (t778 / 0.2e1 + t779 / 0.2e1) * t853) * qJD(2), t76 * qJD(1) + t1 * qJD(2) + t789 * qJD(3) + t29 * qJD(4) + t5 * qJD(5) + ((-t521 * t547 - t522 * t545 + t979) * t1179 + (-t438 * t475 - t440 * t474 + t980) * t1177 + (t820 + (-t697 * t787 + t699 * t785) * t723) * t1180) * t1182, qJD(1) * t205 + qJD(2) * t25 + qJD(3) * t29 + qJD(5) * t100, t1257 + t2 * qJD(2) + t5 * qJD(3) + t100 * qJD(4) + (m(6) * (t407 * t442 + t408 * t443 - t415 * t540 - t417 * t541) - t355 * t766 + ((t308 / 0.2e1 + t263 / 0.2e1) * t787 + (t307 / 0.2e1 + t262 / 0.2e1) * t785) * t765) * qJD(5); (t788 - t1222 - (t735 + t734) * t786 / 0.2e1) * qJD(1) + t11 * qJD(2) + t8 * qJD(3) + t26 * qJD(4) + t3 * qJD(5) + (-t1156 / 0.4e1 - t1149 / 0.4e1 - t1133 / 0.4e1 - t1111 / 0.4e1) * t1185, t11 * qJD(1) + (m(6) * (t251 * t305 - t426 * t427 - t428 * t429) + m(5) * (t324 * t378 - t513 * t514 - t515 * t516) + m(4) * (-t1207 * t642 - t1208 * t640 + t401 * t486) + m(3) * ((t785 * (rSges(3,1) * t1010 - t935) + t787 * (rSges(3,1) * t1009 + t785 * rSges(3,3) - t756)) * (-t785 * t711 - t712 * t787) + t934 * t739 * t737) + (t779 * t705 + (-t787 * t706 + t1189) * t785) * t1158 + (t778 * t706 + (-t785 * t705 + t1189) * t787) * t1160 + t830) * qJD(2) + t18 * qJD(3) + t104 * qJD(4) + t27 * qJD(5), t8 * qJD(1) + t18 * qJD(2) + t22 * qJD(4) + t15 * qJD(5) + ((t110 + t103) * t1177 + (t196 + t169) * t1179 + (t299 + t363) * t1180) * t1182 + (t830 - t1262) * qJD(3), qJD(2) * t104 + qJD(3) * t22 + t1084 + t1259, t3 * qJD(1) + t27 * qJD(2) + t15 * qJD(3) + t56 + (t797 + m(6) * (t392 * t251 - t443 * t426 - t442 * t428 + t814)) * qJD(5); t788 * qJD(1) + t9 * qJD(2) + t795 * qJD(3) + t30 * qJD(4) + t6 * qJD(5) + (-t1148 / 0.4e1 - t1132 / 0.4e1 - t1110 / 0.4e1) * t1185, t9 * qJD(1) + t19 + t23 * qJD(4) + t16 * qJD(5) + ((t371 * t378 - t514 * t545 - t516 * t547 + t169) * t1179 + (t511 * t486 + (-t640 * t785 - t642 * t787) * t723 + t299) * t1180 + (t285 * t305 - t427 * t438 - t429 * t440 + t103) * t1177) * t1184 + (t830 - t1261) * qJD(2), qJD(1) * t795 + qJD(2) * t20 + qJD(4) * t111 + qJD(5) * t31 + t19, qJD(2) * t23 + qJD(3) * t111 + t1084 + t1258, t6 * qJD(1) + t16 * qJD(2) + t31 * qJD(3) + t56 + (t797 + m(6) * (t392 * t285 - t443 * t438 - t442 * t440 + t814)) * qJD(5); -t26 * qJD(2) - t30 * qJD(3) + t99 * qJD(5) + (-t1109 / 0.4e1 - t1130 / 0.4e1) * t1185 + 0.2e1 * (t398 * t1177 + t436 * t1179 + (-t408 * t1177 - t447 * t1179) * t1036) * qJD(1), -t1259 + (m(6) * (-t766 * t305 + (t427 * t785 + t429 * t787) * t765 + t909) + m(5) * (-t766 * t378 + (t514 * t785 + t516 * t787) * t765 + t908) - t104) * qJD(2) + t21 * qJD(3) + t1083, -t1258 + t21 * qJD(2) + (-t111 + t1125 + t1143) * qJD(3) + t1083, t933 * t448, t99 * qJD(1) + m(6) * (-t392 * t766 + (t442 * t787 + t443 * t785) * t765) * qJD(5) + t933 * t58; t4 * qJD(2) + t7 * qJD(3) + t101 * qJD(4) - t1257, t4 * qJD(1) + ((t305 * t380 + t415 * t429 - t417 * t427 - t131 + t911) * m(6) + t810) * qJD(2) + t10 * qJD(3) + t1101, t7 * qJD(1) + t10 * qJD(2) + ((t37 - t147) * m(6) + t810) * qJD(3) + t1101, qJD(1) * t101 + t60 * t933, (m(6) * (t380 * t392 + t415 * t442 - t417 * t443) + t766 ^ 2 * t355 / 0.2e1 + (t91 * t1157 + t90 * t1160 + (t262 * t785 + t263 * t787) * t1162) * t765) * qJD(5) + t933 * t14;];
Cq = t12;