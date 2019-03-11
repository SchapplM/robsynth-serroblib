% Calculate vector of inverse dynamics joint torques for
% S6RPRPRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPRPRR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR2_invdynJ_fixb_slag_vp1: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_invdynJ_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_invdynJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR2_invdynJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR2_invdynJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:36:59
% EndTime: 2019-03-09 03:38:56
% DurationCPUTime: 111.65s
% Computational Cost: add. (59201->1484), mult. (48706->1928), div. (0->0), fcn. (45117->12), ass. (0->691)
t1135 = Icges(4,3) + Icges(5,3);
t599 = qJ(3) + pkin(11);
t589 = sin(t599);
t591 = cos(t599);
t468 = Icges(5,5) * t591 - Icges(5,6) * t589;
t604 = sin(qJ(3));
t607 = cos(qJ(3));
t530 = Icges(4,5) * t607 - Icges(4,6) * t604;
t1128 = t468 + t530;
t600 = qJ(1) + pkin(10);
t592 = cos(t600);
t1134 = t1135 * t592;
t590 = sin(t600);
t926 = t590 * t607;
t928 = t590 * t604;
t932 = t590 * t591;
t936 = t589 * t590;
t1123 = -Icges(4,5) * t926 - Icges(5,5) * t932 + Icges(4,6) * t928 + Icges(5,6) * t936 + t1134;
t1125 = t1128 * t592 + t1135 * t590;
t961 = Icges(5,6) * t592;
t339 = Icges(5,4) * t932 - Icges(5,2) * t936 - t961;
t962 = Icges(4,6) * t592;
t380 = Icges(4,4) * t926 - Icges(4,2) * t928 - t962;
t1133 = t339 * t589 + t380 * t604;
t975 = Icges(5,4) * t589;
t472 = Icges(5,1) * t591 - t975;
t342 = Icges(5,5) * t590 + t472 * t592;
t976 = Icges(4,4) * t604;
t534 = Icges(4,1) * t607 - t976;
t385 = Icges(4,5) * t590 + t534 * t592;
t1132 = -t342 * t932 - t385 * t926;
t506 = Icges(5,4) * t936;
t967 = Icges(5,5) * t592;
t341 = Icges(5,1) * t932 - t506 - t967;
t542 = Icges(4,4) * t928;
t968 = Icges(4,5) * t592;
t384 = Icges(4,1) * t926 - t542 - t968;
t1118 = -t341 * t591 - t384 * t607 + t1133;
t469 = Icges(5,2) * t591 + t975;
t570 = Icges(5,4) * t591;
t471 = Icges(5,1) * t589 + t570;
t531 = Icges(4,2) * t607 + t976;
t595 = Icges(4,4) * t607;
t533 = Icges(4,1) * t604 + t595;
t1127 = t469 * t589 - t471 * t591 + t531 * t604 - t533 * t607;
t1131 = t1125 * t592 + t1132;
t914 = t592 * t607;
t925 = t591 * t592;
t1130 = t1123 * t590 - t341 * t925 - t384 * t914;
t1098 = t1125 * t590 + t342 * t925 + t385 * t914;
t1085 = -t1118 * t590 + t1123 * t592;
t723 = -Icges(5,2) * t589 + t570;
t340 = Icges(5,6) * t590 + t592 * t723;
t724 = -Icges(4,2) * t604 + t595;
t381 = Icges(4,6) * t590 + t592 * t724;
t1084 = -t340 * t936 - t381 * t928 - t1131;
t916 = t592 * t604;
t935 = t589 * t592;
t1083 = -t339 * t935 - t380 * t916 - t1130;
t1082 = -t340 * t935 - t381 * t916 + t1098;
t467 = Icges(5,5) * t589 + Icges(5,6) * t591;
t529 = Icges(4,5) * t604 + Icges(4,6) * t607;
t1129 = t467 + t529;
t1126 = -t469 * t591 - t471 * t589 - t531 * t607 - t533 * t604;
t938 = t529 * t592;
t940 = t467 * t592;
t1103 = -t1127 * t590 - t938 - t940;
t939 = t529 * t590;
t941 = t467 * t590;
t1102 = -t1127 * t592 + t939 + t941;
t1124 = t340 * t589 + t381 * t604;
t605 = sin(qJ(1));
t1013 = t605 * pkin(1);
t602 = -qJ(4) - pkin(7);
t554 = t592 * t602;
t596 = t607 * pkin(3);
t585 = t596 + pkin(2);
t874 = -t590 * t585 - t554;
t763 = t874 - t1013;
t1081 = t339 * t591 + t341 * t589 + t380 * t607 + t384 * t604;
t1080 = t340 * t591 + t342 * t589 + t381 * t607 + t385 * t604;
t452 = t723 * qJD(3);
t453 = t472 * qJD(3);
t494 = t724 * qJD(3);
t495 = t534 * qJD(3);
t1122 = qJD(1) * t1129 + qJD(3) * t1126 - t452 * t589 + t453 * t591 - t494 * t604 + t495 * t607;
t1121 = t1129 * qJD(3);
t839 = rSges(5,1) * t932;
t1120 = -t839 + t763;
t1119 = t342 * t591 + t385 * t607 - t1124;
t1117 = t1082 * t590 - t1083 * t592;
t1116 = t1084 * t590 - t1085 * t592;
t1115 = t1127 * qJD(1) + qJD(3) * t1128;
t608 = cos(qJ(1));
t597 = t608 * pkin(1);
t1114 = t1102 * qJD(1);
t1113 = t1103 * qJD(1);
t1112 = t1125 * qJD(1);
t477 = rSges(3,1) * t590 + rSges(3,2) * t592;
t448 = -t477 - t1013;
t606 = cos(qJ(5));
t915 = t592 * t606;
t603 = sin(qJ(5));
t922 = t591 * t603;
t426 = t590 * t922 + t915;
t917 = t592 * t603;
t921 = t591 * t606;
t427 = t590 * t921 - t917;
t754 = t427 * rSges(6,1) - t426 * rSges(6,2);
t236 = -rSges(6,3) * t936 - t754;
t1006 = rSges(6,1) * t606;
t753 = -rSges(6,2) * t603 + t1006;
t386 = -rSges(6,3) * t591 + t589 * t753;
t857 = qJD(5) * t589;
t859 = qJD(3) * t592;
t443 = -t590 * t857 + t859;
t856 = qJD(5) * t591;
t546 = qJD(1) - t856;
t561 = qJD(4) * t590;
t1018 = pkin(3) * t604;
t480 = pkin(4) * t589 - pkin(8) * t591;
t794 = -t480 - t1018;
t658 = t794 * t859 + t561;
t1111 = t236 * t546 - t386 * t443 + t658;
t406 = Icges(6,4) * t427;
t227 = -Icges(6,2) * t426 + Icges(6,6) * t936 + t406;
t405 = Icges(6,4) * t426;
t231 = -Icges(6,1) * t427 - Icges(6,5) * t936 + t405;
t1092 = t227 * t603 + t231 * t606;
t224 = Icges(6,5) * t427 - Icges(6,6) * t426 + Icges(6,3) * t936;
t94 = -t1092 * t589 - t224 * t591;
t601 = qJ(5) + qJ(6);
t593 = sin(t601);
t919 = t592 * t593;
t594 = cos(t601);
t923 = t591 * t594;
t391 = t590 * t923 - t919;
t364 = Icges(7,4) * t391;
t918 = t592 * t594;
t924 = t591 * t593;
t390 = t590 * t924 + t918;
t188 = -Icges(7,2) * t390 + Icges(7,6) * t936 + t364;
t363 = Icges(7,4) * t390;
t192 = -Icges(7,1) * t391 - Icges(7,5) * t936 + t363;
t1093 = t188 * t593 + t192 * t594;
t185 = Icges(7,5) * t391 - Icges(7,6) * t390 + Icges(7,3) * t936;
t87 = -t1093 * t589 - t185 * t591;
t1110 = t592 ^ 2;
t611 = qJD(1) ^ 2;
t840 = t611 * t597;
t1109 = qJD(3) * t1116 + t1113;
t1108 = qJD(3) * t1117 + t1114;
t667 = qJD(3) * t469;
t219 = qJD(1) * t340 - t590 * t667;
t669 = qJD(3) * t471;
t221 = qJD(1) * t342 - t590 * t669;
t668 = qJD(3) * t531;
t267 = qJD(1) * t381 - t590 * t668;
t670 = qJD(3) * t533;
t270 = qJD(1) * t385 - t590 * t670;
t1107 = qJD(3) * t1118 - t219 * t591 - t221 * t589 - t267 * t607 - t270 * t604;
t218 = -t592 * t667 + (-t590 * t723 + t961) * qJD(1);
t220 = -t592 * t669 + (-t472 * t590 + t967) * qJD(1);
t266 = -t592 * t668 + (-t590 * t724 + t962) * qJD(1);
t269 = -t592 * t670 + (-t534 * t590 + t968) * qJD(1);
t1106 = qJD(3) * t1119 + t218 * t591 + t220 * t589 + t266 * t607 + t269 * t604;
t1105 = t1115 * t590 + t1122 * t592;
t1104 = -t1115 * t592 + t1122 * t590;
t1101 = -qJD(3) * t1080 - t218 * t589 + t220 * t591 - t266 * t604 + t269 * t607 + t1112;
t1100 = t1123 * qJD(1) + qJD(3) * t1081 + t219 * t589 - t221 * t591 + t267 * t604 - t270 * t607;
t1099 = t1123 + t1124;
t563 = qJD(3) * t590;
t1097 = qJD(1) * t1118 - t1121 * t590 + t1112;
t1096 = -t1121 * t592 + (-t1128 * t590 - t1119 + t1134) * qJD(1);
t598 = qJD(5) + qJD(6);
t780 = t591 * t598;
t1059 = qJD(1) - t780;
t1069 = -t391 * rSges(7,1) + t390 * rSges(7,2);
t196 = -rSges(7,3) * t936 + t1069;
t609 = -pkin(9) - pkin(8);
t1008 = pkin(8) + t609;
t1063 = t1008 * t589;
t581 = t591 * pkin(4);
t584 = pkin(5) * t606 + pkin(4);
t520 = t591 * t584;
t841 = pkin(5) * t917;
t767 = -t590 * t520 + t841;
t214 = (t581 + t1063) * t590 + t767;
t1009 = pkin(4) - t584;
t1061 = t1009 * t589;
t1062 = t1008 * t591;
t331 = -t1061 + t1062;
t855 = qJD(6) * t589;
t336 = -t590 * t855 + t443;
t1002 = rSges(7,2) * t593;
t1005 = rSges(7,1) * t594;
t750 = -t1002 + t1005;
t359 = -rSges(7,3) * t591 + t589 * t750;
t1095 = t1059 * t196 + t546 * t214 - t443 * t331 - t336 * t359 + t658;
t719 = Icges(7,5) * t594 - Icges(7,6) * t593;
t353 = -Icges(7,3) * t591 + t589 * t719;
t969 = Icges(7,4) * t594;
t721 = -Icges(7,2) * t593 + t969;
t355 = -Icges(7,6) * t591 + t589 * t721;
t970 = Icges(7,4) * t593;
t725 = Icges(7,1) * t594 - t970;
t357 = -Icges(7,5) * t591 + t589 * t725;
t120 = t353 * t936 - t355 * t390 + t357 * t391;
t442 = t592 * t857 + t563;
t335 = t592 * t855 + t442;
t75 = t185 * t936 - t188 * t390 - t192 * t391;
t930 = t590 * t594;
t392 = -t591 * t919 + t930;
t931 = t590 * t593;
t393 = t591 * t918 + t931;
t187 = Icges(7,5) * t393 + Icges(7,6) * t392 + Icges(7,3) * t935;
t971 = Icges(7,4) * t393;
t190 = Icges(7,2) * t392 + Icges(7,6) * t935 + t971;
t365 = Icges(7,4) * t392;
t193 = Icges(7,1) * t393 + Icges(7,5) * t935 + t365;
t76 = t187 * t936 - t390 * t190 + t391 * t193;
t36 = t1059 * t120 + t335 * t76 - t336 * t75;
t720 = Icges(6,5) * t606 - Icges(6,6) * t603;
t374 = -Icges(6,3) * t591 + t589 * t720;
t972 = Icges(6,4) * t606;
t722 = -Icges(6,2) * t603 + t972;
t378 = -Icges(6,6) * t591 + t589 * t722;
t973 = Icges(6,4) * t603;
t726 = Icges(6,1) * t606 - t973;
t382 = -Icges(6,5) * t591 + t589 * t726;
t128 = t374 * t936 - t378 * t426 + t382 * t427;
t81 = t224 * t936 - t227 * t426 - t231 * t427;
t927 = t590 * t606;
t428 = -t591 * t917 + t927;
t929 = t590 * t603;
t429 = t591 * t915 + t929;
t226 = Icges(6,5) * t429 + Icges(6,6) * t428 + Icges(6,3) * t935;
t974 = Icges(6,4) * t429;
t229 = Icges(6,2) * t428 + Icges(6,6) * t935 + t974;
t407 = Icges(6,4) * t428;
t232 = Icges(6,1) * t429 + Icges(6,5) * t935 + t407;
t82 = t226 * t936 - t426 * t229 + t427 * t232;
t41 = t128 * t546 + t442 * t82 - t443 * t81;
t129 = t374 * t935 + t378 * t428 + t382 * t429;
t83 = t224 * t935 + t428 * t227 - t231 * t429;
t84 = t226 * t935 + t428 * t229 + t429 * t232;
t42 = t129 * t546 + t442 * t84 - t443 * t83;
t121 = t353 * t935 + t355 * t392 + t357 * t393;
t77 = t185 * t935 + t392 * t188 - t192 * t393;
t78 = t187 * t935 + t392 * t190 + t393 * t193;
t37 = t1059 * t121 + t335 * t78 - t336 * t77;
t1079 = t1097 * t1110 + (t1101 * t590 + (-t1096 + t1100) * t592) * t590;
t1078 = t1100 * t1110 + (t1096 * t590 + (-t1097 + t1101) * t592) * t590;
t799 = t859 / 0.2e1;
t864 = qJD(1) * t590;
t822 = t589 * t864;
t1077 = t591 * t799 - t822 / 0.2e1;
t801 = t563 / 0.2e1;
t862 = qJD(1) * t592;
t804 = t862 / 0.2e1;
t1076 = t589 * t804 + t591 * t801;
t849 = qJD(3) * qJD(5);
t1075 = qJDD(5) * t589 + t591 * t849;
t860 = qJD(3) * t591;
t661 = t589 * t862 + t590 * t860;
t197 = t393 * rSges(7,1) + t392 * rSges(7,2) + rSges(7,3) * t935;
t518 = pkin(4) * t925;
t419 = pkin(8) * t935 + t518;
t547 = pkin(5) * t929;
t933 = t589 * t609;
t690 = t584 * t925 - t592 * t933 + t547;
t215 = t690 - t419;
t1066 = t589 * pkin(8) + t581;
t417 = t1066 * t590;
t582 = t592 * pkin(7);
t481 = pkin(2) * t590 - t582;
t329 = t481 + t874;
t580 = t590 * pkin(7);
t483 = t592 * pkin(2) + t580;
t784 = t592 * t585 - t590 * t602;
t330 = t784 - t483;
t810 = -t329 * t563 + t330 * t859 + qJD(2);
t729 = t417 * t563 + t419 * t859 + t810;
t48 = -t196 * t335 + t197 * t336 - t214 * t442 + t215 * t443 + t729;
t908 = -t196 - t214;
t1073 = t48 * t908;
t465 = qJD(1) * t481;
t1070 = qJD(1) * t329 - t465;
t574 = t590 * rSges(4,3);
t389 = rSges(4,1) * t914 - rSges(4,2) * t916 + t574;
t796 = t483 + t597;
t297 = t389 + t796;
t1058 = t591 * rSges(5,1) - rSges(5,2) * t589;
t1068 = t1058 + t596;
t479 = t592 * rSges(3,1) - rSges(3,2) * t590;
t449 = t479 + t597;
t1067 = -rSges(5,2) * t936 - t592 * rSges(5,3);
t1014 = g(2) * t590;
t1057 = g(1) * t592 + t1014;
t1050 = t1070 + (-t1013 - t417) * qJD(1);
t861 = qJD(3) * t589;
t648 = t1059 * t594 + t593 * t861;
t863 = qJD(1) * t591;
t781 = -t598 + t863;
t169 = t592 * t648 + t781 * t931;
t647 = t1059 * t593 - t594 * t861;
t170 = t592 * t647 - t781 * t930;
t817 = t591 * t859;
t827 = t170 * rSges(7,1) + t169 * rSges(7,2) + rSges(7,3) * t817;
t104 = -rSges(7,3) * t822 + t827;
t171 = t590 * t648 - t781 * t919;
t172 = t590 * t647 + t781 * t918;
t752 = rSges(7,1) * t172 + rSges(7,2) * t171;
t105 = rSges(7,3) * t661 + t752;
t492 = pkin(8) * t817;
t1001 = pkin(5) * qJD(5);
t830 = t606 * t1001;
t824 = qJD(1) * t841 + t590 * t830 + t609 * t822;
t831 = t603 * t1001;
t124 = -t492 + (pkin(8) * t864 + t1009 * t859) * t589 + ((-qJD(3) * t609 - t831) * t592 + t1009 * t864) * t591 + t824;
t843 = pkin(4) * t936;
t491 = qJD(3) * t843;
t657 = -t1009 * t591 - t1063;
t778 = t592 * t830;
t779 = t591 * t831;
t937 = t584 * t589;
t125 = -t778 + t491 + (-t779 + (-t937 - t1062) * qJD(3)) * t590 + (t592 * t657 + t547) * qJD(1);
t851 = qJD(1) * qJD(3);
t457 = qJDD(3) * t590 + t592 * t851;
t774 = t1075 * t592 + t457;
t153 = qJD(6) * t817 + (qJDD(6) * t592 - t598 * t864) * t589 + t774;
t458 = -qJDD(3) * t592 + t590 * t851;
t807 = qJD(1) * t857;
t255 = t1075 * t590 + t592 * t807 + t458;
t154 = qJD(6) * t661 + qJDD(6) * t936 + t255;
t254 = -t590 * t807 + t774;
t819 = t589 * t859;
t662 = -t590 * t863 - t819;
t273 = pkin(4) * t662 - pkin(8) * t822 + t492;
t274 = pkin(8) * t661 + qJD(1) * t518 - t491;
t1010 = pkin(2) - t585;
t558 = pkin(7) * t862;
t813 = t604 * t859;
t692 = -pkin(3) * t813 + t561;
t243 = -t558 + (t1010 * t590 - t554) * qJD(1) + t692;
t846 = pkin(3) * t928;
t872 = qJD(3) * t846 + qJD(4) * t592;
t823 = t602 * t864 + t872;
t244 = (-t1010 * t592 - t580) * qJD(1) - t823;
t769 = t243 * t859 + t244 * t563 - t457 * t329 + qJDD(2);
t895 = -t330 - t419;
t633 = t273 * t859 + t274 * t563 + t457 * t417 + t458 * t895 + t769;
t11 = t104 * t336 + t105 * t335 + t124 * t443 + t125 * t442 - t153 * t196 - t154 * t197 - t214 * t254 - t215 * t255 + t633;
t907 = t197 + t215;
t911 = t104 + t124;
t1049 = -t11 * t907 - t48 * t911;
t1048 = t590 * (-t469 * t592 + t342) - t592 * (-Icges(5,2) * t932 + t341 - t506);
t885 = -Icges(4,2) * t926 + t384 - t542;
t888 = t533 * t590 + t380;
t1047 = -t604 * t885 - t607 * t888;
t354 = Icges(7,3) * t589 + t591 * t719;
t705 = -t355 * t593 + t357 * t594;
t713 = -t190 * t593 + t193 * t594;
t1046 = t335 * (-t353 * t592 - t713) - t336 * (-t353 * t590 + t1093) + t1059 * (t354 - t705);
t421 = (-Icges(7,2) * t594 - t970) * t589;
t1045 = t335 * (-Icges(7,2) * t393 + t193 + t365) - t336 * (-Icges(7,2) * t391 - t192 - t363) + t1059 * (t357 + t421);
t375 = Icges(6,3) * t589 + t591 * t720;
t704 = -t378 * t603 + t382 * t606;
t710 = -t229 * t603 + t232 * t606;
t1044 = t442 * (-t374 * t592 - t710) - t443 * (-t374 * t590 + t1092) + t546 * (t375 - t704);
t436 = (-Icges(6,2) * t606 - t973) * t589;
t1043 = t442 * (-Icges(6,2) * t429 + t232 + t407) - t443 * (-Icges(6,2) * t427 - t231 - t405) + t546 * (t382 + t436);
t1042 = t153 / 0.2e1;
t1041 = t154 / 0.2e1;
t1040 = t254 / 0.2e1;
t1039 = t255 / 0.2e1;
t847 = t589 * t849 + qJDD(1);
t328 = qJD(3) * t855 + (-qJDD(5) - qJDD(6)) * t591 + t847;
t1038 = t328 / 0.2e1;
t1037 = -t335 / 0.2e1;
t1036 = t335 / 0.2e1;
t1035 = -t336 / 0.2e1;
t1034 = t336 / 0.2e1;
t1033 = -t442 / 0.2e1;
t1032 = t442 / 0.2e1;
t1031 = -t443 / 0.2e1;
t1030 = t443 / 0.2e1;
t450 = -qJDD(5) * t591 + t847;
t1029 = t450 / 0.2e1;
t1028 = t457 / 0.2e1;
t1027 = t458 / 0.2e1;
t1026 = -t1059 / 0.2e1;
t1025 = t1059 / 0.2e1;
t1024 = -t546 / 0.2e1;
t1023 = t546 / 0.2e1;
t1022 = t590 / 0.2e1;
t1021 = -t591 / 0.2e1;
t1020 = -t592 / 0.2e1;
t1019 = -rSges(6,3) - pkin(8);
t1017 = pkin(5) * t603;
t1016 = g(1) * t590;
t1007 = rSges(4,1) * t607;
t101 = Icges(7,4) * t172 + Icges(7,2) * t171 + Icges(7,6) * t661;
t103 = Icges(7,1) * t172 + Icges(7,4) * t171 + Icges(7,5) * t661;
t99 = Icges(7,5) * t172 + Icges(7,6) * t171 + Icges(7,3) * t661;
t22 = (-qJD(3) * t1093 - t99) * t591 + (qJD(3) * t185 + (-t188 * t598 + t103) * t594 + (t192 * t598 - t101) * t593) * t589;
t1000 = t22 * t336;
t660 = t817 - t822;
t100 = Icges(7,4) * t170 + Icges(7,2) * t169 + Icges(7,6) * t660;
t102 = Icges(7,1) * t170 + Icges(7,4) * t169 + Icges(7,5) * t660;
t98 = Icges(7,5) * t170 + Icges(7,6) * t169 + Icges(7,3) * t660;
t23 = (qJD(3) * t713 - t98) * t591 + (qJD(3) * t187 + (-t190 * t598 + t102) * t594 + (-t193 * t598 - t100) * t593) * t589;
t999 = t23 * t335;
t646 = t546 * t606 + t603 * t861;
t777 = -qJD(5) + t863;
t207 = t590 * t646 - t777 * t917;
t645 = t546 * t603 - t606 * t861;
t208 = t590 * t645 + t777 * t915;
t113 = Icges(6,5) * t208 + Icges(6,6) * t207 + Icges(6,3) * t661;
t115 = Icges(6,4) * t208 + Icges(6,2) * t207 + Icges(6,6) * t661;
t117 = Icges(6,1) * t208 + Icges(6,4) * t207 + Icges(6,5) * t661;
t28 = (-qJD(3) * t1092 - t113) * t591 + (qJD(3) * t224 - t115 * t603 + t117 * t606 + (-t227 * t606 + t231 * t603) * qJD(5)) * t589;
t998 = t28 * t443;
t205 = t592 * t646 + t777 * t929;
t206 = t592 * t645 - t777 * t927;
t112 = Icges(6,5) * t206 + Icges(6,6) * t205 + Icges(6,3) * t660;
t114 = Icges(6,4) * t206 + Icges(6,2) * t205 + Icges(6,6) * t660;
t116 = Icges(6,1) * t206 + Icges(6,4) * t205 + Icges(6,5) * t660;
t29 = (qJD(3) * t710 - t112) * t591 + (qJD(3) * t226 - t114 * t603 + t116 * t606 + (-t229 * t606 - t232 * t603) * qJD(5)) * t589;
t997 = t29 * t442;
t755 = rSges(6,1) * t208 + rSges(6,2) * t207;
t119 = rSges(6,3) * t661 + t755;
t444 = (-rSges(6,1) * t603 - rSges(6,2) * t606) * t589;
t572 = t589 * rSges(6,3);
t275 = qJD(5) * t444 + (t591 * t753 + t572) * qJD(3);
t455 = t1066 * qJD(3);
t850 = qJD(1) * qJD(4);
t675 = qJDD(4) * t590 + t458 * t1018 + t592 * t850 - t840;
t797 = -t481 - t1013;
t766 = t329 + t797;
t695 = -t417 + t766;
t913 = t607 * qJD(3) ^ 2;
t844 = pkin(3) * t913;
t456 = t483 * qJD(1);
t903 = -t244 - t456;
t614 = t458 * t480 + (-t274 + t903) * qJD(1) + (-qJD(3) * t455 - t844) * t592 + t695 * qJDD(1) + t675;
t38 = -t119 * t546 + t236 * t450 + t255 * t386 - t275 * t443 + t614;
t994 = t38 * t590;
t826 = t206 * rSges(6,1) + t205 * rSges(6,2) + rSges(6,3) * t817;
t118 = -rSges(6,3) * t822 + t826;
t237 = t429 * rSges(6,1) + t428 * rSges(6,2) + rSges(6,3) * t935;
t768 = qJDD(1) * t597 - t1013 * t611;
t674 = qJD(1) * (-pkin(2) * t864 + t558) + qJDD(1) * t483 + t768;
t617 = qJD(1) * t243 + qJDD(1) * t330 + t590 * t850 + (-t457 * t604 - t590 * t913) * pkin(3) - qJDD(4) * t592 + t674;
t613 = qJD(1) * t273 + qJDD(1) * t419 - t455 * t563 - t457 * t480 + t617;
t39 = t118 * t546 + t237 * t450 - t254 * t386 - t275 * t442 + t613;
t993 = t39 * t592;
t571 = t589 * rSges(7,3);
t573 = t590 * rSges(5,3);
t664 = t695 * qJD(1);
t60 = t1095 + t664;
t990 = t590 * t60;
t765 = t330 + t796;
t634 = (t419 + t765) * qJD(1) - t480 * t563 - t872;
t86 = t237 * t546 - t386 * t442 + t634;
t989 = t590 * t86;
t61 = t1059 * t197 + t215 * t546 - t331 * t442 - t335 * t359 + t634;
t988 = t592 * t61;
t987 = t61 * t215;
t986 = t87 * t154;
t88 = -t187 * t591 + t589 * t713;
t985 = t88 * t153;
t984 = t94 * t255;
t95 = -t226 * t591 + t589 * t710;
t983 = t95 * t254;
t981 = -rSges(7,3) + t609;
t140 = -t353 * t591 + t589 * t705;
t420 = (-Icges(7,5) * t593 - Icges(7,6) * t594) * t589;
t202 = qJD(3) * t354 + t420 * t598;
t356 = Icges(7,6) * t589 + t591 * t721;
t203 = qJD(3) * t356 + t421 * t598;
t358 = Icges(7,5) * t589 + t591 * t725;
t422 = (-Icges(7,1) * t593 - t969) * t589;
t204 = qJD(3) * t358 + t422 * t598;
t53 = (qJD(3) * t705 - t202) * t591 + (qJD(3) * t353 + (-t355 * t598 + t204) * t594 + (-t357 * t598 - t203) * t593) * t589;
t980 = t1059 * t53 + t140 * t328;
t142 = -t374 * t591 + t589 * t704;
t433 = (-Icges(6,5) * t603 - Icges(6,6) * t606) * t589;
t262 = qJD(3) * t375 + qJD(5) * t433;
t379 = Icges(6,6) * t589 + t591 * t722;
t265 = qJD(3) * t379 + qJD(5) * t436;
t383 = Icges(6,5) * t589 + t591 * t726;
t439 = (-Icges(6,1) * t603 - t972) * t589;
t268 = qJD(3) * t383 + qJD(5) * t439;
t63 = (qJD(3) * t704 - t262) * t591 + (qJD(3) * t374 - t265 * t603 + t268 * t606 + (-t378 * t606 - t382 * t603) * qJD(5)) * t589;
t979 = t142 * t450 + t63 * t546;
t476 = rSges(5,1) * t589 + rSges(5,2) * t591;
t663 = -t476 - t1018;
t659 = t663 * t859 + t561;
t344 = t839 + t1067;
t696 = -t344 + t766;
t133 = qJD(1) * t696 + t659;
t952 = t133 * t476;
t869 = rSges(4,2) * t928 + t592 * rSges(4,3);
t387 = rSges(4,1) * t926 - t869;
t764 = -t387 + t797;
t536 = rSges(4,1) * t604 + rSges(4,2) * t607;
t814 = t536 * t859;
t198 = qJD(1) * t764 - t814;
t950 = t198 * t590;
t949 = t198 * t592;
t199 = qJD(1) * t297 - t536 * t563;
t446 = t536 * t592;
t948 = t199 * t446;
t934 = t589 * t603;
t920 = t591 * t609;
t749 = -rSges(7,1) * t593 - rSges(7,2) * t594;
t423 = t749 * t589;
t211 = t598 * t423 + (t591 * t750 + t571) * qJD(3);
t279 = qJD(3) * t657 - t589 * t831;
t906 = -t211 - t279;
t899 = -t590 * t329 + t592 * t330;
t345 = rSges(5,1) * t925 - rSges(5,2) * t935 + t573;
t896 = -t330 - t345;
t894 = t331 + t359;
t515 = pkin(8) * t932;
t416 = t515 - t843;
t517 = pkin(8) * t925;
t418 = -pkin(4) * t935 + t517;
t891 = t416 * t563 + t418 * t859;
t251 = -t390 * rSges(7,1) - t391 * rSges(7,2);
t252 = t392 * rSges(7,1) - t393 * rSges(7,2);
t887 = -t533 * t592 - t381;
t884 = -t531 * t592 + t385;
t820 = t604 * t864;
t526 = pkin(3) * t820;
t882 = t480 * t864 + t526;
t835 = t589 * t1002;
t881 = rSges(7,3) * t932 + t590 * t835;
t880 = rSges(7,3) * t925 + t592 * t835;
t879 = -t469 + t472;
t878 = t471 + t723;
t836 = rSges(6,2) * t934;
t877 = rSges(6,3) * t932 + t590 * t836;
t876 = rSges(6,3) * t925 + t592 * t836;
t875 = rSges(5,2) * t822 + rSges(5,3) * t862;
t873 = rSges(4,2) * t820 + rSges(4,3) * t862;
t871 = -t531 + t534;
t870 = t533 + t724;
t866 = qJD(1) * t468;
t865 = qJD(1) * t530;
t858 = qJD(3) * t607;
t852 = -m(5) - m(6) - m(7);
t845 = pkin(3) * t916;
t842 = pkin(5) * t934;
t838 = t589 * t1006;
t837 = t589 * t1005;
t832 = pkin(3) * t858;
t825 = t592 * t243 + t590 * t244 - t329 * t862;
t812 = t936 / 0.2e1;
t811 = t935 / 0.2e1;
t809 = -pkin(2) - t1007;
t803 = t861 / 0.2e1;
t802 = -t563 / 0.2e1;
t800 = -t859 / 0.2e1;
t793 = t604 * (-t590 ^ 2 - t1110);
t792 = t335 * t251 + t252 * t336;
t791 = t1059 * t252 - t335 * t423;
t790 = -t1059 * t251 - t336 * t423;
t785 = t520 - t933;
t783 = -qJD(1) * t416 + t526;
t775 = t590 * t417 + t592 * t419 + t899;
t762 = -t386 + t794;
t454 = t1058 * qJD(3);
t761 = -t454 - t832;
t760 = -t455 - t832;
t758 = t597 + t784;
t539 = rSges(2,1) * t608 - rSges(2,2) * t605;
t537 = rSges(2,1) * t605 + rSges(2,2) * t608;
t538 = -rSges(4,2) * t604 + t1007;
t748 = -t920 - t1018;
t745 = -t590 * t61 - t592 * t60;
t742 = t590 * t76 - t592 * t75;
t741 = t590 * t75 + t592 * t76;
t740 = t590 * t78 - t592 * t77;
t739 = t590 * t77 + t592 * t78;
t738 = t590 * t82 - t592 * t81;
t737 = t590 * t81 + t592 * t82;
t736 = t590 * t84 - t592 * t83;
t735 = t590 * t83 + t592 * t84;
t85 = t1111 + t664;
t734 = t592 * t85 + t989;
t733 = t590 * t88 - t592 * t87;
t732 = t590 * t87 + t592 * t88;
t731 = t590 * t95 - t592 * t94;
t730 = t590 * t94 + t592 * t95;
t712 = -t199 * t590 - t949;
t709 = -t236 * t592 - t237 * t590;
t276 = -rSges(4,2) * t592 * t858 + (-t607 * t864 - t813) * rSges(4,1) + t873;
t445 = t536 * t590;
t277 = -qJD(3) * t445 + (t538 * t592 + t574) * qJD(1);
t708 = t276 * t592 + t277 * t590;
t701 = t387 * t590 + t389 * t592;
t409 = t428 * pkin(5);
t694 = t794 - t894;
t693 = (t418 - t845) * qJD(1);
t388 = rSges(6,1) * t921 - rSges(6,2) * t922 + t572;
t360 = rSges(7,1) * t923 - rSges(7,2) * t924 + t571;
t691 = -t275 + t760;
t688 = t1019 * t589 - t581;
t687 = t592 * t273 + t590 * t274 + t417 * t862 + t825;
t686 = -t520 - t585 - t571;
t410 = t476 * t590;
t673 = -t920 + t1061;
t672 = t60 * t331 - t48 * t907;
t671 = t760 + t906;
t408 = t426 * pkin(5);
t656 = -t224 * t443 + t226 * t442 + t374 * t546;
t655 = (-Icges(7,5) * t390 - Icges(7,6) * t391) * t336 - (Icges(7,5) * t392 - Icges(7,6) * t393) * t335 - t420 * t1059;
t654 = (-Icges(6,5) * t426 - Icges(6,6) * t427) * t443 - (Icges(6,5) * t428 - Icges(6,6) * t429) * t442 - t433 * t546;
t653 = t339 * t592 - t340 * t590;
t652 = -t604 * t884 + t607 * t887;
t651 = t592 * t663;
t650 = t589 * t655;
t649 = t589 * t654;
t638 = (-t589 * t878 + t591 * t879) * qJD(1);
t637 = (-t604 * t870 + t607 * t871) * qJD(1);
t632 = (Icges(7,1) * t392 - t190 - t971) * t335 - (-Icges(7,1) * t390 - t188 - t364) * t336 + (-t355 + t422) * t1059;
t630 = (Icges(6,1) * t428 - t229 - t974) * t442 - (-Icges(6,1) * t426 - t227 - t406) * t443 + (-t378 + t439) * t546;
t16 = t101 * t392 + t103 * t393 + t169 * t188 - t170 * t192 + t185 * t660 + t935 * t99;
t17 = t100 * t392 + t102 * t393 + t169 * t190 + t170 * t193 + t187 * t660 + t935 * t98;
t18 = -t101 * t390 + t103 * t391 + t171 * t188 - t172 * t192 + t185 * t661 + t936 * t99;
t19 = -t100 * t390 + t102 * t391 + t171 * t190 + t172 * t193 + t187 * t661 + t936 * t98;
t43 = t1059 * t140 + t335 * t88 - t336 * t87;
t46 = t169 * t355 + t170 * t357 + t202 * t935 + t203 * t392 + t204 * t393 + t353 * t660;
t47 = t171 * t355 + t172 * t357 + t202 * t936 - t203 * t390 + t204 * t391 + t353 * t661;
t5 = t1059 * t46 + t121 * t328 + t153 * t78 + t154 * t77 - t16 * t336 + t17 * t335;
t6 = t1059 * t47 + t120 * t328 + t153 * t76 + t154 * t75 - t18 * t336 + t19 * t335;
t627 = ((qJD(3) * t739 - t46) * t591 + (-qJD(1) * t740 + qJD(3) * t121 + t16 * t590 + t17 * t592) * t589) * t1036 + (t1045 * t392 + t632 * t393 - t592 * t650) * t1037 + (-t1045 * t390 + t391 * t632 - t590 * t650) * t1034 + ((qJD(3) * t741 - t47) * t591 + (-qJD(1) * t742 + qJD(3) * t120 + t18 * t590 + t19 * t592) * t589) * t1035 + (t655 * t591 + (-t1045 * t593 + t594 * t632) * t589) * t1026 + t6 * t812 + (-t121 * t591 + t589 * t739) * t1042 + (-t120 * t591 + t589 * t741) * t1041 + t5 * t811 + t43 * t803 + (-t140 * t591 + t589 * t732) * t1038 + ((qJD(3) * t732 - t53) * t591 + (-qJD(1) * t733 + qJD(3) * t140 + t22 * t590 + t23 * t592) * t589) * t1025 + (t980 + t985 + t986 + t999 - t1000) * t1021 + t1077 * t37 + t1076 * t36;
t72 = -t236 * t442 + t237 * t443 + t729;
t626 = t72 * t709 + (t590 * t85 - t592 * t86) * t386;
t13 = -t105 * t1059 - t546 * t125 + t154 * t359 + t196 * t328 - t336 * t211 + t450 * t214 + t255 * t331 - t443 * t279 + t614;
t618 = -t11 * t196 * t935 + t13 * (-t196 * t591 + t359 * t936) + t48 * (t105 * t935 - t196 * t817) + t60 * (t591 * t105 + t211 * t936 + t359 * t661) + t61 * (t197 * t861 + t359 * t822);
t616 = -t1048 * t589 + t653 * t591;
t615 = t1044 * t589;
t612 = (t1059 * t353 - t185 * t336 + t187 * t335) * t591 + t1046 * t589;
t496 = t538 * qJD(3);
t466 = t589 * t598;
t411 = t476 * t592;
t398 = t592 * t780;
t397 = t590 * t780;
t332 = t785 - t1066;
t311 = -t592 * t838 + t876;
t310 = -t590 * t838 + t877;
t309 = t382 * t592;
t308 = t382 * t590;
t307 = t378 * t592;
t306 = t378 * t590;
t295 = -t592 * t837 + t880;
t294 = -t590 * t837 + t881;
t293 = t357 * t592;
t292 = t357 * t590;
t291 = t355 * t592;
t290 = t355 * t590;
t287 = rSges(6,1) * t428 - rSges(6,2) * t429;
t286 = -rSges(6,1) * t426 - rSges(6,2) * t427;
t261 = t592 * t673 - t517;
t260 = t590 * t673 - t515;
t223 = -qJD(3) * t410 + (t1058 * t592 + t573) * qJD(1);
t222 = rSges(5,1) * t662 - rSges(5,2) * t817 + t875;
t194 = qJD(3) * t701 + qJD(2);
t134 = -t476 * t563 + (t345 + t765) * qJD(1) - t872;
t111 = (t344 * t590 + t345 * t592) * qJD(3) + t810;
t109 = qJD(1) * t276 + qJDD(1) * t389 - t457 * t536 - t496 * t563 + t674;
t108 = -t840 - t496 * t859 + t458 * t536 + (-t277 - t456) * qJD(1) + t764 * qJDD(1);
t93 = qJD(3) * t708 + t387 * t457 - t389 * t458 + qJDD(2);
t69 = qJD(1) * t222 + qJDD(1) * t345 - t454 * t563 - t457 * t476 + t617;
t68 = t458 * t476 + (-qJD(3) * t454 - t844) * t592 + (-t223 + t903) * qJD(1) + t696 * qJDD(1) + t675;
t51 = t344 * t457 + t896 * t458 + (t222 * t592 + t223 * t590) * qJD(3) + t769;
t50 = t207 * t378 + t208 * t382 + t262 * t936 - t265 * t426 + t268 * t427 + t374 * t661;
t49 = t205 * t378 + t206 * t382 + t262 * t935 + t265 * t428 + t268 * t429 + t374 * t660;
t45 = t142 * t546 + t442 * t95 - t443 * t94;
t27 = t112 * t936 - t114 * t426 + t116 * t427 + t207 * t229 + t208 * t232 + t226 * t661;
t26 = t113 * t936 - t115 * t426 + t117 * t427 + t207 * t227 - t208 * t231 + t224 * t661;
t25 = t112 * t935 + t114 * t428 + t116 * t429 + t205 * t229 + t206 * t232 + t226 * t660;
t24 = t113 * t935 + t115 * t428 + t117 * t429 + t205 * t227 - t206 * t231 + t224 * t660;
t20 = t118 * t443 + t119 * t442 - t236 * t254 - t237 * t255 + t633;
t12 = t104 * t1059 + t546 * t124 - t153 * t359 + t328 * t197 - t335 * t211 + t450 * t215 - t254 * t331 - t442 * t279 + t613;
t10 = t128 * t450 + t254 * t82 + t255 * t81 - t26 * t443 + t27 * t442 + t50 * t546;
t9 = t129 * t450 - t24 * t443 + t25 * t442 + t254 * t84 + t255 * t83 + t49 * t546;
t1 = [(t60 * (t590 * t779 - t752 + t778 + t823) + t61 * (-t592 * t779 + t561 + t824 + t827) + ((t748 - t937) * t988 + (t591 * t981 + t937) * t990) * qJD(3) + ((-t60 * t608 - t605 * t61) * pkin(1) + (-t1017 * t60 + t61 * t686) * t590 + (t60 * (t686 + t933) - t61 * t602) * t592) * qJD(1) - (t1050 + t1095 - t60) * t61 + (t12 - g(2)) * (t690 + t758 + t197) + (t13 - g(1)) * (t936 * t981 + t1069 + t763 + t767)) * m(7) + (-(-t814 - t198 - t465 + (-t387 - t1013) * qJD(1)) * t199 + t199 * (t558 + t873) + (t536 * t950 - t948) * qJD(3) + ((-t198 * t608 - t199 * t605) * pkin(1) + (-pkin(2) - t538) * t949 + (t198 * (-rSges(4,3) - pkin(7)) + t199 * t809) * t590) * qJD(1) + (t109 - g(2)) * t297 + (t108 - g(1)) * (t590 * t809 - t1013 + t582 + t869)) * m(4) + (t1102 + t1080) * t1028 + (t1103 + t1081) * t1027 + (-(t1050 + t1111 - t85) * t86 + t85 * (t491 - t755 + t823) + t86 * (-pkin(4) * t819 + t492 + t692 + t826) + (t1019 * t85 * t860 + t38 * t688) * t590 + ((-t605 * t86 - t608 * t85) * pkin(1) + (t85 * (-t585 - t1066 - t572) - t86 * t602) * t592 + (-t585 + t688) * t989) * qJD(1) - t688 * t1016 + (t39 - g(2)) * (t758 + t237 + t419) + (t38 - g(1)) * (-t754 + t763)) * m(6) + (t1105 + t1106) * t801 + (t50 + t42) * t1031 + (t1104 - t1107 + t1108) * t800 + t37 * t1034 + t42 * t1030 + ((-t477 * t611 - g(2) + t768) * t449 + (-t840 + (-0.2e1 * t479 - t597 + t449) * t611 - g(1)) * t448) * m(3) + (t952 * t563 + (t69 - g(2)) * (t345 + t758) + (t68 - g(1)) * (-t1067 + t1120) + (t823 + (-t573 - t597 + (-t585 - t1058) * t592) * qJD(1)) * t133 + (t651 * qJD(3) - t1070 + t133 + t561 - t659 + t875 + (t1013 + t344 + t1120) * qJD(1)) * t134) * m(5) + ((t1098 * t590 + ((t1125 + t1133) * t592 + t1084 + t1130 + t1132) * t592) * qJD(3) + t1114) * t799 + t984 / 0.2e1 + t985 / 0.2e1 + t983 / 0.2e1 + (((t1099 * t592 + t1082 - t1098) * t592 + (t1099 * t590 + t1083 + t1131) * t590) * qJD(3) + t1109 - t1113) * t802 + (t47 + t37) * t1035 + t979 + t980 + t986 / 0.2e1 + t999 / 0.2e1 + t997 / 0.2e1 - t998 / 0.2e1 - m(2) * (-g(1) * t537 + g(2) * t539) - t1000 / 0.2e1 + (m(3) * (t448 ^ 2 + t479 * t449) + m(2) * (t537 ^ 2 + t539 ^ 2) + Icges(2,3) + Icges(3,3) - t1126) * qJDD(1) + (-t1127 * qJD(3) + t452 * t591 + t453 * t589 + t494 * t607 + t495 * t604) * qJD(1) + t49 * t1032 + t46 * t1036 + t128 * t1039 + t129 * t1040 + t120 * t1041 + t121 * t1042; m(3) * qJDD(2) + (-m(3) - m(4) + t852) * g(3) + m(4) * t93 + m(5) * t51 + m(6) * t20 + m(7) * t11; (t93 * t701 + t194 * ((t387 * t592 - t389 * t590) * qJD(1) + t708) + t712 * t496 + (-t108 * t592 - t109 * t590 + (-t199 * t592 + t950) * qJD(1)) * t536 - (t198 * t445 - t948) * qJD(1) - (t194 * (-t445 * t590 - t446 * t592) + t712 * t538) * qJD(3) + g(1) * t446 + g(2) * t445 - g(3) * t538) * m(4) + (t85 * t882 + t20 * t775 + t72 * t687 + (t38 * t762 + t85 * t691 + t20 * t237 + t72 * t118 + (-t236 * t72 + t762 * t86) * qJD(1)) * t592 + (t39 * t762 + t86 * t691 - t20 * t236 + t72 * t119 + (t85 * t386 + t72 * (-t237 + t895)) * qJD(1)) * t590 - t85 * (-t310 * t546 - t388 * t443 + t783) - t86 * (t311 * t546 - t388 * t442 + t693) - t72 * (t310 * t442 + t311 * t443 + t891) - ((t236 * t85 + t237 * t86) * t589 + t626 * t591) * qJD(5) - (-t734 * t1066 + (-t607 * t734 + t72 * t793) * pkin(3)) * qJD(3) - g(1) * (t517 - t845 + t876) - g(2) * (t515 - t846 + t877) - g(3) * (t388 + t596 + t1066) - t1057 * t589 * (-pkin(4) - t1006)) * m(6) + (-t60 * (-t1059 * t294 + t196 * t466 - t260 * t546 - t332 * t443 - t336 * t360 + t359 * t397 + t783) - t61 * (t1059 * t295 + t197 * t466 + t261 * t546 - t332 * t442 - t335 * t360 - t359 * t398 + t693) - t48 * (-t196 * t398 - t197 * t397 + t260 * t442 + t261 * t443 + t294 * t335 + t295 * t336 + t891) - ((t214 * t60 + t987) * t589 + (t48 * (-t214 * t592 - t215 * t590) + (-t988 + t990) * t331) * t591) * qJD(5) - (t745 * t1066 + (t48 * t793 + t607 * t745) * pkin(3)) * qJD(3) + t60 * t882 + t11 * t775 + t48 * t687 + (t13 * t694 + t60 * t671 + (t61 * t694 + t1073) * qJD(1) - t1049) * t592 + (t12 * t694 + t61 * t671 + t11 * t908 + t48 * (t105 + t125) + (t60 * t894 + t48 * (t895 - t907)) * qJD(1)) * t590 - g(1) * t880 - g(2) * t881 - g(3) * (t360 + t596 + t785) - t1057 * ((-t584 - t1005) * t589 + t748)) * m(7) - ((t1048 * t591 + t653 * t589 + (t590 * t884 - t592 * t885) * t607 + (t590 * t887 + t592 * t888) * t604) * qJD(3) + (t589 * t879 + t591 * t878 + t604 * t871 + t607 * t870) * qJD(1)) * qJD(1) / 0.2e1 + ((-t563 * t940 + t866) * t590 + (t638 + (t590 * t941 + t616) * qJD(3)) * t592 + (-t563 * t938 + t865) * t590 + (t637 + (-t1047 * t592 + (t939 + t652) * t590) * qJD(3)) * t592) * t802 + ((-t859 * t941 - t866) * t592 + (t638 + (t592 * t940 + t616) * qJD(3)) * t590 + (-t859 * t939 - t865) * t592 + (t637 + (t652 * t590 + (-t1047 + t938) * t592) * qJD(3)) * t590) * t799 - (t41 * t590 + t42 * t592) * t856 / 0.2e1 + (qJD(1) * t1104 + t1079 * qJD(3) + qJDD(1) * t1103 + t1084 * t457 + t1085 * t458 + t10 + t6) * t1020 + (qJD(1) * t1105 + t1078 * qJD(3) + qJDD(1) * t1102 + t1082 * t457 + t1083 * t458 + t5 + t9) * t1022 + (-(t133 * t410 + t134 * (-t411 - t845)) * qJD(1) - (t111 * pkin(3) * t793 + (-t1068 * t133 - t111 * t411) * t592 + (-t1068 * t134 - t111 * t410) * t590) * qJD(3) - g(1) * t651 - g(3) * t1068 - t663 * t1014 + t51 * t899 + t111 * t825 + (t68 * t663 + t133 * t761 + t51 * t345 + t111 * t222 + (t111 * t344 + t134 * t663) * qJD(1)) * t592 + (t69 * t663 + t134 * t761 + t51 * t344 + t111 * t223 + (t111 * t896 + t952) * qJD(1)) * t590) * m(5) + (t1080 * t590 - t1081 * t592) * qJDD(1) / 0.2e1 + ((t1082 * t592 + t1083 * t590) * qJD(1) + t1078) * t801 + ((t1084 * t592 + t1085 * t590) * qJD(1) + t1079) * t800 + (t42 + t37 + t1108) * t804 + (t41 + t36 + t1109) * t864 / 0.2e1 + t1116 * t1027 + t1117 * t1028 + ((-t291 * t392 - t293 * t393) * t335 + t78 * t398 - (-t290 * t392 - t292 * t393) * t336 + t77 * t397 + (t356 * t392 + t358 * t393) * t1059 + t121 * t466 + t612 * t592) * t1037 + (t140 * t466 + t87 * t397 + t88 * t398 - t1046 * t591 + ((t291 * t593 - t293 * t594 + t187) * t335 - (t290 * t593 - t292 * t594 + t185) * t336 + (-t356 * t593 + t358 * t594 + t353) * t1059) * t589) * t1026 + ((t291 * t390 - t293 * t391) * t335 + t76 * t398 - (t290 * t390 - t292 * t391) * t336 + t75 * t397 + (-t356 * t390 + t358 * t391) * t1059 + t120 * t466 + t612 * t590) * t1034 - t45 * t857 / 0.2e1 + (((t307 * t603 - t309 * t606 + t226) * t442 - (t306 * t603 - t308 * t606 + t224) * t443 + (-t379 * t603 + t383 * t606 + t374) * t546 + t142 * qJD(5)) * t589 + (qJD(5) * t730 - t1044) * t591) * t1024 + ((t307 * t426 - t309 * t427) * t442 - (t306 * t426 - t308 * t427) * t443 + (-t379 * t426 + t383 * t427) * t546 + (t128 * t589 + t82 * t925) * qJD(5) + ((qJD(5) * t81 + t656) * t591 + t615) * t590) * t1030 + ((-t307 * t428 - t309 * t429) * t442 - (-t306 * t428 - t308 * t429) * t443 + (t379 * t428 + t383 * t429) * t546 + (t129 * t589 + t83 * t932) * qJD(5) + ((qJD(5) * t84 + t656) * t591 + t615) * t592) * t1033 - t466 * t43 / 0.2e1 + (qJD(1) * t730 - t28 * t592 + t29 * t590) * t1023 + (qJD(1) * t732 - t22 * t592 + t23 * t590) * t1025 + (t1107 * t592 + t1106 * t590 + (t1080 * t592 + t1081 * t590) * qJD(1)) * qJD(1) / 0.2e1 + t731 * t1029 + (qJD(1) * t737 - t26 * t592 + t27 * t590) * t1031 + (qJD(1) * t735 - t24 * t592 + t25 * t590) * t1032 + (qJD(1) * t741 - t18 * t592 + t19 * t590) * t1035 + (qJD(1) * t739 - t16 * t592 + t17 * t590) * t1036 + t733 * t1038 + t738 * t1039 + t736 * t1040 + t742 * t1041 + t740 * t1042 - t397 * t36 / 0.2e1 - t398 * t37 / 0.2e1; t852 * (-g(2) * t592 + t1016) + 0.2e1 * (t1020 * t12 + t1022 * t13) * m(7) + 0.2e1 * (t994 / 0.2e1 - t993 / 0.2e1) * m(6) + 0.2e1 * (t1020 * t69 + t1022 * t68) * m(5); (t654 * t591 + (-t1043 * t603 + t606 * t630) * t589) * t1024 + (-t129 * t591 + t589 * t735) * t1040 + (t979 + t983 + t984 + t997 - t998) * t1021 + t9 * t811 + (-t1043 * t426 + t427 * t630 - t590 * t649) * t1030 + ((qJD(3) * t730 - t63) * t591 + (-qJD(1) * t731 + qJD(3) * t142 + t28 * t590 + t29 * t592) * t589) * t1023 + (-t142 * t591 + t589 * t730) * t1029 + t627 + ((qJD(3) * t735 - t49) * t591 + (-qJD(1) * t736 + qJD(3) * t129 + t24 * t590 + t25 * t592) * t589) * t1032 + (-t128 * t591 + t589 * t737) * t1039 + t45 * t803 + t10 * t812 + ((qJD(3) * t737 - t50) * t591 + (-qJD(1) * t738 + qJD(3) * t128 + t26 * t590 + t27 * t592) * t589) * t1031 + (t1043 * t428 + t630 * t429 - t592 * t649) * t1033 + t1077 * t42 + t1076 * t41 + (-t60 * (t408 * t546 + t443 * t842 + t790) - t61 * (t409 * t546 + t442 * t842 + t791) - t48 * (-t408 * t442 + t409 * t443 + t792) + (-t13 * t214 + t60 * t125 - t12 * t907 - t61 * t911 + ((-t48 * t214 - t61 * t894) * t592 + t672 * t590) * qJD(3)) * t591 + t618 - g(1) * (t409 + t252) - g(2) * (-t408 + t251) + ((-t60 * t908 + t987) * qJD(3) + (qJD(1) * t672 - t11 * t214 - t12 * t894 + t48 * t125 + t61 * t906) * t592 + (t13 * t331 + t60 * t279 + (t61 * t331 - t1073) * qJD(1) + t1049) * t590 - g(3) * (t749 - t1017)) * t589) * m(7) + (-g(1) * t287 - g(2) * t286 - g(3) * t444 - t85 * (-t286 * t546 - t443 * t444) - t86 * (t287 * t546 - t442 * t444) - t72 * (t286 * t442 + t287 * t443) + (qJD(3) * t626 - t86 * t118 + t85 * t119 - t236 * t38 - t39 * t237) * t591 + (t85 * (qJD(3) * t236 + t275 * t590) + t86 * (qJD(3) * t237 - t275 * t592) + t20 * t709 + t72 * (-t118 * t590 + t119 * t592 + t236 * t864 - t237 * t862) + (qJD(1) * t734 - t993 + t994) * t386) * t589) * m(6); t627 + (-t48 * t792 - t60 * t790 - t61 * t791 + (t61 * (-t359 * t859 - t104) + (-t48 * t563 - t12) * t197) * t591 + (t60 * qJD(3) * t196 + (-qJD(1) * t197 * t48 - t12 * t359 - t211 * t61) * t592 + (-t11 * t197 + t48 * (qJD(1) * t196 - t104)) * t590) * t589 + t618 - g(1) * t252 - g(2) * t251 - g(3) * t423) * m(7);];
tau  = t1;
