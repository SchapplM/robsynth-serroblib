% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP5_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP5_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:14:32
% EndTime: 2019-03-09 03:15:59
% DurationCPUTime: 75.93s
% Computational Cost: add. (44955->1261), mult. (49184->1607), div. (0->0), fcn. (45942->10), ass. (0->605)
t1047 = Icges(6,1) + Icges(7,1);
t1067 = Icges(7,4) + Icges(6,5);
t1066 = Icges(6,6) - Icges(7,6);
t546 = pkin(10) + qJ(5);
t530 = sin(t546);
t1081 = (Icges(6,4) - Icges(7,5)) * t530;
t1046 = Icges(6,2) + Icges(7,3);
t1080 = Icges(7,2) + Icges(6,3);
t532 = cos(t546);
t1079 = -t1066 * t530 + t1067 * t532;
t1078 = t1047 * t532 - t1081;
t547 = pkin(9) + qJ(3);
t531 = sin(t547);
t533 = cos(t547);
t889 = Icges(6,4) * t532;
t664 = -Icges(6,2) * t530 + t889;
t1077 = -t1066 * t533 + t531 * t664;
t555 = cos(qJ(1));
t851 = t555 * t530;
t554 = sin(qJ(1));
t852 = t554 * t532;
t417 = t533 * t851 - t852;
t856 = t533 * t555;
t780 = t532 * t856;
t418 = t530 * t554 + t780;
t858 = t531 * t555;
t204 = Icges(6,5) * t418 - Icges(6,6) * t417 + Icges(6,3) * t858;
t207 = Icges(7,4) * t418 + Icges(7,2) * t858 + Icges(7,6) * t417;
t1044 = t204 + t207;
t884 = Icges(7,5) * t417;
t213 = Icges(7,1) * t418 + Icges(7,4) * t858 + t884;
t399 = Icges(6,4) * t417;
t216 = Icges(6,1) * t418 + Icges(6,5) * t858 - t399;
t1055 = t213 + t216;
t396 = Icges(7,5) * t418;
t201 = Icges(7,6) * t858 + Icges(7,3) * t417 + t396;
t891 = Icges(6,4) * t418;
t210 = -Icges(6,2) * t417 + Icges(6,6) * t858 + t891;
t1057 = t201 - t210;
t1018 = t1044 * t858 + t1055 * t418 + t1057 * t417;
t857 = t533 * t554;
t415 = t530 * t857 + t532 * t555;
t416 = t533 * t852 - t851;
t859 = t531 * t554;
t202 = Icges(6,5) * t416 - Icges(6,6) * t415 + Icges(6,3) * t859;
t205 = Icges(7,4) * t416 + Icges(7,2) * t859 + Icges(7,6) * t415;
t1011 = t202 + t205;
t395 = Icges(7,5) * t416;
t200 = -Icges(7,6) * t859 - Icges(7,3) * t415 - t395;
t398 = Icges(6,4) * t416;
t208 = -Icges(6,2) * t415 + Icges(6,6) * t859 + t398;
t1012 = t200 + t208;
t394 = Icges(7,5) * t415;
t211 = Icges(7,1) * t416 + Icges(7,4) * t859 + t394;
t397 = Icges(6,4) * t415;
t215 = -Icges(6,1) * t416 - Icges(6,5) * t859 + t397;
t1056 = t211 - t215;
t1019 = t1011 * t858 - t1012 * t417 + t1056 * t418;
t790 = qJD(5) * t555;
t798 = qJD(3) * t554;
t453 = t531 * t790 + t798;
t791 = qJD(5) * t554;
t797 = qJD(3) * t555;
t454 = -t531 * t791 + t797;
t792 = qJD(5) * t533;
t508 = qJD(1) - t792;
t1042 = -t1067 * t533 + t1078 * t531;
t1043 = t1079 * t531 - t1080 * t533;
t860 = t531 * t532;
t474 = Icges(7,5) * t860;
t861 = t530 * t531;
t1054 = Icges(7,3) * t861 - t1077 + t474;
t980 = t1042 * t418 + t1043 * t858 + t1054 * t417;
t991 = t1018 * t453 - t1019 * t454 + t508 * t980;
t1020 = t1044 * t859 + t1055 * t416 + t1057 * t415;
t1021 = t1011 * t859 - t1012 * t415 + t1056 * t416;
t981 = t1042 * t416 + t1043 * t859 + t1054 * t415;
t992 = t1020 * t453 - t1021 * t454 + t981 * t508;
t882 = Icges(7,5) * t532;
t660 = Icges(7,3) * t530 + t882;
t974 = (t660 - t664) * t533 - t1066 * t531;
t1076 = t1079 * t533 + t1080 * t531;
t973 = t1067 * t531 + t1078 * t533;
t1075 = (t1046 * t532 + t1081) * t531;
t1074 = (-t1066 * t532 - t1067 * t530) * t531;
t759 = t531 * t797;
t180 = qJD(1) * t415 - qJD(5) * t780 + (t759 - t791) * t530;
t644 = t530 * t508;
t803 = qJD(1) * t533;
t716 = -qJD(5) + t803;
t181 = t555 * t644 + (-t554 * t716 - t759) * t532;
t757 = t533 * t797;
t802 = qJD(1) * t554;
t763 = t531 * t802;
t424 = t757 - t763;
t101 = Icges(7,4) * t181 + Icges(7,2) * t424 - Icges(7,6) * t180;
t99 = Icges(6,5) * t181 + Icges(6,6) * t180 + Icges(6,3) * t424;
t1073 = t101 + t99;
t103 = Icges(6,4) * t181 + Icges(6,2) * t180 + Icges(6,6) * t424;
t97 = Icges(7,5) * t181 + Icges(7,6) * t424 - Icges(7,3) * t180;
t1072 = -t103 + t97;
t760 = t531 * t798;
t801 = qJD(1) * t555;
t182 = (t533 * t791 - t802) * t532 + (t533 * t801 - t760 - t790) * t530;
t183 = t554 * t644 + (t555 * t716 - t760) * t532;
t758 = t533 * t798;
t423 = t531 * t801 + t758;
t104 = Icges(6,4) * t183 - Icges(6,2) * t182 + Icges(6,6) * t423;
t98 = Icges(7,5) * t183 + Icges(7,6) * t423 + Icges(7,3) * t182;
t1071 = -t104 + t98;
t550 = cos(pkin(10));
t853 = t550 * t555;
t548 = sin(pkin(10));
t855 = t548 * t554;
t443 = t533 * t855 + t853;
t850 = t555 * t548;
t854 = t550 * t554;
t444 = t533 * t854 - t850;
t244 = Icges(5,5) * t444 - Icges(5,6) * t443 + Icges(5,3) * t859;
t875 = Icges(4,3) * t555;
t373 = Icges(4,5) * t857 - Icges(4,6) * t859 - t875;
t500 = Icges(4,4) * t859;
t887 = Icges(4,5) * t555;
t377 = Icges(4,1) * t857 - t500 - t887;
t879 = Icges(4,6) * t555;
t375 = Icges(4,4) * t857 - Icges(4,2) * t859 - t879;
t865 = t375 * t531;
t647 = -t377 * t533 + t865;
t247 = Icges(5,4) * t444 - Icges(5,2) * t443 + Icges(5,6) * t859;
t250 = Icges(5,1) * t444 - Icges(5,4) * t443 + Icges(5,5) * t859;
t653 = -t247 * t443 + t250 * t444;
t1070 = t244 * t859 - t373 * t555 - t554 * t647 + t653;
t520 = Icges(4,4) * t533;
t666 = -Icges(4,2) * t531 + t520;
t376 = Icges(4,6) * t554 + t555 * t666;
t892 = Icges(4,4) * t531;
t462 = Icges(4,1) * t533 - t892;
t378 = Icges(4,5) * t554 + t462 * t555;
t331 = t378 * t857;
t458 = Icges(4,5) * t533 - Icges(4,6) * t531;
t374 = Icges(4,3) * t554 + t458 * t555;
t721 = t374 * t555 - t331;
t142 = -t376 * t859 - t721;
t445 = -t533 * t850 + t854;
t446 = t533 * t853 + t855;
t246 = Icges(5,5) * t446 + Icges(5,6) * t445 + Icges(5,3) * t858;
t249 = Icges(5,4) * t446 + Icges(5,2) * t445 + Icges(5,6) * t858;
t252 = Icges(5,1) * t446 + Icges(5,4) * t445 + Icges(5,5) * t858;
t79 = t246 * t859 - t443 * t249 + t444 * t252;
t1069 = t142 + t79;
t833 = t554 * t374 + t378 * t856;
t144 = -t376 * t858 + t833;
t81 = t246 * t858 + t445 * t249 + t446 * t252;
t1068 = t144 + t81;
t100 = Icges(6,5) * t183 - Icges(6,6) * t182 + Icges(6,3) * t423;
t102 = Icges(7,4) * t183 + Icges(7,2) * t423 + Icges(7,6) * t182;
t1065 = t100 + t102;
t105 = Icges(7,1) * t181 + Icges(7,4) * t424 - Icges(7,5) * t180;
t107 = Icges(6,1) * t181 + Icges(6,4) * t180 + Icges(6,5) * t424;
t1064 = t105 + t107;
t106 = Icges(7,1) * t183 + Icges(7,4) * t423 + Icges(7,5) * t182;
t108 = Icges(6,1) * t183 - Icges(6,4) * t182 + Icges(6,5) * t423;
t1063 = t106 + t108;
t662 = Icges(5,5) * t550 - Icges(5,6) * t548;
t361 = -Icges(5,3) * t533 + t531 * t662;
t665 = Icges(5,4) * t550 - Icges(5,2) * t548;
t363 = -Icges(5,6) * t533 + t531 * t665;
t669 = Icges(5,1) * t550 - Icges(5,4) * t548;
t365 = -Icges(5,5) * t533 + t531 * t669;
t459 = Icges(4,2) * t533 + t892;
t461 = Icges(4,1) * t531 + t520;
t645 = t459 * t531 - t461 * t533;
t457 = Icges(4,5) * t531 + Icges(4,6) * t533;
t862 = t457 * t555;
t1062 = t361 * t859 - t363 * t443 + t365 * t444 - t554 * t645 - t862;
t863 = t457 * t554;
t1061 = t361 * t858 + t363 * t445 + t365 * t446 - t555 * t645 + t863;
t1060 = qJD(3) * t974 + qJD(5) * t1075;
t1059 = qJD(3) * t1076 + qJD(5) * t1074;
t384 = (-Icges(6,1) * t530 - t889) * t531;
t793 = qJD(5) * t531;
t1058 = (-Icges(7,1) * t530 + t882) * t793 + qJD(5) * t384 + t973 * qJD(3);
t1053 = t1042 * t532 + t1054 * t530;
t1017 = rSges(7,3) + qJ(6);
t1030 = rSges(7,1) + pkin(5);
t1052 = -t1017 * t415 - t1030 * t416;
t655 = -t210 * t530 + t216 * t532;
t656 = -t208 * t530 - t215 * t532;
t657 = t201 * t530 + t213 * t532;
t658 = -t200 * t530 + t211 * t532;
t1051 = -(-t1043 * t554 - t656 - t658) * t454 + (-t1043 * t555 - t655 - t657) * t453 + (-t1053 + t1076) * t508;
t999 = t1011 * t424 + t1012 * t180 + t1056 * t181 + t1063 * t418 + t1065 * t858 + t1071 * t417;
t997 = t1044 * t424 + t1055 * t181 - t1057 * t180 + t1064 * t418 + t1072 * t417 + t1073 * t858;
t996 = t1011 * t423 - t1012 * t182 + t1056 * t183 + t1063 * t416 + t1065 * t859 + t1071 * t415;
t995 = t1044 * t423 + t1055 * t183 + t1057 * t182 + t1064 * t416 + t1072 * t415 + t1073 * t859;
t1027 = t1042 * t181 + t1043 * t424 - t1054 * t180 + t1058 * t418 + t1059 * t858 + t1060 * t417;
t1026 = t1042 * t183 + t1043 * t423 + t1054 * t182 + t1058 * t416 + t1059 * t859 + t1060 * t415;
t1050 = (qJD(3) * t1053 - t1059) * t533 + (t1058 * t532 + t1060 * t530 + (-t1042 * t530 + t1054 * t532) * qJD(5) + t1043 * qJD(3)) * t531;
t74 = -t205 * t533 + t531 * t658;
t76 = -t202 * t533 + t531 * t656;
t1049 = t74 + t76;
t75 = -t207 * t533 + t531 * t657;
t77 = -t204 * t533 + t531 * t655;
t1048 = t75 + t77;
t553 = -pkin(7) - qJ(2);
t526 = t555 * t553;
t1045 = t1061 * qJD(1);
t979 = -t1043 * t533 + t1053 * t531;
t685 = pkin(5) * t532 + qJ(6) * t530;
t686 = rSges(7,1) * t532 + rSges(7,3) * t530;
t830 = -rSges(7,2) * t533 + (t685 + t686) * t531;
t1041 = -t531 * t660 + t1077;
t1040 = -t1074 * t508 + (-t1066 * t416 - t1067 * t415) * t454 + (t1066 * t418 + t1067 * t417) * t453;
t834 = -t554 * t373 - t377 * t856;
t143 = -t375 * t858 - t834;
t1004 = t445 * t247 + t446 * t250;
t80 = t244 * t858 + t1004;
t903 = t555 * t80;
t1039 = (t1068 * t554 - t143 * t555 - t903) * qJD(3);
t1038 = (t1069 * t554 - t1070 * t555) * qJD(3);
t1034 = t1062 * qJD(1);
t386 = qJD(6) * t417;
t847 = rSges(7,2) * t859 - t1052;
t1031 = t508 * t847 - t386;
t1029 = t1034 + t1038;
t1028 = t1039 + t1045;
t1003 = t247 * t548 - t250 * t550;
t328 = qJD(1) * t445 + t548 * t760;
t329 = qJD(1) * t446 - t550 * t760;
t134 = Icges(5,5) * t329 + Icges(5,6) * t328 + Icges(5,3) * t423;
t136 = Icges(5,4) * t329 + Icges(5,2) * t328 + Icges(5,6) * t423;
t138 = Icges(5,1) * t329 + Icges(5,4) * t328 + Icges(5,5) * t423;
t614 = qJD(3) * t459;
t256 = qJD(1) * t376 - t554 * t614;
t615 = qJD(3) * t461;
t258 = qJD(1) * t378 - t554 * t615;
t1025 = (t134 - t256) * t533 + (t136 * t548 - t138 * t550 - t258) * t531 + (t1003 * t533 - t244 * t531 + t647) * qJD(3);
t326 = qJD(1) * t443 + t548 * t759;
t327 = -qJD(1) * t444 - t550 * t759;
t133 = Icges(5,5) * t327 + Icges(5,6) * t326 + Icges(5,3) * t424;
t135 = Icges(5,4) * t327 + Icges(5,2) * t326 + Icges(5,6) * t424;
t137 = Icges(5,1) * t327 + Icges(5,4) * t326 + Icges(5,5) * t424;
t255 = -t555 * t614 + (-t554 * t666 + t879) * qJD(1);
t257 = -t555 * t615 + (-t462 * t554 + t887) * qJD(1);
t864 = t376 * t531;
t646 = -t378 * t533 + t864;
t651 = -t249 * t548 + t252 * t550;
t1024 = (-t133 + t255) * t533 + (-t135 * t548 + t137 * t550 + t257) * t531 + (t246 * t531 + t533 * t651 - t646) * qJD(3);
t362 = Icges(5,3) * t531 + t533 * t662;
t334 = t362 * qJD(3);
t364 = Icges(5,6) * t531 + t533 * t665;
t335 = t364 * qJD(3);
t366 = Icges(5,5) * t531 + t533 * t669;
t336 = t366 * qJD(3);
t448 = t666 * qJD(3);
t449 = t462 * qJD(3);
t562 = qJD(1) * t457 - t448 * t531 + t449 * t533 + (-t459 * t533 - t461 * t531) * qJD(3);
t945 = t645 * qJD(1) + t458 * qJD(3);
t1023 = t326 * t363 + t327 * t365 + t334 * t858 + t335 * t445 + t336 * t446 + t361 * t424 + t554 * t945 + t562 * t555;
t1022 = t328 * t363 + t329 * t365 + t334 * t859 - t335 * t443 + t336 * t444 + t361 * t423 + t562 * t554 - t555 * t945;
t1016 = t143 + t80;
t194 = t375 * t533 + t377 * t531;
t1015 = -t1003 * t531 - t244 * t533 + t194;
t195 = t376 * t533 + t378 * t531;
t1014 = -t246 * t533 + t531 * t651 + t195;
t385 = qJD(6) * t415;
t509 = pkin(3) * t856;
t439 = qJ(4) * t858 + t509;
t518 = pkin(4) * t855;
t522 = pkin(4) * t550 + pkin(3);
t552 = -pkin(8) - qJ(4);
t642 = t522 * t856 - t552 * t858 + t518;
t243 = t642 - t439;
t923 = pkin(3) - t522;
t746 = t923 * t531;
t849 = qJ(4) + t552;
t597 = -t533 * t849 + t746;
t465 = pkin(3) * t531 - qJ(4) * t533;
t795 = qJD(4) * t554;
t488 = t531 * t795;
t535 = qJD(2) * t555;
t807 = t535 - t488;
t765 = -t465 * t798 - t807;
t551 = cos(pkin(9));
t523 = pkin(2) * t551 + pkin(1);
t502 = t555 * t523;
t719 = -t553 * t554 + t502;
t769 = t439 + t719;
t578 = qJD(1) * (t243 + t769) + t597 * t798 + t765;
t845 = rSges(7,2) * t858 + t1017 * t417 + t1030 * t418;
t48 = -t453 * t830 + t845 * t508 + t385 + t578;
t729 = t48 * t830;
t789 = qJD(6) * t530;
t473 = t531 * t789;
t783 = pkin(4) * t850;
t702 = -t522 * t857 + t783;
t928 = pkin(3) * t533;
t984 = t531 * t849;
t242 = (t928 + t984) * t554 + t702;
t872 = qJ(4) * t531;
t467 = t872 + t928;
t436 = t467 * t554;
t796 = qJD(4) * t533;
t696 = t436 * t798 + t439 * t797 - t796;
t616 = -t242 * t798 + t243 * t797 + t696;
t42 = t453 * t847 + t454 * t845 + t473 + t616;
t731 = t42 * t847;
t1008 = t729 - t731;
t747 = t923 * t533;
t1007 = t747 + t872;
t689 = -rSges(6,1) * t416 + rSges(6,2) * t415;
t221 = rSges(6,3) * t859 - t689;
t688 = rSges(6,1) * t532 - rSges(6,2) * t530;
t352 = -rSges(6,3) * t533 + t531 * t688;
t1006 = -t221 * t508 - t352 * t454;
t648 = -t363 * t548 + t365 * t550;
t1005 = qJD(3) * (-t1003 * t555 - t554 * t651) + (t362 - t648) * qJD(1);
t1002 = 0.2e1 * qJD(3);
t784 = qJD(3) * qJD(5);
t741 = t533 * t784;
t357 = qJD(1) * t453 + t554 * t741;
t358 = qJD(1) * t454 + t555 * t741;
t742 = t531 * t784;
t1001 = t1018 * t358 + t1019 * t357 + t1027 * t508 + t453 * t997 - t454 * t999 + t742 * t980;
t1000 = t1020 * t358 + t1021 * t357 + t1026 * t508 + t453 * t995 - t454 * t996 + t742 * t981;
t479 = qJ(4) * t757;
t815 = qJD(1) * t783 + t552 * t763;
t131 = -t479 + (-t533 * t552 + t746) * t797 + t1007 * t802 + t815;
t484 = pkin(3) * t760;
t766 = -t423 * t552 - t522 * t760;
t132 = -qJ(4) * t758 + t484 + (-t1007 * t555 + t518) * qJD(1) + t766;
t794 = qJD(4) * t555;
t490 = t531 * t794;
t606 = -t533 * t802 - t759;
t238 = pkin(3) * t606 - qJ(4) * t763 + t479 + t490;
t239 = qJ(4) * t423 + qJD(1) * t509 - t484 + t488;
t786 = qJD(1) * qJD(3);
t744 = t555 * t786;
t785 = qJD(3) * qJD(4);
t714 = t238 * t797 + t239 * t798 + t436 * t744 + t531 * t785;
t610 = t131 * t797 + t132 * t798 - t242 * t744 + t714;
t844 = -t243 - t439;
t698 = t844 * t802;
t752 = t533 * t789;
t753 = t532 * t793;
t989 = t1017 * t182 + t1030 * t183 + t385;
t896 = rSges(7,2) * t423 + t989;
t986 = rSges(7,2) * t757 - t1017 * t180 + t1030 * t181 + t386;
t897 = -rSges(7,2) * t763 + t986;
t5 = qJD(6) * t753 + t897 * t454 + t896 * t453 + t847 * t358 - t845 * t357 + (t698 + t752) * qJD(3) + t610;
t998 = t5 * t845;
t17 = (qJD(3) * t658 - t102) * t533 + (qJD(3) * t205 + t106 * t532 + t530 * t98 + (-t200 * t532 - t211 * t530) * qJD(5)) * t531;
t19 = (qJD(3) * t656 - t100) * t533 + (qJD(3) * t202 - t104 * t530 + t108 * t532 + (-t208 * t532 + t215 * t530) * qJD(5)) * t531;
t994 = t17 + t19;
t18 = (qJD(3) * t657 - t101) * t533 + (qJD(3) * t207 + t105 * t532 + t530 * t97 + (t201 * t532 - t213 * t530) * qJD(5)) * t531;
t20 = (qJD(3) * t655 - t99) * t533 + (qJD(3) * t204 - t103 * t530 + t107 * t532 + (-t210 * t532 - t216 * t530) * qJD(5)) * t531;
t993 = t18 + t20;
t990 = t1048 * t453 - t1049 * t454 + t508 * t979;
t978 = t1041 * t554;
t977 = t1041 * t555;
t976 = t1042 * t554;
t975 = t1042 * t555;
t836 = t597 - t465;
t537 = t555 * qJ(2);
t493 = pkin(1) * t554 - t537;
t812 = -t554 * t523 - t526;
t371 = t493 + t812;
t367 = qJD(1) * t371;
t475 = qJD(1) * t493;
t972 = t367 - t475;
t541 = t554 * rSges(4,3);
t411 = rSges(4,1) * t856 - rSges(4,2) * t858 + t541;
t971 = t411 + t719;
t325 = -t747 - t984;
t400 = qJD(3) * t467 - t796;
t840 = -t325 * qJD(3) - t400;
t919 = rSges(7,2) * t531;
t353 = t533 * t686 + t919;
t402 = (-rSges(7,1) * t530 + rSges(7,3) * t532) * t531;
t799 = qJD(3) * t533;
t605 = t530 * t799 + t753;
t848 = t473 + t605 * qJ(6) + (-t530 * t793 + t532 * t799) * pkin(5) + qJD(3) * t353 + qJD(5) * t402;
t970 = t473 + t840 - t848;
t536 = t554 * qJ(2);
t495 = t555 * pkin(1) + t536;
t920 = rSges(3,2) * sin(pkin(9));
t922 = rSges(3,1) * t551;
t643 = t554 * rSges(3,3) + (-t920 + t922) * t555;
t969 = t495 + t643;
t816 = t461 + t666;
t817 = -t459 + t462;
t968 = t1005 * t531 + t361 * t803 + (-t531 * t816 + t533 * t817) * qJD(1);
t967 = t1051 * t531;
t596 = t375 * t555 - t376 * t554;
t943 = t554 * (-t459 * t555 + t378) - t555 * (-Icges(4,2) * t857 + t377 - t500);
t966 = -t531 * t943 + (-t244 * t555 + t246 * t554 + t596) * t533;
t965 = (-t1042 + t1075) * t508 + (-t1046 * t416 + t1056 + t394 - t397) * t454 + (t1046 * t418 - t1055 + t399 - t884) * t453;
t964 = (-Icges(7,1) * t861 + t1054 + t384 + t474) * t508 + (t1047 * t415 + t1012 - t395 + t398) * t454 + (-t1047 * t417 + t1057 + t396 - t891) * t453;
t963 = t1040 * t531;
t962 = -t1011 * t454 + t1043 * t508 + t1044 * t453;
t288 = t597 * t554;
t289 = t597 * t555;
t434 = t465 * t554;
t437 = t465 * t555;
t519 = qJD(4) * t531;
t768 = -t434 * t798 - t437 * t797 + t519;
t776 = t555 * t238 + t554 * t239 + t436 * t801;
t961 = t555 * t131 + t554 * t132 - t242 * t801 - t288 * t798 - t289 * t797 - t768 + t776;
t960 = t1048 * t555 + t1049 * t554;
t959 = t1048 * t554 - t1049 * t555;
t958 = t1018 * t555 + t1019 * t554;
t957 = t1018 * t554 - t1019 * t555;
t956 = t1020 * t555 + t1021 * t554;
t955 = t1020 * t554 - t1021 * t555;
t950 = t1050 * t508 + t742 * t979;
t691 = rSges(5,1) * t550 - rSges(5,2) * t548;
t368 = -rSges(5,3) * t533 + t531 * t691;
t613 = qJD(3) * t457;
t947 = -t555 * t613 + (-t458 * t554 + t646 + t875) * qJD(1);
t805 = qJD(1) * t374;
t946 = qJD(1) * t647 - t554 * t613 + t805;
t942 = m(5) / 0.2e1;
t941 = m(6) / 0.2e1;
t940 = m(7) / 0.2e1;
t939 = t357 / 0.2e1;
t938 = t358 / 0.2e1;
t937 = -t453 / 0.2e1;
t936 = t453 / 0.2e1;
t935 = -t454 / 0.2e1;
t934 = t454 / 0.2e1;
t933 = -t508 / 0.2e1;
t932 = t508 / 0.2e1;
t930 = t554 / 0.2e1;
t929 = -t555 / 0.2e1;
t927 = pkin(4) * t548;
t924 = pkin(1) - t523;
t921 = rSges(4,1) * t533;
t917 = rSges(5,3) * t531;
t915 = rSges(6,3) * t531;
t913 = t17 * t454;
t912 = t18 * t453;
t911 = t19 * t454;
t910 = t20 * t453;
t779 = t181 * rSges(6,1) + t180 * rSges(6,2) + rSges(6,3) * t757;
t110 = -rSges(6,3) * t763 + t779;
t354 = t533 * t688 + t915;
t403 = (-rSges(6,1) * t530 - rSges(6,2) * t532) * t531;
t193 = qJD(3) * t354 + qJD(5) * t403;
t225 = t418 * rSges(6,1) - t417 * rSges(6,2) + rSges(6,3) * t858;
t725 = t836 * t555;
t589 = qJD(1) * t725 + t554 * t840;
t743 = t533 * t785;
t527 = qJ(2) * t801;
t787 = qJD(1) * qJD(2);
t534 = qJD(2) * t554;
t809 = t527 + t534;
t818 = qJD(1) * (-pkin(1) * t802 + t809) + t554 * t787;
t771 = qJD(1) * (-t527 + (t554 * t924 - t526) * qJD(1)) + t818;
t640 = t554 * t743 + t771 + (t238 + t490) * qJD(1);
t609 = qJD(1) * t131 + t640;
t21 = t110 * t508 - t193 * t453 - t352 * t358 + (t225 * t793 + t589) * qJD(3) + t609;
t909 = t21 * t555;
t690 = rSges(6,1) * t183 - rSges(6,2) * t182;
t112 = rSges(6,3) * t423 + t690;
t745 = t554 * t786;
t525 = t555 * t787;
t767 = t465 * t745 + t555 * t743 + t525;
t450 = qJD(1) * t495 - t535;
t515 = t553 * t802;
t831 = t515 - (-t555 * t924 - t536) * qJD(1) - t450;
t775 = -t239 + t831;
t566 = -t597 * t745 + (-t132 - t488 + t775) * qJD(1) + t767;
t726 = t840 * t555;
t22 = -t112 * t508 - t193 * t454 + t352 * t357 + (-t221 * t793 + t726) * qJD(3) + t566;
t908 = t22 * t554;
t902 = t74 * t357;
t901 = t75 * t358;
t900 = t76 * t357;
t899 = t77 * t358;
t898 = -rSges(5,3) - qJ(4);
t410 = rSges(4,1) * t857 - rSges(4,2) * t859 - t555 * rSges(4,3);
t466 = rSges(4,1) * t531 + rSges(4,2) * t533;
t754 = t466 * t797;
t697 = t534 - t754;
t826 = t371 - t493;
t147 = (-t410 + t826) * qJD(1) + t697;
t869 = t147 * t554;
t761 = t466 * t798;
t148 = qJD(1) * t971 - t535 - t761;
t438 = t466 * t555;
t868 = t148 * t438;
t261 = t446 * rSges(5,1) + t445 * rSges(5,2) + rSges(5,3) * t858;
t843 = -t261 - t439;
t842 = t1017 * t416 - t1030 * t415;
t841 = t1017 * t418 - t1030 * t417;
t442 = t465 * t802;
t839 = -t597 * t802 + t442;
t838 = t830 * t554;
t837 = t830 * t555;
t835 = -t325 - t467;
t369 = t533 * t691 + t917;
t832 = -t369 * qJD(3) - t400;
t829 = -t685 * t533 - t353;
t828 = -t368 - t465;
t827 = -t369 - t467;
t821 = -(-pkin(5) * t530 + qJ(6) * t532) * t531 - t402;
t820 = t554 * t436 + t555 * t439;
t819 = -qJD(1) * t437 + t533 * t795;
t814 = rSges(4,2) * t763 + rSges(4,3) * t801;
t813 = t490 + t534;
t516 = t554 * t920;
t811 = rSges(3,3) * t801 + qJD(1) * t516;
t810 = t555 * rSges(3,3) + t516;
t808 = t534 - t475;
t804 = qJD(1) * t458;
t800 = qJD(3) * t531;
t782 = t554 * t922;
t777 = -t193 + t840;
t774 = qJD(1) * t289 + t819;
t773 = t327 * rSges(5,1) + t326 * rSges(5,2) + rSges(5,3) * t757;
t772 = -t352 + t836;
t770 = -t436 + t826;
t764 = t515 + t807;
t748 = -pkin(1) - t922;
t739 = t801 / 0.2e1;
t738 = t800 / 0.2e1;
t737 = -t798 / 0.2e1;
t736 = t798 / 0.2e1;
t734 = t797 / 0.2e1;
t608 = qJD(3) * t725 + t813;
t571 = (t242 + t770) * qJD(1) + t608;
t47 = -t454 * t830 - t1031 + t571;
t730 = t47 * t830;
t87 = -t368 * t798 + (t261 + t769) * qJD(1) + t765;
t728 = t87 * t828;
t724 = t835 * t554;
t723 = t835 * t555;
t722 = t828 * t555;
t720 = -t373 + t864;
t718 = -t522 * t533 - t523;
t717 = qJD(1) * t434 + t533 * t794;
t712 = -t554 * t242 + t555 * t243 + t820;
t710 = -t830 + t836;
t699 = qJD(5) * t738;
t694 = -rSges(4,2) * t531 + t921;
t693 = rSges(5,1) * t329 + rSges(5,2) * t328;
t692 = rSges(5,1) * t444 - rSges(5,2) * t443;
t419 = qJD(1) * t436;
t671 = t367 + t490 - t419 + t808;
t659 = -t147 * t555 - t148 * t554;
t654 = t221 * t555 - t225 * t554;
t639 = t764 - t766;
t638 = t718 - t919;
t637 = t718 - t915;
t623 = t702 + t812;
t435 = t466 * t554;
t621 = -qJD(1) * t288 + t717;
t226 = (t410 * t554 + t411 * t555) * qJD(3);
t607 = -t467 - t917;
t604 = t42 * t896 + t5 * t847;
t599 = -t47 * t847 + t48 * t845;
t598 = -t42 * t845 + t730;
t259 = rSges(5,3) * t859 + t692;
t595 = t642 + t719;
t572 = -t522 * t759 - t552 * t757 + t813 + t815;
t59 = t221 * t453 + t225 * t454 + t616;
t62 = t1006 + t571;
t63 = t225 * t508 - t352 * t453 + t578;
t565 = t59 * t654 + (t554 * t62 - t555 * t63) * t352;
t564 = qJD(1) * t373 - qJD(3) * t194 - t256 * t531 + t258 * t533;
t563 = -qJD(3) * t195 - t255 * t531 + t257 * t533 + t805;
t561 = -t1008 * t555 + t598 * t554;
t451 = t694 * qJD(3);
t427 = t782 - t810;
t422 = (t554 ^ 2 + t555 ^ 2) * t800;
t323 = t368 * t555;
t322 = t368 * t554;
t321 = t365 * t555;
t320 = t365 * t554;
t319 = t363 * t555;
t318 = t363 * t554;
t313 = t352 * t555;
t311 = t352 * t554;
t307 = qJD(1) * t969 - t535;
t306 = t534 + (-t427 - t493) * qJD(1);
t285 = -rSges(6,1) * t417 - rSges(6,2) * t418;
t280 = -rSges(6,1) * t415 - rSges(6,2) * t416;
t265 = -qJD(3) * t435 + (t555 * t694 + t541) * qJD(1);
t264 = rSges(4,1) * t606 - rSges(4,2) * t757 + t814;
t241 = t525 + (-qJD(1) * t643 - t450) * qJD(1);
t240 = qJD(1) * (-qJD(1) * t782 + t811) + t818;
t237 = qJD(1) * t242;
t146 = rSges(5,3) * t423 + t693;
t145 = -rSges(5,3) * t763 + t773;
t114 = -t451 * t797 + t525 + (-t265 + t761 + t831) * qJD(1);
t113 = -t451 * t798 + (t264 - t754) * qJD(1) + t771;
t92 = (t259 * t554 + t261 * t555) * qJD(3) + t696;
t86 = qJD(3) * t722 + (-t259 + t770) * qJD(1) + t813;
t52 = t832 * t797 + (-t146 + (qJD(3) * t368 - t519) * t554 + t775) * qJD(1) + t767;
t51 = qJD(1) * t145 + (qJD(1) * t722 + t554 * t832) * qJD(3) + t640;
t39 = (t145 * t555 + t146 * t554 + (t259 * t555 + t554 * t843) * qJD(1)) * qJD(3) + t714;
t8 = qJD(3) * t698 + t110 * t454 + t112 * t453 + t221 * t358 - t225 * t357 + t610;
t7 = qJD(6) * t182 + t897 * t508 - t848 * t453 - t830 * t358 + (t793 * t845 + t589) * qJD(3) + t609;
t6 = -qJD(6) * t180 - t896 * t508 - t848 * t454 + t830 * t357 + (-t793 * t847 + t726) * qJD(3) + t566;
t1 = [(((t555 * t720 + t144 + t653 - t833) * t555 + (t554 * t720 - t1004 + t1016 + t721 - t79) * t554) * qJD(3) + t1029 - t1034) * t737 + (-(-qJD(1) * t427 - t306 + t808) * t307 + t241 * (t554 * t748 + t537 + t810) + t306 * t535 + t240 * t969 + t307 * (t809 + t811) + (t306 * (t748 + t920) * t555 + (t306 * (-rSges(3,3) - qJ(2)) + t307 * t748) * t554) * qJD(1)) * m(3) + (((t142 - t331 + (t374 + t865) * t555 + t834) * t555 - t903 + (t81 + t833) * t554) * qJD(3) + t1045) * t734 + (t1011 * t533 + (t1012 * t530 - t1056 * t532) * t531 + t1049) * t453 * t933 + t991 * t934 + (t114 * (-t410 + t812) + t147 * (t515 + t535) + t113 * t971 + t148 * (t534 + t814) + (t466 * t869 - t868) * qJD(3) + ((-t147 * rSges(4,3) + t148 * (-t523 - t921)) * t554 + (t147 * (-t523 - t694) - t148 * t553) * t555) * qJD(1) - (-qJD(1) * t410 - t147 + t697 + t972) * t148) * m(4) + t901 / 0.2e1 + t902 / 0.2e1 + (t22 * ((-rSges(6,3) + t552) * t859 + t623 + t689) + t62 * (-rSges(6,3) * t758 + t639 - t690) + t21 * (t595 + t225) + t63 * (t572 + t779) + ((-t63 * t553 + t62 * t637) * t555 + (-t62 * t927 + t63 * t637) * t554) * qJD(1) - (t1006 + t237 - t419 + t608 - t62 + t972) * t63) * m(6) - t913 / 0.2e1 + t899 / 0.2e1 + t900 / 0.2e1 + t912 / 0.2e1 + t910 / 0.2e1 - t911 / 0.2e1 + ((t1014 + t1061) * t555 + (t1015 + t1062) * t554) * t786 / 0.2e1 + t980 * t938 + t981 * t939 + ((t448 - t334) * t533 + (-t335 * t548 + t336 * t550 + t449) * t531 + (t361 * t531 + t533 * t648 - t645) * qJD(3)) * qJD(1) + t950 + (t6 * ((-rSges(7,2) + t552) * t859 + t623 + t1052) + t7 * (t595 + t845) + t729 * t454 + (-rSges(7,2) * t758 + t639 - t989 + (-t554 * t927 + t555 * t638) * qJD(1)) * t47 + (t572 + t986 + (t554 * t638 - t526) * qJD(1) - t797 * t836 - t237 - t671 + t47 + t1031) * t48) * m(7) + (t1023 + t1024) * t736 + t1027 * t936 - (t1022 - t1025 + t1028) * t797 / 0.2e1 + (t1026 + t991) * t935 + (-t728 * t797 + t52 * (-t692 + t812) + t51 * (t502 - t843) + (-t51 * t553 + t52 * t607) * t554 + (t484 - t693 + t764 + t799 * t898 * t554 + (-t523 + t607) * t801) * t86 + (-pkin(3) * t759 + t479 - t671 + t773 + t813 + t86 + (t259 + (t531 * t898 - t523 - t928) * t554 - t526) * qJD(1)) * t87) * m(5); 0.2e1 * (t6 * t930 + t7 * t929) * m(7) + 0.2e1 * (-t909 / 0.2e1 + t908 / 0.2e1) * m(6) + 0.2e1 * (t51 * t929 + t52 * t930) * m(5) + 0.2e1 * (t113 * t929 + t114 * t930) * m(4) + 0.2e1 * (t240 * t929 + t241 * t930) * m(3); t955 * t939 + t957 * t938 + ((t1019 * t857 + t980 * t531) * qJD(5) + ((qJD(5) * t1018 + t962) * t533 + t967) * t555 + (t417 * t974 + t418 * t973) * t508 + (-t417 * t978 + t418 * t976) * t454 + (t417 * t977 - t418 * t975) * t453) * t937 + (qJD(1) * t958 + t554 * t997 - t555 * t999) * t936 + (qJD(1) * t956 + t554 * t995 - t555 * t996) * t935 + ((t1020 * t856 + t981 * t531) * qJD(5) + ((qJD(5) * t1021 + t962) * t533 + t967) * t554 + (t415 * t974 + t416 * t973) * t508 + (-t415 * t978 + t416 * t976) * t454 + (t415 * t977 - t416 * t975) * t453) * t934 + ((t960 * qJD(5) - t1051) * t533 + ((t530 * t974 + t532 * t973 + t1043) * t508 + (-t530 * t978 + t532 * t976 - t1011) * t454 + (t530 * t977 - t532 * t975 + t1044) * t453 + t979 * qJD(5)) * t531) * t933 + (qJD(1) * t960 + t554 * t993 - t555 * t994) * t932 - (-t1005 * t533 + ((-t364 * t548 + t366 * t550 + t361) * qJD(1) + ((t319 * t548 - t321 * t550 + t246) * t554 - (t318 * t548 - t320 * t550 + t244) * t555) * qJD(3)) * t531 + (t531 * t817 + t533 * t816) * qJD(1) + (t596 * t531 + t533 * t943) * qJD(3)) * qJD(1) / 0.2e1 + (t1025 * t555 + t1024 * t554 + (t1014 * t555 + t1015 * t554) * qJD(1)) * qJD(1) / 0.2e1 + ((-t319 * t445 - t321 * t446) * t798 + (t364 * t445 + t366 * t446) * qJD(1) + (-t798 * t862 + t804) * t554 + ((t318 * t445 + t320 * t446 + t554 * t863 + t966) * qJD(3) + t968) * t555) * t737 + (-(t318 * t443 - t320 * t444) * t797 + (-t364 * t443 + t366 * t444) * qJD(1) + (-t797 * t863 - t804) * t555 + ((t319 * t443 - t321 * t444 + t555 * t862 + t966) * qJD(3) + t968) * t554) * t734 - t990 * t793 / 0.2e1 + t959 * t699 + (-(t531 * t599 + t533 * t561) * qJD(5) + t5 * t712 + (qJD(1) * t731 + t6 * t710 + t998) * t555 + (qJD(1) * t730 + t7 * t710 + t604) * t554 + (-t724 * qJD(3) - t829 * t453 + t837 * t508 + t554 * t970 + t710 * t801 - t774) * t48 + (-t752 + t837 * t454 + t838 * t453 + t897 * t555 + (t844 - t845) * t802 + t961) * t42 + (-t723 * qJD(3) - t829 * t454 - t838 * t508 + t555 * t970 - t621 + t839) * t47) * m(7) + (t62 * t839 + t8 * t712 + (t8 * t225 + t62 * t777 + (qJD(1) * t63 + t22) * t772) * t555 + (qJD(1) * t352 * t62 + t21 * t772 + t8 * t221 + t63 * t777) * t554 - t62 * (t311 * t508 - t354 * t454 + t621) - t63 * (-t313 * t508 - t354 * t453 + t774) - (t62 * t723 + t63 * t724) * qJD(3) - ((-t221 * t62 + t225 * t63) * t531 + t565 * t533) * qJD(5) + ((qJD(1) * t221 + t110) * t555 + (t112 + (-t225 + t844) * qJD(1)) * t554 + t311 * t453 + t313 * t454 + t961) * t59) * m(6) + (t86 * t442 + t39 * t820 + t92 * t776 + (t52 * t828 + t86 * t832 + t39 * t261 + t92 * t145 + (t92 * t259 + t728) * qJD(1)) * t555 + (t51 * t828 + t87 * t832 + t39 * t259 + t92 * t146 + (t86 * t368 + t843 * t92) * qJD(1)) * t554 - t86 * (qJD(1) * t322 + t717) - t87 * (-qJD(1) * t323 + t819) - t92 * t768 - ((-t92 * t323 + t827 * t86) * t555 + (-t92 * t322 + t827 * t87) * t554) * qJD(3)) * m(5) + (-(t147 * t435 - t868) * qJD(1) - (t226 * (-t435 * t554 - t438 * t555) + t659 * t694) * qJD(3) + 0.2e1 * t226 * (t264 * t555 + t265 * t554 + (t410 * t555 - t411 * t554) * qJD(1)) + t659 * t451 + (-t113 * t554 - t114 * t555 + (-t148 * t555 + t869) * qJD(1)) * t466) * m(4) + (t1023 * qJD(1) + t1001 + ((-t134 * t858 - t136 * t445 - t138 * t446 - t244 * t424 - t247 * t326 - t250 * t327 - t564 * t555) * t555 + (t133 * t858 + t135 * t445 + t137 * t446 + t246 * t424 + t249 * t326 + t252 * t327 + t947 * t554 + (t563 - t946) * t555) * t554 + (t1016 * t554 + t1068 * t555) * qJD(1)) * t1002) * t930 + (t1022 * qJD(1) + t1000 + ((-t134 * t859 + t136 * t443 - t138 * t444 - t244 * t423 - t247 * t328 - t250 * t329 + t555 * t946) * t555 + (t133 * t859 - t135 * t443 + t137 * t444 + t246 * t423 + t249 * t328 + t252 * t329 + t563 * t554 + (-t564 - t947) * t555) * t554 + (t1069 * t555 + t1070 * t554) * qJD(1)) * t1002) * t929 + (t992 + t1029 + t1038) * t802 / 0.2e1 + (t991 + t1028 + t1039) * t739 - (t554 * t992 + t555 * t991) * t792 / 0.2e1; -m(5) * (t422 * t92 + t423 * t87 + t424 * t86) - m(6) * (t422 * t59 + t423 * t63 + t424 * t62) - m(7) * (t42 * t422 + t423 * t48 + t424 * t47) + 0.2e1 * ((t797 * t86 + t798 * t87 - t39) * t942 + (t62 * t797 + t63 * t798 - t8) * t941 + (t47 * t797 + t48 * t798 - t5) * t940) * t533 + 0.2e1 * ((qJD(3) * t92 + t51 * t554 + t52 * t555 + t801 * t87 - t802 * t86) * t942 + (qJD(3) * t59 + t21 * t554 + t22 * t555 - t62 * t802 + t63 * t801) * t941 + (qJD(3) * t42 - t47 * t802 + t48 * t801 + t554 * t7 + t555 * t6) * t940) * t531; (t531 * t956 - t533 * t981) * t939 + (t531 * t958 - t533 * t980) * t938 + (t417 * t965 + t418 * t964 - t555 * t963) * t937 + ((qJD(3) * t958 - t1027) * t533 + (-qJD(1) * t957 + qJD(3) * t980 + t554 * t999 + t555 * t997) * t531) * t936 + ((qJD(3) * t956 - t1026) * t533 + (-qJD(1) * t955 + qJD(3) * t981 + t554 * t996 + t555 * t995) * t531) * t935 + (t415 * t965 + t416 * t964 - t554 * t963) * t934 + (t1040 * t533 + (t530 * t965 + t532 * t964) * t531) * t933 + ((qJD(3) * t960 - t1050) * t533 + (-qJD(1) * t959 + qJD(3) * t979 + t554 * t994 + t555 * t993) * t531) * t932 - (t901 + t902 + t912 - t913 + t899 + t900 + t910 - t911 + t950) * t533 / 0.2e1 + t1000 * t859 / 0.2e1 + t1001 * t858 / 0.2e1 + t990 * t738 + (t531 * t960 - t533 * t979) * t699 + (-(t416 * t48 + t418 * t47 + t42 * t860) * qJD(6) - (-t47 * t842 + t48 * t841) * t508 - (t42 * t841 + t47 * t821) * t454 - (t42 * t842 + t48 * t821) * t453 + (qJD(3) * t561 + t47 * t896 - t48 * t897 + t6 * t847 - t7 * t845) * t533 + (t599 * qJD(3) + (qJD(1) * t598 - t48 * t848 - t7 * t830 + t604) * t555 + (qJD(1) * t1008 - t42 * t897 + t47 * t848 + t6 * t830 - t998) * t554) * t531) * m(7) + ((qJD(3) * t565 - t63 * t110 + t62 * t112 - t21 * t225 + t22 * t221) * t533 + (t62 * (-qJD(3) * t221 + t193 * t554) + t63 * (qJD(3) * t225 - t193 * t555) + t8 * t654 + t59 * (-t110 * t554 + t112 * t555 - t221 * t802 - t225 * t801) + (-t909 + t908 + (t554 * t63 + t555 * t62) * qJD(1)) * t352) * t531 - t62 * (-t280 * t508 - t403 * t454) - t63 * (t285 * t508 - t403 * t453) - t59 * (t280 * t453 + t285 * t454)) * m(6) + t991 * (t533 * t734 - t763 / 0.2e1) + t992 * (t531 * t739 + t533 * t736); (t415 * t7 + t417 * t6 + t5 * t861 + (-t417 * t508 + t453 * t861 + t182) * t48 + (t415 * t508 + t454 * t861 - t180) * t47 + (-t415 * t453 - t417 * t454 + t605) * t42) * m(7);];
tauc  = t1(:);