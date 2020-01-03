% Calculate vector of inverse dynamics joint torques for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPR8_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR8_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR8_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR8_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR8_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:19:08
% EndTime: 2019-12-31 21:20:25
% DurationCPUTime: 66.05s
% Computational Cost: add. (28386->1257), mult. (39127->1595), div. (0->0), fcn. (34758->8), ass. (0->625)
t1098 = Icges(5,4) - Icges(4,5);
t1097 = Icges(5,5) - Icges(4,6);
t1096 = Icges(5,1) + Icges(4,3);
t575 = qJ(2) + qJ(3);
t554 = cos(t575);
t532 = Icges(4,4) * t554;
t553 = sin(t575);
t439 = Icges(4,1) * t553 + t532;
t913 = Icges(5,6) * t554;
t683 = Icges(5,2) * t553 + t913;
t1090 = t683 + t439;
t578 = sin(qJ(1));
t882 = t554 * t578;
t885 = t553 * t578;
t581 = cos(qJ(1));
t915 = Icges(4,6) * t581;
t329 = Icges(4,4) * t882 - Icges(4,2) * t885 - t915;
t920 = Icges(5,5) * t581;
t334 = Icges(5,6) * t882 - Icges(5,3) * t885 + t920;
t1095 = t329 + t334;
t691 = -Icges(4,2) * t553 + t532;
t330 = Icges(4,6) * t578 + t581 * t691;
t682 = -Icges(5,3) * t553 + t913;
t333 = Icges(5,5) * t578 - t581 * t682;
t1089 = t330 - t333;
t503 = Icges(4,4) * t885;
t921 = Icges(4,5) * t581;
t331 = Icges(4,1) * t882 - t503 - t921;
t498 = Icges(5,6) * t885;
t926 = Icges(5,4) * t581;
t336 = Icges(5,2) * t882 - t498 + t926;
t1077 = t331 + t336;
t928 = Icges(4,4) * t553;
t440 = Icges(4,1) * t554 - t928;
t332 = Icges(4,5) * t578 + t440 * t581;
t884 = t553 * t581;
t499 = Icges(5,6) * t884;
t880 = t554 * t581;
t927 = Icges(5,4) * t578;
t335 = -Icges(5,2) * t880 + t499 + t927;
t1094 = t332 - t335;
t1074 = t1097 * t553 - t1098 * t554;
t914 = Icges(5,6) * t553;
t684 = Icges(5,2) * t554 - t914;
t1093 = -t440 - t684;
t437 = Icges(4,2) * t554 + t928;
t681 = Icges(5,3) * t554 + t914;
t1092 = t681 + t437;
t1091 = t682 + t691;
t433 = Icges(4,5) * t553 + Icges(4,6) * t554;
t689 = Icges(5,4) * t553 + Icges(5,5) * t554;
t1088 = t433 - t689;
t1087 = t1096 * t581;
t1086 = t437 * t581 - t1094;
t1085 = -t1090 - t1091;
t1084 = -t1092 - t1093;
t1083 = t1090 * t578 + t1095;
t1082 = -t1090 * t581 - t1089;
t1078 = -t1097 * t885 + t1098 * t882 + t1087;
t1081 = t1074 * t581 + t1096 * t578;
t663 = t437 * t553 - t439 * t554;
t1080 = -t553 * t681 + t554 * t683 - t663;
t1073 = -t1077 * t554 + t1095 * t553;
t1079 = -t332 * t882 - t333 * t885;
t1076 = t1088 * t578;
t387 = t689 * t581;
t890 = t433 * t581;
t1075 = t387 - t890;
t572 = qJD(2) + qJD(3);
t1072 = t1084 * t572;
t1071 = t1085 * t572;
t1070 = t1082 * t572 + (t1093 * t578 + t921 - t926) * qJD(1);
t1069 = t1083 * t572 + (-t581 * t684 - t332 + t927) * qJD(1);
t641 = t572 * t681;
t1068 = t581 * t641 + t1086 * t572 + (t1091 * t578 - t915 + t920) * qJD(1);
t879 = t572 * t578;
t1067 = t1089 * qJD(1) + t1077 * t572 - t437 * t879 - t578 * t641;
t1066 = t330 * t553 + t335 * t554;
t151 = -t578 * t663 - t890;
t154 = t681 * t885 - t683 * t882 - t387;
t1065 = t151 - t154;
t1029 = t1080 * t581 + t1076;
t1064 = t1073 * t578;
t1063 = t1081 * t581 + t1079;
t1062 = t1081 * qJD(1);
t1008 = -t1081 * t578 - t332 * t880 - t333 * t884;
t1061 = t1078 * t578 - t331 * t880 + t334 * t884;
t1060 = t1078 * t581;
t579 = cos(qJ(5));
t874 = t578 * t579;
t576 = sin(qJ(5));
t877 = t576 * t581;
t412 = t553 * t874 + t877;
t872 = t579 * t581;
t878 = t576 * t578;
t413 = -t553 * t878 + t872;
t836 = t413 * rSges(6,1) - t412 * rSges(6,2);
t219 = rSges(6,3) * t882 - t836;
t949 = rSges(6,1) * t576;
t709 = rSges(6,2) * t579 + t949;
t318 = rSges(6,3) * t553 - t554 * t709;
t473 = t572 * t581;
t801 = qJD(5) * t578;
t371 = -t554 * t801 + t473;
t804 = qJD(5) * t553;
t517 = qJD(1) + t804;
t1059 = t219 * t517 + t371 * t318;
t1033 = t1060 - t1064;
t1032 = -t330 * t885 - t335 * t882 - t1063;
t1031 = -t329 * t884 + t336 * t880 - t1061;
t1030 = -t330 * t884 - t335 * t880 - t1008;
t1058 = t1088 * qJD(1) + t1071 * t553 + t1072 * t554;
t1057 = t1078 * qJD(1) + t1067 * t553 + t1069 * t554;
t1056 = t1068 * t553 + t1070 * t554 + t1062;
t552 = qJD(2) * t578;
t472 = qJD(3) * t578 + t552;
t1055 = t1084 * qJD(1) + t1082 * t472 + t1083 * t473;
t1054 = (-t498 - t503 + (-Icges(4,2) - Icges(5,3)) * t882 + t1077) * t473 + (Icges(5,3) * t880 + t1086 + t499) * t472 + t1085 * qJD(1);
t1053 = t1080 * qJD(1) - t1074 * t572;
t1052 = -t1075 * t572 + (t1074 * t578 + t332 * t554 + t333 * t553 - t1066 - t1087) * qJD(1);
t1051 = t1073 * qJD(1) - t1076 * t572 + t1062;
t1050 = t1029 * qJD(1);
t1049 = t1065 * qJD(1);
t410 = t553 * t872 - t878;
t411 = t553 * t877 + t874;
t685 = Icges(6,5) * t576 + Icges(6,6) * t579;
t613 = -Icges(6,3) * t553 + t554 * t685;
t924 = Icges(6,4) * t576;
t687 = Icges(6,2) * t579 + t924;
t614 = -Icges(6,6) * t553 + t554 * t687;
t923 = Icges(6,4) * t579;
t693 = Icges(6,1) * t576 + t923;
t615 = -Icges(6,5) * t553 + t554 * t693;
t123 = -t410 * t614 - t411 * t615 - t613 * t880;
t800 = qJD(5) * t581;
t370 = t554 * t800 + t472;
t208 = Icges(6,5) * t411 + Icges(6,6) * t410 + Icges(6,3) * t880;
t925 = Icges(6,4) * t411;
t211 = Icges(6,2) * t410 + Icges(6,6) * t880 + t925;
t375 = Icges(6,4) * t410;
t214 = Icges(6,1) * t411 + Icges(6,5) * t880 + t375;
t91 = t208 * t880 + t410 * t211 + t411 * t214;
t210 = -Icges(6,5) * t413 + Icges(6,6) * t412 + Icges(6,3) * t882;
t377 = Icges(6,4) * t413;
t213 = Icges(6,2) * t412 + Icges(6,6) * t882 - t377;
t376 = Icges(6,4) * t412;
t215 = Icges(6,1) * t413 - Icges(6,5) * t882 - t376;
t92 = t210 * t880 + t410 * t213 - t215 * t411;
t33 = t123 * t517 + t370 * t91 - t371 * t92;
t124 = -t412 * t614 + t413 * t615 - t613 * t882;
t93 = t208 * t882 + t412 * t211 - t413 * t214;
t94 = t210 * t882 + t213 * t412 + t215 * t413;
t34 = t124 * t517 + t370 * t93 - t371 * t94;
t676 = t213 * t579 - t215 * t576;
t100 = t210 * t553 - t554 * t676;
t1047 = -t1051 * t578 + t1057 * t581;
t1046 = -t1052 * t578 + t1056 * t581;
t1045 = t1051 * t581 + t1057 * t578;
t1044 = t1052 * t581 + t1056 * t578;
t1043 = t1032 * t472 - t1033 * t473 + t1049;
t1042 = t1030 * t472 - t1031 * t473 + t1050;
t1041 = -t1053 * t578 + t1058 * t581;
t1040 = t1053 * t581 + t1058 * t578;
t1039 = -t1067 * t554 + t1069 * t553;
t1038 = -t1068 * t554 + t1070 * t553;
t577 = sin(qJ(2));
t1037 = rSges(3,2) * t577;
t1028 = t1077 * t553 + t1095 * t554;
t1027 = t1089 * t554 + t1094 * t553;
t217 = t411 * rSges(6,1) + t410 * rSges(6,2) + rSges(6,3) * t880;
t568 = t578 * pkin(4);
t427 = pkin(8) * t880 + t568;
t1026 = t217 + t427;
t428 = pkin(4) * t581 - pkin(8) * t882;
t1025 = t219 - t428;
t881 = t554 * t579;
t1009 = rSges(6,2) * t881 + t554 * t949;
t825 = t1009 * t578;
t266 = -rSges(6,3) * t885 + t825;
t824 = t1009 * t581;
t267 = -rSges(6,3) * t884 + t824;
t808 = qJD(1) * t581;
t768 = t553 * t808;
t315 = t472 * t554 + t768;
t317 = t554 * rSges(6,3) + t553 * t709;
t496 = qJ(4) * t880;
t403 = -pkin(3) * t884 + t496;
t531 = t553 * qJ(4);
t994 = t554 * pkin(3) + t531;
t401 = t994 * t578;
t570 = t581 * pkin(6);
t492 = pkin(1) * t578 - t570;
t580 = cos(qJ(2));
t569 = t580 * pkin(2);
t542 = t569 + pkin(1);
t582 = -pkin(7) - pkin(6);
t546 = t581 * t582;
t819 = -t578 * t542 - t546;
t325 = t492 + t819;
t567 = t578 * pkin(6);
t493 = t581 * pkin(1) + t567;
t505 = t581 * t542;
t728 = -t578 * t582 + t505;
t326 = t728 - t493;
t807 = qJD(2) * t581;
t857 = -t325 * t552 + t326 * t807;
t649 = -qJD(4) * t554 + t472 * t401 + t857;
t495 = qJ(4) * t884;
t406 = pkin(3) * t880 + t495;
t833 = t406 + t427;
t63 = t217 * t371 + t219 * t370 - t428 * t472 + t473 * t833 + t649;
t959 = -rSges(6,3) - pkin(3);
t792 = -pkin(8) + t959;
t725 = t581 * t792;
t726 = t578 * t792;
t759 = t553 * t801;
t847 = t326 + t493;
t773 = t406 + t847;
t806 = qJD(4) * t578;
t957 = pkin(3) * t553;
t441 = -qJ(4) * t554 + t957;
t942 = pkin(2) * qJD(2);
t786 = t577 * t942;
t522 = t578 * t786;
t849 = -t472 * t441 - t522;
t79 = t217 * t517 - t318 * t370 + (-pkin(8) * t472 + t806) * t553 + (t427 + t773) * qJD(1) + t849;
t956 = pkin(8) * t553;
t790 = t473 * t956;
t803 = qJD(5) * t554;
t489 = t554 * t806;
t840 = qJD(1) * t403 + t489;
t494 = qJ(4) * t882;
t398 = -pkin(3) * t885 + t494;
t530 = qJD(4) * t553;
t859 = t472 * t398 + t530;
t889 = t472 * t578;
t1017 = -t63 * (t267 * t371 + t403 * t473 + t217 * t759 + t370 * t266 - t581 * t790 + (-pkin(8) * t889 - t219 * t800) * t553 + t859) - t79 * (t553 * t318 * t800 - pkin(8) * t315 + t217 * t803 + t517 * t267 - t317 * t370 - t472 * t994 + t840) - (g(1) * t725 + g(2) * t726) * t553;
t1016 = t1054 * t553 + t1055 * t554;
t1015 = t1074 * qJD(1) + t1075 * t472 + t1076 * t473;
t1012 = t1066 + t1078;
t809 = qJD(1) * t578;
t767 = t554 * t809;
t783 = t553 * t473;
t1011 = t767 + t783;
t408 = t441 * t809;
t805 = qJD(4) * t581;
t491 = t554 * t805;
t820 = -t473 * t994 + t491;
t714 = -qJD(1) * t398 + t820;
t781 = t473 * t554;
t1010 = pkin(8) * t781 + t266 * t517 + t371 * t317 + t318 * t809 + t408 - t714;
t1007 = t33 * t581 + t34 * t578;
t766 = t554 * t808;
t784 = t553 * t879;
t629 = t766 - t784;
t723 = qJD(1) * t553 + qJD(5);
t656 = t723 * t581;
t657 = t517 * t576;
t183 = t579 * t656 + (t572 * t881 - t657) * t578;
t658 = t579 * t517;
t782 = t554 * t879;
t184 = t578 * t658 + (t656 + t782) * t576;
t712 = rSges(6,1) * t184 + rSges(6,2) * t183;
t115 = rSges(6,3) * t629 + t712;
t710 = -rSges(6,1) * t579 + rSges(6,2) * t576;
t886 = t553 * t572;
t943 = rSges(6,3) * t572;
t166 = t709 * t886 + (qJD(5) * t710 + t943) * t554;
t791 = -qJDD(2) - qJDD(3);
t793 = qJDD(5) * t554;
t802 = qJD(5) * t572;
t796 = qJD(1) * qJD(2);
t538 = t578 * t796;
t795 = qJD(1) * qJD(3);
t815 = t578 * t795 + t538;
t173 = (-t553 * t802 + t793) * t578 + (qJD(1) * t803 + t791) * t581 + t815;
t282 = pkin(4) * t809 + pkin(8) * t629;
t347 = t581 * t791 + t815;
t372 = qJDD(5) * t553 + t554 * t802 + qJDD(1);
t295 = qJ(4) * t886 + (pkin(3) * t572 - qJD(4)) * t554;
t461 = -qJDD(2) * t581 + t538;
t870 = t580 * qJD(2) ^ 2;
t958 = pkin(2) * t577;
t660 = -pkin(2) * t581 * t870 + t461 * t958;
t794 = qJDD(4) * t553;
t617 = -t473 * t295 + t347 * t441 + t572 * t491 + t581 * t794 + t660;
t458 = pkin(3) * t784;
t760 = t553 * t806;
t175 = pkin(3) * t766 + t760 - t458 + (t768 + t782) * qJ(4);
t817 = t582 * t809 + t522;
t951 = pkin(1) - t542;
t252 = (-t581 * t951 - t567) * qJD(1) - t817;
t459 = t493 * qJD(1);
t864 = -t252 - t459;
t635 = -t175 - t760 + t864;
t848 = t325 - t492;
t774 = -t401 + t848;
t719 = t428 + t774;
t883 = t554 * t572;
t24 = -t517 * t115 - t371 * t166 + t173 * t318 - t372 * t219 + (t347 * t553 - t473 * t883) * pkin(8) + t719 * qJDD(1) + (-t282 + t635) * qJD(1) + t617;
t1006 = (qJD(1) * t79 + t24) * t581;
t340 = rSges(4,1) * t882 - rSges(4,2) * t885 - t581 * rSges(4,3);
t562 = t578 * rSges(4,3);
t341 = rSges(4,1) * t880 - rSges(4,2) * t884 + t562;
t122 = t340 * t472 + t341 * t473 + t857;
t443 = rSges(4,1) * t553 + rSges(4,2) * t554;
t762 = t577 * t807;
t724 = pkin(2) * t762;
t631 = -t443 * t473 - t724;
t775 = -t340 + t848;
t127 = qJD(1) * t775 + t631;
t128 = -t443 * t472 - t522 + (t341 + t847) * qJD(1);
t400 = t443 * t578;
t405 = t443 * t581;
t535 = t554 * rSges(4,1);
t995 = -rSges(4,2) * t553 + t535;
t1005 = -t127 * (qJD(1) * t400 - t473 * t995) - t122 * (-t472 * t400 - t405 * t473) - t128 * (-qJD(1) * t405 - t472 * t995);
t616 = -t578 * t723 + t781;
t185 = t579 * t616 - t581 * t657;
t186 = t576 * t616 + t581 * t658;
t869 = t186 * rSges(6,1) + t185 * rSges(6,2);
t116 = -rSges(6,3) * t1011 + t869;
t460 = qJDD(2) * t578 + t581 * t796;
t346 = qJDD(3) * t578 + t581 * t795 + t460;
t172 = -qJD(5) * t1011 + t581 * t793 + t346;
t548 = pkin(4) * t808;
t283 = -pkin(8) * t1011 + t548;
t769 = t553 * t809;
t452 = t572 * t496;
t490 = t553 * t805;
t826 = t452 + t490;
t174 = -pkin(3) * t1011 - qJ(4) * t769 + t826;
t547 = pkin(6) * t808;
t251 = -t724 - t547 + (t578 * t951 - t546) * qJD(1);
t831 = qJD(1) * (-pkin(1) * t809 + t547) + qJDD(1) * t493;
t607 = qJD(1) * t251 + qJDD(1) * t326 + (-t460 * t577 - t578 * t870) * pkin(2) + t831;
t593 = qJDD(1) * t406 + t572 * t489 + t578 * t794 + t607 + (t174 + t490) * qJD(1);
t25 = t517 * t116 - t472 * t295 - t346 * t441 + qJDD(1) * t427 + t372 * t217 - t370 * t166 + t593 + (-t346 * t553 - t472 * t883) * pkin(8) + qJD(1) * t283 - t172 * t318;
t1003 = t25 * t578;
t720 = rSges(5,1) * t808 + rSges(5,2) * t1011 + rSges(5,3) * t781;
t204 = -rSges(5,3) * t769 + t720;
t564 = t578 * rSges(5,1);
t342 = -rSges(5,2) * t880 + rSges(5,3) * t884 + t564;
t946 = rSges(5,2) * t553;
t707 = rSges(5,3) * t554 + t946;
t828 = -t441 + t707;
t533 = t553 * rSges(5,3);
t945 = rSges(5,2) * t554;
t996 = t533 - t945;
t362 = t996 * t572;
t855 = -t295 - t362;
t40 = qJD(1) * t204 + qJDD(1) * t342 + t346 * t828 + t472 * t855 + t593;
t1002 = t40 * t578;
t102 = t760 + t707 * t472 + (t342 + t773) * qJD(1) + t849;
t203 = t707 * t879 + (t581 * t996 + t564) * qJD(1);
t748 = rSges(5,1) * t581 - rSges(5,3) * t885;
t343 = rSges(5,2) * t882 + t748;
t722 = t343 + t774;
t39 = -t347 * t707 - t362 * t473 + t722 * qJDD(1) + (-t203 + t635) * qJD(1) + t617;
t1000 = qJD(1) * t102 + t39;
t467 = qJD(1) * t492;
t997 = qJD(1) * t325 - t467;
t561 = Icges(3,4) * t580;
t692 = -Icges(3,2) * t577 + t561;
t480 = Icges(3,1) * t577 + t561;
t992 = g(1) * t581 + g(2) * t578;
t873 = t578 * t580;
t876 = t577 * t578;
t911 = Icges(3,3) * t581;
t364 = Icges(3,5) * t873 - Icges(3,6) * t876 - t911;
t527 = Icges(3,4) * t876;
t922 = Icges(3,5) * t581;
t368 = Icges(3,1) * t873 - t527 - t922;
t916 = Icges(3,6) * t581;
t366 = Icges(3,4) * t873 - Icges(3,2) * t876 - t916;
t894 = t366 * t577;
t667 = -t368 * t580 + t894;
t138 = -t364 * t581 - t578 * t667;
t477 = Icges(3,5) * t580 - Icges(3,6) * t577;
t476 = Icges(3,5) * t577 + Icges(3,6) * t580;
t632 = qJD(2) * t476;
t929 = Icges(3,4) * t577;
t481 = Icges(3,1) * t580 - t929;
t369 = Icges(3,5) * t578 + t481 * t581;
t367 = Icges(3,6) * t578 + t581 * t692;
t893 = t367 * t577;
t666 = -t369 * t580 + t893;
t987 = -t581 * t632 + (-t477 * t578 + t666 + t911) * qJD(1);
t365 = Icges(3,3) * t578 + t477 * t581;
t811 = qJD(1) * t365;
t986 = qJD(1) * t667 - t578 * t632 + t811;
t478 = Icges(3,2) * t580 + t929;
t661 = t478 * t577 - t480 * t580;
t983 = t661 * qJD(1) + t477 * qJD(2);
t982 = t578 * (-t478 * t581 + t369) - t581 * (-Icges(3,2) * t873 + t368 - t527);
t688 = Icges(6,2) * t576 - t923;
t605 = t370 * (-Icges(6,2) * t411 + t214 + t375) - t371 * (Icges(6,2) * t413 - t215 + t376) + t517 * (t688 * t554 - t615);
t694 = -Icges(6,1) * t579 + t924;
t606 = t370 * (-Icges(6,1) * t410 + t211 + t925) - t371 * (-Icges(6,1) * t412 + t213 - t377) + t517 * (-t694 * t554 - t614);
t979 = m(5) / 0.2e1;
t978 = m(6) / 0.2e1;
t977 = t172 / 0.2e1;
t976 = t173 / 0.2e1;
t975 = t346 / 0.2e1;
t974 = t347 / 0.2e1;
t973 = -t370 / 0.2e1;
t972 = t370 / 0.2e1;
t971 = -t371 / 0.2e1;
t970 = t371 / 0.2e1;
t969 = t460 / 0.2e1;
t968 = t461 / 0.2e1;
t967 = -t472 / 0.2e1;
t966 = t472 / 0.2e1;
t965 = -t473 / 0.2e1;
t964 = t473 / 0.2e1;
t963 = -t517 / 0.2e1;
t962 = t517 / 0.2e1;
t961 = t578 / 0.2e1;
t960 = -t581 / 0.2e1;
t953 = -qJD(1) / 0.2e1;
t952 = qJD(1) / 0.2e1;
t950 = rSges(3,1) * t580;
t948 = rSges(3,2) * t580;
t110 = Icges(6,5) * t186 + Icges(6,6) * t185 - Icges(6,3) * t1011;
t112 = Icges(6,4) * t186 + Icges(6,2) * t185 - Icges(6,6) * t1011;
t114 = Icges(6,1) * t186 + Icges(6,4) * t185 - Icges(6,5) * t1011;
t677 = t211 * t579 + t214 * t576;
t26 = (t572 * t677 + t110) * t553 + (-t112 * t579 - t114 * t576 + t208 * t572 + (t211 * t576 - t214 * t579) * qJD(5)) * t554;
t941 = t26 * t370;
t109 = Icges(6,5) * t184 + Icges(6,6) * t183 + Icges(6,3) * t629;
t111 = Icges(6,4) * t184 + Icges(6,2) * t183 + Icges(6,6) * t629;
t113 = Icges(6,1) * t184 + Icges(6,4) * t183 + Icges(6,5) * t629;
t27 = (t572 * t676 + t109) * t553 + (-t111 * t579 - t113 * t576 + t210 * t572 + (t213 * t576 + t215 * t579) * qJD(5)) * t554;
t940 = t27 * t371;
t563 = t578 * rSges(3,3);
t99 = t208 * t553 - t554 * t677;
t935 = t99 * t172;
t934 = qJDD(1) / 0.2e1;
t672 = -t576 * t615 - t579 * t614;
t126 = -t553 * t613 - t554 * t672;
t638 = t685 * t553;
t686 = -Icges(6,5) * t579 + Icges(6,6) * t576;
t162 = t572 * t638 + (Icges(6,3) * t572 + qJD(5) * t686) * t554;
t639 = t687 * t553;
t163 = t572 * t639 + (Icges(6,6) * t572 + qJD(5) * t688) * t554;
t640 = t693 * t553;
t164 = t572 * t640 + (Icges(6,5) * t572 + qJD(5) * t694) * t554;
t51 = (t572 * t672 + t162) * t553 + (-t163 * t579 - t164 * t576 - t613 * t572 + (-t576 * t614 + t579 * t615) * qJD(5)) * t554;
t933 = t126 * t372 + t51 * t517;
t908 = qJD(1) * t63;
t842 = t342 + t406;
t90 = -t343 * t472 + t473 * t842 + t649;
t906 = qJD(1) * t90;
t905 = t100 * t173;
t903 = t127 * t578;
t482 = rSges(3,1) * t577 + t948;
t763 = t482 * t807;
t816 = rSges(3,2) * t876 + t581 * rSges(3,3);
t373 = rSges(3,1) * t873 - t816;
t837 = -t373 - t492;
t221 = qJD(1) * t837 - t763;
t902 = t221 * t578;
t901 = t221 * t581;
t871 = t580 * t581;
t875 = t577 * t581;
t374 = rSges(3,1) * t871 - rSges(3,2) * t875 + t563;
t292 = t374 + t493;
t222 = qJD(1) * t292 - t482 * t552;
t425 = t482 * t581;
t900 = t222 * t425;
t888 = t476 * t578;
t887 = t476 * t581;
t856 = -t578 * t325 + t581 * t326;
t854 = t578 * t340 + t581 * t341;
t853 = -t578 * t364 - t368 * t871;
t852 = t578 * t365 + t369 * t871;
t841 = t578 * t401 + t581 * t406;
t399 = rSges(5,2) * t885 + rSges(5,3) * t882;
t835 = -t398 - t399;
t404 = rSges(5,2) * t884 + rSges(5,3) * t880;
t834 = t403 + t404;
t832 = -t707 * t809 + t408;
t827 = -t994 - t996;
t823 = rSges(4,2) * t769 + rSges(4,3) * t808;
t822 = -t478 + t481;
t821 = t480 + t692;
t818 = rSges(3,3) * t808 + t1037 * t809;
t810 = qJD(1) * t477;
t225 = -t578 * t661 - t887;
t797 = t225 * qJD(1);
t785 = t580 * t942;
t780 = t581 * t174 + t578 * t175 + t401 * t808;
t201 = -rSges(4,1) * t1011 - rSges(4,2) * t781 + t823;
t202 = -t572 * t400 + (t581 * t995 + t562) * qJD(1);
t779 = t581 * t201 + t578 * t202 + t340 * t808;
t778 = -t217 - t833;
t777 = t581 * t251 + t578 * t252 - t325 * t808;
t772 = t458 + t817;
t771 = t494 + t825;
t770 = t496 + t824;
t764 = t577 * t808;
t757 = -pkin(1) - t950;
t756 = t809 / 0.2e1;
t755 = t808 / 0.2e1;
t754 = -t552 / 0.2e1;
t753 = t552 / 0.2e1;
t752 = -t807 / 0.2e1;
t751 = t807 / 0.2e1;
t630 = -t443 - t958;
t749 = -t441 - t956;
t747 = (-t578 ^ 2 - t581 ^ 2) * t577;
t306 = t369 * t873;
t736 = t365 * t581 - t306;
t729 = -t364 + t893;
t721 = t581 * t342 - t578 * t343 + t841;
t718 = t828 - t958;
t717 = -t318 + t749;
t363 = t995 * t572;
t716 = -t363 - t785;
t715 = -t957 - t958;
t485 = rSges(2,1) * t581 - rSges(2,2) * t578;
t483 = rSges(2,1) * t578 + rSges(2,2) * t581;
t484 = t950 - t1037;
t224 = t367 * t580 + t369 * t577;
t633 = qJD(2) * t478;
t236 = -t581 * t633 + (-t578 * t692 + t916) * qJD(1);
t634 = qJD(2) * t480;
t238 = -t581 * t634 + (-t481 * t578 + t922) * qJD(1);
t595 = -qJD(2) * t224 - t236 * t577 + t238 * t580 + t811;
t223 = t366 * t580 + t368 * t577;
t237 = qJD(1) * t367 - t578 * t633;
t239 = qJD(1) * t369 - t578 * t634;
t596 = qJD(1) * t364 - qJD(2) * t223 - t237 * t577 + t239 * t580;
t705 = t578 * (t578 * t987 + t595 * t581) - t581 * (t578 * t986 + t596 * t581);
t704 = t578 * (t595 * t578 - t581 * t987) - t581 * (t596 * t578 - t581 * t986);
t655 = t490 - t724;
t643 = -t473 * t441 + t655;
t78 = qJD(1) * t719 - t1059 + t643 - t790;
t703 = t578 * t79 + t581 * t78;
t702 = t578 * t92 + t581 * t91;
t701 = t578 * t91 - t581 * t92;
t700 = t578 * t94 + t581 * t93;
t699 = t578 * t93 - t581 * t94;
t698 = t100 * t581 - t578 * t99;
t697 = t100 * t578 + t581 * t99;
t680 = -t127 * t581 - t128 * t578;
t139 = -t367 * t876 - t736;
t679 = -t138 * t581 + t139 * t578;
t140 = -t366 * t875 - t853;
t141 = -t367 * t875 + t852;
t678 = -t140 * t581 + t141 * t578;
t675 = t217 * t578 - t219 * t581;
t674 = -t222 * t578 - t901;
t240 = -t807 * t948 + (-t580 * t809 - t762) * rSges(3,1) + t818;
t424 = t482 * t578;
t241 = -qJD(2) * t424 + (t484 * t581 + t563) * qJD(1);
t673 = t240 * t581 + t241 * t578;
t665 = t373 * t578 + t374 * t581;
t662 = t478 * t580 + t480 * t577;
t659 = -pkin(8) * t883 - t166 - t295;
t654 = -t785 + t855;
t652 = t251 * t807 + t252 * t552 - t460 * t325 - t326 * t461;
t650 = t728 + t406;
t648 = t578 * t203 + t581 * t204 - t343 * t808 + t780;
t647 = t1025 * t578 + t1026 * t581 + t841;
t646 = t554 * pkin(8) + t317 + t994;
t642 = t717 - t958;
t627 = -t208 * t370 + t210 * t371 + t517 * t613;
t626 = (Icges(6,5) * t410 - Icges(6,6) * t411) * t370 - (Icges(6,5) * t412 + Icges(6,6) * t413) * t371 + t686 * t554 * t517;
t623 = t366 * t581 - t367 * t578;
t622 = -t542 - t994 - t533;
t621 = -qJD(1) * t401 + t655 + t997;
t620 = t659 - t785;
t619 = t780 + t1025 * t808 + (t116 + t283) * t581 + (t115 + t282) * t578;
t618 = t554 * t626;
t612 = (-t577 * t821 + t580 * t822) * qJD(1);
t609 = -qJDD(4) * t554 + t472 * t175 + t346 * t401 + t572 * t530 + t652;
t450 = t692 * qJD(2);
t451 = t481 * qJD(2);
t594 = qJD(1) * t476 - qJD(2) * t662 - t450 * t577 + t451 * t580;
t592 = t102 * (qJD(1) * t404 + t472 * t827 + t840) + t90 * (t472 * t399 + t473 * t834 + t859);
t591 = -t577 * t982 + t623 * t580;
t590 = (t613 * t581 + t677) * t370 - (t613 * t578 + t676) * t371 + (Icges(6,3) * t554 + t638 + t672) * t517;
t589 = t590 * t554;
t20 = t110 * t882 + t112 * t412 - t114 * t413 + t183 * t211 + t184 * t214 + t208 * t629;
t21 = t109 * t882 + t111 * t412 - t113 * t413 + t183 * t213 - t184 * t215 + t210 * t629;
t22 = -t1011 * t208 + t110 * t880 + t112 * t410 + t114 * t411 + t185 * t211 + t186 * t214;
t23 = -t1011 * t210 + t109 * t880 + t111 * t410 + t113 * t411 + t185 * t213 - t186 * t215;
t262 = t614 * t578;
t263 = t614 * t581;
t264 = t615 * t578;
t265 = t615 * t581;
t41 = t162 * t882 + t163 * t412 - t164 * t413 - t183 * t614 - t184 * t615 - t613 * t629;
t3 = t124 * t372 + t172 * t93 + t173 * t94 + t20 * t370 - t21 * t371 + t41 * t517;
t311 = Icges(6,6) * t554 + t639;
t313 = Icges(6,5) * t554 + t640;
t38 = -t100 * t371 + t126 * t517 + t370 * t99;
t42 = t1011 * t613 + t162 * t880 + t163 * t410 + t164 * t411 - t185 * t614 - t186 * t615;
t4 = t123 * t372 + t172 * t91 + t173 * t92 + t22 * t370 - t23 * t371 + t42 * t517;
t586 = t699 * t976 + t701 * t977 + ((t263 * t412 - t265 * t413) * t370 - (t262 * t412 - t264 * t413) * t371 + (t311 * t412 - t313 * t413) * t517 + (t124 * t554 - t884 * t93) * qJD(5) + ((-qJD(5) * t94 + t627) * t553 + t589) * t578) * t970 + (qJD(1) * t700 + t20 * t578 - t21 * t581) * t971 + (qJD(1) * t702 + t22 * t578 - t23 * t581) * t972 + ((t263 * t410 + t265 * t411) * t370 - (t262 * t410 + t264 * t411) * t371 + (t311 * t410 + t313 * t411) * t517 + (t123 * t554 - t885 * t92) * qJD(5) + ((-qJD(5) * t91 + t627) * t553 + t589) * t581) * t973 + (qJD(1) * t697 + t26 * t578 - t27 * t581) * t962 + (((-t263 * t579 - t265 * t576 + t208) * t370 - (-t262 * t579 - t264 * t576 + t210) * t371 + (-t311 * t579 - t313 * t576 - t613) * t517 + t126 * qJD(5)) * t554 + (-qJD(5) * t697 + t590) * t553) * t963 - t372 * t698 / 0.2e1 - t38 * t803 / 0.2e1 + (t1030 * t578 - t1031 * t581) * t975 + (t1032 * t578 - t1033 * t581) * t974 + (t1015 * t578 + t1016 * t581) * t967 + (t1047 * t581 + t1046 * t578 + (t1030 * t581 + t1031 * t578) * qJD(1)) * t966 + (t1045 * t581 + t1044 * t578 + (t1032 * t581 + t1033 * t578) * qJD(1)) * t965 + (-t1015 * t581 + t1016 * t578) * t964 + (-t1054 * t554 + t1055 * t553) * t953 + (t1039 * t581 + t1038 * t578 + (t1027 * t581 + t1028 * t578) * qJD(1)) * t952 + (t1027 * t578 - t1028 * t581) * t934 + t1007 * t804 / 0.2e1 + (t1041 * qJD(1) + t1029 * qJDD(1) + t1030 * t346 + t1031 * t347 + t1046 * t472 + t1047 * t473 + t4) * t961 + (t1040 * qJD(1) + t1065 * qJDD(1) + t1032 * t346 + t1033 * t347 + t1044 * t472 + t1045 * t473 + t3) * t960 + (t34 + t1043) * t756 + (t33 + t1042) * t755;
t455 = t484 * qJD(2);
t402 = t710 * t554;
t323 = t473 * t996;
t316 = -t769 + t781;
t256 = (t473 * t581 + t889) * t553;
t250 = rSges(6,1) * t412 + rSges(6,2) * t413;
t249 = rSges(6,1) * t410 - rSges(6,2) * t411;
t226 = -t581 * t661 + t888;
t220 = t226 * qJD(1);
t207 = t665 * qJD(2);
t119 = qJD(1) * t240 + qJDD(1) * t374 - t455 * t552 - t460 * t482 + t831;
t118 = -t455 * t807 + t461 * t482 + t837 * qJDD(1) + (-t241 - t459) * qJD(1);
t106 = t594 * t578 - t581 * t983;
t105 = t578 * t983 + t594 * t581;
t104 = -qJD(2) * t666 + t236 * t580 + t238 * t577;
t103 = -t667 * qJD(2) + t237 * t580 + t239 * t577;
t101 = qJD(1) * t722 + t473 * t707 + t643;
t85 = qJD(2) * t678 + t220;
t84 = qJD(2) * t679 + t797;
t77 = qJD(1) * t201 + qJDD(1) * t341 - t346 * t443 - t363 * t472 + t607;
t76 = t347 * t443 - t363 * t473 + t775 * qJDD(1) + (-t202 + t864) * qJD(1) + t660;
t50 = t201 * t473 + t202 * t472 + t340 * t346 - t341 * t347 + t652;
t28 = t203 * t472 - t343 * t346 + (t174 + t204) * t473 - t842 * t347 + t609;
t15 = (t174 + t283) * t473 + t115 * t370 + t116 * t371 + t172 * t219 - t173 * t217 + t282 * t472 - t346 * t428 - t833 * t347 + t609;
t1 = [t33 * t970 + (t101 * t772 + t102 * (-pkin(3) * t783 + t452 + t655 + t720) + t101 * (-t530 + (-t946 + (-rSges(5,3) - qJ(4)) * t554) * t572) * t578 + ((-t101 * rSges(5,1) + t102 * t622) * t578 + (t101 * (t622 + t945) - t102 * t582) * t581) * qJD(1) - (qJD(1) * t343 + t473 * t828 - t101 + t621) * t102 + (t40 - g(2)) * (t342 + t650) + (t39 - g(1)) * ((-t531 + (rSges(5,2) - pkin(3)) * t554) * t578 + t748 + t819)) * m(5) + (t222 * (t547 + t818) + (t482 * t902 - t900) * qJD(2) + ((-pkin(1) - t484) * t901 + (t221 * (-rSges(3,3) - pkin(6)) + t222 * t757) * t578) * qJD(1) - (-qJD(1) * t373 - t221 - t467 - t763) * t222 + (t119 - g(2)) * t292 + (t118 - g(1)) * (t757 * t578 + t570 + t816)) * m(3) - m(2) * (-g(1) * t483 + g(2) * t485) + t124 * t976 + t123 * t977 + t42 * t972 + (t84 - t797 + ((t581 * t729 + t141 - t852) * t581 + (t578 * t729 + t140 + t736) * t578) * qJD(2)) * t754 + ((-t712 + t772 + (-qJ(4) * t883 + (pkin(8) * t572 - qJD(4) + t943) * t553) * t578 + (t725 * t554 - t495 - t505 - t568) * qJD(1)) * t78 + (-t473 * t749 - t621 + t78 + t548 + t826 + t869 + (t792 * t886 - t786) * t581 + (-qJ(4) * t885 + t726 * t554 - t428 + t819) * qJD(1) + t1059) * t79 + (t25 - g(2)) * (t650 + t1026) + (t24 - g(1)) * (t836 + t428 + t819 + (t554 * t959 - t531) * t578)) * m(6) + t935 / 0.2e1 + (t104 + t105) * t753 + (t225 + t223) * t968 + (t41 + t33) * t971 + ((t1012 * t581 + t1008 + t1030 - t1064) * t473 + ((t1073 + t1081) * t581 + t1012 * t578 + t1031 + t1079) * t472 + t1043 - t1049) * t967 - t347 * t154 / 0.2e1 + (t103 + t106 + t85) * t752 + t905 / 0.2e1 + (t151 + t1028) * t974 + (t1027 + t1029) * t975 + (-qJD(2) * t661 - t1071 * t554 + t1072 * t553 + t450 * t580 + t451 * t577) * qJD(1) + t941 / 0.2e1 - t940 / 0.2e1 + (t220 + ((t139 - t306 + (t365 + t894) * t581 + t853) * t581 + t852 * t578) * qJD(2)) * t751 + t933 + (((t335 * t578 - t336 * t581) * t554 + (t329 * t581 + t330 * t578) * t553 + t1032 + t1061 + t1063) * t473 + (t334 * t885 + (-t335 * t581 - t336 * t578) * t554 - t331 * t882 + (t329 * t578 - t330 * t581) * t553 - t1008 + t1033 - t1060) * t472 + t1050) * t964 + (t1038 + t1041) * t966 + (-t1039 + t1040 + t1042) * t965 + (m(2) * (t483 ^ 2 + t485 ^ 2) + Icges(2,3) + t662 + t1092 * t554 + t1090 * t553) * qJDD(1) + (t226 + t224) * t969 + (t127 * t817 + t128 * (-t724 + t823) + (-t128 * t405 + t443 * t903) * t572 + ((-t127 * rSges(4,3) + t128 * (-t542 - t535)) * t578 + (t127 * (-t542 - t995) - t128 * t582) * t581) * qJD(1) - (-qJD(1) * t340 - t127 + t631 + t997) * t128 + (t77 - g(2)) * (t341 + t728) + (t76 - g(1)) * (-t340 + t819)) * m(4); ((t138 * t578 + t139 * t581) * qJD(1) + t704) * t752 + ((t140 * t578 + t141 * t581) * qJD(1) + t705) * t753 + ((-t807 * t888 - t810) * t581 + (t612 + (t581 * t887 + t591) * qJD(2)) * t578) * t751 + ((-t552 * t887 + t810) * t578 + (t612 + (t578 * t888 + t591) * qJD(2)) * t581) * t754 + t679 * t968 + t678 * t969 + (qJD(1) * t106 + qJD(2) * t704 + qJDD(1) * t225 + t138 * t461 + t139 * t460) * t960 + (qJD(1) * t105 + qJD(2) * t705 + qJDD(1) * t226 + t140 * t461 + t141 * t460) * t961 + (-t103 * t581 + t104 * t578 + (t223 * t578 + t224 * t581) * qJD(1)) * t952 + (-t223 * t581 + t224 * t578) * t934 + ((t577 * t822 + t580 * t821) * qJD(1) + (t623 * t577 + t580 * t982) * qJD(2)) * t953 + t85 * t755 + t586 + t84 * t756 + (-(-t79 * t764 + (-t580 * t703 + t63 * t747) * qJD(2)) * pkin(2) - g(1) * (-pkin(2) * t875 + t770) - g(2) * (-pkin(2) * t876 + t771) - g(3) * (t569 + t646) + t15 * (t647 + t856) + t63 * (t619 + t777) + t642 * t1006 + (t25 * t642 + t79 * t620 + (-t326 + t778) * t908) * t578 + (t219 * t803 + t318 * t759 + t581 * t620 + t1010) * t78 + t1017) * m(6) + (-g(1) * (t581 * t715 + t404 + t496) - g(2) * (t578 * t715 + t399 + t494) - g(3) * (t569 - t827) - t101 * (-qJD(1) * t399 - t323 + t714) - (-t102 * t764 + ((-t101 * t581 - t102 * t578) * t580 + t90 * t747) * qJD(2)) * pkin(2) - t592 + t101 * t832 + t28 * (t721 + t856) + t90 * (t648 + t777) + (t1000 * t718 + t101 * t654) * t581 + (t40 * t718 + t102 * t654 + (-t326 - t842) * t906) * t578) * m(5) + (-(-t128 * t764 + (t122 * t747 + t580 * t680) * qJD(2)) * pkin(2) - g(3) * (t995 + t569) - t992 * t630 + t50 * (t854 + t856) + t122 * (t777 + t779) + (t127 * t716 + (qJD(1) * t128 + t76) * t630) * t581 + (t77 * t630 + t128 * t716 + (t127 * t443 + t122 * (-t326 - t341)) * qJD(1)) * t578 + t1005) * m(4) + (g(1) * t425 + g(2) * t424 - g(3) * t484 - (t221 * t424 - t900) * qJD(1) - (t207 * (-t424 * t578 - t425 * t581) + t674 * t484) * qJD(2) + (qJD(2) * t673 + t373 * t460 - t374 * t461) * t665 + t207 * ((t373 * t581 - t374 * t578) * qJD(1) + t673) + t674 * t455 + (-t118 * t581 - t119 * t578 + (-t222 * t581 + t902) * qJD(1)) * t482) * m(3); t586 + (t15 * t647 + t63 * t619 + (t659 * t79 + t778 * t908) * t578 - g(1) * t770 - g(2) * t771 - g(3) * t646 + (t1003 + t1006) * t717 + (-(-t219 * t554 - t318 * t885) * qJD(5) + t659 * t581 + t1010) * t78 + t1017) * m(6) + (-t592 + t28 * t721 + t90 * t648 + (t102 * t855 - t842 * t906) * t578 - g(1) * t834 + g(2) * t835 + g(3) * t827 + (t1000 * t581 + t1002) * t828 + (-qJD(1) * t835 + t581 * t855 + t323 - t820 + t832) * t101) * m(5) + (g(1) * t405 + g(2) * t400 - g(3) * t995 + t50 * t854 + t122 * (-t341 * t809 + t779) + t680 * t363 + (-t77 * t578 - t76 * t581 + (-t128 * t581 + t903) * qJD(1)) * t443 + t1005) * m(4); (-m(5) - m(6)) * (-g(3) * t554 + t553 * t992) - m(5) * (t101 * t316 + t102 * t315 + t256 * t90) - m(6) * (t256 * t63 + t315 * t79 + t316 * t78) + 0.2e1 * ((t101 * t473 + t102 * t879 - t28) * t979 + (t473 * t78 + t79 * t879 - t15) * t978) * t554 + 0.2e1 * ((-t101 * t809 + t102 * t808 + t39 * t581 + t572 * t90 + t1002) * t979 + (t24 * t581 + t572 * t63 - t78 * t809 + t79 * t808 + t1003) * t978) * t553; -t33 * t767 / 0.2e1 + t4 * t880 / 0.2e1 + (t123 * t553 + t554 * t702) * t977 + ((-t572 * t702 + t42) * t553 + (-qJD(1) * t701 + t123 * t572 + t22 * t581 + t23 * t578) * t554) * t972 + t554 * t34 * t755 + t3 * t882 / 0.2e1 + (t124 * t553 + t554 * t700) * t976 + ((-t572 * t700 + t41) * t553 + (-qJD(1) * t699 + t124 * t572 + t20 * t581 + t21 * t578) * t554) * t971 + t38 * t883 / 0.2e1 + t553 * (t905 + t933 + t935 - t940 + t941) / 0.2e1 + t372 * (t126 * t553 + t554 * t697) / 0.2e1 + ((-t572 * t697 + t51) * t553 + (qJD(1) * t698 + t126 * t572 + t26 * t581 + t27 * t578) * t554) * t962 + (t605 * t410 - t411 * t606 + t581 * t618) * t973 + (t412 * t605 + t413 * t606 + t578 * t618) * t970 + (t626 * t553 + (t606 * t576 - t579 * t605) * t554) * t963 - t1007 * t886 / 0.2e1 + ((-t78 * t115 + t79 * t116 + t25 * t217 - t24 * t219 + (t63 * t675 + (-t578 * t78 + t581 * t79) * t318) * t572) * t553 + (t78 * (t166 * t578 - t219 * t572) + t79 * (-t166 * t581 + t217 * t572) - t15 * t675 + t63 * (t115 * t581 - t116 * t578 - t217 * t808 - t219 * t809) + (qJD(1) * t703 + t24 * t578 - t25 * t581) * t318) * t554 - t78 * (-t250 * t517 - t371 * t402) - t79 * (t249 * t517 - t370 * t402) - t63 * (t249 * t371 + t250 * t370) - g(1) * t249 - g(2) * t250 - g(3) * t402) * m(6);];
tau = t1;