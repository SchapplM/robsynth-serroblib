% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRR13_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR13_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:21:44
% EndTime: 2019-03-09 04:24:24
% DurationCPUTime: 141.89s
% Computational Cost: add. (122021->1575), mult. (358970->2054), div. (0->0), fcn. (436060->14), ass. (0->601)
t1020 = Icges(5,5) - Icges(4,6);
t1021 = Icges(5,4) - Icges(4,5);
t1072 = Icges(5,1) + Icges(4,3);
t887 = sin(pkin(7));
t888 = sin(pkin(6));
t774 = t888 * t887;
t913 = cos(qJ(3));
t736 = t913 * t774;
t914 = cos(qJ(1));
t889 = cos(pkin(12));
t891 = cos(pkin(6));
t778 = t891 * t889;
t886 = sin(pkin(12));
t911 = sin(qJ(1));
t690 = -t778 * t914 + t886 * t911;
t890 = cos(pkin(7));
t954 = t690 * t890;
t1000 = -t914 * t736 - t913 * t954;
t776 = t891 * t886;
t615 = t776 * t914 + t889 * t911;
t910 = sin(qJ(3));
t553 = t615 * t910 - t1000;
t735 = t910 * t774;
t952 = t914 * t735 + t910 * t954;
t556 = -t615 * t913 + t952;
t775 = t890 * t888;
t591 = t690 * t887 - t914 * t775;
t1009 = t1020 * t553 + t1021 * t556 + t1072 * t591;
t1019 = Icges(4,2) + Icges(5,3);
t1071 = Icges(4,4) + Icges(5,6);
t1081 = t1071 * t556;
t1011 = t1019 * t553 + t1020 * t591 + t1081;
t1023 = Icges(4,1) + Icges(5,2);
t1084 = t1071 * t553;
t1099 = t1021 * t591 + t1023 * t556 + t1084;
t773 = t888 * t886;
t988 = t889 * t775 + t891 * t887;
t587 = t910 * t773 - t913 * t988;
t588 = t773 * t913 + t910 * t988;
t614 = -t774 * t889 + t890 * t891;
t975 = t1009 * t614 + t1011 * t587 - t1099 * t588;
t692 = -t776 * t911 + t889 * t914;
t691 = t778 * t911 + t886 * t914;
t682 = t691 * t890;
t950 = t913 * t682 - t911 * t736;
t557 = t692 * t910 + t950;
t951 = -t910 * t682 + t911 * t735;
t558 = t692 * t913 + t951;
t1106 = t1011 * t557 - t1099 * t558;
t1105 = t1011 * t553 + t1099 * t556;
t681 = t691 * t887;
t733 = t911 * t775;
t592 = t681 + t733;
t1101 = t1009 * t592;
t649 = sin(qJ(5));
t912 = cos(qJ(5));
t466 = t553 * t912 - t591 * t649;
t997 = qJD(3) * t591;
t473 = -qJD(5) * t556 + t997;
t301 = qJD(6) * t466 - t473;
t470 = -t557 * t912 + t592 * t649;
t584 = qJD(3) * t592;
t474 = qJD(5) * t558 + t584;
t303 = qJD(6) * t470 + t474;
t551 = -t587 * t912 + t614 * t649;
t601 = qJD(3) * t614 + qJD(1);
t559 = qJD(5) * t588 + t601;
t440 = qJD(6) * t551 + t559;
t469 = t553 * t649 + t591 * t912;
t648 = sin(qJ(6));
t650 = cos(qJ(6));
t333 = t469 * t648 + t556 * t650;
t336 = t469 * t650 - t556 * t648;
t160 = Icges(7,5) * t336 - Icges(7,6) * t333 - Icges(7,3) * t466;
t879 = Icges(7,4) * t336;
t163 = -Icges(7,2) * t333 - Icges(7,6) * t466 + t879;
t329 = Icges(7,4) * t333;
t166 = Icges(7,1) * t336 - Icges(7,5) * t466 - t329;
t471 = t557 * t649 + t592 * t912;
t337 = -t471 * t648 + t558 * t650;
t338 = t471 * t650 + t558 * t648;
t57 = t470 * t160 + t337 * t163 + t338 * t166;
t161 = Icges(7,5) * t338 + Icges(7,6) * t337 + Icges(7,3) * t470;
t878 = Icges(7,4) * t338;
t164 = Icges(7,2) * t337 + Icges(7,6) * t470 + t878;
t330 = Icges(7,4) * t337;
t167 = Icges(7,1) * t338 + Icges(7,5) * t470 + t330;
t58 = t470 * t161 + t337 * t164 + t338 * t167;
t552 = t587 * t649 + t614 * t912;
t458 = -t552 * t648 + t588 * t650;
t459 = t552 * t650 + t588 * t648;
t243 = Icges(7,5) * t459 + Icges(7,6) * t458 + Icges(7,3) * t551;
t877 = Icges(7,4) * t459;
t244 = Icges(7,2) * t458 + Icges(7,6) * t551 + t877;
t457 = Icges(7,4) * t458;
t245 = Icges(7,1) * t459 + Icges(7,5) * t551 + t457;
t80 = t243 * t470 + t244 * t337 + t245 * t338;
t17 = -t301 * t57 + t303 * t58 + t80 * t440;
t1074 = t243 * t466 + t244 * t333 - t245 * t336;
t55 = -t160 * t466 - t163 * t333 + t166 * t336;
t56 = -t161 * t466 - t164 * t333 + t336 * t167;
t16 = -t1074 * t440 - t301 * t55 + t303 * t56;
t61 = t160 * t551 + t163 * t458 + t166 * t459;
t995 = t1009 * t591;
t1018 = t1105 + t995;
t1008 = t1020 * t557 - t1021 * t558 + t1072 * t592;
t1039 = t1008 * t592;
t1095 = t1018 + t1039;
t767 = -t336 * rSges(7,1) + rSges(7,2) * t333;
t169 = -rSges(7,3) * t466 - t767;
t246 = rSges(7,1) * t459 + rSges(7,2) * t458 + rSges(7,3) * t551;
t906 = t469 * pkin(5);
t294 = -pkin(11) * t466 + t906;
t418 = pkin(5) * t552 + pkin(11) * t551;
t1094 = -t169 * t440 - t246 * t301 - t294 * t559 + t418 * t473;
t248 = Icges(6,5) * t469 + Icges(6,6) * t466 - Icges(6,3) * t556;
t882 = Icges(6,4) * t469;
t251 = Icges(6,2) * t466 - Icges(6,6) * t556 + t882;
t460 = Icges(6,4) * t466;
t254 = Icges(6,1) * t469 - Icges(6,5) * t556 + t460;
t88 = -t248 * t556 + t251 * t466 + t254 * t469;
t249 = Icges(6,5) * t471 - Icges(6,6) * t470 + Icges(6,3) * t558;
t881 = Icges(6,4) * t471;
t252 = -Icges(6,2) * t470 + Icges(6,6) * t558 + t881;
t461 = Icges(6,4) * t470;
t255 = Icges(6,1) * t471 + Icges(6,5) * t558 - t461;
t89 = -t249 * t556 + t252 * t466 + t469 * t255;
t340 = Icges(6,5) * t552 - Icges(6,6) * t551 + Icges(6,3) * t588;
t880 = Icges(6,4) * t552;
t341 = -Icges(6,2) * t551 + Icges(6,6) * t588 + t880;
t528 = Icges(6,4) * t551;
t342 = Icges(6,1) * t552 + Icges(6,5) * t588 - t528;
t99 = -t340 * t556 + t341 * t466 + t342 * t469;
t30 = t473 * t88 + t474 * t89 + t559 * t99;
t100 = t340 * t558 - t341 * t470 + t342 * t471;
t90 = t558 * t248 - t470 * t251 + t471 * t254;
t91 = t558 * t249 - t470 * t252 + t471 * t255;
t31 = t100 * t559 + t473 * t90 + t474 * t91;
t96 = t248 * t588 - t251 * t551 + t254 * t552;
t1083 = t1071 * t557;
t1082 = t1071 * t558;
t769 = -t469 * rSges(6,1) - rSges(6,2) * t466;
t257 = -rSges(6,3) * t556 - t769;
t343 = rSges(6,1) * t552 - rSges(6,2) * t551 + rSges(6,3) * t588;
t1077 = -t257 * t559 + t343 * t473;
t1010 = t1019 * t557 + t1020 * t592 - t1082;
t1043 = t1021 * t592 - t1023 * t558 + t1083;
t1076 = t1010 * t553 + t1043 * t556;
t1005 = t1020 * t587 - t1021 * t588 + t1072 * t614;
t1070 = t1071 * t587;
t1042 = t1021 * t614 - t1023 * t588 + t1070;
t1069 = t1071 * t588;
t991 = t1019 * t587 + t1020 * t614 - t1069;
t973 = t1005 * t591 + t1042 * t556 + t553 * t991;
t1073 = t601 * t973;
t1040 = qJD(4) * t553;
t439 = t558 * pkin(3) + t557 * qJ(4);
t521 = pkin(3) * t588 + qJ(4) * t587;
t631 = pkin(9) * t733;
t675 = pkin(9) * t681;
t666 = t692 * pkin(2) + t631 + t675;
t780 = t914 * t888;
t641 = qJD(2) * t780;
t647 = t914 * pkin(1);
t779 = t888 * t911;
t712 = qJ(2) * t779 + t647;
t949 = qJD(1) * t712 - t641;
t814 = qJD(1) * t666 + t949;
t1068 = t601 * t439 - t521 * t584 + t1040 + t814;
t477 = t592 * pkin(4) + t558 * pkin(10);
t560 = pkin(4) * t614 + pkin(10) * t588;
t1054 = t601 * t477 - t560 * t584 + t1068;
t367 = rSges(4,1) * t556 + t553 * rSges(4,2) - t591 * rSges(4,3);
t489 = rSges(4,1) * t588 - rSges(4,2) * t587 + rSges(4,3) * t614;
t1051 = t367 * t601 + t489 * t997;
t1049 = t31 / 0.2e1;
t978 = t1010 * t557 - t1043 * t558 + t1039;
t1047 = t1008 * t591 + t1076;
t1017 = t1101 + t1106;
t608 = t615 * qJD(1);
t400 = qJD(1) * t1000 + qJD(3) * t558 - t608 * t910;
t401 = qJD(1) * t952 - qJD(3) * t557 - t608 * t913;
t575 = t591 * qJD(1);
t1015 = t1019 * t400 - t1020 * t575 - t1071 * t401;
t609 = t692 * qJD(1);
t402 = qJD(1) * t950 - qJD(3) * t556 + t609 * t910;
t403 = qJD(1) * t951 - qJD(3) * t553 + t609 * t913;
t576 = t592 * qJD(1);
t1014 = t1019 * t402 + t1020 * t576 - t1071 * t403;
t1046 = t1020 * t400 - t1021 * t401 - t1072 * t575;
t1013 = -t1021 * t575 - t1023 * t401 + t1071 * t400;
t1012 = t1021 * t576 - t1023 * t403 + t1071 * t402;
t1045 = t1020 * t402 - t1021 * t403 + t1072 * t576;
t571 = t587 * qJD(3);
t572 = t588 * qJD(3);
t1004 = t1019 * t572 + t1071 * t571;
t1003 = t1023 * t571 + t1071 * t572;
t1002 = t1020 * t572 + t1021 * t571;
t538 = t553 * qJ(4);
t905 = pkin(3) * t556;
t434 = t538 - t905;
t170 = t338 * rSges(7,1) + t337 * rSges(7,2) + t470 * rSges(7,3);
t296 = t471 * pkin(5) + t470 * pkin(11);
t42 = t170 * t440 - t246 * t303 + t296 * t559 - t418 * t474 + t1054;
t258 = t471 * rSges(6,1) - t470 * rSges(6,2) + t558 * rSges(6,3);
t85 = t258 * t559 - t343 * t474 + t1054;
t371 = t592 * rSges(5,1) - t558 * rSges(5,2) + t557 * rSges(5,3);
t490 = rSges(5,1) * t614 - rSges(5,2) * t588 + rSges(5,3) * t587;
t105 = t371 * t601 - t490 * t584 + t1068;
t799 = t997 / 0.2e1;
t972 = t1005 * t592 - t1042 * t558 + t557 * t991;
t771 = -t591 * rSges(5,1) - t553 * rSges(5,3);
t370 = rSges(5,2) * t556 - t771;
t1038 = t1011 * t572 - t1012 * t588 + t1014 * t587 + t1045 * t614 + t1099 * t571;
t1037 = t1010 * t572 - t1013 * t588 + t1015 * t587 + t1043 * t571 + t1046 * t614;
t1016 = t1002 * t614 - t1003 * t588 + t1004 * t587 + t1042 * t571 + t572 * t991;
t974 = t1008 * t614 + t1010 * t587 - t1043 * t588;
t1001 = t591 * pkin(9);
t999 = t1017 * t591 + t592 * t978;
t998 = t1018 * t591 + t1047 * t592;
t476 = t591 * pkin(4) - pkin(10) * t556;
t993 = t1045 * t591 + t1046 * t592;
t992 = t972 * t601;
t989 = -t615 * pkin(2) - t1001;
t239 = qJD(5) * t471 - t400 * t912 - t575 * t649;
t567 = qJD(3) * t575;
t325 = qJD(5) * t401 - t567;
t157 = qJD(6) * t239 + t325;
t241 = qJD(5) * t469 - t402 * t912 + t576 * t649;
t568 = qJD(3) * t576;
t326 = qJD(5) * t403 + t568;
t158 = qJD(6) * t241 + t326;
t455 = -qJD(5) * t551 + t572 * t649;
t259 = -qJD(6) * t459 - t455 * t648 - t571 * t650;
t260 = qJD(6) * t458 + t455 * t650 - t571 * t648;
t456 = qJD(5) * t552 - t572 * t912;
t126 = Icges(7,5) * t260 + Icges(7,6) * t259 + Icges(7,3) * t456;
t127 = Icges(7,4) * t260 + Icges(7,2) * t259 + Icges(7,6) * t456;
t128 = Icges(7,1) * t260 + Icges(7,4) * t259 + Icges(7,5) * t456;
t240 = -qJD(5) * t470 + t400 * t649 - t575 * t912;
t142 = -qJD(6) * t338 - t240 * t648 + t401 * t650;
t143 = qJD(6) * t337 + t240 * t650 + t401 * t648;
t19 = t126 * t470 + t127 * t337 + t128 * t338 + t142 * t244 + t143 * t245 + t239 * t243;
t831 = qJD(5) * t571;
t377 = qJD(6) * t456 - t831;
t242 = qJD(5) * t466 + t402 * t649 + t576 * t912;
t144 = -qJD(6) * t336 - t242 * t648 + t403 * t650;
t145 = -qJD(6) * t333 + t242 * t650 + t403 * t648;
t64 = Icges(7,5) * t145 + Icges(7,6) * t144 + Icges(7,3) * t241;
t66 = Icges(7,4) * t145 + Icges(7,2) * t144 + Icges(7,6) * t241;
t68 = Icges(7,1) * t145 + Icges(7,4) * t144 + Icges(7,5) * t241;
t8 = t142 * t163 + t143 * t166 + t160 * t239 + t337 * t66 + t338 * t68 + t470 * t64;
t63 = Icges(7,5) * t143 + Icges(7,6) * t142 + Icges(7,3) * t239;
t65 = Icges(7,4) * t143 + Icges(7,2) * t142 + Icges(7,6) * t239;
t67 = Icges(7,1) * t143 + Icges(7,4) * t142 + Icges(7,5) * t239;
t9 = t142 * t164 + t143 * t167 + t161 * t239 + t337 * t65 + t338 * t67 + t470 * t63;
t1 = t157 * t58 + t158 * t57 + t19 * t440 - t301 * t8 + t303 * t9 + t377 * t80;
t107 = Icges(6,5) * t242 - Icges(6,6) * t241 + Icges(6,3) * t403;
t109 = Icges(6,4) * t242 - Icges(6,2) * t241 + Icges(6,6) * t403;
t111 = Icges(6,1) * t242 - Icges(6,4) * t241 + Icges(6,5) * t403;
t22 = t107 * t558 - t109 * t470 + t111 * t471 - t239 * t251 + t240 * t254 + t248 * t401;
t106 = Icges(6,5) * t240 - Icges(6,6) * t239 + Icges(6,3) * t401;
t108 = Icges(6,4) * t240 - Icges(6,2) * t239 + Icges(6,6) * t401;
t110 = Icges(6,1) * t240 - Icges(6,4) * t239 + Icges(6,5) * t401;
t23 = t106 * t558 - t108 * t470 + t110 * t471 - t239 * t252 + t240 * t255 + t249 * t401;
t261 = Icges(6,5) * t455 - Icges(6,6) * t456 - Icges(6,3) * t571;
t262 = Icges(6,4) * t455 - Icges(6,2) * t456 - Icges(6,6) * t571;
t263 = Icges(6,1) * t455 - Icges(6,4) * t456 - Icges(6,5) * t571;
t36 = -t239 * t341 + t240 * t342 + t261 * t558 - t262 * t470 + t263 * t471 + t340 * t401;
t987 = -t100 * t831 + t22 * t473 + t23 * t474 + t325 * t91 + t326 * t90 + t36 * t559 + t1;
t10 = t144 * t163 + t145 * t166 + t160 * t241 - t333 * t66 + t336 * t68 - t466 * t64;
t11 = t144 * t164 + t145 * t167 + t161 * t241 - t333 * t65 + t336 * t67 - t466 * t63;
t20 = -t126 * t466 - t127 * t333 + t128 * t336 + t144 * t244 + t145 * t245 + t241 * t243;
t2 = -t10 * t301 - t1074 * t377 + t11 * t303 + t157 * t56 + t158 * t55 + t20 * t440;
t24 = -t107 * t556 + t109 * t466 + t111 * t469 - t241 * t251 + t242 * t254 + t248 * t403;
t25 = -t106 * t556 + t108 * t466 + t110 * t469 - t241 * t252 + t242 * t255 + t249 * t403;
t37 = -t241 * t341 + t242 * t342 - t261 * t556 + t262 * t466 + t263 * t469 + t340 * t403;
t986 = t24 * t473 + t25 * t474 + t325 * t89 + t326 * t88 + t37 * t559 - t831 * t99 + t2;
t62 = t161 * t551 + t164 * t458 + t167 * t459;
t897 = t62 * t157;
t898 = t61 * t158;
t13 = t161 * t456 + t164 * t259 + t167 * t260 + t458 * t65 + t459 * t67 + t551 * t63;
t901 = t13 * t303;
t12 = t160 * t456 + t163 * t259 + t166 * t260 + t458 * t66 + t459 * t68 + t551 * t64;
t902 = t12 * t301;
t27 = t126 * t551 + t127 * t458 + t128 * t459 + t243 * t456 + t244 * t259 + t245 * t260;
t86 = t243 * t551 + t244 * t458 + t245 * t459;
t903 = t27 * t440 + t86 * t377;
t3 = t897 + t898 + t901 - t902 + t903;
t117 = t340 * t588 - t341 * t551 + t342 * t552;
t60 = t261 * t588 - t262 * t551 + t263 * t552 - t340 * t571 - t341 * t456 + t342 * t455;
t892 = -t117 * t831 + t60 * t559;
t97 = t249 * t588 - t252 * t551 + t255 * t552;
t895 = t97 * t325;
t896 = t96 * t326;
t29 = t106 * t588 - t108 * t551 + t110 * t552 - t249 * t571 - t252 * t456 + t255 * t455;
t899 = t29 * t474;
t28 = t107 * t588 - t109 * t551 + t111 * t552 - t248 * t571 - t251 * t456 + t254 * t455;
t900 = t28 * t473;
t985 = t892 + t895 + t896 + t899 + t900 + t3;
t984 = t30 + t16;
t983 = t31 + t17;
t982 = qJD(3) * t998 + t1073;
t981 = qJD(3) * t999 + t992;
t980 = t1002 * t592 - t1003 * t558 + t1004 * t557 - t1005 * t575 - t1042 * t401 + t400 * t991;
t979 = t1002 * t591 + t1003 * t556 + t1004 * t553 + t1005 * t576 - t1042 * t403 + t402 * t991;
t976 = t1016 * t601;
t968 = t1020 * t588 + t1021 * t587;
t685 = t690 * rSges(3,2);
t967 = -t615 * rSges(3,1) + rSges(3,3) * t780 + t685;
t966 = t1037 * t592 + t1038 * t591 - t575 * t974 + t576 * t975;
t965 = (t1011 * t400 - t1012 * t558 + t1014 * t557 - t1099 * t401) * t591 + t1017 * t576 - (t995 + t978) * t575 + (-t1008 * t575 + t1010 * t400 - t1013 * t558 + t1015 * t557 - t1043 * t401 + t993) * t592;
t964 = (t1010 * t402 + t1013 * t556 + t1015 * t553 - t1043 * t403) * t592 + t1095 * t576 - t1047 * t575 + (t1009 * t576 + t1011 * t402 + t1012 * t556 + t1014 * t553 - t1099 * t403 + t993) * t591;
t963 = (t1019 * t558 + t1043 + t1083) * t592 + (-t1019 * t556 + t1084 + t1099) * t591;
t962 = (-t1023 * t557 + t1010 - t1082) * t592 + (-t1023 * t553 + t1011 + t1081) * t591;
t961 = (t1020 * t558 + t1021 * t557) * t592 + (-t1020 * t556 + t1021 * t553) * t591;
t960 = -t557 * t601 + t587 * t584;
t959 = t553 * t601 - t587 * t997;
t958 = t1019 * t588 + t1042 + t1070;
t957 = -t1023 * t587 - t1069 + t991;
t956 = t1017 - t1101;
t824 = t911 * pkin(1);
t623 = -qJ(2) * t780 + t824;
t953 = -t623 + t967;
t948 = -t17 / 0.2e1;
t947 = t157 / 0.2e1;
t946 = t158 / 0.2e1;
t945 = t301 / 0.2e1;
t944 = -t301 / 0.2e1;
t943 = -t303 / 0.2e1;
t942 = t303 / 0.2e1;
t941 = t325 / 0.2e1;
t940 = t326 / 0.2e1;
t939 = t377 / 0.2e1;
t936 = -t440 / 0.2e1;
t935 = t440 / 0.2e1;
t934 = -t473 / 0.2e1;
t933 = t473 / 0.2e1;
t932 = -t474 / 0.2e1;
t931 = t474 / 0.2e1;
t928 = -t559 / 0.2e1;
t927 = t559 / 0.2e1;
t917 = rSges(5,2) - pkin(3);
t916 = -rSges(6,3) - pkin(3);
t915 = -rSges(7,3) - pkin(11);
t909 = pkin(10) * t571;
t908 = t242 * pkin(5);
t907 = t403 * pkin(3);
t150 = t240 * pkin(5) + t239 * pkin(11);
t69 = t143 * rSges(7,1) + t142 * rSges(7,2) + t239 * rSges(7,3);
t894 = t150 + t69;
t151 = t241 * pkin(11) + t908;
t768 = -t145 * rSges(7,1) - t144 * rSges(7,2);
t70 = t241 * rSges(7,3) - t768;
t893 = t151 + t70;
t873 = t553 * t592;
t870 = t648 * t649;
t869 = t649 * t650;
t129 = rSges(7,1) * t260 + rSges(7,2) * t259 + rSges(7,3) * t456;
t277 = pkin(5) * t455 + pkin(11) * t456;
t868 = t129 + t277;
t867 = t169 + t294;
t866 = t170 + t296;
t850 = -t402 * qJ(4) - t1040;
t207 = -t850 + t907;
t865 = t592 * t207 - t575 * t434;
t536 = qJD(4) * t557;
t206 = t401 * pkin(3) + t400 * qJ(4) + t536;
t327 = -t575 * pkin(4) + t401 * pkin(10);
t864 = -t206 - t327;
t772 = -t576 * rSges(5,1) - t402 * rSges(5,3);
t234 = -t403 * rSges(5,2) - t772;
t863 = -t207 - t234;
t328 = t576 * pkin(4) + t403 * pkin(10);
t862 = -t207 - t328;
t861 = t246 + t418;
t860 = -(-t557 * t591 + t873) * qJD(3) + t572;
t432 = -pkin(3) * t553 - qJ(4) * t556;
t859 = qJD(4) * t588 + t432 * t584;
t858 = -t370 - t434;
t857 = -t371 - t439;
t372 = t592 * t434;
t856 = t592 * t476 + t372;
t855 = t400 + t959;
t854 = t402 + t960;
t578 = qJD(4) * t587;
t443 = -pkin(3) * t571 + qJ(4) * t572 + t578;
t853 = t591 * t443 + t576 * t521;
t437 = -pkin(3) * t557 + qJ(4) * t558;
t852 = -qJD(4) * t556 + t601 * t437;
t381 = t614 * t439;
t851 = t614 * t477 + t381;
t849 = -t434 - t476;
t848 = -t439 - t477;
t519 = -pkin(3) * t587 + qJ(4) * t588;
t847 = qJD(4) * t558 + t519 * t997;
t454 = t591 * t521;
t846 = t591 * t560 + t454;
t841 = -t490 - t521;
t840 = -t521 - t560;
t719 = -t623 + t989;
t640 = qJD(2) * t779;
t747 = qJD(1) * t780;
t740 = qJ(2) * t747 - qJD(1) * t824 + t640;
t838 = (t640 + t740) * qJD(1);
t836 = -qJD(1) * t623 + t640;
t833 = qJD(5) * t553;
t832 = qJD(5) * t557;
t830 = qJD(5) * t587;
t829 = 2 * m(3);
t828 = 2 * m(4);
t827 = 2 * m(5);
t826 = 2 * m(6);
t825 = 2 * m(7);
t522 = t592 * t909;
t823 = pkin(10) * t997;
t822 = qJD(4) * t572 + t207 * t584 - t434 * t567;
t189 = t614 * t206;
t821 = t614 * t327 + t189 + t522;
t112 = t240 * rSges(6,1) - t239 * rSges(6,2) + t401 * rSges(6,3);
t820 = -t258 + t848;
t819 = -t343 + t840;
t646 = qJD(2) * t891;
t818 = t434 * t584 + t578 + t646;
t817 = t576 * t560 + t853;
t231 = t401 * rSges(4,1) - t400 * rSges(4,2) - t575 * rSges(4,3);
t665 = -t608 * pkin(2) - qJD(1) * t1001;
t816 = qJD(1) * t665 + t838;
t368 = t558 * rSges(4,1) - t557 * rSges(4,2) + t592 * rSges(4,3);
t815 = qJD(1) * t989 + t836;
t233 = -t575 * rSges(5,1) - t401 * rSges(5,2) + t400 * rSges(5,3);
t813 = -t608 * rSges(3,1) + rSges(3,3) * t747 + qJD(1) * t685;
t611 = t691 * rSges(3,2);
t812 = t692 * rSges(3,1) + rSges(3,3) * t779 - t611;
t811 = t556 * t912;
t810 = t558 * t912;
t809 = t588 * t912;
t804 = qJD(6) * t912;
t803 = -t567 / 0.2e1;
t802 = t568 / 0.2e1;
t800 = -t997 / 0.2e1;
t798 = -t584 / 0.2e1;
t797 = t584 / 0.2e1;
t795 = -t831 / 0.2e1;
t794 = t609 * pkin(2) + qJD(1) * t675;
t792 = -t828 / 0.2e1;
t791 = t828 / 0.2e1;
t790 = -t827 / 0.2e1;
t789 = t827 / 0.2e1;
t788 = -t826 / 0.2e1;
t787 = t826 / 0.2e1;
t786 = -t825 / 0.2e1;
t785 = t825 / 0.2e1;
t783 = t848 - t866;
t782 = t592 * t328 - t575 * t476 + t865;
t781 = t840 - t861;
t770 = -t242 * rSges(6,1) + t241 * rSges(6,2);
t766 = -rSges(7,1) * t650 + rSges(7,2) * t648;
t765 = qJD(4) * t402 + t601 * t206 + t816;
t764 = -t434 * t601 + t536 + t815;
t762 = -Icges(7,1) * t650 + Icges(7,4) * t648;
t761 = -Icges(7,4) * t650 + Icges(7,2) * t648;
t760 = -Icges(7,5) * t650 + Icges(7,6) * t648;
t755 = -t367 * t592 - t368 * t591;
t504 = -rSges(4,1) * t571 - rSges(4,2) * t572;
t752 = t489 * t575 - t504 * t592;
t751 = t489 * t576 + t504 * t591;
t503 = rSges(5,2) * t571 + rSges(5,3) * t572;
t750 = t490 * t576 + t503 * t591;
t637 = qJD(1) * t641;
t749 = t637 + (-qJD(1) * t631 - t794 - t949) * qJD(1);
t748 = qJD(1) * t719 + t640;
t746 = qJD(1) * t779;
t745 = pkin(5) * t649 - pkin(11) * t912;
t744 = rSges(6,1) * t649 + rSges(6,2) * t912;
t743 = Icges(6,1) * t649 + Icges(6,4) * t912;
t742 = Icges(6,4) * t649 + Icges(6,2) * t912;
t741 = Icges(6,5) * t649 + Icges(6,6) * t912;
t739 = t370 * t592 + t591 * t857;
t728 = t521 * t997 + t536 + t748;
t232 = t403 * rSges(4,1) - t402 * rSges(4,2) + t576 * rSges(4,3);
t727 = -(-Icges(7,5) * t333 - Icges(7,6) * t336) * t301 + (Icges(7,5) * t337 - Icges(7,6) * t338) * t303 + (Icges(7,5) * t458 - Icges(7,6) * t459) * t440;
t726 = (Icges(6,5) * t466 - Icges(6,6) * t469) * t473 + (-Icges(6,5) * t470 - Icges(6,6) * t471) * t474 + (-Icges(6,5) * t551 - Icges(6,6) * t552) * t559;
t721 = (-t443 - t503) * t592 - t841 * t575;
t720 = qJD(4) * t400 + t443 * t997 + t521 * t568 + t749;
t714 = t476 * t584 + t848 * t997 + t818;
t711 = -t538 + t719;
t704 = pkin(10) * t960 - t519 * t584 + t852;
t703 = pkin(10) * t959 - t432 * t601 + t847;
t702 = -t231 * t591 + t232 * t592 + t367 * t575 - t368 * t576;
t701 = t557 * t823 + (-pkin(10) * t873 - t437 * t591) * qJD(3) + t859;
t700 = -t476 * t601 - t840 * t997 + t764;
t697 = (Icges(7,1) * t337 - t164 - t878) * t303 - (-Icges(7,1) * t333 - t163 - t879) * t301 + (Icges(7,1) * t458 - t244 - t877) * t440;
t696 = (-Icges(7,2) * t338 + t167 + t330) * t303 - (-Icges(7,2) * t336 + t166 - t329) * t301 + (-Icges(7,2) * t459 + t245 + t457) * t440;
t695 = (-Icges(6,1) * t470 - t252 - t881) * t474 + (Icges(6,1) * t466 - t251 - t882) * t473 + (-Icges(6,1) * t551 - t341 - t880) * t559;
t694 = (Icges(6,2) * t471 - t255 + t461) * t474 + (Icges(6,2) * t469 - t254 - t460) * t473 + (Icges(6,2) * t552 - t342 + t528) * t559;
t693 = t560 * t997 + t601 * t849 + t728;
t689 = -t476 + t711;
t687 = t234 * t592 - t370 * t575 + (-t206 - t233) * t591 + t857 * t576;
t686 = t560 * t568 + t601 * t862 + t720;
t683 = t328 * t584 - t476 * t567 + (t576 * t848 + t591 * t864) * qJD(3) + t822;
t676 = t601 * t327 + t765 + (-t443 * t592 - t575 * t840 + t522) * qJD(3);
t664 = t666 + t712;
t662 = t641 + (-t631 - t712) * qJD(1) - t794;
t661 = t439 + t664;
t660 = (Icges(7,3) * t471 + t164 * t648 - t167 * t650 + t470 * t760) * t303 - (Icges(7,3) * t469 + t163 * t648 - t166 * t650 - t466 * t760) * t301 + (Icges(7,3) * t552 + t244 * t648 - t245 * t650 + t551 * t760) * t440;
t659 = t609 * rSges(3,1) + rSges(3,3) * t746 - qJD(1) * t611;
t658 = t662 + t850;
t657 = t665 + t740;
t656 = t477 + t661;
t655 = (-Icges(6,3) * t557 + t252 * t912 + t255 * t649 + t558 * t741) * t474 + (-Icges(6,3) * t553 + t251 * t912 + t254 * t649 - t556 * t741) * t473 + (-Icges(6,3) * t587 + t341 * t912 + t342 * t649 + t588 * t741) * t559;
t654 = -t328 + t658;
t653 = t206 + t657;
t652 = t327 + t653;
t520 = -rSges(4,1) * t587 - rSges(4,2) * t588;
t518 = rSges(5,2) * t587 + rSges(5,3) * t588;
t510 = qJD(1) * t812 + t949;
t509 = qJD(1) * t953 + t640;
t508 = -t588 * t804 - t830;
t507 = t745 * t588;
t506 = -t587 * t648 + t588 * t869;
t505 = -t587 * t650 - t588 * t870;
t480 = t637 + (-t949 - t659) * qJD(1);
t479 = qJD(1) * t813 + t838;
t448 = -t587 * rSges(6,3) + t588 * t744;
t446 = -Icges(6,5) * t587 + t588 * t743;
t445 = -Icges(6,6) * t587 + t588 * t742;
t438 = -rSges(4,1) * t557 - rSges(4,2) * t558;
t436 = rSges(5,2) * t557 + rSges(5,3) * t558;
t433 = -rSges(4,1) * t553 + rSges(4,2) * t556;
t431 = rSges(5,2) * t553 - rSges(5,3) * t556;
t417 = -pkin(5) * t551 + pkin(11) * t552;
t416 = -rSges(6,1) * t551 - rSges(6,2) * t552;
t412 = -t558 * t804 - t832;
t411 = t556 * t804 - t833;
t410 = t745 * t558;
t409 = t745 * t556;
t407 = -t557 * t648 + t558 * t869;
t406 = -t557 * t650 - t558 * t870;
t405 = -t553 * t648 - t556 * t869;
t404 = -t553 * t650 + t556 * t870;
t323 = -t557 * rSges(6,3) + t558 * t744;
t322 = -t553 * rSges(6,3) - t556 * t744;
t320 = -Icges(6,5) * t557 + t558 * t743;
t319 = -Icges(6,5) * t553 - t556 * t743;
t318 = -Icges(6,6) * t557 + t558 * t742;
t317 = -Icges(6,6) * t553 - t556 * t742;
t308 = rSges(7,3) * t552 + t551 * t766;
t307 = Icges(7,5) * t552 + t551 * t762;
t306 = Icges(7,6) * t552 + t551 * t761;
t300 = t506 * rSges(7,1) + t505 * rSges(7,2) - rSges(7,3) * t809;
t299 = Icges(7,1) * t506 + Icges(7,4) * t505 - Icges(7,5) * t809;
t298 = Icges(7,4) * t506 + Icges(7,2) * t505 - Icges(7,6) * t809;
t297 = Icges(7,5) * t506 + Icges(7,6) * t505 - Icges(7,3) * t809;
t295 = -pkin(5) * t470 + pkin(11) * t471;
t293 = pkin(5) * t466 + pkin(11) * t469;
t291 = -rSges(6,1) * t470 - rSges(6,2) * t471;
t290 = rSges(6,1) * t466 - rSges(6,2) * t469;
t281 = rSges(7,1) * t458 - rSges(7,2) * t459;
t264 = rSges(6,1) * t455 - rSges(6,2) * t456 - rSges(6,3) * t571;
t218 = rSges(7,3) * t471 + t470 * t766;
t217 = rSges(7,3) * t469 - t466 * t766;
t216 = Icges(7,5) * t471 + t470 * t762;
t215 = Icges(7,5) * t469 - t466 * t762;
t214 = Icges(7,6) * t471 + t470 * t761;
t213 = Icges(7,6) * t469 - t466 * t761;
t205 = t407 * rSges(7,1) + t406 * rSges(7,2) - rSges(7,3) * t810;
t204 = t405 * rSges(7,1) + t404 * rSges(7,2) + rSges(7,3) * t811;
t203 = Icges(7,1) * t407 + Icges(7,4) * t406 - Icges(7,5) * t810;
t202 = Icges(7,1) * t405 + Icges(7,4) * t404 + Icges(7,5) * t811;
t201 = Icges(7,4) * t407 + Icges(7,2) * t406 - Icges(7,6) * t810;
t200 = Icges(7,4) * t405 + Icges(7,2) * t404 + Icges(7,6) * t811;
t199 = Icges(7,5) * t407 + Icges(7,6) * t406 - Icges(7,3) * t810;
t198 = Icges(7,5) * t405 + Icges(7,6) * t404 + Icges(7,3) * t811;
t197 = rSges(7,1) * t337 - rSges(7,2) * t338;
t196 = -rSges(7,1) * t333 - rSges(7,2) * t336;
t179 = t368 * t601 - t489 * t584 + t814;
t178 = t1051 + t748;
t177 = qJD(3) * t755 + t646;
t116 = qJD(3) * t751 - t232 * t601 + t749;
t115 = qJD(3) * t752 + t231 * t601 + t816;
t113 = t403 * rSges(6,3) - t770;
t104 = t490 * t997 + t601 * t858 + t728;
t101 = qJD(3) * t739 + t818;
t87 = t702 * qJD(3);
t84 = t1077 + t693;
t78 = t257 * t474 - t258 * t473 + t714;
t72 = qJD(3) * t750 + t601 * t863 + t720;
t71 = qJD(3) * t721 + t233 * t601 + t765;
t41 = t1094 + t693;
t38 = qJD(3) * t687 + t822;
t35 = t169 * t303 + t170 * t301 + t294 * t474 - t296 * t473 + t714;
t34 = t117 * t559 + t473 * t96 + t474 * t97;
t33 = -t113 * t559 + t264 * t473 + t326 * t343 + (qJD(5) * t257 - t823) * t571 + t686;
t32 = t112 * t559 - t258 * t831 - t264 * t474 - t325 * t343 + t676;
t21 = -t112 * t473 + t113 * t474 + t257 * t325 - t258 * t326 + t683;
t18 = -t301 * t61 + t303 * t62 + t440 * t86;
t15 = -t301 * t129 - t559 * t151 + t158 * t246 - t377 * t169 + t473 * t277 + t326 * t418 - t440 * t70 + (qJD(5) * t294 - t823) * t571 + t686;
t14 = -t129 * t303 + t150 * t559 - t157 * t246 + t170 * t377 - t277 * t474 - t296 * t831 - t325 * t418 + t440 * t69 + t676;
t7 = -t150 * t473 + t151 * t474 + t157 * t169 - t158 * t170 + t294 * t325 - t296 * t326 + t301 * t69 + t303 * t70 + t683;
t4 = [(t979 + t1038 + t981) * t799 - (qJD(1) * t967 - t509 + t836) * t510 * t829 / 0.2e1 + (t972 + t974) * t803 + (t973 + t975) * t802 + (t980 + t1037) * t797 + t31 * t934 + t17 * t945 + (t116 * (t367 + t719) + t178 * (-t232 + t662) + t115 * (t664 + t368) + t179 * (t657 + t231)) * t791 + t901 / 0.2e1 - t902 / 0.2e1 + (t33 * (-t556 * t916 + t689 + t769) + t84 * (t403 * t916 + t654 + t770) + t32 * (t656 + t258) + t85 * (t652 + t112)) * t787 + (t72 * (-t556 * t917 + t711 + t771) + t104 * (t403 * t917 + t658 + t772) + t71 * (t661 + t371) + t105 * (t653 + t233)) * t789 + t897 / 0.2e1 + t80 * t947 + t301 * t948 + (t15 * (-t466 * t915 + t689 + t767 + t905 - t906) + t41 * (t241 * t915 + t654 + t768 - t907 - t908) + t14 * (t656 + t866) + t42 * (t652 + t894)) * t785 + t473 * t1049 + t896 / 0.2e1 + (t1051 - t178 + t815) * t179 * t792 + t36 * t931 + t37 * t933 + (((t1018 - t1105 + t978) * t592 + t956 * t591) * qJD(3) + t992) * t800 + (-t1073 + ((-t1076 - t1106 + t956) * t592 - t1095 * t591) * qJD(3) + t982) * t798 + t892 + t895 / 0.2e1 + (t1077 + t700 - t84) * t85 * t788 + t898 / 0.2e1 + (t480 * t953 + t509 * (-qJ(2) * t746 - qJD(1) * t647 + t641 - t659) + t479 * (t812 + t712) + t510 * (t740 + t813)) * t829 / 0.2e1 + (t1094 - t41 + t700) * t42 * t786 + t99 * t940 + t100 * t941 + t19 * t942 + t20 * t944 + t976 + t899 / 0.2e1 + t900 / 0.2e1 + t903 + (-t370 * t601 - t841 * t997 - t104 + t764) * t105 * t790 - t1074 * t946; (-t479 * t780 + t480 * t779) * m(3) + (-t14 * t780 + t15 * t779 + t7 * t891) * m(7) + (t21 * t891 - t32 * t780 + t33 * t779) * m(6) + (t38 * t891 - t71 * t780 + t72 * t779) * m(5) + (-t115 * t780 + t116 * t779 + t87 * t891) * m(4); m(7) * (t14 * t851 + t15 * t846 + t35 * t782 + t41 * t817 + t42 * t821 + t7 * t856 + (t35 * t783 + t41 * t861) * t576 - (t35 * t867 + t42 * t781) * t575 + (t14 * t781 + t42 * (-t443 - t868) + t7 * t867 + t35 * t893) * t592 + (t15 * (t849 - t867) + t41 * (t862 - t893) + t14 * t866 + t42 * t894) * t614 + (t15 * t861 + t41 * (t868 - t909) + t7 * t783 + t35 * (t864 - t894)) * t591) + (t28 * t591 + t29 * t592 - t575 * t97 + t576 * t96 + t60 * t614) * t927 + (t22 * t591 + t23 * t592 + t36 * t614 - t575 * t91 + t576 * t90) * t931 + (t24 * t591 + t25 * t592 + t37 * t614 - t575 * t89 + t576 * t88) * t933 + (t12 * t591 + t13 * t592 + t27 * t614 - t575 * t62 + t576 * t61) * t935 + (t33 * (t343 * t591 + (-t257 + t849) * t614 + t846) + t84 * (t343 * t576 + (t264 - t909) * t591 + (-t113 + t862) * t614 + t817) + t32 * (t258 * t614 + t592 * t819 + t851) + t85 * (t112 * t614 + (-t264 - t443) * t592 - t819 * t575 + t821) + t21 * (t257 * t592 + t591 * t820 + t856) + t78 * (t113 * t592 - t257 * t575 + (-t112 + t864) * t591 + t820 * t576 + t782)) * t787 + (t19 * t614 + t57 * t576 - t575 * t58 + t591 * t8 + t592 * t9) * t942 + (t10 * t591 + t11 * t592 + t20 * t614 + t55 * t576 - t56 * t575) * t944 + (t178 * (-t433 * t601 + t520 * t997) + t179 * (t438 * t601 - t520 * t584) + (t433 * t592 - t438 * t591) * qJD(3) * t177) * t792 + (t1016 * t614 + t966) * t601 / 0.2e1 + (t966 * qJD(3) + t976 + t985) * t614 / 0.2e1 + (t964 * qJD(3) + t979 * t601 + t986) * t591 / 0.2e1 + (t104 * (t518 * t997 + (-t431 - t432) * t601 + t847) + t105 * (t436 * t601 + (-t518 - t519) * t584 + t852) + t101 * ((t431 * t592 + (-t436 - t437) * t591) * qJD(3) + t859)) * t790 + t832 * t1049 - ((t587 * t958 + t588 * t957 + t614 * t968) * t601 + (t587 * t963 + t588 * t962 + t614 * t961) * qJD(3)) * t601 / 0.2e1 + (t614 * t979 + t964) * t799 + (t614 * t980 + t965) * t797 - (t981 + t983) * t575 / 0.2e1 + (t982 + t984) * t576 / 0.2e1 + (t116 * (t367 * t614 + t489 * t591) + t178 * (-t232 * t614 + t751) + t115 * (t368 * t614 - t489 * t592) + t179 * (t231 * t614 + t752) + t87 * t755 + t177 * t702) * t791 + ((-t553 * t249 + t318 * t466 + t469 * t320) * t474 + (-t553 * t248 + t317 * t466 + t469 * t319) * t473 + (-t553 * t340 + t445 * t466 + t469 * t446) * t559 + (-t553 * t88 - t557 * t89 - t587 * t99) * qJD(5) - t655 * t556) * t934 + ((t553 * t958 - t556 * t957 + t591 * t968) * t601 + (t553 * t963 - t556 * t962 + t591 * t961) * qJD(3)) * t800 + (t57 * t591 + t58 * t592 + t614 * t80) * t947 + t412 * t948 + (t41 * (-t169 * t508 - t204 * t440 + t246 * t411 - t300 * t301 + t409 * t559 + t473 * t507 + (t294 * t587 - t418 * t553) * qJD(5) + t703) + t42 * (t170 * t508 + t205 * t440 - t246 * t412 - t300 * t303 + t410 * t559 - t474 * t507 + (-t296 * t587 + t418 * t557) * qJD(5) + t704) + t35 * (t169 * t412 - t170 * t411 + t204 * t303 + t205 * t301 - t409 * t474 - t410 * t473 + (-t294 * t557 + t296 * t553) * qJD(5) + t701)) * t786 + ((-t161 * t809 + t505 * t164 + t506 * t167 + t551 * t199 + t458 * t201 + t459 * t203) * t303 + t62 * t412 - (-t160 * t809 + t505 * t163 + t506 * t166 + t551 * t198 + t458 * t200 + t459 * t202) * t301 + t61 * t411 + (-t243 * t809 + t505 * t244 + t506 * t245 + t551 * t297 + t458 * t298 + t459 * t299) * t440 + t86 * t508) * t936 + ((-t161 * t810 + t406 * t164 + t407 * t167 + t470 * t199 + t337 * t201 + t338 * t203) * t303 + t58 * t412 - (-t160 * t810 + t406 * t163 + t407 * t166 + t470 * t198 + t337 * t200 + t338 * t202) * t301 + t57 * t411 + (-t243 * t810 + t406 * t244 + t407 * t245 + t470 * t297 + t337 * t298 + t338 * t299) * t440 + t80 * t508) * t943 + ((t557 * t958 + t558 * t957 + t592 * t968) * t601 + (t557 * t963 + t558 * t962 + t592 * t961) * qJD(3)) * t798 + ((-t587 * t249 - t551 * t318 + t552 * t320) * t474 + (-t587 * t248 - t551 * t317 + t552 * t319) * t473 + (-t587 * t340 - t551 * t445 + t552 * t446) * t559 + (-t117 * t587 - t553 * t96 - t557 * t97) * qJD(5) + t655 * t588) * t928 + ((-t557 * t249 - t470 * t318 + t471 * t320) * t474 + (-t557 * t248 - t470 * t317 + t471 * t319) * t473 + (-t557 * t340 - t470 * t445 + t471 * t446) * t559 + (-t100 * t587 - t553 * t90 - t557 * t91) * qJD(5) + t655 * t558) * t932 + (t614 * t973 + t998) * t802 + (t614 * t972 + t999) * t803 - t411 * t16 / 0.2e1 + t34 * t830 / 0.2e1 + t30 * t833 / 0.2e1 + (t591 * t61 + t592 * t62 + t614 * t86) * t939 + (t591 * t88 + t592 * t89 + t614 * t99) * t940 + (t100 * t614 + t591 * t90 + t592 * t91) * t941 + (t117 * t614 + t591 * t96 + t592 * t97) * t795 - t508 * t18 / 0.2e1 + (t72 * (t490 * t591 + t614 * t858 + t454) + t104 * (t614 * t863 + t750 + t853) + t71 * (t371 * t614 + t592 * t841 + t381) + t105 * (t233 * t614 + t189 + t721) + t38 * (t372 + t739) + t101 * (t687 + t865)) * t789 + (t965 * qJD(3) + t980 * t601 + t987) * t592 / 0.2e1 + (t84 * (-t322 * t559 + t448 * t473 + (t257 * t587 - t343 * t553) * qJD(5) + t703) + t85 * (t323 * t559 - t448 * t474 + (-t258 * t587 + t343 * t557) * qJD(5) + t704) + t78 * (t322 * t474 - t323 * t473 + (-t257 * t557 + t258 * t553) * qJD(5) + t701)) * t788 + ((t161 * t811 + t404 * t164 + t405 * t167 - t199 * t466 - t201 * t333 + t336 * t203) * t303 + t56 * t412 - (t160 * t811 + t404 * t163 + t405 * t166 - t198 * t466 - t200 * t333 + t336 * t202) * t301 + t55 * t411 + (t243 * t811 + t404 * t244 + t405 * t245 - t297 * t466 - t298 * t333 + t336 * t299) * t440 - t1074 * t508) * t945 + (-t1074 * t614 + t55 * t591 + t56 * t592) * t946; (t14 * t553 + t15 * t557 + t35 * t860 + t41 * t855 + t42 * t854 + t587 * t7) * m(7) + (t21 * t587 + t32 * t553 + t33 * t557 + t78 * t860 + t84 * t855 + t85 * t854) * m(6) + (t101 * t860 + t104 * t855 + t105 * t854 + t38 * t587 + t553 * t71 + t557 * t72) * m(5); -(t16 * t469 + t17 * t471 + t18 * t552) * qJD(6) / 0.2e1 - (t34 + t18) * t571 / 0.2e1 + t985 * t588 / 0.2e1 + t983 * t401 / 0.2e1 + t984 * t403 / 0.2e1 + (-t466 * t694 + t469 * t695 - t556 * t726) * t934 - t986 * t556 / 0.2e1 + (t33 * (-t257 * t588 - t343 * t556) + t84 * (-t113 * t588 + t257 * t571 - t264 * t556 + t343 * t403) + t32 * (t258 * t588 - t343 * t558) + t85 * (t112 * t588 - t258 * t571 - t264 * t558 - t343 * t401) + t21 * (t257 * t558 + t258 * t556) + t78 * (t112 * t556 + t113 * t558 + t257 * t401 - t258 * t403)) * t787 + (-t556 * t57 + t558 * t58 + t588 * t80) * t947 + (-t117 * t571 - t28 * t556 + t29 * t558 + t401 * t97 + t403 * t96 + t588 * t60) * t927 + (-t100 * t571 - t22 * t556 + t23 * t558 + t36 * t588 + t401 * t91 + t403 * t90) * t931 + (-t24 * t556 + t25 * t558 + t37 * t588 + t401 * t89 + t403 * t88 - t571 * t99) * t933 + (-t12 * t556 + t13 * t558 + t27 * t588 + t401 * t62 + t403 * t61 - t571 * t86) * t935 + (-t556 * t61 + t558 * t62 + t588 * t86) * t939 + (-t556 * t88 + t558 * t89 + t588 * t99) * t940 + (t100 * t588 - t556 * t90 + t558 * t91) * t941 + (t19 * t588 + t401 * t58 + t403 * t57 - t556 * t8 + t558 * t9 - t571 * t80) * t942 + (t117 * t588 - t556 * t96 + t558 * t97) * t795 + (t15 * (-t556 * t861 - t588 * t867) + t41 * (t403 * t861 - t556 * t868 + t571 * t867 - t588 * t893) + t14 * (-t558 * t861 + t588 * t866) + t42 * (-t401 * t861 - t558 * t868 - t571 * t866 + t588 * t894) + t7 * (t556 * t866 + t558 * t867) + t35 * (t401 * t867 - t403 * t866 + t556 * t894 + t558 * t893)) * t785 + ((t161 * t552 + t214 * t458 + t216 * t459) * t303 - (t160 * t552 + t213 * t458 + t215 * t459) * t301 + (t243 * t552 + t306 * t458 + t307 * t459) * t440 + (t469 * t61 + t471 * t62 + t552 * t86) * qJD(6) + t660 * t551) * t936 + ((t161 * t471 + t214 * t337 + t216 * t338) * t303 - (t160 * t471 + t213 * t337 + t215 * t338) * t301 + (t243 * t471 + t306 * t337 + t307 * t338) * t440 + (t469 * t57 + t471 * t58 + t552 * t80) * qJD(6) + t660 * t470) * t943 + (t41 * (-t217 * t440 - t293 * t559 - t301 * t308 + t417 * t473 + (-t169 * t552 + t246 * t469) * qJD(6)) + t42 * (t218 * t440 + t295 * t559 - t303 * t308 - t417 * t474 + (t170 * t552 - t246 * t471) * qJD(6)) + t35 * (t217 * t303 + t218 * t301 + t293 * t474 - t295 * t473 + (t169 * t471 - t170 * t469) * qJD(6))) * t786 + (t551 * t694 + t552 * t695 + t588 * t726) * t928 + (t470 * t694 + t471 * t695 + t558 * t726) * t932 + t987 * t558 / 0.2e1 + (t84 * (-t290 * t559 + t416 * t473) + t85 * (t291 * t559 - t416 * t474) + t78 * (t290 * t474 - t291 * t473)) * t788 + (-t10 * t556 + t1074 * t571 + t11 * t558 + t20 * t588 + t401 * t56 + t403 * t55) * t944 + (-t1074 * t588 - t55 * t556 + t558 * t56) * t946 + ((t161 * t469 - t214 * t333 + t216 * t336) * t303 - (t160 * t469 - t213 * t333 + t215 * t336) * t301 + (t243 * t469 - t306 * t333 + t307 * t336) * t440 + (-t1074 * t552 + t469 * t55 + t471 * t56) * qJD(6) - t660 * t466) * t945; t239 * t17 / 0.2e1 + t470 * t1 / 0.2e1 + (-t466 * t57 + t470 * t58 + t551 * t80) * t947 + (t19 * t551 + t239 * t58 + t241 * t57 + t456 * t80 - t466 * t8 + t470 * t9) * t942 + t241 * t16 / 0.2e1 - t466 * t2 / 0.2e1 + (-t1074 * t551 - t466 * t55 + t470 * t56) * t946 + (-t10 * t466 - t1074 * t456 + t11 * t470 + t20 * t551 + t239 * t56 + t241 * t55) * t944 + t456 * t18 / 0.2e1 + t551 * t3 / 0.2e1 + (-t466 * t61 + t470 * t62 + t551 * t86) * t939 + (-t12 * t466 + t13 * t470 + t239 * t62 + t241 * t61 + t27 * t551 + t456 * t86) * t935 + (t15 * (-t169 * t551 - t246 * t466) + t41 * (-t129 * t466 - t169 * t456 + t241 * t246 - t551 * t70) + t14 * (t170 * t551 - t246 * t470) + t42 * (-t129 * t470 + t170 * t456 - t239 * t246 + t551 * t69) + t7 * (t169 * t470 + t170 * t466) + t35 * (t169 * t239 - t170 * t241 + t466 * t69 + t470 * t70)) * t785 + (t337 * t696 + t338 * t697 + t470 * t727) * t943 + (-t333 * t696 + t336 * t697 - t466 * t727) * t945 + (t458 * t696 + t459 * t697 + t551 * t727) * t936 + (t41 * (-t196 * t440 - t281 * t301) + t42 * (t197 * t440 - t281 * t303) + t35 * (t196 * t303 + t197 * t301)) * t786;];
tauc  = t4(:);
