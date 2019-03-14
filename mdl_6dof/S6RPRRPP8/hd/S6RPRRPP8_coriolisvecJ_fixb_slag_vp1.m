% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPP8_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP8_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:53:54
% EndTime: 2019-03-09 04:55:47
% DurationCPUTime: 104.81s
% Computational Cost: add. (21699->1231), mult. (57458->1547), div. (0->0), fcn. (54881->6), ass. (0->568)
t1126 = -Icges(6,4) + Icges(5,5) + Icges(7,5);
t1051 = Icges(7,4) + Icges(6,5) - Icges(5,6);
t557 = sin(qJ(4));
t1130 = (-Icges(5,4) + Icges(7,6)) * t557;
t560 = cos(qJ(4));
t1129 = (-Icges(5,4) - Icges(6,6)) * t560;
t1128 = t1130 + (Icges(5,1) + Icges(7,3)) * t560;
t1127 = -t1129 + (-Icges(5,2) - Icges(6,3)) * t557;
t558 = sin(qJ(3));
t561 = cos(qJ(3));
t1125 = t1126 * t558 + t1128 * t561;
t1124 = -t1051 * t558 + t1127 * t561;
t864 = t560 * t561;
t529 = Icges(7,6) * t864;
t872 = t557 * t561;
t1099 = Icges(7,2) * t872 - t1124 + t529;
t530 = Icges(6,6) * t872;
t1092 = -Icges(6,2) * t864 - t1125 + t530;
t1121 = Icges(5,3) + Icges(6,1) + Icges(7,1);
t1120 = t1051 * t557 + t1126 * t560;
t1119 = -Icges(5,1) - Icges(6,2);
t1117 = -Icges(5,2) - Icges(7,2);
t559 = sin(qJ(1));
t866 = t559 * t560;
t562 = cos(qJ(1));
t870 = t558 * t562;
t438 = t557 * t870 + t866;
t416 = Icges(7,6) * t438;
t862 = t562 * t560;
t867 = t559 * t557;
t439 = t558 * t862 - t867;
t863 = t561 * t562;
t206 = Icges(7,5) * t863 - Icges(7,3) * t439 - t416;
t419 = Icges(6,6) * t438;
t216 = -Icges(6,4) * t863 - Icges(6,2) * t439 + t419;
t422 = Icges(5,4) * t438;
t230 = Icges(5,1) * t439 - Icges(5,5) * t863 - t422;
t1011 = -t230 + t216 + t206;
t420 = Icges(6,6) * t439;
t209 = Icges(6,5) * t863 - Icges(6,3) * t438 + t420;
t417 = Icges(7,6) * t439;
t213 = -Icges(7,4) * t863 + Icges(7,2) * t438 + t417;
t423 = Icges(5,4) * t439;
t228 = Icges(5,2) * t438 + Icges(5,6) * t863 - t423;
t1012 = -t228 + t209 - t213;
t218 = Icges(7,1) * t863 - Icges(7,4) * t438 - Icges(7,5) * t439;
t221 = Icges(6,1) * t863 + Icges(6,4) * t439 - Icges(6,5) * t438;
t225 = -Icges(5,5) * t439 + Icges(5,6) * t438 + Icges(5,3) * t863;
t1013 = t218 + t221 + t225;
t1064 = -t1011 * t439 - t1012 * t438 + t1013 * t863;
t1115 = -t1120 * t558 + t1121 * t561;
t1108 = t1120 * t561 + t1121 * t558;
t1098 = t1092 * t560 - t1099 * t557;
t436 = t558 * t867 - t862;
t437 = t557 * t562 + t558 * t866;
t865 = t559 * t561;
t217 = -Icges(7,1) * t865 + Icges(7,4) * t436 + Icges(7,5) * t437;
t220 = -Icges(6,1) * t865 - Icges(6,4) * t437 + Icges(6,5) * t436;
t223 = Icges(5,5) * t437 - Icges(5,6) * t436 - Icges(5,3) * t865;
t1083 = t217 + t220 + t223;
t891 = Icges(6,6) * t437;
t208 = -Icges(6,5) * t865 + Icges(6,3) * t436 - t891;
t415 = Icges(7,6) * t437;
t211 = -Icges(7,4) * t865 + Icges(7,2) * t436 + t415;
t901 = Icges(5,4) * t437;
t226 = -Icges(5,2) * t436 - Icges(5,6) * t865 + t901;
t1084 = t208 + t211 - t226;
t888 = Icges(7,6) * t436;
t205 = -Icges(7,5) * t865 + Icges(7,3) * t437 + t888;
t418 = Icges(6,6) * t436;
t214 = -Icges(6,4) * t865 - Icges(6,2) * t437 + t418;
t421 = Icges(5,4) * t436;
t229 = Icges(5,1) * t437 - Icges(5,5) * t865 - t421;
t1085 = t205 - t214 + t229;
t1017 = -t1083 * t865 + t1084 * t436 + t1085 * t437;
t1114 = t1017 - t1064;
t1015 = t1083 * t863 - t1084 * t438 - t1085 * t439;
t787 = qJD(4) * t561;
t789 = qJD(3) * t562;
t478 = -t559 * t787 + t789;
t788 = qJD(4) * t558;
t534 = qJD(1) + t788;
t978 = t1092 * t439 - t1099 * t438 + t1108 * t863;
t1113 = -t1015 * t478 - t978 * t534;
t1111 = (t1117 * t560 + t1130) * t561;
t1110 = (t1119 * t557 + t1129) * t561;
t886 = Icges(7,6) * t560;
t665 = Icges(7,2) * t557 + t886;
t1006 = -t1051 * t561 + (t665 - t1127) * t558;
t890 = Icges(6,6) * t557;
t666 = Icges(6,2) * t560 - t890;
t1005 = t1126 * t561 + (-t666 - t1128) * t558;
t1107 = (t1051 * t560 - t1126 * t557) * t561;
t1065 = t1011 * t437 + t1012 * t436 - t1013 * t865;
t786 = qJD(4) * t562;
t791 = qJD(3) * t559;
t477 = t561 * t786 + t791;
t653 = t228 * t557 + t230 * t560;
t655 = t226 * t557 - t229 * t560;
t656 = t209 * t557 + t216 * t560;
t658 = t208 * t557 - t214 * t560;
t659 = t206 * t560 - t213 * t557;
t661 = t205 * t560 + t211 * t557;
t1101 = (-t1098 - t1115) * t534 + (-t1108 * t559 - t655 + t658 + t661) * t478 + (t1108 * t562 - t653 + t656 + t659) * t477;
t752 = t558 * t789;
t794 = qJD(1) * t561;
t1036 = t559 * t794 + t752;
t721 = qJD(1) * t558 + qJD(4);
t748 = t558 * t786;
t750 = t561 * t789;
t793 = qJD(1) * t562;
t198 = -t557 * t750 + t721 * t867 + (-t748 - t793) * t560;
t199 = qJD(1) * t437 + qJD(4) * t438 - t560 * t750;
t102 = -Icges(7,5) * t1036 + Icges(7,6) * t198 + Icges(7,3) * t199;
t108 = -Icges(6,4) * t1036 - Icges(6,2) * t199 + Icges(6,6) * t198;
t118 = Icges(5,1) * t199 - Icges(5,4) * t198 - Icges(5,5) * t1036;
t1091 = t102 - t108 + t118;
t790 = qJD(3) * t561;
t751 = t559 * t790;
t587 = t562 * t721 + t751;
t200 = t534 * t866 + t557 * t587;
t749 = t559 * t788;
t795 = qJD(1) * t559;
t201 = t560 * t587 + (-t749 - t795) * t557;
t753 = t558 * t791;
t755 = t561 * t793;
t612 = t753 - t755;
t103 = Icges(7,5) * t612 + Icges(7,6) * t200 + Icges(7,3) * t201;
t109 = Icges(6,4) * t612 - Icges(6,2) * t201 + Icges(6,6) * t200;
t119 = Icges(5,1) * t201 - Icges(5,4) * t200 + Icges(5,5) * t612;
t1090 = t103 - t109 + t119;
t104 = -Icges(6,5) * t1036 - Icges(6,6) * t199 + Icges(6,3) * t198;
t106 = -Icges(7,4) * t1036 + Icges(7,2) * t198 + Icges(7,6) * t199;
t116 = Icges(5,4) * t199 - Icges(5,2) * t198 - Icges(5,6) * t1036;
t1089 = t104 + t106 - t116;
t105 = Icges(6,5) * t612 - Icges(6,6) * t201 + Icges(6,3) * t200;
t107 = Icges(7,4) * t612 + Icges(7,2) * t200 + Icges(7,6) * t201;
t117 = Icges(5,4) * t201 - Icges(5,2) * t200 + Icges(5,6) * t612;
t1088 = t105 + t107 - t117;
t110 = -Icges(7,1) * t1036 + Icges(7,4) * t198 + Icges(7,5) * t199;
t112 = -Icges(6,1) * t1036 - Icges(6,4) * t199 + Icges(6,5) * t198;
t114 = Icges(5,5) * t199 - Icges(5,6) * t198 - Icges(5,3) * t1036;
t1087 = t110 + t112 + t114;
t111 = Icges(7,1) * t612 + Icges(7,4) * t200 + Icges(7,5) * t201;
t113 = Icges(6,1) * t612 - Icges(6,4) * t201 + Icges(6,5) * t200;
t115 = Icges(5,5) * t201 - Icges(5,6) * t200 + Icges(5,3) * t612;
t1086 = t111 + t113 + t115;
t1082 = qJD(3) * t1115 + t1107 * qJD(4);
t1081 = -(Icges(6,3) * t560 + t890) * t787 + t1111 * qJD(4) + t1006 * qJD(3);
t1080 = (-Icges(7,3) * t557 + t886) * t787 + t1110 * qJD(4) + t1005 * qJD(3);
t1074 = Icges(7,3) - t1119;
t1073 = Icges(6,3) - t1117;
t979 = -t1092 * t437 + t1099 * t436 - t1108 * t865;
t1069 = t1101 * t561;
t944 = rSges(7,1) + pkin(5);
t1068 = t562 * pkin(1);
t799 = t559 * qJ(2) + t1068;
t489 = qJD(1) * t799;
t546 = qJD(2) * t562;
t424 = -t546 + t489;
t725 = -rSges(3,2) * t562 + t559 * rSges(3,3);
t319 = qJD(1) * t725 + t424;
t1026 = -t1036 * t1083 + t1084 * t198 + t1085 * t199 + t1086 * t863 - t1088 * t438 - t1090 * t439;
t1025 = t1011 * t199 + t1012 * t198 - t1013 * t1036 + t1087 * t863 - t1089 * t438 - t1091 * t439;
t1024 = t1083 * t612 + t1084 * t200 + t1085 * t201 - t1086 * t865 + t1088 * t436 + t1090 * t437;
t1023 = t1011 * t201 + t1012 * t200 + t1013 * t612 - t1087 * t865 + t1089 * t436 + t1091 * t437;
t1019 = -t1036 * t1108 - t1080 * t439 + t1081 * t438 + t1082 * t863 - t1092 * t199 + t1099 * t198;
t1018 = t1080 * t437 - t1081 * t436 - t1082 * t865 - t1092 * t201 + t1099 * t200 + t1108 * t612;
t84 = t223 * t558 - t561 * t655;
t86 = t217 * t558 + t561 * t661;
t88 = t220 * t558 + t561 * t658;
t1063 = t84 + t86 + t88;
t85 = t225 * t558 - t561 * t653;
t87 = t218 * t558 + t561 * t659;
t89 = t221 * t558 + t561 * t656;
t1062 = t85 + t87 + t89;
t977 = -t1098 * t561 + t1108 * t558;
t1059 = -t561 * t665 + t1124;
t1058 = t561 * t666 + t1125;
t1057 = (t1074 * t436 - t1084 - t415 + t891 + t901) * t478 + (-t1074 * t438 - t1012 + t417 - t420 - t423) * t477 - (-Icges(7,3) * t872 + t1099 + t1110 + t529) * t534;
t1056 = (-t1073 * t437 + t1085 - t418 - t421 + t888) * t478 + (t1073 * t439 + t1011 - t416 + t419 + t422) * t477 - (Icges(6,3) * t864 + t1092 - t1111 + t530) * t534;
t1041 = rSges(7,3) + qJ(6);
t871 = t558 * t559;
t535 = pkin(3) * t871;
t761 = -pkin(3) * t750 - pkin(8) * t1036;
t312 = qJD(1) * t535 + t761;
t509 = pkin(3) * t561 + pkin(8) * t558;
t469 = t509 * t562;
t868 = t559 * t509;
t1055 = t562 * t312 + t469 * t789 + t791 * t868;
t1053 = (t1080 * t560 - t1081 * t557 + (t1092 * t557 + t1099 * t560) * qJD(4) + t1108 * qJD(3)) * t561 + (qJD(3) * t1098 + t1082) * t558;
t1043 = t1065 * t477 + t534 * t979;
t1042 = 0.2e1 * qJD(3);
t251 = -rSges(5,1) * t439 + rSges(5,2) * t438 + rSges(5,3) * t863;
t1040 = t251 * t534;
t698 = rSges(7,2) * t557 + rSges(7,3) * t560;
t814 = qJ(6) * t864 + t558 * t944 + t561 * t698;
t784 = qJD(5) * t561;
t525 = t557 * t784;
t792 = qJD(3) * t558;
t609 = -t557 * t792 + t560 * t787;
t610 = -t557 * t787 - t560 * t792;
t242 = pkin(4) * t610 + qJ(5) * t609 + t525;
t697 = pkin(4) * t560 + qJ(5) * t557;
t463 = t697 * t561;
t391 = t463 * t793;
t458 = t697 * t558;
t941 = pkin(3) * t558;
t508 = pkin(8) * t561 - t941;
t724 = qJD(1) * t469 + t508 * t791;
t487 = qJD(3) * t508;
t806 = t559 * t487 + t509 * t793;
t1038 = t477 * t458 + t391 - t724 + t806 + (t242 - t525) * t559;
t1037 = t558 * t793 + t751;
t304 = t437 * pkin(4) + qJ(5) * t436;
t466 = -pkin(8) * t865 + t535;
t801 = pkin(7) * t793 - t546;
t764 = t489 + t801;
t677 = qJD(1) * t466 - t509 * t789 + t764;
t780 = t438 * qJD(5);
t1035 = -t534 * t304 + t478 * t463 - t677 + t780;
t474 = t509 * t795;
t869 = t559 * t463;
t1034 = -t304 * t787 + t463 * t795 - t534 * t869 + t474;
t852 = t198 * qJ(5) - t780;
t940 = pkin(4) * t199;
t120 = t852 + t940;
t414 = t438 * qJ(5);
t309 = pkin(4) * t439 + t414;
t365 = t697 * t863;
t785 = qJD(5) * t557;
t746 = t558 * t785;
t1033 = t562 * t120 - t304 * t748 + t309 * t749 + t478 * t365 + t1055 + t746;
t408 = t439 * qJD(6);
t1032 = t478 * t814 + t1035 + t408;
t531 = Icges(4,4) * t865;
t896 = Icges(4,5) * t562;
t376 = Icges(4,1) * t871 + t531 + t896;
t902 = Icges(4,4) * t561;
t674 = Icges(4,1) * t558 + t902;
t377 = -Icges(4,5) * t559 + t562 * t674;
t498 = -Icges(4,2) * t558 + t902;
t454 = t498 * t562;
t595 = t559 * (t377 + t454) - t562 * (-Icges(4,2) * t871 + t376 + t531);
t903 = Icges(4,4) * t558;
t672 = Icges(4,2) * t561 + t903;
t372 = Icges(4,6) * t562 + t559 * t672;
t373 = -Icges(4,6) * t559 + t562 * t672;
t500 = Icges(4,1) * t561 - t903;
t456 = t500 * t559;
t457 = t500 * t562;
t596 = t559 * (t373 - t457) - t562 * (t372 - t456);
t1031 = -t596 * t558 + t595 * t561;
t804 = t498 + t674;
t805 = -t672 + t500;
t1029 = (t558 * t804 - t561 * t805) * qJD(1);
t777 = qJD(3) * qJD(4);
t740 = t558 * t777;
t356 = -qJD(1) * t477 + t559 * t740;
t357 = qJD(1) * t478 - t562 * t740;
t739 = t561 * t777;
t1028 = t1015 * t356 + t1019 * t534 + t1025 * t477 + t1026 * t478 + t1064 * t357 + t739 * t978;
t1027 = t1017 * t356 + t1018 * t534 + t1023 * t477 + t1024 * t478 + t1065 * t357 + t739 * t979;
t25 = (qJD(3) * t655 + t115) * t558 + (qJD(3) * t223 - t117 * t557 + t119 * t560 + (-t226 * t560 - t229 * t557) * qJD(4)) * t561;
t27 = (-qJD(3) * t661 + t111) * t558 + (qJD(3) * t217 + t103 * t560 + t107 * t557 + (-t205 * t557 + t211 * t560) * qJD(4)) * t561;
t29 = (-qJD(3) * t658 + t113) * t558 + (qJD(3) * t220 + t105 * t557 - t109 * t560 + (t208 * t560 + t214 * t557) * qJD(4)) * t561;
t1022 = t25 + t27 + t29;
t26 = (qJD(3) * t653 + t114) * t558 + (qJD(3) * t225 - t116 * t557 + t118 * t560 + (-t228 * t560 + t230 * t557) * qJD(4)) * t561;
t28 = (-qJD(3) * t659 + t110) * t558 + (qJD(3) * t218 + t102 * t560 + t106 * t557 + (-t206 * t557 - t213 * t560) * qJD(4)) * t561;
t30 = (-qJD(3) * t656 + t112) * t558 + (qJD(3) * t221 + t104 * t557 - t108 * t560 + (t209 * t560 - t216 * t557) * qJD(4)) * t561;
t1021 = t26 + t28 + t30;
t981 = t1017 * t478 + t1043;
t980 = t1064 * t477 - t1113;
t1020 = t1062 * t477 + t1063 * t478 + t534 * t977;
t1010 = t1059 * t559;
t1009 = t1059 * t562;
t1008 = t1058 * t559;
t1007 = t1058 * t562;
t1004 = -t1013 * t477 - t1083 * t478 - t1108 * t534;
t972 = -t1062 * t562 + t1063 * t559;
t1003 = t1062 * t559 + t1063 * t562;
t971 = t1015 * t559 - t1064 * t562;
t970 = t1017 * t559 - t1065 * t562;
t1002 = t1015 * t562 + t1064 * t559;
t1001 = t1017 * t562 + t1065 * t559;
t982 = t436 * rSges(7,2) + t1041 * t437;
t837 = -t865 * t944 + t982;
t53 = t837 * t534 - t1032;
t765 = t463 + t814;
t999 = t53 * t765;
t699 = rSges(6,2) * t560 - rSges(6,3) * t557;
t964 = -rSges(6,1) * t558 + t561 * t699;
t813 = -t964 + t463;
t809 = -t437 * rSges(6,2) + t436 * rSges(6,3);
t244 = -rSges(6,1) * t865 + t809;
t983 = -t244 * t534 - t478 * t964 + t1035;
t998 = t983 * t813;
t648 = t373 * t561 + t377 * t558;
t993 = t648 * t562;
t645 = t436 * t477 + t438 * t478;
t989 = t645 + t609;
t460 = (-pkin(4) * t557 + qJ(5) * t560) * t561;
t988 = -qJD(5) * t437 + t242 * t863 - t477 * t460;
t632 = -t436 * t534 + t478 * t872;
t987 = t198 + t632;
t930 = rSges(7,2) * t438;
t835 = -t1041 * t439 + t863 * t944 - t930;
t985 = -rSges(7,2) * t198 + t408;
t548 = t562 * qJ(2);
t502 = pkin(1) * t559 - t548;
t937 = pkin(7) * t559;
t728 = -t502 - t937;
t409 = qJD(6) * t437;
t984 = t200 * rSges(7,2) + t1041 * t201 + t753 * t944 + t409;
t973 = t1107 * t534 + (t1051 * t437 - t1126 * t436) * t478 + (-t1051 * t439 + t1126 * t438) * t477;
t969 = t1053 * t534 + t739 * t977;
t701 = rSges(5,1) * t560 - rSges(5,2) * t557;
t394 = rSges(5,3) * t558 + t561 * t701;
t197 = t373 * t558 - t377 * t561;
t256 = qJD(1) * t372 - qJD(3) * t454;
t259 = -qJD(3) * t457 + (t559 * t674 + t896) * qJD(1);
t668 = Icges(4,5) * t558 + Icges(4,6) * t561;
t369 = -Icges(4,3) * t559 + t562 * t668;
t797 = qJD(1) * t369;
t963 = qJD(3) * t197 + t256 * t561 + t259 * t558 + t797;
t480 = t672 * qJD(3);
t481 = t674 * qJD(3);
t496 = Icges(4,5) * t561 - Icges(4,6) * t558;
t962 = qJD(1) * t496 + qJD(3) * (t498 * t558 - t500 * t561) + t480 * t561 + t481 * t558;
t257 = qJD(1) * t373 + t498 * t791;
t260 = qJD(1) * t377 + qJD(3) * t456;
t649 = t372 * t558 - t376 * t561;
t368 = Icges(4,3) * t562 + t559 * t668;
t798 = qJD(1) * t368;
t961 = qJD(3) * t649 - t257 * t561 - t260 * t558 + t798;
t957 = -pkin(1) - pkin(7);
t956 = t356 / 0.2e1;
t955 = t357 / 0.2e1;
t954 = -t477 / 0.2e1;
t953 = t477 / 0.2e1;
t952 = -t478 / 0.2e1;
t951 = t478 / 0.2e1;
t950 = -t534 / 0.2e1;
t949 = t534 / 0.2e1;
t947 = t559 / 0.2e1;
t946 = -t562 / 0.2e1;
t943 = rSges(3,2) - pkin(1);
t942 = -rSges(4,3) - pkin(1);
t938 = pkin(5) * t561;
t555 = t562 * pkin(7);
t935 = rSges(6,1) * t561;
t933 = rSges(7,1) * t561;
t929 = rSges(3,3) * t562;
t927 = rSges(5,3) * t561;
t232 = t438 * t534 + t477 * t872;
t410 = qJD(5) * t436;
t536 = pkin(8) * t863;
t468 = pkin(3) * t870 - t536;
t472 = t509 * t791;
t545 = qJD(2) * t559;
t591 = t472 + t545 + (t468 + t728) * qJD(1);
t584 = t477 * t463 + t410 + t591;
t766 = t309 - t835;
t52 = t477 * t814 + t534 * t766 + t409 + t584;
t926 = t232 * t52;
t925 = t25 * t478;
t924 = t26 * t477;
t923 = t27 * t478;
t922 = t28 * t477;
t921 = t29 * t478;
t920 = t30 * t477;
t554 = t562 * rSges(4,3);
t916 = t983 * t964;
t915 = t84 * t356;
t914 = t85 * t357;
t913 = t86 * t356;
t912 = t87 * t357;
t911 = t88 * t356;
t910 = t89 * t357;
t807 = t437 * rSges(5,1) - t436 * rSges(5,2);
t249 = -rSges(5,3) * t865 + t807;
t97 = t249 * t534 - t394 * t478 + t677;
t909 = t97 * t394;
t763 = rSges(4,1) * t1037 + rSges(4,2) * t755;
t272 = (-rSges(4,2) * t792 - rSges(4,3) * qJD(1)) * t559 + t763;
t505 = rSges(4,1) * t561 - rSges(4,2) * t558;
t470 = t505 * t791;
t703 = rSges(4,1) * t558 + rSges(4,2) * t561;
t483 = t703 * qJD(3);
t563 = qJD(1) ^ 2;
t779 = qJD(1) * qJD(2);
t803 = qJ(2) * t793 + t545;
t811 = qJD(1) * (-pkin(1) * t795 + t803) + t559 * t779;
t641 = -t563 * t937 + t811;
t128 = t483 * t789 + (t272 + t470) * qJD(1) + t641;
t883 = t128 * t562;
t465 = t505 * t562;
t271 = -qJD(3) * t465 + (t559 * t703 + t554) * qJD(1);
t538 = t562 * t779;
t710 = -t555 * t563 + t538;
t754 = t505 * t789;
t129 = -t483 * t791 + (-t271 - t424 + t754) * qJD(1) + t710;
t882 = t129 * t559;
t552 = t559 * rSges(4,3);
t395 = t562 * t703 - t552;
t173 = t470 + t545 + (t395 + t728) * qJD(1);
t878 = t173 * t562;
t874 = t487 * t562;
t873 = t496 * t562;
t448 = t559 * t496;
t700 = -rSges(6,2) * t199 + rSges(6,3) * t198;
t123 = -rSges(6,1) * t1036 + t700;
t860 = -t120 - t123;
t121 = t201 * pkin(4) + qJ(5) * t200 + t410;
t773 = rSges(6,1) * t753 - t201 * rSges(6,2) + t200 * rSges(6,3);
t125 = -rSges(6,1) * t755 + t773;
t859 = -t121 - t125;
t858 = -t1036 * t944 + t1041 * t199 - t985;
t857 = -t755 * t944 + t984;
t306 = pkin(4) * t438 - qJ(5) * t439;
t526 = t560 * t784;
t853 = t478 * t306 + t526;
t302 = -pkin(4) * t436 + qJ(5) * t437;
t839 = -qJD(5) * t439 + t534 * t302;
t399 = t558 * t699 + t935;
t461 = (rSges(6,2) * t557 + rSges(6,3) * t560) * t561;
t274 = qJD(3) * t399 + qJD(4) * t461;
t838 = -t242 - t274;
t836 = -t244 - t304;
t247 = rSges(6,1) * t863 + rSges(6,2) * t439 - rSges(6,3) * t438;
t833 = -t247 + t309;
t832 = -t249 - t466;
t398 = -t558 * t698 + t933;
t464 = (rSges(7,2) * t560 - rSges(7,3) * t557) * t561;
t782 = qJD(6) * t561;
t524 = t560 * t782;
t831 = pkin(5) * t790 + qJ(6) * t610 + qJD(3) * t398 + qJD(4) * t464 + t524;
t830 = t558 * t304 + t463 * t865;
t404 = t562 * t468;
t829 = -t309 * t562 - t404;
t827 = t814 * t559;
t826 = t814 * t562;
t812 = -qJ(6) * t558 * t560 + t398 + t938;
t810 = t869 + t868;
t802 = rSges(3,2) * t795 + rSges(3,3) * t793;
t800 = -qJD(1) * t502 + t545;
t796 = qJD(1) * t668;
t783 = qJD(6) * t560;
t644 = t498 * t561 + t500 * t558;
t203 = t562 * t644 - t448;
t781 = t203 * qJD(1);
t778 = qJD(1) * qJD(3);
t775 = -t120 - t858;
t774 = -t121 - t857;
t771 = t201 * rSges(5,1) - t200 * rSges(5,2) + rSges(5,3) * t753;
t770 = t1036 * t304 - t309 * t753;
t769 = -t242 - t831;
t768 = -t304 - t837;
t767 = -t466 + t836;
t144 = t562 * t368 + t372 * t865 + t376 * t871;
t145 = -t562 * t369 - t373 * t865 - t377 * t871;
t762 = pkin(3) * t1037 + pkin(8) * t753;
t393 = rSges(4,1) * t871 + rSges(4,2) * t865 + t554;
t760 = t555 + t799;
t759 = t561 * (-rSges(5,3) - pkin(8));
t758 = (-rSges(6,1) - pkin(8)) * t561;
t745 = t557 * t782;
t741 = t509 * t778;
t735 = -t791 / 0.2e1;
t734 = t791 / 0.2e1;
t733 = t790 / 0.2e1;
t732 = -t789 / 0.2e1;
t729 = -qJ(2) - t941;
t727 = t53 * t814;
t723 = (-pkin(8) - t944) * t561;
t278 = t312 * t789;
t722 = qJD(4) * t526 + t478 * t120 - t309 * t356 + t278;
t720 = t558 * t121 + t242 * t865 + t304 * t790 + t561 * t391;
t718 = -t466 + t768;
t717 = t535 + t760;
t708 = qJD(4) * t733;
t707 = -t466 * t791 - t468 * t789;
t434 = qJD(1) * t868;
t705 = -t508 * t789 + t434;
t702 = rSges(5,1) * t199 - rSges(5,2) * t198;
t583 = t487 * t791 + t562 * t741 + (-t312 - t424) * qJD(1) + t710;
t572 = qJD(5) * t200 + t477 * t242 + t357 * t463 + t583;
t12 = t274 * t477 - t357 * t964 + t534 * t860 + t739 * t833 + t572;
t67 = -t477 * t964 + t534 * t833 + t584;
t696 = -t12 * t964 + t67 * t274;
t676 = t762 + t803;
t675 = -t761 - t801;
t174 = qJD(1) * t393 - t754 + t764;
t662 = t173 * t559 - t174 * t562;
t652 = t249 * t562 + t251 * t559;
t650 = t372 * t561 + t376 * t558;
t642 = t468 + t548;
t640 = -t783 - t785;
t462 = (-rSges(5,1) * t557 - rSges(5,2) * t560) * t561;
t621 = -t244 * t983 + t67 * t833;
t618 = -t309 * t478 + t525 + t707;
t64 = t247 * t478 + t477 * t836 + t618;
t620 = t64 * t244 - t67 * t813;
t617 = (t144 * t562 + t145 * t559) * qJD(3);
t358 = t559 * t368;
t147 = -t562 * t650 + t358;
t148 = -t369 * t559 + t993;
t616 = (t147 * t562 + t148 * t559) * qJD(3);
t176 = (-t393 * t559 - t395 * t562) * qJD(3);
t614 = t675 - t852;
t313 = -pkin(8) * t755 + t762;
t613 = qJD(1) * t313 + t559 * t741 + t641;
t608 = -pkin(7) * t795 + qJD(1) * t468 + t472 + t800;
t607 = t304 + t717;
t9 = qJD(6) * t201 + t357 * t814 + t477 * t831 + t534 * t775 + t739 * t766 + t572;
t600 = t52 * t831 + t814 * t9;
t594 = -qJD(1) * t648 - qJD(3) * t873 + t798;
t593 = qJD(1) * t650 + qJD(3) * t448 + t797;
t592 = qJD(1) * t644 - t668 * qJD(3);
t590 = t52 * t766 + t53 * t837;
t49 = t477 * t768 + t478 * t835 + t524 + t618;
t589 = t49 * t837 - t52 * t765;
t392 = -t558 * t701 + t927;
t588 = t534 * t309 + t410 + t608;
t586 = -t313 * t559 + (-t466 * t562 + t468 * t559) * qJD(1);
t585 = qJD(5) * t198 + t534 * t121 + t304 * t739 + t613;
t582 = t121 + t676;
t83 = -t249 * t477 + t251 * t478 + t707;
t96 = t394 * t477 - t1040 + t591;
t569 = t83 * t652 + (-t559 * t97 - t562 * t96) * t394;
t568 = t620 * t562 + (t64 * t247 + t998) * t559;
t567 = t589 * t562 + (t49 * t835 - t999) * t559;
t503 = rSges(3,2) * t559 + t929;
t459 = t505 * t559;
t363 = t463 * t863;
t348 = t964 * t562;
t346 = t964 * t559;
t344 = t394 * t562;
t343 = t394 * t559;
t318 = t545 + (-t502 + t503) * qJD(1);
t311 = -rSges(7,2) * t439 + rSges(7,3) * t438;
t308 = -rSges(6,2) * t438 - rSges(6,3) * t439;
t307 = rSges(5,1) * t438 + rSges(5,2) * t439;
t305 = rSges(7,2) * t437 - rSges(7,3) * t436;
t303 = -rSges(5,1) * t436 - rSges(5,2) * t437;
t301 = rSges(6,2) * t436 + rSges(6,3) * t437;
t270 = qJD(3) * t392 + qJD(4) * t462;
t269 = -qJD(1) * t319 + t538;
t268 = qJD(1) * t802 + t811;
t202 = t559 * t644 + t873;
t185 = t202 * qJD(1);
t180 = t477 * t309;
t127 = -rSges(5,3) * t755 + t771;
t126 = -rSges(5,3) * t1036 + t702;
t95 = -t559 * t962 + t592 * t562;
t94 = t592 * t559 + t562 * t962;
t93 = qJD(3) * t648 - t256 * t558 + t259 * t561;
t92 = -qJD(3) * t650 - t257 * t558 + t260 * t561;
t70 = t616 - t781;
t69 = t185 + t617;
t42 = -t126 * t534 - t251 * t739 + t270 * t477 + t357 * t394 + t583;
t41 = t127 * t534 - t270 * t478 - t356 * t394 + (t249 * t787 - t874) * qJD(3) + t613;
t31 = qJD(3) * t586 + t126 * t478 - t127 * t477 - t249 * t357 + t251 * t356 + t278;
t11 = t125 * t534 + t838 * t478 - t813 * t356 + (t244 * t787 - t874) * qJD(3) + t585;
t10 = t123 * t478 + t247 * t356 + t859 * t477 + t836 * t357 + (t586 - t746) * qJD(3) + t722;
t8 = qJD(6) * t199 + t857 * t534 + t769 * t478 - t765 * t356 + (t787 * t837 - t874) * qJD(3) + t585;
t7 = -qJD(4) * t745 + t858 * t478 + t835 * t356 + t774 * t477 + t768 * t357 + (t558 * t640 + t586) * qJD(3) + t722;
t1 = [(t93 + t94 + t69) * t734 + (-t52 * t1032 - t53 * (t409 + t588) - t49 * t180 - (-t52 * t837 - t53 * t835) * t534 - (-t309 * t49 + t999) * t477 + t9 * (t414 + t548 - t536 + t930) + t52 * (t614 + t985) + t8 * (t607 + t982) + t53 * (t582 + t984) + (t9 * (-t561 * t944 + t941) + t52 * t944 * t792) * t562 + (t8 * t723 + t9 * t957) * t559 + ((-t52 * pkin(1) + t53 * t723) * t562 + (t52 * (t729 + t933 + t938) + t53 * t957) * t559) * qJD(1) + (-t199 * t52 + t439 * t9) * (pkin(4) + t1041)) * m(7) + t914 / 0.2e1 + t915 / 0.2e1 + t912 / 0.2e1 + t913 / 0.2e1 + t910 / 0.2e1 + t911 / 0.2e1 + (t129 * (-t552 + t728) - t173 * t801 + t128 * (t760 + t393) + t174 * (-rSges(4,2) * t753 + t763 + t803) + (qJD(3) * t173 * t505 + t129 * t703) * t562 + (t942 * t878 + (t173 * (-qJ(2) - t703) + t174 * (-pkin(7) + t942)) * t559) * qJD(1) - (-t173 + t470 + (t395 - t937) * qJD(1) + t800) * t174) * m(4) + (qJD(1) * t197 + t92 + t95) * t789 / 0.2e1 + (t781 + (t369 * t559 ^ 2 + (-t358 + t145 + (t369 + t650) * t562) * t562) * qJD(3) + t70) * t732 + (t185 + ((-t147 + t358 + t145) * t559 + (t148 - t993 + (t369 - t650) * t559 + t144) * t562) * qJD(3)) * t735 + ((-t1011 * t560 - t1012 * t557) * t561 - t1013 * t558 + t1062) * t478 * t950 + (-t64 * t180 - (-t309 * t64 - t998) * t477 + t12 * (t642 + t833) + t11 * (t607 + t809) + (t11 * t758 + t12 * t957) * t559 + (rSges(6,1) * t752 + t614 - t700 - t940 + (-t1068 + (t729 + t935) * t559) * qJD(1)) * t67 + (-t67 - t247 * t534 + t588 - t582 - t773 + (-t957 * t559 - t562 * t758) * qJD(1)) * t983) * m(6) - ((-t649 + t202) * t559 + t562 * t203) * t778 / 0.2e1 + t978 * t955 + t979 * t956 + t1018 * t951 + t924 / 0.2e1 + t925 / 0.2e1 + t922 / 0.2e1 + t923 / 0.2e1 + (-(qJD(1) * t503 - t318 + t800) * t319 + t269 * (t559 * t943 + t548 + t929) + t318 * t546 + t268 * (t725 + t799) + t319 * (t802 + t803) + (t318 * t943 * t562 + (t318 * (-rSges(3,3) - qJ(2)) - t319 * pkin(1)) * t559) * qJD(1)) * m(3) + (-qJD(3) * t644 + t480 * t558 - t481 * t561) * qJD(1) + t969 + t920 / 0.2e1 + t921 / 0.2e1 + (t1019 + t981) * t953 + ((t1064 + t1114) * t478 + t1043) * t954 + ((-t1017 + t1114) * t477 + t980 + t1113) * t952 + (t42 * (-t251 + t642) + t96 * (rSges(5,3) * t752 + t675 - t702) + t41 * (t717 + t807) + t97 * (t676 + t771) + (t41 * t759 + t42 * t957) * t559 + ((-t96 * pkin(1) + t759 * t97) * t562 + (t96 * (t729 + t927) + t97 * t957) * t559) * qJD(1) - (t608 - t96 - t1040) * t97 - t909 * t477) * m(5); 0.2e1 * (t8 * t946 + t9 * t947) * m(7) + 0.2e1 * (t11 * t946 + t12 * t947) * m(6) + 0.2e1 * (t41 * t946 + t42 * t947) * m(5) + 0.2e1 * (-t883 / 0.2e1 + t882 / 0.2e1) * m(4) + 0.2e1 * (t268 * t946 + t269 * t947) * m(3); ((-t791 * t873 - t796) * t559 + (t1029 + (t559 * t448 + t1031) * qJD(3)) * t562) * t735 + ((t448 * t789 - t796) * t562 + (-t1029 + (-t562 * t873 - t1031) * qJD(3)) * t559) * t732 - qJD(1) * ((-t558 * t805 - t561 * t804) * qJD(1) + (t558 * t595 + t561 * t596) * qJD(3)) / 0.2e1 + qJD(1) * (t93 * t559 + t92 * t562 + (t197 * t562 + t559 * t649) * qJD(1)) / 0.2e1 + (t9 * t810 + t7 * t829 + (qJD(1) * t727 + t7 * t718 + t600) * t559 + (t8 * (-t509 - t765) + t7 * t835) * t562 - (t558 * t567 + t561 * t590) * qJD(4) + (-t434 - t827 * t534 - (t458 - t812) * t478 + (-t561 * t640 + t769) * t562 + t1034) * t53 + (-t524 * t559 - (t365 + t826) * t534 - t812 * t477 + t814 * t793 + t1038) * t52 + (t558 * t783 + t826 * t478 - (-t869 - t827) * t477 + (-t313 + t774 + (t468 + t766) * qJD(1)) * t559 + (qJD(1) * t718 + t858) * t562 + t1033) * t49) * m(7) + (t12 * t810 + t10 * t829 + (qJD(1) * t916 + t10 * t767 + t696) * t559 + (t11 * (-t509 - t813) + t10 * t247) * t562 - (t558 * t568 + t561 * t621) * qJD(4) - (t346 * t534 - t705 - (-t399 + t458) * t478 + (t525 - t487 + t838) * t562 + t1034) * t983 + (-t399 * t477 - (-t348 + t365) * t534 - t964 * t793 + t1038) * t67 + (-t348 * t478 - (t346 - t869) * t477 + (-t313 + t859 + (t468 + t833) * qJD(1)) * t559 + (qJD(1) * t767 + t123) * t562 + t1033) * t64) * m(6) + (-t96 * (t344 * t534 + t392 * t477 + t724) - t97 * (t343 * t534 - t392 * t478 + t705) - ((t249 * t97 - t251 * t96) * t561 + t569 * t558) * qJD(4) + t42 * t868 + t96 * t806 + t97 * t474 - t31 * t404 + (qJD(1) * t909 + t96 * t270 + t31 * t832 + t42 * t394) * t559 + (t41 * (-t394 - t509) + t97 * (-t270 - t487) + t31 * t251 + t96 * t394 * qJD(1)) * t562 + ((-t127 - t313 + (-t251 + t468) * qJD(1)) * t559 + (qJD(1) * t832 + t126) * t562 + t343 * t477 + t344 * t478 + t1055) * t83) * m(5) + (-(t173 * t465 + t174 * t459) * qJD(1) - (t176 * (-t459 * t559 - t465 * t562) - t662 * t703) * qJD(3) + 0.2e1 * t176 * (t271 * t562 - t272 * t559 + (-t393 * t562 + t395 * t559) * qJD(1)) - t662 * t483 + (-t883 + t882 + (t174 * t559 + t878) * qJD(1)) * t505) * m(4) + t1001 * t956 + t1002 * t955 + ((t1015 * t871 + t561 * t978) * qJD(4) + ((-qJD(4) * t1064 + t1004) * t558 - t1069) * t562 + (-t1005 * t439 + t1006 * t438) * t534 + (-t1008 * t439 + t1010 * t438) * t478 + (t1007 * t439 - t1009 * t438) * t477) * t954 + (-qJD(1) * t971 + t1025 * t559 + t1026 * t562) * t953 + ((-t1065 * t870 + t561 * t979) * qJD(4) + ((qJD(4) * t1017 - t1004) * t558 + t1069) * t559 + (t1005 * t437 - t1006 * t436) * t534 + (t1008 * t437 - t1010 * t436) * t478 + (-t1007 * t437 + t1009 * t436) * t477) * t952 + (-qJD(1) * t970 + t1023 * t559 + t1024 * t562) * t951 + (((t1005 * t560 - t1006 * t557 + t1108) * t534 + (t1008 * t560 - t1010 * t557 + t1083) * t478 + (-t1007 * t560 + t1009 * t557 + t1013) * t477 + t977 * qJD(4)) * t561 + (qJD(4) * t972 - t1101) * t558) * t950 + (-qJD(1) * t972 + t1021 * t559 + t1022 * t562) * t949 - t1020 * t787 / 0.2e1 - t981 * t749 / 0.2e1 + t980 * t748 / 0.2e1 + t1003 * t708 + (qJD(1) * t94 + (t559 * (t594 * t559 - t562 * t963) + t562 * (t593 * t559 + t562 * t961) + (-t147 * t559 + t148 * t562) * qJD(1)) * t1042 + t1028) * t947 + (qJD(1) * t95 + (t559 * (t559 * t963 + t594 * t562) + t562 * (-t559 * t961 + t593 * t562) + (-t144 * t559 + t145 * t562) * qJD(1)) * t1042 + t1027) * t562 / 0.2e1 - (t617 + t69 + t981) * t795 / 0.2e1 + (t616 + t70 + t980) * t793 / 0.2e1; (-t53 * (qJD(6) * t438 + t305 * t534 + t839 + (-t460 - t464) * t478) - t49 * (t311 * t478 - t745 + t853 + (-t302 - t305) * t477) - (t49 * t645 + t53 * t632 - t926) * qJ(6) + t9 * t363 + t8 * t830 + t53 * t720 + t49 * t770 + (qJD(3) * t567 + t53 * t857 + t766 * t9 + t8 * t837) * t558 + (t590 * qJD(3) + (t7 * t768 + t49 * t774 + (t49 * t766 + t727) * qJD(1) + t600) * t562 + (qJD(1) * t589 + t49 * t775 + t53 * t831 + t7 * t766 + t8 * t814) * t559) * t561 + (qJD(6) * t436 - t464 * t477 + t775 * t558 + (t306 + t311) * t534 + t988) * t52) * m(7) + (t12 * t363 + t11 * t830 - t983 * t720 + t64 * t770 + (qJD(3) * t568 + t11 * t244 + t12 * t833 - t125 * t983) * t558 + (t621 * qJD(3) + (t10 * t836 + t64 * t859 + (t64 * t833 + t916) * qJD(1) + t696) * t562 + (qJD(1) * t620 + t10 * t833 - t11 * t964 - t274 * t983 + t64 * t860) * t559) * t561 + t983 * (t301 * t534 + (-t460 - t461) * t478 + t839) - t64 * (t308 * t478 + (-t301 - t302) * t477 + t853) + (t860 * t558 - t461 * t477 - (-t306 - t308) * t534 + t988) * t67) * m(6) + (-t96 * (-t307 * t534 + t462 * t477) - t97 * (t303 * t534 - t462 * t478) - t83 * (-t303 * t477 + t307 * t478) + (qJD(3) * t569 - t96 * t126 + t97 * t127 + t41 * t249 - t42 * t251) * t558 + (t96 * (-qJD(3) * t251 + t270 * t562) + t97 * (qJD(3) * t249 + t270 * t559) - t31 * t652 + t83 * (-t126 * t559 - t127 * t562 + t249 * t795 - t251 * t793) + (t41 * t559 + t42 * t562 + (-t559 * t96 + t562 * t97) * qJD(1)) * t394) * t561) * m(5) + (t558 * t979 - t561 * t970) * t956 + (t558 * t978 - t561 * t971) * t955 + (t1056 * t438 + t1057 * t439 + t973 * t863) * t954 + ((-qJD(1) * t1002 + t978 * qJD(3) + t1025 * t562 - t1026 * t559) * t561 + (qJD(3) * t971 + t1019) * t558) * t953 + (-t1056 * t436 - t1057 * t437 - t865 * t973) * t952 + ((-qJD(1) * t1001 + t979 * qJD(3) + t1023 * t562 - t1024 * t559) * t561 + (qJD(3) * t970 + t1018) * t558) * t951 + ((-t1056 * t557 - t1057 * t560) * t561 + t973 * t558) * t950 + ((-qJD(1) * t1003 + t977 * qJD(3) + t1021 * t562 - t1022 * t559) * t561 + (qJD(3) * t972 + t1053) * t558) * t949 + (t914 + t915 + t924 + t925 + t912 + t913 + t922 + t923 + t910 + t911 + t920 + t921 + t969) * t558 / 0.2e1 - t1027 * t865 / 0.2e1 + t1028 * t863 / 0.2e1 + t1020 * t733 + t981 * t558 * t734 + t980 * t558 * t732 + (t558 * t977 - t561 * t972) * t708 - (t559 * t980 + t562 * t981) * t794 / 0.2e1; (t200 * t52 + t436 * t9 - t438 * t8 + t49 * t989 + t53 * t987 + t7 * t872 - t926) * m(7) + (t10 * t872 - t11 * t438 + t12 * t436 - t987 * t983 + (t200 - t232) * t67 + t989 * t64) * m(6); (t437 * t9 - t439 * t8 + t7 * t864 + (-t437 * t534 + t478 * t864 + t199) * t53 + (-t439 * t534 - t477 * t864 + t201) * t52 + (t437 * t477 + t439 * t478 + t610) * t49) * m(7);];
tauc  = t1(:);