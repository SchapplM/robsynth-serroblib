% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP8_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP8_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP8_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:31
% EndTime: 2019-12-31 21:09:14
% DurationCPUTime: 93.22s
% Computational Cost: add. (20544->1157), mult. (54768->1459), div. (0->0), fcn. (52841->6), ass. (0->545)
t1033 = Icges(5,4) - Icges(4,5) - Icges(6,5);
t1032 = Icges(6,4) + Icges(5,5) - Icges(4,6);
t548 = cos(qJ(3));
t1079 = (-Icges(4,4) - Icges(5,6)) * t548;
t545 = sin(qJ(3));
t1078 = (Icges(4,4) - Icges(6,6)) * t545;
t1077 = -Icges(4,1) - Icges(5,2);
t1074 = Icges(4,2) + Icges(6,2);
t1073 = -t1078 + (Icges(4,1) + Icges(6,3)) * t548;
t1072 = -t1079 + (-Icges(4,2) - Icges(5,3)) * t545;
t546 = sin(qJ(2));
t549 = cos(qJ(2));
t1071 = t1033 * t549 + t1073 * t546;
t1070 = t1032 * t549 + t1072 * t546;
t547 = sin(qJ(1));
t550 = cos(qJ(1));
t839 = t550 * t545;
t434 = -t547 * t548 + t549 * t839;
t840 = t549 * t550;
t435 = t545 * t547 + t548 * t840;
t843 = t546 * t550;
t217 = Icges(6,1) * t843 + Icges(6,4) * t434 + Icges(6,5) * t435;
t220 = Icges(5,1) * t843 - Icges(5,4) * t435 + Icges(5,5) * t434;
t223 = Icges(4,5) * t435 - Icges(4,6) * t434 + Icges(4,3) * t843;
t1025 = t217 + t220 + t223;
t867 = Icges(5,6) * t435;
t208 = Icges(5,5) * t843 + Icges(5,3) * t434 - t867;
t415 = Icges(6,6) * t435;
t211 = Icges(6,4) * t843 + Icges(6,2) * t434 + t415;
t879 = Icges(4,4) * t435;
t226 = -Icges(4,2) * t434 + Icges(4,6) * t843 + t879;
t1041 = t208 + t211 - t226;
t864 = Icges(6,6) * t434;
t205 = Icges(6,5) * t843 + Icges(6,3) * t435 + t864;
t418 = Icges(5,6) * t434;
t214 = Icges(5,4) * t843 - Icges(5,2) * t435 + t418;
t421 = Icges(4,4) * t434;
t229 = Icges(4,1) * t435 + Icges(4,5) * t843 - t421;
t1043 = t205 - t214 + t229;
t1009 = t1025 * t843 + t1041 * t434 + t1043 * t435;
t842 = t547 * t549;
t432 = t545 * t842 + t548 * t550;
t841 = t548 * t549;
t433 = t547 * t841 - t839;
t845 = t546 * t547;
t215 = Icges(6,1) * t845 + Icges(6,4) * t432 + Icges(6,5) * t433;
t218 = Icges(5,1) * t845 - Icges(5,4) * t433 + Icges(5,5) * t432;
t221 = Icges(4,5) * t433 - Icges(4,6) * t432 + Icges(4,3) * t845;
t1040 = t215 + t218 + t221;
t417 = Icges(5,6) * t433;
t206 = Icges(5,5) * t845 + Icges(5,3) * t432 - t417;
t414 = Icges(6,6) * t433;
t210 = -Icges(6,4) * t845 - Icges(6,2) * t432 - t414;
t420 = Icges(4,4) * t433;
t224 = -Icges(4,2) * t432 + Icges(4,6) * t845 + t420;
t1042 = t206 - t210 - t224;
t413 = Icges(6,6) * t432;
t203 = Icges(6,5) * t845 + Icges(6,3) * t433 + t413;
t416 = Icges(5,6) * t432;
t213 = -Icges(5,4) * t845 + Icges(5,2) * t433 - t416;
t419 = Icges(4,4) * t432;
t228 = -Icges(4,1) * t433 - Icges(4,5) * t845 + t419;
t1044 = t203 + t213 - t228;
t1010 = t1040 * t843 + t1042 * t434 + t1044 * t435;
t788 = qJD(3) * t550;
t795 = qJD(2) * t547;
t467 = t546 * t788 + t795;
t791 = qJD(3) * t547;
t793 = qJD(2) * t550;
t468 = -t546 * t791 + t793;
t789 = qJD(3) * t549;
t528 = qJD(1) - t789;
t1068 = t1032 * t545 - t1033 * t548;
t1069 = Icges(5,1) + Icges(6,1) + Icges(4,3);
t1024 = t1068 * t546 - t1069 * t549;
t846 = t545 * t546;
t519 = Icges(5,6) * t846;
t844 = t546 * t548;
t1051 = -Icges(5,2) * t844 - t1071 + t519;
t518 = Icges(6,6) * t844;
t1065 = Icges(6,2) * t846 - t1070 + t518;
t948 = t1024 * t843 - t1051 * t435 + t1065 * t434;
t979 = t1009 * t467 - t1010 * t468 + t948 * t528;
t1011 = t1025 * t845 + t1041 * t432 + t1043 * t433;
t1012 = t1040 * t845 + t1042 * t432 + t1044 * t433;
t949 = t1024 * t845 - t1051 * t433 + t1065 * t432;
t980 = t1011 * t467 - t1012 * t468 + t949 * t528;
t1064 = (t1074 * t548 + t1078) * t546;
t1063 = (t1077 * t545 + t1079) * t546;
t1058 = t1068 * t549 + t1069 * t546;
t862 = Icges(6,6) * t548;
t667 = Icges(6,2) * t545 + t862;
t973 = (t667 - t1072) * t549 + t1032 * t546;
t866 = Icges(5,6) * t545;
t668 = Icges(5,2) * t548 - t866;
t972 = (t668 + t1073) * t549 - t1033 * t546;
t1055 = (-t1032 * t548 - t1033 * t545) * t546;
t748 = t548 * t788;
t750 = t545 * t791;
t757 = t546 * t793;
t197 = qJD(1) * t432 + t545 * t757 - t549 * t748 - t750;
t790 = qJD(3) * t548;
t797 = qJD(1) * t550;
t608 = t545 * t797 + t547 * t790;
t798 = qJD(1) * t547;
t611 = t549 * t798 + t757;
t749 = t545 * t788;
t198 = t548 * t611 + t549 * t749 - t608;
t755 = t549 * t793;
t763 = t546 * t798;
t609 = t755 - t763;
t104 = Icges(5,4) * t609 + Icges(5,2) * t198 - Icges(5,6) * t197;
t114 = -Icges(4,1) * t198 + Icges(4,4) * t197 + Icges(4,5) * t609;
t98 = Icges(6,5) * t609 - Icges(6,6) * t197 - Icges(6,3) * t198;
t1050 = -t104 + t114 + t98;
t758 = t546 * t795;
t199 = -t545 * t758 - t548 * t798 + t549 * t608 - t749;
t200 = qJD(1) * t435 - t548 * t758 - t549 * t750 - t748;
t794 = qJD(2) * t549;
t610 = t546 * t797 + t547 * t794;
t105 = Icges(5,4) * t610 - Icges(5,2) * t200 + Icges(5,6) * t199;
t115 = Icges(4,1) * t200 - Icges(4,4) * t199 + Icges(4,5) * t610;
t99 = Icges(6,5) * t610 + Icges(6,6) * t199 + Icges(6,3) * t200;
t1049 = -t105 + t115 + t99;
t100 = Icges(5,5) * t609 + Icges(5,6) * t198 - Icges(5,3) * t197;
t102 = Icges(6,4) * t609 - Icges(6,2) * t197 - Icges(6,6) * t198;
t112 = -Icges(4,4) * t198 + Icges(4,2) * t197 + Icges(4,6) * t609;
t1048 = t100 + t102 - t112;
t101 = Icges(5,5) * t610 - Icges(5,6) * t200 + Icges(5,3) * t199;
t103 = Icges(6,4) * t610 + Icges(6,2) * t199 + Icges(6,6) * t200;
t113 = Icges(4,4) * t200 - Icges(4,2) * t199 + Icges(4,6) * t610;
t1047 = t101 + t103 - t113;
t106 = Icges(6,1) * t609 - Icges(6,4) * t197 - Icges(6,5) * t198;
t108 = Icges(5,1) * t609 + Icges(5,4) * t198 - Icges(5,5) * t197;
t110 = -Icges(4,5) * t198 + Icges(4,6) * t197 + Icges(4,3) * t609;
t1046 = t106 + t108 + t110;
t107 = Icges(6,1) * t610 + Icges(6,4) * t199 + Icges(6,5) * t200;
t109 = Icges(5,1) * t610 - Icges(5,4) * t200 + Icges(5,5) * t199;
t111 = Icges(4,5) * t200 - Icges(4,6) * t199 + Icges(4,3) * t610;
t1045 = t107 + t109 + t111;
t1039 = t1058 * qJD(2) - t1055 * qJD(3);
t792 = qJD(3) * t546;
t1038 = (Icges(5,3) * t548 + t866) * t792 + t1064 * qJD(3) + t973 * qJD(2);
t1037 = (-Icges(6,3) * t545 + t862) * t792 + t1063 * qJD(3) + t972 * qJD(2);
t1034 = -t1051 * t548 + t1065 * t545;
t1002 = rSges(6,1) + pkin(4);
t543 = t550 * pkin(6);
t495 = pkin(1) * t547 - t543;
t478 = qJD(1) * t495;
t505 = pkin(7) * t755;
t535 = pkin(6) * t797;
t494 = pkin(2) * t546 - pkin(7) * t549;
t759 = t494 * t793;
t918 = pkin(2) * t549;
t496 = pkin(7) * t546 + t918;
t461 = t496 * t547;
t952 = qJD(1) * t461;
t1029 = -pkin(2) * t757 + t478 + t505 + t535 + t759 + t952;
t986 = t1040 * t609 - t1042 * t197 - t1044 * t198 + t1045 * t843 + t1047 * t434 + t1049 * t435;
t985 = t1025 * t609 - t1041 * t197 - t1043 * t198 + t1046 * t843 + t1048 * t434 + t1050 * t435;
t984 = t1040 * t610 + t1042 * t199 + t1044 * t200 + t1045 * t845 + t1047 * t432 + t1049 * t433;
t983 = t1025 * t610 + t1041 * t199 + t1043 * t200 + t1046 * t845 + t1048 * t432 + t1050 * t433;
t1015 = t1024 * t610 + t1037 * t433 + t1038 * t432 + t1039 * t845 - t1051 * t200 + t1065 * t199;
t1014 = t1024 * t609 + t1037 * t435 + t1038 * t434 + t1039 * t843 + t1051 * t198 - t1065 * t197;
t659 = -t224 * t545 - t228 * t548;
t993 = t221 * t549;
t80 = t546 * t659 - t993;
t663 = t203 * t548 - t210 * t545;
t999 = t215 * t549;
t82 = t546 * t663 - t999;
t661 = t206 * t545 + t213 * t548;
t996 = t218 * t549;
t84 = t546 * t661 - t996;
t1028 = t80 + t82 + t84;
t658 = -t226 * t545 + t229 * t548;
t81 = -t223 * t549 + t546 * t658;
t662 = t205 * t548 + t211 * t545;
t83 = -t217 * t549 + t546 * t662;
t660 = t208 * t545 - t214 * t548;
t85 = -t220 * t549 + t546 * t660;
t1027 = t81 + t83 + t85;
t947 = -t1024 * t549 + t1034 * t546;
t1023 = -t546 * t667 + t1070;
t1022 = t546 * t668 + t1071;
t1021 = (-t1034 + t1058) * t528 + (t1024 * t547 + t659 + t661 + t663) * t468 + (-t1024 * t550 - t658 - t660 - t662) * t467;
t1001 = rSges(6,3) + qJ(5);
t1018 = t547 * t550;
t698 = rSges(6,2) * t545 + rSges(6,3) * t548;
t990 = -qJ(5) * t844 + t1002 * t549 - t546 * t698;
t410 = qJD(4) * t434;
t116 = -t198 * pkin(3) - qJ(4) * t197 + t410;
t412 = t432 * qJ(4);
t302 = pkin(3) * t433 + t412;
t1017 = t302 * t528 + t1029 + t116 - t410;
t490 = rSges(3,1) * t546 + rSges(3,2) * t549;
t801 = t550 * pkin(1) + t547 * pkin(6);
t989 = qJD(1) * t801;
t1016 = t490 * t795 - t989;
t1013 = (t1034 * qJD(2) - t1039) * t549 + (t1037 * t548 + t1038 * t545 + (t1051 * t545 + t1065 * t548) * qJD(3) + t1024 * qJD(2)) * t546;
t1008 = Icges(6,3) - t1077;
t1007 = Icges(5,3) + t1074;
t1004 = t1055 * t528 + (t1032 * t433 + t1033 * t432) * t468 + (-t1032 * t435 - t1033 * t434) * t467;
t308 = t435 * pkin(3) + qJ(4) * t434;
t531 = pkin(2) * t840;
t463 = pkin(7) * t843 + t531;
t766 = qJD(1) * t463 - t494 * t795 + t989;
t697 = pkin(3) * t548 + qJ(4) * t545;
t962 = t697 * t546;
t677 = t528 * t308 - t467 * t962 + t766;
t406 = qJD(5) * t433;
t787 = qJD(4) * t432;
t807 = -t787 - t406;
t829 = t434 * rSges(6,2) + t1001 * t435 + t1002 * t843;
t51 = t467 * t990 + t829 * t528 + t677 - t807;
t1003 = 0.2e1 * qJD(2);
t246 = t435 * rSges(4,1) - t434 * rSges(4,2) + rSges(4,3) * t843;
t702 = rSges(4,1) * t548 - rSges(4,2) * t545;
t386 = -rSges(4,3) * t549 + t546 * t702;
t92 = t246 * t528 - t386 * t467 + t766;
t1000 = qJD(1) * t92;
t359 = t547 * t962;
t781 = qJD(2) * qJD(3);
t741 = t549 * t781;
t349 = qJD(1) * t467 + t547 * t741;
t350 = qJD(1) * t468 + t550 * t741;
t742 = t546 * t781;
t988 = t1009 * t350 + t1010 * t349 + t1014 * t528 + t467 * t985 - t468 * t986 + t742 * t948;
t987 = t1011 * t350 + t1012 * t349 + t1015 * t528 + t467 * t983 - t468 * t984 + t742 * t949;
t25 = (qJD(2) * t659 - t111) * t549 + (qJD(2) * t221 - t113 * t545 + t115 * t548 + (-t224 * t548 + t228 * t545) * qJD(3)) * t546;
t27 = (qJD(2) * t663 - t107) * t549 + (qJD(2) * t215 + t103 * t545 + t548 * t99 + (-t203 * t545 - t210 * t548) * qJD(3)) * t546;
t29 = (qJD(2) * t661 - t109) * t549 + (qJD(2) * t218 + t101 * t545 - t105 * t548 + (t206 * t548 - t213 * t545) * qJD(3)) * t546;
t982 = t25 + t27 + t29;
t26 = (qJD(2) * t658 - t110) * t549 + (qJD(2) * t223 - t112 * t545 + t114 * t548 + (-t226 * t548 - t229 * t545) * qJD(3)) * t546;
t28 = (qJD(2) * t662 - t106) * t549 + (qJD(2) * t217 + t102 * t545 + t548 * t98 + (-t205 * t545 + t211 * t548) * qJD(3)) * t546;
t30 = (qJD(2) * t660 - t108) * t549 + (qJD(2) * t220 + t100 * t545 - t104 * t548 + (t208 * t548 + t214 * t545) * qJD(3)) * t546;
t981 = t26 + t28 + t30;
t978 = t1027 * t467 - t1028 * t468 + t528 * t947;
t977 = t1023 * t547;
t976 = t1023 * t550;
t975 = t1022 * t547;
t974 = t1022 * t550;
t971 = t1021 * t546;
t970 = t1024 * t528 + t1025 * t467 - t1040 * t468;
t943 = t1027 * t550 + t1028 * t547;
t968 = t1027 * t547 - t1028 * t550;
t942 = t1009 * t550 + t1010 * t547;
t941 = t1011 * t550 + t1012 * t547;
t967 = t1009 * t547 - t1010 * t550;
t966 = t1011 * t547 - t1012 * t550;
t703 = -rSges(4,1) * t433 + rSges(4,2) * t432;
t244 = rSges(4,3) * t845 - t703;
t965 = -t244 * t528 - t386 * t468;
t964 = pkin(6) * qJD(1);
t639 = t432 * t528 + t468 * t846;
t958 = -t197 + t639;
t638 = -t434 * t528 + t467 * t846;
t957 = t199 + t638;
t911 = rSges(6,2) * t432;
t832 = t1001 * t433 + t1002 * t845 + t911;
t539 = Icges(3,4) * t549;
t673 = -Icges(3,2) * t546 + t539;
t488 = Icges(3,1) * t546 + t539;
t407 = qJD(5) * t435;
t956 = -t197 * rSges(6,2) - t1001 * t198 + t1002 * t755 + t407;
t764 = -pkin(1) - t918;
t780 = -pkin(7) - t1002;
t612 = t780 * t546 + t764;
t243 = rSges(5,1) * t843 - t435 * rSges(5,2) + t434 * rSges(5,3);
t699 = rSges(5,2) * t548 - rSges(5,3) * t545;
t586 = rSges(5,1) * t549 + t546 * t699;
t64 = t243 * t528 + t467 * t586 + t677 + t787;
t946 = (Icges(5,3) * t844 + t1051 + t1064 + t519) * t528 + (-t1007 * t433 + t1044 + t413 - t416 - t419) * t468 + (t1007 * t435 - t1043 + t418 + t421 - t864) * t467;
t945 = (-Icges(6,3) * t846 + t1063 + t1065 + t518) * t528 + (t1008 * t432 - t1042 - t414 + t417 + t420) * t468 + (-t1008 * t434 + t1041 + t415 - t867 - t879) * t467;
t944 = t1004 * t546;
t940 = t1013 * t528 + t742 * t947;
t606 = t545 * t794 + t546 * t790;
t861 = Icges(3,3) * t550;
t363 = Icges(3,5) * t842 - Icges(3,6) * t845 - t861;
t522 = Icges(3,4) * t845;
t874 = Icges(3,5) * t550;
t371 = Icges(3,1) * t842 - t522 - t874;
t869 = Icges(3,6) * t550;
t367 = Icges(3,4) * t842 - Icges(3,2) * t845 - t869;
t852 = t367 * t546;
t655 = -t371 * t549 + t852;
t141 = -t363 * t550 - t547 * t655;
t485 = Icges(3,5) * t549 - Icges(3,6) * t546;
t484 = Icges(3,5) * t546 + Icges(3,6) * t549;
t615 = qJD(2) * t484;
t880 = Icges(3,4) * t546;
t489 = Icges(3,1) * t549 - t880;
t372 = Icges(3,5) * t547 + t489 * t550;
t368 = Icges(3,6) * t547 + t550 * t673;
t851 = t368 * t546;
t654 = -t372 * t549 + t851;
t939 = -t550 * t615 + (-t485 * t547 + t654 + t861) * qJD(1);
t364 = Icges(3,3) * t547 + t485 * t550;
t938 = -t547 * t615 + (t364 + t655) * qJD(1);
t486 = Icges(3,2) * t549 + t880;
t651 = t486 * t546 - t488 * t549;
t937 = qJD(1) * t651 + t485 * qJD(2);
t187 = t199 * qJ(4);
t117 = pkin(3) * t200 + t187 + t787;
t785 = qJD(4) * t546;
t513 = t548 * t785;
t786 = qJD(4) * t545;
t514 = t549 * t786;
t936 = qJD(2) * t514 + qJD(3) * t513 + t467 * t117 + t350 * t302;
t935 = t547 * (-t486 * t550 + t372) - t550 * (-Icges(3,2) * t842 + t371 - t522);
t931 = t349 / 0.2e1;
t930 = t350 / 0.2e1;
t929 = -t467 / 0.2e1;
t928 = t467 / 0.2e1;
t927 = -t468 / 0.2e1;
t926 = t468 / 0.2e1;
t925 = -t528 / 0.2e1;
t924 = t528 / 0.2e1;
t920 = -rSges(5,1) - pkin(7);
t919 = -rSges(4,3) - pkin(7);
t916 = rSges(3,1) * t549;
t915 = rSges(5,1) * t546;
t912 = rSges(6,2) * t199;
t910 = rSges(4,3) * t546;
t148 = t432 * t467 + t434 * t468;
t784 = qJD(5) * t548;
t511 = t546 * t784;
t512 = t545 * t785;
t810 = t461 * t795 + t463 * t793;
t723 = t467 * t302 + t512 + t810;
t771 = t308 + t829;
t49 = t467 * t832 + t468 * t771 + t511 + t723;
t908 = t148 * t49;
t907 = t25 * t468;
t906 = t26 * t467;
t905 = t27 * t468;
t904 = t28 * t467;
t903 = t29 * t468;
t902 = t30 * t467;
t540 = t547 * rSges(3,3);
t700 = rSges(5,2) * t433 - rSges(5,3) * t432;
t240 = rSges(5,1) * t845 - t700;
t828 = t243 + t308;
t62 = t240 * t467 + t468 * t828 + t723;
t895 = t62 * t240;
t595 = (-t461 - t495) * qJD(1) - t759;
t573 = -t468 * t962 + t410 + t595;
t830 = -t240 - t302;
t63 = t468 * t586 + t528 * t830 + t573;
t894 = t63 * t586;
t893 = t80 * t349;
t892 = t81 * t350;
t891 = t82 * t349;
t890 = t83 * t350;
t889 = t84 * t349;
t888 = t85 * t350;
t884 = t117 * t843 + t302 * t755;
t802 = rSges(3,2) * t845 + t550 * rSges(3,3);
t387 = rSges(3,1) * t842 - t802;
t760 = t490 * t793;
t185 = -t760 + (-t387 - t495) * qJD(1);
t856 = t185 * t547;
t855 = t185 * t550;
t389 = rSges(3,1) * t840 - rSges(3,2) * t843 + t540;
t186 = qJD(1) * t389 - t1016;
t459 = t490 * t550;
t854 = t186 * t459;
t476 = qJD(2) * t496;
t849 = t476 * t550;
t848 = t484 * t547;
t847 = t484 * t550;
t276 = t549 * t302;
t775 = rSges(5,1) * t755 + t198 * rSges(5,2) - t197 * rSges(5,3);
t119 = -rSges(5,1) * t763 + t775;
t838 = t116 + t119;
t837 = -t1002 * t763 + t956;
t836 = t1001 * t200 + t1002 * t610 + t406 + t912;
t300 = -pkin(3) * t432 + qJ(4) * t433;
t835 = t467 * t300 + t513;
t306 = -pkin(3) * t434 + qJ(4) * t435;
t834 = qJD(4) * t433 + t528 * t306;
t607 = -t545 * t792 + t548 * t794;
t237 = pkin(3) * t607 + qJ(4) * t606 + t512;
t391 = -t549 * t699 + t915;
t453 = (rSges(5,2) * t545 + rSges(5,3) * t548) * t546;
t268 = qJD(2) * t391 + qJD(3) * t453;
t833 = -t237 - t268;
t388 = t549 * t702 + t910;
t454 = (-rSges(4,1) * t545 - rSges(4,2) * t548) * t546;
t264 = qJD(2) * t388 + qJD(3) * t454;
t827 = -t264 - t476;
t390 = rSges(6,1) * t546 + t549 * t698;
t456 = (rSges(6,2) * t548 - rSges(6,3) * t545) * t546;
t796 = qJD(2) * t546;
t826 = pkin(4) * t796 + qJ(5) * t607 + qJD(2) * t390 + qJD(3) * t456 + t511;
t385 = t962 * t798;
t825 = t308 * t796 + t546 * t385;
t824 = t845 * t962 + t276;
t310 = -pkin(2) * t611 - pkin(7) * t763 + t505;
t464 = qJD(1) * (-pkin(1) * t798 + t535);
t823 = qJD(1) * t310 + t464;
t822 = t990 * t547;
t821 = t990 * t550;
t820 = -t547 * t363 - t371 * t840;
t819 = t547 * t364 + t372 * t840;
t466 = t494 * t798;
t816 = t385 + t466;
t815 = -t386 - t494;
t814 = -pkin(4) * t546 - qJ(5) * t841 - t390;
t812 = t586 - t962;
t460 = t494 * t547;
t462 = t494 * t550;
t811 = -t460 * t795 - t462 * t793;
t808 = t547 * t461 + t550 * t463;
t452 = (-pkin(3) * t545 + qJ(4) * t548) * t546;
t806 = qJD(4) * t435 - t468 * t452;
t805 = -t486 + t489;
t804 = t488 + t673;
t803 = rSges(3,2) * t763 + rSges(3,3) * t797;
t799 = qJD(1) * t485;
t201 = -t547 * t651 - t847;
t783 = t201 * qJD(1);
t782 = qJD(1) * qJD(2);
t779 = -pkin(3) - t1001;
t778 = t116 + t837;
t777 = t468 * t302;
t776 = -t198 * rSges(4,1) + t197 * rSges(4,2) + rSges(4,3) * t755;
t774 = -t237 - t826;
t773 = -t476 + t833;
t772 = -t302 - t832;
t504 = pkin(2) * t758;
t311 = pkin(7) * t610 + qJD(1) * t531 - t504;
t770 = t311 * t795 + (t310 + t952) * t793;
t769 = t550 * t310 + t547 * t311 + t461 * t797;
t768 = -t962 + t990;
t767 = -t494 + t812;
t752 = qJD(5) * t846;
t751 = t549 * t784;
t744 = -pkin(1) - t916;
t743 = t547 * t782;
t739 = t797 / 0.2e1;
t738 = t796 / 0.2e1;
t737 = -t795 / 0.2e1;
t736 = t795 / 0.2e1;
t734 = t793 / 0.2e1;
t731 = t49 * t832;
t50 = t468 * t990 + t528 * t772 + t407 + t573;
t730 = t50 * t990;
t729 = t64 * t812;
t343 = t372 * t842;
t727 = t364 * t550 - t343;
t726 = -t363 + t851;
t724 = t549 * t117 + t237 * t845 + t610 * t962;
t722 = -t476 + t774;
t721 = t547 * t302 + t550 * t308 + t808;
t720 = -t494 + t768;
t719 = t801 + t463;
t712 = t51 * t768;
t709 = qJD(3) * t738;
t708 = -qJD(1) * t462 - t496 * t795;
t707 = qJD(4) * t199 + t528 * t116 + t308 * t742 + t823;
t705 = -rSges(3,2) * t546 + t916;
t704 = rSges(4,1) * t200 - rSges(4,2) * t199;
t701 = -rSges(5,2) * t200 + rSges(5,3) * t199;
t121 = rSges(5,1) * t610 + t701;
t596 = -t463 * t743 + t770;
t10 = t121 * t467 + t240 * t350 - t349 * t828 + t468 * t838 + t596 + t936;
t696 = t10 * t240 + t62 * t121;
t676 = t788 * t276 - t467 * t359 + t514 + t811;
t664 = -t186 * t547 - t855;
t657 = t244 * t550 - t246 * t547;
t195 = t367 * t549 + t371 * t546;
t196 = t368 * t549 + t372 * t546;
t650 = t494 * t743 + (-t311 - t989) * qJD(1);
t649 = -pkin(1) - t496;
t647 = t550 * t116 + t547 * t117 + t302 * t797 + t769;
t628 = -t476 * t547 - t494 * t797;
t457 = t490 * t547;
t625 = -t62 * t828 - t894;
t624 = t64 * t243 + t63 * t830;
t623 = qJD(1) * t460 - t496 * t793;
t620 = t546 * t920 + t764;
t619 = t546 * t919 + t764;
t360 = t550 * t962;
t618 = t308 * t792 - t528 * t360 + t708;
t617 = qJD(2) * t488;
t616 = qJD(2) * t486;
t142 = -t368 * t845 - t727;
t614 = (-t141 * t550 + t142 * t547) * qJD(2);
t143 = -t367 * t843 - t820;
t144 = -t368 * t843 + t819;
t613 = (-t143 * t550 + t144 * t547) * qJD(2);
t168 = (t387 * t547 + t389 * t550) * qJD(2);
t605 = t308 + t719;
t7 = -qJD(3) * t752 + t836 * t467 + t832 * t350 + (-t463 * t798 + t751) * qJD(2) + t778 * t468 - t771 * t349 + t770 + t936;
t597 = t49 * t836 + t7 * t832;
t594 = t367 * t550 - t368 * t547;
t593 = -qJD(4) * t197 - t468 * t237 + t349 * t962 + t650;
t458 = t697 * t549;
t592 = t359 * t789 - t468 * t458 + t623;
t591 = -t49 * t771 - t730;
t590 = t50 * t772 + t51 * t829;
t574 = (-t546 * t804 + t549 * t805) * qJD(1);
t79 = t244 * t467 + t246 * t468 + t810;
t91 = t595 + t965;
t561 = t79 * t657 + (t91 * t547 - t550 * t92) * t386;
t470 = t673 * qJD(2);
t471 = t489 * qJD(2);
t558 = qJD(1) * t484 - t470 * t546 + t471 * t549 + (-t486 * t549 - t488 * t546) * qJD(2);
t557 = (t729 + t895) * t550 + t625 * t547;
t556 = -t546 * t935 + t594 * t549;
t555 = (t712 + t731) * t550 + t591 * t547;
t472 = t705 * qJD(2);
t342 = t586 * t550;
t340 = t586 * t547;
t338 = t386 * t550;
t337 = t386 * t547;
t309 = rSges(6,2) * t435 - rSges(6,3) * t434;
t307 = -rSges(4,1) * t434 - rSges(4,2) * t435;
t305 = rSges(5,2) * t434 + rSges(5,3) * t435;
t303 = rSges(6,2) * t433 - rSges(6,3) * t432;
t301 = -rSges(4,1) * t432 - rSges(4,2) * t433;
t299 = rSges(5,2) * t432 + rSges(5,3) * t433;
t266 = -qJD(2) * t457 + (t550 * t705 + t540) * qJD(1);
t265 = -rSges(3,1) * t611 - rSges(3,2) * t755 + t803;
t263 = t302 * t843;
t202 = -t550 * t651 + t848;
t182 = t202 * qJD(1);
t134 = -t472 * t793 + (-t266 + t1016) * qJD(1);
t133 = -t472 * t795 + t464 + (t265 - t760) * qJD(1);
t123 = rSges(4,3) * t610 + t704;
t122 = -rSges(4,3) * t763 + t776;
t89 = t558 * t547 - t550 * t937;
t88 = t547 * t937 + t558 * t550;
t87 = -qJD(2) * t654 + (-t550 * t616 + (-t547 * t673 + t869) * qJD(1)) * t549 + (-t550 * t617 + (-t489 * t547 + t874) * qJD(1)) * t546;
t86 = -qJD(2) * t655 + (qJD(1) * t368 - t547 * t616) * t549 + (qJD(1) * t372 - t547 * t617) * t546;
t66 = t182 + t613;
t65 = t614 + t783;
t42 = -t123 * t528 - t264 * t468 + t349 * t386 + (-t244 * t792 - t849) * qJD(2) + t650;
t41 = t122 * t528 - t264 * t467 - t350 * t386 + (t246 * t792 + t628) * qJD(2) + t823;
t31 = t122 * t468 + t123 * t467 + t244 * t350 - t246 * t349 + t596;
t24 = -t268 * t468 - t349 * t586 + (-t117 - t121) * t528 + (t792 * t830 - t849) * qJD(2) + t593;
t23 = t119 * t528 + t833 * t467 + t812 * t350 + (t243 * t792 + t628) * qJD(2) + t707;
t9 = -qJD(5) * t198 - t826 * t468 - t990 * t349 + (-t117 - t836) * t528 + (t772 * t792 - t849) * qJD(2) + t593;
t8 = qJD(5) * t200 + t837 * t528 + t774 * t467 + t768 * t350 + (t792 * t829 + t628) * qJD(2) + t707;
t1 = [(-t62 * t777 - (-t302 * t62 + t729) * t468 + t23 * (t605 + t243) + (-t117 + t504 - t701 + t620 * t797 + (t794 * t920 - t964) * t547) * t63 + (t620 * t547 - t302 + t543 + t700) * t24 + (t240 * t528 + t63 + t775 + (t649 - t915) * t798 + t1017) * t64) * m(5) + (t8 * (t605 + t829) + (t433 * t779 + t547 * t612 - t412 + t543 - t911) * t9 - t49 * t777 - (-t49 * t302 + t712) * t468 + (t528 * t832 + t612 * t798 + t1017 - t407 + t956) * t51 + (-t187 + t504 + t807 - t912 + t779 * t200 + t612 * t797 + (t780 * t794 - t964) * t547 + t51) * t50) * m(6) + (t1015 + t979) * t927 + t948 * t930 + t949 * t931 + t1014 * t928 + (-(-qJD(1) * t387 - t185 - t478 - t760) * t186 + t134 * (t547 * t744 + t543 + t802) + t133 * (t389 + t801) + t186 * (t535 + t803) + (t490 * t856 - t854) * qJD(2) + ((-pkin(1) - t705) * t855 + (t185 * (-rSges(3,3) - pkin(6)) + t186 * t744) * t547) * qJD(1)) * m(3) + (t993 + t996 + t999 + (-t1042 * t545 - t1044 * t548) * t546 + t1028) * t467 * t925 + (-qJD(2) * t651 + t470 * t549 + t471 * t546) * qJD(1) + ((t196 + t202) * t550 + (t195 + t201) * t547) * t782 / 0.2e1 + t979 * t926 - (t86 + t89 + t66) * t793 / 0.2e1 + (t182 + ((t142 - t343 + (t364 + t852) * t550 + t820) * t550 + t819 * t547) * qJD(2)) * t734 + (t87 + t88) * t736 + (t42 * (t543 + t703) + t41 * (t719 + t246) + (t42 * t619 + (t649 - t910) * t1000) * t547 + (t504 - t704 + t619 * t797 + (t794 * t919 - t964) * t547) * t91 + (t776 + t91 - t965 + t1029) * t92) * m(4) + t906 / 0.2e1 - t907 / 0.2e1 + t904 / 0.2e1 - t905 / 0.2e1 + t902 / 0.2e1 - t903 / 0.2e1 + t940 + t893 / 0.2e1 + t891 / 0.2e1 + t892 / 0.2e1 + t889 / 0.2e1 + t890 / 0.2e1 + (-t783 + ((t550 * t726 + t144 - t819) * t550 + (t547 * t726 + t143 + t727) * t547) * qJD(2) + t65) * t737 + t888 / 0.2e1; qJD(1) * (t547 * t87 - t550 * t86 + (t195 * t547 + t196 * t550) * qJD(1)) / 0.2e1 + ((-t793 * t848 - t799) * t550 + (t574 + (t550 * t847 + t556) * qJD(2)) * t547) * t734 + ((-t795 * t847 + t799) * t547 + (t574 + (t547 * t848 + t556) * qJD(2)) * t550) * t737 - qJD(1) * ((t546 * t805 + t549 * t804) * qJD(1) + (t594 * t546 + t549 * t935) * qJD(2)) / 0.2e1 + (t50 * t816 + t7 * t721 + t49 * t647 + (t9 * t720 + t50 * t722 + t7 * t829 + t49 * t837 + (t51 * t720 + t731) * qJD(1)) * t550 + (t8 * t720 + t51 * t722 + (-t730 + t49 * (-t463 - t771)) * qJD(1) + t597) * t547 - t50 * t592 - t51 * t618 - t49 * (t676 + t751) - (t50 * (t359 - t822) + t51 * t821) * t528 - (t50 * t814 + t49 * (-t360 + t821)) * t468 - (t51 * (-t458 + t814) + t49 * t822) * t467 - (t546 * t590 + t549 * t555) * qJD(3) - (t50 * t550 + t51 * t547) * t546 * (-t784 - t786)) * m(6) + (-t64 * (t342 * t528 - t512 * t547 + t618 + (-t391 - t458) * t467) - t62 * (t340 * t467 + t676 + (t342 - t360) * t468) - (t546 * t624 + t549 * t557) * qJD(3) + t10 * t721 + t62 * t647 + (t24 * t767 + t10 * t243 + t62 * t119 + (t64 * t767 + t895) * qJD(1)) * t550 + (t23 * t767 + t64 * t773 + (-t894 + t62 * (-t463 - t828)) * qJD(1) + t696) * t547 + (t391 * t468 - t592 - (-t340 + t359) * t528 + t816 + (t512 + t773) * t550) * t63) * m(5) + (t91 * t466 + t31 * t808 + (t31 * t246 + t91 * t827 + (t1000 + t42) * t815) * t550 + (t91 * t386 * qJD(1) + t31 * t244 + t41 * t815 + t92 * t827) * t547 - t91 * (t337 * t528 - t388 * t468 + t623) - t92 * (-t338 * t528 - t388 * t467 + t708) - ((-t244 * t91 + t246 * t92) * t546 + t561 * t549) * qJD(3) + (t769 + (t244 * qJD(1) + t122) * t550 + (t123 + (-t246 - t463) * qJD(1)) * t547 + t337 * t467 + t338 * t468 - t811) * t79) * m(4) + (0.2e1 * t168 * (t265 * t550 + t266 * t547 + (t387 * t550 - t389 * t547) * qJD(1)) + t664 * t472 + (-t133 * t547 - t134 * t550 + (-t186 * t550 + t856) * qJD(1)) * t490 - (t185 * t457 - t854) * qJD(1) - (t168 * (-t457 * t547 - t459 * t550) + t664 * t705) * qJD(2)) * m(3) + t966 * t931 + t967 * t930 + ((t1010 * t842 + t948 * t546) * qJD(3) + ((qJD(3) * t1009 + t970) * t549 + t971) * t550 + (t434 * t973 + t435 * t972) * t528 + (-t434 * t977 + t435 * t975) * t468 + (t434 * t976 - t435 * t974) * t467) * t929 + (qJD(1) * t942 + t547 * t985 - t550 * t986) * t928 + (qJD(1) * t941 + t547 * t983 - t550 * t984) * t927 + ((t1011 * t840 + t949 * t546) * qJD(3) + ((qJD(3) * t1012 + t970) * t549 + t971) * t547 + (t432 * t973 + t433 * t972) * t528 + (-t432 * t977 + t433 * t975) * t468 + (t432 * t976 - t433 * t974) * t467) * t926 + ((qJD(3) * t943 - t1021) * t549 + ((t545 * t973 + t548 * t972 + t1024) * t528 + (-t545 * t977 + t548 * t975 - t1040) * t468 + (t545 * t976 - t548 * t974 + t1025) * t467 + t947 * qJD(3)) * t546) * t925 + (qJD(1) * t943 + t547 * t981 - t550 * t982) * t924 - t978 * t792 / 0.2e1 + t968 * t709 + (qJD(1) * t88 + (t547 ^ 2 * t939 - t938 * t1018 + (t143 * t547 + t144 * t550) * qJD(1)) * t1003 + t988) * t547 / 0.2e1 - (qJD(1) * t89 + (-t939 * t1018 + t550 ^ 2 * t938 + (t141 * t547 + t142 * t550) * qJD(1)) * t1003 + t987) * t550 / 0.2e1 + (t614 + t65 + t980) * t798 / 0.2e1 + (t613 + t66 + t979) * t739 - (t547 * t980 + t550 * t979) * t789 / 0.2e1; (t9 * t824 + t50 * t724 + t51 * t825 + t7 * t263 + t49 * t884 + (qJD(2) * t555 + t50 * t836 - t51 * t778 - t771 * t8 + t832 * t9) * t549 + (t590 * qJD(2) + (qJD(1) * t591 + t51 * t774 + t768 * t8 + t597) * t550 + (-t9 * t990 + t50 * t826 - t7 * t771 - t49 * t778 + (t49 * t772 - t51 * t990) * qJD(1)) * t547) * t546 - t50 * (-qJD(5) * t434 - t456 * t468 + t806 + (-t300 - t303) * t528) - t51 * (-qJD(5) * t432 + t309 * t528 + t834 + (-t452 - t456) * t467) - t49 * (t303 * t467 - t752 + t835 + (t306 + t309) * t468) - (t50 * t639 + t51 * t638 - t908) * qJ(5)) * m(6) + (t24 * t824 + t63 * t724 + t64 * t825 + t10 * t263 + t62 * t884 + (qJD(2) * t557 + t63 * t121 - t23 * t828 + t24 * t240 - t64 * t838) * t549 + (t624 * qJD(2) + (qJD(1) * t625 + t23 * t812 + t64 * t833 + t696) * t550 + (-t24 * t586 + t63 * t268 - t10 * t828 - t62 * t838 + (-t586 * t64 + t62 * t830) * qJD(1)) * t547) * t546 - t63 * (-t453 * t468 + (-t299 - t300) * t528 + t806) - t64 * (t305 * t528 + (-t452 - t453) * t467 + t834) - t62 * (t299 * t467 + (t305 + t306) * t468 + t835)) * m(5) + ((qJD(2) * t561 - t92 * t122 + t91 * t123 + t42 * t244 - t41 * t246) * t549 + (t91 * (-qJD(2) * t244 + t264 * t547) + t92 * (qJD(2) * t246 - t264 * t550) + t31 * t657 + t79 * (-t122 * t547 + t123 * t550 - t244 * t798 - t246 * t797) + (-t41 * t550 + t42 * t547 + (t92 * t547 + t91 * t550) * qJD(1)) * t386) * t546 - t91 * (-t301 * t528 - t454 * t468) - t92 * (t307 * t528 - t454 * t467) - t79 * (t301 * t467 + t307 * t468)) * m(4) + (t546 * t941 - t549 * t949) * t931 + (t546 * t942 - t549 * t948) * t930 + (t434 * t946 + t435 * t945 - t550 * t944) * t929 + ((qJD(2) * t942 - t1014) * t549 + (-qJD(1) * t967 + t948 * qJD(2) + t547 * t986 + t550 * t985) * t546) * t928 + ((qJD(2) * t941 - t1015) * t549 + (-qJD(1) * t966 + t949 * qJD(2) + t547 * t984 + t550 * t983) * t546) * t927 + (t432 * t946 + t433 * t945 - t547 * t944) * t926 + (t1004 * t549 + (t545 * t946 + t548 * t945) * t546) * t925 + ((qJD(2) * t943 - t1013) * t549 + (-qJD(1) * t968 + t947 * qJD(2) + t547 * t982 + t550 * t981) * t546) * t924 - (t892 + t893 + t906 - t907 + t890 + t891 + t904 - t905 + t888 + t889 + t902 - t903 + t940) * t549 / 0.2e1 + t987 * t845 / 0.2e1 + t988 * t843 / 0.2e1 + t978 * t738 + (t546 * t943 - t549 * t947) * t709 + (t546 * t739 + t549 * t736) * t980 + (-t763 / 0.2e1 + t549 * t734) * t979; (t432 * t8 + t434 * t9 + t49 * t606 + t50 * t958 + t51 * t957 + t7 * t846 - t908) * m(6) + (t10 * t846 + t23 * t432 + t24 * t434 + t957 * t64 + t958 * t63 + (-t148 + t606) * t62) * m(5); (t433 * t8 + t435 * t9 + t7 * t844 + (-t435 * t528 + t467 * t844 + t200) * t51 + (t433 * t528 + t468 * t844 - t198) * t50 + (-t433 * t467 - t435 * t468 + t607) * t49) * m(6);];
tauc = t1(:);
