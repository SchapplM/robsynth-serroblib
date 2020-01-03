% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRRPP7
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
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRRPP7_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP7_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP7_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP7_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:03:57
% EndTime: 2019-12-31 21:05:45
% DurationCPUTime: 97.87s
% Computational Cost: add. (20475->1195), mult. (54561->1489), div. (0->0), fcn. (52562->6), ass. (0->570)
t1085 = Icges(6,4) + Icges(5,5);
t1026 = -Icges(5,1) - Icges(6,1);
t1018 = Icges(4,1) - t1026;
t1082 = Icges(6,2) + Icges(5,3);
t1046 = Icges(5,4) + Icges(4,5) - Icges(6,5);
t1045 = Icges(4,6) - Icges(5,6) + Icges(6,6);
t541 = cos(qJ(3));
t1084 = t1085 * t541;
t538 = sin(qJ(3));
t1083 = (Icges(4,4) - t1085) * t538;
t1081 = -t1082 * t538 - t1084;
t1017 = Icges(4,2) + t1082;
t1080 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t1079 = t1045 * t538 - t1046 * t541;
t1078 = t1018 * t541 - t1083;
t543 = cos(qJ(1));
t840 = t541 * t543;
t540 = sin(qJ(1));
t542 = cos(qJ(2));
t842 = t540 * t542;
t429 = t538 * t842 + t840;
t838 = t543 * t538;
t841 = t541 * t542;
t430 = t540 * t841 - t838;
t539 = sin(qJ(2));
t845 = t539 * t540;
t207 = Icges(4,5) * t430 - Icges(4,6) * t429 + Icges(4,3) * t845;
t213 = Icges(5,4) * t430 + Icges(5,2) * t845 + Icges(5,6) * t429;
t1077 = t207 + t213;
t844 = t539 * t541;
t1076 = t1085 * t844;
t880 = Icges(4,4) * t541;
t659 = -Icges(4,2) * t538 + t880;
t1075 = -t1045 * t542 + t539 * t659;
t410 = Icges(5,5) * t430;
t205 = -Icges(5,6) * t845 - Icges(5,3) * t429 - t410;
t416 = Icges(4,4) * t430;
t216 = -Icges(4,2) * t429 + Icges(4,6) * t845 + t416;
t1074 = t205 + t216;
t409 = Icges(5,5) * t429;
t222 = Icges(5,1) * t430 + Icges(5,4) * t845 + t409;
t415 = Icges(4,4) * t429;
t226 = -Icges(4,1) * t430 - Icges(4,5) * t845 + t415;
t1073 = t222 - t226;
t1048 = t1079 * t539 + t1080 * t542;
t1072 = -t1079 * t542 + t1080 * t539;
t846 = t538 * t539;
t1044 = t1082 * t846 - t1075 + t1076;
t980 = (-t659 - t1081) * t542 - t1045 * t539;
t1033 = -t1046 * t542 + t1078 * t539;
t979 = t1046 * t539 + t1078 * t542;
t1071 = (t1045 * t541 + t1046 * t538) * t539;
t1070 = (t1017 * t541 + t1083) * t539;
t839 = t542 * t543;
t432 = t538 * t540 + t541 * t839;
t411 = Icges(5,5) * t432;
t431 = -t540 * t541 + t542 * t838;
t843 = t539 * t543;
t206 = Icges(5,6) * t843 + Icges(5,3) * t431 + t411;
t882 = Icges(4,4) * t432;
t218 = -Icges(4,2) * t431 + Icges(4,6) * t843 + t882;
t1069 = t206 - t218;
t209 = Icges(4,5) * t432 - Icges(4,6) * t431 + Icges(4,3) * t843;
t215 = Icges(5,4) * t432 + Icges(5,2) * t843 + Icges(5,6) * t431;
t1068 = t209 + t215;
t873 = Icges(5,5) * t431;
t224 = Icges(5,1) * t432 + Icges(5,4) * t843 + t873;
t417 = Icges(4,4) * t431;
t227 = Icges(4,1) * t432 + Icges(4,5) * t843 - t417;
t1067 = t224 + t227;
t1066 = t1077 * t845;
t1065 = t1077 * t843;
t1064 = t1073 * t430 - t1074 * t429 + t1066;
t1063 = -t1073 * t432 + t1074 * t431 - t1065;
t1030 = t1067 * t430 + t1068 * t845 + t1069 * t429;
t203 = Icges(6,5) * t432 + Icges(6,6) * t431 - Icges(6,3) * t843;
t414 = Icges(6,4) * t432;
t212 = Icges(6,2) * t431 - Icges(6,6) * t843 + t414;
t878 = Icges(6,4) * t431;
t221 = Icges(6,1) * t432 - Icges(6,5) * t843 + t878;
t834 = -t429 * t212 - t430 * t221;
t68 = -t203 * t845 - t834;
t1042 = t68 + t1030;
t1029 = t1067 * t432 + t1068 * t843 + t1069 * t431;
t832 = t431 * t212 + t432 * t221;
t74 = -t203 * t843 + t832;
t1041 = t74 + t1029;
t778 = qJD(3) * t543;
t742 = t542 * t778;
t783 = qJD(2) * t543;
t747 = t539 * t783;
t781 = qJD(3) * t540;
t195 = qJD(1) * t429 - t541 * t742 + (t747 - t781) * t538;
t779 = qJD(3) * t542;
t521 = qJD(1) - t779;
t636 = t521 * t538;
t788 = qJD(1) * t542;
t713 = -qJD(3) + t788;
t196 = t543 * t636 + (-t540 * t713 - t747) * t541;
t745 = t542 * t783;
t789 = qJD(1) * t540;
t753 = t539 * t789;
t600 = -t745 + t753;
t102 = Icges(4,5) * t196 + Icges(4,6) * t195 - Icges(4,3) * t600;
t106 = Icges(5,4) * t196 - Icges(5,2) * t600 - Icges(5,6) * t195;
t98 = Icges(6,5) * t196 - Icges(6,6) * t195 + Icges(6,3) * t600;
t1062 = t102 + t106 - t98;
t785 = qJD(2) * t540;
t748 = t539 * t785;
t780 = qJD(3) * t541;
t787 = qJD(1) * t543;
t197 = t780 * t842 - t541 * t789 + (t542 * t787 - t748 - t778) * t538;
t786 = qJD(2) * t539;
t198 = t713 * t840 + (-t541 * t786 + t636) * t540;
t784 = qJD(2) * t542;
t601 = t539 * t787 + t540 * t784;
t103 = Icges(4,5) * t198 - Icges(4,6) * t197 + Icges(4,3) * t601;
t107 = Icges(5,4) * t198 + Icges(5,2) * t601 + Icges(5,6) * t197;
t99 = Icges(6,5) * t198 + Icges(6,6) * t197 - Icges(6,3) * t601;
t1061 = t103 + t107 - t99;
t100 = Icges(5,5) * t196 - Icges(5,6) * t600 - Icges(5,3) * t195;
t104 = Icges(6,4) * t196 - Icges(6,2) * t195 + Icges(6,6) * t600;
t108 = Icges(4,4) * t196 + Icges(4,2) * t195 - Icges(4,6) * t600;
t1060 = t100 + t104 - t108;
t101 = Icges(5,5) * t198 + Icges(5,6) * t601 + Icges(5,3) * t197;
t105 = Icges(6,4) * t198 + Icges(6,2) * t197 - Icges(6,6) * t601;
t109 = Icges(4,4) * t198 - Icges(4,2) * t197 + Icges(4,6) * t601;
t1059 = t101 + t105 - t109;
t110 = Icges(6,1) * t196 - Icges(6,4) * t195 + Icges(6,5) * t600;
t112 = Icges(5,1) * t196 - Icges(5,4) * t600 - Icges(5,5) * t195;
t114 = Icges(4,1) * t196 + Icges(4,4) * t195 - Icges(4,5) * t600;
t1058 = t110 + t112 + t114;
t111 = Icges(6,1) * t198 + Icges(6,4) * t197 - Icges(6,5) * t601;
t113 = Icges(5,1) * t198 + Icges(5,4) * t601 + Icges(5,5) * t197;
t115 = Icges(4,1) * t198 - Icges(4,4) * t197 + Icges(4,5) * t601;
t1057 = t111 + t113 + t115;
t950 = t1033 * t432 + t1044 * t431 - t1048 * t843;
t202 = -Icges(6,5) * t430 - Icges(6,6) * t429 + Icges(6,3) * t845;
t1056 = t202 + t1077;
t1055 = t203 - t1068;
t413 = Icges(6,4) * t430;
t211 = -Icges(6,2) * t429 + Icges(6,6) * t845 - t413;
t1016 = t211 + t1074;
t1054 = t212 + t1069;
t412 = Icges(6,4) * t429;
t219 = Icges(6,1) * t430 - Icges(6,5) * t845 + t412;
t1053 = t219 + t1073;
t1052 = t221 + t1067;
t1051 = t1072 * qJD(2) - t1071 * qJD(3);
t1050 = t980 * qJD(2) + t1070 * qJD(3);
t446 = (-Icges(4,1) * t538 - t880) * t539;
t782 = qJD(3) * t539;
t1049 = qJD(3) * t446 + (t1026 * t538 + t1084) * t782 + t979 * qJD(2);
t951 = -t1033 * t430 - t1044 * t429 + t1048 * t845;
t1047 = t1033 * t541 + t1044 * t538;
t1028 = rSges(6,3) + qJ(5);
t536 = t543 * pkin(6);
t490 = pkin(1) * t540 - t536;
t475 = qJD(1) * t490;
t501 = pkin(7) * t745;
t527 = pkin(6) * t787;
t489 = pkin(2) * t539 - pkin(7) * t542;
t749 = t489 * t783;
t917 = pkin(2) * t542;
t491 = pkin(7) * t539 + t917;
t458 = t491 * t540;
t954 = qJD(1) * t458;
t1043 = t501 + t527 + t475 + t749 + t954;
t993 = t1016 * t195 + t1053 * t196 - t1056 * t600 + t1057 * t432 + t1059 * t431 + t1061 * t843;
t992 = t1052 * t196 - t1054 * t195 + t1055 * t600 + t1058 * t432 + t1060 * t431 + t1062 * t843;
t991 = -t1016 * t197 + t1053 * t198 + t1056 * t601 + t1057 * t430 + t1059 * t429 + t1061 * t845;
t990 = t1052 * t198 + t1054 * t197 - t1055 * t601 + t1058 * t430 + t1060 * t429 + t1062 * t845;
t1023 = t1033 * t196 - t1044 * t195 + t1048 * t600 + t1049 * t432 + t1050 * t431 + t1051 * t843;
t1022 = t1033 * t198 + t1044 * t197 - t1048 * t601 + t1049 * t430 + t1050 * t429 + t1051 * t845;
t650 = -t211 * t429 + t219 * t430;
t67 = t202 * t845 + t650;
t1020 = t67 + t1064;
t833 = t211 * t431 - t432 * t219;
t73 = t202 * t843 - t833;
t1019 = t73 - t1063;
t1003 = t542 * t202;
t649 = -t211 * t538 + t219 * t541;
t80 = t539 * t649 - t1003;
t1006 = t213 * t542;
t652 = -t205 * t538 + t222 * t541;
t82 = t539 * t652 - t1006;
t1009 = t207 * t542;
t647 = -t216 * t538 - t226 * t541;
t84 = t539 * t647 - t1009;
t1040 = t80 + t82 + t84;
t648 = t212 * t538 + t221 * t541;
t81 = t203 * t542 + t539 * t648;
t651 = t206 * t538 + t224 * t541;
t83 = -t215 * t542 + t539 * t651;
t646 = -t218 * t538 + t227 * t541;
t85 = -t209 * t542 + t539 * t646;
t1039 = t81 + t83 + t85;
t465 = -t539 * t781 + t783;
t1038 = t1063 * t465 + t950 * t521;
t949 = t1047 * t539 + t1048 * t542;
t1034 = t1081 * t539 + t1075;
t464 = t539 * t778 + t785;
t1032 = (-t1047 + t1072) * t521 + (-t1048 * t540 + t647 + t649 + t652) * t465 + (t1048 * t543 - t646 - t648 - t651) * t464;
t1031 = t1064 * t465 + t951 * t521;
t1011 = rSges(6,1) + pkin(4);
t823 = t431 * rSges(6,2) + t1011 * t432 - t1028 * t843;
t1027 = t540 * t543;
t686 = rSges(6,1) * t541 + rSges(6,2) * t538;
t1001 = -pkin(4) * t844 - t1028 * t542 - t539 * t686;
t406 = qJD(4) * t431;
t116 = t196 * pkin(3) - qJ(4) * t195 + t406;
t408 = t429 * qJ(4);
t301 = pkin(3) * t430 + t408;
t1025 = t301 * t521 + t1043 + t116 - t406;
t485 = rSges(3,1) * t539 + rSges(3,2) * t542;
t793 = t543 * pkin(1) + t540 * pkin(6);
t996 = qJD(1) * t793;
t1024 = t485 * t785 - t996;
t987 = t1042 * t464 - t465 * t67 - t1031;
t986 = t1041 * t464 - t465 * t73 + t1038;
t1021 = (t1047 * qJD(2) - t1051) * t542 + (t1049 * t541 + t1050 * t538 + (-t1033 * t538 + t1044 * t541) * qJD(3) - t1048 * qJD(2)) * t539;
t1014 = t1071 * t521 + (-t1045 * t430 - t1046 * t429) * t465 + (t1045 * t432 + t1046 * t431) * t464;
t307 = t432 * pkin(3) + qJ(4) * t431;
t523 = pkin(2) * t839;
t460 = pkin(7) * t843 + t523;
t756 = qJD(1) * t460 - t489 * t785 + t996;
t777 = qJD(4) * t429;
t685 = pkin(3) * t541 + qJ(4) * t538;
t968 = t685 * t539;
t1000 = -t521 * t307 + t464 * t968 - t756 - t777;
t775 = qJD(5) * t539;
t504 = t540 * t775;
t54 = t1001 * t464 + t823 * t521 - t1000 - t504;
t1013 = -t539 * (t202 * t540 + t203 * t543) + t832;
t1012 = 0.2e1 * qJD(2);
t244 = t432 * rSges(4,1) - t431 * rSges(4,2) + rSges(4,3) * t843;
t690 = rSges(4,1) * t541 - rSges(4,2) * t538;
t386 = -rSges(4,3) * t542 + t539 * t690;
t92 = t244 * t521 - t386 * t464 + t756;
t1010 = qJD(1) * t92;
t357 = t540 * t968;
t999 = rSges(6,2) * t429 - t1028 * t845;
t383 = t968 * t789;
t455 = t685 * t542;
t463 = t489 * t789;
t998 = -t357 * t779 + t465 * t455 + t383 + t463;
t772 = qJD(2) * qJD(3);
t734 = t542 * t772;
t347 = qJD(1) * t464 + t540 * t734;
t348 = qJD(1) * t465 + t543 * t734;
t735 = t539 * t772;
t995 = t1019 * t347 + t1023 * t521 + t1041 * t348 + t464 * t992 - t465 * t993 + t735 * t950;
t994 = t1020 * t347 + t1022 * t521 + t1042 * t348 + t464 * t990 - t465 * t991 - t735 * t951;
t25 = (qJD(2) * t649 + t99) * t542 + (qJD(2) * t202 + t105 * t538 + t111 * t541 + (-t211 * t541 - t219 * t538) * qJD(3)) * t539;
t27 = (qJD(2) * t652 - t107) * t542 + (qJD(2) * t213 + t101 * t538 + t113 * t541 + (-t205 * t541 - t222 * t538) * qJD(3)) * t539;
t29 = (qJD(2) * t647 - t103) * t542 + (qJD(2) * t207 - t109 * t538 + t115 * t541 + (-t216 * t541 + t226 * t538) * qJD(3)) * t539;
t989 = t25 + t27 + t29;
t26 = (qJD(2) * t648 + t98) * t542 + (-qJD(2) * t203 + t104 * t538 + t110 * t541 + (t212 * t541 - t221 * t538) * qJD(3)) * t539;
t28 = (qJD(2) * t651 - t106) * t542 + (qJD(2) * t215 + t100 * t538 + t112 * t541 + (t206 * t541 - t224 * t538) * qJD(3)) * t539;
t30 = (qJD(2) * t646 - t102) * t542 + (qJD(2) * t209 - t108 * t538 + t114 * t541 + (-t218 * t541 - t227 * t538) * qJD(3)) * t539;
t988 = t26 + t28 + t30;
t985 = t1039 * t464 - t1040 * t465 + t521 * t949;
t984 = t1034 * t540;
t983 = t1034 * t543;
t982 = t1033 * t540;
t981 = t1033 * t543;
t978 = t1032 * t539;
t977 = -t1048 * t521 - t1055 * t464 - t1056 * t465;
t944 = t1039 * t543 + t1040 * t540;
t976 = t1039 * t540 - t1040 * t543;
t943 = t1019 * t540 + t1041 * t543;
t975 = -t1019 * t543 + t1041 * t540;
t942 = t1020 * t540 + t1042 * t543;
t974 = -t1020 * t543 + t1042 * t540;
t691 = -rSges(4,1) * t430 + rSges(4,2) * t429;
t238 = rSges(4,3) * t845 - t691;
t972 = -t238 * t521 - t386 * t465;
t971 = t834 + (-t202 * t543 + t203 * t540) * t539;
t970 = pkin(6) * qJD(1);
t626 = t429 * t521 + t465 * t846;
t965 = -t195 + t626;
t625 = -t431 * t521 + t464 * t846;
t964 = t197 + t625;
t826 = t1011 * t430 + t999;
t358 = t543 * t968;
t963 = t307 * t782 - t521 * t358;
t961 = -t465 * t968 + t406;
t473 = qJD(2) * t491;
t528 = qJD(5) * t542;
t716 = -t473 - t528;
t959 = -rSges(6,2) * t197 + t504;
t532 = Icges(3,4) * t542;
t660 = -Icges(3,2) * t539 + t532;
t483 = Icges(3,1) * t539 + t532;
t957 = -t195 * rSges(6,2) + t1011 * t196 + t1028 * t753;
t243 = t432 * rSges(5,1) + rSges(5,2) * t843 + t431 * rSges(5,3);
t687 = rSges(5,1) * t541 + rSges(5,3) * t538;
t385 = -rSges(5,2) * t542 + t539 * t687;
t64 = t243 * t521 - t385 * t464 - t1000;
t776 = qJD(4) * t539;
t505 = t538 * t776;
t598 = t538 * t784 + t539 * t780;
t599 = -t538 * t782 + t541 * t784;
t235 = pkin(3) * t599 + qJ(4) * t598 + t505;
t948 = -qJD(4) * t195 - t465 * t235 + t347 * t968;
t947 = (-t1033 + t1070) * t521 + (-t1017 * t430 + t1053 + t409 + t412 - t415) * t465 + (t1017 * t432 - t1052 + t417 - t873 - t878) * t464;
t946 = (t1026 * t846 + t1044 + t1076 + t446) * t521 + (t1018 * t429 + t1016 - t410 - t413 + t416) * t465 + (-t1018 * t431 + t1054 + t411 + t414 - t882) * t464;
t945 = t1014 * t539;
t941 = t1021 * t521 + t735 * t949;
t863 = Icges(3,3) * t543;
t365 = Icges(3,5) * t842 - Icges(3,6) * t845 - t863;
t515 = Icges(3,4) * t845;
t875 = Icges(3,5) * t543;
t381 = Icges(3,1) * t842 - t515 - t875;
t867 = Icges(3,6) * t543;
t373 = Icges(3,4) * t842 - Icges(3,2) * t845 - t867;
t852 = t373 * t539;
t641 = -t381 * t542 + t852;
t141 = -t365 * t543 - t540 * t641;
t480 = Icges(3,5) * t542 - Icges(3,6) * t539;
t479 = Icges(3,5) * t539 + Icges(3,6) * t542;
t605 = qJD(2) * t479;
t883 = Icges(3,4) * t539;
t484 = Icges(3,1) * t542 - t883;
t382 = Icges(3,5) * t540 + t484 * t543;
t374 = Icges(3,6) * t540 + t543 * t660;
t851 = t374 * t539;
t640 = -t382 * t542 + t851;
t940 = -t543 * t605 + (-t480 * t540 + t640 + t863) * qJD(1);
t366 = Icges(3,3) * t540 + t480 * t543;
t939 = -t540 * t605 + (t366 + t641) * qJD(1);
t481 = Icges(3,2) * t542 + t883;
t639 = t481 * t539 - t483 * t542;
t938 = qJD(1) * t639 + t480 * qJD(2);
t830 = t197 * qJ(4) + t777;
t117 = pkin(3) * t198 + t830;
t506 = t541 * t776;
t507 = qJD(4) * t542 * t538;
t937 = qJD(2) * t507 + qJD(3) * t506 + t464 * t117 + t348 * t301;
t936 = t540 * (-t481 * t543 + t382) - t543 * (-Icges(3,2) * t842 + t381 - t515);
t932 = -m(6) / 0.2e1;
t931 = m(6) / 0.2e1;
t930 = t347 / 0.2e1;
t929 = t348 / 0.2e1;
t928 = -t464 / 0.2e1;
t927 = t464 / 0.2e1;
t926 = -t465 / 0.2e1;
t925 = t465 / 0.2e1;
t924 = -t521 / 0.2e1;
t923 = t521 / 0.2e1;
t919 = -rSges(5,2) - pkin(7);
t918 = -rSges(4,3) - pkin(7);
t915 = rSges(3,1) * t542;
t914 = rSges(5,2) * t539;
t911 = rSges(4,3) * t539;
t148 = t429 * t464 + t431 * t465;
t803 = t458 * t785 + t460 * t783;
t712 = t464 * t301 + t505 + t803;
t761 = t307 + t823;
t49 = t464 * t826 + t465 * t761 + t528 + t712;
t908 = t148 * t49;
t907 = t25 * t465;
t906 = t26 * t464;
t905 = t27 * t465;
t904 = t28 * t464;
t903 = t29 * t465;
t902 = t30 * t464;
t533 = t540 * rSges(3,3);
t688 = rSges(5,1) * t430 + rSges(5,3) * t429;
t237 = rSges(5,2) * t845 + t688;
t822 = t243 + t307;
t62 = t237 * t464 + t465 * t822 + t712;
t897 = t62 * t237;
t717 = (-t458 - t490) * qJD(1);
t587 = t717 - t749;
t825 = -t237 - t301;
t63 = -t385 * t465 + t521 * t825 + t587 + t961;
t896 = t63 * t385;
t895 = t80 * t347;
t894 = t81 * t348;
t893 = t82 * t347;
t892 = t83 * t348;
t891 = t84 * t347;
t890 = t85 * t348;
t885 = t117 * t843 + t301 * t745;
t794 = rSges(3,2) * t845 + t543 * rSges(3,3);
t387 = rSges(3,1) * t842 - t794;
t750 = t485 * t783;
t182 = -t750 + (-t387 - t490) * qJD(1);
t858 = t182 * t540;
t857 = t182 * t543;
t391 = rSges(3,1) * t839 - rSges(3,2) * t843 + t533;
t183 = qJD(1) * t391 - t1024;
t456 = t485 * t543;
t856 = t183 * t456;
t849 = t473 * t543;
t848 = t479 * t540;
t847 = t479 * t543;
t765 = t196 * rSges(5,1) + rSges(5,2) * t745 - t195 * rSges(5,3);
t119 = -rSges(5,2) * t753 + t765;
t837 = t116 + t119;
t836 = -rSges(6,3) * t745 + (-qJ(5) * t784 - t775) * t543 + t957;
t835 = t1011 * t198 - t1028 * t601 - t959;
t297 = -pkin(3) * t429 + qJ(4) * t430;
t831 = t464 * t297 + t506;
t303 = -pkin(3) * t431 + qJ(4) * t432;
t828 = qJD(4) * t430 + t521 * t303;
t389 = t542 * t687 + t914;
t451 = (-rSges(5,1) * t538 + rSges(5,3) * t541) * t539;
t263 = qJD(2) * t389 + qJD(3) * t451;
t827 = -t235 - t263;
t388 = -rSges(6,3) * t539 + t542 * t686;
t450 = (-rSges(6,1) * t538 + rSges(6,2) * t541) * t539;
t821 = pkin(4) * t599 - qJ(5) * t786 + qJD(2) * t388 + qJD(3) * t450 + t528;
t390 = t542 * t690 + t911;
t452 = (-rSges(4,1) * t538 - rSges(4,2) * t541) * t539;
t264 = qJD(2) * t390 + qJD(3) * t452;
t820 = -t264 - t473;
t819 = t307 * t786 + t539 * t383;
t818 = t542 * t301 + t845 * t968;
t602 = -t540 * t788 - t747;
t308 = pkin(2) * t602 - pkin(7) * t753 + t501;
t461 = qJD(1) * (-pkin(1) * t789 + t527);
t817 = qJD(1) * t308 + t461;
t500 = pkin(2) * t748;
t309 = pkin(7) * t601 + qJD(1) * t523 - t500;
t816 = -t309 - t996;
t815 = t1001 * t540;
t814 = t1001 * t543;
t813 = -t540 * t365 - t381 * t839;
t812 = t540 * t366 + t382 * t839;
t807 = -t385 - t968;
t806 = -t386 - t489;
t805 = -pkin(4) * t841 + qJ(5) * t539 - t388;
t457 = t489 * t540;
t459 = t489 * t543;
t804 = -t457 * t785 - t459 * t783;
t801 = t540 * t458 + t543 * t460;
t449 = (-pkin(3) * t538 + qJ(4) * t541) * t539;
t800 = qJD(4) * t432 - t465 * t449;
t798 = -t481 + t484;
t797 = t483 + t660;
t796 = rSges(3,2) * t753 + rSges(3,3) * t787;
t791 = qJD(1) * t457;
t790 = qJD(1) * t480;
t199 = -t540 * t639 - t847;
t774 = t199 * qJD(1);
t773 = qJD(1) * qJD(2);
t771 = -pkin(3) - t1011;
t770 = -pkin(7) + t1028;
t769 = t116 + t836;
t768 = t465 * t301;
t767 = t196 * rSges(4,1) + t195 * rSges(4,2) + rSges(4,3) * t745;
t764 = -t235 - t821;
t763 = -t473 + t827;
t762 = -t301 - t826;
t760 = t309 * t785 + (t308 + t954) * t783;
t759 = t543 * t308 + t540 * t309 + t458 * t787;
t758 = -t968 + t1001;
t757 = -t489 + t807;
t755 = -pkin(1) - t917;
t754 = t489 * t787;
t741 = t543 * t775;
t737 = -pkin(1) - t915;
t736 = t540 * t773;
t732 = t787 / 0.2e1;
t731 = t786 / 0.2e1;
t730 = -t785 / 0.2e1;
t729 = t785 / 0.2e1;
t727 = t783 / 0.2e1;
t724 = t49 * t826;
t53 = (-qJD(2) * t489 - t775) * t543 + t1001 * t465 + t717 + t762 * t521 + t961;
t723 = t53 * t1001;
t722 = t64 * t807;
t341 = t382 * t842;
t719 = t366 * t543 - t341;
t718 = -t365 + t851;
t715 = pkin(2) * t747;
t714 = t542 * t117 + t235 * t845 + t601 * t968;
t711 = -t473 + t764;
t710 = t540 * t301 + t543 * t307 + t801;
t709 = -t489 + t758;
t708 = t793 + t460;
t701 = t54 * t758;
t698 = qJD(3) * t731;
t427 = qJD(1) * t459;
t697 = -t491 * t785 - t427;
t689 = rSges(5,1) * t198 + rSges(5,3) * t197;
t122 = rSges(5,2) * t601 + t689;
t588 = -t460 * t736 + t760;
t8 = t122 * t464 + t237 * t348 - t347 * t822 + t465 * t837 + t588 + t937;
t696 = t62 * t122 + t8 * t237;
t695 = qJD(4) * t197 + t521 * t116 + t307 * t735 + t817;
t693 = -rSges(3,2) * t539 + t915;
t692 = rSges(4,1) * t198 - rSges(4,2) * t197;
t665 = t301 * t742 - t464 * t357 + t507 + t804;
t653 = -t183 * t540 - t857;
t645 = t238 * t543 - t244 * t540;
t193 = t373 * t542 + t381 * t539;
t194 = t374 * t542 + t382 * t539;
t433 = t489 * t736;
t638 = qJD(1) * t816 + t433;
t637 = -pkin(1) - t491;
t635 = t543 * t116 + t540 * t117 + t301 * t787 + t759;
t615 = -t473 * t540 - t754;
t454 = t485 * t540;
t613 = -t62 * t822 + t896;
t612 = t64 * t243 + t63 * t825;
t611 = -t491 * t783 + t791;
t609 = t539 * t919 + t755;
t608 = t539 * t918 + t755;
t607 = qJD(2) * t483;
t606 = qJD(2) * t481;
t142 = -t374 * t845 - t719;
t604 = (-t141 * t543 + t142 * t540) * qJD(2);
t143 = -t373 * t843 - t813;
t144 = -t374 * t843 + t812;
t603 = (-t143 * t543 + t144 * t540) * qJD(2);
t168 = (t387 * t540 + t391 * t543) * qJD(2);
t597 = t307 + t708;
t7 = t835 * t464 + t826 * t348 + (-t460 * t789 - t775) * qJD(2) + t769 * t465 - t761 * t347 + t760 + t937;
t589 = t49 * t835 + t7 * t826;
t586 = t373 * t543 - t374 * t540;
t583 = -t49 * t761 - t723;
t582 = t53 * t762 + t54 * t823;
t581 = -t505 + t716;
t565 = (-t539 * t797 + t542 * t798) * qJD(1);
t79 = t238 * t464 + t244 * t465 + t803;
t91 = t587 + t972;
t554 = t79 * t645 + (t540 * t91 - t543 * t92) * t386;
t467 = t660 * qJD(2);
t468 = t484 * qJD(2);
t551 = qJD(1) * t479 - t467 * t539 + t468 * t542 + (-t481 * t542 - t483 * t539) * qJD(2);
t550 = (t722 + t897) * t543 + t613 * t540;
t549 = -t539 * t936 + t586 * t542;
t548 = (t701 + t724) * t543 + t583 * t540;
t469 = t693 * qJD(2);
t340 = t386 * t543;
t339 = t385 * t543;
t337 = t386 * t540;
t336 = t385 * t540;
t306 = -rSges(4,1) * t431 - rSges(4,2) * t432;
t305 = -rSges(5,1) * t431 + rSges(5,3) * t432;
t304 = -rSges(6,1) * t431 + rSges(6,2) * t432;
t300 = -rSges(4,1) * t429 - rSges(4,2) * t430;
t299 = -rSges(5,1) * t429 + rSges(5,3) * t430;
t298 = -rSges(6,1) * t429 + rSges(6,2) * t430;
t266 = -qJD(2) * t454 + (t543 * t693 + t533) * qJD(1);
t265 = rSges(3,1) * t602 - rSges(3,2) * t745 + t796;
t261 = t301 * t843;
t200 = -t543 * t639 + t848;
t179 = t200 * qJD(1);
t134 = -t469 * t783 + (-t266 + t1024) * qJD(1);
t133 = -t469 * t785 + t461 + (t265 - t750) * qJD(1);
t123 = rSges(4,3) * t601 + t692;
t120 = -rSges(4,3) * t753 + t767;
t89 = t551 * t540 - t543 * t938;
t88 = t540 * t938 + t551 * t543;
t87 = -qJD(2) * t640 + (-t543 * t606 + (-t540 * t660 + t867) * qJD(1)) * t542 + (-t543 * t607 + (-t484 * t540 + t875) * qJD(1)) * t539;
t86 = -qJD(2) * t641 + (qJD(1) * t374 - t540 * t606) * t542 + (qJD(1) * t382 - t540 * t607) * t539;
t66 = t179 + t603;
t65 = t604 + t774;
t42 = -t123 * t521 - t264 * t465 + t347 * t386 + (-t238 * t782 - t849) * qJD(2) + t638;
t41 = t120 * t521 - t264 * t464 - t348 * t386 + (t244 * t782 + t615) * qJD(2) + t817;
t31 = t120 * t465 + t123 * t464 + t238 * t348 - t244 * t347 + t588;
t24 = -t263 * t465 + t347 * t385 + (-t117 - t122) * t521 + (t782 * t825 - t849) * qJD(2) + t638 + t948;
t23 = t119 * t521 + t827 * t464 + t807 * t348 + (t243 * t782 + t615) * qJD(2) + t695;
t10 = -qJD(1) * t741 + t836 * t521 + t764 * t464 + t758 * t348 + (t540 * t716 + t782 * t823 - t754) * qJD(2) + t695;
t9 = t433 - t821 * t465 - t1001 * t347 + (-t117 - t835) * t521 + (t504 + t816) * qJD(1) + (t543 * t716 + t762 * t782) * qJD(2) + t948;
t1 = [(t134 * (t540 * t737 + t536 + t794) + t133 * (t391 + t793) + t183 * (t527 + t796) + (t485 * t858 - t856) * qJD(2) + ((-pkin(1) - t693) * t857 + (t182 * (-rSges(3,3) - pkin(6)) + t183 * t737) * t540) * qJD(1) - (-qJD(1) * t387 - t182 - t475 - t750) * t183) * m(3) + (t1003 + t1006 + t1009 + (t1016 * t538 - t1053 * t541) * t539 + t1040) * t464 * t924 + (t65 - t774 + ((t543 * t718 + t144 - t812) * t543 + (t540 * t718 + t143 + t719) * t540) * qJD(2)) * t730 + (-t62 * t768 - (-t301 * t62 + t722) * t465 + t23 * (t597 + t243) + (-t117 + t500 - t689 + t609 * t787 + (t784 * t919 - t970) * t540) * t63 + (t609 * t540 - t301 + t536 - t688) * t24 + (t237 * t521 + t63 - t715 + t765 + (t637 - t914) * t789 + t1025) * t64) * m(5) + ((t68 + t833 + t971) * t465 + (t1016 * t429 - t1053 * t430 + t1013 + t1020 + t1029 - t1066) * t464 + t1038) * t925 + ((t194 + t200) * t543 + (t193 + t199) * t540) * t773 / 0.2e1 - (t86 + t89 + t66) * t783 / 0.2e1 + (-t49 * t768 - (-t49 * t301 + t701) * t465 + (t430 * t771 + t540 * t637 - t408 + t536 - t999) * t9 + (t597 + t823) * t10 + (t741 + t826 * t521 + t957 + (-t1028 * t784 + (-pkin(2) * qJD(2) - qJD(5)) * t539) * t543 + t637 * t789 + t1025) * t54 + (t500 - t830 + t771 * t198 + (t770 * t784 - t970) * t540 + (t539 * t770 + t755) * t787 + t959 + t54) * t53) * m(6) + (t87 + t88) * t729 - t951 * t930 + (t42 * (t536 + t691) + t41 * (t708 + t244) + (t42 * t608 + (t637 - t911) * t1010) * t540 + (t500 - t692 + t608 * t787 + (t784 * t918 - t970) * t540) * t91 + (t91 - t972 - t715 + t767 + t1043) * t92) * m(4) + t941 + ((-t1013 + t650 + t74) * t465 + (t1016 * t431 - t1053 * t432 + t1019 - t1030 - t1065 + t971) * t464 + t987 + t1031) * t928 + (t1022 + t986) * t926 + t906 / 0.2e1 - t907 / 0.2e1 - t903 / 0.2e1 + t904 / 0.2e1 - t905 / 0.2e1 + (t179 + ((t142 - t341 + (t366 + t852) * t543 + t813) * t543 + t812 * t540) * qJD(2)) * t727 + t950 * t929 + t1023 * t927 + t895 / 0.2e1 + t902 / 0.2e1 + t893 / 0.2e1 + t894 / 0.2e1 + t891 / 0.2e1 + t892 / 0.2e1 + t890 / 0.2e1 + (-qJD(2) * t639 + t467 * t542 + t468 * t539) * qJD(1); ((-t783 * t848 - t790) * t543 + (t565 + (t543 * t847 + t549) * qJD(2)) * t540) * t727 + ((-t785 * t847 + t790) * t540 + (t565 + (t540 * t848 + t549) * qJD(2)) * t543) * t730 + qJD(1) * (t87 * t540 - t86 * t543 + (t193 * t540 + t194 * t543) * qJD(1)) / 0.2e1 - qJD(1) * ((t539 * t798 + t542 * t797) * qJD(1) + (t586 * t539 + t542 * t936) * qJD(2)) / 0.2e1 + (t7 * t710 + t49 * t635 + (t9 * t709 + t7 * t823 + t49 * t836 + (t54 * t709 + t724) * qJD(1)) * t543 + (t10 * t709 + t54 * t711 + (-t723 + t49 * (-t460 - t761)) * qJD(1) + t589) * t540 - (t539 * t582 + t542 * t548) * qJD(3) - (t665 - t775 + (-t358 + t814) * t465 + t815 * t464) * t49 - (-t427 + t581 * t540 + t814 * t521 + (-t455 + t805) * t464 + t963) * t54 + (-t791 - (t357 - t815) * t521 - t805 * t465 + (t711 - t581) * t543 + t998) * t53) * m(6) + (t8 * t710 + t62 * t635 + (t24 * t757 + t8 * t243 + t62 * t119 + (t64 * t757 + t897) * qJD(1)) * t543 + (t23 * t757 + t64 * t763 + (t896 + t62 * (-t460 - t822)) * qJD(1) + t696) * t540 - t64 * (-t339 * t521 - t505 * t540 + t697 + (-t389 - t455) * t464 + t963) - t62 * (-t336 * t464 + t665 + (-t339 - t358) * t465) - (t539 * t612 + t542 * t550) * qJD(3) + (t389 * t465 - t611 - (t336 + t357) * t521 + (t763 + t505) * t543 + t998) * t63) * m(5) + (-t91 * (t337 * t521 - t390 * t465 + t611) - t92 * (-t340 * t521 - t390 * t464 + t697) - ((-t238 * t91 + t244 * t92) * t539 + t554 * t542) * qJD(3) + t91 * t463 + t31 * t801 + (t31 * t244 + t91 * t820 + (t1010 + t42) * t806) * t543 + (t91 * t386 * qJD(1) + t31 * t238 + t41 * t806 + t92 * t820) * t540 + (t337 * t464 + t340 * t465 - t804 + t759 + (t238 * qJD(1) + t120) * t543 + (t123 + (-t244 - t460) * qJD(1)) * t540) * t79) * m(4) + (-(t182 * t454 - t856) * qJD(1) - (t168 * (-t454 * t540 - t456 * t543) + t653 * t693) * qJD(2) + 0.2e1 * t168 * (t265 * t543 + t266 * t540 + (t387 * t543 - t391 * t540) * qJD(1)) + t653 * t469 + (-t133 * t540 - t134 * t543 + (-t183 * t543 + t858) * qJD(1)) * t485) * m(3) + t974 * t930 + t975 * t929 + ((t1019 * t842 + t950 * t539) * qJD(3) + ((qJD(3) * t1041 + t977) * t542 + t978) * t543 + (t431 * t980 + t432 * t979) * t521 + (-t431 * t984 + t432 * t982) * t465 + (t431 * t983 - t432 * t981) * t464) * t928 + (qJD(1) * t943 + t540 * t992 - t543 * t993) * t927 + (qJD(1) * t942 + t540 * t990 - t543 * t991) * t926 + ((t1042 * t839 - t951 * t539) * qJD(3) + ((qJD(3) * t1020 + t977) * t542 + t978) * t540 + (t429 * t980 + t430 * t979) * t521 + (-t429 * t984 + t430 * t982) * t465 + (t429 * t983 - t430 * t981) * t464) * t925 + ((qJD(3) * t944 - t1032) * t542 + ((t538 * t980 + t541 * t979 - t1048) * t521 + (-t538 * t984 + t541 * t982 - t1056) * t465 + (t538 * t983 - t541 * t981 - t1055) * t464 + t949 * qJD(3)) * t539) * t924 + (qJD(1) * t944 + t540 * t988 - t543 * t989) * t923 - t985 * t782 / 0.2e1 + t976 * t698 + (qJD(1) * t88 + (t540 ^ 2 * t940 - t939 * t1027 + (t143 * t540 + t144 * t543) * qJD(1)) * t1012 + t995) * t540 / 0.2e1 - (qJD(1) * t89 + (-t940 * t1027 + t543 ^ 2 * t939 + (t141 * t540 + t142 * t543) * qJD(1)) * t1012 + t994) * t543 / 0.2e1 + (t65 + t604 + t987) * t789 / 0.2e1 + (t603 + t66 + t986) * t732 - (t540 * t987 + t543 * t986) * t779 / 0.2e1; (-t53 * (-t450 * t465 + t800 + (-t297 - t298) * t521) - t54 * (t304 * t521 + t828 + (-t449 - t450) * t464) - t49 * (t298 * t464 + t831 + (t303 + t304) * t465) - (t53 * t626 + t54 * t625 - t908) * pkin(4) + t9 * t818 + t53 * t714 + t54 * t819 + t7 * t261 + t49 * t885 + (qJD(2) * t548 - t10 * t761 + t53 * t835 - t54 * t769 + t826 * t9) * t542 + (t582 * qJD(2) + (qJD(1) * t583 + t10 * t758 + t54 * t764 + t589) * t543 + (-t9 * t1001 + t53 * t821 - t7 * t761 - t49 * t769 + (-t1001 * t54 + t49 * t762) * qJD(1)) * t540) * t539) * m(6) + (t24 * t818 + t63 * t714 + t64 * t819 + t8 * t261 + t62 * t885 + (qJD(2) * t550 + t63 * t122 - t23 * t822 + t24 * t237 - t64 * t837) * t542 + (t612 * qJD(2) + (qJD(1) * t613 + t23 * t807 + t64 * t827 + t696) * t543 + (t24 * t385 + t63 * t263 - t8 * t822 - t62 * t837 + (t64 * t385 + t62 * t825) * qJD(1)) * t540) * t539 - t63 * (-t451 * t465 + (-t297 - t299) * t521 + t800) - t64 * (t305 * t521 + (-t449 - t451) * t464 + t828) - t62 * (t299 * t464 + (t303 + t305) * t465 + t831)) * m(5) + ((qJD(2) * t554 - t92 * t120 + t91 * t123 + t42 * t238 - t41 * t244) * t542 + (t91 * (-qJD(2) * t238 + t264 * t540) + t92 * (qJD(2) * t244 - t264 * t543) + t31 * t645 + t79 * (-t120 * t540 + t123 * t543 - t238 * t789 - t244 * t787) + (-t41 * t543 + t42 * t540 + (t540 * t92 + t543 * t91) * qJD(1)) * t386) * t539 - t91 * (-t300 * t521 - t452 * t465) - t92 * (t306 * t521 - t452 * t464) - t79 * (t300 * t464 + t306 * t465)) * m(4) + (t539 * t942 + t542 * t951) * t930 + (t539 * t943 - t542 * t950) * t929 + (t431 * t947 + t432 * t946 - t543 * t945) * t928 + ((qJD(2) * t943 - t1023) * t542 + (-qJD(1) * t975 + t950 * qJD(2) + t540 * t993 + t543 * t992) * t539) * t927 + ((qJD(2) * t942 - t1022) * t542 + (-qJD(1) * t974 - t951 * qJD(2) + t540 * t991 + t543 * t990) * t539) * t926 + (t429 * t947 + t430 * t946 - t540 * t945) * t925 + (t1014 * t542 + (t538 * t947 + t541 * t946) * t539) * t924 + ((qJD(2) * t944 - t1021) * t542 + (-qJD(1) * t976 + t949 * qJD(2) + t540 * t989 + t543 * t988) * t539) * t923 - (t894 + t895 + t906 - t907 + t892 + t893 + t904 - t905 + t890 + t891 + t902 - t903 + t941) * t542 / 0.2e1 + t994 * t845 / 0.2e1 + t995 * t843 / 0.2e1 + t985 * t731 + (t539 * t944 - t542 * t949) * t698 + (t539 * t732 + t542 * t729) * t987 + (-t753 / 0.2e1 + t542 * t727) * t986; (t10 * t429 + t431 * t9 + t49 * t598 + t53 * t965 + t54 * t964 + t7 * t846 - t908) * m(6) + (t8 * t846 + t23 * t429 + t24 * t431 + t964 * t64 + t965 * t63 + (-t148 + t598) * t62) * m(5); 0.2e1 * ((-t53 * t783 - t54 * t785 + t7) * t931 + (-t464 * t54 - t465 * t53) * t932) * t542 + 0.2e1 * ((-qJD(2) * t49 - t10 * t540 + t53 * t789 - t54 * t787 - t543 * t9) * t931 + (t49 * (-t464 * t540 - t465 * t543) + (t53 * t540 - t54 * t543) * t521) * t932) * t539;];
tauc = t1(:);
