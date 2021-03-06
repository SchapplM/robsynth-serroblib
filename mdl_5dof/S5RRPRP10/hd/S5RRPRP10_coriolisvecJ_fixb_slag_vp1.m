% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RRPRP10_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP10_coriolisvecJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRP10_coriolisvecJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRP10_coriolisvecJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:09:25
% EndTime: 2019-12-31 20:10:54
% DurationCPUTime: 79.18s
% Computational Cost: add. (15835->1088), mult. (41231->1385), div. (0->0), fcn. (38059->6), ass. (0->529)
t1090 = Icges(4,4) - Icges(3,5);
t1089 = Icges(5,4) + Icges(6,4);
t1088 = Icges(4,5) - Icges(3,6);
t1087 = Icges(4,1) + Icges(3,3);
t1086 = Icges(5,1) + Icges(6,1);
t1023 = Icges(5,5) + Icges(6,5);
t1022 = Icges(5,2) + Icges(6,2);
t1083 = Icges(5,6) + Icges(6,6);
t482 = sin(qJ(2));
t485 = cos(qJ(2));
t1020 = t1088 * t482 - t1090 * t485;
t484 = cos(qJ(4));
t1085 = t1089 * t484;
t481 = sin(qJ(4));
t1084 = t1089 * t481;
t1082 = Icges(5,3) + Icges(6,3);
t1081 = t1023 * t481 + t1083 * t484;
t1080 = t1022 * t484 + t1084;
t1079 = t1086 * t481 + t1085;
t486 = cos(qJ(1));
t1078 = t1087 * t486;
t483 = sin(qJ(1));
t790 = t483 * t484;
t341 = t481 * t486 + t482 * t790;
t792 = t482 * t483;
t702 = t481 * t792;
t788 = t484 * t486;
t342 = -t702 + t788;
t789 = t483 * t485;
t165 = -Icges(6,5) * t342 + Icges(6,6) * t341 + Icges(6,3) * t789;
t168 = -Icges(5,5) * t342 + Icges(5,6) * t341 + Icges(5,3) * t789;
t1005 = t165 + t168;
t339 = -t481 * t483 + t482 * t788;
t791 = t482 * t486;
t703 = t481 * t791;
t340 = t703 + t790;
t787 = t485 * t486;
t315 = Icges(6,4) * t341;
t176 = Icges(6,1) * t342 - Icges(6,5) * t789 - t315;
t318 = Icges(5,4) * t341;
t179 = Icges(5,1) * t342 - Icges(5,5) * t789 - t318;
t973 = t176 + t179;
t316 = Icges(6,4) * t342;
t171 = Icges(6,2) * t341 + Icges(6,6) * t789 - t316;
t319 = Icges(5,4) * t342;
t174 = Icges(5,2) * t341 + Icges(5,6) * t789 - t319;
t974 = t171 + t174;
t1013 = t1005 * t787 + t339 * t974 - t340 * t973;
t314 = Icges(6,4) * t339;
t175 = Icges(6,1) * t340 + Icges(6,5) * t787 + t314;
t317 = Icges(5,4) * t339;
t178 = Icges(5,1) * t340 + Icges(5,5) * t787 + t317;
t1039 = t175 + t178;
t823 = Icges(6,4) * t340;
t169 = Icges(6,2) * t339 + Icges(6,6) * t787 + t823;
t826 = Icges(5,4) * t340;
t172 = Icges(5,2) * t339 + Icges(5,6) * t787 + t826;
t1040 = t169 + t172;
t163 = Icges(6,5) * t340 + Icges(6,6) * t339 + Icges(6,3) * t787;
t166 = Icges(5,5) * t340 + Icges(5,6) * t339 + Icges(5,3) * t787;
t1041 = t163 + t166;
t1014 = t1039 * t340 + t1040 * t339 + t1041 * t787;
t717 = qJD(4) * t485;
t724 = qJD(2) * t483;
t377 = t486 * t717 + t724;
t722 = qJD(2) * t486;
t378 = -t483 * t717 + t722;
t719 = qJD(4) * t482;
t455 = qJD(1) + t719;
t1001 = -t1023 * t482 + t1079 * t485;
t1002 = t1080 * t485 - t1083 * t482;
t1003 = t1081 * t485 - t1082 * t482;
t916 = -t1001 * t340 - t1002 * t339 - t1003 * t787;
t947 = -t1013 * t378 + t1014 * t377 + t916 * t455;
t1011 = t1005 * t789 + t341 * t974 + t342 * t973;
t1012 = -t1039 * t342 + t1040 * t341 + t1041 * t789;
t908 = t1001 * t342 - t1002 * t341 - t1003 * t789;
t946 = -t1011 * t378 + t1012 * t377 + t455 * t908;
t1066 = -t1088 * t792 + t1090 * t789 + t1078;
t1035 = t1020 * t486 + t1087 * t483;
t829 = Icges(3,4) * t482;
t408 = Icges(3,1) * t485 - t829;
t814 = Icges(4,6) * t482;
t592 = Icges(4,2) * t485 - t814;
t1077 = t408 + t592;
t813 = Icges(4,6) * t485;
t590 = -Icges(4,3) * t482 + t813;
t470 = Icges(3,4) * t485;
t599 = -Icges(3,2) * t482 + t470;
t1076 = t590 + t599;
t401 = Icges(3,5) * t482 + Icges(3,6) * t485;
t597 = Icges(4,4) * t482 + Icges(4,5) * t485;
t1056 = -t597 + t401;
t815 = Icges(3,6) * t486;
t280 = Icges(3,4) * t789 - Icges(3,2) * t792 - t815;
t443 = Icges(4,6) * t792;
t827 = Icges(4,4) * t486;
t291 = Icges(4,2) * t789 - t443 + t827;
t1075 = t280 * t482 - t291 * t485;
t287 = Icges(3,5) * t483 + t408 * t486;
t288 = Icges(4,5) * t483 - t486 * t590;
t1068 = -t287 * t789 - t288 * t792;
t1067 = t1081 * t482 + t1082 * t485;
t933 = t1080 * t482 + t1083 * t485;
t932 = t1023 * t485 + t1079 * t482;
t1065 = (-t1023 * t484 + t1083 * t481) * t485;
t1064 = (t1022 * t481 - t1085) * t485;
t1063 = (-t1086 * t484 + t1084) * t485;
t405 = Icges(3,2) * t485 + t829;
t589 = Icges(4,3) * t485 + t814;
t1059 = -t405 - t589;
t407 = Icges(3,1) * t482 + t470;
t591 = Icges(4,2) * t482 + t813;
t1062 = t407 + t591;
t574 = t405 * t482 - t407 * t485;
t1057 = -t482 * t589 + t485 * t591 - t574;
t448 = Icges(3,4) * t792;
t820 = Icges(3,5) * t486;
t286 = Icges(3,1) * t789 - t448 - t820;
t819 = Icges(4,5) * t486;
t289 = Icges(4,6) * t789 - Icges(4,3) * t792 + t819;
t1017 = -t286 * t485 + t289 * t482 + t1075;
t725 = qJD(2) * t482;
t723 = qJD(2) * t485;
t1061 = t1076 * qJD(2);
t1060 = t1077 * qJD(2);
t1055 = t1035 * t486 + t1068;
t965 = t1035 * t483 + t287 * t787 + t288 * t791;
t1054 = t1066 * t483 - t286 * t787 + t289 * t791;
t918 = t1056 * t483;
t281 = Icges(3,6) * t483 + t486 * t599;
t1033 = t281 - t288;
t1051 = t1033 * t485;
t1050 = -t1059 * t723 + t1062 * t725;
t648 = qJD(1) * t482 + qJD(4);
t570 = t648 * t486;
t571 = t455 * t481;
t152 = t484 * t570 + (t484 * t723 - t571) * t483;
t572 = t484 * t455;
t689 = t483 * t723;
t153 = t483 * t572 + (t570 + t689) * t481;
t691 = t482 * t724;
t726 = qJD(1) * t486;
t694 = t485 * t726;
t540 = -t691 + t694;
t78 = Icges(6,5) * t153 + Icges(6,6) * t152 + Icges(6,3) * t540;
t80 = Icges(5,5) * t153 + Icges(5,6) * t152 + Icges(5,3) * t540;
t1049 = t78 + t80;
t688 = t485 * t722;
t523 = -t483 * t648 + t688;
t154 = t484 * t523 - t486 * t571;
t155 = t481 * t523 + t486 * t572;
t690 = t482 * t722;
t727 = qJD(1) * t485;
t695 = t483 * t727;
t921 = t690 + t695;
t79 = Icges(6,5) * t155 + Icges(6,6) * t154 - Icges(6,3) * t921;
t81 = Icges(5,5) * t155 + Icges(5,6) * t154 - Icges(5,3) * t921;
t1048 = t79 + t81;
t82 = Icges(6,4) * t153 + Icges(6,2) * t152 + Icges(6,6) * t540;
t84 = Icges(5,4) * t153 + Icges(5,2) * t152 + Icges(5,6) * t540;
t1047 = t82 + t84;
t83 = Icges(6,4) * t155 + Icges(6,2) * t154 - Icges(6,6) * t921;
t85 = Icges(5,4) * t155 + Icges(5,2) * t154 - Icges(5,6) * t921;
t1046 = t83 + t85;
t86 = Icges(6,1) * t153 + Icges(6,4) * t152 + Icges(6,5) * t540;
t88 = Icges(5,1) * t153 + Icges(5,4) * t152 + Icges(5,5) * t540;
t1045 = t86 + t88;
t87 = Icges(6,1) * t155 + Icges(6,4) * t154 - Icges(6,5) * t921;
t89 = Icges(5,1) * t155 + Icges(5,4) * t154 - Icges(5,5) * t921;
t1044 = t87 + t89;
t1043 = -t1017 * t483 + t1066 * t486;
t444 = Icges(4,6) * t791;
t828 = Icges(4,4) * t483;
t290 = -Icges(4,2) * t787 + t444 + t828;
t979 = -t281 * t792 - t290 * t789 - t1055;
t978 = -t280 * t791 + t291 * t787 - t1054;
t977 = -t281 * t791 - t290 * t787 + t965;
t1042 = t1057 * t486 + t918;
t1038 = t1067 * qJD(2) + t1065 * qJD(4);
t1037 = t933 * qJD(2) + t1064 * qJD(4);
t1036 = t932 * qJD(2) + t1063 * qJD(4);
t1034 = t280 + t289;
t1032 = t286 + t291;
t1031 = t287 - t290;
t1027 = -t1001 * t481 - t1002 * t484;
t1025 = t281 * t482 + t290 * t485;
t1019 = t1060 * t485 - t1061 * t482 + (t1059 * t485 - t1062 * t482) * qJD(2) + t1056 * qJD(1);
t1018 = t1056 * qJD(2);
t1016 = t287 * t485 + t288 * t482 - t1025;
t1015 = t1057 * qJD(1) - qJD(2) * t1020;
t953 = t1039 * t153 + t1040 * t152 + t1041 * t540 - t1044 * t342 + t1046 * t341 + t1048 * t789;
t952 = t1005 * t540 - t1045 * t342 + t1047 * t341 + t1049 * t789 + t152 * t974 - t153 * t973;
t951 = t1039 * t155 + t1040 * t154 - t1041 * t921 + t1044 * t340 + t1046 * t339 + t1048 * t787;
t950 = -t1005 * t921 + t1045 * t340 + t1047 * t339 + t1049 * t787 + t154 * t974 - t155 * t973;
t992 = -t1001 * t153 - t1002 * t152 - t1003 * t540 - t1036 * t342 + t1037 * t341 + t1038 * t789;
t991 = -t1001 * t155 - t1002 * t154 + t1003 * t921 + t1036 * t340 + t1037 * t339 + t1038 * t787;
t587 = t169 * t484 + t175 * t481;
t66 = t163 * t482 - t485 * t587;
t585 = t172 * t484 + t178 * t481;
t68 = t166 * t482 - t485 * t585;
t1010 = t66 + t68;
t586 = t171 * t484 - t176 * t481;
t944 = t165 * t482;
t67 = -t485 * t586 + t944;
t584 = t174 * t484 - t179 * t481;
t941 = t168 * t482;
t69 = -t485 * t584 + t941;
t1009 = t67 + t69;
t1008 = t486 * (-t1018 * t483 + (t1017 + t1035) * qJD(1));
t907 = -t1003 * t482 - t1027 * t485;
t1007 = t1042 * qJD(1);
t1004 = (t1027 + t1067) * t455 + (-t1003 * t483 - t584 - t586) * t378 + (t1003 * t486 + t585 + t587) * t377;
t1000 = -t1033 * t483 + t1034 * t486;
t999 = (-t1043 * t486 + t483 * t979) * qJD(2);
t998 = (t483 * t977 - t486 * t978) * qJD(2);
t997 = t1076 + t1062;
t996 = t1077 + t1059;
t995 = (t443 + t448 + (Icges(3,2) + Icges(4,3)) * t789 - t1032) * t486 + (-Icges(4,3) * t787 - t405 * t486 + t1031 - t444) * t483;
t901 = t377 * (-t1086 * t339 + t1040 + t823 + t826) - t378 * (-t1086 * t341 - t316 - t319 + t974) + t455 * (-t1002 - t1063);
t902 = t377 * (-t1022 * t340 + t1039 + t314 + t317) - t378 * (t1022 * t342 + t315 + t318 - t973) + t455 * (-t1001 + t1064);
t990 = (-t1037 * t484 - t1036 * t481 + (t1001 * t484 - t1002 * t481) * qJD(4) - t1003 * qJD(2)) * t485 + (qJD(2) * t1027 + t1038) * t482;
t354 = t597 * t486;
t713 = (t589 * t792 - t591 * t789 - t354) * qJD(1);
t794 = t401 * t486;
t156 = -t483 * t574 - t794;
t714 = t156 * qJD(1);
t989 = t714 - t713 + t999;
t988 = t998 + t1007;
t983 = t1017 * qJD(2) + t1050 * t483 + (-t1051 + (-t486 * t592 - t287 + t828) * t482) * qJD(1);
t982 = t1016 * qJD(2) - t1050 * t486 + ((t815 - t819) * t485 + (t820 - t827) * t482 + (-t1076 * t485 - t1077 * t482) * t483) * qJD(1);
t981 = -t1015 * t483 + t1019 * t486;
t980 = t1015 * t486 + t1019 * t483;
t976 = t1032 * t482 + t1034 * t485;
t975 = t1031 * t482 + t1051;
t970 = t1065 * t455 + (-t1023 * t341 - t1083 * t342) * t378 + (t1023 * t339 - t1083 * t340) * t377;
t460 = pkin(4) * t484 + pkin(3);
t480 = -qJ(5) - pkin(7);
t969 = rSges(6,1) * t342 - rSges(6,2) * t341 + t486 * t460 + t480 * t789;
t968 = t1025 + t1066;
t960 = -t1018 * t486 + (-t1020 * t483 - t1016 + t1078) * qJD(1);
t715 = qJD(5) * t486;
t436 = t485 * t715;
t456 = pkin(7) * t789;
t390 = pkin(3) * t486 - t456;
t774 = rSges(6,3) * t789 + pkin(4) * t702 + t390 - t969;
t958 = t455 * t774 - t436;
t957 = 0.2e1 * qJD(2);
t710 = qJD(2) * qJD(4);
t675 = t482 * t710;
t261 = qJD(1) * t377 - t483 * t675;
t262 = qJD(1) * t378 - t486 * t675;
t674 = t485 * t710;
t956 = t1011 * t261 + t1012 * t262 + t377 * t953 - t952 * t378 + t992 * t455 + t908 * t674;
t955 = t1013 * t261 + t1014 * t262 + t377 * t951 - t378 * t950 + t455 * t991 + t674 * t916;
t954 = rSges(6,1) + pkin(4);
t17 = (qJD(2) * t587 + t79) * t482 + (qJD(2) * t163 - t481 * t87 - t484 * t83 + (t169 * t481 - t175 * t484) * qJD(4)) * t485;
t19 = (qJD(2) * t585 + t81) * t482 + (qJD(2) * t166 - t481 * t89 - t484 * t85 + (t172 * t481 - t178 * t484) * qJD(4)) * t485;
t949 = t17 + t19;
t18 = (qJD(2) * t586 + t78) * t482 + (qJD(2) * t165 - t481 * t86 - t484 * t82 + (t171 * t481 + t176 * t484) * qJD(4)) * t485;
t20 = (qJD(2) * t584 + t80) * t482 + (qJD(2) * t168 - t481 * t88 - t484 * t84 + (t174 * t481 + t179 * t484) * qJD(4)) * t485;
t948 = t18 + t20;
t945 = -t1009 * t378 + t1010 * t377 + t455 * t907;
t937 = t1002 * t483;
t936 = t1002 * t486;
t935 = t1001 * t483;
t934 = t1001 * t486;
t863 = pkin(4) * t481;
t708 = t482 * t863;
t859 = pkin(7) + t480;
t312 = -t485 * t859 + t708;
t618 = rSges(6,1) * t481 + rSges(6,2) * t484;
t849 = rSges(6,3) * t485;
t931 = t482 * t618 + t312 + t849;
t930 = t1004 * t485;
t929 = t1000 * t485 - t482 * t995;
t928 = (-t482 * t997 + t485 * t996) * qJD(1);
t862 = pkin(4) * t485;
t927 = t481 * t862 + t485 * t618 + (-rSges(6,3) + t859) * t482;
t926 = t1003 * t455 + t1005 * t378 - t1041 * t377;
t898 = t1009 * t483 + t1010 * t486;
t925 = -t1009 * t486 + t1010 * t483;
t897 = t1011 * t483 + t1012 * t486;
t924 = -t1011 * t486 + t1012 * t483;
t896 = t1013 * t483 + t1014 * t486;
t923 = -t1013 * t486 + t1014 * t483;
t465 = qJD(5) * t482;
t457 = pkin(7) * t787;
t389 = t483 * pkin(3) + t457;
t793 = t482 * qJ(3);
t414 = pkin(2) * t485 + t793;
t366 = t414 * t483;
t372 = pkin(2) * t787 + qJ(3) * t791;
t721 = qJD(3) * t485;
t627 = t366 * t724 + t372 * t722 - t721;
t555 = t389 * t722 - t390 * t724 + t627;
t906 = t340 * rSges(6,1) + t339 * rSges(6,2) + rSges(6,3) * t787 + pkin(4) * t703 + t483 * t460;
t775 = -t480 * t787 - t389 + t906;
t33 = t377 * t774 + t378 * t775 + t465 + t555;
t720 = qJD(3) * t486;
t439 = t482 * t720;
t410 = pkin(2) * t482 - qJ(3) * t485;
t662 = -pkin(7) * t482 - t410;
t626 = t662 * t486;
t538 = qJD(2) * t626 + t439;
t476 = t486 * pkin(6);
t419 = pkin(1) * t483 - t476;
t747 = -t366 - t419;
t508 = (t390 + t747) * qJD(1) + t538;
t47 = t378 * t927 + t508 - t958;
t660 = t47 * t927;
t922 = -t33 * t775 - t660;
t920 = t1020 * qJD(1);
t919 = t794 - t354;
t622 = rSges(5,1) * t342 - rSges(5,2) * t341;
t188 = rSges(5,3) * t789 - t622;
t621 = rSges(5,1) * t481 + rSges(5,2) * t484;
t525 = -rSges(5,3) * t482 + t485 * t621;
t917 = -t188 * t455 + t378 * t525;
t915 = t33 * t774;
t373 = t410 * t724;
t433 = pkin(7) * t691;
t420 = t486 * pkin(1) + t483 * pkin(6);
t642 = t372 + t420;
t903 = t389 + t642;
t542 = qJD(1) * t903 - t373 - t433;
t466 = qJD(3) * t482;
t716 = qJD(5) * t485;
t910 = t483 * (t716 + t466);
t48 = t377 * t927 + t455 * t775 + t542 + t910;
t913 = t48 * t927;
t855 = rSges(4,2) * t482;
t617 = rSges(4,3) * t485 + t855;
t909 = (-qJD(2) * t617 - t466) * t483;
t663 = -qJ(3) - t863;
t472 = t483 * rSges(4,1);
t307 = -rSges(4,2) * t787 + rSges(4,3) * t791 + t472;
t905 = t307 + t642;
t391 = qJD(1) * t419;
t904 = -qJD(1) * t366 - t391;
t900 = t970 * t485;
t718 = qJD(4) * t484;
t682 = t482 * t718;
t899 = t486 * pkin(4) * t682 + t155 * rSges(6,1) + t154 * rSges(6,2) + t460 * t726 + t480 * t921 + t688 * t863 + t436;
t709 = -rSges(5,3) - pkin(2) - pkin(7);
t893 = t709 * t485 - pkin(1) - t793;
t892 = t455 * t990 + t674 * t907;
t64 = t508 + t917;
t184 = t340 * rSges(5,1) + t339 * rSges(5,2) + rSges(5,3) * t787;
t685 = t483 * t466;
t65 = t184 * t455 + t377 * t525 + t542 + t685;
t891 = t65 * t483 + t64 * t486;
t267 = qJD(1) * t389 - t433;
t464 = pkin(3) * t726;
t268 = -pkin(7) * t921 + t464;
t728 = qJD(1) * t483;
t696 = t482 * t728;
t736 = qJ(3) * t688 + t439;
t181 = -pkin(2) * t921 - qJ(3) * t696 + t736;
t337 = t482 * t726 + t689;
t434 = pkin(2) * t691;
t182 = pkin(2) * t694 + qJ(3) * t337 - t434 + t685;
t712 = qJD(1) * qJD(2);
t677 = t486 * t712;
t711 = qJD(2) * qJD(3);
t647 = t181 * t722 + t182 * t724 + t366 * t677 + t482 * t711;
t544 = t267 * t724 + t268 * t722 - t390 * t677 + t647;
t746 = t372 + t389;
t629 = t746 * t728;
t857 = -rSges(6,3) * t921 + pkin(7) * t690 - t464 + (pkin(7) * t727 - t648 * t863) * t483 + t899;
t620 = rSges(6,1) * t153 + rSges(6,2) * t152;
t706 = qJD(4) * t863;
t858 = rSges(6,3) * t540 + t620 + t433 + (qJD(1) * t312 + t706) * t486 + (t480 * t725 + t716 + (-pkin(3) + t460) * qJD(1) + (t481 * t723 + t682) * pkin(4)) * t483;
t5 = t857 * t378 + t858 * t377 + t774 * t262 - t775 * t261 + (-t629 + t716) * qJD(2) + t544;
t882 = t33 * t857 + t5 * t775;
t879 = m(4) / 0.2e1;
t878 = m(5) / 0.2e1;
t877 = -m(6) / 0.2e1;
t876 = m(6) / 0.2e1;
t875 = t261 / 0.2e1;
t874 = t262 / 0.2e1;
t873 = -t377 / 0.2e1;
t872 = t377 / 0.2e1;
t871 = -t378 / 0.2e1;
t870 = t378 / 0.2e1;
t869 = -t455 / 0.2e1;
t868 = t455 / 0.2e1;
t864 = -rSges(6,3) - pkin(2);
t856 = rSges(3,1) * t485;
t854 = rSges(4,2) * t485;
t853 = rSges(4,3) * t482;
t851 = rSges(5,3) * t485;
t848 = t17 * t377;
t847 = t18 * t378;
t846 = t19 * t377;
t845 = t20 * t378;
t471 = t483 * rSges(3,3);
t838 = t66 * t262;
t837 = t67 * t261;
t836 = t68 * t262;
t835 = t69 * t261;
t834 = -rSges(4,3) - qJ(3);
t734 = rSges(3,2) * t792 + t486 * rSges(3,3);
t303 = rSges(3,1) * t789 - t734;
t412 = rSges(3,1) * t482 + rSges(3,2) * t485;
t692 = t412 * t722;
t140 = -t692 + (-t303 - t419) * qJD(1);
t805 = t140 * t483;
t804 = t140 * t486;
t693 = t412 * t724;
t306 = rSges(3,1) * t787 - rSges(3,2) * t791 + t471;
t753 = t306 + t420;
t141 = qJD(1) * t753 - t693;
t371 = t412 * t486;
t803 = t141 * t371;
t785 = t155 * rSges(5,1) + t154 * rSges(5,2);
t388 = t420 * qJD(1);
t776 = -t182 - t388;
t367 = (-rSges(6,1) * t484 + rSges(6,2) * t481) * t485;
t705 = pkin(4) * t718;
t773 = qJD(2) * t931 + qJD(4) * t367 - t485 * t705 + t465;
t772 = -rSges(6,2) * t340 + t339 * t954;
t771 = rSges(6,2) * t342 + t341 * t954;
t770 = t927 * t483;
t752 = -t307 - t372;
t751 = t483 * t366 + t486 * t372;
t320 = qJD(2) * t414 - t721;
t415 = t853 - t854;
t750 = -t415 * qJD(2) - t320;
t369 = t410 * t486;
t749 = -qJD(1) * t369 + t483 * t721;
t676 = t485 * t711;
t748 = qJD(1) * t373 + t486 * t676;
t374 = t410 * t728;
t435 = pkin(7) * t696;
t745 = t374 + t435;
t740 = -t410 + t617;
t739 = -t414 - t415;
t735 = rSges(3,2) * t696 + rSges(3,3) * t726;
t733 = t483 ^ 2 + t486 ^ 2;
t704 = t48 * t726;
t701 = t65 * t726;
t700 = t486 * t181 + t483 * t182 + t366 * t726;
t699 = -t267 + t776;
t363 = t410 * t483;
t698 = -t363 * t724 - t369 * t722 + t466;
t463 = pkin(6) * t726;
t697 = t463 + t736;
t684 = t184 * t717;
t683 = t188 * t717;
t678 = -pkin(1) - t856;
t672 = t726 / 0.2e1;
t671 = -t724 / 0.2e1;
t669 = t723 / 0.2e1;
t668 = -t722 / 0.2e1;
t667 = t722 / 0.2e1;
t658 = rSges(4,1) * t486 - rSges(4,3) * t792;
t657 = t486 * t740;
t652 = qJD(2) * t733;
t651 = qJD(1) * t363 + t485 * t720;
t650 = -t320 - t465;
t375 = qJD(1) * (-pkin(1) * t728 + t463);
t646 = t483 * t676 + t375 + (t181 + t439) * qJD(1);
t644 = t486 * t389 - t483 * t390 + t751;
t643 = rSges(4,1) * t726 + rSges(4,2) * t921 + rSges(4,3) * t688;
t635 = t484 * t862 - t367;
t634 = t525 + t662;
t631 = qJD(4) * t669;
t630 = -pkin(7) * t723 - t320;
t628 = t435 + t651;
t624 = -rSges(3,2) * t482 + t856;
t623 = rSges(5,1) * t153 + rSges(5,2) * t152;
t616 = -t47 * t486 - t48 * t483;
t603 = qJD(1) * t268 + t646;
t588 = -t141 * t483 - t804;
t583 = t184 * t483 - t188 * t486;
t573 = t662 + t927;
t302 = t482 * t621 + t851;
t368 = (-rSges(5,1) * t484 + rSges(5,2) * t481) * t485;
t208 = qJD(2) * t302 + qJD(4) * t368;
t567 = -t208 + t630;
t566 = -pkin(1) - t414;
t565 = qJD(1) * t626;
t564 = t483 * t267 + t486 * t268 - t390 * t726 + t700;
t558 = qJD(2) * t657 + t439;
t365 = t412 * t483;
t364 = t617 * t483;
t557 = t630 - t773;
t131 = (t303 * t483 + t306 * t486) * qJD(2);
t487 = qJD(2) ^ 2;
t541 = qJD(1) * t433 - t457 * t487 + t748;
t537 = t33 * t858 + t5 * t774;
t532 = t47 * t722 + t48 * t724 - t5;
t531 = -t47 * t774 + t48 * t775;
t530 = -t913 - t915;
t514 = qJD(1) * t390 + t538 + t904;
t55 = t184 * t378 + t188 * t377 + t555;
t503 = t55 * t583 - (-t64 * t483 + t486 * t65) * t525;
t6 = -t858 * t455 - t773 * t378 - t927 * t261 + (t486 * t650 - t717 * t774) * qJD(2) + (t699 - t910) * qJD(1) + t541;
t7 = (-pkin(7) * t483 * t487 + qJD(1) * t715) * t485 + t857 * t455 - t773 * t377 + t927 * t262 + (t483 * t650 + t717 * t775 + t565) * qJD(2) + t603;
t496 = (qJD(2) * t33 - t47 * t728 + t483 * t7 + t486 * t6 + t704) * t876;
t495 = -t483 * t922 + t530 * t486;
t386 = t624 * qJD(2);
t370 = t617 * t486;
t338 = t688 - t696;
t336 = t733 * t725;
t308 = rSges(4,2) * t789 + t658;
t249 = t525 * t486;
t247 = t525 * t483;
t231 = rSges(5,1) * t341 + rSges(5,2) * t342;
t229 = rSges(5,1) * t339 - rSges(5,2) * t340;
t212 = -rSges(4,3) * t696 + t643;
t211 = qJD(2) * t364 + (t415 * t486 + t472) * qJD(1);
t210 = -qJD(2) * t365 + (t486 * t624 + t471) * qJD(1);
t209 = -rSges(3,1) * t921 - rSges(3,2) * t688 + t735;
t106 = qJD(1) * t905 - t373 - t909;
t105 = (t308 + t747) * qJD(1) + t558;
t104 = -t386 * t722 + (-t210 - t388 + t693) * qJD(1);
t103 = -t386 * t724 + t375 + (t209 - t692) * qJD(1);
t102 = (t307 * t486 - t308 * t483) * qJD(2) + t627;
t93 = -rSges(5,3) * t921 + t785;
t91 = rSges(5,3) * t540 + t623;
t54 = t750 * t722 + (-t211 + t776 + t909) * qJD(1) + t748;
t53 = qJD(1) * t212 + (qJD(1) * t657 + t750 * t483) * qJD(2) + t646;
t34 = (t211 * t483 + t212 * t486 + (-t308 * t486 + t483 * t752) * qJD(1)) * qJD(2) + t647;
t26 = -t208 * t378 - t261 * t525 - t455 * t91 + (-t320 * t486 - t683) * qJD(2) + (-t685 + t699) * qJD(1) + t541;
t25 = -t487 * t456 - t208 * t377 + t262 * t525 + t455 * t93 + (-t320 * t483 + t565 + t684) * qJD(2) + t603;
t16 = -qJD(2) * t629 - t184 * t261 + t188 * t262 + t377 * t91 + t378 * t93 + t544;
t1 = [(qJD(2) * t1057 + t1060 * t482 + t1061 * t485) * qJD(1) + (0.2e1 * t713 - t714 + ((t486 * t968 - t965 + t977) * t486 + (t483 * t968 + t1055 + t978) * t483) * qJD(2) + t989) * t671 + (-(t514 - t64 + t917) * t65 + t26 * (-t456 + t476 + t622) + t64 * (t433 + t434 - t623) + t25 * (t184 + t903) + t65 * (t464 + t697 + t785) + (qJD(1) * t64 * t893 + t65 * t709 * t725 + t26 * pkin(3)) * t486 + (t26 * (t566 - t851) + t64 * (rSges(5,3) * t725 - qJ(3) * t723 - t466) + (t64 * (-pkin(3) - pkin(6)) + t893 * t65) * qJD(1)) * t483) * m(5) + (-(-qJD(1) * t303 - t140 - t391 - t692) * t141 + t104 * (t483 * t678 + t476 + t734) + t103 * t753 + t141 * (t463 + t735) + (t412 * t805 - t803) * qJD(2) + ((-pkin(1) - t624) * t804 + (t140 * (-rSges(3,3) - pkin(6)) + t141 * t678) * t483) * qJD(1)) * m(3) + (t53 * t905 + (t434 + (-pkin(1) + (rSges(4,2) - pkin(2)) * t485 + t834 * t482) * t726 + (-t466 + (t485 * t834 - t855) * qJD(2) + (-rSges(4,1) - pkin(6)) * qJD(1)) * t483) * t105 + (t476 + t658 + (t566 + t854) * t483) * t54 + (-pkin(2) * t690 + t105 - t558 + t643 + t697 - t904 + (-t308 + (t566 - t853) * t483) * qJD(1)) * t106) * m(4) + (-t941 - t944 + (-t481 * t973 + t484 * t974) * t485 + t1009) * t377 * t869 + t837 / 0.2e1 + t838 / 0.2e1 + t835 / 0.2e1 + t836 / 0.2e1 + ((t965 * t483 + ((t1035 + t1075) * t486 + t979 + t1054 + t1068) * t486) * qJD(2) + t1007) * t667 - t847 / 0.2e1 + t848 / 0.2e1 - t845 / 0.2e1 + t846 / 0.2e1 + t991 * t872 + (t992 + t947) * t871 + t892 + t947 * t870 + (t981 + t982) * t724 / 0.2e1 + (t980 - t983 + t988) * t668 + t916 * t874 + (t6 * (t476 + t969) + t47 * (t434 - t620) + t7 * (t642 + t906) + (-t7 * t480 * t485 + (-t706 + (-pkin(1) + t663 * t482 + (t480 + t864) * t485) * qJD(1)) * t47) * t486 + (-t6 * pkin(1) + (t6 * t864 + t47 * (qJD(2) * t663 - qJD(5))) * t485 + (t6 * t663 + t47 * (-qJD(3) - t705 + (rSges(6,3) - t480) * qJD(2))) * t482 + t47 * (-pkin(6) - t460) * qJD(1)) * t483 - t913 * t378 + (t697 + t899 + t864 * t725 * t486 + (-t706 + (t566 - t708 - t849) * qJD(1)) * t483 + t47 - t514 + t958) * t48) * m(6) + t908 * t875 + ((t975 + t1042) * t486 + (t156 + t976) * t483) * t712 / 0.2e1; t924 * t875 + t923 * t874 + ((-t1013 * t792 + t916 * t485) * qJD(4) + ((-qJD(4) * t1014 + t926) * t482 + t930) * t486 + (t339 * t933 + t340 * t932) * t455 + (-t339 * t937 - t340 * t935) * t378 + (t339 * t936 + t340 * t934) * t377) * t873 + (qJD(1) * t896 + t483 * t951 - t486 * t950) * t872 + (qJD(1) * t897 + t483 * t953 - t486 * t952) * t871 + ((-t1012 * t791 + t908 * t485) * qJD(4) + ((-qJD(4) * t1011 + t926) * t482 + t930) * t483 + (t341 * t933 - t342 * t932) * t455 + (-t341 * t937 + t342 * t935) * t378 + (t341 * t936 - t342 * t934) * t377) * t870 + (((-t481 * t932 - t484 * t933 - t1003) * t455 + (t481 * t935 + t484 * t937 - t1005) * t378 + (-t481 * t934 - t484 * t936 + t1041) * t377 + t907 * qJD(4)) * t485 + (-qJD(4) * t898 + t1004) * t482) * t869 + (qJD(1) * t898 + t483 * t949 - t486 * t948) * t868 - ((t1000 * t482 + t485 * t995) * qJD(2) + (t996 * t482 + t997 * t485) * qJD(1)) * qJD(1) / 0.2e1 + (t983 * t486 + t982 * t483 + (t483 * t976 + t486 * t975) * qJD(1)) * qJD(1) / 0.2e1 + ((-t724 * t919 + t920) * t483 + ((t483 * t918 + t929) * qJD(2) + t928) * t486) * t671 + ((-t722 * t918 - t920) * t486 + ((t486 * t919 + t929) * qJD(2) + t928) * t483) * t667 - t945 * t717 / 0.2e1 + t925 * t631 + (t5 * t644 + t33 * t564 + (t7 * t573 + t48 * t557 + (-t660 + t33 * (-t746 - t775)) * qJD(1) + t537) * t483 - t48 * (-t414 * t724 + t749) - t33 * t698 - (t33 * t770 - t48 * t931) * t377 - (pkin(7) * qJD(2) * t616 + qJD(4) * t531 + t33 * qJD(5)) * t485 - (t616 * qJD(5) + (-t33 * t652 - t704) * pkin(7) + t495 * qJD(4)) * t482 + (t378 * t931 + t414 * t722 + t770 * t455 - t628 + t745) * t47 + (t6 * t573 + (t48 * t573 + t915) * qJD(1) + t882 + (-t33 * t378 - t48 * t455) * t927 + t47 * t557) * t486) * m(6) + (t64 * t745 + t16 * t644 + t55 * t564 + (t26 * t634 + t64 * t567 + t16 * t184 + t55 * t93 + (t55 * t188 + t634 * t65) * qJD(1)) * t486 + (t25 * t634 + t65 * t567 + t16 * t188 + t55 * t91 + (-t64 * t525 + t55 * (-t184 - t746)) * qJD(1)) * t483 - t64 * (-t247 * t455 - t302 * t378 + t628 - t683) - t65 * (t249 * t455 - t302 * t377 + t684 + t749) - t55 * (t247 * t377 + t249 * t378 + t698) - ((-t55 * t652 - t701) * pkin(7) + t503 * qJD(4)) * t482 - t891 * qJD(2) * (-pkin(7) * t485 - t414)) * m(5) + (-t105 * (-qJD(1) * t364 + t651) - t106 * (qJD(1) * t370 + t749) - t102 * t698 - ((t102 * t370 + t105 * t739) * t486 + (t102 * t364 + t106 * t739) * t483) * qJD(2) + t105 * t374 + t34 * t751 + t102 * t700 + (t54 * t740 + t105 * t750 + t34 * t307 + t102 * t212 + (-t102 * t308 + t106 * t740) * qJD(1)) * t486 + (t53 * t740 + t106 * t750 - t34 * t308 + t102 * t211 + (t102 * t752 - t105 * t617) * qJD(1)) * t483) * m(4) + (-(t140 * t365 - t803) * qJD(1) - (t131 * (-t365 * t483 - t371 * t486) + t588 * t624) * qJD(2) + 0.2e1 * t131 * (t209 * t486 + t210 * t483 + (t303 * t486 - t306 * t483) * qJD(1)) + t588 * t386 + (-t103 * t483 - t104 * t486 + (-t141 * t486 + t805) * qJD(1)) * t412) * m(3) + (t981 * qJD(1) + t955 + (t977 * t726 + (t978 * qJD(1) + t960 * t483 - t1008) * t483) * t957) * t483 / 0.2e1 - (t980 * qJD(1) + t956 + ((qJD(1) * t979 + t1008) * t486 + (qJD(1) * t1043 - t960 * t486) * t483) * t957) * t486 / 0.2e1 + (t946 + t989 + t999) * t728 / 0.2e1 + (t947 + t988 + t998) * t672 + (t483 * t946 + t486 * t947) * t719 / 0.2e1; -m(4) * (t102 * t336 + t105 * t338 + t106 * t337) - m(5) * (t336 * t55 + t337 * t65 + t338 * t64) - m(6) * (t33 * t336 + t337 * t48 + t338 * t47) + 0.2e1 * ((t105 * t722 + t106 * t724 - t34) * t879 + (t64 * t722 + t65 * t724 - t16) * t878 + t532 * t876) * t485 + 0.2e1 * ((qJD(2) * t102 - t105 * t728 + t106 * t726 + t483 * t53 + t486 * t54) * t879 + (qJD(2) * t55 + t25 * t483 + t26 * t486 - t64 * t728 + t701) * t878 + t496) * t482; (t482 * t908 + t485 * t897) * t875 + (t482 * t916 + t896 * t485) * t874 + (t902 * t339 - t340 * t901 + t900 * t486) * t873 + ((-qJD(1) * t923 + t916 * qJD(2) + t483 * t950 + t486 * t951) * t485 + (-qJD(2) * t896 + t991) * t482) * t872 + ((-qJD(1) * t924 + t908 * qJD(2) + t483 * t952 + t486 * t953) * t485 + (-qJD(2) * t897 + t992) * t482) * t871 + (t341 * t902 + t342 * t901 + t483 * t900) * t870 + ((t901 * t481 - t484 * t902) * t485 + t970 * t482) * t869 + ((-qJD(1) * t925 + t907 * qJD(2) + t483 * t948 + t486 * t949) * t485 + (-qJD(2) * t898 + t990) * t482) * t868 + (t837 + t838 - t847 + t848 + t835 + t836 - t845 + t846 + t892) * t482 / 0.2e1 + t956 * t789 / 0.2e1 + t955 * t787 / 0.2e1 + t945 * t669 + (t482 * t907 + t485 * t898) * t631 + ((qJD(2) * t495 - t47 * t858 + t48 * t857 - t6 * t774 + t7 * t775) * t482 + (t531 * qJD(2) + (qJD(1) * t922 - t48 * t773 + t7 * t927 + t537) * t486 + (qJD(1) * t530 + t47 * t773 - t6 * t927 - t882) * t483) * t485 - (-t47 * t771 + t48 * t772) * t455 - (t33 * t772 + t47 * t635) * t378 - (t33 * t771 + t48 * t635) * t377) * m(6) + (-t64 * (-t231 * t455 - t368 * t378) - t65 * (t229 * t455 - t368 * t377) - t55 * (t229 * t378 + t231 * t377) + (qJD(2) * t503 + t25 * t184 - t26 * t188 - t64 * t91 + t65 * t93) * t482 + (t64 * (-qJD(2) * t188 + t208 * t483) + t65 * (qJD(2) * t184 - t208 * t486) - t16 * t583 + t55 * (-t184 * t726 - t188 * t728 - t483 * t93 + t486 * t91) - (qJD(1) * t891 - t25 * t486 + t26 * t483) * t525) * t485) * m(5) + t947 * (t482 * t668 - t695 / 0.2e1) + t946 * (t482 * t671 + t485 * t672); 0.2e1 * (-t377 * t48 - t378 * t47 + t532) * t877 * t482 + 0.2e1 * (t496 + (t33 * (t377 * t483 + t378 * t486) + (-t47 * t483 + t48 * t486) * t455) * t877) * t485;];
tauc = t1(:);
