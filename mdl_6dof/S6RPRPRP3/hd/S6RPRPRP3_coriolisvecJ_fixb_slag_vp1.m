% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPRP3
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
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPRP3_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRP3_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:07:49
% EndTime: 2019-03-09 03:09:17
% DurationCPUTime: 80.22s
% Computational Cost: add. (45562->1228), mult. (47151->1576), div. (0->0), fcn. (44323->10), ass. (0->588)
t1004 = Icges(6,1) + Icges(7,1);
t1038 = Icges(7,4) + Icges(6,5);
t1037 = Icges(6,6) - Icges(7,6);
t517 = pkin(10) + qJ(5);
t510 = sin(t517);
t1060 = (Icges(6,4) - Icges(7,5)) * t510;
t1003 = Icges(6,2) + Icges(7,3);
t1059 = Icges(7,2) + Icges(6,3);
t512 = cos(t517);
t1058 = -t1037 * t510 + t1038 * t512;
t1057 = t1004 * t512 - t1060;
t522 = sin(qJ(3));
t524 = cos(qJ(3));
t862 = Icges(6,4) * t512;
t636 = -Icges(6,2) * t510 + t862;
t1056 = -t1037 * t524 + t522 * t636;
t518 = qJ(1) + pkin(9);
t513 = cos(t518);
t827 = t513 * t510;
t511 = sin(t518);
t830 = t511 * t524;
t361 = t512 * t830 - t827;
t341 = Icges(7,5) * t361;
t829 = t512 * t513;
t360 = t510 * t830 + t829;
t831 = t511 * t522;
t180 = -Icges(7,6) * t831 - Icges(7,3) * t360 - t341;
t344 = Icges(6,4) * t361;
t188 = -Icges(6,2) * t360 + Icges(6,6) * t831 + t344;
t1036 = t180 + t188;
t824 = t513 * t524;
t758 = t512 * t824;
t363 = t510 * t511 + t758;
t342 = Icges(7,5) * t363;
t362 = t510 * t824 - t511 * t512;
t825 = t513 * t522;
t181 = Icges(7,6) * t825 + Icges(7,3) * t362 + t342;
t864 = Icges(6,4) * t363;
t190 = -Icges(6,2) * t362 + Icges(6,6) * t825 + t864;
t1035 = t181 - t190;
t182 = Icges(6,5) * t361 - Icges(6,6) * t360 + Icges(6,3) * t831;
t185 = Icges(7,4) * t361 + Icges(7,2) * t831 + Icges(7,6) * t360;
t1034 = t182 + t185;
t184 = Icges(6,5) * t363 - Icges(6,6) * t362 + Icges(6,3) * t825;
t187 = Icges(7,4) * t363 + Icges(7,2) * t825 + Icges(7,6) * t362;
t1001 = t184 + t187;
t340 = Icges(7,5) * t360;
t191 = Icges(7,1) * t361 + Icges(7,4) * t831 + t340;
t343 = Icges(6,4) * t360;
t195 = -Icges(6,1) * t361 - Icges(6,5) * t831 + t343;
t1033 = t191 - t195;
t857 = Icges(7,5) * t362;
t193 = Icges(7,1) * t363 + Icges(7,4) * t825 + t857;
t345 = Icges(6,4) * t362;
t196 = Icges(6,1) * t363 + Icges(6,5) * t825 - t345;
t1032 = t193 + t196;
t828 = t512 * t522;
t487 = Icges(7,5) * t828;
t833 = t510 * t522;
t1031 = Icges(7,3) * t833 - t1056 + t487;
t855 = Icges(7,5) * t512;
t632 = Icges(7,3) * t510 + t855;
t943 = (t632 - t636) * t524 - t1037 * t522;
t1000 = t1058 * t522 - t1059 * t524;
t1055 = t1058 * t524 + t1059 * t522;
t942 = t1038 * t522 + t1057 * t524;
t999 = -t1038 * t524 + t1057 * t522;
t1054 = (t1003 * t512 + t1060) * t522;
t1053 = (-t1037 * t512 - t1038 * t510) * t522;
t520 = cos(pkin(10));
t519 = sin(pkin(10));
t823 = t519 * t524;
t387 = t511 * t823 + t513 * t520;
t822 = t520 * t524;
t826 = t513 * t519;
t388 = t511 * t822 - t826;
t219 = Icges(5,5) * t388 - Icges(5,6) * t387 + Icges(5,3) * t831;
t848 = Icges(4,3) * t513;
t325 = Icges(4,5) * t830 - Icges(4,6) * t831 - t848;
t488 = Icges(4,4) * t831;
t860 = Icges(4,5) * t513;
t329 = Icges(4,1) * t830 - t488 - t860;
t852 = Icges(4,6) * t513;
t327 = Icges(4,4) * t830 - Icges(4,2) * t831 - t852;
t838 = t327 * t522;
t622 = -t329 * t524 + t838;
t222 = Icges(5,4) * t388 - Icges(5,2) * t387 + Icges(5,6) * t831;
t225 = Icges(5,1) * t388 - Icges(5,4) * t387 + Icges(5,5) * t831;
t625 = -t222 * t387 + t225 * t388;
t1052 = t219 * t831 - t325 * t513 - t511 * t622 + t625;
t515 = Icges(4,4) * t524;
t638 = -Icges(4,2) * t522 + t515;
t328 = Icges(4,6) * t511 + t513 * t638;
t865 = Icges(4,4) * t522;
t470 = Icges(4,1) * t524 - t865;
t330 = Icges(4,5) * t511 + t470 * t513;
t309 = t330 * t830;
t466 = Icges(4,5) * t524 - Icges(4,6) * t522;
t326 = Icges(4,3) * t511 + t466 * t513;
t692 = t326 * t513 - t309;
t125 = -t328 * t831 - t692;
t389 = t511 * t520 - t513 * t823;
t832 = t511 * t519;
t390 = t513 * t822 + t832;
t221 = Icges(5,5) * t390 + Icges(5,6) * t389 + Icges(5,3) * t825;
t224 = Icges(5,4) * t390 + Icges(5,2) * t389 + Icges(5,6) * t825;
t227 = Icges(5,1) * t390 + Icges(5,4) * t389 + Icges(5,5) * t825;
t71 = t221 * t831 - t387 * t224 + t388 * t227;
t1051 = t125 + t71;
t634 = Icges(5,5) * t520 - Icges(5,6) * t519;
t391 = -Icges(5,3) * t524 + t522 * t634;
t637 = Icges(5,4) * t520 - Icges(5,2) * t519;
t393 = -Icges(5,6) * t524 + t522 * t637;
t641 = Icges(5,1) * t520 - Icges(5,4) * t519;
t395 = -Icges(5,5) * t524 + t522 * t641;
t467 = Icges(4,2) * t524 + t865;
t469 = Icges(4,1) * t522 + t515;
t616 = t467 * t522 - t469 * t524;
t465 = Icges(4,5) * t522 + Icges(4,6) * t524;
t834 = t465 * t513;
t1050 = -t387 * t393 + t388 * t395 + t391 * t831 - t511 * t616 - t834;
t835 = t465 * t511;
t1049 = t389 * t393 + t390 * t395 + t391 * t825 - t513 * t616 + t835;
t1009 = t1033 * t363 + t1034 * t825 - t1036 * t362;
t1008 = t1001 * t825 + t1032 * t363 + t1035 * t362;
t949 = t1000 * t825 + t1031 * t362 + t363 * t999;
t1048 = qJD(3) * t943 + qJD(5) * t1054;
t1047 = qJD(3) * t1055 + qJD(5) * t1053;
t411 = (-Icges(6,1) * t510 - t862) * t522;
t768 = qJD(5) * t522;
t1046 = (-Icges(7,1) * t510 + t855) * t768 + qJD(5) * t411 + t942 * qJD(3);
t1011 = t1033 * t361 + t1034 * t831 - t1036 * t360;
t1010 = t1001 * t831 + t1032 * t361 + t1035 * t360;
t772 = qJD(3) * t522;
t736 = t513 * t772;
t769 = qJD(5) * t511;
t169 = qJD(1) * t360 - qJD(5) * t758 + (t736 - t769) * t510;
t767 = qJD(5) * t524;
t498 = qJD(1) - t767;
t613 = t498 * t510;
t775 = qJD(1) * t524;
t689 = -qJD(5) + t775;
t170 = t513 * t613 + (-t511 * t689 - t736) * t512;
t771 = qJD(3) * t524;
t735 = t513 * t771;
t776 = qJD(1) * t522;
t744 = t511 * t776;
t399 = t735 - t744;
t86 = Icges(7,5) * t170 + Icges(7,6) * t399 - Icges(7,3) * t169;
t92 = Icges(6,4) * t170 + Icges(6,2) * t169 + Icges(6,6) * t399;
t1045 = t86 - t92;
t741 = t511 * t772;
t777 = qJD(1) * t513;
t778 = qJD(1) * t511;
t171 = -t510 * t741 - qJD(5) * t827 - t512 * t778 + (t510 * t777 + t512 * t769) * t524;
t172 = t689 * t829 + (-t512 * t772 + t613) * t511;
t740 = t511 * t771;
t398 = t513 * t776 + t740;
t87 = Icges(7,5) * t172 + Icges(7,6) * t398 + Icges(7,3) * t171;
t93 = Icges(6,4) * t172 - Icges(6,2) * t171 + Icges(6,6) * t398;
t1044 = t87 - t93;
t88 = Icges(6,5) * t170 + Icges(6,6) * t169 + Icges(6,3) * t399;
t90 = Icges(7,4) * t170 + Icges(7,2) * t399 - Icges(7,6) * t169;
t1043 = t88 + t90;
t89 = Icges(6,5) * t172 - Icges(6,6) * t171 + Icges(6,3) * t398;
t91 = Icges(7,4) * t172 + Icges(7,2) * t398 + Icges(7,6) * t171;
t1042 = t89 + t91;
t94 = Icges(7,1) * t170 + Icges(7,4) * t399 - Icges(7,5) * t169;
t96 = Icges(6,1) * t170 + Icges(6,4) * t169 + Icges(6,5) * t399;
t1041 = t94 + t96;
t95 = Icges(7,1) * t172 + Icges(7,4) * t398 + Icges(7,5) * t171;
t97 = Icges(6,1) * t172 - Icges(6,4) * t171 + Icges(6,5) * t398;
t1040 = t95 + t97;
t973 = t389 * t222 + t390 * t225;
t72 = t219 * t825 + t973;
t73 = t221 * t825 + t389 * t224 + t390 * t227;
t589 = (t511 * t73 - t513 * t72) * qJD(3);
t1039 = qJD(1) * t1049 + t589;
t950 = t1000 * t831 + t1031 * t360 + t361 * t999;
t1030 = (t1051 * t511 - t1052 * t513) * qJD(3);
t1029 = t1031 * t510 + t512 * t999;
t1028 = t1050 * qJD(1);
t992 = t188 * t510 + t195 * t512;
t993 = t180 * t510 - t191 * t512;
t1027 = t992 + t993;
t991 = rSges(7,1) + pkin(5);
t525 = cos(qJ(1));
t516 = t525 * pkin(1);
t774 = qJD(3) * t511;
t427 = t513 * t768 + t774;
t773 = qJD(3) * t513;
t428 = -t511 * t768 + t773;
t960 = t1008 * t427 - t1009 * t428 + t949 * t498;
t1026 = t1028 + t1030;
t808 = -t511 * t325 - t329 * t824;
t126 = -t327 * t825 - t808;
t807 = t511 * t326 + t330 * t824;
t127 = -t328 * t825 + t807;
t583 = (-t126 * t513 + t127 * t511) * qJD(3);
t1025 = t583 + t1039;
t734 = t519 * t772;
t301 = qJD(1) * t389 + t511 * t734;
t733 = t520 * t772;
t302 = qJD(1) * t390 - t511 * t733;
t131 = Icges(5,5) * t302 + Icges(5,6) * t301 + Icges(5,3) * t398;
t133 = Icges(5,4) * t302 + Icges(5,2) * t301 + Icges(5,6) * t398;
t135 = Icges(5,1) * t302 + Icges(5,4) * t301 + Icges(5,5) * t398;
t586 = qJD(3) * t467;
t587 = qJD(3) * t469;
t972 = t222 * t519 - t225 * t520;
t1024 = (-qJD(1) * t328 + t511 * t586 + t131) * t524 + (-qJD(1) * t330 + t133 * t519 - t135 * t520 + t511 * t587) * t522 + (-t219 * t522 + t524 * t972 + t622) * qJD(3);
t299 = qJD(1) * t387 + t513 * t734;
t300 = -qJD(1) * t388 - t513 * t733;
t130 = Icges(5,5) * t300 + Icges(5,6) * t299 + Icges(5,3) * t399;
t132 = Icges(5,4) * t300 + Icges(5,2) * t299 + Icges(5,6) * t399;
t134 = Icges(5,1) * t300 + Icges(5,4) * t299 + Icges(5,5) * t399;
t837 = t328 * t522;
t621 = -t330 * t524 + t837;
t623 = -t224 * t519 + t227 * t520;
t1023 = (-t130 - t513 * t586 + (-t511 * t638 + t852) * qJD(1)) * t524 + (-t132 * t519 + t134 * t520 - t513 * t587 + (-t470 * t511 + t860) * qJD(1)) * t522 + (t221 * t522 + t524 * t623 - t621) * qJD(3);
t1014 = t1000 * t399 - t1031 * t169 + t1046 * t363 + t1047 * t825 + t1048 * t362 + t170 * t999;
t1013 = t1000 * t398 + t1031 * t171 + t1046 * t361 + t1047 * t831 + t1048 * t360 + t172 * t999;
t986 = t185 * t524;
t74 = -t522 * t993 - t986;
t989 = t182 * t524;
t76 = -t522 * t992 - t989;
t1007 = t74 + t76;
t990 = rSges(7,3) + qJ(6);
t392 = Icges(5,3) * t522 + t524 * t634;
t357 = t392 * qJD(3);
t394 = Icges(5,6) * t522 + t524 * t637;
t358 = t394 * qJD(3);
t396 = Icges(5,5) * t522 + t524 * t641;
t359 = t396 * qJD(3);
t443 = t638 * qJD(3);
t444 = t470 * qJD(3);
t532 = qJD(1) * t465 - t443 * t522 + t444 * t524 + (-t467 * t524 - t469 * t522) * qJD(3);
t916 = qJD(1) * t616 + t466 * qJD(3);
t1022 = t299 * t393 + t300 * t395 + t357 * t825 + t358 * t389 + t359 * t390 + t391 * t399 + t511 * t916 + t532 * t513;
t1021 = t301 * t393 + t302 * t395 + t357 * t831 - t358 * t387 + t359 * t388 + t391 * t398 + t532 * t511 - t513 * t916;
t1020 = t126 + t72;
t1019 = (-t219 + t327) * t524 + (t329 - t972) * t522;
t1018 = (-t221 + t328) * t524 + (t330 + t623) * t522;
t332 = qJD(6) * t362;
t982 = -t990 * t360 - t361 * t991;
t819 = rSges(7,2) * t831 - t982;
t1017 = -t498 * t819 + t332;
t331 = qJD(6) * t360;
t656 = pkin(5) * t512 + qJ(6) * t510;
t657 = rSges(7,1) * t512 + rSges(7,3) * t510;
t795 = -rSges(7,2) * t524 + (t656 + t657) * t522;
t817 = rSges(7,2) * t825 + t362 * t990 + t363 * t991;
t1016 = -t427 * t795 + t498 * t817 + t331;
t627 = -t190 * t510 + t196 * t512;
t629 = t181 * t510 + t193 * t512;
t1015 = t427 * (-t1000 * t513 - t627 - t629) - t428 * (-t1000 * t511 + t1027) + t498 * (-t1029 + t1055);
t968 = t1033 * t170 + t1034 * t399 + t1036 * t169 + t1040 * t363 + t1042 * t825 + t1044 * t362;
t966 = t1001 * t399 + t1032 * t170 - t1035 * t169 + t1041 * t363 + t1043 * t825 + t1045 * t362;
t965 = t1033 * t172 + t1034 * t398 - t1036 * t171 + t1040 * t361 + t1042 * t831 + t1044 * t360;
t964 = t1001 * t398 + t1032 * t172 + t1035 * t171 + t1041 * t361 + t1043 * t831 + t1045 * t360;
t1012 = (qJD(3) * t1029 - t1047) * t524 + (t1046 * t512 + t1048 * t510 + (t1031 * t512 - t510 * t999) * qJD(5) + t1000 * qJD(3)) * t522;
t75 = -t187 * t524 + t522 * t629;
t77 = -t184 * t524 + t522 * t627;
t1006 = t75 + t77;
t585 = qJD(3) * t465;
t1005 = t513 * (-t511 * t585 + (t326 + t622) * qJD(1));
t948 = -t1000 * t524 + t1029 * t522;
t998 = -t522 * t632 + t1056;
t997 = -t1053 * t498 + (-t1037 * t361 - t1038 * t360) * t428 + (t1037 * t363 + t1038 * t362) * t427;
t526 = qJD(1) ^ 2;
t961 = t1010 * t427 - t1011 * t428 + t950 * t498;
t502 = t511 * rSges(4,3);
t351 = rSges(4,1) * t824 - rSges(4,2) * t825 + t502;
t441 = t513 * pkin(2) + t511 * pkin(7);
t939 = t516 + t441;
t981 = t351 + t939;
t494 = pkin(3) * t824;
t423 = qJ(4) * t825 + t494;
t483 = pkin(4) * t832;
t508 = pkin(4) * t520 + pkin(3);
t521 = -pkin(8) - qJ(4);
t821 = t521 * t522;
t612 = t508 * t824 - t513 * t821 + t483;
t233 = t612 - t423;
t673 = t423 + t939;
t614 = t233 + t673;
t895 = pkin(3) - t508;
t724 = t895 * t522;
t820 = qJ(4) + t521;
t567 = -t524 * t820 + t724;
t514 = qJD(4) * t522;
t473 = t511 * t514;
t477 = pkin(3) * t522 - qJ(4) * t524;
t785 = -t477 * t774 + t473;
t748 = t567 * t774 + t785;
t547 = qJD(1) * t614 + t748;
t44 = t1016 + t547;
t701 = t44 * t795;
t766 = qJD(6) * t510;
t471 = t522 * t766;
t761 = pkin(4) * t826;
t675 = -t508 * t830 + t761;
t898 = pkin(3) * t524;
t955 = t522 * t820;
t232 = (t898 + t955) * t511 + t675;
t845 = qJ(4) * t522;
t480 = t845 + t898;
t416 = t480 * t511;
t770 = qJD(4) * t524;
t609 = t416 * t774 + t423 * t773 + qJD(2) - t770;
t570 = -t232 * t774 + t233 * t773 + t609;
t38 = t427 * t819 + t428 * t817 + t471 + t570;
t703 = t38 * t819;
t979 = t701 - t703;
t725 = t895 * t524;
t978 = t725 + t845;
t661 = rSges(6,1) * t361 - rSges(6,2) * t360;
t199 = rSges(6,3) * t831 + t661;
t660 = rSges(6,1) * t512 - rSges(6,2) * t510;
t378 = -rSges(6,3) * t524 + t522 * t660;
t975 = -t199 * t498 - t378 * t428;
t617 = -t393 * t519 + t395 * t520;
t974 = qJD(3) * (-t511 * t623 - t513 * t972) + (t392 - t617) * qJD(1);
t971 = 0.2e1 * qJD(3);
t762 = qJD(3) * qJD(5);
t719 = t524 * t762;
t319 = qJD(1) * t427 + t511 * t719;
t320 = qJD(1) * t428 + t513 * t719;
t720 = t522 * t762;
t970 = t1008 * t320 + t1009 * t319 + t1014 * t498 + t427 * t966 - t428 * t968 + t720 * t949;
t969 = t1010 * t320 + t1011 * t319 + t1013 * t498 + t427 * t964 - t428 * t965 + t720 * t950;
t456 = qJ(4) * t735;
t743 = t521 * t776;
t789 = qJD(1) * t761 + t511 * t743;
t122 = -t456 + (-t521 * t524 + t724) * t773 + t978 * t778 + t789;
t462 = pkin(3) * t741;
t732 = t521 * t771;
t749 = -t508 * t741 - t511 * t732 - t513 * t743;
t123 = -qJ(4) * t740 + t462 + (-t513 * t978 + t483) * qJD(1) + t749;
t475 = t513 * t514;
t578 = -t511 * t775 - t736;
t247 = pkin(3) * t578 - qJ(4) * t744 + t456 + t475;
t784 = t473 - t462;
t248 = qJ(4) * t398 + qJD(1) * t494 + t784;
t764 = qJD(1) * qJD(3);
t722 = t513 * t764;
t763 = qJD(3) * qJD(4);
t687 = t247 * t773 + t248 * t774 + t416 * t722 + t522 * t763;
t581 = t122 * t773 + t123 * t774 - t232 * t722 + t687;
t815 = -t233 - t423;
t669 = t815 * t778;
t730 = t524 * t766;
t731 = t512 * t768;
t958 = -t990 * t171 - t172 * t991 - t331;
t869 = rSges(7,2) * t398 - t958;
t957 = rSges(7,2) * t735 - t990 * t169 + t170 * t991 + t332;
t894 = -rSges(7,2) * t744 + t957;
t5 = qJD(6) * t731 + t894 * t428 + t869 * t427 + t819 * t320 - t817 * t319 + (t669 + t730) * qJD(3) + t581;
t967 = t5 * t817;
t17 = (-qJD(3) * t993 - t91) * t524 + (qJD(3) * t185 + t510 * t87 + t512 * t95 + (-t180 * t512 - t191 * t510) * qJD(5)) * t522;
t19 = (-qJD(3) * t992 - t89) * t524 + (qJD(3) * t182 - t510 * t93 + t512 * t97 + (-t188 * t512 + t195 * t510) * qJD(5)) * t522;
t963 = t17 + t19;
t18 = (qJD(3) * t629 - t90) * t524 + (qJD(3) * t187 + t510 * t86 + t512 * t94 + (t181 * t512 - t193 * t510) * qJD(5)) * t522;
t20 = (qJD(3) * t627 - t88) * t524 + (qJD(3) * t184 - t510 * t92 + t512 * t96 + (-t190 * t512 - t196 * t510) * qJD(5)) * t522;
t962 = t18 + t20;
t959 = t1006 * t427 - t1007 * t428 + t498 * t948;
t947 = t998 * t511;
t946 = t998 * t513;
t945 = t999 * t511;
t944 = t999 * t513;
t799 = t567 - t477;
t506 = t513 * pkin(7);
t440 = pkin(2) * t511 - t506;
t434 = qJD(1) * t440;
t941 = -qJD(1) * t416 - t434;
t354 = -t725 - t955;
t429 = qJD(3) * t480 - t770;
t805 = -t354 * qJD(3) - t429;
t891 = rSges(7,2) * t522;
t379 = t524 * t657 + t891;
t418 = (-rSges(7,1) * t510 + rSges(7,3) * t512) * t522;
t577 = t510 * t771 + t731;
t814 = t471 + t577 * qJ(6) + (-t510 * t768 + t512 * t771) * pkin(5) + qJD(3) * t379 + qJD(5) * t418;
t940 = t471 + t805 - t814;
t699 = t513 * rSges(3,1) - rSges(3,2) * t511;
t938 = t516 + t699;
t937 = t1015 * t522;
t936 = (-t999 + t1054) * t498 + (-t1003 * t361 + t1033 + t340 - t343) * t428 + (t1003 * t363 - t1032 + t345 - t857) * t427;
t935 = (-Icges(7,1) * t833 + t1031 + t411 + t487) * t498 + (t1004 * t360 + t1036 - t341 + t344) * t428 + (-t1004 * t362 + t1035 + t342 - t864) * t427;
t934 = t997 * t522;
t933 = t1000 * t498 + t1001 * t427 - t1034 * t428;
t279 = t567 * t511;
t280 = t567 * t513;
t414 = t477 * t511;
t421 = t477 * t513;
t751 = -t414 * t774 - t421 * t773 + t514;
t755 = t513 * t247 + t511 * t248 + t416 * t777;
t932 = t513 * t122 + t511 * t123 - t232 * t777 - t279 * t774 - t280 * t773 - t751 + t755;
t931 = t1006 * t513 + t1007 * t511;
t930 = t1006 * t511 - t1007 * t513;
t929 = t1008 * t513 + t1009 * t511;
t928 = t1008 * t511 - t1009 * t513;
t927 = t1010 * t513 + t1011 * t511;
t926 = t1010 * t511 - t1011 * t513;
t923 = t1012 * t498 + t720 * t948;
t663 = rSges(5,1) * t520 - rSges(5,2) * t519;
t400 = -rSges(5,3) * t524 + t522 * t663;
t918 = -t513 * t585 + (-t466 * t511 + t621 + t848) * qJD(1);
t802 = -Icges(4,2) * t830 + t329 - t488;
t804 = t469 * t511 + t327;
t914 = -t522 * t802 - t524 * t804;
t913 = m(5) / 0.2e1;
t912 = m(6) / 0.2e1;
t911 = m(7) / 0.2e1;
t910 = t319 / 0.2e1;
t909 = t320 / 0.2e1;
t908 = -t427 / 0.2e1;
t907 = t427 / 0.2e1;
t906 = -t428 / 0.2e1;
t905 = t428 / 0.2e1;
t904 = -t498 / 0.2e1;
t903 = t498 / 0.2e1;
t523 = sin(qJ(1));
t899 = pkin(1) * t523;
t893 = rSges(4,1) * t524;
t889 = rSges(5,3) * t522;
t887 = rSges(6,3) * t522;
t885 = t17 * t428;
t884 = t18 * t427;
t883 = t19 * t428;
t882 = t20 * t427;
t697 = t513 * t799;
t594 = qJD(3) * t697 + t475;
t708 = -t440 - t899;
t674 = -t416 + t708;
t538 = (t232 + t674) * qJD(1) + t594;
t60 = t538 + t975;
t876 = t513 * t60;
t874 = t74 * t319;
t873 = t75 * t320;
t872 = t76 * t319;
t871 = t77 * t320;
t870 = -rSges(5,3) - qJ(4);
t782 = rSges(4,2) * t831 + t513 * rSges(4,3);
t350 = rSges(4,1) * t830 - t782;
t478 = rSges(4,1) * t522 + rSges(4,2) * t524;
t737 = t478 * t773;
t174 = -t737 + (-t350 + t708) * qJD(1);
t842 = t174 * t511;
t841 = t174 * t513;
t742 = t478 * t774;
t175 = qJD(1) * t981 - t742;
t422 = t478 * t513;
t840 = t175 * t422;
t230 = t390 * rSges(5,1) + t389 * rSges(5,2) + rSges(5,3) * t825;
t816 = -t230 - t423;
t430 = t441 * qJD(1);
t813 = -t248 - t430;
t812 = -t360 * t991 + t361 * t990;
t811 = -t362 * t991 + t363 * t990;
t810 = t795 * t511;
t809 = t795 * t513;
t432 = t477 * t778;
t806 = -t567 * t778 + t432;
t803 = -t469 * t513 - t328;
t801 = -t467 * t513 + t330;
t798 = -t354 - t480;
t797 = t511 * t416 + t513 * t423;
t401 = t524 * t663 + t889;
t796 = -t401 * qJD(3) - t429;
t794 = -t656 * t524 - t379;
t793 = -qJD(1) * t421 + t511 * t770;
t792 = -t400 - t477;
t791 = -t401 - t480;
t790 = -(-pkin(5) * t510 + qJ(6) * t512) * t522 - t418;
t788 = rSges(4,2) * t744 + rSges(4,3) * t777;
t787 = -t467 + t470;
t786 = t469 + t638;
t497 = pkin(7) * t777;
t783 = t475 + t497;
t779 = qJD(1) * t466;
t760 = t526 * t899;
t759 = t526 * t516;
t757 = t170 * rSges(6,1) + t169 * rSges(6,2) + rSges(6,3) * t735;
t380 = t524 * t660 + t887;
t419 = (-rSges(6,1) * t510 - rSges(6,2) * t512) * t522;
t277 = qJD(3) * t380 + qJD(5) * t419;
t754 = -t277 + t805;
t753 = qJD(1) * t280 + t793;
t752 = t300 * rSges(5,1) + t299 * rSges(5,2) + rSges(5,3) * t735;
t203 = t363 * rSges(6,1) - t362 * rSges(6,2) + rSges(6,3) * t825;
t750 = -t378 + t799;
t747 = t475 + t941;
t745 = -pkin(4) * t519 - pkin(7);
t726 = -pkin(2) - t893;
t723 = t511 * t764;
t721 = t524 * t763;
t717 = t777 / 0.2e1;
t716 = -t774 / 0.2e1;
t713 = t773 / 0.2e1;
t712 = t772 / 0.2e1;
t711 = t771 / 0.2e1;
t706 = t506 - t899;
t705 = -t508 * t524 - pkin(2);
t704 = -pkin(2) + t821;
t43 = -t428 * t795 + t1017 + t538;
t702 = t43 * t795;
t105 = -t400 * t774 + (t230 + t673) * qJD(1) + t785;
t698 = t105 * t792;
t696 = t513 * t792;
t695 = t805 * t513;
t694 = t798 * t511;
t693 = t798 * t513;
t691 = -t325 + t837;
t690 = qJD(1) * t414 + t513 * t770;
t686 = -t511 * t232 + t513 * t233 + t797;
t684 = -t795 + t799;
t683 = -t473 - t749;
t676 = qJD(1) * (-pkin(2) * t778 + t497) - t760;
t670 = qJD(5) * t712;
t439 = rSges(3,1) * t511 + rSges(3,2) * t513;
t666 = -rSges(4,2) * t522 + t893;
t665 = rSges(5,1) * t302 + rSges(5,2) * t301;
t664 = rSges(5,1) * t388 - rSges(5,2) * t387;
t662 = rSges(6,1) * t172 - rSges(6,2) * t171;
t631 = -t175 * t511 - t841;
t626 = t199 * t513 - t203 * t511;
t620 = t350 * t511 + t351 * t513;
t615 = t477 * t723 + t513 * t721 - t759;
t611 = t705 - t891;
t610 = t705 - t887;
t415 = t478 * t511;
t592 = -qJD(1) * t279 + t690;
t588 = (-t219 * t513 + t221 * t511) * t524;
t582 = t522 * t870 - pkin(2) - t898;
t579 = t511 * t721 + t676 + (t247 + t475) * qJD(1);
t576 = t675 + t706;
t575 = t38 * t869 + t5 * t819;
t569 = -t43 * t819 + t44 * t817;
t568 = -t38 * t817 + t702;
t566 = -t522 * t801 + t524 * t803;
t228 = rSges(5,3) * t831 + t664;
t565 = qJD(1) * t122 + t579;
t564 = t612 + t939;
t558 = qJD(1) * t697 + t511 * t805;
t548 = (-t522 * t786 + t524 * t787) * qJD(1);
t268 = rSges(4,1) * t578 - rSges(4,2) * t735 + t788;
t269 = -qJD(3) * t415 + (t513 * t666 + t502) * qJD(1);
t542 = t268 * t513 + t269 * t511 + (t350 * t513 - t351 * t511) * qJD(1);
t541 = -t508 * t736 - t513 * t732 + t783 + t789;
t51 = t199 * t427 + t203 * t428 + t570;
t61 = t203 * t498 - t378 * t427 + t547;
t535 = t51 * t626 + (t511 * t60 - t513 * t61) * t378;
t531 = -t567 * t723 + (-t123 - t473 + t813) * qJD(1) + t615;
t530 = t568 * t511 - t513 * t979;
t527 = t391 * t775 + t522 * t974;
t448 = t666 * qJD(3);
t397 = (t511 ^ 2 + t513 ^ 2) * t772;
t313 = t400 * t513;
t312 = t400 * t511;
t308 = t395 * t513;
t307 = t395 * t511;
t306 = t393 * t513;
t305 = t393 * t511;
t296 = t378 * t513;
t294 = t378 * t511;
t256 = -rSges(6,1) * t362 - rSges(6,2) * t363;
t251 = -rSges(6,1) * t360 - rSges(6,2) * t361;
t217 = qJD(1) * t232;
t173 = qJD(3) * t620 + qJD(2);
t137 = rSges(5,3) * t398 + t665;
t136 = -rSges(5,3) * t744 + t752;
t117 = -t759 - t448 * t773 + (-t269 - t430 + t742) * qJD(1);
t116 = -t448 * t774 + (t268 - t737) * qJD(1) + t676;
t104 = t475 + qJD(3) * t696 + (-t228 + t674) * qJD(1);
t101 = rSges(6,3) * t398 + t662;
t99 = -rSges(6,3) * t744 + t757;
t79 = (t228 * t511 + t230 * t513) * qJD(3) + t609;
t78 = t542 * qJD(3);
t55 = t796 * t773 + (-t137 + (qJD(3) * t400 - t514) * t511 + t813) * qJD(1) + t615;
t54 = qJD(1) * t136 + (qJD(1) * t696 + t511 * t796) * qJD(3) + t579;
t35 = (t136 * t513 + t137 * t511 + (t228 * t513 + t511 * t816) * qJD(1)) * qJD(3) + t687;
t22 = -t101 * t498 - t277 * t428 + t319 * t378 + (-t199 * t768 + t695) * qJD(3) + t531;
t21 = -t277 * t427 - t320 * t378 + t498 * t99 + (t203 * t768 + t558) * qJD(3) + t565;
t8 = qJD(3) * t669 + t101 * t427 + t199 * t320 - t203 * t319 + t428 * t99 + t581;
t7 = -qJD(6) * t169 - t869 * t498 - t814 * t428 + t795 * t319 + (-t768 * t819 + t695) * qJD(3) + t531;
t6 = qJD(6) * t171 + t894 * t498 - t814 * t427 - t795 * t320 + (t768 * t817 + t558) * qJD(3) + t565;
t1 = [(t22 * (t576 - t661) + t60 * (-t662 + t683) + t21 * (t564 + t203) + t61 * (t541 + t757) + (t22 * (t704 - t887) - t60 * rSges(6,3) * t771) * t511 + ((-t523 * t61 - t525 * t60) * pkin(1) + t610 * t876 + (t60 * t745 + t61 * t610) * t511) * qJD(1) - (-qJD(1) * t899 + t217 + t594 - t60 + t941 + t975) * t61) * m(6) + m(3) * ((-t439 * t526 - t760) * t938 + (-t759 + (-0.2e1 * t699 - t516 + t938) * t526) * (-t439 - t899)) + (-(-t104 + (-t228 - t899) * qJD(1) + t747) * t105 - t698 * t773 + t55 * (-t664 + t706) + t104 * (-t665 - t784) + t54 * (t939 - t816) + t105 * (-pkin(3) * t736 + t456 + t752 + t783) + (t104 * t771 * t870 + t55 * t582) * t511 + ((-t104 * t525 - t105 * t523) * pkin(1) + t104 * t582 * t513 + (-t104 * pkin(7) + t105 * (-pkin(2) - t480 - t889)) * t511) * qJD(1)) * m(5) + ((t1018 + t1049) * t513 + (t1019 + t1050) * t511) * t764 / 0.2e1 + t1014 * t907 + (t1022 + t1023) * t774 / 0.2e1 - (t1021 - t1024 + t1025) * t773 / 0.2e1 + (((t513 * t691 + t127 + t625 - t807) * t513 + (t511 * t691 + t1020 + t692 - t71 - t973) * t511) * qJD(3) + t1026 - t1028) * t716 + t872 / 0.2e1 + t873 / 0.2e1 + (((t125 - t309 + (t326 + t838) * t513 + t808) * t513 + t807 * t511) * qJD(3) + t1039) * t713 + (t117 * (t511 * t726 + t706 + t782) + t116 * t981 + t175 * (t497 + t788) + (t478 * t842 - t840) * qJD(3) + ((-t174 * t525 - t175 * t523) * pkin(1) + (-pkin(2) - t666) * t841 + (t174 * (-rSges(4,3) - pkin(7)) + t175 * t726) * t511) * qJD(1) - (-t737 - t174 - t434 + (-t350 - t899) * qJD(1)) * t175) * m(4) + t874 / 0.2e1 + t923 + (t1013 + t960) * t906 + t960 * t905 + (t6 * (t564 + t817) + t701 * t428 + (t576 + t982 + (t704 - t891) * t511) * t7 + (t611 * t778 - t773 * t799 - t1017 - t217 + t541 - t747 + t957) * t44 + (-rSges(7,2) * t740 + t683 + t748 + t958 + (t511 * t745 + t513 * t611 - t516 + t614) * qJD(1) + t1016) * t43) * m(7) + (t1027 * t522 + t1007 + t986 + t989) * t427 * t904 + ((t443 - t357) * t524 + (-t358 * t519 + t359 * t520 + t444) * t522 + (t391 * t522 + t524 * t617 - t616) * qJD(3)) * qJD(1) + t871 / 0.2e1 - t885 / 0.2e1 - t883 / 0.2e1 + t884 / 0.2e1 + t882 / 0.2e1 + t949 * t909 + t950 * t910; m(4) * t78 + m(5) * t35 + m(6) * t8 + m(7) * t5; t926 * t910 + t928 * t909 + ((t1009 * t830 + t949 * t522) * qJD(5) + ((qJD(5) * t1008 + t933) * t524 + t937) * t513 + (t362 * t943 + t363 * t942) * t498 + (-t362 * t947 + t363 * t945) * t428 + (t362 * t946 - t363 * t944) * t427) * t908 + (qJD(1) * t929 + t511 * t966 - t513 * t968) * t907 + (qJD(1) * t927 + t511 * t964 - t513 * t965) * t906 + ((t1010 * t824 + t950 * t522) * qJD(5) + ((qJD(5) * t1011 + t933) * t524 + t937) * t511 + (t360 * t943 + t361 * t942) * t498 + (-t360 * t947 + t361 * t945) * t428 + (t360 * t946 - t361 * t944) * t427) * t905 + ((t931 * qJD(5) - t1015) * t524 + ((t510 * t943 + t512 * t942 + t1000) * t498 + (-t510 * t947 + t512 * t945 - t1034) * t428 + (t510 * t946 - t512 * t944 + t1001) * t427 + t948 * qJD(5)) * t522) * t904 + (qJD(1) * t931 + t511 * t962 - t513 * t963) * t903 - (-t974 * t524 + ((-t394 * t519 + t396 * t520 + t391) * qJD(1) + ((t306 * t519 - t308 * t520 + t221) * t511 - (t305 * t519 - t307 * t520 + t219) * t513) * qJD(3)) * t522 + (t522 * t787 + t524 * t786) * qJD(1) + ((t511 * t801 - t513 * t802) * t524 + (t511 * t803 + t513 * t804) * t522) * qJD(3)) * qJD(1) / 0.2e1 + (t1024 * t513 + t1023 * t511 + (t1018 * t513 + t1019 * t511) * qJD(1)) * qJD(1) / 0.2e1 + ((-t306 * t389 - t308 * t390) * t774 + (t389 * t394 + t390 * t396) * qJD(1) + ((t305 * t389 + t307 * t390 + t588) * qJD(3) + t527) * t513 + (-t774 * t834 + t779) * t511 + (t548 + (-t914 * t513 + (t835 + t566) * t511) * qJD(3)) * t513) * t716 + (-(t305 * t387 - t307 * t388) * t773 + (-t387 * t394 + t388 * t396) * qJD(1) + ((t306 * t387 - t308 * t388 + t588) * qJD(3) + t527) * t511 + (-t773 * t835 - t779) * t513 + (t548 + (t566 * t511 + (t834 - t914) * t513) * qJD(3)) * t511) * t713 - t959 * t768 / 0.2e1 + t930 * t670 + (t5 * t686 + (qJD(1) * t703 + t7 * t684 + t967) * t513 + (qJD(1) * t702 + t6 * t684 + t575) * t511 - (t522 * t569 + t524 * t530) * qJD(5) + (-t694 * qJD(3) - t794 * t427 + t809 * t498 + t511 * t940 + t684 * t777 - t753) * t44 + (t894 * t513 + (t815 - t817) * t778 - t730 + t809 * t428 + t810 * t427 + t932) * t38 + (-t693 * qJD(3) - t794 * t428 - t810 * t498 + t513 * t940 - t592 + t806) * t43) * m(7) + (t60 * t806 + t8 * t686 + (t8 * t203 + t60 * t754 + (qJD(1) * t61 + t22) * t750) * t513 + (qJD(1) * t378 * t60 + t8 * t199 + t21 * t750 + t61 * t754) * t511 - t60 * (t294 * t498 - t380 * t428 + t592) - t61 * (-t296 * t498 - t380 * t427 + t753) - (t60 * t693 + t61 * t694) * qJD(3) - ((-t199 * t60 + t203 * t61) * t522 + t535 * t524) * qJD(5) + ((qJD(1) * t199 + t99) * t513 + (t101 + (-t203 + t815) * qJD(1)) * t511 + t294 * t427 + t296 * t428 + t932) * t51) * m(6) + (t104 * t432 + t35 * t797 + t79 * t755 + (t55 * t792 + t104 * t796 + t35 * t230 + t79 * t136 + (t79 * t228 + t698) * qJD(1)) * t513 + (t54 * t792 + t105 * t796 + t35 * t228 + t79 * t137 + (t104 * t400 + t79 * t816) * qJD(1)) * t511 - t104 * (qJD(1) * t312 + t690) - t105 * (-qJD(1) * t313 + t793) - t79 * t751 - ((t104 * t791 - t79 * t313) * t513 + (t105 * t791 - t79 * t312) * t511) * qJD(3)) * m(5) + (-(t174 * t415 - t840) * qJD(1) - (t173 * (-t415 * t511 - t422 * t513) + t631 * t666) * qJD(3) + t78 * t620 + t173 * t542 + t631 * t448 + (-t116 * t511 - t117 * t513 + (-t175 * t513 + t842) * qJD(1)) * t478) * m(4) + (t1022 * qJD(1) + t970 + ((-t131 * t825 - t133 * t389 - t135 * t390 - t219 * t399 - t222 * t299 - t225 * t300) * t513 + (t130 * t825 + t132 * t389 + t134 * t390 + t221 * t399 + t224 * t299 + t227 * t300 + t918 * t511 - t1005) * t511 + ((t127 + t73) * t513 + t1020 * t511) * qJD(1)) * t971) * t511 / 0.2e1 - (t1021 * qJD(1) + t969 + ((-t131 * t831 + t133 * t387 - t135 * t388 - t219 * t398 - t222 * t301 - t225 * t302 + t1005) * t513 + (t130 * t831 - t132 * t387 + t134 * t388 + t221 * t398 + t224 * t301 + t227 * t302 - t513 * t918) * t511 + (t1051 * t513 + t1052 * t511) * qJD(1)) * t971) * t513 / 0.2e1 + (t961 + t1026 + t1030) * t778 / 0.2e1 + (t589 + t583 + t960 + t1025) * t717 - (t511 * t961 + t513 * t960) * t767 / 0.2e1; -m(5) * (t104 * t399 + t105 * t398 + t397 * t79) - m(6) * (t397 * t51 + t398 * t61 + t399 * t60) - m(7) * (t38 * t397 + t398 * t44 + t399 * t43) + 0.2e1 * ((t104 * t773 + t105 * t774 - t35) * t913 + (t60 * t773 + t61 * t774 - t8) * t912 + (t43 * t773 + t44 * t774 - t5) * t911) * t524 + 0.2e1 * ((qJD(3) * t79 - t104 * t778 + t105 * t777 + t511 * t54 + t513 * t55) * t913 + (qJD(3) * t51 + t21 * t511 + t22 * t513 - t60 * t778 + t61 * t777) * t912 + (qJD(3) * t38 - t43 * t778 + t44 * t777 + t511 * t6 + t513 * t7) * t911) * t522; (t522 * t927 - t524 * t950) * t910 + (t522 * t929 - t524 * t949) * t909 + (t362 * t936 + t363 * t935 - t513 * t934) * t908 + ((qJD(3) * t929 - t1014) * t524 + (-qJD(1) * t928 + qJD(3) * t949 + t511 * t968 + t513 * t966) * t522) * t907 + ((qJD(3) * t927 - t1013) * t524 + (-qJD(1) * t926 + qJD(3) * t950 + t511 * t965 + t513 * t964) * t522) * t906 + (t360 * t936 + t361 * t935 - t511 * t934) * t905 + (t997 * t524 + (t510 * t936 + t512 * t935) * t522) * t904 + ((qJD(3) * t931 - t1012) * t524 + (-qJD(1) * t930 + qJD(3) * t948 + t511 * t963 + t513 * t962) * t522) * t903 - (t873 + t874 + t884 - t885 + t871 + t872 + t882 - t883 + t923) * t524 / 0.2e1 + t969 * t831 / 0.2e1 + t970 * t825 / 0.2e1 + t959 * t712 + (t522 * t931 - t524 * t948) * t670 + ((qJD(3) * t530 + t43 * t869 - t44 * t894 - t6 * t817 + t7 * t819) * t524 + (t569 * qJD(3) + (qJD(1) * t568 - t44 * t814 - t6 * t795 + t575) * t513 + (qJD(1) * t979 - t38 * t894 + t43 * t814 + t7 * t795 - t967) * t511) * t522 - (t361 * t44 + t363 * t43 + t38 * t828) * qJD(6) - (-t43 * t812 + t44 * t811) * t498 - (t38 * t811 + t43 * t790) * t428 - (t38 * t812 + t44 * t790) * t427) * m(7) + (-t60 * (-t251 * t498 - t419 * t428) - t61 * (t256 * t498 - t419 * t427) - t51 * (t251 * t427 + t256 * t428) + (qJD(3) * t535 + t60 * t101 + t22 * t199 - t21 * t203 - t61 * t99) * t524 + (t60 * (-qJD(3) * t199 + t277 * t511) + t61 * (qJD(3) * t203 - t277 * t513) + t8 * t626 + t51 * (t101 * t513 - t199 * t778 - t203 * t777 - t511 * t99) + (-t21 * t513 + t22 * t511 + (t511 * t61 + t876) * qJD(1)) * t378) * t522) * m(6) + t960 * (t513 * t711 - t744 / 0.2e1) + t961 * (t511 * t711 + t522 * t717); (t360 * t6 + t362 * t7 + t5 * t833 + (-t362 * t498 + t427 * t833 + t171) * t44 + (t360 * t498 + t428 * t833 - t169) * t43 + (-t360 * t427 - t362 * t428 + t577) * t38) * m(7);];
tauc  = t1(:);