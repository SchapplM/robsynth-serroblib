% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRRR4_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR4_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:04:48
% EndTime: 2019-03-09 07:05:46
% DurationCPUTime: 45.07s
% Computational Cost: add. (61854->1301), mult. (48102->1697), div. (0->0), fcn. (42697->12), ass. (0->687)
t579 = pkin(11) + qJ(3);
t550 = qJ(4) + t579;
t540 = qJ(5) + t550;
t523 = cos(t540);
t586 = sin(qJ(6));
t589 = cos(qJ(1));
t876 = t589 * t586;
t587 = sin(qJ(1));
t588 = cos(qJ(6));
t878 = t587 * t588;
t419 = t523 * t878 - t876;
t408 = Icges(7,4) * t419;
t877 = t588 * t589;
t879 = t586 * t587;
t418 = t523 * t879 + t877;
t522 = sin(t540);
t895 = t522 * t587;
t236 = -Icges(7,2) * t418 + Icges(7,6) * t895 + t408;
t407 = Icges(7,4) * t418;
t240 = -Icges(7,1) * t419 - Icges(7,5) * t895 + t407;
t1049 = t236 * t586 + t240 * t588;
t233 = Icges(7,5) * t419 - Icges(7,6) * t418 + Icges(7,3) * t895;
t105 = -t1049 * t522 - t233 * t523;
t420 = -t523 * t876 + t878;
t421 = t523 * t877 + t879;
t894 = t522 * t589;
t235 = Icges(7,5) * t421 + Icges(7,6) * t420 + Icges(7,3) * t894;
t944 = Icges(7,4) * t421;
t238 = Icges(7,2) * t420 + Icges(7,6) * t894 + t944;
t409 = Icges(7,4) * t420;
t241 = Icges(7,1) * t421 + Icges(7,5) * t894 + t409;
t100 = t235 * t895 - t418 * t238 + t419 * t241;
t698 = Icges(7,5) * t588 - Icges(7,6) * t586;
t306 = -Icges(7,3) * t523 + t522 * t698;
t942 = Icges(7,4) * t588;
t700 = -Icges(7,2) * t586 + t942;
t308 = -Icges(7,6) * t523 + t522 * t700;
t943 = Icges(7,4) * t586;
t705 = Icges(7,1) * t588 - t943;
t310 = -Icges(7,5) * t523 + t522 * t705;
t125 = t306 * t895 - t308 * t418 + t310 * t419;
t553 = qJD(3) * t587;
t503 = qJD(4) * t587 + t553;
t465 = qJD(5) * t587 + t503;
t813 = qJD(6) * t589;
t353 = t522 * t813 + t465;
t580 = qJD(3) + qJD(4);
t549 = qJD(5) + t580;
t466 = t549 * t589;
t814 = qJD(6) * t587;
t354 = -t522 * t814 + t466;
t815 = qJD(6) * t523;
t491 = qJD(1) - t815;
t99 = t233 * t895 - t236 * t418 - t240 * t419;
t33 = t100 * t353 + t125 * t491 - t354 * t99;
t101 = t233 * t894 + t420 * t236 - t240 * t421;
t102 = t235 * t894 + t420 * t238 + t421 * t241;
t126 = t306 * t894 + t308 * t420 + t310 * t421;
t34 = -t101 * t354 + t102 * t353 + t126 * t491;
t715 = rSges(7,1) * t419 - rSges(7,2) * t418;
t242 = rSges(7,3) * t895 + t715;
t973 = pkin(5) * t523;
t444 = pkin(10) * t522 + t973;
t386 = t444 * t587;
t1045 = t242 + t386;
t244 = t421 * rSges(7,1) + t420 * rSges(7,2) + rSges(7,3) * t894;
t892 = t523 * t589;
t388 = pkin(5) * t892 + pkin(10) * t894;
t1044 = t244 + t388;
t714 = rSges(7,1) * t588 - rSges(7,2) * t586;
t312 = -rSges(7,3) * t523 + t522 * t714;
t974 = pkin(5) * t522;
t443 = -pkin(10) * t523 + t974;
t1043 = t312 + t443;
t893 = t523 * t587;
t932 = Icges(6,6) * t589;
t333 = Icges(6,4) * t893 - Icges(6,2) * t895 - t932;
t516 = Icges(6,4) * t523;
t436 = Icges(6,1) * t522 + t516;
t1042 = -t436 * t587 - t333;
t701 = -Icges(6,2) * t522 + t516;
t334 = Icges(6,6) * t587 + t589 * t701;
t1041 = -t436 * t589 - t334;
t945 = Icges(6,4) * t522;
t437 = Icges(6,1) * t523 - t945;
t336 = Icges(6,5) * t587 + t437 * t589;
t434 = Icges(6,2) * t523 + t945;
t1040 = -t434 * t589 + t336;
t530 = cos(t550);
t888 = t530 * t587;
t529 = sin(t550);
t891 = t529 * t587;
t933 = Icges(5,6) * t589;
t357 = Icges(5,4) * t888 - Icges(5,2) * t891 - t933;
t519 = Icges(5,4) * t530;
t457 = Icges(5,1) * t529 + t519;
t1039 = -t457 * t587 - t357;
t702 = -Icges(5,2) * t529 + t519;
t358 = Icges(5,6) * t587 + t589 * t702;
t1038 = -t457 * t589 - t358;
t946 = Icges(5,4) * t529;
t458 = Icges(5,1) * t530 - t946;
t360 = Icges(5,5) * t587 + t458 * t589;
t455 = Icges(5,2) * t530 + t946;
t1037 = -t455 * t589 + t360;
t1036 = -t434 + t437;
t1035 = t436 + t701;
t504 = t580 * t589;
t975 = pkin(4) * t530;
t447 = t504 * t975;
t819 = qJD(1) * t587;
t784 = t529 * t819;
t487 = pkin(4) * t784;
t1034 = t447 + t487;
t1033 = t457 + t702;
t288 = t312 * t587;
t661 = t714 * t523;
t962 = rSges(7,3) * t522;
t313 = t661 + t962;
t385 = t443 * t587;
t776 = t523 * t814;
t816 = qJD(6) * t522;
t1032 = -qJD(1) * t385 + t1043 * t819 + t242 * t816 - t288 * t491 - t312 * t776 + t354 * t313 + t466 * t444;
t1031 = t33 * t587 + t34 * t589;
t818 = qJD(1) * t589;
t1030 = -t242 * t491 - t312 * t354 - t443 * t466;
t1029 = 0.2e1 * qJD(3);
t547 = sin(t579);
t1028 = rSges(4,2) * t547;
t799 = t522 * t466;
t644 = -t523 * t819 - t799;
t798 = t523 * t466;
t785 = t522 * t819;
t836 = rSges(6,2) * t785 + rSges(6,3) * t818;
t189 = rSges(6,1) * t644 - rSges(6,2) * t798 + t836;
t966 = rSges(6,1) * t523;
t442 = -rSges(6,2) * t522 + t966;
t340 = t442 * t549;
t808 = qJD(1) * qJD(3);
t537 = t589 * t808;
t807 = qJD(1) * qJD(4);
t478 = t589 * t807 + t537;
t806 = qJD(1) * qJD(5);
t439 = t589 * t806 + t478;
t441 = rSges(6,1) * t522 + rSges(6,2) * t523;
t960 = pkin(3) * qJD(3);
t802 = t547 * t960;
t976 = pkin(4) * t529;
t423 = -t580 * t976 - t802;
t400 = t589 * t423;
t817 = qJD(3) * t589;
t780 = t547 * t817;
t733 = pkin(3) * t780;
t585 = -pkin(7) - qJ(2);
t578 = -pkin(8) + t585;
t556 = -pkin(9) + t578;
t826 = -t556 + t578;
t584 = cos(pkin(11));
t531 = t584 * pkin(2) + pkin(1);
t548 = cos(t579);
t479 = pkin(3) * t548 + t531;
t417 = t479 + t975;
t843 = t417 - t479;
t145 = t733 + t400 + (-t587 * t843 + t589 * t826) * qJD(1);
t825 = -t578 + t585;
t835 = t479 - t531;
t245 = -t733 + (-t587 * t835 + t589 * t825) * qJD(1);
t541 = qJ(2) * t818;
t880 = t585 * t589;
t970 = pkin(1) - t531;
t660 = t587 * t970 - t880;
t809 = qJD(1) * qJD(2);
t554 = qJD(2) * t587;
t827 = t541 + t554;
t841 = qJD(1) * (-pkin(1) * t819 + t827) + t587 * t809;
t787 = qJD(1) * (qJD(1) * t660 - t541) + t841;
t882 = t548 * qJD(3) ^ 2;
t603 = qJD(1) * t245 + (-t537 * t547 - t587 * t882) * pkin(3) + t787;
t889 = t530 * t580;
t596 = qJD(1) * t145 + (-t478 * t529 - t503 * t889) * pkin(4) + t603;
t55 = qJD(1) * t189 - t340 * t465 - t439 * t441 + t596;
t1027 = t55 * t587;
t559 = t589 * qJ(2);
t365 = -t559 + t660;
t505 = pkin(1) * t587 - t559;
t494 = qJD(1) * t505;
t1025 = qJD(1) * t365 - t494;
t571 = t587 * rSges(4,3);
t883 = t548 * t589;
t885 = t547 * t589;
t394 = rSges(4,1) * t883 - rSges(4,2) * t885 + t571;
t512 = t589 * t531;
t736 = -t585 * t587 + t512;
t1024 = t394 + t736;
t557 = t587 * qJ(2);
t507 = t589 * pkin(1) + t557;
t965 = rSges(3,2) * sin(pkin(11));
t969 = rSges(3,1) * t584;
t674 = t587 * rSges(3,3) + (-t965 + t969) * t589;
t1023 = t507 + t674;
t527 = Icges(4,4) * t548;
t703 = -Icges(4,2) * t547 + t527;
t472 = Icges(4,1) * t547 + t527;
t1022 = t1032 + t1034;
t383 = t441 * t587;
t735 = qJD(1) * t383 - t466 * t442;
t1021 = t441 * t819 + t1034 - t735;
t362 = rSges(5,1) * t888 - rSges(5,2) * t891 - t589 * rSges(5,3);
t570 = t587 * rSges(5,3);
t887 = t530 * t589;
t890 = t529 * t589;
t363 = rSges(5,1) * t887 - rSges(5,2) * t890 + t570;
t839 = t587 * t479 + t589 * t578;
t295 = t531 * t587 - t839 + t880;
t464 = t589 * t479;
t296 = t587 * t825 + t464 - t512;
t868 = -t295 * t553 + t296 * t817;
t114 = t362 * t503 + t363 * t504 + t868;
t675 = t554 - t733;
t850 = t365 - t505;
t790 = t295 + t850;
t963 = rSges(5,2) * t530;
t460 = rSges(5,1) * t529 + t963;
t900 = t460 * t504;
t129 = -t900 + (-t362 + t790) * qJD(1) + t675;
t789 = t296 + t736;
t499 = t587 * t802;
t555 = qJD(2) * t589;
t831 = t499 + t555;
t130 = -t460 * t503 + (t363 + t789) * qJD(1) - t831;
t410 = t460 * t587;
t411 = t460 * t589;
t967 = rSges(5,1) * t530;
t461 = -rSges(5,2) * t529 + t967;
t1018 = -t129 * (qJD(1) * t410 - t504 * t461) - t114 * (-t503 * t410 - t411 * t504) - t130 * (-qJD(1) * t411 - t461 * t503);
t884 = t548 * t587;
t886 = t547 * t587;
t930 = Icges(4,3) * t589;
t376 = Icges(4,5) * t884 - Icges(4,6) * t886 - t930;
t511 = Icges(4,4) * t886;
t941 = Icges(4,5) * t589;
t380 = Icges(4,1) * t884 - t511 - t941;
t934 = Icges(4,6) * t589;
t378 = Icges(4,4) * t884 - Icges(4,2) * t886 - t934;
t909 = t378 * t547;
t682 = -t380 * t548 + t909;
t149 = -t376 * t589 - t587 * t682;
t732 = qJD(1) * t523 - qJD(6);
t1017 = t587 * t732 + t799;
t433 = Icges(6,5) * t523 - Icges(6,6) * t522;
t432 = Icges(6,5) * t522 + Icges(6,6) * t523;
t906 = t432 * t589;
t913 = t334 * t522;
t928 = Icges(6,3) * t589;
t1016 = -t549 * t906 + (-t336 * t523 - t433 * t587 + t913 + t928) * qJD(1);
t486 = Icges(6,4) * t895;
t939 = Icges(6,5) * t589;
t335 = Icges(6,1) * t893 - t486 - t939;
t686 = t333 * t522 - t335 * t523;
t332 = Icges(6,3) * t587 + t433 * t589;
t824 = qJD(1) * t332;
t907 = t432 * t587;
t1015 = qJD(1) * t686 - t549 * t907 + t824;
t454 = Icges(5,5) * t530 - Icges(5,6) * t529;
t453 = Icges(5,5) * t529 + Icges(5,6) * t530;
t902 = t453 * t589;
t911 = t358 * t529;
t929 = Icges(5,3) * t589;
t1014 = -t580 * t902 + (-t360 * t530 - t454 * t587 + t911 + t929) * qJD(1);
t498 = Icges(5,4) * t891;
t940 = Icges(5,5) * t589;
t359 = Icges(5,1) * t888 - t498 - t940;
t684 = t357 * t529 - t359 * t530;
t356 = Icges(5,3) * t587 + t454 * t589;
t822 = qJD(1) * t356;
t903 = t453 * t587;
t1013 = qJD(1) * t684 - t580 * t903 + t822;
t469 = Icges(4,5) * t548 - Icges(4,6) * t547;
t468 = Icges(4,5) * t547 + Icges(4,6) * t548;
t648 = qJD(3) * t468;
t947 = Icges(4,4) * t547;
t473 = Icges(4,1) * t548 - t947;
t381 = Icges(4,5) * t587 + t473 * t589;
t379 = Icges(4,6) * t587 + t589 * t703;
t908 = t379 * t547;
t681 = -t381 * t548 + t908;
t1012 = -t589 * t648 + (-t469 * t587 + t681 + t930) * qJD(1);
t377 = Icges(4,3) * t587 + t469 * t589;
t821 = qJD(1) * t377;
t1011 = qJD(1) * t682 - t587 * t648 + t821;
t680 = t434 * t522 - t436 * t523;
t1010 = qJD(1) * t680 + t433 * t549;
t679 = t455 * t529 - t457 * t530;
t1009 = qJD(1) * t679 + t454 * t580;
t470 = Icges(4,2) * t548 + t947;
t678 = t470 * t547 - t472 * t548;
t1008 = t678 * qJD(1) + t469 * qJD(3);
t1007 = t587 * (-t470 * t589 + t381) - t589 * (-Icges(4,2) * t884 + t380 - t511);
t676 = t491 * t589;
t199 = t1017 * t586 + t588 * t676;
t200 = -t1017 * t588 + t586 * t676;
t795 = t200 * rSges(7,1) + t199 * rSges(7,2) + rSges(7,3) * t798;
t121 = -rSges(7,3) * t785 + t795;
t713 = -rSges(7,1) * t586 - rSges(7,2) * t588;
t170 = t549 * t661 + (rSges(7,3) * t549 + qJD(6) * t713) * t522;
t462 = pkin(10) * t798;
t214 = pkin(5) * t644 - pkin(10) * t785 + t462;
t642 = -t785 + t798;
t259 = qJD(6) * t642 + t439;
t344 = t444 * t549;
t777 = t549 * t816;
t26 = qJD(1) * t214 + t491 * t121 - t353 * t170 + t244 * t777 - t259 * t312 - t465 * t344 - t439 * t443 + t596;
t881 = t549 * t587;
t643 = t522 * t818 + t523 * t881;
t800 = t522 * t881;
t201 = t491 * t878 + (-t589 * t732 + t800) * t586;
t896 = t522 * t549;
t202 = t732 * t877 + (t491 * t586 - t588 * t896) * t587;
t716 = rSges(7,1) * t202 + rSges(7,2) * t201;
t122 = rSges(7,3) * t643 + t716;
t215 = t643 * pkin(10) + (t523 * t818 - t800) * pkin(5);
t536 = t587 * t808;
t477 = t587 * t807 + t536;
t438 = t587 * t806 + t477;
t258 = qJD(6) * t643 + t438;
t539 = t589 * t809;
t662 = -pkin(3) * t589 * t882 + qJD(1) * t499 + t539;
t632 = -t580 * t447 + t477 * t976 + t662;
t515 = t556 * t819;
t737 = t423 * t587 - t515;
t518 = t578 * t819;
t832 = t499 + t518;
t146 = t818 * t843 + t737 + t832;
t521 = t585 * t819;
t246 = t818 * t835 + t521 - t832;
t451 = qJD(1) * t507 - t555;
t856 = t521 - (-t589 * t970 - t557) * qJD(1) - t451;
t791 = -t246 + t856;
t731 = -t146 + t791;
t27 = -t242 * t777 - t122 * t491 - t170 * t354 + t258 * t312 - t344 * t466 + t438 * t443 + (-t215 + t731) * qJD(1) + t632;
t397 = t589 * t417;
t248 = t587 * t826 + t397 - t464;
t728 = t248 + t789;
t445 = t503 * t976;
t786 = -t445 - t831;
t67 = t244 * t491 - t312 * t353 - t443 * t465 + (t388 + t728) * qJD(1) + t786;
t925 = qJD(1) * t67;
t1006 = (t27 + t925) * t589 + t26 * t587;
t569 = t587 * rSges(6,3);
t190 = -t549 * t383 + (t442 * t589 + t569) * qJD(1);
t56 = -t340 * t466 + t438 * t441 + (-t190 + t731) * qJD(1) + t632;
t342 = rSges(6,1) * t892 - rSges(6,2) * t894 + t569;
t92 = -t441 * t465 + (t342 + t728) * qJD(1) + t786;
t1005 = (qJD(1) * t92 + t56) * t589 + t1027;
t654 = t698 * t523;
t687 = -t308 * t586 + t310 * t588;
t689 = -t238 * t586 + t241 * t588;
t1004 = t353 * (-t306 * t589 - t689) - t354 * (-t306 * t587 + t1049) + t491 * (Icges(7,3) * t522 + t654 - t687);
t699 = -Icges(7,2) * t588 - t943;
t1003 = t353 * (-Icges(7,2) * t421 + t241 + t409) - t354 * (-Icges(7,2) * t419 - t240 - t407) + t491 * (t699 * t522 + t310);
t1002 = qJD(1) * t1035 + t465 * t1040 - t466 * (-Icges(6,2) * t893 + t335 - t486);
t1001 = qJD(1) * t1033 + t503 * t1037 - t504 * (-Icges(5,2) * t888 + t359 - t498);
t1000 = t258 / 0.2e1;
t999 = t259 / 0.2e1;
t998 = -t353 / 0.2e1;
t997 = t353 / 0.2e1;
t996 = -t354 / 0.2e1;
t995 = t354 / 0.2e1;
t994 = t438 / 0.2e1;
t993 = t439 / 0.2e1;
t992 = -t465 / 0.2e1;
t991 = t465 / 0.2e1;
t990 = -t466 / 0.2e1;
t989 = t466 / 0.2e1;
t988 = t477 / 0.2e1;
t987 = t478 / 0.2e1;
t986 = -t491 / 0.2e1;
t985 = t491 / 0.2e1;
t984 = -t503 / 0.2e1;
t983 = t503 / 0.2e1;
t982 = -t504 / 0.2e1;
t981 = t504 / 0.2e1;
t980 = t587 / 0.2e1;
t979 = -t589 / 0.2e1;
t978 = -rSges(7,3) - pkin(10);
t977 = pkin(3) * t547;
t972 = -qJD(1) / 0.2e1;
t971 = qJD(1) / 0.2e1;
t968 = rSges(4,1) * t548;
t964 = rSges(4,2) * t548;
t116 = Icges(7,5) * t202 + Icges(7,6) * t201 + Icges(7,3) * t643;
t118 = Icges(7,4) * t202 + Icges(7,2) * t201 + Icges(7,6) * t643;
t120 = Icges(7,1) * t202 + Icges(7,4) * t201 + Icges(7,5) * t643;
t24 = (-t1049 * t549 - t116) * t523 + (-t118 * t586 + t120 * t588 + t233 * t549 + (-t236 * t588 + t240 * t586) * qJD(6)) * t522;
t959 = t24 * t354;
t115 = Icges(7,5) * t200 + Icges(7,6) * t199 + Icges(7,3) * t642;
t117 = Icges(7,4) * t200 + Icges(7,2) * t199 + Icges(7,6) * t642;
t119 = Icges(7,1) * t200 + Icges(7,4) * t199 + Icges(7,5) * t642;
t25 = (t549 * t689 - t115) * t523 + (-t117 * t586 + t119 * t588 + t235 * t549 + (-t238 * t588 - t241 * t586) * qJD(6)) * t522;
t958 = t25 * t353;
t957 = t26 * t589;
t956 = t27 * t587;
t341 = rSges(6,1) * t893 - rSges(6,2) * t895 - t589 * rSges(6,3);
t446 = t504 * t976;
t659 = -t446 + t675;
t845 = -t587 * t417 - t589 * t556;
t247 = t839 + t845;
t729 = t247 + t790;
t905 = t441 * t466;
t91 = -t905 + (-t341 + t729) * qJD(1) + t659;
t952 = t587 * t91;
t132 = -t306 * t523 + t522 * t687;
t697 = -Icges(7,5) * t586 - Icges(7,6) * t588;
t165 = t549 * t654 + (Icges(7,3) * t549 + qJD(6) * t697) * t522;
t655 = t700 * t523;
t166 = t549 * t655 + (Icges(7,6) * t549 + qJD(6) * t699) * t522;
t656 = t705 * t523;
t704 = -Icges(7,1) * t586 - t942;
t167 = t549 * t656 + (Icges(7,5) * t549 + qJD(6) * t704) * t522;
t46 = (t549 * t687 - t165) * t523 + (-t166 * t586 + t167 * t588 + t306 * t549 + (-t308 * t588 - t310 * t586) * qJD(6)) * t522;
t951 = t132 * t777 + t46 * t491;
t672 = -t503 * t247 + t248 * t504 + t868;
t49 = t242 * t353 + t244 * t354 + t386 * t465 + t388 * t466 + t672;
t926 = qJD(1) * t49;
t78 = t341 * t465 + t342 * t466 + t672;
t924 = qJD(1) * t78;
t921 = t105 * t258;
t106 = -t235 * t523 + t522 * t689;
t920 = t106 * t259;
t918 = t129 * t587;
t830 = rSges(4,2) * t886 + t589 * rSges(4,3);
t393 = rSges(4,1) * t884 - t830;
t474 = rSges(4,1) * t547 + t964;
t778 = t474 * t817;
t720 = t554 - t778;
t153 = (-t393 + t850) * qJD(1) + t720;
t917 = t153 * t587;
t781 = t474 * t553;
t154 = qJD(1) * t1024 - t555 - t781;
t431 = t474 * t589;
t916 = t154 * t431;
t331 = Icges(6,5) * t893 - Icges(6,6) * t895 - t928;
t914 = t331 * t589;
t355 = Icges(5,5) * t888 - Icges(5,6) * t891 - t929;
t912 = t355 * t589;
t901 = t455 * t580;
t899 = t468 * t587;
t898 = t468 * t589;
t897 = t503 * t530;
t875 = -t170 - t344;
t874 = -t587 * t247 + t589 * t248;
t870 = -t248 - t342;
t867 = -t587 * t295 + t589 * t296;
t866 = -t587 * t331 - t335 * t892;
t865 = t587 * t332 + t336 * t892;
t863 = -t587 * t355 - t359 * t887;
t862 = t587 * t356 + t360 * t887;
t859 = -t587 * t376 - t380 * t883;
t858 = t587 * t377 + t381 * t883;
t857 = t587 * t341 + t589 * t342;
t855 = t587 * t362 + t589 * t363;
t844 = t400 + t554;
t838 = -t470 + t473;
t837 = t472 + t703;
t834 = rSges(5,2) * t784 + rSges(5,3) * t818;
t833 = rSges(4,3) * t818 + t1028 * t819;
t524 = t587 * t965;
t829 = rSges(3,3) * t818 + qJD(1) * t524;
t828 = t589 * rSges(3,3) + t524;
t467 = -t976 - t977;
t763 = t467 + t977;
t345 = t763 * t587;
t823 = qJD(1) * t345;
t820 = qJD(1) * t469;
t156 = -t587 * t680 - t906;
t812 = t156 * qJD(1);
t171 = -t587 * t679 - t902;
t811 = t171 * qJD(1);
t197 = -t587 * t678 - t898;
t810 = t197 * qJD(1);
t805 = pkin(4) * t889;
t804 = t587 * t969;
t801 = t548 * t960;
t797 = t589 * t145 + t587 * t146 - t247 * t818;
t796 = t589 * t189 + t587 * t190 + t341 * t818;
t228 = -t504 * t963 + (-t504 * t529 - t530 * t819) * rSges(5,1) + t834;
t229 = -t580 * t410 + (t461 * t589 + t570) * qJD(1);
t794 = t589 * t228 + t587 * t229 + t362 * t818;
t793 = t589 * t245 + t587 * t246 - t295 * t818;
t792 = -t248 - t1044;
t782 = t547 * t818;
t775 = t523 * t813;
t772 = -pkin(1) - t969;
t771 = t819 / 0.2e1;
t770 = t818 / 0.2e1;
t769 = -t553 / 0.2e1;
t766 = t817 / 0.2e1;
t764 = -t460 - t977;
t760 = -t531 - t968;
t759 = t547 * (-t587 ^ 2 - t589 ^ 2);
t758 = (-t587 * t701 + t932) * qJD(1) + t1040 * t549;
t757 = qJD(1) * t334 + t335 * t549 - t434 * t881;
t756 = (-t437 * t587 + t939) * qJD(1) + t1041 * t549;
t755 = qJD(1) * t336 + t1042 * t549;
t754 = (-t587 * t702 + t933) * qJD(1) + t1037 * t580;
t753 = qJD(1) * t358 + t359 * t580 - t587 * t901;
t752 = (-t458 * t587 + t940) * qJD(1) + t1038 * t580;
t751 = qJD(1) * t360 + t1039 * t580;
t290 = t336 * t893;
t750 = t332 * t589 - t290;
t384 = t441 * t589;
t749 = -t465 * t383 - t384 * t466;
t346 = t763 * t589;
t748 = t503 * t345 + t346 * t504;
t301 = t360 * t888;
t747 = t356 * t589 - t301;
t317 = t381 * t884;
t746 = t377 * t589 - t317;
t745 = -t331 + t913;
t744 = t1035 * t549;
t743 = t1036 * t549;
t742 = -t355 + t911;
t741 = -qJD(1) * t384 - t442 * t465;
t740 = -t376 + t908;
t739 = t1033 * t580;
t738 = t458 * t580 - t901;
t730 = t857 + t874;
t727 = t1044 * t589 + t1045 * t587;
t726 = -pkin(4) * t897 + qJD(1) * t346;
t725 = -t340 - t805;
t398 = t461 * t580;
t721 = -t398 - t801;
t717 = t968 - t1028;
t66 = (-t386 + t729) * qJD(1) + t659 + t1030;
t712 = t587 * t67 + t589 * t66;
t711 = -t587 * t92 - t589 * t91;
t710 = t100 * t589 + t587 * t99;
t709 = t100 * t587 - t589 * t99;
t696 = t101 * t589 - t102 * t587;
t695 = t101 * t587 + t102 * t589;
t694 = t105 * t589 - t106 * t587;
t693 = t105 * t587 + t106 * t589;
t692 = -t129 * t589 - t130 * t587;
t691 = -t153 * t589 - t154 * t587;
t688 = t242 * t589 - t244 * t587;
t168 = t333 * t523 + t335 * t522;
t178 = t357 * t530 + t359 * t529;
t209 = t378 * t548 + t380 * t547;
t210 = t379 * t548 + t381 * t547;
t677 = -t805 + t875;
t670 = t1045 * t818 + (t121 + t214) * t589 + (t122 + t215) * t587;
t669 = t796 + t797;
t668 = t727 + t874;
t430 = t474 * t587;
t653 = t686 * t587;
t652 = t684 * t587;
t651 = -t801 - t805;
t650 = qJD(3) * t472;
t649 = qJD(3) * t470;
t150 = -t379 * t886 - t746;
t647 = (-t149 * t589 + t150 * t587) * qJD(3);
t151 = -t378 * t885 - t859;
t152 = -t379 * t885 + t858;
t646 = (-t151 * t589 + t152 * t587) * qJD(3);
t230 = (t393 * t587 + t394 * t589) * qJD(3);
t645 = -t444 - t962;
t641 = (-t529 * t818 - t897) * pkin(4);
t640 = -t233 * t354 + t235 * t353 + t306 * t491;
t639 = (-Icges(7,5) * t418 - Icges(7,6) * t419) * t354 - (Icges(7,5) * t420 - Icges(7,6) * t421) * t353 - t697 * t522 * t491;
t638 = -t340 + t651;
t637 = qJD(1) * t433 - t465 * t906 + t466 * t907;
t636 = qJD(1) * t454 - t503 * t902 + t504 * t903;
t635 = t245 * t817 + t246 * t553 - t295 * t537 - t296 * t536;
t289 = t312 * t589;
t387 = t443 * t589;
t634 = t242 * t775 - t353 * t288 - t289 * t354 - t465 * t385 - t387 * t466;
t633 = t378 * t589 - t379 * t587;
t631 = qJD(1) * t295 + t1025 + t675;
t630 = t651 + t875;
t629 = t670 + t797;
t628 = t522 * t639;
t623 = (-t547 * t837 + t548 * t838) * qJD(1);
t620 = (Icges(7,1) * t420 - t238 - t944) * t353 - (-Icges(7,1) * t418 - t236 - t408) * t354 + (t704 * t522 - t308) * t491;
t617 = qJD(1) * t1036 + t1041 * t465 - t1042 * t466;
t616 = t1038 * t503 - t1039 * t504 + (-t455 + t458) * qJD(1);
t615 = qJD(1) * t247 - t446 + t631;
t614 = -t244 * t776 + t634;
t613 = -qJD(1) * t387 + t244 * t816 - t491 * t289 - t312 * t775 - t313 * t353 - t444 * t465;
t609 = qJD(1) * t331 - t522 * t757 + t523 * t755;
t608 = -t522 * t758 + t523 * t756 + t824;
t607 = qJD(1) * t355 - t529 * t753 + t530 * t751;
t606 = -t529 * t754 + t530 * t752 + t822;
t605 = qJD(1) * t432 - t522 * t744 + t523 * t743;
t604 = qJD(1) * t453 - t529 * t739 + t530 * t738;
t602 = t145 * t504 + t503 * t146 - t478 * t247 - t248 * t477 + t635;
t255 = qJD(1) * t379 - t587 * t649;
t257 = qJD(1) * t381 - t587 * t650;
t601 = qJD(1) * t376 - qJD(3) * t209 - t255 * t547 + t257 * t548;
t254 = -t589 * t649 + (-t587 * t703 + t934) * qJD(1);
t256 = -t589 * t650 + (-t473 * t587 + t941) * qJD(1);
t600 = -qJD(3) * t210 - t254 * t547 + t256 * t548 + t821;
t449 = t703 * qJD(3);
t450 = t473 * qJD(3);
t599 = qJD(1) * t468 - t449 * t547 + t450 * t548 + (-t470 * t548 - t472 * t547) * qJD(3);
t598 = -t1007 * t547 + t633 * t548;
t133 = -t653 - t914;
t134 = -t334 * t895 - t750;
t135 = -t333 * t894 - t866;
t136 = -t334 * t894 + t865;
t169 = t334 * t523 + t336 * t522;
t20 = t116 * t894 + t118 * t420 + t120 * t421 + t199 * t236 - t200 * t240 + t233 * t642;
t21 = t115 * t894 + t117 * t420 + t119 * t421 + t199 * t238 + t200 * t241 + t235 * t642;
t22 = t116 * t895 - t118 * t418 + t120 * t419 + t201 * t236 - t202 * t240 + t233 * t643;
t23 = t115 * t895 - t117 * t418 + t119 * t419 + t201 * t238 + t202 * t241 + t235 * t643;
t283 = t308 * t587;
t284 = t308 * t589;
t285 = t310 * t587;
t286 = t310 * t589;
t39 = t165 * t894 + t166 * t420 + t167 * t421 + t199 * t308 + t200 * t310 + t306 * t642;
t3 = t101 * t258 + t102 * t259 + t126 * t777 - t20 * t354 + t21 * t353 + t39 * t491;
t309 = Icges(7,6) * t522 + t655;
t311 = Icges(7,5) * t522 + t656;
t37 = -t105 * t354 + t106 * t353 + t132 * t491;
t40 = t165 * t895 - t166 * t418 + t167 * t419 + t201 * t308 + t202 * t310 + t306 * t643;
t4 = t100 * t259 + t125 * t777 - t22 * t354 + t23 * t353 + t258 * t99 + t40 * t491;
t50 = t1015 * t587 + t609 * t589;
t51 = t1016 * t587 + t608 * t589;
t52 = -t1015 * t589 + t609 * t587;
t53 = -t1016 * t589 + t608 * t587;
t594 = -t1002 * t522 + t617 * t523;
t595 = t1004 * t522;
t70 = -t133 * t466 + t134 * t465 + t812;
t157 = -t589 * t680 + t907;
t155 = t157 * qJD(1);
t71 = -t135 * t466 + t136 * t465 + t155;
t81 = t1010 * t587 + t605 * t589;
t82 = -t1010 * t589 + t605 * t587;
t93 = t522 * t755 + t523 * t757;
t94 = t522 * t756 + t523 * t758;
t597 = (((t284 * t586 - t286 * t588 + t235) * t353 - (t283 * t586 - t285 * t588 + t233) * t354 + (-t309 * t586 + t311 * t588 + t306) * t491 + t132 * qJD(6)) * t522 + (qJD(6) * t693 - t1004) * t523) * t986 - t37 * t816 / 0.2e1 - t259 * t696 / 0.2e1 + (t1002 * t523 + t617 * t522) * t972 - t694 * t777 / 0.2e1 + (t587 * t94 - t589 * t93 + (t168 * t587 + t169 * t589) * qJD(1)) * t971 + (qJD(1) * t693 - t24 * t589 + t25 * t587) * t985 + (t587 * t594 - t589 * t637) * t989 + (-t52 * t589 + t53 * t587 + (t133 * t587 + t134 * t589) * qJD(1)) * t990 + (-t50 * t589 + t51 * t587 + (t135 * t587 + t136 * t589) * qJD(1)) * t991 + (t587 * t637 + t589 * t594) * t992 + (-t135 * t589 + t136 * t587) * t993 + (-t133 * t589 + t134 * t587) * t994 + ((t284 * t418 - t286 * t419) * t353 - (t283 * t418 - t285 * t419) * t354 + (-t309 * t418 + t311 * t419) * t491 + (t100 * t892 + t125 * t522) * qJD(6) + ((qJD(6) * t99 + t640) * t523 + t595) * t587) * t995 + (qJD(1) * t710 - t22 * t589 + t23 * t587) * t996 + (qJD(1) * t695 - t20 * t589 + t21 * t587) * t997 + ((-t284 * t420 - t286 * t421) * t353 - (-t283 * t420 - t285 * t421) * t354 + (t309 * t420 + t311 * t421) * t491 + (t101 * t893 + t126 * t522) * qJD(6) + ((qJD(6) * t102 + t640) * t523 + t595) * t589) * t998 + t709 * t1000 + (qJD(1) * t81 + t135 * t438 + t136 * t439 + t465 * t51 - t466 * t50 + t3) * t980 + (qJD(1) * t82 + t133 * t438 + t134 * t439 + t465 * t53 - t466 * t52 + t4) * t979 + (t70 + t33) * t771 + (t71 + t34) * t770 - t1031 * t815 / 0.2e1;
t593 = -t1001 * t529 + t616 * t530;
t103 = t529 * t751 + t530 * t753;
t104 = t529 * t752 + t530 * t754;
t138 = -t652 - t912;
t139 = -t358 * t891 - t747;
t140 = -t357 * t890 - t863;
t141 = -t358 * t890 + t862;
t179 = t358 * t530 + t360 * t529;
t59 = t1013 * t587 + t607 * t589;
t60 = t1014 * t587 + t606 * t589;
t61 = -t1013 * t589 + t607 * t587;
t62 = -t1014 * t589 + t606 * t587;
t79 = -t138 * t504 + t139 * t503 + t811;
t172 = -t589 * t679 + t903;
t164 = t172 * qJD(1);
t80 = -t140 * t504 + t141 * t503 + t164;
t95 = t1009 * t587 + t604 * t589;
t96 = -t1009 * t589 + t604 * t587;
t592 = t597 + (t587 * t62 - t589 * t61 + (t138 * t587 + t139 * t589) * qJD(1)) * t982 + (qJD(1) * t96 + t138 * t477 + t139 * t478 + t503 * t62 - t504 * t61) * t979 + t80 * t770 + (t587 * t636 + t589 * t593) * t984 + (-t140 * t589 + t141 * t587) * t987 + (t1001 * t530 + t616 * t529) * t972 + (t587 * t593 - t589 * t636) * t981 + (qJD(1) * t95 + t140 * t477 + t141 * t478 + t503 * t60 - t504 * t59) * t980 + (t587 * t60 - t589 * t59 + (t140 * t587 + t141 * t589) * qJD(1)) * t983 + (-t103 * t589 + t104 * t587 + (t178 * t587 + t179 * t589) * qJD(1)) * t971 + t79 * t771 + (-t138 * t589 + t139 * t587) * t988;
t452 = t717 * qJD(3);
t422 = t804 - t828;
t416 = t589 * t446;
t382 = t713 * t522;
t298 = qJD(1) * t1023 - t555;
t297 = t554 + (-t422 - t505) * qJD(1);
t271 = rSges(7,1) * t420 - rSges(7,2) * t421;
t270 = -rSges(7,1) * t418 - rSges(7,2) * t419;
t261 = -qJD(3) * t430 + (t589 * t717 + t571) * qJD(1);
t260 = -t817 * t964 + (-t548 * t819 - t780) * rSges(4,1) + t833;
t251 = t539 + (-qJD(1) * t674 - t451) * qJD(1);
t250 = qJD(1) * (-qJD(1) * t804 + t829) + t841;
t198 = -t589 * t678 + t899;
t193 = t198 * qJD(1);
t128 = -t452 * t817 + t539 + (-t261 + t781 + t856) * qJD(1);
t127 = -t452 * t553 + (t260 - t778) * qJD(1) + t787;
t110 = -qJD(3) * t681 + t254 * t548 + t256 * t547;
t109 = -qJD(3) * t682 + t255 * t548 + t257 * t547;
t108 = -t1008 * t589 + t599 * t587;
t107 = t1008 * t587 + t599 * t589;
t86 = t193 + t646;
t85 = t647 + t810;
t84 = -t398 * t504 + t460 * t477 + (-t229 + t791) * qJD(1) + t662;
t83 = qJD(1) * t228 - t398 * t503 - t460 * t478 + t603;
t54 = t228 * t504 + t229 * t503 + t362 * t478 - t363 * t477 + t635;
t28 = t189 * t466 + t190 * t465 + t341 * t439 - t342 * t438 + t602;
t13 = t121 * t354 + t122 * t353 + t214 * t466 + t215 * t465 + t242 * t259 - t244 * t258 + t386 * t439 - t388 * t438 + t602;
t1 = [(-(-qJD(1) * t341 + t615 - t905 - t91) * t92 + t56 * (-t341 + t845) + t91 * (t555 - t737) + t55 * (-t556 * t587 + t342 + t397) + t92 * (t836 + t844) + (-t384 * t92 + t441 * t952) * t549 + ((-t91 * rSges(6,3) + t92 * (-t417 - t966)) * t587 + (t91 * (-t417 - t442) - t92 * t556) * t589) * qJD(1)) * m(6) + (t85 - t810 + ((t589 * t740 + t152 - t858) * t589 + (t587 * t740 + t151 + t746) * t587) * qJD(3)) * t769 + (t157 + t169) * t993 + (t27 * (-t715 + t845) + t26 * (t397 + t1044) + (t27 * t645 - t26 * t556 + (t522 * t978 - t417 - t973) * t925) * t587 + (t515 + t555 - t716 + (-t417 + t645) * t818 + (-t423 + (t523 * t978 + t974) * t549) * t587) * t66 + (qJD(1) * t386 - t615 + t66 + t462 + t795 + t844 + (-pkin(5) * t896 - qJD(1) * t556) * t589 - t1030) * t67) * m(7) + (t94 + t81) * t991 - (t109 + t108 + t86) * t817 / 0.2e1 + (t110 + t107) * t553 / 0.2e1 + (t172 + t179) * t987 + (-t811 + (t141 - t652 - t862) * t504 + (t587 * t742 + t140 - t301) * t503 + ((t356 + t684) * t503 + t742 * t504) * t589 + t79) * t984 + (t251 * (t587 * t772 + t559 + t828) + t297 * t555 + t250 * t1023 + t298 * (t827 + t829) + (t297 * (t772 + t965) * t589 + (t297 * (-rSges(3,3) - qJ(2)) + t298 * t772) * t587) * qJD(1) - (-qJD(1) * t422 - t297 - t494 + t554) * t298) * m(3) + (t104 + t95) * t983 + (t156 + t168) * t994 + (t40 + t34) * t996 + (t193 + ((t150 - t317 + (t377 + t909) * t589 + t859) * t589 + t858 * t587) * qJD(3)) * t766 + ((t209 + t197) * t587 + (t210 + t198) * t589) * t808 / 0.2e1 + t921 / 0.2e1 + t920 / 0.2e1 + t958 / 0.2e1 - t959 / 0.2e1 + (-(-qJD(1) * t362 - t129 + t631 - t900) * t130 + t84 * (-t362 - t839) + t129 * (t518 + t831) + t83 * (-t578 * t587 + t363 + t464) + t130 * (t675 + t834) + (-t130 * t411 + t460 * t918) * t580 + ((-t129 * rSges(5,3) + t130 * (-t479 - t967)) * t587 + (t129 * (-t461 - t479) - t130 * t578) * t589) * qJD(1)) * m(5) + t951 + (t103 + t96 + t80) * t982 + (t171 + t178) * t988 + t34 * t995 + (-t812 + (t136 - t653 - t865) * t466 + (t745 * t587 + t135 - t290) * t465 + ((t332 + t686) * t465 + t745 * t466) * t589 + t70) * t992 + (-qJD(3) * t678 + t449 * t548 + t450 * t547 + t522 * t743 + t523 * t744 + t529 * t738 + t530 * t739) * qJD(1) + (t93 + t82 + t71) * t990 + (-(-qJD(1) * t393 + t1025 - t153 + t720) * t154 + t128 * (t587 * t760 + t830 - t880) + t153 * (t521 + t555) + t127 * t1024 + t154 * (t554 + t833) + (t474 * t917 - t916) * qJD(3) + ((-t153 * rSges(4,3) + t154 * t760) * t587 + (t153 * (-t531 - t717) - t154 * t585) * t589) * qJD(1)) * m(4) + (t164 + (t139 + (t357 * t589 + t358 * t587) * t529 + t747 + t863) * t504 + (-t359 * t888 + t912 + t138 + (t357 * t587 - t358 * t589) * t529 + t862) * t503) * t981 + (t155 + (t134 + (t333 * t589 + t334 * t587) * t522 + t750 + t866) * t466 + (-t335 * t893 + t914 + t133 + (t333 * t587 - t334 * t589) * t522 + t865) * t465) * t989 + t39 * t997 + t126 * t999 + t125 * t1000; 0.2e1 * (-t957 / 0.2e1 + t956 / 0.2e1) * m(7) + 0.2e1 * (t55 * t979 + t56 * t980) * m(6) + 0.2e1 * (t83 * t979 + t84 * t980) * m(5) + 0.2e1 * (t127 * t979 + t128 * t980) * m(4) + 0.2e1 * (t250 * t979 + t251 * t980) * m(3); t592 + ((-t817 * t899 - t820) * t589 + (t623 + (t589 * t898 + t598) * qJD(3)) * t587) * t766 + ((t547 * t838 + t548 * t837) * qJD(1) + (t1007 * t548 + t633 * t547) * qJD(3)) * t972 + (-t109 * t589 + t110 * t587 + (t209 * t587 + t210 * t589) * qJD(1)) * t971 + ((-t553 * t898 + t820) * t587 + (t623 + (t587 * t899 + t598) * qJD(3)) * t589) * t769 + (qJD(1) * t107 + (t587 * (t1012 * t587 + t600 * t589) - t589 * (t1011 * t587 + t601 * t589) + (t151 * t587 + t152 * t589) * qJD(1)) * t1029) * t980 + (qJD(1) * t108 + (t587 * (-t1012 * t589 + t600 * t587) - t589 * (-t1011 * t589 + t601 * t587) + (t149 * t587 + t150 * t589) * qJD(1)) * t1029) * t979 + (t647 + t85) * t771 + (t86 + t646) * t770 + (-t67 * (t613 + t726) - t49 * (t614 + t748) - (-t67 * t782 + (t49 * t759 - t548 * t712) * qJD(3)) * pkin(3) + t13 * (t668 + t867) + t49 * (t629 + t793) + (t67 * t630 + (-t296 + t792) * t926) * t587 + t1006 * (t467 - t1043) + (t589 * t630 + t1022 + t823) * t66) * m(7) + (-t92 * (t726 + t741) - t78 * (t748 + t749) - (-t92 * t782 + (t548 * t711 + t759 * t78) * qJD(3)) * pkin(3) + t28 * (t730 + t867) + t78 * (t669 + t793) + (t92 * t638 + (-t296 + t870) * t924) * t587 + t1005 * (-t441 + t467) + (t589 * t638 + t1021 + t823) * t91) * m(6) + (-(-t130 * t782 + (t114 * t759 + t548 * t692) * qJD(3)) * pkin(3) + t54 * (t855 + t867) + t114 * (t793 + t794) + (t129 * t721 + (qJD(1) * t130 + t84) * t764) * t589 + (t83 * t764 + t130 * t721 + (t129 * t460 + t114 * (-t296 - t363)) * qJD(1)) * t587 + t1018) * m(5) + (0.2e1 * t230 * (t260 * t589 + t261 * t587 + (t393 * t589 - t394 * t587) * qJD(1)) + t691 * t452 + (-t127 * t587 - t128 * t589 + (-t154 * t589 + t917) * qJD(1)) * t474 - (t153 * t430 - t916) * qJD(1) - (t230 * (-t430 * t587 - t431 * t589) + t691 * t717) * qJD(3)) * m(4); t592 + (t13 * t668 + (t67 * t677 + t792 * t926) * t587 - t67 * (t641 + t613) + t1006 * (-t1043 - t976) + (t589 * t677 + t1022 - t487) * t66 + (t629 + t416 - (-t244 * t815 - t445) * t587 - t634) * t49) * m(7) + (t28 * t730 + (t725 * t92 + t870 * t924) * t587 - t92 * (t641 + t741) + t1005 * (-t441 - t976) + (t589 * t725 + t1021 - t487) * t91 + (t445 * t587 + t416 + t669 - t749) * t78) * m(6) + (t54 * t855 + t114 * (-t363 * t819 + t794) + t692 * t398 + (-t83 * t587 - t84 * t589 + (-t130 * t589 + t918) * qJD(1)) * t460 + t1018) * m(5); t597 + (-t613 * t67 + t13 * t727 + (-t1044 * t926 + t67 * t875) * t587 - t1006 * t1043 + (t589 * t875 + t1032) * t66 + (-t614 + t670) * t49) * m(7) + (-t735 * t91 - t741 * t92 + t28 * t857 + t711 * t340 + (-t1027 - t56 * t589 + (-t589 * t92 + t952) * qJD(1)) * t441 + (-t342 * t819 - t749 + t796) * t78) * m(6); -t34 * t785 / 0.2e1 + t3 * t894 / 0.2e1 + (-t126 * t523 + t522 * t695) * t999 + ((t549 * t695 - t39) * t523 + (qJD(1) * t696 + t126 * t549 + t20 * t587 + t21 * t589) * t522) * t997 + t522 * t33 * t770 + t4 * t895 / 0.2e1 + (-t125 * t523 + t522 * t710) * t1000 + ((t549 * t710 - t40) * t523 + (-qJD(1) * t709 + t125 * t549 + t22 * t587 + t23 * t589) * t522) * t996 - t523 * (t920 + t921 + t951 + t958 - t959) / 0.2e1 + ((t549 * t693 - t46) * t523 + (qJD(1) * t694 + t132 * t549 + t24 * t587 + t25 * t589) * t522) * t985 + (t1003 * t420 + t620 * t421 - t589 * t628) * t998 + (-t1003 * t418 + t419 * t620 - t587 * t628) * t995 + (t639 * t523 + (-t1003 * t586 + t588 * t620) * t522) * t986 + (t37 + qJD(6) * (-t132 * t523 + t522 * t693)) * t896 / 0.2e1 + t1031 * t523 * t549 / 0.2e1 + ((-t67 * t121 + t66 * t122 + t27 * t242 - t26 * t244 + (t49 * t688 + (t587 * t66 - t589 * t67) * t312) * t549) * t523 + (t66 * (t170 * t587 - t242 * t549) + t67 * (-t170 * t589 + t244 * t549) + t13 * t688 + t49 * (-t121 * t587 + t122 * t589 - t242 * t819 - t244 * t818) + (qJD(1) * t712 + t956 - t957) * t312) * t522 - t66 * (-t270 * t491 - t354 * t382) - t67 * (t271 * t491 - t353 * t382) - t49 * (t270 * t353 + t271 * t354)) * m(7);];
tauc  = t1(:);
