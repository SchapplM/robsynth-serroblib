% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRPPR8_coriolisvecJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPPR8_coriolisvecJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:58:47
% EndTime: 2019-03-09 02:59:48
% DurationCPUTime: 56.79s
% Computational Cost: add. (11412->1112), mult. (29898->1372), div. (0->0), fcn. (25850->6), ass. (0->558)
t480 = cos(qJ(3));
t481 = cos(qJ(1));
t757 = t480 * t481;
t477 = sin(qJ(3));
t762 = t477 * t481;
t478 = sin(qJ(1));
t778 = Icges(5,6) * t478;
t236 = Icges(5,5) * t762 - Icges(5,3) * t757 - t778;
t791 = Icges(4,4) * t477;
t587 = Icges(4,2) * t480 + t791;
t246 = -Icges(4,6) * t478 + t481 * t587;
t969 = t236 - t246;
t423 = Icges(5,5) * t757;
t789 = Icges(5,4) * t478;
t252 = Icges(5,1) * t762 - t423 - t789;
t790 = Icges(4,4) * t480;
t591 = Icges(4,1) * t477 + t790;
t254 = -Icges(4,5) * t478 + t481 * t591;
t968 = t252 + t254;
t783 = Icges(5,5) * t477;
t371 = Icges(5,1) * t480 + t783;
t373 = Icges(4,1) * t480 - t791;
t464 = Icges(6,4) * t477;
t859 = -Icges(6,2) * t480 + t464;
t964 = -t373 + t859;
t937 = t371 - t964;
t583 = Icges(4,5) * t477 + Icges(4,6) * t480;
t237 = Icges(4,3) * t481 + t478 * t583;
t586 = Icges(5,4) * t477 - Icges(5,6) * t480;
t243 = Icges(5,2) * t481 + t478 * t586;
t967 = -t237 - t243;
t966 = (Icges(5,4) + Icges(4,5)) * t480 + (-Icges(4,6) + Icges(5,6)) * t477;
t788 = Icges(6,4) * t480;
t368 = Icges(6,1) * t477 - t788;
t462 = Icges(5,5) * t480;
t860 = Icges(5,3) * t477 + t462;
t965 = t368 + t860;
t963 = t968 * t477 - t480 * t969;
t363 = Icges(6,2) * t477 + t788;
t590 = Icges(5,1) * t477 - t462;
t962 = -t363 - t590 - t591;
t356 = Icges(6,5) * t477 - Icges(6,6) * t480;
t961 = t356 - t966;
t581 = -Icges(5,3) * t480 + t783;
t589 = Icges(6,1) * t480 + t464;
t960 = t581 - t587 - t589;
t763 = t477 * t478;
t422 = Icges(5,5) * t763;
t759 = t478 * t480;
t777 = Icges(5,6) * t481;
t235 = -Icges(5,3) * t759 + t422 + t777;
t245 = Icges(4,6) * t481 + t478 * t587;
t249 = -Icges(6,5) * t481 - t478 * t589;
t951 = t249 - t245;
t959 = -t235 - t951;
t250 = -Icges(6,5) * t478 + t481 * t589;
t952 = t246 + t250;
t958 = -t236 + t952;
t238 = -Icges(4,3) * t478 + t481 * t583;
t780 = Icges(5,2) * t478;
t244 = Icges(5,4) * t762 - Icges(5,6) * t757 - t780;
t957 = -t238 - t244;
t956 = t937 * t481;
t241 = -Icges(6,6) * t481 - t363 * t478;
t745 = t241 * t762 + t249 * t757;
t904 = -t235 * t757 + t478 * t967;
t955 = t745 - t904;
t367 = -Icges(4,2) * t477 + t790;
t560 = t367 * t480 + t373 * t477;
t561 = t368 * t480 + t477 * t859;
t936 = t371 * t477 - t480 * t860 + t560 - t561;
t908 = t957 * t481 + t759 * t969 - t968 * t763;
t242 = -Icges(6,6) * t478 + t363 * t481;
t568 = t242 * t477 + t250 * t480;
t357 = Icges(6,5) * t480 + Icges(6,6) * t477;
t234 = -Icges(6,3) * t478 + t357 * t481;
t769 = t234 * t481;
t97 = -t478 * t568 - t769;
t954 = t97 + t908;
t953 = -t242 - t254;
t251 = Icges(5,4) * t481 + t478 * t590;
t950 = t251 - t241;
t949 = t965 * t478;
t948 = t966 * t478;
t947 = t964 * t478;
t945 = (t367 - t368) * t481;
t944 = t963 * t481;
t424 = Icges(4,4) * t759;
t784 = Icges(4,5) * t481;
t253 = Icges(4,1) * t763 + t424 + t784;
t943 = -t253 - t950;
t942 = -t252 + t953;
t941 = t960 * qJD(3);
t940 = t962 * qJD(3);
t939 = -t357 - t583 - t586;
t938 = -t367 + t965;
t886 = t961 * t481;
t935 = t371 * t762 + t481 * t560 - t757 * t860 - t948;
t567 = t245 * t480 + t253 * t477;
t902 = t251 * t477 + t567;
t934 = (-Icges(5,1) * t759 - t422 + t947 + t959) * t481 + (t956 - t958) * t478;
t746 = t481 * t243 + t251 * t763;
t92 = -t235 * t759 + t746;
t94 = t481 * t237 + t245 * t759 + t253 * t763;
t233 = -Icges(6,3) * t481 - t357 * t478;
t570 = t241 * t477 + t249 * t480;
t536 = t570 * t478;
t96 = -t233 * t481 - t536;
t933 = t92 + t94 + t96;
t932 = -t233 * t478 - t251 * t762 - t481 * t567 + t955;
t893 = t481 * t568 + t944 + (-t234 + t957) * t478;
t931 = t478 * t936 - t886;
t476 = sin(qJ(6));
t479 = cos(qJ(6));
t760 = t478 * t479;
t286 = t476 * t757 + t760;
t287 = -t478 * t476 + t479 * t757;
t143 = -Icges(7,5) * t287 + Icges(7,6) * t286 + Icges(7,3) * t762;
t269 = Icges(7,4) * t287;
t145 = -Icges(7,2) * t286 - Icges(7,6) * t762 + t269;
t268 = Icges(7,4) * t286;
t149 = -Icges(7,1) * t287 + Icges(7,5) * t762 + t268;
t909 = t145 * t476 + t149 * t479;
t54 = -t480 * t143 - t477 * t909;
t688 = qJD(3) * t481;
t929 = -t860 * t688 + t945 * qJD(3) + (t478 * t581 + t777 + t951) * qJD(1);
t690 = qJD(3) * t478;
t928 = -t367 * t690 + t949 * qJD(3) + (t481 * t581 - t778 - t952) * qJD(1);
t927 = -t956 * qJD(3) + (t478 * t591 + t784 + t950) * qJD(1);
t926 = -t371 * t690 + t947 * qJD(3) + (-t481 * t590 + t789 + t953) * qJD(1);
t289 = t356 * t478;
t892 = -t289 + t948;
t925 = t936 * qJD(1) + qJD(3) * t939;
t924 = (Icges(4,2) * t763 - t424 + t943 + t949) * t481 + (-Icges(5,3) * t762 - t423 - t942 + t945) * t478;
t923 = t568 + t963;
t768 = t235 * t480;
t922 = t570 + t768 - t902;
t889 = t477 * t958 + t942 * t480;
t888 = t477 * t959 + t943 * t480;
t921 = -t938 - t962;
t920 = -t937 - t960;
t919 = t941 * t480 + t940 * t477 + (t477 * t938 + t480 * t937) * qJD(3) + t961 * qJD(1);
t918 = (t234 + t238) * qJD(1);
t917 = (t477 * t921 + t480 * t920) * qJD(1);
t916 = t931 * qJD(1);
t915 = (t478 * t893 + t481 * t932) * qJD(3);
t914 = (t478 * t954 + t933 * t481) * qJD(3);
t138 = t481 * t561 - t289;
t913 = (t138 - t935) * qJD(1);
t912 = (-t233 - t967) * qJD(1);
t684 = qJD(6) * t477;
t337 = t478 * t684 + t688;
t683 = qJD(6) * t480;
t431 = qJD(1) + t683;
t758 = t479 * t481;
t284 = t476 * t759 - t758;
t285 = t476 * t481 + t479 * t759;
t141 = -Icges(7,5) * t285 + Icges(7,6) * t284 + Icges(7,3) * t763;
t787 = Icges(7,4) * t285;
t144 = Icges(7,2) * t284 + Icges(7,6) * t763 - t787;
t267 = Icges(7,4) * t284;
t147 = -Icges(7,1) * t285 + Icges(7,5) * t763 + t267;
t49 = t141 * t763 + t284 * t144 - t285 * t147;
t910 = -t284 * t145 - t149 * t285;
t50 = -t143 * t763 - t910;
t682 = qJD(6) * t481;
t543 = t477 * t682 - t690;
t580 = Icges(7,5) * t479 - Icges(7,6) * t476;
t231 = Icges(7,3) * t480 + t477 * t580;
t785 = Icges(7,4) * t479;
t584 = -Icges(7,2) * t476 + t785;
t239 = Icges(7,6) * t480 + t477 * t584;
t786 = Icges(7,4) * t476;
t588 = Icges(7,1) * t479 - t786;
t247 = Icges(7,5) * t480 + t477 * t588;
t86 = t231 * t763 + t239 * t284 - t247 * t285;
t12 = t337 * t49 + t86 * t431 - t50 * t543;
t911 = t477 * t934 + t924 * t480;
t52 = t143 * t762 - t145 * t286 - t149 * t287;
t601 = rSges(7,1) * t479 - rSges(7,2) * t476;
t257 = rSges(7,3) * t480 + t477 * t601;
t906 = t257 * t543;
t470 = t481 * rSges(4,3);
t259 = rSges(4,1) * t763 + rSges(4,2) * t759 + t470;
t387 = t481 * pkin(1) + t478 * qJ(2);
t472 = t481 * pkin(7);
t858 = t472 + t387;
t903 = t259 + t858;
t689 = qJD(3) * t480;
t643 = t478 * t689;
t692 = qJD(1) * t481;
t649 = t477 * t692;
t528 = t643 + t649;
t644 = t477 * t688;
t693 = qJD(1) * t480;
t901 = t478 * t693 + t644;
t900 = t914 + t916;
t899 = t913 + t915;
t898 = qJD(3) * t922 + t477 * t928 - t480 * t926;
t897 = qJD(3) * t923 + t477 * t929 + t480 * t927;
t896 = t478 * t925 - t481 * t919;
t895 = t478 * t919 + t481 * t925;
t891 = t918 + t892 * qJD(3) + (t481 * t586 - t780 - t922) * qJD(1);
t890 = -qJD(1) * t923 + qJD(3) * t886 + t912;
t887 = t939 * qJD(1);
t884 = qJD(1) * t244 + qJD(3) * t889 + t477 * t927 - t480 * t929 + t918;
t883 = qJD(3) * t888 + t477 * t926 + t480 * t928 + t912;
t153 = rSges(7,1) * t287 - rSges(7,2) * t286 - rSges(7,3) * t762;
t819 = pkin(5) * t477;
t389 = pkin(8) * t480 + t819;
t882 = -t153 * t431 + t389 * t690;
t499 = t231 * t762 + t239 * t286 - t247 * t287;
t881 = t499 * t431 + t543 * t52;
t879 = 0.2e1 * qJD(3);
t436 = pkin(8) * t762;
t323 = pkin(5) * t757 - t436;
t406 = pkin(4) * t643;
t450 = qJD(4) * t480;
t641 = t478 * t450;
t685 = qJD(5) * t481;
t542 = -t641 - t685;
t383 = pkin(3) * t480 + qJ(4) * t477;
t451 = qJD(2) * t478;
t862 = t383 * t690 + t451;
t497 = t406 + t542 + t862;
t353 = pkin(4) * t762 + qJ(5) * t478;
t419 = qJ(4) * t757;
t315 = pkin(3) * t762 - t419;
t454 = t481 * qJ(2);
t380 = pkin(1) * t478 - t454;
t817 = pkin(7) * t478;
t628 = -t380 - t817;
t612 = t315 + t628;
t556 = t353 + t612;
t47 = -t906 + (-t323 + t556) * qJD(1) + t497 + t882;
t878 = t47 * t478;
t152 = -t285 * rSges(7,1) + t284 * rSges(7,2) + rSges(7,3) * t763;
t432 = pkin(8) * t763;
t321 = -pkin(5) * t759 + t432;
t433 = pkin(4) * t763;
t352 = -qJ(5) * t481 + t433;
t434 = pkin(3) * t763;
t310 = -qJ(4) * t759 + t434;
t613 = t310 + t858;
t557 = t352 + t613;
t448 = qJD(5) * t478;
t642 = t480 * t688;
t710 = -pkin(4) * t642 - t448;
t452 = qJD(2) * t481;
t720 = -t383 * t688 - t452;
t618 = t710 + t720;
t868 = (qJD(3) * t389 - t450) * t481;
t48 = t152 * t431 - t257 * t337 - t868 + (t321 + t557) * qJD(1) + t618;
t877 = t48 * t481;
t377 = rSges(6,1) * t477 - rSges(6,2) * t480;
t316 = t377 * t481;
t232 = -Icges(7,3) * t477 + t480 * t580;
t572 = -t239 * t476 + t247 * t479;
t578 = -t144 * t476 + t147 * t479;
t484 = t543 * (t231 * t481 + t909) - t337 * (-t231 * t478 - t578) - t431 * (-t232 - t572);
t874 = t484 * t477;
t385 = rSges(5,1) * t480 + rSges(5,3) * t477;
t867 = (qJD(3) * t385 - t450) * t481;
t554 = qJD(3) * t377 - t450;
t866 = t554 * t481;
t864 = t286 * t144 - t287 * t147;
t355 = qJD(1) * t380;
t863 = qJD(1) * t315 - t355;
t625 = -rSges(3,2) * t481 + t478 * rSges(3,3);
t861 = t387 + t625;
t708 = -t478 * rSges(5,2) - rSges(5,3) * t757;
t694 = qJD(1) * t478;
t650 = t477 * t694;
t658 = pkin(3) * t642 + qJ(4) * t901;
t686 = qJD(4) * t481;
t150 = pkin(3) * t650 + t480 * t686 - t658;
t139 = t481 * t150;
t653 = qJ(5) * t692 - t710;
t202 = pkin(4) * t650 - t653;
t317 = t383 * t481;
t728 = -t317 * t688 + t450;
t857 = t481 * t202 + t139 - t728;
t727 = t257 + t389;
t855 = qJD(1) * t727;
t621 = qJD(6) + t693;
t854 = t478 * t621 + t644;
t295 = (-Icges(7,2) * t479 - t786) * t477;
t494 = t543 * (-Icges(7,2) * t287 - t149 - t268) - t337 * (Icges(7,2) * t285 + t147 + t267) - t431 * (t247 + t295);
t302 = (-Icges(7,1) * t476 - t785) * t477;
t495 = t543 * (Icges(7,1) * t286 + t145 + t269) - t337 * (-Icges(7,1) * t284 + t144 - t787) - t431 * (t239 - t302);
t474 = t478 ^ 2;
t840 = m(5) / 0.2e1;
t839 = m(6) / 0.2e1;
t838 = m(7) / 0.2e1;
t837 = -pkin(1) - pkin(7);
t836 = pkin(3) + pkin(4);
t676 = qJD(3) * qJD(6);
t637 = t480 * t676;
t221 = qJD(1) * t543 + t478 * t637;
t835 = t221 / 0.2e1;
t222 = qJD(1) * t337 - t481 * t637;
t834 = t222 / 0.2e1;
t833 = t543 / 0.2e1;
t832 = -t543 / 0.2e1;
t831 = -t337 / 0.2e1;
t830 = t337 / 0.2e1;
t829 = -t431 / 0.2e1;
t828 = t431 / 0.2e1;
t827 = t478 / 0.2e1;
t826 = t480 / 0.2e1;
t824 = rSges(3,2) - pkin(1);
t823 = rSges(7,3) + pkin(8);
t822 = pkin(3) * t477;
t821 = pkin(4) * t477;
t820 = pkin(4) * t480;
t818 = pkin(5) * t480;
t645 = t477 * t690;
t131 = t431 * t760 + (t481 * t621 - t645) * t476;
t691 = qJD(3) * t477;
t132 = -t621 * t758 + (t431 * t476 + t479 * t691) * t478;
t74 = Icges(7,5) * t132 + Icges(7,6) * t131 + Icges(7,3) * t528;
t76 = Icges(7,4) * t132 + Icges(7,2) * t131 + Icges(7,6) * t528;
t78 = Icges(7,1) * t132 + Icges(7,4) * t131 + Icges(7,5) * t528;
t8 = (qJD(3) * t578 + t74) * t480 + (-qJD(3) * t141 - t476 * t76 + t479 * t78 + (-t144 * t479 - t147 * t476) * qJD(6)) * t477;
t816 = t8 * t337;
t555 = t481 * t431;
t129 = t476 * t854 - t479 * t555;
t130 = -t476 * t555 - t479 * t854;
t527 = -t642 + t650;
t73 = Icges(7,5) * t130 + Icges(7,6) * t129 + Icges(7,3) * t527;
t75 = Icges(7,4) * t130 + Icges(7,2) * t129 + Icges(7,6) * t527;
t77 = Icges(7,1) * t130 + Icges(7,4) * t129 + Icges(7,5) * t527;
t9 = (-qJD(3) * t909 + t73) * t480 + (qJD(3) * t143 - t476 * t75 + t479 * t77 + (-t145 * t479 + t149 * t476) * qJD(6)) * t477;
t815 = t9 * t543;
t812 = -pkin(5) - qJ(4);
t811 = rSges(6,1) * t480;
t809 = rSges(6,2) * t477;
t808 = rSges(3,3) * t481;
t261 = -rSges(7,3) * t477 + t480 * t601;
t309 = (-rSges(7,1) * t476 - rSges(7,2) * t479) * t477;
t178 = qJD(3) * t261 + qJD(6) * t309;
t647 = t480 * t692;
t655 = pkin(5) * t645 + pkin(8) * t528;
t194 = -pkin(5) * t647 + t655;
t390 = -pkin(8) * t477 + t818;
t354 = t390 * qJD(3);
t654 = pkin(4) * t649 + qJ(5) * t694 + t406;
t203 = t654 - t685;
t659 = pkin(3) * t528 + qJ(4) * t645;
t687 = qJD(4) * t478;
t151 = (-qJ(4) * t692 - t687) * t480 + t659;
t483 = qJD(1) ^ 2;
t678 = qJD(1) * qJD(2);
t707 = qJ(2) * t692 + t451;
t724 = qJD(1) * (-pkin(1) * t694 + t707) + t478 * t678;
t558 = -t483 * t817 + t724;
t677 = qJD(1) * qJD(3);
t639 = t383 * t677;
t530 = qJD(1) * t151 + t478 * t639 + t558;
t620 = t677 * t820;
t482 = qJD(3) ^ 2;
t671 = t482 * t821;
t496 = qJD(1) * t203 + t478 * t620 + t481 * t671 + t530;
t771 = qJ(4) * t480;
t376 = t771 - t822;
t449 = qJD(4) * t477;
t270 = qJD(3) * t376 + t449;
t622 = -t270 - t449;
t80 = t132 * rSges(7,1) + t131 * rSges(7,2) + rSges(7,3) * t528;
t10 = -t178 * t337 - t221 * t257 + t431 * t80 + (t194 + t542) * qJD(1) + (t389 * t694 - t152 * t684 + (-t354 + t622) * t481) * qJD(3) + t496;
t806 = t10 * t481;
t193 = -pkin(5) * t901 + pkin(8) * t527;
t418 = t477 * t687;
t440 = t481 * t678;
t614 = -t472 * t483 + t440;
t529 = qJD(3) * t418 + t270 * t690 + t481 * t639 + t614;
t498 = qJD(1) * t448 + t481 * t620 + t529;
t271 = qJD(1) * t387 - t452;
t750 = -t150 - t271;
t668 = -t202 + t750;
t602 = rSges(7,1) * t130 + rSges(7,2) * t129;
t79 = rSges(7,3) * t527 + t602;
t11 = -t482 * t433 - t543 * t178 + t222 * t257 - t431 * t79 + (t153 * t684 + t354 * t478) * qJD(3) + (-t193 + t668 + t868) * qJD(1) + t498;
t805 = t11 * t478;
t712 = rSges(6,1) * t645 + rSges(6,3) * t694;
t184 = -rSges(6,1) * t647 - rSges(6,2) * t528 + t712;
t384 = t809 + t811;
t350 = t384 * qJD(3);
t33 = (-t350 + t622) * t688 + (t478 * t554 + t184 - t685) * qJD(1) + t496;
t804 = t33 * t481;
t260 = -rSges(6,3) * t481 - t384 * t478;
t181 = qJD(1) * t260 - qJD(3) * t316;
t34 = (qJD(3) * t350 - t671) * t478 + (-t181 + t668 + t866) * qJD(1) + t498;
t801 = t34 * t478;
t657 = rSges(5,1) * t528 + rSges(5,3) * t645;
t182 = qJD(1) * t708 + t657;
t603 = rSges(5,1) * t477 - rSges(5,3) * t480;
t348 = t603 * qJD(3);
t42 = (t182 - t641) * qJD(1) + (t385 * t694 + (t348 + t622) * t481) * qJD(3) + t530;
t800 = t42 * t481;
t318 = t385 * t481;
t471 = t481 * rSges(5,2);
t179 = -qJD(3) * t318 + (t478 * t603 + t471) * qJD(1);
t43 = -t348 * t690 + (-t179 + t750 + t867) * qJD(1) + t529;
t799 = t43 * t478;
t53 = t141 * t480 + t477 * t578;
t798 = t53 * t221;
t797 = t54 * t222;
t656 = rSges(4,1) * t528 + rSges(4,2) * t647;
t183 = (-rSges(4,2) * t691 - rSges(4,3) * qJD(1)) * t478 + t656;
t386 = rSges(4,1) * t480 - rSges(4,2) * t477;
t328 = t386 * t690;
t604 = rSges(4,1) * t477 + rSges(4,2) * t480;
t349 = t604 * qJD(3);
t81 = t349 * t688 + (t183 + t328) * qJD(1) + t558;
t796 = t81 * t481;
t319 = t386 * t481;
t180 = -qJD(3) * t319 + (t478 * t604 + t470) * qJD(1);
t646 = t386 * t688;
t82 = -t349 * t690 + (-t180 - t271 + t646) * qJD(1) + t614;
t795 = t82 * t478;
t794 = -rSges(5,3) - qJ(4);
t467 = t478 * rSges(4,3);
t263 = t481 * t604 - t467;
t112 = t328 + t451 + (t263 + t628) * qJD(1);
t770 = t112 * t481;
t761 = t478 * t383;
t756 = qJD(3) * t450 + t150 * t688;
t749 = -t151 - t182;
t748 = -t151 - t203;
t747 = -t178 - t354;
t743 = t478 * t270 + t383 * t692;
t709 = rSges(5,1) * t763 + t471;
t258 = -rSges(5,3) * t759 + t709;
t726 = -t258 - t310;
t262 = rSges(5,1) * t762 + t708;
t725 = t262 + t315;
t266 = t481 * t315;
t723 = -t481 * t353 - t266;
t722 = -t310 - t352;
t721 = t315 + t353;
t719 = pkin(4) * t759 + t761;
t711 = t901 * pkin(4);
t706 = rSges(3,2) * t694 + rSges(3,3) * t692;
t705 = -t315 * t688 + t449;
t704 = t451 - t355;
t703 = t454 - t419;
t702 = t481 ^ 2 + t474;
t675 = -rSges(5,2) + t837;
t674 = rSges(6,2) - t836;
t673 = -rSges(4,3) + t837;
t672 = qJ(5) + t837;
t670 = pkin(4) * t691;
t669 = t202 * t688 + t756;
t667 = -t184 + t748;
t666 = -t194 + t748;
t416 = pkin(4) * t647;
t665 = t416 + t743;
t664 = -t260 + t722;
t466 = t478 * rSges(6,3);
t264 = t384 * t481 - t466;
t663 = -t264 + t721;
t662 = -t321 + t722;
t661 = -t323 + t721;
t330 = t383 * t694;
t660 = t330 + t711;
t652 = -t353 * t688 + t705;
t640 = t763 / 0.2e1;
t638 = t477 * t676;
t635 = t692 / 0.2e1;
t634 = -t691 / 0.2e1;
t633 = -t690 / 0.2e1;
t632 = t690 / 0.2e1;
t631 = -t688 / 0.2e1;
t629 = -t383 - t820;
t624 = t244 - t768;
t623 = -rSges(6,2) * qJD(3) - qJD(4);
t619 = -t152 + t662;
t617 = t452 + t658;
t616 = t434 + t858;
t615 = t702 * t820;
t608 = qJD(6) * t634;
t607 = qJD(1) * t761 - t477 * t686;
t606 = qJD(1) * t317 + t376 * t690 + t418;
t600 = -t10 * t478 - t11 * t481;
t599 = t478 * t50 + t481 * t49;
t598 = t478 * t49 - t481 * t50;
t51 = -t141 * t762 - t864;
t597 = t478 * t52 + t481 * t51;
t596 = t478 * t51 - t481 * t52;
t595 = t478 * t54 + t481 * t53;
t594 = t478 * t53 - t481 * t54;
t593 = t659 + t707;
t592 = t433 + t616;
t113 = qJD(1) * t903 - t452 - t646;
t579 = t112 * t478 - t113 * t481;
t576 = t152 * t481 + t153 * t478;
t288 = (-Icges(7,5) * t476 - Icges(7,6) * t479) * t477;
t155 = qJD(3) * t232 + qJD(6) * t288;
t240 = -Icges(7,6) * t477 + t480 * t584;
t162 = qJD(3) * t240 + qJD(6) * t295;
t248 = -Icges(7,5) * t477 + t480 * t588;
t169 = qJD(3) * t248 + qJD(6) * t302;
t20 = (qJD(3) * t572 + t155) * t480 + (-qJD(3) * t231 - t162 * t476 + t169 * t479 + (-t239 * t479 - t247 * t476) * qJD(6)) * t477;
t91 = t231 * t480 + t477 * t572;
t551 = t20 * t431 - t638 * t91;
t550 = t416 + t606;
t544 = -t233 - t568;
t535 = t607 + t711;
t534 = t385 * t690 - t641 + t862;
t114 = (-t259 * t478 - t263 * t481) * qJD(3);
t531 = t617 + t653;
t526 = t141 * t337 + t143 * t543 + t231 * t431;
t525 = (Icges(7,5) * t284 + Icges(7,6) * t285) * t337 - (-Icges(7,5) * t286 - Icges(7,6) * t287) * t543 + t288 * t431;
t515 = t593 + t654;
t493 = qJD(1) * t353 + t497 + t863;
t35 = t152 * t543 + t153 * t337 + (t323 * t481 + t478 * t662) * qJD(3) + t652;
t485 = t35 * t576 + (-t47 * t481 - t478 * t48) * t257;
t381 = rSges(3,2) * t478 + t808;
t325 = t377 * t690;
t322 = t389 * t481;
t320 = t389 * t478;
t314 = t386 * t478;
t313 = t385 * t478;
t281 = t702 * t689;
t280 = t645 - t647;
t211 = t257 * t481;
t210 = t257 * t478;
t209 = t247 * t481;
t208 = t247 * t478;
t207 = t239 * t481;
t206 = t239 * t478;
t201 = qJD(1) * t861 - t452;
t200 = t451 + (-t380 + t381) * qJD(1);
t192 = -rSges(7,1) * t286 - rSges(7,2) * t287;
t191 = rSges(7,1) * t284 + rSges(7,2) * t285;
t177 = t440 + (-qJD(1) * t625 - t271) * qJD(1);
t176 = qJD(1) * t706 + t724;
t90 = (-t262 * t481 + t478 * t726) * qJD(3) + t705;
t89 = -t867 + (t258 + t613) * qJD(1) + t720;
t88 = (t262 + t612) * qJD(1) + t534;
t72 = -t866 + (t260 + t557) * qJD(1) + t618;
t71 = t325 + (-t264 + t556) * qJD(1) + t497;
t70 = (t264 * t481 + t478 * t664) * qJD(3) + t652;
t18 = (t179 * t481 + t749 * t478 + (t478 * t725 + t481 * t726) * qJD(1)) * qJD(3) + t756;
t17 = t131 * t239 + t132 * t247 + t155 * t763 + t162 * t284 - t169 * t285 + t231 * t528;
t16 = t129 * t239 + t130 * t247 - t155 * t762 - t162 * t286 + t169 * t287 + t231 * t527;
t15 = (t181 * t481 + t667 * t478 + (t478 * t663 + t481 * t664) * qJD(1)) * qJD(3) + t669;
t14 = t337 * t53 + t431 * t91 - t54 * t543;
t13 = t337 * t51 - t881;
t7 = t131 * t145 - t132 * t149 - t143 * t528 + t284 * t75 - t285 * t77 + t73 * t763;
t6 = t131 * t144 + t132 * t147 + t141 * t528 + t284 * t76 - t285 * t78 + t74 * t763;
t5 = t129 * t145 - t130 * t149 - t143 * t527 - t286 * t75 + t287 * t77 - t73 * t762;
t4 = t129 * t144 + t130 * t147 + t141 * t527 - t286 * t76 + t287 * t78 - t74 * t762;
t3 = -t152 * t222 + t153 * t221 + t543 * t80 + t337 * t79 + (t193 * t481 + t666 * t478 + (t478 * t661 + t481 * t662) * qJD(1)) * qJD(3) + t669;
t2 = t17 * t431 + t221 * t49 + t222 * t50 + t337 * t6 - t543 * t7 - t638 * t86;
t1 = t16 * t431 + t221 * t51 + t222 * t52 + t337 * t4 + t499 * t638 - t5 * t543;
t19 = [(-t936 * qJD(3) - t941 * t477 + t940 * t480) * qJD(1) + ((t50 + (t141 * t481 + t143 * t478) * t477 + t864 + t910) * t337 + t13 + t881) * t831 + (-(t325 - t71 + (-t264 - t817) * qJD(1) + t493) * t72 + t34 * (t466 + t703) + t71 * t531 + t33 * t592 + t72 * (t515 + t712) + (t34 * t672 - t33 * t809 + (t33 * (-rSges(6,1) - qJ(4)) + t72 * t623) * t480 + (t72 * t837 + (t477 * t674 - qJ(2) + t811) * t71) * qJD(1)) * t478 + (t33 * (-rSges(6,3) - qJ(5)) - t72 * qJD(5) + (-t34 * rSges(6,1) + t623 * t71) * t480 + (t71 * rSges(6,1) * qJD(3) - t34 * t674) * t477 + (t71 * (rSges(6,3) + t837) + t72 * (-t384 - t771)) * qJD(1)) * t481) * m(6) - t499 * t834 + (t177 * (t478 * t824 + t454 + t808) + t200 * t452 + t176 * t861 + t201 * (t706 + t707) + (t200 * t824 * t481 + (t200 * (-rSges(3,3) - qJ(2)) - t201 * pkin(1)) * t478) * qJD(1) - (qJD(1) * t381 - t200 + t704) * t201) * m(3) + (((t481 * t544 - t536 + t746 + t893 + t94 - t944) * t481 + (t478 * t544 + (t624 + t238 - t902) * t481 - t932 + t908 + t955) * t478) * qJD(3) + t916) * t633 + ((-t745 + t97 + t769) * t688 + (-t96 + (t234 - t570) * t478) * t690 + (t238 * t474 + (t478 * t624 + t746 - t92) * t478 + ((t902 - t957) * t481 + t904 + t908) * t481) * qJD(3) + t899 - t913) * t631 + (-(-t47 + (-t323 - t817) * qJD(1) + t493 + t882) * t48 + t48 * t906 + t11 * (-t153 + t436 + t703) + t47 * (t531 - t602) + t10 * (t432 + t592 + t152) + t48 * (t515 + t80 + t655) + (t11 * t672 + (-t48 * qJD(4) + t10 * t812) * t480) * t478 + (t11 * (-t818 + t821 + t822) - t10 * qJ(5) - t48 * qJD(5) + (-t450 + (t480 * t823 + t819) * qJD(3)) * t47) * t481 + ((t48 * t480 * t812 + t47 * t837) * t481 + (t48 * t837 + (-qJ(2) + t818 + (-t823 - t836) * t477) * t47) * t478) * qJD(1)) * m(7) + t12 * t833 + t86 * t835 + t17 * t830 + (t896 + t897 + t900) * t632 + t797 / 0.2e1 + t798 / 0.2e1 + (t82 * (-t467 + t628) + t112 * t452 + t81 * t903 + t113 * (-rSges(4,2) * t645 + t656 + t707) + (qJD(3) * t112 * t386 + t604 * t82) * t481 + (t673 * t770 + (t112 * (-qJ(2) - t604) + t113 * t673) * t478) * qJD(1) - (-t112 + t328 + (t263 - t817) * qJD(1) + t704) * t113) * m(4) + (t16 + t12) * t832 + t551 - t815 / 0.2e1 + t816 / 0.2e1 + (-(-t88 + (t262 - t817) * qJD(1) + t534 + t863) * t89 + t43 * (t703 + t708) + t88 * t617 + t42 * (t616 + t709) + t89 * (t593 + t657) + (t43 * t837 + (-t89 * qJD(4) + t42 * t794) * t480) * t478 + (t88 * (rSges(5,1) * qJD(3) - qJD(4)) * t480 + (t43 * (rSges(5,1) + pkin(3)) + t88 * rSges(5,3) * qJD(3)) * t477) * t481 + ((t480 * t794 * t89 + t675 * t88) * t481 + (t88 * (-qJ(2) - t603 - t822) + t89 * t675) * t478) * qJD(1)) * m(5) + (t895 + t898 + (t138 + t889) * qJD(1)) * t688 / 0.2e1 - (t935 * t481 + (-t888 + t931) * t478) * t677 / 0.2e1; 0.2e1 * (-t806 / 0.2e1 + t805 / 0.2e1) * m(7) + 0.2e1 * (-t804 / 0.2e1 + t801 / 0.2e1) * m(6) + 0.2e1 * (-t800 / 0.2e1 + t799 / 0.2e1) * m(5) + 0.2e1 * (t795 / 0.2e1 - t796 / 0.2e1) * m(4) + 0.2e1 * (-t176 * t481 / 0.2e1 + t177 * t827) * m(3); t14 * t684 / 0.2e1 + (-qJD(1) * t596 + t4 * t481 + t478 * t5) * t832 + ((-t206 * t286 + t208 * t287) * t337 - (t207 * t286 - t209 * t287) * t543 + (-t240 * t286 + t248 * t287) * t431 + (t477 * t499 + t51 * t759) * qJD(6) + ((-qJD(6) * t52 - t526) * t480 - t874) * t481) * t833 + t597 * t834 + t599 * t835 + (-qJD(1) * t594 + t478 * t9 + t481 * t8) * t828 + (((-t206 * t476 + t208 * t479 - t141) * t337 - (t207 * t476 - t209 * t479 + t143) * t543 + (-t240 * t476 + t248 * t479 - t231) * t431 - t91 * qJD(6)) * t477 + (qJD(6) * t594 + t484) * t480) * t829 + (-qJD(1) * t598 + t478 * t7 + t481 * t6) * t830 + ((t206 * t284 - t208 * t285) * t337 - (-t207 * t284 + t209 * t285) * t543 + (t240 * t284 - t248 * t285) * t431 + (-t477 * t86 - t50 * t757) * qJD(6) + ((qJD(6) * t49 + t526) * t480 + t874) * t478) * t831 + t595 * t608 - t478 * t12 * t683 / 0.2e1 + t13 * t682 * t826 + (t11 * t719 + t47 * t665 + t48 * t660 + t3 * t723 + (t10 * (t629 - t727) + t48 * (-t270 + t747) + t3 * (t153 + t323) + t47 * t855) * t481 + (t11 * t727 + t47 * (-t670 - t747) + t3 * t619 + t48 * t855) * t478 - t47 * (qJD(1) * t322 + t211 * t431 - t261 * t543 + t550) - t48 * (qJD(1) * t320 + t210 * t431 - t261 * t337 + t535) - ((-t152 * t48 + t153 * t47) * t477 + t485 * t480) * qJD(6) - ((-t376 - t390) * t877 + (t390 - t821) * t878) * qJD(3) + ((qJD(1) * t619 + t193 + t79) * t481 + (-t80 + t666 + (-t153 + t661) * qJD(1)) * t478 - t210 * t543 + t211 * t337 - (-t615 - t322 * t481 + (-t761 - t320) * t478) * qJD(3) + t857) * t35) * m(7) + (t34 * t719 + t15 * t723 + (t33 * (-t377 + t629) + t15 * t264) * t481 + (t15 * t664 + t34 * t377) * t478 + (t660 - t535 + (-t270 - t350 - (-t376 - t384) * qJD(3)) * t481) * t72 + (t665 - t550 + (t350 - t670 - (t384 - t821) * qJD(3)) * t478) * t71 + ((qJD(1) * t664 + t181) * t481 + (qJD(1) * t663 + t667) * t478 - (-t615 - t316 * t481 + (-t377 * t478 - t761) * t478) * qJD(3) + t857) * t70) * m(6) + (t43 * t761 + t88 * t743 + t89 * t330 - t18 * t266 + t90 * t139 + (t43 * t385 - t88 * t348 + t18 * t726 + t90 * t749 + (t89 * t385 + t725 * t90) * qJD(1)) * t478 + (t42 * (-t383 - t385) + t89 * (-t270 + t348) - t18 * t262 + t90 * t179 + (t88 * t385 + t726 * t90) * qJD(1)) * t481 - t88 * (qJD(1) * t318 + t606) - t89 * (qJD(1) * t313 + t607) - t90 * t728 - ((t89 * (-t376 + t603) - t90 * t318) * t481 + (-t88 * t603 + t90 * (-t761 - t313)) * t478) * qJD(3)) * m(5) + (0.2e1 * t114 * (t180 * t481 - t183 * t478 + (-t259 * t481 + t263 * t478) * qJD(1)) - t579 * t349 + (t795 - t796 + (t113 * t478 + t770) * qJD(1)) * t386 - (t112 * t319 + t113 * t314) * qJD(1) - (t114 * (-t314 * t478 - t319 * t481) - t579 * t604) * qJD(3)) * m(4) - ((t924 * t477 - t480 * t934) * qJD(3) + (t920 * t477 - t921 * t480) * qJD(1)) * qJD(1) / 0.2e1 + (t898 * t481 + t897 * t478 + (t888 * t478 + t889 * t481) * qJD(1)) * qJD(1) / 0.2e1 + ((t690 * t886 + t887) * t478 + ((t478 * t892 + t911) * qJD(3) + t917) * t481) * t633 + ((t688 * t892 + t887) * t481 + ((t481 * t886 - t911) * qJD(3) - t917) * t478) * t631 + (t1 + t896 * qJD(1) + ((t893 * qJD(1) + t883 * t481) * t481 + (t890 * t478 - t932 * qJD(1) + (-t884 + t891) * t481) * t478) * t879) * t827 + (t2 + t895 * qJD(1) + ((t954 * qJD(1) + t891 * t481) * t481 + (t884 * t478 - t933 * qJD(1) + (-t883 + t890) * t481) * t478) * t879) * t481 / 0.2e1 - (t12 + t900 + t914) * t694 / 0.2e1 + (t13 + t899 + t915) * t635; -m(5) * (t280 * t88 + t281 * t90 - t89 * t901) - m(6) * (t280 * t71 + t281 * t70 - t72 * t901) - m(7) * (t280 * t47 + t281 * t35 - t48 * t901) + 0.2e1 * ((-t688 * t89 + t690 * t88 + t18) * t840 + (-t688 * t72 + t690 * t71 + t15) * t839 + (t47 * t690 - t48 * t688 + t3) * t838) * t477 + 0.2e1 * ((qJD(3) * t90 - t692 * t88 - t694 * t89 - t799 + t800) * t840 + (qJD(3) * t70 - t692 * t71 - t694 * t72 - t801 + t804) * t839 + (qJD(3) * t35 - t47 * t692 - t48 * t694 - t805 + t806) * t838) * t480; m(6) * (-t33 * t478 - t34 * t481) + m(7) * t600; t2 * t640 + (t477 * t598 + t480 * t86) * t835 + ((qJD(3) * t598 + t17) * t480 + (qJD(1) * t599 - qJD(3) * t86 + t478 * t6 - t481 * t7) * t477) * t830 - t1 * t762 / 0.2e1 + (t477 * t596 - t480 * t499) * t834 + ((qJD(3) * t596 + t16) * t480 + (qJD(1) * t597 + qJD(3) * t499 + t4 * t478 - t481 * t5) * t477) * t832 + t14 * t634 + (t551 + t797 + t798 - t815 + t816) * t826 + (t477 * t594 + t480 * t91) * t608 + ((qJD(3) * t594 + t20) * t480 + (qJD(1) * t595 - qJD(3) * t91 + t478 * t8 - t481 * t9) * t477) * t828 + (-t284 * t494 - t285 * t495 + t525 * t763) * t831 + (t286 * t494 + t287 * t495 - t525 * t762) * t833 + (t525 * t480 + (t476 * t494 + t479 * t495) * t477) * t829 + (qJD(1) * t640 + t480 * t631) * t13 + (t477 * t635 + t480 * t632) * t12 + ((qJD(3) * t485 + t10 * t152 - t11 * t153 - t47 * t79 + t48 * t80) * t480 + (t47 * (qJD(3) * t153 - t178 * t481) + t48 * (-qJD(3) * t152 - t178 * t478) + t3 * t576 + t35 * (-t152 * t694 + t153 * t692 + t478 * t79 + t481 * t80) + ((-t877 + t878) * qJD(1) + t600) * t257) * t477 - t47 * (-t192 * t431 - t309 * t543) - t48 * (t191 * t431 - t309 * t337) - t35 * (t191 * t543 + t192 * t337)) * m(7);];
tauc  = t19(:);
