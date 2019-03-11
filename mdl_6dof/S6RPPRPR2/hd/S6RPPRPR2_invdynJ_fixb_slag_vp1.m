% Calculate vector of inverse dynamics joint torques for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPRPR2_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR2_invdynJ_fixb_slag_vp1: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR2_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_invdynJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_invdynJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR2_invdynJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR2_invdynJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:41:44
% EndTime: 2019-03-09 01:42:35
% DurationCPUTime: 51.29s
% Computational Cost: add. (28323->1009), mult. (25828->1259), div. (0->0), fcn. (22861->10), ass. (0->485)
t850 = Icges(6,4) - Icges(5,5);
t849 = Icges(6,5) - Icges(5,6);
t848 = Icges(6,1) + Icges(5,3);
t429 = pkin(10) + qJ(4);
t424 = sin(t429);
t426 = cos(t429);
t816 = t849 * t424 - t426 * t850;
t430 = qJ(1) + pkin(9);
t427 = cos(t430);
t847 = t848 * t427;
t425 = sin(t430);
t670 = t425 * t426;
t672 = t424 * t425;
t817 = t670 * t850 - t849 * t672 + t847;
t824 = t425 * t848 + t816 * t427;
t687 = Icges(5,6) * t427;
t203 = Icges(5,4) * t670 - Icges(5,2) * t672 - t687;
t357 = Icges(6,6) * t672;
t695 = Icges(6,4) * t427;
t210 = Icges(6,2) * t670 - t357 + t695;
t846 = t203 * t424 - t210 * t426;
t697 = Icges(5,4) * t424;
t317 = Icges(5,1) * t426 - t697;
t206 = Icges(5,5) * t425 + t317 * t427;
t685 = Icges(6,6) * t426;
t513 = -Icges(6,3) * t424 + t685;
t207 = Icges(6,5) * t425 - t427 * t513;
t845 = -t206 * t670 - t207 * t672;
t310 = Icges(5,5) * t424 + Icges(5,6) * t426;
t518 = Icges(6,4) * t424 + Icges(6,5) * t426;
t832 = t310 - t518;
t314 = Icges(5,2) * t426 + t697;
t686 = Icges(6,6) * t424;
t512 = Icges(6,3) * t426 + t686;
t844 = -t314 - t512;
t408 = Icges(5,4) * t426;
t316 = Icges(5,1) * t424 + t408;
t514 = Icges(6,2) * t424 + t685;
t843 = t316 + t514;
t515 = Icges(6,2) * t426 - t686;
t842 = t317 + t515;
t362 = Icges(5,4) * t672;
t691 = Icges(5,5) * t427;
t205 = Icges(5,1) * t670 - t362 - t691;
t690 = Icges(6,5) * t427;
t208 = Icges(6,6) * t670 - Icges(6,3) * t672 + t690;
t812 = -t205 * t426 + t208 * t424 + t846;
t520 = -Icges(5,2) * t424 + t408;
t841 = t513 + t520;
t497 = t314 * t424 - t316 * t426;
t828 = -t424 * t512 + t426 * t514 - t497;
t840 = t427 * t824 + t845;
t667 = t426 * t427;
t671 = t424 * t427;
t788 = t206 * t667 + t207 * t671 + t425 * t824;
t839 = -t205 * t667 + t208 * t671 + t425 * t817;
t783 = -t425 * t812 + t427 * t817;
t204 = Icges(5,6) * t425 + t427 * t520;
t358 = Icges(6,6) * t671;
t696 = Icges(6,4) * t425;
t209 = -Icges(6,2) * t667 + t358 + t696;
t782 = -t204 * t672 - t209 * t670 - t840;
t781 = -t203 * t671 + t210 * t667 - t839;
t780 = -t204 * t671 - t209 * t667 + t788;
t838 = t203 + t208;
t837 = t204 - t207;
t836 = t205 + t210;
t835 = -t206 + t209;
t834 = t841 * qJD(4);
t833 = t842 * qJD(4);
t830 = t844 * qJD(4);
t829 = t843 * qJD(4);
t827 = -t424 * t843 + t426 * t844;
t764 = t832 * t425;
t673 = t310 * t427;
t109 = -t425 * t497 - t673;
t253 = t518 * t427;
t112 = t512 * t672 - t514 * t670 - t253;
t825 = t109 - t112;
t795 = t427 * t828 + t764;
t433 = -pkin(7) - qJ(3);
t391 = t427 * t433;
t432 = cos(pkin(10));
t420 = pkin(3) * t432 + pkin(2);
t626 = -t425 * t420 - t391;
t435 = sin(qJ(1));
t728 = pkin(1) * t435;
t548 = t626 - t728;
t823 = t204 * t424 + t209 * t426;
t779 = t424 * t836 + t426 * t838;
t822 = t424 * t835 - t426 * t837;
t821 = t830 * t427 + (-t425 * t841 + t687 - t690) * qJD(1);
t820 = qJD(1) * t837 + t425 * t830;
t819 = -t829 * t427 + (-t425 * t842 + t691 - t695) * qJD(1);
t818 = t829 * t425 + (-t427 * t515 - t206 + t696) * qJD(1);
t815 = qJD(1) * t832 + qJD(4) * t827 - t424 * t834 + t426 * t833;
t814 = t832 * qJD(4);
t597 = rSges(5,1) * t670;
t813 = -t597 + t548;
t811 = t206 * t426 + t207 * t424 - t823;
t810 = t425 * t780 - t427 * t781;
t809 = t782 * t425 - t427 * t783;
t808 = t828 * qJD(1) - qJD(4) * t816;
t437 = cos(qJ(1));
t428 = t437 * pkin(1);
t807 = t795 * qJD(1);
t326 = rSges(3,1) * t425 + rSges(3,2) * t427;
t291 = -t326 - t728;
t806 = t825 * qJD(1);
t805 = t824 * qJD(1);
t434 = sin(qJ(6));
t665 = t427 * t434;
t436 = cos(qJ(6));
t668 = t425 * t436;
t279 = t424 * t668 + t665;
t664 = t427 * t436;
t669 = t425 * t434;
t280 = -t424 * t669 + t664;
t638 = t280 * rSges(7,1) - t279 * rSges(7,2);
t153 = rSges(7,3) * t670 - t638;
t541 = rSges(7,1) * t434 + rSges(7,2) * t436;
t767 = t426 * t541;
t238 = rSges(7,3) * t424 - t767;
t607 = qJD(6) * t426;
t611 = qJD(4) * t427;
t287 = -t425 * t607 + t611;
t608 = qJD(6) * t424;
t385 = qJD(1) + t608;
t322 = pkin(4) * t424 - qJ(5) * t426;
t567 = -pkin(8) * t424 - t322;
t609 = qJD(5) * t427;
t351 = t424 * t609;
t400 = qJD(3) * t425;
t627 = t351 + t400;
t804 = t153 * t385 + t287 * t238 - t567 * t611 - t627;
t803 = t427 ^ 2;
t439 = qJD(1) ^ 2;
t802 = t439 * t428;
t801 = qJD(4) * t809 + t806;
t800 = qJD(4) * t810 + t807;
t799 = qJD(4) * t812 + t424 * t818 - t426 * t820;
t798 = qJD(4) * t811 + t424 * t819 + t426 * t821;
t797 = -t425 * t808 + t427 * t815;
t796 = t425 * t815 + t427 * t808;
t794 = t817 * qJD(1) + qJD(4) * t779 + t820 * t424 + t818 * t426;
t793 = qJD(4) * t822 - t424 * t821 + t426 * t819 + t805;
t792 = -t425 * t837 + t427 * t838;
t791 = t817 + t823;
t790 = t841 + t843;
t789 = t842 + t844;
t787 = (t357 + t362 + (Icges(5,2) + Icges(6,3)) * t670 - t836) * t427 + (-Icges(6,3) * t667 - t314 * t427 - t358 - t835) * t425;
t786 = -t814 * t427 + (-t425 * t816 - t811 + t847) * qJD(1);
t785 = qJD(1) * t812 - t425 * t814 + t805;
t413 = t426 * rSges(7,3);
t237 = t424 * t541 + t413;
t290 = pkin(5) * t427 - pkin(8) * t670;
t784 = t290 + t548 + t638;
t613 = qJD(4) * t425;
t286 = t427 * t607 + t613;
t277 = t424 * t664 - t669;
t278 = t424 * t665 + t668;
t142 = Icges(7,5) * t278 + Icges(7,6) * t277 + Icges(7,3) * t667;
t694 = Icges(7,4) * t278;
t145 = Icges(7,2) * t277 + Icges(7,6) * t667 + t694;
t258 = Icges(7,4) * t277;
t148 = Icges(7,1) * t278 + Icges(7,5) * t667 + t258;
t41 = t142 * t667 + t277 * t145 + t278 * t148;
t144 = -Icges(7,5) * t280 + Icges(7,6) * t279 + Icges(7,3) * t670;
t260 = Icges(7,4) * t280;
t147 = Icges(7,2) * t279 + Icges(7,6) * t670 - t260;
t259 = Icges(7,4) * t279;
t149 = Icges(7,1) * t280 - Icges(7,5) * t670 - t259;
t42 = t144 * t667 + t277 * t147 - t149 * t278;
t516 = Icges(7,5) * t434 + Icges(7,6) * t436;
t459 = -Icges(7,3) * t424 + t426 * t516;
t693 = Icges(7,4) * t434;
t517 = Icges(7,2) * t436 + t693;
t460 = -Icges(7,6) * t424 + t426 * t517;
t692 = Icges(7,4) * t436;
t521 = Icges(7,1) * t434 + t692;
t461 = -Icges(7,5) * t424 + t426 * t521;
t70 = -t277 * t460 - t278 * t461 - t459 * t667;
t12 = t286 * t41 - t287 * t42 + t70 * t385;
t43 = t142 * t670 + t279 * t145 - t280 * t148;
t44 = t144 * t670 + t147 * t279 + t149 * t280;
t71 = -t279 * t460 + t280 * t461 - t459 * t670;
t13 = t286 * t43 - t287 * t44 + t385 * t71;
t509 = t147 * t436 - t149 * t434;
t55 = t144 * t424 - t426 * t509;
t773 = -t424 * t787 + t426 * t792;
t772 = (-t424 * t790 + t426 * t789) * qJD(1);
t771 = t794 * t803 + (t786 * t425 + (-t785 + t793) * t427) * t425;
t770 = t785 * t803 + (t793 * t425 + (-t786 + t794) * t427) * t425;
t602 = qJD(4) * qJD(5);
t769 = qJDD(5) * t424 + t426 * t602;
t582 = t424 * t611;
t616 = qJD(1) * t425;
t585 = t426 * t616;
t768 = t582 + t585;
t766 = t816 * qJD(1);
t765 = t673 - t253;
t417 = t425 * pkin(5);
t404 = t427 * qJ(3);
t325 = pkin(2) * t425 - t404;
t197 = t325 + t626;
t305 = qJD(1) * t325;
t762 = qJD(1) * t197 - t305;
t431 = sin(pkin(10));
t715 = rSges(4,2) * t431;
t718 = rSges(4,1) * t432;
t229 = t425 * rSges(4,3) + (-t715 + t718) * t427;
t403 = t425 * qJ(3);
t330 = t427 * pkin(2) + t403;
t568 = t330 + t428;
t173 = t229 + t568;
t402 = t424 * qJ(5);
t759 = t426 * pkin(4) + t402;
t761 = -pkin(8) * t426 - t759;
t331 = t427 * rSges(3,1) - rSges(3,2) * t425;
t292 = t331 + t428;
t760 = -rSges(5,2) * t672 - t427 * rSges(5,3);
t409 = t424 * rSges(6,3);
t714 = rSges(6,2) * t426;
t540 = t409 - t714;
t289 = pkin(8) * t667 + t417;
t729 = -rSges(7,3) - pkin(4);
t599 = -pkin(8) + t729;
t758 = t599 * t426 - t402 - t420;
t266 = t759 * t425;
t569 = -t325 - t728;
t551 = t197 + t569;
t494 = -t266 + t551;
t484 = t290 + t494;
t39 = qJD(1) * t484 - t804;
t151 = t278 * rSges(7,1) + t277 * rSges(7,2) + rSges(7,3) * t667;
t583 = t424 * t613;
t342 = pkin(8) * t583;
t270 = pkin(4) * t667 + qJ(5) * t671;
t559 = t427 * t420 - t425 * t433;
t198 = t559 - t330;
t550 = t198 + t568;
t493 = t270 + t550;
t399 = qJD(5) * t424;
t579 = t425 * t399;
t401 = qJD(3) * t427;
t637 = -t322 * t613 - t401;
t40 = t579 + t151 * t385 - t238 * t286 - t342 + (t289 + t493) * qJD(1) + t637;
t757 = t39 * t427 + t40 * t425;
t754 = -qJD(1) * t266 + t762;
t284 = (Icges(7,2) * t434 - t692) * t426;
t452 = t286 * (-Icges(7,2) * t278 + t148 + t258) - t287 * (Icges(7,2) * t280 - t149 + t259) + t385 * (-t461 + t284);
t285 = (-Icges(7,1) * t436 + t693) * t426;
t453 = t286 * (-Icges(7,1) * t277 + t145 + t694) - t287 * (-Icges(7,1) * t279 + t147 - t260) + t385 * (-t460 - t285);
t745 = m(6) / 0.2e1;
t744 = m(7) / 0.2e1;
t743 = -m(6) - m(7);
t603 = qJD(1) * qJD(4);
t302 = qJDD(4) * t425 + t427 * t603;
t600 = qJDD(6) * t426;
t156 = -qJD(6) * t768 + t427 * t600 + t302;
t742 = t156 / 0.2e1;
t303 = -qJDD(4) * t427 + t425 * t603;
t615 = qJD(1) * t427;
t584 = t426 * t615;
t474 = -t583 + t584;
t157 = qJD(6) * t474 + t425 * t600 + t303;
t741 = t157 / 0.2e1;
t740 = -t286 / 0.2e1;
t739 = t286 / 0.2e1;
t738 = -t287 / 0.2e1;
t737 = t287 / 0.2e1;
t293 = qJD(4) * t607 + qJDD(6) * t424 + qJDD(1);
t736 = t293 / 0.2e1;
t735 = t302 / 0.2e1;
t734 = t303 / 0.2e1;
t733 = -t385 / 0.2e1;
t732 = t385 / 0.2e1;
t731 = t425 / 0.2e1;
t730 = -t427 / 0.2e1;
t726 = g(1) * t425;
t725 = g(2) * t425;
t510 = t145 * t436 + t148 * t434;
t612 = qJD(4) * t426;
t462 = -t385 * t434 + t436 * t612;
t554 = qJD(1) * t424 + qJD(6);
t492 = t425 * t554;
t124 = t427 * t462 - t436 * t492;
t463 = t385 * t436 + t434 * t612;
t125 = t427 * t463 - t434 * t492;
t61 = Icges(7,5) * t125 + Icges(7,6) * t124 - Icges(7,3) * t768;
t63 = Icges(7,4) * t125 + Icges(7,2) * t124 - Icges(7,6) * t768;
t65 = Icges(7,1) * t125 + Icges(7,4) * t124 - Icges(7,5) * t768;
t8 = (qJD(4) * t510 + t61) * t424 + (qJD(4) * t142 - t434 * t65 - t436 * t63 + (t145 * t434 - t148 * t436) * qJD(6)) * t426;
t724 = t8 * t286;
t491 = t427 * t554;
t122 = t425 * t462 + t436 * t491;
t123 = t425 * t463 + t434 * t491;
t60 = Icges(7,5) * t123 + Icges(7,6) * t122 + Icges(7,3) * t474;
t62 = Icges(7,4) * t123 + Icges(7,2) * t122 + Icges(7,6) * t474;
t64 = Icges(7,1) * t123 + Icges(7,4) * t122 + Icges(7,5) * t474;
t9 = (qJD(4) * t509 + t60) * t424 + (qJD(4) * t144 - t434 * t64 - t436 * t62 + (t147 * t434 + t149 * t436) * qJD(6)) * t426;
t723 = t9 * t287;
t720 = pkin(2) - t420;
t230 = Icges(7,3) * t426 + t424 * t516;
t283 = (-Icges(7,5) * t436 + Icges(7,6) * t434) * t426;
t160 = qJD(4) * t230 + qJD(6) * t283;
t232 = Icges(7,6) * t426 + t424 * t517;
t161 = qJD(4) * t232 + qJD(6) * t284;
t234 = Icges(7,5) * t426 + t424 * t521;
t162 = qJD(4) * t234 + qJD(6) * t285;
t501 = -t434 * t461 - t436 * t460;
t29 = (qJD(4) * t501 + t160) * t424 + (-qJD(4) * t459 - t161 * t436 - t162 * t434 + (-t434 * t460 + t436 * t461) * qJD(6)) * t426;
t88 = -t424 * t459 - t426 * t501;
t719 = t29 * t385 + t88 * t293;
t288 = (-rSges(7,1) * t436 + rSges(7,2) * t434) * t426;
t163 = qJD(4) * t237 + qJD(6) * t288;
t396 = pkin(5) * t615;
t200 = -pkin(8) * t768 + t396;
t610 = qJD(5) * t426;
t261 = qJD(4) * t759 - t610;
t581 = t426 * t611;
t336 = qJ(5) * t581;
t586 = t424 * t616;
t120 = -pkin(4) * t768 - qJ(5) * t586 + t336 + t351;
t392 = qJ(3) * t615;
t552 = qJDD(1) * t428 - t439 * t728;
t604 = qJD(1) * qJD(3);
t622 = t392 + t400;
t456 = -qJDD(3) * t427 + qJD(1) * (-pkin(2) * t616 + t622) + qJDD(1) * t330 + t425 * t604 + t552;
t454 = qJD(1) * (-t392 + (t425 * t720 - t391) * qJD(1)) + qJDD(1) * t198 + t456;
t444 = qJDD(1) * t270 + t454 + t769 * t425 + (t120 + t351) * qJD(1);
t666 = t426 * qJD(4) ^ 2;
t662 = t125 * rSges(7,1) + t124 * rSges(7,2);
t67 = -rSges(7,3) * t768 + t662;
t10 = qJD(1) * t200 + (-t302 * t424 - t425 * t666) * pkin(8) + t444 + qJDD(1) * t289 + t293 * t151 - t286 * t163 - t302 * t322 - t261 * t613 + t385 * t67 - t156 * t238;
t712 = t10 * t427;
t199 = qJD(1) * t289 - t342;
t496 = qJDD(3) * t425 + t427 * t604 - t802;
t465 = t303 * t322 + t427 * t769 + t496;
t244 = t424 * t615 + t425 * t612;
t343 = pkin(4) * t583;
t121 = pkin(4) * t584 + qJ(5) * t244 - t343 + t579;
t262 = qJD(1) * t330 - t401;
t378 = t433 * t616;
t651 = t378 - (-t427 * t720 - t403) * qJD(1) - t262;
t482 = -t121 - t579 + t651;
t543 = rSges(7,1) * t123 + rSges(7,2) * t122;
t66 = rSges(7,3) * t474 + t543;
t11 = -t261 * t611 - t293 * t153 + t157 * t238 - t287 * t163 - t385 * t66 + (t303 * t424 - t427 * t666) * pkin(8) + t484 * qJDD(1) + (-t199 + t482) * qJD(1) + t465;
t711 = t11 * t425;
t324 = rSges(5,1) * t424 + rSges(5,2) * t426;
t269 = t324 * t427;
t410 = t425 * rSges(5,3);
t222 = rSges(5,1) * t667 - rSges(5,2) * t671 + t410;
t87 = -t324 * t613 - t401 + (t222 + t550) * qJD(1);
t710 = t269 * t87;
t704 = t424 * rSges(6,2);
t412 = t425 * rSges(6,1);
t221 = t597 + t760;
t495 = -t221 + t551;
t545 = -t324 * t611 + t400;
t86 = qJD(1) * t495 + t545;
t703 = t425 * t86;
t54 = t142 * t424 - t426 * t510;
t702 = t54 * t156;
t701 = t55 * t157;
t488 = t266 * t613 + t270 * t611 + qJD(2) - t610;
t36 = t151 * t287 + t153 * t286 + (t289 * t427 - t290 * t425) * qJD(4) + t488;
t681 = qJD(4) * t36;
t657 = t151 + t289;
t656 = t153 - t290;
t223 = -rSges(6,2) * t667 + rSges(6,3) * t671 + t412;
t645 = -t223 - t270;
t644 = t425 * t266 + t427 * t270;
t355 = qJ(5) * t667;
t267 = -pkin(4) * t671 + t355;
t641 = qJD(1) * t267 + t425 * t610;
t640 = -t540 * qJD(4) - t261;
t639 = -t270 - t289;
t632 = t767 * t425;
t631 = t767 * t427;
t539 = rSges(6,3) * t426 + t704;
t630 = -t322 + t539;
t629 = -t759 - t540;
t628 = rSges(5,2) * t586 + rSges(5,3) * t615;
t264 = rSges(6,2) * t672 + rSges(6,3) * t670;
t268 = rSges(6,2) * t671 + rSges(6,3) * t667;
t379 = t425 * t715;
t625 = rSges(4,3) * t615 + qJD(1) * t379;
t624 = t378 + t401;
t623 = t427 * rSges(4,3) + t379;
t621 = t425 ^ 2 + t803;
t614 = qJD(4) * t424;
t598 = t425 * t718;
t592 = -m(4) - m(5) + t743;
t591 = t40 * t615;
t590 = t427 * t120 + t425 * t121 + t266 * t615;
t353 = qJ(5) * t670;
t263 = -pkin(4) * t672 + t353;
t589 = t263 * t613 + t267 * t611 + t399;
t588 = t336 + t627;
t587 = t343 + t624;
t578 = -pkin(2) - t718;
t575 = t615 / 0.2e1;
t574 = -t613 / 0.2e1;
t573 = t613 / 0.2e1;
t572 = -t611 / 0.2e1;
t571 = t611 / 0.2e1;
t565 = rSges(6,1) * t427 - rSges(6,3) * t672;
t558 = qJD(4) * t640;
t557 = -qJD(1) * t263 + t426 * t609;
t556 = t427 * t599;
t553 = rSges(6,1) * t615 + rSges(6,2) * t768 + rSges(6,3) * t581;
t228 = t598 - t623;
t549 = -t228 + t569;
t547 = -t238 + t567;
t546 = t428 + t559;
t382 = rSges(2,1) * t437 - rSges(2,2) * t435;
t381 = rSges(2,1) * t435 + rSges(2,2) * t437;
t329 = rSges(5,1) * t426 - rSges(5,2) * t424;
t533 = t41 * t427 + t42 * t425;
t532 = t41 * t425 - t42 * t427;
t531 = t425 * t44 + t427 * t43;
t530 = t425 * t43 - t427 * t44;
t529 = t425 * t55 + t427 * t54;
t528 = t425 * t54 - t427 * t55;
t523 = -t425 * t87 - t427 * t86;
t138 = -rSges(5,1) * t768 - rSges(5,2) * t581 + t628;
t265 = t324 * t425;
t139 = -qJD(4) * t265 + (t329 * t427 + t410) * qJD(1);
t511 = t138 * t427 + t139 * t425;
t508 = t151 * t425 - t153 * t427;
t502 = t221 * t425 + t222 * t427;
t490 = -pkin(8) * t612 - t163 - t261;
t224 = rSges(6,2) * t670 + t565;
t485 = t224 + t494;
t483 = t546 + t270;
t475 = t611 * t630 + t627;
t471 = -t142 * t286 + t144 * t287 + t385 * t459;
t470 = (Icges(7,5) * t277 - Icges(7,6) * t278) * t286 - (Icges(7,5) * t279 + Icges(7,6) * t280) * t287 + t283 * t385;
t469 = -qJDD(5) * t426 + t120 * t611 + t121 * t613 + t302 * t266 + t424 * t602 + qJDD(2);
t466 = -t420 - t759 - t409;
t464 = t426 * t470;
t451 = t36 * t508 + (-t39 * t425 + t40 * t427) * t238;
t441 = (t459 * t427 + t510) * t286 - (t459 * t425 + t509) * t287 + (t230 + t501) * t385;
t440 = t441 * t426;
t301 = t329 * qJD(4);
t282 = t322 * t616;
t245 = t581 - t586;
t243 = t621 * t614;
t189 = -rSges(7,3) * t671 + t631;
t188 = -rSges(7,3) * t672 + t632;
t186 = t461 * t427;
t185 = t461 * t425;
t184 = t460 * t427;
t183 = t460 * t425;
t171 = rSges(7,1) * t279 + rSges(7,2) * t280;
t170 = rSges(7,1) * t277 - rSges(7,2) * t278;
t155 = qJD(1) * t173 - t401;
t154 = qJD(1) * t549 + t400;
t141 = -rSges(6,3) * t586 + t553;
t140 = t539 * t613 + (t427 * t540 + t412) * qJD(1);
t97 = qJD(4) * t502 + qJD(2);
t73 = qJDD(1) * t229 + qJD(1) * (-qJD(1) * t598 + t625) + t456;
t72 = t549 * qJDD(1) + (-qJD(1) * t229 - t262) * qJD(1) + t496;
t69 = (t223 * t427 - t224 * t425) * qJD(4) + t488;
t59 = (qJD(4) * t539 + t399) * t425 + (t223 + t493) * qJD(1) + t637;
t58 = qJD(1) * t485 + t475;
t45 = qJD(4) * t511 + t221 * t302 - t222 * t303 + qJDD(2);
t35 = qJD(1) * t138 + qJDD(1) * t222 - t301 * t613 - t302 * t324 + t454;
t34 = -t301 * t611 + t303 * t324 + (-t139 + t651) * qJD(1) + t495 * qJDD(1) + t496;
t19 = -t124 * t460 - t125 * t461 + t160 * t667 + t161 * t277 + t162 * t278 + t459 * t768;
t18 = -t122 * t460 - t123 * t461 + t160 * t670 + t161 * t279 - t162 * t280 - t459 * t474;
t17 = -t224 * t302 + t645 * t303 + (t140 * t425 + t141 * t427) * qJD(4) + t469;
t16 = qJD(1) * t141 + qJDD(1) * t223 + t302 * t630 + t425 * t558 + t444;
t15 = -t303 * t539 + t427 * t558 + t485 * qJDD(1) + (-t140 + t482) * qJD(1) + t465;
t14 = t286 * t54 - t287 * t55 + t385 * t88;
t7 = t124 * t147 - t125 * t149 - t144 * t768 + t277 * t62 + t278 * t64 + t60 * t667;
t6 = t124 * t145 + t125 * t148 - t142 * t768 + t277 * t63 + t278 * t65 + t61 * t667;
t5 = t122 * t147 - t123 * t149 + t144 * t474 + t279 * t62 - t280 * t64 + t60 * t670;
t4 = t122 * t145 + t123 * t148 + t142 * t474 + t279 * t63 - t280 * t65 + t61 * t670;
t3 = -t151 * t157 + t153 * t156 + t286 * t66 + t287 * t67 - t290 * t302 + t639 * t303 + (t199 * t425 + t200 * t427) * qJD(4) + t469;
t2 = t156 * t41 + t157 * t42 + t19 * t385 + t286 * t6 - t287 * t7 + t293 * t70;
t1 = t156 * t43 + t157 * t44 + t18 * t385 + t286 * t4 - t287 * t5 + t293 * t71;
t20 = [(-(-t58 + (t224 - t728) * qJD(1) + t475 + t754) * t59 + t58 * t587 + t59 * (-pkin(4) * t582 + t553 + t588) + t58 * (-t399 + (-t704 + (-rSges(6,3) - qJ(5)) * t426) * qJD(4)) * t425 + ((-t435 * t59 - t437 * t58) * pkin(1) + (-t58 * rSges(6,1) + t466 * t59) * t425 + (t58 * (t466 + t714) - t59 * t433) * t427) * qJD(1) + (t16 - g(2)) * (t223 + t483) + (t15 - g(1)) * ((-t402 + (rSges(6,2) - pkin(4)) * t426) * t425 + t548 + t565)) * m(6) + (-(-t154 - t305 + t400 + (-t228 - t728) * qJD(1)) * t155 + t154 * t401 + t155 * (t622 + t625) + ((-t154 * t437 - t155 * t435) * pkin(1) + t154 * (t578 + t715) * t427 + (t154 * (-rSges(4,3) - qJ(3)) + t155 * t578) * t425) * qJD(1) + (-g(2) + t73) * t173 + (-g(1) + t72) * (t425 * t578 + t404 + t623 - t728)) * m(4) + (t795 - t822) * t735 + ((t324 * t703 - t710) * qJD(4) + (t624 + (-t410 - t428 + (-t329 - t420) * t427) * qJD(1)) * t86 + (t400 + t628 + t86 - t545 - t762 + (t221 + t728 + t813) * qJD(1)) * t87 + (-g(2) + t35) * (t222 + t546) + (-g(1) + t34) * (-t760 + t813)) * m(5) + (t797 + t798) * t573 + (t796 - t799 + t800) * t572 - t723 / 0.2e1 + t724 / 0.2e1 + t719 + (((-t759 - t413) * t425 + t784) * t11 + (t342 - t543 + t587 + (rSges(7,3) * t614 - qJ(5) * t612 - t399) * t425 + (t427 * t758 - t417 - t428) * qJD(1)) * t39 - g(1) * t784 - (t426 * t729 - t402) * t726 + (-g(2) + t10) * (t483 + t657) + (t556 * t614 + t39 + t396 + t588 + t662 - t754 + (t425 * t758 - t290 - t391) * qJD(1) + t804) * t40) * m(7) + (t828 * qJD(4) + t833 * t424 + t834 * t426) * qJD(1) + ((-t326 * t439 - g(2) + t552) * t292 + (-t802 + (-0.2e1 * t331 - t428 + t292) * t439 - g(1)) * t291) * m(3) + (t109 + t779) * t734 - t303 * t112 / 0.2e1 + (m(2) * (t381 ^ 2 + t382 ^ 2) + m(3) * (t291 ^ 2 + t331 * t292) + Icges(2,3) + Icges(3,3) + Icges(4,2) * t432 ^ 2 + (Icges(4,1) * t431 + 0.2e1 * Icges(4,4) * t432) * t431 - t827) * qJDD(1) + t701 / 0.2e1 + t702 / 0.2e1 + (t738 + t737) * t12 + (((t427 * t791 + t780 - t788) * t427 + (t425 * t791 + t781 + t840) * t425) * qJD(4) + t801 - t806) * t574 - m(2) * (-g(1) * t381 + g(2) * t382) + t18 * t738 + t19 * t739 + t71 * t741 + t70 * t742 + ((t788 * t425 + ((t824 + t846) * t427 + t782 + t839 + t845) * t427) * qJD(4) + t807) * t571; (m(3) + m(4)) * qJDD(2) + m(5) * t45 + m(6) * t17 + m(7) * t3 + (-m(3) + t592) * g(3); t592 * (-g(2) * t427 + t726) + 0.2e1 * (-t712 / 0.2e1 + t711 / 0.2e1) * m(7) + 0.2e1 * (t15 * t731 + t16 * t730) * m(6) + 0.2e1 * (t34 * t731 + t35 * t730) * m(5) + 0.2e1 * (t72 * t731 + t73 * t730) * m(4); (t12 * t427 + t13 * t425) * t608 / 0.2e1 + (t58 * t282 + t17 * t644 + t69 * t590 + (t15 * t630 + t58 * t640 + t17 * t223 + t69 * t141 + (-t69 * t224 + t59 * t630) * qJD(1)) * t427 + (t16 * t630 + t59 * t640 - t17 * t224 + t69 * t140 + (-t539 * t58 + t645 * t69) * qJD(1)) * t425 - g(1) * (t267 + t268) - g(2) * (t263 + t264) + g(3) * t629 - t58 * (-qJD(1) * t264 + t557) - t59 * (qJD(1) * t268 + t641) - t69 * t589 - ((t69 * t268 + t58 * t629) * t427 + (t69 * t264 + t59 * t629) * t425) * qJD(4)) * m(6) + (t799 * t427 + t798 * t425 + (t425 * t779 - t427 * t822) * qJD(1)) * qJD(1) / 0.2e1 + (-t425 * t822 - t427 * t779) * qJDD(1) / 0.2e1 + t809 * t734 + t810 * t735 + (t3 * t644 + (t10 * t547 + t3 * t656 + t40 * t490) * t425 + (t3 * t657 + (qJD(1) * t40 + t11) * t547) * t427 - g(1) * (t355 + t631) - g(2) * (t353 + t632) - g(3) * (t237 - t761) - t40 * (t151 * t607 + t189 * t385 - t237 * t286 + t641) - t757 * qJD(4) * t761 + (-g(1) * t556 - t599 * t725 - (-t621 * t681 - t591) * pkin(8) - t451 * qJD(6)) * t424 + (t153 * t607 + t188 * t385 + t237 * t287 + t238 * t616 + t427 * t490 + t282 - t557) * t39 + (t590 + (t199 + t66 + (-t151 + t639) * qJD(1)) * t425 + (qJD(1) * t656 + t200 + t67) * t427 - t188 * t286 - t189 * t287 - t589) * t36) * m(7) + (qJD(1) * t797 + t771 * qJD(4) + qJDD(1) * t795 + t780 * t302 + t781 * t303 + t2) * t731 + (t12 + t800) * t575 + (t13 + t801) * t616 / 0.2e1 - ((t424 * t792 + t426 * t787) * qJD(4) + (t424 * t789 + t426 * t790) * qJD(1)) * qJD(1) / 0.2e1 + ((t425 * t783 + t427 * t782) * qJD(1) + t770) * t572 - t14 * t607 / 0.2e1 + (((-t184 * t436 - t186 * t434 + t142) * t286 - (-t183 * t436 - t185 * t434 + t144) * t287 + (-t232 * t436 - t234 * t434 - t459) * t385 + t88 * qJD(6)) * t426 + (-qJD(6) * t529 + t441) * t424) * t733 + (-(t265 * t86 - t710) * qJD(1) - (t97 * (-t265 * t425 - t269 * t427) + t523 * t329) * qJD(4) + t45 * t502 + t97 * ((t221 * t427 - t222 * t425) * qJD(1) + t511) + t523 * t301 + (-t34 * t427 - t35 * t425 + (-t427 * t87 + t703) * qJD(1)) * t324 + g(1) * t269 + g(2) * t265 - g(3) * t329) * m(5) + (qJD(1) * t529 + t425 * t8 - t427 * t9) * t732 + t528 * t736 + ((t184 * t279 - t186 * t280) * t286 - (t183 * t279 - t185 * t280) * t287 + (t232 * t279 - t234 * t280) * t385 + (t426 * t71 - t43 * t671) * qJD(6) + ((-qJD(6) * t44 + t471) * t424 + t440) * t425) * t737 + (qJD(1) * t531 + t4 * t425 - t427 * t5) * t738 + (qJD(1) * t533 + t425 * t6 - t427 * t7) * t739 + ((t184 * t277 + t186 * t278) * t286 - (t183 * t277 + t185 * t278) * t287 + (t232 * t277 + t234 * t278) * t385 + (-t42 * t672 + t426 * t70) * qJD(6) + ((-qJD(6) * t41 + t471) * t424 + t440) * t427) * t740 + t530 * t741 + t532 * t742 + (t796 * qJD(1) + t770 * qJD(4) + qJDD(1) * t825 + t782 * t302 + t783 * t303 + t1) * t730 + ((-t611 * t764 - t766) * t427 + ((t427 * t765 + t773) * qJD(4) + t772) * t425) * t571 + ((-t613 * t765 + t766) * t425 + ((t425 * t764 + t773) * qJD(4) + t772) * t427) * t574 + ((t425 * t781 + t427 * t780) * qJD(1) + t771) * t573; t743 * (-g(3) * t426 + (g(1) * t427 + t725) * t424) - m(6) * (t243 * t69 + t244 * t59 + t245 * t58) - m(7) * (t243 * t36 + t244 * t40 + t245 * t39) + 0.2e1 * ((t58 * t611 + t59 * t613 - t17) * t745 + (t39 * t611 + t40 * t613 - t3) * t744) * t426 + 0.2e1 * ((qJD(4) * t69 + t15 * t427 + t16 * t425 - t58 * t616 + t59 * t615) * t745 + (t10 * t425 + t11 * t427 - t39 * t616 + t591 + t681) * t744) * t424; t2 * t667 / 0.2e1 + (t424 * t70 + t426 * t533) * t742 + ((-qJD(4) * t533 + t19) * t424 + (-qJD(1) * t532 + qJD(4) * t70 + t425 * t7 + t427 * t6) * t426) * t739 + t1 * t670 / 0.2e1 + (t424 * t71 + t426 * t531) * t741 + ((-qJD(4) * t531 + t18) * t424 + (-qJD(1) * t530 + qJD(4) * t71 + t4 * t427 + t425 * t5) * t426) * t738 + t14 * t612 / 0.2e1 + t424 * (t701 + t702 + t719 - t723 + t724) / 0.2e1 + (t424 * t88 + t426 * t529) * t736 + ((-qJD(4) * t529 + t29) * t424 + (-qJD(1) * t528 + qJD(4) * t88 + t425 * t9 + t427 * t8) * t426) * t732 + (t452 * t277 - t278 * t453 + t427 * t464) * t740 + (t279 * t452 + t280 * t453 + t425 * t464) * t737 + (t470 * t424 + (t453 * t434 - t436 * t452) * t426) * t733 + (t424 * t574 + t426 * t575) * t13 + (-t585 / 0.2e1 + t424 * t572) * t12 + ((qJD(4) * t451 + t10 * t151 - t11 * t153 - t39 * t66 + t40 * t67) * t424 + (t39 * (-qJD(4) * t153 + t163 * t425) + t40 * (qJD(4) * t151 - t163 * t427) - t3 * t508 + t36 * (-t151 * t615 - t153 * t616 - t425 * t67 + t427 * t66) + (qJD(1) * t757 + t711 - t712) * t238) * t426 - t39 * (-t171 * t385 - t287 * t288) - t40 * (t170 * t385 - t286 * t288) - t36 * (t170 * t287 + t171 * t286) - g(1) * t170 - g(2) * t171 - g(3) * t288) * m(7);];
tau  = t20;
