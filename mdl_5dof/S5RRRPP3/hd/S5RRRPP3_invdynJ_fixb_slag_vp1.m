% Calculate vector of inverse dynamics joint torques for
% S5RRRPP3
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRPP3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_invdynJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:08
% EndTime: 2019-12-31 20:53:40
% DurationCPUTime: 26.93s
% Computational Cost: add. (14811->693), mult. (16372->782), div. (0->0), fcn. (12713->6), ass. (0->388)
t836 = Icges(5,4) - Icges(4,5);
t835 = Icges(5,5) - Icges(4,6);
t425 = sin(qJ(3));
t427 = cos(qJ(3));
t834 = t835 * t425 - t836 * t427;
t833 = Icges(5,1) + Icges(4,3);
t666 = Icges(4,4) * t425;
t334 = Icges(4,2) * t427 + t666;
t415 = Icges(6,6) * t425;
t325 = -Icges(6,2) * t427 + t415;
t663 = Icges(5,6) * t425;
t508 = Icges(5,3) * t427 + t663;
t825 = t325 - t508;
t832 = -t334 + t825;
t831 = Icges(5,4) - Icges(6,5);
t830 = Icges(6,4) + Icges(5,5);
t424 = qJ(1) + qJ(2);
t412 = sin(t424);
t330 = Icges(6,4) * t425 + Icges(6,5) * t427;
t413 = cos(t424);
t478 = t330 * t413;
t208 = Icges(6,1) * t412 + t478;
t828 = t834 * t413;
t795 = t833 * t412 + t828;
t829 = t208 + t795;
t638 = t412 * t425;
t378 = Icges(4,4) * t638;
t637 = t412 * t427;
t198 = Icges(4,1) * t637 - Icges(4,5) * t413 - t378;
t367 = Icges(6,6) * t638;
t201 = Icges(6,5) * t413 - Icges(6,3) * t637 - t367;
t827 = -t198 + t201;
t337 = Icges(4,1) * t427 - t666;
t481 = t337 * t413;
t199 = Icges(4,5) * t412 + t481;
t506 = Icges(6,3) * t427 + t415;
t475 = t506 * t413;
t200 = Icges(6,5) * t412 + t475;
t816 = t199 + t200;
t662 = Icges(5,6) * t427;
t509 = -Icges(5,3) * t425 + t662;
t476 = t509 * t413;
t202 = Icges(5,5) * t412 - t476;
t634 = t413 * t427;
t369 = Icges(6,6) * t634;
t635 = t413 * t425;
t204 = Icges(6,4) * t412 + Icges(6,2) * t635 + t369;
t815 = t202 + t204;
t203 = Icges(5,5) * t413 + Icges(5,6) * t637 - Icges(5,3) * t638;
t368 = Icges(6,6) * t637;
t205 = Icges(6,4) * t413 - Icges(6,2) * t638 - t368;
t826 = t203 + t205;
t794 = t833 * t413 + t836 * t637 - t835 * t638;
t661 = Icges(6,6) * t427;
t507 = -Icges(6,3) * t425 + t661;
t511 = Icges(5,2) * t425 + t662;
t824 = t507 - t511;
t416 = Icges(4,4) * t427;
t336 = Icges(4,1) * t425 + t416;
t783 = t511 + t336;
t823 = t507 - t783;
t196 = Icges(4,4) * t637 - Icges(4,2) * t638 - Icges(4,6) * t413;
t370 = Icges(5,6) * t638;
t207 = Icges(5,4) * t413 + Icges(5,2) * t637 - t370;
t822 = t196 * t425 - t207 * t427;
t821 = -t199 * t637 - t202 * t638;
t819 = t198 * t427 - t203 * t425 - t822;
t818 = t831 * t425 + t830 * t427;
t324 = Icges(6,2) * t425 + t661;
t516 = -Icges(4,2) * t425 + t416;
t807 = t324 - t509 - t516;
t328 = Icges(4,5) * t425 + Icges(4,6) * t427;
t806 = t328 - t818;
t512 = Icges(5,2) * t427 - t663;
t817 = t337 + t506 + t512;
t529 = -t200 * t637 - t204 * t638 + t208 * t413;
t480 = t516 * t413;
t197 = Icges(4,6) * t412 + t480;
t371 = Icges(5,6) * t635;
t206 = Icges(5,4) * t412 - Icges(5,2) * t634 + t371;
t814 = -t795 * t413 - t821;
t766 = -t197 * t638 - t206 * t637 + t814;
t714 = -t529 + t766;
t813 = t829 * t412 + t816 * t634 + t815 * t635;
t209 = Icges(6,1) * t413 - Icges(6,4) * t638 - Icges(6,5) * t637;
t186 = t412 * t209;
t812 = t794 * t412 + t827 * t634 + t826 * t635 + t186;
t811 = t196 + t826;
t810 = t197 - t815;
t809 = t207 - t827;
t808 = -t206 + t816;
t423 = qJD(1) + qJD(2);
t804 = (Icges(4,6) - t830) * t423 + t832 * qJD(3);
t803 = (Icges(4,5) - t831) * t423 + t823 * qJD(3);
t490 = t334 * t425 - t336 * t427;
t798 = t825 * t425 - t824 * t427 - t490;
t802 = t197 * t425 + t206 * t427;
t801 = t818 * t413;
t716 = -t196 * t635 + t207 * t634 - t812;
t715 = -t197 * t635 - t206 * t634 + t813;
t647 = t209 * t413;
t501 = t201 * t427 + t205 * t425;
t733 = t412 * t501;
t84 = t647 - t733;
t713 = t819 * t412 + t794 * t413 + t84;
t800 = t817 * qJD(3);
t799 = t807 * qJD(3);
t788 = t330 + t834;
t797 = t823 * t425 + t832 * t427;
t709 = t806 * t412;
t765 = rSges(6,3) + qJ(5);
t643 = t328 * t413;
t155 = -t412 * t490 - t643;
t760 = t824 * t637 - t638 * t825 - t801;
t796 = t155 - t760;
t793 = -t325 + t334;
t712 = t809 * t425 + t811 * t427;
t711 = t808 * t425 + t810 * t427;
t639 = t412 * t423;
t792 = t804 * t413 + t807 * t639;
t636 = t413 * t423;
t791 = -t324 * t636 + (t476 + t480) * t423 + t804 * t412;
t790 = t803 * t413 - t639 * t817;
t789 = t512 * t636 + (t475 + t481) * t423 + t803 * t412;
t748 = t798 * t413 + t709;
t787 = (Icges(6,1) + t833) * t423 - t806 * qJD(3);
t786 = -t501 + t819;
t785 = t815 * t425 + t816 * t427 - t802;
t679 = rSges(6,1) + pkin(4);
t759 = -rSges(6,2) * t638 + t413 * t679 - t765 * t637;
t784 = t759 * t423;
t782 = t797 * qJD(3) + t806 * t423 + t799 * t425 + t800 * t427;
t781 = t715 * t412 - t716 * t413;
t780 = t714 * t412 - t713 * t413;
t779 = -t788 * qJD(3) + t798 * t423;
t405 = t413 * pkin(7);
t287 = pkin(2) * t412 - t405;
t342 = pkin(7) * t636;
t778 = t423 * t287 + t342;
t776 = t748 * t423;
t775 = -t711 * qJD(3) + t829 * t423 - t792 * t425 + t790 * t427;
t774 = -t789 * t427 + t791 * t425 + (t209 + t794) * t423 + t712 * qJD(3);
t773 = -t807 - t823;
t772 = -t508 - t793 + t817;
t771 = Icges(5,3) * t634 + t413 * t793 + t371 - t808;
t770 = -Icges(6,3) * t635 - t413 * t783 + t369 - t810;
t769 = t796 * t423;
t768 = -t413 * t787 + t423 * t785 + t639 * t788;
t767 = (t478 - t786 + t828) * t423 + t787 * t412;
t758 = t794 + t802;
t578 = qJD(3) * qJD(4);
t757 = qJDD(4) * t425 + t427 * t578;
t582 = qJD(3) * t427;
t557 = t413 * t582;
t633 = t423 * t425;
t573 = t412 * t633;
t234 = t557 - t573;
t583 = qJD(3) * t425;
t558 = t413 * t583;
t632 = t423 * t427;
t756 = t412 * t632 + t558;
t560 = t412 * t582;
t571 = t413 * t633;
t233 = t560 + t571;
t414 = t425 * qJ(4);
t420 = t427 * pkin(3);
t725 = t420 + t414;
t261 = t725 * t412;
t230 = t423 * t261;
t411 = qJD(4) * t425;
t358 = t413 * t411;
t601 = qJ(4) * t557 + t358;
t755 = t601 + t230 + t778;
t426 = sin(qJ(1));
t673 = pkin(1) * qJD(1);
t574 = t426 * t673;
t285 = rSges(3,1) * t412 + rSges(3,2) * t413;
t645 = t285 * t423;
t225 = -t574 - t645;
t754 = qJD(3) * t781 + t776;
t753 = qJD(3) * t780 + t769;
t752 = qJD(3) * t785 + t425 * t790 + t427 * t792;
t751 = qJD(3) * t786 + t425 * t789 + t427 * t791;
t750 = -t412 * t779 + t413 * t782;
t749 = t412 * t782 + t413 * t779;
t747 = -t647 + t813;
t561 = t412 * t583;
t311 = qJ(5) * t561;
t579 = qJD(5) * t427;
t746 = t412 * t579 - t311;
t745 = qJ(5) * t634 + t412 * t679;
t584 = qJD(3) * t423;
t269 = qJDD(3) * t412 + t413 * t584;
t356 = t413 * t579;
t422 = qJDD(1) + qJDD(2);
t580 = qJD(5) * t425;
t581 = qJD(4) * t427;
t278 = qJD(3) * t725 - t581;
t418 = t425 * rSges(6,2);
t726 = t427 * rSges(6,3) + t418;
t604 = -qJD(3) * t726 - t278;
t444 = (-qJ(5) * qJD(3) ^ 2 + qJDD(5)) * t427 + (-0.2e1 * t580 + t604) * qJD(3);
t146 = -pkin(3) * t756 - qJ(4) * t573 + t601;
t266 = pkin(3) * t634 + t413 * t414;
t288 = t413 * pkin(2) + t412 * pkin(7);
t428 = cos(qJ(1));
t421 = t428 * pkin(1);
t430 = qJD(1) ^ 2;
t678 = pkin(1) * t426;
t533 = qJDD(1) * t421 - t430 * t678;
t483 = t423 * (-pkin(2) * t639 + t342) + t422 * t288 + t533;
t446 = t422 * t266 + t483 + (t146 + t358) * t423 + t757 * t412;
t344 = pkin(3) * t425 - qJ(4) * t427;
t350 = -rSges(6,2) * t427 + rSges(6,3) * t425;
t659 = qJ(5) * t425;
t545 = t350 + t659;
t531 = -t344 - t545;
t610 = rSges(6,2) * t635 + rSges(6,3) * t634 + t745;
t721 = rSges(6,2) * t557 + t679 * t636 + t356;
t630 = -rSges(6,2) * t573 - t756 * t765 + t721;
t8 = t610 * t422 + (t356 + t630) * t423 + t531 * t269 + t444 * t412 + t446;
t744 = t8 - g(2);
t270 = -qJDD(3) * t413 + t412 * t584;
t474 = (-qJDD(1) * t426 - t428 * t430) * pkin(1);
t448 = t270 * t344 + t757 * t413 + t474;
t605 = -t261 - t287;
t538 = t605 + t759;
t319 = pkin(3) * t561;
t556 = t412 * t411;
t570 = t413 * t632;
t147 = pkin(3) * t570 + qJ(4) * t233 - t319 + t556;
t238 = t288 * t423;
t629 = -t147 - t238;
t313 = rSges(6,3) * t561;
t631 = rSges(6,2) * t560 - t313 + (t413 * t726 + t745) * t423 + t746;
t734 = t412 * (t579 + t411);
t9 = t545 * t270 + t538 * t422 + t444 * t413 + (t629 - t631 - t734) * t423 + t448;
t743 = t9 - g(1);
t564 = rSges(4,1) * t561 + t233 * rSges(4,2);
t728 = rSges(4,1) * t634 + t412 * rSges(4,3);
t141 = t423 * t728 - t564;
t676 = rSges(4,1) * t427;
t352 = -rSges(4,2) * t425 + t676;
t306 = t352 * qJD(3);
t347 = rSges(4,1) * t425 + rSges(4,2) * t427;
t585 = qJD(3) * t413;
t588 = rSges(4,2) * t638 + t413 * rSges(4,3);
t218 = rSges(4,1) * t637 - t588;
t611 = -t218 - t287;
t39 = -t306 * t585 + t270 * t347 + (-t141 - t238) * t423 + t611 * t422 + t474;
t742 = -g(1) + t39;
t482 = -t234 * rSges(4,2) + rSges(4,3) * t636;
t140 = -rSges(4,1) * t756 + t482;
t219 = -rSges(4,2) * t635 + t728;
t586 = qJD(3) * t412;
t40 = t140 * t423 + t219 * t422 - t269 * t347 - t306 * t586 + t483;
t741 = -g(2) + t40;
t303 = rSges(5,2) * t570;
t675 = rSges(5,2) * t425;
t527 = rSges(5,3) * t427 + t675;
t143 = rSges(5,3) * t571 - t303 + (rSges(5,1) * t423 + qJD(3) * t527) * t412;
t417 = t425 * rSges(5,3);
t727 = -rSges(5,2) * t427 + t417;
t603 = -qJD(3) * t727 - t278;
t540 = qJD(3) * t603;
t546 = rSges(5,1) * t413 - rSges(5,3) * t638;
t223 = rSges(5,2) * t637 + t546;
t565 = t223 + t605;
t23 = -t270 * t527 + t413 * t540 + t565 * t422 + (-t143 - t556 + t629) * t423 + t448;
t740 = t23 - g(1);
t537 = rSges(5,1) * t636 + t756 * rSges(5,2) + rSges(5,3) * t557;
t145 = -rSges(5,3) * t573 + t537;
t221 = t412 * rSges(5,1) - rSges(5,2) * t634 + rSges(5,3) * t635;
t592 = -t344 + t527;
t24 = t145 * t423 + t221 * t422 + t269 * t592 + t412 * t540 + t446;
t739 = t24 - g(2);
t236 = rSges(3,1) * t636 - rSges(3,2) * t639;
t738 = -t236 * t423 - t285 * t422 - g(1) + t474;
t286 = t413 * rSges(3,1) - rSges(3,2) * t412;
t737 = t286 * t422 - t423 * t645 - g(2) + t533;
t469 = t531 * t585;
t590 = t356 + t358;
t66 = t423 * t538 + t469 - t574 + t590;
t736 = t413 * t66;
t488 = -pkin(2) - t725;
t668 = -rSges(5,3) - qJ(4);
t484 = t585 * t592 + t358;
t461 = t484 - t574;
t73 = t423 * t565 + t461;
t670 = t423 * t73;
t279 = t344 * t586;
t486 = t527 * t586 - t279 + t556;
t575 = t428 * t673;
t536 = t266 + t288;
t732 = t221 + t536;
t74 = t423 * t732 + t486 + t575;
t735 = (t73 * (-t411 + (t427 * t668 - t675) * qJD(3)) + (t73 * (-rSges(5,1) - pkin(7)) + t74 * (t488 - t417)) * t423) * t412 + (-t74 * pkin(3) * t583 + (t425 * t668 - pkin(2) - t420) * t670) * t413;
t528 = t261 * t586 + t266 * t585 - t581;
t53 = t580 + (-t412 * t759 + t610 * t413) * qJD(3) + t528;
t657 = qJD(3) * t53;
t361 = qJ(4) * t637;
t258 = -pkin(3) * t638 + t361;
t259 = rSges(5,2) * t638 + rSges(5,3) * t637;
t730 = t258 + t259;
t729 = -t350 * t586 - t279;
t576 = -pkin(3) - t765;
t534 = t576 * t427;
t535 = t413 * t576;
t547 = -pkin(2) - t414;
t669 = -rSges(6,2) - qJ(4);
t567 = -t266 - t610;
t539 = t288 - t567;
t67 = t423 * t539 - t311 + t575 + t729 + t734;
t724 = ((t425 * t669 - pkin(2) + t534) * t736 + (t66 * (-pkin(7) - t679) + (t547 - t418 + t534) * t67) * t412) * t423 + (t425 * t535 * t67 + t637 * t66 * t669) * qJD(3) - t469 * t67;
t723 = -t590 + t721 + t755 - t784;
t722 = -qJ(5) * t582 + t604;
t720 = -t767 * t412 + t774 * t413;
t719 = -t768 * t412 + t775 * t413;
t718 = t774 * t412 + t767 * t413;
t717 = t775 * t412 + t768 * t413;
t189 = t423 * t218;
t710 = -rSges(4,1) * t558 + t189 + t482 + t778;
t708 = -t643 + t801;
t706 = (-t773 * t425 + t772 * t427) * t423;
t705 = t771 * t425 + t770 * t427;
t704 = t367 - t370 - t378 + (-Icges(4,2) - Icges(6,2) - Icges(5,3)) * t637 + t809;
t703 = Icges(6,3) * t638 + t412 * t783 - t368 + t811;
t702 = t788 * t423;
t166 = t219 + t288;
t695 = t166 * t423 - t347 * t586;
t691 = -t423 * t223 + t537 + t755;
t690 = t425 * t704 + t427 * t703;
t689 = m(5) / 0.2e1;
t688 = m(6) / 0.2e1;
t687 = t269 / 0.2e1;
t685 = t270 / 0.2e1;
t677 = g(2) * t412;
t671 = t423 * t66;
t562 = t347 * t585;
t470 = -t562 - t574;
t103 = t423 * t611 + t470;
t656 = t103 * t413;
t609 = -t221 - t266;
t607 = t412 * t261 + t413 * t266;
t364 = qJ(4) * t634;
t263 = -pkin(3) * t635 + t364;
t606 = t423 * t263 + t412 * t581;
t602 = t303 + t319;
t591 = -t725 - t727;
t264 = rSges(5,2) * t635 + rSges(5,3) * t634;
t587 = t412 ^ 2 + t413 ^ 2;
t569 = t412 * t147 + (t146 + t230) * t413;
t568 = t258 * t586 + t263 * t585 + t411;
t553 = -pkin(2) - t676;
t551 = -t586 / 0.2e1;
t550 = t586 / 0.2e1;
t549 = -t585 / 0.2e1;
t548 = t585 / 0.2e1;
t532 = g(1) * t413 + t677;
t353 = rSges(2,1) * t428 - rSges(2,2) * t426;
t348 = rSges(2,1) * t426 + rSges(2,2) * t428;
t104 = t575 + t695;
t505 = -t104 * t412 - t656;
t496 = t218 * t412 + t219 * t413;
t473 = -qJDD(4) * t427 + t146 * t585 + t147 * t586 + t269 * t261 + t425 * t578;
t165 = t412 * t553 + t405 + t588;
t464 = -t556 - t746;
t102 = t536 + t610;
t447 = t313 + t319 + t464;
t151 = t405 + ((rSges(5,2) - pkin(3)) * t427 + t547) * t412 + t546;
t101 = t412 * t488 + t405 + t759;
t434 = (t553 * t656 + (t103 * (-rSges(4,3) - pkin(7)) + t104 * t553) * t412) * t423;
t433 = -t760 * t270 / 0.2e1 + (((t84 + t733 + t747) * t412 + ((t795 + t822) * t413 + t766 + t812 + t821) * t413) * qJD(3) + t776) * t548 + (t798 * qJD(3) + t800 * t425 - t799 * t427) * t423 + (t155 + t712) * t685 + (Icges(3,3) - t797) * t422 + (t748 + t711) * t687 + (t750 + t752) * t550 + (((t758 * t413 + t715 - t747) * t413 + (t758 * t412 + t186 + t529 + t716 - t814) * t412) * qJD(3) + t753 - t769) * t551 + (t749 + t751 + t754) * t549;
t391 = rSges(6,2) * t634;
t385 = rSges(6,2) * t637;
t359 = t413 * t581;
t268 = t344 * t639;
t267 = -rSges(6,3) * t635 + t391;
t265 = t347 * t413;
t262 = -rSges(6,3) * t638 + t385;
t260 = t347 * t412;
t239 = t587 * t583;
t226 = t286 * t423 + t575;
t108 = t496 * qJD(3);
t72 = (t221 * t413 - t223 * t412) * qJD(3) + t528;
t10 = -t223 * t269 + t609 * t270 + (t143 * t412 + t145 * t413) * qJD(3) + t473;
t7 = qJDD(5) * t425 - t759 * t269 + t567 * t270 + (t631 * t412 + t630 * t413 + t579) * qJD(3) + t473;
t1 = [Icges(2,3) * qJDD(1) + t433 + (t737 * (t286 + t421) + t738 * (-t285 - t678) + (-t236 - t575 + t226) * t225) * m(3) + (g(1) * t348 - g(2) * t353 + (t348 ^ 2 + t353 ^ 2) * qJDD(1)) * m(2) + (t66 * (t447 - t575) + t744 * (t421 + t102) + t743 * (t101 - t678) + (t66 + t723) * t67 + t724) * m(6) + (t73 * (-t575 + t602) + (-t574 + t73 - t461 + t691) * t74 + t739 * (t732 + t421) + t740 * (t151 - t678) + t735) * m(5) + (t103 * (t564 - t575) + t434 + t741 * (t166 + t421) + t742 * (t165 - t678) + (t103 - t470 - t574 + t710) * t104) * m(4); t433 + (t539 * t671 + t723 * t67 + (t447 - t464 + t729) * t66 + t744 * t102 + t743 * t101 + t724) * m(6) + ((-t484 + t691) * t74 + (t602 + t486) * t73 + t740 * t151 + t735 + (t670 + t739) * t732) * m(5) + (t434 + t741 * t166 + t742 * t165 + (t562 + t710) * t104 + (t564 + t695) * t103) * m(4) + (-t225 * t236 - t226 * t645 + (t225 * t423 + t737) * t286 + (t226 * t423 - t738) * t285) * m(3); t781 * t687 + t780 * t685 + (t750 * t423 + t748 * t422 + t716 * t270 + t715 * t269 + (t719 * t412 + t720 * t413) * qJD(3)) * t412 / 0.2e1 - (t749 * t423 + t796 * t422 + t713 * t270 + t714 * t269 + (t717 * t412 + t413 * t718) * qJD(3)) * t413 / 0.2e1 + (t412 * t711 - t413 * t712) * t422 / 0.2e1 - ((t772 * t425 + t773 * t427) * t423 + ((-t771 * t412 - t704 * t413) * t427 + (t770 * t412 + t703 * t413) * t425) * qJD(3)) * t423 / 0.2e1 + ((t423 * t711 - t751) * t413 + (t423 * t712 + t752) * t412) * t423 / 0.2e1 + t753 * t639 / 0.2e1 + t754 * t636 / 0.2e1 + ((t586 * t708 + t702) * t412 + ((t690 * t413 + (t705 + t709) * t412) * qJD(3) + t706) * t413) * t551 + ((t423 * t715 + t720) * t413 + (t423 * t716 + t719) * t412) * t550 + ((t423 * t714 + t718) * t413 + (t423 * t713 + t717) * t412) * t549 + ((-t585 * t709 - t702) * t413 + ((t705 * t412 + (t690 - t708) * t413) * qJD(3) + t706) * t412) * t548 + (-g(1) * (t364 + t391) - g(2) * (t361 + t385) - (g(1) * t535 + t576 * t677) * t425 + t7 * t607 + (t8 * t531 - t7 * t759) * t412 + (t9 * t531 + t7 * t610) * t413 + (-t606 + (qJ(5) * t635 + t531 * t413 - t267) * t423 + t722 * t412) * t67 + (-t359 + t268 + (t350 * t412 + t258 + t262) * t423 + t722 * t413) * t66 + (-t568 - t579 - (t262 * t412 + t267 * t413 - t587 * t659) * qJD(3) + t569 + (t567 * t423 + t631) * t412 + (t630 - t784) * t413) * t53 + (-g(3) + (t412 * t67 + t736) * qJD(3)) * (qJ(5) * t427 + t725 + t726)) * m(6) + (-t73 * (-t423 * t730 + t359) - t74 * (t264 * t423 + t606) - t72 * t568 - ((t72 * t264 + t591 * t73) * t413 + (t72 * t259 + t591 * t74) * t412) * qJD(3) + t73 * t268 + t10 * t607 + t72 * t569 + (t23 * t592 + t73 * t603 + t10 * t221 + t72 * t145 + (-t72 * t223 + t592 * t74) * t423) * t413 + (t24 * t592 + t74 * t603 - t10 * t223 + t72 * t143 + (-t527 * t73 + t609 * t72) * t423) * t412 - g(1) * (t263 + t264) - g(2) * t730 + g(3) * t591) * m(5) + (-(t103 * t260 - t104 * t265) * t423 - (t108 * (-t260 * t412 - t265 * t413) + t505 * t352) * qJD(3) + (t218 * t269 - t219 * t270 + (t140 * t413 + t141 * t412) * qJD(3)) * t496 + t108 * ((t140 + t189) * t413 + (-t219 * t423 + t141) * t412) + t505 * t306 + ((-t104 * t423 - t39) * t413 + (t103 * t423 - t40) * t412) * t347 + g(1) * t265 + g(2) * t260 - g(3) * t352) * m(4); (-m(5) - m(6)) * (-g(3) * t427 + t425 * t532) - m(5) * (t233 * t74 + t234 * t73 + t239 * t72) - m(6) * (t233 * t67 + t234 * t66 + t239 * t53) + 0.2e1 * ((t585 * t73 + t586 * t74 - t10) * t689 + (t585 * t66 + t586 * t67 - t7) * t688) * t427 + 0.2e1 * ((qJD(3) * t72 + t23 * t413 + t24 * t412 + t636 * t74 - t639 * t73) * t689 + (t412 * t8 + t413 * t9 + t636 * t67 - t639 * t66 + t657) * t688) * t425; (-(-t412 * t66 + t413 * t67) * t632 + (t7 - g(3)) * t425 + (t657 + (t423 * t67 + t9) * t413 + (t8 - t671) * t412 - t587 * t657 - t532) * t427) * m(6);];
tau = t1;
