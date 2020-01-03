% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPP4_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP4_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP4_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP4_coriolismatJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP4_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP4_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP4_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:54:54
% EndTime: 2019-12-31 20:55:19
% DurationCPUTime: 16.41s
% Computational Cost: add. (42316->596), mult. (38773->764), div. (0->0), fcn. (36100->8), ass. (0->367)
t544 = qJ(2) + qJ(3);
t520 = pkin(8) + t544;
t516 = sin(t520);
t517 = cos(t520);
t521 = sin(t544);
t522 = cos(t544);
t885 = Icges(4,5) * t522 + Icges(5,5) * t517 - Icges(4,6) * t521 - Icges(5,6) * t516;
t884 = Icges(4,3) + Icges(5,3);
t548 = cos(qJ(1));
t883 = t885 * t548;
t546 = sin(qJ(1));
t705 = t522 * t546;
t709 = t521 * t546;
t713 = t517 * t546;
t719 = t516 * t546;
t875 = Icges(4,5) * t705 + Icges(5,5) * t713 - Icges(4,6) * t709 - Icges(5,6) * t719 - t884 * t548;
t882 = t884 * t546 + t883;
t738 = Icges(5,4) * t516;
t452 = Icges(5,1) * t517 - t738;
t348 = Icges(5,5) * t546 + t452 * t548;
t506 = Icges(6,5) * t516;
t832 = Icges(6,1) * t517 + t506;
t881 = Icges(6,4) * t546 + t548 * t832 + t348;
t507 = Icges(5,4) * t517;
t451 = Icges(5,1) * t516 + t507;
t737 = Icges(6,5) * t517;
t880 = Icges(6,1) * t516 + t451 - t737;
t574 = Icges(6,3) * t517 - t506;
t873 = Icges(5,2) * t517 + t574 + t738;
t739 = Icges(4,4) * t521;
t475 = Icges(4,1) * t522 - t739;
t400 = Icges(4,5) * t546 + t475 * t548;
t879 = -t348 * t713 - t400 * t705;
t708 = t521 * t548;
t610 = -rSges(4,2) * t708 + rSges(4,3) * t546;
t629 = rSges(4,2) * t709 + t548 * rSges(4,3);
t704 = t522 * t548;
t262 = t546 * (rSges(4,1) * t705 - t629) + t548 * (rSges(4,1) * t704 + t610);
t549 = -pkin(7) - pkin(6);
t519 = t546 * t549;
t541 = t548 * pkin(6);
t679 = t548 * t549;
t547 = cos(qJ(2));
t540 = t547 * pkin(2);
t518 = t540 + pkin(1);
t759 = pkin(1) - t518;
t655 = -t546 * (t759 * t546 - t541 - t679) + t548 * (-pkin(6) * t546 - t759 * t548 - t519);
t178 = t262 + t655;
t476 = rSges(4,1) * t521 + rSges(4,2) * t522;
t441 = t476 * t546;
t442 = t476 * t548;
t281 = -t546 * t441 - t548 * t442;
t754 = rSges(4,1) * t522;
t477 = -rSges(4,2) * t521 + t754;
t545 = sin(qJ(2));
t762 = pkin(2) * t545;
t558 = t476 + t762;
t835 = t558 * t548;
t836 = t558 * t546;
t104 = t178 * t281 + (t546 * t836 + t548 * t835) * t477;
t456 = pkin(4) * t517 + qJ(5) * t516;
t457 = rSges(6,1) * t517 + rSges(6,3) * t516;
t539 = t548 * rSges(6,2);
t543 = t548 ^ 2;
t760 = pkin(3) * t522;
t478 = t518 + t760;
t625 = -qJ(4) + t549;
t513 = t546 * t625;
t633 = -t546 * t478 - t548 * t625;
t661 = -t546 * (t546 * t518 + t633 + t679) + t548 * (-t513 + t519 + (t478 - t518) * t548);
t752 = rSges(6,2) * t546;
t860 = t456 + t457;
t130 = t548 * (t457 * t548 + t752) + t543 * t456 + t661 + (t546 * t860 - t539) * t546;
t101 = t130 + t655;
t542 = t546 ^ 2;
t627 = t542 + t543;
t761 = pkin(3) * t521;
t593 = t627 * t761;
t786 = rSges(6,1) + pkin(4);
t631 = t786 * t719;
t712 = t517 * t548;
t745 = rSges(6,3) + qJ(5);
t632 = t745 * t712;
t718 = t516 * t548;
t594 = (-t786 * t718 + t632) * t548 + (t745 * t713 - t631) * t546;
t166 = -t593 + t594;
t589 = t761 + t762;
t612 = t745 * t517;
t617 = t786 * t516;
t834 = t617 - t612;
t562 = t589 + t834;
t258 = t562 * t546;
t260 = t562 * t548;
t590 = -t760 - t860;
t270 = t590 * t546;
t272 = t590 * t548;
t671 = -t258 * t270 - t260 * t272;
t46 = t101 * t166 + t671;
t354 = rSges(5,1) * t713 - rSges(5,2) * t719 - t548 * rSges(5,3);
t609 = -rSges(5,2) * t718 + rSges(5,3) * t546;
t147 = t546 * t354 + t548 * (rSges(5,1) * t712 + t609) + t661;
t122 = t147 + t655;
t455 = rSges(5,1) * t516 + rSges(5,2) * t517;
t630 = rSges(5,1) * t719 + rSges(5,2) * t713;
t643 = -t543 * t455 - t546 * t630;
t233 = -t593 + t643;
t566 = t455 + t589;
t290 = t566 * t546;
t292 = t566 * t548;
t753 = rSges(5,1) * t517;
t458 = -rSges(5,2) * t516 + t753;
t613 = -t458 - t760;
t326 = t613 * t546;
t328 = t613 * t548;
t667 = -t290 * t326 - t292 * t328;
t68 = t122 * t233 + t667;
t878 = m(4) * t104 + m(5) * t68 + m(6) * t46;
t591 = t761 + t834;
t269 = t591 * t546;
t264 = t546 * t269;
t271 = t591 * t548;
t614 = t455 + t761;
t325 = t614 * t546;
t817 = m(6) / 0.2e1;
t818 = m(5) / 0.2e1;
t586 = t614 * t548;
t837 = t548 * t586;
t672 = (-t271 * t548 - t264) * t817 + (-t546 * t325 - t837) * t818;
t266 = (-t612 + t761) * t546 + t631;
t267 = (-t617 - t761) * t548 + t632;
t309 = pkin(3) * t709 + t630;
t674 = (t266 * t546 - t267 * t548) * t817 + (t309 * t546 + t837) * t818;
t53 = t674 - t672;
t877 = qJD(1) * t53;
t247 = t546 * t258;
t676 = (-t260 * t548 - t247) * t817 + (-t546 * t290 - t292 * t548) * t818;
t470 = t548 * t589;
t256 = -t548 * t617 - t470 + t632;
t257 = (t589 - t612) * t546 + t631;
t294 = t546 * t589 + t630;
t295 = -t455 * t548 - t470;
t678 = (-t256 * t548 + t257 * t546) * t817 + (t294 * t546 - t295 * t548) * t818;
t40 = t678 - t676;
t876 = t40 * qJD(1);
t448 = -Icges(5,2) * t516 + t507;
t874 = t448 + t880;
t872 = -t882 * t548 - t879;
t871 = -t873 * t548 + t881;
t345 = -Icges(6,4) * t548 + t546 * t832;
t484 = Icges(5,4) * t719;
t347 = Icges(5,1) * t713 - Icges(5,5) * t548 - t484;
t870 = -Icges(5,2) * t713 - t574 * t546 + t345 + t347 - t484;
t483 = Icges(6,5) * t712;
t338 = Icges(6,6) * t546 + Icges(6,3) * t718 + t483;
t344 = Icges(5,6) * t546 + t448 * t548;
t869 = -Icges(6,1) * t718 - t451 * t548 + t338 - t344 + t483;
t444 = Icges(6,3) * t516 + t737;
t337 = -Icges(6,6) * t548 + t444 * t546;
t343 = Icges(5,4) * t713 - Icges(5,2) * t719 - Icges(5,6) * t548;
t868 = t880 * t546 - t337 + t343;
t503 = Icges(4,4) * t709;
t399 = Icges(4,1) * t705 - Icges(4,5) * t548 - t503;
t867 = -t347 * t712 - t399 * t704 - t875 * t546;
t866 = t832 + t452;
t446 = Icges(6,4) * t517 + Icges(6,6) * t516;
t727 = t446 * t548;
t865 = t338 * t718 + t400 * t704 + t881 * t712 + (Icges(6,2) * t546 + t727 + t882) * t546;
t864 = Icges(4,5) * t521 + Icges(4,6) * t522 + (Icges(5,6) - Icges(6,6)) * t517 + (Icges(6,4) + Icges(5,5)) * t516;
t515 = Icges(4,4) * t522;
t473 = -Icges(4,2) * t521 + t515;
t398 = Icges(4,6) * t546 + t473 * t548;
t863 = -t344 * t719 - t398 * t709 + t872;
t397 = Icges(4,4) * t705 - Icges(4,2) * t709 - Icges(4,6) * t548;
t862 = t343 * t718 + t397 * t708 + t867;
t474 = Icges(4,1) * t521 + t515;
t859 = t473 + t474;
t858 = t344 * t516 + t398 * t521 - t875;
t857 = t343 * t516 + t397 * t521;
t856 = -t344 * t718 - t398 * t708 + t865;
t838 = (t337 * t516 + t345 * t517) * t546;
t855 = t838 + t865;
t819 = m(4) / 0.2e1;
t789 = t546 / 0.2e1;
t787 = -t548 / 0.2e1;
t850 = t548 / 0.2e1;
t728 = t446 * t546;
t315 = t546 * (-Icges(6,2) * t548 + t728);
t189 = t337 * t718 + t345 * t712 + t315;
t849 = t189 * t548;
t740 = Icges(3,4) * t545;
t493 = Icges(3,2) * t547 + t740;
t496 = Icges(3,1) * t547 - t740;
t848 = (t496 / 0.2e1 - t493 / 0.2e1) * t545;
t611 = t518 + t754;
t287 = -t546 * t611 + t629 - t679;
t288 = t548 * t611 - t519 + t610;
t560 = (-t287 * t548 - t288 * t546) * t477;
t273 = -t354 + t633;
t274 = -t513 + (t478 + t753) * t548 + t609;
t668 = t328 * t273 + t326 * t274;
t829 = -t745 * t516 - t786 * t517;
t231 = t546 * t829 + t539 + t633;
t232 = t752 - t513 + (t478 - t829) * t548;
t675 = t272 * t231 + t270 * t232;
t622 = (t560 + (t546 * t835 - t548 * t836) * t476) * t819 + (-t256 * t269 - t257 * t271 + t675) * t817 + (-t294 * t586 - t295 * t325 + t668) * t818;
t623 = (-t441 * t835 + t442 * t836 + t560) * t819 + (-t258 * t267 - t260 * t266 + t675) * t817 + (t290 * t586 - t292 * t309 + t668) * t818;
t847 = t622 - t623;
t846 = t864 * t546;
t845 = t864 * t548;
t472 = Icges(4,2) * t522 + t739;
t844 = (-t472 + t475) * t522 - t859 * t521 + (t866 - t873) * t517 + (t444 - t874) * t516;
t639 = -t472 * t548 + t400;
t640 = -Icges(4,2) * t705 + t399 - t503;
t641 = -t474 * t548 - t398;
t642 = t474 * t546 + t397;
t842 = (-t639 * t546 + t548 * t640) * t521 + (t641 * t546 + t548 * t642) * t522 + (t869 * t546 + t868 * t548) * t517 + (-t871 * t546 + t870 * t548) * t516;
t535 = Icges(3,4) * t547;
t494 = -Icges(3,2) * t545 + t535;
t495 = Icges(3,1) * t545 + t535;
t432 = Icges(3,5) * t546 + t496 * t548;
t635 = -t493 * t548 + t432;
t702 = t545 * t546;
t510 = Icges(3,4) * t702;
t694 = t546 * t547;
t431 = Icges(3,1) * t694 - Icges(3,5) * t548 - t510;
t636 = -Icges(3,2) * t694 + t431 - t510;
t430 = Icges(3,6) * t546 + t494 * t548;
t637 = -t495 * t548 - t430;
t429 = Icges(3,4) * t694 - Icges(3,2) * t702 - Icges(3,6) * t548;
t638 = t495 * t546 + t429;
t828 = (-t635 * t546 + t548 * t636) * t545 + (t637 * t546 + t548 * t638) * t547;
t554 = (t856 * t546 + t862 * t548 - t849) * t850 + (-t849 + ((t857 + t882) * t548 + t863 + t867 + t879) * t548 + (-t838 + t855) * t546) * t787 + ((t546 * t858 + t189 - t315 - t862 + t863 - t872) * t546 + ((t858 + t875) * t548 + (-t347 * t517 - t399 * t522 + t857) * t546 - t855 + t856) * t548) * t789;
t553 = t859 * t522 / 0.2e1 + (-t472 / 0.2e1 + t475 / 0.2e1) * t521 + (-t444 / 0.2e1 + t874 / 0.2e1) * t517 + (-t873 / 0.2e1 + t866 / 0.2e1) * t516;
t824 = 0.4e1 * qJD(1);
t823 = 2 * qJD(2);
t821 = 2 * qJD(3);
t820 = 4 * qJD(3);
t701 = t545 * t548;
t644 = t542 * (-t589 + t762) + t548 * (pkin(2) * t701 - t470);
t193 = t643 + t644;
t665 = -t325 * t326 - t328 * t586;
t810 = m(5) * (t147 * t193 + t665);
t433 = t627 * t516;
t666 = -t517 * t247 - t260 * t712;
t67 = t101 * t433 + t666;
t663 = -t517 * t264 - t271 * t712;
t77 = t130 * t433 + t663;
t806 = m(6) * (t77 + t67);
t619 = t516 * t101 + t666;
t664 = t270 * t719 + t272 * t718;
t732 = t166 * t517;
t804 = m(6) * (t619 + t664 - t732);
t145 = t594 + t644;
t585 = t516 * t130 + t663 + t664;
t803 = m(6) * (-t145 * t517 + t585);
t670 = -t269 * t270 - t271 * t272;
t800 = m(6) * (t130 * t145 + t670);
t669 = t231 * t712 + t232 * t713;
t795 = m(6) * ((t256 * t546 + t257 * t548) * t516 + t669);
t794 = m(6) * (-t258 * t718 + t260 * t719 + t669);
t793 = m(6) * ((t266 * t548 + t267 * t546) * t516 + t669);
t792 = m(6) * (-t269 * t718 + t271 * t719 + t669);
t791 = m(6) * (t231 * t257 + t232 * t256);
t790 = -t546 / 0.2e1;
t755 = rSges(3,1) * t547;
t616 = pkin(1) + t755;
t628 = rSges(3,2) * t702 + t548 * rSges(3,3);
t329 = -t546 * t616 + t541 + t628;
t512 = rSges(3,2) * t701;
t330 = -t512 + t616 * t548 + (rSges(3,3) + pkin(6)) * t546;
t497 = rSges(3,1) * t545 + rSges(3,2) * t547;
t468 = t497 * t546;
t469 = t497 * t548;
t785 = m(3) * (t329 * t468 - t330 * t469);
t140 = t627 * t476 * t477 + t262 * t281;
t138 = m(4) * t140;
t781 = m(4) * (t287 * t836 - t288 * t835);
t780 = m(4) * (t287 * t441 - t288 * t442);
t779 = m(5) * (t273 * t294 + t274 * t295);
t778 = m(5) * (t273 * t309 - t274 * t586);
t777 = m(5) * (t273 * t548 + t274 * t546);
t770 = m(6) * (t231 * t266 + t232 * t267);
t769 = m(6) * (t231 * t548 + t232 * t546);
t758 = m(6) * qJD(2);
t757 = m(6) * qJD(3);
t756 = m(6) * qJD(5);
t720 = t516 * t517;
t703 = t545 * t429;
t693 = t547 * t548;
t427 = Icges(3,5) * t694 - Icges(3,6) * t702 - Icges(3,3) * t548;
t646 = -t546 * t427 - t431 * t693;
t578 = Icges(3,5) * t547 - Icges(3,6) * t545;
t428 = Icges(3,3) * t546 + t548 * t578;
t645 = t546 * t428 + t432 * t693;
t634 = t627 * t720;
t275 = m(6) * t433;
t626 = t275 * qJD(1);
t141 = -t231 * t719 + t232 * t718;
t621 = m(6) * t141 * qJD(1);
t620 = t130 * t166 + t670;
t618 = t147 * t233 + t665;
t615 = -t477 - t540;
t355 = t432 * t694;
t597 = t548 * t428 - t355;
t595 = t545 * t430 - t427;
t592 = t627 * t762;
t588 = -t540 - t760;
t584 = ((t846 * t546 + t842) * t548 - t845 * t542) * t789 + ((t845 * t548 + t842) * t546 - t846 * t543) * t787;
t577 = -Icges(3,5) * t545 - Icges(3,6) * t547;
t565 = -t458 + t588;
t563 = t138 + t584;
t561 = t588 - t860;
t559 = -t592 + t644;
t551 = -t554 + (t869 * t516 + t871 * t517 + t521 * t641 + t522 * t639 + t885 * t546 + t844 * t548 + t728) * t789 + (-t868 * t516 + t870 * t517 - t521 * t642 + t522 * t640 + t844 * t546 - t727 - t883) * t787;
t550 = -t553 + (t789 + t790) * (t397 * t522 + t399 * t521 + (t337 + t343) * t517 + (-t345 + t347) * t516);
t499 = -rSges(3,2) * t545 + t755;
t463 = t577 * t548;
t462 = t577 * t546;
t393 = t615 * t548;
t391 = t615 * t546;
t293 = t565 * t548;
t291 = t565 * t546;
t280 = t634 - t720;
t268 = t280 * t756;
t261 = t561 * t548;
t259 = t561 * t546;
t255 = -t592 + t281;
t220 = -t430 * t701 + t645;
t219 = -t429 * t701 - t646;
t218 = -t430 * t702 - t597;
t170 = (-t433 * t517 - t280 + t634) * t756;
t168 = t559 + t643;
t143 = -t219 * t548 + t220 * t546;
t142 = -(-t546 * (-t547 * t431 + t703) - t548 * t427) * t548 + t218 * t546;
t137 = t559 + t594;
t135 = m(6) * (-t270 * t548 + t272 * t546) + m(5) * (-t326 * t548 + t328 * t546);
t134 = t135 * qJD(3);
t98 = t769 + t777;
t94 = t792 / 0.2e1;
t92 = t793 / 0.2e1;
t85 = t794 / 0.2e1;
t83 = t795 / 0.2e1;
t59 = (t218 - t355 + (t428 + t703) * t548 + t646) * t548 + t645 * t546;
t58 = (t548 * t595 + t220 - t645) * t548 + (t546 * t595 + t219 + t597) * t546;
t54 = t672 + t674;
t43 = t803 / 0.2e1;
t42 = t676 + t678;
t34 = t804 / 0.2e1;
t23 = t806 / 0.2e1;
t22 = t553 + t770 + t778 + t780;
t21 = t94 - t793 / 0.2e1;
t20 = t94 + t92;
t19 = t92 - t792 / 0.2e1;
t16 = t85 - t795 / 0.2e1;
t15 = t85 + t83;
t14 = t83 - t794 / 0.2e1;
t13 = t785 + t781 + t779 + t791 + (t495 / 0.2e1 + t494 / 0.2e1) * t547 + t848 + t553;
t10 = t23 + t43 - t804 / 0.2e1;
t9 = t23 + t34 - t803 / 0.2e1;
t8 = t34 + t43 - t806 / 0.2e1;
t7 = t563 + t800 + t810;
t6 = t584 + t878;
t4 = t554 + (t58 / 0.2e1 + t142 / 0.2e1) * t546 + (t143 / 0.2e1 - t59 / 0.2e1) * t548;
t3 = t554 + t847;
t2 = t554 - t847;
t1 = t551 + t622 + t623;
t5 = [t13 * qJD(2) + t22 * qJD(3) + t98 * qJD(4) + t141 * t756, t13 * qJD(1) + t1 * qJD(3) + t42 * qJD(4) + t15 * qJD(5) + (m(3) * ((-t329 * t548 - t330 * t546) * t499 + (-t468 * t548 + t469 * t546) * t497) / 0.2e1 + (t287 * t393 + t288 * t391) * t819 + (t273 * t293 + t274 * t291 - t290 * t295 - t292 * t294) * t818 + (t231 * t261 + t232 * t259 - t256 * t258 - t257 * t260) * t817) * t823 + (t551 + t59 * t850 + (t545 * t637 + t547 * t635) * t789 + (t58 + t142) * t790 + (-t545 * t638 + t547 * t636 + t143) * t787 + (t543 / 0.2e1 + t542 / 0.2e1) * t578) * qJD(2), t22 * qJD(1) + t1 * qJD(2) + t551 * qJD(3) + t54 * qJD(4) + t20 * qJD(5) + ((t560 + (-t441 * t548 + t442 * t546) * t476) * t819 + (t668 + (-t309 + t325) * t586) * t818 + (-t266 * t271 - t267 * t269 + t675) * t817) * t821, qJD(1) * t98 + qJD(2) * t42 + qJD(3) * t54, t15 * qJD(2) + t20 * qJD(3) + t621; (t550 - (t495 + t494) * t547 / 0.2e1 - t848) * qJD(1) + t4 * qJD(2) + t2 * qJD(3) - t40 * qJD(4) + t16 * qJD(5) + (-t785 / 0.4e1 - t781 / 0.4e1 - t779 / 0.4e1 - t791 / 0.4e1) * t824, t4 * qJD(1) + (m(5) * (t122 * t168 - t290 * t291 - t292 * t293) + m(4) * (t178 * t255 - t391 * t836 - t393 * t835) + (t542 * t463 + (-t546 * t462 + t828) * t548) * t789 + (t543 * t462 + (-t548 * t463 + t828) * t546) * t787 + m(3) * ((t546 * (rSges(3,1) * t694 - t628) + t548 * (rSges(3,1) * t693 + rSges(3,3) * t546 - t512)) * (-t468 * t546 - t469 * t548) + t627 * t499 * t497) + m(6) * (t101 * t137 - t258 * t259 - t260 * t261) + t584) * qJD(2) + t6 * qJD(3) + t67 * t756, t2 * qJD(1) + t6 * qJD(2) + t584 * qJD(3) + t9 * qJD(5) + (-t810 / 0.4e1 - t800 / 0.4e1 - t138 / 0.4e1) * t820 + ((t104 + t140) * t819 + (t620 + t46) * t817 + (t618 + t68) * t818) * t821, -t876, t16 * qJD(1) + t9 * qJD(3) + t67 * t758 + t170; t550 * qJD(1) + t3 * qJD(2) + t554 * qJD(3) - t53 * qJD(4) + t21 * qJD(5) + (-t780 / 0.4e1 - t778 / 0.4e1 - t770 / 0.4e1) * t824, t3 * qJD(1) + t7 * qJD(3) + t10 * qJD(5) + ((t122 * t193 + t147 * t168 - t291 * t325 - t293 * t586 + t667) * t818 + (t255 * t262 + (-t391 * t546 - t393 * t548) * t476 + t104) * t819 + (t101 * t145 + t130 * t137 - t259 * t269 - t261 * t271 + t671) * t817) * t823 + (t584 - t878) * qJD(2), t554 * qJD(1) + t7 * qJD(2) + t563 * qJD(3) + (m(5) * t618 / 0.4e1 + m(6) * t620 / 0.4e1) * t820 + t77 * t756, -t877, t21 * qJD(1) + t10 * qJD(2) + t757 * t77 + t170; t40 * qJD(2) + t53 * qJD(3) - t275 * qJD(5) + (-t777 / 0.4e1 - t769 / 0.4e1) * t824, t876 + ((-t259 * t548 + t261 * t546) * t817 + (-t291 * t548 + t293 * t546) * t818) * t823 + t134, qJD(2) * t135 + t134 + t877, 0, -t626; t14 * qJD(2) + t19 * qJD(3) + t275 * qJD(4) - t621, t14 * qJD(1) + (-t137 * t517 + (t259 * t546 + t261 * t548) * t516 - t67 + t619) * t758 + t8 * qJD(3) + t268, t19 * qJD(1) + t8 * qJD(2) + (t585 - t77 - t732) * t757 + t268, t626, 0.4e1 * (qJD(2) / 0.4e1 + qJD(3) / 0.4e1) * t280 * m(6);];
Cq = t5;
