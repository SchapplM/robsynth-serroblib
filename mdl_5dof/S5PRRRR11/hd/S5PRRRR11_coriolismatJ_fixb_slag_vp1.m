% Calculate matrix of centrifugal and coriolis load on the joints for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% m [6x1]
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
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5PRRRR11_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_coriolismatJ_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR11_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR11_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:14
% EndTime: 2024-09-27 21:45:32
% DurationCPUTime: 9.61s
% Computational Cost: add. (227128->697), mult. (126653->945), div. (0->0), fcn. (119634->20), ass. (0->458)
t867 = -qJD(2) / 0.4e1;
t866 = m(6) * t867;
t526 = pkin(10) + qJ(2);
t523 = cos(t526);
t749 = cos(pkin(5));
t765 = cos(qJ(3));
t598 = t749 * t765;
t577 = t523 * t598;
t522 = sin(t526);
t763 = sin(qJ(3));
t647 = t522 * t763;
t486 = t647 - t577;
t597 = t749 * t763;
t648 = t522 * t765;
t487 = t523 * t597 + t648;
t527 = sin(pkin(5));
t713 = t523 * t527;
t512 = rSges(4,3) * t713;
t402 = t487 * rSges(4,1) - t486 * rSges(4,2) - t512;
t513 = pkin(7) * t713;
t360 = -t522 * pkin(2) - t402 + t513;
t646 = t523 * t765;
t489 = -t522 * t597 + t646;
t548 = t522 * t598 + t523 * t763;
t590 = t489 * rSges(4,1) - rSges(4,2) * t548;
t715 = t522 * t527;
t362 = t523 * pkin(2) + (rSges(4,3) + pkin(7)) * t715 + t590;
t417 = -rSges(4,1) * t486 - rSges(4,2) * t487;
t418 = -rSges(4,1) * t548 - rSges(4,2) * t489;
t709 = qJ(3) + qJ(4);
t633 = pkin(5) + t709;
t574 = sin(t633) / 0.2e1;
t634 = pkin(5) - t709;
t599 = sin(t634);
t502 = t574 + t599 / 0.2e1;
t575 = cos(t633) / 0.2e1;
t600 = cos(t634);
t504 = t575 - t600 / 0.2e1;
t747 = Icges(5,4) * t504;
t429 = Icges(5,2) * t502 + Icges(5,6) * t749 - t747;
t452 = Icges(5,1) * t502 + t747;
t603 = -qJ(5) + t634;
t579 = sin(t603);
t516 = t579 / 0.2e1;
t602 = qJ(5) + t633;
t578 = sin(t602);
t565 = t578 / 0.2e1;
t494 = t565 + t516;
t566 = cos(t602) / 0.2e1;
t580 = cos(t603);
t495 = t566 - t580 / 0.2e1;
t746 = Icges(6,4) * t495;
t422 = Icges(6,2) * t494 + Icges(6,6) * t749 - t746;
t433 = Icges(6,1) * t494 + t746;
t431 = Icges(6,5) * t494 + Icges(6,6) * t495;
t625 = t749 * t431;
t491 = Icges(6,4) * t494;
t423 = -Icges(6,1) * t495 + Icges(6,5) * t749 + t491;
t680 = Icges(6,2) * t495 + t423 + t491;
t854 = t680 * t494;
t588 = t625 / 0.2e1 + t854 / 0.2e1 + (t422 / 0.2e1 - t433 / 0.2e1) * t495;
t450 = Icges(5,5) * t502 + Icges(5,6) * t504;
t624 = t749 * t450;
t497 = Icges(5,4) * t502;
t430 = -Icges(5,1) * t504 + Icges(5,5) * t749 + t497;
t678 = Icges(5,2) * t504 + t430 + t497;
t855 = t678 * t502;
t547 = t588 + t624 / 0.2e1 + t855 / 0.2e1 + (-t452 / 0.2e1 + t429 / 0.2e1) * t504;
t498 = (Icges(4,5) * t765 - Icges(4,6) * t763) * t527;
t622 = t749 * t498;
t864 = -t547 - t622 / 0.2e1 - m(4) * (-t360 * t417 + t362 * t418);
t863 = m(5) * t867;
t425 = -t495 * rSges(6,1) + t494 * rSges(6,2) + rSges(6,3) * t749;
t760 = pkin(4) * t527;
t762 = sin(qJ(4));
t764 = cos(qJ(4));
t679 = t749 * pkin(9) + (t765 * t762 + t763 * t764) * t760 + t425;
t583 = pkin(3) * t597;
t783 = pkin(7) + pkin(8);
t550 = -t527 * t783 + t583;
t545 = t523 * t550;
t664 = t765 * pkin(3);
t606 = t664 + pkin(2);
t536 = t522 * t606 + t545;
t560 = (pkin(4) * t764 + pkin(3)) * t749;
t596 = t749 * t762;
t529 = t763 * t560 + t765 * pkin(4) * t596 - t527 * (pkin(9) + t783);
t635 = cos(t709);
t552 = pkin(4) * t635 + t606;
t676 = -t522 * t552 - t523 * t529;
t301 = t536 + t676;
t541 = t580 / 0.2e1 + t566;
t654 = qJ(5) + t709;
t604 = sin(t654);
t446 = t522 * t604 - t523 * t541;
t540 = t565 - t579 / 0.2e1;
t605 = cos(t654);
t573 = t522 * t605;
t447 = t523 * t540 + t573;
t315 = t447 * rSges(6,1) - t446 * rSges(6,2) - rSges(6,3) * t713;
t853 = t749 * t315;
t860 = t749 * t301 - t853;
t183 = -t679 * t713 + t860;
t410 = -pkin(3) * t648 - t513 - t545;
t404 = t749 * t410;
t645 = t527 * t763;
t492 = pkin(3) * t645 + pkin(8) * t749;
t649 = t492 + t679;
t148 = -t649 * t713 + t404 + t860;
t530 = -t522 * t541 - t523 * t604;
t572 = t523 * t605;
t531 = -t522 * t540 + t572;
t308 = Icges(6,5) * t531 + Icges(6,6) * t530 + Icges(6,3) * t715;
t862 = t308 * t713;
t544 = t600 / 0.2e1 + t575;
t524 = sin(t709);
t714 = t523 * t524;
t465 = -t522 * t544 - t714;
t503 = t574 - t599 / 0.2e1;
t592 = t523 * t635;
t558 = -t522 * t503 + t592;
t326 = Icges(5,5) * t558 + Icges(5,6) * t465 + Icges(5,3) * t715;
t861 = t326 * t713;
t716 = t522 * t524;
t462 = -t523 * t544 + t716;
t463 = t523 * t503 + t522 * t635;
t325 = -Icges(5,5) * t463 + Icges(5,6) * t462 + Icges(5,3) * t713;
t859 = t325 * t713;
t307 = -Icges(6,5) * t447 + Icges(6,6) * t446 + Icges(6,3) * t713;
t858 = t307 * t713;
t435 = -t504 * rSges(5,1) + t502 * rSges(5,2) + rSges(5,3) * t749;
t335 = t463 * rSges(5,1) - t462 * rSges(5,2) - rSges(5,3) * t713;
t852 = t749 * t335;
t251 = -t435 * t713 - t852;
t677 = t435 + t492;
t216 = -t677 * t713 + t404 - t852;
t243 = -t425 * t713 - t853;
t245 = -t315 + t676;
t836 = Icges(5,4) * t463;
t328 = Icges(5,2) * t462 + Icges(5,6) * t713 - t836;
t851 = -Icges(5,1) * t462 + t328 - t836;
t835 = Icges(5,4) * t558;
t329 = Icges(5,2) * t465 + Icges(5,6) * t715 + t835;
t850 = -Icges(5,1) * t465 + t329 + t835;
t849 = t429 - t452;
t390 = -Icges(4,5) * t487 + Icges(4,6) * t486 + Icges(4,3) * t713;
t848 = t390 * t713;
t310 = -Icges(6,4) * t447 + Icges(6,2) * t446 + Icges(6,6) * t713;
t438 = Icges(6,4) * t446;
t313 = -Icges(6,1) * t447 + Icges(6,5) * t713 + t438;
t847 = t310 * t530 + t313 * t531;
t454 = Icges(5,4) * t462;
t331 = -Icges(5,1) * t463 + Icges(5,5) * t713 + t454;
t846 = t465 * t328 + t331 * t558;
t485 = t749 * rSges(4,3) + (rSges(4,1) * t763 + rSges(4,2) * t765) * t527;
t293 = -t402 * t749 - t485 * t713;
t474 = Icges(4,4) * t487;
t392 = -Icges(4,2) * t486 - Icges(4,6) * t713 + t474;
t473 = Icges(4,4) * t486;
t396 = -Icges(4,1) * t487 + Icges(4,5) * t713 + t473;
t845 = t392 * t548 + t489 * t396;
t843 = -t851 * t713 - t850 * t715;
t842 = m(5) / 0.2e1;
t786 = m(5) / 0.4e1;
t841 = m(6) / 0.2e1;
t784 = m(6) / 0.4e1;
t159 = -t307 * t715 - t847;
t732 = t308 * t523;
t733 = t307 * t522;
t821 = t527 * t523 ^ 2;
t838 = (t307 * t821 - t858 * t523 + (t159 + t847 + (t732 + t733) * t527 - t862) * t522) * t527;
t171 = -t325 * t715 - t846;
t730 = t326 * t523;
t731 = t325 * t522;
t837 = (t325 * t821 - t859 * t523 + (t171 + t846 + (t730 + t731) * t527 - t861) * t522) * t527;
t434 = rSges(6,1) * t494 + rSges(6,2) * t495;
t801 = (-t434 - (-t762 * t763 + t764 * t765) * t760) * t527;
t823 = t522 * t801;
t822 = t523 * t801;
t391 = Icges(4,5) * t489 - Icges(4,6) * t548 + Icges(4,3) * t715;
t817 = t391 * t713;
t528 = rSges(6,1) * t531 + rSges(6,2) * t530 + rSges(6,3) * t715;
t298 = t749 * t528;
t444 = t522 * t529;
t534 = pkin(4) * t592 + t522 * t550 - t444;
t696 = -t749 * t534 - t298;
t181 = t679 * t715 + t696;
t671 = t516 - t578 / 0.2e1;
t449 = -t522 * t671 - t572;
t356 = t530 * rSges(6,1) + t449 * rSges(6,2);
t323 = t749 * t356;
t477 = t548 * pkin(3);
t562 = t763 * t596;
t471 = -pkin(4) * t562 + t560 * t765;
t509 = -pkin(3) * t763 - pkin(4) * t524;
t612 = -t522 * t471 + t523 * t509;
t365 = t612 + t477;
t692 = t749 * t365 + t323;
t199 = t692 + t823;
t514 = pkin(3) * t647;
t717 = t522 * t509;
t364 = t717 + t514 + (-pkin(3) * t598 + t471) * t523;
t448 = t523 * t671 - t573;
t589 = rSges(6,1) * t446 - rSges(6,2) * t448;
t627 = t749 * t589;
t556 = -t364 * t749 + t627;
t200 = t556 + t822;
t247 = t523 * t552 - t444 + t528;
t750 = t465 * rSges(5,2);
t538 = rSges(5,1) * t558 + rSges(5,3) * t715 + t750;
t322 = t749 * t538;
t249 = t435 * t715 - t322;
t262 = -t523 * t471 + t589 - t717;
t263 = t612 + t356;
t372 = -t462 * rSges(5,1) - rSges(5,2) * t463;
t476 = pkin(3) * t577 - t514;
t295 = -t372 - t476;
t373 = t465 * rSges(5,1) - rSges(5,2) * t558;
t296 = -t477 + t373;
t363 = t749 * t373;
t453 = rSges(5,1) * t502 + rSges(5,2) * t504;
t267 = -t453 * t715 + t363;
t279 = t750 + (rSges(5,1) * t635 + t606) * t523 + (-t583 - t503 * rSges(5,1) + (rSges(5,3) + t783) * t527) * t522;
t202 = 0.2e1 * t267 * t279;
t626 = t749 * t372;
t268 = -t453 * t713 - t626;
t277 = -t536 - t335;
t203 = 0.2e1 * t268 * t277;
t705 = t202 + t203;
t757 = (-t181 * t263 + t183 * t262 + t199 * t247 + t200 * t245) * t841 + (-0.2e1 * t249 * t296 + 0.2e1 * t251 * t295 + t705) * t786;
t539 = pkin(3) * t646 - t522 * (-t527 * pkin(8) + t583);
t405 = t749 * t539;
t146 = t649 * t715 - t405 + t696;
t484 = (t598 * t764 - t562) * pkin(4);
t437 = -pkin(4) * t714 - t484 * t522;
t219 = t437 * t749 + t323 + t823;
t793 = 0.2e1 * t247;
t155 = t219 * t793;
t436 = -pkin(4) * t716 + t484 * t523;
t220 = -t436 * t749 + t627 + t822;
t794 = 0.2e1 * t245;
t156 = t220 * t794;
t214 = t677 * t715 - t322 - t405;
t270 = t437 + t356;
t791 = -0.2e1 * t270;
t269 = -t436 + t589;
t792 = 0.2e1 * t269;
t758 = (t146 * t791 + t148 * t792 + t155 + t156) * t784 + (-0.2e1 * t214 * t373 - 0.2e1 * t216 * t372 + t705) * t786;
t816 = t757 - t758;
t813 = 0.4e1 * m(6);
t642 = t763 * Icges(4,4);
t481 = Icges(4,6) * t749 + (Icges(4,2) * t765 + t642) * t527;
t521 = Icges(4,4) * t527 * t765;
t482 = Icges(4,1) * t645 + Icges(4,5) * t749 + t521;
t266 = (Icges(4,3) * t749 + (Icges(4,5) * t763 + Icges(4,6) * t765) * t527) * t713 + t486 * t481 - t487 * t482;
t811 = -t266 / 0.2e1;
t810 = t527 / 0.2e1;
t657 = t749 / 0.2e1;
t694 = Icges(6,2) * t448 - t313 - t438;
t439 = Icges(6,4) * t530;
t314 = Icges(6,1) * t531 + Icges(6,5) * t715 + t439;
t693 = Icges(6,2) * t449 + t314 + t439;
t237 = t356 * t713 - t589 * t715;
t239 = t372 * t715 + t373 * t713;
t684 = -Icges(4,1) * t486 - t392 - t474;
t748 = Icges(4,4) * t489;
t394 = -Icges(4,2) * t548 + Icges(4,6) * t715 + t748;
t683 = -Icges(4,1) * t548 - t394 - t748;
t682 = -Icges(4,2) * t487 - t396 - t473;
t475 = Icges(4,4) * t548;
t397 = Icges(4,1) * t489 + Icges(4,5) * t715 - t475;
t681 = -Icges(4,2) * t489 + t397 - t475;
t674 = t476 * t715 - t477 * t713;
t225 = t315 * t715 + t528 * t713;
t242 = t425 * t715 - t298;
t253 = -t434 * t715 + t323;
t254 = -t434 * t713 + t627;
t804 = t225 * t237 - t242 * t253 + t243 * t254;
t233 = t335 * t715 + t538 * t713;
t685 = -t410 * t715 + t539 * t713;
t178 = t685 + t233;
t591 = t178 * t239 - t214 * t267 + t216 * t268;
t142 = -t301 * t715 + t534 * t713 + t225;
t803 = t142 * t237 - t181 * t253 + t183 * t254;
t798 = 0.2e1 * t225;
t796 = -0.2e1 * t242;
t795 = 0.2e1 * t243;
t790 = 0.2e1 * t589;
t789 = -0.2e1 * t356;
t652 = 0.2e1 * t803;
t121 = t142 + t685;
t745 = t121 * t237;
t104 = 0.2e1 * t745;
t743 = t146 * t253;
t123 = -0.2e1 * t743;
t742 = t148 * t254;
t124 = 0.2e1 * t742;
t653 = t104 + t123 + t124;
t777 = m(6) * (t652 + t653);
t726 = t437 * t523;
t727 = t436 * t522;
t563 = (t726 + t727) * t527;
t197 = t563 + t237;
t651 = t197 * t798 + t219 * t796 + t220 * t795;
t775 = m(6) * (t651 + t653);
t343 = t364 * t715;
t344 = t365 * t713;
t166 = t343 + t344 + t237;
t774 = m(6) * (t166 * t798 + t199 * t796 + t200 * t795 + t652);
t194 = t253 * t793;
t195 = t254 * t794;
t706 = t194 + t195;
t771 = m(6) * (t146 * t789 + t148 * t790 + t706);
t770 = t803 * t813;
t769 = m(6) * (t181 * t789 + t183 * t790 + t706);
t768 = m(6) * (t262 * t795 + t263 * t796 + t706);
t767 = m(6) * (t242 * t791 + t243 * t792 + t706);
t766 = t804 * t813;
t236 = 0.2e1 * t237;
t761 = m(6) * t236;
t755 = m(5) * qJD(3);
t754 = m(5) * qJD(4);
t752 = m(6) * qJD(3);
t751 = m(6) * qJD(4);
t729 = t390 * t522;
t728 = t391 * t523;
t349 = -Icges(6,5) * t446 + Icges(6,6) * t448;
t712 = t527 * t349;
t711 = t527 * t431;
t691 = -Icges(5,2) * t463 - t331 - t454;
t455 = Icges(5,4) * t465;
t332 = Icges(5,1) * t558 + Icges(5,5) * t715 + t455;
t690 = -Icges(5,2) * t558 + t332 + t455;
t688 = 0.2e1 * t239;
t675 = 0.2e1 * t674;
t500 = (Icges(4,1) * t765 - t642) * t527;
t673 = -t481 + t500;
t499 = -Icges(4,2) * t645 + t521;
t672 = t482 + t499;
t333 = 0.2e1 * t343;
t334 = 0.2e1 * t344;
t109 = (t333 / 0.4e1 + t334 / 0.4e1 + (-t727 / 0.2e1 - t726 / 0.2e1) * t527) * m(6);
t670 = t109 * qJD(1);
t353 = -Icges(6,1) * t446 + Icges(6,4) * t448;
t127 = t749 * t349 + (-t310 - t353) * t495 + t694 * t494;
t311 = Icges(6,4) * t531 + Icges(6,2) * t530 + Icges(6,6) * t715;
t350 = Icges(6,5) * t530 + Icges(6,6) * t449;
t354 = Icges(6,1) * t530 + Icges(6,4) * t449;
t128 = t749 * t350 + (t311 - t354) * t495 + t693 * t494;
t149 = t448 * t422 + t447 * t433 - t446 * t680 - t523 * t711;
t150 = t449 * t422 + t433 * t531 + t522 * t711 + t680 * t530;
t180 = (t625 + (t422 - t433) * t495 + t854) * t749;
t637 = -t713 / 0.2e1;
t638 = t715 / 0.2e1;
t665 = ((t449 * t311 + t350 * t715 + t354 * t531) * t715 - (-t310 * t449 + t353 * t531 + t522 * t712) * t713 + t150 * t749 + (t693 * t715 - t694 * t713) * t530) * t638 + ((t448 * t311 - t350 * t713 + t447 * t354 - t693 * t446) * t715 - (-t310 * t448 + t447 * t353 - t694 * t446 - t523 * t712) * t713 + t149 * t749) * t637 + (t180 + (-t127 * t523 + t128 * t522) * t527) * t657;
t208 = -t390 * t715 - t845;
t209 = t391 * t715 - t394 * t548 + t489 * t397;
t656 = (-t208 * t523 + t209 * t522) * t810 - ((t209 + t848) * t522 + ((t728 - t729) * t527 - t208 - t817) * t523) * t527 / 0.2e1;
t655 = t266 * t657 + t749 * t811 + (t390 * t821 - t848 * t523 + (t208 + (t728 + t729) * t527 - t817 + t845) * t522) * t810;
t650 = 0.2e1 * t804;
t160 = t308 * t715 + t530 * t311 + t531 * t314;
t172 = t326 * t715 + t465 * t329 + t558 * t332;
t644 = t765 * t482;
t643 = t765 * t499;
t641 = t763 * t481;
t640 = t763 * t500;
t639 = -t715 / 0.2e1;
t636 = t713 / 0.2e1;
t623 = t749 * t476;
t616 = t774 / 0.4e1 + t665;
t615 = -t766 / 0.4e1 + t665;
t611 = t527 ^ 2 * t664;
t610 = t333 + t334 + t236;
t210 = ((-Icges(6,5) * t495 + Icges(6,6) * t494 + Icges(6,3) * t749) * t715 + t422 * t530 + t423 * t531) * t749;
t42 = t210 + ((t160 + t858) * t522 + ((t732 - t733) * t527 - t159 - t862) * t523) * t527;
t88 = t210 + (-t159 * t523 + t160 * t522) * t527;
t601 = t42 * t637 + t88 * t636 + t838 * t638;
t366 = -Icges(5,5) * t462 - Icges(5,6) * t463;
t134 = t749 * t366 + t691 * t502 - t851 * t504;
t367 = Icges(5,5) * t465 - Icges(5,6) * t558;
t135 = t749 * t367 + t690 * t502 + t850 * t504;
t167 = -t450 * t713 - t462 * t678 - t849 * t463;
t168 = t450 * t715 + t465 * t678 - t849 * t558;
t198 = (t849 * t504 + t624 + t855) * t749;
t587 = ((t367 * t715 + t465 * t690) * t715 - (t366 * t715 + t465 * t691) * t713 + t168 * t749 + t843 * t558) * t638 + ((-t367 * t713 - t462 * t690) * t715 - (-t366 * t713 - t462 * t691) * t713 + t167 * t749 + t843 * t463) * t637 + (t198 + (-t134 * t523 + t135 * t522) * t527) * t657 + t665;
t559 = t233 * t239 - t249 * t267 + t251 * t268;
t101 = 0.4e1 * t559;
t571 = t101 * t786 + t587;
t564 = (t417 * t522 + t418 * t523) * t527;
t561 = -t453 * t527 - t611;
t555 = -t611 + t801;
t221 = ((-Icges(5,5) * t504 + Icges(5,6) * t502 + Icges(5,3) * t749) * t715 + t465 * t429 + t430 * t558) * t749;
t48 = t221 + ((t172 + t859) * t522 + ((t730 - t731) * t527 - t171 - t861) * t523) * t527;
t99 = t221 + (-t171 * t523 + t172 * t522) * t527;
t554 = t48 * t637 + t99 * t636 + t837 * t638 + t601;
t549 = t42 * t636 + t180 + t838 * t639 + (t128 + t150) * t638 + (t88 + t127 + t149) * t637;
t542 = m(6) * (t236 + 0.2e1 * t563);
t537 = t644 / 0.2e1 + t640 / 0.2e1 - t641 / 0.2e1 + t643 / 0.2e1;
t532 = t48 * t636 + t198 + t549 + t837 * t639 + (t135 + t168) * t638 + (t99 + t134 + t167) * t637;
t501 = (rSges(4,1) * t765 - rSges(4,2) * t763) * t527;
t470 = t749 * t477;
t412 = -Icges(4,5) * t548 - Icges(4,6) * t489;
t411 = -Icges(4,5) * t486 - Icges(4,6) * t487;
t358 = -t417 * t749 - t501 * t713;
t357 = t418 * t749 - t501 * t715;
t292 = t485 * t715 - t749 * (rSges(4,3) * t715 + t590);
t255 = (t622 + (t644 + t643 - t641 + t640) * t527) * t749;
t241 = t523 * t561 - t623 - t626;
t240 = t522 * t561 + t363 - t470;
t238 = t688 * t842;
t235 = qJD(5) * t761 / 0.2e1;
t229 = t489 * t673 + t498 * t715 - t548 * t672;
t228 = -t486 * t672 + t487 * t673 - t498 * t713;
t211 = t674 + t239;
t189 = t523 * t555 + t556 - t623;
t188 = t522 * t555 - t470 + t692;
t185 = t749 * t412 + (t681 * t765 + t683 * t763) * t527;
t184 = t749 * t411 + (t682 * t765 + t684 * t763) * t527;
t152 = -0.4e1 * t277 * t372 + 0.4e1 * t279 * t373;
t140 = t166 + t674;
t138 = 0.4e1 * t277 * t295 + 0.4e1 * t279 * t296;
t133 = 0.4e1 * t245 * t589 + 0.4e1 * t247 * t356;
t119 = 0.4e1 * t245 * t269 + 0.4e1 * t247 * t270;
t118 = 0.4e1 * t245 * t262 + 0.4e1 * t247 * t263;
t103 = t133 * t784 + t588;
t102 = t238 + t542 / 0.4e1 + t610 * t784;
t76 = t767 / 0.4e1;
t74 = t768 / 0.4e1;
t73 = 0.4e1 * t591;
t69 = t769 / 0.4e1;
t65 = t119 * t784 + t152 * t786 + t547;
t62 = t771 / 0.4e1;
t53 = 0.4e1 * t742 - 0.4e1 * t743 + 0.4e1 * t745;
t52 = 0.4e1 * t142 * t166 - 0.4e1 * t181 * t199 + 0.4e1 * t183 * t200;
t51 = t118 * t784 + t138 * t786 + t537 * t527 - t864;
t50 = 0.4e1 * t121 * t197 - 0.4e1 * t146 * t219 + 0.4e1 * t148 * t220;
t30 = t775 / 0.4e1;
t26 = t777 / 0.4e1;
t21 = t766 / 0.4e1 + t665;
t20 = t21 * qJD(5);
t19 = t770 / 0.4e1 + t665;
t18 = t53 * t784 + t665;
t16 = t52 * t784 + t571;
t15 = t50 * t784 + t73 * t786 + t587;
t14 = t76 - t769 / 0.4e1 + t601;
t13 = t69 - t767 / 0.4e1 + t601;
t12 = t30 - t777 / 0.4e1 + t616;
t11 = t26 - t775 / 0.4e1 + t616;
t10 = t26 + t30 - t774 / 0.4e1 + t665;
t9 = t74 - t771 / 0.4e1 + t601;
t8 = t62 - t768 / 0.4e1 + t601;
t7 = t69 + t76 + t549;
t6 = t62 + t74 + t549;
t4 = t554 + (t522 * t655 + t523 * t656) * t527;
t3 = t554 + t816;
t2 = t554 - t816;
t1 = t532 + t757 + t758;
t5 = [0, 0, ((t675 + t688) * t842 + (t610 + t675) * t841 + m(4) * t564) * qJD(3) + t102 * qJD(4) + t235, t102 * qJD(3) + (t238 + t542 / 0.2e1) * qJD(4) + t235, t235 + (qJD(3) / 0.2e1 + qJD(4) / 0.2e1) * t761; 0, t51 * qJD(3) + t65 * qJD(4) + t103 * qJD(5), t51 * qJD(2) + t1 * qJD(4) + t6 * qJD(5) + (-t146 * t263 + t148 * t262 + t188 * t247 + t189 * t245) * t752 + (-t214 * t296 + t216 * t295 + t240 * t279 + t241 * t277) * t755 + (t532 + t255 + (-t292 * t418 - t293 * t417 + t357 * t362 + t358 * t360) * m(4) + ((-t184 / 0.2e1 - t228 / 0.2e1 - t656) * t523 + (t185 / 0.2e1 + t229 / 0.2e1 - t655) * t522) * t527) * qJD(3), t65 * qJD(2) + t1 * qJD(3) + t532 * qJD(4) + t7 * qJD(5) + (-t181 * t270 + t183 * t269 + t155 / 0.2e1 + t156 / 0.2e1) * t751 + (-t249 * t373 - t251 * t372 + t202 / 0.2e1 + t203 / 0.2e1) * t754, t103 * qJD(2) + t6 * qJD(3) + t7 * qJD(4) + (t549 + (-t242 * t356 + t243 * t589 + t194 / 0.2e1 + t195 / 0.2e1) * m(6)) * qJD(5); -qJD(4) * t109, t4 * qJD(3) + t2 * qJD(4) + t8 * qJD(5) + t118 * t866 + t138 * t863 + (((t811 + t266 / 0.2e1) * t522 - t537) * t527 + t864) * qJD(2), t4 * qJD(2) + (((t412 * t715 + t489 * t683 - t548 * t681) * t715 - (t411 * t715 + t489 * t684 - t548 * t682) * t713 + t229 * t749) * t638 + ((-t412 * t713 - t486 * t681 + t487 * t683) * t715 - (-t411 * t713 - t486 * t682 + t487 * t684) * t713 + t228 * t749) * t637 + (t121 * t140 - t146 * t188 + t148 * t189) * m(6) + (t178 * t211 - t214 * t240 + t216 * t241) * m(5) + (t255 + (-t184 * t523 + t185 * t522) * t527) * t657 + ((t523 * t590 + (t402 + t512) * t522) * t527 * t564 - t292 * t357 + t293 * t358) * m(4) + t587) * qJD(3) + t15 * qJD(4) + t18 * qJD(5), -t670 + t2 * qJD(2) + t15 * qJD(3) + t587 * qJD(4) + t10 * qJD(5) + (-t101 / 0.4e1 + t559 + t591) * t754 + (-t52 / 0.4e1 + (t148 + t183) * t220 + (-t146 - t181) * t219 + (t121 + t142) * t197) * t751, t8 * qJD(2) + t18 * qJD(3) + t10 * qJD(4) + ((t650 + t653) * t841 + t615) * qJD(5); qJD(3) * t109, -t547 * qJD(2) + t3 * qJD(3) + qJD(4) * t554 + t13 * qJD(5) + t119 * t866 + t152 * t863, t670 + t3 * qJD(2) + t587 * qJD(3) + t16 * qJD(4) + t11 * qJD(5) + (t211 * t233 - t240 * t249 + t241 * t251 - t73 / 0.4e1 + t591) * t755 + (t121 * t166 + t140 * t142 - t146 * t199 + t148 * t200 - t181 * t188 + t183 * t189 - t50 / 0.4e1) * t752, t554 * qJD(2) + t16 * qJD(3) + ((t142 * t197 - t181 * t219 + t183 * t220) * m(6) + t571) * qJD(4) + t19 * qJD(5), t13 * qJD(2) + t11 * qJD(3) + t19 * qJD(4) + ((t650 + t652) * t841 + t615) * qJD(5); 0, -t588 * qJD(2) + t9 * qJD(3) + t14 * qJD(4) + qJD(5) * t601 + t133 * t866, t9 * qJD(2) + t665 * qJD(3) + t12 * qJD(4) + t20 + (t140 * t225 - t188 * t242 + t189 * t243 + t104 / 0.2e1 + t123 / 0.2e1 + t124 / 0.2e1 - t53 / 0.4e1) * t752, t14 * qJD(2) + t12 * qJD(3) + ((t651 + t652) * t841 - t770 / 0.4e1 + t665) * qJD(4) + t20, qJD(2) * t601 + t20 + (qJD(3) + qJD(4)) * t21;];
Cq = t5;
