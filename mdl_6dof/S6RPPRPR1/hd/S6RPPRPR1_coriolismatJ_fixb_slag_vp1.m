% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
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
% Cq [6x6]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:40
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPPRPR1_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_coriolismatJ_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_coriolismatJ_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_coriolismatJ_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR1_coriolismatJ_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRPR1_coriolismatJ_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRPR1_coriolismatJ_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:39:09
% EndTime: 2019-03-09 01:39:35
% DurationCPUTime: 22.04s
% Computational Cost: add. (72796->728), mult. (50639->1045), div. (0->0), fcn. (54118->11), ass. (0->423)
t460 = pkin(10) + qJ(4);
t454 = sin(t460);
t461 = qJ(1) + pkin(9);
t455 = sin(t461);
t615 = t454 * t455;
t719 = m(7) / 0.2e1;
t721 = m(6) / 0.2e1;
t550 = t719 + t721;
t451 = t455 ^ 2;
t458 = cos(t461);
t452 = t458 ^ 2;
t561 = t451 + t452;
t409 = t561 * t454;
t765 = m(6) + m(7);
t567 = t765 * t409 / 0.2e1;
t283 = -t454 * t550 + t567;
t457 = cos(t460);
t459 = pkin(11) + qJ(6);
t453 = sin(t459);
t456 = cos(t459);
t658 = rSges(7,1) * t456;
t517 = -rSges(7,2) * t453 + t658;
t360 = -rSges(7,3) * t457 + t454 * t517;
t629 = qJ(5) * t457;
t666 = pkin(4) * t454;
t422 = -t629 + t666;
t465 = -pkin(8) - qJ(5);
t597 = qJ(5) + t465;
t463 = cos(pkin(11));
t449 = pkin(5) * t463 + pkin(4);
t663 = -pkin(4) + t449;
t542 = t454 * t663 + t457 * t597 + t360 + t422;
t247 = t542 * t455;
t249 = t542 * t458;
t462 = sin(pkin(11));
t659 = rSges(6,1) * t463;
t520 = -rSges(6,2) * t462 + t659;
t486 = t520 * t454;
t754 = rSges(6,3) * t457 - t486;
t569 = t422 - t754;
t316 = t569 * t455;
t318 = t569 * t458;
t600 = t458 * t462;
t610 = t455 * t463;
t407 = -t457 * t600 + t610;
t599 = t458 * t463;
t611 = t455 * t462;
t408 = t457 * t599 + t611;
t521 = t408 * rSges(6,1) + t407 * rSges(6,2);
t662 = -pkin(7) - qJ(3);
t667 = cos(qJ(1)) * pkin(1);
t532 = -t455 * t662 + t667;
t450 = cos(pkin(10)) * pkin(3) + pkin(2);
t650 = rSges(6,3) + qJ(5);
t665 = t457 * pkin(4);
t728 = t454 * t650 + t450 + t665;
t231 = t458 * t728 + t521 + t532;
t604 = t457 * t458;
t612 = t455 * t457;
t664 = sin(qJ(1)) * pkin(1);
t495 = -t458 * t662 - t664;
t405 = t457 * t611 + t599;
t406 = t457 * t610 - t600;
t737 = -t406 * rSges(6,1) + t405 * rSges(6,2);
t756 = -t455 * t728 + t495 + t737;
t591 = t231 * t612 + t604 * t756;
t601 = t458 * t453;
t384 = t455 * t456 - t457 * t601;
t385 = t453 * t455 + t456 * t604;
t518 = t385 * rSges(7,1) + t384 * rSges(7,2);
t613 = t454 * t458;
t524 = pkin(5) * t611 - t465 * t613;
t618 = t449 * t457;
t656 = rSges(7,3) * t454;
t200 = (t450 + t618 + t656) * t458 + t518 + t524 + t532;
t412 = t449 * t612;
t649 = -rSges(7,3) + t465;
t607 = t456 * t458;
t382 = t453 * t612 + t607;
t383 = t456 * t612 - t601;
t738 = -t383 * rSges(7,1) + t382 * rSges(7,2);
t757 = (t454 * t649 - t450) * t455 + (pkin(5) * t462 - t662) * t458 - t412 - t664 + t738;
t594 = t200 * t612 + t604 * t757;
t647 = (-t316 * t613 + t318 * t615 + t591) * t721 + (-t247 * t613 + t249 * t615 + t594) * t719;
t262 = (t649 * t457 + (t449 + t517) * t454) * t455;
t565 = t454 * rSges(7,2) * t601 + rSges(7,3) * t604;
t603 = t457 * t465;
t263 = (-t603 + (-t449 - t658) * t454) * t458 + t565;
t435 = pkin(4) * t615;
t284 = t435 + (-t457 * t650 + t486) * t455;
t427 = qJ(5) * t604;
t563 = t454 * rSges(6,2) * t600 + rSges(6,3) * t604;
t285 = t427 + (-pkin(4) - t659) * t613 + t563;
t648 = ((t284 * t458 + t285 * t455) * t454 + t591) * t721 + ((t262 * t458 + t263 * t455) * t454 + t594) * t719;
t755 = t647 - t648;
t778 = -t755 * qJD(1) - t283 * qJD(2);
t508 = Icges(7,5) * t456 - Icges(7,6) * t453;
t348 = -Icges(7,3) * t457 + t454 * t508;
t637 = Icges(7,4) * t456;
t512 = -Icges(7,2) * t453 + t637;
t352 = -Icges(7,6) * t457 + t454 * t512;
t638 = Icges(7,4) * t453;
t514 = Icges(7,1) * t456 - t638;
t356 = -Icges(7,5) * t457 + t454 * t514;
t184 = t348 * t615 - t352 * t382 + t356 * t383;
t270 = Icges(7,5) * t383 - Icges(7,6) * t382 + Icges(7,3) * t615;
t370 = Icges(7,4) * t383;
t273 = -Icges(7,2) * t382 + Icges(7,6) * t615 + t370;
t369 = Icges(7,4) * t382;
t277 = -Icges(7,1) * t383 - Icges(7,5) * t615 + t369;
t136 = t270 * t615 - t273 * t382 - t277 * t383;
t272 = Icges(7,5) * t385 + Icges(7,6) * t384 + Icges(7,3) * t613;
t639 = Icges(7,4) * t385;
t275 = Icges(7,2) * t384 + Icges(7,6) * t613 + t639;
t371 = Icges(7,4) * t384;
t278 = Icges(7,1) * t385 + Icges(7,5) * t613 + t371;
t137 = t272 * t615 - t382 * t275 + t383 * t278;
t507 = t136 * t455 + t137 * t458;
t17 = t457 * t184 - t454 * t507;
t279 = rSges(7,3) * t615 - t738;
t362 = t457 * t517 + t656;
t166 = (t362 * t455 - t279) * t454;
t281 = rSges(7,3) * t613 + t518;
t328 = -rSges(7,1) * t454 * t607 + t565;
t167 = (-t360 * t458 - t328) * t457 + (-t362 * t458 + t281) * t454;
t104 = t166 * t455 - t167 * t458;
t423 = rSges(5,1) * t454 + rSges(5,2) * t457;
t402 = t423 * t455;
t404 = t423 * t458;
t302 = t402 * t455 + t404 * t458;
t722 = m(5) / 0.2e1;
t768 = t302 * t722;
t543 = (t262 * t455 - t263 * t458) * t719 + (t284 * t455 - t285 * t458) * t721 + t768;
t734 = t247 * t455 + t249 * t458;
t742 = t561 * t423;
t544 = -t734 * m(7) / 0.2e1 + (-t316 * t455 - t318 * t458) * t721 - m(5) * t742 / 0.2e1;
t38 = t544 - t543;
t661 = m(7) * qJD(6);
t777 = t38 * qJD(1) - t104 * t661 / 0.2e1;
t312 = -rSges(7,1) * t382 - rSges(7,2) * t383;
t313 = rSges(7,1) * t384 - rSges(7,2) * t385;
t776 = m(7) * (t200 * t313 - t312 * t757);
t186 = t348 * t613 + t352 * t384 + t356 * t385;
t138 = t270 * t613 + t384 * t273 - t277 * t385;
t139 = t272 * t613 + t384 * t275 + t385 * t278;
t506 = t138 * t455 + t139 * t458;
t774 = -t457 * t186 + t454 * t506;
t291 = Icges(6,5) * t406 - Icges(6,6) * t405 + Icges(6,3) * t615;
t294 = Icges(6,4) * t406 - Icges(6,2) * t405 + Icges(6,6) * t615;
t297 = Icges(6,1) * t406 - Icges(6,4) * t405 + Icges(6,5) * t615;
t752 = t407 * t294 + t408 * t297;
t156 = t291 * t613 + t752;
t293 = Icges(6,5) * t408 + Icges(6,6) * t407 + Icges(6,3) * t613;
t296 = Icges(6,4) * t408 + Icges(6,2) * t407 + Icges(6,6) * t613;
t299 = Icges(6,1) * t408 + Icges(6,4) * t407 + Icges(6,5) * t613;
t770 = -t138 * t458 + t139 * t455;
t773 = (t293 * t613 + t407 * t296 + t408 * t299) * t455 - t156 * t458 + t770;
t75 = -t136 * t458 + t137 * t455;
t226 = t281 * t457 + t360 * t613;
t769 = t279 * t457 + t360 * t615;
t740 = t226 * t455 - t458 * t769;
t624 = t270 * t457;
t767 = t273 * t453 + t277 * t456;
t161 = t454 * t767 + t624;
t741 = -t200 * t455 - t458 * t757;
t696 = -t457 / 0.2e1;
t501 = -t294 * t405 + t297 * t406;
t751 = -0.2e1 * t409;
t700 = t454 / 0.2e1;
t746 = t360 * t455;
t509 = Icges(6,5) * t463 - Icges(6,6) * t462;
t745 = t454 * t509;
t744 = t454 * t597;
t196 = (t312 * t458 - t313 * t455) * t454;
t551 = m(6) / 0.4e1 + m(7) / 0.4e1;
t614 = t454 * t457;
t564 = t561 * t614;
t743 = t551 * (t564 - t614);
t660 = rSges(5,1) * t457;
t530 = t450 + t660;
t562 = rSges(5,2) * t615 + t458 * rSges(5,3);
t314 = -t455 * t530 + t495 + t562;
t529 = -rSges(5,2) * t613 + t455 * rSges(5,3);
t315 = t458 * t530 + t529 + t532;
t739 = t314 * t458 + t315 * t455;
t631 = Icges(6,3) * t457;
t477 = t631 - t745;
t490 = t294 * t462 - t297 * t463 + t477 * t455;
t430 = Icges(5,4) * t615;
t358 = Icges(5,1) * t612 - Icges(5,5) * t458 - t430;
t575 = -Icges(5,2) * t612 + t358 - t430;
t736 = -t490 + t575;
t732 = t293 * t615 - t405 * t296 + t406 * t299;
t513 = Icges(6,4) * t463 - Icges(6,2) * t462;
t730 = -Icges(6,6) * t457 + t454 * t513;
t515 = Icges(6,1) * t463 - Icges(6,4) * t462;
t729 = -Icges(6,5) * t457 + t454 * t515;
t447 = Icges(5,4) * t457;
t634 = Icges(5,2) * t454;
t355 = Icges(5,6) * t455 + (t447 - t634) * t458;
t641 = Icges(5,1) * t454;
t516 = -t447 - t641;
t577 = t516 * t458 - t355;
t354 = Icges(5,4) * t612 - Icges(5,2) * t615 - Icges(5,6) * t458;
t578 = -t516 * t455 + t354;
t727 = (t577 * t455 + t458 * t578) * t457;
t395 = (-Icges(7,2) * t456 - t638) * t454;
t398 = (-Icges(7,1) * t453 - t637) * t454;
t726 = -t453 * (t356 / 0.2e1 + t395 / 0.2e1) + t456 * (t398 / 0.2e1 - t352 / 0.2e1);
t724 = 0.4e1 * qJD(1);
t723 = 2 * qJD(4);
t717 = -t774 / 0.2e1;
t716 = t75 / 0.2e1;
t715 = t770 / 0.2e1;
t502 = t279 * t458 - t281 * t455;
t140 = t502 * t457 + (-t328 * t455 - t458 * t746) * t454;
t182 = t502 * t454;
t592 = -t226 * t612 + t604 * t769;
t712 = m(7) * (-t140 * t457 + (t166 * t458 + t167 * t455 + t182) * t454 + t592);
t710 = m(7) * (t166 * t757 + t167 * t200 - t226 * t263 + t262 * t769);
t401 = (-rSges(7,1) * t453 - rSges(7,2) * t456) * t454;
t708 = m(7) * (-t247 * t313 + t249 * t312 + t401 * t741);
t703 = m(7) * (t182 * t409 + t592);
t306 = -Icges(7,5) * t382 - Icges(7,6) * t383;
t585 = -Icges(7,2) * t383 - t277 - t369;
t587 = -Icges(7,1) * t382 - t273 - t370;
t119 = -t306 * t457 + (-t453 * t585 + t456 * t587) * t454;
t702 = t119 / 0.2e1;
t698 = t455 / 0.2e1;
t697 = t455 / 0.4e1;
t695 = -t458 / 0.2e1;
t694 = -t458 / 0.4e1;
t693 = t458 / 0.2e1;
t651 = rSges(4,3) + qJ(3);
t692 = m(4) * ((t651 * t458 - t664) * t458 + (t651 * t455 + t667) * t455);
t691 = m(5) * (t314 * t402 - t315 * t404);
t690 = m(5) * t739;
t685 = m(6) * (t231 * t285 + t284 * t756);
t221 = t231 * t613;
t683 = m(6) * (-t615 * t756 + t221);
t682 = m(6) * (t231 * t455 + t458 * t756);
t678 = m(7) * (t200 * t263 + t262 * t757);
t191 = t200 * t613;
t677 = m(7) * (-t615 * t757 + t191);
t676 = m(7) * t741;
t675 = m(7) * (-t226 * t613 - t615 * t769);
t674 = m(7) * t740;
t204 = t312 * t455 + t313 * t458;
t671 = m(7) * (-t204 * t457 - t401 * t409);
t669 = m(7) * t196;
t668 = m(7) * t204;
t654 = t455 * t17;
t652 = t458 * t774;
t640 = Icges(5,4) * t454;
t323 = t352 * t455;
t325 = t356 * t455;
t492 = -t348 * t455 + t767;
t110 = -t492 * t457 + (t323 * t453 - t325 * t456 + t270) * t454;
t628 = t110 * t458;
t324 = t352 * t458;
t326 = t356 * t458;
t503 = -t275 * t453 + t278 * t456;
t491 = -t348 * t458 - t503;
t111 = -t491 * t457 + (t324 * t453 - t326 * t456 + t272) * t454;
t627 = t111 * t455;
t623 = t272 * t457;
t621 = t348 * t457;
t619 = t354 * t454;
t617 = t453 * t352;
t353 = Icges(7,6) * t454 + t457 * t512;
t616 = t453 * t353;
t609 = t456 * t356;
t357 = Icges(7,5) * t454 + t457 * t514;
t608 = t456 * t357;
t392 = (-Icges(7,5) * t453 - Icges(7,6) * t456) * t454;
t605 = t457 * t392;
t588 = -t247 * t612 - t249 * t604;
t586 = Icges(7,1) * t384 - t275 - t639;
t584 = -Icges(7,2) * t385 + t278 + t371;
t582 = -t316 * t612 - t318 * t604;
t350 = Icges(5,5) * t612 - Icges(5,6) * t615 - Icges(5,3) * t458;
t581 = -t455 * t350 - t358 * t604;
t511 = Icges(5,5) * t457 - Icges(5,6) * t454;
t351 = Icges(5,3) * t455 + t458 * t511;
t421 = Icges(5,1) * t457 - t640;
t359 = Icges(5,5) * t455 + t421 * t458;
t580 = t455 * t351 + t359 * t604;
t579 = -t352 + t398;
t576 = t356 + t395;
t418 = Icges(5,2) * t457 + t640;
t574 = -t418 * t458 + t359;
t573 = t455 * (qJ(5) * t612 - t435) + t458 * (-pkin(4) * t613 + t427);
t424 = qJ(5) * t454 + t665;
t572 = t561 * t424;
t568 = -rSges(6,3) * t454 - t457 * t520 - t424;
t560 = qJD(1) * t454;
t559 = qJD(1) * t457;
t558 = qJD(6) * t454;
t557 = qJD(6) * t457;
t556 = t104 * qJD(3);
t108 = 0.2e1 * (t140 / 0.4e1 - t204 / 0.4e1) * m(7);
t555 = t108 * qJD(2);
t236 = t550 * t751;
t553 = t236 * qJD(1);
t546 = t717 + t774 / 0.2e1;
t541 = -t457 * t663 - t362 - t424 + t744;
t540 = t615 / 0.4e1;
t527 = t574 * t455;
t339 = t359 * t612;
t526 = t351 * t458 - t339;
t525 = t355 * t454 - t350;
t123 = (t279 - pkin(5) * t600 + t412 - (t665 + t744) * t455) * t455 + (t281 + (-t424 + t618) * t458 + t524) * t458 + t572;
t173 = t455 * (rSges(6,3) * t615 - t737) + t458 * (rSges(6,3) * t613 + t521) + t572;
t44 = m(6) * (t173 * t409 + t582) + m(7) * (t123 * t409 + t588);
t307 = Icges(7,5) * t384 - Icges(7,6) * t385;
t120 = -t307 * t457 + (-t453 * t584 + t456 * t586) * t454;
t144 = -t382 * t576 + t383 * t579 + t392 * t615;
t145 = t384 * t576 + t385 * t579 + t392 * t613;
t523 = t708 / 0.2e1 + (t120 + t145) * t697 + (t119 + t144) * t694;
t510 = -Icges(5,5) * t454 - Icges(5,6) * t457;
t162 = t454 * t503 - t623;
t505 = -t161 * t455 + t162 * t458;
t499 = t609 - t617;
t497 = -t449 * t454 - t603;
t100 = t306 * t615 - t382 * t585 + t383 * t587;
t101 = t307 * t615 - t382 * t584 + t383 * t586;
t47 = -t100 * t458 + t101 * t455;
t102 = t306 * t613 + t384 * t585 + t385 * t587;
t103 = t307 * t613 + t384 * t584 + t385 * t586;
t48 = -t102 * t458 + t103 * t455;
t496 = t47 * t695 + t48 * t698;
t489 = t296 * t462 - t299 * t463 + t477 * t458;
t349 = Icges(7,3) * t454 + t457 * t508;
t488 = t349 - t499;
t484 = (-t291 * t458 + t293 * t455) * t457;
t483 = t489 * t455;
t482 = t17 * t697 + t774 * t694 - t654 / 0.4e1 + t652 / 0.4e1 + (t540 - t615 / 0.4e1) * t770;
t475 = t454 * t492 + t624;
t474 = t454 * t491 + t623;
t473 = t454 * t488 + t621;
t130 = -t353 * t382 + t357 * t383 + t455 * t473;
t131 = t353 * t384 + t357 * t385 + t458 * t473;
t158 = -t488 * t457 + (t348 + t608 - t616) * t454;
t214 = t454 * t499 - t621;
t472 = t158 * t696 + t214 * t700 + t710 / 0.2e1 + (t110 + t130) * t540 + (t111 + t131) * t613 / 0.4e1 + (-t161 + t184) * t612 / 0.4e1 + (t162 + t186) * t604 / 0.4e1;
t377 = Icges(6,6) * t454 + t457 * t513;
t379 = Icges(6,5) * t454 + t457 * t515;
t469 = t608 / 0.2e1 - t616 / 0.2e1 + t348 / 0.2e1 + t421 / 0.2e1 - t418 / 0.2e1 + t463 * t379 / 0.2e1 - t462 * t377 / 0.2e1 + t745 / 0.2e1 - t631 / 0.2e1;
t468 = -t463 * t729 / 0.2e1 + t462 * t730 / 0.2e1 - t609 / 0.2e1 + t617 / 0.2e1 + t349 / 0.2e1 - t447 - t641 / 0.2e1 + t634 / 0.2e1 + Icges(6,3) * t700 + t457 * t509 / 0.2e1;
t425 = -rSges(5,2) * t454 + t660;
t394 = t510 * t458;
t393 = t510 * t455;
t338 = t729 * t458;
t337 = t729 * t455;
t336 = t730 * t458;
t335 = t730 * t455;
t319 = t568 * t458;
t317 = t568 * t455;
t282 = t700 * t765 + t567;
t251 = 0.4e1 * t743;
t250 = t541 * t458;
t248 = t541 * t455;
t240 = -t313 * t457 - t401 * t613;
t239 = t312 * t457 + t401 * t615;
t235 = t551 * t751 + t567;
t218 = -t355 * t613 + t580;
t217 = -t354 * t613 - t581;
t216 = -t355 * t615 - t526;
t202 = -t668 / 0.2e1;
t192 = -t669 / 0.2e1;
t187 = t458 * (-rSges(6,1) * t454 * t599 + t563) + t754 * t451 + t573;
t175 = -t605 + (-t453 * t576 + t456 * t579) * t454;
t168 = t671 / 0.2e1;
t149 = -t674 / 0.2e1;
t148 = -t217 * t458 + t218 * t455;
t147 = -(-t455 * (-t358 * t457 + t619) - t350 * t458) * t458 + t216 * t455;
t143 = t675 / 0.2e1;
t133 = (-t427 + t328 + (t497 + t666) * t458) * t458 + (t435 - t746 + (t497 - t629) * t455) * t455 + t573;
t109 = (t140 + t204) * t719;
t99 = t104 * qJD(4) * t719;
t97 = t703 / 0.2e1;
t96 = -t324 * t384 - t326 * t385 + t458 * t474;
t95 = -t323 * t384 - t325 * t385 + t458 * t475;
t94 = t324 * t382 - t326 * t383 + t455 * t474;
t93 = t323 * t382 - t325 * t383 + t455 * t475;
t84 = -(t291 * t615 + t501) * t458 + t455 * t732;
t81 = t677 + t683;
t70 = t776 - t605 / 0.2e1 + t726 * t454;
t67 = -t214 * t457 + t454 * t505;
t60 = t123 * t204 + t401 * t734;
t55 = (t216 - t339 + (t351 + t619) * t458 + t581) * t458 + t580 * t455;
t54 = (t458 * t525 + t218 - t580) * t458 + (t455 * t525 + t217 + t526) * t455;
t53 = -t676 + t682 + t690 + t692;
t51 = t149 + t668 / 0.2e1;
t50 = t202 + t149;
t49 = t202 + t674 / 0.2e1;
t46 = t455 * t96 - t458 * t95;
t45 = t455 * t94 - t458 * t93;
t41 = t143 + t669 / 0.2e1;
t40 = t192 + t143;
t39 = t192 - t675 / 0.2e1;
t36 = t543 + t544;
t34 = t140 * t182 + t166 * t769 - t167 * t226;
t32 = t712 / 0.2e1;
t31 = -t145 * t457 + (t102 * t455 + t103 * t458) * t454;
t30 = -t144 * t457 + (t100 * t455 + t101 * t458) * t454;
t29 = t454 * t469 - t457 * t468 + t678 + t685 + t691;
t27 = t501 * t458 + (t156 - t732 - t752) * t455;
t22 = (-t158 + t505) * t457 + (t110 * t455 + t111 * t458 + t214) * t454;
t21 = t97 + t32 - t671 / 0.2e1;
t20 = t168 + t97 - t712 / 0.2e1;
t19 = t168 + t32 - t703 / 0.2e1;
t12 = (-t131 + t506) * t457 + (t455 * t95 + t458 * t96 + t186) * t454;
t11 = (-t130 + t507) * t457 + (t455 * t93 + t458 * t94 + t184) * t454;
t10 = m(7) * t60 + t496;
t8 = t647 + t648;
t6 = t546 * t615;
t5 = m(7) * t34 + (t652 / 0.2e1 - t654 / 0.2e1 - t22 / 0.2e1) * t457 + (t12 * t693 + t11 * t698 + t67 / 0.2e1) * t454;
t4 = (t715 - t770 / 0.2e1 + t148 / 0.2e1 - t55 / 0.2e1) * t458 + (-t75 / 0.2e1 + t716 + t54 / 0.2e1 + t147 / 0.2e1 + t84 / 0.2e1 + t27 / 0.2e1) * t455;
t3 = t472 + t523;
t2 = (t158 / 0.2e1 + (-t186 / 0.4e1 - t162 / 0.4e1) * t458 + (-t184 / 0.4e1 + t161 / 0.4e1) * t455) * t457 + t482 - t710 / 0.2e1 + (-t214 / 0.2e1 + (-t131 / 0.4e1 - t111 / 0.4e1) * t458 + (-t110 / 0.4e1 - t130 / 0.4e1) * t455) * t454 + t523;
t1 = t472 + t482 + (-t145 / 0.4e1 - t120 / 0.4e1) * t455 + (t144 / 0.4e1 + t119 / 0.4e1) * t458 - t708 / 0.2e1;
t7 = [t53 * qJD(3) + t29 * qJD(4) + t81 * qJD(5) + t70 * qJD(6), 0, qJD(1) * t53 + qJD(4) * t36 + qJD(5) * t235 + qJD(6) * t50, t29 * qJD(1) + t36 * qJD(3) + t8 * qJD(5) + t3 * qJD(6) + ((-t739 * t425 + (-t402 * t458 + t404 * t455) * t423) * t722 + (t231 * t317 - t284 * t318 - t285 * t316 + t319 * t756) * t721 + (t200 * t248 - t247 * t263 - t249 * t262 + t250 * t757) * t719) * t723 + (t627 / 0.2e1 + (t452 / 0.2e1 + t451 / 0.2e1) * t511 - t628 / 0.2e1 + (t55 + t773) * t693 + (t377 * t407 + t379 * t408 + t131 + (t574 - t489) * t457 + (t336 * t462 - t338 * t463 + t293 + t577) * t454) * t698 - (t54 + t147 + t84 + t27) * t455 / 0.2e1 + (-t377 * t405 + t379 * t406 + t130 + t148 + t736 * t457 + (t335 * t462 - t337 * t463 + t291 - t578) * t454 + t773) * t695) * qJD(4), qJD(1) * t81 + qJD(3) * t235 + qJD(4) * t8 + qJD(6) * t40, t70 * qJD(1) + t50 * qJD(3) + t3 * qJD(4) + t40 * qJD(5) - t175 * t557 + (t200 * t240 - t226 * t313 + t239 * t757 - t312 * t769) * t661 + ((t120 / 0.2e1 + t145 / 0.2e1) * t458 + (t144 / 0.2e1 + t702 - t546) * t455) * t558; 0, 0, 0 (t133 * t719 + t187 * t721 - t768) * t723 + t282 * qJD(5) + t109 * qJD(6), t282 * qJD(4), t109 * qJD(4) + t196 * t661; -t38 * qJD(4) + t236 * qJD(5) + t49 * qJD(6) + (t676 / 0.4e1 - t682 / 0.4e1 - t690 / 0.4e1 - t692 / 0.4e1) * t724, 0, 0 ((-t317 * t458 + t319 * t455) * t721 + (-t248 * t458 + t250 * t455) * t719) * t723 - t777, t553, t49 * qJD(1) + t99 + (t239 * t455 - t240 * t458) * t661; t38 * qJD(3) + t4 * qJD(4) + t755 * qJD(5) + t2 * qJD(6) + (-t691 / 0.4e1 - t685 / 0.4e1 - t678 / 0.4e1) * t724 + t468 * t559 - t469 * t560, qJD(5) * t283 - qJD(6) * t108, t777, t4 * qJD(1) + t44 * qJD(5) + t10 * qJD(6) + (m(7) * (t123 * t133 - t247 * t248 - t249 * t250) + m(6) * (t173 * t187 - t316 * t317 - t318 * t319) + m(5) * (t425 * t742 - (t455 * (rSges(5,1) * t612 - t562) + t458 * (rSges(5,1) * t604 + t529)) * t302) + ((-t336 * t407 - t338 * t408) * t455 + (t407 * t335 + t408 * t337 + t484 + (-t458 * t490 + t483) * t454) * t458 + t46 + t451 * t394 + (-t455 * t393 + t727 + (t458 * t575 - t527) * t454) * t458) * t698 + (t45 + t452 * t393 - (t335 * t405 - t337 * t406) * t458 + (-t458 * t394 + t727 + t405 * t336 - t406 * t338 + t484 + (t458 * t736 + t483 - t527) * t454) * t455) * t695) * qJD(4), t44 * qJD(4) + (-0.4e1 * t743 + 0.2e1 * t550 * (-t409 * t457 + t564)) * qJD(5) + t20 * qJD(6) - t778, t2 * qJD(1) - t555 + t10 * qJD(4) + t20 * qJD(5) + (t30 * t695 + t31 * t698) * qJD(6) + (t22 / 0.2e1 + (t702 + t717) * t458 + (-t120 / 0.2e1 + t17 / 0.2e1) * t455) * t557 + (-t556 / 0.2e1 + (t123 * t196 + t182 * t204 - t239 * t249 - t240 * t247 + t401 * t740 - t34) * qJD(6)) * m(7) + (-t67 / 0.2e1 + (t48 / 0.2e1 - t12 / 0.2e1) * t458 + (t47 / 0.2e1 - t11 / 0.2e1) * t455) * t558; -t236 * qJD(3) - t755 * qJD(4) + t39 * qJD(6) + (-t677 / 0.4e1 - t683 / 0.4e1) * t724 + 0.2e1 * (t191 * t719 + t221 * t721 + (-t200 * t719 - t231 * t721) * t613) * qJD(1), -t283 * qJD(4), -t553 (m(7) * (-t133 * t457 + t588) + m(6) * (-t187 * t457 + t582) + 0.2e1 * ((t248 * t455 + t250 * t458 + t123) * t719 + (t317 * t455 + t319 * t458 + t173) * t721) * t454 - t44) * qJD(4) + t251 * qJD(5) + t19 * qJD(6) + t778, t251 * qJD(4), t39 * qJD(1) + t19 * qJD(4) + (-t196 * t457 + (t239 * t458 + t240 * t455) * t454) * t661; t392 * t559 / 0.2e1 + t51 * qJD(3) + t1 * qJD(4) + t41 * qJD(5) + t6 * qJD(6) - qJD(1) * t776 - t726 * t560, qJD(4) * t108, qJD(1) * t51 + t99, t1 * qJD(1) + t555 + (t12 * t698 + t11 * t695 + (t161 * t458 + t162 * t455) * t700 + (t627 - t628) * t696 + t612 * t716 + t45 * t615 / 0.2e1 + t604 * t715 + t46 * t613 / 0.2e1 - t496) * qJD(4) + t21 * qJD(5) + t5 * qJD(6) + (t556 / 0.2e1 + (t123 * t140 + t133 * t182 - t166 * t249 - t167 * t247 - t226 * t248 + t250 * t769 - t60) * qJD(4)) * m(7), qJD(1) * t41 + qJD(4) * t21, t6 * qJD(1) + t5 * qJD(4) + (t457 ^ 2 * t175 / 0.2e1 + m(7) * (t182 * t196 - t226 * t240 + t239 * t769) + (t31 * t693 + t30 * t698 + (t119 * t455 + t120 * t458) * t696) * t454) * qJD(6);];
Cq  = t7;