% Calculate vector of inverse dynamics joint torques for
% S5PRPRP4
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-05 15:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP4_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP4_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP4_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP4_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP4_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP4_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP4_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:35:04
% EndTime: 2019-12-05 15:36:04
% DurationCPUTime: 45.49s
% Computational Cost: add. (15556->748), mult. (23383->1044), div. (0->0), fcn. (22407->8), ass. (0->392)
t700 = Icges(6,4) + Icges(5,5);
t699 = Icges(5,6) - Icges(6,6);
t678 = Icges(5,2) + Icges(6,3);
t721 = Icges(6,2) + Icges(5,3);
t335 = qJ(2) + pkin(8);
t330 = sin(t335);
t339 = sin(qJ(4));
t341 = cos(qJ(4));
t720 = (-t700 * t339 - t699 * t341) * t330;
t336 = sin(pkin(7));
t331 = cos(t335);
t512 = qJD(4) * t341;
t495 = t331 * t512;
t337 = cos(pkin(7));
t514 = qJD(4) * t337;
t522 = qJD(2) * t336;
t160 = -t336 * t495 + (t330 * t522 + t514) * t339;
t541 = t337 * t341;
t544 = t336 * t339;
t271 = t331 * t544 + t541;
t518 = qJD(2) * t341;
t498 = t330 * t518;
t161 = -qJD(4) * t271 - t336 * t498;
t497 = t331 * t522;
t719 = t699 * t160 + t700 * t161 + t721 * t497;
t542 = t337 * t339;
t500 = t330 * t542;
t513 = qJD(4) * t339;
t162 = qJD(2) * t500 - t336 * t513 - t337 * t495;
t543 = t336 * t341;
t273 = t331 * t542 - t543;
t163 = -qJD(4) * t273 - t337 * t498;
t520 = qJD(2) * t337;
t496 = t331 * t520;
t718 = t699 * t162 + t700 * t163 + t721 * t496;
t717 = Icges(5,1) + Icges(6,1);
t716 = Icges(5,4) - Icges(6,5);
t436 = Icges(5,5) * t341 - Icges(5,6) * t339;
t179 = Icges(5,3) * t330 + t331 * t436;
t441 = Icges(6,4) * t341 + Icges(6,6) * t339;
t181 = Icges(6,2) * t330 + t331 * t441;
t715 = -t720 * qJD(4) + (-t179 - t181) * qJD(2);
t272 = t331 * t543 - t542;
t552 = t330 * t336;
t110 = Icges(5,5) * t272 - Icges(5,6) * t271 + Icges(5,3) * t552;
t112 = Icges(6,4) * t272 + Icges(6,2) * t552 + Icges(6,6) * t271;
t686 = t110 + t112;
t274 = t331 * t541 + t544;
t551 = t330 * t337;
t111 = Icges(5,5) * t274 - Icges(5,6) * t273 + Icges(5,3) * t551;
t113 = Icges(6,4) * t274 + Icges(6,2) * t551 + Icges(6,6) * t273;
t685 = t111 + t113;
t559 = Icges(6,5) * t341;
t435 = Icges(6,3) * t339 + t559;
t177 = Icges(6,6) * t330 + t331 * t435;
t565 = Icges(5,4) * t341;
t442 = -Icges(5,2) * t339 + t565;
t183 = Icges(5,6) * t330 + t331 * t442;
t645 = -t177 + t183;
t178 = -Icges(5,3) * t331 + t330 * t436;
t180 = -Icges(6,2) * t331 + t330 * t441;
t704 = -t180 - t178;
t560 = Icges(6,5) * t339;
t447 = Icges(6,1) * t341 + t560;
t185 = Icges(6,4) * t330 + t331 * t447;
t566 = Icges(5,4) * t339;
t448 = Icges(5,1) * t341 - t566;
t187 = Icges(5,5) * t330 + t331 * t448;
t644 = -t185 - t187;
t714 = (-t678 * t341 + t560 - t566) * t330;
t340 = sin(qJ(2));
t342 = cos(qJ(2));
t439 = -Icges(3,5) * t340 - Icges(3,6) * t342;
t682 = -Icges(4,5) * t330 - Icges(4,6) * t331 + t439;
t713 = -t678 * t160 - t716 * t161 - t699 * t497;
t712 = -t678 * t162 - t716 * t163 - t699 * t496;
t711 = t716 * t160 + t717 * t161 + t700 * t497;
t710 = t716 * t162 + t717 * t163 + t700 * t496;
t709 = t645 * qJD(2) + t714 * qJD(4);
t255 = (-Icges(5,1) * t339 - t565) * t330;
t516 = qJD(4) * t330;
t708 = -(-Icges(6,1) * t339 + t559) * t516 - qJD(4) * t255 + t644 * qJD(2);
t238 = Icges(6,5) * t272;
t108 = Icges(6,6) * t552 + Icges(6,3) * t271 + t238;
t568 = Icges(5,4) * t272;
t114 = -Icges(5,2) * t271 + Icges(5,6) * t552 + t568;
t688 = t108 - t114;
t239 = Icges(6,5) * t274;
t109 = Icges(6,6) * t551 + Icges(6,3) * t273 + t239;
t567 = Icges(5,4) * t274;
t115 = -Icges(5,2) * t273 + Icges(5,6) * t551 + t567;
t687 = t109 - t115;
t562 = Icges(6,5) * t271;
t116 = Icges(6,1) * t272 + Icges(6,4) * t552 + t562;
t240 = Icges(5,4) * t271;
t118 = Icges(5,1) * t272 + Icges(5,5) * t552 - t240;
t707 = t116 + t118;
t561 = Icges(6,5) * t273;
t117 = Icges(6,1) * t274 + Icges(6,4) * t551 + t561;
t241 = Icges(5,4) * t273;
t119 = Icges(5,1) * t274 + Icges(5,5) * t551 - t241;
t706 = t117 + t119;
t549 = t330 * t341;
t313 = Icges(6,5) * t549;
t550 = t330 * t339;
t556 = Icges(6,6) * t331;
t176 = Icges(6,3) * t550 + t313 - t556;
t182 = -Icges(5,6) * t331 + t330 * t442;
t705 = t176 - t182;
t184 = -Icges(6,4) * t331 + t330 * t447;
t186 = -Icges(5,5) * t331 + t330 * t448;
t676 = t184 + t186;
t523 = qJD(2) * t331;
t703 = t715 * t330 + t704 * t523;
t702 = t718 * t330 + t685 * t523;
t701 = t719 * t330 + t686 * t523;
t698 = Icges(3,3) + Icges(4,3);
t697 = t682 * qJD(2);
t696 = Icges(3,5) * t342 + Icges(4,5) * t331 - Icges(3,6) * t340 - Icges(4,6) * t330;
t659 = -t688 * t160 + t707 * t161 + t713 * t271 + t711 * t272 + t701 * t336;
t658 = -t687 * t160 + t706 * t161 + t712 * t271 + t710 * t272 + t702 * t336;
t657 = -t688 * t162 + t707 * t163 + t713 * t273 + t711 * t274 + t701 * t337;
t656 = -t687 * t162 + t706 * t163 + t712 * t273 + t710 * t274 + t702 * t337;
t695 = t705 * t160 - t676 * t161 + t709 * t271 + t708 * t272 + t703 * t336;
t694 = t705 * t162 - t676 * t163 + t709 * t273 + t708 * t274 + t703 * t337;
t681 = t688 * t271 + t707 * t272 + t686 * t552;
t693 = t687 * t271 + t706 * t272 + t685 * t552;
t692 = t688 * t273 + t707 * t274 + t686 * t551;
t653 = t687 * t273 + t706 * t274 + t685 * t551;
t434 = t108 * t339 + t116 * t341;
t48 = -t112 * t331 + t330 * t434;
t430 = -t114 * t339 + t118 * t341;
t50 = -t110 * t331 + t330 * t430;
t652 = t48 + t50;
t433 = t109 * t339 + t117 * t341;
t49 = -t113 * t331 + t330 * t433;
t429 = -t115 * t339 + t119 * t341;
t51 = -t111 * t331 + t330 * t429;
t651 = t49 + t51;
t691 = t705 * t271 + t676 * t272 - t704 * t552;
t690 = t705 * t273 + t676 * t274 - t704 * t551;
t427 = t176 * t339 + t184 * t341;
t553 = t180 * t331;
t67 = t330 * t427 - t553;
t426 = -t182 * t339 + t186 * t341;
t554 = t178 * t331;
t68 = t330 * t426 - t554;
t625 = t67 + t68;
t689 = t337 * t336;
t650 = rSges(6,3) + qJ(5);
t662 = rSges(6,1) + pkin(4);
t528 = t650 * t549 - t550 * t662;
t334 = t337 ^ 2;
t307 = rSges(3,1) * t340 + rSges(3,2) * t342;
t333 = t336 ^ 2;
t615 = t333 + t334;
t680 = t307 * t615;
t675 = t696 * t336 - t698 * t337;
t674 = t698 * t336 + t696 * t337;
t673 = t697 * t336;
t672 = t697 * t337;
t571 = Icges(3,4) * t342;
t446 = -Icges(3,2) * t340 + t571;
t226 = -Icges(3,6) * t337 + t336 * t446;
t572 = Icges(3,4) * t340;
t452 = Icges(3,1) * t342 - t572;
t228 = -Icges(3,5) * t337 + t336 * t452;
t445 = -Icges(3,2) * t342 - t572;
t451 = -Icges(3,1) * t340 - t571;
t569 = Icges(4,4) * t331;
t570 = Icges(4,4) * t330;
t444 = -Icges(4,2) * t330 + t569;
t198 = -Icges(4,6) * t337 + t336 * t444;
t450 = Icges(4,1) * t331 - t570;
t200 = -Icges(4,5) * t337 + t336 * t450;
t630 = t198 * t331 + t200 * t330;
t631 = qJD(2) * t330;
t670 = ((-Icges(4,2) * t331 - t570) * t631 - (-Icges(4,1) * t330 - t569) * t523 - (-t340 * t445 + t342 * t451) * qJD(2)) * t336 + (t226 * t342 + t228 * t340 + t630) * qJD(2);
t286 = t330 * t514 + t522;
t287 = t336 * t516 - t520;
t515 = qJD(4) * t331;
t669 = t720 * t515 + (t700 * t271 + t699 * t272) * t287 + (t700 * t273 + t699 * t274) * t286;
t667 = t198 * t330 - t200 * t331 + t226 * t340 - t228 * t342;
t666 = 0.2e1 * qJD(2);
t665 = 2 * qJDD(2);
t507 = qJD(2) * qJD(4);
t385 = qJDD(4) * t330 + t331 * t507;
t506 = qJDD(2) * t336;
t172 = t337 * t385 + t506;
t505 = qJDD(2) * t337;
t173 = t336 * t385 - t505;
t279 = -qJDD(4) * t331 + t330 * t507;
t664 = t693 * t172 + t681 * t173 + t691 * t279 + t658 * t286 + t659 * t287 + t695 * t515;
t663 = t653 * t172 + t692 * t173 + t690 * t279 + t656 * t286 + t657 * t287 + t694 * t515;
t628 = ((t430 + t434) * qJD(2) - t719) * t331 + (t711 * t341 + t713 * t339 + (-t339 * t707 + t688 * t341) * qJD(4) + t686 * qJD(2)) * t330;
t627 = ((t429 + t433) * qJD(2) - t718) * t331 + (t710 * t341 + t712 * t339 + (-t339 * t706 + t687 * t341) * qJD(4) + t685 * qJD(2)) * t330;
t661 = t693 * t286 + t681 * t287 - t691 * t515;
t660 = t653 * t286 + t692 * t287 - t690 * t515;
t655 = t651 * t286 + t652 * t287 - t625 * t515;
t363 = -t330 * t435 + t556;
t144 = t363 * t336;
t150 = t182 * t336;
t649 = t144 + t150;
t145 = t363 * t337;
t151 = t182 * t337;
t648 = t145 + t151;
t152 = t184 * t336;
t154 = t186 * t336;
t647 = -t152 - t154;
t153 = t184 * t337;
t155 = t186 * t337;
t646 = -t153 - t155;
t643 = t682 * t336;
t642 = t682 * t337;
t405 = t179 - t426;
t406 = -t181 + t427;
t604 = (t180 * t337 + t433) * t286 + (t180 * t336 + t434) * t287;
t605 = -(-t178 * t337 - t429) * t286 - (-t178 * t336 - t430) * t287;
t641 = (-t604 - t605 + (-t405 + t406) * t515) * t330;
t466 = rSges(6,1) * t341 + rSges(6,3) * t339;
t640 = (pkin(4) * t341 + qJ(5) * t339 + t466) * t330;
t639 = t685 * t286 + t686 * t287;
t638 = -t553 - t554;
t637 = t690 * t330;
t636 = t691 * t330;
t635 = t693 * t337;
t634 = t692 * t336;
t626 = ((-t426 - t427) * qJD(2) - t715) * t331 + (t708 * t341 + t709 * t339 + (t339 * t676 - t341 * t705) * qJD(4) + t704 * qJD(2)) * t330;
t332 = t342 * pkin(2);
t617 = t331 * rSges(4,1) - rSges(4,2) * t330;
t619 = t617 + t332;
t616 = t331 * pkin(3) + t330 * pkin(6);
t614 = (t676 + t714) * t515 + (t272 * t678 + t240 - t562 - t707) * t287 + (t274 * t678 + t241 - t561 - t706) * t286;
t613 = (Icges(6,1) * t550 - t255 - t313 - t705) * t515 + (-t271 * t717 + t238 - t568 + t688) * t287 + (-t273 * t717 + t239 - t567 + t687) * t286;
t612 = t669 * t330;
t457 = t50 * t336 + t51 * t337;
t458 = t48 * t336 + t49 * t337;
t611 = t457 + t458;
t610 = t337 * t653 + t634;
t609 = t336 * t681 + t635;
t608 = g(1) * t337 + g(2) * t336;
t607 = t330 * t608;
t509 = qJD(5) * t339;
t306 = t330 * t509;
t248 = t616 * t336;
t249 = t616 * t337;
t194 = -qJ(3) * t337 + t332 * t336;
t195 = qJ(3) * t336 + t332 * t337;
t492 = t194 * t522 + t195 * t520 + qJD(1);
t453 = t248 * t522 + t249 * t520 + t492;
t539 = rSges(6,2) * t551 + t273 * t650 + t274 * t662;
t540 = rSges(6,2) * t552 + t271 * t650 + t272 * t662;
t27 = t286 * t540 - t287 * t539 + t306 + t453;
t294 = pkin(3) * t330 - pkin(6) * t331;
t404 = qJD(2) * t294;
t222 = t336 * t404;
t223 = t337 * t404;
t582 = pkin(2) * qJD(2);
t502 = t340 * t582;
t517 = qJD(3) * t337;
t290 = -t336 * t502 - t517;
t329 = qJD(3) * t336;
t291 = -t337 * t502 + t329;
t418 = t194 * t506 + t195 * t505 + t290 * t522 + t291 * t520 + qJDD(1);
t368 = -t222 * t522 - t223 * t520 + t248 * t506 + t249 * t505 + t418;
t377 = t330 * t512 + t339 * t523;
t510 = qJD(5) * t273;
t585 = rSges(6,2) * t496 - t162 * t650 + t163 * t662 + t510;
t511 = qJD(5) * t271;
t586 = rSges(6,2) * t497 - t160 * t650 + t161 * t662 + t511;
t5 = qJD(5) * t377 + qJDD(5) * t550 + t172 * t540 - t173 * t539 + t286 * t586 - t287 * t585 + t368;
t606 = t27 * t585 + t5 * t539;
t603 = -t340 * (t445 * t336 + t228) - t342 * (-t451 * t336 + t226);
t343 = qJD(2) ^ 2;
t602 = t172 / 0.2e1;
t601 = t173 / 0.2e1;
t600 = t279 / 0.2e1;
t599 = -t286 / 0.2e1;
t598 = t286 / 0.2e1;
t597 = -t287 / 0.2e1;
t596 = t287 / 0.2e1;
t591 = pkin(2) * t340;
t584 = rSges(5,1) * t341;
t321 = t330 * rSges(6,2);
t320 = t330 * rSges(5,3);
t573 = (-rSges(6,1) * t339 + rSges(6,3) * t341) * t516 + (t331 * t466 + t321) * qJD(2) + t306 + t377 * qJ(5) + (-t330 * t513 + t331 * t518) * pkin(4);
t548 = t331 * t336;
t547 = t331 * t337;
t546 = t331 * t339;
t545 = t331 * t341;
t538 = -t271 * t662 + t272 * t650;
t537 = t273 * t662 - t274 * t650;
t303 = rSges(6,2) * t548;
t536 = -t336 * t640 + t303;
t305 = rSges(6,2) * t547;
t535 = t337 * t640 - t305;
t534 = t336 * t194 + t337 * t195;
t533 = -rSges(6,2) * t331 + t640;
t532 = t545 * t662 + t546 * t650 + t321;
t310 = pkin(6) * t548;
t312 = pkin(6) * t547;
t531 = (-pkin(3) * t552 + t310) * t522 + (-pkin(3) * t551 + t312) * t520;
t527 = t336 * t290 + t337 * t291;
t526 = t330 * rSges(5,2) * t544 + rSges(5,3) * t548;
t525 = rSges(5,2) * t500 + rSges(5,3) * t547;
t508 = -m(4) - m(5) - m(6);
t504 = qJDD(3) * t337;
t503 = rSges(5,1) * t549;
t501 = t342 * t582;
t499 = t332 + t616;
t487 = t520 / 0.2e1;
t485 = -t515 / 0.2e1;
t484 = t515 / 0.2e1;
t292 = rSges(4,1) * t330 + rSges(4,2) * t331;
t378 = -t292 - t591;
t482 = -t294 - t591;
t481 = t615 * t340;
t480 = t336 * t248 + t337 * t249 + t534;
t479 = -t336 * t222 - t337 * t223 + t527;
t474 = -t336 * t591 + t310;
t473 = -t337 * t591 + t312;
t467 = -rSges(5,2) * t339 + t584;
t191 = -rSges(5,3) * t331 + t330 * t467;
t472 = -t191 + t482;
t277 = t617 * qJD(2);
t471 = -t277 - t501;
t278 = t616 * qJD(2);
t470 = -t278 - t501;
t308 = rSges(3,1) * t342 - rSges(3,2) * t340;
t260 = (-rSges(5,1) * t339 - rSges(5,2) * t341) * t330;
t107 = qJD(4) * t260 + (t331 * t467 + t320) * qJD(2);
t123 = rSges(5,1) * t274 - rSges(5,2) * t273 + rSges(5,3) * t551;
t390 = (-qJDD(2) * t340 - t342 * t343) * pkin(2);
t354 = -qJD(2) * t278 - qJDD(2) * t294 + t390;
t348 = t336 * t354 - t504;
t91 = rSges(5,1) * t163 + rSges(5,2) * t162 + rSges(5,3) * t496;
t34 = -t107 * t286 + t123 * t279 - t172 * t191 - t515 * t91 + t348;
t121 = rSges(5,1) * t272 - rSges(5,2) * t271 + rSges(5,3) * t552;
t328 = qJDD(3) * t336;
t349 = t337 * t354 + t328;
t89 = rSges(5,1) * t161 + rSges(5,2) * t160 + rSges(5,3) * t497;
t35 = t107 * t287 - t121 * t279 + t173 * t191 + t515 * t89 + t349;
t464 = t336 * t35 - t337 * t34;
t454 = qJD(2) * t482;
t369 = t336 * t454 - t517;
t37 = -t286 * t533 - t515 * t539 + t369 + t511;
t376 = t337 * t454 + t329;
t38 = t287 * t533 + t515 * t540 + t376 + t510;
t463 = -t336 * t37 - t337 * t38;
t56 = -t123 * t515 - t191 * t286 + t369;
t57 = t121 * t515 + t191 * t287 + t376;
t456 = -t336 * t56 - t337 * t57;
t455 = qJD(2) * t378;
t428 = t121 * t337 - t123 * t336;
t206 = -rSges(4,3) * t337 + t336 * t617;
t207 = rSges(4,3) * t336 + t337 * t617;
t423 = t206 * t336 + t207 * t337;
t420 = t615 * t308;
t419 = qJD(2) * t680;
t417 = t482 - t533;
t416 = rSges(5,1) * t545 - rSges(5,2) * t546 + t320;
t415 = -t107 + t470;
t66 = qJD(2) * t423 + t492;
t403 = t66 * t292;
t400 = t470 - t573;
t396 = qJD(2) * t292;
t36 = t121 * t286 - t123 * t287 + t453;
t395 = t36 * t428;
t375 = t27 * t586 + t5 * t540;
t374 = t37 * t539 - t38 * t540;
t373 = t439 * t336 + ((-t446 + t451) * t342 + (-t445 - t452) * t340) * t337;
t355 = -qJD(2) * t277 - qJDD(2) * t292 + t390;
t347 = (t27 * t540 - t37 * t533) * t337 + (-t27 * t539 + t38 * t533) * t336;
t346 = t336 * (-(Icges(4,6) * t336 + t337 * t444) * t331 - (Icges(4,5) * t336 + t337 * t450) * t330) + t337 * t630;
t289 = t307 * t337;
t288 = t307 * t336;
t221 = t337 * t396;
t220 = t336 * t396;
t165 = t337 * t455 + t329;
t164 = t336 * t455 - t517;
t159 = -t337 * t503 + t525;
t157 = -t336 * t503 + t526;
t142 = -rSges(5,1) * t273 - rSges(5,2) * t274;
t138 = -rSges(5,1) * t271 - rSges(5,2) * t272;
t97 = t337 * t355 + t328;
t96 = t336 * t355 - t504;
t69 = -qJD(2) * t419 + qJDD(2) * t420 + qJDD(1);
t39 = t423 * qJDD(2) + (-t220 * t336 - t221 * t337) * qJD(2) + t418;
t16 = t121 * t172 - t123 * t173 + t286 * t89 - t287 * t91 + t368;
t7 = -qJD(5) * t162 + qJDD(5) * t273 + t173 * t533 - t279 * t540 + t287 * t573 + t515 * t586 + t349;
t6 = -qJD(5) * t160 + qJDD(5) * t271 - t172 * t533 + t279 * t539 - t286 * t573 - t515 * t585 + t348;
t1 = [m(2) * qJDD(1) + (-m(2) - m(3) + t508) * g(3) + m(3) * t69 + m(4) * t39 + m(5) * t16 + m(6) * t5; (t653 * t336 - t692 * t337) * t602 + (t693 * t336 - t337 * t681) * t601 + (t336 * t651 - t337 * t652) * t600 + (((t273 * t645 + t274 * t644 + t634) * t331 + t637) * qJD(4) + (((t638 + t653) * qJD(4) + t639) * t331 + t641) * t337 + (t273 * t649 + t274 * t647) * t287 + (t273 * t648 + t274 * t646) * t286) * t599 + (t336 * t656 - t337 * t657) * t598 + (((t271 * t645 + t272 * t644 + t635) * t331 + t636) * qJD(4) + (((t638 + t681) * qJD(4) + t639) * t331 + t641) * t336 + (t271 * t649 + t272 * t647) * t287 + (t271 * t648 + t272 * t646) * t286) * t597 + (t336 * t658 - t337 * t659) * t596 - (t642 * qJD(2) * t333 + (-t603 * t337 + t346 + (t373 - t643) * t336) * t520) * t522 / 0.2e1 + ((t373 * t336 + t346 + (-t603 - t642) * t337) * t522 + t643 * qJD(2) * t334) * t487 - t655 * t516 / 0.2e1 + (((t145 * t339 - t153 * t341 + t113) * t286 + (t144 * t339 - t152 * t341 + t112) * t287 + t67 * qJD(4)) * t330 + ((-t406 * t331 + (-t177 * t339 - t185 * t341 - t180) * t330 + t458) * qJD(4) + t604) * t331 + ((t151 * t339 - t155 * t341 + t111) * t286 + (t150 * t339 - t154 * t341 + t110) * t287 + t68 * qJD(4)) * t330 + ((t405 * t331 + (t183 * t339 - t187 * t341 - t178) * t330 + t457) * qJD(4) + t605) * t331) * t484 + (-t27 * t531 - (t27 * t331 + t330 * t463) * t509 - (t27 * t535 + t38 * t532) * t287 - (t27 * t536 - t37 * t532) * t286 - (t463 * t616 + (-t27 * t481 + t342 * t463) * pkin(2)) * qJD(2) - (t374 * t330 + (t37 * t535 + t38 * t536 + t347) * t331) * qJD(4) + t5 * t480 + t27 * t479 + (t38 * t400 + t417 * t7 + t606) * t337 + (t37 * t400 + t417 * t6 + t375) * t336 - g(1) * (t305 + t473) - g(2) * (t303 + t474) - g(3) * (t499 + t532) - (-t339 * t650 - t341 * t662 - pkin(3)) * t607) * m(6) + (-t36 * (t157 * t286 - t159 * t287 + t531) - (-t286 * t56 + t287 * t57) * t416 - (t456 * t616 + (t342 * t456 - t36 * t481) * pkin(2)) * qJD(2) - ((-t121 * t57 + t123 * t56) * t330 + (t57 * (t191 * t336 + t157) + t56 * (-t191 * t337 - t159) + t395) * t331) * qJD(4) + t16 * t480 + t36 * t479 + (t16 * t123 + t35 * t472 + t36 * t91 + t415 * t57) * t337 + (t16 * t121 + t34 * t472 + t36 * t89 + t415 * t56) * t336 - g(1) * (t473 + t525) - g(2) * (t474 + t526) - g(3) * (t416 + t499) - (-pkin(3) - t584) * t607) * m(5) + (-(-t66 * pkin(2) * t481 + (-t165 * t619 - t337 * t403) * t337 + (-t164 * t619 - t336 * t403) * t336) * qJD(2) + t39 * t534 + t66 * t527 + (t165 * t471 + t39 * t207 - t66 * t221 + t378 * t97) * t337 + (t164 * t471 + t39 * t206 - t66 * t220 + t378 * t96) * t336 - g(3) * t619 - t378 * t608) * m(4) + (g(1) * t289 + g(2) * t288 - g(3) * t308 + t69 * t420 + (-t419 - (-t288 * t336 - t289 * t337) * qJD(2)) * (qJD(2) * t420 + qJD(1)) + (qJDD(2) * t307 + t308 * t343) * t680) * m(3) + (t663 + (t670 * t334 + (t672 * t336 - t337 * t673) * t336) * t666 + (t667 * t334 + (t674 * t336 - t337 * t675) * t336) * t665) * t336 / 0.2e1 - (t664 + (t673 * t334 + (t670 - t672) * t689) * t666 + (t675 * t334 + (t667 - t674) * t689) * t665) * t337 / 0.2e1 + ((-t628 + t660) * t337 + (t627 + t661) * t336) * t485; t508 * (g(1) * t336 - g(2) * t337) + m(4) * (t336 * t97 - t337 * t96) + m(5) * t464 + m(6) * (t336 * t7 - t337 * t6); (t610 * t330 - t690 * t331) * t602 + (t609 * t330 - t691 * t331) * t601 + (t330 * t611 - t331 * t625) * t600 + (t273 * t614 + t274 * t613 - t337 * t612) * t599 + (t694 * t331 + (t336 * t657 + t337 * t656) * t330 + (t331 * t610 + t637) * qJD(2)) * t598 + (t271 * t614 + t272 * t613 - t336 * t612) * t597 + (t695 * t331 + (t336 * t659 + t337 * t658) * t330 + (t331 * t609 + t636) * qJD(2)) * t596 - (t172 * t651 + t173 * t652 + t625 * t279 + t627 * t286 + t628 * t287 + t626 * t515) * t331 / 0.2e1 + t664 * t552 / 0.2e1 + t663 * t551 / 0.2e1 + t655 * t631 / 0.2e1 + (t626 * t331 + (t336 * t628 + t337 * t627) * t330 + (t330 * t625 + t331 * t611) * qJD(2)) * t485 + (t669 * t331 + (t339 * t614 + t341 * t613) * t330) * t484 + t661 * t497 / 0.2e1 + t660 * t331 * t487 + ((qJD(2) * t347 - t37 * t585 + t38 * t586 - t539 * t6 + t540 * t7) * t331 + (t374 * qJD(2) + (-t37 * t573 - t533 * t6 + t375) * t337 + (t38 * t573 + t533 * t7 - t606) * t336) * t330 - (t27 * t549 + t272 * t37 + t274 * t38) * qJD(5) - (t27 * t537 + t38 * t528) * t287 - (t27 * t538 - t37 * t528) * t286 - (t37 * t537 + t38 * t538) * t515 + g(1) * t537 - g(2) * t538 - g(3) * t528) * m(6) + ((t35 * t121 - t34 * t123 - t56 * t91 + t57 * t89 + (t395 + (t336 * t57 - t337 * t56) * t191) * qJD(2)) * t331 + (t57 * (-qJD(2) * t121 + t107 * t336) + t56 * (qJD(2) * t123 - t107 * t337) + t16 * t428 + t36 * (-t336 * t91 + t337 * t89) + t464 * t191) * t330 - t57 * (t138 * t515 + t260 * t287) - t56 * (-t142 * t515 - t260 * t286) - t36 * (t138 * t286 - t142 * t287) - g(1) * t142 - g(2) * t138 - g(3) * t260) * m(5); (-t160 * t37 - t162 * t38 + (t37 * t515 - g(1) + t7) * t273 + (t286 * t37 - t287 * t38 - g(3) + t5) * t550 + (-t38 * t515 - g(2) + t6) * t271 + (-t271 * t286 + t273 * t287 + t377) * t27) * m(6);];
tau = t1;
