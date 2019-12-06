% Calculate vector of inverse dynamics joint torques for
% S5PRPRP3
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
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPRP3_invdynJ_fixb_slag_vp1(qJ, qJD, qJDD, g, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_slag_vp1: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_invdynJ_fixb_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_invdynJ_fixb_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRP3_invdynJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRP3_invdynJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRP3_invdynJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:32:39
% EndTime: 2019-12-05 15:33:40
% DurationCPUTime: 46.84s
% Computational Cost: add. (15958->757), mult. (23293->1045), div. (0->0), fcn. (22052->8), ass. (0->399)
t710 = Icges(5,5) + Icges(6,5);
t708 = Icges(5,6) + Icges(6,6);
t687 = Icges(5,1) + Icges(6,1);
t729 = Icges(5,2) + Icges(6,2);
t734 = Icges(5,3) + Icges(6,3);
t328 = qJ(2) + pkin(8);
t323 = sin(t328);
t333 = sin(qJ(4));
t335 = cos(qJ(4));
t733 = (-t710 * t333 - t708 * t335) * t323;
t324 = cos(t328);
t330 = cos(pkin(7));
t545 = t330 * t333;
t329 = sin(pkin(7));
t546 = t329 * t335;
t269 = t324 * t546 - t545;
t520 = qJD(2) * t323;
t495 = t333 * t520;
t160 = -qJD(4) * t269 + t329 * t495;
t544 = t330 * t335;
t547 = t329 * t333;
t268 = -t324 * t547 - t544;
t494 = t335 * t520;
t161 = qJD(4) * t268 - t329 * t494;
t518 = qJD(2) * t329;
t493 = t324 * t518;
t732 = t708 * t160 + t710 * t161 + t734 * t493;
t271 = t324 * t544 + t547;
t162 = -qJD(4) * t271 + t330 * t495;
t646 = t324 * t545 - t546;
t163 = -qJD(4) * t646 - t330 * t494;
t516 = qJD(2) * t330;
t492 = t324 * t516;
t731 = t708 * t162 + t710 * t163 + t734 * t492;
t730 = Icges(5,4) + Icges(6,4);
t430 = Icges(6,5) * t335 - Icges(6,6) * t333;
t179 = Icges(6,3) * t323 + t324 * t430;
t431 = Icges(5,5) * t335 - Icges(5,6) * t333;
t181 = Icges(5,3) * t323 + t324 * t431;
t728 = -t733 * qJD(4) + (-t179 - t181) * qJD(2);
t554 = t323 * t329;
t110 = Icges(6,5) * t269 + Icges(6,6) * t268 + Icges(6,3) * t554;
t112 = Icges(5,5) * t269 + Icges(5,6) * t268 + Icges(5,3) * t554;
t697 = t112 + t110;
t553 = t323 * t330;
t111 = Icges(6,5) * t271 - Icges(6,6) * t646 + Icges(6,3) * t553;
t113 = Icges(5,5) * t271 - Icges(5,6) * t646 + Icges(5,3) * t553;
t696 = t113 + t111;
t178 = -Icges(6,3) * t324 + t323 * t430;
t180 = -Icges(5,3) * t324 + t323 * t431;
t714 = t178 + t180;
t563 = Icges(6,4) * t335;
t436 = -Icges(6,2) * t333 + t563;
t183 = Icges(6,6) * t323 + t324 * t436;
t567 = Icges(5,4) * t335;
t437 = -Icges(5,2) * t333 + t567;
t185 = Icges(5,6) * t323 + t324 * t437;
t727 = t183 + t185;
t564 = Icges(6,4) * t333;
t442 = Icges(6,1) * t335 - t564;
t187 = Icges(6,5) * t323 + t324 * t442;
t568 = Icges(5,4) * t333;
t443 = Icges(5,1) * t335 - t568;
t189 = Icges(5,5) * t323 + t324 * t443;
t652 = -t187 - t189;
t726 = (-t729 * t335 - t564 - t568) * t323;
t725 = (t687 * t333 + t563 + t567) * t323;
t334 = sin(qJ(2));
t336 = cos(qJ(2));
t434 = -Icges(3,5) * t334 - Icges(3,6) * t336;
t690 = -Icges(4,5) * t323 - Icges(4,6) * t324 + t434;
t724 = t729 * t160 + t730 * t161 + t708 * t493;
t723 = t729 * t162 + t730 * t163 + t708 * t492;
t722 = t730 * t160 + t687 * t161 + t710 * t493;
t721 = t730 * t162 + t687 * t163 + t710 * t492;
t720 = t727 * qJD(2) + t726 * qJD(4);
t719 = t652 * qJD(2) + t725 * qJD(4);
t566 = Icges(6,4) * t269;
t114 = Icges(6,2) * t268 + Icges(6,6) * t554 + t566;
t570 = Icges(5,4) * t269;
t116 = Icges(5,2) * t268 + Icges(5,6) * t554 + t570;
t718 = t114 + t116;
t565 = Icges(6,4) * t271;
t115 = -Icges(6,2) * t646 + Icges(6,6) * t553 + t565;
t569 = Icges(5,4) * t271;
t117 = -Icges(5,2) * t646 + Icges(5,6) * t553 + t569;
t717 = t115 + t117;
t238 = Icges(6,4) * t268;
t118 = Icges(6,1) * t269 + Icges(6,5) * t554 + t238;
t240 = Icges(5,4) * t268;
t120 = Icges(5,1) * t269 + Icges(5,5) * t554 + t240;
t716 = t118 + t120;
t239 = Icges(6,4) * t646;
t119 = Icges(6,1) * t271 + Icges(6,5) * t553 - t239;
t241 = Icges(5,4) * t646;
t121 = Icges(5,1) * t271 + Icges(5,5) * t553 - t241;
t715 = t119 + t121;
t182 = -Icges(6,6) * t324 + t323 * t436;
t184 = -Icges(5,6) * t324 + t323 * t437;
t686 = t182 + t184;
t186 = -Icges(6,5) * t324 + t323 * t442;
t188 = -Icges(5,5) * t324 + t323 * t443;
t691 = t186 + t188;
t519 = qJD(2) * t324;
t713 = t728 * t323 - t714 * t519;
t712 = t731 * t323 + t696 * t519;
t711 = t732 * t323 + t697 * t519;
t669 = rSges(6,1) + pkin(4);
t707 = Icges(3,3) + Icges(4,3);
t706 = t690 * qJD(2);
t705 = Icges(3,5) * t336 + Icges(4,5) * t324 - Icges(3,6) * t334 - Icges(4,6) * t323;
t666 = t718 * t160 + t716 * t161 + t724 * t268 + t722 * t269 + t711 * t329;
t665 = t717 * t160 + t715 * t161 + t723 * t268 + t721 * t269 + t712 * t329;
t664 = t718 * t162 + t716 * t163 + t722 * t271 + t711 * t330 - t724 * t646;
t663 = t717 * t162 + t715 * t163 + t721 * t271 + t712 * t330 - t723 * t646;
t704 = -t686 * t160 - t691 * t161 - t720 * t268 + t719 * t269 + t713 * t329;
t703 = -t686 * t162 - t691 * t163 + t719 * t271 + t713 * t330 + t720 * t646;
t689 = t718 * t268 + t716 * t269 + t697 * t554;
t702 = t717 * t268 + t715 * t269 + t696 * t554;
t701 = t716 * t271 + t697 * t553 - t718 * t646;
t660 = t715 * t271 + t696 * t553 - t717 * t646;
t427 = -t114 * t333 + t118 * t335;
t48 = -t110 * t324 + t323 * t427;
t425 = -t116 * t333 + t120 * t335;
t50 = -t112 * t324 + t323 * t425;
t659 = t48 + t50;
t426 = -t115 * t333 + t119 * t335;
t49 = -t111 * t324 + t323 * t426;
t424 = -t117 * t333 + t121 * t335;
t51 = -t113 * t324 + t323 * t424;
t658 = t49 + t51;
t700 = t686 * t268 + t691 * t269 + t714 * t554;
t699 = t691 * t271 + t714 * t553 - t686 * t646;
t422 = -t182 * t333 + t186 * t335;
t556 = t178 * t324;
t67 = t323 * t422 - t556;
t421 = -t184 * t333 + t188 * t335;
t555 = t180 * t324;
t68 = t323 * t421 - t555;
t633 = t67 + t68;
t698 = t329 * t330;
t513 = qJD(4) * t323;
t490 = t330 * t513;
t283 = t490 + t518;
t491 = t329 * t513;
t284 = t491 - t516;
t512 = qJD(4) * t324;
t621 = t283 * (-t271 * t729 - t239 - t241 + t715) + t284 * (-t269 * t729 + t238 + t240 + t716) - t512 * (t691 + t726);
t327 = t330 ^ 2;
t306 = rSges(3,1) * t334 + rSges(3,2) * t336;
t326 = t329 ^ 2;
t622 = t326 + t327;
t688 = t306 * t622;
t685 = t329 * t705 - t330 * t707;
t684 = t329 * t707 + t330 * t705;
t683 = t706 * t329;
t682 = t706 * t330;
t573 = Icges(3,4) * t336;
t441 = -Icges(3,2) * t334 + t573;
t226 = -Icges(3,6) * t330 + t329 * t441;
t574 = Icges(3,4) * t334;
t447 = Icges(3,1) * t336 - t574;
t228 = -Icges(3,5) * t330 + t329 * t447;
t440 = -Icges(3,2) * t336 - t574;
t446 = -Icges(3,1) * t334 - t573;
t571 = Icges(4,4) * t324;
t572 = Icges(4,4) * t323;
t439 = -Icges(4,2) * t323 + t571;
t200 = -Icges(4,6) * t330 + t329 * t439;
t445 = Icges(4,1) * t324 - t572;
t202 = -Icges(4,5) * t330 + t329 * t445;
t638 = t200 * t324 + t202 * t323;
t680 = ((-Icges(4,2) * t324 - t572) * t520 - (-Icges(4,1) * t323 - t571) * t519 - (-t334 * t440 + t336 * t446) * qJD(2)) * t329 + (t226 * t336 + t228 * t334 + t638) * qJD(2);
t679 = t733 * t512 + (-t268 * t710 + t269 * t708) * t284 + (t271 * t708 + t646 * t710) * t283;
t677 = t200 * t323 - t202 * t324 + t226 * t334 - t228 * t336;
t673 = 0.2e1 * qJD(2);
t672 = 2 * qJDD(2);
t508 = qJD(2) * qJD(4);
t381 = qJDD(4) * t323 + t324 * t508;
t507 = qJDD(2) * t329;
t174 = t330 * t381 + t507;
t506 = qJDD(2) * t330;
t175 = t329 * t381 - t506;
t276 = -qJDD(4) * t324 + t323 * t508;
t671 = t702 * t174 + t689 * t175 + t700 * t276 + t665 * t283 + t284 * t666 + t704 * t512;
t670 = t174 * t660 + t175 * t701 + t276 * t699 + t283 * t663 + t284 * t664 + t512 * t703;
t636 = ((t425 + t427) * qJD(2) - t732) * t324 + (t722 * t335 - t724 * t333 + (-t333 * t716 - t335 * t718) * qJD(4) + t697 * qJD(2)) * t323;
t635 = ((t424 + t426) * qJD(2) - t731) * t324 + (t721 * t335 - t723 * t333 + (-t333 * t715 - t335 * t717) * qJD(4) + t696 * qJD(2)) * t323;
t668 = t283 * t702 + t284 * t689 - t512 * t700;
t667 = t283 * t660 + t284 * t701 - t512 * t699;
t662 = t283 * t658 + t284 * t659 - t512 * t633;
t148 = t182 * t329;
t150 = t184 * t329;
t657 = -t148 - t150;
t149 = t182 * t330;
t151 = t184 * t330;
t656 = -t149 - t151;
t152 = t186 * t329;
t154 = t188 * t329;
t655 = -t152 - t154;
t153 = t186 * t330;
t155 = t188 * t330;
t654 = -t153 - t155;
t651 = t690 * t329;
t650 = t690 * t330;
t400 = t181 - t421;
t401 = t179 - t422;
t612 = -(-t178 * t330 - t426) * t283 - (-t178 * t329 - t427) * t284;
t613 = -(-t180 * t330 - t424) * t283 - (-t180 * t329 - t425) * t284;
t649 = (-t612 - t613 + (-t400 - t401) * t512) * t323;
t648 = t283 * t696 + t284 * t697;
t319 = pkin(4) * t335 + pkin(3);
t596 = pkin(3) - t319;
t484 = t596 * t323;
t331 = -qJ(5) - pkin(6);
t594 = pkin(6) + t331;
t647 = -t324 * t594 + t484;
t645 = -t555 - t556;
t644 = t699 * t323;
t643 = t700 * t323;
t642 = t702 * t330;
t641 = t701 * t329;
t634 = ((-t421 - t422) * qJD(2) - t728) * t324 + (t719 * t335 + t720 * t333 + (t333 * t691 + t686 * t335) * qJD(4) - t714 * qJD(2)) * t323;
t313 = t323 * rSges(6,3);
t548 = t324 * t335;
t549 = t324 * t333;
t626 = rSges(6,1) * t548 - rSges(6,2) * t549 + t324 * t319 - t323 * t331 + t313;
t325 = t336 * pkin(2);
t624 = t324 * rSges(4,1) - rSges(4,2) * t323;
t625 = t624 + t325;
t623 = t324 * pkin(3) + t323 * pkin(6);
t620 = (t686 + t725) * t512 + (t268 * t687 - t566 - t570 - t718) * t284 + (-t646 * t687 - t565 - t569 - t717) * t283;
t619 = t679 * t323;
t452 = t50 * t329 + t51 * t330;
t453 = t48 * t329 + t49 * t330;
t618 = t452 + t453;
t617 = t330 * t660 + t641;
t616 = t329 * t689 + t642;
t615 = g(1) * t330 + g(2) * t329;
t248 = t623 * t329;
t249 = t623 * t330;
t196 = -qJ(3) * t330 + t325 * t329;
t197 = qJ(3) * t329 + t325 * t330;
t487 = t196 * t518 + t197 * t516 + qJD(1);
t448 = t248 * t518 + t249 * t516 + t487;
t510 = qJD(5) * t324;
t371 = -t323 * t594 - t324 * t596;
t543 = rSges(6,1) * t271 - rSges(6,2) * t646 + rSges(6,3) * t553 + pkin(4) * t547 + t330 * t371;
t575 = rSges(6,1) * t269 + rSges(6,2) * t268 + rSges(6,3) * t554 - pkin(4) * t545 + t329 * t371;
t27 = t283 * t575 - t284 * t543 + t448 - t510;
t291 = pkin(3) * t323 - pkin(6) * t324;
t399 = qJD(2) * t291;
t222 = t329 * t399;
t223 = t330 * t399;
t587 = pkin(2) * qJD(2);
t504 = t334 * t587;
t514 = qJD(3) * t330;
t287 = -t329 * t504 - t514;
t322 = qJD(3) * t329;
t288 = -t330 * t504 + t322;
t413 = t196 * t507 + t197 * t506 + t287 * t518 + t288 * t516 + qJDD(1);
t365 = -t222 * t518 - t223 * t516 + t248 * t507 + t249 * t506 + t413;
t511 = qJD(5) * t323;
t300 = t330 * t511;
t586 = pkin(4) * qJD(4);
t502 = t333 * t586;
t343 = qJD(2) * t647 - t324 * t502;
t501 = t335 * t586;
t592 = rSges(6,1) * t163 + rSges(6,2) * t162 + rSges(6,3) * t492 + t329 * t501 + t330 * t343 + t300;
t299 = t329 * t511;
t593 = rSges(6,1) * t161 + rSges(6,2) * t160 + rSges(6,3) * t493 + t329 * t343 - t330 * t501 + t299;
t5 = qJD(2) * t511 - qJDD(5) * t324 + t174 * t575 - t175 * t543 + t283 * t593 - t284 * t592 + t365;
t614 = t27 * t592 + t5 * t543;
t611 = -t334 * (t440 * t329 + t228) - t336 * (-t446 * t329 + t226);
t337 = qJD(2) ^ 2;
t610 = t174 / 0.2e1;
t609 = t175 / 0.2e1;
t608 = t276 / 0.2e1;
t607 = -t283 / 0.2e1;
t606 = t283 / 0.2e1;
t605 = -t284 / 0.2e1;
t604 = t284 / 0.2e1;
t600 = pkin(2) * t334;
t591 = rSges(5,1) * t335;
t590 = rSges(6,1) * t335;
t588 = rSges(6,2) * t335;
t585 = t27 * t323;
t314 = t323 * rSges(5,3);
t258 = (-rSges(6,1) * t333 - t588) * t323;
t460 = -rSges(6,2) * t333 + t590;
t576 = qJD(4) * t258 - t323 * t502 - t510 + (t324 * t460 + t313 + t371) * qJD(2);
t552 = t324 * t329;
t551 = t324 * t330;
t550 = t324 * t331;
t308 = pkin(6) * t552;
t394 = t484 - t550;
t499 = t323 * t546;
t500 = t323 * t547;
t525 = rSges(6,2) * t500 + rSges(6,3) * t552;
t538 = -rSges(6,1) * t499 + t329 * t394 - t308 + t525;
t309 = pkin(6) * t551;
t497 = t323 * t544;
t498 = t323 * t545;
t523 = rSges(6,2) * t498 + rSges(6,3) * t551;
t537 = rSges(6,1) * t497 - t330 * t394 + t309 - t523;
t536 = -t269 * rSges(6,2) + t268 * t669;
t535 = t271 * rSges(6,2) + t646 * t669;
t534 = -rSges(6,3) * t324 + t323 * t460 - t647;
t533 = -t623 + t626;
t532 = t329 * t196 + t330 * t197;
t529 = (-pkin(3) * t554 + t308) * t518 + (-pkin(3) * t553 + t309) * t516;
t526 = t329 * t287 + t330 * t288;
t524 = rSges(5,2) * t500 + rSges(5,3) * t552;
t522 = rSges(5,2) * t498 + rSges(5,3) * t551;
t509 = -m(4) - m(5) - m(6);
t505 = qJDD(3) * t330;
t503 = t336 * t587;
t481 = t516 / 0.2e1;
t479 = -t512 / 0.2e1;
t478 = t512 / 0.2e1;
t289 = rSges(4,1) * t323 + rSges(4,2) * t324;
t374 = -t289 - t600;
t476 = -t291 - t600;
t475 = t622 * t334;
t473 = t329 * t248 + t330 * t249 + t532;
t472 = -t329 * t222 - t330 * t223 + t526;
t467 = pkin(4) * t323 * t333 - t258;
t461 = -rSges(5,2) * t333 + t591;
t193 = -rSges(5,3) * t324 + t323 * t461;
t466 = -t193 + t476;
t274 = t624 * qJD(2);
t465 = -t274 - t503;
t275 = t623 * qJD(2);
t464 = -t275 - t503;
t307 = rSges(3,1) * t336 - rSges(3,2) * t334;
t259 = (-rSges(5,1) * t333 - rSges(5,2) * t335) * t323;
t109 = qJD(4) * t259 + (t324 * t461 + t314) * qJD(2);
t125 = rSges(5,1) * t271 - rSges(5,2) * t646 + rSges(5,3) * t553;
t386 = (-qJDD(2) * t334 - t336 * t337) * pkin(2);
t357 = -qJDD(2) * t291 + t386;
t348 = -qJD(2) * t275 + t357;
t91 = rSges(5,1) * t163 + rSges(5,2) * t162 + rSges(5,3) * t492;
t34 = -t109 * t283 + t125 * t276 - t174 * t193 + t329 * t348 - t512 * t91 - t505;
t123 = rSges(5,1) * t269 + rSges(5,2) * t268 + rSges(5,3) * t554;
t321 = qJDD(3) * t329;
t89 = rSges(5,1) * t161 + rSges(5,2) * t160 + rSges(5,3) * t493;
t35 = t109 * t284 - t123 * t276 + t175 * t193 + t330 * t348 + t512 * t89 + t321;
t459 = t329 * t35 - t330 * t34;
t449 = qJD(2) * t476;
t366 = t329 * t449 - t514;
t36 = -t283 * t534 - t512 * t543 + t299 + t366;
t372 = t330 * t449 + t322;
t37 = t284 * t534 + t512 * t575 + t300 + t372;
t458 = t329 * t36 + t330 * t37;
t56 = -t125 * t512 - t193 * t283 + t366;
t57 = t123 * t512 + t193 * t284 + t372;
t451 = -t329 * t56 - t330 * t57;
t450 = qJD(2) * t374;
t423 = t123 * t330 - t125 * t329;
t206 = -rSges(4,3) * t330 + t329 * t624;
t207 = rSges(4,3) * t329 + t330 * t624;
t418 = t206 * t329 + t207 * t330;
t415 = t622 * t307;
t414 = qJD(2) * t688;
t412 = t476 - t534;
t411 = rSges(5,1) * t548 - rSges(5,2) * t549 + t314;
t410 = -t109 + t464;
t66 = qJD(2) * t418 + t487;
t398 = t66 * t289;
t395 = t464 - t576;
t392 = qJD(2) * t289;
t38 = t123 * t283 - t125 * t284 + t448;
t391 = t38 * t423;
t373 = t27 * t593 + t5 * t575;
t370 = t36 * t543 - t37 * t575;
t369 = t434 * t329 + ((-t441 + t446) * t336 + (-t440 - t447) * t334) * t330;
t364 = t324 * t458 + t585;
t349 = -qJD(2) * t274 - qJDD(2) * t289 + t386;
t342 = qJDD(5) * t323 + (-t275 + t510) * qJD(2) + t357;
t341 = (t27 * t575 - t36 * t534) * t330 + (-t27 * t543 + t37 * t534) * t329;
t340 = t329 * (-(Icges(4,6) * t329 + t330 * t439) * t324 - (Icges(4,5) * t329 + t330 * t445) * t323) + t330 * t638;
t286 = t306 * t330;
t285 = t306 * t329;
t221 = t330 * t392;
t220 = t329 * t392;
t167 = t330 * t450 + t322;
t166 = t329 * t450 - t514;
t159 = -rSges(5,1) * t497 + t522;
t157 = -rSges(5,1) * t499 + t524;
t143 = -rSges(5,1) * t646 - rSges(5,2) * t271;
t141 = rSges(5,1) * t268 - rSges(5,2) * t269;
t97 = t330 * t349 + t321;
t96 = t329 * t349 - t505;
t69 = -qJD(2) * t414 + qJDD(2) * t415 + qJDD(1);
t39 = t418 * qJDD(2) + (-t220 * t329 - t221 * t330) * qJD(2) + t413;
t16 = t123 * t174 - t125 * t175 + t283 * t89 - t284 * t91 + t365;
t7 = t175 * t534 - t276 * t575 + t284 * t576 + t330 * t342 + t512 * t593 + t321;
t6 = -t174 * t534 + t276 * t543 - t283 * t576 + t329 * t342 - t512 * t592 - t505;
t1 = [m(2) * qJDD(1) + (-m(2) - m(3) + t509) * g(3) + m(3) * t69 + m(4) * t39 + m(5) * t16 + m(6) * t5; (t660 * t329 - t330 * t701) * t610 + (t329 * t702 - t330 * t689) * t609 + (t329 * t658 - t330 * t659) * t608 + (((t271 * t652 + t646 * t727 + t641) * t324 + t644) * qJD(4) + (((t645 + t660) * qJD(4) + t648) * t324 + t649) * t330 + (t271 * t655 - t646 * t657) * t284 + (t271 * t654 - t646 * t656) * t283) * t607 + (t329 * t663 - t330 * t664) * t606 + (((-t268 * t727 + t269 * t652 + t642) * t324 + t643) * qJD(4) + (((t645 + t689) * qJD(4) + t648) * t324 + t649) * t329 + (t268 * t657 + t269 * t655) * t284 + (t268 * t656 + t269 * t654) * t283) * t605 + (t329 * t665 - t330 * t666) * t604 - (t650 * qJD(2) * t326 + (-t611 * t330 + t340 + (t369 - t651) * t329) * t516) * t518 / 0.2e1 + ((t369 * t329 + t340 + (-t611 - t650) * t330) * t518 + t651 * qJD(2) * t327) * t481 - t662 * t513 / 0.2e1 + (((t151 * t333 - t155 * t335 + t113) * t283 + (t150 * t333 - t154 * t335 + t112) * t284 + t68 * qJD(4)) * t323 + ((t400 * t324 + (t185 * t333 - t189 * t335 - t180) * t323 + t452) * qJD(4) + t613) * t324 + ((t149 * t333 - t153 * t335 + t111) * t283 + (t148 * t333 - t152 * t335 + t110) * t284 + t67 * qJD(4)) * t323 + ((t401 * t324 + (t183 * t333 - t187 * t335 - t178) * t323 + t453) * qJD(4) + t612) * t324) * t478 + (-g(1) * t523 - g(2) * t525 - g(3) * (t325 + t626) - t615 * (-t600 - t550 + (-t319 - t590) * t323) + t5 * t473 + t27 * t472 + (t37 * t395 + t412 * t7 + t614) * t330 + (t36 * t395 + t412 * t6 + t373) * t329 - t27 * t529 - t364 * qJD(5) - (t27 * t537 + t37 * t533) * t284 - (t27 * t538 - t36 * t533) * t283 - (-t458 * t623 + (-t27 * t475 - t336 * t458) * pkin(2)) * qJD(2) - (t370 * t323 + (t36 * t537 + t37 * t538 + t341) * t324) * qJD(4)) * m(6) + (t16 * t473 + t38 * t472 + (t16 * t125 + t35 * t466 + t38 * t91 + t410 * t57) * t330 + (t16 * t123 + t34 * t466 + t38 * t89 + t410 * t56) * t329 - g(1) * (-t330 * t600 + t309 + t522) - g(2) * (-t329 * t600 + t308 + t524) - g(3) * (t325 + t411 + t623) - t615 * t323 * (-pkin(3) - t591) - t38 * (t157 * t283 - t159 * t284 + t529) - (-t283 * t56 + t284 * t57) * t411 - (t451 * t623 + (t336 * t451 - t38 * t475) * pkin(2)) * qJD(2) - ((-t123 * t57 + t125 * t56) * t323 + (t57 * (t193 * t329 + t157) + t56 * (-t193 * t330 - t159) + t391) * t324) * qJD(4)) * m(5) + (-(-t66 * pkin(2) * t475 + (-t167 * t625 - t330 * t398) * t330 + (-t166 * t625 - t329 * t398) * t329) * qJD(2) + t39 * t532 + t66 * t526 + (t167 * t465 + t39 * t207 - t66 * t221 + t374 * t97) * t330 + (t166 * t465 + t39 * t206 - t66 * t220 + t374 * t96) * t329 - g(3) * t625 - t374 * t615) * m(4) + (g(1) * t286 + g(2) * t285 - g(3) * t307 + t69 * t415 + (qJDD(2) * t306 + t307 * t337) * t688 + (-(-t285 * t329 - t286 * t330) * qJD(2) - t414) * (qJD(2) * t415 + qJD(1))) * m(3) + (t670 + (t680 * t327 + (t682 * t329 - t330 * t683) * t329) * t673 + (t677 * t327 + (t684 * t329 - t330 * t685) * t329) * t672) * t329 / 0.2e1 - (t671 + (t683 * t327 + (t680 - t682) * t698) * t673 + (t685 * t327 + (t677 - t684) * t698) * t672) * t330 / 0.2e1 + ((-t636 + t667) * t330 + (t635 + t668) * t329) * t479; t509 * (g(1) * t329 - g(2) * t330) + m(4) * (t329 * t97 - t330 * t96) + m(5) * t459 + m(6) * (t329 * t7 - t330 * t6); (t617 * t323 - t324 * t699) * t610 + (t616 * t323 - t324 * t700) * t609 + (t323 * t618 - t324 * t633) * t608 + (t271 * t620 - t330 * t619 - t621 * t646) * t607 + (t703 * t324 + (t329 * t664 + t330 * t663) * t323 + (t324 * t617 + t644) * qJD(2)) * t606 + (t268 * t621 + t269 * t620 - t329 * t619) * t605 + (t704 * t324 + (t329 * t666 + t330 * t665) * t323 + (t324 * t616 + t643) * qJD(2)) * t604 - (t174 * t658 + t175 * t659 + t633 * t276 + t635 * t283 + t636 * t284 + t634 * t512) * t324 / 0.2e1 + t671 * t554 / 0.2e1 + t670 * t553 / 0.2e1 + t662 * t520 / 0.2e1 + (t634 * t324 + (t329 * t636 + t330 * t635) * t323 + (t323 * t633 + t324 * t618) * qJD(2)) * t479 + (t679 * t324 + (-t333 * t621 + t620 * t335) * t323) * t478 + t668 * t493 / 0.2e1 + t667 * t324 * t481 + ((qJD(2) * t341 - t36 * t592 + t37 * t593 - t543 * t6 + t575 * t7) * t324 + g(1) * t535 - g(2) * t536 - (t27 * t535 - t37 * t467) * t284 - (t27 * t536 + t36 * t467) * t283 - (t36 * t535 + t37 * t536) * t512 + (t370 * qJD(2) + (-t36 * t576 - t534 * t6 + t373) * t330 + (t37 * t576 + t534 * t7 - t614) * t329 - g(3) * (-t333 * t669 - t588)) * t323) * m(6) + ((t35 * t123 - t34 * t125 - t56 * t91 + t57 * t89 + (t391 + (t329 * t57 - t330 * t56) * t193) * qJD(2)) * t324 + (t57 * (-qJD(2) * t123 + t109 * t329) + t56 * (qJD(2) * t125 - t109 * t330) + t16 * t423 + t38 * (-t329 * t91 + t330 * t89) + t459 * t193) * t323 - t57 * (t141 * t512 + t259 * t284) - t56 * (-t143 * t512 - t259 * t283) - t38 * (t141 * t283 - t143 * t284) - g(1) * t143 - g(2) * t141 - g(3) * t259) * m(5); (t364 * qJD(2) - (t283 * t329 - t284 * t330) * t585 + (t329 * t6 + t330 * t7 - t615) * t323 + (-t5 - t37 * (-t284 + t491) - t36 * (t283 - t490) + g(3)) * t324) * m(6);];
tau = t1;
