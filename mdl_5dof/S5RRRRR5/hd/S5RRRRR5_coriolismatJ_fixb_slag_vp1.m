% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRRR5_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR5_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:58:10
% EndTime: 2019-12-05 18:58:24
% DurationCPUTime: 8.57s
% Computational Cost: add. (76859->499), mult. (43211->597), div. (0->0), fcn. (38666->10), ass. (0->321)
t456 = qJ(1) + qJ(2);
t454 = qJ(3) + t456;
t446 = sin(t454);
t459 = cos(qJ(4));
t603 = pkin(4) * t459;
t448 = pkin(3) + t603;
t447 = cos(t454);
t455 = qJ(4) + qJ(5);
t449 = sin(t455);
t552 = t446 * t449;
t507 = -rSges(6,2) * t552 - t447 * rSges(6,3);
t461 = -pkin(9) - pkin(8);
t541 = t447 * t461;
t451 = cos(t455);
t599 = rSges(6,1) * t451;
t300 = t541 + (t448 + t599) * t446 + t507;
t450 = sin(t456);
t606 = pkin(2) * t450;
t290 = t300 + t606;
t602 = sin(qJ(1)) * pkin(1);
t284 = t290 + t602;
t544 = t447 * t451;
t545 = t447 * t449;
t474 = rSges(6,1) * t544 - rSges(6,2) * t545 + rSges(6,3) * t446;
t503 = t446 * t461 - t447 * t448;
t301 = -t474 + t503;
t452 = cos(t456);
t605 = pkin(2) * t452;
t291 = t301 - t605;
t607 = cos(qJ(1)) * pkin(1);
t285 = t291 - t607;
t146 = -t301 * t284 + t285 * t300;
t147 = -t301 * t290 + t291 * t300;
t441 = t447 * pkin(8);
t600 = rSges(5,1) * t459;
t495 = pkin(3) + t600;
t457 = sin(qJ(4));
t550 = t446 * t457;
t506 = -rSges(5,2) * t550 - t447 * rSges(5,3);
t313 = t446 * t495 - t441 + t506;
t303 = t313 + t606;
t298 = t303 + t602;
t543 = t447 * t457;
t426 = rSges(5,2) * t543;
t314 = t426 - t495 * t447 + (-rSges(5,3) - pkin(8)) * t446;
t304 = t314 - t605;
t299 = t304 - t607;
t155 = -t314 * t298 + t299 * t313;
t157 = -t314 * t303 + t304 * t313;
t393 = rSges(4,1) * t446 + rSges(4,2) * t447;
t373 = t393 + t606;
t363 = t373 + t602;
t394 = -rSges(4,1) * t447 + t446 * rSges(4,2);
t374 = t394 - t605;
t364 = t374 - t607;
t242 = -t394 * t363 + t364 * t393;
t262 = -t394 * t373 + t374 * t393;
t672 = m(6) / 0.2e1;
t673 = m(5) / 0.2e1;
t674 = m(4) / 0.2e1;
t499 = (-t262 + t242) * t674 + (-t147 + t146) * t672 + (-t157 + t155) * t673;
t500 = (t262 + t242) * t674 + (t147 + t146) * t672 + (t157 + t155) * t673;
t12 = t500 - t499;
t716 = t12 * qJD(1);
t338 = Icges(6,5) * t544 - Icges(6,6) * t545 + Icges(6,3) * t446;
t442 = Icges(6,4) * t451;
t694 = Icges(6,2) * t449 - t442;
t339 = Icges(6,6) * t447 + t694 * t446;
t403 = Icges(6,5) * t451 - Icges(6,6) * t449;
t560 = t403 * t446;
t522 = t447 * (Icges(6,3) * t447 - t560) + t339 * t552;
t340 = Icges(6,4) * t544 - Icges(6,2) * t545 + Icges(6,6) * t446;
t415 = Icges(6,4) * t545;
t342 = Icges(6,1) * t544 + Icges(6,5) * t446 - t415;
t698 = (t340 * t449 - t342 * t451) * t447;
t715 = t338 * t446 + t522 - t698;
t714 = t339 * t545;
t433 = rSges(5,1) * t457 + rSges(5,2) * t459;
t382 = t433 * t447;
t408 = rSges(6,1) * t449 + rSges(6,2) * t451;
t604 = pkin(4) * t457;
t487 = (t408 + t604) * t447;
t248 = t487 * t301;
t427 = pkin(4) * t550;
t558 = t408 * t446;
t350 = t427 + t558;
t527 = t350 * t300 - t248;
t567 = t487 * t291;
t529 = t350 * t290 - t567;
t381 = t433 * t446;
t573 = t313 * t381;
t576 = t303 * t381;
t594 = (-t576 + t573 + (t304 - t314) * t382) * t673 + (t527 - t529) * t672;
t551 = t446 * t451;
t371 = rSges(6,1) * t552 + rSges(6,2) * t551;
t333 = t427 + t371;
t156 = -t290 * t333 + t567;
t160 = -t300 * t333 + t248;
t180 = t304 * t382 - t576;
t184 = t314 * t382 - t573;
t702 = (t184 + t180) * t673 + (t160 + t156) * t672;
t713 = t594 - t702;
t226 = t487 * t285;
t531 = -t350 * t284 + t226;
t695 = -t298 * t381 + t299 * t382;
t595 = (t695 - t184) * t673 + (t527 + t531) * t672;
t153 = -t284 * t333 + t226;
t703 = (t184 + t695) * t673 + (t160 + t153) * t672;
t712 = t595 - t703;
t711 = -t446 / 0.2e1;
t656 = t446 / 0.2e1;
t654 = t447 / 0.2e1;
t653 = m(3) * (-t607 * (rSges(3,1) * t450 + rSges(3,2) * t452) - (-rSges(3,1) * t452 + t450 * rSges(3,2)) * t602);
t372 = t408 * t447;
t244 = t291 * t372;
t162 = -t290 * t371 + t244;
t710 = m(6) * t162;
t589 = Icges(6,4) * t449;
t404 = Icges(6,2) * t451 + t589;
t407 = Icges(6,1) * t451 - t589;
t697 = Icges(6,1) * t449 + t442;
t491 = (-t694 / 0.2e1 + t697 / 0.2e1) * t451 + (-t404 / 0.2e1 + t407 / 0.2e1) * t449;
t585 = t285 * t372;
t159 = -t284 * t371 + t585;
t257 = t301 * t372;
t171 = -t300 * t371 + t257;
t657 = m(6) * (t171 + t159);
t709 = t491 + t657 / 0.2e1;
t628 = m(6) * (t171 + t162);
t708 = t491 + t628 / 0.2e1;
t443 = t446 ^ 2;
t444 = t447 ^ 2;
t502 = t443 + t444;
t655 = -t447 / 0.2e1;
t705 = t654 + t655;
t704 = m(6) * t171;
t320 = t447 * t474;
t343 = -rSges(6,1) * t551 - t507;
t166 = -t447 * (pkin(3) * t447 + pkin(8) * t446 + t503) + t320 + (t541 - t343 + t441 + (-pkin(3) + t448) * t446) * t446;
t287 = -t371 * t446 - t447 * t372;
t410 = -rSges(6,2) * t449 + t599;
t556 = t410 * t447;
t557 = t410 * t446;
t480 = Icges(6,5) * t449 + Icges(6,6) * t451;
t365 = t446 * t480;
t366 = t480 * t447;
t414 = Icges(6,4) * t552;
t341 = -Icges(6,1) * t551 + Icges(6,5) * t447 + t414;
t515 = Icges(6,2) * t551 + t341 + t414;
t517 = -t446 * t697 + t339;
t686 = -t449 * t515 - t451 * t517;
t514 = -Icges(6,2) * t544 + t342 - t415;
t516 = t447 * t697 + t340;
t687 = t449 * t514 + t451 * t516;
t601 = (t444 * t365 + (t687 * t446 + (-t366 - t686) * t447) * t446) * t654 + (-t443 * t366 + (t686 * t447 + (t365 - t687) * t446) * t447) * t656;
t29 = t601 + m(6) * (t166 * t287 + t350 * t557 + t487 * t556);
t701 = t29 * qJD(5);
t435 = -rSges(5,2) * t457 + t600;
t700 = t435 * t673;
t542 = t447 * t459;
t357 = Icges(5,4) * t542 - Icges(5,2) * t543 + Icges(5,6) * t446;
t424 = Icges(5,4) * t543;
t359 = Icges(5,1) * t542 + Icges(5,5) * t446 - t424;
t699 = (t357 * t457 - t359 * t459) * t447;
t139 = -t291 * t284 + t285 * t290;
t152 = -t304 * t298 + t299 * t303;
t220 = -t374 * t363 + t364 * t373;
t453 = Icges(5,4) * t459;
t696 = Icges(5,1) * t457 + t453;
t693 = Icges(5,2) * t457 - t453;
t692 = -m(5) * t184 - m(6) * t160;
t691 = qJD(1) + qJD(2) + qJD(3);
t596 = (t695 - t180) * t673 + (t529 + t531) * t672;
t689 = (t180 + t695) * t673 + (t156 + t153) * t672;
t688 = (-t694 + t697) * t449 + (t404 - t407) * t451;
t590 = Icges(5,4) * t457;
t429 = Icges(5,2) * t459 + t590;
t432 = Icges(5,1) * t459 - t590;
t685 = (-t693 + t696) * t457 + (t429 - t432) * t459;
t510 = -Icges(5,2) * t542 + t359 - t424;
t512 = t447 * t696 + t357;
t684 = t457 * t510 + t459 * t512;
t423 = Icges(5,4) * t550;
t549 = t446 * t459;
t358 = -Icges(5,1) * t549 + Icges(5,5) * t447 + t423;
t511 = Icges(5,2) * t549 + t358 + t423;
t356 = Icges(5,6) * t447 + t693 * t446;
t513 = -t446 * t696 + t356;
t683 = -t457 * t511 - t459 * t513;
t682 = (t696 / 0.2e1 - t693 / 0.2e1) * t459 + (t432 / 0.2e1 - t429 / 0.2e1) * t457;
t185 = -t341 * t551 + t522;
t186 = t447 * t338 + t340 * t552 - t342 * t551;
t493 = -t341 * t451 - t338;
t569 = t339 * t449;
t19 = (t185 * t447 + t186 * t446) * t711 + ((t698 + t715) * t447 + ((t493 - t569) * t447 + t186 + t714) * t446) * t656 + ((t186 + (-t338 + t569) * t447 - t714) * t447 + (t446 * t493 - t185 + t715) * t446) * t654;
t680 = 0.4e1 * qJD(1);
t678 = 0.4e1 * qJD(2);
t677 = 0.2e1 * qJD(3);
t676 = 2 * qJD(4);
t659 = m(6) * (t162 + t159);
t649 = m(4) * t220;
t647 = m(4) * t242;
t646 = m(4) * t262;
t638 = m(5) * t152;
t636 = m(5) * t155;
t635 = m(5) * t157;
t633 = m(5) * t695;
t632 = m(5) * t180;
t497 = t290 * t558;
t526 = -t284 * t558 + t585;
t627 = m(6) * (t497 + t526 - t244);
t267 = t300 * t558;
t626 = m(6) * (t267 + t526 - t257);
t625 = m(6) * (-t497 + t267 + (t291 - t301) * t372);
t250 = t285 * t557;
t523 = t350 * t372 - t371 * t487;
t624 = m(6) * (t284 * t556 + t250 + t523);
t297 = t487 * t558;
t570 = t333 * t408;
t586 = t284 * t410;
t623 = m(6) * (t250 + t297 + (-t570 + t586) * t447);
t255 = t291 * t557;
t622 = m(6) * (t290 * t556 + t255 + t523);
t584 = t290 * t410;
t621 = m(6) * (t255 + t297 + (-t570 + t584) * t447);
t268 = t301 * t557;
t620 = m(6) * (t300 * t556 + t268 + t523);
t579 = t300 * t410;
t619 = m(6) * (t268 + t297 + (-t570 + t579) * t447);
t618 = m(6) * t139;
t616 = m(6) * t146;
t615 = m(6) * t147;
t613 = m(6) * t153;
t612 = m(6) * t156;
t611 = m(6) * t159;
t566 = t356 * t457;
t565 = t371 * t408;
t428 = Icges(5,5) * t459 - Icges(5,6) * t457;
t554 = t428 * t446;
t524 = (-t333 + t350) * t487;
t354 = Icges(5,3) * t447 - t554;
t519 = t447 * t354 + t356 * t550;
t518 = t446 * t354 + t358 * t542;
t494 = t410 + t603;
t355 = Icges(5,5) * t542 - Icges(5,6) * t543 + Icges(5,3) * t446;
t492 = -t358 * t459 - t355;
t486 = t659 / 0.2e1 + t491;
t481 = Icges(5,5) * t457 + Icges(5,6) * t459;
t207 = t447 * t355 + t357 * t550 - t359 * t549;
t375 = t446 * t481;
t473 = -t19 + (-t688 * t447 - t449 * t516 + t451 * t514 + t560) * t656 + (t403 * t447 + t688 * t446 - t449 * t517 + t451 * t515) * t654;
t472 = -t491 + t705 * (t451 * t340 + t449 * t342);
t471 = t491 + t682;
t467 = t471 + t689;
t464 = t472 - t682 + t705 * (t459 * t357 + t359 * t457);
t206 = -t358 * t549 + t519;
t144 = t206 * t447 + t207 * t446;
t208 = -t356 * t543 + t518;
t209 = t355 * t446 - t699;
t145 = t208 * t447 + t209 * t446;
t46 = (t209 + t519 + t699) * t447 + (-t208 + (t492 - t566) * t447 + t207 + t518) * t446;
t47 = (t207 + (-t355 + t566) * t447 - t518) * t447 + (t446 * t492 - t206 + t519) * t446;
t462 = (t46 * t711 + t473 + (t47 + t145) * t655 + (t428 * t447 + t685 * t446 - t457 * t513 + t459 * t511) * t654 + (-t685 * t447 - t457 * t512 + t459 * t510 + t144 + t554) * t656) * qJD(4);
t376 = t481 * t447;
t353 = t494 * t447;
t351 = t494 * t446;
t307 = t372 * t558;
t259 = -t446 * t343 + t320;
t252 = -t502 * t604 + t287;
t142 = t491 + t704;
t136 = t491 + t710;
t132 = t619 / 0.2e1;
t130 = t491 + t611;
t129 = t620 / 0.2e1;
t122 = t621 / 0.2e1;
t120 = t622 / 0.2e1;
t111 = t623 / 0.2e1;
t110 = t624 / 0.2e1;
t105 = t625 / 0.2e1;
t103 = t626 / 0.2e1;
t99 = t627 / 0.2e1;
t68 = t471 - t692;
t57 = t471 + t612 + t632;
t56 = t471 + t613 + t633;
t55 = t615 + t635 + t646;
t54 = t616 + t636 + t647;
t48 = t618 + t638 + t649 + t653;
t37 = -t625 / 0.2e1 + t708;
t36 = t105 + t708;
t35 = -t626 / 0.2e1 + t709;
t34 = t103 + t709;
t33 = -t627 / 0.2e1 + t486;
t32 = t99 + t486;
t31 = m(6) * (t502 * t408 * t410 + t259 * t287) + t601;
t30 = t31 * qJD(5);
t28 = t105 - t628 / 0.2e1 + t472;
t27 = t103 - t657 / 0.2e1 + t472;
t26 = t99 - t659 / 0.2e1 + t472;
t25 = t471 - t713;
t24 = t471 + t594 + t702;
t23 = t471 + t595 + t703;
t22 = t471 - t712;
t21 = t467 + t596;
t20 = t467 - t596;
t18 = t19 * qJD(5);
t17 = t464 + t713;
t16 = t464 + t712;
t15 = t464 + t596 - t689;
t13 = t499 + t500;
t11 = t129 - t619 / 0.2e1 + t19;
t10 = t132 - t620 / 0.2e1 + t19;
t9 = t120 - t621 / 0.2e1 + t19;
t8 = t122 - t622 / 0.2e1 + t19;
t7 = t110 - t623 / 0.2e1 + t19;
t6 = t111 - t624 / 0.2e1 + t19;
t5 = t129 + t132 + t473;
t4 = t120 + t122 + t473;
t3 = t110 + t111 + t473;
t2 = (t145 / 0.2e1 + t47 / 0.2e1) * t447 + (t46 / 0.2e1 - t144 / 0.2e1) * t446 + t19;
t1 = t2 * qJD(4);
t14 = [qJD(2) * t48 + qJD(3) * t54 + qJD(4) * t56 + qJD(5) * t130, t48 * qJD(1) + t13 * qJD(3) + t21 * qJD(4) + t32 * qJD(5) + 0.2e1 * (t653 / 0.2e1 + t139 * t672 + t152 * t673 + t220 * t674) * qJD(2), t54 * qJD(1) + t13 * qJD(2) + t23 * qJD(4) + t34 * qJD(5) + (t146 * t672 + t155 * t673 + t242 * t674) * t677, t56 * qJD(1) + t21 * qJD(2) + t23 * qJD(3) + t462 + t3 * qJD(5) + ((t298 * t447 + t299 * t446) * t700 + (t284 * t353 + t285 * t351 + t524) * t672) * t676, t130 * qJD(1) + t32 * qJD(2) + t34 * qJD(3) + t3 * qJD(4) + ((t250 + t307 + (-t565 + t586) * t447) * m(6) + t473) * qJD(5); t12 * qJD(3) + t20 * qJD(4) + t33 * qJD(5) + (-t653 / 0.4e1 - t649 / 0.4e1 - t638 / 0.4e1 - t618 / 0.4e1) * t680, qJD(3) * t55 + qJD(4) * t57 + qJD(5) * t136, t716 + t55 * qJD(2) + t24 * qJD(4) + t36 * qJD(5) + (t147 * t672 + t157 * t673 + t262 * t674) * t677, t20 * qJD(1) + t57 * qJD(2) + t24 * qJD(3) + t462 + t4 * qJD(5) + ((t290 * t353 + t291 * t351 + t524) * t672 + (t303 * t447 + t304 * t446) * t700) * t676, t33 * qJD(1) + t136 * qJD(2) + t36 * qJD(3) + t4 * qJD(4) + ((t255 + t307 + (-t565 + t584) * t447) * m(6) + t473) * qJD(5); -t12 * qJD(2) + t22 * qJD(4) + t35 * qJD(5) + (-t647 / 0.4e1 - t636 / 0.4e1 - t616 / 0.4e1) * t680, -t716 + t25 * qJD(4) + t37 * qJD(5) + (-t615 / 0.4e1 - t635 / 0.4e1 - t646 / 0.4e1) * t678, qJD(4) * t68 + qJD(5) * t142, t22 * qJD(1) + t25 * qJD(2) + t68 * qJD(3) + t462 + t5 * qJD(5) + ((t313 * t447 + t314 * t446) * t700 + (t300 * t353 + t301 * t351 + t524) * t672) * t676, t35 * qJD(1) + t37 * qJD(2) + t142 * qJD(3) + t5 * qJD(4) + ((t268 + t307 + (-t565 + t579) * t447) * m(6) + t473) * qJD(5); t15 * qJD(2) + t16 * qJD(3) + t1 + t7 * qJD(5) + (-t613 / 0.4e1 - t633 / 0.4e1) * t680 + t464 * qJD(1), t15 * qJD(1) + t464 * qJD(2) + t17 * qJD(3) + t1 + t9 * qJD(5) + (-t612 / 0.4e1 - t632 / 0.4e1) * t678, t16 * qJD(1) + t17 * qJD(2) + t11 * qJD(5) + t1 + (t464 + t692) * qJD(3), (m(5) * ((t447 * (rSges(5,1) * t542 + rSges(5,3) * t446 - t426) - t446 * (-rSges(5,1) * t549 - t506)) * (-t381 * t446 - t382 * t447) + t502 * t435 * t433) + (t444 * t375 + (t684 * t446 + (-t376 - t683) * t447) * t446) * t654 + (-t443 * t376 + (t683 * t447 + (t375 - t684) * t446) * t447) * t656 + m(6) * (t166 * t252 + t350 * t351 + t353 * t487) + t601) * qJD(4) + t701 + t691 * t2, t7 * qJD(1) + t9 * qJD(2) + t11 * qJD(3) + t29 * qJD(4) + t701; (t472 - t611) * qJD(1) + t26 * qJD(2) + t27 * qJD(3) + t6 * qJD(4) + t18, t26 * qJD(1) + (t472 - t710) * qJD(2) + t28 * qJD(3) + t8 * qJD(4) + t18, t27 * qJD(1) + t28 * qJD(2) + (t472 - t704) * qJD(3) + t10 * qJD(4) + t18, t6 * qJD(1) + t8 * qJD(2) + t10 * qJD(3) + ((t252 * t259 + (t351 * t446 + t353 * t447) * t408) * m(6) + t601) * qJD(4) + t30, qJD(4) * t31 + t691 * t19 + t30;];
Cq = t14;
