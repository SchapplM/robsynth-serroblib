% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% Cq [5x5]
%   matrix of coriolis and centrifugal joint torques.
%   Gives coriolis joint torques when multiplied with joint velocities

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRRPP3_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP3_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_coriolismatJ_fixb_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP3_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPP3_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPP3_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:53:09
% EndTime: 2019-12-31 20:53:24
% DurationCPUTime: 10.79s
% Computational Cost: add. (27824->525), mult. (28705->635), div. (0->0), fcn. (25883->6), ass. (0->316)
t473 = cos(qJ(3));
t467 = Icges(4,4) * t473;
t471 = sin(qJ(3));
t607 = Icges(5,6) * t473;
t733 = t467 + t607 + (Icges(4,1) + Icges(5,2)) * t471;
t470 = qJ(1) + qJ(2);
t464 = sin(t470);
t465 = cos(t470);
t619 = rSges(6,3) + qJ(5);
t527 = pkin(3) + t619;
t621 = rSges(6,2) + qJ(4);
t659 = rSges(6,1) + pkin(4);
t226 = (pkin(7) + t659) * t464 + (t471 * t621 + t473 * t527 + pkin(2)) * t465;
t631 = cos(qJ(1)) * pkin(1);
t224 = t226 + t631;
t583 = t465 * t471;
t211 = t224 * t583;
t582 = t465 * t473;
t451 = rSges(5,2) * t582;
t620 = rSges(5,3) + qJ(4);
t630 = pkin(3) * t473;
t260 = -t451 + (rSges(5,1) + pkin(7)) * t464 + (t471 * t620 + pkin(2) + t630) * t465;
t249 = t260 + t631;
t231 = t249 * t583;
t460 = t465 * pkin(7);
t589 = t464 * t471;
t518 = rSges(5,1) * t465 - rSges(5,3) * t589;
t604 = qJ(4) * t471;
t259 = t460 + (-t604 - pkin(2) + (rSges(5,2) - pkin(3)) * t473) * t464 + t518;
t632 = sin(qJ(1)) * pkin(1);
t248 = t259 - t632;
t422 = t604 + t630;
t588 = t464 * t473;
t427 = qJ(5) * t588;
t520 = -rSges(6,2) * t589 - rSges(6,3) * t588 - t427;
t225 = t460 + t659 * t465 + (-pkin(2) - t422) * t464 + t520;
t223 = t225 - t632;
t517 = t464 * (-t223 - t225);
t521 = t260 * t583;
t523 = t226 * t583;
t689 = m(6) / 0.2e1;
t690 = m(5) / 0.2e1;
t628 = (t471 * t517 + t211 + t523) * t689 + (t231 + t521 + (-t248 - t259) * t589) * t690;
t723 = t223 - t225;
t516 = t464 * t723;
t722 = t248 - t259;
t629 = (-t471 * t516 + t211 - t523) * t689 + (-t589 * t722 + t231 - t521) * t690;
t6 = t629 - t628;
t732 = t6 * qJD(1);
t325 = Icges(4,4) * t588 - Icges(4,2) * t589 - Icges(4,6) * t465;
t332 = Icges(5,5) * t465 + Icges(5,6) * t588 - Icges(5,3) * t589;
t731 = t325 + t332;
t442 = Icges(4,4) * t589;
t327 = Icges(4,1) * t588 - Icges(4,5) * t465 - t442;
t434 = Icges(5,6) * t589;
t336 = Icges(5,4) * t465 + Icges(5,2) * t588 - t434;
t730 = t327 + t336;
t405 = Icges(5,3) * t471 - t607;
t606 = Icges(6,6) * t473;
t729 = Icges(6,2) * t471 + t405 + t606;
t611 = Icges(4,4) * t471;
t413 = Icges(4,2) * t473 + t611;
t608 = Icges(5,6) * t471;
t728 = -Icges(5,2) * t473 + t413 + t608;
t414 = -Icges(4,2) * t471 + t467;
t726 = Icges(6,3) * t471 + t414 - t606 + t733;
t466 = Icges(6,6) * t471;
t402 = Icges(6,3) * t473 + t466;
t407 = -Icges(6,2) * t473 + t466;
t416 = Icges(4,1) * t473 - t611;
t725 = -Icges(5,3) * t473 + t402 + t407 + t416 - t608;
t417 = pkin(3) * t471 - qJ(4) * t473;
t623 = rSges(5,2) * t471;
t506 = rSges(5,3) * t473 + t623;
t535 = t417 - t506;
t316 = t535 * t464;
t318 = t535 * t465;
t509 = -t316 * t583 + t318 * t589;
t505 = rSges(6,2) * t473 - rSges(6,3) * t471;
t603 = qJ(5) * t471;
t508 = t417 - t505 + t603;
t282 = t508 * t464;
t284 = t508 * t465;
t510 = -t282 * t583 + t284 * t589;
t560 = t259 * t582 + t260 * t588;
t563 = t225 * t582 + t226 * t588;
t615 = (t509 + t560) * t690 + (t510 + t563) * t689;
t452 = pkin(3) * t589;
t292 = t452 + (-t473 * t620 - t623) * t464;
t428 = qJ(4) * t582;
t511 = -pkin(3) * t583 + t428;
t531 = rSges(5,2) * t583 + rSges(5,3) * t582;
t293 = t511 + t531;
t556 = t292 * t583 + t293 * t589;
t272 = t452 + (t471 * t619 - t473 * t621) * t464;
t450 = rSges(6,2) * t582;
t273 = -t527 * t583 + t428 + t450;
t559 = t272 * t583 + t273 * t589;
t617 = (t556 + t560) * t690 + (t559 + t563) * t689;
t12 = t617 - t615;
t561 = t248 * t582 + t249 * t588;
t564 = t223 * t582 + t224 * t588;
t616 = (t509 + t561) * t690 + (t510 + t564) * t689;
t618 = (t556 + t561) * t690 + (t559 + t564) * t689;
t9 = t618 - t616;
t724 = t9 * qJD(1) + t12 * qJD(2);
t462 = t464 ^ 2;
t463 = t465 ^ 2;
t529 = t462 + t463;
t719 = (Icges(6,4) + Icges(5,5) - Icges(4,6)) * t473 + (Icges(5,4) - Icges(4,5) - Icges(6,5)) * t471;
t431 = Icges(6,6) * t589;
t330 = Icges(6,5) * t465 - Icges(6,3) * t588 - t431;
t718 = -t330 + t431 - t434 - t442 + t730 + (-Icges(4,2) - Icges(6,2) - Icges(5,3)) * t588;
t328 = Icges(4,5) * t464 + t416 * t465;
t329 = Icges(6,5) * t464 + t402 * t465;
t435 = Icges(5,6) * t583;
t335 = Icges(5,4) * t464 - Icges(5,2) * t582 + t435;
t717 = Icges(5,3) * t582 - t328 - t329 + t335 + t435 + (-t407 + t413) * t465;
t432 = Icges(6,6) * t588;
t334 = Icges(6,4) * t465 - Icges(6,2) * t589 - t432;
t716 = Icges(6,3) * t589 + t464 * t733 + t334 - t432 + t731;
t326 = Icges(4,6) * t464 + t414 * t465;
t331 = Icges(5,5) * t464 + t405 * t465;
t433 = Icges(6,6) * t582;
t333 = Icges(6,4) * t464 + Icges(6,2) * t583 + t433;
t715 = -Icges(6,3) * t583 - t465 * t733 - t326 + t331 + t333 + t433;
t691 = m(4) / 0.2e1;
t658 = m(3) * (t631 * (-rSges(3,1) * t464 - rSges(3,2) * t465) + (rSges(3,1) * t465 - rSges(3,2) * t464) * t632);
t522 = t226 * t582;
t140 = -t225 * t588 + t522;
t710 = t140 * m(6) * qJD(2);
t338 = Icges(6,1) * t465 - Icges(6,4) * t589 - Icges(6,5) * t588;
t313 = t464 * t338;
t185 = -t330 * t582 - t334 * t583 - t313;
t709 = t185 * t465;
t374 = t529 * t473;
t320 = (t374 - t473) * t471;
t625 = m(6) * qJD(5);
t708 = t320 * t625;
t707 = (t725 - t728) * t473 + (-t726 + t729) * t471;
t704 = (t330 * t473 + t334 * t471) * t464;
t703 = (m(5) / 0.4e1 + m(6) / 0.4e1) * t320;
t68 = -t223 * t226 + t224 * t225;
t102 = -t248 * t260 + t249 * t259;
t624 = rSges(4,1) * t473;
t519 = pkin(2) + t624;
t532 = rSges(4,2) * t589 + rSges(4,3) * t465;
t290 = -t464 * t519 + t460 + t532;
t278 = t290 - t632;
t448 = rSges(4,2) * t583;
t291 = -t448 + t519 * t465 + (rSges(4,3) + pkin(7)) * t464;
t279 = t291 + t631;
t136 = -t278 * t291 + t279 * t290;
t702 = t719 * t465;
t701 = t719 * t464;
t700 = t471 * t718 + t473 * t716;
t699 = t471 * t717 + t473 * t715;
t420 = rSges(4,1) * t471 + rSges(4,2) * t473;
t526 = ((-t279 + t291) * t465 + (t278 - t290) * t464) * t420 * t691 + ((-t224 + t226) * t284 + t723 * t282) * t689 + ((-t249 + t260) * t318 + t722 * t316) * t690;
t101 = t225 * t272 + t226 * t273;
t129 = t248 * t292 + t249 * t293;
t131 = t259 * t292 + t260 * t293;
t368 = t420 * t464;
t370 = t420 * t465;
t155 = t278 * t368 - t279 * t370;
t165 = t290 * t368 - t291 * t370;
t98 = t223 * t272 + t224 * t273;
t698 = (t165 + t155) * t691 + (t101 + t98) * t689 + (t131 + t129) * t690;
t478 = (-t729 / 0.2e1 + t726 / 0.2e1) * t473 + (-t728 / 0.2e1 + t725 / 0.2e1) * t471;
t696 = 0.4e1 * qJD(1);
t694 = 0.4e1 * qJD(2);
t693 = 2 * qJD(3);
t212 = t224 * t582;
t678 = m(6) * (-t473 * t516 + t212 - t522);
t677 = m(6) * (t473 * t517 + t212 + t522);
t676 = m(6) * t68;
t558 = t272 * t582 + t273 * t588;
t674 = m(6) * ((-t223 * t465 - t224 * t464) * t471 + t558);
t672 = m(6) * ((-t225 * t465 - t226 * t464) * t471 + t558);
t271 = t282 * t582;
t600 = t284 * t473;
t670 = m(6) * (-t223 * t583 - t271 + (-t224 * t471 + t600) * t464);
t668 = m(6) * (-t225 * t583 - t271 + (-t226 * t471 + t600) * t464);
t664 = m(6) * t98;
t663 = -t464 / 0.2e1;
t662 = t464 / 0.2e1;
t661 = -t465 / 0.2e1;
t654 = m(4) * t136;
t652 = m(4) * t155;
t651 = m(4) * t165;
t650 = m(5) * t102;
t644 = m(5) * t129;
t642 = m(5) * t131;
t641 = m(5) * (-t248 * t589 + t231);
t640 = m(5) * (-t259 * t589 + t521);
t639 = m(6) * t101;
t638 = m(6) * (-t223 * t589 + t211);
t636 = m(6) * (-t225 * t589 + t523);
t468 = t471 ^ 2;
t469 = t473 ^ 2;
t530 = t529 * t469;
t635 = m(6) * (-t469 + (0.1e1 - t529) * t468 + t530);
t634 = m(6) * (-t374 * t473 - t468 * t529);
t373 = t529 * t471;
t633 = m(6) * (t373 * t471 + t530);
t626 = m(6) * qJD(3);
t622 = rSges(6,2) * t471;
t598 = t325 * t471;
t597 = t336 * t473;
t410 = Icges(4,5) * t473 - Icges(4,6) * t471;
t595 = t410 * t465;
t411 = Icges(6,4) * t471 + Icges(6,5) * t473;
t594 = t411 * t465;
t412 = -Icges(5,4) * t473 + Icges(5,5) * t471;
t593 = t412 * t465;
t586 = t465 * t338;
t565 = -t272 * t284 - t273 * t282;
t562 = -t292 * t318 - t293 * t316;
t557 = -t282 * t588 - t284 * t582;
t555 = -t316 * t588 - t318 * t582;
t339 = Icges(5,1) * t464 + t593;
t554 = t331 * t583 + t339 * t464;
t340 = Icges(5,1) * t465 + Icges(5,4) * t588 - Icges(5,5) * t589;
t553 = t332 * t583 + t464 * t340;
t323 = Icges(4,5) * t588 - Icges(4,6) * t589 - Icges(4,3) * t465;
t552 = -t464 * t323 - t327 * t582;
t324 = Icges(4,3) * t464 + t595;
t551 = t324 * t464 + t328 * t582;
t538 = t464 * (qJ(4) * t588 - t452) + t465 * t511;
t537 = t529 * t422;
t534 = -rSges(6,3) * t473 - t422 - t622;
t533 = rSges(5,2) * t473 - rSges(5,3) * t471 - t422;
t133 = -t223 * t588 + t212;
t524 = m(6) * t133 * qJD(1);
t337 = Icges(6,1) * t464 + t594;
t184 = t329 * t582 + t333 * t583 + t337 * t464;
t296 = t331 * t589;
t515 = t465 * t339 - t296;
t298 = t328 * t588;
t514 = t465 * t324 - t298;
t513 = t326 * t471 - t323;
t512 = t335 * t473 + t340;
t141 = -t520 * t464 + t537 + (t473 * t619 + t622) * t463;
t177 = -t464 * (rSges(5,2) * t588 + t518) + t465 * (t464 * rSges(5,1) + rSges(5,3) * t583 - t451) + t537;
t43 = m(5) * (t177 * t373 + t555) + m(6) * (t141 * t373 + t557);
t507 = -t329 * t588 - t333 * t589 + t337 * t465;
t499 = t282 * t464 + t284 * t465;
t491 = t184 - t586;
t283 = t464 * t534 - t427;
t285 = (-qJ(5) * t473 + t534) * t465;
t489 = t283 * t464 + t285 * t465 + t141;
t488 = (-t368 * t465 + t370 * t464) * t420;
t477 = t478 + t698;
t476 = -t478 + ((-t334 + t731) * t473 + (t330 + t730) * t471) * (t662 + t663);
t181 = -t326 * t589 - t514;
t119 = -(-t464 * (-t327 * t473 + t598) - t465 * t323) * t465 + t181 * t464;
t182 = -t325 * t583 - t552;
t183 = -t326 * t583 + t551;
t120 = -t182 * t465 + t183 * t464;
t121 = t184 * t464 - t709;
t186 = -t335 * t582 + t554;
t187 = t336 * t582 - t553;
t122 = t186 * t464 - t187 * t465;
t189 = t586 - t704;
t123 = -t189 * t465 - t464 * t507;
t190 = -t335 * t588 - t515;
t124 = t190 * t464 - (-t464 * (t332 * t471 - t597) + t465 * t340) * t465;
t35 = (t465 * t513 + t183 - t551) * t465 + (t464 * t513 + t182 + t514) * t464;
t36 = (t181 - t298 + (t324 + t598) * t465 + t552) * t465 + t551 * t464;
t37 = -t709 + (t189 + t491 + t704) * t464;
t38 = (t190 - t296 + (t339 - t597) * t465 + t553) * t465 + t554 * t464;
t39 = (t184 - t491) * t465 + (t185 + t507 + t313) * t464;
t40 = (t465 * t512 + t186 - t554) * t465 + (t464 * t512 + t187 + t515) * t464;
t475 = ((t38 + t37 + t36) * t465 / 0.2e1 + (t124 + t123 + t119 + t40 + t39 + t35) * t663 + (-t717 * t473 + t715 * t471 + t707 * t465 + (t410 + t411 + t412) * t464) * t662 + (t707 * t464 - t471 * t716 + t473 * t718 + t120 + t121 + t122 - t593 - t594 - t595) * t661) * qJD(3);
t425 = -rSges(4,2) * t471 + t624;
t319 = t533 * t465;
t317 = t533 * t464;
t281 = t633 / 0.2e1;
t280 = t634 / 0.2e1;
t276 = t635 / 0.2e1;
t238 = 0.4e1 * t703;
t194 = t462 * t506 + t465 * t531 + t538;
t164 = t465 * (-rSges(6,3) * t583 + t450) + t505 * t462 - t529 * t603 + t538;
t154 = t281 + t276 - t634 / 0.2e1;
t153 = t280 + t281 - t635 / 0.2e1;
t152 = t280 + t276 - t633 / 0.2e1;
t88 = t668 / 0.2e1;
t82 = t141 * t374 + t471 * t499;
t78 = t670 / 0.2e1;
t74 = t672 / 0.2e1;
t70 = t674 / 0.2e1;
t67 = t636 + t640;
t66 = t638 + t641;
t61 = t677 / 0.2e1;
t60 = t678 / 0.2e1;
t28 = t478 + t639 + t642 + t651;
t25 = t478 + t644 + t652 + t664;
t24 = t650 + t654 + t658 + t676;
t23 = t88 - t672 / 0.2e1;
t22 = t88 + t74;
t21 = t74 - t668 / 0.2e1;
t20 = t78 - t674 / 0.2e1;
t19 = t78 + t70;
t18 = t70 - t670 / 0.2e1;
t17 = t61 - t678 / 0.2e1;
t16 = t61 + t60;
t15 = t60 - t677 / 0.2e1;
t14 = t615 + t617;
t10 = t616 + t618;
t7 = t628 + t629;
t5 = t477 + t526;
t4 = t477 - t526;
t3 = t476 + t526 - t698;
t2 = (-t38 / 0.2e1 + t120 / 0.2e1 - t37 / 0.2e1 - t36 / 0.2e1 + t121 / 0.2e1 + t122 / 0.2e1) * t465 + (t124 / 0.2e1 + t35 / 0.2e1 + t123 / 0.2e1 + t119 / 0.2e1 + t39 / 0.2e1 + t40 / 0.2e1) * t464;
t1 = t2 * qJD(3);
t8 = [qJD(2) * t24 + qJD(3) * t25 + qJD(4) * t66 + t133 * t625, t24 * qJD(1) + t5 * qJD(3) + t7 * qJD(4) + t16 * qJD(5) + 0.2e1 * (t658 / 0.2e1 + t102 * t690 + t136 * t691 + t68 * t689) * qJD(2), t25 * qJD(1) + t5 * qJD(2) + t475 + t10 * qJD(4) + t19 * qJD(5) + ((t248 * t319 + t249 * t317 + t562) * t690 + (t223 * t285 + t224 * t283 + t565) * t689 + ((-t278 * t465 - t279 * t464) * t425 + t488) * t691) * t693, qJD(1) * t66 + qJD(2) * t7 + qJD(3) * t10, qJD(2) * t16 + qJD(3) * t19 + t524; t4 * qJD(3) - t6 * qJD(4) + t17 * qJD(5) + (-t676 / 0.4e1 - t650 / 0.4e1 - t654 / 0.4e1 - t658 / 0.4e1) * t696, qJD(3) * t28 + qJD(4) * t67 + t140 * t625, t4 * qJD(1) + t28 * qJD(2) + t475 + t14 * qJD(4) + t22 * qJD(5) + ((t259 * t319 + t260 * t317 + t562) * t690 + (t225 * t285 + t226 * t283 + t565) * t689 + ((-t290 * t465 - t291 * t464) * t425 + t488) * t691) * t693, qJD(2) * t67 + qJD(3) * t14 - t732, qJD(1) * t17 + qJD(3) * t22 + t710; t476 * qJD(1) + t3 * qJD(2) + t1 - t9 * qJD(4) + t20 * qJD(5) + (-t652 / 0.4e1 - t644 / 0.4e1 - t664 / 0.4e1) * t696, t3 * qJD(1) + t476 * qJD(2) + t1 - t12 * qJD(4) + t23 * qJD(5) + (-t651 / 0.4e1 - t642 / 0.4e1 - t639 / 0.4e1) * t694, t43 * qJD(4) + t82 * t625 + (qJD(1) + qJD(2)) * t2 + (m(6) * (t141 * t164 - t282 * t283 - t284 * t285) + m(5) * (t177 * t194 - t316 * t317 - t318 * t319) + m(4) * ((t464 * (rSges(4,1) * t588 - t532) + t465 * (rSges(4,1) * t582 + t464 * rSges(4,3) - t448)) * (-t368 * t464 - t370 * t465) + t529 * t425 * t420) + (t702 * t462 + (t700 * t465 + (t699 - t701) * t464) * t465) * t662 + (t701 * t463 + (t699 * t464 + (t700 - t702) * t465) * t464) * t661) * qJD(3), qJD(3) * t43 - 0.4e1 * qJD(4) * t703 + qJD(5) * t153 - t724, t20 * qJD(1) + t23 * qJD(2) + t153 * qJD(4) + t82 * t626 + t708; t6 * qJD(2) + t9 * qJD(3) + (-t638 / 0.4e1 - t641 / 0.4e1) * t696, t732 + t12 * qJD(3) + (-t640 / 0.4e1 - t636 / 0.4e1) * t694, (m(6) * (-t164 * t473 + t557) + m(5) * (-t194 * t473 + t555) + 0.2e1 * (t489 * t689 + (t317 * t464 + t319 * t465 + t177) * t690) * t471 - t43) * qJD(3) + t238 * qJD(4) + t152 * qJD(5) + t724, t238 * qJD(3), t152 * qJD(3); t15 * qJD(2) + t18 * qJD(3) - t524, qJD(1) * t15 + qJD(3) * t21 - t710, t18 * qJD(1) + t21 * qJD(2) + (t489 * t473 + (t164 + t499) * t471 - t82) * t626 + t154 * qJD(4) - t708, t154 * qJD(3), -t320 * t626;];
Cq = t8;
