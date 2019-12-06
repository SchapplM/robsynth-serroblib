% Calculate matrix of centrifugal and coriolis load on the joints for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S5RRPRR6_coriolismatJ_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR6_coriolismatJ_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:35:43
% EndTime: 2019-12-05 18:36:04
% DurationCPUTime: 15.04s
% Computational Cost: add. (80557->527), mult. (74468->709), div. (0->0), fcn. (80559->10), ass. (0->334)
t477 = qJ(1) + qJ(2);
t475 = cos(t477);
t470 = t475 * qJ(3);
t474 = sin(t477);
t479 = cos(pkin(9));
t482 = cos(qJ(4));
t519 = t479 * (pkin(4) * t482 + pkin(3));
t478 = sin(pkin(9));
t522 = -t478 * rSges(6,3) - pkin(2);
t480 = sin(qJ(4));
t601 = t475 * t480;
t467 = pkin(4) * t601;
t540 = t478 * (-pkin(8) - pkin(7));
t545 = t474 * t540 + t467;
t591 = qJ(4) + qJ(5);
t521 = cos(t591);
t516 = t479 * t521;
t473 = sin(t591);
t603 = t475 * t473;
t423 = t474 * t516 - t603;
t609 = t474 * t473;
t496 = t475 * t521 + t479 * t609;
t553 = t423 * rSges(6,1) - t496 * rSges(6,2);
t278 = -t470 + (t519 - t522) * t474 - t545 + t553;
t364 = t496 * rSges(6,1) + t423 * rSges(6,2);
t424 = -t474 * t521 + t479 * t603;
t425 = t475 * t516 + t609;
t365 = -t424 * rSges(6,1) - t425 * rSges(6,2);
t607 = t474 * t480;
t466 = pkin(4) * t607;
t535 = -t466 + (-t519 + t540) * t475;
t610 = t474 * qJ(3);
t703 = -t425 * rSges(6,1) + t424 * rSges(6,2);
t716 = t475 * t522 + t535 - t610 + t703;
t141 = -t278 * t364 - t365 * t716;
t749 = m(6) * t141;
t486 = t478 * (-t479 * rSges(6,3) + (rSges(6,1) * t521 - rSges(6,2) * t473) * t478);
t388 = t475 * t486;
t598 = t478 * t482;
t498 = t478 * (pkin(4) * t598 - t479 * pkin(8));
t554 = t475 * t498 + t388;
t602 = t475 * t478;
t340 = rSges(6,3) * t602 - t703;
t630 = t479 * pkin(3);
t631 = t478 * pkin(7);
t354 = (t630 + t631) * t475 + t535;
t740 = t340 - t354;
t201 = t740 * t479 + t554;
t328 = t479 * t340;
t202 = -t354 * t479 + t328 + t554;
t746 = m(6) * (t201 - t202);
t629 = sin(qJ(1)) * pkin(1);
t275 = t278 + t629;
t266 = t475 * t275;
t593 = t479 * t480;
t443 = t474 * t593 + t475 * t482;
t592 = t479 * t482;
t444 = t474 * t592 - t601;
t550 = t444 * rSges(5,1) - t443 * rSges(5,2);
t695 = (rSges(5,3) + pkin(7)) * t478 + pkin(2) + t630;
t308 = t695 * t474 - t470 + t550;
t296 = t308 + t629;
t287 = t475 * t296;
t523 = rSges(4,1) * t479 + pkin(2);
t608 = t474 * t478;
t379 = -rSges(4,2) * t608 - t475 * rSges(4,3) + t474 * t523 - t470;
t377 = t379 + t629;
t366 = t475 * t377;
t380 = rSges(4,2) * t602 - t523 * t475 + (-rSges(4,3) - qJ(3)) * t474;
t628 = cos(qJ(1)) * pkin(1);
t378 = t380 - t628;
t604 = t475 * t379;
t605 = t475 * t308;
t606 = t475 * t278;
t685 = m(6) / 0.2e1;
t686 = m(5) / 0.2e1;
t445 = -t474 * t482 + t475 * t593;
t446 = t475 * t592 + t607;
t701 = -t446 * rSges(5,1) + t445 * rSges(5,2);
t717 = -t695 * t475 - t610 + t701;
t725 = m(4) / 0.2e1;
t728 = t717 - t628;
t729 = t716 - t628;
t537 = (-t266 - t606 + (-t729 - t716) * t474) * t685 + (-t287 - t605 + (-t728 - t717) * t474) * t686 + (-t366 - t604 + (-t378 - t380) * t474) * t725;
t167 = -t474 * t729 - t266;
t204 = -t474 * t728 - t287;
t571 = t474 * t717 + t605;
t574 = t474 * t716 + t606;
t538 = (t167 + t574) * t685 + (t204 + t571) * t686 + (-t366 + t604 + (-t378 + t380) * t474) * t725;
t18 = t538 - t537;
t748 = t18 * qJD(1);
t138 = -t275 * t364 - t365 * t729;
t747 = m(6) * t138;
t476 = t478 ^ 2;
t745 = t328 + t388;
t387 = t474 * t486;
t339 = -rSges(6,3) * t608 - t553;
t562 = t339 + (-pkin(4) * t592 + t631) * t474 + t545;
t200 = -t474 * t498 + t479 * t562 - t387;
t743 = t200 * t729;
t331 = Icges(6,5) * t425 - Icges(6,6) * t424 + Icges(6,3) * t602;
t739 = t331 * t602;
t114 = -t275 * t716 + t278 * t729;
t136 = -t296 * t717 + t308 * t728;
t356 = rSges(5,3) * t602 - t701;
t492 = t478 * (-t479 * rSges(5,3) + (rSges(5,1) * t482 - rSges(5,2) * t480) * t478);
t738 = t356 * t479 + t475 * t492;
t737 = -t478 / 0.2e1;
t736 = t200 * t716;
t618 = Icges(6,4) * t473;
t418 = -Icges(6,5) * t479 + (Icges(6,1) * t521 - t618) * t478;
t702 = t418 + (-Icges(6,2) * t521 - t618) * t478;
t731 = t473 * t702;
t621 = Icges(5,4) * t480;
t439 = -Icges(5,5) * t479 + (Icges(5,1) * t482 - t621) * t478;
t547 = t439 + (-Icges(5,2) * t482 - t621) * t478;
t730 = t480 * t547;
t407 = Icges(6,4) * t425;
t334 = -Icges(6,2) * t424 + Icges(6,6) * t602 + t407;
t406 = Icges(6,4) * t424;
t338 = -Icges(6,1) * t425 - Icges(6,5) * t602 + t406;
t568 = t496 * t334 + t423 * t338;
t428 = Icges(5,4) * t446;
t348 = -Icges(5,2) * t445 + Icges(5,6) * t602 + t428;
t427 = Icges(5,4) * t445;
t352 = -Icges(5,1) * t446 - Icges(5,5) * t602 + t427;
t566 = t443 * t348 + t352 * t444;
t510 = -t445 * t348 - t446 * t352;
t665 = m(3) * (-t628 * (rSges(3,1) * t474 + rSges(3,2) * t475) - (-t475 * rSges(3,1) + t474 * rSges(3,2)) * t629);
t661 = m(4) * (-t380 * t377 + t378 * t379);
t162 = -t331 * t608 + t568;
t330 = -Icges(6,5) * t423 + Icges(6,6) * t496 - Icges(6,3) * t608;
t619 = Icges(6,4) * t423;
t333 = Icges(6,2) * t496 - Icges(6,6) * t608 - t619;
t405 = Icges(6,4) * t496;
t336 = -Icges(6,1) * t423 - Icges(6,5) * t608 + t405;
t567 = -t496 * t333 + t423 * t336;
t507 = t330 * t608 + t567;
t694 = t474 ^ 2;
t724 = (-t331 * t694 * t478 + (-t162 + t568) * t474 + (t507 - t567 + (-t330 * t474 - t331 * t475) * t478 + t739) * t475) * t478;
t437 = -Icges(5,3) * t479 + (Icges(5,5) * t482 - Icges(5,6) * t480) * t478;
t620 = Icges(5,4) * t482;
t438 = -Icges(5,6) * t479 + (-Icges(5,2) * t480 + t620) * t478;
t722 = (t437 * t602 - t438 * t445 + t439 * t446) * t479;
t430 = t445 * pkin(4);
t323 = -t365 + t430;
t355 = -rSges(5,3) * t608 - t550;
t298 = t479 * t355 - t474 * t492;
t590 = (-t200 * t475 - t202 * t474) * t685 + (-t298 * t475 - t474 * t738) * t686;
t623 = t474 * t746 / 0.2e1;
t714 = t590 - t623;
t588 = t202 * t278 - t736;
t508 = (-t201 * t278 + t588 + t736) * t685;
t589 = -t201 * t275 + t743;
t624 = (t588 + t589) * t685 + ((-t296 + t308) * t738 + (-t717 + t728) * t298) * t686;
t713 = t624 - t508;
t450 = (-Icges(5,5) * t480 - Icges(5,6) * t482) * t478;
t594 = t479 * t450;
t712 = -t594 / 0.2e1 + t730 * t737;
t440 = (-Icges(6,5) * t473 - Icges(6,6) * t521) * t478;
t595 = t479 * t440;
t711 = -t595 / 0.2e1 + t731 * t737;
t155 = (-t740 * t474 - t562 * t475) * t478;
t263 = (-t364 * t475 - t365 * t474) * t478;
t447 = (-rSges(6,1) * t473 - rSges(6,2) * t521) * t478;
t401 = t447 * t608;
t310 = -t364 * t479 + t401;
t311 = t479 * t365 + t447 * t602;
t358 = Icges(6,5) * t496 + Icges(6,6) * t423;
t557 = -Icges(6,1) * t496 + t333 - t619;
t705 = Icges(6,2) * t423 + t336 + t405;
t142 = -t479 * t358 + (-t473 * t705 - t557 * t521) * t478;
t359 = -Icges(6,5) * t424 - Icges(6,6) * t425;
t556 = Icges(6,1) * t424 + t334 + t407;
t704 = -Icges(6,2) * t425 - t338 - t406;
t143 = -t479 * t359 + (-t473 * t704 - t556 * t521) * t478;
t515 = t521 * Icges(6,4);
t417 = -Icges(6,6) * t479 + (-Icges(6,2) * t473 + t515) * t478;
t442 = (-Icges(6,1) * t473 - t515) * t478;
t546 = -t442 + t417;
t175 = t423 * t546 - t440 * t608 + t702 * t496;
t176 = -t424 * t702 - t425 * t546 + t440 * t602;
t528 = t602 / 0.2e1;
t531 = -t608 / 0.2e1;
t617 = (-t595 + (-t546 * t521 - t731) * t478) * t479;
t668 = -t479 / 0.2e1;
t543 = (-t617 + (-t142 * t474 + t143 * t475) * t478) * t668 + (-t175 * t479 - (-t358 * t608 + t423 * t557) * t608 + (-t359 * t608 + t423 * t556) * t602 + (t704 * t602 - t705 * t608) * t496) * t531 + (-t176 * t479 - (t358 * t602 - t424 * t705 - t425 * t557) * t608 + (t359 * t602 - t424 * t704 - t425 * t556) * t602) * t528;
t14 = t543 + m(6) * (t155 * t263 - t200 * t310 + t202 * t311);
t709 = t14 * qJD(5);
t544 = qJD(1) + qJD(2);
t344 = -Icges(5,5) * t444 + Icges(5,6) * t443 - Icges(5,3) * t608;
t622 = Icges(5,4) * t444;
t347 = Icges(5,2) * t443 - Icges(5,6) * t608 - t622;
t426 = Icges(5,4) * t443;
t350 = -Icges(5,1) * t444 - Icges(5,5) * t608 + t426;
t698 = t344 * t602 - t445 * t347 + t446 * t350;
t375 = rSges(5,1) * t443 + rSges(5,2) * t444;
t376 = -rSges(5,1) * t445 - rSges(5,2) * t446;
t555 = -t443 * pkin(4) - t364;
t577 = (t323 * t475 - t474 * t555) * t685 + (t375 * t474 - t376 * t475) * t686;
t129 = t275 * t555 + t323 * t729;
t134 = t278 * t555 + t323 * t716;
t154 = -t296 * t375 - t376 * t728;
t160 = -t308 * t375 - t376 * t717;
t697 = (t160 + t154) * t686 + (t134 + t129) * t685;
t696 = t200 * t746;
t692 = 0.4e1 * qJD(1);
t691 = 0.2e1 * qJD(2);
t690 = 0.4e1 * qJD(2);
t689 = 2 * qJD(4);
t283 = t339 * t479 - t387;
t682 = t283 * t746;
t583 = -t310 * t275 + t311 * t729;
t586 = -t200 * t364 - t202 * t365;
t678 = m(6) * (t583 + t586);
t582 = -t310 * t278 + t311 * t716;
t677 = m(6) * (t582 + t586);
t676 = m(6) * ((-t275 + t278) * t745 + (-t716 + t729) * t283);
t578 = t283 * t555 + t323 * t745;
t673 = m(6) * (t578 + t583);
t671 = m(6) * (t578 + t582);
t669 = m(6) * (t141 + t138);
t345 = Icges(5,5) * t446 - Icges(5,6) * t445 + Icges(5,3) * t602;
t181 = -t345 * t608 + t566;
t565 = -t443 * t347 + t444 * t350;
t506 = t344 * t608 + t565;
t614 = (-t437 * t608 + t438 * t443 - t439 * t444) * t479;
t100 = -t614 + (t181 * t475 + t474 * t506) * t478;
t183 = t345 * t602 + t510;
t101 = -t722 + (t183 * t475 - t698 * t474) * t478;
t41 = -t614 + (t566 * t475 + (-t183 + t506 + t510) * t474) * t478;
t42 = t722 + (-(-t566 - t698) * t474 + t506 * t475 + (-t510 - t565) * t475 - t474 * t181 + (-t345 * t694 + (-t344 * t474 - t345 * t475) * t475) * t478) * t478;
t616 = (-(-Icges(6,3) * t479 + (Icges(6,5) * t521 - Icges(6,6) * t473) * t478) * t608 + t417 * t496 - t423 * t418) * t479;
t37 = -t616 + (t568 * t475 + (t507 - t739) * t474) * t478;
t529 = -t602 / 0.2e1;
t93 = -t616 + (t162 * t475 + t474 * t507) * t478;
t520 = t37 * t528 + t93 * t529 + t724 * t531;
t2 = ((-t100 / 0.2e1 + t41 / 0.2e1) * t475 + (-t101 / 0.2e1 - t42 / 0.2e1) * t474) * t478 + t520 - t696;
t29 = t590 + t623 - t577;
t666 = t29 * qJD(3) + t2 * qJD(4);
t659 = m(4) * (-t378 * t474 - t366);
t658 = m(4) * (-t380 * t474 - t604);
t653 = m(5) * t136;
t651 = m(5) * t154;
t650 = m(5) * t160;
t649 = m(5) * t204;
t647 = m(5) * t571;
t643 = m(6) * t114;
t639 = m(6) * t129;
t638 = m(6) * t134;
t637 = m(6) * t167;
t636 = m(6) * t574;
t635 = m(6) * (-t283 * t475 - t474 * t745);
t634 = m(6) * (t310 * t474 + t311 * t475);
t632 = m(6) * (t364 * t474 - t365 * t475);
t179 = t635 / 0.2e1;
t74 = t179 - t632 / 0.2e1;
t627 = t74 * qJD(3) + qJD(5) * t520;
t27 = t577 - t714;
t265 = t632 / 0.2e1;
t72 = t265 - t635 / 0.2e1;
t626 = t27 * qJD(4) + t72 * qJD(5);
t28 = t577 + t714;
t73 = t265 + t179;
t625 = t28 * qJD(4) + t73 * qJD(5);
t452 = (-Icges(5,1) * t480 - t620) * t478;
t548 = t438 - t452;
t615 = (-t594 + (-t482 * t548 - t730) * t478) * t479;
t587 = t200 * t555 + t202 * t323;
t576 = -t283 * t364 - t365 * t745;
t575 = -t298 * t375 - t376 * t738;
t561 = -Icges(5,1) * t443 + t347 - t622;
t560 = Icges(5,1) * t445 + t348 + t428;
t559 = Icges(5,2) * t444 + t350 + t426;
t558 = -Icges(5,2) * t446 - t352 - t427;
t530 = t608 / 0.2e1;
t525 = -t598 / 0.2e1;
t524 = t598 / 0.2e1;
t517 = t478 * t521;
t504 = t517 / 0.2e1;
t505 = -t517 / 0.2e1;
t512 = t417 * t505 + t442 * t504 + t711;
t509 = -t682 / 0.2e1 + t520;
t503 = t669 / 0.2e1 + t512;
t495 = t417 * t504 + t442 * t505 - t711;
t494 = t438 * t525 + t452 * t524 + t512 + t712;
t491 = t37 * t529 - t617 + (t142 + t175) * t531 + t724 * t530 + (t93 + t143 + t176) * t528;
t489 = t494 + t697;
t488 = t682 / 0.2e1 + t491;
t485 = t438 * t524 + t452 * t525 + t495 - t712;
t367 = Icges(5,5) * t443 + Icges(5,6) * t444;
t152 = -t479 * t367 + (-t480 * t559 - t482 * t561) * t478;
t368 = -Icges(5,5) * t445 - Icges(5,6) * t446;
t153 = -t479 * t368 + (-t480 * t558 - t482 * t560) * t478;
t215 = t443 * t547 + t444 * t548 - t450 * t608;
t216 = -t445 * t547 - t446 * t548 + t450 * t602;
t484 = t28 * qJD(3) + (t41 * t529 + t491 - t615 + t696 + (t100 + t153 + t216) * t528 + (t152 + t215) * t531 + (t101 + t42) * t530) * qJD(4);
t453 = (-rSges(5,1) * t480 - rSges(5,2) * t482) * t478;
t325 = t376 * t479 + t453 * t602;
t324 = -t375 * t479 + t453 * t608;
t262 = -t430 * t479 - t467 * t476 + t311;
t261 = -t466 * t476 + t479 * t555 + t401;
t247 = (-t339 * t475 - t340 * t474) * t478;
t219 = (t323 * t474 + t555 * t475) * t478;
t209 = qJD(5) * t634;
t107 = t512 + t749;
t106 = t512 + t747;
t105 = -t636 - t647 + t658;
t102 = t637 + t649 + t659;
t78 = t671 / 0.2e1;
t75 = t673 / 0.2e1;
t70 = t73 * qJD(3);
t61 = t676 / 0.2e1;
t59 = t677 / 0.2e1;
t57 = t678 / 0.2e1;
t56 = t494 + t638 + t650;
t55 = t494 + t639 + t651;
t49 = t643 + t653 + t661 + t665;
t22 = -t676 / 0.2e1 + t503;
t21 = t61 + t503;
t20 = t537 + t538;
t17 = t61 - t669 / 0.2e1 + t495;
t16 = m(6) * (t247 * t263 - t283 * t310 + t311 * t745) + t543;
t15 = t16 * qJD(5);
t13 = t489 + t713;
t12 = t489 - t713;
t9 = t485 + t508 + t624 - t697;
t8 = t78 - t677 / 0.2e1 + t509;
t7 = t59 - t671 / 0.2e1 + t509;
t6 = t75 - t678 / 0.2e1 + t509;
t5 = t57 - t673 / 0.2e1 + t509;
t4 = t488 + t59 + t78;
t3 = t488 + t57 + t75;
t1 = [t49 * qJD(2) + t102 * qJD(3) + t55 * qJD(4) + t106 * qJD(5), t49 * qJD(1) + t20 * qJD(3) + t13 * qJD(4) + t21 * qJD(5) + (t665 / 0.2e1 + t661 / 0.2e1 + t136 * t686 + t114 * t685) * t691, qJD(1) * t102 + qJD(2) * t20 + t625, t55 * qJD(1) + t13 * qJD(2) + t3 * qJD(5) + ((-t296 * t324 + t325 * t728 + t575) * t686 + (-t261 * t275 + t262 * t729 + t587) * t685) * t689 + t484, t106 * qJD(1) + t21 * qJD(2) + t70 + t3 * qJD(4) + (m(6) * (t576 + t583) + t491) * qJD(5); -t18 * qJD(3) + t12 * qJD(4) + t22 * qJD(5) + (-t665 / 0.4e1 - t661 / 0.4e1 - t653 / 0.4e1 - t643 / 0.4e1) * t692, qJD(3) * t105 + qJD(4) * t56 + qJD(5) * t107, qJD(2) * t105 + t625 - t748, t12 * qJD(1) + t56 * qJD(2) + t4 * qJD(5) + ((-t261 * t278 + t262 * t716 + t587) * t685 + (-t308 * t324 + t325 * t717 + t575) * t686) * t689 + t484, t22 * qJD(1) + t107 * qJD(2) + t70 + t4 * qJD(4) + (m(6) * (t576 + t582) + t491) * qJD(5); t18 * qJD(2) + (-t659 / 0.4e1 - t649 / 0.4e1 - t637 / 0.4e1) * t692 + t626, t748 + (t636 / 0.4e1 + t647 / 0.4e1 - t658 / 0.4e1) * t690 + t626, 0, ((t324 * t474 + t325 * t475) * t686 + (t261 * t474 + t262 * t475) * t685) * t689 + t209 + t544 * t27, qJD(4) * t634 + t544 * t72 + t209; t9 * qJD(2) + t5 * qJD(5) + (-t651 / 0.4e1 - t639 / 0.4e1) * t692 + t666 + (t485 + 0.2e1 * (t202 * t275 + t589 - t743) * t685) * qJD(1), t9 * qJD(1) + t485 * qJD(2) + t7 * qJD(5) + (-t638 / 0.4e1 - t650 / 0.4e1) * t690 + t508 * t691 + t666, t544 * t29, (m(5) * (-t298 * t324 + t738 * t325 + (-t355 * t475 - t356 * t474) * t476 * (-t375 * t475 - t376 * t474)) + (-t615 + (-t152 * t474 + t153 * t475) * t478) * t668 + (-t215 * t479 - (-t367 * t608 + t443 * t559 + t444 * t561) * t608 + (-t368 * t608 + t443 * t558 + t444 * t560) * t602) * t531 + (-t216 * t479 - (t367 * t602 - t445 * t559 - t446 * t561) * t608 + (t368 * t602 - t445 * t558 - t446 * t560) * t602) * t528 + m(6) * (t155 * t219 - t200 * t261 + t202 * t262) + t543) * qJD(4) + t709 + t544 * t2, t5 * qJD(1) + t7 * qJD(2) + t14 * qJD(4) + t709; (t495 - t747) * qJD(1) + t17 * qJD(2) + t6 * qJD(4) + t627, t17 * qJD(1) + (t495 - t749) * qJD(2) + t8 * qJD(4) + t627, t544 * t74, t6 * qJD(1) + t8 * qJD(2) + ((t219 * t247 - t261 * t283 + t262 * t745) * m(6) + t543) * qJD(4) + t15, qJD(4) * t16 + t520 * t544 + t15;];
Cq = t1;
