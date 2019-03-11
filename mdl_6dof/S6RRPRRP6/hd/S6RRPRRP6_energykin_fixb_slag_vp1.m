% Calculate kinetic energy for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:07:44
% EndTime: 2019-03-09 12:07:47
% DurationCPUTime: 3.23s
% Computational Cost: add. (3500->326), mult. (8808->490), div. (0->0), fcn. (11199->12), ass. (0->155)
t623 = Icges(6,1) + Icges(7,1);
t622 = -Icges(6,4) + Icges(7,5);
t621 = Icges(7,4) + Icges(6,5);
t620 = Icges(6,2) + Icges(7,3);
t619 = Icges(7,2) + Icges(6,3);
t618 = -Icges(6,6) + Icges(7,6);
t617 = rSges(7,1) + pkin(5);
t616 = rSges(7,3) + qJ(6);
t551 = sin(pkin(11));
t556 = sin(qJ(2));
t558 = cos(qJ(2));
t594 = cos(pkin(11));
t536 = -t558 * t551 - t556 * t594;
t557 = sin(qJ(1));
t559 = cos(qJ(1));
t537 = -t556 * t551 + t558 * t594;
t553 = cos(pkin(6));
t566 = t553 * t537;
t513 = t557 * t536 + t559 * t566;
t529 = t536 * t553;
t514 = -t529 * t559 + t537 * t557;
t552 = sin(pkin(6));
t592 = t552 * t559;
t464 = Icges(4,5) * t514 + Icges(4,6) * t513 - Icges(4,3) * t592;
t515 = t536 * t559 - t557 * t566;
t516 = t529 * t557 + t537 * t559;
t593 = t552 * t557;
t465 = Icges(4,5) * t516 + Icges(4,6) * t515 + Icges(4,3) * t593;
t587 = t558 * t559;
t590 = t556 * t557;
t531 = t553 * t587 - t590;
t588 = t557 * t558;
t589 = t556 * t559;
t532 = t553 * t589 + t588;
t533 = -t553 * t588 - t589;
t534 = -t553 * t590 + t587;
t567 = (Icges(3,5) * t532 + Icges(3,6) * t531 - Icges(3,3) * t592) * t559 - (Icges(3,5) * t534 + Icges(3,6) * t533 + Icges(3,3) * t593) * t557;
t615 = t552 * (t464 * t559 - t465 * t557 + t567);
t555 = sin(qJ(4));
t598 = cos(qJ(4));
t486 = t514 * t598 - t555 * t592;
t554 = sin(qJ(5));
t597 = cos(qJ(5));
t458 = t486 * t554 + t513 * t597;
t459 = t486 * t597 - t513 * t554;
t578 = t552 * t598;
t485 = t514 * t555 + t559 * t578;
t614 = t620 * t458 + t622 * t459 + t618 * t485;
t488 = t516 * t598 + t555 * t593;
t460 = t488 * t554 + t515 * t597;
t461 = t488 * t597 - t515 * t554;
t487 = t516 * t555 - t557 * t578;
t613 = t620 * t460 + t622 * t461 + t618 * t487;
t612 = t618 * t458 + t621 * t459 + t619 * t485;
t611 = t618 * t460 + t621 * t461 + t619 * t487;
t610 = t622 * t458 + t623 * t459 + t621 * t485;
t609 = t622 * t460 + t623 * t461 + t621 * t487;
t528 = t536 * t552;
t518 = -t528 * t598 + t553 * t555;
t527 = t537 * t552;
t483 = t518 * t554 + t527 * t597;
t484 = t518 * t597 - t527 * t554;
t517 = -t528 * t555 - t553 * t598;
t608 = t620 * t483 + t622 * t484 + t618 * t517;
t607 = t618 * t483 + t621 * t484 + t619 * t517;
t606 = t622 * t483 + t623 * t484 + t621 * t517;
t605 = -Icges(4,5) * t528 + Icges(4,6) * t527 + (Icges(3,5) * t556 + Icges(3,6) * t558) * t552 + (Icges(4,3) + Icges(3,3)) * t553;
t596 = pkin(2) * t556;
t595 = pkin(2) * t558;
t586 = rSges(7,2) * t485 + t616 * t458 + t459 * t617;
t585 = rSges(7,2) * t487 + t616 * t460 + t461 * t617;
t584 = rSges(7,2) * t517 + t616 * t483 + t484 * t617;
t576 = -qJ(3) * t552 + t553 * t596;
t512 = -t557 * t576 + t559 * t595;
t535 = qJD(1) * (pkin(1) * t559 + pkin(8) * t593);
t546 = qJD(2) * t553 + qJD(1);
t583 = t546 * t512 + t535;
t581 = qJD(2) * t552;
t545 = t557 * t581;
t489 = -qJD(4) * t515 + t545;
t582 = qJD(1) * (pkin(1) * t557 - pkin(8) * t592);
t580 = qJD(3) * t559;
t511 = t557 * t595 + t559 * t576;
t577 = t559 * t581;
t579 = qJD(3) * t553 + t511 * t545 + t512 * t577;
t519 = -qJD(4) * t527 + t546;
t538 = qJ(3) * t553 + t552 * t596;
t574 = qJD(2) * (rSges(4,1) * t528 - rSges(4,2) * t527 - rSges(4,3) * t553 - t538);
t573 = qJD(2) * (pkin(3) * t528 + pkin(9) * t527 - t538);
t572 = qJD(3) * t593 - t582;
t490 = -qJD(4) * t513 - t577;
t479 = pkin(3) * t514 - pkin(9) * t513;
t480 = pkin(3) * t516 - pkin(9) * t515;
t569 = t479 * t545 + t480 * t577 + t579;
t453 = pkin(4) * t486 + pkin(10) * t485;
t454 = pkin(4) * t488 + pkin(10) * t487;
t565 = t489 * t453 - t454 * t490 + t569;
t564 = t546 * t480 + (t557 * t573 - t580) * t552 + t583;
t563 = (-t479 - t511) * t546 + t573 * t592 + t572;
t481 = pkin(4) * t518 + pkin(10) * t517;
t562 = t519 * t454 - t481 * t489 + t564;
t561 = -t453 * t519 + t490 * t481 + t563;
t541 = rSges(2,1) * t559 - rSges(2,2) * t557;
t540 = rSges(2,1) * t557 + rSges(2,2) * t559;
t526 = rSges(3,3) * t553 + (rSges(3,1) * t556 + rSges(3,2) * t558) * t552;
t525 = Icges(3,5) * t553 + (Icges(3,1) * t556 + Icges(3,4) * t558) * t552;
t524 = Icges(3,6) * t553 + (Icges(3,4) * t556 + Icges(3,2) * t558) * t552;
t507 = rSges(3,1) * t534 + rSges(3,2) * t533 + rSges(3,3) * t593;
t506 = rSges(3,1) * t532 + rSges(3,2) * t531 - rSges(3,3) * t592;
t503 = Icges(3,1) * t534 + Icges(3,4) * t533 + Icges(3,5) * t593;
t502 = Icges(3,1) * t532 + Icges(3,4) * t531 - Icges(3,5) * t592;
t501 = Icges(3,4) * t534 + Icges(3,2) * t533 + Icges(3,6) * t593;
t500 = Icges(3,4) * t532 + Icges(3,2) * t531 - Icges(3,6) * t592;
t496 = -Icges(4,1) * t528 + Icges(4,4) * t527 + Icges(4,5) * t553;
t495 = -Icges(4,4) * t528 + Icges(4,2) * t527 + Icges(4,6) * t553;
t482 = qJD(5) * t517 + t519;
t478 = rSges(5,1) * t518 - rSges(5,2) * t517 - rSges(5,3) * t527;
t477 = Icges(5,1) * t518 - Icges(5,4) * t517 - Icges(5,5) * t527;
t476 = Icges(5,4) * t518 - Icges(5,2) * t517 - Icges(5,6) * t527;
t475 = Icges(5,5) * t518 - Icges(5,6) * t517 - Icges(5,3) * t527;
t472 = rSges(4,1) * t516 + rSges(4,2) * t515 + rSges(4,3) * t593;
t471 = rSges(4,1) * t514 + rSges(4,2) * t513 - rSges(4,3) * t592;
t469 = Icges(4,1) * t516 + Icges(4,4) * t515 + Icges(4,5) * t593;
t468 = Icges(4,1) * t514 + Icges(4,4) * t513 - Icges(4,5) * t592;
t467 = Icges(4,4) * t516 + Icges(4,2) * t515 + Icges(4,6) * t593;
t466 = Icges(4,4) * t514 + Icges(4,2) * t513 - Icges(4,6) * t592;
t463 = t507 * t546 - t526 * t545 + t535;
t462 = -t506 * t546 - t526 * t577 - t582;
t457 = (t506 * t557 + t507 * t559) * t581;
t456 = qJD(5) * t485 + t490;
t455 = qJD(5) * t487 + t489;
t449 = rSges(6,1) * t484 - rSges(6,2) * t483 + rSges(6,3) * t517;
t441 = rSges(5,1) * t488 - rSges(5,2) * t487 - rSges(5,3) * t515;
t440 = rSges(5,1) * t486 - rSges(5,2) * t485 - rSges(5,3) * t513;
t439 = Icges(5,1) * t488 - Icges(5,4) * t487 - Icges(5,5) * t515;
t438 = Icges(5,1) * t486 - Icges(5,4) * t485 - Icges(5,5) * t513;
t437 = Icges(5,4) * t488 - Icges(5,2) * t487 - Icges(5,6) * t515;
t436 = Icges(5,4) * t486 - Icges(5,2) * t485 - Icges(5,6) * t513;
t435 = Icges(5,5) * t488 - Icges(5,6) * t487 - Icges(5,3) * t515;
t434 = Icges(5,5) * t486 - Icges(5,6) * t485 - Icges(5,3) * t513;
t430 = t472 * t546 + (t557 * t574 - t580) * t552 + t583;
t429 = (-t471 - t511) * t546 + t574 * t592 + t572;
t428 = rSges(6,1) * t461 - rSges(6,2) * t460 + rSges(6,3) * t487;
t426 = rSges(6,1) * t459 - rSges(6,2) * t458 + rSges(6,3) * t485;
t412 = (t471 * t557 + t472 * t559) * t581 + t579;
t411 = t441 * t519 - t478 * t489 + t564;
t410 = -t440 * t519 + t478 * t490 + t563;
t409 = t440 * t489 - t441 * t490 + t569;
t408 = t428 * t482 - t449 * t455 + t562;
t407 = -t426 * t482 + t449 * t456 + t561;
t406 = t426 * t455 - t428 * t456 + t565;
t405 = qJD(6) * t458 - t455 * t584 + t482 * t585 + t562;
t404 = qJD(6) * t460 + t456 * t584 - t482 * t586 + t561;
t403 = qJD(6) * t483 + t455 * t586 - t456 * t585 + t565;
t1 = m(5) * (t409 ^ 2 + t410 ^ 2 + t411 ^ 2) / 0.2e1 + m(3) * (t457 ^ 2 + t462 ^ 2 + t463 ^ 2) / 0.2e1 + m(4) * (t412 ^ 2 + t429 ^ 2 + t430 ^ 2) / 0.2e1 + m(7) * (t403 ^ 2 + t404 ^ 2 + t405 ^ 2) / 0.2e1 + m(6) * (t406 ^ 2 + t407 ^ 2 + t408 ^ 2) / 0.2e1 + t489 * ((-t435 * t515 - t437 * t487 + t439 * t488) * t489 + (-t434 * t515 - t436 * t487 + t438 * t488) * t490 + (-t475 * t515 - t476 * t487 + t477 * t488) * t519) / 0.2e1 + t490 * ((-t435 * t513 - t437 * t485 + t439 * t486) * t489 + (-t434 * t513 - t436 * t485 + t438 * t486) * t490 + (-t475 * t513 - t476 * t485 + t477 * t486) * t519) / 0.2e1 + t519 * ((-t435 * t527 - t437 * t517 + t439 * t518) * t489 + (-t434 * t527 - t436 * t517 + t438 * t518) * t490 + (-t475 * t527 - t476 * t517 + t477 * t518) * t519) / 0.2e1 + ((t460 * t608 + t461 * t606 + t487 * t607) * t482 + (t460 * t614 + t610 * t461 + t612 * t487) * t456 + (t613 * t460 + t609 * t461 + t611 * t487) * t455) * t455 / 0.2e1 + ((t458 * t608 + t459 * t606 + t485 * t607) * t482 + (t614 * t458 + t610 * t459 + t612 * t485) * t456 + (t458 * t613 + t459 * t609 + t485 * t611) * t455) * t456 / 0.2e1 + ((t608 * t483 + t606 * t484 + t607 * t517) * t482 + (t483 * t614 + t610 * t484 + t612 * t517) * t456 + (t483 * t613 + t484 * t609 + t517 * t611) * t455) * t482 / 0.2e1 + (((t465 * t553 + t467 * t527 - t469 * t528) * t557 - (t464 * t553 + t466 * t527 - t468 * t528) * t559 + ((t501 * t558 + t503 * t556) * t557 - (t500 * t558 + t502 * t556) * t559) * t552 - t567 * t553) * t581 + (t495 * t527 - t496 * t528 + (t524 * t558 + t525 * t556) * t552 + t605 * t553) * t546) * t546 / 0.2e1 + (Icges(2,3) + m(2) * (t540 ^ 2 + t541 ^ 2)) * qJD(1) ^ 2 / 0.2e1 - (((-t513 * t466 - t514 * t468 - t531 * t500 - t532 * t502 + t615) * t559 + (t467 * t513 + t469 * t514 + t501 * t531 + t503 * t532) * t557) * t581 + (t495 * t513 + t496 * t514 + t524 * t531 + t525 * t532 - t592 * t605) * t546) * t577 / 0.2e1 + (((-t466 * t515 - t468 * t516 - t500 * t533 - t502 * t534) * t559 + (t515 * t467 + t516 * t469 + t533 * t501 + t534 * t503 - t615) * t557) * t581 + (t495 * t515 + t496 * t516 + t524 * t533 + t525 * t534 + t593 * t605) * t546) * t545 / 0.2e1;
T  = t1;
