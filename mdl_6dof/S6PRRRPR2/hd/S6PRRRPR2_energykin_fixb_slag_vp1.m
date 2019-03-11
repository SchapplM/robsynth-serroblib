% Calculate kinetic energy for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:05:20
% EndTime: 2019-03-08 23:05:23
% DurationCPUTime: 3.01s
% Computational Cost: add. (3344->376), mult. (5851->570), div. (0->0), fcn. (7074->14), ass. (0->170)
t575 = Icges(5,2) + Icges(6,3);
t524 = sin(pkin(11));
t527 = cos(pkin(11));
t533 = cos(qJ(2));
t528 = cos(pkin(6));
t531 = sin(qJ(2));
t555 = t528 * t531;
t505 = t524 * t533 + t527 * t555;
t553 = qJ(3) + qJ(4);
t521 = sin(t553);
t544 = cos(t553);
t525 = sin(pkin(6));
t561 = t525 * t527;
t483 = t505 * t544 - t521 * t561;
t554 = t528 * t533;
t504 = t524 * t531 - t527 * t554;
t523 = sin(pkin(12));
t526 = cos(pkin(12));
t446 = -t483 * t523 + t504 * t526;
t564 = t504 * t523;
t447 = t483 * t526 + t564;
t543 = t525 * t544;
t482 = t505 * t521 + t527 * t543;
t574 = -Icges(5,4) * t483 + Icges(6,5) * t447 - Icges(5,6) * t504 + Icges(6,6) * t446 + t482 * t575;
t507 = -t524 * t555 + t527 * t533;
t562 = t524 * t525;
t485 = t507 * t544 + t521 * t562;
t506 = t524 * t554 + t527 * t531;
t448 = -t485 * t523 + t506 * t526;
t563 = t506 * t523;
t449 = t485 * t526 + t563;
t484 = t507 * t521 - t524 * t543;
t573 = -Icges(5,4) * t485 + Icges(6,5) * t449 - Icges(5,6) * t506 + Icges(6,6) * t448 + t484 * t575;
t499 = t528 * t521 + t531 * t543;
t557 = t525 * t533;
t480 = -t499 * t523 - t526 * t557;
t549 = t523 * t557;
t481 = t499 * t526 - t549;
t559 = t525 * t531;
t498 = t521 * t559 - t528 * t544;
t572 = -Icges(5,4) * t499 + Icges(6,5) * t481 + Icges(5,6) * t557 + Icges(6,6) * t480 + t498 * t575;
t571 = qJD(2) ^ 2;
t532 = cos(qJ(3));
t567 = pkin(3) * t532;
t566 = pkin(5) * t526;
t530 = sin(qJ(3));
t560 = t525 * t530;
t558 = t525 * t532;
t556 = t528 * t530;
t551 = qJD(2) * t525;
t515 = t524 * t551;
t490 = qJD(3) * t506 + t515;
t518 = qJD(2) * t528;
t548 = t524 * t560;
t547 = t527 * t560;
t454 = qJD(4) * t506 + t490;
t546 = t527 * t551;
t476 = t505 * pkin(2) + t504 * pkin(8);
t477 = t507 * pkin(2) + t506 * pkin(8);
t545 = t476 * t515 + t477 * t546 + qJD(1);
t491 = qJD(3) * t504 - t546;
t510 = (pkin(2) * t531 - pkin(8) * t533) * t525;
t542 = t477 * t518 - t510 * t515;
t455 = qJD(4) * t504 + t491;
t493 = t518 + (-qJD(3) - qJD(4)) * t557;
t425 = -pkin(3) * t547 + pkin(9) * t504 + t505 * t567;
t426 = pkin(3) * t548 + pkin(9) * t506 + t507 * t567;
t541 = t490 * t425 - t426 * t491 + t545;
t540 = (-t476 * t528 - t510 * t561) * qJD(2);
t439 = pkin(4) * t483 + qJ(5) * t482;
t539 = qJD(5) * t498 + t454 * t439 + t541;
t471 = pkin(3) * t556 + (-pkin(9) * t533 + t531 * t567) * t525;
t511 = -qJD(3) * t557 + t518;
t538 = t511 * t426 - t471 * t490 + t542;
t440 = pkin(4) * t485 + qJ(5) * t484;
t537 = qJD(5) * t482 + t493 * t440 + t538;
t536 = -t425 * t511 + t491 * t471 + t540;
t469 = pkin(4) * t499 + qJ(5) * t498;
t535 = qJD(5) * t484 + t455 * t469 + t536;
t522 = pkin(12) + qJ(6);
t520 = cos(t522);
t519 = sin(t522);
t509 = t531 * t558 + t556;
t508 = t528 * t532 - t530 * t559;
t497 = rSges(3,3) * t528 + (rSges(3,1) * t531 + rSges(3,2) * t533) * t525;
t496 = Icges(3,5) * t528 + (Icges(3,1) * t531 + Icges(3,4) * t533) * t525;
t495 = Icges(3,6) * t528 + (Icges(3,4) * t531 + Icges(3,2) * t533) * t525;
t494 = Icges(3,3) * t528 + (Icges(3,5) * t531 + Icges(3,6) * t533) * t525;
t489 = t507 * t532 + t548;
t488 = -t507 * t530 + t524 * t558;
t487 = t505 * t532 - t547;
t486 = -t505 * t530 - t527 * t558;
t475 = t499 * t520 - t519 * t557;
t474 = -t499 * t519 - t520 * t557;
t472 = qJD(6) * t498 + t493;
t470 = rSges(4,1) * t509 + rSges(4,2) * t508 - rSges(4,3) * t557;
t468 = Icges(4,1) * t509 + Icges(4,4) * t508 - Icges(4,5) * t557;
t467 = Icges(4,4) * t509 + Icges(4,2) * t508 - Icges(4,6) * t557;
t466 = Icges(4,5) * t509 + Icges(4,6) * t508 - Icges(4,3) * t557;
t463 = rSges(3,1) * t507 - rSges(3,2) * t506 + rSges(3,3) * t562;
t462 = rSges(3,1) * t505 - rSges(3,2) * t504 - rSges(3,3) * t561;
t461 = Icges(3,1) * t507 - Icges(3,4) * t506 + Icges(3,5) * t562;
t460 = Icges(3,1) * t505 - Icges(3,4) * t504 - Icges(3,5) * t561;
t459 = Icges(3,4) * t507 - Icges(3,2) * t506 + Icges(3,6) * t562;
t458 = Icges(3,4) * t505 - Icges(3,2) * t504 - Icges(3,6) * t561;
t457 = Icges(3,5) * t507 - Icges(3,6) * t506 + Icges(3,3) * t562;
t456 = Icges(3,5) * t505 - Icges(3,6) * t504 - Icges(3,3) * t561;
t453 = rSges(5,1) * t499 - rSges(5,2) * t498 - rSges(5,3) * t557;
t452 = Icges(5,1) * t499 - Icges(5,4) * t498 - Icges(5,5) * t557;
t450 = Icges(5,5) * t499 - Icges(5,6) * t498 - Icges(5,3) * t557;
t445 = t485 * t520 + t506 * t519;
t444 = -t485 * t519 + t506 * t520;
t443 = t483 * t520 + t504 * t519;
t442 = -t483 * t519 + t504 * t520;
t438 = (-t462 * t528 - t497 * t561) * qJD(2);
t437 = (t463 * t528 - t497 * t562) * qJD(2);
t436 = qJD(6) * t482 + t455;
t435 = qJD(6) * t484 + t454;
t434 = rSges(4,1) * t489 + rSges(4,2) * t488 + rSges(4,3) * t506;
t433 = rSges(4,1) * t487 + rSges(4,2) * t486 + rSges(4,3) * t504;
t432 = Icges(4,1) * t489 + Icges(4,4) * t488 + Icges(4,5) * t506;
t431 = Icges(4,1) * t487 + Icges(4,4) * t486 + Icges(4,5) * t504;
t430 = Icges(4,4) * t489 + Icges(4,2) * t488 + Icges(4,6) * t506;
t429 = Icges(4,4) * t487 + Icges(4,2) * t486 + Icges(4,6) * t504;
t428 = Icges(4,5) * t489 + Icges(4,6) * t488 + Icges(4,3) * t506;
t427 = Icges(4,5) * t487 + Icges(4,6) * t486 + Icges(4,3) * t504;
t424 = rSges(5,1) * t485 - rSges(5,2) * t484 + rSges(5,3) * t506;
t423 = rSges(5,1) * t483 - rSges(5,2) * t482 + rSges(5,3) * t504;
t422 = Icges(5,1) * t485 - Icges(5,4) * t484 + Icges(5,5) * t506;
t421 = Icges(5,1) * t483 - Icges(5,4) * t482 + Icges(5,5) * t504;
t418 = Icges(5,5) * t485 - Icges(5,6) * t484 + Icges(5,3) * t506;
t417 = Icges(5,5) * t483 - Icges(5,6) * t482 + Icges(5,3) * t504;
t415 = rSges(6,1) * t481 + rSges(6,2) * t480 + rSges(6,3) * t498;
t413 = Icges(6,1) * t481 + Icges(6,4) * t480 + Icges(6,5) * t498;
t412 = Icges(6,4) * t481 + Icges(6,2) * t480 + Icges(6,6) * t498;
t410 = rSges(7,1) * t475 + rSges(7,2) * t474 + rSges(7,3) * t498;
t409 = Icges(7,1) * t475 + Icges(7,4) * t474 + Icges(7,5) * t498;
t408 = Icges(7,4) * t475 + Icges(7,2) * t474 + Icges(7,6) * t498;
t407 = Icges(7,5) * t475 + Icges(7,6) * t474 + Icges(7,3) * t498;
t405 = -pkin(5) * t549 + pkin(10) * t498 + t499 * t566;
t404 = qJD(1) + (t462 * t524 + t463 * t527) * t551;
t401 = rSges(6,1) * t449 + rSges(6,2) * t448 + rSges(6,3) * t484;
t400 = rSges(6,1) * t447 + rSges(6,2) * t446 + rSges(6,3) * t482;
t399 = Icges(6,1) * t449 + Icges(6,4) * t448 + Icges(6,5) * t484;
t398 = Icges(6,1) * t447 + Icges(6,4) * t446 + Icges(6,5) * t482;
t397 = Icges(6,4) * t449 + Icges(6,2) * t448 + Icges(6,6) * t484;
t396 = Icges(6,4) * t447 + Icges(6,2) * t446 + Icges(6,6) * t482;
t393 = rSges(7,1) * t445 + rSges(7,2) * t444 + rSges(7,3) * t484;
t392 = rSges(7,1) * t443 + rSges(7,2) * t442 + rSges(7,3) * t482;
t391 = Icges(7,1) * t445 + Icges(7,4) * t444 + Icges(7,5) * t484;
t390 = Icges(7,1) * t443 + Icges(7,4) * t442 + Icges(7,5) * t482;
t389 = Icges(7,4) * t445 + Icges(7,2) * t444 + Icges(7,6) * t484;
t388 = Icges(7,4) * t443 + Icges(7,2) * t442 + Icges(7,6) * t482;
t387 = Icges(7,5) * t445 + Icges(7,6) * t444 + Icges(7,3) * t484;
t386 = Icges(7,5) * t443 + Icges(7,6) * t442 + Icges(7,3) * t482;
t385 = pkin(5) * t563 + pkin(10) * t484 + t485 * t566;
t384 = pkin(5) * t564 + pkin(10) * t482 + t483 * t566;
t383 = -t433 * t511 + t470 * t491 + t540;
t382 = t434 * t511 - t470 * t490 + t542;
t381 = t433 * t490 - t434 * t491 + t545;
t380 = -t423 * t493 + t453 * t455 + t536;
t379 = t424 * t493 - t453 * t454 + t538;
t378 = t423 * t454 - t424 * t455 + t541;
t377 = t415 * t455 + (-t400 - t439) * t493 + t535;
t376 = t401 * t493 + (-t415 - t469) * t454 + t537;
t375 = t400 * t454 + (-t401 - t440) * t455 + t539;
t374 = -t392 * t472 + t405 * t455 + t410 * t436 + (-t384 - t439) * t493 + t535;
t373 = t385 * t493 + t393 * t472 - t410 * t435 + (-t405 - t469) * t454 + t537;
t372 = t384 * t454 + t392 * t435 - t393 * t436 + (-t385 - t440) * t455 + t539;
t1 = -t571 * ((-t457 * t561 - t459 * t504 + t461 * t505) * t562 - (-t456 * t561 - t458 * t504 + t460 * t505) * t561 + (-t494 * t561 - t495 * t504 + t496 * t505) * t528) * t561 / 0.2e1 + t490 * ((t428 * t506 + t430 * t488 + t432 * t489) * t490 + (t427 * t506 + t429 * t488 + t431 * t489) * t491 + (t466 * t506 + t467 * t488 + t468 * t489) * t511) / 0.2e1 + t511 * ((-t428 * t557 + t430 * t508 + t432 * t509) * t490 + (-t427 * t557 + t429 * t508 + t431 * t509) * t491 + (-t466 * t557 + t467 * t508 + t468 * t509) * t511) / 0.2e1 + t491 * ((t428 * t504 + t430 * t486 + t432 * t487) * t490 + (t427 * t504 + t429 * t486 + t431 * t487) * t491 + (t466 * t504 + t467 * t486 + t468 * t487) * t511) / 0.2e1 + t472 * ((t387 * t498 + t389 * t474 + t391 * t475) * t435 + (t386 * t498 + t388 * t474 + t390 * t475) * t436 + (t407 * t498 + t408 * t474 + t409 * t475) * t472) / 0.2e1 + t436 * ((t387 * t482 + t389 * t442 + t391 * t443) * t435 + (t482 * t386 + t442 * t388 + t443 * t390) * t436 + (t407 * t482 + t408 * t442 + t409 * t443) * t472) / 0.2e1 + t435 * ((t484 * t387 + t444 * t389 + t445 * t391) * t435 + (t386 * t484 + t388 * t444 + t390 * t445) * t436 + (t407 * t484 + t408 * t444 + t409 * t445) * t472) / 0.2e1 + m(5) * (t378 ^ 2 + t379 ^ 2 + t380 ^ 2) / 0.2e1 + m(4) * (t381 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + m(3) * (t404 ^ 2 + t437 ^ 2 + t438 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(7) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(6) * (t375 ^ 2 + t376 ^ 2 + t377 ^ 2) / 0.2e1 + ((t412 * t448 + t413 * t449 + t450 * t506 + t452 * t485 + t572 * t484) * t493 + (t396 * t448 + t398 * t449 + t417 * t506 + t421 * t485 + t574 * t484) * t455 + (t448 * t397 + t449 * t399 + t506 * t418 + t485 * t422 + t573 * t484) * t454) * t454 / 0.2e1 + ((t412 * t446 + t413 * t447 + t450 * t504 + t452 * t483 + t572 * t482) * t493 + (t446 * t396 + t447 * t398 + t504 * t417 + t483 * t421 + t574 * t482) * t455 + (t397 * t446 + t399 * t447 + t418 * t504 + t422 * t483 + t573 * t482) * t454) * t455 / 0.2e1 + ((t412 * t480 + t413 * t481 - t450 * t557 + t452 * t499 + t572 * t498) * t493 + (t396 * t480 + t398 * t481 - t417 * t557 + t421 * t499 + t574 * t498) * t455 + (t397 * t480 + t399 * t481 - t418 * t557 + t422 * t499 + t573 * t498) * t454) * t493 / 0.2e1 + (((t457 * t562 - t459 * t506 + t461 * t507) * t562 - (t456 * t562 - t458 * t506 + t460 * t507) * t561 + (t494 * t562 - t495 * t506 + t496 * t507) * t528) * t562 + t528 * (t528 ^ 2 * t494 + (((t459 * t533 + t461 * t531) * t524 - (t458 * t533 + t460 * t531) * t527) * t525 + (-t456 * t527 + t457 * t524 + t495 * t533 + t496 * t531) * t528) * t525)) * t571 / 0.2e1;
T  = t1;
