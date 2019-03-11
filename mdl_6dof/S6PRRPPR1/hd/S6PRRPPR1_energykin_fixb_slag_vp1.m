% Calculate kinetic energy for
% S6PRRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4,theta5]';
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
% Datum: 2019-03-08 21:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPPR1_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPPR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:58:30
% EndTime: 2019-03-08 20:58:32
% DurationCPUTime: 2.73s
% Computational Cost: add. (3239->370), mult. (5641->540), div. (0->0), fcn. (6822->14), ass. (0->168)
t583 = Icges(5,2) + Icges(6,3);
t582 = Icges(4,3) + Icges(5,3);
t525 = sin(pkin(10));
t528 = cos(pkin(10));
t535 = cos(qJ(2));
t529 = cos(pkin(6));
t533 = sin(qJ(2));
t559 = t529 * t533;
t506 = t525 * t535 + t528 * t559;
t551 = qJ(3) + pkin(11);
t521 = sin(t551);
t545 = cos(t551);
t526 = sin(pkin(6));
t565 = t526 * t528;
t485 = t506 * t545 - t521 * t565;
t558 = t529 * t535;
t505 = t525 * t533 - t528 * t558;
t524 = sin(pkin(12));
t527 = cos(pkin(12));
t450 = -t485 * t524 + t505 * t527;
t568 = t505 * t524;
t451 = t485 * t527 + t568;
t544 = t526 * t545;
t484 = t506 * t521 + t528 * t544;
t581 = -Icges(5,4) * t485 + Icges(6,5) * t451 - Icges(5,6) * t505 + Icges(6,6) * t450 + t583 * t484;
t508 = -t525 * t559 + t528 * t535;
t566 = t525 * t526;
t487 = t508 * t545 + t521 * t566;
t507 = t525 * t558 + t528 * t533;
t452 = -t487 * t524 + t507 * t527;
t567 = t507 * t524;
t453 = t487 * t527 + t567;
t486 = t508 * t521 - t525 * t544;
t580 = -Icges(5,4) * t487 + Icges(6,5) * t453 - Icges(5,6) * t507 + Icges(6,6) * t452 + t583 * t486;
t499 = t529 * t521 + t533 * t544;
t561 = t526 * t535;
t482 = -t499 * t524 - t527 * t561;
t550 = t524 * t561;
t483 = t499 * t527 - t550;
t563 = t526 * t533;
t498 = t521 * t563 - t529 * t545;
t579 = -Icges(5,4) * t499 + Icges(6,5) * t483 + Icges(5,6) * t561 + Icges(6,6) * t482 + t583 * t498;
t532 = sin(qJ(3));
t534 = cos(qJ(3));
t562 = t526 * t534;
t488 = -t506 * t532 - t528 * t562;
t564 = t526 * t532;
t548 = t528 * t564;
t489 = t506 * t534 - t548;
t578 = Icges(4,5) * t489 + Icges(5,5) * t485 + Icges(4,6) * t488 - Icges(5,6) * t484 + t582 * t505;
t490 = -t508 * t532 + t525 * t562;
t549 = t525 * t564;
t491 = t508 * t534 + t549;
t577 = Icges(4,5) * t491 + Icges(5,5) * t487 + Icges(4,6) * t490 - Icges(5,6) * t486 + t582 * t507;
t509 = t529 * t534 - t532 * t563;
t560 = t529 * t532;
t510 = t533 * t562 + t560;
t576 = Icges(4,5) * t510 + Icges(5,5) * t499 + Icges(4,6) * t509 - Icges(5,6) * t498 - t582 * t561;
t575 = qJD(2) ^ 2;
t571 = pkin(3) * t534;
t570 = pkin(5) * t527;
t427 = -pkin(3) * t548 + qJ(4) * t505 + t506 * t571;
t441 = pkin(4) * t485 + qJ(5) * t484;
t556 = -t427 - t441;
t428 = pkin(3) * t549 + qJ(4) * t507 + t508 * t571;
t442 = pkin(4) * t487 + qJ(5) * t486;
t555 = -t428 - t442;
t458 = t499 * pkin(4) + t498 * qJ(5);
t472 = pkin(3) * t560 + (-qJ(4) * t535 + t533 * t571) * t526;
t554 = -t458 - t472;
t553 = qJD(2) * t526;
t516 = t525 * t553;
t492 = qJD(3) * t507 + t516;
t519 = qJD(2) * t529;
t547 = t528 * t553;
t477 = pkin(2) * t506 + pkin(8) * t505;
t478 = pkin(2) * t508 + pkin(8) * t507;
t546 = t477 * t516 + t478 * t547 + qJD(1);
t493 = qJD(3) * t505 - t547;
t512 = -qJD(3) * t561 + t519;
t511 = (pkin(2) * t533 - pkin(8) * t535) * t526;
t543 = t478 * t519 - t511 * t516;
t542 = qJD(4) * t505 + t512 * t428 + t543;
t541 = (-t477 * t529 - t511 * t565) * qJD(2);
t540 = -qJD(4) * t561 + t492 * t427 + t546;
t539 = qJD(5) * t484 + t512 * t442 + t542;
t538 = qJD(5) * t498 + t492 * t441 + t540;
t537 = qJD(4) * t507 + t493 * t472 + t541;
t536 = qJD(5) * t486 + t493 * t458 + t537;
t523 = pkin(12) + qJ(6);
t522 = cos(t523);
t520 = sin(t523);
t500 = t529 * rSges(3,3) + (rSges(3,1) * t533 + rSges(3,2) * t535) * t526;
t497 = Icges(3,5) * t529 + (Icges(3,1) * t533 + Icges(3,4) * t535) * t526;
t496 = Icges(3,6) * t529 + (Icges(3,4) * t533 + Icges(3,2) * t535) * t526;
t495 = Icges(3,3) * t529 + (Icges(3,5) * t533 + Icges(3,6) * t535) * t526;
t481 = qJD(6) * t498 + t512;
t476 = t499 * t522 - t520 * t561;
t475 = -t499 * t520 - t522 * t561;
t473 = t510 * rSges(4,1) + t509 * rSges(4,2) - rSges(4,3) * t561;
t471 = Icges(4,1) * t510 + Icges(4,4) * t509 - Icges(4,5) * t561;
t470 = Icges(4,4) * t510 + Icges(4,2) * t509 - Icges(4,6) * t561;
t466 = rSges(3,1) * t508 - rSges(3,2) * t507 + rSges(3,3) * t566;
t465 = rSges(3,1) * t506 - rSges(3,2) * t505 - rSges(3,3) * t565;
t464 = Icges(3,1) * t508 - Icges(3,4) * t507 + Icges(3,5) * t566;
t463 = Icges(3,1) * t506 - Icges(3,4) * t505 - Icges(3,5) * t565;
t462 = Icges(3,4) * t508 - Icges(3,2) * t507 + Icges(3,6) * t566;
t461 = Icges(3,4) * t506 - Icges(3,2) * t505 - Icges(3,6) * t565;
t460 = Icges(3,5) * t508 - Icges(3,6) * t507 + Icges(3,3) * t566;
t459 = Icges(3,5) * t506 - Icges(3,6) * t505 - Icges(3,3) * t565;
t457 = t499 * rSges(5,1) - t498 * rSges(5,2) - rSges(5,3) * t561;
t456 = Icges(5,1) * t499 - Icges(5,4) * t498 - Icges(5,5) * t561;
t449 = t487 * t522 + t507 * t520;
t448 = -t487 * t520 + t507 * t522;
t447 = t485 * t522 + t505 * t520;
t446 = -t485 * t520 + t505 * t522;
t445 = qJD(6) * t484 + t493;
t444 = qJD(6) * t486 + t492;
t439 = (-t465 * t529 - t500 * t565) * qJD(2);
t438 = (t466 * t529 - t500 * t566) * qJD(2);
t437 = rSges(4,1) * t491 + rSges(4,2) * t490 + rSges(4,3) * t507;
t436 = rSges(4,1) * t489 + rSges(4,2) * t488 + rSges(4,3) * t505;
t435 = Icges(4,1) * t491 + Icges(4,4) * t490 + Icges(4,5) * t507;
t434 = Icges(4,1) * t489 + Icges(4,4) * t488 + Icges(4,5) * t505;
t433 = Icges(4,4) * t491 + Icges(4,2) * t490 + Icges(4,6) * t507;
t432 = Icges(4,4) * t489 + Icges(4,2) * t488 + Icges(4,6) * t505;
t426 = rSges(5,1) * t487 - rSges(5,2) * t486 + rSges(5,3) * t507;
t425 = rSges(5,1) * t485 - rSges(5,2) * t484 + rSges(5,3) * t505;
t424 = Icges(5,1) * t487 - Icges(5,4) * t486 + Icges(5,5) * t507;
t423 = Icges(5,1) * t485 - Icges(5,4) * t484 + Icges(5,5) * t505;
t418 = rSges(6,1) * t483 + rSges(6,2) * t482 + rSges(6,3) * t498;
t417 = Icges(6,1) * t483 + Icges(6,4) * t482 + Icges(6,5) * t498;
t416 = Icges(6,4) * t483 + Icges(6,2) * t482 + Icges(6,6) * t498;
t413 = rSges(7,1) * t476 + rSges(7,2) * t475 + rSges(7,3) * t498;
t412 = Icges(7,1) * t476 + Icges(7,4) * t475 + Icges(7,5) * t498;
t411 = Icges(7,4) * t476 + Icges(7,2) * t475 + Icges(7,6) * t498;
t410 = Icges(7,5) * t476 + Icges(7,6) * t475 + Icges(7,3) * t498;
t408 = -pkin(5) * t550 + pkin(9) * t498 + t499 * t570;
t407 = qJD(1) + (t465 * t525 + t466 * t528) * t553;
t405 = rSges(6,1) * t453 + rSges(6,2) * t452 + rSges(6,3) * t486;
t404 = rSges(6,1) * t451 + rSges(6,2) * t450 + rSges(6,3) * t484;
t403 = Icges(6,1) * t453 + Icges(6,4) * t452 + Icges(6,5) * t486;
t402 = Icges(6,1) * t451 + Icges(6,4) * t450 + Icges(6,5) * t484;
t401 = Icges(6,4) * t453 + Icges(6,2) * t452 + Icges(6,6) * t486;
t400 = Icges(6,4) * t451 + Icges(6,2) * t450 + Icges(6,6) * t484;
t397 = rSges(7,1) * t449 + rSges(7,2) * t448 + rSges(7,3) * t486;
t396 = rSges(7,1) * t447 + rSges(7,2) * t446 + rSges(7,3) * t484;
t395 = Icges(7,1) * t449 + Icges(7,4) * t448 + Icges(7,5) * t486;
t394 = Icges(7,1) * t447 + Icges(7,4) * t446 + Icges(7,5) * t484;
t393 = Icges(7,4) * t449 + Icges(7,2) * t448 + Icges(7,6) * t486;
t392 = Icges(7,4) * t447 + Icges(7,2) * t446 + Icges(7,6) * t484;
t391 = Icges(7,5) * t449 + Icges(7,6) * t448 + Icges(7,3) * t486;
t390 = Icges(7,5) * t447 + Icges(7,6) * t446 + Icges(7,3) * t484;
t389 = pkin(5) * t567 + pkin(9) * t486 + t487 * t570;
t388 = pkin(5) * t568 + pkin(9) * t484 + t485 * t570;
t387 = -t436 * t512 + t473 * t493 + t541;
t386 = t437 * t512 - t473 * t492 + t543;
t385 = t436 * t492 - t437 * t493 + t546;
t384 = t457 * t493 + (-t425 - t427) * t512 + t537;
t383 = t426 * t512 + (-t457 - t472) * t492 + t542;
t382 = t492 * t425 + (-t426 - t428) * t493 + t540;
t381 = t418 * t493 + (-t404 + t556) * t512 + t536;
t380 = t405 * t512 + (-t418 + t554) * t492 + t539;
t379 = t492 * t404 + (-t405 + t555) * t493 + t538;
t378 = -t396 * t481 + t408 * t493 + t413 * t445 + (-t388 + t556) * t512 + t536;
t377 = t389 * t512 + t397 * t481 - t413 * t444 + (-t408 + t554) * t492 + t539;
t376 = t492 * t388 + t444 * t396 - t445 * t397 + (-t389 + t555) * t493 + t538;
t1 = m(4) * (t385 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + m(3) * (t407 ^ 2 + t438 ^ 2 + t439 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(7) * (t376 ^ 2 + t377 ^ 2 + t378 ^ 2) / 0.2e1 + m(5) * (t382 ^ 2 + t383 ^ 2 + t384 ^ 2) / 0.2e1 - t575 * ((-t460 * t565 - t462 * t505 + t464 * t506) * t566 - (-t459 * t565 - t461 * t505 + t463 * t506) * t565 + (-t495 * t565 - t496 * t505 + t497 * t506) * t529) * t565 / 0.2e1 + m(6) * (t379 ^ 2 + t380 ^ 2 + t381 ^ 2) / 0.2e1 + t444 * ((t486 * t391 + t448 * t393 + t449 * t395) * t444 + (t390 * t486 + t392 * t448 + t394 * t449) * t445 + (t410 * t486 + t411 * t448 + t412 * t449) * t481) / 0.2e1 + t445 * ((t391 * t484 + t393 * t446 + t395 * t447) * t444 + (t484 * t390 + t446 * t392 + t447 * t394) * t445 + (t410 * t484 + t411 * t446 + t412 * t447) * t481) / 0.2e1 + t481 * ((t391 * t498 + t393 * t475 + t395 * t476) * t444 + (t390 * t498 + t392 * t475 + t394 * t476) * t445 + (t498 * t410 + t475 * t411 + t476 * t412) * t481) / 0.2e1 + (t529 * (t529 ^ 2 * t495 + (((t462 * t535 + t464 * t533) * t525 - (t461 * t535 + t463 * t533) * t528) * t526 + (-t459 * t528 + t460 * t525 + t496 * t535 + t497 * t533) * t529) * t526) + ((t460 * t566 - t462 * t507 + t464 * t508) * t566 - (t459 * t566 - t461 * t507 + t463 * t508) * t565 + (t495 * t566 - t496 * t507 + t497 * t508) * t529) * t566) * t575 / 0.2e1 + ((t416 * t452 + t417 * t453 + t456 * t487 + t470 * t490 + t471 * t491 + t486 * t579 + t507 * t576) * t512 + (t400 * t452 + t402 * t453 + t423 * t487 + t432 * t490 + t434 * t491 + t486 * t581 + t578 * t507) * t493 + (t401 * t452 + t403 * t453 + t424 * t487 + t433 * t490 + t435 * t491 + t580 * t486 + t577 * t507) * t492) * t492 / 0.2e1 + ((t416 * t450 + t417 * t451 + t456 * t485 + t470 * t488 + t471 * t489 + t484 * t579 + t505 * t576) * t512 + (t400 * t450 + t402 * t451 + t423 * t485 + t432 * t488 + t434 * t489 + t581 * t484 + t578 * t505) * t493 + (t401 * t450 + t403 * t451 + t424 * t485 + t433 * t488 + t435 * t489 + t484 * t580 + t505 * t577) * t492) * t493 / 0.2e1 + ((t416 * t482 + t417 * t483 + t499 * t456 + t509 * t470 + t510 * t471 + t579 * t498 - t576 * t561) * t512 + (t400 * t482 + t402 * t483 + t499 * t423 + t509 * t432 + t510 * t434 + t498 * t581 - t578 * t561) * t493 + (t401 * t482 + t403 * t483 + t499 * t424 + t509 * t433 + t510 * t435 + t498 * t580 - t561 * t577) * t492) * t512 / 0.2e1;
T  = t1;
