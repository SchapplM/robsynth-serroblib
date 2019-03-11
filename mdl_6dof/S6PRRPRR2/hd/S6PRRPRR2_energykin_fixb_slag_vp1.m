% Calculate kinetic energy for
% S6PRRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR2_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:55:53
% EndTime: 2019-03-08 21:55:56
% DurationCPUTime: 3.03s
% Computational Cost: add. (3407->375), mult. (5857->570), div. (0->0), fcn. (7098->14), ass. (0->170)
t575 = Icges(4,3) + Icges(5,3);
t523 = sin(pkin(11));
t525 = cos(pkin(11));
t533 = cos(qJ(2));
t526 = cos(pkin(6));
t530 = sin(qJ(2));
t554 = t526 * t530;
t505 = t523 * t533 + t525 * t554;
t550 = qJ(3) + pkin(12);
t519 = sin(t550);
t524 = sin(pkin(6));
t544 = cos(t550);
t543 = t524 * t544;
t481 = t505 * t519 + t525 * t543;
t560 = t524 * t525;
t482 = t505 * t544 - t519 * t560;
t529 = sin(qJ(3));
t532 = cos(qJ(3));
t557 = t524 * t532;
t487 = -t505 * t529 - t525 * t557;
t559 = t524 * t529;
t548 = t525 * t559;
t488 = t505 * t532 - t548;
t553 = t526 * t533;
t504 = t523 * t530 - t525 * t553;
t574 = Icges(4,5) * t488 + Icges(5,5) * t482 + Icges(4,6) * t487 - Icges(5,6) * t481 + t504 * t575;
t507 = -t523 * t554 + t525 * t533;
t483 = t507 * t519 - t523 * t543;
t561 = t523 * t524;
t484 = t507 * t544 + t519 * t561;
t489 = -t507 * t529 + t523 * t557;
t549 = t523 * t559;
t490 = t507 * t532 + t549;
t506 = t523 * t553 + t525 * t530;
t573 = Icges(4,5) * t490 + Icges(5,5) * t484 + Icges(4,6) * t489 - Icges(5,6) * t483 + t506 * t575;
t558 = t524 * t530;
t497 = t519 * t558 - t526 * t544;
t498 = t526 * t519 + t530 * t543;
t508 = t526 * t532 - t529 * t558;
t555 = t526 * t529;
t509 = t530 * t557 + t555;
t556 = t524 * t533;
t572 = Icges(4,5) * t509 + Icges(5,5) * t498 + Icges(4,6) * t508 - Icges(5,6) * t497 - t556 * t575;
t571 = qJD(2) ^ 2;
t567 = pkin(3) * t532;
t531 = cos(qJ(5));
t566 = pkin(5) * t531;
t528 = sin(qJ(5));
t563 = t504 * t528;
t562 = t506 * t528;
t552 = qJD(2) * t524;
t515 = t523 * t552;
t491 = qJD(3) * t506 + t515;
t518 = qJD(2) * t526;
t547 = t528 * t556;
t442 = qJD(5) * t483 + t491;
t546 = t525 * t552;
t476 = pkin(2) * t505 + pkin(8) * t504;
t477 = pkin(2) * t507 + pkin(8) * t506;
t545 = t476 * t515 + t477 * t546 + qJD(1);
t492 = qJD(3) * t504 - t546;
t511 = -qJD(3) * t556 + t518;
t510 = (pkin(2) * t530 - pkin(8) * t533) * t524;
t542 = t477 * t518 - t510 * t515;
t443 = qJD(5) * t481 + t492;
t480 = qJD(5) * t497 + t511;
t426 = pkin(3) * t549 + qJ(4) * t506 + t507 * t567;
t541 = qJD(4) * t504 + t511 * t426 + t542;
t540 = (-t476 * t526 - t510 * t560) * qJD(2);
t425 = -pkin(3) * t548 + qJ(4) * t504 + t505 * t567;
t539 = -qJD(4) * t556 + t491 * t425 + t545;
t471 = pkin(3) * t555 + (-qJ(4) * t533 + t530 * t567) * t524;
t538 = qJD(4) * t506 + t492 * t471 + t540;
t440 = t484 * pkin(4) + t483 * pkin(9);
t463 = t498 * pkin(4) + t497 * pkin(9);
t537 = t511 * t440 + (-t463 - t471) * t491 + t541;
t439 = t482 * pkin(4) + t481 * pkin(9);
t536 = t491 * t439 + (-t426 - t440) * t492 + t539;
t535 = t492 * t463 + (-t425 - t439) * t511 + t538;
t522 = qJ(5) + qJ(6);
t521 = cos(t522);
t520 = sin(t522);
t499 = rSges(3,3) * t526 + (rSges(3,1) * t530 + rSges(3,2) * t533) * t524;
t496 = Icges(3,5) * t526 + (Icges(3,1) * t530 + Icges(3,4) * t533) * t524;
t495 = Icges(3,6) * t526 + (Icges(3,4) * t530 + Icges(3,2) * t533) * t524;
t494 = Icges(3,3) * t526 + (Icges(3,5) * t530 + Icges(3,6) * t533) * t524;
t486 = t498 * t531 - t547;
t485 = -t498 * t528 - t531 * t556;
t475 = t498 * t521 - t520 * t556;
t474 = -t498 * t520 - t521 * t556;
t472 = rSges(4,1) * t509 + rSges(4,2) * t508 - rSges(4,3) * t556;
t470 = Icges(4,1) * t509 + Icges(4,4) * t508 - Icges(4,5) * t556;
t469 = Icges(4,4) * t509 + Icges(4,2) * t508 - Icges(4,6) * t556;
t465 = rSges(3,1) * t507 - rSges(3,2) * t506 + rSges(3,3) * t561;
t464 = rSges(3,1) * t505 - rSges(3,2) * t504 - rSges(3,3) * t560;
t462 = Icges(3,1) * t507 - Icges(3,4) * t506 + Icges(3,5) * t561;
t461 = Icges(3,1) * t505 - Icges(3,4) * t504 - Icges(3,5) * t560;
t460 = Icges(3,4) * t507 - Icges(3,2) * t506 + Icges(3,6) * t561;
t459 = Icges(3,4) * t505 - Icges(3,2) * t504 - Icges(3,6) * t560;
t458 = Icges(3,5) * t507 - Icges(3,6) * t506 + Icges(3,3) * t561;
t457 = Icges(3,5) * t505 - Icges(3,6) * t504 - Icges(3,3) * t560;
t456 = rSges(5,1) * t498 - rSges(5,2) * t497 - rSges(5,3) * t556;
t455 = Icges(5,1) * t498 - Icges(5,4) * t497 - Icges(5,5) * t556;
t454 = Icges(5,4) * t498 - Icges(5,2) * t497 - Icges(5,6) * t556;
t452 = t484 * t531 + t562;
t451 = -t484 * t528 + t506 * t531;
t450 = t482 * t531 + t563;
t449 = -t482 * t528 + t504 * t531;
t448 = qJD(6) * t497 + t480;
t447 = t484 * t521 + t506 * t520;
t446 = -t484 * t520 + t506 * t521;
t445 = t482 * t521 + t504 * t520;
t444 = -t482 * t520 + t504 * t521;
t437 = (-t464 * t526 - t499 * t560) * qJD(2);
t436 = (t465 * t526 - t499 * t561) * qJD(2);
t435 = rSges(4,1) * t490 + rSges(4,2) * t489 + rSges(4,3) * t506;
t434 = rSges(4,1) * t488 + rSges(4,2) * t487 + rSges(4,3) * t504;
t433 = Icges(4,1) * t490 + Icges(4,4) * t489 + Icges(4,5) * t506;
t432 = Icges(4,1) * t488 + Icges(4,4) * t487 + Icges(4,5) * t504;
t431 = Icges(4,4) * t490 + Icges(4,2) * t489 + Icges(4,6) * t506;
t430 = Icges(4,4) * t488 + Icges(4,2) * t487 + Icges(4,6) * t504;
t424 = rSges(5,1) * t484 - rSges(5,2) * t483 + rSges(5,3) * t506;
t423 = rSges(5,1) * t482 - rSges(5,2) * t481 + rSges(5,3) * t504;
t422 = Icges(5,1) * t484 - Icges(5,4) * t483 + Icges(5,5) * t506;
t421 = Icges(5,1) * t482 - Icges(5,4) * t481 + Icges(5,5) * t504;
t420 = Icges(5,4) * t484 - Icges(5,2) * t483 + Icges(5,6) * t506;
t419 = Icges(5,4) * t482 - Icges(5,2) * t481 + Icges(5,6) * t504;
t416 = rSges(6,1) * t486 + rSges(6,2) * t485 + rSges(6,3) * t497;
t415 = Icges(6,1) * t486 + Icges(6,4) * t485 + Icges(6,5) * t497;
t414 = Icges(6,4) * t486 + Icges(6,2) * t485 + Icges(6,6) * t497;
t413 = Icges(6,5) * t486 + Icges(6,6) * t485 + Icges(6,3) * t497;
t412 = qJD(6) * t481 + t443;
t411 = qJD(6) * t483 + t442;
t409 = rSges(7,1) * t475 + rSges(7,2) * t474 + rSges(7,3) * t497;
t408 = Icges(7,1) * t475 + Icges(7,4) * t474 + Icges(7,5) * t497;
t407 = Icges(7,4) * t475 + Icges(7,2) * t474 + Icges(7,6) * t497;
t406 = Icges(7,5) * t475 + Icges(7,6) * t474 + Icges(7,3) * t497;
t404 = -pkin(5) * t547 + pkin(10) * t497 + t498 * t566;
t403 = qJD(1) + (t464 * t523 + t465 * t525) * t552;
t401 = rSges(6,1) * t452 + rSges(6,2) * t451 + rSges(6,3) * t483;
t400 = rSges(6,1) * t450 + rSges(6,2) * t449 + rSges(6,3) * t481;
t399 = Icges(6,1) * t452 + Icges(6,4) * t451 + Icges(6,5) * t483;
t398 = Icges(6,1) * t450 + Icges(6,4) * t449 + Icges(6,5) * t481;
t397 = Icges(6,4) * t452 + Icges(6,2) * t451 + Icges(6,6) * t483;
t396 = Icges(6,4) * t450 + Icges(6,2) * t449 + Icges(6,6) * t481;
t395 = Icges(6,5) * t452 + Icges(6,6) * t451 + Icges(6,3) * t483;
t394 = Icges(6,5) * t450 + Icges(6,6) * t449 + Icges(6,3) * t481;
t393 = rSges(7,1) * t447 + rSges(7,2) * t446 + rSges(7,3) * t483;
t392 = rSges(7,1) * t445 + rSges(7,2) * t444 + rSges(7,3) * t481;
t391 = Icges(7,1) * t447 + Icges(7,4) * t446 + Icges(7,5) * t483;
t390 = Icges(7,1) * t445 + Icges(7,4) * t444 + Icges(7,5) * t481;
t389 = Icges(7,4) * t447 + Icges(7,2) * t446 + Icges(7,6) * t483;
t388 = Icges(7,4) * t445 + Icges(7,2) * t444 + Icges(7,6) * t481;
t387 = Icges(7,5) * t447 + Icges(7,6) * t446 + Icges(7,3) * t483;
t386 = Icges(7,5) * t445 + Icges(7,6) * t444 + Icges(7,3) * t481;
t385 = pkin(5) * t562 + pkin(10) * t483 + t484 * t566;
t384 = pkin(5) * t563 + pkin(10) * t481 + t482 * t566;
t383 = -t434 * t511 + t472 * t492 + t540;
t382 = t435 * t511 - t472 * t491 + t542;
t381 = t434 * t491 - t435 * t492 + t545;
t380 = t456 * t492 + (-t423 - t425) * t511 + t538;
t379 = t424 * t511 + (-t456 - t471) * t491 + t541;
t378 = t423 * t491 + (-t424 - t426) * t492 + t539;
t377 = -t400 * t480 + t416 * t443 + t535;
t376 = t401 * t480 - t416 * t442 + t537;
t375 = t400 * t442 - t401 * t443 + t536;
t374 = -t384 * t480 - t392 * t448 + t404 * t443 + t409 * t412 + t535;
t373 = t385 * t480 + t393 * t448 - t404 * t442 - t409 * t411 + t537;
t372 = t384 * t442 - t385 * t443 + t392 * t411 - t393 * t412 + t536;
t1 = -t571 * ((-t458 * t560 - t460 * t504 + t462 * t505) * t561 - (-t457 * t560 - t459 * t504 + t461 * t505) * t560 + (-t494 * t560 - t495 * t504 + t496 * t505) * t526) * t560 / 0.2e1 + t411 * ((t483 * t387 + t446 * t389 + t447 * t391) * t411 + (t386 * t483 + t388 * t446 + t390 * t447) * t412 + (t406 * t483 + t407 * t446 + t408 * t447) * t448) / 0.2e1 + t412 * ((t387 * t481 + t389 * t444 + t391 * t445) * t411 + (t481 * t386 + t444 * t388 + t445 * t390) * t412 + (t406 * t481 + t407 * t444 + t408 * t445) * t448) / 0.2e1 + t448 * ((t387 * t497 + t389 * t474 + t391 * t475) * t411 + (t386 * t497 + t388 * t474 + t390 * t475) * t412 + (t497 * t406 + t474 * t407 + t475 * t408) * t448) / 0.2e1 + t480 * ((t395 * t497 + t397 * t485 + t399 * t486) * t442 + (t394 * t497 + t396 * t485 + t398 * t486) * t443 + (t497 * t413 + t485 * t414 + t486 * t415) * t480) / 0.2e1 + t442 * ((t483 * t395 + t451 * t397 + t452 * t399) * t442 + (t394 * t483 + t396 * t451 + t398 * t452) * t443 + (t413 * t483 + t414 * t451 + t415 * t452) * t480) / 0.2e1 + t443 * ((t395 * t481 + t397 * t449 + t399 * t450) * t442 + (t481 * t394 + t449 * t396 + t450 * t398) * t443 + (t413 * t481 + t414 * t449 + t415 * t450) * t480) / 0.2e1 + m(7) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(6) * (t375 ^ 2 + t376 ^ 2 + t377 ^ 2) / 0.2e1 + m(5) * (t378 ^ 2 + t379 ^ 2 + t380 ^ 2) / 0.2e1 + m(4) * (t381 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + m(3) * (t403 ^ 2 + t436 ^ 2 + t437 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + ((-t454 * t483 + t455 * t484 + t469 * t489 + t470 * t490 + t506 * t572) * t511 + (-t419 * t483 + t421 * t484 + t430 * t489 + t432 * t490 + t506 * t574) * t492 + (-t420 * t483 + t422 * t484 + t431 * t489 + t433 * t490 + t506 * t573) * t491) * t491 / 0.2e1 + ((-t454 * t481 + t455 * t482 + t469 * t487 + t470 * t488 + t504 * t572) * t511 + (-t419 * t481 + t421 * t482 + t430 * t487 + t432 * t488 + t504 * t574) * t492 + (-t420 * t481 + t422 * t482 + t431 * t487 + t433 * t488 + t504 * t573) * t491) * t492 / 0.2e1 + ((-t454 * t497 + t455 * t498 + t469 * t508 + t470 * t509 - t556 * t572) * t511 + (-t419 * t497 + t421 * t498 + t430 * t508 + t432 * t509 - t556 * t574) * t492 + (-t420 * t497 + t422 * t498 + t431 * t508 + t433 * t509 - t556 * t573) * t491) * t511 / 0.2e1 + (((t458 * t561 - t460 * t506 + t462 * t507) * t561 - (t457 * t561 - t459 * t506 + t461 * t507) * t560 + (t494 * t561 - t495 * t506 + t496 * t507) * t526) * t561 + t526 * (t526 ^ 2 * t494 + (((t460 * t533 + t462 * t530) * t523 - (t459 * t533 + t461 * t530) * t525) * t524 + (-t457 * t525 + t458 * t523 + t495 * t533 + t496 * t530) * t526) * t524)) * t571 / 0.2e1;
T  = t1;
