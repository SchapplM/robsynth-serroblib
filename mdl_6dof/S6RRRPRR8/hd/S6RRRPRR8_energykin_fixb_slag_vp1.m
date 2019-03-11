% Calculate kinetic energy for
% S6RRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:55
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR8_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:44:05
% EndTime: 2019-03-09 18:44:08
% DurationCPUTime: 3.17s
% Computational Cost: add. (3476->383), mult. (5909->578), div. (0->0), fcn. (7132->14), ass. (0->173)
t580 = Icges(4,3) + Icges(5,3);
t530 = sin(qJ(2));
t531 = sin(qJ(1));
t534 = cos(qJ(2));
t535 = cos(qJ(1));
t567 = cos(pkin(6));
t550 = t535 * t567;
t504 = t530 * t531 - t534 * t550;
t505 = t530 * t550 + t531 * t534;
t526 = sin(pkin(6));
t560 = t526 * t535;
t461 = Icges(3,5) * t505 - Icges(3,6) * t504 - Icges(3,3) * t560;
t551 = t531 * t567;
t506 = t535 * t530 + t534 * t551;
t507 = -t530 * t551 + t535 * t534;
t563 = t526 * t531;
t462 = Icges(3,5) * t507 - Icges(3,6) * t506 + Icges(3,3) * t563;
t579 = t526 * (t461 * t535 - t462 * t531);
t556 = qJ(3) + pkin(12);
t522 = sin(t556);
t548 = cos(t556);
t547 = t526 * t548;
t481 = t505 * t522 + t535 * t547;
t482 = t505 * t548 - t522 * t560;
t529 = sin(qJ(3));
t533 = cos(qJ(3));
t485 = -t505 * t529 - t533 * t560;
t553 = t529 * t560;
t486 = t505 * t533 - t553;
t578 = Icges(4,5) * t486 + Icges(5,5) * t482 + Icges(4,6) * t485 - Icges(5,6) * t481 + t580 * t504;
t483 = t507 * t522 - t531 * t547;
t484 = t507 * t548 + t522 * t563;
t562 = t526 * t533;
t487 = -t507 * t529 + t531 * t562;
t554 = t529 * t563;
t488 = t507 * t533 + t554;
t577 = Icges(4,5) * t488 + Icges(5,5) * t484 + Icges(4,6) * t487 - Icges(5,6) * t483 + t580 * t506;
t564 = t526 * t530;
t495 = t522 * t564 - t548 * t567;
t496 = t522 * t567 + t530 * t547;
t502 = -t529 * t564 + t533 * t567;
t549 = t567 * t529;
t503 = t530 * t562 + t549;
t561 = t526 * t534;
t576 = Icges(4,5) * t503 + Icges(5,5) * t496 + Icges(4,6) * t502 - Icges(5,6) * t495 - t580 * t561;
t571 = pkin(3) * t533;
t532 = cos(qJ(5));
t570 = pkin(5) * t532;
t528 = sin(qJ(5));
t566 = t504 * t528;
t565 = t506 * t528;
t477 = pkin(2) * t505 + pkin(9) * t504;
t478 = pkin(2) * t507 + pkin(9) * t506;
t557 = qJD(2) * t526;
t517 = t531 * t557;
t552 = t535 * t557;
t559 = t477 * t517 + t478 * t552;
t489 = qJD(3) * t506 + t517;
t558 = qJD(1) * (pkin(1) * t531 - pkin(8) * t560);
t518 = qJD(2) * t567 + qJD(1);
t555 = t528 * t561;
t440 = qJD(5) * t483 + t489;
t490 = qJD(3) * t504 - t552;
t441 = qJD(5) * t481 + t490;
t509 = -qJD(3) * t561 + t518;
t508 = (pkin(2) * t530 - pkin(9) * t534) * t526;
t510 = qJD(1) * (pkin(1) * t535 + pkin(8) * t563);
t545 = t518 * t478 - t508 * t517 + t510;
t425 = -pkin(3) * t553 + qJ(4) * t504 + t505 * t571;
t544 = -qJD(4) * t561 + t489 * t425 + t559;
t472 = qJD(5) * t495 + t509;
t426 = pkin(3) * t554 + qJ(4) * t506 + t507 * t571;
t543 = qJD(4) * t504 + t509 * t426 + t545;
t542 = -t477 * t518 - t508 * t552 - t558;
t438 = pkin(4) * t482 + pkin(10) * t481;
t439 = pkin(4) * t484 + pkin(10) * t483;
t541 = t489 * t438 + (-t426 - t439) * t490 + t544;
t459 = pkin(3) * t549 + (-qJ(4) * t534 + t530 * t571) * t526;
t540 = qJD(4) * t506 + t490 * t459 + t542;
t455 = pkin(4) * t496 + pkin(10) * t495;
t539 = t509 * t439 + (-t455 - t459) * t489 + t543;
t538 = t490 * t455 + (-t425 - t438) * t509 + t540;
t525 = qJ(5) + qJ(6);
t524 = cos(t525);
t523 = sin(t525);
t514 = rSges(2,1) * t535 - rSges(2,2) * t531;
t513 = rSges(2,1) * t531 + rSges(2,2) * t535;
t497 = t567 * rSges(3,3) + (rSges(3,1) * t530 + rSges(3,2) * t534) * t526;
t494 = Icges(3,5) * t567 + (Icges(3,1) * t530 + Icges(3,4) * t534) * t526;
t493 = Icges(3,6) * t567 + (Icges(3,4) * t530 + Icges(3,2) * t534) * t526;
t492 = Icges(3,3) * t567 + (Icges(3,5) * t530 + Icges(3,6) * t534) * t526;
t480 = t496 * t532 - t555;
t479 = -t496 * t528 - t532 * t561;
t474 = t496 * t524 - t523 * t561;
t473 = -t496 * t523 - t524 * t561;
t469 = rSges(3,1) * t507 - rSges(3,2) * t506 + rSges(3,3) * t563;
t468 = rSges(3,1) * t505 - rSges(3,2) * t504 - rSges(3,3) * t560;
t466 = Icges(3,1) * t507 - Icges(3,4) * t506 + Icges(3,5) * t563;
t465 = Icges(3,1) * t505 - Icges(3,4) * t504 - Icges(3,5) * t560;
t464 = Icges(3,4) * t507 - Icges(3,2) * t506 + Icges(3,6) * t563;
t463 = Icges(3,4) * t505 - Icges(3,2) * t504 - Icges(3,6) * t560;
t460 = rSges(4,1) * t503 + rSges(4,2) * t502 - rSges(4,3) * t561;
t458 = Icges(4,1) * t503 + Icges(4,4) * t502 - Icges(4,5) * t561;
t457 = Icges(4,4) * t503 + Icges(4,2) * t502 - Icges(4,6) * t561;
t454 = rSges(5,1) * t496 - rSges(5,2) * t495 - rSges(5,3) * t561;
t453 = Icges(5,1) * t496 - Icges(5,4) * t495 - Icges(5,5) * t561;
t452 = Icges(5,4) * t496 - Icges(5,2) * t495 - Icges(5,6) * t561;
t450 = t484 * t532 + t565;
t449 = -t484 * t528 + t506 * t532;
t448 = t482 * t532 + t566;
t447 = -t482 * t528 + t504 * t532;
t446 = t484 * t524 + t506 * t523;
t445 = -t484 * t523 + t506 * t524;
t444 = t482 * t524 + t504 * t523;
t443 = -t482 * t523 + t504 * t524;
t442 = qJD(6) * t495 + t472;
t435 = rSges(4,1) * t488 + rSges(4,2) * t487 + rSges(4,3) * t506;
t434 = rSges(4,1) * t486 + rSges(4,2) * t485 + rSges(4,3) * t504;
t433 = Icges(4,1) * t488 + Icges(4,4) * t487 + Icges(4,5) * t506;
t432 = Icges(4,1) * t486 + Icges(4,4) * t485 + Icges(4,5) * t504;
t431 = Icges(4,4) * t488 + Icges(4,2) * t487 + Icges(4,6) * t506;
t430 = Icges(4,4) * t486 + Icges(4,2) * t485 + Icges(4,6) * t504;
t424 = rSges(5,1) * t484 - rSges(5,2) * t483 + rSges(5,3) * t506;
t423 = rSges(5,1) * t482 - rSges(5,2) * t481 + rSges(5,3) * t504;
t422 = Icges(5,1) * t484 - Icges(5,4) * t483 + Icges(5,5) * t506;
t421 = Icges(5,1) * t482 - Icges(5,4) * t481 + Icges(5,5) * t504;
t420 = Icges(5,4) * t484 - Icges(5,2) * t483 + Icges(5,6) * t506;
t419 = Icges(5,4) * t482 - Icges(5,2) * t481 + Icges(5,6) * t504;
t416 = t469 * t518 - t497 * t517 + t510;
t415 = -t468 * t518 - t497 * t552 - t558;
t414 = rSges(6,1) * t480 + rSges(6,2) * t479 + rSges(6,3) * t495;
t413 = Icges(6,1) * t480 + Icges(6,4) * t479 + Icges(6,5) * t495;
t412 = Icges(6,4) * t480 + Icges(6,2) * t479 + Icges(6,6) * t495;
t411 = Icges(6,5) * t480 + Icges(6,6) * t479 + Icges(6,3) * t495;
t410 = qJD(6) * t481 + t441;
t409 = qJD(6) * t483 + t440;
t407 = rSges(7,1) * t474 + rSges(7,2) * t473 + rSges(7,3) * t495;
t405 = (t468 * t531 + t469 * t535) * t557;
t404 = Icges(7,1) * t474 + Icges(7,4) * t473 + Icges(7,5) * t495;
t403 = Icges(7,4) * t474 + Icges(7,2) * t473 + Icges(7,6) * t495;
t402 = Icges(7,5) * t474 + Icges(7,6) * t473 + Icges(7,3) * t495;
t401 = -pkin(5) * t555 + pkin(11) * t495 + t496 * t570;
t399 = rSges(6,1) * t450 + rSges(6,2) * t449 + rSges(6,3) * t483;
t398 = rSges(6,1) * t448 + rSges(6,2) * t447 + rSges(6,3) * t481;
t397 = Icges(6,1) * t450 + Icges(6,4) * t449 + Icges(6,5) * t483;
t396 = Icges(6,1) * t448 + Icges(6,4) * t447 + Icges(6,5) * t481;
t395 = Icges(6,4) * t450 + Icges(6,2) * t449 + Icges(6,6) * t483;
t394 = Icges(6,4) * t448 + Icges(6,2) * t447 + Icges(6,6) * t481;
t393 = Icges(6,5) * t450 + Icges(6,6) * t449 + Icges(6,3) * t483;
t392 = Icges(6,5) * t448 + Icges(6,6) * t447 + Icges(6,3) * t481;
t391 = rSges(7,1) * t446 + rSges(7,2) * t445 + rSges(7,3) * t483;
t390 = rSges(7,1) * t444 + rSges(7,2) * t443 + rSges(7,3) * t481;
t389 = Icges(7,1) * t446 + Icges(7,4) * t445 + Icges(7,5) * t483;
t388 = Icges(7,1) * t444 + Icges(7,4) * t443 + Icges(7,5) * t481;
t387 = Icges(7,4) * t446 + Icges(7,2) * t445 + Icges(7,6) * t483;
t386 = Icges(7,4) * t444 + Icges(7,2) * t443 + Icges(7,6) * t481;
t385 = Icges(7,5) * t446 + Icges(7,6) * t445 + Icges(7,3) * t483;
t384 = Icges(7,5) * t444 + Icges(7,6) * t443 + Icges(7,3) * t481;
t383 = pkin(5) * t565 + pkin(11) * t483 + t484 * t570;
t382 = pkin(5) * t566 + pkin(11) * t481 + t482 * t570;
t381 = t435 * t509 - t460 * t489 + t545;
t380 = -t434 * t509 + t460 * t490 + t542;
t379 = t434 * t489 - t435 * t490 + t559;
t378 = t424 * t509 + (-t454 - t459) * t489 + t543;
t377 = t454 * t490 + (-t423 - t425) * t509 + t540;
t376 = t423 * t489 + (-t424 - t426) * t490 + t544;
t375 = t399 * t472 - t414 * t440 + t539;
t374 = -t398 * t472 + t414 * t441 + t538;
t373 = t398 * t440 - t399 * t441 + t541;
t372 = t383 * t472 + t391 * t442 - t401 * t440 - t407 * t409 + t539;
t371 = -t382 * t472 - t390 * t442 + t401 * t441 + t407 * t410 + t538;
t370 = t382 * t440 - t383 * t441 + t390 * t409 - t391 * t410 + t541;
t1 = t410 * ((t385 * t481 + t387 * t443 + t389 * t444) * t409 + (t481 * t384 + t443 * t386 + t444 * t388) * t410 + (t402 * t481 + t403 * t443 + t404 * t444) * t442) / 0.2e1 + t409 * ((t483 * t385 + t445 * t387 + t446 * t389) * t409 + (t384 * t483 + t386 * t445 + t388 * t446) * t410 + (t402 * t483 + t403 * t445 + t404 * t446) * t442) / 0.2e1 + t442 * ((t385 * t495 + t387 * t473 + t389 * t474) * t409 + (t384 * t495 + t386 * t473 + t388 * t474) * t410 + (t495 * t402 + t473 * t403 + t474 * t404) * t442) / 0.2e1 + t441 * ((t393 * t481 + t395 * t447 + t397 * t448) * t440 + (t481 * t392 + t447 * t394 + t448 * t396) * t441 + (t411 * t481 + t412 * t447 + t413 * t448) * t472) / 0.2e1 + t440 * ((t483 * t393 + t449 * t395 + t450 * t397) * t440 + (t392 * t483 + t394 * t449 + t396 * t450) * t441 + (t411 * t483 + t412 * t449 + t413 * t450) * t472) / 0.2e1 + t472 * ((t393 * t495 + t395 * t479 + t397 * t480) * t440 + (t392 * t495 + t394 * t479 + t396 * t480) * t441 + (t495 * t411 + t479 * t412 + t480 * t413) * t472) / 0.2e1 + t518 * ((t567 * t462 + (t464 * t534 + t466 * t530) * t526) * t517 - (t567 * t461 + (t463 * t534 + t465 * t530) * t526) * t552 + (t567 * t492 + (t493 * t534 + t494 * t530) * t526) * t518) / 0.2e1 + m(7) * (t370 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + m(6) * (t373 ^ 2 + t374 ^ 2 + t375 ^ 2) / 0.2e1 + m(5) * (t376 ^ 2 + t377 ^ 2 + t378 ^ 2) / 0.2e1 + m(4) * (t379 ^ 2 + t380 ^ 2 + t381 ^ 2) / 0.2e1 + m(3) * (t405 ^ 2 + t415 ^ 2 + t416 ^ 2) / 0.2e1 + ((t492 * t563 - t493 * t506 + t494 * t507) * t518 + (-(-t463 * t506 + t465 * t507) * t535 + (-t506 * t464 + t507 * t466 - t579) * t531) * t557) * t517 / 0.2e1 - ((-t492 * t560 - t493 * t504 + t494 * t505) * t518 + ((-t464 * t504 + t466 * t505) * t531 + (t504 * t463 - t505 * t465 + t579) * t535) * t557) * t552 / 0.2e1 + ((-t452 * t483 + t453 * t484 + t457 * t487 + t458 * t488 + t506 * t576) * t509 + (-t419 * t483 + t421 * t484 + t430 * t487 + t432 * t488 + t506 * t578) * t490 + (-t483 * t420 + t484 * t422 + t487 * t431 + t488 * t433 + t577 * t506) * t489) * t489 / 0.2e1 + ((-t452 * t481 + t453 * t482 + t457 * t485 + t458 * t486 + t504 * t576) * t509 + (-t481 * t419 + t482 * t421 + t485 * t430 + t486 * t432 + t578 * t504) * t490 + (-t420 * t481 + t422 * t482 + t431 * t485 + t433 * t486 + t504 * t577) * t489) * t490 / 0.2e1 + ((-t495 * t452 + t496 * t453 + t502 * t457 + t503 * t458 - t576 * t561) * t509 + (-t419 * t495 + t421 * t496 + t430 * t502 + t432 * t503 - t561 * t578) * t490 + (-t420 * t495 + t422 * t496 + t431 * t502 + t433 * t503 - t561 * t577) * t489) * t509 / 0.2e1 + (m(2) * (t513 ^ 2 + t514 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
