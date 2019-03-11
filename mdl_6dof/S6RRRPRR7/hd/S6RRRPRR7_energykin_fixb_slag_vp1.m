% Calculate kinetic energy for
% S6RRRPRR7
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
% Datum: 2019-03-09 18:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:32:46
% EndTime: 2019-03-09 18:32:49
% DurationCPUTime: 3.30s
% Computational Cost: add. (3418->383), mult. (5345->577), div. (0->0), fcn. (6312->14), ass. (0->172)
t586 = Icges(4,3) + Icges(5,3);
t532 = sin(pkin(6));
t533 = cos(pkin(6));
t541 = cos(qJ(2));
t542 = cos(qJ(1));
t565 = t541 * t542;
t537 = sin(qJ(2));
t538 = sin(qJ(1));
t568 = t537 * t538;
t507 = -t533 * t565 + t568;
t566 = t538 * t541;
t567 = t537 * t542;
t508 = t533 * t567 + t566;
t509 = t533 * t566 + t567;
t510 = -t533 * t568 + t565;
t570 = t532 * t542;
t573 = t532 * t538;
t553 = (Icges(3,5) * t508 - Icges(3,6) * t507 - Icges(3,3) * t570) * t542 - (Icges(3,5) * t510 - Icges(3,6) * t509 + Icges(3,3) * t573) * t538;
t585 = t532 * t553;
t531 = qJ(3) + pkin(12);
t527 = sin(t531);
t528 = cos(t531);
t480 = -t508 * t527 - t528 * t570;
t481 = t508 * t528 - t527 * t570;
t536 = sin(qJ(3));
t540 = cos(qJ(3));
t484 = -t508 * t536 - t540 * t570;
t558 = t536 * t570;
t485 = t508 * t540 - t558;
t584 = Icges(4,5) * t485 + Icges(5,5) * t481 + Icges(4,6) * t484 + Icges(5,6) * t480 + t507 * t586;
t482 = -t510 * t527 + t528 * t573;
t483 = t510 * t528 + t527 * t573;
t572 = t532 * t540;
t486 = -t510 * t536 + t538 * t572;
t559 = t536 * t573;
t487 = t510 * t540 + t559;
t583 = Icges(4,5) * t487 + Icges(5,5) * t483 + Icges(4,6) * t486 + Icges(5,6) * t482 + t509 * t586;
t574 = t532 * t537;
t496 = -t527 * t574 + t528 * t533;
t497 = t527 * t533 + t528 * t574;
t505 = t533 * t540 - t536 * t574;
t569 = t533 * t536;
t506 = t537 * t572 + t569;
t571 = t532 * t541;
t582 = Icges(4,5) * t506 + Icges(5,5) * t497 + Icges(4,6) * t505 + Icges(5,6) * t496 - t571 * t586;
t576 = t540 * pkin(3);
t478 = pkin(2) * t508 + pkin(9) * t507;
t479 = pkin(2) * t510 + pkin(9) * t509;
t560 = qJD(2) * t532;
t522 = t538 * t560;
t557 = t542 * t560;
t564 = t478 * t522 + t479 * t557;
t563 = pkin(4) * t528;
t488 = qJD(3) * t509 + t522;
t561 = qJD(1) * (pkin(1) * t538 - pkin(8) * t570);
t523 = qJD(2) * t533 + qJD(1);
t457 = qJD(5) * t509 + t488;
t556 = qJ(5) + t531;
t555 = pkin(4) * t527;
t554 = cos(t556);
t489 = qJD(3) * t507 - t557;
t458 = qJD(5) * t507 + t489;
t552 = t532 * t554;
t511 = (pkin(2) * t537 - pkin(9) * t541) * t532;
t513 = qJD(1) * (pkin(1) * t542 + pkin(8) * t573);
t551 = t523 * t479 - t511 * t522 + t513;
t424 = -pkin(3) * t558 + qJ(4) * t507 + t508 * t576;
t550 = -qJD(4) * t571 + t488 * t424 + t564;
t492 = (-qJD(3) - qJD(5)) * t571 + t523;
t425 = pkin(3) * t559 + qJ(4) * t509 + t510 * t576;
t512 = -qJD(3) * t571 + t523;
t549 = qJD(4) * t507 + t512 * t425 + t551;
t548 = -t478 * t523 - t511 * t557 - t561;
t396 = pkin(10) * t507 + t508 * t563 - t555 * t570;
t397 = pkin(10) * t509 + t510 * t563 + t555 * t573;
t547 = t488 * t396 + (-t397 - t425) * t489 + t550;
t459 = pkin(3) * t569 + (-qJ(4) * t541 + t537 * t576) * t532;
t546 = qJD(4) * t509 + t489 * t459 + t548;
t438 = t555 * t533 + (-pkin(10) * t541 + t537 * t563) * t532;
t545 = t512 * t397 + (-t438 - t459) * t488 + t549;
t544 = t489 * t438 + (-t396 - t424) * t512 + t546;
t539 = cos(qJ(6));
t535 = sin(qJ(6));
t524 = sin(t556);
t520 = rSges(2,1) * t542 - rSges(2,2) * t538;
t519 = rSges(2,1) * t538 + rSges(2,2) * t542;
t498 = rSges(3,3) * t533 + (rSges(3,1) * t537 + rSges(3,2) * t541) * t532;
t495 = Icges(3,5) * t533 + (Icges(3,1) * t537 + Icges(3,4) * t541) * t532;
t494 = Icges(3,6) * t533 + (Icges(3,4) * t537 + Icges(3,2) * t541) * t532;
t493 = Icges(3,3) * t533 + (Icges(3,5) * t537 + Icges(3,6) * t541) * t532;
t491 = t533 * t524 + t537 * t552;
t490 = t524 * t574 - t533 * t554;
t477 = t510 * t554 + t524 * t573;
t476 = t510 * t524 - t538 * t552;
t475 = t508 * t554 - t524 * t570;
t474 = t508 * t524 + t542 * t552;
t473 = t491 * t539 - t535 * t571;
t472 = -t491 * t535 - t539 * t571;
t469 = rSges(3,1) * t510 - rSges(3,2) * t509 + rSges(3,3) * t573;
t468 = rSges(3,1) * t508 - rSges(3,2) * t507 - rSges(3,3) * t570;
t466 = Icges(3,1) * t510 - Icges(3,4) * t509 + Icges(3,5) * t573;
t465 = Icges(3,1) * t508 - Icges(3,4) * t507 - Icges(3,5) * t570;
t464 = Icges(3,4) * t510 - Icges(3,2) * t509 + Icges(3,6) * t573;
t463 = Icges(3,4) * t508 - Icges(3,2) * t507 - Icges(3,6) * t570;
t460 = rSges(4,1) * t506 + rSges(4,2) * t505 - rSges(4,3) * t571;
t456 = Icges(4,1) * t506 + Icges(4,4) * t505 - Icges(4,5) * t571;
t455 = Icges(4,4) * t506 + Icges(4,2) * t505 - Icges(4,6) * t571;
t453 = qJD(6) * t490 + t492;
t452 = pkin(5) * t491 + pkin(11) * t490;
t451 = rSges(5,1) * t497 + rSges(5,2) * t496 - rSges(5,3) * t571;
t450 = Icges(5,1) * t497 + Icges(5,4) * t496 - Icges(5,5) * t571;
t449 = Icges(5,4) * t497 + Icges(5,2) * t496 - Icges(5,6) * t571;
t447 = rSges(6,1) * t491 - rSges(6,2) * t490 - rSges(6,3) * t571;
t446 = Icges(6,1) * t491 - Icges(6,4) * t490 - Icges(6,5) * t571;
t445 = Icges(6,4) * t491 - Icges(6,2) * t490 - Icges(6,6) * t571;
t444 = Icges(6,5) * t491 - Icges(6,6) * t490 - Icges(6,3) * t571;
t443 = t477 * t539 + t509 * t535;
t442 = -t477 * t535 + t509 * t539;
t441 = t475 * t539 + t507 * t535;
t440 = -t475 * t535 + t507 * t539;
t437 = pkin(5) * t477 + pkin(11) * t476;
t436 = pkin(5) * t475 + pkin(11) * t474;
t435 = rSges(4,1) * t487 + rSges(4,2) * t486 + rSges(4,3) * t509;
t434 = rSges(4,1) * t485 + rSges(4,2) * t484 + rSges(4,3) * t507;
t433 = Icges(4,1) * t487 + Icges(4,4) * t486 + Icges(4,5) * t509;
t432 = Icges(4,1) * t485 + Icges(4,4) * t484 + Icges(4,5) * t507;
t431 = Icges(4,4) * t487 + Icges(4,2) * t486 + Icges(4,6) * t509;
t430 = Icges(4,4) * t485 + Icges(4,2) * t484 + Icges(4,6) * t507;
t427 = qJD(6) * t474 + t458;
t426 = qJD(6) * t476 + t457;
t423 = rSges(5,1) * t483 + rSges(5,2) * t482 + rSges(5,3) * t509;
t422 = rSges(5,1) * t481 + rSges(5,2) * t480 + rSges(5,3) * t507;
t421 = Icges(5,1) * t483 + Icges(5,4) * t482 + Icges(5,5) * t509;
t420 = Icges(5,1) * t481 + Icges(5,4) * t480 + Icges(5,5) * t507;
t419 = Icges(5,4) * t483 + Icges(5,2) * t482 + Icges(5,6) * t509;
t418 = Icges(5,4) * t481 + Icges(5,2) * t480 + Icges(5,6) * t507;
t415 = rSges(6,1) * t477 - rSges(6,2) * t476 + rSges(6,3) * t509;
t414 = rSges(6,1) * t475 - rSges(6,2) * t474 + rSges(6,3) * t507;
t413 = Icges(6,1) * t477 - Icges(6,4) * t476 + Icges(6,5) * t509;
t412 = Icges(6,1) * t475 - Icges(6,4) * t474 + Icges(6,5) * t507;
t411 = Icges(6,4) * t477 - Icges(6,2) * t476 + Icges(6,6) * t509;
t410 = Icges(6,4) * t475 - Icges(6,2) * t474 + Icges(6,6) * t507;
t409 = Icges(6,5) * t477 - Icges(6,6) * t476 + Icges(6,3) * t509;
t408 = Icges(6,5) * t475 - Icges(6,6) * t474 + Icges(6,3) * t507;
t407 = t469 * t523 - t498 * t522 + t513;
t406 = -t468 * t523 - t498 * t557 - t561;
t404 = rSges(7,1) * t473 + rSges(7,2) * t472 + rSges(7,3) * t490;
t403 = Icges(7,1) * t473 + Icges(7,4) * t472 + Icges(7,5) * t490;
t402 = Icges(7,4) * t473 + Icges(7,2) * t472 + Icges(7,6) * t490;
t401 = Icges(7,5) * t473 + Icges(7,6) * t472 + Icges(7,3) * t490;
t399 = (t468 * t538 + t469 * t542) * t560;
t393 = rSges(7,1) * t443 + rSges(7,2) * t442 + rSges(7,3) * t476;
t392 = rSges(7,1) * t441 + rSges(7,2) * t440 + rSges(7,3) * t474;
t391 = Icges(7,1) * t443 + Icges(7,4) * t442 + Icges(7,5) * t476;
t390 = Icges(7,1) * t441 + Icges(7,4) * t440 + Icges(7,5) * t474;
t389 = Icges(7,4) * t443 + Icges(7,2) * t442 + Icges(7,6) * t476;
t388 = Icges(7,4) * t441 + Icges(7,2) * t440 + Icges(7,6) * t474;
t387 = Icges(7,5) * t443 + Icges(7,6) * t442 + Icges(7,3) * t476;
t386 = Icges(7,5) * t441 + Icges(7,6) * t440 + Icges(7,3) * t474;
t385 = t435 * t512 - t460 * t488 + t551;
t384 = -t434 * t512 + t460 * t489 + t548;
t383 = t434 * t488 - t435 * t489 + t564;
t382 = t423 * t512 + (-t451 - t459) * t488 + t549;
t381 = t451 * t489 + (-t422 - t424) * t512 + t546;
t380 = t422 * t488 + (-t423 - t425) * t489 + t550;
t379 = t415 * t492 - t447 * t457 + t545;
t378 = -t414 * t492 + t447 * t458 + t544;
t377 = t414 * t457 - t415 * t458 + t547;
t376 = t393 * t453 - t404 * t426 + t437 * t492 - t452 * t457 + t545;
t375 = -t392 * t453 + t404 * t427 - t436 * t492 + t452 * t458 + t544;
t374 = t392 * t426 - t393 * t427 + t436 * t457 - t437 * t458 + t547;
t1 = m(7) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + m(6) * (t377 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + m(5) * (t380 ^ 2 + t381 ^ 2 + t382 ^ 2) / 0.2e1 + m(4) * (t383 ^ 2 + t384 ^ 2 + t385 ^ 2) / 0.2e1 + m(3) * (t399 ^ 2 + t406 ^ 2 + t407 ^ 2) / 0.2e1 + t523 * ((t493 * t533 + (t494 * t541 + t495 * t537) * t532) * t523 + (((t464 * t541 + t466 * t537) * t538 - (t463 * t541 + t465 * t537) * t542) * t532 - t553 * t533) * t560) / 0.2e1 + t453 * ((t387 * t490 + t389 * t472 + t391 * t473) * t426 + (t386 * t490 + t388 * t472 + t390 * t473) * t427 + (t401 * t490 + t402 * t472 + t403 * t473) * t453) / 0.2e1 + t426 * ((t387 * t476 + t389 * t442 + t391 * t443) * t426 + (t386 * t476 + t388 * t442 + t390 * t443) * t427 + (t401 * t476 + t402 * t442 + t403 * t443) * t453) / 0.2e1 + t427 * ((t387 * t474 + t389 * t440 + t391 * t441) * t426 + (t386 * t474 + t388 * t440 + t390 * t441) * t427 + (t401 * t474 + t402 * t440 + t403 * t441) * t453) / 0.2e1 + t458 * ((t409 * t507 - t411 * t474 + t413 * t475) * t457 + (t507 * t408 - t474 * t410 + t475 * t412) * t458 + (t444 * t507 - t445 * t474 + t446 * t475) * t492) / 0.2e1 + t492 * ((-t409 * t571 - t411 * t490 + t413 * t491) * t457 + (-t408 * t571 - t410 * t490 + t412 * t491) * t458 + (-t444 * t571 - t445 * t490 + t446 * t491) * t492) / 0.2e1 + t457 * ((t409 * t509 - t411 * t476 + t413 * t477) * t457 + (t408 * t509 - t410 * t476 + t412 * t477) * t458 + (t444 * t509 - t445 * t476 + t446 * t477) * t492) / 0.2e1 - ((-t493 * t570 - t494 * t507 + t495 * t508) * t523 + ((-t464 * t507 + t466 * t508) * t538 + (t463 * t507 - t465 * t508 + t585) * t542) * t560) * t557 / 0.2e1 + ((t493 * t573 - t494 * t509 + t495 * t510) * t523 + (-(-t463 * t509 + t465 * t510) * t542 + (-t464 * t509 + t466 * t510 - t585) * t538) * t560) * t522 / 0.2e1 + ((t449 * t482 + t450 * t483 + t455 * t486 + t456 * t487 + t582 * t509) * t512 + (t418 * t482 + t420 * t483 + t430 * t486 + t432 * t487 + t584 * t509) * t489 + (t419 * t482 + t421 * t483 + t431 * t486 + t433 * t487 + t583 * t509) * t488) * t488 / 0.2e1 + ((t449 * t480 + t450 * t481 + t455 * t484 + t456 * t485 + t582 * t507) * t512 + (t418 * t480 + t420 * t481 + t430 * t484 + t432 * t485 + t584 * t507) * t489 + (t419 * t480 + t421 * t481 + t431 * t484 + t433 * t485 + t583 * t507) * t488) * t489 / 0.2e1 + ((t449 * t496 + t450 * t497 + t455 * t505 + t456 * t506 - t582 * t571) * t512 + (t418 * t496 + t420 * t497 + t430 * t505 + t432 * t506 - t584 * t571) * t489 + (t419 * t496 + t421 * t497 + t431 * t505 + t433 * t506 - t583 * t571) * t488) * t512 / 0.2e1 + (Icges(2,3) + m(2) * (t519 ^ 2 + t520 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
