% Calculate kinetic energy for
% S6RRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 22:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 22:26:30
% EndTime: 2019-03-09 22:26:33
% DurationCPUTime: 3.32s
% Computational Cost: add. (3463->383), mult. (5435->577), div. (0->0), fcn. (6420->14), ass. (0->172)
t586 = Icges(5,3) + Icges(6,3);
t532 = sin(pkin(6));
t533 = cos(pkin(6));
t540 = cos(qJ(2));
t541 = cos(qJ(1));
t565 = t540 * t541;
t536 = sin(qJ(2));
t537 = sin(qJ(1));
t568 = t536 * t537;
t507 = -t533 * t565 + t568;
t566 = t537 * t540;
t567 = t536 * t541;
t508 = t533 * t567 + t566;
t509 = t533 * t566 + t567;
t510 = -t533 * t568 + t565;
t570 = t532 * t541;
t573 = t532 * t537;
t553 = (Icges(3,5) * t508 - Icges(3,6) * t507 - Icges(3,3) * t570) * t541 - (Icges(3,5) * t510 - Icges(3,6) * t509 + Icges(3,3) * t573) * t537;
t585 = t532 * t553;
t531 = qJ(3) + qJ(4);
t556 = pkin(12) + t531;
t524 = sin(t556);
t554 = cos(t556);
t551 = t532 * t554;
t474 = t508 * t524 + t541 * t551;
t475 = t508 * t554 - t524 * t570;
t527 = sin(t531);
t528 = cos(t531);
t480 = -t508 * t527 - t528 * t570;
t481 = t508 * t528 - t527 * t570;
t584 = Icges(5,5) * t481 + Icges(6,5) * t475 + Icges(5,6) * t480 - Icges(6,6) * t474 + t507 * t586;
t476 = t510 * t524 - t537 * t551;
t477 = t510 * t554 + t524 * t573;
t482 = -t510 * t527 + t528 * t573;
t483 = t510 * t528 + t527 * t573;
t583 = Icges(5,5) * t483 + Icges(6,5) * t477 + Icges(5,6) * t482 - Icges(6,6) * t476 + t509 * t586;
t574 = t532 * t536;
t490 = t524 * t574 - t533 * t554;
t491 = t533 * t524 + t536 * t551;
t497 = -t527 * t574 + t528 * t533;
t498 = t527 * t533 + t528 * t574;
t571 = t532 * t540;
t582 = Icges(5,5) * t498 + Icges(6,5) * t491 + Icges(5,6) * t497 - Icges(6,6) * t490 - t571 * t586;
t539 = cos(qJ(3));
t576 = t539 * pkin(3);
t572 = t532 * t539;
t535 = sin(qJ(3));
t569 = t533 * t535;
t478 = pkin(2) * t508 + pkin(9) * t507;
t479 = pkin(2) * t510 + pkin(9) * t509;
t560 = qJD(2) * t532;
t522 = t537 * t560;
t557 = t541 * t560;
t564 = t478 * t522 + t479 * t557;
t563 = pkin(4) * t528;
t488 = qJD(3) * t509 + t522;
t561 = qJD(1) * (pkin(1) * t537 - pkin(8) * t570);
t523 = qJD(2) * t533 + qJD(1);
t559 = t535 * t573;
t558 = t535 * t570;
t457 = qJD(4) * t509 + t488;
t555 = pkin(4) * t527;
t489 = qJD(3) * t507 - t557;
t424 = -pkin(3) * t558 + pkin(10) * t507 + t508 * t576;
t425 = pkin(3) * t559 + pkin(10) * t509 + t510 * t576;
t552 = t488 * t424 - t425 * t489 + t564;
t458 = qJD(4) * t507 + t489;
t511 = (pkin(2) * t536 - pkin(9) * t540) * t532;
t513 = qJD(1) * (pkin(1) * t541 + pkin(8) * t573);
t550 = t523 * t479 - t511 * t522 + t513;
t492 = (-qJD(3) - qJD(4)) * t571 + t523;
t549 = -t478 * t523 - t511 * t557 - t561;
t466 = pkin(3) * t569 + (-pkin(10) * t540 + t536 * t576) * t532;
t512 = -qJD(3) * t571 + t523;
t548 = t512 * t425 - t466 * t488 + t550;
t396 = qJ(5) * t507 + t508 * t563 - t555 * t570;
t547 = -qJD(5) * t571 + t457 * t396 + t552;
t397 = qJ(5) * t509 + t510 * t563 + t555 * t573;
t546 = qJD(5) * t507 + t492 * t397 + t548;
t545 = -t424 * t512 + t489 * t466 + t549;
t438 = t555 * t533 + (-qJ(5) * t540 + t536 * t563) * t532;
t544 = qJD(5) * t509 + t458 * t438 + t545;
t538 = cos(qJ(6));
t534 = sin(qJ(6));
t520 = rSges(2,1) * t541 - rSges(2,2) * t537;
t519 = rSges(2,1) * t537 + rSges(2,2) * t541;
t506 = t536 * t572 + t569;
t505 = t533 * t539 - t535 * t574;
t496 = rSges(3,3) * t533 + (rSges(3,1) * t536 + rSges(3,2) * t540) * t532;
t495 = Icges(3,5) * t533 + (Icges(3,1) * t536 + Icges(3,4) * t540) * t532;
t494 = Icges(3,6) * t533 + (Icges(3,4) * t536 + Icges(3,2) * t540) * t532;
t493 = Icges(3,3) * t533 + (Icges(3,5) * t536 + Icges(3,6) * t540) * t532;
t487 = t510 * t539 + t559;
t486 = -t510 * t535 + t537 * t572;
t485 = t508 * t539 - t558;
t484 = -t508 * t535 - t539 * t570;
t473 = t491 * t538 - t534 * t571;
t472 = -t491 * t534 - t538 * t571;
t469 = rSges(3,1) * t510 - rSges(3,2) * t509 + rSges(3,3) * t573;
t468 = rSges(3,1) * t508 - rSges(3,2) * t507 - rSges(3,3) * t570;
t465 = Icges(3,1) * t510 - Icges(3,4) * t509 + Icges(3,5) * t573;
t464 = Icges(3,1) * t508 - Icges(3,4) * t507 - Icges(3,5) * t570;
t463 = Icges(3,4) * t510 - Icges(3,2) * t509 + Icges(3,6) * t573;
t462 = Icges(3,4) * t508 - Icges(3,2) * t507 - Icges(3,6) * t570;
t459 = rSges(4,1) * t506 + rSges(4,2) * t505 - rSges(4,3) * t571;
t456 = Icges(4,1) * t506 + Icges(4,4) * t505 - Icges(4,5) * t571;
t455 = Icges(4,4) * t506 + Icges(4,2) * t505 - Icges(4,6) * t571;
t454 = Icges(4,5) * t506 + Icges(4,6) * t505 - Icges(4,3) * t571;
t453 = qJD(6) * t490 + t492;
t452 = pkin(5) * t491 + pkin(11) * t490;
t451 = rSges(5,1) * t498 + rSges(5,2) * t497 - rSges(5,3) * t571;
t450 = Icges(5,1) * t498 + Icges(5,4) * t497 - Icges(5,5) * t571;
t449 = Icges(5,4) * t498 + Icges(5,2) * t497 - Icges(5,6) * t571;
t447 = rSges(6,1) * t491 - rSges(6,2) * t490 - rSges(6,3) * t571;
t446 = Icges(6,1) * t491 - Icges(6,4) * t490 - Icges(6,5) * t571;
t445 = Icges(6,4) * t491 - Icges(6,2) * t490 - Icges(6,6) * t571;
t443 = t477 * t538 + t509 * t534;
t442 = -t477 * t534 + t509 * t538;
t441 = t475 * t538 + t507 * t534;
t440 = -t475 * t534 + t507 * t538;
t437 = pkin(5) * t477 + pkin(11) * t476;
t436 = pkin(5) * t475 + pkin(11) * t474;
t435 = rSges(4,1) * t487 + rSges(4,2) * t486 + rSges(4,3) * t509;
t434 = rSges(4,1) * t485 + rSges(4,2) * t484 + rSges(4,3) * t507;
t433 = Icges(4,1) * t487 + Icges(4,4) * t486 + Icges(4,5) * t509;
t432 = Icges(4,1) * t485 + Icges(4,4) * t484 + Icges(4,5) * t507;
t431 = Icges(4,4) * t487 + Icges(4,2) * t486 + Icges(4,6) * t509;
t430 = Icges(4,4) * t485 + Icges(4,2) * t484 + Icges(4,6) * t507;
t429 = Icges(4,5) * t487 + Icges(4,6) * t486 + Icges(4,3) * t509;
t428 = Icges(4,5) * t485 + Icges(4,6) * t484 + Icges(4,3) * t507;
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
t407 = t469 * t523 - t496 * t522 + t513;
t406 = -t468 * t523 - t496 * t557 - t561;
t405 = rSges(7,1) * t473 + rSges(7,2) * t472 + rSges(7,3) * t490;
t403 = Icges(7,1) * t473 + Icges(7,4) * t472 + Icges(7,5) * t490;
t402 = Icges(7,4) * t473 + Icges(7,2) * t472 + Icges(7,6) * t490;
t401 = Icges(7,5) * t473 + Icges(7,6) * t472 + Icges(7,3) * t490;
t400 = (t468 * t537 + t469 * t541) * t560;
t394 = rSges(7,1) * t443 + rSges(7,2) * t442 + rSges(7,3) * t476;
t393 = rSges(7,1) * t441 + rSges(7,2) * t440 + rSges(7,3) * t474;
t392 = Icges(7,1) * t443 + Icges(7,4) * t442 + Icges(7,5) * t476;
t391 = Icges(7,1) * t441 + Icges(7,4) * t440 + Icges(7,5) * t474;
t390 = Icges(7,4) * t443 + Icges(7,2) * t442 + Icges(7,6) * t476;
t389 = Icges(7,4) * t441 + Icges(7,2) * t440 + Icges(7,6) * t474;
t388 = Icges(7,5) * t443 + Icges(7,6) * t442 + Icges(7,3) * t476;
t387 = Icges(7,5) * t441 + Icges(7,6) * t440 + Icges(7,3) * t474;
t385 = t435 * t512 - t459 * t488 + t550;
t384 = -t434 * t512 + t459 * t489 + t549;
t383 = t434 * t488 - t435 * t489 + t564;
t382 = t423 * t492 - t451 * t457 + t548;
t381 = -t422 * t492 + t451 * t458 + t545;
t380 = t422 * t457 - t423 * t458 + t552;
t379 = t415 * t492 + (-t438 - t447) * t457 + t546;
t378 = t447 * t458 + (-t396 - t414) * t492 + t544;
t377 = t414 * t457 + (-t397 - t415) * t458 + t547;
t376 = t394 * t453 - t405 * t426 + t437 * t492 + (-t438 - t452) * t457 + t546;
t375 = -t393 * t453 + t405 * t427 + t452 * t458 + (-t396 - t436) * t492 + t544;
t374 = t393 * t426 - t394 * t427 + t436 * t457 + (-t397 - t437) * t458 + t547;
t1 = ((t493 * t573 - t494 * t509 + t495 * t510) * t523 + (-(-t462 * t509 + t464 * t510) * t541 + (-t509 * t463 + t510 * t465 - t585) * t537) * t560) * t522 / 0.2e1 - ((-t493 * t570 - t494 * t507 + t495 * t508) * t523 + ((-t463 * t507 + t465 * t508) * t537 + (t507 * t462 - t508 * t464 + t585) * t541) * t560) * t557 / 0.2e1 + t523 * ((t533 * t493 + (t494 * t540 + t495 * t536) * t532) * t523 + (((t463 * t540 + t465 * t536) * t537 - (t462 * t540 + t464 * t536) * t541) * t532 - t553 * t533) * t560) / 0.2e1 + t427 * ((t388 * t474 + t390 * t440 + t392 * t441) * t426 + (t387 * t474 + t440 * t389 + t441 * t391) * t427 + (t401 * t474 + t402 * t440 + t403 * t441) * t453) / 0.2e1 + t453 * ((t388 * t490 + t390 * t472 + t392 * t473) * t426 + (t387 * t490 + t389 * t472 + t391 * t473) * t427 + (t490 * t401 + t472 * t402 + t473 * t403) * t453) / 0.2e1 + t426 * ((t388 * t476 + t442 * t390 + t443 * t392) * t426 + (t387 * t476 + t389 * t442 + t391 * t443) * t427 + (t401 * t476 + t402 * t442 + t403 * t443) * t453) / 0.2e1 + t488 * ((t509 * t429 + t486 * t431 + t487 * t433) * t488 + (t428 * t509 + t430 * t486 + t432 * t487) * t489 + (t454 * t509 + t455 * t486 + t456 * t487) * t512) / 0.2e1 + t489 * ((t429 * t507 + t431 * t484 + t433 * t485) * t488 + (t507 * t428 + t484 * t430 + t485 * t432) * t489 + (t454 * t507 + t455 * t484 + t456 * t485) * t512) / 0.2e1 + t512 * ((-t429 * t571 + t431 * t505 + t433 * t506) * t488 + (-t428 * t571 + t430 * t505 + t432 * t506) * t489 + (-t454 * t571 + t505 * t455 + t506 * t456) * t512) / 0.2e1 + m(7) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + m(5) * (t380 ^ 2 + t381 ^ 2 + t382 ^ 2) / 0.2e1 + m(6) * (t377 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + m(4) * (t383 ^ 2 + t384 ^ 2 + t385 ^ 2) / 0.2e1 + m(3) * (t400 ^ 2 + t406 ^ 2 + t407 ^ 2) / 0.2e1 + ((-t445 * t476 + t446 * t477 + t449 * t482 + t450 * t483 + t582 * t509) * t492 + (-t410 * t476 + t412 * t477 + t482 * t418 + t483 * t420 + t584 * t509) * t458 + (-t476 * t411 + t477 * t413 + t482 * t419 + t483 * t421 + t583 * t509) * t457) * t457 / 0.2e1 + ((-t445 * t474 + t446 * t475 + t449 * t480 + t450 * t481 + t582 * t507) * t492 + (-t474 * t410 + t475 * t412 + t480 * t418 + t481 * t420 + t584 * t507) * t458 + (-t411 * t474 + t413 * t475 + t419 * t480 + t421 * t481 + t583 * t507) * t457) * t458 / 0.2e1 + ((-t490 * t445 + t491 * t446 + t497 * t449 + t498 * t450 - t582 * t571) * t492 + (-t410 * t490 + t412 * t491 + t418 * t497 + t420 * t498 - t584 * t571) * t458 + (-t411 * t490 + t413 * t491 + t419 * t497 + t421 * t498 - t583 * t571) * t457) * t492 / 0.2e1 + (Icges(2,3) + m(2) * (t519 ^ 2 + t520 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
