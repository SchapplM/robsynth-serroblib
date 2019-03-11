% Calculate kinetic energy for
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR9_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR9_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR9_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR9_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:55:33
% EndTime: 2019-03-09 10:55:36
% DurationCPUTime: 3.19s
% Computational Cost: add. (3253->382), mult. (5583->562), div. (0->0), fcn. (6724->14), ass. (0->173)
t587 = Icges(5,2) + Icges(6,3);
t532 = sin(qJ(2));
t533 = sin(qJ(1));
t534 = cos(qJ(2));
t535 = cos(qJ(1));
t572 = cos(pkin(6));
t552 = t535 * t572;
t503 = t532 * t533 - t534 * t552;
t504 = t532 * t552 + t533 * t534;
t527 = sin(pkin(6));
t566 = t527 * t535;
t460 = Icges(3,5) * t504 - Icges(3,6) * t503 - Icges(3,3) * t566;
t553 = t533 * t572;
t505 = t535 * t532 + t534 * t553;
t506 = -t532 * t553 + t535 * t534;
t568 = t527 * t533;
t461 = Icges(3,5) * t506 - Icges(3,6) * t505 + Icges(3,3) * t568;
t586 = t527 * (t460 * t535 - t461 * t533);
t559 = pkin(11) + qJ(4);
t522 = sin(t559);
t550 = cos(t559);
t481 = t504 * t550 - t522 * t566;
t525 = sin(pkin(12));
t528 = cos(pkin(12));
t446 = -t481 * t525 + t503 * t528;
t571 = t503 * t525;
t447 = t481 * t528 + t571;
t544 = t527 * t550;
t480 = t504 * t522 + t535 * t544;
t585 = -Icges(5,4) * t481 + Icges(6,5) * t447 - Icges(5,6) * t503 + Icges(6,6) * t446 + t587 * t480;
t483 = t506 * t550 + t522 * t568;
t448 = -t483 * t525 + t505 * t528;
t570 = t505 * t525;
t449 = t483 * t528 + t570;
t482 = t506 * t522 - t533 * t544;
t584 = -Icges(5,4) * t483 + Icges(6,5) * t449 - Icges(5,6) * t505 + Icges(6,6) * t448 + t587 * t482;
t495 = t522 * t572 + t532 * t544;
t567 = t527 * t534;
t478 = -t495 * t525 - t528 * t567;
t558 = t525 * t567;
t479 = t495 * t528 - t558;
t569 = t527 * t532;
t494 = t522 * t569 - t550 * t572;
t583 = -Icges(5,4) * t495 + Icges(6,5) * t479 + Icges(5,6) * t567 + Icges(6,6) * t478 + t587 * t494;
t526 = sin(pkin(11));
t529 = cos(pkin(11));
t484 = -t504 * t526 - t529 * t566;
t556 = t526 * t566;
t485 = t504 * t529 - t556;
t429 = Icges(4,5) * t485 + Icges(4,6) * t484 + Icges(4,3) * t503;
t462 = Icges(3,4) * t504 - Icges(3,2) * t503 - Icges(3,6) * t566;
t582 = -t429 + t462;
t486 = -t506 * t526 + t529 * t568;
t557 = t526 * t568;
t487 = t506 * t529 + t557;
t430 = Icges(4,5) * t487 + Icges(4,6) * t486 + Icges(4,3) * t505;
t463 = Icges(3,4) * t506 - Icges(3,2) * t505 + Icges(3,6) * t568;
t581 = t430 - t463;
t501 = -t526 * t569 + t529 * t572;
t551 = t572 * t526;
t502 = t529 * t569 + t551;
t454 = Icges(4,5) * t502 + Icges(4,6) * t501 - Icges(4,3) * t567;
t492 = Icges(3,6) * t572 + (Icges(3,4) * t532 + Icges(3,2) * t534) * t527;
t580 = t454 - t492;
t574 = pkin(3) * t529;
t573 = pkin(5) * t528;
t474 = pkin(2) * t504 + qJ(3) * t503;
t475 = pkin(2) * t506 + qJ(3) * t505;
t561 = qJD(2) * t527;
t516 = t533 * t561;
t554 = t535 * t561;
t563 = t474 * t516 + t475 * t554;
t488 = qJD(4) * t505 + t516;
t562 = qJD(1) * (pkin(1) * t533 - pkin(8) * t566);
t560 = qJD(3) * t534;
t517 = qJD(2) * t572 + qJD(1);
t509 = qJD(1) * (pkin(1) * t535 + pkin(8) * t568);
t555 = qJD(3) * t503 + t517 * t475 + t509;
t549 = qJD(3) * t505 - t562;
t507 = (pkin(2) * t532 - qJ(3) * t534) * t527;
t546 = (-rSges(4,1) * t502 - rSges(4,2) * t501 + rSges(4,3) * t567 - t507) * t561;
t545 = (-pkin(3) * t551 - (-pkin(9) * t534 + t532 * t574) * t527 - t507) * t561;
t489 = qJD(4) * t503 - t554;
t508 = -qJD(4) * t567 + t517;
t426 = -pkin(3) * t556 + pkin(9) * t503 + t504 * t574;
t427 = pkin(3) * t557 + pkin(9) * t505 + t506 * t574;
t542 = t426 * t516 + t427 * t554 - t527 * t560 + t563;
t438 = pkin(4) * t481 + qJ(5) * t480;
t541 = qJD(5) * t494 + t488 * t438 + t542;
t540 = t517 * t427 + t533 * t545 + t555;
t439 = pkin(4) * t483 + qJ(5) * t482;
t539 = qJD(5) * t480 + t508 * t439 + t540;
t538 = (-t426 - t474) * t517 + t535 * t545 + t549;
t457 = pkin(4) * t495 + qJ(5) * t494;
t537 = qJD(5) * t482 + t489 * t457 + t538;
t524 = pkin(12) + qJ(6);
t523 = cos(t524);
t521 = sin(t524);
t513 = rSges(2,1) * t535 - rSges(2,2) * t533;
t512 = rSges(2,1) * t533 + rSges(2,2) * t535;
t496 = t572 * rSges(3,3) + (rSges(3,1) * t532 + rSges(3,2) * t534) * t527;
t493 = Icges(3,5) * t572 + (Icges(3,1) * t532 + Icges(3,4) * t534) * t527;
t491 = Icges(3,3) * t572 + (Icges(3,5) * t532 + Icges(3,6) * t534) * t527;
t473 = qJD(6) * t494 + t508;
t472 = t495 * t523 - t521 * t567;
t471 = -t495 * t521 - t523 * t567;
t470 = rSges(3,1) * t506 - rSges(3,2) * t505 + rSges(3,3) * t568;
t469 = rSges(3,1) * t504 - rSges(3,2) * t503 - rSges(3,3) * t566;
t465 = Icges(3,1) * t506 - Icges(3,4) * t505 + Icges(3,5) * t568;
t464 = Icges(3,1) * t504 - Icges(3,4) * t503 - Icges(3,5) * t566;
t456 = Icges(4,1) * t502 + Icges(4,4) * t501 - Icges(4,5) * t567;
t455 = Icges(4,4) * t502 + Icges(4,2) * t501 - Icges(4,6) * t567;
t453 = rSges(5,1) * t495 - rSges(5,2) * t494 - rSges(5,3) * t567;
t452 = Icges(5,1) * t495 - Icges(5,4) * t494 - Icges(5,5) * t567;
t450 = Icges(5,5) * t495 - Icges(5,6) * t494 - Icges(5,3) * t567;
t445 = t483 * t523 + t505 * t521;
t444 = -t483 * t521 + t505 * t523;
t443 = t481 * t523 + t503 * t521;
t442 = -t481 * t521 + t503 * t523;
t441 = qJD(6) * t480 + t489;
t440 = qJD(6) * t482 + t488;
t436 = rSges(4,1) * t487 + rSges(4,2) * t486 + rSges(4,3) * t505;
t435 = rSges(4,1) * t485 + rSges(4,2) * t484 + rSges(4,3) * t503;
t434 = Icges(4,1) * t487 + Icges(4,4) * t486 + Icges(4,5) * t505;
t433 = Icges(4,1) * t485 + Icges(4,4) * t484 + Icges(4,5) * t503;
t432 = Icges(4,4) * t487 + Icges(4,2) * t486 + Icges(4,6) * t505;
t431 = Icges(4,4) * t485 + Icges(4,2) * t484 + Icges(4,6) * t503;
t425 = rSges(5,1) * t483 - rSges(5,2) * t482 + rSges(5,3) * t505;
t424 = rSges(5,1) * t481 - rSges(5,2) * t480 + rSges(5,3) * t503;
t423 = Icges(5,1) * t483 - Icges(5,4) * t482 + Icges(5,5) * t505;
t422 = Icges(5,1) * t481 - Icges(5,4) * t480 + Icges(5,5) * t503;
t419 = Icges(5,5) * t483 - Icges(5,6) * t482 + Icges(5,3) * t505;
t418 = Icges(5,5) * t481 - Icges(5,6) * t480 + Icges(5,3) * t503;
t417 = t470 * t517 - t496 * t516 + t509;
t416 = -t469 * t517 - t496 * t554 - t562;
t413 = rSges(6,1) * t479 + rSges(6,2) * t478 + rSges(6,3) * t494;
t411 = Icges(6,1) * t479 + Icges(6,4) * t478 + Icges(6,5) * t494;
t410 = Icges(6,4) * t479 + Icges(6,2) * t478 + Icges(6,6) * t494;
t407 = (t469 * t533 + t470 * t535) * t561;
t406 = rSges(7,1) * t472 + rSges(7,2) * t471 + rSges(7,3) * t494;
t405 = Icges(7,1) * t472 + Icges(7,4) * t471 + Icges(7,5) * t494;
t404 = Icges(7,4) * t472 + Icges(7,2) * t471 + Icges(7,6) * t494;
t403 = Icges(7,5) * t472 + Icges(7,6) * t471 + Icges(7,3) * t494;
t402 = -pkin(5) * t558 + pkin(10) * t494 + t495 * t573;
t401 = rSges(6,1) * t449 + rSges(6,2) * t448 + rSges(6,3) * t482;
t400 = rSges(6,1) * t447 + rSges(6,2) * t446 + rSges(6,3) * t480;
t399 = Icges(6,1) * t449 + Icges(6,4) * t448 + Icges(6,5) * t482;
t398 = Icges(6,1) * t447 + Icges(6,4) * t446 + Icges(6,5) * t480;
t397 = Icges(6,4) * t449 + Icges(6,2) * t448 + Icges(6,6) * t482;
t396 = Icges(6,4) * t447 + Icges(6,2) * t446 + Icges(6,6) * t480;
t393 = rSges(7,1) * t445 + rSges(7,2) * t444 + rSges(7,3) * t482;
t392 = rSges(7,1) * t443 + rSges(7,2) * t442 + rSges(7,3) * t480;
t391 = Icges(7,1) * t445 + Icges(7,4) * t444 + Icges(7,5) * t482;
t390 = Icges(7,1) * t443 + Icges(7,4) * t442 + Icges(7,5) * t480;
t389 = Icges(7,4) * t445 + Icges(7,2) * t444 + Icges(7,6) * t482;
t388 = Icges(7,4) * t443 + Icges(7,2) * t442 + Icges(7,6) * t480;
t387 = Icges(7,5) * t445 + Icges(7,6) * t444 + Icges(7,3) * t482;
t386 = Icges(7,5) * t443 + Icges(7,6) * t442 + Icges(7,3) * t480;
t385 = pkin(5) * t570 + pkin(10) * t482 + t483 * t573;
t384 = pkin(5) * t571 + pkin(10) * t480 + t481 * t573;
t383 = t436 * t517 + t533 * t546 + t555;
t382 = (-t435 - t474) * t517 + t535 * t546 + t549;
t381 = (-t560 + (t435 * t533 + t436 * t535) * qJD(2)) * t527 + t563;
t380 = t425 * t508 - t453 * t488 + t540;
t379 = -t424 * t508 + t453 * t489 + t538;
t378 = t424 * t488 - t425 * t489 + t542;
t377 = t401 * t508 + (-t413 - t457) * t488 + t539;
t376 = t413 * t489 + (-t400 - t438) * t508 + t537;
t375 = t400 * t488 + (-t401 - t439) * t489 + t541;
t374 = t385 * t508 + t393 * t473 - t406 * t440 + (-t402 - t457) * t488 + t539;
t373 = -t392 * t473 + t402 * t489 + t406 * t441 + (-t384 - t438) * t508 + t537;
t372 = t384 * t488 + t392 * t440 - t393 * t441 + (-t385 - t439) * t489 + t541;
t1 = m(5) * (t378 ^ 2 + t379 ^ 2 + t380 ^ 2) / 0.2e1 + m(4) * (t381 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + m(3) * (t407 ^ 2 + t416 ^ 2 + t417 ^ 2) / 0.2e1 + m(7) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(6) * (t375 ^ 2 + t376 ^ 2 + t377 ^ 2) / 0.2e1 + t473 * ((t387 * t494 + t389 * t471 + t391 * t472) * t440 + (t386 * t494 + t388 * t471 + t390 * t472) * t441 + (t403 * t494 + t404 * t471 + t405 * t472) * t473) / 0.2e1 + t441 * ((t387 * t480 + t389 * t442 + t391 * t443) * t440 + (t480 * t386 + t442 * t388 + t443 * t390) * t441 + (t403 * t480 + t404 * t442 + t405 * t443) * t473) / 0.2e1 + t440 * ((t482 * t387 + t444 * t389 + t445 * t391) * t440 + (t386 * t482 + t388 * t444 + t390 * t445) * t441 + (t403 * t482 + t404 * t444 + t405 * t445) * t473) / 0.2e1 + ((t410 * t448 + t411 * t449 + t450 * t505 + t452 * t483 + t482 * t583) * t508 + (t396 * t448 + t398 * t449 + t418 * t505 + t422 * t483 + t482 * t585) * t489 + (t397 * t448 + t399 * t449 + t419 * t505 + t423 * t483 + t584 * t482) * t488) * t488 / 0.2e1 + ((t410 * t446 + t411 * t447 + t450 * t503 + t452 * t481 + t480 * t583) * t508 + (t396 * t446 + t398 * t447 + t418 * t503 + t422 * t481 + t585 * t480) * t489 + (t397 * t446 + t399 * t447 + t419 * t503 + t423 * t481 + t480 * t584) * t488) * t489 / 0.2e1 + ((t410 * t478 + t411 * t479 - t450 * t567 + t452 * t495 + t494 * t583) * t508 + (t396 * t478 + t398 * t479 - t418 * t567 + t422 * t495 + t494 * t585) * t489 + (t397 * t478 + t399 * t479 - t419 * t567 + t423 * t495 + t494 * t584) * t488) * t508 / 0.2e1 + ((t572 * t461 + (t463 * t534 + t465 * t532) * t527) * t516 - (t572 * t460 + (t462 * t534 + t464 * t532) * t527) * t554 + ((t432 * t501 + t434 * t502) * t533 - (t431 * t501 + t433 * t502) * t535 + (t429 * t535 - t430 * t533) * t567) * t561 + (t572 * t491 + (t492 * t534 + t493 * t532) * t527 - t454 * t567 + t455 * t501 + t456 * t502) * t517) * t517 / 0.2e1 + (Icges(2,3) + m(2) * (t512 ^ 2 + t513 ^ 2)) * qJD(1) ^ 2 / 0.2e1 - (((-t431 * t484 - t433 * t485 - t504 * t464 + t503 * t582 + t586) * t535 + (t432 * t484 + t434 * t485 + t465 * t504 + t503 * t581) * t533) * t561 + (t455 * t484 + t456 * t485 - t491 * t566 + t493 * t504 + t503 * t580) * t517) * t554 / 0.2e1 + (((-t431 * t486 - t433 * t487 - t464 * t506 + t505 * t582) * t535 + (t432 * t486 + t434 * t487 + t506 * t465 + t505 * t581 - t586) * t533) * t561 + (t455 * t486 + t456 * t487 + t491 * t568 + t493 * t506 + t505 * t580) * t517) * t516 / 0.2e1;
T  = t1;
