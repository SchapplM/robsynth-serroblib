% Calculate kinetic energy for
% S6PRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR3_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:31:40
% EndTime: 2019-03-08 20:31:43
% DurationCPUTime: 2.95s
% Computational Cost: add. (3294->382), mult. (5183->596), div. (0->0), fcn. (6146->14), ass. (0->170)
t522 = sin(pkin(12));
t525 = cos(pkin(12));
t527 = cos(pkin(6));
t524 = sin(pkin(6));
t530 = sin(qJ(2));
t563 = t524 * t530;
t500 = -t522 * t563 + t525 * t527;
t565 = t522 * t527;
t501 = t525 * t563 + t565;
t532 = cos(qJ(2));
t562 = t524 * t532;
t451 = Icges(4,5) * t501 + Icges(4,6) * t500 - Icges(4,3) * t562;
t489 = Icges(3,6) * t527 + (Icges(3,4) * t530 + Icges(3,2) * t532) * t524;
t569 = t451 - t489;
t568 = qJD(2) ^ 2;
t566 = t525 * pkin(3);
t523 = sin(pkin(11));
t564 = t523 * t524;
t526 = cos(pkin(11));
t561 = t526 * t524;
t560 = t527 * t530;
t559 = t527 * t532;
t504 = t523 * t559 + t526 * t530;
t505 = -t523 * t560 + t526 * t532;
t468 = pkin(2) * t505 + qJ(3) * t504;
t502 = t523 * t530 - t526 * t559;
t516 = qJD(2) * t527;
t557 = qJD(3) * t502 + t468 * t516;
t521 = pkin(12) + qJ(4);
t518 = cos(t521);
t556 = pkin(4) * t518;
t554 = qJD(2) * t524;
t513 = t523 * t554;
t483 = qJD(4) * t504 + t513;
t553 = qJD(3) * t532;
t551 = t522 * t564;
t550 = t522 * t561;
t449 = qJD(5) * t504 + t483;
t549 = t526 * t554;
t503 = t523 * t532 + t526 * t560;
t467 = pkin(2) * t503 + qJ(3) * t502;
t548 = t467 * t513 + t468 * t549 + qJD(1);
t546 = qJ(5) + t521;
t517 = sin(t521);
t545 = pkin(4) * t517;
t506 = (pkin(2) * t530 - qJ(3) * t532) * t524;
t544 = (-rSges(4,1) * t501 - rSges(4,2) * t500 + rSges(4,3) * t562 - t506) * t524;
t543 = (-pkin(3) * t565 - (-pkin(8) * t532 + t530 * t566) * t524 - t506) * t524;
t540 = cos(t546);
t484 = qJD(4) * t502 - t549;
t450 = qJD(5) * t502 + t484;
t539 = t524 * t540;
t487 = t516 + (-qJD(4) - qJD(5)) * t562;
t418 = -pkin(3) * t550 + pkin(8) * t502 + t503 * t566;
t419 = pkin(3) * t551 + pkin(8) * t504 + t505 * t566;
t538 = t418 * t513 + t419 * t549 - t524 * t553 + t548;
t537 = qJD(2) * t523 * t543 + t419 * t516 + t557;
t391 = pkin(9) * t502 + t503 * t556 - t545 * t561;
t392 = pkin(9) * t504 + t505 * t556 + t545 * t564;
t536 = t391 * t483 - t392 * t484 + t538;
t499 = qJD(3) * t504;
t535 = t499 + ((-t418 - t467) * t527 + t526 * t543) * qJD(2);
t434 = t545 * t527 + (-pkin(9) * t532 + t530 * t556) * t524;
t507 = -qJD(4) * t562 + t516;
t534 = t392 * t507 - t434 * t483 + t537;
t533 = -t391 * t507 + t434 * t484 + t535;
t531 = cos(qJ(6));
t529 = sin(qJ(6));
t514 = sin(t546);
t493 = t527 * rSges(3,3) + (rSges(3,1) * t530 + rSges(3,2) * t532) * t524;
t492 = t517 * t527 + t518 * t563;
t491 = -t517 * t563 + t518 * t527;
t490 = Icges(3,5) * t527 + (Icges(3,1) * t530 + Icges(3,4) * t532) * t524;
t488 = Icges(3,3) * t527 + (Icges(3,5) * t530 + Icges(3,6) * t532) * t524;
t486 = t514 * t527 + t530 * t539;
t485 = t514 * t563 - t527 * t540;
t482 = t505 * t525 + t551;
t481 = -t505 * t522 + t525 * t564;
t480 = t503 * t525 - t550;
t479 = -t503 * t522 - t525 * t561;
t478 = t505 * t518 + t517 * t564;
t477 = -t505 * t517 + t518 * t564;
t476 = t503 * t518 - t517 * t561;
t475 = -t503 * t517 - t518 * t561;
t474 = t486 * t531 - t529 * t562;
t473 = -t486 * t529 - t531 * t562;
t472 = t505 * t540 + t514 * t564;
t471 = t505 * t514 - t523 * t539;
t470 = t503 * t540 - t514 * t561;
t469 = t503 * t514 + t526 * t539;
t465 = rSges(3,1) * t505 - rSges(3,2) * t504 + rSges(3,3) * t564;
t464 = rSges(3,1) * t503 - rSges(3,2) * t502 - rSges(3,3) * t561;
t459 = Icges(3,1) * t505 - Icges(3,4) * t504 + Icges(3,5) * t564;
t458 = Icges(3,1) * t503 - Icges(3,4) * t502 - Icges(3,5) * t561;
t457 = Icges(3,4) * t505 - Icges(3,2) * t504 + Icges(3,6) * t564;
t456 = Icges(3,4) * t503 - Icges(3,2) * t502 - Icges(3,6) * t561;
t455 = Icges(3,5) * t505 - Icges(3,6) * t504 + Icges(3,3) * t564;
t454 = Icges(3,5) * t503 - Icges(3,6) * t502 - Icges(3,3) * t561;
t453 = Icges(4,1) * t501 + Icges(4,4) * t500 - Icges(4,5) * t562;
t452 = Icges(4,4) * t501 + Icges(4,2) * t500 - Icges(4,6) * t562;
t448 = qJD(6) * t485 + t487;
t447 = pkin(5) * t486 + pkin(10) * t485;
t446 = rSges(5,1) * t492 + rSges(5,2) * t491 - rSges(5,3) * t562;
t445 = Icges(5,1) * t492 + Icges(5,4) * t491 - Icges(5,5) * t562;
t444 = Icges(5,4) * t492 + Icges(5,2) * t491 - Icges(5,6) * t562;
t443 = Icges(5,5) * t492 + Icges(5,6) * t491 - Icges(5,3) * t562;
t442 = rSges(6,1) * t486 - rSges(6,2) * t485 - rSges(6,3) * t562;
t441 = Icges(6,1) * t486 - Icges(6,4) * t485 - Icges(6,5) * t562;
t440 = Icges(6,4) * t486 - Icges(6,2) * t485 - Icges(6,6) * t562;
t439 = Icges(6,5) * t486 - Icges(6,6) * t485 - Icges(6,3) * t562;
t438 = t472 * t531 + t504 * t529;
t437 = -t472 * t529 + t504 * t531;
t436 = t470 * t531 + t502 * t529;
t435 = -t470 * t529 + t502 * t531;
t433 = pkin(5) * t472 + pkin(10) * t471;
t432 = pkin(5) * t470 + pkin(10) * t469;
t431 = (-t464 * t527 - t493 * t561) * qJD(2);
t430 = (t465 * t527 - t493 * t564) * qJD(2);
t429 = rSges(4,1) * t482 + rSges(4,2) * t481 + rSges(4,3) * t504;
t428 = rSges(4,1) * t480 + rSges(4,2) * t479 + rSges(4,3) * t502;
t427 = Icges(4,1) * t482 + Icges(4,4) * t481 + Icges(4,5) * t504;
t426 = Icges(4,1) * t480 + Icges(4,4) * t479 + Icges(4,5) * t502;
t425 = Icges(4,4) * t482 + Icges(4,2) * t481 + Icges(4,6) * t504;
t424 = Icges(4,4) * t480 + Icges(4,2) * t479 + Icges(4,6) * t502;
t423 = Icges(4,5) * t482 + Icges(4,6) * t481 + Icges(4,3) * t504;
t422 = Icges(4,5) * t480 + Icges(4,6) * t479 + Icges(4,3) * t502;
t421 = qJD(6) * t469 + t450;
t420 = qJD(6) * t471 + t449;
t417 = rSges(5,1) * t478 + rSges(5,2) * t477 + rSges(5,3) * t504;
t416 = rSges(5,1) * t476 + rSges(5,2) * t475 + rSges(5,3) * t502;
t415 = Icges(5,1) * t478 + Icges(5,4) * t477 + Icges(5,5) * t504;
t414 = Icges(5,1) * t476 + Icges(5,4) * t475 + Icges(5,5) * t502;
t413 = Icges(5,4) * t478 + Icges(5,2) * t477 + Icges(5,6) * t504;
t412 = Icges(5,4) * t476 + Icges(5,2) * t475 + Icges(5,6) * t502;
t411 = Icges(5,5) * t478 + Icges(5,6) * t477 + Icges(5,3) * t504;
t410 = Icges(5,5) * t476 + Icges(5,6) * t475 + Icges(5,3) * t502;
t408 = rSges(6,1) * t472 - rSges(6,2) * t471 + rSges(6,3) * t504;
t407 = rSges(6,1) * t470 - rSges(6,2) * t469 + rSges(6,3) * t502;
t406 = Icges(6,1) * t472 - Icges(6,4) * t471 + Icges(6,5) * t504;
t405 = Icges(6,1) * t470 - Icges(6,4) * t469 + Icges(6,5) * t502;
t404 = Icges(6,4) * t472 - Icges(6,2) * t471 + Icges(6,6) * t504;
t403 = Icges(6,4) * t470 - Icges(6,2) * t469 + Icges(6,6) * t502;
t402 = Icges(6,5) * t472 - Icges(6,6) * t471 + Icges(6,3) * t504;
t401 = Icges(6,5) * t470 - Icges(6,6) * t469 + Icges(6,3) * t502;
t398 = rSges(7,1) * t474 + rSges(7,2) * t473 + rSges(7,3) * t485;
t397 = Icges(7,1) * t474 + Icges(7,4) * t473 + Icges(7,5) * t485;
t396 = Icges(7,4) * t474 + Icges(7,2) * t473 + Icges(7,6) * t485;
t395 = Icges(7,5) * t474 + Icges(7,6) * t473 + Icges(7,3) * t485;
t393 = qJD(1) + (t464 * t523 + t465 * t526) * t554;
t388 = rSges(7,1) * t438 + rSges(7,2) * t437 + rSges(7,3) * t471;
t387 = rSges(7,1) * t436 + rSges(7,2) * t435 + rSges(7,3) * t469;
t386 = Icges(7,1) * t438 + Icges(7,4) * t437 + Icges(7,5) * t471;
t385 = Icges(7,1) * t436 + Icges(7,4) * t435 + Icges(7,5) * t469;
t384 = Icges(7,4) * t438 + Icges(7,2) * t437 + Icges(7,6) * t471;
t383 = Icges(7,4) * t436 + Icges(7,2) * t435 + Icges(7,6) * t469;
t382 = Icges(7,5) * t438 + Icges(7,6) * t437 + Icges(7,3) * t471;
t381 = Icges(7,5) * t436 + Icges(7,6) * t435 + Icges(7,3) * t469;
t380 = t499 + ((-t428 - t467) * t527 + t526 * t544) * qJD(2);
t379 = (t429 * t527 + t523 * t544) * qJD(2) + t557;
t378 = (-t553 + (t428 * t523 + t429 * t526) * qJD(2)) * t524 + t548;
t377 = -t416 * t507 + t446 * t484 + t535;
t376 = t417 * t507 - t446 * t483 + t537;
t375 = t416 * t483 - t417 * t484 + t538;
t374 = -t407 * t487 + t442 * t450 + t533;
t373 = t408 * t487 - t442 * t449 + t534;
t372 = t407 * t449 - t408 * t450 + t536;
t371 = -t387 * t448 + t398 * t421 - t432 * t487 + t447 * t450 + t533;
t370 = t388 * t448 - t398 * t420 + t433 * t487 - t447 * t449 + t534;
t369 = t387 * t420 - t388 * t421 + t432 * t449 - t433 * t450 + t536;
t1 = m(4) * (t378 ^ 2 + t379 ^ 2 + t380 ^ 2) / 0.2e1 + m(3) * (t393 ^ 2 + t430 ^ 2 + t431 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t421 * ((t382 * t469 + t384 * t435 + t386 * t436) * t420 + (t469 * t381 + t435 * t383 + t436 * t385) * t421 + (t395 * t469 + t396 * t435 + t397 * t436) * t448) / 0.2e1 + t448 * ((t382 * t485 + t384 * t473 + t386 * t474) * t420 + (t381 * t485 + t383 * t473 + t385 * t474) * t421 + (t485 * t395 + t473 * t396 + t474 * t397) * t448) / 0.2e1 + t420 * ((t471 * t382 + t437 * t384 + t438 * t386) * t420 + (t381 * t471 + t383 * t437 + t385 * t438) * t421 + (t395 * t471 + t396 * t437 + t397 * t438) * t448) / 0.2e1 + t449 * ((t504 * t402 - t471 * t404 + t472 * t406) * t449 + (t401 * t504 - t403 * t471 + t405 * t472) * t450 + (t439 * t504 - t440 * t471 + t441 * t472) * t487) / 0.2e1 + t487 * ((-t402 * t562 - t404 * t485 + t406 * t486) * t449 + (-t401 * t562 - t403 * t485 + t405 * t486) * t450 + (-t439 * t562 - t485 * t440 + t486 * t441) * t487) / 0.2e1 + t450 * ((t402 * t502 - t404 * t469 + t406 * t470) * t449 + (t502 * t401 - t469 * t403 + t470 * t405) * t450 + (t439 * t502 - t440 * t469 + t441 * t470) * t487) / 0.2e1 + t483 * ((t504 * t411 + t477 * t413 + t478 * t415) * t483 + (t410 * t504 + t412 * t477 + t414 * t478) * t484 + (t443 * t504 + t444 * t477 + t445 * t478) * t507) / 0.2e1 + t484 * ((t411 * t502 + t413 * t475 + t415 * t476) * t483 + (t502 * t410 + t475 * t412 + t476 * t414) * t484 + (t443 * t502 + t444 * t475 + t445 * t476) * t507) / 0.2e1 + t507 * ((-t411 * t562 + t413 * t491 + t415 * t492) * t483 + (-t410 * t562 + t412 * t491 + t414 * t492) * t484 + (-t443 * t562 + t491 * t444 + t492 * t445) * t507) / 0.2e1 + m(7) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + m(6) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(5) * (t375 ^ 2 + t376 ^ 2 + t377 ^ 2) / 0.2e1 - ((-t455 * t561 - t457 * t502 + t459 * t503) * t564 - (-t454 * t561 - t456 * t502 + t458 * t503) * t561 + ((t423 * t502 + t425 * t479 + t427 * t480) * t523 - (t422 * t502 + t424 * t479 + t426 * t480) * t526) * t524 + (t452 * t479 + t453 * t480 - t488 * t561 + t490 * t503 + t502 * t569) * t527) * t568 * t561 / 0.2e1 + (((-t423 * t562 + t425 * t500 + t427 * t501) * t564 - (-t422 * t562 + t424 * t500 + t426 * t501) * t561 + ((t457 * t532 + t459 * t530) * t523 - (t456 * t532 + t458 * t530) * t526) * t524 ^ 2 + (-t451 * t562 + t500 * t452 + t501 * t453 + (-t454 * t526 + t455 * t523 + t489 * t532 + t490 * t530) * t524 + t527 * t488) * t527) * t527 + ((t455 * t564 - t457 * t504 + t459 * t505) * t564 - (t454 * t564 - t456 * t504 + t458 * t505) * t561 + ((t423 * t504 + t425 * t481 + t427 * t482) * t523 - (t422 * t504 + t424 * t481 + t426 * t482) * t526) * t524 + (t452 * t481 + t453 * t482 + t488 * t564 + t490 * t505 + t504 * t569) * t527) * t564) * t568 / 0.2e1;
T  = t1;
