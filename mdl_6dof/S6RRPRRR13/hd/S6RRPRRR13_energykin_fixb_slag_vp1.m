% Calculate kinetic energy for
% S6RRPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR13_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR13_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR13_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR13_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:44:56
% EndTime: 2019-03-09 14:44:59
% DurationCPUTime: 3.15s
% Computational Cost: add. (2295->340), mult. (5289->521), div. (0->0), fcn. (6342->12), ass. (0->159)
t570 = Icges(3,1) + Icges(4,2);
t569 = Icges(3,4) + Icges(4,6);
t568 = Icges(3,5) - Icges(4,4);
t567 = Icges(3,2) + Icges(4,3);
t566 = Icges(3,6) - Icges(4,5);
t565 = Icges(3,3) + Icges(4,1);
t510 = sin(pkin(6));
t511 = cos(pkin(6));
t517 = cos(qJ(2));
t518 = cos(qJ(1));
t542 = t517 * t518;
t514 = sin(qJ(2));
t515 = sin(qJ(1));
t545 = t514 * t515;
t489 = -t511 * t542 + t545;
t543 = t515 * t517;
t544 = t514 * t518;
t490 = t511 * t544 + t543;
t491 = t511 * t543 + t544;
t492 = -t511 * t545 + t542;
t546 = t510 * t518;
t547 = t510 * t515;
t556 = (-t566 * t489 + t568 * t490 - t565 * t546) * t518 + (t566 * t491 - t568 * t492 - t565 * t547) * t515;
t564 = t556 * t510;
t563 = t567 * t491 - t569 * t492 - t566 * t547;
t562 = t567 * t489 - t569 * t490 + t566 * t546;
t561 = -t569 * t491 + t570 * t492 + t568 * t547;
t560 = t569 * t489 - t570 * t490 + t568 * t546;
t559 = t565 * t511 + (t568 * t514 + t566 * t517) * t510;
t558 = t566 * t511 + (t569 * t514 + t567 * t517) * t510;
t557 = t568 * t511 + (t570 * t514 + t569 * t517) * t510;
t553 = cos(qJ(4));
t516 = cos(qJ(5));
t552 = pkin(5) * t516;
t512 = sin(qJ(5));
t550 = t490 * t512;
t549 = t492 * t512;
t548 = t510 * t514;
t453 = pkin(2) * t490 + qJ(3) * t489;
t454 = pkin(2) * t492 + qJ(3) * t491;
t539 = qJD(2) * t510;
t503 = t515 * t539;
t534 = t518 * t539;
t541 = t453 * t503 + t454 * t534;
t469 = qJD(4) * t492 + t503;
t540 = qJD(1) * (pkin(1) * t515 - pkin(8) * t546);
t538 = qJD(3) * t517;
t504 = qJD(2) * t511 + qJD(1);
t537 = t512 * t548;
t495 = qJD(1) * (pkin(1) * t518 + pkin(8) * t547);
t536 = qJD(3) * t489 + t504 * t454 + t495;
t513 = sin(qJ(4));
t465 = -t491 * t553 + t513 * t547;
t418 = qJD(5) * t465 + t469;
t535 = t510 * t553;
t494 = qJD(4) * t548 + t504;
t533 = qJD(3) * t491 - t540;
t487 = t511 * t513 + t517 * t535;
t455 = qJD(5) * t487 + t494;
t493 = (pkin(2) * t514 - qJ(3) * t517) * t510;
t530 = (-rSges(4,1) * t511 - (-rSges(4,2) * t514 - rSges(4,3) * t517) * t510 - t493) * t539;
t529 = (-pkin(3) * t511 - pkin(9) * t548 - t493) * t539;
t470 = qJD(4) * t490 - t534;
t467 = t489 * t553 + t513 * t546;
t419 = -qJD(5) * t467 + t470;
t471 = pkin(3) * t547 + pkin(9) * t492;
t472 = -pkin(3) * t546 + pkin(9) * t490;
t526 = t471 * t534 + t472 * t503 - t510 * t538 + t541;
t525 = t504 * t471 + t515 * t529 + t536;
t466 = t491 * t513 + t515 * t535;
t416 = pkin(4) * t466 + pkin(10) * t465;
t468 = t489 * t513 - t518 * t535;
t417 = pkin(4) * t468 - pkin(10) * t467;
t524 = -t416 * t470 + t469 * t417 + t526;
t488 = -t510 * t517 * t513 + t511 * t553;
t452 = pkin(4) * t488 + pkin(10) * t487;
t523 = t494 * t416 - t452 * t469 + t525;
t522 = (-t453 - t472) * t504 + t518 * t529 + t533;
t521 = -t417 * t494 + t470 * t452 + t522;
t509 = qJ(5) + qJ(6);
t508 = cos(t509);
t507 = sin(t509);
t499 = rSges(2,1) * t518 - rSges(2,2) * t515;
t498 = rSges(2,1) * t515 + rSges(2,2) * t518;
t479 = rSges(3,3) * t511 + (rSges(3,1) * t514 + rSges(3,2) * t517) * t510;
t464 = t488 * t516 + t537;
t463 = -t488 * t512 + t516 * t548;
t458 = t488 * t508 + t507 * t548;
t457 = -t488 * t507 + t508 * t548;
t451 = rSges(3,1) * t492 - rSges(3,2) * t491 + rSges(3,3) * t547;
t450 = rSges(3,1) * t490 - rSges(3,2) * t489 - rSges(3,3) * t546;
t449 = -rSges(4,1) * t546 - rSges(4,2) * t490 + rSges(4,3) * t489;
t448 = rSges(4,1) * t547 - rSges(4,2) * t492 + rSges(4,3) * t491;
t432 = rSges(5,1) * t488 - rSges(5,2) * t487 + rSges(5,3) * t548;
t431 = Icges(5,1) * t488 - Icges(5,4) * t487 + Icges(5,5) * t548;
t430 = Icges(5,4) * t488 - Icges(5,2) * t487 + Icges(5,6) * t548;
t429 = Icges(5,5) * t488 - Icges(5,6) * t487 + Icges(5,3) * t548;
t428 = qJD(6) * t487 + t455;
t427 = t468 * t516 + t550;
t426 = -t468 * t512 + t490 * t516;
t425 = t466 * t516 + t549;
t424 = -t466 * t512 + t492 * t516;
t423 = t468 * t508 + t490 * t507;
t422 = -t468 * t507 + t490 * t508;
t421 = t466 * t508 + t492 * t507;
t420 = -t466 * t507 + t492 * t508;
t413 = rSges(5,1) * t468 + rSges(5,2) * t467 + rSges(5,3) * t490;
t412 = rSges(5,1) * t466 - rSges(5,2) * t465 + rSges(5,3) * t492;
t411 = Icges(5,1) * t468 + Icges(5,4) * t467 + Icges(5,5) * t490;
t410 = Icges(5,1) * t466 - Icges(5,4) * t465 + Icges(5,5) * t492;
t409 = Icges(5,4) * t468 + Icges(5,2) * t467 + Icges(5,6) * t490;
t408 = Icges(5,4) * t466 - Icges(5,2) * t465 + Icges(5,6) * t492;
t407 = Icges(5,5) * t468 + Icges(5,6) * t467 + Icges(5,3) * t490;
t406 = Icges(5,5) * t466 - Icges(5,6) * t465 + Icges(5,3) * t492;
t405 = rSges(6,1) * t464 + rSges(6,2) * t463 + rSges(6,3) * t487;
t404 = Icges(6,1) * t464 + Icges(6,4) * t463 + Icges(6,5) * t487;
t403 = Icges(6,4) * t464 + Icges(6,2) * t463 + Icges(6,6) * t487;
t402 = Icges(6,5) * t464 + Icges(6,6) * t463 + Icges(6,3) * t487;
t401 = pkin(5) * t537 + pkin(11) * t487 + t488 * t552;
t400 = rSges(7,1) * t458 + rSges(7,2) * t457 + rSges(7,3) * t487;
t399 = -qJD(6) * t467 + t419;
t398 = qJD(6) * t465 + t418;
t396 = Icges(7,1) * t458 + Icges(7,4) * t457 + Icges(7,5) * t487;
t395 = Icges(7,4) * t458 + Icges(7,2) * t457 + Icges(7,6) * t487;
t394 = Icges(7,5) * t458 + Icges(7,6) * t457 + Icges(7,3) * t487;
t393 = t451 * t504 - t479 * t503 + t495;
t392 = -t450 * t504 - t479 * t534 - t540;
t391 = (t450 * t515 + t451 * t518) * t539;
t390 = rSges(6,1) * t427 + rSges(6,2) * t426 - rSges(6,3) * t467;
t389 = rSges(6,1) * t425 + rSges(6,2) * t424 + rSges(6,3) * t465;
t388 = Icges(6,1) * t427 + Icges(6,4) * t426 - Icges(6,5) * t467;
t387 = Icges(6,1) * t425 + Icges(6,4) * t424 + Icges(6,5) * t465;
t386 = Icges(6,4) * t427 + Icges(6,2) * t426 - Icges(6,6) * t467;
t385 = Icges(6,4) * t425 + Icges(6,2) * t424 + Icges(6,6) * t465;
t384 = Icges(6,5) * t427 + Icges(6,6) * t426 - Icges(6,3) * t467;
t383 = Icges(6,5) * t425 + Icges(6,6) * t424 + Icges(6,3) * t465;
t382 = rSges(7,1) * t423 + rSges(7,2) * t422 - rSges(7,3) * t467;
t381 = rSges(7,1) * t421 + rSges(7,2) * t420 + rSges(7,3) * t465;
t380 = Icges(7,1) * t423 + Icges(7,4) * t422 - Icges(7,5) * t467;
t379 = Icges(7,1) * t421 + Icges(7,4) * t420 + Icges(7,5) * t465;
t378 = Icges(7,4) * t423 + Icges(7,2) * t422 - Icges(7,6) * t467;
t377 = Icges(7,4) * t421 + Icges(7,2) * t420 + Icges(7,6) * t465;
t376 = Icges(7,5) * t423 + Icges(7,6) * t422 - Icges(7,3) * t467;
t375 = Icges(7,5) * t421 + Icges(7,6) * t420 + Icges(7,3) * t465;
t374 = pkin(5) * t550 - pkin(11) * t467 + t468 * t552;
t373 = pkin(5) * t549 + pkin(11) * t465 + t466 * t552;
t372 = t448 * t504 + t515 * t530 + t536;
t371 = (-t449 - t453) * t504 + t518 * t530 + t533;
t370 = (-t538 + (t448 * t518 + t449 * t515) * qJD(2)) * t510 + t541;
t369 = t412 * t494 - t432 * t469 + t525;
t368 = -t413 * t494 + t432 * t470 + t522;
t367 = -t412 * t470 + t413 * t469 + t526;
t366 = t389 * t455 - t405 * t418 + t523;
t365 = -t390 * t455 + t405 * t419 + t521;
t364 = -t389 * t419 + t390 * t418 + t524;
t363 = t373 * t455 + t381 * t428 - t398 * t400 - t401 * t418 + t523;
t362 = -t374 * t455 - t382 * t428 + t399 * t400 + t401 * t419 + t521;
t361 = -t373 * t419 + t374 * t418 - t381 * t399 + t382 * t398 + t524;
t1 = m(7) * (t361 ^ 2 + t362 ^ 2 + t363 ^ 2) / 0.2e1 + m(6) * (t364 ^ 2 + t365 ^ 2 + t366 ^ 2) / 0.2e1 + m(5) * (t367 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + m(4) * (t370 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + m(3) * (t391 ^ 2 + t392 ^ 2 + t393 ^ 2) / 0.2e1 + t398 * ((t465 * t375 + t420 * t377 + t421 * t379) * t398 + (t376 * t465 + t378 * t420 + t380 * t421) * t399 + (t394 * t465 + t395 * t420 + t396 * t421) * t428) / 0.2e1 + t455 * ((t383 * t487 + t385 * t463 + t387 * t464) * t418 + (t384 * t487 + t386 * t463 + t388 * t464) * t419 + (t487 * t402 + t463 * t403 + t464 * t404) * t455) / 0.2e1 + t399 * ((-t375 * t467 + t377 * t422 + t379 * t423) * t398 + (-t467 * t376 + t422 * t378 + t423 * t380) * t399 + (-t394 * t467 + t395 * t422 + t396 * t423) * t428) / 0.2e1 + t428 * ((t375 * t487 + t377 * t457 + t379 * t458) * t398 + (t376 * t487 + t378 * t457 + t380 * t458) * t399 + (t487 * t394 + t457 * t395 + t458 * t396) * t428) / 0.2e1 + t419 * ((-t383 * t467 + t385 * t426 + t387 * t427) * t418 + (-t467 * t384 + t426 * t386 + t427 * t388) * t419 + (-t402 * t467 + t403 * t426 + t404 * t427) * t455) / 0.2e1 + t418 * ((t465 * t383 + t424 * t385 + t425 * t387) * t418 + (t384 * t465 + t386 * t424 + t388 * t425) * t419 + (t402 * t465 + t403 * t424 + t404 * t425) * t455) / 0.2e1 + t470 * ((t406 * t490 + t408 * t467 + t410 * t468) * t469 + (t490 * t407 + t467 * t409 + t468 * t411) * t470 + (t429 * t490 + t430 * t467 + t431 * t468) * t494) / 0.2e1 + t494 * ((t406 * t548 - t408 * t487 + t410 * t488) * t469 + (t407 * t548 - t409 * t487 + t411 * t488) * t470 + (t429 * t548 - t487 * t430 + t488 * t431) * t494) / 0.2e1 + t469 * ((t492 * t406 - t465 * t408 + t466 * t410) * t469 + (t407 * t492 - t409 * t465 + t411 * t466) * t470 + (t429 * t492 - t430 * t465 + t431 * t466) * t494) / 0.2e1 + ((-t556 * t511 + ((t514 * t560 + t517 * t562) * t518 + (t561 * t514 - t563 * t517) * t515) * t510) * t539 + (t559 * t511 + (t514 * t557 + t517 * t558) * t510) * t504) * t504 / 0.2e1 + (Icges(2,3) + m(2) * (t498 ^ 2 + t499 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((-t491 * t562 + t492 * t560) * t518 + (t491 * t563 + t561 * t492 - t564) * t515) * t539 + (-t491 * t558 + t557 * t492 + t559 * t547) * t504) * t503 / 0.2e1 - (((-t562 * t489 + t490 * t560 + t564) * t518 + (t489 * t563 + t561 * t490) * t515) * t539 + (-t489 * t558 + t490 * t557 - t546 * t559) * t504) * t534 / 0.2e1;
T  = t1;
