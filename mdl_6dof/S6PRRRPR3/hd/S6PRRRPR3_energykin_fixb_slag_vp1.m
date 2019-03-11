% Calculate kinetic energy for
% S6PRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR3_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:10:56
% EndTime: 2019-03-08 23:10:59
% DurationCPUTime: 2.86s
% Computational Cost: add. (2814->335), mult. (5255->512), div. (0->0), fcn. (6266->12), ass. (0->156)
t565 = Icges(5,1) + Icges(6,2);
t564 = Icges(6,1) + Icges(5,3);
t563 = -Icges(5,4) - Icges(6,6);
t562 = Icges(6,4) - Icges(5,5);
t561 = Icges(6,5) - Icges(5,6);
t560 = Icges(5,2) + Icges(6,3);
t507 = sin(pkin(11));
t509 = cos(pkin(11));
t516 = cos(qJ(2));
t510 = cos(pkin(6));
t513 = sin(qJ(2));
t538 = t510 * t513;
t493 = t507 * t516 + t509 * t538;
t508 = sin(pkin(6));
t536 = qJ(3) + qJ(4);
t529 = cos(t536);
t527 = t508 * t529;
t528 = sin(t536);
t467 = t493 * t528 + t509 * t527;
t526 = t508 * t528;
t468 = t493 * t529 - t509 * t526;
t537 = t510 * t516;
t492 = t507 * t513 - t509 * t537;
t559 = t560 * t467 + t563 * t468 + t561 * t492;
t495 = -t507 * t538 + t509 * t516;
t469 = t495 * t528 - t507 * t527;
t470 = t495 * t529 + t507 * t526;
t494 = t507 * t537 + t509 * t513;
t558 = t560 * t469 + t563 * t470 + t561 * t494;
t557 = t561 * t467 - t562 * t468 + t564 * t492;
t556 = t561 * t469 - t562 * t470 + t564 * t494;
t555 = t563 * t467 + t565 * t468 - t562 * t492;
t554 = t563 * t469 + t565 * t470 - t562 * t494;
t486 = -t510 * t529 + t513 * t526;
t487 = t510 * t528 + t513 * t527;
t540 = t508 * t516;
t553 = t560 * t486 + t563 * t487 - t561 * t540;
t552 = t563 * t486 + t565 * t487 + t562 * t540;
t551 = t561 * t486 - t562 * t487 - t564 * t540;
t550 = qJD(2) ^ 2;
t515 = cos(qJ(3));
t546 = pkin(3) * t515;
t544 = t507 * t508;
t543 = t508 * t509;
t512 = sin(qJ(3));
t542 = t508 * t512;
t541 = t508 * t515;
t539 = t510 * t512;
t535 = qJD(2) * t508;
t504 = t507 * t535;
t478 = qJD(3) * t494 + t504;
t506 = qJD(2) * t510;
t533 = t507 * t542;
t532 = t509 * t542;
t443 = qJD(4) * t494 + t478;
t531 = t509 * t535;
t463 = t493 * pkin(2) + t492 * pkin(8);
t464 = t495 * pkin(2) + t494 * pkin(8);
t530 = t463 * t504 + t464 * t531 + qJD(1);
t479 = qJD(3) * t492 - t531;
t498 = (pkin(2) * t513 - pkin(8) * t516) * t508;
t525 = t464 * t506 - t498 * t504;
t444 = qJD(4) * t492 + t479;
t481 = t506 + (-qJD(3) - qJD(4)) * t540;
t412 = -pkin(3) * t532 + pkin(9) * t492 + t493 * t546;
t413 = pkin(3) * t533 + pkin(9) * t494 + t495 * t546;
t524 = t478 * t412 - t413 * t479 + t530;
t523 = (-t463 * t510 - t498 * t543) * qJD(2);
t426 = pkin(4) * t468 + qJ(5) * t467;
t522 = qJD(5) * t486 + t443 * t426 + t524;
t460 = pkin(3) * t539 + (-pkin(9) * t516 + t513 * t546) * t508;
t499 = -qJD(3) * t540 + t506;
t521 = t499 * t413 - t460 * t478 + t525;
t427 = pkin(4) * t470 + qJ(5) * t469;
t520 = qJD(5) * t467 + t481 * t427 + t521;
t519 = -t412 * t499 + t479 * t460 + t523;
t458 = pkin(4) * t487 + qJ(5) * t486;
t518 = qJD(5) * t469 + t444 * t458 + t519;
t514 = cos(qJ(6));
t511 = sin(qJ(6));
t497 = t513 * t541 + t539;
t496 = t510 * t515 - t513 * t542;
t485 = rSges(3,3) * t510 + (rSges(3,1) * t513 + rSges(3,2) * t516) * t508;
t484 = Icges(3,5) * t510 + (Icges(3,1) * t513 + Icges(3,4) * t516) * t508;
t483 = Icges(3,6) * t510 + (Icges(3,4) * t513 + Icges(3,2) * t516) * t508;
t482 = Icges(3,3) * t510 + (Icges(3,5) * t513 + Icges(3,6) * t516) * t508;
t477 = -pkin(5) * t540 + pkin(10) * t487;
t476 = t495 * t515 + t533;
t475 = -t495 * t512 + t507 * t541;
t474 = t493 * t515 - t532;
t473 = -t493 * t512 - t509 * t541;
t472 = t486 * t511 - t514 * t540;
t471 = t486 * t514 + t511 * t540;
t461 = qJD(6) * t487 + t481;
t459 = rSges(4,1) * t497 + rSges(4,2) * t496 - rSges(4,3) * t540;
t457 = Icges(4,1) * t497 + Icges(4,4) * t496 - Icges(4,5) * t540;
t456 = Icges(4,4) * t497 + Icges(4,2) * t496 - Icges(4,6) * t540;
t455 = Icges(4,5) * t497 + Icges(4,6) * t496 - Icges(4,3) * t540;
t452 = rSges(3,1) * t495 - rSges(3,2) * t494 + rSges(3,3) * t544;
t451 = rSges(3,1) * t493 - rSges(3,2) * t492 - rSges(3,3) * t543;
t450 = Icges(3,1) * t495 - Icges(3,4) * t494 + Icges(3,5) * t544;
t449 = Icges(3,1) * t493 - Icges(3,4) * t492 - Icges(3,5) * t543;
t448 = Icges(3,4) * t495 - Icges(3,2) * t494 + Icges(3,6) * t544;
t447 = Icges(3,4) * t493 - Icges(3,2) * t492 - Icges(3,6) * t543;
t446 = Icges(3,5) * t495 - Icges(3,6) * t494 + Icges(3,3) * t544;
t445 = Icges(3,5) * t493 - Icges(3,6) * t492 - Icges(3,3) * t543;
t442 = rSges(5,1) * t487 - rSges(5,2) * t486 - rSges(5,3) * t540;
t441 = -rSges(6,1) * t540 - rSges(6,2) * t487 + rSges(6,3) * t486;
t434 = pkin(5) * t494 + pkin(10) * t470;
t433 = pkin(5) * t492 + pkin(10) * t468;
t432 = t469 * t511 + t494 * t514;
t431 = t469 * t514 - t494 * t511;
t430 = t467 * t511 + t492 * t514;
t429 = t467 * t514 - t492 * t511;
t425 = (-t451 * t510 - t485 * t543) * qJD(2);
t424 = (t452 * t510 - t485 * t544) * qJD(2);
t423 = qJD(6) * t468 + t444;
t422 = qJD(6) * t470 + t443;
t421 = rSges(4,1) * t476 + rSges(4,2) * t475 + rSges(4,3) * t494;
t420 = rSges(4,1) * t474 + rSges(4,2) * t473 + rSges(4,3) * t492;
t419 = Icges(4,1) * t476 + Icges(4,4) * t475 + Icges(4,5) * t494;
t418 = Icges(4,1) * t474 + Icges(4,4) * t473 + Icges(4,5) * t492;
t417 = Icges(4,4) * t476 + Icges(4,2) * t475 + Icges(4,6) * t494;
t416 = Icges(4,4) * t474 + Icges(4,2) * t473 + Icges(4,6) * t492;
t415 = Icges(4,5) * t476 + Icges(4,6) * t475 + Icges(4,3) * t494;
t414 = Icges(4,5) * t474 + Icges(4,6) * t473 + Icges(4,3) * t492;
t411 = rSges(5,1) * t470 - rSges(5,2) * t469 + rSges(5,3) * t494;
t410 = rSges(5,1) * t468 - rSges(5,2) * t467 + rSges(5,3) * t492;
t409 = rSges(6,1) * t494 - rSges(6,2) * t470 + rSges(6,3) * t469;
t408 = rSges(6,1) * t492 - rSges(6,2) * t468 + rSges(6,3) * t467;
t394 = rSges(7,1) * t472 + rSges(7,2) * t471 + rSges(7,3) * t487;
t393 = Icges(7,1) * t472 + Icges(7,4) * t471 + Icges(7,5) * t487;
t392 = Icges(7,4) * t472 + Icges(7,2) * t471 + Icges(7,6) * t487;
t391 = Icges(7,5) * t472 + Icges(7,6) * t471 + Icges(7,3) * t487;
t388 = qJD(1) + (t451 * t507 + t452 * t509) * t535;
t385 = rSges(7,1) * t432 + rSges(7,2) * t431 + rSges(7,3) * t470;
t384 = rSges(7,1) * t430 + rSges(7,2) * t429 + rSges(7,3) * t468;
t383 = Icges(7,1) * t432 + Icges(7,4) * t431 + Icges(7,5) * t470;
t382 = Icges(7,1) * t430 + Icges(7,4) * t429 + Icges(7,5) * t468;
t381 = Icges(7,4) * t432 + Icges(7,2) * t431 + Icges(7,6) * t470;
t380 = Icges(7,4) * t430 + Icges(7,2) * t429 + Icges(7,6) * t468;
t379 = Icges(7,5) * t432 + Icges(7,6) * t431 + Icges(7,3) * t470;
t378 = Icges(7,5) * t430 + Icges(7,6) * t429 + Icges(7,3) * t468;
t377 = -t420 * t499 + t459 * t479 + t523;
t376 = t421 * t499 - t459 * t478 + t525;
t375 = t420 * t478 - t421 * t479 + t530;
t374 = -t410 * t481 + t442 * t444 + t519;
t373 = t411 * t481 - t442 * t443 + t521;
t372 = t410 * t443 - t411 * t444 + t524;
t371 = t441 * t444 + (-t408 - t426) * t481 + t518;
t370 = t409 * t481 + (-t441 - t458) * t443 + t520;
t369 = t408 * t443 + (-t409 - t427) * t444 + t522;
t368 = -t384 * t461 + t394 * t423 + t444 * t477 + (-t426 - t433) * t481 + t518;
t367 = t385 * t461 - t394 * t422 + t434 * t481 + (-t458 - t477) * t443 + t520;
t366 = t384 * t422 - t385 * t423 + t433 * t443 + (-t427 - t434) * t444 + t522;
t1 = m(5) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(6) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + m(7) * (t366 ^ 2 + t367 ^ 2 + t368 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t388 ^ 2 + t424 ^ 2 + t425 ^ 2) / 0.2e1 + t478 * ((t415 * t494 + t417 * t475 + t419 * t476) * t478 + (t414 * t494 + t416 * t475 + t418 * t476) * t479 + (t455 * t494 + t456 * t475 + t457 * t476) * t499) / 0.2e1 + t479 * ((t415 * t492 + t417 * t473 + t419 * t474) * t478 + (t414 * t492 + t416 * t473 + t418 * t474) * t479 + (t455 * t492 + t456 * t473 + t457 * t474) * t499) / 0.2e1 + t499 * ((-t415 * t540 + t417 * t496 + t419 * t497) * t478 + (-t414 * t540 + t416 * t496 + t418 * t497) * t479 + (-t455 * t540 + t456 * t496 + t457 * t497) * t499) / 0.2e1 + t422 * ((t470 * t379 + t431 * t381 + t432 * t383) * t422 + (t378 * t470 + t380 * t431 + t382 * t432) * t423 + (t391 * t470 + t392 * t431 + t393 * t432) * t461) / 0.2e1 + t423 * ((t379 * t468 + t381 * t429 + t383 * t430) * t422 + (t468 * t378 + t429 * t380 + t430 * t382) * t423 + (t391 * t468 + t392 * t429 + t393 * t430) * t461) / 0.2e1 + t461 * ((t379 * t487 + t381 * t471 + t383 * t472) * t422 + (t378 * t487 + t380 * t471 + t382 * t472) * t423 + (t391 * t487 + t392 * t471 + t393 * t472) * t461) / 0.2e1 + m(4) * (t375 ^ 2 + t376 ^ 2 + t377 ^ 2) / 0.2e1 - t550 * ((-t446 * t543 - t448 * t492 + t450 * t493) * t544 - (-t445 * t543 - t447 * t492 + t449 * t493) * t543 + (-t482 * t543 - t483 * t492 + t484 * t493) * t510) * t543 / 0.2e1 + ((t469 * t553 + t470 * t552 + t494 * t551) * t481 + (t469 * t559 + t555 * t470 + t557 * t494) * t444 + (t469 * t558 + t470 * t554 + t494 * t556) * t443) * t443 / 0.2e1 + ((t467 * t553 + t468 * t552 + t492 * t551) * t481 + (t467 * t559 + t555 * t468 + t557 * t492) * t444 + (t467 * t558 + t468 * t554 + t492 * t556) * t443) * t444 / 0.2e1 + ((t486 * t553 + t487 * t552 - t540 * t551) * t481 + (t486 * t559 + t555 * t487 - t557 * t540) * t444 + (t486 * t558 + t487 * t554 - t540 * t556) * t443) * t481 / 0.2e1 + (t510 * (t510 ^ 2 * t482 + (((t448 * t516 + t450 * t513) * t507 - (t447 * t516 + t449 * t513) * t509) * t508 + (-t445 * t509 + t446 * t507 + t483 * t516 + t484 * t513) * t510) * t508) + ((t446 * t544 - t448 * t494 + t450 * t495) * t544 - (t445 * t544 - t447 * t494 + t449 * t495) * t543 + (t482 * t544 - t483 * t494 + t484 * t495) * t510) * t544) * t550 / 0.2e1;
T  = t1;
