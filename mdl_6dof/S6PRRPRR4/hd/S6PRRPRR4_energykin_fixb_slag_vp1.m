% Calculate kinetic energy for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:10:16
% EndTime: 2019-03-08 22:10:19
% DurationCPUTime: 2.92s
% Computational Cost: add. (2632->330), mult. (6747->506), div. (0->0), fcn. (8389->12), ass. (0->148)
t568 = Icges(4,1) + Icges(5,1);
t567 = -Icges(4,4) + Icges(5,5);
t566 = Icges(5,4) + Icges(4,5);
t565 = Icges(4,2) + Icges(5,3);
t564 = Icges(5,2) + Icges(4,3);
t563 = Icges(4,6) - Icges(5,6);
t520 = sin(pkin(11));
t522 = cos(pkin(11));
t528 = cos(qJ(2));
t526 = sin(qJ(2));
t548 = cos(pkin(6));
t538 = t526 * t548;
t507 = t520 * t528 + t522 * t538;
t525 = sin(qJ(3));
t521 = sin(pkin(6));
t550 = cos(qJ(3));
t541 = t521 * t550;
t488 = t507 * t525 + t522 * t541;
t544 = t522 * t521;
t489 = t507 * t550 - t525 * t544;
t537 = t528 * t548;
t506 = t520 * t526 - t522 * t537;
t562 = t488 * t565 + t489 * t567 - t506 * t563;
t509 = -t520 * t538 + t522 * t528;
t490 = t509 * t525 - t520 * t541;
t546 = t521 * t525;
t491 = t509 * t550 + t520 * t546;
t508 = t520 * t537 + t522 * t526;
t561 = t490 * t565 + t491 * t567 - t508 * t563;
t560 = -t488 * t563 + t489 * t566 + t506 * t564;
t559 = -t490 * t563 + t491 * t566 + t508 * t564;
t558 = t567 * t488 + t489 * t568 + t566 * t506;
t557 = t567 * t490 + t491 * t568 + t566 * t508;
t510 = t526 * t546 - t548 * t550;
t511 = t525 * t548 + t526 * t541;
t545 = t521 * t528;
t556 = t510 * t565 + t511 * t567 + t545 * t563;
t555 = -t510 * t563 + t511 * t566 - t545 * t564;
t554 = t567 * t510 + t511 * t568 - t566 * t545;
t549 = cos(qJ(5));
t547 = t520 * t521;
t543 = qJD(2) * t521;
t516 = t520 * t543;
t492 = qJD(3) * t508 + t516;
t458 = -qJD(5) * t508 + t492;
t540 = t522 * t543;
t481 = pkin(2) * t507 + pkin(8) * t506;
t482 = pkin(2) * t509 + pkin(8) * t508;
t539 = t481 * t516 + t482 * t540 + qJD(1);
t519 = qJD(2) * t548;
t493 = qJD(3) * t506 - t540;
t513 = -qJD(3) * t545 + t519;
t450 = pkin(3) * t489 + qJ(4) * t488;
t536 = qJD(4) * t510 + t492 * t450 + t539;
t512 = (pkin(2) * t526 - pkin(8) * t528) * t521;
t535 = t482 * t519 - t512 * t516;
t459 = -qJD(5) * t506 + t493;
t495 = qJD(5) * t545 + t513;
t451 = pkin(3) * t491 + qJ(4) * t490;
t534 = qJD(4) * t488 + t513 * t451 + t535;
t533 = (-t481 * t548 - t512 * t544) * qJD(2);
t456 = pkin(4) * t489 - pkin(9) * t506;
t457 = pkin(4) * t491 - pkin(9) * t508;
t532 = t492 * t456 + (-t451 - t457) * t493 + t536;
t483 = pkin(3) * t511 + qJ(4) * t510;
t531 = qJD(4) * t490 + t493 * t483 + t533;
t494 = t511 * pkin(4) + pkin(9) * t545;
t530 = t513 * t457 + (-t483 - t494) * t492 + t534;
t529 = t493 * t494 + (-t450 - t456) * t513 + t531;
t527 = cos(qJ(6));
t524 = sin(qJ(5));
t523 = sin(qJ(6));
t499 = t548 * rSges(3,3) + (rSges(3,1) * t526 + rSges(3,2) * t528) * t521;
t498 = Icges(3,5) * t548 + (Icges(3,1) * t526 + Icges(3,4) * t528) * t521;
t497 = Icges(3,6) * t548 + (Icges(3,4) * t526 + Icges(3,2) * t528) * t521;
t496 = Icges(3,3) * t548 + (Icges(3,5) * t526 + Icges(3,6) * t528) * t521;
t480 = t510 * t524 + t511 * t549;
t479 = -t510 * t549 + t511 * t524;
t477 = t511 * rSges(4,1) - t510 * rSges(4,2) - rSges(4,3) * t545;
t476 = t511 * rSges(5,1) - rSges(5,2) * t545 + t510 * rSges(5,3);
t467 = rSges(3,1) * t509 - rSges(3,2) * t508 + rSges(3,3) * t547;
t466 = rSges(3,1) * t507 - rSges(3,2) * t506 - rSges(3,3) * t544;
t465 = Icges(3,1) * t509 - Icges(3,4) * t508 + Icges(3,5) * t547;
t464 = Icges(3,1) * t507 - Icges(3,4) * t506 - Icges(3,5) * t544;
t463 = Icges(3,4) * t509 - Icges(3,2) * t508 + Icges(3,6) * t547;
t462 = Icges(3,4) * t507 - Icges(3,2) * t506 - Icges(3,6) * t544;
t461 = Icges(3,5) * t509 - Icges(3,6) * t508 + Icges(3,3) * t547;
t460 = Icges(3,5) * t507 - Icges(3,6) * t506 - Icges(3,3) * t544;
t455 = t480 * t527 + t523 * t545;
t454 = -t480 * t523 + t527 * t545;
t452 = qJD(6) * t479 + t495;
t447 = t490 * t524 + t491 * t549;
t446 = -t490 * t549 + t491 * t524;
t445 = t488 * t524 + t489 * t549;
t444 = -t488 * t549 + t489 * t524;
t441 = (-t466 * t548 - t499 * t544) * qJD(2);
t440 = (t467 * t548 - t499 * t547) * qJD(2);
t439 = rSges(4,1) * t491 - rSges(4,2) * t490 + rSges(4,3) * t508;
t438 = rSges(5,1) * t491 + rSges(5,2) * t508 + rSges(5,3) * t490;
t437 = rSges(4,1) * t489 - rSges(4,2) * t488 + rSges(4,3) * t506;
t436 = rSges(5,1) * t489 + rSges(5,2) * t506 + rSges(5,3) * t488;
t423 = pkin(5) * t480 + pkin(10) * t479;
t422 = t480 * rSges(6,1) - t479 * rSges(6,2) + rSges(6,3) * t545;
t421 = t447 * t527 - t508 * t523;
t420 = -t447 * t523 - t508 * t527;
t419 = t445 * t527 - t506 * t523;
t418 = -t445 * t523 - t506 * t527;
t417 = Icges(6,1) * t480 - Icges(6,4) * t479 + Icges(6,5) * t545;
t416 = Icges(6,4) * t480 - Icges(6,2) * t479 + Icges(6,6) * t545;
t415 = Icges(6,5) * t480 - Icges(6,6) * t479 + Icges(6,3) * t545;
t413 = qJD(1) + (t466 * t520 + t467 * t522) * t543;
t412 = qJD(6) * t444 + t459;
t411 = qJD(6) * t446 + t458;
t410 = pkin(5) * t447 + pkin(10) * t446;
t409 = pkin(5) * t445 + pkin(10) * t444;
t408 = rSges(7,1) * t455 + rSges(7,2) * t454 + rSges(7,3) * t479;
t407 = Icges(7,1) * t455 + Icges(7,4) * t454 + Icges(7,5) * t479;
t406 = Icges(7,4) * t455 + Icges(7,2) * t454 + Icges(7,6) * t479;
t405 = Icges(7,5) * t455 + Icges(7,6) * t454 + Icges(7,3) * t479;
t404 = rSges(6,1) * t447 - rSges(6,2) * t446 - rSges(6,3) * t508;
t403 = rSges(6,1) * t445 - rSges(6,2) * t444 - rSges(6,3) * t506;
t402 = Icges(6,1) * t447 - Icges(6,4) * t446 - Icges(6,5) * t508;
t401 = Icges(6,1) * t445 - Icges(6,4) * t444 - Icges(6,5) * t506;
t400 = Icges(6,4) * t447 - Icges(6,2) * t446 - Icges(6,6) * t508;
t399 = Icges(6,4) * t445 - Icges(6,2) * t444 - Icges(6,6) * t506;
t398 = Icges(6,5) * t447 - Icges(6,6) * t446 - Icges(6,3) * t508;
t397 = Icges(6,5) * t445 - Icges(6,6) * t444 - Icges(6,3) * t506;
t396 = -t513 * t437 + t493 * t477 + t533;
t395 = t439 * t513 - t477 * t492 + t535;
t394 = rSges(7,1) * t421 + rSges(7,2) * t420 + rSges(7,3) * t446;
t393 = rSges(7,1) * t419 + rSges(7,2) * t418 + rSges(7,3) * t444;
t392 = Icges(7,1) * t421 + Icges(7,4) * t420 + Icges(7,5) * t446;
t391 = Icges(7,1) * t419 + Icges(7,4) * t418 + Icges(7,5) * t444;
t390 = Icges(7,4) * t421 + Icges(7,2) * t420 + Icges(7,6) * t446;
t389 = Icges(7,4) * t419 + Icges(7,2) * t418 + Icges(7,6) * t444;
t388 = Icges(7,5) * t421 + Icges(7,6) * t420 + Icges(7,3) * t446;
t387 = Icges(7,5) * t419 + Icges(7,6) * t418 + Icges(7,3) * t444;
t386 = t437 * t492 - t439 * t493 + t539;
t385 = t493 * t476 + (-t436 - t450) * t513 + t531;
t384 = t438 * t513 + (-t476 - t483) * t492 + t534;
t383 = t436 * t492 + (-t438 - t451) * t493 + t536;
t382 = -t495 * t403 + t459 * t422 + t529;
t381 = t404 * t495 - t422 * t458 + t530;
t380 = t403 * t458 - t404 * t459 + t532;
t379 = -t452 * t393 + t412 * t408 - t495 * t409 + t459 * t423 + t529;
t378 = t394 * t452 - t408 * t411 + t410 * t495 - t423 * t458 + t530;
t377 = t393 * t411 - t394 * t412 + t409 * t458 - t410 * t459 + t532;
t1 = m(6) * (t380 ^ 2 + t381 ^ 2 + t382 ^ 2) / 0.2e1 + m(7) * (t377 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t413 ^ 2 + t440 ^ 2 + t441 ^ 2) / 0.2e1 + m(4) * (t386 ^ 2 + t395 ^ 2 + t396 ^ 2) / 0.2e1 + m(5) * (t383 ^ 2 + t384 ^ 2 + t385 ^ 2) / 0.2e1 + ((t548 * t461 + (t463 * t528 + t465 * t526) * t521) * t516 - (t548 * t460 + (t462 * t528 + t464 * t526) * t521) * t540 + (t548 * t496 + (t497 * t528 + t498 * t526) * t521) * t519) * t519 / 0.2e1 + t458 * ((-t398 * t508 - t400 * t446 + t402 * t447) * t458 + (-t397 * t508 - t399 * t446 + t401 * t447) * t459 + (-t415 * t508 - t416 * t446 + t417 * t447) * t495) / 0.2e1 + t459 * ((-t398 * t506 - t400 * t444 + t402 * t445) * t458 + (-t397 * t506 - t399 * t444 + t401 * t445) * t459 + (-t415 * t506 - t416 * t444 + t417 * t445) * t495) / 0.2e1 + t495 * ((t398 * t545 - t479 * t400 + t480 * t402) * t458 + (t397 * t545 - t479 * t399 + t480 * t401) * t459 + (t415 * t545 - t479 * t416 + t480 * t417) * t495) / 0.2e1 + t411 * ((t388 * t446 + t390 * t420 + t392 * t421) * t411 + (t387 * t446 + t389 * t420 + t391 * t421) * t412 + (t405 * t446 + t406 * t420 + t407 * t421) * t452) / 0.2e1 + t412 * ((t388 * t444 + t390 * t418 + t392 * t419) * t411 + (t387 * t444 + t389 * t418 + t391 * t419) * t412 + (t405 * t444 + t406 * t418 + t407 * t419) * t452) / 0.2e1 + t452 * ((t479 * t388 + t454 * t390 + t455 * t392) * t411 + (t479 * t387 + t454 * t389 + t455 * t391) * t412 + (t479 * t405 + t454 * t406 + t455 * t407) * t452) / 0.2e1 + ((t490 * t556 + t491 * t554 + t508 * t555) * t513 + (t490 * t562 + t558 * t491 + t560 * t508) * t493 + (t561 * t490 + t557 * t491 + t559 * t508) * t492) * t492 / 0.2e1 + ((t488 * t556 + t489 * t554 + t506 * t555) * t513 + (t562 * t488 + t558 * t489 + t560 * t506) * t493 + (t488 * t561 + t489 * t557 + t506 * t559) * t492) * t493 / 0.2e1 + ((t556 * t510 + t554 * t511 - t555 * t545) * t513 + (t510 * t562 + t558 * t511 - t560 * t545) * t493 + (t510 * t561 + t511 * t557 - t545 * t559) * t492) * t513 / 0.2e1 + (-t522 * ((-t461 * t544 - t463 * t506 + t465 * t507) * t547 - (-t460 * t544 - t462 * t506 + t464 * t507) * t544 + (-t496 * t544 - t497 * t506 + t498 * t507) * t548) / 0.2e1 + t520 * ((t461 * t547 - t463 * t508 + t465 * t509) * t547 - (t460 * t547 - t462 * t508 + t464 * t509) * t544 + (t496 * t547 - t497 * t508 + t498 * t509) * t548) / 0.2e1) * t521 * qJD(2) ^ 2;
T  = t1;
