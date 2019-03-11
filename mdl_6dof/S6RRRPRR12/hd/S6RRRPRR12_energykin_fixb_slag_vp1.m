% Calculate kinetic energy for
% S6RRRPRR12
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
% Datum: 2019-03-09 19:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR12_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR12_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR12_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR12_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR12_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR12_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:38:09
% EndTime: 2019-03-09 19:38:12
% DurationCPUTime: 3.22s
% Computational Cost: add. (3197->383), mult. (6571->577), div. (0->0), fcn. (8049->14), ass. (0->169)
t576 = Icges(4,2) + Icges(5,3);
t533 = sin(qJ(2));
t534 = sin(qJ(1));
t535 = cos(qJ(2));
t536 = cos(qJ(1));
t564 = cos(pkin(6));
t547 = t536 * t564;
t501 = t533 * t534 - t535 * t547;
t502 = t533 * t547 + t534 * t535;
t529 = sin(pkin(6));
t559 = t529 * t536;
t460 = Icges(3,5) * t502 - Icges(3,6) * t501 - Icges(3,3) * t559;
t548 = t534 * t564;
t503 = t536 * t533 + t535 * t548;
t504 = -t533 * t548 + t536 * t535;
t561 = t529 * t534;
t461 = Icges(3,5) * t504 - Icges(3,6) * t503 + Icges(3,3) * t561;
t575 = t529 * (t460 * t536 - t461 * t534);
t532 = sin(qJ(3));
t567 = cos(qJ(3));
t486 = t502 * t567 - t532 * t559;
t528 = sin(pkin(12));
t530 = cos(pkin(12));
t451 = -t486 * t528 + t501 * t530;
t563 = t501 * t528;
t452 = t486 * t530 + t563;
t551 = t529 * t567;
t485 = t502 * t532 + t536 * t551;
t574 = -Icges(4,4) * t486 + Icges(5,5) * t452 - Icges(4,6) * t501 + Icges(5,6) * t451 + t485 * t576;
t488 = t504 * t567 + t532 * t561;
t453 = -t488 * t528 + t503 * t530;
t562 = t503 * t528;
t454 = t488 * t530 + t562;
t487 = t504 * t532 - t534 * t551;
t573 = -Icges(4,4) * t488 + Icges(5,5) * t454 - Icges(4,6) * t503 + Icges(5,6) * t453 + t487 * t576;
t500 = t532 * t564 + t533 * t551;
t560 = t529 * t535;
t483 = -t500 * t528 - t530 * t560;
t552 = t528 * t560;
t484 = t500 * t530 - t552;
t499 = t529 * t532 * t533 - t564 * t567;
t572 = -Icges(4,4) * t500 + Icges(5,5) * t484 + Icges(4,6) * t560 + Icges(5,6) * t483 + t499 * t576;
t565 = t530 * pkin(4);
t474 = pkin(2) * t502 + pkin(9) * t501;
t475 = pkin(2) * t504 + pkin(9) * t503;
t553 = qJD(2) * t529;
t515 = t534 * t553;
t550 = t536 * t553;
t557 = t474 * t515 + t475 * t550;
t527 = pkin(12) + qJ(5);
t523 = cos(t527);
t556 = pkin(5) * t523;
t489 = qJD(3) * t503 + t515;
t554 = qJD(1) * (pkin(1) * t534 - pkin(8) * t559);
t516 = qJD(2) * t564 + qJD(1);
t445 = qJD(5) * t487 + t489;
t522 = sin(t527);
t549 = pkin(5) * t522;
t439 = pkin(3) * t486 + qJ(4) * t485;
t546 = qJD(4) * t499 + t489 * t439 + t557;
t490 = qJD(3) * t501 - t550;
t446 = qJD(5) * t485 + t490;
t506 = -qJD(3) * t560 + t516;
t505 = (pkin(2) * t533 - pkin(9) * t535) * t529;
t507 = qJD(1) * (pkin(1) * t536 + pkin(8) * t561);
t544 = t516 * t475 - t505 * t515 + t507;
t478 = qJD(5) * t499 + t506;
t440 = pkin(3) * t488 + qJ(4) * t487;
t543 = qJD(4) * t485 + t506 * t440 + t544;
t390 = pkin(4) * t563 + pkin(10) * t485 + t486 * t565;
t391 = pkin(4) * t562 + pkin(10) * t487 + t488 * t565;
t542 = t489 * t390 + (-t391 - t440) * t490 + t546;
t541 = -t474 * t516 - t505 * t550 - t554;
t471 = pkin(3) * t500 + qJ(4) * t499;
t540 = qJD(4) * t487 + t490 * t471 + t541;
t421 = -pkin(4) * t552 + pkin(10) * t499 + t500 * t565;
t539 = t506 * t391 + (-t421 - t471) * t489 + t543;
t538 = t490 * t421 + (-t390 - t439) * t506 + t540;
t524 = qJ(6) + t527;
t519 = cos(t524);
t518 = sin(t524);
t512 = rSges(2,1) * t536 - rSges(2,2) * t534;
t511 = rSges(2,1) * t534 + rSges(2,2) * t536;
t494 = t564 * rSges(3,3) + (rSges(3,1) * t533 + rSges(3,2) * t535) * t529;
t493 = Icges(3,5) * t564 + (Icges(3,1) * t533 + Icges(3,4) * t535) * t529;
t492 = Icges(3,6) * t564 + (Icges(3,4) * t533 + Icges(3,2) * t535) * t529;
t491 = Icges(3,3) * t564 + (Icges(3,5) * t533 + Icges(3,6) * t535) * t529;
t477 = t500 * t523 - t522 * t560;
t476 = -t500 * t522 - t523 * t560;
t473 = t500 * t519 - t518 * t560;
t472 = -t500 * t518 - t519 * t560;
t468 = rSges(3,1) * t504 - rSges(3,2) * t503 + rSges(3,3) * t561;
t467 = rSges(3,1) * t502 - rSges(3,2) * t501 - rSges(3,3) * t559;
t465 = Icges(3,1) * t504 - Icges(3,4) * t503 + Icges(3,5) * t561;
t464 = Icges(3,1) * t502 - Icges(3,4) * t501 - Icges(3,5) * t559;
t463 = Icges(3,4) * t504 - Icges(3,2) * t503 + Icges(3,6) * t561;
t462 = Icges(3,4) * t502 - Icges(3,2) * t501 - Icges(3,6) * t559;
t459 = rSges(4,1) * t500 - rSges(4,2) * t499 - rSges(4,3) * t560;
t458 = Icges(4,1) * t500 - Icges(4,4) * t499 - Icges(4,5) * t560;
t456 = Icges(4,5) * t500 - Icges(4,6) * t499 - Icges(4,3) * t560;
t455 = qJD(6) * t499 + t478;
t450 = t488 * t523 + t503 * t522;
t449 = -t488 * t522 + t503 * t523;
t448 = t486 * t523 + t501 * t522;
t447 = -t486 * t522 + t501 * t523;
t444 = t488 * t519 + t503 * t518;
t443 = -t488 * t518 + t503 * t519;
t442 = t486 * t519 + t501 * t518;
t441 = -t486 * t518 + t501 * t519;
t436 = rSges(4,1) * t488 - rSges(4,2) * t487 + rSges(4,3) * t503;
t435 = rSges(4,1) * t486 - rSges(4,2) * t485 + rSges(4,3) * t501;
t434 = Icges(4,1) * t488 - Icges(4,4) * t487 + Icges(4,5) * t503;
t433 = Icges(4,1) * t486 - Icges(4,4) * t485 + Icges(4,5) * t501;
t430 = Icges(4,5) * t488 - Icges(4,6) * t487 + Icges(4,3) * t503;
t429 = Icges(4,5) * t486 - Icges(4,6) * t485 + Icges(4,3) * t501;
t428 = rSges(5,1) * t484 + rSges(5,2) * t483 + rSges(5,3) * t499;
t427 = Icges(5,1) * t484 + Icges(5,4) * t483 + Icges(5,5) * t499;
t426 = Icges(5,4) * t484 + Icges(5,2) * t483 + Icges(5,6) * t499;
t424 = qJD(6) * t485 + t446;
t423 = qJD(6) * t487 + t445;
t420 = rSges(6,1) * t477 + rSges(6,2) * t476 + rSges(6,3) * t499;
t419 = Icges(6,1) * t477 + Icges(6,4) * t476 + Icges(6,5) * t499;
t418 = Icges(6,4) * t477 + Icges(6,2) * t476 + Icges(6,6) * t499;
t417 = Icges(6,5) * t477 + Icges(6,6) * t476 + Icges(6,3) * t499;
t416 = t468 * t516 - t494 * t515 + t507;
t415 = -t467 * t516 - t494 * t550 - t554;
t414 = rSges(7,1) * t473 + rSges(7,2) * t472 + rSges(7,3) * t499;
t413 = Icges(7,1) * t473 + Icges(7,4) * t472 + Icges(7,5) * t499;
t412 = Icges(7,4) * t473 + Icges(7,2) * t472 + Icges(7,6) * t499;
t411 = Icges(7,5) * t473 + Icges(7,6) * t472 + Icges(7,3) * t499;
t410 = (t467 * t534 + t468 * t536) * t553;
t408 = pkin(11) * t499 + t500 * t556 - t549 * t560;
t407 = rSges(5,1) * t454 + rSges(5,2) * t453 + rSges(5,3) * t487;
t406 = rSges(5,1) * t452 + rSges(5,2) * t451 + rSges(5,3) * t485;
t405 = Icges(5,1) * t454 + Icges(5,4) * t453 + Icges(5,5) * t487;
t404 = Icges(5,1) * t452 + Icges(5,4) * t451 + Icges(5,5) * t485;
t403 = Icges(5,4) * t454 + Icges(5,2) * t453 + Icges(5,6) * t487;
t402 = Icges(5,4) * t452 + Icges(5,2) * t451 + Icges(5,6) * t485;
t399 = rSges(6,1) * t450 + rSges(6,2) * t449 + rSges(6,3) * t487;
t398 = rSges(6,1) * t448 + rSges(6,2) * t447 + rSges(6,3) * t485;
t397 = Icges(6,1) * t450 + Icges(6,4) * t449 + Icges(6,5) * t487;
t396 = Icges(6,1) * t448 + Icges(6,4) * t447 + Icges(6,5) * t485;
t395 = Icges(6,4) * t450 + Icges(6,2) * t449 + Icges(6,6) * t487;
t394 = Icges(6,4) * t448 + Icges(6,2) * t447 + Icges(6,6) * t485;
t393 = Icges(6,5) * t450 + Icges(6,6) * t449 + Icges(6,3) * t487;
t392 = Icges(6,5) * t448 + Icges(6,6) * t447 + Icges(6,3) * t485;
t389 = rSges(7,1) * t444 + rSges(7,2) * t443 + rSges(7,3) * t487;
t388 = rSges(7,1) * t442 + rSges(7,2) * t441 + rSges(7,3) * t485;
t387 = Icges(7,1) * t444 + Icges(7,4) * t443 + Icges(7,5) * t487;
t386 = Icges(7,1) * t442 + Icges(7,4) * t441 + Icges(7,5) * t485;
t385 = Icges(7,4) * t444 + Icges(7,2) * t443 + Icges(7,6) * t487;
t384 = Icges(7,4) * t442 + Icges(7,2) * t441 + Icges(7,6) * t485;
t383 = Icges(7,5) * t444 + Icges(7,6) * t443 + Icges(7,3) * t487;
t382 = Icges(7,5) * t442 + Icges(7,6) * t441 + Icges(7,3) * t485;
t379 = pkin(11) * t487 + t488 * t556 + t503 * t549;
t378 = pkin(11) * t485 + t486 * t556 + t501 * t549;
t377 = t436 * t506 - t459 * t489 + t544;
t376 = -t435 * t506 + t459 * t490 + t541;
t375 = t435 * t489 - t436 * t490 + t557;
t374 = t407 * t506 + (-t428 - t471) * t489 + t543;
t373 = t428 * t490 + (-t406 - t439) * t506 + t540;
t372 = t406 * t489 + (-t407 - t440) * t490 + t546;
t371 = t399 * t478 - t420 * t445 + t539;
t370 = -t398 * t478 + t420 * t446 + t538;
t369 = t398 * t445 - t399 * t446 + t542;
t368 = t379 * t478 + t389 * t455 - t408 * t445 - t414 * t423 + t539;
t367 = -t378 * t478 - t388 * t455 + t408 * t446 + t414 * t424 + t538;
t366 = t378 * t445 - t379 * t446 + t388 * t423 - t389 * t424 + t542;
t1 = m(3) * (t410 ^ 2 + t415 ^ 2 + t416 ^ 2) / 0.2e1 + m(4) * (t375 ^ 2 + t376 ^ 2 + t377 ^ 2) / 0.2e1 + m(5) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(6) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + m(7) * (t366 ^ 2 + t367 ^ 2 + t368 ^ 2) / 0.2e1 + t516 * ((t564 * t461 + (t463 * t535 + t465 * t533) * t529) * t515 - (t564 * t460 + (t462 * t535 + t464 * t533) * t529) * t550 + (t564 * t491 + (t492 * t535 + t493 * t533) * t529) * t516) / 0.2e1 + t445 * ((t393 * t487 + t395 * t449 + t397 * t450) * t445 + (t392 * t487 + t394 * t449 + t396 * t450) * t446 + (t417 * t487 + t418 * t449 + t419 * t450) * t478) / 0.2e1 + t446 * ((t393 * t485 + t395 * t447 + t397 * t448) * t445 + (t392 * t485 + t394 * t447 + t396 * t448) * t446 + (t417 * t485 + t418 * t447 + t419 * t448) * t478) / 0.2e1 + t478 * ((t393 * t499 + t395 * t476 + t397 * t477) * t445 + (t392 * t499 + t394 * t476 + t396 * t477) * t446 + (t417 * t499 + t418 * t476 + t419 * t477) * t478) / 0.2e1 + t423 * ((t383 * t487 + t385 * t443 + t387 * t444) * t423 + (t382 * t487 + t384 * t443 + t386 * t444) * t424 + (t411 * t487 + t412 * t443 + t413 * t444) * t455) / 0.2e1 + t424 * ((t383 * t485 + t385 * t441 + t387 * t442) * t423 + (t382 * t485 + t384 * t441 + t386 * t442) * t424 + (t411 * t485 + t412 * t441 + t413 * t442) * t455) / 0.2e1 + t455 * ((t383 * t499 + t385 * t472 + t387 * t473) * t423 + (t382 * t499 + t384 * t472 + t386 * t473) * t424 + (t411 * t499 + t412 * t472 + t413 * t473) * t455) / 0.2e1 - ((-t491 * t559 - t492 * t501 + t493 * t502) * t516 + ((-t463 * t501 + t465 * t502) * t534 + (t501 * t462 - t502 * t464 + t575) * t536) * t553) * t550 / 0.2e1 + ((t491 * t561 - t492 * t503 + t493 * t504) * t516 + (-(-t462 * t503 + t464 * t504) * t536 + (-t503 * t463 + t504 * t465 - t575) * t534) * t553) * t515 / 0.2e1 + ((t426 * t453 + t427 * t454 + t456 * t503 + t458 * t488 + t487 * t572) * t506 + (t402 * t453 + t404 * t454 + t429 * t503 + t433 * t488 + t487 * t574) * t490 + (t403 * t453 + t405 * t454 + t430 * t503 + t434 * t488 + t487 * t573) * t489) * t489 / 0.2e1 + ((t426 * t451 + t427 * t452 + t456 * t501 + t458 * t486 + t485 * t572) * t506 + (t402 * t451 + t404 * t452 + t429 * t501 + t433 * t486 + t574 * t485) * t490 + (t403 * t451 + t405 * t452 + t430 * t501 + t434 * t486 + t485 * t573) * t489) * t490 / 0.2e1 + ((t426 * t483 + t427 * t484 - t456 * t560 + t458 * t500 + t499 * t572) * t506 + (t402 * t483 + t404 * t484 - t429 * t560 + t433 * t500 + t499 * t574) * t490 + (t403 * t483 + t405 * t484 - t430 * t560 + t434 * t500 + t499 * t573) * t489) * t506 / 0.2e1 + (m(2) * (t511 ^ 2 + t512 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
