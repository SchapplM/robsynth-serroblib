% Calculate kinetic energy for
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:11:39
% EndTime: 2019-03-10 04:11:42
% DurationCPUTime: 2.55s
% Computational Cost: add. (3581->389), mult. (6119->608), div. (0->0), fcn. (7384->14), ass. (0->178)
t528 = sin(qJ(2));
t529 = sin(qJ(1));
t532 = cos(qJ(2));
t533 = cos(qJ(1));
t566 = cos(pkin(6));
t548 = t533 * t566;
t503 = t528 * t529 - t532 * t548;
t504 = t528 * t548 + t529 * t532;
t525 = sin(pkin(6));
t559 = t525 * t533;
t459 = Icges(3,5) * t504 - Icges(3,6) * t503 - Icges(3,3) * t559;
t549 = t529 * t566;
t505 = t533 * t528 + t532 * t549;
t506 = -t528 * t549 + t533 * t532;
t562 = t525 * t529;
t460 = Icges(3,5) * t506 - Icges(3,6) * t505 + Icges(3,3) * t562;
t572 = t525 * (t459 * t533 - t460 * t529);
t531 = cos(qJ(3));
t570 = pkin(3) * t531;
t530 = cos(qJ(5));
t569 = pkin(5) * t530;
t526 = sin(qJ(5));
t565 = t503 * t526;
t564 = t505 * t526;
t563 = t525 * t528;
t561 = t525 * t531;
t560 = t525 * t532;
t558 = qJ(3) + qJ(4);
t473 = pkin(2) * t504 + pkin(9) * t503;
t474 = pkin(2) * t506 + pkin(9) * t505;
t555 = qJD(2) * t525;
t516 = t529 * t555;
t551 = t533 * t555;
t557 = t473 * t516 + t474 * t551;
t487 = qJD(3) * t505 + t516;
t556 = qJD(1) * (pkin(1) * t529 - pkin(8) * t559);
t517 = qJD(2) * t566 + qJD(1);
t554 = t526 * t560;
t527 = sin(qJ(3));
t553 = t527 * t562;
t552 = t527 * t559;
t455 = qJD(4) * t505 + t487;
t550 = cos(t558);
t547 = t566 * t527;
t522 = sin(t558);
t546 = t525 * t550;
t481 = t506 * t522 - t529 * t546;
t433 = qJD(5) * t481 + t455;
t488 = qJD(3) * t503 - t551;
t423 = -pkin(3) * t552 + pkin(10) * t503 + t504 * t570;
t424 = pkin(3) * t553 + pkin(10) * t505 + t506 * t570;
t544 = t487 * t423 - t424 * t488 + t557;
t456 = qJD(4) * t503 + t488;
t507 = (pkin(2) * t528 - pkin(9) * t532) * t525;
t509 = qJD(1) * (pkin(1) * t533 + pkin(8) * t562);
t543 = t517 * t474 - t507 * t516 + t509;
t479 = t504 * t522 + t533 * t546;
t434 = qJD(5) * t479 + t456;
t489 = (-qJD(3) - qJD(4)) * t560 + t517;
t495 = t522 * t563 - t550 * t566;
t451 = qJD(5) * t495 + t489;
t480 = t504 * t550 - t522 * t559;
t436 = pkin(4) * t480 + pkin(11) * t479;
t482 = t506 * t550 + t522 * t562;
t437 = pkin(4) * t482 + pkin(11) * t481;
t542 = t455 * t436 - t437 * t456 + t544;
t541 = -t473 * t517 - t507 * t551 - t556;
t465 = pkin(3) * t547 + (-pkin(10) * t532 + t528 * t570) * t525;
t508 = -qJD(3) * t560 + t517;
t540 = t508 * t424 - t465 * t487 + t543;
t539 = -t423 * t508 + t488 * t465 + t541;
t496 = t522 * t566 + t528 * t546;
t458 = pkin(4) * t496 + pkin(11) * t495;
t538 = t489 * t437 - t455 * t458 + t540;
t537 = -t436 * t489 + t456 * t458 + t539;
t524 = qJ(5) + qJ(6);
t523 = cos(t524);
t521 = sin(t524);
t513 = rSges(2,1) * t533 - rSges(2,2) * t529;
t512 = rSges(2,1) * t529 + rSges(2,2) * t533;
t502 = t528 * t561 + t547;
t501 = -t527 * t563 + t531 * t566;
t494 = t566 * rSges(3,3) + (rSges(3,1) * t528 + rSges(3,2) * t532) * t525;
t493 = Icges(3,5) * t566 + (Icges(3,1) * t528 + Icges(3,4) * t532) * t525;
t492 = Icges(3,6) * t566 + (Icges(3,4) * t528 + Icges(3,2) * t532) * t525;
t491 = Icges(3,3) * t566 + (Icges(3,5) * t528 + Icges(3,6) * t532) * t525;
t486 = t506 * t531 + t553;
t485 = -t506 * t527 + t529 * t561;
t484 = t504 * t531 - t552;
t483 = -t504 * t527 - t531 * t559;
t478 = t496 * t530 - t554;
t477 = -t496 * t526 - t530 * t560;
t472 = t496 * t523 - t521 * t560;
t471 = -t496 * t521 - t523 * t560;
t468 = rSges(3,1) * t506 - rSges(3,2) * t505 + rSges(3,3) * t562;
t467 = rSges(3,1) * t504 - rSges(3,2) * t503 - rSges(3,3) * t559;
t464 = Icges(3,1) * t506 - Icges(3,4) * t505 + Icges(3,5) * t562;
t463 = Icges(3,1) * t504 - Icges(3,4) * t503 - Icges(3,5) * t559;
t462 = Icges(3,4) * t506 - Icges(3,2) * t505 + Icges(3,6) * t562;
t461 = Icges(3,4) * t504 - Icges(3,2) * t503 - Icges(3,6) * t559;
t457 = rSges(4,1) * t502 + rSges(4,2) * t501 - rSges(4,3) * t560;
t454 = Icges(4,1) * t502 + Icges(4,4) * t501 - Icges(4,5) * t560;
t453 = Icges(4,4) * t502 + Icges(4,2) * t501 - Icges(4,6) * t560;
t452 = Icges(4,5) * t502 + Icges(4,6) * t501 - Icges(4,3) * t560;
t450 = rSges(5,1) * t496 - rSges(5,2) * t495 - rSges(5,3) * t560;
t449 = Icges(5,1) * t496 - Icges(5,4) * t495 - Icges(5,5) * t560;
t448 = Icges(5,4) * t496 - Icges(5,2) * t495 - Icges(5,6) * t560;
t447 = Icges(5,5) * t496 - Icges(5,6) * t495 - Icges(5,3) * t560;
t446 = t482 * t530 + t564;
t445 = -t482 * t526 + t505 * t530;
t444 = t480 * t530 + t565;
t443 = -t480 * t526 + t503 * t530;
t442 = t482 * t523 + t505 * t521;
t441 = -t482 * t521 + t505 * t523;
t440 = t480 * t523 + t503 * t521;
t439 = -t480 * t521 + t503 * t523;
t438 = qJD(6) * t495 + t451;
t432 = rSges(4,1) * t486 + rSges(4,2) * t485 + rSges(4,3) * t505;
t431 = rSges(4,1) * t484 + rSges(4,2) * t483 + rSges(4,3) * t503;
t430 = Icges(4,1) * t486 + Icges(4,4) * t485 + Icges(4,5) * t505;
t429 = Icges(4,1) * t484 + Icges(4,4) * t483 + Icges(4,5) * t503;
t428 = Icges(4,4) * t486 + Icges(4,2) * t485 + Icges(4,6) * t505;
t427 = Icges(4,4) * t484 + Icges(4,2) * t483 + Icges(4,6) * t503;
t426 = Icges(4,5) * t486 + Icges(4,6) * t485 + Icges(4,3) * t505;
t425 = Icges(4,5) * t484 + Icges(4,6) * t483 + Icges(4,3) * t503;
t422 = rSges(5,1) * t482 - rSges(5,2) * t481 + rSges(5,3) * t505;
t421 = rSges(5,1) * t480 - rSges(5,2) * t479 + rSges(5,3) * t503;
t420 = Icges(5,1) * t482 - Icges(5,4) * t481 + Icges(5,5) * t505;
t419 = Icges(5,1) * t480 - Icges(5,4) * t479 + Icges(5,5) * t503;
t418 = Icges(5,4) * t482 - Icges(5,2) * t481 + Icges(5,6) * t505;
t417 = Icges(5,4) * t480 - Icges(5,2) * t479 + Icges(5,6) * t503;
t416 = Icges(5,5) * t482 - Icges(5,6) * t481 + Icges(5,3) * t505;
t415 = Icges(5,5) * t480 - Icges(5,6) * t479 + Icges(5,3) * t503;
t413 = rSges(6,1) * t478 + rSges(6,2) * t477 + rSges(6,3) * t495;
t412 = Icges(6,1) * t478 + Icges(6,4) * t477 + Icges(6,5) * t495;
t411 = Icges(6,4) * t478 + Icges(6,2) * t477 + Icges(6,6) * t495;
t410 = Icges(6,5) * t478 + Icges(6,6) * t477 + Icges(6,3) * t495;
t408 = t468 * t517 - t494 * t516 + t509;
t407 = -t467 * t517 - t494 * t551 - t556;
t406 = rSges(7,1) * t472 + rSges(7,2) * t471 + rSges(7,3) * t495;
t405 = Icges(7,1) * t472 + Icges(7,4) * t471 + Icges(7,5) * t495;
t404 = Icges(7,4) * t472 + Icges(7,2) * t471 + Icges(7,6) * t495;
t403 = Icges(7,5) * t472 + Icges(7,6) * t471 + Icges(7,3) * t495;
t401 = (t467 * t529 + t468 * t533) * t555;
t400 = -pkin(5) * t554 + pkin(12) * t495 + t496 * t569;
t398 = qJD(6) * t479 + t434;
t397 = qJD(6) * t481 + t433;
t395 = rSges(6,1) * t446 + rSges(6,2) * t445 + rSges(6,3) * t481;
t394 = rSges(6,1) * t444 + rSges(6,2) * t443 + rSges(6,3) * t479;
t393 = Icges(6,1) * t446 + Icges(6,4) * t445 + Icges(6,5) * t481;
t392 = Icges(6,1) * t444 + Icges(6,4) * t443 + Icges(6,5) * t479;
t391 = Icges(6,4) * t446 + Icges(6,2) * t445 + Icges(6,6) * t481;
t390 = Icges(6,4) * t444 + Icges(6,2) * t443 + Icges(6,6) * t479;
t389 = Icges(6,5) * t446 + Icges(6,6) * t445 + Icges(6,3) * t481;
t388 = Icges(6,5) * t444 + Icges(6,6) * t443 + Icges(6,3) * t479;
t387 = rSges(7,1) * t442 + rSges(7,2) * t441 + rSges(7,3) * t481;
t386 = rSges(7,1) * t440 + rSges(7,2) * t439 + rSges(7,3) * t479;
t385 = Icges(7,1) * t442 + Icges(7,4) * t441 + Icges(7,5) * t481;
t384 = Icges(7,1) * t440 + Icges(7,4) * t439 + Icges(7,5) * t479;
t383 = Icges(7,4) * t442 + Icges(7,2) * t441 + Icges(7,6) * t481;
t382 = Icges(7,4) * t440 + Icges(7,2) * t439 + Icges(7,6) * t479;
t381 = Icges(7,5) * t442 + Icges(7,6) * t441 + Icges(7,3) * t481;
t380 = Icges(7,5) * t440 + Icges(7,6) * t439 + Icges(7,3) * t479;
t379 = pkin(5) * t564 + pkin(12) * t481 + t482 * t569;
t378 = pkin(5) * t565 + pkin(12) * t479 + t480 * t569;
t377 = t432 * t508 - t457 * t487 + t543;
t376 = -t431 * t508 + t457 * t488 + t541;
t375 = t431 * t487 - t432 * t488 + t557;
t374 = t422 * t489 - t450 * t455 + t540;
t373 = -t421 * t489 + t450 * t456 + t539;
t372 = t421 * t455 - t422 * t456 + t544;
t371 = t395 * t451 - t413 * t433 + t538;
t370 = -t394 * t451 + t413 * t434 + t537;
t369 = t394 * t433 - t395 * t434 + t542;
t368 = t379 * t451 + t387 * t438 - t397 * t406 - t400 * t433 + t538;
t367 = -t378 * t451 - t386 * t438 + t398 * t406 + t400 * t434 + t537;
t366 = t378 * t433 - t379 * t434 + t386 * t397 - t387 * t398 + t542;
t1 = t455 * ((t505 * t416 - t481 * t418 + t482 * t420) * t455 + (t415 * t505 - t417 * t481 + t419 * t482) * t456 + (t447 * t505 - t448 * t481 + t449 * t482) * t489) / 0.2e1 + t456 * ((t416 * t503 - t418 * t479 + t420 * t480) * t455 + (t503 * t415 - t479 * t417 + t480 * t419) * t456 + (t447 * t503 - t448 * t479 + t449 * t480) * t489) / 0.2e1 + t489 * ((-t416 * t560 - t418 * t495 + t420 * t496) * t455 + (-t415 * t560 - t417 * t495 + t419 * t496) * t456 + (-t447 * t560 - t495 * t448 + t496 * t449) * t489) / 0.2e1 + t488 * ((t426 * t503 + t428 * t483 + t430 * t484) * t487 + (t503 * t425 + t483 * t427 + t484 * t429) * t488 + (t452 * t503 + t453 * t483 + t454 * t484) * t508) / 0.2e1 + t508 * ((-t426 * t560 + t428 * t501 + t430 * t502) * t487 + (-t425 * t560 + t427 * t501 + t429 * t502) * t488 + (-t452 * t560 + t501 * t453 + t502 * t454) * t508) / 0.2e1 + t487 * ((t505 * t426 + t485 * t428 + t486 * t430) * t487 + (t425 * t505 + t427 * t485 + t429 * t486) * t488 + (t452 * t505 + t453 * t485 + t454 * t486) * t508) / 0.2e1 + t517 * ((t566 * t460 + (t462 * t532 + t464 * t528) * t525) * t516 - (t566 * t459 + (t461 * t532 + t463 * t528) * t525) * t551 + (t566 * t491 + (t492 * t532 + t493 * t528) * t525) * t517) / 0.2e1 + m(7) * (t366 ^ 2 + t367 ^ 2 + t368 ^ 2) / 0.2e1 + m(6) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + m(5) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(4) * (t375 ^ 2 + t376 ^ 2 + t377 ^ 2) / 0.2e1 + m(3) * (t401 ^ 2 + t407 ^ 2 + t408 ^ 2) / 0.2e1 + t398 * ((t381 * t479 + t383 * t439 + t385 * t440) * t397 + (t479 * t380 + t439 * t382 + t440 * t384) * t398 + (t403 * t479 + t404 * t439 + t405 * t440) * t438) / 0.2e1 + t397 * ((t481 * t381 + t441 * t383 + t442 * t385) * t397 + (t380 * t481 + t382 * t441 + t384 * t442) * t398 + (t403 * t481 + t404 * t441 + t405 * t442) * t438) / 0.2e1 + t433 * ((t481 * t389 + t445 * t391 + t446 * t393) * t433 + (t388 * t481 + t390 * t445 + t392 * t446) * t434 + (t410 * t481 + t411 * t445 + t412 * t446) * t451) / 0.2e1 + t438 * ((t381 * t495 + t383 * t471 + t385 * t472) * t397 + (t380 * t495 + t382 * t471 + t384 * t472) * t398 + (t495 * t403 + t471 * t404 + t472 * t405) * t438) / 0.2e1 + t434 * ((t389 * t479 + t391 * t443 + t393 * t444) * t433 + (t479 * t388 + t443 * t390 + t444 * t392) * t434 + (t410 * t479 + t411 * t443 + t412 * t444) * t451) / 0.2e1 + t451 * ((t389 * t495 + t391 * t477 + t393 * t478) * t433 + (t388 * t495 + t390 * t477 + t392 * t478) * t434 + (t495 * t410 + t477 * t411 + t478 * t412) * t451) / 0.2e1 - ((-t491 * t559 - t492 * t503 + t493 * t504) * t517 + ((-t462 * t503 + t464 * t504) * t529 + (t503 * t461 - t504 * t463 + t572) * t533) * t555) * t551 / 0.2e1 + ((t491 * t562 - t492 * t505 + t493 * t506) * t517 + (-(-t461 * t505 + t463 * t506) * t533 + (-t505 * t462 + t506 * t464 - t572) * t529) * t555) * t516 / 0.2e1 + (Icges(2,3) + m(2) * (t512 ^ 2 + t513 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
