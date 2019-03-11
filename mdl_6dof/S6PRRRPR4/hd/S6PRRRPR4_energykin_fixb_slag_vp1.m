% Calculate kinetic energy for
% S6PRRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR4_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:15:37
% EndTime: 2019-03-08 23:15:40
% DurationCPUTime: 2.93s
% Computational Cost: add. (3200->375), mult. (6681->568), div. (0->0), fcn. (8222->14), ass. (0->166)
t571 = Icges(5,3) + Icges(6,3);
t525 = sin(pkin(11));
t527 = cos(pkin(11));
t534 = cos(qJ(2));
t528 = cos(pkin(6));
t532 = sin(qJ(2));
t553 = t528 * t532;
t502 = t525 * t534 + t527 * t553;
t531 = sin(qJ(3));
t526 = sin(pkin(6));
t554 = t527 * t526;
t563 = cos(qJ(3));
t486 = t502 * t563 - t531 * t554;
t552 = t528 * t534;
t501 = t525 * t532 - t527 * t552;
t524 = qJ(4) + pkin(12);
t519 = sin(t524);
t520 = cos(t524);
t449 = -t486 * t519 + t501 * t520;
t450 = t486 * t520 + t501 * t519;
t530 = sin(qJ(4));
t533 = cos(qJ(4));
t453 = -t486 * t530 + t501 * t533;
t559 = t501 * t530;
t454 = t486 * t533 + t559;
t546 = t526 * t563;
t485 = t502 * t531 + t527 * t546;
t570 = Icges(5,5) * t454 + Icges(6,5) * t450 + Icges(5,6) * t453 + Icges(6,6) * t449 + t571 * t485;
t504 = -t525 * t553 + t527 * t534;
t556 = t526 * t531;
t488 = t504 * t563 + t525 * t556;
t503 = t525 * t552 + t527 * t532;
t451 = -t488 * t519 + t503 * t520;
t452 = t488 * t520 + t503 * t519;
t455 = -t488 * t530 + t503 * t533;
t558 = t503 * t530;
t456 = t488 * t533 + t558;
t487 = t504 * t531 - t525 * t546;
t569 = Icges(5,5) * t456 + Icges(6,5) * t452 + Icges(5,6) * t455 + Icges(6,6) * t451 + t571 * t487;
t506 = t528 * t531 + t532 * t546;
t555 = t526 * t534;
t478 = -t506 * t519 - t520 * t555;
t479 = t506 * t520 - t519 * t555;
t489 = -t506 * t530 - t533 * t555;
t547 = t530 * t555;
t490 = t506 * t533 - t547;
t505 = -t528 * t563 + t532 * t556;
t568 = Icges(5,5) * t490 + Icges(6,5) * t479 + Icges(5,6) * t489 + Icges(6,6) * t478 + t571 * t505;
t567 = qJD(2) ^ 2;
t561 = t533 * pkin(4);
t557 = t525 * t526;
t551 = pkin(5) * t520;
t549 = qJD(2) * t526;
t513 = t525 * t549;
t491 = qJD(3) * t503 + t513;
t518 = qJD(2) * t528;
t447 = qJD(4) * t487 + t491;
t545 = t527 * t549;
t473 = pkin(2) * t502 + pkin(8) * t501;
t474 = pkin(2) * t504 + pkin(8) * t503;
t544 = t473 * t513 + t474 * t545 + qJD(1);
t543 = pkin(5) * t519;
t492 = qJD(3) * t501 - t545;
t508 = -qJD(3) * t555 + t518;
t507 = (pkin(2) * t532 - pkin(8) * t534) * t526;
t542 = t474 * t518 - t507 * t513;
t448 = qJD(4) * t485 + t492;
t484 = qJD(4) * t505 + t508;
t441 = pkin(3) * t486 + pkin(9) * t485;
t442 = pkin(3) * t488 + pkin(9) * t487;
t541 = t491 * t441 - t442 * t492 + t544;
t540 = (-t473 * t528 - t507 * t554) * qJD(2);
t392 = pkin(4) * t559 + qJ(5) * t485 + t486 * t561;
t539 = qJD(5) * t505 + t447 * t392 + t541;
t477 = t506 * pkin(3) + t505 * pkin(9);
t538 = t508 * t442 - t477 * t491 + t542;
t393 = pkin(4) * t558 + qJ(5) * t487 + t488 * t561;
t537 = qJD(5) * t485 + t484 * t393 + t538;
t536 = -t441 * t508 + t492 * t477 + t540;
t424 = -pkin(4) * t547 + qJ(5) * t505 + t506 * t561;
t535 = qJD(5) * t487 + t448 * t424 + t536;
t521 = qJ(6) + t524;
t516 = cos(t521);
t515 = sin(t521);
t496 = t528 * rSges(3,3) + (rSges(3,1) * t532 + rSges(3,2) * t534) * t526;
t495 = Icges(3,5) * t528 + (Icges(3,1) * t532 + Icges(3,4) * t534) * t526;
t494 = Icges(3,6) * t528 + (Icges(3,4) * t532 + Icges(3,2) * t534) * t526;
t493 = Icges(3,3) * t528 + (Icges(3,5) * t532 + Icges(3,6) * t534) * t526;
t476 = t506 * t516 - t515 * t555;
t475 = -t506 * t515 - t516 * t555;
t471 = t506 * rSges(4,1) - t505 * rSges(4,2) - rSges(4,3) * t555;
t470 = Icges(4,1) * t506 - Icges(4,4) * t505 - Icges(4,5) * t555;
t469 = Icges(4,4) * t506 - Icges(4,2) * t505 - Icges(4,6) * t555;
t468 = Icges(4,5) * t506 - Icges(4,6) * t505 - Icges(4,3) * t555;
t465 = rSges(3,1) * t504 - rSges(3,2) * t503 + rSges(3,3) * t557;
t464 = rSges(3,1) * t502 - rSges(3,2) * t501 - rSges(3,3) * t554;
t463 = Icges(3,1) * t504 - Icges(3,4) * t503 + Icges(3,5) * t557;
t462 = Icges(3,1) * t502 - Icges(3,4) * t501 - Icges(3,5) * t554;
t461 = Icges(3,4) * t504 - Icges(3,2) * t503 + Icges(3,6) * t557;
t460 = Icges(3,4) * t502 - Icges(3,2) * t501 - Icges(3,6) * t554;
t459 = Icges(3,5) * t504 - Icges(3,6) * t503 + Icges(3,3) * t557;
t458 = Icges(3,5) * t502 - Icges(3,6) * t501 - Icges(3,3) * t554;
t457 = qJD(6) * t505 + t484;
t446 = t488 * t516 + t503 * t515;
t445 = -t488 * t515 + t503 * t516;
t444 = t486 * t516 + t501 * t515;
t443 = -t486 * t515 + t501 * t516;
t438 = (-t464 * t528 - t496 * t554) * qJD(2);
t437 = (t465 * t528 - t496 * t557) * qJD(2);
t436 = rSges(5,1) * t490 + rSges(5,2) * t489 + rSges(5,3) * t505;
t435 = Icges(5,1) * t490 + Icges(5,4) * t489 + Icges(5,5) * t505;
t434 = Icges(5,4) * t490 + Icges(5,2) * t489 + Icges(5,6) * t505;
t432 = rSges(4,1) * t488 - rSges(4,2) * t487 + rSges(4,3) * t503;
t431 = rSges(4,1) * t486 - rSges(4,2) * t485 + rSges(4,3) * t501;
t430 = Icges(4,1) * t488 - Icges(4,4) * t487 + Icges(4,5) * t503;
t429 = Icges(4,1) * t486 - Icges(4,4) * t485 + Icges(4,5) * t501;
t428 = Icges(4,4) * t488 - Icges(4,2) * t487 + Icges(4,6) * t503;
t427 = Icges(4,4) * t486 - Icges(4,2) * t485 + Icges(4,6) * t501;
t426 = Icges(4,5) * t488 - Icges(4,6) * t487 + Icges(4,3) * t503;
t425 = Icges(4,5) * t486 - Icges(4,6) * t485 + Icges(4,3) * t501;
t423 = rSges(6,1) * t479 + rSges(6,2) * t478 + rSges(6,3) * t505;
t422 = Icges(6,1) * t479 + Icges(6,4) * t478 + Icges(6,5) * t505;
t421 = Icges(6,4) * t479 + Icges(6,2) * t478 + Icges(6,6) * t505;
t419 = qJD(6) * t485 + t448;
t418 = qJD(6) * t487 + t447;
t416 = rSges(7,1) * t476 + rSges(7,2) * t475 + rSges(7,3) * t505;
t415 = Icges(7,1) * t476 + Icges(7,4) * t475 + Icges(7,5) * t505;
t414 = Icges(7,4) * t476 + Icges(7,2) * t475 + Icges(7,6) * t505;
t413 = Icges(7,5) * t476 + Icges(7,6) * t475 + Icges(7,3) * t505;
t412 = qJD(1) + (t464 * t525 + t465 * t527) * t549;
t411 = pkin(10) * t505 + t506 * t551 - t543 * t555;
t410 = rSges(5,1) * t456 + rSges(5,2) * t455 + rSges(5,3) * t487;
t409 = rSges(5,1) * t454 + rSges(5,2) * t453 + rSges(5,3) * t485;
t408 = Icges(5,1) * t456 + Icges(5,4) * t455 + Icges(5,5) * t487;
t407 = Icges(5,1) * t454 + Icges(5,4) * t453 + Icges(5,5) * t485;
t406 = Icges(5,4) * t456 + Icges(5,2) * t455 + Icges(5,6) * t487;
t405 = Icges(5,4) * t454 + Icges(5,2) * t453 + Icges(5,6) * t485;
t401 = rSges(6,1) * t452 + rSges(6,2) * t451 + rSges(6,3) * t487;
t400 = rSges(6,1) * t450 + rSges(6,2) * t449 + rSges(6,3) * t485;
t399 = Icges(6,1) * t452 + Icges(6,4) * t451 + Icges(6,5) * t487;
t398 = Icges(6,1) * t450 + Icges(6,4) * t449 + Icges(6,5) * t485;
t397 = Icges(6,4) * t452 + Icges(6,2) * t451 + Icges(6,6) * t487;
t396 = Icges(6,4) * t450 + Icges(6,2) * t449 + Icges(6,6) * t485;
t391 = rSges(7,1) * t446 + rSges(7,2) * t445 + rSges(7,3) * t487;
t390 = rSges(7,1) * t444 + rSges(7,2) * t443 + rSges(7,3) * t485;
t389 = Icges(7,1) * t446 + Icges(7,4) * t445 + Icges(7,5) * t487;
t388 = Icges(7,1) * t444 + Icges(7,4) * t443 + Icges(7,5) * t485;
t387 = Icges(7,4) * t446 + Icges(7,2) * t445 + Icges(7,6) * t487;
t386 = Icges(7,4) * t444 + Icges(7,2) * t443 + Icges(7,6) * t485;
t385 = Icges(7,5) * t446 + Icges(7,6) * t445 + Icges(7,3) * t487;
t384 = Icges(7,5) * t444 + Icges(7,6) * t443 + Icges(7,3) * t485;
t382 = pkin(10) * t487 + t488 * t551 + t503 * t543;
t381 = pkin(10) * t485 + t486 * t551 + t501 * t543;
t379 = -t431 * t508 + t471 * t492 + t540;
t378 = t432 * t508 - t471 * t491 + t542;
t377 = t431 * t491 - t432 * t492 + t544;
t376 = -t409 * t484 + t436 * t448 + t536;
t375 = t410 * t484 - t436 * t447 + t538;
t374 = t409 * t447 - t410 * t448 + t541;
t373 = t423 * t448 + (-t392 - t400) * t484 + t535;
t372 = t401 * t484 + (-t423 - t424) * t447 + t537;
t371 = t400 * t447 + (-t393 - t401) * t448 + t539;
t370 = -t390 * t457 + t411 * t448 + t416 * t419 + (-t381 - t392) * t484 + t535;
t369 = t382 * t484 + t391 * t457 - t416 * t418 + (-t411 - t424) * t447 + t537;
t368 = t381 * t447 + t390 * t418 - t391 * t419 + (-t382 - t393) * t448 + t539;
t1 = -t567 * ((-t459 * t554 - t461 * t501 + t463 * t502) * t557 - (-t458 * t554 - t460 * t501 + t462 * t502) * t554 + (-t493 * t554 - t494 * t501 + t495 * t502) * t528) * t554 / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t412 ^ 2 + t437 ^ 2 + t438 ^ 2) / 0.2e1 + m(4) * (t377 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + m(5) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + m(6) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + m(7) * (t368 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + t491 * ((t503 * t426 - t487 * t428 + t488 * t430) * t491 + (t425 * t503 - t427 * t487 + t429 * t488) * t492 + (t468 * t503 - t469 * t487 + t470 * t488) * t508) / 0.2e1 + t492 * ((t426 * t501 - t428 * t485 + t430 * t486) * t491 + (t501 * t425 - t485 * t427 + t486 * t429) * t492 + (t468 * t501 - t469 * t485 + t470 * t486) * t508) / 0.2e1 + t508 * ((-t426 * t555 - t505 * t428 + t506 * t430) * t491 + (-t425 * t555 - t505 * t427 + t506 * t429) * t492 + (-t468 * t555 - t505 * t469 + t506 * t470) * t508) / 0.2e1 + t418 * ((t487 * t385 + t445 * t387 + t446 * t389) * t418 + (t384 * t487 + t386 * t445 + t388 * t446) * t419 + (t413 * t487 + t414 * t445 + t415 * t446) * t457) / 0.2e1 + t419 * ((t385 * t485 + t387 * t443 + t389 * t444) * t418 + (t485 * t384 + t443 * t386 + t444 * t388) * t419 + (t413 * t485 + t414 * t443 + t415 * t444) * t457) / 0.2e1 + t457 * ((t385 * t505 + t387 * t475 + t389 * t476) * t418 + (t384 * t505 + t386 * t475 + t388 * t476) * t419 + (t505 * t413 + t475 * t414 + t476 * t415) * t457) / 0.2e1 + ((t421 * t451 + t422 * t452 + t434 * t455 + t435 * t456 + t568 * t487) * t484 + (t396 * t451 + t398 * t452 + t405 * t455 + t407 * t456 + t570 * t487) * t448 + (t451 * t397 + t452 * t399 + t455 * t406 + t456 * t408 + t569 * t487) * t447) * t447 / 0.2e1 + ((t421 * t449 + t422 * t450 + t434 * t453 + t435 * t454 + t568 * t485) * t484 + (t449 * t396 + t450 * t398 + t453 * t405 + t454 * t407 + t570 * t485) * t448 + (t397 * t449 + t399 * t450 + t406 * t453 + t408 * t454 + t569 * t485) * t447) * t448 / 0.2e1 + ((t478 * t421 + t479 * t422 + t489 * t434 + t490 * t435 + t568 * t505) * t484 + (t396 * t478 + t398 * t479 + t405 * t489 + t407 * t490 + t570 * t505) * t448 + (t397 * t478 + t399 * t479 + t406 * t489 + t408 * t490 + t569 * t505) * t447) * t484 / 0.2e1 + (((t459 * t557 - t461 * t503 + t463 * t504) * t557 - (t458 * t557 - t460 * t503 + t462 * t504) * t554 + (t493 * t557 - t494 * t503 + t495 * t504) * t528) * t557 + t528 * (t528 ^ 2 * t493 + (((t461 * t534 + t463 * t532) * t525 - (t460 * t534 + t462 * t532) * t527) * t526 + (-t458 * t527 + t459 * t525 + t494 * t534 + t495 * t532) * t528) * t526)) * t567 / 0.2e1;
T  = t1;
