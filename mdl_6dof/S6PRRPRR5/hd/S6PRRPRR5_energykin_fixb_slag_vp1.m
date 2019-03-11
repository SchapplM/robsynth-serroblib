% Calculate kinetic energy for
% S6PRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:16:17
% EndTime: 2019-03-08 22:16:19
% DurationCPUTime: 2.98s
% Computational Cost: add. (3128->375), mult. (6519->568), div. (0->0), fcn. (8015->14), ass. (0->166)
t571 = Icges(4,2) + Icges(5,3);
t526 = sin(pkin(11));
t529 = cos(pkin(11));
t534 = cos(qJ(2));
t530 = cos(pkin(6));
t533 = sin(qJ(2));
t554 = t530 * t533;
t502 = t526 * t534 + t529 * t554;
t532 = sin(qJ(3));
t527 = sin(pkin(6));
t555 = t529 * t527;
t563 = cos(qJ(3));
t488 = t502 * t563 - t532 * t555;
t553 = t530 * t534;
t501 = t526 * t533 - t529 * t553;
t525 = sin(pkin(12));
t528 = cos(pkin(12));
t453 = -t488 * t525 + t501 * t528;
t560 = t501 * t525;
t454 = t488 * t528 + t560;
t546 = t527 * t563;
t487 = t502 * t532 + t529 * t546;
t570 = -Icges(4,4) * t488 + Icges(5,5) * t454 - Icges(4,6) * t501 + Icges(5,6) * t453 + t571 * t487;
t504 = -t526 * t554 + t529 * t534;
t557 = t527 * t532;
t490 = t504 * t563 + t526 * t557;
t503 = t526 * t553 + t529 * t533;
t455 = -t490 * t525 + t503 * t528;
t559 = t503 * t525;
t456 = t490 * t528 + t559;
t489 = t504 * t532 - t526 * t546;
t569 = -Icges(4,4) * t490 + Icges(5,5) * t456 - Icges(4,6) * t503 + Icges(5,6) * t455 + t571 * t489;
t506 = t530 * t532 + t533 * t546;
t556 = t527 * t534;
t485 = -t506 * t525 - t528 * t556;
t547 = t525 * t556;
t486 = t506 * t528 - t547;
t505 = -t530 * t563 + t533 * t557;
t568 = -Icges(4,4) * t506 + Icges(5,5) * t486 + Icges(4,6) * t556 + Icges(5,6) * t485 + t571 * t505;
t567 = qJD(2) ^ 2;
t561 = t528 * pkin(4);
t558 = t526 * t527;
t524 = pkin(12) + qJ(5);
t520 = cos(t524);
t551 = pkin(5) * t520;
t549 = qJD(2) * t527;
t513 = t526 * t549;
t491 = qJD(3) * t503 + t513;
t518 = qJD(2) * t530;
t447 = qJD(5) * t489 + t491;
t545 = t529 * t549;
t473 = pkin(2) * t502 + pkin(8) * t501;
t474 = pkin(2) * t504 + pkin(8) * t503;
t544 = t473 * t513 + t474 * t545 + qJD(1);
t519 = sin(t524);
t543 = pkin(5) * t519;
t492 = qJD(3) * t501 - t545;
t508 = -qJD(3) * t556 + t518;
t441 = pkin(3) * t488 + qJ(4) * t487;
t542 = qJD(4) * t505 + t491 * t441 + t544;
t507 = (pkin(2) * t533 - pkin(8) * t534) * t527;
t541 = t474 * t518 - t507 * t513;
t448 = qJD(5) * t487 + t492;
t484 = qJD(5) * t505 + t508;
t442 = pkin(3) * t490 + qJ(4) * t489;
t540 = qJD(4) * t487 + t508 * t442 + t541;
t539 = (-t473 * t530 - t507 * t555) * qJD(2);
t392 = pkin(4) * t560 + pkin(9) * t487 + t488 * t561;
t393 = pkin(4) * t559 + pkin(9) * t489 + t490 * t561;
t538 = t491 * t392 + (-t393 - t442) * t492 + t542;
t475 = t506 * pkin(3) + t505 * qJ(4);
t537 = qJD(4) * t489 + t492 * t475 + t539;
t424 = -pkin(4) * t547 + pkin(9) * t505 + t506 * t561;
t536 = t508 * t393 + (-t424 - t475) * t491 + t540;
t535 = t492 * t424 + (-t392 - t441) * t508 + t537;
t521 = qJ(6) + t524;
t516 = cos(t521);
t515 = sin(t521);
t496 = t530 * rSges(3,3) + (rSges(3,1) * t533 + rSges(3,2) * t534) * t527;
t495 = Icges(3,5) * t530 + (Icges(3,1) * t533 + Icges(3,4) * t534) * t527;
t494 = Icges(3,6) * t530 + (Icges(3,4) * t533 + Icges(3,2) * t534) * t527;
t493 = Icges(3,3) * t530 + (Icges(3,5) * t533 + Icges(3,6) * t534) * t527;
t479 = t506 * t520 - t519 * t556;
t478 = -t506 * t519 - t520 * t556;
t477 = t506 * t516 - t515 * t556;
t476 = -t506 * t515 - t516 * t556;
t471 = t506 * rSges(4,1) - t505 * rSges(4,2) - rSges(4,3) * t556;
t470 = Icges(4,1) * t506 - Icges(4,4) * t505 - Icges(4,5) * t556;
t468 = Icges(4,5) * t506 - Icges(4,6) * t505 - Icges(4,3) * t556;
t465 = rSges(3,1) * t504 - rSges(3,2) * t503 + rSges(3,3) * t558;
t464 = rSges(3,1) * t502 - rSges(3,2) * t501 - rSges(3,3) * t555;
t463 = Icges(3,1) * t504 - Icges(3,4) * t503 + Icges(3,5) * t558;
t462 = Icges(3,1) * t502 - Icges(3,4) * t501 - Icges(3,5) * t555;
t461 = Icges(3,4) * t504 - Icges(3,2) * t503 + Icges(3,6) * t558;
t460 = Icges(3,4) * t502 - Icges(3,2) * t501 - Icges(3,6) * t555;
t459 = Icges(3,5) * t504 - Icges(3,6) * t503 + Icges(3,3) * t558;
t458 = Icges(3,5) * t502 - Icges(3,6) * t501 - Icges(3,3) * t555;
t457 = qJD(6) * t505 + t484;
t452 = t490 * t520 + t503 * t519;
t451 = -t490 * t519 + t503 * t520;
t450 = t488 * t520 + t501 * t519;
t449 = -t488 * t519 + t501 * t520;
t446 = t490 * t516 + t503 * t515;
t445 = -t490 * t515 + t503 * t516;
t444 = t488 * t516 + t501 * t515;
t443 = -t488 * t515 + t501 * t516;
t438 = (-t464 * t530 - t496 * t555) * qJD(2);
t437 = (t465 * t530 - t496 * t558) * qJD(2);
t436 = rSges(5,1) * t486 + rSges(5,2) * t485 + rSges(5,3) * t505;
t435 = rSges(4,1) * t490 - rSges(4,2) * t489 + rSges(4,3) * t503;
t434 = rSges(4,1) * t488 - rSges(4,2) * t487 + rSges(4,3) * t501;
t433 = Icges(5,1) * t486 + Icges(5,4) * t485 + Icges(5,5) * t505;
t432 = Icges(5,4) * t486 + Icges(5,2) * t485 + Icges(5,6) * t505;
t430 = Icges(4,1) * t490 - Icges(4,4) * t489 + Icges(4,5) * t503;
t429 = Icges(4,1) * t488 - Icges(4,4) * t487 + Icges(4,5) * t501;
t426 = Icges(4,5) * t490 - Icges(4,6) * t489 + Icges(4,3) * t503;
t425 = Icges(4,5) * t488 - Icges(4,6) * t487 + Icges(4,3) * t501;
t423 = rSges(6,1) * t479 + rSges(6,2) * t478 + rSges(6,3) * t505;
t422 = Icges(6,1) * t479 + Icges(6,4) * t478 + Icges(6,5) * t505;
t421 = Icges(6,4) * t479 + Icges(6,2) * t478 + Icges(6,6) * t505;
t420 = Icges(6,5) * t479 + Icges(6,6) * t478 + Icges(6,3) * t505;
t419 = qJD(6) * t487 + t448;
t418 = qJD(6) * t489 + t447;
t416 = rSges(7,1) * t477 + rSges(7,2) * t476 + rSges(7,3) * t505;
t415 = Icges(7,1) * t477 + Icges(7,4) * t476 + Icges(7,5) * t505;
t414 = Icges(7,4) * t477 + Icges(7,2) * t476 + Icges(7,6) * t505;
t413 = Icges(7,5) * t477 + Icges(7,6) * t476 + Icges(7,3) * t505;
t412 = qJD(1) + (t464 * t526 + t465 * t529) * t549;
t410 = pkin(10) * t505 + t506 * t551 - t543 * t556;
t409 = rSges(5,1) * t456 + rSges(5,2) * t455 + rSges(5,3) * t489;
t408 = rSges(5,1) * t454 + rSges(5,2) * t453 + rSges(5,3) * t487;
t407 = Icges(5,1) * t456 + Icges(5,4) * t455 + Icges(5,5) * t489;
t406 = Icges(5,1) * t454 + Icges(5,4) * t453 + Icges(5,5) * t487;
t405 = Icges(5,4) * t456 + Icges(5,2) * t455 + Icges(5,6) * t489;
t404 = Icges(5,4) * t454 + Icges(5,2) * t453 + Icges(5,6) * t487;
t401 = rSges(6,1) * t452 + rSges(6,2) * t451 + rSges(6,3) * t489;
t400 = rSges(6,1) * t450 + rSges(6,2) * t449 + rSges(6,3) * t487;
t399 = Icges(6,1) * t452 + Icges(6,4) * t451 + Icges(6,5) * t489;
t398 = Icges(6,1) * t450 + Icges(6,4) * t449 + Icges(6,5) * t487;
t397 = Icges(6,4) * t452 + Icges(6,2) * t451 + Icges(6,6) * t489;
t396 = Icges(6,4) * t450 + Icges(6,2) * t449 + Icges(6,6) * t487;
t395 = Icges(6,5) * t452 + Icges(6,6) * t451 + Icges(6,3) * t489;
t394 = Icges(6,5) * t450 + Icges(6,6) * t449 + Icges(6,3) * t487;
t391 = rSges(7,1) * t446 + rSges(7,2) * t445 + rSges(7,3) * t489;
t390 = rSges(7,1) * t444 + rSges(7,2) * t443 + rSges(7,3) * t487;
t389 = Icges(7,1) * t446 + Icges(7,4) * t445 + Icges(7,5) * t489;
t388 = Icges(7,1) * t444 + Icges(7,4) * t443 + Icges(7,5) * t487;
t387 = Icges(7,4) * t446 + Icges(7,2) * t445 + Icges(7,6) * t489;
t386 = Icges(7,4) * t444 + Icges(7,2) * t443 + Icges(7,6) * t487;
t385 = Icges(7,5) * t446 + Icges(7,6) * t445 + Icges(7,3) * t489;
t384 = Icges(7,5) * t444 + Icges(7,6) * t443 + Icges(7,3) * t487;
t381 = pkin(10) * t489 + t490 * t551 + t503 * t543;
t380 = pkin(10) * t487 + t488 * t551 + t501 * t543;
t379 = -t434 * t508 + t471 * t492 + t539;
t378 = t435 * t508 - t471 * t491 + t541;
t377 = t434 * t491 - t435 * t492 + t544;
t376 = t436 * t492 + (-t408 - t441) * t508 + t537;
t375 = t409 * t508 + (-t436 - t475) * t491 + t540;
t374 = t408 * t491 + (-t409 - t442) * t492 + t542;
t373 = -t400 * t484 + t423 * t448 + t535;
t372 = t401 * t484 - t423 * t447 + t536;
t371 = t400 * t447 - t401 * t448 + t538;
t370 = -t380 * t484 - t390 * t457 + t410 * t448 + t416 * t419 + t535;
t369 = t381 * t484 + t391 * t457 - t410 * t447 - t416 * t418 + t536;
t368 = t380 * t447 - t381 * t448 + t390 * t418 - t391 * t419 + t538;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t412 ^ 2 + t437 ^ 2 + t438 ^ 2) / 0.2e1 + m(4) * (t377 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + m(5) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + m(6) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + m(7) * (t368 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + t447 * ((t489 * t395 + t451 * t397 + t452 * t399) * t447 + (t394 * t489 + t396 * t451 + t398 * t452) * t448 + (t420 * t489 + t421 * t451 + t422 * t452) * t484) / 0.2e1 + t448 * ((t395 * t487 + t397 * t449 + t399 * t450) * t447 + (t487 * t394 + t449 * t396 + t450 * t398) * t448 + (t420 * t487 + t421 * t449 + t422 * t450) * t484) / 0.2e1 + t484 * ((t395 * t505 + t397 * t478 + t399 * t479) * t447 + (t394 * t505 + t396 * t478 + t398 * t479) * t448 + (t505 * t420 + t478 * t421 + t479 * t422) * t484) / 0.2e1 + t418 * ((t489 * t385 + t445 * t387 + t446 * t389) * t418 + (t384 * t489 + t386 * t445 + t388 * t446) * t419 + (t413 * t489 + t414 * t445 + t415 * t446) * t457) / 0.2e1 + t419 * ((t385 * t487 + t387 * t443 + t389 * t444) * t418 + (t487 * t384 + t443 * t386 + t444 * t388) * t419 + (t413 * t487 + t414 * t443 + t415 * t444) * t457) / 0.2e1 + t457 * ((t385 * t505 + t387 * t476 + t389 * t477) * t418 + (t384 * t505 + t386 * t476 + t388 * t477) * t419 + (t505 * t413 + t476 * t414 + t477 * t415) * t457) / 0.2e1 - t567 * ((-t459 * t555 - t461 * t501 + t463 * t502) * t558 - (-t458 * t555 - t460 * t501 + t462 * t502) * t555 + (-t493 * t555 - t494 * t501 + t495 * t502) * t530) * t555 / 0.2e1 + ((t432 * t455 + t433 * t456 + t468 * t503 + t470 * t490 + t568 * t489) * t508 + (t404 * t455 + t406 * t456 + t425 * t503 + t429 * t490 + t570 * t489) * t492 + (t455 * t405 + t456 * t407 + t503 * t426 + t490 * t430 + t569 * t489) * t491) * t491 / 0.2e1 + ((t432 * t453 + t433 * t454 + t468 * t501 + t470 * t488 + t568 * t487) * t508 + (t453 * t404 + t454 * t406 + t501 * t425 + t488 * t429 + t570 * t487) * t492 + (t405 * t453 + t407 * t454 + t426 * t501 + t430 * t488 + t569 * t487) * t491) * t492 / 0.2e1 + ((t485 * t432 + t486 * t433 - t468 * t556 + t506 * t470 + t568 * t505) * t508 + (t404 * t485 + t406 * t486 - t425 * t556 + t506 * t429 + t570 * t505) * t492 + (t405 * t485 + t407 * t486 - t426 * t556 + t506 * t430 + t569 * t505) * t491) * t508 / 0.2e1 + (t530 * (t530 ^ 2 * t493 + (((t461 * t534 + t463 * t533) * t526 - (t460 * t534 + t462 * t533) * t529) * t527 + (-t458 * t529 + t459 * t526 + t494 * t534 + t495 * t533) * t530) * t527) + ((t459 * t558 - t461 * t503 + t463 * t504) * t558 - (t458 * t558 - t460 * t503 + t462 * t504) * t555 + (t493 * t558 - t494 * t503 + t495 * t504) * t530) * t558) * t567 / 0.2e1;
T  = t1;
