% Calculate kinetic energy for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:38:15
% EndTime: 2019-03-09 15:38:17
% DurationCPUTime: 2.95s
% Computational Cost: add. (3308->378), mult. (5693->548), div. (0->0), fcn. (6856->14), ass. (0->171)
t588 = Icges(5,2) + Icges(6,3);
t587 = Icges(4,3) + Icges(5,3);
t533 = sin(qJ(2));
t534 = sin(qJ(1));
t536 = cos(qJ(2));
t537 = cos(qJ(1));
t572 = cos(pkin(6));
t551 = t537 * t572;
t505 = t533 * t534 - t536 * t551;
t506 = t533 * t551 + t534 * t536;
t528 = sin(pkin(6));
t565 = t528 * t537;
t462 = Icges(3,5) * t506 - Icges(3,6) * t505 - Icges(3,3) * t565;
t552 = t534 * t572;
t507 = t537 * t533 + t536 * t552;
t508 = -t533 * t552 + t537 * t536;
t568 = t528 * t534;
t463 = Icges(3,5) * t508 - Icges(3,6) * t507 + Icges(3,3) * t568;
t586 = t528 * (t462 * t537 - t463 * t534);
t557 = qJ(3) + pkin(11);
t524 = sin(t557);
t549 = cos(t557);
t483 = t506 * t549 - t524 * t565;
t527 = sin(pkin(12));
t529 = cos(pkin(12));
t448 = -t483 * t527 + t505 * t529;
t571 = t505 * t527;
t449 = t483 * t529 + t571;
t548 = t528 * t549;
t482 = t506 * t524 + t537 * t548;
t585 = -Icges(5,4) * t483 + Icges(6,5) * t449 - Icges(5,6) * t505 + Icges(6,6) * t448 + t588 * t482;
t485 = t508 * t549 + t524 * t568;
t450 = -t485 * t527 + t507 * t529;
t570 = t507 * t527;
t451 = t485 * t529 + t570;
t484 = t508 * t524 - t534 * t548;
t584 = -Icges(5,4) * t485 + Icges(6,5) * t451 - Icges(5,6) * t507 + Icges(6,6) * t450 + t588 * t484;
t497 = t524 * t572 + t533 * t548;
t566 = t528 * t536;
t480 = -t497 * t527 - t529 * t566;
t556 = t527 * t566;
t481 = t497 * t529 - t556;
t569 = t528 * t533;
t496 = t524 * t569 - t549 * t572;
t583 = -Icges(5,4) * t497 + Icges(6,5) * t481 + Icges(5,6) * t566 + Icges(6,6) * t480 + t588 * t496;
t532 = sin(qJ(3));
t535 = cos(qJ(3));
t486 = -t506 * t532 - t535 * t565;
t554 = t532 * t565;
t487 = t506 * t535 - t554;
t582 = Icges(4,5) * t487 + Icges(5,5) * t483 + Icges(4,6) * t486 - Icges(5,6) * t482 + t587 * t505;
t567 = t528 * t535;
t488 = -t508 * t532 + t534 * t567;
t555 = t532 * t568;
t489 = t508 * t535 + t555;
t581 = Icges(4,5) * t489 + Icges(5,5) * t485 + Icges(4,6) * t488 - Icges(5,6) * t484 + t587 * t507;
t503 = -t532 * t569 + t535 * t572;
t550 = t572 * t532;
t504 = t533 * t567 + t550;
t580 = Icges(4,5) * t504 + Icges(5,5) * t497 + Icges(4,6) * t503 - Icges(5,6) * t496 - t587 * t566;
t575 = pkin(3) * t535;
t574 = pkin(5) * t529;
t427 = -pkin(3) * t554 + qJ(4) * t505 + t506 * t575;
t439 = pkin(4) * t483 + qJ(5) * t482;
t563 = -t427 - t439;
t428 = pkin(3) * t555 + qJ(4) * t507 + t508 * t575;
t440 = pkin(4) * t485 + qJ(5) * t484;
t562 = -t428 - t440;
t456 = pkin(4) * t497 + qJ(5) * t496;
t460 = pkin(3) * t550 + (-qJ(4) * t536 + t533 * t575) * t528;
t561 = -t456 - t460;
t478 = pkin(2) * t506 + pkin(9) * t505;
t479 = pkin(2) * t508 + pkin(9) * t507;
t558 = qJD(2) * t528;
t518 = t534 * t558;
t553 = t537 * t558;
t560 = t478 * t518 + t479 * t553;
t490 = qJD(3) * t507 + t518;
t559 = qJD(1) * (pkin(1) * t534 - pkin(8) * t565);
t519 = qJD(2) * t572 + qJD(1);
t491 = qJD(3) * t505 - t553;
t510 = -qJD(3) * t566 + t519;
t509 = (pkin(2) * t533 - pkin(9) * t536) * t528;
t511 = qJD(1) * (pkin(1) * t537 + pkin(8) * t568);
t546 = t519 * t479 - t509 * t518 + t511;
t545 = -qJD(4) * t566 + t490 * t427 + t560;
t544 = qJD(4) * t505 + t510 * t428 + t546;
t543 = qJD(5) * t496 + t490 * t439 + t545;
t542 = -t478 * t519 - t509 * t553 - t559;
t541 = qJD(5) * t482 + t510 * t440 + t544;
t540 = qJD(4) * t507 + t491 * t460 + t542;
t539 = qJD(5) * t484 + t491 * t456 + t540;
t526 = pkin(12) + qJ(6);
t525 = cos(t526);
t523 = sin(t526);
t515 = rSges(2,1) * t537 - rSges(2,2) * t534;
t514 = rSges(2,1) * t534 + rSges(2,2) * t537;
t498 = t572 * rSges(3,3) + (rSges(3,1) * t533 + rSges(3,2) * t536) * t528;
t495 = Icges(3,5) * t572 + (Icges(3,1) * t533 + Icges(3,4) * t536) * t528;
t494 = Icges(3,6) * t572 + (Icges(3,4) * t533 + Icges(3,2) * t536) * t528;
t493 = Icges(3,3) * t572 + (Icges(3,5) * t533 + Icges(3,6) * t536) * t528;
t475 = qJD(6) * t496 + t510;
t474 = t497 * t525 - t523 * t566;
t473 = -t497 * t523 - t525 * t566;
t470 = rSges(3,1) * t508 - rSges(3,2) * t507 + rSges(3,3) * t568;
t469 = rSges(3,1) * t506 - rSges(3,2) * t505 - rSges(3,3) * t565;
t467 = Icges(3,1) * t508 - Icges(3,4) * t507 + Icges(3,5) * t568;
t466 = Icges(3,1) * t506 - Icges(3,4) * t505 - Icges(3,5) * t565;
t465 = Icges(3,4) * t508 - Icges(3,2) * t507 + Icges(3,6) * t568;
t464 = Icges(3,4) * t506 - Icges(3,2) * t505 - Icges(3,6) * t565;
t461 = rSges(4,1) * t504 + rSges(4,2) * t503 - rSges(4,3) * t566;
t459 = Icges(4,1) * t504 + Icges(4,4) * t503 - Icges(4,5) * t566;
t458 = Icges(4,4) * t504 + Icges(4,2) * t503 - Icges(4,6) * t566;
t455 = rSges(5,1) * t497 - rSges(5,2) * t496 - rSges(5,3) * t566;
t454 = Icges(5,1) * t497 - Icges(5,4) * t496 - Icges(5,5) * t566;
t447 = t485 * t525 + t507 * t523;
t446 = -t485 * t523 + t507 * t525;
t445 = t483 * t525 + t505 * t523;
t444 = -t483 * t523 + t505 * t525;
t443 = qJD(6) * t482 + t491;
t442 = qJD(6) * t484 + t490;
t437 = rSges(4,1) * t489 + rSges(4,2) * t488 + rSges(4,3) * t507;
t436 = rSges(4,1) * t487 + rSges(4,2) * t486 + rSges(4,3) * t505;
t435 = Icges(4,1) * t489 + Icges(4,4) * t488 + Icges(4,5) * t507;
t434 = Icges(4,1) * t487 + Icges(4,4) * t486 + Icges(4,5) * t505;
t433 = Icges(4,4) * t489 + Icges(4,2) * t488 + Icges(4,6) * t507;
t432 = Icges(4,4) * t487 + Icges(4,2) * t486 + Icges(4,6) * t505;
t426 = rSges(5,1) * t485 - rSges(5,2) * t484 + rSges(5,3) * t507;
t425 = rSges(5,1) * t483 - rSges(5,2) * t482 + rSges(5,3) * t505;
t424 = Icges(5,1) * t485 - Icges(5,4) * t484 + Icges(5,5) * t507;
t423 = Icges(5,1) * t483 - Icges(5,4) * t482 + Icges(5,5) * t505;
t418 = t470 * t519 - t498 * t518 + t511;
t417 = -t469 * t519 - t498 * t553 - t559;
t416 = rSges(6,1) * t481 + rSges(6,2) * t480 + rSges(6,3) * t496;
t415 = Icges(6,1) * t481 + Icges(6,4) * t480 + Icges(6,5) * t496;
t414 = Icges(6,4) * t481 + Icges(6,2) * t480 + Icges(6,6) * t496;
t410 = (t469 * t534 + t470 * t537) * t558;
t409 = rSges(7,1) * t474 + rSges(7,2) * t473 + rSges(7,3) * t496;
t408 = Icges(7,1) * t474 + Icges(7,4) * t473 + Icges(7,5) * t496;
t407 = Icges(7,4) * t474 + Icges(7,2) * t473 + Icges(7,6) * t496;
t406 = Icges(7,5) * t474 + Icges(7,6) * t473 + Icges(7,3) * t496;
t405 = -pkin(5) * t556 + pkin(10) * t496 + t497 * t574;
t403 = rSges(6,1) * t451 + rSges(6,2) * t450 + rSges(6,3) * t484;
t402 = rSges(6,1) * t449 + rSges(6,2) * t448 + rSges(6,3) * t482;
t401 = Icges(6,1) * t451 + Icges(6,4) * t450 + Icges(6,5) * t484;
t400 = Icges(6,1) * t449 + Icges(6,4) * t448 + Icges(6,5) * t482;
t399 = Icges(6,4) * t451 + Icges(6,2) * t450 + Icges(6,6) * t484;
t398 = Icges(6,4) * t449 + Icges(6,2) * t448 + Icges(6,6) * t482;
t395 = rSges(7,1) * t447 + rSges(7,2) * t446 + rSges(7,3) * t484;
t394 = rSges(7,1) * t445 + rSges(7,2) * t444 + rSges(7,3) * t482;
t393 = Icges(7,1) * t447 + Icges(7,4) * t446 + Icges(7,5) * t484;
t392 = Icges(7,1) * t445 + Icges(7,4) * t444 + Icges(7,5) * t482;
t391 = Icges(7,4) * t447 + Icges(7,2) * t446 + Icges(7,6) * t484;
t390 = Icges(7,4) * t445 + Icges(7,2) * t444 + Icges(7,6) * t482;
t389 = Icges(7,5) * t447 + Icges(7,6) * t446 + Icges(7,3) * t484;
t388 = Icges(7,5) * t445 + Icges(7,6) * t444 + Icges(7,3) * t482;
t387 = pkin(5) * t570 + pkin(10) * t484 + t485 * t574;
t386 = pkin(5) * t571 + pkin(10) * t482 + t483 * t574;
t385 = t437 * t510 - t461 * t490 + t546;
t384 = -t436 * t510 + t461 * t491 + t542;
t383 = t436 * t490 - t437 * t491 + t560;
t382 = t426 * t510 + (-t455 - t460) * t490 + t544;
t381 = t455 * t491 + (-t425 - t427) * t510 + t540;
t380 = t425 * t490 + (-t426 - t428) * t491 + t545;
t379 = t403 * t510 + (-t416 + t561) * t490 + t541;
t378 = t416 * t491 + (-t402 + t563) * t510 + t539;
t377 = t402 * t490 + (-t403 + t562) * t491 + t543;
t376 = t387 * t510 + t395 * t475 - t409 * t442 + (-t405 + t561) * t490 + t541;
t375 = -t394 * t475 + t405 * t491 + t409 * t443 + (-t386 + t563) * t510 + t539;
t374 = t386 * t490 + t394 * t442 - t395 * t443 + (-t387 + t562) * t491 + t543;
t1 = ((t493 * t568 - t494 * t507 + t495 * t508) * t519 + (-(-t464 * t507 + t466 * t508) * t537 + (-t507 * t465 + t508 * t467 - t586) * t534) * t558) * t518 / 0.2e1 - ((-t493 * t565 - t494 * t505 + t495 * t506) * t519 + ((-t465 * t505 + t467 * t506) * t534 + (t505 * t464 - t506 * t466 + t586) * t537) * t558) * t553 / 0.2e1 + t475 * ((t389 * t496 + t391 * t473 + t393 * t474) * t442 + (t388 * t496 + t390 * t473 + t392 * t474) * t443 + (t496 * t406 + t473 * t407 + t474 * t408) * t475) / 0.2e1 + t442 * ((t484 * t389 + t446 * t391 + t447 * t393) * t442 + (t388 * t484 + t390 * t446 + t392 * t447) * t443 + (t406 * t484 + t407 * t446 + t408 * t447) * t475) / 0.2e1 + t443 * ((t389 * t482 + t391 * t444 + t393 * t445) * t442 + (t482 * t388 + t444 * t390 + t445 * t392) * t443 + (t406 * t482 + t407 * t444 + t408 * t445) * t475) / 0.2e1 + t519 * ((t572 * t463 + (t465 * t536 + t467 * t533) * t528) * t518 - (t572 * t462 + (t464 * t536 + t466 * t533) * t528) * t553 + (t572 * t493 + (t494 * t536 + t495 * t533) * t528) * t519) / 0.2e1 + m(6) * (t377 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + m(7) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + m(4) * (t383 ^ 2 + t384 ^ 2 + t385 ^ 2) / 0.2e1 + m(5) * (t380 ^ 2 + t381 ^ 2 + t382 ^ 2) / 0.2e1 + m(3) * (t410 ^ 2 + t417 ^ 2 + t418 ^ 2) / 0.2e1 + (Icges(2,3) + m(2) * (t514 ^ 2 + t515 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t414 * t450 + t415 * t451 + t454 * t485 + t458 * t488 + t459 * t489 + t484 * t583 + t507 * t580) * t510 + (t398 * t450 + t400 * t451 + t423 * t485 + t432 * t488 + t434 * t489 + t484 * t585 + t507 * t582) * t491 + (t399 * t450 + t401 * t451 + t424 * t485 + t433 * t488 + t435 * t489 + t484 * t584 + t507 * t581) * t490) * t490 / 0.2e1 + ((t414 * t448 + t415 * t449 + t454 * t483 + t458 * t486 + t459 * t487 + t482 * t583 + t505 * t580) * t510 + (t398 * t448 + t400 * t449 + t423 * t483 + t432 * t486 + t434 * t487 + t482 * t585 + t505 * t582) * t491 + (t399 * t448 + t401 * t449 + t424 * t483 + t433 * t486 + t435 * t487 + t482 * t584 + t505 * t581) * t490) * t491 / 0.2e1 + ((t414 * t480 + t415 * t481 + t454 * t497 + t458 * t503 + t459 * t504 + t496 * t583 - t566 * t580) * t510 + (t398 * t480 + t400 * t481 + t423 * t497 + t432 * t503 + t434 * t504 + t496 * t585 - t566 * t582) * t491 + (t399 * t480 + t401 * t481 + t424 * t497 + t433 * t503 + t435 * t504 + t496 * t584 - t566 * t581) * t490) * t510 / 0.2e1;
T  = t1;
