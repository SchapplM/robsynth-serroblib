% Calculate kinetic energy for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:42:56
% EndTime: 2019-03-08 19:42:59
% DurationCPUTime: 2.97s
% Computational Cost: add. (2654->335), mult. (4935->507), div. (0->0), fcn. (5882->12), ass. (0->154)
t565 = Icges(5,1) + Icges(6,2);
t564 = Icges(6,1) + Icges(5,3);
t563 = -Icges(5,4) - Icges(6,6);
t562 = Icges(6,4) - Icges(5,5);
t561 = Icges(6,5) - Icges(5,6);
t560 = Icges(5,2) + Icges(6,3);
t503 = sin(pkin(10));
t506 = cos(pkin(10));
t512 = cos(qJ(2));
t507 = cos(pkin(6));
t510 = sin(qJ(2));
t539 = t507 * t510;
t490 = t503 * t512 + t506 * t539;
t504 = sin(pkin(6));
t532 = pkin(11) + qJ(4);
t526 = cos(t532);
t520 = t504 * t526;
t525 = sin(t532);
t463 = t490 * t525 + t506 * t520;
t519 = t504 * t525;
t464 = t490 * t526 - t506 * t519;
t538 = t507 * t512;
t489 = t503 * t510 - t506 * t538;
t559 = t560 * t463 + t563 * t464 + t561 * t489;
t492 = -t503 * t539 + t506 * t512;
t465 = t492 * t525 - t503 * t520;
t466 = t492 * t526 + t503 * t519;
t491 = t503 * t538 + t506 * t510;
t558 = t560 * t465 + t563 * t466 + t561 * t491;
t557 = t561 * t463 - t562 * t464 + t564 * t489;
t556 = t561 * t465 - t562 * t466 + t564 * t491;
t555 = t563 * t463 + t565 * t464 - t562 * t489;
t554 = t563 * t465 + t565 * t466 - t562 * t491;
t480 = -t507 * t526 + t510 * t519;
t481 = t507 * t525 + t510 * t520;
t540 = t504 * t512;
t553 = t560 * t480 + t563 * t481 - t561 * t540;
t552 = t563 * t480 + t565 * t481 + t562 * t540;
t551 = t561 * t480 - t562 * t481 - t564 * t540;
t502 = sin(pkin(11));
t505 = cos(pkin(11));
t541 = t504 * t510;
t487 = -t502 * t541 + t505 * t507;
t544 = t502 * t507;
t488 = t505 * t541 + t544;
t441 = Icges(4,5) * t488 + Icges(4,6) * t487 - Icges(4,3) * t540;
t478 = Icges(3,6) * t507 + (Icges(3,4) * t510 + Icges(3,2) * t512) * t504;
t550 = t441 - t478;
t549 = qJD(2) ^ 2;
t545 = pkin(3) * t505;
t543 = t503 * t504;
t542 = t504 * t506;
t459 = t492 * pkin(2) + t491 * qJ(3);
t501 = qJD(2) * t507;
t536 = qJD(3) * t489 + t459 * t501;
t535 = qJD(2) * t504;
t499 = t503 * t535;
t474 = qJD(4) * t491 + t499;
t534 = qJD(3) * t512;
t531 = t502 * t543;
t530 = t502 * t542;
t529 = t506 * t535;
t458 = t490 * pkin(2) + t489 * qJ(3);
t528 = t458 * t499 + t459 * t529 + qJD(1);
t493 = (pkin(2) * t510 - qJ(3) * t512) * t504;
t524 = (-t488 * rSges(4,1) - t487 * rSges(4,2) + rSges(4,3) * t540 - t493) * t504;
t523 = (-pkin(3) * t544 - (-pkin(8) * t512 + t510 * t545) * t504 - t493) * t504;
t475 = qJD(4) * t489 - t529;
t494 = -qJD(4) * t540 + t501;
t409 = -pkin(3) * t530 + pkin(8) * t489 + t490 * t545;
t410 = pkin(3) * t531 + pkin(8) * t491 + t492 * t545;
t518 = t409 * t499 + t410 * t529 - t504 * t534 + t528;
t517 = qJD(2) * t503 * t523 + t410 * t501 + t536;
t423 = pkin(4) * t464 + qJ(5) * t463;
t516 = qJD(5) * t480 + t474 * t423 + t518;
t424 = pkin(4) * t466 + qJ(5) * t465;
t515 = qJD(5) * t463 + t494 * t424 + t517;
t486 = qJD(3) * t491;
t514 = t486 + ((-t409 - t458) * t507 + t506 * t523) * qJD(2);
t444 = pkin(4) * t481 + qJ(5) * t480;
t513 = qJD(5) * t465 + t475 * t444 + t514;
t511 = cos(qJ(6));
t509 = sin(qJ(6));
t482 = t507 * rSges(3,3) + (rSges(3,1) * t510 + rSges(3,2) * t512) * t504;
t479 = Icges(3,5) * t507 + (Icges(3,1) * t510 + Icges(3,4) * t512) * t504;
t477 = Icges(3,3) * t507 + (Icges(3,5) * t510 + Icges(3,6) * t512) * t504;
t473 = -pkin(5) * t540 + t481 * pkin(9);
t472 = t492 * t505 + t531;
t471 = -t492 * t502 + t505 * t543;
t470 = t490 * t505 - t530;
t469 = -t490 * t502 - t505 * t542;
t468 = t480 * t509 - t511 * t540;
t467 = t480 * t511 + t509 * t540;
t462 = qJD(6) * t481 + t494;
t456 = rSges(3,1) * t492 - rSges(3,2) * t491 + rSges(3,3) * t543;
t455 = rSges(3,1) * t490 - rSges(3,2) * t489 - rSges(3,3) * t542;
t450 = Icges(3,1) * t492 - Icges(3,4) * t491 + Icges(3,5) * t543;
t449 = Icges(3,1) * t490 - Icges(3,4) * t489 - Icges(3,5) * t542;
t448 = Icges(3,4) * t492 - Icges(3,2) * t491 + Icges(3,6) * t543;
t447 = Icges(3,4) * t490 - Icges(3,2) * t489 - Icges(3,6) * t542;
t446 = Icges(3,5) * t492 - Icges(3,6) * t491 + Icges(3,3) * t543;
t445 = Icges(3,5) * t490 - Icges(3,6) * t489 - Icges(3,3) * t542;
t443 = Icges(4,1) * t488 + Icges(4,4) * t487 - Icges(4,5) * t540;
t442 = Icges(4,4) * t488 + Icges(4,2) * t487 - Icges(4,6) * t540;
t440 = t481 * rSges(5,1) - t480 * rSges(5,2) - rSges(5,3) * t540;
t439 = -rSges(6,1) * t540 - t481 * rSges(6,2) + t480 * rSges(6,3);
t432 = pkin(5) * t491 + pkin(9) * t466;
t431 = pkin(5) * t489 + pkin(9) * t464;
t430 = t465 * t509 + t491 * t511;
t429 = t465 * t511 - t491 * t509;
t428 = t463 * t509 + t489 * t511;
t427 = t463 * t511 - t489 * t509;
t426 = qJD(6) * t464 + t475;
t425 = qJD(6) * t466 + t474;
t421 = (-t455 * t507 - t482 * t542) * qJD(2);
t420 = (t456 * t507 - t482 * t543) * qJD(2);
t419 = rSges(4,1) * t472 + rSges(4,2) * t471 + rSges(4,3) * t491;
t418 = rSges(4,1) * t470 + rSges(4,2) * t469 + rSges(4,3) * t489;
t417 = Icges(4,1) * t472 + Icges(4,4) * t471 + Icges(4,5) * t491;
t416 = Icges(4,1) * t470 + Icges(4,4) * t469 + Icges(4,5) * t489;
t415 = Icges(4,4) * t472 + Icges(4,2) * t471 + Icges(4,6) * t491;
t414 = Icges(4,4) * t470 + Icges(4,2) * t469 + Icges(4,6) * t489;
t413 = Icges(4,5) * t472 + Icges(4,6) * t471 + Icges(4,3) * t491;
t412 = Icges(4,5) * t470 + Icges(4,6) * t469 + Icges(4,3) * t489;
t408 = rSges(5,1) * t466 - rSges(5,2) * t465 + rSges(5,3) * t491;
t407 = rSges(5,1) * t464 - rSges(5,2) * t463 + rSges(5,3) * t489;
t406 = rSges(6,1) * t491 - rSges(6,2) * t466 + rSges(6,3) * t465;
t405 = rSges(6,1) * t489 - rSges(6,2) * t464 + rSges(6,3) * t463;
t392 = rSges(7,1) * t468 + rSges(7,2) * t467 + rSges(7,3) * t481;
t391 = Icges(7,1) * t468 + Icges(7,4) * t467 + Icges(7,5) * t481;
t390 = Icges(7,4) * t468 + Icges(7,2) * t467 + Icges(7,6) * t481;
t389 = Icges(7,5) * t468 + Icges(7,6) * t467 + Icges(7,3) * t481;
t384 = qJD(1) + (t455 * t503 + t456 * t506) * t535;
t383 = rSges(7,1) * t430 + rSges(7,2) * t429 + rSges(7,3) * t466;
t382 = rSges(7,1) * t428 + rSges(7,2) * t427 + rSges(7,3) * t464;
t381 = Icges(7,1) * t430 + Icges(7,4) * t429 + Icges(7,5) * t466;
t380 = Icges(7,1) * t428 + Icges(7,4) * t427 + Icges(7,5) * t464;
t379 = Icges(7,4) * t430 + Icges(7,2) * t429 + Icges(7,6) * t466;
t378 = Icges(7,4) * t428 + Icges(7,2) * t427 + Icges(7,6) * t464;
t377 = Icges(7,5) * t430 + Icges(7,6) * t429 + Icges(7,3) * t466;
t376 = Icges(7,5) * t428 + Icges(7,6) * t427 + Icges(7,3) * t464;
t375 = t486 + ((-t418 - t458) * t507 + t506 * t524) * qJD(2);
t374 = (t419 * t507 + t503 * t524) * qJD(2) + t536;
t373 = (-t534 + (t418 * t503 + t419 * t506) * qJD(2)) * t504 + t528;
t372 = -t494 * t407 + t475 * t440 + t514;
t371 = t494 * t408 - t474 * t440 + t517;
t370 = t474 * t407 - t475 * t408 + t518;
t369 = t475 * t439 + (-t405 - t423) * t494 + t513;
t368 = t494 * t406 + (-t439 - t444) * t474 + t515;
t367 = t474 * t405 + (-t406 - t424) * t475 + t516;
t366 = -t462 * t382 + t426 * t392 + t475 * t473 + (-t423 - t431) * t494 + t513;
t365 = t462 * t383 - t425 * t392 + t494 * t432 + (-t444 - t473) * t474 + t515;
t364 = t425 * t382 - t426 * t383 + t474 * t431 + (-t424 - t432) * t475 + t516;
t1 = t425 * ((t466 * t377 + t429 * t379 + t430 * t381) * t425 + (t376 * t466 + t378 * t429 + t380 * t430) * t426 + (t389 * t466 + t390 * t429 + t391 * t430) * t462) / 0.2e1 + t426 * ((t377 * t464 + t379 * t427 + t381 * t428) * t425 + (t464 * t376 + t427 * t378 + t428 * t380) * t426 + (t389 * t464 + t390 * t427 + t391 * t428) * t462) / 0.2e1 + t462 * ((t377 * t481 + t379 * t467 + t381 * t468) * t425 + (t376 * t481 + t378 * t467 + t380 * t468) * t426 + (t481 * t389 + t467 * t390 + t468 * t391) * t462) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t384 ^ 2 + t420 ^ 2 + t421 ^ 2) / 0.2e1 + m(4) * (t373 ^ 2 + t374 ^ 2 + t375 ^ 2) / 0.2e1 + m(5) * (t370 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 + m(6) * (t367 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + m(7) * (t364 ^ 2 + t365 ^ 2 + t366 ^ 2) / 0.2e1 + ((t553 * t465 + t552 * t466 + t551 * t491) * t494 + (t559 * t465 + t555 * t466 + t557 * t491) * t475 + (t558 * t465 + t554 * t466 + t556 * t491) * t474) * t474 / 0.2e1 + ((t553 * t463 + t552 * t464 + t551 * t489) * t494 + (t559 * t463 + t555 * t464 + t557 * t489) * t475 + (t558 * t463 + t554 * t464 + t556 * t489) * t474) * t475 / 0.2e1 + ((t553 * t480 + t552 * t481 - t551 * t540) * t494 + (t559 * t480 + t555 * t481 - t557 * t540) * t475 + (t558 * t480 + t554 * t481 - t556 * t540) * t474) * t494 / 0.2e1 - (((t413 * t489 + t415 * t469 + t417 * t470) * t503 - (t412 * t489 + t414 * t469 + t416 * t470) * t506) * t504 + (-t446 * t542 - t448 * t489 + t450 * t490) * t543 - (-t445 * t542 - t447 * t489 + t449 * t490) * t542 + (t442 * t469 + t443 * t470 - t477 * t542 + t479 * t490 + t550 * t489) * t507) * t549 * t542 / 0.2e1 + ((((t448 * t512 + t450 * t510) * t503 - (t447 * t512 + t449 * t510) * t506) * t504 ^ 2 + (-t413 * t540 + t487 * t415 + t488 * t417) * t543 - (-t412 * t540 + t487 * t414 + t488 * t416) * t542 + ((-t445 * t506 + t446 * t503 + t478 * t512 + t479 * t510) * t504 - t441 * t540 + t487 * t442 + t488 * t443 + t507 * t477) * t507) * t507 + ((t446 * t543 - t448 * t491 + t450 * t492) * t543 - (t445 * t543 - t447 * t491 + t449 * t492) * t542 + ((t413 * t491 + t415 * t471 + t417 * t472) * t503 - (t412 * t491 + t414 * t471 + t416 * t472) * t506) * t504 + (t442 * t471 + t443 * t472 + t477 * t543 + t479 * t492 + t550 * t491) * t507) * t543) * t549 / 0.2e1;
T  = t1;
