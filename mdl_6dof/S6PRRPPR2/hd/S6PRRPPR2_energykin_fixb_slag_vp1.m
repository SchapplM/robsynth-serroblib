% Calculate kinetic energy for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPPR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPPR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:04:05
% EndTime: 2019-03-08 21:04:08
% DurationCPUTime: 2.50s
% Computational Cost: add. (2709->328), mult. (5045->482), div. (0->0), fcn. (6014->12), ass. (0->153)
t565 = Icges(5,1) + Icges(6,2);
t564 = -Icges(5,4) - Icges(6,6);
t563 = Icges(6,4) - Icges(5,5);
t562 = Icges(6,5) - Icges(5,6);
t561 = Icges(5,2) + Icges(6,3);
t560 = Icges(6,1) + Icges(4,3) + Icges(5,3);
t504 = sin(pkin(10));
t506 = cos(pkin(10));
t514 = cos(qJ(2));
t507 = cos(pkin(6));
t511 = sin(qJ(2));
t538 = t507 * t511;
t490 = t504 * t514 + t506 * t538;
t505 = sin(pkin(6));
t531 = qJ(3) + pkin(11);
t526 = cos(t531);
t524 = t505 * t526;
t525 = sin(t531);
t465 = t490 * t525 + t506 * t524;
t523 = t505 * t525;
t466 = t490 * t526 - t506 * t523;
t537 = t507 * t514;
t489 = t504 * t511 - t506 * t537;
t559 = t561 * t465 + t564 * t466 + t562 * t489;
t492 = -t504 * t538 + t506 * t514;
t467 = t492 * t525 - t504 * t524;
t468 = t492 * t526 + t504 * t523;
t491 = t504 * t537 + t506 * t511;
t558 = t561 * t467 + t564 * t468 + t562 * t491;
t557 = t564 * t465 + t565 * t466 - t563 * t489;
t556 = t564 * t467 + t565 * t468 - t563 * t491;
t482 = -t507 * t526 + t511 * t523;
t483 = t507 * t525 + t511 * t524;
t540 = t505 * t514;
t555 = t561 * t482 + t564 * t483 - t562 * t540;
t554 = t564 * t482 + t565 * t483 + t563 * t540;
t510 = sin(qJ(3));
t513 = cos(qJ(3));
t541 = t505 * t513;
t472 = -t490 * t510 - t506 * t541;
t542 = t505 * t510;
t529 = t506 * t542;
t473 = t490 * t513 - t529;
t553 = Icges(4,5) * t473 + Icges(4,6) * t472 + t562 * t465 - t563 * t466 + t560 * t489;
t474 = -t492 * t510 + t504 * t541;
t530 = t504 * t542;
t475 = t492 * t513 + t530;
t552 = Icges(4,5) * t475 + Icges(4,6) * t474 + t562 * t467 - t563 * t468 + t560 * t491;
t493 = t507 * t513 - t511 * t542;
t539 = t507 * t510;
t494 = t511 * t541 + t539;
t551 = Icges(4,5) * t494 + Icges(4,6) * t493 + t562 * t482 - t563 * t483 - t560 * t540;
t550 = qJD(2) ^ 2;
t546 = pkin(3) * t513;
t544 = t504 * t505;
t543 = t505 * t506;
t410 = -pkin(3) * t529 + qJ(4) * t489 + t490 * t546;
t424 = pkin(4) * t466 + qJ(5) * t465;
t536 = -t410 - t424;
t411 = pkin(3) * t530 + qJ(4) * t491 + t492 * t546;
t425 = pkin(4) * t468 + qJ(5) * t467;
t535 = -t411 - t425;
t443 = pkin(4) * t483 + qJ(5) * t482;
t457 = pkin(3) * t539 + (-qJ(4) * t514 + t511 * t546) * t505;
t534 = -t443 - t457;
t533 = qJD(2) * t505;
t501 = t504 * t533;
t476 = qJD(3) * t491 + t501;
t503 = qJD(2) * t507;
t528 = t506 * t533;
t460 = pkin(2) * t490 + pkin(8) * t489;
t461 = pkin(2) * t492 + pkin(8) * t491;
t527 = t460 * t501 + t461 * t528 + qJD(1);
t477 = qJD(3) * t489 - t528;
t496 = -qJD(3) * t540 + t503;
t495 = (pkin(2) * t511 - pkin(8) * t514) * t505;
t522 = t461 * t503 - t495 * t501;
t521 = qJD(4) * t489 + t496 * t411 + t522;
t520 = (-t460 * t507 - t495 * t543) * qJD(2);
t519 = -qJD(4) * t540 + t476 * t410 + t527;
t518 = qJD(5) * t465 + t496 * t425 + t521;
t517 = qJD(5) * t482 + t476 * t424 + t519;
t516 = qJD(4) * t491 + t477 * t457 + t520;
t515 = qJD(5) * t467 + t477 * t443 + t516;
t512 = cos(qJ(6));
t509 = sin(qJ(6));
t484 = t507 * rSges(3,3) + (rSges(3,1) * t511 + rSges(3,2) * t514) * t505;
t481 = Icges(3,5) * t507 + (Icges(3,1) * t511 + Icges(3,4) * t514) * t505;
t480 = Icges(3,6) * t507 + (Icges(3,4) * t511 + Icges(3,2) * t514) * t505;
t479 = Icges(3,3) * t507 + (Icges(3,5) * t511 + Icges(3,6) * t514) * t505;
t471 = -pkin(5) * t540 + t483 * pkin(9);
t470 = t482 * t509 - t512 * t540;
t469 = t482 * t512 + t509 * t540;
t464 = qJD(6) * t483 + t496;
t458 = t494 * rSges(4,1) + t493 * rSges(4,2) - rSges(4,3) * t540;
t456 = Icges(4,1) * t494 + Icges(4,4) * t493 - Icges(4,5) * t540;
t455 = Icges(4,4) * t494 + Icges(4,2) * t493 - Icges(4,6) * t540;
t451 = rSges(3,1) * t492 - rSges(3,2) * t491 + rSges(3,3) * t544;
t450 = rSges(3,1) * t490 - rSges(3,2) * t489 - rSges(3,3) * t543;
t449 = Icges(3,1) * t492 - Icges(3,4) * t491 + Icges(3,5) * t544;
t448 = Icges(3,1) * t490 - Icges(3,4) * t489 - Icges(3,5) * t543;
t447 = Icges(3,4) * t492 - Icges(3,2) * t491 + Icges(3,6) * t544;
t446 = Icges(3,4) * t490 - Icges(3,2) * t489 - Icges(3,6) * t543;
t445 = Icges(3,5) * t492 - Icges(3,6) * t491 + Icges(3,3) * t544;
t444 = Icges(3,5) * t490 - Icges(3,6) * t489 - Icges(3,3) * t543;
t442 = t483 * rSges(5,1) - t482 * rSges(5,2) - rSges(5,3) * t540;
t441 = -rSges(6,1) * t540 - t483 * rSges(6,2) + t482 * rSges(6,3);
t434 = pkin(5) * t491 + pkin(9) * t468;
t433 = pkin(5) * t489 + pkin(9) * t466;
t432 = t467 * t509 + t491 * t512;
t431 = t467 * t512 - t491 * t509;
t430 = t465 * t509 + t489 * t512;
t429 = t465 * t512 - t489 * t509;
t428 = qJD(6) * t466 + t477;
t427 = qJD(6) * t468 + t476;
t422 = (-t450 * t507 - t484 * t543) * qJD(2);
t421 = (t451 * t507 - t484 * t544) * qJD(2);
t420 = rSges(4,1) * t475 + rSges(4,2) * t474 + rSges(4,3) * t491;
t419 = rSges(4,1) * t473 + rSges(4,2) * t472 + rSges(4,3) * t489;
t418 = Icges(4,1) * t475 + Icges(4,4) * t474 + Icges(4,5) * t491;
t417 = Icges(4,1) * t473 + Icges(4,4) * t472 + Icges(4,5) * t489;
t416 = Icges(4,4) * t475 + Icges(4,2) * t474 + Icges(4,6) * t491;
t415 = Icges(4,4) * t473 + Icges(4,2) * t472 + Icges(4,6) * t489;
t409 = rSges(5,1) * t468 - rSges(5,2) * t467 + rSges(5,3) * t491;
t408 = rSges(5,1) * t466 - rSges(5,2) * t465 + rSges(5,3) * t489;
t407 = rSges(6,1) * t491 - rSges(6,2) * t468 + rSges(6,3) * t467;
t406 = rSges(6,1) * t489 - rSges(6,2) * t466 + rSges(6,3) * t465;
t393 = rSges(7,1) * t470 + rSges(7,2) * t469 + rSges(7,3) * t483;
t392 = Icges(7,1) * t470 + Icges(7,4) * t469 + Icges(7,5) * t483;
t391 = Icges(7,4) * t470 + Icges(7,2) * t469 + Icges(7,6) * t483;
t390 = Icges(7,5) * t470 + Icges(7,6) * t469 + Icges(7,3) * t483;
t387 = qJD(1) + (t450 * t504 + t451 * t506) * t533;
t385 = rSges(7,1) * t432 + rSges(7,2) * t431 + rSges(7,3) * t468;
t384 = rSges(7,1) * t430 + rSges(7,2) * t429 + rSges(7,3) * t466;
t383 = Icges(7,1) * t432 + Icges(7,4) * t431 + Icges(7,5) * t468;
t382 = Icges(7,1) * t430 + Icges(7,4) * t429 + Icges(7,5) * t466;
t381 = Icges(7,4) * t432 + Icges(7,2) * t431 + Icges(7,6) * t468;
t380 = Icges(7,4) * t430 + Icges(7,2) * t429 + Icges(7,6) * t466;
t379 = Icges(7,5) * t432 + Icges(7,6) * t431 + Icges(7,3) * t468;
t378 = Icges(7,5) * t430 + Icges(7,6) * t429 + Icges(7,3) * t466;
t377 = -t419 * t496 + t458 * t477 + t520;
t376 = t420 * t496 - t458 * t476 + t522;
t375 = t419 * t476 - t420 * t477 + t527;
t374 = t442 * t477 + (-t408 - t410) * t496 + t516;
t373 = t409 * t496 + (-t442 - t457) * t476 + t521;
t372 = t476 * t408 + (-t409 - t411) * t477 + t519;
t371 = t441 * t477 + (-t406 + t536) * t496 + t515;
t370 = t407 * t496 + (-t441 + t534) * t476 + t518;
t369 = t476 * t406 + (-t407 + t535) * t477 + t517;
t368 = -t384 * t464 + t393 * t428 + t471 * t477 + (-t433 + t536) * t496 + t515;
t367 = t385 * t464 - t393 * t427 + t434 * t496 + (-t471 + t534) * t476 + t518;
t366 = t427 * t384 - t428 * t385 + t476 * t433 + (-t434 + t535) * t477 + t517;
t1 = m(3) * (t387 ^ 2 + t421 ^ 2 + t422 ^ 2) / 0.2e1 + m(4) * (t375 ^ 2 + t376 ^ 2 + t377 ^ 2) / 0.2e1 + m(5) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(6) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + m(7) * (t366 ^ 2 + t367 ^ 2 + t368 ^ 2) / 0.2e1 + t427 * ((t379 * t468 + t381 * t431 + t383 * t432) * t427 + (t378 * t468 + t380 * t431 + t382 * t432) * t428 + (t390 * t468 + t391 * t431 + t392 * t432) * t464) / 0.2e1 + t428 * ((t379 * t466 + t381 * t429 + t383 * t430) * t427 + (t378 * t466 + t380 * t429 + t382 * t430) * t428 + (t390 * t466 + t391 * t429 + t392 * t430) * t464) / 0.2e1 + t464 * ((t379 * t483 + t381 * t469 + t383 * t470) * t427 + (t378 * t483 + t380 * t469 + t382 * t470) * t428 + (t390 * t483 + t391 * t469 + t392 * t470) * t464) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 - t550 * ((-t445 * t543 - t447 * t489 + t449 * t490) * t544 - (-t444 * t543 - t446 * t489 + t448 * t490) * t543 + (-t479 * t543 - t480 * t489 + t481 * t490) * t507) * t543 / 0.2e1 + (t507 * (t507 ^ 2 * t479 + (((t447 * t514 + t449 * t511) * t504 - (t446 * t514 + t448 * t511) * t506) * t505 + (-t444 * t506 + t445 * t504 + t480 * t514 + t481 * t511) * t507) * t505) + ((t445 * t544 - t447 * t491 + t449 * t492) * t544 - (t444 * t544 - t446 * t491 + t448 * t492) * t543 + (t479 * t544 - t480 * t491 + t481 * t492) * t507) * t544) * t550 / 0.2e1 + ((t455 * t474 + t456 * t475 + t467 * t555 + t468 * t554 + t491 * t551) * t496 + (t415 * t474 + t417 * t475 + t467 * t559 + t557 * t468 + t553 * t491) * t477 + (t416 * t474 + t418 * t475 + t467 * t558 + t468 * t556 + t491 * t552) * t476) * t476 / 0.2e1 + ((t455 * t472 + t456 * t473 + t465 * t555 + t466 * t554 + t489 * t551) * t496 + (t415 * t472 + t417 * t473 + t465 * t559 + t557 * t466 + t553 * t489) * t477 + (t416 * t472 + t418 * t473 + t465 * t558 + t466 * t556 + t489 * t552) * t476) * t477 / 0.2e1 + ((t493 * t455 + t494 * t456 + t482 * t555 + t483 * t554 - t540 * t551) * t496 + (t493 * t415 + t494 * t417 + t482 * t559 + t557 * t483 - t553 * t540) * t477 + (t493 * t416 + t494 * t418 + t482 * t558 + t483 * t556 - t552 * t540) * t476) * t496 / 0.2e1;
T  = t1;
