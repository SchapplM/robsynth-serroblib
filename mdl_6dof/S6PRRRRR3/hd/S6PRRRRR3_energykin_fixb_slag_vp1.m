% Calculate kinetic energy for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:47:50
% EndTime: 2019-03-09 00:47:52
% DurationCPUTime: 2.12s
% Computational Cost: add. (3296->380), mult. (6897->598), div. (0->0), fcn. (8498->14), ass. (0->171)
t563 = qJD(2) ^ 2;
t562 = cos(qJ(3));
t531 = cos(qJ(4));
t560 = t531 * pkin(4);
t524 = sin(pkin(12));
t526 = cos(pkin(12));
t530 = sin(qJ(2));
t527 = cos(pkin(6));
t532 = cos(qJ(2));
t551 = t527 * t532;
t500 = t524 * t530 - t526 * t551;
t528 = sin(qJ(4));
t558 = t500 * t528;
t502 = t524 * t551 + t526 * t530;
t557 = t502 * t528;
t525 = sin(pkin(6));
t556 = t524 * t525;
t555 = t525 * t526;
t529 = sin(qJ(3));
t554 = t525 * t529;
t553 = t525 * t532;
t552 = t527 * t530;
t523 = qJ(4) + qJ(5);
t519 = cos(t523);
t550 = pkin(5) * t519;
t548 = qJD(2) * t525;
t512 = t524 * t548;
t490 = qJD(3) * t502 + t512;
t517 = qJD(2) * t527;
t546 = t528 * t553;
t503 = -t524 * t552 + t526 * t532;
t545 = t525 * t562;
t486 = t503 * t529 - t524 * t545;
t446 = qJD(4) * t486 + t490;
t544 = t526 * t548;
t501 = t524 * t532 + t526 * t552;
t472 = pkin(2) * t501 + pkin(8) * t500;
t473 = pkin(2) * t503 + pkin(8) * t502;
t543 = t472 * t512 + t473 * t544 + qJD(1);
t518 = sin(t523);
t542 = pkin(5) * t518;
t416 = qJD(5) * t486 + t446;
t491 = qJD(3) * t500 - t544;
t507 = -qJD(3) * t553 + t517;
t506 = (pkin(2) * t530 - pkin(8) * t532) * t525;
t541 = t473 * t517 - t506 * t512;
t484 = t501 * t529 + t526 * t545;
t447 = qJD(4) * t484 + t491;
t504 = -t527 * t562 + t530 * t554;
t483 = qJD(4) * t504 + t507;
t417 = qJD(5) * t484 + t447;
t456 = qJD(5) * t504 + t483;
t485 = t501 * t562 - t526 * t554;
t440 = t485 * pkin(3) + t484 * pkin(9);
t487 = t503 * t562 + t524 * t554;
t441 = t487 * pkin(3) + t486 * pkin(9);
t540 = t490 * t440 - t441 * t491 + t543;
t539 = (-t472 * t527 - t506 * t555) * qJD(2);
t505 = t527 * t529 + t530 * t545;
t474 = t505 * pkin(3) + t504 * pkin(9);
t538 = t507 * t441 - t474 * t490 + t541;
t388 = pkin(4) * t558 + pkin(10) * t484 + t485 * t560;
t389 = pkin(4) * t557 + pkin(10) * t486 + t487 * t560;
t537 = t446 * t388 - t389 * t447 + t540;
t536 = -t507 * t440 + t491 * t474 + t539;
t422 = -pkin(4) * t546 + pkin(10) * t504 + t505 * t560;
t535 = t483 * t389 - t422 * t446 + t538;
t534 = -t483 * t388 + t447 * t422 + t536;
t520 = qJ(6) + t523;
t515 = cos(t520);
t514 = sin(t520);
t495 = t527 * rSges(3,3) + (rSges(3,1) * t530 + rSges(3,2) * t532) * t525;
t494 = Icges(3,5) * t527 + (Icges(3,1) * t530 + Icges(3,4) * t532) * t525;
t493 = Icges(3,6) * t527 + (Icges(3,4) * t530 + Icges(3,2) * t532) * t525;
t492 = Icges(3,3) * t527 + (Icges(3,5) * t530 + Icges(3,6) * t532) * t525;
t489 = t505 * t531 - t546;
t488 = -t505 * t528 - t531 * t553;
t478 = t505 * t519 - t518 * t553;
t477 = -t505 * t518 - t519 * t553;
t476 = t505 * t515 - t514 * t553;
t475 = -t505 * t514 - t515 * t553;
t470 = rSges(4,1) * t505 - rSges(4,2) * t504 - rSges(4,3) * t553;
t469 = Icges(4,1) * t505 - Icges(4,4) * t504 - Icges(4,5) * t553;
t468 = Icges(4,4) * t505 - Icges(4,2) * t504 - Icges(4,6) * t553;
t467 = Icges(4,5) * t505 - Icges(4,6) * t504 - Icges(4,3) * t553;
t464 = rSges(3,1) * t503 - rSges(3,2) * t502 + rSges(3,3) * t556;
t463 = rSges(3,1) * t501 - rSges(3,2) * t500 - rSges(3,3) * t555;
t462 = Icges(3,1) * t503 - Icges(3,4) * t502 + Icges(3,5) * t556;
t461 = Icges(3,1) * t501 - Icges(3,4) * t500 - Icges(3,5) * t555;
t460 = Icges(3,4) * t503 - Icges(3,2) * t502 + Icges(3,6) * t556;
t459 = Icges(3,4) * t501 - Icges(3,2) * t500 - Icges(3,6) * t555;
t458 = Icges(3,5) * t503 - Icges(3,6) * t502 + Icges(3,3) * t556;
t457 = Icges(3,5) * t501 - Icges(3,6) * t500 - Icges(3,3) * t555;
t455 = t487 * t531 + t557;
t454 = -t487 * t528 + t502 * t531;
t453 = t485 * t531 + t558;
t452 = -t485 * t528 + t500 * t531;
t451 = t487 * t519 + t502 * t518;
t450 = -t487 * t518 + t502 * t519;
t449 = t485 * t519 + t500 * t518;
t448 = -t485 * t518 + t500 * t519;
t445 = t487 * t515 + t502 * t514;
t444 = -t487 * t514 + t502 * t515;
t443 = t485 * t515 + t500 * t514;
t442 = -t485 * t514 + t500 * t515;
t438 = qJD(6) * t504 + t456;
t436 = (-t463 * t527 - t495 * t555) * qJD(2);
t435 = (t464 * t527 - t495 * t556) * qJD(2);
t434 = rSges(5,1) * t489 + rSges(5,2) * t488 + rSges(5,3) * t504;
t433 = Icges(5,1) * t489 + Icges(5,4) * t488 + Icges(5,5) * t504;
t432 = Icges(5,4) * t489 + Icges(5,2) * t488 + Icges(5,6) * t504;
t431 = Icges(5,5) * t489 + Icges(5,6) * t488 + Icges(5,3) * t504;
t430 = rSges(4,1) * t487 - rSges(4,2) * t486 + rSges(4,3) * t502;
t429 = rSges(4,1) * t485 - rSges(4,2) * t484 + rSges(4,3) * t500;
t428 = Icges(4,1) * t487 - Icges(4,4) * t486 + Icges(4,5) * t502;
t427 = Icges(4,1) * t485 - Icges(4,4) * t484 + Icges(4,5) * t500;
t426 = Icges(4,4) * t487 - Icges(4,2) * t486 + Icges(4,6) * t502;
t425 = Icges(4,4) * t485 - Icges(4,2) * t484 + Icges(4,6) * t500;
t424 = Icges(4,5) * t487 - Icges(4,6) * t486 + Icges(4,3) * t502;
t423 = Icges(4,5) * t485 - Icges(4,6) * t484 + Icges(4,3) * t500;
t421 = rSges(6,1) * t478 + rSges(6,2) * t477 + rSges(6,3) * t504;
t420 = Icges(6,1) * t478 + Icges(6,4) * t477 + Icges(6,5) * t504;
t419 = Icges(6,4) * t478 + Icges(6,2) * t477 + Icges(6,6) * t504;
t418 = Icges(6,5) * t478 + Icges(6,6) * t477 + Icges(6,3) * t504;
t414 = rSges(7,1) * t476 + rSges(7,2) * t475 + rSges(7,3) * t504;
t413 = Icges(7,1) * t476 + Icges(7,4) * t475 + Icges(7,5) * t504;
t412 = Icges(7,4) * t476 + Icges(7,2) * t475 + Icges(7,6) * t504;
t411 = Icges(7,5) * t476 + Icges(7,6) * t475 + Icges(7,3) * t504;
t410 = qJD(1) + (t463 * t524 + t464 * t526) * t548;
t409 = pkin(11) * t504 + t505 * t550 - t542 * t553;
t408 = qJD(6) * t484 + t417;
t407 = qJD(6) * t486 + t416;
t406 = rSges(5,1) * t455 + rSges(5,2) * t454 + rSges(5,3) * t486;
t405 = rSges(5,1) * t453 + rSges(5,2) * t452 + rSges(5,3) * t484;
t404 = Icges(5,1) * t455 + Icges(5,4) * t454 + Icges(5,5) * t486;
t403 = Icges(5,1) * t453 + Icges(5,4) * t452 + Icges(5,5) * t484;
t402 = Icges(5,4) * t455 + Icges(5,2) * t454 + Icges(5,6) * t486;
t401 = Icges(5,4) * t453 + Icges(5,2) * t452 + Icges(5,6) * t484;
t400 = Icges(5,5) * t455 + Icges(5,6) * t454 + Icges(5,3) * t486;
t399 = Icges(5,5) * t453 + Icges(5,6) * t452 + Icges(5,3) * t484;
t397 = rSges(6,1) * t451 + rSges(6,2) * t450 + rSges(6,3) * t486;
t396 = rSges(6,1) * t449 + rSges(6,2) * t448 + rSges(6,3) * t484;
t395 = Icges(6,1) * t451 + Icges(6,4) * t450 + Icges(6,5) * t486;
t394 = Icges(6,1) * t449 + Icges(6,4) * t448 + Icges(6,5) * t484;
t393 = Icges(6,4) * t451 + Icges(6,2) * t450 + Icges(6,6) * t486;
t392 = Icges(6,4) * t449 + Icges(6,2) * t448 + Icges(6,6) * t484;
t391 = Icges(6,5) * t451 + Icges(6,6) * t450 + Icges(6,3) * t486;
t390 = Icges(6,5) * t449 + Icges(6,6) * t448 + Icges(6,3) * t484;
t387 = rSges(7,1) * t445 + rSges(7,2) * t444 + rSges(7,3) * t486;
t386 = rSges(7,1) * t443 + rSges(7,2) * t442 + rSges(7,3) * t484;
t385 = Icges(7,1) * t445 + Icges(7,4) * t444 + Icges(7,5) * t486;
t384 = Icges(7,1) * t443 + Icges(7,4) * t442 + Icges(7,5) * t484;
t383 = Icges(7,4) * t445 + Icges(7,2) * t444 + Icges(7,6) * t486;
t382 = Icges(7,4) * t443 + Icges(7,2) * t442 + Icges(7,6) * t484;
t381 = Icges(7,5) * t445 + Icges(7,6) * t444 + Icges(7,3) * t486;
t380 = Icges(7,5) * t443 + Icges(7,6) * t442 + Icges(7,3) * t484;
t378 = pkin(11) * t486 + t487 * t550 + t502 * t542;
t377 = pkin(11) * t484 + t485 * t550 + t500 * t542;
t375 = -t507 * t429 + t491 * t470 + t539;
t374 = t430 * t507 - t470 * t490 + t541;
t373 = t429 * t490 - t430 * t491 + t543;
t372 = -t483 * t405 + t447 * t434 + t536;
t371 = t406 * t483 - t434 * t446 + t538;
t370 = t405 * t446 - t406 * t447 + t540;
t369 = -t456 * t396 + t417 * t421 + t534;
t368 = t397 * t456 - t416 * t421 + t535;
t367 = t396 * t416 - t397 * t417 + t537;
t366 = -t456 * t377 - t438 * t386 + t408 * t414 + t417 * t409 + t534;
t365 = t378 * t456 + t387 * t438 - t407 * t414 - t409 * t416 + t535;
t364 = t377 * t416 - t378 * t417 + t386 * t407 - t387 * t408 + t537;
t1 = m(6) * (t367 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + m(7) * (t364 ^ 2 + t365 ^ 2 + t366 ^ 2) / 0.2e1 + t490 * ((t502 * t424 - t486 * t426 + t487 * t428) * t490 + (t423 * t502 - t425 * t486 + t427 * t487) * t491 + (t467 * t502 - t468 * t486 + t469 * t487) * t507) / 0.2e1 + t491 * ((t424 * t500 - t426 * t484 + t428 * t485) * t490 + (t500 * t423 - t484 * t425 + t485 * t427) * t491 + (t467 * t500 - t468 * t484 + t469 * t485) * t507) / 0.2e1 + t507 * ((-t424 * t553 - t426 * t504 + t428 * t505) * t490 + (-t423 * t553 - t425 * t504 + t427 * t505) * t491 + (-t467 * t553 - t504 * t468 + t505 * t469) * t507) / 0.2e1 + t446 * ((t486 * t400 + t454 * t402 + t455 * t404) * t446 + (t399 * t486 + t401 * t454 + t403 * t455) * t447 + (t431 * t486 + t432 * t454 + t433 * t455) * t483) / 0.2e1 + t447 * ((t400 * t484 + t402 * t452 + t404 * t453) * t446 + (t484 * t399 + t452 * t401 + t453 * t403) * t447 + (t431 * t484 + t432 * t452 + t433 * t453) * t483) / 0.2e1 + t483 * ((t400 * t504 + t402 * t488 + t404 * t489) * t446 + (t399 * t504 + t401 * t488 + t403 * t489) * t447 + (t504 * t431 + t488 * t432 + t489 * t433) * t483) / 0.2e1 + t416 * ((t486 * t391 + t450 * t393 + t451 * t395) * t416 + (t390 * t486 + t392 * t450 + t394 * t451) * t417 + (t418 * t486 + t419 * t450 + t420 * t451) * t456) / 0.2e1 + t417 * ((t391 * t484 + t393 * t448 + t395 * t449) * t416 + (t484 * t390 + t448 * t392 + t449 * t394) * t417 + (t418 * t484 + t419 * t448 + t420 * t449) * t456) / 0.2e1 + t456 * ((t391 * t504 + t393 * t477 + t395 * t478) * t416 + (t390 * t504 + t392 * t477 + t394 * t478) * t417 + (t504 * t418 + t477 * t419 + t478 * t420) * t456) / 0.2e1 + t407 * ((t486 * t381 + t444 * t383 + t445 * t385) * t407 + (t380 * t486 + t382 * t444 + t384 * t445) * t408 + (t411 * t486 + t412 * t444 + t413 * t445) * t438) / 0.2e1 + t408 * ((t381 * t484 + t383 * t442 + t385 * t443) * t407 + (t484 * t380 + t442 * t382 + t443 * t384) * t408 + (t411 * t484 + t412 * t442 + t413 * t443) * t438) / 0.2e1 + t438 * ((t381 * t504 + t383 * t475 + t385 * t476) * t407 + (t380 * t504 + t382 * t475 + t384 * t476) * t408 + (t504 * t411 + t475 * t412 + t476 * t413) * t438) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t410 ^ 2 + t435 ^ 2 + t436 ^ 2) / 0.2e1 + m(4) * (t373 ^ 2 + t374 ^ 2 + t375 ^ 2) / 0.2e1 + m(5) * (t370 ^ 2 + t371 ^ 2 + t372 ^ 2) / 0.2e1 - t563 * ((-t458 * t555 - t460 * t500 + t462 * t501) * t556 - (-t457 * t555 - t459 * t500 + t461 * t501) * t555 + (-t492 * t555 - t493 * t500 + t494 * t501) * t527) * t555 / 0.2e1 + (t527 * (t527 ^ 2 * t492 + (((t460 * t532 + t462 * t530) * t524 - (t459 * t532 + t461 * t530) * t526) * t525 + (-t457 * t526 + t458 * t524 + t493 * t532 + t494 * t530) * t527) * t525) + ((t458 * t556 - t460 * t502 + t462 * t503) * t556 - (t457 * t556 - t459 * t502 + t461 * t503) * t555 + (t492 * t556 - t493 * t502 + t494 * t503) * t527) * t556) * t563 / 0.2e1;
T  = t1;
