% Calculate kinetic energy for
% S6PRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR4_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:35:54
% EndTime: 2019-03-08 20:35:57
% DurationCPUTime: 2.81s
% Computational Cost: add. (3352->381), mult. (5747->594), div. (0->0), fcn. (6966->14), ass. (0->172)
t516 = sin(pkin(12));
t519 = cos(pkin(12));
t521 = cos(pkin(6));
t518 = sin(pkin(6));
t524 = sin(qJ(2));
t555 = t518 * t524;
t497 = -t516 * t555 + t519 * t521;
t558 = t516 * t521;
t498 = t519 * t555 + t558;
t526 = cos(qJ(2));
t554 = t518 * t526;
t450 = Icges(4,5) * t498 + Icges(4,6) * t497 - Icges(4,3) * t554;
t488 = Icges(3,6) * t521 + (Icges(3,4) * t524 + Icges(3,2) * t526) * t518;
t565 = t450 - t488;
t564 = qJD(2) ^ 2;
t563 = pkin(3) * t519;
t525 = cos(qJ(5));
t562 = pkin(5) * t525;
t517 = sin(pkin(11));
t520 = cos(pkin(11));
t552 = t521 * t526;
t499 = t517 * t524 - t520 * t552;
t523 = sin(qJ(5));
t560 = t499 * t523;
t501 = t517 * t552 + t520 * t524;
t559 = t501 * t523;
t557 = t517 * t518;
t556 = t518 * t520;
t553 = t521 * t524;
t502 = -t517 * t553 + t520 * t526;
t470 = pkin(2) * t502 + qJ(3) * t501;
t511 = qJD(2) * t521;
t550 = qJD(3) * t499 + t470 * t511;
t549 = qJD(2) * t518;
t508 = t517 * t549;
t484 = qJD(4) * t501 + t508;
t548 = qJD(3) * t526;
t546 = pkin(12) + qJ(4);
t545 = t516 * t557;
t544 = t516 * t556;
t543 = t523 * t554;
t512 = sin(t546);
t539 = cos(t546);
t534 = t518 * t539;
t476 = t502 * t512 - t517 * t534;
t435 = qJD(5) * t476 + t484;
t542 = t520 * t549;
t500 = t517 * t526 + t520 * t553;
t469 = pkin(2) * t500 + qJ(3) * t499;
t541 = t469 * t508 + t470 * t542 + qJD(1);
t503 = (pkin(2) * t524 - qJ(3) * t526) * t518;
t538 = (-rSges(4,1) * t498 - rSges(4,2) * t497 + rSges(4,3) * t554 - t503) * t518;
t537 = (-pkin(3) * t558 - (-pkin(8) * t526 + t524 * t563) * t518 - t503) * t518;
t485 = qJD(4) * t499 - t542;
t504 = -qJD(4) * t554 + t511;
t474 = t500 * t512 + t520 * t534;
t436 = qJD(5) * t474 + t485;
t490 = t512 * t555 - t521 * t539;
t473 = qJD(5) * t490 + t504;
t419 = -pkin(3) * t544 + pkin(8) * t499 + t500 * t563;
t420 = pkin(3) * t545 + pkin(8) * t501 + t502 * t563;
t533 = t419 * t508 + t420 * t542 - t518 * t548 + t541;
t532 = qJD(2) * t517 * t537 + t420 * t511 + t550;
t475 = t500 * t539 - t512 * t556;
t433 = t475 * pkin(4) + t474 * pkin(9);
t477 = t502 * t539 + t512 * t557;
t434 = t477 * pkin(4) + t476 * pkin(9);
t531 = t484 * t433 - t434 * t485 + t533;
t496 = qJD(3) * t501;
t530 = t496 + ((-t419 - t469) * t521 + t520 * t537) * qJD(2);
t491 = t521 * t512 + t524 * t534;
t460 = t491 * pkin(4) + t490 * pkin(9);
t529 = t504 * t434 - t460 * t484 + t532;
t528 = -t433 * t504 + t485 * t460 + t530;
t515 = qJ(5) + qJ(6);
t514 = cos(t515);
t513 = sin(t515);
t492 = rSges(3,3) * t521 + (rSges(3,1) * t524 + rSges(3,2) * t526) * t518;
t489 = Icges(3,5) * t521 + (Icges(3,1) * t524 + Icges(3,4) * t526) * t518;
t487 = Icges(3,3) * t521 + (Icges(3,5) * t524 + Icges(3,6) * t526) * t518;
t483 = t502 * t519 + t545;
t482 = -t502 * t516 + t519 * t557;
t481 = t500 * t519 - t544;
t480 = -t500 * t516 - t519 * t556;
t479 = t491 * t525 - t543;
t478 = -t491 * t523 - t525 * t554;
t468 = t491 * t514 - t513 * t554;
t467 = -t491 * t513 - t514 * t554;
t465 = rSges(3,1) * t502 - rSges(3,2) * t501 + rSges(3,3) * t557;
t464 = rSges(3,1) * t500 - rSges(3,2) * t499 - rSges(3,3) * t556;
t458 = Icges(3,1) * t502 - Icges(3,4) * t501 + Icges(3,5) * t557;
t457 = Icges(3,1) * t500 - Icges(3,4) * t499 - Icges(3,5) * t556;
t456 = Icges(3,4) * t502 - Icges(3,2) * t501 + Icges(3,6) * t557;
t455 = Icges(3,4) * t500 - Icges(3,2) * t499 - Icges(3,6) * t556;
t454 = Icges(3,5) * t502 - Icges(3,6) * t501 + Icges(3,3) * t557;
t453 = Icges(3,5) * t500 - Icges(3,6) * t499 - Icges(3,3) * t556;
t452 = Icges(4,1) * t498 + Icges(4,4) * t497 - Icges(4,5) * t554;
t451 = Icges(4,4) * t498 + Icges(4,2) * t497 - Icges(4,6) * t554;
t449 = rSges(5,1) * t491 - rSges(5,2) * t490 - rSges(5,3) * t554;
t448 = Icges(5,1) * t491 - Icges(5,4) * t490 - Icges(5,5) * t554;
t447 = Icges(5,4) * t491 - Icges(5,2) * t490 - Icges(5,6) * t554;
t446 = Icges(5,5) * t491 - Icges(5,6) * t490 - Icges(5,3) * t554;
t445 = t477 * t525 + t559;
t444 = -t477 * t523 + t501 * t525;
t443 = t475 * t525 + t560;
t442 = -t475 * t523 + t499 * t525;
t441 = qJD(6) * t490 + t473;
t440 = t477 * t514 + t501 * t513;
t439 = -t477 * t513 + t501 * t514;
t438 = t475 * t514 + t499 * t513;
t437 = -t475 * t513 + t499 * t514;
t431 = (-t464 * t521 - t492 * t556) * qJD(2);
t430 = (t465 * t521 - t492 * t557) * qJD(2);
t429 = rSges(4,1) * t483 + rSges(4,2) * t482 + rSges(4,3) * t501;
t428 = rSges(4,1) * t481 + rSges(4,2) * t480 + rSges(4,3) * t499;
t427 = Icges(4,1) * t483 + Icges(4,4) * t482 + Icges(4,5) * t501;
t426 = Icges(4,1) * t481 + Icges(4,4) * t480 + Icges(4,5) * t499;
t425 = Icges(4,4) * t483 + Icges(4,2) * t482 + Icges(4,6) * t501;
t424 = Icges(4,4) * t481 + Icges(4,2) * t480 + Icges(4,6) * t499;
t423 = Icges(4,5) * t483 + Icges(4,6) * t482 + Icges(4,3) * t501;
t422 = Icges(4,5) * t481 + Icges(4,6) * t480 + Icges(4,3) * t499;
t418 = rSges(5,1) * t477 - rSges(5,2) * t476 + rSges(5,3) * t501;
t417 = rSges(5,1) * t475 - rSges(5,2) * t474 + rSges(5,3) * t499;
t416 = Icges(5,1) * t477 - Icges(5,4) * t476 + Icges(5,5) * t501;
t415 = Icges(5,1) * t475 - Icges(5,4) * t474 + Icges(5,5) * t499;
t414 = Icges(5,4) * t477 - Icges(5,2) * t476 + Icges(5,6) * t501;
t413 = Icges(5,4) * t475 - Icges(5,2) * t474 + Icges(5,6) * t499;
t412 = Icges(5,5) * t477 - Icges(5,6) * t476 + Icges(5,3) * t501;
t411 = Icges(5,5) * t475 - Icges(5,6) * t474 + Icges(5,3) * t499;
t410 = rSges(6,1) * t479 + rSges(6,2) * t478 + rSges(6,3) * t490;
t409 = Icges(6,1) * t479 + Icges(6,4) * t478 + Icges(6,5) * t490;
t408 = Icges(6,4) * t479 + Icges(6,2) * t478 + Icges(6,6) * t490;
t407 = Icges(6,5) * t479 + Icges(6,6) * t478 + Icges(6,3) * t490;
t403 = qJD(6) * t474 + t436;
t402 = qJD(6) * t476 + t435;
t400 = rSges(7,1) * t468 + rSges(7,2) * t467 + rSges(7,3) * t490;
t399 = Icges(7,1) * t468 + Icges(7,4) * t467 + Icges(7,5) * t490;
t398 = Icges(7,4) * t468 + Icges(7,2) * t467 + Icges(7,6) * t490;
t397 = Icges(7,5) * t468 + Icges(7,6) * t467 + Icges(7,3) * t490;
t396 = -pkin(5) * t543 + pkin(10) * t490 + t491 * t562;
t395 = qJD(1) + (t464 * t517 + t465 * t520) * t549;
t394 = rSges(6,1) * t445 + rSges(6,2) * t444 + rSges(6,3) * t476;
t393 = rSges(6,1) * t443 + rSges(6,2) * t442 + rSges(6,3) * t474;
t392 = Icges(6,1) * t445 + Icges(6,4) * t444 + Icges(6,5) * t476;
t391 = Icges(6,1) * t443 + Icges(6,4) * t442 + Icges(6,5) * t474;
t390 = Icges(6,4) * t445 + Icges(6,2) * t444 + Icges(6,6) * t476;
t389 = Icges(6,4) * t443 + Icges(6,2) * t442 + Icges(6,6) * t474;
t388 = Icges(6,5) * t445 + Icges(6,6) * t444 + Icges(6,3) * t476;
t387 = Icges(6,5) * t443 + Icges(6,6) * t442 + Icges(6,3) * t474;
t386 = rSges(7,1) * t440 + rSges(7,2) * t439 + rSges(7,3) * t476;
t385 = rSges(7,1) * t438 + rSges(7,2) * t437 + rSges(7,3) * t474;
t384 = Icges(7,1) * t440 + Icges(7,4) * t439 + Icges(7,5) * t476;
t383 = Icges(7,1) * t438 + Icges(7,4) * t437 + Icges(7,5) * t474;
t382 = Icges(7,4) * t440 + Icges(7,2) * t439 + Icges(7,6) * t476;
t381 = Icges(7,4) * t438 + Icges(7,2) * t437 + Icges(7,6) * t474;
t380 = Icges(7,5) * t440 + Icges(7,6) * t439 + Icges(7,3) * t476;
t379 = Icges(7,5) * t438 + Icges(7,6) * t437 + Icges(7,3) * t474;
t378 = pkin(5) * t559 + pkin(10) * t476 + t477 * t562;
t377 = pkin(5) * t560 + pkin(10) * t474 + t475 * t562;
t376 = t496 + ((-t428 - t469) * t521 + t520 * t538) * qJD(2);
t375 = (t429 * t521 + t517 * t538) * qJD(2) + t550;
t374 = (-t548 + (t428 * t517 + t429 * t520) * qJD(2)) * t518 + t541;
t373 = -t417 * t504 + t449 * t485 + t530;
t372 = t418 * t504 - t449 * t484 + t532;
t371 = t417 * t484 - t418 * t485 + t533;
t370 = -t393 * t473 + t410 * t436 + t528;
t369 = t394 * t473 - t410 * t435 + t529;
t368 = t393 * t435 - t394 * t436 + t531;
t367 = -t377 * t473 - t385 * t441 + t396 * t436 + t403 * t400 + t528;
t366 = t378 * t473 + t386 * t441 - t396 * t435 - t402 * t400 + t529;
t365 = t377 * t435 - t378 * t436 + t402 * t385 - t403 * t386 + t531;
t1 = t441 * ((t380 * t490 + t382 * t467 + t384 * t468) * t402 + (t379 * t490 + t381 * t467 + t383 * t468) * t403 + (t490 * t397 + t467 * t398 + t468 * t399) * t441) / 0.2e1 + t403 * ((t380 * t474 + t382 * t437 + t384 * t438) * t402 + (t474 * t379 + t437 * t381 + t438 * t383) * t403 + (t397 * t474 + t398 * t437 + t399 * t438) * t441) / 0.2e1 + t402 * ((t476 * t380 + t439 * t382 + t440 * t384) * t402 + (t379 * t476 + t381 * t439 + t383 * t440) * t403 + (t397 * t476 + t398 * t439 + t399 * t440) * t441) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(5) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + m(7) * (t365 ^ 2 + t366 ^ 2 + t367 ^ 2) / 0.2e1 + m(6) * (t368 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + m(4) * (t374 ^ 2 + t375 ^ 2 + t376 ^ 2) / 0.2e1 + m(3) * (t395 ^ 2 + t430 ^ 2 + t431 ^ 2) / 0.2e1 + t436 * ((t388 * t474 + t390 * t442 + t392 * t443) * t435 + (t474 * t387 + t442 * t389 + t443 * t391) * t436 + (t407 * t474 + t408 * t442 + t409 * t443) * t473) / 0.2e1 + t473 * ((t388 * t490 + t390 * t478 + t392 * t479) * t435 + (t387 * t490 + t389 * t478 + t391 * t479) * t436 + (t407 * t490 + t408 * t478 + t409 * t479) * t473) / 0.2e1 + t435 * ((t476 * t388 + t444 * t390 + t445 * t392) * t435 + (t387 * t476 + t389 * t444 + t391 * t445) * t436 + (t407 * t476 + t408 * t444 + t409 * t445) * t473) / 0.2e1 + t504 * ((-t412 * t554 - t414 * t490 + t416 * t491) * t484 + (-t411 * t554 - t413 * t490 + t415 * t491) * t485 + (-t446 * t554 - t447 * t490 + t448 * t491) * t504) / 0.2e1 + t484 * ((t412 * t501 - t414 * t476 + t416 * t477) * t484 + (t411 * t501 - t413 * t476 + t415 * t477) * t485 + (t446 * t501 - t447 * t476 + t448 * t477) * t504) / 0.2e1 + t485 * ((t412 * t499 - t414 * t474 + t416 * t475) * t484 + (t411 * t499 - t413 * t474 + t415 * t475) * t485 + (t446 * t499 - t447 * t474 + t448 * t475) * t504) / 0.2e1 - (((t423 * t499 + t425 * t480 + t427 * t481) * t517 - (t422 * t499 + t424 * t480 + t426 * t481) * t520) * t518 + (-t454 * t556 - t456 * t499 + t458 * t500) * t557 - (-t453 * t556 - t455 * t499 + t457 * t500) * t556 + (t451 * t480 + t452 * t481 - t487 * t556 + t489 * t500 + t499 * t565) * t521) * t564 * t556 / 0.2e1 + (((-t423 * t554 + t425 * t497 + t427 * t498) * t557 - (-t422 * t554 + t424 * t497 + t426 * t498) * t556 + ((t456 * t526 + t458 * t524) * t517 - (t455 * t526 + t457 * t524) * t520) * t518 ^ 2 + (-t450 * t554 + t451 * t497 + t452 * t498 + (-t453 * t520 + t454 * t517 + t488 * t526 + t489 * t524) * t518 + t487 * t521) * t521) * t521 + (((t423 * t501 + t425 * t482 + t427 * t483) * t517 - (t422 * t501 + t424 * t482 + t426 * t483) * t520) * t518 + (t454 * t557 - t456 * t501 + t458 * t502) * t557 - (t453 * t557 - t455 * t501 + t457 * t502) * t556 + (t451 * t482 + t452 * t483 + t487 * t557 + t489 * t502 + t501 * t565) * t521) * t557) * t564 / 0.2e1;
T  = t1;
