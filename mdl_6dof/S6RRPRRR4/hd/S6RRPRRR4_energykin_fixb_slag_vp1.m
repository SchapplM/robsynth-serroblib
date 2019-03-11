% Calculate kinetic energy for
% S6RRPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR4_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:27:46
% EndTime: 2019-03-09 13:27:49
% DurationCPUTime: 3.49s
% Computational Cost: add. (3723->385), mult. (8045->593), div. (0->0), fcn. (10104->14), ass. (0->175)
t553 = sin(qJ(2));
t595 = sin(pkin(12));
t596 = cos(pkin(12));
t601 = cos(qJ(2));
t533 = -t553 * t596 - t595 * t601;
t554 = sin(qJ(1));
t557 = cos(qJ(1));
t597 = cos(pkin(6));
t569 = t597 * t595;
t570 = t597 * t596;
t563 = -t553 * t569 + t570 * t601;
t508 = t554 * t533 + t557 * t563;
t526 = t553 * t570 + t569 * t601;
t534 = -t553 * t595 + t601 * t596;
t509 = t526 * t557 + t534 * t554;
t550 = sin(pkin(6));
t593 = t550 * t557;
t454 = Icges(4,5) * t509 + Icges(4,6) * t508 - Icges(4,3) * t593;
t510 = t533 * t557 - t554 * t563;
t511 = -t526 * t554 + t534 * t557;
t594 = t550 * t554;
t455 = Icges(4,5) * t511 + Icges(4,6) * t510 + Icges(4,3) * t594;
t573 = t597 * t601;
t528 = -t554 * t553 + t557 * t573;
t581 = t553 * t597;
t529 = t554 * t601 + t557 * t581;
t493 = Icges(3,5) * t529 + Icges(3,6) * t528 - Icges(3,3) * t593;
t530 = -t557 * t553 - t554 * t573;
t531 = -t554 * t581 + t557 * t601;
t494 = Icges(3,5) * t531 + Icges(3,6) * t530 + Icges(3,3) * t594;
t606 = t550 * ((t454 + t493) * t557 + (-t455 - t494) * t554);
t524 = t534 * t550;
t525 = t533 * t550;
t605 = -Icges(4,5) * t525 + Icges(4,6) * t524 + (Icges(3,5) * t553 + Icges(3,6) * t601) * t550 + (Icges(4,3) + Icges(3,3)) * t597;
t600 = pkin(2) * t601;
t556 = cos(qJ(4));
t599 = pkin(4) * t556;
t592 = qJ(4) + qJ(5);
t582 = pkin(2) * t581 - qJ(3) * t550;
t507 = -t554 * t582 + t557 * t600;
t532 = qJD(1) * (pkin(1) * t557 + pkin(8) * t594);
t544 = qJD(2) * t597 + qJD(1);
t591 = t544 * t507 + t532;
t589 = qJD(2) * t550;
t543 = t554 * t589;
t483 = -qJD(4) * t510 + t543;
t590 = qJD(1) * (pkin(1) * t554 - pkin(8) * t593);
t588 = qJD(3) * t557;
t552 = sin(qJ(4));
t587 = t552 * t594;
t586 = t552 * t593;
t506 = t554 * t600 + t557 * t582;
t584 = t557 * t589;
t585 = qJD(3) * t597 + t506 * t543 + t507 * t584;
t451 = -qJD(5) * t510 + t483;
t516 = -qJD(4) * t524 + t544;
t583 = cos(t592);
t579 = t597 * t552;
t535 = t550 * t553 * pkin(2) + qJ(3) * t597;
t578 = qJD(2) * (t525 * rSges(4,1) - t524 * rSges(4,2) - rSges(4,3) * t597 - t535);
t577 = qJD(2) * (t525 * pkin(3) + t524 * pkin(9) - t535);
t576 = qJD(3) * t594 - t590;
t485 = -qJD(5) * t524 + t516;
t571 = t550 * t583;
t484 = -qJD(4) * t508 - t584;
t470 = pkin(3) * t509 - pkin(9) * t508;
t471 = pkin(3) * t511 - pkin(9) * t510;
t568 = t470 * t543 + t471 * t584 + t585;
t452 = -qJD(5) * t508 + t484;
t413 = -pkin(4) * t586 - pkin(10) * t508 + t509 * t599;
t414 = pkin(4) * t587 - pkin(10) * t510 + t511 * t599;
t565 = t483 * t413 - t414 * t484 + t568;
t564 = t544 * t471 + (t554 * t577 - t588) * t550 + t591;
t562 = (-t470 - t506) * t544 + t577 * t593 + t576;
t445 = pkin(4) * t579 - pkin(10) * t524 - t525 * t599;
t561 = t516 * t414 - t445 * t483 + t564;
t560 = -t413 * t516 + t484 * t445 + t562;
t555 = cos(qJ(6));
t551 = sin(qJ(6));
t549 = sin(t592);
t539 = rSges(2,1) * t557 - rSges(2,2) * t554;
t538 = rSges(2,1) * t554 + rSges(2,2) * t557;
t523 = t597 * rSges(3,3) + (rSges(3,1) * t553 + rSges(3,2) * t601) * t550;
t522 = Icges(3,5) * t597 + (Icges(3,1) * t553 + Icges(3,4) * t601) * t550;
t521 = Icges(3,6) * t597 + (Icges(3,4) * t553 + Icges(3,2) * t601) * t550;
t515 = -t525 * t556 + t579;
t514 = t525 * t552 + t556 * t597;
t513 = -t525 * t583 + t549 * t597;
t512 = -t525 * t549 - t583 * t597;
t500 = rSges(3,1) * t531 + rSges(3,2) * t530 + rSges(3,3) * t594;
t499 = rSges(3,1) * t529 + rSges(3,2) * t528 - rSges(3,3) * t593;
t498 = Icges(3,1) * t531 + Icges(3,4) * t530 + Icges(3,5) * t594;
t497 = Icges(3,1) * t529 + Icges(3,4) * t528 - Icges(3,5) * t593;
t496 = Icges(3,4) * t531 + Icges(3,2) * t530 + Icges(3,6) * t594;
t495 = Icges(3,4) * t529 + Icges(3,2) * t528 - Icges(3,6) * t593;
t491 = -Icges(4,1) * t525 + Icges(4,4) * t524 + Icges(4,5) * t597;
t490 = -Icges(4,4) * t525 + Icges(4,2) * t524 + Icges(4,6) * t597;
t482 = t511 * t556 + t587;
t481 = -t511 * t552 + t556 * t594;
t480 = t509 * t556 - t586;
t479 = -t509 * t552 - t556 * t593;
t478 = t511 * t583 + t549 * t594;
t477 = t511 * t549 - t554 * t571;
t476 = t509 * t583 - t549 * t593;
t475 = t509 * t549 + t557 * t571;
t474 = t513 * t555 - t524 * t551;
t473 = -t513 * t551 - t524 * t555;
t472 = pkin(5) * t513 + pkin(11) * t512;
t469 = rSges(5,1) * t515 + rSges(5,2) * t514 - rSges(5,3) * t524;
t468 = Icges(5,1) * t515 + Icges(5,4) * t514 - Icges(5,5) * t524;
t467 = Icges(5,4) * t515 + Icges(5,2) * t514 - Icges(5,6) * t524;
t466 = Icges(5,5) * t515 + Icges(5,6) * t514 - Icges(5,3) * t524;
t465 = qJD(6) * t512 + t485;
t462 = rSges(4,1) * t511 + rSges(4,2) * t510 + rSges(4,3) * t594;
t461 = rSges(4,1) * t509 + rSges(4,2) * t508 - rSges(4,3) * t593;
t459 = Icges(4,1) * t511 + Icges(4,4) * t510 + Icges(4,5) * t594;
t458 = Icges(4,1) * t509 + Icges(4,4) * t508 - Icges(4,5) * t593;
t457 = Icges(4,4) * t511 + Icges(4,2) * t510 + Icges(4,6) * t594;
t456 = Icges(4,4) * t509 + Icges(4,2) * t508 - Icges(4,6) * t593;
t453 = rSges(6,1) * t513 - rSges(6,2) * t512 - rSges(6,3) * t524;
t450 = Icges(6,1) * t513 - Icges(6,4) * t512 - Icges(6,5) * t524;
t449 = Icges(6,4) * t513 - Icges(6,2) * t512 - Icges(6,6) * t524;
t448 = Icges(6,5) * t513 - Icges(6,6) * t512 - Icges(6,3) * t524;
t447 = t500 * t544 - t523 * t543 + t532;
t446 = -t499 * t544 - t523 * t584 - t590;
t444 = (t499 * t554 + t500 * t557) * t589;
t443 = t478 * t555 - t510 * t551;
t442 = -t478 * t551 - t510 * t555;
t441 = t476 * t555 - t508 * t551;
t440 = -t476 * t551 - t508 * t555;
t439 = pkin(5) * t478 + pkin(11) * t477;
t438 = pkin(5) * t476 + pkin(11) * t475;
t436 = qJD(6) * t475 + t452;
t435 = qJD(6) * t477 + t451;
t434 = rSges(5,1) * t482 + rSges(5,2) * t481 - rSges(5,3) * t510;
t433 = rSges(5,1) * t480 + rSges(5,2) * t479 - rSges(5,3) * t508;
t432 = Icges(5,1) * t482 + Icges(5,4) * t481 - Icges(5,5) * t510;
t431 = Icges(5,1) * t480 + Icges(5,4) * t479 - Icges(5,5) * t508;
t430 = Icges(5,4) * t482 + Icges(5,2) * t481 - Icges(5,6) * t510;
t429 = Icges(5,4) * t480 + Icges(5,2) * t479 - Icges(5,6) * t508;
t428 = Icges(5,5) * t482 + Icges(5,6) * t481 - Icges(5,3) * t510;
t427 = Icges(5,5) * t480 + Icges(5,6) * t479 - Icges(5,3) * t508;
t426 = rSges(6,1) * t478 - rSges(6,2) * t477 - rSges(6,3) * t510;
t425 = rSges(6,1) * t476 - rSges(6,2) * t475 - rSges(6,3) * t508;
t424 = rSges(7,1) * t474 + rSges(7,2) * t473 + rSges(7,3) * t512;
t423 = Icges(6,1) * t478 - Icges(6,4) * t477 - Icges(6,5) * t510;
t422 = Icges(6,1) * t476 - Icges(6,4) * t475 - Icges(6,5) * t508;
t421 = Icges(6,4) * t478 - Icges(6,2) * t477 - Icges(6,6) * t510;
t420 = Icges(6,4) * t476 - Icges(6,2) * t475 - Icges(6,6) * t508;
t419 = Icges(6,5) * t478 - Icges(6,6) * t477 - Icges(6,3) * t510;
t418 = Icges(6,5) * t476 - Icges(6,6) * t475 - Icges(6,3) * t508;
t417 = Icges(7,1) * t474 + Icges(7,4) * t473 + Icges(7,5) * t512;
t416 = Icges(7,4) * t474 + Icges(7,2) * t473 + Icges(7,6) * t512;
t415 = Icges(7,5) * t474 + Icges(7,6) * t473 + Icges(7,3) * t512;
t410 = t462 * t544 + (t554 * t578 - t588) * t550 + t591;
t409 = (-t461 - t506) * t544 + t578 * t593 + t576;
t408 = rSges(7,1) * t443 + rSges(7,2) * t442 + rSges(7,3) * t477;
t407 = rSges(7,1) * t441 + rSges(7,2) * t440 + rSges(7,3) * t475;
t406 = Icges(7,1) * t443 + Icges(7,4) * t442 + Icges(7,5) * t477;
t405 = Icges(7,1) * t441 + Icges(7,4) * t440 + Icges(7,5) * t475;
t404 = Icges(7,4) * t443 + Icges(7,2) * t442 + Icges(7,6) * t477;
t403 = Icges(7,4) * t441 + Icges(7,2) * t440 + Icges(7,6) * t475;
t402 = Icges(7,5) * t443 + Icges(7,6) * t442 + Icges(7,3) * t477;
t401 = Icges(7,5) * t441 + Icges(7,6) * t440 + Icges(7,3) * t475;
t400 = (t461 * t554 + t462 * t557) * t589 + t585;
t399 = t434 * t516 - t469 * t483 + t564;
t398 = -t433 * t516 + t469 * t484 + t562;
t397 = t433 * t483 - t434 * t484 + t568;
t396 = t426 * t485 - t451 * t453 + t561;
t395 = -t425 * t485 + t452 * t453 + t560;
t394 = t425 * t451 - t426 * t452 + t565;
t393 = t408 * t465 - t424 * t435 + t439 * t485 - t451 * t472 + t561;
t392 = -t407 * t465 + t424 * t436 - t438 * t485 + t452 * t472 + t560;
t391 = t407 * t435 - t408 * t436 + t438 * t451 - t439 * t452 + t565;
t1 = t436 * ((t402 * t475 + t404 * t440 + t406 * t441) * t435 + (t401 * t475 + t403 * t440 + t405 * t441) * t436 + (t415 * t475 + t416 * t440 + t417 * t441) * t465) / 0.2e1 + t465 * ((t402 * t512 + t404 * t473 + t406 * t474) * t435 + (t401 * t512 + t403 * t473 + t405 * t474) * t436 + (t415 * t512 + t416 * t473 + t417 * t474) * t465) / 0.2e1 + t485 * ((-t419 * t524 - t421 * t512 + t423 * t513) * t451 + (-t418 * t524 - t420 * t512 + t422 * t513) * t452 + (-t448 * t524 - t449 * t512 + t450 * t513) * t485) / 0.2e1 + t452 * ((-t419 * t508 - t421 * t475 + t423 * t476) * t451 + (-t418 * t508 - t420 * t475 + t422 * t476) * t452 + (-t448 * t508 - t449 * t475 + t450 * t476) * t485) / 0.2e1 + t451 * ((-t419 * t510 - t421 * t477 + t423 * t478) * t451 + (-t418 * t510 - t420 * t477 + t422 * t478) * t452 + (-t448 * t510 - t449 * t477 + t450 * t478) * t485) / 0.2e1 + t483 * ((-t428 * t510 + t430 * t481 + t432 * t482) * t483 + (-t427 * t510 + t429 * t481 + t431 * t482) * t484 + (-t466 * t510 + t467 * t481 + t468 * t482) * t516) / 0.2e1 + t484 * ((-t428 * t508 + t430 * t479 + t432 * t480) * t483 + (-t427 * t508 + t429 * t479 + t431 * t480) * t484 + (-t466 * t508 + t467 * t479 + t468 * t480) * t516) / 0.2e1 + t516 * ((-t428 * t524 + t430 * t514 + t432 * t515) * t483 + (-t427 * t524 + t429 * t514 + t431 * t515) * t484 + (-t466 * t524 + t467 * t514 + t468 * t515) * t516) / 0.2e1 + m(7) * (t391 ^ 2 + t392 ^ 2 + t393 ^ 2) / 0.2e1 + m(6) * (t394 ^ 2 + t395 ^ 2 + t396 ^ 2) / 0.2e1 + m(5) * (t397 ^ 2 + t398 ^ 2 + t399 ^ 2) / 0.2e1 + m(4) * (t400 ^ 2 + t409 ^ 2 + t410 ^ 2) / 0.2e1 + m(3) * (t444 ^ 2 + t446 ^ 2 + t447 ^ 2) / 0.2e1 + t435 * ((t402 * t477 + t404 * t442 + t406 * t443) * t435 + (t401 * t477 + t403 * t442 + t405 * t443) * t436 + (t415 * t477 + t416 * t442 + t417 * t443) * t465) / 0.2e1 + ((t597 * t494 + (t496 * t601 + t498 * t553) * t550) * t543 - (t597 * t493 + (t495 * t601 + t497 * t553) * t550) * t584 + ((t455 * t597 + t524 * t457 - t525 * t459) * t554 - (t454 * t597 + t524 * t456 - t525 * t458) * t557) * t589 + ((t521 * t601 + t522 * t553) * t550 + t524 * t490 - t525 * t491 + t605 * t597) * t544) * t544 / 0.2e1 + (m(2) * (t538 ^ 2 + t539 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t456 * t510 - t458 * t511 - t495 * t530 - t497 * t531) * t557 + (t457 * t510 + t459 * t511 + t496 * t530 + t498 * t531 - t606) * t554) * t589 + (t490 * t510 + t491 * t511 + t521 * t530 + t522 * t531 + t594 * t605) * t544) * t543 / 0.2e1 - (((-t456 * t508 - t458 * t509 - t495 * t528 - t497 * t529 + t606) * t557 + (t457 * t508 + t459 * t509 + t496 * t528 + t498 * t529) * t554) * t589 + (t490 * t508 + t491 * t509 + t521 * t528 + t522 * t529 - t593 * t605) * t544) * t584 / 0.2e1;
T  = t1;
