% Calculate kinetic energy for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:54:30
% EndTime: 2019-03-09 08:54:34
% DurationCPUTime: 3.57s
% Computational Cost: add. (3513->379), mult. (7562->566), div. (0->0), fcn. (9474->14), ass. (0->171)
t615 = -Icges(4,2) - Icges(5,3);
t557 = sin(qJ(2));
t601 = sin(pkin(11));
t602 = cos(pkin(11));
t606 = cos(qJ(2));
t535 = -t557 * t602 - t601 * t606;
t558 = sin(qJ(1));
t560 = cos(qJ(1));
t603 = cos(pkin(6));
t570 = t603 * t601;
t571 = t603 * t602;
t564 = -t557 * t570 + t571 * t606;
t510 = t558 * t535 + t560 * t564;
t528 = t557 * t571 + t570 * t606;
t536 = -t557 * t601 + t606 * t602;
t511 = t528 * t560 + t536 * t558;
t553 = sin(pkin(6));
t599 = t553 * t560;
t457 = Icges(4,5) * t511 + Icges(4,6) * t510 - Icges(4,3) * t599;
t512 = t535 * t560 - t558 * t564;
t513 = -t528 * t558 + t536 * t560;
t600 = t553 * t558;
t458 = Icges(4,5) * t513 + Icges(4,6) * t512 + Icges(4,3) * t600;
t576 = t603 * t606;
t530 = -t558 * t557 + t560 * t576;
t585 = t557 * t603;
t531 = t558 * t606 + t560 * t585;
t495 = Icges(3,5) * t531 + Icges(3,6) * t530 - Icges(3,3) * t599;
t532 = -t560 * t557 - t558 * t576;
t533 = -t558 * t585 + t560 * t606;
t496 = Icges(3,5) * t533 + Icges(3,6) * t532 + Icges(3,3) * t600;
t614 = ((t457 + t495) * t560 + (-t458 - t496) * t558) * t553;
t552 = sin(pkin(12));
t554 = cos(pkin(12));
t482 = -t511 * t552 - t554 * t599;
t589 = t552 * t599;
t483 = t511 * t554 - t589;
t613 = -Icges(4,4) * t511 + Icges(5,5) * t483 + Icges(4,6) * t599 + Icges(5,6) * t482 + t615 * t510;
t484 = -t513 * t552 + t554 * t600;
t590 = t552 * t600;
t485 = t513 * t554 + t590;
t612 = -Icges(4,4) * t513 + Icges(5,5) * t485 - Icges(4,6) * t600 + Icges(5,6) * t484 + t615 * t512;
t527 = t535 * t553;
t516 = t527 * t552 + t554 * t603;
t583 = t603 * t552;
t517 = -t527 * t554 + t583;
t526 = t536 * t553;
t611 = Icges(4,4) * t527 + Icges(5,5) * t517 - Icges(4,6) * t603 + Icges(5,6) * t516 + t615 * t526;
t610 = -Icges(4,5) * t527 + Icges(4,6) * t526 + (Icges(3,5) * t557 + Icges(3,6) * t606) * t553 + (Icges(4,3) + Icges(3,3)) * t603;
t605 = pkin(2) * t606;
t604 = pkin(4) * t554;
t472 = pkin(3) * t511 - qJ(4) * t510;
t586 = pkin(2) * t585 - qJ(3) * t553;
t508 = t558 * t605 + t560 * t586;
t597 = -t472 - t508;
t509 = -t558 * t586 + t560 * t605;
t534 = qJD(1) * (pkin(1) * t560 + pkin(8) * t600);
t546 = qJD(2) * t603 + qJD(1);
t596 = t546 * t509 + t534;
t537 = t553 * t557 * pkin(2) + qJ(3) * t603;
t595 = t527 * pkin(3) + t526 * qJ(4) - t537;
t593 = qJD(2) * t553;
t545 = t558 * t593;
t486 = -qJD(5) * t512 + t545;
t594 = qJD(1) * (pkin(1) * t558 - pkin(8) * t599);
t592 = qJD(3) * t560;
t591 = pkin(12) + qJ(5);
t587 = t560 * t593;
t588 = qJD(3) * t603 + t508 * t545 + t509 * t587;
t518 = -qJD(5) * t526 + t546;
t582 = cos(t591);
t581 = qJD(2) * (t527 * rSges(4,1) - t526 * rSges(4,2) - rSges(4,3) * t603 - t537);
t580 = qJD(3) * t600 - t594;
t473 = pkin(3) * t513 - qJ(4) * t512;
t579 = -qJD(4) * t510 + t546 * t473 + t596;
t574 = qJD(2) * (-pkin(4) * t583 + pkin(9) * t526 + t527 * t604 + t595);
t573 = qJD(2) * (-rSges(5,1) * t517 - rSges(5,2) * t516 + rSges(5,3) * t526 + t595);
t487 = -qJD(5) * t510 - t587;
t572 = -qJD(4) * t512 + t580;
t569 = t553 * t582;
t566 = -qJD(4) * t526 + t472 * t545 + t473 * t587 + t588;
t419 = -pkin(4) * t589 - pkin(9) * t510 + t511 * t604;
t420 = pkin(4) * t590 - pkin(9) * t512 + t513 * t604;
t565 = t419 * t545 + t420 * t587 + t566;
t563 = t546 * t420 + (t558 * t574 - t592) * t553 + t579;
t562 = (-t419 + t597) * t546 + t574 * t599 + t572;
t559 = cos(qJ(6));
t556 = sin(qJ(6));
t551 = sin(t591);
t541 = rSges(2,1) * t560 - rSges(2,2) * t558;
t540 = rSges(2,1) * t558 + rSges(2,2) * t560;
t525 = t603 * rSges(3,3) + (rSges(3,1) * t557 + rSges(3,2) * t606) * t553;
t524 = Icges(3,5) * t603 + (Icges(3,1) * t557 + Icges(3,4) * t606) * t553;
t523 = Icges(3,6) * t603 + (Icges(3,4) * t557 + Icges(3,2) * t606) * t553;
t515 = -t527 * t582 + t551 * t603;
t514 = -t527 * t551 - t582 * t603;
t502 = rSges(3,1) * t533 + rSges(3,2) * t532 + rSges(3,3) * t600;
t501 = rSges(3,1) * t531 + rSges(3,2) * t530 - rSges(3,3) * t599;
t500 = Icges(3,1) * t533 + Icges(3,4) * t532 + Icges(3,5) * t600;
t499 = Icges(3,1) * t531 + Icges(3,4) * t530 - Icges(3,5) * t599;
t498 = Icges(3,4) * t533 + Icges(3,2) * t532 + Icges(3,6) * t600;
t497 = Icges(3,4) * t531 + Icges(3,2) * t530 - Icges(3,6) * t599;
t493 = -Icges(4,1) * t527 + Icges(4,4) * t526 + Icges(4,5) * t603;
t481 = t513 * t582 + t551 * t600;
t480 = t513 * t551 - t558 * t569;
t479 = t511 * t582 - t551 * t599;
t478 = t511 * t551 + t560 * t569;
t477 = t515 * t559 - t526 * t556;
t476 = -t515 * t556 - t526 * t559;
t475 = qJD(6) * t514 + t518;
t474 = pkin(5) * t515 + pkin(10) * t514;
t470 = Icges(5,1) * t517 + Icges(5,4) * t516 - Icges(5,5) * t526;
t469 = Icges(5,4) * t517 + Icges(5,2) * t516 - Icges(5,6) * t526;
t467 = rSges(4,1) * t513 + rSges(4,2) * t512 + rSges(4,3) * t600;
t466 = rSges(4,1) * t511 + rSges(4,2) * t510 - rSges(4,3) * t599;
t462 = Icges(4,1) * t513 + Icges(4,4) * t512 + Icges(4,5) * t600;
t461 = Icges(4,1) * t511 + Icges(4,4) * t510 - Icges(4,5) * t599;
t456 = t502 * t546 - t525 * t545 + t534;
t455 = -t501 * t546 - t525 * t587 - t594;
t454 = rSges(6,1) * t515 - rSges(6,2) * t514 - rSges(6,3) * t526;
t453 = Icges(6,1) * t515 - Icges(6,4) * t514 - Icges(6,5) * t526;
t452 = Icges(6,4) * t515 - Icges(6,2) * t514 - Icges(6,6) * t526;
t451 = Icges(6,5) * t515 - Icges(6,6) * t514 - Icges(6,3) * t526;
t449 = (t501 * t558 + t502 * t560) * t593;
t448 = t481 * t559 - t512 * t556;
t447 = -t481 * t556 - t512 * t559;
t446 = t479 * t559 - t510 * t556;
t445 = -t479 * t556 - t510 * t559;
t444 = qJD(6) * t478 + t487;
t443 = qJD(6) * t480 + t486;
t442 = pkin(5) * t481 + pkin(10) * t480;
t441 = pkin(5) * t479 + pkin(10) * t478;
t440 = rSges(5,1) * t485 + rSges(5,2) * t484 - rSges(5,3) * t512;
t439 = rSges(5,1) * t483 + rSges(5,2) * t482 - rSges(5,3) * t510;
t438 = Icges(5,1) * t485 + Icges(5,4) * t484 - Icges(5,5) * t512;
t437 = Icges(5,1) * t483 + Icges(5,4) * t482 - Icges(5,5) * t510;
t436 = Icges(5,4) * t485 + Icges(5,2) * t484 - Icges(5,6) * t512;
t435 = Icges(5,4) * t483 + Icges(5,2) * t482 - Icges(5,6) * t510;
t432 = rSges(6,1) * t481 - rSges(6,2) * t480 - rSges(6,3) * t512;
t431 = rSges(6,1) * t479 - rSges(6,2) * t478 - rSges(6,3) * t510;
t430 = Icges(6,1) * t481 - Icges(6,4) * t480 - Icges(6,5) * t512;
t429 = Icges(6,1) * t479 - Icges(6,4) * t478 - Icges(6,5) * t510;
t428 = Icges(6,4) * t481 - Icges(6,2) * t480 - Icges(6,6) * t512;
t427 = Icges(6,4) * t479 - Icges(6,2) * t478 - Icges(6,6) * t510;
t426 = Icges(6,5) * t481 - Icges(6,6) * t480 - Icges(6,3) * t512;
t425 = Icges(6,5) * t479 - Icges(6,6) * t478 - Icges(6,3) * t510;
t424 = rSges(7,1) * t477 + rSges(7,2) * t476 + rSges(7,3) * t514;
t423 = Icges(7,1) * t477 + Icges(7,4) * t476 + Icges(7,5) * t514;
t422 = Icges(7,4) * t477 + Icges(7,2) * t476 + Icges(7,6) * t514;
t421 = Icges(7,5) * t477 + Icges(7,6) * t476 + Icges(7,3) * t514;
t415 = t467 * t546 + (t558 * t581 - t592) * t553 + t596;
t414 = (-t466 - t508) * t546 + t581 * t599 + t580;
t413 = rSges(7,1) * t448 + rSges(7,2) * t447 + rSges(7,3) * t480;
t412 = rSges(7,1) * t446 + rSges(7,2) * t445 + rSges(7,3) * t478;
t411 = Icges(7,1) * t448 + Icges(7,4) * t447 + Icges(7,5) * t480;
t410 = Icges(7,1) * t446 + Icges(7,4) * t445 + Icges(7,5) * t478;
t409 = Icges(7,4) * t448 + Icges(7,2) * t447 + Icges(7,6) * t480;
t408 = Icges(7,4) * t446 + Icges(7,2) * t445 + Icges(7,6) * t478;
t407 = Icges(7,5) * t448 + Icges(7,6) * t447 + Icges(7,3) * t480;
t406 = Icges(7,5) * t446 + Icges(7,6) * t445 + Icges(7,3) * t478;
t405 = (t466 * t558 + t467 * t560) * t593 + t588;
t404 = t440 * t546 + (t558 * t573 - t592) * t553 + t579;
t403 = (-t439 + t597) * t546 + t573 * t599 + t572;
t402 = (t439 * t558 + t440 * t560) * t593 + t566;
t401 = t432 * t518 - t454 * t486 + t563;
t400 = -t431 * t518 + t454 * t487 + t562;
t399 = t431 * t486 - t432 * t487 + t565;
t398 = t413 * t475 - t424 * t443 + t442 * t518 - t474 * t486 + t563;
t397 = -t412 * t475 + t424 * t444 - t441 * t518 + t474 * t487 + t562;
t396 = t412 * t443 - t413 * t444 + t441 * t486 - t442 * t487 + t565;
t1 = t444 * ((t407 * t478 + t409 * t445 + t411 * t446) * t443 + (t478 * t406 + t445 * t408 + t446 * t410) * t444 + (t421 * t478 + t422 * t445 + t423 * t446) * t475) / 0.2e1 + t475 * ((t407 * t514 + t409 * t476 + t411 * t477) * t443 + (t406 * t514 + t408 * t476 + t410 * t477) * t444 + (t421 * t514 + t422 * t476 + t423 * t477) * t475) / 0.2e1 + t443 * ((t480 * t407 + t447 * t409 + t448 * t411) * t443 + (t406 * t480 + t408 * t447 + t410 * t448) * t444 + (t421 * t480 + t422 * t447 + t423 * t448) * t475) / 0.2e1 + t486 * ((-t426 * t512 - t428 * t480 + t430 * t481) * t486 + (-t425 * t512 - t427 * t480 + t429 * t481) * t487 + (-t451 * t512 - t452 * t480 + t453 * t481) * t518) / 0.2e1 + t487 * ((-t426 * t510 - t428 * t478 + t430 * t479) * t486 + (-t425 * t510 - t427 * t478 + t429 * t479) * t487 + (-t451 * t510 - t452 * t478 + t453 * t479) * t518) / 0.2e1 + t518 * ((-t426 * t526 - t428 * t514 + t430 * t515) * t486 + (-t425 * t526 - t427 * t514 + t429 * t515) * t487 + (-t451 * t526 - t452 * t514 + t453 * t515) * t518) / 0.2e1 + m(7) * (t396 ^ 2 + t397 ^ 2 + t398 ^ 2) / 0.2e1 + m(6) * (t399 ^ 2 + t400 ^ 2 + t401 ^ 2) / 0.2e1 + m(5) * (t402 ^ 2 + t403 ^ 2 + t404 ^ 2) / 0.2e1 + m(4) * (t405 ^ 2 + t414 ^ 2 + t415 ^ 2) / 0.2e1 + m(3) * (t449 ^ 2 + t455 ^ 2 + t456 ^ 2) / 0.2e1 + (Icges(2,3) + m(2) * (t540 ^ 2 + t541 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t603 * t496 + (t498 * t606 + t500 * t557) * t553) * t545 - (t603 * t495 + (t497 * t606 + t499 * t557) * t553) * t587 + ((-t435 * t516 - t437 * t517 - t457 * t603 + t527 * t461 + t526 * t613) * t560 + (t436 * t516 + t438 * t517 + t458 * t603 - t527 * t462 - t526 * t612) * t558) * t593 + ((t523 * t606 + t524 * t557) * t553 - t527 * t493 + t469 * t516 + t470 * t517 + t610 * t603 - t611 * t526) * t546) * t546 / 0.2e1 + (((-t435 * t484 - t437 * t485 - t461 * t513 - t497 * t532 - t499 * t533 + t512 * t613) * t560 + (t436 * t484 + t438 * t485 + t513 * t462 + t532 * t498 + t533 * t500 - t512 * t612 - t614) * t558) * t593 + (t469 * t484 + t470 * t485 + t493 * t513 - t512 * t611 + t523 * t532 + t524 * t533 + t600 * t610) * t546) * t545 / 0.2e1 - (((-t435 * t482 - t437 * t483 - t511 * t461 - t530 * t497 - t531 * t499 + t510 * t613 + t614) * t560 + (t436 * t482 + t438 * t483 + t462 * t511 + t498 * t530 + t500 * t531 - t510 * t612) * t558) * t593 + (t469 * t482 + t470 * t483 + t493 * t511 - t510 * t611 + t523 * t530 + t524 * t531 - t599 * t610) * t546) * t587 / 0.2e1;
T  = t1;
