% Calculate kinetic energy for
% S6PRPRPR3
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
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:34:27
% EndTime: 2019-03-08 19:34:29
% DurationCPUTime: 3.10s
% Computational Cost: add. (2962->328), mult. (7553->497), div. (0->0), fcn. (9514->12), ass. (0->153)
t588 = Icges(5,1) + Icges(6,2);
t587 = -Icges(6,1) - Icges(5,3);
t586 = -Icges(5,4) - Icges(6,6);
t585 = Icges(6,4) - Icges(5,5);
t584 = Icges(6,5) - Icges(5,6);
t583 = Icges(5,2) + Icges(6,3);
t527 = sin(pkin(10));
t529 = cos(pkin(10));
t532 = sin(qJ(2));
t560 = sin(pkin(11));
t561 = cos(pkin(11));
t565 = cos(qJ(2));
t516 = -t532 * t560 + t565 * t561;
t530 = cos(pkin(6));
t538 = t530 * t516;
t539 = t532 * t561 + t560 * t565;
t493 = -t527 * t539 + t529 * t538;
t509 = t539 * t530;
t494 = t509 * t529 + t516 * t527;
t528 = sin(pkin(6));
t558 = t529 * t528;
t437 = Icges(4,5) * t494 + Icges(4,6) * t493 - Icges(4,3) * t558;
t495 = -t527 * t538 - t529 * t539;
t496 = -t509 * t527 + t516 * t529;
t559 = t527 * t528;
t438 = Icges(4,5) * t496 + Icges(4,6) * t495 + Icges(4,3) * t559;
t507 = t516 * t528;
t508 = t539 * t528;
t476 = Icges(4,5) * t508 + Icges(4,6) * t507 + Icges(4,3) * t530;
t550 = t530 * t565;
t511 = -t527 * t532 + t529 * t550;
t557 = t530 * t532;
t512 = t527 * t565 + t529 * t557;
t480 = Icges(3,5) * t512 + Icges(3,6) * t511 - Icges(3,3) * t558;
t513 = -t527 * t550 - t529 * t532;
t514 = -t527 * t557 + t529 * t565;
t481 = Icges(3,5) * t514 + Icges(3,6) * t513 + Icges(3,3) * t559;
t503 = Icges(3,3) * t530 + (Icges(3,5) * t532 + Icges(3,6) * t565) * t528;
t582 = -(t476 + t503) * t530 + (t480 + t437) * t558 - (t481 + t438) * t559;
t564 = cos(qJ(4));
t552 = t528 * t564;
t563 = sin(qJ(4));
t466 = t494 * t563 + t529 * t552;
t551 = t528 * t563;
t467 = t494 * t564 - t529 * t551;
t581 = t466 * t583 + t467 * t586 - t493 * t584;
t468 = t496 * t563 - t527 * t552;
t469 = t496 * t564 + t527 * t551;
t580 = t468 * t583 + t469 * t586 - t495 * t584;
t579 = t466 * t584 - t467 * t585 + t493 * t587;
t578 = t468 * t584 - t469 * t585 + t495 * t587;
t577 = t586 * t466 + t467 * t588 + t585 * t493;
t576 = t586 * t468 + t469 * t588 + t585 * t495;
t498 = t508 * t563 - t530 * t564;
t499 = t508 * t564 + t530 * t563;
t575 = t498 * t583 + t499 * t586 - t507 * t584;
t574 = t498 * t584 - t499 * t585 + t507 * t587;
t573 = t586 * t498 + t499 * t588 + t585 * t507;
t569 = qJD(2) ^ 2;
t562 = pkin(2) * t565;
t517 = pkin(2) * t528 * t532 + qJ(3) * t530;
t556 = -pkin(3) * t508 + pkin(8) * t507 - t517;
t555 = qJD(2) * t528;
t522 = t527 * t555;
t471 = -qJD(4) * t495 + t522;
t526 = qJD(2) * t530;
t500 = -qJD(4) * t507 + t526;
t554 = qJD(3) * t529;
t549 = t529 * t555;
t548 = pkin(2) * t557 - qJ(3) * t528;
t546 = (-rSges(4,1) * t508 - rSges(4,2) * t507 - rSges(4,3) * t530 - t517) * t528;
t488 = t527 * t562 + t529 * t548;
t489 = -t527 * t548 + t529 * t562;
t545 = qJD(3) * t530 + t488 * t522 + t489 * t549 + qJD(1);
t472 = -qJD(4) * t493 - t549;
t448 = pkin(3) * t494 - pkin(8) * t493;
t449 = pkin(3) * t496 - pkin(8) * t495;
t541 = t448 * t522 + t449 * t549 + t545;
t426 = pkin(4) * t467 + qJ(5) * t466;
t540 = qJD(5) * t498 + t471 * t426 + t541;
t475 = t489 * t526;
t537 = t449 * t526 + t475 + (qJD(2) * t527 * t556 - t554) * t528;
t521 = qJD(3) * t559;
t536 = t521 + ((-t448 - t488) * t530 + t556 * t558) * qJD(2);
t427 = pkin(4) * t469 + qJ(5) * t468;
t535 = qJD(5) * t466 + t500 * t427 + t537;
t460 = pkin(4) * t499 + qJ(5) * t498;
t534 = qJD(5) * t468 + t472 * t460 + t536;
t533 = cos(qJ(6));
t531 = sin(qJ(6));
t506 = t530 * rSges(3,3) + (rSges(3,1) * t532 + rSges(3,2) * t565) * t528;
t505 = Icges(3,5) * t530 + (Icges(3,1) * t532 + Icges(3,4) * t565) * t528;
t504 = Icges(3,6) * t530 + (Icges(3,4) * t532 + Icges(3,2) * t565) * t528;
t487 = rSges(3,1) * t514 + rSges(3,2) * t513 + rSges(3,3) * t559;
t486 = rSges(3,1) * t512 + rSges(3,2) * t511 - rSges(3,3) * t558;
t485 = Icges(3,1) * t514 + Icges(3,4) * t513 + Icges(3,5) * t559;
t484 = Icges(3,1) * t512 + Icges(3,4) * t511 - Icges(3,5) * t558;
t483 = Icges(3,4) * t514 + Icges(3,2) * t513 + Icges(3,6) * t559;
t482 = Icges(3,4) * t512 + Icges(3,2) * t511 - Icges(3,6) * t558;
t478 = Icges(4,1) * t508 + Icges(4,4) * t507 + Icges(4,5) * t530;
t477 = Icges(4,4) * t508 + Icges(4,2) * t507 + Icges(4,6) * t530;
t470 = -pkin(5) * t507 + pkin(9) * t499;
t463 = t498 * t531 - t507 * t533;
t462 = t498 * t533 + t507 * t531;
t461 = qJD(6) * t499 + t500;
t459 = (-t486 * t530 - t506 * t558) * qJD(2);
t458 = (t487 * t530 - t506 * t559) * qJD(2);
t457 = rSges(5,1) * t499 - rSges(5,2) * t498 - rSges(5,3) * t507;
t456 = -rSges(6,1) * t507 - rSges(6,2) * t499 + rSges(6,3) * t498;
t444 = rSges(4,1) * t496 + rSges(4,2) * t495 + rSges(4,3) * t559;
t443 = rSges(4,1) * t494 + rSges(4,2) * t493 - rSges(4,3) * t558;
t442 = Icges(4,1) * t496 + Icges(4,4) * t495 + Icges(4,5) * t559;
t441 = Icges(4,1) * t494 + Icges(4,4) * t493 - Icges(4,5) * t558;
t440 = Icges(4,4) * t496 + Icges(4,2) * t495 + Icges(4,6) * t559;
t439 = Icges(4,4) * t494 + Icges(4,2) * t493 - Icges(4,6) * t558;
t436 = -pkin(5) * t495 + pkin(9) * t469;
t435 = -pkin(5) * t493 + pkin(9) * t467;
t434 = t468 * t531 - t495 * t533;
t433 = t468 * t533 + t495 * t531;
t432 = t466 * t531 - t493 * t533;
t431 = t466 * t533 + t493 * t531;
t430 = qJD(1) + (t486 * t527 + t487 * t529) * t555;
t429 = qJD(6) * t467 + t472;
t428 = qJD(6) * t469 + t471;
t423 = rSges(7,1) * t463 + rSges(7,2) * t462 + rSges(7,3) * t499;
t422 = Icges(7,1) * t463 + Icges(7,4) * t462 + Icges(7,5) * t499;
t421 = Icges(7,4) * t463 + Icges(7,2) * t462 + Icges(7,6) * t499;
t420 = Icges(7,5) * t463 + Icges(7,6) * t462 + Icges(7,3) * t499;
t419 = rSges(5,1) * t469 - rSges(5,2) * t468 - rSges(5,3) * t495;
t418 = rSges(5,1) * t467 - rSges(5,2) * t466 - rSges(5,3) * t493;
t417 = -rSges(6,1) * t495 - rSges(6,2) * t469 + rSges(6,3) * t468;
t416 = -rSges(6,1) * t493 - rSges(6,2) * t467 + rSges(6,3) * t466;
t402 = t521 + ((-t443 - t488) * t530 + t529 * t546) * qJD(2);
t401 = -t528 * t554 + t475 + (t444 * t530 + t527 * t546) * qJD(2);
t400 = rSges(7,1) * t434 + rSges(7,2) * t433 + rSges(7,3) * t469;
t399 = rSges(7,1) * t432 + rSges(7,2) * t431 + rSges(7,3) * t467;
t398 = Icges(7,1) * t434 + Icges(7,4) * t433 + Icges(7,5) * t469;
t397 = Icges(7,1) * t432 + Icges(7,4) * t431 + Icges(7,5) * t467;
t396 = Icges(7,4) * t434 + Icges(7,2) * t433 + Icges(7,6) * t469;
t395 = Icges(7,4) * t432 + Icges(7,2) * t431 + Icges(7,6) * t467;
t394 = Icges(7,5) * t434 + Icges(7,6) * t433 + Icges(7,3) * t469;
t393 = Icges(7,5) * t432 + Icges(7,6) * t431 + Icges(7,3) * t467;
t392 = (t443 * t527 + t444 * t529) * t555 + t545;
t391 = -t418 * t500 + t457 * t472 + t536;
t390 = t419 * t500 - t457 * t471 + t537;
t389 = t418 * t471 - t419 * t472 + t541;
t388 = t456 * t472 + (-t416 - t426) * t500 + t534;
t387 = t417 * t500 + (-t456 - t460) * t471 + t535;
t386 = t416 * t471 + (-t417 - t427) * t472 + t540;
t385 = -t399 * t461 + t423 * t429 + t470 * t472 + (-t426 - t435) * t500 + t534;
t384 = t400 * t461 - t423 * t428 + t436 * t500 + (-t460 - t470) * t471 + t535;
t383 = t399 * t428 - t400 * t429 + t435 * t471 + (-t427 - t436) * t472 + t540;
t1 = m(5) * (t389 ^ 2 + t390 ^ 2 + t391 ^ 2) / 0.2e1 + m(6) * (t386 ^ 2 + t387 ^ 2 + t388 ^ 2) / 0.2e1 + m(7) * (t383 ^ 2 + t384 ^ 2 + t385 ^ 2) / 0.2e1 + ((t530 * t481 + (t483 * t565 + t485 * t532) * t528) * t522 - (t530 * t480 + (t482 * t565 + t484 * t532) * t528) * t549 + (t530 * t503 + (t504 * t565 + t505 * t532) * t528) * t526) * t526 / 0.2e1 + t428 * ((t469 * t394 + t433 * t396 + t434 * t398) * t428 + (t393 * t469 + t395 * t433 + t397 * t434) * t429 + (t420 * t469 + t421 * t433 + t422 * t434) * t461) / 0.2e1 + t429 * ((t394 * t467 + t396 * t431 + t398 * t432) * t428 + (t467 * t393 + t431 * t395 + t432 * t397) * t429 + (t420 * t467 + t421 * t431 + t422 * t432) * t461) / 0.2e1 + t461 * ((t394 * t499 + t396 * t462 + t398 * t463) * t428 + (t393 * t499 + t395 * t462 + t397 * t463) * t429 + (t499 * t420 + t462 * t421 + t463 * t422) * t461) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t430 ^ 2 + t458 ^ 2 + t459 ^ 2) / 0.2e1 + m(4) * (t392 ^ 2 + t401 ^ 2 + t402 ^ 2) / 0.2e1 + ((t468 * t575 + t469 * t573 - t495 * t574) * t500 + (t468 * t581 + t577 * t469 - t579 * t495) * t472 + (t580 * t468 + t576 * t469 - t578 * t495) * t471) * t471 / 0.2e1 + ((t466 * t575 + t467 * t573 - t493 * t574) * t500 + (t581 * t466 + t577 * t467 - t579 * t493) * t472 + (t466 * t580 + t467 * t576 - t493 * t578) * t471) * t472 / 0.2e1 + ((t575 * t498 + t573 * t499 - t574 * t507) * t500 + (t498 * t581 + t577 * t499 - t579 * t507) * t472 + (t498 * t580 + t499 * t576 - t507 * t578) * t471) * t500 / 0.2e1 - ((t440 * t493 + t442 * t494 + t483 * t511 + t485 * t512) * t559 + (t477 * t493 + t478 * t494 + t504 * t511 + t505 * t512) * t530 + (-t439 * t493 - t441 * t494 - t482 * t511 - t484 * t512 + t582) * t558) * t569 * t558 / 0.2e1 + (t530 * ((t476 * t530 + t477 * t507 + t478 * t508) * t530 + ((t438 * t530 + t440 * t507 + t442 * t508) * t527 - (t437 * t530 + t439 * t507 + t441 * t508) * t529) * t528) + ((-t439 * t495 - t441 * t496 - t482 * t513 - t484 * t514) * t558 + (t477 * t495 + t478 * t496 + t504 * t513 + t505 * t514) * t530 + (t440 * t495 + t442 * t496 + t483 * t513 + t485 * t514 - t582) * t559) * t559) * t569 / 0.2e1;
T  = t1;
