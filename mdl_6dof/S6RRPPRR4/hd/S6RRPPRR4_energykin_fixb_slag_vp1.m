% Calculate kinetic energy for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:01:09
% EndTime: 2019-03-09 09:01:13
% DurationCPUTime: 3.66s
% Computational Cost: add. (2707->335), mult. (6786->501), div. (0->0), fcn. (8432->12), ass. (0->157)
t592 = Icges(5,1) + Icges(4,3);
t589 = -Icges(5,4) + Icges(4,5);
t588 = Icges(5,5) - Icges(4,6);
t591 = Icges(4,1) + Icges(5,2);
t590 = Icges(4,4) + Icges(5,6);
t587 = -Icges(4,2) - Icges(5,3);
t528 = sin(qJ(1));
t531 = cos(qJ(1));
t527 = sin(qJ(2));
t530 = cos(qJ(2));
t569 = sin(pkin(11));
t570 = cos(pkin(11));
t511 = -t527 * t569 + t530 * t570;
t524 = cos(pkin(6));
t536 = t524 * t511;
t538 = t527 * t570 + t530 * t569;
t487 = -t528 * t536 - t531 * t538;
t535 = t524 * t538;
t488 = t511 * t531 - t528 * t535;
t523 = sin(pkin(6));
t568 = t523 * t528;
t586 = -t588 * t487 + t589 * t488 + t592 * t568;
t485 = -t528 * t538 + t531 * t536;
t486 = t528 * t511 + t531 * t535;
t567 = t523 * t531;
t585 = -t588 * t485 + t589 * t486 - t592 * t567;
t563 = t530 * t531;
t566 = t527 * t528;
t505 = t524 * t563 - t566;
t564 = t528 * t530;
t565 = t527 * t531;
t506 = t524 * t565 + t564;
t507 = -t524 * t564 - t565;
t508 = -t524 * t566 + t563;
t540 = (Icges(3,5) * t506 + Icges(3,6) * t505 - Icges(3,3) * t567) * t531 - (Icges(3,5) * t508 + Icges(3,6) * t507 + Icges(3,3) * t568) * t528;
t584 = (-t586 * t528 + t585 * t531 + t540) * t523;
t583 = t487 * t587 - t488 * t590 + t568 * t588;
t582 = t485 * t587 - t486 * t590 - t567 * t588;
t581 = t590 * t487 + t488 * t591 + t589 * t568;
t580 = -t590 * t485 - t486 * t591 + t589 * t567;
t502 = t511 * t523;
t503 = t538 * t523;
t579 = t587 * t502 - t503 * t590 + t524 * t588;
t578 = t590 * t502 + t503 * t591 + t589 * t524;
t577 = (Icges(3,5) * t527 + Icges(3,6) * t530) * t523 + t589 * t503 - t588 * t502 + (Icges(3,3) + t592) * t524;
t573 = cos(qJ(5));
t572 = pkin(2) * t527;
t571 = pkin(2) * t530;
t441 = pkin(3) * t486 - qJ(4) * t485;
t553 = -qJ(3) * t523 + t524 * t572;
t483 = t528 * t571 + t531 * t553;
t562 = -t441 - t483;
t484 = -t528 * t553 + t531 * t571;
t509 = qJD(1) * (pkin(1) * t531 + pkin(8) * t568);
t519 = qJD(2) * t524 + qJD(1);
t561 = t519 * t484 + t509;
t512 = qJ(3) * t524 + t523 * t572;
t560 = -pkin(3) * t503 + qJ(4) * t502 - t512;
t558 = qJD(2) * t523;
t518 = t528 * t558;
t454 = qJD(5) * t488 + t518;
t559 = qJD(1) * (pkin(1) * t528 - pkin(8) * t567);
t557 = qJD(3) * t531;
t554 = t531 * t558;
t556 = qJD(3) * t524 + t483 * t518 + t484 * t554;
t555 = t523 * t573;
t491 = qJD(5) * t503 + t519;
t550 = qJD(2) * (-rSges(4,1) * t503 - rSges(4,2) * t502 - rSges(4,3) * t524 - t512);
t549 = qJD(3) * t568 - t559;
t442 = pkin(3) * t488 - qJ(4) * t487;
t548 = -qJD(4) * t485 + t519 * t442 + t561;
t545 = qJD(2) * (-rSges(5,1) * t524 + rSges(5,2) * t503 + rSges(5,3) * t502 + t560);
t544 = qJD(2) * (-pkin(4) * t524 - pkin(9) * t503 + t560);
t455 = qJD(5) * t486 - t554;
t543 = -qJD(4) * t487 + t549;
t539 = -qJD(4) * t502 + t441 * t518 + t442 * t554 + t556;
t459 = pkin(4) * t568 + pkin(9) * t488;
t460 = -pkin(4) * t567 + pkin(9) * t486;
t537 = t459 * t554 + t460 * t518 + t539;
t534 = t519 * t459 + (t528 * t544 - t557) * t523 + t548;
t533 = (-t460 + t562) * t519 + t544 * t567 + t543;
t529 = cos(qJ(6));
t526 = sin(qJ(5));
t525 = sin(qJ(6));
t515 = rSges(2,1) * t531 - rSges(2,2) * t528;
t514 = rSges(2,1) * t528 + rSges(2,2) * t531;
t501 = rSges(3,3) * t524 + (rSges(3,1) * t527 + rSges(3,2) * t530) * t523;
t500 = Icges(3,5) * t524 + (Icges(3,1) * t527 + Icges(3,4) * t530) * t523;
t499 = Icges(3,6) * t524 + (Icges(3,4) * t527 + Icges(3,2) * t530) * t523;
t490 = -t502 * t526 + t524 * t573;
t489 = t502 * t573 + t524 * t526;
t477 = rSges(3,1) * t508 + rSges(3,2) * t507 + rSges(3,3) * t568;
t476 = rSges(3,1) * t506 + rSges(3,2) * t505 - rSges(3,3) * t567;
t474 = Icges(3,1) * t508 + Icges(3,4) * t507 + Icges(3,5) * t568;
t473 = Icges(3,1) * t506 + Icges(3,4) * t505 - Icges(3,5) * t567;
t472 = Icges(3,4) * t508 + Icges(3,2) * t507 + Icges(3,6) * t568;
t471 = Icges(3,4) * t506 + Icges(3,2) * t505 - Icges(3,6) * t567;
t453 = -t485 * t526 - t531 * t555;
t452 = -t485 * t573 + t526 * t567;
t451 = -t487 * t526 + t528 * t555;
t450 = t487 * t573 + t526 * t568;
t449 = t490 * t529 + t503 * t525;
t448 = -t490 * t525 + t503 * t529;
t444 = qJD(6) * t489 + t491;
t443 = pkin(5) * t490 + pkin(10) * t489;
t440 = rSges(6,1) * t490 - rSges(6,2) * t489 + rSges(6,3) * t503;
t439 = Icges(6,1) * t490 - Icges(6,4) * t489 + Icges(6,5) * t503;
t438 = Icges(6,4) * t490 - Icges(6,2) * t489 + Icges(6,6) * t503;
t437 = Icges(6,5) * t490 - Icges(6,6) * t489 + Icges(6,3) * t503;
t436 = rSges(4,1) * t488 + rSges(4,2) * t487 + rSges(4,3) * t568;
t435 = rSges(4,1) * t486 + rSges(4,2) * t485 - rSges(4,3) * t567;
t434 = -rSges(5,1) * t567 - rSges(5,2) * t486 - rSges(5,3) * t485;
t433 = rSges(5,1) * t568 - rSges(5,2) * t488 - rSges(5,3) * t487;
t417 = t477 * t519 - t501 * t518 + t509;
t416 = -t476 * t519 - t501 * t554 - t559;
t415 = t453 * t529 + t486 * t525;
t414 = -t453 * t525 + t486 * t529;
t413 = t451 * t529 + t488 * t525;
t412 = -t451 * t525 + t488 * t529;
t411 = (t476 * t528 + t477 * t531) * t558;
t410 = -qJD(6) * t452 + t455;
t409 = qJD(6) * t450 + t454;
t408 = pkin(5) * t453 - pkin(10) * t452;
t407 = pkin(5) * t451 + pkin(10) * t450;
t406 = rSges(7,1) * t449 + rSges(7,2) * t448 + rSges(7,3) * t489;
t405 = Icges(7,1) * t449 + Icges(7,4) * t448 + Icges(7,5) * t489;
t404 = Icges(7,4) * t449 + Icges(7,2) * t448 + Icges(7,6) * t489;
t403 = Icges(7,5) * t449 + Icges(7,6) * t448 + Icges(7,3) * t489;
t402 = rSges(6,1) * t453 + rSges(6,2) * t452 + rSges(6,3) * t486;
t401 = rSges(6,1) * t451 - rSges(6,2) * t450 + rSges(6,3) * t488;
t400 = Icges(6,1) * t453 + Icges(6,4) * t452 + Icges(6,5) * t486;
t399 = Icges(6,1) * t451 - Icges(6,4) * t450 + Icges(6,5) * t488;
t398 = Icges(6,4) * t453 + Icges(6,2) * t452 + Icges(6,6) * t486;
t397 = Icges(6,4) * t451 - Icges(6,2) * t450 + Icges(6,6) * t488;
t396 = Icges(6,5) * t453 + Icges(6,6) * t452 + Icges(6,3) * t486;
t395 = Icges(6,5) * t451 - Icges(6,6) * t450 + Icges(6,3) * t488;
t394 = t436 * t519 + (t528 * t550 - t557) * t523 + t561;
t393 = (-t435 - t483) * t519 + t550 * t567 + t549;
t392 = rSges(7,1) * t415 + rSges(7,2) * t414 - rSges(7,3) * t452;
t391 = rSges(7,1) * t413 + rSges(7,2) * t412 + rSges(7,3) * t450;
t390 = Icges(7,1) * t415 + Icges(7,4) * t414 - Icges(7,5) * t452;
t389 = Icges(7,1) * t413 + Icges(7,4) * t412 + Icges(7,5) * t450;
t388 = Icges(7,4) * t415 + Icges(7,2) * t414 - Icges(7,6) * t452;
t387 = Icges(7,4) * t413 + Icges(7,2) * t412 + Icges(7,6) * t450;
t386 = Icges(7,5) * t415 + Icges(7,6) * t414 - Icges(7,3) * t452;
t385 = Icges(7,5) * t413 + Icges(7,6) * t412 + Icges(7,3) * t450;
t384 = (t435 * t528 + t436 * t531) * t558 + t556;
t383 = t433 * t519 + (t528 * t545 - t557) * t523 + t548;
t382 = (-t434 + t562) * t519 + t545 * t567 + t543;
t381 = (t433 * t531 + t434 * t528) * t558 + t539;
t380 = t401 * t491 - t440 * t454 + t534;
t379 = -t402 * t491 + t440 * t455 + t533;
t378 = -t401 * t455 + t402 * t454 + t537;
t377 = t391 * t444 - t406 * t409 + t407 * t491 - t443 * t454 + t534;
t376 = -t392 * t444 + t406 * t410 - t408 * t491 + t443 * t455 + t533;
t375 = -t391 * t410 + t392 * t409 - t407 * t455 + t408 * t454 + t537;
t1 = m(7) * (t375 ^ 2 + t376 ^ 2 + t377 ^ 2) / 0.2e1 + m(6) * (t378 ^ 2 + t379 ^ 2 + t380 ^ 2) / 0.2e1 + m(5) * (t381 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + m(4) * (t384 ^ 2 + t393 ^ 2 + t394 ^ 2) / 0.2e1 + m(3) * (t411 ^ 2 + t416 ^ 2 + t417 ^ 2) / 0.2e1 + t409 * ((t450 * t385 + t412 * t387 + t413 * t389) * t409 + (t386 * t450 + t388 * t412 + t390 * t413) * t410 + (t403 * t450 + t404 * t412 + t405 * t413) * t444) / 0.2e1 + t410 * ((-t385 * t452 + t387 * t414 + t389 * t415) * t409 + (-t452 * t386 + t414 * t388 + t415 * t390) * t410 + (-t403 * t452 + t404 * t414 + t405 * t415) * t444) / 0.2e1 + t444 * ((t385 * t489 + t387 * t448 + t389 * t449) * t409 + (t386 * t489 + t388 * t448 + t390 * t449) * t410 + (t489 * t403 + t448 * t404 + t449 * t405) * t444) / 0.2e1 + t455 * ((t395 * t486 + t397 * t452 + t399 * t453) * t454 + (t486 * t396 + t452 * t398 + t453 * t400) * t455 + (t437 * t486 + t438 * t452 + t439 * t453) * t491) / 0.2e1 + t454 * ((t488 * t395 - t450 * t397 + t451 * t399) * t454 + (t396 * t488 - t398 * t450 + t400 * t451) * t455 + (t437 * t488 - t438 * t450 + t439 * t451) * t491) / 0.2e1 + t491 * ((t395 * t503 - t397 * t489 + t399 * t490) * t454 + (t396 * t503 - t398 * t489 + t400 * t490) * t455 + (t503 * t437 - t489 * t438 + t490 * t439) * t491) / 0.2e1 + (Icges(2,3) + m(2) * (t514 ^ 2 + t515 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((-t540 * t524 + (-(t471 * t530 + t473 * t527) * t523 - t585 * t524 + t580 * t503 + t582 * t502) * t531 + ((t472 * t530 + t474 * t527) * t523 + t586 * t524 + t581 * t503 - t583 * t502) * t528) * t558 + ((t499 * t530 + t500 * t527) * t523 + t578 * t503 - t579 * t502 + t577 * t524) * t519) * t519 / 0.2e1 + (((-t471 * t507 - t473 * t508 + t487 * t582 + t488 * t580) * t531 + (t507 * t472 + t508 * t474 - t487 * t583 + t581 * t488 - t584) * t528) * t558 + (-t487 * t579 + t488 * t578 + t499 * t507 + t500 * t508 + t568 * t577) * t519) * t518 / 0.2e1 - (((-t505 * t471 - t506 * t473 + t485 * t582 + t486 * t580 + t584) * t531 + (t472 * t505 + t474 * t506 - t485 * t583 + t581 * t486) * t528) * t558 + (-t485 * t579 + t486 * t578 + t499 * t505 + t500 * t506 - t567 * t577) * t519) * t554 / 0.2e1;
T  = t1;
