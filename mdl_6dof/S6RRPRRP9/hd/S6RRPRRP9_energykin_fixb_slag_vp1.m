% Calculate kinetic energy for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP9_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP9_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP9_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP9_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:27:29
% EndTime: 2019-03-09 12:27:32
% DurationCPUTime: 2.96s
% Computational Cost: add. (3183->336), mult. (5709->502), div. (0->0), fcn. (6885->12), ass. (0->161)
t599 = Icges(6,1) + Icges(7,1);
t598 = Icges(6,4) + Icges(7,4);
t597 = Icges(6,5) + Icges(7,5);
t596 = Icges(6,2) + Icges(7,2);
t595 = Icges(6,6) + Icges(7,6);
t594 = Icges(6,3) + Icges(7,3);
t593 = rSges(7,3) + qJ(6);
t528 = sin(qJ(2));
t529 = sin(qJ(1));
t531 = cos(qJ(2));
t532 = cos(qJ(1));
t571 = cos(pkin(6));
t549 = t532 * t571;
t503 = t528 * t529 - t531 * t549;
t504 = t528 * t549 + t529 * t531;
t523 = sin(pkin(6));
t565 = t523 * t532;
t465 = Icges(3,5) * t504 - Icges(3,6) * t503 - Icges(3,3) * t565;
t550 = t529 * t571;
t505 = t532 * t528 + t531 * t550;
t506 = -t528 * t550 + t532 * t531;
t567 = t523 * t529;
t466 = Icges(3,5) * t506 - Icges(3,6) * t505 + Icges(3,3) * t567;
t592 = (t465 * t532 - t466 * t529) * t523;
t556 = pkin(11) + qJ(4);
t521 = sin(t556);
t547 = cos(t556);
t482 = t504 * t547 - t521 * t565;
t527 = sin(qJ(5));
t530 = cos(qJ(5));
t451 = -t482 * t527 + t503 * t530;
t570 = t503 * t527;
t452 = t482 * t530 + t570;
t541 = t523 * t547;
t481 = t504 * t521 + t532 * t541;
t591 = t595 * t451 + t597 * t452 + t594 * t481;
t484 = t506 * t547 + t521 * t567;
t453 = -t484 * t527 + t505 * t530;
t569 = t505 * t527;
t454 = t484 * t530 + t569;
t483 = t506 * t521 - t529 * t541;
t590 = t595 * t453 + t597 * t454 + t594 * t483;
t589 = t596 * t451 + t598 * t452 + t595 * t481;
t588 = t596 * t453 + t598 * t454 + t595 * t483;
t587 = t598 * t451 + t599 * t452 + t597 * t481;
t586 = t598 * t453 + t599 * t454 + t597 * t483;
t495 = t521 * t571 + t528 * t541;
t566 = t523 * t531;
t479 = -t495 * t527 - t530 * t566;
t553 = t527 * t566;
t480 = t495 * t530 - t553;
t568 = t523 * t528;
t494 = t521 * t568 - t547 * t571;
t585 = t595 * t479 + t597 * t480 + t594 * t494;
t584 = t596 * t479 + t598 * t480 + t595 * t494;
t583 = t598 * t479 + t599 * t480 + t597 * t494;
t522 = sin(pkin(11));
t524 = cos(pkin(11));
t485 = -t504 * t522 - t524 * t565;
t554 = t522 * t565;
t486 = t504 * t524 - t554;
t438 = Icges(4,5) * t486 + Icges(4,6) * t485 + Icges(4,3) * t503;
t467 = Icges(3,4) * t504 - Icges(3,2) * t503 - Icges(3,6) * t565;
t582 = -t438 + t467;
t487 = -t506 * t522 + t524 * t567;
t555 = t522 * t567;
t488 = t506 * t524 + t555;
t439 = Icges(4,5) * t488 + Icges(4,6) * t487 + Icges(4,3) * t505;
t468 = Icges(3,4) * t506 - Icges(3,2) * t505 + Icges(3,6) * t567;
t581 = t439 - t468;
t501 = -t522 * t568 + t524 * t571;
t548 = t571 * t522;
t502 = t524 * t568 + t548;
t459 = Icges(4,5) * t502 + Icges(4,6) * t501 - Icges(4,3) * t566;
t492 = Icges(3,6) * t571 + (Icges(3,4) * t528 + Icges(3,2) * t531) * t523;
t580 = t459 - t492;
t574 = pkin(3) * t524;
t573 = pkin(5) * t530;
t563 = rSges(7,1) * t452 + rSges(7,2) * t451 + pkin(5) * t570 + t593 * t481 + t482 * t573;
t562 = rSges(7,1) * t454 + rSges(7,2) * t453 + pkin(5) * t569 + t593 * t483 + t484 * t573;
t561 = rSges(7,1) * t480 + rSges(7,2) * t479 - pkin(5) * t553 + t593 * t494 + t495 * t573;
t477 = pkin(2) * t504 + qJ(3) * t503;
t478 = pkin(2) * t506 + qJ(3) * t505;
t558 = qJD(2) * t523;
t516 = t529 * t558;
t551 = t532 * t558;
t560 = t477 * t516 + t478 * t551;
t489 = qJD(4) * t505 + t516;
t559 = qJD(1) * (pkin(1) * t529 - pkin(8) * t565);
t557 = qJD(3) * t531;
t517 = qJD(2) * t571 + qJD(1);
t509 = qJD(1) * (pkin(1) * t532 + pkin(8) * t567);
t552 = qJD(3) * t503 + t517 * t478 + t509;
t546 = qJD(3) * t505 - t559;
t507 = (pkin(2) * t528 - qJ(3) * t531) * t523;
t543 = (-rSges(4,1) * t502 - rSges(4,2) * t501 + rSges(4,3) * t566 - t507) * t558;
t542 = (-pkin(3) * t548 - (-pkin(9) * t531 + t528 * t574) * t523 - t507) * t558;
t490 = qJD(4) * t503 - t551;
t508 = -qJD(4) * t566 + t517;
t435 = -pkin(3) * t554 + pkin(9) * t503 + t504 * t574;
t436 = pkin(3) * t555 + pkin(9) * t505 + t506 * t574;
t539 = t435 * t516 + t436 * t551 - t523 * t557 + t560;
t538 = t517 * t436 + t529 * t542 + t552;
t447 = pkin(4) * t482 + pkin(10) * t481;
t448 = pkin(4) * t484 + pkin(10) * t483;
t537 = t489 * t447 - t448 * t490 + t539;
t463 = pkin(4) * t495 + pkin(10) * t494;
t536 = t508 * t448 - t463 * t489 + t538;
t535 = (-t435 - t477) * t517 + t532 * t542 + t546;
t534 = -t447 * t508 + t490 * t463 + t535;
t513 = rSges(2,1) * t532 - rSges(2,2) * t529;
t512 = rSges(2,1) * t529 + rSges(2,2) * t532;
t496 = t571 * rSges(3,3) + (rSges(3,1) * t528 + rSges(3,2) * t531) * t523;
t493 = Icges(3,5) * t571 + (Icges(3,1) * t528 + Icges(3,4) * t531) * t523;
t491 = Icges(3,3) * t571 + (Icges(3,5) * t528 + Icges(3,6) * t531) * t523;
t476 = qJD(5) * t494 + t508;
t475 = rSges(3,1) * t506 - rSges(3,2) * t505 + rSges(3,3) * t567;
t474 = rSges(3,1) * t504 - rSges(3,2) * t503 - rSges(3,3) * t565;
t470 = Icges(3,1) * t506 - Icges(3,4) * t505 + Icges(3,5) * t567;
t469 = Icges(3,1) * t504 - Icges(3,4) * t503 - Icges(3,5) * t565;
t461 = Icges(4,1) * t502 + Icges(4,4) * t501 - Icges(4,5) * t566;
t460 = Icges(4,4) * t502 + Icges(4,2) * t501 - Icges(4,6) * t566;
t458 = rSges(5,1) * t495 - rSges(5,2) * t494 - rSges(5,3) * t566;
t457 = Icges(5,1) * t495 - Icges(5,4) * t494 - Icges(5,5) * t566;
t456 = Icges(5,4) * t495 - Icges(5,2) * t494 - Icges(5,6) * t566;
t455 = Icges(5,5) * t495 - Icges(5,6) * t494 - Icges(5,3) * t566;
t450 = qJD(5) * t481 + t490;
t449 = qJD(5) * t483 + t489;
t445 = rSges(4,1) * t488 + rSges(4,2) * t487 + rSges(4,3) * t505;
t444 = rSges(4,1) * t486 + rSges(4,2) * t485 + rSges(4,3) * t503;
t443 = Icges(4,1) * t488 + Icges(4,4) * t487 + Icges(4,5) * t505;
t442 = Icges(4,1) * t486 + Icges(4,4) * t485 + Icges(4,5) * t503;
t441 = Icges(4,4) * t488 + Icges(4,2) * t487 + Icges(4,6) * t505;
t440 = Icges(4,4) * t486 + Icges(4,2) * t485 + Icges(4,6) * t503;
t434 = rSges(5,1) * t484 - rSges(5,2) * t483 + rSges(5,3) * t505;
t433 = rSges(5,1) * t482 - rSges(5,2) * t481 + rSges(5,3) * t503;
t432 = Icges(5,1) * t484 - Icges(5,4) * t483 + Icges(5,5) * t505;
t431 = Icges(5,1) * t482 - Icges(5,4) * t481 + Icges(5,5) * t503;
t430 = Icges(5,4) * t484 - Icges(5,2) * t483 + Icges(5,6) * t505;
t429 = Icges(5,4) * t482 - Icges(5,2) * t481 + Icges(5,6) * t503;
t428 = Icges(5,5) * t484 - Icges(5,6) * t483 + Icges(5,3) * t505;
t427 = Icges(5,5) * t482 - Icges(5,6) * t481 + Icges(5,3) * t503;
t426 = t475 * t517 - t496 * t516 + t509;
t425 = -t474 * t517 - t496 * t551 - t559;
t424 = rSges(6,1) * t480 + rSges(6,2) * t479 + rSges(6,3) * t494;
t412 = (t474 * t529 + t475 * t532) * t558;
t410 = rSges(6,1) * t454 + rSges(6,2) * t453 + rSges(6,3) * t483;
t408 = rSges(6,1) * t452 + rSges(6,2) * t451 + rSges(6,3) * t481;
t392 = t445 * t517 + t529 * t543 + t552;
t391 = (-t444 - t477) * t517 + t532 * t543 + t546;
t390 = (-t557 + (t444 * t529 + t445 * t532) * qJD(2)) * t523 + t560;
t389 = t434 * t508 - t458 * t489 + t538;
t388 = -t433 * t508 + t458 * t490 + t535;
t387 = t433 * t489 - t434 * t490 + t539;
t386 = t410 * t476 - t424 * t449 + t536;
t385 = -t408 * t476 + t424 * t450 + t534;
t384 = t408 * t449 - t410 * t450 + t537;
t383 = qJD(6) * t481 - t449 * t561 + t476 * t562 + t536;
t382 = qJD(6) * t483 + t450 * t561 - t476 * t563 + t534;
t381 = qJD(6) * t494 + t449 * t563 - t450 * t562 + t537;
t1 = m(6) * (t384 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + m(7) * (t381 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + m(5) * (t387 ^ 2 + t388 ^ 2 + t389 ^ 2) / 0.2e1 + m(4) * (t390 ^ 2 + t391 ^ 2 + t392 ^ 2) / 0.2e1 + m(3) * (t412 ^ 2 + t425 ^ 2 + t426 ^ 2) / 0.2e1 + t489 * ((t505 * t428 - t483 * t430 + t484 * t432) * t489 + (t427 * t505 - t429 * t483 + t431 * t484) * t490 + (t455 * t505 - t456 * t483 + t457 * t484) * t508) / 0.2e1 + t490 * ((t428 * t503 - t430 * t481 + t432 * t482) * t489 + (t503 * t427 - t481 * t429 + t482 * t431) * t490 + (t455 * t503 - t456 * t481 + t457 * t482) * t508) / 0.2e1 + t508 * ((-t428 * t566 - t430 * t494 + t432 * t495) * t489 + (-t427 * t566 - t429 * t494 + t431 * t495) * t490 + (-t455 * t566 - t456 * t494 + t457 * t495) * t508) / 0.2e1 + ((t453 * t584 + t454 * t583 + t483 * t585) * t476 + (t453 * t589 + t454 * t587 + t483 * t591) * t450 + (t588 * t453 + t586 * t454 + t590 * t483) * t449) * t449 / 0.2e1 + ((t451 * t584 + t452 * t583 + t481 * t585) * t476 + (t589 * t451 + t587 * t452 + t591 * t481) * t450 + (t451 * t588 + t452 * t586 + t481 * t590) * t449) * t450 / 0.2e1 + ((t584 * t479 + t583 * t480 + t585 * t494) * t476 + (t479 * t589 + t480 * t587 + t494 * t591) * t450 + (t479 * t588 + t480 * t586 + t494 * t590) * t449) * t476 / 0.2e1 + (((t441 * t501 + t443 * t502) * t529 - (t440 * t501 + t442 * t502) * t532 + (t438 * t532 - t439 * t529) * t566) * t558 + (t571 * t466 + (t468 * t531 + t470 * t528) * t523) * t516 - (t571 * t465 + (t467 * t531 + t469 * t528) * t523) * t551 + (-t459 * t566 + t460 * t501 + t461 * t502 + t571 * t491 + (t492 * t531 + t493 * t528) * t523) * t517) * t517 / 0.2e1 + (m(2) * (t512 ^ 2 + t513 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t440 * t487 - t442 * t488 - t469 * t506 + t505 * t582) * t532 + (t441 * t487 + t443 * t488 + t506 * t470 + t505 * t581 - t592) * t529) * t558 + (t460 * t487 + t461 * t488 + t491 * t567 + t493 * t506 + t505 * t580) * t517) * t516 / 0.2e1 - (((-t440 * t485 - t442 * t486 - t504 * t469 + t503 * t582 + t592) * t532 + (t441 * t485 + t443 * t486 + t470 * t504 + t503 * t581) * t529) * t558 + (t460 * t485 + t461 * t486 - t491 * t565 + t493 * t504 + t503 * t580) * t517) * t551 / 0.2e1;
T  = t1;
