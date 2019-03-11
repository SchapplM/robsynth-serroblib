% Calculate kinetic energy for
% S6RRRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP7_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:34:28
% EndTime: 2019-03-10 01:34:31
% DurationCPUTime: 3.08s
% Computational Cost: add. (3343->338), mult. (6029->518), div. (0->0), fcn. (7269->12), ass. (0->161)
t592 = Icges(6,1) + Icges(7,1);
t591 = Icges(6,4) + Icges(7,4);
t590 = Icges(6,5) + Icges(7,5);
t589 = Icges(6,2) + Icges(7,2);
t588 = Icges(6,6) + Icges(7,6);
t587 = Icges(6,3) + Icges(7,3);
t586 = rSges(7,3) + qJ(6);
t527 = sin(qJ(2));
t528 = sin(qJ(1));
t531 = cos(qJ(2));
t532 = cos(qJ(1));
t567 = cos(pkin(6));
t546 = t532 * t567;
t504 = t527 * t528 - t531 * t546;
t505 = t527 * t546 + t528 * t531;
t523 = sin(pkin(6));
t560 = t523 * t532;
t465 = Icges(3,5) * t505 - Icges(3,6) * t504 - Icges(3,3) * t560;
t547 = t528 * t567;
t506 = t532 * t527 + t531 * t547;
t507 = -t527 * t547 + t532 * t531;
t563 = t523 * t528;
t466 = Icges(3,5) * t507 - Icges(3,6) * t506 + Icges(3,3) * t563;
t585 = t523 * (t465 * t532 - t466 * t528);
t559 = qJ(3) + qJ(4);
t522 = sin(t559);
t548 = cos(t559);
t482 = t505 * t548 - t522 * t560;
t525 = sin(qJ(5));
t529 = cos(qJ(5));
t449 = -t482 * t525 + t504 * t529;
t566 = t504 * t525;
t450 = t482 * t529 + t566;
t544 = t523 * t548;
t481 = t505 * t522 + t532 * t544;
t584 = t588 * t449 + t590 * t450 + t587 * t481;
t484 = t507 * t548 + t522 * t563;
t451 = -t484 * t525 + t506 * t529;
t565 = t506 * t525;
t452 = t484 * t529 + t565;
t483 = t507 * t522 - t528 * t544;
t583 = t588 * t451 + t590 * t452 + t587 * t483;
t582 = t589 * t449 + t591 * t450 + t588 * t481;
t581 = t589 * t451 + t591 * t452 + t588 * t483;
t580 = t591 * t449 + t592 * t450 + t590 * t481;
t579 = t591 * t451 + t592 * t452 + t590 * t483;
t497 = t522 * t567 + t527 * t544;
t561 = t523 * t531;
t479 = -t497 * t525 - t529 * t561;
t552 = t525 * t561;
t480 = t497 * t529 - t552;
t564 = t523 * t527;
t496 = t522 * t564 - t548 * t567;
t578 = t588 * t479 + t590 * t480 + t587 * t496;
t577 = t589 * t479 + t591 * t480 + t588 * t496;
t576 = t591 * t479 + t592 * t480 + t590 * t496;
t530 = cos(qJ(3));
t571 = pkin(3) * t530;
t570 = pkin(5) * t529;
t562 = t523 * t530;
t558 = rSges(7,1) * t450 + rSges(7,2) * t449 + pkin(5) * t566 + t586 * t481 + t482 * t570;
t557 = rSges(7,1) * t452 + rSges(7,2) * t451 + pkin(5) * t565 + t586 * t483 + t484 * t570;
t556 = rSges(7,1) * t480 + rSges(7,2) * t479 - pkin(5) * t552 + t586 * t496 + t497 * t570;
t477 = pkin(2) * t505 + pkin(9) * t504;
t478 = pkin(2) * t507 + pkin(9) * t506;
t553 = qJD(2) * t523;
t517 = t528 * t553;
t549 = t532 * t553;
t555 = t477 * t517 + t478 * t549;
t489 = qJD(3) * t506 + t517;
t554 = qJD(1) * (pkin(1) * t528 - pkin(8) * t560);
t518 = qJD(2) * t567 + qJD(1);
t526 = sin(qJ(3));
t551 = t526 * t563;
t550 = t526 * t560;
t461 = qJD(4) * t506 + t489;
t545 = t567 * t526;
t490 = qJD(3) * t504 - t549;
t434 = -pkin(3) * t550 + pkin(10) * t504 + t505 * t571;
t435 = pkin(3) * t551 + pkin(10) * t506 + t507 * t571;
t542 = t489 * t434 - t435 * t490 + t555;
t462 = qJD(4) * t504 + t490;
t508 = (pkin(2) * t527 - pkin(9) * t531) * t523;
t510 = qJD(1) * (pkin(1) * t532 + pkin(8) * t563);
t541 = t518 * t478 - t508 * t517 + t510;
t491 = (-qJD(3) - qJD(4)) * t561 + t518;
t447 = pkin(4) * t482 + pkin(11) * t481;
t448 = pkin(4) * t484 + pkin(11) * t483;
t540 = t461 * t447 - t448 * t462 + t542;
t539 = -t477 * t518 - t508 * t549 - t554;
t471 = pkin(3) * t545 + (-pkin(10) * t531 + t527 * t571) * t523;
t509 = -qJD(3) * t561 + t518;
t538 = t509 * t435 - t471 * t489 + t541;
t537 = -t434 * t509 + t490 * t471 + t539;
t464 = pkin(4) * t497 + pkin(11) * t496;
t536 = t491 * t448 - t461 * t464 + t538;
t535 = -t447 * t491 + t462 * t464 + t537;
t514 = rSges(2,1) * t532 - rSges(2,2) * t528;
t513 = rSges(2,1) * t528 + rSges(2,2) * t532;
t503 = t527 * t562 + t545;
t502 = -t526 * t564 + t530 * t567;
t495 = t567 * rSges(3,3) + (rSges(3,1) * t527 + rSges(3,2) * t531) * t523;
t494 = Icges(3,5) * t567 + (Icges(3,1) * t527 + Icges(3,4) * t531) * t523;
t493 = Icges(3,6) * t567 + (Icges(3,4) * t527 + Icges(3,2) * t531) * t523;
t492 = Icges(3,3) * t567 + (Icges(3,5) * t527 + Icges(3,6) * t531) * t523;
t488 = t507 * t530 + t551;
t487 = -t507 * t526 + t528 * t562;
t486 = t505 * t530 - t550;
t485 = -t505 * t526 - t530 * t560;
t474 = rSges(3,1) * t507 - rSges(3,2) * t506 + rSges(3,3) * t563;
t473 = rSges(3,1) * t505 - rSges(3,2) * t504 - rSges(3,3) * t560;
t470 = Icges(3,1) * t507 - Icges(3,4) * t506 + Icges(3,5) * t563;
t469 = Icges(3,1) * t505 - Icges(3,4) * t504 - Icges(3,5) * t560;
t468 = Icges(3,4) * t507 - Icges(3,2) * t506 + Icges(3,6) * t563;
t467 = Icges(3,4) * t505 - Icges(3,2) * t504 - Icges(3,6) * t560;
t463 = rSges(4,1) * t503 + rSges(4,2) * t502 - rSges(4,3) * t561;
t460 = Icges(4,1) * t503 + Icges(4,4) * t502 - Icges(4,5) * t561;
t459 = Icges(4,4) * t503 + Icges(4,2) * t502 - Icges(4,6) * t561;
t458 = Icges(4,5) * t503 + Icges(4,6) * t502 - Icges(4,3) * t561;
t457 = qJD(5) * t496 + t491;
t456 = rSges(5,1) * t497 - rSges(5,2) * t496 - rSges(5,3) * t561;
t455 = Icges(5,1) * t497 - Icges(5,4) * t496 - Icges(5,5) * t561;
t454 = Icges(5,4) * t497 - Icges(5,2) * t496 - Icges(5,6) * t561;
t453 = Icges(5,5) * t497 - Icges(5,6) * t496 - Icges(5,3) * t561;
t445 = qJD(5) * t481 + t462;
t444 = qJD(5) * t483 + t461;
t443 = rSges(4,1) * t488 + rSges(4,2) * t487 + rSges(4,3) * t506;
t442 = rSges(4,1) * t486 + rSges(4,2) * t485 + rSges(4,3) * t504;
t441 = Icges(4,1) * t488 + Icges(4,4) * t487 + Icges(4,5) * t506;
t440 = Icges(4,1) * t486 + Icges(4,4) * t485 + Icges(4,5) * t504;
t439 = Icges(4,4) * t488 + Icges(4,2) * t487 + Icges(4,6) * t506;
t438 = Icges(4,4) * t486 + Icges(4,2) * t485 + Icges(4,6) * t504;
t437 = Icges(4,5) * t488 + Icges(4,6) * t487 + Icges(4,3) * t506;
t436 = Icges(4,5) * t486 + Icges(4,6) * t485 + Icges(4,3) * t504;
t433 = rSges(5,1) * t484 - rSges(5,2) * t483 + rSges(5,3) * t506;
t432 = rSges(5,1) * t482 - rSges(5,2) * t481 + rSges(5,3) * t504;
t431 = Icges(5,1) * t484 - Icges(5,4) * t483 + Icges(5,5) * t506;
t430 = Icges(5,1) * t482 - Icges(5,4) * t481 + Icges(5,5) * t504;
t429 = Icges(5,4) * t484 - Icges(5,2) * t483 + Icges(5,6) * t506;
t428 = Icges(5,4) * t482 - Icges(5,2) * t481 + Icges(5,6) * t504;
t427 = Icges(5,5) * t484 - Icges(5,6) * t483 + Icges(5,3) * t506;
t426 = Icges(5,5) * t482 - Icges(5,6) * t481 + Icges(5,3) * t504;
t424 = rSges(6,1) * t480 + rSges(6,2) * t479 + rSges(6,3) * t496;
t415 = t474 * t518 - t495 * t517 + t510;
t414 = -t473 * t518 - t495 * t549 - t554;
t412 = (t473 * t528 + t474 * t532) * t553;
t408 = rSges(6,1) * t452 + rSges(6,2) * t451 + rSges(6,3) * t483;
t406 = rSges(6,1) * t450 + rSges(6,2) * t449 + rSges(6,3) * t481;
t390 = t443 * t509 - t463 * t489 + t541;
t389 = -t442 * t509 + t463 * t490 + t539;
t388 = t442 * t489 - t443 * t490 + t555;
t387 = t433 * t491 - t456 * t461 + t538;
t386 = -t432 * t491 + t456 * t462 + t537;
t385 = t432 * t461 - t433 * t462 + t542;
t384 = t408 * t457 - t424 * t444 + t536;
t383 = -t406 * t457 + t424 * t445 + t535;
t382 = t406 * t444 - t408 * t445 + t540;
t381 = qJD(6) * t481 - t444 * t556 + t457 * t557 + t536;
t380 = qJD(6) * t483 + t445 * t556 - t457 * t558 + t535;
t379 = qJD(6) * t496 + t444 * t558 - t445 * t557 + t540;
t1 = m(5) * (t385 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + m(4) * (t388 ^ 2 + t389 ^ 2 + t390 ^ 2) / 0.2e1 + m(3) * (t412 ^ 2 + t414 ^ 2 + t415 ^ 2) / 0.2e1 + t491 * ((-t427 * t561 - t496 * t429 + t497 * t431) * t461 + (-t426 * t561 - t428 * t496 + t430 * t497) * t462 + (-t453 * t561 - t496 * t454 + t497 * t455) * t491) / 0.2e1 + t461 * ((t506 * t427 - t483 * t429 + t484 * t431) * t461 + (t426 * t506 - t428 * t483 + t430 * t484) * t462 + (t453 * t506 - t454 * t483 + t455 * t484) * t491) / 0.2e1 + m(7) * (t379 ^ 2 + t380 ^ 2 + t381 ^ 2) / 0.2e1 + m(6) * (t382 ^ 2 + t383 ^ 2 + t384 ^ 2) / 0.2e1 + t462 * ((t427 * t504 - t429 * t481 + t431 * t482) * t461 + (t504 * t426 - t481 * t428 + t482 * t430) * t462 + (t453 * t504 - t454 * t481 + t455 * t482) * t491) / 0.2e1 + t490 * ((t437 * t504 + t439 * t485 + t441 * t486) * t489 + (t504 * t436 + t485 * t438 + t486 * t440) * t490 + (t458 * t504 + t459 * t485 + t460 * t486) * t509) / 0.2e1 + t509 * ((-t437 * t561 + t439 * t502 + t441 * t503) * t489 + (-t436 * t561 + t438 * t502 + t440 * t503) * t490 + (-t458 * t561 + t502 * t459 + t503 * t460) * t509) / 0.2e1 + t489 * ((t506 * t437 + t487 * t439 + t488 * t441) * t489 + (t436 * t506 + t438 * t487 + t440 * t488) * t490 + (t458 * t506 + t459 * t487 + t460 * t488) * t509) / 0.2e1 + t518 * ((t567 * t466 + (t468 * t531 + t470 * t527) * t523) * t517 - (t567 * t465 + (t467 * t531 + t469 * t527) * t523) * t549 + (t567 * t492 + (t493 * t531 + t494 * t527) * t523) * t518) / 0.2e1 - ((-t492 * t560 - t493 * t504 + t494 * t505) * t518 + ((-t468 * t504 + t470 * t505) * t528 + (t504 * t467 - t505 * t469 + t585) * t532) * t553) * t549 / 0.2e1 + ((t492 * t563 - t493 * t506 + t494 * t507) * t518 + (-(-t467 * t506 + t469 * t507) * t532 + (-t506 * t468 + t507 * t470 - t585) * t528) * t553) * t517 / 0.2e1 + ((t451 * t577 + t452 * t576 + t483 * t578) * t457 + (t451 * t582 + t452 * t580 + t483 * t584) * t445 + (t451 * t581 + t452 * t579 + t483 * t583) * t444) * t444 / 0.2e1 + ((t449 * t577 + t450 * t576 + t481 * t578) * t457 + (t449 * t582 + t450 * t580 + t481 * t584) * t445 + (t449 * t581 + t450 * t579 + t481 * t583) * t444) * t445 / 0.2e1 + ((t479 * t577 + t480 * t576 + t496 * t578) * t457 + (t479 * t582 + t480 * t580 + t496 * t584) * t445 + (t479 * t581 + t480 * t579 + t496 * t583) * t444) * t457 / 0.2e1 + (m(2) * (t513 ^ 2 + t514 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
