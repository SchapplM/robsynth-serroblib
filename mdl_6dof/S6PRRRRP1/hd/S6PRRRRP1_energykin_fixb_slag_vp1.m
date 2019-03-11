% Calculate kinetic energy for
% S6PRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP1_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:55:43
% EndTime: 2019-03-08 23:55:46
% DurationCPUTime: 2.75s
% Computational Cost: add. (3274->330), mult. (5977->510), div. (0->0), fcn. (7235->12), ass. (0->158)
t587 = Icges(6,1) + Icges(7,1);
t586 = Icges(6,4) + Icges(7,4);
t585 = Icges(6,5) + Icges(7,5);
t584 = Icges(6,2) + Icges(7,2);
t583 = Icges(6,6) + Icges(7,6);
t582 = Icges(6,3) + Icges(7,3);
t581 = rSges(7,3) + qJ(6);
t520 = sin(pkin(11));
t522 = cos(pkin(11));
t530 = cos(qJ(2));
t523 = cos(pkin(6));
t527 = sin(qJ(2));
t554 = t523 * t527;
t505 = t520 * t530 + t522 * t554;
t552 = qJ(3) + qJ(4);
t519 = sin(t552);
t541 = cos(t552);
t521 = sin(pkin(6));
t560 = t521 * t522;
t482 = t505 * t541 - t519 * t560;
t553 = t523 * t530;
t504 = t520 * t527 - t522 * t553;
t525 = sin(qJ(5));
t528 = cos(qJ(5));
t451 = -t482 * t525 + t504 * t528;
t563 = t504 * t525;
t452 = t482 * t528 + t563;
t540 = t521 * t541;
t481 = t505 * t519 + t522 * t540;
t580 = t583 * t451 + t585 * t452 + t582 * t481;
t507 = -t520 * t554 + t522 * t530;
t561 = t520 * t521;
t484 = t507 * t541 + t519 * t561;
t506 = t520 * t553 + t522 * t527;
t453 = -t484 * t525 + t506 * t528;
t562 = t506 * t525;
t454 = t484 * t528 + t562;
t483 = t507 * t519 - t520 * t540;
t579 = t583 * t453 + t585 * t454 + t582 * t483;
t578 = t584 * t451 + t586 * t452 + t583 * t481;
t577 = t584 * t453 + t586 * t454 + t583 * t483;
t576 = t586 * t451 + t587 * t452 + t585 * t481;
t575 = t586 * t453 + t587 * t454 + t585 * t483;
t499 = t523 * t519 + t527 * t540;
t556 = t521 * t530;
t485 = -t499 * t525 - t528 * t556;
t544 = t525 * t556;
t486 = t499 * t528 - t544;
t558 = t521 * t527;
t498 = t519 * t558 - t523 * t541;
t574 = t583 * t485 + t585 * t486 + t582 * t498;
t573 = t584 * t485 + t586 * t486 + t583 * t498;
t572 = t586 * t485 + t587 * t486 + t585 * t498;
t571 = qJD(2) ^ 2;
t529 = cos(qJ(3));
t567 = pkin(3) * t529;
t566 = pkin(5) * t528;
t526 = sin(qJ(3));
t559 = t521 * t526;
t557 = t521 * t529;
t555 = t523 * t526;
t551 = rSges(7,1) * t452 + rSges(7,2) * t451 + pkin(5) * t563 + t481 * t581 + t482 * t566;
t550 = rSges(7,1) * t454 + rSges(7,2) * t453 + pkin(5) * t562 + t483 * t581 + t484 * t566;
t549 = rSges(7,1) * t486 + rSges(7,2) * t485 - pkin(5) * t544 + t498 * t581 + t499 * t566;
t548 = qJD(2) * t521;
t515 = t520 * t548;
t491 = qJD(3) * t506 + t515;
t518 = qJD(2) * t523;
t546 = t520 * t559;
t545 = t522 * t559;
t459 = qJD(4) * t506 + t491;
t543 = t522 * t548;
t479 = t505 * pkin(2) + t504 * pkin(8);
t480 = t507 * pkin(2) + t506 * pkin(8);
t542 = t479 * t515 + t480 * t543 + qJD(1);
t492 = qJD(3) * t504 - t543;
t510 = (pkin(2) * t527 - pkin(8) * t530) * t521;
t539 = t480 * t518 - t510 * t515;
t460 = qJD(4) * t504 + t492;
t493 = t518 + (-qJD(3) - qJD(4)) * t556;
t434 = -pkin(3) * t545 + pkin(9) * t504 + t505 * t567;
t435 = pkin(3) * t546 + pkin(9) * t506 + t507 * t567;
t538 = t491 * t434 - t435 * t492 + t542;
t537 = (-t479 * t523 - t510 * t560) * qJD(2);
t476 = pkin(3) * t555 + (-pkin(9) * t530 + t527 * t567) * t521;
t511 = -qJD(3) * t556 + t518;
t536 = t511 * t435 - t476 * t491 + t539;
t449 = pkin(4) * t482 + pkin(10) * t481;
t450 = pkin(4) * t484 + pkin(10) * t483;
t535 = t459 * t449 - t450 * t460 + t538;
t534 = -t434 * t511 + t492 * t476 + t537;
t475 = pkin(4) * t499 + pkin(10) * t498;
t533 = t493 * t450 - t459 * t475 + t536;
t532 = -t449 * t493 + t460 * t475 + t534;
t509 = t527 * t557 + t555;
t508 = t523 * t529 - t526 * t558;
t497 = rSges(3,3) * t523 + (rSges(3,1) * t527 + rSges(3,2) * t530) * t521;
t496 = Icges(3,5) * t523 + (Icges(3,1) * t527 + Icges(3,4) * t530) * t521;
t495 = Icges(3,6) * t523 + (Icges(3,4) * t527 + Icges(3,2) * t530) * t521;
t494 = Icges(3,3) * t523 + (Icges(3,5) * t527 + Icges(3,6) * t530) * t521;
t490 = t507 * t529 + t546;
t489 = -t507 * t526 + t520 * t557;
t488 = t505 * t529 - t545;
t487 = -t505 * t526 - t522 * t557;
t477 = qJD(5) * t498 + t493;
t474 = rSges(4,1) * t509 + rSges(4,2) * t508 - rSges(4,3) * t556;
t473 = Icges(4,1) * t509 + Icges(4,4) * t508 - Icges(4,5) * t556;
t472 = Icges(4,4) * t509 + Icges(4,2) * t508 - Icges(4,6) * t556;
t471 = Icges(4,5) * t509 + Icges(4,6) * t508 - Icges(4,3) * t556;
t468 = rSges(3,1) * t507 - rSges(3,2) * t506 + rSges(3,3) * t561;
t467 = rSges(3,1) * t505 - rSges(3,2) * t504 - rSges(3,3) * t560;
t466 = Icges(3,1) * t507 - Icges(3,4) * t506 + Icges(3,5) * t561;
t465 = Icges(3,1) * t505 - Icges(3,4) * t504 - Icges(3,5) * t560;
t464 = Icges(3,4) * t507 - Icges(3,2) * t506 + Icges(3,6) * t561;
t463 = Icges(3,4) * t505 - Icges(3,2) * t504 - Icges(3,6) * t560;
t462 = Icges(3,5) * t507 - Icges(3,6) * t506 + Icges(3,3) * t561;
t461 = Icges(3,5) * t505 - Icges(3,6) * t504 - Icges(3,3) * t560;
t458 = rSges(5,1) * t499 - rSges(5,2) * t498 - rSges(5,3) * t556;
t457 = Icges(5,1) * t499 - Icges(5,4) * t498 - Icges(5,5) * t556;
t456 = Icges(5,4) * t499 - Icges(5,2) * t498 - Icges(5,6) * t556;
t455 = Icges(5,5) * t499 - Icges(5,6) * t498 - Icges(5,3) * t556;
t447 = (-t467 * t523 - t497 * t560) * qJD(2);
t446 = (t468 * t523 - t497 * t561) * qJD(2);
t445 = qJD(5) * t481 + t460;
t444 = qJD(5) * t483 + t459;
t443 = rSges(4,1) * t490 + rSges(4,2) * t489 + rSges(4,3) * t506;
t442 = rSges(4,1) * t488 + rSges(4,2) * t487 + rSges(4,3) * t504;
t441 = Icges(4,1) * t490 + Icges(4,4) * t489 + Icges(4,5) * t506;
t440 = Icges(4,1) * t488 + Icges(4,4) * t487 + Icges(4,5) * t504;
t439 = Icges(4,4) * t490 + Icges(4,2) * t489 + Icges(4,6) * t506;
t438 = Icges(4,4) * t488 + Icges(4,2) * t487 + Icges(4,6) * t504;
t437 = Icges(4,5) * t490 + Icges(4,6) * t489 + Icges(4,3) * t506;
t436 = Icges(4,5) * t488 + Icges(4,6) * t487 + Icges(4,3) * t504;
t433 = rSges(5,1) * t484 - rSges(5,2) * t483 + rSges(5,3) * t506;
t432 = rSges(5,1) * t482 - rSges(5,2) * t481 + rSges(5,3) * t504;
t430 = Icges(5,1) * t484 - Icges(5,4) * t483 + Icges(5,5) * t506;
t429 = Icges(5,1) * t482 - Icges(5,4) * t481 + Icges(5,5) * t504;
t428 = Icges(5,4) * t484 - Icges(5,2) * t483 + Icges(5,6) * t506;
t427 = Icges(5,4) * t482 - Icges(5,2) * t481 + Icges(5,6) * t504;
t426 = Icges(5,5) * t484 - Icges(5,6) * t483 + Icges(5,3) * t506;
t425 = Icges(5,5) * t482 - Icges(5,6) * t481 + Icges(5,3) * t504;
t424 = rSges(6,1) * t486 + rSges(6,2) * t485 + rSges(6,3) * t498;
t413 = qJD(1) + (t467 * t520 + t468 * t522) * t548;
t410 = rSges(6,1) * t454 + rSges(6,2) * t453 + rSges(6,3) * t483;
t408 = rSges(6,1) * t452 + rSges(6,2) * t451 + rSges(6,3) * t481;
t392 = -t442 * t511 + t474 * t492 + t537;
t391 = t443 * t511 - t474 * t491 + t539;
t390 = t442 * t491 - t443 * t492 + t542;
t389 = -t432 * t493 + t458 * t460 + t534;
t388 = t433 * t493 - t458 * t459 + t536;
t387 = t432 * t459 - t433 * t460 + t538;
t386 = -t408 * t477 + t424 * t445 + t532;
t385 = t410 * t477 - t424 * t444 + t533;
t384 = t408 * t444 - t410 * t445 + t535;
t383 = qJD(6) * t483 + t445 * t549 - t477 * t551 + t532;
t382 = qJD(6) * t481 - t444 * t549 + t477 * t550 + t533;
t381 = qJD(6) * t498 + t444 * t551 - t445 * t550 + t535;
t1 = -t571 * ((-t462 * t560 - t464 * t504 + t466 * t505) * t561 - (-t461 * t560 - t463 * t504 + t465 * t505) * t560 + (-t494 * t560 - t495 * t504 + t496 * t505) * t523) * t560 / 0.2e1 + m(6) * (t384 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + m(5) * (t387 ^ 2 + t388 ^ 2 + t389 ^ 2) / 0.2e1 + t492 * ((t437 * t504 + t439 * t487 + t441 * t488) * t491 + (t504 * t436 + t487 * t438 + t488 * t440) * t492 + (t471 * t504 + t472 * t487 + t473 * t488) * t511) / 0.2e1 + t511 * ((-t437 * t556 + t439 * t508 + t441 * t509) * t491 + (-t436 * t556 + t438 * t508 + t440 * t509) * t492 + (-t471 * t556 + t472 * t508 + t473 * t509) * t511) / 0.2e1 + t493 * ((-t426 * t556 - t428 * t498 + t430 * t499) * t459 + (-t425 * t556 - t427 * t498 + t429 * t499) * t460 + (-t455 * t556 - t456 * t498 + t457 * t499) * t493) / 0.2e1 + t459 * ((t506 * t426 - t483 * t428 + t484 * t430) * t459 + (t425 * t506 - t427 * t483 + t429 * t484) * t460 + (t455 * t506 - t456 * t483 + t457 * t484) * t493) / 0.2e1 + t460 * ((t426 * t504 - t428 * t481 + t430 * t482) * t459 + (t504 * t425 - t481 * t427 + t482 * t429) * t460 + (t455 * t504 - t456 * t481 + t457 * t482) * t493) / 0.2e1 + t491 * ((t506 * t437 + t489 * t439 + t490 * t441) * t491 + (t436 * t506 + t438 * t489 + t440 * t490) * t492 + (t471 * t506 + t472 * t489 + t473 * t490) * t511) / 0.2e1 + m(4) * (t390 ^ 2 + t391 ^ 2 + t392 ^ 2) / 0.2e1 + m(3) * (t413 ^ 2 + t446 ^ 2 + t447 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(7) * (t381 ^ 2 + t382 ^ 2 + t383 ^ 2) / 0.2e1 + ((t453 * t573 + t454 * t572 + t483 * t574) * t477 + (t578 * t453 + t576 * t454 + t483 * t580) * t445 + (t577 * t453 + t575 * t454 + t579 * t483) * t444) * t444 / 0.2e1 + ((t451 * t573 + t452 * t572 + t481 * t574) * t477 + (t578 * t451 + t576 * t452 + t580 * t481) * t445 + (t451 * t577 + t452 * t575 + t481 * t579) * t444) * t445 / 0.2e1 + ((t573 * t485 + t572 * t486 + t574 * t498) * t477 + (t578 * t485 + t576 * t486 + t498 * t580) * t445 + (t485 * t577 + t486 * t575 + t498 * t579) * t444) * t477 / 0.2e1 + (((t462 * t561 - t464 * t506 + t466 * t507) * t561 - (t461 * t561 - t463 * t506 + t465 * t507) * t560 + (t494 * t561 - t495 * t506 + t496 * t507) * t523) * t561 + t523 * (t523 ^ 2 * t494 + (((t464 * t530 + t466 * t527) * t520 - (t463 * t530 + t465 * t527) * t522) * t521 + (-t461 * t522 + t462 * t520 + t495 * t530 + t496 * t527) * t523) * t521)) * t571 / 0.2e1;
T  = t1;
