% Calculate kinetic energy for
% S6PRRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRP1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:22:37
% EndTime: 2019-03-08 21:22:40
% DurationCPUTime: 2.66s
% Computational Cost: add. (3169->324), mult. (5767->480), div. (0->0), fcn. (6983->12), ass. (0->153)
t595 = Icges(6,1) + Icges(7,1);
t594 = Icges(6,4) + Icges(7,4);
t593 = Icges(6,5) + Icges(7,5);
t592 = Icges(6,2) + Icges(7,2);
t591 = Icges(6,6) + Icges(7,6);
t590 = Icges(4,3) + Icges(5,3);
t589 = Icges(6,3) + Icges(7,3);
t588 = rSges(7,3) + qJ(6);
t521 = sin(pkin(10));
t523 = cos(pkin(10));
t532 = cos(qJ(2));
t524 = cos(pkin(6));
t529 = sin(qJ(2));
t555 = t524 * t529;
t506 = t521 * t532 + t523 * t555;
t548 = qJ(3) + pkin(11);
t520 = sin(t548);
t542 = cos(t548);
t522 = sin(pkin(6));
t561 = t522 * t523;
t484 = t506 * t542 - t520 * t561;
t554 = t524 * t532;
t505 = t521 * t529 - t523 * t554;
t527 = sin(qJ(5));
t530 = cos(qJ(5));
t455 = -t484 * t527 + t505 * t530;
t564 = t505 * t527;
t456 = t484 * t530 + t564;
t541 = t522 * t542;
t483 = t506 * t520 + t523 * t541;
t587 = t591 * t455 + t593 * t456 + t589 * t483;
t508 = -t521 * t555 + t523 * t532;
t562 = t521 * t522;
t486 = t508 * t542 + t520 * t562;
t507 = t521 * t554 + t523 * t529;
t457 = -t486 * t527 + t507 * t530;
t563 = t507 * t527;
t458 = t486 * t530 + t563;
t485 = t508 * t520 - t521 * t541;
t586 = t591 * t457 + t593 * t458 + t589 * t485;
t585 = t592 * t455 + t594 * t456 + t591 * t483;
t584 = t592 * t457 + t594 * t458 + t591 * t485;
t583 = t594 * t455 + t595 * t456 + t593 * t483;
t582 = t594 * t457 + t595 * t458 + t593 * t485;
t499 = t524 * t520 + t529 * t541;
t557 = t522 * t532;
t487 = -t499 * t527 - t530 * t557;
t545 = t527 * t557;
t488 = t499 * t530 - t545;
t559 = t522 * t529;
t498 = t520 * t559 - t524 * t542;
t581 = t591 * t487 + t593 * t488 + t589 * t498;
t580 = t592 * t487 + t594 * t488 + t591 * t498;
t579 = t594 * t487 + t595 * t488 + t593 * t498;
t528 = sin(qJ(3));
t531 = cos(qJ(3));
t558 = t522 * t531;
t489 = -t506 * t528 - t523 * t558;
t560 = t522 * t528;
t546 = t523 * t560;
t490 = t506 * t531 - t546;
t578 = Icges(4,5) * t490 + Icges(5,5) * t484 + Icges(4,6) * t489 - Icges(5,6) * t483 + t590 * t505;
t491 = -t508 * t528 + t521 * t558;
t547 = t521 * t560;
t492 = t508 * t531 + t547;
t577 = Icges(4,5) * t492 + Icges(5,5) * t486 + Icges(4,6) * t491 - Icges(5,6) * t485 + t590 * t507;
t509 = t524 * t531 - t528 * t559;
t556 = t524 * t528;
t510 = t529 * t558 + t556;
t576 = Icges(4,5) * t510 + Icges(5,5) * t499 + Icges(4,6) * t509 - Icges(5,6) * t498 - t590 * t557;
t575 = qJD(2) ^ 2;
t568 = pkin(3) * t531;
t567 = pkin(5) * t530;
t553 = rSges(7,1) * t456 + rSges(7,2) * t455 + pkin(5) * t564 + t483 * t588 + t484 * t567;
t552 = rSges(7,1) * t458 + rSges(7,2) * t457 + pkin(5) * t563 + t485 * t588 + t486 * t567;
t551 = rSges(7,1) * t488 + rSges(7,2) * t487 - pkin(5) * t545 + t498 * t588 + t499 * t567;
t550 = qJD(2) * t522;
t516 = t521 * t550;
t493 = qJD(3) * t507 + t516;
t519 = qJD(2) * t524;
t544 = t523 * t550;
t480 = pkin(2) * t506 + pkin(8) * t505;
t481 = pkin(2) * t508 + pkin(8) * t507;
t543 = t480 * t516 + t481 * t544 + qJD(1);
t494 = qJD(3) * t505 - t544;
t512 = -qJD(3) * t557 + t519;
t511 = (pkin(2) * t529 - pkin(8) * t532) * t522;
t540 = t481 * t519 - t511 * t516;
t437 = pkin(3) * t547 + qJ(4) * t507 + t508 * t568;
t539 = qJD(4) * t505 + t512 * t437 + t540;
t538 = (-t480 * t524 - t511 * t561) * qJD(2);
t436 = -pkin(3) * t546 + qJ(4) * t505 + t506 * t568;
t537 = -qJD(4) * t557 + t493 * t436 + t543;
t477 = pkin(3) * t556 + (-qJ(4) * t532 + t529 * t568) * t522;
t536 = qJD(4) * t507 + t494 * t477 + t538;
t451 = pkin(4) * t486 + pkin(9) * t485;
t469 = t499 * pkin(4) + t498 * pkin(9);
t535 = t512 * t451 + (-t469 - t477) * t493 + t539;
t450 = pkin(4) * t484 + pkin(9) * t483;
t534 = t493 * t450 + (-t437 - t451) * t494 + t537;
t533 = t494 * t469 + (-t436 - t450) * t512 + t536;
t500 = t524 * rSges(3,3) + (rSges(3,1) * t529 + rSges(3,2) * t532) * t522;
t497 = Icges(3,5) * t524 + (Icges(3,1) * t529 + Icges(3,4) * t532) * t522;
t496 = Icges(3,6) * t524 + (Icges(3,4) * t529 + Icges(3,2) * t532) * t522;
t495 = Icges(3,3) * t524 + (Icges(3,5) * t529 + Icges(3,6) * t532) * t522;
t482 = qJD(5) * t498 + t512;
t478 = t510 * rSges(4,1) + t509 * rSges(4,2) - rSges(4,3) * t557;
t476 = Icges(4,1) * t510 + Icges(4,4) * t509 - Icges(4,5) * t557;
t475 = Icges(4,4) * t510 + Icges(4,2) * t509 - Icges(4,6) * t557;
t471 = rSges(3,1) * t508 - rSges(3,2) * t507 + rSges(3,3) * t562;
t470 = rSges(3,1) * t506 - rSges(3,2) * t505 - rSges(3,3) * t561;
t468 = Icges(3,1) * t508 - Icges(3,4) * t507 + Icges(3,5) * t562;
t467 = Icges(3,1) * t506 - Icges(3,4) * t505 - Icges(3,5) * t561;
t466 = Icges(3,4) * t508 - Icges(3,2) * t507 + Icges(3,6) * t562;
t465 = Icges(3,4) * t506 - Icges(3,2) * t505 - Icges(3,6) * t561;
t464 = Icges(3,5) * t508 - Icges(3,6) * t507 + Icges(3,3) * t562;
t463 = Icges(3,5) * t506 - Icges(3,6) * t505 - Icges(3,3) * t561;
t462 = t499 * rSges(5,1) - t498 * rSges(5,2) - rSges(5,3) * t557;
t461 = Icges(5,1) * t499 - Icges(5,4) * t498 - Icges(5,5) * t557;
t460 = Icges(5,4) * t499 - Icges(5,2) * t498 - Icges(5,6) * t557;
t454 = qJD(5) * t483 + t494;
t453 = qJD(5) * t485 + t493;
t448 = (-t470 * t524 - t500 * t561) * qJD(2);
t447 = (t471 * t524 - t500 * t562) * qJD(2);
t446 = rSges(4,1) * t492 + rSges(4,2) * t491 + rSges(4,3) * t507;
t445 = rSges(4,1) * t490 + rSges(4,2) * t489 + rSges(4,3) * t505;
t444 = Icges(4,1) * t492 + Icges(4,4) * t491 + Icges(4,5) * t507;
t443 = Icges(4,1) * t490 + Icges(4,4) * t489 + Icges(4,5) * t505;
t442 = Icges(4,4) * t492 + Icges(4,2) * t491 + Icges(4,6) * t507;
t441 = Icges(4,4) * t490 + Icges(4,2) * t489 + Icges(4,6) * t505;
t435 = rSges(5,1) * t486 - rSges(5,2) * t485 + rSges(5,3) * t507;
t434 = rSges(5,1) * t484 - rSges(5,2) * t483 + rSges(5,3) * t505;
t433 = Icges(5,1) * t486 - Icges(5,4) * t485 + Icges(5,5) * t507;
t432 = Icges(5,1) * t484 - Icges(5,4) * t483 + Icges(5,5) * t505;
t431 = Icges(5,4) * t486 - Icges(5,2) * t485 + Icges(5,6) * t507;
t430 = Icges(5,4) * t484 - Icges(5,2) * t483 + Icges(5,6) * t505;
t427 = rSges(6,1) * t488 + rSges(6,2) * t487 + rSges(6,3) * t498;
t416 = qJD(1) + (t470 * t521 + t471 * t523) * t550;
t414 = rSges(6,1) * t458 + rSges(6,2) * t457 + rSges(6,3) * t485;
t412 = rSges(6,1) * t456 + rSges(6,2) * t455 + rSges(6,3) * t483;
t396 = -t445 * t512 + t478 * t494 + t538;
t395 = t446 * t512 - t478 * t493 + t540;
t394 = t445 * t493 - t446 * t494 + t543;
t393 = t462 * t494 + (-t434 - t436) * t512 + t536;
t392 = t435 * t512 + (-t462 - t477) * t493 + t539;
t391 = t493 * t434 + (-t435 - t437) * t494 + t537;
t390 = -t412 * t482 + t427 * t454 + t533;
t389 = t414 * t482 - t427 * t453 + t535;
t388 = t453 * t412 - t454 * t414 + t534;
t387 = qJD(6) * t485 + t454 * t551 - t482 * t553 + t533;
t386 = qJD(6) * t483 - t453 * t551 + t482 * t552 + t535;
t385 = qJD(6) * t498 + t453 * t553 - t454 * t552 + t534;
t1 = -t575 * ((-t464 * t561 - t466 * t505 + t468 * t506) * t562 - (-t463 * t561 - t465 * t505 + t467 * t506) * t561 + (-t495 * t561 - t496 * t505 + t497 * t506) * t524) * t561 / 0.2e1 + m(6) * (t388 ^ 2 + t389 ^ 2 + t390 ^ 2) / 0.2e1 + m(7) * (t385 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t416 ^ 2 + t447 ^ 2 + t448 ^ 2) / 0.2e1 + m(4) * (t394 ^ 2 + t395 ^ 2 + t396 ^ 2) / 0.2e1 + m(5) * (t391 ^ 2 + t392 ^ 2 + t393 ^ 2) / 0.2e1 + ((t457 * t580 + t458 * t579 + t485 * t581) * t482 + (t585 * t457 + t583 * t458 + t485 * t587) * t454 + (t457 * t584 + t458 * t582 + t485 * t586) * t453) * t453 / 0.2e1 + ((t455 * t580 + t456 * t579 + t483 * t581) * t482 + (t585 * t455 + t583 * t456 + t483 * t587) * t454 + (t455 * t584 + t456 * t582 + t483 * t586) * t453) * t454 / 0.2e1 + ((t487 * t580 + t579 * t488 + t581 * t498) * t482 + (t585 * t487 + t583 * t488 + t498 * t587) * t454 + (t487 * t584 + t488 * t582 + t498 * t586) * t453) * t482 / 0.2e1 + ((-t460 * t485 + t461 * t486 + t475 * t491 + t476 * t492 + t507 * t576) * t512 + (-t430 * t485 + t432 * t486 + t441 * t491 + t443 * t492 + t507 * t578) * t494 + (-t431 * t485 + t433 * t486 + t442 * t491 + t444 * t492 + t507 * t577) * t493) * t493 / 0.2e1 + ((-t460 * t483 + t461 * t484 + t475 * t489 + t476 * t490 + t505 * t576) * t512 + (-t483 * t430 + t432 * t484 + t441 * t489 + t443 * t490 + t505 * t578) * t494 + (-t431 * t483 + t433 * t484 + t442 * t489 + t444 * t490 + t505 * t577) * t493) * t494 / 0.2e1 + ((-t498 * t460 + t499 * t461 + t509 * t475 + t510 * t476 - t557 * t576) * t512 + (-t498 * t430 + t499 * t432 + t509 * t441 + t510 * t443 - t557 * t578) * t494 + (-t498 * t431 + t499 * t433 + t509 * t442 + t510 * t444 - t557 * t577) * t493) * t512 / 0.2e1 + (((t464 * t562 - t466 * t507 + t468 * t508) * t562 - (t463 * t562 - t465 * t507 + t467 * t508) * t561 + (t495 * t562 - t496 * t507 + t497 * t508) * t524) * t562 + t524 * (t524 ^ 2 * t495 + (((t466 * t532 + t468 * t529) * t521 - (t465 * t532 + t467 * t529) * t523) * t522 + (-t463 * t523 + t464 * t521 + t496 * t532 + t497 * t529) * t524) * t522)) * t575 / 0.2e1;
T  = t1;
