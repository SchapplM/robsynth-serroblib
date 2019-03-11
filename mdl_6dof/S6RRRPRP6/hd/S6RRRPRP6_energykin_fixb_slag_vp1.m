% Calculate kinetic energy for
% S6RRRPRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 17:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP6_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRP6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRP6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:53:51
% EndTime: 2019-03-09 16:53:54
% DurationCPUTime: 2.96s
% Computational Cost: add. (3238->332), mult. (5819->488), div. (0->0), fcn. (7017->12), ass. (0->156)
t600 = Icges(6,1) + Icges(7,1);
t599 = Icges(6,4) + Icges(7,4);
t598 = Icges(6,5) + Icges(7,5);
t597 = Icges(6,2) + Icges(7,2);
t596 = Icges(6,6) + Icges(7,6);
t595 = Icges(4,3) + Icges(5,3);
t594 = Icges(6,3) + Icges(7,3);
t593 = rSges(7,3) + qJ(6);
t529 = sin(qJ(2));
t530 = sin(qJ(1));
t533 = cos(qJ(2));
t534 = cos(qJ(1));
t568 = cos(pkin(6));
t548 = t534 * t568;
t505 = t529 * t530 - t533 * t548;
t506 = t529 * t548 + t530 * t533;
t524 = sin(pkin(6));
t561 = t524 * t534;
t467 = Icges(3,5) * t506 - Icges(3,6) * t505 - Icges(3,3) * t561;
t549 = t530 * t568;
t507 = t534 * t529 + t533 * t549;
t508 = -t529 * t549 + t534 * t533;
t564 = t524 * t530;
t468 = Icges(3,5) * t508 - Icges(3,6) * t507 + Icges(3,3) * t564;
t592 = t524 * (t467 * t534 - t468 * t530);
t554 = qJ(3) + pkin(11);
t523 = sin(t554);
t546 = cos(t554);
t484 = t506 * t546 - t523 * t561;
t527 = sin(qJ(5));
t531 = cos(qJ(5));
t453 = -t484 * t527 + t505 * t531;
t567 = t505 * t527;
t454 = t484 * t531 + t567;
t545 = t524 * t546;
t483 = t506 * t523 + t534 * t545;
t591 = t596 * t453 + t598 * t454 + t594 * t483;
t486 = t508 * t546 + t523 * t564;
t455 = -t486 * t527 + t507 * t531;
t566 = t507 * t527;
t456 = t486 * t531 + t566;
t485 = t508 * t523 - t530 * t545;
t590 = t596 * t455 + t598 * t456 + t594 * t485;
t589 = t597 * t453 + t599 * t454 + t596 * t483;
t588 = t597 * t455 + t599 * t456 + t596 * t485;
t587 = t599 * t453 + t600 * t454 + t598 * t483;
t586 = t599 * t455 + t600 * t456 + t598 * t485;
t497 = t523 * t568 + t529 * t545;
t562 = t524 * t533;
t481 = -t497 * t527 - t531 * t562;
t553 = t527 * t562;
t482 = t497 * t531 - t553;
t565 = t524 * t529;
t496 = t523 * t565 - t546 * t568;
t585 = t596 * t481 + t598 * t482 + t594 * t496;
t584 = t597 * t481 + t599 * t482 + t596 * t496;
t583 = t599 * t481 + t600 * t482 + t598 * t496;
t528 = sin(qJ(3));
t532 = cos(qJ(3));
t487 = -t506 * t528 - t532 * t561;
t551 = t528 * t561;
t488 = t506 * t532 - t551;
t582 = Icges(4,5) * t488 + Icges(5,5) * t484 + Icges(4,6) * t487 - Icges(5,6) * t483 + t595 * t505;
t563 = t524 * t532;
t489 = -t508 * t528 + t530 * t563;
t552 = t528 * t564;
t490 = t508 * t532 + t552;
t581 = Icges(4,5) * t490 + Icges(5,5) * t486 + Icges(4,6) * t489 - Icges(5,6) * t485 + t595 * t507;
t503 = -t528 * t565 + t532 * t568;
t547 = t568 * t528;
t504 = t529 * t563 + t547;
t580 = Icges(4,5) * t504 + Icges(5,5) * t497 + Icges(4,6) * t503 - Icges(5,6) * t496 - t595 * t562;
t572 = pkin(3) * t532;
t571 = pkin(5) * t531;
t560 = rSges(7,1) * t454 + rSges(7,2) * t453 + pkin(5) * t567 + t593 * t483 + t484 * t571;
t559 = rSges(7,1) * t456 + rSges(7,2) * t455 + pkin(5) * t566 + t593 * t485 + t486 * t571;
t558 = rSges(7,1) * t482 + rSges(7,2) * t481 - pkin(5) * t553 + t593 * t496 + t497 * t571;
t479 = pkin(2) * t506 + pkin(9) * t505;
t480 = pkin(2) * t508 + pkin(9) * t507;
t555 = qJD(2) * t524;
t518 = t530 * t555;
t550 = t534 * t555;
t557 = t479 * t518 + t480 * t550;
t491 = qJD(3) * t507 + t518;
t556 = qJD(1) * (pkin(1) * t530 - pkin(8) * t561);
t519 = qJD(2) * t568 + qJD(1);
t492 = qJD(3) * t505 - t550;
t510 = -qJD(3) * t562 + t519;
t509 = (pkin(2) * t529 - pkin(9) * t533) * t524;
t511 = qJD(1) * (pkin(1) * t534 + pkin(8) * t564);
t543 = t519 * t480 - t509 * t518 + t511;
t436 = -pkin(3) * t551 + qJ(4) * t505 + t506 * t572;
t542 = -qJD(4) * t562 + t491 * t436 + t557;
t437 = pkin(3) * t552 + qJ(4) * t507 + t508 * t572;
t541 = qJD(4) * t505 + t510 * t437 + t543;
t540 = -t479 * t519 - t509 * t550 - t556;
t449 = pkin(4) * t484 + pkin(10) * t483;
t450 = pkin(4) * t486 + pkin(10) * t485;
t539 = t491 * t449 + (-t437 - t450) * t492 + t542;
t465 = pkin(3) * t547 + (-qJ(4) * t533 + t529 * t572) * t524;
t538 = qJD(4) * t507 + t492 * t465 + t540;
t461 = pkin(4) * t497 + pkin(10) * t496;
t537 = t510 * t450 + (-t461 - t465) * t491 + t541;
t536 = t492 * t461 + (-t436 - t449) * t510 + t538;
t515 = rSges(2,1) * t534 - rSges(2,2) * t530;
t514 = rSges(2,1) * t530 + rSges(2,2) * t534;
t498 = t568 * rSges(3,3) + (rSges(3,1) * t529 + rSges(3,2) * t533) * t524;
t495 = Icges(3,5) * t568 + (Icges(3,1) * t529 + Icges(3,4) * t533) * t524;
t494 = Icges(3,6) * t568 + (Icges(3,4) * t529 + Icges(3,2) * t533) * t524;
t493 = Icges(3,3) * t568 + (Icges(3,5) * t529 + Icges(3,6) * t533) * t524;
t478 = qJD(5) * t496 + t510;
t475 = rSges(3,1) * t508 - rSges(3,2) * t507 + rSges(3,3) * t564;
t474 = rSges(3,1) * t506 - rSges(3,2) * t505 - rSges(3,3) * t561;
t472 = Icges(3,1) * t508 - Icges(3,4) * t507 + Icges(3,5) * t564;
t471 = Icges(3,1) * t506 - Icges(3,4) * t505 - Icges(3,5) * t561;
t470 = Icges(3,4) * t508 - Icges(3,2) * t507 + Icges(3,6) * t564;
t469 = Icges(3,4) * t506 - Icges(3,2) * t505 - Icges(3,6) * t561;
t466 = rSges(4,1) * t504 + rSges(4,2) * t503 - rSges(4,3) * t562;
t464 = Icges(4,1) * t504 + Icges(4,4) * t503 - Icges(4,5) * t562;
t463 = Icges(4,4) * t504 + Icges(4,2) * t503 - Icges(4,6) * t562;
t460 = rSges(5,1) * t497 - rSges(5,2) * t496 - rSges(5,3) * t562;
t459 = Icges(5,1) * t497 - Icges(5,4) * t496 - Icges(5,5) * t562;
t458 = Icges(5,4) * t497 - Icges(5,2) * t496 - Icges(5,6) * t562;
t452 = qJD(5) * t483 + t492;
t451 = qJD(5) * t485 + t491;
t446 = rSges(4,1) * t490 + rSges(4,2) * t489 + rSges(4,3) * t507;
t445 = rSges(4,1) * t488 + rSges(4,2) * t487 + rSges(4,3) * t505;
t444 = Icges(4,1) * t490 + Icges(4,4) * t489 + Icges(4,5) * t507;
t443 = Icges(4,1) * t488 + Icges(4,4) * t487 + Icges(4,5) * t505;
t442 = Icges(4,4) * t490 + Icges(4,2) * t489 + Icges(4,6) * t507;
t441 = Icges(4,4) * t488 + Icges(4,2) * t487 + Icges(4,6) * t505;
t435 = rSges(5,1) * t486 - rSges(5,2) * t485 + rSges(5,3) * t507;
t434 = rSges(5,1) * t484 - rSges(5,2) * t483 + rSges(5,3) * t505;
t433 = Icges(5,1) * t486 - Icges(5,4) * t485 + Icges(5,5) * t507;
t432 = Icges(5,1) * t484 - Icges(5,4) * t483 + Icges(5,5) * t505;
t431 = Icges(5,4) * t486 - Icges(5,2) * t485 + Icges(5,6) * t507;
t430 = Icges(5,4) * t484 - Icges(5,2) * t483 + Icges(5,6) * t505;
t427 = t475 * t519 - t498 * t518 + t511;
t426 = -t474 * t519 - t498 * t550 - t556;
t425 = rSges(6,1) * t482 + rSges(6,2) * t481 + rSges(6,3) * t496;
t415 = (t474 * t530 + t475 * t534) * t555;
t412 = rSges(6,1) * t456 + rSges(6,2) * t455 + rSges(6,3) * t485;
t410 = rSges(6,1) * t454 + rSges(6,2) * t453 + rSges(6,3) * t483;
t394 = t446 * t510 - t466 * t491 + t543;
t393 = -t445 * t510 + t466 * t492 + t540;
t392 = t445 * t491 - t446 * t492 + t557;
t391 = t435 * t510 + (-t460 - t465) * t491 + t541;
t390 = t460 * t492 + (-t434 - t436) * t510 + t538;
t389 = t434 * t491 + (-t435 - t437) * t492 + t542;
t388 = t412 * t478 - t425 * t451 + t537;
t387 = -t410 * t478 + t425 * t452 + t536;
t386 = t410 * t451 - t412 * t452 + t539;
t385 = qJD(6) * t483 - t451 * t558 + t478 * t559 + t537;
t384 = qJD(6) * t485 + t452 * t558 - t478 * t560 + t536;
t383 = qJD(6) * t496 + t451 * t560 - t452 * t559 + t539;
t1 = ((t493 * t564 - t494 * t507 + t495 * t508) * t519 + (-(-t469 * t507 + t471 * t508) * t534 + (-t507 * t470 + t508 * t472 - t592) * t530) * t555) * t518 / 0.2e1 - ((-t493 * t561 - t494 * t505 + t495 * t506) * t519 + ((-t470 * t505 + t472 * t506) * t530 + (t505 * t469 - t506 * t471 + t592) * t534) * t555) * t550 / 0.2e1 + t519 * ((t568 * t468 + (t470 * t533 + t472 * t529) * t524) * t518 - (t568 * t467 + (t469 * t533 + t471 * t529) * t524) * t550 + (t568 * t493 + (t494 * t533 + t495 * t529) * t524) * t519) / 0.2e1 + m(6) * (t386 ^ 2 + t387 ^ 2 + t388 ^ 2) / 0.2e1 + m(7) * (t383 ^ 2 + t384 ^ 2 + t385 ^ 2) / 0.2e1 + m(4) * (t392 ^ 2 + t393 ^ 2 + t394 ^ 2) / 0.2e1 + m(5) * (t389 ^ 2 + t390 ^ 2 + t391 ^ 2) / 0.2e1 + m(3) * (t415 ^ 2 + t426 ^ 2 + t427 ^ 2) / 0.2e1 + ((t455 * t584 + t456 * t583 + t485 * t585) * t478 + (t455 * t589 + t456 * t587 + t485 * t591) * t452 + (t588 * t455 + t586 * t456 + t590 * t485) * t451) * t451 / 0.2e1 + ((t453 * t584 + t454 * t583 + t483 * t585) * t478 + (t589 * t453 + t587 * t454 + t591 * t483) * t452 + (t453 * t588 + t454 * t586 + t483 * t590) * t451) * t452 / 0.2e1 + ((t584 * t481 + t583 * t482 + t585 * t496) * t478 + (t481 * t589 + t482 * t587 + t496 * t591) * t452 + (t481 * t588 + t482 * t586 + t496 * t590) * t451) * t478 / 0.2e1 + ((-t485 * t458 + t459 * t486 + t463 * t489 + t464 * t490 + t507 * t580) * t510 + (-t430 * t485 + t432 * t486 + t441 * t489 + t443 * t490 + t507 * t582) * t492 + (-t431 * t485 + t433 * t486 + t442 * t489 + t444 * t490 + t581 * t507) * t491) * t491 / 0.2e1 + ((-t458 * t483 + t459 * t484 + t463 * t487 + t464 * t488 + t505 * t580) * t510 + (-t430 * t483 + t432 * t484 + t441 * t487 + t443 * t488 + t582 * t505) * t492 + (-t431 * t483 + t484 * t433 + t442 * t487 + t444 * t488 + t505 * t581) * t491) * t492 / 0.2e1 + ((-t458 * t496 + t459 * t497 + t463 * t503 + t464 * t504 - t580 * t562) * t510 + (-t430 * t496 + t432 * t497 + t441 * t503 + t443 * t504 - t562 * t582) * t492 + (-t431 * t496 + t433 * t497 + t442 * t503 + t444 * t504 - t562 * t581) * t491) * t510 / 0.2e1 + (Icges(2,3) + m(2) * (t514 ^ 2 + t515 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
