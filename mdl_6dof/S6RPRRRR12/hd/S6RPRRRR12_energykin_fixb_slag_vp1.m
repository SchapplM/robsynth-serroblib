% Calculate kinetic energy for
% S6RPRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR12_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_energykin_fixb_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR12_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR12_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:50:24
% EndTime: 2019-03-09 07:50:26
% DurationCPUTime: 2.31s
% Computational Cost: add. (9212->349), mult. (25933->555), div. (0->0), fcn. (34233->18), ass. (0->162)
t603 = cos(pkin(6));
t638 = t603 ^ 2;
t636 = cos(qJ(4));
t635 = cos(qJ(5));
t634 = cos(pkin(8));
t633 = sin(pkin(8));
t599 = sin(pkin(7));
t632 = t599 * t603;
t600 = sin(pkin(6));
t608 = sin(qJ(1));
t631 = t600 * t608;
t611 = cos(qJ(1));
t630 = t600 * t611;
t601 = cos(pkin(14));
t602 = cos(pkin(7));
t629 = t601 * t602;
t628 = t603 * t608;
t627 = t603 * t611;
t598 = sin(pkin(14));
t589 = t598 * t627 + t601 * t608;
t607 = sin(qJ(3));
t610 = cos(qJ(3));
t588 = -t598 * t608 + t601 * t627;
t619 = t588 * t602 - t599 * t630;
t569 = -t589 * t607 + t610 * t619;
t582 = -t588 * t599 - t602 * t630;
t554 = -t569 * t633 + t582 * t634;
t578 = qJD(3) * t582;
t544 = qJD(4) * t554 + t578;
t591 = -t598 * t628 + t601 * t611;
t590 = -t598 * t611 - t601 * t628;
t618 = t590 * t602 + t599 * t631;
t571 = -t591 * t607 + t610 * t618;
t583 = -t590 * t599 + t602 * t631;
t555 = -t571 * t633 + t583 * t634;
t579 = qJD(3) * t583;
t545 = qJD(4) * t555 + t579;
t626 = qJD(2) * t600;
t587 = -t599 * t600 * t601 + t602 * t603;
t585 = qJD(3) * t587 + qJD(1);
t570 = t589 * t610 + t607 * t619;
t606 = sin(qJ(4));
t624 = t636 * t633;
t625 = t634 * t636;
t530 = -t569 * t625 + t570 * t606 - t582 * t624;
t510 = qJD(5) * t530 + t544;
t572 = t591 * t610 + t607 * t618;
t532 = -t571 * t625 + t572 * t606 - t583 * t624;
t511 = qJD(5) * t532 + t545;
t580 = t610 * t632 + (-t598 * t607 + t610 * t629) * t600;
t568 = -t580 * t633 + t587 * t634;
t561 = qJD(4) * t568 + t585;
t581 = t607 * t632 + (t598 * t610 + t607 * t629) * t600;
t548 = -t580 * t625 + t581 * t606 - t587 * t624;
t524 = qJD(5) * t548 + t561;
t623 = -t611 * t626 + qJD(1) * (pkin(1) * t611 + qJ(2) * t631);
t593 = pkin(1) * t608 - qJ(2) * t630;
t596 = t608 * t626;
t622 = t596 + (-pkin(2) * t589 - pkin(10) * t582 - t593) * qJD(1);
t621 = qJD(1) * (pkin(2) * t591 + pkin(10) * t583) + t623;
t534 = t570 * pkin(3) + pkin(11) * t554;
t535 = t572 * pkin(3) + pkin(11) * t555;
t597 = qJD(2) * t603;
t620 = t534 * t579 - t535 * t578 + t597;
t556 = t581 * pkin(3) + pkin(11) * t568;
t617 = -t534 * t585 + t556 * t578 + t622;
t531 = t570 * t636 + (t569 * t634 + t582 * t633) * t606;
t508 = pkin(4) * t531 + pkin(12) * t530;
t533 = t572 * t636 + (t571 * t634 + t583 * t633) * t606;
t509 = pkin(4) * t533 + pkin(12) * t532;
t616 = t545 * t508 - t509 * t544 + t620;
t615 = t585 * t535 - t556 * t579 + t621;
t549 = t581 * t636 + (t580 * t634 + t587 * t633) * t606;
t522 = pkin(4) * t549 + pkin(12) * t548;
t614 = -t508 * t561 + t544 * t522 + t617;
t613 = t561 * t509 - t522 * t545 + t615;
t609 = cos(qJ(6));
t605 = sin(qJ(5));
t604 = sin(qJ(6));
t595 = rSges(2,1) * t611 - rSges(2,2) * t608;
t594 = rSges(2,1) * t608 + rSges(2,2) * t611;
t563 = qJD(1) * (rSges(3,1) * t591 + rSges(3,2) * t590 + rSges(3,3) * t631) + t623;
t562 = t596 + (-rSges(3,1) * t589 - rSges(3,2) * t588 + rSges(3,3) * t630 - t593) * qJD(1);
t560 = rSges(4,1) * t581 + rSges(4,2) * t580 + rSges(4,3) * t587;
t559 = Icges(4,1) * t581 + Icges(4,4) * t580 + Icges(4,5) * t587;
t558 = Icges(4,4) * t581 + Icges(4,2) * t580 + Icges(4,6) * t587;
t557 = Icges(4,5) * t581 + Icges(4,6) * t580 + Icges(4,3) * t587;
t543 = rSges(4,1) * t572 + rSges(4,2) * t571 + rSges(4,3) * t583;
t542 = rSges(4,1) * t570 + rSges(4,2) * t569 + rSges(4,3) * t582;
t541 = Icges(4,1) * t572 + Icges(4,4) * t571 + Icges(4,5) * t583;
t540 = Icges(4,1) * t570 + Icges(4,4) * t569 + Icges(4,5) * t582;
t539 = Icges(4,4) * t572 + Icges(4,2) * t571 + Icges(4,6) * t583;
t538 = Icges(4,4) * t570 + Icges(4,2) * t569 + Icges(4,6) * t582;
t537 = Icges(4,5) * t572 + Icges(4,6) * t571 + Icges(4,3) * t583;
t536 = Icges(4,5) * t570 + Icges(4,6) * t569 + Icges(4,3) * t582;
t527 = t549 * t635 + t568 * t605;
t526 = t549 * t605 - t568 * t635;
t521 = rSges(5,1) * t549 - rSges(5,2) * t548 + rSges(5,3) * t568;
t520 = t533 * t635 + t555 * t605;
t519 = t533 * t605 - t555 * t635;
t518 = t531 * t635 + t554 * t605;
t517 = t531 * t605 - t554 * t635;
t516 = Icges(5,1) * t549 - Icges(5,4) * t548 + Icges(5,5) * t568;
t515 = Icges(5,4) * t549 - Icges(5,2) * t548 + Icges(5,6) * t568;
t514 = Icges(5,5) * t549 - Icges(5,6) * t548 + Icges(5,3) * t568;
t513 = t527 * t609 + t548 * t604;
t512 = -t527 * t604 + t548 * t609;
t507 = pkin(5) * t527 + pkin(13) * t526;
t505 = qJD(6) * t526 + t524;
t504 = t543 * t585 - t560 * t579 + t621;
t503 = -t542 * t585 + t560 * t578 + t622;
t502 = t597 + (t542 * t583 - t543 * t582) * qJD(3);
t500 = rSges(5,1) * t533 - rSges(5,2) * t532 + rSges(5,3) * t555;
t499 = rSges(5,1) * t531 - rSges(5,2) * t530 + rSges(5,3) * t554;
t498 = Icges(5,1) * t533 - Icges(5,4) * t532 + Icges(5,5) * t555;
t497 = Icges(5,1) * t531 - Icges(5,4) * t530 + Icges(5,5) * t554;
t496 = Icges(5,4) * t533 - Icges(5,2) * t532 + Icges(5,6) * t555;
t495 = Icges(5,4) * t531 - Icges(5,2) * t530 + Icges(5,6) * t554;
t494 = Icges(5,5) * t533 - Icges(5,6) * t532 + Icges(5,3) * t555;
t493 = Icges(5,5) * t531 - Icges(5,6) * t530 + Icges(5,3) * t554;
t492 = t520 * t609 + t532 * t604;
t491 = -t520 * t604 + t532 * t609;
t490 = t518 * t609 + t530 * t604;
t489 = -t518 * t604 + t530 * t609;
t488 = rSges(6,1) * t527 - rSges(6,2) * t526 + rSges(6,3) * t548;
t487 = Icges(6,1) * t527 - Icges(6,4) * t526 + Icges(6,5) * t548;
t486 = Icges(6,4) * t527 - Icges(6,2) * t526 + Icges(6,6) * t548;
t485 = Icges(6,5) * t527 - Icges(6,6) * t526 + Icges(6,3) * t548;
t483 = pkin(5) * t520 + pkin(13) * t519;
t482 = pkin(5) * t518 + pkin(13) * t517;
t481 = qJD(6) * t519 + t511;
t480 = qJD(6) * t517 + t510;
t479 = rSges(6,1) * t520 - rSges(6,2) * t519 + rSges(6,3) * t532;
t478 = rSges(6,1) * t518 - rSges(6,2) * t517 + rSges(6,3) * t530;
t477 = Icges(6,1) * t520 - Icges(6,4) * t519 + Icges(6,5) * t532;
t476 = Icges(6,1) * t518 - Icges(6,4) * t517 + Icges(6,5) * t530;
t475 = Icges(6,4) * t520 - Icges(6,2) * t519 + Icges(6,6) * t532;
t474 = Icges(6,4) * t518 - Icges(6,2) * t517 + Icges(6,6) * t530;
t473 = Icges(6,5) * t520 - Icges(6,6) * t519 + Icges(6,3) * t532;
t472 = Icges(6,5) * t518 - Icges(6,6) * t517 + Icges(6,3) * t530;
t471 = rSges(7,1) * t513 + rSges(7,2) * t512 + rSges(7,3) * t526;
t470 = Icges(7,1) * t513 + Icges(7,4) * t512 + Icges(7,5) * t526;
t469 = Icges(7,4) * t513 + Icges(7,2) * t512 + Icges(7,6) * t526;
t468 = Icges(7,5) * t513 + Icges(7,6) * t512 + Icges(7,3) * t526;
t467 = rSges(7,1) * t492 + rSges(7,2) * t491 + rSges(7,3) * t519;
t466 = rSges(7,1) * t490 + rSges(7,2) * t489 + rSges(7,3) * t517;
t465 = Icges(7,1) * t492 + Icges(7,4) * t491 + Icges(7,5) * t519;
t464 = Icges(7,1) * t490 + Icges(7,4) * t489 + Icges(7,5) * t517;
t463 = Icges(7,4) * t492 + Icges(7,2) * t491 + Icges(7,6) * t519;
t462 = Icges(7,4) * t490 + Icges(7,2) * t489 + Icges(7,6) * t517;
t461 = Icges(7,5) * t492 + Icges(7,6) * t491 + Icges(7,3) * t519;
t460 = Icges(7,5) * t490 + Icges(7,6) * t489 + Icges(7,3) * t517;
t459 = t500 * t561 - t521 * t545 + t615;
t458 = -t499 * t561 + t521 * t544 + t617;
t457 = t499 * t545 - t500 * t544 + t620;
t456 = t479 * t524 - t488 * t511 + t613;
t455 = -t478 * t524 + t488 * t510 + t614;
t454 = t478 * t511 - t479 * t510 + t616;
t453 = t467 * t505 - t471 * t481 + t483 * t524 - t507 * t511 + t613;
t452 = -t466 * t505 + t471 * t480 - t482 * t524 + t507 * t510 + t614;
t451 = t466 * t481 - t467 * t480 + t482 * t511 - t483 * t510 + t616;
t1 = m(7) * (t451 ^ 2 + t452 ^ 2 + t453 ^ 2) / 0.2e1 + m(6) * (t454 ^ 2 + t455 ^ 2 + t456 ^ 2) / 0.2e1 + m(5) * (t457 ^ 2 + t458 ^ 2 + t459 ^ 2) / 0.2e1 + m(4) * (t502 ^ 2 + t503 ^ 2 + t504 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 * t638 + t562 ^ 2 + t563 ^ 2) / 0.2e1 + t545 * ((t555 * t494 - t532 * t496 + t533 * t498) * t545 + (t493 * t555 - t495 * t532 + t497 * t533) * t544 + (t514 * t555 - t515 * t532 + t516 * t533) * t561) / 0.2e1 + t544 * ((t494 * t554 - t496 * t530 + t498 * t531) * t545 + (t554 * t493 - t530 * t495 + t531 * t497) * t544 + (t514 * t554 - t515 * t530 + t516 * t531) * t561) / 0.2e1 + t561 * ((t494 * t568 - t496 * t548 + t498 * t549) * t545 + (t493 * t568 - t495 * t548 + t497 * t549) * t544 + (t568 * t514 - t548 * t515 + t549 * t516) * t561) / 0.2e1 + ((t557 * t582 + t558 * t569 + t559 * t570) * t585 + ((t537 * t582 + t539 * t569 + t541 * t570) * t583 + (t536 * t582 + t538 * t569 + t540 * t570) * t582) * qJD(3)) * t578 / 0.2e1 + t585 * ((t587 * t557 + t580 * t558 + t581 * t559) * t585 + ((t537 * t587 + t539 * t580 + t541 * t581) * t583 + (t536 * t587 + t538 * t580 + t540 * t581) * t582) * qJD(3)) / 0.2e1 + t481 * ((t519 * t461 + t491 * t463 + t492 * t465) * t481 + (t460 * t519 + t462 * t491 + t464 * t492) * t480 + (t468 * t519 + t469 * t491 + t470 * t492) * t505) / 0.2e1 + t505 * ((t461 * t526 + t463 * t512 + t465 * t513) * t481 + (t460 * t526 + t462 * t512 + t464 * t513) * t480 + (t526 * t468 + t512 * t469 + t513 * t470) * t505) / 0.2e1 + t480 * ((t461 * t517 + t463 * t489 + t465 * t490) * t481 + (t517 * t460 + t489 * t462 + t490 * t464) * t480 + (t468 * t517 + t469 * t489 + t470 * t490) * t505) / 0.2e1 + t511 * ((t532 * t473 - t519 * t475 + t520 * t477) * t511 + (t472 * t532 - t474 * t519 + t476 * t520) * t510 + (t485 * t532 - t486 * t519 + t487 * t520) * t524) / 0.2e1 + t510 * ((t473 * t530 - t475 * t517 + t477 * t518) * t511 + (t530 * t472 - t517 * t474 + t518 * t476) * t510 + (t485 * t530 - t486 * t517 + t487 * t518) * t524) / 0.2e1 + t524 * ((t473 * t548 - t526 * t475 + t527 * t477) * t511 + (t472 * t548 - t474 * t526 + t476 * t527) * t510 + (t548 * t485 - t526 * t486 + t527 * t487) * t524) / 0.2e1 + ((t557 * t583 + t558 * t571 + t559 * t572) * t585 + ((t537 * t583 + t539 * t571 + t541 * t572) * t583 + (t536 * t583 + t538 * t571 + t540 * t572) * t582) * qJD(3)) * t579 / 0.2e1 + (Icges(3,3) * t638 + ((Icges(3,2) * t601 ^ 2 + (Icges(3,1) * t598 + 0.2e1 * Icges(3,4) * t601) * t598) * t600 + 0.2e1 * t603 * (Icges(3,5) * t598 + Icges(3,6) * t601)) * t600 + m(2) * (t594 ^ 2 + t595 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
