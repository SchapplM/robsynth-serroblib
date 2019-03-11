% Calculate kinetic energy for
% S6PRRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1]';
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
% Datum: 2019-03-08 23:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR8_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:48:11
% EndTime: 2019-03-08 23:48:13
% DurationCPUTime: 2.97s
% Computational Cost: add. (4578->339), mult. (12553->525), div. (0->0), fcn. (16090->14), ass. (0->157)
t617 = Icges(5,1) + Icges(6,2);
t616 = Icges(6,1) + Icges(5,3);
t615 = -Icges(5,4) - Icges(6,6);
t614 = -Icges(6,4) + Icges(5,5);
t613 = Icges(6,5) - Icges(5,6);
t612 = Icges(5,2) + Icges(6,3);
t563 = sin(pkin(12));
t565 = cos(pkin(12));
t570 = sin(qJ(2));
t566 = cos(pkin(6));
t572 = cos(qJ(2));
t590 = t566 * t572;
t554 = -t563 * t570 + t565 * t590;
t591 = t566 * t570;
t555 = t563 * t572 + t565 * t591;
t569 = sin(qJ(3));
t564 = sin(pkin(6));
t595 = sin(pkin(7));
t584 = t564 * t595;
t596 = cos(pkin(7));
t598 = cos(qJ(3));
t516 = t555 * t598 + (t554 * t596 - t565 * t584) * t569;
t585 = t564 * t596;
t541 = -t554 * t595 - t565 * t585;
t568 = sin(qJ(4));
t597 = cos(qJ(4));
t497 = t516 * t568 - t541 * t597;
t498 = t516 * t597 + t541 * t568;
t582 = t598 * t595;
t581 = t564 * t582;
t583 = t596 * t598;
t515 = -t554 * t583 + t555 * t569 + t565 * t581;
t611 = t497 * t612 + t498 * t615 + t515 * t613;
t556 = -t563 * t590 - t565 * t570;
t557 = -t563 * t591 + t565 * t572;
t518 = t557 * t598 + (t556 * t596 + t563 * t584) * t569;
t542 = -t556 * t595 + t563 * t585;
t499 = t518 * t568 - t542 * t597;
t500 = t518 * t597 + t542 * t568;
t517 = -t556 * t583 + t557 * t569 - t563 * t581;
t610 = t499 * t612 + t500 * t615 + t517 * t613;
t609 = t497 * t613 + t498 * t614 + t515 * t616;
t608 = t499 * t613 + t500 * t614 + t517 * t616;
t607 = t615 * t497 + t498 * t617 + t614 * t515;
t606 = t615 * t499 + t500 * t617 + t614 * t517;
t540 = t566 * t595 * t569 + (t569 * t572 * t596 + t570 * t598) * t564;
t553 = t566 * t596 - t572 * t584;
t520 = t540 * t568 - t553 * t597;
t521 = t540 * t597 + t553 * t568;
t592 = t564 * t570;
t539 = -t564 * t572 * t583 - t566 * t582 + t569 * t592;
t605 = t520 * t612 + t521 * t615 + t539 * t613;
t604 = t520 * t613 + t521 * t614 + t539 * t616;
t603 = t615 * t520 + t521 * t617 + t614 * t539;
t602 = qJD(2) ^ 2;
t594 = t563 * t564;
t593 = t564 * t565;
t589 = qJD(2) * t564;
t561 = t563 * t589;
t532 = qJD(3) * t542 + t561;
t562 = qJD(2) * t566;
t544 = qJD(3) * t553 + t562;
t491 = qJD(4) * t517 + t532;
t511 = qJD(4) * t539 + t544;
t587 = t565 * t589;
t522 = t555 * pkin(2) + pkin(9) * t541;
t523 = t557 * pkin(2) + pkin(9) * t542;
t586 = t522 * t561 + t523 * t587 + qJD(1);
t533 = qJD(3) * t541 - t587;
t543 = pkin(2) * t592 + pkin(9) * t553;
t580 = t523 * t562 - t543 * t561;
t492 = qJD(4) * t515 + t533;
t486 = pkin(3) * t516 + pkin(10) * t515;
t487 = pkin(3) * t518 + pkin(10) * t517;
t579 = t532 * t486 - t487 * t533 + t586;
t578 = (-t522 * t566 - t543 * t593) * qJD(2);
t459 = pkin(4) * t498 + qJ(5) * t497;
t577 = qJD(5) * t520 + t491 * t459 + t579;
t508 = pkin(3) * t540 + pkin(10) * t539;
t576 = t544 * t487 - t508 * t532 + t580;
t460 = pkin(4) * t500 + qJ(5) * t499;
t575 = qJD(5) * t497 + t511 * t460 + t576;
t574 = -t486 * t544 + t533 * t508 + t578;
t488 = pkin(4) * t521 + qJ(5) * t520;
t573 = qJD(5) * t499 + t492 * t488 + t574;
t571 = cos(qJ(6));
t567 = sin(qJ(6));
t550 = t566 * rSges(3,3) + (rSges(3,1) * t570 + rSges(3,2) * t572) * t564;
t549 = Icges(3,5) * t566 + (Icges(3,1) * t570 + Icges(3,4) * t572) * t564;
t548 = Icges(3,6) * t566 + (Icges(3,4) * t570 + Icges(3,2) * t572) * t564;
t547 = Icges(3,3) * t566 + (Icges(3,5) * t570 + Icges(3,6) * t572) * t564;
t531 = rSges(3,1) * t557 + rSges(3,2) * t556 + rSges(3,3) * t594;
t530 = rSges(3,1) * t555 + rSges(3,2) * t554 - rSges(3,3) * t593;
t529 = Icges(3,1) * t557 + Icges(3,4) * t556 + Icges(3,5) * t594;
t528 = Icges(3,1) * t555 + Icges(3,4) * t554 - Icges(3,5) * t593;
t527 = Icges(3,4) * t557 + Icges(3,2) * t556 + Icges(3,6) * t594;
t526 = Icges(3,4) * t555 + Icges(3,2) * t554 - Icges(3,6) * t593;
t525 = Icges(3,5) * t557 + Icges(3,6) * t556 + Icges(3,3) * t594;
t524 = Icges(3,5) * t555 + Icges(3,6) * t554 - Icges(3,3) * t593;
t507 = (-t530 * t566 - t550 * t593) * qJD(2);
t506 = (t531 * t566 - t550 * t594) * qJD(2);
t505 = rSges(4,1) * t540 - rSges(4,2) * t539 + rSges(4,3) * t553;
t504 = Icges(4,1) * t540 - Icges(4,4) * t539 + Icges(4,5) * t553;
t503 = Icges(4,4) * t540 - Icges(4,2) * t539 + Icges(4,6) * t553;
t502 = Icges(4,5) * t540 - Icges(4,6) * t539 + Icges(4,3) * t553;
t501 = pkin(5) * t539 + pkin(11) * t521;
t496 = t520 * t567 + t539 * t571;
t495 = t520 * t571 - t539 * t567;
t490 = qJD(1) + (t530 * t563 + t531 * t565) * t589;
t485 = qJD(6) * t521 + t511;
t483 = rSges(5,1) * t521 - rSges(5,2) * t520 + rSges(5,3) * t539;
t482 = rSges(6,1) * t539 - rSges(6,2) * t521 + rSges(6,3) * t520;
t481 = rSges(4,1) * t518 - rSges(4,2) * t517 + rSges(4,3) * t542;
t480 = rSges(4,1) * t516 - rSges(4,2) * t515 + rSges(4,3) * t541;
t473 = Icges(4,1) * t518 - Icges(4,4) * t517 + Icges(4,5) * t542;
t472 = Icges(4,1) * t516 - Icges(4,4) * t515 + Icges(4,5) * t541;
t471 = Icges(4,4) * t518 - Icges(4,2) * t517 + Icges(4,6) * t542;
t470 = Icges(4,4) * t516 - Icges(4,2) * t515 + Icges(4,6) * t541;
t469 = Icges(4,5) * t518 - Icges(4,6) * t517 + Icges(4,3) * t542;
t468 = Icges(4,5) * t516 - Icges(4,6) * t515 + Icges(4,3) * t541;
t467 = pkin(5) * t517 + pkin(11) * t500;
t466 = pkin(5) * t515 + pkin(11) * t498;
t465 = t499 * t567 + t517 * t571;
t464 = t499 * t571 - t517 * t567;
t463 = t497 * t567 + t515 * t571;
t462 = t497 * t571 - t515 * t567;
t458 = qJD(6) * t498 + t492;
t457 = qJD(6) * t500 + t491;
t455 = rSges(5,1) * t500 - rSges(5,2) * t499 + rSges(5,3) * t517;
t454 = rSges(5,1) * t498 - rSges(5,2) * t497 + rSges(5,3) * t515;
t453 = rSges(7,1) * t496 + rSges(7,2) * t495 + rSges(7,3) * t521;
t452 = rSges(6,1) * t517 - rSges(6,2) * t500 + rSges(6,3) * t499;
t451 = rSges(6,1) * t515 - rSges(6,2) * t498 + rSges(6,3) * t497;
t446 = Icges(7,1) * t496 + Icges(7,4) * t495 + Icges(7,5) * t521;
t441 = Icges(7,4) * t496 + Icges(7,2) * t495 + Icges(7,6) * t521;
t436 = Icges(7,5) * t496 + Icges(7,6) * t495 + Icges(7,3) * t521;
t433 = -t480 * t544 + t505 * t533 + t578;
t432 = t481 * t544 - t505 * t532 + t580;
t431 = rSges(7,1) * t465 + rSges(7,2) * t464 + rSges(7,3) * t500;
t430 = rSges(7,1) * t463 + rSges(7,2) * t462 + rSges(7,3) * t498;
t429 = Icges(7,1) * t465 + Icges(7,4) * t464 + Icges(7,5) * t500;
t428 = Icges(7,1) * t463 + Icges(7,4) * t462 + Icges(7,5) * t498;
t427 = Icges(7,4) * t465 + Icges(7,2) * t464 + Icges(7,6) * t500;
t426 = Icges(7,4) * t463 + Icges(7,2) * t462 + Icges(7,6) * t498;
t425 = Icges(7,5) * t465 + Icges(7,6) * t464 + Icges(7,3) * t500;
t424 = Icges(7,5) * t463 + Icges(7,6) * t462 + Icges(7,3) * t498;
t423 = t480 * t532 - t481 * t533 + t586;
t422 = -t454 * t511 + t483 * t492 + t574;
t421 = t455 * t511 - t483 * t491 + t576;
t420 = t454 * t491 - t455 * t492 + t579;
t419 = t482 * t492 + (-t451 - t459) * t511 + t573;
t418 = t452 * t511 + (-t482 - t488) * t491 + t575;
t417 = t451 * t491 + (-t452 - t460) * t492 + t577;
t416 = -t430 * t485 + t453 * t458 + t492 * t501 + (-t459 - t466) * t511 + t573;
t415 = t431 * t485 - t453 * t457 + t467 * t511 + (-t488 - t501) * t491 + t575;
t414 = t430 * t457 - t431 * t458 + t466 * t491 + (-t460 - t467) * t492 + t577;
t1 = -t602 * ((-t525 * t593 + t527 * t554 + t529 * t555) * t594 - (-t524 * t593 + t526 * t554 + t528 * t555) * t593 + (-t547 * t593 + t548 * t554 + t549 * t555) * t566) * t593 / 0.2e1 + m(7) * (t414 ^ 2 + t415 ^ 2 + t416 ^ 2) / 0.2e1 + m(5) * (t420 ^ 2 + t421 ^ 2 + t422 ^ 2) / 0.2e1 + m(6) * (t417 ^ 2 + t418 ^ 2 + t419 ^ 2) / 0.2e1 + m(4) * (t423 ^ 2 + t432 ^ 2 + t433 ^ 2) / 0.2e1 + m(3) * (t490 ^ 2 + t506 ^ 2 + t507 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t457 * ((t425 * t500 + t427 * t464 + t429 * t465) * t457 + (t424 * t500 + t426 * t464 + t428 * t465) * t458 + (t436 * t500 + t441 * t464 + t446 * t465) * t485) / 0.2e1 + t458 * ((t425 * t498 + t427 * t462 + t429 * t463) * t457 + (t424 * t498 + t426 * t462 + t428 * t463) * t458 + (t436 * t498 + t441 * t462 + t446 * t463) * t485) / 0.2e1 + t485 * ((t425 * t521 + t427 * t495 + t429 * t496) * t457 + (t424 * t521 + t426 * t495 + t428 * t496) * t458 + (t436 * t521 + t441 * t495 + t446 * t496) * t485) / 0.2e1 + t532 * ((t542 * t469 - t517 * t471 + t518 * t473) * t532 + (t468 * t542 - t470 * t517 + t472 * t518) * t533 + (t502 * t542 - t503 * t517 + t504 * t518) * t544) / 0.2e1 + t544 * ((t469 * t553 - t471 * t539 + t473 * t540) * t532 + (t468 * t553 - t470 * t539 + t472 * t540) * t533 + (t553 * t502 - t539 * t503 + t540 * t504) * t544) / 0.2e1 + t533 * ((t469 * t541 - t471 * t515 + t473 * t516) * t532 + (t541 * t468 - t515 * t470 + t516 * t472) * t533 + (t502 * t541 - t503 * t515 + t504 * t516) * t544) / 0.2e1 + ((t605 * t499 + t603 * t500 + t604 * t517) * t511 + (t611 * t499 + t607 * t500 + t609 * t517) * t492 + (t610 * t499 + t606 * t500 + t608 * t517) * t491) * t491 / 0.2e1 + ((t605 * t497 + t603 * t498 + t604 * t515) * t511 + (t611 * t497 + t607 * t498 + t609 * t515) * t492 + (t610 * t497 + t606 * t498 + t608 * t515) * t491) * t492 / 0.2e1 + ((t605 * t520 + t603 * t521 + t604 * t539) * t511 + (t611 * t520 + t607 * t521 + t609 * t539) * t492 + (t610 * t520 + t606 * t521 + t608 * t539) * t491) * t511 / 0.2e1 + (((t525 * t594 + t527 * t556 + t529 * t557) * t594 - (t524 * t594 + t526 * t556 + t528 * t557) * t593 + (t547 * t594 + t548 * t556 + t549 * t557) * t566) * t594 + t566 * (t566 ^ 2 * t547 + (((t527 * t572 + t529 * t570) * t563 - (t526 * t572 + t528 * t570) * t565) * t564 + (-t524 * t565 + t525 * t563 + t548 * t572 + t549 * t570) * t566) * t564)) * t602 / 0.2e1;
T  = t1;
