% Calculate kinetic energy for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR15_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR15_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR15_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:23:55
% EndTime: 2024-09-27 22:23:56
% DurationCPUTime: 0.45s
% Computational Cost: add. (3298->363), mult. (3705->547), div. (0->0), fcn. (3982->26), ass. (0->201)
t572 = sin(pkin(5));
t574 = cos(pkin(5));
t583 = cos(qJ(2));
t584 = cos(qJ(1));
t610 = t583 * t584;
t578 = sin(qJ(2));
t579 = sin(qJ(1));
t615 = t578 * t579;
t523 = t574 * t610 - t615;
t612 = t579 * t583;
t614 = t578 * t584;
t524 = t574 * t614 + t612;
t525 = -t574 * t612 - t614;
t526 = -t574 * t615 + t610;
t618 = t572 * t584;
t619 = t572 * t579;
t592 = (Icges(3,5) * t524 + Icges(3,6) * t523 - Icges(3,3) * t618) * t584 - (Icges(3,5) * t526 + Icges(3,6) * t525 + Icges(3,3) * t619) * t579;
t633 = t572 * t592;
t585 = pkin(8) + pkin(9);
t632 = -t579 / 0.2e1;
t631 = t584 / 0.2e1;
t629 = pkin(2) * t578;
t577 = sin(qJ(3));
t628 = pkin(3) * t577;
t571 = sin(pkin(6));
t627 = pkin(11) * t571;
t559 = t583 * pkin(2) + pkin(1);
t626 = -pkin(1) + t559;
t570 = qJ(2) + qJ(3);
t567 = qJ(4) + t570;
t556 = sin(t567);
t625 = t556 * t579;
t624 = t556 * t584;
t557 = cos(t567);
t573 = cos(pkin(6));
t623 = t557 * t573;
t622 = t557 * t579;
t621 = t557 * t584;
t620 = t571 * t574;
t575 = sin(qJ(5));
t617 = t575 * t579;
t616 = t575 * t584;
t580 = cos(qJ(5));
t613 = t579 * t580;
t611 = t580 * t584;
t576 = sin(qJ(4));
t581 = cos(qJ(4));
t531 = pkin(4) * t581 + t576 * t627 + pkin(3);
t535 = -pkin(4) * t576 + t581 * t627;
t582 = cos(qJ(3));
t497 = t531 * t582 + t535 * t577 + pkin(2);
t498 = -t531 * t577 + t535 * t582;
t558 = pkin(3) * t582 + pkin(2);
t569 = pkin(10) + t585;
t601 = t583 * t628;
t503 = -t569 * t572 + (t558 * t578 + t601) * t574;
t542 = pkin(11) * t573 + t569;
t609 = t542 * t572 + (-t497 * t578 + t498 * t583) * t574 + t503;
t566 = cos(t570);
t537 = pkin(3) * t566 + t559;
t608 = t497 * t583 + t498 * t578 + pkin(1) - t537;
t530 = -t572 * t585 + t574 * t629;
t593 = pkin(8) * t572 + t530;
t487 = t579 * t626 + t584 * t593;
t488 = -t579 * t593 + t584 * t626;
t603 = qJD(2) * t572;
t545 = t579 * t603;
t594 = t584 * t603;
t607 = t487 * t545 + t488 * t594;
t606 = t503 - t530;
t605 = t537 - t559;
t521 = qJD(3) * t619 + t545;
t604 = qJD(1) * (pkin(1) * t579 - pkin(8) * t618);
t548 = qJD(2) * t574 + qJD(1);
t602 = -qJD(2) - qJD(3);
t600 = t571 * t619;
t599 = t571 * t618;
t598 = t574 * t617;
t597 = t574 * t616;
t596 = t574 * t613;
t595 = t574 * t611;
t508 = qJD(4) * t619 + t521;
t534 = qJD(3) * t574 + t548;
t564 = pkin(5) - t570;
t563 = pkin(5) + t570;
t516 = qJD(4) * t574 + t534;
t444 = t579 * t605 + t584 * t606;
t445 = -t579 * t606 + t584 * t605;
t522 = t602 * t618;
t591 = t521 * t444 - t445 * t522 + t607;
t509 = (-qJD(4) + t602) * t618;
t515 = t572 * t629 + (-pkin(8) + t585) * t574;
t529 = qJD(1) * (pkin(1) * t584 + pkin(8) * t619);
t590 = t548 * t488 - t515 * t545 + t529;
t589 = -t487 * t548 - t515 * t594 - t604;
t482 = (t569 - t585) * t574 + (t601 + (-pkin(2) + t558) * t578) * t572;
t588 = t534 * t445 - t482 * t521 + t590;
t587 = -t444 * t534 + t522 * t482 + t589;
t565 = sin(t570);
t555 = -qJ(4) + t564;
t554 = qJ(4) + t563;
t552 = cos(t564);
t551 = cos(t563);
t550 = sin(t564);
t549 = sin(t563);
t547 = cos(t554);
t546 = sin(t555);
t541 = cos(t555) / 0.2e1;
t540 = sin(t554) / 0.2e1;
t539 = rSges(2,1) * t584 - rSges(2,2) * t579;
t538 = rSges(2,1) * t579 + rSges(2,2) * t584;
t533 = t552 + t551;
t532 = t549 - t550;
t528 = t552 / 0.2e1 - t551 / 0.2e1;
t527 = t549 / 0.2e1 + t550 / 0.2e1;
t520 = t541 - t547 / 0.2e1;
t519 = t541 + t547 / 0.2e1;
t518 = t540 - t546 / 0.2e1;
t517 = t540 + t546 / 0.2e1;
t514 = -t557 * t571 * t572 + t573 * t574;
t513 = rSges(3,3) * t574 + (rSges(3,1) * t578 + rSges(3,2) * t583) * t572;
t512 = Icges(3,5) * t574 + (Icges(3,1) * t578 + Icges(3,4) * t583) * t572;
t511 = Icges(3,6) * t574 + (Icges(3,4) * t578 + Icges(3,2) * t583) * t572;
t510 = Icges(3,3) * t574 + (Icges(3,5) * t578 + Icges(3,6) * t583) * t572;
t507 = t532 * t632 + t584 * t566;
t506 = t533 * t632 - t584 * t565;
t505 = t532 * t631 + t579 * t566;
t504 = t533 * t631 - t579 * t565;
t502 = -t518 * t579 + t621;
t501 = -t519 * t579 - t624;
t500 = t518 * t584 + t622;
t499 = t519 * t584 - t625;
t496 = t575 * t620 + (t556 * t580 + t575 * t623) * t572;
t495 = t580 * t620 + (-t556 * t575 + t580 * t623) * t572;
t494 = -t573 * t618 + (-t574 * t621 + t625) * t571;
t493 = t573 * t619 + (t574 * t622 + t624) * t571;
t492 = rSges(4,1) * t528 + rSges(4,2) * t527 + rSges(4,3) * t574;
t491 = Icges(4,1) * t528 + Icges(4,4) * t527 + Icges(4,5) * t574;
t490 = Icges(4,4) * t528 + Icges(4,2) * t527 + Icges(4,6) * t574;
t489 = Icges(4,5) * t528 + Icges(4,6) * t527 + Icges(4,3) * t574;
t486 = qJD(5) * t514 + t516;
t485 = rSges(5,1) * t520 + rSges(5,2) * t517 + rSges(5,3) * t574;
t484 = rSges(3,1) * t526 + rSges(3,2) * t525 + rSges(3,3) * t619;
t483 = rSges(3,1) * t524 + rSges(3,2) * t523 - rSges(3,3) * t618;
t481 = Icges(5,1) * t520 + Icges(5,4) * t517 + Icges(5,5) * t574;
t480 = Icges(5,4) * t520 + Icges(5,2) * t517 + Icges(5,6) * t574;
t479 = Icges(5,5) * t520 + Icges(5,6) * t517 + Icges(5,3) * t574;
t478 = Icges(3,1) * t526 + Icges(3,4) * t525 + Icges(3,5) * t619;
t477 = Icges(3,1) * t524 + Icges(3,4) * t523 - Icges(3,5) * t618;
t476 = Icges(3,4) * t526 + Icges(3,2) * t525 + Icges(3,6) * t619;
t475 = Icges(3,4) * t524 + Icges(3,2) * t523 - Icges(3,6) * t618;
t468 = rSges(4,1) * t507 + rSges(4,2) * t506 + rSges(4,3) * t619;
t467 = rSges(4,1) * t505 + rSges(4,2) * t504 - rSges(4,3) * t618;
t466 = (t573 * t595 - t617) * t557 + (-t573 * t613 - t597) * t556 - t580 * t599;
t465 = (t573 * t597 + t613) * t557 + (-t573 * t617 + t595) * t556 - t575 * t599;
t464 = (-t573 * t596 - t616) * t557 + (-t573 * t611 + t598) * t556 + t580 * t600;
t463 = (-t573 * t598 + t611) * t557 + (-t573 * t616 - t596) * t556 + t575 * t600;
t462 = qJD(5) * t494 + t509;
t461 = qJD(5) * t493 + t508;
t460 = Icges(4,1) * t507 + Icges(4,4) * t506 + Icges(4,5) * t619;
t459 = Icges(4,1) * t505 + Icges(4,4) * t504 - Icges(4,5) * t618;
t458 = Icges(4,4) * t507 + Icges(4,2) * t506 + Icges(4,6) * t619;
t457 = Icges(4,4) * t505 + Icges(4,2) * t504 - Icges(4,6) * t618;
t456 = Icges(4,5) * t507 + Icges(4,6) * t506 + Icges(4,3) * t619;
t455 = Icges(4,5) * t505 + Icges(4,6) * t504 - Icges(4,3) * t618;
t454 = rSges(5,1) * t502 + rSges(5,2) * t501 + rSges(5,3) * t619;
t453 = rSges(5,1) * t500 + rSges(5,2) * t499 - rSges(5,3) * t618;
t452 = Icges(5,1) * t502 + Icges(5,4) * t501 + Icges(5,5) * t619;
t451 = Icges(5,1) * t500 + Icges(5,4) * t499 - Icges(5,5) * t618;
t450 = Icges(5,4) * t502 + Icges(5,2) * t501 + Icges(5,6) * t619;
t449 = Icges(5,4) * t500 + Icges(5,2) * t499 - Icges(5,6) * t618;
t448 = Icges(5,5) * t502 + Icges(5,6) * t501 + Icges(5,3) * t619;
t447 = Icges(5,5) * t500 + Icges(5,6) * t499 - Icges(5,3) * t618;
t440 = t484 * t548 - t513 * t545 + t529;
t439 = -t483 * t548 - t513 * t594 - t604;
t438 = rSges(6,1) * t496 + rSges(6,2) * t495 + rSges(6,3) * t514;
t437 = Icges(6,1) * t496 + Icges(6,4) * t495 + Icges(6,5) * t514;
t436 = Icges(6,4) * t496 + Icges(6,2) * t495 + Icges(6,6) * t514;
t435 = Icges(6,5) * t496 + Icges(6,6) * t495 + Icges(6,3) * t514;
t434 = (t483 * t579 + t484 * t584) * t603;
t433 = (t542 - t569) * t574 + ((-t498 - t628) * t583 + (t497 - t558) * t578) * t572;
t432 = rSges(6,1) * t465 + rSges(6,2) * t466 + rSges(6,3) * t494;
t431 = rSges(6,1) * t463 + rSges(6,2) * t464 + rSges(6,3) * t493;
t430 = Icges(6,1) * t465 + Icges(6,4) * t466 + Icges(6,5) * t494;
t429 = Icges(6,1) * t463 + Icges(6,4) * t464 + Icges(6,5) * t493;
t428 = Icges(6,4) * t465 + Icges(6,2) * t466 + Icges(6,6) * t494;
t427 = Icges(6,4) * t463 + Icges(6,2) * t464 + Icges(6,6) * t493;
t426 = Icges(6,5) * t465 + Icges(6,6) * t466 + Icges(6,3) * t494;
t425 = Icges(6,5) * t463 + Icges(6,6) * t464 + Icges(6,3) * t493;
t424 = t468 * t534 - t492 * t521 + t590;
t423 = -t467 * t534 + t492 * t522 + t589;
t422 = t579 * t609 + t584 * t608;
t421 = t579 * t608 - t584 * t609;
t420 = t467 * t521 - t468 * t522 + t607;
t419 = t454 * t516 - t485 * t508 + t588;
t418 = -t453 * t516 + t485 * t509 + t587;
t417 = t453 * t508 - t454 * t509 + t591;
t416 = t422 * t516 + t431 * t486 - t433 * t508 - t438 * t461 + t588;
t415 = -t421 * t516 - t432 * t486 + t433 * t509 + t438 * t462 + t587;
t414 = t421 * t508 - t422 * t509 - t431 * t462 + t432 * t461 + t591;
t1 = m(3) * (t434 ^ 2 + t439 ^ 2 + t440 ^ 2) / 0.2e1 + ((t510 * t619 + t511 * t525 + t512 * t526) * t548 + (-(t475 * t525 + t477 * t526) * t584 + (t525 * t476 + t526 * t478 - t633) * t579) * t603) * t545 / 0.2e1 - ((-t510 * t618 + t511 * t523 + t512 * t524) * t548 + ((t476 * t523 + t478 * t524) * t579 + (-t523 * t475 - t524 * t477 + t633) * t584) * t603) * t594 / 0.2e1 + t548 * ((t574 * t510 + (t511 * t583 + t512 * t578) * t572) * t548 + (((t476 * t583 + t478 * t578) * t579 - (t475 * t583 + t477 * t578) * t584) * t572 - t592 * t574) * t603) / 0.2e1 + m(4) * (t420 ^ 2 + t423 ^ 2 + t424 ^ 2) / 0.2e1 + t521 * ((t456 * t619 + t506 * t458 + t507 * t460) * t521 + (t455 * t619 + t457 * t506 + t459 * t507) * t522 + (t489 * t619 + t490 * t506 + t491 * t507) * t534) / 0.2e1 + t522 * ((-t456 * t618 + t458 * t504 + t460 * t505) * t521 + (-t455 * t618 + t504 * t457 + t505 * t459) * t522 + (-t489 * t618 + t490 * t504 + t491 * t505) * t534) / 0.2e1 + t534 * ((t456 * t574 + t458 * t527 + t460 * t528) * t521 + (t455 * t574 + t457 * t527 + t459 * t528) * t522 + (t574 * t489 + t527 * t490 + t528 * t491) * t534) / 0.2e1 + m(5) * (t417 ^ 2 + t418 ^ 2 + t419 ^ 2) / 0.2e1 + t508 * ((t448 * t619 + t501 * t450 + t502 * t452) * t508 + (t447 * t619 + t449 * t501 + t451 * t502) * t509 + (t479 * t619 + t480 * t501 + t481 * t502) * t516) / 0.2e1 + t509 * ((-t448 * t618 + t450 * t499 + t452 * t500) * t508 + (-t447 * t618 + t499 * t449 + t500 * t451) * t509 + (-t479 * t618 + t480 * t499 + t481 * t500) * t516) / 0.2e1 + t516 * ((t448 * t574 + t450 * t517 + t452 * t520) * t508 + (t447 * t574 + t449 * t517 + t451 * t520) * t509 + (t574 * t479 + t517 * t480 + t520 * t481) * t516) / 0.2e1 + m(6) * (t414 ^ 2 + t415 ^ 2 + t416 ^ 2) / 0.2e1 + t461 * ((t493 * t425 + t464 * t427 + t463 * t429) * t461 + (t426 * t493 + t428 * t464 + t430 * t463) * t462 + (t435 * t493 + t436 * t464 + t437 * t463) * t486) / 0.2e1 + t462 * ((t425 * t494 + t427 * t466 + t429 * t465) * t461 + (t494 * t426 + t466 * t428 + t465 * t430) * t462 + (t435 * t494 + t436 * t466 + t437 * t465) * t486) / 0.2e1 + t486 * ((t425 * t514 + t427 * t495 + t429 * t496) * t461 + (t426 * t514 + t428 * t495 + t430 * t496) * t462 + (t514 * t435 + t495 * t436 + t496 * t437) * t486) / 0.2e1 + (m(2) * (t538 ^ 2 + t539 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T = t1;
