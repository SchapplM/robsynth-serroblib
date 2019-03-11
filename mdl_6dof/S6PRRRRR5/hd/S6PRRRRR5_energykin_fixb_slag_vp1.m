% Calculate kinetic energy for
% S6PRRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:03:53
% EndTime: 2019-03-09 01:03:55
% DurationCPUTime: 2.38s
% Computational Cost: add. (5612->385), mult. (14885->613), div. (0->0), fcn. (19222->16), ass. (0->176)
t615 = qJD(2) ^ 2;
t614 = cos(qJ(3));
t613 = cos(qJ(4));
t578 = cos(qJ(5));
t612 = pkin(5) * t578;
t610 = cos(pkin(7));
t609 = sin(pkin(7));
t570 = sin(pkin(13));
t572 = cos(pkin(13));
t577 = sin(qJ(2));
t573 = cos(pkin(6));
t579 = cos(qJ(2));
t601 = t573 * t579;
t557 = -t570 * t577 + t572 * t601;
t602 = t573 * t577;
t558 = t570 * t579 + t572 * t602;
t576 = sin(qJ(3));
t571 = sin(pkin(6));
t593 = t614 * t609;
t592 = t571 * t593;
t594 = t610 * t614;
t521 = -t557 * t594 + t558 * t576 + t572 * t592;
t574 = sin(qJ(5));
t608 = t521 * t574;
t559 = -t570 * t601 - t572 * t577;
t560 = -t570 * t602 + t572 * t579;
t523 = -t559 * t594 + t560 * t576 - t570 * t592;
t607 = t523 * t574;
t604 = t571 * t577;
t545 = -t571 * t579 * t594 - t573 * t593 + t576 * t604;
t606 = t545 * t574;
t605 = t570 * t571;
t603 = t572 * t571;
t600 = qJD(2) * t571;
t564 = t570 * t600;
t596 = t571 * t610;
t587 = -t559 * t609 + t570 * t596;
t538 = qJD(3) * t587 + t564;
t566 = qJD(2) * t573;
t595 = t571 * t609;
t586 = t573 * t610 - t579 * t595;
t548 = qJD(3) * t586 + t566;
t496 = qJD(4) * t523 + t538;
t517 = qJD(4) * t545 + t548;
t598 = t572 * t600;
t588 = -t557 * t609 - t572 * t596;
t528 = t558 * pkin(2) + pkin(9) * t588;
t529 = t560 * pkin(2) + pkin(9) * t587;
t597 = t528 * t564 + t529 * t598 + qJD(1);
t524 = t560 * t614 + (t559 * t610 + t570 * t595) * t576;
t575 = sin(qJ(4));
t506 = t524 * t575 - t587 * t613;
t463 = qJD(5) * t506 + t496;
t546 = t573 * t609 * t576 + (t576 * t579 * t610 + t577 * t614) * t571;
t526 = t546 * t575 - t586 * t613;
t490 = qJD(5) * t526 + t517;
t539 = qJD(3) * t588 - t598;
t547 = pkin(2) * t604 + pkin(9) * t586;
t591 = t529 * t566 - t547 * t564;
t497 = qJD(4) * t521 + t539;
t522 = t558 * t614 + (t557 * t610 - t572 * t595) * t576;
t504 = t522 * t575 - t588 * t613;
t464 = qJD(5) * t504 + t497;
t491 = pkin(3) * t522 + pkin(10) * t521;
t492 = pkin(3) * t524 + pkin(10) * t523;
t590 = t538 * t491 - t492 * t539 + t597;
t589 = (-t528 * t573 - t547 * t603) * qJD(2);
t514 = pkin(3) * t546 + pkin(10) * t545;
t585 = t548 * t492 - t514 * t538 + t591;
t505 = t522 * t613 + t575 * t588;
t465 = t505 * pkin(4) + t504 * pkin(11);
t507 = t524 * t613 + t575 * t587;
t466 = t507 * pkin(4) + t506 * pkin(11);
t584 = t496 * t465 - t466 * t497 + t590;
t583 = -t491 * t548 + t539 * t514 + t589;
t527 = t546 * t613 + t575 * t586;
t493 = t527 * pkin(4) + t526 * pkin(11);
t582 = t517 * t466 - t493 * t496 + t585;
t581 = -t465 * t517 + t497 * t493 + t583;
t569 = qJ(5) + qJ(6);
t568 = cos(t569);
t567 = sin(t569);
t554 = rSges(3,3) * t573 + (rSges(3,1) * t577 + rSges(3,2) * t579) * t571;
t553 = Icges(3,5) * t573 + (Icges(3,1) * t577 + Icges(3,4) * t579) * t571;
t552 = Icges(3,6) * t573 + (Icges(3,4) * t577 + Icges(3,2) * t579) * t571;
t551 = Icges(3,3) * t573 + (Icges(3,5) * t577 + Icges(3,6) * t579) * t571;
t537 = rSges(3,1) * t560 + rSges(3,2) * t559 + rSges(3,3) * t605;
t536 = rSges(3,1) * t558 + rSges(3,2) * t557 - rSges(3,3) * t603;
t535 = Icges(3,1) * t560 + Icges(3,4) * t559 + Icges(3,5) * t605;
t534 = Icges(3,1) * t558 + Icges(3,4) * t557 - Icges(3,5) * t603;
t533 = Icges(3,4) * t560 + Icges(3,2) * t559 + Icges(3,6) * t605;
t532 = Icges(3,4) * t558 + Icges(3,2) * t557 - Icges(3,6) * t603;
t531 = Icges(3,5) * t560 + Icges(3,6) * t559 + Icges(3,3) * t605;
t530 = Icges(3,5) * t558 + Icges(3,6) * t557 - Icges(3,3) * t603;
t513 = (-t536 * t573 - t554 * t603) * qJD(2);
t512 = (t537 * t573 - t554 * t605) * qJD(2);
t511 = t546 * rSges(4,1) - t545 * rSges(4,2) + rSges(4,3) * t586;
t510 = Icges(4,1) * t546 - Icges(4,4) * t545 + Icges(4,5) * t586;
t509 = Icges(4,4) * t546 - Icges(4,2) * t545 + Icges(4,6) * t586;
t508 = Icges(4,5) * t546 - Icges(4,6) * t545 + Icges(4,3) * t586;
t503 = t527 * t578 + t606;
t502 = -t527 * t574 + t545 * t578;
t499 = t527 * t568 + t545 * t567;
t498 = -t527 * t567 + t545 * t568;
t495 = qJD(1) + (t536 * t570 + t537 * t572) * t600;
t488 = rSges(5,1) * t527 - rSges(5,2) * t526 + rSges(5,3) * t545;
t487 = t524 * rSges(4,1) - t523 * rSges(4,2) + rSges(4,3) * t587;
t486 = t522 * rSges(4,1) - t521 * rSges(4,2) + rSges(4,3) * t588;
t485 = Icges(5,1) * t527 - Icges(5,4) * t526 + Icges(5,5) * t545;
t484 = Icges(5,4) * t527 - Icges(5,2) * t526 + Icges(5,6) * t545;
t483 = Icges(5,5) * t527 - Icges(5,6) * t526 + Icges(5,3) * t545;
t482 = Icges(4,1) * t524 - Icges(4,4) * t523 + Icges(4,5) * t587;
t481 = Icges(4,1) * t522 - Icges(4,4) * t521 + Icges(4,5) * t588;
t480 = Icges(4,4) * t524 - Icges(4,2) * t523 + Icges(4,6) * t587;
t479 = Icges(4,4) * t522 - Icges(4,2) * t521 + Icges(4,6) * t588;
t478 = Icges(4,5) * t524 - Icges(4,6) * t523 + Icges(4,3) * t587;
t477 = Icges(4,5) * t522 - Icges(4,6) * t521 + Icges(4,3) * t588;
t476 = t507 * t578 + t607;
t475 = -t507 * t574 + t523 * t578;
t474 = t505 * t578 + t608;
t473 = -t505 * t574 + t521 * t578;
t471 = t507 * t568 + t523 * t567;
t470 = -t507 * t567 + t523 * t568;
t469 = t505 * t568 + t521 * t567;
t468 = -t505 * t567 + t521 * t568;
t467 = qJD(6) * t526 + t490;
t461 = rSges(6,1) * t503 + rSges(6,2) * t502 + rSges(6,3) * t526;
t460 = rSges(5,1) * t507 - rSges(5,2) * t506 + rSges(5,3) * t523;
t459 = rSges(5,1) * t505 - rSges(5,2) * t504 + rSges(5,3) * t521;
t458 = Icges(5,1) * t507 - Icges(5,4) * t506 + Icges(5,5) * t523;
t457 = Icges(5,1) * t505 - Icges(5,4) * t504 + Icges(5,5) * t521;
t456 = Icges(6,1) * t503 + Icges(6,4) * t502 + Icges(6,5) * t526;
t455 = Icges(5,4) * t507 - Icges(5,2) * t506 + Icges(5,6) * t523;
t454 = Icges(5,4) * t505 - Icges(5,2) * t504 + Icges(5,6) * t521;
t453 = Icges(6,4) * t503 + Icges(6,2) * t502 + Icges(6,6) * t526;
t452 = Icges(5,5) * t507 - Icges(5,6) * t506 + Icges(5,3) * t523;
t451 = Icges(5,5) * t505 - Icges(5,6) * t504 + Icges(5,3) * t521;
t450 = Icges(6,5) * t503 + Icges(6,6) * t502 + Icges(6,3) * t526;
t448 = rSges(7,1) * t499 + rSges(7,2) * t498 + rSges(7,3) * t526;
t447 = Icges(7,1) * t499 + Icges(7,4) * t498 + Icges(7,5) * t526;
t446 = Icges(7,4) * t499 + Icges(7,2) * t498 + Icges(7,6) * t526;
t445 = Icges(7,5) * t499 + Icges(7,6) * t498 + Icges(7,3) * t526;
t444 = pkin(5) * t606 + pkin(12) * t526 + t527 * t612;
t443 = qJD(6) * t504 + t464;
t442 = qJD(6) * t506 + t463;
t440 = -t486 * t548 + t511 * t539 + t589;
t439 = t487 * t548 - t511 * t538 + t591;
t438 = rSges(6,1) * t476 + rSges(6,2) * t475 + rSges(6,3) * t506;
t437 = rSges(6,1) * t474 + rSges(6,2) * t473 + rSges(6,3) * t504;
t436 = Icges(6,1) * t476 + Icges(6,4) * t475 + Icges(6,5) * t506;
t435 = Icges(6,1) * t474 + Icges(6,4) * t473 + Icges(6,5) * t504;
t434 = Icges(6,4) * t476 + Icges(6,2) * t475 + Icges(6,6) * t506;
t433 = Icges(6,4) * t474 + Icges(6,2) * t473 + Icges(6,6) * t504;
t432 = Icges(6,5) * t476 + Icges(6,6) * t475 + Icges(6,3) * t506;
t431 = Icges(6,5) * t474 + Icges(6,6) * t473 + Icges(6,3) * t504;
t430 = rSges(7,1) * t471 + rSges(7,2) * t470 + rSges(7,3) * t506;
t429 = rSges(7,1) * t469 + rSges(7,2) * t468 + rSges(7,3) * t504;
t428 = Icges(7,1) * t471 + Icges(7,4) * t470 + Icges(7,5) * t506;
t427 = Icges(7,1) * t469 + Icges(7,4) * t468 + Icges(7,5) * t504;
t426 = Icges(7,4) * t471 + Icges(7,2) * t470 + Icges(7,6) * t506;
t425 = Icges(7,4) * t469 + Icges(7,2) * t468 + Icges(7,6) * t504;
t424 = Icges(7,5) * t471 + Icges(7,6) * t470 + Icges(7,3) * t506;
t423 = Icges(7,5) * t469 + Icges(7,6) * t468 + Icges(7,3) * t504;
t422 = pkin(5) * t607 + pkin(12) * t506 + t507 * t612;
t421 = pkin(5) * t608 + pkin(12) * t504 + t505 * t612;
t420 = t486 * t538 - t487 * t539 + t597;
t419 = -t459 * t517 + t488 * t497 + t583;
t418 = t460 * t517 - t488 * t496 + t585;
t417 = t459 * t496 - t460 * t497 + t590;
t416 = -t437 * t490 + t461 * t464 + t581;
t415 = t438 * t490 - t461 * t463 + t582;
t414 = t437 * t463 - t438 * t464 + t584;
t413 = -t421 * t490 - t429 * t467 + t443 * t448 + t444 * t464 + t581;
t412 = t422 * t490 + t430 * t467 - t442 * t448 - t444 * t463 + t582;
t411 = t421 * t463 - t422 * t464 + t429 * t442 - t430 * t443 + t584;
t1 = -t615 * ((-t531 * t603 + t533 * t557 + t535 * t558) * t605 - (-t530 * t603 + t532 * t557 + t534 * t558) * t603 + (-t551 * t603 + t552 * t557 + t553 * t558) * t573) * t603 / 0.2e1 + t443 * ((t424 * t504 + t426 * t468 + t428 * t469) * t442 + (t504 * t423 + t468 * t425 + t469 * t427) * t443 + (t445 * t504 + t446 * t468 + t447 * t469) * t467) / 0.2e1 + t467 * ((t424 * t526 + t426 * t498 + t428 * t499) * t442 + (t423 * t526 + t425 * t498 + t427 * t499) * t443 + (t526 * t445 + t498 * t446 + t499 * t447) * t467) / 0.2e1 + t464 * ((t432 * t504 + t434 * t473 + t436 * t474) * t463 + (t504 * t431 + t473 * t433 + t474 * t435) * t464 + (t450 * t504 + t453 * t473 + t456 * t474) * t490) / 0.2e1 + t442 * ((t506 * t424 + t470 * t426 + t471 * t428) * t442 + (t423 * t506 + t425 * t470 + t427 * t471) * t443 + (t445 * t506 + t446 * t470 + t447 * t471) * t467) / 0.2e1 + t463 * ((t506 * t432 + t475 * t434 + t476 * t436) * t463 + (t431 * t506 + t433 * t475 + t435 * t476) * t464 + (t450 * t506 + t453 * t475 + t456 * t476) * t490) / 0.2e1 + t490 * ((t432 * t526 + t434 * t502 + t436 * t503) * t463 + (t431 * t526 + t433 * t502 + t435 * t503) * t464 + (t526 * t450 + t502 * t453 + t503 * t456) * t490) / 0.2e1 + t496 * ((t523 * t452 - t506 * t455 + t507 * t458) * t496 + (t451 * t523 - t454 * t506 + t457 * t507) * t497 + (t483 * t523 - t484 * t506 + t485 * t507) * t517) / 0.2e1 + t497 * ((t452 * t521 - t455 * t504 + t458 * t505) * t496 + (t521 * t451 - t504 * t454 + t505 * t457) * t497 + (t483 * t521 - t484 * t504 + t485 * t505) * t517) / 0.2e1 + t517 * ((t452 * t545 - t455 * t526 + t458 * t527) * t496 + (t451 * t545 - t454 * t526 + t457 * t527) * t497 + (t545 * t483 - t526 * t484 + t527 * t485) * t517) / 0.2e1 + t539 * ((t478 * t588 - t521 * t480 + t522 * t482) * t538 + (t588 * t477 - t521 * t479 + t522 * t481) * t539 + (t508 * t588 - t521 * t509 + t522 * t510) * t548) / 0.2e1 + t538 * ((t587 * t478 - t523 * t480 + t524 * t482) * t538 + (t477 * t587 - t523 * t479 + t524 * t481) * t539 + (t508 * t587 - t523 * t509 + t524 * t510) * t548) / 0.2e1 + t548 * ((t478 * t586 - t545 * t480 + t546 * t482) * t538 + (t477 * t586 - t545 * t479 + t546 * t481) * t539 + (t586 * t508 - t545 * t509 + t546 * t510) * t548) / 0.2e1 + m(7) * (t411 ^ 2 + t412 ^ 2 + t413 ^ 2) / 0.2e1 + m(6) * (t414 ^ 2 + t415 ^ 2 + t416 ^ 2) / 0.2e1 + m(5) * (t417 ^ 2 + t418 ^ 2 + t419 ^ 2) / 0.2e1 + m(4) * (t420 ^ 2 + t439 ^ 2 + t440 ^ 2) / 0.2e1 + m(3) * (t495 ^ 2 + t512 ^ 2 + t513 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + (((t531 * t605 + t533 * t559 + t535 * t560) * t605 - (t530 * t605 + t532 * t559 + t534 * t560) * t603 + (t551 * t605 + t552 * t559 + t553 * t560) * t573) * t605 + t573 * (t573 ^ 2 * t551 + (((t533 * t579 + t535 * t577) * t570 - (t532 * t579 + t534 * t577) * t572) * t571 + (-t530 * t572 + t531 * t570 + t552 * t579 + t553 * t577) * t573) * t571)) * t615 / 0.2e1;
T  = t1;
