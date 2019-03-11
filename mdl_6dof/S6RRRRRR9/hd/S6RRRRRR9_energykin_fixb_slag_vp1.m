% Calculate kinetic energy for
% S6RRRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 05:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR9_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR9_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR9_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR9_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR9_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR9_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:16:42
% EndTime: 2019-03-10 05:16:45
% DurationCPUTime: 2.88s
% Computational Cost: add. (5681->392), mult. (14937->619), div. (0->0), fcn. (19256->16), ass. (0->181)
t576 = sin(pkin(6));
t577 = cos(pkin(6));
t584 = cos(qJ(2));
t585 = cos(qJ(1));
t609 = t584 * t585;
t581 = sin(qJ(2));
t582 = sin(qJ(1));
t612 = t581 * t582;
t558 = t577 * t609 - t612;
t610 = t582 * t584;
t611 = t581 * t585;
t559 = t577 * t611 + t610;
t560 = -t577 * t610 - t611;
t561 = -t577 * t612 + t609;
t613 = t576 * t585;
t614 = t576 * t582;
t600 = (Icges(3,5) * t559 + Icges(3,6) * t558 - Icges(3,3) * t613) * t585 - (Icges(3,5) * t561 + Icges(3,6) * t560 + Icges(3,3) * t614) * t582;
t626 = t576 * t600;
t624 = cos(qJ(3));
t623 = cos(qJ(4));
t583 = cos(qJ(5));
t622 = pkin(5) * t583;
t620 = cos(pkin(7));
t619 = sin(pkin(7));
t580 = sin(qJ(3));
t601 = t624 * t619;
t599 = t576 * t601;
t602 = t620 * t624;
t525 = -t558 * t602 + t559 * t580 + t585 * t599;
t578 = sin(qJ(5));
t618 = t525 * t578;
t527 = -t560 * t602 + t561 * t580 - t582 * t599;
t617 = t527 * t578;
t615 = t576 * t581;
t546 = -t576 * t584 * t602 - t577 * t601 + t580 * t615;
t616 = t546 * t578;
t604 = t576 * t620;
t596 = -t558 * t619 - t585 * t604;
t529 = t559 * pkin(2) + pkin(10) * t596;
t595 = -t560 * t619 + t582 * t604;
t530 = t561 * pkin(2) + pkin(10) * t595;
t606 = qJD(2) * t576;
t569 = t582 * t606;
t605 = t585 * t606;
t608 = t529 * t569 + t530 * t605;
t539 = qJD(3) * t595 + t569;
t607 = qJD(1) * (pkin(1) * t582 - pkin(9) * t613);
t570 = qJD(2) * t577 + qJD(1);
t499 = qJD(4) * t527 + t539;
t603 = t576 * t619;
t594 = t577 * t620 - t584 * t603;
t548 = qJD(3) * t594 + t570;
t528 = t561 * t624 + (t560 * t620 + t582 * t603) * t580;
t579 = sin(qJ(4));
t509 = t528 * t579 - t595 * t623;
t464 = qJD(5) * t509 + t499;
t516 = qJD(4) * t546 + t548;
t540 = qJD(3) * t596 - t605;
t547 = t577 * t619 * t580 + (t580 * t584 * t620 + t581 * t624) * t576;
t523 = t547 * t579 - t594 * t623;
t491 = qJD(5) * t523 + t516;
t526 = t559 * t624 + (t558 * t620 - t585 * t603) * t580;
t493 = pkin(3) * t526 + pkin(11) * t525;
t494 = pkin(3) * t528 + pkin(11) * t527;
t598 = t539 * t493 - t494 * t540 + t608;
t500 = qJD(4) * t525 + t540;
t549 = pkin(2) * t615 + pkin(10) * t594;
t562 = qJD(1) * (pkin(1) * t585 + pkin(9) * t614);
t597 = t570 * t530 - t549 * t569 + t562;
t507 = t526 * t579 - t596 * t623;
t465 = qJD(5) * t507 + t500;
t508 = t526 * t623 + t579 * t596;
t466 = pkin(4) * t508 + pkin(12) * t507;
t510 = t528 * t623 + t579 * t595;
t467 = pkin(4) * t510 + pkin(12) * t509;
t593 = t499 * t466 - t467 * t500 + t598;
t592 = -t529 * t570 - t549 * t605 - t607;
t515 = pkin(3) * t547 + pkin(11) * t546;
t591 = t548 * t494 - t515 * t539 + t597;
t590 = -t493 * t548 + t540 * t515 + t592;
t524 = t547 * t623 + t579 * t594;
t492 = pkin(4) * t524 + pkin(12) * t523;
t589 = t516 * t467 - t492 * t499 + t591;
t588 = -t466 * t516 + t500 * t492 + t590;
t575 = qJ(5) + qJ(6);
t574 = cos(t575);
t573 = sin(t575);
t567 = rSges(2,1) * t585 - rSges(2,2) * t582;
t566 = rSges(2,1) * t582 + rSges(2,2) * t585;
t555 = rSges(3,3) * t577 + (rSges(3,1) * t581 + rSges(3,2) * t584) * t576;
t554 = Icges(3,5) * t577 + (Icges(3,1) * t581 + Icges(3,4) * t584) * t576;
t553 = Icges(3,6) * t577 + (Icges(3,4) * t581 + Icges(3,2) * t584) * t576;
t552 = Icges(3,3) * t577 + (Icges(3,5) * t581 + Icges(3,6) * t584) * t576;
t538 = rSges(3,1) * t561 + rSges(3,2) * t560 + rSges(3,3) * t614;
t537 = rSges(3,1) * t559 + rSges(3,2) * t558 - rSges(3,3) * t613;
t536 = Icges(3,1) * t561 + Icges(3,4) * t560 + Icges(3,5) * t614;
t535 = Icges(3,1) * t559 + Icges(3,4) * t558 - Icges(3,5) * t613;
t534 = Icges(3,4) * t561 + Icges(3,2) * t560 + Icges(3,6) * t614;
t533 = Icges(3,4) * t559 + Icges(3,2) * t558 - Icges(3,6) * t613;
t514 = t547 * rSges(4,1) - t546 * rSges(4,2) + rSges(4,3) * t594;
t513 = Icges(4,1) * t547 - Icges(4,4) * t546 + Icges(4,5) * t594;
t512 = Icges(4,4) * t547 - Icges(4,2) * t546 + Icges(4,6) * t594;
t511 = Icges(4,5) * t547 - Icges(4,6) * t546 + Icges(4,3) * t594;
t504 = t524 * t583 + t616;
t503 = -t524 * t578 + t546 * t583;
t502 = t538 * t570 - t555 * t569 + t562;
t501 = -t537 * t570 - t555 * t605 - t607;
t498 = t524 * t574 + t546 * t573;
t497 = -t524 * t573 + t546 * t574;
t496 = (t537 * t582 + t538 * t585) * t606;
t489 = t528 * rSges(4,1) - t527 * rSges(4,2) + rSges(4,3) * t595;
t488 = t526 * rSges(4,1) - t525 * rSges(4,2) + rSges(4,3) * t596;
t487 = Icges(4,1) * t528 - Icges(4,4) * t527 + Icges(4,5) * t595;
t486 = Icges(4,1) * t526 - Icges(4,4) * t525 + Icges(4,5) * t596;
t485 = Icges(4,4) * t528 - Icges(4,2) * t527 + Icges(4,6) * t595;
t484 = Icges(4,4) * t526 - Icges(4,2) * t525 + Icges(4,6) * t596;
t483 = Icges(4,5) * t528 - Icges(4,6) * t527 + Icges(4,3) * t595;
t482 = Icges(4,5) * t526 - Icges(4,6) * t525 + Icges(4,3) * t596;
t481 = rSges(5,1) * t524 - rSges(5,2) * t523 + rSges(5,3) * t546;
t480 = Icges(5,1) * t524 - Icges(5,4) * t523 + Icges(5,5) * t546;
t479 = Icges(5,4) * t524 - Icges(5,2) * t523 + Icges(5,6) * t546;
t478 = Icges(5,5) * t524 - Icges(5,6) * t523 + Icges(5,3) * t546;
t477 = t510 * t583 + t617;
t476 = -t510 * t578 + t527 * t583;
t475 = t508 * t583 + t618;
t474 = -t508 * t578 + t525 * t583;
t472 = t510 * t574 + t527 * t573;
t471 = -t510 * t573 + t527 * t574;
t470 = t508 * t574 + t525 * t573;
t469 = -t508 * t573 + t525 * t574;
t468 = qJD(6) * t523 + t491;
t462 = rSges(5,1) * t510 - rSges(5,2) * t509 + rSges(5,3) * t527;
t461 = rSges(5,1) * t508 - rSges(5,2) * t507 + rSges(5,3) * t525;
t460 = Icges(5,1) * t510 - Icges(5,4) * t509 + Icges(5,5) * t527;
t459 = Icges(5,1) * t508 - Icges(5,4) * t507 + Icges(5,5) * t525;
t458 = Icges(5,4) * t510 - Icges(5,2) * t509 + Icges(5,6) * t527;
t457 = Icges(5,4) * t508 - Icges(5,2) * t507 + Icges(5,6) * t525;
t456 = Icges(5,5) * t510 - Icges(5,6) * t509 + Icges(5,3) * t527;
t455 = Icges(5,5) * t508 - Icges(5,6) * t507 + Icges(5,3) * t525;
t453 = rSges(6,1) * t504 + rSges(6,2) * t503 + rSges(6,3) * t523;
t452 = Icges(6,1) * t504 + Icges(6,4) * t503 + Icges(6,5) * t523;
t451 = Icges(6,4) * t504 + Icges(6,2) * t503 + Icges(6,6) * t523;
t450 = Icges(6,5) * t504 + Icges(6,6) * t503 + Icges(6,3) * t523;
t449 = rSges(7,1) * t498 + rSges(7,2) * t497 + rSges(7,3) * t523;
t448 = Icges(7,1) * t498 + Icges(7,4) * t497 + Icges(7,5) * t523;
t447 = Icges(7,4) * t498 + Icges(7,2) * t497 + Icges(7,6) * t523;
t446 = Icges(7,5) * t498 + Icges(7,6) * t497 + Icges(7,3) * t523;
t445 = pkin(5) * t616 + pkin(13) * t523 + t524 * t622;
t444 = qJD(6) * t507 + t465;
t443 = qJD(6) * t509 + t464;
t441 = rSges(6,1) * t477 + rSges(6,2) * t476 + rSges(6,3) * t509;
t440 = rSges(6,1) * t475 + rSges(6,2) * t474 + rSges(6,3) * t507;
t439 = t489 * t548 - t514 * t539 + t597;
t438 = -t488 * t548 + t514 * t540 + t592;
t437 = Icges(6,1) * t477 + Icges(6,4) * t476 + Icges(6,5) * t509;
t436 = Icges(6,1) * t475 + Icges(6,4) * t474 + Icges(6,5) * t507;
t435 = Icges(6,4) * t477 + Icges(6,2) * t476 + Icges(6,6) * t509;
t434 = Icges(6,4) * t475 + Icges(6,2) * t474 + Icges(6,6) * t507;
t433 = Icges(6,5) * t477 + Icges(6,6) * t476 + Icges(6,3) * t509;
t432 = Icges(6,5) * t475 + Icges(6,6) * t474 + Icges(6,3) * t507;
t431 = rSges(7,1) * t472 + rSges(7,2) * t471 + rSges(7,3) * t509;
t430 = rSges(7,1) * t470 + rSges(7,2) * t469 + rSges(7,3) * t507;
t429 = Icges(7,1) * t472 + Icges(7,4) * t471 + Icges(7,5) * t509;
t428 = Icges(7,1) * t470 + Icges(7,4) * t469 + Icges(7,5) * t507;
t427 = Icges(7,4) * t472 + Icges(7,2) * t471 + Icges(7,6) * t509;
t426 = Icges(7,4) * t470 + Icges(7,2) * t469 + Icges(7,6) * t507;
t425 = Icges(7,5) * t472 + Icges(7,6) * t471 + Icges(7,3) * t509;
t424 = Icges(7,5) * t470 + Icges(7,6) * t469 + Icges(7,3) * t507;
t423 = pkin(5) * t617 + pkin(13) * t509 + t510 * t622;
t422 = pkin(5) * t618 + pkin(13) * t507 + t508 * t622;
t421 = t488 * t539 - t489 * t540 + t608;
t420 = t462 * t516 - t481 * t499 + t591;
t419 = -t461 * t516 + t481 * t500 + t590;
t418 = t461 * t499 - t462 * t500 + t598;
t417 = t441 * t491 - t453 * t464 + t589;
t416 = -t440 * t491 + t453 * t465 + t588;
t415 = t440 * t464 - t441 * t465 + t593;
t414 = t423 * t491 + t431 * t468 - t443 * t449 - t445 * t464 + t589;
t413 = -t422 * t491 - t430 * t468 + t444 * t449 + t445 * t465 + t588;
t412 = t422 * t464 - t423 * t465 + t430 * t443 - t431 * t444 + t593;
t1 = m(4) * (t421 ^ 2 + t438 ^ 2 + t439 ^ 2) / 0.2e1 + m(3) * (t496 ^ 2 + t501 ^ 2 + t502 ^ 2) / 0.2e1 + t444 * ((t425 * t507 + t427 * t469 + t429 * t470) * t443 + (t507 * t424 + t469 * t426 + t470 * t428) * t444 + (t446 * t507 + t447 * t469 + t448 * t470) * t468) / 0.2e1 + t468 * ((t425 * t523 + t427 * t497 + t429 * t498) * t443 + (t424 * t523 + t426 * t497 + t428 * t498) * t444 + (t523 * t446 + t497 * t447 + t498 * t448) * t468) / 0.2e1 + t443 * ((t509 * t425 + t471 * t427 + t472 * t429) * t443 + (t424 * t509 + t426 * t471 + t428 * t472) * t444 + (t446 * t509 + t447 * t471 + t448 * t472) * t468) / 0.2e1 + t491 * ((t433 * t523 + t435 * t503 + t437 * t504) * t464 + (t432 * t523 + t434 * t503 + t436 * t504) * t465 + (t523 * t450 + t503 * t451 + t504 * t452) * t491) / 0.2e1 + t465 * ((t433 * t507 + t435 * t474 + t437 * t475) * t464 + (t507 * t432 + t474 * t434 + t475 * t436) * t465 + (t450 * t507 + t451 * t474 + t452 * t475) * t491) / 0.2e1 + t464 * ((t509 * t433 + t476 * t435 + t477 * t437) * t464 + (t432 * t509 + t434 * t476 + t436 * t477) * t465 + (t450 * t509 + t451 * t476 + t452 * t477) * t491) / 0.2e1 + t500 * ((t456 * t525 - t458 * t507 + t460 * t508) * t499 + (t525 * t455 - t507 * t457 + t508 * t459) * t500 + (t478 * t525 - t479 * t507 + t480 * t508) * t516) / 0.2e1 + t516 * ((t456 * t546 - t458 * t523 + t460 * t524) * t499 + (t455 * t546 - t457 * t523 + t459 * t524) * t500 + (t546 * t478 - t523 * t479 + t524 * t480) * t516) / 0.2e1 + t499 * ((t527 * t456 - t509 * t458 + t510 * t460) * t499 + (t455 * t527 - t457 * t509 + t459 * t510) * t500 + (t478 * t527 - t479 * t509 + t480 * t510) * t516) / 0.2e1 + t539 * ((t595 * t483 - t527 * t485 + t528 * t487) * t539 + (t482 * t595 - t527 * t484 + t528 * t486) * t540 + (t511 * t595 - t527 * t512 + t528 * t513) * t548) / 0.2e1 + t540 * ((t483 * t596 - t525 * t485 + t526 * t487) * t539 + (t596 * t482 - t525 * t484 + t526 * t486) * t540 + (t511 * t596 - t525 * t512 + t526 * t513) * t548) / 0.2e1 + t548 * ((t483 * t594 - t546 * t485 + t547 * t487) * t539 + (t482 * t594 - t546 * t484 + t547 * t486) * t540 + (t594 * t511 - t546 * t512 + t547 * t513) * t548) / 0.2e1 + t570 * ((t577 * t552 + (t553 * t584 + t554 * t581) * t576) * t570 + (((t534 * t584 + t536 * t581) * t582 - (t533 * t584 + t535 * t581) * t585) * t576 - t600 * t577) * t606) / 0.2e1 + m(7) * (t412 ^ 2 + t413 ^ 2 + t414 ^ 2) / 0.2e1 + m(6) * (t415 ^ 2 + t416 ^ 2 + t417 ^ 2) / 0.2e1 + m(5) * (t418 ^ 2 + t419 ^ 2 + t420 ^ 2) / 0.2e1 - ((-t552 * t613 + t553 * t558 + t554 * t559) * t570 + ((t534 * t558 + t536 * t559) * t582 + (-t558 * t533 - t559 * t535 + t626) * t585) * t606) * t605 / 0.2e1 + ((t552 * t614 + t553 * t560 + t554 * t561) * t570 + (-(t533 * t560 + t535 * t561) * t585 + (t560 * t534 + t561 * t536 - t626) * t582) * t606) * t569 / 0.2e1 + (m(2) * (t566 ^ 2 + t567 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
