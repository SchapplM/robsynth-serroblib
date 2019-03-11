% Calculate kinetic energy for
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:38:58
% EndTime: 2019-03-08 23:39:01
% DurationCPUTime: 3.38s
% Computational Cost: add. (5360->380), mult. (14225->583), div. (0->0), fcn. (18346->16), ass. (0->171)
t623 = Icges(5,2) + Icges(6,3);
t572 = sin(pkin(12));
t575 = cos(pkin(12));
t580 = sin(qJ(2));
t576 = cos(pkin(6));
t581 = cos(qJ(2));
t603 = t576 * t581;
t558 = -t572 * t580 + t575 * t603;
t604 = t576 * t580;
t559 = t572 * t581 + t575 * t604;
t579 = sin(qJ(3));
t573 = sin(pkin(6));
t611 = sin(pkin(7));
t596 = t573 * t611;
t612 = cos(pkin(7));
t615 = cos(qJ(3));
t523 = t559 * t615 + (t558 * t612 - t575 * t596) * t579;
t578 = sin(qJ(4));
t597 = t573 * t612;
t588 = -t558 * t611 - t575 * t597;
t614 = cos(qJ(4));
t506 = t523 * t614 + t578 * t588;
t594 = t615 * t611;
t593 = t573 * t594;
t595 = t612 * t615;
t522 = -t558 * t595 + t559 * t579 + t575 * t593;
t571 = sin(pkin(13));
t574 = cos(pkin(13));
t473 = -t506 * t571 + t522 * t574;
t610 = t522 * t571;
t474 = t506 * t574 + t610;
t505 = t523 * t578 - t588 * t614;
t622 = -Icges(5,4) * t506 + Icges(6,5) * t474 - Icges(5,6) * t522 + Icges(6,6) * t473 + t505 * t623;
t560 = -t572 * t603 - t575 * t580;
t561 = -t572 * t604 + t575 * t581;
t525 = t561 * t615 + (t560 * t612 + t572 * t596) * t579;
t587 = -t560 * t611 + t572 * t597;
t508 = t525 * t614 + t578 * t587;
t524 = -t560 * t595 + t561 * t579 - t572 * t593;
t475 = -t508 * t571 + t524 * t574;
t609 = t524 * t571;
t476 = t508 * t574 + t609;
t507 = t525 * t578 - t587 * t614;
t621 = -Icges(5,4) * t508 + Icges(6,5) * t476 - Icges(5,6) * t524 + Icges(6,6) * t475 + t507 * t623;
t547 = t576 * t611 * t579 + (t579 * t581 * t612 + t580 * t615) * t573;
t586 = t576 * t612 - t581 * t596;
t528 = t547 * t614 + t578 * t586;
t606 = t573 * t580;
t546 = -t573 * t581 * t595 - t576 * t594 + t579 * t606;
t503 = -t528 * t571 + t546 * t574;
t608 = t546 * t571;
t504 = t528 * t574 + t608;
t527 = t547 * t578 - t586 * t614;
t620 = -Icges(5,4) * t528 + Icges(6,5) * t504 - Icges(5,6) * t546 + Icges(6,6) * t503 + t527 * t623;
t619 = qJD(2) ^ 2;
t613 = pkin(5) * t574;
t607 = t572 * t573;
t605 = t575 * t573;
t601 = qJD(2) * t573;
t565 = t572 * t601;
t539 = qJD(3) * t587 + t565;
t567 = qJD(2) * t576;
t549 = qJD(3) * t586 + t567;
t497 = qJD(4) * t524 + t539;
t518 = qJD(4) * t546 + t549;
t599 = t575 * t601;
t529 = t559 * pkin(2) + pkin(9) * t588;
t530 = t561 * pkin(2) + pkin(9) * t587;
t598 = t529 * t565 + t530 * t599 + qJD(1);
t540 = qJD(3) * t588 - t599;
t548 = pkin(2) * t606 + pkin(9) * t586;
t592 = t530 * t567 - t548 * t565;
t498 = qJD(4) * t522 + t540;
t492 = pkin(3) * t523 + pkin(10) * t522;
t493 = pkin(3) * t525 + pkin(10) * t524;
t591 = t539 * t492 - t493 * t540 + t598;
t590 = (-t529 * t576 - t548 * t605) * qJD(2);
t467 = pkin(4) * t506 + qJ(5) * t505;
t589 = qJD(5) * t527 + t497 * t467 + t591;
t515 = pkin(3) * t547 + pkin(10) * t546;
t585 = t549 * t493 - t515 * t539 + t592;
t468 = pkin(4) * t508 + qJ(5) * t507;
t584 = qJD(5) * t505 + t518 * t468 + t585;
t583 = -t492 * t549 + t540 * t515 + t590;
t494 = pkin(4) * t528 + qJ(5) * t527;
t582 = qJD(5) * t507 + t498 * t494 + t583;
t570 = pkin(13) + qJ(6);
t569 = cos(t570);
t568 = sin(t570);
t555 = t576 * rSges(3,3) + (rSges(3,1) * t580 + rSges(3,2) * t581) * t573;
t554 = Icges(3,5) * t576 + (Icges(3,1) * t580 + Icges(3,4) * t581) * t573;
t553 = Icges(3,6) * t576 + (Icges(3,4) * t580 + Icges(3,2) * t581) * t573;
t552 = Icges(3,3) * t576 + (Icges(3,5) * t580 + Icges(3,6) * t581) * t573;
t538 = rSges(3,1) * t561 + rSges(3,2) * t560 + rSges(3,3) * t607;
t537 = rSges(3,1) * t559 + rSges(3,2) * t558 - rSges(3,3) * t605;
t536 = Icges(3,1) * t561 + Icges(3,4) * t560 + Icges(3,5) * t607;
t535 = Icges(3,1) * t559 + Icges(3,4) * t558 - Icges(3,5) * t605;
t534 = Icges(3,4) * t561 + Icges(3,2) * t560 + Icges(3,6) * t607;
t533 = Icges(3,4) * t559 + Icges(3,2) * t558 - Icges(3,6) * t605;
t532 = Icges(3,5) * t561 + Icges(3,6) * t560 + Icges(3,3) * t607;
t531 = Icges(3,5) * t559 + Icges(3,6) * t558 - Icges(3,3) * t605;
t514 = (-t537 * t576 - t555 * t605) * qJD(2);
t513 = (t538 * t576 - t555 * t607) * qJD(2);
t512 = t547 * rSges(4,1) - t546 * rSges(4,2) + rSges(4,3) * t586;
t511 = Icges(4,1) * t547 - Icges(4,4) * t546 + Icges(4,5) * t586;
t510 = Icges(4,4) * t547 - Icges(4,2) * t546 + Icges(4,6) * t586;
t509 = Icges(4,5) * t547 - Icges(4,6) * t546 + Icges(4,3) * t586;
t500 = t528 * t569 + t546 * t568;
t499 = -t528 * t568 + t546 * t569;
t496 = qJD(1) + (t537 * t572 + t538 * t575) * t601;
t491 = qJD(6) * t527 + t518;
t489 = rSges(5,1) * t528 - rSges(5,2) * t527 + rSges(5,3) * t546;
t488 = t525 * rSges(4,1) - t524 * rSges(4,2) + rSges(4,3) * t587;
t487 = t523 * rSges(4,1) - t522 * rSges(4,2) + rSges(4,3) * t588;
t486 = Icges(5,1) * t528 - Icges(5,4) * t527 + Icges(5,5) * t546;
t484 = Icges(5,5) * t528 - Icges(5,6) * t527 + Icges(5,3) * t546;
t483 = Icges(4,1) * t525 - Icges(4,4) * t524 + Icges(4,5) * t587;
t482 = Icges(4,1) * t523 - Icges(4,4) * t522 + Icges(4,5) * t588;
t481 = Icges(4,4) * t525 - Icges(4,2) * t524 + Icges(4,6) * t587;
t480 = Icges(4,4) * t523 - Icges(4,2) * t522 + Icges(4,6) * t588;
t479 = Icges(4,5) * t525 - Icges(4,6) * t524 + Icges(4,3) * t587;
t478 = Icges(4,5) * t523 - Icges(4,6) * t522 + Icges(4,3) * t588;
t472 = t508 * t569 + t524 * t568;
t471 = -t508 * t568 + t524 * t569;
t470 = t506 * t569 + t522 * t568;
t469 = -t506 * t568 + t522 * t569;
t466 = qJD(6) * t505 + t498;
t465 = qJD(6) * t507 + t497;
t463 = rSges(5,1) * t508 - rSges(5,2) * t507 + rSges(5,3) * t524;
t462 = rSges(5,1) * t506 - rSges(5,2) * t505 + rSges(5,3) * t522;
t461 = Icges(5,1) * t508 - Icges(5,4) * t507 + Icges(5,5) * t524;
t460 = Icges(5,1) * t506 - Icges(5,4) * t505 + Icges(5,5) * t522;
t457 = Icges(5,5) * t508 - Icges(5,6) * t507 + Icges(5,3) * t524;
t456 = Icges(5,5) * t506 - Icges(5,6) * t505 + Icges(5,3) * t522;
t454 = rSges(6,1) * t504 + rSges(6,2) * t503 + rSges(6,3) * t527;
t453 = Icges(6,1) * t504 + Icges(6,4) * t503 + Icges(6,5) * t527;
t452 = Icges(6,4) * t504 + Icges(6,2) * t503 + Icges(6,6) * t527;
t450 = rSges(7,1) * t500 + rSges(7,2) * t499 + rSges(7,3) * t527;
t449 = Icges(7,1) * t500 + Icges(7,4) * t499 + Icges(7,5) * t527;
t448 = Icges(7,4) * t500 + Icges(7,2) * t499 + Icges(7,6) * t527;
t447 = Icges(7,5) * t500 + Icges(7,6) * t499 + Icges(7,3) * t527;
t446 = pkin(5) * t608 + pkin(11) * t527 + t528 * t613;
t444 = -t487 * t549 + t512 * t540 + t590;
t443 = t488 * t549 - t512 * t539 + t592;
t442 = rSges(6,1) * t476 + rSges(6,2) * t475 + rSges(6,3) * t507;
t441 = rSges(6,1) * t474 + rSges(6,2) * t473 + rSges(6,3) * t505;
t440 = Icges(6,1) * t476 + Icges(6,4) * t475 + Icges(6,5) * t507;
t439 = Icges(6,1) * t474 + Icges(6,4) * t473 + Icges(6,5) * t505;
t438 = Icges(6,4) * t476 + Icges(6,2) * t475 + Icges(6,6) * t507;
t437 = Icges(6,4) * t474 + Icges(6,2) * t473 + Icges(6,6) * t505;
t434 = rSges(7,1) * t472 + rSges(7,2) * t471 + rSges(7,3) * t507;
t433 = rSges(7,1) * t470 + rSges(7,2) * t469 + rSges(7,3) * t505;
t432 = Icges(7,1) * t472 + Icges(7,4) * t471 + Icges(7,5) * t507;
t431 = Icges(7,1) * t470 + Icges(7,4) * t469 + Icges(7,5) * t505;
t430 = Icges(7,4) * t472 + Icges(7,2) * t471 + Icges(7,6) * t507;
t429 = Icges(7,4) * t470 + Icges(7,2) * t469 + Icges(7,6) * t505;
t428 = Icges(7,5) * t472 + Icges(7,6) * t471 + Icges(7,3) * t507;
t427 = Icges(7,5) * t470 + Icges(7,6) * t469 + Icges(7,3) * t505;
t426 = pkin(5) * t609 + pkin(11) * t507 + t508 * t613;
t425 = pkin(5) * t610 + pkin(11) * t505 + t506 * t613;
t424 = t487 * t539 - t488 * t540 + t598;
t423 = -t462 * t518 + t489 * t498 + t583;
t422 = t463 * t518 - t489 * t497 + t585;
t421 = t462 * t497 - t463 * t498 + t591;
t420 = t454 * t498 + (-t441 - t467) * t518 + t582;
t419 = t442 * t518 + (-t454 - t494) * t497 + t584;
t418 = t441 * t497 + (-t442 - t468) * t498 + t589;
t417 = -t433 * t491 + t446 * t498 + t450 * t466 + (-t425 - t467) * t518 + t582;
t416 = t426 * t518 + t434 * t491 - t450 * t465 + (-t446 - t494) * t497 + t584;
t415 = t425 * t497 + t433 * t465 - t434 * t466 + (-t426 - t468) * t498 + t589;
t1 = t465 * ((t507 * t428 + t471 * t430 + t472 * t432) * t465 + (t427 * t507 + t429 * t471 + t431 * t472) * t466 + (t447 * t507 + t448 * t471 + t449 * t472) * t491) / 0.2e1 + t466 * ((t428 * t505 + t430 * t469 + t432 * t470) * t465 + (t505 * t427 + t469 * t429 + t470 * t431) * t466 + (t447 * t505 + t448 * t469 + t449 * t470) * t491) / 0.2e1 + t491 * ((t428 * t527 + t430 * t499 + t432 * t500) * t465 + (t427 * t527 + t429 * t499 + t431 * t500) * t466 + (t527 * t447 + t499 * t448 + t500 * t449) * t491) / 0.2e1 + m(7) * (t415 ^ 2 + t416 ^ 2 + t417 ^ 2) / 0.2e1 + m(6) * (t418 ^ 2 + t419 ^ 2 + t420 ^ 2) / 0.2e1 + m(5) * (t421 ^ 2 + t422 ^ 2 + t423 ^ 2) / 0.2e1 + m(4) * (t424 ^ 2 + t443 ^ 2 + t444 ^ 2) / 0.2e1 + m(3) * (t496 ^ 2 + t513 ^ 2 + t514 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t539 * ((t587 * t479 - t524 * t481 + t525 * t483) * t539 + (t478 * t587 - t524 * t480 + t525 * t482) * t540 + (t509 * t587 - t524 * t510 + t525 * t511) * t549) / 0.2e1 + t540 * ((t479 * t588 - t522 * t481 + t523 * t483) * t539 + (t588 * t478 - t522 * t480 + t523 * t482) * t540 + (t509 * t588 - t522 * t510 + t523 * t511) * t549) / 0.2e1 + t549 * ((t479 * t586 - t546 * t481 + t547 * t483) * t539 + (t478 * t586 - t546 * t480 + t547 * t482) * t540 + (t586 * t509 - t546 * t510 + t547 * t511) * t549) / 0.2e1 - t619 * ((-t532 * t605 + t534 * t558 + t536 * t559) * t607 - (-t531 * t605 + t533 * t558 + t535 * t559) * t605 + (-t552 * t605 + t553 * t558 + t554 * t559) * t576) * t605 / 0.2e1 + ((t452 * t475 + t453 * t476 + t484 * t524 + t486 * t508 + t507 * t620) * t518 + (t437 * t475 + t439 * t476 + t456 * t524 + t460 * t508 + t507 * t622) * t498 + (t475 * t438 + t476 * t440 + t524 * t457 + t508 * t461 + t507 * t621) * t497) * t497 / 0.2e1 + ((t452 * t473 + t453 * t474 + t484 * t522 + t486 * t506 + t505 * t620) * t518 + (t473 * t437 + t474 * t439 + t522 * t456 + t506 * t460 + t505 * t622) * t498 + (t438 * t473 + t440 * t474 + t457 * t522 + t461 * t506 + t505 * t621) * t497) * t498 / 0.2e1 + ((t503 * t452 + t504 * t453 + t546 * t484 + t528 * t486 + t527 * t620) * t518 + (t437 * t503 + t439 * t504 + t456 * t546 + t460 * t528 + t527 * t622) * t498 + (t438 * t503 + t504 * t440 + t457 * t546 + t461 * t528 + t527 * t621) * t497) * t518 / 0.2e1 + (t576 * (t576 ^ 2 * t552 + (((t534 * t581 + t536 * t580) * t572 - (t533 * t581 + t535 * t580) * t575) * t573 + (-t531 * t575 + t532 * t572 + t553 * t581 + t554 * t580) * t576) * t573) + ((t532 * t607 + t534 * t560 + t536 * t561) * t607 - (t531 * t607 + t533 * t560 + t535 * t561) * t605 + (t552 * t607 + t553 * t560 + t554 * t561) * t576) * t607) * t619 / 0.2e1;
T  = t1;
