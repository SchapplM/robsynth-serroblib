% Calculate kinetic energy for
% S6RRRRPR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-10 00:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR14_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR14_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR14_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR14_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR14_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR14_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR14_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:08:41
% EndTime: 2019-03-10 00:08:44
% DurationCPUTime: 3.59s
% Computational Cost: add. (5429->387), mult. (14277->589), div. (0->0), fcn. (18380->16), ass. (0->176)
t634 = Icges(5,2) + Icges(6,3);
t578 = sin(pkin(6));
t580 = cos(pkin(6));
t586 = cos(qJ(2));
t587 = cos(qJ(1));
t611 = t586 * t587;
t584 = sin(qJ(2));
t585 = sin(qJ(1));
t614 = t584 * t585;
t559 = t580 * t611 - t614;
t612 = t585 * t586;
t613 = t584 * t587;
t560 = t580 * t613 + t612;
t561 = -t580 * t612 - t613;
t562 = -t580 * t614 + t611;
t615 = t578 * t587;
t616 = t578 * t585;
t601 = (Icges(3,5) * t560 + Icges(3,6) * t559 - Icges(3,3) * t615) * t587 - (Icges(3,5) * t562 + Icges(3,6) * t561 + Icges(3,3) * t616) * t585;
t633 = t578 * t601;
t583 = sin(qJ(3));
t621 = sin(pkin(7));
t604 = t578 * t621;
t622 = cos(pkin(7));
t625 = cos(qJ(3));
t527 = t560 * t625 + (t559 * t622 - t587 * t604) * t583;
t582 = sin(qJ(4));
t605 = t578 * t622;
t596 = -t559 * t621 - t587 * t605;
t624 = cos(qJ(4));
t509 = t527 * t624 + t582 * t596;
t602 = t625 * t621;
t600 = t578 * t602;
t603 = t622 * t625;
t526 = -t559 * t603 + t560 * t583 + t587 * t600;
t577 = sin(pkin(13));
t579 = cos(pkin(13));
t474 = -t509 * t577 + t526 * t579;
t620 = t526 * t577;
t475 = t509 * t579 + t620;
t508 = t527 * t582 - t596 * t624;
t632 = -Icges(5,4) * t509 + Icges(6,5) * t475 - Icges(5,6) * t526 + Icges(6,6) * t474 + t634 * t508;
t529 = t562 * t625 + (t561 * t622 + t585 * t604) * t583;
t595 = -t561 * t621 + t585 * t605;
t511 = t529 * t624 + t582 * t595;
t528 = -t561 * t603 + t562 * t583 - t585 * t600;
t476 = -t511 * t577 + t528 * t579;
t619 = t528 * t577;
t477 = t511 * t579 + t619;
t510 = t529 * t582 - t595 * t624;
t631 = -Icges(5,4) * t511 + Icges(6,5) * t477 - Icges(5,6) * t528 + Icges(6,6) * t476 + t634 * t510;
t548 = t580 * t621 * t583 + (t583 * t586 * t622 + t584 * t625) * t578;
t594 = t580 * t622 - t586 * t604;
t525 = t548 * t624 + t582 * t594;
t617 = t578 * t584;
t547 = -t578 * t586 * t603 - t580 * t602 + t583 * t617;
t502 = -t525 * t577 + t547 * t579;
t618 = t547 * t577;
t503 = t525 * t579 + t618;
t524 = t548 * t582 - t594 * t624;
t630 = -Icges(5,4) * t525 + Icges(6,5) * t503 - Icges(5,6) * t547 + Icges(6,6) * t502 + t634 * t524;
t623 = pkin(5) * t579;
t530 = t560 * pkin(2) + pkin(10) * t596;
t531 = t562 * pkin(2) + pkin(10) * t595;
t607 = qJD(2) * t578;
t570 = t585 * t607;
t606 = t587 * t607;
t609 = t530 * t570 + t531 * t606;
t540 = qJD(3) * t595 + t570;
t608 = qJD(1) * (pkin(1) * t585 - pkin(9) * t615);
t571 = qJD(2) * t580 + qJD(1);
t500 = qJD(4) * t528 + t540;
t549 = qJD(3) * t594 + t571;
t517 = qJD(4) * t547 + t549;
t541 = qJD(3) * t596 - t606;
t494 = pkin(3) * t527 + pkin(11) * t526;
t495 = pkin(3) * t529 + pkin(11) * t528;
t599 = t540 * t494 - t495 * t541 + t609;
t501 = qJD(4) * t526 + t541;
t550 = pkin(2) * t617 + pkin(10) * t594;
t563 = qJD(1) * (pkin(1) * t587 + pkin(9) * t616);
t598 = t571 * t531 - t550 * t570 + t563;
t468 = pkin(4) * t509 + qJ(5) * t508;
t597 = qJD(5) * t524 + t500 * t468 + t599;
t593 = -t530 * t571 - t550 * t606 - t608;
t516 = pkin(3) * t548 + pkin(11) * t547;
t592 = t549 * t495 - t516 * t540 + t598;
t469 = pkin(4) * t511 + qJ(5) * t510;
t591 = qJD(5) * t508 + t517 * t469 + t592;
t590 = -t494 * t549 + t541 * t516 + t593;
t493 = pkin(4) * t525 + qJ(5) * t524;
t589 = qJD(5) * t510 + t501 * t493 + t590;
t576 = pkin(13) + qJ(6);
t575 = cos(t576);
t574 = sin(t576);
t568 = rSges(2,1) * t587 - rSges(2,2) * t585;
t567 = rSges(2,1) * t585 + rSges(2,2) * t587;
t556 = rSges(3,3) * t580 + (rSges(3,1) * t584 + rSges(3,2) * t586) * t578;
t555 = Icges(3,5) * t580 + (Icges(3,1) * t584 + Icges(3,4) * t586) * t578;
t554 = Icges(3,6) * t580 + (Icges(3,4) * t584 + Icges(3,2) * t586) * t578;
t553 = Icges(3,3) * t580 + (Icges(3,5) * t584 + Icges(3,6) * t586) * t578;
t539 = rSges(3,1) * t562 + rSges(3,2) * t561 + rSges(3,3) * t616;
t538 = rSges(3,1) * t560 + rSges(3,2) * t559 - rSges(3,3) * t615;
t537 = Icges(3,1) * t562 + Icges(3,4) * t561 + Icges(3,5) * t616;
t536 = Icges(3,1) * t560 + Icges(3,4) * t559 - Icges(3,5) * t615;
t535 = Icges(3,4) * t562 + Icges(3,2) * t561 + Icges(3,6) * t616;
t534 = Icges(3,4) * t560 + Icges(3,2) * t559 - Icges(3,6) * t615;
t515 = t548 * rSges(4,1) - t547 * rSges(4,2) + rSges(4,3) * t594;
t514 = Icges(4,1) * t548 - Icges(4,4) * t547 + Icges(4,5) * t594;
t513 = Icges(4,4) * t548 - Icges(4,2) * t547 + Icges(4,6) * t594;
t512 = Icges(4,5) * t548 - Icges(4,6) * t547 + Icges(4,3) * t594;
t505 = t539 * t571 - t556 * t570 + t563;
t504 = -t538 * t571 - t556 * t606 - t608;
t499 = t525 * t575 + t547 * t574;
t498 = -t525 * t574 + t547 * t575;
t497 = (t538 * t585 + t539 * t587) * t607;
t492 = qJD(6) * t524 + t517;
t490 = t529 * rSges(4,1) - t528 * rSges(4,2) + rSges(4,3) * t595;
t489 = t527 * rSges(4,1) - t526 * rSges(4,2) + rSges(4,3) * t596;
t488 = Icges(4,1) * t529 - Icges(4,4) * t528 + Icges(4,5) * t595;
t487 = Icges(4,1) * t527 - Icges(4,4) * t526 + Icges(4,5) * t596;
t486 = Icges(4,4) * t529 - Icges(4,2) * t528 + Icges(4,6) * t595;
t485 = Icges(4,4) * t527 - Icges(4,2) * t526 + Icges(4,6) * t596;
t484 = Icges(4,5) * t529 - Icges(4,6) * t528 + Icges(4,3) * t595;
t483 = Icges(4,5) * t527 - Icges(4,6) * t526 + Icges(4,3) * t596;
t482 = rSges(5,1) * t525 - rSges(5,2) * t524 + rSges(5,3) * t547;
t481 = Icges(5,1) * t525 - Icges(5,4) * t524 + Icges(5,5) * t547;
t479 = Icges(5,5) * t525 - Icges(5,6) * t524 + Icges(5,3) * t547;
t473 = t511 * t575 + t528 * t574;
t472 = -t511 * t574 + t528 * t575;
t471 = t509 * t575 + t526 * t574;
t470 = -t509 * t574 + t526 * t575;
t467 = qJD(6) * t508 + t501;
t466 = qJD(6) * t510 + t500;
t464 = rSges(5,1) * t511 - rSges(5,2) * t510 + rSges(5,3) * t528;
t463 = rSges(5,1) * t509 - rSges(5,2) * t508 + rSges(5,3) * t526;
t462 = Icges(5,1) * t511 - Icges(5,4) * t510 + Icges(5,5) * t528;
t461 = Icges(5,1) * t509 - Icges(5,4) * t508 + Icges(5,5) * t526;
t458 = Icges(5,5) * t511 - Icges(5,6) * t510 + Icges(5,3) * t528;
t457 = Icges(5,5) * t509 - Icges(5,6) * t508 + Icges(5,3) * t526;
t455 = rSges(6,1) * t503 + rSges(6,2) * t502 + rSges(6,3) * t524;
t454 = Icges(6,1) * t503 + Icges(6,4) * t502 + Icges(6,5) * t524;
t453 = Icges(6,4) * t503 + Icges(6,2) * t502 + Icges(6,6) * t524;
t451 = rSges(7,1) * t499 + rSges(7,2) * t498 + rSges(7,3) * t524;
t450 = Icges(7,1) * t499 + Icges(7,4) * t498 + Icges(7,5) * t524;
t449 = Icges(7,4) * t499 + Icges(7,2) * t498 + Icges(7,6) * t524;
t448 = Icges(7,5) * t499 + Icges(7,6) * t498 + Icges(7,3) * t524;
t447 = pkin(5) * t618 + pkin(12) * t524 + t525 * t623;
t445 = t490 * t549 - t515 * t540 + t598;
t444 = -t489 * t549 + t515 * t541 + t593;
t443 = rSges(6,1) * t477 + rSges(6,2) * t476 + rSges(6,3) * t510;
t442 = rSges(6,1) * t475 + rSges(6,2) * t474 + rSges(6,3) * t508;
t441 = Icges(6,1) * t477 + Icges(6,4) * t476 + Icges(6,5) * t510;
t440 = Icges(6,1) * t475 + Icges(6,4) * t474 + Icges(6,5) * t508;
t439 = Icges(6,4) * t477 + Icges(6,2) * t476 + Icges(6,6) * t510;
t438 = Icges(6,4) * t475 + Icges(6,2) * t474 + Icges(6,6) * t508;
t435 = rSges(7,1) * t473 + rSges(7,2) * t472 + rSges(7,3) * t510;
t434 = rSges(7,1) * t471 + rSges(7,2) * t470 + rSges(7,3) * t508;
t433 = Icges(7,1) * t473 + Icges(7,4) * t472 + Icges(7,5) * t510;
t432 = Icges(7,1) * t471 + Icges(7,4) * t470 + Icges(7,5) * t508;
t431 = Icges(7,4) * t473 + Icges(7,2) * t472 + Icges(7,6) * t510;
t430 = Icges(7,4) * t471 + Icges(7,2) * t470 + Icges(7,6) * t508;
t429 = Icges(7,5) * t473 + Icges(7,6) * t472 + Icges(7,3) * t510;
t428 = Icges(7,5) * t471 + Icges(7,6) * t470 + Icges(7,3) * t508;
t427 = pkin(5) * t619 + pkin(12) * t510 + t511 * t623;
t426 = pkin(5) * t620 + pkin(12) * t508 + t509 * t623;
t425 = t489 * t540 - t490 * t541 + t609;
t424 = t464 * t517 - t482 * t500 + t592;
t423 = -t463 * t517 + t482 * t501 + t590;
t422 = t463 * t500 - t464 * t501 + t599;
t421 = t443 * t517 + (-t455 - t493) * t500 + t591;
t420 = t455 * t501 + (-t442 - t468) * t517 + t589;
t419 = t442 * t500 + (-t443 - t469) * t501 + t597;
t418 = t427 * t517 + t435 * t492 - t451 * t466 + (-t447 - t493) * t500 + t591;
t417 = -t434 * t492 + t447 * t501 + t451 * t467 + (-t426 - t468) * t517 + t589;
t416 = t426 * t500 + t434 * t466 - t435 * t467 + (-t427 - t469) * t501 + t597;
t1 = m(7) * (t416 ^ 2 + t417 ^ 2 + t418 ^ 2) / 0.2e1 + m(6) * (t419 ^ 2 + t420 ^ 2 + t421 ^ 2) / 0.2e1 + m(5) * (t422 ^ 2 + t423 ^ 2 + t424 ^ 2) / 0.2e1 + m(4) * (t425 ^ 2 + t444 ^ 2 + t445 ^ 2) / 0.2e1 + m(3) * (t497 ^ 2 + t504 ^ 2 + t505 ^ 2) / 0.2e1 + t466 * ((t510 * t429 + t472 * t431 + t473 * t433) * t466 + (t428 * t510 + t430 * t472 + t432 * t473) * t467 + (t448 * t510 + t449 * t472 + t450 * t473) * t492) / 0.2e1 + t467 * ((t429 * t508 + t431 * t470 + t433 * t471) * t466 + (t508 * t428 + t470 * t430 + t471 * t432) * t467 + (t448 * t508 + t449 * t470 + t450 * t471) * t492) / 0.2e1 + t492 * ((t429 * t524 + t431 * t498 + t433 * t499) * t466 + (t428 * t524 + t430 * t498 + t432 * t499) * t467 + (t448 * t524 + t449 * t498 + t450 * t499) * t492) / 0.2e1 + t549 * ((t484 * t594 - t547 * t486 + t548 * t488) * t540 + (t483 * t594 - t547 * t485 + t548 * t487) * t541 + (t594 * t512 - t547 * t513 + t548 * t514) * t549) / 0.2e1 + t541 * ((t484 * t596 - t526 * t486 + t527 * t488) * t540 + (t596 * t483 - t526 * t485 + t527 * t487) * t541 + (t512 * t596 - t526 * t513 + t527 * t514) * t549) / 0.2e1 + t540 * ((t595 * t484 - t528 * t486 + t529 * t488) * t540 + (t483 * t595 - t528 * t485 + t529 * t487) * t541 + (t512 * t595 - t528 * t513 + t529 * t514) * t549) / 0.2e1 + t571 * ((t580 * t553 + (t554 * t586 + t555 * t584) * t578) * t571 + (((t535 * t586 + t537 * t584) * t585 - (t534 * t586 + t536 * t584) * t587) * t578 - t601 * t580) * t607) / 0.2e1 - ((-t553 * t615 + t554 * t559 + t555 * t560) * t571 + ((t535 * t559 + t537 * t560) * t585 + (-t559 * t534 - t560 * t536 + t633) * t587) * t607) * t606 / 0.2e1 + ((t553 * t616 + t554 * t561 + t555 * t562) * t571 + (-(t534 * t561 + t536 * t562) * t587 + (t561 * t535 + t562 * t537 - t633) * t585) * t607) * t570 / 0.2e1 + ((t453 * t476 + t454 * t477 + t479 * t528 + t481 * t511 + t630 * t510) * t517 + (t438 * t476 + t440 * t477 + t457 * t528 + t461 * t511 + t632 * t510) * t501 + (t439 * t476 + t441 * t477 + t458 * t528 + t462 * t511 + t631 * t510) * t500) * t500 / 0.2e1 + ((t453 * t474 + t454 * t475 + t479 * t526 + t481 * t509 + t630 * t508) * t517 + (t438 * t474 + t440 * t475 + t457 * t526 + t461 * t509 + t632 * t508) * t501 + (t439 * t474 + t441 * t475 + t458 * t526 + t462 * t509 + t631 * t508) * t500) * t501 / 0.2e1 + ((t453 * t502 + t454 * t503 + t479 * t547 + t481 * t525 + t630 * t524) * t517 + (t438 * t502 + t440 * t503 + t457 * t547 + t461 * t525 + t632 * t524) * t501 + (t439 * t502 + t441 * t503 + t458 * t547 + t462 * t525 + t631 * t524) * t500) * t517 / 0.2e1 + (m(2) * (t567 ^ 2 + t568 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
