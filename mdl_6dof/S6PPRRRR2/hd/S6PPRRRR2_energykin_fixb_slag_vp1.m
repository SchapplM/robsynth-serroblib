% Calculate kinetic energy for
% S6PPRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRRR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:03:14
% EndTime: 2019-03-08 19:03:16
% DurationCPUTime: 1.72s
% Computational Cost: add. (5260->321), mult. (14021->506), div. (0->0), fcn. (18311->16), ass. (0->153)
t587 = cos(qJ(3));
t586 = cos(qJ(4));
t553 = cos(qJ(5));
t585 = pkin(5) * t553;
t583 = cos(pkin(6));
t582 = cos(pkin(7));
t581 = cos(pkin(13));
t580 = sin(pkin(7));
t579 = sin(pkin(13));
t547 = sin(pkin(12));
t549 = cos(pkin(12));
t568 = t583 * t579;
t534 = t547 * t581 + t549 * t568;
t552 = sin(qJ(3));
t570 = t583 * t581;
t561 = t547 * t579 - t549 * t570;
t557 = t561 * t582;
t548 = sin(pkin(6));
t571 = t548 * t580;
t566 = t587 * t571;
t513 = t534 * t552 + t549 * t566 + t557 * t587;
t550 = sin(qJ(5));
t578 = t513 * t550;
t535 = -t547 * t568 + t549 * t581;
t560 = t547 * t570 + t549 * t579;
t556 = t560 * t582;
t515 = t535 * t552 - t547 * t566 + t556 * t587;
t577 = t515 * t550;
t567 = t582 * t581;
t569 = t583 * t580;
t525 = -t569 * t587 + (t552 * t579 - t567 * t587) * t548;
t576 = t525 * t550;
t572 = t548 * t582;
t527 = -t549 * t572 + t561 * t580;
t523 = qJD(3) * t527;
t503 = qJD(4) * t513 + t523;
t528 = t547 * t572 + t560 * t580;
t524 = qJD(3) * t528;
t504 = qJD(4) * t515 + t524;
t533 = -t571 * t581 + t582 * t583;
t532 = qJD(3) * t533;
t519 = qJD(4) * t525 + t532;
t575 = qJD(2) * t548;
t540 = qJD(2) * t583 + qJD(1);
t514 = t534 * t587 + (-t549 * t571 - t557) * t552;
t551 = sin(qJ(4));
t497 = t514 * t551 - t527 * t586;
t463 = qJD(5) * t497 + t503;
t516 = t535 * t587 + (t547 * t571 - t556) * t552;
t499 = t516 * t551 - t528 * t586;
t464 = qJD(5) * t499 + t504;
t526 = t552 * t569 + (t552 * t567 + t579 * t587) * t548;
t517 = t526 * t551 - t533 * t586;
t491 = qJD(5) * t517 + t519;
t573 = t549 * t575;
t488 = pkin(3) * t514 + pkin(9) * t513;
t509 = pkin(3) * t526 + pkin(9) * t525;
t539 = t547 * t575;
t565 = -t488 * t532 + t509 * t523 + t539;
t489 = pkin(3) * t516 + pkin(9) * t515;
t564 = t488 * t524 - t489 * t523 + t540;
t563 = t489 * t532 - t509 * t524 - t573;
t498 = t514 * t586 + t527 * t551;
t461 = pkin(4) * t498 + pkin(10) * t497;
t518 = t526 * t586 + t533 * t551;
t490 = pkin(4) * t518 + pkin(10) * t517;
t562 = -t461 * t519 + t503 * t490 + t565;
t500 = t516 * t586 + t528 * t551;
t462 = pkin(4) * t500 + pkin(10) * t499;
t559 = t504 * t461 - t462 * t503 + t564;
t558 = t519 * t462 - t490 * t504 + t563;
t546 = qJ(5) + qJ(6);
t544 = cos(t546);
t543 = sin(t546);
t508 = rSges(4,1) * t526 - rSges(4,2) * t525 + rSges(4,3) * t533;
t507 = Icges(4,1) * t526 - Icges(4,4) * t525 + Icges(4,5) * t533;
t506 = Icges(4,4) * t526 - Icges(4,2) * t525 + Icges(4,6) * t533;
t505 = Icges(4,5) * t526 - Icges(4,6) * t525 + Icges(4,3) * t533;
t502 = t518 * t553 + t576;
t501 = -t518 * t550 + t525 * t553;
t494 = t518 * t544 + t525 * t543;
t493 = -t518 * t543 + t525 * t544;
t486 = rSges(5,1) * t518 - rSges(5,2) * t517 + rSges(5,3) * t525;
t485 = Icges(5,1) * t518 - Icges(5,4) * t517 + Icges(5,5) * t525;
t484 = Icges(5,4) * t518 - Icges(5,2) * t517 + Icges(5,6) * t525;
t483 = Icges(5,5) * t518 - Icges(5,6) * t517 + Icges(5,3) * t525;
t481 = rSges(4,1) * t516 - rSges(4,2) * t515 + rSges(4,3) * t528;
t480 = rSges(4,1) * t514 - rSges(4,2) * t513 + rSges(4,3) * t527;
t479 = Icges(4,1) * t516 - Icges(4,4) * t515 + Icges(4,5) * t528;
t478 = Icges(4,1) * t514 - Icges(4,4) * t513 + Icges(4,5) * t527;
t477 = Icges(4,4) * t516 - Icges(4,2) * t515 + Icges(4,6) * t528;
t476 = Icges(4,4) * t514 - Icges(4,2) * t513 + Icges(4,6) * t527;
t475 = Icges(4,5) * t516 - Icges(4,6) * t515 + Icges(4,3) * t528;
t474 = Icges(4,5) * t514 - Icges(4,6) * t513 + Icges(4,3) * t527;
t473 = t500 * t553 + t577;
t472 = -t500 * t550 + t515 * t553;
t471 = t498 * t553 + t578;
t470 = -t498 * t550 + t513 * t553;
t469 = t500 * t544 + t515 * t543;
t468 = -t500 * t543 + t515 * t544;
t467 = t498 * t544 + t513 * t543;
t466 = -t498 * t543 + t513 * t544;
t465 = qJD(6) * t517 + t491;
t458 = rSges(6,1) * t502 + rSges(6,2) * t501 + rSges(6,3) * t517;
t457 = Icges(6,1) * t502 + Icges(6,4) * t501 + Icges(6,5) * t517;
t456 = Icges(6,4) * t502 + Icges(6,2) * t501 + Icges(6,6) * t517;
t455 = Icges(6,5) * t502 + Icges(6,6) * t501 + Icges(6,3) * t517;
t454 = rSges(5,1) * t500 - rSges(5,2) * t499 + rSges(5,3) * t515;
t453 = rSges(5,1) * t498 - rSges(5,2) * t497 + rSges(5,3) * t513;
t452 = Icges(5,1) * t500 - Icges(5,4) * t499 + Icges(5,5) * t515;
t451 = Icges(5,1) * t498 - Icges(5,4) * t497 + Icges(5,5) * t513;
t450 = Icges(5,4) * t500 - Icges(5,2) * t499 + Icges(5,6) * t515;
t449 = Icges(5,4) * t498 - Icges(5,2) * t497 + Icges(5,6) * t513;
t448 = Icges(5,5) * t500 - Icges(5,6) * t499 + Icges(5,3) * t515;
t447 = Icges(5,5) * t498 - Icges(5,6) * t497 + Icges(5,3) * t513;
t446 = rSges(7,1) * t494 + rSges(7,2) * t493 + rSges(7,3) * t517;
t445 = Icges(7,1) * t494 + Icges(7,4) * t493 + Icges(7,5) * t517;
t444 = Icges(7,4) * t494 + Icges(7,2) * t493 + Icges(7,6) * t517;
t443 = Icges(7,5) * t494 + Icges(7,6) * t493 + Icges(7,3) * t517;
t442 = pkin(5) * t576 + pkin(11) * t517 + t518 * t585;
t441 = qJD(6) * t499 + t464;
t440 = qJD(6) * t497 + t463;
t438 = -t573 + (t481 * t533 - t508 * t528) * qJD(3);
t437 = t539 + (-t480 * t533 + t508 * t527) * qJD(3);
t436 = (t480 * t528 - t481 * t527) * qJD(3) + t540;
t435 = rSges(6,1) * t473 + rSges(6,2) * t472 + rSges(6,3) * t499;
t434 = rSges(6,1) * t471 + rSges(6,2) * t470 + rSges(6,3) * t497;
t433 = Icges(6,1) * t473 + Icges(6,4) * t472 + Icges(6,5) * t499;
t432 = Icges(6,1) * t471 + Icges(6,4) * t470 + Icges(6,5) * t497;
t431 = Icges(6,4) * t473 + Icges(6,2) * t472 + Icges(6,6) * t499;
t430 = Icges(6,4) * t471 + Icges(6,2) * t470 + Icges(6,6) * t497;
t429 = Icges(6,5) * t473 + Icges(6,6) * t472 + Icges(6,3) * t499;
t428 = Icges(6,5) * t471 + Icges(6,6) * t470 + Icges(6,3) * t497;
t427 = rSges(7,1) * t469 + rSges(7,2) * t468 + rSges(7,3) * t499;
t426 = rSges(7,1) * t467 + rSges(7,2) * t466 + rSges(7,3) * t497;
t425 = Icges(7,1) * t469 + Icges(7,4) * t468 + Icges(7,5) * t499;
t424 = Icges(7,1) * t467 + Icges(7,4) * t466 + Icges(7,5) * t497;
t423 = Icges(7,4) * t469 + Icges(7,2) * t468 + Icges(7,6) * t499;
t422 = Icges(7,4) * t467 + Icges(7,2) * t466 + Icges(7,6) * t497;
t421 = Icges(7,5) * t469 + Icges(7,6) * t468 + Icges(7,3) * t499;
t420 = Icges(7,5) * t467 + Icges(7,6) * t466 + Icges(7,3) * t497;
t419 = pkin(5) * t577 + pkin(11) * t499 + t500 * t585;
t418 = pkin(5) * t578 + pkin(11) * t497 + t498 * t585;
t417 = t454 * t519 - t486 * t504 + t563;
t416 = -t453 * t519 + t486 * t503 + t565;
t415 = t453 * t504 - t454 * t503 + t564;
t414 = t435 * t491 - t458 * t464 + t558;
t413 = -t434 * t491 + t458 * t463 + t562;
t412 = t434 * t464 - t435 * t463 + t559;
t411 = t419 * t491 + t427 * t465 - t441 * t446 - t442 * t464 + t558;
t410 = -t418 * t491 - t426 * t465 + t440 * t446 + t442 * t463 + t562;
t409 = t418 * t464 - t419 * t463 + t426 * t441 - t427 * t440 + t559;
t1 = m(7) * (t409 ^ 2 + t410 ^ 2 + t411 ^ 2) / 0.2e1 + m(6) * (t412 ^ 2 + t413 ^ 2 + t414 ^ 2) / 0.2e1 + m(5) * (t415 ^ 2 + t416 ^ 2 + t417 ^ 2) / 0.2e1 + m(4) * (t436 ^ 2 + t437 ^ 2 + t438 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t540 ^ 2 + (t547 ^ 2 + t549 ^ 2) * qJD(2) ^ 2 * t548 ^ 2) / 0.2e1 + t441 * ((t499 * t421 + t468 * t423 + t469 * t425) * t441 + (t420 * t499 + t422 * t468 + t424 * t469) * t440 + (t443 * t499 + t444 * t468 + t445 * t469) * t465) / 0.2e1 + t440 * ((t421 * t497 + t423 * t466 + t425 * t467) * t441 + (t497 * t420 + t466 * t422 + t467 * t424) * t440 + (t443 * t497 + t444 * t466 + t445 * t467) * t465) / 0.2e1 + t465 * ((t421 * t517 + t423 * t493 + t425 * t494) * t441 + (t420 * t517 + t422 * t493 + t424 * t494) * t440 + (t517 * t443 + t493 * t444 + t494 * t445) * t465) / 0.2e1 + t464 * ((t499 * t429 + t472 * t431 + t473 * t433) * t464 + (t428 * t499 + t430 * t472 + t432 * t473) * t463 + (t455 * t499 + t456 * t472 + t457 * t473) * t491) / 0.2e1 + t463 * ((t429 * t497 + t431 * t470 + t433 * t471) * t464 + (t497 * t428 + t470 * t430 + t471 * t432) * t463 + (t455 * t497 + t456 * t470 + t457 * t471) * t491) / 0.2e1 + t491 * ((t429 * t517 + t431 * t501 + t433 * t502) * t464 + (t428 * t517 + t430 * t501 + t432 * t502) * t463 + (t517 * t455 + t501 * t456 + t502 * t457) * t491) / 0.2e1 + t504 * ((t515 * t448 - t499 * t450 + t500 * t452) * t504 + (t447 * t515 - t449 * t499 + t451 * t500) * t503 + (t483 * t515 - t484 * t499 + t485 * t500) * t519) / 0.2e1 + t503 * ((t448 * t513 - t450 * t497 + t452 * t498) * t504 + (t513 * t447 - t497 * t449 + t498 * t451) * t503 + (t483 * t513 - t484 * t497 + t485 * t498) * t519) / 0.2e1 + t519 * ((t448 * t525 - t450 * t517 + t452 * t518) * t504 + (t447 * t525 - t449 * t517 + t451 * t518) * t503 + (t525 * t483 - t517 * t484 + t518 * t485) * t519) / 0.2e1 + (t528 * ((t475 * t528 - t477 * t515 + t479 * t516) * t528 + (t474 * t528 - t476 * t515 + t478 * t516) * t527 + (t505 * t528 - t506 * t515 + t507 * t516) * t533) + t527 * ((t475 * t527 - t477 * t513 + t479 * t514) * t528 + (t474 * t527 - t476 * t513 + t478 * t514) * t527 + (t505 * t527 - t506 * t513 + t507 * t514) * t533) + t533 * ((t475 * t533 - t477 * t525 + t479 * t526) * t528 + (t474 * t533 - t476 * t525 + t478 * t526) * t527 + (t505 * t533 - t506 * t525 + t507 * t526) * t533)) * qJD(3) ^ 2 / 0.2e1;
T  = t1;
