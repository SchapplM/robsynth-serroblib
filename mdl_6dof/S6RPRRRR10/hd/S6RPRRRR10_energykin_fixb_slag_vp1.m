% Calculate kinetic energy for
% S6RPRRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR10_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR10_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:27:40
% EndTime: 2019-03-09 07:27:42
% DurationCPUTime: 2.20s
% Computational Cost: add. (5072->348), mult. (12495->546), div. (0->0), fcn. (16133->16), ass. (0->163)
t562 = cos(pkin(6));
t603 = cos(pkin(13));
t607 = sin(qJ(1));
t587 = t607 * t603;
t568 = cos(qJ(1));
t601 = sin(pkin(13));
t589 = t568 * t601;
t577 = t562 * t587 + t589;
t561 = sin(pkin(6));
t602 = sin(pkin(7));
t592 = t561 * t602;
t604 = cos(pkin(7));
t610 = t577 * t604 - t607 * t592;
t608 = cos(qJ(3));
t567 = cos(qJ(4));
t606 = pkin(4) * t567;
t586 = t607 * t601;
t590 = t568 * t603;
t576 = -t562 * t590 + t586;
t593 = t561 * t604;
t540 = -t568 * t593 + t576 * t602;
t564 = sin(qJ(4));
t600 = t540 * t564;
t541 = t577 * t602 + t593 * t607;
t599 = t541 * t564;
t547 = t562 * t604 - t592 * t603;
t598 = t547 * t564;
t597 = t561 * t568;
t596 = qJ(4) + qJ(5);
t548 = t562 * t589 + t587;
t565 = sin(qJ(3));
t572 = t576 * t604;
t588 = t608 * t602;
t525 = t548 * t565 + t572 * t608 + t588 * t597;
t536 = qJD(3) * t540;
t507 = qJD(4) * t525 + t536;
t549 = -t562 * t586 + t590;
t527 = t549 * t565 + t610 * t608;
t537 = qJD(3) * t541;
t508 = qJD(4) * t527 + t537;
t545 = qJD(3) * t547 + qJD(1);
t481 = qJD(5) * t525 + t507;
t482 = qJD(5) * t527 + t508;
t595 = t561 * t607;
t584 = t604 * t603;
t591 = t561 * t601;
t538 = -t561 * t584 * t608 - t562 * t588 + t565 * t591;
t529 = qJD(4) * t538 + t545;
t594 = cos(t596);
t509 = qJD(5) * t538 + t529;
t585 = -qJD(2) * t597 + qJD(1) * (t568 * pkin(1) + qJ(2) * t595);
t551 = pkin(1) * t607 - qJ(2) * t597;
t557 = qJD(2) * t595;
t583 = t557 + (-t548 * pkin(2) - pkin(9) * t540 - t551) * qJD(1);
t581 = qJD(1) * (t549 * pkin(2) + pkin(9) * t541) + t585;
t526 = t548 * t608 + (-t568 * t592 - t572) * t565;
t494 = t526 * pkin(3) + t525 * pkin(10);
t528 = t549 * t608 - t610 * t565;
t495 = t528 * pkin(3) + t527 * pkin(10);
t559 = qJD(2) * t562;
t580 = t494 * t537 - t495 * t536 + t559;
t539 = t562 * t602 * t565 + (t565 * t584 + t601 * t608) * t561;
t516 = t539 * pkin(3) + t538 * pkin(10);
t579 = -t545 * t494 + t516 * t536 + t583;
t440 = pkin(4) * t600 + pkin(11) * t525 + t526 * t606;
t441 = pkin(4) * t599 + pkin(11) * t527 + t528 * t606;
t578 = t508 * t440 - t441 * t507 + t580;
t575 = t545 * t495 - t516 * t537 + t581;
t467 = pkin(4) * t598 + pkin(11) * t538 + t539 * t606;
t574 = -t529 * t440 + t507 * t467 + t579;
t571 = t529 * t441 - t467 * t508 + t575;
t566 = cos(qJ(6));
t563 = sin(qJ(6));
t560 = sin(t596);
t555 = t568 * rSges(2,1) - rSges(2,2) * t607;
t554 = rSges(2,1) * t607 + t568 * rSges(2,2);
t524 = t539 * t567 + t598;
t523 = -t539 * t564 + t547 * t567;
t518 = t539 * t594 + t547 * t560;
t517 = t539 * t560 - t547 * t594;
t515 = qJD(1) * (t549 * rSges(3,1) - rSges(3,2) * t577 + rSges(3,3) * t595) + t585;
t514 = t557 + (-t548 * rSges(3,1) + rSges(3,2) * t576 + rSges(3,3) * t597 - t551) * qJD(1);
t513 = rSges(4,1) * t539 - rSges(4,2) * t538 + rSges(4,3) * t547;
t512 = Icges(4,1) * t539 - Icges(4,4) * t538 + Icges(4,5) * t547;
t511 = Icges(4,4) * t539 - Icges(4,2) * t538 + Icges(4,6) * t547;
t510 = Icges(4,5) * t539 - Icges(4,6) * t538 + Icges(4,3) * t547;
t506 = t528 * t567 + t599;
t505 = -t528 * t564 + t541 * t567;
t504 = t526 * t567 + t600;
t503 = -t526 * t564 + t540 * t567;
t502 = t528 * t594 + t541 * t560;
t501 = t528 * t560 - t541 * t594;
t500 = t526 * t594 + t540 * t560;
t499 = t526 * t560 - t540 * t594;
t498 = t518 * t566 + t538 * t563;
t497 = -t518 * t563 + t538 * t566;
t493 = pkin(5) * t518 + pkin(12) * t517;
t490 = rSges(4,1) * t528 - rSges(4,2) * t527 + rSges(4,3) * t541;
t489 = rSges(4,1) * t526 - rSges(4,2) * t525 + rSges(4,3) * t540;
t488 = Icges(4,1) * t528 - Icges(4,4) * t527 + Icges(4,5) * t541;
t487 = Icges(4,1) * t526 - Icges(4,4) * t525 + Icges(4,5) * t540;
t486 = Icges(4,4) * t528 - Icges(4,2) * t527 + Icges(4,6) * t541;
t485 = Icges(4,4) * t526 - Icges(4,2) * t525 + Icges(4,6) * t540;
t484 = Icges(4,5) * t528 - Icges(4,6) * t527 + Icges(4,3) * t541;
t483 = Icges(4,5) * t526 - Icges(4,6) * t525 + Icges(4,3) * t540;
t480 = rSges(5,1) * t524 + rSges(5,2) * t523 + rSges(5,3) * t538;
t479 = Icges(5,1) * t524 + Icges(5,4) * t523 + Icges(5,5) * t538;
t478 = Icges(5,4) * t524 + Icges(5,2) * t523 + Icges(5,6) * t538;
t477 = Icges(5,5) * t524 + Icges(5,6) * t523 + Icges(5,3) * t538;
t476 = qJD(6) * t517 + t509;
t475 = rSges(6,1) * t518 - rSges(6,2) * t517 + rSges(6,3) * t538;
t474 = Icges(6,1) * t518 - Icges(6,4) * t517 + Icges(6,5) * t538;
t473 = Icges(6,4) * t518 - Icges(6,2) * t517 + Icges(6,6) * t538;
t472 = Icges(6,5) * t518 - Icges(6,6) * t517 + Icges(6,3) * t538;
t471 = t502 * t566 + t527 * t563;
t470 = -t502 * t563 + t527 * t566;
t469 = t500 * t566 + t525 * t563;
t468 = -t500 * t563 + t525 * t566;
t466 = pkin(5) * t502 + pkin(12) * t501;
t465 = pkin(5) * t500 + pkin(12) * t499;
t464 = qJD(6) * t501 + t482;
t463 = qJD(6) * t499 + t481;
t462 = rSges(5,1) * t506 + rSges(5,2) * t505 + rSges(5,3) * t527;
t461 = rSges(5,1) * t504 + rSges(5,2) * t503 + rSges(5,3) * t525;
t460 = Icges(5,1) * t506 + Icges(5,4) * t505 + Icges(5,5) * t527;
t459 = Icges(5,1) * t504 + Icges(5,4) * t503 + Icges(5,5) * t525;
t458 = Icges(5,4) * t506 + Icges(5,2) * t505 + Icges(5,6) * t527;
t457 = Icges(5,4) * t504 + Icges(5,2) * t503 + Icges(5,6) * t525;
t456 = Icges(5,5) * t506 + Icges(5,6) * t505 + Icges(5,3) * t527;
t455 = Icges(5,5) * t504 + Icges(5,6) * t503 + Icges(5,3) * t525;
t453 = rSges(6,1) * t502 - rSges(6,2) * t501 + rSges(6,3) * t527;
t452 = rSges(6,1) * t500 - rSges(6,2) * t499 + rSges(6,3) * t525;
t451 = Icges(6,1) * t502 - Icges(6,4) * t501 + Icges(6,5) * t527;
t450 = Icges(6,1) * t500 - Icges(6,4) * t499 + Icges(6,5) * t525;
t449 = Icges(6,4) * t502 - Icges(6,2) * t501 + Icges(6,6) * t527;
t448 = Icges(6,4) * t500 - Icges(6,2) * t499 + Icges(6,6) * t525;
t447 = Icges(6,5) * t502 - Icges(6,6) * t501 + Icges(6,3) * t527;
t446 = Icges(6,5) * t500 - Icges(6,6) * t499 + Icges(6,3) * t525;
t445 = rSges(7,1) * t498 + rSges(7,2) * t497 + rSges(7,3) * t517;
t444 = Icges(7,1) * t498 + Icges(7,4) * t497 + Icges(7,5) * t517;
t443 = Icges(7,4) * t498 + Icges(7,2) * t497 + Icges(7,6) * t517;
t442 = Icges(7,5) * t498 + Icges(7,6) * t497 + Icges(7,3) * t517;
t437 = t490 * t545 - t513 * t537 + t581;
t436 = -t545 * t489 + t513 * t536 + t583;
t435 = t559 + (t489 * t541 - t490 * t540) * qJD(3);
t434 = rSges(7,1) * t471 + rSges(7,2) * t470 + rSges(7,3) * t501;
t433 = rSges(7,1) * t469 + rSges(7,2) * t468 + rSges(7,3) * t499;
t432 = Icges(7,1) * t471 + Icges(7,4) * t470 + Icges(7,5) * t501;
t431 = Icges(7,1) * t469 + Icges(7,4) * t468 + Icges(7,5) * t499;
t430 = Icges(7,4) * t471 + Icges(7,2) * t470 + Icges(7,6) * t501;
t429 = Icges(7,4) * t469 + Icges(7,2) * t468 + Icges(7,6) * t499;
t428 = Icges(7,5) * t471 + Icges(7,6) * t470 + Icges(7,3) * t501;
t427 = Icges(7,5) * t469 + Icges(7,6) * t468 + Icges(7,3) * t499;
t426 = t462 * t529 - t480 * t508 + t575;
t425 = -t529 * t461 + t507 * t480 + t579;
t424 = t461 * t508 - t462 * t507 + t580;
t423 = t453 * t509 - t475 * t482 + t571;
t422 = -t509 * t452 + t481 * t475 + t574;
t421 = t452 * t482 - t453 * t481 + t578;
t420 = t434 * t476 - t445 * t464 + t466 * t509 - t482 * t493 + t571;
t419 = -t476 * t433 + t463 * t445 - t509 * t465 + t481 * t493 + t574;
t418 = t433 * t464 - t434 * t463 + t465 * t482 - t466 * t481 + t578;
t1 = m(7) * (t418 ^ 2 + t419 ^ 2 + t420 ^ 2) / 0.2e1 + m(6) * (t421 ^ 2 + t422 ^ 2 + t423 ^ 2) / 0.2e1 + m(5) * (t424 ^ 2 + t425 ^ 2 + t426 ^ 2) / 0.2e1 + m(4) * (t435 ^ 2 + t436 ^ 2 + t437 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 * t562 ^ 2 + t514 ^ 2 + t515 ^ 2) / 0.2e1 + ((t510 * t540 - t511 * t525 + t512 * t526) * t545 + ((t484 * t540 - t486 * t525 + t488 * t526) * t541 + (t483 * t540 - t485 * t525 + t487 * t526) * t540) * qJD(3)) * t536 / 0.2e1 + t476 * ((t428 * t517 + t430 * t497 + t432 * t498) * t464 + (t427 * t517 + t429 * t497 + t431 * t498) * t463 + (t517 * t442 + t497 * t443 + t498 * t444) * t476) / 0.2e1 + t464 * ((t501 * t428 + t470 * t430 + t471 * t432) * t464 + (t427 * t501 + t429 * t470 + t431 * t471) * t463 + (t442 * t501 + t443 * t470 + t444 * t471) * t476) / 0.2e1 + t463 * ((t428 * t499 + t430 * t468 + t432 * t469) * t464 + (t499 * t427 + t468 * t429 + t469 * t431) * t463 + (t442 * t499 + t443 * t468 + t444 * t469) * t476) / 0.2e1 + t482 * ((t527 * t447 - t501 * t449 + t502 * t451) * t482 + (t446 * t527 - t448 * t501 + t450 * t502) * t481 + (t472 * t527 - t473 * t501 + t474 * t502) * t509) / 0.2e1 + t481 * ((t447 * t525 - t449 * t499 + t451 * t500) * t482 + (t525 * t446 - t499 * t448 + t500 * t450) * t481 + (t472 * t525 - t473 * t499 + t474 * t500) * t509) / 0.2e1 + t509 * ((t447 * t538 - t449 * t517 + t451 * t518) * t482 + (t446 * t538 - t448 * t517 + t450 * t518) * t481 + (t538 * t472 - t517 * t473 + t518 * t474) * t509) / 0.2e1 + t507 * ((t456 * t525 + t458 * t503 + t460 * t504) * t508 + (t525 * t455 + t503 * t457 + t504 * t459) * t507 + (t477 * t525 + t478 * t503 + t479 * t504) * t529) / 0.2e1 + t508 * ((t527 * t456 + t505 * t458 + t506 * t460) * t508 + (t455 * t527 + t457 * t505 + t459 * t506) * t507 + (t477 * t527 + t478 * t505 + t479 * t506) * t529) / 0.2e1 + t529 * ((t456 * t538 + t458 * t523 + t460 * t524) * t508 + (t455 * t538 + t457 * t523 + t459 * t524) * t507 + (t538 * t477 + t523 * t478 + t524 * t479) * t529) / 0.2e1 + t545 * ((t547 * t510 - t538 * t511 + t539 * t512) * t545 + ((t484 * t547 - t486 * t538 + t488 * t539) * t541 + (t483 * t547 - t485 * t538 + t487 * t539) * t540) * qJD(3)) / 0.2e1 + ((t510 * t541 - t511 * t527 + t512 * t528) * t545 + ((t484 * t541 - t486 * t527 + t488 * t528) * t541 + (t483 * t541 - t485 * t527 + t487 * t528) * t540) * qJD(3)) * t537 / 0.2e1 + (m(2) * (t554 ^ 2 + t555 ^ 2) + (Icges(3,5) * t562 + (Icges(3,1) * t601 + Icges(3,4) * t603) * t561) * t591 + t561 * t603 * (Icges(3,6) * t562 + (Icges(3,4) * t601 + Icges(3,2) * t603) * t561) + t562 * (Icges(3,3) * t562 + (Icges(3,5) * t601 + Icges(3,6) * t603) * t561) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
