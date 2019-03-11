% Calculate kinetic energy for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:55:13
% EndTime: 2019-03-08 19:55:16
% DurationCPUTime: 3.23s
% Computational Cost: add. (3465->322), mult. (8808->488), div. (0->0), fcn. (11221->12), ass. (0->154)
t600 = Icges(6,1) + Icges(7,1);
t599 = Icges(6,4) + Icges(7,4);
t598 = Icges(6,5) + Icges(7,5);
t597 = Icges(6,2) + Icges(7,2);
t596 = Icges(6,6) + Icges(7,6);
t595 = Icges(6,3) + Icges(7,3);
t526 = sin(pkin(10));
t528 = cos(pkin(10));
t533 = sin(qJ(2));
t535 = cos(qJ(2));
t569 = sin(pkin(11));
t570 = cos(pkin(11));
t515 = -t533 * t569 + t535 * t570;
t529 = cos(pkin(6));
t541 = t529 * t515;
t542 = t533 * t570 + t535 * t569;
t493 = -t526 * t542 + t528 * t541;
t508 = t542 * t529;
t494 = t508 * t528 + t515 * t526;
t527 = sin(pkin(6));
t564 = t527 * t528;
t444 = Icges(4,5) * t494 + Icges(4,6) * t493 - Icges(4,3) * t564;
t495 = -t526 * t541 - t528 * t542;
t496 = -t508 * t526 + t515 * t528;
t565 = t526 * t527;
t445 = Icges(4,5) * t496 + Icges(4,6) * t495 + Icges(4,3) * t565;
t561 = t529 * t535;
t510 = -t526 * t533 + t528 * t561;
t562 = t529 * t533;
t511 = t526 * t535 + t528 * t562;
t480 = Icges(3,5) * t511 + Icges(3,6) * t510 - Icges(3,3) * t564;
t512 = -t526 * t561 - t528 * t533;
t513 = -t526 * t562 + t528 * t535;
t481 = Icges(3,5) * t513 + Icges(3,6) * t512 + Icges(3,3) * t565;
t506 = t515 * t527;
t507 = t542 * t527;
t579 = (Icges(4,5) * t507 + Icges(4,6) * t506 + (Icges(3,5) * t533 + Icges(3,6) * t535) * t527 + (Icges(4,3) + Icges(3,3)) * t529) * t529;
t594 = -t579 + (t444 + t480) * t564 - (t481 + t445) * t565;
t593 = rSges(7,3) + qJ(6);
t532 = sin(qJ(4));
t563 = t527 * t532;
t574 = cos(qJ(4));
t468 = t494 * t574 - t528 * t563;
t531 = sin(qJ(5));
t534 = cos(qJ(5));
t440 = -t468 * t531 - t493 * t534;
t568 = t493 * t531;
t441 = t468 * t534 - t568;
t553 = t527 * t574;
t467 = t494 * t532 + t528 * t553;
t591 = t440 * t596 + t441 * t598 + t467 * t595;
t470 = t496 * t574 + t526 * t563;
t442 = -t470 * t531 - t495 * t534;
t567 = t495 * t531;
t443 = t470 * t534 - t567;
t469 = t496 * t532 - t526 * t553;
t590 = t442 * t596 + t443 * t598 + t469 * t595;
t589 = t440 * t597 + t441 * t599 + t467 * t596;
t588 = t442 * t597 + t443 * t599 + t469 * t596;
t587 = t599 * t440 + t441 * t600 + t598 * t467;
t586 = t599 * t442 + t443 * t600 + t598 * t469;
t498 = t507 * t574 + t529 * t532;
t465 = -t498 * t531 - t506 * t534;
t566 = t506 * t531;
t466 = t498 * t534 - t566;
t497 = t507 * t532 - t529 * t574;
t585 = t465 * t596 + t466 * t598 + t497 * t595;
t584 = t465 * t597 + t466 * t599 + t497 * t596;
t583 = t599 * t465 + t466 * t600 + t598 * t497;
t578 = qJD(2) ^ 2;
t573 = pkin(2) * t535;
t572 = pkin(5) * t534;
t560 = rSges(7,1) * t441 + rSges(7,2) * t440 - pkin(5) * t568 + t467 * t593 + t468 * t572;
t559 = rSges(7,1) * t443 + rSges(7,2) * t442 - pkin(5) * t567 + t469 * t593 + t470 * t572;
t558 = rSges(7,1) * t466 + rSges(7,2) * t465 - pkin(5) * t566 + t497 * t593 + t498 * t572;
t516 = pkin(2) * t527 * t533 + qJ(3) * t529;
t557 = -pkin(3) * t507 + pkin(8) * t506 - t516;
t556 = qJD(2) * t527;
t520 = t526 * t556;
t471 = -qJD(4) * t495 + t520;
t525 = qJD(2) * t529;
t499 = -qJD(4) * t506 + t525;
t555 = qJD(3) * t528;
t552 = t528 * t556;
t550 = pkin(2) * t562 - qJ(3) * t527;
t547 = (-rSges(4,1) * t507 - rSges(4,2) * t506 - rSges(4,3) * t529 - t516) * t527;
t488 = t526 * t573 + t528 * t550;
t489 = -t526 * t550 + t528 * t573;
t546 = qJD(3) * t529 + t488 * t520 + t489 * t552 + qJD(1);
t472 = -qJD(4) * t493 - t552;
t455 = pkin(3) * t494 - pkin(8) * t493;
t456 = pkin(3) * t496 - pkin(8) * t495;
t543 = t455 * t520 + t456 * t552 + t546;
t435 = pkin(4) * t468 + pkin(9) * t467;
t436 = pkin(4) * t470 + pkin(9) * t469;
t540 = t471 * t435 - t436 * t472 + t543;
t475 = t489 * t525;
t539 = t456 * t525 + t475 + (qJD(2) * t526 * t557 - t555) * t527;
t519 = qJD(3) * t565;
t538 = t519 + ((-t455 - t488) * t529 + t557 * t564) * qJD(2);
t463 = pkin(4) * t498 + pkin(9) * t497;
t537 = t499 * t436 - t463 * t471 + t539;
t536 = -t435 * t499 + t472 * t463 + t538;
t505 = t529 * rSges(3,3) + (rSges(3,1) * t533 + rSges(3,2) * t535) * t527;
t504 = Icges(3,5) * t529 + (Icges(3,1) * t533 + Icges(3,4) * t535) * t527;
t503 = Icges(3,6) * t529 + (Icges(3,4) * t533 + Icges(3,2) * t535) * t527;
t487 = rSges(3,1) * t513 + rSges(3,2) * t512 + rSges(3,3) * t565;
t486 = rSges(3,1) * t511 + rSges(3,2) * t510 - rSges(3,3) * t564;
t485 = Icges(3,1) * t513 + Icges(3,4) * t512 + Icges(3,5) * t565;
t484 = Icges(3,1) * t511 + Icges(3,4) * t510 - Icges(3,5) * t564;
t483 = Icges(3,4) * t513 + Icges(3,2) * t512 + Icges(3,6) * t565;
t482 = Icges(3,4) * t511 + Icges(3,2) * t510 - Icges(3,6) * t564;
t478 = Icges(4,1) * t507 + Icges(4,4) * t506 + Icges(4,5) * t529;
t477 = Icges(4,4) * t507 + Icges(4,2) * t506 + Icges(4,6) * t529;
t464 = qJD(5) * t497 + t499;
t462 = (-t486 * t529 - t505 * t564) * qJD(2);
t461 = (t487 * t529 - t505 * t565) * qJD(2);
t460 = rSges(5,1) * t498 - rSges(5,2) * t497 - rSges(5,3) * t506;
t459 = Icges(5,1) * t498 - Icges(5,4) * t497 - Icges(5,5) * t506;
t458 = Icges(5,4) * t498 - Icges(5,2) * t497 - Icges(5,6) * t506;
t457 = Icges(5,5) * t498 - Icges(5,6) * t497 - Icges(5,3) * t506;
t451 = rSges(4,1) * t496 + rSges(4,2) * t495 + rSges(4,3) * t565;
t450 = rSges(4,1) * t494 + rSges(4,2) * t493 - rSges(4,3) * t564;
t449 = Icges(4,1) * t496 + Icges(4,4) * t495 + Icges(4,5) * t565;
t448 = Icges(4,1) * t494 + Icges(4,4) * t493 - Icges(4,5) * t564;
t447 = Icges(4,4) * t496 + Icges(4,2) * t495 + Icges(4,6) * t565;
t446 = Icges(4,4) * t494 + Icges(4,2) * t493 - Icges(4,6) * t564;
t439 = qJD(1) + (t486 * t526 + t487 * t528) * t556;
t438 = qJD(5) * t467 + t472;
t437 = qJD(5) * t469 + t471;
t432 = rSges(6,1) * t466 + rSges(6,2) * t465 + rSges(6,3) * t497;
t424 = rSges(5,1) * t470 - rSges(5,2) * t469 - rSges(5,3) * t495;
t423 = rSges(5,1) * t468 - rSges(5,2) * t467 - rSges(5,3) * t493;
t422 = Icges(5,1) * t470 - Icges(5,4) * t469 - Icges(5,5) * t495;
t421 = Icges(5,1) * t468 - Icges(5,4) * t467 - Icges(5,5) * t493;
t420 = Icges(5,4) * t470 - Icges(5,2) * t469 - Icges(5,6) * t495;
t419 = Icges(5,4) * t468 - Icges(5,2) * t467 - Icges(5,6) * t493;
t418 = Icges(5,5) * t470 - Icges(5,6) * t469 - Icges(5,3) * t495;
t417 = Icges(5,5) * t468 - Icges(5,6) * t467 - Icges(5,3) * t493;
t414 = t519 + ((-t450 - t488) * t529 + t528 * t547) * qJD(2);
t413 = -t527 * t555 + t475 + (t451 * t529 + t526 * t547) * qJD(2);
t412 = rSges(6,1) * t443 + rSges(6,2) * t442 + rSges(6,3) * t469;
t410 = rSges(6,1) * t441 + rSges(6,2) * t440 + rSges(6,3) * t467;
t394 = (t450 * t526 + t451 * t528) * t556 + t546;
t393 = -t423 * t499 + t460 * t472 + t538;
t392 = t424 * t499 - t460 * t471 + t539;
t391 = t423 * t471 - t424 * t472 + t543;
t390 = -t410 * t464 + t432 * t438 + t536;
t389 = t412 * t464 - t432 * t437 + t537;
t388 = t410 * t437 - t412 * t438 + t540;
t387 = qJD(6) * t469 + t438 * t558 - t464 * t560 + t536;
t386 = qJD(6) * t467 - t437 * t558 + t464 * t559 + t537;
t385 = qJD(6) * t497 + t437 * t560 - t438 * t559 + t540;
t1 = t471 * ((-t495 * t418 - t469 * t420 + t470 * t422) * t471 + (-t417 * t495 - t419 * t469 + t421 * t470) * t472 + (-t457 * t495 - t458 * t469 + t459 * t470) * t499) / 0.2e1 + t472 * ((-t418 * t493 - t420 * t467 + t422 * t468) * t471 + (-t493 * t417 - t467 * t419 + t468 * t421) * t472 + (-t457 * t493 - t458 * t467 + t459 * t468) * t499) / 0.2e1 + t499 * ((-t418 * t506 - t420 * t497 + t422 * t498) * t471 + (-t417 * t506 - t419 * t497 + t421 * t498) * t472 + (-t506 * t457 - t497 * t458 + t498 * t459) * t499) / 0.2e1 + m(5) * (t391 ^ 2 + t392 ^ 2 + t393 ^ 2) / 0.2e1 + m(3) * (t439 ^ 2 + t461 ^ 2 + t462 ^ 2) / 0.2e1 + m(4) * (t394 ^ 2 + t413 ^ 2 + t414 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(7) * (t385 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + m(6) * (t388 ^ 2 + t389 ^ 2 + t390 ^ 2) / 0.2e1 + ((t442 * t584 + t443 * t583 + t469 * t585) * t464 + (t589 * t442 + t587 * t443 + t469 * t591) * t438 + (t588 * t442 + t586 * t443 + t590 * t469) * t437) * t437 / 0.2e1 + ((t440 * t584 + t441 * t583 + t467 * t585) * t464 + (t589 * t440 + t587 * t441 + t591 * t467) * t438 + (t440 * t588 + t441 * t586 + t467 * t590) * t437) * t438 / 0.2e1 + ((t584 * t465 + t583 * t466 + t585 * t497) * t464 + (t589 * t465 + t587 * t466 + t497 * t591) * t438 + (t465 * t588 + t466 * t586 + t497 * t590) * t437) * t464 / 0.2e1 - ((t447 * t493 + t449 * t494 + t483 * t510 + t485 * t511) * t565 + (t477 * t493 + t478 * t494 + t503 * t510 + t504 * t511) * t529 + (-t446 * t493 - t448 * t494 - t482 * t510 - t484 * t511 + t594) * t564) * t578 * t564 / 0.2e1 + (((t506 * t477 + t507 * t478 + t579) * t529 + ((t445 * t529 + t447 * t506 + t449 * t507) * t526 - (t444 * t529 + t446 * t506 + t448 * t507) * t528 + (-t480 * t528 + t481 * t526 + t503 * t535 + t504 * t533) * t529 + ((t483 * t535 + t485 * t533) * t526 - (t482 * t535 + t484 * t533) * t528) * t527) * t527) * t529 + ((-t446 * t495 - t448 * t496 - t482 * t512 - t484 * t513) * t564 + (t477 * t495 + t478 * t496 + t503 * t512 + t504 * t513) * t529 + (t447 * t495 + t449 * t496 + t483 * t512 + t485 * t513 - t594) * t565) * t565) * t578 / 0.2e1;
T  = t1;
