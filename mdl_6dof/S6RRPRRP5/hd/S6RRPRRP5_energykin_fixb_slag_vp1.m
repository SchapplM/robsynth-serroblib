% Calculate kinetic energy for
% S6RRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:57:46
% EndTime: 2019-03-09 11:57:49
% DurationCPUTime: 3.33s
% Computational Cost: add. (3534->334), mult. (8860->503), div. (0->0), fcn. (11255->12), ass. (0->159)
t614 = Icges(6,1) + Icges(7,1);
t613 = Icges(6,4) + Icges(7,4);
t612 = Icges(6,5) + Icges(7,5);
t611 = Icges(6,2) + Icges(7,2);
t610 = Icges(6,6) + Icges(7,6);
t609 = Icges(6,3) + Icges(7,3);
t608 = rSges(7,3) + qJ(6);
t541 = sin(qJ(1));
t544 = cos(qJ(1));
t540 = sin(qJ(2));
t543 = cos(qJ(2));
t583 = sin(pkin(11));
t585 = cos(pkin(6));
t556 = t585 * t583;
t584 = cos(pkin(11));
t557 = t585 * t584;
t550 = -t540 * t556 + t543 * t557;
t552 = t540 * t584 + t543 * t583;
t498 = -t541 * t552 + t544 * t550;
t513 = t540 * t557 + t543 * t556;
t521 = -t540 * t583 + t543 * t584;
t499 = t513 * t544 + t521 * t541;
t536 = sin(pkin(6));
t578 = t536 * t544;
t451 = Icges(4,5) * t499 + Icges(4,6) * t498 - Icges(4,3) * t578;
t500 = -t541 * t550 - t544 * t552;
t501 = -t513 * t541 + t521 * t544;
t579 = t536 * t541;
t452 = Icges(4,5) * t501 + Icges(4,6) * t500 + Icges(4,3) * t579;
t563 = t544 * t585;
t515 = -t541 * t540 + t543 * t563;
t516 = t540 * t563 + t541 * t543;
t485 = Icges(3,5) * t516 + Icges(3,6) * t515 - Icges(3,3) * t578;
t565 = t541 * t585;
t517 = -t544 * t540 - t543 * t565;
t518 = -t540 * t565 + t544 * t543;
t486 = Icges(3,5) * t518 + Icges(3,6) * t517 + Icges(3,3) * t579;
t607 = t536 * ((t451 + t485) * t544 + (-t452 - t486) * t541);
t539 = sin(qJ(4));
t590 = cos(qJ(4));
t473 = t499 * t590 - t539 * t578;
t538 = sin(qJ(5));
t542 = cos(qJ(5));
t445 = -t473 * t538 - t498 * t542;
t582 = t498 * t538;
t446 = t473 * t542 - t582;
t569 = t536 * t590;
t472 = t499 * t539 + t544 * t569;
t606 = t445 * t610 + t446 * t612 + t472 * t609;
t475 = t501 * t590 + t539 * t579;
t447 = -t475 * t538 - t500 * t542;
t581 = t500 * t538;
t448 = t475 * t542 - t581;
t474 = t501 * t539 - t541 * t569;
t605 = t447 * t610 + t448 * t612 + t474 * t609;
t604 = t445 * t611 + t446 * t613 + t472 * t610;
t603 = t447 * t611 + t448 * t613 + t474 * t610;
t602 = t613 * t445 + t446 * t614 + t612 * t472;
t601 = t613 * t447 + t448 * t614 + t612 * t474;
t512 = t552 * t536;
t503 = t512 * t590 + t539 * t585;
t511 = t521 * t536;
t470 = -t503 * t538 - t511 * t542;
t580 = t511 * t538;
t471 = t503 * t542 - t580;
t502 = t512 * t539 - t585 * t590;
t600 = t470 * t610 + t471 * t612 + t502 * t609;
t599 = t470 * t611 + t471 * t613 + t502 * t610;
t598 = t613 * t470 + t471 * t614 + t612 * t502;
t597 = Icges(4,5) * t512 + Icges(4,6) * t511 + (Icges(3,5) * t540 + Icges(3,6) * t543) * t536 + (Icges(4,3) + Icges(3,3)) * t585;
t589 = pkin(2) * t540;
t588 = pkin(2) * t543;
t587 = pkin(5) * t542;
t577 = rSges(7,1) * t446 + rSges(7,2) * t445 - pkin(5) * t582 + t472 * t608 + t473 * t587;
t576 = rSges(7,1) * t448 + rSges(7,2) * t447 - pkin(5) * t581 + t474 * t608 + t475 * t587;
t575 = rSges(7,1) * t471 + rSges(7,2) * t470 - pkin(5) * t580 + t502 * t608 + t503 * t587;
t567 = -qJ(3) * t536 + t585 * t589;
t497 = -t541 * t567 + t544 * t588;
t519 = qJD(1) * (pkin(1) * t544 + pkin(8) * t579);
t530 = qJD(2) * t585 + qJD(1);
t574 = t530 * t497 + t519;
t572 = qJD(2) * t536;
t529 = t541 * t572;
t476 = -qJD(4) * t500 + t529;
t573 = qJD(1) * (pkin(1) * t541 - pkin(8) * t578);
t571 = qJD(3) * t544;
t496 = t541 * t588 + t544 * t567;
t568 = t544 * t572;
t570 = qJD(3) * t585 + t496 * t529 + t497 * t568;
t504 = -qJD(4) * t511 + t530;
t522 = qJ(3) * t585 + t536 * t589;
t562 = qJD(2) * (-t512 * rSges(4,1) - t511 * rSges(4,2) - rSges(4,3) * t585 - t522);
t561 = qJD(2) * (-pkin(3) * t512 + pkin(9) * t511 - t522);
t560 = qJD(3) * t579 - t573;
t477 = -qJD(4) * t498 - t568;
t466 = pkin(3) * t499 - pkin(9) * t498;
t467 = pkin(3) * t501 - pkin(9) * t500;
t555 = t466 * t529 + t467 * t568 + t570;
t440 = pkin(4) * t473 + pkin(10) * t472;
t441 = pkin(4) * t475 + pkin(10) * t474;
t551 = t476 * t440 - t441 * t477 + t555;
t549 = t530 * t467 + (t541 * t561 - t571) * t536 + t574;
t548 = (-t466 - t496) * t530 + t561 * t578 + t560;
t468 = pkin(4) * t503 + pkin(10) * t502;
t547 = t504 * t441 - t468 * t476 + t549;
t546 = -t440 * t504 + t477 * t468 + t548;
t525 = rSges(2,1) * t544 - rSges(2,2) * t541;
t524 = rSges(2,1) * t541 + rSges(2,2) * t544;
t510 = t585 * rSges(3,3) + (rSges(3,1) * t540 + rSges(3,2) * t543) * t536;
t509 = Icges(3,5) * t585 + (Icges(3,1) * t540 + Icges(3,4) * t543) * t536;
t508 = Icges(3,6) * t585 + (Icges(3,4) * t540 + Icges(3,2) * t543) * t536;
t492 = rSges(3,1) * t518 + rSges(3,2) * t517 + rSges(3,3) * t579;
t491 = rSges(3,1) * t516 + rSges(3,2) * t515 - rSges(3,3) * t578;
t490 = Icges(3,1) * t518 + Icges(3,4) * t517 + Icges(3,5) * t579;
t489 = Icges(3,1) * t516 + Icges(3,4) * t515 - Icges(3,5) * t578;
t488 = Icges(3,4) * t518 + Icges(3,2) * t517 + Icges(3,6) * t579;
t487 = Icges(3,4) * t516 + Icges(3,2) * t515 - Icges(3,6) * t578;
t483 = Icges(4,1) * t512 + Icges(4,4) * t511 + Icges(4,5) * t585;
t482 = Icges(4,4) * t512 + Icges(4,2) * t511 + Icges(4,6) * t585;
t469 = qJD(5) * t502 + t504;
t465 = rSges(5,1) * t503 - rSges(5,2) * t502 - rSges(5,3) * t511;
t464 = Icges(5,1) * t503 - Icges(5,4) * t502 - Icges(5,5) * t511;
t463 = Icges(5,4) * t503 - Icges(5,2) * t502 - Icges(5,6) * t511;
t462 = Icges(5,5) * t503 - Icges(5,6) * t502 - Icges(5,3) * t511;
t459 = rSges(4,1) * t501 + rSges(4,2) * t500 + rSges(4,3) * t579;
t458 = rSges(4,1) * t499 + rSges(4,2) * t498 - rSges(4,3) * t578;
t456 = Icges(4,1) * t501 + Icges(4,4) * t500 + Icges(4,5) * t579;
t455 = Icges(4,1) * t499 + Icges(4,4) * t498 - Icges(4,5) * t578;
t454 = Icges(4,4) * t501 + Icges(4,2) * t500 + Icges(4,6) * t579;
t453 = Icges(4,4) * t499 + Icges(4,2) * t498 - Icges(4,6) * t578;
t450 = t492 * t530 - t510 * t529 + t519;
t449 = -t491 * t530 - t510 * t568 - t573;
t444 = (t491 * t541 + t492 * t544) * t572;
t443 = qJD(5) * t472 + t477;
t442 = qJD(5) * t474 + t476;
t437 = rSges(6,1) * t471 + rSges(6,2) * t470 + rSges(6,3) * t502;
t429 = rSges(5,1) * t475 - rSges(5,2) * t474 - rSges(5,3) * t500;
t428 = rSges(5,1) * t473 - rSges(5,2) * t472 - rSges(5,3) * t498;
t427 = Icges(5,1) * t475 - Icges(5,4) * t474 - Icges(5,5) * t500;
t426 = Icges(5,1) * t473 - Icges(5,4) * t472 - Icges(5,5) * t498;
t425 = Icges(5,4) * t475 - Icges(5,2) * t474 - Icges(5,6) * t500;
t424 = Icges(5,4) * t473 - Icges(5,2) * t472 - Icges(5,6) * t498;
t423 = Icges(5,5) * t475 - Icges(5,6) * t474 - Icges(5,3) * t500;
t422 = Icges(5,5) * t473 - Icges(5,6) * t472 - Icges(5,3) * t498;
t419 = t459 * t530 + (t541 * t562 - t571) * t536 + t574;
t418 = (-t458 - t496) * t530 + t562 * t578 + t560;
t417 = rSges(6,1) * t448 + rSges(6,2) * t447 + rSges(6,3) * t474;
t415 = rSges(6,1) * t446 + rSges(6,2) * t445 + rSges(6,3) * t472;
t399 = (t458 * t541 + t459 * t544) * t572 + t570;
t398 = t429 * t504 - t465 * t476 + t549;
t397 = -t428 * t504 + t465 * t477 + t548;
t396 = t428 * t476 - t429 * t477 + t555;
t395 = t417 * t469 - t437 * t442 + t547;
t394 = -t415 * t469 + t437 * t443 + t546;
t393 = t415 * t442 - t417 * t443 + t551;
t392 = qJD(6) * t472 - t442 * t575 + t469 * t576 + t547;
t391 = qJD(6) * t474 + t443 * t575 - t469 * t577 + t546;
t390 = qJD(6) * t502 + t442 * t577 - t443 * t576 + t551;
t1 = m(4) * (t399 ^ 2 + t418 ^ 2 + t419 ^ 2) / 0.2e1 + t504 * ((-t423 * t511 - t425 * t502 + t427 * t503) * t476 + (-t422 * t511 - t424 * t502 + t426 * t503) * t477 + (-t462 * t511 - t463 * t502 + t464 * t503) * t504) / 0.2e1 + t476 * ((-t423 * t500 - t425 * t474 + t427 * t475) * t476 + (-t422 * t500 - t424 * t474 + t426 * t475) * t477 + (-t462 * t500 - t463 * t474 + t464 * t475) * t504) / 0.2e1 + t477 * ((-t423 * t498 - t425 * t472 + t427 * t473) * t476 + (-t422 * t498 - t424 * t472 + t426 * t473) * t477 + (-t462 * t498 - t463 * t472 + t464 * t473) * t504) / 0.2e1 + m(3) * (t444 ^ 2 + t449 ^ 2 + t450 ^ 2) / 0.2e1 + m(7) * (t390 ^ 2 + t391 ^ 2 + t392 ^ 2) / 0.2e1 + m(5) * (t396 ^ 2 + t397 ^ 2 + t398 ^ 2) / 0.2e1 + m(6) * (t393 ^ 2 + t394 ^ 2 + t395 ^ 2) / 0.2e1 + ((t447 * t599 + t448 * t598 + t474 * t600) * t469 + (t604 * t447 + t602 * t448 + t474 * t606) * t443 + (t603 * t447 + t601 * t448 + t605 * t474) * t442) * t442 / 0.2e1 + ((t445 * t599 + t446 * t598 + t472 * t600) * t469 + (t604 * t445 + t602 * t446 + t606 * t472) * t443 + (t445 * t603 + t446 * t601 + t472 * t605) * t442) * t443 / 0.2e1 + ((t599 * t470 + t598 * t471 + t600 * t502) * t469 + (t604 * t470 + t602 * t471 + t502 * t606) * t443 + (t470 * t603 + t471 * t601 + t502 * t605) * t442) * t469 / 0.2e1 + ((t585 * t486 + (t488 * t543 + t490 * t540) * t536) * t529 - (t585 * t485 + (t487 * t543 + t489 * t540) * t536) * t568 + ((t452 * t585 + t511 * t454 + t512 * t456) * t541 - (t451 * t585 + t511 * t453 + t512 * t455) * t544) * t572 + ((t508 * t543 + t509 * t540) * t536 + t511 * t482 + t512 * t483 + t597 * t585) * t530) * t530 / 0.2e1 + (m(2) * (t524 ^ 2 + t525 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 - (((-t498 * t453 - t499 * t455 - t515 * t487 - t516 * t489 + t607) * t544 + (t454 * t498 + t456 * t499 + t488 * t515 + t490 * t516) * t541) * t572 + (t482 * t498 + t483 * t499 + t508 * t515 + t509 * t516 - t578 * t597) * t530) * t568 / 0.2e1 + (((-t453 * t500 - t455 * t501 - t487 * t517 - t489 * t518) * t544 + (t500 * t454 + t501 * t456 + t517 * t488 + t518 * t490 - t607) * t541) * t572 + (t482 * t500 + t483 * t501 + t508 * t517 + t509 * t518 + t579 * t597) * t530) * t529 / 0.2e1;
T  = t1;
