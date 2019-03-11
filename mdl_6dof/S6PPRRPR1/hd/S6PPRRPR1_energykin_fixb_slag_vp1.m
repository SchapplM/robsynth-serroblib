% Calculate kinetic energy for
% S6PPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-03-08 18:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_energykin_fixb_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRPR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRPR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:45:17
% EndTime: 2019-03-08 18:45:19
% DurationCPUTime: 2.24s
% Computational Cost: add. (5008->316), mult. (13361->476), div. (0->0), fcn. (17435->16), ass. (0->148)
t595 = Icges(5,2) + Icges(6,3);
t549 = sin(pkin(11));
t552 = cos(pkin(11));
t581 = sin(pkin(12));
t585 = cos(pkin(6));
t569 = t585 * t581;
t583 = cos(pkin(12));
t535 = t549 * t583 + t552 * t569;
t555 = sin(qJ(3));
t571 = t585 * t583;
t561 = t549 * t581 - t552 * t571;
t584 = cos(pkin(7));
t558 = t561 * t584;
t550 = sin(pkin(6));
t582 = sin(pkin(7));
t572 = t550 * t582;
t588 = cos(qJ(3));
t515 = t535 * t588 + (-t552 * t572 - t558) * t555;
t573 = t550 * t584;
t528 = -t552 * t573 + t561 * t582;
t554 = sin(qJ(4));
t587 = cos(qJ(4));
t501 = t515 * t587 + t528 * t554;
t567 = t588 * t572;
t514 = t535 * t555 + t552 * t567 + t558 * t588;
t548 = sin(pkin(13));
t551 = cos(pkin(13));
t471 = -t501 * t548 + t514 * t551;
t580 = t514 * t548;
t472 = t501 * t551 + t580;
t500 = t515 * t554 - t528 * t587;
t594 = -Icges(5,4) * t501 + Icges(6,5) * t472 - Icges(5,6) * t514 + Icges(6,6) * t471 + t595 * t500;
t536 = -t549 * t569 + t552 * t583;
t560 = t549 * t571 + t552 * t581;
t557 = t560 * t584;
t517 = t536 * t588 + (t549 * t572 - t557) * t555;
t529 = t549 * t573 + t560 * t582;
t503 = t517 * t587 + t529 * t554;
t516 = t536 * t555 - t549 * t567 + t557 * t588;
t473 = -t503 * t548 + t516 * t551;
t579 = t516 * t548;
t474 = t503 * t551 + t579;
t502 = t517 * t554 - t529 * t587;
t593 = -Icges(5,4) * t503 + Icges(6,5) * t474 - Icges(5,6) * t516 + Icges(6,6) * t473 + t595 * t502;
t568 = t584 * t583;
t570 = t585 * t582;
t527 = t555 * t570 + (t555 * t568 + t581 * t588) * t550;
t534 = -t572 * t583 + t584 * t585;
t519 = t527 * t587 + t534 * t554;
t526 = -t570 * t588 + (t555 * t581 - t568 * t588) * t550;
t498 = -t519 * t548 + t526 * t551;
t578 = t526 * t548;
t499 = t519 * t551 + t578;
t518 = t527 * t554 - t534 * t587;
t592 = -Icges(5,4) * t519 + Icges(6,5) * t499 - Icges(5,6) * t526 + Icges(6,6) * t498 + t595 * t518;
t586 = pkin(5) * t551;
t524 = qJD(3) * t528;
t504 = qJD(4) * t514 + t524;
t525 = qJD(3) * t529;
t505 = qJD(4) * t516 + t525;
t533 = qJD(3) * t534;
t520 = qJD(4) * t526 + t533;
t576 = qJD(2) * t550;
t541 = qJD(2) * t585 + qJD(1);
t574 = t552 * t576;
t489 = pkin(3) * t515 + pkin(9) * t514;
t510 = pkin(3) * t527 + pkin(9) * t526;
t540 = t549 * t576;
t566 = -t489 * t533 + t510 * t524 + t540;
t490 = pkin(3) * t517 + pkin(9) * t516;
t565 = t489 * t525 - t490 * t524 + t541;
t491 = t519 * pkin(4) + t518 * qJ(5);
t564 = qJD(5) * t502 + t504 * t491 + t566;
t563 = t490 * t533 - t510 * t525 - t574;
t463 = t501 * pkin(4) + t500 * qJ(5);
t562 = qJD(5) * t518 + t505 * t463 + t565;
t464 = t503 * pkin(4) + t502 * qJ(5);
t559 = qJD(5) * t500 + t520 * t464 + t563;
t547 = pkin(13) + qJ(6);
t545 = cos(t547);
t544 = sin(t547);
t509 = rSges(4,1) * t527 - rSges(4,2) * t526 + rSges(4,3) * t534;
t508 = Icges(4,1) * t527 - Icges(4,4) * t526 + Icges(4,5) * t534;
t507 = Icges(4,4) * t527 - Icges(4,2) * t526 + Icges(4,6) * t534;
t506 = Icges(4,5) * t527 - Icges(4,6) * t526 + Icges(4,3) * t534;
t495 = t519 * t545 + t526 * t544;
t494 = -t519 * t544 + t526 * t545;
t492 = qJD(6) * t518 + t520;
t487 = rSges(5,1) * t519 - rSges(5,2) * t518 + rSges(5,3) * t526;
t486 = Icges(5,1) * t519 - Icges(5,4) * t518 + Icges(5,5) * t526;
t484 = Icges(5,5) * t519 - Icges(5,6) * t518 + Icges(5,3) * t526;
t482 = rSges(4,1) * t517 - rSges(4,2) * t516 + rSges(4,3) * t529;
t481 = rSges(4,1) * t515 - rSges(4,2) * t514 + rSges(4,3) * t528;
t480 = Icges(4,1) * t517 - Icges(4,4) * t516 + Icges(4,5) * t529;
t479 = Icges(4,1) * t515 - Icges(4,4) * t514 + Icges(4,5) * t528;
t478 = Icges(4,4) * t517 - Icges(4,2) * t516 + Icges(4,6) * t529;
t477 = Icges(4,4) * t515 - Icges(4,2) * t514 + Icges(4,6) * t528;
t476 = Icges(4,5) * t517 - Icges(4,6) * t516 + Icges(4,3) * t529;
t475 = Icges(4,5) * t515 - Icges(4,6) * t514 + Icges(4,3) * t528;
t470 = t503 * t545 + t516 * t544;
t469 = -t503 * t544 + t516 * t545;
t468 = t501 * t545 + t514 * t544;
t467 = -t501 * t544 + t514 * t545;
t466 = qJD(6) * t502 + t505;
t465 = qJD(6) * t500 + t504;
t460 = rSges(5,1) * t503 - rSges(5,2) * t502 + rSges(5,3) * t516;
t459 = rSges(5,1) * t501 - rSges(5,2) * t500 + rSges(5,3) * t514;
t458 = rSges(6,1) * t499 + rSges(6,2) * t498 + rSges(6,3) * t518;
t457 = Icges(5,1) * t503 - Icges(5,4) * t502 + Icges(5,5) * t516;
t456 = Icges(5,1) * t501 - Icges(5,4) * t500 + Icges(5,5) * t514;
t453 = Icges(5,5) * t503 - Icges(5,6) * t502 + Icges(5,3) * t516;
t452 = Icges(5,5) * t501 - Icges(5,6) * t500 + Icges(5,3) * t514;
t451 = Icges(6,1) * t499 + Icges(6,4) * t498 + Icges(6,5) * t518;
t450 = Icges(6,4) * t499 + Icges(6,2) * t498 + Icges(6,6) * t518;
t448 = rSges(7,1) * t495 + rSges(7,2) * t494 + rSges(7,3) * t518;
t447 = Icges(7,1) * t495 + Icges(7,4) * t494 + Icges(7,5) * t518;
t446 = Icges(7,4) * t495 + Icges(7,2) * t494 + Icges(7,6) * t518;
t445 = Icges(7,5) * t495 + Icges(7,6) * t494 + Icges(7,3) * t518;
t444 = pkin(5) * t578 + pkin(10) * t518 + t519 * t586;
t442 = -t574 + (t482 * t534 - t509 * t529) * qJD(3);
t441 = t540 + (-t481 * t534 + t509 * t528) * qJD(3);
t440 = (t481 * t529 - t482 * t528) * qJD(3) + t541;
t439 = rSges(6,1) * t474 + rSges(6,2) * t473 + rSges(6,3) * t502;
t438 = rSges(6,1) * t472 + rSges(6,2) * t471 + rSges(6,3) * t500;
t437 = Icges(6,1) * t474 + Icges(6,4) * t473 + Icges(6,5) * t502;
t436 = Icges(6,1) * t472 + Icges(6,4) * t471 + Icges(6,5) * t500;
t435 = Icges(6,4) * t474 + Icges(6,2) * t473 + Icges(6,6) * t502;
t434 = Icges(6,4) * t472 + Icges(6,2) * t471 + Icges(6,6) * t500;
t431 = rSges(7,1) * t470 + rSges(7,2) * t469 + rSges(7,3) * t502;
t430 = rSges(7,1) * t468 + rSges(7,2) * t467 + rSges(7,3) * t500;
t429 = Icges(7,1) * t470 + Icges(7,4) * t469 + Icges(7,5) * t502;
t428 = Icges(7,1) * t468 + Icges(7,4) * t467 + Icges(7,5) * t500;
t427 = Icges(7,4) * t470 + Icges(7,2) * t469 + Icges(7,6) * t502;
t426 = Icges(7,4) * t468 + Icges(7,2) * t467 + Icges(7,6) * t500;
t425 = Icges(7,5) * t470 + Icges(7,6) * t469 + Icges(7,3) * t502;
t424 = Icges(7,5) * t468 + Icges(7,6) * t467 + Icges(7,3) * t500;
t423 = pkin(5) * t579 + pkin(10) * t502 + t503 * t586;
t422 = pkin(5) * t580 + pkin(10) * t500 + t501 * t586;
t421 = t460 * t520 - t487 * t505 + t563;
t420 = -t459 * t520 + t487 * t504 + t566;
t419 = t459 * t505 - t460 * t504 + t565;
t418 = t520 * t439 + (-t458 - t491) * t505 + t559;
t417 = t504 * t458 + (-t438 - t463) * t520 + t564;
t416 = t505 * t438 + (-t439 - t464) * t504 + t562;
t415 = t520 * t423 + t492 * t431 - t466 * t448 + (-t444 - t491) * t505 + t559;
t414 = -t492 * t430 + t504 * t444 + t465 * t448 + (-t422 - t463) * t520 + t564;
t413 = t505 * t422 + t466 * t430 - t465 * t431 + (-t423 - t464) * t504 + t562;
t1 = m(6) * (t416 ^ 2 + t417 ^ 2 + t418 ^ 2) / 0.2e1 + m(7) * (t413 ^ 2 + t414 ^ 2 + t415 ^ 2) / 0.2e1 + m(5) * (t419 ^ 2 + t420 ^ 2 + t421 ^ 2) / 0.2e1 + m(4) * (t440 ^ 2 + t441 ^ 2 + t442 ^ 2) / 0.2e1 + m(3) * (t541 ^ 2 + (t549 ^ 2 + t552 ^ 2) * qJD(2) ^ 2 * t550 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t466 * ((t502 * t425 + t469 * t427 + t470 * t429) * t466 + (t424 * t502 + t426 * t469 + t428 * t470) * t465 + (t445 * t502 + t446 * t469 + t447 * t470) * t492) / 0.2e1 + t465 * ((t425 * t500 + t427 * t467 + t429 * t468) * t466 + (t500 * t424 + t467 * t426 + t468 * t428) * t465 + (t445 * t500 + t446 * t467 + t447 * t468) * t492) / 0.2e1 + t492 * ((t425 * t518 + t427 * t494 + t429 * t495) * t466 + (t424 * t518 + t426 * t494 + t428 * t495) * t465 + (t518 * t445 + t494 * t446 + t495 * t447) * t492) / 0.2e1 + ((t450 * t471 + t451 * t472 + t484 * t514 + t486 * t501 + t500 * t592) * t520 + (t435 * t471 + t437 * t472 + t453 * t514 + t457 * t501 + t500 * t593) * t505 + (t471 * t434 + t472 * t436 + t514 * t452 + t501 * t456 + t500 * t594) * t504) * t504 / 0.2e1 + ((t450 * t473 + t451 * t474 + t484 * t516 + t486 * t503 + t502 * t592) * t520 + (t473 * t435 + t474 * t437 + t516 * t453 + t503 * t457 + t502 * t593) * t505 + (t434 * t473 + t436 * t474 + t452 * t516 + t456 * t503 + t502 * t594) * t504) * t505 / 0.2e1 + ((t498 * t450 + t499 * t451 + t526 * t484 + t519 * t486 + t592 * t518) * t520 + (t435 * t498 + t437 * t499 + t453 * t526 + t457 * t519 + t518 * t593) * t505 + (t434 * t498 + t436 * t499 + t452 * t526 + t456 * t519 + t518 * t594) * t504) * t520 / 0.2e1 + (t529 * ((t529 * t476 - t516 * t478 + t517 * t480) * t529 + (t475 * t529 - t477 * t516 + t479 * t517) * t528 + (t506 * t529 - t507 * t516 + t508 * t517) * t534) + t528 * ((t476 * t528 - t478 * t514 + t480 * t515) * t529 + (t528 * t475 - t514 * t477 + t515 * t479) * t528 + (t506 * t528 - t507 * t514 + t508 * t515) * t534) + t534 * ((t476 * t534 - t478 * t526 + t480 * t527) * t529 + (t475 * t534 - t477 * t526 + t479 * t527) * t528 + (t534 * t506 - t526 * t507 + t527 * t508) * t534)) * qJD(3) ^ 2 / 0.2e1;
T  = t1;
