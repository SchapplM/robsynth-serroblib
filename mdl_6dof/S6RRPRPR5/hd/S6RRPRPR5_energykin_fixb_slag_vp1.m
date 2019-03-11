% Calculate kinetic energy for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:28:57
% EndTime: 2019-03-09 10:29:00
% DurationCPUTime: 3.51s
% Computational Cost: add. (3611->380), mult. (8643->563), div. (0->0), fcn. (10968->14), ass. (0->170)
t603 = Icges(5,2) + Icges(6,3);
t547 = sin(qJ(1));
t548 = cos(qJ(1));
t546 = sin(qJ(2));
t585 = sin(pkin(11));
t587 = cos(pkin(6));
t560 = t587 * t585;
t586 = cos(pkin(11));
t561 = t587 * t586;
t591 = cos(qJ(2));
t553 = -t546 * t560 + t561 * t591;
t555 = t546 * t586 + t585 * t591;
t499 = -t547 * t555 + t548 * t553;
t515 = t546 * t561 + t560 * t591;
t523 = -t546 * t585 + t591 * t586;
t500 = t515 * t548 + t523 * t547;
t542 = sin(pkin(6));
t580 = t542 * t548;
t448 = Icges(4,5) * t500 + Icges(4,6) * t499 - Icges(4,3) * t580;
t501 = -t547 * t553 - t548 * t555;
t502 = -t515 * t547 + t523 * t548;
t581 = t542 * t547;
t449 = Icges(4,5) * t502 + Icges(4,6) * t501 + Icges(4,3) * t581;
t563 = t587 * t591;
t517 = -t547 * t546 + t548 * t563;
t570 = t546 * t587;
t518 = t547 * t591 + t548 * t570;
t486 = Icges(3,5) * t518 + Icges(3,6) * t517 - Icges(3,3) * t580;
t519 = -t548 * t546 - t547 * t563;
t520 = -t547 * t570 + t548 * t591;
t487 = Icges(3,5) * t520 + Icges(3,6) * t519 + Icges(3,3) * t581;
t602 = ((t448 + t486) * t548 + (-t449 - t487) * t547) * t542;
t545 = sin(qJ(4));
t590 = cos(qJ(4));
t474 = t500 * t590 - t545 * t580;
t541 = sin(pkin(12));
t543 = cos(pkin(12));
t442 = -t474 * t541 - t499 * t543;
t584 = t499 * t541;
t443 = t474 * t543 - t584;
t573 = t542 * t590;
t473 = t500 * t545 + t548 * t573;
t601 = -Icges(5,4) * t474 + Icges(6,5) * t443 + Icges(5,6) * t499 + Icges(6,6) * t442 + t603 * t473;
t476 = t502 * t590 + t545 * t581;
t444 = -t476 * t541 - t501 * t543;
t583 = t501 * t541;
t445 = t476 * t543 - t583;
t475 = t502 * t545 - t547 * t573;
t600 = -Icges(5,4) * t476 + Icges(6,5) * t445 + Icges(5,6) * t501 + Icges(6,6) * t444 + t603 * t475;
t514 = t555 * t542;
t505 = t514 * t590 + t545 * t587;
t513 = t523 * t542;
t469 = -t505 * t541 - t513 * t543;
t582 = t513 * t541;
t470 = t505 * t543 - t582;
t504 = t514 * t545 - t587 * t590;
t599 = -Icges(5,4) * t505 + Icges(6,5) * t470 + Icges(5,6) * t513 + Icges(6,6) * t469 + t603 * t504;
t598 = Icges(4,5) * t514 + Icges(4,6) * t513 + (Icges(3,5) * t546 + Icges(3,6) * t591) * t542 + (Icges(4,3) + Icges(3,3)) * t587;
t589 = pkin(2) * t591;
t588 = pkin(5) * t543;
t571 = pkin(2) * t570 - qJ(3) * t542;
t498 = -t547 * t571 + t548 * t589;
t521 = qJD(1) * (pkin(1) * t548 + pkin(8) * t581);
t532 = qJD(2) * t587 + qJD(1);
t578 = t532 * t498 + t521;
t576 = qJD(2) * t542;
t531 = t547 * t576;
t477 = -qJD(4) * t501 + t531;
t577 = qJD(1) * (pkin(1) * t547 - pkin(8) * t580);
t575 = qJD(3) * t548;
t497 = t547 * t589 + t548 * t571;
t572 = t548 * t576;
t574 = qJD(3) * t587 + t497 * t531 + t498 * t572;
t506 = -qJD(4) * t513 + t532;
t524 = t542 * t546 * pkin(2) + qJ(3) * t587;
t568 = qJD(2) * (-t514 * rSges(4,1) - t513 * rSges(4,2) - rSges(4,3) * t587 - t524);
t567 = qJD(2) * (-pkin(3) * t514 + pkin(9) * t513 - t524);
t566 = qJD(3) * t581 - t577;
t478 = -qJD(4) * t499 - t572;
t463 = pkin(3) * t500 - pkin(9) * t499;
t464 = pkin(3) * t502 - pkin(9) * t501;
t559 = t463 * t531 + t464 * t572 + t574;
t433 = pkin(4) * t474 + qJ(5) * t473;
t556 = qJD(5) * t504 + t477 * t433 + t559;
t554 = t532 * t464 + (t547 * t567 - t575) * t542 + t578;
t552 = (-t463 - t497) * t532 + t567 * t580 + t566;
t434 = pkin(4) * t476 + qJ(5) * t475;
t551 = qJD(5) * t473 + t506 * t434 + t554;
t465 = pkin(4) * t505 + qJ(5) * t504;
t550 = qJD(5) * t475 + t478 * t465 + t552;
t540 = pkin(12) + qJ(6);
t539 = cos(t540);
t538 = sin(t540);
t527 = rSges(2,1) * t548 - rSges(2,2) * t547;
t526 = rSges(2,1) * t547 + rSges(2,2) * t548;
t512 = t587 * rSges(3,3) + (rSges(3,1) * t546 + rSges(3,2) * t591) * t542;
t511 = Icges(3,5) * t587 + (Icges(3,1) * t546 + Icges(3,4) * t591) * t542;
t510 = Icges(3,6) * t587 + (Icges(3,4) * t546 + Icges(3,2) * t591) * t542;
t493 = rSges(3,1) * t520 + rSges(3,2) * t519 + rSges(3,3) * t581;
t492 = rSges(3,1) * t518 + rSges(3,2) * t517 - rSges(3,3) * t580;
t491 = Icges(3,1) * t520 + Icges(3,4) * t519 + Icges(3,5) * t581;
t490 = Icges(3,1) * t518 + Icges(3,4) * t517 - Icges(3,5) * t580;
t489 = Icges(3,4) * t520 + Icges(3,2) * t519 + Icges(3,6) * t581;
t488 = Icges(3,4) * t518 + Icges(3,2) * t517 - Icges(3,6) * t580;
t484 = Icges(4,1) * t514 + Icges(4,4) * t513 + Icges(4,5) * t587;
t483 = Icges(4,4) * t514 + Icges(4,2) * t513 + Icges(4,6) * t587;
t468 = qJD(6) * t504 + t506;
t467 = t505 * t539 - t513 * t538;
t466 = -t505 * t538 - t513 * t539;
t462 = rSges(5,1) * t505 - rSges(5,2) * t504 - rSges(5,3) * t513;
t461 = Icges(5,1) * t505 - Icges(5,4) * t504 - Icges(5,5) * t513;
t459 = Icges(5,5) * t505 - Icges(5,6) * t504 - Icges(5,3) * t513;
t456 = rSges(4,1) * t502 + rSges(4,2) * t501 + rSges(4,3) * t581;
t455 = rSges(4,1) * t500 + rSges(4,2) * t499 - rSges(4,3) * t580;
t453 = Icges(4,1) * t502 + Icges(4,4) * t501 + Icges(4,5) * t581;
t452 = Icges(4,1) * t500 + Icges(4,4) * t499 - Icges(4,5) * t580;
t451 = Icges(4,4) * t502 + Icges(4,2) * t501 + Icges(4,6) * t581;
t450 = Icges(4,4) * t500 + Icges(4,2) * t499 - Icges(4,6) * t580;
t447 = t493 * t532 - t512 * t531 + t521;
t446 = -t492 * t532 - t512 * t572 - t577;
t441 = (t492 * t547 + t493 * t548) * t576;
t440 = t476 * t539 - t501 * t538;
t439 = -t476 * t538 - t501 * t539;
t438 = t474 * t539 - t499 * t538;
t437 = -t474 * t538 - t499 * t539;
t436 = qJD(6) * t473 + t478;
t435 = qJD(6) * t475 + t477;
t430 = rSges(6,1) * t470 + rSges(6,2) * t469 + rSges(6,3) * t504;
t429 = Icges(6,1) * t470 + Icges(6,4) * t469 + Icges(6,5) * t504;
t428 = Icges(6,4) * t470 + Icges(6,2) * t469 + Icges(6,6) * t504;
t426 = rSges(5,1) * t476 - rSges(5,2) * t475 - rSges(5,3) * t501;
t425 = rSges(5,1) * t474 - rSges(5,2) * t473 - rSges(5,3) * t499;
t424 = Icges(5,1) * t476 - Icges(5,4) * t475 - Icges(5,5) * t501;
t423 = Icges(5,1) * t474 - Icges(5,4) * t473 - Icges(5,5) * t499;
t420 = Icges(5,5) * t476 - Icges(5,6) * t475 - Icges(5,3) * t501;
t419 = Icges(5,5) * t474 - Icges(5,6) * t473 - Icges(5,3) * t499;
t418 = -pkin(5) * t582 + pkin(10) * t504 + t505 * t588;
t417 = rSges(7,1) * t467 + rSges(7,2) * t466 + rSges(7,3) * t504;
t416 = Icges(7,1) * t467 + Icges(7,4) * t466 + Icges(7,5) * t504;
t415 = Icges(7,4) * t467 + Icges(7,2) * t466 + Icges(7,6) * t504;
t414 = Icges(7,5) * t467 + Icges(7,6) * t466 + Icges(7,3) * t504;
t412 = t456 * t532 + (t547 * t568 - t575) * t542 + t578;
t411 = (-t455 - t497) * t532 + t568 * t580 + t566;
t410 = rSges(6,1) * t445 + rSges(6,2) * t444 + rSges(6,3) * t475;
t409 = rSges(6,1) * t443 + rSges(6,2) * t442 + rSges(6,3) * t473;
t408 = Icges(6,1) * t445 + Icges(6,4) * t444 + Icges(6,5) * t475;
t407 = Icges(6,1) * t443 + Icges(6,4) * t442 + Icges(6,5) * t473;
t406 = Icges(6,4) * t445 + Icges(6,2) * t444 + Icges(6,6) * t475;
t405 = Icges(6,4) * t443 + Icges(6,2) * t442 + Icges(6,6) * t473;
t402 = rSges(7,1) * t440 + rSges(7,2) * t439 + rSges(7,3) * t475;
t401 = rSges(7,1) * t438 + rSges(7,2) * t437 + rSges(7,3) * t473;
t400 = Icges(7,1) * t440 + Icges(7,4) * t439 + Icges(7,5) * t475;
t399 = Icges(7,1) * t438 + Icges(7,4) * t437 + Icges(7,5) * t473;
t398 = Icges(7,4) * t440 + Icges(7,2) * t439 + Icges(7,6) * t475;
t397 = Icges(7,4) * t438 + Icges(7,2) * t437 + Icges(7,6) * t473;
t396 = Icges(7,5) * t440 + Icges(7,6) * t439 + Icges(7,3) * t475;
t395 = Icges(7,5) * t438 + Icges(7,6) * t437 + Icges(7,3) * t473;
t394 = -pkin(5) * t583 + pkin(10) * t475 + t476 * t588;
t393 = -pkin(5) * t584 + pkin(10) * t473 + t474 * t588;
t392 = (t455 * t547 + t456 * t548) * t576 + t574;
t391 = t426 * t506 - t462 * t477 + t554;
t390 = -t425 * t506 + t462 * t478 + t552;
t389 = t425 * t477 - t426 * t478 + t559;
t388 = t410 * t506 + (-t430 - t465) * t477 + t551;
t387 = t430 * t478 + (-t409 - t433) * t506 + t550;
t386 = t409 * t477 + (-t410 - t434) * t478 + t556;
t385 = t394 * t506 + t402 * t468 - t417 * t435 + (-t418 - t465) * t477 + t551;
t384 = -t401 * t468 + t417 * t436 + t418 * t478 + (-t393 - t433) * t506 + t550;
t383 = t393 * t477 + t401 * t435 - t402 * t436 + (-t394 - t434) * t478 + t556;
t1 = m(5) * (t389 ^ 2 + t390 ^ 2 + t391 ^ 2) / 0.2e1 + m(4) * (t392 ^ 2 + t411 ^ 2 + t412 ^ 2) / 0.2e1 + m(3) * (t441 ^ 2 + t446 ^ 2 + t447 ^ 2) / 0.2e1 + m(7) * (t383 ^ 2 + t384 ^ 2 + t385 ^ 2) / 0.2e1 + m(6) * (t386 ^ 2 + t387 ^ 2 + t388 ^ 2) / 0.2e1 + t468 * ((t396 * t504 + t398 * t466 + t400 * t467) * t435 + (t395 * t504 + t397 * t466 + t399 * t467) * t436 + (t504 * t414 + t466 * t415 + t467 * t416) * t468) / 0.2e1 + t435 * ((t475 * t396 + t439 * t398 + t440 * t400) * t435 + (t395 * t475 + t397 * t439 + t399 * t440) * t436 + (t414 * t475 + t415 * t439 + t416 * t440) * t468) / 0.2e1 + t436 * ((t396 * t473 + t398 * t437 + t400 * t438) * t435 + (t473 * t395 + t437 * t397 + t438 * t399) * t436 + (t414 * t473 + t415 * t437 + t416 * t438) * t468) / 0.2e1 + ((t428 * t444 + t429 * t445 - t459 * t501 + t461 * t476 + t599 * t475) * t506 + (t405 * t444 + t407 * t445 - t419 * t501 + t423 * t476 + t601 * t475) * t478 + (t444 * t406 + t445 * t408 - t501 * t420 + t476 * t424 + t600 * t475) * t477) * t477 / 0.2e1 + ((t428 * t442 + t429 * t443 - t459 * t499 + t461 * t474 + t599 * t473) * t506 + (t442 * t405 + t443 * t407 - t499 * t419 + t474 * t423 + t601 * t473) * t478 + (t406 * t442 + t408 * t443 - t420 * t499 + t424 * t474 + t600 * t473) * t477) * t478 / 0.2e1 + ((t469 * t428 + t470 * t429 - t513 * t459 + t505 * t461 + t599 * t504) * t506 + (t405 * t469 + t407 * t470 - t419 * t513 + t423 * t505 + t601 * t504) * t478 + (t406 * t469 + t408 * t470 - t420 * t513 + t424 * t505 + t600 * t504) * t477) * t506 / 0.2e1 + (((t449 * t587 + t513 * t451 + t514 * t453) * t547 - (t448 * t587 + t513 * t450 + t514 * t452) * t548) * t576 + (t587 * t487 + (t489 * t591 + t491 * t546) * t542) * t531 - (t587 * t486 + (t488 * t591 + t490 * t546) * t542) * t572 + (t513 * t483 + t514 * t484 + (t510 * t591 + t511 * t546) * t542 + t598 * t587) * t532) * t532 / 0.2e1 + (m(2) * (t526 ^ 2 + t527 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t450 * t501 - t452 * t502 - t488 * t519 - t490 * t520) * t548 + (t501 * t451 + t502 * t453 + t519 * t489 + t520 * t491 - t602) * t547) * t576 + (t483 * t501 + t484 * t502 + t510 * t519 + t511 * t520 + t598 * t581) * t532) * t531 / 0.2e1 - (((-t499 * t450 - t500 * t452 - t517 * t488 - t518 * t490 + t602) * t548 + (t499 * t451 + t500 * t453 + t489 * t517 + t491 * t518) * t547) * t576 + (t483 * t499 + t484 * t500 + t510 * t517 + t511 * t518 - t598 * t580) * t532) * t572 / 0.2e1;
T  = t1;
