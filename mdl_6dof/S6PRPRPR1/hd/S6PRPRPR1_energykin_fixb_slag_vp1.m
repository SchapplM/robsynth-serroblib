% Calculate kinetic energy for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:24:40
% EndTime: 2019-03-08 19:24:43
% DurationCPUTime: 3.41s
% Computational Cost: add. (3534->369), mult. (7717->554), div. (0->0), fcn. (9710->14), ass. (0->166)
t601 = -Icges(5,3) - Icges(6,3);
t552 = sin(qJ(2));
t584 = sin(pkin(11));
t585 = cos(pkin(11));
t589 = cos(qJ(2));
t532 = -t552 * t585 - t584 * t589;
t545 = sin(pkin(10));
t547 = cos(pkin(10));
t533 = -t552 * t584 + t589 * t585;
t548 = cos(pkin(6));
t559 = t548 * t533;
t508 = t545 * t532 + t547 * t559;
t526 = t532 * t548;
t509 = -t526 * t547 + t533 * t545;
t546 = sin(pkin(6));
t582 = t547 * t546;
t452 = Icges(4,5) * t509 + Icges(4,6) * t508 - Icges(4,3) * t582;
t510 = t532 * t547 - t545 * t559;
t511 = t526 * t545 + t533 * t547;
t583 = t545 * t546;
t453 = Icges(4,5) * t511 + Icges(4,6) * t510 + Icges(4,3) * t583;
t524 = t533 * t546;
t525 = t532 * t546;
t489 = -Icges(4,5) * t525 + Icges(4,6) * t524 + Icges(4,3) * t548;
t572 = t548 * t589;
t528 = -t545 * t552 + t547 * t572;
t580 = t548 * t552;
t529 = t545 * t589 + t547 * t580;
t493 = Icges(3,5) * t529 + Icges(3,6) * t528 - Icges(3,3) * t582;
t530 = -t545 * t572 - t547 * t552;
t531 = -t545 * t580 + t547 * t589;
t494 = Icges(3,5) * t531 + Icges(3,6) * t530 + Icges(3,3) * t583;
t520 = Icges(3,3) * t548 + (Icges(3,5) * t552 + Icges(3,6) * t589) * t546;
t600 = -(t489 + t520) * t548 + (t452 + t493) * t582 - (t494 + t453) * t583;
t575 = qJ(4) + pkin(12);
t544 = sin(t575);
t568 = cos(t575);
t562 = t546 * t568;
t476 = t509 * t544 + t547 * t562;
t477 = t509 * t568 - t544 * t582;
t551 = sin(qJ(4));
t554 = cos(qJ(4));
t480 = -t509 * t551 - t554 * t582;
t573 = t551 * t582;
t481 = t509 * t554 - t573;
t599 = Icges(5,5) * t481 + Icges(6,5) * t477 + Icges(5,6) * t480 - Icges(6,6) * t476 + t601 * t508;
t478 = t511 * t544 - t545 * t562;
t479 = t511 * t568 + t544 * t583;
t482 = -t511 * t551 + t554 * t583;
t574 = t551 * t583;
t483 = t511 * t554 + t574;
t598 = Icges(5,5) * t483 + Icges(6,5) * t479 + Icges(5,6) * t482 - Icges(6,6) * t478 + t601 * t510;
t512 = -t525 * t544 - t548 * t568;
t513 = -t525 * t568 + t548 * t544;
t514 = t525 * t551 + t548 * t554;
t581 = t548 * t551;
t515 = -t525 * t554 + t581;
t597 = Icges(5,5) * t515 + Icges(6,5) * t513 + Icges(5,6) * t514 - Icges(6,6) * t512 + t601 * t524;
t593 = qJD(2) ^ 2;
t588 = pkin(2) * t589;
t587 = pkin(4) * t554;
t534 = pkin(2) * t546 * t552 + qJ(3) * t548;
t579 = pkin(3) * t525 + pkin(8) * t524 - t534;
t578 = qJD(2) * t546;
t539 = t545 * t578;
t484 = -qJD(4) * t510 + t539;
t543 = qJD(2) * t548;
t516 = -qJD(4) * t524 + t543;
t577 = qJD(3) * t547;
t571 = t547 * t578;
t570 = pkin(2) * t580 - qJ(3) * t546;
t567 = (rSges(4,1) * t525 - rSges(4,2) * t524 - rSges(4,3) * t548 - t534) * t546;
t501 = t545 * t588 + t547 * t570;
t502 = -t545 * t570 + t547 * t588;
t566 = qJD(3) * t548 + t501 * t539 + t502 * t571 + qJD(1);
t485 = -qJD(4) * t508 - t571;
t464 = pkin(3) * t509 - pkin(8) * t508;
t465 = pkin(3) * t511 - pkin(8) * t510;
t561 = t464 * t539 + t465 * t571 + t566;
t416 = -pkin(4) * t573 - qJ(5) * t508 + t509 * t587;
t560 = -qJD(5) * t524 + t484 * t416 + t561;
t488 = t502 * t543;
t558 = t465 * t543 + t488 + (qJD(2) * t545 * t579 - t577) * t546;
t538 = qJD(3) * t583;
t557 = t538 + ((-t464 - t501) * t548 + t579 * t582) * qJD(2);
t417 = pkin(4) * t574 - qJ(5) * t510 + t511 * t587;
t556 = -qJD(5) * t508 + t516 * t417 + t558;
t448 = pkin(4) * t581 - qJ(5) * t524 - t525 * t587;
t555 = -qJD(5) * t510 + t485 * t448 + t557;
t553 = cos(qJ(6));
t550 = sin(qJ(6));
t523 = t548 * rSges(3,3) + (rSges(3,1) * t552 + rSges(3,2) * t589) * t546;
t522 = Icges(3,5) * t548 + (Icges(3,1) * t552 + Icges(3,4) * t589) * t546;
t521 = Icges(3,6) * t548 + (Icges(3,4) * t552 + Icges(3,2) * t589) * t546;
t500 = rSges(3,1) * t531 + rSges(3,2) * t530 + rSges(3,3) * t583;
t499 = rSges(3,1) * t529 + rSges(3,2) * t528 - rSges(3,3) * t582;
t498 = Icges(3,1) * t531 + Icges(3,4) * t530 + Icges(3,5) * t583;
t497 = Icges(3,1) * t529 + Icges(3,4) * t528 - Icges(3,5) * t582;
t496 = Icges(3,4) * t531 + Icges(3,2) * t530 + Icges(3,6) * t583;
t495 = Icges(3,4) * t529 + Icges(3,2) * t528 - Icges(3,6) * t582;
t491 = -Icges(4,1) * t525 + Icges(4,4) * t524 + Icges(4,5) * t548;
t490 = -Icges(4,4) * t525 + Icges(4,2) * t524 + Icges(4,6) * t548;
t475 = t513 * t553 - t524 * t550;
t474 = -t513 * t550 - t524 * t553;
t473 = qJD(6) * t512 + t516;
t472 = pkin(5) * t513 + pkin(9) * t512;
t471 = (-t499 * t548 - t523 * t582) * qJD(2);
t470 = (t500 * t548 - t523 * t583) * qJD(2);
t469 = rSges(5,1) * t515 + rSges(5,2) * t514 - rSges(5,3) * t524;
t468 = Icges(5,1) * t515 + Icges(5,4) * t514 - Icges(5,5) * t524;
t467 = Icges(5,4) * t515 + Icges(5,2) * t514 - Icges(5,6) * t524;
t460 = rSges(4,1) * t511 + rSges(4,2) * t510 + rSges(4,3) * t583;
t459 = rSges(4,1) * t509 + rSges(4,2) * t508 - rSges(4,3) * t582;
t458 = rSges(6,1) * t513 - rSges(6,2) * t512 - rSges(6,3) * t524;
t457 = Icges(4,1) * t511 + Icges(4,4) * t510 + Icges(4,5) * t583;
t456 = Icges(4,1) * t509 + Icges(4,4) * t508 - Icges(4,5) * t582;
t455 = Icges(4,4) * t511 + Icges(4,2) * t510 + Icges(4,6) * t583;
t454 = Icges(4,4) * t509 + Icges(4,2) * t508 - Icges(4,6) * t582;
t451 = Icges(6,1) * t513 - Icges(6,4) * t512 - Icges(6,5) * t524;
t450 = Icges(6,4) * t513 - Icges(6,2) * t512 - Icges(6,6) * t524;
t447 = t479 * t553 - t510 * t550;
t446 = -t479 * t550 - t510 * t553;
t445 = t477 * t553 - t508 * t550;
t444 = -t477 * t550 - t508 * t553;
t443 = qJD(1) + (t499 * t545 + t500 * t547) * t578;
t442 = qJD(6) * t476 + t485;
t441 = qJD(6) * t478 + t484;
t440 = pkin(5) * t479 + pkin(9) * t478;
t439 = pkin(5) * t477 + pkin(9) * t476;
t437 = rSges(5,1) * t483 + rSges(5,2) * t482 - rSges(5,3) * t510;
t436 = rSges(5,1) * t481 + rSges(5,2) * t480 - rSges(5,3) * t508;
t435 = Icges(5,1) * t483 + Icges(5,4) * t482 - Icges(5,5) * t510;
t434 = Icges(5,1) * t481 + Icges(5,4) * t480 - Icges(5,5) * t508;
t433 = Icges(5,4) * t483 + Icges(5,2) * t482 - Icges(5,6) * t510;
t432 = Icges(5,4) * t481 + Icges(5,2) * t480 - Icges(5,6) * t508;
t429 = rSges(7,1) * t475 + rSges(7,2) * t474 + rSges(7,3) * t512;
t428 = Icges(7,1) * t475 + Icges(7,4) * t474 + Icges(7,5) * t512;
t427 = Icges(7,4) * t475 + Icges(7,2) * t474 + Icges(7,6) * t512;
t426 = Icges(7,5) * t475 + Icges(7,6) * t474 + Icges(7,3) * t512;
t425 = rSges(6,1) * t479 - rSges(6,2) * t478 - rSges(6,3) * t510;
t424 = rSges(6,1) * t477 - rSges(6,2) * t476 - rSges(6,3) * t508;
t423 = Icges(6,1) * t479 - Icges(6,4) * t478 - Icges(6,5) * t510;
t422 = Icges(6,1) * t477 - Icges(6,4) * t476 - Icges(6,5) * t508;
t421 = Icges(6,4) * t479 - Icges(6,2) * t478 - Icges(6,6) * t510;
t420 = Icges(6,4) * t477 - Icges(6,2) * t476 - Icges(6,6) * t508;
t414 = t538 + ((-t459 - t501) * t548 + t547 * t567) * qJD(2);
t413 = -t546 * t577 + t488 + (t460 * t548 + t545 * t567) * qJD(2);
t411 = rSges(7,1) * t447 + rSges(7,2) * t446 + rSges(7,3) * t478;
t410 = rSges(7,1) * t445 + rSges(7,2) * t444 + rSges(7,3) * t476;
t409 = Icges(7,1) * t447 + Icges(7,4) * t446 + Icges(7,5) * t478;
t408 = Icges(7,1) * t445 + Icges(7,4) * t444 + Icges(7,5) * t476;
t407 = Icges(7,4) * t447 + Icges(7,2) * t446 + Icges(7,6) * t478;
t406 = Icges(7,4) * t445 + Icges(7,2) * t444 + Icges(7,6) * t476;
t405 = Icges(7,5) * t447 + Icges(7,6) * t446 + Icges(7,3) * t478;
t404 = Icges(7,5) * t445 + Icges(7,6) * t444 + Icges(7,3) * t476;
t403 = (t459 * t545 + t460 * t547) * t578 + t566;
t402 = -t436 * t516 + t469 * t485 + t557;
t401 = t437 * t516 - t469 * t484 + t558;
t400 = t436 * t484 - t437 * t485 + t561;
t399 = t458 * t485 + (-t416 - t424) * t516 + t555;
t398 = t425 * t516 + (-t448 - t458) * t484 + t556;
t397 = t424 * t484 + (-t417 - t425) * t485 + t560;
t396 = -t410 * t473 + t429 * t442 + t472 * t485 + (-t416 - t439) * t516 + t555;
t395 = t411 * t473 - t429 * t441 + t440 * t516 + (-t448 - t472) * t484 + t556;
t394 = t410 * t441 - t411 * t442 + t439 * t484 + (-t417 - t440) * t485 + t560;
t1 = ((t548 * t494 + (t496 * t589 + t498 * t552) * t546) * t539 - (t548 * t493 + (t495 * t589 + t497 * t552) * t546) * t571 + (t548 * t520 + (t521 * t589 + t522 * t552) * t546) * t543) * t543 / 0.2e1 + t441 * ((t478 * t405 + t446 * t407 + t447 * t409) * t441 + (t404 * t478 + t406 * t446 + t408 * t447) * t442 + (t426 * t478 + t427 * t446 + t428 * t447) * t473) / 0.2e1 + t442 * ((t405 * t476 + t407 * t444 + t409 * t445) * t441 + (t476 * t404 + t444 * t406 + t445 * t408) * t442 + (t426 * t476 + t427 * t444 + t428 * t445) * t473) / 0.2e1 + t473 * ((t405 * t512 + t407 * t474 + t409 * t475) * t441 + (t404 * t512 + t406 * t474 + t408 * t475) * t442 + (t512 * t426 + t474 * t427 + t475 * t428) * t473) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t443 ^ 2 + t470 ^ 2 + t471 ^ 2) / 0.2e1 + m(4) * (t403 ^ 2 + t413 ^ 2 + t414 ^ 2) / 0.2e1 + m(5) * (t400 ^ 2 + t401 ^ 2 + t402 ^ 2) / 0.2e1 + m(6) * (t397 ^ 2 + t398 ^ 2 + t399 ^ 2) / 0.2e1 + m(7) * (t394 ^ 2 + t395 ^ 2 + t396 ^ 2) / 0.2e1 + ((-t450 * t478 + t451 * t479 + t467 * t482 + t468 * t483 - t510 * t597) * t516 + (-t420 * t478 + t422 * t479 + t432 * t482 + t434 * t483 - t510 * t599) * t485 + (-t478 * t421 + t479 * t423 + t482 * t433 + t483 * t435 - t598 * t510) * t484) * t484 / 0.2e1 + ((-t450 * t476 + t451 * t477 + t467 * t480 + t468 * t481 - t508 * t597) * t516 + (-t476 * t420 + t477 * t422 + t480 * t432 + t481 * t434 - t599 * t508) * t485 + (-t421 * t476 + t423 * t477 + t433 * t480 + t435 * t481 - t508 * t598) * t484) * t485 / 0.2e1 + ((-t512 * t450 + t513 * t451 + t514 * t467 + t515 * t468 - t597 * t524) * t516 + (-t420 * t512 + t422 * t513 + t432 * t514 + t434 * t515 - t524 * t599) * t485 + (-t421 * t512 + t423 * t513 + t433 * t514 + t435 * t515 - t524 * t598) * t484) * t516 / 0.2e1 - ((t455 * t508 + t457 * t509 + t496 * t528 + t498 * t529) * t583 + (t490 * t508 + t491 * t509 + t521 * t528 + t522 * t529) * t548 + (-t454 * t508 - t456 * t509 - t495 * t528 - t497 * t529 + t600) * t582) * t593 * t582 / 0.2e1 + (t548 * ((t489 * t548 + t490 * t524 - t491 * t525) * t548 + ((t453 * t548 + t455 * t524 - t457 * t525) * t545 - (t452 * t548 + t454 * t524 - t456 * t525) * t547) * t546) + ((-t454 * t510 - t456 * t511 - t495 * t530 - t497 * t531) * t582 + (t490 * t510 + t491 * t511 + t521 * t530 + t522 * t531) * t548 + (t455 * t510 + t457 * t511 + t496 * t530 + t498 * t531 - t600) * t583) * t583) * t593 / 0.2e1;
T  = t1;
