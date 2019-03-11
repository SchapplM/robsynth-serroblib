% Calculate kinetic energy for
% S6PRPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:22:21
% EndTime: 2019-03-08 20:22:24
% DurationCPUTime: 3.29s
% Computational Cost: add. (3654->374), mult. (7993->586), div. (0->0), fcn. (10070->14), ass. (0->173)
t552 = sin(qJ(2));
t587 = sin(pkin(12));
t588 = cos(pkin(12));
t592 = cos(qJ(2));
t533 = -t552 * t588 - t587 * t592;
t546 = sin(pkin(11));
t548 = cos(pkin(11));
t534 = -t552 * t587 + t592 * t588;
t549 = cos(pkin(6));
t560 = t549 * t534;
t509 = t546 * t533 + t548 * t560;
t527 = t533 * t549;
t510 = -t527 * t548 + t534 * t546;
t547 = sin(pkin(6));
t585 = t547 * t548;
t449 = Icges(4,5) * t510 + Icges(4,6) * t509 - Icges(4,3) * t585;
t511 = t533 * t548 - t546 * t560;
t512 = t527 * t546 + t534 * t548;
t586 = t546 * t547;
t450 = Icges(4,5) * t512 + Icges(4,6) * t511 + Icges(4,3) * t586;
t525 = t534 * t547;
t526 = t533 * t547;
t490 = -Icges(4,5) * t526 + Icges(4,6) * t525 + Icges(4,3) * t549;
t573 = t549 * t592;
t529 = -t546 * t552 + t548 * t573;
t581 = t549 * t552;
t530 = t546 * t592 + t548 * t581;
t494 = Icges(3,5) * t530 + Icges(3,6) * t529 - Icges(3,3) * t585;
t531 = -t546 * t573 - t548 * t552;
t532 = -t546 * t581 + t548 * t592;
t495 = Icges(3,5) * t532 + Icges(3,6) * t531 + Icges(3,3) * t586;
t521 = Icges(3,3) * t549 + (Icges(3,5) * t552 + Icges(3,6) * t592) * t547;
t597 = -(t490 + t521) * t549 + (t494 + t449) * t585 - (t495 + t450) * t586;
t593 = qJD(2) ^ 2;
t591 = pkin(2) * t592;
t554 = cos(qJ(4));
t590 = pkin(4) * t554;
t551 = sin(qJ(4));
t584 = t547 * t551;
t583 = t547 * t554;
t582 = t549 * t551;
t580 = qJ(4) + qJ(5);
t535 = pkin(2) * t547 * t552 + qJ(3) * t549;
t579 = t526 * pkin(3) + t525 * pkin(8) - t535;
t578 = qJD(2) * t547;
t540 = t546 * t578;
t484 = -qJD(4) * t511 + t540;
t544 = qJD(2) * t549;
t517 = -qJD(4) * t525 + t544;
t577 = qJD(3) * t548;
t575 = t546 * t584;
t574 = t548 * t584;
t447 = -qJD(5) * t511 + t484;
t488 = -qJD(5) * t525 + t517;
t572 = t548 * t578;
t571 = cos(t580);
t570 = pkin(2) * t581 - qJ(3) * t547;
t568 = (rSges(4,1) * t526 - rSges(4,2) * t525 - rSges(4,3) * t549 - t535) * t547;
t502 = t546 * t591 + t548 * t570;
t503 = -t546 * t570 + t548 * t591;
t567 = qJD(3) * t549 + t502 * t540 + t503 * t572 + qJD(1);
t563 = t547 * t571;
t485 = -qJD(4) * t509 - t572;
t448 = -qJD(5) * t509 + t485;
t464 = t510 * pkin(3) - t509 * pkin(8);
t465 = t512 * pkin(3) - t511 * pkin(8);
t562 = t464 * t540 + t465 * t572 + t567;
t414 = -pkin(4) * t574 - pkin(9) * t509 + t510 * t590;
t415 = pkin(4) * t575 - pkin(9) * t511 + t512 * t590;
t561 = t484 * t414 - t415 * t485 + t562;
t489 = t503 * t544;
t559 = t465 * t544 + t489 + (qJD(2) * t546 * t579 - t577) * t547;
t539 = qJD(3) * t586;
t558 = t539 + ((-t464 - t502) * t549 + t579 * t585) * qJD(2);
t446 = pkin(4) * t582 - pkin(9) * t525 - t526 * t590;
t557 = t517 * t415 - t446 * t484 + t559;
t556 = -t414 * t517 + t485 * t446 + t558;
t553 = cos(qJ(6));
t550 = sin(qJ(6));
t545 = sin(t580);
t524 = t549 * rSges(3,3) + (rSges(3,1) * t552 + rSges(3,2) * t592) * t547;
t523 = Icges(3,5) * t549 + (Icges(3,1) * t552 + Icges(3,4) * t592) * t547;
t522 = Icges(3,6) * t549 + (Icges(3,4) * t552 + Icges(3,2) * t592) * t547;
t516 = -t526 * t554 + t582;
t515 = t526 * t551 + t549 * t554;
t514 = -t526 * t571 + t549 * t545;
t513 = -t526 * t545 - t549 * t571;
t501 = rSges(3,1) * t532 + rSges(3,2) * t531 + rSges(3,3) * t586;
t500 = rSges(3,1) * t530 + rSges(3,2) * t529 - rSges(3,3) * t585;
t499 = Icges(3,1) * t532 + Icges(3,4) * t531 + Icges(3,5) * t586;
t498 = Icges(3,1) * t530 + Icges(3,4) * t529 - Icges(3,5) * t585;
t497 = Icges(3,4) * t532 + Icges(3,2) * t531 + Icges(3,6) * t586;
t496 = Icges(3,4) * t530 + Icges(3,2) * t529 - Icges(3,6) * t585;
t492 = -Icges(4,1) * t526 + Icges(4,4) * t525 + Icges(4,5) * t549;
t491 = -Icges(4,4) * t526 + Icges(4,2) * t525 + Icges(4,6) * t549;
t483 = t512 * t554 + t575;
t482 = -t512 * t551 + t546 * t583;
t481 = t510 * t554 - t574;
t480 = -t510 * t551 - t548 * t583;
t479 = t512 * t571 + t545 * t586;
t478 = t512 * t545 - t546 * t563;
t477 = t510 * t571 - t545 * t585;
t476 = t510 * t545 + t548 * t563;
t475 = t514 * t553 - t525 * t550;
t474 = -t514 * t550 - t525 * t553;
t473 = pkin(5) * t514 + pkin(10) * t513;
t472 = (-t500 * t549 - t524 * t585) * qJD(2);
t471 = (t501 * t549 - t524 * t586) * qJD(2);
t470 = qJD(6) * t513 + t488;
t469 = rSges(5,1) * t516 + rSges(5,2) * t515 - rSges(5,3) * t525;
t468 = Icges(5,1) * t516 + Icges(5,4) * t515 - Icges(5,5) * t525;
t467 = Icges(5,4) * t516 + Icges(5,2) * t515 - Icges(5,6) * t525;
t466 = Icges(5,5) * t516 + Icges(5,6) * t515 - Icges(5,3) * t525;
t462 = rSges(6,1) * t514 - rSges(6,2) * t513 - rSges(6,3) * t525;
t461 = Icges(6,1) * t514 - Icges(6,4) * t513 - Icges(6,5) * t525;
t460 = Icges(6,4) * t514 - Icges(6,2) * t513 - Icges(6,6) * t525;
t459 = Icges(6,5) * t514 - Icges(6,6) * t513 - Icges(6,3) * t525;
t456 = rSges(4,1) * t512 + rSges(4,2) * t511 + rSges(4,3) * t586;
t455 = rSges(4,1) * t510 + rSges(4,2) * t509 - rSges(4,3) * t585;
t454 = Icges(4,1) * t512 + Icges(4,4) * t511 + Icges(4,5) * t586;
t453 = Icges(4,1) * t510 + Icges(4,4) * t509 - Icges(4,5) * t585;
t452 = Icges(4,4) * t512 + Icges(4,2) * t511 + Icges(4,6) * t586;
t451 = Icges(4,4) * t510 + Icges(4,2) * t509 - Icges(4,6) * t585;
t445 = t479 * t553 - t511 * t550;
t444 = -t479 * t550 - t511 * t553;
t443 = t477 * t553 - t509 * t550;
t442 = -t477 * t550 - t509 * t553;
t441 = qJD(1) + (t500 * t546 + t501 * t548) * t578;
t440 = pkin(5) * t479 + pkin(10) * t478;
t439 = pkin(5) * t477 + pkin(10) * t476;
t437 = qJD(6) * t476 + t448;
t436 = qJD(6) * t478 + t447;
t435 = rSges(5,1) * t483 + rSges(5,2) * t482 - rSges(5,3) * t511;
t434 = rSges(5,1) * t481 + rSges(5,2) * t480 - rSges(5,3) * t509;
t433 = Icges(5,1) * t483 + Icges(5,4) * t482 - Icges(5,5) * t511;
t432 = Icges(5,1) * t481 + Icges(5,4) * t480 - Icges(5,5) * t509;
t431 = Icges(5,4) * t483 + Icges(5,2) * t482 - Icges(5,6) * t511;
t430 = Icges(5,4) * t481 + Icges(5,2) * t480 - Icges(5,6) * t509;
t429 = Icges(5,5) * t483 + Icges(5,6) * t482 - Icges(5,3) * t511;
t428 = Icges(5,5) * t481 + Icges(5,6) * t480 - Icges(5,3) * t509;
t427 = rSges(7,1) * t475 + rSges(7,2) * t474 + rSges(7,3) * t513;
t426 = Icges(7,1) * t475 + Icges(7,4) * t474 + Icges(7,5) * t513;
t425 = Icges(7,4) * t475 + Icges(7,2) * t474 + Icges(7,6) * t513;
t424 = Icges(7,5) * t475 + Icges(7,6) * t474 + Icges(7,3) * t513;
t423 = rSges(6,1) * t479 - rSges(6,2) * t478 - rSges(6,3) * t511;
t422 = rSges(6,1) * t477 - rSges(6,2) * t476 - rSges(6,3) * t509;
t421 = Icges(6,1) * t479 - Icges(6,4) * t478 - Icges(6,5) * t511;
t420 = Icges(6,1) * t477 - Icges(6,4) * t476 - Icges(6,5) * t509;
t419 = Icges(6,4) * t479 - Icges(6,2) * t478 - Icges(6,6) * t511;
t418 = Icges(6,4) * t477 - Icges(6,2) * t476 - Icges(6,6) * t509;
t417 = Icges(6,5) * t479 - Icges(6,6) * t478 - Icges(6,3) * t511;
t416 = Icges(6,5) * t477 - Icges(6,6) * t476 - Icges(6,3) * t509;
t412 = t539 + ((-t455 - t502) * t549 + t548 * t568) * qJD(2);
t411 = -t547 * t577 + t489 + (t456 * t549 + t546 * t568) * qJD(2);
t409 = rSges(7,1) * t445 + rSges(7,2) * t444 + rSges(7,3) * t478;
t408 = rSges(7,1) * t443 + rSges(7,2) * t442 + rSges(7,3) * t476;
t407 = Icges(7,1) * t445 + Icges(7,4) * t444 + Icges(7,5) * t478;
t406 = Icges(7,1) * t443 + Icges(7,4) * t442 + Icges(7,5) * t476;
t405 = Icges(7,4) * t445 + Icges(7,2) * t444 + Icges(7,6) * t478;
t404 = Icges(7,4) * t443 + Icges(7,2) * t442 + Icges(7,6) * t476;
t403 = Icges(7,5) * t445 + Icges(7,6) * t444 + Icges(7,3) * t478;
t402 = Icges(7,5) * t443 + Icges(7,6) * t442 + Icges(7,3) * t476;
t401 = (t455 * t546 + t456 * t548) * t578 + t567;
t400 = -t434 * t517 + t469 * t485 + t558;
t399 = t435 * t517 - t469 * t484 + t559;
t398 = t434 * t484 - t435 * t485 + t562;
t397 = -t422 * t488 + t448 * t462 + t556;
t396 = t423 * t488 - t447 * t462 + t557;
t395 = t422 * t447 - t423 * t448 + t561;
t394 = -t408 * t470 + t427 * t437 - t439 * t488 + t448 * t473 + t556;
t393 = t409 * t470 - t427 * t436 + t440 * t488 - t447 * t473 + t557;
t392 = t408 * t436 - t409 * t437 + t439 * t447 - t440 * t448 + t561;
t1 = m(7) * (t392 ^ 2 + t393 ^ 2 + t394 ^ 2) / 0.2e1 + t436 * ((t478 * t403 + t444 * t405 + t445 * t407) * t436 + (t402 * t478 + t404 * t444 + t406 * t445) * t437 + (t424 * t478 + t425 * t444 + t426 * t445) * t470) / 0.2e1 + t437 * ((t403 * t476 + t405 * t442 + t407 * t443) * t436 + (t476 * t402 + t442 * t404 + t443 * t406) * t437 + (t424 * t476 + t425 * t442 + t426 * t443) * t470) / 0.2e1 + t470 * ((t403 * t513 + t405 * t474 + t407 * t475) * t436 + (t402 * t513 + t404 * t474 + t406 * t475) * t437 + (t513 * t424 + t474 * t425 + t475 * t426) * t470) / 0.2e1 + ((t549 * t495 + (t497 * t592 + t499 * t552) * t547) * t540 - (t549 * t494 + (t496 * t592 + t498 * t552) * t547) * t572 + (t549 * t521 + (t522 * t592 + t523 * t552) * t547) * t544) * t544 / 0.2e1 + t484 * ((-t511 * t429 + t482 * t431 + t483 * t433) * t484 + (-t428 * t511 + t430 * t482 + t432 * t483) * t485 + (-t466 * t511 + t467 * t482 + t468 * t483) * t517) / 0.2e1 + t485 * ((-t429 * t509 + t431 * t480 + t433 * t481) * t484 + (-t509 * t428 + t480 * t430 + t481 * t432) * t485 + (-t466 * t509 + t467 * t480 + t468 * t481) * t517) / 0.2e1 + t517 * ((-t429 * t525 + t431 * t515 + t433 * t516) * t484 + (-t428 * t525 + t430 * t515 + t432 * t516) * t485 + (-t466 * t525 + t467 * t515 + t468 * t516) * t517) / 0.2e1 + t447 * ((-t511 * t417 - t478 * t419 + t479 * t421) * t447 + (-t416 * t511 - t418 * t478 + t420 * t479) * t448 + (-t459 * t511 - t460 * t478 + t461 * t479) * t488) / 0.2e1 + t448 * ((-t417 * t509 - t419 * t476 + t421 * t477) * t447 + (-t509 * t416 - t476 * t418 + t477 * t420) * t448 + (-t459 * t509 - t460 * t476 + t461 * t477) * t488) / 0.2e1 + t488 * ((-t417 * t525 - t419 * t513 + t421 * t514) * t447 + (-t416 * t525 - t418 * t513 + t420 * t514) * t448 + (-t525 * t459 - t513 * t460 + t514 * t461) * t488) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t441 ^ 2 + t471 ^ 2 + t472 ^ 2) / 0.2e1 + m(4) * (t401 ^ 2 + t411 ^ 2 + t412 ^ 2) / 0.2e1 + m(5) * (t398 ^ 2 + t399 ^ 2 + t400 ^ 2) / 0.2e1 + m(6) * (t395 ^ 2 + t396 ^ 2 + t397 ^ 2) / 0.2e1 - ((t452 * t509 + t454 * t510 + t497 * t529 + t499 * t530) * t586 + (t491 * t509 + t492 * t510 + t522 * t529 + t523 * t530) * t549 + (-t451 * t509 - t453 * t510 - t496 * t529 - t498 * t530 + t597) * t585) * t593 * t585 / 0.2e1 + (t549 * ((t490 * t549 + t491 * t525 - t492 * t526) * t549 + ((t450 * t549 + t452 * t525 - t454 * t526) * t546 - (t449 * t549 + t451 * t525 - t453 * t526) * t548) * t547) + ((-t451 * t511 - t453 * t512 - t496 * t531 - t498 * t532) * t585 + (t491 * t511 + t492 * t512 + t522 * t531 + t523 * t532) * t549 + (t452 * t511 + t454 * t512 + t497 * t531 + t499 * t532 - t597) * t586) * t586) * t593 / 0.2e1;
T  = t1;
