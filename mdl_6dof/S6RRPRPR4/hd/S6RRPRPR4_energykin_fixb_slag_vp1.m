% Calculate kinetic energy for
% S6RRPRPR4
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
% Datum: 2019-03-09 10:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR4_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRPR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:21:16
% EndTime: 2019-03-09 10:21:19
% DurationCPUTime: 3.68s
% Computational Cost: add. (3603->380), mult. (7769->563), div. (0->0), fcn. (9744->14), ass. (0->170)
t614 = -Icges(5,3) - Icges(6,3);
t555 = sin(qJ(2));
t596 = sin(pkin(11));
t597 = cos(pkin(11));
t602 = cos(qJ(2));
t534 = -t555 * t597 - t596 * t602;
t556 = sin(qJ(1));
t559 = cos(qJ(1));
t598 = cos(pkin(6));
t571 = t598 * t596;
t572 = t598 * t597;
t564 = -t555 * t571 + t572 * t602;
t509 = t556 * t534 + t559 * t564;
t527 = t555 * t572 + t571 * t602;
t535 = -t555 * t596 + t602 * t597;
t510 = t527 * t559 + t535 * t556;
t551 = sin(pkin(6));
t594 = t551 * t559;
t456 = Icges(4,5) * t510 + Icges(4,6) * t509 - Icges(4,3) * t594;
t511 = t534 * t559 - t556 * t564;
t512 = -t527 * t556 + t535 * t559;
t595 = t551 * t556;
t457 = Icges(4,5) * t512 + Icges(4,6) * t511 + Icges(4,3) * t595;
t574 = t598 * t602;
t529 = -t556 * t555 + t559 * t574;
t583 = t555 * t598;
t530 = t556 * t602 + t559 * t583;
t494 = Icges(3,5) * t530 + Icges(3,6) * t529 - Icges(3,3) * t594;
t531 = -t559 * t555 - t556 * t574;
t532 = -t556 * t583 + t559 * t602;
t495 = Icges(3,5) * t532 + Icges(3,6) * t531 + Icges(3,3) * t595;
t613 = ((t456 + t494) * t559 + (-t457 - t495) * t556) * t551;
t589 = qJ(4) + pkin(12);
t550 = sin(t589);
t580 = cos(t589);
t570 = t551 * t580;
t477 = t510 * t550 + t559 * t570;
t478 = t510 * t580 - t550 * t594;
t554 = sin(qJ(4));
t558 = cos(qJ(4));
t481 = -t510 * t554 - t558 * t594;
t587 = t554 * t594;
t482 = t510 * t558 - t587;
t612 = Icges(5,5) * t482 + Icges(6,5) * t478 + Icges(5,6) * t481 - Icges(6,6) * t477 + t614 * t509;
t479 = t512 * t550 - t556 * t570;
t480 = t512 * t580 + t550 * t595;
t483 = -t512 * t554 + t558 * t595;
t588 = t554 * t595;
t484 = t512 * t558 + t588;
t611 = Icges(5,5) * t484 + Icges(6,5) * t480 + Icges(5,6) * t483 - Icges(6,6) * t479 + t614 * t511;
t526 = t534 * t551;
t513 = -t526 * t550 - t580 * t598;
t514 = -t526 * t580 + t550 * t598;
t515 = t526 * t554 + t558 * t598;
t581 = t598 * t554;
t516 = -t526 * t558 + t581;
t525 = t535 * t551;
t610 = Icges(5,5) * t516 + Icges(6,5) * t514 + Icges(5,6) * t515 - Icges(6,6) * t513 + t614 * t525;
t609 = -Icges(4,5) * t526 + Icges(4,6) * t525 + (Icges(3,5) * t555 + Icges(3,6) * t602) * t551 + (Icges(4,3) + Icges(3,3)) * t598;
t601 = pkin(2) * t602;
t600 = pkin(4) * t558;
t584 = pkin(2) * t583 - qJ(3) * t551;
t508 = -t556 * t584 + t559 * t601;
t533 = qJD(1) * (pkin(1) * t559 + pkin(8) * t595);
t545 = qJD(2) * t598 + qJD(1);
t593 = t545 * t508 + t533;
t591 = qJD(2) * t551;
t544 = t556 * t591;
t485 = -qJD(4) * t511 + t544;
t592 = qJD(1) * (pkin(1) * t556 - pkin(8) * t594);
t590 = qJD(3) * t559;
t507 = t556 * t601 + t559 * t584;
t585 = t559 * t591;
t586 = qJD(3) * t598 + t507 * t544 + t508 * t585;
t517 = -qJD(4) * t525 + t545;
t536 = t551 * t555 * pkin(2) + qJ(3) * t598;
t579 = qJD(2) * (t526 * rSges(4,1) - t525 * rSges(4,2) - rSges(4,3) * t598 - t536);
t578 = qJD(2) * (t526 * pkin(3) + t525 * pkin(9) - t536);
t577 = qJD(3) * t595 - t592;
t486 = -qJD(4) * t509 - t585;
t471 = pkin(3) * t510 - pkin(9) * t509;
t472 = pkin(3) * t512 - pkin(9) * t511;
t569 = t471 * t544 + t472 * t585 + t586;
t417 = -pkin(4) * t587 - qJ(5) * t509 + t510 * t600;
t566 = -qJD(5) * t525 + t485 * t417 + t569;
t565 = t545 * t472 + (t556 * t578 - t590) * t551 + t593;
t563 = (-t471 - t507) * t545 + t578 * t594 + t577;
t418 = pkin(4) * t588 - qJ(5) * t511 + t512 * t600;
t562 = -qJD(5) * t509 + t517 * t418 + t565;
t449 = pkin(4) * t581 - qJ(5) * t525 - t526 * t600;
t561 = -qJD(5) * t511 + t486 * t449 + t563;
t557 = cos(qJ(6));
t553 = sin(qJ(6));
t540 = rSges(2,1) * t559 - rSges(2,2) * t556;
t539 = rSges(2,1) * t556 + rSges(2,2) * t559;
t524 = t598 * rSges(3,3) + (rSges(3,1) * t555 + rSges(3,2) * t602) * t551;
t523 = Icges(3,5) * t598 + (Icges(3,1) * t555 + Icges(3,4) * t602) * t551;
t522 = Icges(3,6) * t598 + (Icges(3,4) * t555 + Icges(3,2) * t602) * t551;
t501 = rSges(3,1) * t532 + rSges(3,2) * t531 + rSges(3,3) * t595;
t500 = rSges(3,1) * t530 + rSges(3,2) * t529 - rSges(3,3) * t594;
t499 = Icges(3,1) * t532 + Icges(3,4) * t531 + Icges(3,5) * t595;
t498 = Icges(3,1) * t530 + Icges(3,4) * t529 - Icges(3,5) * t594;
t497 = Icges(3,4) * t532 + Icges(3,2) * t531 + Icges(3,6) * t595;
t496 = Icges(3,4) * t530 + Icges(3,2) * t529 - Icges(3,6) * t594;
t492 = -Icges(4,1) * t526 + Icges(4,4) * t525 + Icges(4,5) * t598;
t491 = -Icges(4,4) * t526 + Icges(4,2) * t525 + Icges(4,6) * t598;
t476 = t514 * t557 - t525 * t553;
t475 = -t514 * t553 - t525 * t557;
t474 = qJD(6) * t513 + t517;
t473 = pkin(5) * t514 + pkin(10) * t513;
t470 = rSges(5,1) * t516 + rSges(5,2) * t515 - rSges(5,3) * t525;
t469 = Icges(5,1) * t516 + Icges(5,4) * t515 - Icges(5,5) * t525;
t468 = Icges(5,4) * t516 + Icges(5,2) * t515 - Icges(5,6) * t525;
t464 = rSges(4,1) * t512 + rSges(4,2) * t511 + rSges(4,3) * t595;
t463 = rSges(4,1) * t510 + rSges(4,2) * t509 - rSges(4,3) * t594;
t461 = Icges(4,1) * t512 + Icges(4,4) * t511 + Icges(4,5) * t595;
t460 = Icges(4,1) * t510 + Icges(4,4) * t509 - Icges(4,5) * t594;
t459 = Icges(4,4) * t512 + Icges(4,2) * t511 + Icges(4,6) * t595;
t458 = Icges(4,4) * t510 + Icges(4,2) * t509 - Icges(4,6) * t594;
t455 = t501 * t545 - t524 * t544 + t533;
t454 = -t500 * t545 - t524 * t585 - t592;
t453 = rSges(6,1) * t514 - rSges(6,2) * t513 - rSges(6,3) * t525;
t452 = Icges(6,1) * t514 - Icges(6,4) * t513 - Icges(6,5) * t525;
t451 = Icges(6,4) * t514 - Icges(6,2) * t513 - Icges(6,6) * t525;
t448 = (t500 * t556 + t501 * t559) * t591;
t447 = t480 * t557 - t511 * t553;
t446 = -t480 * t553 - t511 * t557;
t445 = t478 * t557 - t509 * t553;
t444 = -t478 * t553 - t509 * t557;
t443 = qJD(6) * t477 + t486;
t442 = qJD(6) * t479 + t485;
t441 = pkin(5) * t480 + pkin(10) * t479;
t440 = pkin(5) * t478 + pkin(10) * t477;
t438 = rSges(5,1) * t484 + rSges(5,2) * t483 - rSges(5,3) * t511;
t437 = rSges(5,1) * t482 + rSges(5,2) * t481 - rSges(5,3) * t509;
t436 = Icges(5,1) * t484 + Icges(5,4) * t483 - Icges(5,5) * t511;
t435 = Icges(5,1) * t482 + Icges(5,4) * t481 - Icges(5,5) * t509;
t434 = Icges(5,4) * t484 + Icges(5,2) * t483 - Icges(5,6) * t511;
t433 = Icges(5,4) * t482 + Icges(5,2) * t481 - Icges(5,6) * t509;
t430 = rSges(6,1) * t480 - rSges(6,2) * t479 - rSges(6,3) * t511;
t429 = rSges(6,1) * t478 - rSges(6,2) * t477 - rSges(6,3) * t509;
t428 = Icges(6,1) * t480 - Icges(6,4) * t479 - Icges(6,5) * t511;
t427 = Icges(6,1) * t478 - Icges(6,4) * t477 - Icges(6,5) * t509;
t426 = Icges(6,4) * t480 - Icges(6,2) * t479 - Icges(6,6) * t511;
t425 = Icges(6,4) * t478 - Icges(6,2) * t477 - Icges(6,6) * t509;
t422 = rSges(7,1) * t476 + rSges(7,2) * t475 + rSges(7,3) * t513;
t421 = Icges(7,1) * t476 + Icges(7,4) * t475 + Icges(7,5) * t513;
t420 = Icges(7,4) * t476 + Icges(7,2) * t475 + Icges(7,6) * t513;
t419 = Icges(7,5) * t476 + Icges(7,6) * t475 + Icges(7,3) * t513;
t414 = t464 * t545 + (t556 * t579 - t590) * t551 + t593;
t413 = (-t463 - t507) * t545 + t579 * t594 + t577;
t412 = rSges(7,1) * t447 + rSges(7,2) * t446 + rSges(7,3) * t479;
t411 = rSges(7,1) * t445 + rSges(7,2) * t444 + rSges(7,3) * t477;
t410 = Icges(7,1) * t447 + Icges(7,4) * t446 + Icges(7,5) * t479;
t409 = Icges(7,1) * t445 + Icges(7,4) * t444 + Icges(7,5) * t477;
t408 = Icges(7,4) * t447 + Icges(7,2) * t446 + Icges(7,6) * t479;
t407 = Icges(7,4) * t445 + Icges(7,2) * t444 + Icges(7,6) * t477;
t406 = Icges(7,5) * t447 + Icges(7,6) * t446 + Icges(7,3) * t479;
t405 = Icges(7,5) * t445 + Icges(7,6) * t444 + Icges(7,3) * t477;
t404 = (t463 * t556 + t464 * t559) * t591 + t586;
t403 = t438 * t517 - t470 * t485 + t565;
t402 = -t437 * t517 + t470 * t486 + t563;
t401 = t437 * t485 - t438 * t486 + t569;
t400 = t430 * t517 + (-t449 - t453) * t485 + t562;
t399 = t453 * t486 + (-t417 - t429) * t517 + t561;
t398 = t429 * t485 + (-t418 - t430) * t486 + t566;
t397 = t412 * t474 - t422 * t442 + t441 * t517 + (-t449 - t473) * t485 + t562;
t396 = -t411 * t474 + t422 * t443 + t473 * t486 + (-t417 - t440) * t517 + t561;
t395 = t411 * t442 - t412 * t443 + t440 * t485 + (-t418 - t441) * t486 + t566;
t1 = m(4) * (t404 ^ 2 + t413 ^ 2 + t414 ^ 2) / 0.2e1 + m(3) * (t448 ^ 2 + t454 ^ 2 + t455 ^ 2) / 0.2e1 + m(7) * (t395 ^ 2 + t396 ^ 2 + t397 ^ 2) / 0.2e1 + m(5) * (t401 ^ 2 + t402 ^ 2 + t403 ^ 2) / 0.2e1 + m(6) * (t398 ^ 2 + t399 ^ 2 + t400 ^ 2) / 0.2e1 + t442 * ((t406 * t479 + t408 * t446 + t410 * t447) * t442 + (t405 * t479 + t407 * t446 + t409 * t447) * t443 + (t419 * t479 + t420 * t446 + t421 * t447) * t474) / 0.2e1 + t443 * ((t406 * t477 + t408 * t444 + t410 * t445) * t442 + (t405 * t477 + t407 * t444 + t409 * t445) * t443 + (t419 * t477 + t420 * t444 + t421 * t445) * t474) / 0.2e1 + t474 * ((t406 * t513 + t408 * t475 + t410 * t476) * t442 + (t405 * t513 + t407 * t475 + t409 * t476) * t443 + (t419 * t513 + t420 * t475 + t421 * t476) * t474) / 0.2e1 + ((-t451 * t479 + t452 * t480 + t468 * t483 + t469 * t484 - t610 * t511) * t517 + (-t425 * t479 + t427 * t480 + t433 * t483 + t435 * t484 - t612 * t511) * t486 + (-t426 * t479 + t428 * t480 + t434 * t483 + t436 * t484 - t611 * t511) * t485) * t485 / 0.2e1 + ((-t451 * t477 + t452 * t478 + t468 * t481 + t469 * t482 - t610 * t509) * t517 + (-t425 * t477 + t427 * t478 + t433 * t481 + t435 * t482 - t612 * t509) * t486 + (-t426 * t477 + t428 * t478 + t434 * t481 + t436 * t482 - t611 * t509) * t485) * t486 / 0.2e1 + ((-t451 * t513 + t514 * t452 + t468 * t515 + t469 * t516 - t610 * t525) * t517 + (-t425 * t513 + t427 * t514 + t433 * t515 + t435 * t516 - t612 * t525) * t486 + (-t426 * t513 + t428 * t514 + t434 * t515 + t436 * t516 - t611 * t525) * t485) * t517 / 0.2e1 + (((t457 * t598 + t525 * t459 - t526 * t461) * t556 - (t456 * t598 + t525 * t458 - t526 * t460) * t559) * t591 + (t598 * t495 + (t497 * t602 + t499 * t555) * t551) * t544 - (t598 * t494 + (t496 * t602 + t498 * t555) * t551) * t585 + (t525 * t491 - t526 * t492 + (t522 * t602 + t523 * t555) * t551 + t609 * t598) * t545) * t545 / 0.2e1 + (m(2) * (t539 ^ 2 + t540 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t458 * t511 - t460 * t512 - t496 * t531 - t498 * t532) * t559 + (t511 * t459 + t512 * t461 + t531 * t497 + t532 * t499 - t613) * t556) * t591 + (t491 * t511 + t492 * t512 + t522 * t531 + t523 * t532 + t609 * t595) * t545) * t544 / 0.2e1 - (((-t509 * t458 - t510 * t460 - t529 * t496 - t530 * t498 + t613) * t559 + (t459 * t509 + t461 * t510 + t497 * t529 + t499 * t530) * t556) * t591 + (t491 * t509 + t492 * t510 + t522 * t529 + t523 * t530 - t609 * t594) * t545) * t585 / 0.2e1;
T  = t1;
