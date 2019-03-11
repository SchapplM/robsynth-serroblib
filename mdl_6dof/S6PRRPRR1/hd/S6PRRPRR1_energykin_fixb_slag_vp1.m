% Calculate kinetic energy for
% S6PRRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR1_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:51:01
% EndTime: 2019-03-08 21:51:04
% DurationCPUTime: 3.00s
% Computational Cost: add. (3349->376), mult. (5293->571), div. (0->0), fcn. (6278->14), ass. (0->167)
t578 = Icges(4,3) + Icges(5,3);
t529 = sin(pkin(11));
t531 = cos(pkin(11));
t539 = cos(qJ(2));
t532 = cos(pkin(6));
t536 = sin(qJ(2));
t561 = t532 * t536;
t508 = t529 * t539 + t531 * t561;
t528 = qJ(3) + pkin(12);
t524 = sin(t528);
t525 = cos(t528);
t530 = sin(pkin(6));
t563 = t531 * t530;
t482 = -t508 * t524 - t525 * t563;
t483 = t508 * t525 - t524 * t563;
t535 = sin(qJ(3));
t538 = cos(qJ(3));
t486 = -t508 * t535 - t538 * t563;
t554 = t535 * t563;
t487 = t508 * t538 - t554;
t560 = t532 * t539;
t507 = t529 * t536 - t531 * t560;
t577 = Icges(4,5) * t487 + Icges(5,5) * t483 + Icges(4,6) * t486 + Icges(5,6) * t482 + t578 * t507;
t510 = -t529 * t561 + t531 * t539;
t567 = t529 * t530;
t484 = -t510 * t524 + t525 * t567;
t485 = t510 * t525 + t524 * t567;
t565 = t530 * t538;
t488 = -t510 * t535 + t529 * t565;
t555 = t535 * t567;
t489 = t510 * t538 + t555;
t509 = t529 * t560 + t531 * t536;
t576 = Icges(4,5) * t489 + Icges(5,5) * t485 + Icges(4,6) * t488 + Icges(5,6) * t484 + t578 * t509;
t566 = t530 * t536;
t498 = -t524 * t566 + t525 * t532;
t499 = t524 * t532 + t525 * t566;
t511 = t532 * t538 - t535 * t566;
t562 = t532 * t535;
t512 = t536 * t565 + t562;
t564 = t530 * t539;
t575 = Icges(4,5) * t512 + Icges(5,5) * t499 + Icges(4,6) * t511 + Icges(5,6) * t498 - t578 * t564;
t574 = qJD(2) ^ 2;
t569 = t538 * pkin(3);
t559 = pkin(4) * t525;
t557 = qJD(2) * t530;
t520 = t529 * t557;
t490 = qJD(3) * t509 + t520;
t523 = qJD(2) * t532;
t456 = qJD(5) * t509 + t490;
t553 = t531 * t557;
t478 = pkin(2) * t508 + pkin(8) * t507;
t479 = pkin(2) * t510 + pkin(8) * t509;
t552 = t478 * t520 + t479 * t553 + qJD(1);
t551 = qJ(5) + t528;
t550 = pkin(4) * t524;
t549 = cos(t551);
t491 = qJD(3) * t507 - t553;
t513 = (pkin(2) * t536 - pkin(8) * t539) * t530;
t548 = t479 * t523 - t513 * t520;
t457 = qJD(5) * t507 + t491;
t547 = t530 * t549;
t494 = t523 + (-qJD(3) - qJD(5)) * t564;
t425 = pkin(3) * t555 + qJ(4) * t509 + t510 * t569;
t514 = -qJD(3) * t564 + t523;
t546 = qJD(4) * t507 + t514 * t425 + t548;
t545 = (-t478 * t532 - t513 * t563) * qJD(2);
t424 = -pkin(3) * t554 + qJ(4) * t507 + t508 * t569;
t544 = -qJD(4) * t564 + t490 * t424 + t552;
t471 = pkin(3) * t562 + (-qJ(4) * t539 + t536 * t569) * t530;
t543 = qJD(4) * t509 + t491 * t471 + t545;
t399 = pkin(9) * t509 + t510 * t559 + t550 * t567;
t440 = t550 * t532 + (-pkin(9) * t539 + t536 * t559) * t530;
t542 = t514 * t399 + (-t440 - t471) * t490 + t546;
t398 = pkin(9) * t507 + t508 * t559 - t550 * t563;
t541 = t490 * t398 + (-t399 - t425) * t491 + t544;
t540 = t491 * t440 + (-t398 - t424) * t514 + t543;
t537 = cos(qJ(6));
t534 = sin(qJ(6));
t521 = sin(t551);
t500 = t532 * rSges(3,3) + (rSges(3,1) * t536 + rSges(3,2) * t539) * t530;
t497 = Icges(3,5) * t532 + (Icges(3,1) * t536 + Icges(3,4) * t539) * t530;
t496 = Icges(3,6) * t532 + (Icges(3,4) * t536 + Icges(3,2) * t539) * t530;
t495 = Icges(3,3) * t532 + (Icges(3,5) * t536 + Icges(3,6) * t539) * t530;
t493 = t532 * t521 + t536 * t547;
t492 = t521 * t566 - t532 * t549;
t481 = t493 * t537 - t534 * t564;
t480 = -t493 * t534 - t537 * t564;
t477 = t510 * t549 + t521 * t567;
t476 = t510 * t521 - t529 * t547;
t475 = t508 * t549 - t521 * t563;
t474 = t508 * t521 + t531 * t547;
t472 = t512 * rSges(4,1) + t511 * rSges(4,2) - rSges(4,3) * t564;
t470 = Icges(4,1) * t512 + Icges(4,4) * t511 - Icges(4,5) * t564;
t469 = Icges(4,4) * t512 + Icges(4,2) * t511 - Icges(4,6) * t564;
t465 = rSges(3,1) * t510 - rSges(3,2) * t509 + rSges(3,3) * t567;
t464 = rSges(3,1) * t508 - rSges(3,2) * t507 - rSges(3,3) * t563;
t463 = Icges(3,1) * t510 - Icges(3,4) * t509 + Icges(3,5) * t567;
t462 = Icges(3,1) * t508 - Icges(3,4) * t507 - Icges(3,5) * t563;
t461 = Icges(3,4) * t510 - Icges(3,2) * t509 + Icges(3,6) * t567;
t460 = Icges(3,4) * t508 - Icges(3,2) * t507 - Icges(3,6) * t563;
t459 = Icges(3,5) * t510 - Icges(3,6) * t509 + Icges(3,3) * t567;
t458 = Icges(3,5) * t508 - Icges(3,6) * t507 - Icges(3,3) * t563;
t455 = qJD(6) * t492 + t494;
t454 = pkin(5) * t493 + pkin(10) * t492;
t453 = t499 * rSges(5,1) + t498 * rSges(5,2) - rSges(5,3) * t564;
t452 = Icges(5,1) * t499 + Icges(5,4) * t498 - Icges(5,5) * t564;
t451 = Icges(5,4) * t499 + Icges(5,2) * t498 - Icges(5,6) * t564;
t449 = t493 * rSges(6,1) - t492 * rSges(6,2) - rSges(6,3) * t564;
t448 = Icges(6,1) * t493 - Icges(6,4) * t492 - Icges(6,5) * t564;
t447 = Icges(6,4) * t493 - Icges(6,2) * t492 - Icges(6,6) * t564;
t446 = Icges(6,5) * t493 - Icges(6,6) * t492 - Icges(6,3) * t564;
t445 = t477 * t537 + t509 * t534;
t444 = -t477 * t534 + t509 * t537;
t443 = t475 * t537 + t507 * t534;
t442 = -t475 * t534 + t507 * t537;
t439 = pkin(5) * t477 + pkin(10) * t476;
t438 = pkin(5) * t475 + pkin(10) * t474;
t437 = (-t464 * t532 - t500 * t563) * qJD(2);
t436 = (t465 * t532 - t500 * t567) * qJD(2);
t435 = rSges(4,1) * t489 + rSges(4,2) * t488 + rSges(4,3) * t509;
t434 = rSges(4,1) * t487 + rSges(4,2) * t486 + rSges(4,3) * t507;
t433 = Icges(4,1) * t489 + Icges(4,4) * t488 + Icges(4,5) * t509;
t432 = Icges(4,1) * t487 + Icges(4,4) * t486 + Icges(4,5) * t507;
t431 = Icges(4,4) * t489 + Icges(4,2) * t488 + Icges(4,6) * t509;
t430 = Icges(4,4) * t487 + Icges(4,2) * t486 + Icges(4,6) * t507;
t427 = qJD(6) * t474 + t457;
t426 = qJD(6) * t476 + t456;
t423 = rSges(5,1) * t485 + rSges(5,2) * t484 + rSges(5,3) * t509;
t422 = rSges(5,1) * t483 + rSges(5,2) * t482 + rSges(5,3) * t507;
t421 = Icges(5,1) * t485 + Icges(5,4) * t484 + Icges(5,5) * t509;
t420 = Icges(5,1) * t483 + Icges(5,4) * t482 + Icges(5,5) * t507;
t419 = Icges(5,4) * t485 + Icges(5,2) * t484 + Icges(5,6) * t509;
t418 = Icges(5,4) * t483 + Icges(5,2) * t482 + Icges(5,6) * t507;
t415 = rSges(6,1) * t477 - rSges(6,2) * t476 + rSges(6,3) * t509;
t414 = rSges(6,1) * t475 - rSges(6,2) * t474 + rSges(6,3) * t507;
t413 = Icges(6,1) * t477 - Icges(6,4) * t476 + Icges(6,5) * t509;
t412 = Icges(6,1) * t475 - Icges(6,4) * t474 + Icges(6,5) * t507;
t411 = Icges(6,4) * t477 - Icges(6,2) * t476 + Icges(6,6) * t509;
t410 = Icges(6,4) * t475 - Icges(6,2) * t474 + Icges(6,6) * t507;
t409 = Icges(6,5) * t477 - Icges(6,6) * t476 + Icges(6,3) * t509;
t408 = Icges(6,5) * t475 - Icges(6,6) * t474 + Icges(6,3) * t507;
t406 = rSges(7,1) * t481 + rSges(7,2) * t480 + rSges(7,3) * t492;
t405 = Icges(7,1) * t481 + Icges(7,4) * t480 + Icges(7,5) * t492;
t404 = Icges(7,4) * t481 + Icges(7,2) * t480 + Icges(7,6) * t492;
t403 = Icges(7,5) * t481 + Icges(7,6) * t480 + Icges(7,3) * t492;
t401 = qJD(1) + (t464 * t529 + t465 * t531) * t557;
t395 = rSges(7,1) * t445 + rSges(7,2) * t444 + rSges(7,3) * t476;
t394 = rSges(7,1) * t443 + rSges(7,2) * t442 + rSges(7,3) * t474;
t393 = Icges(7,1) * t445 + Icges(7,4) * t444 + Icges(7,5) * t476;
t392 = Icges(7,1) * t443 + Icges(7,4) * t442 + Icges(7,5) * t474;
t391 = Icges(7,4) * t445 + Icges(7,2) * t444 + Icges(7,6) * t476;
t390 = Icges(7,4) * t443 + Icges(7,2) * t442 + Icges(7,6) * t474;
t389 = Icges(7,5) * t445 + Icges(7,6) * t444 + Icges(7,3) * t476;
t388 = Icges(7,5) * t443 + Icges(7,6) * t442 + Icges(7,3) * t474;
t387 = -t434 * t514 + t472 * t491 + t545;
t386 = t435 * t514 - t472 * t490 + t548;
t385 = t434 * t490 - t435 * t491 + t552;
t384 = t453 * t491 + (-t422 - t424) * t514 + t543;
t383 = t423 * t514 + (-t453 - t471) * t490 + t546;
t382 = t490 * t422 + (-t423 - t425) * t491 + t544;
t381 = -t414 * t494 + t449 * t457 + t540;
t380 = t415 * t494 - t449 * t456 + t542;
t379 = t456 * t414 - t457 * t415 + t541;
t378 = -t394 * t455 + t406 * t427 - t438 * t494 + t454 * t457 + t540;
t377 = t395 * t455 - t406 * t426 + t439 * t494 - t454 * t456 + t542;
t376 = t426 * t394 - t427 * t395 + t456 * t438 - t457 * t439 + t541;
t1 = -t574 * ((-t459 * t563 - t461 * t507 + t463 * t508) * t567 - (-t458 * t563 - t460 * t507 + t462 * t508) * t563 + (-t495 * t563 - t496 * t507 + t497 * t508) * t532) * t563 / 0.2e1 + t427 * ((t389 * t474 + t391 * t442 + t393 * t443) * t426 + (t474 * t388 + t442 * t390 + t443 * t392) * t427 + (t403 * t474 + t404 * t442 + t405 * t443) * t455) / 0.2e1 + t455 * ((t389 * t492 + t391 * t480 + t393 * t481) * t426 + (t388 * t492 + t390 * t480 + t392 * t481) * t427 + (t492 * t403 + t480 * t404 + t481 * t405) * t455) / 0.2e1 + t426 * ((t476 * t389 + t444 * t391 + t445 * t393) * t426 + (t388 * t476 + t390 * t444 + t392 * t445) * t427 + (t403 * t476 + t404 * t444 + t405 * t445) * t455) / 0.2e1 + t456 * ((t509 * t409 - t476 * t411 + t477 * t413) * t456 + (t408 * t509 - t410 * t476 + t412 * t477) * t457 + (t446 * t509 - t447 * t476 + t448 * t477) * t494) / 0.2e1 + t457 * ((t409 * t507 - t411 * t474 + t413 * t475) * t456 + (t507 * t408 - t474 * t410 + t475 * t412) * t457 + (t446 * t507 - t447 * t474 + t448 * t475) * t494) / 0.2e1 + t494 * ((-t409 * t564 - t492 * t411 + t493 * t413) * t456 + (-t408 * t564 - t492 * t410 + t493 * t412) * t457 + (-t446 * t564 - t492 * t447 + t493 * t448) * t494) / 0.2e1 + m(7) * (t376 ^ 2 + t377 ^ 2 + t378 ^ 2) / 0.2e1 + m(6) * (t379 ^ 2 + t380 ^ 2 + t381 ^ 2) / 0.2e1 + m(4) * (t385 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + m(5) * (t382 ^ 2 + t383 ^ 2 + t384 ^ 2) / 0.2e1 + m(3) * (t401 ^ 2 + t436 ^ 2 + t437 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + ((t451 * t484 + t452 * t485 + t469 * t488 + t489 * t470 + t509 * t575) * t514 + (t418 * t484 + t420 * t485 + t430 * t488 + t432 * t489 + t509 * t577) * t491 + (t484 * t419 + t485 * t421 + t488 * t431 + t489 * t433 + t576 * t509) * t490) * t490 / 0.2e1 + ((t451 * t482 + t452 * t483 + t469 * t486 + t470 * t487 + t507 * t575) * t514 + (t482 * t418 + t483 * t420 + t486 * t430 + t487 * t432 + t577 * t507) * t491 + (t419 * t482 + t421 * t483 + t431 * t486 + t487 * t433 + t507 * t576) * t490) * t491 / 0.2e1 + ((t498 * t451 + t499 * t452 + t511 * t469 + t512 * t470 - t575 * t564) * t514 + (t498 * t418 + t499 * t420 + t511 * t430 + t512 * t432 - t564 * t577) * t491 + (t498 * t419 + t499 * t421 + t511 * t431 + t512 * t433 - t564 * t576) * t490) * t514 / 0.2e1 + (((t459 * t567 - t461 * t509 + t463 * t510) * t567 - (t458 * t567 - t460 * t509 + t462 * t510) * t563 + (t495 * t567 - t496 * t509 + t497 * t510) * t532) * t567 + t532 * (t532 ^ 2 * t495 + (((t461 * t539 + t463 * t536) * t529 - (t460 * t539 + t462 * t536) * t531) * t530 + (-t458 * t531 + t459 * t529 + t496 * t539 + t497 * t536) * t532) * t530)) * t574 / 0.2e1;
T  = t1;
