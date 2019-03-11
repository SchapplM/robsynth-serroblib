% Calculate kinetic energy for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR1_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR1_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR1_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:00:10
% EndTime: 2019-03-08 23:00:13
% DurationCPUTime: 3.05s
% Computational Cost: add. (3394->376), mult. (5383->572), div. (0->0), fcn. (6386->14), ass. (0->168)
t579 = Icges(5,3) + Icges(6,3);
t529 = sin(pkin(11));
t531 = cos(pkin(11));
t538 = cos(qJ(2));
t532 = cos(pkin(6));
t535 = sin(qJ(2));
t561 = t532 * t535;
t508 = t529 * t538 + t531 * t561;
t528 = qJ(3) + qJ(4);
t551 = pkin(12) + t528;
t521 = sin(t551);
t530 = sin(pkin(6));
t549 = cos(t551);
t547 = t530 * t549;
t474 = t508 * t521 + t531 * t547;
t567 = t530 * t531;
t475 = t508 * t549 - t521 * t567;
t524 = sin(t528);
t525 = cos(t528);
t482 = -t508 * t524 - t525 * t567;
t483 = t508 * t525 - t524 * t567;
t560 = t532 * t538;
t507 = t529 * t535 - t531 * t560;
t578 = Icges(5,5) * t483 + Icges(6,5) * t475 + Icges(5,6) * t482 - Icges(6,6) * t474 + t507 * t579;
t510 = -t529 * t561 + t531 * t538;
t476 = t510 * t521 - t529 * t547;
t568 = t529 * t530;
t477 = t510 * t549 + t521 * t568;
t484 = -t510 * t524 + t525 * t568;
t485 = t510 * t525 + t524 * t568;
t509 = t529 * t560 + t531 * t535;
t577 = Icges(5,5) * t485 + Icges(6,5) * t477 + Icges(5,6) * t484 - Icges(6,6) * t476 + t509 * t579;
t565 = t530 * t535;
t492 = t521 * t565 - t532 * t549;
t493 = t532 * t521 + t535 * t547;
t499 = -t524 * t565 + t525 * t532;
t500 = t524 * t532 + t525 * t565;
t563 = t530 * t538;
t576 = Icges(5,5) * t500 + Icges(6,5) * t493 + Icges(5,6) * t499 - Icges(6,6) * t492 - t563 * t579;
t575 = qJD(2) ^ 2;
t537 = cos(qJ(3));
t570 = t537 * pkin(3);
t534 = sin(qJ(3));
t566 = t530 * t534;
t564 = t530 * t537;
t562 = t532 * t534;
t559 = pkin(4) * t525;
t557 = qJD(2) * t530;
t520 = t529 * t557;
t490 = qJD(3) * t509 + t520;
t523 = qJD(2) * t532;
t555 = t529 * t566;
t554 = t531 * t566;
t456 = qJD(4) * t509 + t490;
t553 = t531 * t557;
t478 = t508 * pkin(2) + t507 * pkin(8);
t479 = t510 * pkin(2) + t509 * pkin(8);
t552 = t478 * t520 + t479 * t553 + qJD(1);
t550 = pkin(4) * t524;
t491 = qJD(3) * t507 - t553;
t513 = (pkin(2) * t535 - pkin(8) * t538) * t530;
t548 = t479 * t523 - t513 * t520;
t457 = qJD(4) * t507 + t491;
t494 = t523 + (-qJD(3) - qJD(4)) * t563;
t424 = -pkin(3) * t554 + pkin(9) * t507 + t508 * t570;
t425 = pkin(3) * t555 + pkin(9) * t509 + t510 * t570;
t546 = t490 * t424 - t425 * t491 + t552;
t545 = (-t478 * t532 - t513 * t567) * qJD(2);
t472 = pkin(3) * t562 + (-pkin(9) * t538 + t535 * t570) * t530;
t514 = -qJD(3) * t563 + t523;
t544 = t514 * t425 - t472 * t490 + t548;
t399 = qJ(5) * t509 + t510 * t559 + t550 * t568;
t543 = qJD(5) * t507 + t494 * t399 + t544;
t398 = qJ(5) * t507 + t508 * t559 - t550 * t567;
t542 = -qJD(5) * t563 + t456 * t398 + t546;
t541 = -t424 * t514 + t491 * t472 + t545;
t440 = t550 * t532 + (-qJ(5) * t538 + t535 * t559) * t530;
t540 = qJD(5) * t509 + t457 * t440 + t541;
t536 = cos(qJ(6));
t533 = sin(qJ(6));
t512 = t535 * t564 + t562;
t511 = t532 * t537 - t534 * t565;
t498 = rSges(3,3) * t532 + (rSges(3,1) * t535 + rSges(3,2) * t538) * t530;
t497 = Icges(3,5) * t532 + (Icges(3,1) * t535 + Icges(3,4) * t538) * t530;
t496 = Icges(3,6) * t532 + (Icges(3,4) * t535 + Icges(3,2) * t538) * t530;
t495 = Icges(3,3) * t532 + (Icges(3,5) * t535 + Icges(3,6) * t538) * t530;
t489 = t510 * t537 + t555;
t488 = -t510 * t534 + t529 * t564;
t487 = t508 * t537 - t554;
t486 = -t508 * t534 - t531 * t564;
t481 = t493 * t536 - t533 * t563;
t480 = -t493 * t533 - t536 * t563;
t471 = rSges(4,1) * t512 + rSges(4,2) * t511 - rSges(4,3) * t563;
t470 = Icges(4,1) * t512 + Icges(4,4) * t511 - Icges(4,5) * t563;
t469 = Icges(4,4) * t512 + Icges(4,2) * t511 - Icges(4,6) * t563;
t468 = Icges(4,5) * t512 + Icges(4,6) * t511 - Icges(4,3) * t563;
t465 = rSges(3,1) * t510 - rSges(3,2) * t509 + rSges(3,3) * t568;
t464 = rSges(3,1) * t508 - rSges(3,2) * t507 - rSges(3,3) * t567;
t463 = Icges(3,1) * t510 - Icges(3,4) * t509 + Icges(3,5) * t568;
t462 = Icges(3,1) * t508 - Icges(3,4) * t507 - Icges(3,5) * t567;
t461 = Icges(3,4) * t510 - Icges(3,2) * t509 + Icges(3,6) * t568;
t460 = Icges(3,4) * t508 - Icges(3,2) * t507 - Icges(3,6) * t567;
t459 = Icges(3,5) * t510 - Icges(3,6) * t509 + Icges(3,3) * t568;
t458 = Icges(3,5) * t508 - Icges(3,6) * t507 - Icges(3,3) * t567;
t455 = qJD(6) * t492 + t494;
t454 = pkin(5) * t493 + pkin(10) * t492;
t453 = rSges(5,1) * t500 + rSges(5,2) * t499 - rSges(5,3) * t563;
t452 = Icges(5,1) * t500 + Icges(5,4) * t499 - Icges(5,5) * t563;
t451 = Icges(5,4) * t500 + Icges(5,2) * t499 - Icges(5,6) * t563;
t449 = rSges(6,1) * t493 - rSges(6,2) * t492 - rSges(6,3) * t563;
t448 = Icges(6,1) * t493 - Icges(6,4) * t492 - Icges(6,5) * t563;
t447 = Icges(6,4) * t493 - Icges(6,2) * t492 - Icges(6,6) * t563;
t445 = t477 * t536 + t509 * t533;
t444 = -t477 * t533 + t509 * t536;
t443 = t475 * t536 + t507 * t533;
t442 = -t475 * t533 + t507 * t536;
t439 = pkin(5) * t477 + pkin(10) * t476;
t438 = pkin(5) * t475 + pkin(10) * t474;
t437 = (-t464 * t532 - t498 * t567) * qJD(2);
t436 = (t465 * t532 - t498 * t568) * qJD(2);
t435 = rSges(4,1) * t489 + rSges(4,2) * t488 + rSges(4,3) * t509;
t434 = rSges(4,1) * t487 + rSges(4,2) * t486 + rSges(4,3) * t507;
t433 = Icges(4,1) * t489 + Icges(4,4) * t488 + Icges(4,5) * t509;
t432 = Icges(4,1) * t487 + Icges(4,4) * t486 + Icges(4,5) * t507;
t431 = Icges(4,4) * t489 + Icges(4,2) * t488 + Icges(4,6) * t509;
t430 = Icges(4,4) * t487 + Icges(4,2) * t486 + Icges(4,6) * t507;
t429 = Icges(4,5) * t489 + Icges(4,6) * t488 + Icges(4,3) * t509;
t428 = Icges(4,5) * t487 + Icges(4,6) * t486 + Icges(4,3) * t507;
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
t407 = rSges(7,1) * t481 + rSges(7,2) * t480 + rSges(7,3) * t492;
t406 = Icges(7,1) * t481 + Icges(7,4) * t480 + Icges(7,5) * t492;
t405 = Icges(7,4) * t481 + Icges(7,2) * t480 + Icges(7,6) * t492;
t404 = Icges(7,5) * t481 + Icges(7,6) * t480 + Icges(7,3) * t492;
t402 = qJD(1) + (t464 * t529 + t465 * t531) * t557;
t396 = rSges(7,1) * t445 + rSges(7,2) * t444 + rSges(7,3) * t476;
t395 = rSges(7,1) * t443 + rSges(7,2) * t442 + rSges(7,3) * t474;
t394 = Icges(7,1) * t445 + Icges(7,4) * t444 + Icges(7,5) * t476;
t393 = Icges(7,1) * t443 + Icges(7,4) * t442 + Icges(7,5) * t474;
t392 = Icges(7,4) * t445 + Icges(7,2) * t444 + Icges(7,6) * t476;
t391 = Icges(7,4) * t443 + Icges(7,2) * t442 + Icges(7,6) * t474;
t390 = Icges(7,5) * t445 + Icges(7,6) * t444 + Icges(7,3) * t476;
t389 = Icges(7,5) * t443 + Icges(7,6) * t442 + Icges(7,3) * t474;
t387 = -t434 * t514 + t471 * t491 + t545;
t386 = t435 * t514 - t471 * t490 + t548;
t385 = t434 * t490 - t435 * t491 + t552;
t384 = -t422 * t494 + t453 * t457 + t541;
t383 = t423 * t494 - t453 * t456 + t544;
t382 = t422 * t456 - t423 * t457 + t546;
t381 = t449 * t457 + (-t398 - t414) * t494 + t540;
t380 = t415 * t494 + (-t440 - t449) * t456 + t543;
t379 = t414 * t456 + (-t399 - t415) * t457 + t542;
t378 = -t395 * t455 + t407 * t427 + t454 * t457 + (-t398 - t438) * t494 + t540;
t377 = t396 * t455 - t407 * t426 + t439 * t494 + (-t440 - t454) * t456 + t543;
t376 = t395 * t426 - t396 * t427 + t438 * t456 + (-t399 - t439) * t457 + t542;
t1 = -t575 * ((-t459 * t567 - t461 * t507 + t463 * t508) * t568 - (-t458 * t567 - t460 * t507 + t462 * t508) * t567 + (-t495 * t567 - t496 * t507 + t497 * t508) * t532) * t567 / 0.2e1 + t491 * ((t429 * t507 + t431 * t486 + t433 * t487) * t490 + (t428 * t507 + t430 * t486 + t432 * t487) * t491 + (t468 * t507 + t469 * t486 + t470 * t487) * t514) / 0.2e1 + t514 * ((-t429 * t563 + t431 * t511 + t433 * t512) * t490 + (-t428 * t563 + t430 * t511 + t432 * t512) * t491 + (-t468 * t563 + t469 * t511 + t470 * t512) * t514) / 0.2e1 + t490 * ((t429 * t509 + t431 * t488 + t433 * t489) * t490 + (t428 * t509 + t430 * t488 + t432 * t489) * t491 + (t468 * t509 + t469 * t488 + t470 * t489) * t514) / 0.2e1 + m(7) * (t376 ^ 2 + t377 ^ 2 + t378 ^ 2) / 0.2e1 + m(6) * (t379 ^ 2 + t380 ^ 2 + t381 ^ 2) / 0.2e1 + m(5) * (t382 ^ 2 + t383 ^ 2 + t384 ^ 2) / 0.2e1 + m(4) * (t385 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + m(3) * (t402 ^ 2 + t436 ^ 2 + t437 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t427 * ((t390 * t474 + t392 * t442 + t394 * t443) * t426 + (t474 * t389 + t442 * t391 + t443 * t393) * t427 + (t404 * t474 + t405 * t442 + t406 * t443) * t455) / 0.2e1 + t455 * ((t390 * t492 + t392 * t480 + t394 * t481) * t426 + (t389 * t492 + t391 * t480 + t393 * t481) * t427 + (t492 * t404 + t480 * t405 + t481 * t406) * t455) / 0.2e1 + t426 * ((t476 * t390 + t444 * t392 + t445 * t394) * t426 + (t389 * t476 + t391 * t444 + t393 * t445) * t427 + (t404 * t476 + t405 * t444 + t406 * t445) * t455) / 0.2e1 + ((-t447 * t476 + t448 * t477 + t451 * t484 + t452 * t485 + t509 * t576) * t494 + (-t410 * t476 + t412 * t477 + t418 * t484 + t420 * t485 + t509 * t578) * t457 + (-t476 * t411 + t477 * t413 + t484 * t419 + t485 * t421 + t509 * t577) * t456) * t456 / 0.2e1 + ((-t447 * t474 + t448 * t475 + t451 * t482 + t452 * t483 + t507 * t576) * t494 + (-t474 * t410 + t475 * t412 + t482 * t418 + t483 * t420 + t507 * t578) * t457 + (-t411 * t474 + t413 * t475 + t419 * t482 + t421 * t483 + t507 * t577) * t456) * t457 / 0.2e1 + ((-t447 * t492 + t448 * t493 + t451 * t499 + t452 * t500 - t563 * t576) * t494 + (-t410 * t492 + t412 * t493 + t418 * t499 + t420 * t500 - t563 * t578) * t457 + (-t411 * t492 + t413 * t493 + t419 * t499 + t421 * t500 - t563 * t577) * t456) * t494 / 0.2e1 + (((t459 * t568 - t461 * t509 + t463 * t510) * t568 - (t458 * t568 - t460 * t509 + t462 * t510) * t567 + (t495 * t568 - t496 * t509 + t497 * t510) * t532) * t568 + t532 * (t532 ^ 2 * t495 + (((t461 * t538 + t463 * t535) * t529 - (t460 * t538 + t462 * t535) * t531) * t530 + (-t458 * t531 + t459 * t529 + t496 * t538 + t497 * t535) * t532) * t530)) * t575 / 0.2e1;
T  = t1;
