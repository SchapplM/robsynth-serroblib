% Calculate kinetic energy for
% S6PRPRPR4
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
% Datum: 2019-03-08 19:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR4_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR4_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR4_energykin_fixb_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR4_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR4_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR4_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:39:01
% EndTime: 2019-03-08 19:39:04
% DurationCPUTime: 3.23s
% Computational Cost: add. (3184->376), mult. (5531->564), div. (0->0), fcn. (6690->14), ass. (0->167)
t573 = Icges(5,2) + Icges(6,3);
t519 = sin(pkin(10));
t523 = cos(pkin(10));
t528 = cos(qJ(2));
t524 = cos(pkin(6));
t527 = sin(qJ(2));
t555 = t524 * t527;
t501 = t519 * t528 + t523 * t555;
t547 = pkin(11) + qJ(4);
t514 = sin(t547);
t540 = cos(t547);
t520 = sin(pkin(6));
t558 = t520 * t523;
t478 = t501 * t540 - t514 * t558;
t554 = t524 * t528;
t500 = t519 * t527 - t523 * t554;
t517 = sin(pkin(12));
t521 = cos(pkin(12));
t443 = -t478 * t517 + t500 * t521;
t562 = t500 * t517;
t444 = t478 * t521 + t562;
t535 = t520 * t540;
t477 = t501 * t514 + t523 * t535;
t572 = -Icges(5,4) * t478 + Icges(6,5) * t444 - Icges(5,6) * t500 + Icges(6,6) * t443 + t573 * t477;
t503 = -t519 * t555 + t523 * t528;
t559 = t519 * t520;
t480 = t503 * t540 + t514 * t559;
t502 = t519 * t554 + t523 * t527;
t445 = -t480 * t517 + t502 * t521;
t561 = t502 * t517;
t446 = t480 * t521 + t561;
t479 = t503 * t514 - t519 * t535;
t571 = -Icges(5,4) * t480 + Icges(6,5) * t446 - Icges(5,6) * t502 + Icges(6,6) * t445 + t573 * t479;
t492 = t524 * t514 + t527 * t535;
t556 = t520 * t528;
t475 = -t492 * t517 - t521 * t556;
t546 = t517 * t556;
t476 = t492 * t521 - t546;
t557 = t520 * t527;
t491 = t514 * t557 - t524 * t540;
t570 = -Icges(5,4) * t492 + Icges(6,5) * t476 + Icges(5,6) * t556 + Icges(6,6) * t475 + t573 * t491;
t518 = sin(pkin(11));
t522 = cos(pkin(11));
t498 = -t518 * t557 + t522 * t524;
t560 = t518 * t524;
t499 = t522 * t557 + t560;
t451 = Icges(4,5) * t499 + Icges(4,6) * t498 - Icges(4,3) * t556;
t489 = Icges(3,6) * t524 + (Icges(3,4) * t527 + Icges(3,2) * t528) * t520;
t569 = t451 - t489;
t568 = qJD(2) ^ 2;
t564 = pkin(3) * t522;
t563 = pkin(5) * t521;
t471 = pkin(2) * t503 + qJ(3) * t502;
t512 = qJD(2) * t524;
t551 = qJD(3) * t500 + t471 * t512;
t550 = qJD(2) * t520;
t509 = t519 * t550;
t485 = qJD(4) * t502 + t509;
t549 = qJD(3) * t528;
t545 = t518 * t559;
t544 = t518 * t558;
t543 = t523 * t550;
t470 = pkin(2) * t501 + qJ(3) * t500;
t542 = t470 * t509 + t471 * t543 + qJD(1);
t504 = (pkin(2) * t527 - qJ(3) * t528) * t520;
t539 = (-t499 * rSges(4,1) - t498 * rSges(4,2) + rSges(4,3) * t556 - t504) * t520;
t538 = (-pkin(3) * t560 - (-pkin(8) * t528 + t527 * t564) * t520 - t504) * t520;
t486 = qJD(4) * t500 - t543;
t505 = -qJD(4) * t556 + t512;
t421 = -pkin(3) * t544 + pkin(8) * t500 + t501 * t564;
t422 = pkin(3) * t545 + pkin(8) * t502 + t503 * t564;
t534 = t421 * t509 + t422 * t543 - t520 * t549 + t542;
t533 = qJD(2) * t519 * t538 + t422 * t512 + t551;
t435 = pkin(4) * t478 + qJ(5) * t477;
t532 = qJD(5) * t491 + t485 * t435 + t534;
t436 = pkin(4) * t480 + qJ(5) * t479;
t531 = qJD(5) * t477 + t505 * t436 + t533;
t497 = qJD(3) * t502;
t530 = t497 + ((-t421 - t470) * t524 + t523 * t538) * qJD(2);
t454 = t492 * pkin(4) + t491 * qJ(5);
t529 = qJD(5) * t479 + t486 * t454 + t530;
t516 = pkin(12) + qJ(6);
t515 = cos(t516);
t513 = sin(t516);
t493 = t524 * rSges(3,3) + (rSges(3,1) * t527 + rSges(3,2) * t528) * t520;
t490 = Icges(3,5) * t524 + (Icges(3,1) * t527 + Icges(3,4) * t528) * t520;
t488 = Icges(3,3) * t524 + (Icges(3,5) * t527 + Icges(3,6) * t528) * t520;
t484 = t503 * t522 + t545;
t483 = -t503 * t518 + t522 * t559;
t482 = t501 * t522 - t544;
t481 = -t501 * t518 - t522 * t558;
t474 = qJD(6) * t491 + t505;
t469 = t492 * t515 - t513 * t556;
t468 = -t492 * t513 - t515 * t556;
t466 = rSges(3,1) * t503 - rSges(3,2) * t502 + rSges(3,3) * t559;
t465 = rSges(3,1) * t501 - rSges(3,2) * t500 - rSges(3,3) * t558;
t460 = Icges(3,1) * t503 - Icges(3,4) * t502 + Icges(3,5) * t559;
t459 = Icges(3,1) * t501 - Icges(3,4) * t500 - Icges(3,5) * t558;
t458 = Icges(3,4) * t503 - Icges(3,2) * t502 + Icges(3,6) * t559;
t457 = Icges(3,4) * t501 - Icges(3,2) * t500 - Icges(3,6) * t558;
t456 = Icges(3,5) * t503 - Icges(3,6) * t502 + Icges(3,3) * t559;
t455 = Icges(3,5) * t501 - Icges(3,6) * t500 - Icges(3,3) * t558;
t453 = Icges(4,1) * t499 + Icges(4,4) * t498 - Icges(4,5) * t556;
t452 = Icges(4,4) * t499 + Icges(4,2) * t498 - Icges(4,6) * t556;
t450 = t492 * rSges(5,1) - t491 * rSges(5,2) - rSges(5,3) * t556;
t449 = Icges(5,1) * t492 - Icges(5,4) * t491 - Icges(5,5) * t556;
t447 = Icges(5,5) * t492 - Icges(5,6) * t491 - Icges(5,3) * t556;
t442 = t480 * t515 + t502 * t513;
t441 = -t480 * t513 + t502 * t515;
t440 = t478 * t515 + t500 * t513;
t439 = -t478 * t513 + t500 * t515;
t438 = qJD(6) * t477 + t486;
t437 = qJD(6) * t479 + t485;
t433 = (-t465 * t524 - t493 * t558) * qJD(2);
t432 = (t466 * t524 - t493 * t559) * qJD(2);
t431 = rSges(4,1) * t484 + rSges(4,2) * t483 + rSges(4,3) * t502;
t430 = rSges(4,1) * t482 + rSges(4,2) * t481 + rSges(4,3) * t500;
t429 = Icges(4,1) * t484 + Icges(4,4) * t483 + Icges(4,5) * t502;
t428 = Icges(4,1) * t482 + Icges(4,4) * t481 + Icges(4,5) * t500;
t427 = Icges(4,4) * t484 + Icges(4,2) * t483 + Icges(4,6) * t502;
t426 = Icges(4,4) * t482 + Icges(4,2) * t481 + Icges(4,6) * t500;
t425 = Icges(4,5) * t484 + Icges(4,6) * t483 + Icges(4,3) * t502;
t424 = Icges(4,5) * t482 + Icges(4,6) * t481 + Icges(4,3) * t500;
t420 = rSges(5,1) * t480 - rSges(5,2) * t479 + rSges(5,3) * t502;
t419 = rSges(5,1) * t478 - rSges(5,2) * t477 + rSges(5,3) * t500;
t418 = Icges(5,1) * t480 - Icges(5,4) * t479 + Icges(5,5) * t502;
t417 = Icges(5,1) * t478 - Icges(5,4) * t477 + Icges(5,5) * t500;
t414 = Icges(5,5) * t480 - Icges(5,6) * t479 + Icges(5,3) * t502;
t413 = Icges(5,5) * t478 - Icges(5,6) * t477 + Icges(5,3) * t500;
t412 = rSges(6,1) * t476 + rSges(6,2) * t475 + rSges(6,3) * t491;
t410 = Icges(6,1) * t476 + Icges(6,4) * t475 + Icges(6,5) * t491;
t409 = Icges(6,4) * t476 + Icges(6,2) * t475 + Icges(6,6) * t491;
t404 = rSges(7,1) * t469 + rSges(7,2) * t468 + rSges(7,3) * t491;
t403 = Icges(7,1) * t469 + Icges(7,4) * t468 + Icges(7,5) * t491;
t402 = Icges(7,4) * t469 + Icges(7,2) * t468 + Icges(7,6) * t491;
t401 = Icges(7,5) * t469 + Icges(7,6) * t468 + Icges(7,3) * t491;
t400 = -pkin(5) * t546 + pkin(9) * t491 + t492 * t563;
t399 = qJD(1) + (t465 * t519 + t466 * t523) * t550;
t398 = rSges(6,1) * t446 + rSges(6,2) * t445 + rSges(6,3) * t479;
t397 = rSges(6,1) * t444 + rSges(6,2) * t443 + rSges(6,3) * t477;
t396 = Icges(6,1) * t446 + Icges(6,4) * t445 + Icges(6,5) * t479;
t395 = Icges(6,1) * t444 + Icges(6,4) * t443 + Icges(6,5) * t477;
t394 = Icges(6,4) * t446 + Icges(6,2) * t445 + Icges(6,6) * t479;
t393 = Icges(6,4) * t444 + Icges(6,2) * t443 + Icges(6,6) * t477;
t390 = rSges(7,1) * t442 + rSges(7,2) * t441 + rSges(7,3) * t479;
t389 = rSges(7,1) * t440 + rSges(7,2) * t439 + rSges(7,3) * t477;
t388 = Icges(7,1) * t442 + Icges(7,4) * t441 + Icges(7,5) * t479;
t387 = Icges(7,1) * t440 + Icges(7,4) * t439 + Icges(7,5) * t477;
t386 = Icges(7,4) * t442 + Icges(7,2) * t441 + Icges(7,6) * t479;
t385 = Icges(7,4) * t440 + Icges(7,2) * t439 + Icges(7,6) * t477;
t384 = Icges(7,5) * t442 + Icges(7,6) * t441 + Icges(7,3) * t479;
t383 = Icges(7,5) * t440 + Icges(7,6) * t439 + Icges(7,3) * t477;
t382 = pkin(5) * t561 + pkin(9) * t479 + t480 * t563;
t381 = pkin(5) * t562 + pkin(9) * t477 + t478 * t563;
t380 = t497 + ((-t430 - t470) * t524 + t523 * t539) * qJD(2);
t379 = (t431 * t524 + t519 * t539) * qJD(2) + t551;
t378 = (-t549 + (t430 * t519 + t431 * t523) * qJD(2)) * t520 + t542;
t377 = -t419 * t505 + t450 * t486 + t530;
t376 = t420 * t505 - t450 * t485 + t533;
t375 = t485 * t419 - t486 * t420 + t534;
t374 = t412 * t486 + (-t397 - t435) * t505 + t529;
t373 = t398 * t505 + (-t412 - t454) * t485 + t531;
t372 = t485 * t397 + (-t398 - t436) * t486 + t532;
t371 = -t389 * t474 + t400 * t486 + t404 * t438 + (-t381 - t435) * t505 + t529;
t370 = t382 * t505 + t390 * t474 - t404 * t437 + (-t400 - t454) * t485 + t531;
t369 = t485 * t381 + t437 * t389 - t438 * t390 + (-t382 - t436) * t486 + t532;
t1 = m(7) * (t369 ^ 2 + t370 ^ 2 + t371 ^ 2) / 0.2e1 + m(6) * (t372 ^ 2 + t373 ^ 2 + t374 ^ 2) / 0.2e1 + m(5) * (t375 ^ 2 + t376 ^ 2 + t377 ^ 2) / 0.2e1 + m(4) * (t378 ^ 2 + t379 ^ 2 + t380 ^ 2) / 0.2e1 + m(3) * (t399 ^ 2 + t432 ^ 2 + t433 ^ 2) / 0.2e1 + m(2) * qJD(1) ^ 2 / 0.2e1 + t438 * ((t384 * t477 + t386 * t439 + t388 * t440) * t437 + (t477 * t383 + t439 * t385 + t440 * t387) * t438 + (t401 * t477 + t402 * t439 + t403 * t440) * t474) / 0.2e1 + t474 * ((t384 * t491 + t386 * t468 + t388 * t469) * t437 + (t383 * t491 + t385 * t468 + t387 * t469) * t438 + (t491 * t401 + t468 * t402 + t469 * t403) * t474) / 0.2e1 + t437 * ((t479 * t384 + t441 * t386 + t442 * t388) * t437 + (t383 * t479 + t385 * t441 + t387 * t442) * t438 + (t401 * t479 + t402 * t441 + t403 * t442) * t474) / 0.2e1 + ((t409 * t445 + t410 * t446 + t447 * t502 + t449 * t480 + t479 * t570) * t505 + (t393 * t445 + t395 * t446 + t413 * t502 + t417 * t480 + t479 * t572) * t486 + (t394 * t445 + t396 * t446 + t414 * t502 + t418 * t480 + t571 * t479) * t485) * t485 / 0.2e1 + ((t409 * t443 + t410 * t444 + t447 * t500 + t449 * t478 + t477 * t570) * t505 + (t393 * t443 + t395 * t444 + t413 * t500 + t417 * t478 + t572 * t477) * t486 + (t394 * t443 + t396 * t444 + t414 * t500 + t418 * t478 + t477 * t571) * t485) * t486 / 0.2e1 + ((t409 * t475 + t410 * t476 - t447 * t556 + t492 * t449 + t570 * t491) * t505 + (t393 * t475 + t395 * t476 - t413 * t556 + t492 * t417 + t491 * t572) * t486 + (t394 * t475 + t396 * t476 - t414 * t556 + t492 * t418 + t491 * t571) * t485) * t505 / 0.2e1 - (((t425 * t500 + t427 * t481 + t429 * t482) * t519 - (t424 * t500 + t426 * t481 + t428 * t482) * t523) * t520 + (-t456 * t558 - t458 * t500 + t460 * t501) * t559 - (-t455 * t558 - t457 * t500 + t459 * t501) * t558 + (t452 * t481 + t453 * t482 - t488 * t558 + t490 * t501 + t500 * t569) * t524) * t568 * t558 / 0.2e1 + (((-t425 * t556 + t498 * t427 + t499 * t429) * t559 - (-t424 * t556 + t498 * t426 + t499 * t428) * t558 + ((t458 * t528 + t460 * t527) * t519 - (t457 * t528 + t459 * t527) * t523) * t520 ^ 2 + (-t451 * t556 + t498 * t452 + t499 * t453 + (-t455 * t523 + t456 * t519 + t489 * t528 + t490 * t527) * t520 + t524 * t488) * t524) * t524 + (((t425 * t502 + t427 * t483 + t429 * t484) * t519 - (t424 * t502 + t426 * t483 + t428 * t484) * t523) * t520 + (t456 * t559 - t458 * t502 + t460 * t503) * t559 - (t455 * t559 - t457 * t502 + t459 * t503) * t558 + (t452 * t483 + t453 * t484 + t488 * t559 + t490 * t503 + t502 * t569) * t524) * t559) * t568 / 0.2e1;
T  = t1;
