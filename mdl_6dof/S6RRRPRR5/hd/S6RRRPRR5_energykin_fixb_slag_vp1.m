% Calculate kinetic energy for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR5_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR5_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR5_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:20:35
% EndTime: 2019-03-09 18:20:38
% DurationCPUTime: 3.73s
% Computational Cost: add. (1594->290), mult. (1914->456), div. (0->0), fcn. (1815->10), ass. (0->163)
t568 = Icges(4,4) + Icges(5,6);
t567 = Icges(4,1) + Icges(5,2);
t566 = -Icges(4,2) - Icges(5,3);
t477 = qJ(2) + qJ(3);
t475 = cos(t477);
t565 = t568 * t475;
t473 = sin(t477);
t564 = t568 * t473;
t563 = Icges(5,4) - Icges(4,5);
t562 = Icges(5,5) - Icges(4,6);
t561 = t473 * t566 + t565;
t560 = t475 * t567 - t564;
t559 = Icges(5,1) + Icges(4,3);
t480 = sin(qJ(1));
t483 = cos(qJ(1));
t558 = t561 * t480 + t562 * t483;
t557 = -t562 * t480 + t561 * t483;
t556 = t560 * t480 + t563 * t483;
t555 = -t563 * t480 + t560 * t483;
t554 = t475 * t566 - t564;
t553 = t473 * t567 + t565;
t552 = t562 * t473 - t563 * t475;
t471 = qJD(2) * t480;
t453 = qJD(3) * t480 + t471;
t454 = (-qJD(2) - qJD(3)) * t483;
t551 = (-t558 * t473 + t556 * t475) * t454 + (-t557 * t473 + t555 * t475) * t453 + (t554 * t473 + t553 * t475) * qJD(1);
t550 = (t552 * t480 - t559 * t483) * t454 + (t559 * t480 + t552 * t483) * t453 + (-t563 * t473 - t562 * t475) * qJD(1);
t478 = sin(qJ(5));
t546 = pkin(5) * t478;
t545 = pkin(9) * t473;
t482 = cos(qJ(2));
t543 = pkin(2) * t482;
t481 = cos(qJ(5));
t542 = pkin(5) * t481;
t479 = sin(qJ(2));
t539 = Icges(3,4) * t479;
t538 = Icges(3,4) * t482;
t476 = qJ(5) + qJ(6);
t472 = sin(t476);
t533 = t472 * t480;
t532 = t472 * t483;
t474 = cos(t476);
t531 = t474 * t480;
t530 = t474 * t483;
t529 = t475 * t480;
t528 = t475 * t483;
t527 = t478 * t480;
t526 = t478 * t483;
t525 = t480 * t481;
t524 = t481 * t483;
t398 = -pkin(8) * t483 + t480 * t543;
t399 = pkin(8) * t480 + t483 * t543;
t521 = qJD(2) * t483;
t523 = t398 * t471 + t399 * t521;
t465 = pkin(1) * t480 - pkin(7) * t483;
t522 = -t398 - t465;
t520 = qJD(4) * t473;
t519 = qJD(5) * t475;
t518 = qJD(6) * t475;
t466 = qJD(5) * t473 + qJD(1);
t517 = pkin(2) * qJD(2) * t479;
t511 = pkin(3) * t475 + qJ(4) * t473;
t432 = t511 * t480;
t516 = -t432 + t522;
t428 = t483 * t519 + t453;
t515 = t483 * t517;
t514 = rSges(3,1) * t482 - rSges(3,2) * t479;
t513 = rSges(4,1) * t475 - rSges(4,2) * t473;
t512 = -rSges(5,2) * t475 + rSges(5,3) * t473;
t510 = Icges(3,1) * t482 - t539;
t508 = -Icges(3,2) * t479 + t538;
t505 = Icges(3,5) * t482 - Icges(3,6) * t479;
t424 = -Icges(3,6) * t483 + t480 * t508;
t426 = -Icges(3,5) * t483 + t480 * t510;
t501 = t424 * t479 - t426 * t482;
t425 = Icges(3,6) * t480 + t483 * t508;
t427 = Icges(3,5) * t480 + t483 * t510;
t500 = -t425 * t479 + t427 * t482;
t456 = Icges(3,2) * t482 + t539;
t457 = Icges(3,1) * t479 + t538;
t499 = -t456 * t479 + t457 * t482;
t429 = t480 * t519 + t454;
t498 = -qJD(4) * t475 + t453 * t432 + t523;
t452 = qJD(1) * (pkin(1) * t483 + pkin(7) * t480);
t497 = qJD(1) * t399 - t480 * t517 + t452;
t449 = pkin(3) * t473 - qJ(4) * t475;
t496 = t454 * t449 + t483 * t520 - t515;
t495 = pkin(10) * t475 + t473 * t546;
t433 = t511 * t483;
t492 = qJD(1) * t433 + t480 * t520 + t497;
t441 = pkin(4) * t480 + pkin(9) * t528;
t442 = -pkin(4) * t483 + pkin(9) * t529;
t491 = t453 * t442 + (-t433 - t441) * t454 + t498;
t490 = t454 * t545 + (-t442 + t516) * qJD(1) + t496;
t489 = qJD(1) * t441 + (-t449 - t545) * t453 + t492;
t460 = rSges(2,1) * t483 - rSges(2,2) * t480;
t459 = rSges(2,1) * t480 + rSges(2,2) * t483;
t458 = rSges(3,1) * t479 + rSges(3,2) * t482;
t455 = Icges(3,5) * t479 + Icges(3,6) * t482;
t451 = rSges(4,1) * t473 + rSges(4,2) * t475;
t450 = -rSges(5,2) * t473 - rSges(5,3) * t475;
t440 = qJD(6) * t473 + t466;
t437 = t473 * t527 - t524;
t436 = t473 * t525 + t526;
t435 = t473 * t526 + t525;
t434 = t473 * t524 - t527;
t431 = rSges(3,3) * t480 + t483 * t514;
t430 = -rSges(3,3) * t483 + t480 * t514;
t423 = Icges(3,3) * t480 + t483 * t505;
t422 = -Icges(3,3) * t483 + t480 * t505;
t420 = t473 * t533 - t530;
t419 = t473 * t531 + t532;
t418 = t473 * t532 + t531;
t417 = t473 * t530 - t533;
t416 = -rSges(5,1) * t483 + t480 * t512;
t415 = rSges(5,1) * t480 + t483 * t512;
t414 = rSges(4,3) * t480 + t483 * t513;
t413 = -rSges(4,3) * t483 + t480 * t513;
t412 = pkin(10) * t473 - t475 * t546;
t395 = rSges(6,3) * t473 + (-rSges(6,1) * t478 - rSges(6,2) * t481) * t475;
t394 = Icges(6,5) * t473 + (-Icges(6,1) * t478 - Icges(6,4) * t481) * t475;
t393 = Icges(6,6) * t473 + (-Icges(6,4) * t478 - Icges(6,2) * t481) * t475;
t392 = Icges(6,3) * t473 + (-Icges(6,5) * t478 - Icges(6,6) * t481) * t475;
t388 = rSges(7,3) * t473 + (-rSges(7,1) * t472 - rSges(7,2) * t474) * t475;
t387 = t480 * t518 + t429;
t386 = t483 * t518 + t428;
t385 = Icges(7,5) * t473 + (-Icges(7,1) * t472 - Icges(7,4) * t474) * t475;
t384 = Icges(7,6) * t473 + (-Icges(7,4) * t472 - Icges(7,2) * t474) * t475;
t383 = Icges(7,3) * t473 + (-Icges(7,5) * t472 - Icges(7,6) * t474) * t475;
t381 = t480 * t495 - t483 * t542;
t380 = t480 * t542 + t483 * t495;
t379 = qJD(1) * t431 - t458 * t471 + t452;
t378 = -t458 * t521 + (-t430 - t465) * qJD(1);
t377 = rSges(6,1) * t437 + rSges(6,2) * t436 + rSges(6,3) * t529;
t376 = rSges(6,1) * t435 + rSges(6,2) * t434 + rSges(6,3) * t528;
t375 = Icges(6,1) * t437 + Icges(6,4) * t436 + Icges(6,5) * t529;
t374 = Icges(6,1) * t435 + Icges(6,4) * t434 + Icges(6,5) * t528;
t373 = Icges(6,4) * t437 + Icges(6,2) * t436 + Icges(6,6) * t529;
t372 = Icges(6,4) * t435 + Icges(6,2) * t434 + Icges(6,6) * t528;
t371 = Icges(6,5) * t437 + Icges(6,6) * t436 + Icges(6,3) * t529;
t370 = Icges(6,5) * t435 + Icges(6,6) * t434 + Icges(6,3) * t528;
t369 = (t430 * t480 + t431 * t483) * qJD(2);
t368 = rSges(7,1) * t420 + rSges(7,2) * t419 + rSges(7,3) * t529;
t367 = rSges(7,1) * t418 + rSges(7,2) * t417 + rSges(7,3) * t528;
t366 = Icges(7,1) * t420 + Icges(7,4) * t419 + Icges(7,5) * t529;
t365 = Icges(7,1) * t418 + Icges(7,4) * t417 + Icges(7,5) * t528;
t364 = Icges(7,4) * t420 + Icges(7,2) * t419 + Icges(7,6) * t529;
t363 = Icges(7,4) * t418 + Icges(7,2) * t417 + Icges(7,6) * t528;
t362 = Icges(7,5) * t420 + Icges(7,6) * t419 + Icges(7,3) * t529;
t361 = Icges(7,5) * t418 + Icges(7,6) * t417 + Icges(7,3) * t528;
t360 = qJD(1) * t414 - t451 * t453 + t497;
t359 = -t515 + t451 * t454 + (-t413 + t522) * qJD(1);
t358 = t413 * t453 - t414 * t454 + t523;
t357 = qJD(1) * t415 + (-t449 - t450) * t453 + t492;
t356 = t450 * t454 + (-t416 + t516) * qJD(1) + t496;
t355 = t416 * t453 + (-t415 - t433) * t454 + t498;
t354 = t376 * t466 - t395 * t428 + t489;
t353 = -t377 * t466 + t395 * t429 + t490;
t352 = -t376 * t429 + t377 * t428 + t491;
t351 = t367 * t440 + t380 * t466 - t386 * t388 - t412 * t428 + t489;
t350 = -t368 * t440 - t381 * t466 + t387 * t388 + t412 * t429 + t490;
t349 = -t367 * t387 + t368 * t386 - t380 * t429 + t381 * t428 + t491;
t1 = -((-t483 * t455 + t480 * t499) * qJD(1) + (t483 ^ 2 * t422 + (t500 * t480 + (-t423 + t501) * t483) * t480) * qJD(2)) * t521 / 0.2e1 + ((t480 * t455 + t483 * t499) * qJD(1) + (t480 ^ 2 * t423 + (t501 * t483 + (-t422 + t500) * t480) * t483) * qJD(2)) * t471 / 0.2e1 + t387 * ((t361 * t529 + t363 * t419 + t365 * t420) * t386 + (t362 * t529 + t419 * t364 + t420 * t366) * t387 + (t383 * t529 + t384 * t419 + t385 * t420) * t440) / 0.2e1 + t440 * ((t361 * t386 + t362 * t387 + t383 * t440) * t473 + ((-t363 * t474 - t365 * t472) * t386 + (-t364 * t474 - t366 * t472) * t387 + (-t384 * t474 - t385 * t472) * t440) * t475) / 0.2e1 + t386 * ((t361 * t528 + t417 * t363 + t418 * t365) * t386 + (t362 * t528 + t364 * t417 + t366 * t418) * t387 + (t383 * t528 + t384 * t417 + t385 * t418) * t440) / 0.2e1 + t428 * ((t370 * t528 + t434 * t372 + t435 * t374) * t428 + (t371 * t528 + t373 * t434 + t375 * t435) * t429 + (t392 * t528 + t393 * t434 + t394 * t435) * t466) / 0.2e1 + t429 * ((t370 * t529 + t372 * t436 + t374 * t437) * t428 + (t371 * t529 + t436 * t373 + t437 * t375) * t429 + (t392 * t529 + t393 * t436 + t394 * t437) * t466) / 0.2e1 + t466 * ((t370 * t428 + t371 * t429 + t392 * t466) * t473 + ((-t372 * t481 - t374 * t478) * t428 + (-t373 * t481 - t375 * t478) * t429 + (-t393 * t481 - t394 * t478) * t466) * t475) / 0.2e1 + m(7) * (t349 ^ 2 + t350 ^ 2 + t351 ^ 2) / 0.2e1 + m(6) * (t352 ^ 2 + t353 ^ 2 + t354 ^ 2) / 0.2e1 + m(5) * (t355 ^ 2 + t356 ^ 2 + t357 ^ 2) / 0.2e1 + m(4) * (t358 ^ 2 + t359 ^ 2 + t360 ^ 2) / 0.2e1 + m(3) * (t369 ^ 2 + t378 ^ 2 + t379 ^ 2) / 0.2e1 + (t550 * t480 + t551 * t483) * t453 / 0.2e1 + (t551 * t480 - t550 * t483) * t454 / 0.2e1 + (Icges(2,3) + m(2) * (t459 ^ 2 + t460 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + (((t425 * t482 + t427 * t479) * t480 - (t424 * t482 + t426 * t479) * t483) * qJD(2) + (t556 * t473 + t558 * t475) * t454 + (t555 * t473 + t557 * t475) * t453 + (t482 * t456 + t479 * t457 + t553 * t473 - t554 * t475) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
