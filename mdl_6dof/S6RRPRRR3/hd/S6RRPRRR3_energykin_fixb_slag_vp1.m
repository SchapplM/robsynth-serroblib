% Calculate kinetic energy for
% S6RRPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR3_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:21:44
% EndTime: 2019-03-09 13:21:47
% DurationCPUTime: 3.08s
% Computational Cost: add. (2207->333), mult. (2259->530), div. (0->0), fcn. (2216->12), ass. (0->175)
t563 = Icges(3,3) + Icges(4,3);
t480 = qJ(2) + pkin(11);
t472 = sin(t480);
t473 = cos(t480);
t485 = sin(qJ(2));
t488 = cos(qJ(2));
t562 = Icges(3,5) * t488 + Icges(4,5) * t473 - Icges(3,6) * t485 - Icges(4,6) * t472;
t486 = sin(qJ(1));
t489 = cos(qJ(1));
t561 = t562 * t486 - t563 * t489;
t560 = t563 * t486 + t562 * t489;
t559 = Icges(3,5) * t485 + Icges(4,5) * t472 + Icges(3,6) * t488 + Icges(4,6) * t473;
t543 = Icges(4,4) * t472;
t449 = Icges(4,2) * t473 + t543;
t542 = Icges(4,4) * t473;
t450 = Icges(4,1) * t472 + t542;
t545 = Icges(3,4) * t485;
t457 = Icges(3,2) * t488 + t545;
t544 = Icges(3,4) * t488;
t458 = Icges(3,1) * t485 + t544;
t558 = -t449 * t472 + t450 * t473 - t457 * t485 + t458 * t488;
t507 = -Icges(4,2) * t472 + t542;
t414 = Icges(4,6) * t486 + t489 * t507;
t509 = Icges(4,1) * t473 - t543;
t416 = Icges(4,5) * t486 + t489 * t509;
t508 = -Icges(3,2) * t485 + t544;
t434 = Icges(3,6) * t486 + t489 * t508;
t510 = Icges(3,1) * t488 - t545;
t436 = Icges(3,5) * t486 + t489 * t510;
t557 = -t414 * t472 + t416 * t473 - t434 * t485 + t436 * t488;
t413 = -Icges(4,6) * t489 + t486 * t507;
t415 = -Icges(4,5) * t489 + t486 * t509;
t433 = -Icges(3,6) * t489 + t486 * t508;
t435 = -Icges(3,5) * t489 + t486 * t510;
t556 = t413 * t472 - t415 * t473 + t433 * t485 - t435 * t488;
t552 = pkin(2) * t485;
t549 = pkin(2) * t488;
t487 = cos(qJ(4));
t548 = t487 * pkin(4);
t541 = t472 * t486;
t540 = t472 * t489;
t539 = t473 * t486;
t538 = t473 * t489;
t482 = qJ(4) + qJ(5);
t476 = sin(t482);
t537 = t476 * t486;
t536 = t476 * t489;
t477 = cos(t482);
t535 = t477 * t486;
t534 = t477 * t489;
t484 = sin(qJ(4));
t533 = t484 * t486;
t532 = t484 * t489;
t531 = t486 * t487;
t530 = t487 * t489;
t409 = -qJ(3) * t489 + t486 * t549;
t410 = qJ(3) * t486 + t489 * t549;
t475 = qJD(2) * t486;
t525 = qJD(2) * t489;
t529 = t409 * t475 + t410 * t525;
t466 = pkin(1) * t486 - pkin(7) * t489;
t528 = -t409 - t466;
t527 = pkin(5) * t477;
t524 = qJD(4) * t472;
t446 = t489 * t524 + t475;
t523 = qJD(5) * t472;
t522 = qJD(6) * t472;
t521 = -qJD(4) - qJD(5);
t407 = t489 * t523 + t446;
t518 = pkin(5) * t476;
t447 = t486 * t524 - t525;
t516 = pkin(3) * t473 + pkin(8) * t472;
t439 = t516 * t486;
t440 = t516 * t489;
t517 = t439 * t475 + t440 * t525 + t529;
t408 = t486 * t523 + t447;
t454 = qJD(1) * (pkin(1) * t489 + pkin(7) * t486);
t515 = qJD(1) * t410 - qJD(3) * t489 + t454;
t514 = rSges(3,1) * t488 - rSges(3,2) * t485;
t513 = rSges(4,1) * t473 - rSges(4,2) * t472;
t512 = qJD(2) * (-rSges(4,1) * t472 - rSges(4,2) * t473 - t552);
t511 = qJD(2) * (-pkin(3) * t472 + pkin(8) * t473 - t552);
t497 = pkin(9) * t472 + t473 * t548;
t374 = -pkin(4) * t532 + t486 * t497;
t375 = pkin(4) * t533 + t489 * t497;
t498 = t446 * t374 - t375 * t447 + t517;
t496 = pkin(10) * t472 + t473 * t527;
t495 = qJD(1) * t440 + t486 * t511 + t515;
t474 = qJD(3) * t486;
t494 = t474 + (-t439 + t528) * qJD(1) + t489 * t511;
t391 = -pkin(9) * t473 + t472 * t548;
t467 = -qJD(4) * t473 + qJD(1);
t493 = t467 * t375 - t391 * t446 + t495;
t492 = -t374 * t467 + t447 * t391 + t494;
t478 = qJ(6) + t482;
t469 = cos(t478);
t468 = sin(t478);
t465 = rSges(2,1) * t489 - rSges(2,2) * t486;
t464 = rSges(2,1) * t486 + rSges(2,2) * t489;
t463 = rSges(3,1) * t485 + rSges(3,2) * t488;
t445 = t473 * t521 + qJD(1);
t444 = t473 * t530 + t533;
t443 = -t473 * t532 + t531;
t442 = t473 * t531 - t532;
t441 = -t473 * t533 - t530;
t438 = rSges(3,3) * t486 + t489 * t514;
t437 = -rSges(3,3) * t489 + t486 * t514;
t429 = t473 * t534 + t537;
t428 = -t473 * t536 + t535;
t427 = t473 * t535 - t536;
t426 = -t473 * t537 - t534;
t425 = qJD(1) + (-qJD(6) + t521) * t473;
t422 = rSges(4,3) * t486 + t489 * t513;
t421 = -rSges(4,3) * t489 + t486 * t513;
t420 = t468 * t486 + t469 * t538;
t419 = -t468 * t538 + t469 * t486;
t418 = -t468 * t489 + t469 * t539;
t417 = -t468 * t539 - t469 * t489;
t406 = -rSges(5,3) * t473 + (rSges(5,1) * t487 - rSges(5,2) * t484) * t472;
t404 = -Icges(5,5) * t473 + (Icges(5,1) * t487 - Icges(5,4) * t484) * t472;
t403 = -Icges(5,6) * t473 + (Icges(5,4) * t487 - Icges(5,2) * t484) * t472;
t402 = -Icges(5,3) * t473 + (Icges(5,5) * t487 - Icges(5,6) * t484) * t472;
t399 = -rSges(6,3) * t473 + (rSges(6,1) * t477 - rSges(6,2) * t476) * t472;
t398 = -Icges(6,5) * t473 + (Icges(6,1) * t477 - Icges(6,4) * t476) * t472;
t397 = -Icges(6,6) * t473 + (Icges(6,4) * t477 - Icges(6,2) * t476) * t472;
t396 = -Icges(6,3) * t473 + (Icges(6,5) * t477 - Icges(6,6) * t476) * t472;
t395 = -rSges(7,3) * t473 + (rSges(7,1) * t469 - rSges(7,2) * t468) * t472;
t394 = -Icges(7,5) * t473 + (Icges(7,1) * t469 - Icges(7,4) * t468) * t472;
t393 = -Icges(7,6) * t473 + (Icges(7,4) * t469 - Icges(7,2) * t468) * t472;
t392 = -Icges(7,3) * t473 + (Icges(7,5) * t469 - Icges(7,6) * t468) * t472;
t390 = t486 * t522 + t408;
t389 = t489 * t522 + t407;
t388 = -pkin(10) * t473 + t472 * t527;
t387 = qJD(1) * t438 - t463 * t475 + t454;
t386 = -t463 * t525 + (-t437 - t466) * qJD(1);
t385 = (t437 * t486 + t438 * t489) * qJD(2);
t384 = rSges(5,1) * t444 + rSges(5,2) * t443 + rSges(5,3) * t540;
t383 = rSges(5,1) * t442 + rSges(5,2) * t441 + rSges(5,3) * t541;
t381 = Icges(5,1) * t444 + Icges(5,4) * t443 + Icges(5,5) * t540;
t380 = Icges(5,1) * t442 + Icges(5,4) * t441 + Icges(5,5) * t541;
t379 = Icges(5,4) * t444 + Icges(5,2) * t443 + Icges(5,6) * t540;
t378 = Icges(5,4) * t442 + Icges(5,2) * t441 + Icges(5,6) * t541;
t377 = Icges(5,5) * t444 + Icges(5,6) * t443 + Icges(5,3) * t540;
t376 = Icges(5,5) * t442 + Icges(5,6) * t441 + Icges(5,3) * t541;
t373 = rSges(6,1) * t429 + rSges(6,2) * t428 + rSges(6,3) * t540;
t372 = rSges(6,1) * t427 + rSges(6,2) * t426 + rSges(6,3) * t541;
t371 = Icges(6,1) * t429 + Icges(6,4) * t428 + Icges(6,5) * t540;
t370 = Icges(6,1) * t427 + Icges(6,4) * t426 + Icges(6,5) * t541;
t369 = Icges(6,4) * t429 + Icges(6,2) * t428 + Icges(6,6) * t540;
t368 = Icges(6,4) * t427 + Icges(6,2) * t426 + Icges(6,6) * t541;
t367 = Icges(6,5) * t429 + Icges(6,6) * t428 + Icges(6,3) * t540;
t366 = Icges(6,5) * t427 + Icges(6,6) * t426 + Icges(6,3) * t541;
t364 = rSges(7,1) * t420 + rSges(7,2) * t419 + rSges(7,3) * t540;
t363 = rSges(7,1) * t418 + rSges(7,2) * t417 + rSges(7,3) * t541;
t362 = Icges(7,1) * t420 + Icges(7,4) * t419 + Icges(7,5) * t540;
t361 = Icges(7,1) * t418 + Icges(7,4) * t417 + Icges(7,5) * t541;
t360 = Icges(7,4) * t420 + Icges(7,2) * t419 + Icges(7,6) * t540;
t359 = Icges(7,4) * t418 + Icges(7,2) * t417 + Icges(7,6) * t541;
t358 = Icges(7,5) * t420 + Icges(7,6) * t419 + Icges(7,3) * t540;
t357 = Icges(7,5) * t418 + Icges(7,6) * t417 + Icges(7,3) * t541;
t355 = t486 * t518 + t489 * t496;
t354 = t486 * t496 - t489 * t518;
t353 = qJD(1) * t422 + t486 * t512 + t515;
t352 = t474 + t489 * t512 + (-t421 + t528) * qJD(1);
t351 = (t421 * t486 + t422 * t489) * qJD(2) + t529;
t350 = t384 * t467 - t406 * t446 + t495;
t349 = -t383 * t467 + t406 * t447 + t494;
t348 = t383 * t446 - t384 * t447 + t517;
t347 = t373 * t445 - t399 * t407 + t493;
t346 = -t372 * t445 + t399 * t408 + t492;
t345 = t372 * t407 - t373 * t408 + t498;
t344 = t355 * t445 + t364 * t425 - t388 * t407 - t389 * t395 + t493;
t343 = -t354 * t445 - t363 * t425 + t388 * t408 + t390 * t395 + t492;
t342 = t354 * t407 - t355 * t408 + t363 * t389 - t364 * t390 + t498;
t1 = t467 * ((-t376 * t447 - t377 * t446 - t402 * t467) * t473 + ((-t379 * t484 + t381 * t487) * t446 + (-t378 * t484 + t380 * t487) * t447 + (-t403 * t484 + t404 * t487) * t467) * t472) / 0.2e1 + m(7) * (t342 ^ 2 + t343 ^ 2 + t344 ^ 2) / 0.2e1 + m(6) * (t345 ^ 2 + t346 ^ 2 + t347 ^ 2) / 0.2e1 + m(5) * (t348 ^ 2 + t349 ^ 2 + t350 ^ 2) / 0.2e1 + m(3) * (t385 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + m(4) * (t351 ^ 2 + t352 ^ 2 + t353 ^ 2) / 0.2e1 + t389 * ((t358 * t540 + t419 * t360 + t420 * t362) * t389 + (t357 * t540 + t359 * t419 + t361 * t420) * t390 + (t392 * t540 + t393 * t419 + t394 * t420) * t425) / 0.2e1 + t390 * ((t358 * t541 + t360 * t417 + t362 * t418) * t389 + (t357 * t541 + t417 * t359 + t418 * t361) * t390 + (t392 * t541 + t393 * t417 + t394 * t418) * t425) / 0.2e1 + t425 * ((-t357 * t390 - t358 * t389 - t392 * t425) * t473 + ((-t360 * t468 + t362 * t469) * t389 + (-t359 * t468 + t361 * t469) * t390 + (-t393 * t468 + t394 * t469) * t425) * t472) / 0.2e1 + t407 * ((t367 * t540 + t428 * t369 + t429 * t371) * t407 + (t366 * t540 + t368 * t428 + t370 * t429) * t408 + (t396 * t540 + t397 * t428 + t398 * t429) * t445) / 0.2e1 + t408 * ((t367 * t541 + t369 * t426 + t371 * t427) * t407 + (t366 * t541 + t426 * t368 + t427 * t370) * t408 + (t396 * t541 + t397 * t426 + t398 * t427) * t445) / 0.2e1 + t445 * ((-t366 * t408 - t367 * t407 - t396 * t445) * t473 + ((-t369 * t476 + t371 * t477) * t407 + (-t368 * t476 + t370 * t477) * t408 + (-t397 * t476 + t398 * t477) * t445) * t472) / 0.2e1 + t446 * ((t377 * t540 + t443 * t379 + t444 * t381) * t446 + (t376 * t540 + t378 * t443 + t380 * t444) * t447 + (t402 * t540 + t403 * t443 + t404 * t444) * t467) / 0.2e1 + t447 * ((t377 * t541 + t379 * t441 + t381 * t442) * t446 + (t376 * t541 + t441 * t378 + t442 * t380) * t447 + (t402 * t541 + t403 * t441 + t404 * t442) * t467) / 0.2e1 + (m(2) * (t464 ^ 2 + t465 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((-t413 * t473 - t472 * t415 - t433 * t488 - t435 * t485) * t489 + (t414 * t473 + t416 * t472 + t434 * t488 + t436 * t485) * t486) * qJD(2) + (t473 * t449 + t472 * t450 + t488 * t457 + t485 * t458) * qJD(1)) * qJD(1) / 0.2e1 + ((t560 * t486 ^ 2 + (t556 * t489 + (t557 - t561) * t486) * t489) * qJD(2) + (t559 * t486 + t558 * t489) * qJD(1)) * t475 / 0.2e1 - ((t561 * t489 ^ 2 + (t557 * t486 + (t556 - t560) * t489) * t486) * qJD(2) + (t558 * t486 - t559 * t489) * qJD(1)) * t525 / 0.2e1;
T  = t1;
