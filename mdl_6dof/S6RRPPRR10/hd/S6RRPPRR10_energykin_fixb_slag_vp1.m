% Calculate kinetic energy for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
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
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPRR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR10_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR10_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPPRR10_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:34:05
% EndTime: 2019-03-09 09:34:09
% DurationCPUTime: 4.17s
% Computational Cost: add. (1378->313), mult. (2128->481), div. (0->0), fcn. (2085->10), ass. (0->157)
t573 = Icges(3,4) + Icges(4,6);
t572 = Icges(3,1) + Icges(4,2);
t571 = Icges(3,2) + Icges(4,3);
t492 = cos(qJ(2));
t570 = t573 * t492;
t490 = sin(qJ(2));
t569 = t573 * t490;
t568 = -Icges(4,4) + Icges(3,5);
t567 = Icges(4,5) - Icges(3,6);
t566 = t571 * t490 - t570;
t565 = -t572 * t492 + t569;
t564 = Icges(4,1) + Icges(3,3);
t491 = sin(qJ(1));
t493 = cos(qJ(1));
t563 = -t491 * t566 + t493 * t567;
t562 = t491 * t567 + t493 * t566;
t561 = t491 * t568 - t565 * t493;
t560 = t565 * t491 + t493 * t568;
t559 = -t571 * t492 - t569;
t558 = t572 * t490 + t570;
t557 = t567 * t490 + t492 * t568;
t556 = t491 * t564 + t557 * t493;
t555 = t557 * t491 - t493 * t564;
t554 = t490 * t568 - t567 * t492;
t553 = t490 * t559 + t492 * t558;
t552 = t490 * t563 + t492 * t560;
t551 = t490 * t562 + t492 * t561;
t487 = sin(pkin(10));
t547 = pkin(4) * t487;
t488 = cos(pkin(10));
t545 = t488 * pkin(4);
t540 = t490 * t491;
t539 = t490 * t493;
t538 = t491 * t492;
t537 = t492 * t493;
t516 = pkin(2) * t492 + qJ(3) * t490;
t445 = t516 * t491;
t467 = pkin(1) * t491 - pkin(7) * t493;
t535 = -t445 - t467;
t486 = pkin(10) + qJ(5);
t479 = cos(t486);
t534 = pkin(5) * t479;
t530 = qJD(3) * t490;
t473 = t493 * t530;
t529 = qJD(4) * t492;
t533 = t493 * t529 + t473;
t483 = qJD(2) * t491;
t528 = qJD(5) * t492;
t448 = t493 * t528 + t483;
t531 = qJD(2) * t493;
t527 = qJD(6) * t492;
t474 = qJD(5) * t490 + qJD(1);
t446 = t516 * t493;
t455 = qJD(1) * (pkin(1) * t493 + pkin(7) * t491);
t526 = qJD(1) * t446 + t491 * t530 + t455;
t452 = -pkin(3) * t493 + qJ(4) * t538;
t525 = -t452 + t535;
t478 = sin(t486);
t522 = pkin(5) * t478;
t462 = pkin(2) * t490 - qJ(3) * t492;
t521 = -qJ(4) * t490 - t462;
t520 = qJD(2) * (rSges(4,2) * t490 + rSges(4,3) * t492 - t462);
t449 = t491 * t528 - t531;
t519 = -qJD(3) * t492 + t445 * t483 + t446 * t531;
t518 = rSges(3,1) * t492 - rSges(3,2) * t490;
t517 = -rSges(4,2) * t492 + rSges(4,3) * t490;
t451 = pkin(3) * t491 + qJ(4) * t537;
t515 = qJD(1) * t451 + t491 * t529 + t526;
t502 = qJD(2) * (-rSges(5,3) * t490 - (-rSges(5,1) * t487 - rSges(5,2) * t488) * t492 + t521);
t501 = qJD(2) * (-pkin(8) * t490 + t492 * t547 + t521);
t500 = pkin(8) * t492 + t490 * t547;
t499 = qJD(4) * t490 + t451 * t531 + t452 * t483 + t519;
t498 = pkin(9) * t492 + t490 * t522;
t394 = t491 * t545 + t493 * t500;
t395 = t491 * t500 - t493 * t545;
t497 = t394 * t531 + t395 * t483 + t499;
t496 = qJD(1) * t394 + t491 * t501 + t515;
t495 = (-t395 + t525) * qJD(1) + t493 * t501 + t533;
t480 = qJ(6) + t486;
t476 = cos(t480);
t475 = sin(t480);
t466 = rSges(2,1) * t493 - rSges(2,2) * t491;
t465 = rSges(2,1) * t491 + rSges(2,2) * t493;
t464 = rSges(3,1) * t490 + rSges(3,2) * t492;
t454 = qJD(6) * t490 + t474;
t444 = t487 * t540 - t488 * t493;
t443 = t487 * t493 + t488 * t540;
t442 = t487 * t539 + t488 * t491;
t441 = -t487 * t491 + t488 * t539;
t436 = -rSges(4,1) * t493 + t491 * t517;
t435 = rSges(4,1) * t491 + t493 * t517;
t434 = rSges(3,3) * t491 + t493 * t518;
t433 = -rSges(3,3) * t493 + t491 * t518;
t418 = t478 * t540 - t479 * t493;
t417 = t478 * t493 + t479 * t540;
t416 = t478 * t539 + t479 * t491;
t415 = -t478 * t491 + t479 * t539;
t413 = t491 * t527 + t449;
t412 = t493 * t527 + t448;
t411 = Icges(5,5) * t490 + (-Icges(5,1) * t487 - Icges(5,4) * t488) * t492;
t410 = Icges(5,6) * t490 + (-Icges(5,4) * t487 - Icges(5,2) * t488) * t492;
t409 = Icges(5,3) * t490 + (-Icges(5,5) * t487 - Icges(5,6) * t488) * t492;
t408 = t475 * t540 - t476 * t493;
t407 = t475 * t493 + t476 * t540;
t406 = t475 * t539 + t476 * t491;
t405 = -t475 * t491 + t476 * t539;
t404 = rSges(6,3) * t490 + (-rSges(6,1) * t478 - rSges(6,2) * t479) * t492;
t403 = Icges(6,5) * t490 + (-Icges(6,1) * t478 - Icges(6,4) * t479) * t492;
t402 = Icges(6,6) * t490 + (-Icges(6,4) * t478 - Icges(6,2) * t479) * t492;
t401 = Icges(6,3) * t490 + (-Icges(6,5) * t478 - Icges(6,6) * t479) * t492;
t400 = rSges(7,3) * t490 + (-rSges(7,1) * t475 - rSges(7,2) * t476) * t492;
t399 = Icges(7,5) * t490 + (-Icges(7,1) * t475 - Icges(7,4) * t476) * t492;
t398 = Icges(7,6) * t490 + (-Icges(7,4) * t475 - Icges(7,2) * t476) * t492;
t397 = Icges(7,3) * t490 + (-Icges(7,5) * t475 - Icges(7,6) * t476) * t492;
t396 = pkin(9) * t490 - t492 * t522;
t392 = rSges(5,1) * t444 + rSges(5,2) * t443 + rSges(5,3) * t538;
t391 = rSges(5,1) * t442 + rSges(5,2) * t441 + rSges(5,3) * t537;
t390 = Icges(5,1) * t444 + Icges(5,4) * t443 + Icges(5,5) * t538;
t389 = Icges(5,1) * t442 + Icges(5,4) * t441 + Icges(5,5) * t537;
t388 = Icges(5,4) * t444 + Icges(5,2) * t443 + Icges(5,6) * t538;
t387 = Icges(5,4) * t442 + Icges(5,2) * t441 + Icges(5,6) * t537;
t386 = Icges(5,5) * t444 + Icges(5,6) * t443 + Icges(5,3) * t538;
t385 = Icges(5,5) * t442 + Icges(5,6) * t441 + Icges(5,3) * t537;
t382 = qJD(1) * t434 - t464 * t483 + t455;
t381 = -t464 * t531 + (-t433 - t467) * qJD(1);
t380 = (t433 * t491 + t434 * t493) * qJD(2);
t379 = rSges(6,1) * t418 + rSges(6,2) * t417 + rSges(6,3) * t538;
t378 = rSges(6,1) * t416 + rSges(6,2) * t415 + rSges(6,3) * t537;
t377 = Icges(6,1) * t418 + Icges(6,4) * t417 + Icges(6,5) * t538;
t376 = Icges(6,1) * t416 + Icges(6,4) * t415 + Icges(6,5) * t537;
t375 = Icges(6,4) * t418 + Icges(6,2) * t417 + Icges(6,6) * t538;
t374 = Icges(6,4) * t416 + Icges(6,2) * t415 + Icges(6,6) * t537;
t373 = Icges(6,5) * t418 + Icges(6,6) * t417 + Icges(6,3) * t538;
t372 = Icges(6,5) * t416 + Icges(6,6) * t415 + Icges(6,3) * t537;
t371 = rSges(7,1) * t408 + rSges(7,2) * t407 + rSges(7,3) * t538;
t370 = rSges(7,1) * t406 + rSges(7,2) * t405 + rSges(7,3) * t537;
t369 = Icges(7,1) * t408 + Icges(7,4) * t407 + Icges(7,5) * t538;
t368 = Icges(7,1) * t406 + Icges(7,4) * t405 + Icges(7,5) * t537;
t367 = Icges(7,4) * t408 + Icges(7,2) * t407 + Icges(7,6) * t538;
t366 = Icges(7,4) * t406 + Icges(7,2) * t405 + Icges(7,6) * t537;
t365 = Icges(7,5) * t408 + Icges(7,6) * t407 + Icges(7,3) * t538;
t364 = Icges(7,5) * t406 + Icges(7,6) * t405 + Icges(7,3) * t537;
t363 = t491 * t498 - t493 * t534;
t362 = t491 * t534 + t493 * t498;
t361 = qJD(1) * t435 + t491 * t520 + t526;
t360 = t473 + t493 * t520 + (-t436 + t535) * qJD(1);
t359 = (t435 * t493 + t436 * t491) * qJD(2) + t519;
t358 = qJD(1) * t391 + t491 * t502 + t515;
t357 = t493 * t502 + (-t392 + t525) * qJD(1) + t533;
t356 = (t391 * t493 + t392 * t491) * qJD(2) + t499;
t355 = t378 * t474 - t404 * t448 + t496;
t354 = -t379 * t474 + t404 * t449 + t495;
t353 = -t378 * t449 + t379 * t448 + t497;
t352 = t362 * t474 + t370 * t454 - t396 * t448 - t400 * t412 + t496;
t351 = -t363 * t474 - t371 * t454 + t396 * t449 + t400 * t413 + t495;
t350 = -t362 * t449 + t363 * t448 - t370 * t413 + t371 * t412 + t497;
t1 = t413 * ((t364 * t538 + t366 * t407 + t368 * t408) * t412 + (t365 * t538 + t407 * t367 + t408 * t369) * t413 + (t397 * t538 + t398 * t407 + t399 * t408) * t454) / 0.2e1 + t454 * ((t364 * t412 + t365 * t413 + t397 * t454) * t490 + ((-t366 * t476 - t368 * t475) * t412 + (-t367 * t476 - t369 * t475) * t413 + (-t398 * t476 - t399 * t475) * t454) * t492) / 0.2e1 + t449 * ((t372 * t538 + t374 * t417 + t376 * t418) * t448 + (t373 * t538 + t375 * t417 + t377 * t418) * t449 + (t401 * t538 + t402 * t417 + t403 * t418) * t474) / 0.2e1 + t474 * ((t372 * t448 + t373 * t449 + t401 * t474) * t490 + ((-t374 * t479 - t376 * t478) * t448 + (-t375 * t479 - t377 * t478) * t449 + (-t402 * t479 - t403 * t478) * t474) * t492) / 0.2e1 + t448 * ((t372 * t537 + t415 * t374 + t416 * t376) * t448 + (t373 * t537 + t375 * t415 + t377 * t416) * t449 + (t401 * t537 + t402 * t415 + t403 * t416) * t474) / 0.2e1 + t412 * ((t364 * t537 + t405 * t366 + t406 * t368) * t412 + (t365 * t537 + t367 * t405 + t369 * t406) * t413 + (t397 * t537 + t398 * t405 + t399 * t406) * t454) / 0.2e1 + m(7) * (t350 ^ 2 + t351 ^ 2 + t352 ^ 2) / 0.2e1 + m(6) * (t353 ^ 2 + t354 ^ 2 + t355 ^ 2) / 0.2e1 + m(4) * (t359 ^ 2 + t360 ^ 2 + t361 ^ 2) / 0.2e1 + m(5) * (t356 ^ 2 + t357 ^ 2 + t358 ^ 2) / 0.2e1 + m(3) * (t380 ^ 2 + t381 ^ 2 + t382 ^ 2) / 0.2e1 + (m(2) * (t465 ^ 2 + t466 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((((t388 * t488 + t390 * t487 - t563) * t492 + (-t386 + t560) * t490) * t493 + ((-t387 * t488 - t389 * t487 - t562) * t492 + (t385 + t561) * t490) * t491) * qJD(2) + ((-t410 * t488 - t411 * t487 - t559) * t492 + (t409 + t558) * t490) * qJD(1)) * qJD(1) / 0.2e1 + (((-t386 * t537 - t388 * t441 - t390 * t442 + t552 * t493) * t493 + (t385 * t537 + t387 * t441 + t389 * t442 + (t551 - t555) * t493 + t556 * t491) * t491) * qJD(2) + (t409 * t537 + t410 * t441 + t411 * t442 + t491 * t554 + t493 * t553) * qJD(1)) * t483 / 0.2e1 - (((t385 * t538 + t387 * t443 + t389 * t444 + t551 * t491) * t491 + (-t386 * t538 - t388 * t443 - t390 * t444 + (t552 - t556) * t491 + t555 * t493) * t493) * qJD(2) + (t409 * t538 + t410 * t443 + t411 * t444 + t553 * t491 - t554 * t493) * qJD(1)) * t531 / 0.2e1;
T  = t1;
