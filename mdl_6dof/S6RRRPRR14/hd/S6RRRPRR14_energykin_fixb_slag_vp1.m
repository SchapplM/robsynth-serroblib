% Calculate kinetic energy for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR14_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR14_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR14_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR14_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR14_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:12:21
% EndTime: 2019-03-09 20:12:24
% DurationCPUTime: 2.96s
% Computational Cost: add. (2501->342), mult. (5795->518), div. (0->0), fcn. (7021->12), ass. (0->158)
t568 = Icges(4,1) + Icges(5,2);
t567 = Icges(5,1) + Icges(4,3);
t566 = -Icges(4,4) - Icges(5,6);
t565 = Icges(5,4) - Icges(4,5);
t564 = Icges(5,5) - Icges(4,6);
t563 = Icges(4,2) + Icges(5,3);
t514 = sin(qJ(2));
t515 = sin(qJ(1));
t517 = cos(qJ(2));
t518 = cos(qJ(1));
t544 = cos(pkin(6));
t530 = t518 * t544;
t491 = t514 * t515 - t517 * t530;
t492 = t514 * t530 + t515 * t517;
t512 = sin(pkin(6));
t538 = t512 * t518;
t451 = Icges(3,5) * t492 - Icges(3,6) * t491 - Icges(3,3) * t538;
t531 = t515 * t544;
t493 = t518 * t514 + t517 * t531;
t494 = -t514 * t531 + t518 * t517;
t540 = t512 * t515;
t452 = Icges(3,5) * t494 - Icges(3,6) * t493 + Icges(3,3) * t540;
t562 = (t451 * t518 - t452 * t515) * t512;
t548 = cos(qJ(3));
t534 = t512 * t548;
t547 = sin(qJ(3));
t474 = t492 * t547 + t518 * t534;
t533 = t512 * t547;
t475 = t492 * t548 - t518 * t533;
t561 = t563 * t474 + t566 * t475 + t564 * t491;
t476 = t494 * t547 - t515 * t534;
t477 = t494 * t548 + t515 * t533;
t560 = t563 * t476 + t566 * t477 + t564 * t493;
t559 = t564 * t474 - t565 * t475 + t567 * t491;
t558 = t564 * t476 - t565 * t477 + t567 * t493;
t557 = t566 * t474 + t568 * t475 - t565 * t491;
t556 = t566 * t476 + t568 * t477 - t565 * t493;
t489 = t514 * t533 - t544 * t548;
t490 = t514 * t534 + t544 * t547;
t539 = t512 * t517;
t555 = t563 * t489 + t566 * t490 - t564 * t539;
t554 = t566 * t489 + t568 * t490 + t565 * t539;
t553 = t564 * t489 - t565 * t490 - t567 * t539;
t516 = cos(qJ(5));
t546 = pkin(5) * t516;
t513 = sin(qJ(5));
t543 = t474 * t513;
t542 = t476 * t513;
t541 = t489 * t513;
t463 = pkin(2) * t492 + pkin(9) * t491;
t464 = pkin(2) * t494 + pkin(9) * t493;
t535 = qJD(2) * t512;
t504 = t515 * t535;
t532 = t518 * t535;
t537 = t463 * t504 + t464 * t532;
t478 = qJD(3) * t493 + t504;
t536 = qJD(1) * (pkin(1) * t515 - pkin(8) * t538);
t505 = qJD(2) * t544 + qJD(1);
t429 = qJD(5) * t477 + t478;
t427 = pkin(3) * t475 + qJ(4) * t474;
t529 = qJD(4) * t489 + t478 * t427 + t537;
t479 = qJD(3) * t491 - t532;
t430 = qJD(5) * t475 + t479;
t496 = -qJD(3) * t539 + t505;
t495 = (pkin(2) * t514 - pkin(9) * t517) * t512;
t497 = qJD(1) * (pkin(1) * t518 + pkin(8) * t540);
t527 = t505 * t464 - t495 * t504 + t497;
t465 = qJD(5) * t490 + t496;
t428 = pkin(3) * t477 + qJ(4) * t476;
t526 = qJD(4) * t474 + t496 * t428 + t527;
t441 = pkin(4) * t491 + pkin(10) * t475;
t442 = pkin(4) * t493 + pkin(10) * t477;
t525 = t478 * t441 + (-t428 - t442) * t479 + t529;
t524 = -t463 * t505 - t495 * t532 - t536;
t462 = pkin(3) * t490 + qJ(4) * t489;
t523 = qJD(4) * t476 + t479 * t462 + t524;
t480 = -pkin(4) * t539 + pkin(10) * t490;
t522 = t496 * t442 + (-t462 - t480) * t478 + t526;
t521 = t479 * t480 + (-t427 - t441) * t496 + t523;
t511 = qJ(5) + qJ(6);
t510 = cos(t511);
t509 = sin(t511);
t500 = rSges(2,1) * t518 - rSges(2,2) * t515;
t499 = rSges(2,1) * t515 + rSges(2,2) * t518;
t484 = t544 * rSges(3,3) + (rSges(3,1) * t514 + rSges(3,2) * t517) * t512;
t483 = Icges(3,5) * t544 + (Icges(3,1) * t514 + Icges(3,4) * t517) * t512;
t482 = Icges(3,6) * t544 + (Icges(3,4) * t514 + Icges(3,2) * t517) * t512;
t481 = Icges(3,3) * t544 + (Icges(3,5) * t514 + Icges(3,6) * t517) * t512;
t473 = -t516 * t539 + t541;
t472 = t489 * t516 + t513 * t539;
t467 = t489 * t509 - t510 * t539;
t466 = t489 * t510 + t509 * t539;
t459 = rSges(3,1) * t494 - rSges(3,2) * t493 + rSges(3,3) * t540;
t458 = rSges(3,1) * t492 - rSges(3,2) * t491 - rSges(3,3) * t538;
t456 = Icges(3,1) * t494 - Icges(3,4) * t493 + Icges(3,5) * t540;
t455 = Icges(3,1) * t492 - Icges(3,4) * t491 - Icges(3,5) * t538;
t454 = Icges(3,4) * t494 - Icges(3,2) * t493 + Icges(3,6) * t540;
t453 = Icges(3,4) * t492 - Icges(3,2) * t491 - Icges(3,6) * t538;
t450 = rSges(4,1) * t490 - rSges(4,2) * t489 - rSges(4,3) * t539;
t449 = -rSges(5,1) * t539 - rSges(5,2) * t490 + rSges(5,3) * t489;
t440 = qJD(6) * t490 + t465;
t439 = t493 * t516 + t542;
t438 = t476 * t516 - t493 * t513;
t437 = t491 * t516 + t543;
t436 = t474 * t516 - t491 * t513;
t434 = t476 * t509 + t493 * t510;
t433 = t476 * t510 - t493 * t509;
t432 = t474 * t509 + t491 * t510;
t431 = t474 * t510 - t491 * t509;
t423 = pkin(5) * t541 + pkin(11) * t490 - t539 * t546;
t421 = rSges(4,1) * t477 - rSges(4,2) * t476 + rSges(4,3) * t493;
t420 = rSges(4,1) * t475 - rSges(4,2) * t474 + rSges(4,3) * t491;
t419 = rSges(5,1) * t493 - rSges(5,2) * t477 + rSges(5,3) * t476;
t418 = rSges(5,1) * t491 - rSges(5,2) * t475 + rSges(5,3) * t474;
t405 = rSges(6,1) * t473 + rSges(6,2) * t472 + rSges(6,3) * t490;
t404 = Icges(6,1) * t473 + Icges(6,4) * t472 + Icges(6,5) * t490;
t403 = Icges(6,4) * t473 + Icges(6,2) * t472 + Icges(6,6) * t490;
t402 = Icges(6,5) * t473 + Icges(6,6) * t472 + Icges(6,3) * t490;
t401 = rSges(7,1) * t467 + rSges(7,2) * t466 + rSges(7,3) * t490;
t400 = qJD(6) * t475 + t430;
t399 = qJD(6) * t477 + t429;
t398 = Icges(7,1) * t467 + Icges(7,4) * t466 + Icges(7,5) * t490;
t397 = Icges(7,4) * t467 + Icges(7,2) * t466 + Icges(7,6) * t490;
t396 = Icges(7,5) * t467 + Icges(7,6) * t466 + Icges(7,3) * t490;
t394 = t459 * t505 - t484 * t504 + t497;
t393 = -t458 * t505 - t484 * t532 - t536;
t392 = (t458 * t515 + t459 * t518) * t535;
t391 = pkin(5) * t542 + pkin(11) * t477 + t493 * t546;
t390 = pkin(5) * t543 + pkin(11) * t475 + t491 * t546;
t389 = rSges(6,1) * t439 + rSges(6,2) * t438 + rSges(6,3) * t477;
t388 = rSges(6,1) * t437 + rSges(6,2) * t436 + rSges(6,3) * t475;
t387 = Icges(6,1) * t439 + Icges(6,4) * t438 + Icges(6,5) * t477;
t386 = Icges(6,1) * t437 + Icges(6,4) * t436 + Icges(6,5) * t475;
t385 = Icges(6,4) * t439 + Icges(6,2) * t438 + Icges(6,6) * t477;
t384 = Icges(6,4) * t437 + Icges(6,2) * t436 + Icges(6,6) * t475;
t383 = Icges(6,5) * t439 + Icges(6,6) * t438 + Icges(6,3) * t477;
t382 = Icges(6,5) * t437 + Icges(6,6) * t436 + Icges(6,3) * t475;
t381 = rSges(7,1) * t434 + rSges(7,2) * t433 + rSges(7,3) * t477;
t380 = rSges(7,1) * t432 + rSges(7,2) * t431 + rSges(7,3) * t475;
t379 = Icges(7,1) * t434 + Icges(7,4) * t433 + Icges(7,5) * t477;
t378 = Icges(7,1) * t432 + Icges(7,4) * t431 + Icges(7,5) * t475;
t377 = Icges(7,4) * t434 + Icges(7,2) * t433 + Icges(7,6) * t477;
t376 = Icges(7,4) * t432 + Icges(7,2) * t431 + Icges(7,6) * t475;
t375 = Icges(7,5) * t434 + Icges(7,6) * t433 + Icges(7,3) * t477;
t374 = Icges(7,5) * t432 + Icges(7,6) * t431 + Icges(7,3) * t475;
t373 = t421 * t496 - t450 * t478 + t527;
t372 = -t420 * t496 + t450 * t479 + t524;
t371 = t420 * t478 - t421 * t479 + t537;
t370 = t419 * t496 + (-t449 - t462) * t478 + t526;
t369 = t449 * t479 + (-t418 - t427) * t496 + t523;
t368 = t418 * t478 + (-t419 - t428) * t479 + t529;
t367 = t389 * t465 - t405 * t429 + t522;
t366 = -t388 * t465 + t405 * t430 + t521;
t365 = t388 * t429 - t389 * t430 + t525;
t364 = t381 * t440 + t391 * t465 - t399 * t401 - t423 * t429 + t522;
t363 = -t380 * t440 - t390 * t465 + t400 * t401 + t423 * t430 + t521;
t362 = t380 * t399 - t381 * t400 + t390 * t429 - t391 * t430 + t525;
t1 = m(7) * (t362 ^ 2 + t363 ^ 2 + t364 ^ 2) / 0.2e1 + m(6) * (t365 ^ 2 + t366 ^ 2 + t367 ^ 2) / 0.2e1 + m(5) * (t368 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + m(4) * (t371 ^ 2 + t372 ^ 2 + t373 ^ 2) / 0.2e1 + m(3) * (t392 ^ 2 + t393 ^ 2 + t394 ^ 2) / 0.2e1 + t505 * ((t544 * t452 + (t454 * t517 + t456 * t514) * t512) * t504 - (t544 * t451 + (t453 * t517 + t455 * t514) * t512) * t532 + (t544 * t481 + (t482 * t517 + t483 * t514) * t512) * t505) / 0.2e1 + t400 * ((t375 * t475 + t377 * t431 + t379 * t432) * t399 + (t374 * t475 + t376 * t431 + t378 * t432) * t400 + (t396 * t475 + t397 * t431 + t398 * t432) * t440) / 0.2e1 + t429 * ((t383 * t477 + t385 * t438 + t387 * t439) * t429 + (t382 * t477 + t384 * t438 + t386 * t439) * t430 + (t402 * t477 + t403 * t438 + t404 * t439) * t465) / 0.2e1 + t440 * ((t375 * t490 + t377 * t466 + t379 * t467) * t399 + (t374 * t490 + t376 * t466 + t378 * t467) * t400 + (t396 * t490 + t466 * t397 + t467 * t398) * t440) / 0.2e1 + t399 * ((t375 * t477 + t377 * t433 + t379 * t434) * t399 + (t374 * t477 + t376 * t433 + t378 * t434) * t400 + (t396 * t477 + t397 * t433 + t398 * t434) * t440) / 0.2e1 + t465 * ((t383 * t490 + t385 * t472 + t387 * t473) * t429 + (t382 * t490 + t384 * t472 + t386 * t473) * t430 + (t402 * t490 + t403 * t472 + t404 * t473) * t465) / 0.2e1 + t430 * ((t383 * t475 + t385 * t436 + t387 * t437) * t429 + (t382 * t475 + t384 * t436 + t386 * t437) * t430 + (t402 * t475 + t403 * t436 + t404 * t437) * t465) / 0.2e1 - ((-t481 * t538 - t482 * t491 + t483 * t492) * t505 + ((-t454 * t491 + t456 * t492) * t515 + (t491 * t453 - t492 * t455 + t562) * t518) * t535) * t532 / 0.2e1 + ((t481 * t540 - t482 * t493 + t483 * t494) * t505 + (-(-t453 * t493 + t455 * t494) * t518 + (-t493 * t454 + t494 * t456 - t562) * t515) * t535) * t504 / 0.2e1 + ((t476 * t555 + t477 * t554 + t493 * t553) * t496 + (t476 * t561 + t477 * t557 + t493 * t559) * t479 + (t560 * t476 + t556 * t477 + t558 * t493) * t478) * t478 / 0.2e1 + ((t474 * t555 + t475 * t554 + t491 * t553) * t496 + (t561 * t474 + t557 * t475 + t559 * t491) * t479 + (t474 * t560 + t475 * t556 + t491 * t558) * t478) * t479 / 0.2e1 + ((t555 * t489 + t554 * t490 - t553 * t539) * t496 + (t489 * t561 + t490 * t557 - t539 * t559) * t479 + (t489 * t560 + t490 * t556 - t539 * t558) * t478) * t496 / 0.2e1 + (m(2) * (t499 ^ 2 + t500 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
