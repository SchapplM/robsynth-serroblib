% Calculate kinetic energy for
% S6RRRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
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
% Datum: 2019-03-09 16:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:03:34
% EndTime: 2019-03-09 16:03:37
% DurationCPUTime: 2.50s
% Computational Cost: add. (2001->290), mult. (4999->430), div. (0->0), fcn. (5959->10), ass. (0->136)
t538 = Icges(4,1) + Icges(5,1) + Icges(6,2);
t537 = Icges(6,1) + Icges(4,2) + Icges(5,3);
t536 = -Icges(4,4) - Icges(6,4) + Icges(5,5);
t535 = Icges(5,4) + Icges(4,5) + Icges(6,6);
t534 = Icges(6,5) + Icges(4,6) - Icges(5,6);
t533 = Icges(5,2) + Icges(4,3) + Icges(6,3);
t487 = sin(qJ(2));
t488 = sin(qJ(1));
t490 = cos(qJ(2));
t491 = cos(qJ(1));
t516 = cos(pkin(6));
t502 = t491 * t516;
t467 = t487 * t488 - t490 * t502;
t468 = t487 * t502 + t488 * t490;
t485 = sin(pkin(6));
t513 = t485 * t491;
t429 = Icges(3,5) * t468 - Icges(3,6) * t467 - Icges(3,3) * t513;
t503 = t488 * t516;
t469 = t491 * t487 + t490 * t503;
t470 = -t487 * t503 + t491 * t490;
t515 = t485 * t488;
t430 = Icges(3,5) * t470 - Icges(3,6) * t469 + Icges(3,3) * t515;
t532 = (t429 * t491 - t430 * t488) * t485;
t518 = cos(qJ(3));
t506 = t485 * t518;
t517 = sin(qJ(3));
t449 = t468 * t517 + t491 * t506;
t505 = t485 * t517;
t450 = t468 * t518 - t491 * t505;
t531 = t449 * t534 - t450 * t535 - t467 * t533;
t451 = t470 * t517 - t488 * t506;
t452 = t470 * t518 + t488 * t505;
t530 = t451 * t534 - t452 * t535 - t469 * t533;
t529 = t449 * t537 + t450 * t536 - t467 * t534;
t528 = t451 * t537 + t452 * t536 - t469 * t534;
t527 = t536 * t449 + t450 * t538 + t535 * t467;
t526 = t536 * t451 + t452 * t538 + t535 * t469;
t465 = t487 * t505 - t516 * t518;
t466 = t487 * t506 + t516 * t517;
t514 = t485 * t490;
t525 = t465 * t534 - t466 * t535 + t514 * t533;
t524 = t465 * t537 + t466 * t536 + t514 * t534;
t523 = t536 * t465 + t466 * t538 - t535 * t514;
t404 = pkin(3) * t450 + qJ(4) * t449;
t415 = pkin(4) * t450 - qJ(5) * t467;
t512 = -t404 - t415;
t405 = pkin(3) * t452 + qJ(4) * t451;
t416 = pkin(4) * t452 - qJ(5) * t469;
t511 = -t405 - t416;
t442 = pkin(2) * t468 + pkin(9) * t467;
t443 = pkin(2) * t470 + pkin(9) * t469;
t507 = qJD(2) * t485;
t481 = t488 * t507;
t504 = t491 * t507;
t510 = t442 * t481 + t443 * t504;
t440 = pkin(3) * t466 + qJ(4) * t465;
t455 = pkin(4) * t466 + qJ(5) * t514;
t509 = -t440 - t455;
t453 = qJD(3) * t469 + t481;
t508 = qJD(1) * (pkin(1) * t488 - pkin(8) * t513);
t482 = qJD(2) * t516 + qJD(1);
t501 = qJD(4) * t465 + t453 * t404 + t510;
t454 = qJD(3) * t467 - t504;
t472 = -qJD(3) * t514 + t482;
t499 = qJD(5) * t514 + t453 * t415 + t501;
t471 = (pkin(2) * t487 - pkin(9) * t490) * t485;
t473 = qJD(1) * (pkin(1) * t491 + pkin(8) * t515);
t498 = t482 * t443 - t471 * t481 + t473;
t497 = qJD(4) * t449 + t472 * t405 + t498;
t496 = -t442 * t482 - t471 * t504 - t508;
t495 = -qJD(5) * t467 + t472 * t416 + t497;
t494 = qJD(4) * t451 + t454 * t440 + t496;
t493 = -qJD(5) * t469 + t454 * t455 + t494;
t489 = cos(qJ(6));
t486 = sin(qJ(6));
t476 = rSges(2,1) * t491 - rSges(2,2) * t488;
t475 = rSges(2,1) * t488 + rSges(2,2) * t491;
t459 = t516 * rSges(3,3) + (rSges(3,1) * t487 + rSges(3,2) * t490) * t485;
t458 = Icges(3,5) * t516 + (Icges(3,1) * t487 + Icges(3,4) * t490) * t485;
t457 = Icges(3,6) * t516 + (Icges(3,4) * t487 + Icges(3,2) * t490) * t485;
t456 = Icges(3,3) * t516 + (Icges(3,5) * t487 + Icges(3,6) * t490) * t485;
t448 = t465 * t489 + t486 * t514;
t447 = -t465 * t486 + t489 * t514;
t444 = qJD(6) * t466 + t472;
t441 = pkin(5) * t465 + pkin(10) * t466;
t437 = rSges(3,1) * t470 - rSges(3,2) * t469 + rSges(3,3) * t515;
t436 = rSges(3,1) * t468 - rSges(3,2) * t467 - rSges(3,3) * t513;
t434 = Icges(3,1) * t470 - Icges(3,4) * t469 + Icges(3,5) * t515;
t433 = Icges(3,1) * t468 - Icges(3,4) * t467 - Icges(3,5) * t513;
t432 = Icges(3,4) * t470 - Icges(3,2) * t469 + Icges(3,6) * t515;
t431 = Icges(3,4) * t468 - Icges(3,2) * t467 - Icges(3,6) * t513;
t428 = rSges(4,1) * t466 - rSges(4,2) * t465 - rSges(4,3) * t514;
t427 = rSges(5,1) * t466 - rSges(5,2) * t514 + rSges(5,3) * t465;
t426 = rSges(6,1) * t465 - rSges(6,2) * t466 + rSges(6,3) * t514;
t414 = t451 * t489 - t469 * t486;
t413 = -t451 * t486 - t469 * t489;
t412 = t449 * t489 - t467 * t486;
t411 = -t449 * t486 - t467 * t489;
t409 = qJD(6) * t450 + t454;
t408 = qJD(6) * t452 + t453;
t407 = pkin(5) * t451 + pkin(10) * t452;
t406 = pkin(5) * t449 + pkin(10) * t450;
t399 = rSges(4,1) * t452 - rSges(4,2) * t451 + rSges(4,3) * t469;
t398 = rSges(5,1) * t452 + rSges(5,2) * t469 + rSges(5,3) * t451;
t397 = rSges(6,1) * t451 - rSges(6,2) * t452 - rSges(6,3) * t469;
t396 = rSges(4,1) * t450 - rSges(4,2) * t449 + rSges(4,3) * t467;
t395 = rSges(5,1) * t450 + rSges(5,2) * t467 + rSges(5,3) * t449;
t394 = rSges(6,1) * t449 - rSges(6,2) * t450 - rSges(6,3) * t467;
t375 = rSges(7,1) * t448 + rSges(7,2) * t447 + rSges(7,3) * t466;
t374 = Icges(7,1) * t448 + Icges(7,4) * t447 + Icges(7,5) * t466;
t373 = Icges(7,4) * t448 + Icges(7,2) * t447 + Icges(7,6) * t466;
t372 = Icges(7,5) * t448 + Icges(7,6) * t447 + Icges(7,3) * t466;
t370 = t437 * t482 - t459 * t481 + t473;
t369 = -t436 * t482 - t459 * t504 - t508;
t368 = (t436 * t488 + t437 * t491) * t507;
t367 = rSges(7,1) * t414 + rSges(7,2) * t413 + rSges(7,3) * t452;
t366 = rSges(7,1) * t412 + rSges(7,2) * t411 + rSges(7,3) * t450;
t365 = Icges(7,1) * t414 + Icges(7,4) * t413 + Icges(7,5) * t452;
t364 = Icges(7,1) * t412 + Icges(7,4) * t411 + Icges(7,5) * t450;
t363 = Icges(7,4) * t414 + Icges(7,2) * t413 + Icges(7,6) * t452;
t362 = Icges(7,4) * t412 + Icges(7,2) * t411 + Icges(7,6) * t450;
t361 = Icges(7,5) * t414 + Icges(7,6) * t413 + Icges(7,3) * t452;
t360 = Icges(7,5) * t412 + Icges(7,6) * t411 + Icges(7,3) * t450;
t359 = t399 * t472 - t428 * t453 + t498;
t358 = -t396 * t472 + t428 * t454 + t496;
t357 = t396 * t453 - t399 * t454 + t510;
t356 = t398 * t472 + (-t427 - t440) * t453 + t497;
t355 = t427 * t454 + (-t395 - t404) * t472 + t494;
t354 = t395 * t453 + (-t398 - t405) * t454 + t501;
t353 = t397 * t472 + (-t426 + t509) * t453 + t495;
t352 = t426 * t454 + (-t394 + t512) * t472 + t493;
t351 = t394 * t453 + (-t397 + t511) * t454 + t499;
t350 = t367 * t444 - t375 * t408 + t407 * t472 + (-t441 + t509) * t453 + t495;
t349 = -t366 * t444 + t375 * t409 + t441 * t454 + (-t406 + t512) * t472 + t493;
t348 = t366 * t408 - t367 * t409 + t406 * t453 + (-t407 + t511) * t454 + t499;
t1 = ((t456 * t515 - t457 * t469 + t458 * t470) * t482 + (-(-t431 * t469 + t433 * t470) * t491 + (-t469 * t432 + t470 * t434 - t532) * t488) * t507) * t481 / 0.2e1 - ((-t456 * t513 - t457 * t467 + t458 * t468) * t482 + ((-t432 * t467 + t434 * t468) * t488 + (t431 * t467 - t433 * t468 + t532) * t491) * t507) * t504 / 0.2e1 + m(4) * (t357 ^ 2 + t358 ^ 2 + t359 ^ 2) / 0.2e1 + m(5) * (t354 ^ 2 + t355 ^ 2 + t356 ^ 2) / 0.2e1 + m(6) * (t351 ^ 2 + t352 ^ 2 + t353 ^ 2) / 0.2e1 + m(7) * (t348 ^ 2 + t349 ^ 2 + t350 ^ 2) / 0.2e1 + t482 * ((t516 * t430 + (t432 * t490 + t434 * t487) * t485) * t481 - (t516 * t429 + (t431 * t490 + t433 * t487) * t485) * t504 + (t516 * t456 + (t457 * t490 + t458 * t487) * t485) * t482) / 0.2e1 + t408 * ((t452 * t361 + t413 * t363 + t414 * t365) * t408 + (t360 * t452 + t362 * t413 + t364 * t414) * t409 + (t372 * t452 + t373 * t413 + t374 * t414) * t444) / 0.2e1 + t409 * ((t361 * t450 + t363 * t411 + t365 * t412) * t408 + (t450 * t360 + t411 * t362 + t412 * t364) * t409 + (t372 * t450 + t373 * t411 + t374 * t412) * t444) / 0.2e1 + t444 * ((t361 * t466 + t363 * t447 + t365 * t448) * t408 + (t360 * t466 + t362 * t447 + t364 * t448) * t409 + (t372 * t466 + t373 * t447 + t374 * t448) * t444) / 0.2e1 + m(3) * (t368 ^ 2 + t369 ^ 2 + t370 ^ 2) / 0.2e1 + (m(2) * (t475 ^ 2 + t476 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((t451 * t524 + t452 * t523 - t469 * t525) * t472 + (t451 * t529 + t452 * t527 - t469 * t531) * t454 + (t528 * t451 + t526 * t452 - t530 * t469) * t453) * t453 / 0.2e1 + ((t449 * t524 + t450 * t523 - t467 * t525) * t472 + (t529 * t449 + t527 * t450 - t531 * t467) * t454 + (t449 * t528 + t450 * t526 - t467 * t530) * t453) * t454 / 0.2e1 + ((t524 * t465 + t523 * t466 + t525 * t514) * t472 + (t465 * t529 + t466 * t527 + t514 * t531) * t454 + (t465 * t528 + t466 * t526 + t514 * t530) * t453) * t472 / 0.2e1;
T  = t1;
