% Calculate kinetic energy for
% S6RRRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPP6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:11:07
% EndTime: 2019-03-09 21:11:10
% DurationCPUTime: 2.57s
% Computational Cost: add. (1628->261), mult. (2573->405), div. (0->0), fcn. (2612->8), ass. (0->135)
t512 = Icges(5,1) + Icges(6,2) + Icges(7,3);
t511 = -Icges(5,4) - Icges(6,6) + Icges(7,6);
t510 = -Icges(5,5) - Icges(7,5) + Icges(6,4);
t509 = Icges(5,2) + Icges(7,2) + Icges(6,3);
t508 = Icges(5,6) - Icges(6,5) - Icges(7,4);
t507 = -Icges(5,3) - Icges(7,1) - Icges(6,1);
t506 = rSges(7,1) + pkin(5);
t505 = rSges(7,3) + qJ(6);
t447 = qJ(3) + qJ(4);
t445 = sin(t447);
t446 = cos(t447);
t453 = cos(qJ(1));
t450 = sin(qJ(1));
t452 = cos(qJ(2));
t482 = t452 * t450;
t410 = t445 * t482 + t446 * t453;
t480 = t453 * t445;
t411 = t446 * t482 - t480;
t449 = sin(qJ(2));
t484 = t449 * t450;
t504 = t511 * t410 + t512 * t411 - t510 * t484;
t412 = -t450 * t446 + t452 * t480;
t481 = t452 * t453;
t413 = t445 * t450 + t446 * t481;
t483 = t449 * t453;
t503 = t511 * t412 + t512 * t413 - t510 * t483;
t502 = t509 * t410 + t511 * t411 - t508 * t484;
t501 = t509 * t412 + t511 * t413 - t508 * t483;
t500 = -t508 * t410 - t510 * t411 - t507 * t484;
t499 = -t508 * t412 - t510 * t413 - t507 * t483;
t498 = t507 * t452 + (-t508 * t445 - t510 * t446) * t449;
t497 = t508 * t452 + (t509 * t445 + t511 * t446) * t449;
t496 = t510 * t452 + (t511 * t445 + t512 * t446) * t449;
t451 = cos(qJ(3));
t491 = pkin(3) * t451;
t489 = Icges(3,4) * t449;
t488 = Icges(3,4) * t452;
t487 = t446 * t449;
t448 = sin(qJ(3));
t486 = t448 * t450;
t485 = t448 * t453;
t479 = rSges(7,2) * t410 + t505 * t411 + t506 * t484;
t478 = rSges(7,2) * t412 + t505 * t413 + t506 * t483;
t477 = (rSges(7,2) * t445 + rSges(7,3) * t446) * t449 + qJ(6) * t487 - t506 * t452;
t472 = pkin(2) * t452 + pkin(8) * t449;
t422 = t472 * t450;
t423 = t472 * t453;
t444 = qJD(2) * t450;
t475 = qJD(2) * t453;
t476 = t422 * t444 + t423 * t475;
t474 = qJD(3) * t449;
t425 = t453 * t474 + t444;
t473 = qJD(4) * t449;
t426 = t450 * t474 - t475;
t471 = rSges(3,1) * t452 - rSges(3,2) * t449;
t470 = Icges(3,1) * t452 - t489;
t469 = -Icges(3,2) * t449 + t488;
t468 = Icges(3,5) * t452 - Icges(3,6) * t449;
t402 = -Icges(3,6) * t453 + t450 * t469;
t405 = -Icges(3,5) * t453 + t450 * t470;
t467 = t402 * t449 - t405 * t452;
t403 = Icges(3,6) * t450 + t453 * t469;
t406 = Icges(3,5) * t450 + t453 * t470;
t466 = -t403 * t449 + t406 * t452;
t430 = Icges(3,2) * t452 + t489;
t431 = Icges(3,1) * t449 + t488;
t465 = -t430 * t449 + t431 * t452;
t461 = pkin(9) * t449 + t452 * t491;
t375 = -pkin(3) * t485 + t450 * t461;
t376 = pkin(3) * t486 + t453 * t461;
t464 = t425 * t375 - t376 * t426 + t476;
t428 = qJD(1) * (pkin(1) * t453 + pkin(7) * t450);
t436 = pkin(2) * t449 - pkin(8) * t452;
t463 = qJD(1) * t423 - t436 * t444 + t428;
t365 = pkin(4) * t411 + qJ(5) * t410;
t396 = t453 * t473 + t425;
t462 = qJD(5) * t449 * t445 + t396 * t365 + t464;
t437 = pkin(1) * t450 - pkin(7) * t453;
t460 = (-t422 - t437) * qJD(1) - t436 * t475;
t381 = -pkin(9) * t452 + t449 * t491;
t442 = -qJD(3) * t452 + qJD(1);
t459 = t442 * t376 - t381 * t425 + t463;
t366 = pkin(4) * t413 + qJ(5) * t412;
t427 = qJD(1) + (-qJD(3) - qJD(4)) * t452;
t458 = qJD(5) * t410 + t427 * t366 + t459;
t457 = -t375 * t442 + t426 * t381 + t460;
t397 = t450 * t473 + t426;
t416 = (pkin(4) * t446 + qJ(5) * t445) * t449;
t456 = qJD(5) * t412 + t397 * t416 + t457;
t434 = rSges(2,1) * t453 - rSges(2,2) * t450;
t433 = rSges(2,1) * t450 + rSges(2,2) * t453;
t432 = rSges(3,1) * t449 + rSges(3,2) * t452;
t429 = Icges(3,5) * t449 + Icges(3,6) * t452;
t421 = t451 * t481 + t486;
t420 = -t448 * t481 + t450 * t451;
t419 = t451 * t482 - t485;
t418 = -t448 * t482 - t451 * t453;
t409 = rSges(3,3) * t450 + t453 * t471;
t408 = -rSges(3,3) * t453 + t450 * t471;
t407 = -rSges(4,3) * t452 + (rSges(4,1) * t451 - rSges(4,2) * t448) * t449;
t404 = -Icges(4,5) * t452 + (Icges(4,1) * t451 - Icges(4,4) * t448) * t449;
t401 = -Icges(4,6) * t452 + (Icges(4,4) * t451 - Icges(4,2) * t448) * t449;
t400 = Icges(3,3) * t450 + t453 * t468;
t399 = -Icges(3,3) * t453 + t450 * t468;
t398 = -Icges(4,3) * t452 + (Icges(4,5) * t451 - Icges(4,6) * t448) * t449;
t393 = -rSges(6,1) * t452 + (-rSges(6,2) * t446 + rSges(6,3) * t445) * t449;
t391 = -rSges(5,3) * t452 + (rSges(5,1) * t446 - rSges(5,2) * t445) * t449;
t374 = rSges(4,1) * t421 + rSges(4,2) * t420 + rSges(4,3) * t483;
t373 = rSges(4,1) * t419 + rSges(4,2) * t418 + rSges(4,3) * t484;
t372 = Icges(4,1) * t421 + Icges(4,4) * t420 + Icges(4,5) * t483;
t371 = Icges(4,1) * t419 + Icges(4,4) * t418 + Icges(4,5) * t484;
t370 = Icges(4,4) * t421 + Icges(4,2) * t420 + Icges(4,6) * t483;
t369 = Icges(4,4) * t419 + Icges(4,2) * t418 + Icges(4,6) * t484;
t368 = Icges(4,5) * t421 + Icges(4,6) * t420 + Icges(4,3) * t483;
t367 = Icges(4,5) * t419 + Icges(4,6) * t418 + Icges(4,3) * t484;
t364 = qJD(1) * t409 - t432 * t444 + t428;
t363 = -t432 * t475 + (-t408 - t437) * qJD(1);
t362 = (t408 * t450 + t409 * t453) * qJD(2);
t360 = rSges(5,1) * t413 - rSges(5,2) * t412 + rSges(5,3) * t483;
t359 = rSges(5,1) * t411 - rSges(5,2) * t410 + rSges(5,3) * t484;
t358 = rSges(6,1) * t483 - rSges(6,2) * t413 + rSges(6,3) * t412;
t356 = rSges(6,1) * t484 - rSges(6,2) * t411 + rSges(6,3) * t410;
t333 = t374 * t442 - t407 * t425 + t463;
t332 = -t373 * t442 + t407 * t426 + t460;
t331 = t373 * t425 - t374 * t426 + t476;
t330 = t360 * t427 - t391 * t396 + t459;
t329 = -t359 * t427 + t391 * t397 + t457;
t328 = t359 * t396 - t360 * t397 + t464;
t327 = t358 * t427 + (-t393 - t416) * t396 + t458;
t326 = t393 * t397 + (-t356 - t365) * t427 + t456;
t325 = t356 * t396 + (-t358 - t366) * t397 + t462;
t324 = qJD(6) * t411 + t478 * t427 + (-t416 - t477) * t396 + t458;
t323 = qJD(6) * t413 + t477 * t397 + (-t365 - t479) * t427 + t456;
t322 = qJD(6) * t487 + t479 * t396 + (-t366 - t478) * t397 + t462;
t1 = ((t450 * t429 + t453 * t465) * qJD(1) + (t450 ^ 2 * t400 + (t467 * t453 + (-t399 + t466) * t450) * t453) * qJD(2)) * t444 / 0.2e1 - ((-t453 * t429 + t450 * t465) * qJD(1) + (t453 ^ 2 * t399 + (t466 * t450 + (-t400 + t467) * t453) * t450) * qJD(2)) * t475 / 0.2e1 + m(3) * (t362 ^ 2 + t363 ^ 2 + t364 ^ 2) / 0.2e1 + m(4) * (t331 ^ 2 + t332 ^ 2 + t333 ^ 2) / 0.2e1 + m(5) * (t328 ^ 2 + t329 ^ 2 + t330 ^ 2) / 0.2e1 + m(6) * (t325 ^ 2 + t326 ^ 2 + t327 ^ 2) / 0.2e1 + m(7) * (t322 ^ 2 + t323 ^ 2 + t324 ^ 2) / 0.2e1 + qJD(1) * ((t430 * t452 + t431 * t449) * qJD(1) + ((t403 * t452 + t406 * t449) * t450 - (t402 * t452 + t405 * t449) * t453) * qJD(2)) / 0.2e1 + t425 * ((t368 * t483 + t420 * t370 + t421 * t372) * t425 + (t367 * t483 + t369 * t420 + t371 * t421) * t426 + (t398 * t483 + t401 * t420 + t404 * t421) * t442) / 0.2e1 + t426 * ((t368 * t484 + t370 * t418 + t372 * t419) * t425 + (t367 * t484 + t418 * t369 + t419 * t371) * t426 + (t398 * t484 + t401 * t418 + t404 * t419) * t442) / 0.2e1 + t442 * ((-t367 * t426 - t368 * t425 - t398 * t442) * t452 + ((-t370 * t448 + t372 * t451) * t425 + (-t369 * t448 + t371 * t451) * t426 + (-t401 * t448 + t404 * t451) * t442) * t449) / 0.2e1 + (Icges(2,3) + m(2) * (t433 ^ 2 + t434 ^ 2)) * qJD(1) ^ 2 / 0.2e1 + ((t497 * t412 + t496 * t413 + t498 * t483) * t427 + (t502 * t412 + t504 * t413 + t500 * t483) * t397 + (t501 * t412 + t503 * t413 + t499 * t483) * t396) * t396 / 0.2e1 + ((t497 * t410 + t496 * t411 + t498 * t484) * t427 + (t502 * t410 + t504 * t411 + t500 * t484) * t397 + (t501 * t410 + t503 * t411 + t499 * t484) * t396) * t397 / 0.2e1 + ((-t499 * t396 - t500 * t397 - t498 * t427) * t452 + ((t497 * t445 + t496 * t446) * t427 + (t502 * t445 + t504 * t446) * t397 + (t501 * t445 + t503 * t446) * t396) * t449) * t427 / 0.2e1;
T  = t1;
