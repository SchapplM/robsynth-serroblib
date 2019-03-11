% Calculate kinetic energy for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
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
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP6_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP6_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPP6_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPP6_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:46:20
% EndTime: 2019-03-09 04:46:22
% DurationCPUTime: 1.91s
% Computational Cost: add. (1116->238), mult. (1850->353), div. (0->0), fcn. (1850->8), ass. (0->123)
t493 = Icges(6,1) + Icges(7,1);
t492 = Icges(6,4) - Icges(7,5);
t491 = Icges(7,4) + Icges(6,5);
t490 = Icges(6,2) + Icges(7,3);
t489 = Icges(7,6) - Icges(6,6);
t488 = Icges(6,3) + Icges(7,2) + Icges(5,3);
t487 = rSges(7,1) + pkin(5);
t486 = rSges(7,3) + qJ(6);
t427 = qJ(4) + pkin(9);
t423 = sin(t427);
t424 = cos(t427);
t434 = cos(qJ(1));
t459 = t434 * t424;
t430 = sin(qJ(3));
t431 = sin(qJ(1));
t465 = t430 * t431;
t381 = t423 * t465 - t459;
t382 = t423 * t434 + t424 * t465;
t433 = cos(qJ(3));
t462 = t431 * t433;
t485 = t490 * t381 - t492 * t382 - t489 * t462;
t464 = t430 * t434;
t383 = t423 * t464 + t424 * t431;
t384 = t423 * t431 - t430 * t459;
t460 = t433 * t434;
t484 = -t490 * t383 - t492 * t384 + t489 * t460;
t483 = -t492 * t381 + t493 * t382 - t491 * t462;
t482 = t492 * t383 + t493 * t384 + t491 * t460;
t481 = (t490 * t423 - t492 * t424) * t433 + t489 * t430;
t480 = (-t492 * t423 + t493 * t424) * t433 + t491 * t430;
t429 = sin(qJ(4));
t432 = cos(qJ(4));
t461 = t432 * t434;
t400 = -t429 * t465 + t461;
t463 = t431 * t432;
t466 = t429 * t434;
t401 = t430 * t463 + t466;
t479 = Icges(5,5) * t401 + Icges(5,6) * t400 + t489 * t381 + t491 * t382 - t488 * t462;
t402 = t429 * t464 + t463;
t467 = t429 * t431;
t403 = -t430 * t461 + t467;
t478 = Icges(5,5) * t403 + Icges(5,6) * t402 - t489 * t383 + t491 * t384 + t488 * t460;
t477 = (Icges(5,5) * t432 - Icges(5,6) * t429 + t489 * t423 + t491 * t424) * t433 + t488 * t430;
t471 = pkin(4) * t432;
t476 = -qJ(5) * t433 + t430 * t471;
t469 = Icges(4,4) * t430;
t468 = Icges(4,4) * t433;
t458 = -rSges(7,2) * t462 + t486 * t381 + t487 * t382;
t457 = rSges(7,2) * t460 - t486 * t383 + t487 * t384;
t456 = rSges(7,2) * t430 + (t486 * t423 + t487 * t424) * t433;
t409 = qJD(1) * (pkin(1) * t434 + qJ(2) * t431);
t455 = qJD(1) * t434 * pkin(7) + t409;
t454 = qJD(3) * t431;
t453 = qJD(3) * t434;
t452 = qJD(4) * t433;
t451 = qJD(5) * t433;
t413 = pkin(1) * t431 - qJ(2) * t434;
t450 = -pkin(7) * t431 - t413;
t449 = pkin(3) * t430 - pkin(8) * t433;
t404 = t449 * t431;
t405 = t449 * t434;
t448 = -t404 * t454 - t405 * t453;
t447 = rSges(4,1) * t430 + rSges(4,2) * t433;
t446 = Icges(4,1) * t430 + t468;
t445 = Icges(4,2) * t433 + t469;
t444 = Icges(4,5) * t430 + Icges(4,6) * t433;
t389 = Icges(4,6) * t434 + t431 * t445;
t392 = Icges(4,5) * t434 + t431 * t446;
t443 = -t389 * t433 - t392 * t430;
t390 = Icges(4,6) * t431 - t434 * t445;
t393 = Icges(4,5) * t431 - t434 * t446;
t442 = t390 * t433 + t393 * t430;
t411 = -Icges(4,2) * t430 + t468;
t412 = Icges(4,1) * t433 - t469;
t441 = t411 * t433 + t412 * t430;
t366 = pkin(4) * t467 - t434 * t476;
t408 = -t431 * t452 + t453;
t440 = qJD(5) * t430 + t408 * t366 + t448;
t417 = pkin(3) * t433 + pkin(8) * t430;
t426 = qJD(2) * t431;
t439 = t417 * t454 + t426 + (t405 + t450) * qJD(1);
t438 = qJD(1) * t404 + (-qJD(3) * t417 - qJD(2)) * t434 + t455;
t365 = pkin(4) * t466 + t431 * t476;
t420 = qJD(4) * t430 + qJD(1);
t437 = t420 * t365 + t434 * t451 + t438;
t372 = qJ(5) * t430 + t433 * t471;
t407 = t434 * t452 + t454;
t436 = t407 * t372 - t431 * t451 + t439;
t416 = rSges(2,1) * t434 - rSges(2,2) * t431;
t415 = rSges(4,1) * t433 - rSges(4,2) * t430;
t414 = rSges(2,1) * t431 + rSges(2,2) * t434;
t410 = Icges(4,5) * t433 - Icges(4,6) * t430;
t397 = rSges(4,3) * t431 - t434 * t447;
t396 = rSges(5,3) * t430 + (rSges(5,1) * t432 - rSges(5,2) * t429) * t433;
t395 = rSges(4,3) * t434 + t431 * t447;
t391 = Icges(5,5) * t430 + (Icges(5,1) * t432 - Icges(5,4) * t429) * t433;
t388 = Icges(5,6) * t430 + (Icges(5,4) * t432 - Icges(5,2) * t429) * t433;
t387 = Icges(4,3) * t431 - t434 * t444;
t386 = Icges(4,3) * t434 + t431 * t444;
t380 = rSges(6,3) * t430 + (rSges(6,1) * t424 - rSges(6,2) * t423) * t433;
t371 = t409 - qJD(2) * t434 + qJD(1) * (-rSges(3,2) * t434 + rSges(3,3) * t431);
t370 = t426 + (rSges(3,2) * t431 + rSges(3,3) * t434 - t413) * qJD(1);
t368 = rSges(5,1) * t403 + rSges(5,2) * t402 + rSges(5,3) * t460;
t367 = rSges(5,1) * t401 + rSges(5,2) * t400 - rSges(5,3) * t462;
t364 = Icges(5,1) * t403 + Icges(5,4) * t402 + Icges(5,5) * t460;
t363 = Icges(5,1) * t401 + Icges(5,4) * t400 - Icges(5,5) * t462;
t362 = Icges(5,4) * t403 + Icges(5,2) * t402 + Icges(5,6) * t460;
t361 = Icges(5,4) * t401 + Icges(5,2) * t400 - Icges(5,6) * t462;
t356 = (-t395 * t431 + t397 * t434) * qJD(3);
t354 = rSges(6,1) * t384 + rSges(6,2) * t383 + rSges(6,3) * t460;
t352 = rSges(6,1) * t382 - rSges(6,2) * t381 - rSges(6,3) * t462;
t337 = qJD(1) * t395 + (-qJD(3) * t415 - qJD(2)) * t434 + t455;
t336 = t415 * t454 + t426 + (-t397 + t450) * qJD(1);
t335 = t367 * t420 - t396 * t408 + t438;
t334 = -t368 * t420 + t396 * t407 + t439;
t333 = -t367 * t407 + t368 * t408 + t448;
t332 = t352 * t420 + (-t372 - t380) * t408 + t437;
t331 = t380 * t407 + (-t354 - t366) * t420 + t436;
t330 = t354 * t408 + (-t352 - t365) * t407 + t440;
t329 = -qJD(6) * t383 + t458 * t420 + (-t372 - t456) * t408 + t437;
t328 = qJD(6) * t381 + t456 * t407 + (-t366 - t457) * t420 + t436;
t327 = qJD(6) * t423 * t433 + t457 * t408 + (-t365 - t458) * t407 + t440;
t1 = ((t434 * t410 + t431 * t441) * qJD(1) + (t434 ^ 2 * t386 + (t442 * t431 + (t387 - t443) * t434) * t431) * qJD(3)) * t453 / 0.2e1 + ((t431 * t410 - t434 * t441) * qJD(1) + (t431 ^ 2 * t387 + (t443 * t434 + (t386 - t442) * t431) * t434) * qJD(3)) * t454 / 0.2e1 + qJD(1) * ((-t430 * t411 + t433 * t412) * qJD(1) + ((-t389 * t430 + t392 * t433) * t434 + (-t390 * t430 + t393 * t433) * t431) * qJD(3)) / 0.2e1 + m(3) * (t370 ^ 2 + t371 ^ 2) / 0.2e1 + m(4) * (t336 ^ 2 + t337 ^ 2 + t356 ^ 2) / 0.2e1 + m(5) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + m(6) * (t330 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + m(7) * (t327 ^ 2 + t328 ^ 2 + t329 ^ 2) / 0.2e1 + ((-t383 * t481 + t384 * t480 + t388 * t402 + t391 * t403 + t460 * t477) * t420 + (t361 * t402 + t363 * t403 - t383 * t485 + t483 * t384 + t479 * t460) * t408 + (t402 * t362 + t403 * t364 - t383 * t484 + t384 * t482 + t460 * t478) * t407) * t407 / 0.2e1 + ((t381 * t481 + t382 * t480 + t388 * t400 + t391 * t401 - t462 * t477) * t420 + (t400 * t361 + t401 * t363 + t381 * t485 + t483 * t382 - t479 * t462) * t408 + (t362 * t400 + t364 * t401 + t381 * t484 + t382 * t482 - t462 * t478) * t407) * t408 / 0.2e1 + (((-t388 * t429 + t391 * t432 + t423 * t481 + t424 * t480) * t420 + (-t361 * t429 + t363 * t432 + t423 * t485 + t483 * t424) * t408 + (-t362 * t429 + t364 * t432 + t423 * t484 + t424 * t482) * t407) * t433 + (t407 * t478 + t408 * t479 + t420 * t477) * t430) * t420 / 0.2e1 + (Icges(3,1) + m(2) * (t414 ^ 2 + t416 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
