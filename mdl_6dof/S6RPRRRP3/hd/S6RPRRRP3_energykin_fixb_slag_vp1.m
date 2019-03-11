% Calculate kinetic energy for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP3_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP3_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP3_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP3_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP3_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:02:52
% EndTime: 2019-03-09 06:02:54
% DurationCPUTime: 2.17s
% Computational Cost: add. (1895->243), mult. (1899->384), div. (0->0), fcn. (1898->10), ass. (0->130)
t496 = Icges(6,1) + Icges(7,1);
t495 = -Icges(6,4) + Icges(7,5);
t494 = Icges(7,4) + Icges(6,5);
t493 = Icges(6,2) + Icges(7,3);
t492 = -Icges(7,6) + Icges(6,6);
t491 = -Icges(6,3) - Icges(7,2);
t490 = rSges(7,1) + pkin(5);
t489 = rSges(7,3) + qJ(6);
t430 = qJ(1) + pkin(10);
t426 = sin(t430);
t427 = cos(t430);
t431 = qJ(4) + qJ(5);
t429 = cos(t431);
t428 = sin(t431);
t436 = cos(qJ(3));
t466 = t428 * t436;
t381 = t426 * t466 + t427 * t429;
t465 = t429 * t436;
t382 = t426 * t465 - t427 * t428;
t433 = sin(qJ(3));
t469 = t426 * t433;
t488 = t493 * t381 + t495 * t382 - t492 * t469;
t383 = -t426 * t429 + t427 * t466;
t384 = t426 * t428 + t427 * t465;
t467 = t427 * t433;
t487 = t493 * t383 + t495 * t384 - t492 * t467;
t486 = -t492 * t381 + t494 * t382 - t491 * t469;
t485 = -t492 * t383 + t494 * t384 - t491 * t467;
t484 = t495 * t381 + t496 * t382 + t494 * t469;
t483 = t495 * t383 + t496 * t384 + t494 * t467;
t482 = t492 * t436 + (t493 * t428 + t495 * t429) * t433;
t481 = t491 * t436 + (-t492 * t428 + t494 * t429) * t433;
t480 = -t494 * t436 + (t495 * t428 + t496 * t429) * t433;
t434 = sin(qJ(1));
t475 = pkin(1) * t434;
t435 = cos(qJ(4));
t474 = pkin(4) * t435;
t472 = Icges(4,4) * t433;
t471 = Icges(4,4) * t436;
t432 = sin(qJ(4));
t470 = t426 * t432;
t468 = t427 * t432;
t464 = t432 * t436;
t463 = t435 * t436;
t462 = rSges(7,2) * t469 + t489 * t381 + t382 * t490;
t461 = rSges(7,2) * t467 + t489 * t383 + t384 * t490;
t460 = -rSges(7,2) * t436 + (t489 * t428 + t429 * t490) * t433;
t437 = cos(qJ(1));
t425 = qJD(1) * t437 * pkin(1);
t459 = qJD(1) * (pkin(2) * t427 + pkin(7) * t426) + t425;
t422 = qJD(3) * t426;
t457 = qJD(4) * t433;
t407 = t427 * t457 + t422;
t458 = qJD(3) * t427;
t456 = qJD(5) * t433;
t453 = pkin(3) * t436 + pkin(8) * t433;
t404 = t453 * t426;
t405 = t453 * t427;
t455 = t404 * t422 + t405 * t458 + qJD(2);
t454 = -pkin(2) * t426 + pkin(7) * t427 - t475;
t408 = t426 * t457 - t458;
t452 = rSges(4,1) * t436 - rSges(4,2) * t433;
t451 = Icges(4,1) * t436 - t472;
t450 = -Icges(4,2) * t433 + t471;
t449 = Icges(4,5) * t436 - Icges(4,6) * t433;
t372 = -Icges(4,6) * t427 + t426 * t450;
t374 = -Icges(4,5) * t427 + t426 * t451;
t448 = t372 * t433 - t374 * t436;
t373 = Icges(4,6) * t426 + t427 * t450;
t375 = Icges(4,5) * t426 + t427 * t451;
t447 = -t373 * t433 + t375 * t436;
t414 = Icges(4,2) * t436 + t472;
t415 = Icges(4,1) * t433 + t471;
t446 = -t414 * t433 + t415 * t436;
t421 = pkin(3) * t433 - pkin(8) * t436;
t445 = qJD(1) * t405 - t421 * t422 + t459;
t443 = pkin(9) * t433 + t436 * t474;
t363 = -pkin(4) * t468 + t426 * t443;
t364 = pkin(4) * t470 + t427 * t443;
t444 = t407 * t363 - t364 * t408 + t455;
t442 = (-t404 + t454) * qJD(1) - t421 * t458;
t380 = -pkin(9) * t436 + t433 * t474;
t423 = -qJD(4) * t436 + qJD(1);
t441 = t423 * t364 - t380 * t407 + t445;
t440 = -t363 * t423 + t408 * t380 + t442;
t420 = rSges(2,1) * t437 - rSges(2,2) * t434;
t419 = rSges(2,1) * t434 + rSges(2,2) * t437;
t418 = rSges(4,1) * t433 + rSges(4,2) * t436;
t413 = Icges(4,5) * t433 + Icges(4,6) * t436;
t411 = qJD(1) + (-qJD(4) - qJD(5)) * t436;
t403 = -rSges(5,3) * t436 + (rSges(5,1) * t435 - rSges(5,2) * t432) * t433;
t402 = -Icges(5,5) * t436 + (Icges(5,1) * t435 - Icges(5,4) * t432) * t433;
t401 = -Icges(5,6) * t436 + (Icges(5,4) * t435 - Icges(5,2) * t432) * t433;
t400 = -Icges(5,3) * t436 + (Icges(5,5) * t435 - Icges(5,6) * t432) * t433;
t399 = t427 * t463 + t470;
t398 = t426 * t435 - t427 * t464;
t397 = t426 * t463 - t468;
t396 = -t426 * t464 - t427 * t435;
t394 = t425 + qJD(1) * (rSges(3,1) * t427 - rSges(3,2) * t426);
t393 = (-rSges(3,1) * t426 - rSges(3,2) * t427 - t475) * qJD(1);
t392 = -rSges(6,3) * t436 + (rSges(6,1) * t429 - rSges(6,2) * t428) * t433;
t377 = rSges(4,3) * t426 + t427 * t452;
t376 = -rSges(4,3) * t427 + t426 * t452;
t371 = Icges(4,3) * t426 + t427 * t449;
t370 = -Icges(4,3) * t427 + t426 * t449;
t369 = t426 * t456 + t408;
t368 = t427 * t456 + t407;
t362 = rSges(5,1) * t399 + rSges(5,2) * t398 + rSges(5,3) * t467;
t361 = rSges(5,1) * t397 + rSges(5,2) * t396 + rSges(5,3) * t469;
t360 = Icges(5,1) * t399 + Icges(5,4) * t398 + Icges(5,5) * t467;
t359 = Icges(5,1) * t397 + Icges(5,4) * t396 + Icges(5,5) * t469;
t358 = Icges(5,4) * t399 + Icges(5,2) * t398 + Icges(5,6) * t467;
t357 = Icges(5,4) * t397 + Icges(5,2) * t396 + Icges(5,6) * t469;
t356 = Icges(5,5) * t399 + Icges(5,6) * t398 + Icges(5,3) * t467;
t355 = Icges(5,5) * t397 + Icges(5,6) * t396 + Icges(5,3) * t469;
t353 = rSges(6,1) * t384 - rSges(6,2) * t383 + rSges(6,3) * t467;
t351 = rSges(6,1) * t382 - rSges(6,2) * t381 + rSges(6,3) * t469;
t337 = qJD(1) * t377 - t418 * t422 + t459;
t336 = -t418 * t458 + (-t376 + t454) * qJD(1);
t335 = qJD(2) + (t376 * t426 + t377 * t427) * qJD(3);
t333 = t362 * t423 - t403 * t407 + t445;
t332 = -t361 * t423 + t403 * t408 + t442;
t331 = t361 * t407 - t362 * t408 + t455;
t330 = t353 * t411 - t368 * t392 + t441;
t329 = -t351 * t411 + t369 * t392 + t440;
t328 = t351 * t368 - t353 * t369 + t444;
t327 = qJD(6) * t381 - t368 * t460 + t411 * t461 + t441;
t326 = qJD(6) * t383 + t369 * t460 - t411 * t462 + t440;
t325 = qJD(6) * t428 * t433 + t368 * t462 - t369 * t461 + t444;
t1 = m(7) * (t325 ^ 2 + t326 ^ 2 + t327 ^ 2) / 0.2e1 + m(5) * (t331 ^ 2 + t332 ^ 2 + t333 ^ 2) / 0.2e1 + m(4) * (t335 ^ 2 + t336 ^ 2 + t337 ^ 2) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t393 ^ 2 + t394 ^ 2) / 0.2e1 + t423 * ((-t355 * t408 - t356 * t407 - t400 * t423) * t436 + ((-t358 * t432 + t360 * t435) * t407 + (-t357 * t432 + t359 * t435) * t408 + (-t401 * t432 + t402 * t435) * t423) * t433) / 0.2e1 + t407 * ((t356 * t467 + t398 * t358 + t399 * t360) * t407 + (t355 * t467 + t357 * t398 + t359 * t399) * t408 + (t398 * t401 + t399 * t402 + t400 * t467) * t423) / 0.2e1 + t408 * ((t356 * t469 + t358 * t396 + t360 * t397) * t407 + (t355 * t469 + t396 * t357 + t397 * t359) * t408 + (t396 * t401 + t397 * t402 + t400 * t469) * t423) / 0.2e1 - ((-t427 * t413 + t426 * t446) * qJD(1) + (t427 ^ 2 * t370 + (t447 * t426 + (-t371 + t448) * t427) * t426) * qJD(3)) * t458 / 0.2e1 + qJD(1) * ((t436 * t414 + t433 * t415) * qJD(1) + ((t373 * t436 + t375 * t433) * t426 - (t372 * t436 + t374 * t433) * t427) * qJD(3)) / 0.2e1 + m(6) * (t328 ^ 2 + t329 ^ 2 + t330 ^ 2) / 0.2e1 + ((t426 * t413 + t446 * t427) * qJD(1) + (t426 ^ 2 * t371 + (t448 * t427 + (-t370 + t447) * t426) * t427) * qJD(3)) * t422 / 0.2e1 + ((t482 * t383 + t480 * t384 + t481 * t467) * t411 + (t488 * t383 + t484 * t384 + t486 * t467) * t369 + (t487 * t383 + t483 * t384 + t485 * t467) * t368) * t368 / 0.2e1 + ((t482 * t381 + t480 * t382 + t481 * t469) * t411 + (t488 * t381 + t484 * t382 + t486 * t469) * t369 + (t487 * t381 + t483 * t382 + t485 * t469) * t368) * t369 / 0.2e1 + ((-t485 * t368 - t486 * t369 - t481 * t411) * t436 + ((t482 * t428 + t480 * t429) * t411 + (t488 * t428 + t484 * t429) * t369 + (t487 * t428 + t483 * t429) * t368) * t433) * t411 / 0.2e1 + (m(2) * (t419 ^ 2 + t420 ^ 2) + Icges(2,3) + Icges(3,3)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
