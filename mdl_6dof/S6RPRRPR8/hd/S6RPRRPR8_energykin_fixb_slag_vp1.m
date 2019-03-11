% Calculate kinetic energy for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR8_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR8_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:22:34
% EndTime: 2019-03-09 05:22:37
% DurationCPUTime: 2.38s
% Computational Cost: add. (1306->292), mult. (1908->450), div. (0->0), fcn. (1898->10), ass. (0->140)
t490 = Icges(6,3) + Icges(5,3);
t438 = qJ(4) + pkin(10);
t428 = sin(t438);
t429 = cos(t438);
t445 = cos(qJ(1));
t441 = sin(qJ(3));
t442 = sin(qJ(1));
t473 = t441 * t442;
t382 = -t428 * t473 + t429 * t445;
t383 = t428 * t445 + t429 * t473;
t440 = sin(qJ(4));
t443 = cos(qJ(4));
t469 = t443 * t445;
t400 = -t440 * t473 + t469;
t471 = t442 * t443;
t474 = t440 * t445;
t401 = t441 * t471 + t474;
t444 = cos(qJ(3));
t470 = t442 * t444;
t489 = Icges(5,5) * t401 + Icges(6,5) * t383 + Icges(5,6) * t400 + Icges(6,6) * t382 - t490 * t470;
t472 = t441 * t445;
t384 = t428 * t472 + t429 * t442;
t385 = t428 * t442 - t429 * t472;
t402 = t440 * t472 + t471;
t475 = t440 * t442;
t403 = -t441 * t469 + t475;
t468 = t444 * t445;
t488 = Icges(5,5) * t403 + Icges(6,5) * t385 + Icges(5,6) * t402 + Icges(6,6) * t384 + t490 * t468;
t487 = (Icges(5,5) * t443 + Icges(6,5) * t429 - Icges(5,6) * t440 - Icges(6,6) * t428) * t444 + t490 * t441;
t467 = pkin(5) * t429;
t486 = -pkin(9) * t444 + t441 * t467;
t479 = t443 * pkin(4);
t485 = -qJ(5) * t444 + t441 * t479;
t477 = Icges(4,4) * t441;
t476 = Icges(4,4) * t444;
t411 = qJD(1) * (pkin(1) * t445 + qJ(2) * t442);
t466 = qJD(1) * t445 * pkin(7) + t411;
t433 = qJD(3) * t442;
t464 = qJD(4) * t444;
t407 = t445 * t464 + t433;
t434 = qJD(3) * t445;
t463 = qJD(5) * t444;
t423 = qJD(4) * t441 + qJD(1);
t462 = pkin(5) * t428;
t416 = pkin(1) * t442 - qJ(2) * t445;
t461 = -pkin(7) * t442 - t416;
t460 = pkin(3) * t441 - pkin(8) * t444;
t404 = t460 * t442;
t405 = t460 * t445;
t459 = -t404 * t433 - t405 * t434;
t458 = rSges(4,1) * t441 + rSges(4,2) * t444;
t457 = Icges(4,1) * t441 + t476;
t456 = Icges(4,2) * t444 + t477;
t455 = Icges(4,5) * t441 + Icges(4,6) * t444;
t390 = Icges(4,6) * t445 + t442 * t456;
t393 = Icges(4,5) * t445 + t442 * t457;
t454 = -t390 * t444 - t393 * t441;
t391 = Icges(4,6) * t442 - t445 * t456;
t394 = Icges(4,5) * t442 - t445 * t457;
t453 = t391 * t444 + t394 * t441;
t414 = -Icges(4,2) * t441 + t476;
t415 = Icges(4,1) * t444 - t477;
t452 = t414 * t444 + t415 * t441;
t360 = pkin(4) * t475 - t485 * t445;
t408 = -t442 * t464 + t434;
t451 = qJD(5) * t441 + t408 * t360 + t459;
t420 = pkin(3) * t444 + pkin(8) * t441;
t435 = qJD(2) * t442;
t450 = t420 * t433 + t435 + (t405 + t461) * qJD(1);
t449 = qJD(1) * t404 + (-qJD(3) * t420 - qJD(2)) * t445 + t466;
t359 = pkin(4) * t474 + t485 * t442;
t448 = t423 * t359 + t445 * t463 + t449;
t371 = qJ(5) * t441 + t444 * t479;
t447 = t407 * t371 - t442 * t463 + t450;
t430 = qJ(6) + t438;
t425 = cos(t430);
t424 = sin(t430);
t419 = rSges(2,1) * t445 - rSges(2,2) * t442;
t418 = rSges(4,1) * t444 - rSges(4,2) * t441;
t417 = rSges(2,1) * t442 + rSges(2,2) * t445;
t413 = Icges(4,5) * t444 - Icges(4,6) * t441;
t410 = qJD(6) * t441 + t423;
t397 = rSges(4,3) * t442 - t445 * t458;
t396 = rSges(5,3) * t441 + (rSges(5,1) * t443 - rSges(5,2) * t440) * t444;
t395 = rSges(4,3) * t445 + t442 * t458;
t392 = Icges(5,5) * t441 + (Icges(5,1) * t443 - Icges(5,4) * t440) * t444;
t389 = Icges(5,6) * t441 + (Icges(5,4) * t443 - Icges(5,2) * t440) * t444;
t388 = Icges(4,3) * t442 - t445 * t455;
t387 = Icges(4,3) * t445 + t442 * t455;
t381 = t434 + (-qJD(4) - qJD(6)) * t470;
t380 = qJD(6) * t468 + t407;
t379 = t424 * t442 - t425 * t472;
t378 = t424 * t472 + t425 * t442;
t377 = t424 * t445 + t425 * t473;
t376 = -t424 * t473 + t425 * t445;
t375 = rSges(6,3) * t441 + (rSges(6,1) * t429 - rSges(6,2) * t428) * t444;
t374 = Icges(6,5) * t441 + (Icges(6,1) * t429 - Icges(6,4) * t428) * t444;
t373 = Icges(6,6) * t441 + (Icges(6,4) * t429 - Icges(6,2) * t428) * t444;
t370 = rSges(7,3) * t441 + (rSges(7,1) * t425 - rSges(7,2) * t424) * t444;
t369 = Icges(7,5) * t441 + (Icges(7,1) * t425 - Icges(7,4) * t424) * t444;
t368 = Icges(7,6) * t441 + (Icges(7,4) * t425 - Icges(7,2) * t424) * t444;
t367 = Icges(7,3) * t441 + (Icges(7,5) * t425 - Icges(7,6) * t424) * t444;
t366 = t411 - qJD(2) * t445 + qJD(1) * (-rSges(3,2) * t445 + rSges(3,3) * t442);
t365 = t435 + (rSges(3,2) * t442 + rSges(3,3) * t445 - t416) * qJD(1);
t363 = pkin(9) * t441 + t444 * t467;
t362 = rSges(5,1) * t403 + rSges(5,2) * t402 + rSges(5,3) * t468;
t361 = rSges(5,1) * t401 + rSges(5,2) * t400 - rSges(5,3) * t470;
t358 = Icges(5,1) * t403 + Icges(5,4) * t402 + Icges(5,5) * t468;
t357 = Icges(5,1) * t401 + Icges(5,4) * t400 - Icges(5,5) * t470;
t356 = Icges(5,4) * t403 + Icges(5,2) * t402 + Icges(5,6) * t468;
t355 = Icges(5,4) * t401 + Icges(5,2) * t400 - Icges(5,6) * t470;
t352 = (-t395 * t442 + t397 * t445) * qJD(3);
t350 = rSges(6,1) * t385 + rSges(6,2) * t384 + rSges(6,3) * t468;
t349 = rSges(6,1) * t383 + rSges(6,2) * t382 - rSges(6,3) * t470;
t348 = Icges(6,1) * t385 + Icges(6,4) * t384 + Icges(6,5) * t468;
t347 = Icges(6,1) * t383 + Icges(6,4) * t382 - Icges(6,5) * t470;
t346 = Icges(6,4) * t385 + Icges(6,2) * t384 + Icges(6,6) * t468;
t345 = Icges(6,4) * t383 + Icges(6,2) * t382 - Icges(6,6) * t470;
t341 = rSges(7,1) * t379 + rSges(7,2) * t378 + rSges(7,3) * t468;
t340 = rSges(7,1) * t377 + rSges(7,2) * t376 - rSges(7,3) * t470;
t339 = Icges(7,1) * t379 + Icges(7,4) * t378 + Icges(7,5) * t468;
t338 = Icges(7,1) * t377 + Icges(7,4) * t376 - Icges(7,5) * t470;
t337 = Icges(7,4) * t379 + Icges(7,2) * t378 + Icges(7,6) * t468;
t336 = Icges(7,4) * t377 + Icges(7,2) * t376 - Icges(7,6) * t470;
t335 = Icges(7,5) * t379 + Icges(7,6) * t378 + Icges(7,3) * t468;
t334 = Icges(7,5) * t377 + Icges(7,6) * t376 - Icges(7,3) * t470;
t333 = qJD(1) * t395 + (-qJD(3) * t418 - qJD(2)) * t445 + t466;
t332 = t418 * t433 + t435 + (-t397 + t461) * qJD(1);
t331 = t462 * t442 - t486 * t445;
t330 = t486 * t442 + t462 * t445;
t329 = t361 * t423 - t396 * t408 + t449;
t328 = -t362 * t423 + t396 * t407 + t450;
t327 = -t361 * t407 + t362 * t408 + t459;
t326 = t349 * t423 + (-t371 - t375) * t408 + t448;
t325 = t375 * t407 + (-t350 - t360) * t423 + t447;
t324 = t350 * t408 + (-t349 - t359) * t407 + t451;
t323 = t330 * t423 + t340 * t410 - t370 * t381 + (-t363 - t371) * t408 + t448;
t322 = -t341 * t410 + t363 * t407 + t370 * t380 + (-t331 - t360) * t423 + t447;
t321 = t331 * t408 - t340 * t380 + t341 * t381 + (-t330 - t359) * t407 + t451;
t1 = ((t445 * t413 + t442 * t452) * qJD(1) + (t445 ^ 2 * t387 + (t453 * t442 + (t388 - t454) * t445) * t442) * qJD(3)) * t434 / 0.2e1 + ((t442 * t413 - t445 * t452) * qJD(1) + (t442 ^ 2 * t388 + (t454 * t445 + (t387 - t453) * t442) * t445) * qJD(3)) * t433 / 0.2e1 + qJD(1) * ((-t441 * t414 + t444 * t415) * qJD(1) + ((-t390 * t441 + t393 * t444) * t445 + (-t391 * t441 + t394 * t444) * t442) * qJD(3)) / 0.2e1 + m(7) * (t321 ^ 2 + t322 ^ 2 + t323 ^ 2) / 0.2e1 + m(5) * (t327 ^ 2 + t328 ^ 2 + t329 ^ 2) / 0.2e1 + m(6) * (t324 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + m(4) * (t332 ^ 2 + t333 ^ 2 + t352 ^ 2) / 0.2e1 + m(3) * (t365 ^ 2 + t366 ^ 2) / 0.2e1 + t410 * ((t334 * t381 + t335 * t380 + t367 * t410) * t441 + ((-t336 * t424 + t338 * t425) * t381 + (-t337 * t424 + t339 * t425) * t380 + (-t368 * t424 + t369 * t425) * t410) * t444) / 0.2e1 + t381 * ((-t334 * t470 + t376 * t336 + t377 * t338) * t381 + (-t335 * t470 + t337 * t376 + t339 * t377) * t380 + (-t367 * t470 + t368 * t376 + t369 * t377) * t410) / 0.2e1 + t380 * ((t334 * t468 + t336 * t378 + t338 * t379) * t381 + (t335 * t468 + t378 * t337 + t379 * t339) * t380 + (t367 * t468 + t368 * t378 + t369 * t379) * t410) / 0.2e1 + ((t373 * t384 + t374 * t385 + t389 * t402 + t392 * t403 + t487 * t468) * t423 + (t345 * t384 + t347 * t385 + t355 * t402 + t357 * t403 + t489 * t468) * t408 + (t384 * t346 + t385 * t348 + t402 * t356 + t403 * t358 + t488 * t468) * t407) * t407 / 0.2e1 + ((t373 * t382 + t374 * t383 + t389 * t400 + t392 * t401 - t487 * t470) * t423 + (t382 * t345 + t383 * t347 + t400 * t355 + t401 * t357 - t489 * t470) * t408 + (t346 * t382 + t348 * t383 + t356 * t400 + t358 * t401 - t488 * t470) * t407) * t408 / 0.2e1 + (((-t373 * t428 + t374 * t429 - t389 * t440 + t392 * t443) * t423 + (-t345 * t428 + t347 * t429 - t355 * t440 + t357 * t443) * t408 + (-t346 * t428 + t348 * t429 - t356 * t440 + t358 * t443) * t407) * t444 + (t488 * t407 + t489 * t408 + t487 * t423) * t441) * t423 / 0.2e1 + (Icges(2,3) + m(2) * (t417 ^ 2 + t419 ^ 2) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
