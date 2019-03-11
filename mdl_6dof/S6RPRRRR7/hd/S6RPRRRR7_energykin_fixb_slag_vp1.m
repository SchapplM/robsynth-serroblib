% Calculate kinetic energy for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR7_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR7_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR7_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR7_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:16:55
% EndTime: 2019-03-09 07:16:56
% DurationCPUTime: 1.54s
% Computational Cost: add. (1228->248), mult. (1354->400), div. (0->0), fcn. (1232->10), ass. (0->141)
t420 = sin(qJ(1));
t423 = cos(qJ(1));
t417 = qJ(3) + qJ(4);
t415 = qJ(5) + t417;
t404 = cos(t415);
t403 = sin(t415);
t471 = Icges(6,4) * t403;
t444 = Icges(6,2) * t404 + t471;
t344 = Icges(6,6) * t423 + t444 * t420;
t345 = Icges(6,6) * t420 - t444 * t423;
t470 = Icges(6,4) * t404;
t447 = Icges(6,1) * t403 + t470;
t346 = Icges(6,5) * t423 + t447 * t420;
t347 = Icges(6,5) * t420 - t447 * t423;
t380 = -Icges(6,2) * t403 + t470;
t381 = Icges(6,1) * t404 - t471;
t410 = qJD(3) * t420;
t392 = qJD(4) * t420 + t410;
t384 = qJD(5) * t420 + t392;
t411 = qJD(3) * t423;
t393 = qJD(4) * t423 + t411;
t385 = qJD(5) * t423 + t393;
t482 = (t344 * t404 + t346 * t403) * t385 + (t345 * t404 + t347 * t403) * t384 + (t380 * t404 + t381 * t403) * qJD(1);
t414 = cos(t417);
t413 = sin(t417);
t473 = Icges(5,4) * t413;
t445 = Icges(5,2) * t414 + t473;
t354 = Icges(5,6) * t423 + t445 * t420;
t355 = Icges(5,6) * t420 - t445 * t423;
t472 = Icges(5,4) * t414;
t448 = Icges(5,1) * t413 + t472;
t356 = Icges(5,5) * t423 + t448 * t420;
t357 = Icges(5,5) * t420 - t448 * t423;
t387 = -Icges(5,2) * t413 + t472;
t388 = Icges(5,1) * t414 - t473;
t481 = (t354 * t414 + t356 * t413) * t393 + (t355 * t414 + t357 * t413) * t392 + (t387 * t414 + t388 * t413) * qJD(1);
t419 = sin(qJ(3));
t479 = pkin(3) * t419;
t478 = pkin(4) * t414;
t475 = Icges(4,4) * t419;
t422 = cos(qJ(3));
t474 = Icges(4,4) * t422;
t469 = t404 * t420;
t468 = t404 * t423;
t418 = sin(qJ(6));
t467 = t418 * t420;
t466 = t418 * t423;
t421 = cos(qJ(6));
t465 = t420 * t421;
t464 = t421 * t423;
t390 = qJD(1) * (pkin(1) * t423 + qJ(2) * t420);
t463 = qJD(1) * t423 * pkin(7) + t390;
t412 = qJD(2) * t420;
t459 = pkin(3) * qJD(3) * t422;
t462 = t420 * t459 + t412;
t460 = qJD(6) * t404;
t458 = t392 * t478 + t462;
t457 = pkin(4) * t413;
t398 = pkin(1) * t420 - qJ(2) * t423;
t456 = -pkin(7) * t420 - t398;
t376 = pkin(8) * t420 - t423 * t479;
t455 = -t376 + t456;
t454 = pkin(5) * t403 - pkin(10) * t404;
t377 = pkin(8) * t423 + t420 * t479;
t453 = t376 * t411 - t377 * t410;
t452 = rSges(4,1) * t419 + rSges(4,2) * t422;
t451 = rSges(5,1) * t413 + rSges(5,2) * t414;
t450 = rSges(6,1) * t403 + rSges(6,2) * t404;
t449 = Icges(4,1) * t419 + t474;
t446 = Icges(4,2) * t422 + t475;
t443 = Icges(4,5) * t419 + Icges(4,6) * t422;
t442 = Icges(5,5) * t413 + Icges(5,6) * t414;
t441 = Icges(6,5) * t403 + Icges(6,6) * t404;
t365 = Icges(4,6) * t423 + t446 * t420;
t367 = Icges(4,5) * t423 + t449 * t420;
t436 = -t365 * t422 - t367 * t419;
t366 = Icges(4,6) * t420 - t446 * t423;
t368 = Icges(4,5) * t420 - t449 * t423;
t435 = t366 * t422 + t368 * t419;
t396 = -Icges(4,2) * t419 + t474;
t397 = Icges(4,1) * t422 - t475;
t432 = t396 * t422 + t397 * t419;
t334 = pkin(9) * t420 - t457 * t423;
t431 = -t334 + t455;
t430 = (Icges(6,5) * t404 - Icges(6,6) * t403) * qJD(1) + (Icges(6,3) * t423 + t441 * t420) * t385 + (Icges(6,3) * t420 - t441 * t423) * t384;
t429 = (Icges(5,5) * t414 - Icges(5,6) * t413) * qJD(1) + (Icges(5,3) * t423 + t442 * t420) * t393 + (Icges(5,3) * t420 - t442 * t423) * t392;
t335 = pkin(9) * t423 + t457 * t420;
t428 = t393 * t334 - t335 * t392 + t453;
t427 = qJD(1) * t377 + (-qJD(2) - t459) * t423 + t463;
t426 = qJD(1) * t335 - t393 * t478 + t427;
t401 = rSges(2,1) * t423 - rSges(2,2) * t420;
t400 = rSges(4,1) * t422 - rSges(4,2) * t419;
t399 = rSges(2,1) * t420 + rSges(2,2) * t423;
t395 = Icges(4,5) * t422 - Icges(4,6) * t419;
t394 = qJD(6) * t403 + qJD(1);
t389 = rSges(5,1) * t414 - rSges(5,2) * t413;
t383 = pkin(5) * t404 + pkin(10) * t403;
t382 = rSges(6,1) * t404 - rSges(6,2) * t403;
t375 = -t403 * t464 + t467;
t374 = t403 * t466 + t465;
t373 = t403 * t465 + t466;
t372 = -t403 * t467 + t464;
t371 = rSges(4,3) * t420 - t452 * t423;
t370 = rSges(4,3) * t423 + t452 * t420;
t364 = Icges(4,3) * t420 - t443 * t423;
t363 = Icges(4,3) * t423 + t443 * t420;
t362 = t454 * t423;
t361 = t454 * t420;
t359 = rSges(5,3) * t420 - t451 * t423;
t358 = rSges(5,3) * t423 + t451 * t420;
t351 = rSges(6,3) * t420 - t450 * t423;
t350 = rSges(6,3) * t423 + t450 * t420;
t349 = -t420 * t460 + t385;
t348 = t423 * t460 + t384;
t341 = t390 - qJD(2) * t423 + qJD(1) * (-rSges(3,2) * t423 + rSges(3,3) * t420);
t340 = t412 + (rSges(3,2) * t420 + rSges(3,3) * t423 - t398) * qJD(1);
t339 = rSges(7,3) * t403 + (rSges(7,1) * t421 - rSges(7,2) * t418) * t404;
t338 = Icges(7,5) * t403 + (Icges(7,1) * t421 - Icges(7,4) * t418) * t404;
t337 = Icges(7,6) * t403 + (Icges(7,4) * t421 - Icges(7,2) * t418) * t404;
t336 = Icges(7,3) * t403 + (Icges(7,5) * t421 - Icges(7,6) * t418) * t404;
t331 = (-t370 * t420 + t371 * t423) * qJD(3);
t330 = rSges(7,1) * t375 + rSges(7,2) * t374 + rSges(7,3) * t468;
t329 = rSges(7,1) * t373 + rSges(7,2) * t372 - rSges(7,3) * t469;
t328 = Icges(7,1) * t375 + Icges(7,4) * t374 + Icges(7,5) * t468;
t327 = Icges(7,1) * t373 + Icges(7,4) * t372 - Icges(7,5) * t469;
t326 = Icges(7,4) * t375 + Icges(7,2) * t374 + Icges(7,6) * t468;
t325 = Icges(7,4) * t373 + Icges(7,2) * t372 - Icges(7,6) * t469;
t324 = Icges(7,5) * t375 + Icges(7,6) * t374 + Icges(7,3) * t468;
t323 = Icges(7,5) * t373 + Icges(7,6) * t372 - Icges(7,3) * t469;
t322 = qJD(1) * t370 + (-qJD(3) * t400 - qJD(2)) * t423 + t463;
t321 = t400 * t410 + t412 + (-t371 + t456) * qJD(1);
t320 = qJD(1) * t358 - t389 * t393 + t427;
t319 = t389 * t392 + (-t359 + t455) * qJD(1) + t462;
t318 = -t358 * t392 + t359 * t393 + t453;
t317 = qJD(1) * t350 - t382 * t385 + t426;
t316 = t382 * t384 + (-t351 + t431) * qJD(1) + t458;
t315 = -t350 * t384 + t351 * t385 + t428;
t314 = qJD(1) * t361 + t329 * t394 - t339 * t349 - t383 * t385 + t426;
t313 = -t330 * t394 + t339 * t348 + t383 * t384 + (t362 + t431) * qJD(1) + t458;
t312 = -t329 * t348 + t330 * t349 - t361 * t384 - t362 * t385 + t428;
t1 = m(6) * (t315 ^ 2 + t316 ^ 2 + t317 ^ 2) / 0.2e1 + m(5) * (t318 ^ 2 + t319 ^ 2 + t320 ^ 2) / 0.2e1 + m(4) * (t321 ^ 2 + t322 ^ 2 + t331 ^ 2) / 0.2e1 + m(3) * (t340 ^ 2 + t341 ^ 2) / 0.2e1 + t394 * ((t323 * t349 + t324 * t348 + t336 * t394) * t403 + ((-t325 * t418 + t327 * t421) * t349 + (-t326 * t418 + t328 * t421) * t348 + (-t337 * t418 + t338 * t421) * t394) * t404) / 0.2e1 + t349 * ((-t323 * t469 + t372 * t325 + t373 * t327) * t349 + (-t324 * t469 + t326 * t372 + t328 * t373) * t348 + (-t336 * t469 + t337 * t372 + t338 * t373) * t394) / 0.2e1 + t348 * ((t323 * t468 + t325 * t374 + t327 * t375) * t349 + (t324 * t468 + t374 * t326 + t375 * t328) * t348 + (t336 * t468 + t337 * t374 + t338 * t375) * t394) / 0.2e1 + t385 * (t482 * t420 + t430 * t423) / 0.2e1 + t384 * (t430 * t420 - t482 * t423) / 0.2e1 + t393 * (t481 * t420 + t429 * t423) / 0.2e1 + t392 * (t429 * t420 - t481 * t423) / 0.2e1 + ((t420 * t395 - t432 * t423) * qJD(1) + (t420 ^ 2 * t364 + (t436 * t423 + (t363 - t435) * t420) * t423) * qJD(3)) * t410 / 0.2e1 + m(7) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + ((t423 * t395 + t432 * t420) * qJD(1) + (t423 ^ 2 * t363 + (t435 * t420 + (t364 - t436) * t423) * t420) * qJD(3)) * t411 / 0.2e1 + (m(2) * (t399 ^ 2 + t401 ^ 2) + Icges(2,3) + Icges(3,1)) * qJD(1) ^ 2 / 0.2e1 + ((-t344 * t403 + t346 * t404) * t385 + (-t345 * t403 + t347 * t404) * t384 + (-t354 * t413 + t356 * t414) * t393 + (-t355 * t413 + t357 * t414) * t392 + ((-t365 * t419 + t367 * t422) * t423 + (-t366 * t419 + t368 * t422) * t420) * qJD(3) + (-t403 * t380 + t404 * t381 - t413 * t387 + t414 * t388 - t419 * t396 + t422 * t397) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
