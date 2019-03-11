% Calculate kinetic energy for
% S6RPRPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR2_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_energykin_fixb_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_energykin_fixb_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_energykin_fixb_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_energykin_fixb_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR2_energykin_fixb_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR2_energykin_fixb_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:36:56
% EndTime: 2019-03-09 03:36:58
% DurationCPUTime: 2.13s
% Computational Cost: add. (1954->267), mult. (1563->419), div. (0->0), fcn. (1496->12), ass. (0->142)
t499 = Icges(4,3) + Icges(5,3);
t427 = qJ(3) + pkin(11);
t421 = sin(t427);
t423 = cos(t427);
t432 = sin(qJ(3));
t435 = cos(qJ(3));
t498 = Icges(4,5) * t435 + Icges(5,5) * t423 - Icges(4,6) * t432 - Icges(5,6) * t421;
t428 = qJ(1) + pkin(10);
t422 = sin(t428);
t424 = cos(t428);
t497 = t499 * t422 + t498 * t424;
t496 = t498 * t422 - t499 * t424;
t495 = Icges(4,5) * t432 + Icges(5,5) * t421 + Icges(4,6) * t435 + Icges(5,6) * t423;
t479 = Icges(5,4) * t421;
t402 = Icges(5,2) * t423 + t479;
t478 = Icges(5,4) * t423;
t403 = Icges(5,1) * t421 + t478;
t481 = Icges(4,4) * t432;
t410 = Icges(4,2) * t435 + t481;
t480 = Icges(4,4) * t435;
t411 = Icges(4,1) * t432 + t480;
t494 = -t402 * t421 + t403 * t423 - t410 * t432 + t411 * t435;
t451 = -Icges(5,2) * t421 + t478;
t361 = Icges(5,6) * t422 + t424 * t451;
t453 = Icges(5,1) * t423 - t479;
t363 = Icges(5,5) * t422 + t424 * t453;
t452 = -Icges(4,2) * t432 + t480;
t377 = Icges(4,6) * t422 + t424 * t452;
t454 = Icges(4,1) * t435 - t481;
t380 = Icges(4,5) * t422 + t424 * t454;
t493 = -t361 * t421 + t363 * t423 - t377 * t432 + t380 * t435;
t360 = -Icges(5,6) * t424 + t422 * t451;
t362 = -Icges(5,5) * t424 + t422 * t453;
t376 = -Icges(4,6) * t424 + t422 * t452;
t379 = -Icges(4,5) * t424 + t422 * t454;
t492 = t360 * t421 - t362 * t423 + t376 * t432 - t379 * t435;
t433 = sin(qJ(1));
t488 = pkin(1) * t433;
t487 = pkin(3) * t432;
t485 = pkin(3) * t435;
t434 = cos(qJ(5));
t484 = pkin(5) * t434;
t477 = t421 * t422;
t476 = t421 * t424;
t431 = sin(qJ(5));
t475 = t422 * t431;
t474 = t423 * t422;
t429 = qJ(5) + qJ(6);
t425 = sin(t429);
t473 = t424 * t425;
t426 = cos(t429);
t472 = t424 * t426;
t471 = t424 * t431;
t470 = t424 * t434;
t436 = cos(qJ(1));
t420 = qJD(1) * t436 * pkin(1);
t469 = qJD(1) * (pkin(2) * t424 + pkin(7) * t422) + t420;
t417 = qJD(3) * t422;
t467 = qJD(5) * t421;
t397 = t424 * t467 + t417;
t468 = qJD(3) * t424;
t466 = qJD(6) * t421;
t353 = -qJ(4) * t424 + t422 * t485;
t354 = qJ(4) * t422 + t424 * t485;
t465 = t353 * t417 + t354 * t468 + qJD(2);
t462 = -pkin(2) * t422 + pkin(7) * t424 - t488;
t398 = t422 * t467 - t468;
t461 = -t353 + t462;
t460 = pkin(4) * t423 + pkin(8) * t421;
t459 = rSges(4,1) * t435 - rSges(4,2) * t432;
t458 = rSges(5,1) * t423 - rSges(5,2) * t421;
t457 = qJD(3) * (-rSges(5,1) * t421 - rSges(5,2) * t423 - t487);
t456 = qJD(3) * (-pkin(4) * t421 + pkin(8) * t423 - t487);
t389 = t460 * t422;
t390 = t460 * t424;
t455 = t389 * t417 + t390 * t468 + t465;
t442 = qJD(1) * t354 - qJD(4) * t424 + t469;
t441 = pkin(9) * t421 + t423 * t484;
t440 = qJD(1) * t390 + t422 * t456 + t442;
t416 = qJD(4) * t422;
t439 = t416 + (-t389 + t461) * qJD(1) + t424 * t456;
t415 = -qJD(5) * t423 + qJD(1);
t414 = rSges(2,1) * t436 - rSges(2,2) * t433;
t413 = rSges(2,1) * t433 + rSges(2,2) * t436;
t412 = rSges(4,1) * t432 + rSges(4,2) * t435;
t399 = qJD(1) + (-qJD(5) - qJD(6)) * t423;
t396 = t423 * t470 + t475;
t395 = t422 * t434 - t423 * t471;
t394 = t434 * t474 - t471;
t393 = -t431 * t474 - t470;
t392 = t420 + qJD(1) * (rSges(3,1) * t424 - rSges(3,2) * t422);
t391 = (-rSges(3,1) * t422 - rSges(3,2) * t424 - t488) * qJD(1);
t387 = t422 * t425 + t423 * t472;
t386 = t422 * t426 - t423 * t473;
t385 = t426 * t474 - t473;
t384 = -t425 * t474 - t472;
t383 = rSges(4,3) * t422 + t424 * t459;
t382 = -rSges(4,3) * t424 + t422 * t459;
t381 = -rSges(6,3) * t423 + (rSges(6,1) * t434 - rSges(6,2) * t431) * t421;
t378 = -Icges(6,5) * t423 + (Icges(6,1) * t434 - Icges(6,4) * t431) * t421;
t375 = -Icges(6,6) * t423 + (Icges(6,4) * t434 - Icges(6,2) * t431) * t421;
t372 = -Icges(6,3) * t423 + (Icges(6,5) * t434 - Icges(6,6) * t431) * t421;
t371 = -rSges(7,3) * t423 + (rSges(7,1) * t426 - rSges(7,2) * t425) * t421;
t370 = -Icges(7,5) * t423 + (Icges(7,1) * t426 - Icges(7,4) * t425) * t421;
t369 = -Icges(7,6) * t423 + (Icges(7,4) * t426 - Icges(7,2) * t425) * t421;
t368 = -Icges(7,3) * t423 + (Icges(7,5) * t426 - Icges(7,6) * t425) * t421;
t365 = rSges(5,3) * t422 + t424 * t458;
t364 = -rSges(5,3) * t424 + t422 * t458;
t357 = t422 * t466 + t398;
t356 = t424 * t466 + t397;
t355 = -pkin(9) * t423 + t421 * t484;
t349 = rSges(6,1) * t396 + rSges(6,2) * t395 + rSges(6,3) * t476;
t348 = rSges(6,1) * t394 + rSges(6,2) * t393 + rSges(6,3) * t477;
t347 = Icges(6,1) * t396 + Icges(6,4) * t395 + Icges(6,5) * t476;
t346 = Icges(6,1) * t394 + Icges(6,4) * t393 + Icges(6,5) * t477;
t345 = Icges(6,4) * t396 + Icges(6,2) * t395 + Icges(6,6) * t476;
t344 = Icges(6,4) * t394 + Icges(6,2) * t393 + Icges(6,6) * t477;
t343 = Icges(6,5) * t396 + Icges(6,6) * t395 + Icges(6,3) * t476;
t342 = Icges(6,5) * t394 + Icges(6,6) * t393 + Icges(6,3) * t477;
t341 = pkin(5) * t475 + t424 * t441;
t340 = -pkin(5) * t471 + t422 * t441;
t339 = qJD(1) * t383 - t412 * t417 + t469;
t338 = -t412 * t468 + (-t382 + t462) * qJD(1);
t337 = rSges(7,1) * t387 + rSges(7,2) * t386 + rSges(7,3) * t476;
t336 = rSges(7,1) * t385 + rSges(7,2) * t384 + rSges(7,3) * t477;
t335 = qJD(2) + (t382 * t422 + t383 * t424) * qJD(3);
t334 = Icges(7,1) * t387 + Icges(7,4) * t386 + Icges(7,5) * t476;
t333 = Icges(7,1) * t385 + Icges(7,4) * t384 + Icges(7,5) * t477;
t332 = Icges(7,4) * t387 + Icges(7,2) * t386 + Icges(7,6) * t476;
t331 = Icges(7,4) * t385 + Icges(7,2) * t384 + Icges(7,6) * t477;
t330 = Icges(7,5) * t387 + Icges(7,6) * t386 + Icges(7,3) * t476;
t329 = Icges(7,5) * t385 + Icges(7,6) * t384 + Icges(7,3) * t477;
t328 = qJD(1) * t365 + t422 * t457 + t442;
t327 = t416 + t424 * t457 + (-t364 + t461) * qJD(1);
t326 = (t364 * t422 + t365 * t424) * qJD(3) + t465;
t325 = t349 * t415 - t381 * t397 + t440;
t324 = -t348 * t415 + t381 * t398 + t439;
t323 = t348 * t397 - t349 * t398 + t455;
t322 = t337 * t399 + t341 * t415 - t355 * t397 - t356 * t371 + t440;
t321 = -t336 * t399 - t340 * t415 + t355 * t398 + t357 * t371 + t439;
t320 = t336 * t356 - t337 * t357 + t340 * t397 - t341 * t398 + t455;
t1 = t398 * ((t343 * t477 + t345 * t393 + t347 * t394) * t397 + (t342 * t477 + t393 * t344 + t394 * t346) * t398 + (t372 * t477 + t375 * t393 + t378 * t394) * t415) / 0.2e1 + t415 * ((-t342 * t398 - t343 * t397 - t372 * t415) * t423 + ((-t345 * t431 + t347 * t434) * t397 + (-t344 * t431 + t346 * t434) * t398 + (-t375 * t431 + t378 * t434) * t415) * t421) / 0.2e1 + t356 * ((t330 * t476 + t386 * t332 + t387 * t334) * t356 + (t329 * t476 + t331 * t386 + t333 * t387) * t357 + (t368 * t476 + t369 * t386 + t370 * t387) * t399) / 0.2e1 + t357 * ((t330 * t477 + t332 * t384 + t334 * t385) * t356 + (t329 * t477 + t384 * t331 + t385 * t333) * t357 + (t368 * t477 + t369 * t384 + t370 * t385) * t399) / 0.2e1 + t399 * ((-t329 * t357 - t330 * t356 - t368 * t399) * t423 + ((-t332 * t425 + t334 * t426) * t356 + (-t331 * t425 + t333 * t426) * t357 + (-t369 * t425 + t370 * t426) * t399) * t421) / 0.2e1 + m(5) * (t326 ^ 2 + t327 ^ 2 + t328 ^ 2) / 0.2e1 + m(6) * (t323 ^ 2 + t324 ^ 2 + t325 ^ 2) / 0.2e1 + m(7) * (t320 ^ 2 + t321 ^ 2 + t322 ^ 2) / 0.2e1 + t397 * ((t343 * t476 + t395 * t345 + t396 * t347) * t397 + (t342 * t476 + t344 * t395 + t346 * t396) * t398 + (t372 * t476 + t375 * t395 + t378 * t396) * t415) / 0.2e1 + m(3) * (qJD(2) ^ 2 + t391 ^ 2 + t392 ^ 2) / 0.2e1 + m(4) * (t335 ^ 2 + t338 ^ 2 + t339 ^ 2) / 0.2e1 + (((-t360 * t423 - t362 * t421 - t376 * t435 - t379 * t432) * t424 + (t361 * t423 + t363 * t421 + t377 * t435 + t380 * t432) * t422) * qJD(3) + (t423 * t402 + t421 * t403 + t435 * t410 + t432 * t411) * qJD(1)) * qJD(1) / 0.2e1 + ((t497 * t422 ^ 2 + (t492 * t424 + (t493 - t496) * t422) * t424) * qJD(3) + (t422 * t495 + t424 * t494) * qJD(1)) * t417 / 0.2e1 - ((t496 * t424 ^ 2 + (t493 * t422 + (t492 - t497) * t424) * t422) * qJD(3) + (t422 * t494 - t424 * t495) * qJD(1)) * t468 / 0.2e1 + (Icges(2,3) + Icges(3,3) + m(2) * (t413 ^ 2 + t414 ^ 2)) * qJD(1) ^ 2 / 0.2e1;
T  = t1;
