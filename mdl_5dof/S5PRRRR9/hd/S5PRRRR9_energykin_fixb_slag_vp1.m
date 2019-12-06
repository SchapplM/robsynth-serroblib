% Calculate kinetic energy for
% S5PRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR9_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_energykin_fixb_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR9_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR9_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR9_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:19:01
% EndTime: 2019-12-05 17:19:02
% DurationCPUTime: 1.62s
% Computational Cost: add. (1982->298), mult. (4597->477), div. (0->0), fcn. (5634->12), ass. (0->136)
t445 = qJD(2) ^ 2;
t444 = cos(qJ(3));
t420 = cos(qJ(4));
t443 = pkin(4) * t420;
t413 = sin(pkin(10));
t415 = cos(pkin(10));
t419 = sin(qJ(2));
t416 = cos(pkin(5));
t421 = cos(qJ(2));
t434 = t416 * t421;
t396 = t413 * t419 - t415 * t434;
t417 = sin(qJ(4));
t441 = t396 * t417;
t398 = t413 * t434 + t415 * t419;
t440 = t398 * t417;
t414 = sin(pkin(5));
t439 = t413 * t414;
t418 = sin(qJ(3));
t438 = t414 * t418;
t437 = t414 * t421;
t436 = t415 * t414;
t435 = t416 * t419;
t433 = qJD(2) * t414;
t406 = t413 * t433;
t387 = qJD(3) * t398 + t406;
t409 = qJD(2) * t416;
t431 = t417 * t437;
t399 = -t413 * t435 + t415 * t421;
t430 = t414 * t444;
t383 = t399 * t418 - t413 * t430;
t347 = qJD(4) * t383 + t387;
t429 = t415 * t433;
t397 = t413 * t421 + t415 * t435;
t373 = pkin(2) * t397 + pkin(7) * t396;
t374 = pkin(2) * t399 + pkin(7) * t398;
t428 = t373 * t406 + t374 * t429 + qJD(1);
t388 = qJD(3) * t396 - t429;
t403 = -qJD(3) * t437 + t409;
t402 = (pkin(2) * t419 - pkin(7) * t421) * t414;
t427 = t374 * t409 - t402 * t406;
t381 = t397 * t418 + t415 * t430;
t348 = qJD(4) * t381 + t388;
t400 = -t416 * t444 + t419 * t438;
t380 = qJD(4) * t400 + t403;
t382 = t397 * t444 - t418 * t436;
t345 = t382 * pkin(3) + t381 * pkin(8);
t384 = t399 * t444 + t413 * t438;
t346 = t384 * pkin(3) + t383 * pkin(8);
t426 = t387 * t345 - t346 * t388 + t428;
t425 = (-t373 * t416 - t402 * t436) * qJD(2);
t401 = t416 * t418 + t419 * t430;
t375 = t401 * pkin(3) + t400 * pkin(8);
t424 = t403 * t346 - t375 * t387 + t427;
t423 = -t345 * t403 + t388 * t375 + t425;
t412 = qJ(4) + qJ(5);
t411 = cos(t412);
t410 = sin(t412);
t392 = rSges(3,3) * t416 + (rSges(3,1) * t419 + rSges(3,2) * t421) * t414;
t391 = Icges(3,5) * t416 + (Icges(3,1) * t419 + Icges(3,4) * t421) * t414;
t390 = Icges(3,6) * t416 + (Icges(3,4) * t419 + Icges(3,2) * t421) * t414;
t389 = Icges(3,3) * t416 + (Icges(3,5) * t419 + Icges(3,6) * t421) * t414;
t386 = t401 * t420 - t431;
t385 = -t401 * t417 - t420 * t437;
t377 = t401 * t411 - t410 * t437;
t376 = -t401 * t410 - t411 * t437;
t371 = rSges(4,1) * t401 - rSges(4,2) * t400 - rSges(4,3) * t437;
t370 = Icges(4,1) * t401 - Icges(4,4) * t400 - Icges(4,5) * t437;
t369 = Icges(4,4) * t401 - Icges(4,2) * t400 - Icges(4,6) * t437;
t368 = Icges(4,5) * t401 - Icges(4,6) * t400 - Icges(4,3) * t437;
t365 = rSges(3,1) * t399 - rSges(3,2) * t398 + rSges(3,3) * t439;
t364 = rSges(3,1) * t397 - rSges(3,2) * t396 - rSges(3,3) * t436;
t363 = Icges(3,1) * t399 - Icges(3,4) * t398 + Icges(3,5) * t439;
t362 = Icges(3,1) * t397 - Icges(3,4) * t396 - Icges(3,5) * t436;
t361 = Icges(3,4) * t399 - Icges(3,2) * t398 + Icges(3,6) * t439;
t360 = Icges(3,4) * t397 - Icges(3,2) * t396 - Icges(3,6) * t436;
t359 = Icges(3,5) * t399 - Icges(3,6) * t398 + Icges(3,3) * t439;
t358 = Icges(3,5) * t397 - Icges(3,6) * t396 - Icges(3,3) * t436;
t357 = qJD(5) * t400 + t380;
t356 = t384 * t420 + t440;
t355 = -t384 * t417 + t398 * t420;
t354 = t382 * t420 + t441;
t353 = -t382 * t417 + t396 * t420;
t352 = t384 * t411 + t398 * t410;
t351 = -t384 * t410 + t398 * t411;
t350 = t382 * t411 + t396 * t410;
t349 = -t382 * t410 + t396 * t411;
t342 = (-t364 * t416 - t392 * t436) * qJD(2);
t341 = (t365 * t416 - t392 * t439) * qJD(2);
t340 = rSges(5,1) * t386 + rSges(5,2) * t385 + rSges(5,3) * t400;
t339 = Icges(5,1) * t386 + Icges(5,4) * t385 + Icges(5,5) * t400;
t338 = Icges(5,4) * t386 + Icges(5,2) * t385 + Icges(5,6) * t400;
t337 = Icges(5,5) * t386 + Icges(5,6) * t385 + Icges(5,3) * t400;
t336 = rSges(4,1) * t384 - rSges(4,2) * t383 + rSges(4,3) * t398;
t335 = rSges(4,1) * t382 - rSges(4,2) * t381 + rSges(4,3) * t396;
t334 = Icges(4,1) * t384 - Icges(4,4) * t383 + Icges(4,5) * t398;
t333 = Icges(4,1) * t382 - Icges(4,4) * t381 + Icges(4,5) * t396;
t332 = Icges(4,4) * t384 - Icges(4,2) * t383 + Icges(4,6) * t398;
t331 = Icges(4,4) * t382 - Icges(4,2) * t381 + Icges(4,6) * t396;
t330 = Icges(4,5) * t384 - Icges(4,6) * t383 + Icges(4,3) * t398;
t329 = Icges(4,5) * t382 - Icges(4,6) * t381 + Icges(4,3) * t396;
t328 = -pkin(4) * t431 + pkin(9) * t400 + t443 * t401;
t327 = rSges(6,1) * t377 + rSges(6,2) * t376 + rSges(6,3) * t400;
t326 = Icges(6,1) * t377 + Icges(6,4) * t376 + Icges(6,5) * t400;
t325 = Icges(6,4) * t377 + Icges(6,2) * t376 + Icges(6,6) * t400;
t324 = Icges(6,5) * t377 + Icges(6,6) * t376 + Icges(6,3) * t400;
t323 = qJD(5) * t381 + t348;
t322 = qJD(5) * t383 + t347;
t320 = qJD(1) + (t364 * t413 + t365 * t415) * t433;
t319 = rSges(5,1) * t356 + rSges(5,2) * t355 + rSges(5,3) * t383;
t318 = rSges(5,1) * t354 + rSges(5,2) * t353 + rSges(5,3) * t381;
t317 = Icges(5,1) * t356 + Icges(5,4) * t355 + Icges(5,5) * t383;
t316 = Icges(5,1) * t354 + Icges(5,4) * t353 + Icges(5,5) * t381;
t315 = Icges(5,4) * t356 + Icges(5,2) * t355 + Icges(5,6) * t383;
t314 = Icges(5,4) * t354 + Icges(5,2) * t353 + Icges(5,6) * t381;
t313 = Icges(5,5) * t356 + Icges(5,6) * t355 + Icges(5,3) * t383;
t312 = Icges(5,5) * t354 + Icges(5,6) * t353 + Icges(5,3) * t381;
t311 = rSges(6,1) * t352 + rSges(6,2) * t351 + rSges(6,3) * t383;
t310 = rSges(6,1) * t350 + rSges(6,2) * t349 + rSges(6,3) * t381;
t309 = Icges(6,1) * t352 + Icges(6,4) * t351 + Icges(6,5) * t383;
t308 = Icges(6,1) * t350 + Icges(6,4) * t349 + Icges(6,5) * t381;
t307 = Icges(6,4) * t352 + Icges(6,2) * t351 + Icges(6,6) * t383;
t306 = Icges(6,4) * t350 + Icges(6,2) * t349 + Icges(6,6) * t381;
t305 = Icges(6,5) * t352 + Icges(6,6) * t351 + Icges(6,3) * t383;
t304 = Icges(6,5) * t350 + Icges(6,6) * t349 + Icges(6,3) * t381;
t303 = pkin(4) * t440 + pkin(9) * t383 + t443 * t384;
t302 = pkin(4) * t441 + pkin(9) * t381 + t443 * t382;
t301 = -t335 * t403 + t371 * t388 + t425;
t300 = t336 * t403 - t371 * t387 + t427;
t299 = t335 * t387 - t336 * t388 + t428;
t298 = -t318 * t380 + t340 * t348 + t423;
t297 = t319 * t380 - t340 * t347 + t424;
t296 = t318 * t347 - t319 * t348 + t426;
t295 = -t302 * t380 - t310 * t357 + t323 * t327 + t328 * t348 + t423;
t294 = t303 * t380 + t311 * t357 - t322 * t327 - t328 * t347 + t424;
t293 = t302 * t347 - t303 * t348 + t310 * t322 - t311 * t323 + t426;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t320 ^ 2 + t341 ^ 2 + t342 ^ 2) / 0.2e1 - t445 * ((-t359 * t436 - t396 * t361 + t363 * t397) * t439 - (-t358 * t436 - t360 * t396 + t362 * t397) * t436 + (-t389 * t436 - t390 * t396 + t391 * t397) * t416) * t436 / 0.2e1 + m(4) * (t299 ^ 2 + t300 ^ 2 + t301 ^ 2) / 0.2e1 + t387 * ((t398 * t330 - t383 * t332 + t384 * t334) * t387 + (t329 * t398 - t331 * t383 + t333 * t384) * t388 + (t368 * t398 - t369 * t383 + t370 * t384) * t403) / 0.2e1 + t388 * ((t330 * t396 - t332 * t381 + t334 * t382) * t387 + (t396 * t329 - t381 * t331 + t382 * t333) * t388 + (t368 * t396 - t369 * t381 + t370 * t382) * t403) / 0.2e1 + t403 * ((-t330 * t437 - t332 * t400 + t334 * t401) * t387 + (-t329 * t437 - t331 * t400 + t333 * t401) * t388 + (-t368 * t437 - t400 * t369 + t401 * t370) * t403) / 0.2e1 + m(5) * (t296 ^ 2 + t297 ^ 2 + t298 ^ 2) / 0.2e1 + t347 * ((t383 * t313 + t355 * t315 + t356 * t317) * t347 + (t312 * t383 + t314 * t355 + t316 * t356) * t348 + (t337 * t383 + t338 * t355 + t339 * t356) * t380) / 0.2e1 + t348 * ((t313 * t381 + t315 * t353 + t317 * t354) * t347 + (t381 * t312 + t353 * t314 + t354 * t316) * t348 + (t337 * t381 + t338 * t353 + t339 * t354) * t380) / 0.2e1 + t380 * ((t313 * t400 + t315 * t385 + t317 * t386) * t347 + (t312 * t400 + t314 * t385 + t316 * t386) * t348 + (t400 * t337 + t385 * t338 + t386 * t339) * t380) / 0.2e1 + m(6) * (t293 ^ 2 + t294 ^ 2 + t295 ^ 2) / 0.2e1 + t322 * ((t383 * t305 + t351 * t307 + t352 * t309) * t322 + (t304 * t383 + t306 * t351 + t308 * t352) * t323 + (t324 * t383 + t325 * t351 + t326 * t352) * t357) / 0.2e1 + t323 * ((t305 * t381 + t307 * t349 + t309 * t350) * t322 + (t381 * t304 + t349 * t306 + t350 * t308) * t323 + (t324 * t381 + t325 * t349 + t326 * t350) * t357) / 0.2e1 + t357 * ((t305 * t400 + t307 * t376 + t309 * t377) * t322 + (t304 * t400 + t306 * t376 + t308 * t377) * t323 + (t400 * t324 + t376 * t325 + t377 * t326) * t357) / 0.2e1 + (((t359 * t439 - t361 * t398 + t363 * t399) * t439 - (t358 * t439 - t398 * t360 + t362 * t399) * t436 + (t389 * t439 - t390 * t398 + t391 * t399) * t416) * t439 + t416 * (t416 ^ 2 * t389 + (((t361 * t421 + t363 * t419) * t413 - (t360 * t421 + t362 * t419) * t415) * t414 + (-t358 * t415 + t359 * t413 + t390 * t421 + t391 * t419) * t416) * t414)) * t445 / 0.2e1;
T = t1;
