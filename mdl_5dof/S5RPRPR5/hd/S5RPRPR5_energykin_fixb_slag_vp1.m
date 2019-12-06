% Calculate kinetic energy for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:55:50
% EndTime: 2019-12-05 17:55:52
% DurationCPUTime: 1.87s
% Computational Cost: add. (1145->232), mult. (1594->367), div. (0->0), fcn. (1632->10), ass. (0->117)
t428 = -Icges(5,3) - Icges(4,3);
t373 = sin(pkin(8));
t374 = cos(pkin(8));
t372 = qJ(3) + pkin(9);
t366 = cos(t372);
t379 = cos(qJ(1));
t407 = t379 * t366;
t365 = sin(t372);
t377 = sin(qJ(1));
t415 = t377 * t365;
t341 = t374 * t415 + t407;
t408 = t379 * t365;
t414 = t377 * t366;
t342 = -t374 * t414 + t408;
t343 = -t374 * t408 + t414;
t344 = t374 * t407 + t415;
t378 = cos(qJ(3));
t404 = t379 * t378;
t376 = sin(qJ(3));
t412 = t377 * t376;
t348 = t374 * t412 + t404;
t405 = t379 * t376;
t411 = t377 * t378;
t349 = -t374 * t411 + t405;
t350 = -t374 * t405 + t411;
t351 = t374 * t404 + t412;
t406 = t379 * t373;
t413 = t377 * t373;
t425 = (-Icges(4,5) * t351 - Icges(5,5) * t344 - Icges(4,6) * t350 - Icges(5,6) * t343 + t428 * t406) * t379 + (Icges(4,5) * t349 + Icges(5,5) * t342 + Icges(4,6) * t348 + Icges(5,6) * t341 + t428 * t413) * t377;
t427 = t425 * t373;
t426 = -t428 * t374 + (-Icges(4,5) * t378 - Icges(5,5) * t366 + Icges(4,6) * t376 + Icges(5,6) * t365) * t373;
t401 = pkin(4) * t366;
t424 = pkin(7) * t373 + t374 * t401;
t419 = t378 * pkin(3);
t423 = qJ(4) * t373 + t374 * t419;
t367 = qJ(5) + t372;
t362 = sin(t367);
t417 = t377 * t362;
t363 = cos(t367);
t416 = t377 * t363;
t410 = t379 * t362;
t409 = t379 * t363;
t316 = pkin(3) * t412 + t423 * t379;
t390 = pkin(4) * t365;
t403 = -t390 * t377 - t424 * t379 - t316;
t402 = -t344 * rSges(5,1) - t343 * rSges(5,2) - rSges(5,3) * t406 - t316;
t400 = qJD(1) * (-t377 * pkin(1) + t379 * qJ(2)) + qJD(2) * t377;
t398 = qJD(3) * t373;
t397 = qJD(3) * t379;
t396 = qJD(4) * t374;
t395 = qJD(4) * t377;
t394 = qJD(3) + qJD(5);
t387 = pkin(2) * t374 + pkin(6) * t373;
t393 = -qJD(1) * t387 * t377 + t400;
t392 = t377 * t398;
t391 = t373 * t397;
t386 = -rSges(3,1) * t374 + rSges(3,2) * t373;
t358 = t379 * pkin(1) + t377 * qJ(2);
t369 = qJD(2) * t379;
t383 = t369 + (-t387 * t379 - t358) * qJD(1);
t315 = pkin(3) * t405 - t423 * t377;
t328 = -qJ(4) * t374 + t373 * t419;
t361 = -qJD(3) * t374 + qJD(1);
t382 = qJD(4) * t406 + t361 * t315 + t328 * t392 + t393;
t381 = t328 * t391 + t383;
t359 = t379 * rSges(2,1) - t377 * rSges(2,2);
t357 = -t377 * rSges(2,1) - t379 * rSges(2,2);
t354 = -t374 * t394 + qJD(1);
t347 = t394 * t406;
t346 = t394 * t413;
t340 = -t374 * rSges(4,3) + (rSges(4,1) * t378 - rSges(4,2) * t376) * t373;
t339 = -Icges(4,5) * t374 + (Icges(4,1) * t378 - Icges(4,4) * t376) * t373;
t338 = -Icges(4,6) * t374 + (Icges(4,4) * t378 - Icges(4,2) * t376) * t373;
t336 = t374 * t409 + t417;
t335 = -t374 * t410 + t416;
t334 = -t374 * t416 + t410;
t333 = t374 * t417 + t409;
t332 = -t374 * rSges(5,3) + (rSges(5,1) * t366 - rSges(5,2) * t365) * t373;
t331 = -Icges(5,5) * t374 + (Icges(5,1) * t366 - Icges(5,4) * t365) * t373;
t330 = -Icges(5,6) * t374 + (Icges(5,4) * t366 - Icges(5,2) * t365) * t373;
t327 = -t374 * rSges(6,3) + (rSges(6,1) * t363 - rSges(6,2) * t362) * t373;
t326 = -Icges(6,5) * t374 + (Icges(6,1) * t363 - Icges(6,4) * t362) * t373;
t325 = -Icges(6,6) * t374 + (Icges(6,4) * t363 - Icges(6,2) * t362) * t373;
t324 = -Icges(6,3) * t374 + (Icges(6,5) * t363 - Icges(6,6) * t362) * t373;
t321 = t369 + (-t377 * rSges(3,3) + t379 * t386 - t358) * qJD(1);
t320 = qJD(1) * (t379 * rSges(3,3) + t377 * t386) + t400;
t319 = -pkin(7) * t374 + t373 * t401;
t318 = t351 * rSges(4,1) + t350 * rSges(4,2) + rSges(4,3) * t406;
t317 = t349 * rSges(4,1) + t348 * rSges(4,2) - rSges(4,3) * t413;
t314 = Icges(4,1) * t351 + Icges(4,4) * t350 + Icges(4,5) * t406;
t313 = Icges(4,1) * t349 + Icges(4,4) * t348 - Icges(4,5) * t413;
t312 = Icges(4,4) * t351 + Icges(4,2) * t350 + Icges(4,6) * t406;
t311 = Icges(4,4) * t349 + Icges(4,2) * t348 - Icges(4,6) * t413;
t306 = t342 * rSges(5,1) + t341 * rSges(5,2) - rSges(5,3) * t413;
t305 = Icges(5,1) * t344 + Icges(5,4) * t343 + Icges(5,5) * t406;
t304 = Icges(5,1) * t342 + Icges(5,4) * t341 - Icges(5,5) * t413;
t303 = Icges(5,4) * t344 + Icges(5,2) * t343 + Icges(5,6) * t406;
t302 = Icges(5,4) * t342 + Icges(5,2) * t341 - Icges(5,6) * t413;
t299 = t336 * rSges(6,1) + t335 * rSges(6,2) + rSges(6,3) * t406;
t298 = t334 * rSges(6,1) + t333 * rSges(6,2) - rSges(6,3) * t413;
t297 = Icges(6,1) * t336 + Icges(6,4) * t335 + Icges(6,5) * t406;
t296 = Icges(6,1) * t334 + Icges(6,4) * t333 - Icges(6,5) * t413;
t295 = Icges(6,4) * t336 + Icges(6,2) * t335 + Icges(6,6) * t406;
t294 = Icges(6,4) * t334 + Icges(6,2) * t333 - Icges(6,6) * t413;
t293 = Icges(6,5) * t336 + Icges(6,6) * t335 + Icges(6,3) * t406;
t292 = Icges(6,5) * t334 + Icges(6,6) * t333 - Icges(6,3) * t413;
t290 = -t424 * t377 + t390 * t379;
t289 = (-t317 * t379 - t318 * t377) * t398;
t288 = -t361 * t318 + t340 * t391 + t383;
t287 = t361 * t317 + t340 * t392 + t393;
t286 = (t332 * t397 - t395) * t373 + t402 * t361 + t381;
t285 = t361 * t306 + t332 * t392 + t382;
t284 = -t396 + ((-t306 - t315) * t379 + t402 * t377) * t398;
t283 = -t354 * t299 + t347 * t327 + (t319 * t397 - t395) * t373 + t403 * t361 + t381;
t282 = t361 * t290 + t354 * t298 + t319 * t392 + t346 * t327 + t382;
t281 = -t396 - t347 * t298 - t346 * t299 + ((-t290 - t315) * t379 + t403 * t377) * t398;
t1 = m(3) * (t320 ^ 2 + t321 ^ 2) / 0.2e1 + m(4) * (t287 ^ 2 + t288 ^ 2 + t289 ^ 2) / 0.2e1 + m(5) * (t284 ^ 2 + t285 ^ 2 + t286 ^ 2) / 0.2e1 + m(6) * (t281 ^ 2 + t282 ^ 2 + t283 ^ 2) / 0.2e1 + t354 * ((t292 * t346 - t293 * t347 - t324 * t354) * t374 + ((-t325 * t362 + t326 * t363) * t354 - (-t294 * t362 + t296 * t363) * t346 + (-t295 * t362 + t297 * t363) * t347) * t373) / 0.2e1 - t346 * ((-t324 * t413 + t333 * t325 + t334 * t326) * t354 - (-t292 * t413 + t333 * t294 + t334 * t296) * t346 + (-t293 * t413 + t333 * t295 + t334 * t297) * t347) / 0.2e1 + t347 * ((t324 * t406 + t335 * t325 + t336 * t326) * t354 - (t292 * t406 + t335 * t294 + t336 * t296) * t346 + (t293 * t406 + t335 * t295 + t336 * t297) * t347) / 0.2e1 + ((t425 * t374 + ((-t303 * t365 + t305 * t366 - t312 * t376 + t314 * t378) * t379 + (t302 * t365 - t304 * t366 + t311 * t376 - t313 * t378) * t377) * t373) * t398 + (t426 * t374 + (-t330 * t365 + t331 * t366 - t338 * t376 + t339 * t378) * t373) * t361) * t361 / 0.2e1 - (((t341 * t303 + t342 * t305 + t348 * t312 + t349 * t314) * t379 + (-t341 * t302 - t342 * t304 - t348 * t311 - t349 * t313 + t427) * t377) * t398 + (t341 * t330 + t342 * t331 + t348 * t338 + t349 * t339 + t426 * t413) * t361) * t392 / 0.2e1 + (((t343 * t303 + t344 * t305 + t350 * t312 + t351 * t314 - t427) * t379 + (-t343 * t302 - t344 * t304 - t350 * t311 - t351 * t313) * t377) * t398 + (t343 * t330 + t344 * t331 + t350 * t338 + t351 * t339 - t426 * t406) * t361) * t391 / 0.2e1 + (m(2) * (t357 ^ 2 + t359 ^ 2) + Icges(2,3) + Icges(3,2) * t374 ^ 2 + (Icges(3,1) * t373 + 0.2e1 * Icges(3,4) * t374) * t373) * qJD(1) ^ 2 / 0.2e1;
T = t1;
