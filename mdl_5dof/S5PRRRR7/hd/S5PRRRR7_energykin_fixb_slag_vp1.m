% Calculate kinetic energy for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR7_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR7_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR7_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:11:45
% EndTime: 2019-12-05 17:11:47
% DurationCPUTime: 2.25s
% Computational Cost: add. (1229->238), mult. (1788->402), div. (0->0), fcn. (1824->10), ass. (0->122)
t372 = cos(pkin(9));
t415 = t372 ^ 2;
t371 = sin(pkin(9));
t416 = t371 ^ 2;
t417 = t415 + t416;
t419 = qJD(2) * t417;
t418 = t371 * t372;
t414 = qJD(2) ^ 2;
t375 = cos(qJ(3));
t412 = t375 * pkin(3);
t373 = sin(qJ(3));
t410 = t371 * t373;
t374 = sin(qJ(2));
t409 = t371 * t374;
t376 = cos(qJ(2));
t408 = t371 * t376;
t407 = t372 * t373;
t406 = t372 * t374;
t405 = t372 * t376;
t404 = t373 * t376;
t403 = t375 * t376;
t370 = qJ(3) + qJ(4);
t366 = cos(t370);
t402 = pkin(4) * t366;
t364 = qJD(2) * t371;
t399 = qJD(3) * t374;
t349 = t372 * t399 + t364;
t400 = qJD(2) * t372;
t398 = qJD(3) * t376;
t397 = qJD(4) * t374;
t396 = qJD(5) * t374;
t395 = -qJD(3) - qJD(4);
t326 = t372 * t397 + t349;
t356 = pkin(2) * t374 - pkin(6) * t376;
t394 = t356 * t364;
t393 = t356 * t400;
t392 = qJD(1) + (pkin(2) * t376 + pkin(6) * t374) * t419;
t365 = sin(t370);
t391 = pkin(4) * t365;
t350 = t371 * t399 - t400;
t327 = t371 * t397 + t350;
t388 = Icges(3,5) * t376 - Icges(3,6) * t374;
t381 = pkin(7) * t374 + t376 * t412;
t307 = -pkin(3) * t407 + t371 * t381;
t317 = -pkin(7) * t376 + t374 * t412;
t385 = t307 * t398 + t350 * t317 - t393;
t308 = pkin(3) * t410 + t372 * t381;
t382 = t349 * t307 - t308 * t350 + t392;
t380 = pkin(8) * t374 + t376 * t402;
t379 = -t308 * t398 - t317 * t349 - t394;
t367 = qJ(5) + t370;
t362 = cos(t367);
t361 = sin(t367);
t355 = rSges(3,1) * t374 + rSges(3,2) * t376;
t353 = t395 * t376;
t348 = (-qJD(5) + t395) * t376;
t347 = t372 * t403 + t410;
t346 = t371 * t375 - t372 * t404;
t345 = t371 * t403 - t407;
t344 = -t371 * t404 - t372 * t375;
t343 = -rSges(4,3) * t376 + (rSges(4,1) * t375 - rSges(4,2) * t373) * t374;
t342 = -Icges(4,5) * t376 + (Icges(4,1) * t375 - Icges(4,4) * t373) * t374;
t341 = -Icges(4,6) * t376 + (Icges(4,4) * t375 - Icges(4,2) * t373) * t374;
t340 = -Icges(4,3) * t376 + (Icges(4,5) * t375 - Icges(4,6) * t373) * t374;
t339 = t365 * t371 + t366 * t405;
t338 = -t365 * t405 + t366 * t371;
t337 = -t365 * t372 + t366 * t408;
t336 = -t365 * t408 - t366 * t372;
t329 = Icges(3,3) * t371 + t372 * t388;
t328 = -Icges(3,3) * t372 + t371 * t388;
t325 = t361 * t371 + t362 * t405;
t324 = -t361 * t405 + t362 * t371;
t323 = -t361 * t372 + t362 * t408;
t322 = -t361 * t408 - t362 * t372;
t321 = -rSges(5,3) * t376 + (rSges(5,1) * t366 - rSges(5,2) * t365) * t374;
t320 = -Icges(5,5) * t376 + (Icges(5,1) * t366 - Icges(5,4) * t365) * t374;
t319 = -Icges(5,6) * t376 + (Icges(5,4) * t366 - Icges(5,2) * t365) * t374;
t318 = -Icges(5,3) * t376 + (Icges(5,5) * t366 - Icges(5,6) * t365) * t374;
t316 = -rSges(6,3) * t376 + (rSges(6,1) * t362 - rSges(6,2) * t361) * t374;
t315 = -Icges(6,5) * t376 + (Icges(6,1) * t362 - Icges(6,4) * t361) * t374;
t314 = -Icges(6,6) * t376 + (Icges(6,4) * t362 - Icges(6,2) * t361) * t374;
t313 = -Icges(6,3) * t376 + (Icges(6,5) * t362 - Icges(6,6) * t361) * t374;
t312 = t371 * t396 + t327;
t311 = t372 * t396 + t326;
t310 = -pkin(8) * t376 + t374 * t402;
t306 = rSges(4,1) * t347 + rSges(4,2) * t346 + rSges(4,3) * t406;
t305 = rSges(4,1) * t345 + rSges(4,2) * t344 + rSges(4,3) * t409;
t304 = Icges(4,1) * t347 + Icges(4,4) * t346 + Icges(4,5) * t406;
t303 = Icges(4,1) * t345 + Icges(4,4) * t344 + Icges(4,5) * t409;
t302 = Icges(4,4) * t347 + Icges(4,2) * t346 + Icges(4,6) * t406;
t301 = Icges(4,4) * t345 + Icges(4,2) * t344 + Icges(4,6) * t409;
t300 = Icges(4,5) * t347 + Icges(4,6) * t346 + Icges(4,3) * t406;
t299 = Icges(4,5) * t345 + Icges(4,6) * t344 + Icges(4,3) * t409;
t297 = rSges(5,1) * t339 + rSges(5,2) * t338 + rSges(5,3) * t406;
t296 = rSges(5,1) * t337 + rSges(5,2) * t336 + rSges(5,3) * t409;
t295 = Icges(5,1) * t339 + Icges(5,4) * t338 + Icges(5,5) * t406;
t294 = Icges(5,1) * t337 + Icges(5,4) * t336 + Icges(5,5) * t409;
t293 = Icges(5,4) * t339 + Icges(5,2) * t338 + Icges(5,6) * t406;
t292 = Icges(5,4) * t337 + Icges(5,2) * t336 + Icges(5,6) * t409;
t291 = Icges(5,5) * t339 + Icges(5,6) * t338 + Icges(5,3) * t406;
t290 = Icges(5,5) * t337 + Icges(5,6) * t336 + Icges(5,3) * t409;
t289 = qJD(1) + (rSges(3,1) * t376 - rSges(3,2) * t374) * t419;
t288 = rSges(6,1) * t325 + rSges(6,2) * t324 + rSges(6,3) * t406;
t287 = rSges(6,1) * t323 + rSges(6,2) * t322 + rSges(6,3) * t409;
t286 = Icges(6,1) * t325 + Icges(6,4) * t324 + Icges(6,5) * t406;
t285 = Icges(6,1) * t323 + Icges(6,4) * t322 + Icges(6,5) * t409;
t284 = Icges(6,4) * t325 + Icges(6,2) * t324 + Icges(6,6) * t406;
t283 = Icges(6,4) * t323 + Icges(6,2) * t322 + Icges(6,6) * t409;
t282 = Icges(6,5) * t325 + Icges(6,6) * t324 + Icges(6,3) * t406;
t281 = Icges(6,5) * t323 + Icges(6,6) * t322 + Icges(6,3) * t409;
t279 = t371 * t391 + t372 * t380;
t278 = t371 * t380 - t372 * t391;
t277 = t305 * t398 + t343 * t350 - t393;
t276 = -t306 * t398 - t343 * t349 - t394;
t275 = t305 * t349 - t306 * t350 + t392;
t274 = -t296 * t353 + t321 * t327 + t385;
t273 = t297 * t353 - t321 * t326 + t379;
t272 = t296 * t326 - t297 * t327 + t382;
t271 = -t278 * t353 - t287 * t348 + t310 * t327 + t312 * t316 + t385;
t270 = t279 * t353 + t288 * t348 - t310 * t326 - t311 * t316 + t379;
t269 = t278 * t326 - t279 * t327 + t287 * t311 - t288 * t312 + t382;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t417 * t414 * t355 ^ 2 + t289 ^ 2) / 0.2e1 + t414 * t371 * (-t328 * t418 + t416 * t329) / 0.2e1 - t414 * t372 * (t415 * t328 - t329 * t418) / 0.2e1 + m(4) * (t275 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + t349 * ((t300 * t406 + t302 * t346 + t304 * t347) * t349 + (t299 * t406 + t301 * t346 + t303 * t347) * t350 - (t340 * t406 + t341 * t346 + t342 * t347) * t398) / 0.2e1 + t350 * ((t300 * t409 + t302 * t344 + t304 * t345) * t349 + (t299 * t409 + t301 * t344 + t303 * t345) * t350 - (t340 * t409 + t341 * t344 + t342 * t345) * t398) / 0.2e1 - ((-t299 * t350 - t300 * t349 + t340 * t398) * t376 + ((-t302 * t373 + t304 * t375) * t349 + (-t301 * t373 + t303 * t375) * t350 - (-t341 * t373 + t342 * t375) * t398) * t374) * t398 / 0.2e1 + m(5) * (t272 ^ 2 + t273 ^ 2 + t274 ^ 2) / 0.2e1 + t326 * ((t291 * t406 + t293 * t338 + t295 * t339) * t326 + (t290 * t406 + t292 * t338 + t294 * t339) * t327 + (t318 * t406 + t319 * t338 + t320 * t339) * t353) / 0.2e1 + t327 * ((t291 * t409 + t293 * t336 + t295 * t337) * t326 + (t290 * t409 + t292 * t336 + t294 * t337) * t327 + (t318 * t409 + t319 * t336 + t320 * t337) * t353) / 0.2e1 + t353 * ((-t290 * t327 - t291 * t326 - t318 * t353) * t376 + ((-t293 * t365 + t295 * t366) * t326 + (-t292 * t365 + t294 * t366) * t327 + (-t319 * t365 + t320 * t366) * t353) * t374) / 0.2e1 + m(6) * (t269 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + t311 * ((t282 * t406 + t284 * t324 + t286 * t325) * t311 + (t281 * t406 + t283 * t324 + t285 * t325) * t312 + (t313 * t406 + t314 * t324 + t315 * t325) * t348) / 0.2e1 + t312 * ((t282 * t409 + t284 * t322 + t286 * t323) * t311 + (t281 * t409 + t283 * t322 + t285 * t323) * t312 + (t313 * t409 + t314 * t322 + t315 * t323) * t348) / 0.2e1 + t348 * ((-t281 * t312 - t282 * t311 - t313 * t348) * t376 + ((-t284 * t361 + t286 * t362) * t311 + (-t283 * t361 + t285 * t362) * t312 + (-t314 * t361 + t315 * t362) * t348) * t374) / 0.2e1;
T = t1;
