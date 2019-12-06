% Calculate kinetic energy for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:24
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR4_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR4_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR4_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRPR4_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRPR4_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:21:26
% EndTime: 2019-12-05 16:21:29
% DurationCPUTime: 2.33s
% Computational Cost: add. (1193->231), mult. (1728->369), div. (0->0), fcn. (1764->10), ass. (0->116)
t427 = -Icges(5,3) - Icges(4,3);
t372 = cos(pkin(8));
t419 = t372 ^ 2;
t371 = sin(pkin(8));
t420 = t371 ^ 2;
t421 = t419 + t420;
t426 = qJD(2) * t421;
t425 = t371 * t372;
t370 = qJ(3) + pkin(9);
t365 = sin(t370);
t366 = cos(t370);
t377 = cos(qJ(2));
t410 = t371 * t377;
t335 = -t365 * t410 - t366 * t372;
t336 = -t365 * t372 + t366 * t410;
t376 = cos(qJ(3));
t374 = sin(qJ(3));
t406 = t374 * t377;
t345 = -t371 * t406 - t372 * t376;
t405 = t376 * t377;
t409 = t372 * t374;
t346 = t371 * t405 - t409;
t375 = sin(qJ(2));
t411 = t371 * t375;
t424 = Icges(4,5) * t346 + Icges(5,5) * t336 + Icges(4,6) * t345 + Icges(5,6) * t335 - t427 * t411;
t407 = t372 * t377;
t337 = -t365 * t407 + t366 * t371;
t338 = t365 * t371 + t366 * t407;
t347 = t371 * t376 - t372 * t406;
t412 = t371 * t374;
t348 = t372 * t405 + t412;
t408 = t372 * t375;
t423 = Icges(4,5) * t348 + Icges(5,5) * t338 + Icges(4,6) * t347 + Icges(5,6) * t337 - t427 * t408;
t422 = t427 * t377 + (Icges(4,5) * t376 + Icges(5,5) * t366 - Icges(4,6) * t374 - Icges(5,6) * t365) * t375;
t418 = qJD(2) ^ 2;
t414 = t376 * pkin(3);
t380 = qJ(4) * t375 + t377 * t414;
t309 = pkin(3) * t412 + t372 * t380;
t402 = pkin(4) * t366;
t379 = pkin(7) * t375 + t377 * t402;
t391 = pkin(4) * t365;
t404 = -t371 * t391 - t372 * t379 - t309;
t403 = -rSges(5,1) * t338 - rSges(5,2) * t337 - rSges(5,3) * t408 - t309;
t364 = qJD(2) * t371;
t399 = qJD(3) * t375;
t349 = t372 * t399 + t364;
t400 = qJD(2) * t372;
t398 = qJD(3) * t377;
t397 = qJD(4) * t375;
t396 = qJD(5) * t375;
t356 = pkin(2) * t375 - pkin(6) * t377;
t395 = t356 * t364;
t394 = t356 * t400;
t393 = qJD(1) + (pkin(2) * t377 + pkin(6) * t375) * t426;
t350 = t371 * t399 - t400;
t390 = t371 * t397 - t395;
t387 = Icges(3,5) * t377 - Icges(3,6) * t375;
t308 = -pkin(3) * t409 + t371 * t380;
t318 = -qJ(4) * t377 + t375 * t414;
t382 = t308 * t398 + t350 * t318 + t372 * t397 - t394;
t381 = -qJD(4) * t377 + t349 * t308 + t393;
t367 = qJ(5) + t370;
t362 = cos(t367);
t361 = sin(t367);
t355 = rSges(3,1) * t375 + rSges(3,2) * t377;
t353 = (-qJD(3) - qJD(5)) * t377;
t344 = -rSges(4,3) * t377 + (rSges(4,1) * t376 - rSges(4,2) * t374) * t375;
t343 = -Icges(4,5) * t377 + (Icges(4,1) * t376 - Icges(4,4) * t374) * t375;
t342 = -Icges(4,6) * t377 + (Icges(4,4) * t376 - Icges(4,2) * t374) * t375;
t330 = Icges(3,3) * t371 + t372 * t387;
t329 = -Icges(3,3) * t372 + t371 * t387;
t328 = t371 * t396 + t350;
t327 = t372 * t396 + t349;
t326 = t361 * t371 + t362 * t407;
t325 = -t361 * t407 + t362 * t371;
t324 = -t361 * t372 + t362 * t410;
t323 = -t361 * t410 - t362 * t372;
t322 = -rSges(5,3) * t377 + (rSges(5,1) * t366 - rSges(5,2) * t365) * t375;
t321 = -Icges(5,5) * t377 + (Icges(5,1) * t366 - Icges(5,4) * t365) * t375;
t320 = -Icges(5,6) * t377 + (Icges(5,4) * t366 - Icges(5,2) * t365) * t375;
t317 = -rSges(6,3) * t377 + (rSges(6,1) * t362 - rSges(6,2) * t361) * t375;
t316 = -Icges(6,5) * t377 + (Icges(6,1) * t362 - Icges(6,4) * t361) * t375;
t315 = -Icges(6,6) * t377 + (Icges(6,4) * t362 - Icges(6,2) * t361) * t375;
t314 = -Icges(6,3) * t377 + (Icges(6,5) * t362 - Icges(6,6) * t361) * t375;
t312 = -pkin(7) * t377 + t375 * t402;
t311 = rSges(4,1) * t348 + rSges(4,2) * t347 + rSges(4,3) * t408;
t310 = rSges(4,1) * t346 + rSges(4,2) * t345 + rSges(4,3) * t411;
t307 = Icges(4,1) * t348 + Icges(4,4) * t347 + Icges(4,5) * t408;
t306 = Icges(4,1) * t346 + Icges(4,4) * t345 + Icges(4,5) * t411;
t305 = Icges(4,4) * t348 + Icges(4,2) * t347 + Icges(4,6) * t408;
t304 = Icges(4,4) * t346 + Icges(4,2) * t345 + Icges(4,6) * t411;
t300 = qJD(1) + (rSges(3,1) * t377 - rSges(3,2) * t375) * t426;
t298 = rSges(5,1) * t336 + rSges(5,2) * t335 + rSges(5,3) * t411;
t297 = Icges(5,1) * t338 + Icges(5,4) * t337 + Icges(5,5) * t408;
t296 = Icges(5,1) * t336 + Icges(5,4) * t335 + Icges(5,5) * t411;
t295 = Icges(5,4) * t338 + Icges(5,2) * t337 + Icges(5,6) * t408;
t294 = Icges(5,4) * t336 + Icges(5,2) * t335 + Icges(5,6) * t411;
t291 = rSges(6,1) * t326 + rSges(6,2) * t325 + rSges(6,3) * t408;
t290 = rSges(6,1) * t324 + rSges(6,2) * t323 + rSges(6,3) * t411;
t288 = Icges(6,1) * t326 + Icges(6,4) * t325 + Icges(6,5) * t408;
t287 = Icges(6,1) * t324 + Icges(6,4) * t323 + Icges(6,5) * t411;
t286 = Icges(6,4) * t326 + Icges(6,2) * t325 + Icges(6,6) * t408;
t285 = Icges(6,4) * t324 + Icges(6,2) * t323 + Icges(6,6) * t411;
t284 = Icges(6,5) * t326 + Icges(6,6) * t325 + Icges(6,3) * t408;
t283 = Icges(6,5) * t324 + Icges(6,6) * t323 + Icges(6,3) * t411;
t281 = t371 * t379 - t372 * t391;
t280 = t310 * t398 + t344 * t350 - t394;
t279 = -t311 * t398 - t344 * t349 - t395;
t278 = t310 * t349 - t311 * t350 + t393;
t277 = t298 * t398 + t322 * t350 + t382;
t276 = (-t318 - t322) * t349 + t403 * t398 + t390;
t275 = t298 * t349 + t350 * t403 + t381;
t274 = t281 * t398 - t290 * t353 + t312 * t350 + t317 * t328 + t382;
t273 = t291 * t353 - t317 * t327 + (-t312 - t318) * t349 + t404 * t398 + t390;
t272 = t281 * t349 + t290 * t327 - t291 * t328 + t350 * t404 + t381;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t421 * t418 * t355 ^ 2 + t300 ^ 2) / 0.2e1 + t418 * t371 * (-t329 * t425 + t420 * t330) / 0.2e1 - t418 * t372 * (t419 * t329 - t330 * t425) / 0.2e1 + m(4) * (t278 ^ 2 + t279 ^ 2 + t280 ^ 2) / 0.2e1 + m(5) * (t275 ^ 2 + t276 ^ 2 + t277 ^ 2) / 0.2e1 + m(6) * (t272 ^ 2 + t273 ^ 2 + t274 ^ 2) / 0.2e1 + t327 * ((t284 * t408 + t325 * t286 + t326 * t288) * t327 + (t283 * t408 + t285 * t325 + t287 * t326) * t328 + (t314 * t408 + t315 * t325 + t316 * t326) * t353) / 0.2e1 + t328 * ((t284 * t411 + t286 * t323 + t288 * t324) * t327 + (t283 * t411 + t323 * t285 + t324 * t287) * t328 + (t314 * t411 + t315 * t323 + t316 * t324) * t353) / 0.2e1 + t353 * ((-t283 * t328 - t284 * t327 - t314 * t353) * t377 + ((-t286 * t361 + t288 * t362) * t327 + (-t285 * t361 + t287 * t362) * t328 + (-t315 * t361 + t316 * t362) * t353) * t375) / 0.2e1 + ((-t320 * t337 - t338 * t321 - t342 * t347 - t343 * t348 - t422 * t408) * t398 + (t294 * t337 + t296 * t338 + t304 * t347 + t306 * t348 + t424 * t408) * t350 + (t337 * t295 + t338 * t297 + t347 * t305 + t348 * t307 + t423 * t408) * t349) * t349 / 0.2e1 + ((-t320 * t335 - t321 * t336 - t342 * t345 - t343 * t346 - t422 * t411) * t398 + (t335 * t294 + t336 * t296 + t345 * t304 + t346 * t306 + t424 * t411) * t350 + (t295 * t335 + t336 * t297 + t305 * t345 + t307 * t346 + t423 * t411) * t349) * t350 / 0.2e1 - ((-t423 * t349 - t424 * t350 + t422 * t398) * t377 + ((t320 * t365 - t321 * t366 + t342 * t374 - t343 * t376) * t398 + (-t294 * t365 + t296 * t366 - t304 * t374 + t306 * t376) * t350 + (-t295 * t365 + t297 * t366 - t305 * t374 + t307 * t376) * t349) * t375) * t398 / 0.2e1;
T = t1;
