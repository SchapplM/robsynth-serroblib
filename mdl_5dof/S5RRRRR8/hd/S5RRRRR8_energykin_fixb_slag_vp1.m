% Calculate kinetic energy for
% S5RRRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR8_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR8_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR8_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:24:21
% EndTime: 2019-12-31 22:24:23
% DurationCPUTime: 1.74s
% Computational Cost: add. (1318->258), mult. (1571->426), div. (0->0), fcn. (1524->10), ass. (0->139)
t370 = cos(qJ(2));
t416 = pkin(2) * t370;
t369 = cos(qJ(4));
t415 = pkin(4) * t369;
t367 = sin(qJ(2));
t412 = Icges(3,4) * t367;
t411 = Icges(3,4) * t370;
t365 = qJ(2) + qJ(3);
t361 = sin(t365);
t410 = Icges(4,4) * t361;
t363 = cos(t365);
t409 = Icges(4,4) * t363;
t368 = sin(qJ(1));
t408 = t361 * t368;
t371 = cos(qJ(1));
t407 = t361 * t371;
t406 = t363 * t368;
t405 = t363 * t371;
t366 = sin(qJ(4));
t404 = t366 * t368;
t403 = t366 * t371;
t402 = t368 * t369;
t401 = t369 * t371;
t305 = -pkin(7) * t371 + t368 * t416;
t306 = pkin(7) * t368 + t371 * t416;
t359 = qJD(2) * t368;
t398 = qJD(2) * t371;
t400 = t305 * t359 + t306 * t398;
t354 = pkin(1) * t368 - pkin(6) * t371;
t399 = -t305 - t354;
t344 = qJD(3) * t368 + t359;
t397 = qJD(4) * t361;
t396 = qJD(5) * t361;
t395 = pkin(2) * qJD(2) * t367;
t327 = t371 * t397 + t344;
t345 = (-qJD(2) - qJD(3)) * t371;
t394 = t371 * t395;
t393 = pkin(3) * t363 + pkin(8) * t361;
t392 = rSges(3,1) * t370 - rSges(3,2) * t367;
t391 = rSges(4,1) * t363 - rSges(4,2) * t361;
t390 = Icges(3,1) * t370 - t412;
t389 = Icges(4,1) * t363 - t410;
t388 = -Icges(3,2) * t367 + t411;
t387 = -Icges(4,2) * t361 + t409;
t386 = Icges(3,5) * t370 - Icges(3,6) * t367;
t385 = Icges(4,5) * t363 - Icges(4,6) * t361;
t323 = -Icges(3,6) * t371 + t368 * t388;
t325 = -Icges(3,5) * t371 + t368 * t390;
t384 = t323 * t367 - t325 * t370;
t324 = Icges(3,6) * t368 + t371 * t388;
t326 = Icges(3,5) * t368 + t371 * t390;
t383 = -t324 * t367 + t326 * t370;
t347 = Icges(3,2) * t370 + t412;
t348 = Icges(3,1) * t367 + t411;
t382 = -t347 * t367 + t348 * t370;
t328 = t368 * t397 + t345;
t331 = t393 * t368;
t332 = t393 * t371;
t381 = t344 * t331 - t332 * t345 + t400;
t343 = qJD(1) * (pkin(1) * t371 + pkin(6) * t368);
t380 = qJD(1) * t306 - t368 * t395 + t343;
t379 = pkin(9) * t361 + t363 * t415;
t378 = (Icges(4,5) * t361 + Icges(4,6) * t363) * qJD(1) + (-Icges(4,3) * t371 + t368 * t385) * t345 + (Icges(4,3) * t368 + t371 * t385) * t344;
t342 = pkin(3) * t361 - pkin(8) * t363;
t377 = qJD(1) * t332 - t342 * t344 + t380;
t376 = t345 * t342 + (-t331 + t399) * qJD(1) - t394;
t310 = -Icges(4,6) * t371 + t368 * t387;
t311 = Icges(4,6) * t368 + t371 * t387;
t312 = -Icges(4,5) * t371 + t368 * t389;
t313 = Icges(4,5) * t368 + t371 * t389;
t339 = Icges(4,2) * t363 + t410;
t340 = Icges(4,1) * t361 + t409;
t375 = (-t311 * t361 + t313 * t363) * t344 + (-t310 * t361 + t312 * t363) * t345 + (-t339 * t361 + t340 * t363) * qJD(1);
t364 = qJ(4) + qJ(5);
t362 = cos(t364);
t360 = sin(t364);
t355 = -qJD(4) * t363 + qJD(1);
t351 = rSges(2,1) * t371 - rSges(2,2) * t368;
t350 = rSges(2,1) * t368 + rSges(2,2) * t371;
t349 = rSges(3,1) * t367 + rSges(3,2) * t370;
t346 = Icges(3,5) * t367 + Icges(3,6) * t370;
t341 = rSges(4,1) * t361 + rSges(4,2) * t363;
t337 = qJD(1) + (-qJD(4) - qJD(5)) * t363;
t336 = t363 * t401 + t404;
t335 = -t363 * t403 + t402;
t334 = t363 * t402 - t403;
t333 = -t363 * t404 - t401;
t330 = rSges(3,3) * t368 + t371 * t392;
t329 = -rSges(3,3) * t371 + t368 * t392;
t322 = Icges(3,3) * t368 + t371 * t386;
t321 = -Icges(3,3) * t371 + t368 * t386;
t319 = t360 * t368 + t362 * t405;
t318 = -t360 * t405 + t362 * t368;
t317 = -t360 * t371 + t362 * t406;
t316 = -t360 * t406 - t362 * t371;
t315 = rSges(4,3) * t368 + t371 * t391;
t314 = -rSges(4,3) * t371 + t368 * t391;
t304 = -rSges(5,3) * t363 + (rSges(5,1) * t369 - rSges(5,2) * t366) * t361;
t303 = -Icges(5,5) * t363 + (Icges(5,1) * t369 - Icges(5,4) * t366) * t361;
t302 = -Icges(5,6) * t363 + (Icges(5,4) * t369 - Icges(5,2) * t366) * t361;
t301 = -Icges(5,3) * t363 + (Icges(5,5) * t369 - Icges(5,6) * t366) * t361;
t297 = -rSges(6,3) * t363 + (rSges(6,1) * t362 - rSges(6,2) * t360) * t361;
t296 = t368 * t396 + t328;
t295 = t371 * t396 + t327;
t293 = -Icges(6,5) * t363 + (Icges(6,1) * t362 - Icges(6,4) * t360) * t361;
t292 = -Icges(6,6) * t363 + (Icges(6,4) * t362 - Icges(6,2) * t360) * t361;
t291 = -Icges(6,3) * t363 + (Icges(6,5) * t362 - Icges(6,6) * t360) * t361;
t290 = -pkin(9) * t363 + t361 * t415;
t289 = qJD(1) * t330 - t349 * t359 + t343;
t288 = -t349 * t398 + (-t329 - t354) * qJD(1);
t287 = rSges(5,1) * t336 + rSges(5,2) * t335 + rSges(5,3) * t407;
t286 = rSges(5,1) * t334 + rSges(5,2) * t333 + rSges(5,3) * t408;
t285 = Icges(5,1) * t336 + Icges(5,4) * t335 + Icges(5,5) * t407;
t284 = Icges(5,1) * t334 + Icges(5,4) * t333 + Icges(5,5) * t408;
t283 = Icges(5,4) * t336 + Icges(5,2) * t335 + Icges(5,6) * t407;
t282 = Icges(5,4) * t334 + Icges(5,2) * t333 + Icges(5,6) * t408;
t281 = Icges(5,5) * t336 + Icges(5,6) * t335 + Icges(5,3) * t407;
t280 = Icges(5,5) * t334 + Icges(5,6) * t333 + Icges(5,3) * t408;
t279 = (t329 * t368 + t330 * t371) * qJD(2);
t278 = pkin(4) * t404 + t371 * t379;
t277 = -pkin(4) * t403 + t368 * t379;
t276 = rSges(6,1) * t319 + rSges(6,2) * t318 + rSges(6,3) * t407;
t275 = rSges(6,1) * t317 + rSges(6,2) * t316 + rSges(6,3) * t408;
t274 = Icges(6,1) * t319 + Icges(6,4) * t318 + Icges(6,5) * t407;
t273 = Icges(6,1) * t317 + Icges(6,4) * t316 + Icges(6,5) * t408;
t272 = Icges(6,4) * t319 + Icges(6,2) * t318 + Icges(6,6) * t407;
t271 = Icges(6,4) * t317 + Icges(6,2) * t316 + Icges(6,6) * t408;
t270 = Icges(6,5) * t319 + Icges(6,6) * t318 + Icges(6,3) * t407;
t269 = Icges(6,5) * t317 + Icges(6,6) * t316 + Icges(6,3) * t408;
t268 = qJD(1) * t315 - t341 * t344 + t380;
t267 = -t394 + t341 * t345 + (-t314 + t399) * qJD(1);
t266 = t314 * t344 - t315 * t345 + t400;
t265 = t287 * t355 - t304 * t327 + t377;
t264 = -t286 * t355 + t304 * t328 + t376;
t263 = t286 * t327 - t287 * t328 + t381;
t262 = t276 * t337 + t278 * t355 - t290 * t327 - t295 * t297 + t377;
t261 = -t275 * t337 - t277 * t355 + t290 * t328 + t296 * t297 + t376;
t260 = t275 * t295 - t276 * t296 + t277 * t327 - t278 * t328 + t381;
t1 = m(3) * (t279 ^ 2 + t288 ^ 2 + t289 ^ 2) / 0.2e1 + ((t368 * t346 + t371 * t382) * qJD(1) + (t368 ^ 2 * t322 + (t384 * t371 + (-t321 + t383) * t368) * t371) * qJD(2)) * t359 / 0.2e1 - ((-t371 * t346 + t368 * t382) * qJD(1) + (t371 ^ 2 * t321 + (t383 * t368 + (-t322 + t384) * t371) * t368) * qJD(2)) * t398 / 0.2e1 + m(4) * (t266 ^ 2 + t267 ^ 2 + t268 ^ 2) / 0.2e1 + t344 * (t368 * t378 + t371 * t375) / 0.2e1 + t345 * (t368 * t375 - t378 * t371) / 0.2e1 + m(5) * (t263 ^ 2 + t264 ^ 2 + t265 ^ 2) / 0.2e1 + t327 * ((t281 * t407 + t335 * t283 + t336 * t285) * t327 + (t280 * t407 + t282 * t335 + t284 * t336) * t328 + (t301 * t407 + t302 * t335 + t303 * t336) * t355) / 0.2e1 + t328 * ((t281 * t408 + t283 * t333 + t285 * t334) * t327 + (t280 * t408 + t333 * t282 + t334 * t284) * t328 + (t301 * t408 + t302 * t333 + t303 * t334) * t355) / 0.2e1 + t355 * ((-t280 * t328 - t281 * t327 - t301 * t355) * t363 + ((-t283 * t366 + t285 * t369) * t327 + (-t282 * t366 + t284 * t369) * t328 + (-t302 * t366 + t303 * t369) * t355) * t361) / 0.2e1 + m(6) * (t260 ^ 2 + t261 ^ 2 + t262 ^ 2) / 0.2e1 + t295 * ((t270 * t407 + t318 * t272 + t319 * t274) * t295 + (t269 * t407 + t271 * t318 + t273 * t319) * t296 + (t291 * t407 + t292 * t318 + t293 * t319) * t337) / 0.2e1 + t296 * ((t270 * t408 + t272 * t316 + t274 * t317) * t295 + (t269 * t408 + t316 * t271 + t317 * t273) * t296 + (t291 * t408 + t292 * t316 + t293 * t317) * t337) / 0.2e1 + t337 * ((-t269 * t296 - t270 * t295 - t291 * t337) * t363 + ((-t272 * t360 + t274 * t362) * t295 + (-t271 * t360 + t273 * t362) * t296 + (-t292 * t360 + t293 * t362) * t337) * t361) / 0.2e1 + (m(2) * (t350 ^ 2 + t351 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t324 * t370 + t326 * t367) * t368 - (t323 * t370 + t325 * t367) * t371) * qJD(2) + (t311 * t363 + t313 * t361) * t344 + (t310 * t363 + t312 * t361) * t345 + (t339 * t363 + t340 * t361 + t370 * t347 + t367 * t348) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
