% Calculate kinetic energy for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR7_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR7_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR7_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR7_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRPR7_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:16:05
% EndTime: 2019-12-31 21:16:06
% DurationCPUTime: 1.75s
% Computational Cost: add. (1246->260), mult. (1511->415), div. (0->0), fcn. (1464->10), ass. (0->137)
t367 = cos(qJ(2));
t413 = pkin(2) * t367;
t363 = cos(pkin(9));
t412 = pkin(4) * t363;
t365 = sin(qJ(2));
t410 = Icges(3,4) * t365;
t409 = Icges(3,4) * t367;
t361 = qJ(2) + qJ(3);
t358 = sin(t361);
t408 = Icges(4,4) * t358;
t359 = cos(t361);
t407 = Icges(4,4) * t359;
t366 = sin(qJ(1));
t406 = t358 * t366;
t368 = cos(qJ(1));
t405 = t358 * t368;
t404 = t359 * t366;
t403 = t359 * t368;
t362 = sin(pkin(9));
t402 = t362 * t366;
t401 = t362 * t368;
t400 = t363 * t366;
t399 = t363 * t368;
t303 = -pkin(7) * t368 + t366 * t413;
t304 = pkin(7) * t366 + t368 * t413;
t357 = qJD(2) * t366;
t395 = qJD(2) * t368;
t397 = t303 * t357 + t304 * t395;
t350 = pkin(1) * t366 - pkin(6) * t368;
t396 = -t303 - t350;
t340 = qJD(3) * t366 + t357;
t394 = qJD(4) * t358;
t393 = qJD(5) * t358;
t392 = pkin(2) * qJD(2) * t365;
t387 = pkin(3) * t359 + qJ(4) * t358;
t328 = t387 * t366;
t391 = -t328 + t396;
t341 = (-qJD(2) - qJD(3)) * t368;
t390 = t368 * t392;
t389 = rSges(3,1) * t367 - rSges(3,2) * t365;
t388 = rSges(4,1) * t359 - rSges(4,2) * t358;
t386 = Icges(3,1) * t367 - t410;
t385 = Icges(4,1) * t359 - t408;
t384 = -Icges(3,2) * t365 + t409;
t383 = -Icges(4,2) * t358 + t407;
t382 = Icges(3,5) * t367 - Icges(3,6) * t365;
t381 = Icges(4,5) * t359 - Icges(4,6) * t358;
t320 = -Icges(3,6) * t368 + t366 * t384;
t322 = -Icges(3,5) * t368 + t366 * t386;
t380 = t320 * t365 - t322 * t367;
t321 = Icges(3,6) * t366 + t368 * t384;
t323 = Icges(3,5) * t366 + t368 * t386;
t379 = -t321 * t365 + t323 * t367;
t343 = Icges(3,2) * t367 + t410;
t344 = Icges(3,1) * t365 + t409;
t378 = -t343 * t365 + t344 * t367;
t377 = -qJD(4) * t359 + t340 * t328 + t397;
t339 = qJD(1) * (pkin(1) * t368 + pkin(6) * t366);
t376 = qJD(1) * t304 - t366 * t392 + t339;
t337 = pkin(3) * t358 - qJ(4) * t359;
t375 = t341 * t337 + t368 * t394 - t390;
t374 = (Icges(4,5) * t358 + Icges(4,6) * t359) * qJD(1) + (-Icges(4,3) * t368 + t366 * t381) * t341 + (Icges(4,3) * t366 + t368 * t381) * t340;
t373 = pkin(8) * t358 + t359 * t412;
t329 = t387 * t368;
t372 = qJD(1) * t329 + t366 * t394 + t376;
t307 = -Icges(4,6) * t368 + t366 * t383;
t308 = Icges(4,6) * t366 + t368 * t383;
t309 = -Icges(4,5) * t368 + t366 * t385;
t310 = Icges(4,5) * t366 + t368 * t385;
t335 = Icges(4,2) * t359 + t408;
t336 = Icges(4,1) * t358 + t407;
t371 = (-t308 * t358 + t310 * t359) * t340 + (-t307 * t358 + t309 * t359) * t341 + (-t335 * t358 + t336 * t359) * qJD(1);
t360 = pkin(9) + qJ(5);
t355 = cos(t360);
t354 = sin(t360);
t351 = -qJD(5) * t359 + qJD(1);
t347 = rSges(2,1) * t368 - rSges(2,2) * t366;
t346 = rSges(2,1) * t366 + rSges(2,2) * t368;
t345 = rSges(3,1) * t365 + rSges(3,2) * t367;
t342 = Icges(3,5) * t365 + Icges(3,6) * t367;
t338 = rSges(4,1) * t358 + rSges(4,2) * t359;
t333 = t359 * t399 + t402;
t332 = -t359 * t401 + t400;
t331 = t359 * t400 - t401;
t330 = -t359 * t402 - t399;
t327 = rSges(3,3) * t366 + t368 * t389;
t326 = -rSges(3,3) * t368 + t366 * t389;
t325 = t366 * t393 + t341;
t324 = t368 * t393 + t340;
t319 = Icges(3,3) * t366 + t368 * t382;
t318 = -Icges(3,3) * t368 + t366 * t382;
t316 = t354 * t366 + t355 * t403;
t315 = -t354 * t403 + t355 * t366;
t314 = -t354 * t368 + t355 * t404;
t313 = -t354 * t404 - t355 * t368;
t312 = rSges(4,3) * t366 + t368 * t388;
t311 = -rSges(4,3) * t368 + t366 * t388;
t301 = -rSges(5,3) * t359 + (rSges(5,1) * t363 - rSges(5,2) * t362) * t358;
t299 = -Icges(5,5) * t359 + (Icges(5,1) * t363 - Icges(5,4) * t362) * t358;
t298 = -Icges(5,6) * t359 + (Icges(5,4) * t363 - Icges(5,2) * t362) * t358;
t297 = -Icges(5,3) * t359 + (Icges(5,5) * t363 - Icges(5,6) * t362) * t358;
t293 = -rSges(6,3) * t359 + (rSges(6,1) * t355 - rSges(6,2) * t354) * t358;
t292 = -Icges(6,5) * t359 + (Icges(6,1) * t355 - Icges(6,4) * t354) * t358;
t291 = -Icges(6,6) * t359 + (Icges(6,4) * t355 - Icges(6,2) * t354) * t358;
t290 = -Icges(6,3) * t359 + (Icges(6,5) * t355 - Icges(6,6) * t354) * t358;
t289 = -pkin(8) * t359 + t358 * t412;
t288 = qJD(1) * t327 - t345 * t357 + t339;
t287 = -t345 * t395 + (-t326 - t350) * qJD(1);
t286 = (t326 * t366 + t327 * t368) * qJD(2);
t285 = rSges(5,1) * t333 + rSges(5,2) * t332 + rSges(5,3) * t405;
t284 = rSges(5,1) * t331 + rSges(5,2) * t330 + rSges(5,3) * t406;
t283 = Icges(5,1) * t333 + Icges(5,4) * t332 + Icges(5,5) * t405;
t282 = Icges(5,1) * t331 + Icges(5,4) * t330 + Icges(5,5) * t406;
t281 = Icges(5,4) * t333 + Icges(5,2) * t332 + Icges(5,6) * t405;
t280 = Icges(5,4) * t331 + Icges(5,2) * t330 + Icges(5,6) * t406;
t279 = Icges(5,5) * t333 + Icges(5,6) * t332 + Icges(5,3) * t405;
t278 = Icges(5,5) * t331 + Icges(5,6) * t330 + Icges(5,3) * t406;
t277 = pkin(4) * t402 + t368 * t373;
t276 = -pkin(4) * t401 + t366 * t373;
t275 = rSges(6,1) * t316 + rSges(6,2) * t315 + rSges(6,3) * t405;
t274 = rSges(6,1) * t314 + rSges(6,2) * t313 + rSges(6,3) * t406;
t273 = Icges(6,1) * t316 + Icges(6,4) * t315 + Icges(6,5) * t405;
t272 = Icges(6,1) * t314 + Icges(6,4) * t313 + Icges(6,5) * t406;
t271 = Icges(6,4) * t316 + Icges(6,2) * t315 + Icges(6,6) * t405;
t270 = Icges(6,4) * t314 + Icges(6,2) * t313 + Icges(6,6) * t406;
t269 = Icges(6,5) * t316 + Icges(6,6) * t315 + Icges(6,3) * t405;
t268 = Icges(6,5) * t314 + Icges(6,6) * t313 + Icges(6,3) * t406;
t267 = qJD(1) * t312 - t338 * t340 + t376;
t266 = -t390 + t338 * t341 + (-t311 + t396) * qJD(1);
t265 = t311 * t340 - t312 * t341 + t397;
t264 = qJD(1) * t285 + (-t301 - t337) * t340 + t372;
t263 = t301 * t341 + (-t284 + t391) * qJD(1) + t375;
t262 = t284 * t340 + (-t285 - t329) * t341 + t377;
t261 = qJD(1) * t277 + t275 * t351 - t293 * t324 + (-t289 - t337) * t340 + t372;
t260 = -t274 * t351 + t289 * t341 + t293 * t325 + (-t276 + t391) * qJD(1) + t375;
t259 = t274 * t324 - t275 * t325 + t276 * t340 + (-t277 - t329) * t341 + t377;
t1 = m(3) * (t286 ^ 2 + t287 ^ 2 + t288 ^ 2) / 0.2e1 + ((t366 * t342 + t368 * t378) * qJD(1) + (t366 ^ 2 * t319 + (t380 * t368 + (-t318 + t379) * t366) * t368) * qJD(2)) * t357 / 0.2e1 - ((-t368 * t342 + t366 * t378) * qJD(1) + (t368 ^ 2 * t318 + (t379 * t366 + (-t319 + t380) * t368) * t366) * qJD(2)) * t395 / 0.2e1 + m(4) * (t265 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + m(5) * (t262 ^ 2 + t263 ^ 2 + t264 ^ 2) / 0.2e1 + m(6) * (t259 ^ 2 + t260 ^ 2 + t261 ^ 2) / 0.2e1 + t324 * ((t269 * t405 + t271 * t315 + t273 * t316) * t324 + (t268 * t405 + t270 * t315 + t272 * t316) * t325 + (t290 * t405 + t291 * t315 + t292 * t316) * t351) / 0.2e1 + t325 * ((t269 * t406 + t271 * t313 + t273 * t314) * t324 + (t268 * t406 + t270 * t313 + t272 * t314) * t325 + (t290 * t406 + t291 * t313 + t292 * t314) * t351) / 0.2e1 + t351 * ((-t268 * t325 - t269 * t324 - t290 * t351) * t359 + ((-t271 * t354 + t273 * t355) * t324 + (-t270 * t354 + t272 * t355) * t325 + (-t291 * t354 + t292 * t355) * t351) * t358) / 0.2e1 + (t374 * t366 + t371 * t368 + (t279 * t405 + t281 * t332 + t283 * t333) * t340 + (t278 * t405 + t280 * t332 + t282 * t333) * t341 + (t297 * t405 + t298 * t332 + t299 * t333) * qJD(1)) * t340 / 0.2e1 + (t371 * t366 - t374 * t368 + (t279 * t406 + t281 * t330 + t283 * t331) * t340 + (t278 * t406 + t280 * t330 + t282 * t331) * t341 + (t297 * t406 + t298 * t330 + t299 * t331) * qJD(1)) * t341 / 0.2e1 + (m(2) * (t346 ^ 2 + t347 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + (((t321 * t367 + t323 * t365) * t366 - (t320 * t367 + t322 * t365) * t368) * qJD(2) + (t308 * t359 + t310 * t358) * t340 + (t307 * t359 + t309 * t358) * t341 + (-t278 * t341 - t279 * t340) * t359 + ((-t281 * t362 + t283 * t363) * t340 + (-t280 * t362 + t282 * t363) * t341) * t358 + (t367 * t343 + t365 * t344 + (t335 - t297) * t359 + (-t298 * t362 + t299 * t363 + t336) * t358) * qJD(1)) * qJD(1) / 0.2e1;
T = t1;
