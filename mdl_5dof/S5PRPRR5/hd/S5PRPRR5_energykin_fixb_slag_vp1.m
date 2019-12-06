% Calculate kinetic energy for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPRR5_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR5_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRPRR5_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRPRR5_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:53:24
% EndTime: 2019-12-05 15:53:25
% DurationCPUTime: 1.66s
% Computational Cost: add. (1077->213), mult. (1464->356), div. (0->0), fcn. (1493->10), ass. (0->108)
t356 = cos(pkin(8));
t396 = t356 ^ 2;
t354 = sin(pkin(8));
t397 = t354 ^ 2;
t398 = t396 + t397;
t399 = t398 * qJD(2);
t395 = qJD(2) ^ 2;
t355 = cos(pkin(9));
t393 = t355 * pkin(3);
t353 = sin(pkin(9));
t392 = t353 * t354;
t391 = t353 * t356;
t358 = sin(qJ(2));
t390 = t354 * t358;
t359 = cos(qJ(2));
t389 = t354 * t359;
t388 = t356 * t358;
t387 = t356 * t359;
t352 = pkin(9) + qJ(4);
t348 = cos(t352);
t385 = pkin(4) * t348;
t346 = qJD(2) * t354;
t381 = qJD(4) * t358;
t331 = t356 * t381 + t346;
t383 = qJD(2) * t356;
t382 = qJD(3) * t358;
t380 = qJD(4) * t359;
t379 = qJD(5) * t358;
t347 = sin(t352);
t376 = pkin(4) * t347;
t337 = pkin(2) * t358 - qJ(3) * t359;
t375 = qJD(2) * (pkin(6) * t359 - t393 * t358 - t337);
t374 = qJD(2) * (rSges(4,3) * t359 - (rSges(4,1) * t355 - rSges(4,2) * t353) * t358 - t337);
t332 = t354 * t381 - t383;
t371 = Icges(3,5) * t359 - Icges(3,6) * t358;
t368 = -qJD(3) * t359 + qJD(1) + (pkin(2) * t359 + qJ(3) * t358) * t399;
t341 = t354 * t382;
t367 = t354 * t375 + t341;
t342 = t356 * t382;
t366 = t356 * t375 + t342;
t362 = pkin(6) * t358 + t393 * t359;
t363 = (-pkin(3) * t391 + t362 * t354) * t346 + (pkin(3) * t392 + t362 * t356) * t383 + t368;
t361 = pkin(7) * t358 + t385 * t359;
t349 = qJ(5) + t352;
t344 = cos(t349);
t343 = sin(t349);
t338 = rSges(3,1) * t358 + rSges(3,2) * t359;
t335 = (-qJD(4) - qJD(5)) * t359;
t330 = t355 * t387 + t392;
t329 = -t353 * t387 + t354 * t355;
t328 = t355 * t389 - t391;
t327 = -t353 * t389 - t355 * t356;
t323 = t347 * t354 + t348 * t387;
t322 = -t347 * t387 + t348 * t354;
t321 = -t347 * t356 + t348 * t389;
t320 = -t347 * t389 - t348 * t356;
t315 = Icges(3,3) * t354 + t371 * t356;
t314 = -Icges(3,3) * t356 + t371 * t354;
t313 = t354 * t379 + t332;
t312 = t356 * t379 + t331;
t311 = t343 * t354 + t344 * t387;
t310 = -t343 * t387 + t344 * t354;
t309 = -t343 * t356 + t344 * t389;
t308 = -t343 * t389 - t344 * t356;
t307 = -rSges(5,3) * t359 + (rSges(5,1) * t348 - rSges(5,2) * t347) * t358;
t306 = -Icges(5,5) * t359 + (Icges(5,1) * t348 - Icges(5,4) * t347) * t358;
t305 = -Icges(5,6) * t359 + (Icges(5,4) * t348 - Icges(5,2) * t347) * t358;
t304 = -Icges(5,3) * t359 + (Icges(5,5) * t348 - Icges(5,6) * t347) * t358;
t302 = -rSges(6,3) * t359 + (rSges(6,1) * t344 - rSges(6,2) * t343) * t358;
t301 = -Icges(6,5) * t359 + (Icges(6,1) * t344 - Icges(6,4) * t343) * t358;
t300 = -Icges(6,6) * t359 + (Icges(6,4) * t344 - Icges(6,2) * t343) * t358;
t299 = -Icges(6,3) * t359 + (Icges(6,5) * t344 - Icges(6,6) * t343) * t358;
t298 = -pkin(7) * t359 + t385 * t358;
t297 = Icges(4,1) * t330 + Icges(4,4) * t329 + Icges(4,5) * t388;
t296 = Icges(4,1) * t328 + Icges(4,4) * t327 + Icges(4,5) * t390;
t295 = Icges(4,4) * t330 + Icges(4,2) * t329 + Icges(4,6) * t388;
t294 = Icges(4,4) * t328 + Icges(4,2) * t327 + Icges(4,6) * t390;
t293 = Icges(4,5) * t330 + Icges(4,6) * t329 + Icges(4,3) * t388;
t292 = Icges(4,5) * t328 + Icges(4,6) * t327 + Icges(4,3) * t390;
t291 = t356 * t374 + t342;
t290 = t354 * t374 + t341;
t287 = qJD(1) + (rSges(3,1) * t359 - rSges(3,2) * t358) * t399;
t286 = rSges(5,1) * t323 + rSges(5,2) * t322 + rSges(5,3) * t388;
t285 = rSges(5,1) * t321 + rSges(5,2) * t320 + rSges(5,3) * t390;
t284 = Icges(5,1) * t323 + Icges(5,4) * t322 + Icges(5,5) * t388;
t283 = Icges(5,1) * t321 + Icges(5,4) * t320 + Icges(5,5) * t390;
t282 = Icges(5,4) * t323 + Icges(5,2) * t322 + Icges(5,6) * t388;
t281 = Icges(5,4) * t321 + Icges(5,2) * t320 + Icges(5,6) * t390;
t280 = Icges(5,5) * t323 + Icges(5,6) * t322 + Icges(5,3) * t388;
t279 = Icges(5,5) * t321 + Icges(5,6) * t320 + Icges(5,3) * t390;
t278 = rSges(6,1) * t311 + rSges(6,2) * t310 + rSges(6,3) * t388;
t277 = rSges(6,1) * t309 + rSges(6,2) * t308 + rSges(6,3) * t390;
t276 = Icges(6,1) * t311 + Icges(6,4) * t310 + Icges(6,5) * t388;
t275 = Icges(6,1) * t309 + Icges(6,4) * t308 + Icges(6,5) * t390;
t274 = Icges(6,4) * t311 + Icges(6,2) * t310 + Icges(6,6) * t388;
t273 = Icges(6,4) * t309 + Icges(6,2) * t308 + Icges(6,6) * t390;
t272 = Icges(6,5) * t311 + Icges(6,6) * t310 + Icges(6,3) * t388;
t271 = Icges(6,5) * t309 + Icges(6,6) * t308 + Icges(6,3) * t390;
t270 = t376 * t354 + t361 * t356;
t269 = t361 * t354 - t376 * t356;
t268 = (t354 * (rSges(4,1) * t328 + rSges(4,2) * t327 + rSges(4,3) * t390) + t356 * (rSges(4,1) * t330 + rSges(4,2) * t329 + rSges(4,3) * t388)) * qJD(2) + t368;
t267 = t285 * t380 + t307 * t332 + t366;
t266 = -t286 * t380 - t307 * t331 + t367;
t265 = t285 * t331 - t286 * t332 + t363;
t264 = t269 * t380 - t277 * t335 + t298 * t332 + t302 * t313 + t366;
t263 = -t270 * t380 + t278 * t335 - t298 * t331 - t302 * t312 + t367;
t262 = t269 * t331 - t270 * t332 + t277 * t312 - t278 * t313 + t363;
t1 = m(2) * qJD(1) ^ 2 / 0.2e1 + m(3) * (t398 * t395 * t338 ^ 2 + t287 ^ 2) / 0.2e1 + m(4) * (t268 ^ 2 + t290 ^ 2 + t291 ^ 2) / 0.2e1 + m(5) * (t265 ^ 2 + t266 ^ 2 + t267 ^ 2) / 0.2e1 + t331 * ((t280 * t388 + t282 * t322 + t284 * t323) * t331 + (t279 * t388 + t281 * t322 + t283 * t323) * t332 - (t304 * t388 + t305 * t322 + t306 * t323) * t380) / 0.2e1 + t332 * ((t280 * t390 + t282 * t320 + t284 * t321) * t331 + (t279 * t390 + t281 * t320 + t283 * t321) * t332 - (t304 * t390 + t305 * t320 + t306 * t321) * t380) / 0.2e1 - ((-t279 * t332 - t280 * t331 + t304 * t380) * t359 + ((-t282 * t347 + t284 * t348) * t331 + (-t281 * t347 + t283 * t348) * t332 - (-t305 * t347 + t306 * t348) * t380) * t358) * t380 / 0.2e1 + m(6) * (t262 ^ 2 + t263 ^ 2 + t264 ^ 2) / 0.2e1 + t312 * ((t272 * t388 + t274 * t310 + t276 * t311) * t312 + (t271 * t388 + t273 * t310 + t275 * t311) * t313 + (t299 * t388 + t300 * t310 + t301 * t311) * t335) / 0.2e1 + t313 * ((t272 * t390 + t274 * t308 + t276 * t309) * t312 + (t271 * t390 + t273 * t308 + t275 * t309) * t313 + (t299 * t390 + t300 * t308 + t301 * t309) * t335) / 0.2e1 + t335 * ((-t271 * t313 - t272 * t312 - t299 * t335) * t359 + ((-t274 * t343 + t276 * t344) * t312 + (-t273 * t343 + t275 * t344) * t313 + (-t300 * t343 + t301 * t344) * t335) * t358) / 0.2e1 + (t397 * t315 + (t293 * t388 + t295 * t329 + t330 * t297) * t354 + (-t292 * t388 - t294 * t329 - t296 * t330 - t354 * t314) * t356) * t354 * t395 / 0.2e1 - (t396 * t314 - (t292 * t390 + t294 * t327 + t296 * t328) * t356 + (t293 * t390 + t295 * t327 + t297 * t328 - t356 * t315) * t354) * t356 * t395 / 0.2e1;
T = t1;
