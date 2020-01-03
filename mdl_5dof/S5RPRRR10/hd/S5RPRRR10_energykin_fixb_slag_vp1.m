% Calculate kinetic energy for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR10_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_energykin_fixb_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR10_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRRR10_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRRR10_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:09:42
% EndTime: 2019-12-31 19:09:43
% DurationCPUTime: 1.30s
% Computational Cost: add. (1166->224), mult. (1302->369), div. (0->0), fcn. (1284->10), ass. (0->118)
t359 = cos(pkin(9));
t400 = pkin(2) * t359;
t363 = cos(qJ(4));
t399 = pkin(4) * t363;
t356 = pkin(9) + qJ(3);
t350 = sin(t356);
t397 = Icges(4,4) * t350;
t351 = cos(t356);
t396 = Icges(4,4) * t351;
t362 = sin(qJ(1));
t395 = t350 * t362;
t364 = cos(qJ(1));
t394 = t350 * t364;
t357 = qJ(4) + qJ(5);
t354 = sin(t357);
t393 = t354 * t362;
t392 = t354 * t364;
t355 = cos(t357);
t391 = t355 * t362;
t390 = t355 * t364;
t361 = sin(qJ(4));
t389 = t361 * t362;
t388 = t361 * t364;
t387 = t362 * t363;
t386 = t363 * t364;
t344 = pkin(1) * t362 - qJ(2) * t364;
t384 = pkin(6) * t364 - t362 * t400 - t344;
t379 = pkin(3) * t351 + pkin(7) * t350;
t327 = t379 * t362;
t328 = t379 * t364;
t352 = qJD(3) * t362;
t382 = qJD(3) * t364;
t383 = t327 * t352 + t328 * t382;
t381 = qJD(4) * t350;
t334 = t364 * t381 + t352;
t380 = qJD(5) * t350;
t335 = t362 * t381 - t382;
t341 = qJD(1) * (pkin(1) * t364 + qJ(2) * t362);
t378 = -qJD(2) * t364 + qJD(1) * (pkin(6) * t362 + t364 * t400) + t341;
t358 = sin(pkin(9));
t377 = rSges(3,1) * t359 - rSges(3,2) * t358;
t376 = rSges(4,1) * t351 - rSges(4,2) * t350;
t375 = Icges(4,1) * t351 - t397;
t374 = -Icges(4,2) * t350 + t396;
t373 = Icges(4,5) * t351 - Icges(4,6) * t350;
t314 = -Icges(4,6) * t364 + t362 * t374;
t316 = -Icges(4,5) * t364 + t362 * t375;
t372 = t314 * t350 - t316 * t351;
t315 = Icges(4,6) * t362 + t364 * t374;
t317 = Icges(4,5) * t362 + t364 * t375;
t371 = -t315 * t350 + t317 * t351;
t337 = Icges(4,2) * t351 + t397;
t338 = Icges(4,1) * t350 + t396;
t370 = -t337 * t350 + t338 * t351;
t369 = pkin(8) * t350 + t351 * t399;
t340 = pkin(3) * t350 - pkin(7) * t351;
t368 = qJD(1) * t328 - t340 * t352 + t378;
t353 = qJD(2) * t362;
t367 = t353 + (-t327 + t384) * qJD(1) - t340 * t382;
t347 = -qJD(4) * t351 + qJD(1);
t346 = rSges(2,1) * t364 - rSges(2,2) * t362;
t345 = rSges(2,1) * t362 + rSges(2,2) * t364;
t339 = rSges(4,1) * t350 + rSges(4,2) * t351;
t336 = Icges(4,5) * t350 + Icges(4,6) * t351;
t333 = qJD(1) + (-qJD(4) - qJD(5)) * t351;
t332 = t351 * t386 + t389;
t331 = -t351 * t388 + t387;
t330 = t351 * t387 - t388;
t329 = -t351 * t389 - t386;
t325 = t351 * t390 + t393;
t324 = -t351 * t392 + t391;
t323 = t351 * t391 - t392;
t322 = -t351 * t393 - t390;
t319 = rSges(4,3) * t362 + t364 * t376;
t318 = -rSges(4,3) * t364 + t362 * t376;
t313 = Icges(4,3) * t362 + t364 * t373;
t312 = -Icges(4,3) * t364 + t362 * t373;
t311 = t362 * t380 + t335;
t310 = t364 * t380 + t334;
t308 = -rSges(5,3) * t351 + (rSges(5,1) * t363 - rSges(5,2) * t361) * t350;
t307 = -Icges(5,5) * t351 + (Icges(5,1) * t363 - Icges(5,4) * t361) * t350;
t306 = -Icges(5,6) * t351 + (Icges(5,4) * t363 - Icges(5,2) * t361) * t350;
t305 = -Icges(5,3) * t351 + (Icges(5,5) * t363 - Icges(5,6) * t361) * t350;
t303 = -rSges(6,3) * t351 + (rSges(6,1) * t355 - rSges(6,2) * t354) * t350;
t302 = -Icges(6,5) * t351 + (Icges(6,1) * t355 - Icges(6,4) * t354) * t350;
t301 = -Icges(6,6) * t351 + (Icges(6,4) * t355 - Icges(6,2) * t354) * t350;
t300 = -Icges(6,3) * t351 + (Icges(6,5) * t355 - Icges(6,6) * t354) * t350;
t299 = -pkin(8) * t351 + t350 * t399;
t298 = qJD(1) * t362 * rSges(3,3) + t341 + (qJD(1) * t377 - qJD(2)) * t364;
t297 = t353 + (t364 * rSges(3,3) - t362 * t377 - t344) * qJD(1);
t296 = rSges(5,1) * t332 + rSges(5,2) * t331 + rSges(5,3) * t394;
t295 = rSges(5,1) * t330 + rSges(5,2) * t329 + rSges(5,3) * t395;
t294 = Icges(5,1) * t332 + Icges(5,4) * t331 + Icges(5,5) * t394;
t293 = Icges(5,1) * t330 + Icges(5,4) * t329 + Icges(5,5) * t395;
t292 = Icges(5,4) * t332 + Icges(5,2) * t331 + Icges(5,6) * t394;
t291 = Icges(5,4) * t330 + Icges(5,2) * t329 + Icges(5,6) * t395;
t290 = Icges(5,5) * t332 + Icges(5,6) * t331 + Icges(5,3) * t394;
t289 = Icges(5,5) * t330 + Icges(5,6) * t329 + Icges(5,3) * t395;
t288 = pkin(4) * t389 + t364 * t369;
t287 = -pkin(4) * t388 + t362 * t369;
t286 = rSges(6,1) * t325 + rSges(6,2) * t324 + rSges(6,3) * t394;
t285 = rSges(6,1) * t323 + rSges(6,2) * t322 + rSges(6,3) * t395;
t284 = Icges(6,1) * t325 + Icges(6,4) * t324 + Icges(6,5) * t394;
t283 = Icges(6,1) * t323 + Icges(6,4) * t322 + Icges(6,5) * t395;
t282 = Icges(6,4) * t325 + Icges(6,2) * t324 + Icges(6,6) * t394;
t281 = Icges(6,4) * t323 + Icges(6,2) * t322 + Icges(6,6) * t395;
t280 = Icges(6,5) * t325 + Icges(6,6) * t324 + Icges(6,3) * t394;
t279 = Icges(6,5) * t323 + Icges(6,6) * t322 + Icges(6,3) * t395;
t278 = (t318 * t362 + t319 * t364) * qJD(3);
t277 = qJD(1) * t319 - t339 * t352 + t378;
t276 = -t339 * t382 + t353 + (-t318 + t384) * qJD(1);
t275 = t295 * t334 - t296 * t335 + t383;
t274 = t296 * t347 - t308 * t334 + t368;
t273 = -t295 * t347 + t308 * t335 + t367;
t272 = t286 * t333 + t288 * t347 - t299 * t334 - t303 * t310 + t368;
t271 = -t285 * t333 - t287 * t347 + t299 * t335 + t303 * t311 + t367;
t270 = t285 * t310 - t286 * t311 + t287 * t334 - t288 * t335 + t383;
t1 = m(3) * (t297 ^ 2 + t298 ^ 2) / 0.2e1 + m(4) * (t276 ^ 2 + t277 ^ 2 + t278 ^ 2) / 0.2e1 + ((t362 * t336 + t364 * t370) * qJD(1) + (t362 ^ 2 * t313 + (t372 * t364 + (-t312 + t371) * t362) * t364) * qJD(3)) * t352 / 0.2e1 - ((-t364 * t336 + t362 * t370) * qJD(1) + (t364 ^ 2 * t312 + (t371 * t362 + (-t313 + t372) * t364) * t362) * qJD(3)) * t382 / 0.2e1 + qJD(1) * ((t351 * t337 + t350 * t338) * qJD(1) + ((t315 * t351 + t317 * t350) * t362 - (t314 * t351 + t316 * t350) * t364) * qJD(3)) / 0.2e1 + m(5) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + t334 * ((t290 * t394 + t331 * t292 + t332 * t294) * t334 + (t289 * t394 + t291 * t331 + t293 * t332) * t335 + (t305 * t394 + t306 * t331 + t307 * t332) * t347) / 0.2e1 + t335 * ((t290 * t395 + t292 * t329 + t294 * t330) * t334 + (t289 * t395 + t329 * t291 + t330 * t293) * t335 + (t305 * t395 + t306 * t329 + t307 * t330) * t347) / 0.2e1 + t347 * ((-t289 * t335 - t290 * t334 - t305 * t347) * t351 + ((-t292 * t361 + t294 * t363) * t334 + (-t291 * t361 + t293 * t363) * t335 + (-t306 * t361 + t307 * t363) * t347) * t350) / 0.2e1 + m(6) * (t270 ^ 2 + t271 ^ 2 + t272 ^ 2) / 0.2e1 + t310 * ((t280 * t394 + t324 * t282 + t325 * t284) * t310 + (t279 * t394 + t281 * t324 + t283 * t325) * t311 + (t300 * t394 + t301 * t324 + t302 * t325) * t333) / 0.2e1 + t311 * ((t280 * t395 + t282 * t322 + t284 * t323) * t310 + (t279 * t395 + t322 * t281 + t323 * t283) * t311 + (t300 * t395 + t301 * t322 + t302 * t323) * t333) / 0.2e1 + t333 * ((-t279 * t311 - t280 * t310 - t300 * t333) * t351 + ((-t282 * t354 + t284 * t355) * t310 + (-t281 * t354 + t283 * t355) * t311 + (-t301 * t354 + t302 * t355) * t333) * t350) / 0.2e1 + (m(2) * (t345 ^ 2 + t346 ^ 2) + Icges(2,3) + Icges(3,2) * t359 ^ 2 + (Icges(3,1) * t358 + 0.2e1 * Icges(3,4) * t359) * t358) * qJD(1) ^ 2 / 0.2e1;
T = t1;
