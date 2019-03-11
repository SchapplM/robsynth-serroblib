% Calculate kinetic energy for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
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
% Datum: 2019-03-08 18:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR1_energykin_fixb_slag_vp1(qJ, qJD, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_energykin_fixb_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_energykin_fixb_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_energykin_fixb_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_energykin_fixb_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR1_energykin_fixb_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR1_energykin_fixb_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:36:54
% EndTime: 2019-03-08 18:36:55
% DurationCPUTime: 1.30s
% Computational Cost: add. (1157->227), mult. (1273->378), div. (0->0), fcn. (1170->10), ass. (0->137)
t339 = sin(qJ(1));
t392 = pkin(1) * t339;
t336 = qJ(2) + qJ(3);
t332 = sin(t336);
t391 = pkin(3) * t332;
t341 = cos(qJ(2));
t389 = t341 * pkin(2);
t338 = sin(qJ(2));
t388 = Icges(3,4) * t338;
t387 = Icges(3,4) * t341;
t386 = Icges(4,4) * t332;
t333 = cos(t336);
t385 = Icges(4,4) * t333;
t334 = qJ(4) + t336;
t325 = sin(t334);
t384 = Icges(5,4) * t325;
t326 = cos(t334);
t383 = Icges(5,4) * t326;
t382 = t325 * t339;
t342 = cos(qJ(1));
t381 = t325 * t342;
t337 = sin(qJ(5));
t380 = t337 * t339;
t379 = t337 * t342;
t340 = cos(qJ(5));
t378 = t339 * t340;
t377 = t340 * t342;
t376 = pkin(3) * t333;
t331 = qJD(2) * t342;
t317 = qJD(3) * t342 + t331;
t375 = qJD(2) * t339;
t374 = qJD(5) * t325;
t373 = -qJD(2) - qJD(3);
t372 = pkin(2) * qJD(2) * t338;
t308 = qJD(4) * t342 + t317;
t309 = t389 * t339;
t371 = -t309 - t392;
t370 = t342 * t372;
t276 = t376 * t339;
t369 = -t276 + t371;
t368 = pkin(4) * t326 + pkin(6) * t325;
t367 = rSges(3,1) * t341 - rSges(3,2) * t338;
t366 = rSges(4,1) * t333 - rSges(4,2) * t332;
t365 = rSges(5,1) * t326 - rSges(5,2) * t325;
t307 = (-qJD(4) + t373) * t339;
t364 = Icges(3,1) * t341 - t388;
t363 = Icges(4,1) * t333 - t386;
t362 = Icges(5,1) * t326 - t384;
t361 = -Icges(3,2) * t338 + t387;
t360 = -Icges(4,2) * t332 + t385;
t359 = -Icges(5,2) * t325 + t383;
t358 = Icges(3,5) * t341 - Icges(3,6) * t338;
t357 = Icges(4,5) * t333 - Icges(4,6) * t332;
t356 = Icges(5,5) * t326 - Icges(5,6) * t325;
t290 = Icges(3,6) * t342 + t361 * t339;
t292 = Icges(3,5) * t342 + t364 * t339;
t355 = -t290 * t338 + t292 * t341;
t291 = -Icges(3,6) * t339 + t361 * t342;
t293 = -Icges(3,5) * t339 + t364 * t342;
t354 = t291 * t338 - t293 * t341;
t320 = -Icges(3,2) * t341 - t388;
t321 = -Icges(3,1) * t338 - t387;
t353 = -t320 * t338 + t321 * t341;
t310 = t389 * t342;
t328 = qJD(1) * t342 * pkin(1);
t352 = qJD(1) * t310 - t339 * t372 + t328;
t351 = (-t309 * t339 - t310 * t342) * qJD(2);
t350 = -t317 * t391 - t370;
t349 = (-Icges(5,5) * t325 - Icges(5,6) * t326) * qJD(1) + (Icges(5,3) * t342 + t356 * t339) * t308 + (-Icges(5,3) * t339 + t356 * t342) * t307;
t316 = t373 * t339;
t348 = (-Icges(4,5) * t332 - Icges(4,6) * t333) * qJD(1) + (Icges(4,3) * t342 + t357 * t339) * t317 + (-Icges(4,3) * t339 + t357 * t342) * t316;
t277 = t376 * t342;
t347 = qJD(1) * t277 + t316 * t391 + t352;
t346 = t316 * t276 - t277 * t317 + t351;
t268 = Icges(5,6) * t342 + t359 * t339;
t269 = -Icges(5,6) * t339 + t359 * t342;
t270 = Icges(5,5) * t342 + t362 * t339;
t271 = -Icges(5,5) * t339 + t362 * t342;
t302 = -Icges(5,2) * t326 - t384;
t303 = -Icges(5,1) * t325 - t383;
t345 = (-t269 * t325 + t271 * t326) * t307 + (-t268 * t325 + t270 * t326) * t308 + (-t302 * t325 + t303 * t326) * qJD(1);
t280 = Icges(4,6) * t342 + t360 * t339;
t281 = -Icges(4,6) * t339 + t360 * t342;
t282 = Icges(4,5) * t342 + t363 * t339;
t283 = -Icges(4,5) * t339 + t363 * t342;
t312 = -Icges(4,2) * t333 - t386;
t313 = -Icges(4,1) * t332 - t385;
t344 = (-t281 * t332 + t283 * t333) * t316 + (-t280 * t332 + t282 * t333) * t317 + (-t312 * t332 + t313 * t333) * qJD(1);
t324 = rSges(2,1) * t342 - rSges(2,2) * t339;
t323 = rSges(2,1) * t339 + rSges(2,2) * t342;
t322 = -rSges(3,1) * t338 - rSges(3,2) * t341;
t319 = -Icges(3,5) * t338 - Icges(3,6) * t341;
t318 = qJD(5) * t326 + qJD(1);
t314 = -rSges(4,1) * t332 - rSges(4,2) * t333;
t306 = -pkin(4) * t325 + pkin(6) * t326;
t305 = -rSges(5,1) * t325 - rSges(5,2) * t326;
t299 = t326 * t377 - t380;
t298 = -t326 * t379 - t378;
t297 = t326 * t378 + t379;
t296 = -t326 * t380 + t377;
t295 = -rSges(3,3) * t339 + t367 * t342;
t294 = rSges(3,3) * t342 + t367 * t339;
t289 = -Icges(3,3) * t339 + t358 * t342;
t288 = Icges(3,3) * t342 + t358 * t339;
t287 = t368 * t342;
t286 = t368 * t339;
t285 = -rSges(4,3) * t339 + t366 * t342;
t284 = rSges(4,3) * t342 + t366 * t339;
t275 = -rSges(5,3) * t339 + t365 * t342;
t274 = rSges(5,3) * t342 + t365 * t339;
t273 = t339 * t374 + t308;
t272 = t342 * t374 + t307;
t264 = rSges(6,3) * t326 + (-rSges(6,1) * t340 + rSges(6,2) * t337) * t325;
t263 = Icges(6,5) * t326 + (-Icges(6,1) * t340 + Icges(6,4) * t337) * t325;
t262 = Icges(6,6) * t326 + (-Icges(6,4) * t340 + Icges(6,2) * t337) * t325;
t261 = Icges(6,3) * t326 + (-Icges(6,5) * t340 + Icges(6,6) * t337) * t325;
t259 = qJD(1) * t295 + t322 * t375 + t328;
t258 = t322 * t331 + (-t294 - t392) * qJD(1);
t257 = (-t294 * t339 - t295 * t342) * qJD(2);
t256 = rSges(6,1) * t299 + rSges(6,2) * t298 + rSges(6,3) * t381;
t255 = rSges(6,1) * t297 + rSges(6,2) * t296 + rSges(6,3) * t382;
t254 = Icges(6,1) * t299 + Icges(6,4) * t298 + Icges(6,5) * t381;
t253 = Icges(6,1) * t297 + Icges(6,4) * t296 + Icges(6,5) * t382;
t252 = Icges(6,4) * t299 + Icges(6,2) * t298 + Icges(6,6) * t381;
t251 = Icges(6,4) * t297 + Icges(6,2) * t296 + Icges(6,6) * t382;
t250 = Icges(6,5) * t299 + Icges(6,6) * t298 + Icges(6,3) * t381;
t249 = Icges(6,5) * t297 + Icges(6,6) * t296 + Icges(6,3) * t382;
t248 = qJD(1) * t285 - t314 * t316 + t352;
t247 = -t370 + t314 * t317 + (-t284 + t371) * qJD(1);
t246 = t284 * t316 - t285 * t317 + t351;
t245 = qJD(1) * t275 - t305 * t307 + t347;
t244 = t305 * t308 + (-t274 + t369) * qJD(1) + t350;
t243 = t274 * t307 - t275 * t308 + t346;
t242 = qJD(1) * t287 + t256 * t318 - t264 * t272 - t306 * t307 + t347;
t241 = -t255 * t318 + t264 * t273 + t306 * t308 + (-t286 + t369) * qJD(1) + t350;
t240 = t255 * t272 - t256 * t273 + t286 * t307 - t287 * t308 + t346;
t1 = m(3) * (t257 ^ 2 + t258 ^ 2 + t259 ^ 2) / 0.2e1 - ((-t339 * t319 + t353 * t342) * qJD(1) + (t339 ^ 2 * t289 + (t355 * t342 + (-t288 + t354) * t339) * t342) * qJD(2)) * t375 / 0.2e1 + ((t342 * t319 + t353 * t339) * qJD(1) + (t342 ^ 2 * t288 + (t354 * t339 + (-t289 + t355) * t342) * t339) * qJD(2)) * t331 / 0.2e1 + m(4) * (t246 ^ 2 + t247 ^ 2 + t248 ^ 2) / 0.2e1 + t316 * (-t348 * t339 + t344 * t342) / 0.2e1 + t317 * (t344 * t339 + t348 * t342) / 0.2e1 + m(5) * (t243 ^ 2 + t244 ^ 2 + t245 ^ 2) / 0.2e1 + t307 * (-t349 * t339 + t345 * t342) / 0.2e1 + t308 * (t345 * t339 + t349 * t342) / 0.2e1 + m(6) * (t240 ^ 2 + t241 ^ 2 + t242 ^ 2) / 0.2e1 + t272 * ((t250 * t381 + t298 * t252 + t299 * t254) * t272 + (t249 * t381 + t251 * t298 + t253 * t299) * t273 + (t261 * t381 + t262 * t298 + t263 * t299) * t318) / 0.2e1 + t273 * ((t250 * t382 + t252 * t296 + t254 * t297) * t272 + (t249 * t382 + t296 * t251 + t297 * t253) * t273 + (t261 * t382 + t262 * t296 + t263 * t297) * t318) / 0.2e1 + t318 * ((t249 * t273 + t250 * t272 + t261 * t318) * t326 + ((t252 * t337 - t254 * t340) * t272 + (t251 * t337 - t253 * t340) * t273 + (t262 * t337 - t263 * t340) * t318) * t325) / 0.2e1 + (m(2) * (t323 ^ 2 + t324 ^ 2) + Icges(2,3)) * qJD(1) ^ 2 / 0.2e1 + ((-(-t291 * t341 - t293 * t338) * t339 + (-t290 * t341 - t292 * t338) * t342) * qJD(2) + (-t281 * t333 - t283 * t332) * t316 + (-t280 * t333 - t282 * t332) * t317 + (-t269 * t326 - t271 * t325) * t307 + (-t268 * t326 - t270 * t325) * t308 + (-t326 * t302 - t325 * t303 - t333 * t312 - t332 * t313 - t341 * t320 - t338 * t321) * qJD(1)) * qJD(1) / 0.2e1;
T  = t1;
