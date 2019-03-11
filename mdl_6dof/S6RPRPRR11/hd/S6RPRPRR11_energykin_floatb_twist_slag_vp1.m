% Calculate kinetic energy for
% S6RPRPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 04:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR11_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR11_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRR11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR11_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR11_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR11_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR11_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:11:11
% EndTime: 2019-03-09 04:11:16
% DurationCPUTime: 5.23s
% Computational Cost: add. (5245->447), mult. (12423->633), div. (0->0), fcn. (15815->16), ass. (0->193)
t424 = Icges(4,2) + Icges(5,3);
t369 = cos(pkin(12));
t366 = sin(pkin(12));
t374 = sin(qJ(1));
t403 = t374 * t366;
t370 = cos(pkin(6));
t376 = cos(qJ(1));
t404 = t370 * t376;
t336 = t369 * t404 - t403;
t402 = t374 * t369;
t337 = t366 * t404 + t402;
t373 = sin(qJ(3));
t367 = sin(pkin(6));
t413 = sin(pkin(7));
t393 = t367 * t413;
t414 = cos(pkin(7));
t416 = cos(qJ(3));
t294 = t337 * t416 + (t336 * t414 - t376 * t393) * t373;
t394 = t367 * t414;
t318 = -t336 * t413 - t376 * t394;
t365 = sin(pkin(13));
t368 = cos(pkin(13));
t272 = -t294 * t365 + t318 * t368;
t410 = t318 * t365;
t273 = t294 * t368 + t410;
t389 = t416 * t413;
t387 = t367 * t389;
t390 = t414 * t416;
t293 = -t336 * t390 + t337 * t373 + t376 * t387;
t423 = -Icges(4,4) * t294 + Icges(5,5) * t273 - Icges(4,6) * t318 + Icges(5,6) * t272 + t293 * t424;
t338 = -t366 * t376 - t370 * t402;
t339 = t369 * t376 - t370 * t403;
t296 = t339 * t416 + (t338 * t414 + t374 * t393) * t373;
t319 = -t338 * t413 + t374 * t394;
t274 = -t296 * t365 + t319 * t368;
t409 = t319 * t365;
t275 = t296 * t368 + t409;
t295 = -t338 * t390 + t339 * t373 - t374 * t387;
t422 = -Icges(4,4) * t296 + Icges(5,5) * t275 - Icges(4,6) * t319 + Icges(5,6) * t274 + t295 * t424;
t317 = t370 * t413 * t373 + (t369 * t373 * t414 + t366 * t416) * t367;
t335 = -t369 * t393 + t370 * t414;
t290 = -t317 * t365 + t335 * t368;
t408 = t335 * t365;
t291 = t317 * t368 + t408;
t407 = t366 * t367;
t316 = -t367 * t369 * t390 - t370 * t389 + t373 * t407;
t421 = -Icges(4,4) * t317 + Icges(5,5) * t291 - Icges(4,6) * t335 + Icges(5,6) * t290 + t316 * t424;
t415 = pkin(4) * t368;
t412 = Icges(2,4) * t374;
t411 = qJ(2) * t370;
t406 = t367 * t374;
t405 = t367 * t376;
t400 = qJD(2) * t367;
t399 = pkin(13) + qJ(5);
t398 = V_base(5) * pkin(8) + V_base(1);
t311 = qJD(3) * t319 + V_base(4);
t310 = qJD(3) * t318 + V_base(5);
t362 = V_base(6) + qJD(1);
t395 = -pkin(8) - t411;
t392 = cos(t399);
t341 = t374 * pkin(1) - qJ(2) * t405;
t391 = qJD(2) * t370 + V_base(4) * t341 + V_base(3);
t277 = qJD(5) * t295 + t311;
t276 = qJD(5) * t293 + t310;
t325 = qJD(3) * t335 + t362;
t388 = t374 * t400 + V_base(5) * t411 + t398;
t285 = qJD(5) * t316 + t325;
t342 = pkin(1) * t376 + qJ(2) * t406;
t386 = t362 * t342 - t376 * t400 + V_base(2);
t298 = t337 * pkin(2) + pkin(9) * t318;
t322 = pkin(2) * t407 + pkin(9) * t335;
t385 = V_base(5) * t322 + (-t298 - t341) * t362 + t388;
t299 = t339 * pkin(2) + pkin(9) * t319;
t384 = V_base(4) * t298 + (-t299 - t342) * V_base(5) + t391;
t282 = pkin(3) * t317 + qJ(4) * t316;
t383 = qJD(4) * t295 + t310 * t282 + t385;
t261 = pkin(3) * t294 + qJ(4) * t293;
t382 = qJD(4) * t316 + t311 * t261 + t384;
t381 = t362 * t299 + (-t322 + t395) * V_base(4) + t386;
t262 = pkin(3) * t296 + qJ(4) * t295;
t380 = qJD(4) * t293 + t325 * t262 + t381;
t208 = pkin(4) * t410 + pkin(10) * t293 + t294 * t415;
t235 = pkin(4) * t408 + pkin(10) * t316 + t317 * t415;
t379 = t310 * t235 + (-t208 - t261) * t325 + t383;
t209 = pkin(4) * t409 + pkin(10) * t295 + t296 * t415;
t378 = t311 * t208 + (-t209 - t262) * t310 + t382;
t377 = t325 * t209 + (-t235 - t282) * t311 + t380;
t375 = cos(qJ(6));
t372 = sin(qJ(6));
t363 = Icges(2,4) * t376;
t361 = sin(t399);
t355 = rSges(2,1) * t376 - t374 * rSges(2,2);
t354 = t374 * rSges(2,1) + rSges(2,2) * t376;
t353 = Icges(2,1) * t376 - t412;
t352 = Icges(2,1) * t374 + t363;
t351 = -Icges(2,2) * t374 + t363;
t350 = Icges(2,2) * t376 + t412;
t349 = Icges(2,5) * t376 - Icges(2,6) * t374;
t348 = Icges(2,5) * t374 + Icges(2,6) * t376;
t347 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t346 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t345 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t332 = rSges(3,3) * t370 + (rSges(3,1) * t366 + rSges(3,2) * t369) * t367;
t331 = Icges(3,5) * t370 + (Icges(3,1) * t366 + Icges(3,4) * t369) * t367;
t330 = Icges(3,6) * t370 + (Icges(3,4) * t366 + Icges(3,2) * t369) * t367;
t329 = Icges(3,3) * t370 + (Icges(3,5) * t366 + Icges(3,6) * t369) * t367;
t324 = V_base(5) * rSges(2,3) - t354 * t362 + t398;
t323 = t355 * t362 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t321 = t354 * V_base(4) - t355 * V_base(5) + V_base(3);
t307 = rSges(3,1) * t339 + rSges(3,2) * t338 + rSges(3,3) * t406;
t306 = t337 * rSges(3,1) + t336 * rSges(3,2) - rSges(3,3) * t405;
t305 = Icges(3,1) * t339 + Icges(3,4) * t338 + Icges(3,5) * t406;
t304 = Icges(3,1) * t337 + Icges(3,4) * t336 - Icges(3,5) * t405;
t303 = Icges(3,4) * t339 + Icges(3,2) * t338 + Icges(3,6) * t406;
t302 = Icges(3,4) * t337 + Icges(3,2) * t336 - Icges(3,6) * t405;
t301 = Icges(3,5) * t339 + Icges(3,6) * t338 + Icges(3,3) * t406;
t300 = Icges(3,5) * t337 + Icges(3,6) * t336 - Icges(3,3) * t405;
t284 = t317 * t392 + t335 * t361;
t283 = t317 * t361 - t335 * t392;
t281 = rSges(4,1) * t317 - rSges(4,2) * t316 + rSges(4,3) * t335;
t280 = Icges(4,1) * t317 - Icges(4,4) * t316 + Icges(4,5) * t335;
t278 = Icges(4,5) * t317 - Icges(4,6) * t316 + Icges(4,3) * t335;
t271 = t296 * t392 + t319 * t361;
t270 = t296 * t361 - t319 * t392;
t269 = t294 * t392 + t318 * t361;
t268 = t294 * t361 - t318 * t392;
t267 = t284 * t375 + t316 * t372;
t266 = -t284 * t372 + t316 * t375;
t264 = t332 * V_base(5) + (-t306 - t341) * t362 + t388;
t263 = t362 * t307 + (-t332 + t395) * V_base(4) + t386;
t260 = qJD(6) * t283 + t285;
t259 = pkin(5) * t284 + pkin(11) * t283;
t258 = t306 * V_base(4) + (-t307 - t342) * V_base(5) + t391;
t256 = rSges(4,1) * t296 - rSges(4,2) * t295 + rSges(4,3) * t319;
t255 = rSges(4,1) * t294 - rSges(4,2) * t293 + rSges(4,3) * t318;
t254 = Icges(4,1) * t296 - Icges(4,4) * t295 + Icges(4,5) * t319;
t253 = Icges(4,1) * t294 - Icges(4,4) * t293 + Icges(4,5) * t318;
t250 = Icges(4,5) * t296 - Icges(4,6) * t295 + Icges(4,3) * t319;
t249 = Icges(4,5) * t294 - Icges(4,6) * t293 + Icges(4,3) * t318;
t247 = rSges(5,1) * t291 + rSges(5,2) * t290 + rSges(5,3) * t316;
t246 = Icges(5,1) * t291 + Icges(5,4) * t290 + Icges(5,5) * t316;
t245 = Icges(5,4) * t291 + Icges(5,2) * t290 + Icges(5,6) * t316;
t243 = rSges(6,1) * t284 - rSges(6,2) * t283 + rSges(6,3) * t316;
t242 = t271 * t375 + t295 * t372;
t241 = -t271 * t372 + t295 * t375;
t240 = t269 * t375 + t293 * t372;
t239 = -t269 * t372 + t293 * t375;
t238 = Icges(6,1) * t284 - Icges(6,4) * t283 + Icges(6,5) * t316;
t237 = Icges(6,4) * t284 - Icges(6,2) * t283 + Icges(6,6) * t316;
t236 = Icges(6,5) * t284 - Icges(6,6) * t283 + Icges(6,3) * t316;
t234 = qJD(6) * t270 + t277;
t233 = qJD(6) * t268 + t276;
t232 = pkin(5) * t271 + pkin(11) * t270;
t231 = pkin(5) * t269 + pkin(11) * t268;
t229 = rSges(5,1) * t275 + rSges(5,2) * t274 + rSges(5,3) * t295;
t228 = rSges(5,1) * t273 + rSges(5,2) * t272 + rSges(5,3) * t293;
t227 = Icges(5,1) * t275 + Icges(5,4) * t274 + Icges(5,5) * t295;
t226 = Icges(5,1) * t273 + Icges(5,4) * t272 + Icges(5,5) * t293;
t225 = Icges(5,4) * t275 + Icges(5,2) * t274 + Icges(5,6) * t295;
t224 = Icges(5,4) * t273 + Icges(5,2) * t272 + Icges(5,6) * t293;
t221 = rSges(6,1) * t271 - rSges(6,2) * t270 + rSges(6,3) * t295;
t220 = rSges(6,1) * t269 - rSges(6,2) * t268 + rSges(6,3) * t293;
t219 = Icges(6,1) * t271 - Icges(6,4) * t270 + Icges(6,5) * t295;
t218 = Icges(6,1) * t269 - Icges(6,4) * t268 + Icges(6,5) * t293;
t217 = Icges(6,4) * t271 - Icges(6,2) * t270 + Icges(6,6) * t295;
t216 = Icges(6,4) * t269 - Icges(6,2) * t268 + Icges(6,6) * t293;
t215 = Icges(6,5) * t271 - Icges(6,6) * t270 + Icges(6,3) * t295;
t214 = Icges(6,5) * t269 - Icges(6,6) * t268 + Icges(6,3) * t293;
t213 = rSges(7,1) * t267 + rSges(7,2) * t266 + rSges(7,3) * t283;
t212 = Icges(7,1) * t267 + Icges(7,4) * t266 + Icges(7,5) * t283;
t211 = Icges(7,4) * t267 + Icges(7,2) * t266 + Icges(7,6) * t283;
t210 = Icges(7,5) * t267 + Icges(7,6) * t266 + Icges(7,3) * t283;
t205 = -t255 * t325 + t281 * t310 + t385;
t204 = t325 * t256 - t311 * t281 + t381;
t203 = rSges(7,1) * t242 + rSges(7,2) * t241 + rSges(7,3) * t270;
t202 = rSges(7,1) * t240 + rSges(7,2) * t239 + rSges(7,3) * t268;
t201 = Icges(7,1) * t242 + Icges(7,4) * t241 + Icges(7,5) * t270;
t200 = Icges(7,1) * t240 + Icges(7,4) * t239 + Icges(7,5) * t268;
t199 = Icges(7,4) * t242 + Icges(7,2) * t241 + Icges(7,6) * t270;
t198 = Icges(7,4) * t240 + Icges(7,2) * t239 + Icges(7,6) * t268;
t197 = Icges(7,5) * t242 + Icges(7,6) * t241 + Icges(7,3) * t270;
t196 = Icges(7,5) * t240 + Icges(7,6) * t239 + Icges(7,3) * t268;
t195 = t255 * t311 - t256 * t310 + t384;
t194 = t247 * t310 + (-t228 - t261) * t325 + t383;
t193 = t325 * t229 + (-t247 - t282) * t311 + t380;
t192 = t228 * t311 + (-t229 - t262) * t310 + t382;
t191 = -t220 * t285 + t243 * t276 + t379;
t190 = t285 * t221 - t277 * t243 + t377;
t189 = t220 * t277 - t221 * t276 + t378;
t188 = -t202 * t260 + t213 * t233 - t231 * t285 + t259 * t276 + t379;
t187 = t260 * t203 - t234 * t213 + t285 * t232 - t277 * t259 + t377;
t186 = t202 * t234 - t203 * t233 + t231 * t277 - t232 * t276 + t378;
t1 = m(1) * (t345 ^ 2 + t346 ^ 2 + t347 ^ 2) / 0.2e1 + m(2) * (t321 ^ 2 + t323 ^ 2 + t324 ^ 2) / 0.2e1 + t285 * ((t215 * t316 - t217 * t283 + t219 * t284) * t277 + (t214 * t316 - t216 * t283 + t218 * t284) * t276 + (t236 * t316 - t237 * t283 + t238 * t284) * t285) / 0.2e1 + t277 * ((t215 * t295 - t217 * t270 + t219 * t271) * t277 + (t214 * t295 - t216 * t270 + t218 * t271) * t276 + (t236 * t295 - t237 * t270 + t238 * t271) * t285) / 0.2e1 + t276 * ((t215 * t293 - t217 * t268 + t219 * t269) * t277 + (t214 * t293 - t216 * t268 + t218 * t269) * t276 + (t236 * t293 - t237 * t268 + t238 * t269) * t285) / 0.2e1 + t260 * ((t197 * t283 + t199 * t266 + t201 * t267) * t234 + (t196 * t283 + t198 * t266 + t200 * t267) * t233 + (t210 * t283 + t211 * t266 + t212 * t267) * t260) / 0.2e1 + m(7) * (t186 ^ 2 + t187 ^ 2 + t188 ^ 2) / 0.2e1 + m(6) * (t189 ^ 2 + t190 ^ 2 + t191 ^ 2) / 0.2e1 + m(5) * (t192 ^ 2 + t193 ^ 2 + t194 ^ 2) / 0.2e1 + m(4) * (t195 ^ 2 + t204 ^ 2 + t205 ^ 2) / 0.2e1 + m(3) * (t258 ^ 2 + t263 ^ 2 + t264 ^ 2) / 0.2e1 + t233 * ((t197 * t268 + t199 * t239 + t201 * t240) * t234 + (t196 * t268 + t198 * t239 + t200 * t240) * t233 + (t210 * t268 + t211 * t239 + t212 * t240) * t260) / 0.2e1 + t234 * ((t197 * t270 + t199 * t241 + t201 * t242) * t234 + (t196 * t270 + t198 * t241 + t200 * t242) * t233 + (t210 * t270 + t211 * t241 + t212 * t242) * t260) / 0.2e1 + ((t245 * t272 + t246 * t273 + t278 * t318 + t280 * t294 + t293 * t421) * t325 + (t225 * t272 + t227 * t273 + t250 * t318 + t254 * t294 + t293 * t422) * t311 + (t224 * t272 + t226 * t273 + t249 * t318 + t253 * t294 + t293 * t423) * t310) * t310 / 0.2e1 + ((t245 * t274 + t246 * t275 + t278 * t319 + t280 * t296 + t295 * t421) * t325 + (t225 * t274 + t227 * t275 + t250 * t319 + t254 * t296 + t295 * t422) * t311 + (t224 * t274 + t226 * t275 + t249 * t319 + t253 * t296 + t295 * t423) * t310) * t311 / 0.2e1 + ((t245 * t290 + t246 * t291 + t278 * t335 + t280 * t317 + t316 * t421) * t325 + (t225 * t290 + t227 * t291 + t250 * t335 + t254 * t317 + t316 * t422) * t311 + (t224 * t290 + t226 * t291 + t249 * t335 + t253 * t317 + t316 * t423) * t310) * t325 / 0.2e1 + ((t300 * V_base(5) + t301 * V_base(4) + t329 * t362) * t370 + ((t303 * t369 + t305 * t366) * V_base(4) + (t302 * t369 + t304 * t366) * V_base(5) + (t330 * t369 + t331 * t366) * t362) * t367 + Icges(2,3) * t362 + t348 * V_base(5) + t349 * V_base(4)) * t362 / 0.2e1 + ((t329 * t406 + t330 * t338 + t331 * t339 + t349) * t362 + (t300 * t406 + t302 * t338 + t304 * t339 - t374 * t350 + t352 * t376 + Icges(1,4)) * V_base(5) + (t301 * t406 + t303 * t338 + t305 * t339 - t374 * t351 + t353 * t376 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t329 * t405 + t336 * t330 + t337 * t331 + t348) * t362 + (-t300 * t405 + t336 * t302 + t337 * t304 + t350 * t376 + t374 * t352 + Icges(1,2)) * V_base(5) + (-t301 * t405 + t336 * t303 + t337 * t305 + t351 * t376 + t374 * t353 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
