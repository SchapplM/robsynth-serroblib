% Calculate kinetic energy for
% S6PPRRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PPRRRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRRR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:59:38
% EndTime: 2019-03-08 18:59:42
% DurationCPUTime: 3.63s
% Computational Cost: add. (5458->452), mult. (13137->657), div. (0->0), fcn. (16739->16), ass. (0->197)
t367 = sin(pkin(12));
t370 = cos(pkin(12));
t420 = Icges(2,5) * t370 - Icges(2,6) * t367 + Icges(1,5);
t419 = Icges(2,5) * t367 + Icges(2,6) * t370 + Icges(1,6);
t418 = cos(qJ(3));
t376 = cos(qJ(4));
t417 = pkin(4) * t376;
t415 = cos(pkin(7));
t414 = sin(pkin(7));
t413 = Icges(2,4) * t367;
t371 = cos(pkin(6));
t412 = qJ(2) * t371;
t366 = sin(pkin(13));
t369 = cos(pkin(13));
t404 = t370 * t371;
t337 = -t366 * t367 + t369 * t404;
t368 = sin(pkin(6));
t395 = t368 * t415;
t321 = -t337 * t414 - t370 * t395;
t373 = sin(qJ(4));
t411 = t321 * t373;
t406 = t367 * t371;
t339 = -t366 * t370 - t369 * t406;
t322 = -t339 * t414 + t367 * t395;
t410 = t322 * t373;
t394 = t368 * t414;
t336 = -t369 * t394 + t371 * t415;
t409 = t336 * t373;
t408 = t366 * t368;
t407 = t367 * t368;
t405 = t368 * t370;
t403 = qJ(4) + qJ(5);
t402 = qJD(2) * t368;
t401 = V_base(5) * qJ(1) + V_base(1);
t397 = qJD(1) + V_base(3);
t311 = qJD(3) * t322 + V_base(4);
t310 = qJD(3) * t321 + V_base(5);
t329 = qJD(3) * t336 + V_base(6);
t396 = cos(t403);
t393 = -qJ(1) - t412;
t340 = -t366 * t406 + t369 * t370;
t374 = sin(qJ(3));
t391 = t418 * t414;
t388 = t368 * t391;
t392 = t415 * t418;
t293 = -t339 * t392 + t340 * t374 - t367 * t388;
t274 = qJD(4) * t293 + t311;
t338 = t366 * t404 + t367 * t369;
t291 = -t337 * t392 + t338 * t374 + t370 * t388;
t273 = qJD(4) * t291 + t310;
t318 = -t368 * t369 * t392 - t371 * t391 + t374 * t408;
t295 = qJD(4) * t318 + t329;
t390 = t367 * t402 + V_base(5) * t412 + t401;
t343 = pkin(1) * t367 - qJ(2) * t405;
t389 = qJD(2) * t371 + V_base(4) * t343 + t397;
t245 = qJD(5) * t293 + t274;
t244 = qJD(5) * t291 + t273;
t279 = qJD(5) * t318 + t295;
t344 = pkin(1) * t370 + qJ(2) * t407;
t387 = V_base(6) * t344 - t370 * t402 + V_base(2);
t300 = t338 * pkin(2) + pkin(8) * t321;
t324 = pkin(2) * t408 + pkin(8) * t336;
t386 = V_base(5) * t324 + (-t300 - t343) * V_base(6) + t390;
t301 = t340 * pkin(2) + pkin(8) * t322;
t385 = V_base(4) * t300 + (-t301 - t344) * V_base(5) + t389;
t292 = t338 * t418 + (t337 * t415 - t370 * t394) * t374;
t262 = t292 * pkin(3) + t291 * pkin(9);
t319 = t371 * t414 * t374 + (t415 * t369 * t374 + t366 * t418) * t368;
t284 = t319 * pkin(3) + t318 * pkin(9);
t384 = -t262 * t329 + t310 * t284 + t386;
t294 = t340 * t418 + (t339 * t415 + t367 * t394) * t374;
t263 = t294 * pkin(3) + t293 * pkin(9);
t383 = t311 * t262 - t263 * t310 + t385;
t382 = V_base(6) * t301 + (-t324 + t393) * V_base(4) + t387;
t207 = pkin(4) * t411 + pkin(10) * t291 + t292 * t417;
t238 = pkin(4) * t409 + pkin(10) * t318 + t319 * t417;
t381 = -t207 * t295 + t273 * t238 + t384;
t208 = pkin(4) * t410 + pkin(10) * t293 + t294 * t417;
t380 = t274 * t207 - t208 * t273 + t383;
t379 = t329 * t263 - t284 * t311 + t382;
t378 = t295 * t208 - t238 * t274 + t379;
t375 = cos(qJ(6));
t372 = sin(qJ(6));
t364 = sin(t403);
t363 = Icges(2,4) * t370;
t357 = rSges(2,1) * t370 - rSges(2,2) * t367;
t356 = rSges(2,1) * t367 + rSges(2,2) * t370;
t355 = Icges(2,1) * t370 - t413;
t354 = Icges(2,1) * t367 + t363;
t353 = -Icges(2,2) * t367 + t363;
t352 = Icges(2,2) * t370 + t413;
t349 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t348 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t347 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t334 = rSges(3,3) * t371 + (rSges(3,1) * t366 + rSges(3,2) * t369) * t368;
t333 = Icges(3,5) * t371 + (Icges(3,1) * t366 + Icges(3,4) * t369) * t368;
t332 = Icges(3,6) * t371 + (Icges(3,4) * t366 + Icges(3,2) * t369) * t368;
t331 = Icges(3,3) * t371 + (Icges(3,5) * t366 + Icges(3,6) * t369) * t368;
t326 = V_base(5) * rSges(2,3) - t356 * V_base(6) + t401;
t325 = t357 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t320 = t356 * V_base(4) - t357 * V_base(5) + t397;
t309 = rSges(3,1) * t340 + rSges(3,2) * t339 + rSges(3,3) * t407;
t308 = rSges(3,1) * t338 + rSges(3,2) * t337 - rSges(3,3) * t405;
t307 = Icges(3,1) * t340 + Icges(3,4) * t339 + Icges(3,5) * t407;
t306 = Icges(3,1) * t338 + Icges(3,4) * t337 - Icges(3,5) * t405;
t305 = Icges(3,4) * t340 + Icges(3,2) * t339 + Icges(3,6) * t407;
t304 = Icges(3,4) * t338 + Icges(3,2) * t337 - Icges(3,6) * t405;
t303 = Icges(3,5) * t340 + Icges(3,6) * t339 + Icges(3,3) * t407;
t302 = Icges(3,5) * t338 + Icges(3,6) * t337 - Icges(3,3) * t405;
t297 = t319 * t376 + t409;
t296 = -t319 * t373 + t336 * t376;
t286 = t319 * t396 + t336 * t364;
t285 = t319 * t364 - t336 * t396;
t283 = rSges(4,1) * t319 - rSges(4,2) * t318 + rSges(4,3) * t336;
t282 = Icges(4,1) * t319 - Icges(4,4) * t318 + Icges(4,5) * t336;
t281 = Icges(4,4) * t319 - Icges(4,2) * t318 + Icges(4,6) * t336;
t280 = Icges(4,5) * t319 - Icges(4,6) * t318 + Icges(4,3) * t336;
t278 = t294 * t376 + t410;
t277 = -t294 * t373 + t322 * t376;
t276 = t292 * t376 + t411;
t275 = -t292 * t373 + t321 * t376;
t272 = t286 * t375 + t318 * t372;
t271 = -t286 * t372 + t318 * t375;
t270 = t294 * t396 + t322 * t364;
t269 = t294 * t364 - t322 * t396;
t268 = t292 * t396 + t321 * t364;
t267 = t292 * t364 - t321 * t396;
t265 = t334 * V_base(5) + (-t308 - t343) * V_base(6) + t390;
t264 = t309 * V_base(6) + (-t334 + t393) * V_base(4) + t387;
t261 = pkin(5) * t286 + pkin(11) * t285;
t259 = t308 * V_base(4) + (-t309 - t344) * V_base(5) + t389;
t258 = rSges(5,1) * t297 + rSges(5,2) * t296 + rSges(5,3) * t318;
t257 = Icges(5,1) * t297 + Icges(5,4) * t296 + Icges(5,5) * t318;
t256 = Icges(5,4) * t297 + Icges(5,2) * t296 + Icges(5,6) * t318;
t255 = Icges(5,5) * t297 + Icges(5,6) * t296 + Icges(5,3) * t318;
t254 = rSges(4,1) * t294 - rSges(4,2) * t293 + rSges(4,3) * t322;
t253 = rSges(4,1) * t292 - rSges(4,2) * t291 + rSges(4,3) * t321;
t252 = Icges(4,1) * t294 - Icges(4,4) * t293 + Icges(4,5) * t322;
t251 = Icges(4,1) * t292 - Icges(4,4) * t291 + Icges(4,5) * t321;
t250 = Icges(4,4) * t294 - Icges(4,2) * t293 + Icges(4,6) * t322;
t249 = Icges(4,4) * t292 - Icges(4,2) * t291 + Icges(4,6) * t321;
t248 = Icges(4,5) * t294 - Icges(4,6) * t293 + Icges(4,3) * t322;
t247 = Icges(4,5) * t292 - Icges(4,6) * t291 + Icges(4,3) * t321;
t246 = qJD(6) * t285 + t279;
t242 = rSges(6,1) * t286 - rSges(6,2) * t285 + rSges(6,3) * t318;
t241 = Icges(6,1) * t286 - Icges(6,4) * t285 + Icges(6,5) * t318;
t240 = Icges(6,4) * t286 - Icges(6,2) * t285 + Icges(6,6) * t318;
t239 = Icges(6,5) * t286 - Icges(6,6) * t285 + Icges(6,3) * t318;
t237 = t270 * t375 + t293 * t372;
t236 = -t270 * t372 + t293 * t375;
t235 = t268 * t375 + t291 * t372;
t234 = -t268 * t372 + t291 * t375;
t233 = pkin(5) * t270 + pkin(11) * t269;
t232 = pkin(5) * t268 + pkin(11) * t267;
t231 = rSges(5,1) * t278 + rSges(5,2) * t277 + rSges(5,3) * t293;
t230 = rSges(5,1) * t276 + rSges(5,2) * t275 + rSges(5,3) * t291;
t229 = Icges(5,1) * t278 + Icges(5,4) * t277 + Icges(5,5) * t293;
t228 = Icges(5,1) * t276 + Icges(5,4) * t275 + Icges(5,5) * t291;
t227 = Icges(5,4) * t278 + Icges(5,2) * t277 + Icges(5,6) * t293;
t226 = Icges(5,4) * t276 + Icges(5,2) * t275 + Icges(5,6) * t291;
t225 = Icges(5,5) * t278 + Icges(5,6) * t277 + Icges(5,3) * t293;
t224 = Icges(5,5) * t276 + Icges(5,6) * t275 + Icges(5,3) * t291;
t223 = qJD(6) * t269 + t245;
t222 = qJD(6) * t267 + t244;
t220 = rSges(6,1) * t270 - rSges(6,2) * t269 + rSges(6,3) * t293;
t219 = rSges(6,1) * t268 - rSges(6,2) * t267 + rSges(6,3) * t291;
t218 = Icges(6,1) * t270 - Icges(6,4) * t269 + Icges(6,5) * t293;
t217 = Icges(6,1) * t268 - Icges(6,4) * t267 + Icges(6,5) * t291;
t216 = Icges(6,4) * t270 - Icges(6,2) * t269 + Icges(6,6) * t293;
t215 = Icges(6,4) * t268 - Icges(6,2) * t267 + Icges(6,6) * t291;
t214 = Icges(6,5) * t270 - Icges(6,6) * t269 + Icges(6,3) * t293;
t213 = Icges(6,5) * t268 - Icges(6,6) * t267 + Icges(6,3) * t291;
t212 = rSges(7,1) * t272 + rSges(7,2) * t271 + rSges(7,3) * t285;
t211 = Icges(7,1) * t272 + Icges(7,4) * t271 + Icges(7,5) * t285;
t210 = Icges(7,4) * t272 + Icges(7,2) * t271 + Icges(7,6) * t285;
t209 = Icges(7,5) * t272 + Icges(7,6) * t271 + Icges(7,3) * t285;
t204 = -t253 * t329 + t283 * t310 + t386;
t203 = t254 * t329 - t283 * t311 + t382;
t202 = rSges(7,1) * t237 + rSges(7,2) * t236 + rSges(7,3) * t269;
t201 = rSges(7,1) * t235 + rSges(7,2) * t234 + rSges(7,3) * t267;
t200 = Icges(7,1) * t237 + Icges(7,4) * t236 + Icges(7,5) * t269;
t199 = Icges(7,1) * t235 + Icges(7,4) * t234 + Icges(7,5) * t267;
t198 = Icges(7,4) * t237 + Icges(7,2) * t236 + Icges(7,6) * t269;
t197 = Icges(7,4) * t235 + Icges(7,2) * t234 + Icges(7,6) * t267;
t196 = Icges(7,5) * t237 + Icges(7,6) * t236 + Icges(7,3) * t269;
t195 = Icges(7,5) * t235 + Icges(7,6) * t234 + Icges(7,3) * t267;
t194 = t253 * t311 - t254 * t310 + t385;
t193 = -t230 * t295 + t258 * t273 + t384;
t192 = t231 * t295 - t258 * t274 + t379;
t191 = t230 * t274 - t231 * t273 + t383;
t190 = -t219 * t279 + t242 * t244 + t381;
t189 = t220 * t279 - t242 * t245 + t378;
t188 = t219 * t245 - t220 * t244 + t380;
t187 = -t201 * t246 + t212 * t222 - t232 * t279 + t244 * t261 + t381;
t186 = t202 * t246 - t212 * t223 + t233 * t279 - t245 * t261 + t378;
t185 = t201 * t223 - t202 * t222 + t232 * t245 - t233 * t244 + t380;
t1 = m(1) * (t347 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 + t329 * ((t248 * t336 - t250 * t318 + t252 * t319) * t311 + (t247 * t336 - t249 * t318 + t251 * t319) * t310 + (t280 * t336 - t281 * t318 + t282 * t319) * t329) / 0.2e1 + m(2) * (t320 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + t311 * ((t248 * t322 - t250 * t293 + t252 * t294) * t311 + (t247 * t322 - t249 * t293 + t251 * t294) * t310 + (t280 * t322 - t281 * t293 + t282 * t294) * t329) / 0.2e1 + t310 * ((t248 * t321 - t250 * t291 + t252 * t292) * t311 + (t247 * t321 - t249 * t291 + t251 * t292) * t310 + (t280 * t321 - t281 * t291 + t282 * t292) * t329) / 0.2e1 + t295 * ((t225 * t318 + t227 * t296 + t229 * t297) * t274 + (t224 * t318 + t226 * t296 + t228 * t297) * t273 + (t255 * t318 + t256 * t296 + t257 * t297) * t295) / 0.2e1 + t279 * ((t214 * t318 - t216 * t285 + t218 * t286) * t245 + (t213 * t318 - t215 * t285 + t217 * t286) * t244 + (t239 * t318 - t240 * t285 + t241 * t286) * t279) / 0.2e1 + t274 * ((t225 * t293 + t227 * t277 + t229 * t278) * t274 + (t224 * t293 + t226 * t277 + t228 * t278) * t273 + (t255 * t293 + t256 * t277 + t257 * t278) * t295) / 0.2e1 + t273 * ((t225 * t291 + t227 * t275 + t229 * t276) * t274 + (t291 * t224 + t275 * t226 + t276 * t228) * t273 + (t255 * t291 + t256 * t275 + t257 * t276) * t295) / 0.2e1 + t245 * ((t293 * t214 - t269 * t216 + t270 * t218) * t245 + (t213 * t293 - t215 * t269 + t217 * t270) * t244 + (t239 * t293 - t240 * t269 + t241 * t270) * t279) / 0.2e1 + t244 * ((t214 * t291 - t216 * t267 + t218 * t268) * t245 + (t291 * t213 - t267 * t215 + t268 * t217) * t244 + (t239 * t291 - t240 * t267 + t241 * t268) * t279) / 0.2e1 + t246 * ((t196 * t285 + t198 * t271 + t200 * t272) * t223 + (t195 * t285 + t197 * t271 + t199 * t272) * t222 + (t285 * t209 + t271 * t210 + t272 * t211) * t246) / 0.2e1 + t223 * ((t269 * t196 + t236 * t198 + t237 * t200) * t223 + (t195 * t269 + t197 * t236 + t199 * t237) * t222 + (t209 * t269 + t210 * t236 + t211 * t237) * t246) / 0.2e1 + m(7) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(6) * (t188 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + m(5) * (t191 ^ 2 + t192 ^ 2 + t193 ^ 2) / 0.2e1 + m(4) * (t194 ^ 2 + t203 ^ 2 + t204 ^ 2) / 0.2e1 + m(3) * (t259 ^ 2 + t264 ^ 2 + t265 ^ 2) / 0.2e1 + t222 * ((t196 * t267 + t198 * t234 + t200 * t235) * t223 + (t267 * t195 + t234 * t197 + t235 * t199) * t222 + (t209 * t267 + t210 * t234 + t211 * t235) * t246) / 0.2e1 + ((t331 * t407 + t332 * t339 + t333 * t340 + t420) * V_base(6) + (t302 * t407 + t304 * t339 + t306 * t340 - t352 * t367 + t354 * t370 + Icges(1,4)) * V_base(5) + (t303 * t407 + t305 * t339 + t307 * t340 - t353 * t367 + t355 * t370 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t331 * t405 + t332 * t337 + t333 * t338 + t419) * V_base(6) + (-t302 * t405 + t304 * t337 + t306 * t338 + t352 * t370 + t354 * t367 + Icges(1,2)) * V_base(5) + (-t303 * t405 + t305 * t337 + t307 * t338 + t353 * t370 + t355 * t367 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t371 * t331 + (t369 * t332 + t333 * t366) * t368 + Icges(2,3) + Icges(1,3)) * V_base(6) + (t302 * t371 + (t304 * t369 + t306 * t366) * t368 + t419) * V_base(5) + (t303 * t371 + (t305 * t369 + t307 * t366) * t368 + t420) * V_base(4)) * V_base(6) / 0.2e1;
T  = t1;
