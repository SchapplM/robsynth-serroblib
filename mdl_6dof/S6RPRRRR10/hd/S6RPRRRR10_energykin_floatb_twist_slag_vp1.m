% Calculate kinetic energy for
% S6RPRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR10_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR10_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR10_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR10_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR10_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR10_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:27:35
% EndTime: 2019-03-09 07:27:40
% DurationCPUTime: 4.67s
% Computational Cost: add. (5518->452), mult. (13137->664), div. (0->0), fcn. (16739->16), ass. (0->197)
t416 = cos(qJ(3));
t375 = cos(qJ(4));
t415 = pkin(4) * t375;
t413 = cos(pkin(7));
t412 = sin(pkin(7));
t373 = sin(qJ(1));
t411 = Icges(2,4) * t373;
t369 = cos(pkin(6));
t410 = qJ(2) * t369;
t366 = sin(pkin(13));
t368 = cos(pkin(13));
t376 = cos(qJ(1));
t402 = t369 * t376;
t337 = -t366 * t373 + t368 * t402;
t367 = sin(pkin(6));
t394 = t367 * t413;
t319 = -t337 * t412 - t376 * t394;
t371 = sin(qJ(4));
t409 = t319 * t371;
t403 = t369 * t373;
t339 = -t366 * t376 - t368 * t403;
t320 = -t339 * t412 + t373 * t394;
t408 = t320 * t371;
t393 = t367 * t412;
t336 = -t368 * t393 + t369 * t413;
t407 = t336 * t371;
t406 = t366 * t367;
t405 = t367 * t373;
t404 = t367 * t376;
t401 = qJ(4) + qJ(5);
t400 = qJD(2) * t367;
t399 = V_base(5) * pkin(8) + V_base(1);
t310 = qJD(3) * t320 + V_base(4);
t309 = qJD(3) * t319 + V_base(5);
t396 = cos(t401);
t362 = V_base(6) + qJD(1);
t395 = -pkin(8) - t410;
t342 = pkin(1) * t373 - qJ(2) * t404;
t392 = qJD(2) * t369 + V_base(4) * t342 + V_base(3);
t340 = -t366 * t403 + t368 * t376;
t372 = sin(qJ(3));
t390 = t416 * t412;
t388 = t367 * t390;
t391 = t413 * t416;
t296 = -t339 * t391 + t340 * t372 - t373 * t388;
t273 = qJD(4) * t296 + t310;
t338 = t366 * t402 + t368 * t373;
t294 = -t337 * t391 + t338 * t372 + t376 * t388;
t272 = qJD(4) * t294 + t309;
t326 = qJD(3) * t336 + t362;
t389 = t373 * t400 + V_base(5) * t410 + t399;
t248 = qJD(5) * t296 + t273;
t247 = qJD(5) * t294 + t272;
t317 = -t367 * t368 * t391 - t369 * t390 + t372 * t406;
t286 = qJD(4) * t317 + t326;
t278 = qJD(5) * t317 + t286;
t343 = pkin(1) * t376 + qJ(2) * t405;
t387 = t362 * t343 - t376 * t400 + V_base(2);
t299 = t338 * pkin(2) + pkin(9) * t319;
t323 = pkin(2) * t406 + pkin(9) * t336;
t386 = V_base(5) * t323 + (-t299 - t342) * t362 + t389;
t300 = t340 * pkin(2) + pkin(9) * t320;
t385 = V_base(4) * t299 + (-t300 - t343) * V_base(5) + t392;
t295 = t338 * t416 + (t337 * t413 - t376 * t393) * t372;
t263 = t295 * pkin(3) + t294 * pkin(10);
t318 = t369 * t412 * t372 + (t368 * t372 * t413 + t366 * t416) * t367;
t283 = t318 * pkin(3) + t317 * pkin(10);
t384 = -t263 * t326 + t309 * t283 + t386;
t297 = t340 * t416 + (t339 * t413 + t373 * t393) * t372;
t264 = t297 * pkin(3) + t296 * pkin(10);
t383 = t310 * t263 - t264 * t309 + t385;
t382 = t362 * t300 + (-t323 + t395) * V_base(4) + t387;
t206 = pkin(4) * t409 + pkin(11) * t294 + t295 * t415;
t233 = pkin(4) * t407 + pkin(11) * t317 + t318 * t415;
t381 = -t206 * t286 + t272 * t233 + t384;
t207 = pkin(4) * t408 + pkin(11) * t296 + t297 * t415;
t380 = t273 * t206 - t207 * t272 + t383;
t379 = t326 * t264 - t283 * t310 + t382;
t378 = t286 * t207 - t233 * t273 + t379;
t374 = cos(qJ(6));
t370 = sin(qJ(6));
t364 = Icges(2,4) * t376;
t363 = sin(t401);
t356 = rSges(2,1) * t376 - rSges(2,2) * t373;
t355 = rSges(2,1) * t373 + rSges(2,2) * t376;
t354 = Icges(2,1) * t376 - t411;
t353 = Icges(2,1) * t373 + t364;
t352 = -Icges(2,2) * t373 + t364;
t351 = Icges(2,2) * t376 + t411;
t350 = Icges(2,5) * t376 - Icges(2,6) * t373;
t349 = Icges(2,5) * t373 + Icges(2,6) * t376;
t348 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t347 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t346 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t333 = rSges(3,3) * t369 + (rSges(3,1) * t366 + rSges(3,2) * t368) * t367;
t332 = Icges(3,5) * t369 + (Icges(3,1) * t366 + Icges(3,4) * t368) * t367;
t331 = Icges(3,6) * t369 + (Icges(3,4) * t366 + Icges(3,2) * t368) * t367;
t330 = Icges(3,3) * t369 + (Icges(3,5) * t366 + Icges(3,6) * t368) * t367;
t325 = V_base(5) * rSges(2,3) - t355 * t362 + t399;
t324 = t356 * t362 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t322 = t355 * V_base(4) - t356 * V_base(5) + V_base(3);
t308 = rSges(3,1) * t340 + rSges(3,2) * t339 + rSges(3,3) * t405;
t307 = rSges(3,1) * t338 + rSges(3,2) * t337 - rSges(3,3) * t404;
t306 = Icges(3,1) * t340 + Icges(3,4) * t339 + Icges(3,5) * t405;
t305 = Icges(3,1) * t338 + Icges(3,4) * t337 - Icges(3,5) * t404;
t304 = Icges(3,4) * t340 + Icges(3,2) * t339 + Icges(3,6) * t405;
t303 = Icges(3,4) * t338 + Icges(3,2) * t337 - Icges(3,6) * t404;
t302 = Icges(3,5) * t340 + Icges(3,6) * t339 + Icges(3,3) * t405;
t301 = Icges(3,5) * t338 + Icges(3,6) * t337 - Icges(3,3) * t404;
t292 = t318 * t375 + t407;
t291 = -t318 * t371 + t336 * t375;
t285 = t318 * t396 + t336 * t363;
t284 = t318 * t363 - t336 * t396;
t282 = rSges(4,1) * t318 - rSges(4,2) * t317 + rSges(4,3) * t336;
t281 = Icges(4,1) * t318 - Icges(4,4) * t317 + Icges(4,5) * t336;
t280 = Icges(4,4) * t318 - Icges(4,2) * t317 + Icges(4,6) * t336;
t279 = Icges(4,5) * t318 - Icges(4,6) * t317 + Icges(4,3) * t336;
t277 = t297 * t375 + t408;
t276 = -t297 * t371 + t320 * t375;
t275 = t295 * t375 + t409;
t274 = -t295 * t371 + t319 * t375;
t271 = t297 * t396 + t320 * t363;
t270 = t297 * t363 - t320 * t396;
t269 = t295 * t396 + t319 * t363;
t268 = t295 * t363 - t319 * t396;
t267 = t285 * t374 + t317 * t370;
t266 = -t285 * t370 + t317 * t374;
t262 = t333 * V_base(5) + (-t307 - t342) * t362 + t389;
t261 = t308 * t362 + (-t333 + t395) * V_base(4) + t387;
t260 = pkin(5) * t285 + pkin(12) * t284;
t259 = t307 * V_base(4) + (-t308 - t343) * V_base(5) + t392;
t257 = rSges(4,1) * t297 - rSges(4,2) * t296 + rSges(4,3) * t320;
t256 = rSges(4,1) * t295 - rSges(4,2) * t294 + rSges(4,3) * t319;
t255 = Icges(4,1) * t297 - Icges(4,4) * t296 + Icges(4,5) * t320;
t254 = Icges(4,1) * t295 - Icges(4,4) * t294 + Icges(4,5) * t319;
t253 = Icges(4,4) * t297 - Icges(4,2) * t296 + Icges(4,6) * t320;
t252 = Icges(4,4) * t295 - Icges(4,2) * t294 + Icges(4,6) * t319;
t251 = Icges(4,5) * t297 - Icges(4,6) * t296 + Icges(4,3) * t320;
t250 = Icges(4,5) * t295 - Icges(4,6) * t294 + Icges(4,3) * t319;
t249 = rSges(5,1) * t292 + rSges(5,2) * t291 + rSges(5,3) * t317;
t245 = Icges(5,1) * t292 + Icges(5,4) * t291 + Icges(5,5) * t317;
t244 = Icges(5,4) * t292 + Icges(5,2) * t291 + Icges(5,6) * t317;
t243 = Icges(5,5) * t292 + Icges(5,6) * t291 + Icges(5,3) * t317;
t242 = qJD(6) * t284 + t278;
t241 = rSges(6,1) * t285 - rSges(6,2) * t284 + rSges(6,3) * t317;
t240 = Icges(6,1) * t285 - Icges(6,4) * t284 + Icges(6,5) * t317;
t239 = Icges(6,4) * t285 - Icges(6,2) * t284 + Icges(6,6) * t317;
t238 = Icges(6,5) * t285 - Icges(6,6) * t284 + Icges(6,3) * t317;
t237 = t271 * t374 + t296 * t370;
t236 = -t271 * t370 + t296 * t374;
t235 = t269 * t374 + t294 * t370;
t234 = -t269 * t370 + t294 * t374;
t232 = pkin(5) * t271 + pkin(12) * t270;
t231 = pkin(5) * t269 + pkin(12) * t268;
t230 = rSges(5,1) * t277 + rSges(5,2) * t276 + rSges(5,3) * t296;
t229 = rSges(5,1) * t275 + rSges(5,2) * t274 + rSges(5,3) * t294;
t228 = Icges(5,1) * t277 + Icges(5,4) * t276 + Icges(5,5) * t296;
t227 = Icges(5,1) * t275 + Icges(5,4) * t274 + Icges(5,5) * t294;
t226 = Icges(5,4) * t277 + Icges(5,2) * t276 + Icges(5,6) * t296;
t225 = Icges(5,4) * t275 + Icges(5,2) * t274 + Icges(5,6) * t294;
t224 = Icges(5,5) * t277 + Icges(5,6) * t276 + Icges(5,3) * t296;
t223 = Icges(5,5) * t275 + Icges(5,6) * t274 + Icges(5,3) * t294;
t222 = qJD(6) * t270 + t248;
t221 = qJD(6) * t268 + t247;
t220 = rSges(6,1) * t271 - rSges(6,2) * t270 + rSges(6,3) * t296;
t219 = rSges(6,1) * t269 - rSges(6,2) * t268 + rSges(6,3) * t294;
t218 = Icges(6,1) * t271 - Icges(6,4) * t270 + Icges(6,5) * t296;
t217 = Icges(6,1) * t269 - Icges(6,4) * t268 + Icges(6,5) * t294;
t216 = Icges(6,4) * t271 - Icges(6,2) * t270 + Icges(6,6) * t296;
t215 = Icges(6,4) * t269 - Icges(6,2) * t268 + Icges(6,6) * t294;
t214 = Icges(6,5) * t271 - Icges(6,6) * t270 + Icges(6,3) * t296;
t213 = Icges(6,5) * t269 - Icges(6,6) * t268 + Icges(6,3) * t294;
t211 = rSges(7,1) * t267 + rSges(7,2) * t266 + rSges(7,3) * t284;
t210 = Icges(7,1) * t267 + Icges(7,4) * t266 + Icges(7,5) * t284;
t209 = Icges(7,4) * t267 + Icges(7,2) * t266 + Icges(7,6) * t284;
t208 = Icges(7,5) * t267 + Icges(7,6) * t266 + Icges(7,3) * t284;
t203 = -t256 * t326 + t282 * t309 + t386;
t202 = t257 * t326 - t282 * t310 + t382;
t201 = rSges(7,1) * t237 + rSges(7,2) * t236 + rSges(7,3) * t270;
t200 = rSges(7,1) * t235 + rSges(7,2) * t234 + rSges(7,3) * t268;
t199 = Icges(7,1) * t237 + Icges(7,4) * t236 + Icges(7,5) * t270;
t198 = Icges(7,1) * t235 + Icges(7,4) * t234 + Icges(7,5) * t268;
t197 = Icges(7,4) * t237 + Icges(7,2) * t236 + Icges(7,6) * t270;
t196 = Icges(7,4) * t235 + Icges(7,2) * t234 + Icges(7,6) * t268;
t195 = Icges(7,5) * t237 + Icges(7,6) * t236 + Icges(7,3) * t270;
t194 = Icges(7,5) * t235 + Icges(7,6) * t234 + Icges(7,3) * t268;
t193 = t256 * t310 - t257 * t309 + t385;
t192 = -t229 * t286 + t249 * t272 + t384;
t191 = t230 * t286 - t249 * t273 + t379;
t190 = t229 * t273 - t230 * t272 + t383;
t189 = -t219 * t278 + t241 * t247 + t381;
t188 = t220 * t278 - t241 * t248 + t378;
t187 = t219 * t248 - t220 * t247 + t380;
t186 = -t200 * t242 + t211 * t221 - t231 * t278 + t247 * t260 + t381;
t185 = t201 * t242 - t211 * t222 + t232 * t278 - t248 * t260 + t378;
t184 = t200 * t222 - t201 * t221 + t231 * t248 - t232 * t247 + t380;
t1 = m(1) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + t326 * ((t251 * t336 - t253 * t317 + t255 * t318) * t310 + (t250 * t336 - t252 * t317 + t254 * t318) * t309 + (t279 * t336 - t280 * t317 + t281 * t318) * t326) / 0.2e1 + t309 * ((t251 * t319 - t253 * t294 + t255 * t295) * t310 + (t250 * t319 - t252 * t294 + t254 * t295) * t309 + (t279 * t319 - t280 * t294 + t281 * t295) * t326) / 0.2e1 + t310 * ((t251 * t320 - t253 * t296 + t255 * t297) * t310 + (t250 * t320 - t252 * t296 + t254 * t297) * t309 + (t279 * t320 - t280 * t296 + t281 * t297) * t326) / 0.2e1 + m(2) * (t322 ^ 2 + t324 ^ 2 + t325 ^ 2) / 0.2e1 + t286 * ((t224 * t317 + t226 * t291 + t228 * t292) * t273 + (t223 * t317 + t225 * t291 + t227 * t292) * t272 + (t243 * t317 + t244 * t291 + t245 * t292) * t286) / 0.2e1 + t278 * ((t214 * t317 - t216 * t284 + t218 * t285) * t248 + (t213 * t317 - t215 * t284 + t217 * t285) * t247 + (t238 * t317 - t239 * t284 + t240 * t285) * t278) / 0.2e1 + t248 * ((t214 * t296 - t216 * t270 + t218 * t271) * t248 + (t213 * t296 - t215 * t270 + t217 * t271) * t247 + (t238 * t296 - t239 * t270 + t240 * t271) * t278) / 0.2e1 + t273 * ((t224 * t296 + t226 * t276 + t228 * t277) * t273 + (t223 * t296 + t225 * t276 + t227 * t277) * t272 + (t243 * t296 + t244 * t276 + t245 * t277) * t286) / 0.2e1 + t272 * ((t224 * t294 + t226 * t274 + t228 * t275) * t273 + (t223 * t294 + t225 * t274 + t227 * t275) * t272 + (t243 * t294 + t244 * t274 + t245 * t275) * t286) / 0.2e1 + t247 * ((t214 * t294 - t216 * t268 + t218 * t269) * t248 + (t213 * t294 - t215 * t268 + t217 * t269) * t247 + (t238 * t294 - t239 * t268 + t240 * t269) * t278) / 0.2e1 + t242 * ((t195 * t284 + t197 * t266 + t199 * t267) * t222 + (t194 * t284 + t196 * t266 + t198 * t267) * t221 + (t208 * t284 + t209 * t266 + t210 * t267) * t242) / 0.2e1 + t222 * ((t195 * t270 + t197 * t236 + t199 * t237) * t222 + (t194 * t270 + t196 * t236 + t198 * t237) * t221 + (t208 * t270 + t209 * t236 + t210 * t237) * t242) / 0.2e1 + m(7) * (t184 ^ 2 + t185 ^ 2 + t186 ^ 2) / 0.2e1 + m(6) * (t187 ^ 2 + t188 ^ 2 + t189 ^ 2) / 0.2e1 + m(5) * (t190 ^ 2 + t191 ^ 2 + t192 ^ 2) / 0.2e1 + m(4) * (t193 ^ 2 + t202 ^ 2 + t203 ^ 2) / 0.2e1 + m(3) * (t259 ^ 2 + t261 ^ 2 + t262 ^ 2) / 0.2e1 + t221 * ((t195 * t268 + t197 * t234 + t199 * t235) * t222 + (t194 * t268 + t196 * t234 + t198 * t235) * t221 + (t208 * t268 + t209 * t234 + t210 * t235) * t242) / 0.2e1 + ((t301 * V_base(5) + t302 * V_base(4) + t330 * t362) * t369 + ((t304 * t368 + t306 * t366) * V_base(4) + (t303 * t368 + t305 * t366) * V_base(5) + (t331 * t368 + t332 * t366) * t362) * t367 + Icges(2,3) * t362 + t349 * V_base(5) + t350 * V_base(4)) * t362 / 0.2e1 + ((t330 * t405 + t331 * t339 + t332 * t340 + t350) * t362 + (t301 * t405 + t303 * t339 + t305 * t340 - t351 * t373 + t353 * t376 + Icges(1,4)) * V_base(5) + (t302 * t405 + t304 * t339 + t306 * t340 - t352 * t373 + t354 * t376 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t330 * t404 + t331 * t337 + t332 * t338 + t349) * t362 + (-t301 * t404 + t303 * t337 + t305 * t338 + t351 * t376 + t353 * t373 + Icges(1,2)) * V_base(5) + (-t302 * t404 + t304 * t337 + t306 * t338 + t352 * t376 + t354 * t373 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
