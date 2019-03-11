% Calculate kinetic energy for
% S6RRRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 20:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR15_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRR15_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR15_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR15_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR15_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:24:26
% EndTime: 2019-03-09 20:24:31
% DurationCPUTime: 4.75s
% Computational Cost: add. (4426->401), mult. (11404->586), div. (0->0), fcn. (14382->14), ass. (0->182)
t417 = Icges(4,1) + Icges(5,2);
t416 = Icges(5,1) + Icges(4,3);
t415 = -Icges(4,4) - Icges(5,6);
t414 = -Icges(5,4) + Icges(4,5);
t413 = Icges(5,5) - Icges(4,6);
t412 = Icges(4,2) + Icges(5,3);
t355 = cos(pkin(6));
t362 = cos(qJ(2));
t363 = cos(qJ(1));
t387 = t362 * t363;
t359 = sin(qJ(2));
t360 = sin(qJ(1));
t389 = t360 * t359;
t323 = t355 * t387 - t389;
t388 = t360 * t362;
t390 = t359 * t363;
t324 = t355 * t390 + t388;
t358 = sin(qJ(3));
t354 = sin(pkin(6));
t395 = sin(pkin(7));
t399 = cos(qJ(3));
t377 = t399 * t395;
t376 = t354 * t377;
t396 = cos(pkin(7));
t378 = t396 * t399;
t282 = -t323 * t378 + t324 * t358 + t363 * t376;
t379 = t358 * t395;
t380 = t358 * t396;
t391 = t354 * t363;
t283 = t323 * t380 + t324 * t399 - t379 * t391;
t382 = t354 * t396;
t307 = -t323 * t395 - t363 * t382;
t411 = t282 * t412 + t283 * t415 + t307 * t413;
t325 = -t355 * t388 - t390;
t326 = -t355 * t389 + t387;
t284 = -t325 * t378 + t326 * t358 - t360 * t376;
t381 = t354 * t395;
t285 = t326 * t399 + (t325 * t396 + t360 * t381) * t358;
t308 = -t325 * t395 + t360 * t382;
t410 = t284 * t412 + t285 * t415 + t308 * t413;
t409 = t282 * t413 + t283 * t414 + t307 * t416;
t408 = t284 * t413 + t285 * t414 + t308 * t416;
t407 = t415 * t282 + t283 * t417 + t414 * t307;
t406 = t415 * t284 + t285 * t417 + t414 * t308;
t393 = t354 * t359;
t305 = -t354 * t362 * t378 - t355 * t377 + t358 * t393;
t306 = t355 * t379 + (t359 * t399 + t362 * t380) * t354;
t322 = t355 * t396 - t362 * t381;
t405 = t305 * t412 + t306 * t415 + t322 * t413;
t404 = t305 * t413 + t306 * t414 + t322 * t416;
t403 = t415 * t305 + t306 * t417 + t414 * t322;
t398 = cos(qJ(5));
t397 = pkin(9) * t355;
t394 = Icges(2,4) * t360;
t392 = t354 * t360;
t386 = qJD(2) * t354;
t385 = V_base(5) * pkin(8) + V_base(1);
t337 = t360 * t386 + V_base(4);
t351 = V_base(6) + qJD(1);
t298 = qJD(3) * t308 + t337;
t338 = qJD(2) * t355 + t351;
t252 = qJD(5) * t285 + t298;
t309 = qJD(3) * t322 + t338;
t336 = -t363 * t386 + V_base(5);
t328 = t360 * pkin(1) - pkin(9) * t391;
t375 = -t328 * t351 + V_base(5) * t397 + t385;
t271 = qJD(5) * t306 + t309;
t329 = pkin(1) * t363 + pkin(9) * t392;
t374 = V_base(4) * t328 - t329 * V_base(5) + V_base(3);
t297 = qJD(3) * t307 + t336;
t251 = qJD(5) * t283 + t297;
t373 = t351 * t329 + V_base(2) + (-pkin(8) - t397) * V_base(4);
t287 = t324 * pkin(2) + pkin(10) * t307;
t313 = pkin(2) * t393 + pkin(10) * t322;
t372 = -t287 * t338 + t336 * t313 + t375;
t288 = t326 * pkin(2) + pkin(10) * t308;
t371 = t337 * t287 - t288 * t336 + t374;
t270 = pkin(3) * t306 + qJ(4) * t305;
t370 = qJD(4) * t284 + t297 * t270 + t372;
t248 = pkin(3) * t283 + qJ(4) * t282;
t369 = qJD(4) * t305 + t298 * t248 + t371;
t368 = t338 * t288 - t313 * t337 + t373;
t249 = pkin(3) * t285 + qJ(4) * t284;
t367 = qJD(4) * t282 + t309 * t249 + t368;
t260 = pkin(4) * t307 + pkin(11) * t283;
t286 = pkin(4) * t322 + pkin(11) * t306;
t366 = t297 * t286 + (-t248 - t260) * t309 + t370;
t261 = pkin(4) * t308 + pkin(11) * t285;
t365 = t298 * t260 + (-t249 - t261) * t297 + t369;
t364 = t309 * t261 + (-t270 - t286) * t298 + t367;
t361 = cos(qJ(6));
t357 = sin(qJ(5));
t356 = sin(qJ(6));
t352 = Icges(2,4) * t363;
t346 = rSges(2,1) * t363 - t360 * rSges(2,2);
t345 = t360 * rSges(2,1) + rSges(2,2) * t363;
t344 = Icges(2,1) * t363 - t394;
t343 = Icges(2,1) * t360 + t352;
t342 = -Icges(2,2) * t360 + t352;
t341 = Icges(2,2) * t363 + t394;
t335 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t334 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t333 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t319 = rSges(3,3) * t355 + (rSges(3,1) * t359 + rSges(3,2) * t362) * t354;
t318 = Icges(3,5) * t355 + (Icges(3,1) * t359 + Icges(3,4) * t362) * t354;
t317 = Icges(3,6) * t355 + (Icges(3,4) * t359 + Icges(3,2) * t362) * t354;
t316 = Icges(3,3) * t355 + (Icges(3,5) * t359 + Icges(3,6) * t362) * t354;
t312 = V_base(5) * rSges(2,3) - t345 * t351 + t385;
t311 = t346 * t351 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t310 = t345 * V_base(4) - t346 * V_base(5) + V_base(3);
t296 = rSges(3,1) * t326 + rSges(3,2) * t325 + rSges(3,3) * t392;
t295 = t324 * rSges(3,1) + t323 * rSges(3,2) - rSges(3,3) * t391;
t294 = Icges(3,1) * t326 + Icges(3,4) * t325 + Icges(3,5) * t392;
t293 = Icges(3,1) * t324 + Icges(3,4) * t323 - Icges(3,5) * t391;
t292 = Icges(3,4) * t326 + Icges(3,2) * t325 + Icges(3,6) * t392;
t291 = Icges(3,4) * t324 + Icges(3,2) * t323 - Icges(3,6) * t391;
t290 = Icges(3,5) * t326 + Icges(3,6) * t325 + Icges(3,3) * t392;
t289 = Icges(3,5) * t324 + Icges(3,6) * t323 - Icges(3,3) * t391;
t281 = t305 * t357 + t322 * t398;
t280 = -t305 * t398 + t322 * t357;
t269 = rSges(5,1) * t322 - rSges(5,2) * t306 + rSges(5,3) * t305;
t268 = rSges(4,1) * t306 - rSges(4,2) * t305 + rSges(4,3) * t322;
t259 = t284 * t357 + t308 * t398;
t258 = -t284 * t398 + t308 * t357;
t257 = t282 * t357 + t307 * t398;
t256 = -t282 * t398 + t307 * t357;
t255 = t281 * t361 + t306 * t356;
t254 = -t281 * t356 + t306 * t361;
t247 = pkin(5) * t281 + pkin(12) * t280;
t245 = -t295 * t338 + t319 * t336 + t375;
t244 = t296 * t338 - t319 * t337 + t373;
t243 = qJD(6) * t280 + t271;
t240 = rSges(5,1) * t308 - rSges(5,2) * t285 + rSges(5,3) * t284;
t239 = rSges(5,1) * t307 - rSges(5,2) * t283 + rSges(5,3) * t282;
t238 = rSges(4,1) * t285 - rSges(4,2) * t284 + rSges(4,3) * t308;
t237 = rSges(4,1) * t283 - rSges(4,2) * t282 + rSges(4,3) * t307;
t236 = t295 * t337 - t296 * t336 + t374;
t223 = rSges(6,1) * t281 - rSges(6,2) * t280 + rSges(6,3) * t306;
t222 = Icges(6,1) * t281 - Icges(6,4) * t280 + Icges(6,5) * t306;
t221 = Icges(6,4) * t281 - Icges(6,2) * t280 + Icges(6,6) * t306;
t220 = Icges(6,5) * t281 - Icges(6,6) * t280 + Icges(6,3) * t306;
t219 = t259 * t361 + t285 * t356;
t218 = -t259 * t356 + t285 * t361;
t217 = t257 * t361 + t283 * t356;
t216 = -t257 * t356 + t283 * t361;
t214 = pkin(5) * t259 + pkin(12) * t258;
t213 = pkin(5) * t257 + pkin(12) * t256;
t212 = qJD(6) * t258 + t252;
t211 = qJD(6) * t256 + t251;
t210 = rSges(6,1) * t259 - rSges(6,2) * t258 + rSges(6,3) * t285;
t209 = rSges(6,1) * t257 - rSges(6,2) * t256 + rSges(6,3) * t283;
t208 = Icges(6,1) * t259 - Icges(6,4) * t258 + Icges(6,5) * t285;
t207 = Icges(6,1) * t257 - Icges(6,4) * t256 + Icges(6,5) * t283;
t206 = Icges(6,4) * t259 - Icges(6,2) * t258 + Icges(6,6) * t285;
t205 = Icges(6,4) * t257 - Icges(6,2) * t256 + Icges(6,6) * t283;
t204 = Icges(6,5) * t259 - Icges(6,6) * t258 + Icges(6,3) * t285;
t203 = Icges(6,5) * t257 - Icges(6,6) * t256 + Icges(6,3) * t283;
t202 = rSges(7,1) * t255 + rSges(7,2) * t254 + rSges(7,3) * t280;
t201 = Icges(7,1) * t255 + Icges(7,4) * t254 + Icges(7,5) * t280;
t200 = Icges(7,4) * t255 + Icges(7,2) * t254 + Icges(7,6) * t280;
t199 = Icges(7,5) * t255 + Icges(7,6) * t254 + Icges(7,3) * t280;
t198 = rSges(7,1) * t219 + rSges(7,2) * t218 + rSges(7,3) * t258;
t197 = rSges(7,1) * t217 + rSges(7,2) * t216 + rSges(7,3) * t256;
t196 = Icges(7,1) * t219 + Icges(7,4) * t218 + Icges(7,5) * t258;
t195 = Icges(7,1) * t217 + Icges(7,4) * t216 + Icges(7,5) * t256;
t194 = Icges(7,4) * t219 + Icges(7,2) * t218 + Icges(7,6) * t258;
t193 = Icges(7,4) * t217 + Icges(7,2) * t216 + Icges(7,6) * t256;
t192 = Icges(7,5) * t219 + Icges(7,6) * t218 + Icges(7,3) * t258;
t191 = Icges(7,5) * t217 + Icges(7,6) * t216 + Icges(7,3) * t256;
t190 = -t237 * t309 + t268 * t297 + t372;
t189 = t238 * t309 - t268 * t298 + t368;
t188 = t237 * t298 - t238 * t297 + t371;
t187 = t269 * t297 + (-t239 - t248) * t309 + t370;
t186 = t240 * t309 + (-t269 - t270) * t298 + t367;
t185 = t239 * t298 + (-t240 - t249) * t297 + t369;
t184 = -t209 * t271 + t223 * t251 + t366;
t183 = t210 * t271 - t223 * t252 + t364;
t182 = t209 * t252 - t210 * t251 + t365;
t181 = -t197 * t243 + t202 * t211 - t213 * t271 + t247 * t251 + t366;
t180 = t198 * t243 - t202 * t212 + t214 * t271 - t247 * t252 + t364;
t179 = t197 * t212 - t198 * t211 + t213 * t252 - t214 * t251 + t365;
t1 = t338 * ((t289 * t336 + t290 * t337 + t316 * t338) * t355 + ((t292 * t362 + t294 * t359) * t337 + (t291 * t362 + t293 * t359) * t336 + (t317 * t362 + t318 * t359) * t338) * t354) / 0.2e1 + m(1) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + m(2) * (t310 ^ 2 + t311 ^ 2 + t312 ^ 2) / 0.2e1 + t271 * ((t204 * t306 - t206 * t280 + t208 * t281) * t252 + (t203 * t306 - t205 * t280 + t207 * t281) * t251 + (t220 * t306 - t221 * t280 + t222 * t281) * t271) / 0.2e1 + t252 * ((t204 * t285 - t206 * t258 + t208 * t259) * t252 + (t203 * t285 - t205 * t258 + t207 * t259) * t251 + (t220 * t285 - t221 * t258 + t222 * t259) * t271) / 0.2e1 + t251 * ((t204 * t283 - t206 * t256 + t208 * t257) * t252 + (t203 * t283 - t205 * t256 + t207 * t257) * t251 + (t220 * t283 - t221 * t256 + t222 * t257) * t271) / 0.2e1 + t243 * ((t192 * t280 + t194 * t254 + t196 * t255) * t212 + (t191 * t280 + t193 * t254 + t195 * t255) * t211 + (t199 * t280 + t200 * t254 + t201 * t255) * t243) / 0.2e1 + t212 * ((t192 * t258 + t218 * t194 + t219 * t196) * t212 + (t191 * t258 + t193 * t218 + t195 * t219) * t211 + (t199 * t258 + t200 * t218 + t201 * t219) * t243) / 0.2e1 + t211 * ((t192 * t256 + t194 * t216 + t196 * t217) * t212 + (t191 * t256 + t216 * t193 + t217 * t195) * t211 + (t199 * t256 + t200 * t216 + t201 * t217) * t243) / 0.2e1 + m(3) * (t236 ^ 2 + t244 ^ 2 + t245 ^ 2) / 0.2e1 + t336 * ((-t290 * t391 + t323 * t292 + t324 * t294) * t337 + (-t289 * t391 + t323 * t291 + t324 * t293) * t336 + (-t316 * t391 + t323 * t317 + t324 * t318) * t338) / 0.2e1 + t337 * ((t290 * t392 + t292 * t325 + t294 * t326) * t337 + (t289 * t392 + t291 * t325 + t293 * t326) * t336 + (t316 * t392 + t317 * t325 + t318 * t326) * t338) / 0.2e1 + m(4) * (t188 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + m(5) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(6) * (t182 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + m(7) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + ((t282 * t405 + t283 * t403 + t307 * t404) * t309 + (t282 * t410 + t283 * t406 + t307 * t408) * t298 + (t282 * t411 + t407 * t283 + t409 * t307) * t297) * t297 / 0.2e1 + ((t284 * t405 + t285 * t403 + t308 * t404) * t309 + (t284 * t410 + t285 * t406 + t308 * t408) * t298 + (t284 * t411 + t407 * t285 + t409 * t308) * t297) * t298 / 0.2e1 + ((t305 * t405 + t306 * t403 + t322 * t404) * t309 + (t305 * t410 + t306 * t406 + t322 * t408) * t298 + (t305 * t411 + t407 * t306 + t409 * t322) * t297) * t309 / 0.2e1 + ((-t360 * t341 + t343 * t363 + Icges(1,4)) * V_base(5) + (-t360 * t342 + t344 * t363 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t341 * t363 + t360 * t343 + Icges(1,2)) * V_base(5) + (t342 * t363 + t360 * t344 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t360 + Icges(2,6) * t363) * V_base(5) + (Icges(2,5) * t363 - Icges(2,6) * t360) * V_base(4) + Icges(2,3) * t351 / 0.2e1) * t351;
T  = t1;
