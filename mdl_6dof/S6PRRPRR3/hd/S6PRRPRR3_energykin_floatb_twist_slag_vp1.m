% Calculate kinetic energy for
% S6PRRPRR3
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:02:08
% EndTime: 2019-03-08 22:02:13
% DurationCPUTime: 5.44s
% Computational Cost: add. (6009->454), mult. (15638->663), div. (0->0), fcn. (20097->16), ass. (0->201)
t432 = Icges(4,3) + Icges(5,3);
t370 = sin(pkin(12));
t373 = cos(pkin(12));
t379 = sin(qJ(2));
t375 = cos(pkin(6));
t382 = cos(qJ(2));
t406 = t375 * t382;
t342 = -t370 * t406 - t373 * t379;
t371 = sin(pkin(7));
t372 = sin(pkin(6));
t374 = cos(pkin(7));
t411 = t372 * t374;
t321 = -t342 * t371 + t370 * t411;
t431 = pkin(9) * t321;
t340 = -t370 * t379 + t373 * t406;
t399 = t340 * t371 + t373 * t411;
t430 = pkin(9) * t399;
t378 = sin(qJ(3));
t381 = cos(qJ(3));
t416 = sin(pkin(13));
t417 = cos(pkin(13));
t347 = -t378 * t416 + t381 * t417;
t333 = t347 * t374;
t407 = t375 * t379;
t341 = t370 * t382 + t373 * t407;
t346 = -t378 * t417 - t381 * t416;
t393 = t371 * t347;
t392 = t372 * t393;
t278 = t333 * t340 + t341 * t346 - t373 * t392;
t332 = t346 * t371;
t334 = t346 * t374;
t412 = t372 * t373;
t279 = t332 * t412 - t334 * t340 + t341 * t347;
t398 = t340 * t374 - t371 * t412;
t294 = -t341 * t378 + t381 * t398;
t295 = t341 * t381 + t378 * t398;
t429 = Icges(4,5) * t295 + Icges(5,5) * t279 + Icges(4,6) * t294 + Icges(5,6) * t278 - t432 * t399;
t343 = -t370 * t407 + t373 * t382;
t280 = t342 * t333 + t343 * t346 + t370 * t392;
t414 = t370 * t372;
t281 = -t332 * t414 - t334 * t342 + t343 * t347;
t397 = t342 * t374 + t371 * t414;
t296 = -t343 * t378 + t381 * t397;
t297 = t343 * t381 + t378 * t397;
t428 = Icges(4,5) * t297 + Icges(5,5) * t281 + Icges(4,6) * t296 + Icges(5,6) * t280 + t432 * t321;
t409 = t372 * t382;
t410 = t372 * t379;
t290 = t333 * t409 + t346 * t410 + t375 * t393;
t291 = -t375 * t332 + (-t334 * t382 + t347 * t379) * t372;
t408 = t374 * t382;
t318 = t375 * t371 * t381 + (-t378 * t379 + t381 * t408) * t372;
t413 = t371 * t378;
t319 = t375 * t413 + (t378 * t408 + t379 * t381) * t372;
t339 = -t371 * t409 + t375 * t374;
t427 = Icges(4,5) * t319 + Icges(5,5) * t291 + Icges(4,6) * t318 + Icges(5,6) * t290 + t432 * t339;
t421 = cos(qJ(5));
t420 = pkin(8) * t375;
t419 = pkin(3) * t381;
t418 = pkin(9) + qJ(4);
t415 = Icges(2,4) * t370;
t405 = qJD(2) * t372;
t404 = V_base(5) * qJ(1) + V_base(1);
t400 = qJD(1) + V_base(3);
t354 = t370 * t405 + V_base(4);
t363 = qJD(2) * t375 + V_base(6);
t309 = qJD(3) * t321 + t354;
t322 = qJD(3) * t339 + t363;
t265 = -qJD(5) * t280 + t309;
t275 = -qJD(5) * t290 + t322;
t353 = -t373 * t405 + V_base(5);
t308 = -qJD(3) * t399 + t353;
t348 = pkin(1) * t370 - pkin(8) * t412;
t396 = -t348 * V_base(6) + V_base(5) * t420 + t404;
t349 = pkin(1) * t373 + pkin(8) * t414;
t395 = V_base(4) * t348 - t349 * V_base(5) + t400;
t264 = -qJD(5) * t278 + t308;
t394 = V_base(6) * t349 + V_base(2) + (-qJ(1) - t420) * V_base(4);
t298 = pkin(2) * t341 - t430;
t324 = pkin(2) * t410 + pkin(9) * t339;
t391 = -t298 * t363 + t353 * t324 + t396;
t299 = pkin(2) * t343 + t431;
t390 = t354 * t298 - t299 * t353 + t395;
t389 = t363 * t299 - t324 * t354 + t394;
t337 = pkin(3) * t413 + t374 * t418;
t338 = pkin(3) * t374 * t378 - t371 * t418;
t288 = (-pkin(9) * t374 + t337) * t375 + ((pkin(9) * t371 + t338) * t382 + t419 * t379) * t372;
t388 = qJD(4) * t321 + t308 * t288 + t391;
t266 = -t337 * t412 + t338 * t340 + t341 * t419 + t430;
t387 = qJD(4) * t339 + t309 * t266 + t390;
t267 = t337 * t414 + t338 * t342 + t343 * t419 - t431;
t386 = -qJD(4) * t399 + t322 * t267 + t389;
t242 = pkin(4) * t279 - pkin(10) * t278;
t263 = pkin(4) * t291 - pkin(10) * t290;
t385 = t308 * t263 + (-t242 - t266) * t322 + t388;
t243 = pkin(4) * t281 - pkin(10) * t280;
t384 = t309 * t242 + (-t243 - t267) * t308 + t387;
t383 = t322 * t243 + (-t263 - t288) * t309 + t386;
t380 = cos(qJ(6));
t377 = sin(qJ(5));
t376 = sin(qJ(6));
t368 = Icges(2,4) * t373;
t362 = rSges(2,1) * t373 - rSges(2,2) * t370;
t361 = rSges(2,1) * t370 + rSges(2,2) * t373;
t360 = Icges(2,1) * t373 - t415;
t359 = Icges(2,1) * t370 + t368;
t358 = -Icges(2,2) * t370 + t368;
t357 = Icges(2,2) * t373 + t415;
t352 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t351 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t350 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t330 = t375 * rSges(3,3) + (rSges(3,1) * t379 + rSges(3,2) * t382) * t372;
t329 = Icges(3,5) * t375 + (Icges(3,1) * t379 + Icges(3,4) * t382) * t372;
t328 = Icges(3,6) * t375 + (Icges(3,4) * t379 + Icges(3,2) * t382) * t372;
t327 = Icges(3,3) * t375 + (Icges(3,5) * t379 + Icges(3,6) * t382) * t372;
t326 = V_base(5) * rSges(2,3) - t361 * V_base(6) + t404;
t325 = t362 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t317 = t361 * V_base(4) - t362 * V_base(5) + t400;
t307 = rSges(3,1) * t343 + rSges(3,2) * t342 + rSges(3,3) * t414;
t306 = rSges(3,1) * t341 + rSges(3,2) * t340 - rSges(3,3) * t412;
t305 = Icges(3,1) * t343 + Icges(3,4) * t342 + Icges(3,5) * t414;
t304 = Icges(3,1) * t341 + Icges(3,4) * t340 - Icges(3,5) * t412;
t303 = Icges(3,4) * t343 + Icges(3,2) * t342 + Icges(3,6) * t414;
t302 = Icges(3,4) * t341 + Icges(3,2) * t340 - Icges(3,6) * t412;
t301 = Icges(3,5) * t343 + Icges(3,6) * t342 + Icges(3,3) * t414;
t300 = Icges(3,5) * t341 + Icges(3,6) * t340 - Icges(3,3) * t412;
t287 = rSges(4,1) * t319 + rSges(4,2) * t318 + rSges(4,3) * t339;
t286 = Icges(4,1) * t319 + Icges(4,4) * t318 + Icges(4,5) * t339;
t285 = Icges(4,4) * t319 + Icges(4,2) * t318 + Icges(4,6) * t339;
t283 = t291 * t421 + t339 * t377;
t282 = t291 * t377 - t339 * t421;
t274 = -t306 * t363 + t330 * t353 + t396;
t273 = t307 * t363 - t330 * t354 + t394;
t271 = t281 * t421 + t321 * t377;
t270 = t281 * t377 - t321 * t421;
t269 = t279 * t421 - t377 * t399;
t268 = t279 * t377 + t399 * t421;
t262 = rSges(4,1) * t297 + rSges(4,2) * t296 + rSges(4,3) * t321;
t261 = rSges(4,1) * t295 + rSges(4,2) * t294 - rSges(4,3) * t399;
t260 = Icges(4,1) * t297 + Icges(4,4) * t296 + Icges(4,5) * t321;
t259 = Icges(4,1) * t295 + Icges(4,4) * t294 - Icges(4,5) * t399;
t258 = Icges(4,4) * t297 + Icges(4,2) * t296 + Icges(4,6) * t321;
t257 = Icges(4,4) * t295 + Icges(4,2) * t294 - Icges(4,6) * t399;
t254 = t306 * t354 - t307 * t353 + t395;
t252 = rSges(5,1) * t291 + rSges(5,2) * t290 + rSges(5,3) * t339;
t251 = Icges(5,1) * t291 + Icges(5,4) * t290 + Icges(5,5) * t339;
t250 = Icges(5,4) * t291 + Icges(5,2) * t290 + Icges(5,6) * t339;
t248 = t283 * t380 - t290 * t376;
t247 = -t283 * t376 - t290 * t380;
t244 = pkin(5) * t283 + pkin(11) * t282;
t241 = qJD(6) * t282 + t275;
t239 = rSges(5,1) * t281 + rSges(5,2) * t280 + rSges(5,3) * t321;
t238 = rSges(5,1) * t279 + rSges(5,2) * t278 - rSges(5,3) * t399;
t237 = Icges(5,1) * t281 + Icges(5,4) * t280 + Icges(5,5) * t321;
t236 = Icges(5,1) * t279 + Icges(5,4) * t278 - Icges(5,5) * t399;
t235 = Icges(5,4) * t281 + Icges(5,2) * t280 + Icges(5,6) * t321;
t234 = Icges(5,4) * t279 + Icges(5,2) * t278 - Icges(5,6) * t399;
t231 = t271 * t380 - t280 * t376;
t230 = -t271 * t376 - t280 * t380;
t229 = t269 * t380 - t278 * t376;
t228 = -t269 * t376 - t278 * t380;
t226 = rSges(6,1) * t283 - rSges(6,2) * t282 - rSges(6,3) * t290;
t225 = Icges(6,1) * t283 - Icges(6,4) * t282 - Icges(6,5) * t290;
t224 = Icges(6,4) * t283 - Icges(6,2) * t282 - Icges(6,6) * t290;
t223 = Icges(6,5) * t283 - Icges(6,6) * t282 - Icges(6,3) * t290;
t222 = pkin(5) * t271 + pkin(11) * t270;
t221 = pkin(5) * t269 + pkin(11) * t268;
t220 = qJD(6) * t270 + t265;
t219 = qJD(6) * t268 + t264;
t218 = rSges(6,1) * t271 - rSges(6,2) * t270 - rSges(6,3) * t280;
t217 = rSges(6,1) * t269 - rSges(6,2) * t268 - rSges(6,3) * t278;
t216 = Icges(6,1) * t271 - Icges(6,4) * t270 - Icges(6,5) * t280;
t215 = Icges(6,1) * t269 - Icges(6,4) * t268 - Icges(6,5) * t278;
t214 = Icges(6,4) * t271 - Icges(6,2) * t270 - Icges(6,6) * t280;
t213 = Icges(6,4) * t269 - Icges(6,2) * t268 - Icges(6,6) * t278;
t212 = Icges(6,5) * t271 - Icges(6,6) * t270 - Icges(6,3) * t280;
t211 = Icges(6,5) * t269 - Icges(6,6) * t268 - Icges(6,3) * t278;
t210 = -t261 * t322 + t287 * t308 + t391;
t209 = t262 * t322 - t287 * t309 + t389;
t208 = rSges(7,1) * t248 + rSges(7,2) * t247 + rSges(7,3) * t282;
t207 = Icges(7,1) * t248 + Icges(7,4) * t247 + Icges(7,5) * t282;
t206 = Icges(7,4) * t248 + Icges(7,2) * t247 + Icges(7,6) * t282;
t205 = Icges(7,5) * t248 + Icges(7,6) * t247 + Icges(7,3) * t282;
t204 = t261 * t309 - t262 * t308 + t390;
t203 = rSges(7,1) * t231 + rSges(7,2) * t230 + rSges(7,3) * t270;
t202 = rSges(7,1) * t229 + rSges(7,2) * t228 + rSges(7,3) * t268;
t201 = Icges(7,1) * t231 + Icges(7,4) * t230 + Icges(7,5) * t270;
t200 = Icges(7,1) * t229 + Icges(7,4) * t228 + Icges(7,5) * t268;
t199 = Icges(7,4) * t231 + Icges(7,2) * t230 + Icges(7,6) * t270;
t198 = Icges(7,4) * t229 + Icges(7,2) * t228 + Icges(7,6) * t268;
t197 = Icges(7,5) * t231 + Icges(7,6) * t230 + Icges(7,3) * t270;
t196 = Icges(7,5) * t229 + Icges(7,6) * t228 + Icges(7,3) * t268;
t195 = t252 * t308 + (-t238 - t266) * t322 + t388;
t194 = t239 * t322 + (-t252 - t288) * t309 + t386;
t193 = t238 * t309 + (-t239 - t267) * t308 + t387;
t192 = -t217 * t275 + t226 * t264 + t385;
t191 = t218 * t275 - t226 * t265 + t383;
t190 = t217 * t265 - t218 * t264 + t384;
t189 = -t202 * t241 + t208 * t219 - t221 * t275 + t244 * t264 + t385;
t188 = t203 * t241 - t208 * t220 + t222 * t275 - t244 * t265 + t383;
t187 = t202 * t220 - t203 * t219 + t221 * t265 - t222 * t264 + t384;
t1 = m(6) * (t190 ^ 2 + t191 ^ 2 + t192 ^ 2) / 0.2e1 + m(7) * (t187 ^ 2 + t188 ^ 2 + t189 ^ 2) / 0.2e1 + m(4) * (t204 ^ 2 + t209 ^ 2 + t210 ^ 2) / 0.2e1 + m(5) * (t193 ^ 2 + t194 ^ 2 + t195 ^ 2) / 0.2e1 + m(1) * (t350 ^ 2 + t351 ^ 2 + t352 ^ 2) / 0.2e1 + t264 * ((-t212 * t278 - t214 * t268 + t216 * t269) * t265 + (-t211 * t278 - t213 * t268 + t215 * t269) * t264 + (-t223 * t278 - t224 * t268 + t225 * t269) * t275) / 0.2e1 + t265 * ((-t212 * t280 - t214 * t270 + t216 * t271) * t265 + (-t211 * t280 - t213 * t270 + t215 * t271) * t264 + (-t223 * t280 - t224 * t270 + t225 * t271) * t275) / 0.2e1 + t220 * ((t197 * t270 + t199 * t230 + t201 * t231) * t220 + (t196 * t270 + t198 * t230 + t200 * t231) * t219 + (t205 * t270 + t206 * t230 + t207 * t231) * t241) / 0.2e1 + m(3) * (t254 ^ 2 + t273 ^ 2 + t274 ^ 2) / 0.2e1 + t219 * ((t197 * t268 + t199 * t228 + t201 * t229) * t220 + (t196 * t268 + t198 * t228 + t200 * t229) * t219 + (t205 * t268 + t206 * t228 + t207 * t229) * t241) / 0.2e1 + m(2) * (t317 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + t275 * ((-t212 * t290 - t214 * t282 + t216 * t283) * t265 + (-t211 * t290 - t213 * t282 + t215 * t283) * t264 + (-t223 * t290 - t224 * t282 + t225 * t283) * t275) / 0.2e1 + t241 * ((t197 * t282 + t199 * t247 + t201 * t248) * t220 + (t196 * t282 + t198 * t247 + t200 * t248) * t219 + (t205 * t282 + t206 * t247 + t207 * t248) * t241) / 0.2e1 + t353 * ((-t301 * t412 + t303 * t340 + t305 * t341) * t354 + (-t300 * t412 + t302 * t340 + t304 * t341) * t353 + (-t327 * t412 + t328 * t340 + t329 * t341) * t363) / 0.2e1 + t354 * ((t301 * t414 + t303 * t342 + t305 * t343) * t354 + (t300 * t414 + t302 * t342 + t304 * t343) * t353 + (t327 * t414 + t328 * t342 + t329 * t343) * t363) / 0.2e1 + t363 * ((t300 * t353 + t301 * t354 + t327 * t363) * t375 + ((t303 * t382 + t305 * t379) * t354 + (t302 * t382 + t304 * t379) * t353 + (t328 * t382 + t329 * t379) * t363) * t372) / 0.2e1 + ((t250 * t278 + t251 * t279 + t285 * t294 + t286 * t295 - t399 * t427) * t322 + (t235 * t278 + t237 * t279 + t258 * t294 + t260 * t295 - t399 * t428) * t309 + (t234 * t278 + t236 * t279 + t257 * t294 + t259 * t295 - t399 * t429) * t308) * t308 / 0.2e1 + ((t250 * t280 + t251 * t281 + t285 * t296 + t286 * t297 + t321 * t427) * t322 + (t235 * t280 + t237 * t281 + t258 * t296 + t260 * t297 + t428 * t321) * t309 + (t234 * t280 + t236 * t281 + t257 * t296 + t259 * t297 + t321 * t429) * t308) * t309 / 0.2e1 + ((t250 * t290 + t251 * t291 + t285 * t318 + t286 * t319 + t427 * t339) * t322 + (t235 * t290 + t237 * t291 + t258 * t318 + t260 * t319 + t339 * t428) * t309 + (t234 * t290 + t236 * t291 + t257 * t318 + t259 * t319 + t339 * t429) * t308) * t322 / 0.2e1 + ((-t357 * t370 + t359 * t373 + Icges(1,4)) * V_base(5) + (-t358 * t370 + t360 * t373 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t357 * t373 + t359 * t370 + Icges(1,2)) * V_base(5) + (t358 * t373 + t360 * t370 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t373 - Icges(2,6) * t370 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t370 + Icges(2,6) * t373 + Icges(1,6)) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
