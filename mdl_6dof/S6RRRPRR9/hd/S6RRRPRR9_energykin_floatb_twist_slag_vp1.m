% Calculate kinetic energy for
% S6RRRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 19:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR9_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR9_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR9_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR9_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR9_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR9_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:56:27
% EndTime: 2019-03-09 18:56:32
% DurationCPUTime: 5.34s
% Computational Cost: add. (6069->454), mult. (15638->664), div. (0->0), fcn. (20097->16), ass. (0->202)
t432 = Icges(4,3) + Icges(5,3);
t376 = cos(pkin(6));
t381 = sin(qJ(1));
t384 = cos(qJ(2));
t408 = t381 * t384;
t380 = sin(qJ(2));
t385 = cos(qJ(1));
t410 = t380 * t385;
t345 = -t376 * t408 - t410;
t373 = sin(pkin(7));
t375 = cos(pkin(7));
t374 = sin(pkin(6));
t414 = t374 * t381;
t322 = -t345 * t373 + t375 * t414;
t431 = pkin(10) * t322;
t407 = t384 * t385;
t409 = t381 * t380;
t343 = t376 * t407 - t409;
t412 = t374 * t385;
t402 = t343 * t373 + t375 * t412;
t430 = pkin(10) * t402;
t379 = sin(qJ(3));
t383 = cos(qJ(3));
t418 = sin(pkin(13));
t419 = cos(pkin(13));
t349 = -t379 * t418 + t383 * t419;
t335 = t349 * t375;
t344 = t376 * t410 + t408;
t348 = -t379 * t419 - t383 * t418;
t396 = t373 * t349;
t393 = t374 * t396;
t282 = t335 * t343 + t344 * t348 - t385 * t393;
t334 = t348 * t373;
t336 = t348 * t375;
t283 = t334 * t412 - t343 * t336 + t344 * t349;
t401 = t343 * t375 - t373 * t412;
t296 = -t344 * t379 + t383 * t401;
t297 = t344 * t383 + t379 * t401;
t429 = Icges(4,5) * t297 + Icges(5,5) * t283 + Icges(4,6) * t296 + Icges(5,6) * t282 - t432 * t402;
t346 = -t376 * t409 + t407;
t284 = t345 * t335 + t346 * t348 + t381 * t393;
t285 = -t334 * t414 - t336 * t345 + t346 * t349;
t400 = t345 * t375 + t373 * t414;
t298 = -t346 * t379 + t383 * t400;
t299 = t346 * t383 + t379 * t400;
t428 = Icges(4,5) * t299 + Icges(5,5) * t285 + Icges(4,6) * t298 + Icges(5,6) * t284 + t432 * t322;
t413 = t374 * t384;
t415 = t374 * t380;
t292 = t335 * t413 + t348 * t415 + t376 * t396;
t293 = -t334 * t376 + (-t336 * t384 + t349 * t380) * t374;
t411 = t375 * t384;
t319 = t373 * t376 * t383 + (-t379 * t380 + t383 * t411) * t374;
t416 = t373 * t379;
t320 = t376 * t416 + (t379 * t411 + t380 * t383) * t374;
t342 = -t373 * t413 + t375 * t376;
t427 = Icges(4,5) * t320 + Icges(5,5) * t293 + Icges(4,6) * t319 + Icges(5,6) * t292 + t432 * t342;
t423 = cos(qJ(5));
t422 = pkin(9) * t376;
t421 = pkin(3) * t383;
t420 = pkin(10) + qJ(4);
t417 = Icges(2,4) * t381;
t406 = qJD(2) * t374;
t405 = V_base(5) * pkin(8) + V_base(1);
t356 = t381 * t406 + V_base(4);
t370 = V_base(6) + qJD(1);
t311 = qJD(3) * t322 + t356;
t357 = qJD(2) * t376 + t370;
t267 = -qJD(5) * t284 + t311;
t323 = qJD(3) * t342 + t357;
t355 = -t385 * t406 + V_base(5);
t350 = t381 * pkin(1) - pkin(9) * t412;
t399 = -t350 * t370 + V_base(5) * t422 + t405;
t277 = -qJD(5) * t292 + t323;
t351 = pkin(1) * t385 + pkin(9) * t414;
t398 = V_base(4) * t350 - t351 * V_base(5) + V_base(3);
t310 = -qJD(3) * t402 + t355;
t266 = -qJD(5) * t282 + t310;
t397 = t370 * t351 + V_base(2) + (-pkin(8) - t422) * V_base(4);
t300 = t344 * pkin(2) - t430;
t327 = pkin(2) * t415 + pkin(10) * t342;
t395 = -t300 * t357 + t355 * t327 + t399;
t301 = pkin(2) * t346 + t431;
t394 = t356 * t300 - t301 * t355 + t398;
t340 = pkin(3) * t416 + t375 * t420;
t341 = pkin(3) * t375 * t379 - t373 * t420;
t290 = (-pkin(10) * t375 + t340) * t376 + ((pkin(10) * t373 + t341) * t384 + t421 * t380) * t374;
t392 = qJD(4) * t322 + t310 * t290 + t395;
t272 = -t340 * t412 + t343 * t341 + t344 * t421 + t430;
t391 = qJD(4) * t342 + t311 * t272 + t394;
t390 = t357 * t301 - t327 * t356 + t397;
t273 = t340 * t414 + t341 * t345 + t346 * t421 - t431;
t389 = -qJD(4) * t402 + t323 * t273 + t390;
t245 = pkin(4) * t283 - pkin(11) * t282;
t265 = pkin(4) * t293 - pkin(11) * t292;
t388 = t310 * t265 + (-t245 - t272) * t323 + t392;
t246 = pkin(4) * t285 - pkin(11) * t284;
t387 = t311 * t245 + (-t246 - t273) * t310 + t391;
t386 = t323 * t246 + (-t265 - t290) * t311 + t389;
t382 = cos(qJ(6));
t378 = sin(qJ(5));
t377 = sin(qJ(6));
t371 = Icges(2,4) * t385;
t365 = rSges(2,1) * t385 - t381 * rSges(2,2);
t364 = t381 * rSges(2,1) + rSges(2,2) * t385;
t363 = Icges(2,1) * t385 - t417;
t362 = Icges(2,1) * t381 + t371;
t361 = -Icges(2,2) * t381 + t371;
t360 = Icges(2,2) * t385 + t417;
t354 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t353 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t352 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t332 = rSges(3,3) * t376 + (rSges(3,1) * t380 + rSges(3,2) * t384) * t374;
t331 = Icges(3,5) * t376 + (Icges(3,1) * t380 + Icges(3,4) * t384) * t374;
t330 = Icges(3,6) * t376 + (Icges(3,4) * t380 + Icges(3,2) * t384) * t374;
t329 = Icges(3,3) * t376 + (Icges(3,5) * t380 + Icges(3,6) * t384) * t374;
t326 = V_base(5) * rSges(2,3) - t364 * t370 + t405;
t325 = t365 * t370 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t324 = t364 * V_base(4) - t365 * V_base(5) + V_base(3);
t309 = rSges(3,1) * t346 + rSges(3,2) * t345 + rSges(3,3) * t414;
t308 = t344 * rSges(3,1) + t343 * rSges(3,2) - rSges(3,3) * t412;
t307 = Icges(3,1) * t346 + Icges(3,4) * t345 + Icges(3,5) * t414;
t306 = Icges(3,1) * t344 + Icges(3,4) * t343 - Icges(3,5) * t412;
t305 = Icges(3,4) * t346 + Icges(3,2) * t345 + Icges(3,6) * t414;
t304 = Icges(3,4) * t344 + Icges(3,2) * t343 - Icges(3,6) * t412;
t303 = Icges(3,5) * t346 + Icges(3,6) * t345 + Icges(3,3) * t414;
t302 = Icges(3,5) * t344 + Icges(3,6) * t343 - Icges(3,3) * t412;
t289 = rSges(4,1) * t320 + rSges(4,2) * t319 + rSges(4,3) * t342;
t288 = Icges(4,1) * t320 + Icges(4,4) * t319 + Icges(4,5) * t342;
t287 = Icges(4,4) * t320 + Icges(4,2) * t319 + Icges(4,6) * t342;
t281 = t293 * t423 + t342 * t378;
t280 = t293 * t378 - t342 * t423;
t276 = -t308 * t357 + t332 * t355 + t399;
t275 = t309 * t357 - t332 * t356 + t397;
t271 = t285 * t423 + t322 * t378;
t270 = t285 * t378 - t322 * t423;
t269 = t283 * t423 - t378 * t402;
t268 = t283 * t378 + t402 * t423;
t264 = rSges(4,1) * t299 + rSges(4,2) * t298 + rSges(4,3) * t322;
t263 = rSges(4,1) * t297 + rSges(4,2) * t296 - rSges(4,3) * t402;
t262 = t308 * t356 - t309 * t355 + t398;
t261 = Icges(4,1) * t299 + Icges(4,4) * t298 + Icges(4,5) * t322;
t260 = Icges(4,1) * t297 + Icges(4,4) * t296 - Icges(4,5) * t402;
t259 = Icges(4,4) * t299 + Icges(4,2) * t298 + Icges(4,6) * t322;
t258 = Icges(4,4) * t297 + Icges(4,2) * t296 - Icges(4,6) * t402;
t254 = rSges(5,1) * t293 + rSges(5,2) * t292 + rSges(5,3) * t342;
t253 = Icges(5,1) * t293 + Icges(5,4) * t292 + Icges(5,5) * t342;
t252 = Icges(5,4) * t293 + Icges(5,2) * t292 + Icges(5,6) * t342;
t250 = t281 * t382 - t292 * t377;
t249 = -t281 * t377 - t292 * t382;
t244 = pkin(5) * t281 + pkin(12) * t280;
t243 = qJD(6) * t280 + t277;
t241 = rSges(5,1) * t285 + rSges(5,2) * t284 + rSges(5,3) * t322;
t240 = rSges(5,1) * t283 + rSges(5,2) * t282 - rSges(5,3) * t402;
t239 = Icges(5,1) * t285 + Icges(5,4) * t284 + Icges(5,5) * t322;
t238 = Icges(5,1) * t283 + Icges(5,4) * t282 - Icges(5,5) * t402;
t237 = Icges(5,4) * t285 + Icges(5,2) * t284 + Icges(5,6) * t322;
t236 = Icges(5,4) * t283 + Icges(5,2) * t282 - Icges(5,6) * t402;
t233 = t271 * t382 - t284 * t377;
t232 = -t271 * t377 - t284 * t382;
t231 = t269 * t382 - t282 * t377;
t230 = -t269 * t377 - t282 * t382;
t228 = pkin(5) * t271 + pkin(12) * t270;
t227 = pkin(5) * t269 + pkin(12) * t268;
t226 = rSges(6,1) * t281 - rSges(6,2) * t280 - rSges(6,3) * t292;
t225 = Icges(6,1) * t281 - Icges(6,4) * t280 - Icges(6,5) * t292;
t224 = Icges(6,4) * t281 - Icges(6,2) * t280 - Icges(6,6) * t292;
t223 = Icges(6,5) * t281 - Icges(6,6) * t280 - Icges(6,3) * t292;
t222 = qJD(6) * t270 + t267;
t221 = qJD(6) * t268 + t266;
t220 = rSges(6,1) * t271 - rSges(6,2) * t270 - rSges(6,3) * t284;
t219 = rSges(6,1) * t269 - rSges(6,2) * t268 - rSges(6,3) * t282;
t218 = Icges(6,1) * t271 - Icges(6,4) * t270 - Icges(6,5) * t284;
t217 = Icges(6,1) * t269 - Icges(6,4) * t268 - Icges(6,5) * t282;
t216 = Icges(6,4) * t271 - Icges(6,2) * t270 - Icges(6,6) * t284;
t215 = Icges(6,4) * t269 - Icges(6,2) * t268 - Icges(6,6) * t282;
t214 = Icges(6,5) * t271 - Icges(6,6) * t270 - Icges(6,3) * t284;
t213 = Icges(6,5) * t269 - Icges(6,6) * t268 - Icges(6,3) * t282;
t212 = -t263 * t323 + t289 * t310 + t395;
t211 = t264 * t323 - t289 * t311 + t390;
t210 = rSges(7,1) * t250 + rSges(7,2) * t249 + rSges(7,3) * t280;
t209 = Icges(7,1) * t250 + Icges(7,4) * t249 + Icges(7,5) * t280;
t208 = Icges(7,4) * t250 + Icges(7,2) * t249 + Icges(7,6) * t280;
t207 = Icges(7,5) * t250 + Icges(7,6) * t249 + Icges(7,3) * t280;
t206 = t263 * t311 - t264 * t310 + t394;
t205 = rSges(7,1) * t233 + rSges(7,2) * t232 + rSges(7,3) * t270;
t204 = rSges(7,1) * t231 + rSges(7,2) * t230 + rSges(7,3) * t268;
t203 = Icges(7,1) * t233 + Icges(7,4) * t232 + Icges(7,5) * t270;
t202 = Icges(7,1) * t231 + Icges(7,4) * t230 + Icges(7,5) * t268;
t201 = Icges(7,4) * t233 + Icges(7,2) * t232 + Icges(7,6) * t270;
t200 = Icges(7,4) * t231 + Icges(7,2) * t230 + Icges(7,6) * t268;
t199 = Icges(7,5) * t233 + Icges(7,6) * t232 + Icges(7,3) * t270;
t198 = Icges(7,5) * t231 + Icges(7,6) * t230 + Icges(7,3) * t268;
t197 = t254 * t310 + (-t240 - t272) * t323 + t392;
t196 = t241 * t323 + (-t254 - t290) * t311 + t389;
t195 = t240 * t311 + (-t241 - t273) * t310 + t391;
t194 = -t219 * t277 + t226 * t266 + t388;
t193 = t220 * t277 - t226 * t267 + t386;
t192 = t219 * t267 - t220 * t266 + t387;
t191 = -t204 * t243 + t210 * t221 - t227 * t277 + t244 * t266 + t388;
t190 = t205 * t243 - t210 * t222 + t228 * t277 - t244 * t267 + t386;
t189 = t204 * t222 - t205 * t221 + t227 * t267 - t228 * t266 + t387;
t1 = m(4) * (t206 ^ 2 + t211 ^ 2 + t212 ^ 2) / 0.2e1 + m(5) * (t195 ^ 2 + t196 ^ 2 + t197 ^ 2) / 0.2e1 + m(6) * (t192 ^ 2 + t193 ^ 2 + t194 ^ 2) / 0.2e1 + m(7) * (t189 ^ 2 + t190 ^ 2 + t191 ^ 2) / 0.2e1 + m(1) * (t352 ^ 2 + t353 ^ 2 + t354 ^ 2) / 0.2e1 + m(2) * (t324 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + t277 * ((-t214 * t292 - t216 * t280 + t218 * t281) * t267 + (-t213 * t292 - t215 * t280 + t217 * t281) * t266 + (-t223 * t292 - t280 * t224 + t281 * t225) * t277) / 0.2e1 + t266 * ((-t214 * t282 - t216 * t268 + t218 * t269) * t267 + (-t282 * t213 - t268 * t215 + t269 * t217) * t266 + (-t223 * t282 - t224 * t268 + t225 * t269) * t277) / 0.2e1 + t267 * ((-t284 * t214 - t216 * t270 + t271 * t218) * t267 + (-t213 * t284 - t215 * t270 + t217 * t271) * t266 + (-t223 * t284 - t224 * t270 + t225 * t271) * t277) / 0.2e1 + t243 * ((t199 * t280 + t201 * t249 + t203 * t250) * t222 + (t198 * t280 + t200 * t249 + t202 * t250) * t221 + (t280 * t207 + t249 * t208 + t250 * t209) * t243) / 0.2e1 + m(3) * (t262 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + t222 * ((t270 * t199 + t232 * t201 + t233 * t203) * t222 + (t198 * t270 + t200 * t232 + t202 * t233) * t221 + (t207 * t270 + t208 * t232 + t209 * t233) * t243) / 0.2e1 + t221 * ((t199 * t268 + t201 * t230 + t203 * t231) * t222 + (t198 * t268 + t230 * t200 + t231 * t202) * t221 + (t207 * t268 + t208 * t230 + t209 * t231) * t243) / 0.2e1 + t355 * ((-t303 * t412 + t343 * t305 + t344 * t307) * t356 + (-t302 * t412 + t343 * t304 + t344 * t306) * t355 + (-t329 * t412 + t343 * t330 + t344 * t331) * t357) / 0.2e1 + t356 * ((t303 * t414 + t305 * t345 + t307 * t346) * t356 + (t302 * t414 + t304 * t345 + t306 * t346) * t355 + (t329 * t414 + t330 * t345 + t331 * t346) * t357) / 0.2e1 + t357 * ((t302 * t355 + t303 * t356 + t329 * t357) * t376 + ((t305 * t384 + t307 * t380) * t356 + (t304 * t384 + t306 * t380) * t355 + (t330 * t384 + t331 * t380) * t357) * t374) / 0.2e1 + ((t252 * t282 + t253 * t283 + t287 * t296 + t288 * t297 - t402 * t427) * t323 + (t237 * t282 + t239 * t283 + t259 * t296 + t261 * t297 - t402 * t428) * t311 + (t282 * t236 + t283 * t238 + t296 * t258 + t297 * t260 - t402 * t429) * t310) * t310 / 0.2e1 + ((t252 * t284 + t253 * t285 + t287 * t298 + t288 * t299 + t322 * t427) * t323 + (t237 * t284 + t285 * t239 + t298 * t259 + t299 * t261 + t428 * t322) * t311 + (t236 * t284 + t238 * t285 + t258 * t298 + t260 * t299 + t322 * t429) * t310) * t311 / 0.2e1 + ((t252 * t292 + t253 * t293 + t287 * t319 + t288 * t320 + t427 * t342) * t323 + (t237 * t292 + t239 * t293 + t259 * t319 + t261 * t320 + t342 * t428) * t311 + (t236 * t292 + t238 * t293 + t258 * t319 + t260 * t320 + t342 * t429) * t310) * t323 / 0.2e1 + ((-t381 * t360 + t362 * t385 + Icges(1,4)) * V_base(5) + (-t381 * t361 + t363 * t385 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t360 * t385 + t381 * t362 + Icges(1,2)) * V_base(5) + (t361 * t385 + t381 * t363 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t381 + Icges(2,6) * t385) * V_base(5) + (Icges(2,5) * t385 - Icges(2,6) * t381) * V_base(4) + Icges(2,3) * t370 / 0.2e1) * t370;
T  = t1;
