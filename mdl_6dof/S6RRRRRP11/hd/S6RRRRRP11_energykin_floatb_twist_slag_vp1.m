% Calculate kinetic energy for
% S6RRRRRP11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 02:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP11_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP11_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRP11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP11_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP11_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP11_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP11_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:32:29
% EndTime: 2019-03-10 02:32:34
% DurationCPUTime: 5.14s
% Computational Cost: add. (5682->396), mult. (14811->585), div. (0->0), fcn. (18967->14), ass. (0->182)
t422 = Icges(6,1) + Icges(7,1);
t421 = Icges(6,4) + Icges(7,4);
t420 = Icges(6,5) + Icges(7,5);
t419 = Icges(6,2) + Icges(7,2);
t418 = Icges(6,6) + Icges(7,6);
t417 = Icges(6,3) + Icges(7,3);
t416 = rSges(7,3) + qJ(6);
t349 = cos(pkin(6));
t357 = cos(qJ(2));
t358 = cos(qJ(1));
t386 = t357 * t358;
t354 = sin(qJ(2));
t355 = sin(qJ(1));
t388 = t355 * t354;
t317 = t349 * t386 - t388;
t387 = t355 * t357;
t389 = t354 * t358;
t318 = t349 * t389 + t387;
t353 = sin(qJ(3));
t348 = sin(pkin(6));
t397 = sin(pkin(7));
t377 = t348 * t397;
t398 = cos(pkin(7));
t403 = cos(qJ(3));
t280 = t318 * t403 + (t317 * t398 - t358 * t377) * t353;
t352 = sin(qJ(4));
t378 = t348 * t398;
t370 = -t317 * t397 - t358 * t378;
t402 = cos(qJ(4));
t264 = t280 * t402 + t352 * t370;
t375 = t403 * t397;
t374 = t348 * t375;
t376 = t398 * t403;
t279 = -t317 * t376 + t318 * t353 + t358 * t374;
t351 = sin(qJ(5));
t356 = cos(qJ(5));
t234 = -t264 * t351 + t279 * t356;
t395 = t279 * t351;
t235 = t264 * t356 + t395;
t263 = t280 * t352 - t370 * t402;
t415 = t234 * t418 + t235 * t420 + t263 * t417;
t319 = -t349 * t387 - t389;
t320 = -t349 * t388 + t386;
t282 = t320 * t403 + (t319 * t398 + t355 * t377) * t353;
t369 = -t319 * t397 + t355 * t378;
t266 = t282 * t402 + t352 * t369;
t281 = -t319 * t376 + t320 * t353 - t355 * t374;
t236 = -t266 * t351 + t281 * t356;
t394 = t281 * t351;
t237 = t266 * t356 + t394;
t265 = t282 * t352 - t369 * t402;
t414 = t236 * t418 + t237 * t420 + t265 * t417;
t413 = t234 * t419 + t235 * t421 + t263 * t418;
t412 = t236 * t419 + t237 * t421 + t265 * t418;
t411 = t421 * t234 + t235 * t422 + t420 * t263;
t410 = t421 * t236 + t237 * t422 + t420 * t265;
t302 = t349 * t397 * t353 + (t353 * t357 * t398 + t354 * t403) * t348;
t368 = t349 * t398 - t357 * t377;
t278 = t302 * t402 + t352 * t368;
t392 = t348 * t354;
t301 = -t348 * t357 * t376 - t349 * t375 + t353 * t392;
t261 = -t278 * t351 + t301 * t356;
t393 = t301 * t351;
t262 = t278 * t356 + t393;
t277 = t302 * t352 - t368 * t402;
t409 = t261 * t418 + t262 * t420 + t277 * t417;
t408 = t261 * t419 + t262 * t421 + t277 * t418;
t407 = t421 * t261 + t262 * t422 + t420 * t277;
t401 = pkin(9) * t349;
t400 = pkin(5) * t356;
t396 = Icges(2,4) * t355;
t391 = t348 * t355;
t390 = t348 * t358;
t385 = rSges(7,1) * t235 + rSges(7,2) * t234 + pkin(5) * t395 + t263 * t416 + t264 * t400;
t384 = rSges(7,1) * t237 + rSges(7,2) * t236 + pkin(5) * t394 + t265 * t416 + t266 * t400;
t383 = rSges(7,1) * t262 + rSges(7,2) * t261 + pkin(5) * t393 + t277 * t416 + t278 * t400;
t382 = qJD(2) * t348;
t381 = V_base(5) * pkin(8) + V_base(1);
t330 = t355 * t382 + V_base(4);
t345 = V_base(6) + qJD(1);
t294 = qJD(3) * t369 + t330;
t331 = qJD(2) * t349 + t345;
t260 = qJD(4) * t281 + t294;
t303 = qJD(3) * t368 + t331;
t329 = -t358 * t382 + V_base(5);
t322 = t355 * pkin(1) - pkin(9) * t390;
t373 = -t322 * t345 + V_base(5) * t401 + t381;
t272 = qJD(4) * t301 + t303;
t323 = pkin(1) * t358 + pkin(9) * t391;
t372 = V_base(4) * t322 - t323 * V_base(5) + V_base(3);
t293 = qJD(3) * t370 + t329;
t259 = qJD(4) * t279 + t293;
t371 = t345 * t323 + V_base(2) + (-pkin(8) - t401) * V_base(4);
t283 = t318 * pkin(2) + pkin(10) * t370;
t307 = pkin(2) * t392 + pkin(10) * t368;
t367 = -t283 * t331 + t329 * t307 + t373;
t284 = t320 * pkin(2) + pkin(10) * t369;
t366 = t330 * t283 - t284 * t329 + t372;
t365 = t331 * t284 - t307 * t330 + t371;
t256 = pkin(3) * t280 + pkin(11) * t279;
t271 = pkin(3) * t302 + pkin(11) * t301;
t364 = -t256 * t303 + t293 * t271 + t367;
t257 = pkin(3) * t282 + pkin(11) * t281;
t363 = t294 * t256 - t257 * t293 + t366;
t362 = t303 * t257 - t271 * t294 + t365;
t231 = pkin(4) * t264 + pkin(12) * t263;
t255 = pkin(4) * t278 + pkin(12) * t277;
t361 = -t231 * t272 + t259 * t255 + t364;
t232 = pkin(4) * t266 + pkin(12) * t265;
t360 = t260 * t231 - t232 * t259 + t363;
t359 = t272 * t232 - t255 * t260 + t362;
t346 = Icges(2,4) * t358;
t339 = rSges(2,1) * t358 - t355 * rSges(2,2);
t338 = t355 * rSges(2,1) + rSges(2,2) * t358;
t337 = Icges(2,1) * t358 - t396;
t336 = Icges(2,1) * t355 + t346;
t335 = -Icges(2,2) * t355 + t346;
t334 = Icges(2,2) * t358 + t396;
t328 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t327 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t326 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t313 = rSges(3,3) * t349 + (rSges(3,1) * t354 + rSges(3,2) * t357) * t348;
t312 = Icges(3,5) * t349 + (Icges(3,1) * t354 + Icges(3,4) * t357) * t348;
t311 = Icges(3,6) * t349 + (Icges(3,4) * t354 + Icges(3,2) * t357) * t348;
t310 = Icges(3,3) * t349 + (Icges(3,5) * t354 + Icges(3,6) * t357) * t348;
t306 = V_base(5) * rSges(2,3) - t338 * t345 + t381;
t305 = t339 * t345 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t304 = t338 * V_base(4) - t339 * V_base(5) + V_base(3);
t292 = rSges(3,1) * t320 + rSges(3,2) * t319 + rSges(3,3) * t391;
t291 = t318 * rSges(3,1) + t317 * rSges(3,2) - rSges(3,3) * t390;
t290 = Icges(3,1) * t320 + Icges(3,4) * t319 + Icges(3,5) * t391;
t289 = Icges(3,1) * t318 + Icges(3,4) * t317 - Icges(3,5) * t390;
t288 = Icges(3,4) * t320 + Icges(3,2) * t319 + Icges(3,6) * t391;
t287 = Icges(3,4) * t318 + Icges(3,2) * t317 - Icges(3,6) * t390;
t286 = Icges(3,5) * t320 + Icges(3,6) * t319 + Icges(3,3) * t391;
t285 = Icges(3,5) * t318 + Icges(3,6) * t317 - Icges(3,3) * t390;
t270 = t302 * rSges(4,1) - t301 * rSges(4,2) + rSges(4,3) * t368;
t269 = Icges(4,1) * t302 - Icges(4,4) * t301 + Icges(4,5) * t368;
t268 = Icges(4,4) * t302 - Icges(4,2) * t301 + Icges(4,6) * t368;
t267 = Icges(4,5) * t302 - Icges(4,6) * t301 + Icges(4,3) * t368;
t254 = -t291 * t331 + t313 * t329 + t373;
t253 = t292 * t331 - t313 * t330 + t371;
t252 = qJD(5) * t277 + t272;
t250 = t282 * rSges(4,1) - t281 * rSges(4,2) + rSges(4,3) * t369;
t249 = t280 * rSges(4,1) - t279 * rSges(4,2) + rSges(4,3) * t370;
t248 = t291 * t330 - t292 * t329 + t372;
t247 = Icges(4,1) * t282 - Icges(4,4) * t281 + Icges(4,5) * t369;
t246 = Icges(4,1) * t280 - Icges(4,4) * t279 + Icges(4,5) * t370;
t245 = Icges(4,4) * t282 - Icges(4,2) * t281 + Icges(4,6) * t369;
t244 = Icges(4,4) * t280 - Icges(4,2) * t279 + Icges(4,6) * t370;
t243 = Icges(4,5) * t282 - Icges(4,6) * t281 + Icges(4,3) * t369;
t242 = Icges(4,5) * t280 - Icges(4,6) * t279 + Icges(4,3) * t370;
t241 = rSges(5,1) * t278 - rSges(5,2) * t277 + rSges(5,3) * t301;
t240 = Icges(5,1) * t278 - Icges(5,4) * t277 + Icges(5,5) * t301;
t239 = Icges(5,4) * t278 - Icges(5,2) * t277 + Icges(5,6) * t301;
t238 = Icges(5,5) * t278 - Icges(5,6) * t277 + Icges(5,3) * t301;
t230 = qJD(5) * t265 + t260;
t229 = qJD(5) * t263 + t259;
t227 = rSges(5,1) * t266 - rSges(5,2) * t265 + rSges(5,3) * t281;
t226 = rSges(5,1) * t264 - rSges(5,2) * t263 + rSges(5,3) * t279;
t225 = Icges(5,1) * t266 - Icges(5,4) * t265 + Icges(5,5) * t281;
t224 = Icges(5,1) * t264 - Icges(5,4) * t263 + Icges(5,5) * t279;
t223 = Icges(5,4) * t266 - Icges(5,2) * t265 + Icges(5,6) * t281;
t222 = Icges(5,4) * t264 - Icges(5,2) * t263 + Icges(5,6) * t279;
t221 = Icges(5,5) * t266 - Icges(5,6) * t265 + Icges(5,3) * t281;
t220 = Icges(5,5) * t264 - Icges(5,6) * t263 + Icges(5,3) * t279;
t218 = rSges(6,1) * t262 + rSges(6,2) * t261 + rSges(6,3) * t277;
t208 = rSges(6,1) * t237 + rSges(6,2) * t236 + rSges(6,3) * t265;
t206 = rSges(6,1) * t235 + rSges(6,2) * t234 + rSges(6,3) * t263;
t192 = -t249 * t303 + t270 * t293 + t367;
t191 = t250 * t303 - t270 * t294 + t365;
t188 = t249 * t294 - t250 * t293 + t366;
t187 = -t226 * t272 + t241 * t259 + t364;
t186 = t227 * t272 - t241 * t260 + t362;
t185 = t226 * t260 - t227 * t259 + t363;
t184 = -t206 * t252 + t218 * t229 + t361;
t183 = t208 * t252 - t218 * t230 + t359;
t182 = t206 * t230 - t208 * t229 + t360;
t181 = qJD(6) * t265 + t229 * t383 - t252 * t385 + t361;
t180 = qJD(6) * t263 - t230 * t383 + t252 * t384 + t359;
t179 = qJD(6) * t277 - t229 * t384 + t230 * t385 + t360;
t1 = t331 * ((t285 * t329 + t286 * t330 + t310 * t331) * t349 + ((t288 * t357 + t290 * t354) * t330 + (t287 * t357 + t289 * t354) * t329 + (t311 * t357 + t312 * t354) * t331) * t348) / 0.2e1 + t303 * ((t243 * t368 - t301 * t245 + t302 * t247) * t294 + (t242 * t368 - t301 * t244 + t302 * t246) * t293 + (t267 * t368 - t301 * t268 + t302 * t269) * t303) / 0.2e1 + t294 * ((t243 * t369 - t281 * t245 + t282 * t247) * t294 + (t242 * t369 - t281 * t244 + t282 * t246) * t293 + (t267 * t369 - t281 * t268 + t282 * t269) * t303) / 0.2e1 + t293 * ((t243 * t370 - t279 * t245 + t280 * t247) * t294 + (t242 * t370 - t279 * t244 + t280 * t246) * t293 + (t267 * t370 - t279 * t268 + t280 * t269) * t303) / 0.2e1 + m(1) * (t326 ^ 2 + t327 ^ 2 + t328 ^ 2) / 0.2e1 + t272 * ((t221 * t301 - t223 * t277 + t225 * t278) * t260 + (t220 * t301 - t222 * t277 + t224 * t278) * t259 + (t238 * t301 - t239 * t277 + t240 * t278) * t272) / 0.2e1 + m(2) * (t304 ^ 2 + t305 ^ 2 + t306 ^ 2) / 0.2e1 + t259 * ((t221 * t279 - t223 * t263 + t225 * t264) * t260 + (t220 * t279 - t222 * t263 + t224 * t264) * t259 + (t238 * t279 - t239 * t263 + t240 * t264) * t272) / 0.2e1 + t260 * ((t221 * t281 - t223 * t265 + t225 * t266) * t260 + (t220 * t281 - t222 * t265 + t224 * t266) * t259 + (t238 * t281 - t239 * t265 + t240 * t266) * t272) / 0.2e1 + m(3) * (t248 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + m(4) * (t188 ^ 2 + t191 ^ 2 + t192 ^ 2) / 0.2e1 + m(5) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(6) * (t182 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + m(7) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + t329 * ((-t286 * t390 + t317 * t288 + t318 * t290) * t330 + (-t285 * t390 + t317 * t287 + t318 * t289) * t329 + (-t310 * t390 + t317 * t311 + t318 * t312) * t331) / 0.2e1 + t330 * ((t286 * t391 + t288 * t319 + t290 * t320) * t330 + (t285 * t391 + t287 * t319 + t289 * t320) * t329 + (t310 * t391 + t311 * t319 + t312 * t320) * t331) / 0.2e1 + ((t234 * t408 + t235 * t407 + t263 * t409) * t252 + (t234 * t412 + t235 * t410 + t263 * t414) * t230 + (t413 * t234 + t411 * t235 + t415 * t263) * t229) * t229 / 0.2e1 + ((t236 * t408 + t237 * t407 + t265 * t409) * t252 + (t412 * t236 + t410 * t237 + t414 * t265) * t230 + (t413 * t236 + t411 * t237 + t265 * t415) * t229) * t230 / 0.2e1 + ((t408 * t261 + t407 * t262 + t409 * t277) * t252 + (t261 * t412 + t262 * t410 + t277 * t414) * t230 + (t413 * t261 + t411 * t262 + t277 * t415) * t229) * t252 / 0.2e1 + ((-t355 * t334 + t336 * t358 + Icges(1,4)) * V_base(5) + (-t355 * t335 + t337 * t358 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t334 * t358 + t355 * t336 + Icges(1,2)) * V_base(5) + (t335 * t358 + t355 * t337 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t355 + Icges(2,6) * t358) * V_base(5) + (Icges(2,5) * t358 - Icges(2,6) * t355) * V_base(4) + Icges(2,3) * t345 / 0.2e1) * t345;
T  = t1;
