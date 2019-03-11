% Calculate kinetic energy for
% S6RPRRRP12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP12_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRP12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP12_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP12_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:43:53
% EndTime: 2019-03-09 06:43:58
% DurationCPUTime: 4.87s
% Computational Cost: add. (5490->396), mult. (14471->566), div. (0->0), fcn. (18600->14), ass. (0->178)
t433 = Icges(6,1) + Icges(7,1);
t432 = -Icges(6,4) + Icges(7,5);
t431 = Icges(7,4) + Icges(6,5);
t430 = Icges(6,2) + Icges(7,3);
t429 = Icges(7,2) + Icges(6,3);
t428 = -Icges(6,6) + Icges(7,6);
t427 = rSges(7,1) + pkin(5);
t426 = rSges(7,3) + qJ(6);
t368 = cos(pkin(12));
t366 = sin(pkin(12));
t373 = sin(qJ(1));
t401 = t373 * t366;
t369 = cos(pkin(6));
t374 = cos(qJ(1));
t402 = t369 * t374;
t339 = t368 * t402 - t401;
t400 = t373 * t368;
t340 = t366 * t402 + t400;
t372 = sin(qJ(3));
t367 = sin(pkin(6));
t408 = sin(pkin(7));
t390 = t367 * t408;
t409 = cos(pkin(7));
t412 = cos(qJ(3));
t297 = t340 * t412 + (t339 * t409 - t374 * t390) * t372;
t391 = t367 * t409;
t321 = -t339 * t408 - t374 * t391;
t371 = sin(qJ(4));
t411 = cos(qJ(4));
t280 = t297 * t411 + t321 * t371;
t387 = t412 * t408;
t385 = t367 * t387;
t388 = t409 * t412;
t296 = -t339 * t388 + t340 * t372 + t374 * t385;
t370 = sin(qJ(5));
t410 = cos(qJ(5));
t249 = t280 * t370 - t296 * t410;
t250 = t280 * t410 + t296 * t370;
t279 = t297 * t371 - t321 * t411;
t425 = t430 * t249 + t432 * t250 + t428 * t279;
t341 = -t366 * t374 - t369 * t400;
t342 = t368 * t374 - t369 * t401;
t299 = t342 * t412 + (t341 * t409 + t373 * t390) * t372;
t322 = -t341 * t408 + t373 * t391;
t282 = t299 * t411 + t322 * t371;
t298 = -t341 * t388 + t342 * t372 - t373 * t385;
t251 = t282 * t370 - t298 * t410;
t252 = t282 * t410 + t298 * t370;
t281 = t299 * t371 - t322 * t411;
t424 = t430 * t251 + t432 * t252 + t428 * t281;
t423 = t428 * t249 + t431 * t250 + t429 * t279;
t422 = t428 * t251 + t431 * t252 + t429 * t281;
t421 = t432 * t249 + t433 * t250 + t431 * t279;
t420 = t432 * t251 + t433 * t252 + t431 * t281;
t320 = t369 * t408 * t372 + (t368 * t372 * t409 + t366 * t412) * t367;
t338 = -t368 * t390 + t369 * t409;
t294 = t320 * t411 + t338 * t371;
t405 = t366 * t367;
t319 = -t367 * t368 * t388 - t369 * t387 + t372 * t405;
t275 = t294 * t370 - t319 * t410;
t276 = t294 * t410 + t319 * t370;
t293 = t320 * t371 - t338 * t411;
t419 = t430 * t275 + t432 * t276 + t428 * t293;
t418 = t428 * t275 + t431 * t276 + t429 * t293;
t417 = t432 * t275 + t433 * t276 + t431 * t293;
t407 = Icges(2,4) * t373;
t406 = qJ(2) * t369;
t404 = t367 * t373;
t403 = t367 * t374;
t399 = rSges(7,2) * t279 + t426 * t249 + t250 * t427;
t398 = rSges(7,2) * t281 + t426 * t251 + t252 * t427;
t397 = rSges(7,2) * t293 + t426 * t275 + t276 * t427;
t396 = qJD(2) * t367;
t395 = V_base(5) * pkin(8) + V_base(1);
t312 = qJD(3) * t322 + V_base(4);
t311 = qJD(3) * t321 + V_base(5);
t363 = V_base(6) + qJD(1);
t392 = -pkin(8) - t406;
t344 = t373 * pkin(1) - qJ(2) * t403;
t389 = qJD(2) * t369 + V_base(4) * t344 + V_base(3);
t278 = qJD(4) * t298 + t312;
t277 = qJD(4) * t296 + t311;
t328 = qJD(3) * t338 + t363;
t386 = t373 * t396 + V_base(5) * t406 + t395;
t290 = qJD(4) * t319 + t328;
t345 = pkin(1) * t374 + qJ(2) * t404;
t384 = t363 * t345 - t374 * t396 + V_base(2);
t301 = t340 * pkin(2) + pkin(9) * t321;
t325 = pkin(2) * t405 + pkin(9) * t338;
t383 = V_base(5) * t325 + (-t301 - t344) * t363 + t386;
t302 = t342 * pkin(2) + pkin(9) * t322;
t382 = V_base(4) * t301 + (-t302 - t345) * V_base(5) + t389;
t272 = pkin(3) * t297 + pkin(10) * t296;
t287 = pkin(3) * t320 + pkin(10) * t319;
t381 = -t272 * t328 + t311 * t287 + t383;
t273 = pkin(3) * t299 + pkin(10) * t298;
t380 = t312 * t272 - t273 * t311 + t382;
t379 = t363 * t302 + (-t325 + t392) * V_base(4) + t384;
t245 = pkin(4) * t280 + pkin(11) * t279;
t268 = pkin(4) * t294 + pkin(11) * t293;
t378 = -t245 * t290 + t277 * t268 + t381;
t246 = pkin(4) * t282 + pkin(11) * t281;
t377 = t278 * t245 - t246 * t277 + t380;
t376 = t328 * t273 - t312 * t287 + t379;
t375 = t290 * t246 - t278 * t268 + t376;
t364 = Icges(2,4) * t374;
t358 = rSges(2,1) * t374 - t373 * rSges(2,2);
t357 = t373 * rSges(2,1) + rSges(2,2) * t374;
t356 = Icges(2,1) * t374 - t407;
t355 = Icges(2,1) * t373 + t364;
t354 = -Icges(2,2) * t373 + t364;
t353 = Icges(2,2) * t374 + t407;
t352 = Icges(2,5) * t374 - Icges(2,6) * t373;
t351 = Icges(2,5) * t373 + Icges(2,6) * t374;
t350 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t349 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t348 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t334 = rSges(3,3) * t369 + (rSges(3,1) * t366 + rSges(3,2) * t368) * t367;
t333 = Icges(3,5) * t369 + (Icges(3,1) * t366 + Icges(3,4) * t368) * t367;
t332 = Icges(3,6) * t369 + (Icges(3,4) * t366 + Icges(3,2) * t368) * t367;
t331 = Icges(3,3) * t369 + (Icges(3,5) * t366 + Icges(3,6) * t368) * t367;
t327 = V_base(5) * rSges(2,3) - t357 * t363 + t395;
t326 = t358 * t363 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t324 = t357 * V_base(4) - t358 * V_base(5) + V_base(3);
t310 = rSges(3,1) * t342 + rSges(3,2) * t341 + rSges(3,3) * t404;
t309 = t340 * rSges(3,1) + t339 * rSges(3,2) - rSges(3,3) * t403;
t308 = Icges(3,1) * t342 + Icges(3,4) * t341 + Icges(3,5) * t404;
t307 = Icges(3,1) * t340 + Icges(3,4) * t339 - Icges(3,5) * t403;
t306 = Icges(3,4) * t342 + Icges(3,2) * t341 + Icges(3,6) * t404;
t305 = Icges(3,4) * t340 + Icges(3,2) * t339 - Icges(3,6) * t403;
t304 = Icges(3,5) * t342 + Icges(3,6) * t341 + Icges(3,3) * t404;
t303 = Icges(3,5) * t340 + Icges(3,6) * t339 - Icges(3,3) * t403;
t286 = rSges(4,1) * t320 - rSges(4,2) * t319 + rSges(4,3) * t338;
t285 = Icges(4,1) * t320 - Icges(4,4) * t319 + Icges(4,5) * t338;
t284 = Icges(4,4) * t320 - Icges(4,2) * t319 + Icges(4,6) * t338;
t283 = Icges(4,5) * t320 - Icges(4,6) * t319 + Icges(4,3) * t338;
t271 = t334 * V_base(5) + (-t309 - t344) * t363 + t386;
t270 = t363 * t310 + (-t334 + t392) * V_base(4) + t384;
t269 = qJD(5) * t293 + t290;
t267 = t309 * V_base(4) + (-t310 - t345) * V_base(5) + t389;
t265 = rSges(4,1) * t299 - rSges(4,2) * t298 + rSges(4,3) * t322;
t264 = rSges(4,1) * t297 - rSges(4,2) * t296 + rSges(4,3) * t321;
t263 = Icges(4,1) * t299 - Icges(4,4) * t298 + Icges(4,5) * t322;
t262 = Icges(4,1) * t297 - Icges(4,4) * t296 + Icges(4,5) * t321;
t261 = Icges(4,4) * t299 - Icges(4,2) * t298 + Icges(4,6) * t322;
t260 = Icges(4,4) * t297 - Icges(4,2) * t296 + Icges(4,6) * t321;
t259 = Icges(4,5) * t299 - Icges(4,6) * t298 + Icges(4,3) * t322;
t258 = Icges(4,5) * t297 - Icges(4,6) * t296 + Icges(4,3) * t321;
t257 = rSges(5,1) * t294 - rSges(5,2) * t293 + rSges(5,3) * t319;
t255 = Icges(5,1) * t294 - Icges(5,4) * t293 + Icges(5,5) * t319;
t254 = Icges(5,4) * t294 - Icges(5,2) * t293 + Icges(5,6) * t319;
t253 = Icges(5,5) * t294 - Icges(5,6) * t293 + Icges(5,3) * t319;
t248 = qJD(5) * t281 + t278;
t247 = qJD(5) * t279 + t277;
t242 = rSges(5,1) * t282 - rSges(5,2) * t281 + rSges(5,3) * t298;
t241 = rSges(5,1) * t280 - rSges(5,2) * t279 + rSges(5,3) * t296;
t240 = Icges(5,1) * t282 - Icges(5,4) * t281 + Icges(5,5) * t298;
t239 = Icges(5,1) * t280 - Icges(5,4) * t279 + Icges(5,5) * t296;
t238 = Icges(5,4) * t282 - Icges(5,2) * t281 + Icges(5,6) * t298;
t237 = Icges(5,4) * t280 - Icges(5,2) * t279 + Icges(5,6) * t296;
t236 = Icges(5,5) * t282 - Icges(5,6) * t281 + Icges(5,3) * t298;
t235 = Icges(5,5) * t280 - Icges(5,6) * t279 + Icges(5,3) * t296;
t233 = rSges(6,1) * t276 - rSges(6,2) * t275 + rSges(6,3) * t293;
t222 = rSges(6,1) * t252 - rSges(6,2) * t251 + rSges(6,3) * t281;
t220 = rSges(6,1) * t250 - rSges(6,2) * t249 + rSges(6,3) * t279;
t206 = -t264 * t328 + t286 * t311 + t383;
t205 = t328 * t265 - t312 * t286 + t379;
t204 = t264 * t312 - t265 * t311 + t382;
t203 = -t241 * t290 + t257 * t277 + t381;
t202 = t290 * t242 - t278 * t257 + t376;
t201 = t241 * t278 - t242 * t277 + t380;
t200 = -t220 * t269 + t233 * t247 + t378;
t199 = t269 * t222 - t248 * t233 + t375;
t198 = t220 * t248 - t222 * t247 + t377;
t197 = qJD(6) * t251 + t247 * t397 - t269 * t399 + t378;
t196 = qJD(6) * t249 - t248 * t397 + t269 * t398 + t375;
t195 = qJD(6) * t275 - t247 * t398 + t248 * t399 + t377;
t1 = t311 * ((t259 * t321 - t261 * t296 + t263 * t297) * t312 + (t258 * t321 - t260 * t296 + t262 * t297) * t311 + (t283 * t321 - t284 * t296 + t285 * t297) * t328) / 0.2e1 + t312 * ((t259 * t322 - t261 * t298 + t263 * t299) * t312 + (t258 * t322 - t260 * t298 + t262 * t299) * t311 + (t283 * t322 - t284 * t298 + t285 * t299) * t328) / 0.2e1 + t328 * ((t259 * t338 - t261 * t319 + t263 * t320) * t312 + (t258 * t338 - t260 * t319 + t262 * t320) * t311 + (t283 * t338 - t284 * t319 + t285 * t320) * t328) / 0.2e1 + m(2) * (t324 ^ 2 + t326 ^ 2 + t327 ^ 2) / 0.2e1 + t290 * ((t236 * t319 - t238 * t293 + t240 * t294) * t278 + (t235 * t319 - t237 * t293 + t239 * t294) * t277 + (t253 * t319 - t254 * t293 + t255 * t294) * t290) / 0.2e1 + t277 * ((t236 * t296 - t238 * t279 + t240 * t280) * t278 + (t235 * t296 - t237 * t279 + t239 * t280) * t277 + (t253 * t296 - t254 * t279 + t255 * t280) * t290) / 0.2e1 + t278 * ((t236 * t298 - t238 * t281 + t240 * t282) * t278 + (t235 * t298 - t237 * t281 + t239 * t282) * t277 + (t253 * t298 - t254 * t281 + t255 * t282) * t290) / 0.2e1 + m(3) * (t267 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + m(1) * (t348 ^ 2 + t349 ^ 2 + t350 ^ 2) / 0.2e1 + m(4) * (t204 ^ 2 + t205 ^ 2 + t206 ^ 2) / 0.2e1 + m(5) * (t201 ^ 2 + t202 ^ 2 + t203 ^ 2) / 0.2e1 + m(6) * (t198 ^ 2 + t199 ^ 2 + t200 ^ 2) / 0.2e1 + m(7) * (t195 ^ 2 + t196 ^ 2 + t197 ^ 2) / 0.2e1 + ((t249 * t419 + t250 * t417 + t279 * t418) * t269 + (t249 * t424 + t250 * t420 + t279 * t422) * t248 + (t425 * t249 + t421 * t250 + t423 * t279) * t247) * t247 / 0.2e1 + ((t251 * t419 + t252 * t417 + t281 * t418) * t269 + (t424 * t251 + t420 * t252 + t422 * t281) * t248 + (t251 * t425 + t421 * t252 + t423 * t281) * t247) * t248 / 0.2e1 + ((t419 * t275 + t417 * t276 + t418 * t293) * t269 + (t275 * t424 + t276 * t420 + t293 * t422) * t248 + (t275 * t425 + t421 * t276 + t423 * t293) * t247) * t269 / 0.2e1 + ((t303 * V_base(5) + t304 * V_base(4) + t331 * t363) * t369 + ((t306 * t368 + t308 * t366) * V_base(4) + (t305 * t368 + t307 * t366) * V_base(5) + (t332 * t368 + t333 * t366) * t363) * t367 + Icges(2,3) * t363 + t351 * V_base(5) + t352 * V_base(4)) * t363 / 0.2e1 + ((t331 * t404 + t332 * t341 + t333 * t342 + t352) * t363 + (t303 * t404 + t305 * t341 + t307 * t342 - t373 * t353 + t355 * t374 + Icges(1,4)) * V_base(5) + (t304 * t404 + t306 * t341 + t308 * t342 - t373 * t354 + t356 * t374 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t331 * t403 + t339 * t332 + t340 * t333 + t351) * t363 + (-t303 * t403 + t339 * t305 + t340 * t307 + t353 * t374 + t373 * t355 + Icges(1,2)) * V_base(5) + (-t304 * t403 + t339 * t306 + t340 * t308 + t354 * t374 + t373 * t356 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
