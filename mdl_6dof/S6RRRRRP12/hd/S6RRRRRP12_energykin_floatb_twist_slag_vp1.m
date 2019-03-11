% Calculate kinetic energy for
% S6RRRRRP12
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
% Datum: 2019-03-10 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP12_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP12_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRP12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRP12_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP12_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP12_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP12_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 02:59:32
% EndTime: 2019-03-10 02:59:37
% DurationCPUTime: 4.90s
% Computational Cost: add. (5625->391), mult. (14696->578), div. (0->0), fcn. (18825->14), ass. (0->179)
t435 = Icges(6,1) + Icges(7,1);
t434 = -Icges(6,4) + Icges(7,5);
t433 = Icges(7,4) + Icges(6,5);
t432 = Icges(6,2) + Icges(7,3);
t431 = Icges(7,2) + Icges(6,3);
t430 = -Icges(6,6) + Icges(7,6);
t429 = rSges(7,1) + pkin(5);
t428 = rSges(7,3) + qJ(6);
t370 = cos(pkin(6));
t376 = cos(qJ(2));
t377 = cos(qJ(1));
t402 = t376 * t377;
t374 = sin(qJ(2));
t375 = sin(qJ(1));
t404 = t375 * t374;
t339 = t370 * t402 - t404;
t403 = t375 * t376;
t405 = t374 * t377;
t340 = t370 * t405 + t403;
t373 = sin(qJ(3));
t369 = sin(pkin(6));
t410 = sin(pkin(7));
t393 = t369 * t410;
t411 = cos(pkin(7));
t415 = cos(qJ(3));
t298 = t340 * t415 + (t339 * t411 - t377 * t393) * t373;
t394 = t369 * t411;
t322 = -t339 * t410 - t377 * t394;
t372 = sin(qJ(4));
t414 = cos(qJ(4));
t280 = t298 * t414 + t322 * t372;
t391 = t415 * t410;
t390 = t369 * t391;
t392 = t411 * t415;
t297 = -t339 * t392 + t340 * t373 + t377 * t390;
t371 = sin(qJ(5));
t413 = cos(qJ(5));
t250 = t280 * t371 - t297 * t413;
t251 = t280 * t413 + t297 * t371;
t279 = t298 * t372 - t322 * t414;
t427 = t250 * t432 + t251 * t434 + t279 * t430;
t341 = -t370 * t403 - t405;
t342 = -t370 * t404 + t402;
t300 = t342 * t415 + (t341 * t411 + t375 * t393) * t373;
t323 = -t341 * t410 + t375 * t394;
t282 = t300 * t414 + t323 * t372;
t299 = -t341 * t392 + t342 * t373 - t375 * t390;
t252 = t282 * t371 - t299 * t413;
t253 = t282 * t413 + t299 * t371;
t281 = t300 * t372 - t323 * t414;
t426 = t252 * t432 + t253 * t434 + t281 * t430;
t425 = t250 * t430 + t251 * t433 + t279 * t431;
t424 = t252 * t430 + t253 * t433 + t281 * t431;
t423 = t434 * t250 + t251 * t435 + t433 * t279;
t422 = t434 * t252 + t253 * t435 + t433 * t281;
t321 = t370 * t410 * t373 + (t373 * t376 * t411 + t374 * t415) * t369;
t338 = t370 * t411 - t376 * t393;
t296 = t321 * t414 + t338 * t372;
t408 = t369 * t374;
t320 = -t369 * t376 * t392 - t370 * t391 + t373 * t408;
t277 = t296 * t371 - t320 * t413;
t278 = t296 * t413 + t320 * t371;
t295 = t321 * t372 - t338 * t414;
t421 = t277 * t432 + t278 * t434 + t295 * t430;
t420 = t277 * t430 + t278 * t433 + t295 * t431;
t419 = t434 * t277 + t278 * t435 + t433 * t295;
t412 = pkin(9) * t370;
t409 = Icges(2,4) * t375;
t407 = t369 * t375;
t406 = t369 * t377;
t401 = rSges(7,2) * t279 + t428 * t250 + t251 * t429;
t400 = rSges(7,2) * t281 + t428 * t252 + t253 * t429;
t399 = rSges(7,2) * t295 + t428 * t277 + t278 * t429;
t398 = qJD(2) * t369;
t397 = V_base(5) * pkin(8) + V_base(1);
t352 = t375 * t398 + V_base(4);
t366 = V_base(6) + qJD(1);
t312 = qJD(3) * t323 + t352;
t353 = qJD(2) * t370 + t366;
t276 = qJD(4) * t299 + t312;
t324 = qJD(3) * t338 + t353;
t351 = -t377 * t398 + V_base(5);
t344 = t375 * pkin(1) - pkin(9) * t406;
t389 = -t344 * t366 + V_base(5) * t412 + t397;
t288 = qJD(4) * t320 + t324;
t345 = pkin(1) * t377 + pkin(9) * t407;
t388 = V_base(4) * t344 - t345 * V_base(5) + V_base(3);
t311 = qJD(3) * t322 + t351;
t275 = qJD(4) * t297 + t311;
t387 = t366 * t345 + V_base(2) + (-pkin(8) - t412) * V_base(4);
t301 = t340 * pkin(2) + pkin(10) * t322;
t328 = pkin(2) * t408 + pkin(10) * t338;
t386 = -t301 * t353 + t351 * t328 + t389;
t302 = t342 * pkin(2) + pkin(10) * t323;
t385 = t352 * t301 - t302 * t351 + t388;
t384 = t353 * t302 - t328 * t352 + t387;
t272 = pkin(3) * t298 + pkin(11) * t297;
t287 = pkin(3) * t321 + pkin(11) * t320;
t383 = -t272 * t324 + t311 * t287 + t386;
t273 = pkin(3) * t300 + pkin(11) * t299;
t382 = t312 * t272 - t273 * t311 + t385;
t381 = t324 * t273 - t287 * t312 + t384;
t247 = pkin(4) * t280 + pkin(12) * t279;
t271 = pkin(4) * t296 + pkin(12) * t295;
t380 = -t247 * t288 + t275 * t271 + t383;
t248 = pkin(4) * t282 + pkin(12) * t281;
t379 = t276 * t247 - t248 * t275 + t382;
t378 = t288 * t248 - t271 * t276 + t381;
t367 = Icges(2,4) * t377;
t361 = rSges(2,1) * t377 - t375 * rSges(2,2);
t360 = t375 * rSges(2,1) + rSges(2,2) * t377;
t359 = Icges(2,1) * t377 - t409;
t358 = Icges(2,1) * t375 + t367;
t357 = -Icges(2,2) * t375 + t367;
t356 = Icges(2,2) * t377 + t409;
t350 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t349 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t348 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t334 = rSges(3,3) * t370 + (rSges(3,1) * t374 + rSges(3,2) * t376) * t369;
t333 = Icges(3,5) * t370 + (Icges(3,1) * t374 + Icges(3,4) * t376) * t369;
t332 = Icges(3,6) * t370 + (Icges(3,4) * t374 + Icges(3,2) * t376) * t369;
t331 = Icges(3,3) * t370 + (Icges(3,5) * t374 + Icges(3,6) * t376) * t369;
t327 = V_base(5) * rSges(2,3) - t360 * t366 + t397;
t326 = t361 * t366 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t325 = t360 * V_base(4) - t361 * V_base(5) + V_base(3);
t310 = rSges(3,1) * t342 + rSges(3,2) * t341 + rSges(3,3) * t407;
t309 = t340 * rSges(3,1) + t339 * rSges(3,2) - rSges(3,3) * t406;
t308 = Icges(3,1) * t342 + Icges(3,4) * t341 + Icges(3,5) * t407;
t307 = Icges(3,1) * t340 + Icges(3,4) * t339 - Icges(3,5) * t406;
t306 = Icges(3,4) * t342 + Icges(3,2) * t341 + Icges(3,6) * t407;
t305 = Icges(3,4) * t340 + Icges(3,2) * t339 - Icges(3,6) * t406;
t304 = Icges(3,5) * t342 + Icges(3,6) * t341 + Icges(3,3) * t407;
t303 = Icges(3,5) * t340 + Icges(3,6) * t339 - Icges(3,3) * t406;
t286 = rSges(4,1) * t321 - rSges(4,2) * t320 + rSges(4,3) * t338;
t285 = Icges(4,1) * t321 - Icges(4,4) * t320 + Icges(4,5) * t338;
t284 = Icges(4,4) * t321 - Icges(4,2) * t320 + Icges(4,6) * t338;
t283 = Icges(4,5) * t321 - Icges(4,6) * t320 + Icges(4,3) * t338;
t270 = -t309 * t353 + t334 * t351 + t389;
t269 = t310 * t353 - t334 * t352 + t387;
t268 = qJD(5) * t295 + t288;
t266 = rSges(4,1) * t300 - rSges(4,2) * t299 + rSges(4,3) * t323;
t265 = rSges(4,1) * t298 - rSges(4,2) * t297 + rSges(4,3) * t322;
t264 = t309 * t352 - t310 * t351 + t388;
t263 = Icges(4,1) * t300 - Icges(4,4) * t299 + Icges(4,5) * t323;
t262 = Icges(4,1) * t298 - Icges(4,4) * t297 + Icges(4,5) * t322;
t261 = Icges(4,4) * t300 - Icges(4,2) * t299 + Icges(4,6) * t323;
t260 = Icges(4,4) * t298 - Icges(4,2) * t297 + Icges(4,6) * t322;
t259 = Icges(4,5) * t300 - Icges(4,6) * t299 + Icges(4,3) * t323;
t258 = Icges(4,5) * t298 - Icges(4,6) * t297 + Icges(4,3) * t322;
t257 = rSges(5,1) * t296 - rSges(5,2) * t295 + rSges(5,3) * t320;
t256 = Icges(5,1) * t296 - Icges(5,4) * t295 + Icges(5,5) * t320;
t255 = Icges(5,4) * t296 - Icges(5,2) * t295 + Icges(5,6) * t320;
t254 = Icges(5,5) * t296 - Icges(5,6) * t295 + Icges(5,3) * t320;
t246 = qJD(5) * t281 + t276;
t245 = qJD(5) * t279 + t275;
t242 = rSges(5,1) * t282 - rSges(5,2) * t281 + rSges(5,3) * t299;
t241 = rSges(5,1) * t280 - rSges(5,2) * t279 + rSges(5,3) * t297;
t240 = Icges(5,1) * t282 - Icges(5,4) * t281 + Icges(5,5) * t299;
t239 = Icges(5,1) * t280 - Icges(5,4) * t279 + Icges(5,5) * t297;
t238 = Icges(5,4) * t282 - Icges(5,2) * t281 + Icges(5,6) * t299;
t237 = Icges(5,4) * t280 - Icges(5,2) * t279 + Icges(5,6) * t297;
t236 = Icges(5,5) * t282 - Icges(5,6) * t281 + Icges(5,3) * t299;
t235 = Icges(5,5) * t280 - Icges(5,6) * t279 + Icges(5,3) * t297;
t233 = rSges(6,1) * t278 - rSges(6,2) * t277 + rSges(6,3) * t295;
t222 = rSges(6,1) * t253 - rSges(6,2) * t252 + rSges(6,3) * t281;
t220 = rSges(6,1) * t251 - rSges(6,2) * t250 + rSges(6,3) * t279;
t206 = -t265 * t324 + t286 * t311 + t386;
t205 = t266 * t324 - t286 * t312 + t384;
t204 = t265 * t312 - t266 * t311 + t385;
t203 = -t241 * t288 + t257 * t275 + t383;
t202 = t242 * t288 - t257 * t276 + t381;
t201 = t241 * t276 - t242 * t275 + t382;
t200 = -t220 * t268 + t233 * t245 + t380;
t199 = t222 * t268 - t233 * t246 + t378;
t198 = t220 * t246 - t222 * t245 + t379;
t197 = qJD(6) * t252 + t245 * t399 - t268 * t401 + t380;
t196 = qJD(6) * t250 - t246 * t399 + t268 * t400 + t378;
t195 = qJD(6) * t277 - t245 * t400 + t246 * t401 + t379;
t1 = m(7) * (t195 ^ 2 + t196 ^ 2 + t197 ^ 2) / 0.2e1 + t276 * ((t236 * t299 - t238 * t281 + t240 * t282) * t276 + (t235 * t299 - t237 * t281 + t239 * t282) * t275 + (t254 * t299 - t255 * t281 + t256 * t282) * t288) / 0.2e1 + t275 * ((t236 * t297 - t238 * t279 + t240 * t280) * t276 + (t235 * t297 - t237 * t279 + t239 * t280) * t275 + (t254 * t297 - t255 * t279 + t256 * t280) * t288) / 0.2e1 + m(5) * (t201 ^ 2 + t202 ^ 2 + t203 ^ 2) / 0.2e1 + m(3) * (t264 ^ 2 + t269 ^ 2 + t270 ^ 2) / 0.2e1 + m(6) * (t198 ^ 2 + t199 ^ 2 + t200 ^ 2) / 0.2e1 + t324 * ((t259 * t338 - t261 * t320 + t263 * t321) * t312 + (t258 * t338 - t260 * t320 + t262 * t321) * t311 + (t283 * t338 - t320 * t284 + t321 * t285) * t324) / 0.2e1 + t353 * ((t303 * t351 + t304 * t352 + t331 * t353) * t370 + ((t306 * t376 + t308 * t374) * t352 + (t305 * t376 + t307 * t374) * t351 + (t332 * t376 + t333 * t374) * t353) * t369) / 0.2e1 + t352 * ((t304 * t407 + t306 * t341 + t308 * t342) * t352 + (t303 * t407 + t305 * t341 + t307 * t342) * t351 + (t331 * t407 + t332 * t341 + t333 * t342) * t353) / 0.2e1 + t351 * ((-t304 * t406 + t339 * t306 + t340 * t308) * t352 + (-t303 * t406 + t339 * t305 + t340 * t307) * t351 + (-t331 * t406 + t339 * t332 + t340 * t333) * t353) / 0.2e1 + t311 * ((t259 * t322 - t261 * t297 + t263 * t298) * t312 + (t258 * t322 - t260 * t297 + t262 * t298) * t311 + (t283 * t322 - t284 * t297 + t285 * t298) * t324) / 0.2e1 + t312 * ((t259 * t323 - t261 * t299 + t263 * t300) * t312 + (t258 * t323 - t260 * t299 + t262 * t300) * t311 + (t283 * t323 - t284 * t299 + t285 * t300) * t324) / 0.2e1 + m(2) * (t325 ^ 2 + t326 ^ 2 + t327 ^ 2) / 0.2e1 + t288 * ((t236 * t320 - t238 * t295 + t240 * t296) * t276 + (t235 * t320 - t237 * t295 + t239 * t296) * t275 + (t254 * t320 - t255 * t295 + t256 * t296) * t288) / 0.2e1 + m(1) * (t348 ^ 2 + t349 ^ 2 + t350 ^ 2) / 0.2e1 + m(4) * (t204 ^ 2 + t205 ^ 2 + t206 ^ 2) / 0.2e1 + ((t250 * t421 + t251 * t419 + t279 * t420) * t268 + (t250 * t426 + t251 * t422 + t279 * t424) * t246 + (t250 * t427 + t423 * t251 + t425 * t279) * t245) * t245 / 0.2e1 + ((t252 * t421 + t253 * t419 + t281 * t420) * t268 + (t252 * t426 + t253 * t422 + t281 * t424) * t246 + (t252 * t427 + t423 * t253 + t425 * t281) * t245) * t246 / 0.2e1 + ((t277 * t421 + t278 * t419 + t295 * t420) * t268 + (t277 * t426 + t278 * t422 + t295 * t424) * t246 + (t277 * t427 + t423 * t278 + t425 * t295) * t245) * t268 / 0.2e1 + ((-t375 * t356 + t358 * t377 + Icges(1,4)) * V_base(5) + (-t375 * t357 + t359 * t377 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t356 * t377 + t375 * t358 + Icges(1,2)) * V_base(5) + (t357 * t377 + t375 * t359 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t375 + Icges(2,6) * t377) * V_base(5) + (Icges(2,5) * t377 - Icges(2,6) * t375) * V_base(4) + Icges(2,3) * t366 / 0.2e1) * t366;
T  = t1;
