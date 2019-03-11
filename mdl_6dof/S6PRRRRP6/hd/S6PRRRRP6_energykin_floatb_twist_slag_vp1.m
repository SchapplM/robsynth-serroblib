% Calculate kinetic energy for
% S6PRRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRRP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP6_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:28:28
% EndTime: 2019-03-09 00:28:33
% DurationCPUTime: 5.02s
% Computational Cost: add. (5565->391), mult. (14696->576), div. (0->0), fcn. (18825->14), ass. (0->177)
t434 = Icges(6,1) + Icges(7,1);
t433 = -Icges(6,4) + Icges(7,5);
t432 = Icges(7,4) + Icges(6,5);
t431 = Icges(6,2) + Icges(7,3);
t430 = Icges(7,2) + Icges(6,3);
t429 = -Icges(6,6) + Icges(7,6);
t428 = rSges(7,1) + pkin(5);
t427 = rSges(7,3) + qJ(6);
t366 = sin(pkin(12));
t368 = cos(pkin(12));
t373 = sin(qJ(2));
t369 = cos(pkin(6));
t374 = cos(qJ(2));
t401 = t369 * t374;
t336 = -t366 * t373 + t368 * t401;
t402 = t369 * t373;
t337 = t366 * t374 + t368 * t402;
t372 = sin(qJ(3));
t367 = sin(pkin(6));
t407 = sin(pkin(7));
t390 = t367 * t407;
t408 = cos(pkin(7));
t412 = cos(qJ(3));
t294 = t337 * t412 + (t336 * t408 - t368 * t390) * t372;
t391 = t367 * t408;
t321 = -t336 * t407 - t368 * t391;
t371 = sin(qJ(4));
t411 = cos(qJ(4));
t278 = t294 * t411 + t321 * t371;
t388 = t412 * t407;
t387 = t367 * t388;
t389 = t408 * t412;
t293 = -t336 * t389 + t337 * t372 + t368 * t387;
t370 = sin(qJ(5));
t410 = cos(qJ(5));
t248 = t278 * t370 - t293 * t410;
t249 = t278 * t410 + t293 * t370;
t277 = t294 * t371 - t321 * t411;
t426 = t248 * t431 + t249 * t433 + t277 * t429;
t338 = -t366 * t401 - t368 * t373;
t339 = -t366 * t402 + t368 * t374;
t296 = t339 * t412 + (t338 * t408 + t366 * t390) * t372;
t322 = -t338 * t407 + t366 * t391;
t280 = t296 * t411 + t322 * t371;
t295 = -t338 * t389 + t339 * t372 - t366 * t387;
t250 = t280 * t370 - t295 * t410;
t251 = t280 * t410 + t295 * t370;
t279 = t296 * t371 - t322 * t411;
t425 = t250 * t431 + t251 * t433 + t279 * t429;
t424 = t248 * t429 + t249 * t432 + t277 * t430;
t423 = t250 * t429 + t251 * t432 + t279 * t430;
t422 = t433 * t248 + t249 * t434 + t432 * t277;
t421 = t433 * t250 + t251 * t434 + t432 * t279;
t320 = t369 * t407 * t372 + (t372 * t374 * t408 + t373 * t412) * t367;
t335 = t369 * t408 - t374 * t390;
t298 = t320 * t411 + t335 * t371;
t403 = t367 * t373;
t319 = -t367 * t374 * t389 - t369 * t388 + t372 * t403;
t275 = t298 * t370 - t319 * t410;
t276 = t298 * t410 + t319 * t370;
t297 = t320 * t371 - t335 * t411;
t420 = t275 * t431 + t276 * t433 + t297 * t429;
t419 = t275 * t429 + t276 * t432 + t297 * t430;
t418 = t433 * t275 + t276 * t434 + t432 * t297;
t409 = pkin(8) * t369;
t406 = Icges(2,4) * t366;
t405 = t366 * t367;
t404 = t367 * t368;
t400 = rSges(7,2) * t277 + t427 * t248 + t249 * t428;
t399 = rSges(7,2) * t279 + t427 * t250 + t251 * t428;
t398 = rSges(7,2) * t297 + t427 * t275 + t276 * t428;
t397 = qJD(2) * t367;
t396 = V_base(5) * qJ(1) + V_base(1);
t392 = qJD(1) + V_base(3);
t350 = t366 * t397 + V_base(4);
t360 = qJD(2) * t369 + V_base(6);
t310 = qJD(3) * t322 + t350;
t323 = qJD(3) * t335 + t360;
t274 = qJD(4) * t295 + t310;
t287 = qJD(4) * t319 + t323;
t349 = -t368 * t397 + V_base(5);
t309 = qJD(3) * t321 + t349;
t342 = pkin(1) * t366 - pkin(8) * t404;
t386 = -t342 * V_base(6) + V_base(5) * t409 + t396;
t343 = pkin(1) * t368 + pkin(8) * t405;
t385 = V_base(4) * t342 - t343 * V_base(5) + t392;
t273 = qJD(4) * t293 + t309;
t384 = V_base(6) * t343 + V_base(2) + (-qJ(1) - t409) * V_base(4);
t299 = t337 * pkin(2) + pkin(9) * t321;
t324 = pkin(2) * t403 + pkin(9) * t335;
t383 = -t299 * t360 + t349 * t324 + t386;
t300 = t339 * pkin(2) + pkin(9) * t322;
t382 = t350 * t299 - t300 * t349 + t385;
t381 = t360 * t300 - t324 * t350 + t384;
t269 = pkin(3) * t294 + pkin(10) * t293;
t285 = pkin(3) * t320 + pkin(10) * t319;
t380 = -t269 * t323 + t309 * t285 + t383;
t270 = pkin(3) * t296 + pkin(10) * t295;
t379 = t310 * t269 - t270 * t309 + t382;
t378 = t323 * t270 - t285 * t310 + t381;
t245 = pkin(4) * t278 + pkin(11) * t277;
t271 = pkin(4) * t298 + pkin(11) * t297;
t377 = -t245 * t287 + t273 * t271 + t380;
t246 = pkin(4) * t280 + pkin(11) * t279;
t376 = t274 * t245 - t246 * t273 + t379;
t375 = t287 * t246 - t271 * t274 + t378;
t364 = Icges(2,4) * t368;
t358 = rSges(2,1) * t368 - rSges(2,2) * t366;
t357 = rSges(2,1) * t366 + rSges(2,2) * t368;
t356 = Icges(2,1) * t368 - t406;
t355 = Icges(2,1) * t366 + t364;
t354 = -Icges(2,2) * t366 + t364;
t353 = Icges(2,2) * t368 + t406;
t348 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t347 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t346 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t332 = t369 * rSges(3,3) + (rSges(3,1) * t373 + rSges(3,2) * t374) * t367;
t331 = Icges(3,5) * t369 + (Icges(3,1) * t373 + Icges(3,4) * t374) * t367;
t330 = Icges(3,6) * t369 + (Icges(3,4) * t373 + Icges(3,2) * t374) * t367;
t329 = Icges(3,3) * t369 + (Icges(3,5) * t373 + Icges(3,6) * t374) * t367;
t326 = V_base(5) * rSges(2,3) - t357 * V_base(6) + t396;
t325 = t358 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t318 = t357 * V_base(4) - t358 * V_base(5) + t392;
t308 = rSges(3,1) * t339 + rSges(3,2) * t338 + rSges(3,3) * t405;
t307 = rSges(3,1) * t337 + rSges(3,2) * t336 - rSges(3,3) * t404;
t306 = Icges(3,1) * t339 + Icges(3,4) * t338 + Icges(3,5) * t405;
t305 = Icges(3,1) * t337 + Icges(3,4) * t336 - Icges(3,5) * t404;
t304 = Icges(3,4) * t339 + Icges(3,2) * t338 + Icges(3,6) * t405;
t303 = Icges(3,4) * t337 + Icges(3,2) * t336 - Icges(3,6) * t404;
t302 = Icges(3,5) * t339 + Icges(3,6) * t338 + Icges(3,3) * t405;
t301 = Icges(3,5) * t337 + Icges(3,6) * t336 - Icges(3,3) * t404;
t284 = rSges(4,1) * t320 - rSges(4,2) * t319 + rSges(4,3) * t335;
t283 = Icges(4,1) * t320 - Icges(4,4) * t319 + Icges(4,5) * t335;
t282 = Icges(4,4) * t320 - Icges(4,2) * t319 + Icges(4,6) * t335;
t281 = Icges(4,5) * t320 - Icges(4,6) * t319 + Icges(4,3) * t335;
t268 = -t307 * t360 + t332 * t349 + t386;
t267 = t308 * t360 - t332 * t350 + t384;
t266 = qJD(5) * t297 + t287;
t264 = rSges(5,1) * t298 - rSges(5,2) * t297 + rSges(5,3) * t319;
t263 = rSges(4,1) * t296 - rSges(4,2) * t295 + rSges(4,3) * t322;
t262 = rSges(4,1) * t294 - rSges(4,2) * t293 + rSges(4,3) * t321;
t261 = Icges(5,1) * t298 - Icges(5,4) * t297 + Icges(5,5) * t319;
t260 = Icges(5,4) * t298 - Icges(5,2) * t297 + Icges(5,6) * t319;
t259 = Icges(5,5) * t298 - Icges(5,6) * t297 + Icges(5,3) * t319;
t258 = Icges(4,1) * t296 - Icges(4,4) * t295 + Icges(4,5) * t322;
t257 = Icges(4,1) * t294 - Icges(4,4) * t293 + Icges(4,5) * t321;
t256 = Icges(4,4) * t296 - Icges(4,2) * t295 + Icges(4,6) * t322;
t255 = Icges(4,4) * t294 - Icges(4,2) * t293 + Icges(4,6) * t321;
t254 = Icges(4,5) * t296 - Icges(4,6) * t295 + Icges(4,3) * t322;
t253 = Icges(4,5) * t294 - Icges(4,6) * t293 + Icges(4,3) * t321;
t252 = t307 * t350 - t308 * t349 + t385;
t243 = qJD(5) * t279 + t274;
t242 = qJD(5) * t277 + t273;
t240 = rSges(6,1) * t276 - rSges(6,2) * t275 + rSges(6,3) * t297;
t238 = rSges(5,1) * t280 - rSges(5,2) * t279 + rSges(5,3) * t295;
t237 = rSges(5,1) * t278 - rSges(5,2) * t277 + rSges(5,3) * t293;
t236 = Icges(5,1) * t280 - Icges(5,4) * t279 + Icges(5,5) * t295;
t235 = Icges(5,1) * t278 - Icges(5,4) * t277 + Icges(5,5) * t293;
t232 = Icges(5,4) * t280 - Icges(5,2) * t279 + Icges(5,6) * t295;
t231 = Icges(5,4) * t278 - Icges(5,2) * t277 + Icges(5,6) * t293;
t228 = Icges(5,5) * t280 - Icges(5,6) * t279 + Icges(5,3) * t295;
t227 = Icges(5,5) * t278 - Icges(5,6) * t277 + Icges(5,3) * t293;
t220 = rSges(6,1) * t251 - rSges(6,2) * t250 + rSges(6,3) * t279;
t218 = rSges(6,1) * t249 - rSges(6,2) * t248 + rSges(6,3) * t277;
t204 = -t262 * t323 + t284 * t309 + t383;
t203 = t263 * t323 - t284 * t310 + t381;
t202 = t262 * t310 - t263 * t309 + t382;
t201 = -t237 * t287 + t264 * t273 + t380;
t200 = t238 * t287 - t264 * t274 + t378;
t199 = t237 * t274 - t238 * t273 + t379;
t198 = -t218 * t266 + t240 * t242 + t377;
t197 = t220 * t266 - t240 * t243 + t375;
t196 = t218 * t243 - t220 * t242 + t376;
t195 = qJD(6) * t250 + t242 * t398 - t266 * t400 + t377;
t194 = qJD(6) * t248 - t243 * t398 + t266 * t399 + t375;
t193 = qJD(6) * t275 - t242 * t399 + t243 * t400 + t376;
t1 = t309 * ((t254 * t321 - t256 * t293 + t258 * t294) * t310 + (t253 * t321 - t255 * t293 + t257 * t294) * t309 + (t281 * t321 - t282 * t293 + t283 * t294) * t323) / 0.2e1 + t310 * ((t254 * t322 - t256 * t295 + t258 * t296) * t310 + (t253 * t322 - t255 * t295 + t257 * t296) * t309 + (t281 * t322 - t282 * t295 + t283 * t296) * t323) / 0.2e1 + m(2) * (t318 ^ 2 + t325 ^ 2 + t326 ^ 2) / 0.2e1 + t287 * ((t228 * t319 - t232 * t297 + t236 * t298) * t274 + (t227 * t319 - t231 * t297 + t235 * t298) * t273 + (t259 * t319 - t260 * t297 + t261 * t298) * t287) / 0.2e1 + t323 * ((t254 * t335 - t256 * t319 + t258 * t320) * t310 + (t253 * t335 - t255 * t319 + t257 * t320) * t309 + (t281 * t335 - t282 * t319 + t283 * t320) * t323) / 0.2e1 + m(1) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + m(4) * (t202 ^ 2 + t203 ^ 2 + t204 ^ 2) / 0.2e1 + m(5) * (t199 ^ 2 + t200 ^ 2 + t201 ^ 2) / 0.2e1 + m(6) * (t196 ^ 2 + t197 ^ 2 + t198 ^ 2) / 0.2e1 + t273 * ((t228 * t293 - t232 * t277 + t236 * t278) * t274 + (t227 * t293 - t231 * t277 + t235 * t278) * t273 + (t259 * t293 - t260 * t277 + t261 * t278) * t287) / 0.2e1 + t274 * ((t228 * t295 - t232 * t279 + t236 * t280) * t274 + (t227 * t295 - t231 * t279 + t235 * t280) * t273 + (t259 * t295 - t260 * t279 + t261 * t280) * t287) / 0.2e1 + t360 * ((t301 * t349 + t302 * t350 + t329 * t360) * t369 + ((t304 * t374 + t306 * t373) * t350 + (t303 * t374 + t305 * t373) * t349 + (t330 * t374 + t331 * t373) * t360) * t367) / 0.2e1 + m(7) * (t193 ^ 2 + t194 ^ 2 + t195 ^ 2) / 0.2e1 + m(3) * (t252 ^ 2 + t267 ^ 2 + t268 ^ 2) / 0.2e1 + t349 * ((-t302 * t404 + t304 * t336 + t306 * t337) * t350 + (-t301 * t404 + t303 * t336 + t305 * t337) * t349 + (-t329 * t404 + t330 * t336 + t331 * t337) * t360) / 0.2e1 + t350 * ((t302 * t405 + t304 * t338 + t306 * t339) * t350 + (t301 * t405 + t303 * t338 + t305 * t339) * t349 + (t329 * t405 + t330 * t338 + t331 * t339) * t360) / 0.2e1 + ((t248 * t420 + t249 * t418 + t277 * t419) * t266 + (t248 * t425 + t249 * t421 + t277 * t423) * t243 + (t248 * t426 + t422 * t249 + t424 * t277) * t242) * t242 / 0.2e1 + ((t250 * t420 + t251 * t418 + t279 * t419) * t266 + (t250 * t425 + t251 * t421 + t279 * t423) * t243 + (t250 * t426 + t422 * t251 + t424 * t279) * t242) * t243 / 0.2e1 + ((t275 * t420 + t276 * t418 + t297 * t419) * t266 + (t275 * t425 + t276 * t421 + t297 * t423) * t243 + (t275 * t426 + t422 * t276 + t424 * t297) * t242) * t266 / 0.2e1 + ((-t353 * t366 + t355 * t368 + Icges(1,4)) * V_base(5) + (-t354 * t366 + t356 * t368 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t353 * t368 + t355 * t366 + Icges(1,2)) * V_base(5) + (t354 * t368 + t356 * t366 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t368 - Icges(2,6) * t366 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t366 + Icges(2,6) * t368 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
