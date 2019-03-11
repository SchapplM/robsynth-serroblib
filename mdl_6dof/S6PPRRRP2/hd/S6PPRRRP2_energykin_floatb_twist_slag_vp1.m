% Calculate kinetic energy for
% S6PPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PPRRRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRP2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRRP2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:55:56
% EndTime: 2019-03-08 18:56:00
% DurationCPUTime: 4.07s
% Computational Cost: add. (5430->396), mult. (14471->560), div. (0->0), fcn. (18600->14), ass. (0->177)
t432 = Icges(6,1) + Icges(7,1);
t431 = -Icges(6,4) + Icges(7,5);
t430 = Icges(7,4) + Icges(6,5);
t429 = Icges(6,2) + Icges(7,3);
t428 = Icges(7,2) + Icges(6,3);
t427 = -Icges(6,6) + Icges(7,6);
t426 = rSges(7,1) + pkin(5);
t425 = rSges(7,3) + qJ(6);
t364 = sin(pkin(11));
t367 = cos(pkin(11));
t424 = Icges(2,5) * t367 - Icges(2,6) * t364 + Icges(1,5);
t423 = Icges(2,5) * t364 + Icges(2,6) * t367 + Icges(1,6);
t363 = sin(pkin(12));
t366 = cos(pkin(12));
t368 = cos(pkin(6));
t399 = t367 * t368;
t336 = -t363 * t364 + t366 * t399;
t337 = t363 * t399 + t364 * t366;
t371 = sin(qJ(3));
t365 = sin(pkin(6));
t406 = sin(pkin(7));
t388 = t365 * t406;
t407 = cos(pkin(7));
t410 = cos(qJ(3));
t291 = t337 * t410 + (t336 * t407 - t367 * t388) * t371;
t389 = t365 * t407;
t320 = -t336 * t406 - t367 * t389;
t370 = sin(qJ(4));
t409 = cos(qJ(4));
t276 = t291 * t409 + t320 * t370;
t385 = t410 * t406;
t382 = t365 * t385;
t386 = t407 * t410;
t290 = -t336 * t386 + t337 * t371 + t367 * t382;
t369 = sin(qJ(5));
t408 = cos(qJ(5));
t247 = t276 * t369 - t290 * t408;
t248 = t276 * t408 + t290 * t369;
t275 = t291 * t370 - t320 * t409;
t422 = t429 * t247 + t431 * t248 + t427 * t275;
t401 = t364 * t368;
t338 = -t363 * t367 - t366 * t401;
t339 = -t363 * t401 + t366 * t367;
t293 = t339 * t410 + (t338 * t407 + t364 * t388) * t371;
t321 = -t338 * t406 + t364 * t389;
t278 = t293 * t409 + t321 * t370;
t292 = -t338 * t386 + t339 * t371 - t364 * t382;
t249 = t278 * t369 - t292 * t408;
t250 = t278 * t408 + t292 * t369;
t277 = t293 * t370 - t321 * t409;
t421 = t429 * t249 + t431 * t250 + t427 * t277;
t420 = t427 * t247 + t430 * t248 + t428 * t275;
t419 = t427 * t249 + t430 * t250 + t428 * t277;
t418 = t431 * t247 + t432 * t248 + t430 * t275;
t417 = t431 * t249 + t432 * t250 + t430 * t277;
t318 = t368 * t406 * t371 + (t366 * t371 * t407 + t363 * t410) * t365;
t335 = -t366 * t388 + t368 * t407;
t296 = t318 * t409 + t335 * t370;
t403 = t363 * t365;
t317 = -t365 * t366 * t386 - t368 * t385 + t371 * t403;
t279 = t296 * t369 - t317 * t408;
t280 = t296 * t408 + t317 * t369;
t295 = t318 * t370 - t335 * t409;
t416 = t429 * t279 + t431 * t280 + t427 * t295;
t415 = t427 * t279 + t430 * t280 + t428 * t295;
t414 = t431 * t279 + t432 * t280 + t430 * t295;
t405 = Icges(2,4) * t364;
t404 = qJ(2) * t368;
t402 = t364 * t365;
t400 = t365 * t367;
t398 = rSges(7,2) * t275 + t425 * t247 + t248 * t426;
t397 = rSges(7,2) * t277 + t425 * t249 + t250 * t426;
t396 = rSges(7,2) * t295 + t425 * t279 + t280 * t426;
t395 = qJD(2) * t365;
t394 = V_base(5) * qJ(1) + V_base(1);
t390 = qJD(1) + V_base(3);
t310 = qJD(3) * t321 + V_base(4);
t309 = qJD(3) * t320 + V_base(5);
t328 = qJD(3) * t335 + V_base(6);
t387 = -qJ(1) - t404;
t274 = qJD(4) * t292 + t310;
t273 = qJD(4) * t290 + t309;
t294 = qJD(4) * t317 + t328;
t384 = t364 * t395 + V_base(5) * t404 + t394;
t342 = pkin(1) * t364 - qJ(2) * t400;
t383 = qJD(2) * t368 + V_base(4) * t342 + t390;
t343 = pkin(1) * t367 + qJ(2) * t402;
t381 = V_base(6) * t343 - t367 * t395 + V_base(2);
t299 = t337 * pkin(2) + pkin(8) * t320;
t323 = pkin(2) * t403 + pkin(8) * t335;
t380 = V_base(5) * t323 + (-t299 - t342) * V_base(6) + t384;
t300 = t339 * pkin(2) + pkin(8) * t321;
t379 = V_base(4) * t299 + (-t300 - t343) * V_base(5) + t383;
t266 = pkin(3) * t291 + pkin(9) * t290;
t285 = pkin(3) * t318 + pkin(9) * t317;
t378 = -t266 * t328 + t309 * t285 + t380;
t267 = pkin(3) * t293 + pkin(9) * t292;
t377 = t310 * t266 - t267 * t309 + t379;
t376 = V_base(6) * t300 + (-t323 + t387) * V_base(4) + t381;
t242 = pkin(4) * t276 + pkin(10) * t275;
t268 = pkin(4) * t296 + pkin(10) * t295;
t375 = -t242 * t294 + t273 * t268 + t378;
t243 = pkin(4) * t278 + pkin(10) * t277;
t374 = t274 * t242 - t243 * t273 + t377;
t373 = t328 * t267 - t285 * t310 + t376;
t372 = t294 * t243 - t268 * t274 + t373;
t361 = Icges(2,4) * t367;
t356 = rSges(2,1) * t367 - rSges(2,2) * t364;
t355 = rSges(2,1) * t364 + rSges(2,2) * t367;
t354 = Icges(2,1) * t367 - t405;
t353 = Icges(2,1) * t364 + t361;
t352 = -Icges(2,2) * t364 + t361;
t351 = Icges(2,2) * t367 + t405;
t348 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t347 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t346 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t332 = rSges(3,3) * t368 + (rSges(3,1) * t363 + rSges(3,2) * t366) * t365;
t331 = Icges(3,5) * t368 + (Icges(3,1) * t363 + Icges(3,4) * t366) * t365;
t330 = Icges(3,6) * t368 + (Icges(3,4) * t363 + Icges(3,2) * t366) * t365;
t329 = Icges(3,3) * t368 + (Icges(3,5) * t363 + Icges(3,6) * t366) * t365;
t325 = V_base(5) * rSges(2,3) - t355 * V_base(6) + t394;
t324 = t356 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t319 = t355 * V_base(4) - t356 * V_base(5) + t390;
t308 = rSges(3,1) * t339 + rSges(3,2) * t338 + rSges(3,3) * t402;
t307 = rSges(3,1) * t337 + rSges(3,2) * t336 - rSges(3,3) * t400;
t306 = Icges(3,1) * t339 + Icges(3,4) * t338 + Icges(3,5) * t402;
t305 = Icges(3,1) * t337 + Icges(3,4) * t336 - Icges(3,5) * t400;
t304 = Icges(3,4) * t339 + Icges(3,2) * t338 + Icges(3,6) * t402;
t303 = Icges(3,4) * t337 + Icges(3,2) * t336 - Icges(3,6) * t400;
t302 = Icges(3,5) * t339 + Icges(3,6) * t338 + Icges(3,3) * t402;
t301 = Icges(3,5) * t337 + Icges(3,6) * t336 - Icges(3,3) * t400;
t284 = rSges(4,1) * t318 - rSges(4,2) * t317 + rSges(4,3) * t335;
t283 = Icges(4,1) * t318 - Icges(4,4) * t317 + Icges(4,5) * t335;
t282 = Icges(4,4) * t318 - Icges(4,2) * t317 + Icges(4,6) * t335;
t281 = Icges(4,5) * t318 - Icges(4,6) * t317 + Icges(4,3) * t335;
t271 = t332 * V_base(5) + (-t307 - t342) * V_base(6) + t384;
t270 = t308 * V_base(6) + (-t332 + t387) * V_base(4) + t381;
t269 = qJD(5) * t295 + t294;
t264 = t307 * V_base(4) + (-t308 - t343) * V_base(5) + t383;
t263 = rSges(5,1) * t296 - rSges(5,2) * t295 + rSges(5,3) * t317;
t262 = Icges(5,1) * t296 - Icges(5,4) * t295 + Icges(5,5) * t317;
t261 = Icges(5,4) * t296 - Icges(5,2) * t295 + Icges(5,6) * t317;
t260 = Icges(5,5) * t296 - Icges(5,6) * t295 + Icges(5,3) * t317;
t259 = rSges(4,1) * t293 - rSges(4,2) * t292 + rSges(4,3) * t321;
t258 = rSges(4,1) * t291 - rSges(4,2) * t290 + rSges(4,3) * t320;
t257 = Icges(4,1) * t293 - Icges(4,4) * t292 + Icges(4,5) * t321;
t256 = Icges(4,1) * t291 - Icges(4,4) * t290 + Icges(4,5) * t320;
t255 = Icges(4,4) * t293 - Icges(4,2) * t292 + Icges(4,6) * t321;
t254 = Icges(4,4) * t291 - Icges(4,2) * t290 + Icges(4,6) * t320;
t253 = Icges(4,5) * t293 - Icges(4,6) * t292 + Icges(4,3) * t321;
t252 = Icges(4,5) * t291 - Icges(4,6) * t290 + Icges(4,3) * t320;
t246 = qJD(5) * t277 + t274;
t245 = qJD(5) * t275 + t273;
t239 = rSges(6,1) * t280 - rSges(6,2) * t279 + rSges(6,3) * t295;
t231 = rSges(5,1) * t278 - rSges(5,2) * t277 + rSges(5,3) * t292;
t230 = rSges(5,1) * t276 - rSges(5,2) * t275 + rSges(5,3) * t290;
t229 = Icges(5,1) * t278 - Icges(5,4) * t277 + Icges(5,5) * t292;
t228 = Icges(5,1) * t276 - Icges(5,4) * t275 + Icges(5,5) * t290;
t227 = Icges(5,4) * t278 - Icges(5,2) * t277 + Icges(5,6) * t292;
t226 = Icges(5,4) * t276 - Icges(5,2) * t275 + Icges(5,6) * t290;
t225 = Icges(5,5) * t278 - Icges(5,6) * t277 + Icges(5,3) * t292;
t224 = Icges(5,5) * t276 - Icges(5,6) * t275 + Icges(5,3) * t290;
t220 = rSges(6,1) * t250 - rSges(6,2) * t249 + rSges(6,3) * t277;
t218 = rSges(6,1) * t248 - rSges(6,2) * t247 + rSges(6,3) * t275;
t204 = -t258 * t328 + t284 * t309 + t380;
t203 = t259 * t328 - t284 * t310 + t376;
t202 = t258 * t310 - t259 * t309 + t379;
t201 = -t230 * t294 + t263 * t273 + t378;
t200 = t231 * t294 - t263 * t274 + t373;
t199 = t230 * t274 - t231 * t273 + t377;
t198 = -t218 * t269 + t239 * t245 + t375;
t197 = t220 * t269 - t239 * t246 + t372;
t196 = t218 * t246 - t220 * t245 + t374;
t195 = qJD(6) * t249 + t245 * t396 - t269 * t398 + t375;
t194 = qJD(6) * t247 - t246 * t396 + t269 * t397 + t372;
t193 = qJD(6) * t279 - t245 * t397 + t246 * t398 + t374;
t1 = m(1) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + t328 * ((t253 * t335 - t255 * t317 + t257 * t318) * t310 + (t252 * t335 - t254 * t317 + t256 * t318) * t309 + (t281 * t335 - t282 * t317 + t283 * t318) * t328) / 0.2e1 + m(2) * (t319 ^ 2 + t324 ^ 2 + t325 ^ 2) / 0.2e1 + t309 * ((t253 * t320 - t255 * t290 + t257 * t291) * t310 + (t252 * t320 - t254 * t290 + t256 * t291) * t309 + (t281 * t320 - t282 * t290 + t283 * t291) * t328) / 0.2e1 + t310 * ((t253 * t321 - t255 * t292 + t257 * t293) * t310 + (t252 * t321 - t254 * t292 + t256 * t293) * t309 + (t281 * t321 - t282 * t292 + t283 * t293) * t328) / 0.2e1 + t294 * ((t225 * t317 - t227 * t295 + t229 * t296) * t274 + (t224 * t317 - t226 * t295 + t228 * t296) * t273 + (t260 * t317 - t261 * t295 + t262 * t296) * t294) / 0.2e1 + t274 * ((t225 * t292 - t277 * t227 + t278 * t229) * t274 + (t224 * t292 - t226 * t277 + t228 * t278) * t273 + (t260 * t292 - t261 * t277 + t262 * t278) * t294) / 0.2e1 + t273 * ((t225 * t290 - t227 * t275 + t229 * t276) * t274 + (t224 * t290 - t226 * t275 + t228 * t276) * t273 + (t260 * t290 - t261 * t275 + t262 * t276) * t294) / 0.2e1 + m(3) * (t264 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + m(4) * (t202 ^ 2 + t203 ^ 2 + t204 ^ 2) / 0.2e1 + m(5) * (t199 ^ 2 + t200 ^ 2 + t201 ^ 2) / 0.2e1 + m(6) * (t196 ^ 2 + t197 ^ 2 + t198 ^ 2) / 0.2e1 + m(7) * (t193 ^ 2 + t194 ^ 2 + t195 ^ 2) / 0.2e1 + ((t247 * t416 + t248 * t414 + t275 * t415) * t269 + (t247 * t421 + t248 * t417 + t275 * t419) * t246 + (t422 * t247 + t418 * t248 + t420 * t275) * t245) * t245 / 0.2e1 + ((t249 * t416 + t250 * t414 + t277 * t415) * t269 + (t421 * t249 + t417 * t250 + t419 * t277) * t246 + (t249 * t422 + t250 * t418 + t277 * t420) * t245) * t246 / 0.2e1 + ((t416 * t279 + t414 * t280 + t415 * t295) * t269 + (t279 * t421 + t280 * t417 + t295 * t419) * t246 + (t279 * t422 + t280 * t418 + t295 * t420) * t245) * t269 / 0.2e1 + ((t329 * t402 + t330 * t338 + t331 * t339 + t424) * V_base(6) + (t301 * t402 + t303 * t338 + t305 * t339 - t351 * t364 + t353 * t367 + Icges(1,4)) * V_base(5) + (t302 * t402 + t304 * t338 + t306 * t339 - t352 * t364 + t354 * t367 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t329 * t400 + t330 * t336 + t331 * t337 + t423) * V_base(6) + (-t301 * t400 + t303 * t336 + t305 * t337 + t351 * t367 + t353 * t364 + Icges(1,2)) * V_base(5) + (-t302 * t400 + t304 * t336 + t306 * t337 + t352 * t367 + t354 * t364 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t329 * t368 + (t330 * t366 + t331 * t363) * t365 + Icges(2,3) + Icges(1,3)) * V_base(6) + (t301 * t368 + (t303 * t366 + t305 * t363) * t365 + t423) * V_base(5) + (t302 * t368 + (t304 * t366 + t306 * t363) * t365 + t424) * V_base(4)) * V_base(6) / 0.2e1;
T  = t1;
