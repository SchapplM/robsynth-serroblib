% Calculate kinetic energy for
% S6PPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PPRRRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPRRRR3_energykin_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRRR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:07:12
% EndTime: 2019-03-08 19:07:16
% DurationCPUTime: 4.06s
% Computational Cost: add. (9598->455), mult. (26575->665), div. (0->0), fcn. (34839->18), ass. (0->199)
t401 = sin(pkin(13));
t405 = cos(pkin(13));
t452 = Icges(2,5) * t405 - Icges(2,6) * t401 + Icges(1,5);
t451 = Icges(2,5) * t401 + Icges(2,6) * t405 + Icges(1,6);
t450 = cos(qJ(4));
t449 = cos(qJ(5));
t448 = cos(pkin(8));
t447 = sin(pkin(8));
t446 = Icges(2,4) * t401;
t407 = cos(pkin(6));
t445 = qJ(2) * t407;
t403 = sin(pkin(6));
t444 = t401 * t403;
t443 = t401 * t407;
t402 = sin(pkin(7));
t442 = t402 * t403;
t441 = t402 * t407;
t440 = t403 * t405;
t406 = cos(pkin(7));
t439 = t403 * t406;
t404 = cos(pkin(14));
t438 = t404 * t406;
t437 = t405 * t407;
t436 = qJD(2) * t403;
t435 = V_base(5) * qJ(1) + V_base(1);
t431 = qJD(1) + V_base(3);
t400 = sin(pkin(14));
t378 = -t400 * t405 - t404 * t443;
t363 = -t378 * t402 + t401 * t439;
t356 = qJD(3) * t363 + V_base(4);
t376 = -t400 * t401 + t404 * t437;
t362 = -t376 * t402 - t405 * t439;
t355 = qJD(3) * t362 + V_base(5);
t375 = -t404 * t442 + t406 * t407;
t369 = qJD(3) * t375 + V_base(6);
t430 = -qJ(1) - t445;
t379 = -t400 * t443 + t404 * t405;
t411 = sin(qJ(3));
t413 = cos(qJ(3));
t424 = t378 * t406 + t401 * t442;
t337 = -t379 * t411 + t413 * t424;
t324 = -t337 * t447 + t363 * t448;
t312 = qJD(4) * t324 + t356;
t377 = t400 * t437 + t401 * t404;
t425 = t376 * t406 - t402 * t440;
t335 = -t377 * t411 + t413 * t425;
t323 = -t335 * t447 + t362 * t448;
t311 = qJD(4) * t323 + t355;
t359 = t413 * t441 + (-t400 * t411 + t413 * t438) * t403;
t339 = -t359 * t447 + t375 * t448;
t330 = qJD(4) * t339 + t369;
t429 = t448 * t450;
t428 = t450 * t447;
t427 = t401 * t436 + V_base(5) * t445 + t435;
t382 = pkin(1) * t401 - qJ(2) * t440;
t426 = qJD(2) * t407 + V_base(4) * t382 + t431;
t338 = t379 * t413 + t411 * t424;
t410 = sin(qJ(4));
t296 = -t337 * t429 + t338 * t410 - t363 * t428;
t277 = qJD(5) * t296 + t312;
t336 = t377 * t413 + t411 * t425;
t294 = -t335 * t429 + t336 * t410 - t362 * t428;
t276 = qJD(5) * t294 + t311;
t360 = t411 * t441 + (t400 * t413 + t411 * t438) * t403;
t321 = -t359 * t429 + t360 * t410 - t375 * t428;
t290 = qJD(5) * t321 + t330;
t383 = pkin(1) * t405 + qJ(2) * t444;
t423 = V_base(6) * t383 - t405 * t436 + V_base(2);
t342 = pkin(2) * t377 + pkin(9) * t362;
t365 = pkin(2) * t400 * t403 + pkin(9) * t375;
t422 = V_base(5) * t365 + (-t342 - t382) * V_base(6) + t427;
t343 = pkin(2) * t379 + pkin(9) * t363;
t421 = V_base(4) * t342 + (-t343 - t383) * V_base(5) + t426;
t300 = t336 * pkin(3) + pkin(10) * t323;
t325 = t360 * pkin(3) + pkin(10) * t339;
t420 = -t300 * t369 + t355 * t325 + t422;
t301 = t338 * pkin(3) + pkin(10) * t324;
t419 = t356 * t300 - t301 * t355 + t421;
t418 = V_base(6) * t343 + (-t365 + t430) * V_base(4) + t423;
t295 = t336 * t450 + (t335 * t448 + t362 * t447) * t410;
t273 = pkin(4) * t295 + pkin(11) * t294;
t322 = t360 * t450 + (t359 * t448 + t375 * t447) * t410;
t288 = pkin(4) * t322 + pkin(11) * t321;
t417 = -t273 * t330 + t311 * t288 + t420;
t297 = t338 * t450 + (t337 * t448 + t363 * t447) * t410;
t274 = pkin(4) * t297 + pkin(11) * t296;
t416 = t312 * t273 - t274 * t311 + t419;
t415 = t369 * t301 - t325 * t356 + t418;
t414 = t330 * t274 - t288 * t312 + t415;
t412 = cos(qJ(6));
t409 = sin(qJ(5));
t408 = sin(qJ(6));
t398 = Icges(2,4) * t405;
t394 = rSges(2,1) * t405 - rSges(2,2) * t401;
t393 = rSges(2,1) * t401 + rSges(2,2) * t405;
t392 = Icges(2,1) * t405 - t446;
t391 = Icges(2,1) * t401 + t398;
t390 = -Icges(2,2) * t401 + t398;
t389 = Icges(2,2) * t405 + t446;
t386 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t385 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t384 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t373 = rSges(3,3) * t407 + (rSges(3,1) * t400 + rSges(3,2) * t404) * t403;
t372 = Icges(3,5) * t407 + (Icges(3,1) * t400 + Icges(3,4) * t404) * t403;
t371 = Icges(3,6) * t407 + (Icges(3,4) * t400 + Icges(3,2) * t404) * t403;
t370 = Icges(3,3) * t407 + (Icges(3,5) * t400 + Icges(3,6) * t404) * t403;
t367 = V_base(5) * rSges(2,3) - t393 * V_base(6) + t435;
t366 = t394 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t361 = t393 * V_base(4) - t394 * V_base(5) + t431;
t351 = rSges(3,1) * t379 + rSges(3,2) * t378 + rSges(3,3) * t444;
t350 = rSges(3,1) * t377 + rSges(3,2) * t376 - rSges(3,3) * t440;
t349 = Icges(3,1) * t379 + Icges(3,4) * t378 + Icges(3,5) * t444;
t348 = Icges(3,1) * t377 + Icges(3,4) * t376 - Icges(3,5) * t440;
t347 = Icges(3,4) * t379 + Icges(3,2) * t378 + Icges(3,6) * t444;
t346 = Icges(3,4) * t377 + Icges(3,2) * t376 - Icges(3,6) * t440;
t345 = Icges(3,5) * t379 + Icges(3,6) * t378 + Icges(3,3) * t444;
t344 = Icges(3,5) * t377 + Icges(3,6) * t376 - Icges(3,3) * t440;
t329 = rSges(4,1) * t360 + rSges(4,2) * t359 + rSges(4,3) * t375;
t328 = Icges(4,1) * t360 + Icges(4,4) * t359 + Icges(4,5) * t375;
t327 = Icges(4,4) * t360 + Icges(4,2) * t359 + Icges(4,6) * t375;
t326 = Icges(4,5) * t360 + Icges(4,6) * t359 + Icges(4,3) * t375;
t315 = t373 * V_base(5) + (-t350 - t382) * V_base(6) + t427;
t314 = t351 * V_base(6) + (-t373 + t430) * V_base(4) + t423;
t310 = t350 * V_base(4) + (-t351 - t383) * V_base(5) + t426;
t309 = rSges(4,1) * t338 + rSges(4,2) * t337 + rSges(4,3) * t363;
t308 = rSges(4,1) * t336 + rSges(4,2) * t335 + rSges(4,3) * t362;
t307 = Icges(4,1) * t338 + Icges(4,4) * t337 + Icges(4,5) * t363;
t306 = Icges(4,1) * t336 + Icges(4,4) * t335 + Icges(4,5) * t362;
t305 = Icges(4,4) * t338 + Icges(4,2) * t337 + Icges(4,6) * t363;
t304 = Icges(4,4) * t336 + Icges(4,2) * t335 + Icges(4,6) * t362;
t303 = Icges(4,5) * t338 + Icges(4,6) * t337 + Icges(4,3) * t363;
t302 = Icges(4,5) * t336 + Icges(4,6) * t335 + Icges(4,3) * t362;
t299 = t322 * t449 + t339 * t409;
t298 = t322 * t409 - t339 * t449;
t287 = rSges(5,1) * t322 - rSges(5,2) * t321 + rSges(5,3) * t339;
t286 = Icges(5,1) * t322 - Icges(5,4) * t321 + Icges(5,5) * t339;
t285 = Icges(5,4) * t322 - Icges(5,2) * t321 + Icges(5,6) * t339;
t284 = Icges(5,5) * t322 - Icges(5,6) * t321 + Icges(5,3) * t339;
t283 = t297 * t449 + t324 * t409;
t282 = t297 * t409 - t324 * t449;
t281 = t295 * t449 + t323 * t409;
t280 = t295 * t409 - t323 * t449;
t279 = t299 * t412 + t321 * t408;
t278 = -t299 * t408 + t321 * t412;
t275 = pkin(5) * t299 + pkin(12) * t298;
t271 = qJD(6) * t298 + t290;
t269 = -t308 * t369 + t329 * t355 + t422;
t268 = t309 * t369 - t329 * t356 + t418;
t267 = rSges(5,1) * t297 - rSges(5,2) * t296 + rSges(5,3) * t324;
t266 = rSges(5,1) * t295 - rSges(5,2) * t294 + rSges(5,3) * t323;
t265 = Icges(5,1) * t297 - Icges(5,4) * t296 + Icges(5,5) * t324;
t264 = Icges(5,1) * t295 - Icges(5,4) * t294 + Icges(5,5) * t323;
t263 = Icges(5,4) * t297 - Icges(5,2) * t296 + Icges(5,6) * t324;
t262 = Icges(5,4) * t295 - Icges(5,2) * t294 + Icges(5,6) * t323;
t261 = Icges(5,5) * t297 - Icges(5,6) * t296 + Icges(5,3) * t324;
t260 = Icges(5,5) * t295 - Icges(5,6) * t294 + Icges(5,3) * t323;
t259 = rSges(6,1) * t299 - rSges(6,2) * t298 + rSges(6,3) * t321;
t258 = Icges(6,1) * t299 - Icges(6,4) * t298 + Icges(6,5) * t321;
t257 = Icges(6,4) * t299 - Icges(6,2) * t298 + Icges(6,6) * t321;
t256 = Icges(6,5) * t299 - Icges(6,6) * t298 + Icges(6,3) * t321;
t255 = t283 * t412 + t296 * t408;
t254 = -t283 * t408 + t296 * t412;
t253 = t281 * t412 + t294 * t408;
t252 = -t281 * t408 + t294 * t412;
t250 = t308 * t356 - t309 * t355 + t421;
t249 = pkin(5) * t283 + pkin(12) * t282;
t248 = pkin(5) * t281 + pkin(12) * t280;
t247 = qJD(6) * t282 + t277;
t246 = qJD(6) * t280 + t276;
t245 = rSges(6,1) * t283 - rSges(6,2) * t282 + rSges(6,3) * t296;
t244 = rSges(6,1) * t281 - rSges(6,2) * t280 + rSges(6,3) * t294;
t243 = Icges(6,1) * t283 - Icges(6,4) * t282 + Icges(6,5) * t296;
t242 = Icges(6,1) * t281 - Icges(6,4) * t280 + Icges(6,5) * t294;
t241 = Icges(6,4) * t283 - Icges(6,2) * t282 + Icges(6,6) * t296;
t240 = Icges(6,4) * t281 - Icges(6,2) * t280 + Icges(6,6) * t294;
t239 = Icges(6,5) * t283 - Icges(6,6) * t282 + Icges(6,3) * t296;
t238 = Icges(6,5) * t281 - Icges(6,6) * t280 + Icges(6,3) * t294;
t237 = rSges(7,1) * t279 + rSges(7,2) * t278 + rSges(7,3) * t298;
t236 = Icges(7,1) * t279 + Icges(7,4) * t278 + Icges(7,5) * t298;
t235 = Icges(7,4) * t279 + Icges(7,2) * t278 + Icges(7,6) * t298;
t234 = Icges(7,5) * t279 + Icges(7,6) * t278 + Icges(7,3) * t298;
t233 = rSges(7,1) * t255 + rSges(7,2) * t254 + rSges(7,3) * t282;
t232 = rSges(7,1) * t253 + rSges(7,2) * t252 + rSges(7,3) * t280;
t231 = Icges(7,1) * t255 + Icges(7,4) * t254 + Icges(7,5) * t282;
t230 = Icges(7,1) * t253 + Icges(7,4) * t252 + Icges(7,5) * t280;
t229 = Icges(7,4) * t255 + Icges(7,2) * t254 + Icges(7,6) * t282;
t228 = Icges(7,4) * t253 + Icges(7,2) * t252 + Icges(7,6) * t280;
t227 = Icges(7,5) * t255 + Icges(7,6) * t254 + Icges(7,3) * t282;
t226 = Icges(7,5) * t253 + Icges(7,6) * t252 + Icges(7,3) * t280;
t225 = -t266 * t330 + t287 * t311 + t420;
t224 = t267 * t330 - t287 * t312 + t415;
t223 = t266 * t312 - t267 * t311 + t419;
t222 = -t244 * t290 + t259 * t276 + t417;
t221 = t245 * t290 - t259 * t277 + t414;
t220 = t244 * t277 - t245 * t276 + t416;
t219 = -t232 * t271 + t237 * t246 - t248 * t290 + t275 * t276 + t417;
t218 = t233 * t271 - t237 * t247 + t249 * t290 - t275 * t277 + t414;
t217 = t232 * t247 - t233 * t246 + t248 * t277 - t249 * t276 + t416;
t1 = m(7) * (t217 ^ 2 + t218 ^ 2 + t219 ^ 2) / 0.2e1 + t247 * ((t282 * t227 + t254 * t229 + t255 * t231) * t247 + (t226 * t282 + t228 * t254 + t230 * t255) * t246 + (t234 * t282 + t235 * t254 + t236 * t255) * t271) / 0.2e1 + t311 * ((t261 * t323 - t263 * t294 + t265 * t295) * t312 + (t260 * t323 - t262 * t294 + t264 * t295) * t311 + (t284 * t323 - t285 * t294 + t286 * t295) * t330) / 0.2e1 + t312 * ((t261 * t324 - t263 * t296 + t265 * t297) * t312 + (t260 * t324 - t262 * t296 + t264 * t297) * t311 + (t284 * t324 - t285 * t296 + t286 * t297) * t330) / 0.2e1 + m(3) * (t310 ^ 2 + t314 ^ 2 + t315 ^ 2) / 0.2e1 + t290 * ((t239 * t321 - t241 * t298 + t243 * t299) * t277 + (t238 * t321 - t240 * t298 + t242 * t299) * t276 + (t256 * t321 - t257 * t298 + t258 * t299) * t290) / 0.2e1 + m(1) * (t384 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + t369 * ((t303 * t375 + t305 * t359 + t307 * t360) * t356 + (t302 * t375 + t304 * t359 + t306 * t360) * t355 + (t326 * t375 + t327 * t359 + t328 * t360) * t369) / 0.2e1 + m(2) * (t361 ^ 2 + t366 ^ 2 + t367 ^ 2) / 0.2e1 + t355 * ((t303 * t362 + t305 * t335 + t307 * t336) * t356 + (t302 * t362 + t304 * t335 + t306 * t336) * t355 + (t326 * t362 + t327 * t335 + t328 * t336) * t369) / 0.2e1 + t356 * ((t303 * t363 + t305 * t337 + t307 * t338) * t356 + (t302 * t363 + t304 * t337 + t306 * t338) * t355 + (t326 * t363 + t327 * t337 + t328 * t338) * t369) / 0.2e1 + t246 * ((t227 * t280 + t229 * t252 + t231 * t253) * t247 + (t280 * t226 + t252 * t228 + t253 * t230) * t246 + (t234 * t280 + t235 * t252 + t236 * t253) * t271) / 0.2e1 + m(6) * (t220 ^ 2 + t221 ^ 2 + t222 ^ 2) / 0.2e1 + t277 * ((t296 * t239 - t282 * t241 + t283 * t243) * t277 + (t238 * t296 - t240 * t282 + t242 * t283) * t276 + (t256 * t296 - t257 * t282 + t258 * t283) * t290) / 0.2e1 + t271 * ((t227 * t298 + t229 * t278 + t231 * t279) * t247 + (t226 * t298 + t228 * t278 + t230 * t279) * t246 + (t298 * t234 + t278 * t235 + t279 * t236) * t271) / 0.2e1 + t276 * ((t239 * t294 - t241 * t280 + t243 * t281) * t277 + (t294 * t238 - t280 * t240 + t281 * t242) * t276 + (t256 * t294 - t257 * t280 + t258 * t281) * t290) / 0.2e1 + m(5) * (t223 ^ 2 + t224 ^ 2 + t225 ^ 2) / 0.2e1 + t330 * ((t261 * t339 - t263 * t321 + t265 * t322) * t312 + (t260 * t339 - t262 * t321 + t264 * t322) * t311 + (t284 * t339 - t285 * t321 + t286 * t322) * t330) / 0.2e1 + m(4) * (t250 ^ 2 + t268 ^ 2 + t269 ^ 2) / 0.2e1 + ((t370 * t444 + t371 * t378 + t372 * t379 + t452) * V_base(6) + (t344 * t444 + t346 * t378 + t348 * t379 - t389 * t401 + t391 * t405 + Icges(1,4)) * V_base(5) + (t345 * t444 + t347 * t378 + t349 * t379 - t390 * t401 + t392 * t405 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t370 * t440 + t371 * t376 + t372 * t377 + t451) * V_base(6) + (-t344 * t440 + t346 * t376 + t348 * t377 + t389 * t405 + t391 * t401 + Icges(1,2)) * V_base(5) + (-t345 * t440 + t347 * t376 + t349 * t377 + t390 * t405 + t392 * t401 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t370 * t407 + (t371 * t404 + t372 * t400) * t403 + Icges(2,3) + Icges(1,3)) * V_base(6) + (t344 * t407 + (t346 * t404 + t348 * t400) * t403 + t451) * V_base(5) + (t345 * t407 + (t347 * t404 + t349 * t400) * t403 + t452) * V_base(4)) * V_base(6) / 0.2e1;
T  = t1;
