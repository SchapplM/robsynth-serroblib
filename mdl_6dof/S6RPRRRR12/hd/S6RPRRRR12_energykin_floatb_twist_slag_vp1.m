% Calculate kinetic energy for
% S6RPRRRR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 08:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR12_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR12_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRR12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RPRRRR12_energykin_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR12_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR12_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR12_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:50:19
% EndTime: 2019-03-09 07:50:24
% DurationCPUTime: 5.15s
% Computational Cost: add. (9658->455), mult. (26575->670), div. (0->0), fcn. (34839->18), ass. (0->198)
t449 = cos(qJ(4));
t448 = cos(qJ(5));
t447 = cos(pkin(8));
t446 = sin(pkin(8));
t412 = sin(qJ(1));
t445 = Icges(2,4) * t412;
t407 = cos(pkin(6));
t444 = qJ(2) * t407;
t403 = sin(pkin(7));
t443 = t403 * t407;
t404 = sin(pkin(6));
t442 = t404 * t412;
t415 = cos(qJ(1));
t441 = t404 * t415;
t405 = cos(pkin(14));
t406 = cos(pkin(7));
t440 = t405 * t406;
t439 = t407 * t415;
t402 = sin(pkin(14));
t438 = t412 * t402;
t437 = t412 * t405;
t436 = qJD(2) * t404;
t435 = V_base(5) * pkin(9) + V_base(1);
t380 = -t402 * t415 - t407 * t437;
t363 = -t380 * t403 + t406 * t442;
t357 = qJD(3) * t363 + V_base(4);
t378 = t405 * t439 - t438;
t362 = -t378 * t403 - t406 * t441;
t356 = qJD(3) * t362 + V_base(5);
t399 = V_base(6) + qJD(1);
t432 = -pkin(9) - t444;
t383 = pkin(1) * t412 - qJ(2) * t441;
t431 = qJD(2) * t407 + t383 * V_base(4) + V_base(3);
t381 = t405 * t415 - t407 * t438;
t411 = sin(qJ(3));
t414 = cos(qJ(3));
t426 = t380 * t406 + t403 * t442;
t340 = -t381 * t411 + t414 * t426;
t325 = -t340 * t446 + t363 * t447;
t313 = qJD(4) * t325 + t357;
t379 = t402 * t439 + t437;
t427 = t378 * t406 - t403 * t441;
t338 = -t379 * t411 + t414 * t427;
t324 = -t338 * t446 + t362 * t447;
t312 = qJD(4) * t324 + t356;
t377 = -t403 * t404 * t405 + t406 * t407;
t369 = qJD(3) * t377 + t399;
t430 = t447 * t449;
t429 = t449 * t446;
t428 = t412 * t436 + t444 * V_base(5) + t435;
t341 = t381 * t414 + t411 * t426;
t410 = sin(qJ(4));
t299 = -t340 * t430 + t341 * t410 - t363 * t429;
t278 = qJD(5) * t299 + t313;
t339 = t379 * t414 + t411 * t427;
t297 = -t338 * t430 + t339 * t410 - t362 * t429;
t277 = qJD(5) * t297 + t312;
t360 = t414 * t443 + (-t402 * t411 + t414 * t440) * t404;
t336 = -t360 * t446 + t377 * t447;
t331 = qJD(4) * t336 + t369;
t361 = t411 * t443 + (t402 * t414 + t411 * t440) * t404;
t318 = -t360 * t430 + t361 * t410 - t377 * t429;
t291 = qJD(5) * t318 + t331;
t384 = pkin(1) * t415 + qJ(2) * t442;
t425 = t384 * t399 - t415 * t436 + V_base(2);
t343 = pkin(2) * t379 + pkin(10) * t362;
t366 = pkin(2) * t402 * t404 + pkin(10) * t377;
t424 = V_base(5) * t366 + (-t343 - t383) * t399 + t428;
t344 = pkin(2) * t381 + pkin(10) * t363;
t423 = V_base(4) * t343 + (-t344 - t384) * V_base(5) + t431;
t301 = pkin(3) * t339 + pkin(11) * t324;
t326 = pkin(3) * t361 + pkin(11) * t336;
t422 = -t301 * t369 + t326 * t356 + t424;
t302 = pkin(3) * t341 + pkin(11) * t325;
t421 = t301 * t357 - t302 * t356 + t423;
t420 = t399 * t344 + (-t366 + t432) * V_base(4) + t425;
t298 = t339 * t449 + (t338 * t447 + t362 * t446) * t410;
t275 = pkin(4) * t298 + pkin(12) * t297;
t319 = t361 * t449 + (t360 * t447 + t377 * t446) * t410;
t289 = pkin(4) * t319 + pkin(12) * t318;
t419 = -t275 * t331 + t289 * t312 + t422;
t300 = t341 * t449 + (t340 * t447 + t363 * t446) * t410;
t276 = pkin(4) * t300 + pkin(12) * t299;
t418 = t275 * t313 - t276 * t312 + t421;
t417 = t302 * t369 - t357 * t326 + t420;
t416 = t276 * t331 - t313 * t289 + t417;
t413 = cos(qJ(6));
t409 = sin(qJ(5));
t408 = sin(qJ(6));
t400 = Icges(2,4) * t415;
t395 = rSges(2,1) * t415 - rSges(2,2) * t412;
t394 = rSges(2,1) * t412 + rSges(2,2) * t415;
t393 = Icges(2,1) * t415 - t445;
t392 = Icges(2,1) * t412 + t400;
t391 = -Icges(2,2) * t412 + t400;
t390 = Icges(2,2) * t415 + t445;
t389 = Icges(2,5) * t415 - Icges(2,6) * t412;
t388 = Icges(2,5) * t412 + Icges(2,6) * t415;
t387 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t386 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t385 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t374 = rSges(3,3) * t407 + (rSges(3,1) * t402 + rSges(3,2) * t405) * t404;
t373 = Icges(3,5) * t407 + (Icges(3,1) * t402 + Icges(3,4) * t405) * t404;
t372 = Icges(3,6) * t407 + (Icges(3,4) * t402 + Icges(3,2) * t405) * t404;
t371 = Icges(3,3) * t407 + (Icges(3,5) * t402 + Icges(3,6) * t405) * t404;
t368 = V_base(5) * rSges(2,3) - t394 * t399 + t435;
t367 = t395 * t399 + V_base(2) + (-rSges(2,3) - pkin(9)) * V_base(4);
t365 = t394 * V_base(4) - t395 * V_base(5) + V_base(3);
t352 = rSges(3,1) * t381 + rSges(3,2) * t380 + rSges(3,3) * t442;
t351 = rSges(3,1) * t379 + rSges(3,2) * t378 - rSges(3,3) * t441;
t350 = Icges(3,1) * t381 + Icges(3,4) * t380 + Icges(3,5) * t442;
t349 = Icges(3,1) * t379 + Icges(3,4) * t378 - Icges(3,5) * t441;
t348 = Icges(3,4) * t381 + Icges(3,2) * t380 + Icges(3,6) * t442;
t347 = Icges(3,4) * t379 + Icges(3,2) * t378 - Icges(3,6) * t441;
t346 = Icges(3,5) * t381 + Icges(3,6) * t380 + Icges(3,3) * t442;
t345 = Icges(3,5) * t379 + Icges(3,6) * t378 - Icges(3,3) * t441;
t330 = rSges(4,1) * t361 + rSges(4,2) * t360 + rSges(4,3) * t377;
t329 = Icges(4,1) * t361 + Icges(4,4) * t360 + Icges(4,5) * t377;
t328 = Icges(4,4) * t361 + Icges(4,2) * t360 + Icges(4,6) * t377;
t327 = Icges(4,5) * t361 + Icges(4,6) * t360 + Icges(4,3) * t377;
t316 = t374 * V_base(5) + (-t351 - t383) * t399 + t428;
t315 = t399 * t352 + (-t374 + t432) * V_base(4) + t425;
t311 = t351 * V_base(4) + (-t352 - t384) * V_base(5) + t431;
t310 = rSges(4,1) * t341 + rSges(4,2) * t340 + rSges(4,3) * t363;
t309 = rSges(4,1) * t339 + rSges(4,2) * t338 + rSges(4,3) * t362;
t308 = Icges(4,1) * t341 + Icges(4,4) * t340 + Icges(4,5) * t363;
t307 = Icges(4,1) * t339 + Icges(4,4) * t338 + Icges(4,5) * t362;
t306 = Icges(4,4) * t341 + Icges(4,2) * t340 + Icges(4,6) * t363;
t305 = Icges(4,4) * t339 + Icges(4,2) * t338 + Icges(4,6) * t362;
t304 = Icges(4,5) * t341 + Icges(4,6) * t340 + Icges(4,3) * t363;
t303 = Icges(4,5) * t339 + Icges(4,6) * t338 + Icges(4,3) * t362;
t294 = t319 * t448 + t336 * t409;
t293 = t319 * t409 - t336 * t448;
t288 = rSges(5,1) * t319 - rSges(5,2) * t318 + rSges(5,3) * t336;
t287 = t300 * t448 + t325 * t409;
t286 = t300 * t409 - t325 * t448;
t285 = t298 * t448 + t324 * t409;
t284 = t298 * t409 - t324 * t448;
t283 = Icges(5,1) * t319 - Icges(5,4) * t318 + Icges(5,5) * t336;
t282 = Icges(5,4) * t319 - Icges(5,2) * t318 + Icges(5,6) * t336;
t281 = Icges(5,5) * t319 - Icges(5,6) * t318 + Icges(5,3) * t336;
t280 = t294 * t413 + t318 * t408;
t279 = -t294 * t408 + t318 * t413;
t274 = pkin(5) * t294 + pkin(13) * t293;
t272 = qJD(6) * t293 + t291;
t270 = rSges(5,1) * t300 - rSges(5,2) * t299 + rSges(5,3) * t325;
t269 = rSges(5,1) * t298 - rSges(5,2) * t297 + rSges(5,3) * t324;
t268 = Icges(5,1) * t300 - Icges(5,4) * t299 + Icges(5,5) * t325;
t267 = Icges(5,1) * t298 - Icges(5,4) * t297 + Icges(5,5) * t324;
t266 = Icges(5,4) * t300 - Icges(5,2) * t299 + Icges(5,6) * t325;
t265 = Icges(5,4) * t298 - Icges(5,2) * t297 + Icges(5,6) * t324;
t264 = Icges(5,5) * t300 - Icges(5,6) * t299 + Icges(5,3) * t325;
t263 = Icges(5,5) * t298 - Icges(5,6) * t297 + Icges(5,3) * t324;
t262 = -t309 * t369 + t330 * t356 + t424;
t261 = t310 * t369 - t330 * t357 + t420;
t260 = t287 * t413 + t299 * t408;
t259 = -t287 * t408 + t299 * t413;
t258 = t285 * t413 + t297 * t408;
t257 = -t285 * t408 + t297 * t413;
t256 = rSges(6,1) * t294 - rSges(6,2) * t293 + rSges(6,3) * t318;
t255 = Icges(6,1) * t294 - Icges(6,4) * t293 + Icges(6,5) * t318;
t254 = Icges(6,4) * t294 - Icges(6,2) * t293 + Icges(6,6) * t318;
t253 = Icges(6,5) * t294 - Icges(6,6) * t293 + Icges(6,3) * t318;
t251 = t309 * t357 - t310 * t356 + t423;
t250 = pkin(5) * t287 + pkin(13) * t286;
t249 = pkin(5) * t285 + pkin(13) * t284;
t248 = qJD(6) * t286 + t278;
t247 = qJD(6) * t284 + t277;
t246 = rSges(6,1) * t287 - rSges(6,2) * t286 + rSges(6,3) * t299;
t245 = rSges(6,1) * t285 - rSges(6,2) * t284 + rSges(6,3) * t297;
t244 = Icges(6,1) * t287 - Icges(6,4) * t286 + Icges(6,5) * t299;
t243 = Icges(6,1) * t285 - Icges(6,4) * t284 + Icges(6,5) * t297;
t242 = Icges(6,4) * t287 - Icges(6,2) * t286 + Icges(6,6) * t299;
t241 = Icges(6,4) * t285 - Icges(6,2) * t284 + Icges(6,6) * t297;
t240 = Icges(6,5) * t287 - Icges(6,6) * t286 + Icges(6,3) * t299;
t239 = Icges(6,5) * t285 - Icges(6,6) * t284 + Icges(6,3) * t297;
t238 = rSges(7,1) * t280 + rSges(7,2) * t279 + rSges(7,3) * t293;
t237 = Icges(7,1) * t280 + Icges(7,4) * t279 + Icges(7,5) * t293;
t236 = Icges(7,4) * t280 + Icges(7,2) * t279 + Icges(7,6) * t293;
t235 = Icges(7,5) * t280 + Icges(7,6) * t279 + Icges(7,3) * t293;
t234 = rSges(7,1) * t260 + rSges(7,2) * t259 + rSges(7,3) * t286;
t233 = rSges(7,1) * t258 + rSges(7,2) * t257 + rSges(7,3) * t284;
t232 = Icges(7,1) * t260 + Icges(7,4) * t259 + Icges(7,5) * t286;
t231 = Icges(7,1) * t258 + Icges(7,4) * t257 + Icges(7,5) * t284;
t230 = Icges(7,4) * t260 + Icges(7,2) * t259 + Icges(7,6) * t286;
t229 = Icges(7,4) * t258 + Icges(7,2) * t257 + Icges(7,6) * t284;
t228 = Icges(7,5) * t260 + Icges(7,6) * t259 + Icges(7,3) * t286;
t227 = Icges(7,5) * t258 + Icges(7,6) * t257 + Icges(7,3) * t284;
t226 = -t269 * t331 + t288 * t312 + t422;
t225 = t331 * t270 - t313 * t288 + t417;
t224 = t269 * t313 - t270 * t312 + t421;
t223 = -t245 * t291 + t256 * t277 + t419;
t222 = t291 * t246 - t278 * t256 + t416;
t221 = t245 * t278 - t246 * t277 + t418;
t220 = -t233 * t272 + t238 * t247 - t249 * t291 + t274 * t277 + t419;
t219 = t234 * t272 - t238 * t248 + t250 * t291 - t274 * t278 + t416;
t218 = t233 * t248 - t234 * t247 + t249 * t278 - t250 * t277 + t418;
t1 = m(5) * (t224 ^ 2 + t225 ^ 2 + t226 ^ 2) / 0.2e1 + m(6) * (t221 ^ 2 + t222 ^ 2 + t223 ^ 2) / 0.2e1 + t331 * ((t264 * t336 - t266 * t318 + t268 * t319) * t313 + (t263 * t336 - t265 * t318 + t267 * t319) * t312 + (t281 * t336 - t282 * t318 + t283 * t319) * t331) / 0.2e1 + t312 * ((t264 * t324 - t266 * t297 + t268 * t298) * t313 + (t263 * t324 - t265 * t297 + t267 * t298) * t312 + (t281 * t324 - t282 * t297 + t283 * t298) * t331) / 0.2e1 + t313 * ((t264 * t325 - t266 * t299 + t268 * t300) * t313 + (t263 * t325 - t265 * t299 + t267 * t300) * t312 + (t281 * t325 - t282 * t299 + t283 * t300) * t331) / 0.2e1 + t248 * ((t286 * t228 + t259 * t230 + t260 * t232) * t248 + (t227 * t286 + t229 * t259 + t231 * t260) * t247 + (t235 * t286 + t236 * t259 + t237 * t260) * t272) / 0.2e1 + t247 * ((t228 * t284 + t230 * t257 + t232 * t258) * t248 + (t284 * t227 + t257 * t229 + t258 * t231) * t247 + (t235 * t284 + t236 * t257 + t237 * t258) * t272) / 0.2e1 + t291 * ((t240 * t318 - t242 * t293 + t244 * t294) * t278 + (t239 * t318 - t241 * t293 + t243 * t294) * t277 + (t253 * t318 - t254 * t293 + t255 * t294) * t291) / 0.2e1 + m(3) * (t311 ^ 2 + t315 ^ 2 + t316 ^ 2) / 0.2e1 + t278 * ((t299 * t240 - t286 * t242 + t287 * t244) * t278 + (t239 * t299 - t241 * t286 + t243 * t287) * t277 + (t253 * t299 - t254 * t286 + t255 * t287) * t291) / 0.2e1 + m(7) * (t218 ^ 2 + t219 ^ 2 + t220 ^ 2) / 0.2e1 + m(2) * (t365 ^ 2 + t367 ^ 2 + t368 ^ 2) / 0.2e1 + t356 * ((t304 * t362 + t306 * t338 + t308 * t339) * t357 + (t303 * t362 + t305 * t338 + t307 * t339) * t356 + (t327 * t362 + t328 * t338 + t329 * t339) * t369) / 0.2e1 + t357 * ((t304 * t363 + t306 * t340 + t308 * t341) * t357 + (t303 * t363 + t305 * t340 + t307 * t341) * t356 + (t327 * t363 + t328 * t340 + t329 * t341) * t369) / 0.2e1 + t369 * ((t304 * t377 + t306 * t360 + t308 * t361) * t357 + (t303 * t377 + t305 * t360 + t307 * t361) * t356 + (t327 * t377 + t328 * t360 + t329 * t361) * t369) / 0.2e1 + m(1) * (t385 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + m(4) * (t251 ^ 2 + t261 ^ 2 + t262 ^ 2) / 0.2e1 + t277 * ((t240 * t297 - t242 * t284 + t244 * t285) * t278 + (t297 * t239 - t284 * t241 + t285 * t243) * t277 + (t253 * t297 - t254 * t284 + t255 * t285) * t291) / 0.2e1 + t272 * ((t228 * t293 + t230 * t279 + t232 * t280) * t248 + (t227 * t293 + t229 * t279 + t231 * t280) * t247 + (t293 * t235 + t279 * t236 + t280 * t237) * t272) / 0.2e1 + ((t345 * V_base(5) + t346 * V_base(4) + t371 * t399) * t407 + ((t348 * t405 + t350 * t402) * V_base(4) + (t347 * t405 + t349 * t402) * V_base(5) + (t372 * t405 + t373 * t402) * t399) * t404 + Icges(2,3) * t399 + t388 * V_base(5) + t389 * V_base(4)) * t399 / 0.2e1 + ((t371 * t442 + t372 * t380 + t373 * t381 + t389) * t399 + (t345 * t442 + t347 * t380 + t349 * t381 - t390 * t412 + t392 * t415 + Icges(1,4)) * V_base(5) + (t346 * t442 + t348 * t380 + t350 * t381 - t412 * t391 + t393 * t415 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t371 * t441 + t372 * t378 + t373 * t379 + t388) * t399 + (-t345 * t441 + t378 * t347 + t379 * t349 + t390 * t415 + t412 * t392 + Icges(1,2)) * V_base(5) + (-t346 * t441 + t348 * t378 + t350 * t379 + t391 * t415 + t393 * t412 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
