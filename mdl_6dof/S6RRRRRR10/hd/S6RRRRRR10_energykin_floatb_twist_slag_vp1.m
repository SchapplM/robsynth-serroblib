% Calculate kinetic energy for
% S6RRRRRR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 06:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_energykin_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR10_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR10_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR10_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 05:49:36
% EndTime: 2019-03-10 05:49:41
% DurationCPUTime: 4.93s
% Computational Cost: add. (9793->450), mult. (26800->682), div. (0->0), fcn. (35064->18), ass. (0->199)
t453 = cos(qJ(4));
t452 = cos(qJ(5));
t409 = cos(pkin(6));
t451 = pkin(10) * t409;
t450 = cos(pkin(8));
t449 = sin(pkin(8));
t415 = sin(qJ(1));
t448 = Icges(2,4) * t415;
t406 = sin(pkin(7));
t447 = t406 * t409;
t407 = sin(pkin(6));
t446 = t407 * t415;
t419 = cos(qJ(1));
t445 = t407 * t419;
t408 = cos(pkin(7));
t418 = cos(qJ(2));
t444 = t408 * t418;
t414 = sin(qJ(2));
t443 = t414 * t419;
t442 = t415 * t414;
t441 = t415 * t418;
t440 = t418 * t419;
t439 = qJD(2) * t407;
t438 = V_base(5) * pkin(9) + V_base(1);
t390 = t415 * t439 + V_base(4);
t403 = V_base(6) + qJD(1);
t381 = -t409 * t441 - t443;
t365 = -t381 * t406 + t408 * t446;
t355 = qJD(3) * t365 + t390;
t391 = qJD(2) * t409 + t403;
t435 = t450 * t453;
t434 = t453 * t449;
t382 = -t409 * t442 + t440;
t413 = sin(qJ(3));
t417 = cos(qJ(3));
t432 = t381 * t408 + t406 * t446;
t342 = -t382 * t413 + t417 * t432;
t326 = -t342 * t449 + t365 * t450;
t314 = qJD(4) * t326 + t355;
t378 = -t406 * t407 * t418 + t408 * t409;
t366 = qJD(3) * t378 + t391;
t389 = -t419 * t439 + V_base(5);
t343 = t382 * t417 + t413 * t432;
t412 = sin(qJ(4));
t300 = -t342 * t435 + t343 * t412 - t365 * t434;
t279 = qJD(5) * t300 + t314;
t379 = t409 * t440 - t442;
t364 = -t379 * t406 - t408 * t445;
t433 = t379 * t408 - t406 * t445;
t384 = t415 * pkin(1) - pkin(10) * t445;
t431 = -t384 * t403 + V_base(5) * t451 + t438;
t362 = t417 * t447 + (-t413 * t414 + t417 * t444) * t407;
t339 = -t362 * t449 + t378 * t450;
t328 = qJD(4) * t339 + t366;
t385 = pkin(1) * t419 + pkin(10) * t446;
t430 = V_base(4) * t384 - t385 * V_base(5) + V_base(3);
t354 = qJD(3) * t364 + t389;
t363 = t413 * t447 + (t413 * t444 + t414 * t417) * t407;
t319 = -t362 * t435 + t363 * t412 - t378 * t434;
t292 = qJD(5) * t319 + t328;
t380 = t409 * t443 + t441;
t340 = -t380 * t413 + t417 * t433;
t325 = -t340 * t449 + t364 * t450;
t313 = qJD(4) * t325 + t354;
t429 = t403 * t385 + V_base(2) + (-pkin(9) - t451) * V_base(4);
t341 = t380 * t417 + t413 * t433;
t298 = -t340 * t435 + t341 * t412 - t364 * t434;
t278 = qJD(5) * t298 + t313;
t344 = t380 * pkin(2) + pkin(11) * t364;
t370 = pkin(2) * t407 * t414 + pkin(11) * t378;
t428 = -t344 * t391 + t389 * t370 + t431;
t345 = pkin(2) * t382 + pkin(11) * t365;
t427 = t390 * t344 - t345 * t389 + t430;
t426 = t391 * t345 - t370 * t390 + t429;
t302 = t341 * pkin(3) + pkin(12) * t325;
t327 = t363 * pkin(3) + pkin(12) * t339;
t425 = -t302 * t366 + t354 * t327 + t428;
t303 = t343 * pkin(3) + pkin(12) * t326;
t424 = t355 * t302 - t303 * t354 + t427;
t423 = t366 * t303 - t327 * t355 + t426;
t299 = t341 * t453 + (t340 * t450 + t364 * t449) * t412;
t276 = pkin(4) * t299 + pkin(13) * t298;
t320 = t363 * t453 + (t362 * t450 + t378 * t449) * t412;
t290 = pkin(4) * t320 + pkin(13) * t319;
t422 = -t276 * t328 + t313 * t290 + t425;
t301 = t343 * t453 + (t342 * t450 + t365 * t449) * t412;
t277 = pkin(4) * t301 + pkin(13) * t300;
t421 = t314 * t276 - t277 * t313 + t424;
t420 = t328 * t277 - t290 * t314 + t423;
t416 = cos(qJ(6));
t411 = sin(qJ(5));
t410 = sin(qJ(6));
t404 = Icges(2,4) * t419;
t399 = rSges(2,1) * t419 - t415 * rSges(2,2);
t398 = t415 * rSges(2,1) + rSges(2,2) * t419;
t397 = Icges(2,1) * t419 - t448;
t396 = Icges(2,1) * t415 + t404;
t395 = -Icges(2,2) * t415 + t404;
t394 = Icges(2,2) * t419 + t448;
t388 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t387 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t386 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t375 = rSges(3,3) * t409 + (rSges(3,1) * t414 + rSges(3,2) * t418) * t407;
t374 = Icges(3,5) * t409 + (Icges(3,1) * t414 + Icges(3,4) * t418) * t407;
t373 = Icges(3,6) * t409 + (Icges(3,4) * t414 + Icges(3,2) * t418) * t407;
t372 = Icges(3,3) * t409 + (Icges(3,5) * t414 + Icges(3,6) * t418) * t407;
t369 = V_base(5) * rSges(2,3) - t398 * t403 + t438;
t368 = t399 * t403 + V_base(2) + (-rSges(2,3) - pkin(9)) * V_base(4);
t367 = t398 * V_base(4) - t399 * V_base(5) + V_base(3);
t353 = rSges(3,1) * t382 + rSges(3,2) * t381 + rSges(3,3) * t446;
t352 = t380 * rSges(3,1) + t379 * rSges(3,2) - rSges(3,3) * t445;
t351 = Icges(3,1) * t382 + Icges(3,4) * t381 + Icges(3,5) * t446;
t350 = Icges(3,1) * t380 + Icges(3,4) * t379 - Icges(3,5) * t445;
t349 = Icges(3,4) * t382 + Icges(3,2) * t381 + Icges(3,6) * t446;
t348 = Icges(3,4) * t380 + Icges(3,2) * t379 - Icges(3,6) * t445;
t347 = Icges(3,5) * t382 + Icges(3,6) * t381 + Icges(3,3) * t446;
t346 = Icges(3,5) * t380 + Icges(3,6) * t379 - Icges(3,3) * t445;
t332 = rSges(4,1) * t363 + rSges(4,2) * t362 + rSges(4,3) * t378;
t331 = Icges(4,1) * t363 + Icges(4,4) * t362 + Icges(4,5) * t378;
t330 = Icges(4,4) * t363 + Icges(4,2) * t362 + Icges(4,6) * t378;
t329 = Icges(4,5) * t363 + Icges(4,6) * t362 + Icges(4,3) * t378;
t317 = -t352 * t391 + t375 * t389 + t431;
t316 = t353 * t391 - t375 * t390 + t429;
t312 = rSges(4,1) * t343 + rSges(4,2) * t342 + rSges(4,3) * t365;
t311 = rSges(4,1) * t341 + rSges(4,2) * t340 + rSges(4,3) * t364;
t310 = t352 * t390 - t353 * t389 + t430;
t309 = Icges(4,1) * t343 + Icges(4,4) * t342 + Icges(4,5) * t365;
t308 = Icges(4,1) * t341 + Icges(4,4) * t340 + Icges(4,5) * t364;
t307 = Icges(4,4) * t343 + Icges(4,2) * t342 + Icges(4,6) * t365;
t306 = Icges(4,4) * t341 + Icges(4,2) * t340 + Icges(4,6) * t364;
t305 = Icges(4,5) * t343 + Icges(4,6) * t342 + Icges(4,3) * t365;
t304 = Icges(4,5) * t341 + Icges(4,6) * t340 + Icges(4,3) * t364;
t295 = t320 * t452 + t339 * t411;
t294 = t320 * t411 - t339 * t452;
t289 = t301 * t452 + t326 * t411;
t288 = t301 * t411 - t326 * t452;
t287 = t299 * t452 + t325 * t411;
t286 = t299 * t411 - t325 * t452;
t285 = rSges(5,1) * t320 - rSges(5,2) * t319 + rSges(5,3) * t339;
t284 = Icges(5,1) * t320 - Icges(5,4) * t319 + Icges(5,5) * t339;
t283 = Icges(5,4) * t320 - Icges(5,2) * t319 + Icges(5,6) * t339;
t282 = Icges(5,5) * t320 - Icges(5,6) * t319 + Icges(5,3) * t339;
t281 = t295 * t416 + t319 * t410;
t280 = -t295 * t410 + t319 * t416;
t275 = pkin(5) * t295 + pkin(14) * t294;
t273 = qJD(6) * t294 + t292;
t271 = rSges(5,1) * t301 - rSges(5,2) * t300 + rSges(5,3) * t326;
t270 = rSges(5,1) * t299 - rSges(5,2) * t298 + rSges(5,3) * t325;
t269 = Icges(5,1) * t301 - Icges(5,4) * t300 + Icges(5,5) * t326;
t268 = Icges(5,1) * t299 - Icges(5,4) * t298 + Icges(5,5) * t325;
t267 = Icges(5,4) * t301 - Icges(5,2) * t300 + Icges(5,6) * t326;
t266 = Icges(5,4) * t299 - Icges(5,2) * t298 + Icges(5,6) * t325;
t265 = Icges(5,5) * t301 - Icges(5,6) * t300 + Icges(5,3) * t326;
t264 = Icges(5,5) * t299 - Icges(5,6) * t298 + Icges(5,3) * t325;
t263 = t289 * t416 + t300 * t410;
t262 = -t289 * t410 + t300 * t416;
t261 = t287 * t416 + t298 * t410;
t260 = -t287 * t410 + t298 * t416;
t259 = rSges(6,1) * t295 - rSges(6,2) * t294 + rSges(6,3) * t319;
t258 = Icges(6,1) * t295 - Icges(6,4) * t294 + Icges(6,5) * t319;
t257 = Icges(6,4) * t295 - Icges(6,2) * t294 + Icges(6,6) * t319;
t256 = Icges(6,5) * t295 - Icges(6,6) * t294 + Icges(6,3) * t319;
t255 = -t311 * t366 + t332 * t354 + t428;
t254 = t312 * t366 - t332 * t355 + t426;
t252 = pkin(5) * t289 + pkin(14) * t288;
t251 = pkin(5) * t287 + pkin(14) * t286;
t250 = t311 * t355 - t312 * t354 + t427;
t249 = qJD(6) * t288 + t279;
t248 = qJD(6) * t286 + t278;
t247 = rSges(6,1) * t289 - rSges(6,2) * t288 + rSges(6,3) * t300;
t246 = rSges(6,1) * t287 - rSges(6,2) * t286 + rSges(6,3) * t298;
t245 = Icges(6,1) * t289 - Icges(6,4) * t288 + Icges(6,5) * t300;
t244 = Icges(6,1) * t287 - Icges(6,4) * t286 + Icges(6,5) * t298;
t243 = Icges(6,4) * t289 - Icges(6,2) * t288 + Icges(6,6) * t300;
t242 = Icges(6,4) * t287 - Icges(6,2) * t286 + Icges(6,6) * t298;
t241 = Icges(6,5) * t289 - Icges(6,6) * t288 + Icges(6,3) * t300;
t240 = Icges(6,5) * t287 - Icges(6,6) * t286 + Icges(6,3) * t298;
t239 = rSges(7,1) * t281 + rSges(7,2) * t280 + rSges(7,3) * t294;
t238 = Icges(7,1) * t281 + Icges(7,4) * t280 + Icges(7,5) * t294;
t237 = Icges(7,4) * t281 + Icges(7,2) * t280 + Icges(7,6) * t294;
t236 = Icges(7,5) * t281 + Icges(7,6) * t280 + Icges(7,3) * t294;
t235 = rSges(7,1) * t263 + rSges(7,2) * t262 + rSges(7,3) * t288;
t234 = rSges(7,1) * t261 + rSges(7,2) * t260 + rSges(7,3) * t286;
t233 = Icges(7,1) * t263 + Icges(7,4) * t262 + Icges(7,5) * t288;
t232 = Icges(7,1) * t261 + Icges(7,4) * t260 + Icges(7,5) * t286;
t231 = Icges(7,4) * t263 + Icges(7,2) * t262 + Icges(7,6) * t288;
t230 = Icges(7,4) * t261 + Icges(7,2) * t260 + Icges(7,6) * t286;
t229 = Icges(7,5) * t263 + Icges(7,6) * t262 + Icges(7,3) * t288;
t228 = Icges(7,5) * t261 + Icges(7,6) * t260 + Icges(7,3) * t286;
t227 = -t270 * t328 + t285 * t313 + t425;
t226 = t271 * t328 - t285 * t314 + t423;
t225 = t270 * t314 - t271 * t313 + t424;
t224 = -t246 * t292 + t259 * t278 + t422;
t223 = t247 * t292 - t259 * t279 + t420;
t222 = t246 * t279 - t247 * t278 + t421;
t221 = -t234 * t273 + t239 * t248 - t251 * t292 + t275 * t278 + t422;
t220 = t235 * t273 - t239 * t249 + t252 * t292 - t275 * t279 + t420;
t219 = t234 * t249 - t235 * t248 + t251 * t279 - t252 * t278 + t421;
t1 = m(4) * (t250 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + t389 * ((-t347 * t445 + t379 * t349 + t380 * t351) * t390 + (-t346 * t445 + t379 * t348 + t380 * t350) * t389 + (-t372 * t445 + t379 * t373 + t380 * t374) * t391) / 0.2e1 + t390 * ((t347 * t446 + t349 * t381 + t351 * t382) * t390 + (t346 * t446 + t348 * t381 + t350 * t382) * t389 + (t372 * t446 + t373 * t381 + t374 * t382) * t391) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + t328 * ((t265 * t339 - t267 * t319 + t269 * t320) * t314 + (t264 * t339 - t266 * t319 + t268 * t320) * t313 + (t282 * t339 - t283 * t319 + t284 * t320) * t328) / 0.2e1 + t313 * ((t265 * t325 - t267 * t298 + t269 * t299) * t314 + (t264 * t325 - t266 * t298 + t268 * t299) * t313 + (t282 * t325 - t283 * t298 + t284 * t299) * t328) / 0.2e1 + t314 * ((t265 * t326 - t267 * t300 + t269 * t301) * t314 + (t264 * t326 - t266 * t300 + t268 * t301) * t313 + (t282 * t326 - t283 * t300 + t284 * t301) * t328) / 0.2e1 + m(3) * (t310 ^ 2 + t316 ^ 2 + t317 ^ 2) / 0.2e1 + t292 * ((t241 * t319 - t243 * t294 + t245 * t295) * t279 + (t240 * t319 - t242 * t294 + t244 * t295) * t278 + (t256 * t319 - t257 * t294 + t258 * t295) * t292) / 0.2e1 + ((Icges(2,5) * t415 + Icges(2,6) * t419) * V_base(5) + (Icges(2,5) * t419 - Icges(2,6) * t415) * V_base(4) + Icges(2,3) * t403 / 0.2e1) * t403 + m(1) * (t386 ^ 2 + t387 ^ 2 + t388 ^ 2) / 0.2e1 + m(6) * (t222 ^ 2 + t223 ^ 2 + t224 ^ 2) / 0.2e1 + t278 * ((t241 * t298 - t243 * t286 + t245 * t287) * t279 + (t298 * t240 - t286 * t242 + t287 * t244) * t278 + (t256 * t298 - t257 * t286 + t258 * t287) * t292) / 0.2e1 + t279 * ((t300 * t241 - t288 * t243 + t289 * t245) * t279 + (t240 * t300 - t242 * t288 + t244 * t289) * t278 + (t256 * t300 - t257 * t288 + t258 * t289) * t292) / 0.2e1 + m(7) * (t219 ^ 2 + t220 ^ 2 + t221 ^ 2) / 0.2e1 + t366 * ((t305 * t378 + t307 * t362 + t309 * t363) * t355 + (t304 * t378 + t306 * t362 + t308 * t363) * t354 + (t329 * t378 + t330 * t362 + t331 * t363) * t366) / 0.2e1 + t354 * ((t305 * t364 + t307 * t340 + t309 * t341) * t355 + (t304 * t364 + t306 * t340 + t308 * t341) * t354 + (t329 * t364 + t330 * t340 + t331 * t341) * t366) / 0.2e1 + t355 * ((t305 * t365 + t307 * t342 + t309 * t343) * t355 + (t304 * t365 + t306 * t342 + t308 * t343) * t354 + (t329 * t365 + t330 * t342 + t331 * t343) * t366) / 0.2e1 + m(2) * (t367 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + t391 * ((t346 * t389 + t347 * t390 + t372 * t391) * t409 + ((t349 * t418 + t351 * t414) * t390 + (t348 * t418 + t350 * t414) * t389 + (t373 * t418 + t374 * t414) * t391) * t407) / 0.2e1 + m(5) * (t225 ^ 2 + t226 ^ 2 + t227 ^ 2) / 0.2e1 + t273 * ((t229 * t294 + t231 * t280 + t233 * t281) * t249 + (t228 * t294 + t230 * t280 + t232 * t281) * t248 + (t294 * t236 + t280 * t237 + t281 * t238) * t273) / 0.2e1 + t249 * ((t288 * t229 + t262 * t231 + t263 * t233) * t249 + (t228 * t288 + t230 * t262 + t232 * t263) * t248 + (t236 * t288 + t237 * t262 + t238 * t263) * t273) / 0.2e1 + t248 * ((t229 * t286 + t231 * t260 + t233 * t261) * t249 + (t286 * t228 + t260 * t230 + t261 * t232) * t248 + (t236 * t286 + t237 * t260 + t238 * t261) * t273) / 0.2e1 + ((-t415 * t394 + t396 * t419 + Icges(1,4)) * V_base(5) + (-t415 * t395 + t397 * t419 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t394 * t419 + t415 * t396 + Icges(1,2)) * V_base(5) + (t395 * t419 + t415 * t397 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1;
T  = t1;
