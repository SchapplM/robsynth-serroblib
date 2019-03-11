% Calculate kinetic energy for
% S6RRPRRR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 15:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR14_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR14_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:57:11
% EndTime: 2019-03-09 14:57:16
% DurationCPUTime: 5.17s
% Computational Cost: add. (9537->453), mult. (26224->666), div. (0->0), fcn. (34328->18), ass. (0->197)
t453 = cos(qJ(4));
t452 = cos(qJ(5));
t410 = cos(pkin(6));
t451 = pkin(10) * t410;
t450 = cos(pkin(8));
t449 = sin(pkin(8));
t415 = sin(qJ(1));
t448 = Icges(2,4) * t415;
t406 = sin(pkin(7));
t447 = t406 * t410;
t407 = sin(pkin(6));
t414 = sin(qJ(2));
t446 = t407 * t414;
t445 = t407 * t415;
t418 = cos(qJ(1));
t444 = t407 * t418;
t409 = cos(pkin(7));
t417 = cos(qJ(2));
t443 = t409 * t417;
t442 = t414 * t418;
t441 = t415 * t414;
t440 = t415 * t417;
t439 = t417 * t418;
t438 = qJD(2) * t407;
t437 = V_base(5) * pkin(9) + V_base(1);
t389 = t415 * t438 + V_base(4);
t402 = V_base(6) + qJD(1);
t381 = -t410 * t441 + t439;
t405 = sin(pkin(14));
t408 = cos(pkin(14));
t380 = -t410 * t440 - t442;
t431 = t380 * t409 + t406 * t445;
t344 = -t381 * t405 + t408 * t431;
t365 = -t380 * t406 + t409 * t445;
t329 = -t344 * t449 + t365 * t450;
t318 = qJD(4) * t329 + t389;
t390 = qJD(2) * t410 + t402;
t434 = t450 * t453;
t433 = t453 * t449;
t345 = t381 * t408 + t405 * t431;
t413 = sin(qJ(4));
t302 = -t344 * t434 + t345 * t413 - t365 * t433;
t283 = qJD(5) * t302 + t318;
t360 = t408 * t447 + (-t405 * t414 + t408 * t443) * t407;
t377 = -t406 * t407 * t417 + t409 * t410;
t341 = -t360 * t449 + t377 * t450;
t334 = qJD(4) * t341 + t390;
t388 = -t418 * t438 + V_base(5);
t378 = t410 * t439 - t441;
t364 = -t378 * t406 - t409 * t444;
t432 = t378 * t409 - t406 * t444;
t383 = pkin(1) * t415 - pkin(10) * t444;
t430 = -t383 * t402 + t451 * V_base(5) + t437;
t361 = t408 * t446 + (t407 * t443 + t447) * t405;
t321 = -t360 * t434 + t361 * t413 - t377 * t433;
t293 = qJD(5) * t321 + t334;
t384 = pkin(1) * t418 + pkin(10) * t445;
t429 = t383 * V_base(4) - t384 * V_base(5) + V_base(3);
t379 = t410 * t442 + t440;
t342 = -t379 * t405 + t408 * t432;
t328 = -t342 * t449 + t364 * t450;
t317 = qJD(4) * t328 + t388;
t343 = t379 * t408 + t405 * t432;
t300 = -t342 * t434 + t343 * t413 - t364 * t433;
t282 = qJD(5) * t300 + t317;
t367 = pkin(2) * t446 + qJ(3) * t377;
t428 = qJD(3) * t365 + t367 * t388 + t430;
t346 = pkin(2) * t379 + qJ(3) * t364;
t427 = qJD(3) * t377 + t346 * t389 + t429;
t426 = t402 * t384 + V_base(2) + (-pkin(9) - t451) * V_base(4);
t347 = pkin(2) * t381 + qJ(3) * t365;
t425 = qJD(3) * t364 + t347 * t390 + t426;
t304 = pkin(3) * t343 + pkin(11) * t328;
t327 = pkin(3) * t361 + pkin(11) * t341;
t424 = t388 * t327 + (-t304 - t346) * t390 + t428;
t305 = pkin(3) * t345 + pkin(11) * t329;
t423 = t389 * t304 + (-t305 - t347) * t388 + t427;
t422 = t390 * t305 + (-t327 - t367) * t389 + t425;
t301 = t343 * t453 + (t342 * t450 + t364 * t449) * t413;
t277 = pkin(4) * t301 + pkin(12) * t300;
t322 = t361 * t453 + (t360 * t450 + t377 * t449) * t413;
t292 = pkin(4) * t322 + pkin(12) * t321;
t421 = -t277 * t334 + t292 * t317 + t424;
t303 = t345 * t453 + (t344 * t450 + t365 * t449) * t413;
t278 = pkin(4) * t303 + pkin(12) * t302;
t420 = t277 * t318 - t278 * t317 + t423;
t419 = t278 * t334 - t292 * t318 + t422;
t416 = cos(qJ(6));
t412 = sin(qJ(5));
t411 = sin(qJ(6));
t403 = Icges(2,4) * t418;
t398 = rSges(2,1) * t418 - rSges(2,2) * t415;
t397 = rSges(2,1) * t415 + rSges(2,2) * t418;
t396 = Icges(2,1) * t418 - t448;
t395 = Icges(2,1) * t415 + t403;
t394 = -Icges(2,2) * t415 + t403;
t393 = Icges(2,2) * t418 + t448;
t387 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t386 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t385 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t374 = rSges(3,3) * t410 + (rSges(3,1) * t414 + rSges(3,2) * t417) * t407;
t373 = Icges(3,5) * t410 + (Icges(3,1) * t414 + Icges(3,4) * t417) * t407;
t372 = Icges(3,6) * t410 + (Icges(3,4) * t414 + Icges(3,2) * t417) * t407;
t371 = Icges(3,3) * t410 + (Icges(3,5) * t414 + Icges(3,6) * t417) * t407;
t369 = V_base(5) * rSges(2,3) - t397 * t402 + t437;
t368 = t398 * t402 + V_base(2) + (-rSges(2,3) - pkin(9)) * V_base(4);
t366 = t397 * V_base(4) - t398 * V_base(5) + V_base(3);
t355 = rSges(3,1) * t381 + rSges(3,2) * t380 + rSges(3,3) * t445;
t354 = rSges(3,1) * t379 + rSges(3,2) * t378 - rSges(3,3) * t444;
t353 = Icges(3,1) * t381 + Icges(3,4) * t380 + Icges(3,5) * t445;
t352 = Icges(3,1) * t379 + Icges(3,4) * t378 - Icges(3,5) * t444;
t351 = Icges(3,4) * t381 + Icges(3,2) * t380 + Icges(3,6) * t445;
t350 = Icges(3,4) * t379 + Icges(3,2) * t378 - Icges(3,6) * t444;
t349 = Icges(3,5) * t381 + Icges(3,6) * t380 + Icges(3,3) * t445;
t348 = Icges(3,5) * t379 + Icges(3,6) * t378 - Icges(3,3) * t444;
t333 = rSges(4,1) * t361 + rSges(4,2) * t360 + rSges(4,3) * t377;
t332 = Icges(4,1) * t361 + Icges(4,4) * t360 + Icges(4,5) * t377;
t331 = Icges(4,4) * t361 + Icges(4,2) * t360 + Icges(4,6) * t377;
t330 = Icges(4,5) * t361 + Icges(4,6) * t360 + Icges(4,3) * t377;
t316 = -t354 * t390 + t374 * t388 + t430;
t315 = t355 * t390 - t374 * t389 + t426;
t314 = t354 * t389 - t355 * t388 + t429;
t313 = rSges(4,1) * t345 + rSges(4,2) * t344 + rSges(4,3) * t365;
t312 = rSges(4,1) * t343 + rSges(4,2) * t342 + rSges(4,3) * t364;
t311 = Icges(4,1) * t345 + Icges(4,4) * t344 + Icges(4,5) * t365;
t310 = Icges(4,1) * t343 + Icges(4,4) * t342 + Icges(4,5) * t364;
t309 = Icges(4,4) * t345 + Icges(4,2) * t344 + Icges(4,6) * t365;
t308 = Icges(4,4) * t343 + Icges(4,2) * t342 + Icges(4,6) * t364;
t307 = Icges(4,5) * t345 + Icges(4,6) * t344 + Icges(4,3) * t365;
t306 = Icges(4,5) * t343 + Icges(4,6) * t342 + Icges(4,3) * t364;
t296 = t322 * t452 + t341 * t412;
t295 = t322 * t412 - t341 * t452;
t291 = t303 * t452 + t329 * t412;
t290 = t303 * t412 - t329 * t452;
t289 = t301 * t452 + t328 * t412;
t288 = t301 * t412 - t328 * t452;
t287 = rSges(5,1) * t322 - rSges(5,2) * t321 + rSges(5,3) * t341;
t286 = Icges(5,1) * t322 - Icges(5,4) * t321 + Icges(5,5) * t341;
t285 = Icges(5,4) * t322 - Icges(5,2) * t321 + Icges(5,6) * t341;
t284 = Icges(5,5) * t322 - Icges(5,6) * t321 + Icges(5,3) * t341;
t281 = t296 * t416 + t321 * t411;
t280 = -t296 * t411 + t321 * t416;
t276 = pkin(5) * t296 + pkin(13) * t295;
t275 = qJD(6) * t295 + t293;
t273 = rSges(5,1) * t303 - rSges(5,2) * t302 + rSges(5,3) * t329;
t272 = rSges(5,1) * t301 - rSges(5,2) * t300 + rSges(5,3) * t328;
t271 = Icges(5,1) * t303 - Icges(5,4) * t302 + Icges(5,5) * t329;
t270 = Icges(5,1) * t301 - Icges(5,4) * t300 + Icges(5,5) * t328;
t269 = Icges(5,4) * t303 - Icges(5,2) * t302 + Icges(5,6) * t329;
t268 = Icges(5,4) * t301 - Icges(5,2) * t300 + Icges(5,6) * t328;
t267 = Icges(5,5) * t303 - Icges(5,6) * t302 + Icges(5,3) * t329;
t266 = Icges(5,5) * t301 - Icges(5,6) * t300 + Icges(5,3) * t328;
t265 = t333 * t388 + (-t312 - t346) * t390 + t428;
t264 = t313 * t390 + (-t333 - t367) * t389 + t425;
t263 = t291 * t416 + t302 * t411;
t262 = -t291 * t411 + t302 * t416;
t261 = t289 * t416 + t300 * t411;
t260 = -t289 * t411 + t300 * t416;
t258 = rSges(6,1) * t296 - rSges(6,2) * t295 + rSges(6,3) * t321;
t257 = Icges(6,1) * t296 - Icges(6,4) * t295 + Icges(6,5) * t321;
t256 = Icges(6,4) * t296 - Icges(6,2) * t295 + Icges(6,6) * t321;
t255 = Icges(6,5) * t296 - Icges(6,6) * t295 + Icges(6,3) * t321;
t254 = t312 * t389 + (-t313 - t347) * t388 + t427;
t253 = pkin(5) * t291 + pkin(13) * t290;
t252 = pkin(5) * t289 + pkin(13) * t288;
t251 = qJD(6) * t290 + t283;
t250 = qJD(6) * t288 + t282;
t249 = rSges(6,1) * t291 - rSges(6,2) * t290 + rSges(6,3) * t302;
t248 = rSges(6,1) * t289 - rSges(6,2) * t288 + rSges(6,3) * t300;
t247 = Icges(6,1) * t291 - Icges(6,4) * t290 + Icges(6,5) * t302;
t246 = Icges(6,1) * t289 - Icges(6,4) * t288 + Icges(6,5) * t300;
t245 = Icges(6,4) * t291 - Icges(6,2) * t290 + Icges(6,6) * t302;
t244 = Icges(6,4) * t289 - Icges(6,2) * t288 + Icges(6,6) * t300;
t243 = Icges(6,5) * t291 - Icges(6,6) * t290 + Icges(6,3) * t302;
t242 = Icges(6,5) * t289 - Icges(6,6) * t288 + Icges(6,3) * t300;
t241 = rSges(7,1) * t281 + rSges(7,2) * t280 + rSges(7,3) * t295;
t240 = Icges(7,1) * t281 + Icges(7,4) * t280 + Icges(7,5) * t295;
t239 = Icges(7,4) * t281 + Icges(7,2) * t280 + Icges(7,6) * t295;
t238 = Icges(7,5) * t281 + Icges(7,6) * t280 + Icges(7,3) * t295;
t237 = rSges(7,1) * t263 + rSges(7,2) * t262 + rSges(7,3) * t290;
t236 = rSges(7,1) * t261 + rSges(7,2) * t260 + rSges(7,3) * t288;
t235 = Icges(7,1) * t263 + Icges(7,4) * t262 + Icges(7,5) * t290;
t234 = Icges(7,1) * t261 + Icges(7,4) * t260 + Icges(7,5) * t288;
t233 = Icges(7,4) * t263 + Icges(7,2) * t262 + Icges(7,6) * t290;
t232 = Icges(7,4) * t261 + Icges(7,2) * t260 + Icges(7,6) * t288;
t231 = Icges(7,5) * t263 + Icges(7,6) * t262 + Icges(7,3) * t290;
t230 = Icges(7,5) * t261 + Icges(7,6) * t260 + Icges(7,3) * t288;
t229 = -t272 * t334 + t287 * t317 + t424;
t228 = t273 * t334 - t287 * t318 + t422;
t227 = t272 * t318 - t273 * t317 + t423;
t226 = -t248 * t293 + t258 * t282 + t421;
t225 = t249 * t293 - t258 * t283 + t419;
t224 = t248 * t283 - t249 * t282 + t420;
t223 = -t236 * t275 + t241 * t250 - t252 * t293 + t276 * t282 + t421;
t222 = t237 * t275 - t241 * t251 + t253 * t293 - t276 * t283 + t419;
t221 = t236 * t251 - t237 * t250 + t252 * t283 - t253 * t282 + t420;
t1 = t282 * ((t243 * t300 - t245 * t288 + t247 * t289) * t283 + (t300 * t242 - t288 * t244 + t289 * t246) * t282 + (t255 * t300 - t256 * t288 + t257 * t289) * t293) / 0.2e1 + t283 * ((t302 * t243 - t290 * t245 + t291 * t247) * t283 + (t242 * t302 - t244 * t290 + t246 * t291) * t282 + (t255 * t302 - t256 * t290 + t257 * t291) * t293) / 0.2e1 + t275 * ((t231 * t295 + t233 * t280 + t235 * t281) * t251 + (t230 * t295 + t232 * t280 + t234 * t281) * t250 + (t295 * t238 + t280 * t239 + t281 * t240) * t275) / 0.2e1 + t251 * ((t290 * t231 + t262 * t233 + t263 * t235) * t251 + (t230 * t290 + t232 * t262 + t234 * t263) * t250 + (t238 * t290 + t239 * t262 + t240 * t263) * t275) / 0.2e1 + t250 * ((t231 * t288 + t233 * t260 + t235 * t261) * t251 + (t288 * t230 + t260 * t232 + t261 * t234) * t250 + (t238 * t288 + t239 * t260 + t240 * t261) * t275) / 0.2e1 + m(4) * (t254 ^ 2 + t264 ^ 2 + t265 ^ 2) / 0.2e1 + m(5) * (t227 ^ 2 + t228 ^ 2 + t229 ^ 2) / 0.2e1 + m(6) * (t224 ^ 2 + t225 ^ 2 + t226 ^ 2) / 0.2e1 + m(3) * (t314 ^ 2 + t315 ^ 2 + t316 ^ 2) / 0.2e1 + t293 * ((t243 * t321 - t245 * t295 + t247 * t296) * t283 + (t242 * t321 - t244 * t295 + t246 * t296) * t282 + (t321 * t255 - t295 * t256 + t296 * t257) * t293) / 0.2e1 + t317 * ((t267 * t328 - t269 * t300 + t271 * t301) * t318 + (t266 * t328 - t268 * t300 + t270 * t301) * t317 + (t284 * t328 - t285 * t300 + t286 * t301) * t334) / 0.2e1 + t318 * ((t267 * t329 - t269 * t302 + t271 * t303) * t318 + (t266 * t329 - t268 * t302 + t270 * t303) * t317 + (t284 * t329 - t285 * t302 + t286 * t303) * t334) / 0.2e1 + t334 * ((t267 * t341 - t269 * t321 + t271 * t322) * t318 + (t266 * t341 - t268 * t321 + t270 * t322) * t317 + (t284 * t341 - t285 * t321 + t286 * t322) * t334) / 0.2e1 + m(2) * (t366 ^ 2 + t368 ^ 2 + t369 ^ 2) / 0.2e1 + m(7) * (t221 ^ 2 + t222 ^ 2 + t223 ^ 2) / 0.2e1 + m(1) * (t385 ^ 2 + t386 ^ 2 + t387 ^ 2) / 0.2e1 + ((t330 * t364 + t331 * t342 + t332 * t343 - t371 * t444 + t372 * t378 + t373 * t379) * t390 + (t307 * t364 + t309 * t342 + t311 * t343 - t349 * t444 + t351 * t378 + t353 * t379) * t389 + (t306 * t364 + t308 * t342 + t310 * t343 - t348 * t444 + t378 * t350 + t379 * t352) * t388) * t388 / 0.2e1 + ((t330 * t365 + t331 * t344 + t332 * t345 + t371 * t445 + t372 * t380 + t373 * t381) * t390 + (t307 * t365 + t309 * t344 + t311 * t345 + t349 * t445 + t351 * t380 + t353 * t381) * t389 + (t306 * t365 + t308 * t344 + t310 * t345 + t348 * t445 + t350 * t380 + t352 * t381) * t388) * t389 / 0.2e1 + ((t348 * t388 + t349 * t389 + t371 * t390) * t410 + ((t351 * t417 + t353 * t414) * t389 + (t350 * t417 + t352 * t414) * t388 + (t372 * t417 + t373 * t414) * t390) * t407 + (t307 * t377 + t309 * t360 + t311 * t361) * t389 + (t306 * t377 + t308 * t360 + t310 * t361) * t388 + (t330 * t377 + t331 * t360 + t332 * t361) * t390) * t390 / 0.2e1 + ((-t393 * t415 + t395 * t418 + Icges(1,4)) * V_base(5) + (-t415 * t394 + t396 * t418 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t393 * t418 + t415 * t395 + Icges(1,2)) * V_base(5) + (t394 * t418 + t396 * t415 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t415 + Icges(2,6) * t418) * V_base(5) + (Icges(2,5) * t418 - Icges(2,6) * t415) * V_base(4) + Icges(2,3) * t402 / 0.2e1) * t402;
T  = t1;
