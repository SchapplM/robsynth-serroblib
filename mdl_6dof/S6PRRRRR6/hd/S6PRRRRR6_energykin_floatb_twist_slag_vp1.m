% Calculate kinetic energy for
% S6PRRRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_energykin_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:13:59
% EndTime: 2019-03-09 01:14:04
% DurationCPUTime: 4.90s
% Computational Cost: add. (9733->450), mult. (26800->681), div. (0->0), fcn. (35064->18), ass. (0->199)
t452 = cos(qJ(4));
t451 = cos(qJ(5));
t408 = cos(pkin(6));
t450 = pkin(9) * t408;
t449 = cos(pkin(8));
t448 = sin(pkin(8));
t403 = sin(pkin(14));
t447 = Icges(2,4) * t403;
t405 = sin(pkin(6));
t446 = t403 * t405;
t404 = sin(pkin(7));
t445 = t404 * t405;
t444 = t404 * t408;
t406 = cos(pkin(14));
t443 = t405 * t406;
t407 = cos(pkin(7));
t442 = t405 * t407;
t416 = cos(qJ(2));
t441 = t407 * t416;
t413 = sin(qJ(2));
t440 = t408 * t413;
t439 = t408 * t416;
t438 = qJD(2) * t405;
t437 = V_base(5) * qJ(1) + V_base(1);
t433 = qJD(1) + V_base(3);
t388 = t403 * t438 + V_base(4);
t397 = qJD(2) * t408 + V_base(6);
t378 = -t403 * t439 - t406 * t413;
t364 = -t378 * t404 + t403 * t442;
t353 = qJD(3) * t364 + t388;
t375 = t407 * t408 - t416 * t445;
t365 = qJD(3) * t375 + t397;
t432 = t449 * t452;
t431 = t452 * t448;
t379 = -t403 * t440 + t406 * t416;
t412 = sin(qJ(3));
t415 = cos(qJ(3));
t429 = t378 * t407 + t403 * t445;
t340 = -t379 * t412 + t415 * t429;
t324 = -t340 * t448 + t364 * t449;
t312 = qJD(4) * t324 + t353;
t361 = t415 * t444 + (-t412 * t413 + t415 * t441) * t405;
t337 = -t361 * t448 + t375 * t449;
t326 = qJD(4) * t337 + t365;
t387 = -t406 * t438 + V_base(5);
t341 = t379 * t415 + t412 * t429;
t411 = sin(qJ(4));
t298 = -t340 * t432 + t341 * t411 - t364 * t431;
t277 = qJD(5) * t298 + t312;
t362 = t412 * t444 + (t412 * t441 + t413 * t415) * t405;
t321 = -t361 * t432 + t362 * t411 - t375 * t431;
t290 = qJD(5) * t321 + t326;
t376 = -t403 * t413 + t406 * t439;
t363 = -t376 * t404 - t406 * t442;
t430 = t376 * t407 - t404 * t443;
t352 = qJD(3) * t363 + t387;
t382 = pkin(1) * t403 - pkin(9) * t443;
t428 = -t382 * V_base(6) + t450 * V_base(5) + t437;
t383 = pkin(1) * t406 + pkin(9) * t446;
t427 = t382 * V_base(4) - t383 * V_base(5) + t433;
t377 = t403 * t416 + t406 * t440;
t338 = -t377 * t412 + t415 * t430;
t323 = -t338 * t448 + t363 * t449;
t311 = qJD(4) * t323 + t352;
t339 = t377 * t415 + t412 * t430;
t296 = -t338 * t432 + t339 * t411 - t363 * t431;
t276 = qJD(5) * t296 + t311;
t426 = V_base(6) * t383 + V_base(2) + (-qJ(1) - t450) * V_base(4);
t342 = pkin(2) * t377 + pkin(10) * t363;
t366 = pkin(2) * t405 * t413 + pkin(10) * t375;
t425 = -t342 * t397 + t366 * t387 + t428;
t343 = pkin(2) * t379 + pkin(10) * t364;
t424 = t342 * t388 - t343 * t387 + t427;
t423 = t343 * t397 - t366 * t388 + t426;
t300 = pkin(3) * t339 + pkin(11) * t323;
t325 = pkin(3) * t362 + pkin(11) * t337;
t422 = -t300 * t365 + t325 * t352 + t425;
t301 = pkin(3) * t341 + pkin(11) * t324;
t421 = t300 * t353 - t301 * t352 + t424;
t420 = t301 * t365 - t325 * t353 + t423;
t297 = t339 * t452 + (t338 * t449 + t363 * t448) * t411;
t274 = pkin(4) * t297 + pkin(12) * t296;
t322 = t362 * t452 + (t361 * t449 + t375 * t448) * t411;
t289 = pkin(4) * t322 + pkin(12) * t321;
t419 = -t274 * t326 + t289 * t311 + t422;
t299 = t341 * t452 + (t340 * t449 + t364 * t448) * t411;
t275 = pkin(4) * t299 + pkin(12) * t298;
t418 = t274 * t312 - t275 * t311 + t421;
t417 = t275 * t326 - t289 * t312 + t420;
t414 = cos(qJ(6));
t410 = sin(qJ(5));
t409 = sin(qJ(6));
t401 = Icges(2,4) * t406;
t396 = rSges(2,1) * t406 - rSges(2,2) * t403;
t395 = rSges(2,1) * t403 + rSges(2,2) * t406;
t394 = Icges(2,1) * t406 - t447;
t393 = Icges(2,1) * t403 + t401;
t392 = -Icges(2,2) * t403 + t401;
t391 = Icges(2,2) * t406 + t447;
t386 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t385 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t384 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t373 = t408 * rSges(3,3) + (rSges(3,1) * t413 + rSges(3,2) * t416) * t405;
t372 = Icges(3,5) * t408 + (Icges(3,1) * t413 + Icges(3,4) * t416) * t405;
t371 = Icges(3,6) * t408 + (Icges(3,4) * t413 + Icges(3,2) * t416) * t405;
t370 = Icges(3,3) * t408 + (Icges(3,5) * t413 + Icges(3,6) * t416) * t405;
t368 = V_base(5) * rSges(2,3) - t395 * V_base(6) + t437;
t367 = t396 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t360 = t395 * V_base(4) - t396 * V_base(5) + t433;
t351 = rSges(3,1) * t379 + rSges(3,2) * t378 + rSges(3,3) * t446;
t350 = rSges(3,1) * t377 + rSges(3,2) * t376 - rSges(3,3) * t443;
t349 = Icges(3,1) * t379 + Icges(3,4) * t378 + Icges(3,5) * t446;
t348 = Icges(3,1) * t377 + Icges(3,4) * t376 - Icges(3,5) * t443;
t347 = Icges(3,4) * t379 + Icges(3,2) * t378 + Icges(3,6) * t446;
t346 = Icges(3,4) * t377 + Icges(3,2) * t376 - Icges(3,6) * t443;
t345 = Icges(3,5) * t379 + Icges(3,6) * t378 + Icges(3,3) * t446;
t344 = Icges(3,5) * t377 + Icges(3,6) * t376 - Icges(3,3) * t443;
t330 = rSges(4,1) * t362 + rSges(4,2) * t361 + rSges(4,3) * t375;
t329 = Icges(4,1) * t362 + Icges(4,4) * t361 + Icges(4,5) * t375;
t328 = Icges(4,4) * t362 + Icges(4,2) * t361 + Icges(4,6) * t375;
t327 = Icges(4,5) * t362 + Icges(4,6) * t361 + Icges(4,3) * t375;
t315 = -t350 * t397 + t373 * t387 + t428;
t314 = t351 * t397 - t373 * t388 + t426;
t310 = rSges(4,1) * t341 + rSges(4,2) * t340 + rSges(4,3) * t364;
t309 = rSges(4,1) * t339 + rSges(4,2) * t338 + rSges(4,3) * t363;
t308 = Icges(4,1) * t341 + Icges(4,4) * t340 + Icges(4,5) * t364;
t307 = Icges(4,1) * t339 + Icges(4,4) * t338 + Icges(4,5) * t363;
t306 = Icges(4,4) * t341 + Icges(4,2) * t340 + Icges(4,6) * t364;
t305 = Icges(4,4) * t339 + Icges(4,2) * t338 + Icges(4,6) * t363;
t304 = Icges(4,5) * t341 + Icges(4,6) * t340 + Icges(4,3) * t364;
t303 = Icges(4,5) * t339 + Icges(4,6) * t338 + Icges(4,3) * t363;
t302 = t350 * t388 - t351 * t387 + t427;
t295 = t322 * t451 + t337 * t410;
t294 = t322 * t410 - t337 * t451;
t287 = rSges(5,1) * t322 - rSges(5,2) * t321 + rSges(5,3) * t337;
t286 = Icges(5,1) * t322 - Icges(5,4) * t321 + Icges(5,5) * t337;
t285 = Icges(5,4) * t322 - Icges(5,2) * t321 + Icges(5,6) * t337;
t284 = Icges(5,5) * t322 - Icges(5,6) * t321 + Icges(5,3) * t337;
t283 = t299 * t451 + t324 * t410;
t282 = t299 * t410 - t324 * t451;
t281 = t297 * t451 + t323 * t410;
t280 = t297 * t410 - t323 * t451;
t279 = t295 * t414 + t321 * t409;
t278 = -t295 * t409 + t321 * t414;
t273 = pkin(5) * t295 + pkin(13) * t294;
t271 = qJD(6) * t294 + t290;
t269 = rSges(5,1) * t299 - rSges(5,2) * t298 + rSges(5,3) * t324;
t268 = rSges(5,1) * t297 - rSges(5,2) * t296 + rSges(5,3) * t323;
t267 = Icges(5,1) * t299 - Icges(5,4) * t298 + Icges(5,5) * t324;
t266 = Icges(5,1) * t297 - Icges(5,4) * t296 + Icges(5,5) * t323;
t265 = Icges(5,4) * t299 - Icges(5,2) * t298 + Icges(5,6) * t324;
t264 = Icges(5,4) * t297 - Icges(5,2) * t296 + Icges(5,6) * t323;
t263 = Icges(5,5) * t299 - Icges(5,6) * t298 + Icges(5,3) * t324;
t262 = Icges(5,5) * t297 - Icges(5,6) * t296 + Icges(5,3) * t323;
t261 = rSges(6,1) * t295 - rSges(6,2) * t294 + rSges(6,3) * t321;
t260 = Icges(6,1) * t295 - Icges(6,4) * t294 + Icges(6,5) * t321;
t259 = Icges(6,4) * t295 - Icges(6,2) * t294 + Icges(6,6) * t321;
t258 = Icges(6,5) * t295 - Icges(6,6) * t294 + Icges(6,3) * t321;
t257 = t283 * t414 + t298 * t409;
t256 = -t283 * t409 + t298 * t414;
t255 = t281 * t414 + t296 * t409;
t254 = -t281 * t409 + t296 * t414;
t253 = -t309 * t365 + t330 * t352 + t425;
t252 = t310 * t365 - t330 * t353 + t423;
t250 = pkin(5) * t283 + pkin(13) * t282;
t249 = pkin(5) * t281 + pkin(13) * t280;
t248 = t309 * t353 - t310 * t352 + t424;
t247 = qJD(6) * t282 + t277;
t246 = qJD(6) * t280 + t276;
t245 = rSges(6,1) * t283 - rSges(6,2) * t282 + rSges(6,3) * t298;
t244 = rSges(6,1) * t281 - rSges(6,2) * t280 + rSges(6,3) * t296;
t243 = Icges(6,1) * t283 - Icges(6,4) * t282 + Icges(6,5) * t298;
t242 = Icges(6,1) * t281 - Icges(6,4) * t280 + Icges(6,5) * t296;
t241 = Icges(6,4) * t283 - Icges(6,2) * t282 + Icges(6,6) * t298;
t240 = Icges(6,4) * t281 - Icges(6,2) * t280 + Icges(6,6) * t296;
t239 = Icges(6,5) * t283 - Icges(6,6) * t282 + Icges(6,3) * t298;
t238 = Icges(6,5) * t281 - Icges(6,6) * t280 + Icges(6,3) * t296;
t237 = rSges(7,1) * t279 + rSges(7,2) * t278 + rSges(7,3) * t294;
t236 = Icges(7,1) * t279 + Icges(7,4) * t278 + Icges(7,5) * t294;
t235 = Icges(7,4) * t279 + Icges(7,2) * t278 + Icges(7,6) * t294;
t234 = Icges(7,5) * t279 + Icges(7,6) * t278 + Icges(7,3) * t294;
t233 = rSges(7,1) * t257 + rSges(7,2) * t256 + rSges(7,3) * t282;
t232 = rSges(7,1) * t255 + rSges(7,2) * t254 + rSges(7,3) * t280;
t231 = Icges(7,1) * t257 + Icges(7,4) * t256 + Icges(7,5) * t282;
t230 = Icges(7,1) * t255 + Icges(7,4) * t254 + Icges(7,5) * t280;
t229 = Icges(7,4) * t257 + Icges(7,2) * t256 + Icges(7,6) * t282;
t228 = Icges(7,4) * t255 + Icges(7,2) * t254 + Icges(7,6) * t280;
t227 = Icges(7,5) * t257 + Icges(7,6) * t256 + Icges(7,3) * t282;
t226 = Icges(7,5) * t255 + Icges(7,6) * t254 + Icges(7,3) * t280;
t225 = -t268 * t326 + t287 * t311 + t422;
t224 = t269 * t326 - t287 * t312 + t420;
t223 = t268 * t312 - t269 * t311 + t421;
t222 = -t244 * t290 + t261 * t276 + t419;
t221 = t245 * t290 - t261 * t277 + t417;
t220 = t244 * t277 - t245 * t276 + t418;
t219 = -t232 * t271 + t237 * t246 - t249 * t290 + t273 * t276 + t419;
t218 = t233 * t271 - t237 * t247 + t250 * t290 - t273 * t277 + t417;
t217 = t232 * t247 - t233 * t246 + t249 * t277 - t250 * t276 + t418;
t1 = t326 * ((t263 * t337 - t265 * t321 + t267 * t322) * t312 + (t262 * t337 - t264 * t321 + t266 * t322) * t311 + (t284 * t337 - t285 * t321 + t286 * t322) * t326) / 0.2e1 + t311 * ((t263 * t323 - t265 * t296 + t267 * t297) * t312 + (t262 * t323 - t264 * t296 + t266 * t297) * t311 + (t284 * t323 - t285 * t296 + t286 * t297) * t326) / 0.2e1 + t312 * ((t263 * t324 - t298 * t265 + t299 * t267) * t312 + (t262 * t324 - t264 * t298 + t266 * t299) * t311 + (t284 * t324 - t285 * t298 + t286 * t299) * t326) / 0.2e1 + ((-t391 * t403 + t393 * t406 + Icges(1,4)) * V_base(5) + (-t392 * t403 + t394 * t406 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t391 * t406 + t393 * t403 + Icges(1,2)) * V_base(5) + (t392 * t406 + t394 * t403 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + m(3) * (t302 ^ 2 + t314 ^ 2 + t315 ^ 2) / 0.2e1 + t290 * ((t239 * t321 - t241 * t294 + t243 * t295) * t277 + (t238 * t321 - t240 * t294 + t242 * t295) * t276 + (t258 * t321 - t259 * t294 + t260 * t295) * t290) / 0.2e1 + m(4) * (t248 ^ 2 + t252 ^ 2 + t253 ^ 2) / 0.2e1 + ((Icges(2,5) * t403 + Icges(2,6) * t406 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t406 - Icges(2,6) * t403 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + m(1) * (t384 ^ 2 + t385 ^ 2 + t386 ^ 2) / 0.2e1 + m(5) * (t223 ^ 2 + t224 ^ 2 + t225 ^ 2) / 0.2e1 + m(6) * (t220 ^ 2 + t221 ^ 2 + t222 ^ 2) / 0.2e1 + t247 * ((t282 * t227 + t256 * t229 + t257 * t231) * t247 + (t226 * t282 + t228 * t256 + t230 * t257) * t246 + (t234 * t282 + t235 * t256 + t236 * t257) * t271) / 0.2e1 + t271 * ((t227 * t294 + t229 * t278 + t231 * t279) * t247 + (t226 * t294 + t228 * t278 + t230 * t279) * t246 + (t294 * t234 + t278 * t235 + t279 * t236) * t271) / 0.2e1 + t276 * ((t239 * t296 - t241 * t280 + t243 * t281) * t277 + (t296 * t238 - t280 * t240 + t281 * t242) * t276 + (t258 * t296 - t259 * t280 + t260 * t281) * t290) / 0.2e1 + t277 * ((t298 * t239 - t282 * t241 + t283 * t243) * t277 + (t238 * t298 - t240 * t282 + t242 * t283) * t276 + (t258 * t298 - t259 * t282 + t260 * t283) * t290) / 0.2e1 + t365 * ((t304 * t375 + t306 * t361 + t308 * t362) * t353 + (t303 * t375 + t305 * t361 + t307 * t362) * t352 + (t327 * t375 + t328 * t361 + t329 * t362) * t365) / 0.2e1 + t353 * ((t304 * t364 + t306 * t340 + t308 * t341) * t353 + (t303 * t364 + t305 * t340 + t307 * t341) * t352 + (t327 * t364 + t328 * t340 + t329 * t341) * t365) / 0.2e1 + t352 * ((t304 * t363 + t306 * t338 + t308 * t339) * t353 + (t303 * t363 + t305 * t338 + t307 * t339) * t352 + (t327 * t363 + t328 * t338 + t329 * t339) * t365) / 0.2e1 + m(2) * (t360 ^ 2 + t367 ^ 2 + t368 ^ 2) / 0.2e1 + m(7) * (t217 ^ 2 + t218 ^ 2 + t219 ^ 2) / 0.2e1 + t246 * ((t227 * t280 + t229 * t254 + t231 * t255) * t247 + (t280 * t226 + t254 * t228 + t255 * t230) * t246 + (t234 * t280 + t235 * t254 + t236 * t255) * t271) / 0.2e1 + t387 * ((-t345 * t443 + t347 * t376 + t349 * t377) * t388 + (-t344 * t443 + t346 * t376 + t348 * t377) * t387 + (-t370 * t443 + t371 * t376 + t372 * t377) * t397) / 0.2e1 + t388 * ((t345 * t446 + t347 * t378 + t349 * t379) * t388 + (t344 * t446 + t346 * t378 + t348 * t379) * t387 + (t370 * t446 + t371 * t378 + t372 * t379) * t397) / 0.2e1 + t397 * ((t344 * t387 + t345 * t388 + t370 * t397) * t408 + ((t347 * t416 + t349 * t413) * t388 + (t346 * t416 + t348 * t413) * t387 + (t371 * t416 + t372 * t413) * t397) * t405) / 0.2e1;
T  = t1;
