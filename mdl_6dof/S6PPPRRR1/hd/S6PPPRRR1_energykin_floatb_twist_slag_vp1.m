% Calculate kinetic energy for
% S6PPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
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
% Datum: 2019-03-08 18:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPPRRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PPPRRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_energykin_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPPRRR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPPRRR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:39:03
% EndTime: 2019-03-08 18:39:08
% DurationCPUTime: 4.67s
% Computational Cost: add. (9342->458), mult. (25999->646), div. (0->0), fcn. (34103->18), ass. (0->200)
t401 = sin(pkin(12));
t406 = cos(pkin(12));
t455 = Icges(2,5) * t406 - Icges(2,6) * t401 + Icges(1,5);
t454 = Icges(2,5) * t401 + Icges(2,6) * t406 + Icges(1,6);
t453 = cos(qJ(4));
t452 = cos(qJ(5));
t451 = cos(pkin(8));
t450 = sin(pkin(8));
t449 = Icges(2,4) * t401;
t408 = cos(pkin(6));
t448 = qJ(2) * t408;
t400 = sin(pkin(13));
t403 = sin(pkin(6));
t447 = t400 * t403;
t446 = t401 * t403;
t445 = t401 * t408;
t402 = sin(pkin(7));
t444 = t402 * t403;
t443 = t402 * t408;
t442 = t403 * t406;
t407 = cos(pkin(7));
t441 = t403 * t407;
t405 = cos(pkin(13));
t440 = t405 * t407;
t439 = t406 * t408;
t375 = -t400 * t401 + t405 * t439;
t362 = -t375 * t402 - t406 * t441;
t376 = t400 * t439 + t401 * t405;
t344 = pkin(2) * t376 + qJ(3) * t362;
t381 = pkin(1) * t401 - qJ(2) * t442;
t438 = -t344 - t381;
t377 = -t400 * t406 - t405 * t445;
t363 = -t377 * t402 + t401 * t441;
t378 = -t400 * t445 + t405 * t406;
t345 = pkin(2) * t378 + qJ(3) * t363;
t382 = pkin(1) * t406 + qJ(2) * t446;
t437 = -t345 - t382;
t436 = qJD(2) * t403;
t435 = V_base(5) * qJ(1) + V_base(1);
t431 = qJD(1) + V_base(3);
t399 = sin(pkin(14));
t404 = cos(pkin(14));
t423 = t377 * t407 + t401 * t444;
t339 = -t378 * t399 + t404 * t423;
t326 = -t339 * t450 + t363 * t451;
t316 = qJD(4) * t326 + V_base(4);
t424 = t375 * t407 - t402 * t442;
t337 = -t376 * t399 + t404 * t424;
t325 = -t337 * t450 + t362 * t451;
t315 = qJD(4) * t325 + V_base(5);
t359 = t404 * t443 + (-t399 * t400 + t404 * t440) * t403;
t374 = -t405 * t444 + t407 * t408;
t341 = -t359 * t450 + t374 * t451;
t334 = qJD(4) * t341 + V_base(6);
t430 = -qJ(1) - t448;
t340 = t378 * t404 + t399 * t423;
t411 = sin(qJ(4));
t428 = t453 * t450;
t429 = t451 * t453;
t298 = -t339 * t429 + t340 * t411 - t363 * t428;
t281 = qJD(5) * t298 + t316;
t338 = t376 * t404 + t399 * t424;
t296 = -t337 * t429 + t338 * t411 - t362 * t428;
t280 = qJD(5) * t296 + t315;
t360 = t404 * t447 + (t403 * t440 + t443) * t399;
t320 = -t359 * t429 + t360 * t411 - t374 * t428;
t293 = qJD(5) * t320 + t334;
t427 = t401 * t436 + V_base(5) * t448 + t435;
t365 = pkin(2) * t447 + qJ(3) * t374;
t426 = -t365 + t430;
t425 = qJD(2) * t408 + V_base(4) * t381 + t431;
t422 = V_base(6) * t382 - t406 * t436 + V_base(2);
t421 = qJD(3) * t363 + V_base(5) * t365 + t427;
t420 = qJD(3) * t374 + V_base(4) * t344 + t425;
t419 = qJD(3) * t362 + V_base(6) * t345 + t422;
t302 = t338 * pkin(3) + pkin(9) * t325;
t327 = t360 * pkin(3) + pkin(9) * t341;
t418 = V_base(5) * t327 + (-t302 + t438) * V_base(6) + t421;
t303 = t340 * pkin(3) + pkin(9) * t326;
t417 = V_base(4) * t302 + (-t303 + t437) * V_base(5) + t420;
t297 = t338 * t453 + (t337 * t451 + t362 * t450) * t411;
t275 = pkin(4) * t297 + pkin(10) * t296;
t321 = t360 * t453 + (t359 * t451 + t374 * t450) * t411;
t290 = pkin(4) * t321 + pkin(10) * t320;
t416 = -t275 * t334 + t315 * t290 + t418;
t299 = t340 * t453 + (t339 * t451 + t363 * t450) * t411;
t276 = pkin(4) * t299 + pkin(10) * t298;
t415 = t316 * t275 - t276 * t315 + t417;
t414 = V_base(6) * t303 + (-t327 + t426) * V_base(4) + t419;
t413 = t334 * t276 - t290 * t316 + t414;
t412 = cos(qJ(6));
t410 = sin(qJ(5));
t409 = sin(qJ(6));
t397 = Icges(2,4) * t406;
t393 = rSges(2,1) * t406 - rSges(2,2) * t401;
t392 = rSges(2,1) * t401 + rSges(2,2) * t406;
t391 = Icges(2,1) * t406 - t449;
t390 = Icges(2,1) * t401 + t397;
t389 = -Icges(2,2) * t401 + t397;
t388 = Icges(2,2) * t406 + t449;
t385 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t384 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t383 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t372 = rSges(3,3) * t408 + (rSges(3,1) * t400 + rSges(3,2) * t405) * t403;
t371 = Icges(3,5) * t408 + (Icges(3,1) * t400 + Icges(3,4) * t405) * t403;
t370 = Icges(3,6) * t408 + (Icges(3,4) * t400 + Icges(3,2) * t405) * t403;
t369 = Icges(3,3) * t408 + (Icges(3,5) * t400 + Icges(3,6) * t405) * t403;
t367 = V_base(5) * rSges(2,3) - t392 * V_base(6) + t435;
t366 = t393 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t361 = t392 * V_base(4) - t393 * V_base(5) + t431;
t353 = rSges(3,1) * t378 + rSges(3,2) * t377 + rSges(3,3) * t446;
t352 = rSges(3,1) * t376 + rSges(3,2) * t375 - rSges(3,3) * t442;
t351 = Icges(3,1) * t378 + Icges(3,4) * t377 + Icges(3,5) * t446;
t350 = Icges(3,1) * t376 + Icges(3,4) * t375 - Icges(3,5) * t442;
t349 = Icges(3,4) * t378 + Icges(3,2) * t377 + Icges(3,6) * t446;
t348 = Icges(3,4) * t376 + Icges(3,2) * t375 - Icges(3,6) * t442;
t347 = Icges(3,5) * t378 + Icges(3,6) * t377 + Icges(3,3) * t446;
t346 = Icges(3,5) * t376 + Icges(3,6) * t375 - Icges(3,3) * t442;
t331 = rSges(4,1) * t360 + rSges(4,2) * t359 + rSges(4,3) * t374;
t330 = Icges(4,1) * t360 + Icges(4,4) * t359 + Icges(4,5) * t374;
t329 = Icges(4,4) * t360 + Icges(4,2) * t359 + Icges(4,6) * t374;
t328 = Icges(4,5) * t360 + Icges(4,6) * t359 + Icges(4,3) * t374;
t314 = t372 * V_base(5) + (-t352 - t381) * V_base(6) + t427;
t313 = t353 * V_base(6) + (-t372 + t430) * V_base(4) + t422;
t312 = t352 * V_base(4) + (-t353 - t382) * V_base(5) + t425;
t311 = rSges(4,1) * t340 + rSges(4,2) * t339 + rSges(4,3) * t363;
t310 = rSges(4,1) * t338 + rSges(4,2) * t337 + rSges(4,3) * t362;
t309 = Icges(4,1) * t340 + Icges(4,4) * t339 + Icges(4,5) * t363;
t308 = Icges(4,1) * t338 + Icges(4,4) * t337 + Icges(4,5) * t362;
t307 = Icges(4,4) * t340 + Icges(4,2) * t339 + Icges(4,6) * t363;
t306 = Icges(4,4) * t338 + Icges(4,2) * t337 + Icges(4,6) * t362;
t305 = Icges(4,5) * t340 + Icges(4,6) * t339 + Icges(4,3) * t363;
t304 = Icges(4,5) * t338 + Icges(4,6) * t337 + Icges(4,3) * t362;
t295 = t321 * t452 + t341 * t410;
t294 = t321 * t410 - t341 * t452;
t289 = rSges(5,1) * t321 - rSges(5,2) * t320 + rSges(5,3) * t341;
t288 = Icges(5,1) * t321 - Icges(5,4) * t320 + Icges(5,5) * t341;
t287 = Icges(5,4) * t321 - Icges(5,2) * t320 + Icges(5,6) * t341;
t286 = Icges(5,5) * t321 - Icges(5,6) * t320 + Icges(5,3) * t341;
t285 = t299 * t452 + t326 * t410;
t284 = t299 * t410 - t326 * t452;
t283 = t297 * t452 + t325 * t410;
t282 = t297 * t410 - t325 * t452;
t279 = t295 * t412 + t320 * t409;
t278 = -t295 * t409 + t320 * t412;
t274 = qJD(6) * t294 + t293;
t273 = pkin(5) * t295 + pkin(11) * t294;
t271 = t331 * V_base(5) + (-t310 + t438) * V_base(6) + t421;
t270 = t311 * V_base(6) + (-t331 + t426) * V_base(4) + t419;
t269 = rSges(5,1) * t299 - rSges(5,2) * t298 + rSges(5,3) * t326;
t268 = rSges(5,1) * t297 - rSges(5,2) * t296 + rSges(5,3) * t325;
t267 = Icges(5,1) * t299 - Icges(5,4) * t298 + Icges(5,5) * t326;
t266 = Icges(5,1) * t297 - Icges(5,4) * t296 + Icges(5,5) * t325;
t265 = Icges(5,4) * t299 - Icges(5,2) * t298 + Icges(5,6) * t326;
t264 = Icges(5,4) * t297 - Icges(5,2) * t296 + Icges(5,6) * t325;
t263 = Icges(5,5) * t299 - Icges(5,6) * t298 + Icges(5,3) * t326;
t262 = Icges(5,5) * t297 - Icges(5,6) * t296 + Icges(5,3) * t325;
t261 = rSges(6,1) * t295 - rSges(6,2) * t294 + rSges(6,3) * t320;
t259 = Icges(6,1) * t295 - Icges(6,4) * t294 + Icges(6,5) * t320;
t258 = Icges(6,4) * t295 - Icges(6,2) * t294 + Icges(6,6) * t320;
t257 = Icges(6,5) * t295 - Icges(6,6) * t294 + Icges(6,3) * t320;
t256 = t285 * t412 + t298 * t409;
t255 = -t285 * t409 + t298 * t412;
t254 = t283 * t412 + t296 * t409;
t253 = -t283 * t409 + t296 * t412;
t252 = t310 * V_base(4) + (-t311 + t437) * V_base(5) + t420;
t251 = qJD(6) * t284 + t281;
t250 = qJD(6) * t282 + t280;
t249 = pkin(5) * t285 + pkin(11) * t284;
t248 = pkin(5) * t283 + pkin(11) * t282;
t247 = rSges(6,1) * t285 - rSges(6,2) * t284 + rSges(6,3) * t298;
t246 = rSges(6,1) * t283 - rSges(6,2) * t282 + rSges(6,3) * t296;
t245 = Icges(6,1) * t285 - Icges(6,4) * t284 + Icges(6,5) * t298;
t244 = Icges(6,1) * t283 - Icges(6,4) * t282 + Icges(6,5) * t296;
t243 = Icges(6,4) * t285 - Icges(6,2) * t284 + Icges(6,6) * t298;
t242 = Icges(6,4) * t283 - Icges(6,2) * t282 + Icges(6,6) * t296;
t241 = Icges(6,5) * t285 - Icges(6,6) * t284 + Icges(6,3) * t298;
t240 = Icges(6,5) * t283 - Icges(6,6) * t282 + Icges(6,3) * t296;
t239 = rSges(7,1) * t279 + rSges(7,2) * t278 + rSges(7,3) * t294;
t238 = Icges(7,1) * t279 + Icges(7,4) * t278 + Icges(7,5) * t294;
t237 = Icges(7,4) * t279 + Icges(7,2) * t278 + Icges(7,6) * t294;
t236 = Icges(7,5) * t279 + Icges(7,6) * t278 + Icges(7,3) * t294;
t235 = rSges(7,1) * t256 + rSges(7,2) * t255 + rSges(7,3) * t284;
t234 = rSges(7,1) * t254 + rSges(7,2) * t253 + rSges(7,3) * t282;
t233 = Icges(7,1) * t256 + Icges(7,4) * t255 + Icges(7,5) * t284;
t232 = Icges(7,1) * t254 + Icges(7,4) * t253 + Icges(7,5) * t282;
t231 = Icges(7,4) * t256 + Icges(7,2) * t255 + Icges(7,6) * t284;
t230 = Icges(7,4) * t254 + Icges(7,2) * t253 + Icges(7,6) * t282;
t229 = Icges(7,5) * t256 + Icges(7,6) * t255 + Icges(7,3) * t284;
t228 = Icges(7,5) * t254 + Icges(7,6) * t253 + Icges(7,3) * t282;
t227 = -t268 * t334 + t289 * t315 + t418;
t226 = t269 * t334 - t289 * t316 + t414;
t225 = t268 * t316 - t269 * t315 + t417;
t224 = -t246 * t293 + t261 * t280 + t416;
t223 = t247 * t293 - t261 * t281 + t413;
t222 = t246 * t281 - t247 * t280 + t415;
t221 = -t234 * t274 + t239 * t250 - t248 * t293 + t273 * t280 + t416;
t220 = t235 * t274 - t239 * t251 + t249 * t293 - t273 * t281 + t413;
t219 = t234 * t251 - t235 * t250 + t248 * t281 - t249 * t280 + t415;
t1 = m(7) * (t219 ^ 2 + t220 ^ 2 + t221 ^ 2) / 0.2e1 + m(5) * (t225 ^ 2 + t226 ^ 2 + t227 ^ 2) / 0.2e1 + m(6) * (t222 ^ 2 + t223 ^ 2 + t224 ^ 2) / 0.2e1 + m(1) * (t383 ^ 2 + t384 ^ 2 + t385 ^ 2) / 0.2e1 + m(2) * (t361 ^ 2 + t366 ^ 2 + t367 ^ 2) / 0.2e1 + t334 * ((t263 * t341 - t265 * t320 + t267 * t321) * t316 + (t262 * t341 - t264 * t320 + t266 * t321) * t315 + (t286 * t341 - t287 * t320 + t288 * t321) * t334) / 0.2e1 + t315 * ((t263 * t325 - t265 * t296 + t267 * t297) * t316 + (t262 * t325 - t264 * t296 + t266 * t297) * t315 + (t286 * t325 - t287 * t296 + t288 * t297) * t334) / 0.2e1 + t316 * ((t263 * t326 - t265 * t298 + t267 * t299) * t316 + (t262 * t326 - t264 * t298 + t266 * t299) * t315 + (t286 * t326 - t287 * t298 + t288 * t299) * t334) / 0.2e1 + t293 * ((t241 * t320 - t243 * t294 + t245 * t295) * t281 + (t240 * t320 - t242 * t294 + t244 * t295) * t280 + (t320 * t257 - t294 * t258 + t295 * t259) * t293) / 0.2e1 + m(3) * (t312 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + t281 * ((t298 * t241 - t284 * t243 + t285 * t245) * t281 + (t240 * t298 - t242 * t284 + t244 * t285) * t280 + (t257 * t298 - t258 * t284 + t259 * t285) * t293) / 0.2e1 + t274 * ((t229 * t294 + t231 * t278 + t233 * t279) * t251 + (t228 * t294 + t230 * t278 + t232 * t279) * t250 + (t294 * t236 + t278 * t237 + t279 * t238) * t274) / 0.2e1 + t280 * ((t241 * t296 - t243 * t282 + t245 * t283) * t281 + (t296 * t240 - t282 * t242 + t283 * t244) * t280 + (t257 * t296 - t258 * t282 + t259 * t283) * t293) / 0.2e1 + t251 * ((t284 * t229 + t255 * t231 + t256 * t233) * t251 + (t228 * t284 + t230 * t255 + t232 * t256) * t250 + (t236 * t284 + t237 * t255 + t238 * t256) * t274) / 0.2e1 + t250 * ((t229 * t282 + t231 * t253 + t233 * t254) * t251 + (t282 * t228 + t253 * t230 + t254 * t232) * t250 + (t236 * t282 + t237 * t253 + t238 * t254) * t274) / 0.2e1 + m(4) * (t252 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + ((t328 * t363 + t329 * t339 + t330 * t340 + t369 * t446 + t370 * t377 + t371 * t378 + t455) * V_base(6) + (t304 * t363 + t306 * t339 + t308 * t340 + t346 * t446 + t348 * t377 + t350 * t378 - t388 * t401 + t390 * t406 + Icges(1,4)) * V_base(5) + (t305 * t363 + t307 * t339 + t309 * t340 + t347 * t446 + t349 * t377 + t351 * t378 - t389 * t401 + t391 * t406 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t328 * t362 + t329 * t337 + t330 * t338 - t369 * t442 + t370 * t375 + t371 * t376 + t454) * V_base(6) + (t304 * t362 + t306 * t337 + t308 * t338 - t346 * t442 + t348 * t375 + t350 * t376 + t388 * t406 + t390 * t401 + Icges(1,2)) * V_base(5) + (t305 * t362 + t307 * t337 + t309 * t338 - t347 * t442 + t349 * t375 + t351 * t376 + t389 * t406 + t391 * t401 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t369 * t408 + (t370 * t405 + t371 * t400) * t403 + t328 * t374 + t329 * t359 + t330 * t360 + Icges(2,3) + Icges(1,3)) * V_base(6) + (t346 * t408 + (t348 * t405 + t350 * t400) * t403 + t304 * t374 + t306 * t359 + t308 * t360 + t454) * V_base(5) + (t347 * t408 + (t349 * t405 + t351 * t400) * t403 + t305 * t374 + t307 * t359 + t309 * t360 + t455) * V_base(4)) * V_base(6) / 0.2e1;
T  = t1;
