% Calculate kinetic energy for
% S6PRPRRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRRR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRPRRR7_energykin_floatb_twist_slag_vp1: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:50:30
% EndTime: 2019-03-08 20:50:35
% DurationCPUTime: 5.14s
% Computational Cost: add. (9477->453), mult. (26224->665), div. (0->0), fcn. (34328->18), ass. (0->197)
t452 = cos(qJ(4));
t451 = cos(qJ(5));
t409 = cos(pkin(6));
t450 = pkin(9) * t409;
t449 = cos(pkin(8));
t448 = sin(pkin(8));
t403 = sin(pkin(13));
t447 = Icges(2,4) * t403;
t405 = sin(pkin(6));
t446 = t403 * t405;
t404 = sin(pkin(7));
t445 = t404 * t405;
t444 = t404 * t409;
t407 = cos(pkin(13));
t443 = t405 * t407;
t408 = cos(pkin(7));
t442 = t405 * t408;
t413 = sin(qJ(2));
t441 = t405 * t413;
t415 = cos(qJ(2));
t440 = t408 * t415;
t439 = t409 * t413;
t438 = t409 * t415;
t437 = qJD(2) * t405;
t436 = V_base(5) * qJ(1) + V_base(1);
t432 = qJD(1) + V_base(3);
t387 = t403 * t437 + V_base(4);
t396 = qJD(2) * t409 + V_base(6);
t378 = -t403 * t439 + t407 * t415;
t402 = sin(pkin(14));
t406 = cos(pkin(14));
t377 = -t403 * t438 - t407 * t413;
t428 = t377 * t408 + t403 * t445;
t342 = -t378 * t402 + t406 * t428;
t364 = -t377 * t404 + t403 * t442;
t326 = -t342 * t448 + t364 * t449;
t316 = qJD(4) * t326 + t387;
t360 = t406 * t444 + (-t402 * t413 + t406 * t440) * t405;
t374 = t409 * t408 - t415 * t445;
t339 = -t360 * t448 + t374 * t449;
t332 = qJD(4) * t339 + t396;
t431 = t449 * t452;
t430 = t452 * t448;
t343 = t378 * t406 + t402 * t428;
t412 = sin(qJ(4));
t300 = -t342 * t431 + t343 * t412 - t364 * t430;
t279 = qJD(5) * t300 + t316;
t361 = t406 * t441 + (t405 * t440 + t444) * t402;
t321 = -t360 * t431 + t361 * t412 - t374 * t430;
t291 = qJD(5) * t321 + t332;
t386 = -t407 * t437 + V_base(5);
t375 = -t403 * t413 + t407 * t438;
t363 = -t375 * t404 - t407 * t442;
t429 = t375 * t408 - t404 * t443;
t376 = t403 * t415 + t407 * t439;
t340 = -t376 * t402 + t406 * t429;
t325 = -t340 * t448 + t363 * t449;
t315 = qJD(4) * t325 + t386;
t381 = pkin(1) * t403 - pkin(9) * t443;
t427 = -t381 * V_base(6) + V_base(5) * t450 + t436;
t382 = pkin(1) * t407 + pkin(9) * t446;
t426 = V_base(4) * t381 - t382 * V_base(5) + t432;
t341 = t376 * t406 + t402 * t429;
t298 = -t340 * t431 + t341 * t412 - t363 * t430;
t278 = qJD(5) * t298 + t315;
t425 = V_base(6) * t382 + V_base(2) + (-qJ(1) - t450) * V_base(4);
t365 = pkin(2) * t441 + qJ(3) * t374;
t424 = qJD(3) * t364 + t386 * t365 + t427;
t344 = pkin(2) * t376 + qJ(3) * t363;
t423 = qJD(3) * t374 + t387 * t344 + t426;
t345 = pkin(2) * t378 + qJ(3) * t364;
t422 = qJD(3) * t363 + t396 * t345 + t425;
t302 = t341 * pkin(3) + pkin(10) * t325;
t327 = t361 * pkin(3) + pkin(10) * t339;
t421 = t386 * t327 + (-t302 - t344) * t396 + t424;
t303 = t343 * pkin(3) + pkin(10) * t326;
t420 = t387 * t302 + (-t303 - t345) * t386 + t423;
t419 = t396 * t303 + (-t327 - t365) * t387 + t422;
t299 = t341 * t452 + (t340 * t449 + t363 * t448) * t412;
t275 = pkin(4) * t299 + pkin(11) * t298;
t322 = t361 * t452 + (t360 * t449 + t374 * t448) * t412;
t290 = pkin(4) * t322 + pkin(11) * t321;
t418 = -t275 * t332 + t315 * t290 + t421;
t301 = t343 * t452 + (t342 * t449 + t364 * t448) * t412;
t276 = pkin(4) * t301 + pkin(11) * t300;
t417 = t316 * t275 - t276 * t315 + t420;
t416 = t332 * t276 - t290 * t316 + t419;
t414 = cos(qJ(6));
t411 = sin(qJ(5));
t410 = sin(qJ(6));
t400 = Icges(2,4) * t407;
t395 = rSges(2,1) * t407 - rSges(2,2) * t403;
t394 = rSges(2,1) * t403 + rSges(2,2) * t407;
t393 = Icges(2,1) * t407 - t447;
t392 = Icges(2,1) * t403 + t400;
t391 = -Icges(2,2) * t403 + t400;
t390 = Icges(2,2) * t407 + t447;
t385 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t384 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t383 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t372 = t409 * rSges(3,3) + (rSges(3,1) * t413 + rSges(3,2) * t415) * t405;
t371 = Icges(3,5) * t409 + (Icges(3,1) * t413 + Icges(3,4) * t415) * t405;
t370 = Icges(3,6) * t409 + (Icges(3,4) * t413 + Icges(3,2) * t415) * t405;
t369 = Icges(3,3) * t409 + (Icges(3,5) * t413 + Icges(3,6) * t415) * t405;
t367 = V_base(5) * rSges(2,3) - t394 * V_base(6) + t436;
t366 = t395 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t362 = t394 * V_base(4) - t395 * V_base(5) + t432;
t353 = rSges(3,1) * t378 + rSges(3,2) * t377 + rSges(3,3) * t446;
t352 = rSges(3,1) * t376 + rSges(3,2) * t375 - rSges(3,3) * t443;
t351 = Icges(3,1) * t378 + Icges(3,4) * t377 + Icges(3,5) * t446;
t350 = Icges(3,1) * t376 + Icges(3,4) * t375 - Icges(3,5) * t443;
t349 = Icges(3,4) * t378 + Icges(3,2) * t377 + Icges(3,6) * t446;
t348 = Icges(3,4) * t376 + Icges(3,2) * t375 - Icges(3,6) * t443;
t347 = Icges(3,5) * t378 + Icges(3,6) * t377 + Icges(3,3) * t446;
t346 = Icges(3,5) * t376 + Icges(3,6) * t375 - Icges(3,3) * t443;
t331 = rSges(4,1) * t361 + rSges(4,2) * t360 + rSges(4,3) * t374;
t330 = Icges(4,1) * t361 + Icges(4,4) * t360 + Icges(4,5) * t374;
t329 = Icges(4,4) * t361 + Icges(4,2) * t360 + Icges(4,6) * t374;
t328 = Icges(4,5) * t361 + Icges(4,6) * t360 + Icges(4,3) * t374;
t314 = -t352 * t396 + t372 * t386 + t427;
t313 = t353 * t396 - t372 * t387 + t425;
t312 = rSges(4,1) * t343 + rSges(4,2) * t342 + rSges(4,3) * t364;
t311 = rSges(4,1) * t341 + rSges(4,2) * t340 + rSges(4,3) * t363;
t310 = Icges(4,1) * t343 + Icges(4,4) * t342 + Icges(4,5) * t364;
t309 = Icges(4,1) * t341 + Icges(4,4) * t340 + Icges(4,5) * t363;
t308 = Icges(4,4) * t343 + Icges(4,2) * t342 + Icges(4,6) * t364;
t307 = Icges(4,4) * t341 + Icges(4,2) * t340 + Icges(4,6) * t363;
t306 = Icges(4,5) * t343 + Icges(4,6) * t342 + Icges(4,3) * t364;
t305 = Icges(4,5) * t341 + Icges(4,6) * t340 + Icges(4,3) * t363;
t304 = t352 * t387 - t353 * t386 + t426;
t297 = t322 * t451 + t339 * t411;
t296 = t322 * t411 - t339 * t451;
t289 = rSges(5,1) * t322 - rSges(5,2) * t321 + rSges(5,3) * t339;
t288 = Icges(5,1) * t322 - Icges(5,4) * t321 + Icges(5,5) * t339;
t287 = Icges(5,4) * t322 - Icges(5,2) * t321 + Icges(5,6) * t339;
t286 = Icges(5,5) * t322 - Icges(5,6) * t321 + Icges(5,3) * t339;
t285 = t301 * t451 + t326 * t411;
t284 = t301 * t411 - t326 * t451;
t283 = t299 * t451 + t325 * t411;
t282 = t299 * t411 - t325 * t451;
t281 = t297 * t414 + t321 * t410;
t280 = -t297 * t410 + t321 * t414;
t274 = pkin(5) * t297 + pkin(12) * t296;
t273 = qJD(6) * t296 + t291;
t271 = t331 * t386 + (-t311 - t344) * t396 + t424;
t270 = t312 * t396 + (-t331 - t365) * t387 + t422;
t269 = rSges(5,1) * t301 - rSges(5,2) * t300 + rSges(5,3) * t326;
t268 = rSges(5,1) * t299 - rSges(5,2) * t298 + rSges(5,3) * t325;
t267 = Icges(5,1) * t301 - Icges(5,4) * t300 + Icges(5,5) * t326;
t266 = Icges(5,1) * t299 - Icges(5,4) * t298 + Icges(5,5) * t325;
t265 = Icges(5,4) * t301 - Icges(5,2) * t300 + Icges(5,6) * t326;
t264 = Icges(5,4) * t299 - Icges(5,2) * t298 + Icges(5,6) * t325;
t263 = Icges(5,5) * t301 - Icges(5,6) * t300 + Icges(5,3) * t326;
t262 = Icges(5,5) * t299 - Icges(5,6) * t298 + Icges(5,3) * t325;
t261 = rSges(6,1) * t297 - rSges(6,2) * t296 + rSges(6,3) * t321;
t260 = Icges(6,1) * t297 - Icges(6,4) * t296 + Icges(6,5) * t321;
t259 = Icges(6,4) * t297 - Icges(6,2) * t296 + Icges(6,6) * t321;
t258 = Icges(6,5) * t297 - Icges(6,6) * t296 + Icges(6,3) * t321;
t257 = t285 * t414 + t300 * t410;
t256 = -t285 * t410 + t300 * t414;
t255 = t283 * t414 + t298 * t410;
t254 = -t283 * t410 + t298 * t414;
t252 = t311 * t387 + (-t312 - t345) * t386 + t423;
t251 = pkin(5) * t285 + pkin(12) * t284;
t250 = pkin(5) * t283 + pkin(12) * t282;
t249 = qJD(6) * t284 + t279;
t248 = qJD(6) * t282 + t278;
t247 = rSges(6,1) * t285 - rSges(6,2) * t284 + rSges(6,3) * t300;
t246 = rSges(6,1) * t283 - rSges(6,2) * t282 + rSges(6,3) * t298;
t245 = Icges(6,1) * t285 - Icges(6,4) * t284 + Icges(6,5) * t300;
t244 = Icges(6,1) * t283 - Icges(6,4) * t282 + Icges(6,5) * t298;
t243 = Icges(6,4) * t285 - Icges(6,2) * t284 + Icges(6,6) * t300;
t242 = Icges(6,4) * t283 - Icges(6,2) * t282 + Icges(6,6) * t298;
t241 = Icges(6,5) * t285 - Icges(6,6) * t284 + Icges(6,3) * t300;
t240 = Icges(6,5) * t283 - Icges(6,6) * t282 + Icges(6,3) * t298;
t239 = rSges(7,1) * t281 + rSges(7,2) * t280 + rSges(7,3) * t296;
t238 = Icges(7,1) * t281 + Icges(7,4) * t280 + Icges(7,5) * t296;
t237 = Icges(7,4) * t281 + Icges(7,2) * t280 + Icges(7,6) * t296;
t236 = Icges(7,5) * t281 + Icges(7,6) * t280 + Icges(7,3) * t296;
t235 = rSges(7,1) * t257 + rSges(7,2) * t256 + rSges(7,3) * t284;
t234 = rSges(7,1) * t255 + rSges(7,2) * t254 + rSges(7,3) * t282;
t233 = Icges(7,1) * t257 + Icges(7,4) * t256 + Icges(7,5) * t284;
t232 = Icges(7,1) * t255 + Icges(7,4) * t254 + Icges(7,5) * t282;
t231 = Icges(7,4) * t257 + Icges(7,2) * t256 + Icges(7,6) * t284;
t230 = Icges(7,4) * t255 + Icges(7,2) * t254 + Icges(7,6) * t282;
t229 = Icges(7,5) * t257 + Icges(7,6) * t256 + Icges(7,3) * t284;
t228 = Icges(7,5) * t255 + Icges(7,6) * t254 + Icges(7,3) * t282;
t227 = -t268 * t332 + t289 * t315 + t421;
t226 = t269 * t332 - t289 * t316 + t419;
t225 = t268 * t316 - t269 * t315 + t420;
t224 = -t246 * t291 + t261 * t278 + t418;
t223 = t247 * t291 - t261 * t279 + t416;
t222 = t246 * t279 - t247 * t278 + t417;
t221 = -t234 * t273 + t239 * t248 - t250 * t291 + t274 * t278 + t418;
t220 = t235 * t273 - t239 * t249 + t251 * t291 - t274 * t279 + t416;
t219 = t234 * t249 - t235 * t248 + t250 * t279 - t251 * t278 + t417;
t1 = t315 * ((t263 * t325 - t265 * t298 + t267 * t299) * t316 + (t262 * t325 - t264 * t298 + t266 * t299) * t315 + (t286 * t325 - t287 * t298 + t288 * t299) * t332) / 0.2e1 + t316 * ((t263 * t326 - t265 * t300 + t267 * t301) * t316 + (t262 * t326 - t264 * t300 + t266 * t301) * t315 + (t286 * t326 - t287 * t300 + t288 * t301) * t332) / 0.2e1 + t332 * ((t263 * t339 - t265 * t321 + t267 * t322) * t316 + (t262 * t339 - t264 * t321 + t266 * t322) * t315 + (t286 * t339 - t287 * t321 + t288 * t322) * t332) / 0.2e1 + m(5) * (t225 ^ 2 + t226 ^ 2 + t227 ^ 2) / 0.2e1 + m(4) * (t252 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + t291 * ((t241 * t321 - t243 * t296 + t245 * t297) * t279 + (t240 * t321 - t242 * t296 + t244 * t297) * t278 + (t321 * t258 - t296 * t259 + t297 * t260) * t291) / 0.2e1 + m(3) * (t304 ^ 2 + t313 ^ 2 + t314 ^ 2) / 0.2e1 + t278 * ((t241 * t298 - t243 * t282 + t245 * t283) * t279 + (t298 * t240 - t282 * t242 + t283 * t244) * t278 + (t258 * t298 - t259 * t282 + t260 * t283) * t291) / 0.2e1 + t279 * ((t300 * t241 - t284 * t243 + t285 * t245) * t279 + (t240 * t300 - t242 * t284 + t244 * t285) * t278 + (t258 * t300 - t259 * t284 + t260 * t285) * t291) / 0.2e1 + t273 * ((t229 * t296 + t231 * t280 + t233 * t281) * t249 + (t228 * t296 + t230 * t280 + t232 * t281) * t248 + (t296 * t236 + t280 * t237 + t281 * t238) * t273) / 0.2e1 + m(1) * (t383 ^ 2 + t384 ^ 2 + t385 ^ 2) / 0.2e1 + m(6) * (t222 ^ 2 + t223 ^ 2 + t224 ^ 2) / 0.2e1 + m(7) * (t219 ^ 2 + t220 ^ 2 + t221 ^ 2) / 0.2e1 + m(2) * (t362 ^ 2 + t366 ^ 2 + t367 ^ 2) / 0.2e1 + t249 * ((t284 * t229 + t256 * t231 + t257 * t233) * t249 + (t228 * t284 + t230 * t256 + t232 * t257) * t248 + (t236 * t284 + t237 * t256 + t238 * t257) * t273) / 0.2e1 + t248 * ((t229 * t282 + t231 * t254 + t233 * t255) * t249 + (t282 * t228 + t254 * t230 + t255 * t232) * t248 + (t236 * t282 + t237 * t254 + t238 * t255) * t273) / 0.2e1 + ((t328 * t363 + t329 * t340 + t330 * t341 - t369 * t443 + t370 * t375 + t371 * t376) * t396 + (t306 * t363 + t308 * t340 + t310 * t341 - t347 * t443 + t349 * t375 + t351 * t376) * t387 + (t305 * t363 + t307 * t340 + t309 * t341 - t346 * t443 + t348 * t375 + t350 * t376) * t386) * t386 / 0.2e1 + ((t328 * t364 + t329 * t342 + t330 * t343 + t369 * t446 + t370 * t377 + t371 * t378) * t396 + (t306 * t364 + t308 * t342 + t310 * t343 + t347 * t446 + t349 * t377 + t351 * t378) * t387 + (t305 * t364 + t307 * t342 + t309 * t343 + t346 * t446 + t348 * t377 + t350 * t378) * t386) * t387 / 0.2e1 + ((t346 * t386 + t347 * t387 + t369 * t396) * t409 + ((t349 * t415 + t351 * t413) * t387 + (t348 * t415 + t350 * t413) * t386 + (t370 * t415 + t371 * t413) * t396) * t405 + (t306 * t374 + t308 * t360 + t310 * t361) * t387 + (t305 * t374 + t307 * t360 + t309 * t361) * t386 + (t328 * t374 + t329 * t360 + t330 * t361) * t396) * t396 / 0.2e1 + ((-t390 * t403 + t392 * t407 + Icges(1,4)) * V_base(5) + (-t391 * t403 + t393 * t407 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t390 * t407 + t392 * t403 + Icges(1,2)) * V_base(5) + (t391 * t407 + t393 * t403 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t407 - Icges(2,6) * t403 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t403 + Icges(2,6) * t407 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
