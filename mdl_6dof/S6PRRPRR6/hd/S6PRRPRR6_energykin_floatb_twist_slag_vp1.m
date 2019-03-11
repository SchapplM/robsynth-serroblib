% Calculate kinetic energy for
% S6PRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:22:45
% EndTime: 2019-03-08 22:22:50
% DurationCPUTime: 5.02s
% Computational Cost: add. (5320->442), mult. (12648->643), div. (0->0), fcn. (16040->16), ass. (0->192)
t426 = Icges(4,2) + Icges(5,3);
t367 = sin(pkin(12));
t370 = cos(pkin(12));
t375 = sin(qJ(2));
t371 = cos(pkin(6));
t377 = cos(qJ(2));
t404 = t371 * t377;
t334 = -t367 * t375 + t370 * t404;
t405 = t371 * t375;
t335 = t367 * t377 + t370 * t405;
t374 = sin(qJ(3));
t368 = sin(pkin(6));
t413 = sin(pkin(7));
t394 = t368 * t413;
t414 = cos(pkin(7));
t417 = cos(qJ(3));
t294 = t335 * t417 + (t334 * t414 - t370 * t394) * t374;
t395 = t368 * t414;
t319 = -t334 * t413 - t370 * t395;
t366 = sin(pkin(13));
t369 = cos(pkin(13));
t273 = -t294 * t366 + t319 * t369;
t411 = t319 * t366;
t274 = t294 * t369 + t411;
t391 = t417 * t413;
t390 = t368 * t391;
t392 = t414 * t417;
t293 = -t334 * t392 + t335 * t374 + t370 * t390;
t425 = -Icges(4,4) * t294 + Icges(5,5) * t274 - Icges(4,6) * t319 + Icges(5,6) * t273 + t426 * t293;
t336 = -t367 * t404 - t370 * t375;
t337 = -t367 * t405 + t370 * t377;
t296 = t337 * t417 + (t336 * t414 + t367 * t394) * t374;
t320 = -t336 * t413 + t367 * t395;
t275 = -t296 * t366 + t320 * t369;
t410 = t320 * t366;
t276 = t296 * t369 + t410;
t295 = -t336 * t392 + t337 * t374 - t367 * t390;
t424 = -Icges(4,4) * t296 + Icges(5,5) * t276 - Icges(4,6) * t320 + Icges(5,6) * t275 + t426 * t295;
t318 = t371 * t413 * t374 + (t374 * t377 * t414 + t375 * t417) * t368;
t333 = t371 * t414 - t377 * t394;
t291 = -t318 * t366 + t333 * t369;
t409 = t333 * t366;
t292 = t318 * t369 + t409;
t406 = t368 * t375;
t317 = -t368 * t377 * t392 - t371 * t391 + t374 * t406;
t423 = -Icges(4,4) * t318 + Icges(5,5) * t292 - Icges(4,6) * t333 + Icges(5,6) * t291 + t426 * t317;
t416 = pkin(8) * t371;
t415 = pkin(4) * t369;
t412 = Icges(2,4) * t367;
t408 = t367 * t368;
t407 = t368 * t370;
t402 = qJD(2) * t368;
t401 = pkin(13) + qJ(5);
t400 = V_base(5) * qJ(1) + V_base(1);
t396 = qJD(1) + V_base(3);
t348 = t367 * t402 + V_base(4);
t358 = qJD(2) * t371 + V_base(6);
t393 = cos(t401);
t308 = qJD(3) * t320 + t348;
t321 = qJD(3) * t333 + t358;
t266 = qJD(5) * t295 + t308;
t283 = qJD(5) * t317 + t321;
t347 = -t370 * t402 + V_base(5);
t307 = qJD(3) * t319 + t347;
t340 = pkin(1) * t367 - pkin(8) * t407;
t389 = -t340 * V_base(6) + V_base(5) * t416 + t400;
t341 = pkin(1) * t370 + pkin(8) * t408;
t388 = V_base(4) * t340 - t341 * V_base(5) + t396;
t265 = qJD(5) * t293 + t307;
t387 = V_base(6) * t341 + V_base(2) + (-qJ(1) - t416) * V_base(4);
t297 = t335 * pkin(2) + pkin(9) * t319;
t322 = pkin(2) * t406 + pkin(9) * t333;
t386 = -t297 * t358 + t347 * t322 + t389;
t298 = t337 * pkin(2) + pkin(9) * t320;
t385 = t348 * t297 - t298 * t347 + t388;
t384 = t358 * t298 - t322 * t348 + t387;
t281 = pkin(3) * t318 + qJ(4) * t317;
t383 = qJD(4) * t295 + t307 * t281 + t386;
t262 = pkin(3) * t294 + qJ(4) * t293;
t382 = qJD(4) * t317 + t308 * t262 + t385;
t263 = pkin(3) * t296 + qJ(4) * t295;
t381 = qJD(4) * t293 + t321 * t263 + t384;
t207 = pkin(4) * t411 + pkin(10) * t293 + t294 * t415;
t238 = pkin(4) * t409 + pkin(10) * t317 + t318 * t415;
t380 = t307 * t238 + (-t207 - t262) * t321 + t383;
t208 = pkin(4) * t410 + pkin(10) * t295 + t296 * t415;
t379 = t308 * t207 + (-t208 - t263) * t307 + t382;
t378 = t321 * t208 + (-t238 - t281) * t308 + t381;
t376 = cos(qJ(6));
t373 = sin(qJ(6));
t364 = Icges(2,4) * t370;
t363 = sin(t401);
t356 = rSges(2,1) * t370 - rSges(2,2) * t367;
t355 = rSges(2,1) * t367 + rSges(2,2) * t370;
t354 = Icges(2,1) * t370 - t412;
t353 = Icges(2,1) * t367 + t364;
t352 = -Icges(2,2) * t367 + t364;
t351 = Icges(2,2) * t370 + t412;
t346 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t345 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t344 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t331 = t371 * rSges(3,3) + (rSges(3,1) * t375 + rSges(3,2) * t377) * t368;
t330 = Icges(3,5) * t371 + (Icges(3,1) * t375 + Icges(3,4) * t377) * t368;
t329 = Icges(3,6) * t371 + (Icges(3,4) * t375 + Icges(3,2) * t377) * t368;
t328 = Icges(3,3) * t371 + (Icges(3,5) * t375 + Icges(3,6) * t377) * t368;
t324 = V_base(5) * rSges(2,3) - t355 * V_base(6) + t400;
t323 = t356 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t316 = t355 * V_base(4) - t356 * V_base(5) + t396;
t306 = rSges(3,1) * t337 + rSges(3,2) * t336 + rSges(3,3) * t408;
t305 = rSges(3,1) * t335 + rSges(3,2) * t334 - rSges(3,3) * t407;
t304 = Icges(3,1) * t337 + Icges(3,4) * t336 + Icges(3,5) * t408;
t303 = Icges(3,1) * t335 + Icges(3,4) * t334 - Icges(3,5) * t407;
t302 = Icges(3,4) * t337 + Icges(3,2) * t336 + Icges(3,6) * t408;
t301 = Icges(3,4) * t335 + Icges(3,2) * t334 - Icges(3,6) * t407;
t300 = Icges(3,5) * t337 + Icges(3,6) * t336 + Icges(3,3) * t408;
t299 = Icges(3,5) * t335 + Icges(3,6) * t334 - Icges(3,3) * t407;
t285 = t318 * t393 + t333 * t363;
t284 = t318 * t363 - t333 * t393;
t280 = rSges(4,1) * t318 - rSges(4,2) * t317 + rSges(4,3) * t333;
t279 = Icges(4,1) * t318 - Icges(4,4) * t317 + Icges(4,5) * t333;
t277 = Icges(4,5) * t318 - Icges(4,6) * t317 + Icges(4,3) * t333;
t272 = t296 * t393 + t320 * t363;
t271 = t296 * t363 - t320 * t393;
t270 = t294 * t393 + t319 * t363;
t269 = t294 * t363 - t319 * t393;
t268 = t285 * t376 + t317 * t373;
t267 = -t285 * t373 + t317 * t376;
t261 = -t305 * t358 + t331 * t347 + t389;
t260 = t306 * t358 - t331 * t348 + t387;
t259 = qJD(6) * t284 + t283;
t258 = pkin(5) * t285 + pkin(11) * t284;
t256 = rSges(4,1) * t296 - rSges(4,2) * t295 + rSges(4,3) * t320;
t255 = rSges(4,1) * t294 - rSges(4,2) * t293 + rSges(4,3) * t319;
t254 = Icges(4,1) * t296 - Icges(4,4) * t295 + Icges(4,5) * t320;
t253 = Icges(4,1) * t294 - Icges(4,4) * t293 + Icges(4,5) * t319;
t250 = Icges(4,5) * t296 - Icges(4,6) * t295 + Icges(4,3) * t320;
t249 = Icges(4,5) * t294 - Icges(4,6) * t293 + Icges(4,3) * t319;
t248 = rSges(5,1) * t292 + rSges(5,2) * t291 + rSges(5,3) * t317;
t247 = Icges(5,1) * t292 + Icges(5,4) * t291 + Icges(5,5) * t317;
t246 = Icges(5,4) * t292 + Icges(5,2) * t291 + Icges(5,6) * t317;
t244 = t305 * t348 - t306 * t347 + t388;
t243 = rSges(6,1) * t285 - rSges(6,2) * t284 + rSges(6,3) * t317;
t242 = Icges(6,1) * t285 - Icges(6,4) * t284 + Icges(6,5) * t317;
t241 = Icges(6,4) * t285 - Icges(6,2) * t284 + Icges(6,6) * t317;
t240 = Icges(6,5) * t285 - Icges(6,6) * t284 + Icges(6,3) * t317;
t237 = t272 * t376 + t295 * t373;
t236 = -t272 * t373 + t295 * t376;
t235 = t270 * t376 + t293 * t373;
t234 = -t270 * t373 + t293 * t376;
t233 = pkin(5) * t272 + pkin(11) * t271;
t232 = pkin(5) * t270 + pkin(11) * t269;
t231 = qJD(6) * t271 + t266;
t230 = qJD(6) * t269 + t265;
t228 = rSges(5,1) * t276 + rSges(5,2) * t275 + rSges(5,3) * t295;
t227 = rSges(5,1) * t274 + rSges(5,2) * t273 + rSges(5,3) * t293;
t226 = Icges(5,1) * t276 + Icges(5,4) * t275 + Icges(5,5) * t295;
t225 = Icges(5,1) * t274 + Icges(5,4) * t273 + Icges(5,5) * t293;
t224 = Icges(5,4) * t276 + Icges(5,2) * t275 + Icges(5,6) * t295;
t223 = Icges(5,4) * t274 + Icges(5,2) * t273 + Icges(5,6) * t293;
t220 = rSges(6,1) * t272 - rSges(6,2) * t271 + rSges(6,3) * t295;
t219 = rSges(6,1) * t270 - rSges(6,2) * t269 + rSges(6,3) * t293;
t218 = Icges(6,1) * t272 - Icges(6,4) * t271 + Icges(6,5) * t295;
t217 = Icges(6,1) * t270 - Icges(6,4) * t269 + Icges(6,5) * t293;
t216 = Icges(6,4) * t272 - Icges(6,2) * t271 + Icges(6,6) * t295;
t215 = Icges(6,4) * t270 - Icges(6,2) * t269 + Icges(6,6) * t293;
t214 = Icges(6,5) * t272 - Icges(6,6) * t271 + Icges(6,3) * t295;
t213 = Icges(6,5) * t270 - Icges(6,6) * t269 + Icges(6,3) * t293;
t212 = rSges(7,1) * t268 + rSges(7,2) * t267 + rSges(7,3) * t284;
t211 = Icges(7,1) * t268 + Icges(7,4) * t267 + Icges(7,5) * t284;
t210 = Icges(7,4) * t268 + Icges(7,2) * t267 + Icges(7,6) * t284;
t209 = Icges(7,5) * t268 + Icges(7,6) * t267 + Icges(7,3) * t284;
t204 = rSges(7,1) * t237 + rSges(7,2) * t236 + rSges(7,3) * t271;
t203 = rSges(7,1) * t235 + rSges(7,2) * t234 + rSges(7,3) * t269;
t202 = Icges(7,1) * t237 + Icges(7,4) * t236 + Icges(7,5) * t271;
t201 = Icges(7,1) * t235 + Icges(7,4) * t234 + Icges(7,5) * t269;
t200 = Icges(7,4) * t237 + Icges(7,2) * t236 + Icges(7,6) * t271;
t199 = Icges(7,4) * t235 + Icges(7,2) * t234 + Icges(7,6) * t269;
t198 = Icges(7,5) * t237 + Icges(7,6) * t236 + Icges(7,3) * t271;
t197 = Icges(7,5) * t235 + Icges(7,6) * t234 + Icges(7,3) * t269;
t196 = -t255 * t321 + t280 * t307 + t386;
t195 = t256 * t321 - t280 * t308 + t384;
t194 = t255 * t308 - t256 * t307 + t385;
t193 = t248 * t307 + (-t227 - t262) * t321 + t383;
t192 = t228 * t321 + (-t248 - t281) * t308 + t381;
t191 = t227 * t308 + (-t228 - t263) * t307 + t382;
t190 = -t219 * t283 + t243 * t265 + t380;
t189 = t220 * t283 - t243 * t266 + t378;
t188 = t219 * t266 - t220 * t265 + t379;
t187 = -t203 * t259 + t212 * t230 - t232 * t283 + t258 * t265 + t380;
t186 = t204 * t259 - t212 * t231 + t233 * t283 - t258 * t266 + t378;
t185 = t203 * t231 - t204 * t230 + t232 * t266 - t233 * t265 + t379;
t1 = t358 * ((t299 * t347 + t300 * t348 + t328 * t358) * t371 + ((t302 * t377 + t304 * t375) * t348 + (t301 * t377 + t303 * t375) * t347 + (t329 * t377 + t330 * t375) * t358) * t368) / 0.2e1 + t347 * ((-t300 * t407 + t302 * t334 + t304 * t335) * t348 + (-t299 * t407 + t301 * t334 + t303 * t335) * t347 + (-t328 * t407 + t329 * t334 + t330 * t335) * t358) / 0.2e1 + t348 * ((t300 * t408 + t302 * t336 + t304 * t337) * t348 + (t299 * t408 + t301 * t336 + t303 * t337) * t347 + (t328 * t408 + t329 * t336 + t330 * t337) * t358) / 0.2e1 + m(1) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + m(2) * (t316 ^ 2 + t323 ^ 2 + t324 ^ 2) / 0.2e1 + t283 * ((t214 * t317 - t216 * t284 + t218 * t285) * t266 + (t213 * t317 - t215 * t284 + t217 * t285) * t265 + (t317 * t240 - t284 * t241 + t285 * t242) * t283) / 0.2e1 + m(7) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(6) * (t188 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + m(5) * (t191 ^ 2 + t192 ^ 2 + t193 ^ 2) / 0.2e1 + m(4) * (t194 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + m(3) * (t244 ^ 2 + t260 ^ 2 + t261 ^ 2) / 0.2e1 + t230 * ((t198 * t269 + t200 * t234 + t202 * t235) * t231 + (t197 * t269 + t234 * t199 + t235 * t201) * t230 + (t209 * t269 + t210 * t234 + t211 * t235) * t259) / 0.2e1 + t231 * ((t198 * t271 + t236 * t200 + t237 * t202) * t231 + (t197 * t271 + t199 * t236 + t201 * t237) * t230 + (t209 * t271 + t210 * t236 + t211 * t237) * t259) / 0.2e1 + t259 * ((t198 * t284 + t200 * t267 + t202 * t268) * t231 + (t197 * t284 + t199 * t267 + t201 * t268) * t230 + (t209 * t284 + t210 * t267 + t211 * t268) * t259) / 0.2e1 + t265 * ((t214 * t293 - t216 * t269 + t218 * t270) * t266 + (t213 * t293 - t215 * t269 + t217 * t270) * t265 + (t240 * t293 - t241 * t269 + t242 * t270) * t283) / 0.2e1 + t266 * ((t214 * t295 - t216 * t271 + t218 * t272) * t266 + (t213 * t295 - t215 * t271 + t217 * t272) * t265 + (t240 * t295 - t241 * t271 + t242 * t272) * t283) / 0.2e1 + ((t246 * t273 + t247 * t274 + t277 * t319 + t279 * t294 + t293 * t423) * t321 + (t224 * t273 + t226 * t274 + t250 * t319 + t254 * t294 + t293 * t424) * t308 + (t223 * t273 + t225 * t274 + t249 * t319 + t253 * t294 + t425 * t293) * t307) * t307 / 0.2e1 + ((t246 * t275 + t247 * t276 + t277 * t320 + t279 * t296 + t295 * t423) * t321 + (t224 * t275 + t226 * t276 + t250 * t320 + t254 * t296 + t424 * t295) * t308 + (t223 * t275 + t225 * t276 + t249 * t320 + t253 * t296 + t295 * t425) * t307) * t308 / 0.2e1 + ((t246 * t291 + t247 * t292 + t277 * t333 + t279 * t318 + t423 * t317) * t321 + (t224 * t291 + t226 * t292 + t250 * t333 + t254 * t318 + t317 * t424) * t308 + (t223 * t291 + t225 * t292 + t249 * t333 + t253 * t318 + t317 * t425) * t307) * t321 / 0.2e1 + ((-t351 * t367 + t353 * t370 + Icges(1,4)) * V_base(5) + (-t352 * t367 + t354 * t370 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t351 * t370 + t353 * t367 + Icges(1,2)) * V_base(5) + (t352 * t370 + t354 * t367 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t370 - Icges(2,6) * t367 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t367 + Icges(2,6) * t370 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
