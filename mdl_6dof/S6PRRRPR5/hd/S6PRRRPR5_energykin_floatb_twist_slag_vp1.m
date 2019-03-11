% Calculate kinetic energy for
% S6PRRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRPR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:22:54
% EndTime: 2019-03-08 23:22:59
% DurationCPUTime: 5.23s
% Computational Cost: add. (5437->442), mult. (12954->643), div. (0->0), fcn. (16436->16), ass. (0->192)
t426 = Icges(5,3) + Icges(6,3);
t366 = sin(pkin(12));
t368 = cos(pkin(12));
t374 = sin(qJ(2));
t369 = cos(pkin(6));
t377 = cos(qJ(2));
t403 = t369 * t377;
t334 = -t366 * t374 + t368 * t403;
t404 = t369 * t374;
t335 = t366 * t377 + t368 * t404;
t373 = sin(qJ(3));
t367 = sin(pkin(6));
t412 = sin(pkin(7));
t394 = t367 * t412;
t413 = cos(pkin(7));
t417 = cos(qJ(3));
t292 = t335 * t417 + (t334 * t413 - t368 * t394) * t373;
t395 = t367 * t413;
t319 = -t334 * t412 - t368 * t395;
t401 = qJ(4) + pkin(13);
t363 = sin(t401);
t393 = cos(t401);
t269 = t292 * t363 - t319 * t393;
t270 = t292 * t393 + t319 * t363;
t372 = sin(qJ(4));
t376 = cos(qJ(4));
t273 = -t292 * t372 + t319 * t376;
t410 = t319 * t372;
t274 = t292 * t376 + t410;
t391 = t417 * t412;
t390 = t367 * t391;
t392 = t413 * t417;
t291 = -t334 * t392 + t335 * t373 + t368 * t390;
t424 = Icges(5,5) * t274 + Icges(6,5) * t270 + Icges(5,6) * t273 - Icges(6,6) * t269 + t291 * t426;
t336 = -t366 * t403 - t368 * t374;
t337 = -t366 * t404 + t368 * t377;
t294 = t337 * t417 + (t336 * t413 + t366 * t394) * t373;
t320 = -t336 * t412 + t366 * t395;
t271 = t294 * t363 - t320 * t393;
t272 = t294 * t393 + t320 * t363;
t275 = -t294 * t372 + t320 * t376;
t409 = t320 * t372;
t276 = t294 * t376 + t409;
t293 = -t336 * t392 + t337 * t373 - t366 * t390;
t423 = Icges(5,5) * t276 + Icges(6,5) * t272 + Icges(5,6) * t275 - Icges(6,6) * t271 + t293 * t426;
t318 = t369 * t412 * t373 + (t373 * t377 * t413 + t374 * t417) * t367;
t333 = t369 * t413 - t377 * t394;
t284 = t318 * t363 - t333 * t393;
t285 = t318 * t393 + t333 * t363;
t295 = -t318 * t372 + t333 * t376;
t408 = t333 * t372;
t296 = t318 * t376 + t408;
t405 = t367 * t374;
t317 = -t367 * t377 * t392 - t369 * t391 + t373 * t405;
t422 = Icges(5,5) * t296 + Icges(6,5) * t285 + Icges(5,6) * t295 - Icges(6,6) * t284 + t317 * t426;
t416 = pkin(8) * t369;
t415 = pkin(4) * t376;
t411 = Icges(2,4) * t366;
t407 = t366 * t367;
t406 = t367 * t368;
t402 = qJD(2) * t367;
t400 = V_base(5) * qJ(1) + V_base(1);
t396 = qJD(1) + V_base(3);
t348 = t366 * t402 + V_base(4);
t358 = qJD(2) * t369 + V_base(6);
t308 = qJD(3) * t320 + t348;
t321 = qJD(3) * t333 + t358;
t266 = qJD(4) * t293 + t308;
t283 = qJD(4) * t317 + t321;
t347 = -t368 * t402 + V_base(5);
t307 = qJD(3) * t319 + t347;
t340 = pkin(1) * t366 - pkin(8) * t406;
t389 = -t340 * V_base(6) + V_base(5) * t416 + t400;
t341 = pkin(1) * t368 + pkin(8) * t407;
t388 = V_base(4) * t340 - t341 * V_base(5) + t396;
t265 = qJD(4) * t291 + t307;
t387 = V_base(6) * t341 + V_base(2) + (-qJ(1) - t416) * V_base(4);
t297 = t335 * pkin(2) + pkin(9) * t319;
t322 = pkin(2) * t405 + pkin(9) * t333;
t386 = -t297 * t358 + t347 * t322 + t389;
t298 = t337 * pkin(2) + pkin(9) * t320;
t385 = t348 * t297 - t298 * t347 + t388;
t384 = t358 * t298 - t322 * t348 + t387;
t262 = pkin(3) * t292 + pkin(10) * t291;
t281 = pkin(3) * t318 + pkin(10) * t317;
t383 = -t262 * t321 + t307 * t281 + t386;
t263 = pkin(3) * t294 + pkin(10) * t293;
t382 = t308 * t262 - t263 * t307 + t385;
t381 = t321 * t263 - t281 * t308 + t384;
t239 = pkin(4) * t408 + qJ(5) * t317 + t318 * t415;
t380 = qJD(5) * t293 + t265 * t239 + t383;
t207 = pkin(4) * t410 + qJ(5) * t291 + t292 * t415;
t379 = qJD(5) * t317 + t266 * t207 + t382;
t208 = pkin(4) * t409 + qJ(5) * t293 + t294 * t415;
t378 = qJD(5) * t291 + t283 * t208 + t381;
t375 = cos(qJ(6));
t371 = sin(qJ(6));
t364 = Icges(2,4) * t368;
t356 = rSges(2,1) * t368 - rSges(2,2) * t366;
t355 = rSges(2,1) * t366 + rSges(2,2) * t368;
t354 = Icges(2,1) * t368 - t411;
t353 = Icges(2,1) * t366 + t364;
t352 = -Icges(2,2) * t366 + t364;
t351 = Icges(2,2) * t368 + t411;
t346 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t345 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t344 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t331 = t369 * rSges(3,3) + (rSges(3,1) * t374 + rSges(3,2) * t377) * t367;
t330 = Icges(3,5) * t369 + (Icges(3,1) * t374 + Icges(3,4) * t377) * t367;
t329 = Icges(3,6) * t369 + (Icges(3,4) * t374 + Icges(3,2) * t377) * t367;
t328 = Icges(3,3) * t369 + (Icges(3,5) * t374 + Icges(3,6) * t377) * t367;
t324 = V_base(5) * rSges(2,3) - t355 * V_base(6) + t400;
t323 = t356 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t316 = t355 * V_base(4) - t356 * V_base(5) + t396;
t306 = rSges(3,1) * t337 + rSges(3,2) * t336 + rSges(3,3) * t407;
t305 = rSges(3,1) * t335 + rSges(3,2) * t334 - rSges(3,3) * t406;
t304 = Icges(3,1) * t337 + Icges(3,4) * t336 + Icges(3,5) * t407;
t303 = Icges(3,1) * t335 + Icges(3,4) * t334 - Icges(3,5) * t406;
t302 = Icges(3,4) * t337 + Icges(3,2) * t336 + Icges(3,6) * t407;
t301 = Icges(3,4) * t335 + Icges(3,2) * t334 - Icges(3,6) * t406;
t300 = Icges(3,5) * t337 + Icges(3,6) * t336 + Icges(3,3) * t407;
t299 = Icges(3,5) * t335 + Icges(3,6) * t334 - Icges(3,3) * t406;
t280 = rSges(4,1) * t318 - rSges(4,2) * t317 + rSges(4,3) * t333;
t279 = Icges(4,1) * t318 - Icges(4,4) * t317 + Icges(4,5) * t333;
t278 = Icges(4,4) * t318 - Icges(4,2) * t317 + Icges(4,6) * t333;
t277 = Icges(4,5) * t318 - Icges(4,6) * t317 + Icges(4,3) * t333;
t268 = t285 * t375 + t317 * t371;
t267 = -t285 * t371 + t317 * t375;
t261 = -t305 * t358 + t331 * t347 + t389;
t260 = t306 * t358 - t331 * t348 + t387;
t259 = qJD(6) * t284 + t283;
t258 = pkin(5) * t285 + pkin(11) * t284;
t256 = rSges(5,1) * t296 + rSges(5,2) * t295 + rSges(5,3) * t317;
t255 = rSges(4,1) * t294 - rSges(4,2) * t293 + rSges(4,3) * t320;
t254 = rSges(4,1) * t292 - rSges(4,2) * t291 + rSges(4,3) * t319;
t253 = Icges(5,1) * t296 + Icges(5,4) * t295 + Icges(5,5) * t317;
t252 = Icges(5,4) * t296 + Icges(5,2) * t295 + Icges(5,6) * t317;
t250 = Icges(4,1) * t294 - Icges(4,4) * t293 + Icges(4,5) * t320;
t249 = Icges(4,1) * t292 - Icges(4,4) * t291 + Icges(4,5) * t319;
t248 = Icges(4,4) * t294 - Icges(4,2) * t293 + Icges(4,6) * t320;
t247 = Icges(4,4) * t292 - Icges(4,2) * t291 + Icges(4,6) * t319;
t246 = Icges(4,5) * t294 - Icges(4,6) * t293 + Icges(4,3) * t320;
t245 = Icges(4,5) * t292 - Icges(4,6) * t291 + Icges(4,3) * t319;
t244 = t305 * t348 - t306 * t347 + t388;
t243 = rSges(6,1) * t285 - rSges(6,2) * t284 + rSges(6,3) * t317;
t242 = Icges(6,1) * t285 - Icges(6,4) * t284 + Icges(6,5) * t317;
t241 = Icges(6,4) * t285 - Icges(6,2) * t284 + Icges(6,6) * t317;
t237 = t272 * t375 + t293 * t371;
t236 = -t272 * t371 + t293 * t375;
t235 = t270 * t375 + t291 * t371;
t234 = -t270 * t371 + t291 * t375;
t233 = pkin(5) * t272 + pkin(11) * t271;
t232 = pkin(5) * t270 + pkin(11) * t269;
t231 = qJD(6) * t271 + t266;
t230 = qJD(6) * t269 + t265;
t229 = rSges(5,1) * t276 + rSges(5,2) * t275 + rSges(5,3) * t293;
t228 = rSges(5,1) * t274 + rSges(5,2) * t273 + rSges(5,3) * t291;
t227 = Icges(5,1) * t276 + Icges(5,4) * t275 + Icges(5,5) * t293;
t226 = Icges(5,1) * t274 + Icges(5,4) * t273 + Icges(5,5) * t291;
t225 = Icges(5,4) * t276 + Icges(5,2) * t275 + Icges(5,6) * t293;
t224 = Icges(5,4) * t274 + Icges(5,2) * t273 + Icges(5,6) * t291;
t221 = rSges(6,1) * t272 - rSges(6,2) * t271 + rSges(6,3) * t293;
t220 = rSges(6,1) * t270 - rSges(6,2) * t269 + rSges(6,3) * t291;
t219 = Icges(6,1) * t272 - Icges(6,4) * t271 + Icges(6,5) * t293;
t218 = Icges(6,1) * t270 - Icges(6,4) * t269 + Icges(6,5) * t291;
t217 = Icges(6,4) * t272 - Icges(6,2) * t271 + Icges(6,6) * t293;
t216 = Icges(6,4) * t270 - Icges(6,2) * t269 + Icges(6,6) * t291;
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
t196 = -t254 * t321 + t280 * t307 + t386;
t195 = t255 * t321 - t280 * t308 + t384;
t194 = t254 * t308 - t255 * t307 + t385;
t193 = -t228 * t283 + t256 * t265 + t383;
t192 = t229 * t283 - t256 * t266 + t381;
t191 = t228 * t266 - t229 * t265 + t382;
t190 = t243 * t265 + (-t207 - t220) * t283 + t380;
t189 = t221 * t283 + (-t239 - t243) * t266 + t378;
t188 = t220 * t266 + (-t208 - t221) * t265 + t379;
t187 = -t203 * t259 + t380 + t212 * t230 + t258 * t265 + (-t207 - t232) * t283;
t186 = t204 * t259 - t212 * t231 + t233 * t283 + (-t239 - t258) * t266 + t378;
t185 = t203 * t231 - t204 * t230 + t232 * t266 + (-t208 - t233) * t265 + t379;
t1 = t347 * ((-t300 * t406 + t302 * t334 + t304 * t335) * t348 + (-t299 * t406 + t301 * t334 + t303 * t335) * t347 + (-t328 * t406 + t329 * t334 + t330 * t335) * t358) / 0.2e1 + t348 * ((t300 * t407 + t302 * t336 + t304 * t337) * t348 + (t299 * t407 + t301 * t336 + t303 * t337) * t347 + (t328 * t407 + t329 * t336 + t330 * t337) * t358) / 0.2e1 + t358 * ((t299 * t347 + t300 * t348 + t328 * t358) * t369 + ((t302 * t377 + t304 * t374) * t348 + (t301 * t377 + t303 * t374) * t347 + (t329 * t377 + t330 * t374) * t358) * t367) / 0.2e1 + m(1) * (t344 ^ 2 + t345 ^ 2 + t346 ^ 2) / 0.2e1 + t321 * ((t246 * t333 - t248 * t317 + t250 * t318) * t308 + (t245 * t333 - t247 * t317 + t249 * t318) * t307 + (t277 * t333 - t278 * t317 + t279 * t318) * t321) / 0.2e1 + t307 * ((t246 * t319 - t248 * t291 + t250 * t292) * t308 + (t245 * t319 - t247 * t291 + t249 * t292) * t307 + (t277 * t319 - t278 * t291 + t279 * t292) * t321) / 0.2e1 + t308 * ((t246 * t320 - t248 * t293 + t250 * t294) * t308 + (t245 * t320 - t247 * t293 + t249 * t294) * t307 + (t277 * t320 - t278 * t293 + t279 * t294) * t321) / 0.2e1 + m(2) * (t316 ^ 2 + t323 ^ 2 + t324 ^ 2) / 0.2e1 + m(7) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(6) * (t188 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + m(5) * (t191 ^ 2 + t192 ^ 2 + t193 ^ 2) / 0.2e1 + m(4) * (t194 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + m(3) * (t244 ^ 2 + t260 ^ 2 + t261 ^ 2) / 0.2e1 + t230 * ((t198 * t269 + t200 * t234 + t202 * t235) * t231 + (t197 * t269 + t199 * t234 + t201 * t235) * t230 + (t209 * t269 + t210 * t234 + t211 * t235) * t259) / 0.2e1 + t231 * ((t198 * t271 + t200 * t236 + t202 * t237) * t231 + (t197 * t271 + t199 * t236 + t201 * t237) * t230 + (t209 * t271 + t210 * t236 + t211 * t237) * t259) / 0.2e1 + t259 * ((t198 * t284 + t200 * t267 + t202 * t268) * t231 + (t197 * t284 + t199 * t267 + t201 * t268) * t230 + (t209 * t284 + t210 * t267 + t211 * t268) * t259) / 0.2e1 + ((-t241 * t269 + t242 * t270 + t252 * t273 + t253 * t274 + t291 * t422) * t283 + (-t217 * t269 + t219 * t270 + t225 * t273 + t227 * t274 + t291 * t423) * t266 + (-t216 * t269 + t218 * t270 + t224 * t273 + t226 * t274 + t291 * t424) * t265) * t265 / 0.2e1 + ((-t241 * t271 + t242 * t272 + t252 * t275 + t253 * t276 + t293 * t422) * t283 + (-t217 * t271 + t219 * t272 + t225 * t275 + t227 * t276 + t293 * t423) * t266 + (-t216 * t271 + t218 * t272 + t224 * t275 + t226 * t276 + t293 * t424) * t265) * t266 / 0.2e1 + ((-t241 * t284 + t242 * t285 + t252 * t295 + t253 * t296 + t317 * t422) * t283 + (-t217 * t284 + t219 * t285 + t225 * t295 + t227 * t296 + t317 * t423) * t266 + (-t216 * t284 + t218 * t285 + t224 * t295 + t226 * t296 + t317 * t424) * t265) * t283 / 0.2e1 + ((-t351 * t366 + t353 * t368 + Icges(1,4)) * V_base(5) + (-t352 * t366 + t354 * t368 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t351 * t368 + t353 * t366 + Icges(1,2)) * V_base(5) + (t352 * t368 + t354 * t366 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t368 - Icges(2,6) * t366 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t366 + Icges(2,6) * t368 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
