% Calculate kinetic energy for
% S6RPRRPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR9_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR9_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR9_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR9_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR9_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR9_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:26:44
% EndTime: 2019-03-09 05:26:49
% DurationCPUTime: 5.24s
% Computational Cost: add. (5362->447), mult. (12729->633), div. (0->0), fcn. (16211->16), ass. (0->193)
t424 = Icges(5,3) + Icges(6,3);
t367 = cos(pkin(12));
t365 = sin(pkin(12));
t373 = sin(qJ(1));
t402 = t373 * t365;
t368 = cos(pkin(6));
t376 = cos(qJ(1));
t403 = t368 * t376;
t336 = t367 * t403 - t402;
t401 = t373 * t367;
t337 = t365 * t403 + t401;
t372 = sin(qJ(3));
t366 = sin(pkin(6));
t412 = sin(pkin(7));
t393 = t366 * t412;
t413 = cos(pkin(7));
t416 = cos(qJ(3));
t294 = t337 * t416 + (t336 * t413 - t376 * t393) * t372;
t394 = t366 * t413;
t318 = -t336 * t412 - t376 * t394;
t399 = qJ(4) + pkin(13);
t361 = sin(t399);
t392 = cos(t399);
t268 = t294 * t361 - t318 * t392;
t269 = t294 * t392 + t318 * t361;
t371 = sin(qJ(4));
t375 = cos(qJ(4));
t274 = -t294 * t371 + t318 * t375;
t409 = t318 * t371;
t275 = t294 * t375 + t409;
t389 = t416 * t412;
t387 = t366 * t389;
t390 = t413 * t416;
t293 = -t336 * t390 + t337 * t372 + t376 * t387;
t423 = Icges(5,5) * t275 + Icges(6,5) * t269 + Icges(5,6) * t274 - Icges(6,6) * t268 + t293 * t424;
t338 = -t365 * t376 - t368 * t401;
t339 = t367 * t376 - t368 * t402;
t296 = t339 * t416 + (t338 * t413 + t373 * t393) * t372;
t319 = -t338 * t412 + t373 * t394;
t270 = t296 * t361 - t319 * t392;
t271 = t296 * t392 + t319 * t361;
t276 = -t296 * t371 + t319 * t375;
t408 = t319 * t371;
t277 = t296 * t375 + t408;
t295 = -t338 * t390 + t339 * t372 - t373 * t387;
t422 = Icges(5,5) * t277 + Icges(6,5) * t271 + Icges(5,6) * t276 - Icges(6,6) * t270 + t295 * t424;
t317 = t368 * t412 * t372 + (t367 * t372 * t413 + t365 * t416) * t366;
t335 = -t367 * t393 + t368 * t413;
t283 = t317 * t361 - t335 * t392;
t284 = t317 * t392 + t335 * t361;
t290 = -t317 * t371 + t335 * t375;
t407 = t335 * t371;
t291 = t317 * t375 + t407;
t406 = t365 * t366;
t316 = -t366 * t367 * t390 - t368 * t389 + t372 * t406;
t421 = Icges(5,5) * t291 + Icges(6,5) * t284 + Icges(5,6) * t290 - Icges(6,6) * t283 + t316 * t424;
t415 = pkin(4) * t375;
t411 = Icges(2,4) * t373;
t410 = qJ(2) * t368;
t405 = t366 * t373;
t404 = t366 * t376;
t400 = qJD(2) * t366;
t398 = V_base(5) * pkin(8) + V_base(1);
t311 = qJD(3) * t319 + V_base(4);
t310 = qJD(3) * t318 + V_base(5);
t362 = V_base(6) + qJD(1);
t395 = -pkin(8) - t410;
t341 = t373 * pkin(1) - qJ(2) * t404;
t391 = qJD(2) * t368 + V_base(4) * t341 + V_base(3);
t273 = qJD(4) * t295 + t311;
t272 = qJD(4) * t293 + t310;
t325 = qJD(3) * t335 + t362;
t388 = t373 * t400 + V_base(5) * t410 + t398;
t285 = qJD(4) * t316 + t325;
t342 = pkin(1) * t376 + qJ(2) * t405;
t386 = t362 * t342 - t376 * t400 + V_base(2);
t298 = t337 * pkin(2) + pkin(9) * t318;
t322 = pkin(2) * t406 + pkin(9) * t335;
t385 = V_base(5) * t322 + (-t298 - t341) * t362 + t388;
t299 = t339 * pkin(2) + pkin(9) * t319;
t384 = V_base(4) * t298 + (-t299 - t342) * V_base(5) + t391;
t263 = pkin(3) * t294 + pkin(10) * t293;
t282 = pkin(3) * t317 + pkin(10) * t316;
t383 = -t263 * t325 + t310 * t282 + t385;
t264 = pkin(3) * t296 + pkin(10) * t295;
t382 = t311 * t263 - t264 * t310 + t384;
t381 = t362 * t299 + (-t322 + t395) * V_base(4) + t386;
t235 = pkin(4) * t407 + qJ(5) * t316 + t317 * t415;
t380 = qJD(5) * t295 + t272 * t235 + t383;
t208 = pkin(4) * t409 + qJ(5) * t293 + t294 * t415;
t379 = qJD(5) * t316 + t273 * t208 + t382;
t378 = t325 * t264 - t311 * t282 + t381;
t209 = pkin(4) * t408 + qJ(5) * t295 + t296 * t415;
t377 = qJD(5) * t293 + t285 * t209 + t378;
t374 = cos(qJ(6));
t370 = sin(qJ(6));
t363 = Icges(2,4) * t376;
t355 = rSges(2,1) * t376 - t373 * rSges(2,2);
t354 = t373 * rSges(2,1) + rSges(2,2) * t376;
t353 = Icges(2,1) * t376 - t411;
t352 = Icges(2,1) * t373 + t363;
t351 = -Icges(2,2) * t373 + t363;
t350 = Icges(2,2) * t376 + t411;
t349 = Icges(2,5) * t376 - Icges(2,6) * t373;
t348 = Icges(2,5) * t373 + Icges(2,6) * t376;
t347 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t346 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t345 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t332 = rSges(3,3) * t368 + (rSges(3,1) * t365 + rSges(3,2) * t367) * t366;
t331 = Icges(3,5) * t368 + (Icges(3,1) * t365 + Icges(3,4) * t367) * t366;
t330 = Icges(3,6) * t368 + (Icges(3,4) * t365 + Icges(3,2) * t367) * t366;
t329 = Icges(3,3) * t368 + (Icges(3,5) * t365 + Icges(3,6) * t367) * t366;
t324 = V_base(5) * rSges(2,3) - t354 * t362 + t398;
t323 = t355 * t362 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t321 = t354 * V_base(4) - t355 * V_base(5) + V_base(3);
t307 = rSges(3,1) * t339 + rSges(3,2) * t338 + rSges(3,3) * t405;
t306 = t337 * rSges(3,1) + t336 * rSges(3,2) - rSges(3,3) * t404;
t305 = Icges(3,1) * t339 + Icges(3,4) * t338 + Icges(3,5) * t405;
t304 = Icges(3,1) * t337 + Icges(3,4) * t336 - Icges(3,5) * t404;
t303 = Icges(3,4) * t339 + Icges(3,2) * t338 + Icges(3,6) * t405;
t302 = Icges(3,4) * t337 + Icges(3,2) * t336 - Icges(3,6) * t404;
t301 = Icges(3,5) * t339 + Icges(3,6) * t338 + Icges(3,3) * t405;
t300 = Icges(3,5) * t337 + Icges(3,6) * t336 - Icges(3,3) * t404;
t281 = rSges(4,1) * t317 - rSges(4,2) * t316 + rSges(4,3) * t335;
t280 = Icges(4,1) * t317 - Icges(4,4) * t316 + Icges(4,5) * t335;
t279 = Icges(4,4) * t317 - Icges(4,2) * t316 + Icges(4,6) * t335;
t278 = Icges(4,5) * t317 - Icges(4,6) * t316 + Icges(4,3) * t335;
t267 = t284 * t374 + t316 * t370;
t266 = -t284 * t370 + t316 * t374;
t262 = t332 * V_base(5) + (-t306 - t341) * t362 + t388;
t261 = t362 * t307 + (-t332 + t395) * V_base(4) + t386;
t260 = qJD(6) * t283 + t285;
t259 = pkin(5) * t284 + pkin(11) * t283;
t258 = t306 * V_base(4) + (-t307 - t342) * V_base(5) + t391;
t256 = rSges(4,1) * t296 - rSges(4,2) * t295 + rSges(4,3) * t319;
t255 = rSges(4,1) * t294 - rSges(4,2) * t293 + rSges(4,3) * t318;
t254 = Icges(4,1) * t296 - Icges(4,4) * t295 + Icges(4,5) * t319;
t253 = Icges(4,1) * t294 - Icges(4,4) * t293 + Icges(4,5) * t318;
t252 = Icges(4,4) * t296 - Icges(4,2) * t295 + Icges(4,6) * t319;
t251 = Icges(4,4) * t294 - Icges(4,2) * t293 + Icges(4,6) * t318;
t250 = Icges(4,5) * t296 - Icges(4,6) * t295 + Icges(4,3) * t319;
t249 = Icges(4,5) * t294 - Icges(4,6) * t293 + Icges(4,3) * t318;
t248 = rSges(5,1) * t291 + rSges(5,2) * t290 + rSges(5,3) * t316;
t246 = Icges(5,1) * t291 + Icges(5,4) * t290 + Icges(5,5) * t316;
t245 = Icges(5,4) * t291 + Icges(5,2) * t290 + Icges(5,6) * t316;
t243 = rSges(6,1) * t284 - rSges(6,2) * t283 + rSges(6,3) * t316;
t242 = t271 * t374 + t295 * t370;
t241 = -t271 * t370 + t295 * t374;
t240 = t269 * t374 + t293 * t370;
t239 = -t269 * t370 + t293 * t374;
t238 = Icges(6,1) * t284 - Icges(6,4) * t283 + Icges(6,5) * t316;
t237 = Icges(6,4) * t284 - Icges(6,2) * t283 + Icges(6,6) * t316;
t234 = qJD(6) * t270 + t273;
t233 = qJD(6) * t268 + t272;
t232 = pkin(5) * t271 + pkin(11) * t270;
t231 = pkin(5) * t269 + pkin(11) * t268;
t230 = rSges(5,1) * t277 + rSges(5,2) * t276 + rSges(5,3) * t295;
t229 = rSges(5,1) * t275 + rSges(5,2) * t274 + rSges(5,3) * t293;
t228 = Icges(5,1) * t277 + Icges(5,4) * t276 + Icges(5,5) * t295;
t227 = Icges(5,1) * t275 + Icges(5,4) * t274 + Icges(5,5) * t293;
t226 = Icges(5,4) * t277 + Icges(5,2) * t276 + Icges(5,6) * t295;
t225 = Icges(5,4) * t275 + Icges(5,2) * t274 + Icges(5,6) * t293;
t221 = rSges(6,1) * t271 - rSges(6,2) * t270 + rSges(6,3) * t295;
t220 = rSges(6,1) * t269 - rSges(6,2) * t268 + rSges(6,3) * t293;
t219 = Icges(6,1) * t271 - Icges(6,4) * t270 + Icges(6,5) * t295;
t218 = Icges(6,1) * t269 - Icges(6,4) * t268 + Icges(6,5) * t293;
t217 = Icges(6,4) * t271 - Icges(6,2) * t270 + Icges(6,6) * t295;
t216 = Icges(6,4) * t269 - Icges(6,2) * t268 + Icges(6,6) * t293;
t213 = rSges(7,1) * t267 + rSges(7,2) * t266 + rSges(7,3) * t283;
t212 = Icges(7,1) * t267 + Icges(7,4) * t266 + Icges(7,5) * t283;
t211 = Icges(7,4) * t267 + Icges(7,2) * t266 + Icges(7,6) * t283;
t210 = Icges(7,5) * t267 + Icges(7,6) * t266 + Icges(7,3) * t283;
t205 = -t255 * t325 + t281 * t310 + t385;
t204 = t325 * t256 - t311 * t281 + t381;
t203 = rSges(7,1) * t242 + rSges(7,2) * t241 + rSges(7,3) * t270;
t202 = rSges(7,1) * t240 + rSges(7,2) * t239 + rSges(7,3) * t268;
t201 = Icges(7,1) * t242 + Icges(7,4) * t241 + Icges(7,5) * t270;
t200 = Icges(7,1) * t240 + Icges(7,4) * t239 + Icges(7,5) * t268;
t199 = Icges(7,4) * t242 + Icges(7,2) * t241 + Icges(7,6) * t270;
t198 = Icges(7,4) * t240 + Icges(7,2) * t239 + Icges(7,6) * t268;
t197 = Icges(7,5) * t242 + Icges(7,6) * t241 + Icges(7,3) * t270;
t196 = Icges(7,5) * t240 + Icges(7,6) * t239 + Icges(7,3) * t268;
t195 = t255 * t311 - t256 * t310 + t384;
t194 = -t229 * t285 + t248 * t272 + t383;
t193 = t285 * t230 - t273 * t248 + t378;
t192 = t229 * t273 - t230 * t272 + t382;
t191 = t243 * t272 + (-t208 - t220) * t285 + t380;
t190 = t285 * t221 + (-t235 - t243) * t273 + t377;
t189 = t220 * t273 + (-t209 - t221) * t272 + t379;
t188 = -t202 * t260 + t213 * t233 + t259 * t272 + t380 + (-t208 - t231) * t285;
t187 = t260 * t203 - t234 * t213 + t285 * t232 + (-t235 - t259) * t273 + t377;
t186 = t202 * t234 - t203 * t233 + t231 * t273 + (-t209 - t232) * t272 + t379;
t1 = m(1) * (t345 ^ 2 + t346 ^ 2 + t347 ^ 2) / 0.2e1 + t325 * ((t250 * t335 - t252 * t316 + t254 * t317) * t311 + (t249 * t335 - t251 * t316 + t253 * t317) * t310 + (t278 * t335 - t279 * t316 + t280 * t317) * t325) / 0.2e1 + m(2) * (t321 ^ 2 + t323 ^ 2 + t324 ^ 2) / 0.2e1 + t310 * ((t250 * t318 - t252 * t293 + t254 * t294) * t311 + (t249 * t318 - t251 * t293 + t253 * t294) * t310 + (t278 * t318 - t279 * t293 + t280 * t294) * t325) / 0.2e1 + t311 * ((t250 * t319 - t252 * t295 + t254 * t296) * t311 + (t249 * t319 - t251 * t295 + t253 * t296) * t310 + (t278 * t319 - t279 * t295 + t280 * t296) * t325) / 0.2e1 + t260 * ((t197 * t283 + t199 * t266 + t201 * t267) * t234 + (t196 * t283 + t198 * t266 + t200 * t267) * t233 + (t210 * t283 + t211 * t266 + t212 * t267) * t260) / 0.2e1 + m(7) * (t186 ^ 2 + t187 ^ 2 + t188 ^ 2) / 0.2e1 + m(6) * (t189 ^ 2 + t190 ^ 2 + t191 ^ 2) / 0.2e1 + m(5) * (t192 ^ 2 + t193 ^ 2 + t194 ^ 2) / 0.2e1 + m(4) * (t195 ^ 2 + t204 ^ 2 + t205 ^ 2) / 0.2e1 + m(3) * (t258 ^ 2 + t261 ^ 2 + t262 ^ 2) / 0.2e1 + t233 * ((t197 * t268 + t199 * t239 + t201 * t240) * t234 + (t268 * t196 + t239 * t198 + t240 * t200) * t233 + (t210 * t268 + t211 * t239 + t212 * t240) * t260) / 0.2e1 + t234 * ((t270 * t197 + t241 * t199 + t242 * t201) * t234 + (t196 * t270 + t198 * t241 + t200 * t242) * t233 + (t210 * t270 + t211 * t241 + t212 * t242) * t260) / 0.2e1 + ((-t237 * t268 + t238 * t269 + t245 * t274 + t246 * t275 + t293 * t421) * t285 + (-t217 * t268 + t219 * t269 + t226 * t274 + t228 * t275 + t293 * t422) * t273 + (-t216 * t268 + t218 * t269 + t225 * t274 + t227 * t275 + t423 * t293) * t272) * t272 / 0.2e1 + ((-t237 * t270 + t238 * t271 + t245 * t276 + t246 * t277 + t295 * t421) * t285 + (-t217 * t270 + t219 * t271 + t226 * t276 + t228 * t277 + t422 * t295) * t273 + (-t216 * t270 + t218 * t271 + t225 * t276 + t227 * t277 + t295 * t423) * t272) * t273 / 0.2e1 + ((-t237 * t283 + t238 * t284 + t245 * t290 + t246 * t291 + t421 * t316) * t285 + (-t217 * t283 + t219 * t284 + t226 * t290 + t228 * t291 + t316 * t422) * t273 + (-t216 * t283 + t218 * t284 + t225 * t290 + t227 * t291 + t316 * t423) * t272) * t285 / 0.2e1 + ((t300 * V_base(5) + t301 * V_base(4) + t329 * t362) * t368 + ((t303 * t367 + t305 * t365) * V_base(4) + (t302 * t367 + t304 * t365) * V_base(5) + (t330 * t367 + t331 * t365) * t362) * t366 + Icges(2,3) * t362 + t348 * V_base(5) + t349 * V_base(4)) * t362 / 0.2e1 + ((t329 * t405 + t330 * t338 + t331 * t339 + t349) * t362 + (t300 * t405 + t302 * t338 + t304 * t339 - t373 * t350 + t352 * t376 + Icges(1,4)) * V_base(5) + (t301 * t405 + t303 * t338 + t305 * t339 - t373 * t351 + t353 * t376 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t329 * t404 + t336 * t330 + t337 * t331 + t348) * t362 + (-t300 * t404 + t336 * t302 + t337 * t304 + t350 * t376 + t373 * t352 + Icges(1,2)) * V_base(5) + (-t301 * t404 + t336 * t303 + t337 * t305 + t351 * t376 + t373 * t353 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
