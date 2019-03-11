% Calculate kinetic energy for
% S6RRRRPR12
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 23:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR12_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR12_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR12_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR12_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRPR12_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR12_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR12_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR12_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 23:29:04
% EndTime: 2019-03-09 23:29:10
% DurationCPUTime: 5.28s
% Computational Cost: add. (5497->442), mult. (12954->645), div. (0->0), fcn. (16436->16), ass. (0->194)
t427 = Icges(5,3) + Icges(6,3);
t370 = cos(pkin(6));
t379 = cos(qJ(2));
t380 = cos(qJ(1));
t404 = t379 * t380;
t375 = sin(qJ(2));
t376 = sin(qJ(1));
t406 = t376 * t375;
t337 = t370 * t404 - t406;
t405 = t376 * t379;
t407 = t375 * t380;
t338 = t370 * t407 + t405;
t374 = sin(qJ(3));
t369 = sin(pkin(6));
t415 = sin(pkin(7));
t397 = t369 * t415;
t416 = cos(pkin(7));
t420 = cos(qJ(3));
t296 = t338 * t420 + (t337 * t416 - t380 * t397) * t374;
t398 = t369 * t416;
t320 = -t337 * t415 - t380 * t398;
t402 = qJ(4) + pkin(13);
t365 = sin(t402);
t396 = cos(t402);
t271 = t296 * t365 - t320 * t396;
t272 = t296 * t396 + t320 * t365;
t373 = sin(qJ(4));
t378 = cos(qJ(4));
t275 = -t296 * t373 + t320 * t378;
t413 = t320 * t373;
t276 = t296 * t378 + t413;
t394 = t420 * t415;
t393 = t369 * t394;
t395 = t416 * t420;
t295 = -t337 * t395 + t338 * t374 + t380 * t393;
t426 = Icges(5,5) * t276 + Icges(6,5) * t272 + Icges(5,6) * t275 - Icges(6,6) * t271 + t295 * t427;
t339 = -t370 * t405 - t407;
t340 = -t370 * t406 + t404;
t298 = t340 * t420 + (t339 * t416 + t376 * t397) * t374;
t321 = -t339 * t415 + t376 * t398;
t273 = t298 * t365 - t321 * t396;
t274 = t298 * t396 + t321 * t365;
t277 = -t298 * t373 + t321 * t378;
t412 = t321 * t373;
t278 = t298 * t378 + t412;
t297 = -t339 * t395 + t340 * t374 - t376 * t393;
t425 = Icges(5,5) * t278 + Icges(6,5) * t274 + Icges(5,6) * t277 - Icges(6,6) * t273 + t297 * t427;
t319 = t370 * t415 * t374 + (t374 * t379 * t416 + t375 * t420) * t369;
t336 = t370 * t416 - t379 * t397;
t286 = t319 * t365 - t336 * t396;
t287 = t319 * t396 + t336 * t365;
t293 = -t319 * t373 + t336 * t378;
t411 = t336 * t373;
t294 = t319 * t378 + t411;
t410 = t369 * t375;
t318 = -t369 * t379 * t395 - t370 * t394 + t374 * t410;
t424 = Icges(5,5) * t294 + Icges(6,5) * t287 + Icges(5,6) * t293 - Icges(6,6) * t286 + t318 * t427;
t419 = pkin(9) * t370;
t418 = pkin(4) * t378;
t414 = Icges(2,4) * t376;
t409 = t369 * t376;
t408 = t369 * t380;
t403 = qJD(2) * t369;
t401 = V_base(5) * pkin(8) + V_base(1);
t350 = t376 * t403 + V_base(4);
t366 = V_base(6) + qJD(1);
t310 = qJD(3) * t321 + t350;
t351 = qJD(2) * t370 + t366;
t268 = qJD(4) * t297 + t310;
t322 = qJD(3) * t336 + t351;
t349 = -t380 * t403 + V_base(5);
t342 = t376 * pkin(1) - pkin(9) * t408;
t392 = -t342 * t366 + V_base(5) * t419 + t401;
t284 = qJD(4) * t318 + t322;
t343 = pkin(1) * t380 + pkin(9) * t409;
t391 = V_base(4) * t342 - t343 * V_base(5) + V_base(3);
t309 = qJD(3) * t320 + t349;
t267 = qJD(4) * t295 + t309;
t390 = t366 * t343 + V_base(2) + (-pkin(8) - t419) * V_base(4);
t299 = t338 * pkin(2) + pkin(10) * t320;
t326 = pkin(2) * t410 + pkin(10) * t336;
t389 = -t299 * t351 + t349 * t326 + t392;
t300 = t340 * pkin(2) + pkin(10) * t321;
t388 = t350 * t299 - t300 * t349 + t391;
t387 = t351 * t300 - t326 * t350 + t390;
t264 = pkin(3) * t296 + pkin(11) * t295;
t283 = pkin(3) * t319 + pkin(11) * t318;
t386 = -t264 * t322 + t309 * t283 + t389;
t265 = pkin(3) * t298 + pkin(11) * t297;
t385 = t310 * t264 - t265 * t309 + t388;
t236 = pkin(4) * t411 + qJ(5) * t318 + t319 * t418;
t384 = qJD(5) * t297 + t267 * t236 + t386;
t209 = pkin(4) * t413 + qJ(5) * t295 + t296 * t418;
t383 = qJD(5) * t318 + t268 * t209 + t385;
t382 = t322 * t265 - t283 * t310 + t387;
t210 = pkin(4) * t412 + qJ(5) * t297 + t298 * t418;
t381 = qJD(5) * t295 + t284 * t210 + t382;
t377 = cos(qJ(6));
t372 = sin(qJ(6));
t367 = Icges(2,4) * t380;
t359 = rSges(2,1) * t380 - t376 * rSges(2,2);
t358 = t376 * rSges(2,1) + rSges(2,2) * t380;
t357 = Icges(2,1) * t380 - t414;
t356 = Icges(2,1) * t376 + t367;
t355 = -Icges(2,2) * t376 + t367;
t354 = Icges(2,2) * t380 + t414;
t348 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t347 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t346 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t333 = rSges(3,3) * t370 + (rSges(3,1) * t375 + rSges(3,2) * t379) * t369;
t332 = Icges(3,5) * t370 + (Icges(3,1) * t375 + Icges(3,4) * t379) * t369;
t331 = Icges(3,6) * t370 + (Icges(3,4) * t375 + Icges(3,2) * t379) * t369;
t330 = Icges(3,3) * t370 + (Icges(3,5) * t375 + Icges(3,6) * t379) * t369;
t325 = V_base(5) * rSges(2,3) - t358 * t366 + t401;
t324 = t359 * t366 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t323 = t358 * V_base(4) - t359 * V_base(5) + V_base(3);
t308 = rSges(3,1) * t340 + rSges(3,2) * t339 + rSges(3,3) * t409;
t307 = t338 * rSges(3,1) + t337 * rSges(3,2) - rSges(3,3) * t408;
t306 = Icges(3,1) * t340 + Icges(3,4) * t339 + Icges(3,5) * t409;
t305 = Icges(3,1) * t338 + Icges(3,4) * t337 - Icges(3,5) * t408;
t304 = Icges(3,4) * t340 + Icges(3,2) * t339 + Icges(3,6) * t409;
t303 = Icges(3,4) * t338 + Icges(3,2) * t337 - Icges(3,6) * t408;
t302 = Icges(3,5) * t340 + Icges(3,6) * t339 + Icges(3,3) * t409;
t301 = Icges(3,5) * t338 + Icges(3,6) * t337 - Icges(3,3) * t408;
t282 = rSges(4,1) * t319 - rSges(4,2) * t318 + rSges(4,3) * t336;
t281 = Icges(4,1) * t319 - Icges(4,4) * t318 + Icges(4,5) * t336;
t280 = Icges(4,4) * t319 - Icges(4,2) * t318 + Icges(4,6) * t336;
t279 = Icges(4,5) * t319 - Icges(4,6) * t318 + Icges(4,3) * t336;
t270 = t287 * t377 + t318 * t372;
t269 = -t287 * t372 + t318 * t377;
t263 = -t307 * t351 + t333 * t349 + t392;
t262 = t308 * t351 - t333 * t350 + t390;
t261 = pkin(5) * t287 + pkin(12) * t286;
t260 = qJD(6) * t286 + t284;
t258 = rSges(4,1) * t298 - rSges(4,2) * t297 + rSges(4,3) * t321;
t257 = rSges(4,1) * t296 - rSges(4,2) * t295 + rSges(4,3) * t320;
t256 = t307 * t350 - t308 * t349 + t391;
t255 = Icges(4,1) * t298 - Icges(4,4) * t297 + Icges(4,5) * t321;
t254 = Icges(4,1) * t296 - Icges(4,4) * t295 + Icges(4,5) * t320;
t253 = Icges(4,4) * t298 - Icges(4,2) * t297 + Icges(4,6) * t321;
t252 = Icges(4,4) * t296 - Icges(4,2) * t295 + Icges(4,6) * t320;
t251 = Icges(4,5) * t298 - Icges(4,6) * t297 + Icges(4,3) * t321;
t250 = Icges(4,5) * t296 - Icges(4,6) * t295 + Icges(4,3) * t320;
t249 = rSges(5,1) * t294 + rSges(5,2) * t293 + rSges(5,3) * t318;
t248 = Icges(5,1) * t294 + Icges(5,4) * t293 + Icges(5,5) * t318;
t247 = Icges(5,4) * t294 + Icges(5,2) * t293 + Icges(5,6) * t318;
t244 = t274 * t377 + t297 * t372;
t243 = -t274 * t372 + t297 * t377;
t242 = t272 * t377 + t295 * t372;
t241 = -t272 * t372 + t295 * t377;
t240 = rSges(6,1) * t287 - rSges(6,2) * t286 + rSges(6,3) * t318;
t239 = Icges(6,1) * t287 - Icges(6,4) * t286 + Icges(6,5) * t318;
t238 = Icges(6,4) * t287 - Icges(6,2) * t286 + Icges(6,6) * t318;
t235 = pkin(5) * t274 + pkin(12) * t273;
t234 = pkin(5) * t272 + pkin(12) * t271;
t233 = qJD(6) * t273 + t268;
t232 = qJD(6) * t271 + t267;
t231 = rSges(5,1) * t278 + rSges(5,2) * t277 + rSges(5,3) * t297;
t230 = rSges(5,1) * t276 + rSges(5,2) * t275 + rSges(5,3) * t295;
t229 = Icges(5,1) * t278 + Icges(5,4) * t277 + Icges(5,5) * t297;
t228 = Icges(5,1) * t276 + Icges(5,4) * t275 + Icges(5,5) * t295;
t227 = Icges(5,4) * t278 + Icges(5,2) * t277 + Icges(5,6) * t297;
t226 = Icges(5,4) * t276 + Icges(5,2) * t275 + Icges(5,6) * t295;
t223 = rSges(6,1) * t274 - rSges(6,2) * t273 + rSges(6,3) * t297;
t222 = rSges(6,1) * t272 - rSges(6,2) * t271 + rSges(6,3) * t295;
t221 = Icges(6,1) * t274 - Icges(6,4) * t273 + Icges(6,5) * t297;
t220 = Icges(6,1) * t272 - Icges(6,4) * t271 + Icges(6,5) * t295;
t219 = Icges(6,4) * t274 - Icges(6,2) * t273 + Icges(6,6) * t297;
t218 = Icges(6,4) * t272 - Icges(6,2) * t271 + Icges(6,6) * t295;
t214 = rSges(7,1) * t270 + rSges(7,2) * t269 + rSges(7,3) * t286;
t213 = Icges(7,1) * t270 + Icges(7,4) * t269 + Icges(7,5) * t286;
t212 = Icges(7,4) * t270 + Icges(7,2) * t269 + Icges(7,6) * t286;
t211 = Icges(7,5) * t270 + Icges(7,6) * t269 + Icges(7,3) * t286;
t206 = rSges(7,1) * t244 + rSges(7,2) * t243 + rSges(7,3) * t273;
t205 = rSges(7,1) * t242 + rSges(7,2) * t241 + rSges(7,3) * t271;
t204 = Icges(7,1) * t244 + Icges(7,4) * t243 + Icges(7,5) * t273;
t203 = Icges(7,1) * t242 + Icges(7,4) * t241 + Icges(7,5) * t271;
t202 = Icges(7,4) * t244 + Icges(7,2) * t243 + Icges(7,6) * t273;
t201 = Icges(7,4) * t242 + Icges(7,2) * t241 + Icges(7,6) * t271;
t200 = Icges(7,5) * t244 + Icges(7,6) * t243 + Icges(7,3) * t273;
t199 = Icges(7,5) * t242 + Icges(7,6) * t241 + Icges(7,3) * t271;
t198 = -t257 * t322 + t282 * t309 + t389;
t197 = t258 * t322 - t282 * t310 + t387;
t196 = t257 * t310 - t258 * t309 + t388;
t195 = -t230 * t284 + t249 * t267 + t386;
t194 = t231 * t284 - t249 * t268 + t382;
t193 = t230 * t268 - t231 * t267 + t385;
t192 = t240 * t267 + (-t209 - t222) * t284 + t384;
t191 = t223 * t284 + (-t236 - t240) * t268 + t381;
t190 = t222 * t268 + (-t210 - t223) * t267 + t383;
t189 = (-t209 - t234) * t284 - t205 * t260 + t214 * t232 + t261 * t267 + t384;
t188 = t206 * t260 - t214 * t233 + t235 * t284 + (-t236 - t261) * t268 + t381;
t187 = t205 * t233 - t206 * t232 + t234 * t268 + (-t210 - t235) * t267 + t383;
t1 = t351 * ((t301 * t349 + t302 * t350 + t330 * t351) * t370 + ((t304 * t379 + t306 * t375) * t350 + (t303 * t379 + t305 * t375) * t349 + (t331 * t379 + t332 * t375) * t351) * t369) / 0.2e1 + t350 * ((t302 * t409 + t304 * t339 + t306 * t340) * t350 + (t301 * t409 + t303 * t339 + t305 * t340) * t349 + (t330 * t409 + t331 * t339 + t332 * t340) * t351) / 0.2e1 + t349 * ((-t302 * t408 + t337 * t304 + t338 * t306) * t350 + (-t301 * t408 + t337 * t303 + t338 * t305) * t349 + (-t330 * t408 + t337 * t331 + t338 * t332) * t351) / 0.2e1 + m(1) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + t322 * ((t251 * t336 - t253 * t318 + t255 * t319) * t310 + (t250 * t336 - t252 * t318 + t254 * t319) * t309 + (t279 * t336 - t280 * t318 + t281 * t319) * t322) / 0.2e1 + t309 * ((t251 * t320 - t253 * t295 + t255 * t296) * t310 + (t250 * t320 - t252 * t295 + t254 * t296) * t309 + (t279 * t320 - t280 * t295 + t281 * t296) * t322) / 0.2e1 + t310 * ((t251 * t321 - t253 * t297 + t255 * t298) * t310 + (t250 * t321 - t252 * t297 + t254 * t298) * t309 + (t279 * t321 - t280 * t297 + t281 * t298) * t322) / 0.2e1 + m(2) * (t323 ^ 2 + t324 ^ 2 + t325 ^ 2) / 0.2e1 + t260 * ((t200 * t286 + t202 * t269 + t204 * t270) * t233 + (t199 * t286 + t201 * t269 + t203 * t270) * t232 + (t286 * t211 + t269 * t212 + t270 * t213) * t260) / 0.2e1 + t233 * ((t273 * t200 + t243 * t202 + t244 * t204) * t233 + (t199 * t273 + t201 * t243 + t203 * t244) * t232 + (t211 * t273 + t212 * t243 + t213 * t244) * t260) / 0.2e1 + t232 * ((t200 * t271 + t202 * t241 + t204 * t242) * t233 + (t271 * t199 + t241 * t201 + t242 * t203) * t232 + (t211 * t271 + t212 * t241 + t213 * t242) * t260) / 0.2e1 + m(3) * (t256 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + m(4) * (t196 ^ 2 + t197 ^ 2 + t198 ^ 2) / 0.2e1 + m(5) * (t193 ^ 2 + t194 ^ 2 + t195 ^ 2) / 0.2e1 + m(6) * (t190 ^ 2 + t191 ^ 2 + t192 ^ 2) / 0.2e1 + m(7) * (t187 ^ 2 + t188 ^ 2 + t189 ^ 2) / 0.2e1 + ((-t238 * t271 + t239 * t272 + t247 * t275 + t248 * t276 + t295 * t424) * t284 + (-t219 * t271 + t221 * t272 + t227 * t275 + t229 * t276 + t295 * t425) * t268 + (-t271 * t218 + t272 * t220 + t275 * t226 + t276 * t228 + t426 * t295) * t267) * t267 / 0.2e1 + ((-t238 * t273 + t239 * t274 + t247 * t277 + t248 * t278 + t297 * t424) * t284 + (-t273 * t219 + t274 * t221 + t277 * t227 + t278 * t229 + t425 * t297) * t268 + (-t218 * t273 + t220 * t274 + t226 * t277 + t228 * t278 + t297 * t426) * t267) * t268 / 0.2e1 + ((-t238 * t286 + t239 * t287 + t247 * t293 + t248 * t294 + t424 * t318) * t284 + (-t219 * t286 + t221 * t287 + t227 * t293 + t229 * t294 + t318 * t425) * t268 + (-t218 * t286 + t220 * t287 + t226 * t293 + t228 * t294 + t318 * t426) * t267) * t284 / 0.2e1 + ((-t376 * t354 + t356 * t380 + Icges(1,4)) * V_base(5) + (-t376 * t355 + t357 * t380 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t354 * t380 + t376 * t356 + Icges(1,2)) * V_base(5) + (t355 * t380 + t376 * t357 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t376 + Icges(2,6) * t380) * V_base(5) + (Icges(2,5) * t380 - Icges(2,6) * t376) * V_base(4) + Icges(2,3) * t366 / 0.2e1) * t366;
T  = t1;
