% Calculate kinetic energy for
% S6RPRPRR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 04:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR9_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR9_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR9_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:01:09
% EndTime: 2019-03-09 04:01:14
% DurationCPUTime: 5.37s
% Computational Cost: add. (5934->459), mult. (15413->652), div. (0->0), fcn. (19872->16), ass. (0->201)
t426 = Icges(4,3) + Icges(5,3);
t366 = sin(pkin(12));
t371 = cos(pkin(6));
t378 = cos(qJ(1));
t369 = cos(pkin(12));
t375 = sin(qJ(1));
t401 = t375 * t369;
t341 = -t366 * t378 - t371 * t401;
t367 = sin(pkin(7));
t370 = cos(pkin(7));
t368 = sin(pkin(6));
t406 = t368 * t375;
t317 = -t341 * t367 + t370 * t406;
t425 = pkin(9) * t317;
t402 = t375 * t366;
t403 = t371 * t378;
t339 = t369 * t403 - t402;
t405 = t368 * t378;
t393 = t339 * t367 + t370 * t405;
t424 = pkin(9) * t393;
t374 = sin(qJ(3));
t377 = cos(qJ(3));
t412 = sin(pkin(13));
t413 = cos(pkin(13));
t345 = -t374 * t412 + t377 * t413;
t331 = t345 * t370;
t340 = t366 * t403 + t401;
t344 = -t374 * t413 - t377 * t412;
t389 = t367 * t345;
t386 = t368 * t389;
t278 = t331 * t339 + t340 * t344 - t378 * t386;
t330 = t344 * t367;
t332 = t344 * t370;
t279 = t330 * t405 - t339 * t332 + t340 * t345;
t392 = t339 * t370 - t367 * t405;
t291 = -t340 * t374 + t377 * t392;
t292 = t340 * t377 + t374 * t392;
t423 = Icges(4,5) * t292 + Icges(5,5) * t279 + Icges(4,6) * t291 + Icges(5,6) * t278 - t426 * t393;
t342 = t369 * t378 - t371 * t402;
t280 = t341 * t331 + t342 * t344 + t375 * t386;
t281 = -t330 * t406 - t332 * t341 + t342 * t345;
t391 = t341 * t370 + t367 * t406;
t293 = -t342 * t374 + t377 * t391;
t294 = t342 * t377 + t374 * t391;
t422 = Icges(4,5) * t294 + Icges(5,5) * t281 + Icges(4,6) * t293 + Icges(5,6) * t280 + t426 * t317;
t407 = t368 * t369;
t409 = t366 * t368;
t288 = t331 * t407 + t344 * t409 + t371 * t389;
t289 = -t330 * t371 + (-t332 * t369 + t345 * t366) * t368;
t314 = t367 * t371 * t377 + (t369 * t370 * t377 - t366 * t374) * t368;
t404 = t370 * t374;
t408 = t367 * t374;
t315 = t371 * t408 + (t366 * t377 + t369 * t404) * t368;
t338 = -t367 * t407 + t370 * t371;
t421 = Icges(4,5) * t315 + Icges(5,5) * t289 + Icges(4,6) * t314 + Icges(5,6) * t288 + t426 * t338;
t416 = cos(qJ(5));
t415 = pkin(3) * t377;
t414 = pkin(9) + qJ(4);
t411 = Icges(2,4) * t375;
t410 = qJ(2) * t371;
t400 = qJD(2) * t368;
t399 = V_base(5) * pkin(8) + V_base(1);
t307 = qJD(3) * t317 + V_base(4);
t306 = -qJD(3) * t393 + V_base(5);
t363 = V_base(6) + qJD(1);
t396 = -pkin(8) - t410;
t346 = t375 * pkin(1) - qJ(2) * t405;
t395 = qJD(2) * t371 + V_base(4) * t346 + V_base(3);
t265 = -qJD(5) * t280 + t307;
t264 = -qJD(5) * t278 + t306;
t324 = qJD(3) * t338 + t363;
t394 = t375 * t400 + V_base(5) * t410 + t399;
t273 = -qJD(5) * t288 + t324;
t347 = pkin(1) * t378 + qJ(2) * t406;
t390 = t363 * t347 - t378 * t400 + V_base(2);
t296 = t340 * pkin(2) - t424;
t320 = pkin(2) * t409 + pkin(9) * t338;
t388 = V_base(5) * t320 + (-t296 - t346) * t363 + t394;
t297 = pkin(2) * t342 + t425;
t387 = V_base(4) * t296 + (-t297 - t347) * V_base(5) + t395;
t336 = pkin(3) * t408 + t370 * t414;
t337 = pkin(3) * t404 - t367 * t414;
t286 = (-pkin(9) * t370 + t336) * t371 + ((pkin(9) * t367 + t337) * t369 + t415 * t366) * t368;
t385 = qJD(4) * t317 + t306 * t286 + t388;
t262 = -t336 * t405 + t339 * t337 + t340 * t415 + t424;
t384 = qJD(4) * t338 + t307 * t262 + t387;
t383 = t363 * t297 + (-t320 + t396) * V_base(4) + t390;
t263 = t336 * t406 + t337 * t341 + t342 * t415 - t425;
t382 = -qJD(4) * t393 + t324 * t263 + t383;
t241 = pkin(4) * t279 - pkin(10) * t278;
t260 = pkin(4) * t289 - pkin(10) * t288;
t381 = t306 * t260 + (-t241 - t262) * t324 + t385;
t242 = pkin(4) * t281 - pkin(10) * t280;
t380 = t307 * t241 + (-t242 - t263) * t306 + t384;
t379 = t324 * t242 + (-t260 - t286) * t307 + t382;
t376 = cos(qJ(6));
t373 = sin(qJ(5));
t372 = sin(qJ(6));
t364 = Icges(2,4) * t378;
t358 = rSges(2,1) * t378 - t375 * rSges(2,2);
t357 = t375 * rSges(2,1) + rSges(2,2) * t378;
t356 = Icges(2,1) * t378 - t411;
t355 = Icges(2,1) * t375 + t364;
t354 = -Icges(2,2) * t375 + t364;
t353 = Icges(2,2) * t378 + t411;
t352 = Icges(2,5) * t378 - Icges(2,6) * t375;
t351 = Icges(2,5) * t375 + Icges(2,6) * t378;
t350 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t349 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t348 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t328 = rSges(3,3) * t371 + (rSges(3,1) * t366 + rSges(3,2) * t369) * t368;
t327 = Icges(3,5) * t371 + (Icges(3,1) * t366 + Icges(3,4) * t369) * t368;
t326 = Icges(3,6) * t371 + (Icges(3,4) * t366 + Icges(3,2) * t369) * t368;
t325 = Icges(3,3) * t371 + (Icges(3,5) * t366 + Icges(3,6) * t369) * t368;
t322 = V_base(5) * rSges(2,3) - t357 * t363 + t399;
t321 = t358 * t363 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t319 = t357 * V_base(4) - t358 * V_base(5) + V_base(3);
t305 = rSges(3,1) * t342 + rSges(3,2) * t341 + rSges(3,3) * t406;
t304 = t340 * rSges(3,1) + t339 * rSges(3,2) - rSges(3,3) * t405;
t303 = Icges(3,1) * t342 + Icges(3,4) * t341 + Icges(3,5) * t406;
t302 = Icges(3,1) * t340 + Icges(3,4) * t339 - Icges(3,5) * t405;
t301 = Icges(3,4) * t342 + Icges(3,2) * t341 + Icges(3,6) * t406;
t300 = Icges(3,4) * t340 + Icges(3,2) * t339 - Icges(3,6) * t405;
t299 = Icges(3,5) * t342 + Icges(3,6) * t341 + Icges(3,3) * t406;
t298 = Icges(3,5) * t340 + Icges(3,6) * t339 - Icges(3,3) * t405;
t285 = rSges(4,1) * t315 + rSges(4,2) * t314 + rSges(4,3) * t338;
t284 = Icges(4,1) * t315 + Icges(4,4) * t314 + Icges(4,5) * t338;
t283 = Icges(4,4) * t315 + Icges(4,2) * t314 + Icges(4,6) * t338;
t277 = t289 * t416 + t338 * t373;
t276 = t289 * t373 - t338 * t416;
t272 = t328 * V_base(5) + (-t304 - t346) * t363 + t394;
t271 = t363 * t305 + (-t328 + t396) * V_base(4) + t390;
t269 = t281 * t416 + t317 * t373;
t268 = t281 * t373 - t317 * t416;
t267 = t279 * t416 - t373 * t393;
t266 = t279 * t373 + t393 * t416;
t261 = t304 * V_base(4) + (-t305 - t347) * V_base(5) + t395;
t259 = rSges(4,1) * t294 + rSges(4,2) * t293 + rSges(4,3) * t317;
t258 = rSges(4,1) * t292 + rSges(4,2) * t291 - rSges(4,3) * t393;
t257 = Icges(4,1) * t294 + Icges(4,4) * t293 + Icges(4,5) * t317;
t256 = Icges(4,1) * t292 + Icges(4,4) * t291 - Icges(4,5) * t393;
t255 = Icges(4,4) * t294 + Icges(4,2) * t293 + Icges(4,6) * t317;
t254 = Icges(4,4) * t292 + Icges(4,2) * t291 - Icges(4,6) * t393;
t250 = rSges(5,1) * t289 + rSges(5,2) * t288 + rSges(5,3) * t338;
t249 = Icges(5,1) * t289 + Icges(5,4) * t288 + Icges(5,5) * t338;
t248 = Icges(5,4) * t289 + Icges(5,2) * t288 + Icges(5,6) * t338;
t245 = t277 * t376 - t288 * t372;
t244 = -t277 * t372 - t288 * t376;
t240 = qJD(6) * t276 + t273;
t239 = pkin(5) * t277 + pkin(11) * t276;
t237 = rSges(5,1) * t281 + rSges(5,2) * t280 + rSges(5,3) * t317;
t236 = rSges(5,1) * t279 + rSges(5,2) * t278 - rSges(5,3) * t393;
t235 = Icges(5,1) * t281 + Icges(5,4) * t280 + Icges(5,5) * t317;
t234 = Icges(5,1) * t279 + Icges(5,4) * t278 - Icges(5,5) * t393;
t233 = Icges(5,4) * t281 + Icges(5,2) * t280 + Icges(5,6) * t317;
t232 = Icges(5,4) * t279 + Icges(5,2) * t278 - Icges(5,6) * t393;
t228 = t269 * t376 - t280 * t372;
t227 = -t269 * t372 - t280 * t376;
t226 = t267 * t376 - t278 * t372;
t225 = -t267 * t372 - t278 * t376;
t224 = qJD(6) * t268 + t265;
t223 = qJD(6) * t266 + t264;
t222 = rSges(6,1) * t277 - rSges(6,2) * t276 - rSges(6,3) * t288;
t221 = Icges(6,1) * t277 - Icges(6,4) * t276 - Icges(6,5) * t288;
t220 = Icges(6,4) * t277 - Icges(6,2) * t276 - Icges(6,6) * t288;
t219 = Icges(6,5) * t277 - Icges(6,6) * t276 - Icges(6,3) * t288;
t218 = pkin(5) * t269 + pkin(11) * t268;
t217 = pkin(5) * t267 + pkin(11) * t266;
t216 = rSges(6,1) * t269 - rSges(6,2) * t268 - rSges(6,3) * t280;
t215 = rSges(6,1) * t267 - rSges(6,2) * t266 - rSges(6,3) * t278;
t214 = Icges(6,1) * t269 - Icges(6,4) * t268 - Icges(6,5) * t280;
t213 = Icges(6,1) * t267 - Icges(6,4) * t266 - Icges(6,5) * t278;
t212 = Icges(6,4) * t269 - Icges(6,2) * t268 - Icges(6,6) * t280;
t211 = Icges(6,4) * t267 - Icges(6,2) * t266 - Icges(6,6) * t278;
t210 = Icges(6,5) * t269 - Icges(6,6) * t268 - Icges(6,3) * t280;
t209 = Icges(6,5) * t267 - Icges(6,6) * t266 - Icges(6,3) * t278;
t208 = -t258 * t324 + t285 * t306 + t388;
t207 = t324 * t259 - t307 * t285 + t383;
t206 = rSges(7,1) * t245 + rSges(7,2) * t244 + rSges(7,3) * t276;
t205 = Icges(7,1) * t245 + Icges(7,4) * t244 + Icges(7,5) * t276;
t204 = Icges(7,4) * t245 + Icges(7,2) * t244 + Icges(7,6) * t276;
t203 = Icges(7,5) * t245 + Icges(7,6) * t244 + Icges(7,3) * t276;
t202 = t258 * t307 - t259 * t306 + t387;
t201 = rSges(7,1) * t228 + rSges(7,2) * t227 + rSges(7,3) * t268;
t200 = rSges(7,1) * t226 + rSges(7,2) * t225 + rSges(7,3) * t266;
t199 = Icges(7,1) * t228 + Icges(7,4) * t227 + Icges(7,5) * t268;
t198 = Icges(7,1) * t226 + Icges(7,4) * t225 + Icges(7,5) * t266;
t197 = Icges(7,4) * t228 + Icges(7,2) * t227 + Icges(7,6) * t268;
t196 = Icges(7,4) * t226 + Icges(7,2) * t225 + Icges(7,6) * t266;
t195 = Icges(7,5) * t228 + Icges(7,6) * t227 + Icges(7,3) * t268;
t194 = Icges(7,5) * t226 + Icges(7,6) * t225 + Icges(7,3) * t266;
t193 = t250 * t306 + (-t236 - t262) * t324 + t385;
t192 = t324 * t237 + (-t250 - t286) * t307 + t382;
t191 = t236 * t307 + (-t237 - t263) * t306 + t384;
t190 = -t215 * t273 + t222 * t264 + t381;
t189 = t273 * t216 - t265 * t222 + t379;
t188 = t215 * t265 - t216 * t264 + t380;
t187 = -t200 * t240 + t206 * t223 - t217 * t273 + t239 * t264 + t381;
t186 = t240 * t201 - t224 * t206 + t273 * t218 - t265 * t239 + t379;
t185 = t200 * t224 - t201 * t223 + t217 * t265 - t218 * t264 + t380;
t1 = m(4) * (t202 ^ 2 + t207 ^ 2 + t208 ^ 2) / 0.2e1 + t265 * ((-t210 * t280 - t212 * t268 + t214 * t269) * t265 + (-t209 * t280 - t211 * t268 + t213 * t269) * t264 + (-t219 * t280 - t220 * t268 + t221 * t269) * t273) / 0.2e1 + t264 * ((-t210 * t278 - t212 * t266 + t214 * t267) * t265 + (-t209 * t278 - t211 * t266 + t213 * t267) * t264 + (-t219 * t278 - t220 * t266 + t221 * t267) * t273) / 0.2e1 + t240 * ((t195 * t276 + t197 * t244 + t199 * t245) * t224 + (t194 * t276 + t196 * t244 + t198 * t245) * t223 + (t203 * t276 + t204 * t244 + t205 * t245) * t240) / 0.2e1 + t224 * ((t195 * t268 + t197 * t227 + t199 * t228) * t224 + (t194 * t268 + t196 * t227 + t198 * t228) * t223 + (t203 * t268 + t204 * t227 + t205 * t228) * t240) / 0.2e1 + m(3) * (t261 ^ 2 + t271 ^ 2 + t272 ^ 2) / 0.2e1 + t223 * ((t195 * t266 + t197 * t225 + t199 * t226) * t224 + (t194 * t266 + t196 * t225 + t198 * t226) * t223 + (t203 * t266 + t204 * t225 + t205 * t226) * t240) / 0.2e1 + t273 * ((-t210 * t288 - t212 * t276 + t214 * t277) * t265 + (-t209 * t288 - t211 * t276 + t213 * t277) * t264 + (-t219 * t288 - t220 * t276 + t221 * t277) * t273) / 0.2e1 + m(2) * (t319 ^ 2 + t321 ^ 2 + t322 ^ 2) / 0.2e1 + m(5) * (t191 ^ 2 + t192 ^ 2 + t193 ^ 2) / 0.2e1 + m(6) * (t188 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + m(7) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(1) * (t348 ^ 2 + t349 ^ 2 + t350 ^ 2) / 0.2e1 + ((t248 * t278 + t249 * t279 + t283 * t291 + t284 * t292 - t393 * t421) * t324 + (t233 * t278 + t235 * t279 + t255 * t291 + t257 * t292 - t393 * t422) * t307 + (t232 * t278 + t234 * t279 + t254 * t291 + t256 * t292 - t393 * t423) * t306) * t306 / 0.2e1 + ((t248 * t280 + t249 * t281 + t283 * t293 + t284 * t294 + t317 * t421) * t324 + (t233 * t280 + t235 * t281 + t255 * t293 + t257 * t294 + t422 * t317) * t307 + (t232 * t280 + t234 * t281 + t254 * t293 + t256 * t294 + t317 * t423) * t306) * t307 / 0.2e1 + ((t248 * t288 + t249 * t289 + t283 * t314 + t284 * t315 + t421 * t338) * t324 + (t233 * t288 + t235 * t289 + t255 * t314 + t257 * t315 + t338 * t422) * t307 + (t232 * t288 + t234 * t289 + t254 * t314 + t256 * t315 + t338 * t423) * t306) * t324 / 0.2e1 + ((t298 * V_base(5) + t299 * V_base(4) + t325 * t363) * t371 + ((t301 * t369 + t303 * t366) * V_base(4) + (t300 * t369 + t302 * t366) * V_base(5) + (t326 * t369 + t327 * t366) * t363) * t368 + Icges(2,3) * t363 + t351 * V_base(5) + t352 * V_base(4)) * t363 / 0.2e1 + ((t325 * t406 + t326 * t341 + t327 * t342 + t352) * t363 + (t298 * t406 + t300 * t341 + t302 * t342 - t375 * t353 + t355 * t378 + Icges(1,4)) * V_base(5) + (t299 * t406 + t301 * t341 + t303 * t342 - t375 * t354 + t356 * t378 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t325 * t405 + t339 * t326 + t340 * t327 + t351) * t363 + (-t298 * t405 + t339 * t300 + t340 * t302 + t353 * t378 + t375 * t355 + Icges(1,2)) * V_base(5) + (-t299 * t405 + t339 * t301 + t340 * t303 + t354 * t378 + t375 * t356 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
