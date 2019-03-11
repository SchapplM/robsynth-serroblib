% Calculate kinetic energy for
% S6PPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRPRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PPRPRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRPRR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRPRR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:41:56
% EndTime: 2019-03-08 18:42:00
% DurationCPUTime: 4.72s
% Computational Cost: add. (5874->459), mult. (15413->647), div. (0->0), fcn. (19872->16), ass. (0->201)
t428 = Icges(4,3) + Icges(5,3);
t365 = sin(pkin(12));
t369 = cos(pkin(12));
t370 = cos(pkin(11));
t366 = sin(pkin(11));
t372 = cos(pkin(6));
t408 = t366 * t372;
t340 = -t365 * t370 - t369 * t408;
t367 = sin(pkin(7));
t368 = sin(pkin(6));
t371 = cos(pkin(7));
t404 = t368 * t371;
t318 = -t340 * t367 + t366 * t404;
t427 = pkin(8) * t318;
t403 = t370 * t372;
t338 = -t365 * t366 + t369 * t403;
t392 = t338 * t367 + t370 * t404;
t426 = pkin(8) * t392;
t425 = Icges(2,5) * t370 - Icges(2,6) * t366 + Icges(1,5);
t424 = Icges(2,5) * t366 + Icges(2,6) * t370 + Icges(1,6);
t375 = sin(qJ(3));
t377 = cos(qJ(3));
t413 = sin(pkin(13));
t414 = cos(pkin(13));
t345 = -t375 * t413 + t377 * t414;
t331 = t345 * t371;
t339 = t365 * t403 + t366 * t369;
t344 = -t375 * t414 - t377 * t413;
t388 = t367 * t345;
t387 = t368 * t388;
t275 = t331 * t338 + t339 * t344 - t370 * t387;
t330 = t344 * t367;
t332 = t344 * t371;
t405 = t368 * t370;
t276 = t330 * t405 - t332 * t338 + t339 * t345;
t391 = t338 * t371 - t367 * t405;
t290 = -t339 * t375 + t377 * t391;
t291 = t339 * t377 + t375 * t391;
t423 = Icges(4,5) * t291 + Icges(5,5) * t276 + Icges(4,6) * t290 + Icges(5,6) * t275 - t392 * t428;
t341 = -t365 * t408 + t369 * t370;
t277 = t340 * t331 + t341 * t344 + t366 * t387;
t409 = t366 * t368;
t278 = -t330 * t409 - t332 * t340 + t341 * t345;
t390 = t340 * t371 + t367 * t409;
t292 = -t341 * t375 + t377 * t390;
t293 = t341 * t377 + t375 * t390;
t422 = Icges(4,5) * t293 + Icges(5,5) * t278 + Icges(4,6) * t292 + Icges(5,6) * t277 + t318 * t428;
t406 = t368 * t369;
t410 = t365 * t368;
t288 = t331 * t406 + t344 * t410 + t372 * t388;
t289 = -t330 * t372 + (-t332 * t369 + t345 * t365) * t368;
t314 = t367 * t372 * t377 + (t369 * t371 * t377 - t365 * t375) * t368;
t402 = t371 * t375;
t407 = t367 * t375;
t315 = t372 * t407 + (t365 * t377 + t369 * t402) * t368;
t337 = -t367 * t406 + t371 * t372;
t421 = Icges(4,5) * t315 + Icges(5,5) * t289 + Icges(4,6) * t314 + Icges(5,6) * t288 + t337 * t428;
t417 = cos(qJ(5));
t416 = pkin(3) * t377;
t415 = pkin(8) + qJ(4);
t412 = Icges(2,4) * t366;
t411 = qJ(2) * t372;
t401 = qJD(2) * t368;
t400 = V_base(5) * qJ(1) + V_base(1);
t396 = qJD(1) + V_base(3);
t307 = qJD(3) * t318 + V_base(4);
t306 = -qJD(3) * t392 + V_base(5);
t324 = qJD(3) * t337 + V_base(6);
t395 = -qJ(1) - t411;
t265 = -qJD(5) * t277 + t307;
t264 = -qJD(5) * t275 + t306;
t279 = -qJD(5) * t288 + t324;
t394 = t366 * t401 + V_base(5) * t411 + t400;
t346 = pkin(1) * t366 - qJ(2) * t405;
t393 = qJD(2) * t372 + V_base(4) * t346 + t396;
t347 = pkin(1) * t370 + qJ(2) * t409;
t389 = V_base(6) * t347 - t370 * t401 + V_base(2);
t296 = pkin(2) * t339 - t426;
t320 = pkin(2) * t410 + pkin(8) * t337;
t386 = V_base(5) * t320 + (-t296 - t346) * V_base(6) + t394;
t297 = pkin(2) * t341 + t427;
t385 = V_base(4) * t296 + (-t297 - t347) * V_base(5) + t393;
t335 = pkin(3) * t407 + t371 * t415;
t336 = pkin(3) * t402 - t367 * t415;
t286 = (-pkin(8) * t371 + t335) * t372 + ((pkin(8) * t367 + t336) * t369 + t416 * t365) * t368;
t384 = qJD(4) * t318 + t306 * t286 + t386;
t262 = -t335 * t405 + t336 * t338 + t339 * t416 + t426;
t383 = qJD(4) * t337 + t307 * t262 + t385;
t382 = V_base(6) * t297 + (-t320 + t395) * V_base(4) + t389;
t263 = t335 * t409 + t336 * t340 + t341 * t416 - t427;
t381 = -qJD(4) * t392 + t324 * t263 + t382;
t239 = pkin(4) * t276 - pkin(9) * t275;
t260 = pkin(4) * t289 - pkin(9) * t288;
t380 = t306 * t260 + (-t239 - t262) * t324 + t384;
t240 = pkin(4) * t278 - pkin(9) * t277;
t379 = t307 * t239 + (-t240 - t263) * t306 + t383;
t378 = t324 * t240 + (-t260 - t286) * t307 + t381;
t376 = cos(qJ(6));
t374 = sin(qJ(5));
t373 = sin(qJ(6));
t363 = Icges(2,4) * t370;
t358 = rSges(2,1) * t370 - rSges(2,2) * t366;
t357 = rSges(2,1) * t366 + rSges(2,2) * t370;
t356 = Icges(2,1) * t370 - t412;
t355 = Icges(2,1) * t366 + t363;
t354 = -Icges(2,2) * t366 + t363;
t353 = Icges(2,2) * t370 + t412;
t350 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t349 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t348 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t328 = rSges(3,3) * t372 + (rSges(3,1) * t365 + rSges(3,2) * t369) * t368;
t327 = Icges(3,5) * t372 + (Icges(3,1) * t365 + Icges(3,4) * t369) * t368;
t326 = Icges(3,6) * t372 + (Icges(3,4) * t365 + Icges(3,2) * t369) * t368;
t325 = Icges(3,3) * t372 + (Icges(3,5) * t365 + Icges(3,6) * t369) * t368;
t323 = V_base(5) * rSges(2,3) - t357 * V_base(6) + t400;
t322 = t358 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t316 = t357 * V_base(4) - t358 * V_base(5) + t396;
t305 = rSges(3,1) * t341 + rSges(3,2) * t340 + rSges(3,3) * t409;
t304 = rSges(3,1) * t339 + rSges(3,2) * t338 - rSges(3,3) * t405;
t303 = Icges(3,1) * t341 + Icges(3,4) * t340 + Icges(3,5) * t409;
t302 = Icges(3,1) * t339 + Icges(3,4) * t338 - Icges(3,5) * t405;
t301 = Icges(3,4) * t341 + Icges(3,2) * t340 + Icges(3,6) * t409;
t300 = Icges(3,4) * t339 + Icges(3,2) * t338 - Icges(3,6) * t405;
t299 = Icges(3,5) * t341 + Icges(3,6) * t340 + Icges(3,3) * t409;
t298 = Icges(3,5) * t339 + Icges(3,6) * t338 - Icges(3,3) * t405;
t285 = rSges(4,1) * t315 + rSges(4,2) * t314 + rSges(4,3) * t337;
t284 = Icges(4,1) * t315 + Icges(4,4) * t314 + Icges(4,5) * t337;
t283 = Icges(4,4) * t315 + Icges(4,2) * t314 + Icges(4,6) * t337;
t281 = t289 * t417 + t337 * t374;
t280 = t289 * t374 - t337 * t417;
t272 = t328 * V_base(5) + (-t304 - t346) * V_base(6) + t394;
t271 = t305 * V_base(6) + (-t328 + t395) * V_base(4) + t389;
t269 = t278 * t417 + t318 * t374;
t268 = t278 * t374 - t318 * t417;
t267 = t276 * t417 - t374 * t392;
t266 = t276 * t374 + t392 * t417;
t261 = t304 * V_base(4) + (-t305 - t347) * V_base(5) + t393;
t259 = rSges(4,1) * t293 + rSges(4,2) * t292 + rSges(4,3) * t318;
t258 = rSges(4,1) * t291 + rSges(4,2) * t290 - rSges(4,3) * t392;
t257 = Icges(4,1) * t293 + Icges(4,4) * t292 + Icges(4,5) * t318;
t256 = Icges(4,1) * t291 + Icges(4,4) * t290 - Icges(4,5) * t392;
t255 = Icges(4,4) * t293 + Icges(4,2) * t292 + Icges(4,6) * t318;
t254 = Icges(4,4) * t291 + Icges(4,2) * t290 - Icges(4,6) * t392;
t250 = rSges(5,1) * t289 + rSges(5,2) * t288 + rSges(5,3) * t337;
t249 = Icges(5,1) * t289 + Icges(5,4) * t288 + Icges(5,5) * t337;
t248 = Icges(5,4) * t289 + Icges(5,2) * t288 + Icges(5,6) * t337;
t246 = t281 * t376 - t288 * t373;
t245 = -t281 * t373 - t288 * t376;
t242 = qJD(6) * t280 + t279;
t241 = pkin(5) * t281 + pkin(10) * t280;
t237 = rSges(5,1) * t278 + rSges(5,2) * t277 + rSges(5,3) * t318;
t236 = rSges(5,1) * t276 + rSges(5,2) * t275 - rSges(5,3) * t392;
t235 = Icges(5,1) * t278 + Icges(5,4) * t277 + Icges(5,5) * t318;
t234 = Icges(5,1) * t276 + Icges(5,4) * t275 - Icges(5,5) * t392;
t233 = Icges(5,4) * t278 + Icges(5,2) * t277 + Icges(5,6) * t318;
t232 = Icges(5,4) * t276 + Icges(5,2) * t275 - Icges(5,6) * t392;
t228 = t269 * t376 - t277 * t373;
t227 = -t269 * t373 - t277 * t376;
t226 = t267 * t376 - t275 * t373;
t225 = -t267 * t373 - t275 * t376;
t224 = rSges(6,1) * t281 - rSges(6,2) * t280 - rSges(6,3) * t288;
t223 = Icges(6,1) * t281 - Icges(6,4) * t280 - Icges(6,5) * t288;
t222 = Icges(6,4) * t281 - Icges(6,2) * t280 - Icges(6,6) * t288;
t221 = Icges(6,5) * t281 - Icges(6,6) * t280 - Icges(6,3) * t288;
t220 = qJD(6) * t268 + t265;
t219 = qJD(6) * t266 + t264;
t218 = pkin(5) * t269 + pkin(10) * t268;
t217 = pkin(5) * t267 + pkin(10) * t266;
t216 = rSges(6,1) * t269 - rSges(6,2) * t268 - rSges(6,3) * t277;
t215 = rSges(6,1) * t267 - rSges(6,2) * t266 - rSges(6,3) * t275;
t214 = Icges(6,1) * t269 - Icges(6,4) * t268 - Icges(6,5) * t277;
t213 = Icges(6,1) * t267 - Icges(6,4) * t266 - Icges(6,5) * t275;
t212 = Icges(6,4) * t269 - Icges(6,2) * t268 - Icges(6,6) * t277;
t211 = Icges(6,4) * t267 - Icges(6,2) * t266 - Icges(6,6) * t275;
t210 = Icges(6,5) * t269 - Icges(6,6) * t268 - Icges(6,3) * t277;
t209 = Icges(6,5) * t267 - Icges(6,6) * t266 - Icges(6,3) * t275;
t208 = -t258 * t324 + t285 * t306 + t386;
t207 = t259 * t324 - t285 * t307 + t382;
t206 = rSges(7,1) * t246 + rSges(7,2) * t245 + rSges(7,3) * t280;
t205 = Icges(7,1) * t246 + Icges(7,4) * t245 + Icges(7,5) * t280;
t204 = Icges(7,4) * t246 + Icges(7,2) * t245 + Icges(7,6) * t280;
t203 = Icges(7,5) * t246 + Icges(7,6) * t245 + Icges(7,3) * t280;
t202 = t258 * t307 - t259 * t306 + t385;
t201 = rSges(7,1) * t228 + rSges(7,2) * t227 + rSges(7,3) * t268;
t200 = rSges(7,1) * t226 + rSges(7,2) * t225 + rSges(7,3) * t266;
t199 = Icges(7,1) * t228 + Icges(7,4) * t227 + Icges(7,5) * t268;
t198 = Icges(7,1) * t226 + Icges(7,4) * t225 + Icges(7,5) * t266;
t197 = Icges(7,4) * t228 + Icges(7,2) * t227 + Icges(7,6) * t268;
t196 = Icges(7,4) * t226 + Icges(7,2) * t225 + Icges(7,6) * t266;
t195 = Icges(7,5) * t228 + Icges(7,6) * t227 + Icges(7,3) * t268;
t194 = Icges(7,5) * t226 + Icges(7,6) * t225 + Icges(7,3) * t266;
t193 = t250 * t306 + (-t236 - t262) * t324 + t384;
t192 = t237 * t324 + (-t250 - t286) * t307 + t381;
t191 = t236 * t307 + (-t237 - t263) * t306 + t383;
t190 = -t215 * t279 + t224 * t264 + t380;
t189 = t216 * t279 - t224 * t265 + t378;
t188 = t215 * t265 - t216 * t264 + t379;
t187 = -t200 * t242 + t206 * t219 - t217 * t279 + t241 * t264 + t380;
t186 = t201 * t242 - t206 * t220 + t218 * t279 - t241 * t265 + t378;
t185 = t200 * t220 - t201 * t219 + t217 * t265 - t218 * t264 + t379;
t1 = m(1) * (t348 ^ 2 + t349 ^ 2 + t350 ^ 2) / 0.2e1 + m(4) * (t202 ^ 2 + t207 ^ 2 + t208 ^ 2) / 0.2e1 + m(2) * (t316 ^ 2 + t322 ^ 2 + t323 ^ 2) / 0.2e1 + t220 * ((t195 * t268 + t197 * t227 + t199 * t228) * t220 + (t194 * t268 + t196 * t227 + t198 * t228) * t219 + (t203 * t268 + t204 * t227 + t205 * t228) * t242) / 0.2e1 + m(3) * (t261 ^ 2 + t271 ^ 2 + t272 ^ 2) / 0.2e1 + t219 * ((t195 * t266 + t197 * t225 + t199 * t226) * t220 + (t194 * t266 + t196 * t225 + t198 * t226) * t219 + (t203 * t266 + t204 * t225 + t205 * t226) * t242) / 0.2e1 + t279 * ((-t210 * t288 - t212 * t280 + t214 * t281) * t265 + (-t209 * t288 - t211 * t280 + t213 * t281) * t264 + (-t221 * t288 - t222 * t280 + t223 * t281) * t279) / 0.2e1 + t242 * ((t195 * t280 + t197 * t245 + t199 * t246) * t220 + (t194 * t280 + t196 * t245 + t198 * t246) * t219 + (t203 * t280 + t204 * t245 + t205 * t246) * t242) / 0.2e1 + t264 * ((-t210 * t275 - t212 * t266 + t214 * t267) * t265 + (-t209 * t275 - t211 * t266 + t213 * t267) * t264 + (-t221 * t275 - t222 * t266 + t223 * t267) * t279) / 0.2e1 + t265 * ((-t210 * t277 - t212 * t268 + t214 * t269) * t265 + (-t209 * t277 - t211 * t268 + t213 * t269) * t264 + (-t221 * t277 - t222 * t268 + t223 * t269) * t279) / 0.2e1 + m(5) * (t191 ^ 2 + t192 ^ 2 + t193 ^ 2) / 0.2e1 + m(6) * (t188 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + m(7) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + ((t248 * t275 + t249 * t276 + t283 * t290 + t284 * t291 - t392 * t421) * t324 + (t233 * t275 + t235 * t276 + t255 * t290 + t257 * t291 - t392 * t422) * t307 + (t232 * t275 + t234 * t276 + t254 * t290 + t256 * t291 - t392 * t423) * t306) * t306 / 0.2e1 + ((t248 * t277 + t249 * t278 + t283 * t292 + t284 * t293 + t421 * t318) * t324 + (t233 * t277 + t235 * t278 + t255 * t292 + t257 * t293 + t422 * t318) * t307 + (t232 * t277 + t234 * t278 + t254 * t292 + t256 * t293 + t423 * t318) * t306) * t307 / 0.2e1 + ((t248 * t288 + t249 * t289 + t283 * t314 + t284 * t315 + t421 * t337) * t324 + (t233 * t288 + t235 * t289 + t255 * t314 + t257 * t315 + t422 * t337) * t307 + (t232 * t288 + t234 * t289 + t254 * t314 + t256 * t315 + t423 * t337) * t306) * t324 / 0.2e1 + ((t325 * t409 + t326 * t340 + t327 * t341 + t425) * V_base(6) + (t298 * t409 + t300 * t340 + t302 * t341 - t353 * t366 + t355 * t370 + Icges(1,4)) * V_base(5) + (t299 * t409 + t301 * t340 + t303 * t341 - t354 * t366 + t356 * t370 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t325 * t405 + t326 * t338 + t327 * t339 + t424) * V_base(6) + (-t298 * t405 + t300 * t338 + t302 * t339 + t353 * t370 + t355 * t366 + Icges(1,2)) * V_base(5) + (-t299 * t405 + t301 * t338 + t303 * t339 + t354 * t370 + t356 * t366 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t325 * t372 + (t326 * t369 + t327 * t365) * t368 + Icges(2,3) + Icges(1,3)) * V_base(6) + (t298 * t372 + (t300 * t369 + t302 * t365) * t368 + t424) * V_base(5) + (t299 * t372 + (t301 * t369 + t303 * t365) * t368 + t425) * V_base(4)) * V_base(6) / 0.2e1;
T  = t1;
