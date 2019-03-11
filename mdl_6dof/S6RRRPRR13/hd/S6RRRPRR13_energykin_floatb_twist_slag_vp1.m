% Calculate kinetic energy for
% S6RRRPRR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 20:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR13_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR13_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR13_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRR13_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRPRR13_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR13_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR13_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR13_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 19:51:50
% EndTime: 2019-03-09 19:51:55
% DurationCPUTime: 5.06s
% Computational Cost: add. (5380->442), mult. (12648->645), div. (0->0), fcn. (16040->16), ass. (0->194)
t427 = Icges(4,2) + Icges(5,3);
t372 = cos(pkin(6));
t379 = cos(qJ(2));
t380 = cos(qJ(1));
t405 = t379 * t380;
t376 = sin(qJ(2));
t377 = sin(qJ(1));
t407 = t377 * t376;
t337 = t372 * t405 - t407;
t406 = t377 * t379;
t408 = t376 * t380;
t338 = t372 * t408 + t406;
t375 = sin(qJ(3));
t370 = sin(pkin(6));
t416 = sin(pkin(7));
t397 = t370 * t416;
t417 = cos(pkin(7));
t420 = cos(qJ(3));
t296 = t338 * t420 + (t337 * t417 - t380 * t397) * t375;
t398 = t370 * t417;
t320 = -t337 * t416 - t380 * t398;
t369 = sin(pkin(13));
t371 = cos(pkin(13));
t275 = -t296 * t369 + t320 * t371;
t414 = t320 * t369;
t276 = t296 * t371 + t414;
t394 = t420 * t416;
t393 = t370 * t394;
t395 = t417 * t420;
t295 = -t337 * t395 + t338 * t375 + t380 * t393;
t426 = -Icges(4,4) * t296 + Icges(5,5) * t276 - Icges(4,6) * t320 + Icges(5,6) * t275 + t427 * t295;
t339 = -t372 * t406 - t408;
t340 = -t372 * t407 + t405;
t298 = t340 * t420 + (t339 * t417 + t377 * t397) * t375;
t321 = -t339 * t416 + t377 * t398;
t277 = -t298 * t369 + t321 * t371;
t413 = t321 * t369;
t278 = t298 * t371 + t413;
t297 = -t339 * t395 + t340 * t375 - t377 * t393;
t425 = -Icges(4,4) * t298 + Icges(5,5) * t278 - Icges(4,6) * t321 + Icges(5,6) * t277 + t427 * t297;
t319 = t372 * t416 * t375 + (t375 * t379 * t417 + t376 * t420) * t370;
t336 = t372 * t417 - t379 * t397;
t289 = -t319 * t369 + t336 * t371;
t412 = t336 * t369;
t290 = t319 * t371 + t412;
t411 = t370 * t376;
t318 = -t370 * t379 * t395 - t372 * t394 + t375 * t411;
t424 = -Icges(4,4) * t319 + Icges(5,5) * t290 - Icges(4,6) * t336 + Icges(5,6) * t289 + t427 * t318;
t419 = pkin(9) * t372;
t418 = pkin(4) * t371;
t415 = Icges(2,4) * t377;
t410 = t370 * t377;
t409 = t370 * t380;
t403 = qJD(2) * t370;
t402 = pkin(13) + qJ(5);
t401 = V_base(5) * pkin(8) + V_base(1);
t350 = t377 * t403 + V_base(4);
t366 = V_base(6) + qJD(1);
t396 = cos(t402);
t310 = qJD(3) * t321 + t350;
t351 = qJD(2) * t372 + t366;
t268 = qJD(5) * t297 + t310;
t322 = qJD(3) * t336 + t351;
t349 = -t380 * t403 + V_base(5);
t342 = t377 * pkin(1) - pkin(9) * t409;
t392 = -t342 * t366 + V_base(5) * t419 + t401;
t284 = qJD(5) * t318 + t322;
t343 = pkin(1) * t380 + pkin(9) * t410;
t391 = V_base(4) * t342 - t343 * V_base(5) + V_base(3);
t309 = qJD(3) * t320 + t349;
t267 = qJD(5) * t295 + t309;
t390 = t366 * t343 + V_base(2) + (-pkin(8) - t419) * V_base(4);
t299 = t338 * pkin(2) + pkin(10) * t320;
t326 = pkin(2) * t411 + pkin(10) * t336;
t389 = -t299 * t351 + t349 * t326 + t392;
t300 = t340 * pkin(2) + pkin(10) * t321;
t388 = t350 * t299 - t300 * t349 + t391;
t283 = pkin(3) * t319 + qJ(4) * t318;
t387 = qJD(4) * t297 + t309 * t283 + t389;
t264 = pkin(3) * t296 + qJ(4) * t295;
t386 = qJD(4) * t318 + t310 * t264 + t388;
t385 = t351 * t300 - t326 * t350 + t390;
t265 = pkin(3) * t298 + qJ(4) * t297;
t384 = qJD(4) * t295 + t322 * t265 + t385;
t209 = pkin(4) * t414 + pkin(11) * t295 + t296 * t418;
t236 = pkin(4) * t412 + pkin(11) * t318 + t319 * t418;
t383 = t309 * t236 + (-t209 - t264) * t322 + t387;
t210 = pkin(4) * t413 + pkin(11) * t297 + t298 * t418;
t382 = t310 * t209 + (-t210 - t265) * t309 + t386;
t381 = t322 * t210 + (-t236 - t283) * t310 + t384;
t378 = cos(qJ(6));
t374 = sin(qJ(6));
t367 = Icges(2,4) * t380;
t365 = sin(t402);
t359 = rSges(2,1) * t380 - t377 * rSges(2,2);
t358 = t377 * rSges(2,1) + rSges(2,2) * t380;
t357 = Icges(2,1) * t380 - t415;
t356 = Icges(2,1) * t377 + t367;
t355 = -Icges(2,2) * t377 + t367;
t354 = Icges(2,2) * t380 + t415;
t348 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t347 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t346 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t333 = rSges(3,3) * t372 + (rSges(3,1) * t376 + rSges(3,2) * t379) * t370;
t332 = Icges(3,5) * t372 + (Icges(3,1) * t376 + Icges(3,4) * t379) * t370;
t331 = Icges(3,6) * t372 + (Icges(3,4) * t376 + Icges(3,2) * t379) * t370;
t330 = Icges(3,3) * t372 + (Icges(3,5) * t376 + Icges(3,6) * t379) * t370;
t325 = V_base(5) * rSges(2,3) - t358 * t366 + t401;
t324 = t359 * t366 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t323 = t358 * V_base(4) - t359 * V_base(5) + V_base(3);
t308 = rSges(3,1) * t340 + rSges(3,2) * t339 + rSges(3,3) * t410;
t307 = t338 * rSges(3,1) + t337 * rSges(3,2) - rSges(3,3) * t409;
t306 = Icges(3,1) * t340 + Icges(3,4) * t339 + Icges(3,5) * t410;
t305 = Icges(3,1) * t338 + Icges(3,4) * t337 - Icges(3,5) * t409;
t304 = Icges(3,4) * t340 + Icges(3,2) * t339 + Icges(3,6) * t410;
t303 = Icges(3,4) * t338 + Icges(3,2) * t337 - Icges(3,6) * t409;
t302 = Icges(3,5) * t340 + Icges(3,6) * t339 + Icges(3,3) * t410;
t301 = Icges(3,5) * t338 + Icges(3,6) * t337 - Icges(3,3) * t409;
t287 = t319 * t396 + t336 * t365;
t286 = t319 * t365 - t336 * t396;
t282 = rSges(4,1) * t319 - rSges(4,2) * t318 + rSges(4,3) * t336;
t281 = Icges(4,1) * t319 - Icges(4,4) * t318 + Icges(4,5) * t336;
t279 = Icges(4,5) * t319 - Icges(4,6) * t318 + Icges(4,3) * t336;
t274 = t298 * t396 + t321 * t365;
t273 = t298 * t365 - t321 * t396;
t272 = t296 * t396 + t320 * t365;
t271 = t296 * t365 - t320 * t396;
t270 = t287 * t378 + t318 * t374;
t269 = -t287 * t374 + t318 * t378;
t263 = -t307 * t351 + t333 * t349 + t392;
t262 = t308 * t351 - t333 * t350 + t390;
t261 = pkin(5) * t287 + pkin(12) * t286;
t260 = qJD(6) * t286 + t284;
t258 = rSges(4,1) * t298 - rSges(4,2) * t297 + rSges(4,3) * t321;
t257 = rSges(4,1) * t296 - rSges(4,2) * t295 + rSges(4,3) * t320;
t256 = t307 * t350 - t308 * t349 + t391;
t255 = Icges(4,1) * t298 - Icges(4,4) * t297 + Icges(4,5) * t321;
t254 = Icges(4,1) * t296 - Icges(4,4) * t295 + Icges(4,5) * t320;
t251 = Icges(4,5) * t298 - Icges(4,6) * t297 + Icges(4,3) * t321;
t250 = Icges(4,5) * t296 - Icges(4,6) * t295 + Icges(4,3) * t320;
t249 = rSges(5,1) * t290 + rSges(5,2) * t289 + rSges(5,3) * t318;
t248 = Icges(5,1) * t290 + Icges(5,4) * t289 + Icges(5,5) * t318;
t247 = Icges(5,4) * t290 + Icges(5,2) * t289 + Icges(5,6) * t318;
t244 = t274 * t378 + t297 * t374;
t243 = -t274 * t374 + t297 * t378;
t242 = t272 * t378 + t295 * t374;
t241 = -t272 * t374 + t295 * t378;
t240 = rSges(6,1) * t287 - rSges(6,2) * t286 + rSges(6,3) * t318;
t239 = Icges(6,1) * t287 - Icges(6,4) * t286 + Icges(6,5) * t318;
t238 = Icges(6,4) * t287 - Icges(6,2) * t286 + Icges(6,6) * t318;
t237 = Icges(6,5) * t287 - Icges(6,6) * t286 + Icges(6,3) * t318;
t235 = pkin(5) * t274 + pkin(12) * t273;
t234 = pkin(5) * t272 + pkin(12) * t271;
t233 = qJD(6) * t273 + t268;
t232 = qJD(6) * t271 + t267;
t230 = rSges(5,1) * t278 + rSges(5,2) * t277 + rSges(5,3) * t297;
t229 = rSges(5,1) * t276 + rSges(5,2) * t275 + rSges(5,3) * t295;
t228 = Icges(5,1) * t278 + Icges(5,4) * t277 + Icges(5,5) * t297;
t227 = Icges(5,1) * t276 + Icges(5,4) * t275 + Icges(5,5) * t295;
t226 = Icges(5,4) * t278 + Icges(5,2) * t277 + Icges(5,6) * t297;
t225 = Icges(5,4) * t276 + Icges(5,2) * t275 + Icges(5,6) * t295;
t222 = rSges(6,1) * t274 - rSges(6,2) * t273 + rSges(6,3) * t297;
t221 = rSges(6,1) * t272 - rSges(6,2) * t271 + rSges(6,3) * t295;
t220 = Icges(6,1) * t274 - Icges(6,4) * t273 + Icges(6,5) * t297;
t219 = Icges(6,1) * t272 - Icges(6,4) * t271 + Icges(6,5) * t295;
t218 = Icges(6,4) * t274 - Icges(6,2) * t273 + Icges(6,6) * t297;
t217 = Icges(6,4) * t272 - Icges(6,2) * t271 + Icges(6,6) * t295;
t216 = Icges(6,5) * t274 - Icges(6,6) * t273 + Icges(6,3) * t297;
t215 = Icges(6,5) * t272 - Icges(6,6) * t271 + Icges(6,3) * t295;
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
t197 = t258 * t322 - t282 * t310 + t385;
t196 = t257 * t310 - t258 * t309 + t388;
t195 = t249 * t309 + (-t229 - t264) * t322 + t387;
t194 = t230 * t322 + (-t249 - t283) * t310 + t384;
t193 = t229 * t310 + (-t230 - t265) * t309 + t386;
t192 = -t221 * t284 + t240 * t267 + t383;
t191 = t222 * t284 - t240 * t268 + t381;
t190 = t221 * t268 - t222 * t267 + t382;
t189 = -t205 * t260 + t214 * t232 - t234 * t284 + t261 * t267 + t383;
t188 = t206 * t260 - t214 * t233 + t235 * t284 - t261 * t268 + t381;
t187 = t205 * t233 - t206 * t232 + t234 * t268 - t235 * t267 + t382;
t1 = t351 * ((t301 * t349 + t302 * t350 + t330 * t351) * t372 + ((t304 * t379 + t306 * t376) * t350 + (t303 * t379 + t305 * t376) * t349 + (t331 * t379 + t332 * t376) * t351) * t370) / 0.2e1 + t349 * ((-t302 * t409 + t337 * t304 + t338 * t306) * t350 + (-t301 * t409 + t337 * t303 + t338 * t305) * t349 + (-t330 * t409 + t337 * t331 + t338 * t332) * t351) / 0.2e1 + t350 * ((t302 * t410 + t304 * t339 + t306 * t340) * t350 + (t301 * t410 + t303 * t339 + t305 * t340) * t349 + (t330 * t410 + t331 * t339 + t332 * t340) * t351) / 0.2e1 + m(7) * (t187 ^ 2 + t188 ^ 2 + t189 ^ 2) / 0.2e1 + m(6) * (t190 ^ 2 + t191 ^ 2 + t192 ^ 2) / 0.2e1 + m(5) * (t193 ^ 2 + t194 ^ 2 + t195 ^ 2) / 0.2e1 + m(4) * (t196 ^ 2 + t197 ^ 2 + t198 ^ 2) / 0.2e1 + m(3) * (t256 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + t232 * ((t200 * t271 + t202 * t241 + t204 * t242) * t233 + (t271 * t199 + t241 * t201 + t203 * t242) * t232 + (t211 * t271 + t212 * t241 + t213 * t242) * t260) / 0.2e1 + t233 * ((t200 * t273 + t202 * t243 + t244 * t204) * t233 + (t199 * t273 + t201 * t243 + t203 * t244) * t232 + (t211 * t273 + t212 * t243 + t213 * t244) * t260) / 0.2e1 + t260 * ((t200 * t286 + t202 * t269 + t204 * t270) * t233 + (t199 * t286 + t201 * t269 + t203 * t270) * t232 + (t211 * t286 + t212 * t269 + t213 * t270) * t260) / 0.2e1 + t267 * ((t216 * t295 - t218 * t271 + t220 * t272) * t268 + (t215 * t295 - t271 * t217 + t272 * t219) * t267 + (t237 * t295 - t238 * t271 + t239 * t272) * t284) / 0.2e1 + t268 * ((t216 * t297 - t273 * t218 + t274 * t220) * t268 + (t215 * t297 - t217 * t273 + t219 * t274) * t267 + (t237 * t297 - t238 * t273 + t239 * t274) * t284) / 0.2e1 + t284 * ((t216 * t318 - t218 * t286 + t220 * t287) * t268 + (t215 * t318 - t217 * t286 + t219 * t287) * t267 + (t318 * t237 - t286 * t238 + t287 * t239) * t284) / 0.2e1 + m(2) * (t323 ^ 2 + t324 ^ 2 + t325 ^ 2) / 0.2e1 + m(1) * (t346 ^ 2 + t347 ^ 2 + t348 ^ 2) / 0.2e1 + ((t247 * t275 + t248 * t276 + t279 * t320 + t281 * t296 + t424 * t295) * t322 + (t226 * t275 + t228 * t276 + t251 * t320 + t255 * t296 + t425 * t295) * t310 + (t225 * t275 + t227 * t276 + t250 * t320 + t254 * t296 + t426 * t295) * t309) * t309 / 0.2e1 + ((t247 * t277 + t248 * t278 + t279 * t321 + t281 * t298 + t424 * t297) * t322 + (t226 * t277 + t228 * t278 + t251 * t321 + t255 * t298 + t425 * t297) * t310 + (t225 * t277 + t227 * t278 + t250 * t321 + t254 * t298 + t426 * t297) * t309) * t310 / 0.2e1 + ((t247 * t289 + t248 * t290 + t279 * t336 + t281 * t319 + t424 * t318) * t322 + (t226 * t289 + t228 * t290 + t251 * t336 + t255 * t319 + t425 * t318) * t310 + (t225 * t289 + t227 * t290 + t250 * t336 + t254 * t319 + t426 * t318) * t309) * t322 / 0.2e1 + ((-t377 * t354 + t356 * t380 + Icges(1,4)) * V_base(5) + (-t377 * t355 + t357 * t380 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t354 * t380 + t377 * t356 + Icges(1,2)) * V_base(5) + (t355 * t380 + t377 * t357 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t377 + Icges(2,6) * t380) * V_base(5) + (Icges(2,5) * t380 - Icges(2,6) * t377) * V_base(4) + Icges(2,3) * t366 / 0.2e1) * t366;
T  = t1;
