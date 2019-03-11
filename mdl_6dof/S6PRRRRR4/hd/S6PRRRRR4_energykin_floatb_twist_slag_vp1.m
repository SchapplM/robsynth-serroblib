% Calculate kinetic energy for
% S6PRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 01:02
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR4_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:54:45
% EndTime: 2019-03-09 00:54:49
% DurationCPUTime: 4.39s
% Computational Cost: add. (5593->447), mult. (13362->673), div. (0->0), fcn. (16964->16), ass. (0->197)
t420 = cos(qJ(3));
t372 = cos(pkin(6));
t419 = pkin(8) * t372;
t378 = cos(qJ(4));
t418 = pkin(4) * t378;
t416 = cos(pkin(7));
t415 = sin(pkin(7));
t369 = sin(pkin(13));
t414 = Icges(2,4) * t369;
t371 = cos(pkin(13));
t376 = sin(qJ(2));
t379 = cos(qJ(2));
t406 = t372 * t379;
t337 = -t369 * t376 + t371 * t406;
t370 = sin(pkin(6));
t397 = t370 * t416;
t322 = -t337 * t415 - t371 * t397;
t374 = sin(qJ(4));
t413 = t322 * t374;
t339 = -t369 * t406 - t371 * t376;
t323 = -t339 * t415 + t369 * t397;
t412 = t323 * t374;
t396 = t370 * t415;
t336 = t372 * t416 - t379 * t396;
t411 = t336 * t374;
t410 = t369 * t370;
t409 = t370 * t371;
t408 = t370 * t376;
t407 = t372 * t376;
t405 = qJ(4) + qJ(5);
t404 = qJD(2) * t370;
t403 = V_base(5) * qJ(1) + V_base(1);
t399 = qJD(1) + V_base(3);
t351 = t369 * t404 + V_base(4);
t361 = qJD(2) * t372 + V_base(6);
t398 = cos(t405);
t311 = qJD(3) * t323 + t351;
t324 = qJD(3) * t336 + t361;
t395 = t416 * t420;
t394 = t420 * t415;
t340 = -t369 * t407 + t371 * t379;
t375 = sin(qJ(3));
t393 = t370 * t394;
t296 = -t339 * t395 + t340 * t375 - t369 * t393;
t268 = qJD(4) * t296 + t311;
t320 = -t370 * t379 * t395 - t372 * t394 + t375 * t408;
t286 = qJD(4) * t320 + t324;
t350 = -t371 * t404 + V_base(5);
t240 = qJD(5) * t296 + t268;
t275 = qJD(5) * t320 + t286;
t310 = qJD(3) * t322 + t350;
t343 = pkin(1) * t369 - pkin(8) * t409;
t392 = -t343 * V_base(6) + V_base(5) * t419 + t403;
t344 = pkin(1) * t371 + pkin(8) * t410;
t391 = V_base(4) * t343 - t344 * V_base(5) + t399;
t338 = t369 * t379 + t371 * t407;
t294 = -t337 * t395 + t338 * t375 + t371 * t393;
t267 = qJD(4) * t294 + t310;
t239 = qJD(5) * t294 + t267;
t390 = V_base(6) * t344 + V_base(2) + (-qJ(1) - t419) * V_base(4);
t300 = t338 * pkin(2) + pkin(9) * t322;
t325 = pkin(2) * t408 + pkin(9) * t336;
t389 = -t300 * t361 + t350 * t325 + t392;
t301 = t340 * pkin(2) + pkin(9) * t323;
t388 = t351 * t300 - t301 * t350 + t391;
t387 = t361 * t301 - t325 * t351 + t390;
t295 = t338 * t420 + (t337 * t416 - t371 * t396) * t375;
t264 = t295 * pkin(3) + t294 * pkin(10);
t321 = t372 * t415 * t375 + (t416 * t375 * t379 + t376 * t420) * t370;
t284 = t321 * pkin(3) + t320 * pkin(10);
t386 = -t264 * t324 + t310 * t284 + t389;
t297 = t340 * t420 + (t339 * t416 + t369 * t396) * t375;
t265 = t297 * pkin(3) + t296 * pkin(10);
t385 = t311 * t264 - t265 * t310 + t388;
t384 = t324 * t265 - t284 * t311 + t387;
t207 = pkin(4) * t413 + pkin(11) * t294 + t295 * t418;
t241 = pkin(4) * t411 + pkin(11) * t320 + t321 * t418;
t383 = -t207 * t286 + t267 * t241 + t386;
t208 = pkin(4) * t412 + pkin(11) * t296 + t297 * t418;
t382 = t268 * t207 - t208 * t267 + t385;
t381 = t286 * t208 - t241 * t268 + t384;
t377 = cos(qJ(6));
t373 = sin(qJ(6));
t367 = sin(t405);
t366 = Icges(2,4) * t371;
t359 = rSges(2,1) * t371 - rSges(2,2) * t369;
t358 = rSges(2,1) * t369 + rSges(2,2) * t371;
t357 = Icges(2,1) * t371 - t414;
t356 = Icges(2,1) * t369 + t366;
t355 = -Icges(2,2) * t369 + t366;
t354 = Icges(2,2) * t371 + t414;
t349 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t348 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t347 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t334 = rSges(3,3) * t372 + (rSges(3,1) * t376 + rSges(3,2) * t379) * t370;
t333 = Icges(3,5) * t372 + (Icges(3,1) * t376 + Icges(3,4) * t379) * t370;
t332 = Icges(3,6) * t372 + (Icges(3,4) * t376 + Icges(3,2) * t379) * t370;
t331 = Icges(3,3) * t372 + (Icges(3,5) * t376 + Icges(3,6) * t379) * t370;
t327 = V_base(5) * rSges(2,3) - t358 * V_base(6) + t403;
t326 = t359 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t319 = t358 * V_base(4) - t359 * V_base(5) + t399;
t309 = rSges(3,1) * t340 + rSges(3,2) * t339 + rSges(3,3) * t410;
t308 = rSges(3,1) * t338 + rSges(3,2) * t337 - rSges(3,3) * t409;
t307 = Icges(3,1) * t340 + Icges(3,4) * t339 + Icges(3,5) * t410;
t306 = Icges(3,1) * t338 + Icges(3,4) * t337 - Icges(3,5) * t409;
t305 = Icges(3,4) * t340 + Icges(3,2) * t339 + Icges(3,6) * t410;
t304 = Icges(3,4) * t338 + Icges(3,2) * t337 - Icges(3,6) * t409;
t303 = Icges(3,5) * t340 + Icges(3,6) * t339 + Icges(3,3) * t410;
t302 = Icges(3,5) * t338 + Icges(3,6) * t337 - Icges(3,3) * t409;
t299 = t321 * t378 + t411;
t298 = -t321 * t374 + t336 * t378;
t289 = t321 * t398 + t336 * t367;
t288 = t321 * t367 - t336 * t398;
t283 = rSges(4,1) * t321 - rSges(4,2) * t320 + rSges(4,3) * t336;
t282 = Icges(4,1) * t321 - Icges(4,4) * t320 + Icges(4,5) * t336;
t281 = Icges(4,4) * t321 - Icges(4,2) * t320 + Icges(4,6) * t336;
t280 = Icges(4,5) * t321 - Icges(4,6) * t320 + Icges(4,3) * t336;
t279 = t297 * t378 + t412;
t278 = -t297 * t374 + t323 * t378;
t277 = t295 * t378 + t413;
t276 = -t295 * t374 + t322 * t378;
t274 = t297 * t398 + t323 * t367;
t273 = t297 * t367 - t323 * t398;
t272 = t295 * t398 + t322 * t367;
t271 = t295 * t367 - t322 * t398;
t270 = t289 * t377 + t320 * t373;
t269 = -t289 * t373 + t320 * t377;
t263 = -t308 * t361 + t334 * t350 + t392;
t262 = t309 * t361 - t334 * t351 + t390;
t261 = pkin(5) * t289 + pkin(12) * t288;
t259 = rSges(5,1) * t299 + rSges(5,2) * t298 + rSges(5,3) * t320;
t258 = rSges(4,1) * t297 - rSges(4,2) * t296 + rSges(4,3) * t323;
t257 = rSges(4,1) * t295 - rSges(4,2) * t294 + rSges(4,3) * t322;
t256 = Icges(5,1) * t299 + Icges(5,4) * t298 + Icges(5,5) * t320;
t255 = Icges(5,4) * t299 + Icges(5,2) * t298 + Icges(5,6) * t320;
t254 = Icges(5,5) * t299 + Icges(5,6) * t298 + Icges(5,3) * t320;
t253 = Icges(4,1) * t297 - Icges(4,4) * t296 + Icges(4,5) * t323;
t252 = Icges(4,1) * t295 - Icges(4,4) * t294 + Icges(4,5) * t322;
t251 = Icges(4,4) * t297 - Icges(4,2) * t296 + Icges(4,6) * t323;
t250 = Icges(4,4) * t295 - Icges(4,2) * t294 + Icges(4,6) * t322;
t249 = Icges(4,5) * t297 - Icges(4,6) * t296 + Icges(4,3) * t323;
t248 = Icges(4,5) * t295 - Icges(4,6) * t294 + Icges(4,3) * t322;
t247 = t308 * t351 - t309 * t350 + t391;
t246 = qJD(6) * t288 + t275;
t245 = rSges(6,1) * t289 - rSges(6,2) * t288 + rSges(6,3) * t320;
t244 = Icges(6,1) * t289 - Icges(6,4) * t288 + Icges(6,5) * t320;
t243 = Icges(6,4) * t289 - Icges(6,2) * t288 + Icges(6,6) * t320;
t242 = Icges(6,5) * t289 - Icges(6,6) * t288 + Icges(6,3) * t320;
t237 = t274 * t377 + t296 * t373;
t236 = -t274 * t373 + t296 * t377;
t235 = t272 * t377 + t294 * t373;
t234 = -t272 * t373 + t294 * t377;
t233 = pkin(5) * t274 + pkin(12) * t273;
t232 = pkin(5) * t272 + pkin(12) * t271;
t231 = rSges(5,1) * t279 + rSges(5,2) * t278 + rSges(5,3) * t296;
t230 = rSges(5,1) * t277 + rSges(5,2) * t276 + rSges(5,3) * t294;
t229 = Icges(5,1) * t279 + Icges(5,4) * t278 + Icges(5,5) * t296;
t228 = Icges(5,1) * t277 + Icges(5,4) * t276 + Icges(5,5) * t294;
t227 = Icges(5,4) * t279 + Icges(5,2) * t278 + Icges(5,6) * t296;
t226 = Icges(5,4) * t277 + Icges(5,2) * t276 + Icges(5,6) * t294;
t225 = Icges(5,5) * t279 + Icges(5,6) * t278 + Icges(5,3) * t296;
t224 = Icges(5,5) * t277 + Icges(5,6) * t276 + Icges(5,3) * t294;
t223 = rSges(6,1) * t274 - rSges(6,2) * t273 + rSges(6,3) * t296;
t222 = rSges(6,1) * t272 - rSges(6,2) * t271 + rSges(6,3) * t294;
t221 = Icges(6,1) * t274 - Icges(6,4) * t273 + Icges(6,5) * t296;
t220 = Icges(6,1) * t272 - Icges(6,4) * t271 + Icges(6,5) * t294;
t219 = Icges(6,4) * t274 - Icges(6,2) * t273 + Icges(6,6) * t296;
t218 = Icges(6,4) * t272 - Icges(6,2) * t271 + Icges(6,6) * t294;
t217 = Icges(6,5) * t274 - Icges(6,6) * t273 + Icges(6,3) * t296;
t216 = Icges(6,5) * t272 - Icges(6,6) * t271 + Icges(6,3) * t294;
t215 = qJD(6) * t273 + t240;
t214 = qJD(6) * t271 + t239;
t213 = rSges(7,1) * t270 + rSges(7,2) * t269 + rSges(7,3) * t288;
t212 = Icges(7,1) * t270 + Icges(7,4) * t269 + Icges(7,5) * t288;
t211 = Icges(7,4) * t270 + Icges(7,2) * t269 + Icges(7,6) * t288;
t210 = Icges(7,5) * t270 + Icges(7,6) * t269 + Icges(7,3) * t288;
t204 = rSges(7,1) * t237 + rSges(7,2) * t236 + rSges(7,3) * t273;
t203 = rSges(7,1) * t235 + rSges(7,2) * t234 + rSges(7,3) * t271;
t202 = Icges(7,1) * t237 + Icges(7,4) * t236 + Icges(7,5) * t273;
t201 = Icges(7,1) * t235 + Icges(7,4) * t234 + Icges(7,5) * t271;
t200 = Icges(7,4) * t237 + Icges(7,2) * t236 + Icges(7,6) * t273;
t199 = Icges(7,4) * t235 + Icges(7,2) * t234 + Icges(7,6) * t271;
t198 = Icges(7,5) * t237 + Icges(7,6) * t236 + Icges(7,3) * t273;
t197 = Icges(7,5) * t235 + Icges(7,6) * t234 + Icges(7,3) * t271;
t196 = -t257 * t324 + t283 * t310 + t389;
t195 = t258 * t324 - t283 * t311 + t387;
t194 = t257 * t311 - t258 * t310 + t388;
t193 = -t230 * t286 + t259 * t267 + t386;
t192 = t231 * t286 - t259 * t268 + t384;
t191 = t230 * t268 - t231 * t267 + t385;
t190 = -t222 * t275 + t239 * t245 + t383;
t189 = t223 * t275 - t240 * t245 + t381;
t188 = t222 * t240 - t223 * t239 + t382;
t187 = -t203 * t246 + t213 * t214 - t232 * t275 + t239 * t261 + t383;
t186 = t204 * t246 - t213 * t215 + t233 * t275 - t240 * t261 + t381;
t185 = t203 * t215 - t204 * t214 + t232 * t240 - t233 * t239 + t382;
t1 = t361 * ((t302 * t350 + t303 * t351 + t331 * t361) * t372 + ((t305 * t379 + t307 * t376) * t351 + (t304 * t379 + t306 * t376) * t350 + (t332 * t379 + t333 * t376) * t361) * t370) / 0.2e1 + ((-t354 * t369 + t356 * t371 + Icges(1,4)) * V_base(5) + (-t355 * t369 + t357 * t371 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + m(1) * (t347 ^ 2 + t348 ^ 2 + t349 ^ 2) / 0.2e1 + t324 * ((t249 * t336 - t251 * t320 + t253 * t321) * t311 + (t248 * t336 - t250 * t320 + t252 * t321) * t310 + (t280 * t336 - t281 * t320 + t282 * t321) * t324) / 0.2e1 + t311 * ((t249 * t323 - t251 * t296 + t253 * t297) * t311 + (t248 * t323 - t250 * t296 + t252 * t297) * t310 + (t280 * t323 - t281 * t296 + t282 * t297) * t324) / 0.2e1 + t310 * ((t249 * t322 - t251 * t294 + t253 * t295) * t311 + (t248 * t322 - t250 * t294 + t252 * t295) * t310 + (t280 * t322 - t281 * t294 + t282 * t295) * t324) / 0.2e1 + m(2) * (t319 ^ 2 + t326 ^ 2 + t327 ^ 2) / 0.2e1 + t286 * ((t225 * t320 + t227 * t298 + t229 * t299) * t268 + (t224 * t320 + t226 * t298 + t228 * t299) * t267 + (t254 * t320 + t255 * t298 + t256 * t299) * t286) / 0.2e1 + t275 * ((t217 * t320 - t219 * t288 + t221 * t289) * t240 + (t216 * t320 - t218 * t288 + t220 * t289) * t239 + (t242 * t320 - t243 * t288 + t244 * t289) * t275) / 0.2e1 + t268 * ((t296 * t225 + t278 * t227 + t279 * t229) * t268 + (t224 * t296 + t226 * t278 + t228 * t279) * t267 + (t254 * t296 + t255 * t278 + t256 * t279) * t286) / 0.2e1 + t240 * ((t296 * t217 - t273 * t219 + t274 * t221) * t240 + (t216 * t296 - t218 * t273 + t220 * t274) * t239 + (t242 * t296 - t243 * t273 + t244 * t274) * t275) / 0.2e1 + t267 * ((t225 * t294 + t227 * t276 + t229 * t277) * t268 + (t294 * t224 + t276 * t226 + t277 * t228) * t267 + (t254 * t294 + t255 * t276 + t256 * t277) * t286) / 0.2e1 + t239 * ((t217 * t294 - t219 * t271 + t221 * t272) * t240 + (t294 * t216 - t271 * t218 + t272 * t220) * t239 + (t242 * t294 - t243 * t271 + t244 * t272) * t275) / 0.2e1 + t246 * ((t198 * t288 + t200 * t269 + t202 * t270) * t215 + (t197 * t288 + t199 * t269 + t201 * t270) * t214 + (t288 * t210 + t269 * t211 + t270 * t212) * t246) / 0.2e1 + t215 * ((t273 * t198 + t236 * t200 + t237 * t202) * t215 + (t197 * t273 + t199 * t236 + t201 * t237) * t214 + (t210 * t273 + t211 * t236 + t212 * t237) * t246) / 0.2e1 + t214 * ((t198 * t271 + t200 * t234 + t202 * t235) * t215 + (t271 * t197 + t234 * t199 + t235 * t201) * t214 + (t210 * t271 + t211 * t234 + t212 * t235) * t246) / 0.2e1 + m(3) * (t247 ^ 2 + t262 ^ 2 + t263 ^ 2) / 0.2e1 + m(4) * (t194 ^ 2 + t195 ^ 2 + t196 ^ 2) / 0.2e1 + ((Icges(2,5) * t371 - Icges(2,6) * t369 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t369 + Icges(2,6) * t371 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + m(5) * (t191 ^ 2 + t192 ^ 2 + t193 ^ 2) / 0.2e1 + m(6) * (t188 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + m(7) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + ((t354 * t371 + t356 * t369 + Icges(1,2)) * V_base(5) + (t355 * t371 + t357 * t369 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + t350 * ((-t303 * t409 + t305 * t337 + t307 * t338) * t351 + (-t302 * t409 + t304 * t337 + t306 * t338) * t350 + (-t331 * t409 + t332 * t337 + t333 * t338) * t361) / 0.2e1 + t351 * ((t303 * t410 + t305 * t339 + t307 * t340) * t351 + (t302 * t410 + t304 * t339 + t306 * t340) * t350 + (t331 * t410 + t332 * t339 + t333 * t340) * t361) / 0.2e1;
T  = t1;
