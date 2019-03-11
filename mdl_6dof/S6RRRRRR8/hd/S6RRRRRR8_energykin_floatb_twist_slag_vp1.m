% Calculate kinetic energy for
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
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
% Datum: 2019-03-10 05:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRR8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRR8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRR8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 04:48:12
% EndTime: 2019-03-10 04:48:16
% DurationCPUTime: 4.40s
% Computational Cost: add. (5653->447), mult. (13362->675), div. (0->0), fcn. (16964->16), ass. (0->199)
t423 = cos(qJ(3));
t373 = cos(pkin(6));
t422 = pkin(9) * t373;
t380 = cos(qJ(4));
t421 = pkin(4) * t380;
t419 = cos(pkin(7));
t418 = sin(pkin(7));
t378 = sin(qJ(1));
t417 = Icges(2,4) * t378;
t381 = cos(qJ(2));
t382 = cos(qJ(1));
t407 = t381 * t382;
t377 = sin(qJ(2));
t410 = t377 * t378;
t340 = t373 * t407 - t410;
t372 = sin(pkin(6));
t400 = t372 * t419;
t323 = -t340 * t418 - t382 * t400;
t375 = sin(qJ(4));
t416 = t323 * t375;
t408 = t378 * t381;
t409 = t377 * t382;
t342 = -t373 * t408 - t409;
t324 = -t342 * t418 + t378 * t400;
t415 = t324 * t375;
t399 = t372 * t418;
t339 = t373 * t419 - t381 * t399;
t414 = t339 * t375;
t413 = t372 * t377;
t412 = t372 * t378;
t411 = t372 * t382;
t406 = qJ(4) + qJ(5);
t405 = qJD(2) * t372;
t404 = V_base(5) * pkin(8) + V_base(1);
t353 = t378 * t405 + V_base(4);
t401 = cos(t406);
t368 = V_base(6) + qJD(1);
t313 = qJD(3) * t324 + t353;
t354 = qJD(2) * t373 + t368;
t398 = t419 * t423;
t397 = t423 * t418;
t343 = -t373 * t410 + t407;
t376 = sin(qJ(3));
t396 = t372 * t397;
t300 = -t342 * t398 + t343 * t376 - t378 * t396;
t270 = qJD(4) * t300 + t313;
t325 = qJD(3) * t339 + t354;
t352 = -t382 * t405 + V_base(5);
t248 = qJD(5) * t300 + t270;
t345 = pkin(1) * t378 - pkin(9) * t411;
t395 = -t345 * t368 + V_base(5) * t422 + t404;
t321 = -t372 * t381 * t398 - t373 * t397 + t376 * t413;
t287 = qJD(4) * t321 + t325;
t346 = pkin(1) * t382 + pkin(9) * t412;
t394 = V_base(4) * t345 - t346 * V_base(5) + V_base(3);
t312 = qJD(3) * t323 + t352;
t273 = qJD(5) * t321 + t287;
t341 = t373 * t409 + t408;
t298 = -t340 * t398 + t341 * t376 + t382 * t396;
t269 = qJD(4) * t298 + t312;
t393 = t368 * t346 + V_base(2) + (-pkin(8) - t422) * V_base(4);
t247 = qJD(5) * t298 + t269;
t302 = t341 * pkin(2) + pkin(10) * t323;
t329 = pkin(2) * t413 + pkin(10) * t339;
t392 = -t302 * t354 + t352 * t329 + t395;
t303 = t343 * pkin(2) + pkin(10) * t324;
t391 = t353 * t302 - t303 * t352 + t394;
t390 = t354 * t303 - t329 * t353 + t393;
t299 = t341 * t423 + (t340 * t419 - t382 * t399) * t376;
t266 = t299 * pkin(3) + t298 * pkin(11);
t322 = t373 * t418 * t376 + (t419 * t376 * t381 + t377 * t423) * t372;
t286 = t322 * pkin(3) + t321 * pkin(11);
t389 = -t266 * t325 + t312 * t286 + t392;
t301 = t343 * t423 + (t342 * t419 + t378 * t399) * t376;
t267 = t301 * pkin(3) + t300 * pkin(11);
t388 = t313 * t266 - t267 * t312 + t391;
t387 = t325 * t267 - t286 * t313 + t390;
t209 = pkin(4) * t416 + pkin(12) * t298 + t299 * t421;
t236 = pkin(4) * t414 + pkin(12) * t321 + t322 * t421;
t386 = -t209 * t287 + t269 * t236 + t389;
t210 = pkin(4) * t415 + pkin(12) * t300 + t301 * t421;
t385 = t270 * t209 - t210 * t269 + t388;
t384 = t287 * t210 - t236 * t270 + t387;
t379 = cos(qJ(6));
t374 = sin(qJ(6));
t370 = Icges(2,4) * t382;
t369 = sin(t406);
t362 = rSges(2,1) * t382 - rSges(2,2) * t378;
t361 = rSges(2,1) * t378 + rSges(2,2) * t382;
t360 = Icges(2,1) * t382 - t417;
t359 = Icges(2,1) * t378 + t370;
t358 = -Icges(2,2) * t378 + t370;
t357 = Icges(2,2) * t382 + t417;
t351 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t350 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t349 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t336 = rSges(3,3) * t373 + (rSges(3,1) * t377 + rSges(3,2) * t381) * t372;
t335 = Icges(3,5) * t373 + (Icges(3,1) * t377 + Icges(3,4) * t381) * t372;
t334 = Icges(3,6) * t373 + (Icges(3,4) * t377 + Icges(3,2) * t381) * t372;
t333 = Icges(3,3) * t373 + (Icges(3,5) * t377 + Icges(3,6) * t381) * t372;
t328 = V_base(5) * rSges(2,3) - t361 * t368 + t404;
t327 = t362 * t368 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t326 = t361 * V_base(4) - t362 * V_base(5) + V_base(3);
t311 = rSges(3,1) * t343 + rSges(3,2) * t342 + rSges(3,3) * t412;
t310 = rSges(3,1) * t341 + rSges(3,2) * t340 - rSges(3,3) * t411;
t309 = Icges(3,1) * t343 + Icges(3,4) * t342 + Icges(3,5) * t412;
t308 = Icges(3,1) * t341 + Icges(3,4) * t340 - Icges(3,5) * t411;
t307 = Icges(3,4) * t343 + Icges(3,2) * t342 + Icges(3,6) * t412;
t306 = Icges(3,4) * t341 + Icges(3,2) * t340 - Icges(3,6) * t411;
t305 = Icges(3,5) * t343 + Icges(3,6) * t342 + Icges(3,3) * t412;
t304 = Icges(3,5) * t341 + Icges(3,6) * t340 - Icges(3,3) * t411;
t297 = t322 * t380 + t414;
t296 = -t322 * t375 + t339 * t380;
t290 = t322 * t401 + t339 * t369;
t289 = t322 * t369 - t339 * t401;
t285 = rSges(4,1) * t322 - rSges(4,2) * t321 + rSges(4,3) * t339;
t284 = Icges(4,1) * t322 - Icges(4,4) * t321 + Icges(4,5) * t339;
t283 = Icges(4,4) * t322 - Icges(4,2) * t321 + Icges(4,6) * t339;
t282 = Icges(4,5) * t322 - Icges(4,6) * t321 + Icges(4,3) * t339;
t281 = t301 * t380 + t415;
t280 = -t301 * t375 + t324 * t380;
t279 = t299 * t380 + t416;
t278 = -t299 * t375 + t323 * t380;
t277 = t301 * t401 + t324 * t369;
t276 = t301 * t369 - t324 * t401;
t275 = t299 * t401 + t323 * t369;
t274 = t299 * t369 - t323 * t401;
t272 = t290 * t379 + t321 * t374;
t271 = -t290 * t374 + t321 * t379;
t265 = -t310 * t354 + t336 * t352 + t395;
t264 = t311 * t354 - t336 * t353 + t393;
t263 = pkin(5) * t290 + pkin(13) * t289;
t261 = rSges(4,1) * t301 - rSges(4,2) * t300 + rSges(4,3) * t324;
t260 = rSges(4,1) * t299 - rSges(4,2) * t298 + rSges(4,3) * t323;
t259 = t310 * t353 - t311 * t352 + t394;
t258 = Icges(4,1) * t301 - Icges(4,4) * t300 + Icges(4,5) * t324;
t257 = Icges(4,1) * t299 - Icges(4,4) * t298 + Icges(4,5) * t323;
t256 = Icges(4,4) * t301 - Icges(4,2) * t300 + Icges(4,6) * t324;
t255 = Icges(4,4) * t299 - Icges(4,2) * t298 + Icges(4,6) * t323;
t254 = Icges(4,5) * t301 - Icges(4,6) * t300 + Icges(4,3) * t324;
t253 = Icges(4,5) * t299 - Icges(4,6) * t298 + Icges(4,3) * t323;
t252 = rSges(5,1) * t297 + rSges(5,2) * t296 + rSges(5,3) * t321;
t251 = Icges(5,1) * t297 + Icges(5,4) * t296 + Icges(5,5) * t321;
t250 = Icges(5,4) * t297 + Icges(5,2) * t296 + Icges(5,6) * t321;
t249 = Icges(5,5) * t297 + Icges(5,6) * t296 + Icges(5,3) * t321;
t245 = rSges(6,1) * t290 - rSges(6,2) * t289 + rSges(6,3) * t321;
t244 = t277 * t379 + t300 * t374;
t243 = -t277 * t374 + t300 * t379;
t242 = t275 * t379 + t298 * t374;
t241 = -t275 * t374 + t298 * t379;
t240 = Icges(6,1) * t290 - Icges(6,4) * t289 + Icges(6,5) * t321;
t239 = Icges(6,4) * t290 - Icges(6,2) * t289 + Icges(6,6) * t321;
t238 = Icges(6,5) * t290 - Icges(6,6) * t289 + Icges(6,3) * t321;
t237 = qJD(6) * t289 + t273;
t235 = pkin(5) * t277 + pkin(13) * t276;
t234 = pkin(5) * t275 + pkin(13) * t274;
t233 = rSges(5,1) * t281 + rSges(5,2) * t280 + rSges(5,3) * t300;
t232 = rSges(5,1) * t279 + rSges(5,2) * t278 + rSges(5,3) * t298;
t231 = Icges(5,1) * t281 + Icges(5,4) * t280 + Icges(5,5) * t300;
t230 = Icges(5,1) * t279 + Icges(5,4) * t278 + Icges(5,5) * t298;
t229 = Icges(5,4) * t281 + Icges(5,2) * t280 + Icges(5,6) * t300;
t228 = Icges(5,4) * t279 + Icges(5,2) * t278 + Icges(5,6) * t298;
t227 = Icges(5,5) * t281 + Icges(5,6) * t280 + Icges(5,3) * t300;
t226 = Icges(5,5) * t279 + Icges(5,6) * t278 + Icges(5,3) * t298;
t225 = rSges(6,1) * t277 - rSges(6,2) * t276 + rSges(6,3) * t300;
t224 = rSges(6,1) * t275 - rSges(6,2) * t274 + rSges(6,3) * t298;
t223 = Icges(6,1) * t277 - Icges(6,4) * t276 + Icges(6,5) * t300;
t222 = Icges(6,1) * t275 - Icges(6,4) * t274 + Icges(6,5) * t298;
t221 = Icges(6,4) * t277 - Icges(6,2) * t276 + Icges(6,6) * t300;
t220 = Icges(6,4) * t275 - Icges(6,2) * t274 + Icges(6,6) * t298;
t219 = Icges(6,5) * t277 - Icges(6,6) * t276 + Icges(6,3) * t300;
t218 = Icges(6,5) * t275 - Icges(6,6) * t274 + Icges(6,3) * t298;
t217 = qJD(6) * t276 + t248;
t216 = qJD(6) * t274 + t247;
t214 = rSges(7,1) * t272 + rSges(7,2) * t271 + rSges(7,3) * t289;
t213 = Icges(7,1) * t272 + Icges(7,4) * t271 + Icges(7,5) * t289;
t212 = Icges(7,4) * t272 + Icges(7,2) * t271 + Icges(7,6) * t289;
t211 = Icges(7,5) * t272 + Icges(7,6) * t271 + Icges(7,3) * t289;
t206 = rSges(7,1) * t244 + rSges(7,2) * t243 + rSges(7,3) * t276;
t205 = rSges(7,1) * t242 + rSges(7,2) * t241 + rSges(7,3) * t274;
t204 = Icges(7,1) * t244 + Icges(7,4) * t243 + Icges(7,5) * t276;
t203 = Icges(7,1) * t242 + Icges(7,4) * t241 + Icges(7,5) * t274;
t202 = Icges(7,4) * t244 + Icges(7,2) * t243 + Icges(7,6) * t276;
t201 = Icges(7,4) * t242 + Icges(7,2) * t241 + Icges(7,6) * t274;
t200 = Icges(7,5) * t244 + Icges(7,6) * t243 + Icges(7,3) * t276;
t199 = Icges(7,5) * t242 + Icges(7,6) * t241 + Icges(7,3) * t274;
t198 = -t260 * t325 + t285 * t312 + t392;
t197 = t261 * t325 - t285 * t313 + t390;
t196 = t260 * t313 - t261 * t312 + t391;
t195 = -t232 * t287 + t252 * t269 + t389;
t194 = t233 * t287 - t252 * t270 + t387;
t193 = t232 * t270 - t233 * t269 + t388;
t192 = -t224 * t273 + t245 * t247 + t386;
t191 = t225 * t273 - t245 * t248 + t384;
t190 = t224 * t248 - t225 * t247 + t385;
t189 = -t205 * t237 + t214 * t216 - t234 * t273 + t247 * t263 + t386;
t188 = t206 * t237 - t214 * t217 + t235 * t273 - t248 * t263 + t384;
t187 = t205 * t217 - t206 * t216 + t234 * t248 - t235 * t247 + t385;
t1 = t354 * ((t304 * t352 + t305 * t353 + t333 * t354) * t373 + ((t307 * t381 + t309 * t377) * t353 + (t306 * t381 + t308 * t377) * t352 + (t334 * t381 + t335 * t377) * t354) * t372) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + m(1) * (t349 ^ 2 + t350 ^ 2 + t351 ^ 2) / 0.2e1 + t325 * ((t254 * t339 - t256 * t321 + t258 * t322) * t313 + (t253 * t339 - t255 * t321 + t257 * t322) * t312 + (t282 * t339 - t283 * t321 + t284 * t322) * t325) / 0.2e1 + t312 * ((t254 * t323 - t256 * t298 + t258 * t299) * t313 + (t253 * t323 - t255 * t298 + t257 * t299) * t312 + (t282 * t323 - t283 * t298 + t284 * t299) * t325) / 0.2e1 + t313 * ((t254 * t324 - t256 * t300 + t258 * t301) * t313 + (t253 * t324 - t255 * t300 + t257 * t301) * t312 + (t282 * t324 - t283 * t300 + t284 * t301) * t325) / 0.2e1 + m(2) * (t326 ^ 2 + t327 ^ 2 + t328 ^ 2) / 0.2e1 + t287 * ((t227 * t321 + t229 * t296 + t231 * t297) * t270 + (t226 * t321 + t228 * t296 + t230 * t297) * t269 + (t249 * t321 + t250 * t296 + t251 * t297) * t287) / 0.2e1 + t273 * ((t219 * t321 - t221 * t289 + t223 * t290) * t248 + (t218 * t321 - t220 * t289 + t222 * t290) * t247 + (t321 * t238 - t289 * t239 + t290 * t240) * t273) / 0.2e1 + ((-t357 * t378 + t359 * t382 + Icges(1,4)) * V_base(5) + (-t358 * t378 + t360 * t382 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + t270 * ((t300 * t227 + t280 * t229 + t281 * t231) * t270 + (t226 * t300 + t228 * t280 + t230 * t281) * t269 + (t249 * t300 + t250 * t280 + t251 * t281) * t287) / 0.2e1 + t248 * ((t300 * t219 - t276 * t221 + t277 * t223) * t248 + (t218 * t300 - t220 * t276 + t222 * t277) * t247 + (t238 * t300 - t239 * t276 + t240 * t277) * t273) / 0.2e1 + t269 * ((t227 * t298 + t229 * t278 + t231 * t279) * t270 + (t298 * t226 + t278 * t228 + t279 * t230) * t269 + (t249 * t298 + t250 * t278 + t251 * t279) * t287) / 0.2e1 + t247 * ((t219 * t298 - t221 * t274 + t223 * t275) * t248 + (t298 * t218 - t274 * t220 + t275 * t222) * t247 + (t238 * t298 - t239 * t274 + t240 * t275) * t273) / 0.2e1 + t237 * ((t200 * t289 + t202 * t271 + t204 * t272) * t217 + (t199 * t289 + t201 * t271 + t203 * t272) * t216 + (t289 * t211 + t271 * t212 + t272 * t213) * t237) / 0.2e1 + t217 * ((t276 * t200 + t243 * t202 + t244 * t204) * t217 + (t199 * t276 + t201 * t243 + t203 * t244) * t216 + (t211 * t276 + t212 * t243 + t213 * t244) * t237) / 0.2e1 + t216 * ((t200 * t274 + t202 * t241 + t204 * t242) * t217 + (t274 * t199 + t241 * t201 + t242 * t203) * t216 + (t211 * t274 + t212 * t241 + t213 * t242) * t237) / 0.2e1 + m(3) * (t259 ^ 2 + t264 ^ 2 + t265 ^ 2) / 0.2e1 + ((t357 * t382 + t359 * t378 + Icges(1,2)) * V_base(5) + (t358 * t382 + t360 * t378 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + t352 * ((-t305 * t411 + t307 * t340 + t309 * t341) * t353 + (-t304 * t411 + t306 * t340 + t308 * t341) * t352 + (-t333 * t411 + t334 * t340 + t335 * t341) * t354) / 0.2e1 + t353 * ((t305 * t412 + t307 * t342 + t309 * t343) * t353 + (t304 * t412 + t306 * t342 + t308 * t343) * t352 + (t333 * t412 + t334 * t342 + t335 * t343) * t354) / 0.2e1 + m(7) * (t187 ^ 2 + t188 ^ 2 + t189 ^ 2) / 0.2e1 + m(6) * (t190 ^ 2 + t191 ^ 2 + t192 ^ 2) / 0.2e1 + ((Icges(2,5) * t378 + Icges(2,6) * t382) * V_base(5) + (Icges(2,5) * t382 - Icges(2,6) * t378) * V_base(4) + Icges(2,3) * t368 / 0.2e1) * t368 + m(5) * (t193 ^ 2 + t194 ^ 2 + t195 ^ 2) / 0.2e1 + m(4) * (t196 ^ 2 + t197 ^ 2 + t198 ^ 2) / 0.2e1;
T  = t1;
