% Calculate kinetic energy for
% S6RPRRPR11
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
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPR11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR11_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR11_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:38:59
% EndTime: 2019-03-09 05:39:05
% DurationCPUTime: 5.29s
% Computational Cost: add. (5568->447), mult. (14201->633), div. (0->0), fcn. (18231->16), ass. (0->193)
t412 = Icges(5,2) + Icges(6,3);
t356 = cos(pkin(12));
t353 = sin(pkin(12));
t361 = sin(qJ(1));
t390 = t361 * t353;
t357 = cos(pkin(6));
t362 = cos(qJ(1));
t391 = t357 * t362;
t321 = t356 * t391 - t390;
t389 = t361 * t356;
t322 = t353 * t391 + t389;
t360 = sin(qJ(3));
t354 = sin(pkin(6));
t400 = sin(pkin(7));
t381 = t354 * t400;
t401 = cos(pkin(7));
t404 = cos(qJ(3));
t283 = t322 * t404 + (t321 * t401 - t362 * t381) * t360;
t359 = sin(qJ(4));
t382 = t354 * t401;
t374 = -t321 * t400 - t362 * t382;
t403 = cos(qJ(4));
t267 = t283 * t403 + t359 * t374;
t378 = t404 * t400;
t376 = t354 * t378;
t379 = t401 * t404;
t282 = -t321 * t379 + t322 * t360 + t362 * t376;
t352 = sin(pkin(13));
t355 = cos(pkin(13));
t232 = -t267 * t352 + t282 * t355;
t397 = t282 * t352;
t233 = t267 * t355 + t397;
t266 = t283 * t359 - t374 * t403;
t411 = -Icges(5,4) * t267 + Icges(6,5) * t233 - Icges(5,6) * t282 + Icges(6,6) * t232 + t412 * t266;
t323 = -t353 * t362 - t357 * t389;
t324 = t356 * t362 - t357 * t390;
t285 = t324 * t404 + (t323 * t401 + t361 * t381) * t360;
t373 = -t323 * t400 + t361 * t382;
t269 = t285 * t403 + t359 * t373;
t284 = -t323 * t379 + t324 * t360 - t361 * t376;
t234 = -t269 * t352 + t284 * t355;
t396 = t284 * t352;
t235 = t269 * t355 + t396;
t268 = t285 * t359 - t373 * t403;
t410 = -Icges(5,4) * t269 + Icges(6,5) * t235 - Icges(5,6) * t284 + Icges(6,6) * t234 + t412 * t268;
t305 = t357 * t400 * t360 + (t356 * t360 * t401 + t353 * t404) * t354;
t372 = -t356 * t381 + t357 * t401;
t280 = t305 * t403 + t359 * t372;
t394 = t353 * t354;
t304 = -t354 * t356 * t379 - t357 * t378 + t360 * t394;
t260 = -t280 * t352 + t304 * t355;
t395 = t304 * t352;
t261 = t280 * t355 + t395;
t279 = t305 * t359 - t372 * t403;
t409 = -Icges(5,4) * t280 + Icges(6,5) * t261 - Icges(5,6) * t304 + Icges(6,6) * t260 + t412 * t279;
t402 = pkin(5) * t355;
t399 = Icges(2,4) * t361;
t398 = qJ(2) * t357;
t393 = t354 * t361;
t392 = t354 * t362;
t387 = qJD(2) * t354;
t386 = V_base(5) * pkin(8) + V_base(1);
t298 = qJD(3) * t373 + V_base(4);
t297 = qJD(3) * t374 + V_base(5);
t348 = V_base(6) + qJD(1);
t383 = -pkin(8) - t398;
t326 = t361 * pkin(1) - qJ(2) * t392;
t380 = qJD(2) * t357 + V_base(4) * t326 + V_base(3);
t265 = qJD(4) * t284 + t298;
t264 = qJD(4) * t282 + t297;
t311 = qJD(3) * t372 + t348;
t377 = t361 * t387 + V_base(5) * t398 + t386;
t276 = qJD(4) * t304 + t311;
t327 = pkin(1) * t362 + qJ(2) * t393;
t375 = t348 * t327 - t362 * t387 + V_base(2);
t287 = t322 * pkin(2) + pkin(9) * t374;
t308 = pkin(2) * t394 + pkin(9) * t372;
t371 = V_base(5) * t308 + (-t287 - t326) * t348 + t377;
t288 = t324 * pkin(2) + pkin(9) * t373;
t370 = V_base(4) * t287 + (-t288 - t327) * V_base(5) + t380;
t255 = pkin(3) * t283 + pkin(10) * t282;
t274 = pkin(3) * t305 + pkin(10) * t304;
t369 = -t255 * t311 + t297 * t274 + t371;
t256 = pkin(3) * t285 + pkin(10) * t284;
t368 = t298 * t255 - t256 * t297 + t370;
t367 = t348 * t288 + (-t308 + t383) * V_base(4) + t375;
t251 = pkin(4) * t280 + qJ(5) * t279;
t366 = qJD(5) * t268 + t264 * t251 + t369;
t224 = pkin(4) * t267 + qJ(5) * t266;
t365 = qJD(5) * t279 + t265 * t224 + t368;
t364 = t311 * t256 - t298 * t274 + t367;
t225 = pkin(4) * t269 + qJ(5) * t268;
t363 = qJD(5) * t266 + t276 * t225 + t364;
t351 = pkin(13) + qJ(6);
t349 = Icges(2,4) * t362;
t347 = cos(t351);
t346 = sin(t351);
t340 = rSges(2,1) * t362 - t361 * rSges(2,2);
t339 = t361 * rSges(2,1) + rSges(2,2) * t362;
t338 = Icges(2,1) * t362 - t399;
t337 = Icges(2,1) * t361 + t349;
t336 = -Icges(2,2) * t361 + t349;
t335 = Icges(2,2) * t362 + t399;
t334 = Icges(2,5) * t362 - Icges(2,6) * t361;
t333 = Icges(2,5) * t361 + Icges(2,6) * t362;
t332 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t331 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t330 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t317 = rSges(3,3) * t357 + (rSges(3,1) * t353 + rSges(3,2) * t356) * t354;
t316 = Icges(3,5) * t357 + (Icges(3,1) * t353 + Icges(3,4) * t356) * t354;
t315 = Icges(3,6) * t357 + (Icges(3,4) * t353 + Icges(3,2) * t356) * t354;
t314 = Icges(3,3) * t357 + (Icges(3,5) * t353 + Icges(3,6) * t356) * t354;
t310 = V_base(5) * rSges(2,3) - t339 * t348 + t386;
t309 = t340 * t348 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t307 = t339 * V_base(4) - t340 * V_base(5) + V_base(3);
t296 = rSges(3,1) * t324 + rSges(3,2) * t323 + rSges(3,3) * t393;
t295 = t322 * rSges(3,1) + t321 * rSges(3,2) - rSges(3,3) * t392;
t294 = Icges(3,1) * t324 + Icges(3,4) * t323 + Icges(3,5) * t393;
t293 = Icges(3,1) * t322 + Icges(3,4) * t321 - Icges(3,5) * t392;
t292 = Icges(3,4) * t324 + Icges(3,2) * t323 + Icges(3,6) * t393;
t291 = Icges(3,4) * t322 + Icges(3,2) * t321 - Icges(3,6) * t392;
t290 = Icges(3,5) * t324 + Icges(3,6) * t323 + Icges(3,3) * t393;
t289 = Icges(3,5) * t322 + Icges(3,6) * t321 - Icges(3,3) * t392;
t273 = t305 * rSges(4,1) - t304 * rSges(4,2) + rSges(4,3) * t372;
t272 = Icges(4,1) * t305 - Icges(4,4) * t304 + Icges(4,5) * t372;
t271 = Icges(4,4) * t305 - Icges(4,2) * t304 + Icges(4,6) * t372;
t270 = Icges(4,5) * t305 - Icges(4,6) * t304 + Icges(4,3) * t372;
t259 = t280 * t347 + t304 * t346;
t258 = -t280 * t346 + t304 * t347;
t254 = t317 * V_base(5) + (-t295 - t326) * t348 + t377;
t253 = t348 * t296 + (-t317 + t383) * V_base(4) + t375;
t252 = qJD(6) * t279 + t276;
t250 = t295 * V_base(4) + (-t296 - t327) * V_base(5) + t380;
t248 = t285 * rSges(4,1) - t284 * rSges(4,2) + rSges(4,3) * t373;
t247 = t283 * rSges(4,1) - t282 * rSges(4,2) + rSges(4,3) * t374;
t246 = Icges(4,1) * t285 - Icges(4,4) * t284 + Icges(4,5) * t373;
t245 = Icges(4,1) * t283 - Icges(4,4) * t282 + Icges(4,5) * t374;
t244 = Icges(4,4) * t285 - Icges(4,2) * t284 + Icges(4,6) * t373;
t243 = Icges(4,4) * t283 - Icges(4,2) * t282 + Icges(4,6) * t374;
t242 = Icges(4,5) * t285 - Icges(4,6) * t284 + Icges(4,3) * t373;
t241 = Icges(4,5) * t283 - Icges(4,6) * t282 + Icges(4,3) * t374;
t240 = rSges(5,1) * t280 - rSges(5,2) * t279 + rSges(5,3) * t304;
t238 = Icges(5,1) * t280 - Icges(5,4) * t279 + Icges(5,5) * t304;
t236 = Icges(5,5) * t280 - Icges(5,6) * t279 + Icges(5,3) * t304;
t231 = t269 * t347 + t284 * t346;
t230 = -t269 * t346 + t284 * t347;
t229 = t267 * t347 + t282 * t346;
t228 = -t267 * t346 + t282 * t347;
t227 = qJD(6) * t268 + t265;
t226 = qJD(6) * t266 + t264;
t222 = rSges(5,1) * t269 - rSges(5,2) * t268 + rSges(5,3) * t284;
t221 = rSges(5,1) * t267 - rSges(5,2) * t266 + rSges(5,3) * t282;
t220 = Icges(5,1) * t269 - Icges(5,4) * t268 + Icges(5,5) * t284;
t219 = Icges(5,1) * t267 - Icges(5,4) * t266 + Icges(5,5) * t282;
t216 = Icges(5,5) * t269 - Icges(5,6) * t268 + Icges(5,3) * t284;
t215 = Icges(5,5) * t267 - Icges(5,6) * t266 + Icges(5,3) * t282;
t213 = rSges(6,1) * t261 + rSges(6,2) * t260 + rSges(6,3) * t279;
t212 = Icges(6,1) * t261 + Icges(6,4) * t260 + Icges(6,5) * t279;
t211 = Icges(6,4) * t261 + Icges(6,2) * t260 + Icges(6,6) * t279;
t209 = rSges(7,1) * t259 + rSges(7,2) * t258 + rSges(7,3) * t279;
t208 = Icges(7,1) * t259 + Icges(7,4) * t258 + Icges(7,5) * t279;
t207 = Icges(7,4) * t259 + Icges(7,2) * t258 + Icges(7,6) * t279;
t206 = Icges(7,5) * t259 + Icges(7,6) * t258 + Icges(7,3) * t279;
t205 = pkin(5) * t395 + pkin(11) * t279 + t280 * t402;
t203 = rSges(6,1) * t235 + rSges(6,2) * t234 + rSges(6,3) * t268;
t202 = rSges(6,1) * t233 + rSges(6,2) * t232 + rSges(6,3) * t266;
t201 = Icges(6,1) * t235 + Icges(6,4) * t234 + Icges(6,5) * t268;
t200 = Icges(6,1) * t233 + Icges(6,4) * t232 + Icges(6,5) * t266;
t199 = Icges(6,4) * t235 + Icges(6,2) * t234 + Icges(6,6) * t268;
t198 = Icges(6,4) * t233 + Icges(6,2) * t232 + Icges(6,6) * t266;
t195 = rSges(7,1) * t231 + rSges(7,2) * t230 + rSges(7,3) * t268;
t194 = rSges(7,1) * t229 + rSges(7,2) * t228 + rSges(7,3) * t266;
t193 = Icges(7,1) * t231 + Icges(7,4) * t230 + Icges(7,5) * t268;
t192 = Icges(7,1) * t229 + Icges(7,4) * t228 + Icges(7,5) * t266;
t191 = Icges(7,4) * t231 + Icges(7,2) * t230 + Icges(7,6) * t268;
t190 = Icges(7,4) * t229 + Icges(7,2) * t228 + Icges(7,6) * t266;
t189 = Icges(7,5) * t231 + Icges(7,6) * t230 + Icges(7,3) * t268;
t188 = Icges(7,5) * t229 + Icges(7,6) * t228 + Icges(7,3) * t266;
t187 = -t247 * t311 + t273 * t297 + t371;
t186 = t311 * t248 - t298 * t273 + t367;
t185 = pkin(5) * t396 + pkin(11) * t268 + t269 * t402;
t184 = pkin(5) * t397 + pkin(11) * t266 + t267 * t402;
t183 = t247 * t298 - t248 * t297 + t370;
t182 = -t221 * t276 + t240 * t264 + t369;
t181 = t276 * t222 - t265 * t240 + t364;
t180 = t221 * t265 - t222 * t264 + t368;
t179 = t213 * t264 + (-t202 - t224) * t276 + t366;
t178 = t276 * t203 + (-t213 - t251) * t265 + t363;
t177 = t202 * t265 + (-t203 - t225) * t264 + t365;
t176 = (-t184 - t224) * t276 - t194 * t252 + t205 * t264 + t209 * t226 + t366;
t175 = t276 * t185 + t252 * t195 - t227 * t209 + (-t205 - t251) * t265 + t363;
t174 = t184 * t265 + t194 * t227 - t195 * t226 + (-t185 - t225) * t264 + t365;
t1 = t297 * ((t242 * t374 - t282 * t244 + t283 * t246) * t298 + (t241 * t374 - t282 * t243 + t283 * t245) * t297 + (t270 * t374 - t282 * t271 + t283 * t272) * t311) / 0.2e1 + t311 * ((t242 * t372 - t304 * t244 + t305 * t246) * t298 + (t241 * t372 - t304 * t243 + t305 * t245) * t297 + (t270 * t372 - t304 * t271 + t305 * t272) * t311) / 0.2e1 + t298 * ((t242 * t373 - t284 * t244 + t285 * t246) * t298 + (t241 * t373 - t284 * t243 + t285 * t245) * t297 + (t270 * t373 - t284 * t271 + t285 * t272) * t311) / 0.2e1 + m(1) * (t330 ^ 2 + t331 ^ 2 + t332 ^ 2) / 0.2e1 + m(2) * (t307 ^ 2 + t309 ^ 2 + t310 ^ 2) / 0.2e1 + t252 * ((t189 * t279 + t191 * t258 + t193 * t259) * t227 + (t188 * t279 + t190 * t258 + t192 * t259) * t226 + (t206 * t279 + t207 * t258 + t208 * t259) * t252) / 0.2e1 + t227 * ((t189 * t268 + t191 * t230 + t193 * t231) * t227 + (t188 * t268 + t190 * t230 + t192 * t231) * t226 + (t206 * t268 + t207 * t230 + t208 * t231) * t252) / 0.2e1 + t226 * ((t189 * t266 + t191 * t228 + t193 * t229) * t227 + (t266 * t188 + t228 * t190 + t229 * t192) * t226 + (t206 * t266 + t207 * t228 + t208 * t229) * t252) / 0.2e1 + m(3) * (t250 ^ 2 + t253 ^ 2 + t254 ^ 2) / 0.2e1 + m(4) * (t183 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(5) * (t180 ^ 2 + t181 ^ 2 + t182 ^ 2) / 0.2e1 + m(6) * (t177 ^ 2 + t178 ^ 2 + t179 ^ 2) / 0.2e1 + m(7) * (t174 ^ 2 + t175 ^ 2 + t176 ^ 2) / 0.2e1 + ((t211 * t232 + t212 * t233 + t236 * t282 + t238 * t267 + t266 * t409) * t276 + (t199 * t232 + t201 * t233 + t216 * t282 + t220 * t267 + t266 * t410) * t265 + (t198 * t232 + t200 * t233 + t215 * t282 + t219 * t267 + t411 * t266) * t264) * t264 / 0.2e1 + ((t211 * t234 + t212 * t235 + t236 * t284 + t238 * t269 + t268 * t409) * t276 + (t199 * t234 + t201 * t235 + t216 * t284 + t220 * t269 + t410 * t268) * t265 + (t198 * t234 + t200 * t235 + t215 * t284 + t219 * t269 + t268 * t411) * t264) * t265 / 0.2e1 + ((t211 * t260 + t212 * t261 + t236 * t304 + t238 * t280 + t409 * t279) * t276 + (t199 * t260 + t201 * t261 + t216 * t304 + t220 * t280 + t279 * t410) * t265 + (t198 * t260 + t200 * t261 + t215 * t304 + t219 * t280 + t279 * t411) * t264) * t276 / 0.2e1 + ((t289 * V_base(5) + t290 * V_base(4) + t314 * t348) * t357 + ((t292 * t356 + t294 * t353) * V_base(4) + (t291 * t356 + t293 * t353) * V_base(5) + (t315 * t356 + t316 * t353) * t348) * t354 + Icges(2,3) * t348 + t333 * V_base(5) + t334 * V_base(4)) * t348 / 0.2e1 + ((t314 * t393 + t315 * t323 + t316 * t324 + t334) * t348 + (t289 * t393 + t291 * t323 + t293 * t324 - t361 * t335 + t337 * t362 + Icges(1,4)) * V_base(5) + (t290 * t393 + t292 * t323 + t294 * t324 - t361 * t336 + t338 * t362 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t314 * t392 + t321 * t315 + t322 * t316 + t333) * t348 + (-t289 * t392 + t321 * t291 + t322 * t293 + t335 * t362 + t361 * t337 + Icges(1,2)) * V_base(5) + (-t290 * t392 + t321 * t292 + t322 * t294 + t336 * t362 + t361 * t338 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
