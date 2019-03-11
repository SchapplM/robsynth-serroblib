% Calculate kinetic energy for
% S6PPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-03-08 19:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PPRRRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PPRRRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PPRRRR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PPRRRR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:03:10
% EndTime: 2019-03-08 19:03:14
% DurationCPUTime: 3.71s
% Computational Cost: add. (5760->452), mult. (14861->657), div. (0->0), fcn. (19107->16), ass. (0->197)
t354 = sin(pkin(12));
t357 = cos(pkin(12));
t408 = Icges(2,5) * t357 - Icges(2,6) * t354 + Icges(1,5);
t407 = Icges(2,5) * t354 + Icges(2,6) * t357 + Icges(1,6);
t406 = cos(qJ(3));
t405 = cos(qJ(4));
t362 = cos(qJ(5));
t404 = pkin(5) * t362;
t402 = cos(pkin(7));
t401 = sin(pkin(7));
t400 = Icges(2,4) * t354;
t358 = cos(pkin(6));
t399 = qJ(2) * t358;
t353 = sin(pkin(13));
t356 = cos(pkin(13));
t391 = t357 * t358;
t322 = -t353 * t354 + t356 * t391;
t323 = t353 * t391 + t354 * t356;
t361 = sin(qJ(3));
t355 = sin(pkin(6));
t380 = t406 * t401;
t377 = t355 * t380;
t381 = t402 * t406;
t280 = -t322 * t381 + t323 * t361 + t357 * t377;
t359 = sin(qJ(5));
t398 = t280 * t359;
t393 = t354 * t358;
t324 = -t353 * t357 - t356 * t393;
t325 = -t353 * t393 + t356 * t357;
t282 = -t324 * t381 + t325 * t361 - t354 * t377;
t397 = t282 * t359;
t395 = t353 * t355;
t306 = -t355 * t356 * t381 - t358 * t380 + t361 * t395;
t396 = t306 * t359;
t394 = t354 * t355;
t392 = t355 * t357;
t390 = qJD(2) * t355;
t389 = V_base(5) * qJ(1) + V_base(1);
t385 = qJD(1) + V_base(3);
t384 = t355 * t402;
t374 = -t324 * t401 + t354 * t384;
t300 = qJD(3) * t374 + V_base(4);
t375 = -t322 * t401 - t357 * t384;
t299 = qJD(3) * t375 + V_base(5);
t383 = t355 * t401;
t373 = -t356 * t383 + t358 * t402;
t315 = qJD(3) * t373 + V_base(6);
t382 = -qJ(1) - t399;
t265 = qJD(4) * t282 + t300;
t264 = qJD(4) * t280 + t299;
t284 = qJD(4) * t306 + t315;
t379 = t354 * t390 + V_base(5) * t399 + t389;
t328 = pkin(1) * t354 - qJ(2) * t392;
t378 = qJD(2) * t358 + V_base(4) * t328 + t385;
t283 = t325 * t406 + (t324 * t402 + t354 * t383) * t361;
t360 = sin(qJ(4));
t268 = t283 * t360 - t374 * t405;
t228 = qJD(5) * t268 + t265;
t281 = t323 * t406 + (t322 * t402 - t357 * t383) * t361;
t266 = t281 * t360 - t375 * t405;
t227 = qJD(5) * t266 + t264;
t307 = t358 * t401 * t361 + (t402 * t356 * t361 + t353 * t406) * t355;
t285 = t307 * t360 - t373 * t405;
t256 = qJD(5) * t285 + t284;
t329 = pkin(1) * t357 + qJ(2) * t394;
t376 = V_base(6) * t329 - t357 * t390 + V_base(2);
t289 = t323 * pkin(2) + pkin(8) * t375;
t310 = pkin(2) * t395 + pkin(8) * t373;
t372 = V_base(5) * t310 + (-t289 - t328) * V_base(6) + t379;
t290 = t325 * pkin(2) + pkin(8) * t374;
t371 = V_base(4) * t289 + (-t290 - t329) * V_base(5) + t378;
t253 = pkin(3) * t281 + pkin(9) * t280;
t276 = pkin(3) * t307 + pkin(9) * t306;
t370 = -t253 * t315 + t299 * t276 + t372;
t254 = pkin(3) * t283 + pkin(9) * t282;
t369 = t300 * t253 - t254 * t299 + t371;
t368 = V_base(6) * t290 + (-t310 + t382) * V_base(4) + t376;
t267 = t281 * t405 + t360 * t375;
t225 = t267 * pkin(4) + t266 * pkin(10);
t286 = t307 * t405 + t360 * t373;
t255 = t286 * pkin(4) + t285 * pkin(10);
t367 = -t225 * t284 + t264 * t255 + t370;
t269 = t283 * t405 + t360 * t374;
t226 = t269 * pkin(4) + t268 * pkin(10);
t366 = t265 * t225 - t226 * t264 + t369;
t365 = t315 * t254 - t276 * t300 + t368;
t364 = t284 * t226 - t255 * t265 + t365;
t352 = qJ(5) + qJ(6);
t350 = cos(t352);
t349 = sin(t352);
t348 = Icges(2,4) * t357;
t342 = rSges(2,1) * t357 - rSges(2,2) * t354;
t341 = rSges(2,1) * t354 + rSges(2,2) * t357;
t340 = Icges(2,1) * t357 - t400;
t339 = Icges(2,1) * t354 + t348;
t338 = -Icges(2,2) * t354 + t348;
t337 = Icges(2,2) * t357 + t400;
t334 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t333 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t332 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t319 = rSges(3,3) * t358 + (rSges(3,1) * t353 + rSges(3,2) * t356) * t355;
t318 = Icges(3,5) * t358 + (Icges(3,1) * t353 + Icges(3,4) * t356) * t355;
t317 = Icges(3,6) * t358 + (Icges(3,4) * t353 + Icges(3,2) * t356) * t355;
t316 = Icges(3,3) * t358 + (Icges(3,5) * t353 + Icges(3,6) * t356) * t355;
t312 = V_base(5) * rSges(2,3) - t341 * V_base(6) + t389;
t311 = t342 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t308 = t341 * V_base(4) - t342 * V_base(5) + t385;
t298 = rSges(3,1) * t325 + rSges(3,2) * t324 + rSges(3,3) * t394;
t297 = rSges(3,1) * t323 + rSges(3,2) * t322 - rSges(3,3) * t392;
t296 = Icges(3,1) * t325 + Icges(3,4) * t324 + Icges(3,5) * t394;
t295 = Icges(3,1) * t323 + Icges(3,4) * t322 - Icges(3,5) * t392;
t294 = Icges(3,4) * t325 + Icges(3,2) * t324 + Icges(3,6) * t394;
t293 = Icges(3,4) * t323 + Icges(3,2) * t322 - Icges(3,6) * t392;
t292 = Icges(3,5) * t325 + Icges(3,6) * t324 + Icges(3,3) * t394;
t291 = Icges(3,5) * t323 + Icges(3,6) * t322 - Icges(3,3) * t392;
t275 = t307 * rSges(4,1) - t306 * rSges(4,2) + rSges(4,3) * t373;
t274 = Icges(4,1) * t307 - Icges(4,4) * t306 + Icges(4,5) * t373;
t273 = Icges(4,4) * t307 - Icges(4,2) * t306 + Icges(4,6) * t373;
t272 = Icges(4,5) * t307 - Icges(4,6) * t306 + Icges(4,3) * t373;
t271 = t286 * t362 + t396;
t270 = -t286 * t359 + t306 * t362;
t261 = t286 * t350 + t306 * t349;
t260 = -t286 * t349 + t306 * t350;
t258 = t319 * V_base(5) + (-t297 - t328) * V_base(6) + t379;
t257 = t298 * V_base(6) + (-t319 + t382) * V_base(4) + t376;
t251 = t297 * V_base(4) + (-t298 - t329) * V_base(5) + t378;
t250 = rSges(5,1) * t286 - rSges(5,2) * t285 + rSges(5,3) * t306;
t249 = Icges(5,1) * t286 - Icges(5,4) * t285 + Icges(5,5) * t306;
t248 = Icges(5,4) * t286 - Icges(5,2) * t285 + Icges(5,6) * t306;
t247 = Icges(5,5) * t286 - Icges(5,6) * t285 + Icges(5,3) * t306;
t246 = t283 * rSges(4,1) - t282 * rSges(4,2) + rSges(4,3) * t374;
t245 = t281 * rSges(4,1) - t280 * rSges(4,2) + rSges(4,3) * t375;
t244 = Icges(4,1) * t283 - Icges(4,4) * t282 + Icges(4,5) * t374;
t243 = Icges(4,1) * t281 - Icges(4,4) * t280 + Icges(4,5) * t375;
t242 = Icges(4,4) * t283 - Icges(4,2) * t282 + Icges(4,6) * t374;
t241 = Icges(4,4) * t281 - Icges(4,2) * t280 + Icges(4,6) * t375;
t240 = Icges(4,5) * t283 - Icges(4,6) * t282 + Icges(4,3) * t374;
t239 = Icges(4,5) * t281 - Icges(4,6) * t280 + Icges(4,3) * t375;
t237 = t269 * t362 + t397;
t236 = -t269 * t359 + t282 * t362;
t235 = t267 * t362 + t398;
t234 = -t267 * t359 + t280 * t362;
t233 = t269 * t350 + t282 * t349;
t232 = -t269 * t349 + t282 * t350;
t231 = t267 * t350 + t280 * t349;
t230 = -t267 * t349 + t280 * t350;
t229 = qJD(6) * t285 + t256;
t222 = rSges(6,1) * t271 + rSges(6,2) * t270 + rSges(6,3) * t285;
t221 = Icges(6,1) * t271 + Icges(6,4) * t270 + Icges(6,5) * t285;
t220 = Icges(6,4) * t271 + Icges(6,2) * t270 + Icges(6,6) * t285;
t219 = Icges(6,5) * t271 + Icges(6,6) * t270 + Icges(6,3) * t285;
t218 = rSges(5,1) * t269 - rSges(5,2) * t268 + rSges(5,3) * t282;
t217 = rSges(5,1) * t267 - rSges(5,2) * t266 + rSges(5,3) * t280;
t216 = Icges(5,1) * t269 - Icges(5,4) * t268 + Icges(5,5) * t282;
t215 = Icges(5,1) * t267 - Icges(5,4) * t266 + Icges(5,5) * t280;
t214 = Icges(5,4) * t269 - Icges(5,2) * t268 + Icges(5,6) * t282;
t213 = Icges(5,4) * t267 - Icges(5,2) * t266 + Icges(5,6) * t280;
t212 = Icges(5,5) * t269 - Icges(5,6) * t268 + Icges(5,3) * t282;
t211 = Icges(5,5) * t267 - Icges(5,6) * t266 + Icges(5,3) * t280;
t210 = rSges(7,1) * t261 + rSges(7,2) * t260 + rSges(7,3) * t285;
t209 = Icges(7,1) * t261 + Icges(7,4) * t260 + Icges(7,5) * t285;
t208 = Icges(7,4) * t261 + Icges(7,2) * t260 + Icges(7,6) * t285;
t207 = Icges(7,5) * t261 + Icges(7,6) * t260 + Icges(7,3) * t285;
t206 = pkin(5) * t396 + pkin(11) * t285 + t286 * t404;
t205 = qJD(6) * t268 + t228;
t204 = qJD(6) * t266 + t227;
t202 = rSges(6,1) * t237 + rSges(6,2) * t236 + rSges(6,3) * t268;
t201 = rSges(6,1) * t235 + rSges(6,2) * t234 + rSges(6,3) * t266;
t200 = Icges(6,1) * t237 + Icges(6,4) * t236 + Icges(6,5) * t268;
t199 = Icges(6,1) * t235 + Icges(6,4) * t234 + Icges(6,5) * t266;
t198 = Icges(6,4) * t237 + Icges(6,2) * t236 + Icges(6,6) * t268;
t197 = Icges(6,4) * t235 + Icges(6,2) * t234 + Icges(6,6) * t266;
t196 = Icges(6,5) * t237 + Icges(6,6) * t236 + Icges(6,3) * t268;
t195 = Icges(6,5) * t235 + Icges(6,6) * t234 + Icges(6,3) * t266;
t194 = -t245 * t315 + t275 * t299 + t372;
t193 = t246 * t315 - t275 * t300 + t368;
t192 = rSges(7,1) * t233 + rSges(7,2) * t232 + rSges(7,3) * t268;
t191 = rSges(7,1) * t231 + rSges(7,2) * t230 + rSges(7,3) * t266;
t190 = Icges(7,1) * t233 + Icges(7,4) * t232 + Icges(7,5) * t268;
t189 = Icges(7,1) * t231 + Icges(7,4) * t230 + Icges(7,5) * t266;
t188 = Icges(7,4) * t233 + Icges(7,2) * t232 + Icges(7,6) * t268;
t187 = Icges(7,4) * t231 + Icges(7,2) * t230 + Icges(7,6) * t266;
t186 = Icges(7,5) * t233 + Icges(7,6) * t232 + Icges(7,3) * t268;
t185 = Icges(7,5) * t231 + Icges(7,6) * t230 + Icges(7,3) * t266;
t184 = pkin(5) * t397 + pkin(11) * t268 + t269 * t404;
t183 = pkin(5) * t398 + pkin(11) * t266 + t267 * t404;
t182 = t245 * t300 - t246 * t299 + t371;
t181 = -t217 * t284 + t250 * t264 + t370;
t180 = t218 * t284 - t250 * t265 + t365;
t179 = t217 * t265 - t218 * t264 + t369;
t178 = -t201 * t256 + t222 * t227 + t367;
t177 = t202 * t256 - t222 * t228 + t364;
t176 = t201 * t228 - t202 * t227 + t366;
t175 = -t183 * t256 - t191 * t229 + t204 * t210 + t206 * t227 + t367;
t174 = t184 * t256 + t192 * t229 - t205 * t210 - t206 * t228 + t364;
t173 = t183 * t228 - t184 * t227 + t191 * t205 - t192 * t204 + t366;
t1 = m(1) * (t332 ^ 2 + t333 ^ 2 + t334 ^ 2) / 0.2e1 + t315 * ((t240 * t373 - t306 * t242 + t307 * t244) * t300 + (t239 * t373 - t306 * t241 + t307 * t243) * t299 + (t272 * t373 - t306 * t273 + t307 * t274) * t315) / 0.2e1 + t300 * ((t240 * t374 - t282 * t242 + t283 * t244) * t300 + (t239 * t374 - t282 * t241 + t283 * t243) * t299 + (t272 * t374 - t282 * t273 + t283 * t274) * t315) / 0.2e1 + t299 * ((t240 * t375 - t280 * t242 + t281 * t244) * t300 + (t239 * t375 - t280 * t241 + t281 * t243) * t299 + (t272 * t375 - t280 * t273 + t281 * t274) * t315) / 0.2e1 + t284 * ((t212 * t306 - t214 * t285 + t216 * t286) * t265 + (t211 * t306 - t213 * t285 + t215 * t286) * t264 + (t247 * t306 - t248 * t285 + t249 * t286) * t284) / 0.2e1 + m(2) * (t308 ^ 2 + t311 ^ 2 + t312 ^ 2) / 0.2e1 + t264 * ((t212 * t280 - t214 * t266 + t216 * t267) * t265 + (t211 * t280 - t213 * t266 + t215 * t267) * t264 + (t247 * t280 - t248 * t266 + t249 * t267) * t284) / 0.2e1 + t265 * ((t212 * t282 - t214 * t268 + t216 * t269) * t265 + (t211 * t282 - t213 * t268 + t215 * t269) * t264 + (t247 * t282 - t248 * t268 + t249 * t269) * t284) / 0.2e1 + t256 * ((t196 * t285 + t198 * t270 + t200 * t271) * t228 + (t195 * t285 + t197 * t270 + t199 * t271) * t227 + (t219 * t285 + t220 * t270 + t221 * t271) * t256) / 0.2e1 + t229 * ((t186 * t285 + t188 * t260 + t190 * t261) * t205 + (t185 * t285 + t187 * t260 + t189 * t261) * t204 + (t207 * t285 + t208 * t260 + t209 * t261) * t229) / 0.2e1 + t228 * ((t196 * t268 + t198 * t236 + t200 * t237) * t228 + (t195 * t268 + t197 * t236 + t199 * t237) * t227 + (t219 * t268 + t220 * t236 + t221 * t237) * t256) / 0.2e1 + t205 * ((t268 * t186 + t232 * t188 + t233 * t190) * t205 + (t185 * t268 + t187 * t232 + t189 * t233) * t204 + (t207 * t268 + t208 * t232 + t209 * t233) * t229) / 0.2e1 + t227 * ((t196 * t266 + t198 * t234 + t200 * t235) * t228 + (t195 * t266 + t197 * t234 + t199 * t235) * t227 + (t219 * t266 + t220 * t234 + t221 * t235) * t256) / 0.2e1 + t204 * ((t186 * t266 + t188 * t230 + t190 * t231) * t205 + (t266 * t185 + t230 * t187 + t231 * t189) * t204 + (t207 * t266 + t208 * t230 + t209 * t231) * t229) / 0.2e1 + m(3) * (t251 ^ 2 + t257 ^ 2 + t258 ^ 2) / 0.2e1 + m(4) * (t182 ^ 2 + t193 ^ 2 + t194 ^ 2) / 0.2e1 + m(5) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + m(6) * (t176 ^ 2 + t177 ^ 2 + t178 ^ 2) / 0.2e1 + m(7) * (t173 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + ((t316 * t394 + t317 * t324 + t318 * t325 + t408) * V_base(6) + (t291 * t394 + t293 * t324 + t295 * t325 - t337 * t354 + t339 * t357 + Icges(1,4)) * V_base(5) + (t292 * t394 + t294 * t324 + t296 * t325 - t338 * t354 + t340 * t357 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t316 * t392 + t317 * t322 + t318 * t323 + t407) * V_base(6) + (-t291 * t392 + t293 * t322 + t295 * t323 + t337 * t357 + t339 * t354 + Icges(1,2)) * V_base(5) + (-t292 * t392 + t294 * t322 + t296 * t323 + t338 * t357 + t340 * t354 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t316 * t358 + (t317 * t356 + t318 * t353) * t355 + Icges(2,3) + Icges(1,3)) * V_base(6) + (t291 * t358 + (t293 * t356 + t295 * t353) * t355 + t407) * V_base(5) + (t292 * t358 + (t294 * t356 + t296 * t353) * t355 + t408) * V_base(4)) * V_base(6) / 0.2e1;
T  = t1;
