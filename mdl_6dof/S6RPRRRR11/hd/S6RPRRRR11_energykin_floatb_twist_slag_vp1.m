% Calculate kinetic energy for
% S6RPRRRR11
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR11_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR11_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR11_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRR11_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRRR11_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR11_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR11_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR11_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:38:11
% EndTime: 2019-03-09 07:38:16
% DurationCPUTime: 4.86s
% Computational Cost: add. (5820->452), mult. (14861->664), div. (0->0), fcn. (19107->16), ass. (0->197)
t404 = cos(qJ(3));
t403 = cos(qJ(4));
t361 = cos(qJ(5));
t402 = pkin(5) * t361;
t400 = cos(pkin(7));
t399 = sin(pkin(7));
t360 = sin(qJ(1));
t398 = Icges(2,4) * t360;
t356 = cos(pkin(6));
t397 = qJ(2) * t356;
t353 = sin(pkin(13));
t355 = cos(pkin(13));
t362 = cos(qJ(1));
t389 = t356 * t362;
t322 = -t353 * t360 + t355 * t389;
t323 = t353 * t389 + t355 * t360;
t359 = sin(qJ(3));
t354 = sin(pkin(6));
t379 = t404 * t399;
t377 = t354 * t379;
t380 = t400 * t404;
t283 = -t322 * t380 + t323 * t359 + t362 * t377;
t357 = sin(qJ(5));
t396 = t283 * t357;
t390 = t356 * t360;
t324 = -t353 * t362 - t355 * t390;
t325 = -t353 * t390 + t355 * t362;
t285 = -t324 * t380 + t325 * t359 - t360 * t377;
t395 = t285 * t357;
t393 = t353 * t354;
t305 = -t354 * t355 * t380 - t356 * t379 + t359 * t393;
t394 = t305 * t357;
t392 = t354 * t360;
t391 = t354 * t362;
t388 = qJD(2) * t354;
t387 = V_base(5) * pkin(8) + V_base(1);
t383 = t354 * t400;
t374 = -t324 * t399 + t360 * t383;
t299 = qJD(3) * t374 + V_base(4);
t375 = -t322 * t399 - t362 * t383;
t298 = qJD(3) * t375 + V_base(5);
t347 = V_base(6) + qJD(1);
t384 = -pkin(8) - t397;
t382 = t354 * t399;
t327 = pkin(1) * t360 - qJ(2) * t391;
t381 = qJD(2) * t356 + V_base(4) * t327 + V_base(3);
t266 = qJD(4) * t285 + t299;
t265 = qJD(4) * t283 + t298;
t373 = -t355 * t382 + t356 * t400;
t312 = qJD(3) * t373 + t347;
t378 = t360 * t388 + V_base(5) * t397 + t387;
t286 = t325 * t404 + (t324 * t400 + t360 * t382) * t359;
t358 = sin(qJ(4));
t269 = t286 * t358 - t374 * t403;
t227 = qJD(5) * t269 + t266;
t284 = t323 * t404 + (t322 * t400 - t362 * t382) * t359;
t267 = t284 * t358 - t375 * t403;
t226 = qJD(5) * t267 + t265;
t277 = qJD(4) * t305 + t312;
t306 = t356 * t399 * t359 + (t355 * t359 * t400 + t353 * t404) * t354;
t280 = t306 * t358 - t373 * t403;
t253 = qJD(5) * t280 + t277;
t328 = pkin(1) * t362 + qJ(2) * t392;
t376 = t347 * t328 - t362 * t388 + V_base(2);
t288 = t323 * pkin(2) + pkin(9) * t375;
t309 = pkin(2) * t393 + pkin(9) * t373;
t372 = V_base(5) * t309 + (-t288 - t327) * t347 + t378;
t289 = t325 * pkin(2) + pkin(9) * t374;
t371 = V_base(4) * t288 + (-t289 - t328) * V_base(5) + t381;
t256 = pkin(3) * t284 + pkin(10) * t283;
t275 = pkin(3) * t306 + pkin(10) * t305;
t370 = -t256 * t312 + t298 * t275 + t372;
t257 = pkin(3) * t286 + pkin(10) * t285;
t369 = t299 * t256 - t257 * t298 + t371;
t368 = t347 * t289 + (-t309 + t384) * V_base(4) + t376;
t268 = t284 * t403 + t358 * t375;
t224 = t268 * pkin(4) + t267 * pkin(11);
t281 = t306 * t403 + t358 * t373;
t252 = t281 * pkin(4) + t280 * pkin(11);
t367 = -t224 * t277 + t265 * t252 + t370;
t270 = t286 * t403 + t358 * t374;
t225 = t270 * pkin(4) + t269 * pkin(11);
t366 = t266 * t224 - t225 * t265 + t369;
t365 = t312 * t257 - t275 * t299 + t368;
t364 = t277 * t225 - t252 * t266 + t365;
t352 = qJ(5) + qJ(6);
t350 = Icges(2,4) * t362;
t349 = cos(t352);
t348 = sin(t352);
t341 = rSges(2,1) * t362 - rSges(2,2) * t360;
t340 = rSges(2,1) * t360 + rSges(2,2) * t362;
t339 = Icges(2,1) * t362 - t398;
t338 = Icges(2,1) * t360 + t350;
t337 = -Icges(2,2) * t360 + t350;
t336 = Icges(2,2) * t362 + t398;
t335 = Icges(2,5) * t362 - Icges(2,6) * t360;
t334 = Icges(2,5) * t360 + Icges(2,6) * t362;
t333 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t332 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t331 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t318 = rSges(3,3) * t356 + (rSges(3,1) * t353 + rSges(3,2) * t355) * t354;
t317 = Icges(3,5) * t356 + (Icges(3,1) * t353 + Icges(3,4) * t355) * t354;
t316 = Icges(3,6) * t356 + (Icges(3,4) * t353 + Icges(3,2) * t355) * t354;
t315 = Icges(3,3) * t356 + (Icges(3,5) * t353 + Icges(3,6) * t355) * t354;
t311 = V_base(5) * rSges(2,3) - t340 * t347 + t387;
t310 = t341 * t347 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t308 = t340 * V_base(4) - t341 * V_base(5) + V_base(3);
t297 = rSges(3,1) * t325 + rSges(3,2) * t324 + rSges(3,3) * t392;
t296 = rSges(3,1) * t323 + rSges(3,2) * t322 - rSges(3,3) * t391;
t295 = Icges(3,1) * t325 + Icges(3,4) * t324 + Icges(3,5) * t392;
t294 = Icges(3,1) * t323 + Icges(3,4) * t322 - Icges(3,5) * t391;
t293 = Icges(3,4) * t325 + Icges(3,2) * t324 + Icges(3,6) * t392;
t292 = Icges(3,4) * t323 + Icges(3,2) * t322 - Icges(3,6) * t391;
t291 = Icges(3,5) * t325 + Icges(3,6) * t324 + Icges(3,3) * t392;
t290 = Icges(3,5) * t323 + Icges(3,6) * t322 - Icges(3,3) * t391;
t274 = t306 * rSges(4,1) - t305 * rSges(4,2) + rSges(4,3) * t373;
t273 = Icges(4,1) * t306 - Icges(4,4) * t305 + Icges(4,5) * t373;
t272 = Icges(4,4) * t306 - Icges(4,2) * t305 + Icges(4,6) * t373;
t271 = Icges(4,5) * t306 - Icges(4,6) * t305 + Icges(4,3) * t373;
t262 = t281 * t361 + t394;
t261 = -t281 * t357 + t305 * t361;
t260 = t281 * t349 + t305 * t348;
t259 = -t281 * t348 + t305 * t349;
t255 = t318 * V_base(5) + (-t296 - t327) * t347 + t378;
t254 = t297 * t347 + (-t318 + t384) * V_base(4) + t376;
t251 = t296 * V_base(4) + (-t297 - t328) * V_base(5) + t381;
t249 = t286 * rSges(4,1) - t285 * rSges(4,2) + rSges(4,3) * t374;
t248 = t284 * rSges(4,1) - t283 * rSges(4,2) + rSges(4,3) * t375;
t247 = Icges(4,1) * t286 - Icges(4,4) * t285 + Icges(4,5) * t374;
t246 = Icges(4,1) * t284 - Icges(4,4) * t283 + Icges(4,5) * t375;
t245 = Icges(4,4) * t286 - Icges(4,2) * t285 + Icges(4,6) * t374;
t244 = Icges(4,4) * t284 - Icges(4,2) * t283 + Icges(4,6) * t375;
t243 = Icges(4,5) * t286 - Icges(4,6) * t285 + Icges(4,3) * t374;
t242 = Icges(4,5) * t284 - Icges(4,6) * t283 + Icges(4,3) * t375;
t241 = rSges(5,1) * t281 - rSges(5,2) * t280 + rSges(5,3) * t305;
t239 = Icges(5,1) * t281 - Icges(5,4) * t280 + Icges(5,5) * t305;
t238 = Icges(5,4) * t281 - Icges(5,2) * t280 + Icges(5,6) * t305;
t237 = Icges(5,5) * t281 - Icges(5,6) * t280 + Icges(5,3) * t305;
t236 = t270 * t361 + t395;
t235 = -t270 * t357 + t285 * t361;
t234 = t268 * t361 + t396;
t233 = -t268 * t357 + t283 * t361;
t232 = t270 * t349 + t285 * t348;
t231 = -t270 * t348 + t285 * t349;
t230 = t268 * t349 + t283 * t348;
t229 = -t268 * t348 + t283 * t349;
t228 = qJD(6) * t280 + t253;
t222 = rSges(5,1) * t270 - rSges(5,2) * t269 + rSges(5,3) * t285;
t221 = rSges(5,1) * t268 - rSges(5,2) * t267 + rSges(5,3) * t283;
t220 = Icges(5,1) * t270 - Icges(5,4) * t269 + Icges(5,5) * t285;
t219 = Icges(5,1) * t268 - Icges(5,4) * t267 + Icges(5,5) * t283;
t218 = Icges(5,4) * t270 - Icges(5,2) * t269 + Icges(5,6) * t285;
t217 = Icges(5,4) * t268 - Icges(5,2) * t267 + Icges(5,6) * t283;
t216 = Icges(5,5) * t270 - Icges(5,6) * t269 + Icges(5,3) * t285;
t215 = Icges(5,5) * t268 - Icges(5,6) * t267 + Icges(5,3) * t283;
t213 = rSges(6,1) * t262 + rSges(6,2) * t261 + rSges(6,3) * t280;
t212 = Icges(6,1) * t262 + Icges(6,4) * t261 + Icges(6,5) * t280;
t211 = Icges(6,4) * t262 + Icges(6,2) * t261 + Icges(6,6) * t280;
t210 = Icges(6,5) * t262 + Icges(6,6) * t261 + Icges(6,3) * t280;
t209 = rSges(7,1) * t260 + rSges(7,2) * t259 + rSges(7,3) * t280;
t208 = Icges(7,1) * t260 + Icges(7,4) * t259 + Icges(7,5) * t280;
t207 = Icges(7,4) * t260 + Icges(7,2) * t259 + Icges(7,6) * t280;
t206 = Icges(7,5) * t260 + Icges(7,6) * t259 + Icges(7,3) * t280;
t205 = pkin(5) * t394 + pkin(12) * t280 + t281 * t402;
t204 = qJD(6) * t269 + t227;
t203 = qJD(6) * t267 + t226;
t201 = rSges(6,1) * t236 + rSges(6,2) * t235 + rSges(6,3) * t269;
t200 = rSges(6,1) * t234 + rSges(6,2) * t233 + rSges(6,3) * t267;
t199 = Icges(6,1) * t236 + Icges(6,4) * t235 + Icges(6,5) * t269;
t198 = Icges(6,1) * t234 + Icges(6,4) * t233 + Icges(6,5) * t267;
t197 = Icges(6,4) * t236 + Icges(6,2) * t235 + Icges(6,6) * t269;
t196 = Icges(6,4) * t234 + Icges(6,2) * t233 + Icges(6,6) * t267;
t195 = Icges(6,5) * t236 + Icges(6,6) * t235 + Icges(6,3) * t269;
t194 = Icges(6,5) * t234 + Icges(6,6) * t233 + Icges(6,3) * t267;
t193 = rSges(7,1) * t232 + rSges(7,2) * t231 + rSges(7,3) * t269;
t192 = rSges(7,1) * t230 + rSges(7,2) * t229 + rSges(7,3) * t267;
t191 = Icges(7,1) * t232 + Icges(7,4) * t231 + Icges(7,5) * t269;
t190 = Icges(7,1) * t230 + Icges(7,4) * t229 + Icges(7,5) * t267;
t189 = Icges(7,4) * t232 + Icges(7,2) * t231 + Icges(7,6) * t269;
t188 = Icges(7,4) * t230 + Icges(7,2) * t229 + Icges(7,6) * t267;
t187 = Icges(7,5) * t232 + Icges(7,6) * t231 + Icges(7,3) * t269;
t186 = Icges(7,5) * t230 + Icges(7,6) * t229 + Icges(7,3) * t267;
t185 = -t248 * t312 + t274 * t298 + t372;
t184 = t249 * t312 - t274 * t299 + t368;
t183 = pkin(5) * t395 + pkin(12) * t269 + t270 * t402;
t182 = pkin(5) * t396 + pkin(12) * t267 + t268 * t402;
t181 = t248 * t299 - t249 * t298 + t371;
t180 = -t221 * t277 + t241 * t265 + t370;
t179 = t222 * t277 - t241 * t266 + t365;
t178 = t221 * t266 - t222 * t265 + t369;
t177 = -t200 * t253 + t213 * t226 + t367;
t176 = t201 * t253 - t213 * t227 + t364;
t175 = t200 * t227 - t201 * t226 + t366;
t174 = -t182 * t253 - t192 * t228 + t203 * t209 + t205 * t226 + t367;
t173 = t183 * t253 + t193 * t228 - t204 * t209 - t205 * t227 + t364;
t172 = t182 * t227 - t183 * t226 + t192 * t204 - t193 * t203 + t366;
t1 = m(1) * (t331 ^ 2 + t332 ^ 2 + t333 ^ 2) / 0.2e1 + t277 * ((t216 * t305 - t218 * t280 + t220 * t281) * t266 + (t215 * t305 - t217 * t280 + t219 * t281) * t265 + (t237 * t305 - t238 * t280 + t239 * t281) * t277) / 0.2e1 + m(2) * (t308 ^ 2 + t310 ^ 2 + t311 ^ 2) / 0.2e1 + t266 * ((t216 * t285 - t218 * t269 + t220 * t270) * t266 + (t215 * t285 - t217 * t269 + t219 * t270) * t265 + (t237 * t285 - t238 * t269 + t239 * t270) * t277) / 0.2e1 + t265 * ((t216 * t283 - t218 * t267 + t220 * t268) * t266 + (t215 * t283 - t217 * t267 + t219 * t268) * t265 + (t237 * t283 - t238 * t267 + t239 * t268) * t277) / 0.2e1 + t253 * ((t195 * t280 + t197 * t261 + t199 * t262) * t227 + (t194 * t280 + t196 * t261 + t198 * t262) * t226 + (t210 * t280 + t211 * t261 + t212 * t262) * t253) / 0.2e1 + t228 * ((t187 * t280 + t189 * t259 + t191 * t260) * t204 + (t186 * t280 + t188 * t259 + t190 * t260) * t203 + (t280 * t206 + t259 * t207 + t260 * t208) * t228) / 0.2e1 + t226 * ((t195 * t267 + t197 * t233 + t199 * t234) * t227 + (t267 * t194 + t233 * t196 + t234 * t198) * t226 + (t210 * t267 + t211 * t233 + t212 * t234) * t253) / 0.2e1 + t203 * ((t187 * t267 + t189 * t229 + t191 * t230) * t204 + (t267 * t186 + t229 * t188 + t230 * t190) * t203 + (t206 * t267 + t207 * t229 + t208 * t230) * t228) / 0.2e1 + t227 * ((t269 * t195 + t235 * t197 + t236 * t199) * t227 + (t194 * t269 + t196 * t235 + t198 * t236) * t226 + (t210 * t269 + t211 * t235 + t212 * t236) * t253) / 0.2e1 + t204 * ((t269 * t187 + t231 * t189 + t232 * t191) * t204 + (t186 * t269 + t188 * t231 + t190 * t232) * t203 + (t206 * t269 + t207 * t231 + t208 * t232) * t228) / 0.2e1 + m(3) * (t251 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + t312 * ((t243 * t373 - t305 * t245 + t306 * t247) * t299 + (t242 * t373 - t305 * t244 + t306 * t246) * t298 + (t271 * t373 - t305 * t272 + t306 * t273) * t312) / 0.2e1 + t299 * ((t243 * t374 - t285 * t245 + t286 * t247) * t299 + (t242 * t374 - t285 * t244 + t286 * t246) * t298 + (t271 * t374 - t285 * t272 + t286 * t273) * t312) / 0.2e1 + t298 * ((t243 * t375 - t283 * t245 + t284 * t247) * t299 + (t242 * t375 - t283 * t244 + t284 * t246) * t298 + (t271 * t375 - t283 * t272 + t284 * t273) * t312) / 0.2e1 + m(7) * (t172 ^ 2 + t173 ^ 2 + t174 ^ 2) / 0.2e1 + m(6) * (t175 ^ 2 + t176 ^ 2 + t177 ^ 2) / 0.2e1 + m(5) * (t178 ^ 2 + t179 ^ 2 + t180 ^ 2) / 0.2e1 + m(4) * (t181 ^ 2 + t184 ^ 2 + t185 ^ 2) / 0.2e1 + ((t290 * V_base(5) + t291 * V_base(4) + t315 * t347) * t356 + ((t293 * t355 + t295 * t353) * V_base(4) + (t292 * t355 + t294 * t353) * V_base(5) + (t316 * t355 + t317 * t353) * t347) * t354 + Icges(2,3) * t347 + t334 * V_base(5) + t335 * V_base(4)) * t347 / 0.2e1 + ((t315 * t392 + t316 * t324 + t317 * t325 + t335) * t347 + (t290 * t392 + t292 * t324 + t294 * t325 - t336 * t360 + t338 * t362 + Icges(1,4)) * V_base(5) + (t291 * t392 + t293 * t324 + t295 * t325 - t337 * t360 + t339 * t362 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t315 * t391 + t316 * t322 + t317 * t323 + t334) * t347 + (-t290 * t391 + t292 * t322 + t294 * t323 + t336 * t362 + t338 * t360 + Icges(1,2)) * V_base(5) + (-t291 * t391 + t293 * t322 + t295 * t323 + t337 * t362 + t339 * t360 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
