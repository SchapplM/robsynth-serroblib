% Calculate kinetic energy for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% m [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR15_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR15_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR15_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR15_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR15_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR15_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:23:53
% EndTime: 2024-09-27 22:23:55
% DurationCPUTime: 0.78s
% Computational Cost: add. (3517->420), mult. (3840->603), div. (0->0), fcn. (4050->26), ass. (0->219)
t363 = pkin(8) + pkin(9);
t357 = sin(qJ(1));
t413 = -t357 / 0.2e1;
t362 = cos(qJ(1));
t412 = t362 / 0.2e1;
t356 = sin(qJ(2));
t411 = pkin(2) * t356;
t355 = sin(qJ(3));
t410 = pkin(3) * t355;
t349 = sin(pkin(6));
t409 = pkin(11) * t349;
t352 = cos(pkin(5));
t408 = t352 * pkin(8);
t361 = cos(qJ(2));
t334 = t361 * pkin(2) + pkin(1);
t407 = -pkin(1) + t334;
t406 = Icges(2,4) * t357;
t348 = qJ(2) + qJ(3);
t344 = qJ(4) + t348;
t331 = sin(t344);
t405 = t331 * t357;
t404 = t331 * t362;
t332 = cos(t344);
t351 = cos(pkin(6));
t403 = t332 * t351;
t402 = t332 * t357;
t401 = t332 * t362;
t400 = t349 * t352;
t350 = sin(pkin(5));
t399 = t350 * t357;
t398 = t350 * t362;
t353 = sin(qJ(5));
t397 = t353 * t357;
t396 = t353 * t362;
t395 = t356 * t357;
t394 = t356 * t362;
t358 = cos(qJ(5));
t393 = t357 * t358;
t392 = t357 * t361;
t391 = t358 * t362;
t390 = t361 * t362;
t354 = sin(qJ(4));
t359 = cos(qJ(4));
t294 = pkin(4) * t359 + t354 * t409 + pkin(3);
t297 = -pkin(4) * t354 + t359 * t409;
t360 = cos(qJ(3));
t254 = t294 * t360 + t297 * t355 + pkin(2);
t255 = -t294 * t355 + t297 * t360;
t333 = pkin(3) * t360 + pkin(2);
t347 = pkin(10) + t363;
t382 = t361 * t410;
t261 = -t347 * t350 + (t333 * t356 + t382) * t352;
t317 = pkin(11) * t351 + t347;
t389 = t317 * t350 + (-t254 * t356 + t255 * t361) * t352 + t261;
t342 = cos(t348);
t300 = pkin(3) * t342 + t334;
t388 = t254 * t361 + t255 * t356 + pkin(1) - t300;
t293 = -t350 * t363 + t352 * t411;
t387 = t261 - t293;
t386 = t300 - t334;
t385 = qJD(2) * t350;
t384 = -qJD(2) - qJD(3);
t383 = V_base(5) * pkin(7) + V_base(1);
t379 = t349 * t399;
t378 = t349 * t398;
t377 = t352 * t397;
t376 = t352 * t396;
t375 = t352 * t393;
t374 = t352 * t391;
t305 = t357 * t385 + V_base(4);
t340 = V_base(6) + qJD(1);
t339 = pkin(5) - t348;
t338 = pkin(5) + t348;
t373 = pkin(8) * t350 + t293;
t278 = qJD(3) * t399 + t305;
t306 = qJD(2) * t352 + t340;
t264 = qJD(4) * t399 + t278;
t290 = qJD(3) * t352 + t306;
t298 = pkin(1) * t357 - pkin(8) * t398;
t372 = -t298 * t340 + V_base(5) * t408 + t383;
t275 = qJD(4) * t352 + t290;
t299 = pkin(1) * t362 + pkin(8) * t399;
t371 = V_base(4) * t298 - t299 * V_base(5) + V_base(3);
t263 = V_base(5) + (-qJD(4) + t384) * t398;
t370 = t340 * t299 + V_base(2) + (-pkin(7) - t408) * V_base(4);
t244 = t357 * t407 + t362 * t373;
t279 = t350 * t411 + (-pkin(8) + t363) * t352;
t304 = -t362 * t385 + V_base(5);
t369 = -t244 * t306 + t304 * t279 + t372;
t245 = -t357 * t373 + t362 * t407;
t368 = t305 * t244 - t245 * t304 + t371;
t367 = t306 * t245 - t279 * t305 + t370;
t202 = t357 * t386 + t362 * t387;
t239 = (t347 - t363) * t352 + (t382 + (-pkin(2) + t333) * t356) * t350;
t277 = t384 * t398 + V_base(5);
t366 = -t202 * t290 + t277 * t239 + t369;
t203 = -t357 * t387 + t362 * t386;
t365 = t278 * t202 - t203 * t277 + t368;
t364 = t290 * t203 - t239 * t278 + t367;
t343 = Icges(2,4) * t362;
t341 = sin(t348);
t330 = -qJ(4) + t339;
t329 = qJ(4) + t338;
t327 = cos(t339);
t326 = cos(t338);
t325 = sin(t339);
t324 = sin(t338);
t322 = cos(t329);
t321 = sin(t330);
t316 = cos(t330) / 0.2e1;
t315 = sin(t329) / 0.2e1;
t314 = rSges(2,1) * t362 - rSges(2,2) * t357;
t313 = rSges(2,1) * t357 + rSges(2,2) * t362;
t312 = Icges(2,1) * t362 - t406;
t311 = Icges(2,1) * t357 + t343;
t310 = -Icges(2,2) * t357 + t343;
t309 = Icges(2,2) * t362 + t406;
t303 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t302 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t301 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t296 = t327 + t326;
t295 = t324 - t325;
t292 = t327 / 0.2e1 - t326 / 0.2e1;
t291 = t324 / 0.2e1 + t325 / 0.2e1;
t288 = -t352 * t395 + t390;
t287 = -t352 * t392 - t394;
t286 = t352 * t394 + t392;
t285 = t352 * t390 - t395;
t284 = t316 - t322 / 0.2e1;
t283 = t316 + t322 / 0.2e1;
t282 = t315 - t321 / 0.2e1;
t281 = t315 + t321 / 0.2e1;
t276 = -t332 * t349 * t350 + t351 * t352;
t274 = rSges(3,3) * t352 + (rSges(3,1) * t356 + rSges(3,2) * t361) * t350;
t273 = Icges(3,5) * t352 + (Icges(3,1) * t356 + Icges(3,4) * t361) * t350;
t272 = Icges(3,6) * t352 + (Icges(3,4) * t356 + Icges(3,2) * t361) * t350;
t271 = Icges(3,3) * t352 + (Icges(3,5) * t356 + Icges(3,6) * t361) * t350;
t270 = t295 * t413 + t342 * t362;
t269 = t296 * t413 - t341 * t362;
t268 = t295 * t412 + t342 * t357;
t267 = t296 * t412 - t341 * t357;
t266 = V_base(5) * rSges(2,3) - t313 * t340 + t383;
t265 = t314 * t340 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t262 = t313 * V_base(4) - t314 * V_base(5) + V_base(3);
t260 = -t282 * t357 + t401;
t259 = -t283 * t357 - t404;
t258 = t282 * t362 + t402;
t257 = t283 * t362 - t405;
t253 = t353 * t400 + (t331 * t358 + t353 * t403) * t350;
t252 = t358 * t400 + (-t331 * t353 + t358 * t403) * t350;
t251 = -t351 * t398 + (-t352 * t401 + t405) * t349;
t250 = t351 * t399 + (t352 * t402 + t404) * t349;
t249 = rSges(4,1) * t292 + rSges(4,2) * t291 + rSges(4,3) * t352;
t248 = Icges(4,1) * t292 + Icges(4,4) * t291 + Icges(4,5) * t352;
t247 = Icges(4,4) * t292 + Icges(4,2) * t291 + Icges(4,6) * t352;
t246 = Icges(4,5) * t292 + Icges(4,6) * t291 + Icges(4,3) * t352;
t243 = qJD(5) * t276 + t275;
t242 = rSges(5,1) * t284 + rSges(5,2) * t281 + rSges(5,3) * t352;
t241 = rSges(3,1) * t288 + rSges(3,2) * t287 + rSges(3,3) * t399;
t240 = rSges(3,1) * t286 + rSges(3,2) * t285 - rSges(3,3) * t398;
t238 = Icges(5,1) * t284 + Icges(5,4) * t281 + Icges(5,5) * t352;
t237 = Icges(5,4) * t284 + Icges(5,2) * t281 + Icges(5,6) * t352;
t236 = Icges(5,5) * t284 + Icges(5,6) * t281 + Icges(5,3) * t352;
t235 = Icges(3,1) * t288 + Icges(3,4) * t287 + Icges(3,5) * t399;
t234 = Icges(3,1) * t286 + Icges(3,4) * t285 - Icges(3,5) * t398;
t233 = Icges(3,4) * t288 + Icges(3,2) * t287 + Icges(3,6) * t399;
t232 = Icges(3,4) * t286 + Icges(3,2) * t285 - Icges(3,6) * t398;
t231 = Icges(3,5) * t288 + Icges(3,6) * t287 + Icges(3,3) * t399;
t230 = Icges(3,5) * t286 + Icges(3,6) * t285 - Icges(3,3) * t398;
t226 = rSges(4,1) * t270 + rSges(4,2) * t269 + rSges(4,3) * t399;
t225 = rSges(4,1) * t268 + rSges(4,2) * t267 - rSges(4,3) * t398;
t224 = (t351 * t374 - t397) * t332 + (-t351 * t393 - t376) * t331 - t358 * t378;
t223 = (t351 * t376 + t393) * t332 + (-t351 * t397 + t374) * t331 - t353 * t378;
t222 = (-t351 * t375 - t396) * t332 + (-t351 * t391 + t377) * t331 + t358 * t379;
t221 = (-t351 * t377 + t391) * t332 + (-t351 * t396 - t375) * t331 + t353 * t379;
t220 = Icges(4,1) * t270 + Icges(4,4) * t269 + Icges(4,5) * t399;
t219 = Icges(4,1) * t268 + Icges(4,4) * t267 - Icges(4,5) * t398;
t218 = Icges(4,4) * t270 + Icges(4,2) * t269 + Icges(4,6) * t399;
t217 = Icges(4,4) * t268 + Icges(4,2) * t267 - Icges(4,6) * t398;
t216 = Icges(4,5) * t270 + Icges(4,6) * t269 + Icges(4,3) * t399;
t215 = Icges(4,5) * t268 + Icges(4,6) * t267 - Icges(4,3) * t398;
t214 = qJD(5) * t250 + t264;
t213 = qJD(5) * t251 + t263;
t212 = rSges(5,1) * t260 + rSges(5,2) * t259 + rSges(5,3) * t399;
t211 = rSges(5,1) * t258 + rSges(5,2) * t257 - rSges(5,3) * t398;
t210 = Icges(5,1) * t260 + Icges(5,4) * t259 + Icges(5,5) * t399;
t209 = Icges(5,1) * t258 + Icges(5,4) * t257 - Icges(5,5) * t398;
t208 = Icges(5,4) * t260 + Icges(5,2) * t259 + Icges(5,6) * t399;
t207 = Icges(5,4) * t258 + Icges(5,2) * t257 - Icges(5,6) * t398;
t206 = Icges(5,5) * t260 + Icges(5,6) * t259 + Icges(5,3) * t399;
t205 = Icges(5,5) * t258 + Icges(5,6) * t257 - Icges(5,3) * t398;
t198 = rSges(6,1) * t253 + rSges(6,2) * t252 + rSges(6,3) * t276;
t197 = Icges(6,1) * t253 + Icges(6,4) * t252 + Icges(6,5) * t276;
t196 = Icges(6,4) * t253 + Icges(6,2) * t252 + Icges(6,6) * t276;
t195 = Icges(6,5) * t253 + Icges(6,6) * t252 + Icges(6,3) * t276;
t194 = -t240 * t306 + t274 * t304 + t372;
t193 = t241 * t306 - t274 * t305 + t370;
t192 = (t317 - t347) * t352 + ((-t255 - t410) * t361 + (t254 - t333) * t356) * t350;
t191 = t240 * t305 - t241 * t304 + t371;
t190 = rSges(6,1) * t223 + rSges(6,2) * t224 + rSges(6,3) * t251;
t189 = rSges(6,1) * t221 + rSges(6,2) * t222 + rSges(6,3) * t250;
t188 = Icges(6,1) * t223 + Icges(6,4) * t224 + Icges(6,5) * t251;
t187 = Icges(6,1) * t221 + Icges(6,4) * t222 + Icges(6,5) * t250;
t186 = Icges(6,4) * t223 + Icges(6,2) * t224 + Icges(6,6) * t251;
t185 = Icges(6,4) * t221 + Icges(6,2) * t222 + Icges(6,6) * t250;
t184 = Icges(6,5) * t223 + Icges(6,6) * t224 + Icges(6,3) * t251;
t183 = Icges(6,5) * t221 + Icges(6,6) * t222 + Icges(6,3) * t250;
t182 = t357 * t389 + t362 * t388;
t181 = t357 * t388 - t362 * t389;
t180 = -t225 * t290 + t249 * t277 + t369;
t179 = t226 * t290 - t249 * t278 + t367;
t178 = t225 * t278 - t226 * t277 + t368;
t177 = -t211 * t275 + t242 * t263 + t366;
t176 = t212 * t275 - t242 * t264 + t364;
t175 = t211 * t264 - t212 * t263 + t365;
t174 = -t181 * t275 - t190 * t243 + t192 * t263 + t198 * t213 + t366;
t173 = t182 * t275 + t189 * t243 - t192 * t264 - t198 * t214 + t364;
t172 = t181 * t264 - t182 * t263 - t189 * t213 + t190 * t214 + t365;
t1 = m(1) * (t301 ^ 2 + t302 ^ 2 + t303 ^ 2) / 0.2e1 + m(2) * (t262 ^ 2 + t265 ^ 2 + t266 ^ 2) / 0.2e1 + m(3) * (t191 ^ 2 + t193 ^ 2 + t194 ^ 2) / 0.2e1 + t305 * ((t231 * t399 + t233 * t287 + t235 * t288) * t305 + (t230 * t399 + t232 * t287 + t234 * t288) * t304 + (t271 * t399 + t272 * t287 + t273 * t288) * t306) / 0.2e1 + t304 * ((-t231 * t398 + t233 * t285 + t235 * t286) * t305 + (-t230 * t398 + t232 * t285 + t234 * t286) * t304 + (-t271 * t398 + t272 * t285 + t273 * t286) * t306) / 0.2e1 + t306 * ((t230 * t304 + t231 * t305 + t271 * t306) * t352 + ((t233 * t361 + t235 * t356) * t305 + (t232 * t361 + t234 * t356) * t304 + (t272 * t361 + t273 * t356) * t306) * t350) / 0.2e1 + m(4) * (t178 ^ 2 + t179 ^ 2 + t180 ^ 2) / 0.2e1 + t278 * ((t216 * t399 + t218 * t269 + t220 * t270) * t278 + (t215 * t399 + t217 * t269 + t219 * t270) * t277 + (t246 * t399 + t247 * t269 + t248 * t270) * t290) / 0.2e1 + t277 * ((-t216 * t398 + t218 * t267 + t220 * t268) * t278 + (-t215 * t398 + t217 * t267 + t219 * t268) * t277 + (-t246 * t398 + t247 * t267 + t248 * t268) * t290) / 0.2e1 + t290 * ((t216 * t352 + t218 * t291 + t220 * t292) * t278 + (t215 * t352 + t217 * t291 + t219 * t292) * t277 + (t246 * t352 + t247 * t291 + t248 * t292) * t290) / 0.2e1 + m(5) * (t175 ^ 2 + t176 ^ 2 + t177 ^ 2) / 0.2e1 + t264 * ((t206 * t399 + t208 * t259 + t210 * t260) * t264 + (t205 * t399 + t207 * t259 + t209 * t260) * t263 + (t236 * t399 + t237 * t259 + t238 * t260) * t275) / 0.2e1 + t263 * ((-t206 * t398 + t208 * t257 + t210 * t258) * t264 + (-t205 * t398 + t207 * t257 + t209 * t258) * t263 + (-t236 * t398 + t237 * t257 + t238 * t258) * t275) / 0.2e1 + t275 * ((t206 * t352 + t208 * t281 + t210 * t284) * t264 + (t205 * t352 + t207 * t281 + t209 * t284) * t263 + (t236 * t352 + t237 * t281 + t238 * t284) * t275) / 0.2e1 + m(6) * (t172 ^ 2 + t173 ^ 2 + t174 ^ 2) / 0.2e1 + t214 * ((t250 * t183 + t222 * t185 + t221 * t187) * t214 + (t184 * t250 + t186 * t222 + t188 * t221) * t213 + (t195 * t250 + t196 * t222 + t197 * t221) * t243) / 0.2e1 + t213 * ((t183 * t251 + t185 * t224 + t187 * t223) * t214 + (t251 * t184 + t224 * t186 + t223 * t188) * t213 + (t195 * t251 + t196 * t224 + t197 * t223) * t243) / 0.2e1 + t243 * ((t183 * t276 + t185 * t252 + t187 * t253) * t214 + (t184 * t276 + t186 * t252 + t188 * t253) * t213 + (t195 * t276 + t196 * t252 + t197 * t253) * t243) / 0.2e1 + ((-t309 * t357 + t311 * t362 + Icges(1,4)) * V_base(5) + (-t310 * t357 + t312 * t362 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t309 * t362 + t311 * t357 + Icges(1,2)) * V_base(5) + (t310 * t362 + t312 * t357 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t357 + Icges(2,6) * t362) * V_base(5) + (Icges(2,5) * t362 - Icges(2,6) * t357) * V_base(4) + Icges(2,3) * t340 / 0.2e1) * t340;
T = t1;
