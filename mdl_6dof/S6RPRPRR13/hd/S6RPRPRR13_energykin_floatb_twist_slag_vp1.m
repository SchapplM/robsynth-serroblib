% Calculate kinetic energy for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR13_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR13_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRR13_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_energykin_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR13_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR13_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR13_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:21:36
% EndTime: 2019-03-09 04:21:41
% DurationCPUTime: 4.78s
% Computational Cost: add. (4291->406), mult. (11179->574), div. (0->0), fcn. (14157->14), ass. (0->181)
t415 = Icges(4,1) + Icges(5,2);
t414 = Icges(5,1) + Icges(4,3);
t413 = -Icges(4,4) - Icges(5,6);
t412 = -Icges(5,4) + Icges(4,5);
t411 = Icges(5,5) - Icges(4,6);
t410 = Icges(4,2) + Icges(5,3);
t353 = cos(pkin(12));
t351 = sin(pkin(12));
t358 = sin(qJ(1));
t386 = t358 * t351;
t354 = cos(pkin(6));
t360 = cos(qJ(1));
t387 = t354 * t360;
t323 = t353 * t387 - t386;
t385 = t358 * t353;
t324 = t351 * t387 + t385;
t357 = sin(qJ(3));
t352 = sin(pkin(6));
t393 = sin(pkin(7));
t396 = cos(qJ(3));
t373 = t396 * t393;
t371 = t352 * t373;
t394 = cos(pkin(7));
t374 = t394 * t396;
t281 = -t323 * t374 + t324 * t357 + t360 * t371;
t376 = t357 * t393;
t377 = t357 * t394;
t388 = t352 * t360;
t282 = t323 * t377 + t324 * t396 - t376 * t388;
t379 = t352 * t394;
t306 = -t323 * t393 - t360 * t379;
t409 = t410 * t281 + t413 * t282 + t411 * t306;
t325 = -t351 * t360 - t354 * t385;
t326 = t353 * t360 - t354 * t386;
t283 = -t325 * t374 + t326 * t357 - t358 * t371;
t378 = t352 * t393;
t284 = t326 * t396 + (t325 * t394 + t358 * t378) * t357;
t307 = -t325 * t393 + t358 * t379;
t408 = t410 * t283 + t413 * t284 + t411 * t307;
t407 = t411 * t281 + t412 * t282 + t414 * t306;
t406 = t411 * t283 + t412 * t284 + t414 * t307;
t405 = t413 * t281 + t415 * t282 + t412 * t306;
t404 = t413 * t283 + t415 * t284 + t412 * t307;
t390 = t351 * t352;
t304 = -t352 * t353 * t374 - t354 * t373 + t357 * t390;
t305 = t354 * t376 + (t351 * t396 + t353 * t377) * t352;
t322 = -t353 * t378 + t354 * t394;
t403 = t410 * t304 + t413 * t305 + t411 * t322;
t402 = t411 * t304 + t412 * t305 + t414 * t322;
t401 = t413 * t304 + t415 * t305 + t412 * t322;
t395 = cos(qJ(5));
t392 = Icges(2,4) * t358;
t391 = qJ(2) * t354;
t389 = t352 * t358;
t384 = qJD(2) * t352;
t383 = V_base(5) * pkin(8) + V_base(1);
t298 = qJD(3) * t307 + V_base(4);
t297 = qJD(3) * t306 + V_base(5);
t348 = V_base(6) + qJD(1);
t380 = -pkin(8) - t391;
t328 = t358 * pkin(1) - qJ(2) * t388;
t375 = qJD(2) * t354 + V_base(4) * t328 + V_base(3);
t254 = qJD(5) * t284 + t298;
t253 = qJD(5) * t282 + t297;
t313 = qJD(3) * t322 + t348;
t372 = t358 * t384 + V_base(5) * t391 + t383;
t273 = qJD(5) * t305 + t313;
t329 = pkin(1) * t360 + qJ(2) * t389;
t370 = t348 * t329 - t360 * t384 + V_base(2);
t287 = t324 * pkin(2) + pkin(9) * t306;
t310 = pkin(2) * t390 + pkin(9) * t322;
t369 = V_base(5) * t310 + (-t287 - t328) * t348 + t372;
t288 = t326 * pkin(2) + pkin(9) * t307;
t368 = V_base(4) * t287 + (-t288 - t329) * V_base(5) + t375;
t270 = pkin(3) * t305 + qJ(4) * t304;
t367 = qJD(4) * t283 + t297 * t270 + t369;
t245 = pkin(3) * t282 + qJ(4) * t281;
t366 = qJD(4) * t304 + t298 * t245 + t368;
t365 = t348 * t288 + (-t310 + t380) * V_base(4) + t370;
t246 = pkin(3) * t284 + qJ(4) * t283;
t364 = qJD(4) * t281 + t313 * t246 + t365;
t259 = pkin(4) * t306 + pkin(10) * t282;
t286 = pkin(4) * t322 + pkin(10) * t305;
t363 = t297 * t286 + (-t245 - t259) * t313 + t367;
t260 = pkin(4) * t307 + pkin(10) * t284;
t362 = t298 * t259 + (-t246 - t260) * t297 + t366;
t361 = t313 * t260 + (-t270 - t286) * t298 + t364;
t359 = cos(qJ(6));
t356 = sin(qJ(5));
t355 = sin(qJ(6));
t349 = Icges(2,4) * t360;
t343 = rSges(2,1) * t360 - t358 * rSges(2,2);
t342 = t358 * rSges(2,1) + rSges(2,2) * t360;
t341 = Icges(2,1) * t360 - t392;
t340 = Icges(2,1) * t358 + t349;
t339 = -Icges(2,2) * t358 + t349;
t338 = Icges(2,2) * t360 + t392;
t337 = Icges(2,5) * t360 - Icges(2,6) * t358;
t336 = Icges(2,5) * t358 + Icges(2,6) * t360;
t335 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t334 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t333 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t319 = rSges(3,3) * t354 + (rSges(3,1) * t351 + rSges(3,2) * t353) * t352;
t318 = Icges(3,5) * t354 + (Icges(3,1) * t351 + Icges(3,4) * t353) * t352;
t317 = Icges(3,6) * t354 + (Icges(3,4) * t351 + Icges(3,2) * t353) * t352;
t316 = Icges(3,3) * t354 + (Icges(3,5) * t351 + Icges(3,6) * t353) * t352;
t312 = V_base(5) * rSges(2,3) - t342 * t348 + t383;
t311 = t343 * t348 + V_base(2) + (-rSges(2,3) - pkin(8)) * V_base(4);
t309 = t342 * V_base(4) - t343 * V_base(5) + V_base(3);
t296 = rSges(3,1) * t326 + rSges(3,2) * t325 + rSges(3,3) * t389;
t295 = t324 * rSges(3,1) + t323 * rSges(3,2) - rSges(3,3) * t388;
t294 = Icges(3,1) * t326 + Icges(3,4) * t325 + Icges(3,5) * t389;
t293 = Icges(3,1) * t324 + Icges(3,4) * t323 - Icges(3,5) * t388;
t292 = Icges(3,4) * t326 + Icges(3,2) * t325 + Icges(3,6) * t389;
t291 = Icges(3,4) * t324 + Icges(3,2) * t323 - Icges(3,6) * t388;
t290 = Icges(3,5) * t326 + Icges(3,6) * t325 + Icges(3,3) * t389;
t289 = Icges(3,5) * t324 + Icges(3,6) * t323 - Icges(3,3) * t388;
t279 = t304 * t356 + t322 * t395;
t278 = -t304 * t395 + t322 * t356;
t269 = rSges(5,1) * t322 - rSges(5,2) * t305 + rSges(5,3) * t304;
t268 = rSges(4,1) * t305 - rSges(4,2) * t304 + rSges(4,3) * t322;
t258 = t283 * t356 + t307 * t395;
t257 = -t283 * t395 + t307 * t356;
t256 = t281 * t356 + t306 * t395;
t255 = -t281 * t395 + t306 * t356;
t252 = t279 * t359 + t305 * t355;
t251 = -t279 * t355 + t305 * t359;
t248 = t319 * V_base(5) + (-t295 - t328) * t348 + t372;
t247 = t348 * t296 + (-t319 + t380) * V_base(4) + t370;
t244 = qJD(6) * t278 + t273;
t243 = pkin(5) * t279 + pkin(11) * t278;
t241 = t295 * V_base(4) + (-t296 - t329) * V_base(5) + t375;
t239 = rSges(5,1) * t307 - rSges(5,2) * t284 + rSges(5,3) * t283;
t238 = rSges(5,1) * t306 - rSges(5,2) * t282 + rSges(5,3) * t281;
t237 = rSges(4,1) * t284 - rSges(4,2) * t283 + rSges(4,3) * t307;
t236 = rSges(4,1) * t282 - rSges(4,2) * t281 + rSges(4,3) * t306;
t223 = rSges(6,1) * t279 - rSges(6,2) * t278 + rSges(6,3) * t305;
t222 = Icges(6,1) * t279 - Icges(6,4) * t278 + Icges(6,5) * t305;
t221 = Icges(6,4) * t279 - Icges(6,2) * t278 + Icges(6,6) * t305;
t220 = Icges(6,5) * t279 - Icges(6,6) * t278 + Icges(6,3) * t305;
t218 = t258 * t359 + t284 * t355;
t217 = -t258 * t355 + t284 * t359;
t216 = t256 * t359 + t282 * t355;
t215 = -t256 * t355 + t282 * t359;
t214 = qJD(6) * t257 + t254;
t213 = qJD(6) * t255 + t253;
t212 = pkin(5) * t258 + pkin(11) * t257;
t211 = pkin(5) * t256 + pkin(11) * t255;
t210 = rSges(6,1) * t258 - rSges(6,2) * t257 + rSges(6,3) * t284;
t209 = rSges(6,1) * t256 - rSges(6,2) * t255 + rSges(6,3) * t282;
t208 = Icges(6,1) * t258 - Icges(6,4) * t257 + Icges(6,5) * t284;
t207 = Icges(6,1) * t256 - Icges(6,4) * t255 + Icges(6,5) * t282;
t206 = Icges(6,4) * t258 - Icges(6,2) * t257 + Icges(6,6) * t284;
t205 = Icges(6,4) * t256 - Icges(6,2) * t255 + Icges(6,6) * t282;
t204 = Icges(6,5) * t258 - Icges(6,6) * t257 + Icges(6,3) * t284;
t203 = Icges(6,5) * t256 - Icges(6,6) * t255 + Icges(6,3) * t282;
t202 = rSges(7,1) * t252 + rSges(7,2) * t251 + rSges(7,3) * t278;
t201 = Icges(7,1) * t252 + Icges(7,4) * t251 + Icges(7,5) * t278;
t200 = Icges(7,4) * t252 + Icges(7,2) * t251 + Icges(7,6) * t278;
t199 = Icges(7,5) * t252 + Icges(7,6) * t251 + Icges(7,3) * t278;
t198 = rSges(7,1) * t218 + rSges(7,2) * t217 + rSges(7,3) * t257;
t197 = rSges(7,1) * t216 + rSges(7,2) * t215 + rSges(7,3) * t255;
t196 = Icges(7,1) * t218 + Icges(7,4) * t217 + Icges(7,5) * t257;
t195 = Icges(7,1) * t216 + Icges(7,4) * t215 + Icges(7,5) * t255;
t194 = Icges(7,4) * t218 + Icges(7,2) * t217 + Icges(7,6) * t257;
t193 = Icges(7,4) * t216 + Icges(7,2) * t215 + Icges(7,6) * t255;
t192 = Icges(7,5) * t218 + Icges(7,6) * t217 + Icges(7,3) * t257;
t191 = Icges(7,5) * t216 + Icges(7,6) * t215 + Icges(7,3) * t255;
t190 = -t236 * t313 + t268 * t297 + t369;
t189 = t313 * t237 - t298 * t268 + t365;
t188 = t236 * t298 - t237 * t297 + t368;
t187 = t269 * t297 + (-t238 - t245) * t313 + t367;
t186 = t313 * t239 + (-t269 - t270) * t298 + t364;
t185 = t238 * t298 + (-t239 - t246) * t297 + t366;
t184 = -t209 * t273 + t223 * t253 + t363;
t183 = t273 * t210 - t254 * t223 + t361;
t182 = t209 * t254 - t210 * t253 + t362;
t181 = -t197 * t244 + t202 * t213 - t211 * t273 + t243 * t253 + t363;
t180 = t244 * t198 - t214 * t202 + t273 * t212 - t254 * t243 + t361;
t179 = t197 * t214 - t198 * t213 + t211 * t254 - t212 * t253 + t362;
t1 = m(1) * (t333 ^ 2 + t334 ^ 2 + t335 ^ 2) / 0.2e1 + m(2) * (t309 ^ 2 + t311 ^ 2 + t312 ^ 2) / 0.2e1 + t273 * ((t204 * t305 - t206 * t278 + t208 * t279) * t254 + (t203 * t305 - t205 * t278 + t207 * t279) * t253 + (t305 * t220 - t278 * t221 + t279 * t222) * t273) / 0.2e1 + t253 * ((t204 * t282 - t206 * t255 + t208 * t256) * t254 + (t282 * t203 - t255 * t205 + t256 * t207) * t253 + (t220 * t282 - t221 * t255 + t222 * t256) * t273) / 0.2e1 + t254 * ((t284 * t204 - t257 * t206 + t258 * t208) * t254 + (t203 * t284 - t205 * t257 + t207 * t258) * t253 + (t220 * t284 - t221 * t257 + t222 * t258) * t273) / 0.2e1 + t244 * ((t192 * t278 + t194 * t251 + t196 * t252) * t214 + (t191 * t278 + t193 * t251 + t195 * t252) * t213 + (t278 * t199 + t251 * t200 + t252 * t201) * t244) / 0.2e1 + t213 * ((t192 * t255 + t194 * t215 + t196 * t216) * t214 + (t191 * t255 + t215 * t193 + t216 * t195) * t213 + (t199 * t255 + t200 * t215 + t201 * t216) * t244) / 0.2e1 + t214 * ((t192 * t257 + t217 * t194 + t218 * t196) * t214 + (t191 * t257 + t193 * t217 + t195 * t218) * t213 + (t199 * t257 + t200 * t217 + t201 * t218) * t244) / 0.2e1 + m(3) * (t241 ^ 2 + t247 ^ 2 + t248 ^ 2) / 0.2e1 + m(4) * (t188 ^ 2 + t189 ^ 2 + t190 ^ 2) / 0.2e1 + m(5) * (t185 ^ 2 + t186 ^ 2 + t187 ^ 2) / 0.2e1 + m(6) * (t182 ^ 2 + t183 ^ 2 + t184 ^ 2) / 0.2e1 + m(7) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + ((t281 * t403 + t282 * t401 + t306 * t402) * t313 + (t281 * t408 + t282 * t404 + t306 * t406) * t298 + (t409 * t281 + t405 * t282 + t407 * t306) * t297) * t297 / 0.2e1 + ((t283 * t403 + t284 * t401 + t307 * t402) * t313 + (t408 * t283 + t404 * t284 + t406 * t307) * t298 + (t283 * t409 + t405 * t284 + t407 * t307) * t297) * t298 / 0.2e1 + ((t403 * t304 + t401 * t305 + t402 * t322) * t313 + (t304 * t408 + t305 * t404 + t322 * t406) * t298 + (t304 * t409 + t405 * t305 + t407 * t322) * t297) * t313 / 0.2e1 + ((t289 * V_base(5) + t290 * V_base(4) + t316 * t348) * t354 + ((t292 * t353 + t294 * t351) * V_base(4) + (t291 * t353 + t293 * t351) * V_base(5) + (t317 * t353 + t318 * t351) * t348) * t352 + Icges(2,3) * t348 + t336 * V_base(5) + t337 * V_base(4)) * t348 / 0.2e1 + ((t316 * t389 + t317 * t325 + t318 * t326 + t337) * t348 + (t289 * t389 + t291 * t325 + t293 * t326 - t358 * t338 + t340 * t360 + Icges(1,4)) * V_base(5) + (t290 * t389 + t292 * t325 + t294 * t326 - t358 * t339 + t341 * t360 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t316 * t388 + t323 * t317 + t324 * t318 + t336) * t348 + (-t289 * t388 + t323 * t291 + t324 * t293 + t338 * t360 + t358 * t340 + Icges(1,2)) * V_base(5) + (-t290 * t388 + t323 * t292 + t324 * t294 + t339 * t360 + t358 * t341 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
