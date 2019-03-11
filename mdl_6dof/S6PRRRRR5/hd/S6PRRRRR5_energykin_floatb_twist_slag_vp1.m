% Calculate kinetic energy for
% S6PRRRRR5
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
% Datum: 2019-03-09 01:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRRR5_energykin_floatb_twist_slag_vp1: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:03:48
% EndTime: 2019-03-09 01:03:53
% DurationCPUTime: 4.54s
% Computational Cost: add. (5895->447), mult. (15086->673), div. (0->0), fcn. (19332->16), ass. (0->197)
t408 = cos(qJ(3));
t407 = cos(qJ(4));
t359 = cos(pkin(6));
t406 = pkin(8) * t359;
t364 = cos(qJ(5));
t405 = pkin(5) * t364;
t403 = cos(pkin(7));
t402 = sin(pkin(7));
t356 = sin(pkin(13));
t401 = Icges(2,4) * t356;
t358 = cos(pkin(13));
t363 = sin(qJ(2));
t365 = cos(qJ(2));
t393 = t359 * t365;
t322 = -t356 * t363 + t358 * t393;
t394 = t359 * t363;
t323 = t356 * t365 + t358 * t394;
t362 = sin(qJ(3));
t357 = sin(pkin(6));
t383 = t408 * t402;
t382 = t357 * t383;
t384 = t403 * t408;
t283 = -t322 * t384 + t323 * t362 + t358 * t382;
t360 = sin(qJ(5));
t400 = t283 * t360;
t324 = -t356 * t393 - t358 * t363;
t325 = -t356 * t394 + t358 * t365;
t285 = -t324 * t384 + t325 * t362 - t356 * t382;
t399 = t285 * t360;
t395 = t357 * t363;
t308 = -t357 * t365 * t384 - t359 * t383 + t362 * t395;
t398 = t308 * t360;
t397 = t356 * t357;
t396 = t357 * t358;
t392 = qJD(2) * t357;
t391 = V_base(5) * qJ(1) + V_base(1);
t387 = qJD(1) + V_base(3);
t336 = t356 * t392 + V_base(4);
t346 = qJD(2) * t359 + V_base(6);
t386 = t357 * t403;
t385 = t357 * t402;
t378 = -t324 * t402 + t356 * t386;
t300 = qJD(3) * t378 + t336;
t377 = t359 * t403 - t365 * t385;
t310 = qJD(3) * t377 + t346;
t261 = qJD(4) * t285 + t300;
t278 = qJD(4) * t308 + t310;
t335 = -t358 * t392 + V_base(5);
t286 = t325 * t408 + (t324 * t403 + t356 * t385) * t362;
t361 = sin(qJ(4));
t270 = t286 * t361 - t378 * t407;
t226 = qJD(5) * t270 + t261;
t309 = t359 * t402 * t362 + (t362 * t365 * t403 + t363 * t408) * t357;
t287 = t309 * t361 - t377 * t407;
t253 = qJD(5) * t287 + t278;
t379 = -t322 * t402 - t358 * t386;
t299 = qJD(3) * t379 + t335;
t328 = pkin(1) * t356 - pkin(8) * t396;
t381 = -t328 * V_base(6) + t406 * V_base(5) + t391;
t329 = pkin(1) * t358 + pkin(8) * t397;
t380 = t328 * V_base(4) - t329 * V_base(5) + t387;
t260 = qJD(4) * t283 + t299;
t284 = t323 * t408 + (t322 * t403 - t358 * t385) * t362;
t268 = t284 * t361 - t379 * t407;
t225 = qJD(5) * t268 + t260;
t376 = V_base(6) * t329 + V_base(2) + (-qJ(1) - t406) * V_base(4);
t289 = pkin(2) * t323 + pkin(9) * t379;
t311 = pkin(2) * t395 + pkin(9) * t377;
t375 = -t289 * t346 + t311 * t335 + t381;
t290 = pkin(2) * t325 + pkin(9) * t378;
t374 = t289 * t336 - t290 * t335 + t380;
t373 = t290 * t346 - t311 * t336 + t376;
t256 = pkin(3) * t284 + pkin(10) * t283;
t276 = pkin(3) * t309 + pkin(10) * t308;
t372 = -t256 * t310 + t276 * t299 + t375;
t257 = pkin(3) * t286 + pkin(10) * t285;
t371 = t256 * t300 - t257 * t299 + t374;
t370 = t257 * t310 - t276 * t300 + t373;
t269 = t284 * t407 + t361 * t379;
t227 = pkin(4) * t269 + t268 * pkin(11);
t288 = t309 * t407 + t361 * t377;
t258 = pkin(4) * t288 + pkin(11) * t287;
t369 = -t227 * t278 + t258 * t260 + t372;
t271 = t286 * t407 + t361 * t378;
t228 = pkin(4) * t271 + pkin(11) * t270;
t368 = t227 * t261 - t228 * t260 + t371;
t367 = t228 * t278 - t258 * t261 + t370;
t355 = qJ(5) + qJ(6);
t353 = cos(t355);
t352 = sin(t355);
t351 = Icges(2,4) * t358;
t344 = rSges(2,1) * t358 - rSges(2,2) * t356;
t343 = rSges(2,1) * t356 + rSges(2,2) * t358;
t342 = Icges(2,1) * t358 - t401;
t341 = Icges(2,1) * t356 + t351;
t340 = -Icges(2,2) * t356 + t351;
t339 = Icges(2,2) * t358 + t401;
t334 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t333 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t332 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t319 = rSges(3,3) * t359 + (rSges(3,1) * t363 + rSges(3,2) * t365) * t357;
t318 = Icges(3,5) * t359 + (Icges(3,1) * t363 + Icges(3,4) * t365) * t357;
t317 = Icges(3,6) * t359 + (Icges(3,4) * t363 + Icges(3,2) * t365) * t357;
t316 = Icges(3,3) * t359 + (Icges(3,5) * t363 + Icges(3,6) * t365) * t357;
t313 = V_base(5) * rSges(2,3) - t343 * V_base(6) + t391;
t312 = t344 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t307 = t343 * V_base(4) - t344 * V_base(5) + t387;
t298 = rSges(3,1) * t325 + rSges(3,2) * t324 + rSges(3,3) * t397;
t297 = rSges(3,1) * t323 + rSges(3,2) * t322 - rSges(3,3) * t396;
t296 = Icges(3,1) * t325 + Icges(3,4) * t324 + Icges(3,5) * t397;
t295 = Icges(3,1) * t323 + Icges(3,4) * t322 - Icges(3,5) * t396;
t294 = Icges(3,4) * t325 + Icges(3,2) * t324 + Icges(3,6) * t397;
t293 = Icges(3,4) * t323 + Icges(3,2) * t322 - Icges(3,6) * t396;
t292 = Icges(3,5) * t325 + Icges(3,6) * t324 + Icges(3,3) * t397;
t291 = Icges(3,5) * t323 + Icges(3,6) * t322 - Icges(3,3) * t396;
t275 = rSges(4,1) * t309 - rSges(4,2) * t308 + rSges(4,3) * t377;
t274 = Icges(4,1) * t309 - Icges(4,4) * t308 + Icges(4,5) * t377;
t273 = Icges(4,4) * t309 - Icges(4,2) * t308 + Icges(4,6) * t377;
t272 = Icges(4,5) * t309 - Icges(4,6) * t308 + Icges(4,3) * t377;
t267 = t288 * t364 + t398;
t266 = -t288 * t360 + t308 * t364;
t263 = t288 * t353 + t308 * t352;
t262 = -t288 * t352 + t308 * t353;
t255 = -t297 * t346 + t319 * t335 + t381;
t254 = t298 * t346 - t319 * t336 + t376;
t251 = rSges(5,1) * t288 - rSges(5,2) * t287 + rSges(5,3) * t308;
t250 = rSges(4,1) * t286 - rSges(4,2) * t285 + rSges(4,3) * t378;
t249 = rSges(4,1) * t284 - rSges(4,2) * t283 + rSges(4,3) * t379;
t248 = Icges(5,1) * t288 - Icges(5,4) * t287 + Icges(5,5) * t308;
t247 = Icges(5,4) * t288 - Icges(5,2) * t287 + Icges(5,6) * t308;
t246 = Icges(5,5) * t288 - Icges(5,6) * t287 + Icges(5,3) * t308;
t245 = Icges(4,1) * t286 - Icges(4,4) * t285 + Icges(4,5) * t378;
t244 = Icges(4,1) * t284 - Icges(4,4) * t283 + Icges(4,5) * t379;
t243 = Icges(4,4) * t286 - Icges(4,2) * t285 + Icges(4,6) * t378;
t242 = Icges(4,4) * t284 - Icges(4,2) * t283 + Icges(4,6) * t379;
t241 = Icges(4,5) * t286 - Icges(4,6) * t285 + Icges(4,3) * t378;
t240 = Icges(4,5) * t284 - Icges(4,6) * t283 + Icges(4,3) * t379;
t239 = t297 * t336 - t298 * t335 + t380;
t238 = t271 * t364 + t399;
t237 = -t271 * t360 + t285 * t364;
t236 = t269 * t364 + t400;
t235 = -t269 * t360 + t283 * t364;
t233 = t271 * t353 + t285 * t352;
t232 = -t271 * t352 + t285 * t353;
t231 = t269 * t353 + t283 * t352;
t230 = -t269 * t352 + t283 * t353;
t229 = qJD(6) * t287 + t253;
t223 = rSges(6,1) * t267 + rSges(6,2) * t266 + rSges(6,3) * t287;
t222 = rSges(5,1) * t271 - rSges(5,2) * t270 + rSges(5,3) * t285;
t221 = rSges(5,1) * t269 - rSges(5,2) * t268 + rSges(5,3) * t283;
t220 = Icges(5,1) * t271 - Icges(5,4) * t270 + Icges(5,5) * t285;
t219 = Icges(5,1) * t269 - Icges(5,4) * t268 + Icges(5,5) * t283;
t218 = Icges(6,1) * t267 + Icges(6,4) * t266 + Icges(6,5) * t287;
t217 = Icges(5,4) * t271 - Icges(5,2) * t270 + Icges(5,6) * t285;
t216 = Icges(5,4) * t269 - Icges(5,2) * t268 + Icges(5,6) * t283;
t215 = Icges(6,4) * t267 + Icges(6,2) * t266 + Icges(6,6) * t287;
t214 = Icges(5,5) * t271 - Icges(5,6) * t270 + Icges(5,3) * t285;
t213 = Icges(5,5) * t269 - Icges(5,6) * t268 + Icges(5,3) * t283;
t212 = Icges(6,5) * t267 + Icges(6,6) * t266 + Icges(6,3) * t287;
t210 = rSges(7,1) * t263 + rSges(7,2) * t262 + rSges(7,3) * t287;
t209 = Icges(7,1) * t263 + Icges(7,4) * t262 + Icges(7,5) * t287;
t208 = Icges(7,4) * t263 + Icges(7,2) * t262 + Icges(7,6) * t287;
t207 = Icges(7,5) * t263 + Icges(7,6) * t262 + Icges(7,3) * t287;
t206 = pkin(5) * t398 + pkin(12) * t287 + t288 * t405;
t205 = qJD(6) * t270 + t226;
t204 = qJD(6) * t268 + t225;
t202 = rSges(6,1) * t238 + rSges(6,2) * t237 + rSges(6,3) * t270;
t201 = rSges(6,1) * t236 + rSges(6,2) * t235 + rSges(6,3) * t268;
t200 = Icges(6,1) * t238 + Icges(6,4) * t237 + Icges(6,5) * t270;
t199 = Icges(6,1) * t236 + Icges(6,4) * t235 + Icges(6,5) * t268;
t198 = Icges(6,4) * t238 + Icges(6,2) * t237 + Icges(6,6) * t270;
t197 = Icges(6,4) * t236 + Icges(6,2) * t235 + Icges(6,6) * t268;
t196 = Icges(6,5) * t238 + Icges(6,6) * t237 + Icges(6,3) * t270;
t195 = Icges(6,5) * t236 + Icges(6,6) * t235 + Icges(6,3) * t268;
t194 = rSges(7,1) * t233 + rSges(7,2) * t232 + rSges(7,3) * t270;
t193 = rSges(7,1) * t231 + rSges(7,2) * t230 + rSges(7,3) * t268;
t192 = Icges(7,1) * t233 + Icges(7,4) * t232 + Icges(7,5) * t270;
t191 = Icges(7,1) * t231 + Icges(7,4) * t230 + Icges(7,5) * t268;
t190 = Icges(7,4) * t233 + Icges(7,2) * t232 + Icges(7,6) * t270;
t189 = Icges(7,4) * t231 + Icges(7,2) * t230 + Icges(7,6) * t268;
t188 = Icges(7,5) * t233 + Icges(7,6) * t232 + Icges(7,3) * t270;
t187 = Icges(7,5) * t231 + Icges(7,6) * t230 + Icges(7,3) * t268;
t186 = -t249 * t310 + t275 * t299 + t375;
t185 = t250 * t310 - t275 * t300 + t373;
t184 = pkin(5) * t399 + pkin(12) * t270 + t271 * t405;
t183 = pkin(5) * t400 + pkin(12) * t268 + t269 * t405;
t182 = t249 * t300 - t250 * t299 + t374;
t181 = -t221 * t278 + t251 * t260 + t372;
t180 = t222 * t278 - t251 * t261 + t370;
t179 = t221 * t261 - t222 * t260 + t371;
t178 = -t201 * t253 + t223 * t225 + t369;
t177 = t202 * t253 - t223 * t226 + t367;
t176 = t201 * t226 - t202 * t225 + t368;
t175 = -t183 * t253 - t193 * t229 + t204 * t210 + t206 * t225 + t369;
t174 = t184 * t253 + t194 * t229 - t205 * t210 - t206 * t226 + t367;
t173 = t183 * t226 - t184 * t225 + t193 * t205 - t194 * t204 + t368;
t1 = ((t339 * t358 + t341 * t356 + Icges(1,2)) * V_base(5) + (t340 * t358 + t342 * t356 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + t310 * ((t241 * t377 - t308 * t243 + t309 * t245) * t300 + (t240 * t377 - t308 * t242 + t309 * t244) * t299 + (t272 * t377 - t308 * t273 + t309 * t274) * t310) / 0.2e1 + t300 * ((t241 * t378 - t285 * t243 + t286 * t245) * t300 + (t240 * t378 - t285 * t242 + t286 * t244) * t299 + (t272 * t378 - t285 * t273 + t286 * t274) * t310) / 0.2e1 + m(1) * (t332 ^ 2 + t333 ^ 2 + t334 ^ 2) / 0.2e1 + m(2) * (t307 ^ 2 + t312 ^ 2 + t313 ^ 2) / 0.2e1 + t278 * ((t214 * t308 - t217 * t287 + t220 * t288) * t261 + (t213 * t308 - t216 * t287 + t219 * t288) * t260 + (t246 * t308 - t247 * t287 + t248 * t288) * t278) / 0.2e1 + t229 * ((t188 * t287 + t190 * t262 + t192 * t263) * t205 + (t187 * t287 + t189 * t262 + t191 * t263) * t204 + (t287 * t207 + t262 * t208 + t263 * t209) * t229) / 0.2e1 + t253 * ((t196 * t287 + t198 * t266 + t200 * t267) * t226 + (t195 * t287 + t197 * t266 + t199 * t267) * t225 + (t212 * t287 + t215 * t266 + t218 * t267) * t253) / 0.2e1 + t261 * ((t214 * t285 - t217 * t270 + t220 * t271) * t261 + (t213 * t285 - t216 * t270 + t219 * t271) * t260 + (t246 * t285 - t247 * t270 + t248 * t271) * t278) / 0.2e1 + t260 * ((t214 * t283 - t217 * t268 + t220 * t269) * t261 + (t213 * t283 - t216 * t268 + t219 * t269) * t260 + (t246 * t283 - t247 * t268 + t248 * t269) * t278) / 0.2e1 + t226 * ((t270 * t196 + t237 * t198 + t238 * t200) * t226 + (t195 * t270 + t197 * t237 + t199 * t238) * t225 + (t212 * t270 + t215 * t237 + t218 * t238) * t253) / 0.2e1 + t205 * ((t270 * t188 + t232 * t190 + t233 * t192) * t205 + (t187 * t270 + t189 * t232 + t191 * t233) * t204 + (t207 * t270 + t208 * t232 + t209 * t233) * t229) / 0.2e1 + t225 * ((t196 * t268 + t198 * t235 + t200 * t236) * t226 + (t268 * t195 + t235 * t197 + t236 * t199) * t225 + (t212 * t268 + t215 * t235 + t218 * t236) * t253) / 0.2e1 + t204 * ((t188 * t268 + t190 * t230 + t192 * t231) * t205 + (t268 * t187 + t230 * t189 + t231 * t191) * t204 + (t207 * t268 + t208 * t230 + t209 * t231) * t229) / 0.2e1 + m(3) * (t239 ^ 2 + t254 ^ 2 + t255 ^ 2) / 0.2e1 + ((Icges(2,5) * t358 - Icges(2,6) * t356 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t356 + Icges(2,6) * t358 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + ((-t339 * t356 + t341 * t358 + Icges(1,4)) * V_base(5) + (-t340 * t356 + t342 * t358 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + m(4) * (t182 ^ 2 + t185 ^ 2 + t186 ^ 2) / 0.2e1 + m(5) * (t179 ^ 2 + t180 ^ 2 + t181 ^ 2) / 0.2e1 + m(6) * (t176 ^ 2 + t177 ^ 2 + t178 ^ 2) / 0.2e1 + m(7) * (t173 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + t299 * ((t241 * t379 - t283 * t243 + t284 * t245) * t300 + (t240 * t379 - t283 * t242 + t284 * t244) * t299 + (t272 * t379 - t283 * t273 + t284 * t274) * t310) / 0.2e1 + t335 * ((-t292 * t396 + t294 * t322 + t296 * t323) * t336 + (-t291 * t396 + t293 * t322 + t295 * t323) * t335 + (-t316 * t396 + t317 * t322 + t318 * t323) * t346) / 0.2e1 + t336 * ((t292 * t397 + t294 * t324 + t296 * t325) * t336 + (t291 * t397 + t293 * t324 + t295 * t325) * t335 + (t316 * t397 + t317 * t324 + t318 * t325) * t346) / 0.2e1 + t346 * ((t291 * t335 + t292 * t336 + t316 * t346) * t359 + ((t294 * t365 + t296 * t363) * t336 + (t293 * t365 + t295 * t363) * t335 + (t317 * t365 + t318 * t363) * t346) * t357) / 0.2e1;
T  = t1;
