% Calculate kinetic energy for
% S6RRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP6_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:07:40
% EndTime: 2019-03-09 12:07:44
% DurationCPUTime: 4.49s
% Computational Cost: add. (3774->383), mult. (8957->543), div. (0->0), fcn. (11275->12), ass. (0->177)
t402 = Icges(6,1) + Icges(7,1);
t401 = -Icges(6,4) + Icges(7,5);
t400 = Icges(7,4) + Icges(6,5);
t399 = Icges(6,2) + Icges(7,3);
t398 = Icges(7,2) + Icges(6,3);
t397 = -Icges(6,6) + Icges(7,6);
t396 = rSges(7,1) + pkin(5);
t395 = rSges(7,3) + qJ(6);
t330 = sin(pkin(11));
t335 = sin(qJ(2));
t337 = cos(qJ(2));
t371 = cos(pkin(11));
t299 = -t337 * t330 - t335 * t371;
t332 = cos(pkin(6));
t291 = t299 * t332;
t300 = -t335 * t330 + t337 * t371;
t336 = sin(qJ(1));
t338 = cos(qJ(1));
t272 = -t291 * t338 + t336 * t300;
t334 = sin(qJ(4));
t331 = sin(pkin(6));
t368 = t331 * t338;
t376 = cos(qJ(4));
t248 = t272 * t376 - t334 * t368;
t349 = t332 * t300;
t271 = t336 * t299 + t338 * t349;
t333 = sin(qJ(5));
t375 = cos(qJ(5));
t219 = t248 * t333 + t271 * t375;
t220 = t248 * t375 - t271 * t333;
t354 = t331 * t376;
t247 = t272 * t334 + t338 * t354;
t394 = t399 * t219 + t401 * t220 + t397 * t247;
t274 = t336 * t291 + t300 * t338;
t369 = t331 * t336;
t250 = t274 * t376 + t334 * t369;
t273 = t299 * t338 - t336 * t349;
t221 = t250 * t333 + t273 * t375;
t222 = t250 * t375 - t273 * t333;
t249 = t274 * t334 - t336 * t354;
t393 = t399 * t221 + t401 * t222 + t397 * t249;
t392 = t397 * t219 + t400 * t220 + t398 * t247;
t391 = t397 * t221 + t400 * t222 + t398 * t249;
t390 = t401 * t219 + t402 * t220 + t400 * t247;
t389 = t401 * t221 + t402 * t222 + t400 * t249;
t290 = t299 * t331;
t277 = -t290 * t376 + t332 * t334;
t289 = t300 * t331;
t241 = t277 * t333 + t289 * t375;
t242 = t277 * t375 - t289 * t333;
t276 = -t290 * t334 - t332 * t376;
t388 = t399 * t241 + t401 * t242 + t397 * t276;
t387 = t397 * t241 + t400 * t242 + t398 * t276;
t386 = t401 * t241 + t402 * t242 + t400 * t276;
t225 = Icges(4,5) * t272 + Icges(4,6) * t271 - Icges(4,3) * t368;
t363 = t337 * t338;
t365 = t336 * t335;
t294 = t332 * t363 - t365;
t364 = t336 * t337;
t366 = t335 * t338;
t295 = t332 * t366 + t364;
t256 = Icges(3,5) * t295 + Icges(3,6) * t294 - Icges(3,3) * t368;
t385 = t225 + t256;
t226 = Icges(4,5) * t274 + Icges(4,6) * t273 + Icges(4,3) * t369;
t296 = -t332 * t364 - t366;
t297 = -t332 * t365 + t363;
t257 = Icges(3,5) * t297 + Icges(3,6) * t296 + Icges(3,3) * t369;
t384 = t226 + t257;
t252 = -Icges(4,5) * t290 + Icges(4,6) * t289 + Icges(4,3) * t332;
t285 = Icges(3,3) * t332 + (Icges(3,5) * t335 + Icges(3,6) * t337) * t331;
t383 = t252 + t285;
t374 = pkin(2) * t335;
t373 = pkin(8) * t332;
t372 = pkin(2) * t337;
t370 = Icges(2,4) * t336;
t362 = rSges(7,2) * t247 + t395 * t219 + t396 * t220;
t361 = rSges(7,2) * t249 + t395 * t221 + t396 * t222;
t360 = rSges(7,2) * t276 + t395 * t241 + t396 * t242;
t359 = qJD(2) * t331;
t358 = qJD(3) * t331;
t357 = V_base(5) * pkin(7) + V_base(1);
t308 = t336 * t359 + V_base(4);
t327 = V_base(6) + qJD(1);
t353 = -qJ(3) * t331 + t332 * t374;
t246 = -qJD(4) * t273 + t308;
t309 = qJD(2) * t332 + t327;
t275 = -qJD(4) * t289 + t309;
t307 = -t338 * t359 + V_base(5);
t302 = t336 * pkin(1) - pkin(8) * t368;
t351 = -t302 * t327 + V_base(5) * t373 + t357;
t303 = pkin(1) * t338 + pkin(8) * t369;
t350 = V_base(4) * t302 - t303 * V_base(5) + V_base(3);
t245 = -qJD(4) * t271 + t307;
t301 = qJ(3) * t332 + t331 * t374;
t348 = t307 * t301 + t336 * t358 + t351;
t269 = t336 * t372 + t338 * t353;
t347 = qJD(3) * t332 + t308 * t269 + t350;
t346 = t327 * t303 + V_base(2) + (-pkin(7) - t373) * V_base(4);
t237 = pkin(3) * t272 - pkin(9) * t271;
t266 = -pkin(3) * t290 - pkin(9) * t289;
t345 = t307 * t266 + (-t237 - t269) * t309 + t348;
t238 = pkin(3) * t274 - pkin(9) * t273;
t270 = -t336 * t353 + t338 * t372;
t344 = t308 * t237 + (-t238 - t270) * t307 + t347;
t343 = t309 * t270 - t338 * t358 + t346;
t215 = pkin(4) * t248 + pkin(10) * t247;
t239 = pkin(4) * t277 + pkin(10) * t276;
t342 = -t215 * t275 + t245 * t239 + t345;
t216 = pkin(4) * t250 + pkin(10) * t249;
t341 = t246 * t215 - t216 * t245 + t344;
t340 = t309 * t238 + (-t266 - t301) * t308 + t343;
t339 = t275 * t216 - t246 * t239 + t340;
t328 = Icges(2,4) * t338;
t317 = rSges(2,1) * t338 - t336 * rSges(2,2);
t316 = t336 * rSges(2,1) + rSges(2,2) * t338;
t315 = Icges(2,1) * t338 - t370;
t314 = Icges(2,1) * t336 + t328;
t313 = -Icges(2,2) * t336 + t328;
t312 = Icges(2,2) * t338 + t370;
t306 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t305 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t304 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t288 = rSges(3,3) * t332 + (rSges(3,1) * t335 + rSges(3,2) * t337) * t331;
t287 = Icges(3,5) * t332 + (Icges(3,1) * t335 + Icges(3,4) * t337) * t331;
t286 = Icges(3,6) * t332 + (Icges(3,4) * t335 + Icges(3,2) * t337) * t331;
t280 = V_base(5) * rSges(2,3) - t316 * t327 + t357;
t279 = t317 * t327 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t278 = t316 * V_base(4) - t317 * V_base(5) + V_base(3);
t265 = rSges(3,1) * t297 + rSges(3,2) * t296 + rSges(3,3) * t369;
t264 = t295 * rSges(3,1) + t294 * rSges(3,2) - rSges(3,3) * t368;
t261 = Icges(3,1) * t297 + Icges(3,4) * t296 + Icges(3,5) * t369;
t260 = Icges(3,1) * t295 + Icges(3,4) * t294 - Icges(3,5) * t368;
t259 = Icges(3,4) * t297 + Icges(3,2) * t296 + Icges(3,6) * t369;
t258 = Icges(3,4) * t295 + Icges(3,2) * t294 - Icges(3,6) * t368;
t255 = -rSges(4,1) * t290 + rSges(4,2) * t289 + rSges(4,3) * t332;
t254 = -Icges(4,1) * t290 + Icges(4,4) * t289 + Icges(4,5) * t332;
t253 = -Icges(4,4) * t290 + Icges(4,2) * t289 + Icges(4,6) * t332;
t240 = qJD(5) * t276 + t275;
t236 = rSges(5,1) * t277 - rSges(5,2) * t276 - rSges(5,3) * t289;
t235 = Icges(5,1) * t277 - Icges(5,4) * t276 - Icges(5,5) * t289;
t234 = Icges(5,4) * t277 - Icges(5,2) * t276 - Icges(5,6) * t289;
t233 = Icges(5,5) * t277 - Icges(5,6) * t276 - Icges(5,3) * t289;
t232 = rSges(4,1) * t274 + rSges(4,2) * t273 + rSges(4,3) * t369;
t231 = t272 * rSges(4,1) + t271 * rSges(4,2) - rSges(4,3) * t368;
t230 = Icges(4,1) * t274 + Icges(4,4) * t273 + Icges(4,5) * t369;
t229 = Icges(4,1) * t272 + Icges(4,4) * t271 - Icges(4,5) * t368;
t228 = Icges(4,4) * t274 + Icges(4,2) * t273 + Icges(4,6) * t369;
t227 = Icges(4,4) * t272 + Icges(4,2) * t271 - Icges(4,6) * t368;
t218 = qJD(5) * t249 + t246;
t217 = qJD(5) * t247 + t245;
t213 = -t264 * t309 + t288 * t307 + t351;
t212 = t265 * t309 - t288 * t308 + t346;
t209 = t264 * t308 - t265 * t307 + t350;
t208 = rSges(6,1) * t242 - rSges(6,2) * t241 + rSges(6,3) * t276;
t200 = rSges(5,1) * t250 - rSges(5,2) * t249 - rSges(5,3) * t273;
t199 = rSges(5,1) * t248 - rSges(5,2) * t247 - rSges(5,3) * t271;
t198 = Icges(5,1) * t250 - Icges(5,4) * t249 - Icges(5,5) * t273;
t197 = Icges(5,1) * t248 - Icges(5,4) * t247 - Icges(5,5) * t271;
t196 = Icges(5,4) * t250 - Icges(5,2) * t249 - Icges(5,6) * t273;
t195 = Icges(5,4) * t248 - Icges(5,2) * t247 - Icges(5,6) * t271;
t194 = Icges(5,5) * t250 - Icges(5,6) * t249 - Icges(5,3) * t273;
t193 = Icges(5,5) * t248 - Icges(5,6) * t247 - Icges(5,3) * t271;
t189 = rSges(6,1) * t222 - rSges(6,2) * t221 + rSges(6,3) * t249;
t187 = rSges(6,1) * t220 - rSges(6,2) * t219 + rSges(6,3) * t247;
t173 = t255 * t307 + (-t231 - t269) * t309 + t348;
t172 = t309 * t232 + (-t255 - t301) * t308 + t343;
t171 = t231 * t308 + (-t232 - t270) * t307 + t347;
t170 = -t199 * t275 + t236 * t245 + t345;
t169 = t275 * t200 - t246 * t236 + t340;
t168 = t199 * t246 - t200 * t245 + t344;
t167 = -t187 * t240 + t208 * t217 + t342;
t166 = t240 * t189 - t218 * t208 + t339;
t165 = t187 * t218 - t189 * t217 + t341;
t164 = qJD(6) * t221 + t217 * t360 - t240 * t362 + t342;
t163 = qJD(6) * t219 - t218 * t360 + t240 * t361 + t339;
t162 = qJD(6) * t241 - t217 * t361 + t218 * t362 + t341;
t1 = m(7) * (t162 ^ 2 + t163 ^ 2 + t164 ^ 2) / 0.2e1 + m(6) * (t165 ^ 2 + t166 ^ 2 + t167 ^ 2) / 0.2e1 + m(5) * (t168 ^ 2 + t169 ^ 2 + t170 ^ 2) / 0.2e1 + m(4) * (t171 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + m(3) * (t209 ^ 2 + t212 ^ 2 + t213 ^ 2) / 0.2e1 + t245 * ((-t194 * t271 - t196 * t247 + t198 * t248) * t246 + (-t193 * t271 - t195 * t247 + t197 * t248) * t245 + (-t233 * t271 - t234 * t247 + t235 * t248) * t275) / 0.2e1 + t246 * ((-t194 * t273 - t196 * t249 + t198 * t250) * t246 + (-t193 * t273 - t195 * t249 + t197 * t250) * t245 + (-t233 * t273 - t234 * t249 + t235 * t250) * t275) / 0.2e1 + m(2) * (t278 ^ 2 + t279 ^ 2 + t280 ^ 2) / 0.2e1 + t275 * ((-t194 * t289 - t196 * t276 + t198 * t277) * t246 + (-t193 * t289 - t195 * t276 + t197 * t277) * t245 + (-t233 * t289 - t234 * t276 + t235 * t277) * t275) / 0.2e1 + m(1) * (t304 ^ 2 + t305 ^ 2 + t306 ^ 2) / 0.2e1 + ((t219 * t388 + t220 * t386 + t247 * t387) * t240 + (t219 * t393 + t220 * t389 + t247 * t391) * t218 + (t394 * t219 + t390 * t220 + t392 * t247) * t217) * t217 / 0.2e1 + ((t388 * t221 + t222 * t386 + t387 * t249) * t240 + (t393 * t221 + t389 * t222 + t391 * t249) * t218 + (t221 * t394 + t390 * t222 + t392 * t249) * t217) * t218 / 0.2e1 + ((t388 * t241 + t386 * t242 + t387 * t276) * t240 + (t241 * t393 + t242 * t389 + t276 * t391) * t218 + (t241 * t394 + t390 * t242 + t392 * t276) * t217) * t240 / 0.2e1 + ((t271 * t253 + t272 * t254 + t294 * t286 + t295 * t287 - t368 * t383) * t309 + (t271 * t228 + t272 * t230 + t294 * t259 + t295 * t261 - t368 * t384) * t308 + (t271 * t227 + t272 * t229 + t294 * t258 + t295 * t260 - t385 * t368) * t307) * t307 / 0.2e1 + ((t253 * t273 + t254 * t274 + t286 * t296 + t287 * t297 + t369 * t383) * t309 + (t228 * t273 + t230 * t274 + t259 * t296 + t261 * t297 + t384 * t369) * t308 + (t227 * t273 + t229 * t274 + t258 * t296 + t260 * t297 + t369 * t385) * t307) * t308 / 0.2e1 + ((t256 * t307 + t257 * t308 + t285 * t309) * t332 + ((t259 * t337 + t261 * t335) * t308 + (t258 * t337 + t260 * t335) * t307 + (t286 * t337 + t287 * t335) * t309) * t331 + (t226 * t332 + t228 * t289 - t230 * t290) * t308 + (t225 * t332 + t227 * t289 - t229 * t290) * t307 + (t252 * t332 + t253 * t289 - t254 * t290) * t309) * t309 / 0.2e1 + ((-t336 * t312 + t314 * t338 + Icges(1,4)) * V_base(5) + (-t336 * t313 + t315 * t338 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t312 * t338 + t336 * t314 + Icges(1,2)) * V_base(5) + (t313 * t338 + t336 * t315 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t336 + Icges(2,6) * t338) * V_base(5) + (Icges(2,5) * t338 - Icges(2,6) * t336) * V_base(4) + Icges(2,3) * t327 / 0.2e1) * t327;
T  = t1;
