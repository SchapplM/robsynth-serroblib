% Calculate kinetic energy for
% S6PRPRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRP2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:59:51
% EndTime: 2019-03-08 19:59:55
% DurationCPUTime: 4.60s
% Computational Cost: add. (3714->383), mult. (8957->542), div. (0->0), fcn. (11275->12), ass. (0->175)
t399 = Icges(6,1) + Icges(7,1);
t398 = -Icges(6,4) + Icges(7,5);
t397 = Icges(7,4) + Icges(6,5);
t396 = Icges(6,2) + Icges(7,3);
t395 = Icges(7,2) + Icges(6,3);
t394 = -Icges(6,6) + Icges(7,6);
t393 = rSges(7,1) + pkin(5);
t392 = rSges(7,3) + qJ(6);
t325 = sin(pkin(11));
t332 = sin(qJ(2));
t333 = cos(qJ(2));
t367 = cos(pkin(11));
t295 = -t333 * t325 - t332 * t367;
t329 = cos(pkin(6));
t287 = t295 * t329;
t296 = -t332 * t325 + t333 * t367;
t326 = sin(pkin(10));
t328 = cos(pkin(10));
t268 = -t287 * t328 + t296 * t326;
t327 = sin(pkin(6));
t331 = sin(qJ(4));
t363 = t327 * t331;
t371 = cos(qJ(4));
t244 = t268 * t371 - t328 * t363;
t344 = t329 * t296;
t267 = t326 * t295 + t328 * t344;
t330 = sin(qJ(5));
t370 = cos(qJ(5));
t215 = t244 * t330 + t267 * t370;
t216 = t244 * t370 - t267 * t330;
t349 = t327 * t371;
t243 = t268 * t331 + t328 * t349;
t390 = t215 * t396 + t216 * t398 + t243 * t394;
t270 = t287 * t326 + t296 * t328;
t246 = t270 * t371 + t326 * t363;
t269 = t295 * t328 - t326 * t344;
t217 = t246 * t330 + t269 * t370;
t218 = t246 * t370 - t269 * t330;
t245 = t270 * t331 - t326 * t349;
t389 = t217 * t396 + t218 * t398 + t245 * t394;
t388 = t215 * t394 + t216 * t397 + t243 * t395;
t387 = t217 * t394 + t218 * t397 + t245 * t395;
t386 = t398 * t215 + t216 * t399 + t397 * t243;
t385 = t398 * t217 + t218 * t399 + t397 * t245;
t286 = t295 * t327;
t274 = -t286 * t371 + t329 * t331;
t285 = t296 * t327;
t237 = t274 * t330 + t285 * t370;
t238 = t274 * t370 - t285 * t330;
t273 = -t286 * t331 - t329 * t371;
t384 = t237 * t396 + t238 * t398 + t273 * t394;
t383 = t237 * t394 + t238 * t397 + t273 * t395;
t382 = t398 * t237 + t238 * t399 + t397 * t273;
t364 = t327 * t328;
t221 = Icges(4,5) * t268 + Icges(4,6) * t267 - Icges(4,3) * t364;
t361 = t329 * t333;
t289 = -t326 * t332 + t328 * t361;
t362 = t329 * t332;
t290 = t326 * t333 + t328 * t362;
t252 = Icges(3,5) * t290 + Icges(3,6) * t289 - Icges(3,3) * t364;
t381 = t221 + t252;
t365 = t326 * t327;
t222 = Icges(4,5) * t270 + Icges(4,6) * t269 + Icges(4,3) * t365;
t291 = -t326 * t361 - t328 * t332;
t292 = -t326 * t362 + t328 * t333;
t253 = Icges(3,5) * t292 + Icges(3,6) * t291 + Icges(3,3) * t365;
t380 = t222 + t253;
t248 = -Icges(4,5) * t286 + Icges(4,6) * t285 + Icges(4,3) * t329;
t281 = Icges(3,3) * t329 + (Icges(3,5) * t332 + Icges(3,6) * t333) * t327;
t379 = t248 + t281;
t369 = pkin(7) * t329;
t368 = pkin(2) * t333;
t366 = Icges(2,4) * t326;
t359 = rSges(7,2) * t243 + t392 * t215 + t393 * t216;
t358 = rSges(7,2) * t245 + t392 * t217 + t393 * t218;
t357 = rSges(7,2) * t273 + t392 * t237 + t393 * t238;
t356 = qJD(2) * t327;
t355 = qJD(3) * t327;
t354 = V_base(5) * qJ(1) + V_base(1);
t350 = qJD(1) + V_base(3);
t304 = t326 * t356 + V_base(4);
t315 = qJD(2) * t329 + V_base(6);
t348 = pkin(2) * t362 - qJ(3) * t327;
t242 = -qJD(4) * t269 + t304;
t272 = -qJD(4) * t285 + t315;
t303 = -t328 * t356 + V_base(5);
t241 = -qJD(4) * t267 + t303;
t298 = pkin(1) * t326 - pkin(7) * t364;
t346 = -t298 * V_base(6) + V_base(5) * t369 + t354;
t299 = pkin(1) * t328 + pkin(7) * t365;
t345 = V_base(4) * t298 - t299 * V_base(5) + t350;
t343 = V_base(6) * t299 + V_base(2) + (-qJ(1) - t369) * V_base(4);
t297 = pkin(2) * t327 * t332 + qJ(3) * t329;
t342 = t303 * t297 + t326 * t355 + t346;
t262 = t326 * t368 + t328 * t348;
t341 = qJD(3) * t329 + t304 * t262 + t345;
t263 = -t326 * t348 + t328 * t368;
t340 = t315 * t263 - t328 * t355 + t343;
t229 = pkin(3) * t268 - pkin(8) * t267;
t266 = -pkin(3) * t286 - pkin(8) * t285;
t339 = t303 * t266 + (-t229 - t262) * t315 + t342;
t230 = pkin(3) * t270 - pkin(8) * t269;
t338 = t304 * t229 + (-t230 - t263) * t303 + t341;
t211 = pkin(4) * t244 + pkin(9) * t243;
t235 = pkin(4) * t274 + pkin(9) * t273;
t337 = -t211 * t272 + t241 * t235 + t339;
t212 = pkin(4) * t246 + pkin(9) * t245;
t336 = t242 * t211 - t212 * t241 + t338;
t335 = t315 * t230 + (-t266 - t297) * t304 + t340;
t334 = t272 * t212 - t235 * t242 + t335;
t323 = Icges(2,4) * t328;
t312 = rSges(2,1) * t328 - rSges(2,2) * t326;
t311 = rSges(2,1) * t326 + rSges(2,2) * t328;
t310 = Icges(2,1) * t328 - t366;
t309 = Icges(2,1) * t326 + t323;
t308 = -Icges(2,2) * t326 + t323;
t307 = Icges(2,2) * t328 + t366;
t302 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t301 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t300 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t284 = t329 * rSges(3,3) + (rSges(3,1) * t332 + rSges(3,2) * t333) * t327;
t283 = Icges(3,5) * t329 + (Icges(3,1) * t332 + Icges(3,4) * t333) * t327;
t282 = Icges(3,6) * t329 + (Icges(3,4) * t332 + Icges(3,2) * t333) * t327;
t276 = V_base(5) * rSges(2,3) - t311 * V_base(6) + t354;
t275 = t312 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t271 = t311 * V_base(4) - t312 * V_base(5) + t350;
t259 = rSges(3,1) * t292 + rSges(3,2) * t291 + rSges(3,3) * t365;
t258 = rSges(3,1) * t290 + rSges(3,2) * t289 - rSges(3,3) * t364;
t257 = Icges(3,1) * t292 + Icges(3,4) * t291 + Icges(3,5) * t365;
t256 = Icges(3,1) * t290 + Icges(3,4) * t289 - Icges(3,5) * t364;
t255 = Icges(3,4) * t292 + Icges(3,2) * t291 + Icges(3,6) * t365;
t254 = Icges(3,4) * t290 + Icges(3,2) * t289 - Icges(3,6) * t364;
t251 = -rSges(4,1) * t286 + rSges(4,2) * t285 + rSges(4,3) * t329;
t250 = -Icges(4,1) * t286 + Icges(4,4) * t285 + Icges(4,5) * t329;
t249 = -Icges(4,4) * t286 + Icges(4,2) * t285 + Icges(4,6) * t329;
t236 = qJD(5) * t273 + t272;
t234 = rSges(5,1) * t274 - rSges(5,2) * t273 - rSges(5,3) * t285;
t233 = Icges(5,1) * t274 - Icges(5,4) * t273 - Icges(5,5) * t285;
t232 = Icges(5,4) * t274 - Icges(5,2) * t273 - Icges(5,6) * t285;
t231 = Icges(5,5) * t274 - Icges(5,6) * t273 - Icges(5,3) * t285;
t228 = rSges(4,1) * t270 + rSges(4,2) * t269 + rSges(4,3) * t365;
t227 = rSges(4,1) * t268 + rSges(4,2) * t267 - rSges(4,3) * t364;
t226 = Icges(4,1) * t270 + Icges(4,4) * t269 + Icges(4,5) * t365;
t225 = Icges(4,1) * t268 + Icges(4,4) * t267 - Icges(4,5) * t364;
t224 = Icges(4,4) * t270 + Icges(4,2) * t269 + Icges(4,6) * t365;
t223 = Icges(4,4) * t268 + Icges(4,2) * t267 - Icges(4,6) * t364;
t214 = qJD(5) * t245 + t242;
t213 = qJD(5) * t243 + t241;
t210 = -t258 * t315 + t284 * t303 + t346;
t209 = t259 * t315 - t284 * t304 + t343;
t205 = rSges(6,1) * t238 - rSges(6,2) * t237 + rSges(6,3) * t273;
t197 = t258 * t304 - t259 * t303 + t345;
t196 = rSges(5,1) * t246 - rSges(5,2) * t245 - rSges(5,3) * t269;
t195 = rSges(5,1) * t244 - rSges(5,2) * t243 - rSges(5,3) * t267;
t194 = Icges(5,1) * t246 - Icges(5,4) * t245 - Icges(5,5) * t269;
t193 = Icges(5,1) * t244 - Icges(5,4) * t243 - Icges(5,5) * t267;
t192 = Icges(5,4) * t246 - Icges(5,2) * t245 - Icges(5,6) * t269;
t191 = Icges(5,4) * t244 - Icges(5,2) * t243 - Icges(5,6) * t267;
t190 = Icges(5,5) * t246 - Icges(5,6) * t245 - Icges(5,3) * t269;
t189 = Icges(5,5) * t244 - Icges(5,6) * t243 - Icges(5,3) * t267;
t185 = rSges(6,1) * t218 - rSges(6,2) * t217 + rSges(6,3) * t245;
t183 = rSges(6,1) * t216 - rSges(6,2) * t215 + rSges(6,3) * t243;
t169 = t251 * t303 + (-t227 - t262) * t315 + t342;
t168 = t228 * t315 + (-t251 - t297) * t304 + t340;
t167 = t227 * t304 + (-t228 - t263) * t303 + t341;
t166 = -t195 * t272 + t234 * t241 + t339;
t165 = t196 * t272 - t234 * t242 + t335;
t164 = t195 * t242 - t196 * t241 + t338;
t163 = -t183 * t236 + t205 * t213 + t337;
t162 = t185 * t236 - t205 * t214 + t334;
t161 = t183 * t214 - t185 * t213 + t336;
t160 = qJD(6) * t217 + t213 * t357 - t236 * t359 + t337;
t159 = qJD(6) * t215 - t214 * t357 + t236 * t358 + t334;
t158 = qJD(6) * t237 - t213 * t358 + t214 * t359 + t336;
t1 = m(7) * (t158 ^ 2 + t159 ^ 2 + t160 ^ 2) / 0.2e1 + m(6) * (t161 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(5) * (t164 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(4) * (t167 ^ 2 + t168 ^ 2 + t169 ^ 2) / 0.2e1 + m(3) * (t197 ^ 2 + t209 ^ 2 + t210 ^ 2) / 0.2e1 + t241 * ((-t190 * t267 - t192 * t243 + t194 * t244) * t242 + (-t189 * t267 - t191 * t243 + t193 * t244) * t241 + (-t231 * t267 - t232 * t243 + t233 * t244) * t272) / 0.2e1 + t242 * ((-t190 * t269 - t192 * t245 + t194 * t246) * t242 + (-t189 * t269 - t191 * t245 + t193 * t246) * t241 + (-t231 * t269 - t232 * t245 + t233 * t246) * t272) / 0.2e1 + m(2) * (t271 ^ 2 + t275 ^ 2 + t276 ^ 2) / 0.2e1 + t272 * ((-t190 * t285 - t192 * t273 + t194 * t274) * t242 + (-t189 * t285 - t191 * t273 + t193 * t274) * t241 + (-t231 * t285 - t232 * t273 + t233 * t274) * t272) / 0.2e1 + m(1) * (t300 ^ 2 + t301 ^ 2 + t302 ^ 2) / 0.2e1 + ((t215 * t384 + t216 * t382 + t243 * t383) * t236 + (t215 * t389 + t216 * t385 + t243 * t387) * t214 + (t390 * t215 + t386 * t216 + t388 * t243) * t213) * t213 / 0.2e1 + ((t217 * t384 + t218 * t382 + t245 * t383) * t236 + (t389 * t217 + t385 * t218 + t387 * t245) * t214 + (t217 * t390 + t218 * t386 + t245 * t388) * t213) * t214 / 0.2e1 + ((t384 * t237 + t382 * t238 + t383 * t273) * t236 + (t237 * t389 + t238 * t385 + t273 * t387) * t214 + (t390 * t237 + t238 * t386 + t388 * t273) * t213) * t236 / 0.2e1 + ((t249 * t267 + t250 * t268 + t282 * t289 + t283 * t290 - t364 * t379) * t315 + (t224 * t267 + t226 * t268 + t255 * t289 + t257 * t290 - t364 * t380) * t304 + (t223 * t267 + t225 * t268 + t254 * t289 + t256 * t290 - t381 * t364) * t303) * t303 / 0.2e1 + ((t249 * t269 + t250 * t270 + t282 * t291 + t283 * t292 + t365 * t379) * t315 + (t224 * t269 + t226 * t270 + t255 * t291 + t257 * t292 + t380 * t365) * t304 + (t223 * t269 + t225 * t270 + t254 * t291 + t256 * t292 + t365 * t381) * t303) * t304 / 0.2e1 + ((t252 * t303 + t253 * t304 + t281 * t315) * t329 + ((t255 * t333 + t257 * t332) * t304 + (t254 * t333 + t256 * t332) * t303 + (t282 * t333 + t283 * t332) * t315) * t327 + (t222 * t329 + t224 * t285 - t226 * t286) * t304 + (t221 * t329 + t223 * t285 - t225 * t286) * t303 + (t248 * t329 + t249 * t285 - t250 * t286) * t315) * t315 / 0.2e1 + ((-t307 * t326 + t309 * t328 + Icges(1,4)) * V_base(5) + (-t308 * t326 + t310 * t328 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t307 * t328 + t309 * t326 + Icges(1,2)) * V_base(5) + (t308 * t328 + t310 * t326 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t328 - Icges(2,6) * t326 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t326 + Icges(2,6) * t328 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
