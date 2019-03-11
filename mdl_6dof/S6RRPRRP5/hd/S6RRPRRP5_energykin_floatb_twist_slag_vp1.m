% Calculate kinetic energy for
% S6RRPRRP5
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
% Datum: 2019-03-09 12:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRP5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP5_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRP5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRP5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:57:41
% EndTime: 2019-03-09 11:57:46
% DurationCPUTime: 4.64s
% Computational Cost: add. (3808->388), mult. (9009->550), div. (0->0), fcn. (11331->12), ass. (0->180)
t389 = Icges(6,1) + Icges(7,1);
t388 = Icges(6,4) + Icges(7,4);
t387 = Icges(6,5) + Icges(7,5);
t386 = Icges(6,2) + Icges(7,2);
t385 = Icges(6,6) + Icges(7,6);
t384 = Icges(6,3) + Icges(7,3);
t383 = rSges(7,3) + qJ(6);
t312 = cos(pkin(6));
t316 = sin(qJ(2));
t319 = cos(qJ(2));
t357 = sin(pkin(11));
t358 = cos(pkin(11));
t332 = t316 * t358 + t319 * t357;
t271 = t332 * t312;
t280 = -t316 * t357 + t319 * t358;
t317 = sin(qJ(1));
t320 = cos(qJ(1));
t253 = t271 * t320 + t280 * t317;
t315 = sin(qJ(4));
t311 = sin(pkin(6));
t351 = t311 * t320;
t364 = cos(qJ(4));
t231 = t253 * t364 - t315 * t351;
t328 = t312 * t280;
t252 = -t317 * t332 + t320 * t328;
t314 = sin(qJ(5));
t318 = cos(qJ(5));
t202 = -t231 * t314 - t252 * t318;
t355 = t252 * t314;
t203 = t231 * t318 - t355;
t338 = t311 * t364;
t230 = t253 * t315 + t320 * t338;
t382 = t202 * t385 + t203 * t387 + t230 * t384;
t255 = -t271 * t317 + t280 * t320;
t352 = t311 * t317;
t233 = t255 * t364 + t315 * t352;
t254 = -t317 * t328 - t320 * t332;
t204 = -t233 * t314 - t254 * t318;
t354 = t254 * t314;
t205 = t233 * t318 - t354;
t232 = t255 * t315 - t317 * t338;
t381 = t204 * t385 + t205 * t387 + t232 * t384;
t380 = t202 * t386 + t203 * t388 + t230 * t385;
t379 = t204 * t386 + t205 * t388 + t232 * t385;
t378 = t202 * t388 + t203 * t389 + t230 * t387;
t377 = t204 * t388 + t205 * t389 + t232 * t387;
t270 = t332 * t311;
t258 = t270 * t364 + t312 * t315;
t269 = t280 * t311;
t224 = -t258 * t314 - t269 * t318;
t353 = t269 * t314;
t225 = t258 * t318 - t353;
t257 = t270 * t315 - t312 * t364;
t376 = t224 * t385 + t225 * t387 + t257 * t384;
t375 = t224 * t386 + t225 * t388 + t257 * t385;
t374 = t224 * t388 + t225 * t389 + t257 * t387;
t208 = Icges(4,5) * t253 + Icges(4,6) * t252 - Icges(4,3) * t351;
t347 = t319 * t320;
t349 = t317 * t316;
t274 = t312 * t347 - t349;
t348 = t317 * t319;
t350 = t316 * t320;
t275 = t312 * t350 + t348;
t239 = Icges(3,5) * t275 + Icges(3,6) * t274 - Icges(3,3) * t351;
t373 = t208 + t239;
t209 = Icges(4,5) * t255 + Icges(4,6) * t254 + Icges(4,3) * t352;
t276 = -t312 * t348 - t350;
t277 = -t312 * t349 + t347;
t240 = Icges(3,5) * t277 + Icges(3,6) * t276 + Icges(3,3) * t352;
t372 = t209 + t240;
t235 = Icges(4,5) * t270 + Icges(4,6) * t269 + Icges(4,3) * t312;
t265 = Icges(3,3) * t312 + (Icges(3,5) * t316 + Icges(3,6) * t319) * t311;
t371 = t235 + t265;
t363 = pkin(2) * t316;
t362 = pkin(8) * t312;
t361 = pkin(2) * t319;
t360 = pkin(5) * t318;
t356 = Icges(2,4) * t317;
t346 = rSges(7,1) * t203 + rSges(7,2) * t202 - pkin(5) * t355 + t230 * t383 + t231 * t360;
t345 = rSges(7,1) * t205 + rSges(7,2) * t204 - pkin(5) * t354 + t232 * t383 + t233 * t360;
t344 = rSges(7,1) * t225 + rSges(7,2) * t224 - pkin(5) * t353 + t257 * t383 + t258 * t360;
t343 = qJD(2) * t311;
t342 = qJD(3) * t311;
t341 = V_base(5) * pkin(7) + V_base(1);
t288 = t317 * t343 + V_base(4);
t308 = V_base(6) + qJD(1);
t337 = -qJ(3) * t311 + t312 * t363;
t229 = -qJD(4) * t254 + t288;
t289 = qJD(2) * t312 + t308;
t256 = -qJD(4) * t269 + t289;
t287 = -t320 * t343 + V_base(5);
t282 = pkin(1) * t317 - pkin(8) * t351;
t334 = -t282 * t308 + t362 * V_base(5) + t341;
t283 = pkin(1) * t320 + pkin(8) * t352;
t333 = t282 * V_base(4) - t283 * V_base(5) + V_base(3);
t228 = -qJD(4) * t252 + t287;
t281 = qJ(3) * t312 + t311 * t363;
t331 = t281 * t287 + t317 * t342 + t334;
t250 = t317 * t361 + t320 * t337;
t330 = qJD(3) * t312 + t250 * t288 + t333;
t329 = t308 * t283 + V_base(2) + (-pkin(7) - t362) * V_base(4);
t220 = pkin(3) * t253 - pkin(9) * t252;
t247 = pkin(3) * t270 - pkin(9) * t269;
t327 = t287 * t247 + (-t220 - t250) * t289 + t331;
t221 = pkin(3) * t255 - pkin(9) * t254;
t251 = -t317 * t337 + t320 * t361;
t326 = t288 * t220 + (-t221 - t251) * t287 + t330;
t325 = t251 * t289 - t320 * t342 + t329;
t198 = pkin(4) * t231 + pkin(10) * t230;
t222 = pkin(4) * t258 + pkin(10) * t257;
t324 = -t198 * t256 + t222 * t228 + t327;
t199 = pkin(4) * t233 + pkin(10) * t232;
t323 = t198 * t229 - t199 * t228 + t326;
t322 = t289 * t221 + (-t247 - t281) * t288 + t325;
t321 = t199 * t256 - t222 * t229 + t322;
t309 = Icges(2,4) * t320;
t297 = rSges(2,1) * t320 - rSges(2,2) * t317;
t296 = rSges(2,1) * t317 + rSges(2,2) * t320;
t295 = Icges(2,1) * t320 - t356;
t294 = Icges(2,1) * t317 + t309;
t293 = -Icges(2,2) * t317 + t309;
t292 = Icges(2,2) * t320 + t356;
t286 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t285 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t284 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t268 = rSges(3,3) * t312 + (rSges(3,1) * t316 + rSges(3,2) * t319) * t311;
t267 = Icges(3,5) * t312 + (Icges(3,1) * t316 + Icges(3,4) * t319) * t311;
t266 = Icges(3,6) * t312 + (Icges(3,4) * t316 + Icges(3,2) * t319) * t311;
t261 = V_base(5) * rSges(2,3) - t296 * t308 + t341;
t260 = t297 * t308 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t259 = t296 * V_base(4) - t297 * V_base(5) + V_base(3);
t246 = rSges(3,1) * t277 + rSges(3,2) * t276 + rSges(3,3) * t352;
t245 = rSges(3,1) * t275 + rSges(3,2) * t274 - rSges(3,3) * t351;
t244 = Icges(3,1) * t277 + Icges(3,4) * t276 + Icges(3,5) * t352;
t243 = Icges(3,1) * t275 + Icges(3,4) * t274 - Icges(3,5) * t351;
t242 = Icges(3,4) * t277 + Icges(3,2) * t276 + Icges(3,6) * t352;
t241 = Icges(3,4) * t275 + Icges(3,2) * t274 - Icges(3,6) * t351;
t238 = rSges(4,1) * t270 + rSges(4,2) * t269 + rSges(4,3) * t312;
t237 = Icges(4,1) * t270 + Icges(4,4) * t269 + Icges(4,5) * t312;
t236 = Icges(4,4) * t270 + Icges(4,2) * t269 + Icges(4,6) * t312;
t223 = qJD(5) * t257 + t256;
t219 = rSges(5,1) * t258 - rSges(5,2) * t257 - rSges(5,3) * t269;
t218 = Icges(5,1) * t258 - Icges(5,4) * t257 - Icges(5,5) * t269;
t217 = Icges(5,4) * t258 - Icges(5,2) * t257 - Icges(5,6) * t269;
t216 = Icges(5,5) * t258 - Icges(5,6) * t257 - Icges(5,3) * t269;
t215 = rSges(4,1) * t255 + rSges(4,2) * t254 + rSges(4,3) * t352;
t214 = rSges(4,1) * t253 + rSges(4,2) * t252 - rSges(4,3) * t351;
t213 = Icges(4,1) * t255 + Icges(4,4) * t254 + Icges(4,5) * t352;
t212 = Icges(4,1) * t253 + Icges(4,4) * t252 - Icges(4,5) * t351;
t211 = Icges(4,4) * t255 + Icges(4,2) * t254 + Icges(4,6) * t352;
t210 = Icges(4,4) * t253 + Icges(4,2) * t252 - Icges(4,6) * t351;
t201 = qJD(5) * t232 + t229;
t200 = qJD(5) * t230 + t228;
t197 = -t245 * t289 + t268 * t287 + t334;
t196 = t246 * t289 - t268 * t288 + t329;
t193 = t245 * t288 - t246 * t287 + t333;
t192 = rSges(6,1) * t225 + rSges(6,2) * t224 + rSges(6,3) * t257;
t184 = rSges(5,1) * t233 - rSges(5,2) * t232 - rSges(5,3) * t254;
t183 = rSges(5,1) * t231 - rSges(5,2) * t230 - rSges(5,3) * t252;
t182 = Icges(5,1) * t233 - Icges(5,4) * t232 - Icges(5,5) * t254;
t181 = Icges(5,1) * t231 - Icges(5,4) * t230 - Icges(5,5) * t252;
t180 = Icges(5,4) * t233 - Icges(5,2) * t232 - Icges(5,6) * t254;
t179 = Icges(5,4) * t231 - Icges(5,2) * t230 - Icges(5,6) * t252;
t178 = Icges(5,5) * t233 - Icges(5,6) * t232 - Icges(5,3) * t254;
t177 = Icges(5,5) * t231 - Icges(5,6) * t230 - Icges(5,3) * t252;
t174 = rSges(6,1) * t205 + rSges(6,2) * t204 + rSges(6,3) * t232;
t172 = rSges(6,1) * t203 + rSges(6,2) * t202 + rSges(6,3) * t230;
t156 = t238 * t287 + (-t214 - t250) * t289 + t331;
t155 = t289 * t215 + (-t238 - t281) * t288 + t325;
t154 = t214 * t288 + (-t215 - t251) * t287 + t330;
t153 = -t183 * t256 + t219 * t228 + t327;
t152 = t184 * t256 - t219 * t229 + t322;
t151 = t183 * t229 - t184 * t228 + t326;
t150 = -t172 * t223 + t192 * t200 + t324;
t149 = t174 * t223 - t192 * t201 + t321;
t148 = t172 * t201 - t174 * t200 + t323;
t147 = qJD(6) * t232 + t200 * t344 - t223 * t346 + t324;
t146 = qJD(6) * t230 - t201 * t344 + t223 * t345 + t321;
t145 = qJD(6) * t257 - t200 * t345 + t201 * t346 + t323;
t1 = m(6) * (t148 ^ 2 + t149 ^ 2 + t150 ^ 2) / 0.2e1 + m(7) * (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(1) * (t284 ^ 2 + t285 ^ 2 + t286 ^ 2) / 0.2e1 + t256 * ((-t178 * t269 - t180 * t257 + t182 * t258) * t229 + (-t177 * t269 - t179 * t257 + t181 * t258) * t228 + (-t216 * t269 - t217 * t257 + t218 * t258) * t256) / 0.2e1 + m(2) * (t259 ^ 2 + t260 ^ 2 + t261 ^ 2) / 0.2e1 + t229 * ((-t254 * t178 - t232 * t180 + t233 * t182) * t229 + (-t177 * t254 - t179 * t232 + t181 * t233) * t228 + (-t216 * t254 - t217 * t232 + t218 * t233) * t256) / 0.2e1 + t228 * ((-t178 * t252 - t180 * t230 + t182 * t231) * t229 + (-t252 * t177 - t230 * t179 + t231 * t181) * t228 + (-t216 * t252 - t217 * t230 + t218 * t231) * t256) / 0.2e1 + m(3) * (t193 ^ 2 + t196 ^ 2 + t197 ^ 2) / 0.2e1 + m(5) * (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + m(4) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + ((t202 * t375 + t203 * t374 + t230 * t376) * t223 + (t202 * t379 + t203 * t377 + t230 * t381) * t201 + (t380 * t202 + t378 * t203 + t230 * t382) * t200) * t200 / 0.2e1 + ((t204 * t375 + t205 * t374 + t232 * t376) * t223 + (t204 * t379 + t205 * t377 + t232 * t381) * t201 + (t204 * t380 + t205 * t378 + t232 * t382) * t200) * t201 / 0.2e1 + ((t224 * t375 + t225 * t374 + t257 * t376) * t223 + (t224 * t379 + t225 * t377 + t257 * t381) * t201 + (t224 * t380 + t225 * t378 + t257 * t382) * t200) * t223 / 0.2e1 + ((t252 * t236 + t253 * t237 + t274 * t266 + t275 * t267 - t351 * t371) * t289 + (t211 * t252 + t213 * t253 + t242 * t274 + t244 * t275 - t351 * t372) * t288 + (t252 * t210 + t253 * t212 + t274 * t241 + t275 * t243 - t351 * t373) * t287) * t287 / 0.2e1 + ((t236 * t254 + t237 * t255 + t266 * t276 + t267 * t277 + t352 * t371) * t289 + (t211 * t254 + t213 * t255 + t242 * t276 + t244 * t277 + t352 * t372) * t288 + (t210 * t254 + t212 * t255 + t241 * t276 + t243 * t277 + t352 * t373) * t287) * t288 / 0.2e1 + ((t239 * t287 + t240 * t288 + t265 * t289) * t312 + ((t242 * t319 + t244 * t316) * t288 + (t241 * t319 + t243 * t316) * t287 + (t266 * t319 + t267 * t316) * t289) * t311 + (t209 * t312 + t211 * t269 + t213 * t270) * t288 + (t208 * t312 + t210 * t269 + t212 * t270) * t287 + (t235 * t312 + t236 * t269 + t237 * t270) * t289) * t289 / 0.2e1 + ((-t292 * t317 + t294 * t320 + Icges(1,4)) * V_base(5) + (-t317 * t293 + t295 * t320 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t292 * t320 + t317 * t294 + Icges(1,2)) * V_base(5) + (t293 * t320 + t295 * t317 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t317 + Icges(2,6) * t320) * V_base(5) + (Icges(2,5) * t320 - Icges(2,6) * t317) * V_base(4) + Icges(2,3) * t308 / 0.2e1) * t308;
T  = t1;
