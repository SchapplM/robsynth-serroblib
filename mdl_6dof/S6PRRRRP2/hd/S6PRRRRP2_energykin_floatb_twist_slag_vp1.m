% Calculate kinetic energy for
% S6PRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRRRP2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRRRP2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRRRP2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:01:04
% EndTime: 2019-03-09 00:01:08
% DurationCPUTime: 4.51s
% Computational Cost: add. (3504->387), mult. (6139->563), div. (0->0), fcn. (7307->12), ass. (0->176)
t386 = Icges(6,1) + Icges(7,1);
t385 = -Icges(6,4) + Icges(7,5);
t384 = Icges(7,4) + Icges(6,5);
t383 = Icges(6,2) + Icges(7,3);
t382 = Icges(7,2) + Icges(6,3);
t381 = -Icges(6,6) + Icges(7,6);
t380 = rSges(7,1) + pkin(5);
t379 = rSges(7,3) + qJ(6);
t315 = sin(pkin(11));
t317 = cos(pkin(11));
t323 = cos(qJ(2));
t318 = cos(pkin(6));
t321 = sin(qJ(2));
t352 = t318 * t321;
t280 = t315 * t323 + t317 * t352;
t350 = qJ(3) + qJ(4);
t313 = sin(t350);
t338 = cos(t350);
t316 = sin(pkin(6));
t358 = t316 * t317;
t251 = t280 * t338 - t313 * t358;
t351 = t318 * t323;
t279 = t315 * t321 - t317 * t351;
t319 = sin(qJ(5));
t364 = cos(qJ(5));
t221 = t251 * t319 - t279 * t364;
t222 = t251 * t364 + t279 * t319;
t337 = t316 * t338;
t250 = t280 * t313 + t317 * t337;
t377 = t221 * t383 + t222 * t385 + t250 * t381;
t282 = -t315 * t352 + t317 * t323;
t359 = t315 * t316;
t253 = t282 * t338 + t313 * t359;
t281 = t315 * t351 + t317 * t321;
t223 = t253 * t319 - t281 * t364;
t224 = t253 * t364 + t281 * t319;
t252 = t282 * t313 - t315 * t337;
t376 = t223 * t383 + t224 * t385 + t252 * t381;
t375 = t221 * t381 + t222 * t384 + t250 * t382;
t374 = t223 * t381 + t224 * t384 + t252 * t382;
t373 = t385 * t221 + t222 * t386 + t384 * t250;
t372 = t385 * t223 + t224 * t386 + t384 * t252;
t272 = t318 * t313 + t321 * t337;
t354 = t316 * t323;
t254 = t272 * t319 + t354 * t364;
t255 = t272 * t364 - t319 * t354;
t356 = t316 * t321;
t271 = t313 * t356 - t318 * t338;
t371 = t254 * t383 + t255 * t385 + t271 * t381;
t370 = t254 * t381 + t255 * t384 + t271 * t382;
t369 = t385 * t254 + t255 * t386 + t384 * t271;
t363 = pkin(7) * t318;
t322 = cos(qJ(3));
t362 = pkin(3) * t322;
t360 = Icges(2,4) * t315;
t320 = sin(qJ(3));
t357 = t316 * t320;
t355 = t316 * t322;
t353 = t318 * t320;
t349 = rSges(7,2) * t250 + t379 * t221 + t380 * t222;
t348 = rSges(7,2) * t252 + t379 * t223 + t380 * t224;
t347 = rSges(7,2) * t271 + t379 * t254 + t380 * t255;
t346 = qJD(2) * t316;
t345 = V_base(5) * qJ(1) + V_base(1);
t341 = qJD(1) + V_base(3);
t340 = t315 * t357;
t339 = t317 * t357;
t295 = t315 * t346 + V_base(4);
t306 = qJD(2) * t318 + V_base(6);
t258 = qJD(3) * t281 + t295;
t231 = qJD(4) * t281 + t258;
t294 = -t317 * t346 + V_base(5);
t257 = qJD(3) * t279 + t294;
t289 = pkin(1) * t315 - pkin(7) * t358;
t336 = -t289 * V_base(6) + V_base(5) * t363 + t345;
t290 = pkin(1) * t317 + pkin(7) * t359;
t335 = V_base(4) * t289 - t290 * V_base(5) + t341;
t230 = qJD(4) * t279 + t257;
t266 = (-qJD(3) - qJD(4)) * t354 + t306;
t334 = V_base(6) * t290 + V_base(2) + (-qJ(1) - t363) * V_base(4);
t248 = t280 * pkin(2) + t279 * pkin(8);
t288 = (pkin(2) * t321 - pkin(8) * t323) * t316;
t333 = -t248 * t306 + t294 * t288 + t336;
t249 = t282 * pkin(2) + t281 * pkin(8);
t332 = t295 * t248 - t249 * t294 + t335;
t331 = t306 * t249 - t288 * t295 + t334;
t205 = -pkin(3) * t339 + pkin(9) * t279 + t280 * t362;
t247 = pkin(3) * t353 + (-pkin(9) * t323 + t321 * t362) * t316;
t283 = -qJD(3) * t354 + t306;
t330 = -t205 * t283 + t257 * t247 + t333;
t206 = pkin(3) * t340 + pkin(9) * t281 + t282 * t362;
t329 = t258 * t205 - t206 * t257 + t332;
t328 = t283 * t206 - t247 * t258 + t331;
t218 = pkin(4) * t251 + pkin(10) * t250;
t246 = pkin(4) * t272 + pkin(10) * t271;
t327 = -t218 * t266 + t230 * t246 + t330;
t219 = pkin(4) * t253 + pkin(10) * t252;
t326 = t231 * t218 - t219 * t230 + t329;
t325 = t266 * t219 - t231 * t246 + t328;
t312 = Icges(2,4) * t317;
t304 = rSges(2,1) * t317 - rSges(2,2) * t315;
t303 = rSges(2,1) * t315 + rSges(2,2) * t317;
t302 = Icges(2,1) * t317 - t360;
t301 = Icges(2,1) * t315 + t312;
t300 = -Icges(2,2) * t315 + t312;
t299 = Icges(2,2) * t317 + t360;
t293 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t292 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t291 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t287 = t321 * t355 + t353;
t286 = t318 * t322 - t320 * t356;
t270 = rSges(3,3) * t318 + (rSges(3,1) * t321 + rSges(3,2) * t323) * t316;
t269 = Icges(3,5) * t318 + (Icges(3,1) * t321 + Icges(3,4) * t323) * t316;
t268 = Icges(3,6) * t318 + (Icges(3,4) * t321 + Icges(3,2) * t323) * t316;
t267 = Icges(3,3) * t318 + (Icges(3,5) * t321 + Icges(3,6) * t323) * t316;
t265 = V_base(5) * rSges(2,3) - t303 * V_base(6) + t345;
t264 = t304 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t262 = t282 * t322 + t340;
t261 = -t282 * t320 + t315 * t355;
t260 = t280 * t322 - t339;
t259 = -t280 * t320 - t317 * t355;
t256 = t303 * V_base(4) - t304 * V_base(5) + t341;
t245 = rSges(4,1) * t287 + rSges(4,2) * t286 - rSges(4,3) * t354;
t244 = Icges(4,1) * t287 + Icges(4,4) * t286 - Icges(4,5) * t354;
t243 = Icges(4,4) * t287 + Icges(4,2) * t286 - Icges(4,6) * t354;
t242 = Icges(4,5) * t287 + Icges(4,6) * t286 - Icges(4,3) * t354;
t241 = rSges(3,1) * t282 - rSges(3,2) * t281 + rSges(3,3) * t359;
t240 = rSges(3,1) * t280 - rSges(3,2) * t279 - rSges(3,3) * t358;
t239 = Icges(3,1) * t282 - Icges(3,4) * t281 + Icges(3,5) * t359;
t238 = Icges(3,1) * t280 - Icges(3,4) * t279 - Icges(3,5) * t358;
t237 = Icges(3,4) * t282 - Icges(3,2) * t281 + Icges(3,6) * t359;
t236 = Icges(3,4) * t280 - Icges(3,2) * t279 - Icges(3,6) * t358;
t235 = Icges(3,5) * t282 - Icges(3,6) * t281 + Icges(3,3) * t359;
t234 = Icges(3,5) * t280 - Icges(3,6) * t279 - Icges(3,3) * t358;
t232 = qJD(5) * t271 + t266;
t228 = rSges(5,1) * t272 - rSges(5,2) * t271 - rSges(5,3) * t354;
t227 = Icges(5,1) * t272 - Icges(5,4) * t271 - Icges(5,5) * t354;
t226 = Icges(5,4) * t272 - Icges(5,2) * t271 - Icges(5,6) * t354;
t225 = Icges(5,5) * t272 - Icges(5,6) * t271 - Icges(5,3) * t354;
t216 = rSges(4,1) * t262 + rSges(4,2) * t261 + rSges(4,3) * t281;
t215 = rSges(4,1) * t260 + rSges(4,2) * t259 + rSges(4,3) * t279;
t214 = Icges(4,1) * t262 + Icges(4,4) * t261 + Icges(4,5) * t281;
t213 = Icges(4,1) * t260 + Icges(4,4) * t259 + Icges(4,5) * t279;
t212 = Icges(4,4) * t262 + Icges(4,2) * t261 + Icges(4,6) * t281;
t211 = Icges(4,4) * t260 + Icges(4,2) * t259 + Icges(4,6) * t279;
t210 = Icges(4,5) * t262 + Icges(4,6) * t261 + Icges(4,3) * t281;
t209 = Icges(4,5) * t260 + Icges(4,6) * t259 + Icges(4,3) * t279;
t208 = qJD(5) * t252 + t231;
t207 = qJD(5) * t250 + t230;
t204 = rSges(5,1) * t253 - rSges(5,2) * t252 + rSges(5,3) * t281;
t203 = rSges(5,1) * t251 - rSges(5,2) * t250 + rSges(5,3) * t279;
t202 = Icges(5,1) * t253 - Icges(5,4) * t252 + Icges(5,5) * t281;
t201 = Icges(5,1) * t251 - Icges(5,4) * t250 + Icges(5,5) * t279;
t200 = Icges(5,4) * t253 - Icges(5,2) * t252 + Icges(5,6) * t281;
t199 = Icges(5,4) * t251 - Icges(5,2) * t250 + Icges(5,6) * t279;
t198 = Icges(5,5) * t253 - Icges(5,6) * t252 + Icges(5,3) * t281;
t197 = Icges(5,5) * t251 - Icges(5,6) * t250 + Icges(5,3) * t279;
t196 = rSges(6,1) * t255 - rSges(6,2) * t254 + rSges(6,3) * t271;
t185 = -t240 * t306 + t270 * t294 + t336;
t184 = t241 * t306 - t270 * t295 + t334;
t179 = t240 * t295 - t241 * t294 + t335;
t178 = rSges(6,1) * t224 - rSges(6,2) * t223 + rSges(6,3) * t252;
t176 = rSges(6,1) * t222 - rSges(6,2) * t221 + rSges(6,3) * t250;
t162 = -t215 * t283 + t245 * t257 + t333;
t161 = t216 * t283 - t245 * t258 + t331;
t160 = t215 * t258 - t216 * t257 + t332;
t159 = -t203 * t266 + t228 * t230 + t330;
t158 = t204 * t266 - t228 * t231 + t328;
t157 = t203 * t231 - t204 * t230 + t329;
t156 = -t176 * t232 + t196 * t207 + t327;
t155 = t178 * t232 - t196 * t208 + t325;
t154 = t176 * t208 - t178 * t207 + t326;
t153 = qJD(6) * t223 + t207 * t347 - t232 * t349 + t327;
t152 = qJD(6) * t221 - t208 * t347 + t232 * t348 + t325;
t151 = qJD(6) * t254 - t207 * t348 + t208 * t349 + t326;
t1 = t283 * ((-t210 * t354 + t212 * t286 + t214 * t287) * t258 + (-t209 * t354 + t211 * t286 + t213 * t287) * t257 + (-t242 * t354 + t243 * t286 + t244 * t287) * t283) / 0.2e1 + t266 * ((-t198 * t354 - t200 * t271 + t202 * t272) * t231 + (-t197 * t354 - t199 * t271 + t201 * t272) * t230 + (-t225 * t354 - t226 * t271 + t227 * t272) * t266) / 0.2e1 + m(3) * (t179 ^ 2 + t184 ^ 2 + t185 ^ 2) / 0.2e1 + m(2) * (t256 ^ 2 + t264 ^ 2 + t265 ^ 2) / 0.2e1 + t230 * ((t198 * t279 - t200 * t250 + t202 * t251) * t231 + (t197 * t279 - t199 * t250 + t201 * t251) * t230 + (t225 * t279 - t226 * t250 + t227 * t251) * t266) / 0.2e1 + t231 * ((t198 * t281 - t200 * t252 + t202 * t253) * t231 + (t197 * t281 - t199 * t252 + t201 * t253) * t230 + (t225 * t281 - t226 * t252 + t227 * t253) * t266) / 0.2e1 + t257 * ((t210 * t279 + t212 * t259 + t214 * t260) * t258 + (t209 * t279 + t211 * t259 + t213 * t260) * t257 + (t242 * t279 + t243 * t259 + t244 * t260) * t283) / 0.2e1 + t258 * ((t210 * t281 + t212 * t261 + t214 * t262) * t258 + (t209 * t281 + t211 * t261 + t213 * t262) * t257 + (t242 * t281 + t243 * t261 + t244 * t262) * t283) / 0.2e1 + m(1) * (t291 ^ 2 + t292 ^ 2 + t293 ^ 2) / 0.2e1 + t306 * ((t234 * t294 + t235 * t295 + t267 * t306) * t318 + ((t237 * t323 + t239 * t321) * t295 + (t236 * t323 + t238 * t321) * t294 + (t268 * t323 + t269 * t321) * t306) * t316) / 0.2e1 + t294 * ((-t235 * t358 - t237 * t279 + t239 * t280) * t295 + (-t234 * t358 - t236 * t279 + t238 * t280) * t294 + (-t267 * t358 - t268 * t279 + t269 * t280) * t306) / 0.2e1 + t295 * ((t235 * t359 - t237 * t281 + t239 * t282) * t295 + (t234 * t359 - t236 * t281 + t238 * t282) * t294 + (t267 * t359 - t268 * t281 + t269 * t282) * t306) / 0.2e1 + m(7) * (t151 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + m(4) * (t160 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + m(6) * (t154 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + m(5) * (t157 ^ 2 + t158 ^ 2 + t159 ^ 2) / 0.2e1 + ((t221 * t371 + t222 * t369 + t250 * t370) * t232 + (t221 * t376 + t222 * t372 + t250 * t374) * t208 + (t377 * t221 + t373 * t222 + t375 * t250) * t207) * t207 / 0.2e1 + ((t223 * t371 + t224 * t369 + t252 * t370) * t232 + (t376 * t223 + t372 * t224 + t374 * t252) * t208 + (t223 * t377 + t224 * t373 + t252 * t375) * t207) * t208 / 0.2e1 + ((t371 * t254 + t369 * t255 + t370 * t271) * t232 + (t254 * t376 + t255 * t372 + t271 * t374) * t208 + (t254 * t377 + t255 * t373 + t271 * t375) * t207) * t232 / 0.2e1 + ((-t299 * t315 + t301 * t317 + Icges(1,4)) * V_base(5) + (-t300 * t315 + t302 * t317 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t299 * t317 + t301 * t315 + Icges(1,2)) * V_base(5) + (t300 * t317 + t302 * t315 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t317 - Icges(2,6) * t315 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t315 + Icges(2,6) * t317 + Icges(1,6)) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
