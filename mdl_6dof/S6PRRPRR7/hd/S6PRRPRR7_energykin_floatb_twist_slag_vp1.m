% Calculate kinetic energy for
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRR7_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPRR7_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR7_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR7_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPRR7_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:30:09
% EndTime: 2019-03-08 22:30:13
% DurationCPUTime: 4.32s
% Computational Cost: add. (2715->396), mult. (5944->569), div. (0->0), fcn. (7097->12), ass. (0->175)
t362 = Icges(4,1) + Icges(5,2);
t361 = Icges(5,1) + Icges(4,3);
t360 = -Icges(4,4) - Icges(5,6);
t359 = Icges(5,4) - Icges(4,5);
t358 = Icges(5,5) - Icges(4,6);
t357 = Icges(4,2) + Icges(5,3);
t300 = sin(pkin(11));
t302 = cos(pkin(11));
t307 = cos(qJ(2));
t303 = cos(pkin(6));
t305 = sin(qJ(2));
t330 = t303 * t305;
t262 = t300 * t307 + t302 * t330;
t301 = sin(pkin(6));
t342 = cos(qJ(3));
t322 = t301 * t342;
t341 = sin(qJ(3));
t243 = t262 * t341 + t302 * t322;
t321 = t301 * t341;
t244 = t262 * t342 - t302 * t321;
t329 = t303 * t307;
t261 = t300 * t305 - t302 * t329;
t355 = t243 * t357 + t244 * t360 + t261 * t358;
t264 = -t300 * t330 + t302 * t307;
t245 = t264 * t341 - t300 * t322;
t246 = t264 * t342 + t300 * t321;
t263 = t300 * t329 + t302 * t305;
t354 = t245 * t357 + t246 * t360 + t263 * t358;
t353 = t243 * t358 - t244 * t359 + t261 * t361;
t352 = t245 * t358 - t246 * t359 + t263 * t361;
t351 = t360 * t243 + t244 * t362 - t359 * t261;
t350 = t360 * t245 + t246 * t362 - t359 * t263;
t268 = -t303 * t342 + t305 * t321;
t269 = t303 * t341 + t305 * t322;
t331 = t301 * t307;
t349 = t268 * t357 + t269 * t360 - t331 * t358;
t348 = t360 * t268 + t269 * t362 + t359 * t331;
t347 = t268 * t358 - t269 * t359 - t331 * t361;
t340 = pkin(7) * t303;
t306 = cos(qJ(5));
t339 = pkin(5) * t306;
t337 = Icges(2,4) * t300;
t304 = sin(qJ(5));
t336 = t243 * t304;
t335 = t245 * t304;
t334 = t268 * t304;
t333 = t300 * t301;
t332 = t301 * t302;
t328 = qJD(2) * t301;
t327 = V_base(5) * qJ(1) + V_base(1);
t323 = qJD(1) + V_base(3);
t277 = t300 * t328 + V_base(4);
t289 = qJD(2) * t303 + V_base(6);
t242 = qJD(3) * t263 + t277;
t199 = qJD(5) * t246 + t242;
t276 = -t302 * t328 + V_base(5);
t241 = qJD(3) * t261 + t276;
t265 = -qJD(3) * t331 + t289;
t271 = pkin(1) * t300 - pkin(7) * t332;
t320 = -t271 * V_base(6) + V_base(5) * t340 + t327;
t272 = pkin(1) * t302 + pkin(7) * t333;
t319 = V_base(4) * t271 - t272 * V_base(5) + t323;
t198 = qJD(5) * t244 + t241;
t233 = qJD(5) * t269 + t265;
t318 = V_base(6) * t272 + V_base(2) + (-qJ(1) - t340) * V_base(4);
t230 = pkin(2) * t262 + pkin(8) * t261;
t270 = (pkin(2) * t305 - pkin(8) * t307) * t301;
t317 = -t230 * t289 + t276 * t270 + t320;
t231 = pkin(2) * t264 + pkin(8) * t263;
t316 = t277 * t230 - t231 * t276 + t319;
t315 = t289 * t231 - t270 * t277 + t318;
t232 = pkin(3) * t269 + qJ(4) * t268;
t314 = qJD(4) * t245 + t241 * t232 + t317;
t196 = pkin(3) * t244 + qJ(4) * t243;
t313 = qJD(4) * t268 + t242 * t196 + t316;
t197 = pkin(3) * t246 + qJ(4) * t245;
t312 = qJD(4) * t243 + t265 * t197 + t315;
t210 = t261 * pkin(4) + t244 * pkin(9);
t250 = -pkin(4) * t331 + t269 * pkin(9);
t311 = t241 * t250 + (-t196 - t210) * t265 + t314;
t211 = t263 * pkin(4) + t246 * pkin(9);
t310 = t242 * t210 + (-t197 - t211) * t241 + t313;
t309 = t265 * t211 + (-t232 - t250) * t242 + t312;
t299 = qJ(5) + qJ(6);
t297 = cos(t299);
t296 = sin(t299);
t295 = Icges(2,4) * t302;
t285 = rSges(2,1) * t302 - rSges(2,2) * t300;
t284 = rSges(2,1) * t300 + rSges(2,2) * t302;
t283 = Icges(2,1) * t302 - t337;
t282 = Icges(2,1) * t300 + t295;
t281 = -Icges(2,2) * t300 + t295;
t280 = Icges(2,2) * t302 + t337;
t275 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t274 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t273 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t256 = rSges(3,3) * t303 + (rSges(3,1) * t305 + rSges(3,2) * t307) * t301;
t255 = Icges(3,5) * t303 + (Icges(3,1) * t305 + Icges(3,4) * t307) * t301;
t254 = Icges(3,6) * t303 + (Icges(3,4) * t305 + Icges(3,2) * t307) * t301;
t253 = Icges(3,3) * t303 + (Icges(3,5) * t305 + Icges(3,6) * t307) * t301;
t252 = V_base(5) * rSges(2,3) - t284 * V_base(6) + t327;
t251 = t285 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t248 = -t306 * t331 + t334;
t247 = t268 * t306 + t304 * t331;
t240 = t284 * V_base(4) - t285 * V_base(5) + t323;
t235 = t268 * t296 - t297 * t331;
t234 = t268 * t297 + t296 * t331;
t229 = rSges(4,1) * t269 - rSges(4,2) * t268 - rSges(4,3) * t331;
t228 = -rSges(5,1) * t331 - rSges(5,2) * t269 + rSges(5,3) * t268;
t221 = rSges(3,1) * t264 - rSges(3,2) * t263 + rSges(3,3) * t333;
t220 = rSges(3,1) * t262 - rSges(3,2) * t261 - rSges(3,3) * t332;
t219 = Icges(3,1) * t264 - Icges(3,4) * t263 + Icges(3,5) * t333;
t218 = Icges(3,1) * t262 - Icges(3,4) * t261 - Icges(3,5) * t332;
t217 = Icges(3,4) * t264 - Icges(3,2) * t263 + Icges(3,6) * t333;
t216 = Icges(3,4) * t262 - Icges(3,2) * t261 - Icges(3,6) * t332;
t215 = Icges(3,5) * t264 - Icges(3,6) * t263 + Icges(3,3) * t333;
t214 = Icges(3,5) * t262 - Icges(3,6) * t261 - Icges(3,3) * t332;
t209 = qJD(6) * t269 + t233;
t208 = t263 * t306 + t335;
t207 = t245 * t306 - t263 * t304;
t206 = t261 * t306 + t336;
t205 = t243 * t306 - t261 * t304;
t203 = t245 * t296 + t263 * t297;
t202 = t245 * t297 - t263 * t296;
t201 = t243 * t296 + t261 * t297;
t200 = t243 * t297 - t261 * t296;
t193 = pkin(5) * t334 + pkin(10) * t269 - t331 * t339;
t192 = rSges(6,1) * t248 + rSges(6,2) * t247 + rSges(6,3) * t269;
t189 = Icges(6,1) * t248 + Icges(6,4) * t247 + Icges(6,5) * t269;
t188 = Icges(6,4) * t248 + Icges(6,2) * t247 + Icges(6,6) * t269;
t187 = Icges(6,5) * t248 + Icges(6,6) * t247 + Icges(6,3) * t269;
t186 = rSges(4,1) * t246 - rSges(4,2) * t245 + rSges(4,3) * t263;
t185 = rSges(4,1) * t244 - rSges(4,2) * t243 + rSges(4,3) * t261;
t184 = rSges(5,1) * t263 - rSges(5,2) * t246 + rSges(5,3) * t245;
t183 = rSges(5,1) * t261 - rSges(5,2) * t244 + rSges(5,3) * t243;
t170 = rSges(7,1) * t235 + rSges(7,2) * t234 + rSges(7,3) * t269;
t169 = Icges(7,1) * t235 + Icges(7,4) * t234 + Icges(7,5) * t269;
t168 = Icges(7,4) * t235 + Icges(7,2) * t234 + Icges(7,6) * t269;
t167 = Icges(7,5) * t235 + Icges(7,6) * t234 + Icges(7,3) * t269;
t166 = qJD(6) * t246 + t199;
t165 = qJD(6) * t244 + t198;
t163 = -t220 * t289 + t256 * t276 + t320;
t162 = t221 * t289 - t256 * t277 + t318;
t161 = pkin(5) * t335 + pkin(10) * t246 + t263 * t339;
t160 = pkin(5) * t336 + pkin(10) * t244 + t261 * t339;
t159 = rSges(6,1) * t208 + rSges(6,2) * t207 + rSges(6,3) * t246;
t158 = rSges(6,1) * t206 + rSges(6,2) * t205 + rSges(6,3) * t244;
t157 = Icges(6,1) * t208 + Icges(6,4) * t207 + Icges(6,5) * t246;
t156 = Icges(6,1) * t206 + Icges(6,4) * t205 + Icges(6,5) * t244;
t155 = Icges(6,4) * t208 + Icges(6,2) * t207 + Icges(6,6) * t246;
t154 = Icges(6,4) * t206 + Icges(6,2) * t205 + Icges(6,6) * t244;
t153 = Icges(6,5) * t208 + Icges(6,6) * t207 + Icges(6,3) * t246;
t152 = Icges(6,5) * t206 + Icges(6,6) * t205 + Icges(6,3) * t244;
t151 = t220 * t277 - t221 * t276 + t319;
t150 = rSges(7,1) * t203 + rSges(7,2) * t202 + rSges(7,3) * t246;
t149 = rSges(7,1) * t201 + rSges(7,2) * t200 + rSges(7,3) * t244;
t148 = Icges(7,1) * t203 + Icges(7,4) * t202 + Icges(7,5) * t246;
t147 = Icges(7,1) * t201 + Icges(7,4) * t200 + Icges(7,5) * t244;
t146 = Icges(7,4) * t203 + Icges(7,2) * t202 + Icges(7,6) * t246;
t145 = Icges(7,4) * t201 + Icges(7,2) * t200 + Icges(7,6) * t244;
t144 = Icges(7,5) * t203 + Icges(7,6) * t202 + Icges(7,3) * t246;
t143 = Icges(7,5) * t201 + Icges(7,6) * t200 + Icges(7,3) * t244;
t142 = -t185 * t265 + t229 * t241 + t317;
t141 = t186 * t265 - t229 * t242 + t315;
t140 = t185 * t242 - t186 * t241 + t316;
t139 = t228 * t241 + (-t183 - t196) * t265 + t314;
t138 = t184 * t265 + (-t228 - t232) * t242 + t312;
t137 = t183 * t242 + (-t184 - t197) * t241 + t313;
t136 = -t158 * t233 + t192 * t198 + t311;
t135 = t159 * t233 - t192 * t199 + t309;
t134 = t158 * t199 - t159 * t198 + t310;
t133 = -t149 * t209 - t160 * t233 + t165 * t170 + t193 * t198 + t311;
t132 = t150 * t209 + t161 * t233 - t166 * t170 - t193 * t199 + t309;
t131 = t149 * t166 - t150 * t165 + t160 * t199 - t161 * t198 + t310;
t1 = t276 * ((-t215 * t332 - t217 * t261 + t219 * t262) * t277 + (-t214 * t332 - t216 * t261 + t218 * t262) * t276 + (-t253 * t332 - t254 * t261 + t255 * t262) * t289) / 0.2e1 + t277 * ((t215 * t333 - t217 * t263 + t219 * t264) * t277 + (t214 * t333 - t216 * t263 + t218 * t264) * t276 + (t253 * t333 - t254 * t263 + t255 * t264) * t289) / 0.2e1 + t233 * ((t153 * t269 + t155 * t247 + t157 * t248) * t199 + (t152 * t269 + t154 * t247 + t156 * t248) * t198 + (t187 * t269 + t188 * t247 + t189 * t248) * t233) / 0.2e1 + t209 * ((t144 * t269 + t146 * t234 + t148 * t235) * t166 + (t143 * t269 + t145 * t234 + t147 * t235) * t165 + (t167 * t269 + t168 * t234 + t169 * t235) * t209) / 0.2e1 + m(1) * (t273 ^ 2 + t274 ^ 2 + t275 ^ 2) / 0.2e1 + m(2) * (t240 ^ 2 + t251 ^ 2 + t252 ^ 2) / 0.2e1 + t199 * ((t153 * t246 + t155 * t207 + t157 * t208) * t199 + (t152 * t246 + t154 * t207 + t156 * t208) * t198 + (t187 * t246 + t188 * t207 + t189 * t208) * t233) / 0.2e1 + t166 * ((t144 * t246 + t146 * t202 + t148 * t203) * t166 + (t143 * t246 + t145 * t202 + t147 * t203) * t165 + (t167 * t246 + t168 * t202 + t169 * t203) * t209) / 0.2e1 + t198 * ((t153 * t244 + t155 * t205 + t157 * t206) * t199 + (t152 * t244 + t154 * t205 + t156 * t206) * t198 + (t187 * t244 + t188 * t205 + t189 * t206) * t233) / 0.2e1 + t165 * ((t144 * t244 + t146 * t200 + t148 * t201) * t166 + (t143 * t244 + t145 * t200 + t147 * t201) * t165 + (t167 * t244 + t168 * t200 + t169 * t201) * t209) / 0.2e1 + m(3) * (t151 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(4) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(7) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + t289 * ((t214 * t276 + t215 * t277 + t253 * t289) * t303 + ((t217 * t307 + t219 * t305) * t277 + (t216 * t307 + t218 * t305) * t276 + (t254 * t307 + t255 * t305) * t289) * t301) / 0.2e1 + ((t243 * t349 + t244 * t348 + t261 * t347) * t265 + (t243 * t354 + t244 * t350 + t261 * t352) * t242 + (t355 * t243 + t351 * t244 + t353 * t261) * t241) * t241 / 0.2e1 + ((t245 * t349 + t246 * t348 + t263 * t347) * t265 + (t354 * t245 + t350 * t246 + t352 * t263) * t242 + (t245 * t355 + t246 * t351 + t263 * t353) * t241) * t242 / 0.2e1 + ((t349 * t268 + t348 * t269 - t347 * t331) * t265 + (t268 * t354 + t269 * t350 - t331 * t352) * t242 + (t268 * t355 + t269 * t351 - t331 * t353) * t241) * t265 / 0.2e1 + ((-t280 * t300 + t282 * t302 + Icges(1,4)) * V_base(5) + (-t281 * t300 + t283 * t302 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t280 * t302 + t282 * t300 + Icges(1,2)) * V_base(5) + (t281 * t302 + t283 * t300 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t302 - Icges(2,6) * t300 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t300 + Icges(2,6) * t302 + Icges(1,6)) * V_base(5) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
