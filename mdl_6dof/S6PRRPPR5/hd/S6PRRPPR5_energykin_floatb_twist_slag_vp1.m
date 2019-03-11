% Calculate kinetic energy for
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
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
% Datum: 2019-03-08 21:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRRPPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRRPPR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:18:05
% EndTime: 2019-03-08 21:18:09
% DurationCPUTime: 3.67s
% Computational Cost: add. (2619->390), mult. (5728->539), div. (0->0), fcn. (6821->12), ass. (0->172)
t362 = Icges(5,1) + Icges(4,3);
t361 = -Icges(4,4) - Icges(5,6);
t360 = Icges(5,4) - Icges(4,5);
t359 = Icges(5,5) - Icges(4,6);
t358 = Icges(4,2) + Icges(5,3);
t357 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t298 = sin(pkin(10));
t301 = cos(pkin(10));
t305 = cos(qJ(2));
t302 = cos(pkin(6));
t304 = sin(qJ(2));
t331 = t302 * t304;
t259 = t298 * t305 + t301 * t331;
t299 = sin(pkin(6));
t342 = cos(qJ(3));
t319 = t299 * t342;
t341 = sin(qJ(3));
t242 = t259 * t341 + t301 * t319;
t318 = t299 * t341;
t243 = t259 * t342 - t301 * t318;
t330 = t302 * t305;
t258 = t298 * t304 - t301 * t330;
t354 = t358 * t242 + t361 * t243 + t359 * t258;
t261 = -t298 * t331 + t301 * t305;
t244 = t261 * t341 - t298 * t319;
t245 = t261 * t342 + t298 * t318;
t260 = t298 * t330 + t301 * t304;
t353 = t358 * t244 + t361 * t245 + t359 * t260;
t352 = t359 * t242 - t360 * t243 + t362 * t258;
t351 = t359 * t244 - t360 * t245 + t362 * t260;
t265 = -t302 * t342 + t304 * t318;
t266 = t302 * t341 + t304 * t319;
t332 = t299 * t305;
t350 = t358 * t265 + t361 * t266 - t359 * t332;
t349 = t359 * t265 - t360 * t266 - t362 * t332;
t297 = sin(pkin(11));
t300 = cos(pkin(11));
t203 = t242 * t300 - t258 * t297;
t337 = t242 * t297;
t204 = t258 * t300 + t337;
t348 = Icges(6,5) * t204 + Icges(6,6) * t203 + t361 * t242 + t357 * t243 - t360 * t258;
t205 = t244 * t300 - t260 * t297;
t336 = t244 * t297;
t206 = t260 * t300 + t336;
t347 = Icges(6,5) * t206 + Icges(6,6) * t205 + t361 * t244 + t357 * t245 - t360 * t260;
t240 = t265 * t300 + t297 * t332;
t335 = t265 * t297;
t241 = -t300 * t332 + t335;
t346 = Icges(6,5) * t241 + Icges(6,6) * t240 + t361 * t265 + t357 * t266 + t360 * t332;
t340 = pkin(7) * t302;
t339 = pkin(5) * t300;
t338 = Icges(2,4) * t298;
t334 = t298 * t299;
t333 = t299 * t301;
t194 = pkin(3) * t243 + qJ(4) * t242;
t207 = pkin(4) * t258 + qJ(5) * t243;
t328 = -t194 - t207;
t195 = pkin(3) * t245 + qJ(4) * t244;
t208 = pkin(4) * t260 + qJ(5) * t245;
t327 = -t195 - t208;
t229 = pkin(3) * t266 + qJ(4) * t265;
t247 = -pkin(4) * t332 + t266 * qJ(5);
t326 = -t229 - t247;
t325 = qJD(2) * t299;
t324 = V_base(5) * qJ(1) + V_base(1);
t320 = qJD(1) + V_base(3);
t274 = t298 * t325 + V_base(4);
t286 = qJD(2) * t302 + V_base(6);
t239 = qJD(3) * t260 + t274;
t273 = -t301 * t325 + V_base(5);
t238 = qJD(3) * t258 + t273;
t262 = -qJD(3) * t332 + t286;
t268 = pkin(1) * t298 - pkin(7) * t333;
t317 = -t268 * V_base(6) + V_base(5) * t340 + t324;
t269 = pkin(1) * t301 + pkin(7) * t334;
t316 = V_base(4) * t268 - t269 * V_base(5) + t320;
t315 = V_base(6) * t269 + V_base(2) + (-qJ(1) - t340) * V_base(4);
t227 = pkin(2) * t259 + pkin(8) * t258;
t267 = (pkin(2) * t304 - pkin(8) * t305) * t299;
t314 = -t227 * t286 + t273 * t267 + t317;
t228 = pkin(2) * t261 + pkin(8) * t260;
t313 = t274 * t227 - t228 * t273 + t316;
t312 = t286 * t228 - t267 * t274 + t315;
t311 = qJD(4) * t244 + t238 * t229 + t314;
t310 = qJD(4) * t265 + t239 * t194 + t313;
t309 = qJD(4) * t242 + t262 * t195 + t312;
t308 = qJD(5) * t245 + t238 * t247 + t311;
t307 = qJD(5) * t266 + t239 * t207 + t310;
t306 = qJD(5) * t243 + t262 * t208 + t309;
t296 = pkin(11) + qJ(6);
t294 = Icges(2,4) * t301;
t293 = cos(t296);
t292 = sin(t296);
t282 = rSges(2,1) * t301 - rSges(2,2) * t298;
t281 = rSges(2,1) * t298 + rSges(2,2) * t301;
t280 = Icges(2,1) * t301 - t338;
t279 = Icges(2,1) * t298 + t294;
t278 = -Icges(2,2) * t298 + t294;
t277 = Icges(2,2) * t301 + t338;
t272 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t271 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t270 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t253 = t302 * rSges(3,3) + (rSges(3,1) * t304 + rSges(3,2) * t305) * t299;
t252 = Icges(3,5) * t302 + (Icges(3,1) * t304 + Icges(3,4) * t305) * t299;
t251 = Icges(3,6) * t302 + (Icges(3,4) * t304 + Icges(3,2) * t305) * t299;
t250 = Icges(3,3) * t302 + (Icges(3,5) * t304 + Icges(3,6) * t305) * t299;
t249 = V_base(5) * rSges(2,3) - t281 * V_base(6) + t324;
t248 = t282 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t237 = t281 * V_base(4) - t282 * V_base(5) + t320;
t232 = t265 * t292 - t293 * t332;
t231 = t265 * t293 + t292 * t332;
t230 = qJD(6) * t266 + t262;
t226 = t266 * rSges(4,1) - t265 * rSges(4,2) - rSges(4,3) * t332;
t225 = -rSges(5,1) * t332 - t266 * rSges(5,2) + t265 * rSges(5,3);
t218 = rSges(3,1) * t261 - rSges(3,2) * t260 + rSges(3,3) * t334;
t217 = rSges(3,1) * t259 - rSges(3,2) * t258 - rSges(3,3) * t333;
t216 = Icges(3,1) * t261 - Icges(3,4) * t260 + Icges(3,5) * t334;
t215 = Icges(3,1) * t259 - Icges(3,4) * t258 - Icges(3,5) * t333;
t214 = Icges(3,4) * t261 - Icges(3,2) * t260 + Icges(3,6) * t334;
t213 = Icges(3,4) * t259 - Icges(3,2) * t258 - Icges(3,6) * t333;
t212 = Icges(3,5) * t261 - Icges(3,6) * t260 + Icges(3,3) * t334;
t211 = Icges(3,5) * t259 - Icges(3,6) * t258 - Icges(3,3) * t333;
t201 = t244 * t292 + t260 * t293;
t200 = t244 * t293 - t260 * t292;
t199 = t242 * t292 + t258 * t293;
t198 = t242 * t293 - t258 * t292;
t197 = qJD(6) * t245 + t239;
t196 = qJD(6) * t243 + t238;
t191 = pkin(5) * t335 + pkin(9) * t266 - t332 * t339;
t188 = rSges(6,1) * t241 + rSges(6,2) * t240 + rSges(6,3) * t266;
t187 = rSges(4,1) * t245 - rSges(4,2) * t244 + rSges(4,3) * t260;
t186 = rSges(4,1) * t243 - rSges(4,2) * t242 + rSges(4,3) * t258;
t185 = rSges(5,1) * t260 - rSges(5,2) * t245 + rSges(5,3) * t244;
t184 = rSges(5,1) * t258 - rSges(5,2) * t243 + rSges(5,3) * t242;
t183 = Icges(6,1) * t241 + Icges(6,4) * t240 + Icges(6,5) * t266;
t182 = Icges(6,4) * t241 + Icges(6,2) * t240 + Icges(6,6) * t266;
t168 = rSges(7,1) * t232 + rSges(7,2) * t231 + rSges(7,3) * t266;
t167 = Icges(7,1) * t232 + Icges(7,4) * t231 + Icges(7,5) * t266;
t166 = Icges(7,4) * t232 + Icges(7,2) * t231 + Icges(7,6) * t266;
t165 = Icges(7,5) * t232 + Icges(7,6) * t231 + Icges(7,3) * t266;
t163 = -t217 * t286 + t253 * t273 + t317;
t162 = t218 * t286 - t253 * t274 + t315;
t161 = pkin(5) * t336 + pkin(9) * t245 + t260 * t339;
t160 = pkin(5) * t337 + pkin(9) * t243 + t258 * t339;
t159 = rSges(6,1) * t206 + rSges(6,2) * t205 + rSges(6,3) * t245;
t158 = rSges(6,1) * t204 + rSges(6,2) * t203 + rSges(6,3) * t243;
t157 = Icges(6,1) * t206 + Icges(6,4) * t205 + Icges(6,5) * t245;
t156 = Icges(6,1) * t204 + Icges(6,4) * t203 + Icges(6,5) * t243;
t155 = Icges(6,4) * t206 + Icges(6,2) * t205 + Icges(6,6) * t245;
t154 = Icges(6,4) * t204 + Icges(6,2) * t203 + Icges(6,6) * t243;
t151 = t217 * t274 - t218 * t273 + t316;
t150 = rSges(7,1) * t201 + rSges(7,2) * t200 + rSges(7,3) * t245;
t149 = rSges(7,1) * t199 + rSges(7,2) * t198 + rSges(7,3) * t243;
t148 = Icges(7,1) * t201 + Icges(7,4) * t200 + Icges(7,5) * t245;
t147 = Icges(7,1) * t199 + Icges(7,4) * t198 + Icges(7,5) * t243;
t146 = Icges(7,4) * t201 + Icges(7,2) * t200 + Icges(7,6) * t245;
t145 = Icges(7,4) * t199 + Icges(7,2) * t198 + Icges(7,6) * t243;
t144 = Icges(7,5) * t201 + Icges(7,6) * t200 + Icges(7,3) * t245;
t143 = Icges(7,5) * t199 + Icges(7,6) * t198 + Icges(7,3) * t243;
t142 = -t186 * t262 + t226 * t238 + t314;
t141 = t187 * t262 - t226 * t239 + t312;
t140 = t186 * t239 - t187 * t238 + t313;
t139 = t225 * t238 + (-t184 - t194) * t262 + t311;
t138 = t185 * t262 + (-t225 - t229) * t239 + t309;
t137 = t184 * t239 + (-t185 - t195) * t238 + t310;
t136 = t188 * t238 + (-t158 + t328) * t262 + t308;
t135 = t159 * t262 + (-t188 + t326) * t239 + t306;
t134 = t158 * t239 + (-t159 + t327) * t238 + t307;
t133 = -t149 * t230 + t168 * t196 + t191 * t238 + (-t160 + t328) * t262 + t308;
t132 = t150 * t230 + t161 * t262 - t168 * t197 + (-t191 + t326) * t239 + t306;
t131 = t149 * t197 - t150 * t196 + t160 * t239 + (-t161 + t327) * t238 + t307;
t1 = m(5) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(4) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(7) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(2) * (t237 ^ 2 + t248 ^ 2 + t249 ^ 2) / 0.2e1 + t197 * ((t144 * t245 + t146 * t200 + t148 * t201) * t197 + (t143 * t245 + t145 * t200 + t147 * t201) * t196 + (t165 * t245 + t166 * t200 + t167 * t201) * t230) / 0.2e1 + t196 * ((t144 * t243 + t146 * t198 + t148 * t199) * t197 + (t143 * t243 + t145 * t198 + t147 * t199) * t196 + (t165 * t243 + t166 * t198 + t167 * t199) * t230) / 0.2e1 + m(3) * (t151 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + t286 * ((t211 * t273 + t212 * t274 + t250 * t286) * t302 + ((t214 * t305 + t216 * t304) * t274 + (t213 * t305 + t215 * t304) * t273 + (t251 * t305 + t252 * t304) * t286) * t299) / 0.2e1 + t230 * ((t144 * t266 + t146 * t231 + t148 * t232) * t197 + (t143 * t266 + t145 * t231 + t147 * t232) * t196 + (t165 * t266 + t166 * t231 + t167 * t232) * t230) / 0.2e1 + m(1) * (t270 ^ 2 + t271 ^ 2 + t272 ^ 2) / 0.2e1 + t273 * ((-t212 * t333 - t214 * t258 + t216 * t259) * t274 + (-t211 * t333 - t213 * t258 + t215 * t259) * t273 + (-t250 * t333 - t251 * t258 + t252 * t259) * t286) / 0.2e1 + t274 * ((t212 * t334 - t214 * t260 + t216 * t261) * t274 + (t211 * t334 - t213 * t260 + t215 * t261) * t273 + (t250 * t334 - t251 * t260 + t252 * t261) * t286) / 0.2e1 + ((-t277 * t298 + t279 * t301 + Icges(1,4)) * V_base(5) + (-t278 * t298 + t280 * t301 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t277 * t301 + t279 * t298 + Icges(1,2)) * V_base(5) + (t278 * t301 + t280 * t298 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t182 * t203 + t183 * t204 + t350 * t242 + t243 * t346 + t349 * t258) * t262 + (t155 * t203 + t157 * t204 + t242 * t353 + t243 * t347 + t258 * t351) * t239 + (t154 * t203 + t156 * t204 + t354 * t242 + t348 * t243 + t352 * t258) * t238) * t238 / 0.2e1 + ((t182 * t205 + t183 * t206 + t244 * t350 + t245 * t346 + t260 * t349) * t262 + (t155 * t205 + t157 * t206 + t353 * t244 + t347 * t245 + t351 * t260) * t239 + (t154 * t205 + t156 * t206 + t244 * t354 + t245 * t348 + t260 * t352) * t238) * t239 / 0.2e1 + ((t182 * t240 + t183 * t241 + t350 * t265 + t346 * t266 - t349 * t332) * t262 + (t155 * t240 + t157 * t241 + t265 * t353 + t266 * t347 - t332 * t351) * t239 + (t154 * t240 + t156 * t241 + t265 * t354 + t266 * t348 - t332 * t352) * t238) * t262 / 0.2e1 + ((Icges(2,5) * t298 + Icges(2,6) * t301 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t301 - Icges(2,6) * t298 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
