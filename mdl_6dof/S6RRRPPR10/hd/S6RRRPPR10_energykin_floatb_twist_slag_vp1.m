% Calculate kinetic energy for
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
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
% Datum: 2019-03-09 16:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPPR10_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPPR10_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR10_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR10_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPPR10_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:21:23
% EndTime: 2019-03-09 16:21:26
% DurationCPUTime: 3.84s
% Computational Cost: add. (2679->390), mult. (5728->541), div. (0->0), fcn. (6821->12), ass. (0->174)
t361 = Icges(5,1) + Icges(4,3);
t360 = -Icges(4,4) - Icges(5,6);
t359 = Icges(5,4) - Icges(4,5);
t358 = Icges(5,5) - Icges(4,6);
t357 = Icges(4,2) + Icges(5,3);
t356 = Icges(4,1) + Icges(5,2) + Icges(6,3);
t301 = cos(pkin(6));
t304 = sin(qJ(1));
t305 = cos(qJ(2));
t330 = t304 * t305;
t303 = sin(qJ(2));
t306 = cos(qJ(1));
t331 = t303 * t306;
t263 = t301 * t331 + t330;
t299 = sin(pkin(6));
t343 = cos(qJ(3));
t320 = t299 * t343;
t342 = sin(qJ(3));
t241 = t263 * t342 + t306 * t320;
t319 = t299 * t342;
t242 = t263 * t343 - t306 * t319;
t329 = t305 * t306;
t332 = t303 * t304;
t262 = -t301 * t329 + t332;
t355 = t241 * t357 + t242 * t360 + t262 * t358;
t265 = -t301 * t332 + t329;
t243 = t265 * t342 - t304 * t320;
t244 = t265 * t343 + t304 * t319;
t264 = t301 * t330 + t331;
t354 = t243 * t357 + t244 * t360 + t264 * t358;
t353 = t358 * t241 - t359 * t242 + t262 * t361;
t352 = t358 * t243 - t359 * t244 + t264 * t361;
t260 = -t301 * t343 + t303 * t319;
t261 = t301 * t342 + t303 * t320;
t334 = t299 * t305;
t351 = t260 * t357 + t261 * t360 - t334 * t358;
t350 = t358 * t260 - t359 * t261 - t334 * t361;
t298 = sin(pkin(11));
t300 = cos(pkin(11));
t203 = t241 * t300 - t262 * t298;
t338 = t241 * t298;
t204 = t262 * t300 + t338;
t349 = Icges(6,5) * t204 + Icges(6,6) * t203 + t241 * t360 + t242 * t356 - t262 * t359;
t205 = t243 * t300 - t264 * t298;
t337 = t243 * t298;
t206 = t264 * t300 + t337;
t348 = Icges(6,5) * t206 + Icges(6,6) * t205 + t243 * t360 + t244 * t356 - t264 * t359;
t237 = t260 * t300 + t298 * t334;
t336 = t260 * t298;
t238 = -t300 * t334 + t336;
t347 = Icges(6,5) * t238 + Icges(6,6) * t237 + t260 * t360 + t261 * t356 + t334 * t359;
t341 = pkin(8) * t301;
t340 = pkin(5) * t300;
t339 = Icges(2,4) * t304;
t335 = t299 * t304;
t333 = t299 * t306;
t194 = pkin(3) * t242 + qJ(4) * t241;
t207 = pkin(4) * t262 + qJ(5) * t242;
t327 = -t194 - t207;
t195 = pkin(3) * t244 + qJ(4) * t243;
t208 = pkin(4) * t264 + qJ(5) * t244;
t326 = -t195 - t208;
t227 = pkin(3) * t261 + qJ(4) * t260;
t249 = -pkin(4) * t334 + qJ(5) * t261;
t325 = -t227 - t249;
t324 = qJD(2) * t299;
t323 = V_base(5) * pkin(7) + V_base(1);
t274 = t304 * t324 + V_base(4);
t294 = V_base(6) + qJD(1);
t240 = qJD(3) * t264 + t274;
t275 = qJD(2) * t301 + t294;
t273 = -t306 * t324 + V_base(5);
t268 = t304 * pkin(1) - pkin(8) * t333;
t318 = -t268 * t294 + V_base(5) * t341 + t323;
t269 = pkin(1) * t306 + pkin(8) * t335;
t317 = V_base(4) * t268 - t269 * V_base(5) + V_base(3);
t239 = qJD(3) * t262 + t273;
t258 = -qJD(3) * t334 + t275;
t316 = t294 * t269 + V_base(2) + (-pkin(7) - t341) * V_base(4);
t229 = pkin(2) * t263 + pkin(9) * t262;
t267 = (pkin(2) * t303 - pkin(9) * t305) * t299;
t315 = -t229 * t275 + t273 * t267 + t318;
t230 = pkin(2) * t265 + pkin(9) * t264;
t314 = t274 * t229 - t230 * t273 + t317;
t313 = qJD(4) * t243 + t239 * t227 + t315;
t312 = qJD(4) * t260 + t240 * t194 + t314;
t311 = t275 * t230 - t267 * t274 + t316;
t310 = qJD(5) * t244 + t239 * t249 + t313;
t309 = qJD(5) * t261 + t240 * t207 + t312;
t308 = qJD(4) * t241 + t258 * t195 + t311;
t307 = qJD(5) * t242 + t258 * t208 + t308;
t297 = pkin(11) + qJ(6);
t295 = Icges(2,4) * t306;
t293 = cos(t297);
t292 = sin(t297);
t283 = rSges(2,1) * t306 - t304 * rSges(2,2);
t282 = t304 * rSges(2,1) + rSges(2,2) * t306;
t281 = Icges(2,1) * t306 - t339;
t280 = Icges(2,1) * t304 + t295;
t279 = -Icges(2,2) * t304 + t295;
t278 = Icges(2,2) * t306 + t339;
t272 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t271 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t270 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t253 = rSges(3,3) * t301 + (rSges(3,1) * t303 + rSges(3,2) * t305) * t299;
t252 = Icges(3,5) * t301 + (Icges(3,1) * t303 + Icges(3,4) * t305) * t299;
t251 = Icges(3,6) * t301 + (Icges(3,4) * t303 + Icges(3,2) * t305) * t299;
t250 = Icges(3,3) * t301 + (Icges(3,5) * t303 + Icges(3,6) * t305) * t299;
t248 = V_base(5) * rSges(2,3) - t282 * t294 + t323;
t247 = t283 * t294 + V_base(2) + (-rSges(2,3) - pkin(7)) * V_base(4);
t245 = t282 * V_base(4) - t283 * V_base(5) + V_base(3);
t232 = t260 * t292 - t293 * t334;
t231 = t260 * t293 + t292 * t334;
t228 = qJD(6) * t261 + t258;
t226 = rSges(3,1) * t265 - rSges(3,2) * t264 + rSges(3,3) * t335;
t225 = t263 * rSges(3,1) - t262 * rSges(3,2) - rSges(3,3) * t333;
t224 = Icges(3,1) * t265 - Icges(3,4) * t264 + Icges(3,5) * t335;
t223 = Icges(3,1) * t263 - Icges(3,4) * t262 - Icges(3,5) * t333;
t222 = Icges(3,4) * t265 - Icges(3,2) * t264 + Icges(3,6) * t335;
t221 = Icges(3,4) * t263 - Icges(3,2) * t262 - Icges(3,6) * t333;
t220 = Icges(3,5) * t265 - Icges(3,6) * t264 + Icges(3,3) * t335;
t219 = Icges(3,5) * t263 - Icges(3,6) * t262 - Icges(3,3) * t333;
t218 = rSges(4,1) * t261 - rSges(4,2) * t260 - rSges(4,3) * t334;
t217 = -rSges(5,1) * t334 - rSges(5,2) * t261 + rSges(5,3) * t260;
t201 = t243 * t292 + t264 * t293;
t200 = t243 * t293 - t264 * t292;
t199 = t241 * t292 + t262 * t293;
t198 = t241 * t293 - t262 * t292;
t197 = qJD(6) * t244 + t240;
t196 = qJD(6) * t242 + t239;
t190 = pkin(5) * t336 + pkin(10) * t261 - t334 * t340;
t189 = rSges(4,1) * t244 - rSges(4,2) * t243 + rSges(4,3) * t264;
t188 = rSges(4,1) * t242 - rSges(4,2) * t241 + rSges(4,3) * t262;
t187 = rSges(5,1) * t264 - rSges(5,2) * t244 + rSges(5,3) * t243;
t186 = rSges(5,1) * t262 - rSges(5,2) * t242 + rSges(5,3) * t241;
t172 = rSges(6,1) * t238 + rSges(6,2) * t237 + rSges(6,3) * t261;
t171 = Icges(6,1) * t238 + Icges(6,4) * t237 + Icges(6,5) * t261;
t170 = Icges(6,4) * t238 + Icges(6,2) * t237 + Icges(6,6) * t261;
t168 = rSges(7,1) * t232 + rSges(7,2) * t231 + rSges(7,3) * t261;
t167 = Icges(7,1) * t232 + Icges(7,4) * t231 + Icges(7,5) * t261;
t166 = Icges(7,4) * t232 + Icges(7,2) * t231 + Icges(7,6) * t261;
t165 = Icges(7,5) * t232 + Icges(7,6) * t231 + Icges(7,3) * t261;
t163 = -t225 * t275 + t253 * t273 + t318;
t162 = t226 * t275 - t253 * t274 + t316;
t161 = pkin(5) * t337 + pkin(10) * t244 + t264 * t340;
t160 = pkin(5) * t338 + pkin(10) * t242 + t262 * t340;
t159 = rSges(6,1) * t206 + rSges(6,2) * t205 + rSges(6,3) * t244;
t158 = rSges(6,1) * t204 + rSges(6,2) * t203 + rSges(6,3) * t242;
t157 = Icges(6,1) * t206 + Icges(6,4) * t205 + Icges(6,5) * t244;
t156 = Icges(6,1) * t204 + Icges(6,4) * t203 + Icges(6,5) * t242;
t155 = Icges(6,4) * t206 + Icges(6,2) * t205 + Icges(6,6) * t244;
t154 = Icges(6,4) * t204 + Icges(6,2) * t203 + Icges(6,6) * t242;
t151 = t225 * t274 - t226 * t273 + t317;
t150 = rSges(7,1) * t201 + rSges(7,2) * t200 + rSges(7,3) * t244;
t149 = rSges(7,1) * t199 + rSges(7,2) * t198 + rSges(7,3) * t242;
t148 = Icges(7,1) * t201 + Icges(7,4) * t200 + Icges(7,5) * t244;
t147 = Icges(7,1) * t199 + Icges(7,4) * t198 + Icges(7,5) * t242;
t146 = Icges(7,4) * t201 + Icges(7,2) * t200 + Icges(7,6) * t244;
t145 = Icges(7,4) * t199 + Icges(7,2) * t198 + Icges(7,6) * t242;
t144 = Icges(7,5) * t201 + Icges(7,6) * t200 + Icges(7,3) * t244;
t143 = Icges(7,5) * t199 + Icges(7,6) * t198 + Icges(7,3) * t242;
t142 = -t188 * t258 + t218 * t239 + t315;
t141 = t189 * t258 - t218 * t240 + t311;
t140 = t188 * t240 - t189 * t239 + t314;
t139 = t217 * t239 + (-t186 - t194) * t258 + t313;
t138 = t187 * t258 + (-t217 - t227) * t240 + t308;
t137 = t186 * t240 + (-t187 - t195) * t239 + t312;
t136 = t172 * t239 + (-t158 + t327) * t258 + t310;
t135 = t159 * t258 + (-t172 + t325) * t240 + t307;
t134 = t158 * t240 + (-t159 + t326) * t239 + t309;
t133 = t310 - t149 * t228 + t168 * t196 + t190 * t239 + (-t160 + t327) * t258;
t132 = t150 * t228 + t161 * t258 - t168 * t197 + (-t190 + t325) * t240 + t307;
t131 = t149 * t197 - t150 * t196 + t160 * t240 + (-t161 + t326) * t239 + t309;
t1 = m(1) * (t270 ^ 2 + t271 ^ 2 + t272 ^ 2) / 0.2e1 + t228 * ((t144 * t261 + t146 * t231 + t148 * t232) * t197 + (t143 * t261 + t145 * t231 + t147 * t232) * t196 + (t165 * t261 + t166 * t231 + t167 * t232) * t228) / 0.2e1 + t197 * ((t144 * t244 + t146 * t200 + t148 * t201) * t197 + (t143 * t244 + t145 * t200 + t147 * t201) * t196 + (t165 * t244 + t166 * t200 + t167 * t201) * t228) / 0.2e1 + m(2) * (t245 ^ 2 + t247 ^ 2 + t248 ^ 2) / 0.2e1 + t196 * ((t144 * t242 + t146 * t198 + t148 * t199) * t197 + (t143 * t242 + t145 * t198 + t147 * t199) * t196 + (t165 * t242 + t166 * t198 + t167 * t199) * t228) / 0.2e1 + m(3) * (t151 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(5) * (t137 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(4) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(6) * (t134 ^ 2 + t135 ^ 2 + t136 ^ 2) / 0.2e1 + m(7) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + t274 * ((t220 * t335 - t222 * t264 + t224 * t265) * t274 + (t219 * t335 - t221 * t264 + t223 * t265) * t273 + (t250 * t335 - t251 * t264 + t252 * t265) * t275) / 0.2e1 + t273 * ((-t220 * t333 - t262 * t222 + t263 * t224) * t274 + (-t219 * t333 - t262 * t221 + t263 * t223) * t273 + (-t250 * t333 - t262 * t251 + t263 * t252) * t275) / 0.2e1 + t275 * ((t219 * t273 + t220 * t274 + t250 * t275) * t301 + ((t222 * t305 + t224 * t303) * t274 + (t221 * t305 + t223 * t303) * t273 + (t251 * t305 + t252 * t303) * t275) * t299) / 0.2e1 + ((-t304 * t278 + t280 * t306 + Icges(1,4)) * V_base(5) + (-t304 * t279 + t281 * t306 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t278 * t306 + t304 * t280 + Icges(1,2)) * V_base(5) + (t279 * t306 + t304 * t281 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t170 * t203 + t171 * t204 + t241 * t351 + t242 * t347 + t262 * t350) * t258 + (t155 * t203 + t157 * t204 + t241 * t354 + t242 * t348 + t262 * t352) * t240 + (t154 * t203 + t156 * t204 + t355 * t241 + t349 * t242 + t353 * t262) * t239) * t239 / 0.2e1 + ((t170 * t205 + t171 * t206 + t243 * t351 + t244 * t347 + t264 * t350) * t258 + (t155 * t205 + t157 * t206 + t243 * t354 + t244 * t348 + t264 * t352) * t240 + (t154 * t205 + t156 * t206 + t243 * t355 + t349 * t244 + t353 * t264) * t239) * t240 / 0.2e1 + ((t170 * t237 + t171 * t238 + t260 * t351 + t261 * t347 - t334 * t350) * t258 + (t155 * t237 + t157 * t238 + t260 * t354 + t261 * t348 - t334 * t352) * t240 + (t154 * t237 + t156 * t238 + t260 * t355 + t349 * t261 - t353 * t334) * t239) * t258 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t304 + Icges(2,6) * t306) * V_base(5) + (Icges(2,5) * t306 - Icges(2,6) * t304) * V_base(4) + Icges(2,3) * t294 / 0.2e1) * t294;
T  = t1;
