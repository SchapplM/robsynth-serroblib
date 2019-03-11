% Calculate kinetic energy for
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-03-08 20:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRRR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:45:06
% EndTime: 2019-03-08 20:45:10
% DurationCPUTime: 4.36s
% Computational Cost: add. (2509->396), mult. (5438->571), div. (0->0), fcn. (6418->12), ass. (0->176)
t361 = Icges(3,1) + Icges(4,2);
t360 = Icges(3,4) + Icges(4,6);
t359 = Icges(3,5) - Icges(4,4);
t358 = Icges(3,2) + Icges(4,3);
t357 = Icges(3,6) - Icges(4,5);
t356 = Icges(3,3) + Icges(4,1);
t298 = sin(pkin(11));
t300 = cos(pkin(11));
t304 = sin(qJ(2));
t301 = cos(pkin(6));
t306 = cos(qJ(2));
t328 = t301 * t306;
t261 = t298 * t328 + t300 * t304;
t329 = t301 * t304;
t262 = -t298 * t329 + t300 * t306;
t299 = sin(pkin(6));
t334 = t298 * t299;
t353 = t358 * t261 - t360 * t262 - t357 * t334;
t259 = t298 * t304 - t300 * t328;
t260 = t298 * t306 + t300 * t329;
t333 = t299 * t300;
t352 = t358 * t259 - t360 * t260 + t357 * t333;
t351 = -t360 * t261 + t361 * t262 + t359 * t334;
t350 = -t360 * t259 + t361 * t260 - t359 * t333;
t349 = -t357 * t261 + t359 * t262 + t356 * t334;
t348 = -t357 * t259 + t359 * t260 - t356 * t333;
t347 = t356 * t301 + (t359 * t304 + t357 * t306) * t299;
t346 = t357 * t301 + (t360 * t304 + t358 * t306) * t299;
t345 = t359 * t301 + (t361 * t304 + t360 * t306) * t299;
t341 = cos(qJ(4));
t340 = pkin(7) * t301;
t305 = cos(qJ(5));
t339 = pkin(5) * t305;
t337 = Icges(2,4) * t298;
t302 = sin(qJ(5));
t336 = t260 * t302;
t335 = t262 * t302;
t303 = sin(qJ(4));
t332 = t299 * t303;
t331 = t299 * t304;
t330 = t299 * t306;
t327 = qJD(2) * t299;
t326 = V_base(5) * qJ(1) + V_base(1);
t322 = qJD(1) + V_base(3);
t321 = t302 * t331;
t320 = t299 * t341;
t276 = t298 * t327 + V_base(4);
t287 = qJD(2) * t301 + V_base(6);
t232 = qJD(4) * t262 + t276;
t263 = qJD(4) * t331 + t287;
t233 = -t261 * t341 + t298 * t332;
t188 = qJD(5) * t233 + t232;
t266 = t301 * t303 + t306 * t320;
t225 = qJD(5) * t266 + t263;
t275 = -t300 * t327 + V_base(5);
t231 = qJD(4) * t260 + t275;
t269 = pkin(1) * t298 - pkin(7) * t333;
t319 = -t269 * V_base(6) + V_base(5) * t340 + t326;
t270 = pkin(1) * t300 + pkin(7) * t334;
t318 = V_base(4) * t269 - t270 * V_base(5) + t322;
t235 = t259 * t341 + t300 * t332;
t187 = -qJD(5) * t235 + t231;
t317 = V_base(6) * t270 + V_base(2) + (-qJ(1) - t340) * V_base(4);
t268 = (pkin(2) * t304 - qJ(3) * t306) * t299;
t316 = qJD(3) * t261 + t275 * t268 + t319;
t222 = pkin(2) * t262 + qJ(3) * t261;
t315 = qJD(3) * t259 + t287 * t222 + t317;
t221 = pkin(2) * t260 + qJ(3) * t259;
t314 = -qJD(3) * t330 + t276 * t221 + t318;
t241 = -pkin(3) * t333 + pkin(8) * t260;
t271 = pkin(3) * t301 + pkin(8) * t331;
t313 = t275 * t271 + (-t221 - t241) * t287 + t316;
t240 = pkin(3) * t334 + pkin(8) * t262;
t312 = t287 * t240 + (-t268 - t271) * t276 + t315;
t311 = t276 * t241 + (-t222 - t240) * t275 + t314;
t236 = t259 * t303 - t300 * t320;
t186 = t236 * pkin(4) - t235 * pkin(9);
t267 = t301 * t341 - t303 * t330;
t223 = t267 * pkin(4) + t266 * pkin(9);
t310 = -t186 * t263 + t231 * t223 + t313;
t234 = t261 * t303 + t298 * t320;
t185 = t234 * pkin(4) + t233 * pkin(9);
t309 = t263 * t185 - t223 * t232 + t312;
t308 = -t185 * t231 + t232 * t186 + t311;
t297 = qJ(5) + qJ(6);
t295 = cos(t297);
t294 = sin(t297);
t293 = Icges(2,4) * t300;
t284 = rSges(2,1) * t300 - rSges(2,2) * t298;
t283 = rSges(2,1) * t298 + rSges(2,2) * t300;
t282 = Icges(2,1) * t300 - t337;
t281 = Icges(2,1) * t298 + t293;
t280 = -Icges(2,2) * t298 + t293;
t279 = Icges(2,2) * t300 + t337;
t274 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t273 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t272 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t252 = rSges(4,1) * t301 + (-rSges(4,2) * t304 - rSges(4,3) * t306) * t299;
t251 = rSges(3,3) * t301 + (rSges(3,1) * t304 + rSges(3,2) * t306) * t299;
t243 = V_base(5) * rSges(2,3) - t283 * V_base(6) + t326;
t242 = t284 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t238 = t267 * t305 + t321;
t237 = -t267 * t302 + t305 * t331;
t230 = t283 * V_base(4) - t284 * V_base(5) + t322;
t227 = t267 * t295 + t294 * t331;
t226 = -t267 * t294 + t295 * t331;
t219 = rSges(5,1) * t267 - rSges(5,2) * t266 + rSges(5,3) * t331;
t218 = Icges(5,1) * t267 - Icges(5,4) * t266 + Icges(5,5) * t331;
t217 = Icges(5,4) * t267 - Icges(5,2) * t266 + Icges(5,6) * t331;
t216 = Icges(5,5) * t267 - Icges(5,6) * t266 + Icges(5,3) * t331;
t215 = rSges(3,1) * t262 - rSges(3,2) * t261 + rSges(3,3) * t334;
t214 = rSges(3,1) * t260 - rSges(3,2) * t259 - rSges(3,3) * t333;
t213 = -rSges(4,1) * t333 - rSges(4,2) * t260 + rSges(4,3) * t259;
t212 = rSges(4,1) * t334 - rSges(4,2) * t262 + rSges(4,3) * t261;
t197 = qJD(6) * t266 + t225;
t196 = t236 * t305 + t336;
t195 = -t236 * t302 + t260 * t305;
t194 = t234 * t305 + t335;
t193 = -t234 * t302 + t262 * t305;
t192 = t236 * t295 + t260 * t294;
t191 = -t236 * t294 + t260 * t295;
t190 = t234 * t295 + t262 * t294;
t189 = -t234 * t294 + t262 * t295;
t183 = rSges(6,1) * t238 + rSges(6,2) * t237 + rSges(6,3) * t266;
t181 = Icges(6,1) * t238 + Icges(6,4) * t237 + Icges(6,5) * t266;
t180 = Icges(6,4) * t238 + Icges(6,2) * t237 + Icges(6,6) * t266;
t179 = Icges(6,5) * t238 + Icges(6,6) * t237 + Icges(6,3) * t266;
t178 = rSges(5,1) * t236 + rSges(5,2) * t235 + rSges(5,3) * t260;
t177 = rSges(5,1) * t234 - rSges(5,2) * t233 + rSges(5,3) * t262;
t176 = Icges(5,1) * t236 + Icges(5,4) * t235 + Icges(5,5) * t260;
t175 = Icges(5,1) * t234 - Icges(5,4) * t233 + Icges(5,5) * t262;
t174 = Icges(5,4) * t236 + Icges(5,2) * t235 + Icges(5,6) * t260;
t173 = Icges(5,4) * t234 - Icges(5,2) * t233 + Icges(5,6) * t262;
t172 = Icges(5,5) * t236 + Icges(5,6) * t235 + Icges(5,3) * t260;
t171 = Icges(5,5) * t234 - Icges(5,6) * t233 + Icges(5,3) * t262;
t170 = pkin(5) * t321 + pkin(10) * t266 + t267 * t339;
t169 = rSges(7,1) * t227 + rSges(7,2) * t226 + rSges(7,3) * t266;
t168 = Icges(7,1) * t227 + Icges(7,4) * t226 + Icges(7,5) * t266;
t167 = Icges(7,4) * t227 + Icges(7,2) * t226 + Icges(7,6) * t266;
t166 = Icges(7,5) * t227 + Icges(7,6) * t226 + Icges(7,3) * t266;
t165 = qJD(6) * t233 + t188;
t164 = -qJD(6) * t235 + t187;
t162 = -t214 * t287 + t251 * t275 + t319;
t161 = t215 * t287 - t251 * t276 + t317;
t160 = rSges(6,1) * t196 + rSges(6,2) * t195 - rSges(6,3) * t235;
t159 = rSges(6,1) * t194 + rSges(6,2) * t193 + rSges(6,3) * t233;
t158 = Icges(6,1) * t196 + Icges(6,4) * t195 - Icges(6,5) * t235;
t157 = Icges(6,1) * t194 + Icges(6,4) * t193 + Icges(6,5) * t233;
t156 = Icges(6,4) * t196 + Icges(6,2) * t195 - Icges(6,6) * t235;
t155 = Icges(6,4) * t194 + Icges(6,2) * t193 + Icges(6,6) * t233;
t154 = Icges(6,5) * t196 + Icges(6,6) * t195 - Icges(6,3) * t235;
t153 = Icges(6,5) * t194 + Icges(6,6) * t193 + Icges(6,3) * t233;
t152 = t214 * t276 - t215 * t275 + t318;
t151 = rSges(7,1) * t192 + rSges(7,2) * t191 - rSges(7,3) * t235;
t150 = rSges(7,1) * t190 + rSges(7,2) * t189 + rSges(7,3) * t233;
t149 = Icges(7,1) * t192 + Icges(7,4) * t191 - Icges(7,5) * t235;
t148 = Icges(7,1) * t190 + Icges(7,4) * t189 + Icges(7,5) * t233;
t147 = Icges(7,4) * t192 + Icges(7,2) * t191 - Icges(7,6) * t235;
t146 = Icges(7,4) * t190 + Icges(7,2) * t189 + Icges(7,6) * t233;
t145 = Icges(7,5) * t192 + Icges(7,6) * t191 - Icges(7,3) * t235;
t144 = Icges(7,5) * t190 + Icges(7,6) * t189 + Icges(7,3) * t233;
t143 = pkin(5) * t336 - pkin(10) * t235 + t236 * t339;
t142 = pkin(5) * t335 + pkin(10) * t233 + t234 * t339;
t141 = t252 * t275 + (-t213 - t221) * t287 + t316;
t140 = t212 * t287 + (-t252 - t268) * t276 + t315;
t139 = t213 * t276 + (-t212 - t222) * t275 + t314;
t138 = -t178 * t263 + t219 * t231 + t313;
t137 = t177 * t263 - t219 * t232 + t312;
t136 = -t177 * t231 + t178 * t232 + t311;
t135 = -t160 * t225 + t183 * t187 + t310;
t134 = t159 * t225 - t183 * t188 + t309;
t133 = -t159 * t187 + t160 * t188 + t308;
t132 = -t143 * t225 - t151 * t197 + t164 * t169 + t170 * t187 + t310;
t131 = t142 * t225 + t150 * t197 - t165 * t169 - t170 * t188 + t309;
t130 = -t142 * t187 + t143 * t188 - t150 * t164 + t151 * t165 + t308;
t1 = t263 * ((t171 * t331 - t173 * t266 + t175 * t267) * t232 + (t172 * t331 - t174 * t266 + t176 * t267) * t231 + (t216 * t331 - t266 * t217 + t267 * t218) * t263) / 0.2e1 + m(3) * (t152 ^ 2 + t161 ^ 2 + t162 ^ 2) / 0.2e1 + m(4) * (t139 ^ 2 + t140 ^ 2 + t141 ^ 2) / 0.2e1 + m(6) * (t133 ^ 2 + t134 ^ 2 + t135 ^ 2) / 0.2e1 + m(5) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(7) * (t130 ^ 2 + t131 ^ 2 + t132 ^ 2) / 0.2e1 + t165 * ((t233 * t144 + t189 * t146 + t190 * t148) * t165 + (t145 * t233 + t147 * t189 + t149 * t190) * t164 + (t166 * t233 + t167 * t189 + t168 * t190) * t197) / 0.2e1 + t188 * ((t233 * t153 + t193 * t155 + t194 * t157) * t188 + (t154 * t233 + t156 * t193 + t158 * t194) * t187 + (t179 * t233 + t180 * t193 + t181 * t194) * t225) / 0.2e1 + t164 * ((-t144 * t235 + t146 * t191 + t148 * t192) * t165 + (-t235 * t145 + t191 * t147 + t192 * t149) * t164 + (-t166 * t235 + t167 * t191 + t168 * t192) * t197) / 0.2e1 + t187 * ((-t153 * t235 + t155 * t195 + t157 * t196) * t188 + (-t235 * t154 + t195 * t156 + t196 * t158) * t187 + (-t179 * t235 + t180 * t195 + t181 * t196) * t225) / 0.2e1 + m(2) * (t230 ^ 2 + t242 ^ 2 + t243 ^ 2) / 0.2e1 + t231 * ((t171 * t260 + t173 * t235 + t175 * t236) * t232 + (t260 * t172 + t235 * t174 + t236 * t176) * t231 + (t216 * t260 + t217 * t235 + t218 * t236) * t263) / 0.2e1 + t232 * ((t262 * t171 - t233 * t173 + t234 * t175) * t232 + (t172 * t262 - t174 * t233 + t176 * t234) * t231 + (t216 * t262 - t217 * t233 + t218 * t234) * t263) / 0.2e1 + t197 * ((t144 * t266 + t146 * t226 + t148 * t227) * t165 + (t145 * t266 + t147 * t226 + t149 * t227) * t164 + (t266 * t166 + t226 * t167 + t227 * t168) * t197) / 0.2e1 + t225 * ((t153 * t266 + t155 * t237 + t157 * t238) * t188 + (t154 * t266 + t156 * t237 + t158 * t238) * t187 + (t266 * t179 + t237 * t180 + t238 * t181) * t225) / 0.2e1 + m(1) * (t272 ^ 2 + t273 ^ 2 + t274 ^ 2) / 0.2e1 + ((-t259 * t346 + t260 * t345 - t333 * t347) * t287 + (t259 * t353 + t260 * t351 - t333 * t349) * t276 + (t352 * t259 + t350 * t260 - t348 * t333) * t275) * t275 / 0.2e1 + ((-t261 * t346 + t262 * t345 + t334 * t347) * t287 + (t353 * t261 + t351 * t262 + t349 * t334) * t276 + (t261 * t352 + t262 * t350 + t334 * t348) * t275) * t276 / 0.2e1 + ((t275 * t348 + t276 * t349 + t287 * t347) * t301 + ((t304 * t345 + t306 * t346) * t287 + (t304 * t351 - t306 * t353) * t276 + (t304 * t350 - t352 * t306) * t275) * t299) * t287 / 0.2e1 + ((-t279 * t298 + t281 * t300 + Icges(1,4)) * V_base(5) + (-t298 * t280 + t300 * t282 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t300 * t279 + t298 * t281 + Icges(1,2)) * V_base(5) + (t280 * t300 + t282 * t298 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t298 + Icges(2,6) * t300 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t300 - Icges(2,6) * t298 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
