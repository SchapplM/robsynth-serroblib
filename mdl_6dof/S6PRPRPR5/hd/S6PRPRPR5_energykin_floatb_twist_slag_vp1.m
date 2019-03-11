% Calculate kinetic energy for
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6PRPRPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6PRPRPR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:42:52
% EndTime: 2019-03-08 19:42:56
% DurationCPUTime: 4.42s
% Computational Cost: add. (2937->396), mult. (5136->550), div. (0->0), fcn. (5992->12), ass. (0->176)
t363 = Icges(5,1) + Icges(6,2);
t362 = Icges(6,1) + Icges(5,3);
t361 = -Icges(5,4) - Icges(6,6);
t360 = Icges(6,4) - Icges(5,5);
t359 = Icges(6,5) - Icges(5,6);
t358 = Icges(5,2) + Icges(6,3);
t292 = sin(pkin(10));
t295 = cos(pkin(10));
t301 = cos(qJ(2));
t296 = cos(pkin(6));
t299 = sin(qJ(2));
t329 = t296 * t299;
t258 = t292 * t301 + t295 * t329;
t293 = sin(pkin(6));
t325 = pkin(11) + qJ(4);
t317 = cos(t325);
t315 = t293 * t317;
t316 = sin(t325);
t227 = t258 * t316 + t295 * t315;
t314 = t293 * t316;
t228 = t258 * t317 - t295 * t314;
t328 = t296 * t301;
t257 = t292 * t299 - t295 * t328;
t356 = t227 * t358 + t228 * t361 + t257 * t359;
t260 = -t292 * t329 + t295 * t301;
t229 = t260 * t316 - t292 * t315;
t230 = t260 * t317 + t292 * t314;
t259 = t292 * t328 + t295 * t299;
t355 = t229 * t358 + t230 * t361 + t259 * t359;
t354 = t227 * t359 - t228 * t360 + t257 * t362;
t353 = t229 * t359 - t230 * t360 + t259 * t362;
t352 = t361 * t227 + t228 * t363 - t360 * t257;
t351 = t361 * t229 + t230 * t363 - t360 * t259;
t291 = sin(pkin(11));
t294 = cos(pkin(11));
t332 = t293 * t295;
t234 = -t258 * t291 - t294 * t332;
t318 = t291 * t332;
t235 = t258 * t294 - t318;
t178 = Icges(4,5) * t235 + Icges(4,6) * t234 + Icges(4,3) * t257;
t214 = Icges(3,4) * t258 - Icges(3,2) * t257 - Icges(3,6) * t332;
t350 = t178 - t214;
t333 = t292 * t293;
t236 = -t260 * t291 + t294 * t333;
t319 = t291 * t333;
t237 = t260 * t294 + t319;
t179 = Icges(4,5) * t237 + Icges(4,6) * t236 + Icges(4,3) * t259;
t215 = Icges(3,4) * t260 - Icges(3,2) * t259 + Icges(3,6) * t333;
t349 = t179 - t215;
t248 = -t296 * t317 + t299 * t314;
t249 = t296 * t316 + t299 * t315;
t330 = t293 * t301;
t348 = t248 * t358 + t249 * t361 - t330 * t359;
t347 = t361 * t248 + t249 * t363 + t360 * t330;
t346 = t248 * t359 - t249 * t360 - t330 * t362;
t331 = t293 * t299;
t255 = -t291 * t331 + t294 * t296;
t334 = t291 * t296;
t256 = t294 * t331 + t334;
t208 = Icges(4,5) * t256 + Icges(4,6) * t255 - Icges(4,3) * t330;
t246 = Icges(3,6) * t296 + (Icges(3,4) * t299 + Icges(3,2) * t301) * t293;
t345 = t208 - t246;
t337 = pkin(7) * t296;
t336 = pkin(3) * t294;
t335 = Icges(2,4) * t292;
t326 = qJD(2) * t293;
t324 = V_base(5) * qJ(1) + V_base(1);
t320 = qJD(1) + V_base(3);
t273 = t292 * t326 + V_base(4);
t284 = qJD(2) * t296 + V_base(6);
t240 = qJD(4) * t259 + t273;
t272 = -t295 * t326 + V_base(5);
t239 = qJD(4) * t257 + t272;
t261 = -qJD(4) * t330 + t284;
t265 = pkin(1) * t292 - pkin(7) * t332;
t313 = -t265 * V_base(6) + V_base(5) * t337 + t324;
t266 = pkin(1) * t295 + pkin(7) * t333;
t312 = V_base(4) * t265 - V_base(5) * t266 + t320;
t311 = V_base(6) * t266 + V_base(2) + (-qJ(1) - t337) * V_base(4);
t264 = (pkin(2) * t299 - qJ(3) * t301) * t293;
t310 = qJD(3) * t259 + t272 * t264 + t313;
t224 = pkin(2) * t260 + qJ(3) * t259;
t309 = qJD(3) * t257 + t284 * t224 + t311;
t223 = pkin(2) * t258 + qJ(3) * t257;
t308 = -qJD(3) * t330 + t273 * t223 + t312;
t175 = -pkin(3) * t318 + pkin(8) * t257 + t258 * t336;
t219 = pkin(3) * t334 + (-pkin(8) * t301 + t299 * t336) * t293;
t307 = t272 * t219 + (-t175 - t223) * t284 + t310;
t176 = pkin(3) * t319 + pkin(8) * t259 + t260 * t336;
t306 = t284 * t176 + (-t219 - t264) * t273 + t309;
t211 = pkin(4) * t249 + qJ(5) * t248;
t305 = qJD(5) * t229 + t239 * t211 + t307;
t304 = t273 * t175 + (-t176 - t224) * t272 + t308;
t188 = pkin(4) * t230 + qJ(5) * t229;
t303 = qJD(5) * t227 + t261 * t188 + t306;
t187 = pkin(4) * t228 + qJ(5) * t227;
t302 = qJD(5) * t248 + t240 * t187 + t304;
t300 = cos(qJ(6));
t298 = sin(qJ(6));
t289 = Icges(2,4) * t295;
t281 = rSges(2,1) * t295 - rSges(2,2) * t292;
t280 = rSges(2,1) * t292 + rSges(2,2) * t295;
t279 = Icges(2,1) * t295 - t335;
t278 = Icges(2,1) * t292 + t289;
t277 = -Icges(2,2) * t292 + t289;
t276 = Icges(2,2) * t295 + t335;
t271 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t270 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t269 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t250 = t296 * rSges(3,3) + (rSges(3,1) * t299 + rSges(3,2) * t301) * t293;
t247 = Icges(3,5) * t296 + (Icges(3,1) * t299 + Icges(3,4) * t301) * t293;
t245 = Icges(3,3) * t296 + (Icges(3,5) * t299 + Icges(3,6) * t301) * t293;
t243 = V_base(5) * rSges(2,3) - t280 * V_base(6) + t324;
t242 = t281 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t238 = -pkin(5) * t330 + t249 * pkin(9);
t233 = t280 * V_base(4) - t281 * V_base(5) + t320;
t232 = t248 * t298 - t300 * t330;
t231 = t248 * t300 + t298 * t330;
t222 = qJD(6) * t249 + t261;
t221 = rSges(3,1) * t260 - rSges(3,2) * t259 + rSges(3,3) * t333;
t220 = rSges(3,1) * t258 - rSges(3,2) * t257 - rSges(3,3) * t332;
t218 = t256 * rSges(4,1) + t255 * rSges(4,2) - rSges(4,3) * t330;
t217 = Icges(3,1) * t260 - Icges(3,4) * t259 + Icges(3,5) * t333;
t216 = Icges(3,1) * t258 - Icges(3,4) * t257 - Icges(3,5) * t332;
t213 = Icges(3,5) * t260 - Icges(3,6) * t259 + Icges(3,3) * t333;
t212 = Icges(3,5) * t258 - Icges(3,6) * t257 - Icges(3,3) * t332;
t210 = Icges(4,1) * t256 + Icges(4,4) * t255 - Icges(4,5) * t330;
t209 = Icges(4,4) * t256 + Icges(4,2) * t255 - Icges(4,6) * t330;
t205 = t249 * rSges(5,1) - t248 * rSges(5,2) - rSges(5,3) * t330;
t204 = -rSges(6,1) * t330 - t249 * rSges(6,2) + t248 * rSges(6,3);
t197 = pkin(5) * t259 + pkin(9) * t230;
t196 = pkin(5) * t257 + pkin(9) * t228;
t194 = t229 * t298 + t259 * t300;
t193 = t229 * t300 - t259 * t298;
t192 = t227 * t298 + t257 * t300;
t191 = t227 * t300 - t257 * t298;
t190 = qJD(6) * t230 + t240;
t189 = qJD(6) * t228 + t239;
t185 = rSges(4,1) * t237 + rSges(4,2) * t236 + rSges(4,3) * t259;
t184 = rSges(4,1) * t235 + rSges(4,2) * t234 + rSges(4,3) * t257;
t183 = Icges(4,1) * t237 + Icges(4,4) * t236 + Icges(4,5) * t259;
t182 = Icges(4,1) * t235 + Icges(4,4) * t234 + Icges(4,5) * t257;
t181 = Icges(4,4) * t237 + Icges(4,2) * t236 + Icges(4,6) * t259;
t180 = Icges(4,4) * t235 + Icges(4,2) * t234 + Icges(4,6) * t257;
t174 = rSges(5,1) * t230 - rSges(5,2) * t229 + rSges(5,3) * t259;
t173 = rSges(5,1) * t228 - rSges(5,2) * t227 + rSges(5,3) * t257;
t172 = rSges(6,1) * t259 - rSges(6,2) * t230 + rSges(6,3) * t229;
t171 = rSges(6,1) * t257 - rSges(6,2) * t228 + rSges(6,3) * t227;
t158 = rSges(7,1) * t232 + rSges(7,2) * t231 + rSges(7,3) * t249;
t157 = Icges(7,1) * t232 + Icges(7,4) * t231 + Icges(7,5) * t249;
t156 = Icges(7,4) * t232 + Icges(7,2) * t231 + Icges(7,6) * t249;
t155 = Icges(7,5) * t232 + Icges(7,6) * t231 + Icges(7,3) * t249;
t151 = -t220 * t284 + t250 * t272 + t313;
t150 = t221 * t284 - t250 * t273 + t311;
t149 = t220 * t273 - t221 * t272 + t312;
t148 = rSges(7,1) * t194 + rSges(7,2) * t193 + rSges(7,3) * t230;
t147 = rSges(7,1) * t192 + rSges(7,2) * t191 + rSges(7,3) * t228;
t146 = Icges(7,1) * t194 + Icges(7,4) * t193 + Icges(7,5) * t230;
t145 = Icges(7,1) * t192 + Icges(7,4) * t191 + Icges(7,5) * t228;
t144 = Icges(7,4) * t194 + Icges(7,2) * t193 + Icges(7,6) * t230;
t143 = Icges(7,4) * t192 + Icges(7,2) * t191 + Icges(7,6) * t228;
t142 = Icges(7,5) * t194 + Icges(7,6) * t193 + Icges(7,3) * t230;
t141 = Icges(7,5) * t192 + Icges(7,6) * t191 + Icges(7,3) * t228;
t140 = t218 * t272 + (-t184 - t223) * t284 + t310;
t139 = t185 * t284 + (-t218 - t264) * t273 + t309;
t138 = t273 * t184 + (-t185 - t224) * t272 + t308;
t137 = -t173 * t261 + t205 * t239 + t307;
t136 = t174 * t261 - t205 * t240 + t306;
t135 = t240 * t173 - t239 * t174 + t304;
t134 = t204 * t239 + (-t171 - t187) * t261 + t305;
t133 = t172 * t261 + (-t204 - t211) * t240 + t303;
t132 = t240 * t171 + (-t172 - t188) * t239 + t302;
t131 = -t147 * t222 + t158 * t189 + t238 * t239 + (-t187 - t196) * t261 + t305;
t130 = t148 * t222 - t158 * t190 + t197 * t261 + (-t211 - t238) * t240 + t303;
t129 = -t189 * t148 + t190 * t147 + t302 + (-t188 - t197) * t239 + t240 * t196;
t1 = m(3) * (t149 ^ 2 + t150 ^ 2 + t151 ^ 2) / 0.2e1 + m(5) * (t135 ^ 2 + t136 ^ 2 + t137 ^ 2) / 0.2e1 + m(4) * (t138 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(6) * (t132 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(7) * (t129 ^ 2 + t130 ^ 2 + t131 ^ 2) / 0.2e1 + t189 * ((t142 * t228 + t144 * t191 + t146 * t192) * t190 + (t228 * t141 + t191 * t143 + t192 * t145) * t189 + (t155 * t228 + t156 * t191 + t157 * t192) * t222) / 0.2e1 + t190 * ((t230 * t142 + t193 * t144 + t194 * t146) * t190 + (t141 * t230 + t143 * t193 + t145 * t194) * t189 + (t155 * t230 + t156 * t193 + t157 * t194) * t222) / 0.2e1 + m(2) * (t233 ^ 2 + t242 ^ 2 + t243 ^ 2) / 0.2e1 + t222 * ((t142 * t249 + t144 * t231 + t146 * t232) * t190 + (t141 * t249 + t143 * t231 + t145 * t232) * t189 + (t155 * t249 + t156 * t231 + t157 * t232) * t222) / 0.2e1 + m(1) * (t269 ^ 2 + t270 ^ 2 + t271 ^ 2) / 0.2e1 + ((t227 * t348 + t228 * t347 + t257 * t346) * t261 + (t227 * t355 + t228 * t351 + t257 * t353) * t240 + (t356 * t227 + t352 * t228 + t354 * t257) * t239) * t239 / 0.2e1 + ((t229 * t348 + t230 * t347 + t259 * t346) * t261 + (t355 * t229 + t351 * t230 + t353 * t259) * t240 + (t229 * t356 + t230 * t352 + t259 * t354) * t239) * t240 / 0.2e1 + ((t348 * t248 + t347 * t249 - t346 * t330) * t261 + (t248 * t355 + t249 * t351 - t330 * t353) * t240 + (t248 * t356 + t249 * t352 - t330 * t354) * t239) * t261 / 0.2e1 + ((t209 * t234 + t210 * t235 - t245 * t332 + t247 * t258 + t257 * t345) * t284 + (t181 * t234 + t183 * t235 - t213 * t332 + t217 * t258 + t257 * t349) * t273 + (t180 * t234 + t182 * t235 - t212 * t332 + t216 * t258 + t350 * t257) * t272) * t272 / 0.2e1 + ((t209 * t236 + t210 * t237 + t245 * t333 + t247 * t260 + t259 * t345) * t284 + (t181 * t236 + t183 * t237 + t213 * t333 + t217 * t260 + t349 * t259) * t273 + (t180 * t236 + t182 * t237 + t212 * t333 + t216 * t260 + t259 * t350) * t272) * t273 / 0.2e1 + ((t212 * t272 + t213 * t273 + t245 * t284) * t296 + ((t215 * t301 + t217 * t299) * t273 + (t214 * t301 + t216 * t299) * t272 + (t246 * t301 + t247 * t299) * t284) * t293 + (-t179 * t330 + t255 * t181 + t256 * t183) * t273 + (-t178 * t330 + t255 * t180 + t256 * t182) * t272 + (-t208 * t330 + t255 * t209 + t256 * t210) * t284) * t284 / 0.2e1 + ((-t276 * t292 + t278 * t295 + Icges(1,4)) * V_base(5) + (-t277 * t292 + t279 * t295 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t276 * t295 + t278 * t292 + Icges(1,2)) * V_base(5) + (t277 * t295 + t279 * t292 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(2,5) * t295 - Icges(2,6) * t292 + Icges(1,5)) * V_base(4) + (Icges(2,5) * t292 + Icges(2,6) * t295 + Icges(1,6)) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(1,3) / 0.2e1) * V_base(6)) * V_base(6);
T  = t1;
