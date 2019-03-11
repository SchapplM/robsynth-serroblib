% Calculate kinetic energy for
% S6RRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRRP1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRRP1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:55:29
% EndTime: 2019-03-10 00:55:32
% DurationCPUTime: 2.91s
% Computational Cost: add. (2239->314), mult. (2089->461), div. (0->0), fcn. (1927->10), ass. (0->167)
t335 = Icges(6,1) + Icges(7,1);
t334 = Icges(6,4) + Icges(7,4);
t333 = -Icges(7,5) - Icges(6,5);
t332 = Icges(6,2) + Icges(7,2);
t331 = -Icges(7,6) - Icges(6,6);
t330 = -Icges(7,3) - Icges(6,3);
t245 = qJ(2) + qJ(3);
t241 = qJ(4) + t245;
t231 = cos(t241);
t250 = cos(qJ(5));
t252 = cos(qJ(1));
t298 = t250 * t252;
t247 = sin(qJ(5));
t249 = sin(qJ(1));
t301 = t247 * t249;
t189 = -t231 * t301 - t298;
t299 = t249 * t250;
t300 = t247 * t252;
t190 = t231 * t299 - t300;
t230 = sin(t241);
t303 = t230 * t249;
t329 = -t331 * t189 - t333 * t190 - t330 * t303;
t191 = -t231 * t300 + t299;
t192 = t231 * t298 + t301;
t302 = t230 * t252;
t328 = -t331 * t191 - t333 * t192 - t330 * t302;
t327 = t332 * t189 + t334 * t190 - t331 * t303;
t326 = t332 * t191 + t334 * t192 - t331 * t302;
t325 = t334 * t189 + t335 * t190 - t333 * t303;
t324 = t334 * t191 + t335 * t192 - t333 * t302;
t323 = t330 * t231 + (t331 * t247 - t333 * t250) * t230;
t322 = t331 * t231 + (-t332 * t247 + t334 * t250) * t230;
t321 = t333 * t231 + (-t334 * t247 + t335 * t250) * t230;
t248 = sin(qJ(2));
t316 = pkin(2) * t248;
t238 = sin(t245);
t315 = pkin(3) * t238;
t251 = cos(qJ(2));
t314 = t251 * pkin(2);
t313 = pkin(5) * t250;
t310 = Icges(2,4) * t249;
t309 = Icges(3,4) * t248;
t308 = Icges(3,4) * t251;
t307 = Icges(4,4) * t238;
t239 = cos(t245);
t306 = Icges(4,4) * t239;
t305 = Icges(5,4) * t230;
t304 = Icges(5,4) * t231;
t267 = qJ(6) * t230 + t231 * t313;
t297 = rSges(7,1) * t190 + rSges(7,2) * t189 + rSges(7,3) * t303 - pkin(5) * t300 + t249 * t267;
t296 = rSges(7,1) * t192 + rSges(7,2) * t191 + rSges(7,3) * t302 + pkin(5) * t301 + t252 * t267;
t295 = (-qJ(6) - rSges(7,3)) * t231 + (rSges(7,1) * t250 - rSges(7,2) * t247 + t313) * t230;
t168 = -pkin(8) * t252 + t249 * t314;
t226 = t249 * pkin(1) - t252 * pkin(7);
t294 = -t168 - t226;
t293 = pkin(3) * t239;
t291 = qJD(5) * t230;
t290 = qJD(6) * t230;
t289 = -qJD(2) - qJD(3);
t288 = V_base(5) * pkin(6) + V_base(1);
t140 = -pkin(9) * t252 + t249 * t293;
t285 = -t140 + t294;
t229 = qJD(2) * t249 + V_base(4);
t234 = V_base(6) + qJD(1);
t228 = -qJD(2) * t252 + V_base(5);
t284 = t228 * t316 + t288;
t206 = qJD(3) * t249 + t229;
t205 = t252 * t289 + V_base(5);
t283 = t205 * t315 + t284;
t282 = pkin(4) * t231 + pkin(10) * t230;
t281 = rSges(3,1) * t251 - rSges(3,2) * t248;
t280 = rSges(4,1) * t239 - rSges(4,2) * t238;
t279 = rSges(5,1) * t231 - rSges(5,2) * t230;
t194 = qJD(4) * t249 + t206;
t278 = Icges(3,1) * t251 - t309;
t277 = Icges(4,1) * t239 - t307;
t276 = Icges(5,1) * t231 - t305;
t275 = -Icges(3,2) * t248 + t308;
t274 = -Icges(4,2) * t238 + t306;
t273 = -Icges(5,2) * t230 + t304;
t272 = Icges(3,5) * t251 - Icges(3,6) * t248;
t271 = Icges(4,5) * t239 - Icges(4,6) * t238;
t270 = Icges(5,5) * t231 - Icges(5,6) * t230;
t227 = t252 * pkin(1) + t249 * pkin(7);
t269 = -V_base(4) * pkin(6) + t234 * t227 + V_base(2);
t268 = V_base(4) * t226 - t227 * V_base(5) + V_base(3);
t193 = V_base(5) + (-qJD(4) + t289) * t252;
t266 = (-Icges(5,3) * t252 + t249 * t270) * t193 + (Icges(5,3) * t249 + t252 * t270) * t194 + (Icges(5,5) * t230 + Icges(5,6) * t231) * t234;
t265 = (-Icges(4,3) * t252 + t249 * t271) * t205 + (Icges(4,3) * t249 + t252 * t271) * t206 + (Icges(4,5) * t238 + Icges(4,6) * t239) * t234;
t264 = (-Icges(3,3) * t252 + t249 * t272) * t228 + (Icges(3,3) * t249 + t252 * t272) * t229 + (Icges(3,5) * t248 + Icges(3,6) * t251) * t234;
t169 = pkin(8) * t249 + t252 * t314;
t263 = t229 * t168 - t169 * t228 + t268;
t262 = t234 * t169 - t229 * t316 + t269;
t178 = t282 * t249;
t199 = pkin(4) * t230 - pkin(10) * t231;
t261 = t193 * t199 + (-t178 + t285) * t234 + t283;
t141 = pkin(9) * t249 + t252 * t293;
t260 = t206 * t140 - t141 * t205 + t263;
t259 = t234 * t141 - t206 * t315 + t262;
t179 = t282 * t252;
t258 = t194 * t178 - t179 * t193 + t260;
t257 = t234 * t179 - t194 * t199 + t259;
t162 = -Icges(5,6) * t252 + t249 * t273;
t163 = Icges(5,6) * t249 + t252 * t273;
t164 = -Icges(5,5) * t252 + t249 * t276;
t165 = Icges(5,5) * t249 + t252 * t276;
t196 = Icges(5,2) * t231 + t305;
t197 = Icges(5,1) * t230 + t304;
t256 = (-t163 * t230 + t165 * t231) * t194 + (-t162 * t230 + t164 * t231) * t193 + (-t196 * t230 + t197 * t231) * t234;
t172 = -Icges(4,6) * t252 + t249 * t274;
t173 = Icges(4,6) * t249 + t252 * t274;
t174 = -Icges(4,5) * t252 + t249 * t277;
t175 = Icges(4,5) * t249 + t252 * t277;
t202 = Icges(4,2) * t239 + t307;
t203 = Icges(4,1) * t238 + t306;
t255 = (-t173 * t238 + t175 * t239) * t206 + (-t172 * t238 + t174 * t239) * t205 + (-t202 * t238 + t203 * t239) * t234;
t182 = -Icges(3,6) * t252 + t249 * t275;
t183 = Icges(3,6) * t249 + t252 * t275;
t184 = -Icges(3,5) * t252 + t249 * t278;
t185 = Icges(3,5) * t249 + t252 * t278;
t217 = Icges(3,2) * t251 + t309;
t220 = Icges(3,1) * t248 + t308;
t254 = (-t183 * t248 + t185 * t251) * t229 + (-t182 * t248 + t184 * t251) * t228 + (-t217 * t248 + t220 * t251) * t234;
t240 = Icges(2,4) * t252;
t225 = rSges(2,1) * t252 - rSges(2,2) * t249;
t224 = rSges(2,1) * t249 + rSges(2,2) * t252;
t223 = rSges(3,1) * t248 + rSges(3,2) * t251;
t222 = Icges(2,1) * t252 - t310;
t221 = Icges(2,1) * t249 + t240;
t219 = -Icges(2,2) * t249 + t240;
t218 = Icges(2,2) * t252 + t310;
t213 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t212 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t211 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t207 = -qJD(5) * t231 + t234;
t204 = rSges(4,1) * t238 + rSges(4,2) * t239;
t198 = rSges(5,1) * t230 + rSges(5,2) * t231;
t187 = rSges(3,3) * t249 + t252 * t281;
t186 = -rSges(3,3) * t252 + t249 * t281;
t177 = rSges(4,3) * t249 + t252 * t280;
t176 = -rSges(4,3) * t252 + t249 * t280;
t167 = rSges(5,3) * t249 + t252 * t279;
t166 = -rSges(5,3) * t252 + t249 * t279;
t158 = V_base(5) * rSges(2,3) - t224 * t234 + t288;
t157 = t225 * t234 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t156 = t252 * t291 + t194;
t155 = t249 * t291 + t193;
t154 = t224 * V_base(4) - t225 * V_base(5) + V_base(3);
t153 = -rSges(6,3) * t231 + (rSges(6,1) * t250 - rSges(6,2) * t247) * t230;
t136 = rSges(6,1) * t192 + rSges(6,2) * t191 + rSges(6,3) * t302;
t134 = rSges(6,1) * t190 + rSges(6,2) * t189 + rSges(6,3) * t303;
t118 = t223 * t228 + (-t186 - t226) * t234 + t288;
t117 = t187 * t234 - t223 * t229 + t269;
t116 = t186 * t229 - t187 * t228 + t268;
t115 = t204 * t205 + (-t176 + t294) * t234 + t284;
t114 = t177 * t234 - t204 * t206 + t262;
t113 = t176 * t206 - t177 * t205 + t263;
t112 = t193 * t198 + (-t166 + t285) * t234 + t283;
t111 = t167 * t234 - t194 * t198 + t259;
t110 = t166 * t194 - t167 * t193 + t260;
t109 = -t134 * t207 + t153 * t155 + t261;
t108 = t136 * t207 - t153 * t156 + t257;
t107 = t134 * t156 - t136 * t155 + t258;
t106 = t155 * t295 - t207 * t297 + t252 * t290 + t261;
t105 = -t156 * t295 + t207 * t296 + t249 * t290 + t257;
t104 = -qJD(6) * t231 - t155 * t296 + t156 * t297 + t258;
t1 = m(1) * (t211 ^ 2 + t212 ^ 2 + t213 ^ 2) / 0.2e1 + m(2) * (t154 ^ 2 + t157 ^ 2 + t158 ^ 2) / 0.2e1 + m(4) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(3) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(6) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(5) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(7) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + t229 * (t264 * t249 + t254 * t252) / 0.2e1 + t228 * (t254 * t249 - t264 * t252) / 0.2e1 + t206 * (t265 * t249 + t255 * t252) / 0.2e1 + t205 * (t255 * t249 - t265 * t252) / 0.2e1 + t194 * (t266 * t249 + t256 * t252) / 0.2e1 + t193 * (t256 * t249 - t266 * t252) / 0.2e1 + ((t189 * t322 + t190 * t321 + t303 * t323) * t207 + (t189 * t326 + t190 * t324 + t303 * t328) * t156 + (t327 * t189 + t325 * t190 + t329 * t303) * t155) * t155 / 0.2e1 + ((t191 * t322 + t192 * t321 + t302 * t323) * t207 + (t326 * t191 + t324 * t192 + t328 * t302) * t156 + (t327 * t191 + t325 * t192 + t302 * t329) * t155) * t156 / 0.2e1 + ((-t155 * t329 - t328 * t156 - t323 * t207) * t231 + ((-t247 * t322 + t250 * t321) * t207 + (-t247 * t326 + t250 * t324) * t156 + (-t247 * t327 + t250 * t325) * t155) * t230) * t207 / 0.2e1 + ((-t218 * t249 + t221 * t252 + Icges(1,4)) * V_base(5) + (-t219 * t249 + t222 * t252 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t218 * t252 + t221 * t249 + Icges(1,2)) * V_base(5) + (t219 * t252 + t222 * t249 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t183 * t251 + t185 * t248) * t229 + (t182 * t251 + t184 * t248) * t228 + (t173 * t239 + t175 * t238) * t206 + (t172 * t239 + t174 * t238) * t205 + (t163 * t231 + t165 * t230) * t194 + (t162 * t231 + t164 * t230) * t193 + (t196 * t231 + t197 * t230 + t202 * t239 + t203 * t238 + t217 * t251 + t220 * t248 + Icges(2,3)) * t234) * t234 / 0.2e1 + V_base(4) * t234 * (Icges(2,5) * t252 - Icges(2,6) * t249) + V_base(5) * t234 * (Icges(2,5) * t249 + Icges(2,6) * t252) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
