% Calculate kinetic energy for
% S6RPRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRP6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRP6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:14:09
% EndTime: 2019-03-09 06:14:13
% DurationCPUTime: 3.65s
% Computational Cost: add. (2275->341), mult. (2327->494), div. (0->0), fcn. (2221->10), ass. (0->173)
t344 = Icges(6,1) + Icges(7,1);
t343 = Icges(6,4) + Icges(7,4);
t342 = -Icges(7,5) - Icges(6,5);
t341 = Icges(6,2) + Icges(7,2);
t340 = -Icges(7,6) - Icges(6,6);
t339 = -Icges(7,3) - Icges(6,3);
t255 = pkin(10) + qJ(3);
t245 = cos(t255);
t256 = qJ(4) + qJ(5);
t250 = cos(t256);
t263 = cos(qJ(1));
t310 = t250 * t263;
t249 = sin(t256);
t261 = sin(qJ(1));
t313 = t249 * t261;
t189 = -t245 * t313 - t310;
t311 = t250 * t261;
t312 = t249 * t263;
t190 = t245 * t311 - t312;
t244 = sin(t255);
t315 = t244 * t261;
t338 = -t340 * t189 - t342 * t190 - t339 * t315;
t191 = -t245 * t312 + t311;
t192 = t245 * t310 + t313;
t314 = t244 * t263;
t337 = -t340 * t191 - t342 * t192 - t339 * t314;
t336 = t341 * t189 + t343 * t190 - t340 * t315;
t335 = t341 * t191 + t343 * t192 - t340 * t314;
t334 = t343 * t189 + t344 * t190 - t342 * t315;
t333 = t343 * t191 + t344 * t192 - t342 * t314;
t332 = t339 * t245 + (t340 * t249 - t342 * t250) * t244;
t331 = t340 * t245 + (-t341 * t249 + t343 * t250) * t244;
t330 = t342 * t245 + (-t343 * t249 + t344 * t250) * t244;
t257 = sin(pkin(10));
t325 = pkin(2) * t257;
t258 = cos(pkin(10));
t323 = pkin(2) * t258;
t262 = cos(qJ(4));
t322 = t262 * pkin(4);
t320 = Icges(2,4) * t261;
t319 = Icges(3,4) * t257;
t318 = Icges(3,4) * t258;
t317 = Icges(4,4) * t244;
t316 = Icges(4,4) * t245;
t260 = sin(qJ(4));
t309 = t260 * t261;
t308 = t260 * t263;
t307 = t261 * t262;
t306 = t262 * t263;
t300 = pkin(5) * t250;
t276 = qJ(6) * t244 + t245 * t300;
t291 = pkin(5) * t249;
t304 = rSges(7,1) * t190 + rSges(7,2) * t189 + rSges(7,3) * t315 + t261 * t276 - t263 * t291;
t303 = rSges(7,1) * t192 + rSges(7,2) * t191 + rSges(7,3) * t314 + t261 * t291 + t263 * t276;
t302 = (-qJ(6) - rSges(7,3)) * t245 + (rSges(7,1) * t250 - rSges(7,2) * t249 + t300) * t244;
t177 = -pkin(7) * t263 + t261 * t323;
t235 = pkin(1) * t261 - qJ(2) * t263;
t301 = -t177 - t235;
t298 = qJD(4) * t244;
t297 = qJD(5) * t244;
t296 = qJD(6) * t244;
t295 = V_base(4) * t235 + V_base(3);
t294 = V_base(5) * pkin(6) + V_base(1);
t240 = qJD(3) * t261 + V_base(4);
t246 = V_base(6) + qJD(1);
t290 = qJD(2) * t261 + t294;
t205 = t263 * t298 + t240;
t289 = V_base(5) * t325 + t290;
t288 = pkin(3) * t245 + pkin(8) * t244;
t239 = -qJD(3) * t263 + V_base(5);
t287 = rSges(3,1) * t258 - rSges(3,2) * t257;
t286 = rSges(4,1) * t245 - rSges(4,2) * t244;
t285 = Icges(3,1) * t258 - t319;
t284 = Icges(4,1) * t245 - t317;
t283 = -Icges(3,2) * t257 + t318;
t282 = -Icges(4,2) * t244 + t316;
t281 = Icges(3,5) * t258 - Icges(3,6) * t257;
t280 = Icges(4,5) * t245 - Icges(4,6) * t244;
t237 = pkin(1) * t263 + qJ(2) * t261;
t279 = -qJD(2) * t263 + t246 * t237 + V_base(2);
t204 = t261 * t298 + t239;
t278 = pkin(9) * t244 + t245 * t322;
t277 = (-Icges(4,3) * t263 + t261 * t280) * t239 + (Icges(4,3) * t261 + t263 * t280) * t240 + (Icges(4,5) * t244 + Icges(4,6) * t245) * t246;
t178 = pkin(7) * t261 + t263 * t323;
t275 = V_base(4) * t177 + (-t178 - t237) * V_base(5) + t295;
t274 = (-Icges(3,3) * t263 + t261 * t281) * V_base(5) + (Icges(3,3) * t261 + t263 * t281) * V_base(4) + (Icges(3,5) * t257 + Icges(3,6) * t258) * t246;
t201 = t288 * t261;
t215 = t244 * pkin(3) - t245 * pkin(8);
t273 = t239 * t215 + (-t201 + t301) * t246 + t289;
t202 = t288 * t263;
t272 = t240 * t201 - t202 * t239 + t275;
t271 = t246 * t178 + (-pkin(6) - t325) * V_base(4) + t279;
t145 = -pkin(4) * t308 + t261 * t278;
t156 = -pkin(9) * t245 + t244 * t322;
t218 = -qJD(4) * t245 + t246;
t270 = -t145 * t218 + t204 * t156 + t273;
t146 = pkin(4) * t309 + t263 * t278;
t269 = t205 * t145 - t146 * t204 + t272;
t268 = t246 * t202 - t215 * t240 + t271;
t267 = t218 * t146 - t156 * t205 + t268;
t182 = -Icges(4,6) * t263 + t261 * t282;
t183 = Icges(4,6) * t261 + t263 * t282;
t184 = -Icges(4,5) * t263 + t261 * t284;
t185 = Icges(4,5) * t261 + t263 * t284;
t212 = Icges(4,2) * t245 + t317;
t213 = Icges(4,1) * t244 + t316;
t266 = (-t183 * t244 + t185 * t245) * t240 + (-t182 * t244 + t184 * t245) * t239 + (-t212 * t244 + t213 * t245) * t246;
t195 = -Icges(3,6) * t263 + t261 * t283;
t196 = Icges(3,6) * t261 + t263 * t283;
t197 = -Icges(3,5) * t263 + t261 * t285;
t198 = Icges(3,5) * t261 + t263 * t285;
t224 = Icges(3,2) * t258 + t319;
t225 = Icges(3,1) * t257 + t318;
t265 = (-t196 * t257 + t198 * t258) * V_base(4) + (-t195 * t257 + t197 * t258) * V_base(5) + (-t224 * t257 + t225 * t258) * t246;
t251 = Icges(2,4) * t263;
t238 = rSges(2,1) * t263 - rSges(2,2) * t261;
t236 = rSges(2,1) * t261 + rSges(2,2) * t263;
t232 = Icges(2,1) * t263 - t320;
t231 = Icges(2,1) * t261 + t251;
t230 = -Icges(2,2) * t261 + t251;
t229 = Icges(2,2) * t263 + t320;
t228 = Icges(2,5) * t263 - Icges(2,6) * t261;
t227 = Icges(2,5) * t261 + Icges(2,6) * t263;
t226 = rSges(3,1) * t257 + rSges(3,2) * t258;
t221 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t220 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t219 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t214 = rSges(4,1) * t244 + rSges(4,2) * t245;
t209 = t245 * t306 + t309;
t208 = -t245 * t308 + t307;
t207 = t245 * t307 - t308;
t206 = -t245 * t309 - t306;
t203 = (-qJD(4) - qJD(5)) * t245 + t246;
t200 = rSges(3,3) * t261 + t263 * t287;
t199 = -rSges(3,3) * t263 + t261 * t287;
t187 = rSges(4,3) * t261 + t263 * t286;
t186 = -rSges(4,3) * t263 + t261 * t286;
t176 = -rSges(5,3) * t245 + (rSges(5,1) * t262 - rSges(5,2) * t260) * t244;
t175 = V_base(5) * rSges(2,3) - t236 * t246 + t294;
t174 = t238 * t246 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t173 = -Icges(5,5) * t245 + (Icges(5,1) * t262 - Icges(5,4) * t260) * t244;
t172 = -Icges(5,6) * t245 + (Icges(5,4) * t262 - Icges(5,2) * t260) * t244;
t171 = -Icges(5,3) * t245 + (Icges(5,5) * t262 - Icges(5,6) * t260) * t244;
t170 = t263 * t297 + t205;
t169 = t261 * t297 + t204;
t167 = t236 * V_base(4) - t238 * V_base(5) + V_base(3);
t165 = -rSges(6,3) * t245 + (rSges(6,1) * t250 - rSges(6,2) * t249) * t244;
t154 = rSges(5,1) * t209 + rSges(5,2) * t208 + rSges(5,3) * t314;
t153 = rSges(5,1) * t207 + rSges(5,2) * t206 + rSges(5,3) * t315;
t152 = Icges(5,1) * t209 + Icges(5,4) * t208 + Icges(5,5) * t314;
t151 = Icges(5,1) * t207 + Icges(5,4) * t206 + Icges(5,5) * t315;
t150 = Icges(5,4) * t209 + Icges(5,2) * t208 + Icges(5,6) * t314;
t149 = Icges(5,4) * t207 + Icges(5,2) * t206 + Icges(5,6) * t315;
t148 = Icges(5,5) * t209 + Icges(5,6) * t208 + Icges(5,3) * t314;
t147 = Icges(5,5) * t207 + Icges(5,6) * t206 + Icges(5,3) * t315;
t143 = rSges(6,1) * t192 + rSges(6,2) * t191 + rSges(6,3) * t314;
t141 = rSges(6,1) * t190 + rSges(6,2) * t189 + rSges(6,3) * t315;
t126 = t226 * V_base(5) + (-t199 - t235) * t246 + t290;
t125 = t200 * t246 + (-pkin(6) - t226) * V_base(4) + t279;
t123 = t199 * V_base(4) + (-t200 - t237) * V_base(5) + t295;
t120 = t214 * t239 + (-t186 + t301) * t246 + t289;
t119 = t187 * t246 - t214 * t240 + t271;
t118 = t186 * t240 - t187 * t239 + t275;
t117 = -t153 * t218 + t176 * t204 + t273;
t116 = t154 * t218 - t176 * t205 + t268;
t115 = t153 * t205 - t154 * t204 + t272;
t114 = -t141 * t203 + t165 * t169 + t270;
t113 = t143 * t203 - t165 * t170 + t267;
t112 = t141 * t170 - t143 * t169 + t269;
t111 = t169 * t302 - t203 * t304 + t263 * t296 + t270;
t110 = -t170 * t302 + t203 * t303 + t261 * t296 + t267;
t109 = -qJD(6) * t245 - t169 * t303 + t170 * t304 + t269;
t1 = m(2) * (t167 ^ 2 + t174 ^ 2 + t175 ^ 2) / 0.2e1 + m(1) * (t219 ^ 2 + t220 ^ 2 + t221 ^ 2) / 0.2e1 + m(3) * (t123 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(4) * (t118 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(7) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(6) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(5) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + t240 * (t277 * t261 + t266 * t263) / 0.2e1 + t239 * (t266 * t261 - t277 * t263) / 0.2e1 + t205 * ((t148 * t314 + t150 * t208 + t152 * t209) * t205 + (t147 * t314 + t149 * t208 + t151 * t209) * t204 + (t171 * t314 + t172 * t208 + t173 * t209) * t218) / 0.2e1 + t204 * ((t148 * t315 + t150 * t206 + t152 * t207) * t205 + (t147 * t315 + t149 * t206 + t151 * t207) * t204 + (t171 * t315 + t172 * t206 + t173 * t207) * t218) / 0.2e1 + t218 * ((-t147 * t204 - t148 * t205 - t171 * t218) * t245 + ((-t150 * t260 + t152 * t262) * t205 + (-t149 * t260 + t151 * t262) * t204 + (-t172 * t260 + t173 * t262) * t218) * t244) / 0.2e1 + ((t189 * t331 + t190 * t330 + t315 * t332) * t203 + (t189 * t335 + t190 * t333 + t315 * t337) * t170 + (t336 * t189 + t334 * t190 + t338 * t315) * t169) * t169 / 0.2e1 + ((t191 * t331 + t192 * t330 + t314 * t332) * t203 + (t335 * t191 + t333 * t192 + t337 * t314) * t170 + (t336 * t191 + t334 * t192 + t314 * t338) * t169) * t170 / 0.2e1 + ((-t169 * t338 - t337 * t170 - t332 * t203) * t245 + ((-t249 * t331 + t250 * t330) * t203 + (-t249 * t335 + t250 * t333) * t170 + (-t249 * t336 + t250 * t334) * t169) * t244) * t203 / 0.2e1 + ((t183 * t245 + t185 * t244) * t240 + (t182 * t245 + t184 * t244) * t239 + (t195 * t258 + t197 * t257 + t227) * V_base(5) + (t196 * t258 + t198 * t257 + t228) * V_base(4) + (t212 * t245 + t213 * t244 + t224 * t258 + t225 * t257 + Icges(2,3)) * t246) * t246 / 0.2e1 + (t228 * t246 + t274 * t261 + t265 * t263 + (-t229 * t261 + t231 * t263 + Icges(1,4)) * V_base(5) + (-t230 * t261 + t232 * t263 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t227 * t246 + t265 * t261 - t274 * t263 + (t229 * t263 + t231 * t261 + Icges(1,2)) * V_base(5) + (t230 * t263 + t232 * t261 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
