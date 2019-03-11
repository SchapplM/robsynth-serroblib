% Calculate kinetic energy for
% S6RPRRRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR5_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR5_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR5_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR5_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR5_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR5_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:08:07
% EndTime: 2019-03-09 07:08:10
% DurationCPUTime: 3.02s
% Computational Cost: add. (2334->365), mult. (2024->543), div. (0->0), fcn. (1862->12), ass. (0->187)
t262 = sin(pkin(11));
t333 = pkin(2) * t262;
t260 = pkin(11) + qJ(3);
t247 = sin(t260);
t332 = pkin(3) * t247;
t263 = cos(pkin(11));
t331 = t263 * pkin(2);
t267 = cos(qJ(5));
t330 = pkin(5) * t267;
t266 = sin(qJ(1));
t328 = Icges(2,4) * t266;
t327 = Icges(3,4) * t262;
t326 = Icges(3,4) * t263;
t325 = Icges(4,4) * t247;
t248 = cos(t260);
t324 = Icges(4,4) * t248;
t249 = qJ(4) + t260;
t243 = sin(t249);
t323 = Icges(5,4) * t243;
t244 = cos(t249);
t322 = Icges(5,4) * t244;
t321 = t243 * t266;
t268 = cos(qJ(1));
t320 = t243 * t268;
t261 = qJ(5) + qJ(6);
t254 = sin(t261);
t319 = t254 * t266;
t318 = t254 * t268;
t255 = cos(t261);
t317 = t255 * t266;
t316 = t255 * t268;
t265 = sin(qJ(5));
t315 = t265 * t266;
t314 = t265 * t268;
t313 = t266 * t267;
t312 = t267 * t268;
t174 = -pkin(7) * t268 + t266 * t331;
t236 = pkin(1) * t266 - qJ(2) * t268;
t310 = -t174 - t236;
t309 = pkin(3) * t248;
t307 = qJD(5) * t243;
t306 = qJD(6) * t243;
t305 = V_base(4) * t236 + V_base(3);
t304 = V_base(5) * pkin(6) + V_base(1);
t145 = -pkin(8) * t268 + t266 * t309;
t301 = -t145 + t310;
t241 = qJD(3) * t266 + V_base(4);
t250 = V_base(6) + qJD(1);
t300 = qJD(2) * t266 + t304;
t218 = qJD(4) * t266 + t241;
t299 = V_base(5) * t333 + t300;
t298 = pkin(4) * t244 + pkin(9) * t243;
t297 = rSges(3,1) * t263 - rSges(3,2) * t262;
t296 = rSges(4,1) * t248 - rSges(4,2) * t247;
t295 = rSges(5,1) * t244 - rSges(5,2) * t243;
t183 = t268 * t307 + t218;
t294 = Icges(3,1) * t263 - t327;
t293 = Icges(4,1) * t248 - t325;
t292 = Icges(5,1) * t244 - t323;
t291 = -Icges(3,2) * t262 + t326;
t290 = -Icges(4,2) * t247 + t324;
t289 = -Icges(5,2) * t243 + t322;
t288 = Icges(3,5) * t263 - Icges(3,6) * t262;
t287 = Icges(4,5) * t248 - Icges(4,6) * t247;
t286 = Icges(5,5) * t244 - Icges(5,6) * t243;
t238 = pkin(1) * t268 + qJ(2) * t266;
t285 = -qJD(2) * t268 + t250 * t238 + V_base(2);
t240 = -qJD(3) * t268 + V_base(5);
t284 = t240 * t332 + t299;
t217 = V_base(5) + (-qJD(3) - qJD(4)) * t268;
t182 = t266 * t307 + t217;
t283 = pkin(10) * t243 + t244 * t330;
t282 = (-Icges(5,3) * t268 + t266 * t286) * t217 + (Icges(5,3) * t266 + t268 * t286) * t218 + (Icges(5,5) * t243 + Icges(5,6) * t244) * t250;
t281 = (-Icges(4,3) * t268 + t266 * t287) * t240 + (Icges(4,3) * t266 + t268 * t287) * t241 + (Icges(4,5) * t247 + Icges(4,6) * t248) * t250;
t175 = pkin(7) * t266 + t268 * t331;
t280 = V_base(4) * t174 + (-t175 - t238) * V_base(5) + t305;
t279 = (-Icges(3,3) * t268 + t266 * t288) * V_base(5) + (Icges(3,3) * t266 + t268 * t288) * V_base(4) + (Icges(3,5) * t262 + Icges(3,6) * t263) * t250;
t146 = pkin(8) * t266 + t268 * t309;
t278 = t241 * t145 - t146 * t240 + t280;
t277 = t250 * t175 + (-pkin(6) - t333) * V_base(4) + t285;
t190 = t298 * t266;
t209 = t243 * pkin(4) - t244 * pkin(9);
t276 = t217 * t209 + (-t190 + t301) * t250 + t284;
t191 = t298 * t268;
t275 = t218 * t190 - t191 * t217 + t278;
t274 = t250 * t146 - t241 * t332 + t277;
t273 = t250 * t191 - t209 * t218 + t274;
t166 = -Icges(5,6) * t268 + t266 * t289;
t167 = Icges(5,6) * t266 + t268 * t289;
t168 = -Icges(5,5) * t268 + t266 * t292;
t169 = Icges(5,5) * t266 + t268 * t292;
t206 = Icges(5,2) * t244 + t323;
t207 = Icges(5,1) * t243 + t322;
t272 = (-t167 * t243 + t169 * t244) * t218 + (-t166 * t243 + t168 * t244) * t217 + (-t206 * t243 + t207 * t244) * t250;
t178 = -Icges(4,6) * t268 + t266 * t290;
t179 = Icges(4,6) * t266 + t268 * t290;
t180 = -Icges(4,5) * t268 + t266 * t293;
t181 = Icges(4,5) * t266 + t268 * t293;
t213 = Icges(4,2) * t248 + t325;
t214 = Icges(4,1) * t247 + t324;
t271 = (-t179 * t247 + t181 * t248) * t241 + (-t178 * t247 + t180 * t248) * t240 + (-t213 * t247 + t214 * t248) * t250;
t195 = -Icges(3,6) * t268 + t266 * t291;
t196 = Icges(3,6) * t266 + t268 * t291;
t197 = -Icges(3,5) * t268 + t266 * t294;
t198 = Icges(3,5) * t266 + t268 * t294;
t227 = Icges(3,2) * t263 + t327;
t228 = Icges(3,1) * t262 + t326;
t270 = (-t196 * t262 + t198 * t263) * V_base(4) + (-t195 * t262 + t197 * t263) * V_base(5) + (-t227 * t262 + t228 * t263) * t250;
t257 = Icges(2,4) * t268;
t239 = rSges(2,1) * t268 - rSges(2,2) * t266;
t237 = rSges(2,1) * t266 + rSges(2,2) * t268;
t235 = Icges(2,1) * t268 - t328;
t234 = Icges(2,1) * t266 + t257;
t233 = -Icges(2,2) * t266 + t257;
t232 = Icges(2,2) * t268 + t328;
t231 = Icges(2,5) * t268 - Icges(2,6) * t266;
t230 = Icges(2,5) * t266 + Icges(2,6) * t268;
t229 = rSges(3,1) * t262 + rSges(3,2) * t263;
t223 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t222 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t221 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t216 = -qJD(5) * t244 + t250;
t215 = rSges(4,1) * t247 + rSges(4,2) * t248;
t208 = rSges(5,1) * t243 + rSges(5,2) * t244;
t204 = t244 * t312 + t315;
t203 = -t244 * t314 + t313;
t202 = t244 * t313 - t314;
t201 = -t244 * t315 - t312;
t200 = rSges(3,3) * t266 + t268 * t297;
t199 = -rSges(3,3) * t268 + t266 * t297;
t192 = (-qJD(5) - qJD(6)) * t244 + t250;
t189 = t244 * t316 + t319;
t188 = -t244 * t318 + t317;
t187 = t244 * t317 - t318;
t186 = -t244 * t319 - t316;
t185 = rSges(4,3) * t266 + t268 * t296;
t184 = -rSges(4,3) * t268 + t266 * t296;
t173 = rSges(5,3) * t266 + t268 * t295;
t172 = -rSges(5,3) * t268 + t266 * t295;
t171 = V_base(5) * rSges(2,3) - t237 * t250 + t304;
t170 = t239 * t250 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t162 = t237 * V_base(4) - t239 * V_base(5) + V_base(3);
t159 = -rSges(6,3) * t244 + (rSges(6,1) * t267 - rSges(6,2) * t265) * t243;
t158 = -Icges(6,5) * t244 + (Icges(6,1) * t267 - Icges(6,4) * t265) * t243;
t157 = -Icges(6,6) * t244 + (Icges(6,4) * t267 - Icges(6,2) * t265) * t243;
t156 = -Icges(6,3) * t244 + (Icges(6,5) * t267 - Icges(6,6) * t265) * t243;
t154 = -rSges(7,3) * t244 + (rSges(7,1) * t255 - rSges(7,2) * t254) * t243;
t153 = -Icges(7,5) * t244 + (Icges(7,1) * t255 - Icges(7,4) * t254) * t243;
t152 = -Icges(7,6) * t244 + (Icges(7,4) * t255 - Icges(7,2) * t254) * t243;
t151 = -Icges(7,3) * t244 + (Icges(7,5) * t255 - Icges(7,6) * t254) * t243;
t150 = t268 * t306 + t183;
t149 = t266 * t306 + t182;
t147 = -pkin(10) * t244 + t243 * t330;
t142 = rSges(6,1) * t204 + rSges(6,2) * t203 + rSges(6,3) * t320;
t141 = rSges(6,1) * t202 + rSges(6,2) * t201 + rSges(6,3) * t321;
t140 = Icges(6,1) * t204 + Icges(6,4) * t203 + Icges(6,5) * t320;
t139 = Icges(6,1) * t202 + Icges(6,4) * t201 + Icges(6,5) * t321;
t138 = Icges(6,4) * t204 + Icges(6,2) * t203 + Icges(6,6) * t320;
t137 = Icges(6,4) * t202 + Icges(6,2) * t201 + Icges(6,6) * t321;
t136 = Icges(6,5) * t204 + Icges(6,6) * t203 + Icges(6,3) * t320;
t135 = Icges(6,5) * t202 + Icges(6,6) * t201 + Icges(6,3) * t321;
t134 = pkin(5) * t315 + t268 * t283;
t133 = -pkin(5) * t314 + t266 * t283;
t132 = rSges(7,1) * t189 + rSges(7,2) * t188 + rSges(7,3) * t320;
t131 = rSges(7,1) * t187 + rSges(7,2) * t186 + rSges(7,3) * t321;
t130 = Icges(7,1) * t189 + Icges(7,4) * t188 + Icges(7,5) * t320;
t129 = Icges(7,1) * t187 + Icges(7,4) * t186 + Icges(7,5) * t321;
t128 = Icges(7,4) * t189 + Icges(7,2) * t188 + Icges(7,6) * t320;
t127 = Icges(7,4) * t187 + Icges(7,2) * t186 + Icges(7,6) * t321;
t126 = Icges(7,5) * t189 + Icges(7,6) * t188 + Icges(7,3) * t320;
t125 = Icges(7,5) * t187 + Icges(7,6) * t186 + Icges(7,3) * t321;
t124 = t229 * V_base(5) + (-t199 - t236) * t250 + t300;
t123 = t200 * t250 + (-pkin(6) - t229) * V_base(4) + t285;
t122 = t199 * V_base(4) + (-t200 - t238) * V_base(5) + t305;
t121 = t215 * t240 + (-t184 + t310) * t250 + t299;
t120 = t185 * t250 - t215 * t241 + t277;
t119 = t184 * t241 - t185 * t240 + t280;
t118 = t208 * t217 + (-t172 + t301) * t250 + t284;
t117 = t173 * t250 - t208 * t218 + t274;
t116 = t172 * t218 - t173 * t217 + t278;
t115 = -t141 * t216 + t159 * t182 + t276;
t114 = t142 * t216 - t159 * t183 + t273;
t113 = t141 * t183 - t142 * t182 + t275;
t112 = -t131 * t192 - t133 * t216 + t147 * t182 + t149 * t154 + t276;
t111 = t132 * t192 + t134 * t216 - t147 * t183 - t150 * t154 + t273;
t110 = t131 * t150 - t132 * t149 + t133 * t183 - t134 * t182 + t275;
t1 = t241 * (t266 * t281 + t268 * t271) / 0.2e1 + t240 * (t266 * t271 - t268 * t281) / 0.2e1 + t218 * (t266 * t282 + t268 * t272) / 0.2e1 + t217 * (t266 * t272 - t268 * t282) / 0.2e1 + t150 * ((t126 * t320 + t188 * t128 + t189 * t130) * t150 + (t125 * t320 + t127 * t188 + t129 * t189) * t149 + (t151 * t320 + t152 * t188 + t153 * t189) * t192) / 0.2e1 + t182 * ((t136 * t321 + t138 * t201 + t140 * t202) * t183 + (t135 * t321 + t201 * t137 + t202 * t139) * t182 + (t156 * t321 + t157 * t201 + t158 * t202) * t216) / 0.2e1 + t149 * ((t126 * t321 + t128 * t186 + t130 * t187) * t150 + (t125 * t321 + t186 * t127 + t187 * t129) * t149 + (t151 * t321 + t152 * t186 + t153 * t187) * t192) / 0.2e1 + t183 * ((t136 * t320 + t203 * t138 + t204 * t140) * t183 + (t135 * t320 + t137 * t203 + t139 * t204) * t182 + (t156 * t320 + t157 * t203 + t158 * t204) * t216) / 0.2e1 + t216 * ((-t135 * t182 - t136 * t183 - t156 * t216) * t244 + ((-t138 * t265 + t140 * t267) * t183 + (-t137 * t265 + t139 * t267) * t182 + (-t157 * t265 + t158 * t267) * t216) * t243) / 0.2e1 + t192 * ((-t125 * t149 - t126 * t150 - t151 * t192) * t244 + ((-t128 * t254 + t130 * t255) * t150 + (-t127 * t254 + t129 * t255) * t149 + (-t152 * t254 + t153 * t255) * t192) * t243) / 0.2e1 + m(1) * (t221 ^ 2 + t222 ^ 2 + t223 ^ 2) / 0.2e1 + m(2) * (t162 ^ 2 + t170 ^ 2 + t171 ^ 2) / 0.2e1 + m(4) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(3) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(6) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(5) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(7) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + (t231 * t250 + t266 * t279 + t270 * t268 + (-t232 * t266 + t234 * t268 + Icges(1,4)) * V_base(5) + (-t266 * t233 + t268 * t235 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t230 * t250 + t266 * t270 - t268 * t279 + (t268 * t232 + t266 * t234 + Icges(1,2)) * V_base(5) + (t233 * t268 + t235 * t266 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t179 * t248 + t181 * t247) * t241 + (t178 * t248 + t180 * t247) * t240 + (t167 * t244 + t169 * t243) * t218 + (t166 * t244 + t168 * t243) * t217 + (t195 * t263 + t197 * t262 + t230) * V_base(5) + (t196 * t263 + t198 * t262 + t231) * V_base(4) + (t244 * t206 + t243 * t207 + t248 * t213 + t247 * t214 + t263 * t227 + t262 * t228 + Icges(2,3)) * t250) * t250 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
