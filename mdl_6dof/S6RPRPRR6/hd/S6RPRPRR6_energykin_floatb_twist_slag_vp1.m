% Calculate kinetic energy for
% S6RPRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:51:08
% EndTime: 2019-03-09 03:51:11
% DurationCPUTime: 3.56s
% Computational Cost: add. (2319->391), mult. (2247->574), div. (0->0), fcn. (2141->12), ass. (0->188)
t263 = sin(pkin(10));
t329 = pkin(2) * t263;
t265 = cos(pkin(10));
t327 = pkin(2) * t265;
t264 = cos(pkin(11));
t326 = t264 * pkin(4);
t268 = sin(qJ(1));
t325 = Icges(2,4) * t268;
t324 = Icges(3,4) * t263;
t323 = Icges(3,4) * t265;
t261 = pkin(10) + qJ(3);
t249 = sin(t261);
t322 = Icges(4,4) * t249;
t251 = cos(t261);
t321 = Icges(4,4) * t251;
t320 = t249 * t268;
t269 = cos(qJ(1));
t319 = t249 * t269;
t318 = t251 * t269;
t262 = sin(pkin(11));
t317 = t262 * t269;
t316 = t264 * t269;
t260 = pkin(11) + qJ(5);
t252 = qJ(6) + t260;
t244 = sin(t252);
t315 = t268 * t244;
t245 = cos(t252);
t314 = t268 * t245;
t248 = sin(t260);
t313 = t268 * t248;
t250 = cos(t260);
t312 = t268 * t250;
t311 = t268 * t262;
t310 = t268 * t264;
t173 = -pkin(7) * t269 + t268 * t327;
t237 = t268 * pkin(1) - qJ(2) * t269;
t307 = -t173 - t237;
t306 = pkin(5) * t250;
t304 = qJD(4) * t249;
t303 = qJD(5) * t249;
t302 = qJD(6) * t249;
t301 = V_base(4) * t237 + V_base(3);
t300 = V_base(5) * pkin(6) + V_base(1);
t291 = pkin(3) * t251 + qJ(4) * t249;
t201 = t291 * t268;
t297 = -t201 + t307;
t242 = qJD(3) * t268 + V_base(4);
t253 = V_base(6) + qJD(1);
t296 = pkin(5) * t248;
t295 = qJD(2) * t268 + t300;
t209 = t269 * t303 + t242;
t294 = V_base(5) * t329 + t295;
t241 = -qJD(3) * t269 + V_base(5);
t293 = rSges(3,1) * t265 - rSges(3,2) * t263;
t292 = rSges(4,1) * t251 - rSges(4,2) * t249;
t290 = Icges(3,1) * t265 - t324;
t289 = Icges(4,1) * t251 - t322;
t288 = -Icges(3,2) * t263 + t323;
t287 = -Icges(4,2) * t249 + t321;
t286 = Icges(3,5) * t265 - Icges(3,6) * t263;
t285 = Icges(4,5) * t251 - Icges(4,6) * t249;
t239 = pkin(1) * t269 + t268 * qJ(2);
t284 = -qJD(2) * t269 + t253 * t239 + V_base(2);
t208 = t268 * t303 + t241;
t214 = pkin(3) * t249 - qJ(4) * t251;
t283 = t241 * t214 + t269 * t304 + t294;
t282 = (-Icges(4,3) * t269 + t268 * t285) * t241 + (Icges(4,3) * t268 + t269 * t285) * t242 + (Icges(4,5) * t249 + Icges(4,6) * t251) * t253;
t281 = pkin(8) * t249 + t251 * t326;
t280 = pkin(9) * t249 + t251 * t306;
t174 = pkin(7) * t268 + t269 * t327;
t279 = V_base(4) * t173 + (-t174 - t239) * V_base(5) + t301;
t278 = (-Icges(3,3) * t269 + t268 * t286) * V_base(5) + (Icges(3,3) * t268 + t269 * t286) * V_base(4) + (Icges(3,5) * t263 + Icges(3,6) * t265) * t253;
t277 = t253 * t174 + (-pkin(6) - t329) * V_base(4) + t284;
t276 = -qJD(4) * t251 + t242 * t201 + t279;
t140 = -pkin(4) * t317 + t268 * t281;
t152 = -pkin(8) * t251 + t249 * t326;
t275 = t241 * t152 + (-t140 + t297) * t253 + t283;
t202 = t291 * t269;
t274 = t253 * t202 + t268 * t304 + t277;
t141 = pkin(4) * t311 + t269 * t281;
t273 = t242 * t140 + (-t141 - t202) * t241 + t276;
t272 = t253 * t141 + (-t152 - t214) * t242 + t274;
t182 = -Icges(4,6) * t269 + t268 * t287;
t183 = Icges(4,6) * t268 + t269 * t287;
t184 = -Icges(4,5) * t269 + t268 * t289;
t185 = Icges(4,5) * t268 + t269 * t289;
t212 = Icges(4,2) * t251 + t322;
t213 = Icges(4,1) * t249 + t321;
t271 = (-t183 * t249 + t185 * t251) * t242 + (-t182 * t249 + t184 * t251) * t241 + (-t212 * t249 + t213 * t251) * t253;
t195 = -Icges(3,6) * t269 + t268 * t288;
t196 = Icges(3,6) * t268 + t269 * t288;
t197 = -Icges(3,5) * t269 + t268 * t290;
t198 = Icges(3,5) * t268 + t269 * t290;
t224 = Icges(3,2) * t265 + t324;
t225 = Icges(3,1) * t263 + t323;
t270 = (-t196 * t263 + t198 * t265) * V_base(4) + (-t195 * t263 + t197 * t265) * V_base(5) + (-t224 * t263 + t225 * t265) * t253;
t257 = Icges(2,4) * t269;
t240 = rSges(2,1) * t269 - t268 * rSges(2,2);
t238 = t268 * rSges(2,1) + rSges(2,2) * t269;
t232 = Icges(2,1) * t269 - t325;
t231 = Icges(2,1) * t268 + t257;
t230 = -Icges(2,2) * t268 + t257;
t229 = Icges(2,2) * t269 + t325;
t228 = Icges(2,5) * t269 - Icges(2,6) * t268;
t227 = Icges(2,5) * t268 + Icges(2,6) * t269;
t226 = rSges(3,1) * t263 + rSges(3,2) * t265;
t222 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t221 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t220 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t218 = -qJD(5) * t251 + t253;
t215 = rSges(4,1) * t249 + rSges(4,2) * t251;
t207 = t251 * t316 + t311;
t206 = -t251 * t317 + t310;
t205 = t251 * t310 - t317;
t204 = -t251 * t311 - t316;
t203 = (-qJD(5) - qJD(6)) * t251 + t253;
t200 = t268 * rSges(3,3) + t269 * t293;
t199 = -rSges(3,3) * t269 + t268 * t293;
t192 = t250 * t318 + t313;
t191 = -t248 * t318 + t312;
t190 = -t248 * t269 + t251 * t312;
t189 = -t250 * t269 - t251 * t313;
t187 = t268 * rSges(4,3) + t269 * t292;
t186 = -rSges(4,3) * t269 + t268 * t292;
t179 = t245 * t318 + t315;
t178 = -t244 * t318 + t314;
t177 = -t244 * t269 + t251 * t314;
t176 = -t245 * t269 - t251 * t315;
t172 = V_base(5) * rSges(2,3) - t238 * t253 + t300;
t171 = t240 * t253 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t170 = -rSges(5,3) * t251 + (rSges(5,1) * t264 - rSges(5,2) * t262) * t249;
t169 = -Icges(5,5) * t251 + (Icges(5,1) * t264 - Icges(5,4) * t262) * t249;
t168 = -Icges(5,6) * t251 + (Icges(5,4) * t264 - Icges(5,2) * t262) * t249;
t167 = -Icges(5,3) * t251 + (Icges(5,5) * t264 - Icges(5,6) * t262) * t249;
t166 = t269 * t302 + t209;
t165 = t268 * t302 + t208;
t163 = t238 * V_base(4) - t240 * V_base(5) + V_base(3);
t160 = -rSges(6,3) * t251 + (rSges(6,1) * t250 - rSges(6,2) * t248) * t249;
t159 = -Icges(6,5) * t251 + (Icges(6,1) * t250 - Icges(6,4) * t248) * t249;
t158 = -Icges(6,6) * t251 + (Icges(6,4) * t250 - Icges(6,2) * t248) * t249;
t157 = -Icges(6,3) * t251 + (Icges(6,5) * t250 - Icges(6,6) * t248) * t249;
t156 = -rSges(7,3) * t251 + (rSges(7,1) * t245 - rSges(7,2) * t244) * t249;
t155 = -Icges(7,5) * t251 + (Icges(7,1) * t245 - Icges(7,4) * t244) * t249;
t154 = -Icges(7,6) * t251 + (Icges(7,4) * t245 - Icges(7,2) * t244) * t249;
t153 = -Icges(7,3) * t251 + (Icges(7,5) * t245 - Icges(7,6) * t244) * t249;
t150 = -pkin(9) * t251 + t249 * t306;
t149 = t207 * rSges(5,1) + t206 * rSges(5,2) + rSges(5,3) * t319;
t148 = rSges(5,1) * t205 + rSges(5,2) * t204 + rSges(5,3) * t320;
t147 = Icges(5,1) * t207 + Icges(5,4) * t206 + Icges(5,5) * t319;
t146 = Icges(5,1) * t205 + Icges(5,4) * t204 + Icges(5,5) * t320;
t145 = Icges(5,4) * t207 + Icges(5,2) * t206 + Icges(5,6) * t319;
t144 = Icges(5,4) * t205 + Icges(5,2) * t204 + Icges(5,6) * t320;
t143 = Icges(5,5) * t207 + Icges(5,6) * t206 + Icges(5,3) * t319;
t142 = Icges(5,5) * t205 + Icges(5,6) * t204 + Icges(5,3) * t320;
t138 = t192 * rSges(6,1) + t191 * rSges(6,2) + rSges(6,3) * t319;
t137 = rSges(6,1) * t190 + rSges(6,2) * t189 + rSges(6,3) * t320;
t136 = Icges(6,1) * t192 + Icges(6,4) * t191 + Icges(6,5) * t319;
t135 = Icges(6,1) * t190 + Icges(6,4) * t189 + Icges(6,5) * t320;
t134 = Icges(6,4) * t192 + Icges(6,2) * t191 + Icges(6,6) * t319;
t133 = Icges(6,4) * t190 + Icges(6,2) * t189 + Icges(6,6) * t320;
t132 = Icges(6,5) * t192 + Icges(6,6) * t191 + Icges(6,3) * t319;
t131 = Icges(6,5) * t190 + Icges(6,6) * t189 + Icges(6,3) * t320;
t129 = t179 * rSges(7,1) + t178 * rSges(7,2) + rSges(7,3) * t319;
t128 = rSges(7,1) * t177 + rSges(7,2) * t176 + rSges(7,3) * t320;
t127 = Icges(7,1) * t179 + Icges(7,4) * t178 + Icges(7,5) * t319;
t126 = Icges(7,1) * t177 + Icges(7,4) * t176 + Icges(7,5) * t320;
t125 = Icges(7,4) * t179 + Icges(7,2) * t178 + Icges(7,6) * t319;
t124 = Icges(7,4) * t177 + Icges(7,2) * t176 + Icges(7,6) * t320;
t123 = Icges(7,5) * t179 + Icges(7,6) * t178 + Icges(7,3) * t319;
t122 = Icges(7,5) * t177 + Icges(7,6) * t176 + Icges(7,3) * t320;
t121 = t226 * V_base(5) + (-t199 - t237) * t253 + t295;
t120 = t253 * t200 + (-pkin(6) - t226) * V_base(4) + t284;
t119 = t199 * V_base(4) + (-t200 - t239) * V_base(5) + t301;
t118 = t268 * t296 + t269 * t280;
t117 = t268 * t280 - t269 * t296;
t116 = t215 * t241 + (-t186 + t307) * t253 + t294;
t115 = t253 * t187 - t242 * t215 + t277;
t114 = t186 * t242 - t187 * t241 + t279;
t113 = t170 * t241 + (-t148 + t297) * t253 + t283;
t112 = t253 * t149 + (-t170 - t214) * t242 + t274;
t111 = t148 * t242 + (-t149 - t202) * t241 + t276;
t110 = -t137 * t218 + t160 * t208 + t275;
t109 = t218 * t138 - t209 * t160 + t272;
t108 = t137 * t209 - t138 * t208 + t273;
t107 = -t117 * t218 - t128 * t203 + t150 * t208 + t156 * t165 + t275;
t106 = t218 * t118 + t203 * t129 - t209 * t150 - t166 * t156 + t272;
t105 = t117 * t209 - t118 * t208 + t128 * t166 - t129 * t165 + t273;
t1 = t203 * ((-t122 * t165 - t123 * t166 - t153 * t203) * t251 + ((-t125 * t244 + t127 * t245) * t166 + (-t124 * t244 + t126 * t245) * t165 + (-t154 * t244 + t155 * t245) * t203) * t249) / 0.2e1 + m(1) * (t220 ^ 2 + t221 ^ 2 + t222 ^ 2) / 0.2e1 + m(2) * (t163 ^ 2 + t171 ^ 2 + t172 ^ 2) / 0.2e1 + m(3) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(6) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(5) * (t111 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(4) * (t114 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(7) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + t209 * ((t132 * t319 + t191 * t134 + t192 * t136) * t209 + (t131 * t319 + t191 * t133 + t192 * t135) * t208 + (t157 * t319 + t191 * t158 + t192 * t159) * t218) / 0.2e1 + t166 * ((t123 * t319 + t178 * t125 + t179 * t127) * t166 + (t122 * t319 + t178 * t124 + t179 * t126) * t165 + (t153 * t319 + t178 * t154 + t179 * t155) * t203) / 0.2e1 + t208 * ((t132 * t320 + t134 * t189 + t136 * t190) * t209 + (t131 * t320 + t133 * t189 + t135 * t190) * t208 + (t157 * t320 + t158 * t189 + t159 * t190) * t218) / 0.2e1 + t165 * ((t123 * t320 + t125 * t176 + t127 * t177) * t166 + (t122 * t320 + t124 * t176 + t126 * t177) * t165 + (t153 * t320 + t154 * t176 + t155 * t177) * t203) / 0.2e1 + t218 * ((-t131 * t208 - t132 * t209 - t157 * t218) * t251 + ((-t134 * t248 + t136 * t250) * t209 + (-t133 * t248 + t135 * t250) * t208 + (-t158 * t248 + t159 * t250) * t218) * t249) / 0.2e1 + ((t143 * t320 + t145 * t204 + t147 * t205) * t242 + (t142 * t320 + t144 * t204 + t146 * t205) * t241 + (t167 * t320 + t168 * t204 + t169 * t205) * t253 + t268 * t271 - t282 * t269) * t241 / 0.2e1 + ((t143 * t319 + t206 * t145 + t207 * t147) * t242 + (t142 * t319 + t206 * t144 + t207 * t146) * t241 + (t167 * t319 + t206 * t168 + t207 * t169) * t253 + t268 * t282 + t269 * t271) * t242 / 0.2e1 + (t228 * t253 + t268 * t278 + t270 * t269 + (-t268 * t229 + t231 * t269 + Icges(1,4)) * V_base(5) + (-t268 * t230 + t232 * t269 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t227 * t253 + t268 * t270 - t278 * t269 + (t229 * t269 + t268 * t231 + Icges(1,2)) * V_base(5) + (t230 * t269 + t268 * t232 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t142 * t241 - t143 * t242) * t251 + ((-t145 * t262 + t147 * t264) * t242 + (-t144 * t262 + t146 * t264) * t241) * t249 + (t183 * t251 + t185 * t249) * t242 + (t182 * t251 + t184 * t249) * t241 + (t195 * t265 + t197 * t263 + t227) * V_base(5) + (t196 * t265 + t198 * t263 + t228) * V_base(4) + (t224 * t265 + t225 * t263 + Icges(2,3) + (-t167 + t212) * t251 + (-t168 * t262 + t169 * t264 + t213) * t249) * t253) * t253 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
