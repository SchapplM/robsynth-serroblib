% Calculate kinetic energy for
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-03-09 02:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:17:53
% EndTime: 2019-03-09 02:17:55
% DurationCPUTime: 2.56s
% Computational Cost: add. (2109->325), mult. (1486->459), div. (0->0), fcn. (1266->12), ass. (0->171)
t229 = qJ(1) + pkin(10);
t218 = sin(t229);
t220 = cos(t229);
t234 = sin(qJ(1));
t236 = cos(qJ(1));
t301 = Icges(2,5) * t234 + Icges(3,5) * t218 + Icges(2,6) * t236 + Icges(3,6) * t220;
t300 = Icges(2,5) * t236 + Icges(3,5) * t220 - Icges(2,6) * t234 - Icges(3,6) * t218;
t298 = pkin(1) * t234;
t297 = pkin(1) * t236;
t230 = sin(pkin(11));
t296 = pkin(3) * t230;
t228 = pkin(11) + qJ(4);
t217 = sin(t228);
t295 = pkin(4) * t217;
t231 = cos(pkin(11));
t294 = t231 * pkin(3);
t293 = -pkin(6) - qJ(2);
t292 = Icges(2,4) * t234;
t291 = Icges(3,4) * t218;
t290 = Icges(4,4) * t230;
t289 = Icges(4,4) * t231;
t288 = Icges(5,4) * t217;
t219 = cos(t228);
t287 = Icges(5,4) * t219;
t221 = qJ(5) + t228;
t214 = sin(t221);
t286 = Icges(6,4) * t214;
t215 = cos(t221);
t285 = Icges(6,4) * t215;
t284 = t214 * t218;
t283 = t214 * t220;
t233 = sin(qJ(6));
t282 = t218 * t233;
t235 = cos(qJ(6));
t281 = t218 * t235;
t280 = t220 * t233;
t279 = t220 * t235;
t277 = pkin(4) * t219;
t275 = qJD(6) * t214;
t222 = V_base(6) + qJD(1);
t274 = t222 * t297 + V_base(2);
t273 = V_base(5) * pkin(6) + V_base(1);
t197 = qJD(4) * t218 + V_base(4);
t184 = pkin(2) * t218 - qJ(3) * t220;
t270 = -t184 - t298;
t186 = pkin(2) * t220 + qJ(3) * t218;
t269 = -t186 - t297;
t268 = V_base(5) * qJ(2) + t273;
t267 = V_base(4) * t298 + qJD(2) + V_base(3);
t166 = qJD(5) * t218 + t197;
t123 = -pkin(7) * t220 + t294 * t218;
t266 = -t123 + t270;
t265 = qJD(3) * t218 + t268;
t264 = pkin(5) * t215 + pkin(9) * t214;
t263 = V_base(4) * t184 + t267;
t262 = rSges(4,1) * t231 - rSges(4,2) * t230;
t261 = rSges(5,1) * t219 - rSges(5,2) * t217;
t260 = rSges(6,1) * t215 - rSges(6,2) * t214;
t259 = Icges(4,1) * t231 - t290;
t258 = Icges(5,1) * t219 - t288;
t257 = Icges(6,1) * t215 - t286;
t256 = -Icges(4,2) * t230 + t289;
t255 = -Icges(5,2) * t217 + t287;
t254 = -Icges(6,2) * t214 + t285;
t253 = Icges(4,5) * t231 - Icges(4,6) * t230;
t252 = Icges(5,5) * t219 - Icges(5,6) * t217;
t251 = Icges(6,5) * t215 - Icges(6,6) * t214;
t117 = -pkin(8) * t220 + t277 * t218;
t250 = -t117 + t266;
t249 = V_base(5) * t296 + t265;
t248 = -qJD(3) * t220 + t222 * t186 + t274;
t165 = V_base(5) + (-qJD(4) - qJD(5)) * t220;
t196 = -qJD(4) * t220 + V_base(5);
t247 = t196 * t295 + t249;
t246 = (-Icges(6,3) * t220 + t251 * t218) * t165 + (Icges(6,3) * t218 + t251 * t220) * t166 + (Icges(6,5) * t214 + Icges(6,6) * t215) * t222;
t245 = (-Icges(5,3) * t220 + t252 * t218) * t196 + (Icges(5,3) * t218 + t252 * t220) * t197 + (Icges(5,5) * t217 + Icges(5,6) * t219) * t222;
t244 = (-Icges(4,3) * t220 + t253 * t218) * V_base(5) + (Icges(4,3) * t218 + t253 * t220) * V_base(4) + (Icges(4,5) * t230 + Icges(4,6) * t231) * t222;
t124 = pkin(7) * t218 + t294 * t220;
t243 = V_base(4) * t123 + (-t124 + t269) * V_base(5) + t263;
t242 = t222 * t124 + (t293 - t296) * V_base(4) + t248;
t118 = pkin(8) * t218 + t277 * t220;
t241 = t197 * t117 - t196 * t118 + t243;
t240 = t222 * t118 - t197 * t295 + t242;
t127 = -Icges(6,6) * t220 + t254 * t218;
t128 = Icges(6,6) * t218 + t254 * t220;
t129 = -Icges(6,5) * t220 + t257 * t218;
t130 = Icges(6,5) * t218 + t257 * t220;
t169 = Icges(6,2) * t215 + t286;
t170 = Icges(6,1) * t214 + t285;
t239 = (-t128 * t214 + t130 * t215) * t166 + (-t127 * t214 + t129 * t215) * t165 + (-t169 * t214 + t170 * t215) * t222;
t135 = -Icges(5,6) * t220 + t255 * t218;
t136 = Icges(5,6) * t218 + t255 * t220;
t137 = -Icges(5,5) * t220 + t258 * t218;
t138 = Icges(5,5) * t218 + t258 * t220;
t177 = Icges(5,2) * t219 + t288;
t180 = Icges(5,1) * t217 + t287;
t238 = (-t136 * t217 + t138 * t219) * t197 + (-t135 * t217 + t137 * t219) * t196 + (-t177 * t217 + t180 * t219) * t222;
t150 = -Icges(4,6) * t220 + t256 * t218;
t151 = Icges(4,6) * t218 + t256 * t220;
t152 = -Icges(4,5) * t220 + t259 * t218;
t153 = Icges(4,5) * t218 + t259 * t220;
t194 = Icges(4,2) * t231 + t290;
t195 = Icges(4,1) * t230 + t289;
t237 = (-t151 * t230 + t153 * t231) * V_base(4) + (-t150 * t230 + t152 * t231) * V_base(5) + (-t194 * t230 + t195 * t231) * t222;
t225 = Icges(2,4) * t236;
t213 = Icges(3,4) * t220;
t206 = rSges(2,1) * t236 - t234 * rSges(2,2);
t205 = t234 * rSges(2,1) + rSges(2,2) * t236;
t204 = Icges(2,1) * t236 - t292;
t203 = Icges(2,1) * t234 + t225;
t202 = -Icges(2,2) * t234 + t225;
t201 = Icges(2,2) * t236 + t292;
t198 = rSges(4,1) * t230 + rSges(4,2) * t231;
t192 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t191 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t190 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t188 = -qJD(6) * t215 + t222;
t187 = rSges(3,1) * t220 - rSges(3,2) * t218;
t185 = rSges(3,1) * t218 + rSges(3,2) * t220;
t183 = rSges(5,1) * t217 + rSges(5,2) * t219;
t182 = Icges(3,1) * t220 - t291;
t181 = Icges(3,1) * t218 + t213;
t179 = -Icges(3,2) * t218 + t213;
t178 = Icges(3,2) * t220 + t291;
t173 = pkin(5) * t214 - pkin(9) * t215;
t172 = rSges(6,1) * t214 + rSges(6,2) * t215;
t163 = t215 * t279 + t282;
t162 = -t215 * t280 + t281;
t161 = t215 * t281 - t280;
t160 = -t215 * t282 - t279;
t159 = t264 * t220;
t158 = t264 * t218;
t157 = V_base(5) * rSges(2,3) - t205 * t222 + t273;
t156 = t206 * t222 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t155 = rSges(4,3) * t218 + t262 * t220;
t154 = -rSges(4,3) * t220 + t262 * t218;
t147 = t205 * V_base(4) - t206 * V_base(5) + V_base(3);
t146 = -rSges(7,3) * t215 + (rSges(7,1) * t235 - rSges(7,2) * t233) * t214;
t145 = -Icges(7,5) * t215 + (Icges(7,1) * t235 - Icges(7,4) * t233) * t214;
t144 = -Icges(7,6) * t215 + (Icges(7,4) * t235 - Icges(7,2) * t233) * t214;
t143 = -Icges(7,3) * t215 + (Icges(7,5) * t235 - Icges(7,6) * t233) * t214;
t142 = rSges(5,3) * t218 + t261 * t220;
t141 = -rSges(5,3) * t220 + t261 * t218;
t140 = t220 * t275 + t166;
t139 = t218 * t275 + t165;
t132 = rSges(6,3) * t218 + t260 * t220;
t131 = -rSges(6,3) * t220 + t260 * t218;
t120 = V_base(5) * rSges(3,3) + (-t185 - t298) * t222 + t268;
t119 = t187 * t222 + (-rSges(3,3) + t293) * V_base(4) + t274;
t116 = V_base(4) * t185 + (-t187 - t297) * V_base(5) + t267;
t113 = rSges(7,1) * t163 + rSges(7,2) * t162 + rSges(7,3) * t283;
t112 = rSges(7,1) * t161 + rSges(7,2) * t160 + rSges(7,3) * t284;
t111 = Icges(7,1) * t163 + Icges(7,4) * t162 + Icges(7,5) * t283;
t110 = Icges(7,1) * t161 + Icges(7,4) * t160 + Icges(7,5) * t284;
t109 = Icges(7,4) * t163 + Icges(7,2) * t162 + Icges(7,6) * t283;
t108 = Icges(7,4) * t161 + Icges(7,2) * t160 + Icges(7,6) * t284;
t107 = Icges(7,5) * t163 + Icges(7,6) * t162 + Icges(7,3) * t283;
t106 = Icges(7,5) * t161 + Icges(7,6) * t160 + Icges(7,3) * t284;
t105 = t198 * V_base(5) + (-t154 + t270) * t222 + t265;
t104 = t155 * t222 + (-t198 + t293) * V_base(4) + t248;
t103 = V_base(4) * t154 + (-t155 + t269) * V_base(5) + t263;
t102 = t183 * t196 + (-t141 + t266) * t222 + t249;
t101 = t142 * t222 - t183 * t197 + t242;
t100 = t197 * t141 - t196 * t142 + t243;
t99 = t165 * t172 + (-t131 + t250) * t222 + t247;
t98 = t132 * t222 - t166 * t172 + t240;
t97 = t166 * t131 - t165 * t132 + t241;
t96 = -t112 * t188 + t139 * t146 + t165 * t173 + (-t158 + t250) * t222 + t247;
t95 = t113 * t188 - t140 * t146 + t159 * t222 - t166 * t173 + t240;
t94 = t140 * t112 - t139 * t113 + t166 * t158 - t165 * t159 + t241;
t1 = t197 * (t245 * t218 + t238 * t220) / 0.2e1 + t196 * (t238 * t218 - t245 * t220) / 0.2e1 + t166 * (t246 * t218 + t239 * t220) / 0.2e1 + t165 * (t239 * t218 - t246 * t220) / 0.2e1 + t140 * ((t107 * t283 + t162 * t109 + t163 * t111) * t140 + (t106 * t283 + t108 * t162 + t110 * t163) * t139 + (t143 * t283 + t144 * t162 + t145 * t163) * t188) / 0.2e1 + t139 * ((t107 * t284 + t109 * t160 + t111 * t161) * t140 + (t106 * t284 + t160 * t108 + t161 * t110) * t139 + (t143 * t284 + t144 * t160 + t145 * t161) * t188) / 0.2e1 + m(3) * (t116 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(4) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(7) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + m(6) * (t97 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t101 ^ 2 + t102 ^ 2) / 0.2e1 + t188 * ((-t106 * t139 - t107 * t140 - t143 * t188) * t215 + ((-t109 * t233 + t111 * t235) * t140 + (-t108 * t233 + t110 * t235) * t139 + (-t144 * t233 + t145 * t235) * t188) * t214) / 0.2e1 + m(1) * (t190 ^ 2 + t191 ^ 2 + t192 ^ 2) / 0.2e1 + m(2) * (t147 ^ 2 + t156 ^ 2 + t157 ^ 2) / 0.2e1 + (t244 * t218 + t237 * t220 + t300 * t222 + (-t178 * t218 + t181 * t220 - t234 * t201 + t203 * t236 + Icges(1,4)) * V_base(5) + (-t218 * t179 + t220 * t182 - t234 * t202 + t236 * t204 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t237 * t218 - t244 * t220 + t301 * t222 + (t220 * t178 + t218 * t181 + t236 * t201 + t234 * t203 + Icges(1,2)) * V_base(5) + (t179 * t220 + t182 * t218 + t202 * t236 + t234 * t204 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t136 * t219 + t138 * t217) * t197 + (t135 * t219 + t137 * t217) * t196 + (t128 * t215 + t130 * t214) * t166 + (t127 * t215 + t129 * t214) * t165 + (t150 * t231 + t152 * t230 + t301) * V_base(5) + (t151 * t231 + t153 * t230 + t300) * V_base(4) + (t215 * t169 + t214 * t170 + t219 * t177 + t217 * t180 + t231 * t194 + t230 * t195 + Icges(2,3) + Icges(3,3)) * t222) * t222 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
