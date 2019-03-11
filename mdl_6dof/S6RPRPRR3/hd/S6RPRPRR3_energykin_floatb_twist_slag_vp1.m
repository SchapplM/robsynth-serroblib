% Calculate kinetic energy for
% S6RPRPRR3
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
% Datum: 2019-03-09 03:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR3_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:40:28
% EndTime: 2019-03-09 03:40:31
% DurationCPUTime: 3.14s
% Computational Cost: add. (2330->372), mult. (2047->541), div. (0->0), fcn. (1939->12), ass. (0->172)
t252 = sin(qJ(1));
t303 = pkin(1) * t252;
t254 = cos(qJ(1));
t302 = pkin(1) * t254;
t249 = cos(pkin(11));
t300 = t249 * pkin(4);
t299 = -pkin(6) - qJ(2);
t298 = Icges(2,4) * t252;
t247 = qJ(1) + pkin(10);
t236 = sin(t247);
t297 = Icges(3,4) * t236;
t251 = sin(qJ(3));
t296 = Icges(4,4) * t251;
t253 = cos(qJ(3));
t295 = Icges(4,4) * t253;
t248 = sin(pkin(11));
t294 = t236 * t248;
t293 = t236 * t251;
t292 = t236 * t253;
t238 = cos(t247);
t291 = t238 * t248;
t290 = t238 * t251;
t289 = t238 * t253;
t288 = t248 * t253;
t287 = t249 * t253;
t246 = pkin(11) + qJ(5);
t237 = cos(t246);
t285 = pkin(5) * t237;
t283 = qJD(4) * t251;
t282 = qJD(5) * t251;
t281 = qJD(6) * t251;
t240 = V_base(6) + qJD(1);
t280 = t240 * t302 + V_base(2);
t279 = V_base(5) * pkin(6) + V_base(1);
t209 = qJD(3) * t236 + V_base(4);
t201 = pkin(2) * t236 - pkin(7) * t238;
t276 = -t201 - t303;
t235 = sin(t246);
t275 = pkin(5) * t235;
t274 = V_base(5) * qJ(2) + t279;
t273 = V_base(4) * t303 + qJD(2) + V_base(3);
t186 = t238 * t282 + t209;
t270 = pkin(3) * t253 + qJ(4) * t251;
t188 = t270 * t236;
t272 = -t188 + t276;
t208 = -qJD(3) * t238 + V_base(5);
t271 = rSges(4,1) * t253 - rSges(4,2) * t251;
t269 = Icges(4,1) * t253 - t296;
t268 = -Icges(4,2) * t251 + t295;
t267 = Icges(4,5) * t253 - Icges(4,6) * t251;
t224 = pkin(3) * t251 - qJ(4) * t253;
t266 = t208 * t224 + t238 * t283 + t274;
t185 = t236 * t282 + t208;
t265 = (-Icges(4,3) * t238 + t236 * t267) * t208 + (Icges(4,3) * t236 + t238 * t267) * t209 + (Icges(4,5) * t251 + Icges(4,6) * t253) * t240;
t264 = pkin(8) * t251 + t253 * t300;
t263 = pkin(9) * t251 + t253 * t285;
t202 = pkin(2) * t238 + pkin(7) * t236;
t262 = t240 * t202 + t299 * V_base(4) + t280;
t189 = t270 * t238;
t261 = t240 * t189 + t236 * t283 + t262;
t260 = V_base(4) * t201 + (-t202 - t302) * V_base(5) + t273;
t138 = -pkin(4) * t291 + t236 * t264;
t167 = -pkin(8) * t253 + t251 * t300;
t259 = t208 * t167 + (-t138 + t272) * t240 + t266;
t258 = -qJD(4) * t253 + t209 * t188 + t260;
t139 = pkin(4) * t294 + t238 * t264;
t257 = t240 * t139 + (-t167 - t224) * t209 + t261;
t256 = t209 * t138 + (-t139 - t189) * t208 + t258;
t158 = -Icges(4,6) * t238 + t236 * t268;
t159 = Icges(4,6) * t236 + t238 * t268;
t160 = -Icges(4,5) * t238 + t236 * t269;
t161 = Icges(4,5) * t236 + t238 * t269;
t213 = Icges(4,2) * t253 + t296;
t216 = Icges(4,1) * t251 + t295;
t255 = (-t159 * t251 + t161 * t253) * t209 + (-t158 * t251 + t160 * t253) * t208 + (-t213 * t251 + t216 * t253) * t240;
t243 = Icges(2,4) * t254;
t239 = qJ(6) + t246;
t233 = cos(t239);
t232 = sin(t239);
t231 = Icges(3,4) * t238;
t227 = rSges(2,1) * t254 - t252 * rSges(2,2);
t226 = t252 * rSges(2,1) + rSges(2,2) * t254;
t225 = rSges(4,1) * t251 + rSges(4,2) * t253;
t223 = -qJD(5) * t253 + t240;
t218 = Icges(2,1) * t254 - t298;
t217 = Icges(2,1) * t252 + t243;
t215 = -Icges(2,2) * t252 + t243;
t214 = Icges(2,2) * t254 + t298;
t207 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t206 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t205 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t200 = rSges(3,1) * t238 - rSges(3,2) * t236;
t199 = rSges(3,1) * t236 + rSges(3,2) * t238;
t198 = (-qJD(5) - qJD(6)) * t253 + t240;
t197 = Icges(3,1) * t238 - t297;
t196 = Icges(3,1) * t236 + t231;
t195 = -Icges(3,2) * t236 + t231;
t194 = Icges(3,2) * t238 + t297;
t187 = -rSges(5,3) * t253 + (rSges(5,1) * t249 - rSges(5,2) * t248) * t251;
t184 = -Icges(5,5) * t253 + (Icges(5,1) * t249 - Icges(5,4) * t248) * t251;
t183 = -Icges(5,6) * t253 + (Icges(5,4) * t249 - Icges(5,2) * t248) * t251;
t182 = -Icges(5,3) * t253 + (Icges(5,5) * t249 - Icges(5,6) * t248) * t251;
t180 = t238 * t287 + t294;
t179 = t236 * t249 - t238 * t288;
t178 = t236 * t287 - t291;
t177 = -t236 * t288 - t238 * t249;
t176 = -rSges(6,3) * t253 + (rSges(6,1) * t237 - rSges(6,2) * t235) * t251;
t175 = -Icges(6,5) * t253 + (Icges(6,1) * t237 - Icges(6,4) * t235) * t251;
t174 = -Icges(6,6) * t253 + (Icges(6,4) * t237 - Icges(6,2) * t235) * t251;
t173 = -Icges(6,3) * t253 + (Icges(6,5) * t237 - Icges(6,6) * t235) * t251;
t172 = t235 * t236 + t237 * t289;
t171 = -t235 * t289 + t236 * t237;
t170 = -t235 * t238 + t237 * t292;
t169 = -t235 * t292 - t237 * t238;
t166 = rSges(4,3) * t236 + t238 * t271;
t165 = -rSges(4,3) * t238 + t236 * t271;
t164 = -rSges(7,3) * t253 + (rSges(7,1) * t233 - rSges(7,2) * t232) * t251;
t163 = V_base(5) * rSges(2,3) - t226 * t240 + t279;
t162 = t227 * t240 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t155 = -Icges(7,5) * t253 + (Icges(7,1) * t233 - Icges(7,4) * t232) * t251;
t154 = -Icges(7,6) * t253 + (Icges(7,4) * t233 - Icges(7,2) * t232) * t251;
t153 = -Icges(7,3) * t253 + (Icges(7,5) * t233 - Icges(7,6) * t232) * t251;
t152 = t232 * t236 + t233 * t289;
t151 = -t232 * t289 + t233 * t236;
t150 = -t232 * t238 + t233 * t292;
t149 = -t232 * t292 - t233 * t238;
t148 = t226 * V_base(4) - t227 * V_base(5) + V_base(3);
t147 = t238 * t281 + t186;
t146 = t236 * t281 + t185;
t143 = -pkin(9) * t253 + t251 * t285;
t142 = V_base(5) * rSges(3,3) + (-t199 - t303) * t240 + t274;
t141 = t200 * t240 + (-rSges(3,3) + t299) * V_base(4) + t280;
t140 = V_base(4) * t199 + (-t200 - t302) * V_base(5) + t273;
t137 = rSges(5,1) * t180 + rSges(5,2) * t179 + rSges(5,3) * t290;
t136 = rSges(5,1) * t178 + rSges(5,2) * t177 + rSges(5,3) * t293;
t135 = Icges(5,1) * t180 + Icges(5,4) * t179 + Icges(5,5) * t290;
t134 = Icges(5,1) * t178 + Icges(5,4) * t177 + Icges(5,5) * t293;
t133 = Icges(5,4) * t180 + Icges(5,2) * t179 + Icges(5,6) * t290;
t132 = Icges(5,4) * t178 + Icges(5,2) * t177 + Icges(5,6) * t293;
t131 = Icges(5,5) * t180 + Icges(5,6) * t179 + Icges(5,3) * t290;
t130 = Icges(5,5) * t178 + Icges(5,6) * t177 + Icges(5,3) * t293;
t128 = rSges(6,1) * t172 + rSges(6,2) * t171 + rSges(6,3) * t290;
t127 = rSges(6,1) * t170 + rSges(6,2) * t169 + rSges(6,3) * t293;
t126 = Icges(6,1) * t172 + Icges(6,4) * t171 + Icges(6,5) * t290;
t125 = Icges(6,1) * t170 + Icges(6,4) * t169 + Icges(6,5) * t293;
t124 = Icges(6,4) * t172 + Icges(6,2) * t171 + Icges(6,6) * t290;
t123 = Icges(6,4) * t170 + Icges(6,2) * t169 + Icges(6,6) * t293;
t122 = Icges(6,5) * t172 + Icges(6,6) * t171 + Icges(6,3) * t290;
t121 = Icges(6,5) * t170 + Icges(6,6) * t169 + Icges(6,3) * t293;
t119 = rSges(7,1) * t152 + rSges(7,2) * t151 + rSges(7,3) * t290;
t118 = rSges(7,1) * t150 + rSges(7,2) * t149 + rSges(7,3) * t293;
t117 = Icges(7,1) * t152 + Icges(7,4) * t151 + Icges(7,5) * t290;
t116 = Icges(7,1) * t150 + Icges(7,4) * t149 + Icges(7,5) * t293;
t115 = Icges(7,4) * t152 + Icges(7,2) * t151 + Icges(7,6) * t290;
t114 = Icges(7,4) * t150 + Icges(7,2) * t149 + Icges(7,6) * t293;
t113 = Icges(7,5) * t152 + Icges(7,6) * t151 + Icges(7,3) * t290;
t112 = Icges(7,5) * t150 + Icges(7,6) * t149 + Icges(7,3) * t293;
t111 = t236 * t275 + t238 * t263;
t110 = t236 * t263 - t238 * t275;
t109 = t208 * t225 + (-t165 + t276) * t240 + t274;
t108 = t166 * t240 - t209 * t225 + t262;
t107 = t209 * t165 - t208 * t166 + t260;
t106 = t187 * t208 + (-t136 + t272) * t240 + t266;
t105 = t137 * t240 + (-t187 - t224) * t209 + t261;
t104 = t209 * t136 + (-t137 - t189) * t208 + t258;
t103 = -t127 * t223 + t176 * t185 + t259;
t102 = t128 * t223 - t176 * t186 + t257;
t101 = t186 * t127 - t185 * t128 + t256;
t100 = -t110 * t223 - t118 * t198 + t143 * t185 + t146 * t164 + t259;
t99 = t111 * t223 + t119 * t198 - t143 * t186 - t147 * t164 + t257;
t98 = t186 * t110 - t185 * t111 + t147 * t118 - t146 * t119 + t256;
t1 = m(1) * (t205 ^ 2 + t206 ^ 2 + t207 ^ 2) / 0.2e1 + m(2) * (t148 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(3) * (t140 ^ 2 + t141 ^ 2 + t142 ^ 2) / 0.2e1 + m(4) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t101 ^ 2 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + m(5) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + t186 * ((t122 * t290 + t171 * t124 + t172 * t126) * t186 + (t121 * t290 + t123 * t171 + t125 * t172) * t185 + (t171 * t174 + t172 * t175 + t173 * t290) * t223) / 0.2e1 + t147 * ((t113 * t290 + t151 * t115 + t152 * t117) * t147 + (t112 * t290 + t114 * t151 + t116 * t152) * t146 + (t151 * t154 + t152 * t155 + t153 * t290) * t198) / 0.2e1 + t185 * ((t122 * t293 + t124 * t169 + t126 * t170) * t186 + (t121 * t293 + t169 * t123 + t170 * t125) * t185 + (t169 * t174 + t170 * t175 + t173 * t293) * t223) / 0.2e1 + t146 * ((t113 * t293 + t115 * t149 + t117 * t150) * t147 + (t112 * t293 + t149 * t114 + t150 * t116) * t146 + (t149 * t154 + t150 * t155 + t153 * t293) * t198) / 0.2e1 + t223 * ((-t121 * t185 - t122 * t186 - t173 * t223) * t253 + ((-t124 * t235 + t126 * t237) * t186 + (-t123 * t235 + t125 * t237) * t185 + (-t174 * t235 + t175 * t237) * t223) * t251) / 0.2e1 + t198 * ((-t112 * t146 - t113 * t147 - t153 * t198) * t253 + ((-t115 * t232 + t117 * t233) * t147 + (-t114 * t232 + t116 * t233) * t146 + (-t154 * t232 + t155 * t233) * t198) * t251) / 0.2e1 + ((t131 * t293 + t133 * t177 + t135 * t178) * t209 + (t130 * t293 + t177 * t132 + t178 * t134) * t208 + (t177 * t183 + t178 * t184 + t182 * t293) * t240 + t236 * t255 - t265 * t238) * t208 / 0.2e1 + ((t131 * t290 + t179 * t133 + t180 * t135) * t209 + (t130 * t290 + t132 * t179 + t134 * t180) * t208 + (t179 * t183 + t180 * t184 + t182 * t290) * t240 + t236 * t265 + t238 * t255) * t209 / 0.2e1 + ((-t194 * t236 + t196 * t238 - t252 * t214 + t217 * t254 + Icges(1,4)) * V_base(5) + (-t236 * t195 + t238 * t197 - t252 * t215 + t254 * t218 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t238 * t194 + t236 * t196 + t254 * t214 + t252 * t217 + Icges(1,2)) * V_base(5) + (t195 * t238 + t197 * t236 + t215 * t254 + t252 * t218 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t130 * t208 - t131 * t209) * t253 + ((-t133 * t248 + t135 * t249) * t209 + (-t132 * t248 + t134 * t249) * t208) * t251 + (t159 * t253 + t161 * t251) * t209 + (t158 * t253 + t160 * t251) * t208 + (Icges(2,3) + Icges(3,3) + (-t182 + t213) * t253 + (-t183 * t248 + t184 * t249 + t216) * t251) * t240) * t240 / 0.2e1 + t240 * V_base(5) * (Icges(2,5) * t252 + Icges(3,5) * t236 + Icges(2,6) * t254 + Icges(3,6) * t238) + t240 * V_base(4) * (Icges(2,5) * t254 + Icges(3,5) * t238 - Icges(2,6) * t252 - Icges(3,6) * t236) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
