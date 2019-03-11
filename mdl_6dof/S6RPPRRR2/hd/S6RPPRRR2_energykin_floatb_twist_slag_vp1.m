% Calculate kinetic energy for
% S6RPPRRR2
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
% Datum: 2019-03-09 02:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPRRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPRRR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPRRR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:20:14
% EndTime: 2019-03-09 02:20:17
% DurationCPUTime: 2.48s
% Computational Cost: add. (2244->349), mult. (1736->499), div. (0->0), fcn. (1572->12), ass. (0->176)
t245 = qJ(1) + pkin(10);
t235 = sin(t245);
t237 = cos(t245);
t251 = sin(qJ(1));
t253 = cos(qJ(1));
t315 = Icges(2,5) * t251 + Icges(3,5) * t235 + Icges(2,6) * t253 + Icges(3,6) * t237;
t314 = Icges(2,5) * t253 + Icges(3,5) * t237 - Icges(2,6) * t251 - Icges(3,6) * t235;
t312 = pkin(1) * t251;
t311 = pkin(1) * t253;
t247 = sin(pkin(11));
t310 = pkin(3) * t247;
t248 = cos(pkin(11));
t309 = pkin(3) * t248;
t252 = cos(qJ(5));
t308 = pkin(5) * t252;
t307 = -pkin(6) - qJ(2);
t305 = Icges(2,4) * t251;
t304 = Icges(3,4) * t235;
t303 = Icges(4,4) * t247;
t302 = Icges(4,4) * t248;
t244 = pkin(11) + qJ(4);
t234 = sin(t244);
t301 = Icges(5,4) * t234;
t236 = cos(t244);
t300 = Icges(5,4) * t236;
t299 = t234 * t235;
t298 = t234 * t237;
t246 = qJ(5) + qJ(6);
t239 = sin(t246);
t297 = t235 * t239;
t240 = cos(t246);
t296 = t235 * t240;
t250 = sin(qJ(5));
t295 = t235 * t250;
t294 = t235 * t252;
t293 = t237 * t239;
t292 = t237 * t240;
t291 = t237 * t250;
t290 = t237 * t252;
t288 = qJD(5) * t234;
t287 = qJD(6) * t234;
t238 = V_base(6) + qJD(1);
t286 = t238 * t311 + V_base(2);
t285 = V_base(5) * pkin(6) + V_base(1);
t216 = qJD(4) * t235 + V_base(4);
t201 = pkin(2) * t235 - qJ(3) * t237;
t282 = -t201 - t312;
t203 = pkin(2) * t237 + qJ(3) * t235;
t281 = -t203 - t311;
t280 = V_base(5) * qJ(2) + t285;
t279 = V_base(4) * t312 + qJD(2) + V_base(3);
t183 = t237 * t288 + t216;
t142 = -pkin(7) * t237 + t235 * t309;
t278 = -t142 + t282;
t277 = qJD(3) * t235 + t280;
t276 = pkin(4) * t236 + pkin(8) * t234;
t275 = V_base(4) * t201 + t279;
t215 = -qJD(4) * t237 + V_base(5);
t274 = rSges(4,1) * t248 - rSges(4,2) * t247;
t273 = rSges(5,1) * t236 - rSges(5,2) * t234;
t272 = Icges(4,1) * t248 - t303;
t271 = Icges(5,1) * t236 - t301;
t270 = -Icges(4,2) * t247 + t302;
t269 = -Icges(5,2) * t234 + t300;
t268 = Icges(4,5) * t248 - Icges(4,6) * t247;
t267 = Icges(5,5) * t236 - Icges(5,6) * t234;
t266 = V_base(5) * t310 + t277;
t182 = t235 * t288 + t215;
t265 = -qJD(3) * t237 + t238 * t203 + t286;
t264 = pkin(9) * t234 + t236 * t308;
t263 = (-Icges(5,3) * t237 + t235 * t267) * t215 + (Icges(5,3) * t235 + t237 * t267) * t216 + (Icges(5,5) * t234 + Icges(5,6) * t236) * t238;
t262 = (-Icges(4,3) * t237 + t235 * t268) * V_base(5) + (Icges(4,3) * t235 + t237 * t268) * V_base(4) + (Icges(4,5) * t247 + Icges(4,6) * t248) * t238;
t143 = pkin(7) * t235 + t237 * t309;
t261 = V_base(4) * t142 + (-t143 + t281) * V_base(5) + t275;
t180 = t276 * t235;
t205 = t234 * pkin(4) - t236 * pkin(8);
t260 = t215 * t205 + (-t180 + t278) * t238 + t266;
t259 = t238 * t143 + (t307 - t310) * V_base(4) + t265;
t181 = t276 * t237;
t258 = t216 * t180 - t181 * t215 + t261;
t257 = t238 * t181 - t205 * t216 + t259;
t149 = -Icges(5,6) * t237 + t235 * t269;
t150 = Icges(5,6) * t235 + t237 * t269;
t151 = -Icges(5,5) * t237 + t235 * t271;
t152 = Icges(5,5) * t235 + t237 * t271;
t194 = Icges(5,2) * t236 + t301;
t197 = Icges(5,1) * t234 + t300;
t256 = (-t150 * t234 + t152 * t236) * t216 + (-t149 * t234 + t151 * t236) * t215 + (-t194 * t234 + t197 * t236) * t238;
t163 = -Icges(4,6) * t237 + t235 * t270;
t164 = Icges(4,6) * t235 + t237 * t270;
t165 = -Icges(4,5) * t237 + t235 * t272;
t166 = Icges(4,5) * t235 + t237 * t272;
t213 = Icges(4,2) * t248 + t303;
t214 = Icges(4,1) * t247 + t302;
t255 = (-t164 * t247 + t166 * t248) * V_base(4) + (-t163 * t247 + t165 * t248) * V_base(5) + (-t213 * t247 + t214 * t248) * t238;
t242 = Icges(2,4) * t253;
t231 = Icges(3,4) * t237;
t225 = rSges(2,1) * t253 - rSges(2,2) * t251;
t224 = rSges(2,1) * t251 + rSges(2,2) * t253;
t223 = Icges(2,1) * t253 - t305;
t222 = Icges(2,1) * t251 + t242;
t221 = -Icges(2,2) * t251 + t242;
t220 = Icges(2,2) * t253 + t305;
t217 = rSges(4,1) * t247 + rSges(4,2) * t248;
t211 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t210 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t209 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t208 = -qJD(5) * t236 + t238;
t204 = rSges(3,1) * t237 - rSges(3,2) * t235;
t202 = rSges(3,1) * t235 + rSges(3,2) * t237;
t200 = rSges(5,1) * t234 + rSges(5,2) * t236;
t199 = Icges(3,1) * t237 - t304;
t198 = Icges(3,1) * t235 + t231;
t196 = -Icges(3,2) * t235 + t231;
t195 = Icges(3,2) * t237 + t304;
t189 = (-qJD(5) - qJD(6)) * t236 + t238;
t187 = t236 * t290 + t295;
t186 = -t236 * t291 + t294;
t185 = t236 * t294 - t291;
t184 = -t236 * t295 - t290;
t178 = t236 * t292 + t297;
t177 = -t236 * t293 + t296;
t176 = t236 * t296 - t293;
t175 = -t236 * t297 - t292;
t174 = -rSges(6,3) * t236 + (rSges(6,1) * t252 - rSges(6,2) * t250) * t234;
t173 = V_base(5) * rSges(2,3) - t224 * t238 + t285;
t172 = t225 * t238 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t171 = -Icges(6,5) * t236 + (Icges(6,1) * t252 - Icges(6,4) * t250) * t234;
t170 = -Icges(6,6) * t236 + (Icges(6,4) * t252 - Icges(6,2) * t250) * t234;
t169 = -Icges(6,3) * t236 + (Icges(6,5) * t252 - Icges(6,6) * t250) * t234;
t168 = rSges(4,3) * t235 + t237 * t274;
t167 = -rSges(4,3) * t237 + t235 * t274;
t160 = t224 * V_base(4) - t225 * V_base(5) + V_base(3);
t158 = -rSges(7,3) * t236 + (rSges(7,1) * t240 - rSges(7,2) * t239) * t234;
t157 = -Icges(7,5) * t236 + (Icges(7,1) * t240 - Icges(7,4) * t239) * t234;
t156 = -Icges(7,6) * t236 + (Icges(7,4) * t240 - Icges(7,2) * t239) * t234;
t155 = -Icges(7,3) * t236 + (Icges(7,5) * t240 - Icges(7,6) * t239) * t234;
t154 = rSges(5,3) * t235 + t237 * t273;
t153 = -rSges(5,3) * t237 + t235 * t273;
t146 = -pkin(9) * t236 + t234 * t308;
t145 = t237 * t287 + t183;
t144 = t235 * t287 + t182;
t138 = V_base(5) * rSges(3,3) + (-t202 - t312) * t238 + t280;
t137 = t204 * t238 + (-rSges(3,3) + t307) * V_base(4) + t286;
t136 = t202 * V_base(4) + (-t204 - t311) * V_base(5) + t279;
t135 = rSges(6,1) * t187 + rSges(6,2) * t186 + rSges(6,3) * t298;
t134 = rSges(6,1) * t185 + rSges(6,2) * t184 + rSges(6,3) * t299;
t133 = Icges(6,1) * t187 + Icges(6,4) * t186 + Icges(6,5) * t298;
t132 = Icges(6,1) * t185 + Icges(6,4) * t184 + Icges(6,5) * t299;
t131 = Icges(6,4) * t187 + Icges(6,2) * t186 + Icges(6,6) * t298;
t130 = Icges(6,4) * t185 + Icges(6,2) * t184 + Icges(6,6) * t299;
t129 = Icges(6,5) * t187 + Icges(6,6) * t186 + Icges(6,3) * t298;
t128 = Icges(6,5) * t185 + Icges(6,6) * t184 + Icges(6,3) * t299;
t127 = pkin(5) * t295 + t237 * t264;
t126 = -pkin(5) * t291 + t235 * t264;
t125 = rSges(7,1) * t178 + rSges(7,2) * t177 + rSges(7,3) * t298;
t124 = rSges(7,1) * t176 + rSges(7,2) * t175 + rSges(7,3) * t299;
t123 = Icges(7,1) * t178 + Icges(7,4) * t177 + Icges(7,5) * t298;
t122 = Icges(7,1) * t176 + Icges(7,4) * t175 + Icges(7,5) * t299;
t121 = Icges(7,4) * t178 + Icges(7,2) * t177 + Icges(7,6) * t298;
t120 = Icges(7,4) * t176 + Icges(7,2) * t175 + Icges(7,6) * t299;
t119 = Icges(7,5) * t178 + Icges(7,6) * t177 + Icges(7,3) * t298;
t118 = Icges(7,5) * t176 + Icges(7,6) * t175 + Icges(7,3) * t299;
t117 = t217 * V_base(5) + (-t167 + t282) * t238 + t277;
t116 = t168 * t238 + (-t217 + t307) * V_base(4) + t265;
t115 = t167 * V_base(4) + (-t168 + t281) * V_base(5) + t275;
t114 = t200 * t215 + (-t153 + t278) * t238 + t266;
t113 = t154 * t238 - t200 * t216 + t259;
t112 = t153 * t216 - t154 * t215 + t261;
t111 = -t134 * t208 + t174 * t182 + t260;
t110 = t135 * t208 - t174 * t183 + t257;
t109 = t134 * t183 - t135 * t182 + t258;
t108 = -t124 * t189 - t126 * t208 + t144 * t158 + t146 * t182 + t260;
t107 = t125 * t189 + t127 * t208 - t145 * t158 - t146 * t183 + t257;
t106 = t124 * t145 - t125 * t144 + t126 * t183 - t127 * t182 + t258;
t1 = m(3) * (t136 ^ 2 + t137 ^ 2 + t138 ^ 2) / 0.2e1 + m(1) * (t209 ^ 2 + t210 ^ 2 + t211 ^ 2) / 0.2e1 + t182 * ((t129 * t299 + t131 * t184 + t133 * t185) * t183 + (t128 * t299 + t184 * t130 + t185 * t132) * t182 + (t169 * t299 + t170 * t184 + t171 * t185) * t208) / 0.2e1 + t144 * ((t119 * t299 + t121 * t175 + t123 * t176) * t145 + (t118 * t299 + t175 * t120 + t176 * t122) * t144 + (t155 * t299 + t156 * t175 + t157 * t176) * t189) / 0.2e1 + m(4) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(7) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(6) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(5) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + t183 * ((t129 * t298 + t186 * t131 + t187 * t133) * t183 + (t128 * t298 + t130 * t186 + t132 * t187) * t182 + (t169 * t298 + t170 * t186 + t171 * t187) * t208) / 0.2e1 + t145 * ((t119 * t298 + t177 * t121 + t178 * t123) * t145 + (t118 * t298 + t120 * t177 + t122 * t178) * t144 + (t155 * t298 + t156 * t177 + t157 * t178) * t189) / 0.2e1 + t208 * ((-t128 * t182 - t129 * t183 - t169 * t208) * t236 + ((-t131 * t250 + t133 * t252) * t183 + (-t130 * t250 + t132 * t252) * t182 + (-t170 * t250 + t171 * t252) * t208) * t234) / 0.2e1 + t216 * (t235 * t263 + t237 * t256) / 0.2e1 + t215 * (t235 * t256 - t237 * t263) / 0.2e1 + t189 * ((-t118 * t144 - t119 * t145 - t155 * t189) * t236 + ((-t121 * t239 + t123 * t240) * t145 + (-t120 * t239 + t122 * t240) * t144 + (-t156 * t239 + t157 * t240) * t189) * t234) / 0.2e1 + m(2) * (t160 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + ((t150 * t236 + t152 * t234) * t216 + (t149 * t236 + t151 * t234) * t215 + (t163 * t248 + t165 * t247 + t315) * V_base(5) + (t164 * t248 + t166 * t247 + t314) * V_base(4) + (t194 * t236 + t197 * t234 + t213 * t248 + t214 * t247 + Icges(2,3) + Icges(3,3)) * t238) * t238 / 0.2e1 + (t235 * t262 + t237 * t255 + t314 * t238 + (-t195 * t235 + t198 * t237 - t220 * t251 + t222 * t253 + Icges(1,4)) * V_base(5) + (-t196 * t235 + t199 * t237 - t221 * t251 + t223 * t253 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t235 * t255 - t237 * t262 + t315 * t238 + (t195 * t237 + t198 * t235 + t220 * t253 + t222 * t251 + Icges(1,2)) * V_base(5) + (t196 * t237 + t199 * t235 + t221 * t253 + t223 * t251 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
