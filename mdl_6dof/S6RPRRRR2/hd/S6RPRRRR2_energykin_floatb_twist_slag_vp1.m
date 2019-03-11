% Calculate kinetic energy for
% S6RPRRRR2
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
% Datum: 2019-03-09 06:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR2_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRRR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:56:56
% EndTime: 2019-03-09 06:56:59
% DurationCPUTime: 2.58s
% Computational Cost: add. (2372->346), mult. (1800->511), div. (0->0), fcn. (1636->12), ass. (0->172)
t246 = sin(qJ(1));
t305 = pkin(1) * t246;
t249 = cos(qJ(1));
t304 = pkin(1) * t249;
t245 = sin(qJ(3));
t303 = pkin(3) * t245;
t248 = cos(qJ(3));
t302 = pkin(3) * t248;
t247 = cos(qJ(5));
t301 = pkin(5) * t247;
t300 = -pkin(6) - qJ(2);
t297 = Icges(2,4) * t246;
t241 = qJ(1) + pkin(11);
t231 = sin(t241);
t296 = Icges(3,4) * t231;
t295 = Icges(4,4) * t245;
t294 = Icges(4,4) * t248;
t243 = qJ(3) + qJ(4);
t235 = sin(t243);
t293 = Icges(5,4) * t235;
t237 = cos(t243);
t292 = Icges(5,4) * t237;
t291 = t231 * t235;
t244 = sin(qJ(5));
t290 = t231 * t244;
t232 = cos(t241);
t289 = t232 * t235;
t288 = t232 * t244;
t242 = qJ(5) + qJ(6);
t234 = sin(t242);
t287 = t234 * t237;
t236 = cos(t242);
t286 = t236 * t237;
t285 = t237 * t244;
t284 = t237 * t247;
t283 = qJD(5) * t235;
t282 = qJD(6) * t235;
t233 = V_base(6) + qJD(1);
t281 = t233 * t304 + V_base(2);
t280 = V_base(5) * pkin(6) + V_base(1);
t211 = qJD(3) * t231 + V_base(4);
t197 = t231 * pkin(2) - t232 * pkin(7);
t277 = -t197 - t305;
t276 = V_base(5) * qJ(2) + t280;
t275 = V_base(4) * t305 + qJD(2) + V_base(3);
t186 = qJD(4) * t231 + t211;
t139 = -pkin(8) * t232 + t231 * t302;
t274 = -t139 + t277;
t210 = -qJD(3) * t232 + V_base(5);
t273 = t210 * t303 + t276;
t272 = pkin(4) * t237 + pkin(9) * t235;
t271 = rSges(4,1) * t248 - rSges(4,2) * t245;
t270 = rSges(5,1) * t237 - rSges(5,2) * t235;
t155 = t232 * t283 + t186;
t269 = Icges(4,1) * t248 - t295;
t268 = Icges(5,1) * t237 - t293;
t267 = -Icges(4,2) * t245 + t294;
t266 = -Icges(5,2) * t235 + t292;
t265 = Icges(4,5) * t248 - Icges(4,6) * t245;
t264 = Icges(5,5) * t237 - Icges(5,6) * t235;
t185 = V_base(5) + (-qJD(3) - qJD(4)) * t232;
t154 = t231 * t283 + t185;
t263 = pkin(10) * t235 + t237 * t301;
t262 = (-Icges(5,3) * t232 + t231 * t264) * t185 + (Icges(5,3) * t231 + t232 * t264) * t186 + (Icges(5,5) * t235 + Icges(5,6) * t237) * t233;
t261 = (-Icges(4,3) * t232 + t231 * t265) * t210 + (Icges(4,3) * t231 + t232 * t265) * t211 + (Icges(4,5) * t245 + Icges(4,6) * t248) * t233;
t198 = t232 * pkin(2) + t231 * pkin(7);
t260 = t233 * t198 + t300 * V_base(4) + t281;
t259 = V_base(4) * t197 + (-t198 - t304) * V_base(5) + t275;
t140 = pkin(8) * t231 + t232 * t302;
t258 = t233 * t140 - t211 * t303 + t260;
t177 = t272 * t231;
t203 = pkin(4) * t235 - pkin(9) * t237;
t257 = t185 * t203 + (-t177 + t274) * t233 + t273;
t256 = t211 * t139 - t140 * t210 + t259;
t178 = t272 * t232;
t255 = t233 * t178 - t186 * t203 + t258;
t254 = t186 * t177 - t178 * t185 + t256;
t145 = -Icges(5,6) * t232 + t231 * t266;
t146 = Icges(5,6) * t231 + t232 * t266;
t147 = -Icges(5,5) * t232 + t231 * t268;
t148 = Icges(5,5) * t231 + t232 * t268;
t200 = Icges(5,2) * t237 + t293;
t201 = Icges(5,1) * t235 + t292;
t253 = (-t146 * t235 + t148 * t237) * t186 + (-t145 * t235 + t147 * t237) * t185 + (-t200 * t235 + t201 * t237) * t233;
t161 = -Icges(4,6) * t232 + t231 * t267;
t162 = Icges(4,6) * t231 + t232 * t267;
t163 = -Icges(4,5) * t232 + t231 * t269;
t164 = Icges(4,5) * t231 + t232 * t269;
t215 = Icges(4,2) * t248 + t295;
t218 = Icges(4,1) * t245 + t294;
t252 = (-t162 * t245 + t164 * t248) * t211 + (-t161 * t245 + t163 * t248) * t210 + (-t215 * t245 + t218 * t248) * t233;
t239 = Icges(2,4) * t249;
t228 = Icges(3,4) * t232;
t223 = rSges(2,1) * t249 - rSges(2,2) * t246;
t222 = rSges(2,1) * t246 + rSges(2,2) * t249;
t221 = rSges(4,1) * t245 + rSges(4,2) * t248;
t220 = Icges(2,1) * t249 - t297;
t219 = Icges(2,1) * t246 + t239;
t217 = -Icges(2,2) * t246 + t239;
t216 = Icges(2,2) * t249 + t297;
t209 = -qJD(5) * t237 + t233;
t208 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t207 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t206 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t202 = rSges(5,1) * t235 + rSges(5,2) * t237;
t196 = rSges(3,1) * t232 - rSges(3,2) * t231;
t195 = rSges(3,1) * t231 + rSges(3,2) * t232;
t194 = Icges(3,1) * t232 - t296;
t193 = Icges(3,1) * t231 + t228;
t192 = -Icges(3,2) * t231 + t228;
t191 = Icges(3,2) * t232 + t296;
t184 = (-qJD(5) - qJD(6)) * t237 + t233;
t182 = t232 * t284 + t290;
t181 = t231 * t247 - t232 * t285;
t180 = t231 * t284 - t288;
t179 = -t231 * t285 - t232 * t247;
t176 = -rSges(6,3) * t237 + (rSges(6,1) * t247 - rSges(6,2) * t244) * t235;
t175 = -Icges(6,5) * t237 + (Icges(6,1) * t247 - Icges(6,4) * t244) * t235;
t174 = -Icges(6,6) * t237 + (Icges(6,4) * t247 - Icges(6,2) * t244) * t235;
t173 = -Icges(6,3) * t237 + (Icges(6,5) * t247 - Icges(6,6) * t244) * t235;
t172 = t231 * t234 + t232 * t286;
t171 = t231 * t236 - t232 * t287;
t170 = t231 * t286 - t232 * t234;
t169 = -t231 * t287 - t232 * t236;
t168 = rSges(4,3) * t231 + t232 * t271;
t167 = -rSges(4,3) * t232 + t231 * t271;
t166 = V_base(5) * rSges(2,3) - t222 * t233 + t280;
t165 = t223 * t233 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t157 = t222 * V_base(4) - t223 * V_base(5) + V_base(3);
t156 = -rSges(7,3) * t237 + (rSges(7,1) * t236 - rSges(7,2) * t234) * t235;
t153 = -Icges(7,5) * t237 + (Icges(7,1) * t236 - Icges(7,4) * t234) * t235;
t152 = -Icges(7,6) * t237 + (Icges(7,4) * t236 - Icges(7,2) * t234) * t235;
t151 = -Icges(7,3) * t237 + (Icges(7,5) * t236 - Icges(7,6) * t234) * t235;
t150 = rSges(5,3) * t231 + t232 * t270;
t149 = -rSges(5,3) * t232 + t231 * t270;
t141 = -pkin(10) * t237 + t235 * t301;
t137 = t232 * t282 + t155;
t136 = t231 * t282 + t154;
t133 = V_base(5) * rSges(3,3) + (-t195 - t305) * t233 + t276;
t132 = t196 * t233 + (-rSges(3,3) + t300) * V_base(4) + t281;
t131 = t195 * V_base(4) + (-t196 - t304) * V_base(5) + t275;
t130 = rSges(6,1) * t182 + rSges(6,2) * t181 + rSges(6,3) * t289;
t129 = rSges(6,1) * t180 + rSges(6,2) * t179 + rSges(6,3) * t291;
t128 = Icges(6,1) * t182 + Icges(6,4) * t181 + Icges(6,5) * t289;
t127 = Icges(6,1) * t180 + Icges(6,4) * t179 + Icges(6,5) * t291;
t126 = Icges(6,4) * t182 + Icges(6,2) * t181 + Icges(6,6) * t289;
t125 = Icges(6,4) * t180 + Icges(6,2) * t179 + Icges(6,6) * t291;
t124 = Icges(6,5) * t182 + Icges(6,6) * t181 + Icges(6,3) * t289;
t123 = Icges(6,5) * t180 + Icges(6,6) * t179 + Icges(6,3) * t291;
t122 = pkin(5) * t290 + t232 * t263;
t121 = -pkin(5) * t288 + t231 * t263;
t120 = rSges(7,1) * t172 + rSges(7,2) * t171 + rSges(7,3) * t289;
t119 = rSges(7,1) * t170 + rSges(7,2) * t169 + rSges(7,3) * t291;
t118 = Icges(7,1) * t172 + Icges(7,4) * t171 + Icges(7,5) * t289;
t117 = Icges(7,1) * t170 + Icges(7,4) * t169 + Icges(7,5) * t291;
t116 = Icges(7,4) * t172 + Icges(7,2) * t171 + Icges(7,6) * t289;
t115 = Icges(7,4) * t170 + Icges(7,2) * t169 + Icges(7,6) * t291;
t114 = Icges(7,5) * t172 + Icges(7,6) * t171 + Icges(7,3) * t289;
t113 = Icges(7,5) * t170 + Icges(7,6) * t169 + Icges(7,3) * t291;
t112 = t210 * t221 + (-t167 + t277) * t233 + t276;
t111 = t168 * t233 - t211 * t221 + t260;
t110 = t167 * t211 - t168 * t210 + t259;
t109 = t185 * t202 + (-t149 + t274) * t233 + t273;
t108 = t150 * t233 - t186 * t202 + t258;
t107 = t149 * t186 - t150 * t185 + t256;
t106 = -t129 * t209 + t154 * t176 + t257;
t105 = t130 * t209 - t155 * t176 + t255;
t104 = t129 * t155 - t130 * t154 + t254;
t103 = -t119 * t184 - t121 * t209 + t136 * t156 + t141 * t154 + t257;
t102 = t120 * t184 + t122 * t209 - t137 * t156 - t141 * t155 + t255;
t101 = t119 * t137 - t120 * t136 + t121 * t155 - t122 * t154 + t254;
t1 = t154 * ((t124 * t291 + t126 * t179 + t128 * t180) * t155 + (t123 * t291 + t179 * t125 + t180 * t127) * t154 + (t173 * t291 + t174 * t179 + t175 * t180) * t209) / 0.2e1 + t136 * ((t114 * t291 + t116 * t169 + t118 * t170) * t137 + (t113 * t291 + t169 * t115 + t170 * t117) * t136 + (t151 * t291 + t152 * t169 + t153 * t170) * t184) / 0.2e1 + t155 * ((t124 * t289 + t181 * t126 + t182 * t128) * t155 + (t123 * t289 + t125 * t181 + t127 * t182) * t154 + (t173 * t289 + t174 * t181 + t175 * t182) * t209) / 0.2e1 + t137 * ((t114 * t289 + t171 * t116 + t172 * t118) * t137 + (t113 * t289 + t115 * t171 + t117 * t172) * t136 + (t151 * t289 + t152 * t171 + t153 * t172) * t184) / 0.2e1 + m(1) * (t206 ^ 2 + t207 ^ 2 + t208 ^ 2) / 0.2e1 + m(2) * (t157 ^ 2 + t165 ^ 2 + t166 ^ 2) / 0.2e1 + m(3) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(6) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(5) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(4) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(7) * (t101 ^ 2 + t102 ^ 2 + t103 ^ 2) / 0.2e1 + t211 * (t261 * t231 + t252 * t232) / 0.2e1 + t210 * (t252 * t231 - t261 * t232) / 0.2e1 + t186 * (t262 * t231 + t253 * t232) / 0.2e1 + t185 * (t231 * t253 - t232 * t262) / 0.2e1 + t209 * ((-t123 * t154 - t124 * t155 - t173 * t209) * t237 + ((-t126 * t244 + t128 * t247) * t155 + (-t125 * t244 + t127 * t247) * t154 + (-t174 * t244 + t175 * t247) * t209) * t235) / 0.2e1 + t184 * ((-t113 * t136 - t114 * t137 - t151 * t184) * t237 + ((-t116 * t234 + t118 * t236) * t137 + (-t115 * t234 + t117 * t236) * t136 + (-t152 * t234 + t153 * t236) * t184) * t235) / 0.2e1 + ((-t191 * t231 + t193 * t232 - t216 * t246 + t219 * t249 + Icges(1,4)) * V_base(5) + (-t231 * t192 + t232 * t194 - t246 * t217 + t249 * t220 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t232 * t191 + t231 * t193 + t249 * t216 + t246 * t219 + Icges(1,2)) * V_base(5) + (t192 * t232 + t194 * t231 + t217 * t249 + t220 * t246 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t162 * t248 + t164 * t245) * t211 + (t161 * t248 + t163 * t245) * t210 + (t146 * t237 + t148 * t235) * t186 + (t145 * t237 + t147 * t235) * t185 + (t237 * t200 + t235 * t201 + t248 * t215 + t245 * t218 + Icges(2,3) + Icges(3,3)) * t233) * t233 / 0.2e1 + t233 * V_base(5) * (Icges(2,5) * t246 + Icges(3,5) * t231 + Icges(2,6) * t249 + Icges(3,6) * t232) + t233 * V_base(4) * (Icges(2,5) * t249 + Icges(3,5) * t232 - Icges(2,6) * t246 - Icges(3,6) * t231) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
