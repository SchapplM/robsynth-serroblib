% Calculate kinetic energy for
% S6RRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 13:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR1_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:12:22
% EndTime: 2019-03-09 13:12:25
% DurationCPUTime: 2.76s
% Computational Cost: add. (2217->335), mult. (1816->489), div. (0->0), fcn. (1598->12), ass. (0->178)
t324 = Icges(3,3) + Icges(4,3);
t242 = qJ(2) + pkin(11);
t229 = sin(t242);
t230 = cos(t242);
t245 = sin(qJ(2));
t248 = cos(qJ(2));
t323 = Icges(3,5) * t248 + Icges(4,5) * t230 - Icges(3,6) * t245 - Icges(4,6) * t229;
t246 = sin(qJ(1));
t249 = cos(qJ(1));
t308 = Icges(4,4) * t230;
t274 = -Icges(4,2) * t229 + t308;
t158 = -Icges(4,6) * t249 + t246 * t274;
t159 = Icges(4,6) * t246 + t249 * t274;
t309 = Icges(4,4) * t229;
t278 = Icges(4,1) * t230 - t309;
t160 = -Icges(4,5) * t249 + t246 * t278;
t161 = Icges(4,5) * t246 + t249 * t278;
t310 = Icges(3,4) * t248;
t275 = -Icges(3,2) * t245 + t310;
t174 = -Icges(3,6) * t249 + t246 * t275;
t175 = Icges(3,6) * t246 + t249 * t275;
t311 = Icges(3,4) * t245;
t279 = Icges(3,1) * t248 - t311;
t176 = -Icges(3,5) * t249 + t246 * t279;
t177 = Icges(3,5) * t246 + t249 * t279;
t194 = Icges(4,2) * t230 + t309;
t195 = Icges(4,1) * t229 + t308;
t209 = Icges(3,2) * t248 + t311;
t212 = Icges(3,1) * t245 + t310;
t220 = -qJD(2) * t249 + V_base(5);
t221 = qJD(2) * t246 + V_base(4);
t232 = V_base(6) + qJD(1);
t322 = (-t194 * t229 + t195 * t230 - t209 * t245 + t212 * t248) * t232 + (-t159 * t229 + t161 * t230 - t175 * t245 + t177 * t248) * t221 + (-t158 * t229 + t160 * t230 - t174 * t245 + t176 * t248) * t220;
t321 = (Icges(3,5) * t245 + Icges(4,5) * t229 + Icges(3,6) * t248 + Icges(4,6) * t230) * t232 + (t324 * t246 + t323 * t249) * t221 + (t323 * t246 - t324 * t249) * t220;
t317 = pkin(2) * t245;
t316 = pkin(3) * t229;
t231 = qJ(4) + t242;
t225 = sin(t231);
t315 = pkin(4) * t225;
t314 = t248 * pkin(2);
t312 = Icges(2,4) * t246;
t307 = Icges(5,4) * t225;
t226 = cos(t231);
t306 = Icges(5,4) * t226;
t227 = qJ(5) + t231;
t222 = sin(t227);
t305 = Icges(6,4) * t222;
t223 = cos(t227);
t304 = Icges(6,4) * t223;
t303 = t222 * t246;
t302 = t222 * t249;
t244 = sin(qJ(6));
t301 = t244 * t249;
t300 = t246 * t244;
t247 = cos(qJ(6));
t299 = t246 * t247;
t298 = t247 * t249;
t154 = -qJ(3) * t249 + t246 * t314;
t218 = t246 * pkin(1) - pkin(7) * t249;
t297 = -t154 - t218;
t296 = pkin(4) * t226;
t295 = pkin(3) * t230;
t292 = qJD(6) * t222;
t291 = -qJD(2) - qJD(4);
t290 = V_base(5) * pkin(6) + V_base(1);
t125 = -pkin(8) * t249 + t246 * t295;
t287 = -t125 + t297;
t121 = -pkin(9) * t249 + t246 * t296;
t286 = -t121 + t287;
t199 = qJD(4) * t246 + t221;
t285 = qJD(3) * t246 + t220 * t317 + t290;
t284 = pkin(5) * t223 + pkin(10) * t222;
t283 = rSges(3,1) * t248 - rSges(3,2) * t245;
t282 = rSges(4,1) * t230 - rSges(4,2) * t229;
t281 = rSges(5,1) * t226 - rSges(5,2) * t225;
t280 = rSges(6,1) * t223 - rSges(6,2) * t222;
t186 = qJD(5) * t246 + t199;
t277 = Icges(5,1) * t226 - t307;
t276 = Icges(6,1) * t223 - t305;
t273 = -Icges(5,2) * t225 + t306;
t272 = -Icges(6,2) * t222 + t304;
t269 = Icges(5,5) * t226 - Icges(5,6) * t225;
t268 = Icges(6,5) * t223 - Icges(6,6) * t222;
t267 = t220 * t316 + t285;
t219 = pkin(1) * t249 + t246 * pkin(7);
t266 = -V_base(4) * pkin(6) + t232 * t219 + V_base(2);
t265 = V_base(4) * t218 - t219 * V_base(5) + V_base(3);
t198 = t249 * t291 + V_base(5);
t264 = t198 * t315 + t267;
t263 = t221 * t154 + t265;
t185 = V_base(5) + (-qJD(5) + t291) * t249;
t262 = (-Icges(6,3) * t249 + t246 * t268) * t185 + (Icges(6,3) * t246 + t249 * t268) * t186 + (Icges(6,5) * t222 + Icges(6,6) * t223) * t232;
t261 = (-Icges(5,3) * t249 + t246 * t269) * t198 + (Icges(5,3) * t246 + t249 * t269) * t199 + (Icges(5,5) * t225 + Icges(5,6) * t226) * t232;
t155 = qJ(3) * t246 + t249 * t314;
t258 = -qJD(3) * t249 + t232 * t155 + t266;
t126 = pkin(8) * t246 + t249 * t295;
t257 = t221 * t125 + (-t126 - t155) * t220 + t263;
t122 = pkin(9) * t246 + t249 * t296;
t256 = t199 * t121 - t122 * t198 + t257;
t255 = t232 * t126 + (-t316 - t317) * t221 + t258;
t254 = t232 * t122 - t199 * t315 + t255;
t137 = -Icges(6,6) * t249 + t246 * t272;
t138 = Icges(6,6) * t246 + t249 * t272;
t139 = -Icges(6,5) * t249 + t246 * t276;
t140 = Icges(6,5) * t246 + t249 * t276;
t179 = Icges(6,2) * t223 + t305;
t180 = Icges(6,1) * t222 + t304;
t253 = (-t138 * t222 + t140 * t223) * t186 + (-t137 * t222 + t139 * t223) * t185 + (-t179 * t222 + t180 * t223) * t232;
t146 = -Icges(5,6) * t249 + t246 * t273;
t147 = Icges(5,6) * t246 + t249 * t273;
t148 = -Icges(5,5) * t249 + t246 * t277;
t149 = Icges(5,5) * t246 + t249 * t277;
t188 = Icges(5,2) * t226 + t307;
t189 = Icges(5,1) * t225 + t306;
t252 = (-t147 * t225 + t149 * t226) * t199 + (-t146 * t225 + t148 * t226) * t198 + (-t188 * t225 + t189 * t226) * t232;
t238 = Icges(2,4) * t249;
t217 = rSges(2,1) * t249 - t246 * rSges(2,2);
t216 = t246 * rSges(2,1) + rSges(2,2) * t249;
t215 = rSges(3,1) * t245 + rSges(3,2) * t248;
t214 = Icges(2,1) * t249 - t312;
t213 = Icges(2,1) * t246 + t238;
t211 = -Icges(2,2) * t246 + t238;
t210 = Icges(2,2) * t249 + t312;
t205 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t204 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t203 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t197 = rSges(4,1) * t229 + rSges(4,2) * t230;
t196 = -qJD(6) * t223 + t232;
t190 = rSges(5,1) * t225 + rSges(5,2) * t226;
t184 = pkin(5) * t222 - pkin(10) * t223;
t183 = rSges(6,1) * t222 + rSges(6,2) * t223;
t182 = t246 * rSges(3,3) + t249 * t283;
t181 = -rSges(3,3) * t249 + t246 * t283;
t171 = t223 * t298 + t300;
t170 = -t223 * t301 + t299;
t169 = t223 * t299 - t301;
t168 = -t223 * t300 - t298;
t165 = t246 * rSges(4,3) + t249 * t282;
t164 = -rSges(4,3) * t249 + t246 * t282;
t163 = t284 * t249;
t162 = t284 * t246;
t153 = t246 * rSges(5,3) + t249 * t281;
t152 = -rSges(5,3) * t249 + t246 * t281;
t151 = V_base(5) * rSges(2,3) - t216 * t232 + t290;
t150 = t217 * t232 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t143 = t216 * V_base(4) - t217 * V_base(5) + V_base(3);
t142 = t246 * rSges(6,3) + t249 * t280;
t141 = -rSges(6,3) * t249 + t246 * t280;
t134 = t249 * t292 + t186;
t133 = t246 * t292 + t185;
t130 = -rSges(7,3) * t223 + (rSges(7,1) * t247 - rSges(7,2) * t244) * t222;
t129 = -Icges(7,5) * t223 + (Icges(7,1) * t247 - Icges(7,4) * t244) * t222;
t128 = -Icges(7,6) * t223 + (Icges(7,4) * t247 - Icges(7,2) * t244) * t222;
t127 = -Icges(7,3) * t223 + (Icges(7,5) * t247 - Icges(7,6) * t244) * t222;
t120 = t171 * rSges(7,1) + t170 * rSges(7,2) + rSges(7,3) * t302;
t119 = rSges(7,1) * t169 + rSges(7,2) * t168 + rSges(7,3) * t303;
t118 = Icges(7,1) * t171 + Icges(7,4) * t170 + Icges(7,5) * t302;
t117 = Icges(7,1) * t169 + Icges(7,4) * t168 + Icges(7,5) * t303;
t116 = Icges(7,4) * t171 + Icges(7,2) * t170 + Icges(7,6) * t302;
t115 = Icges(7,4) * t169 + Icges(7,2) * t168 + Icges(7,6) * t303;
t114 = Icges(7,5) * t171 + Icges(7,6) * t170 + Icges(7,3) * t302;
t113 = Icges(7,5) * t169 + Icges(7,6) * t168 + Icges(7,3) * t303;
t111 = t215 * t220 + (-t181 - t218) * t232 + t290;
t110 = t182 * t232 - t215 * t221 + t266;
t108 = t181 * t221 - t182 * t220 + t265;
t107 = t197 * t220 + (-t164 + t297) * t232 + t285;
t106 = t232 * t165 + (-t197 - t317) * t221 + t258;
t105 = t164 * t221 + (-t155 - t165) * t220 + t263;
t104 = t190 * t198 + (-t152 + t287) * t232 + t267;
t103 = t232 * t153 - t199 * t190 + t255;
t102 = t152 * t199 - t153 * t198 + t257;
t101 = t183 * t185 + (-t141 + t286) * t232 + t264;
t100 = t232 * t142 - t186 * t183 + t254;
t99 = -t119 * t196 + t130 * t133 + t184 * t185 + (-t162 + t286) * t232 + t264;
t98 = t196 * t120 - t134 * t130 + t232 * t163 - t186 * t184 + t254;
t97 = t141 * t186 - t142 * t185 + t256;
t96 = t119 * t134 - t120 * t133 + t162 * t186 - t163 * t185 + t256;
t1 = t134 * ((t114 * t302 + t170 * t116 + t171 * t118) * t134 + (t113 * t302 + t170 * t115 + t171 * t117) * t133 + (t127 * t302 + t170 * t128 + t171 * t129) * t196) / 0.2e1 + t133 * ((t114 * t303 + t116 * t168 + t118 * t169) * t134 + (t113 * t303 + t168 * t115 + t169 * t117) * t133 + (t127 * t303 + t128 * t168 + t129 * t169) * t196) / 0.2e1 + t199 * (t246 * t261 + t249 * t252) / 0.2e1 + t198 * (t246 * t252 - t249 * t261) / 0.2e1 + t186 * (t246 * t262 + t249 * t253) / 0.2e1 + t185 * (t246 * t253 - t249 * t262) / 0.2e1 + m(2) * (t143 ^ 2 + t150 ^ 2 + t151 ^ 2) / 0.2e1 + t196 * ((-t113 * t133 - t114 * t134 - t127 * t196) * t223 + ((-t116 * t244 + t118 * t247) * t134 + (-t115 * t244 + t117 * t247) * t133 + (-t128 * t244 + t129 * t247) * t196) * t222) / 0.2e1 + m(4) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(3) * (t108 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(7) * (t96 ^ 2 + t98 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t100 ^ 2 + t101 ^ 2 + t97 ^ 2) / 0.2e1 + m(5) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + m(1) * (t203 ^ 2 + t204 ^ 2 + t205 ^ 2) / 0.2e1 + (t322 * t246 - t321 * t249) * t220 / 0.2e1 + (t321 * t246 + t322 * t249) * t221 / 0.2e1 + ((-t246 * t210 + t213 * t249 + Icges(1,4)) * V_base(5) + (-t246 * t211 + t249 * t214 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t249 * t210 + t246 * t213 + Icges(1,2)) * V_base(5) + (t211 * t249 + t246 * t214 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t147 * t226 + t149 * t225) * t199 + (t146 * t226 + t148 * t225) * t198 + (t138 * t223 + t140 * t222) * t186 + (t137 * t223 + t139 * t222) * t185 + (t159 * t230 + t161 * t229 + t175 * t248 + t177 * t245) * t221 + (t158 * t230 + t160 * t229 + t174 * t248 + t176 * t245) * t220 + (t223 * t179 + t222 * t180 + t226 * t188 + t225 * t189 + t230 * t194 + t229 * t195 + t248 * t209 + t245 * t212 + Icges(2,3)) * t232) * t232 / 0.2e1 + t232 * V_base(4) * (Icges(2,5) * t249 - Icges(2,6) * t246) + V_base(5) * t232 * (Icges(2,5) * t246 + Icges(2,6) * t249) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
