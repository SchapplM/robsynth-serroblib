% Calculate kinetic energy for
% S6RRPRRR8
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
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR8_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR8_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR8_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR8_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:02:13
% EndTime: 2019-03-09 14:02:17
% DurationCPUTime: 4.00s
% Computational Cost: add. (2448->412), mult. (2746->624), div. (0->0), fcn. (2696->12), ass. (0->185)
t278 = cos(pkin(11));
t330 = t278 * pkin(3);
t281 = sin(qJ(1));
t329 = Icges(2,4) * t281;
t280 = sin(qJ(2));
t328 = Icges(3,4) * t280;
t282 = cos(qJ(2));
t327 = Icges(3,4) * t282;
t277 = sin(pkin(11));
t283 = cos(qJ(1));
t326 = t277 * t283;
t325 = t280 * t281;
t324 = t280 * t283;
t323 = t281 * t277;
t322 = t281 * t282;
t321 = t282 * t283;
t276 = pkin(11) + qJ(4);
t267 = qJ(5) + t276;
t262 = cos(t267);
t319 = pkin(5) * t262;
t261 = sin(t267);
t318 = pkin(5) * t261;
t302 = pkin(2) * t282 + qJ(3) * t280;
t223 = t302 * t281;
t247 = t281 * pkin(1) - pkin(7) * t283;
t317 = -t223 - t247;
t266 = cos(t276);
t316 = pkin(4) * t266;
t313 = qJD(3) * t280;
t312 = qJD(4) * t280;
t311 = qJD(5) * t280;
t310 = qJD(6) * t280;
t309 = -qJD(4) - qJD(5);
t308 = V_base(5) * pkin(6) + V_base(1);
t250 = qJD(2) * t281 + V_base(4);
t268 = V_base(6) + qJD(1);
t265 = sin(t276);
t305 = pkin(4) * t265;
t222 = t283 * t312 + t250;
t243 = pkin(2) * t280 - qJ(3) * t282;
t249 = -qJD(2) * t283 + V_base(5);
t304 = t249 * t243 + t283 * t313 + t308;
t303 = rSges(3,1) * t282 - rSges(3,2) * t280;
t191 = t283 * t311 + t222;
t301 = Icges(3,1) * t282 - t328;
t300 = -Icges(3,2) * t280 + t327;
t299 = Icges(3,5) * t282 - Icges(3,6) * t280;
t221 = t281 * t312 + t249;
t248 = pkin(1) * t283 + t281 * pkin(7);
t298 = -V_base(4) * pkin(6) + t268 * t248 + V_base(2);
t297 = V_base(4) * t247 - t248 * V_base(5) + V_base(3);
t190 = t281 * t311 + t221;
t296 = (-Icges(3,3) * t283 + t281 * t299) * t249 + (Icges(3,3) * t281 + t283 * t299) * t250 + (Icges(3,5) * t280 + Icges(3,6) * t282) * t268;
t295 = pkin(8) * t280 + t282 * t330;
t224 = t302 * t283;
t294 = t268 * t224 + t281 * t313 + t298;
t293 = pkin(10) * t280 + t282 * t319;
t292 = pkin(9) * t280 + t282 * t316;
t291 = -qJD(3) * t282 + t250 * t223 + t297;
t163 = -pkin(3) * t326 + t281 * t295;
t180 = -pkin(8) * t282 + t280 * t330;
t290 = t249 * t180 + (-t163 + t317) * t268 + t304;
t164 = pkin(3) * t323 + t283 * t295;
t289 = t268 * t164 + (-t180 - t243) * t250 + t294;
t123 = t281 * t292 - t283 * t305;
t165 = -pkin(9) * t282 + t280 * t316;
t242 = -qJD(4) * t282 + t268;
t288 = -t123 * t242 + t221 * t165 + t290;
t287 = t250 * t163 + (-t164 - t224) * t249 + t291;
t124 = t281 * t305 + t283 * t292;
t286 = t242 * t124 - t165 * t222 + t289;
t285 = t222 * t123 - t124 * t221 + t287;
t208 = -Icges(3,6) * t283 + t281 * t300;
t209 = Icges(3,6) * t281 + t283 * t300;
t210 = -Icges(3,5) * t283 + t281 * t301;
t211 = Icges(3,5) * t281 + t283 * t301;
t236 = Icges(3,2) * t282 + t328;
t239 = Icges(3,1) * t280 + t327;
t284 = (-t209 * t280 + t211 * t282) * t250 + (-t208 * t280 + t210 * t282) * t249 + (-t236 * t280 + t239 * t282) * t268;
t273 = Icges(2,4) * t283;
t264 = qJ(6) + t267;
t252 = cos(t264);
t251 = sin(t264);
t246 = rSges(2,1) * t283 - t281 * rSges(2,2);
t245 = t281 * rSges(2,1) + rSges(2,2) * t283;
t244 = rSges(3,1) * t280 + rSges(3,2) * t282;
t241 = Icges(2,1) * t283 - t329;
t240 = Icges(2,1) * t281 + t273;
t238 = -Icges(2,2) * t281 + t273;
t237 = Icges(2,2) * t283 + t329;
t232 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t231 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t230 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t226 = t282 * t309 + t268;
t220 = t278 * t321 + t323;
t219 = -t277 * t321 + t281 * t278;
t218 = t278 * t322 - t326;
t217 = -t277 * t322 - t278 * t283;
t215 = (-qJD(6) + t309) * t282 + t268;
t214 = t281 * rSges(3,3) + t283 * t303;
t213 = -rSges(3,3) * t283 + t281 * t303;
t205 = t281 * t265 + t266 * t321;
t204 = -t265 * t321 + t281 * t266;
t203 = -t265 * t283 + t266 * t322;
t202 = -t265 * t322 - t266 * t283;
t201 = -rSges(4,3) * t282 + (rSges(4,1) * t278 - rSges(4,2) * t277) * t280;
t199 = -Icges(4,5) * t282 + (Icges(4,1) * t278 - Icges(4,4) * t277) * t280;
t198 = -Icges(4,6) * t282 + (Icges(4,4) * t278 - Icges(4,2) * t277) * t280;
t197 = -Icges(4,3) * t282 + (Icges(4,5) * t278 - Icges(4,6) * t277) * t280;
t195 = t281 * t261 + t262 * t321;
t194 = -t261 * t321 + t281 * t262;
t193 = -t261 * t283 + t262 * t322;
t192 = -t261 * t322 - t262 * t283;
t188 = -rSges(5,3) * t282 + (rSges(5,1) * t266 - rSges(5,2) * t265) * t280;
t187 = -Icges(5,5) * t282 + (Icges(5,1) * t266 - Icges(5,4) * t265) * t280;
t186 = -Icges(5,6) * t282 + (Icges(5,4) * t266 - Icges(5,2) * t265) * t280;
t185 = -Icges(5,3) * t282 + (Icges(5,5) * t266 - Icges(5,6) * t265) * t280;
t184 = t281 * t251 + t252 * t321;
t183 = -t251 * t321 + t281 * t252;
t182 = -t251 * t283 + t252 * t322;
t181 = -t251 * t322 - t252 * t283;
t179 = -rSges(6,3) * t282 + (rSges(6,1) * t262 - rSges(6,2) * t261) * t280;
t178 = V_base(5) * rSges(2,3) - t245 * t268 + t308;
t177 = t246 * t268 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t176 = -Icges(6,5) * t282 + (Icges(6,1) * t262 - Icges(6,4) * t261) * t280;
t175 = -Icges(6,6) * t282 + (Icges(6,4) * t262 - Icges(6,2) * t261) * t280;
t174 = -Icges(6,3) * t282 + (Icges(6,5) * t262 - Icges(6,6) * t261) * t280;
t173 = t245 * V_base(4) - t246 * V_base(5) + V_base(3);
t172 = -rSges(7,3) * t282 + (rSges(7,1) * t252 - rSges(7,2) * t251) * t280;
t171 = -Icges(7,5) * t282 + (Icges(7,1) * t252 - Icges(7,4) * t251) * t280;
t170 = -Icges(7,6) * t282 + (Icges(7,4) * t252 - Icges(7,2) * t251) * t280;
t169 = -Icges(7,3) * t282 + (Icges(7,5) * t252 - Icges(7,6) * t251) * t280;
t168 = t283 * t310 + t191;
t167 = t281 * t310 + t190;
t162 = t220 * rSges(4,1) + t219 * rSges(4,2) + rSges(4,3) * t324;
t161 = rSges(4,1) * t218 + rSges(4,2) * t217 + rSges(4,3) * t325;
t160 = Icges(4,1) * t220 + Icges(4,4) * t219 + Icges(4,5) * t324;
t159 = Icges(4,1) * t218 + Icges(4,4) * t217 + Icges(4,5) * t325;
t158 = Icges(4,4) * t220 + Icges(4,2) * t219 + Icges(4,6) * t324;
t157 = Icges(4,4) * t218 + Icges(4,2) * t217 + Icges(4,6) * t325;
t156 = Icges(4,5) * t220 + Icges(4,6) * t219 + Icges(4,3) * t324;
t155 = Icges(4,5) * t218 + Icges(4,6) * t217 + Icges(4,3) * t325;
t153 = t205 * rSges(5,1) + t204 * rSges(5,2) + rSges(5,3) * t324;
t152 = rSges(5,1) * t203 + rSges(5,2) * t202 + rSges(5,3) * t325;
t151 = Icges(5,1) * t205 + Icges(5,4) * t204 + Icges(5,5) * t324;
t150 = Icges(5,1) * t203 + Icges(5,4) * t202 + Icges(5,5) * t325;
t149 = Icges(5,4) * t205 + Icges(5,2) * t204 + Icges(5,6) * t324;
t148 = Icges(5,4) * t203 + Icges(5,2) * t202 + Icges(5,6) * t325;
t147 = Icges(5,5) * t205 + Icges(5,6) * t204 + Icges(5,3) * t324;
t146 = Icges(5,5) * t203 + Icges(5,6) * t202 + Icges(5,3) * t325;
t144 = -pkin(10) * t282 + t280 * t319;
t142 = t195 * rSges(6,1) + t194 * rSges(6,2) + rSges(6,3) * t324;
t141 = rSges(6,1) * t193 + rSges(6,2) * t192 + rSges(6,3) * t325;
t140 = Icges(6,1) * t195 + Icges(6,4) * t194 + Icges(6,5) * t324;
t139 = Icges(6,1) * t193 + Icges(6,4) * t192 + Icges(6,5) * t325;
t138 = Icges(6,4) * t195 + Icges(6,2) * t194 + Icges(6,6) * t324;
t137 = Icges(6,4) * t193 + Icges(6,2) * t192 + Icges(6,6) * t325;
t136 = Icges(6,5) * t195 + Icges(6,6) * t194 + Icges(6,3) * t324;
t135 = Icges(6,5) * t193 + Icges(6,6) * t192 + Icges(6,3) * t325;
t134 = t184 * rSges(7,1) + t183 * rSges(7,2) + rSges(7,3) * t324;
t133 = rSges(7,1) * t182 + rSges(7,2) * t181 + rSges(7,3) * t325;
t132 = Icges(7,1) * t184 + Icges(7,4) * t183 + Icges(7,5) * t324;
t131 = Icges(7,1) * t182 + Icges(7,4) * t181 + Icges(7,5) * t325;
t130 = Icges(7,4) * t184 + Icges(7,2) * t183 + Icges(7,6) * t324;
t129 = Icges(7,4) * t182 + Icges(7,2) * t181 + Icges(7,6) * t325;
t128 = Icges(7,5) * t184 + Icges(7,6) * t183 + Icges(7,3) * t324;
t127 = Icges(7,5) * t182 + Icges(7,6) * t181 + Icges(7,3) * t325;
t126 = t244 * t249 + (-t213 - t247) * t268 + t308;
t125 = t214 * t268 - t244 * t250 + t298;
t121 = t213 * t250 - t214 * t249 + t297;
t119 = t281 * t318 + t283 * t293;
t118 = t281 * t293 - t283 * t318;
t117 = t201 * t249 + (-t161 + t317) * t268 + t304;
t116 = t162 * t268 + (-t201 - t243) * t250 + t294;
t115 = t161 * t250 + (-t162 - t224) * t249 + t291;
t114 = -t152 * t242 + t188 * t221 + t290;
t113 = t153 * t242 - t188 * t222 + t289;
t112 = t152 * t222 - t153 * t221 + t287;
t111 = -t141 * t226 + t179 * t190 + t288;
t110 = t142 * t226 - t179 * t191 + t286;
t109 = t141 * t191 - t142 * t190 + t285;
t108 = -t118 * t226 - t133 * t215 + t144 * t190 + t167 * t172 + t288;
t107 = t119 * t226 + t134 * t215 - t144 * t191 - t168 * t172 + t286;
t106 = t118 * t191 - t119 * t190 + t133 * t168 - t134 * t167 + t285;
t1 = m(1) * (t230 ^ 2 + t231 ^ 2 + t232 ^ 2) / 0.2e1 + m(2) * (t173 ^ 2 + t177 ^ 2 + t178 ^ 2) / 0.2e1 + m(3) * (t121 ^ 2 + t125 ^ 2 + t126 ^ 2) / 0.2e1 + m(4) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(6) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + m(5) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(7) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + t215 * ((-t127 * t167 - t128 * t168 - t169 * t215) * t282 + ((-t130 * t251 + t132 * t252) * t168 + (-t129 * t251 + t131 * t252) * t167 + (-t170 * t251 + t171 * t252) * t215) * t280) / 0.2e1 + t242 * ((-t146 * t221 - t147 * t222 - t185 * t242) * t282 + ((-t149 * t265 + t151 * t266) * t222 + (-t148 * t265 + t150 * t266) * t221 + (-t186 * t265 + t187 * t266) * t242) * t280) / 0.2e1 + t222 * ((t147 * t324 + t204 * t149 + t205 * t151) * t222 + (t146 * t324 + t204 * t148 + t205 * t150) * t221 + (t185 * t324 + t204 * t186 + t205 * t187) * t242) / 0.2e1 + t191 * ((t136 * t324 + t194 * t138 + t195 * t140) * t191 + (t135 * t324 + t194 * t137 + t195 * t139) * t190 + (t174 * t324 + t194 * t175 + t195 * t176) * t226) / 0.2e1 + t168 * ((t128 * t324 + t183 * t130 + t184 * t132) * t168 + (t127 * t324 + t183 * t129 + t184 * t131) * t167 + (t169 * t324 + t183 * t170 + t184 * t171) * t215) / 0.2e1 + t221 * ((t147 * t325 + t149 * t202 + t151 * t203) * t222 + (t146 * t325 + t148 * t202 + t150 * t203) * t221 + (t185 * t325 + t186 * t202 + t187 * t203) * t242) / 0.2e1 + t190 * ((t136 * t325 + t138 * t192 + t140 * t193) * t191 + (t135 * t325 + t137 * t192 + t139 * t193) * t190 + (t174 * t325 + t175 * t192 + t176 * t193) * t226) / 0.2e1 + t167 * ((t128 * t325 + t130 * t181 + t132 * t182) * t168 + (t127 * t325 + t181 * t129 + t182 * t131) * t167 + (t169 * t325 + t170 * t181 + t171 * t182) * t215) / 0.2e1 + t226 * ((-t135 * t190 - t136 * t191 - t174 * t226) * t282 + ((-t138 * t261 + t140 * t262) * t191 + (-t137 * t261 + t139 * t262) * t190 + (-t175 * t261 + t176 * t262) * t226) * t280) / 0.2e1 + ((t156 * t325 + t158 * t217 + t160 * t218) * t250 + (t155 * t325 + t157 * t217 + t159 * t218) * t249 + (t197 * t325 + t198 * t217 + t199 * t218) * t268 + t284 * t281 - t296 * t283) * t249 / 0.2e1 + ((t156 * t324 + t219 * t158 + t220 * t160) * t250 + (t155 * t324 + t219 * t157 + t220 * t159) * t249 + (t197 * t324 + t219 * t198 + t220 * t199) * t268 + t296 * t281 + t284 * t283) * t250 / 0.2e1 + ((-t281 * t237 + t240 * t283 + Icges(1,4)) * V_base(5) + (-t281 * t238 + t241 * t283 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t237 * t283 + t281 * t240 + Icges(1,2)) * V_base(5) + (t238 * t283 + t281 * t241 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t155 * t249 - t156 * t250) * t282 + ((-t158 * t277 + t160 * t278) * t250 + (-t157 * t277 + t159 * t278) * t249) * t280 + (t209 * t282 + t211 * t280) * t250 + (t208 * t282 + t210 * t280) * t249 + (Icges(2,3) + (-t197 + t236) * t282 + (-t198 * t277 + t199 * t278 + t239) * t280) * t268) * t268 / 0.2e1 + t268 * V_base(4) * (Icges(2,5) * t283 - Icges(2,6) * t281) + t268 * V_base(5) * (Icges(2,5) * t281 + Icges(2,6) * t283) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
