% Calculate kinetic energy for
% S6RPRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR6_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR6_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:15:15
% EndTime: 2019-03-09 05:15:19
% DurationCPUTime: 4.12s
% Computational Cost: add. (2373->383), mult. (2292->552), div. (0->0), fcn. (2186->12), ass. (0->186)
t335 = -Icges(6,3) - Icges(5,3);
t258 = pkin(10) + qJ(3);
t248 = cos(t258);
t259 = qJ(4) + pkin(11);
t247 = sin(t259);
t265 = sin(qJ(1));
t310 = t265 * t247;
t249 = cos(t259);
t267 = cos(qJ(1));
t314 = t249 * t267;
t187 = -t248 * t310 - t314;
t309 = t265 * t249;
t188 = -t247 * t267 + t248 * t309;
t266 = cos(qJ(4));
t306 = t266 * t267;
t264 = sin(qJ(4));
t308 = t265 * t264;
t204 = -t248 * t308 - t306;
t307 = t265 * t266;
t313 = t264 * t267;
t205 = t248 * t307 - t313;
t246 = sin(t258);
t317 = t246 * t265;
t334 = Icges(5,5) * t205 + Icges(6,5) * t188 + Icges(5,6) * t204 + Icges(6,6) * t187 - t335 * t317;
t315 = t248 * t267;
t189 = -t247 * t315 + t309;
t190 = t248 * t314 + t310;
t206 = -t248 * t313 + t307;
t207 = t248 * t306 + t308;
t316 = t246 * t267;
t333 = Icges(5,5) * t207 + Icges(6,5) * t190 + Icges(5,6) * t206 + Icges(6,6) * t189 - t335 * t316;
t332 = t335 * t248 + (Icges(5,5) * t266 + Icges(6,5) * t249 - Icges(5,6) * t264 - Icges(6,6) * t247) * t246;
t260 = sin(pkin(10));
t327 = pkin(2) * t260;
t261 = cos(pkin(10));
t325 = pkin(2) * t261;
t324 = t266 * pkin(4);
t322 = Icges(2,4) * t265;
t321 = Icges(3,4) * t260;
t320 = Icges(3,4) * t261;
t319 = Icges(4,4) * t246;
t318 = Icges(4,4) * t248;
t250 = qJ(6) + t259;
t242 = sin(t250);
t312 = t265 * t242;
t243 = cos(t250);
t311 = t265 * t243;
t171 = -pkin(7) * t267 + t265 * t325;
t235 = t265 * pkin(1) - qJ(2) * t267;
t304 = -t171 - t235;
t303 = pkin(5) * t249;
t301 = qJD(4) * t246;
t300 = qJD(5) * t246;
t299 = qJD(6) * t246;
t298 = V_base(4) * t235 + V_base(3);
t297 = V_base(5) * pkin(6) + V_base(1);
t240 = qJD(3) * t265 + V_base(4);
t251 = V_base(6) + qJD(1);
t294 = pkin(5) * t247;
t293 = qJD(2) * t265 + t297;
t203 = t267 * t301 + t240;
t292 = V_base(5) * t327 + t293;
t291 = pkin(3) * t248 + pkin(8) * t246;
t239 = -qJD(3) * t267 + V_base(5);
t290 = rSges(3,1) * t261 - rSges(3,2) * t260;
t289 = rSges(4,1) * t248 - rSges(4,2) * t246;
t288 = Icges(3,1) * t261 - t321;
t287 = Icges(4,1) * t248 - t319;
t286 = -Icges(3,2) * t260 + t320;
t285 = -Icges(4,2) * t246 + t318;
t284 = Icges(3,5) * t261 - Icges(3,6) * t260;
t283 = Icges(4,5) * t248 - Icges(4,6) * t246;
t237 = pkin(1) * t267 + t265 * qJ(2);
t282 = -qJD(2) * t267 + t251 * t237 + V_base(2);
t202 = t265 * t301 + t239;
t281 = qJ(5) * t246 + t248 * t324;
t280 = (-Icges(4,3) * t267 + t265 * t283) * t239 + (Icges(4,3) * t265 + t267 * t283) * t240 + (Icges(4,5) * t246 + Icges(4,6) * t248) * t251;
t279 = pkin(9) * t246 + t248 * t303;
t172 = pkin(7) * t265 + t267 * t325;
t278 = V_base(4) * t171 + (-t172 - t237) * V_base(5) + t298;
t277 = (-Icges(3,3) * t267 + t265 * t284) * V_base(5) + (Icges(3,3) * t265 + t267 * t284) * V_base(4) + (Icges(3,5) * t260 + Icges(3,6) * t261) * t251;
t199 = t291 * t265;
t213 = pkin(3) * t246 - pkin(8) * t248;
t276 = t239 * t213 + (-t199 + t304) * t251 + t292;
t200 = t291 * t267;
t275 = t240 * t199 - t200 * t239 + t278;
t274 = t251 * t172 + (-pkin(6) - t327) * V_base(4) + t282;
t150 = -qJ(5) * t248 + t246 * t324;
t273 = t202 * t150 + t267 * t300 + t276;
t272 = t251 * t200 - t240 * t213 + t274;
t139 = -pkin(4) * t313 + t265 * t281;
t271 = -qJD(5) * t248 + t203 * t139 + t275;
t140 = pkin(4) * t308 + t267 * t281;
t216 = -qJD(4) * t248 + t251;
t270 = t216 * t140 + t265 * t300 + t272;
t180 = -Icges(4,6) * t267 + t265 * t285;
t181 = Icges(4,6) * t265 + t267 * t285;
t182 = -Icges(4,5) * t267 + t265 * t287;
t183 = Icges(4,5) * t265 + t267 * t287;
t210 = Icges(4,2) * t248 + t319;
t211 = Icges(4,1) * t246 + t318;
t269 = (-t181 * t246 + t183 * t248) * t240 + (-t180 * t246 + t182 * t248) * t239 + (-t210 * t246 + t211 * t248) * t251;
t193 = -Icges(3,6) * t267 + t265 * t286;
t194 = Icges(3,6) * t265 + t267 * t286;
t195 = -Icges(3,5) * t267 + t265 * t288;
t196 = Icges(3,5) * t265 + t267 * t288;
t222 = Icges(3,2) * t261 + t321;
t223 = Icges(3,1) * t260 + t320;
t268 = (-t194 * t260 + t196 * t261) * V_base(4) + (-t193 * t260 + t195 * t261) * V_base(5) + (-t222 * t260 + t223 * t261) * t251;
t254 = Icges(2,4) * t267;
t238 = rSges(2,1) * t267 - t265 * rSges(2,2);
t236 = t265 * rSges(2,1) + rSges(2,2) * t267;
t230 = Icges(2,1) * t267 - t322;
t229 = Icges(2,1) * t265 + t254;
t228 = -Icges(2,2) * t265 + t254;
t227 = Icges(2,2) * t267 + t322;
t226 = Icges(2,5) * t267 - Icges(2,6) * t265;
t225 = Icges(2,5) * t265 + Icges(2,6) * t267;
t224 = rSges(3,1) * t260 + rSges(3,2) * t261;
t219 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t218 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t217 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t212 = rSges(4,1) * t246 + rSges(4,2) * t248;
t201 = (-qJD(4) - qJD(6)) * t248 + t251;
t198 = t265 * rSges(3,3) + t267 * t290;
t197 = -rSges(3,3) * t267 + t265 * t290;
t185 = t265 * rSges(4,3) + t267 * t289;
t184 = -rSges(4,3) * t267 + t265 * t289;
t176 = t243 * t315 + t312;
t175 = -t242 * t315 + t311;
t174 = -t242 * t267 + t248 * t311;
t173 = -t243 * t267 - t248 * t312;
t170 = -rSges(5,3) * t248 + (rSges(5,1) * t266 - rSges(5,2) * t264) * t246;
t169 = V_base(5) * rSges(2,3) - t236 * t251 + t297;
t168 = t238 * t251 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t167 = -Icges(5,5) * t248 + (Icges(5,1) * t266 - Icges(5,4) * t264) * t246;
t166 = -Icges(5,6) * t248 + (Icges(5,4) * t266 - Icges(5,2) * t264) * t246;
t164 = t267 * t299 + t203;
t163 = t265 * t299 + t202;
t161 = t236 * V_base(4) - t238 * V_base(5) + V_base(3);
t158 = -rSges(6,3) * t248 + (rSges(6,1) * t249 - rSges(6,2) * t247) * t246;
t157 = -Icges(6,5) * t248 + (Icges(6,1) * t249 - Icges(6,4) * t247) * t246;
t156 = -Icges(6,6) * t248 + (Icges(6,4) * t249 - Icges(6,2) * t247) * t246;
t154 = -rSges(7,3) * t248 + (rSges(7,1) * t243 - rSges(7,2) * t242) * t246;
t153 = -Icges(7,5) * t248 + (Icges(7,1) * t243 - Icges(7,4) * t242) * t246;
t152 = -Icges(7,6) * t248 + (Icges(7,4) * t243 - Icges(7,2) * t242) * t246;
t151 = -Icges(7,3) * t248 + (Icges(7,5) * t243 - Icges(7,6) * t242) * t246;
t149 = -pkin(9) * t248 + t246 * t303;
t148 = t207 * rSges(5,1) + t206 * rSges(5,2) + rSges(5,3) * t316;
t147 = rSges(5,1) * t205 + rSges(5,2) * t204 + rSges(5,3) * t317;
t146 = Icges(5,1) * t207 + Icges(5,4) * t206 + Icges(5,5) * t316;
t145 = Icges(5,1) * t205 + Icges(5,4) * t204 + Icges(5,5) * t317;
t144 = Icges(5,4) * t207 + Icges(5,2) * t206 + Icges(5,6) * t316;
t143 = Icges(5,4) * t205 + Icges(5,2) * t204 + Icges(5,6) * t317;
t137 = t190 * rSges(6,1) + t189 * rSges(6,2) + rSges(6,3) * t316;
t136 = rSges(6,1) * t188 + rSges(6,2) * t187 + rSges(6,3) * t317;
t135 = Icges(6,1) * t190 + Icges(6,4) * t189 + Icges(6,5) * t316;
t134 = Icges(6,1) * t188 + Icges(6,4) * t187 + Icges(6,5) * t317;
t133 = Icges(6,4) * t190 + Icges(6,2) * t189 + Icges(6,6) * t316;
t132 = Icges(6,4) * t188 + Icges(6,2) * t187 + Icges(6,6) * t317;
t129 = t176 * rSges(7,1) + t175 * rSges(7,2) + rSges(7,3) * t316;
t128 = rSges(7,1) * t174 + rSges(7,2) * t173 + rSges(7,3) * t317;
t127 = Icges(7,1) * t176 + Icges(7,4) * t175 + Icges(7,5) * t316;
t126 = Icges(7,1) * t174 + Icges(7,4) * t173 + Icges(7,5) * t317;
t125 = Icges(7,4) * t176 + Icges(7,2) * t175 + Icges(7,6) * t316;
t124 = Icges(7,4) * t174 + Icges(7,2) * t173 + Icges(7,6) * t317;
t123 = Icges(7,5) * t176 + Icges(7,6) * t175 + Icges(7,3) * t316;
t122 = Icges(7,5) * t174 + Icges(7,6) * t173 + Icges(7,3) * t317;
t120 = t224 * V_base(5) + (-t197 - t235) * t251 + t293;
t119 = t251 * t198 + (-pkin(6) - t224) * V_base(4) + t282;
t117 = t197 * V_base(4) + (-t198 - t237) * V_base(5) + t298;
t116 = t265 * t294 + t267 * t279;
t115 = t265 * t279 - t267 * t294;
t114 = t212 * t239 + (-t184 + t304) * t251 + t292;
t113 = t251 * t185 - t240 * t212 + t274;
t112 = t184 * t240 - t185 * t239 + t278;
t111 = -t147 * t216 + t170 * t202 + t276;
t110 = t216 * t148 - t203 * t170 + t272;
t109 = t147 * t203 - t148 * t202 + t275;
t108 = t158 * t202 + (-t136 - t139) * t216 + t273;
t107 = t216 * t137 + (-t150 - t158) * t203 + t270;
t106 = t136 * t203 + (-t137 - t140) * t202 + t271;
t105 = -t128 * t201 + t149 * t202 + t154 * t163 + (-t115 - t139) * t216 + t273;
t104 = t216 * t116 + t201 * t129 - t164 * t154 + (-t149 - t150) * t203 + t270;
t103 = t115 * t203 + t128 * t164 - t129 * t163 + (-t116 - t140) * t202 + t271;
t1 = t240 * (t280 * t265 + t269 * t267) / 0.2e1 + t239 * (t269 * t265 - t280 * t267) / 0.2e1 + m(1) * (t217 ^ 2 + t218 ^ 2 + t219 ^ 2) / 0.2e1 + t164 * ((t123 * t316 + t175 * t125 + t176 * t127) * t164 + (t122 * t316 + t175 * t124 + t176 * t126) * t163 + (t151 * t316 + t175 * t152 + t176 * t153) * t201) / 0.2e1 + t163 * ((t123 * t317 + t125 * t173 + t127 * t174) * t164 + (t122 * t317 + t173 * t124 + t174 * t126) * t163 + (t151 * t317 + t152 * t173 + t153 * t174) * t201) / 0.2e1 + m(2) * (t161 ^ 2 + t168 ^ 2 + t169 ^ 2) / 0.2e1 + m(3) * (t117 ^ 2 + t119 ^ 2 + t120 ^ 2) / 0.2e1 + m(4) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(7) * (t103 ^ 2 + t104 ^ 2 + t105 ^ 2) / 0.2e1 + m(6) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(5) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + t201 * ((-t122 * t163 - t123 * t164 - t151 * t201) * t248 + ((-t125 * t242 + t127 * t243) * t164 + (-t124 * t242 + t126 * t243) * t163 + (-t152 * t242 + t153 * t243) * t201) * t246) / 0.2e1 + ((t156 * t187 + t157 * t188 + t166 * t204 + t167 * t205 + t317 * t332) * t216 + (t133 * t187 + t135 * t188 + t144 * t204 + t146 * t205 + t317 * t333) * t203 + (t187 * t132 + t188 * t134 + t204 * t143 + t205 * t145 + t334 * t317) * t202) * t202 / 0.2e1 + ((t189 * t156 + t190 * t157 + t206 * t166 + t207 * t167 + t316 * t332) * t216 + (t189 * t133 + t190 * t135 + t206 * t144 + t207 * t146 + t333 * t316) * t203 + (t189 * t132 + t190 * t134 + t206 * t143 + t207 * t145 + t316 * t334) * t202) * t203 / 0.2e1 + ((-t202 * t334 - t333 * t203 - t332 * t216) * t248 + ((-t156 * t247 + t157 * t249 - t166 * t264 + t167 * t266) * t216 + (-t133 * t247 + t135 * t249 - t144 * t264 + t146 * t266) * t203 + (-t132 * t247 + t134 * t249 - t143 * t264 + t145 * t266) * t202) * t246) * t216 / 0.2e1 + ((t181 * t248 + t183 * t246) * t240 + (t180 * t248 + t182 * t246) * t239 + (t193 * t261 + t195 * t260 + t225) * V_base(5) + (t194 * t261 + t196 * t260 + t226) * V_base(4) + (t210 * t248 + t211 * t246 + t222 * t261 + t223 * t260 + Icges(2,3)) * t251) * t251 / 0.2e1 + (t226 * t251 + t265 * t277 + t268 * t267 + (-t265 * t227 + t229 * t267 + Icges(1,4)) * V_base(5) + (-t265 * t228 + t230 * t267 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t225 * t251 + t265 * t268 - t277 * t267 + (t227 * t267 + t265 * t229 + Icges(1,2)) * V_base(5) + (t228 * t267 + t265 * t230 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
