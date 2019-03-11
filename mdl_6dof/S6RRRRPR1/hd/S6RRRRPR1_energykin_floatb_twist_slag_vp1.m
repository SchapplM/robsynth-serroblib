% Calculate kinetic energy for
% S6RRRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 21:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR1_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:52:51
% EndTime: 2019-03-09 21:52:53
% DurationCPUTime: 2.77s
% Computational Cost: add. (2257->335), mult. (1856->489), div. (0->0), fcn. (1638->12), ass. (0->178)
t326 = Icges(5,3) + Icges(6,3);
t244 = qJ(2) + qJ(3);
t240 = qJ(4) + t244;
t227 = pkin(11) + t240;
t224 = sin(t227);
t225 = cos(t227);
t228 = sin(t240);
t229 = cos(t240);
t325 = Icges(5,5) * t229 + Icges(6,5) * t225 - Icges(5,6) * t228 - Icges(6,6) * t224;
t247 = sin(qJ(1));
t250 = cos(qJ(1));
t306 = Icges(6,4) * t225;
t273 = -Icges(6,2) * t224 + t306;
t139 = -Icges(6,6) * t250 + t247 * t273;
t140 = Icges(6,6) * t247 + t250 * t273;
t307 = Icges(6,4) * t224;
t277 = Icges(6,1) * t225 - t307;
t141 = -Icges(6,5) * t250 + t247 * t277;
t142 = Icges(6,5) * t247 + t250 * t277;
t308 = Icges(5,4) * t229;
t274 = -Icges(5,2) * t228 + t308;
t150 = -Icges(5,6) * t250 + t247 * t274;
t151 = Icges(5,6) * t247 + t250 * t274;
t309 = Icges(5,4) * t228;
t278 = Icges(5,1) * t229 - t309;
t152 = -Icges(5,5) * t250 + t247 * t278;
t153 = Icges(5,5) * t247 + t250 * t278;
t181 = Icges(6,2) * t225 + t307;
t182 = Icges(6,1) * t224 + t306;
t293 = -qJD(2) - qJD(3);
t188 = V_base(5) + (-qJD(4) + t293) * t250;
t223 = qJD(2) * t247 + V_base(4);
t201 = qJD(3) * t247 + t223;
t189 = qJD(4) * t247 + t201;
t191 = Icges(5,2) * t229 + t309;
t192 = Icges(5,1) * t228 + t308;
t231 = V_base(6) + qJD(1);
t324 = (-t181 * t224 + t182 * t225 - t191 * t228 + t192 * t229) * t231 + (-t140 * t224 + t142 * t225 - t151 * t228 + t153 * t229) * t189 + (-t139 * t224 + t141 * t225 - t150 * t228 + t152 * t229) * t188;
t323 = (Icges(5,5) * t228 + Icges(6,5) * t224 + Icges(5,6) * t229 + Icges(6,6) * t225) * t231 + (t326 * t247 + t325 * t250) * t189 + (t325 * t247 - t326 * t250) * t188;
t246 = sin(qJ(2));
t319 = pkin(2) * t246;
t236 = sin(t244);
t318 = pkin(3) * t236;
t317 = pkin(4) * t228;
t249 = cos(qJ(2));
t316 = t249 * pkin(2);
t314 = Icges(2,4) * t247;
t313 = Icges(3,4) * t246;
t312 = Icges(3,4) * t249;
t311 = Icges(4,4) * t236;
t237 = cos(t244);
t310 = Icges(4,4) * t237;
t305 = t224 * t247;
t304 = t224 * t250;
t245 = sin(qJ(6));
t303 = t245 * t247;
t302 = t245 * t250;
t248 = cos(qJ(6));
t301 = t247 * t248;
t300 = t248 * t250;
t157 = -pkin(8) * t250 + t247 * t316;
t220 = t247 * pkin(1) - t250 * pkin(7);
t299 = -t157 - t220;
t298 = pkin(4) * t229;
t297 = pkin(3) * t237;
t294 = qJD(6) * t224;
t292 = V_base(5) * pkin(6) + V_base(1);
t127 = -pkin(9) * t250 + t247 * t297;
t289 = -t127 + t299;
t222 = -qJD(2) * t250 + V_base(5);
t288 = t222 * t319 + t292;
t123 = -qJ(5) * t250 + t247 * t298;
t287 = -t123 + t289;
t200 = t250 * t293 + V_base(5);
t286 = t200 * t318 + t288;
t285 = pkin(5) * t225 + pkin(10) * t224;
t284 = rSges(3,1) * t249 - rSges(3,2) * t246;
t283 = rSges(4,1) * t237 - rSges(4,2) * t236;
t282 = rSges(5,1) * t229 - rSges(5,2) * t228;
t281 = rSges(6,1) * t225 - rSges(6,2) * t224;
t280 = Icges(3,1) * t249 - t313;
t279 = Icges(4,1) * t237 - t311;
t276 = -Icges(3,2) * t246 + t312;
t275 = -Icges(4,2) * t236 + t310;
t272 = Icges(3,5) * t249 - Icges(3,6) * t246;
t271 = Icges(4,5) * t237 - Icges(4,6) * t236;
t221 = t250 * pkin(1) + t247 * pkin(7);
t268 = -V_base(4) * pkin(6) + t231 * t221 + V_base(2);
t267 = V_base(4) * t220 - t221 * V_base(5) + V_base(3);
t266 = qJD(5) * t247 + t188 * t317 + t286;
t263 = (-Icges(4,3) * t250 + t247 * t271) * t200 + (Icges(4,3) * t247 + t250 * t271) * t201 + (Icges(4,5) * t236 + Icges(4,6) * t237) * t231;
t262 = (-Icges(3,3) * t250 + t247 * t272) * t222 + (Icges(3,3) * t247 + t250 * t272) * t223 + (Icges(3,5) * t246 + Icges(3,6) * t249) * t231;
t158 = pkin(8) * t247 + t250 * t316;
t261 = t223 * t157 - t158 * t222 + t267;
t260 = t231 * t158 - t223 * t319 + t268;
t128 = pkin(9) * t247 + t250 * t297;
t259 = t201 * t127 - t128 * t200 + t261;
t258 = t189 * t123 + t259;
t257 = t231 * t128 - t201 * t318 + t260;
t124 = qJ(5) * t247 + t250 * t298;
t256 = -qJD(5) * t250 + t231 * t124 + t257;
t163 = -Icges(4,6) * t250 + t247 * t275;
t164 = Icges(4,6) * t247 + t250 * t275;
t165 = -Icges(4,5) * t250 + t247 * t279;
t166 = Icges(4,5) * t247 + t250 * t279;
t197 = Icges(4,2) * t237 + t311;
t198 = Icges(4,1) * t236 + t310;
t253 = (-t164 * t236 + t166 * t237) * t201 + (-t163 * t236 + t165 * t237) * t200 + (-t197 * t236 + t198 * t237) * t231;
t175 = -Icges(3,6) * t250 + t247 * t276;
t176 = Icges(3,6) * t247 + t250 * t276;
t177 = -Icges(3,5) * t250 + t247 * t280;
t178 = Icges(3,5) * t247 + t250 * t280;
t211 = Icges(3,2) * t249 + t313;
t214 = Icges(3,1) * t246 + t312;
t252 = (-t176 * t246 + t178 * t249) * t223 + (-t175 * t246 + t177 * t249) * t222 + (-t211 * t246 + t214 * t249) * t231;
t239 = Icges(2,4) * t250;
t219 = rSges(2,1) * t250 - rSges(2,2) * t247;
t218 = rSges(2,1) * t247 + rSges(2,2) * t250;
t217 = rSges(3,1) * t246 + rSges(3,2) * t249;
t216 = Icges(2,1) * t250 - t314;
t215 = Icges(2,1) * t247 + t239;
t213 = -Icges(2,2) * t247 + t239;
t212 = Icges(2,2) * t250 + t314;
t207 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t206 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t205 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t199 = rSges(4,1) * t236 + rSges(4,2) * t237;
t195 = -qJD(6) * t225 + t231;
t193 = rSges(5,1) * t228 + rSges(5,2) * t229;
t187 = pkin(5) * t224 - pkin(10) * t225;
t185 = rSges(6,1) * t224 + rSges(6,2) * t225;
t184 = rSges(3,3) * t247 + t250 * t284;
t183 = -rSges(3,3) * t250 + t247 * t284;
t172 = t225 * t300 + t303;
t171 = -t225 * t302 + t301;
t170 = t225 * t301 - t302;
t169 = -t225 * t303 - t300;
t168 = rSges(4,3) * t247 + t250 * t283;
t167 = -rSges(4,3) * t250 + t247 * t283;
t160 = t285 * t250;
t159 = t285 * t247;
t155 = rSges(5,3) * t247 + t250 * t282;
t154 = -rSges(5,3) * t250 + t247 * t282;
t147 = V_base(5) * rSges(2,3) - t218 * t231 + t292;
t146 = t219 * t231 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t145 = t218 * V_base(4) - t219 * V_base(5) + V_base(3);
t144 = rSges(6,3) * t247 + t250 * t281;
t143 = -rSges(6,3) * t250 + t247 * t281;
t135 = t250 * t294 + t189;
t134 = t247 * t294 + t188;
t132 = -rSges(7,3) * t225 + (rSges(7,1) * t248 - rSges(7,2) * t245) * t224;
t131 = -Icges(7,5) * t225 + (Icges(7,1) * t248 - Icges(7,4) * t245) * t224;
t130 = -Icges(7,6) * t225 + (Icges(7,4) * t248 - Icges(7,2) * t245) * t224;
t129 = -Icges(7,3) * t225 + (Icges(7,5) * t248 - Icges(7,6) * t245) * t224;
t122 = rSges(7,1) * t172 + rSges(7,2) * t171 + rSges(7,3) * t304;
t121 = rSges(7,1) * t170 + rSges(7,2) * t169 + rSges(7,3) * t305;
t120 = Icges(7,1) * t172 + Icges(7,4) * t171 + Icges(7,5) * t304;
t119 = Icges(7,1) * t170 + Icges(7,4) * t169 + Icges(7,5) * t305;
t118 = Icges(7,4) * t172 + Icges(7,2) * t171 + Icges(7,6) * t304;
t117 = Icges(7,4) * t170 + Icges(7,2) * t169 + Icges(7,6) * t305;
t116 = Icges(7,5) * t172 + Icges(7,6) * t171 + Icges(7,3) * t304;
t115 = Icges(7,5) * t170 + Icges(7,6) * t169 + Icges(7,3) * t305;
t113 = t217 * t222 + (-t183 - t220) * t231 + t292;
t112 = t184 * t231 - t217 * t223 + t268;
t110 = t183 * t223 - t184 * t222 + t267;
t109 = t199 * t200 + (-t167 + t299) * t231 + t288;
t108 = t168 * t231 - t199 * t201 + t260;
t107 = t167 * t201 - t168 * t200 + t261;
t106 = t188 * t193 + (-t154 + t289) * t231 + t286;
t105 = t155 * t231 - t189 * t193 + t257;
t104 = t154 * t189 - t155 * t188 + t259;
t103 = t185 * t188 + (-t143 + t287) * t231 + t266;
t102 = t144 * t231 + (-t185 - t317) * t189 + t256;
t101 = -t121 * t195 + t132 * t134 + t187 * t188 + (-t159 + t287) * t231 + t266;
t100 = t122 * t195 - t132 * t135 + t160 * t231 + (-t187 - t317) * t189 + t256;
t99 = t143 * t189 + (-t124 - t144) * t188 + t258;
t98 = t121 * t135 - t122 * t134 + t159 * t189 + (-t124 - t160) * t188 + t258;
t1 = t223 * (t262 * t247 + t252 * t250) / 0.2e1 + t222 * (t247 * t252 - t250 * t262) / 0.2e1 + t201 * (t263 * t247 + t253 * t250) / 0.2e1 + t200 * (t253 * t247 - t263 * t250) / 0.2e1 + t195 * ((-t115 * t134 - t116 * t135 - t129 * t195) * t225 + ((-t118 * t245 + t120 * t248) * t135 + (-t117 * t245 + t119 * t248) * t134 + (-t130 * t245 + t131 * t248) * t195) * t224) / 0.2e1 + m(2) * (t145 ^ 2 + t146 ^ 2 + t147 ^ 2) / 0.2e1 + m(3) * (t110 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(6) * (t102 ^ 2 + t103 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(4) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(1) * (t205 ^ 2 + t206 ^ 2 + t207 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t101 ^ 2 + t98 ^ 2) / 0.2e1 + t135 * ((t116 * t304 + t171 * t118 + t172 * t120) * t135 + (t115 * t304 + t117 * t171 + t119 * t172) * t134 + (t129 * t304 + t130 * t171 + t131 * t172) * t195) / 0.2e1 + t134 * ((t116 * t305 + t118 * t169 + t120 * t170) * t135 + (t115 * t305 + t169 * t117 + t170 * t119) * t134 + (t129 * t305 + t169 * t130 + t170 * t131) * t195) / 0.2e1 + (t324 * t247 - t323 * t250) * t188 / 0.2e1 + (t323 * t247 + t324 * t250) * t189 / 0.2e1 + ((-t212 * t247 + t215 * t250 + Icges(1,4)) * V_base(5) + (-t247 * t213 + t250 * t216 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t250 * t212 + t247 * t215 + Icges(1,2)) * V_base(5) + (t213 * t250 + t216 * t247 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t176 * t249 + t178 * t246) * t223 + (t175 * t249 + t177 * t246) * t222 + (t164 * t237 + t166 * t236) * t201 + (t163 * t237 + t165 * t236) * t200 + (t140 * t225 + t142 * t224 + t151 * t229 + t153 * t228) * t189 + (t139 * t225 + t141 * t224 + t150 * t229 + t152 * t228) * t188 + (t225 * t181 + t224 * t182 + t229 * t191 + t228 * t192 + t237 * t197 + t236 * t198 + t249 * t211 + t246 * t214 + Icges(2,3)) * t231) * t231 / 0.2e1 + V_base(4) * t231 * (Icges(2,5) * t250 - Icges(2,6) * t247) + t231 * V_base(5) * (Icges(2,5) * t247 + Icges(2,6) * t250) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
