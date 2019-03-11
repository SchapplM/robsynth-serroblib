% Calculate kinetic energy for
% S6RRRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:05
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR1_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR1_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR1_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR1_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR1_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR1_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:02:06
% EndTime: 2019-03-09 18:02:09
% DurationCPUTime: 2.77s
% Computational Cost: add. (2239->335), mult. (1838->489), div. (0->0), fcn. (1620->12), ass. (0->178)
t326 = Icges(4,3) + Icges(5,3);
t244 = qJ(2) + qJ(3);
t231 = pkin(11) + t244;
t226 = sin(t231);
t227 = cos(t231);
t237 = sin(t244);
t238 = cos(t244);
t325 = Icges(4,5) * t238 + Icges(5,5) * t227 - Icges(4,6) * t237 - Icges(5,6) * t226;
t247 = sin(qJ(1));
t250 = cos(qJ(1));
t308 = Icges(5,4) * t227;
t275 = -Icges(5,2) * t226 + t308;
t148 = -Icges(5,6) * t250 + t247 * t275;
t149 = Icges(5,6) * t247 + t250 * t275;
t309 = Icges(5,4) * t226;
t279 = Icges(5,1) * t227 - t309;
t150 = -Icges(5,5) * t250 + t247 * t279;
t151 = Icges(5,5) * t247 + t250 * t279;
t310 = Icges(4,4) * t238;
t276 = -Icges(4,2) * t237 + t310;
t162 = -Icges(4,6) * t250 + t247 * t276;
t163 = Icges(4,6) * t247 + t250 * t276;
t311 = Icges(4,4) * t237;
t280 = Icges(4,1) * t238 - t311;
t164 = -Icges(4,5) * t250 + t247 * t280;
t165 = Icges(4,5) * t247 + t250 * t280;
t191 = Icges(5,2) * t227 + t309;
t192 = Icges(5,1) * t226 + t308;
t197 = Icges(4,2) * t238 + t311;
t198 = Icges(4,1) * t237 + t310;
t293 = -qJD(2) - qJD(3);
t200 = t250 * t293 + V_base(5);
t223 = qJD(2) * t247 + V_base(4);
t201 = qJD(3) * t247 + t223;
t232 = V_base(6) + qJD(1);
t324 = (-t191 * t226 + t192 * t227 - t197 * t237 + t198 * t238) * t232 + (-t149 * t226 + t151 * t227 - t163 * t237 + t165 * t238) * t201 + (-t148 * t226 + t150 * t227 - t162 * t237 + t164 * t238) * t200;
t323 = (Icges(4,5) * t237 + Icges(5,5) * t226 + Icges(4,6) * t238 + Icges(5,6) * t227) * t232 + (t326 * t247 + t325 * t250) * t201 + (t325 * t247 - t326 * t250) * t200;
t246 = sin(qJ(2));
t319 = pkin(2) * t246;
t318 = pkin(3) * t237;
t317 = pkin(4) * t226;
t249 = cos(qJ(2));
t316 = t249 * pkin(2);
t314 = Icges(2,4) * t247;
t313 = Icges(3,4) * t246;
t312 = Icges(3,4) * t249;
t229 = qJ(5) + t231;
t224 = sin(t229);
t307 = Icges(6,4) * t224;
t225 = cos(t229);
t306 = Icges(6,4) * t225;
t305 = t224 * t247;
t304 = t224 * t250;
t245 = sin(qJ(6));
t303 = t245 * t247;
t302 = t245 * t250;
t248 = cos(qJ(6));
t301 = t247 * t248;
t300 = t248 * t250;
t156 = -pkin(8) * t250 + t247 * t316;
t220 = t247 * pkin(1) - t250 * pkin(7);
t299 = -t156 - t220;
t298 = pkin(4) * t227;
t297 = pkin(3) * t238;
t294 = qJD(6) * t224;
t292 = V_base(5) * pkin(6) + V_base(1);
t127 = -qJ(4) * t250 + t247 * t297;
t289 = -t127 + t299;
t222 = -qJD(2) * t250 + V_base(5);
t288 = t222 * t319 + t292;
t123 = -pkin(9) * t250 + t247 * t298;
t287 = -t123 + t289;
t286 = pkin(5) * t225 + pkin(10) * t224;
t285 = rSges(3,1) * t249 - rSges(3,2) * t246;
t284 = rSges(4,1) * t238 - rSges(4,2) * t237;
t283 = rSges(5,1) * t227 - rSges(5,2) * t226;
t282 = rSges(6,1) * t225 - rSges(6,2) * t224;
t189 = qJD(5) * t247 + t201;
t281 = Icges(3,1) * t249 - t313;
t278 = Icges(6,1) * t225 - t307;
t277 = -Icges(3,2) * t246 + t312;
t274 = -Icges(6,2) * t224 + t306;
t273 = Icges(3,5) * t249 - Icges(3,6) * t246;
t270 = Icges(6,5) * t225 - Icges(6,6) * t224;
t269 = qJD(4) * t247 + t200 * t318 + t288;
t221 = t250 * pkin(1) + t247 * pkin(7);
t268 = -V_base(4) * pkin(6) + t232 * t221 + V_base(2);
t267 = V_base(4) * t220 - t221 * V_base(5) + V_base(3);
t266 = t200 * t317 + t269;
t188 = V_base(5) + (-qJD(5) + t293) * t250;
t265 = (-Icges(6,3) * t250 + t247 * t270) * t188 + (Icges(6,3) * t247 + t250 * t270) * t189 + (Icges(6,5) * t224 + Icges(6,6) * t225) * t232;
t262 = (-Icges(3,3) * t250 + t247 * t273) * t222 + (Icges(3,3) * t247 + t250 * t273) * t223 + (Icges(3,5) * t246 + Icges(3,6) * t249) * t232;
t157 = pkin(8) * t247 + t250 * t316;
t261 = t223 * t156 - t157 * t222 + t267;
t260 = t232 * t157 - t223 * t319 + t268;
t259 = t201 * t127 + t261;
t128 = qJ(4) * t247 + t250 * t297;
t258 = -qJD(4) * t250 + t232 * t128 + t260;
t124 = pkin(9) * t247 + t250 * t298;
t257 = t201 * t123 + (-t124 - t128) * t200 + t259;
t256 = t232 * t124 + (-t317 - t318) * t201 + t258;
t139 = -Icges(6,6) * t250 + t247 * t274;
t140 = Icges(6,6) * t247 + t250 * t274;
t141 = -Icges(6,5) * t250 + t247 * t278;
t142 = Icges(6,5) * t247 + t250 * t278;
t181 = Icges(6,2) * t225 + t307;
t182 = Icges(6,1) * t224 + t306;
t255 = (-t140 * t224 + t142 * t225) * t189 + (-t139 * t224 + t141 * t225) * t188 + (-t181 * t224 + t182 * t225) * t232;
t176 = -Icges(3,6) * t250 + t247 * t277;
t177 = Icges(3,6) * t247 + t250 * t277;
t178 = -Icges(3,5) * t250 + t247 * t281;
t179 = Icges(3,5) * t247 + t250 * t281;
t211 = Icges(3,2) * t249 + t313;
t214 = Icges(3,1) * t246 + t312;
t252 = (-t177 * t246 + t179 * t249) * t223 + (-t176 * t246 + t178 * t249) * t222 + (-t211 * t246 + t214 * t249) * t232;
t240 = Icges(2,4) * t250;
t219 = rSges(2,1) * t250 - rSges(2,2) * t247;
t218 = rSges(2,1) * t247 + rSges(2,2) * t250;
t217 = rSges(3,1) * t246 + rSges(3,2) * t249;
t216 = Icges(2,1) * t250 - t314;
t215 = Icges(2,1) * t247 + t240;
t213 = -Icges(2,2) * t247 + t240;
t212 = Icges(2,2) * t250 + t314;
t207 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t206 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t205 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t199 = rSges(4,1) * t237 + rSges(4,2) * t238;
t195 = -qJD(6) * t225 + t232;
t193 = rSges(5,1) * t226 + rSges(5,2) * t227;
t187 = pkin(5) * t224 - pkin(10) * t225;
t185 = rSges(6,1) * t224 + rSges(6,2) * t225;
t184 = rSges(3,3) * t247 + t250 * t285;
t183 = -rSges(3,3) * t250 + t247 * t285;
t172 = t225 * t300 + t303;
t171 = -t225 * t302 + t301;
t170 = t225 * t301 - t302;
t169 = -t225 * t303 - t300;
t167 = rSges(4,3) * t247 + t250 * t284;
t166 = -rSges(4,3) * t250 + t247 * t284;
t159 = t286 * t250;
t158 = t286 * t247;
t155 = rSges(5,3) * t247 + t250 * t283;
t154 = -rSges(5,3) * t250 + t247 * t283;
t153 = V_base(5) * rSges(2,3) - t218 * t232 + t292;
t152 = t219 * t232 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t145 = t218 * V_base(4) - t219 * V_base(5) + V_base(3);
t144 = rSges(6,3) * t247 + t250 * t282;
t143 = -rSges(6,3) * t250 + t247 * t282;
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
t113 = t217 * t222 + (-t183 - t220) * t232 + t292;
t112 = t184 * t232 - t217 * t223 + t268;
t110 = t183 * t223 - t184 * t222 + t267;
t109 = t199 * t200 + (-t166 + t299) * t232 + t288;
t108 = t167 * t232 - t199 * t201 + t260;
t107 = t166 * t201 - t167 * t200 + t261;
t106 = t193 * t200 + (-t154 + t289) * t232 + t269;
t105 = t155 * t232 + (-t193 - t318) * t201 + t258;
t104 = t154 * t201 + (-t128 - t155) * t200 + t259;
t103 = t185 * t188 + (-t143 + t287) * t232 + t266;
t102 = t144 * t232 - t185 * t189 + t256;
t101 = -t121 * t195 + t132 * t134 + t187 * t188 + (-t158 + t287) * t232 + t266;
t100 = t122 * t195 - t132 * t135 + t159 * t232 - t187 * t189 + t256;
t99 = t143 * t189 - t144 * t188 + t257;
t98 = t121 * t135 - t122 * t134 + t158 * t189 - t159 * t188 + t257;
t1 = t195 * ((-t115 * t134 - t116 * t135 - t129 * t195) * t225 + ((-t118 * t245 + t120 * t248) * t135 + (-t117 * t245 + t119 * t248) * t134 + (-t130 * t245 + t131 * t248) * t195) * t224) / 0.2e1 + t222 * (t252 * t247 - t262 * t250) / 0.2e1 + t189 * (t265 * t247 + t255 * t250) / 0.2e1 + t188 * (t255 * t247 - t265 * t250) / 0.2e1 + t223 * (t262 * t247 + t252 * t250) / 0.2e1 + m(1) * (t205 ^ 2 + t206 ^ 2 + t207 ^ 2) / 0.2e1 + m(2) * (t145 ^ 2 + t152 ^ 2 + t153 ^ 2) / 0.2e1 + m(5) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(4) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + m(3) * (t110 ^ 2 + t112 ^ 2 + t113 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t101 ^ 2 + t98 ^ 2) / 0.2e1 + m(6) * (t102 ^ 2 + t103 ^ 2 + t99 ^ 2) / 0.2e1 + t135 * ((t116 * t304 + t171 * t118 + t172 * t120) * t135 + (t115 * t304 + t117 * t171 + t119 * t172) * t134 + (t129 * t304 + t130 * t171 + t131 * t172) * t195) / 0.2e1 + t134 * ((t116 * t305 + t118 * t169 + t120 * t170) * t135 + (t115 * t305 + t169 * t117 + t170 * t119) * t134 + (t129 * t305 + t130 * t169 + t131 * t170) * t195) / 0.2e1 + (t324 * t247 - t323 * t250) * t200 / 0.2e1 + (t323 * t247 + t324 * t250) * t201 / 0.2e1 + ((-t212 * t247 + t215 * t250 + Icges(1,4)) * V_base(5) + (-t247 * t213 + t250 * t216 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t250 * t212 + t247 * t215 + Icges(1,2)) * V_base(5) + (t213 * t250 + t216 * t247 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t177 * t249 + t179 * t246) * t223 + (t176 * t249 + t178 * t246) * t222 + (t140 * t225 + t142 * t224) * t189 + (t139 * t225 + t141 * t224) * t188 + (t149 * t227 + t151 * t226 + t163 * t238 + t165 * t237) * t201 + (t148 * t227 + t150 * t226 + t162 * t238 + t164 * t237) * t200 + (t225 * t181 + t224 * t182 + t227 * t191 + t226 * t192 + t238 * t197 + t237 * t198 + t249 * t211 + t246 * t214 + Icges(2,3)) * t232) * t232 / 0.2e1 + t232 * V_base(4) * (Icges(2,5) * t250 - Icges(2,6) * t247) + V_base(5) * t232 * (Icges(2,5) * t247 + Icges(2,6) * t250) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
