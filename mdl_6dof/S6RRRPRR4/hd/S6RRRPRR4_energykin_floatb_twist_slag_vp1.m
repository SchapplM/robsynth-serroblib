% Calculate kinetic energy for
% S6RRRPRR4
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
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:15:18
% EndTime: 2019-03-09 18:15:22
% DurationCPUTime: 3.34s
% Computational Cost: add. (2409->388), mult. (2337->583), div. (0->0), fcn. (2231->12), ass. (0->184)
t261 = sin(qJ(2));
t321 = pkin(2) * t261;
t263 = cos(qJ(2));
t319 = pkin(2) * t263;
t259 = cos(pkin(11));
t318 = t259 * pkin(4);
t262 = sin(qJ(1));
t316 = Icges(2,4) * t262;
t315 = Icges(3,4) * t261;
t314 = Icges(3,4) * t263;
t257 = qJ(2) + qJ(3);
t250 = sin(t257);
t313 = Icges(4,4) * t250;
t251 = cos(t257);
t312 = Icges(4,4) * t251;
t311 = t250 * t262;
t264 = cos(qJ(1));
t310 = t250 * t264;
t309 = t251 * t262;
t308 = t251 * t264;
t258 = sin(pkin(11));
t307 = t258 * t262;
t306 = t258 * t264;
t305 = t259 * t262;
t304 = t259 * t264;
t172 = -pkin(8) * t264 + t262 * t319;
t236 = pkin(1) * t262 - pkin(7) * t264;
t302 = -t172 - t236;
t256 = pkin(11) + qJ(5);
t245 = cos(t256);
t301 = pkin(5) * t245;
t299 = qJD(4) * t250;
t298 = qJD(5) * t250;
t297 = qJD(6) * t250;
t296 = V_base(5) * pkin(6) + V_base(1);
t288 = pkin(3) * t251 + qJ(4) * t250;
t197 = t288 * t262;
t293 = -t197 + t302;
t239 = qJD(2) * t262 + V_base(4);
t247 = V_base(6) + qJD(1);
t244 = sin(t256);
t292 = pkin(5) * t244;
t238 = -qJD(2) * t264 + V_base(5);
t291 = t238 * t321 + t296;
t211 = qJD(3) * t262 + t239;
t290 = rSges(3,1) * t263 - rSges(3,2) * t261;
t289 = rSges(4,1) * t251 - rSges(4,2) * t250;
t184 = t264 * t298 + t211;
t287 = Icges(3,1) * t263 - t315;
t286 = Icges(4,1) * t251 - t313;
t285 = -Icges(3,2) * t261 + t314;
t284 = -Icges(4,2) * t250 + t312;
t283 = Icges(3,5) * t263 - Icges(3,6) * t261;
t282 = Icges(4,5) * t251 - Icges(4,6) * t250;
t208 = pkin(3) * t250 - qJ(4) * t251;
t210 = V_base(5) + (-qJD(2) - qJD(3)) * t264;
t281 = t208 * t210 + t264 * t299 + t291;
t237 = pkin(1) * t264 + pkin(7) * t262;
t280 = -V_base(4) * pkin(6) + t237 * t247 + V_base(2);
t279 = t236 * V_base(4) - t237 * V_base(5) + V_base(3);
t183 = t262 * t298 + t210;
t278 = (-Icges(4,3) * t264 + t262 * t282) * t210 + (Icges(4,3) * t262 + t264 * t282) * t211 + (Icges(4,5) * t250 + Icges(4,6) * t251) * t247;
t277 = (-Icges(3,3) * t264 + t262 * t283) * t238 + (Icges(3,3) * t262 + t264 * t283) * t239 + (Icges(3,5) * t261 + Icges(3,6) * t263) * t247;
t276 = pkin(9) * t250 + t251 * t318;
t275 = pkin(10) * t250 + t251 * t301;
t173 = pkin(8) * t262 + t264 * t319;
t274 = t172 * t239 - t173 * t238 + t279;
t273 = t173 * t247 - t239 * t321 + t280;
t198 = t288 * t264;
t272 = t198 * t247 + t262 * t299 + t273;
t134 = -pkin(4) * t306 + t262 * t276;
t146 = -pkin(9) * t251 + t250 * t318;
t271 = t210 * t146 + (-t134 + t293) * t247 + t281;
t270 = -qJD(4) * t251 + t197 * t211 + t274;
t135 = pkin(4) * t307 + t264 * t276;
t269 = t247 * t135 + (-t146 - t208) * t211 + t272;
t268 = t211 * t134 + (-t135 - t198) * t210 + t270;
t177 = -Icges(4,6) * t264 + t262 * t284;
t178 = Icges(4,6) * t262 + t264 * t284;
t179 = -Icges(4,5) * t264 + t262 * t286;
t180 = Icges(4,5) * t262 + t264 * t286;
t206 = Icges(4,2) * t251 + t313;
t207 = Icges(4,1) * t250 + t312;
t267 = (-t178 * t250 + t180 * t251) * t211 + (-t177 * t250 + t179 * t251) * t210 + (-t206 * t250 + t207 * t251) * t247;
t191 = -Icges(3,6) * t264 + t262 * t285;
t192 = Icges(3,6) * t262 + t264 * t285;
t193 = -Icges(3,5) * t264 + t262 * t287;
t194 = Icges(3,5) * t262 + t264 * t287;
t223 = Icges(3,2) * t263 + t315;
t226 = Icges(3,1) * t261 + t314;
t266 = (-t192 * t261 + t194 * t263) * t239 + (-t191 * t261 + t193 * t263) * t238 + (-t223 * t261 + t226 * t263) * t247;
t253 = Icges(2,4) * t264;
t246 = qJ(6) + t256;
t241 = cos(t246);
t240 = sin(t246);
t231 = rSges(2,1) * t264 - rSges(2,2) * t262;
t230 = rSges(2,1) * t262 + rSges(2,2) * t264;
t229 = rSges(3,1) * t261 + rSges(3,2) * t263;
t228 = Icges(2,1) * t264 - t316;
t227 = Icges(2,1) * t262 + t253;
t225 = -Icges(2,2) * t262 + t253;
t224 = Icges(2,2) * t264 + t316;
t219 = -qJD(5) * t251 + t247;
t218 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t217 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t216 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t209 = rSges(4,1) * t250 + rSges(4,2) * t251;
t203 = t251 * t304 + t307;
t202 = -t251 * t306 + t305;
t201 = t251 * t305 - t306;
t200 = -t251 * t307 - t304;
t199 = (-qJD(5) - qJD(6)) * t251 + t247;
t196 = rSges(3,3) * t262 + t264 * t290;
t195 = -rSges(3,3) * t264 + t262 * t290;
t188 = t244 * t262 + t245 * t308;
t187 = -t244 * t308 + t245 * t262;
t186 = -t244 * t264 + t245 * t309;
t185 = -t244 * t309 - t245 * t264;
t182 = rSges(4,3) * t262 + t264 * t289;
t181 = -rSges(4,3) * t264 + t262 * t289;
t171 = t240 * t262 + t241 * t308;
t170 = -t240 * t308 + t241 * t262;
t169 = -t240 * t264 + t241 * t309;
t168 = -t240 * t309 - t241 * t264;
t167 = -rSges(5,3) * t251 + (rSges(5,1) * t259 - rSges(5,2) * t258) * t250;
t166 = -Icges(5,5) * t251 + (Icges(5,1) * t259 - Icges(5,4) * t258) * t250;
t165 = -Icges(5,6) * t251 + (Icges(5,4) * t259 - Icges(5,2) * t258) * t250;
t164 = -Icges(5,3) * t251 + (Icges(5,5) * t259 - Icges(5,6) * t258) * t250;
t163 = V_base(5) * rSges(2,3) - t230 * t247 + t296;
t162 = t231 * t247 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t160 = t230 * V_base(4) - t231 * V_base(5) + V_base(3);
t158 = -rSges(6,3) * t251 + (rSges(6,1) * t245 - rSges(6,2) * t244) * t250;
t157 = -Icges(6,5) * t251 + (Icges(6,1) * t245 - Icges(6,4) * t244) * t250;
t156 = -Icges(6,6) * t251 + (Icges(6,4) * t245 - Icges(6,2) * t244) * t250;
t155 = -Icges(6,3) * t251 + (Icges(6,5) * t245 - Icges(6,6) * t244) * t250;
t154 = t264 * t297 + t184;
t153 = t262 * t297 + t183;
t152 = -rSges(7,3) * t251 + (rSges(7,1) * t241 - rSges(7,2) * t240) * t250;
t150 = -Icges(7,5) * t251 + (Icges(7,1) * t241 - Icges(7,4) * t240) * t250;
t149 = -Icges(7,6) * t251 + (Icges(7,4) * t241 - Icges(7,2) * t240) * t250;
t148 = -Icges(7,3) * t251 + (Icges(7,5) * t241 - Icges(7,6) * t240) * t250;
t144 = -pkin(10) * t251 + t250 * t301;
t143 = rSges(5,1) * t203 + rSges(5,2) * t202 + rSges(5,3) * t310;
t142 = rSges(5,1) * t201 + rSges(5,2) * t200 + rSges(5,3) * t311;
t141 = Icges(5,1) * t203 + Icges(5,4) * t202 + Icges(5,5) * t310;
t140 = Icges(5,1) * t201 + Icges(5,4) * t200 + Icges(5,5) * t311;
t139 = Icges(5,4) * t203 + Icges(5,2) * t202 + Icges(5,6) * t310;
t138 = Icges(5,4) * t201 + Icges(5,2) * t200 + Icges(5,6) * t311;
t137 = Icges(5,5) * t203 + Icges(5,6) * t202 + Icges(5,3) * t310;
t136 = Icges(5,5) * t201 + Icges(5,6) * t200 + Icges(5,3) * t311;
t132 = rSges(6,1) * t188 + rSges(6,2) * t187 + rSges(6,3) * t310;
t131 = rSges(6,1) * t186 + rSges(6,2) * t185 + rSges(6,3) * t311;
t130 = Icges(6,1) * t188 + Icges(6,4) * t187 + Icges(6,5) * t310;
t129 = Icges(6,1) * t186 + Icges(6,4) * t185 + Icges(6,5) * t311;
t128 = Icges(6,4) * t188 + Icges(6,2) * t187 + Icges(6,6) * t310;
t127 = Icges(6,4) * t186 + Icges(6,2) * t185 + Icges(6,6) * t311;
t126 = Icges(6,5) * t188 + Icges(6,6) * t187 + Icges(6,3) * t310;
t125 = Icges(6,5) * t186 + Icges(6,6) * t185 + Icges(6,3) * t311;
t124 = rSges(7,1) * t171 + rSges(7,2) * t170 + rSges(7,3) * t310;
t123 = rSges(7,1) * t169 + rSges(7,2) * t168 + rSges(7,3) * t311;
t122 = Icges(7,1) * t171 + Icges(7,4) * t170 + Icges(7,5) * t310;
t121 = Icges(7,1) * t169 + Icges(7,4) * t168 + Icges(7,5) * t311;
t120 = Icges(7,4) * t171 + Icges(7,2) * t170 + Icges(7,6) * t310;
t119 = Icges(7,4) * t169 + Icges(7,2) * t168 + Icges(7,6) * t311;
t118 = Icges(7,5) * t171 + Icges(7,6) * t170 + Icges(7,3) * t310;
t117 = Icges(7,5) * t169 + Icges(7,6) * t168 + Icges(7,3) * t311;
t115 = t229 * t238 + (-t195 - t236) * t247 + t296;
t114 = t196 * t247 - t229 * t239 + t280;
t113 = t262 * t292 + t264 * t275;
t112 = t262 * t275 - t264 * t292;
t111 = t195 * t239 - t196 * t238 + t279;
t110 = t209 * t210 + (-t181 + t302) * t247 + t291;
t109 = t182 * t247 - t209 * t211 + t273;
t108 = t181 * t211 - t182 * t210 + t274;
t107 = t167 * t210 + (-t142 + t293) * t247 + t281;
t106 = t143 * t247 + (-t167 - t208) * t211 + t272;
t105 = t142 * t211 + (-t143 - t198) * t210 + t270;
t104 = -t131 * t219 + t158 * t183 + t271;
t103 = t132 * t219 - t158 * t184 + t269;
t102 = t131 * t184 - t132 * t183 + t268;
t101 = -t112 * t219 - t123 * t199 + t144 * t183 + t152 * t153 + t271;
t100 = t113 * t219 + t124 * t199 - t144 * t184 - t152 * t154 + t269;
t99 = t112 * t184 - t113 * t183 + t123 * t154 - t124 * t153 + t268;
t1 = t219 * ((-t125 * t183 - t126 * t184 - t155 * t219) * t251 + ((-t128 * t244 + t130 * t245) * t184 + (-t127 * t244 + t129 * t245) * t183 + (-t156 * t244 + t157 * t245) * t219) * t250) / 0.2e1 + m(2) * (t160 ^ 2 + t162 ^ 2 + t163 ^ 2) / 0.2e1 + m(3) * (t111 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(5) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(4) * (t108 ^ 2 + t109 ^ 2 + t110 ^ 2) / 0.2e1 + m(7) * (t100 ^ 2 + t101 ^ 2 + t99 ^ 2) / 0.2e1 + m(6) * (t102 ^ 2 + t103 ^ 2 + t104 ^ 2) / 0.2e1 + t239 * (t277 * t262 + t266 * t264) / 0.2e1 + t238 * (t266 * t262 - t277 * t264) / 0.2e1 + t153 * ((t118 * t311 + t120 * t168 + t122 * t169) * t154 + (t117 * t311 + t119 * t168 + t121 * t169) * t153 + (t148 * t311 + t149 * t168 + t150 * t169) * t199) / 0.2e1 + t184 * ((t126 * t310 + t128 * t187 + t130 * t188) * t184 + (t125 * t310 + t127 * t187 + t129 * t188) * t183 + (t155 * t310 + t156 * t187 + t157 * t188) * t219) / 0.2e1 + t154 * ((t118 * t310 + t120 * t170 + t122 * t171) * t154 + (t117 * t310 + t119 * t170 + t121 * t171) * t153 + (t148 * t310 + t149 * t170 + t150 * t171) * t199) / 0.2e1 + t183 * ((t126 * t311 + t128 * t185 + t130 * t186) * t184 + (t125 * t311 + t127 * t185 + t129 * t186) * t183 + (t155 * t311 + t156 * t185 + t157 * t186) * t219) / 0.2e1 + t199 * ((-t117 * t153 - t118 * t154 - t148 * t199) * t251 + ((-t120 * t240 + t122 * t241) * t154 + (-t119 * t240 + t121 * t241) * t153 + (-t149 * t240 + t150 * t241) * t199) * t250) / 0.2e1 + m(1) * (t216 ^ 2 + t217 ^ 2 + t218 ^ 2) / 0.2e1 + ((t137 * t311 + t139 * t200 + t141 * t201) * t211 + (t136 * t311 + t138 * t200 + t140 * t201) * t210 + (t164 * t311 + t165 * t200 + t166 * t201) * t247 + t262 * t267 - t264 * t278) * t210 / 0.2e1 + ((t137 * t310 + t139 * t202 + t141 * t203) * t211 + (t136 * t310 + t138 * t202 + t140 * t203) * t210 + (t164 * t310 + t165 * t202 + t166 * t203) * t247 + t262 * t278 + t264 * t267) * t211 / 0.2e1 + ((-t224 * t262 + t227 * t264 + Icges(1,4)) * V_base(5) + (-t225 * t262 + t228 * t264 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t224 * t264 + t227 * t262 + Icges(1,2)) * V_base(5) + (t225 * t264 + t228 * t262 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t136 * t210 - t137 * t211) * t251 + ((-t139 * t258 + t141 * t259) * t211 + (-t138 * t258 + t140 * t259) * t210) * t250 + (t192 * t263 + t194 * t261) * t239 + (t191 * t263 + t193 * t261) * t238 + (t178 * t251 + t180 * t250) * t211 + (t177 * t251 + t179 * t250) * t210 + (t223 * t263 + t226 * t261 + Icges(2,3) + (-t164 + t206) * t251 + (-t165 * t258 + t166 * t259 + t207) * t250) * t247) * t247 / 0.2e1 + t247 * V_base(4) * (Icges(2,5) * t264 - Icges(2,6) * t262) + V_base(5) * t247 * (Icges(2,5) * t262 + Icges(2,6) * t264) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
