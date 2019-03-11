% Calculate kinetic energy for
% S6RPRRPR4
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
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR4_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR4_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRRPR4_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:08:59
% EndTime: 2019-03-09 05:09:01
% DurationCPUTime: 2.96s
% Computational Cost: add. (2226->367), mult. (1964->532), div. (0->0), fcn. (1802->12), ass. (0->185)
t260 = sin(pkin(10));
t331 = pkin(2) * t260;
t258 = pkin(10) + qJ(3);
t245 = sin(t258);
t330 = pkin(3) * t245;
t262 = cos(pkin(10));
t329 = t262 * pkin(2);
t261 = cos(pkin(11));
t328 = pkin(5) * t261;
t265 = sin(qJ(1));
t327 = Icges(2,4) * t265;
t326 = Icges(3,4) * t260;
t325 = Icges(3,4) * t262;
t324 = Icges(4,4) * t245;
t247 = cos(t258);
t323 = Icges(4,4) * t247;
t248 = qJ(4) + t258;
t240 = sin(t248);
t322 = Icges(5,4) * t240;
t241 = cos(t248);
t321 = Icges(5,4) * t241;
t320 = t240 * t265;
t266 = cos(qJ(1));
t319 = t240 * t266;
t257 = pkin(11) + qJ(6);
t244 = sin(t257);
t318 = t244 * t266;
t246 = cos(t257);
t317 = t246 * t266;
t259 = sin(pkin(11));
t316 = t259 * t266;
t315 = t261 * t266;
t314 = t265 * t244;
t313 = t265 * t246;
t312 = t265 * t259;
t311 = t265 * t261;
t172 = -pkin(7) * t266 + t265 * t329;
t233 = t265 * pkin(1) - qJ(2) * t266;
t308 = -t172 - t233;
t307 = pkin(3) * t247;
t305 = qJD(5) * t240;
t304 = qJD(6) * t240;
t303 = V_base(4) * t233 + V_base(3);
t302 = V_base(5) * pkin(6) + V_base(1);
t145 = -pkin(8) * t266 + t265 * t307;
t299 = -t145 + t308;
t238 = qJD(3) * t265 + V_base(4);
t249 = V_base(6) + qJD(1);
t298 = qJD(2) * t265 + t302;
t292 = pkin(4) * t241 + qJ(5) * t240;
t188 = t292 * t265;
t297 = -t188 + t299;
t215 = qJD(4) * t265 + t238;
t296 = V_base(5) * t331 + t298;
t295 = rSges(3,1) * t262 - rSges(3,2) * t260;
t294 = rSges(4,1) * t247 - rSges(4,2) * t245;
t293 = rSges(5,1) * t241 - rSges(5,2) * t240;
t291 = Icges(3,1) * t262 - t326;
t290 = Icges(4,1) * t247 - t324;
t289 = Icges(5,1) * t241 - t322;
t288 = -Icges(3,2) * t260 + t325;
t287 = -Icges(4,2) * t245 + t323;
t286 = -Icges(5,2) * t240 + t321;
t285 = Icges(3,5) * t262 - Icges(3,6) * t260;
t284 = Icges(4,5) * t247 - Icges(4,6) * t245;
t283 = Icges(5,5) * t241 - Icges(5,6) * t240;
t235 = pkin(1) * t266 + t265 * qJ(2);
t282 = -qJD(2) * t266 + t249 * t235 + V_base(2);
t237 = -qJD(3) * t266 + V_base(5);
t281 = t237 * t330 + t296;
t214 = V_base(5) + (-qJD(3) - qJD(4)) * t266;
t205 = pkin(4) * t240 - qJ(5) * t241;
t280 = t214 * t205 + t266 * t305 + t281;
t279 = (-Icges(5,3) * t266 + t265 * t283) * t214 + (Icges(5,3) * t265 + t266 * t283) * t215 + (Icges(5,5) * t240 + Icges(5,6) * t241) * t249;
t278 = (-Icges(4,3) * t266 + t265 * t284) * t237 + (Icges(4,3) * t265 + t266 * t284) * t238 + (Icges(4,5) * t245 + Icges(4,6) * t247) * t249;
t277 = pkin(9) * t240 + t241 * t328;
t173 = pkin(7) * t265 + t266 * t329;
t276 = V_base(4) * t172 + (-t173 - t235) * V_base(5) + t303;
t275 = (-Icges(3,3) * t266 + t265 * t285) * V_base(5) + (Icges(3,3) * t265 + t266 * t285) * V_base(4) + (Icges(3,5) * t260 + Icges(3,6) * t262) * t249;
t146 = pkin(8) * t265 + t266 * t307;
t274 = t238 * t145 - t146 * t237 + t276;
t273 = t249 * t173 + (-pkin(6) - t331) * V_base(4) + t282;
t272 = -qJD(5) * t241 + t215 * t188 + t274;
t271 = t249 * t146 - t238 * t330 + t273;
t189 = t292 * t266;
t270 = t249 * t189 + t265 * t305 + t271;
t164 = -Icges(5,6) * t266 + t265 * t286;
t165 = Icges(5,6) * t265 + t266 * t286;
t166 = -Icges(5,5) * t266 + t265 * t289;
t167 = Icges(5,5) * t265 + t266 * t289;
t203 = Icges(5,2) * t241 + t322;
t204 = Icges(5,1) * t240 + t321;
t269 = (-t165 * t240 + t167 * t241) * t215 + (-t164 * t240 + t166 * t241) * t214 + (-t203 * t240 + t204 * t241) * t249;
t176 = -Icges(4,6) * t266 + t265 * t287;
t177 = Icges(4,6) * t265 + t266 * t287;
t178 = -Icges(4,5) * t266 + t265 * t290;
t179 = Icges(4,5) * t265 + t266 * t290;
t210 = Icges(4,2) * t247 + t324;
t211 = Icges(4,1) * t245 + t323;
t268 = (-t177 * t245 + t179 * t247) * t238 + (-t176 * t245 + t178 * t247) * t237 + (-t210 * t245 + t211 * t247) * t249;
t192 = -Icges(3,6) * t266 + t265 * t288;
t193 = Icges(3,6) * t265 + t266 * t288;
t194 = -Icges(3,5) * t266 + t265 * t291;
t195 = Icges(3,5) * t265 + t266 * t291;
t224 = Icges(3,2) * t262 + t326;
t225 = Icges(3,1) * t260 + t325;
t267 = (-t193 * t260 + t195 * t262) * V_base(4) + (-t192 * t260 + t194 * t262) * V_base(5) + (-t224 * t260 + t225 * t262) * t249;
t254 = Icges(2,4) * t266;
t236 = rSges(2,1) * t266 - t265 * rSges(2,2);
t234 = t265 * rSges(2,1) + rSges(2,2) * t266;
t232 = Icges(2,1) * t266 - t327;
t231 = Icges(2,1) * t265 + t254;
t230 = -Icges(2,2) * t265 + t254;
t229 = Icges(2,2) * t266 + t327;
t228 = Icges(2,5) * t266 - Icges(2,6) * t265;
t227 = Icges(2,5) * t265 + Icges(2,6) * t266;
t226 = rSges(3,1) * t260 + rSges(3,2) * t262;
t220 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t219 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t218 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t213 = -qJD(6) * t241 + t249;
t212 = rSges(4,1) * t245 + rSges(4,2) * t247;
t206 = rSges(5,1) * t240 + rSges(5,2) * t241;
t201 = t241 * t315 + t312;
t200 = -t241 * t316 + t311;
t199 = t241 * t311 - t316;
t198 = -t241 * t312 - t315;
t197 = t265 * rSges(3,3) + t266 * t295;
t196 = -rSges(3,3) * t266 + t265 * t295;
t187 = t241 * t317 + t314;
t186 = -t241 * t318 + t313;
t185 = t241 * t313 - t318;
t184 = -t241 * t314 - t317;
t183 = t265 * rSges(4,3) + t266 * t294;
t182 = -rSges(4,3) * t266 + t265 * t294;
t181 = t266 * t304 + t215;
t180 = t265 * t304 + t214;
t171 = t265 * rSges(5,3) + t266 * t293;
t170 = -rSges(5,3) * t266 + t265 * t293;
t169 = V_base(5) * rSges(2,3) - t234 * t249 + t302;
t168 = t236 * t249 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t160 = t234 * V_base(4) - t236 * V_base(5) + V_base(3);
t157 = -rSges(6,3) * t241 + (rSges(6,1) * t261 - rSges(6,2) * t259) * t240;
t155 = -Icges(6,5) * t241 + (Icges(6,1) * t261 - Icges(6,4) * t259) * t240;
t154 = -Icges(6,6) * t241 + (Icges(6,4) * t261 - Icges(6,2) * t259) * t240;
t153 = -Icges(6,3) * t241 + (Icges(6,5) * t261 - Icges(6,6) * t259) * t240;
t152 = -rSges(7,3) * t241 + (rSges(7,1) * t246 - rSges(7,2) * t244) * t240;
t151 = -Icges(7,5) * t241 + (Icges(7,1) * t246 - Icges(7,4) * t244) * t240;
t150 = -Icges(7,6) * t241 + (Icges(7,4) * t246 - Icges(7,2) * t244) * t240;
t149 = -Icges(7,3) * t241 + (Icges(7,5) * t246 - Icges(7,6) * t244) * t240;
t147 = -pkin(9) * t241 + t240 * t328;
t142 = t201 * rSges(6,1) + t200 * rSges(6,2) + rSges(6,3) * t319;
t141 = rSges(6,1) * t199 + rSges(6,2) * t198 + rSges(6,3) * t320;
t140 = Icges(6,1) * t201 + Icges(6,4) * t200 + Icges(6,5) * t319;
t139 = Icges(6,1) * t199 + Icges(6,4) * t198 + Icges(6,5) * t320;
t138 = Icges(6,4) * t201 + Icges(6,2) * t200 + Icges(6,6) * t319;
t137 = Icges(6,4) * t199 + Icges(6,2) * t198 + Icges(6,6) * t320;
t136 = Icges(6,5) * t201 + Icges(6,6) * t200 + Icges(6,3) * t319;
t135 = Icges(6,5) * t199 + Icges(6,6) * t198 + Icges(6,3) * t320;
t134 = pkin(5) * t312 + t266 * t277;
t133 = -pkin(5) * t316 + t265 * t277;
t132 = t187 * rSges(7,1) + t186 * rSges(7,2) + rSges(7,3) * t319;
t131 = rSges(7,1) * t185 + rSges(7,2) * t184 + rSges(7,3) * t320;
t130 = Icges(7,1) * t187 + Icges(7,4) * t186 + Icges(7,5) * t319;
t129 = Icges(7,1) * t185 + Icges(7,4) * t184 + Icges(7,5) * t320;
t128 = Icges(7,4) * t187 + Icges(7,2) * t186 + Icges(7,6) * t319;
t127 = Icges(7,4) * t185 + Icges(7,2) * t184 + Icges(7,6) * t320;
t126 = Icges(7,5) * t187 + Icges(7,6) * t186 + Icges(7,3) * t319;
t125 = Icges(7,5) * t185 + Icges(7,6) * t184 + Icges(7,3) * t320;
t124 = t226 * V_base(5) + (-t196 - t233) * t249 + t298;
t123 = t249 * t197 + (-pkin(6) - t226) * V_base(4) + t282;
t122 = t196 * V_base(4) + (-t197 - t235) * V_base(5) + t303;
t121 = t212 * t237 + (-t182 + t308) * t249 + t296;
t120 = t249 * t183 - t238 * t212 + t273;
t119 = t182 * t238 - t183 * t237 + t276;
t118 = t206 * t214 + (-t170 + t299) * t249 + t281;
t117 = t249 * t171 - t215 * t206 + t271;
t116 = t170 * t215 - t171 * t214 + t274;
t115 = t157 * t214 + (-t141 + t297) * t249 + t280;
t114 = t249 * t142 + (-t157 - t205) * t215 + t270;
t113 = t141 * t215 + (-t142 - t189) * t214 + t272;
t112 = -t131 * t213 + t147 * t214 + t152 * t180 + (-t133 + t297) * t249 + t280;
t111 = t213 * t132 + t249 * t134 - t181 * t152 + (-t147 - t205) * t215 + t270;
t110 = t131 * t181 - t132 * t180 + t133 * t215 + (-t134 - t189) * t214 + t272;
t1 = t181 * ((t126 * t319 + t186 * t128 + t187 * t130) * t181 + (t125 * t319 + t186 * t127 + t187 * t129) * t180 + (t149 * t319 + t186 * t150 + t187 * t151) * t213) / 0.2e1 + t180 * ((t126 * t320 + t128 * t184 + t130 * t185) * t181 + (t125 * t320 + t184 * t127 + t185 * t129) * t180 + (t149 * t320 + t150 * t184 + t151 * t185) * t213) / 0.2e1 + t238 * (t265 * t278 + t266 * t268) / 0.2e1 + t237 * (t265 * t268 - t266 * t278) / 0.2e1 + m(2) * (t160 ^ 2 + t168 ^ 2 + t169 ^ 2) / 0.2e1 + m(3) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(6) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(5) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(4) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(7) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(1) * (t218 ^ 2 + t219 ^ 2 + t220 ^ 2) / 0.2e1 + t213 * ((-t125 * t180 - t126 * t181 - t149 * t213) * t241 + ((-t128 * t244 + t130 * t246) * t181 + (-t127 * t244 + t129 * t246) * t180 + (-t150 * t244 + t151 * t246) * t213) * t240) / 0.2e1 + ((t136 * t320 + t138 * t198 + t140 * t199) * t215 + (t135 * t320 + t198 * t137 + t199 * t139) * t214 + (t153 * t320 + t154 * t198 + t155 * t199) * t249 + t265 * t269 - t266 * t279) * t214 / 0.2e1 + ((t136 * t319 + t200 * t138 + t201 * t140) * t215 + (t135 * t319 + t200 * t137 + t201 * t139) * t214 + (t153 * t319 + t200 * t154 + t201 * t155) * t249 + t265 * t279 + t266 * t269) * t215 / 0.2e1 + (t228 * t249 + t265 * t275 + t266 * t267 + (-t265 * t229 + t231 * t266 + Icges(1,4)) * V_base(5) + (-t265 * t230 + t266 * t232 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t227 * t249 + t265 * t267 - t266 * t275 + (t266 * t229 + t265 * t231 + Icges(1,2)) * V_base(5) + (t230 * t266 + t265 * t232 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t135 * t214 - t136 * t215) * t241 + ((-t138 * t259 + t140 * t261) * t215 + (-t137 * t259 + t139 * t261) * t214) * t240 + (t177 * t247 + t179 * t245) * t238 + (t176 * t247 + t178 * t245) * t237 + (t165 * t241 + t167 * t240) * t215 + (t164 * t241 + t166 * t240) * t214 + (t192 * t262 + t194 * t260 + t227) * V_base(5) + (t193 * t262 + t195 * t260 + t228) * V_base(4) + (t247 * t210 + t245 * t211 + t262 * t224 + t260 * t225 + Icges(2,3) + (-t153 + t203) * t241 + (-t154 * t259 + t155 * t261 + t204) * t240) * t249) * t249 / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
