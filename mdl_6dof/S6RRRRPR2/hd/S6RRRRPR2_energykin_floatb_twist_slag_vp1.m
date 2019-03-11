% Calculate kinetic energy for
% S6RRRRPR2
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
% Datum: 2019-03-09 22:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRRPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRRPR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:57:19
% EndTime: 2019-03-09 21:57:22
% DurationCPUTime: 2.74s
% Computational Cost: add. (2316->364), mult. (2054->540), div. (0->0), fcn. (1892->12), ass. (0->185)
t264 = sin(qJ(2));
t333 = pkin(2) * t264;
t260 = qJ(2) + qJ(3);
t252 = sin(t260);
t332 = pkin(3) * t252;
t266 = cos(qJ(2));
t331 = t266 * pkin(2);
t262 = cos(pkin(11));
t330 = pkin(5) * t262;
t265 = sin(qJ(1));
t328 = Icges(2,4) * t265;
t327 = Icges(3,4) * t264;
t326 = Icges(3,4) * t266;
t325 = Icges(4,4) * t252;
t253 = cos(t260);
t324 = Icges(4,4) * t253;
t255 = qJ(4) + t260;
t243 = sin(t255);
t323 = Icges(5,4) * t243;
t244 = cos(t255);
t322 = Icges(5,4) * t244;
t321 = t243 * t265;
t267 = cos(qJ(1));
t320 = t243 * t267;
t258 = pkin(11) + qJ(6);
t246 = sin(t258);
t319 = t246 * t265;
t318 = t246 * t267;
t247 = cos(t258);
t317 = t247 * t265;
t316 = t247 * t267;
t261 = sin(pkin(11));
t315 = t261 * t265;
t314 = t261 * t267;
t313 = t262 * t265;
t312 = t262 * t267;
t174 = -pkin(8) * t267 + t265 * t331;
t238 = t265 * pkin(1) - t267 * pkin(7);
t310 = -t174 - t238;
t309 = pkin(3) * t253;
t307 = qJD(5) * t243;
t306 = qJD(6) * t243;
t305 = -qJD(2) - qJD(3);
t304 = V_base(5) * pkin(6) + V_base(1);
t147 = -pkin(9) * t267 + t265 * t309;
t301 = -t147 + t310;
t241 = qJD(2) * t265 + V_base(4);
t248 = V_base(6) + qJD(1);
t240 = -qJD(2) * t267 + V_base(5);
t300 = t240 * t333 + t304;
t294 = pkin(4) * t244 + qJ(5) * t243;
t188 = t294 * t265;
t299 = -t188 + t301;
t216 = qJD(3) * t265 + t241;
t215 = t267 * t305 + V_base(5);
t298 = t215 * t332 + t300;
t297 = rSges(3,1) * t266 - rSges(3,2) * t264;
t296 = rSges(4,1) * t253 - rSges(4,2) * t252;
t295 = rSges(5,1) * t244 - rSges(5,2) * t243;
t204 = qJD(4) * t265 + t216;
t293 = Icges(3,1) * t266 - t327;
t292 = Icges(4,1) * t253 - t325;
t291 = Icges(5,1) * t244 - t323;
t290 = -Icges(3,2) * t264 + t326;
t289 = -Icges(4,2) * t252 + t324;
t288 = -Icges(5,2) * t243 + t322;
t287 = Icges(3,5) * t266 - Icges(3,6) * t264;
t286 = Icges(4,5) * t253 - Icges(4,6) * t252;
t285 = Icges(5,5) * t244 - Icges(5,6) * t243;
t239 = t267 * pkin(1) + t265 * pkin(7);
t284 = -V_base(4) * pkin(6) + t248 * t239 + V_base(2);
t283 = V_base(4) * t238 - t239 * V_base(5) + V_base(3);
t203 = V_base(5) + (-qJD(4) + t305) * t267;
t208 = pkin(4) * t243 - qJ(5) * t244;
t282 = t203 * t208 + t267 * t307 + t298;
t281 = (-Icges(5,3) * t267 + t265 * t285) * t203 + (Icges(5,3) * t265 + t267 * t285) * t204 + (Icges(5,5) * t243 + Icges(5,6) * t244) * t248;
t280 = (-Icges(4,3) * t267 + t265 * t286) * t215 + (Icges(4,3) * t265 + t267 * t286) * t216 + (Icges(4,5) * t252 + Icges(4,6) * t253) * t248;
t279 = (-Icges(3,3) * t267 + t265 * t287) * t240 + (Icges(3,3) * t265 + t267 * t287) * t241 + (Icges(3,5) * t264 + Icges(3,6) * t266) * t248;
t278 = pkin(10) * t243 + t244 * t330;
t175 = pkin(8) * t265 + t267 * t331;
t277 = t241 * t174 - t175 * t240 + t283;
t276 = t248 * t175 - t241 * t333 + t284;
t148 = pkin(9) * t265 + t267 * t309;
t275 = t216 * t147 - t148 * t215 + t277;
t274 = t248 * t148 - t216 * t332 + t276;
t189 = t294 * t267;
t273 = t248 * t189 + t265 * t307 + t274;
t272 = -qJD(5) * t244 + t204 * t188 + t275;
t168 = -Icges(5,6) * t267 + t265 * t288;
t169 = Icges(5,6) * t265 + t267 * t288;
t170 = -Icges(5,5) * t267 + t265 * t291;
t171 = Icges(5,5) * t265 + t267 * t291;
t206 = Icges(5,2) * t244 + t323;
t207 = Icges(5,1) * t243 + t322;
t271 = (-t169 * t243 + t171 * t244) * t204 + (-t168 * t243 + t170 * t244) * t203 + (-t206 * t243 + t207 * t244) * t248;
t178 = -Icges(4,6) * t267 + t265 * t289;
t179 = Icges(4,6) * t265 + t267 * t289;
t180 = -Icges(4,5) * t267 + t265 * t292;
t181 = Icges(4,5) * t265 + t267 * t292;
t212 = Icges(4,2) * t253 + t325;
t213 = Icges(4,1) * t252 + t324;
t270 = (-t179 * t252 + t181 * t253) * t216 + (-t178 * t252 + t180 * t253) * t215 + (-t212 * t252 + t213 * t253) * t248;
t192 = -Icges(3,6) * t267 + t265 * t290;
t193 = Icges(3,6) * t265 + t267 * t290;
t194 = -Icges(3,5) * t267 + t265 * t293;
t195 = Icges(3,5) * t265 + t267 * t293;
t229 = Icges(3,2) * t266 + t327;
t232 = Icges(3,1) * t264 + t326;
t269 = (-t193 * t264 + t195 * t266) * t241 + (-t192 * t264 + t194 * t266) * t240 + (-t229 * t264 + t232 * t266) * t248;
t254 = Icges(2,4) * t267;
t237 = rSges(2,1) * t267 - rSges(2,2) * t265;
t236 = rSges(2,1) * t265 + rSges(2,2) * t267;
t235 = rSges(3,1) * t264 + rSges(3,2) * t266;
t234 = Icges(2,1) * t267 - t328;
t233 = Icges(2,1) * t265 + t254;
t231 = -Icges(2,2) * t265 + t254;
t230 = Icges(2,2) * t267 + t328;
t223 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t222 = rSges(1,1) * V_base(6) - V_base(4) * rSges(1,3) + V_base(2);
t221 = -rSges(1,2) * V_base(6) + V_base(5) * rSges(1,3) + V_base(1);
t217 = -qJD(6) * t244 + t248;
t214 = rSges(4,1) * t252 + rSges(4,2) * t253;
t209 = rSges(5,1) * t243 + rSges(5,2) * t244;
t201 = t244 * t312 + t315;
t200 = -t244 * t314 + t313;
t199 = t244 * t313 - t314;
t198 = -t244 * t315 - t312;
t197 = rSges(3,3) * t265 + t267 * t297;
t196 = -rSges(3,3) * t267 + t265 * t297;
t187 = t244 * t316 + t319;
t186 = -t244 * t318 + t317;
t185 = t244 * t317 - t318;
t184 = -t244 * t319 - t316;
t183 = rSges(4,3) * t265 + t267 * t296;
t182 = -rSges(4,3) * t267 + t265 * t296;
t173 = rSges(5,3) * t265 + t267 * t295;
t172 = -rSges(5,3) * t267 + t265 * t295;
t165 = V_base(5) * rSges(2,3) - t236 * t248 + t304;
t164 = t237 * t248 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t162 = t267 * t306 + t204;
t161 = t265 * t306 + t203;
t160 = t236 * V_base(4) - t237 * V_base(5) + V_base(3);
t159 = -rSges(6,3) * t244 + (rSges(6,1) * t262 - rSges(6,2) * t261) * t243;
t157 = -Icges(6,5) * t244 + (Icges(6,1) * t262 - Icges(6,4) * t261) * t243;
t156 = -Icges(6,6) * t244 + (Icges(6,4) * t262 - Icges(6,2) * t261) * t243;
t155 = -Icges(6,3) * t244 + (Icges(6,5) * t262 - Icges(6,6) * t261) * t243;
t153 = -rSges(7,3) * t244 + (rSges(7,1) * t247 - rSges(7,2) * t246) * t243;
t152 = -Icges(7,5) * t244 + (Icges(7,1) * t247 - Icges(7,4) * t246) * t243;
t151 = -Icges(7,6) * t244 + (Icges(7,4) * t247 - Icges(7,2) * t246) * t243;
t150 = -Icges(7,3) * t244 + (Icges(7,5) * t247 - Icges(7,6) * t246) * t243;
t146 = -pkin(10) * t244 + t243 * t330;
t142 = rSges(6,1) * t201 + rSges(6,2) * t200 + rSges(6,3) * t320;
t141 = rSges(6,1) * t199 + rSges(6,2) * t198 + rSges(6,3) * t321;
t140 = Icges(6,1) * t201 + Icges(6,4) * t200 + Icges(6,5) * t320;
t139 = Icges(6,1) * t199 + Icges(6,4) * t198 + Icges(6,5) * t321;
t138 = Icges(6,4) * t201 + Icges(6,2) * t200 + Icges(6,6) * t320;
t137 = Icges(6,4) * t199 + Icges(6,2) * t198 + Icges(6,6) * t321;
t136 = Icges(6,5) * t201 + Icges(6,6) * t200 + Icges(6,3) * t320;
t135 = Icges(6,5) * t199 + Icges(6,6) * t198 + Icges(6,3) * t321;
t134 = pkin(5) * t315 + t267 * t278;
t133 = -pkin(5) * t314 + t265 * t278;
t132 = rSges(7,1) * t187 + rSges(7,2) * t186 + rSges(7,3) * t320;
t131 = rSges(7,1) * t185 + rSges(7,2) * t184 + rSges(7,3) * t321;
t130 = Icges(7,1) * t187 + Icges(7,4) * t186 + Icges(7,5) * t320;
t129 = Icges(7,1) * t185 + Icges(7,4) * t184 + Icges(7,5) * t321;
t128 = Icges(7,4) * t187 + Icges(7,2) * t186 + Icges(7,6) * t320;
t127 = Icges(7,4) * t185 + Icges(7,2) * t184 + Icges(7,6) * t321;
t126 = Icges(7,5) * t187 + Icges(7,6) * t186 + Icges(7,3) * t320;
t125 = Icges(7,5) * t185 + Icges(7,6) * t184 + Icges(7,3) * t321;
t124 = t235 * t240 + (-t196 - t238) * t248 + t304;
t123 = t197 * t248 - t235 * t241 + t284;
t122 = t196 * t241 - t197 * t240 + t283;
t121 = t214 * t215 + (-t182 + t310) * t248 + t300;
t120 = t183 * t248 - t214 * t216 + t276;
t119 = t182 * t216 - t183 * t215 + t277;
t118 = t203 * t209 + (-t172 + t301) * t248 + t298;
t117 = t173 * t248 - t204 * t209 + t274;
t116 = t172 * t204 - t173 * t203 + t275;
t115 = t159 * t203 + (-t141 + t299) * t248 + t282;
t114 = t142 * t248 + (-t159 - t208) * t204 + t273;
t113 = t141 * t204 + (-t142 - t189) * t203 + t272;
t112 = -t131 * t217 + t146 * t203 + t153 * t161 + (-t133 + t299) * t248 + t282;
t111 = t132 * t217 + t134 * t248 - t153 * t162 + (-t146 - t208) * t204 + t273;
t110 = t131 * t162 - t132 * t161 + t133 * t204 + (-t134 - t189) * t203 + t272;
t1 = t161 * ((t126 * t321 + t128 * t184 + t130 * t185) * t162 + (t125 * t321 + t184 * t127 + t185 * t129) * t161 + (t150 * t321 + t151 * t184 + t152 * t185) * t217) / 0.2e1 + t162 * ((t126 * t320 + t186 * t128 + t187 * t130) * t162 + (t125 * t320 + t127 * t186 + t129 * t187) * t161 + (t150 * t320 + t151 * t186 + t152 * t187) * t217) / 0.2e1 + t241 * (t279 * t265 + t269 * t267) / 0.2e1 + t240 * (t269 * t265 - t279 * t267) / 0.2e1 + t216 * (t280 * t265 + t270 * t267) / 0.2e1 + t215 * (t270 * t265 - t280 * t267) / 0.2e1 + t217 * ((-t125 * t161 - t126 * t162 - t150 * t217) * t244 + ((-t128 * t246 + t130 * t247) * t162 + (-t127 * t246 + t129 * t247) * t161 + (-t151 * t246 + t152 * t247) * t217) * t243) / 0.2e1 + m(1) * (t221 ^ 2 + t222 ^ 2 + t223 ^ 2) / 0.2e1 + m(2) * (t160 ^ 2 + t164 ^ 2 + t165 ^ 2) / 0.2e1 + m(3) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(5) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(4) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(7) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(6) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + ((t136 * t321 + t138 * t198 + t140 * t199) * t204 + (t135 * t321 + t198 * t137 + t199 * t139) * t203 + (t155 * t321 + t156 * t198 + t157 * t199) * t248 + t271 * t265 - t281 * t267) * t203 / 0.2e1 + ((t136 * t320 + t200 * t138 + t201 * t140) * t204 + (t135 * t320 + t137 * t200 + t139 * t201) * t203 + (t155 * t320 + t156 * t200 + t157 * t201) * t248 + t281 * t265 + t271 * t267) * t204 / 0.2e1 + ((-t230 * t265 + t233 * t267 + Icges(1,4)) * V_base(5) + (-t265 * t231 + t267 * t234 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t267 * t230 + t265 * t233 + Icges(1,2)) * V_base(5) + (t231 * t267 + t234 * t265 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((-t135 * t203 - t136 * t204) * t244 + ((-t138 * t261 + t140 * t262) * t204 + (-t137 * t261 + t139 * t262) * t203) * t243 + (t193 * t266 + t195 * t264) * t241 + (t192 * t266 + t194 * t264) * t240 + (t179 * t253 + t181 * t252) * t216 + (t178 * t253 + t180 * t252) * t215 + (t169 * t244 + t171 * t243) * t204 + (t168 * t244 + t170 * t243) * t203 + (t253 * t212 + t252 * t213 + t266 * t229 + t264 * t232 + Icges(2,3) + (-t155 + t206) * t244 + (-t156 * t261 + t157 * t262 + t207) * t243) * t248) * t248 / 0.2e1 + t248 * V_base(4) * (Icges(2,5) * t267 - Icges(2,6) * t265) + t248 * V_base(5) * (Icges(2,5) * t265 + Icges(2,6) * t267) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
