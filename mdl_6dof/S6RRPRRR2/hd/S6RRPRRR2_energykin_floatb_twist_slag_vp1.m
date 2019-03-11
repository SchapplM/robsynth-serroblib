% Calculate kinetic energy for
% S6RRPRRR2
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
% Datum: 2019-03-09 13:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR2_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRPRRR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:16:35
% EndTime: 2019-03-09 13:16:38
% DurationCPUTime: 3.35s
% Computational Cost: add. (2360->359), mult. (2050->529), div. (0->0), fcn. (1888->12), ass. (0->182)
t339 = Icges(3,3) + Icges(4,3);
t260 = qJ(2) + pkin(11);
t247 = sin(t260);
t248 = cos(t260);
t264 = sin(qJ(2));
t267 = cos(qJ(2));
t338 = Icges(3,5) * t267 + Icges(4,5) * t248 - Icges(3,6) * t264 - Icges(4,6) * t247;
t265 = sin(qJ(1));
t268 = cos(qJ(1));
t322 = Icges(4,4) * t248;
t291 = -Icges(4,2) * t247 + t322;
t178 = -Icges(4,6) * t268 + t265 * t291;
t179 = Icges(4,6) * t265 + t268 * t291;
t323 = Icges(4,4) * t247;
t294 = Icges(4,1) * t248 - t323;
t180 = -Icges(4,5) * t268 + t265 * t294;
t181 = Icges(4,5) * t265 + t268 * t294;
t324 = Icges(3,4) * t267;
t292 = -Icges(3,2) * t264 + t324;
t195 = -Icges(3,6) * t268 + t265 * t292;
t196 = Icges(3,6) * t265 + t268 * t292;
t325 = Icges(3,4) * t264;
t295 = Icges(3,1) * t267 - t325;
t197 = -Icges(3,5) * t268 + t265 * t295;
t198 = Icges(3,5) * t265 + t268 * t295;
t213 = Icges(4,2) * t248 + t323;
t214 = Icges(4,1) * t247 + t322;
t230 = Icges(3,2) * t267 + t325;
t233 = Icges(3,1) * t264 + t324;
t241 = -qJD(2) * t268 + V_base(5);
t242 = qJD(2) * t265 + V_base(4);
t250 = V_base(6) + qJD(1);
t337 = (-t213 * t247 + t214 * t248 - t230 * t264 + t233 * t267) * t250 + (-t179 * t247 + t181 * t248 - t196 * t264 + t198 * t267) * t242 + (-t178 * t247 + t180 * t248 - t195 * t264 + t197 * t267) * t241;
t336 = (Icges(3,5) * t264 + Icges(4,5) * t247 + Icges(3,6) * t267 + Icges(4,6) * t248) * t250 + (t265 * t339 + t338 * t268) * t242 + (t338 * t265 - t268 * t339) * t241;
t332 = pkin(2) * t264;
t331 = pkin(3) * t247;
t330 = t267 * pkin(2);
t266 = cos(qJ(5));
t329 = pkin(5) * t266;
t326 = Icges(2,4) * t265;
t249 = qJ(4) + t260;
t243 = sin(t249);
t321 = Icges(5,4) * t243;
t244 = cos(t249);
t320 = Icges(5,4) * t244;
t319 = t243 * t265;
t318 = t243 * t268;
t261 = qJ(5) + qJ(6);
t254 = sin(t261);
t317 = t254 * t265;
t316 = t254 * t268;
t255 = cos(t261);
t315 = t255 * t265;
t314 = t255 * t268;
t263 = sin(qJ(5));
t313 = t263 * t265;
t312 = t263 * t268;
t311 = t265 * t266;
t310 = t266 * t268;
t174 = -qJ(3) * t268 + t265 * t330;
t239 = pkin(1) * t265 - pkin(7) * t268;
t309 = -t174 - t239;
t308 = pkin(3) * t248;
t306 = qJD(5) * t243;
t305 = qJD(6) * t243;
t304 = V_base(5) * pkin(6) + V_base(1);
t145 = -pkin(8) * t268 + t265 * t308;
t301 = -t145 + t309;
t218 = qJD(4) * t265 + t242;
t300 = qJD(3) * t265 + t241 * t332 + t304;
t299 = pkin(4) * t244 + pkin(9) * t243;
t298 = rSges(3,1) * t267 - rSges(3,2) * t264;
t297 = rSges(4,1) * t248 - rSges(4,2) * t247;
t296 = rSges(5,1) * t244 - rSges(5,2) * t243;
t183 = t268 * t306 + t218;
t293 = Icges(5,1) * t244 - t321;
t290 = -Icges(5,2) * t243 + t320;
t287 = Icges(5,5) * t244 - Icges(5,6) * t243;
t286 = t241 * t331 + t300;
t240 = pkin(1) * t268 + pkin(7) * t265;
t285 = -V_base(4) * pkin(6) + t250 * t240 + V_base(2);
t284 = V_base(4) * t239 - t240 * V_base(5) + V_base(3);
t217 = V_base(5) + (-qJD(2) - qJD(4)) * t268;
t283 = t242 * t174 + t284;
t182 = t265 * t306 + t217;
t282 = pkin(10) * t243 + t244 * t329;
t281 = (-Icges(5,3) * t268 + t265 * t287) * t217 + (Icges(5,3) * t265 + t268 * t287) * t218 + (Icges(5,5) * t243 + Icges(5,6) * t244) * t250;
t175 = qJ(3) * t265 + t268 * t330;
t278 = -qJD(3) * t268 + t250 * t175 + t285;
t146 = pkin(8) * t265 + t268 * t308;
t277 = t242 * t145 + (-t146 - t175) * t241 + t283;
t190 = t299 * t265;
t209 = t243 * pkin(4) - t244 * pkin(9);
t276 = t217 * t209 + (-t190 + t301) * t250 + t286;
t191 = t299 * t268;
t275 = t218 * t190 - t191 * t217 + t277;
t274 = t250 * t146 + (-t331 - t332) * t242 + t278;
t273 = t250 * t191 - t209 * t218 + t274;
t166 = -Icges(5,6) * t268 + t265 * t290;
t167 = Icges(5,6) * t265 + t268 * t290;
t168 = -Icges(5,5) * t268 + t265 * t293;
t169 = Icges(5,5) * t265 + t268 * t293;
t206 = Icges(5,2) * t244 + t321;
t207 = Icges(5,1) * t243 + t320;
t272 = (-t167 * t243 + t169 * t244) * t218 + (-t166 * t243 + t168 * t244) * t217 + (-t206 * t243 + t207 * t244) * t250;
t256 = Icges(2,4) * t268;
t238 = rSges(2,1) * t268 - rSges(2,2) * t265;
t237 = rSges(2,1) * t265 + rSges(2,2) * t268;
t236 = rSges(3,1) * t264 + rSges(3,2) * t267;
t235 = Icges(2,1) * t268 - t326;
t234 = Icges(2,1) * t265 + t256;
t232 = -Icges(2,2) * t265 + t256;
t231 = Icges(2,2) * t268 + t326;
t224 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t223 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t222 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t216 = -qJD(5) * t244 + t250;
t215 = rSges(4,1) * t247 + rSges(4,2) * t248;
t208 = rSges(5,1) * t243 + rSges(5,2) * t244;
t204 = t244 * t310 + t313;
t203 = -t244 * t312 + t311;
t202 = t244 * t311 - t312;
t201 = -t244 * t313 - t310;
t200 = rSges(3,3) * t265 + t268 * t298;
t199 = -rSges(3,3) * t268 + t265 * t298;
t192 = (-qJD(5) - qJD(6)) * t244 + t250;
t189 = t244 * t314 + t317;
t188 = -t244 * t316 + t315;
t187 = t244 * t315 - t316;
t186 = -t244 * t317 - t314;
t185 = rSges(4,3) * t265 + t268 * t297;
t184 = -rSges(4,3) * t268 + t265 * t297;
t173 = rSges(5,3) * t265 + t268 * t296;
t172 = -rSges(5,3) * t268 + t265 * t296;
t171 = V_base(5) * rSges(2,3) - t237 * t250 + t304;
t170 = t238 * t250 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t162 = t237 * V_base(4) - t238 * V_base(5) + V_base(3);
t160 = -rSges(6,3) * t244 + (rSges(6,1) * t266 - rSges(6,2) * t263) * t243;
t158 = -Icges(6,5) * t244 + (Icges(6,1) * t266 - Icges(6,4) * t263) * t243;
t157 = -Icges(6,6) * t244 + (Icges(6,4) * t266 - Icges(6,2) * t263) * t243;
t156 = -Icges(6,3) * t244 + (Icges(6,5) * t266 - Icges(6,6) * t263) * t243;
t154 = -rSges(7,3) * t244 + (rSges(7,1) * t255 - rSges(7,2) * t254) * t243;
t153 = -Icges(7,5) * t244 + (Icges(7,1) * t255 - Icges(7,4) * t254) * t243;
t152 = -Icges(7,6) * t244 + (Icges(7,4) * t255 - Icges(7,2) * t254) * t243;
t151 = -Icges(7,3) * t244 + (Icges(7,5) * t255 - Icges(7,6) * t254) * t243;
t150 = t268 * t305 + t183;
t149 = t265 * t305 + t182;
t147 = -pkin(10) * t244 + t243 * t329;
t142 = rSges(6,1) * t204 + rSges(6,2) * t203 + rSges(6,3) * t318;
t141 = rSges(6,1) * t202 + rSges(6,2) * t201 + rSges(6,3) * t319;
t140 = Icges(6,1) * t204 + Icges(6,4) * t203 + Icges(6,5) * t318;
t139 = Icges(6,1) * t202 + Icges(6,4) * t201 + Icges(6,5) * t319;
t138 = Icges(6,4) * t204 + Icges(6,2) * t203 + Icges(6,6) * t318;
t137 = Icges(6,4) * t202 + Icges(6,2) * t201 + Icges(6,6) * t319;
t136 = Icges(6,5) * t204 + Icges(6,6) * t203 + Icges(6,3) * t318;
t135 = Icges(6,5) * t202 + Icges(6,6) * t201 + Icges(6,3) * t319;
t134 = pkin(5) * t313 + t268 * t282;
t133 = -pkin(5) * t312 + t265 * t282;
t132 = rSges(7,1) * t189 + rSges(7,2) * t188 + rSges(7,3) * t318;
t131 = rSges(7,1) * t187 + rSges(7,2) * t186 + rSges(7,3) * t319;
t130 = Icges(7,1) * t189 + Icges(7,4) * t188 + Icges(7,5) * t318;
t129 = Icges(7,1) * t187 + Icges(7,4) * t186 + Icges(7,5) * t319;
t128 = Icges(7,4) * t189 + Icges(7,2) * t188 + Icges(7,6) * t318;
t127 = Icges(7,4) * t187 + Icges(7,2) * t186 + Icges(7,6) * t319;
t126 = Icges(7,5) * t189 + Icges(7,6) * t188 + Icges(7,3) * t318;
t125 = Icges(7,5) * t187 + Icges(7,6) * t186 + Icges(7,3) * t319;
t124 = t236 * t241 + (-t199 - t239) * t250 + t304;
t123 = t200 * t250 - t236 * t242 + t285;
t122 = t199 * t242 - t200 * t241 + t284;
t121 = t215 * t241 + (-t184 + t309) * t250 + t300;
t120 = t185 * t250 + (-t215 - t332) * t242 + t278;
t119 = t184 * t242 + (-t175 - t185) * t241 + t283;
t118 = t208 * t217 + (-t172 + t301) * t250 + t286;
t117 = t173 * t250 - t208 * t218 + t274;
t116 = t172 * t218 - t173 * t217 + t277;
t115 = -t141 * t216 + t160 * t182 + t276;
t114 = t142 * t216 - t160 * t183 + t273;
t113 = t141 * t183 - t142 * t182 + t275;
t112 = -t131 * t192 - t133 * t216 + t147 * t182 + t149 * t154 + t276;
t111 = t132 * t192 + t134 * t216 - t147 * t183 - t150 * t154 + t273;
t110 = t131 * t150 - t132 * t149 + t133 * t183 - t134 * t182 + t275;
t1 = t182 * ((t136 * t319 + t138 * t201 + t140 * t202) * t183 + (t135 * t319 + t201 * t137 + t202 * t139) * t182 + (t156 * t319 + t157 * t201 + t158 * t202) * t216) / 0.2e1 + t149 * ((t126 * t319 + t128 * t186 + t130 * t187) * t150 + (t125 * t319 + t186 * t127 + t187 * t129) * t149 + (t151 * t319 + t152 * t186 + t153 * t187) * t192) / 0.2e1 + t183 * ((t136 * t318 + t203 * t138 + t204 * t140) * t183 + (t135 * t318 + t137 * t203 + t139 * t204) * t182 + (t156 * t318 + t157 * t203 + t158 * t204) * t216) / 0.2e1 + t150 * ((t126 * t318 + t188 * t128 + t189 * t130) * t150 + (t125 * t318 + t127 * t188 + t129 * t189) * t149 + (t151 * t318 + t152 * t188 + t153 * t189) * t192) / 0.2e1 + t218 * (t281 * t265 + t272 * t268) / 0.2e1 + t217 * (t272 * t265 - t281 * t268) / 0.2e1 + t216 * ((-t135 * t182 - t136 * t183 - t156 * t216) * t244 + ((-t138 * t263 + t140 * t266) * t183 + (-t137 * t263 + t139 * t266) * t182 + (-t157 * t263 + t158 * t266) * t216) * t243) / 0.2e1 + t192 * ((-t125 * t149 - t126 * t150 - t151 * t192) * t244 + ((-t128 * t254 + t130 * t255) * t150 + (-t127 * t254 + t129 * t255) * t149 + (-t152 * t254 + t153 * t255) * t192) * t243) / 0.2e1 + m(1) * (t222 ^ 2 + t223 ^ 2 + t224 ^ 2) / 0.2e1 + m(2) * (t162 ^ 2 + t170 ^ 2 + t171 ^ 2) / 0.2e1 + m(5) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + m(4) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(3) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(7) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(6) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + (t337 * t265 - t336 * t268) * t241 / 0.2e1 + (t336 * t265 + t337 * t268) * t242 / 0.2e1 + ((-t231 * t265 + t234 * t268 + Icges(1,4)) * V_base(5) + (-t265 * t232 + t268 * t235 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t268 * t231 + t265 * t234 + Icges(1,2)) * V_base(5) + (t232 * t268 + t235 * t265 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t167 * t244 + t169 * t243) * t218 + (t166 * t244 + t168 * t243) * t217 + (t179 * t248 + t181 * t247 + t196 * t267 + t198 * t264) * t242 + (t178 * t248 + t180 * t247 + t195 * t267 + t197 * t264) * t241 + (t244 * t206 + t243 * t207 + t248 * t213 + t247 * t214 + t267 * t230 + t264 * t233 + Icges(2,3)) * t250) * t250 / 0.2e1 + t250 * V_base(4) * (Icges(2,5) * t268 - Icges(2,6) * t265) + t250 * V_base(5) * (Icges(2,5) * t265 + Icges(2,6) * t268) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
