% Calculate kinetic energy for
% S6RRRPRR2
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
% Datum: 2019-03-09 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR2_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RRRPRR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:06:12
% EndTime: 2019-03-09 18:06:15
% DurationCPUTime: 3.31s
% Computational Cost: add. (2382->359), mult. (2072->529), div. (0->0), fcn. (1910->12), ass. (0->182)
t339 = Icges(4,3) + Icges(5,3);
t261 = qJ(2) + qJ(3);
t247 = pkin(11) + t261;
t243 = sin(t247);
t244 = cos(t247);
t253 = sin(t261);
t255 = cos(t261);
t338 = Icges(4,5) * t255 + Icges(5,5) * t244 - Icges(4,6) * t253 - Icges(5,6) * t243;
t264 = sin(qJ(1));
t267 = cos(qJ(1));
t320 = Icges(5,4) * t244;
t290 = -Icges(5,2) * t243 + t320;
t166 = -Icges(5,6) * t267 + t264 * t290;
t167 = Icges(5,6) * t264 + t267 * t290;
t321 = Icges(5,4) * t243;
t293 = Icges(5,1) * t244 - t321;
t168 = -Icges(5,5) * t267 + t264 * t293;
t169 = Icges(5,5) * t264 + t267 * t293;
t322 = Icges(4,4) * t255;
t291 = -Icges(4,2) * t253 + t322;
t180 = -Icges(4,6) * t267 + t264 * t291;
t181 = Icges(4,6) * t264 + t267 * t291;
t323 = Icges(4,4) * t253;
t294 = Icges(4,1) * t255 - t323;
t182 = -Icges(4,5) * t267 + t264 * t294;
t183 = Icges(4,5) * t264 + t267 * t294;
t207 = Icges(5,2) * t244 + t321;
t208 = Icges(5,1) * t243 + t320;
t213 = Icges(4,2) * t255 + t323;
t214 = Icges(4,1) * t253 + t322;
t217 = V_base(5) + (-qJD(2) - qJD(3)) * t267;
t242 = qJD(2) * t264 + V_base(4);
t218 = qJD(3) * t264 + t242;
t248 = V_base(6) + qJD(1);
t337 = (-t207 * t243 + t208 * t244 - t213 * t253 + t214 * t255) * t248 + (-t167 * t243 + t169 * t244 - t181 * t253 + t183 * t255) * t218 + (-t166 * t243 + t168 * t244 - t180 * t253 + t182 * t255) * t217;
t336 = (Icges(4,5) * t253 + Icges(5,5) * t243 + Icges(4,6) * t255 + Icges(5,6) * t244) * t248 + (t264 * t339 + t338 * t267) * t218 + (t338 * t264 - t267 * t339) * t217;
t263 = sin(qJ(2));
t332 = pkin(2) * t263;
t331 = pkin(3) * t253;
t266 = cos(qJ(2));
t330 = t266 * pkin(2);
t265 = cos(qJ(5));
t329 = pkin(5) * t265;
t326 = Icges(2,4) * t264;
t325 = Icges(3,4) * t263;
t324 = Icges(3,4) * t266;
t319 = t243 * t264;
t318 = t243 * t267;
t260 = qJ(5) + qJ(6);
t252 = sin(t260);
t317 = t252 * t264;
t316 = t252 * t267;
t254 = cos(t260);
t315 = t254 * t264;
t314 = t254 * t267;
t262 = sin(qJ(5));
t313 = t262 * t264;
t312 = t262 * t267;
t311 = t264 * t265;
t310 = t265 * t267;
t174 = -pkin(8) * t267 + t264 * t330;
t239 = t264 * pkin(1) - t267 * pkin(7);
t309 = -t174 - t239;
t308 = pkin(3) * t255;
t306 = qJD(5) * t243;
t305 = qJD(6) * t243;
t304 = V_base(5) * pkin(6) + V_base(1);
t146 = -qJ(4) * t267 + t264 * t308;
t301 = -t146 + t309;
t241 = -qJD(2) * t267 + V_base(5);
t300 = t241 * t332 + t304;
t299 = pkin(4) * t244 + pkin(9) * t243;
t298 = rSges(3,1) * t266 - rSges(3,2) * t263;
t297 = rSges(4,1) * t255 - rSges(4,2) * t253;
t296 = rSges(5,1) * t244 - rSges(5,2) * t243;
t177 = t267 * t306 + t218;
t295 = Icges(3,1) * t266 - t325;
t292 = -Icges(3,2) * t263 + t324;
t289 = Icges(3,5) * t266 - Icges(3,6) * t263;
t286 = qJD(4) * t264 + t217 * t331 + t300;
t240 = t267 * pkin(1) + t264 * pkin(7);
t285 = -V_base(4) * pkin(6) + t248 * t240 + V_base(2);
t284 = V_base(4) * t239 - t240 * V_base(5) + V_base(3);
t176 = t264 * t306 + t217;
t283 = pkin(10) * t243 + t244 * t329;
t280 = (-Icges(3,3) * t267 + t264 * t289) * t241 + (Icges(3,3) * t264 + t267 * t289) * t242 + (Icges(3,5) * t263 + Icges(3,6) * t266) * t248;
t175 = pkin(8) * t264 + t267 * t330;
t279 = t242 * t174 - t175 * t241 + t284;
t278 = t248 * t175 - t242 * t332 + t285;
t277 = t218 * t146 + t279;
t190 = t299 * t264;
t210 = pkin(4) * t243 - pkin(9) * t244;
t276 = t217 * t210 + (-t190 + t301) * t248 + t286;
t147 = qJ(4) * t264 + t267 * t308;
t275 = -qJD(4) * t267 + t248 * t147 + t278;
t191 = t299 * t267;
t274 = t218 * t190 + (-t147 - t191) * t217 + t277;
t273 = t248 * t191 + (-t210 - t331) * t218 + t275;
t195 = -Icges(3,6) * t267 + t264 * t292;
t196 = Icges(3,6) * t264 + t267 * t292;
t197 = -Icges(3,5) * t267 + t264 * t295;
t198 = Icges(3,5) * t264 + t267 * t295;
t230 = Icges(3,2) * t266 + t325;
t233 = Icges(3,1) * t263 + t324;
t270 = (-t196 * t263 + t198 * t266) * t242 + (-t195 * t263 + t197 * t266) * t241 + (-t230 * t263 + t233 * t266) * t248;
t256 = Icges(2,4) * t267;
t238 = rSges(2,1) * t267 - rSges(2,2) * t264;
t237 = rSges(2,1) * t264 + rSges(2,2) * t267;
t236 = rSges(3,1) * t263 + rSges(3,2) * t266;
t235 = Icges(2,1) * t267 - t326;
t234 = Icges(2,1) * t264 + t256;
t232 = -Icges(2,2) * t264 + t256;
t231 = Icges(2,2) * t267 + t326;
t224 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t223 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t222 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t216 = -qJD(5) * t244 + t248;
t215 = rSges(4,1) * t253 + rSges(4,2) * t255;
t209 = rSges(5,1) * t243 + rSges(5,2) * t244;
t205 = t244 * t310 + t313;
t204 = -t244 * t312 + t311;
t203 = t244 * t311 - t312;
t202 = -t244 * t313 - t310;
t200 = rSges(3,3) * t264 + t267 * t298;
t199 = -rSges(3,3) * t267 + t264 * t298;
t192 = (-qJD(5) - qJD(6)) * t244 + t248;
t189 = t244 * t314 + t317;
t188 = -t244 * t316 + t315;
t187 = t244 * t315 - t316;
t186 = -t244 * t317 - t314;
t185 = rSges(4,3) * t264 + t267 * t297;
t184 = -rSges(4,3) * t267 + t264 * t297;
t173 = rSges(5,3) * t264 + t267 * t296;
t172 = -rSges(5,3) * t267 + t264 * t296;
t171 = V_base(5) * rSges(2,3) - t237 * t248 + t304;
t170 = t238 * t248 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t162 = t237 * V_base(4) - t238 * V_base(5) + V_base(3);
t159 = -rSges(6,3) * t244 + (rSges(6,1) * t265 - rSges(6,2) * t262) * t243;
t158 = -Icges(6,5) * t244 + (Icges(6,1) * t265 - Icges(6,4) * t262) * t243;
t157 = -Icges(6,6) * t244 + (Icges(6,4) * t265 - Icges(6,2) * t262) * t243;
t156 = -Icges(6,3) * t244 + (Icges(6,5) * t265 - Icges(6,6) * t262) * t243;
t154 = -rSges(7,3) * t244 + (rSges(7,1) * t254 - rSges(7,2) * t252) * t243;
t153 = -Icges(7,5) * t244 + (Icges(7,1) * t254 - Icges(7,4) * t252) * t243;
t152 = -Icges(7,6) * t244 + (Icges(7,4) * t254 - Icges(7,2) * t252) * t243;
t151 = -Icges(7,3) * t244 + (Icges(7,5) * t254 - Icges(7,6) * t252) * t243;
t150 = t267 * t305 + t177;
t149 = t264 * t305 + t176;
t145 = -pkin(10) * t244 + t243 * t329;
t142 = rSges(6,1) * t205 + rSges(6,2) * t204 + rSges(6,3) * t318;
t141 = rSges(6,1) * t203 + rSges(6,2) * t202 + rSges(6,3) * t319;
t140 = Icges(6,1) * t205 + Icges(6,4) * t204 + Icges(6,5) * t318;
t139 = Icges(6,1) * t203 + Icges(6,4) * t202 + Icges(6,5) * t319;
t138 = Icges(6,4) * t205 + Icges(6,2) * t204 + Icges(6,6) * t318;
t137 = Icges(6,4) * t203 + Icges(6,2) * t202 + Icges(6,6) * t319;
t136 = Icges(6,5) * t205 + Icges(6,6) * t204 + Icges(6,3) * t318;
t135 = Icges(6,5) * t203 + Icges(6,6) * t202 + Icges(6,3) * t319;
t134 = pkin(5) * t313 + t267 * t283;
t133 = -pkin(5) * t312 + t264 * t283;
t132 = rSges(7,1) * t189 + rSges(7,2) * t188 + rSges(7,3) * t318;
t131 = rSges(7,1) * t187 + rSges(7,2) * t186 + rSges(7,3) * t319;
t130 = Icges(7,1) * t189 + Icges(7,4) * t188 + Icges(7,5) * t318;
t129 = Icges(7,1) * t187 + Icges(7,4) * t186 + Icges(7,5) * t319;
t128 = Icges(7,4) * t189 + Icges(7,2) * t188 + Icges(7,6) * t318;
t127 = Icges(7,4) * t187 + Icges(7,2) * t186 + Icges(7,6) * t319;
t126 = Icges(7,5) * t189 + Icges(7,6) * t188 + Icges(7,3) * t318;
t125 = Icges(7,5) * t187 + Icges(7,6) * t186 + Icges(7,3) * t319;
t124 = t236 * t241 + (-t199 - t239) * t248 + t304;
t123 = t200 * t248 - t236 * t242 + t285;
t122 = t199 * t242 - t200 * t241 + t284;
t121 = t215 * t217 + (-t184 + t309) * t248 + t300;
t120 = t185 * t248 - t215 * t218 + t278;
t119 = t184 * t218 - t185 * t217 + t279;
t118 = t209 * t217 + (-t172 + t301) * t248 + t286;
t117 = t173 * t248 + (-t209 - t331) * t218 + t275;
t116 = t172 * t218 + (-t147 - t173) * t217 + t277;
t115 = -t141 * t216 + t159 * t176 + t276;
t114 = t142 * t216 - t159 * t177 + t273;
t113 = t141 * t177 - t142 * t176 + t274;
t112 = -t131 * t192 - t133 * t216 + t145 * t176 + t149 * t154 + t276;
t111 = t132 * t192 + t134 * t216 - t145 * t177 - t150 * t154 + t273;
t110 = t131 * t150 - t132 * t149 + t133 * t177 - t134 * t176 + t274;
t1 = t177 * ((t136 * t318 + t204 * t138 + t205 * t140) * t177 + (t135 * t318 + t137 * t204 + t139 * t205) * t176 + (t156 * t318 + t157 * t204 + t158 * t205) * t216) / 0.2e1 + t150 * ((t126 * t318 + t188 * t128 + t189 * t130) * t150 + (t125 * t318 + t127 * t188 + t129 * t189) * t149 + (t151 * t318 + t152 * t188 + t153 * t189) * t192) / 0.2e1 + t176 * ((t136 * t319 + t138 * t202 + t140 * t203) * t177 + (t135 * t319 + t202 * t137 + t203 * t139) * t176 + (t156 * t319 + t157 * t202 + t158 * t203) * t216) / 0.2e1 + t149 * ((t126 * t319 + t128 * t186 + t130 * t187) * t150 + (t125 * t319 + t186 * t127 + t187 * t129) * t149 + (t151 * t319 + t152 * t186 + t153 * t187) * t192) / 0.2e1 + t242 * (t280 * t264 + t270 * t267) / 0.2e1 + t241 * (t270 * t264 - t280 * t267) / 0.2e1 + t216 * ((-t135 * t176 - t136 * t177 - t156 * t216) * t244 + ((-t138 * t262 + t140 * t265) * t177 + (-t137 * t262 + t139 * t265) * t176 + (-t157 * t262 + t158 * t265) * t216) * t243) / 0.2e1 + t192 * ((-t125 * t149 - t126 * t150 - t151 * t192) * t244 + ((-t128 * t252 + t130 * t254) * t150 + (-t127 * t252 + t129 * t254) * t149 + (-t152 * t252 + t153 * t254) * t192) * t243) / 0.2e1 + m(1) * (t222 ^ 2 + t223 ^ 2 + t224 ^ 2) / 0.2e1 + m(2) * (t162 ^ 2 + t170 ^ 2 + t171 ^ 2) / 0.2e1 + m(4) * (t119 ^ 2 + t120 ^ 2 + t121 ^ 2) / 0.2e1 + m(3) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(7) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(6) * (t113 ^ 2 + t114 ^ 2 + t115 ^ 2) / 0.2e1 + m(5) * (t116 ^ 2 + t117 ^ 2 + t118 ^ 2) / 0.2e1 + (t337 * t264 - t336 * t267) * t217 / 0.2e1 + (t336 * t264 + t337 * t267) * t218 / 0.2e1 + ((-t231 * t264 + t234 * t267 + Icges(1,4)) * V_base(5) + (-t264 * t232 + t267 * t235 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t267 * t231 + t264 * t234 + Icges(1,2)) * V_base(5) + (t232 * t267 + t235 * t264 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t196 * t266 + t198 * t263) * t242 + (t195 * t266 + t197 * t263) * t241 + (t167 * t244 + t169 * t243 + t181 * t255 + t183 * t253) * t218 + (t166 * t244 + t168 * t243 + t180 * t255 + t182 * t253) * t217 + (t244 * t207 + t243 * t208 + t255 * t213 + t253 * t214 + t266 * t230 + t263 * t233 + Icges(2,3)) * t248) * t248 / 0.2e1 + t248 * V_base(4) * (Icges(2,5) * t267 - Icges(2,6) * t264) + t248 * V_base(5) * (Icges(2,5) * t264 + Icges(2,6) * t267) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
