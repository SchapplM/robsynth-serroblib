% Calculate kinetic energy for
% S6RPRPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR2_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR2_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR2_energykin_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR2_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR2_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPRPRR2_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:36:53
% EndTime: 2019-03-09 03:36:56
% DurationCPUTime: 2.91s
% Computational Cost: add. (2288->343), mult. (1758->485), div. (0->0), fcn. (1594->12), ass. (0->170)
t319 = Icges(4,3) + Icges(5,3);
t244 = qJ(3) + pkin(11);
t234 = sin(t244);
t236 = cos(t244);
t249 = sin(qJ(3));
t252 = cos(qJ(3));
t318 = Icges(4,5) * t252 + Icges(5,5) * t236 - Icges(4,6) * t249 - Icges(5,6) * t234;
t245 = qJ(1) + pkin(10);
t235 = sin(t245);
t237 = cos(t245);
t297 = Icges(5,4) * t236;
t270 = -Icges(5,2) * t234 + t297;
t149 = -Icges(5,6) * t237 + t235 * t270;
t150 = Icges(5,6) * t235 + t237 * t270;
t298 = Icges(5,4) * t234;
t272 = Icges(5,1) * t236 - t298;
t151 = -Icges(5,5) * t237 + t235 * t272;
t152 = Icges(5,5) * t235 + t237 * t272;
t299 = Icges(4,4) * t252;
t271 = -Icges(4,2) * t249 + t299;
t165 = -Icges(4,6) * t237 + t235 * t271;
t166 = Icges(4,6) * t235 + t237 * t271;
t300 = Icges(4,4) * t249;
t273 = Icges(4,1) * t252 - t300;
t168 = -Icges(4,5) * t237 + t235 * t273;
t169 = Icges(4,5) * t235 + t237 * t273;
t195 = Icges(5,2) * t236 + t298;
t198 = Icges(5,1) * t234 + t297;
t213 = -qJD(3) * t237 + V_base(5);
t214 = qJD(3) * t235 + V_base(4);
t218 = Icges(4,2) * t252 + t300;
t221 = Icges(4,1) * t249 + t299;
t238 = V_base(6) + qJD(1);
t315 = (-t195 * t234 + t198 * t236 - t218 * t249 + t221 * t252) * t238 + (-t150 * t234 + t152 * t236 - t166 * t249 + t169 * t252) * t214 + (-t149 * t234 + t151 * t236 - t165 * t249 + t168 * t252) * t213;
t314 = (Icges(4,5) * t249 + Icges(5,5) * t234 + Icges(4,6) * t252 + Icges(5,6) * t236) * t238 + (t319 * t235 + t318 * t237) * t214 + (t318 * t235 - t319 * t237) * t213;
t250 = sin(qJ(1));
t310 = pkin(1) * t250;
t253 = cos(qJ(1));
t309 = pkin(1) * t253;
t308 = pkin(3) * t249;
t307 = pkin(3) * t252;
t251 = cos(qJ(5));
t306 = pkin(5) * t251;
t305 = -pkin(6) - qJ(2);
t302 = Icges(2,4) * t250;
t301 = Icges(3,4) * t235;
t296 = t234 * t235;
t295 = t234 * t237;
t246 = qJ(5) + qJ(6);
t239 = sin(t246);
t294 = t235 * t239;
t240 = cos(t246);
t293 = t235 * t240;
t248 = sin(qJ(5));
t292 = t235 * t248;
t291 = t235 * t251;
t290 = t237 * t239;
t289 = t237 * t240;
t288 = t237 * t248;
t287 = t237 * t251;
t286 = qJD(5) * t234;
t285 = qJD(6) * t234;
t284 = t238 * t309 + V_base(2);
t283 = V_base(5) * pkin(6) + V_base(1);
t205 = pkin(2) * t235 - pkin(7) * t237;
t280 = -t205 - t310;
t279 = V_base(5) * qJ(2) + t283;
t278 = V_base(4) * t310 + qJD(2) + V_base(3);
t183 = t237 * t286 + t214;
t144 = -qJ(4) * t237 + t235 * t307;
t277 = -t144 + t280;
t276 = pkin(4) * t236 + pkin(8) * t234;
t275 = rSges(4,1) * t252 - rSges(4,2) * t249;
t274 = rSges(5,1) * t236 - rSges(5,2) * t234;
t267 = qJD(4) * t235 + t213 * t308 + t279;
t182 = t235 * t286 + t213;
t266 = pkin(9) * t234 + t236 * t306;
t206 = pkin(2) * t237 + pkin(7) * t235;
t263 = t238 * t206 + t305 * V_base(4) + t284;
t262 = V_base(4) * t205 + (-t206 - t309) * V_base(5) + t278;
t261 = t214 * t144 + t262;
t145 = qJ(4) * t235 + t237 * t307;
t260 = -qJD(4) * t237 + t238 * t145 + t263;
t180 = t276 * t235;
t204 = t234 * pkin(4) - t236 * pkin(8);
t259 = t213 * t204 + (-t180 + t277) * t238 + t267;
t181 = t276 * t237;
t258 = t214 * t180 + (-t145 - t181) * t213 + t261;
t257 = t238 * t181 + (-t204 - t308) * t214 + t260;
t242 = Icges(2,4) * t253;
t231 = Icges(3,4) * t237;
t226 = rSges(2,1) * t253 - rSges(2,2) * t250;
t225 = rSges(2,1) * t250 + rSges(2,2) * t253;
t224 = rSges(4,1) * t249 + rSges(4,2) * t252;
t223 = Icges(2,1) * t253 - t302;
t222 = Icges(2,1) * t250 + t242;
t220 = -Icges(2,2) * t250 + t242;
t219 = Icges(2,2) * t253 + t302;
t212 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t211 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t210 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t209 = -qJD(5) * t236 + t238;
t203 = rSges(3,1) * t237 - rSges(3,2) * t235;
t202 = rSges(3,1) * t235 + rSges(3,2) * t237;
t201 = rSges(5,1) * t234 + rSges(5,2) * t236;
t200 = Icges(3,1) * t237 - t301;
t199 = Icges(3,1) * t235 + t231;
t197 = -Icges(3,2) * t235 + t231;
t196 = Icges(3,2) * t237 + t301;
t189 = (-qJD(5) - qJD(6)) * t236 + t238;
t187 = t236 * t287 + t292;
t186 = -t236 * t288 + t291;
t185 = t236 * t291 - t288;
t184 = -t236 * t292 - t287;
t178 = t236 * t289 + t294;
t177 = -t236 * t290 + t293;
t176 = t236 * t293 - t290;
t175 = -t236 * t294 - t289;
t174 = rSges(4,3) * t235 + t237 * t275;
t173 = -rSges(4,3) * t237 + t235 * t275;
t172 = -rSges(6,3) * t236 + (rSges(6,1) * t251 - rSges(6,2) * t248) * t234;
t171 = V_base(5) * rSges(2,3) - t225 * t238 + t283;
t170 = t226 * t238 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t167 = -Icges(6,5) * t236 + (Icges(6,1) * t251 - Icges(6,4) * t248) * t234;
t164 = -Icges(6,6) * t236 + (Icges(6,4) * t251 - Icges(6,2) * t248) * t234;
t161 = -Icges(6,3) * t236 + (Icges(6,5) * t251 - Icges(6,6) * t248) * t234;
t160 = t225 * V_base(4) - t226 * V_base(5) + V_base(3);
t158 = -rSges(7,3) * t236 + (rSges(7,1) * t240 - rSges(7,2) * t239) * t234;
t157 = -Icges(7,5) * t236 + (Icges(7,1) * t240 - Icges(7,4) * t239) * t234;
t156 = -Icges(7,6) * t236 + (Icges(7,4) * t240 - Icges(7,2) * t239) * t234;
t155 = -Icges(7,3) * t236 + (Icges(7,5) * t240 - Icges(7,6) * t239) * t234;
t154 = rSges(5,3) * t235 + t237 * t274;
t153 = -rSges(5,3) * t237 + t235 * t274;
t146 = -pkin(9) * t236 + t234 * t306;
t143 = t237 * t285 + t183;
t142 = t235 * t285 + t182;
t139 = V_base(5) * rSges(3,3) + (-t202 - t310) * t238 + t279;
t138 = t203 * t238 + (-rSges(3,3) + t305) * V_base(4) + t284;
t136 = t202 * V_base(4) + (-t203 - t309) * V_base(5) + t278;
t135 = rSges(6,1) * t187 + rSges(6,2) * t186 + rSges(6,3) * t295;
t134 = rSges(6,1) * t185 + rSges(6,2) * t184 + rSges(6,3) * t296;
t133 = Icges(6,1) * t187 + Icges(6,4) * t186 + Icges(6,5) * t295;
t132 = Icges(6,1) * t185 + Icges(6,4) * t184 + Icges(6,5) * t296;
t131 = Icges(6,4) * t187 + Icges(6,2) * t186 + Icges(6,6) * t295;
t130 = Icges(6,4) * t185 + Icges(6,2) * t184 + Icges(6,6) * t296;
t129 = Icges(6,5) * t187 + Icges(6,6) * t186 + Icges(6,3) * t295;
t128 = Icges(6,5) * t185 + Icges(6,6) * t184 + Icges(6,3) * t296;
t127 = pkin(5) * t292 + t237 * t266;
t126 = -pkin(5) * t288 + t235 * t266;
t125 = rSges(7,1) * t178 + rSges(7,2) * t177 + rSges(7,3) * t295;
t124 = rSges(7,1) * t176 + rSges(7,2) * t175 + rSges(7,3) * t296;
t123 = Icges(7,1) * t178 + Icges(7,4) * t177 + Icges(7,5) * t295;
t122 = Icges(7,1) * t176 + Icges(7,4) * t175 + Icges(7,5) * t296;
t121 = Icges(7,4) * t178 + Icges(7,2) * t177 + Icges(7,6) * t295;
t120 = Icges(7,4) * t176 + Icges(7,2) * t175 + Icges(7,6) * t296;
t119 = Icges(7,5) * t178 + Icges(7,6) * t177 + Icges(7,3) * t295;
t118 = Icges(7,5) * t176 + Icges(7,6) * t175 + Icges(7,3) * t296;
t117 = t213 * t224 + (-t173 + t280) * t238 + t279;
t116 = t174 * t238 - t214 * t224 + t263;
t115 = t173 * t214 - t174 * t213 + t262;
t114 = t201 * t213 + (-t153 + t277) * t238 + t267;
t113 = t154 * t238 + (-t201 - t308) * t214 + t260;
t112 = t153 * t214 + (-t145 - t154) * t213 + t261;
t111 = -t134 * t209 + t172 * t182 + t259;
t110 = t135 * t209 - t172 * t183 + t257;
t109 = t134 * t183 - t135 * t182 + t258;
t108 = -t124 * t189 - t126 * t209 + t142 * t158 + t146 * t182 + t259;
t107 = t125 * t189 + t127 * t209 - t143 * t158 - t146 * t183 + t257;
t106 = t124 * t143 - t125 * t142 + t126 * t183 - t127 * t182 + t258;
t1 = t183 * ((t129 * t295 + t186 * t131 + t187 * t133) * t183 + (t128 * t295 + t130 * t186 + t132 * t187) * t182 + (t161 * t295 + t164 * t186 + t167 * t187) * t209) / 0.2e1 + t143 * ((t119 * t295 + t177 * t121 + t178 * t123) * t143 + (t118 * t295 + t120 * t177 + t122 * t178) * t142 + (t155 * t295 + t156 * t177 + t157 * t178) * t189) / 0.2e1 + t182 * ((t129 * t296 + t131 * t184 + t133 * t185) * t183 + (t128 * t296 + t184 * t130 + t185 * t132) * t182 + (t161 * t296 + t164 * t184 + t167 * t185) * t209) / 0.2e1 + t142 * ((t119 * t296 + t121 * t175 + t123 * t176) * t143 + (t118 * t296 + t175 * t120 + t176 * t122) * t142 + (t155 * t296 + t156 * t175 + t157 * t176) * t189) / 0.2e1 + t209 * ((-t128 * t182 - t129 * t183 - t161 * t209) * t236 + ((-t131 * t248 + t133 * t251) * t183 + (-t130 * t248 + t132 * t251) * t182 + (-t164 * t248 + t167 * t251) * t209) * t234) / 0.2e1 + t189 * ((-t118 * t142 - t119 * t143 - t155 * t189) * t236 + ((-t121 * t239 + t123 * t240) * t143 + (-t120 * t239 + t122 * t240) * t142 + (-t156 * t239 + t157 * t240) * t189) * t234) / 0.2e1 + m(1) * (t210 ^ 2 + t211 ^ 2 + t212 ^ 2) / 0.2e1 + m(2) * (t160 ^ 2 + t170 ^ 2 + t171 ^ 2) / 0.2e1 + m(3) * (t136 ^ 2 + t138 ^ 2 + t139 ^ 2) / 0.2e1 + m(5) * (t112 ^ 2 + t113 ^ 2 + t114 ^ 2) / 0.2e1 + m(4) * (t115 ^ 2 + t116 ^ 2 + t117 ^ 2) / 0.2e1 + m(7) * (t106 ^ 2 + t107 ^ 2 + t108 ^ 2) / 0.2e1 + m(6) * (t109 ^ 2 + t110 ^ 2 + t111 ^ 2) / 0.2e1 + (t315 * t235 - t314 * t237) * t213 / 0.2e1 + (t314 * t235 + t315 * t237) * t214 / 0.2e1 + ((-t196 * t235 + t199 * t237 - t219 * t250 + t222 * t253 + Icges(1,4)) * V_base(5) + (-t235 * t197 + t237 * t200 - t250 * t220 + t253 * t223 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t237 * t196 + t235 * t199 + t253 * t219 + t250 * t222 + Icges(1,2)) * V_base(5) + (t197 * t237 + t200 * t235 + t220 * t253 + t223 * t250 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t150 * t236 + t152 * t234 + t166 * t252 + t169 * t249) * t214 + (t149 * t236 + t151 * t234 + t165 * t252 + t168 * t249) * t213 + (t236 * t195 + t234 * t198 + t252 * t218 + t249 * t221 + Icges(2,3) + Icges(3,3)) * t238) * t238 / 0.2e1 + t238 * V_base(5) * (Icges(2,5) * t250 + Icges(3,5) * t235 + Icges(2,6) * t253 + Icges(3,6) * t237) + t238 * V_base(4) * (Icges(2,5) * t253 + Icges(3,5) * t237 - Icges(2,6) * t250 - Icges(3,6) * t235) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
