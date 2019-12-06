% Calculate kinetic energy for
% S5PPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d4,d5,theta1,theta2,theta3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 14:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPPRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPPRR2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPPRR2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPPRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPPRR2_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPPRR2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PPPRR2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PPPRR2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 14:59:22
% EndTime: 2019-12-05 14:59:24
% DurationCPUTime: 2.06s
% Computational Cost: add. (1308->336), mult. (3043->467), div. (0->0), fcn. (3493->10), ass. (0->148)
t235 = sin(pkin(7));
t238 = cos(pkin(7));
t279 = Icges(2,5) * t238 - Icges(2,6) * t235 + Icges(1,5);
t278 = Icges(2,5) * t235 + Icges(2,6) * t238 + Icges(1,6);
t277 = cos(qJ(4));
t276 = Icges(2,4) * t235;
t234 = sin(pkin(8));
t275 = Icges(3,4) * t234;
t237 = cos(pkin(8));
t274 = Icges(3,4) * t237;
t233 = sin(pkin(9));
t273 = t233 * t234;
t272 = t234 * t235;
t271 = t234 * t238;
t240 = sin(qJ(4));
t270 = t234 * t240;
t269 = t235 * t237;
t268 = t237 * t238;
t218 = pkin(2) * t234 - qJ(3) * t237;
t267 = -qJ(1) - t218;
t254 = pkin(2) * t237 + qJ(3) * t234;
t199 = t254 * t235;
t220 = pkin(1) * t235 - qJ(2) * t238;
t266 = -t199 - t220;
t200 = t254 * t238;
t222 = pkin(1) * t238 + qJ(2) * t235;
t265 = -t200 - t222;
t264 = qJD(3) * t234;
t263 = V_base(5) * qJ(1) + V_base(1);
t259 = qJD(1) + V_base(3);
t258 = t234 * t277;
t236 = cos(pkin(9));
t195 = t233 * t268 - t235 * t236;
t175 = qJD(4) * t195 + V_base(4);
t193 = t233 * t269 + t236 * t238;
t174 = qJD(4) * t193 + V_base(5);
t208 = qJD(4) * t273 + V_base(6);
t257 = qJD(2) * t235 + t263;
t256 = V_base(4) * t220 + t259;
t255 = rSges(3,1) * t237 - rSges(3,2) * t234;
t253 = Icges(3,1) * t237 - t275;
t252 = -Icges(3,2) * t234 + t274;
t251 = Icges(3,5) * t237 - Icges(3,6) * t234;
t250 = -qJD(2) * t238 + V_base(6) * t222 + V_base(2);
t249 = V_base(5) * t218 + t238 * t264 + t257;
t248 = V_base(6) * t200 + t235 * t264 + t250;
t247 = -qJD(3) * t237 + V_base(4) * t199 + t256;
t246 = (-Icges(3,3) * t238 + t235 * t251) * V_base(5) + (Icges(3,3) * t235 + t238 * t251) * V_base(4) + (Icges(3,5) * t234 + Icges(3,6) * t237) * V_base(6);
t194 = -t233 * t238 + t236 * t269;
t161 = pkin(3) * t194 + pkin(5) * t193;
t201 = (pkin(3) * t236 + pkin(5) * t233) * t234;
t245 = V_base(5) * t201 + (-t161 + t266) * V_base(6) + t249;
t196 = t233 * t235 + t236 * t268;
t162 = pkin(3) * t196 + pkin(5) * t195;
t244 = V_base(6) * t162 + (-t201 + t267) * V_base(4) + t248;
t243 = V_base(4) * t161 + (-t162 + t265) * V_base(5) + t247;
t180 = -Icges(3,6) * t238 + t235 * t252;
t181 = Icges(3,6) * t235 + t238 * t252;
t183 = -Icges(3,5) * t238 + t235 * t253;
t184 = Icges(3,5) * t235 + t238 * t253;
t212 = Icges(3,2) * t237 + t275;
t215 = Icges(3,1) * t234 + t274;
t242 = (-t181 * t234 + t184 * t237) * V_base(4) + (-t180 * t234 + t183 * t237) * V_base(5) + (-t212 * t234 + t215 * t237) * V_base(6);
t241 = cos(qJ(5));
t239 = sin(qJ(5));
t231 = Icges(2,4) * t238;
t223 = rSges(2,1) * t238 - rSges(2,2) * t235;
t221 = rSges(2,1) * t235 + rSges(2,2) * t238;
t219 = rSges(3,1) * t234 + rSges(3,2) * t237;
t217 = Icges(2,1) * t238 - t276;
t216 = Icges(2,1) * t235 + t231;
t214 = -Icges(2,2) * t235 + t231;
t213 = Icges(2,2) * t238 + t276;
t207 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t206 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t205 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t198 = t236 * t258 - t237 * t240;
t197 = t236 * t270 + t237 * t277;
t187 = rSges(3,3) * t235 + t238 * t255;
t186 = -rSges(3,3) * t238 + t235 * t255;
t185 = -rSges(4,3) * t237 + (rSges(4,1) * t236 - rSges(4,2) * t233) * t234;
t182 = -Icges(4,5) * t237 + (Icges(4,1) * t236 - Icges(4,4) * t233) * t234;
t179 = -Icges(4,6) * t237 + (Icges(4,4) * t236 - Icges(4,2) * t233) * t234;
t176 = -Icges(4,3) * t237 + (Icges(4,5) * t236 - Icges(4,6) * t233) * t234;
t173 = V_base(5) * rSges(2,3) - t221 * V_base(6) + t263;
t172 = t223 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t171 = t198 * t241 + t239 * t273;
t170 = -t198 * t239 + t241 * t273;
t169 = qJD(5) * t197 + t208;
t168 = t196 * t277 + t238 * t270;
t167 = t196 * t240 - t238 * t258;
t166 = t194 * t277 + t235 * t270;
t165 = t194 * t240 - t235 * t258;
t164 = t221 * V_base(4) - t223 * V_base(5) + t259;
t163 = pkin(4) * t198 + pkin(6) * t197;
t158 = rSges(5,1) * t198 - rSges(5,2) * t197 + rSges(5,3) * t273;
t157 = Icges(5,1) * t198 - Icges(5,4) * t197 + Icges(5,5) * t273;
t156 = Icges(5,4) * t198 - Icges(5,2) * t197 + Icges(5,6) * t273;
t155 = Icges(5,5) * t198 - Icges(5,6) * t197 + Icges(5,3) * t273;
t154 = rSges(4,1) * t196 - rSges(4,2) * t195 + rSges(4,3) * t271;
t153 = rSges(4,1) * t194 - rSges(4,2) * t193 + rSges(4,3) * t272;
t152 = Icges(4,1) * t196 - Icges(4,4) * t195 + Icges(4,5) * t271;
t151 = Icges(4,1) * t194 - Icges(4,4) * t193 + Icges(4,5) * t272;
t150 = Icges(4,4) * t196 - Icges(4,2) * t195 + Icges(4,6) * t271;
t149 = Icges(4,4) * t194 - Icges(4,2) * t193 + Icges(4,6) * t272;
t148 = Icges(4,5) * t196 - Icges(4,6) * t195 + Icges(4,3) * t271;
t147 = Icges(4,5) * t194 - Icges(4,6) * t193 + Icges(4,3) * t272;
t146 = t168 * t241 + t195 * t239;
t145 = -t168 * t239 + t195 * t241;
t144 = t166 * t241 + t193 * t239;
t143 = -t166 * t239 + t193 * t241;
t142 = qJD(5) * t167 + t175;
t141 = qJD(5) * t165 + t174;
t140 = t219 * V_base(5) + (-t186 - t220) * V_base(6) + t257;
t139 = t187 * V_base(6) + (-qJ(1) - t219) * V_base(4) + t250;
t138 = pkin(4) * t168 + pkin(6) * t167;
t137 = pkin(4) * t166 + pkin(6) * t165;
t136 = t186 * V_base(4) + (-t187 - t222) * V_base(5) + t256;
t135 = rSges(6,1) * t171 + rSges(6,2) * t170 + rSges(6,3) * t197;
t134 = Icges(6,1) * t171 + Icges(6,4) * t170 + Icges(6,5) * t197;
t133 = Icges(6,4) * t171 + Icges(6,2) * t170 + Icges(6,6) * t197;
t132 = Icges(6,5) * t171 + Icges(6,6) * t170 + Icges(6,3) * t197;
t131 = rSges(5,1) * t168 - rSges(5,2) * t167 + rSges(5,3) * t195;
t130 = rSges(5,1) * t166 - rSges(5,2) * t165 + rSges(5,3) * t193;
t129 = Icges(5,1) * t168 - Icges(5,4) * t167 + Icges(5,5) * t195;
t128 = Icges(5,1) * t166 - Icges(5,4) * t165 + Icges(5,5) * t193;
t127 = Icges(5,4) * t168 - Icges(5,2) * t167 + Icges(5,6) * t195;
t126 = Icges(5,4) * t166 - Icges(5,2) * t165 + Icges(5,6) * t193;
t125 = Icges(5,5) * t168 - Icges(5,6) * t167 + Icges(5,3) * t195;
t124 = Icges(5,5) * t166 - Icges(5,6) * t165 + Icges(5,3) * t193;
t123 = t185 * V_base(5) + (-t153 + t266) * V_base(6) + t249;
t122 = t154 * V_base(6) + (-t185 + t267) * V_base(4) + t248;
t121 = rSges(6,1) * t146 + rSges(6,2) * t145 + rSges(6,3) * t167;
t120 = rSges(6,1) * t144 + rSges(6,2) * t143 + rSges(6,3) * t165;
t119 = Icges(6,1) * t146 + Icges(6,4) * t145 + Icges(6,5) * t167;
t118 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t165;
t117 = Icges(6,4) * t146 + Icges(6,2) * t145 + Icges(6,6) * t167;
t116 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t165;
t115 = Icges(6,5) * t146 + Icges(6,6) * t145 + Icges(6,3) * t167;
t114 = Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t165;
t113 = t153 * V_base(4) + (-t154 + t265) * V_base(5) + t247;
t112 = -t130 * t208 + t158 * t174 + t245;
t111 = t131 * t208 - t158 * t175 + t244;
t110 = t130 * t175 - t131 * t174 + t243;
t109 = -t120 * t169 + t135 * t141 - t137 * t208 + t163 * t174 + t245;
t108 = t121 * t169 - t135 * t142 + t138 * t208 - t163 * t175 + t244;
t107 = t120 * t142 - t121 * t141 + t137 * t175 - t138 * t174 + t243;
t1 = m(1) * (t205 ^ 2 + t206 ^ 2 + t207 ^ 2) / 0.2e1 + m(2) * (t164 ^ 2 + t172 ^ 2 + t173 ^ 2) / 0.2e1 + m(3) * (t136 ^ 2 + t139 ^ 2 + t140 ^ 2) / 0.2e1 + m(4) * (t113 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(5) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + t175 * ((t125 * t195 - t127 * t167 + t129 * t168) * t175 + (t124 * t195 - t126 * t167 + t128 * t168) * t174 + (t155 * t195 - t156 * t167 + t157 * t168) * t208) / 0.2e1 + t174 * ((t125 * t193 - t127 * t165 + t129 * t166) * t175 + (t124 * t193 - t126 * t165 + t128 * t166) * t174 + (t155 * t193 - t156 * t165 + t157 * t166) * t208) / 0.2e1 + t208 * ((t125 * t273 - t127 * t197 + t129 * t198) * t175 + (t124 * t273 - t126 * t197 + t128 * t198) * t174 + (t155 * t273 - t156 * t197 + t157 * t198) * t208) / 0.2e1 + m(6) * (t107 ^ 2 + t108 ^ 2 + t109 ^ 2) / 0.2e1 + t142 * ((t167 * t115 + t145 * t117 + t146 * t119) * t142 + (t114 * t167 + t116 * t145 + t118 * t146) * t141 + (t132 * t167 + t133 * t145 + t134 * t146) * t169) / 0.2e1 + t141 * ((t115 * t165 + t117 * t143 + t119 * t144) * t142 + (t165 * t114 + t143 * t116 + t144 * t118) * t141 + (t132 * t165 + t133 * t143 + t134 * t144) * t169) / 0.2e1 + t169 * ((t115 * t197 + t117 * t170 + t119 * t171) * t142 + (t114 * t197 + t116 * t170 + t118 * t171) * t141 + (t132 * t197 + t133 * t170 + t134 * t171) * t169) / 0.2e1 + (t235 * t246 + t238 * t242 + (t176 * t271 - t179 * t195 + t182 * t196 + t279) * V_base(6) + (t147 * t271 - t149 * t195 + t151 * t196 - t213 * t235 + t216 * t238 + Icges(1,4)) * V_base(5) + (t148 * t271 - t150 * t195 + t152 * t196 - t214 * t235 + t217 * t238 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t235 * t242 - t238 * t246 + (t176 * t272 - t179 * t193 + t182 * t194 + t278) * V_base(6) + (t147 * t272 - t149 * t193 + t151 * t194 + t213 * t238 + t216 * t235 + Icges(1,2)) * V_base(5) + (t148 * t272 - t150 * t193 + t152 * t194 + t214 * t238 + t217 * t235 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((Icges(1,3) + Icges(2,3) + (t212 - t176) * t237 + (-t179 * t233 + t182 * t236 + t215) * t234) * V_base(6) + ((t180 - t147) * t237 + (-t149 * t233 + t151 * t236 + t183) * t234 + t278) * V_base(5) + ((t181 - t148) * t237 + (-t150 * t233 + t152 * t236 + t184) * t234 + t279) * V_base(4)) * V_base(6) / 0.2e1;
T = t1;
