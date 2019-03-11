% Calculate kinetic energy for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_energykin_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_energykin_floatb_twist_slag_vp1: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPPPRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_energykin_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_energykin_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPPPRR3_energykin_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [7 6]), ...
  'S6RPPPRR3_energykin_floatb_twist_slag_vp1: Icges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:24
% EndTime: 2019-03-09 01:33:27
% DurationCPUTime: 2.54s
% Computational Cost: add. (1431->301), mult. (2160->406), div. (0->0), fcn. (2394->10), ass. (0->152)
t287 = Icges(2,4) - Icges(3,5);
t286 = Icges(2,1) + Icges(3,1);
t285 = Icges(3,4) + Icges(2,5);
t284 = Icges(2,2) + Icges(3,3);
t283 = Icges(2,6) - Icges(3,6);
t272 = sin(qJ(1));
t282 = t287 * t272;
t273 = cos(qJ(1));
t281 = t287 * t273;
t280 = -t284 * t273 - t282;
t279 = t284 * t272 - t281;
t278 = t286 * t272 + t281;
t277 = t286 * t273 - t282;
t268 = sin(pkin(9));
t269 = cos(pkin(9));
t171 = -t268 * t272 - t269 * t273;
t172 = t268 * t273 - t269 * t272;
t276 = Icges(4,5) * t172 - Icges(4,6) * t171 + t285 * t272 + t283 * t273;
t275 = Icges(4,5) * t171 + Icges(4,6) * t172 - t283 * t272 + t285 * t273;
t216 = sin(pkin(10));
t271 = pkin(4) * t216;
t217 = cos(pkin(10));
t270 = pkin(4) * t217;
t267 = Icges(4,4) * t171;
t266 = Icges(5,4) * t216;
t265 = Icges(5,4) * t217;
t215 = pkin(10) + qJ(5);
t207 = sin(t215);
t264 = Icges(6,4) * t207;
t208 = cos(t215);
t263 = Icges(6,4) * t208;
t262 = t171 * t207;
t261 = t172 * t207;
t219 = sin(qJ(6));
t260 = t208 * t219;
t220 = cos(qJ(6));
t259 = t208 * t220;
t257 = qJD(6) * t207;
t194 = pkin(1) * t272 - qJ(2) * t273;
t256 = V_base(4) * t194 + V_base(3);
t255 = V_base(5) * pkin(6) + V_base(1);
t252 = t273 * pkin(2);
t251 = t272 * pkin(2);
t159 = qJD(5) * t172 + V_base(4);
t158 = -qJD(5) * t171 + V_base(5);
t209 = V_base(6) + qJD(1);
t248 = qJD(2) * t272 + t255;
t247 = -t194 - t251;
t197 = pkin(1) * t273 + qJ(2) * t272;
t246 = -t197 - t252;
t245 = qJD(4) * t172 + t248;
t244 = -pkin(5) * t208 - pkin(8) * t207;
t243 = V_base(4) * t251 - qJD(3) + t256;
t242 = -rSges(5,1) * t217 + rSges(5,2) * t216;
t241 = -rSges(6,1) * t208 + rSges(6,2) * t207;
t240 = -Icges(5,1) * t217 + t266;
t239 = -Icges(6,1) * t208 + t264;
t238 = Icges(5,2) * t216 - t265;
t237 = Icges(6,2) * t207 - t263;
t236 = -Icges(5,5) * t217 + Icges(5,6) * t216;
t235 = -Icges(6,5) * t208 + Icges(6,6) * t207;
t147 = -pkin(3) * t172 - qJ(4) * t171;
t234 = -t147 + t247;
t149 = -pkin(3) * t171 + qJ(4) * t172;
t233 = -t149 + t246;
t232 = V_base(4) * t147 + t243;
t231 = -qJD(2) * t273 + t209 * t197 + V_base(2);
t110 = -pkin(7) * t171 - t172 * t270;
t230 = -t110 + t234;
t229 = (-Icges(6,3) * t171 + t172 * t235) * t158 + (Icges(6,3) * t172 + t171 * t235) * t159 + (-Icges(6,5) * t207 - Icges(6,6) * t208) * t209;
t228 = V_base(4) * qJ(3) + t209 * t252 + t231;
t227 = (-Icges(5,3) * t171 + t172 * t236) * V_base(5) + (Icges(5,3) * t172 + t171 * t236) * V_base(4) + (-Icges(5,5) * t216 - Icges(5,6) * t217) * t209;
t226 = -qJD(4) * t171 + t209 * t149 + t228;
t225 = (-qJ(3) - t271) * V_base(5) + t245;
t111 = pkin(7) * t172 - t171 * t270;
t224 = t209 * t111 + t226 + (-pkin(6) + t271) * V_base(4);
t223 = V_base(4) * t110 + (-t111 + t233) * V_base(5) + t232;
t114 = -Icges(6,6) * t171 + t172 * t237;
t115 = Icges(6,6) * t172 + t171 * t237;
t116 = -Icges(6,5) * t171 + t172 * t239;
t117 = Icges(6,5) * t172 + t171 * t239;
t167 = -Icges(6,2) * t208 - t264;
t168 = -Icges(6,1) * t207 - t263;
t222 = (t115 * t207 - t117 * t208) * t159 + (t114 * t207 - t116 * t208) * t158 + (t167 * t207 - t168 * t208) * t209;
t123 = -Icges(5,6) * t171 + t172 * t238;
t124 = Icges(5,6) * t172 + t171 * t238;
t125 = -Icges(5,5) * t171 + t172 * t240;
t126 = Icges(5,5) * t172 + t171 * t240;
t179 = -Icges(5,2) * t217 - t266;
t180 = -Icges(5,1) * t216 - t265;
t221 = (t124 * t216 - t126 * t217) * V_base(4) + (t123 * t216 - t125 * t217) * V_base(5) + (t179 * t216 - t180 * t217) * t209;
t199 = rSges(2,1) * t273 - rSges(2,2) * t272;
t198 = rSges(3,1) * t273 + rSges(3,3) * t272;
t196 = rSges(2,1) * t272 + rSges(2,2) * t273;
t195 = rSges(3,1) * t272 - rSges(3,3) * t273;
t181 = -rSges(5,1) * t216 - rSges(5,2) * t217;
t177 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t176 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t175 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t174 = qJD(6) * t208 + t209;
t170 = -pkin(5) * t207 + pkin(8) * t208;
t169 = -rSges(6,1) * t207 - rSges(6,2) * t208;
t164 = Icges(4,4) * t172;
t157 = t208 * rSges(7,3) + (-rSges(7,1) * t220 + rSges(7,2) * t219) * t207;
t156 = V_base(5) * rSges(2,3) - t196 * t209 + t255;
t155 = t199 * t209 + V_base(2) + (-rSges(2,3) - pkin(6)) * V_base(4);
t154 = Icges(7,5) * t208 + (-Icges(7,1) * t220 + Icges(7,4) * t219) * t207;
t153 = Icges(7,6) * t208 + (-Icges(7,4) * t220 + Icges(7,2) * t219) * t207;
t152 = Icges(7,3) * t208 + (-Icges(7,5) * t220 + Icges(7,6) * t219) * t207;
t151 = t196 * V_base(4) - t199 * V_base(5) + V_base(3);
t150 = -rSges(4,1) * t171 - rSges(4,2) * t172;
t148 = -rSges(4,1) * t172 + rSges(4,2) * t171;
t146 = -Icges(4,1) * t171 - t164;
t145 = -Icges(4,1) * t172 + t267;
t144 = -Icges(4,2) * t172 - t267;
t143 = Icges(4,2) * t171 - t164;
t138 = -t171 * t259 + t172 * t219;
t137 = t171 * t260 + t172 * t220;
t136 = -t171 * t219 - t172 * t259;
t135 = -t171 * t220 + t172 * t260;
t134 = -t171 * t257 + t159;
t133 = -t172 * t257 + t158;
t132 = t244 * t171;
t131 = t244 * t172;
t130 = V_base(5) * rSges(3,2) + (-t194 - t195) * t209 + t248;
t129 = t209 * t198 + (-rSges(3,2) - pkin(6)) * V_base(4) + t231;
t128 = rSges(5,3) * t172 + t171 * t242;
t127 = -rSges(5,3) * t171 + t172 * t242;
t120 = t195 * V_base(4) + (-t197 - t198) * V_base(5) + t256;
t119 = rSges(6,3) * t172 + t171 * t241;
t118 = -rSges(6,3) * t171 + t172 * t241;
t107 = (-qJ(3) - rSges(4,3)) * V_base(5) + (-t148 + t247) * t209 + t248;
t106 = t209 * t150 + (rSges(4,3) - pkin(6)) * V_base(4) + t228;
t105 = rSges(7,1) * t138 + rSges(7,2) * t137 - rSges(7,3) * t262;
t104 = rSges(7,1) * t136 + rSges(7,2) * t135 - rSges(7,3) * t261;
t103 = Icges(7,1) * t138 + Icges(7,4) * t137 - Icges(7,5) * t262;
t102 = Icges(7,1) * t136 + Icges(7,4) * t135 - Icges(7,5) * t261;
t101 = Icges(7,4) * t138 + Icges(7,2) * t137 - Icges(7,6) * t262;
t100 = Icges(7,4) * t136 + Icges(7,2) * t135 - Icges(7,6) * t261;
t99 = Icges(7,5) * t138 + Icges(7,6) * t137 - Icges(7,3) * t262;
t98 = Icges(7,5) * t136 + Icges(7,6) * t135 - Icges(7,3) * t261;
t97 = V_base(4) * t148 + (-t150 + t246) * V_base(5) + t243;
t96 = (-qJ(3) + t181) * V_base(5) + (-t127 + t234) * t209 + t245;
t95 = t209 * t128 + (-pkin(6) - t181) * V_base(4) + t226;
t94 = V_base(4) * t127 + (-t128 + t233) * V_base(5) + t232;
t93 = t158 * t169 + (-t118 + t230) * t209 + t225;
t92 = t209 * t119 - t159 * t169 + t224;
t91 = t159 * t118 - t158 * t119 + t223;
t90 = -t174 * t104 + t133 * t157 + t158 * t170 + (-t131 + t230) * t209 + t225;
t89 = t174 * t105 + t209 * t132 - t134 * t157 - t159 * t170 + t224;
t88 = t134 * t104 - t133 * t105 + t159 * t131 - t158 * t132 + t223;
t1 = m(7) * (t88 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(6) * (t91 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(5) * (t94 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + t159 * (t222 * t171 + t229 * t172) / 0.2e1 + t158 * (-t229 * t171 + t222 * t172) / 0.2e1 + m(3) * (t120 ^ 2 + t129 ^ 2 + t130 ^ 2) / 0.2e1 + m(1) * (t175 ^ 2 + t176 ^ 2 + t177 ^ 2) / 0.2e1 + m(4) * (t106 ^ 2 + t107 ^ 2 + t97 ^ 2) / 0.2e1 + m(2) * (t151 ^ 2 + t155 ^ 2 + t156 ^ 2) / 0.2e1 + t133 * ((t101 * t135 + t103 * t136 - t261 * t99) * t134 + (t135 * t100 + t136 * t102 - t98 * t261) * t133 + (t135 * t153 + t136 * t154 - t152 * t261) * t174) / 0.2e1 + t134 * ((t137 * t101 + t138 * t103 - t99 * t262) * t134 + (t100 * t137 + t102 * t138 - t262 * t98) * t133 + (t137 * t153 + t138 * t154 - t152 * t262) * t174) / 0.2e1 + t174 * ((t133 * t98 + t134 * t99 + t152 * t174) * t208 + ((t101 * t219 - t103 * t220) * t134 + (t100 * t219 - t102 * t220) * t133 + (t153 * t219 - t154 * t220) * t174) * t207) / 0.2e1 + ((-t115 * t208 - t117 * t207) * t159 + (-t114 * t208 - t116 * t207) * t158 + (-t123 * t217 - t125 * t216 + t276) * V_base(5) + (-t124 * t217 - t126 * t216 + t275) * V_base(4) + (-t167 * t208 - t168 * t207 - t179 * t217 - t180 * t216 + Icges(3,2) + Icges(2,3) + Icges(4,3)) * t209) * t209 / 0.2e1 + (t171 * t221 + t172 * t227 + t275 * t209 + (-t143 * t172 - t145 * t171 + t272 * t280 + t278 * t273 + Icges(1,4)) * V_base(5) + (-t144 * t172 - t146 * t171 + t272 * t279 + t273 * t277 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (-t227 * t171 + t221 * t172 + t276 * t209 + (t143 * t171 - t145 * t172 + t278 * t272 - t273 * t280 + Icges(1,2)) * V_base(5) + (t144 * t171 - t146 * t172 + t272 * t277 - t273 * t279 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T  = t1;
