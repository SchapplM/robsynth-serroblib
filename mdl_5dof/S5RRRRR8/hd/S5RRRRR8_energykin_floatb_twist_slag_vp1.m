% Calculate kinetic energy for
% S5RRRRR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR8_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR8_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRRR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR8_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR8_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR8_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR8_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:24:19
% EndTime: 2019-12-31 22:24:21
% DurationCPUTime: 2.38s
% Computational Cost: add. (1529->307), mult. (1682->471), div. (0->0), fcn. (1576->10), ass. (0->151)
t206 = sin(qJ(2));
t258 = pkin(2) * t206;
t209 = cos(qJ(2));
t257 = pkin(2) * t209;
t208 = cos(qJ(4));
t256 = pkin(4) * t208;
t207 = sin(qJ(1));
t253 = Icges(2,4) * t207;
t252 = Icges(3,4) * t206;
t251 = Icges(3,4) * t209;
t204 = qJ(2) + qJ(3);
t198 = sin(t204);
t250 = Icges(4,4) * t198;
t200 = cos(t204);
t249 = Icges(4,4) * t200;
t248 = t198 * t207;
t210 = cos(qJ(1));
t247 = t198 * t210;
t246 = t200 * t207;
t245 = t200 * t210;
t205 = sin(qJ(4));
t244 = t205 * t207;
t243 = t205 * t210;
t242 = t207 * t208;
t241 = t208 * t210;
t128 = -pkin(7) * t210 + t207 * t257;
t188 = t207 * pkin(1) - t210 * pkin(6);
t240 = -t128 - t188;
t239 = qJD(4) * t198;
t238 = qJD(5) * t198;
t237 = V_base(5) * pkin(5) + V_base(1);
t191 = qJD(2) * t207 + V_base(4);
t194 = V_base(6) + qJD(1);
t190 = -qJD(2) * t210 + V_base(5);
t234 = t190 * t258 + t237;
t167 = qJD(3) * t207 + t191;
t233 = pkin(3) * t200 + pkin(8) * t198;
t232 = rSges(3,1) * t209 - rSges(3,2) * t206;
t231 = rSges(4,1) * t200 - rSges(4,2) * t198;
t140 = t210 * t239 + t167;
t230 = Icges(3,1) * t209 - t252;
t229 = Icges(4,1) * t200 - t250;
t228 = -Icges(3,2) * t206 + t251;
t227 = -Icges(4,2) * t198 + t249;
t226 = Icges(3,5) * t209 - Icges(3,6) * t206;
t225 = Icges(4,5) * t200 - Icges(4,6) * t198;
t189 = t210 * pkin(1) + t207 * pkin(6);
t224 = -V_base(4) * pkin(5) + t194 * t189 + V_base(2);
t223 = V_base(4) * t188 - t189 * V_base(5) + V_base(3);
t166 = V_base(5) + (-qJD(2) - qJD(3)) * t210;
t139 = t207 * t239 + t166;
t222 = pkin(9) * t198 + t200 * t256;
t221 = (-Icges(4,3) * t210 + t207 * t225) * t166 + (Icges(4,3) * t207 + t210 * t225) * t167 + (Icges(4,5) * t198 + Icges(4,6) * t200) * t194;
t220 = (-Icges(3,3) * t210 + t207 * t226) * t190 + (Icges(3,3) * t207 + t210 * t226) * t191 + (Icges(3,5) * t206 + Icges(3,6) * t209) * t194;
t153 = t233 * t207;
t165 = pkin(3) * t198 - pkin(8) * t200;
t219 = t166 * t165 + (-t153 + t240) * t194 + t234;
t129 = pkin(7) * t207 + t210 * t257;
t218 = t191 * t128 - t129 * t190 + t223;
t217 = t194 * t129 - t191 * t258 + t224;
t154 = t233 * t210;
t216 = t167 * t153 - t154 * t166 + t218;
t215 = t194 * t154 - t165 * t167 + t217;
t133 = -Icges(4,6) * t210 + t207 * t227;
t134 = Icges(4,6) * t207 + t210 * t227;
t135 = -Icges(4,5) * t210 + t207 * t229;
t136 = Icges(4,5) * t207 + t210 * t229;
t162 = Icges(4,2) * t200 + t250;
t163 = Icges(4,1) * t198 + t249;
t214 = (-t134 * t198 + t136 * t200) * t167 + (-t133 * t198 + t135 * t200) * t166 + (-t162 * t198 + t163 * t200) * t194;
t147 = -Icges(3,6) * t210 + t207 * t228;
t148 = Icges(3,6) * t207 + t210 * t228;
t149 = -Icges(3,5) * t210 + t207 * t230;
t150 = Icges(3,5) * t207 + t210 * t230;
t177 = Icges(3,2) * t209 + t252;
t180 = Icges(3,1) * t206 + t251;
t213 = (-t148 * t206 + t150 * t209) * t191 + (-t147 * t206 + t149 * t209) * t190 + (-t177 * t206 + t180 * t209) * t194;
t203 = qJ(4) + qJ(5);
t201 = Icges(2,4) * t210;
t199 = cos(t203);
t197 = sin(t203);
t185 = rSges(2,1) * t210 - rSges(2,2) * t207;
t184 = rSges(2,1) * t207 + rSges(2,2) * t210;
t183 = rSges(3,1) * t206 + rSges(3,2) * t209;
t182 = Icges(2,1) * t210 - t253;
t181 = Icges(2,1) * t207 + t201;
t179 = -Icges(2,2) * t207 + t201;
t178 = Icges(2,2) * t210 + t253;
t173 = -qJD(4) * t200 + t194;
t172 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t171 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t170 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t164 = rSges(4,1) * t198 + rSges(4,2) * t200;
t159 = t200 * t241 + t244;
t158 = -t200 * t243 + t242;
t157 = t200 * t242 - t243;
t156 = -t200 * t244 - t241;
t155 = (-qJD(4) - qJD(5)) * t200 + t194;
t152 = rSges(3,3) * t207 + t210 * t232;
t151 = -rSges(3,3) * t210 + t207 * t232;
t144 = t197 * t207 + t199 * t245;
t143 = -t197 * t245 + t199 * t207;
t142 = -t197 * t210 + t199 * t246;
t141 = -t197 * t246 - t199 * t210;
t138 = rSges(4,3) * t207 + t210 * t231;
t137 = -rSges(4,3) * t210 + t207 * t231;
t127 = -rSges(5,3) * t200 + (rSges(5,1) * t208 - rSges(5,2) * t205) * t198;
t126 = -Icges(5,5) * t200 + (Icges(5,1) * t208 - Icges(5,4) * t205) * t198;
t125 = -Icges(5,6) * t200 + (Icges(5,4) * t208 - Icges(5,2) * t205) * t198;
t124 = -Icges(5,3) * t200 + (Icges(5,5) * t208 - Icges(5,6) * t205) * t198;
t123 = V_base(5) * rSges(2,3) - t184 * t194 + t237;
t122 = t185 * t194 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t120 = t184 * V_base(4) - t185 * V_base(5) + V_base(3);
t119 = -rSges(6,3) * t200 + (rSges(6,1) * t199 - rSges(6,2) * t197) * t198;
t118 = -Icges(6,5) * t200 + (Icges(6,1) * t199 - Icges(6,4) * t197) * t198;
t117 = -Icges(6,6) * t200 + (Icges(6,4) * t199 - Icges(6,2) * t197) * t198;
t116 = -Icges(6,3) * t200 + (Icges(6,5) * t199 - Icges(6,6) * t197) * t198;
t114 = t210 * t238 + t140;
t113 = t207 * t238 + t139;
t111 = -pkin(9) * t200 + t198 * t256;
t109 = rSges(5,1) * t159 + rSges(5,2) * t158 + rSges(5,3) * t247;
t108 = rSges(5,1) * t157 + rSges(5,2) * t156 + rSges(5,3) * t248;
t107 = Icges(5,1) * t159 + Icges(5,4) * t158 + Icges(5,5) * t247;
t106 = Icges(5,1) * t157 + Icges(5,4) * t156 + Icges(5,5) * t248;
t105 = Icges(5,4) * t159 + Icges(5,2) * t158 + Icges(5,6) * t247;
t104 = Icges(5,4) * t157 + Icges(5,2) * t156 + Icges(5,6) * t248;
t103 = Icges(5,5) * t159 + Icges(5,6) * t158 + Icges(5,3) * t247;
t102 = Icges(5,5) * t157 + Icges(5,6) * t156 + Icges(5,3) * t248;
t101 = pkin(4) * t244 + t210 * t222;
t100 = -pkin(4) * t243 + t207 * t222;
t99 = rSges(6,1) * t144 + rSges(6,2) * t143 + rSges(6,3) * t247;
t98 = rSges(6,1) * t142 + rSges(6,2) * t141 + rSges(6,3) * t248;
t97 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t247;
t96 = Icges(6,1) * t142 + Icges(6,4) * t141 + Icges(6,5) * t248;
t95 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t247;
t94 = Icges(6,4) * t142 + Icges(6,2) * t141 + Icges(6,6) * t248;
t93 = Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t247;
t92 = Icges(6,5) * t142 + Icges(6,6) * t141 + Icges(6,3) * t248;
t91 = t183 * t190 + (-t151 - t188) * t194 + t237;
t90 = t152 * t194 - t183 * t191 + t224;
t89 = t151 * t191 - t152 * t190 + t223;
t88 = t164 * t166 + (-t137 + t240) * t194 + t234;
t87 = t138 * t194 - t164 * t167 + t217;
t86 = t137 * t167 - t138 * t166 + t218;
t85 = -t108 * t173 + t127 * t139 + t219;
t84 = t109 * t173 - t127 * t140 + t215;
t83 = t108 * t140 - t109 * t139 + t216;
t82 = -t100 * t173 + t111 * t139 + t113 * t119 - t155 * t98 + t219;
t81 = t101 * t173 - t111 * t140 - t114 * t119 + t155 * t99 + t215;
t80 = t100 * t140 - t101 * t139 - t113 * t99 + t114 * t98 + t216;
t1 = m(1) * (t170 ^ 2 + t171 ^ 2 + t172 ^ 2) / 0.2e1 + m(2) * (t120 ^ 2 + t122 ^ 2 + t123 ^ 2) / 0.2e1 + m(3) * (t89 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t191 * (t207 * t220 + t210 * t213) / 0.2e1 + t190 * (t207 * t213 - t210 * t220) / 0.2e1 + m(4) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + t167 * (t221 * t207 + t214 * t210) / 0.2e1 + t166 * (t214 * t207 - t221 * t210) / 0.2e1 + m(5) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + t140 * ((t103 * t247 + t105 * t158 + t107 * t159) * t140 + (t102 * t247 + t104 * t158 + t106 * t159) * t139 + (t124 * t247 + t125 * t158 + t126 * t159) * t173) / 0.2e1 + t139 * ((t103 * t248 + t105 * t156 + t107 * t157) * t140 + (t102 * t248 + t104 * t156 + t106 * t157) * t139 + (t124 * t248 + t125 * t156 + t126 * t157) * t173) / 0.2e1 + t173 * ((-t102 * t139 - t103 * t140 - t124 * t173) * t200 + ((-t105 * t205 + t107 * t208) * t140 + (-t104 * t205 + t106 * t208) * t139 + (-t125 * t205 + t126 * t208) * t173) * t198) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + t114 * ((t143 * t95 + t144 * t97 + t93 * t247) * t114 + (t143 * t94 + t144 * t96 + t247 * t92) * t113 + (t116 * t247 + t117 * t143 + t118 * t144) * t155) / 0.2e1 + t113 * ((t141 * t95 + t142 * t97 + t248 * t93) * t114 + (t141 * t94 + t142 * t96 + t92 * t248) * t113 + (t116 * t248 + t117 * t141 + t118 * t142) * t155) / 0.2e1 + t155 * ((-t92 * t113 - t93 * t114 - t116 * t155) * t200 + ((-t197 * t95 + t199 * t97) * t114 + (-t197 * t94 + t199 * t96) * t113 + (-t117 * t197 + t118 * t199) * t155) * t198) / 0.2e1 + ((-t178 * t207 + t181 * t210 + Icges(1,4)) * V_base(5) + (-t179 * t207 + t182 * t210 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t178 * t210 + t181 * t207 + Icges(1,2)) * V_base(5) + (t179 * t210 + t182 * t207 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t148 * t209 + t150 * t206) * t191 + (t147 * t209 + t149 * t206) * t190 + (t134 * t200 + t136 * t198) * t167 + (t133 * t200 + t135 * t198) * t166 + (t162 * t200 + t163 * t198 + t177 * t209 + t180 * t206 + Icges(2,3)) * t194) * t194 / 0.2e1 + V_base(4) * t194 * (Icges(2,5) * t210 - Icges(2,6) * t207) + V_base(5) * t194 * (Icges(2,5) * t207 + Icges(2,6) * t210) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
