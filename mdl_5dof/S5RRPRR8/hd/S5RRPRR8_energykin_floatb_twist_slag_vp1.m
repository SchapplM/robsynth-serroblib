% Calculate kinetic energy for
% S5RRPRR8
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR8_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR8_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR8_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR8_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR8_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:16:49
% EndTime: 2019-12-31 20:16:51
% DurationCPUTime: 2.17s
% Computational Cost: add. (1410->280), mult. (1390->407), div. (0->0), fcn. (1228->10), ass. (0->144)
t261 = Icges(3,3) + Icges(4,3);
t195 = qJ(2) + pkin(9);
t184 = sin(t195);
t185 = cos(t195);
t198 = sin(qJ(2));
t201 = cos(qJ(2));
t260 = Icges(3,5) * t201 + Icges(4,5) * t185 - Icges(3,6) * t198 - Icges(4,6) * t184;
t199 = sin(qJ(1));
t202 = cos(qJ(1));
t246 = Icges(4,4) * t185;
t220 = -Icges(4,2) * t184 + t246;
t123 = -Icges(4,6) * t202 + t220 * t199;
t124 = Icges(4,6) * t199 + t220 * t202;
t247 = Icges(4,4) * t184;
t223 = Icges(4,1) * t185 - t247;
t125 = -Icges(4,5) * t202 + t223 * t199;
t126 = Icges(4,5) * t199 + t223 * t202;
t248 = Icges(3,4) * t201;
t221 = -Icges(3,2) * t198 + t248;
t135 = -Icges(3,6) * t202 + t221 * t199;
t136 = Icges(3,6) * t199 + t221 * t202;
t249 = Icges(3,4) * t198;
t224 = Icges(3,1) * t201 - t249;
t137 = -Icges(3,5) * t202 + t224 * t199;
t138 = Icges(3,5) * t199 + t224 * t202;
t153 = Icges(4,2) * t185 + t247;
t154 = Icges(4,1) * t184 + t246;
t168 = Icges(3,2) * t201 + t249;
t171 = Icges(3,1) * t198 + t248;
t179 = -qJD(2) * t202 + V_base(5);
t180 = qJD(2) * t199 + V_base(4);
t187 = V_base(6) + qJD(1);
t259 = (-t153 * t184 + t154 * t185 - t168 * t198 + t171 * t201) * t187 + (-t124 * t184 + t126 * t185 - t136 * t198 + t138 * t201) * t180 + (-t123 * t184 + t125 * t185 - t135 * t198 + t137 * t201) * t179;
t258 = (Icges(3,5) * t198 + Icges(4,5) * t184 + Icges(3,6) * t201 + Icges(4,6) * t185) * t187 + (t261 * t199 + t260 * t202) * t180 + (t260 * t199 - t261 * t202) * t179;
t254 = pkin(2) * t198;
t253 = pkin(3) * t184;
t252 = t201 * pkin(2);
t250 = Icges(2,4) * t199;
t186 = qJ(4) + t195;
t181 = sin(t186);
t245 = Icges(5,4) * t181;
t182 = cos(t186);
t244 = Icges(5,4) * t182;
t243 = t181 * t199;
t242 = t181 * t202;
t197 = sin(qJ(5));
t241 = t197 * t202;
t240 = t199 * t197;
t200 = cos(qJ(5));
t239 = t199 * t200;
t238 = t200 * t202;
t119 = -qJ(3) * t202 + t252 * t199;
t177 = t199 * pkin(1) - pkin(6) * t202;
t237 = -t119 - t177;
t236 = pkin(3) * t185;
t234 = qJD(5) * t181;
t233 = V_base(5) * pkin(5) + V_base(1);
t100 = -pkin(7) * t202 + t236 * t199;
t230 = -t100 + t237;
t158 = qJD(4) * t199 + t180;
t229 = qJD(3) * t199 + t179 * t254 + t233;
t228 = pkin(4) * t182 + pkin(8) * t181;
t227 = rSges(3,1) * t201 - rSges(3,2) * t198;
t226 = rSges(4,1) * t185 - rSges(4,2) * t184;
t225 = rSges(5,1) * t182 - rSges(5,2) * t181;
t222 = Icges(5,1) * t182 - t245;
t219 = -Icges(5,2) * t181 + t244;
t216 = Icges(5,5) * t182 - Icges(5,6) * t181;
t215 = t179 * t253 + t229;
t178 = pkin(1) * t202 + t199 * pkin(6);
t214 = -V_base(4) * pkin(5) + t187 * t178 + V_base(2);
t213 = V_base(4) * t177 - t178 * V_base(5) + V_base(3);
t157 = V_base(5) + (-qJD(2) - qJD(4)) * t202;
t212 = t180 * t119 + t213;
t211 = (-Icges(5,3) * t202 + t216 * t199) * t157 + (Icges(5,3) * t199 + t216 * t202) * t158 + (Icges(5,5) * t181 + Icges(5,6) * t182) * t187;
t120 = qJ(3) * t199 + t252 * t202;
t208 = -qJD(3) * t202 + t187 * t120 + t214;
t101 = pkin(7) * t199 + t236 * t202;
t207 = t180 * t100 + (-t101 - t120) * t179 + t212;
t206 = t187 * t101 + (-t253 - t254) * t180 + t208;
t111 = -Icges(5,6) * t202 + t219 * t199;
t112 = Icges(5,6) * t199 + t219 * t202;
t113 = -Icges(5,5) * t202 + t222 * t199;
t114 = Icges(5,5) * t199 + t222 * t202;
t146 = Icges(5,2) * t182 + t245;
t147 = Icges(5,1) * t181 + t244;
t205 = (-t112 * t181 + t114 * t182) * t158 + (-t111 * t181 + t113 * t182) * t157 + (-t146 * t181 + t147 * t182) * t187;
t191 = Icges(2,4) * t202;
t176 = rSges(2,1) * t202 - t199 * rSges(2,2);
t175 = t199 * rSges(2,1) + rSges(2,2) * t202;
t174 = rSges(3,1) * t198 + rSges(3,2) * t201;
t173 = Icges(2,1) * t202 - t250;
t172 = Icges(2,1) * t199 + t191;
t170 = -Icges(2,2) * t199 + t191;
t169 = Icges(2,2) * t202 + t250;
t164 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t163 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t162 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t156 = -qJD(5) * t182 + t187;
t155 = rSges(4,1) * t184 + rSges(4,2) * t185;
t149 = pkin(4) * t181 - pkin(8) * t182;
t148 = rSges(5,1) * t181 + rSges(5,2) * t182;
t144 = t182 * t238 + t240;
t143 = -t182 * t241 + t239;
t142 = t182 * t239 - t241;
t141 = -t182 * t240 - t238;
t140 = t199 * rSges(3,3) + t227 * t202;
t139 = -rSges(3,3) * t202 + t227 * t199;
t132 = t228 * t202;
t131 = t228 * t199;
t130 = t199 * rSges(4,3) + t226 * t202;
t129 = -rSges(4,3) * t202 + t226 * t199;
t128 = t202 * t234 + t158;
t127 = t199 * t234 + t157;
t118 = t199 * rSges(5,3) + t225 * t202;
t117 = -rSges(5,3) * t202 + t225 * t199;
t116 = V_base(5) * rSges(2,3) - t175 * t187 + t233;
t115 = t176 * t187 + V_base(2) + (-rSges(2,3) - pkin(5)) * V_base(4);
t108 = t175 * V_base(4) - t176 * V_base(5) + V_base(3);
t107 = -rSges(6,3) * t182 + (rSges(6,1) * t200 - rSges(6,2) * t197) * t181;
t105 = -Icges(6,5) * t182 + (Icges(6,1) * t200 - Icges(6,4) * t197) * t181;
t104 = -Icges(6,6) * t182 + (Icges(6,4) * t200 - Icges(6,2) * t197) * t181;
t103 = -Icges(6,3) * t182 + (Icges(6,5) * t200 - Icges(6,6) * t197) * t181;
t97 = t144 * rSges(6,1) + t143 * rSges(6,2) + rSges(6,3) * t242;
t96 = rSges(6,1) * t142 + rSges(6,2) * t141 + rSges(6,3) * t243;
t95 = Icges(6,1) * t144 + Icges(6,4) * t143 + Icges(6,5) * t242;
t94 = Icges(6,1) * t142 + Icges(6,4) * t141 + Icges(6,5) * t243;
t93 = Icges(6,4) * t144 + Icges(6,2) * t143 + Icges(6,6) * t242;
t92 = Icges(6,4) * t142 + Icges(6,2) * t141 + Icges(6,6) * t243;
t91 = Icges(6,5) * t144 + Icges(6,6) * t143 + Icges(6,3) * t242;
t90 = Icges(6,5) * t142 + Icges(6,6) * t141 + Icges(6,3) * t243;
t89 = t174 * t179 + (-t139 - t177) * t187 + t233;
t88 = t140 * t187 - t174 * t180 + t214;
t87 = t139 * t180 - t140 * t179 + t213;
t86 = t155 * t179 + (-t129 + t237) * t187 + t229;
t85 = t187 * t130 + (-t155 - t254) * t180 + t208;
t84 = t129 * t180 + (-t120 - t130) * t179 + t212;
t83 = t148 * t157 + (-t117 + t230) * t187 + t215;
t82 = t187 * t118 - t158 * t148 + t206;
t81 = t117 * t158 - t118 * t157 + t207;
t80 = t107 * t127 + t149 * t157 - t156 * t96 + (-t131 + t230) * t187 + t215;
t79 = -t128 * t107 + t187 * t132 - t158 * t149 + t156 * t97 + t206;
t78 = -t127 * t97 + t128 * t96 + t131 * t158 - t132 * t157 + t207;
t1 = m(1) * (t162 ^ 2 + t163 ^ 2 + t164 ^ 2) / 0.2e1 + m(2) * (t108 ^ 2 + t115 ^ 2 + t116 ^ 2) / 0.2e1 + m(3) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(4) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(5) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + t158 * (t211 * t199 + t205 * t202) / 0.2e1 + t157 * (t205 * t199 - t211 * t202) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + t128 * ((t143 * t93 + t144 * t95 + t91 * t242) * t128 + (t143 * t92 + t144 * t94 + t90 * t242) * t127 + (t103 * t242 + t143 * t104 + t144 * t105) * t156) / 0.2e1 + t127 * ((t141 * t93 + t142 * t95 + t91 * t243) * t128 + (t141 * t92 + t142 * t94 + t90 * t243) * t127 + (t103 * t243 + t104 * t141 + t105 * t142) * t156) / 0.2e1 + t156 * ((-t103 * t156 - t90 * t127 - t91 * t128) * t182 + ((-t197 * t93 + t200 * t95) * t128 + (-t197 * t92 + t200 * t94) * t127 + (-t104 * t197 + t105 * t200) * t156) * t181) / 0.2e1 + (t259 * t199 - t258 * t202) * t179 / 0.2e1 + (t258 * t199 + t259 * t202) * t180 / 0.2e1 + ((-t199 * t169 + t172 * t202 + Icges(1,4)) * V_base(5) + (-t199 * t170 + t173 * t202 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t169 * t202 + t199 * t172 + Icges(1,2)) * V_base(5) + (t170 * t202 + t199 * t173 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t112 * t182 + t114 * t181) * t158 + (t111 * t182 + t113 * t181) * t157 + (t124 * t185 + t126 * t184 + t136 * t201 + t138 * t198) * t180 + (t123 * t185 + t125 * t184 + t135 * t201 + t137 * t198) * t179 + (t146 * t182 + t147 * t181 + t153 * t185 + t154 * t184 + t168 * t201 + t171 * t198 + Icges(2,3)) * t187) * t187 / 0.2e1 + t187 * V_base(4) * (Icges(2,5) * t202 - Icges(2,6) * t199) + t187 * V_base(5) * (Icges(2,5) * t199 + Icges(2,6) * t202) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
