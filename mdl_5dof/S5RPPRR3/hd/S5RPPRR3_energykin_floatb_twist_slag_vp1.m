% Calculate kinetic energy for
% S5RPPRR3
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:29
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPPRR3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:27:23
% EndTime: 2020-01-03 11:27:25
% DurationCPUTime: 2.19s
% Computational Cost: add. (1226->250), mult. (928->338), div. (0->0), fcn. (708->10), ass. (0->133)
t167 = qJ(1) + pkin(8);
t156 = sin(t167);
t158 = cos(t167);
t171 = sin(qJ(1));
t172 = cos(qJ(1));
t236 = Icges(2,5) * t171 + Icges(3,5) * t156 + Icges(2,6) * t172 + Icges(3,6) * t158;
t235 = -Icges(2,5) * t172 - Icges(3,5) * t158 + Icges(2,6) * t171 + Icges(3,6) * t156;
t168 = sin(pkin(9));
t169 = cos(pkin(9));
t219 = Icges(4,4) * t169;
t191 = -Icges(4,2) * t168 + t219;
t100 = -Icges(4,6) * t158 + t156 * t191;
t220 = Icges(4,4) * t168;
t194 = Icges(4,1) * t169 - t220;
t102 = -Icges(4,5) * t158 + t156 * t194;
t221 = Icges(3,4) * t158;
t234 = Icges(3,1) * t156 - t100 * t168 + t102 * t169 + t221;
t101 = -Icges(4,6) * t156 - t158 * t191;
t103 = -Icges(4,5) * t156 - t158 * t194;
t151 = Icges(3,4) * t156;
t233 = -Icges(3,1) * t158 - t101 * t168 + t103 * t169 + t151;
t211 = -qJD(4) - qJD(5);
t109 = t156 * t211 + V_base(6);
t110 = t158 * t211 + V_base(5);
t166 = pkin(9) + qJ(4);
t159 = qJ(5) + t166;
t153 = cos(t159);
t152 = sin(t159);
t216 = Icges(6,4) * t152;
t113 = Icges(6,2) * t153 + t216;
t215 = Icges(6,4) * t153;
t114 = Icges(6,1) * t152 + t215;
t160 = V_base(4) + qJD(1);
t189 = -Icges(6,2) * t152 + t215;
t83 = -Icges(6,6) * t158 + t156 * t189;
t84 = -Icges(6,6) * t156 - t158 * t189;
t192 = Icges(6,1) * t153 - t216;
t85 = -Icges(6,5) * t158 + t156 * t192;
t86 = -Icges(6,5) * t156 - t158 * t192;
t232 = (t113 * t152 - t114 * t153) * t160 + (t152 * t83 - t153 * t85) * t110 + (t152 * t84 - t153 * t86) * t109;
t157 = cos(t166);
t155 = sin(t166);
t218 = Icges(5,4) * t155;
t119 = Icges(5,2) * t157 + t218;
t217 = Icges(5,4) * t157;
t122 = Icges(5,1) * t155 + t217;
t137 = -qJD(4) * t156 + V_base(6);
t138 = -qJD(4) * t158 + V_base(5);
t190 = -Icges(5,2) * t155 + t217;
t91 = -Icges(5,6) * t158 + t156 * t190;
t92 = -Icges(5,6) * t156 - t158 * t190;
t193 = Icges(5,1) * t157 - t218;
t93 = -Icges(5,5) * t158 + t156 * t193;
t94 = -Icges(5,5) * t156 - t158 * t193;
t231 = (t119 * t155 - t122 * t157) * t160 + (t155 * t91 - t157 * t93) * t138 + (t155 * t92 - t157 * t94) * t137;
t229 = pkin(1) * t171;
t228 = pkin(1) * t172;
t227 = pkin(3) * t168;
t226 = pkin(4) * t155;
t225 = t169 * pkin(3);
t224 = -pkin(5) - qJ(2);
t128 = -pkin(2) * t158 - qJ(3) * t156;
t80 = -pkin(6) * t156 - t158 * t225;
t223 = -t128 - t80;
t222 = Icges(2,4) * t172;
t213 = pkin(4) * t157;
t210 = t160 * t229 + V_base(3);
t209 = V_base(6) * pkin(5) + V_base(2);
t206 = qJD(2) + V_base(1);
t205 = t172 * V_base(5);
t126 = pkin(2) * t156 - qJ(3) * t158;
t204 = -t126 - t229;
t203 = V_base(5) * t128 + t206;
t202 = V_base(6) * qJ(2) + t160 * t228 + t209;
t201 = rSges(4,1) * t169 - rSges(4,2) * t168;
t200 = rSges(5,1) * t157 - rSges(5,2) * t155;
t199 = rSges(6,1) * t153 - rSges(6,2) * t152;
t188 = Icges(4,5) * t169 - Icges(4,6) * t168;
t187 = Icges(5,5) * t157 - Icges(5,6) * t155;
t186 = Icges(6,5) * t153 - Icges(6,6) * t152;
t135 = Icges(4,2) * t169 + t220;
t136 = Icges(4,1) * t168 + t219;
t181 = t135 * t168 - t136 * t169;
t180 = -qJD(3) * t156 + t160 * t126 + t210;
t179 = -qJD(3) * t158 + t202;
t178 = -t109 * (-Icges(6,3) * t156 - t158 * t186) - t110 * (-Icges(6,3) * t158 + t156 * t186) - (Icges(6,5) * t152 + Icges(6,6) * t153) * t160;
t177 = -(Icges(5,5) * t155 + Icges(5,6) * t157) * t160 - t137 * (-Icges(5,3) * t156 - t158 * t187) - t138 * (-Icges(5,3) * t158 + t156 * t187);
t176 = V_base(6) * t227 + t179;
t175 = -(Icges(4,5) * t168 + Icges(4,6) * t169) * t160 - (-Icges(4,3) * t158 + t156 * t188) * V_base(5) - (-Icges(4,3) * t156 - t158 * t188) * V_base(6);
t79 = -pkin(6) * t158 + t156 * t225;
t174 = t160 * t79 + (t224 - t227) * V_base(5) + t180;
t173 = V_base(5) * t80 + (t204 - t79) * V_base(6) - pkin(1) * t205 + t203;
t163 = Icges(2,4) * t171;
t147 = -rSges(2,1) * t172 + t171 * rSges(2,2);
t146 = t171 * rSges(2,1) + rSges(2,2) * t172;
t145 = -Icges(2,1) * t172 + t163;
t144 = Icges(2,1) * t171 + t222;
t143 = Icges(2,2) * t171 - t222;
t142 = Icges(2,2) * t172 + t163;
t139 = rSges(4,1) * t168 + rSges(4,2) * t169;
t133 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t132 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t131 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t129 = -rSges(3,1) * t158 + rSges(3,2) * t156;
t127 = rSges(3,1) * t156 + rSges(3,2) * t158;
t125 = rSges(5,1) * t155 + rSges(5,2) * t157;
t121 = Icges(3,2) * t156 - t221;
t120 = Icges(3,2) * t158 + t151;
t115 = rSges(6,1) * t152 + rSges(6,2) * t153;
t107 = V_base(6) * rSges(2,3) - t147 * t160 + t209;
t106 = t146 * t160 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t105 = -rSges(4,3) * t156 - t158 * t201;
t104 = -rSges(4,3) * t158 + t156 * t201;
t97 = -t146 * V_base(6) + t147 * V_base(5) + V_base(1);
t96 = -rSges(5,3) * t156 - t158 * t200;
t95 = -rSges(5,3) * t158 + t156 * t200;
t88 = -rSges(6,3) * t156 - t158 * t199;
t87 = -rSges(6,3) * t158 + t156 * t199;
t76 = V_base(6) * rSges(3,3) - t129 * t160 + t202;
t75 = t127 * t160 + (-rSges(3,3) + t224) * V_base(5) + t210;
t74 = -pkin(7) * t156 - t158 * t213;
t73 = -pkin(7) * t158 + t156 * t213;
t72 = -V_base(6) * t127 + V_base(5) * t129 + (-t171 * V_base(6) - t205) * pkin(1) + t206;
t71 = t139 * V_base(6) + (-t105 - t128) * t160 + t179;
t70 = t104 * t160 + (-t139 + t224) * V_base(5) + t180;
t69 = (t105 - t228) * V_base(5) + (-t104 + t204) * V_base(6) + t203;
t68 = t125 * t137 + (-t96 + t223) * t160 + t176;
t67 = -t125 * t138 + t160 * t95 + t174;
t66 = -t137 * t95 + t138 * t96 + t173;
t65 = t137 * t226 + t109 * t115 + (-t74 - t88 + t223) * t160 + t176;
t64 = -t138 * t226 - t110 * t115 + (t73 + t87) * t160 + t174;
t63 = -t109 * t87 + t110 * t88 - t137 * t73 + t138 * t74 + t173;
t1 = m(1) * (t131 ^ 2 + t132 ^ 2 + t133 ^ 2) / 0.2e1 + m(2) * (t106 ^ 2 + t107 ^ 2 + t97 ^ 2) / 0.2e1 + m(3) * (t72 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + t138 * (-t231 * t156 + t177 * t158) / 0.2e1 + t137 * (t177 * t156 + t231 * t158) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t110 * (-t232 * t156 + t178 * t158) / 0.2e1 + t109 * (t178 * t156 + t232 * t158) / 0.2e1 + (t175 * t158 + (-t181 * t156 + t236) * t160 + (t121 * t158 + t143 * t172 + t171 * t145 + t233 * t156 + Icges(1,6)) * V_base(6) + (t158 * t120 + t172 * t142 + t171 * t144 + t234 * t156 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t175 * t156 + (t181 * t158 + t235) * t160 + (t156 * t121 + t171 * t143 - t172 * t145 - t233 * t158 + Icges(1,3)) * V_base(6) + (t120 * t156 + t171 * t142 - t144 * t172 - t234 * t158 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + ((t155 * t93 + t157 * t91) * t138 + (t155 * t94 + t157 * t92) * t137 + (t152 * t85 + t153 * t83) * t110 + (t152 * t86 + t153 * t84) * t109 + (t101 * t169 + t103 * t168 + t235) * V_base(6) + (t100 * t169 + t102 * t168 + t236) * V_base(5) + (t153 * t113 + t152 * t114 + t157 * t119 + t155 * t122 + t169 * t135 + t168 * t136 + Icges(2,3) + Icges(3,3)) * t160) * t160 / 0.2e1 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
