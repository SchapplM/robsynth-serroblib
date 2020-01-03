% Calculate kinetic energy for
% S5RPRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:38:10
% EndTime: 2020-01-03 11:38:12
% DurationCPUTime: 2.10s
% Computational Cost: add. (1262->241), mult. (946->325), div. (0->0), fcn. (726->10), ass. (0->125)
t238 = Icges(4,3) + Icges(5,3);
t166 = qJ(3) + pkin(9);
t155 = sin(t166);
t157 = cos(t166);
t169 = sin(qJ(3));
t171 = cos(qJ(3));
t237 = Icges(4,5) * t171 + Icges(5,5) * t157 - Icges(4,6) * t169 - Icges(5,6) * t155;
t167 = qJ(1) + pkin(8);
t156 = sin(t167);
t158 = cos(t167);
t215 = Icges(4,4) * t171;
t192 = -Icges(4,2) * t169 + t215;
t100 = -Icges(4,6) * t158 + t156 * t192;
t101 = -Icges(4,6) * t156 - t158 * t192;
t216 = Icges(4,4) * t169;
t195 = Icges(4,1) * t171 - t216;
t102 = -Icges(4,5) * t158 + t156 * t195;
t103 = -Icges(4,5) * t156 - t158 * t195;
t214 = Icges(5,4) * t155;
t120 = Icges(5,2) * t157 + t214;
t213 = Icges(5,4) * t157;
t123 = Icges(5,1) * t155 + t213;
t135 = -qJD(3) * t156 + V_base(6);
t136 = -qJD(3) * t158 + V_base(5);
t140 = Icges(4,2) * t171 + t216;
t143 = Icges(4,1) * t169 + t215;
t160 = V_base(4) + qJD(1);
t191 = -Icges(5,2) * t155 + t213;
t91 = -Icges(5,6) * t158 + t156 * t191;
t92 = -Icges(5,6) * t156 - t158 * t191;
t194 = Icges(5,1) * t157 - t214;
t93 = -Icges(5,5) * t158 + t156 * t194;
t94 = -Icges(5,5) * t156 - t158 * t194;
t236 = t135 * (t101 * t169 - t103 * t171 + t155 * t92 - t157 * t94) + t136 * (t100 * t169 - t102 * t171 + t155 * t91 - t157 * t93) + t160 * (t120 * t155 - t123 * t157 + t140 * t169 - t143 * t171);
t233 = (-Icges(4,5) * t169 - Icges(5,5) * t155 - Icges(4,6) * t171 - Icges(5,6) * t157) * t160 + (-t237 * t156 + t238 * t158) * t136 + (t238 * t156 + t237 * t158) * t135;
t208 = -qJD(3) - qJD(5);
t109 = t156 * t208 + V_base(6);
t110 = t158 * t208 + V_base(5);
t159 = qJ(5) + t166;
t153 = cos(t159);
t152 = sin(t159);
t212 = Icges(6,4) * t152;
t113 = Icges(6,2) * t153 + t212;
t211 = Icges(6,4) * t153;
t114 = Icges(6,1) * t152 + t211;
t190 = -Icges(6,2) * t152 + t211;
t83 = -Icges(6,6) * t158 + t156 * t190;
t84 = -Icges(6,6) * t156 - t158 * t190;
t193 = Icges(6,1) * t153 - t212;
t85 = -Icges(6,5) * t158 + t156 * t193;
t86 = -Icges(6,5) * t156 - t158 * t193;
t229 = (t113 * t152 - t114 * t153) * t160 + (t152 * t83 - t153 * t85) * t110 + (t152 * t84 - t153 * t86) * t109;
t225 = pkin(1) * t160;
t224 = pkin(3) * t169;
t223 = pkin(4) * t155;
t222 = t171 * pkin(3);
t221 = -pkin(5) - qJ(2);
t130 = -pkin(2) * t158 - pkin(6) * t156;
t80 = -qJ(4) * t156 - t158 * t222;
t219 = -t130 - t80;
t172 = cos(qJ(1));
t218 = Icges(2,4) * t172;
t217 = Icges(3,4) * t158;
t210 = pkin(4) * t157;
t170 = sin(qJ(1));
t207 = t170 * t225 + V_base(3);
t206 = V_base(6) * pkin(5) + V_base(2);
t203 = V_base(6) * qJ(2) + t172 * t225 + t206;
t202 = rSges(4,1) * t171 - rSges(4,2) * t169;
t201 = rSges(5,1) * t157 - rSges(5,2) * t155;
t200 = rSges(6,1) * t153 - rSges(6,2) * t152;
t187 = Icges(6,5) * t153 - Icges(6,6) * t152;
t181 = -(-Icges(6,3) * t156 - t158 * t187) * t109 - (-Icges(6,3) * t158 + t156 * t187) * t110 - (Icges(6,5) * t152 + Icges(6,6) * t153) * t160;
t178 = -qJD(4) * t158 + t135 * t224 + t203;
t129 = pkin(2) * t156 - pkin(6) * t158;
t177 = t160 * t129 + t221 * V_base(5) + t207;
t176 = qJD(2) + V_base(1) + (-V_base(6) * t170 - t172 * V_base(5)) * pkin(1);
t79 = -qJ(4) * t158 + t156 * t222;
t175 = -qJD(4) * t156 + t160 * t79 + t177;
t174 = -V_base(6) * t129 + V_base(5) * t130 + t176;
t173 = t136 * t80 + t174;
t162 = Icges(2,4) * t170;
t151 = Icges(3,4) * t156;
t148 = -rSges(2,1) * t172 + t170 * rSges(2,2);
t147 = t170 * rSges(2,1) + rSges(2,2) * t172;
t146 = rSges(4,1) * t169 + rSges(4,2) * t171;
t145 = -Icges(2,1) * t172 + t162;
t144 = Icges(2,1) * t170 + t218;
t142 = Icges(2,2) * t170 - t218;
t141 = Icges(2,2) * t172 + t162;
t134 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t133 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t132 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t128 = -rSges(3,1) * t158 + rSges(3,2) * t156;
t127 = rSges(3,1) * t156 + rSges(3,2) * t158;
t126 = rSges(5,1) * t155 + rSges(5,2) * t157;
t125 = -Icges(3,1) * t158 + t151;
t124 = Icges(3,1) * t156 + t217;
t122 = Icges(3,2) * t156 - t217;
t121 = Icges(3,2) * t158 + t151;
t115 = rSges(6,1) * t152 + rSges(6,2) * t153;
t107 = -rSges(4,3) * t156 - t158 * t202;
t106 = -rSges(4,3) * t158 + t156 * t202;
t105 = V_base(6) * rSges(2,3) - t148 * t160 + t206;
t104 = t147 * t160 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t97 = -t147 * V_base(6) + t148 * V_base(5) + V_base(1);
t96 = -rSges(5,3) * t156 - t158 * t201;
t95 = -rSges(5,3) * t158 + t156 * t201;
t88 = -rSges(6,3) * t156 - t158 * t200;
t87 = -rSges(6,3) * t158 + t156 * t200;
t77 = V_base(6) * rSges(3,3) - t128 * t160 + t203;
t76 = t127 * t160 + (-rSges(3,3) + t221) * V_base(5) + t207;
t74 = -pkin(7) * t156 - t158 * t210;
t73 = -pkin(7) * t158 + t156 * t210;
t72 = -V_base(6) * t127 + V_base(5) * t128 + t176;
t71 = t135 * t146 + (-t107 - t130) * t160 + t203;
t70 = t106 * t160 - t136 * t146 + t177;
t69 = -t135 * t106 + t136 * t107 + t174;
t68 = t126 * t135 + (-t96 + t219) * t160 + t178;
t67 = t160 * t95 + (-t126 - t224) * t136 + t175;
t66 = t136 * t96 + (-t79 - t95) * t135 + t173;
t65 = t135 * t223 + t109 * t115 + (-t74 - t88 + t219) * t160 + t178;
t64 = -t110 * t115 + (t73 + t87) * t160 + (-t223 - t224) * t136 + t175;
t63 = -t109 * t87 + t110 * t88 + t136 * t74 + (-t73 - t79) * t135 + t173;
t1 = m(1) * (t132 ^ 2 + t133 ^ 2 + t134 ^ 2) / 0.2e1 + m(2) * (t104 ^ 2 + t105 ^ 2 + t97 ^ 2) / 0.2e1 + m(3) * (t72 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + t110 * (-t229 * t156 + t181 * t158) / 0.2e1 + t109 * (t181 * t156 + t158 * t229) / 0.2e1 + (t233 * t156 + t236 * t158) * t135 / 0.2e1 + (-t236 * t156 + t233 * t158) * t136 / 0.2e1 + ((t122 * t158 + t125 * t156 + t142 * t172 + t170 * t145 + Icges(1,6)) * V_base(6) + (t121 * t158 + t124 * t156 + t141 * t172 + t170 * t144 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + ((t122 * t156 - t125 * t158 + t170 * t142 - t145 * t172 + Icges(1,3)) * V_base(6) + (t121 * t156 - t124 * t158 + t170 * t141 - t144 * t172 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + ((t152 * t85 + t153 * t83) * t110 + (t152 * t86 + t153 * t84) * t109 + (t100 * t171 + t102 * t169 + t155 * t93 + t157 * t91) * t136 + (t101 * t171 + t103 * t169 + t155 * t94 + t157 * t92) * t135 + (t113 * t153 + t114 * t152 + t120 * t157 + t123 * t155 + t140 * t171 + t143 * t169 + Icges(2,3) + Icges(3,3)) * t160) * t160 / 0.2e1 + t160 * V_base(6) * (-Icges(2,5) * t172 - Icges(3,5) * t158 + Icges(2,6) * t170 + Icges(3,6) * t156) + t160 * V_base(5) * (Icges(2,5) * t170 + Icges(3,5) * t156 + Icges(2,6) * t172 + Icges(3,6) * t158) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
