% Calculate kinetic energy for
% S5RPRPR3
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
% Datum: 2020-01-03 11:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RPRPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:35:53
% EndTime: 2020-01-03 11:35:55
% DurationCPUTime: 2.08s
% Computational Cost: add. (1200->256), mult. (946->342), div. (0->0), fcn. (780->10), ass. (0->123)
t156 = qJ(1) + pkin(8);
t150 = qJ(3) + t156;
t145 = sin(t150);
t157 = sin(pkin(9));
t158 = cos(pkin(9));
t146 = cos(t150);
t194 = Icges(4,4) * t146;
t192 = Icges(5,4) * t158;
t170 = -Icges(5,2) * t157 + t192;
t81 = -Icges(5,6) * t146 + t170 * t145;
t193 = Icges(5,4) * t157;
t171 = Icges(5,1) * t158 - t193;
t83 = -Icges(5,5) * t146 + t171 * t145;
t206 = Icges(4,1) * t145 - t157 * t81 + t158 * t83 + t194;
t143 = Icges(4,4) * t145;
t82 = -Icges(5,6) * t145 - t170 * t146;
t84 = -Icges(5,5) * t145 - t171 * t146;
t205 = -Icges(4,1) * t146 - t157 * t82 + t158 * t84 + t143;
t160 = sin(qJ(1));
t200 = pkin(1) * t160;
t162 = cos(qJ(1));
t199 = pkin(1) * t162;
t151 = V_base(4) + qJD(1);
t198 = pkin(2) * t151;
t197 = -pkin(5) - qJ(2);
t196 = Icges(2,4) * t162;
t149 = cos(t156);
t195 = Icges(3,4) * t149;
t191 = t145 * t157;
t190 = t146 * t157;
t159 = sin(qJ(5));
t189 = t158 * t159;
t161 = cos(qJ(5));
t188 = t158 * t161;
t187 = qJD(5) * t157;
t186 = -pkin(6) + t197;
t185 = t151 * t200 + V_base(3);
t184 = V_base(6) * pkin(5) + V_base(2);
t181 = qJD(2) + V_base(1);
t148 = sin(t156);
t180 = t148 * t198 + t185;
t110 = -pkin(3) * t146 - qJ(4) * t145;
t179 = V_base(5) * t110 + t181;
t147 = qJD(3) + t151;
t178 = V_base(6) * qJ(2) + t151 * t199 + t184;
t177 = -pkin(2) * t148 - t200;
t176 = -pkin(2) * t149 - t199;
t175 = pkin(4) * t158 + pkin(7) * t157;
t174 = rSges(5,1) * t158 - rSges(5,2) * t157;
t169 = Icges(5,5) * t158 - Icges(5,6) * t157;
t129 = Icges(5,2) * t158 + t193;
t130 = Icges(5,1) * t157 + t192;
t168 = t129 * t157 - t130 * t158;
t108 = pkin(3) * t145 - qJ(4) * t146;
t167 = -t108 + t177;
t166 = V_base(6) * pkin(6) + t149 * t198 + t178;
t165 = -qJD(4) * t145 + t147 * t108 + t180;
t164 = -qJD(4) * t146 + t166;
t163 = -(Icges(5,5) * t157 + Icges(5,6) * t158) * t147 - (-Icges(5,3) * t146 + t169 * t145) * V_base(5) - (-Icges(5,3) * t145 - t169 * t146) * V_base(6);
t153 = Icges(2,4) * t160;
t144 = Icges(3,4) * t148;
t140 = -rSges(2,1) * t162 + t160 * rSges(2,2);
t139 = t160 * rSges(2,1) + rSges(2,2) * t162;
t138 = -Icges(2,1) * t162 + t153;
t137 = Icges(2,1) * t160 + t196;
t136 = Icges(2,2) * t160 - t196;
t135 = Icges(2,2) * t162 + t153;
t132 = pkin(4) * t157 - pkin(7) * t158;
t131 = rSges(5,1) * t157 + rSges(5,2) * t158;
t125 = -qJD(5) * t158 + t147;
t124 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t123 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t122 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t121 = -rSges(3,1) * t149 + rSges(3,2) * t148;
t120 = rSges(3,1) * t148 + rSges(3,2) * t149;
t119 = -Icges(3,1) * t149 + t144;
t118 = Icges(3,1) * t148 + t195;
t117 = Icges(3,2) * t148 - t195;
t116 = Icges(3,2) * t149 + t144;
t113 = t145 * t187 + V_base(5);
t112 = -t146 * t187 + V_base(6);
t111 = -rSges(4,1) * t146 + rSges(4,2) * t145;
t109 = rSges(4,1) * t145 + rSges(4,2) * t146;
t105 = Icges(4,2) * t145 - t194;
t104 = Icges(4,2) * t146 + t143;
t103 = -Icges(4,5) * t146 + Icges(4,6) * t145;
t102 = Icges(4,5) * t145 + Icges(4,6) * t146;
t100 = -rSges(6,3) * t158 + (rSges(6,1) * t161 - rSges(6,2) * t159) * t157;
t99 = -Icges(6,5) * t158 + (Icges(6,1) * t161 - Icges(6,4) * t159) * t157;
t98 = -Icges(6,6) * t158 + (Icges(6,4) * t161 - Icges(6,2) * t159) * t157;
t97 = -Icges(6,3) * t158 + (Icges(6,5) * t161 - Icges(6,6) * t159) * t157;
t96 = t175 * t146;
t95 = t175 * t145;
t94 = -t145 * t159 - t146 * t188;
t93 = -t145 * t161 + t146 * t189;
t92 = t145 * t188 - t146 * t159;
t91 = -t145 * t189 - t146 * t161;
t89 = V_base(6) * rSges(2,3) - t140 * t151 + t184;
t88 = t139 * t151 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t87 = -t139 * V_base(6) + t140 * V_base(5) + V_base(1);
t86 = -rSges(5,3) * t145 - t174 * t146;
t85 = -rSges(5,3) * t146 + t174 * t145;
t78 = V_base(6) * rSges(3,3) - t121 * t151 + t178;
t77 = t120 * t151 + (-rSges(3,3) + t197) * V_base(5) + t185;
t76 = -V_base(6) * t120 + V_base(5) * t121 + (-t160 * V_base(6) - t162 * V_base(5)) * pkin(1) + t181;
t75 = rSges(6,1) * t94 + rSges(6,2) * t93 - rSges(6,3) * t190;
t74 = rSges(6,1) * t92 + rSges(6,2) * t91 + rSges(6,3) * t191;
t73 = Icges(6,1) * t94 + Icges(6,4) * t93 - Icges(6,5) * t190;
t72 = Icges(6,1) * t92 + Icges(6,4) * t91 + Icges(6,5) * t191;
t71 = Icges(6,4) * t94 + Icges(6,2) * t93 - Icges(6,6) * t190;
t70 = Icges(6,4) * t92 + Icges(6,2) * t91 + Icges(6,6) * t191;
t69 = Icges(6,5) * t94 + Icges(6,6) * t93 - Icges(6,3) * t190;
t68 = Icges(6,5) * t92 + Icges(6,6) * t91 + Icges(6,3) * t191;
t67 = V_base(6) * rSges(4,3) - t111 * t147 + t166;
t66 = t109 * t147 + (-rSges(4,3) + t186) * V_base(5) + t180;
t65 = (-t109 + t177) * V_base(6) + (t111 + t176) * V_base(5) + t181;
t64 = t131 * V_base(6) + (-t110 - t86) * t147 + t164;
t63 = t147 * t85 + (-t131 + t186) * V_base(5) + t165;
t62 = (t176 + t86) * V_base(5) + (t167 - t85) * V_base(6) + t179;
t61 = t100 * t112 - t125 * t75 + t132 * V_base(6) + (-t110 + t96) * t147 + t164;
t60 = -t100 * t113 + t125 * t74 + t147 * t95 + (-t132 + t186) * V_base(5) + t165;
t59 = -t112 * t74 + t113 * t75 + (t176 - t96) * V_base(5) + (t167 - t95) * V_base(6) + t179;
t1 = m(1) * (t122 ^ 2 + t123 ^ 2 + t124 ^ 2) / 0.2e1 + m(2) * (t87 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(3) * (t76 ^ 2 + t77 ^ 2 + t78 ^ 2) / 0.2e1 + m(4) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(5) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(6) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + t125 * ((-t69 * t112 - t68 * t113 - t97 * t125) * t158 + ((-t159 * t98 + t161 * t99) * t125 + (-t159 * t70 + t161 * t72) * t113 + (-t159 * t71 + t161 * t73) * t112) * t157) / 0.2e1 + t113 * ((t97 * t191 + t91 * t98 + t92 * t99) * t125 + (t68 * t191 + t91 * t70 + t92 * t72) * t113 + (t69 * t191 + t71 * t91 + t73 * t92) * t112) / 0.2e1 + t112 * ((-t97 * t190 + t93 * t98 + t94 * t99) * t125 + (-t68 * t190 + t70 * t93 + t72 * t94) * t113 + (-t69 * t190 + t93 * t71 + t94 * t73) * t112) / 0.2e1 + ((t157 * t84 + t158 * t82 + t103) * V_base(6) + (t157 * t83 + t158 * t81 + t102) * V_base(5) + (t158 * t129 + t157 * t130 + Icges(4,3)) * t147) * t147 / 0.2e1 + (t163 * t146 + (-t145 * t168 + t102) * t147 + (t105 * t146 + t117 * t149 + t119 * t148 + t136 * t162 + t160 * t138 + t205 * t145 + Icges(1,6)) * V_base(6) + (t104 * t146 + t116 * t149 + t118 * t148 + t135 * t162 + t160 * t137 + t206 * t145 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + (t163 * t145 + (t146 * t168 + t103) * t147 + (t105 * t145 + t117 * t148 - t119 * t149 + t160 * t136 - t138 * t162 - t205 * t146 + Icges(1,3)) * V_base(6) + (t104 * t145 + t116 * t148 - t118 * t149 + t160 * t135 - t137 * t162 - t206 * t146 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + ((-Icges(2,5) * t162 - Icges(3,5) * t149 + Icges(2,6) * t160 + Icges(3,6) * t148) * V_base(6) + (Icges(2,5) * t160 + Icges(3,5) * t148 + Icges(2,6) * t162 + Icges(3,6) * t149) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(3,3) / 0.2e1) * t151) * t151 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4);
T = t1;
