% Calculate kinetic energy for
% S4PRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% rSges [5x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [5x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRPR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRPR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR6_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR6_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:17
% EndTime: 2019-12-31 16:24:18
% DurationCPUTime: 1.60s
% Computational Cost: add. (766->254), mult. (1270->373), div. (0->0), fcn. (1214->8), ass. (0->114)
t158 = cos(pkin(7));
t194 = pkin(3) * t158;
t157 = sin(pkin(6));
t193 = Icges(2,4) * t157;
t161 = sin(qJ(2));
t192 = Icges(3,4) * t161;
t162 = cos(qJ(2));
t191 = Icges(3,4) * t162;
t156 = sin(pkin(7));
t190 = t157 * t156;
t189 = t157 * t161;
t188 = t157 * t162;
t159 = cos(pkin(6));
t187 = t159 * t156;
t186 = t159 * t161;
t185 = t159 * t162;
t173 = pkin(2) * t162 + qJ(3) * t161;
t122 = t173 * t157;
t137 = pkin(1) * t157 - pkin(4) * t159;
t183 = -t122 - t137;
t182 = qJD(3) * t161;
t181 = qJD(4) * t161;
t180 = V_base(5) * qJ(1) + V_base(1);
t176 = qJD(1) + V_base(3);
t145 = qJD(2) * t157 + V_base(4);
t142 = t161 * pkin(2) - qJ(3) * t162;
t144 = -qJD(2) * t159 + V_base(5);
t175 = t144 * t142 + t159 * t182 + t180;
t174 = rSges(3,1) * t162 - rSges(3,2) * t161;
t172 = Icges(3,1) * t162 - t192;
t171 = -Icges(3,2) * t161 + t191;
t170 = Icges(3,5) * t162 - Icges(3,6) * t161;
t138 = pkin(1) * t159 + pkin(4) * t157;
t169 = -V_base(4) * qJ(1) + V_base(6) * t138 + V_base(2);
t168 = V_base(4) * t137 - V_base(5) * t138 + t176;
t167 = pkin(5) * t161 + t162 * t194;
t123 = t173 * t159;
t166 = V_base(6) * t123 + t157 * t182 + t169;
t165 = (Icges(3,5) * t161 + Icges(3,6) * t162) * V_base(6) + (-Icges(3,3) * t159 + t157 * t170) * t144 + (Icges(3,3) * t157 + t159 * t170) * t145;
t164 = -qJD(3) * t162 + t145 * t122 + t168;
t100 = -Icges(3,6) * t159 + t157 * t171;
t101 = Icges(3,6) * t157 + t159 * t171;
t102 = -Icges(3,5) * t159 + t157 * t172;
t103 = Icges(3,5) * t157 + t159 * t172;
t140 = Icges(3,2) * t162 + t192;
t141 = Icges(3,1) * t161 + t191;
t163 = (-t101 * t161 + t103 * t162) * t145 + (-t100 * t161 + t102 * t162) * t144 + (-t140 * t161 + t141 * t162) * V_base(6);
t155 = pkin(7) + qJ(4);
t153 = Icges(2,4) * t159;
t152 = cos(t155);
t151 = sin(t155);
t146 = -qJD(4) * t162 + V_base(6);
t143 = t161 * rSges(3,1) + rSges(3,2) * t162;
t136 = rSges(2,1) * t159 - rSges(2,2) * t157;
t135 = rSges(2,1) * t157 + rSges(2,2) * t159;
t134 = Icges(2,1) * t159 - t193;
t133 = Icges(2,1) * t157 + t153;
t132 = -Icges(2,2) * t157 + t153;
t131 = Icges(2,2) * t159 + t193;
t128 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t127 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t126 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t121 = t159 * t181 + t145;
t120 = t157 * t181 + t144;
t119 = t158 * t185 + t190;
t118 = -t156 * t185 + t157 * t158;
t117 = t158 * t188 - t187;
t116 = -t156 * t188 - t159 * t158;
t113 = -rSges(4,3) * t162 + (rSges(4,1) * t158 - rSges(4,2) * t156) * t161;
t112 = t157 * rSges(3,3) + t159 * t174;
t111 = -t159 * rSges(3,3) + t157 * t174;
t110 = -Icges(4,5) * t162 + (Icges(4,1) * t158 - Icges(4,4) * t156) * t161;
t109 = -Icges(4,6) * t162 + (Icges(4,4) * t158 - Icges(4,2) * t156) * t161;
t108 = -Icges(4,3) * t162 + (Icges(4,5) * t158 - Icges(4,6) * t156) * t161;
t107 = t157 * t151 + t152 * t185;
t106 = -t151 * t185 + t157 * t152;
t105 = -t159 * t151 + t152 * t188;
t104 = -t151 * t188 - t159 * t152;
t97 = -rSges(5,3) * t162 + (rSges(5,1) * t152 - rSges(5,2) * t151) * t161;
t96 = -Icges(5,5) * t162 + (Icges(5,1) * t152 - Icges(5,4) * t151) * t161;
t95 = -Icges(5,6) * t162 + (Icges(5,4) * t152 - Icges(5,2) * t151) * t161;
t94 = -Icges(5,3) * t162 + (Icges(5,5) * t152 - Icges(5,6) * t151) * t161;
t92 = V_base(5) * rSges(2,3) - t135 * V_base(6) + t180;
t91 = t136 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t90 = -pkin(5) * t162 + t161 * t194;
t89 = t135 * V_base(4) - t136 * V_base(5) + t176;
t88 = pkin(3) * t190 + t159 * t167;
t87 = -pkin(3) * t187 + t157 * t167;
t86 = rSges(4,1) * t119 + rSges(4,2) * t118 + rSges(4,3) * t186;
t85 = rSges(4,1) * t117 + rSges(4,2) * t116 + rSges(4,3) * t189;
t84 = Icges(4,1) * t119 + Icges(4,4) * t118 + Icges(4,5) * t186;
t83 = Icges(4,1) * t117 + Icges(4,4) * t116 + Icges(4,5) * t189;
t82 = Icges(4,4) * t119 + Icges(4,2) * t118 + Icges(4,6) * t186;
t81 = Icges(4,4) * t117 + Icges(4,2) * t116 + Icges(4,6) * t189;
t80 = Icges(4,5) * t119 + Icges(4,6) * t118 + Icges(4,3) * t186;
t79 = Icges(4,5) * t117 + Icges(4,6) * t116 + Icges(4,3) * t189;
t78 = rSges(5,1) * t107 + rSges(5,2) * t106 + rSges(5,3) * t186;
t77 = rSges(5,1) * t105 + rSges(5,2) * t104 + rSges(5,3) * t189;
t76 = Icges(5,1) * t107 + Icges(5,4) * t106 + Icges(5,5) * t186;
t75 = Icges(5,1) * t105 + Icges(5,4) * t104 + Icges(5,5) * t189;
t74 = Icges(5,4) * t107 + Icges(5,2) * t106 + Icges(5,6) * t186;
t73 = Icges(5,4) * t105 + Icges(5,2) * t104 + Icges(5,6) * t189;
t72 = Icges(5,5) * t107 + Icges(5,6) * t106 + Icges(5,3) * t186;
t71 = Icges(5,5) * t105 + Icges(5,6) * t104 + Icges(5,3) * t189;
t70 = t143 * t144 + (-t111 - t137) * V_base(6) + t180;
t69 = t112 * V_base(6) - t143 * t145 + t169;
t68 = t111 * t145 - t112 * t144 + t168;
t67 = t113 * t144 + (-t85 + t183) * V_base(6) + t175;
t66 = t86 * V_base(6) + (-t113 - t142) * t145 + t166;
t65 = t145 * t85 + (-t123 - t86) * t144 + t164;
t64 = t120 * t97 + t144 * t90 - t146 * t77 + (-t87 + t183) * V_base(6) + t175;
t63 = -t121 * t97 + t146 * t78 + t88 * V_base(6) + (-t142 - t90) * t145 + t166;
t62 = -t120 * t78 + t121 * t77 + t145 * t87 + (-t123 - t88) * t144 + t164;
t1 = m(1) * (t126 ^ 2 + t127 ^ 2 + t128 ^ 2) / 0.2e1 + m(2) * (t89 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(3) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(4) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(5) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + t121 * ((t106 * t74 + t107 * t76 + t72 * t186) * t121 + (t106 * t73 + t107 * t75 + t186 * t71) * t120 + (t106 * t95 + t107 * t96 + t186 * t94) * t146) / 0.2e1 + t120 * ((t104 * t74 + t105 * t76 + t189 * t72) * t121 + (t104 * t73 + t105 * t75 + t71 * t189) * t120 + (t104 * t95 + t105 * t96 + t189 * t94) * t146) / 0.2e1 + t146 * ((-t71 * t120 - t72 * t121 - t94 * t146) * t162 + ((-t151 * t74 + t152 * t76) * t121 + (-t151 * t73 + t152 * t75) * t120 + (-t151 * t95 + t152 * t96) * t146) * t161) / 0.2e1 + (t163 * t157 - t165 * t159 + (t116 * t82 + t117 * t84 + t189 * t80) * t145 + (t116 * t81 + t117 * t83 + t189 * t79) * t144 + (t108 * t189 + t109 * t116 + t110 * t117) * V_base(6)) * t144 / 0.2e1 + (t165 * t157 + t163 * t159 + (t118 * t82 + t119 * t84 + t186 * t80) * t145 + (t118 * t81 + t119 * t83 + t186 * t79) * t144 + (t108 * t186 + t109 * t118 + t110 * t119) * V_base(6)) * t145 / 0.2e1 + ((-t131 * t157 + t133 * t159 + Icges(1,4)) * V_base(5) + (-t132 * t157 + t134 * t159 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t131 * t159 + t133 * t157 + Icges(1,2)) * V_base(5) + (t132 * t159 + t134 * t157 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t101 * t162 + t161 * t103) * t145 + (t100 * t162 + t161 * t102) * t144 + (-t79 * t144 - t80 * t145) * t162 + ((-t156 * t82 + t158 * t84) * t145 + (-t156 * t81 + t158 * t83) * t144) * t161 + (Icges(1,3) + Icges(2,3) + (t140 - t108) * t162 + (-t109 * t156 + t110 * t158 + t141) * t161) * V_base(6)) * V_base(6) / 0.2e1 + V_base(6) * V_base(4) * (Icges(2,5) * t159 - Icges(2,6) * t157 + Icges(1,5)) + V_base(6) * V_base(5) * (Icges(2,5) * t157 + Icges(2,6) * t159 + Icges(1,6));
T = t1;
