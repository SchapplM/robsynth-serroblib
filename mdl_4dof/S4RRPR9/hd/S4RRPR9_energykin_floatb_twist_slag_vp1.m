% Calculate kinetic energy for
% S4RRPR9
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR9_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR9_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRPR9_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR9_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR9_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR9_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:09:18
% EndTime: 2019-12-31 17:09:20
% DurationCPUTime: 1.73s
% Computational Cost: add. (798->254), mult. (1270->378), div. (0->0), fcn. (1214->8), ass. (0->114)
t155 = cos(pkin(7));
t190 = pkin(3) * t155;
t158 = sin(qJ(1));
t189 = Icges(2,4) * t158;
t157 = sin(qJ(2));
t188 = Icges(3,4) * t157;
t159 = cos(qJ(2));
t187 = Icges(3,4) * t159;
t154 = sin(pkin(7));
t160 = cos(qJ(1));
t186 = t154 * t160;
t185 = t157 * t158;
t184 = t157 * t160;
t183 = t158 * t154;
t182 = t158 * t159;
t181 = t159 * t160;
t171 = pkin(2) * t159 + qJ(3) * t157;
t119 = t171 * t158;
t140 = t158 * pkin(1) - pkin(5) * t160;
t179 = -t119 - t140;
t178 = qJD(3) * t157;
t177 = qJD(4) * t157;
t176 = V_base(5) * pkin(4) + V_base(1);
t143 = qJD(2) * t158 + V_base(4);
t149 = V_base(6) + qJD(1);
t136 = pkin(2) * t157 - qJ(3) * t159;
t142 = -qJD(2) * t160 + V_base(5);
t173 = t142 * t136 + t160 * t178 + t176;
t172 = rSges(3,1) * t159 - rSges(3,2) * t157;
t170 = Icges(3,1) * t159 - t188;
t169 = -Icges(3,2) * t157 + t187;
t168 = Icges(3,5) * t159 - Icges(3,6) * t157;
t141 = pkin(1) * t160 + t158 * pkin(5);
t167 = -V_base(4) * pkin(4) + t149 * t141 + V_base(2);
t166 = V_base(4) * t140 - t141 * V_base(5) + V_base(3);
t165 = (-Icges(3,3) * t160 + t158 * t168) * t142 + (Icges(3,3) * t158 + t160 * t168) * t143 + (Icges(3,5) * t157 + Icges(3,6) * t159) * t149;
t120 = t171 * t160;
t164 = t149 * t120 + t158 * t178 + t167;
t163 = pkin(6) * t157 + t159 * t190;
t162 = -qJD(3) * t159 + t143 * t119 + t166;
t106 = -Icges(3,6) * t160 + t158 * t169;
t107 = Icges(3,6) * t158 + t160 * t169;
t108 = -Icges(3,5) * t160 + t158 * t170;
t109 = Icges(3,5) * t158 + t160 * t170;
t129 = Icges(3,2) * t159 + t188;
t132 = Icges(3,1) * t157 + t187;
t161 = (-t107 * t157 + t109 * t159) * t143 + (-t106 * t157 + t108 * t159) * t142 + (-t129 * t157 + t132 * t159) * t149;
t153 = pkin(7) + qJ(4);
t151 = Icges(2,4) * t160;
t148 = cos(t153);
t147 = sin(t153);
t139 = rSges(2,1) * t160 - t158 * rSges(2,2);
t138 = t158 * rSges(2,1) + rSges(2,2) * t160;
t137 = rSges(3,1) * t157 + rSges(3,2) * t159;
t135 = -qJD(4) * t159 + t149;
t134 = Icges(2,1) * t160 - t189;
t133 = Icges(2,1) * t158 + t151;
t131 = -Icges(2,2) * t158 + t151;
t130 = Icges(2,2) * t160 + t189;
t125 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t124 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t123 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t118 = t160 * t177 + t143;
t117 = t158 * t177 + t142;
t116 = t155 * t181 + t183;
t115 = -t154 * t181 + t158 * t155;
t114 = t155 * t182 - t186;
t113 = -t154 * t182 - t155 * t160;
t111 = t158 * rSges(3,3) + t160 * t172;
t110 = -rSges(3,3) * t160 + t158 * t172;
t103 = t158 * t147 + t148 * t181;
t102 = -t147 * t181 + t158 * t148;
t101 = -t147 * t160 + t148 * t182;
t100 = -t147 * t182 - t148 * t160;
t99 = -rSges(4,3) * t159 + (rSges(4,1) * t155 - rSges(4,2) * t154) * t157;
t97 = -Icges(4,5) * t159 + (Icges(4,1) * t155 - Icges(4,4) * t154) * t157;
t96 = -Icges(4,6) * t159 + (Icges(4,4) * t155 - Icges(4,2) * t154) * t157;
t95 = -Icges(4,3) * t159 + (Icges(4,5) * t155 - Icges(4,6) * t154) * t157;
t93 = -rSges(5,3) * t159 + (rSges(5,1) * t148 - rSges(5,2) * t147) * t157;
t92 = -Icges(5,5) * t159 + (Icges(5,1) * t148 - Icges(5,4) * t147) * t157;
t91 = -Icges(5,6) * t159 + (Icges(5,4) * t148 - Icges(5,2) * t147) * t157;
t90 = -Icges(5,3) * t159 + (Icges(5,5) * t148 - Icges(5,6) * t147) * t157;
t89 = -pkin(6) * t159 + t157 * t190;
t88 = V_base(5) * rSges(2,3) - t138 * t149 + t176;
t87 = t139 * t149 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t86 = t138 * V_base(4) - t139 * V_base(5) + V_base(3);
t85 = pkin(3) * t183 + t160 * t163;
t84 = -pkin(3) * t186 + t158 * t163;
t83 = t116 * rSges(4,1) + t115 * rSges(4,2) + rSges(4,3) * t184;
t82 = rSges(4,1) * t114 + rSges(4,2) * t113 + rSges(4,3) * t185;
t81 = Icges(4,1) * t116 + Icges(4,4) * t115 + Icges(4,5) * t184;
t80 = Icges(4,1) * t114 + Icges(4,4) * t113 + Icges(4,5) * t185;
t79 = Icges(4,4) * t116 + Icges(4,2) * t115 + Icges(4,6) * t184;
t78 = Icges(4,4) * t114 + Icges(4,2) * t113 + Icges(4,6) * t185;
t77 = Icges(4,5) * t116 + Icges(4,6) * t115 + Icges(4,3) * t184;
t76 = Icges(4,5) * t114 + Icges(4,6) * t113 + Icges(4,3) * t185;
t75 = t103 * rSges(5,1) + t102 * rSges(5,2) + rSges(5,3) * t184;
t74 = rSges(5,1) * t101 + rSges(5,2) * t100 + rSges(5,3) * t185;
t73 = Icges(5,1) * t103 + Icges(5,4) * t102 + Icges(5,5) * t184;
t72 = Icges(5,1) * t101 + Icges(5,4) * t100 + Icges(5,5) * t185;
t71 = Icges(5,4) * t103 + Icges(5,2) * t102 + Icges(5,6) * t184;
t70 = Icges(5,4) * t101 + Icges(5,2) * t100 + Icges(5,6) * t185;
t69 = Icges(5,5) * t103 + Icges(5,6) * t102 + Icges(5,3) * t184;
t68 = Icges(5,5) * t101 + Icges(5,6) * t100 + Icges(5,3) * t185;
t67 = t137 * t142 + (-t110 - t140) * t149 + t176;
t66 = t111 * t149 - t137 * t143 + t167;
t65 = t110 * t143 - t111 * t142 + t166;
t64 = t142 * t99 + (-t82 + t179) * t149 + t173;
t63 = t149 * t83 + (-t136 - t99) * t143 + t164;
t62 = t143 * t82 + (-t120 - t83) * t142 + t162;
t61 = t117 * t93 - t135 * t74 + t142 * t89 + (-t84 + t179) * t149 + t173;
t60 = -t118 * t93 + t135 * t75 + t149 * t85 + (-t136 - t89) * t143 + t164;
t59 = -t117 * t75 + t118 * t74 + t143 * t84 + (-t120 - t85) * t142 + t162;
t1 = m(1) * (t123 ^ 2 + t124 ^ 2 + t125 ^ 2) / 0.2e1 + m(2) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(3) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(4) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + m(5) * (t59 ^ 2 + t60 ^ 2 + t61 ^ 2) / 0.2e1 + t118 * ((t102 * t71 + t103 * t73 + t69 * t184) * t118 + (t102 * t70 + t103 * t72 + t184 * t68) * t117 + (t102 * t91 + t103 * t92 + t184 * t90) * t135) / 0.2e1 + t117 * ((t100 * t71 + t101 * t73 + t185 * t69) * t118 + (t100 * t70 + t101 * t72 + t68 * t185) * t117 + (t100 * t91 + t101 * t92 + t185 * t90) * t135) / 0.2e1 + t135 * ((-t68 * t117 - t69 * t118 - t90 * t135) * t159 + ((-t147 * t71 + t148 * t73) * t118 + (-t147 * t70 + t148 * t72) * t117 + (-t147 * t91 + t148 * t92) * t135) * t157) / 0.2e1 + (t161 * t158 - t165 * t160 + (t113 * t79 + t114 * t81 + t185 * t77) * t143 + (t113 * t78 + t114 * t80 + t185 * t76) * t142 + (t113 * t96 + t114 * t97 + t185 * t95) * t149) * t142 / 0.2e1 + (t165 * t158 + t161 * t160 + (t115 * t79 + t116 * t81 + t184 * t77) * t143 + (t115 * t78 + t116 * t80 + t184 * t76) * t142 + (t115 * t96 + t116 * t97 + t184 * t95) * t149) * t143 / 0.2e1 + ((-t158 * t130 + t133 * t160 + Icges(1,4)) * V_base(5) + (-t158 * t131 + t134 * t160 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t130 * t160 + t158 * t133 + Icges(1,2)) * V_base(5) + (t131 * t160 + t158 * t134 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t107 * t159 + t109 * t157) * t143 + (t106 * t159 + t108 * t157) * t142 + (-t76 * t142 - t77 * t143) * t159 + ((-t154 * t79 + t155 * t81) * t143 + (-t154 * t78 + t155 * t80) * t142) * t157 + (Icges(2,3) + (t129 - t95) * t159 + (-t154 * t96 + t155 * t97 + t132) * t157) * t149) * t149 / 0.2e1 + V_base(4) * t149 * (Icges(2,5) * t160 - Icges(2,6) * t158) + V_base(5) * t149 * (Icges(2,5) * t158 + Icges(2,6) * t160) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
