% Calculate kinetic energy for
% S5RRPRR3
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
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR3_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRPRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_energykin_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR3_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR3_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRPRR3_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:00:19
% EndTime: 2020-01-03 12:00:20
% DurationCPUTime: 1.42s
% Computational Cost: add. (1001->214), mult. (596->272), div. (0->0), fcn. (372->10), ass. (0->108)
t146 = qJ(1) + qJ(2);
t137 = pkin(9) + t146;
t136 = qJ(4) + t137;
t130 = sin(t136);
t106 = -qJD(5) * t130 + V_base(6);
t131 = cos(t136);
t107 = -qJD(5) * t131 + V_base(5);
t149 = cos(qJ(5));
t147 = sin(qJ(5));
t176 = Icges(6,4) * t147;
t116 = Icges(6,2) * t149 + t176;
t175 = Icges(6,4) * t149;
t119 = Icges(6,1) * t147 + t175;
t138 = V_base(4) + qJD(1);
t135 = qJD(2) + t138;
t129 = qJD(4) + t135;
t158 = -Icges(6,2) * t147 + t175;
t69 = -Icges(6,6) * t131 + t130 * t158;
t70 = -Icges(6,6) * t130 - t131 * t158;
t159 = Icges(6,1) * t149 - t176;
t71 = -Icges(6,5) * t131 + t130 * t159;
t72 = -Icges(6,5) * t130 - t131 * t159;
t189 = (t116 * t147 - t119 * t149) * t129 + (t147 * t69 - t149 * t71) * t107 + (t147 * t70 - t149 * t72) * t106;
t188 = -pkin(5) - pkin(6);
t148 = sin(qJ(1));
t185 = pkin(1) * t148;
t150 = cos(qJ(1));
t184 = pkin(1) * t150;
t139 = sin(t146);
t183 = pkin(2) * t139;
t140 = cos(t146);
t182 = pkin(2) * t140;
t181 = pkin(3) * t135;
t180 = Icges(2,4) * t150;
t179 = Icges(3,4) * t140;
t134 = cos(t137);
t178 = Icges(4,4) * t134;
t177 = Icges(5,4) * t131;
t174 = -qJ(3) + t188;
t173 = t138 * t185 + V_base(3);
t172 = V_base(6) * pkin(5) + V_base(2);
t169 = qJD(3) + V_base(1);
t168 = -pkin(7) + t174;
t167 = t135 * t183 + t173;
t133 = sin(t137);
t166 = t133 * t181 + t167;
t165 = V_base(6) * pkin(6) + t138 * t184 + t172;
t164 = -t183 - t185;
t163 = -t182 - t184;
t162 = rSges(6,1) * t149 - rSges(6,2) * t147;
t157 = Icges(6,5) * t149 - Icges(6,6) * t147;
t155 = V_base(6) * qJ(3) + t135 * t182 + t165;
t154 = -pkin(3) * t133 + t164;
t153 = -pkin(3) * t134 + t163;
t152 = -t106 * (-Icges(6,3) * t130 - t131 * t157) - t107 * (-Icges(6,3) * t131 + t130 * t157) - (Icges(6,5) * t147 + Icges(6,6) * t149) * t129;
t151 = V_base(6) * pkin(7) + t134 * t181 + t155;
t142 = Icges(2,4) * t148;
t132 = Icges(3,4) * t139;
t128 = Icges(4,4) * t133;
t125 = Icges(5,4) * t130;
t124 = -rSges(2,1) * t150 + t148 * rSges(2,2);
t123 = t148 * rSges(2,1) + rSges(2,2) * t150;
t122 = rSges(6,1) * t147 + rSges(6,2) * t149;
t121 = -Icges(2,1) * t150 + t142;
t120 = Icges(2,1) * t148 + t180;
t118 = Icges(2,2) * t148 - t180;
t117 = Icges(2,2) * t150 + t142;
t112 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t111 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t110 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t105 = -rSges(3,1) * t140 + rSges(3,2) * t139;
t104 = rSges(3,1) * t139 + rSges(3,2) * t140;
t103 = -Icges(3,1) * t140 + t132;
t102 = Icges(3,1) * t139 + t179;
t101 = Icges(3,2) * t139 - t179;
t100 = Icges(3,2) * t140 + t132;
t95 = -rSges(4,1) * t134 + rSges(4,2) * t133;
t94 = rSges(4,1) * t133 + rSges(4,2) * t134;
t93 = -Icges(4,1) * t134 + t128;
t92 = Icges(4,1) * t133 + t178;
t91 = Icges(4,2) * t133 - t178;
t90 = Icges(4,2) * t134 + t128;
t87 = -pkin(4) * t131 - pkin(8) * t130;
t86 = pkin(4) * t130 - pkin(8) * t131;
t85 = -rSges(5,1) * t131 + rSges(5,2) * t130;
t84 = rSges(5,1) * t130 + rSges(5,2) * t131;
t83 = -Icges(5,1) * t131 + t125;
t82 = Icges(5,1) * t130 + t177;
t81 = Icges(5,2) * t130 - t177;
t80 = Icges(5,2) * t131 + t125;
t77 = V_base(6) * rSges(2,3) - t124 * t138 + t172;
t76 = t123 * t138 + V_base(3) + (-rSges(2,3) - pkin(5)) * V_base(5);
t75 = -t123 * V_base(6) + t124 * V_base(5) + V_base(1);
t74 = -rSges(6,3) * t130 - t131 * t162;
t73 = -rSges(6,3) * t131 + t130 * t162;
t66 = V_base(6) * rSges(3,3) - t105 * t135 + t165;
t65 = t104 * t135 + (-rSges(3,3) + t188) * V_base(5) + t173;
t64 = -V_base(6) * t104 + V_base(5) * t105 + V_base(1) + (-t148 * V_base(6) - t150 * V_base(5)) * pkin(1);
t63 = V_base(6) * rSges(4,3) - t135 * t95 + t155;
t62 = t135 * t94 + (-rSges(4,3) + t174) * V_base(5) + t167;
t61 = (t164 - t94) * V_base(6) + (t163 + t95) * V_base(5) + t169;
t60 = V_base(6) * rSges(5,3) - t129 * t85 + t151;
t59 = t129 * t84 + (-rSges(5,3) + t168) * V_base(5) + t166;
t58 = (t154 - t84) * V_base(6) + (t153 + t85) * V_base(5) + t169;
t57 = t106 * t122 + (-t74 - t87) * t129 + t151;
t56 = -t107 * t122 + (t73 + t86) * t129 + t168 * V_base(5) + t166;
t55 = -t106 * t73 + t107 * t74 + (t154 - t86) * V_base(6) + (t153 + t87) * V_base(5) + t169;
t1 = m(1) * (t110 ^ 2 + t111 ^ 2 + t112 ^ 2) / 0.2e1 + m(2) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(3) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + m(4) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + m(5) * (t58 ^ 2 + t59 ^ 2 + t60 ^ 2) / 0.2e1 + m(6) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + t107 * (-t189 * t130 + t152 * t131) / 0.2e1 + t106 * (t152 * t130 + t189 * t131) / 0.2e1 + ((t147 * t71 + t149 * t69) * t107 + (t147 * t72 + t149 * t70) * t106 + (t149 * t116 + t147 * t119 + Icges(5,3)) * t129) * t129 / 0.2e1 + ((t101 * t140 + t103 * t139 + t118 * t150 + t148 * t121 + t130 * t83 + t131 * t81 + t133 * t93 + t134 * t91 + Icges(1,6)) * V_base(6) + (t140 * t100 + t139 * t102 + t150 * t117 + t148 * t120 + t130 * t82 + t131 * t80 + t133 * t92 + t134 * t90 + Icges(1,2)) * V_base(5)) * V_base(5) / 0.2e1 + ((t139 * t101 - t140 * t103 + t148 * t118 - t150 * t121 + t130 * t81 - t131 * t83 + t133 * t91 - t134 * t93 + Icges(1,3)) * V_base(6) + (t100 * t139 - t102 * t140 + t148 * t117 - t120 * t150 + t130 * t80 - t131 * t82 + t133 * t90 - t134 * t92 + Icges(1,6)) * V_base(5)) * V_base(6) / 0.2e1 + t129 * V_base(6) * (-Icges(5,5) * t131 + Icges(5,6) * t130) + t129 * V_base(5) * (Icges(5,5) * t130 + Icges(5,6) * t131) + ((-Icges(3,5) * t140 - Icges(4,5) * t134 + Icges(3,6) * t139 + Icges(4,6) * t133) * V_base(6) + (Icges(3,5) * t139 + Icges(4,5) * t133 + Icges(3,6) * t140 + Icges(4,6) * t134) * V_base(5) + (Icges(3,3) / 0.2e1 + Icges(4,3) / 0.2e1) * t135) * t135 + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + Icges(1,1) * V_base(4) / 0.2e1) * V_base(4) + ((Icges(2,5) * t148 + Icges(2,6) * t150) * V_base(5) + (-Icges(2,5) * t150 + Icges(2,6) * t148) * V_base(6) + Icges(2,3) * t138 / 0.2e1) * t138;
T = t1;
