% Calculate kinetic energy for
% S4RRPR6
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
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR6_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRPR6_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR6_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR6_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:20
% EndTime: 2019-12-31 17:04:22
% DurationCPUTime: 1.14s
% Computational Cost: add. (770->204), mult. (838->285), div. (0->0), fcn. (670->8), ass. (0->105)
t189 = Icges(3,3) + Icges(4,3);
t137 = qJ(2) + pkin(7);
t127 = sin(t137);
t128 = cos(t137);
t139 = sin(qJ(2));
t141 = cos(qJ(2));
t188 = Icges(3,5) * t141 + Icges(4,5) * t128 - Icges(3,6) * t139 - Icges(4,6) * t127;
t176 = Icges(3,4) * t139;
t111 = Icges(3,2) * t141 + t176;
t175 = Icges(3,4) * t141;
t114 = Icges(3,1) * t139 + t175;
t142 = cos(qJ(1));
t122 = -qJD(2) * t142 + V_base(5);
t140 = sin(qJ(1));
t123 = qJD(2) * t140 + V_base(4);
t130 = V_base(6) + qJD(1);
t173 = Icges(4,4) * t128;
t157 = -Icges(4,2) * t127 + t173;
t77 = -Icges(4,6) * t142 + t157 * t140;
t78 = Icges(4,6) * t140 + t157 * t142;
t174 = Icges(4,4) * t127;
t160 = Icges(4,1) * t128 - t174;
t79 = -Icges(4,5) * t142 + t160 * t140;
t80 = Icges(4,5) * t140 + t160 * t142;
t158 = -Icges(3,2) * t139 + t175;
t85 = -Icges(3,6) * t142 + t158 * t140;
t86 = Icges(3,6) * t140 + t158 * t142;
t161 = Icges(3,1) * t141 - t176;
t87 = -Icges(3,5) * t142 + t161 * t140;
t88 = Icges(3,5) * t140 + t161 * t142;
t97 = Icges(4,2) * t128 + t174;
t98 = Icges(4,1) * t127 + t173;
t187 = (-t111 * t139 + t114 * t141 - t127 * t97 + t128 * t98) * t130 + (-t127 * t78 + t128 * t80 - t139 * t86 + t141 * t88) * t123 + (-t127 * t77 + t128 * t79 - t139 * t85 + t141 * t87) * t122;
t186 = (Icges(3,5) * t139 + Icges(4,5) * t127 + Icges(3,6) * t141 + Icges(4,6) * t128) * t130 + (t189 * t140 + t188 * t142) * t123 + (t188 * t140 - t189 * t142) * t122;
t182 = pkin(2) * t139;
t181 = pkin(3) * t127;
t180 = t141 * pkin(2);
t120 = t140 * pkin(1) - pkin(5) * t142;
t73 = -qJ(3) * t142 + t180 * t140;
t178 = -t120 - t73;
t177 = Icges(2,4) * t140;
t129 = qJ(4) + t137;
t124 = sin(t129);
t172 = Icges(5,4) * t124;
t125 = cos(t129);
t171 = Icges(5,4) * t125;
t170 = pkin(3) * t128;
t168 = V_base(5) * pkin(4) + V_base(1);
t165 = qJD(3) * t140 + t122 * t182 + t168;
t164 = rSges(3,1) * t141 - rSges(3,2) * t139;
t163 = rSges(4,1) * t128 - rSges(4,2) * t127;
t162 = rSges(5,1) * t125 - rSges(5,2) * t124;
t159 = Icges(5,1) * t125 - t172;
t156 = -Icges(5,2) * t124 + t171;
t153 = Icges(5,5) * t125 - Icges(5,6) * t124;
t121 = pkin(1) * t142 + t140 * pkin(5);
t152 = -V_base(4) * pkin(4) + t130 * t121 + V_base(2);
t151 = V_base(4) * t120 - t121 * V_base(5) + V_base(3);
t150 = t123 * t73 + t151;
t100 = V_base(5) + (-qJD(2) - qJD(4)) * t142;
t101 = qJD(4) * t140 + t123;
t149 = (-Icges(5,3) * t142 + t153 * t140) * t100 + (Icges(5,3) * t140 + t153 * t142) * t101 + (Icges(5,5) * t124 + Icges(5,6) * t125) * t130;
t74 = qJ(3) * t140 + t180 * t142;
t146 = -qJD(3) * t142 + t130 * t74 + t152;
t65 = -Icges(5,6) * t142 + t156 * t140;
t66 = Icges(5,6) * t140 + t156 * t142;
t67 = -Icges(5,5) * t142 + t159 * t140;
t68 = Icges(5,5) * t140 + t159 * t142;
t92 = Icges(5,2) * t125 + t172;
t93 = Icges(5,1) * t124 + t171;
t145 = (-t124 * t66 + t125 * t68) * t101 + (-t124 * t65 + t125 * t67) * t100 + (-t124 * t92 + t125 * t93) * t130;
t133 = Icges(2,4) * t142;
t119 = rSges(2,1) * t142 - t140 * rSges(2,2);
t118 = t140 * rSges(2,1) + rSges(2,2) * t142;
t117 = rSges(3,1) * t139 + rSges(3,2) * t141;
t116 = Icges(2,1) * t142 - t177;
t115 = Icges(2,1) * t140 + t133;
t113 = -Icges(2,2) * t140 + t133;
t112 = Icges(2,2) * t142 + t177;
t107 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t106 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t105 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t99 = rSges(4,1) * t127 + rSges(4,2) * t128;
t94 = rSges(5,1) * t124 + rSges(5,2) * t125;
t90 = t140 * rSges(3,3) + t164 * t142;
t89 = -rSges(3,3) * t142 + t164 * t140;
t82 = t140 * rSges(4,3) + t163 * t142;
t81 = -rSges(4,3) * t142 + t163 * t140;
t72 = t140 * rSges(5,3) + t162 * t142;
t71 = -rSges(5,3) * t142 + t162 * t140;
t70 = V_base(5) * rSges(2,3) - t118 * t130 + t168;
t69 = t119 * t130 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t62 = t118 * V_base(4) - t119 * V_base(5) + V_base(3);
t59 = pkin(6) * t140 + t170 * t142;
t58 = -pkin(6) * t142 + t170 * t140;
t57 = t117 * t122 + (-t120 - t89) * t130 + t168;
t56 = -t117 * t123 + t130 * t90 + t152;
t55 = -t122 * t90 + t123 * t89 + t151;
t54 = t122 * t99 + (-t81 + t178) * t130 + t165;
t53 = t130 * t82 + (-t99 - t182) * t123 + t146;
t52 = t123 * t81 + (-t74 - t82) * t122 + t150;
t51 = t122 * t181 + t100 * t94 + (-t58 - t71 + t178) * t130 + t165;
t50 = -t101 * t94 + (t59 + t72) * t130 + (-t181 - t182) * t123 + t146;
t49 = -t100 * t72 + t101 * t71 + t123 * t58 + (-t59 - t74) * t122 + t150;
t1 = m(1) * (t105 ^ 2 + t106 ^ 2 + t107 ^ 2) / 0.2e1 + m(2) * (t62 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + m(3) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + m(4) * (t52 ^ 2 + t53 ^ 2 + t54 ^ 2) / 0.2e1 + m(5) * (t49 ^ 2 + t50 ^ 2 + t51 ^ 2) / 0.2e1 + t101 * (t149 * t140 + t145 * t142) / 0.2e1 + t100 * (t145 * t140 - t149 * t142) / 0.2e1 + (t187 * t140 - t186 * t142) * t122 / 0.2e1 + (t186 * t140 + t187 * t142) * t123 / 0.2e1 + ((-t140 * t112 + t115 * t142 + Icges(1,4)) * V_base(5) + (-t140 * t113 + t116 * t142 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t112 * t142 + t140 * t115 + Icges(1,2)) * V_base(5) + (t113 * t142 + t140 * t116 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + ((t124 * t68 + t125 * t66) * t101 + (t124 * t67 + t125 * t65) * t100 + (t127 * t80 + t128 * t78 + t139 * t88 + t141 * t86) * t123 + (t127 * t79 + t128 * t77 + t139 * t87 + t141 * t85) * t122 + (t111 * t141 + t114 * t139 + t124 * t93 + t125 * t92 + t127 * t98 + t128 * t97 + Icges(2,3)) * t130) * t130 / 0.2e1 + t130 * V_base(4) * (Icges(2,5) * t142 - Icges(2,6) * t140) + V_base(5) * t130 * (Icges(2,5) * t140 + Icges(2,6) * t142) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
