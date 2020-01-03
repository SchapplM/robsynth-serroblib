% Calculate kinetic energy for
% S4RPRR4
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR4_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR4_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:12
% EndTime: 2019-12-31 16:50:13
% DurationCPUTime: 1.24s
% Computational Cost: add. (795->212), mult. (862->305), div. (0->0), fcn. (748->8), ass. (0->101)
t136 = sin(qJ(1));
t167 = pkin(1) * t136;
t139 = cos(qJ(1));
t166 = pkin(1) * t139;
t165 = -pkin(4) - qJ(2);
t164 = Icges(2,4) * t136;
t133 = qJ(1) + pkin(7);
t127 = sin(t133);
t163 = Icges(3,4) * t127;
t135 = sin(qJ(3));
t162 = Icges(4,4) * t135;
t138 = cos(qJ(3));
t161 = Icges(4,4) * t138;
t160 = t127 * t135;
t128 = cos(t133);
t159 = t128 * t135;
t134 = sin(qJ(4));
t158 = t134 * t138;
t137 = cos(qJ(4));
t157 = t137 * t138;
t156 = qJD(4) * t135;
t129 = V_base(6) + qJD(1);
t155 = t129 * t166 + V_base(2);
t154 = V_base(5) * pkin(4) + V_base(1);
t108 = qJD(3) * t127 + V_base(4);
t102 = pkin(2) * t127 - pkin(5) * t128;
t151 = -t102 - t167;
t150 = V_base(5) * qJ(2) + t154;
t149 = V_base(4) * t167 + qJD(2) + V_base(3);
t148 = pkin(3) * t138 + pkin(6) * t135;
t107 = -qJD(3) * t128 + V_base(5);
t147 = rSges(4,1) * t138 - rSges(4,2) * t135;
t146 = Icges(4,1) * t138 - t162;
t145 = -Icges(4,2) * t135 + t161;
t144 = Icges(4,5) * t138 - Icges(4,6) * t135;
t143 = (-Icges(4,3) * t128 + t127 * t144) * t107 + (Icges(4,3) * t127 + t128 * t144) * t108 + (Icges(4,5) * t135 + Icges(4,6) * t138) * t129;
t103 = pkin(2) * t128 + pkin(5) * t127;
t142 = t129 * t103 + t165 * V_base(4) + t155;
t141 = V_base(4) * t102 + (-t103 - t166) * V_base(5) + t149;
t112 = Icges(4,2) * t138 + t162;
t115 = Icges(4,1) * t135 + t161;
t72 = -Icges(4,6) * t128 + t127 * t145;
t73 = Icges(4,6) * t127 + t128 * t145;
t74 = -Icges(4,5) * t128 + t127 * t146;
t75 = Icges(4,5) * t127 + t128 * t146;
t140 = (-t135 * t73 + t138 * t75) * t108 + (-t135 * t72 + t138 * t74) * t107 + (-t112 * t135 + t115 * t138) * t129;
t131 = Icges(2,4) * t139;
t126 = Icges(3,4) * t128;
t122 = pkin(3) * t135 - pkin(6) * t138;
t121 = rSges(2,1) * t139 - t136 * rSges(2,2);
t120 = t136 * rSges(2,1) + rSges(2,2) * t139;
t119 = rSges(4,1) * t135 + rSges(4,2) * t138;
t118 = -qJD(4) * t138 + t129;
t117 = Icges(2,1) * t139 - t164;
t116 = Icges(2,1) * t136 + t131;
t114 = -Icges(2,2) * t136 + t131;
t113 = Icges(2,2) * t139 + t164;
t106 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t105 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t104 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t101 = rSges(3,1) * t128 - rSges(3,2) * t127;
t100 = rSges(3,1) * t127 + rSges(3,2) * t128;
t99 = Icges(3,1) * t128 - t163;
t98 = Icges(3,1) * t127 + t126;
t97 = -Icges(3,2) * t127 + t126;
t96 = Icges(3,2) * t128 + t163;
t92 = t148 * t128;
t91 = t148 * t127;
t89 = -rSges(5,3) * t138 + (rSges(5,1) * t137 - rSges(5,2) * t134) * t135;
t88 = -Icges(5,5) * t138 + (Icges(5,1) * t137 - Icges(5,4) * t134) * t135;
t87 = -Icges(5,6) * t138 + (Icges(5,4) * t137 - Icges(5,2) * t134) * t135;
t86 = -Icges(5,3) * t138 + (Icges(5,5) * t137 - Icges(5,6) * t134) * t135;
t85 = t127 * t134 + t128 * t157;
t84 = t127 * t137 - t128 * t158;
t83 = t127 * t157 - t128 * t134;
t82 = -t127 * t158 - t128 * t137;
t81 = t128 * t156 + t108;
t80 = t127 * t156 + t107;
t79 = rSges(4,3) * t127 + t128 * t147;
t78 = -rSges(4,3) * t128 + t127 * t147;
t77 = V_base(5) * rSges(2,3) - t120 * t129 + t154;
t76 = t121 * t129 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t69 = t120 * V_base(4) - t121 * V_base(5) + V_base(3);
t68 = V_base(5) * rSges(3,3) + (-t100 - t167) * t129 + t150;
t67 = t101 * t129 + (-rSges(3,3) + t165) * V_base(4) + t155;
t66 = V_base(4) * t100 + (-t101 - t166) * V_base(5) + t149;
t65 = rSges(5,1) * t85 + rSges(5,2) * t84 + rSges(5,3) * t159;
t64 = rSges(5,1) * t83 + rSges(5,2) * t82 + rSges(5,3) * t160;
t63 = Icges(5,1) * t85 + Icges(5,4) * t84 + Icges(5,5) * t159;
t62 = Icges(5,1) * t83 + Icges(5,4) * t82 + Icges(5,5) * t160;
t61 = Icges(5,4) * t85 + Icges(5,2) * t84 + Icges(5,6) * t159;
t60 = Icges(5,4) * t83 + Icges(5,2) * t82 + Icges(5,6) * t160;
t59 = Icges(5,5) * t85 + Icges(5,6) * t84 + Icges(5,3) * t159;
t58 = Icges(5,5) * t83 + Icges(5,6) * t82 + Icges(5,3) * t160;
t57 = t107 * t119 + (t151 - t78) * t129 + t150;
t56 = -t108 * t119 + t129 * t79 + t142;
t55 = -t107 * t79 + t108 * t78 + t141;
t54 = t107 * t122 - t118 * t64 + t80 * t89 + (t151 - t91) * t129 + t150;
t53 = -t108 * t122 + t118 * t65 + t129 * t92 - t81 * t89 + t142;
t52 = -t107 * t92 + t108 * t91 + t81 * t64 - t80 * t65 + t141;
t1 = m(1) * (t104 ^ 2 + t105 ^ 2 + t106 ^ 2) / 0.2e1 + m(2) * (t69 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + m(3) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(4) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + t108 * (t127 * t143 + t128 * t140) / 0.2e1 + t107 * (t127 * t140 - t128 * t143) / 0.2e1 + m(5) * (t52 ^ 2 + t53 ^ 2 + t54 ^ 2) / 0.2e1 + t81 * ((t59 * t159 + t84 * t61 + t85 * t63) * t81 + (t159 * t58 + t60 * t84 + t62 * t85) * t80 + (t159 * t86 + t84 * t87 + t85 * t88) * t118) / 0.2e1 + t80 * ((t160 * t59 + t61 * t82 + t83 * t63) * t81 + (t58 * t160 + t82 * t60 + t83 * t62) * t80 + (t160 * t86 + t82 * t87 + t83 * t88) * t118) / 0.2e1 + t118 * ((-t86 * t118 - t58 * t80 - t59 * t81) * t138 + ((-t134 * t61 + t137 * t63) * t81 + (-t134 * t60 + t137 * t62) * t80 + (-t134 * t87 + t137 * t88) * t118) * t135) / 0.2e1 + ((t135 * t75 + t138 * t73) * t108 + (t135 * t74 + t138 * t72) * t107 + (t138 * t112 + t135 * t115 + Icges(2,3) + Icges(3,3)) * t129) * t129 / 0.2e1 + ((-t136 * t113 + t116 * t139 - t127 * t96 + t128 * t98 + Icges(1,4)) * V_base(5) + (-t136 * t114 + t139 * t117 - t127 * t97 + t128 * t99 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t139 * t113 + t136 * t116 + t127 * t98 + t128 * t96 + Icges(1,2)) * V_base(5) + (t114 * t139 + t136 * t117 + t127 * t99 + t128 * t97 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t129 * (Icges(2,5) * t136 + Icges(3,5) * t127 + Icges(2,6) * t139 + Icges(3,6) * t128) + V_base(4) * t129 * (Icges(2,5) * t139 + Icges(3,5) * t128 - Icges(2,6) * t136 - Icges(3,6) * t127) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
