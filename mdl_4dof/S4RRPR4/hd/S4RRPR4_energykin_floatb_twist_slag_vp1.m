% Calculate kinetic energy for
% S4RRPR4
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
% Datum: 2019-12-31 17:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR4_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR4_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR4_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR4_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR4_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR4_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:02:24
% EndTime: 2019-12-31 17:02:25
% DurationCPUTime: 1.00s
% Computational Cost: add. (725->193), mult. (628->261), div. (0->0), fcn. (458->8), ass. (0->99)
t162 = -pkin(4) - pkin(5);
t126 = sin(qJ(1));
t160 = pkin(1) * t126;
t127 = cos(qJ(1));
t159 = pkin(1) * t127;
t123 = sin(pkin(7));
t158 = pkin(3) * t123;
t124 = cos(pkin(7));
t157 = pkin(3) * t124;
t156 = Icges(2,4) * t126;
t122 = qJ(1) + qJ(2);
t116 = sin(t122);
t155 = Icges(3,4) * t116;
t154 = Icges(4,4) * t123;
t153 = Icges(4,4) * t124;
t121 = pkin(7) + qJ(4);
t113 = sin(t121);
t152 = Icges(5,4) * t113;
t114 = cos(t121);
t151 = Icges(5,4) * t114;
t115 = V_base(6) + qJD(1);
t149 = t115 * t159 + V_base(2);
t148 = V_base(4) * t160 + V_base(3);
t147 = V_base(5) * pkin(4) + V_base(1);
t117 = cos(t122);
t88 = pkin(2) * t117 + qJ(3) * t116;
t144 = -t88 - t159;
t86 = pkin(2) * t116 - qJ(3) * t117;
t143 = V_base(4) * t86 + t148;
t142 = rSges(4,1) * t124 - rSges(4,2) * t123;
t141 = rSges(5,1) * t114 - rSges(5,2) * t113;
t140 = Icges(4,1) * t124 - t154;
t139 = Icges(5,1) * t114 - t152;
t138 = -Icges(4,2) * t123 + t153;
t137 = -Icges(5,2) * t113 + t151;
t136 = Icges(4,5) * t124 - Icges(4,6) * t123;
t135 = Icges(5,5) * t114 - Icges(5,6) * t113;
t111 = qJD(2) + t115;
t134 = -qJD(3) * t117 + t111 * t88 + t149;
t133 = V_base(5) * pkin(5) - t115 * t160 + t147;
t97 = -qJD(4) * t117 + V_base(5);
t98 = qJD(4) * t116 + V_base(4);
t132 = (Icges(5,5) * t113 + Icges(5,6) * t114) * t111 + (-Icges(5,3) * t117 + t116 * t135) * t97 + (Icges(5,3) * t116 + t117 * t135) * t98;
t131 = qJD(3) * t116 + t133;
t130 = (Icges(4,5) * t123 + Icges(4,6) * t124) * t111 + (-Icges(4,3) * t117 + t116 * t136) * V_base(5) + (Icges(4,3) * t116 + t117 * t136) * V_base(4);
t57 = -Icges(5,6) * t117 + t116 * t137;
t58 = Icges(5,6) * t116 + t117 * t137;
t59 = -Icges(5,5) * t117 + t116 * t139;
t60 = Icges(5,5) * t116 + t117 * t139;
t77 = Icges(5,2) * t114 + t152;
t78 = Icges(5,1) * t113 + t151;
t129 = (-t113 * t58 + t114 * t60) * t98 + (-t113 * t57 + t114 * t59) * t97 + (-t113 * t77 + t114 * t78) * t111;
t68 = -Icges(4,6) * t117 + t116 * t138;
t69 = Icges(4,6) * t116 + t117 * t138;
t70 = -Icges(4,5) * t117 + t116 * t140;
t71 = Icges(4,5) * t116 + t117 * t140;
t94 = Icges(4,2) * t124 + t154;
t95 = Icges(4,1) * t123 + t153;
t128 = (-t123 * t69 + t124 * t71) * V_base(4) + (-t123 * t68 + t124 * t70) * V_base(5) + (-t123 * t94 + t124 * t95) * t111;
t118 = Icges(2,4) * t127;
t110 = Icges(3,4) * t117;
t106 = rSges(2,1) * t127 - t126 * rSges(2,2);
t105 = t126 * rSges(2,1) + rSges(2,2) * t127;
t104 = Icges(2,1) * t127 - t156;
t103 = Icges(2,1) * t126 + t118;
t102 = -Icges(2,2) * t126 + t118;
t101 = Icges(2,2) * t127 + t156;
t96 = rSges(4,1) * t123 + rSges(4,2) * t124;
t92 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t91 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t90 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t89 = rSges(3,1) * t117 - rSges(3,2) * t116;
t87 = rSges(3,1) * t116 + rSges(3,2) * t117;
t85 = Icges(3,1) * t117 - t155;
t84 = Icges(3,1) * t116 + t110;
t83 = -Icges(3,2) * t116 + t110;
t82 = Icges(3,2) * t117 + t155;
t81 = Icges(3,5) * t117 - Icges(3,6) * t116;
t80 = Icges(3,5) * t116 + Icges(3,6) * t117;
t79 = rSges(5,1) * t113 + rSges(5,2) * t114;
t73 = rSges(4,3) * t116 + t117 * t142;
t72 = -rSges(4,3) * t117 + t116 * t142;
t65 = V_base(5) * rSges(2,3) - t105 * t115 + t147;
t64 = t106 * t115 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t63 = t105 * V_base(4) - t106 * V_base(5) + V_base(3);
t62 = rSges(5,3) * t116 + t117 * t141;
t61 = -rSges(5,3) * t117 + t116 * t141;
t54 = pkin(6) * t116 + t117 * t157;
t53 = -pkin(6) * t117 + t116 * t157;
t52 = V_base(5) * rSges(3,3) - t111 * t87 + t133;
t51 = t111 * t89 + (-rSges(3,3) + t162) * V_base(4) + t149;
t50 = V_base(4) * t87 + (-t89 - t159) * V_base(5) + t148;
t49 = t96 * V_base(5) + (-t72 - t86) * t111 + t131;
t48 = t111 * t73 + (-t96 + t162) * V_base(4) + t134;
t47 = V_base(4) * t72 + (t144 - t73) * V_base(5) + t143;
t46 = V_base(5) * t158 + t79 * t97 + (-t53 - t61 - t86) * t111 + t131;
t45 = -t79 * t98 + (t54 + t62) * t111 + (-t158 + t162) * V_base(4) + t134;
t44 = V_base(4) * t53 + t98 * t61 - t97 * t62 + (t144 - t54) * V_base(5) + t143;
t1 = m(1) * (t90 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(2) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + m(3) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(4) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + m(5) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + t98 * (t116 * t132 + t117 * t129) / 0.2e1 + t97 * (t116 * t129 - t132 * t117) / 0.2e1 + ((t113 * t60 + t114 * t58) * t98 + (t113 * t59 + t114 * t57) * t97 + (t123 * t70 + t124 * t68 + t80) * V_base(5) + (t123 * t71 + t124 * t69 + t81) * V_base(4) + (t113 * t78 + t114 * t77 + t123 * t95 + t124 * t94 + Icges(3,3)) * t111) * t111 / 0.2e1 + (t81 * t111 + t116 * t130 + t117 * t128 + (-t126 * t101 + t103 * t127 - t116 * t82 + t117 * t84 + Icges(1,4)) * V_base(5) + (-t126 * t102 + t127 * t104 - t116 * t83 + t117 * t85 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + (t80 * t111 + t116 * t128 - t117 * t130 + (t127 * t101 + t126 * t103 + t116 * t84 + t117 * t82 + Icges(1,2)) * V_base(5) + (t102 * t127 + t126 * t104 + t116 * t85 + t117 * t83 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t127 - Icges(2,6) * t126) * V_base(4) + (Icges(2,5) * t126 + Icges(2,6) * t127) * V_base(5) + Icges(2,3) * t115 / 0.2e1) * t115;
T = t1;
