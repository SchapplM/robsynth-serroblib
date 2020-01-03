% Calculate kinetic energy for
% S4RRPR3
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRPR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR3_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR3_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:27
% EndTime: 2019-12-31 17:01:27
% DurationCPUTime: 0.86s
% Computational Cost: add. (658->174), mult. (488->225), div. (0->0), fcn. (316->8), ass. (0->87)
t148 = -pkin(4) - pkin(5);
t118 = sin(qJ(1));
t146 = pkin(1) * t118;
t120 = cos(qJ(1));
t145 = pkin(1) * t120;
t116 = qJ(1) + qJ(2);
t110 = sin(t116);
t144 = pkin(2) * t110;
t111 = cos(t116);
t143 = pkin(2) * t111;
t142 = Icges(2,4) * t118;
t141 = Icges(3,4) * t110;
t108 = pkin(7) + t116;
t105 = sin(t108);
t140 = Icges(4,4) * t105;
t117 = sin(qJ(4));
t139 = Icges(5,4) * t117;
t119 = cos(qJ(4));
t138 = Icges(5,4) * t119;
t137 = -qJ(3) + t148;
t109 = V_base(6) + qJD(1);
t136 = t109 * t145 + V_base(2);
t135 = V_base(4) * t146 + V_base(3);
t134 = V_base(5) * pkin(4) + V_base(1);
t107 = qJD(2) + t109;
t131 = t107 * t143 + t136;
t130 = -t143 - t145;
t129 = V_base(4) * t144 + qJD(3) + t135;
t128 = rSges(5,1) * t119 - rSges(5,2) * t117;
t127 = Icges(5,1) * t119 - t139;
t126 = -Icges(5,2) * t117 + t138;
t125 = Icges(5,5) * t119 - Icges(5,6) * t117;
t124 = V_base(5) * pkin(5) - t109 * t146 + t134;
t106 = cos(t108);
t83 = -qJD(4) * t106 + V_base(5);
t84 = qJD(4) * t105 + V_base(4);
t123 = t107 * (Icges(5,5) * t117 + Icges(5,6) * t119) + (-Icges(5,3) * t106 + t105 * t125) * t83 + (Icges(5,3) * t105 + t106 * t125) * t84;
t122 = V_base(5) * qJ(3) + t124;
t55 = -Icges(5,6) * t106 + t105 * t126;
t56 = Icges(5,6) * t105 + t106 * t126;
t57 = -Icges(5,5) * t106 + t105 * t127;
t58 = Icges(5,5) * t105 + t106 * t127;
t91 = Icges(5,2) * t119 + t139;
t94 = Icges(5,1) * t117 + t138;
t121 = (-t117 * t56 + t119 * t58) * t84 + (-t117 * t55 + t119 * t57) * t83 + (-t117 * t91 + t119 * t94) * t107;
t113 = Icges(2,4) * t120;
t104 = Icges(3,4) * t111;
t102 = Icges(4,4) * t106;
t99 = rSges(2,1) * t120 - t118 * rSges(2,2);
t98 = t118 * rSges(2,1) + rSges(2,2) * t120;
t97 = rSges(5,1) * t117 + rSges(5,2) * t119;
t96 = Icges(2,1) * t120 - t142;
t95 = Icges(2,1) * t118 + t113;
t93 = -Icges(2,2) * t118 + t113;
t92 = Icges(2,2) * t120 + t142;
t87 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t86 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t85 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t81 = rSges(3,1) * t111 - rSges(3,2) * t110;
t80 = rSges(3,1) * t110 + rSges(3,2) * t111;
t79 = Icges(3,1) * t111 - t141;
t78 = Icges(3,1) * t110 + t104;
t77 = -Icges(3,2) * t110 + t104;
t76 = Icges(3,2) * t111 + t141;
t73 = pkin(3) * t106 + pkin(6) * t105;
t72 = pkin(3) * t105 - pkin(6) * t106;
t71 = rSges(4,1) * t106 - rSges(4,2) * t105;
t70 = rSges(4,1) * t105 + rSges(4,2) * t106;
t69 = Icges(4,1) * t106 - t140;
t68 = Icges(4,1) * t105 + t102;
t67 = -Icges(4,2) * t105 + t102;
t66 = Icges(4,2) * t106 + t140;
t63 = V_base(5) * rSges(2,3) - t109 * t98 + t134;
t62 = t109 * t99 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t61 = t98 * V_base(4) - t99 * V_base(5) + V_base(3);
t60 = rSges(5,3) * t105 + t106 * t128;
t59 = -rSges(5,3) * t106 + t105 * t128;
t52 = V_base(5) * rSges(3,3) - t107 * t80 + t124;
t51 = t107 * t81 + (-rSges(3,3) + t148) * V_base(4) + t136;
t50 = V_base(4) * t80 + (-t81 - t145) * V_base(5) + t135;
t49 = V_base(5) * rSges(4,3) + (-t70 - t144) * t107 + t122;
t48 = t107 * t71 + (-rSges(4,3) + t137) * V_base(4) + t131;
t47 = V_base(4) * t70 + (t130 - t71) * V_base(5) + t129;
t46 = t83 * t97 + (-t59 - t72 - t144) * t107 + t122;
t45 = -t84 * t97 + (t60 + t73) * t107 + t137 * V_base(4) + t131;
t44 = t84 * t59 - t83 * t60 + V_base(4) * t72 + (t130 - t73) * V_base(5) + t129;
t1 = m(1) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(2) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + m(3) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(4) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + m(5) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + t84 * (t123 * t105 + t121 * t106) / 0.2e1 + t83 * (t121 * t105 - t123 * t106) / 0.2e1 + ((t117 * t58 + t119 * t56) * t84 + (t117 * t57 + t119 * t55) * t83 + (t117 * t94 + t119 * t91 + Icges(3,3) + Icges(4,3)) * t107) * t107 / 0.2e1 + ((-t105 * t66 + t106 * t68 - t110 * t76 + t111 * t78 - t118 * t92 + t120 * t95 + Icges(1,4)) * V_base(5) + (-t105 * t67 + t106 * t69 - t110 * t77 + t111 * t79 - t118 * t93 + t120 * t96 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t105 * t68 + t106 * t66 + t110 * t78 + t111 * t76 + t118 * t95 + t120 * t92 + Icges(1,2)) * V_base(5) + (t105 * t69 + t106 * t67 + t110 * t79 + t111 * t77 + t118 * t96 + t120 * t93 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t107 * (Icges(3,5) * t110 + Icges(4,5) * t105 + Icges(3,6) * t111 + Icges(4,6) * t106) + V_base(4) * t107 * (Icges(3,5) * t111 + Icges(4,5) * t106 - Icges(3,6) * t110 - Icges(4,6) * t105) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t118 + Icges(2,6) * t120) * V_base(5) + (Icges(2,5) * t120 - Icges(2,6) * t118) * V_base(4) + Icges(2,3) * t109 / 0.2e1) * t109;
T = t1;
