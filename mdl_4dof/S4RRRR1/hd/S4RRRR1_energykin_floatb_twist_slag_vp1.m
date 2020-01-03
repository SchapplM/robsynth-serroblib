% Calculate kinetic energy for
% S4RRRR1
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRR1_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRR1_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:06
% EndTime: 2019-12-31 17:22:07
% DurationCPUTime: 0.78s
% Computational Cost: add. (670->173), mult. (488->231), div. (0->0), fcn. (316->8), ass. (0->88)
t149 = -pkin(4) - pkin(5);
t119 = sin(qJ(1));
t147 = pkin(1) * t119;
t121 = cos(qJ(1));
t146 = pkin(1) * t121;
t117 = qJ(1) + qJ(2);
t110 = sin(t117);
t145 = pkin(2) * t110;
t111 = cos(t117);
t144 = pkin(2) * t111;
t143 = Icges(2,4) * t119;
t142 = Icges(3,4) * t110;
t113 = qJ(3) + t117;
t107 = sin(t113);
t141 = Icges(4,4) * t107;
t118 = sin(qJ(4));
t140 = Icges(5,4) * t118;
t120 = cos(qJ(4));
t139 = Icges(5,4) * t120;
t138 = -pkin(6) + t149;
t109 = V_base(6) + qJD(1);
t137 = t109 * t146 + V_base(2);
t136 = V_base(4) * t147 + V_base(3);
t135 = V_base(5) * pkin(4) + V_base(1);
t106 = qJD(2) + t109;
t132 = t106 * t144 + t137;
t131 = V_base(4) * t145 + t136;
t130 = -t144 - t146;
t129 = rSges(5,1) * t120 - rSges(5,2) * t118;
t128 = Icges(5,1) * t120 - t140;
t127 = -Icges(5,2) * t118 + t139;
t126 = Icges(5,5) * t120 - Icges(5,6) * t118;
t125 = V_base(5) * pkin(5) - t109 * t147 + t135;
t103 = qJD(3) + t106;
t108 = cos(t113);
t86 = -qJD(4) * t108 + V_base(5);
t87 = qJD(4) * t107 + V_base(4);
t124 = (Icges(5,5) * t118 + Icges(5,6) * t120) * t103 + (-Icges(5,3) * t108 + t107 * t126) * t86 + (Icges(5,3) * t107 + t108 * t126) * t87;
t123 = V_base(5) * pkin(6) - t106 * t145 + t125;
t55 = -Icges(5,6) * t108 + t107 * t127;
t56 = Icges(5,6) * t107 + t108 * t127;
t57 = -Icges(5,5) * t108 + t107 * t128;
t58 = Icges(5,5) * t107 + t108 * t128;
t91 = Icges(5,2) * t120 + t140;
t94 = Icges(5,1) * t118 + t139;
t122 = (-t118 * t56 + t120 * t58) * t87 + (-t118 * t55 + t120 * t57) * t86 + (-t118 * t91 + t120 * t94) * t103;
t112 = Icges(2,4) * t121;
t105 = Icges(3,4) * t111;
t102 = Icges(4,4) * t108;
t99 = t121 * rSges(2,1) - t119 * rSges(2,2);
t98 = t119 * rSges(2,1) + t121 * rSges(2,2);
t97 = t118 * rSges(5,1) + t120 * rSges(5,2);
t96 = Icges(2,1) * t121 - t143;
t95 = Icges(2,1) * t119 + t112;
t93 = -Icges(2,2) * t119 + t112;
t92 = Icges(2,2) * t121 + t143;
t85 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t84 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t83 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t81 = t111 * rSges(3,1) - t110 * rSges(3,2);
t80 = t110 * rSges(3,1) + t111 * rSges(3,2);
t79 = Icges(3,1) * t111 - t142;
t78 = Icges(3,1) * t110 + t105;
t77 = -Icges(3,2) * t110 + t105;
t76 = Icges(3,2) * t111 + t142;
t73 = t108 * pkin(3) + t107 * pkin(7);
t72 = t107 * pkin(3) - t108 * pkin(7);
t71 = t108 * rSges(4,1) - t107 * rSges(4,2);
t70 = t107 * rSges(4,1) + t108 * rSges(4,2);
t69 = Icges(4,1) * t108 - t141;
t68 = Icges(4,1) * t107 + t102;
t67 = -Icges(4,2) * t107 + t102;
t66 = Icges(4,2) * t108 + t141;
t63 = V_base(5) * rSges(2,3) - t109 * t98 + t135;
t62 = t109 * t99 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t61 = t98 * V_base(4) - t99 * V_base(5) + V_base(3);
t60 = t107 * rSges(5,3) + t108 * t129;
t59 = -t108 * rSges(5,3) + t107 * t129;
t52 = V_base(5) * rSges(3,3) - t106 * t80 + t125;
t51 = t106 * t81 + (-rSges(3,3) + t149) * V_base(4) + t137;
t50 = V_base(4) * t80 + (-t81 - t146) * V_base(5) + t136;
t49 = V_base(5) * rSges(4,3) - t103 * t70 + t123;
t48 = t103 * t71 + (-rSges(4,3) + t138) * V_base(4) + t132;
t47 = V_base(4) * t70 + (t130 - t71) * V_base(5) + t131;
t46 = t86 * t97 + (-t59 - t72) * t103 + t123;
t45 = -t87 * t97 + (t60 + t73) * t103 + t138 * V_base(4) + t132;
t44 = t87 * t59 - t86 * t60 + V_base(4) * t72 + (t130 - t73) * V_base(5) + t131;
t1 = m(1) * (t83 ^ 2 + t84 ^ 2 + t85 ^ 2) / 0.2e1 + m(2) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + m(3) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(4) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + m(5) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + t87 * (t107 * t124 + t108 * t122) / 0.2e1 + t86 * (t107 * t122 - t124 * t108) / 0.2e1 + ((t118 * t58 + t120 * t56) * t87 + (t118 * t57 + t120 * t55) * t86 + (t118 * t94 + t120 * t91 + Icges(4,3)) * t103) * t103 / 0.2e1 + ((-t107 * t66 + t108 * t68 - t110 * t76 + t111 * t78 - t119 * t92 + t121 * t95 + Icges(1,4)) * V_base(5) + (-t107 * t67 + t108 * t69 - t110 * t77 + t111 * t79 - t119 * t93 + t121 * t96 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t107 * t68 + t108 * t66 + t110 * t78 + t111 * t76 + t119 * t95 + t121 * t92 + Icges(1,2)) * V_base(5) + (t107 * t69 + t108 * t67 + t110 * t79 + t111 * t77 + t119 * t96 + t121 * t93 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t103 * (Icges(4,5) * t108 - Icges(4,6) * t107) + V_base(5) * t103 * (Icges(4,5) * t107 + Icges(4,6) * t108) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t119 + Icges(2,6) * t121) * V_base(5) + (Icges(2,5) * t121 - Icges(2,6) * t119) * V_base(4) + Icges(2,3) * t109 / 0.2e1) * t109 + ((Icges(3,5) * t110 + Icges(3,6) * t111) * V_base(5) + (Icges(3,5) * t111 - Icges(3,6) * t110) * V_base(4) + Icges(3,3) * t106 / 0.2e1) * t106;
T = t1;
