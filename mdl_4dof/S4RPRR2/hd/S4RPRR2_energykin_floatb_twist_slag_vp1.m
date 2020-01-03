% Calculate kinetic energy for
% S4RPRR2
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
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRR2_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRR2_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:06
% EndTime: 2019-12-31 16:48:07
% DurationCPUTime: 0.92s
% Computational Cost: add. (649->174), mult. (488->226), div. (0->0), fcn. (316->8), ass. (0->87)
t118 = sin(qJ(1));
t147 = pkin(1) * t118;
t120 = cos(qJ(1));
t146 = pkin(1) * t120;
t116 = qJ(1) + pkin(7);
t108 = sin(t116);
t145 = pkin(2) * t108;
t109 = cos(t116);
t144 = pkin(2) * t109;
t143 = -pkin(4) - qJ(2);
t142 = Icges(2,4) * t118;
t141 = Icges(3,4) * t108;
t110 = qJ(3) + t116;
t105 = sin(t110);
t140 = Icges(4,4) * t105;
t117 = sin(qJ(4));
t139 = Icges(5,4) * t117;
t119 = cos(qJ(4));
t138 = Icges(5,4) * t119;
t137 = -pkin(5) + t143;
t111 = V_base(6) + qJD(1);
t136 = t111 * t146 + V_base(2);
t135 = V_base(5) * pkin(4) + V_base(1);
t132 = t111 * t144 + t136;
t131 = V_base(5) * qJ(2) + t135;
t130 = V_base(4) * t147 + qJD(2) + V_base(3);
t129 = -t144 - t146;
t128 = V_base(4) * t145 + t130;
t127 = rSges(5,1) * t119 - rSges(5,2) * t117;
t126 = Icges(5,1) * t119 - t139;
t125 = -Icges(5,2) * t117 + t138;
t124 = Icges(5,5) * t119 - Icges(5,6) * t117;
t106 = cos(t110);
t107 = qJD(3) + t111;
t82 = -qJD(4) * t106 + V_base(5);
t83 = qJD(4) * t105 + V_base(4);
t123 = (Icges(5,5) * t117 + Icges(5,6) * t119) * t107 + (-Icges(5,3) * t106 + t105 * t124) * t82 + (Icges(5,3) * t105 + t106 * t124) * t83;
t122 = V_base(5) * pkin(5) + (-t145 - t147) * t111 + t131;
t55 = -Icges(5,6) * t106 + t105 * t125;
t56 = Icges(5,6) * t105 + t106 * t125;
t57 = -Icges(5,5) * t106 + t105 * t126;
t58 = Icges(5,5) * t105 + t106 * t126;
t91 = Icges(5,2) * t119 + t139;
t94 = Icges(5,1) * t117 + t138;
t121 = (-t117 * t56 + t119 * t58) * t83 + (-t117 * t55 + t119 * t57) * t82 + (-t117 * t91 + t119 * t94) * t107;
t113 = Icges(2,4) * t120;
t104 = Icges(3,4) * t109;
t102 = Icges(4,4) * t106;
t99 = rSges(2,1) * t120 - t118 * rSges(2,2);
t98 = t118 * rSges(2,1) + rSges(2,2) * t120;
t97 = rSges(5,1) * t117 + rSges(5,2) * t119;
t96 = Icges(2,1) * t120 - t142;
t95 = Icges(2,1) * t118 + t113;
t93 = -Icges(2,2) * t118 + t113;
t92 = Icges(2,2) * t120 + t142;
t86 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t85 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t84 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t81 = rSges(3,1) * t109 - rSges(3,2) * t108;
t80 = rSges(3,1) * t108 + rSges(3,2) * t109;
t79 = Icges(3,1) * t109 - t141;
t78 = Icges(3,1) * t108 + t104;
t77 = -Icges(3,2) * t108 + t104;
t76 = Icges(3,2) * t109 + t141;
t73 = pkin(3) * t106 + pkin(6) * t105;
t72 = pkin(3) * t105 - pkin(6) * t106;
t71 = rSges(4,1) * t106 - rSges(4,2) * t105;
t70 = rSges(4,1) * t105 + rSges(4,2) * t106;
t69 = Icges(4,1) * t106 - t140;
t68 = Icges(4,1) * t105 + t102;
t67 = -Icges(4,2) * t105 + t102;
t66 = Icges(4,2) * t106 + t140;
t63 = V_base(5) * rSges(2,3) - t111 * t98 + t135;
t62 = t111 * t99 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t61 = t98 * V_base(4) - t99 * V_base(5) + V_base(3);
t60 = rSges(5,3) * t105 + t106 * t127;
t59 = -rSges(5,3) * t106 + t105 * t127;
t52 = V_base(5) * rSges(3,3) + (-t80 - t147) * t111 + t131;
t51 = t111 * t81 + (-rSges(3,3) + t143) * V_base(4) + t136;
t50 = V_base(4) * t80 + (-t81 - t146) * V_base(5) + t130;
t49 = V_base(5) * rSges(4,3) - t107 * t70 + t122;
t48 = t107 * t71 + (-rSges(4,3) + t137) * V_base(4) + t132;
t47 = V_base(4) * t70 + (t129 - t71) * V_base(5) + t128;
t46 = t82 * t97 + (-t59 - t72) * t107 + t122;
t45 = -t83 * t97 + (t60 + t73) * t107 + t137 * V_base(4) + t132;
t44 = t83 * t59 - t82 * t60 + V_base(4) * t72 + (t129 - t73) * V_base(5) + t128;
t1 = m(1) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(2) * (t61 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + m(3) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(4) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + m(5) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + t83 * (t105 * t123 + t106 * t121) / 0.2e1 + t82 * (t105 * t121 - t123 * t106) / 0.2e1 + ((t117 * t58 + t119 * t56) * t83 + (t117 * t57 + t119 * t55) * t82 + (t117 * t94 + t119 * t91 + Icges(4,3)) * t107) * t107 / 0.2e1 + ((-t105 * t66 + t106 * t68 - t108 * t76 + t109 * t78 - t118 * t92 + t120 * t95 + Icges(1,4)) * V_base(5) + (-t105 * t67 + t106 * t69 - t108 * t77 + t109 * t79 - t118 * t93 + t120 * t96 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t105 * t68 + t106 * t66 + t108 * t78 + t109 * t76 + t118 * t95 + t120 * t92 + Icges(1,2)) * V_base(5) + (t105 * t69 + t106 * t67 + t108 * t79 + t109 * t77 + t118 * t96 + t120 * t93 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t107 * (Icges(4,5) * t106 - Icges(4,6) * t105) + V_base(5) * t107 * (Icges(4,5) * t105 + Icges(4,6) * t106) + ((Icges(2,5) * t120 + Icges(3,5) * t109 - Icges(2,6) * t118 - Icges(3,6) * t108) * V_base(4) + (Icges(2,5) * t118 + Icges(3,5) * t108 + Icges(2,6) * t120 + Icges(3,6) * t109) * V_base(5) + (Icges(2,3) / 0.2e1 + Icges(3,3) / 0.2e1) * t111) * t111 + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
