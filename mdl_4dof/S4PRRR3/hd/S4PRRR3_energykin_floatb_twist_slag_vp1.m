% Calculate kinetic energy for
% S4PRRR3
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR3_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR3_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR3_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRRR3_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR3_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR3_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR3_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR3_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:31:32
% EndTime: 2019-12-31 16:31:33
% DurationCPUTime: 0.92s
% Computational Cost: add. (638->173), mult. (488->229), div. (0->0), fcn. (316->8), ass. (0->87)
t118 = cos(pkin(7));
t148 = pkin(1) * t118;
t112 = V_base(6) + qJD(2);
t147 = pkin(2) * t112;
t146 = -pkin(4) - qJ(1);
t117 = sin(pkin(7));
t145 = Icges(2,4) * t117;
t116 = pkin(7) + qJ(2);
t108 = sin(t116);
t144 = Icges(3,4) * t108;
t111 = qJ(3) + t116;
t105 = sin(t111);
t143 = Icges(4,4) * t105;
t119 = sin(qJ(4));
t142 = Icges(5,4) * t119;
t120 = cos(qJ(4));
t141 = Icges(5,4) * t120;
t140 = -pkin(5) + t146;
t133 = pkin(1) * V_base(6);
t139 = t118 * t133 + V_base(2);
t138 = V_base(5) * qJ(1) + V_base(1);
t134 = qJD(1) + V_base(3);
t109 = cos(t116);
t132 = t109 * t147 + t139;
t131 = V_base(4) * t117 * pkin(1) + t134;
t130 = -pkin(2) * t109 - t148;
t129 = V_base(4) * pkin(2) * t108 + t131;
t128 = rSges(5,1) * t120 - rSges(5,2) * t119;
t127 = Icges(5,1) * t120 - t142;
t126 = -Icges(5,2) * t119 + t141;
t125 = Icges(5,5) * t120 - Icges(5,6) * t119;
t106 = cos(t111);
t107 = qJD(3) + t112;
t82 = -qJD(4) * t106 + V_base(5);
t83 = qJD(4) * t105 + V_base(4);
t124 = (Icges(5,5) * t119 + Icges(5,6) * t120) * t107 + (-Icges(5,3) * t106 + t105 * t125) * t82 + (Icges(5,3) * t105 + t106 * t125) * t83;
t123 = V_base(5) * pkin(4) - t117 * t133 + t138;
t122 = V_base(5) * pkin(5) - t108 * t147 + t123;
t56 = -Icges(5,6) * t106 + t105 * t126;
t57 = Icges(5,6) * t105 + t106 * t126;
t58 = -Icges(5,5) * t106 + t105 * t127;
t59 = Icges(5,5) * t105 + t106 * t127;
t97 = Icges(5,2) * t120 + t142;
t98 = Icges(5,1) * t119 + t141;
t121 = (-t119 * t57 + t120 * t59) * t83 + (-t119 * t56 + t120 * t58) * t82 + (-t119 * t97 + t120 * t98) * t107;
t110 = Icges(2,4) * t118;
t104 = Icges(3,4) * t109;
t101 = Icges(4,4) * t106;
t99 = t119 * rSges(5,1) + rSges(5,2) * t120;
t95 = rSges(2,1) * t118 - rSges(2,2) * t117;
t94 = rSges(2,1) * t117 + rSges(2,2) * t118;
t93 = Icges(2,1) * t118 - t145;
t92 = Icges(2,1) * t117 + t110;
t91 = -Icges(2,2) * t117 + t110;
t90 = Icges(2,2) * t118 + t145;
t86 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t85 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t84 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t81 = rSges(3,1) * t109 - rSges(3,2) * t108;
t80 = rSges(3,1) * t108 + rSges(3,2) * t109;
t79 = Icges(3,1) * t109 - t144;
t78 = Icges(3,1) * t108 + t104;
t77 = -Icges(3,2) * t108 + t104;
t76 = Icges(3,2) * t109 + t144;
t73 = pkin(3) * t106 + pkin(6) * t105;
t72 = pkin(3) * t105 - pkin(6) * t106;
t71 = rSges(4,1) * t106 - rSges(4,2) * t105;
t70 = rSges(4,1) * t105 + rSges(4,2) * t106;
t69 = Icges(4,1) * t106 - t143;
t68 = Icges(4,1) * t105 + t101;
t67 = -Icges(4,2) * t105 + t101;
t66 = Icges(4,2) * t106 + t143;
t63 = V_base(5) * rSges(2,3) - t94 * V_base(6) + t138;
t62 = t95 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t61 = t105 * rSges(5,3) + t106 * t128;
t60 = -t106 * rSges(5,3) + t105 * t128;
t53 = t94 * V_base(4) - t95 * V_base(5) + t134;
t52 = V_base(5) * rSges(3,3) - t112 * t80 + t123;
t51 = t112 * t81 + (-rSges(3,3) + t146) * V_base(4) + t139;
t50 = t80 * V_base(4) + (-t81 - t148) * V_base(5) + t131;
t49 = V_base(5) * rSges(4,3) - t107 * t70 + t122;
t48 = t107 * t71 + (-rSges(4,3) + t140) * V_base(4) + t132;
t47 = t70 * V_base(4) + (t130 - t71) * V_base(5) + t129;
t46 = t82 * t99 + (-t60 - t72) * t107 + t122;
t45 = -t83 * t99 + (t61 + t73) * t107 + t140 * V_base(4) + t132;
t44 = t60 * t83 - t61 * t82 + t72 * V_base(4) + (t130 - t73) * V_base(5) + t129;
t1 = m(1) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(2) * (t53 ^ 2 + t62 ^ 2 + t63 ^ 2) / 0.2e1 + m(3) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(4) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + m(5) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + t83 * (t105 * t124 + t106 * t121) / 0.2e1 + t82 * (t105 * t121 - t124 * t106) / 0.2e1 + ((t119 * t59 + t120 * t57) * t83 + (t119 * t58 + t120 * t56) * t82 + (t119 * t98 + t120 * t97 + Icges(4,3)) * t107) * t107 / 0.2e1 + ((-t105 * t66 + t106 * t68 - t108 * t76 + t109 * t78 - t117 * t90 + t118 * t92 + Icges(1,4)) * V_base(5) + (-t105 * t67 + t106 * t69 - t108 * t77 + t109 * t79 - t117 * t91 + t118 * t93 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t105 * t68 + t106 * t66 + t108 * t78 + t109 * t76 + t117 * t92 + t118 * t90 + Icges(1,2)) * V_base(5) + (t105 * t69 + t106 * t67 + t108 * t79 + t109 * t77 + t117 * t93 + t118 * t91 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t107 * (Icges(4,5) * t106 - Icges(4,6) * t105) + V_base(5) * t107 * (Icges(4,5) * t105 + Icges(4,6) * t106) + ((Icges(2,5) * t117 + Icges(2,6) * t118 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t118 - Icges(2,6) * t117 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + ((Icges(3,5) * t108 + Icges(3,6) * t109) * V_base(5) + (Icges(3,5) * t109 - Icges(3,6) * t108) * V_base(4) + Icges(3,3) * t112 / 0.2e1) * t112;
T = t1;
