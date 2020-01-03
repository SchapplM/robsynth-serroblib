% Calculate kinetic energy for
% S4PRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR4_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR4_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:54
% EndTime: 2019-12-31 16:21:55
% DurationCPUTime: 0.98s
% Computational Cost: add. (500->164), mult. (494->208), div. (0->0), fcn. (324->6), ass. (0->81)
t156 = Icges(3,4) + Icges(4,6);
t155 = Icges(3,1) + Icges(4,2);
t154 = -Icges(4,4) + Icges(3,5);
t153 = Icges(4,5) - Icges(3,6);
t152 = Icges(3,2) + Icges(4,3);
t105 = pkin(6) + qJ(2);
t100 = cos(t105);
t151 = t156 * t100;
t99 = sin(t105);
t150 = t156 * t99;
t149 = -t152 * t100 - t150;
t148 = t152 * t99 - t151;
t147 = t155 * t99 + t151;
t146 = t155 * t100 - t150;
t102 = V_base(6) + qJD(2);
t108 = sin(qJ(4));
t109 = cos(qJ(4));
t133 = Icges(5,4) * t108;
t114 = Icges(5,2) * t109 + t133;
t50 = Icges(5,6) * t100 + t114 * t99;
t51 = Icges(5,6) * t99 - t114 * t100;
t132 = Icges(5,4) * t109;
t115 = Icges(5,1) * t108 + t132;
t52 = Icges(5,5) * t100 + t115 * t99;
t53 = Icges(5,5) * t99 - t115 * t100;
t87 = qJD(4) * t99 + V_base(5);
t88 = qJD(4) * t100 + V_base(4);
t92 = -Icges(5,2) * t108 + t132;
t93 = Icges(5,1) * t109 - t133;
t141 = (t108 * t52 + t109 * t50) * t88 + (t108 * t53 + t109 * t51) * t87 + (t108 * t93 + t109 * t92) * t102;
t140 = pkin(5) * t99;
t107 = cos(pkin(6));
t138 = pkin(1) * t107;
t137 = -pkin(4) - qJ(1);
t106 = sin(pkin(6));
t134 = Icges(2,4) * t106;
t124 = pkin(1) * V_base(6);
t130 = t107 * t124 + V_base(2);
t129 = V_base(5) * qJ(1) + V_base(1);
t125 = qJD(1) + V_base(3);
t75 = pkin(2) * t100 + qJ(3) * t99;
t123 = -t75 - t138;
t122 = t102 * t75 + t130;
t121 = V_base(4) * t106 * pkin(1) + t125;
t72 = pkin(2) * t99 - qJ(3) * t100;
t120 = V_base(4) * t72 + t121;
t119 = rSges(5,1) * t108 + rSges(5,2) * t109;
t113 = Icges(5,5) * t108 + Icges(5,6) * t109;
t112 = t102 * (Icges(5,5) * t109 - Icges(5,6) * t108) + (Icges(5,3) * t100 + t113 * t99) * t88 + (Icges(5,3) * t99 - t113 * t100) * t87;
t111 = V_base(5) * pkin(4) - t106 * t124 + t129;
t110 = qJD(3) * t99 + t111;
t101 = Icges(2,4) * t107;
t94 = rSges(5,1) * t109 - t108 * rSges(5,2);
t90 = rSges(2,1) * t107 - rSges(2,2) * t106;
t89 = rSges(2,1) * t106 + rSges(2,2) * t107;
t86 = Icges(2,1) * t107 - t134;
t85 = Icges(2,1) * t106 + t101;
t84 = -Icges(2,2) * t106 + t101;
t83 = Icges(2,2) * t107 + t134;
t80 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t79 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t78 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t77 = rSges(3,1) * t100 - rSges(3,2) * t99;
t76 = -rSges(4,2) * t100 + rSges(4,3) * t99;
t74 = rSges(3,1) * t99 + rSges(3,2) * t100;
t73 = -rSges(4,2) * t99 - rSges(4,3) * t100;
t57 = V_base(5) * rSges(2,3) - t89 * V_base(6) + t129;
t56 = t90 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t55 = t99 * rSges(5,3) - t119 * t100;
t54 = t100 * rSges(5,3) + t119 * t99;
t47 = t89 * V_base(4) - t90 * V_base(5) + t125;
t46 = V_base(5) * rSges(3,3) - t102 * t74 + t111;
t45 = t102 * t77 + (-rSges(3,3) + t137) * V_base(4) + t130;
t44 = t74 * V_base(4) + (-t77 - t138) * V_base(5) + t121;
t43 = V_base(5) * rSges(4,1) + (-t72 - t73) * t102 + t110;
t42 = -qJD(3) * t100 + t102 * t76 + (-rSges(4,1) + t137) * V_base(4) + t122;
t41 = t73 * V_base(4) + (t123 - t76) * V_base(5) + t120;
t40 = V_base(5) * pkin(3) + t87 * t94 + (-t55 - t72 - t140) * t102 + t110;
t39 = t102 * t54 - t88 * t94 + (pkin(5) * t102 - qJD(3)) * t100 + (-pkin(3) + t137) * V_base(4) + t122;
t38 = V_base(4) * t140 - t54 * t87 + t55 * t88 + (-pkin(5) * t100 + t123) * V_base(5) + t120;
t1 = m(1) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + m(2) * (t47 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + m(3) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + m(4) * (t41 ^ 2 + t42 ^ 2 + t43 ^ 2) / 0.2e1 + m(5) * (t38 ^ 2 + t39 ^ 2 + t40 ^ 2) / 0.2e1 + t88 * (t112 * t100 + t141 * t99) / 0.2e1 + t87 * (-t141 * t100 + t112 * t99) / 0.2e1 + ((-t108 * t50 + t109 * t52) * t88 + (-t108 * t51 + t109 * t53) * t87 + (-t108 * t92 + t109 * t93 + Icges(4,1) + Icges(3,3)) * t102) * t102 / 0.2e1 + ((t147 * t100 - t106 * t83 + t107 * t85 + t149 * t99 + Icges(1,4)) * V_base(5) + (t146 * t100 - t106 * t84 + t107 * t86 + t148 * t99 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t149 * t100 + t106 * t85 + t107 * t83 + t147 * t99 + Icges(1,2)) * V_base(5) + (-t148 * t100 + t106 * t86 + t107 * t84 + t146 * t99 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t102 * (-t153 * t100 + t154 * t99) + V_base(4) * t102 * (t154 * t100 + t153 * t99) + ((Icges(2,5) * t106 + Icges(2,6) * t107 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t107 - Icges(2,6) * t106 + Icges(1,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6);
T = t1;
