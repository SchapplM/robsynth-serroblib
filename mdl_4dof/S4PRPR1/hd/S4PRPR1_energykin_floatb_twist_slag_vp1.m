% Calculate kinetic energy for
% S4PRPR1
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4PRPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRPR1_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRPR1_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:42:11
% EndTime: 2018-11-14 13:42:12
% DurationCPUTime: 0.85s
% Computational Cost: add. (479->149), mult. (446->183), div. (0->0), fcn. (316->6), ass. (0->74)
t141 = Icges(3,4) - Icges(4,5);
t140 = Icges(3,1) + Icges(4,1);
t139 = Icges(3,2) + Icges(4,3);
t106 = pkin(6) + qJ(2);
t100 = sin(t106);
t138 = t141 * t100;
t101 = cos(t106);
t137 = t141 * t101;
t136 = -t139 * t101 - t138;
t135 = t139 * t100 - t137;
t134 = t140 * t100 + t137;
t133 = t140 * t101 - t138;
t132 = Icges(4,4) + Icges(3,5);
t131 = Icges(3,6) - Icges(4,6);
t130 = -pkin(5) - rSges(5,3);
t129 = cos(qJ(4));
t128 = sin(qJ(4));
t108 = cos(pkin(6));
t127 = pkin(1) * t108;
t126 = pkin(3) * t100;
t125 = -pkin(4) - qJ(1);
t61 = -t100 * t128 - t101 * t129;
t124 = Icges(5,4) * t61;
t107 = sin(pkin(6));
t123 = Icges(2,4) * t107;
t115 = pkin(1) * V_base(6);
t120 = t108 * t115 + V_base(2);
t119 = V_base(5) * qJ(1) + V_base(1);
t116 = qJD(1) + V_base(3);
t78 = pkin(2) * t101 + qJ(3) * t100;
t114 = -t78 - t127;
t103 = V_base(6) + qJD(2);
t113 = t103 * t78 + t120;
t112 = V_base(4) * t107 * pkin(1) + t116;
t75 = pkin(2) * t100 - qJ(3) * t101;
t111 = V_base(4) * t75 + t112;
t110 = V_base(5) * pkin(4) - t107 * t115 + t119;
t109 = qJD(3) * t100 + t110;
t102 = Icges(2,4) * t108;
t99 = -qJD(4) + t103;
t91 = rSges(2,1) * t108 - t107 * rSges(2,2);
t90 = t107 * rSges(2,1) + rSges(2,2) * t108;
t89 = Icges(2,1) * t108 - t123;
t88 = Icges(2,1) * t107 + t102;
t87 = -Icges(2,2) * t107 + t102;
t86 = Icges(2,2) * t108 + t123;
t83 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t82 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t81 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t80 = rSges(3,1) * t101 - rSges(3,2) * t100;
t79 = rSges(4,1) * t101 + rSges(4,3) * t100;
t77 = rSges(3,1) * t100 + rSges(3,2) * t101;
t76 = rSges(4,1) * t100 - rSges(4,3) * t101;
t62 = t100 * t129 - t101 * t128;
t59 = Icges(5,4) * t62;
t57 = V_base(5) * rSges(2,3) - t90 * V_base(6) + t119;
t56 = t91 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t55 = t90 * V_base(4) - t91 * V_base(5) + t116;
t54 = -rSges(5,1) * t61 + rSges(5,2) * t62;
t53 = rSges(5,1) * t62 + rSges(5,2) * t61;
t52 = -Icges(5,1) * t61 + t59;
t51 = Icges(5,1) * t62 + t124;
t50 = Icges(5,2) * t62 - t124;
t49 = Icges(5,2) * t61 + t59;
t46 = V_base(5) * rSges(3,3) - t103 * t77 + t110;
t45 = t103 * t80 + (-rSges(3,3) + t125) * V_base(4) + t120;
t44 = V_base(4) * t77 + (-t80 - t127) * V_base(5) + t112;
t43 = V_base(5) * rSges(4,2) + (-t75 - t76) * t103 + t109;
t42 = -qJD(3) * t101 + t103 * t79 + (-rSges(4,2) + t125) * V_base(4) + t113;
t41 = V_base(4) * t76 + (t114 - t79) * V_base(5) + t111;
t40 = -t53 * t99 + t130 * V_base(5) + (-t75 - t126) * t103 + t109;
t39 = t54 * t99 + (pkin(3) * t103 - qJD(3)) * t101 + (t125 - t130) * V_base(4) + t113;
t38 = (t53 + t126) * V_base(4) + (-pkin(3) * t101 + t114 - t54) * V_base(5) + t111;
t1 = m(1) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + Icges(1,1) * V_base(4) ^ 2 / 0.2e1 + Icges(1,2) * V_base(5) ^ 2 / 0.2e1 + m(2) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + m(3) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + m(4) * (t41 ^ 2 + t42 ^ 2 + t43 ^ 2) / 0.2e1 + m(5) * (t38 ^ 2 + t39 ^ 2 + t40 ^ 2) / 0.2e1 + V_base(4) * V_base(5) * Icges(1,4) + ((-Icges(5,5) * t62 - Icges(5,6) * t61) * V_base(5) + (Icges(5,5) * t61 - Icges(5,6) * t62) * V_base(4) + Icges(5,3) * t99 / 0.2e1) * t99 + ((Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6) + (Icges(2,5) * t107 + Icges(2,6) * t108 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t108 - Icges(2,6) * t107 + Icges(1,5)) * V_base(4)) * V_base(6) + ((Icges(3,3) / 0.2e1 + Icges(4,2) / 0.2e1) * t103 + (t131 * V_base(5) + t132 * V_base(4)) * t101 + (-t131 * V_base(4) + t132 * V_base(5)) * t100) * t103 + ((t136 * t100 + t134 * t101 - t107 * t86 + t108 * t88 + t49 * t62 - t51 * t61) * V_base(5) + (t135 * t100 + t133 * t101 - t107 * t87 + t108 * t89 + t50 * t62 - t52 * t61) * V_base(4)) * V_base(4) / 0.2e1 + ((t134 * t100 - t136 * t101 + t107 * t88 + t108 * t86 + t49 * t61 + t51 * t62) * V_base(5) + (t133 * t100 - t135 * t101 + t107 * t89 + t108 * t87 + t50 * t61 + t52 * t62) * V_base(4)) * V_base(5) / 0.2e1;
T  = t1;
