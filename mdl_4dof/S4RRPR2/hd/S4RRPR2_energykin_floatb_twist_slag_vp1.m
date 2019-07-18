% Calculate kinetic energy for
% S4RRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR2_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR2_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRPR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPR2_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR2_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR2_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR2_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 18:16:29
% EndTime: 2019-07-18 18:16:29
% DurationCPUTime: 0.86s
% Computational Cost: add. (509->148), mult. (444->186), div. (0->0), fcn. (316->6), ass. (0->73)
t140 = Icges(3,4) - Icges(4,5);
t139 = Icges(3,1) + Icges(4,1);
t138 = Icges(3,2) + Icges(4,3);
t107 = qJ(1) + qJ(2);
t102 = sin(t107);
t137 = t140 * t102;
t103 = cos(t107);
t136 = t140 * t103;
t135 = -t138 * t103 - t137;
t134 = t138 * t102 - t136;
t133 = t139 * t102 + t136;
t132 = t139 * t103 - t137;
t131 = Icges(4,4) + Icges(3,5);
t130 = Icges(3,6) - Icges(4,6);
t129 = -pkin(4) - pkin(5);
t128 = cos(qJ(4));
t127 = sin(qJ(4));
t108 = sin(qJ(1));
t126 = pkin(1) * t108;
t109 = cos(qJ(1));
t125 = pkin(1) * t109;
t124 = pkin(3) * t102;
t61 = -t102 * t127 - t103 * t128;
t123 = Icges(5,4) * t61;
t122 = Icges(2,4) * t108;
t101 = V_base(6) + qJD(1);
t119 = t101 * t125 + V_base(2);
t118 = V_base(4) * t126 + V_base(3);
t117 = V_base(5) * pkin(4) + V_base(1);
t78 = pkin(2) * t103 + qJ(3) * t102;
t114 = -t78 - t125;
t100 = qJD(2) + t101;
t113 = t100 * t78 + t119;
t75 = pkin(2) * t102 - qJ(3) * t103;
t112 = V_base(4) * t75 + t118;
t111 = V_base(5) * pkin(5) - t101 * t126 + t117;
t110 = qJD(3) * t102 + t111;
t104 = Icges(2,4) * t109;
t95 = -qJD(4) + t100;
t91 = rSges(2,1) * t109 - t108 * rSges(2,2);
t90 = t108 * rSges(2,1) + rSges(2,2) * t109;
t89 = Icges(2,1) * t109 - t122;
t88 = Icges(2,1) * t108 + t104;
t87 = -Icges(2,2) * t108 + t104;
t86 = Icges(2,2) * t109 + t122;
t83 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t82 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t81 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t80 = rSges(3,1) * t103 - rSges(3,2) * t102;
t79 = rSges(4,1) * t103 + rSges(4,3) * t102;
t77 = rSges(3,1) * t102 + rSges(3,2) * t103;
t76 = rSges(4,1) * t102 - rSges(4,3) * t103;
t62 = t102 * t128 - t103 * t127;
t59 = Icges(5,4) * t62;
t57 = V_base(5) * rSges(2,3) - t101 * t90 + t117;
t56 = t101 * t91 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t55 = t90 * V_base(4) - t91 * V_base(5) + V_base(3);
t54 = -rSges(5,1) * t61 + rSges(5,2) * t62;
t53 = rSges(5,1) * t62 + rSges(5,2) * t61;
t52 = -Icges(5,1) * t61 + t59;
t51 = Icges(5,1) * t62 + t123;
t50 = Icges(5,2) * t62 - t123;
t49 = Icges(5,2) * t61 + t59;
t46 = V_base(5) * rSges(3,3) - t100 * t77 + t111;
t45 = t100 * t80 + (-rSges(3,3) + t129) * V_base(4) + t119;
t44 = V_base(4) * t77 + (-t80 - t125) * V_base(5) + t118;
t43 = V_base(5) * rSges(4,2) + (-t75 - t76) * t100 + t110;
t42 = -qJD(3) * t103 + t100 * t79 + (-rSges(4,2) + t129) * V_base(4) + t113;
t41 = V_base(4) * t76 + (t114 - t79) * V_base(5) + t112;
t40 = -V_base(5) * rSges(5,3) - t53 * t95 + (-t75 - t124) * t100 + t110;
t39 = t54 * t95 + (pkin(3) * t100 - qJD(3)) * t103 + (rSges(5,3) + t129) * V_base(4) + t113;
t38 = (t53 + t124) * V_base(4) + (-pkin(3) * t103 + t114 - t54) * V_base(5) + t112;
t1 = m(1) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + Icges(1,1) * V_base(4) ^ 2 / 0.2e1 + Icges(1,2) * V_base(5) ^ 2 / 0.2e1 + m(2) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + m(3) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + m(4) * (t41 ^ 2 + t42 ^ 2 + t43 ^ 2) / 0.2e1 + m(5) * (t38 ^ 2 + t39 ^ 2 + t40 ^ 2) / 0.2e1 + V_base(4) * V_base(5) * Icges(1,4) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((-Icges(5,5) * t62 - Icges(5,6) * t61) * V_base(5) + (Icges(5,5) * t61 - Icges(5,6) * t62) * V_base(4) + Icges(5,3) * t95 / 0.2e1) * t95 + ((Icges(2,5) * t108 + Icges(2,6) * t109) * V_base(5) + (Icges(2,5) * t109 - Icges(2,6) * t108) * V_base(4) + Icges(2,3) * t101 / 0.2e1) * t101 + ((Icges(3,3) / 0.2e1 + Icges(4,2) / 0.2e1) * t100 + (t130 * V_base(5) + t131 * V_base(4)) * t103 + (-t130 * V_base(4) + t131 * V_base(5)) * t102) * t100 + ((t135 * t102 + t133 * t103 - t108 * t86 + t109 * t88 + t49 * t62 - t51 * t61) * V_base(5) + (t134 * t102 + t132 * t103 - t108 * t87 + t109 * t89 + t62 * t50 - t61 * t52) * V_base(4)) * V_base(4) / 0.2e1 + ((t133 * t102 - t135 * t103 + t108 * t88 + t109 * t86 + t61 * t49 + t62 * t51) * V_base(5) + (t132 * t102 - t134 * t103 + t108 * t89 + t109 * t87 + t50 * t61 + t52 * t62) * V_base(4)) * V_base(5) / 0.2e1;
T  = t1;
