% Calculate kinetic energy for
% S4RPPR1
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
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2019-03-08 18:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR1_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR1_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR1_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:27:24
% EndTime: 2019-03-08 18:27:25
% DurationCPUTime: 0.77s
% Computational Cost: add. (490->150), mult. (446->178), div. (0->0), fcn. (316->6), ass. (0->76)
t143 = Icges(3,4) - Icges(4,5);
t142 = Icges(3,1) + Icges(4,1);
t141 = Icges(3,2) + Icges(4,3);
t106 = qJ(1) + pkin(6);
t100 = sin(t106);
t140 = t143 * t100;
t101 = cos(t106);
t139 = t143 * t101;
t138 = -t141 * t101 - t140;
t137 = t141 * t100 - t139;
t136 = t142 * t100 + t139;
t135 = t142 * t101 - t140;
t134 = V_base(4) / 0.2e1;
t133 = V_base(5) / 0.2e1;
t132 = Icges(4,4) + Icges(3,5);
t131 = Icges(3,6) - Icges(4,6);
t130 = -pkin(5) - rSges(5,3);
t129 = cos(qJ(4));
t128 = sin(qJ(4));
t107 = sin(qJ(1));
t127 = pkin(1) * t107;
t108 = cos(qJ(1));
t126 = pkin(1) * t108;
t125 = pkin(3) * t100;
t124 = -pkin(4) - qJ(2);
t61 = -t100 * t128 - t101 * t129;
t123 = Icges(5,4) * t61;
t122 = Icges(2,4) * t107;
t102 = V_base(6) + qJD(1);
t119 = t102 * t126 + V_base(2);
t118 = V_base(5) * pkin(4) + V_base(1);
t75 = t100 * pkin(2) - t101 * qJ(3);
t115 = -t75 - t127;
t78 = t101 * pkin(2) + t100 * qJ(3);
t114 = -t78 - t126;
t113 = t102 * t78 + t119;
t112 = V_base(4) * t127 + qJD(2) + V_base(3);
t111 = V_base(5) * qJ(2) + t118;
t110 = V_base(4) * t75 + t112;
t109 = qJD(3) * t100 + t111;
t104 = Icges(2,4) * t108;
t99 = -qJD(4) + t102;
t91 = t108 * rSges(2,1) - t107 * rSges(2,2);
t90 = t107 * rSges(2,1) + t108 * rSges(2,2);
t89 = Icges(2,1) * t108 - t122;
t88 = Icges(2,1) * t107 + t104;
t87 = -Icges(2,2) * t107 + t104;
t86 = Icges(2,2) * t108 + t122;
t83 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t82 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t81 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t80 = t101 * rSges(3,1) - t100 * rSges(3,2);
t79 = t101 * rSges(4,1) + t100 * rSges(4,3);
t77 = t100 * rSges(3,1) + t101 * rSges(3,2);
t76 = t100 * rSges(4,1) - t101 * rSges(4,3);
t62 = t100 * t129 - t101 * t128;
t59 = Icges(5,4) * t62;
t57 = V_base(5) * rSges(2,3) - t102 * t90 + t118;
t56 = t102 * t91 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t55 = V_base(4) * t90 - V_base(5) * t91 + V_base(3);
t54 = -t61 * rSges(5,1) + t62 * rSges(5,2);
t53 = t62 * rSges(5,1) + t61 * rSges(5,2);
t52 = -Icges(5,1) * t61 + t59;
t51 = Icges(5,1) * t62 + t123;
t50 = Icges(5,2) * t62 - t123;
t49 = Icges(5,2) * t61 + t59;
t46 = V_base(5) * rSges(3,3) + (-t77 - t127) * t102 + t111;
t45 = t102 * t80 + (-rSges(3,3) + t124) * V_base(4) + t119;
t44 = V_base(4) * t77 + (-t80 - t126) * V_base(5) + t112;
t43 = V_base(5) * rSges(4,2) + (t115 - t76) * t102 + t109;
t42 = -qJD(3) * t101 + t102 * t79 + (-rSges(4,2) + t124) * V_base(4) + t113;
t41 = V_base(4) * t76 + (t114 - t79) * V_base(5) + t110;
t40 = -t99 * t53 + t130 * V_base(5) + (t115 - t125) * t102 + t109;
t39 = t99 * t54 + (pkin(3) * t102 - qJD(3)) * t101 + (t124 - t130) * V_base(4) + t113;
t38 = (t53 + t125) * V_base(4) + (-pkin(3) * t101 + t114 - t54) * V_base(5) + t110;
t1 = m(1) * (t81 ^ 2 + t82 ^ 2 + t83 ^ 2) / 0.2e1 + m(2) * (t55 ^ 2 + t56 ^ 2 + t57 ^ 2) / 0.2e1 + m(3) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + m(4) * (t41 ^ 2 + t42 ^ 2 + t43 ^ 2) / 0.2e1 + m(5) * (t38 ^ 2 + t39 ^ 2 + t40 ^ 2) / 0.2e1 + Icges(1,3) * V_base(6) ^ 2 / 0.2e1 + Icges(5,3) * t99 ^ 2 / 0.2e1 + (Icges(2,3) / 0.2e1 + Icges(3,3) / 0.2e1 + Icges(4,2) / 0.2e1) * t102 ^ 2 + (Icges(1,6) * V_base(6) + (-Icges(5,5) * t62 - Icges(5,6) * t61) * t99 + (Icges(2,5) * t107 + Icges(2,6) * t108 + t132 * t100 + t131 * t101) * t102 + Icges(1,2) * t133) * V_base(5) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + (Icges(5,5) * t61 - Icges(5,6) * t62) * t99 + (Icges(2,5) * t108 - Icges(2,6) * t107 - t131 * t100 + t132 * t101) * t102 + Icges(1,1) * t134) * V_base(4) + ((t138 * t100 + t136 * t101 - t107 * t86 + t108 * t88 + t62 * t49 - t61 * t51) * V_base(5) + (t137 * t100 + t135 * t101 - t107 * t87 + t108 * t89 + t62 * t50 - t61 * t52) * V_base(4)) * t134 + ((t136 * t100 - t138 * t101 + t107 * t88 + t108 * t86 + t61 * t49 + t62 * t51) * V_base(5) + (t135 * t100 - t137 * t101 + t107 * t89 + t108 * t87 + t61 * t50 + t62 * t52) * V_base(4)) * t133;
T  = t1;
