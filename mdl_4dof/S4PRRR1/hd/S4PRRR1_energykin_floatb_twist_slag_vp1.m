% Calculate kinetic energy for
% S4PRRR1
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
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_energykin_floatb_twist_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR1_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRR1_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRR1_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:25:03
% EndTime: 2019-03-08 18:25:04
% DurationCPUTime: 0.58s
% Computational Cost: add. (517->157), mult. (358->198), div. (0->0), fcn. (184->8), ass. (0->79)
t136 = V_base(4) / 0.2e1;
t135 = pkin(6) + rSges(5,3);
t113 = cos(pkin(7));
t134 = pkin(1) * t113;
t107 = V_base(6) + qJD(2);
t133 = pkin(2) * t107;
t101 = qJD(3) + t107;
t132 = pkin(3) * t101;
t131 = -pkin(4) - qJ(1);
t111 = pkin(7) + qJ(2);
t106 = qJ(3) + t111;
t99 = sin(t106);
t130 = Icges(4,4) * t99;
t102 = qJ(4) + t106;
t94 = sin(t102);
t129 = Icges(5,4) * t94;
t112 = sin(pkin(7));
t128 = Icges(2,4) * t112;
t103 = sin(t111);
t127 = Icges(3,4) * t103;
t120 = pkin(1) * V_base(6);
t126 = t113 * t120 + V_base(2);
t125 = -pkin(5) + t131;
t124 = V_base(5) * qJ(1) + V_base(1);
t121 = qJD(1) + V_base(3);
t104 = cos(t111);
t119 = t104 * t133 + t126;
t118 = V_base(4) * t112 * pkin(1) + t121;
t117 = V_base(4) * pkin(2) * t103 + t118;
t116 = -pkin(2) * t104 - t134;
t115 = V_base(5) * pkin(4) - t112 * t120 + t124;
t114 = V_base(5) * pkin(5) - t103 * t133 + t115;
t105 = Icges(2,4) * t113;
t100 = cos(t106);
t98 = Icges(3,4) * t104;
t95 = cos(t102);
t93 = qJD(4) + t101;
t92 = Icges(4,4) * t100;
t90 = Icges(5,4) * t95;
t89 = rSges(2,1) * t113 - t112 * rSges(2,2);
t88 = t112 * rSges(2,1) + rSges(2,2) * t113;
t87 = Icges(2,1) * t113 - t128;
t86 = Icges(2,1) * t112 + t105;
t85 = -Icges(2,2) * t112 + t105;
t84 = Icges(2,2) * t113 + t128;
t80 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t79 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t78 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t77 = rSges(3,1) * t104 - rSges(3,2) * t103;
t76 = rSges(3,1) * t103 + rSges(3,2) * t104;
t75 = Icges(3,1) * t104 - t127;
t74 = Icges(3,1) * t103 + t98;
t73 = -Icges(3,2) * t103 + t98;
t72 = Icges(3,2) * t104 + t127;
t69 = rSges(4,1) * t100 - rSges(4,2) * t99;
t68 = rSges(4,1) * t99 + rSges(4,2) * t100;
t67 = Icges(4,1) * t100 - t130;
t66 = Icges(4,1) * t99 + t92;
t65 = -Icges(4,2) * t99 + t92;
t64 = Icges(4,2) * t100 + t130;
t61 = rSges(5,1) * t95 - rSges(5,2) * t94;
t60 = rSges(5,1) * t94 + rSges(5,2) * t95;
t59 = Icges(5,1) * t95 - t129;
t58 = Icges(5,1) * t94 + t90;
t57 = -Icges(5,2) * t94 + t90;
t56 = Icges(5,2) * t95 + t129;
t53 = V_base(5) * rSges(2,3) - t88 * V_base(6) + t124;
t52 = t89 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t51 = t88 * V_base(4) - t89 * V_base(5) + t121;
t50 = V_base(5) * rSges(3,3) - t107 * t76 + t115;
t49 = t107 * t77 + (-rSges(3,3) + t131) * V_base(4) + t126;
t48 = V_base(4) * t76 + (-t77 - t134) * V_base(5) + t118;
t47 = V_base(5) * rSges(4,3) - t101 * t68 + t114;
t46 = t101 * t69 + (-rSges(4,3) + t125) * V_base(4) + t119;
t45 = V_base(4) * t68 + (t116 - t69) * V_base(5) + t117;
t44 = -t99 * t132 + t135 * V_base(5) - t60 * t93 + t114;
t43 = t100 * t132 + t61 * t93 + (t125 - t135) * V_base(4) + t119;
t42 = (pkin(3) * t99 + t60) * V_base(4) + (-pkin(3) * t100 + t116 - t61) * V_base(5) + t117;
t1 = m(1) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + Icges(1,2) * V_base(5) ^ 2 / 0.2e1 + m(2) * (t51 ^ 2 + t52 ^ 2 + t53 ^ 2) / 0.2e1 + m(3) * (t48 ^ 2 + t49 ^ 2 + t50 ^ 2) / 0.2e1 + m(4) * (t45 ^ 2 + t46 ^ 2 + t47 ^ 2) / 0.2e1 + m(5) * (t42 ^ 2 + t43 ^ 2 + t44 ^ 2) / 0.2e1 + ((Icges(5,5) * t94 + Icges(5,6) * t95) * V_base(5) + Icges(5,3) * t93 / 0.2e1) * t93 + ((Icges(3,5) * t103 + Icges(3,6) * t104) * V_base(5) + Icges(3,3) * t107 / 0.2e1) * t107 + ((Icges(4,5) * t99 + Icges(4,6) * t100) * V_base(5) + Icges(4,3) * t101 / 0.2e1) * t101 + ((Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6) + (Icges(2,5) * t112 + Icges(2,6) * t113 + Icges(1,6)) * V_base(5)) * V_base(6) + (Icges(1,4) * V_base(5) + (Icges(5,5) * t95 - Icges(5,6) * t94) * t93 + (Icges(3,5) * t104 - Icges(3,6) * t103) * t107 + (Icges(4,5) * t100 - Icges(4,6) * t99) * t101 + (Icges(2,5) * t113 - Icges(2,6) * t112 + Icges(1,5)) * V_base(6) + Icges(1,1) * t136) * V_base(4) + ((t100 * t66 - t103 * t72 + t104 * t74 - t112 * t84 + t113 * t86 - t94 * t56 + t58 * t95 - t64 * t99) * V_base(5) + (t100 * t67 - t103 * t73 + t104 * t75 - t112 * t85 + t113 * t87 - t94 * t57 + t59 * t95 - t65 * t99) * V_base(4)) * t136 + ((t100 * t64 + t103 * t74 + t104 * t72 + t112 * t86 + t113 * t84 + t56 * t95 + t94 * t58 + t66 * t99) * V_base(5) + (t100 * t65 + t103 * t75 + t104 * t73 + t112 * t87 + t113 * t85 + t57 * t95 + t59 * t94 + t67 * t99) * V_base(4)) * V_base(5) / 0.2e1;
T  = t1;
