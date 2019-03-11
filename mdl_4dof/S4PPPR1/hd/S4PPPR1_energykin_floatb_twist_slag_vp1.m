% Calculate kinetic energy for
% S4PPPR1
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
%   pkin=[a2,a3,a4,d4,theta1]';
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
% Datum: 2019-03-08 18:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPPR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPPR1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPPR1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PPPR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPPR1_energykin_floatb_twist_slag_vp1: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPPR1_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPPR1_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPPR1_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:08:52
% EndTime: 2019-03-08 18:08:53
% DurationCPUTime: 0.65s
% Computational Cost: add. (279->133), mult. (454->157), div. (0->0), fcn. (326->4), ass. (0->68)
t137 = Icges(2,4) - Icges(4,5) + Icges(3,6);
t136 = Icges(2,1) + Icges(3,2) + Icges(4,3);
t135 = Icges(4,1) + Icges(2,2) + Icges(3,3);
t99 = cos(pkin(5));
t134 = t137 * t99;
t98 = sin(pkin(5));
t133 = t137 * t98;
t132 = -t135 * t99 - t133;
t131 = t135 * t98 - t134;
t130 = t136 * t98 + t134;
t129 = t136 * t99 - t133;
t128 = V_base(4) / 0.2e1;
t127 = V_base(5) / 0.2e1;
t126 = -Icges(3,4) + Icges(2,5) - Icges(4,6);
t125 = Icges(4,4) - Icges(3,5) + Icges(2,6);
t124 = pkin(3) * t98;
t123 = pkin(3) * t99;
t122 = pkin(4) + rSges(5,3);
t121 = sin(qJ(4));
t120 = -pkin(2) - qJ(1);
t100 = cos(qJ(4));
t56 = t100 * t99 - t121 * t98;
t118 = Icges(5,4) * t56;
t114 = qJ(3) * t98;
t113 = qJ(3) * t99;
t112 = V_base(5) * qJ(1) + V_base(1);
t109 = qJD(1) + V_base(3);
t108 = qJD(2) * t98 + t112;
t79 = pkin(1) * t98 - qJ(2) * t99;
t107 = -t79 - t114;
t83 = pkin(1) * t99 + qJ(2) * t98;
t106 = -t83 - t113;
t105 = V_base(4) * t79 + t109;
t104 = V_base(4) * t114 + t105;
t103 = V_base(5) * pkin(2) + qJD(3) * t99 + t108;
t102 = -qJD(2) * t99 + V_base(6) * t83 + V_base(2);
t101 = qJD(3) * t98 + V_base(6) * t113 + t102;
t95 = V_base(6) + qJD(4);
t86 = rSges(2,1) * t99 - rSges(2,2) * t98;
t85 = -rSges(4,1) * t99 + rSges(4,3) * t98;
t84 = -rSges(3,2) * t99 + rSges(3,3) * t98;
t82 = rSges(2,1) * t98 + rSges(2,2) * t99;
t81 = rSges(4,1) * t98 + rSges(4,3) * t99;
t80 = -rSges(3,2) * t98 - rSges(3,3) * t99;
t60 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t59 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t58 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t57 = t98 * t100 + t121 * t99;
t53 = Icges(5,4) * t57;
t52 = V_base(5) * rSges(2,3) - t82 * V_base(6) + t112;
t51 = t86 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t50 = -rSges(5,1) * t56 + rSges(5,2) * t57;
t49 = rSges(5,1) * t57 + rSges(5,2) * t56;
t48 = t82 * V_base(4) - t86 * V_base(5) + t109;
t47 = -Icges(5,1) * t56 + t53;
t46 = Icges(5,1) * t57 + t118;
t45 = Icges(5,2) * t57 - t118;
t44 = Icges(5,2) * t56 + t53;
t41 = V_base(5) * rSges(3,1) + (-t79 - t80) * V_base(6) + t108;
t40 = t84 * V_base(6) + (-rSges(3,1) - qJ(1)) * V_base(4) + t102;
t39 = t80 * V_base(4) + (-t83 - t84) * V_base(5) + t105;
t38 = -V_base(5) * rSges(4,2) + (t107 - t85) * V_base(6) + t103;
t37 = t81 * V_base(6) + (rSges(4,2) + t120) * V_base(4) + t101;
t36 = t85 * V_base(4) + (t106 - t81) * V_base(5) + t104;
t35 = -t50 * t95 + t122 * V_base(5) + (t107 + t123) * V_base(6) + t103;
t34 = V_base(6) * t124 + t49 * t95 + (t120 - t122) * V_base(4) + t101;
t33 = (t50 - t123) * V_base(4) + (t106 - t49 - t124) * V_base(5) + t104;
t1 = m(1) * (t58 ^ 2 + t59 ^ 2 + t60 ^ 2) / 0.2e1 + m(2) * (t48 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(3) * (t39 ^ 2 + t40 ^ 2 + t41 ^ 2) / 0.2e1 + m(4) * (t36 ^ 2 + t37 ^ 2 + t38 ^ 2) / 0.2e1 + m(5) * (t33 ^ 2 + t34 ^ 2 + t35 ^ 2) / 0.2e1 + Icges(5,3) * t95 ^ 2 / 0.2e1 + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,1) / 0.2e1 + Icges(4,2) / 0.2e1) * V_base(6) ^ 2 + ((-Icges(5,5) * t56 + Icges(5,6) * t57) * t95 + (t125 * t99 + t126 * t98 + Icges(1,6)) * V_base(6) + Icges(1,2) * t127) * V_base(5) + (Icges(1,4) * V_base(5) + (Icges(5,5) * t57 + Icges(5,6) * t56) * t95 + (-t125 * t98 + t126 * t99 + Icges(1,5)) * V_base(6) + Icges(1,1) * t128) * V_base(4) + ((t130 * t99 + t132 * t98 + t45 * t56 + t47 * t57) * V_base(5) + (t129 * t99 + t131 * t98 + t44 * t56 + t46 * t57) * V_base(4)) * t128 + ((t130 * t98 - t132 * t99 + t57 * t45 - t47 * t56) * V_base(5) + (t129 * t98 - t131 * t99 + t44 * t57 - t46 * t56) * V_base(4)) * t127;
T  = t1;
