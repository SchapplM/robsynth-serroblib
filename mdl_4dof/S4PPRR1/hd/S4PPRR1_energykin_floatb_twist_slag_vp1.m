% Calculate kinetic energy for
% S4PPRR1
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
%   pkin=[a2,a3,a4,d3,d4,theta1]';
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
% Datum: 2019-03-08 18:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPRR1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PPRR1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR1_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR1_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PPRR1_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PPRR1_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:15:27
% EndTime: 2019-03-08 18:15:28
% DurationCPUTime: 0.85s
% Computational Cost: add. (421->153), mult. (546->185), div. (0->0), fcn. (458->6), ass. (0->81)
t152 = Icges(2,4) - Icges(3,5);
t151 = Icges(2,1) + Icges(3,1);
t150 = Icges(2,2) + Icges(3,3);
t112 = sin(pkin(6));
t149 = t152 * t112;
t113 = cos(pkin(6));
t148 = t152 * t113;
t147 = -t150 * t113 - t149;
t146 = t150 * t112 - t148;
t145 = t151 * t112 + t148;
t144 = t151 * t113 - t149;
t143 = V_base(4) / 0.2e1;
t142 = V_base(5) / 0.2e1;
t141 = Icges(3,4) + Icges(2,5);
t140 = Icges(2,6) - Icges(3,6);
t139 = rSges(5,3) + pkin(5);
t138 = cos(qJ(3));
t137 = sin(qJ(3));
t136 = t112 * pkin(2);
t135 = t113 * pkin(2);
t134 = t138 * pkin(3);
t124 = t112 * t137;
t75 = -t113 * t138 - t124;
t133 = Icges(4,4) * t75;
t129 = qJ(3) + qJ(4);
t120 = sin(t129);
t121 = cos(t129);
t70 = -t112 * t120 - t113 * t121;
t132 = Icges(5,4) * t70;
t128 = V_base(5) * qJ(1) + V_base(1);
t125 = qJD(1) + V_base(3);
t123 = t113 * t137;
t95 = pkin(1) * t113 + t112 * qJ(2);
t122 = -t95 - t135;
t109 = V_base(6) - qJD(3);
t92 = t112 * pkin(1) - qJ(2) * t113;
t119 = V_base(4) * t92 + t125;
t118 = qJD(2) * t112 + t128;
t117 = V_base(4) * t136 + t119;
t116 = -qJD(2) * t113 + V_base(6) * t95 + V_base(2);
t115 = V_base(4) * pkin(4) + V_base(6) * t135 + t116;
t114 = (-t92 - t136) * V_base(6) + t118;
t104 = -qJD(4) + t109;
t97 = rSges(2,1) * t113 - t112 * rSges(2,2);
t96 = rSges(3,1) * t113 + t112 * rSges(3,3);
t94 = t112 * rSges(2,1) + rSges(2,2) * t113;
t93 = t112 * rSges(3,1) - rSges(3,3) * t113;
t79 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t78 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t77 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t76 = t112 * t138 - t123;
t72 = Icges(4,4) * t76;
t71 = t112 * t121 - t113 * t120;
t69 = Icges(5,4) * t71;
t68 = pkin(3) * t124 + t134 * t113;
t67 = -pkin(3) * t123 + t134 * t112;
t66 = V_base(5) * rSges(2,3) - t94 * V_base(6) + t128;
t65 = V_base(6) * t97 + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t64 = -rSges(4,1) * t75 + rSges(4,2) * t76;
t63 = rSges(4,1) * t76 + rSges(4,2) * t75;
t62 = t94 * V_base(4) - t97 * V_base(5) + t125;
t61 = -Icges(4,1) * t75 + t72;
t60 = Icges(4,1) * t76 + t133;
t59 = Icges(4,2) * t76 - t133;
t58 = Icges(4,2) * t75 + t72;
t55 = -rSges(5,1) * t70 + rSges(5,2) * t71;
t54 = rSges(5,1) * t71 + rSges(5,2) * t70;
t53 = -Icges(5,1) * t70 + t69;
t52 = Icges(5,1) * t71 + t132;
t51 = Icges(5,2) * t71 - t132;
t50 = Icges(5,2) * t70 + t69;
t47 = V_base(5) * rSges(3,2) + (-t92 - t93) * V_base(6) + t118;
t46 = V_base(6) * t96 + (-rSges(3,2) - qJ(1)) * V_base(4) + t116;
t45 = V_base(4) * t93 + (-t95 - t96) * V_base(5) + t119;
t44 = -t109 * t63 + (-pkin(4) - rSges(4,3)) * V_base(5) + t114;
t43 = t109 * t64 + (rSges(4,3) - qJ(1)) * V_base(4) + t115;
t42 = V_base(4) * t63 + (t122 - t64) * V_base(5) + t117;
t41 = -t104 * t54 - t109 * t67 + (-pkin(4) - t139) * V_base(5) + t114;
t40 = t104 * t55 + t109 * t68 + (-qJ(1) + t139) * V_base(4) + t115;
t39 = (t54 + t67) * V_base(4) + (t122 - t55 - t68) * V_base(5) + t117;
t1 = m(1) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(2) * (t62 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + m(3) * (t45 ^ 2 + t46 ^ 2 + t47 ^ 2) / 0.2e1 + m(4) * (t42 ^ 2 + t43 ^ 2 + t44 ^ 2) / 0.2e1 + m(5) * (t39 ^ 2 + t40 ^ 2 + t41 ^ 2) / 0.2e1 + Icges(4,3) * t109 ^ 2 / 0.2e1 + Icges(5,3) * t104 ^ 2 / 0.2e1 + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1 + Icges(3,2) / 0.2e1) * V_base(6) ^ 2 + ((-Icges(4,5) * t76 - Icges(4,6) * t75) * t109 + (-Icges(5,5) * t71 - Icges(5,6) * t70) * t104 + (t141 * t112 + t140 * t113 + Icges(1,6)) * V_base(6) + Icges(1,2) * t142) * V_base(5) + (Icges(1,4) * V_base(5) + (Icges(4,5) * t75 - Icges(4,6) * t76) * t109 + (Icges(5,5) * t70 - Icges(5,6) * t71) * t104 + (-t140 * t112 + t141 * t113 + Icges(1,5)) * V_base(6) + Icges(1,1) * t143) * V_base(4) + ((t147 * t112 + t145 * t113 + t50 * t71 - t52 * t70 + t58 * t76 - t60 * t75) * V_base(5) + (t146 * t112 + t144 * t113 + t71 * t51 - t70 * t53 + t76 * t59 - t75 * t61) * V_base(4)) * t143 + ((t145 * t112 - t147 * t113 + t70 * t50 + t71 * t52 + t75 * t58 + t76 * t60) * V_base(5) + (t144 * t112 - t146 * t113 + t51 * t70 + t53 * t71 + t59 * t75 + t61 * t76) * V_base(4)) * t142;
T  = t1;
