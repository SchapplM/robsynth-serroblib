% Calculate kinetic energy for
% S5PRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% Icges [6x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:05
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR2_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR2_energykin_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR2_energykin_floatb_twist_slag_vp1: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PRRRR2_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5PRRRR2_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR2_energykin_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRR2_energykin_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5PRRRR2_energykin_floatb_twist_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:04:41
% EndTime: 2019-12-05 17:04:42
% DurationCPUTime: 0.92s
% Computational Cost: add. (691->191), mult. (520->242), div. (0->0), fcn. (312->8), ass. (0->95)
t162 = Icges(1,4) + Icges(2,4);
t123 = sin(qJ(2));
t158 = pkin(2) * t123;
t125 = cos(qJ(2));
t157 = pkin(2) * t125;
t121 = qJ(2) + qJ(3);
t113 = sin(t121);
t156 = pkin(3) * t113;
t114 = cos(t121);
t155 = pkin(3) * t114;
t117 = qJ(4) + t121;
t110 = sin(t117);
t154 = pkin(6) * t110;
t111 = cos(t117);
t153 = pkin(6) * t111;
t152 = -pkin(4) - qJ(1);
t151 = Icges(3,4) * t123;
t150 = Icges(4,4) * t113;
t149 = Icges(5,4) * t110;
t122 = sin(qJ(5));
t148 = Icges(6,4) * t122;
t124 = cos(qJ(5));
t147 = Icges(6,4) * t124;
t146 = -pkin(5) + t152;
t145 = V_base(5) * qJ(1) + V_base(1);
t144 = V_base(6) * pkin(1) + V_base(2);
t140 = qJD(1) + V_base(3);
t139 = -pkin(1) - t157;
t112 = V_base(6) + qJD(2);
t138 = t112 * t157 + t144;
t137 = V_base(4) * t158 + t140;
t109 = qJD(3) + t112;
t136 = t109 * t155 + t138;
t135 = V_base(4) * t156 + t137;
t134 = rSges(6,1) * t124 - rSges(6,2) * t122;
t133 = Icges(6,1) * t124 - t148;
t132 = -Icges(6,2) * t122 + t147;
t131 = Icges(6,5) * t124 - Icges(6,6) * t122;
t130 = t139 - t155;
t129 = V_base(5) * pkin(4) - t112 * t158 + t145;
t106 = qJD(4) + t109;
t89 = -qJD(5) * t111 + V_base(5);
t90 = qJD(5) * t110 + V_base(4);
t128 = t106 * (Icges(6,5) * t122 + Icges(6,6) * t124) + (-Icges(6,3) * t111 + t131 * t110) * t89 + (Icges(6,3) * t110 + t131 * t111) * t90;
t127 = V_base(5) * pkin(5) - t109 * t156 + t129;
t59 = -Icges(6,6) * t111 + t132 * t110;
t60 = Icges(6,6) * t110 + t132 * t111;
t61 = -Icges(6,5) * t111 + t133 * t110;
t62 = Icges(6,5) * t110 + t133 * t111;
t94 = Icges(6,2) * t124 + t148;
t97 = Icges(6,1) * t122 + t147;
t126 = (-t122 * t60 + t124 * t62) * t90 + (-t122 * t59 + t124 * t61) * t89 + (-t122 * t94 + t124 * t97) * t106;
t116 = Icges(3,4) * t125;
t108 = Icges(4,4) * t114;
t105 = Icges(5,4) * t111;
t102 = rSges(3,1) * t125 - t123 * rSges(3,2);
t101 = t123 * rSges(3,1) + rSges(3,2) * t125;
t100 = rSges(6,1) * t122 + rSges(6,2) * t124;
t99 = Icges(3,1) * t125 - t151;
t98 = Icges(3,1) * t123 + t116;
t96 = -Icges(3,2) * t123 + t116;
t95 = Icges(3,2) * t125 + t151;
t88 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t87 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t86 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t84 = -V_base(5) * rSges(2,1) + V_base(4) * rSges(2,2) + t140;
t83 = rSges(4,1) * t114 - rSges(4,2) * t113;
t82 = rSges(4,1) * t113 + rSges(4,2) * t114;
t81 = Icges(4,1) * t114 - t150;
t80 = Icges(4,1) * t113 + t108;
t79 = -Icges(4,2) * t113 + t108;
t78 = Icges(4,2) * t114 + t150;
t75 = rSges(5,1) * t111 - rSges(5,2) * t110;
t74 = rSges(5,1) * t110 + rSges(5,2) * t111;
t73 = V_base(6) * rSges(2,1) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t72 = -V_base(6) * rSges(2,2) + V_base(5) * rSges(2,3) + t145;
t71 = Icges(5,1) * t111 - t149;
t70 = Icges(5,1) * t110 + t105;
t69 = -Icges(5,2) * t110 + t105;
t68 = Icges(5,2) * t111 + t149;
t65 = V_base(5) * rSges(3,3) - t101 * t112 + t145;
t64 = rSges(6,3) * t110 + t134 * t111;
t63 = -rSges(6,3) * t111 + t134 * t110;
t56 = t102 * t112 + (-rSges(3,3) - qJ(1)) * V_base(4) + t144;
t55 = t101 * V_base(4) + (-pkin(1) - t102) * V_base(5) + t140;
t54 = V_base(5) * rSges(4,3) - t109 * t82 + t129;
t53 = t109 * t83 + (-rSges(4,3) + t152) * V_base(4) + t138;
t52 = V_base(4) * t82 + (t139 - t83) * V_base(5) + t137;
t51 = V_base(5) * rSges(5,3) - t106 * t74 + t127;
t50 = t106 * t75 + (-rSges(5,3) + t146) * V_base(4) + t136;
t49 = V_base(4) * t74 + (t130 - t75) * V_base(5) + t135;
t48 = t100 * t89 + (-t63 + t153) * t106 + t127;
t47 = -t100 * t90 + (t64 + t154) * t106 + t146 * V_base(4) + t136;
t46 = -V_base(4) * t153 + t90 * t63 - t89 * t64 + (t130 - t154) * V_base(5) + t135;
t1 = m(1) * (t86 ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(2) * (t72 ^ 2 + t73 ^ 2 + t84 ^ 2) / 0.2e1 + m(3) * (t55 ^ 2 + t56 ^ 2 + t65 ^ 2) / 0.2e1 + m(4) * (t52 ^ 2 + t53 ^ 2 + t54 ^ 2) / 0.2e1 + m(5) * (t49 ^ 2 + t50 ^ 2 + t51 ^ 2) / 0.2e1 + m(6) * (t46 ^ 2 + t47 ^ 2 + t48 ^ 2) / 0.2e1 + t90 * (t128 * t110 + t126 * t111) / 0.2e1 + t89 * (t126 * t110 - t128 * t111) / 0.2e1 + ((t122 * t62 + t124 * t60) * t90 + (t122 * t61 + t124 * t59) * t89 + (t122 * t97 + t124 * t94 + Icges(5,3)) * t106) * t106 / 0.2e1 + ((-t110 * t68 + t111 * t70 - t113 * t78 + t114 * t80 - t123 * t95 + t125 * t98 + t162) * V_base(5) + (-t110 * t69 + t111 * t71 - t113 * t79 + t114 * t81 - t123 * t96 + t125 * t99 + Icges(1,1) + Icges(2,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t110 * t70 + t111 * t68 + t113 * t80 + t114 * t78 + t123 * t98 + t125 * t95 + Icges(1,2) + Icges(2,2)) * V_base(5) + (t110 * t71 + t111 * t69 + t113 * t81 + t114 * t79 + t123 * t99 + t125 * t96 + t162) * V_base(4)) * V_base(5) / 0.2e1 + V_base(4) * t106 * (Icges(5,5) * t111 - Icges(5,6) * t110) + V_base(5) * t106 * (Icges(5,5) * t110 + Icges(5,6) * t111) + ((Icges(1,6) + Icges(2,6)) * V_base(5) + (Icges(1,5) + Icges(2,5)) * V_base(4) + (Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + ((Icges(3,5) * t123 + Icges(3,6) * t125) * V_base(5) + (Icges(3,5) * t125 - Icges(3,6) * t123) * V_base(4) + Icges(3,3) * t112 / 0.2e1) * t112 + ((Icges(4,5) * t113 + Icges(4,6) * t114) * V_base(5) + (Icges(4,5) * t114 - Icges(4,6) * t113) * V_base(4) + Icges(4,3) * t109 / 0.2e1) * t109;
T = t1;
