% Calculate kinetic energy for
% S4RPPR4
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR4_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR4_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR4_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPPR4_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR4_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR4_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPPR4_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPPR4_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:38:47
% EndTime: 2019-12-31 16:38:48
% DurationCPUTime: 0.93s
% Computational Cost: add. (511->165), mult. (494->205), div. (0->0), fcn. (324->6), ass. (0->81)
t153 = Icges(3,4) + Icges(4,6);
t152 = Icges(3,1) + Icges(4,2);
t151 = -Icges(4,4) + Icges(3,5);
t150 = Icges(4,5) - Icges(3,6);
t149 = Icges(3,2) + Icges(4,3);
t105 = qJ(1) + pkin(6);
t100 = cos(t105);
t148 = t153 * t100;
t99 = sin(t105);
t147 = t153 * t99;
t146 = -t149 * t100 - t147;
t145 = t149 * t99 - t148;
t144 = t152 * t99 + t148;
t143 = t152 * t100 - t147;
t101 = V_base(6) + qJD(1);
t106 = sin(qJ(4));
t108 = cos(qJ(4));
t131 = Icges(5,4) * t106;
t112 = Icges(5,2) * t108 + t131;
t50 = Icges(5,6) * t100 + t112 * t99;
t51 = Icges(5,6) * t99 - t112 * t100;
t130 = Icges(5,4) * t108;
t113 = Icges(5,1) * t106 + t130;
t52 = Icges(5,5) * t100 + t113 * t99;
t53 = Icges(5,5) * t99 - t113 * t100;
t81 = qJD(4) * t99 + V_base(5);
t82 = qJD(4) * t100 + V_base(4);
t86 = -Icges(5,2) * t106 + t130;
t89 = Icges(5,1) * t108 - t131;
t140 = (t106 * t52 + t108 * t50) * t82 + (t106 * t53 + t108 * t51) * t81 + (t106 * t89 + t108 * t86) * t101;
t139 = pkin(5) * t99;
t107 = sin(qJ(1));
t137 = pkin(1) * t107;
t109 = cos(qJ(1));
t136 = pkin(1) * t109;
t135 = -pkin(4) - qJ(2);
t132 = Icges(2,4) * t107;
t128 = t101 * t136 + V_base(2);
t127 = V_base(5) * pkin(4) + V_base(1);
t72 = pkin(2) * t99 - qJ(3) * t100;
t124 = -t72 - t137;
t75 = pkin(2) * t100 + qJ(3) * t99;
t123 = -t75 - t136;
t122 = t101 * t75 + t128;
t121 = V_base(4) * t137 + qJD(2) + V_base(3);
t120 = V_base(5) * qJ(2) + t127;
t119 = V_base(4) * t72 + t121;
t118 = qJD(3) * t99 + t120;
t117 = rSges(5,1) * t106 + rSges(5,2) * t108;
t111 = Icges(5,5) * t106 + Icges(5,6) * t108;
t110 = t101 * (Icges(5,5) * t108 - Icges(5,6) * t106) + (Icges(5,3) * t100 + t111 * t99) * t82 + (Icges(5,3) * t99 - t111 * t100) * t81;
t103 = Icges(2,4) * t109;
t94 = rSges(2,1) * t109 - t107 * rSges(2,2);
t93 = rSges(5,1) * t108 - rSges(5,2) * t106;
t92 = t107 * rSges(2,1) + rSges(2,2) * t109;
t91 = Icges(2,1) * t109 - t132;
t90 = Icges(2,1) * t107 + t103;
t88 = -Icges(2,2) * t107 + t103;
t87 = Icges(2,2) * t109 + t132;
t80 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t79 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t78 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t77 = rSges(3,1) * t100 - rSges(3,2) * t99;
t76 = -rSges(4,2) * t100 + rSges(4,3) * t99;
t74 = rSges(3,1) * t99 + rSges(3,2) * t100;
t73 = -rSges(4,2) * t99 - rSges(4,3) * t100;
t57 = rSges(5,3) * t99 - t117 * t100;
t56 = rSges(5,3) * t100 + t117 * t99;
t55 = V_base(5) * rSges(2,3) - t101 * t92 + t127;
t54 = t101 * t94 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t47 = t92 * V_base(4) - t94 * V_base(5) + V_base(3);
t46 = V_base(5) * rSges(3,3) + (-t74 - t137) * t101 + t120;
t45 = t101 * t77 + (-rSges(3,3) + t135) * V_base(4) + t128;
t44 = V_base(4) * t74 + (-t77 - t136) * V_base(5) + t121;
t43 = V_base(5) * rSges(4,1) + (t124 - t73) * t101 + t118;
t42 = -qJD(3) * t100 + t101 * t76 + (-rSges(4,1) + t135) * V_base(4) + t122;
t41 = V_base(4) * t73 + (t123 - t76) * V_base(5) + t119;
t40 = V_base(5) * pkin(3) + t81 * t93 + (t124 - t57 - t139) * t101 + t118;
t39 = t101 * t56 - t82 * t93 + (pkin(5) * t101 - qJD(3)) * t100 + (-pkin(3) + t135) * V_base(4) + t122;
t38 = V_base(4) * t139 - t81 * t56 + t82 * t57 + (-pkin(5) * t100 + t123) * V_base(5) + t119;
t1 = m(1) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + m(2) * (t47 ^ 2 + t54 ^ 2 + t55 ^ 2) / 0.2e1 + m(3) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + m(4) * (t41 ^ 2 + t42 ^ 2 + t43 ^ 2) / 0.2e1 + m(5) * (t38 ^ 2 + t39 ^ 2 + t40 ^ 2) / 0.2e1 + t82 * (t110 * t100 + t140 * t99) / 0.2e1 + t81 * (-t140 * t100 + t110 * t99) / 0.2e1 + ((-t106 * t50 + t108 * t52) * t82 + (-t106 * t51 + t108 * t53) * t81 + (-t106 * t86 + t108 * t89 + Icges(4,1) + Icges(2,3) + Icges(3,3)) * t101) * t101 / 0.2e1 + ((t144 * t100 - t107 * t87 + t109 * t90 + t146 * t99 + Icges(1,4)) * V_base(5) + (t143 * t100 - t107 * t88 + t109 * t91 + t145 * t99 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((-t146 * t100 + t107 * t90 + t109 * t87 + t144 * t99 + Icges(1,2)) * V_base(5) + (-t145 * t100 + t107 * t91 + t109 * t88 + t143 * t99 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t101 * (Icges(2,5) * t107 + Icges(2,6) * t109 - t150 * t100 + t151 * t99) + V_base(4) * t101 * (Icges(2,5) * t109 - Icges(2,6) * t107 + t151 * t100 + t150 * t99) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6);
T = t1;
