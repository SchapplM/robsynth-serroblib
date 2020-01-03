% Calculate kinetic energy for
% S4RRPR5
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2019-12-31 17:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPR5_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR5_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR5_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRPR5_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPR5_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR5_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPR5_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPR5_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:03:26
% EndTime: 2019-12-31 17:03:27
% DurationCPUTime: 1.00s
% Computational Cost: add. (532->164), mult. (494->211), div. (0->0), fcn. (324->6), ass. (0->81)
t153 = Icges(3,4) + Icges(4,6);
t152 = Icges(3,1) + Icges(4,2);
t151 = -Icges(4,4) + Icges(3,5);
t150 = Icges(4,5) - Icges(3,6);
t149 = Icges(3,2) + Icges(4,3);
t106 = qJ(1) + qJ(2);
t102 = cos(t106);
t148 = t153 * t102;
t101 = sin(t106);
t147 = t153 * t101;
t146 = -t102 * t149 - t147;
t145 = t101 * t149 - t148;
t144 = t101 * t152 + t148;
t143 = t102 * t152 - t147;
t107 = sin(qJ(4));
t109 = cos(qJ(4));
t132 = Icges(5,4) * t107;
t115 = Icges(5,2) * t109 + t132;
t52 = Icges(5,6) * t102 + t101 * t115;
t53 = Icges(5,6) * t101 - t102 * t115;
t131 = Icges(5,4) * t109;
t116 = Icges(5,1) * t107 + t131;
t54 = Icges(5,5) * t102 + t101 * t116;
t55 = Icges(5,5) * t101 - t102 * t116;
t81 = qJD(4) * t101 + V_base(5);
t82 = qJD(4) * t102 + V_base(4);
t86 = -Icges(5,2) * t107 + t131;
t89 = Icges(5,1) * t109 - t132;
t100 = V_base(6) + qJD(1);
t99 = qJD(2) + t100;
t140 = (t107 * t54 + t109 * t52) * t82 + (t107 * t55 + t109 * t53) * t81 + (t107 * t89 + t109 * t86) * t99;
t138 = -pkin(4) - pkin(5);
t108 = sin(qJ(1));
t137 = pkin(1) * t108;
t110 = cos(qJ(1));
t136 = pkin(1) * t110;
t135 = pkin(6) * t101;
t134 = Icges(2,4) * t108;
t128 = t100 * t136 + V_base(2);
t127 = t137 * V_base(4) + V_base(3);
t126 = V_base(5) * pkin(4) + V_base(1);
t75 = pkin(2) * t102 + qJ(3) * t101;
t123 = -t75 - t136;
t122 = t75 * t99 + t128;
t72 = pkin(2) * t101 - qJ(3) * t102;
t121 = t72 * V_base(4) + t127;
t120 = rSges(5,1) * t107 + rSges(5,2) * t109;
t114 = Icges(5,5) * t107 + Icges(5,6) * t109;
t113 = (Icges(5,3) * t102 + t101 * t114) * t82 + (Icges(5,3) * t101 - t102 * t114) * t81 + (Icges(5,5) * t109 - Icges(5,6) * t107) * t99;
t112 = V_base(5) * pkin(5) - t100 * t137 + t126;
t111 = qJD(3) * t101 + t112;
t103 = Icges(2,4) * t110;
t94 = rSges(2,1) * t110 - rSges(2,2) * t108;
t93 = rSges(5,1) * t109 - rSges(5,2) * t107;
t92 = rSges(2,1) * t108 + rSges(2,2) * t110;
t91 = Icges(2,1) * t110 - t134;
t90 = Icges(2,1) * t108 + t103;
t88 = -Icges(2,2) * t108 + t103;
t87 = Icges(2,2) * t110 + t134;
t80 = -rSges(1,1) * V_base(5) + rSges(1,2) * V_base(4) + V_base(3);
t79 = rSges(1,1) * V_base(6) - rSges(1,3) * V_base(4) + V_base(2);
t78 = -rSges(1,2) * V_base(6) + rSges(1,3) * V_base(5) + V_base(1);
t77 = rSges(3,1) * t102 - rSges(3,2) * t101;
t76 = -rSges(4,2) * t102 + rSges(4,3) * t101;
t74 = rSges(3,1) * t101 + rSges(3,2) * t102;
t73 = -rSges(4,2) * t101 - rSges(4,3) * t102;
t57 = rSges(5,3) * t101 - t102 * t120;
t56 = rSges(5,3) * t102 + t101 * t120;
t49 = V_base(5) * rSges(2,3) - t100 * t92 + t126;
t48 = t100 * t94 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t47 = t92 * V_base(4) - t94 * V_base(5) + V_base(3);
t46 = V_base(5) * rSges(3,3) - t74 * t99 + t112;
t45 = t77 * t99 + (-rSges(3,3) + t138) * V_base(4) + t128;
t44 = V_base(4) * t74 + (-t77 - t136) * V_base(5) + t127;
t43 = V_base(5) * rSges(4,1) + (-t72 - t73) * t99 + t111;
t42 = -qJD(3) * t102 + t76 * t99 + (-rSges(4,1) + t138) * V_base(4) + t122;
t41 = V_base(4) * t73 + (t123 - t76) * V_base(5) + t121;
t40 = V_base(5) * pkin(3) + t81 * t93 + (-t57 - t72 - t135) * t99 + t111;
t39 = t56 * t99 - t82 * t93 + (pkin(6) * t99 - qJD(3)) * t102 + (-pkin(3) + t138) * V_base(4) + t122;
t38 = V_base(4) * t135 - t81 * t56 + t82 * t57 + (-pkin(6) * t102 + t123) * V_base(5) + t121;
t1 = m(1) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + m(2) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + m(3) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + m(4) * (t41 ^ 2 + t42 ^ 2 + t43 ^ 2) / 0.2e1 + m(5) * (t38 ^ 2 + t39 ^ 2 + t40 ^ 2) / 0.2e1 + t82 * (t140 * t101 + t113 * t102) / 0.2e1 + t81 * (t113 * t101 - t140 * t102) / 0.2e1 + ((-t107 * t52 + t109 * t54) * t82 + (-t107 * t53 + t109 * t55) * t81 + (-t107 * t86 + t109 * t89 + Icges(4,1) + Icges(3,3)) * t99) * t99 / 0.2e1 + ((t101 * t146 + t102 * t144 - t108 * t87 + t110 * t90 + Icges(1,4)) * V_base(5) + (t145 * t101 + t143 * t102 - t108 * t88 + t110 * t91 + Icges(1,1)) * V_base(4)) * V_base(4) / 0.2e1 + ((t144 * t101 - t146 * t102 + t108 * t90 + t110 * t87 + Icges(1,2)) * V_base(5) + (t101 * t143 - t102 * t145 + t108 * t91 + t110 * t88 + Icges(1,4)) * V_base(4)) * V_base(5) / 0.2e1 + V_base(5) * t99 * (t101 * t151 - t102 * t150) + V_base(4) * t99 * (t101 * t150 + t102 * t151) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t108 + Icges(2,6) * t110) * V_base(5) + (Icges(2,5) * t110 - Icges(2,6) * t108) * V_base(4) + Icges(2,3) * t100 / 0.2e1) * t100;
T = t1;
