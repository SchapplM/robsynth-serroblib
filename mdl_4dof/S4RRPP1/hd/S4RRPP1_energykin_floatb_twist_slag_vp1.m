% Calculate kinetic energy for
% S4RRPP1
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
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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
% Datum: 2019-03-08 18:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRPP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP1_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRPP1_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRPP1_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:32:57
% EndTime: 2019-03-08 18:32:58
% DurationCPUTime: 0.76s
% Computational Cost: add. (507->144), mult. (362->169), div. (0->0), fcn. (190->6), ass. (0->75)
t146 = Icges(4,4) - Icges(5,5);
t145 = Icges(4,1) + Icges(5,1);
t144 = Icges(4,2) + Icges(5,3);
t107 = qJ(1) + qJ(2);
t99 = pkin(6) + t107;
t96 = sin(t99);
t143 = t146 * t96;
t97 = cos(t99);
t142 = t146 * t97;
t141 = rSges(5,1) + pkin(3);
t140 = -t144 * t97 - t143;
t139 = t144 * t96 - t142;
t138 = t145 * t96 + t142;
t137 = t145 * t97 - t143;
t136 = rSges(5,3) + qJ(4);
t135 = V_base(4) / 0.2e1;
t134 = V_base(5) / 0.2e1;
t133 = Icges(5,4) + Icges(4,5);
t132 = Icges(4,6) - Icges(5,6);
t131 = -pkin(4) - pkin(5);
t108 = sin(qJ(1));
t130 = pkin(1) * t108;
t109 = cos(qJ(1));
t129 = pkin(1) * t109;
t101 = sin(t107);
t128 = pkin(2) * t101;
t102 = cos(t107);
t127 = pkin(2) * t102;
t126 = -t136 * t97 + t141 * t96;
t125 = t136 * t96 + t141 * t97;
t122 = Icges(2,4) * t108;
t121 = Icges(3,4) * t101;
t100 = V_base(6) + qJD(1);
t120 = t100 * t129 + V_base(2);
t119 = V_base(4) * t130 + V_base(3);
t118 = -qJ(3) + t131;
t117 = V_base(5) * pkin(4) + V_base(1);
t98 = qJD(2) + t100;
t114 = t98 * t127 + t120;
t113 = V_base(4) * t128 + qJD(3) + t119;
t112 = -t127 - t129;
t111 = V_base(5) * pkin(5) - t100 * t130 + t117;
t110 = V_base(5) * qJ(3) + t111;
t104 = Icges(2,4) * t109;
t95 = Icges(3,4) * t102;
t89 = rSges(2,1) * t109 - t108 * rSges(2,2);
t88 = t108 * rSges(2,1) + rSges(2,2) * t109;
t87 = Icges(2,1) * t109 - t122;
t86 = Icges(2,1) * t108 + t104;
t85 = -Icges(2,2) * t108 + t104;
t84 = Icges(2,2) * t109 + t122;
t81 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t80 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t79 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t77 = rSges(3,1) * t102 - rSges(3,2) * t101;
t76 = rSges(3,1) * t101 + rSges(3,2) * t102;
t75 = Icges(3,1) * t102 - t121;
t74 = Icges(3,1) * t101 + t95;
t73 = -Icges(3,2) * t101 + t95;
t72 = Icges(3,2) * t102 + t121;
t69 = rSges(4,1) * t97 - rSges(4,2) * t96;
t66 = rSges(4,1) * t96 + rSges(4,2) * t97;
t51 = V_base(5) * rSges(2,3) - t100 * t88 + t117;
t50 = t100 * t89 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t49 = t88 * V_base(4) - t89 * V_base(5) + V_base(3);
t48 = V_base(5) * rSges(3,3) - t76 * t98 + t111;
t47 = t77 * t98 + (-rSges(3,3) + t131) * V_base(4) + t120;
t46 = V_base(4) * t76 + (-t77 - t129) * V_base(5) + t119;
t45 = V_base(5) * rSges(4,3) + (-t66 - t128) * t98 + t110;
t44 = t69 * t98 + (-rSges(4,3) + t118) * V_base(4) + t114;
t43 = V_base(4) * t66 + (t112 - t69) * V_base(5) + t113;
t42 = V_base(5) * rSges(5,2) + qJD(4) * t96 + (-t126 - t128) * t98 + t110;
t41 = -qJD(4) * t97 + t125 * t98 + (-rSges(5,2) + t118) * V_base(4) + t114;
t40 = t126 * V_base(4) + (t112 - t125) * V_base(5) + t113;
t1 = m(1) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + m(2) * (t49 ^ 2 + t50 ^ 2 + t51 ^ 2) / 0.2e1 + m(3) * (t46 ^ 2 + t47 ^ 2 + t48 ^ 2) / 0.2e1 + m(4) * (t43 ^ 2 + t44 ^ 2 + t45 ^ 2) / 0.2e1 + m(5) * (t40 ^ 2 + t41 ^ 2 + t42 ^ 2) / 0.2e1 + Icges(1,3) * V_base(6) ^ 2 / 0.2e1 + Icges(2,3) * t100 ^ 2 / 0.2e1 + (Icges(3,3) / 0.2e1 + Icges(4,3) / 0.2e1 + Icges(5,2) / 0.2e1) * t98 ^ 2 + (Icges(1,6) * V_base(6) + (Icges(2,5) * t108 + Icges(2,6) * t109) * t100 + (Icges(3,5) * t101 + Icges(3,6) * t102 + t132 * t97 + t133 * t96) * t98 + Icges(1,2) * t134) * V_base(5) + (Icges(1,4) * V_base(5) + Icges(1,5) * V_base(6) + (Icges(2,5) * t109 - Icges(2,6) * t108) * t100 + (Icges(3,5) * t102 - Icges(3,6) * t101 - t132 * t96 + t133 * t97) * t98 + Icges(1,1) * t135) * V_base(4) + ((-t101 * t72 + t102 * t74 - t108 * t84 + t109 * t86 + t138 * t97 + t140 * t96) * V_base(5) + (-t101 * t73 + t102 * t75 - t108 * t85 + t109 * t87 + t137 * t97 + t139 * t96) * V_base(4)) * t135 + ((t101 * t74 + t102 * t72 + t108 * t86 + t109 * t84 + t138 * t96 - t140 * t97) * V_base(5) + (t101 * t75 + t102 * t73 + t108 * t87 + t109 * t85 + t137 * t96 - t139 * t97) * V_base(4)) * t134;
T  = t1;
