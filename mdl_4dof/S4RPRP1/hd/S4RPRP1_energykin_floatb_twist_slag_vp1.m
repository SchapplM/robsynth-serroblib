% Calculate kinetic energy for
% S4RPRP1
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
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP1_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RPRP1_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RPRP1_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:29:26
% EndTime: 2019-03-08 18:29:27
% DurationCPUTime: 0.78s
% Computational Cost: add. (498->144), mult. (362->172), div. (0->0), fcn. (190->6), ass. (0->73)
t144 = Icges(4,4) - Icges(5,5);
t143 = Icges(4,1) + Icges(5,1);
t142 = Icges(4,2) + Icges(5,3);
t107 = qJ(1) + pkin(6);
t101 = qJ(3) + t107;
t96 = sin(t101);
t141 = t144 * t96;
t97 = cos(t101);
t140 = t144 * t97;
t139 = rSges(5,1) + pkin(3);
t138 = -t142 * t97 - t141;
t137 = t142 * t96 - t140;
t136 = t143 * t96 + t140;
t135 = t143 * t97 - t141;
t134 = rSges(5,3) + qJ(4);
t133 = Icges(5,4) + Icges(4,5);
t132 = Icges(4,6) - Icges(5,6);
t99 = sin(t107);
t131 = pkin(2) * t99;
t108 = sin(qJ(1));
t130 = pkin(1) * t108;
t109 = cos(qJ(1));
t129 = pkin(1) * t109;
t100 = cos(t107);
t128 = pkin(2) * t100;
t127 = -pkin(4) - qJ(2);
t126 = -t134 * t97 + t139 * t96;
t125 = t134 * t96 + t139 * t97;
t124 = Icges(3,4) * t99;
t121 = Icges(2,4) * t108;
t102 = V_base(6) + qJD(1);
t120 = t102 * t129 + V_base(2);
t119 = -pkin(5) + t127;
t118 = V_base(5) * pkin(4) + V_base(1);
t115 = t102 * t128 + t120;
t114 = V_base(4) * t130 + qJD(2) + V_base(3);
t113 = V_base(5) * qJ(2) + t118;
t112 = V_base(4) * t131 + t114;
t111 = -t128 - t129;
t110 = V_base(5) * pkin(5) + (-t130 - t131) * t102 + t113;
t104 = Icges(2,4) * t109;
t98 = qJD(3) + t102;
t95 = Icges(3,4) * t100;
t89 = t109 * rSges(2,1) - t108 * rSges(2,2);
t88 = t108 * rSges(2,1) + t109 * rSges(2,2);
t87 = Icges(2,1) * t109 - t121;
t86 = Icges(2,1) * t108 + t104;
t85 = -Icges(2,2) * t108 + t104;
t84 = Icges(2,2) * t109 + t121;
t80 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t79 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t78 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t77 = t100 * rSges(3,1) - t99 * rSges(3,2);
t76 = t99 * rSges(3,1) + t100 * rSges(3,2);
t75 = Icges(3,1) * t100 - t124;
t74 = Icges(3,1) * t99 + t95;
t73 = -Icges(3,2) * t99 + t95;
t72 = Icges(3,2) * t100 + t124;
t69 = t97 * rSges(4,1) - t96 * rSges(4,2);
t66 = t96 * rSges(4,1) + t97 * rSges(4,2);
t51 = V_base(5) * rSges(2,3) - t102 * t88 + t118;
t50 = t102 * t89 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t49 = V_base(4) * t88 - V_base(5) * t89 + V_base(3);
t48 = V_base(5) * rSges(3,3) + (-t76 - t130) * t102 + t113;
t47 = t102 * t77 + (-rSges(3,3) + t127) * V_base(4) + t120;
t46 = V_base(4) * t76 + (-t77 - t129) * V_base(5) + t114;
t45 = V_base(5) * rSges(4,3) - t98 * t66 + t110;
t44 = t98 * t69 + (-rSges(4,3) + t119) * V_base(4) + t115;
t43 = V_base(4) * t66 + (t111 - t69) * V_base(5) + t112;
t42 = V_base(5) * rSges(5,2) + qJD(4) * t96 - t126 * t98 + t110;
t41 = -qJD(4) * t97 + t125 * t98 + (-rSges(5,2) + t119) * V_base(4) + t115;
t40 = t126 * V_base(4) + (t111 - t125) * V_base(5) + t112;
t1 = m(1) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + Icges(1,1) * V_base(4) ^ 2 / 0.2e1 + Icges(1,2) * V_base(5) ^ 2 / 0.2e1 + m(2) * (t49 ^ 2 + t50 ^ 2 + t51 ^ 2) / 0.2e1 + m(3) * (t46 ^ 2 + t47 ^ 2 + t48 ^ 2) / 0.2e1 + m(4) * (t43 ^ 2 + t44 ^ 2 + t45 ^ 2) / 0.2e1 + m(5) * (t40 ^ 2 + t41 ^ 2 + t42 ^ 2) / 0.2e1 + V_base(4) * V_base(5) * Icges(1,4) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(4,3) / 0.2e1 + Icges(5,2) / 0.2e1) * t98 + (t132 * V_base(5) + t133 * V_base(4)) * t97 + (-t132 * V_base(4) + t133 * V_base(5)) * t96) * t98 + ((Icges(2,3) / 0.2e1 + Icges(3,3) / 0.2e1) * t102 + (Icges(2,5) * t108 + Icges(3,5) * t99 + Icges(2,6) * t109 + Icges(3,6) * t100) * V_base(5) + (Icges(2,5) * t109 + Icges(3,5) * t100 - Icges(2,6) * t108 - Icges(3,6) * t99) * V_base(4)) * t102 + ((t100 * t74 - t108 * t84 + t109 * t86 + t136 * t97 + t138 * t96 - t99 * t72) * V_base(5) + (t100 * t75 - t108 * t85 + t109 * t87 + t135 * t97 + t137 * t96 - t99 * t73) * V_base(4)) * V_base(4) / 0.2e1 + ((t100 * t72 + t108 * t86 + t109 * t84 + t136 * t96 - t138 * t97 + t99 * t74) * V_base(5) + (t100 * t73 + t108 * t87 + t109 * t85 + t135 * t96 - t137 * t97 + t99 * t75) * V_base(4)) * V_base(5) / 0.2e1;
T  = t1;
