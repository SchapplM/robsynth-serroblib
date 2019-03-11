% Calculate kinetic energy for
% S4PRRP1
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
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2019-03-08 18:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PRRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PRRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP1_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4PRRP1_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4PRRP1_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:22:47
% EndTime: 2019-03-08 18:22:48
% DurationCPUTime: 0.75s
% Computational Cost: add. (487->143), mult. (362->175), div. (0->0), fcn. (190->6), ass. (0->73)
t144 = Icges(4,4) - Icges(5,5);
t143 = Icges(4,1) + Icges(5,1);
t142 = Icges(4,2) + Icges(5,3);
t107 = pkin(6) + qJ(2);
t102 = qJ(3) + t107;
t96 = sin(t102);
t141 = t144 * t96;
t97 = cos(t102);
t140 = t144 * t97;
t139 = rSges(5,1) + pkin(3);
t138 = -t142 * t97 - t141;
t137 = t142 * t96 - t140;
t136 = t143 * t96 + t140;
t135 = t143 * t97 - t141;
t134 = rSges(5,3) + qJ(4);
t133 = Icges(5,4) + Icges(4,5);
t132 = Icges(4,6) - Icges(5,6);
t109 = cos(pkin(6));
t131 = pkin(1) * t109;
t103 = V_base(6) + qJD(2);
t130 = pkin(2) * t103;
t129 = -pkin(4) - qJ(1);
t128 = -t134 * t97 + t139 * t96;
t127 = t134 * t96 + t139 * t97;
t99 = sin(t107);
t126 = Icges(3,4) * t99;
t108 = sin(pkin(6));
t123 = Icges(2,4) * t108;
t116 = pkin(1) * V_base(6);
t122 = t109 * t116 + V_base(2);
t121 = -pkin(5) + t129;
t120 = V_base(5) * qJ(1) + V_base(1);
t117 = qJD(1) + V_base(3);
t100 = cos(t107);
t115 = t100 * t130 + t122;
t114 = V_base(4) * t108 * pkin(1) + t117;
t113 = V_base(4) * pkin(2) * t99 + t114;
t112 = -pkin(2) * t100 - t131;
t111 = V_base(5) * pkin(4) - t108 * t116 + t120;
t110 = V_base(5) * pkin(5) - t99 * t130 + t111;
t101 = Icges(2,4) * t109;
t98 = qJD(3) + t103;
t95 = Icges(3,4) * t100;
t89 = rSges(2,1) * t109 - t108 * rSges(2,2);
t88 = t108 * rSges(2,1) + rSges(2,2) * t109;
t87 = Icges(2,1) * t109 - t123;
t86 = Icges(2,1) * t108 + t101;
t85 = -Icges(2,2) * t108 + t101;
t84 = Icges(2,2) * t109 + t123;
t80 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t79 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t78 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t77 = rSges(3,1) * t100 - rSges(3,2) * t99;
t76 = rSges(3,1) * t99 + rSges(3,2) * t100;
t75 = Icges(3,1) * t100 - t126;
t74 = Icges(3,1) * t99 + t95;
t73 = -Icges(3,2) * t99 + t95;
t72 = Icges(3,2) * t100 + t126;
t69 = rSges(4,1) * t97 - rSges(4,2) * t96;
t66 = rSges(4,1) * t96 + rSges(4,2) * t97;
t51 = V_base(5) * rSges(2,3) - t88 * V_base(6) + t120;
t50 = t89 * V_base(6) + V_base(2) + (-rSges(2,3) - qJ(1)) * V_base(4);
t49 = t88 * V_base(4) - t89 * V_base(5) + t117;
t48 = V_base(5) * rSges(3,3) - t103 * t76 + t111;
t47 = t103 * t77 + (-rSges(3,3) + t129) * V_base(4) + t122;
t46 = V_base(4) * t76 + (-t77 - t131) * V_base(5) + t114;
t45 = V_base(5) * rSges(4,3) - t66 * t98 + t110;
t44 = t69 * t98 + (-rSges(4,3) + t121) * V_base(4) + t115;
t43 = V_base(4) * t66 + (t112 - t69) * V_base(5) + t113;
t42 = V_base(5) * rSges(5,2) + qJD(4) * t96 - t128 * t98 + t110;
t41 = -qJD(4) * t97 + t127 * t98 + (-rSges(5,2) + t121) * V_base(4) + t115;
t40 = t128 * V_base(4) + (t112 - t127) * V_base(5) + t113;
t1 = m(1) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + Icges(1,1) * V_base(4) ^ 2 / 0.2e1 + Icges(1,2) * V_base(5) ^ 2 / 0.2e1 + m(2) * (t49 ^ 2 + t50 ^ 2 + t51 ^ 2) / 0.2e1 + m(3) * (t46 ^ 2 + t47 ^ 2 + t48 ^ 2) / 0.2e1 + m(4) * (t43 ^ 2 + t44 ^ 2 + t45 ^ 2) / 0.2e1 + m(5) * (t40 ^ 2 + t41 ^ 2 + t42 ^ 2) / 0.2e1 + V_base(4) * V_base(5) * Icges(1,4) + ((Icges(3,5) * t99 + Icges(3,6) * t100) * V_base(5) + (Icges(3,5) * t100 - Icges(3,6) * t99) * V_base(4) + Icges(3,3) * t103 / 0.2e1) * t103 + ((Icges(1,3) / 0.2e1 + Icges(2,3) / 0.2e1) * V_base(6) + (Icges(2,5) * t108 + Icges(2,6) * t109 + Icges(1,6)) * V_base(5) + (Icges(2,5) * t109 - Icges(2,6) * t108 + Icges(1,5)) * V_base(4)) * V_base(6) + ((Icges(4,3) / 0.2e1 + Icges(5,2) / 0.2e1) * t98 + (t132 * V_base(5) + t133 * V_base(4)) * t97 + (-t132 * V_base(4) + t133 * V_base(5)) * t96) * t98 + ((t100 * t74 - t108 * t84 + t109 * t86 + t136 * t97 + t138 * t96 - t72 * t99) * V_base(5) + (t100 * t75 - t108 * t85 + t109 * t87 + t135 * t97 + t137 * t96 - t73 * t99) * V_base(4)) * V_base(4) / 0.2e1 + ((t100 * t72 + t108 * t86 + t109 * t84 + t136 * t96 - t138 * t97 + t74 * t99) * V_base(5) + (t100 * t73 + t108 * t87 + t109 * t85 + t135 * t96 - t137 * t97 + t75 * t99) * V_base(4)) * V_base(5) / 0.2e1;
T  = t1;
