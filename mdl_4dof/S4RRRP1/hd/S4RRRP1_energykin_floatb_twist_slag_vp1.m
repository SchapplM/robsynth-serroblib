% Calculate kinetic energy for
% S4RRRP1
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
%   pkin=[a2,a3,a4,d1,d2,d3]';
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

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:55
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4RRRP1_energykin_floatb_twist_slag_vp1(qJ, qJD, V_base, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP1_energykin_floatb_twist_slag_vp1: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP1_energykin_floatb_twist_slag_vp1: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RRRP1_energykin_floatb_twist_slag_vp1: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP1_energykin_floatb_twist_slag_vp1: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRP1_energykin_floatb_twist_slag_vp1: m has to be [5x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [5,3]), ...
  'S4RRRP1_energykin_floatb_twist_slag_vp1: rSges has to be [5x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [5 6]), ...
  'S4RRRP1_energykin_floatb_twist_slag_vp1: Icges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_worldframe_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:54:19
% EndTime: 2018-11-14 13:54:20
% DurationCPUTime: 0.84s
% Computational Cost: add. (504->142), mult. (358->175), div. (0->0), fcn. (184->6), ass. (0->74)
t141 = Icges(4,4) + Icges(5,4);
t140 = Icges(4,1) + Icges(5,1);
t139 = Icges(4,2) + Icges(5,2);
t104 = qJ(1) + qJ(2);
t100 = qJ(3) + t104;
t95 = cos(t100);
t138 = t141 * t95;
t94 = sin(t100);
t137 = t141 * t94;
t136 = rSges(5,1) + pkin(3);
t135 = t139 * t95 + t137;
t134 = -t139 * t94 + t138;
t133 = t140 * t94 + t138;
t132 = t140 * t95 - t137;
t131 = Icges(4,5) + Icges(5,5);
t130 = Icges(4,6) + Icges(5,6);
t129 = -pkin(4) - pkin(5);
t97 = sin(t104);
t128 = pkin(2) * t97;
t98 = cos(t104);
t127 = pkin(2) * t98;
t105 = sin(qJ(1));
t126 = pkin(1) * t105;
t106 = cos(qJ(1));
t125 = pkin(1) * t106;
t124 = Icges(3,4) * t97;
t121 = qJ(4) + rSges(5,3);
t120 = Icges(2,4) * t105;
t119 = -pkin(6) + t129;
t96 = V_base(6) + qJD(1);
t118 = t96 * t125 + V_base(2);
t117 = V_base(4) * t126 + V_base(3);
t116 = V_base(5) * pkin(4) + V_base(1);
t113 = rSges(5,2) * t95 + t136 * t94;
t112 = -rSges(5,2) * t94 + t136 * t95;
t93 = qJD(2) + t96;
t111 = t93 * t127 + t118;
t110 = V_base(4) * t128 + t117;
t109 = -t125 - t127;
t108 = V_base(5) * pkin(5) - t96 * t126 + t116;
t107 = V_base(5) * pkin(6) - t93 * t128 + t108;
t99 = Icges(2,4) * t106;
t92 = Icges(3,4) * t98;
t90 = qJD(3) + t93;
t85 = rSges(2,1) * t106 - t105 * rSges(2,2);
t84 = t105 * rSges(2,1) + rSges(2,2) * t106;
t83 = Icges(2,1) * t106 - t120;
t82 = Icges(2,1) * t105 + t99;
t81 = -Icges(2,2) * t105 + t99;
t80 = Icges(2,2) * t106 + t120;
t77 = -V_base(5) * rSges(1,1) + V_base(4) * rSges(1,2) + V_base(3);
t76 = V_base(6) * rSges(1,1) - V_base(4) * rSges(1,3) + V_base(2);
t75 = -V_base(6) * rSges(1,2) + V_base(5) * rSges(1,3) + V_base(1);
t73 = rSges(3,1) * t98 - rSges(3,2) * t97;
t72 = rSges(3,1) * t97 + rSges(3,2) * t98;
t71 = Icges(3,1) * t98 - t124;
t70 = Icges(3,1) * t97 + t92;
t69 = -Icges(3,2) * t97 + t92;
t68 = Icges(3,2) * t98 + t124;
t65 = rSges(4,1) * t95 - rSges(4,2) * t94;
t63 = rSges(4,1) * t94 + rSges(4,2) * t95;
t49 = V_base(5) * rSges(2,3) - t84 * t96 + t116;
t48 = t85 * t96 + V_base(2) + (-rSges(2,3) - pkin(4)) * V_base(4);
t47 = t84 * V_base(4) - t85 * V_base(5) + V_base(3);
t46 = V_base(5) * rSges(3,3) - t72 * t93 + t108;
t45 = t73 * t93 + (-rSges(3,3) + t129) * V_base(4) + t118;
t44 = V_base(4) * t72 + (-t73 - t125) * V_base(5) + t117;
t43 = V_base(5) * rSges(4,3) - t63 * t90 + t107;
t42 = t65 * t90 + (-rSges(4,3) + t119) * V_base(4) + t111;
t41 = V_base(4) * t63 + (t109 - t65) * V_base(5) + t110;
t40 = -t113 * t90 + t121 * V_base(5) + t107;
t39 = t112 * t90 + (t119 - t121) * V_base(4) + t111;
t38 = qJD(4) + t113 * V_base(4) + (t109 - t112) * V_base(5) + t110;
t1 = m(1) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + Icges(1,1) * V_base(4) ^ 2 / 0.2e1 + Icges(1,2) * V_base(5) ^ 2 / 0.2e1 + m(2) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + m(3) * (t44 ^ 2 + t45 ^ 2 + t46 ^ 2) / 0.2e1 + m(4) * (t41 ^ 2 + t42 ^ 2 + t43 ^ 2) / 0.2e1 + m(5) * (t38 ^ 2 + t39 ^ 2 + t40 ^ 2) / 0.2e1 + V_base(4) * V_base(5) * Icges(1,4) + (Icges(1,5) * V_base(4) + Icges(1,6) * V_base(5) + Icges(1,3) * V_base(6) / 0.2e1) * V_base(6) + ((Icges(2,5) * t105 + Icges(2,6) * t106) * V_base(5) + (Icges(2,5) * t106 - Icges(2,6) * t105) * V_base(4) + Icges(2,3) * t96 / 0.2e1) * t96 + ((Icges(3,5) * t97 + Icges(3,6) * t98) * V_base(5) + (Icges(3,5) * t98 - Icges(3,6) * t97) * V_base(4) + Icges(3,3) * t93 / 0.2e1) * t93 + ((Icges(4,3) / 0.2e1 + Icges(5,3) / 0.2e1) * t90 + (t130 * V_base(5) + t131 * V_base(4)) * t95 + (-t130 * V_base(4) + t131 * V_base(5)) * t94) * t90 + ((-t105 * t80 + t106 * t82 + t133 * t95 - t135 * t94 - t68 * t97 + t70 * t98) * V_base(5) + (-t105 * t81 + t106 * t83 + t132 * t95 - t134 * t94 - t69 * t97 + t71 * t98) * V_base(4)) * V_base(4) / 0.2e1 + ((t105 * t82 + t106 * t80 + t133 * t94 + t135 * t95 + t68 * t98 + t70 * t97) * V_base(5) + (t105 * t83 + t106 * t81 + t132 * t94 + t134 * t95 + t69 * t98 + t71 * t97) * V_base(4)) * V_base(5) / 0.2e1;
T  = t1;
