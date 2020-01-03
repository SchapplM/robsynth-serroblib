% Calculate joint inertia matrix for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Mq [5x5]
%   inertia matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RRRRR5_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRRR5_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RRRRR5_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:03
% EndTime: 2020-01-03 12:13:06
% DurationCPUTime: 0.73s
% Computational Cost: add. (3060->178), mult. (1852->252), div. (0->0), fcn. (1594->10), ass. (0->98)
t101 = qJ(1) + qJ(2);
t97 = qJ(3) + t101;
t90 = sin(t97);
t91 = cos(t97);
t144 = t90 * t91;
t102 = sin(qJ(4));
t104 = cos(qJ(4));
t146 = -rSges(5,1) * t104 + rSges(5,2) * t102;
t100 = qJ(4) + qJ(5);
t93 = sin(t100);
t95 = cos(t100);
t145 = -rSges(6,1) * t95 + rSges(6,2) * t93;
t33 = -t90 * rSges(6,3) + t145 * t91;
t92 = t104 * pkin(4) + pkin(3);
t143 = t91 * t92 - t33;
t83 = t90 ^ 2;
t84 = t91 ^ 2;
t142 = -t90 / 0.2e1;
t141 = -t91 / 0.2e1;
t71 = t102 * rSges(5,1) + t104 * rSges(5,2);
t140 = m(5) * t71;
t56 = t93 * rSges(6,1) + t95 * rSges(6,2);
t139 = m(6) * t56;
t106 = -pkin(9) - pkin(8);
t136 = t91 * t106 + t90 * t92;
t50 = t90 * rSges(4,1) + t91 * rSges(4,2);
t135 = -t91 * pkin(3) - t90 * pkin(8);
t134 = t83 + t84;
t94 = sin(t101);
t96 = cos(t101);
t57 = t94 * rSges(3,1) + t96 * rSges(3,2);
t131 = Icges(6,4) * t93;
t130 = Icges(6,4) * t95;
t129 = Icges(5,4) * t102;
t128 = Icges(5,4) * t104;
t88 = pkin(2) * t94;
t44 = t88 + t50;
t118 = -Icges(6,2) * t93 + t130;
t119 = Icges(6,1) * t95 - t131;
t54 = Icges(6,2) * t95 + t131;
t55 = Icges(6,1) * t93 + t130;
t120 = t54 * t93 - t55 * t95;
t53 = Icges(6,5) * t93 + Icges(6,6) * t95;
t127 = (t120 * t91 + t95 * (-Icges(6,6) * t90 - t118 * t91) + t93 * (-Icges(6,5) * t90 - t119 * t91) - t90 * t53) * t142 + (-t120 * t90 + t95 * (-Icges(6,6) * t91 + t118 * t90) + t93 * (-Icges(6,5) * t91 + t119 * t90) - t91 * t53) * t141;
t58 = t96 * rSges(3,1) - t94 * rSges(3,2);
t51 = t91 * rSges(4,1) - t90 * rSges(4,2);
t126 = pkin(4) * t102 + t56;
t89 = pkin(2) * t96;
t45 = t51 + t89;
t125 = t146 * t90;
t117 = Icges(6,5) * t95 - Icges(6,6) * t93;
t27 = -Icges(6,3) * t91 + t117 * t90;
t28 = -Icges(6,3) * t90 - t117 * t91;
t124 = -t91 * (t28 * t144 + t84 * t27) - t90 * (t27 * t144 + t83 * t28);
t69 = Icges(5,2) * t104 + t129;
t70 = Icges(5,1) * t102 + t128;
t123 = t102 * t70 + t104 * t69 + t95 * t54 + t93 * t55 + Icges(4,3);
t114 = t102 * t69 - t104 * t70;
t113 = -t90 * rSges(5,3) + t146 * t91;
t112 = Icges(5,1) * t104 - t129;
t111 = -Icges(5,2) * t102 + t128;
t110 = Icges(5,5) * t104 - Icges(5,6) * t102;
t109 = Icges(3,3) + t123;
t68 = Icges(5,5) * t102 + Icges(5,6) * t104;
t108 = t127 + (t102 * (-Icges(5,5) * t90 - t112 * t91) + t104 * (-Icges(5,6) * t90 - t111 * t91) + t114 * t91 - t90 * t68) * t142 + (t102 * (-Icges(5,5) * t91 + t112 * t90) + t104 * (-Icges(5,6) * t91 + t111 * t90) - t114 * t90 - t91 * t68) * t141;
t107 = -t91 * rSges(6,3) - t145 * t90;
t25 = -t113 - t135;
t23 = t25 + t89;
t20 = t107 + t136;
t21 = -t90 * t106 + t143;
t81 = t90 * pkin(3);
t24 = t81 + (-rSges(5,3) - pkin(8)) * t91 - t125;
t16 = t20 + t88;
t17 = t21 + t89;
t22 = t88 + t24;
t105 = cos(qJ(1));
t103 = sin(qJ(1));
t99 = t105 * pkin(1);
t98 = t103 * pkin(1);
t73 = t105 * rSges(2,1) - t103 * rSges(2,2);
t72 = t103 * rSges(2,1) + t105 * rSges(2,2);
t49 = t58 + t99;
t48 = t98 + t57;
t43 = t45 + t99;
t42 = t98 + t44;
t37 = -Icges(5,3) * t90 - t110 * t91;
t36 = -Icges(5,3) * t91 + t110 * t90;
t35 = t126 * t91;
t34 = t126 * t90;
t26 = t90 * t107;
t19 = t23 + t99;
t18 = t98 + t22;
t13 = t17 + t99;
t12 = t16 + t98;
t11 = -t91 * t113 + t90 * (-t91 * rSges(5,3) - t125);
t6 = -t91 * t33 + t26;
t3 = t90 * (-t81 + t136) + t26 + ((pkin(8) - t106) * t90 + t135 + t143) * t91;
t1 = [Icges(2,3) + m(2) * (t72 ^ 2 + t73 ^ 2) + m(3) * (t48 ^ 2 + t49 ^ 2) + m(4) * (t42 ^ 2 + t43 ^ 2) + m(5) * (t18 ^ 2 + t19 ^ 2) + m(6) * (t12 ^ 2 + t13 ^ 2) + t109; m(3) * (t57 * t48 + t58 * t49) + m(4) * (t44 * t42 + t45 * t43) + m(5) * (t22 * t18 + t23 * t19) + m(6) * (t16 * t12 + t17 * t13) + t109; m(6) * (t16 ^ 2 + t17 ^ 2) + m(5) * (t22 ^ 2 + t23 ^ 2) + m(4) * (t44 ^ 2 + t45 ^ 2) + m(3) * (t57 ^ 2 + t58 ^ 2) + t109; m(4) * (t50 * t42 + t51 * t43) + m(5) * (t24 * t18 + t25 * t19) + m(6) * (t20 * t12 + t21 * t13) + t123; m(6) * (t20 * t16 + t21 * t17) + m(5) * (t24 * t22 + t25 * t23) + m(4) * (t50 * t44 + t51 * t45) + t123; m(6) * (t20 ^ 2 + t21 ^ 2) + m(5) * (t24 ^ 2 + t25 ^ 2) + m(4) * (t50 ^ 2 + t51 ^ 2) + t123; m(6) * (t35 * t12 - t34 * t13) + (t18 * t91 - t19 * t90) * t140 + t108; m(6) * (t35 * t16 - t34 * t17) + (t22 * t91 - t23 * t90) * t140 + t108; m(6) * (t35 * t20 - t34 * t21) + (t24 * t91 - t25 * t90) * t140 + t108; m(5) * (t134 * t71 ^ 2 + t11 ^ 2) - t91 * (t37 * t144 + t84 * t36) - t90 * (t36 * t144 + t83 * t37) + m(6) * (t3 ^ 2 + t34 ^ 2 + t35 ^ 2) + t124; (t12 * t91 - t13 * t90) * t139 + t127; (t16 * t91 - t17 * t90) * t139 + t127; (t20 * t91 - t21 * t90) * t139 + t127; m(6) * (t6 * t3 + (t34 * t90 + t35 * t91) * t56) + t124; m(6) * (t134 * t56 ^ 2 + t6 ^ 2) + t124;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
