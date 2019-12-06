% Calculate joint inertia matrix for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Mq = S5RPRPR4_inertiaJ_slag_vp1(qJ, ...
  pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_inertiaJ_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_inertiaJ_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR4_inertiaJ_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR4_inertiaJ_slag_vp1: rSges has to be [6x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [6 6]), ...
  'S5RPRPR4_inertiaJ_slag_vp1: Icges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From inertia_joint_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:53:09
% EndTime: 2019-12-05 17:53:13
% DurationCPUTime: 0.91s
% Computational Cost: add. (1929->188), mult. (1340->254), div. (0->0), fcn. (1180->10), ass. (0->96)
t83 = qJ(3) + pkin(9);
t76 = sin(t83);
t78 = cos(t83);
t86 = sin(qJ(3));
t88 = cos(qJ(3));
t141 = Icges(4,5) * t88 + Icges(5,5) * t78 - Icges(4,6) * t86 - Icges(5,6) * t76;
t140 = Icges(4,3) + Icges(5,3);
t84 = qJ(1) + pkin(8);
t77 = sin(t84);
t79 = cos(t84);
t139 = t77 * t79;
t138 = t140 * t79 - t141 * t77;
t137 = t140 * t77 + t141 * t79;
t74 = t77 ^ 2;
t75 = t79 ^ 2;
t136 = t77 / 0.2e1;
t135 = t79 / 0.2e1;
t80 = qJ(5) + t83;
t71 = sin(t80);
t72 = cos(t80);
t90 = Icges(6,5) * t72 - Icges(6,6) * t71;
t22 = Icges(6,3) * t79 - t77 * t90;
t23 = Icges(6,3) * t77 + t79 * t90;
t134 = t77 * (t22 * t139 + t74 * t23) + t79 * (t23 * t139 + t75 * t22);
t133 = pkin(3) * t86;
t87 = sin(qJ(1));
t132 = t87 * pkin(1);
t89 = cos(qJ(1));
t131 = t89 * pkin(1);
t73 = t88 * pkin(3) + pkin(2);
t130 = pkin(2) - t73;
t129 = rSges(4,1) * t88;
t128 = rSges(5,1) * t78;
t127 = rSges(6,1) * t72;
t126 = rSges(5,2) * t76;
t125 = rSges(6,2) * t71;
t124 = t77 * t86;
t69 = t79 * rSges(5,3);
t85 = -qJ(4) - pkin(6);
t123 = t79 * t85;
t122 = t79 * rSges(6,3) + t77 * t125;
t55 = pkin(4) * t78 + t73;
t121 = t55 - t73;
t120 = t77 * t126 + t69;
t119 = rSges(4,2) * t124 + t79 * rSges(4,3);
t118 = t75 + t74;
t117 = Icges(4,4) * t86;
t116 = Icges(4,4) * t88;
t115 = Icges(5,4) * t76;
t114 = Icges(5,4) * t78;
t113 = Icges(6,4) * t71;
t112 = Icges(6,4) * t72;
t46 = Icges(6,5) * t71 + Icges(6,6) * t72;
t93 = -Icges(6,2) * t71 + t112;
t96 = Icges(6,1) * t72 - t113;
t47 = Icges(6,2) * t72 + t113;
t48 = Icges(6,1) * t71 + t112;
t99 = t47 * t71 - t48 * t72;
t111 = (t72 * (Icges(6,6) * t77 + t79 * t93) + t71 * (Icges(6,5) * t77 + t79 * t96) + t77 * t46 - t79 * t99) * t136 + (t72 * (Icges(6,6) * t79 - t77 * t93) + t71 * (Icges(6,5) * t79 - t77 * t96) + t79 * t46 + t77 * t99) * t135;
t49 = t71 * rSges(6,1) + t72 * rSges(6,2);
t110 = pkin(4) * t76 + t49;
t108 = -rSges(4,2) * t86 + t129;
t107 = -t126 + t128;
t106 = -t125 + t127;
t98 = Icges(4,1) * t88 - t117;
t97 = Icges(5,1) * t78 - t115;
t95 = -Icges(4,2) * t86 + t116;
t94 = -Icges(5,2) * t76 + t114;
t82 = -pkin(7) + t85;
t67 = t77 * t85;
t66 = pkin(3) * t124;
t64 = -t89 * rSges(2,1) + t87 * rSges(2,2);
t63 = -t87 * rSges(2,1) - t89 * rSges(2,2);
t62 = t86 * rSges(4,1) + t88 * rSges(4,2);
t53 = t76 * rSges(5,1) + t78 * rSges(5,2);
t44 = -t79 * rSges(3,1) + t77 * rSges(3,2) - t131;
t43 = -t77 * rSges(3,1) - t79 * rSges(3,2) - t132;
t36 = (-t53 - t133) * t79;
t35 = t77 * t53 + t66;
t28 = -t77 * t127 + t122;
t21 = (-pkin(6) - t85) * t79 + t130 * t77;
t20 = t79 * (t77 * rSges(6,3) + t106 * t79);
t19 = t79 * (-t77 * pkin(6) - t130 * t79 - t67);
t18 = -t131 + (-rSges(4,3) - pkin(6)) * t77 + (-pkin(2) - t108) * t79;
t17 = -t132 + t79 * pkin(6) + (-pkin(2) - t129) * t77 + t119;
t16 = (-t110 - t133) * t79;
t15 = t110 * t77 + t66;
t14 = -t77 * rSges(5,3) - t131 + t67 + (-t107 - t73) * t79;
t13 = -t132 - t123 + (-t73 - t128) * t77 + t120;
t12 = -t77 * (-t77 * t129 + t119) + (t77 * rSges(4,3) + t108 * t79) * t79;
t11 = -t131 + (-rSges(6,3) + t82) * t77 + (-t106 - t55) * t79;
t10 = -t132 - t79 * t82 + (-t55 - t127) * t77 + t122;
t9 = -t77 * t28 + t20;
t4 = t19 + t107 * t75 + (t77 * t128 - t120 - t21 + t69) * t77;
t3 = t19 + t20 + (t121 * t79 + t67) * t79 + (t121 * t77 - t123 - t21 - t28) * t77;
t1 = [t72 * t47 + t71 * t48 + t78 * (Icges(5,2) * t78 + t115) + t76 * (Icges(5,1) * t76 + t114) + t88 * (Icges(4,2) * t88 + t117) + t86 * (Icges(4,1) * t86 + t116) + Icges(2,3) + Icges(3,3) + m(2) * (t63 ^ 2 + t64 ^ 2) + m(3) * (t43 ^ 2 + t44 ^ 2) + m(4) * (t17 ^ 2 + t18 ^ 2) + m(5) * (t13 ^ 2 + t14 ^ 2) + m(6) * (t10 ^ 2 + t11 ^ 2); 0; m(3) + m(4) + m(5) + m(6); m(5) * (t36 * t13 + t35 * t14) + m(6) * (t16 * t10 + t15 * t11) + m(4) * (-t17 * t79 + t18 * t77) * t62 + t111 + (t78 * (Icges(5,6) * t77 + t79 * t94) + t76 * (Icges(5,5) * t77 + t79 * t97) + t88 * (Icges(4,6) * t77 + t79 * t95) + t86 * (Icges(4,5) * t77 + t79 * t98)) * t136 + (t78 * (Icges(5,6) * t79 - t77 * t94) + t76 * (Icges(5,5) * t79 - t77 * t97) + t88 * (Icges(4,6) * t79 - t77 * t95) + t86 * (Icges(4,5) * t79 - t77 * t98)) * t135 + (Icges(4,5) * t86 + Icges(5,5) * t76 + Icges(4,6) * t88 + Icges(5,6) * t78) * (t75 / 0.2e1 + t74 / 0.2e1); m(4) * t12 + m(5) * t4 + m(6) * t3; m(6) * (t15 ^ 2 + t16 ^ 2 + t3 ^ 2) + m(5) * (t35 ^ 2 + t36 ^ 2 + t4 ^ 2) + m(4) * (t118 * t62 ^ 2 + t12 ^ 2) + t134 + t138 * t79 * t75 + (t137 * t74 + (t137 * t79 + t138 * t77) * t79) * t77; m(5) * (t77 * t13 + t79 * t14) + m(6) * (t77 * t10 + t79 * t11); 0; m(6) * (t79 * t15 + t77 * t16) + m(5) * (t79 * t35 + t77 * t36); 0.2e1 * (m(5) / 0.2e1 + m(6) / 0.2e1) * t118; m(6) * (-t10 * t79 + t11 * t77) * t49 + t111; m(6) * t9; m(6) * (t9 * t3 + (t15 * t77 - t16 * t79) * t49) + t134; 0; m(6) * (t118 * t49 ^ 2 + t9 ^ 2) + t134;];
%% Postprocessing: Reshape Output
% From vec2symmat_5_matlab.m
res = [t1(1), t1(2), t1(4), t1(7), t1(11); t1(2), t1(3), t1(5), t1(8), t1(12); t1(4), t1(5), t1(6), t1(9), t1(13); t1(7), t1(8), t1(9), t1(10), t1(14); t1(11), t1(12), t1(13), t1(14), t1(15);];
Mq = res;
