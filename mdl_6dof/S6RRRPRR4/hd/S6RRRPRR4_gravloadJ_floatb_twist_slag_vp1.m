% Calculate Gravitation load on the joints for
% S6RRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% rSges [7x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [6x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:15:25
% EndTime: 2019-03-09 18:15:28
% DurationCPUTime: 1.07s
% Computational Cost: add. (598->156), mult. (562->207), div. (0->0), fcn. (517->12), ass. (0->79)
t56 = sin(pkin(11));
t104 = rSges(5,2) * t56;
t55 = qJ(2) + qJ(3);
t49 = sin(t55);
t50 = cos(t55);
t57 = cos(pkin(11));
t85 = -rSges(5,1) * t57 - pkin(3);
t70 = (rSges(5,3) + qJ(4)) * t49 + (-t85 - t104) * t50;
t126 = qJ(4) * t50 + t49 * t104;
t43 = t57 * pkin(4) + pkin(3);
t54 = pkin(11) + qJ(5);
t47 = cos(t54);
t19 = pkin(5) * t47 + t43;
t124 = t49 * rSges(7,3) + t50 * t19;
t123 = t49 * rSges(6,3) + t50 * t43;
t121 = t50 * rSges(4,1) - rSges(4,2) * t49;
t63 = -pkin(8) - pkin(7);
t84 = pkin(4) * t56 - t63;
t60 = sin(qJ(1));
t62 = cos(qJ(1));
t120 = g(1) * t62 + g(2) * t60;
t119 = t120 * t49;
t48 = qJ(6) + t54;
t41 = sin(t48);
t42 = cos(t48);
t99 = t50 * t60;
t6 = t41 * t99 + t42 * t62;
t7 = t41 * t62 - t42 * t99;
t118 = -t6 * rSges(7,1) + t7 * rSges(7,2);
t98 = t50 * t62;
t8 = -t41 * t98 + t42 * t60;
t9 = t41 * t60 + t42 * t98;
t117 = t8 * rSges(7,1) - t9 * rSges(7,2);
t59 = sin(qJ(2));
t116 = pkin(2) * t59;
t46 = sin(t54);
t114 = pkin(5) * t46;
t61 = cos(qJ(2));
t52 = t61 * pkin(2);
t45 = t52 + pkin(1);
t29 = t62 * t45;
t112 = g(2) * t29;
t110 = g(3) * t49;
t109 = rSges(3,3) + pkin(7);
t107 = rSges(6,1) * t47;
t106 = rSges(7,1) * t42;
t103 = rSges(6,2) * t46;
t102 = rSges(7,2) * t41;
t58 = -pkin(9) - qJ(4);
t53 = -pkin(10) + t58;
t101 = t49 * t53;
t100 = t49 * t58;
t97 = rSges(4,3) - t63;
t88 = t49 * t102;
t96 = rSges(7,3) * t99 + t60 * t88;
t95 = rSges(7,3) * t98 + t62 * t88;
t89 = t49 * t103;
t94 = rSges(6,3) * t99 + t60 * t89;
t93 = rSges(6,3) * t98 + t62 * t89;
t92 = t114 + t84;
t87 = rSges(5,3) * t99 + t126 * t60;
t86 = rSges(5,3) * t98 + t126 * t62;
t83 = -t43 - t107;
t82 = -t19 - t106;
t80 = rSges(3,1) * t61 - rSges(3,2) * t59;
t77 = -rSges(4,1) * t49 - rSges(4,2) * t50;
t76 = -rSges(7,1) * t41 - rSges(7,2) * t42;
t75 = pkin(1) + t80;
t12 = -t46 * t98 + t47 * t60;
t10 = t46 * t99 + t47 * t62;
t73 = rSges(5,1) * t56 + rSges(5,2) * t57 - t63;
t72 = t123 + (-t103 + t107) * t50;
t71 = t124 + (-t102 + t106) * t50;
t68 = -t100 + t123;
t67 = -t101 + t124;
t64 = t85 * t119;
t13 = t46 * t60 + t47 * t98;
t11 = t46 * t62 - t47 * t99;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t60 - rSges(2,2) * t62) + g(2) * (rSges(2,1) * t62 - rSges(2,2) * t60)) - m(3) * ((g(1) * t109 + g(2) * t75) * t62 + (-g(1) * t75 + g(2) * t109) * t60) - m(4) * (t112 + (g(1) * t97 + g(2) * t121) * t62 + (g(1) * (-t45 - t121) + g(2) * t97) * t60) - m(5) * (t112 + (g(1) * t73 + t70 * g(2)) * t62 + (g(2) * t73 + (-t45 - t70) * g(1)) * t60) - m(6) * (g(1) * (t11 * rSges(6,1) + t10 * rSges(6,2)) + g(2) * (t13 * rSges(6,1) + t12 * rSges(6,2) + t29) + (g(1) * t84 + g(2) * t68) * t62 + (g(1) * (-t45 - t68) + g(2) * t84) * t60) - m(7) * (g(1) * (t7 * rSges(7,1) + t6 * rSges(7,2)) + g(2) * (t9 * rSges(7,1) + t8 * rSges(7,2) + t29) + (g(1) * t92 + g(2) * t67) * t62 + (g(1) * (-t45 - t67) + g(2) * t92) * t60) -m(3) * (g(3) * t80 + t120 * (-rSges(3,1) * t59 - rSges(3,2) * t61)) - m(4) * (g(3) * (t52 + t121) + t120 * (t77 - t116)) - m(5) * (g(1) * (-t116 * t62 + t86) + g(2) * (-t116 * t60 + t87) + g(3) * (t52 + t70) + t64) - m(6) * (g(1) * t93 + g(2) * t94 + g(3) * (t52 + t72 - t100) + t120 * (t49 * t83 - t50 * t58 - t116)) - m(7) * (g(1) * t95 + g(2) * t96 + g(3) * (t52 + t71 - t101) + t120 * (t49 * t82 - t50 * t53 - t116)) -m(4) * (g(3) * t121 + t120 * t77) - m(5) * (g(1) * t86 + g(2) * t87 + g(3) * t70 + t64) - m(6) * (g(1) * (-t58 * t98 + t93) + g(2) * (-t58 * t99 + t94) + g(3) * t72 + (-g(3) * t58 + t120 * t83) * t49) - m(7) * (g(1) * (-t53 * t98 + t95) + g(2) * (-t53 * t99 + t96) + g(3) * t71 + (-g(3) * t53 + t120 * t82) * t49) (-m(5) - m(6) - m(7)) * (-g(3) * t50 + t119) -m(6) * (g(1) * (rSges(6,1) * t12 - rSges(6,2) * t13) + g(2) * (-rSges(6,1) * t10 + rSges(6,2) * t11)) - m(7) * (g(1) * (pkin(5) * t12 + t117) + g(2) * (-pkin(5) * t10 + t118)) + (-m(6) * (-rSges(6,1) * t46 - rSges(6,2) * t47) - m(7) * (t76 - t114)) * t110, -m(7) * (g(1) * t117 + g(2) * t118 + t76 * t110)];
taug  = t1(:);
