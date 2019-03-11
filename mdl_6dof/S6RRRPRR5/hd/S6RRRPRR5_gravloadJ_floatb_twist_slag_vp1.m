% Calculate Gravitation load on the joints for
% S6RRRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
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
% Datum: 2019-03-09 18:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:20:44
% EndTime: 2019-03-09 18:20:46
% DurationCPUTime: 0.83s
% Computational Cost: add. (503->156), mult. (552->205), div. (0->0), fcn. (503->10), ass. (0->88)
t121 = rSges(6,3) + pkin(9);
t56 = qJ(2) + qJ(3);
t53 = cos(t56);
t132 = t121 * t53;
t55 = qJ(5) + qJ(6);
t52 = cos(t55);
t112 = t52 * t53;
t50 = sin(t55);
t117 = t50 * t53;
t51 = sin(t56);
t63 = -pkin(10) - pkin(9);
t131 = rSges(7,1) * t117 + rSges(7,2) * t112 + t51 * t63;
t64 = -pkin(8) - pkin(7);
t119 = pkin(4) - t64;
t73 = -t53 * rSges(5,2) + t51 * rSges(5,3);
t130 = t53 * rSges(4,1) - rSges(4,2) * t51;
t59 = sin(qJ(1));
t124 = g(2) * t59;
t62 = cos(qJ(1));
t129 = g(1) * t62 + t124;
t57 = sin(qJ(5));
t116 = t51 * t57;
t39 = pkin(5) * t116;
t128 = t39 + (rSges(7,3) - t63) * t53;
t127 = t129 * t51;
t58 = sin(qJ(2));
t126 = pkin(2) * t58;
t123 = g(3) * t53;
t47 = t53 * pkin(3);
t122 = rSges(3,3) + pkin(7);
t120 = -rSges(7,3) - pkin(3);
t115 = t51 * t59;
t114 = t51 * t62;
t111 = t52 * t59;
t110 = t52 * t62;
t60 = cos(qJ(5));
t108 = t53 * t60;
t107 = t53 * t62;
t105 = t57 * t59;
t104 = t57 * t62;
t103 = t59 * t60;
t102 = t60 * t62;
t101 = rSges(5,1) - t64;
t100 = rSges(4,3) - t64;
t41 = t51 * qJ(4);
t99 = t41 + t47;
t98 = pkin(5) * t60 + t119;
t97 = qJ(4) * t53;
t96 = -pkin(3) - t121;
t95 = t59 * t126;
t94 = t62 * t126;
t92 = rSges(6,2) * t108;
t90 = t53 * t105;
t89 = t53 * t104;
t29 = t59 * t97;
t88 = rSges(6,1) * t90 + t59 * t92 + t29;
t31 = t62 * t97;
t87 = rSges(6,1) * t89 + t62 * t92 + t31;
t86 = t59 * t53 * rSges(5,3) + rSges(5,2) * t115 + t29;
t61 = cos(qJ(2));
t54 = t61 * pkin(2);
t49 = t54 + pkin(1);
t32 = t62 * t49;
t85 = pkin(3) * t107 + t62 * t41 + t32;
t84 = rSges(5,2) * t114 + rSges(5,3) * t107 + t31;
t82 = -t49 - t41;
t81 = g(1) * t96;
t80 = pkin(5) * t90 + t131 * t59 + t29;
t79 = pkin(5) * t89 + t131 * t62 + t31;
t78 = -pkin(3) * t51 - t126;
t77 = rSges(3,1) * t61 - rSges(3,2) * t58;
t74 = -rSges(4,1) * t51 - rSges(4,2) * t53;
t72 = t99 + t73;
t71 = t51 * t60 * rSges(6,2) + rSges(6,1) * t116 + t132 + t99;
t70 = pkin(1) + t77;
t10 = t51 * t102 - t105;
t12 = t51 * t103 + t104;
t68 = t128 + t99 + (rSges(7,1) * t50 + rSges(7,2) * t52) * t51;
t67 = t120 * t127;
t66 = (t96 * t124 + t62 * t81) * t51;
t6 = t51 * t110 - t50 * t59;
t7 = t50 * t114 + t111;
t8 = t51 * t111 + t50 * t62;
t9 = -t50 * t115 + t110;
t65 = g(1) * (t6 * rSges(7,1) - t7 * rSges(7,2)) + g(2) * (t8 * rSges(7,1) + t9 * rSges(7,2)) + g(3) * (-rSges(7,1) * t112 + rSges(7,2) * t117);
t13 = -t51 * t105 + t102;
t11 = t51 * t104 + t103;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t59 - rSges(2,2) * t62) + g(2) * (rSges(2,1) * t62 - rSges(2,2) * t59)) - m(3) * ((g(1) * t122 + g(2) * t70) * t62 + (-g(1) * t70 + g(2) * t122) * t59) - m(4) * (g(2) * t32 + (g(1) * t100 + g(2) * t130) * t62 + (g(1) * (-t49 - t130) + g(2) * t100) * t59) - m(5) * (g(2) * t85 + (g(1) * t101 + g(2) * t73) * t62 + (g(1) * (-t73 + t82 - t47) + g(2) * t101) * t59) - m(6) * (g(1) * (t13 * rSges(6,1) - t12 * rSges(6,2)) + g(2) * (t11 * rSges(6,1) + t10 * rSges(6,2) + t85) + (g(1) * t119 + g(2) * t132) * t62 + (g(1) * t82 + g(2) * t119 + t53 * t81) * t59) - m(7) * (g(1) * (t9 * rSges(7,1) - t8 * rSges(7,2)) + g(2) * (t7 * rSges(7,1) + t6 * rSges(7,2) + t85) + (g(1) * t98 + g(2) * t128) * t62 + (g(2) * t98 + (t82 - t39 + (t63 + t120) * t53) * g(1)) * t59) -m(3) * (g(3) * t77 + t129 * (-rSges(3,1) * t58 - rSges(3,2) * t61)) - m(4) * (g(3) * (t54 + t130) + t129 * (t74 - t126)) - m(5) * (g(1) * (t78 * t62 + t84) + g(2) * (t78 * t59 + t86) + g(3) * (t54 + t72)) - m(6) * (g(1) * (t87 - t94) + g(2) * (t88 - t95) + g(3) * (t54 + t71) + t66) - m(7) * (g(1) * (t79 - t94) + g(2) * (t80 - t95) + g(3) * (t54 + t68) + t67) -m(4) * (g(3) * t130 + t129 * t74) - m(5) * (g(1) * (-pkin(3) * t114 + t84) + g(2) * (-pkin(3) * t115 + t86) + g(3) * t72) - m(6) * (g(1) * t87 + g(2) * t88 + g(3) * t71 + t66) - m(7) * (g(1) * t79 + g(2) * t80 + g(3) * t68 + t67) (-m(5) - m(6) - m(7)) * (-t123 + t127) -m(6) * (g(1) * (rSges(6,1) * t10 - rSges(6,2) * t11) + g(2) * (rSges(6,1) * t12 + rSges(6,2) * t13) + (-rSges(6,1) * t60 + rSges(6,2) * t57) * t123) - m(7) * ((g(1) * t10 + g(2) * t12 - g(3) * t108) * pkin(5) + t65) -m(7) * t65];
taug  = t1(:);
