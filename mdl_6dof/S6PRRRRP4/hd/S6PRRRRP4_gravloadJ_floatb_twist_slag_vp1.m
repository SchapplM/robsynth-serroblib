% Calculate Gravitation load on the joints for
% S6PRRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-03-09 00:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:13:09
% EndTime: 2019-03-09 00:13:11
% DurationCPUTime: 1.11s
% Computational Cost: add. (730->163), mult. (1461->239), div. (0->0), fcn. (1761->12), ass. (0->87)
t136 = rSges(7,1) + pkin(5);
t71 = sin(qJ(4));
t128 = pkin(4) * t71;
t139 = pkin(8) + t128;
t135 = rSges(7,3) + qJ(6);
t101 = cos(pkin(11));
t73 = sin(qJ(2));
t76 = cos(qJ(2));
t102 = cos(pkin(6));
t69 = sin(pkin(11));
t97 = t69 * t102;
t53 = t101 * t76 - t73 * t97;
t70 = sin(pkin(6));
t72 = sin(qJ(3));
t75 = cos(qJ(3));
t33 = t69 * t70 * t72 + t53 * t75;
t52 = t101 * t73 + t76 * t97;
t74 = cos(qJ(4));
t138 = -t33 * t71 + t52 * t74;
t85 = t102 * t101;
t51 = t69 * t76 + t73 * t85;
t96 = t70 * t101;
t31 = t51 * t75 - t72 * t96;
t50 = t69 * t73 - t76 * t85;
t137 = -t31 * t71 + t50 * t74;
t134 = g(1) * t52 + g(2) * t50;
t30 = -t51 * t72 - t75 * t96;
t112 = t70 * t75;
t32 = t69 * t112 - t53 * t72;
t113 = t70 * t73;
t54 = t102 * t75 - t72 * t113;
t133 = g(1) * t32 + g(2) * t30 + g(3) * t54;
t111 = t70 * t76;
t132 = (g(3) * t111 - t134) * t72;
t121 = rSges(5,3) + pkin(9);
t83 = rSges(5,1) * t74 - rSges(5,2) * t71 + pkin(3);
t131 = t121 * t72 + t83 * t75;
t68 = qJ(4) + qJ(5);
t66 = sin(t68);
t67 = cos(t68);
t11 = t31 * t66 - t50 * t67;
t12 = t31 * t67 + t50 * t66;
t130 = -t11 * rSges(6,1) - t12 * rSges(6,2);
t13 = t33 * t66 - t52 * t67;
t14 = t33 * t67 + t52 * t66;
t129 = -t13 * rSges(6,1) - t14 * rSges(6,2);
t123 = g(3) * t70;
t122 = rSges(4,3) + pkin(8);
t65 = pkin(4) * t74 + pkin(3);
t116 = t65 * t75;
t115 = t66 * t75;
t114 = t67 * t75;
t110 = t75 * t76;
t77 = -pkin(10) - pkin(9);
t107 = t30 * t65 - t31 * t77;
t106 = t32 * t65 - t33 * t77;
t55 = t102 * t72 + t73 * t112;
t26 = t67 * t111 + t55 * t66;
t27 = -t66 * t111 + t55 * t67;
t105 = -t26 * rSges(6,1) - t27 * rSges(6,2);
t104 = t54 * t65 - t55 * t77;
t103 = pkin(2) * t111 + pkin(8) * t113;
t99 = t70 * t110;
t98 = g(3) * t103;
t95 = t113 * t128 + t65 * t99 + t103;
t94 = t137 * pkin(4);
t93 = t138 * pkin(4);
t92 = rSges(4,1) * t75 - rSges(4,2) * t72;
t91 = t71 * rSges(5,1) + t74 * rSges(5,2);
t90 = rSges(6,1) * t67 - rSges(6,2) * t66;
t89 = -t136 * t11 + t135 * t12;
t48 = t50 * pkin(2);
t88 = -t50 * t116 + t139 * t51 - t48;
t49 = t52 * pkin(2);
t87 = -t52 * t116 + t139 * t53 - t49;
t86 = -t136 * t13 + t135 * t14;
t84 = t135 * t27 - t136 * t26;
t82 = pkin(8) + t91;
t81 = -t74 * t111 - t55 * t71;
t80 = t81 * pkin(4);
t37 = (t67 * t110 + t66 * t73) * t70;
t36 = -t67 * t113 + t66 * t99;
t18 = -t52 * t114 + t53 * t66;
t17 = -t52 * t115 - t53 * t67;
t16 = -t50 * t114 + t51 * t66;
t15 = -t50 * t115 - t51 * t67;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t52 - rSges(3,2) * t53) + g(2) * (-rSges(3,1) * t50 - rSges(3,2) * t51) + (rSges(3,1) * t76 - rSges(3,2) * t73) * t123) - m(4) * (g(1) * (t122 * t53 - t92 * t52 - t49) + g(2) * (t122 * t51 - t92 * t50 - t48) + t98 + (rSges(4,3) * t73 + t92 * t76) * t123) - m(5) * (t98 + (t131 * t76 + t91 * t73) * t123 - t134 * t131 + (t82 * t51 - t48) * g(2) + (t82 * t53 - t49) * g(1)) - m(6) * (g(1) * (t18 * rSges(6,1) - t17 * rSges(6,2) + t87) + g(2) * (t16 * rSges(6,1) - t15 * rSges(6,2) + t88) + g(3) * (t37 * rSges(6,1) - t36 * rSges(6,2) + t95) + (rSges(6,3) - t77) * t132) - m(7) * (g(1) * (t135 * t17 + t136 * t18 + t87) + g(2) * (t135 * t15 + t136 * t16 + t88) + g(3) * (t135 * t36 + t136 * t37 + t95) + (rSges(7,2) - t77) * t132) -m(4) * (g(1) * (rSges(4,1) * t32 - rSges(4,2) * t33) + g(2) * (rSges(4,1) * t30 - rSges(4,2) * t31) + g(3) * (rSges(4,1) * t54 - rSges(4,2) * t55)) - m(5) * (t133 * t83 + (g(1) * t33 + g(2) * t31 + g(3) * t55) * t121) - m(6) * (g(1) * (rSges(6,3) * t33 + t90 * t32 + t106) + g(2) * (rSges(6,3) * t31 + t90 * t30 + t107) + g(3) * (rSges(6,3) * t55 + t90 * t54 + t104)) + (-g(1) * (rSges(7,2) * t33 + t106) - g(2) * (rSges(7,2) * t31 + t107) - g(3) * (rSges(7,2) * t55 + t104) - t133 * (t135 * t66 + t136 * t67)) * m(7), -m(5) * (g(1) * (t138 * rSges(5,1) + (-t33 * t74 - t52 * t71) * rSges(5,2)) + g(2) * (t137 * rSges(5,1) + (-t31 * t74 - t50 * t71) * rSges(5,2)) + g(3) * (t81 * rSges(5,1) + (t71 * t111 - t55 * t74) * rSges(5,2))) - m(6) * (g(1) * (t93 + t129) + g(2) * (t94 + t130) + g(3) * (t80 + t105)) - m(7) * (g(1) * (t86 + t93) + g(2) * (t89 + t94) + g(3) * (t80 + t84)) -m(6) * (g(1) * t129 + g(2) * t130 + g(3) * t105) - m(7) * (g(1) * t86 + g(2) * t89 + g(3) * t84) -m(7) * (g(1) * t13 + g(2) * t11 + g(3) * t26)];
taug  = t1(:);
