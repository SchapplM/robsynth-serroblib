% Calculate Gravitation load on the joints for
% S6RRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
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
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRRPP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRPP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRRPP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:44:03
% EndTime: 2019-03-09 20:44:05
% DurationCPUTime: 1.06s
% Computational Cost: add. (612->165), mult. (637->215), div. (0->0), fcn. (602->10), ass. (0->85)
t126 = rSges(5,3) + pkin(9);
t125 = rSges(7,1) + pkin(5);
t53 = sin(qJ(4));
t58 = cos(qJ(1));
t101 = t53 * t58;
t51 = qJ(2) + qJ(3);
t48 = cos(t51);
t55 = sin(qJ(1));
t56 = cos(qJ(4));
t99 = t55 * t56;
t8 = -t48 * t101 + t99;
t93 = rSges(7,3) + qJ(6);
t47 = sin(t51);
t124 = t48 * rSges(4,1) - rSges(4,2) * t47;
t123 = g(1) * t58 + g(2) * t55;
t67 = t48 * pkin(3) + t126 * t47;
t122 = t123 * t47;
t54 = sin(qJ(2));
t121 = pkin(2) * t54;
t120 = pkin(4) * t53;
t59 = -pkin(8) - pkin(7);
t117 = g(2) * t59;
t52 = -qJ(5) - pkin(9);
t116 = g(3) * t52;
t114 = rSges(3,3) + pkin(7);
t113 = rSges(5,1) * t56;
t111 = rSges(5,2) * t53;
t50 = qJ(4) + pkin(10);
t45 = sin(t50);
t110 = t45 * t47;
t109 = t45 * t48;
t46 = cos(t50);
t108 = t46 * t48;
t39 = t47 * rSges(7,2);
t37 = t47 * rSges(6,3);
t107 = t47 * t58;
t43 = pkin(4) * t56 + pkin(3);
t20 = t48 * t43;
t106 = t48 * t55;
t105 = t48 * t58;
t104 = t52 * t55;
t103 = t52 * t58;
t102 = t53 * t55;
t100 = t55 * t46;
t98 = t56 * t58;
t97 = t58 * t45;
t96 = rSges(4,3) - t59;
t89 = rSges(6,2) * t110;
t95 = rSges(6,3) * t106 + t55 * t89;
t94 = rSges(6,3) * t105 + t58 * t89;
t92 = t55 * t121;
t91 = t58 * t121;
t90 = t47 * t111;
t88 = t48 * t104;
t87 = t48 * t103;
t85 = t126 * t106 + t55 * t90;
t84 = t126 * t105 + t58 * t90;
t82 = -rSges(6,1) * t46 - t43;
t57 = cos(qJ(2));
t49 = t57 * pkin(2);
t44 = t49 + pkin(1);
t80 = -t44 - t20;
t79 = rSges(7,2) * t106 - t88;
t78 = rSges(7,2) * t105 - t87;
t77 = pkin(4) * t101 + t47 * t104 - t58 * t59;
t76 = rSges(3,1) * t57 - rSges(3,2) * t54;
t73 = -rSges(4,1) * t47 - rSges(4,2) * t48;
t72 = t8 * pkin(4);
t71 = t125 * t108 + t93 * t109 + t20 + t39;
t70 = pkin(1) + t76;
t6 = t48 * t102 + t98;
t69 = rSges(6,1) * t108 - rSges(6,2) * t109 + t20 + t37;
t22 = t58 * t44;
t68 = pkin(4) * t102 - t47 * t103 + t43 * t105 + t22;
t66 = t67 + (-t111 + t113) * t48;
t65 = t6 * pkin(4);
t61 = (-pkin(3) - t113) * t122;
t60 = (-t116 + t123 * (-t125 * t46 - t93 * t45 - t43)) * t47;
t9 = t48 * t98 + t102;
t7 = -t48 * t99 + t101;
t5 = t46 * t105 + t45 * t55;
t4 = t48 * t97 - t100;
t3 = t48 * t100 - t97;
t2 = t45 * t106 + t46 * t58;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t55 - rSges(2,2) * t58) + g(2) * (rSges(2,1) * t58 - rSges(2,2) * t55)) - m(3) * ((g(1) * t114 + g(2) * t70) * t58 + (-g(1) * t70 + g(2) * t114) * t55) - m(4) * (g(2) * t22 + (g(1) * t96 + g(2) * t124) * t58 + (g(1) * (-t44 - t124) + g(2) * t96) * t55) - m(5) * (g(1) * (t7 * rSges(5,1) + t6 * rSges(5,2)) + g(2) * (t9 * rSges(5,1) + t8 * rSges(5,2) + t22) + (-g(1) * t59 + g(2) * t67) * t58 + (g(1) * (-t44 - t67) - t117) * t55) - m(6) * (g(1) * (-t3 * rSges(6,1) + t2 * rSges(6,2) + t77) + g(2) * (t5 * rSges(6,1) - t4 * rSges(6,2) + rSges(6,3) * t107 + t68) + (g(1) * (t80 - t37) - t117) * t55) - m(7) * (g(1) * (-t125 * t3 - t93 * t2 + t77) + g(2) * (rSges(7,2) * t107 + t125 * t5 + t93 * t4 + t68) + (g(1) * (t80 - t39) - t117) * t55) -m(3) * (g(3) * t76 + t123 * (-rSges(3,1) * t54 - rSges(3,2) * t57)) - m(4) * (g(3) * (t49 + t124) + t123 * (t73 - t121)) - m(5) * (g(1) * (t84 - t91) + g(2) * (t85 - t92) + g(3) * (t49 + t66) + t61) - m(6) * (g(1) * t94 + g(2) * t95 + g(3) * (-t47 * t52 + t49 + t69) + t123 * (t82 * t47 - t48 * t52 - t121)) - m(7) * (g(1) * (t78 - t91) + g(2) * (t79 - t92) + g(3) * (t49 + t71) + t60) -m(4) * (g(3) * t124 + t123 * t73) - m(5) * (g(1) * t84 + g(2) * t85 + g(3) * t66 + t61) - m(6) * (g(1) * (-t87 + t94) + g(2) * (-t88 + t95) + g(3) * t69 + (t123 * t82 - t116) * t47) - m(7) * (g(1) * t78 + g(2) * t79 + g(3) * t71 + t60) -m(5) * (g(1) * (rSges(5,1) * t8 - rSges(5,2) * t9) + g(2) * (-rSges(5,1) * t6 + rSges(5,2) * t7)) - m(6) * (g(1) * (-rSges(6,1) * t4 - rSges(6,2) * t5 + t72) + g(2) * (-rSges(6,1) * t2 - rSges(6,2) * t3 - t65)) - m(7) * (g(1) * (-t125 * t4 + t93 * t5 + t72) + g(2) * (-t125 * t2 + t93 * t3 - t65)) + (-m(5) * (-rSges(5,1) * t53 - rSges(5,2) * t56) - m(6) * (-rSges(6,1) * t45 - rSges(6,2) * t46 - t120) - m(7) * (-t125 * t45 + t93 * t46 - t120)) * g(3) * t47 (-m(6) - m(7)) * (-g(3) * t48 + t122) -m(7) * (g(1) * t4 + g(2) * t2 + g(3) * t110)];
taug  = t1(:);
