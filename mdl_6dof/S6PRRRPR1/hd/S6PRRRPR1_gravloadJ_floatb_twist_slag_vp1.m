% Calculate Gravitation load on the joints for
% S6PRRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:00:13
% EndTime: 2019-03-08 23:00:15
% DurationCPUTime: 0.81s
% Computational Cost: add. (703->146), mult. (945->210), div. (0->0), fcn. (1077->14), ass. (0->77)
t127 = rSges(7,3) + pkin(10);
t62 = sin(pkin(6));
t65 = sin(qJ(2));
t103 = t62 * t65;
t60 = qJ(3) + qJ(4);
t56 = sin(t60);
t57 = cos(t60);
t94 = cos(pkin(6));
t126 = -t56 * t103 + t94 * t57;
t61 = sin(pkin(11));
t104 = t61 * t62;
t68 = cos(qJ(2));
t90 = t61 * t94;
t93 = cos(pkin(11));
t41 = -t65 * t90 + t68 * t93;
t125 = t57 * t104 - t41 * t56;
t63 = sin(qJ(6));
t66 = cos(qJ(6));
t124 = rSges(7,1) * t66 - rSges(7,2) * t63 + pkin(5);
t101 = t62 * t68;
t123 = g(3) * t101;
t122 = rSges(7,1) * t63 + rSges(7,2) * t66;
t84 = t94 * t93;
t38 = t61 * t65 - t68 * t84;
t115 = g(2) * t38;
t40 = t65 * t93 + t68 * t90;
t121 = g(1) * t40 + t115;
t55 = pkin(12) + t60;
t51 = sin(t55);
t52 = cos(t55);
t120 = t52 * t124 + t127 * t51;
t119 = -m(6) - m(7);
t69 = -pkin(9) - pkin(8);
t39 = t61 * t68 + t65 * t84;
t89 = t62 * t93;
t12 = -t39 * t51 - t52 * t89;
t13 = t39 * t52 - t51 * t89;
t118 = t12 * rSges(6,1) - t13 * rSges(6,2);
t14 = t104 * t52 - t41 * t51;
t15 = t104 * t51 + t41 * t52;
t117 = t14 * rSges(6,1) - t15 * rSges(6,2);
t114 = g(2) * t39;
t67 = cos(qJ(3));
t58 = t67 * pkin(3);
t47 = pkin(4) * t57 + t58;
t45 = pkin(2) + t47;
t113 = t45 * t123;
t112 = g(3) * t62;
t111 = rSges(4,3) + pkin(8);
t102 = t62 * t67;
t100 = rSges(5,3) - t69;
t59 = -qJ(5) + t69;
t99 = -t38 * t45 - t39 * t59;
t98 = -t40 * t45 - t41 * t59;
t64 = sin(qJ(3));
t46 = -pkin(3) * t64 - pkin(4) * t56;
t97 = t47 * t104 + t41 * t46;
t29 = -t103 * t51 + t52 * t94;
t30 = t103 * t52 + t51 * t94;
t96 = t29 * rSges(6,1) - t30 * rSges(6,2);
t95 = t46 * t103 + t94 * t47;
t87 = t125 * pkin(4);
t86 = rSges(6,1) * t52 - rSges(6,2) * t51;
t85 = t126 * pkin(4);
t83 = rSges(4,1) * t67 - rSges(4,2) * t64 + pkin(2);
t81 = t102 * t61 - t41 * t64;
t80 = rSges(5,1) * t57 - rSges(5,2) * t56 + pkin(2) + t58;
t79 = t39 * t46 - t47 * t89;
t78 = -t39 * t56 - t57 * t89;
t77 = -t39 * t64 - t67 * t89;
t76 = -t103 * t64 + t67 * t94;
t75 = t124 * t12 + t127 * t13;
t74 = t124 * t14 + t127 * t15;
t73 = t124 * t29 + t127 * t30;
t72 = t78 * pkin(4);
t71 = g(1) * (t125 * rSges(5,1) + (-t104 * t56 - t41 * t57) * rSges(5,2)) + g(2) * (t78 * rSges(5,1) + (-t39 * t57 + t56 * t89) * rSges(5,2)) + g(3) * (t126 * rSges(5,1) + (-t103 * t57 - t56 * t94) * rSges(5,2));
t1 = [(-m(2) - m(3) - m(4) - m(5) + t119) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t40 - rSges(3,2) * t41) + g(2) * (-rSges(3,1) * t38 - rSges(3,2) * t39) + (rSges(3,1) * t68 - rSges(3,2) * t65) * t112) - m(4) * (g(1) * (t111 * t41 - t40 * t83) + t111 * t114 - t83 * t115 + (t111 * t65 + t68 * t83) * t112) - m(5) * (g(1) * (t100 * t41 - t40 * t80) + t100 * t114 - t80 * t115 + (t100 * t65 + t68 * t80) * t112) - m(6) * (g(1) * (rSges(6,3) * t41 - t40 * t86 + t98) + g(2) * (t39 * rSges(6,3) - t38 * t86 + t99) + t113 + (t86 * t68 + (rSges(6,3) - t59) * t65) * t112) - m(7) * (g(1) * (t122 * t41 + t98) + g(2) * (t122 * t39 + t99) + t113 + ((-t59 + t122) * t65 + t120 * t68) * t112 - t121 * t120) -m(4) * (g(1) * (t81 * rSges(4,1) + (-t104 * t64 - t41 * t67) * rSges(4,2)) + g(2) * (t77 * rSges(4,1) + (-t39 * t67 + t64 * t89) * rSges(4,2)) + g(3) * (t76 * rSges(4,1) + (-t102 * t65 - t64 * t94) * rSges(4,2))) - m(5) * ((g(1) * t81 + g(2) * t77 + g(3) * t76) * pkin(3) + t71) - m(6) * (g(1) * (t97 + t117) + g(2) * (t79 + t118) + g(3) * (t95 + t96)) - m(7) * (g(1) * (t74 + t97) + g(2) * (t75 + t79) + g(3) * (t73 + t95)) -m(5) * t71 - m(6) * (g(1) * (t87 + t117) + g(2) * (t72 + t118) + g(3) * (t85 + t96)) - m(7) * (g(1) * (t74 + t87) + g(2) * (t72 + t75) + g(3) * (t73 + t85)) t119 * (t121 - t123) -m(7) * (g(1) * ((-t15 * t63 + t40 * t66) * rSges(7,1) + (-t15 * t66 - t40 * t63) * rSges(7,2)) + g(2) * ((-t13 * t63 + t38 * t66) * rSges(7,1) + (-t13 * t66 - t38 * t63) * rSges(7,2)) + g(3) * ((-t101 * t66 - t30 * t63) * rSges(7,1) + (t101 * t63 - t30 * t66) * rSges(7,2)))];
taug  = t1(:);
