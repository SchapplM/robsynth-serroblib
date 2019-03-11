% Calculate Gravitation load on the joints for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
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
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 15:38:22
% EndTime: 2019-03-09 15:38:25
% DurationCPUTime: 1.52s
% Computational Cost: add. (791->196), mult. (1306->279), div. (0->0), fcn. (1527->14), ass. (0->80)
t54 = sin(qJ(3));
t57 = cos(qJ(3));
t87 = cos(pkin(6));
t50 = sin(pkin(6));
t55 = sin(qJ(2));
t97 = t50 * t55;
t121 = -t54 * t97 + t87 * t57;
t105 = cos(qJ(2));
t106 = cos(qJ(1));
t56 = sin(qJ(1));
t78 = t55 * t87;
t26 = t106 * t105 - t56 * t78;
t95 = t50 * t57;
t11 = -t26 * t54 + t56 * t95;
t70 = t87 * t105;
t23 = -t106 * t70 + t55 * t56;
t47 = pkin(12) + qJ(6);
t42 = sin(t47);
t44 = cos(t47);
t24 = t56 * t105 + t106 * t78;
t48 = qJ(3) + pkin(11);
t43 = sin(t48);
t45 = cos(t48);
t82 = t50 * t106;
t6 = t24 * t45 - t43 * t82;
t120 = -t23 * t44 + t42 * t6;
t119 = -t23 * t42 - t44 * t6;
t62 = t24 * t54 + t57 * t82;
t118 = t62 * pkin(3);
t89 = pkin(10) + qJ(5) + rSges(7,3);
t88 = qJ(5) + rSges(6,3);
t17 = t43 * t97 - t87 * t45;
t5 = t24 * t43 + t45 * t82;
t96 = t50 * t56;
t9 = t26 * t43 - t45 * t96;
t116 = g(1) * t9 + g(2) * t5 + g(3) * t17;
t25 = t106 * t55 + t56 * t70;
t115 = g(1) * t25 + g(2) * t23;
t112 = -m(6) - m(7);
t109 = g(3) * t50;
t49 = sin(pkin(12));
t108 = t49 * pkin(5);
t107 = -pkin(9) - rSges(4,3);
t104 = rSges(5,2) * t43;
t101 = t23 * t49;
t99 = t25 * t49;
t41 = pkin(3) * t57 + pkin(2);
t52 = -qJ(4) - pkin(9);
t92 = -t23 * t41 - t24 * t52;
t91 = -t25 * t41 - t26 * t52;
t90 = t106 * pkin(1) + pkin(8) * t96;
t84 = t54 * t96;
t83 = t45 * t105;
t81 = t50 * t105;
t80 = -t56 * pkin(1) + pkin(8) * t82;
t34 = t54 * t82;
t79 = -t24 * t57 + t34;
t76 = t43 * t81;
t75 = t45 * t81;
t73 = t11 * pkin(3);
t27 = t41 * t81;
t72 = -t52 * t97 + t27;
t71 = pkin(3) * t84 - t25 * t52 + t26 * t41 + t90;
t69 = -rSges(5,1) * t45 + t104;
t51 = cos(pkin(12));
t68 = t49 * rSges(6,1) + t51 * rSges(6,2);
t67 = t121 * pkin(3);
t65 = -t51 * rSges(6,1) + t49 * rSges(6,2) - pkin(4);
t40 = pkin(5) * t51 + pkin(4);
t64 = -t44 * rSges(7,1) + t42 * rSges(7,2) - t40;
t63 = pkin(3) * t34 + t23 * t52 - t24 * t41 + t80;
t61 = t42 * rSges(7,1) + t44 * rSges(7,2) + t108;
t59 = -t88 * t43 + t65 * t45;
t58 = -t89 * t43 + t64 * t45;
t18 = t87 * t43 + t45 * t97;
t12 = t26 * t57 + t84;
t10 = t26 * t45 + t43 * t96;
t3 = t10 * t44 + t25 * t42;
t2 = -t10 * t42 + t25 * t44;
t1 = [-m(2) * (g(1) * (-t56 * rSges(2,1) - t106 * rSges(2,2)) + g(2) * (t106 * rSges(2,1) - t56 * rSges(2,2))) - m(3) * (g(1) * (-t24 * rSges(3,1) + t23 * rSges(3,2) + rSges(3,3) * t82 + t80) + g(2) * (rSges(3,1) * t26 - rSges(3,2) * t25 + rSges(3,3) * t96 + t90)) - m(4) * (g(1) * (t79 * rSges(4,1) + t62 * rSges(4,2) - t24 * pkin(2) + t107 * t23 + t80) + g(2) * (rSges(4,1) * t12 + rSges(4,2) * t11 + pkin(2) * t26 - t107 * t25 + t90)) - m(5) * (g(1) * (-rSges(5,1) * t6 + rSges(5,2) * t5 - rSges(5,3) * t23 + t63) + g(2) * (rSges(5,1) * t10 - rSges(5,2) * t9 + rSges(5,3) * t25 + t71)) - m(6) * (g(1) * (-t6 * pkin(4) + (-t51 * t6 - t101) * rSges(6,1) + (-t23 * t51 + t49 * t6) * rSges(6,2) - t88 * t5 + t63) + g(2) * (t10 * pkin(4) + (t10 * t51 + t99) * rSges(6,1) + (-t10 * t49 + t25 * t51) * rSges(6,2) + t88 * t9 + t71)) - m(7) * (g(1) * (rSges(7,1) * t119 + rSges(7,2) * t120 - pkin(5) * t101 - t6 * t40 - t89 * t5 + t63) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t2 + pkin(5) * t99 + t10 * t40 + t89 * t9 + t71)) -m(3) * (g(1) * (-rSges(3,1) * t25 - rSges(3,2) * t26) + g(2) * (-rSges(3,1) * t23 - rSges(3,2) * t24) + (t105 * rSges(3,1) - rSges(3,2) * t55) * t109) - m(4) * ((-g(1) * t26 - g(2) * t24 - t55 * t109) * t107 + (t109 * t105 - t115) * (t57 * rSges(4,1) - t54 * rSges(4,2) + pkin(2))) - m(5) * (g(1) * (rSges(5,3) * t26 + t69 * t25 + t91) + g(2) * (rSges(5,3) * t24 + t69 * t23 + t92) + g(3) * t27 + (rSges(5,1) * t83 - t105 * t104 + (rSges(5,3) - t52) * t55) * t109) - m(6) * (g(1) * (t59 * t25 + t68 * t26 + t91) + g(2) * (t59 * t23 + t68 * t24 + t92) + g(3) * (pkin(4) * t75 + t72 + t88 * t76 + ((t49 * t55 + t51 * t83) * rSges(6,1) + (-t49 * t83 + t51 * t55) * rSges(6,2)) * t50)) - m(7) * (g(1) * (t58 * t25 + t61 * t26 + t91) + g(2) * (t58 * t23 + t61 * t24 + t92) + g(3) * (t97 * t108 + t40 * t75 + t72 + t89 * t76 + ((t42 * t55 + t44 * t83) * rSges(7,1) + (-t42 * t83 + t44 * t55) * rSges(7,2)) * t50)) -m(4) * (g(1) * (rSges(4,1) * t11 - rSges(4,2) * t12) + g(2) * (-t62 * rSges(4,1) + t79 * rSges(4,2)) + g(3) * (t121 * rSges(4,1) + (-t87 * t54 - t55 * t95) * rSges(4,2))) - m(5) * (g(1) * (-rSges(5,1) * t9 - rSges(5,2) * t10 + t73) + g(2) * (-t5 * rSges(5,1) - t6 * rSges(5,2) - t118) + g(3) * (-rSges(5,1) * t17 - rSges(5,2) * t18 + t67)) + (-g(1) * (t89 * t10 + t73) - g(2) * (t89 * t6 - t118) - g(3) * (t89 * t18 + t67) - t116 * t64) * m(7) + (-g(1) * (t88 * t10 + t73) - g(2) * (t88 * t6 - t118) - g(3) * (t88 * t18 + t67) - t116 * t65) * m(6) (-m(5) + t112) * (-g(3) * t81 + t115) t112 * t116, -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * (-rSges(7,1) * t120 + rSges(7,2) * t119) + g(3) * ((-t18 * t42 - t44 * t81) * rSges(7,1) + (-t18 * t44 + t42 * t81) * rSges(7,2)))];
taug  = t1(:);
