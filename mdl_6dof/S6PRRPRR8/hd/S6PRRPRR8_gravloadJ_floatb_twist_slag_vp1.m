% Calculate Gravitation load on the joints for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:36:03
% EndTime: 2019-03-08 22:36:06
% DurationCPUTime: 1.08s
% Computational Cost: add. (866->182), mult. (2363->278), div. (0->0), fcn. (2973->14), ass. (0->82)
t113 = rSges(6,3) + pkin(10);
t98 = rSges(7,3) + pkin(11);
t61 = sin(qJ(2));
t64 = cos(qJ(2));
t57 = cos(pkin(12));
t89 = cos(pkin(6));
t80 = t57 * t89;
t87 = sin(pkin(12));
t43 = t61 * t80 + t64 * t87;
t74 = t89 * t87;
t45 = t57 * t64 - t61 * t74;
t112 = g(1) * t45 + g(2) * t43;
t42 = -t61 * t87 + t64 * t80;
t55 = sin(pkin(7));
t60 = sin(qJ(3));
t88 = cos(pkin(7));
t79 = t60 * t88;
t56 = sin(pkin(6));
t94 = t56 * t57;
t97 = cos(qJ(3));
t13 = -t55 * t60 * t94 + t42 * t79 + t43 * t97;
t44 = -t57 * t61 - t64 * t74;
t81 = t56 * t87;
t78 = t55 * t81;
t15 = t45 * t97 + (t44 * t88 + t78) * t60;
t82 = t55 * t89;
t28 = t60 * t82 + (t61 * t97 + t64 * t79) * t56;
t111 = g(1) * t15 + g(2) * t13 + g(3) * t28;
t75 = t88 * t97;
t83 = t56 * t97;
t12 = t55 * t57 * t83 - t42 * t75 + t43 * t60;
t14 = -t44 * t75 + t45 * t60 - t78 * t97;
t92 = t56 * t64;
t93 = t56 * t61;
t27 = t60 * t93 - t75 * t92 - t82 * t97;
t110 = g(1) * t14 + g(2) * t12 + g(3) * t27;
t59 = sin(qJ(5));
t96 = t55 * t59;
t63 = cos(qJ(5));
t95 = t55 * t63;
t84 = t55 * t93;
t91 = pkin(2) * t92 + pkin(9) * t84;
t90 = rSges(5,3) + qJ(4);
t86 = -m(5) - m(6) - m(7);
t85 = g(3) * t93;
t20 = t42 * t60 + t43 * t75;
t21 = t42 * t97 - t43 * t79;
t39 = t42 * pkin(2);
t77 = t21 * pkin(3) + t20 * qJ(4) + t39;
t22 = t44 * t60 + t45 * t75;
t23 = t44 * t97 - t45 * t79;
t40 = t44 * pkin(2);
t76 = t23 * pkin(3) + t22 * qJ(4) + t40;
t37 = (t60 * t64 + t61 * t75) * t56;
t38 = t64 * t83 - t79 * t93;
t73 = t38 * pkin(3) + t37 * qJ(4) + t91;
t58 = sin(qJ(6));
t62 = cos(qJ(6));
t72 = rSges(7,1) * t62 - rSges(7,2) * t58 + pkin(5);
t69 = pkin(4) * t84 + t73;
t68 = t21 * pkin(10) + t77;
t67 = t23 * pkin(10) + t76;
t66 = t112 * (pkin(4) + pkin(9)) * t55;
t41 = -t55 * t92 + t88 * t89;
t30 = -t44 * t55 + t81 * t88;
t29 = -t42 * t55 - t88 * t94;
t26 = t27 * pkin(3);
t25 = t37 * t59 + t63 * t84;
t24 = -t37 * t63 + t59 * t84;
t17 = t27 * t59 + t41 * t63;
t16 = t27 * t63 - t41 * t59;
t11 = t14 * pkin(3);
t10 = t12 * pkin(3);
t9 = t22 * t59 + t45 * t95;
t8 = t22 * t63 - t45 * t96;
t7 = t20 * t59 + t43 * t95;
t6 = t20 * t63 - t43 * t96;
t5 = t14 * t59 + t30 * t63;
t4 = t14 * t63 - t30 * t59;
t3 = t12 * t59 + t29 * t63;
t2 = t12 * t63 - t29 * t59;
t1 = [(-m(2) - m(3) - m(4) + t86) * g(3), -m(3) * (g(1) * (rSges(3,1) * t44 - rSges(3,2) * t45) + g(2) * (rSges(3,1) * t42 - rSges(3,2) * t43) + g(3) * (rSges(3,1) * t64 - rSges(3,2) * t61) * t56) - m(4) * (g(1) * (rSges(4,1) * t23 - rSges(4,2) * t22 + t40) + g(2) * (rSges(4,1) * t21 - rSges(4,2) * t20 + t39) + g(3) * (rSges(4,1) * t38 - rSges(4,2) * t37 + t91) + (rSges(4,3) * t85 + t112 * (rSges(4,3) + pkin(9))) * t55) - m(5) * (g(1) * (-rSges(5,2) * t23 + rSges(5,3) * t22 + t76) + g(2) * (-rSges(5,2) * t21 + rSges(5,3) * t20 + t77) + g(3) * (-rSges(5,2) * t38 + rSges(5,3) * t37 + t73) + (rSges(5,1) * t85 + t112 * (rSges(5,1) + pkin(9))) * t55) - m(6) * (g(1) * (rSges(6,1) * t9 + rSges(6,2) * t8 + rSges(6,3) * t23 + t67) + g(2) * (rSges(6,1) * t7 + rSges(6,2) * t6 + rSges(6,3) * t21 + t68) + g(3) * (rSges(6,1) * t25 - rSges(6,2) * t24 + t113 * t38 + t69) + t66) - m(7) * (g(1) * (t9 * pkin(5) + (t23 * t58 + t62 * t9) * rSges(7,1) + (t23 * t62 - t58 * t9) * rSges(7,2) + t67 - t98 * t8) + g(2) * (t7 * pkin(5) + (t21 * t58 + t62 * t7) * rSges(7,1) + (t21 * t62 - t58 * t7) * rSges(7,2) + t68 - t98 * t6) + g(3) * (t25 * pkin(5) + t38 * pkin(10) + (t25 * t62 + t38 * t58) * rSges(7,1) + (-t25 * t58 + t38 * t62) * rSges(7,2) + t98 * t24 + t69) + t66) -m(4) * (g(1) * (-rSges(4,1) * t14 - rSges(4,2) * t15) + g(2) * (-rSges(4,1) * t12 - rSges(4,2) * t13) + g(3) * (-rSges(4,1) * t27 - rSges(4,2) * t28)) - m(5) * (g(1) * (rSges(5,2) * t14 + t15 * t90 - t11) + g(2) * (rSges(5,2) * t12 + t13 * t90 - t10) + g(3) * (rSges(5,2) * t27 + t28 * t90 - t26)) + (-g(1) * (-t113 * t14 - t11) - g(2) * (-t113 * t12 - t10) - g(3) * (-t113 * t27 - t26) - t111 * (rSges(6,1) * t59 + rSges(6,2) * t63 + qJ(4))) * m(6) + (g(1) * t11 + g(2) * t10 + g(3) * t26 - t110 * (-rSges(7,1) * t58 - rSges(7,2) * t62 - pkin(10)) - t111 * (t59 * t72 - t63 * t98 + qJ(4))) * m(7), t86 * t110, -m(6) * (g(1) * (rSges(6,1) * t4 - rSges(6,2) * t5) + g(2) * (rSges(6,1) * t2 - rSges(6,2) * t3) + g(3) * (rSges(6,1) * t16 - rSges(6,2) * t17)) - m(7) * (g(1) * (t4 * t72 + t5 * t98) + (t16 * t72 + t17 * t98) * g(3) + (t2 * t72 + t3 * t98) * g(2)) -m(7) * (g(1) * ((t15 * t62 - t5 * t58) * rSges(7,1) + (-t15 * t58 - t5 * t62) * rSges(7,2)) + g(2) * ((t13 * t62 - t3 * t58) * rSges(7,1) + (-t13 * t58 - t3 * t62) * rSges(7,2)) + g(3) * ((-t17 * t58 + t28 * t62) * rSges(7,1) + (-t17 * t62 - t28 * t58) * rSges(7,2)))];
taug  = t1(:);
