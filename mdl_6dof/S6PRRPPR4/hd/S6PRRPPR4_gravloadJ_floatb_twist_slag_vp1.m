% Calculate Gravitation load on the joints for
% S6PRRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:12:33
% EndTime: 2019-03-08 21:12:36
% DurationCPUTime: 0.97s
% Computational Cost: add. (548->163), mult. (1402->245), div. (0->0), fcn. (1706->12), ass. (0->77)
t101 = cos(qJ(3));
t59 = sin(qJ(3));
t112 = -pkin(3) * t101 - qJ(4) * t59;
t55 = sin(pkin(11));
t57 = cos(pkin(11));
t111 = -pkin(4) * t57 - qJ(5) * t55;
t60 = sin(qJ(2));
t62 = cos(qJ(2));
t89 = cos(pkin(10));
t90 = cos(pkin(6));
t68 = t90 * t89;
t88 = sin(pkin(10));
t40 = t60 * t68 + t62 * t88;
t56 = sin(pkin(6));
t77 = t56 * t89;
t21 = t101 * t77 + t40 * t59;
t67 = t90 * t88;
t42 = -t60 * t67 + t62 * t89;
t76 = t56 * t88;
t23 = -t101 * t76 + t42 * t59;
t98 = t56 * t60;
t43 = -t101 * t90 + t59 * t98;
t110 = g(1) * t23 + g(2) * t21 + g(3) * t43;
t109 = -m(6) - m(7);
t104 = g(3) * t56;
t103 = rSges(4,3) + pkin(8);
t102 = rSges(7,3) + pkin(9);
t39 = t60 * t88 - t62 * t68;
t100 = t39 * t59;
t41 = t60 * t89 + t62 * t67;
t99 = t41 * t59;
t97 = t56 * t62;
t96 = pkin(2) * t97 + pkin(8) * t98;
t93 = rSges(6,2) + qJ(4);
t92 = rSges(5,3) + qJ(4);
t91 = rSges(6,3) + qJ(5);
t87 = -m(5) + t109;
t86 = t59 * t97;
t17 = t21 * pkin(3);
t84 = t111 * t21 - t17;
t18 = t23 * pkin(3);
t83 = t111 * t23 - t18;
t38 = t43 * pkin(3);
t82 = t111 * t43 - t38;
t80 = t55 * t101;
t79 = t57 * t101;
t78 = t62 * t101;
t74 = t56 * t78;
t75 = pkin(3) * t74 + qJ(4) * t86 + t96;
t27 = (t55 * t60 + t57 * t78) * t56;
t73 = pkin(4) * t27 + t75;
t72 = -rSges(5,1) * t57 + rSges(5,2) * t55;
t71 = -rSges(6,1) * t57 - rSges(6,3) * t55;
t36 = t39 * pkin(2);
t70 = t40 * pkin(8) + t112 * t39 - t36;
t37 = t41 * pkin(2);
t69 = t42 * pkin(8) + t112 * t41 - t37;
t10 = -t39 * t79 + t40 * t55;
t66 = pkin(4) * t10 + t70;
t12 = -t41 * t79 + t42 * t55;
t65 = pkin(4) * t12 + t69;
t64 = rSges(4,1) * t101 - rSges(4,2) * t59;
t61 = cos(qJ(6));
t58 = sin(qJ(6));
t44 = t101 * t98 + t59 * t90;
t26 = t55 * t74 - t57 * t98;
t24 = t101 * t42 + t59 * t76;
t22 = t101 * t40 - t59 * t77;
t20 = t44 * t57 - t55 * t97;
t19 = t44 * t55 + t57 * t97;
t11 = -t41 * t80 - t42 * t57;
t9 = -t39 * t80 - t40 * t57;
t6 = t24 * t57 + t41 * t55;
t5 = t24 * t55 - t41 * t57;
t4 = t22 * t57 + t39 * t55;
t3 = t22 * t55 - t39 * t57;
t1 = [(-m(2) - m(3) - m(4) + t87) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t41 - rSges(3,2) * t42) + g(2) * (-rSges(3,1) * t39 - rSges(3,2) * t40) + (rSges(3,1) * t62 - rSges(3,2) * t60) * t104) - m(4) * (g(1) * (t103 * t42 - t41 * t64 - t37) + g(2) * (t103 * t40 - t39 * t64 - t36) + g(3) * t96 + (rSges(4,3) * t60 + t62 * t64) * t104) - m(5) * (g(1) * (rSges(5,1) * t12 - rSges(5,2) * t11 - rSges(5,3) * t99 + t69) + g(2) * (rSges(5,1) * t10 - rSges(5,2) * t9 - rSges(5,3) * t100 + t70) + g(3) * (rSges(5,1) * t27 - rSges(5,2) * t26 + rSges(5,3) * t86 + t75)) - m(6) * (g(1) * (rSges(6,1) * t12 - rSges(6,2) * t99 + t11 * t91 + t65) + g(2) * (rSges(6,1) * t10 - rSges(6,2) * t100 + t9 * t91 + t66) + g(3) * (t27 * rSges(6,1) + rSges(6,2) * t86 + t26 * t91 + t73)) - m(7) * (g(1) * (t12 * pkin(5) + t11 * qJ(5) + (t11 * t58 + t12 * t61) * rSges(7,1) + (t11 * t61 - t12 * t58) * rSges(7,2) + t65) + g(2) * (t10 * pkin(5) + t9 * qJ(5) + (t10 * t61 + t58 * t9) * rSges(7,1) + (-t10 * t58 + t61 * t9) * rSges(7,2) + t66) + g(3) * (t27 * pkin(5) + t26 * qJ(5) + (t26 * t58 + t27 * t61) * rSges(7,1) + (t26 * t61 - t27 * t58) * rSges(7,2) + t73) + (g(1) * t41 + g(2) * t39 - g(3) * t97) * t59 * t102) -m(4) * (g(1) * (-rSges(4,1) * t23 - rSges(4,2) * t24) + g(2) * (-rSges(4,1) * t21 - rSges(4,2) * t22) + g(3) * (-rSges(4,1) * t43 - rSges(4,2) * t44)) - m(5) * (g(1) * (t23 * t72 + t24 * t92 - t18) + g(2) * (t21 * t72 + t22 * t92 - t17) + g(3) * (t43 * t72 + t44 * t92 - t38)) - m(6) * (g(1) * (t23 * t71 + t24 * t93 + t83) + g(2) * (t21 * t71 + t22 * t93 + t84) + g(3) * (t43 * t71 + t44 * t93 + t82)) + (-g(1) * t83 - g(2) * t84 - g(3) * t82 - (g(1) * t24 + g(2) * t22 + g(3) * t44) * (qJ(4) - t102) - t110 * (-t57 * pkin(5) + (-t55 * t58 - t57 * t61) * rSges(7,1) + (-t55 * t61 + t57 * t58) * rSges(7,2))) * m(7), t87 * t110, t109 * (g(1) * t5 + g(2) * t3 + g(3) * t19) -m(7) * (g(1) * ((t5 * t61 - t58 * t6) * rSges(7,1) + (-t5 * t58 - t6 * t61) * rSges(7,2)) + g(2) * ((t3 * t61 - t4 * t58) * rSges(7,1) + (-t3 * t58 - t4 * t61) * rSges(7,2)) + g(3) * ((t19 * t61 - t20 * t58) * rSges(7,1) + (-t19 * t58 - t20 * t61) * rSges(7,2)))];
taug  = t1(:);
