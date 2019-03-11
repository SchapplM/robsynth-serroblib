% Calculate Gravitation load on the joints for
% S6PRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
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
% Datum: 2019-03-08 22:15
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:10:20
% EndTime: 2019-03-08 22:10:22
% DurationCPUTime: 0.90s
% Computational Cost: add. (602->156), mult. (1524->228), div. (0->0), fcn. (1873->12), ass. (0->80)
t58 = sin(qJ(6));
t62 = cos(qJ(6));
t69 = t62 * rSges(7,1) - t58 * rSges(7,2) + pkin(5);
t57 = cos(pkin(11));
t61 = sin(qJ(2));
t65 = cos(qJ(2));
t89 = sin(pkin(11));
t90 = cos(pkin(6));
t70 = t90 * t89;
t43 = t57 * t65 - t61 * t70;
t60 = sin(qJ(3));
t64 = cos(qJ(3));
t56 = sin(pkin(6));
t82 = t56 * t89;
t24 = t43 * t60 - t64 * t82;
t25 = t43 * t64 + t60 * t82;
t59 = sin(qJ(5));
t63 = cos(qJ(5));
t8 = t24 * t63 - t25 * t59;
t111 = t69 * t8;
t81 = t57 * t90;
t41 = t61 * t81 + t65 * t89;
t97 = t56 * t64;
t22 = t41 * t60 + t57 * t97;
t23 = -t57 * t56 * t60 + t41 * t64;
t4 = t22 * t63 - t23 * t59;
t5 = t22 * t59 + t23 * t63;
t110 = rSges(6,1) * t4 - rSges(6,2) * t5;
t9 = t24 * t59 + t25 * t63;
t109 = rSges(6,1) * t8 - rSges(6,2) * t9;
t98 = t56 * t61;
t44 = t60 * t98 - t64 * t90;
t45 = t60 * t90 + t61 * t97;
t16 = t44 * t63 - t45 * t59;
t108 = (g(2) * t4 + g(3) * t16) * t69;
t17 = t44 * t59 + t45 * t63;
t107 = rSges(6,1) * t16 - rSges(6,2) * t17;
t101 = rSges(7,3) + pkin(10);
t106 = t59 * t64 - t60 * t63;
t105 = g(3) * t56;
t104 = rSges(5,2) + pkin(8);
t103 = rSges(4,3) + pkin(8);
t102 = -rSges(6,3) - pkin(9);
t40 = -t61 * t89 + t65 * t81;
t100 = t40 * t64;
t42 = -t57 * t61 - t65 * t70;
t99 = t42 * t64;
t96 = t56 * t65;
t93 = pkin(2) * t96 + pkin(8) * t98;
t92 = qJ(4) * t60;
t91 = rSges(5,3) + qJ(4);
t88 = -m(5) - m(6) - m(7);
t87 = pkin(8) + t102;
t86 = t60 * t96;
t85 = t64 * t96;
t36 = t40 * pkin(2);
t84 = pkin(3) * t100 + t40 * t92 + t36;
t37 = t42 * pkin(2);
t83 = pkin(3) * t99 + t42 * t92 + t37;
t80 = pkin(4) * t100 + t84;
t79 = pkin(4) * t99 + t83;
t78 = pkin(3) * t85 + qJ(4) * t86 + t93;
t19 = t22 * pkin(3);
t77 = -t22 * pkin(4) + qJ(4) * t23 - t19;
t21 = t24 * pkin(3);
t76 = -t24 * pkin(4) + qJ(4) * t25 - t21;
t39 = t44 * pkin(3);
t75 = -t44 * pkin(4) + qJ(4) * t45 - t39;
t74 = pkin(4) * t85 + t78;
t73 = rSges(4,1) * t64 - rSges(4,2) * t60;
t72 = rSges(5,1) * t64 + rSges(5,3) * t60;
t71 = t59 * t60 + t63 * t64;
t68 = -t58 * rSges(7,1) - t62 * rSges(7,2) + pkin(8) - pkin(9);
t27 = t71 * t96;
t26 = t59 * t85 - t63 * t86;
t13 = t71 * t42;
t12 = t106 * t42;
t11 = t71 * t40;
t10 = t106 * t40;
t1 = [(-m(2) - m(3) - m(4) + t88) * g(3), -m(3) * (g(1) * (rSges(3,1) * t42 - t43 * rSges(3,2)) + g(2) * (rSges(3,1) * t40 - rSges(3,2) * t41) + (rSges(3,1) * t65 - rSges(3,2) * t61) * t105) - m(4) * (g(1) * (t103 * t43 + t42 * t73 + t37) + g(2) * (t103 * t41 + t40 * t73 + t36) + g(3) * t93 + (rSges(4,3) * t61 + t65 * t73) * t105) - m(5) * (g(1) * (t104 * t43 + t42 * t72 + t83) + g(2) * (t104 * t41 + t40 * t72 + t84) + g(3) * t78 + (rSges(5,2) * t61 + t65 * t72) * t105) - m(6) * (g(3) * (rSges(6,1) * t27 - rSges(6,2) * t26 + t102 * t98 + t74) + (rSges(6,1) * t11 - rSges(6,2) * t10 + t41 * t87 + t80) * g(2) + (rSges(6,1) * t13 - rSges(6,2) * t12 + t43 * t87 + t79) * g(1)) - m(7) * (g(1) * (t101 * t12 + t13 * t69 + t43 * t68 + t79) + g(2) * (t10 * t101 + t11 * t69 + t41 * t68 + t80) + g(3) * (t27 * pkin(5) - pkin(9) * t98 + (t27 * t62 - t58 * t98) * rSges(7,1) + (-t27 * t58 - t62 * t98) * rSges(7,2) + t101 * t26 + t74)) -m(4) * (g(1) * (-rSges(4,1) * t24 - rSges(4,2) * t25) + g(2) * (-rSges(4,1) * t22 - rSges(4,2) * t23) + g(3) * (-rSges(4,1) * t44 - rSges(4,2) * t45)) - m(5) * (g(1) * (-rSges(5,1) * t24 + t25 * t91 - t21) + g(2) * (-rSges(5,1) * t22 + t23 * t91 - t19) + g(3) * (-rSges(5,1) * t44 + t45 * t91 - t39)) - m(6) * (g(1) * (-t109 + t76) + g(2) * (-t110 + t77) + g(3) * (-t107 + t75)) - m(7) * (g(2) * (-t5 * t101 + t77) + g(3) * (-t17 * t101 + t75) + (-t9 * t101 - t111 + t76) * g(1) - t108) t88 * (g(1) * t24 + g(2) * t22 + g(3) * t44) -m(6) * (g(1) * t109 + g(2) * t110 + g(3) * t107) - m(7) * (g(1) * t111 + (g(1) * t9 + g(2) * t5 + g(3) * t17) * t101 + t108) -m(7) * (g(1) * ((t42 * t62 - t58 * t9) * rSges(7,1) + (-t42 * t58 - t62 * t9) * rSges(7,2)) + g(2) * ((t40 * t62 - t5 * t58) * rSges(7,1) + (-t40 * t58 - t5 * t62) * rSges(7,2)) + g(3) * ((-t17 * t58 + t62 * t96) * rSges(7,1) + (-t17 * t62 - t58 * t96) * rSges(7,2)))];
taug  = t1(:);
