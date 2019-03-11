% Calculate Gravitation load on the joints for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-03-08 21:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:28:25
% EndTime: 2019-03-08 21:28:27
% DurationCPUTime: 0.84s
% Computational Cost: add. (615->151), mult. (1085->228), div. (0->0), fcn. (1279->12), ass. (0->80)
t53 = sin(pkin(10));
t58 = sin(qJ(2));
t61 = cos(qJ(2));
t84 = cos(pkin(10));
t85 = cos(pkin(6));
t67 = t85 * t84;
t36 = t53 * t61 + t58 * t67;
t52 = qJ(3) + pkin(11);
t50 = sin(t52);
t51 = cos(t52);
t54 = sin(pkin(6));
t77 = t54 * t84;
t12 = -t36 * t50 - t51 * t77;
t57 = sin(qJ(3));
t60 = cos(qJ(3));
t64 = -t36 * t57 - t60 * t77;
t62 = t64 * pkin(3);
t111 = t12 * pkin(4) + t62;
t93 = t54 * t58;
t110 = -t57 * t93 + t85 * t60;
t78 = t53 * t85;
t38 = -t58 * t78 + t61 * t84;
t92 = t54 * t60;
t109 = -t38 * t57 + t53 * t92;
t108 = rSges(7,2) + pkin(9);
t107 = pkin(4) * t51;
t35 = t53 * t58 - t61 * t67;
t106 = g(2) * t35;
t105 = g(3) * t54;
t104 = rSges(7,1) + pkin(5);
t103 = rSges(4,3) + pkin(8);
t102 = rSges(6,3) + pkin(9);
t101 = t35 * t50;
t37 = t58 * t84 + t61 * t78;
t99 = t37 * t50;
t97 = t50 * t61;
t56 = sin(qJ(5));
t96 = t51 * t56;
t59 = cos(qJ(5));
t95 = t51 * t59;
t94 = t53 * t54;
t91 = t54 * t61;
t55 = -qJ(4) - pkin(8);
t90 = t55 * t58;
t89 = t59 * t61;
t49 = pkin(3) * t60 + pkin(2);
t88 = -t35 * t49 - t36 * t55;
t87 = -t37 * t49 - t38 * t55;
t86 = rSges(7,3) + qJ(6);
t83 = -m(5) - m(6) - m(7);
t80 = t56 * t91;
t41 = t49 * t91;
t79 = t41 + (pkin(9) * t50 + t107) * t91;
t75 = -pkin(9) * t101 - t35 * t107 + t88;
t74 = -pkin(9) * t99 - t37 * t107 + t87;
t73 = t109 * pkin(3);
t71 = rSges(5,1) * t51 - rSges(5,2) * t50;
t70 = rSges(6,1) * t59 - rSges(6,2) * t56;
t69 = t110 * pkin(3);
t14 = -t38 * t50 + t51 * t94;
t68 = t14 * pkin(4) + t73;
t66 = rSges(4,1) * t60 - rSges(4,2) * t57 + pkin(2);
t27 = -t50 * t93 + t51 * t85;
t65 = t27 * pkin(4) + t69;
t28 = t50 * t85 + t51 * t93;
t19 = (t51 * t89 + t56 * t58) * t54;
t18 = t51 * t80 - t59 * t93;
t17 = t28 * t59 - t80;
t16 = t28 * t56 + t54 * t89;
t15 = t38 * t51 + t50 * t94;
t13 = t36 * t51 - t50 * t77;
t8 = -t37 * t95 + t38 * t56;
t7 = -t37 * t96 - t38 * t59;
t6 = -t35 * t95 + t36 * t56;
t5 = -t35 * t96 - t36 * t59;
t4 = t15 * t59 + t37 * t56;
t3 = t15 * t56 - t37 * t59;
t2 = t13 * t59 + t35 * t56;
t1 = t13 * t56 - t35 * t59;
t9 = [(-m(2) - m(3) - m(4) + t83) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t37 - rSges(3,2) * t38) + g(2) * (-rSges(3,1) * t35 - rSges(3,2) * t36) + (rSges(3,1) * t61 - rSges(3,2) * t58) * t105) - m(4) * (g(1) * (t103 * t38 - t37 * t66) + g(2) * t103 * t36 - t66 * t106 + (t103 * t58 + t61 * t66) * t105) - m(5) * (g(1) * (rSges(5,3) * t38 - t37 * t71 + t87) + g(2) * (rSges(5,3) * t36 - t35 * t71 + t88) + g(3) * t41 + (t71 * t61 + (rSges(5,3) - t55) * t58) * t105) - m(6) * (g(1) * (rSges(6,1) * t8 - rSges(6,2) * t7 - rSges(6,3) * t99 + t74) + g(2) * (rSges(6,1) * t6 - rSges(6,2) * t5 - rSges(6,3) * t101 + t75) + g(3) * (t19 * rSges(6,1) - t18 * rSges(6,2) + (rSges(6,3) * t97 - t90) * t54 + t79)) - m(7) * (g(1) * (-rSges(7,2) * t99 + t104 * t8 + t7 * t86 + t74) + g(2) * (-rSges(7,2) * t101 + t104 * t6 + t5 * t86 + t75) + g(3) * ((rSges(7,2) * t97 - t90) * t54 + t104 * t19 + t86 * t18 + t79)) -m(4) * (g(1) * (t109 * rSges(4,1) + (-t38 * t60 - t57 * t94) * rSges(4,2)) + g(2) * (t64 * rSges(4,1) + (-t36 * t60 + t57 * t77) * rSges(4,2)) + g(3) * (t110 * rSges(4,1) + (-t57 * t85 - t58 * t92) * rSges(4,2))) - m(5) * (g(1) * (rSges(5,1) * t14 - rSges(5,2) * t15 + t73) + g(2) * (t12 * rSges(5,1) - t13 * rSges(5,2) + t62) + g(3) * (rSges(5,1) * t27 - rSges(5,2) * t28 + t69)) - m(6) * (g(1) * (t102 * t15 + t14 * t70 + t68) + g(2) * (t102 * t13 + t12 * t70 + t111) + g(3) * (t102 * t28 + t27 * t70 + t65)) + (-g(1) * (t108 * t15 + t68) - g(2) * (t108 * t13 + t111) - g(3) * (t108 * t28 + t65) - (g(1) * t14 + g(2) * t12 + g(3) * t27) * (t104 * t59 + t56 * t86)) * m(7), t83 * (g(1) * t37 - g(3) * t91 + t106) -m(6) * (g(1) * (-rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (-rSges(6,1) * t1 - rSges(6,2) * t2) + g(3) * (-rSges(6,1) * t16 - rSges(6,2) * t17)) - m(7) * (g(1) * (-t104 * t3 + t4 * t86) + g(2) * (-t1 * t104 + t2 * t86) + g(3) * (-t104 * t16 + t17 * t86)) -m(7) * (g(1) * t3 + g(2) * t1 + g(3) * t16)];
taug  = t9(:);
