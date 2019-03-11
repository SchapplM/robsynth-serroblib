% Calculate Gravitation load on the joints for
% S6PRRPRP3
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
% Datum: 2019-03-08 21:39
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:34:13
% EndTime: 2019-03-08 21:34:16
% DurationCPUTime: 1.06s
% Computational Cost: add. (588->140), mult. (1190->206), div. (0->0), fcn. (1420->12), ass. (0->76)
t55 = sin(pkin(11));
t106 = pkin(4) * t55;
t114 = pkin(8) + t106;
t60 = sin(qJ(2));
t62 = cos(qJ(2));
t82 = cos(pkin(10));
t83 = cos(pkin(6));
t68 = t83 * t82;
t81 = sin(pkin(10));
t35 = t60 * t68 + t81 * t62;
t59 = sin(qJ(3));
t61 = cos(qJ(3));
t56 = sin(pkin(6));
t76 = t56 * t82;
t16 = t35 * t59 + t61 * t76;
t67 = t83 * t81;
t37 = -t60 * t67 + t82 * t62;
t75 = t56 * t81;
t18 = t37 * t59 - t61 * t75;
t94 = t56 * t60;
t38 = t59 * t94 - t83 * t61;
t113 = -g(1) * t18 - g(2) * t16 - g(3) * t38;
t112 = rSges(7,1) + pkin(5);
t84 = rSges(7,3) + qJ(6);
t34 = t81 * t60 - t62 * t68;
t36 = t82 * t60 + t62 * t67;
t110 = g(1) * t36 + g(2) * t34;
t93 = t56 * t62;
t108 = (g(3) * t93 - t110) * t59;
t57 = cos(pkin(11));
t66 = rSges(5,1) * t57 - rSges(5,2) * t55 + pkin(3);
t85 = rSges(5,3) + qJ(4);
t107 = t85 * t59 + t66 * t61;
t100 = g(3) * t56;
t98 = rSges(4,3) + pkin(8);
t51 = pkin(4) * t57 + pkin(3);
t97 = t51 * t61;
t54 = pkin(11) + qJ(5);
t52 = sin(t54);
t96 = t52 * t61;
t53 = cos(t54);
t95 = t53 * t61;
t92 = t61 * t62;
t17 = t35 * t61 - t59 * t76;
t58 = -pkin(9) - qJ(4);
t89 = -t16 * t51 - t17 * t58;
t19 = t37 * t61 + t59 * t75;
t88 = -t18 * t51 - t19 * t58;
t39 = t83 * t59 + t61 * t94;
t87 = -t38 * t51 - t39 * t58;
t86 = pkin(2) * t93 + pkin(8) * t94;
t80 = -m(5) - m(6) - m(7);
t78 = t52 * t93;
t77 = g(3) * t86;
t74 = t56 * t51 * t92 + t94 * t106 + t86;
t73 = rSges(4,1) * t61 - rSges(4,2) * t59;
t72 = t55 * rSges(5,1) + t57 * rSges(5,2);
t71 = -rSges(6,1) * t53 + rSges(6,2) * t52;
t32 = t34 * pkin(2);
t70 = t114 * t35 - t34 * t97 - t32;
t33 = t36 * pkin(2);
t69 = t114 * t37 - t36 * t97 - t33;
t65 = pkin(8) + t72;
t21 = (t52 * t60 + t53 * t92) * t56;
t20 = -t53 * t94 + t61 * t78;
t13 = t39 * t53 - t78;
t12 = t39 * t52 + t53 * t93;
t9 = -t36 * t95 + t37 * t52;
t8 = -t36 * t96 - t37 * t53;
t7 = -t34 * t95 + t35 * t52;
t6 = -t34 * t96 - t35 * t53;
t5 = t19 * t53 + t36 * t52;
t4 = t19 * t52 - t36 * t53;
t3 = t17 * t53 + t34 * t52;
t2 = t17 * t52 - t34 * t53;
t1 = [(-m(2) - m(3) - m(4) + t80) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t36 - rSges(3,2) * t37) + g(2) * (-rSges(3,1) * t34 - rSges(3,2) * t35) + (rSges(3,1) * t62 - rSges(3,2) * t60) * t100) - m(4) * (g(1) * (-t73 * t36 + t98 * t37 - t33) + g(2) * (-t73 * t34 + t98 * t35 - t32) + t77 + (rSges(4,3) * t60 + t73 * t62) * t100) - m(5) * (t77 + (t107 * t62 + t72 * t60) * t100 - t110 * t107 + (t65 * t35 - t32) * g(2) + (t65 * t37 - t33) * g(1)) - m(6) * (g(1) * (rSges(6,1) * t9 - rSges(6,2) * t8 + t69) + g(2) * (rSges(6,1) * t7 - rSges(6,2) * t6 + t70) + g(3) * (t21 * rSges(6,1) - t20 * rSges(6,2) + t74) + (rSges(6,3) - t58) * t108) - m(7) * (g(1) * (t112 * t9 + t84 * t8 + t69) + g(2) * (t112 * t7 + t84 * t6 + t70) + g(3) * (t112 * t21 + t84 * t20 + t74) + (rSges(7,2) - t58) * t108) -m(4) * (g(1) * (-rSges(4,1) * t18 - rSges(4,2) * t19) + g(2) * (-rSges(4,1) * t16 - rSges(4,2) * t17) + g(3) * (-rSges(4,1) * t38 - rSges(4,2) * t39)) - m(5) * ((g(1) * t19 + g(2) * t17 + g(3) * t39) * t85 + t113 * t66) - m(6) * (g(1) * (rSges(6,3) * t19 + t71 * t18 + t88) + g(2) * (rSges(6,3) * t17 + t71 * t16 + t89) + g(3) * (rSges(6,3) * t39 + t71 * t38 + t87)) + (-g(1) * (rSges(7,2) * t19 + t88) - g(2) * (rSges(7,2) * t17 + t89) - g(3) * (rSges(7,2) * t39 + t87) + t113 * (-t112 * t53 - t84 * t52)) * m(7), -t80 * t113, -m(6) * (g(1) * (-rSges(6,1) * t4 - rSges(6,2) * t5) + g(2) * (-rSges(6,1) * t2 - rSges(6,2) * t3) + g(3) * (-rSges(6,1) * t12 - rSges(6,2) * t13)) - m(7) * (g(1) * (-t112 * t4 + t84 * t5) + g(2) * (-t112 * t2 + t84 * t3) + g(3) * (-t112 * t12 + t84 * t13)) -m(7) * (g(1) * t4 + g(2) * t2 + g(3) * t12)];
taug  = t1(:);
