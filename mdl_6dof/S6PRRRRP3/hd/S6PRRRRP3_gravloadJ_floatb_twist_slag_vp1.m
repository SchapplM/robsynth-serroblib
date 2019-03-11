% Calculate Gravitation load on the joints for
% S6PRRRRP3
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
% Datum: 2019-03-09 00:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:06:23
% EndTime: 2019-03-09 00:06:26
% DurationCPUTime: 0.96s
% Computational Cost: add. (644->157), mult. (1279->234), div. (0->0), fcn. (1517->12), ass. (0->80)
t57 = sin(qJ(2));
t60 = cos(qJ(2));
t83 = cos(pkin(11));
t84 = cos(pkin(6));
t72 = t84 * t83;
t82 = sin(pkin(11));
t31 = t82 * t57 - t60 * t72;
t71 = t84 * t82;
t33 = t83 * t57 + t60 * t71;
t113 = g(1) * t33 + g(2) * t31;
t53 = qJ(4) + qJ(5);
t50 = cos(t53);
t58 = cos(qJ(4));
t51 = t58 * pkin(4);
t39 = pkin(5) * t50 + t51;
t37 = pkin(3) + t39;
t56 = sin(qJ(3));
t59 = cos(qJ(3));
t61 = -pkin(10) - pkin(9);
t87 = rSges(7,3) + qJ(6) - t61;
t112 = t37 * t59 + t87 * t56;
t48 = t51 + pkin(3);
t88 = rSges(6,3) - t61;
t111 = t48 * t59 + t88 * t56;
t32 = t57 * t72 + t82 * t60;
t54 = sin(pkin(6));
t78 = t54 * t83;
t24 = t32 * t59 - t56 * t78;
t34 = -t57 * t71 + t83 * t60;
t77 = t54 * t82;
t26 = t34 * t59 + t56 * t77;
t91 = t54 * t57;
t36 = t84 * t56 + t59 * t91;
t110 = g(1) * t26 + g(2) * t24 + g(3) * t36;
t23 = t32 * t56 + t59 * t78;
t25 = t34 * t56 - t59 * t77;
t35 = t56 * t91 - t84 * t59;
t109 = g(1) * t25 + g(2) * t23 + g(3) * t35;
t55 = sin(qJ(4));
t70 = rSges(5,1) * t58 - rSges(5,2) * t55 + pkin(3);
t96 = rSges(5,3) + pkin(9);
t108 = t96 * t56 + t70 * t59;
t49 = sin(t53);
t10 = -t24 * t50 - t31 * t49;
t9 = -t24 * t49 + t31 * t50;
t107 = t9 * rSges(7,1) + t10 * rSges(7,2);
t11 = -t26 * t49 + t33 * t50;
t12 = -t26 * t50 - t33 * t49;
t106 = t11 * rSges(7,1) + t12 * rSges(7,2);
t105 = pkin(4) * t55;
t98 = g(3) * t54;
t97 = rSges(4,3) + pkin(8);
t93 = t49 * t59;
t92 = t50 * t59;
t90 = t54 * t60;
t89 = t59 * t60;
t21 = -t36 * t49 - t50 * t90;
t22 = -t36 * t50 + t49 * t90;
t86 = t21 * rSges(7,1) + t22 * rSges(7,2);
t85 = pkin(2) * t90 + pkin(8) * t91;
t81 = g(3) * t85;
t29 = t31 * pkin(2);
t80 = t32 * pkin(8) - t29;
t30 = t33 * pkin(2);
t79 = t34 * pkin(8) - t30;
t76 = rSges(4,1) * t59 - rSges(4,2) * t56;
t75 = rSges(5,1) * t55 + rSges(5,2) * t58;
t74 = -t24 * t55 + t31 * t58;
t73 = -t26 * t55 + t33 * t58;
t69 = pkin(8) + t75;
t68 = -t36 * t55 - t58 * t90;
t62 = m(6) * (g(1) * (t11 * rSges(6,1) + t12 * rSges(6,2)) + g(2) * (t9 * rSges(6,1) + t10 * rSges(6,2)) + g(3) * (t21 * rSges(6,1) + t22 * rSges(6,2)));
t38 = pkin(5) * t49 + t105;
t28 = (t49 * t57 + t50 * t89) * t54;
t27 = (-t49 * t89 + t50 * t57) * t54;
t16 = -t33 * t92 + t34 * t49;
t15 = t33 * t93 + t34 * t50;
t14 = -t31 * t92 + t32 * t49;
t13 = t31 * t93 + t32 * t50;
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t33 - rSges(3,2) * t34) + g(2) * (-rSges(3,1) * t31 - rSges(3,2) * t32) + (rSges(3,1) * t60 - rSges(3,2) * t57) * t98) - m(4) * (g(1) * (-t76 * t33 + t97 * t34 - t30) + g(2) * (-t76 * t31 + t97 * t32 - t29) + t81 + (rSges(4,3) * t57 + t76 * t60) * t98) - m(5) * (t81 + (t108 * t60 + t75 * t57) * t98 - t113 * t108 + (t69 * t32 - t29) * g(2) + (t69 * t34 - t30) * g(1)) - m(6) * (g(1) * (t16 * rSges(6,1) + t15 * rSges(6,2) + t34 * t105 + t79) + g(2) * (t14 * rSges(6,1) + t13 * rSges(6,2) + t32 * t105 + t80) + g(3) * (t28 * rSges(6,1) + t27 * rSges(6,2) + t85) + (t57 * t105 + t111 * t60) * t98 - t113 * t111) - m(7) * (g(1) * (rSges(7,1) * t16 + rSges(7,2) * t15 + t34 * t38 + t79) + g(2) * (rSges(7,1) * t14 + rSges(7,2) * t13 + t32 * t38 + t80) + g(3) * (rSges(7,1) * t28 + rSges(7,2) * t27 + t85) + (t112 * t60 + t38 * t57) * t98 - t113 * t112) -m(4) * (g(1) * (-rSges(4,1) * t25 - rSges(4,2) * t26) + g(2) * (-rSges(4,1) * t23 - rSges(4,2) * t24) + g(3) * (-rSges(4,1) * t35 - rSges(4,2) * t36)) - m(5) * (-t109 * t70 + t110 * t96) - m(6) * (t110 * t88 + t109 * (-rSges(6,1) * t50 + rSges(6,2) * t49 - t48)) - m(7) * (t110 * t87 + t109 * (-rSges(7,1) * t50 + rSges(7,2) * t49 - t37)) -m(5) * (g(1) * (t73 * rSges(5,1) + (-t26 * t58 - t33 * t55) * rSges(5,2)) + g(2) * (t74 * rSges(5,1) + (-t24 * t58 - t31 * t55) * rSges(5,2)) + g(3) * (t68 * rSges(5,1) + (-t36 * t58 + t55 * t90) * rSges(5,2))) - t62 - m(7) * (g(1) * (-t26 * t38 + t33 * t39 + t106) + g(2) * (-t24 * t38 + t31 * t39 + t107) + g(3) * (-t36 * t38 - t39 * t90 + t86)) - m(6) * (g(1) * t73 + g(2) * t74 + g(3) * t68) * pkin(4), -t62 + (-g(1) * t106 - g(2) * t107 - g(3) * t86 - (g(1) * t11 + g(2) * t9 + g(3) * t21) * pkin(5)) * m(7), -m(7) * t109];
taug  = t1(:);
