% Calculate Gravitation load on the joints for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:29:35
% EndTime: 2019-03-08 19:29:38
% DurationCPUTime: 0.98s
% Computational Cost: add. (573->133), mult. (1351->212), div. (0->0), fcn. (1681->14), ass. (0->71)
t38 = sin(pkin(11));
t45 = sin(qJ(2));
t73 = cos(pkin(11));
t90 = cos(qJ(2));
t23 = -t90 * t38 - t45 * t73;
t39 = sin(pkin(10));
t41 = cos(pkin(10));
t42 = cos(pkin(6));
t50 = -t45 * t38 + t73 * t90;
t48 = t42 * t50;
t12 = t23 * t41 - t39 * t48;
t69 = t42 * t90;
t51 = -t39 * t69 - t41 * t45;
t49 = t51 * pkin(2);
t104 = t12 * pkin(3) + t49;
t103 = -t39 * t45 + t41 * t69;
t72 = sin(pkin(6));
t55 = t73 * t72;
t64 = t45 * t72;
t20 = t38 * t64 - t55 * t90;
t102 = g(1) * t12 - g(3) * t20;
t58 = t90 * t72;
t21 = t38 * t58 + t45 * t55;
t44 = sin(qJ(4));
t46 = cos(qJ(4));
t15 = t21 * t46 + t42 * t44;
t76 = t23 * t42;
t13 = t39 * t76 + t41 * t50;
t66 = t44 * t72;
t5 = t13 * t46 + t39 * t66;
t101 = g(1) * t5 + g(3) * t15;
t8 = -t39 * t50 + t41 * t76;
t14 = t21 * t44 - t42 * t46;
t63 = t46 * t72;
t2 = t41 * t63 - t44 * t8;
t4 = t13 * t44 - t39 * t63;
t100 = g(1) * t4 + g(2) * t2 + g(3) * t14;
t98 = -m(6) - m(7);
t93 = t46 * pkin(4);
t37 = sin(pkin(12));
t92 = t8 * t37;
t91 = rSges(5,3) + pkin(8);
t89 = t13 * t37;
t88 = t21 * t37;
t36 = pkin(12) + qJ(6);
t34 = sin(t36);
t87 = t34 * t46;
t35 = cos(t36);
t86 = t35 * t46;
t85 = t37 * t46;
t40 = cos(pkin(12));
t82 = t40 * t46;
t79 = t42 * t45;
t33 = pkin(5) * t40 + pkin(4);
t78 = t46 * t33;
t77 = rSges(7,3) + pkin(9) + qJ(5);
t30 = pkin(2) * t58;
t75 = -t20 * pkin(3) + t30;
t74 = rSges(6,3) + qJ(5);
t71 = -m(4) - m(5) + t98;
t70 = g(2) * t77;
t67 = g(2) * t74;
t60 = t103 * pkin(2);
t59 = t21 * pkin(8) + t75;
t57 = rSges(5,1) * t46 - rSges(5,2) * t44;
t9 = t39 * t23 + t41 * t48;
t56 = t9 * pkin(3) + t60;
t52 = -t8 * pkin(8) + t56;
t47 = pkin(8) * t13 + t104;
t3 = -t41 * t66 - t46 * t8;
t1 = [(-m(2) - m(3) + t71) * g(3), -m(3) * (g(1) * (t51 * rSges(3,1) + (t39 * t79 - t41 * t90) * rSges(3,2)) + g(2) * (t103 * rSges(3,1) + (-t39 * t90 - t41 * t79) * rSges(3,2)) + g(3) * (rSges(3,1) * t58 - rSges(3,2) * t64)) - m(4) * (g(1) * (t12 * rSges(4,1) - rSges(4,2) * t13 + t49) + g(2) * (t9 * rSges(4,1) + t8 * rSges(4,2) + t60) + g(3) * (-t20 * rSges(4,1) - t21 * rSges(4,2) + t30)) - m(5) * (g(1) * (t12 * t57 + t13 * t91 + t104) + g(2) * (t57 * t9 - t8 * t91 + t56) + g(3) * (-t20 * t57 + t21 * t91 + t75)) - m(6) * (g(1) * (t12 * t93 + (t12 * t82 + t89) * rSges(6,1) + (-t12 * t85 + t13 * t40) * rSges(6,2) + t47) + g(2) * (t9 * t93 + (t82 * t9 - t92) * rSges(6,1) + (-t8 * t40 - t85 * t9) * rSges(6,2) + t52) + g(3) * (-t20 * t93 + (-t20 * t82 + t88) * rSges(6,1) + (t20 * t85 + t21 * t40) * rSges(6,2) + t59) + (t102 * t74 + t67 * t9) * t44) - m(7) * (g(1) * (t12 * t78 + pkin(5) * t89 + (t12 * t86 + t13 * t34) * rSges(7,1) + (-t12 * t87 + t13 * t35) * rSges(7,2) + t47) + g(2) * (t9 * t78 - pkin(5) * t92 + (-t8 * t34 + t9 * t86) * rSges(7,1) + (-t8 * t35 - t9 * t87) * rSges(7,2) + t52) + g(3) * (-t20 * t78 + pkin(5) * t88 + (-t20 * t86 + t21 * t34) * rSges(7,1) + (t20 * t87 + t21 * t35) * rSges(7,2) + t59) + (t102 * t77 + t9 * t70) * t44) t71 * (g(3) * t42 + (g(1) * t39 - g(2) * t41) * t72) -m(5) * (g(1) * (-rSges(5,1) * t4 - rSges(5,2) * t5) + g(2) * (-rSges(5,1) * t2 - rSges(5,2) * t3) + g(3) * (-rSges(5,1) * t14 - rSges(5,2) * t15)) - m(6) * (t3 * t67 + t101 * t74 + t100 * (-rSges(6,1) * t40 + rSges(6,2) * t37 - pkin(4))) - m(7) * (t3 * t70 + t101 * t77 + t100 * (-rSges(7,1) * t35 + rSges(7,2) * t34 - t33)) t98 * t100, -m(7) * (g(1) * ((-t12 * t35 - t34 * t5) * rSges(7,1) + (t12 * t34 - t35 * t5) * rSges(7,2)) + g(2) * ((-t3 * t34 - t35 * t9) * rSges(7,1) + (-t3 * t35 + t34 * t9) * rSges(7,2)) + g(3) * ((-t15 * t34 + t20 * t35) * rSges(7,1) + (-t15 * t35 - t20 * t34) * rSges(7,2)))];
taug  = t1(:);
