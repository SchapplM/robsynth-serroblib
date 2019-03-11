% Calculate Gravitation load on the joints for
% S6PRRPRP1
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
% Datum: 2019-03-08 21:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:22:45
% EndTime: 2019-03-08 21:22:48
% DurationCPUTime: 0.97s
% Computational Cost: add. (554->147), mult. (969->223), div. (0->0), fcn. (1123->12), ass. (0->74)
t46 = sin(qJ(2));
t49 = cos(qJ(2));
t69 = sin(pkin(10));
t71 = cos(pkin(6));
t56 = t71 * t69;
t70 = cos(pkin(10));
t27 = -t46 * t56 + t70 * t49;
t45 = sin(qJ(3));
t48 = cos(qJ(3));
t41 = sin(pkin(6));
t64 = t41 * t69;
t102 = -t27 * t45 + t48 * t64;
t78 = t41 * t46;
t101 = -t45 * t78 + t71 * t48;
t57 = t71 * t70;
t25 = t46 * t57 + t69 * t49;
t65 = t41 * t70;
t51 = -t25 * t45 - t48 * t65;
t50 = t51 * pkin(3);
t100 = rSges(6,3) + pkin(9);
t74 = rSges(7,3) + qJ(6) + pkin(9);
t26 = t70 * t46 + t49 * t56;
t24 = t69 * t46 - t49 * t57;
t90 = g(2) * t24;
t99 = g(1) * t26 + t90;
t77 = t41 * t49;
t98 = g(3) * t77 - t99;
t40 = qJ(3) + pkin(11);
t38 = sin(t40);
t39 = cos(t40);
t10 = t25 * t38 + t39 * t65;
t12 = t27 * t38 - t39 * t64;
t20 = t38 * t78 - t71 * t39;
t97 = g(1) * t12 + g(2) * t10 + g(3) * t20;
t47 = cos(qJ(5));
t36 = pkin(5) * t47 + pkin(4);
t96 = t36 * t39 + t74 * t38;
t95 = m(7) * t97;
t94 = pkin(4) * t39;
t88 = g(3) * t41;
t87 = rSges(4,3) + pkin(8);
t44 = sin(qJ(5));
t85 = t25 * t44;
t83 = t27 * t44;
t80 = t39 * t44;
t79 = t39 * t47;
t76 = t44 * t49;
t75 = t47 * t49;
t37 = pkin(3) * t48 + pkin(2);
t43 = -qJ(4) - pkin(8);
t73 = -t24 * t37 - t25 * t43;
t72 = -t26 * t37 - t27 * t43;
t68 = -m(5) - m(6) - m(7);
t62 = t102 * pkin(3);
t59 = rSges(5,1) * t39 - rSges(5,2) * t38;
t11 = t25 * t39 - t38 * t65;
t1 = -t11 * t44 + t24 * t47;
t13 = t27 * t39 + t38 * t64;
t3 = -t13 * t44 + t26 * t47;
t58 = t101 * pkin(3);
t55 = rSges(4,1) * t48 - rSges(4,2) * t45 + pkin(2);
t21 = t71 * t38 + t39 * t78;
t14 = -t21 * t44 - t41 * t75;
t28 = t37 * t77;
t17 = (t39 * t75 + t44 * t46) * t41;
t16 = (-t39 * t76 + t46 * t47) * t41;
t15 = -t21 * t47 + t41 * t76;
t8 = -t26 * t79 + t83;
t7 = t26 * t80 + t27 * t47;
t6 = -t24 * t79 + t85;
t5 = t24 * t80 + t25 * t47;
t4 = -t13 * t47 - t26 * t44;
t2 = -t11 * t47 - t24 * t44;
t9 = [(-m(2) - m(3) - m(4) + t68) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t26 - rSges(3,2) * t27) + g(2) * (-rSges(3,1) * t24 - rSges(3,2) * t25) + (rSges(3,1) * t49 - rSges(3,2) * t46) * t88) - m(4) * (g(1) * (-t55 * t26 + t87 * t27) + g(2) * t87 * t25 - t55 * t90 + (t87 * t46 + t55 * t49) * t88) - m(5) * (g(1) * (rSges(5,3) * t27 - t59 * t26 + t72) + g(2) * (rSges(5,3) * t25 - t59 * t24 + t73) + g(3) * t28 + (t59 * t49 + (rSges(5,3) - t43) * t46) * t88) - m(6) * (g(1) * (rSges(6,1) * t8 + rSges(6,2) * t7 - t26 * t94 + t72) + g(2) * (rSges(6,1) * t6 + rSges(6,2) * t5 - t24 * t94 + t73) + g(3) * (t17 * rSges(6,1) + t16 * rSges(6,2) - t43 * t78 + t77 * t94 + t28) + t98 * t38 * t100) - m(7) * (g(1) * (rSges(7,1) * t8 + rSges(7,2) * t7 + pkin(5) * t83 + t72) + g(2) * (rSges(7,1) * t6 + rSges(7,2) * t5 + pkin(5) * t85 + t73) + g(3) * (t17 * rSges(7,1) + t16 * rSges(7,2) + t28) + ((pkin(5) * t44 - t43) * t46 + t96 * t49) * t88 - t99 * t96) -m(4) * (g(1) * (t102 * rSges(4,1) + (-t27 * t48 - t45 * t64) * rSges(4,2)) + g(2) * (t51 * rSges(4,1) + (-t25 * t48 + t45 * t65) * rSges(4,2)) + g(3) * (t101 * rSges(4,1) + (-t71 * t45 - t48 * t78) * rSges(4,2))) - m(5) * (g(1) * (-rSges(5,1) * t12 - rSges(5,2) * t13 + t62) + g(2) * (-t10 * rSges(5,1) - t11 * rSges(5,2) + t50) + g(3) * (-rSges(5,1) * t20 - rSges(5,2) * t21 + t58)) - m(7) * (g(1) * (t13 * t74 + t62) + g(2) * (t11 * t74 + t50) + g(3) * (t21 * t74 + t58)) - (-rSges(7,1) * t47 + rSges(7,2) * t44 - t36) * t95 + (-g(1) * (t100 * t13 + t62) - g(2) * (t100 * t11 + t50) - g(3) * (t100 * t21 + t58) - t97 * (-rSges(6,1) * t47 + rSges(6,2) * t44 - pkin(4))) * m(6), -t68 * t98, -m(6) * (g(1) * (rSges(6,1) * t3 + rSges(6,2) * t4) + g(2) * (rSges(6,1) * t1 + rSges(6,2) * t2) + g(3) * (rSges(6,1) * t14 + rSges(6,2) * t15)) + (-g(1) * (rSges(7,1) * t3 + rSges(7,2) * t4) - g(2) * (rSges(7,1) * t1 + rSges(7,2) * t2) - g(3) * (t14 * rSges(7,1) + t15 * rSges(7,2)) - (g(1) * t3 + g(2) * t1 + g(3) * t14) * pkin(5)) * m(7), -t95];
taug  = t9(:);
