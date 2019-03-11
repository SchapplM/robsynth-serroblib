% Calculate Gravitation load on the joints for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:14:04
% EndTime: 2019-03-08 20:14:05
% DurationCPUTime: 0.66s
% Computational Cost: add. (367->121), mult. (882->188), div. (0->0), fcn. (1029->10), ass. (0->64)
t35 = sin(pkin(10));
t40 = sin(qJ(2));
t43 = cos(qJ(2));
t59 = cos(pkin(10));
t60 = cos(pkin(6));
t51 = t60 * t59;
t23 = t35 * t43 + t40 * t51;
t54 = t35 * t60;
t25 = -t40 * t54 + t59 * t43;
t84 = g(1) * t25 + g(2) * t23;
t24 = t59 * t40 + t43 * t54;
t39 = sin(qJ(4));
t42 = cos(qJ(4));
t36 = sin(pkin(6));
t71 = t35 * t36;
t10 = -t24 * t42 + t39 * t71;
t22 = t35 * t40 - t43 * t51;
t53 = t36 * t59;
t12 = t22 * t42 + t39 * t53;
t69 = t36 * t43;
t26 = t60 * t39 + t42 * t69;
t83 = g(1) * t10 - g(2) * t12 + g(3) * t26;
t11 = t24 * t39 + t42 * t71;
t13 = -t22 * t39 + t42 * t53;
t27 = -t39 * t69 + t60 * t42;
t82 = g(1) * t11 - g(2) * t13 + g(3) * t27;
t75 = g(3) * t36;
t74 = rSges(6,3) + pkin(9);
t38 = sin(qJ(5));
t73 = t22 * t38;
t72 = t24 * t38;
t70 = t36 * t40;
t68 = t38 * t39;
t67 = t38 * t40;
t66 = t38 * t43;
t41 = cos(qJ(5));
t65 = t39 * t41;
t64 = t40 * t41;
t63 = rSges(7,3) + qJ(6) + pkin(9);
t62 = pkin(2) * t69 + qJ(3) * t70;
t61 = rSges(4,3) + qJ(3);
t58 = pkin(8) * t69 + t62;
t57 = -m(4) - m(5) - m(6) - m(7);
t20 = t22 * pkin(2);
t56 = -pkin(8) * t22 - t20;
t21 = t24 * pkin(2);
t55 = -pkin(8) * t24 - t21;
t52 = rSges(5,1) * t39 + rSges(5,2) * t42;
t1 = -t11 * t38 + t25 * t41;
t3 = t13 * t38 + t23 * t41;
t14 = -t27 * t38 + t36 * t64;
t47 = pkin(4) * t39 - t74 * t42;
t34 = pkin(5) * t41 + pkin(4);
t46 = t34 * t39 - t63 * t42;
t17 = (t39 * t64 + t66) * t36;
t16 = (-t39 * t67 + t41 * t43) * t36;
t15 = -t27 * t41 - t36 * t67;
t9 = t25 * t65 - t72;
t8 = -t24 * t41 - t25 * t68;
t7 = t23 * t65 - t73;
t6 = -t22 * t41 - t23 * t68;
t4 = t13 * t41 - t23 * t38;
t2 = -t11 * t41 - t25 * t38;
t5 = [(-m(2) - m(3) + t57) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t24 - rSges(3,2) * t25) + g(2) * (-rSges(3,1) * t22 - rSges(3,2) * t23) + (rSges(3,1) * t43 - rSges(3,2) * t40) * t75) - m(4) * (g(1) * (rSges(4,2) * t24 + t61 * t25 - t21) + g(2) * (rSges(4,2) * t22 + t61 * t23 - t20) + g(3) * ((-rSges(4,2) * t43 + rSges(4,3) * t40) * t36 + t62)) - m(5) * (g(1) * (-rSges(5,3) * t24 + t55) + g(2) * (-rSges(5,3) * t22 + t56) + g(3) * t58 + (rSges(5,3) * t43 + t52 * t40) * t75 + t84 * (qJ(3) + t52)) - m(6) * (g(1) * (rSges(6,1) * t9 + rSges(6,2) * t8 + t55) + g(2) * (rSges(6,1) * t7 + rSges(6,2) * t6 + t56) + t84 * (qJ(3) + t47) + (rSges(6,1) * t17 + rSges(6,2) * t16 + t47 * t70 + t58) * g(3)) - m(7) * (g(1) * (rSges(7,1) * t9 + rSges(7,2) * t8 - pkin(5) * t72 + t55) + g(2) * (rSges(7,1) * t7 + rSges(7,2) * t6 - pkin(5) * t73 + t56) + g(3) * (t17 * rSges(7,1) + t16 * rSges(7,2) + t58) + (pkin(5) * t66 + t46 * t40) * t75 + t84 * (qJ(3) + t46)) t57 * (g(1) * t24 + g(2) * t22 - g(3) * t69) -m(5) * (g(1) * (-rSges(5,1) * t10 - rSges(5,2) * t11) + g(2) * (rSges(5,1) * t12 + rSges(5,2) * t13) + g(3) * (-rSges(5,1) * t26 - rSges(5,2) * t27)) - m(6) * (t82 * t74 - t83 * (rSges(6,1) * t41 - rSges(6,2) * t38 + pkin(4))) - m(7) * (t82 * t63 - t83 * (rSges(7,1) * t41 - rSges(7,2) * t38 + t34)) -m(6) * (g(1) * (rSges(6,1) * t1 + rSges(6,2) * t2) + g(2) * (rSges(6,1) * t3 + rSges(6,2) * t4) + g(3) * (rSges(6,1) * t14 + rSges(6,2) * t15)) + (-g(1) * (rSges(7,1) * t1 + rSges(7,2) * t2) - g(2) * (rSges(7,1) * t3 + rSges(7,2) * t4) - g(3) * (rSges(7,1) * t14 + rSges(7,2) * t15) - (g(1) * t1 + g(2) * t3 + g(3) * t14) * pkin(5)) * m(7), -m(7) * t83];
taug  = t5(:);
