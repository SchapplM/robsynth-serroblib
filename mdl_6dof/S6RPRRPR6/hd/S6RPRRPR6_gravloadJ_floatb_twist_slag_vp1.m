% Calculate Gravitation load on the joints for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:15:36
% EndTime: 2019-03-09 05:15:37
% DurationCPUTime: 0.60s
% Computational Cost: add. (472->133), mult. (450->185), div. (0->0), fcn. (421->12), ass. (0->68)
t71 = rSges(5,3) + pkin(8);
t39 = -qJ(5) - pkin(8);
t59 = rSges(6,3) - t39;
t58 = rSges(7,3) + pkin(9) - t39;
t42 = sin(qJ(1));
t44 = cos(qJ(1));
t79 = g(1) * t44 + g(2) * t42;
t78 = -m(6) - m(7);
t36 = qJ(4) + pkin(11);
t32 = qJ(6) + t36;
t25 = cos(t32);
t35 = pkin(10) + qJ(3);
t30 = cos(t35);
t24 = sin(t32);
t67 = t42 * t24;
t6 = t25 * t44 + t30 * t67;
t66 = t42 * t25;
t7 = t24 * t44 - t30 * t66;
t77 = -t6 * rSges(7,1) + t7 * rSges(7,2);
t70 = t30 * t44;
t8 = -t24 * t70 + t66;
t9 = t25 * t70 + t67;
t76 = t8 * rSges(7,1) - t9 * rSges(7,2);
t41 = sin(qJ(4));
t75 = pkin(4) * t41;
t28 = sin(t35);
t72 = g(3) * t28;
t31 = cos(t36);
t69 = t31 * t44;
t68 = t41 * t44;
t29 = sin(t36);
t65 = t42 * t29;
t64 = t42 * t31;
t63 = t42 * t41;
t43 = cos(qJ(4));
t62 = t42 * t43;
t61 = t43 * t44;
t40 = -pkin(7) - qJ(2);
t60 = rSges(4,3) - t40;
t19 = pkin(5) * t29 + t75;
t57 = t19 - t40;
t33 = t43 * pkin(4);
t20 = pkin(5) * t31 + t33;
t56 = rSges(3,3) + qJ(2);
t55 = -t40 + t75;
t54 = rSges(4,1) * t30 - rSges(4,2) * t28;
t52 = -rSges(7,1) * t24 - rSges(7,2) * t25;
t38 = cos(pkin(10));
t51 = rSges(3,1) * t38 - rSges(3,2) * sin(pkin(10)) + pkin(1);
t50 = rSges(5,1) * t43 - rSges(5,2) * t41 + pkin(3);
t16 = -t30 * t68 + t62;
t14 = t30 * t63 + t61;
t27 = t33 + pkin(3);
t49 = rSges(6,1) * t31 - rSges(6,2) * t29 + t27;
t18 = pkin(3) + t20;
t48 = rSges(7,1) * t25 - rSges(7,2) * t24 + t18;
t47 = t30 * pkin(3) + t71 * t28;
t46 = t30 * t27 + t59 * t28;
t45 = t30 * t18 + t58 * t28;
t26 = pkin(2) * t38 + pkin(1);
t21 = t44 * t26;
t17 = t30 * t61 + t63;
t15 = -t30 * t62 + t68;
t13 = t30 * t69 + t65;
t12 = -t29 * t70 + t64;
t11 = t29 * t44 - t30 * t64;
t10 = t30 * t65 + t69;
t1 = [-m(2) * (g(1) * (-t42 * rSges(2,1) - rSges(2,2) * t44) + g(2) * (rSges(2,1) * t44 - t42 * rSges(2,2))) - m(3) * ((g(1) * t56 + g(2) * t51) * t44 + (-g(1) * t51 + g(2) * t56) * t42) - m(4) * (g(2) * t21 + (g(1) * t60 + g(2) * t54) * t44 + (g(1) * (-t26 - t54) + g(2) * t60) * t42) - m(5) * (g(1) * (t15 * rSges(5,1) + t14 * rSges(5,2)) + g(2) * (t17 * rSges(5,1) + t16 * rSges(5,2) + t21) + (-g(1) * t40 + g(2) * t47) * t44 + (g(1) * (-t26 - t47) - g(2) * t40) * t42) - m(6) * (g(1) * (t11 * rSges(6,1) + t10 * rSges(6,2)) + g(2) * (t13 * rSges(6,1) + t12 * rSges(6,2) + t21) + (g(1) * t55 + g(2) * t46) * t44 + (g(1) * (-t26 - t46) + g(2) * t55) * t42) - m(7) * (g(1) * (t7 * rSges(7,1) + t6 * rSges(7,2)) + g(2) * (t9 * rSges(7,1) + t8 * rSges(7,2) + t21) + (g(1) * t57 + g(2) * t45) * t44 + (g(1) * (-t26 - t45) + g(2) * t57) * t42) (-m(3) - m(4) - m(5) + t78) * (g(1) * t42 - g(2) * t44) -m(4) * (g(3) * t54 + t79 * (-rSges(4,1) * t28 - rSges(4,2) * t30)) - m(5) * ((g(3) * t50 + t79 * t71) * t30 + (g(3) * t71 - t79 * t50) * t28) - m(6) * ((g(3) * t49 + t79 * t59) * t30 + (g(3) * t59 - t79 * t49) * t28) - m(7) * ((g(3) * t48 + t79 * t58) * t30 + (g(3) * t58 - t79 * t48) * t28) -m(5) * (g(1) * (rSges(5,1) * t16 - rSges(5,2) * t17) + g(2) * (-rSges(5,1) * t14 + rSges(5,2) * t15)) - m(6) * (g(1) * (t12 * rSges(6,1) - t13 * rSges(6,2) + t16 * pkin(4)) + g(2) * (-t10 * rSges(6,1) + t11 * rSges(6,2) - t14 * pkin(4))) - m(7) * (g(1) * (-t19 * t70 + t42 * t20 + t76) + g(2) * (-t42 * t30 * t19 - t20 * t44 + t77)) + (-m(5) * (-rSges(5,1) * t41 - rSges(5,2) * t43) - m(6) * (-rSges(6,1) * t29 - rSges(6,2) * t31 - t75) - m(7) * (-t19 + t52)) * t72, t78 * (-g(3) * t30 + t79 * t28) -m(7) * (g(1) * t76 + g(2) * t77 + t52 * t72)];
taug  = t1(:);
