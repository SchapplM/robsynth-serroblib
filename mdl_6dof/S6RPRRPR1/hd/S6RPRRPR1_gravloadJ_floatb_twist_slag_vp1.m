% Calculate Gravitation load on the joints for
% S6RPRRPR1
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
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:57:32
% EndTime: 2019-03-09 04:57:34
% DurationCPUTime: 0.57s
% Computational Cost: add. (451->109), mult. (326->135), div. (0->0), fcn. (277->12), ass. (0->62)
t87 = rSges(7,3) + pkin(9);
t36 = qJ(3) + qJ(4);
t29 = pkin(11) + t36;
t23 = sin(t29);
t24 = cos(t29);
t37 = sin(qJ(6));
t73 = rSges(7,2) * t37;
t86 = t23 * t73 + t24 * t87;
t52 = t24 * rSges(6,1) - t23 * rSges(6,2);
t30 = sin(t36);
t31 = cos(t36);
t57 = t31 * rSges(5,1) - t30 * rSges(5,2);
t35 = qJ(1) + pkin(10);
t27 = sin(t35);
t28 = cos(t35);
t84 = g(1) * t28 + g(2) * t27;
t47 = t24 * pkin(5) + t87 * t23;
t40 = cos(qJ(6));
t74 = rSges(7,1) * t40;
t83 = (-pkin(5) - t74) * t23;
t82 = -m(6) - m(7);
t43 = -pkin(8) - pkin(7);
t81 = pkin(4) * t30;
t38 = sin(qJ(3));
t78 = t38 * pkin(3);
t39 = sin(qJ(1));
t77 = t39 * pkin(1);
t76 = rSges(4,3) + pkin(7);
t41 = cos(qJ(3));
t32 = t41 * pkin(3);
t26 = t32 + pkin(2);
t25 = pkin(4) * t31;
t15 = t25 + t26;
t42 = cos(qJ(1));
t33 = t42 * pkin(1);
t75 = t28 * t15 + t33;
t69 = t27 * t37;
t68 = t27 * t40;
t67 = t28 * t37;
t66 = t28 * t40;
t64 = rSges(5,3) - t43;
t34 = -qJ(5) + t43;
t63 = rSges(6,3) - t34;
t62 = g(1) * t77;
t61 = t86 * t27;
t60 = t86 * t28;
t56 = t25 + t52;
t55 = t41 * rSges(4,1) - t38 * rSges(4,2);
t53 = -rSges(5,1) * t30 - rSges(5,2) * t31;
t51 = -rSges(6,1) * t23 - rSges(6,2) * t24;
t50 = g(2) * t33 - t62;
t49 = pkin(2) + t55;
t48 = t26 + t57;
t45 = t25 + t47 + (-t73 + t74) * t24;
t16 = -t78 - t81;
t7 = t28 * t16;
t6 = t27 * t16;
t4 = t24 * t66 + t69;
t3 = -t24 * t67 + t68;
t2 = -t24 * t68 + t67;
t1 = t24 * t69 + t66;
t5 = [-m(2) * (g(1) * (-t39 * rSges(2,1) - t42 * rSges(2,2)) + g(2) * (t42 * rSges(2,1) - t39 * rSges(2,2))) - m(3) * (g(1) * (-t27 * rSges(3,1) - t28 * rSges(3,2) - t77) + g(2) * (t28 * rSges(3,1) - t27 * rSges(3,2) + t33)) - m(4) * ((g(1) * t76 + g(2) * t49) * t28 + (-g(1) * t49 + g(2) * t76) * t27 + t50) - m(5) * ((g(1) * t64 + g(2) * t48) * t28 + (-g(1) * t48 + g(2) * t64) * t27 + t50) - m(6) * (-t62 + g(2) * t75 + (g(1) * t63 + g(2) * t52) * t28 + (g(1) * (-t15 - t52) + g(2) * t63) * t27) - m(7) * (g(1) * (t2 * rSges(7,1) + t1 * rSges(7,2) - t77) + g(2) * (t4 * rSges(7,1) + t3 * rSges(7,2) + t75) + (-g(1) * t34 + g(2) * t47) * t28 + (g(1) * (-t15 - t47) - g(2) * t34) * t27) (-m(3) - m(4) - m(5) + t82) * g(3), -m(4) * (g(3) * t55 + t84 * (-rSges(4,1) * t38 - rSges(4,2) * t41)) - m(5) * (g(3) * (t32 + t57) + t84 * (t53 - t78)) - m(6) * (g(1) * (t28 * t51 + t7) + g(2) * (t27 * t51 + t6) + g(3) * (t32 + t56)) - m(7) * (g(1) * (t7 + t60) + g(2) * (t6 + t61) + g(3) * (t32 + t45) + t84 * t83) -m(7) * (g(1) * t60 + g(2) * t61) + (-m(5) * t57 - m(6) * t56 - m(7) * t45) * g(3) + t84 * (-m(5) * t53 - m(6) * (t51 - t81) - m(7) * (-t81 + t83)) t82 * (g(1) * t27 - g(2) * t28) -m(7) * (g(1) * (t3 * rSges(7,1) - t4 * rSges(7,2)) + g(2) * (-t1 * rSges(7,1) + t2 * rSges(7,2)) + g(3) * (-rSges(7,1) * t37 - rSges(7,2) * t40) * t23)];
taug  = t5(:);
