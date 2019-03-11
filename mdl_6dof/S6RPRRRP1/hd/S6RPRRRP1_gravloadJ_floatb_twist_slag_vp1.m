% Calculate Gravitation load on the joints for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:55:32
% EndTime: 2019-03-09 05:55:34
% DurationCPUTime: 0.75s
% Computational Cost: add. (523->131), mult. (463->171), div. (0->0), fcn. (429->10), ass. (0->58)
t40 = cos(qJ(5));
t82 = rSges(7,1) + pkin(5);
t83 = -t40 * t82 - pkin(4);
t36 = qJ(3) + qJ(4);
t32 = cos(t36);
t25 = t32 * rSges(5,1);
t31 = sin(t36);
t53 = -rSges(5,2) * t31 + t25;
t26 = t31 * pkin(9);
t63 = t32 * pkin(4) + t26;
t35 = qJ(1) + pkin(10);
t29 = sin(t35);
t30 = cos(t35);
t81 = g(1) * t30 + g(2) * t29;
t61 = rSges(7,3) + qJ(6);
t38 = sin(qJ(3));
t80 = pkin(3) * t38;
t43 = -pkin(8) - pkin(7);
t77 = g(2) * t43;
t39 = sin(qJ(1));
t76 = t39 * pkin(1);
t74 = rSges(4,3) + pkin(7);
t70 = t29 * t32;
t12 = pkin(9) * t70;
t73 = rSges(7,2) * t70 + t12;
t41 = cos(qJ(3));
t33 = t41 * pkin(3);
t28 = t33 + pkin(2);
t42 = cos(qJ(1));
t34 = t42 * pkin(1);
t72 = t30 * t28 + t34;
t69 = t30 * t32;
t24 = t31 * rSges(7,2);
t23 = t31 * rSges(6,3);
t37 = sin(qJ(5));
t68 = t31 * t37;
t67 = t32 * t37;
t66 = t32 * t40;
t65 = rSges(5,3) - t43;
t14 = pkin(9) * t69;
t64 = rSges(7,2) * t69 + t14;
t62 = qJ(6) * t37;
t60 = g(1) * t76;
t59 = rSges(6,2) * t68;
t58 = -rSges(6,1) * t40 - pkin(4);
t57 = pkin(4) * t69 + t30 * t26 + t72;
t56 = -t30 * t43 - t76;
t55 = rSges(4,1) * t41 - rSges(4,2) * t38;
t52 = -t28 - t63;
t51 = pkin(2) + t55;
t50 = rSges(7,3) * t67 + t32 * t62 + t82 * t66 + t24 + t63;
t49 = rSges(6,1) * t66 - rSges(6,2) * t67 + t23 + t63;
t45 = g(1) * (rSges(6,3) * t69 + t30 * t59 + t14) + g(2) * (rSges(6,3) * t70 + t29 * t59 + t12);
t4 = t29 * t37 + t30 * t66;
t3 = -t29 * t40 + t30 * t67;
t2 = t29 * t66 - t30 * t37;
t1 = t29 * t67 + t30 * t40;
t5 = [-m(2) * (g(1) * (-rSges(2,1) * t39 - rSges(2,2) * t42) + g(2) * (rSges(2,1) * t42 - rSges(2,2) * t39)) - m(3) * (g(1) * (-rSges(3,1) * t29 - rSges(3,2) * t30 - t76) + g(2) * (rSges(3,1) * t30 - rSges(3,2) * t29 + t34)) - m(4) * (-t60 + g(2) * t34 + (g(1) * t74 + g(2) * t51) * t30 + (-g(1) * t51 + g(2) * t74) * t29) - m(5) * (-t60 + g(2) * t72 + (g(1) * t65 + g(2) * t53) * t30 + (g(1) * (-t28 - t53) + g(2) * t65) * t29) - m(6) * (g(1) * (-t2 * rSges(6,1) + t1 * rSges(6,2) + t56) + g(2) * (t4 * rSges(6,1) - t3 * rSges(6,2) + t30 * t23 + t57) + (g(1) * (t52 - t23) - t77) * t29) - m(7) * (g(1) * (-t61 * t1 - t82 * t2 + t56) + g(2) * (t30 * t24 + t61 * t3 + t4 * t82 + t57) + (g(1) * (t52 - t24) - t77) * t29) (-m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(4) * (g(3) * t55 + t81 * (-rSges(4,1) * t38 - rSges(4,2) * t41)) - m(5) * (g(3) * (t33 + t53) + t81 * (-rSges(5,1) * t31 - rSges(5,2) * t32 - t80)) - m(6) * (g(3) * (t33 + t49) + t45 + t81 * (t58 * t31 - t80)) - m(7) * (g(1) * (-t30 * t80 + t64) + g(2) * (-t29 * t80 + t73) + g(3) * (t33 + t50) + t81 * t31 * (-t61 * t37 + t83)) -m(5) * (g(3) * t25 + (-g(1) * t69 - g(2) * t70) * rSges(5,2)) - m(6) * (g(3) * t49 + t45) - m(7) * (g(1) * t64 + g(2) * t73 + g(3) * t50) + (m(5) * g(3) * rSges(5,2) + t81 * (m(5) * rSges(5,1) - m(6) * t58 - m(7) * (-rSges(7,3) * t37 - t62 + t83))) * t31, -m(6) * (g(1) * (-rSges(6,1) * t3 - rSges(6,2) * t4) + g(2) * (-rSges(6,1) * t1 - rSges(6,2) * t2)) - m(7) * (g(1) * (-t3 * t82 + t61 * t4) + g(2) * (-t1 * t82 + t61 * t2)) + (-m(6) * (-rSges(6,1) * t37 - rSges(6,2) * t40) - m(7) * (-t82 * t37 + t61 * t40)) * g(3) * t31, -m(7) * (g(1) * t3 + g(2) * t1 + g(3) * t68)];
taug  = t5(:);
