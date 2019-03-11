% Calculate Gravitation load on the joints for
% S6RRPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
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
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPP1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:45:09
% EndTime: 2019-03-09 09:45:11
% DurationCPUTime: 0.85s
% Computational Cost: add. (476->138), mult. (508->178), div. (0->0), fcn. (485->10), ass. (0->64)
t28 = qJ(2) + pkin(9);
t23 = sin(t28);
t69 = rSges(5,3) + pkin(8);
t80 = t69 * t23;
t25 = cos(t28);
t33 = sin(qJ(1));
t34 = cos(qJ(4));
t62 = t33 * t34;
t31 = sin(qJ(4));
t36 = cos(qJ(1));
t66 = t31 * t36;
t8 = -t25 * t66 + t62;
t71 = rSges(7,1) + pkin(5);
t51 = g(1) * t36 + g(2) * t33;
t55 = rSges(7,3) + qJ(6);
t27 = qJ(4) + pkin(10);
t22 = sin(t27);
t24 = cos(t27);
t79 = t55 * t22 + t71 * t24;
t78 = -m(6) - m(7);
t32 = sin(qJ(2));
t77 = pkin(2) * t32;
t76 = pkin(4) * t31;
t30 = -qJ(3) - pkin(7);
t74 = g(2) * t30;
t72 = g(3) * t23;
t70 = rSges(3,3) + pkin(7);
t68 = t23 * t36;
t20 = pkin(4) * t34 + pkin(3);
t11 = t25 * t20;
t67 = t25 * t36;
t65 = t33 * t22;
t64 = t33 * t24;
t63 = t33 * t31;
t61 = t34 * t36;
t60 = t36 * t22;
t29 = -qJ(5) - pkin(8);
t59 = rSges(7,2) - t29;
t58 = rSges(4,3) - t30;
t57 = rSges(6,3) - t29;
t35 = cos(qJ(2));
t26 = t35 * pkin(2);
t56 = t11 + t26;
t21 = t26 + pkin(1);
t53 = -t21 - t11;
t52 = t33 * t23 * t29 + pkin(4) * t66 - t30 * t36;
t50 = rSges(3,1) * t35 - rSges(3,2) * t32;
t48 = rSges(4,1) * t25 - rSges(4,2) * t23;
t47 = rSges(6,1) * t24 - rSges(6,2) * t22;
t46 = t8 * pkin(4);
t45 = pkin(1) + t50;
t44 = rSges(5,1) * t34 - rSges(5,2) * t31 + pkin(3);
t6 = t25 * t63 + t61;
t14 = t36 * t21;
t43 = pkin(4) * t63 + t20 * t67 - t29 * t68 + t14;
t42 = t25 * pkin(3) + t80;
t41 = t6 * pkin(4);
t9 = t25 * t61 + t63;
t7 = -t25 * t62 + t66;
t5 = t24 * t67 + t65;
t4 = t25 * t60 - t64;
t3 = t25 * t64 - t60;
t2 = t24 * t36 + t25 * t65;
t1 = [-m(2) * (g(1) * (-t33 * rSges(2,1) - rSges(2,2) * t36) + g(2) * (rSges(2,1) * t36 - t33 * rSges(2,2))) - m(3) * ((g(1) * t70 + g(2) * t45) * t36 + (-g(1) * t45 + g(2) * t70) * t33) - m(4) * (g(2) * t14 + (g(1) * t58 + g(2) * t48) * t36 + (g(1) * (-t21 - t48) + g(2) * t58) * t33) - m(5) * (g(1) * (t7 * rSges(5,1) + t6 * rSges(5,2)) + g(2) * (t9 * rSges(5,1) + t8 * rSges(5,2) + t14) + (-g(1) * t30 + g(2) * t42) * t36 + (g(1) * (-t21 - t42) - t74) * t33) - m(6) * (g(1) * (-t3 * rSges(6,1) + t2 * rSges(6,2) + t52) + g(2) * (t5 * rSges(6,1) - t4 * rSges(6,2) + rSges(6,3) * t68 + t43) + (g(1) * (-t23 * rSges(6,3) + t53) - t74) * t33) - m(7) * (g(1) * (-t55 * t2 - t71 * t3 + t52) + g(2) * (rSges(7,2) * t68 + t55 * t4 + t71 * t5 + t43) + (g(1) * (-rSges(7,2) * t23 + t53) - t74) * t33) -m(3) * (g(3) * t50 + t51 * (-rSges(3,1) * t32 - rSges(3,2) * t35)) - m(4) * (g(3) * (t26 + t48) + t51 * (-rSges(4,1) * t23 - rSges(4,2) * t25 - t77)) - m(5) * (g(3) * (t44 * t25 + t26 + t80) + t51 * (-t44 * t23 + t69 * t25 - t77)) - m(6) * (g(3) * (t57 * t23 + t47 * t25 + t56) + t51 * (-t77 + t57 * t25 + (-t20 - t47) * t23)) - m(7) * (g(3) * t56 - t51 * t77 + (g(3) * t79 + t51 * t59) * t25 + (g(3) * t59 + t51 * (-t20 - t79)) * t23) (-m(4) - m(5) + t78) * (g(1) * t33 - g(2) * t36) -m(5) * (g(1) * (rSges(5,1) * t8 - rSges(5,2) * t9) + g(2) * (-rSges(5,1) * t6 + rSges(5,2) * t7)) - m(6) * (g(1) * (-t4 * rSges(6,1) - t5 * rSges(6,2) + t46) + g(2) * (-t2 * rSges(6,1) - t3 * rSges(6,2) - t41)) - m(7) * (g(1) * (-t71 * t4 + t55 * t5 + t46) + g(2) * (-t71 * t2 + t55 * t3 - t41)) + (-m(5) * (-rSges(5,1) * t31 - rSges(5,2) * t34) - m(6) * (-rSges(6,1) * t22 - rSges(6,2) * t24 - t76) - m(7) * (-t71 * t22 + t55 * t24 - t76)) * t72, t78 * (-g(3) * t25 + t51 * t23) -m(7) * (g(1) * t4 + g(2) * t2 + t22 * t72)];
taug  = t1(:);
