% Calculate Gravitation load on the joints for
% S6RPRRPR2
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
% Datum: 2019-03-09 05:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:00:15
% EndTime: 2019-03-09 05:00:17
% DurationCPUTime: 0.67s
% Computational Cost: add. (484->132), mult. (432->180), div. (0->0), fcn. (403->12), ass. (0->63)
t43 = cos(qJ(3));
t40 = sin(qJ(3));
t66 = rSges(5,3) + pkin(8);
t57 = t66 * t40;
t79 = t43 * pkin(3) + t57;
t38 = -qJ(5) - pkin(8);
t59 = rSges(7,3) + pkin(9) - t38;
t78 = t59 * t40;
t60 = rSges(6,3) - t38;
t77 = t60 * t40;
t37 = qJ(1) + pkin(10);
t29 = sin(t37);
t31 = cos(t37);
t76 = g(1) * t31 + g(2) * t29;
t75 = -m(6) - m(7);
t36 = qJ(4) + pkin(11);
t32 = qJ(6) + t36;
t25 = sin(t32);
t26 = cos(t32);
t64 = t29 * t43;
t5 = t25 * t64 + t26 * t31;
t6 = t25 * t31 - t26 * t64;
t74 = -t5 * rSges(7,1) + t6 * rSges(7,2);
t63 = t31 * t43;
t7 = -t25 * t63 + t26 * t29;
t8 = t25 * t29 + t26 * t63;
t73 = t7 * rSges(7,1) - t8 * rSges(7,2);
t39 = sin(qJ(4));
t72 = pkin(4) * t39;
t69 = g(3) * t40;
t41 = sin(qJ(1));
t68 = t41 * pkin(1);
t65 = rSges(4,2) * t40;
t62 = t39 * t43;
t42 = cos(qJ(4));
t61 = t42 * t43;
t30 = cos(t36);
t33 = t42 * pkin(4);
t20 = pkin(5) * t30 + t33;
t44 = cos(qJ(1));
t34 = t44 * pkin(1);
t58 = t31 * pkin(2) + t29 * pkin(7) + t34;
t56 = t31 * pkin(7) - t68;
t55 = rSges(4,1) * t43 - t65;
t53 = -rSges(7,1) * t25 - rSges(7,2) * t26;
t52 = rSges(5,1) * t42 - rSges(5,2) * t39 + pkin(3);
t16 = t29 * t42 - t31 * t62;
t14 = t29 * t62 + t31 * t42;
t27 = t33 + pkin(3);
t28 = sin(t36);
t51 = rSges(6,1) * t30 - rSges(6,2) * t28 + t27;
t18 = pkin(3) + t20;
t50 = rSges(7,1) * t26 - rSges(7,2) * t25 + t18;
t49 = t43 * t27 + t77;
t48 = t43 * t18 + t78;
t19 = pkin(5) * t28 + t72;
t17 = t29 * t39 + t31 * t61;
t15 = -t29 * t61 + t31 * t39;
t12 = t28 * t29 + t30 * t63;
t11 = -t28 * t63 + t29 * t30;
t10 = t28 * t31 - t30 * t64;
t9 = t28 * t64 + t30 * t31;
t1 = [-m(2) * (g(1) * (-t41 * rSges(2,1) - rSges(2,2) * t44) + g(2) * (rSges(2,1) * t44 - t41 * rSges(2,2))) - m(3) * (g(1) * (-rSges(3,1) * t29 - rSges(3,2) * t31 - t68) + g(2) * (rSges(3,1) * t31 - rSges(3,2) * t29 + t34)) - m(4) * (g(1) * (t31 * rSges(4,3) + t56) + g(2) * (rSges(4,1) * t63 - t31 * t65 + t58) + (g(1) * (-pkin(2) - t55) + g(2) * rSges(4,3)) * t29) - m(5) * ((t17 * rSges(5,1) + t16 * rSges(5,2) + t79 * t31 + t58) * g(2) + (t15 * rSges(5,1) + t14 * rSges(5,2) + t56 + (-pkin(2) - t79) * t29) * g(1)) - m(6) * (g(1) * (rSges(6,1) * t10 + rSges(6,2) * t9 + t56) + g(2) * (t12 * rSges(6,1) + t11 * rSges(6,2) + t58) + (g(1) * t72 + g(2) * t49) * t31 + (g(1) * (-pkin(2) - t49) + g(2) * t72) * t29) - m(7) * (g(1) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t56) + g(2) * (t8 * rSges(7,1) + t7 * rSges(7,2) + t58) + (g(1) * t19 + g(2) * t48) * t31 + (g(1) * (-pkin(2) - t48) + g(2) * t19) * t29) (-m(3) - m(4) - m(5) + t75) * g(3), -m(4) * (g(3) * t55 + t76 * (-rSges(4,1) * t40 - rSges(4,2) * t43)) - m(5) * (g(3) * (t52 * t43 + t57) + t76 * (-t52 * t40 + t66 * t43)) - m(6) * (g(3) * (t51 * t43 + t77) + t76 * (-t51 * t40 + t60 * t43)) - m(7) * (g(3) * (t50 * t43 + t78) + t76 * (-t50 * t40 + t59 * t43)) -m(5) * (g(1) * (rSges(5,1) * t16 - rSges(5,2) * t17) + g(2) * (-rSges(5,1) * t14 + rSges(5,2) * t15)) - m(6) * (g(1) * (t11 * rSges(6,1) - t12 * rSges(6,2) + t16 * pkin(4)) + g(2) * (-t9 * rSges(6,1) + t10 * rSges(6,2) - t14 * pkin(4))) - m(7) * (g(1) * (-t19 * t63 + t20 * t29 + t73) + g(2) * (-t19 * t64 - t20 * t31 + t74)) + (-m(5) * (-rSges(5,1) * t39 - rSges(5,2) * t42) - m(6) * (-rSges(6,1) * t28 - rSges(6,2) * t30 - t72) - m(7) * (-t19 + t53)) * t69, t75 * (-g(3) * t43 + t76 * t40) -m(7) * (g(1) * t73 + g(2) * t74 + t53 * t69)];
taug  = t1(:);
