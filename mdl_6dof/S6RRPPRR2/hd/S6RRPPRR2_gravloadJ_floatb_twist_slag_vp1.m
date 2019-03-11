% Calculate Gravitation load on the joints for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:49:51
% EndTime: 2019-03-09 08:49:53
% DurationCPUTime: 0.68s
% Computational Cost: add. (465->122), mult. (446->160), div. (0->0), fcn. (413->12), ass. (0->64)
t31 = qJ(2) + pkin(10);
t23 = sin(t31);
t25 = cos(t31);
t32 = sin(pkin(11));
t33 = cos(pkin(11));
t50 = rSges(5,1) * t33 - rSges(5,2) * t32 + pkin(3);
t59 = rSges(5,3) + qJ(4);
t80 = t59 * t23 + t50 * t25;
t35 = -pkin(8) - qJ(4);
t61 = rSges(7,3) + pkin(9) - t35;
t82 = t61 * t23;
t62 = rSges(6,3) - t35;
t81 = t62 * t23;
t34 = -qJ(3) - pkin(7);
t57 = pkin(4) * t32 - t34;
t37 = sin(qJ(1));
t39 = cos(qJ(1));
t79 = g(1) * t39 + g(2) * t37;
t30 = pkin(11) + qJ(5);
t26 = qJ(6) + t30;
t19 = cos(t26);
t18 = sin(t26);
t67 = t37 * t18;
t6 = t19 * t39 + t25 * t67;
t66 = t37 * t19;
t7 = t18 * t39 - t25 * t66;
t78 = -t6 * rSges(7,1) + t7 * rSges(7,2);
t68 = t25 * t39;
t8 = -t18 * t68 + t66;
t9 = t19 * t68 + t67;
t77 = t8 * rSges(7,1) - t9 * rSges(7,2);
t36 = sin(qJ(2));
t76 = pkin(2) * t36;
t22 = sin(t30);
t74 = pkin(5) * t22;
t38 = cos(qJ(2));
t28 = t38 * pkin(2);
t21 = t28 + pkin(1);
t17 = t39 * t21;
t72 = g(2) * t17;
t70 = g(3) * t23;
t69 = rSges(3,3) + pkin(7);
t20 = t33 * pkin(4) + pkin(3);
t65 = t37 * t22;
t24 = cos(t30);
t64 = t37 * t24;
t63 = rSges(4,3) - t34;
t60 = t74 + t57;
t58 = -m(5) - m(6) - m(7);
t55 = rSges(3,1) * t38 - rSges(3,2) * t36;
t53 = rSges(4,1) * t25 - rSges(4,2) * t23;
t52 = -rSges(7,1) * t18 - rSges(7,2) * t19;
t51 = pkin(1) + t55;
t12 = -t22 * t68 + t64;
t10 = t24 * t39 + t25 * t65;
t49 = rSges(5,1) * t32 + rSges(5,2) * t33 - t34;
t48 = rSges(6,1) * t24 - rSges(6,2) * t22 + t20;
t14 = pkin(5) * t24 + t20;
t47 = rSges(7,1) * t19 - rSges(7,2) * t18 + t14;
t44 = t25 * t20 + t81;
t43 = t25 * t14 + t82;
t13 = t24 * t68 + t65;
t11 = t22 * t39 - t25 * t64;
t1 = [-m(2) * (g(1) * (-t37 * rSges(2,1) - rSges(2,2) * t39) + g(2) * (rSges(2,1) * t39 - t37 * rSges(2,2))) - m(3) * ((g(1) * t69 + g(2) * t51) * t39 + (-g(1) * t51 + g(2) * t69) * t37) - m(4) * (t72 + (g(1) * t63 + g(2) * t53) * t39 + (g(1) * (-t21 - t53) + g(2) * t63) * t37) - m(5) * (t72 + (g(1) * t49 + g(2) * t80) * t39 + (g(2) * t49 + (-t21 - t80) * g(1)) * t37) - m(6) * (g(1) * (t11 * rSges(6,1) + t10 * rSges(6,2)) + g(2) * (t13 * rSges(6,1) + t12 * rSges(6,2) + t17) + (g(1) * t57 + g(2) * t44) * t39 + (g(1) * (-t21 - t44) + g(2) * t57) * t37) - m(7) * (g(1) * (t7 * rSges(7,1) + t6 * rSges(7,2)) + g(2) * (t9 * rSges(7,1) + t8 * rSges(7,2) + t17) + (g(1) * t60 + g(2) * t43) * t39 + (g(1) * (-t21 - t43) + g(2) * t60) * t37) -m(3) * (g(3) * t55 + t79 * (-rSges(3,1) * t36 - rSges(3,2) * t38)) - m(4) * (g(3) * (t28 + t53) + t79 * (-rSges(4,1) * t23 - rSges(4,2) * t25 - t76)) - m(5) * (g(3) * (t28 + t80) + t79 * (-t50 * t23 + t59 * t25 - t76)) - m(6) * (g(3) * (t48 * t25 + t28 + t81) + t79 * (-t48 * t23 + t62 * t25 - t76)) - m(7) * (g(3) * (t47 * t25 + t28 + t82) + t79 * (-t47 * t23 + t61 * t25 - t76)) (-m(4) + t58) * (g(1) * t37 - g(2) * t39) t58 * (-g(3) * t25 + t23 * t79) -m(6) * (g(1) * (rSges(6,1) * t12 - rSges(6,2) * t13) + g(2) * (-rSges(6,1) * t10 + rSges(6,2) * t11)) - m(7) * (g(1) * (t12 * pkin(5) + t77) + g(2) * (-t10 * pkin(5) + t78)) + (-m(6) * (-rSges(6,1) * t22 - rSges(6,2) * t24) - m(7) * (t52 - t74)) * t70, -m(7) * (g(1) * t77 + g(2) * t78 + t52 * t70)];
taug  = t1(:);
