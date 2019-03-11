% Calculate Gravitation load on the joints for
% S6RPRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-03-09 03:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRPRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRPRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:51:14
% EndTime: 2019-03-09 03:51:15
% DurationCPUTime: 0.58s
% Computational Cost: add. (449->116), mult. (415->164), div. (0->0), fcn. (385->12), ass. (0->61)
t35 = -pkin(8) - qJ(4);
t57 = rSges(6,3) - t35;
t56 = rSges(7,3) + pkin(9) - t35;
t31 = sin(pkin(11));
t36 = -pkin(7) - qJ(2);
t51 = pkin(4) * t31 - t36;
t37 = sin(qJ(1));
t65 = g(2) * t37;
t38 = cos(qJ(1));
t67 = g(1) * t38;
t72 = t65 + t67;
t53 = rSges(5,3) + qJ(4);
t29 = pkin(11) + qJ(5);
t26 = qJ(6) + t29;
t19 = cos(t26);
t30 = pkin(10) + qJ(3);
t25 = cos(t30);
t18 = sin(t26);
t62 = t37 * t18;
t6 = t19 * t38 + t25 * t62;
t61 = t37 * t19;
t7 = t18 * t38 - t25 * t61;
t71 = -t6 * rSges(7,1) + t7 * rSges(7,2);
t63 = t25 * t38;
t8 = -t18 * t63 + t61;
t9 = t19 * t63 + t62;
t70 = t8 * rSges(7,1) - t9 * rSges(7,2);
t22 = sin(t29);
t68 = pkin(5) * t22;
t34 = cos(pkin(10));
t21 = pkin(2) * t34 + pkin(1);
t16 = t38 * t21;
t66 = g(2) * t16;
t23 = sin(t30);
t64 = g(3) * t23;
t33 = cos(pkin(11));
t20 = t33 * pkin(4) + pkin(3);
t60 = t37 * t22;
t24 = cos(t29);
t59 = t37 * t24;
t58 = rSges(4,3) - t36;
t55 = t68 + t51;
t54 = rSges(3,3) + qJ(2);
t52 = -m(5) - m(6) - m(7);
t50 = g(2) * t53;
t49 = rSges(4,1) * t25 - rSges(4,2) * t23;
t47 = -rSges(7,1) * t18 - rSges(7,2) * t19;
t46 = rSges(3,1) * t34 - rSges(3,2) * sin(pkin(10)) + pkin(1);
t45 = rSges(5,1) * t33 - rSges(5,2) * t31 + pkin(3);
t12 = -t22 * t63 + t59;
t10 = t24 * t38 + t25 * t60;
t44 = rSges(5,1) * t31 + rSges(5,2) * t33 - t36;
t43 = rSges(6,1) * t24 - rSges(6,2) * t22 + t20;
t14 = pkin(5) * t24 + t20;
t42 = rSges(7,1) * t19 - rSges(7,2) * t18 + t14;
t41 = g(1) * t45;
t40 = t25 * t20 + t57 * t23;
t39 = t25 * t14 + t56 * t23;
t13 = t24 * t63 + t60;
t11 = t22 * t38 - t25 * t59;
t1 = [-m(2) * (g(1) * (-t37 * rSges(2,1) - rSges(2,2) * t38) + g(2) * (rSges(2,1) * t38 - t37 * rSges(2,2))) - m(3) * ((g(1) * t54 + g(2) * t46) * t38 + (-g(1) * t46 + g(2) * t54) * t37) - m(4) * (t66 + (g(1) * t58 + g(2) * t49) * t38 + (g(1) * (-t21 - t49) + g(2) * t58) * t37) - m(5) * (t66 + (g(2) * t45 * t25 + g(1) * t44 + t23 * t50) * t38 + (g(1) * (-t53 * t23 - t21) + g(2) * t44 - t25 * t41) * t37) - m(6) * (g(1) * (t11 * rSges(6,1) + t10 * rSges(6,2)) + g(2) * (t13 * rSges(6,1) + t12 * rSges(6,2) + t16) + (g(1) * t51 + g(2) * t40) * t38 + (g(1) * (-t21 - t40) + g(2) * t51) * t37) - m(7) * (g(1) * (t7 * rSges(7,1) + t6 * rSges(7,2)) + g(2) * (t9 * rSges(7,1) + t8 * rSges(7,2) + t16) + (g(1) * t55 + g(2) * t39) * t38 + (g(1) * (-t21 - t39) + g(2) * t55) * t37) (-m(3) - m(4) + t52) * (g(1) * t37 - g(2) * t38) -m(4) * (g(3) * t49 + t72 * (-rSges(4,1) * t23 - rSges(4,2) * t25)) - m(5) * ((g(3) * t45 + t37 * t50 + t53 * t67) * t25 + (g(3) * t53 - t38 * t41 - t45 * t65) * t23) - m(6) * ((g(3) * t43 + t72 * t57) * t25 + (g(3) * t57 - t72 * t43) * t23) - m(7) * ((g(3) * t42 + t72 * t56) * t25 + (g(3) * t56 - t72 * t42) * t23) t52 * (-g(3) * t25 + t72 * t23) -m(6) * (g(1) * (rSges(6,1) * t12 - rSges(6,2) * t13) + g(2) * (-rSges(6,1) * t10 + rSges(6,2) * t11)) - m(7) * (g(1) * (t12 * pkin(5) + t70) + g(2) * (-t10 * pkin(5) + t71)) + (-m(6) * (-rSges(6,1) * t22 - rSges(6,2) * t24) - m(7) * (t47 - t68)) * t64, -m(7) * (g(1) * t70 + g(2) * t71 + t47 * t64)];
taug  = t1(:);
