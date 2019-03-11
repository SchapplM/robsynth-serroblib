% Calculate Gravitation load on the joints for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:15:52
% EndTime: 2019-03-09 10:15:54
% DurationCPUTime: 0.70s
% Computational Cost: add. (488->143), mult. (481->188), div. (0->0), fcn. (449->12), ass. (0->74)
t37 = qJ(2) + pkin(10);
t29 = sin(t37);
t75 = rSges(5,3) + pkin(8);
t88 = t75 * t29;
t38 = -qJ(5) - pkin(8);
t63 = rSges(7,3) + pkin(9) - t38;
t87 = t63 * t29;
t64 = rSges(6,3) - t38;
t86 = t64 * t29;
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t85 = g(1) * t45 + g(2) * t42;
t84 = -m(6) - m(7);
t36 = qJ(4) + pkin(11);
t32 = qJ(6) + t36;
t25 = cos(t32);
t31 = cos(t37);
t24 = sin(t32);
t72 = t42 * t24;
t6 = t25 * t45 + t31 * t72;
t71 = t42 * t25;
t7 = t24 * t45 - t31 * t71;
t83 = -t6 * rSges(7,1) + t7 * rSges(7,2);
t74 = t31 * t45;
t8 = -t24 * t74 + t71;
t9 = t25 * t74 + t72;
t82 = t8 * rSges(7,1) - t9 * rSges(7,2);
t41 = sin(qJ(2));
t81 = pkin(2) * t41;
t40 = sin(qJ(4));
t80 = pkin(4) * t40;
t77 = g(3) * t29;
t76 = rSges(3,3) + pkin(7);
t73 = t40 * t45;
t28 = sin(t36);
t70 = t42 * t28;
t30 = cos(t36);
t69 = t42 * t30;
t68 = t42 * t40;
t43 = cos(qJ(4));
t67 = t42 * t43;
t66 = t43 * t45;
t39 = -qJ(3) - pkin(7);
t65 = rSges(4,3) - t39;
t19 = pkin(5) * t28 + t80;
t62 = t19 - t39;
t33 = t43 * pkin(4);
t20 = pkin(5) * t30 + t33;
t61 = -t39 + t80;
t44 = cos(qJ(2));
t60 = rSges(3,1) * t44 - rSges(3,2) * t41;
t58 = rSges(4,1) * t31 - rSges(4,2) * t29;
t57 = -rSges(7,1) * t24 - rSges(7,2) * t25;
t56 = pkin(1) + t60;
t55 = rSges(5,1) * t43 - rSges(5,2) * t40 + pkin(3);
t16 = -t31 * t73 + t67;
t14 = t31 * t68 + t66;
t26 = t33 + pkin(3);
t54 = rSges(6,1) * t30 - rSges(6,2) * t28 + t26;
t18 = pkin(3) + t20;
t53 = rSges(7,1) * t25 - rSges(7,2) * t24 + t18;
t52 = pkin(3) * t31 + t88;
t50 = t31 * t26 + t86;
t49 = t31 * t18 + t87;
t34 = t44 * pkin(2);
t27 = t34 + pkin(1);
t22 = t45 * t27;
t17 = t31 * t66 + t68;
t15 = -t31 * t67 + t73;
t13 = t30 * t74 + t70;
t12 = -t28 * t74 + t69;
t11 = t28 * t45 - t31 * t69;
t10 = t30 * t45 + t31 * t70;
t1 = [-m(2) * (g(1) * (-t42 * rSges(2,1) - rSges(2,2) * t45) + g(2) * (rSges(2,1) * t45 - t42 * rSges(2,2))) - m(3) * ((g(1) * t76 + g(2) * t56) * t45 + (-g(1) * t56 + g(2) * t76) * t42) - m(4) * (g(2) * t22 + (g(1) * t65 + g(2) * t58) * t45 + (g(1) * (-t27 - t58) + g(2) * t65) * t42) - m(5) * (g(1) * (t15 * rSges(5,1) + t14 * rSges(5,2)) + g(2) * (t17 * rSges(5,1) + t16 * rSges(5,2) + t22) + (-g(1) * t39 + g(2) * t52) * t45 + (g(1) * (-t27 - t52) - g(2) * t39) * t42) - m(6) * (g(1) * (t11 * rSges(6,1) + t10 * rSges(6,2)) + g(2) * (t13 * rSges(6,1) + t12 * rSges(6,2) + t22) + (g(1) * t61 + g(2) * t50) * t45 + (g(1) * (-t27 - t50) + g(2) * t61) * t42) - m(7) * (g(1) * (t7 * rSges(7,1) + t6 * rSges(7,2)) + g(2) * (t9 * rSges(7,1) + t8 * rSges(7,2) + t22) + (g(1) * t62 + g(2) * t49) * t45 + (g(1) * (-t27 - t49) + g(2) * t62) * t42) -m(3) * (g(3) * t60 + t85 * (-rSges(3,1) * t41 - rSges(3,2) * t44)) - m(4) * (g(3) * (t34 + t58) + t85 * (-rSges(4,1) * t29 - rSges(4,2) * t31 - t81)) - m(5) * (g(3) * (t55 * t31 + t34 + t88) + t85 * (-t55 * t29 + t75 * t31 - t81)) - m(6) * (g(3) * (t54 * t31 + t34 + t86) + t85 * (-t54 * t29 + t64 * t31 - t81)) - m(7) * (g(3) * (t53 * t31 + t34 + t87) + t85 * (-t53 * t29 + t63 * t31 - t81)) (-m(4) - m(5) + t84) * (g(1) * t42 - g(2) * t45) -m(5) * (g(1) * (rSges(5,1) * t16 - rSges(5,2) * t17) + g(2) * (-rSges(5,1) * t14 + rSges(5,2) * t15)) - m(6) * (g(1) * (t12 * rSges(6,1) - t13 * rSges(6,2) + t16 * pkin(4)) + g(2) * (-t10 * rSges(6,1) + t11 * rSges(6,2) - t14 * pkin(4))) - m(7) * (g(1) * (-t19 * t74 + t42 * t20 + t82) + g(2) * (-t42 * t31 * t19 - t20 * t45 + t83)) + (-m(5) * (-rSges(5,1) * t40 - rSges(5,2) * t43) - m(6) * (-rSges(6,1) * t28 - rSges(6,2) * t30 - t80) - m(7) * (-t19 + t57)) * t77, t84 * (-g(3) * t31 + t85 * t29) -m(7) * (g(1) * t82 + g(2) * t83 + t57 * t77)];
taug  = t1(:);
