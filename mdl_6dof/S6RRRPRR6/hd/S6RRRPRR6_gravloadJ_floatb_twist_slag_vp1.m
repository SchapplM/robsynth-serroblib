% Calculate Gravitation load on the joints for
% S6RRRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRRPRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRRPRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 18:25:09
% EndTime: 2019-03-09 18:25:11
% DurationCPUTime: 0.80s
% Computational Cost: add. (615->168), mult. (615->230), div. (0->0), fcn. (590->12), ass. (0->70)
t50 = -qJ(4) - pkin(8);
t72 = rSges(5,3) - t50;
t48 = -pkin(9) + t50;
t71 = rSges(6,3) - t48;
t70 = rSges(7,3) + pkin(10) - t48;
t53 = sin(qJ(1));
t56 = cos(qJ(1));
t87 = g(1) * t56 + g(2) * t53;
t52 = sin(qJ(2));
t55 = cos(qJ(2));
t77 = rSges(4,3) + pkin(8);
t86 = t55 * pkin(2) + t52 * t77;
t49 = qJ(3) + pkin(11);
t42 = qJ(5) + t49;
t38 = qJ(6) + t42;
t33 = sin(t38);
t34 = cos(t38);
t74 = t53 * t55;
t5 = t33 * t74 + t34 * t56;
t6 = t33 * t56 - t34 * t74;
t85 = -rSges(7,1) * t5 + rSges(7,2) * t6;
t73 = t55 * t56;
t7 = -t33 * t73 + t34 * t53;
t8 = t33 * t53 + t34 * t73;
t84 = rSges(7,1) * t7 - rSges(7,2) * t8;
t51 = sin(qJ(3));
t83 = pkin(3) * t51;
t36 = sin(t42);
t82 = pkin(5) * t36;
t79 = g(3) * t52;
t37 = cos(t42);
t13 = t36 * t74 + t37 * t56;
t14 = t36 * t56 - t37 * t74;
t76 = -rSges(6,1) * t13 + rSges(6,2) * t14;
t75 = rSges(3,2) * t52;
t15 = -t36 * t73 + t37 * t53;
t16 = t36 * t53 + t37 * t73;
t69 = rSges(6,1) * t15 - rSges(6,2) * t16;
t41 = cos(t49);
t54 = cos(qJ(3));
t45 = t54 * pkin(3);
t31 = pkin(4) * t41 + t45;
t68 = pkin(1) * t56 + pkin(7) * t53;
t23 = pkin(5) * t37 + t31;
t40 = sin(t49);
t30 = pkin(4) * t40 + t83;
t67 = rSges(3,1) * t55 - t75;
t65 = -rSges(6,1) * t36 - rSges(6,2) * t37;
t64 = -rSges(7,1) * t33 - rSges(7,2) * t34;
t63 = rSges(4,1) * t54 - rSges(4,2) * t51 + pkin(2);
t27 = -t51 * t73 + t53 * t54;
t25 = t51 * t74 + t54 * t56;
t39 = t45 + pkin(2);
t62 = rSges(5,1) * t41 - rSges(5,2) * t40 + t39;
t29 = pkin(2) + t31;
t61 = rSges(6,1) * t37 - rSges(6,2) * t36 + t29;
t17 = pkin(2) + t23;
t60 = rSges(7,1) * t34 - rSges(7,2) * t33 + t17;
t59 = t39 * t55 + t52 * t72;
t58 = t29 * t55 + t52 * t71;
t57 = t17 * t55 + t52 * t70;
t46 = t56 * pkin(7);
t28 = t51 * t53 + t54 * t73;
t26 = t51 * t56 - t54 * t74;
t22 = t30 + t82;
t21 = t40 * t53 + t41 * t73;
t20 = -t40 * t73 + t41 * t53;
t19 = t40 * t56 - t41 * t74;
t18 = t40 * t74 + t41 * t56;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t53 - rSges(2,2) * t56) + g(2) * (rSges(2,1) * t56 - rSges(2,2) * t53)) - m(3) * (g(1) * (t56 * rSges(3,3) + t46) + g(2) * (rSges(3,1) * t73 - t56 * t75 + t68) + (g(1) * (-pkin(1) - t67) + g(2) * rSges(3,3)) * t53) - m(4) * ((t28 * rSges(4,1) + t27 * rSges(4,2) + t56 * t86 + t68) * g(2) + (rSges(4,1) * t26 + rSges(4,2) * t25 + t46 + (-pkin(1) - t86) * t53) * g(1)) - m(5) * (g(1) * (t19 * rSges(5,1) + t18 * rSges(5,2) + t46) + g(2) * (t21 * rSges(5,1) + t20 * rSges(5,2) + t68) + (g(1) * t83 + g(2) * t59) * t56 + (g(1) * (-pkin(1) - t59) + g(2) * t83) * t53) - m(6) * (g(1) * (t14 * rSges(6,1) + t13 * rSges(6,2) + t46) + g(2) * (t16 * rSges(6,1) + t15 * rSges(6,2) + t68) + (g(1) * t30 + g(2) * t58) * t56 + (g(1) * (-pkin(1) - t58) + g(2) * t30) * t53) - m(7) * (g(1) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t46) + g(2) * (t8 * rSges(7,1) + t7 * rSges(7,2) + t68) + (g(1) * t22 + g(2) * t57) * t56 + (g(1) * (-pkin(1) - t57) + g(2) * t22) * t53) -m(3) * (g(3) * t67 + t87 * (-rSges(3,1) * t52 - rSges(3,2) * t55)) - m(4) * ((g(3) * t63 + t77 * t87) * t55 + (g(3) * t77 - t63 * t87) * t52) - m(5) * ((g(3) * t62 + t72 * t87) * t55 + (g(3) * t72 - t62 * t87) * t52) - m(6) * ((g(3) * t61 + t71 * t87) * t55 + (g(3) * t71 - t61 * t87) * t52) - m(7) * ((g(3) * t60 + t70 * t87) * t55 + (g(3) * t70 - t60 * t87) * t52) -m(4) * (g(1) * (rSges(4,1) * t27 - rSges(4,2) * t28) + g(2) * (-rSges(4,1) * t25 + rSges(4,2) * t26) + (-rSges(4,1) * t51 - rSges(4,2) * t54) * t79) - m(5) * (g(1) * (t20 * rSges(5,1) - t21 * rSges(5,2)) + g(2) * (-t18 * rSges(5,1) + t19 * rSges(5,2)) + (g(1) * t27 - g(2) * t25) * pkin(3) + (-rSges(5,1) * t40 - rSges(5,2) * t41 - t83) * t79) - m(6) * (g(1) * (-t30 * t73 + t31 * t53 + t69) + g(2) * (-t30 * t74 - t31 * t56 + t76) + (-t30 + t65) * t79) - m(7) * (g(1) * (-t22 * t73 + t23 * t53 + t84) + g(2) * (-t22 * t74 - t23 * t56 + t85) + (-t22 + t64) * t79) (-m(5) - m(6) - m(7)) * (-g(3) * t55 + t52 * t87) -m(6) * (g(1) * t69 + g(2) * t76) - m(7) * (g(1) * (pkin(5) * t15 + t84) + g(2) * (-pkin(5) * t13 + t85)) + (-m(6) * t65 - m(7) * (t64 - t82)) * t79, -m(7) * (g(1) * t84 + g(2) * t85 + t64 * t79)];
taug  = t1(:);
