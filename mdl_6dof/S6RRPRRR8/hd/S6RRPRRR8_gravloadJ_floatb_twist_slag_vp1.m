% Calculate Gravitation load on the joints for
% S6RRPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6,theta3]';
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
% Datum: 2019-03-09 14:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRRR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRRR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:02:21
% EndTime: 2019-03-09 14:02:23
% DurationCPUTime: 0.77s
% Computational Cost: add. (593->160), mult. (575->221), div. (0->0), fcn. (549->12), ass. (0->70)
t50 = -pkin(8) - qJ(3);
t72 = rSges(5,3) - t50;
t46 = -pkin(9) + t50;
t71 = rSges(6,3) - t46;
t70 = rSges(7,3) + pkin(10) - t46;
t54 = cos(qJ(1));
t52 = sin(qJ(1));
t79 = g(2) * t52;
t84 = g(1) * t54 + t79;
t47 = pkin(11) + qJ(4);
t39 = qJ(5) + t47;
t36 = qJ(6) + t39;
t29 = sin(t36);
t30 = cos(t36);
t53 = cos(qJ(2));
t75 = t52 * t53;
t5 = t29 * t75 + t30 * t54;
t6 = t29 * t54 - t30 * t75;
t83 = -t5 * rSges(7,1) + t6 * rSges(7,2);
t74 = t53 * t54;
t7 = -t29 * t74 + t52 * t30;
t8 = t52 * t29 + t30 * t74;
t82 = t7 * rSges(7,1) - t8 * rSges(7,2);
t37 = sin(t47);
t31 = pkin(4) * t37;
t33 = sin(t39);
t81 = pkin(5) * t33;
t51 = sin(qJ(2));
t78 = g(3) * t51;
t48 = sin(pkin(11));
t41 = t48 * pkin(3);
t49 = cos(pkin(11));
t35 = t49 * pkin(3) + pkin(2);
t34 = cos(t39);
t13 = t33 * t75 + t34 * t54;
t14 = t33 * t54 - t34 * t75;
t77 = -t13 * rSges(6,1) + t14 * rSges(6,2);
t76 = rSges(3,2) * t51;
t73 = t54 * t48;
t15 = -t33 * t74 + t52 * t34;
t16 = t52 * t33 + t34 * t74;
t69 = t15 * rSges(6,1) - t16 * rSges(6,2);
t68 = t54 * pkin(1) + t52 * pkin(7);
t67 = rSges(4,3) + qJ(3);
t38 = cos(t47);
t32 = pkin(4) * t38;
t26 = t32 + t35;
t24 = -t31 - t81;
t66 = t54 * t67;
t65 = rSges(3,1) * t53 - t76;
t63 = -rSges(6,1) * t33 - rSges(6,2) * t34;
t62 = -rSges(7,1) * t29 - rSges(7,2) * t30;
t61 = rSges(4,1) * t49 - rSges(4,2) * t48 + pkin(2);
t20 = -t37 * t74 + t52 * t38;
t18 = t37 * t75 + t38 * t54;
t60 = rSges(5,1) * t38 - rSges(5,2) * t37 + t35;
t59 = rSges(6,1) * t34 - rSges(6,2) * t33 + t26;
t28 = pkin(5) * t34;
t17 = t28 + t26;
t58 = rSges(7,1) * t30 - rSges(7,2) * t29 + t17;
t57 = t35 * t53 + t72 * t51;
t56 = t26 * t53 + t71 * t51;
t55 = t17 * t53 + t70 * t51;
t44 = t54 * pkin(7);
t27 = t31 + t41;
t25 = t28 + t32;
t22 = -t24 + t41;
t21 = t52 * t37 + t38 * t74;
t19 = t37 * t54 - t38 * t75;
t1 = [-m(2) * (g(1) * (-t52 * rSges(2,1) - rSges(2,2) * t54) + g(2) * (rSges(2,1) * t54 - t52 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t54 + t44) + g(2) * (rSges(3,1) * t74 - t54 * t76 + t68) + (g(1) * (-pkin(1) - t65) + g(2) * rSges(3,3)) * t52) - m(4) * (g(1) * (-pkin(2) * t75 - t52 * pkin(1) + t44 + (-t49 * t75 + t73) * rSges(4,1) + (t48 * t75 + t49 * t54) * rSges(4,2)) + g(2) * (pkin(2) * t74 + (t52 * t48 + t49 * t74) * rSges(4,1) + (t52 * t49 - t53 * t73) * rSges(4,2) + t68) + (-g(1) * t67 * t52 + g(2) * t66) * t51) - m(5) * (g(1) * (t19 * rSges(5,1) + t18 * rSges(5,2) + t44) + g(2) * (t21 * rSges(5,1) + t20 * rSges(5,2) + t68) + (g(1) * t41 + g(2) * t57) * t54 + (g(1) * (-pkin(1) - t57) + g(2) * t41) * t52) - m(6) * (g(1) * (t14 * rSges(6,1) + t13 * rSges(6,2) + t44) + g(2) * (t16 * rSges(6,1) + t15 * rSges(6,2) + t68) + (g(1) * t27 + g(2) * t56) * t54 + (g(1) * (-pkin(1) - t56) + g(2) * t27) * t52) - m(7) * (g(1) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t44) + g(2) * (t8 * rSges(7,1) + t7 * rSges(7,2) + t68) + (g(1) * t22 + g(2) * t55) * t54 + (g(1) * (-pkin(1) - t55) + g(2) * t22) * t52) -m(3) * (g(3) * t65 + t84 * (-rSges(3,1) * t51 - rSges(3,2) * t53)) - m(4) * ((g(1) * t66 + g(3) * t61 + t67 * t79) * t53 + (g(3) * t67 - t84 * t61) * t51) - m(5) * ((g(3) * t60 + t84 * t72) * t53 + (g(3) * t72 - t84 * t60) * t51) - m(6) * ((g(3) * t59 + t84 * t71) * t53 + (g(3) * t71 - t84 * t59) * t51) - m(7) * ((g(3) * t58 + t84 * t70) * t53 + (g(3) * t70 - t84 * t58) * t51) (-m(4) - m(5) - m(6) - m(7)) * (-g(3) * t53 + t84 * t51) -m(5) * (g(1) * (rSges(5,1) * t20 - rSges(5,2) * t21) + g(2) * (-rSges(5,1) * t18 + rSges(5,2) * t19)) - m(6) * (g(1) * (pkin(4) * t20 + t69) + g(2) * (-pkin(4) * t18 + t77)) - m(7) * (g(1) * (t24 * t74 + t52 * t25 + t82) + g(2) * (t24 * t75 - t25 * t54 + t83)) + (-m(5) * (-rSges(5,1) * t37 - rSges(5,2) * t38) - m(6) * (t63 - t31) - m(7) * (t24 + t62)) * t78, -m(6) * (g(1) * t69 + g(2) * t77) - m(7) * (g(1) * (pkin(5) * t15 + t82) + g(2) * (-pkin(5) * t13 + t83)) + (-m(6) * t63 - m(7) * (t62 - t81)) * t78, -m(7) * (g(1) * t82 + g(2) * t83 + t62 * t78)];
taug  = t1(:);
