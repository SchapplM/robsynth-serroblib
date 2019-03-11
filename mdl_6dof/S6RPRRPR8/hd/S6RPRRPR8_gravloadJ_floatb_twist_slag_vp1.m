% Calculate Gravitation load on the joints for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
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
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:22:41
% EndTime: 2019-03-09 05:22:43
% DurationCPUTime: 0.80s
% Computational Cost: add. (341->132), mult. (444->178), div. (0->0), fcn. (415->10), ass. (0->69)
t38 = sin(qJ(3));
t41 = cos(qJ(3));
t70 = rSges(5,3) + pkin(8);
t83 = t70 * t41;
t86 = t38 * pkin(3) - t83;
t39 = sin(qJ(1));
t42 = cos(qJ(1));
t51 = g(1) * t39 - g(2) * t42;
t36 = -qJ(5) - pkin(8);
t54 = rSges(7,3) + pkin(9) - t36;
t85 = t54 * t41;
t55 = rSges(6,3) - t36;
t84 = t55 * t41;
t79 = -m(6) - m(7);
t78 = -pkin(1) - pkin(7);
t35 = qJ(4) + pkin(10);
t28 = qJ(6) + t35;
t24 = cos(t28);
t60 = t42 * t24;
t23 = sin(t28);
t68 = t39 * t23;
t5 = -t38 * t68 + t60;
t61 = t42 * t23;
t67 = t39 * t24;
t6 = t38 * t67 + t61;
t77 = t5 * rSges(7,1) - t6 * rSges(7,2);
t7 = t38 * t61 + t67;
t8 = t38 * t60 - t68;
t76 = t7 * rSges(7,1) + t8 * rSges(7,2);
t73 = g(3) * t41;
t37 = sin(qJ(4));
t72 = t37 * pkin(4);
t69 = t38 * t42;
t26 = sin(t35);
t66 = t39 * t26;
t27 = cos(t35);
t65 = t39 * t27;
t64 = t39 * t37;
t40 = cos(qJ(4));
t63 = t39 * t40;
t62 = t41 * rSges(4,2);
t59 = t42 * t26;
t58 = t42 * t27;
t57 = t42 * t37;
t56 = t42 * t40;
t31 = t40 * pkin(4);
t20 = pkin(5) * t27 + t31;
t53 = t42 * pkin(1) + t39 * qJ(2);
t52 = t42 * pkin(7) + t53;
t49 = t38 * rSges(4,1) + t62;
t48 = -rSges(7,1) * t23 - rSges(7,2) * t24;
t47 = rSges(5,1) * t40 - rSges(5,2) * t37 + pkin(3);
t16 = t38 * t57 + t63;
t14 = -t38 * t64 + t56;
t25 = t31 + pkin(3);
t46 = rSges(6,1) * t27 - rSges(6,2) * t26 + t25;
t18 = pkin(3) + t20;
t45 = rSges(7,1) * t24 - rSges(7,2) * t23 + t18;
t44 = t38 * t25 - t84;
t43 = t38 * t18 - t85;
t30 = t42 * qJ(2);
t19 = pkin(5) * t26 + t72;
t17 = t38 * t56 - t64;
t15 = t38 * t63 + t57;
t12 = t38 * t58 - t66;
t11 = t38 * t59 + t65;
t10 = t38 * t65 + t59;
t9 = -t38 * t66 + t58;
t1 = [-m(2) * (g(1) * (-t39 * rSges(2,1) - t42 * rSges(2,2)) + g(2) * (t42 * rSges(2,1) - t39 * rSges(2,2))) - m(3) * (g(1) * (t42 * rSges(3,3) + t30 + (rSges(3,2) - pkin(1)) * t39) + g(2) * (-t42 * rSges(3,2) + t39 * rSges(3,3) + t53)) - m(4) * (g(1) * (rSges(4,1) * t69 + t42 * t62 + t30) + g(2) * (t42 * rSges(4,3) + t52) + (g(1) * (-rSges(4,3) + t78) + g(2) * t49) * t39) - m(5) * ((t15 * rSges(5,1) + t14 * rSges(5,2) + t86 * t39 + t52) * g(2) + (t17 * rSges(5,1) - t16 * rSges(5,2) + t78 * t39 + t86 * t42 + t30) * g(1)) - m(6) * (g(1) * (t12 * rSges(6,1) - t11 * rSges(6,2) + t30) + g(2) * (t10 * rSges(6,1) + t9 * rSges(6,2) + t52) + (g(1) * t44 + g(2) * t72) * t42 + (g(1) * (-t72 + t78) + g(2) * t44) * t39) - m(7) * (g(1) * (t8 * rSges(7,1) - t7 * rSges(7,2) + t30) + g(2) * (t6 * rSges(7,1) + t5 * rSges(7,2) + t52) + (g(1) * t43 + g(2) * t19) * t42 + (g(1) * (-t19 + t78) + g(2) * t43) * t39) (-m(3) - m(4) - m(5) + t79) * t51, -m(4) * (-g(3) * t49 + t51 * (rSges(4,1) * t41 - rSges(4,2) * t38)) - m(5) * (g(3) * (-t47 * t38 + t83) + t51 * (t70 * t38 + t47 * t41)) - m(6) * (g(3) * (-t46 * t38 + t84) + t51 * (t55 * t38 + t46 * t41)) - m(7) * (g(3) * (-t45 * t38 + t85) + t51 * (t54 * t38 + t45 * t41)) -m(5) * (g(1) * (t14 * rSges(5,1) - t15 * rSges(5,2)) + g(2) * (t16 * rSges(5,1) + t17 * rSges(5,2))) - m(6) * (g(1) * (t9 * rSges(6,1) - t10 * rSges(6,2) + t14 * pkin(4)) + g(2) * (t11 * rSges(6,1) + t12 * rSges(6,2) + t16 * pkin(4))) - m(7) * (g(1) * (-t39 * t38 * t19 + t42 * t20 + t77) + g(2) * (t19 * t69 + t39 * t20 + t76)) + (-m(5) * (-rSges(5,1) * t37 - rSges(5,2) * t40) - m(6) * (-rSges(6,1) * t26 - rSges(6,2) * t27 - t72) - m(7) * (-t19 + t48)) * t73, t79 * (g(3) * t38 - t51 * t41) -m(7) * (g(1) * t77 + g(2) * t76 + t48 * t73)];
taug  = t1(:);
