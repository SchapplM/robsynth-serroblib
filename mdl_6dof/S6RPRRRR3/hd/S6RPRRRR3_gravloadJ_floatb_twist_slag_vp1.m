% Calculate Gravitation load on the joints for
% S6RPRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:00:23
% EndTime: 2019-03-09 07:00:24
% DurationCPUTime: 0.70s
% Computational Cost: add. (554->139), mult. (480->190), div. (0->0), fcn. (453->12), ass. (0->68)
t45 = cos(qJ(3));
t42 = sin(qJ(3));
t74 = rSges(5,3) + pkin(8);
t61 = t74 * t42;
t87 = t45 * pkin(3) + t61;
t47 = -pkin(9) - pkin(8);
t64 = rSges(7,3) + pkin(10) - t47;
t86 = t64 * t42;
t65 = rSges(6,3) - t47;
t85 = t65 * t42;
t38 = qJ(1) + pkin(11);
t31 = sin(t38);
t32 = cos(t38);
t84 = g(1) * t32 + g(2) * t31;
t40 = qJ(4) + qJ(5);
t35 = qJ(6) + t40;
t28 = sin(t35);
t29 = cos(t35);
t71 = t31 * t45;
t5 = t28 * t71 + t29 * t32;
t6 = t28 * t32 - t29 * t71;
t83 = -t5 * rSges(7,1) + t6 * rSges(7,2);
t70 = t32 * t45;
t7 = -t28 * t70 + t29 * t31;
t8 = t28 * t31 + t29 * t70;
t82 = t7 * rSges(7,1) - t8 * rSges(7,2);
t41 = sin(qJ(4));
t81 = pkin(4) * t41;
t33 = sin(t40);
t80 = pkin(5) * t33;
t77 = g(3) * t42;
t43 = sin(qJ(1));
t76 = t43 * pkin(1);
t34 = cos(t40);
t69 = t33 * t45;
t13 = t31 * t69 + t32 * t34;
t68 = t34 * t45;
t14 = -t31 * t68 + t32 * t33;
t73 = -t13 * rSges(6,1) + t14 * rSges(6,2);
t72 = rSges(4,2) * t42;
t67 = t41 * t45;
t44 = cos(qJ(4));
t66 = t44 * t45;
t15 = t31 * t34 - t32 * t69;
t16 = t31 * t33 + t32 * t68;
t63 = t15 * rSges(6,1) - t16 * rSges(6,2);
t36 = t44 * pkin(4);
t23 = pkin(5) * t34 + t36;
t46 = cos(qJ(1));
t37 = t46 * pkin(1);
t62 = t32 * pkin(2) + t31 * pkin(7) + t37;
t60 = t32 * pkin(7) - t76;
t59 = rSges(4,1) * t45 - t72;
t57 = -rSges(6,1) * t33 - rSges(6,2) * t34;
t56 = -rSges(7,1) * t28 - rSges(7,2) * t29;
t55 = rSges(5,1) * t44 - rSges(5,2) * t41 + pkin(3);
t19 = t31 * t44 - t32 * t67;
t17 = t31 * t67 + t32 * t44;
t30 = t36 + pkin(3);
t54 = rSges(6,1) * t34 - rSges(6,2) * t33 + t30;
t21 = pkin(3) + t23;
t53 = rSges(7,1) * t29 - rSges(7,2) * t28 + t21;
t52 = t45 * t30 + t85;
t51 = t45 * t21 + t86;
t22 = t80 + t81;
t20 = t31 * t41 + t32 * t66;
t18 = -t31 * t66 + t32 * t41;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t43 - rSges(2,2) * t46) + g(2) * (rSges(2,1) * t46 - rSges(2,2) * t43)) - m(3) * (g(1) * (-rSges(3,1) * t31 - rSges(3,2) * t32 - t76) + g(2) * (rSges(3,1) * t32 - rSges(3,2) * t31 + t37)) - m(4) * (g(1) * (rSges(4,3) * t32 + t60) + g(2) * (rSges(4,1) * t70 - t32 * t72 + t62) + (g(1) * (-pkin(2) - t59) + g(2) * rSges(4,3)) * t31) - m(5) * ((rSges(5,1) * t20 + rSges(5,2) * t19 + t87 * t32 + t62) * g(2) + (rSges(5,1) * t18 + rSges(5,2) * t17 + t60 + (-pkin(2) - t87) * t31) * g(1)) - m(6) * (g(1) * (t14 * rSges(6,1) + t13 * rSges(6,2) + t60) + g(2) * (t16 * rSges(6,1) + t15 * rSges(6,2) + t62) + (g(1) * t81 + g(2) * t52) * t32 + (g(1) * (-pkin(2) - t52) + g(2) * t81) * t31) - m(7) * (g(1) * (rSges(7,1) * t6 + rSges(7,2) * t5 + t60) + g(2) * (rSges(7,1) * t8 + rSges(7,2) * t7 + t62) + (g(1) * t22 + g(2) * t51) * t32 + (g(1) * (-pkin(2) - t51) + g(2) * t22) * t31) (-m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(4) * (g(3) * t59 + t84 * (-rSges(4,1) * t42 - rSges(4,2) * t45)) - m(5) * (g(3) * (t55 * t45 + t61) + t84 * (-t55 * t42 + t74 * t45)) - m(6) * (g(3) * (t54 * t45 + t85) + t84 * (-t54 * t42 + t65 * t45)) - m(7) * (g(3) * (t53 * t45 + t86) + t84 * (-t53 * t42 + t64 * t45)) -m(5) * (g(1) * (rSges(5,1) * t19 - rSges(5,2) * t20) + g(2) * (-rSges(5,1) * t17 + rSges(5,2) * t18)) - m(6) * (g(1) * (t19 * pkin(4) + t63) + g(2) * (-t17 * pkin(4) + t73)) - m(7) * (g(1) * (-t22 * t70 + t23 * t31 + t82) + g(2) * (-t22 * t71 - t23 * t32 + t83)) + (-m(5) * (-rSges(5,1) * t41 - rSges(5,2) * t44) - m(6) * (t57 - t81) - m(7) * (-t22 + t56)) * t77, -m(6) * (g(1) * t63 + g(2) * t73) - m(7) * (g(1) * (t15 * pkin(5) + t82) + g(2) * (-t13 * pkin(5) + t83)) + (-m(6) * t57 - m(7) * (t56 - t80)) * t77, -m(7) * (g(1) * t82 + g(2) * t83 + t56 * t77)];
taug  = t1(:);
