% Calculate Gravitation load on the joints for
% S6RPRRRR6
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
% Datum: 2019-03-09 07:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRRR6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRRR6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:12:16
% EndTime: 2019-03-09 07:12:17
% DurationCPUTime: 0.65s
% Computational Cost: add. (535->140), mult. (498->193), div. (0->0), fcn. (471->12), ass. (0->71)
t77 = rSges(5,3) + pkin(8);
t47 = -pkin(9) - pkin(8);
t64 = rSges(6,3) - t47;
t63 = rSges(7,3) + pkin(10) - t47;
t44 = sin(qJ(1));
t46 = cos(qJ(1));
t85 = g(1) * t46 + g(2) * t44;
t39 = qJ(4) + qJ(5);
t35 = qJ(6) + t39;
t28 = sin(t35);
t29 = cos(t35);
t37 = pkin(11) + qJ(3);
t32 = cos(t37);
t75 = t32 * t44;
t5 = t28 * t75 + t29 * t46;
t6 = t28 * t46 - t29 * t75;
t84 = -t5 * rSges(7,1) + t6 * rSges(7,2);
t74 = t32 * t46;
t7 = -t28 * t74 + t29 * t44;
t8 = t28 * t44 + t29 * t74;
t83 = t7 * rSges(7,1) - t8 * rSges(7,2);
t43 = sin(qJ(4));
t82 = pkin(4) * t43;
t33 = sin(t39);
t81 = pkin(5) * t33;
t31 = sin(t37);
t78 = g(3) * t31;
t34 = cos(t39);
t70 = t34 * t46;
t73 = t33 * t44;
t13 = t32 * t73 + t70;
t71 = t34 * t44;
t72 = t33 * t46;
t14 = -t32 * t71 + t72;
t76 = -t13 * rSges(6,1) + t14 * rSges(6,2);
t69 = t43 * t44;
t68 = t43 * t46;
t45 = cos(qJ(4));
t67 = t44 * t45;
t66 = t45 * t46;
t42 = -pkin(7) - qJ(2);
t65 = rSges(4,3) - t42;
t15 = -t32 * t72 + t71;
t16 = t32 * t70 + t73;
t62 = t15 * rSges(6,1) - t16 * rSges(6,2);
t22 = t81 + t82;
t61 = t22 - t42;
t36 = t45 * pkin(4);
t23 = pkin(5) * t34 + t36;
t60 = rSges(3,3) + qJ(2);
t59 = -t42 + t82;
t58 = rSges(4,1) * t32 - t31 * rSges(4,2);
t56 = -rSges(6,1) * t33 - rSges(6,2) * t34;
t55 = -rSges(7,1) * t28 - rSges(7,2) * t29;
t41 = cos(pkin(11));
t54 = rSges(3,1) * t41 - rSges(3,2) * sin(pkin(11)) + pkin(1);
t53 = rSges(5,1) * t45 - rSges(5,2) * t43 + pkin(3);
t19 = -t32 * t68 + t67;
t17 = t32 * t69 + t66;
t30 = t36 + pkin(3);
t52 = rSges(6,1) * t34 - rSges(6,2) * t33 + t30;
t21 = pkin(3) + t23;
t51 = rSges(7,1) * t29 - rSges(7,2) * t28 + t21;
t50 = pkin(3) * t32 + t31 * t77;
t49 = t32 * t30 + t31 * t64;
t48 = t21 * t32 + t31 * t63;
t26 = pkin(2) * t41 + pkin(1);
t24 = t46 * t26;
t20 = t32 * t66 + t69;
t18 = -t32 * t67 + t68;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t44 - rSges(2,2) * t46) + g(2) * (rSges(2,1) * t46 - rSges(2,2) * t44)) - m(3) * ((g(1) * t60 + g(2) * t54) * t46 + (-g(1) * t54 + g(2) * t60) * t44) - m(4) * (g(2) * t24 + (g(1) * t65 + g(2) * t58) * t46 + (g(1) * (-t26 - t58) + g(2) * t65) * t44) - m(5) * (g(1) * (rSges(5,1) * t18 + rSges(5,2) * t17) + g(2) * (rSges(5,1) * t20 + rSges(5,2) * t19 + t24) + (-g(1) * t42 + g(2) * t50) * t46 + (g(1) * (-t26 - t50) - g(2) * t42) * t44) - m(6) * (g(1) * (t14 * rSges(6,1) + t13 * rSges(6,2)) + g(2) * (t16 * rSges(6,1) + t15 * rSges(6,2) + t24) + (g(1) * t59 + g(2) * t49) * t46 + (g(1) * (-t26 - t49) + g(2) * t59) * t44) - m(7) * (g(1) * (rSges(7,1) * t6 + rSges(7,2) * t5) + g(2) * (rSges(7,1) * t8 + rSges(7,2) * t7 + t24) + (g(1) * t61 + g(2) * t48) * t46 + (g(1) * (-t26 - t48) + g(2) * t61) * t44) (-m(3) - m(4) - m(5) - m(6) - m(7)) * (g(1) * t44 - g(2) * t46) -m(4) * (g(3) * t58 + t85 * (-rSges(4,1) * t31 - rSges(4,2) * t32)) - m(5) * ((g(3) * t53 + t77 * t85) * t32 + (g(3) * t77 - t53 * t85) * t31) - m(6) * ((g(3) * t52 + t64 * t85) * t32 + (g(3) * t64 - t52 * t85) * t31) - m(7) * ((g(3) * t51 + t63 * t85) * t32 + (g(3) * t63 - t51 * t85) * t31) -m(5) * (g(1) * (rSges(5,1) * t19 - rSges(5,2) * t20) + g(2) * (-rSges(5,1) * t17 + rSges(5,2) * t18)) - m(6) * (g(1) * (t19 * pkin(4) + t62) + g(2) * (-t17 * pkin(4) + t76)) - m(7) * (g(1) * (-t22 * t74 + t23 * t44 + t83) + g(2) * (-t22 * t75 - t23 * t46 + t84)) + (-m(5) * (-rSges(5,1) * t43 - rSges(5,2) * t45) - m(6) * (t56 - t82) - m(7) * (-t22 + t55)) * t78, -m(6) * (g(1) * t62 + g(2) * t76) - m(7) * (g(1) * (t15 * pkin(5) + t83) + g(2) * (-t13 * pkin(5) + t84)) + (-m(6) * t56 - m(7) * (t55 - t81)) * t78, -m(7) * (g(1) * t83 + g(2) * t84 + t55 * t78)];
taug  = t1(:);
