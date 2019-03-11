% Calculate Gravitation load on the joints for
% S6RPRRPR7
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
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR7_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:19:54
% EndTime: 2019-03-09 05:19:55
% DurationCPUTime: 0.67s
% Computational Cost: add. (344->126), mult. (344->155), div. (0->0), fcn. (295->10), ass. (0->69)
t35 = qJ(3) + qJ(4);
t28 = pkin(10) + t35;
t27 = cos(t28);
t88 = rSges(7,3) + pkin(9);
t89 = t88 * t27;
t26 = sin(t28);
t87 = t26 * t88;
t41 = cos(qJ(1));
t77 = g(2) * t41;
t38 = sin(qJ(1));
t86 = g(1) * t38 - t77;
t85 = g(1) * t41 + g(2) * t38;
t36 = sin(qJ(6));
t68 = rSges(7,2) * t36;
t54 = t27 * t68;
t84 = t54 * t77;
t83 = -m(6) - m(7);
t42 = -pkin(8) - pkin(7);
t29 = sin(t35);
t82 = pkin(4) * t29;
t30 = cos(t35);
t81 = pkin(4) * t30;
t37 = sin(qJ(3));
t76 = t37 * pkin(3);
t40 = cos(qJ(3));
t75 = t40 * pkin(3);
t74 = rSges(4,3) + pkin(7);
t72 = rSges(5,1) * t30;
t71 = rSges(6,1) * t27;
t39 = cos(qJ(6));
t70 = rSges(7,1) * t39;
t69 = rSges(5,2) * t29;
t67 = t26 * t38;
t66 = t26 * t41;
t65 = t27 * t38;
t64 = t30 * t38;
t63 = t38 * t36;
t62 = t38 * t39;
t61 = t41 * t36;
t60 = t41 * t39;
t59 = rSges(5,3) - t42;
t58 = t41 * pkin(1) + t38 * qJ(2);
t10 = t76 + t82;
t32 = t41 * qJ(2);
t34 = -qJ(5) + t42;
t57 = t41 * t10 + t38 * t34 + t32;
t56 = t38 * t10 + t58;
t55 = t27 * t70;
t53 = -pkin(5) - t70;
t52 = rSges(6,1) * t65 - rSges(6,2) * t67;
t50 = t37 * rSges(4,1) + t40 * rSges(4,2);
t49 = -t29 * rSges(5,1) - t30 * rSges(5,2);
t48 = t26 * rSges(6,1) + t27 * rSges(6,2);
t47 = g(1) * t32 + g(2) * t58;
t46 = -t49 + t76;
t45 = -t48 - t82;
t44 = pkin(5) * t65 + t88 * t67 + (-t54 + t55) * t38;
t43 = -t82 + t89 + (t53 + t68) * t26;
t22 = pkin(4) * t64;
t21 = t41 * t69;
t20 = rSges(5,1) * t64;
t15 = rSges(6,2) * t66;
t11 = t75 + t81;
t6 = t38 * t11;
t4 = t26 * t60 - t63;
t3 = t26 * t61 + t62;
t2 = t26 * t62 + t61;
t1 = -t26 * t63 + t60;
t5 = [-m(2) * (g(1) * (-t38 * rSges(2,1) - t41 * rSges(2,2)) + g(2) * (t41 * rSges(2,1) - t38 * rSges(2,2))) - m(3) * (g(1) * (t41 * rSges(3,3) + t32 + (rSges(3,2) - pkin(1)) * t38) + g(2) * (-t41 * rSges(3,2) + t38 * rSges(3,3) + t58)) - m(4) * ((g(1) * t50 + g(2) * t74) * t41 + (g(1) * (-pkin(1) - t74) + g(2) * t50) * t38 + t47) - m(5) * ((g(1) * t46 + g(2) * t59) * t41 + (g(1) * (-pkin(1) - t59) + g(2) * t46) * t38 + t47) - m(6) * (g(1) * t57 + g(2) * t56 + (g(1) * t48 + g(2) * (rSges(6,3) - t34)) * t41 + (g(1) * (-rSges(6,3) - pkin(1)) + g(2) * t48) * t38) - m(7) * (g(1) * (t4 * rSges(7,1) - t3 * rSges(7,2) - t38 * pkin(1) + pkin(5) * t66 + t57) + g(2) * (t2 * rSges(7,1) + t1 * rSges(7,2) + pkin(5) * t67 - t41 * t34 + t56) - t85 * t89) (-m(3) - m(4) - m(5) + t83) * t86, -m(4) * (-g(3) * t50 + t86 * (rSges(4,1) * t40 - rSges(4,2) * t37)) - m(5) * (g(1) * (t20 + (-t69 + t75) * t38) + g(2) * (t21 + (-t72 - t75) * t41) - g(3) * t46) - m(6) * (g(1) * (t52 + t6) + g(2) * (t15 + (-t11 - t71) * t41) + g(3) * (t45 - t76)) - m(7) * (g(1) * (t44 + t6) + t84 + g(3) * (t43 - t76) + (t27 * t53 - t11 - t87) * t77) -m(5) * (g(1) * (-t38 * t69 + t20) + g(2) * t21 + g(3) * t49) - m(6) * (g(1) * (t22 + t52) + g(2) * t15 + g(3) * t45) - m(7) * (g(1) * (t22 + t44) + t84 + g(3) * t43) + (m(5) * t72 - m(6) * (-t71 - t81) - m(7) * (-pkin(5) * t27 - t55 - t81 - t87)) * t77, t83 * t85, -m(7) * (g(1) * (t1 * rSges(7,1) - t2 * rSges(7,2)) + g(2) * (t3 * rSges(7,1) + t4 * rSges(7,2)) + g(3) * (-rSges(7,1) * t36 - rSges(7,2) * t39) * t27)];
taug  = t5(:);
