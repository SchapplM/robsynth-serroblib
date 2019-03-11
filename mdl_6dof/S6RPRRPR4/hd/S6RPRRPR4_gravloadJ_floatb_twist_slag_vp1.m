% Calculate Gravitation load on the joints for
% S6RPRRPR4
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
% Datum: 2019-03-09 05:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR4_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR4_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR4_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR4_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR4_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR4_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:09:04
% EndTime: 2019-03-09 05:09:05
% DurationCPUTime: 0.74s
% Computational Cost: add. (495->115), mult. (402->152), div. (0->0), fcn. (357->12), ass. (0->60)
t40 = pkin(10) + qJ(3);
t36 = qJ(4) + t40;
t28 = sin(t36);
t29 = cos(t36);
t43 = cos(pkin(11));
t66 = -rSges(6,1) * t43 - pkin(4);
t41 = sin(pkin(11));
t84 = rSges(6,2) * t41;
t55 = (rSges(6,3) + qJ(5)) * t28 + (-t66 - t84) * t29;
t97 = qJ(5) * t29 + t28 * t84;
t94 = rSges(5,1) * t29 - rSges(5,2) * t28;
t47 = sin(qJ(1));
t48 = cos(qJ(1));
t93 = g(1) * t48 + g(2) * t47;
t30 = pkin(5) * t43 + pkin(4);
t45 = -pkin(9) - qJ(5);
t52 = t29 * t30 + (rSges(7,3) - t45) * t28;
t92 = t93 * t28;
t35 = cos(t40);
t27 = pkin(3) * t35;
t44 = cos(pkin(10));
t31 = pkin(2) * t44 + pkin(1);
t13 = t27 + t31;
t6 = t48 * t13;
t91 = g(2) * t6;
t90 = -m(6) - m(7);
t33 = sin(t40);
t89 = pkin(3) * t33;
t39 = pkin(11) + qJ(6);
t34 = cos(t39);
t85 = rSges(7,1) * t34;
t32 = sin(t39);
t83 = rSges(7,2) * t32;
t80 = t47 * t29;
t79 = t47 * t32;
t78 = t47 * t34;
t77 = t48 * t29;
t76 = t48 * t32;
t75 = t48 * t34;
t46 = -pkin(7) - qJ(2);
t74 = rSges(4,3) - t46;
t38 = -pkin(8) + t46;
t73 = rSges(5,3) - t38;
t71 = rSges(3,3) + qJ(2);
t70 = rSges(6,3) * t80 + t47 * t97;
t68 = t28 * t83;
t67 = rSges(6,3) * t77 + t48 * t97;
t65 = pkin(5) * t41 - t38;
t64 = -t30 - t85;
t62 = rSges(4,1) * t35 - rSges(4,2) * t33;
t59 = rSges(3,1) * t44 - rSges(3,2) * sin(pkin(10)) + pkin(1);
t57 = t31 + t62;
t56 = rSges(6,1) * t41 + rSges(6,2) * t43 - t38;
t53 = g(1) * (rSges(7,3) * t77 + t48 * t68) + g(2) * (rSges(7,3) * t80 + t47 * t68);
t51 = t52 + (-t83 + t85) * t29;
t5 = t29 * t75 + t79;
t4 = -t29 * t76 + t78;
t3 = -t29 * t78 + t76;
t2 = t29 * t79 + t75;
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t47 - rSges(2,2) * t48) + g(2) * (rSges(2,1) * t48 - rSges(2,2) * t47)) - m(3) * ((g(1) * t71 + g(2) * t59) * t48 + (-g(1) * t59 + g(2) * t71) * t47) - m(4) * ((g(1) * t74 + g(2) * t57) * t48 + (-g(1) * t57 + g(2) * t74) * t47) - m(5) * (t91 + (g(1) * t73 + g(2) * t94) * t48 + (g(1) * (-t13 - t94) + g(2) * t73) * t47) - m(6) * (t91 + (g(1) * t56 + g(2) * t55) * t48 + (g(2) * t56 + (-t13 - t55) * g(1)) * t47) - m(7) * (g(1) * (t3 * rSges(7,1) + t2 * rSges(7,2)) + g(2) * (t5 * rSges(7,1) + t4 * rSges(7,2) + t6) + (g(1) * t65 + g(2) * t52) * t48 + (g(1) * (-t13 - t52) + g(2) * t65) * t47) (-m(3) - m(4) - m(5) + t90) * (g(1) * t47 - g(2) * t48) -m(4) * (g(3) * t62 + t93 * (-rSges(4,1) * t33 - rSges(4,2) * t35)) - m(5) * (g(3) * (t27 + t94) + t93 * (-rSges(5,1) * t28 - rSges(5,2) * t29 - t89)) - m(6) * (g(1) * (-t48 * t89 + t67) + g(2) * (-t47 * t89 + t70) + g(3) * (t27 + t55) + t66 * t92) - m(7) * (g(3) * (t27 + t51) + t53 + t93 * (t28 * t64 - t29 * t45 - t89)) -m(5) * g(3) * t94 - m(6) * (g(1) * t67 + g(2) * t70 + g(3) * t55) - m(7) * (g(3) * t51 + t53) + t93 * ((m(5) * rSges(5,2) + m(7) * t45) * t29 + (m(5) * rSges(5,1) - m(6) * t66 - m(7) * t64) * t28) t90 * (-g(3) * t29 + t92) -m(7) * (g(1) * (rSges(7,1) * t4 - rSges(7,2) * t5) + g(2) * (-rSges(7,1) * t2 + rSges(7,2) * t3) + g(3) * (-rSges(7,1) * t32 - rSges(7,2) * t34) * t28)];
taug  = t1(:);
