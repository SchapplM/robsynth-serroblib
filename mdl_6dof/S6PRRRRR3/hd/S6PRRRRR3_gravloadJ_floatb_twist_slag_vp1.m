% Calculate Gravitation load on the joints for
% S6PRRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,d6,theta1]';
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
% Datum: 2019-03-09 00:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRRRR3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRRRR3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:47:58
% EndTime: 2019-03-09 00:48:00
% DurationCPUTime: 0.91s
% Computational Cost: add. (728->144), mult. (1326->213), div. (0->0), fcn. (1576->14), ass. (0->72)
t43 = sin(pkin(12));
t47 = sin(qJ(2));
t50 = cos(qJ(2));
t78 = cos(pkin(12));
t79 = cos(pkin(6));
t68 = t79 * t78;
t20 = t43 * t50 + t47 * t68;
t76 = t43 * t79;
t22 = -t47 * t76 + t78 * t50;
t106 = g(1) * t22 + g(2) * t20;
t19 = t43 * t47 - t50 * t68;
t21 = t78 * t47 + t50 * t76;
t105 = g(1) * t21 + g(2) * t19;
t46 = sin(qJ(3));
t49 = cos(qJ(3));
t44 = sin(pkin(6));
t75 = t44 * t78;
t13 = -t20 * t46 - t49 * t75;
t83 = t44 * t49;
t15 = -t22 * t46 + t43 * t83;
t84 = t44 * t47;
t23 = -t46 * t84 + t79 * t49;
t104 = g(1) * t15 + g(2) * t13 + g(3) * t23;
t14 = t20 * t49 - t46 * t75;
t16 = t43 * t44 * t46 + t22 * t49;
t24 = t79 * t46 + t47 * t83;
t103 = g(1) * t16 + g(2) * t14 + g(3) * t24;
t42 = qJ(4) + qJ(5);
t38 = cos(t42);
t48 = cos(qJ(4));
t40 = t48 * pkin(4);
t27 = pkin(5) * t38 + t40;
t39 = qJ(6) + t42;
t34 = sin(t39);
t35 = cos(t39);
t62 = rSges(7,1) * t35 - rSges(7,2) * t34 + pkin(3) + t27;
t51 = -pkin(10) - pkin(9);
t80 = rSges(7,3) + pkin(11) - t51;
t102 = t80 * t46 + t62 * t49;
t37 = sin(t42);
t63 = rSges(6,1) * t38 - rSges(6,2) * t37 + pkin(3) + t40;
t81 = rSges(6,3) - t51;
t101 = t81 * t46 + t63 * t49;
t45 = sin(qJ(4));
t67 = rSges(5,1) * t48 - rSges(5,2) * t45 + pkin(3);
t86 = rSges(5,3) + pkin(9);
t100 = t86 * t46 + t67 * t49;
t99 = (-t14 * t34 + t19 * t35) * rSges(7,1) + (-t14 * t35 - t19 * t34) * rSges(7,2);
t98 = (-t16 * t34 + t21 * t35) * rSges(7,1) + (-t16 * t35 - t21 * t34) * rSges(7,2);
t89 = g(3) * t44;
t88 = t45 * pkin(4);
t87 = rSges(4,3) + pkin(8);
t82 = t44 * t50;
t85 = (-t24 * t34 - t35 * t82) * rSges(7,1) + (-t24 * t35 + t34 * t82) * rSges(7,2);
t77 = g(3) * (pkin(2) * t82 + pkin(8) * t84);
t74 = rSges(4,1) * t49 - rSges(4,2) * t46;
t73 = t45 * rSges(5,1) + t48 * rSges(5,2);
t72 = -t14 * t37 + t19 * t38;
t71 = -t14 * t45 + t19 * t48;
t70 = -t16 * t37 + t21 * t38;
t69 = -t16 * t45 + t21 * t48;
t65 = -t24 * t37 - t38 * t82;
t64 = -t24 * t45 - t48 * t82;
t26 = pkin(5) * t37 + t88;
t61 = t34 * rSges(7,1) + t35 * rSges(7,2) + t26;
t59 = t37 * rSges(6,1) + t38 * rSges(6,2) + t88;
t17 = t19 * pkin(2);
t18 = t21 * pkin(2);
t57 = -g(1) * t18 - g(2) * t17 + t77;
t53 = m(7) * (g(1) * t98 + g(2) * t99 + g(3) * t85);
t52 = m(6) * (g(1) * (t70 * rSges(6,1) + (-t16 * t38 - t21 * t37) * rSges(6,2)) + g(2) * (t72 * rSges(6,1) + (-t14 * t38 - t19 * t37) * rSges(6,2)) + g(3) * (t65 * rSges(6,1) + (-t24 * t38 + t37 * t82) * rSges(6,2)));
t1 = [(-m(2) - m(3) - m(4) - m(5) - m(6) - m(7)) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t21 - rSges(3,2) * t22) + g(2) * (-rSges(3,1) * t19 - rSges(3,2) * t20) + (rSges(3,1) * t50 - rSges(3,2) * t47) * t89) - m(4) * (g(1) * (-t74 * t21 + t87 * t22 - t18) + g(2) * (-t74 * t19 + t87 * t20 - t17) + t77 + (rSges(4,3) * t47 + t74 * t50) * t89) - m(5) * ((t100 * t50 + t73 * t47) * t89 + t57 + t106 * (pkin(8) + t73) - t105 * t100) - m(6) * ((t101 * t50 + t59 * t47) * t89 + t57 + t106 * (pkin(8) + t59) - t105 * t101) - m(7) * ((t102 * t50 + t61 * t47) * t89 + t57 + t106 * (pkin(8) + t61) - t105 * t102) -m(4) * (g(1) * (rSges(4,1) * t15 - rSges(4,2) * t16) + g(2) * (rSges(4,1) * t13 - rSges(4,2) * t14) + g(3) * (rSges(4,1) * t23 - rSges(4,2) * t24)) - m(5) * (t103 * t86 + t104 * t67) - m(6) * (t103 * t81 + t104 * t63) - m(7) * (t103 * t80 + t104 * t62) -m(5) * (g(1) * (t69 * rSges(5,1) + (-t16 * t48 - t21 * t45) * rSges(5,2)) + g(2) * (t71 * rSges(5,1) + (-t14 * t48 - t19 * t45) * rSges(5,2)) + g(3) * (t64 * rSges(5,1) + (-t24 * t48 + t45 * t82) * rSges(5,2))) - t52 - m(7) * (g(1) * (-t16 * t26 + t21 * t27 + t98) + g(2) * (-t14 * t26 + t19 * t27 + t99) + g(3) * (-t24 * t26 - t27 * t82 + t85)) - m(6) * (g(1) * t69 + g(2) * t71 + g(3) * t64) * pkin(4), -t52 - t53 - m(7) * (g(1) * t70 + g(2) * t72 + g(3) * t65) * pkin(5), -t53];
taug  = t1(:);
