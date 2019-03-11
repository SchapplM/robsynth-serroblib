% Calculate Gravitation load on the joints for
% S6PRRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta4]';
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
% Datum: 2019-03-08 21:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRRPPR2_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRRPPR2_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:04:08
% EndTime: 2019-03-08 21:04:10
% DurationCPUTime: 0.82s
% Computational Cost: add. (494->141), mult. (857->206), div. (0->0), fcn. (984->12), ass. (0->64)
t38 = sin(pkin(10));
t43 = sin(qJ(2));
t46 = cos(qJ(2));
t69 = cos(pkin(10));
t70 = cos(pkin(6));
t52 = t70 * t69;
t21 = t38 * t46 + t43 * t52;
t42 = sin(qJ(3));
t45 = cos(qJ(3));
t39 = sin(pkin(6));
t64 = t39 * t69;
t48 = -t21 * t42 - t45 * t64;
t47 = t48 * pkin(3);
t37 = qJ(3) + pkin(11);
t35 = sin(t37);
t36 = cos(t37);
t5 = t21 * t35 + t36 * t64;
t95 = -t5 * pkin(4) + t47;
t77 = t39 * t43;
t94 = -t42 * t77 + t70 * t45;
t65 = t38 * t70;
t23 = -t43 * t65 + t46 * t69;
t76 = t39 * t45;
t93 = -t23 * t42 + t38 * t76;
t92 = pkin(4) * t36 + qJ(5) * t35;
t83 = rSges(7,3) + pkin(9);
t91 = t36 * t83;
t22 = t43 * t69 + t46 * t65;
t20 = t38 * t43 - t46 * t52;
t86 = g(2) * t20;
t90 = -g(1) * t22 - t86;
t89 = -m(6) - m(7);
t85 = g(3) * t39;
t84 = rSges(4,3) + pkin(8);
t41 = sin(qJ(6));
t80 = t35 * t41;
t44 = cos(qJ(6));
t79 = t35 * t44;
t78 = t38 * t39;
t75 = t39 * t46;
t34 = pkin(3) * t45 + pkin(2);
t40 = -qJ(4) - pkin(8);
t74 = -t20 * t34 - t21 * t40;
t73 = -t22 * t34 - t23 * t40;
t71 = rSges(6,3) + qJ(5);
t68 = -m(5) + t89;
t62 = -t20 * t92 + t74;
t61 = -t22 * t92 + t73;
t26 = t34 * t75;
t60 = g(3) * (t75 * t92 + t26);
t59 = t93 * pkin(3);
t57 = rSges(5,1) * t36 - rSges(5,2) * t35;
t56 = rSges(7,1) * t41 + rSges(7,2) * t44;
t55 = rSges(6,2) * t36 - rSges(6,3) * t35;
t7 = t23 * t35 - t36 * t78;
t54 = -t7 * pkin(4) + t59;
t53 = t94 * pkin(3);
t51 = rSges(4,1) * t45 - rSges(4,2) * t42 + pkin(2);
t16 = t35 * t77 - t36 * t70;
t49 = -t16 * pkin(4) + t53;
t17 = t35 * t70 + t36 * t77;
t8 = t23 * t36 + t35 * t78;
t6 = t21 * t36 - t35 * t64;
t1 = [(-m(2) - m(3) - m(4) + t68) * g(3), -m(3) * (g(1) * (-rSges(3,1) * t22 - rSges(3,2) * t23) + g(2) * (-rSges(3,1) * t20 - rSges(3,2) * t21) + (rSges(3,1) * t46 - rSges(3,2) * t43) * t85) - m(4) * (g(1) * (-t22 * t51 + t23 * t84) + g(2) * t84 * t21 - t51 * t86 + (t43 * t84 + t46 * t51) * t85) - m(5) * (g(1) * (rSges(5,3) * t23 - t22 * t57 + t73) + g(2) * (rSges(5,3) * t21 - t20 * t57 + t74) + g(3) * t26 + (t57 * t46 + (rSges(5,3) - t40) * t43) * t85) - m(6) * (g(1) * (rSges(6,1) * t23 + t22 * t55 + t61) + g(2) * (rSges(6,1) * t21 + t20 * t55 + t62) + t60 + (-t55 * t46 + (rSges(6,1) - t40) * t43) * t85) - m(7) * (g(1) * (t23 * pkin(5) + (-t22 * t80 + t23 * t44) * rSges(7,1) + (-t22 * t79 - t23 * t41) * rSges(7,2) + t61) + g(2) * (t21 * pkin(5) + (-t20 * t80 + t21 * t44) * rSges(7,1) + (-t20 * t79 - t21 * t41) * rSges(7,2) + t62) + t60 + t90 * t91 + ((t44 * rSges(7,1) - t41 * rSges(7,2) + pkin(5) - t40) * t43 + (t35 * t56 + t91) * t46) * t85) -m(4) * (g(1) * (t93 * rSges(4,1) + (-t23 * t45 - t42 * t78) * rSges(4,2)) + g(2) * (t48 * rSges(4,1) + (-t21 * t45 + t42 * t64) * rSges(4,2)) + g(3) * (t94 * rSges(4,1) + (-t42 * t70 - t43 * t76) * rSges(4,2))) - m(5) * (g(1) * (-rSges(5,1) * t7 - rSges(5,2) * t8 + t59) + g(2) * (-t5 * rSges(5,1) - t6 * rSges(5,2) + t47) + g(3) * (-rSges(5,1) * t16 - rSges(5,2) * t17 + t53)) - m(6) * (g(1) * (rSges(6,2) * t7 + t71 * t8 + t54) + g(2) * (t5 * rSges(6,2) + t6 * t71 + t95) + g(3) * (rSges(6,2) * t16 + t17 * t71 + t49)) + (-g(1) * (-t83 * t7 + t54) - g(2) * (-t83 * t5 + t95) - g(3) * (-t83 * t16 + t49) - (g(1) * t8 + g(2) * t6 + g(3) * t17) * (qJ(5) + t56)) * m(7), t68 * (-g(3) * t75 - t90) t89 * (g(1) * t7 + g(2) * t5 + g(3) * t16) -m(7) * (g(1) * ((-t22 * t41 + t44 * t7) * rSges(7,1) + (-t22 * t44 - t41 * t7) * rSges(7,2)) + g(2) * ((-t20 * t41 + t44 * t5) * rSges(7,1) + (-t20 * t44 - t41 * t5) * rSges(7,2)) + g(3) * ((t16 * t44 + t41 * t75) * rSges(7,1) + (-t16 * t41 + t44 * t75) * rSges(7,2)))];
taug  = t1(:);
