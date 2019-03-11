% Calculate Gravitation load on the joints for
% S6RRPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2019-03-09 11:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR13_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR13_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:24:11
% EndTime: 2019-03-09 11:24:14
% DurationCPUTime: 0.94s
% Computational Cost: add. (562->161), mult. (1215->220), div. (0->0), fcn. (1429->12), ass. (0->63)
t82 = rSges(5,3) + pkin(9);
t44 = sin(qJ(2));
t45 = sin(qJ(1));
t47 = cos(qJ(2));
t48 = cos(qJ(1));
t72 = cos(pkin(6));
t68 = t48 * t72;
t21 = t44 * t68 + t45 * t47;
t69 = t45 * t72;
t23 = -t44 * t69 + t47 * t48;
t97 = g(1) * t23 + g(2) * t21;
t20 = t44 * t45 - t47 * t68;
t22 = t48 * t44 + t47 * t69;
t96 = g(1) * t22 + g(2) * t20;
t43 = sin(qJ(4));
t46 = cos(qJ(4));
t40 = sin(pkin(6));
t79 = t40 * t47;
t18 = t72 * t43 + t46 * t79;
t80 = t40 * t45;
t5 = -t22 * t46 + t43 * t80;
t78 = t40 * t48;
t60 = t20 * t46 + t43 * t78;
t95 = -g(1) * t5 + g(2) * t60 - g(3) * t18;
t19 = -t43 * t79 + t72 * t46;
t6 = t22 * t43 + t46 * t80;
t61 = -t20 * t43 + t46 * t78;
t94 = g(1) * t6 - g(2) * t61 + g(3) * t19;
t91 = -m(6) - m(7);
t39 = sin(pkin(11));
t90 = pkin(5) * t39;
t83 = g(3) * t40;
t81 = t40 * t44;
t77 = rSges(7,3) + pkin(10) + qJ(5);
t76 = pkin(2) * t79 + qJ(3) * t81;
t75 = t48 * pkin(1) + pkin(8) * t80;
t74 = rSges(4,3) + qJ(3);
t73 = qJ(5) + rSges(6,3);
t71 = t23 * pkin(2) + t75;
t70 = -t45 * pkin(1) + pkin(8) * t78;
t67 = g(3) * (pkin(9) * t79 + t76);
t66 = -t21 * pkin(2) + t70;
t65 = rSges(5,1) * t43 + rSges(5,2) * t46;
t41 = cos(pkin(11));
t64 = rSges(6,1) * t39 + rSges(6,2) * t41;
t63 = rSges(6,1) * t41 - rSges(6,2) * t39 + pkin(4);
t34 = pkin(5) * t41 + pkin(4);
t38 = pkin(11) + qJ(6);
t35 = sin(t38);
t36 = cos(t38);
t59 = rSges(7,1) * t36 - rSges(7,2) * t35 + t34;
t57 = pkin(3) * t80 + t22 * qJ(3) + t71;
t56 = rSges(7,1) * t35 + rSges(7,2) * t36 + t90;
t55 = pkin(3) * t78 - t20 * qJ(3) + t66;
t54 = -pkin(9) - t56;
t14 = t20 * pkin(2);
t16 = t22 * pkin(2);
t53 = -g(1) * t16 - g(2) * t14 + t67;
t52 = t63 * t43 - t73 * t46;
t51 = t59 * t43 - t77 * t46;
t3 = t23 * t35 + t36 * t6;
t2 = t23 * t36 - t35 * t6;
t1 = [-m(2) * (g(1) * (-t45 * rSges(2,1) - rSges(2,2) * t48) + g(2) * (rSges(2,1) * t48 - t45 * rSges(2,2))) - m(3) * (g(1) * (-t21 * rSges(3,1) + t20 * rSges(3,2) + rSges(3,3) * t78 + t70) + g(2) * (rSges(3,1) * t23 - rSges(3,2) * t22 + rSges(3,3) * t80 + t75)) - m(4) * (g(1) * (rSges(4,1) * t78 + t21 * rSges(4,2) - t74 * t20 + t66) + g(2) * (rSges(4,1) * t80 - rSges(4,2) * t23 + t74 * t22 + t71)) - m(5) * (g(1) * (rSges(5,1) * t61 - rSges(5,2) * t60 - t82 * t21 + t55) + g(2) * (rSges(5,1) * t6 - rSges(5,2) * t5 + t82 * t23 + t57)) - m(6) * (g(1) * (t61 * pkin(4) - t21 * pkin(9) + (-t21 * t39 + t41 * t61) * rSges(6,1) + (-t21 * t41 - t39 * t61) * rSges(6,2) + t73 * t60 + t55) + g(2) * (t23 * pkin(9) + t6 * pkin(4) + (t23 * t39 + t41 * t6) * rSges(6,1) + (t23 * t41 - t39 * t6) * rSges(6,2) + t73 * t5 + t57)) - m(7) * (g(1) * (t54 * t21 + t59 * t61 + t60 * t77 + t55) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t2 + t34 * t6 + t77 * t5 + (pkin(9) + t90) * t23 + t57)) -m(3) * (g(1) * (-rSges(3,1) * t22 - rSges(3,2) * t23) + g(2) * (-rSges(3,1) * t20 - rSges(3,2) * t21) + (rSges(3,1) * t47 - rSges(3,2) * t44) * t83) - m(4) * (g(1) * (rSges(4,2) * t22 + t74 * t23 - t16) + g(2) * (rSges(4,2) * t20 + t74 * t21 - t14) + g(3) * ((-rSges(4,2) * t47 + rSges(4,3) * t44) * t40 + t76)) - m(5) * (g(1) * (-t82 * t22 - t16) + g(2) * (-t82 * t20 - t14) + t67 + (rSges(5,3) * t47 + t65 * t44) * t83 + t97 * (qJ(3) + t65)) - m(6) * ((t52 * t44 + t64 * t47) * t83 + t53 + t96 * (-pkin(9) - t64) + t97 * (qJ(3) + t52)) - m(7) * ((t51 * t44 + t56 * t47) * t83 + t53 + t96 * t54 + t97 * (qJ(3) + t51)) (-m(4) - m(5) + t91) * (-g(3) * t79 + t96) -m(5) * (g(1) * (-rSges(5,1) * t5 - rSges(5,2) * t6) + g(2) * (rSges(5,1) * t60 + rSges(5,2) * t61) + g(3) * (-rSges(5,1) * t18 - rSges(5,2) * t19)) - m(6) * (t95 * t63 + t94 * t73) - m(7) * (t95 * t59 + t94 * t77) -t91 * t95, -m(7) * (g(1) * (rSges(7,1) * t2 - rSges(7,2) * t3) + g(2) * ((t21 * t36 + t35 * t61) * rSges(7,1) + (-t21 * t35 + t36 * t61) * rSges(7,2)) + g(3) * ((-t19 * t35 + t36 * t81) * rSges(7,1) + (-t19 * t36 - t35 * t81) * rSges(7,2)))];
taug  = t1(:);
