% Calculate Gravitation load on the joints for
% S6PRPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-03-08 20:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRRP6_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRRP6_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:18:14
% EndTime: 2019-03-08 20:18:16
% DurationCPUTime: 0.74s
% Computational Cost: add. (404->125), mult. (998->192), div. (0->0), fcn. (1185->10), ass. (0->68)
t52 = cos(qJ(4));
t46 = sin(pkin(6));
t50 = sin(qJ(2));
t78 = t46 * t50;
t45 = sin(pkin(10));
t47 = cos(pkin(10));
t53 = cos(qJ(2));
t69 = cos(pkin(6));
t62 = t50 * t69;
t31 = t45 * t53 + t47 * t62;
t33 = -t45 * t62 + t47 * t53;
t88 = g(1) * t33 + g(2) * t31;
t91 = t52 * (g(3) * t78 + t88);
t49 = sin(qJ(4));
t86 = pkin(4) * t49;
t90 = qJ(3) + t86;
t82 = rSges(7,1) + pkin(5);
t89 = rSges(7,2) + pkin(9);
t70 = rSges(7,3) + qJ(6);
t83 = g(3) * t46;
t80 = rSges(6,3) + pkin(9);
t79 = t46 * t49;
t77 = t46 * t52;
t76 = t46 * t53;
t48 = sin(qJ(5));
t75 = t48 * t49;
t51 = cos(qJ(5));
t74 = t49 * t51;
t73 = t50 * t51;
t72 = pkin(2) * t76 + qJ(3) * t78;
t71 = rSges(4,3) + qJ(3);
t67 = t48 * t78;
t66 = pkin(8) * t76 + t72;
t65 = -m(4) - m(5) - m(6) - m(7);
t61 = t53 * t69;
t30 = t45 * t50 - t47 * t61;
t27 = t30 * pkin(2);
t64 = -t30 * pkin(8) - t27;
t32 = t45 * t61 + t47 * t50;
t28 = t32 * pkin(2);
t63 = -t32 * pkin(8) - t28;
t60 = t78 * t86 + t66;
t59 = rSges(5,1) * t49 + rSges(5,2) * t52;
t58 = rSges(6,1) * t51 - rSges(6,2) * t48;
t56 = t90 * t31 + t64;
t55 = t90 * t33 + t63;
t35 = -t49 * t76 + t69 * t52;
t34 = -t69 * t49 - t52 * t76;
t29 = t34 * pkin(4);
t19 = (t48 * t53 + t49 * t73) * t46;
t18 = t49 * t67 - t51 * t76;
t17 = t35 * t51 + t67;
t16 = t35 * t48 - t46 * t73;
t15 = -t30 * t49 + t47 * t77;
t14 = t30 * t52 + t47 * t79;
t13 = t32 * t49 + t45 * t77;
t12 = t32 * t52 - t45 * t79;
t11 = t14 * pkin(4);
t10 = t12 * pkin(4);
t9 = -t32 * t48 + t33 * t74;
t8 = t32 * t51 + t33 * t75;
t7 = -t30 * t48 + t31 * t74;
t6 = t30 * t51 + t31 * t75;
t4 = -t15 * t51 + t31 * t48;
t3 = -t15 * t48 - t31 * t51;
t2 = t13 * t51 + t33 * t48;
t1 = t13 * t48 - t33 * t51;
t5 = [(-m(2) - m(3) + t65) * g(3), -m(3) * (g(1) * (-t32 * rSges(3,1) - t33 * rSges(3,2)) + g(2) * (-t30 * rSges(3,1) - t31 * rSges(3,2)) + (rSges(3,1) * t53 - rSges(3,2) * t50) * t83) - m(4) * (g(1) * (t32 * rSges(4,2) + t71 * t33 - t28) + g(2) * (t30 * rSges(4,2) + t71 * t31 - t27) + g(3) * ((-rSges(4,2) * t53 + rSges(4,3) * t50) * t46 + t72)) - m(5) * (g(1) * (-t32 * rSges(5,3) + t63) + g(2) * (-t30 * rSges(5,3) + t64) + g(3) * t66 + (rSges(5,3) * t53 + t59 * t50) * t83 + t88 * (qJ(3) + t59)) - m(6) * (g(1) * (t9 * rSges(6,1) - t8 * rSges(6,2) + t55) + g(2) * (t7 * rSges(6,1) - t6 * rSges(6,2) + t56) + g(3) * (t19 * rSges(6,1) - t18 * rSges(6,2) + t60) - t80 * t91) - m(7) * (g(1) * (t70 * t8 + t82 * t9 + t55) + g(2) * (t70 * t6 + t82 * t7 + t56) + g(3) * (t70 * t18 + t82 * t19 + t60) - t89 * t91) t65 * (g(1) * t32 + g(2) * t30 - g(3) * t76) -m(5) * (g(1) * (rSges(5,1) * t12 - rSges(5,2) * t13) + g(2) * (rSges(5,1) * t14 + rSges(5,2) * t15) + g(3) * (rSges(5,1) * t34 - rSges(5,2) * t35)) - m(6) * (g(1) * (t58 * t12 + t80 * t13 + t10) + g(2) * (t58 * t14 - t80 * t15 + t11) + g(3) * (t58 * t34 + t80 * t35 + t29)) + (-g(1) * (t89 * t13 + t10) - g(2) * (-t89 * t15 + t11) - g(3) * (t89 * t35 + t29) - (g(1) * t12 + g(2) * t14 + g(3) * t34) * (t70 * t48 + t82 * t51)) * m(7), -m(6) * (g(1) * (-rSges(6,1) * t1 - rSges(6,2) * t2) + g(2) * (-rSges(6,1) * t3 - rSges(6,2) * t4) + g(3) * (-rSges(6,1) * t16 - rSges(6,2) * t17)) - m(7) * (g(1) * (-t82 * t1 + t70 * t2) + g(2) * (-t82 * t3 + t70 * t4) + g(3) * (-t82 * t16 + t70 * t17)) -m(7) * (g(1) * t1 + g(2) * t3 + g(3) * t16)];
taug  = t5(:);
