% Calculate Gravitation load on the joints for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RPRRPR10_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR10_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RPRRPR10_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:34:29
% EndTime: 2019-03-09 05:34:31
% DurationCPUTime: 0.98s
% Computational Cost: add. (289->157), mult. (633->211), div. (0->0), fcn. (662->8), ass. (0->66)
t31 = sin(qJ(4));
t33 = sin(qJ(1));
t35 = cos(qJ(4));
t67 = t33 * t35;
t32 = sin(qJ(3));
t37 = cos(qJ(1));
t69 = t32 * t37;
t13 = t31 * t69 + t67;
t63 = t37 * t35;
t68 = t33 * t31;
t14 = t32 * t63 - t68;
t30 = sin(qJ(6));
t34 = cos(qJ(6));
t44 = t13 * t30 + t14 * t34;
t86 = -t13 * t34 + t14 * t30;
t88 = rSges(7,1) * t86 + rSges(7,2) * t44;
t59 = rSges(6,3) + qJ(5);
t75 = -rSges(6,1) - pkin(4);
t87 = -t31 * t59 + t35 * t75 - pkin(3);
t80 = -pkin(7) - pkin(1);
t77 = g(2) * t37;
t79 = g(1) * t33;
t83 = -t77 + t79;
t82 = -m(6) - m(7);
t81 = -pkin(4) - pkin(5);
t36 = cos(qJ(3));
t78 = g(2) * t36;
t76 = g(3) * t36;
t74 = -rSges(6,2) - pkin(8);
t73 = -rSges(5,3) - pkin(8);
t72 = -pkin(9) - rSges(7,3);
t70 = t32 * t33;
t66 = t33 * t36;
t65 = t35 * t36;
t64 = t36 * t37;
t62 = pkin(3) * t66 + pkin(8) * t70;
t26 = t37 * qJ(2);
t61 = pkin(3) * t69 + t26;
t60 = pkin(1) * t37 + qJ(2) * t33;
t58 = -pkin(8) - t72;
t57 = pkin(8) * t64;
t56 = t31 * t66;
t55 = t33 * t65;
t54 = pkin(7) * t37 + t60;
t53 = g(1) * t80;
t51 = pkin(4) * t55 + qJ(5) * t56 + t62;
t50 = pkin(3) * t70 + t54;
t49 = t58 * t37;
t11 = t32 * t68 - t63;
t12 = t31 * t37 + t32 * t67;
t2 = t11 * t34 - t12 * t30;
t3 = t11 * t30 + t12 * t34;
t48 = rSges(7,1) * t2 - rSges(7,2) * t3;
t42 = t30 * t31 + t34 * t35;
t43 = t30 * t35 - t31 * t34;
t47 = (-rSges(7,1) * t43 - rSges(7,2) * t42) * t36;
t45 = rSges(4,1) * t32 + rSges(4,2) * t36;
t41 = pkin(4) * t14 + t13 * qJ(5) + t61;
t40 = -rSges(5,1) * t35 + rSges(5,2) * t31 - pkin(3);
t39 = pkin(4) * t12 + qJ(5) * t11 + t50;
t38 = -pkin(3) + (-rSges(7,1) * t30 - rSges(7,2) * t34 - qJ(5)) * t31 + (-rSges(7,1) * t34 + rSges(7,2) * t30 + t81) * t35;
t27 = t36 * pkin(8);
t18 = qJ(5) * t65;
t9 = t13 * pkin(4);
t7 = t11 * pkin(4);
t1 = [-m(2) * (g(1) * (-rSges(2,1) * t33 - rSges(2,2) * t37) + g(2) * (rSges(2,1) * t37 - rSges(2,2) * t33)) - m(3) * (g(1) * (rSges(3,3) * t37 + t26 + (rSges(3,2) - pkin(1)) * t33) + g(2) * (-rSges(3,2) * t37 + rSges(3,3) * t33 + t60)) - m(4) * (g(1) * (rSges(4,1) * t69 + rSges(4,2) * t64 + t26) + g(2) * (rSges(4,3) * t37 + t54) + (g(1) * (-rSges(4,3) + t80) + g(2) * t45) * t33) - m(5) * (g(1) * (t14 * rSges(5,1) - t13 * rSges(5,2) - rSges(5,3) * t64 - t57 + t61) + g(2) * (rSges(5,1) * t12 - rSges(5,2) * t11 + t50) + (t73 * t78 + t53) * t33) - m(6) * (g(1) * (t14 * rSges(6,1) - rSges(6,2) * t64 + t13 * rSges(6,3) + t41 - t57) + g(2) * (rSges(6,1) * t12 + rSges(6,3) * t11 + t39) + (t74 * t78 + t53) * t33) - m(7) * (g(1) * (rSges(7,1) * t44 - rSges(7,2) * t86 + t14 * pkin(5) + t33 * t80 + t41) + g(2) * (rSges(7,1) * t3 + rSges(7,2) * t2 + pkin(5) * t12 + t39) + (g(2) * t33 * t58 + g(1) * t49) * t36) (-m(3) - m(4) - m(5) + t82) * t83, -m(4) * (-g(3) * t45 + t83 * (rSges(4,1) * t36 - rSges(4,2) * t32)) - m(5) * (g(1) * (rSges(5,1) * t55 - rSges(5,2) * t56 + t62) + g(3) * (rSges(5,3) * t36 + t27) + (rSges(5,3) * t79 + g(3) * t40) * t32 + (t32 * t73 + t36 * t40) * t77) - m(6) * (g(1) * (rSges(6,1) * t55 + rSges(6,3) * t56 + t51) + g(3) * (rSges(6,2) * t36 + t27) + (rSges(6,2) * t79 + g(3) * t87) * t32 + (t74 * t32 + t36 * t87) * t77) - m(7) * (g(1) * t51 + g(3) * t27 + (g(2) * t49 + g(3) * t38 + t72 * t79) * t32 + (g(3) * t72 + (rSges(7,1) * t42 - rSges(7,2) * t43 + t35 * pkin(5)) * t79 + t38 * t77) * t36) -m(5) * (g(1) * (-rSges(5,1) * t11 - rSges(5,2) * t12) + g(2) * (rSges(5,1) * t13 + rSges(5,2) * t14)) - m(6) * (g(1) * (-rSges(6,1) * t11 + t12 * t59 - t7) + g(2) * (rSges(6,1) * t13 - t14 * t59 + t9) + g(3) * t18) - m(7) * (g(1) * (-t11 * pkin(5) + t12 * qJ(5) - t48 - t7) + g(2) * (t13 * pkin(5) - t14 * qJ(5) - t88 + t9) + g(3) * (t18 - t47)) + ((m(5) * rSges(5,2) - m(6) * rSges(6,3)) * t35 + (m(5) * rSges(5,1) - m(6) * t75 - m(7) * t81) * t31) * t76, t82 * (g(1) * t11 - g(2) * t13 + t31 * t76) -m(7) * (g(1) * t48 + g(2) * t88 + g(3) * t47)];
taug  = t1(:);
