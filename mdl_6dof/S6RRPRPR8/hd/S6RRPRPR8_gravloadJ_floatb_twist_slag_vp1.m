% Calculate Gravitation load on the joints for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2019-03-09 10:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPR8_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR8_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPR8_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:49:34
% EndTime: 2019-03-09 10:49:38
% DurationCPUTime: 1.17s
% Computational Cost: add. (530->168), mult. (720->227), div. (0->0), fcn. (750->10), ass. (0->66)
t32 = pkin(10) + qJ(4);
t27 = sin(t32);
t28 = cos(t32);
t41 = cos(qJ(1));
t38 = sin(qJ(1));
t40 = cos(qJ(2));
t73 = t38 * t40;
t11 = t27 * t73 + t28 * t41;
t71 = t41 * t27;
t12 = t28 * t73 - t71;
t36 = sin(qJ(6));
t39 = cos(qJ(6));
t49 = t11 * t36 + t12 * t39;
t89 = t11 * t39 - t12 * t36;
t90 = t89 * rSges(7,1) - t49 * rSges(7,2);
t82 = g(2) * t38;
t87 = g(1) * t41 + t82;
t86 = -m(6) - m(7);
t85 = -pkin(4) - pkin(5);
t84 = g(1) * t38;
t37 = sin(qJ(2));
t81 = g(3) * t37;
t80 = t38 * pkin(1);
t79 = -rSges(6,1) - pkin(4);
t78 = rSges(7,3) + pkin(9);
t77 = rSges(3,2) * t37;
t33 = sin(pkin(10));
t75 = t33 * t41;
t74 = t38 * t33;
t34 = cos(pkin(10));
t26 = pkin(3) * t34 + pkin(2);
t19 = t40 * t26;
t72 = t40 * t41;
t35 = -pkin(8) - qJ(3);
t70 = rSges(6,2) - t35;
t69 = rSges(5,3) - t35;
t68 = t41 * pkin(1) + t38 * pkin(7);
t67 = rSges(4,3) + qJ(3);
t66 = rSges(6,3) + qJ(5);
t65 = -t35 - t78;
t30 = t41 * pkin(7);
t64 = t38 * t37 * t35 + pkin(3) * t75 + t30;
t63 = -pkin(1) - t19;
t61 = t70 * t41;
t60 = t69 * t41;
t59 = t41 * t67;
t58 = pkin(3) * t74 + t26 * t72 + t68;
t57 = g(3) * (t19 + (t28 * pkin(4) + t27 * qJ(5)) * t40);
t56 = t65 * t41;
t14 = t38 * t27 + t28 * t72;
t55 = t14 * pkin(4) + t58;
t13 = -t38 * t28 + t40 * t71;
t2 = t13 * t39 - t14 * t36;
t3 = t13 * t36 + t14 * t39;
t54 = rSges(7,1) * t2 - rSges(7,2) * t3;
t47 = t27 * t36 + t28 * t39;
t48 = t27 * t39 - t28 * t36;
t53 = (rSges(7,1) * t48 - rSges(7,2) * t47) * t37;
t52 = rSges(3,1) * t40 - t77;
t50 = rSges(5,1) * t28 - rSges(5,2) * t27;
t46 = rSges(4,1) * t34 - rSges(4,2) * t33 + pkin(2);
t44 = -t12 * pkin(4) - t11 * qJ(5) + t64;
t17 = t37 * t28 * qJ(5);
t9 = t13 * pkin(4);
t7 = t11 * pkin(4);
t1 = [-m(2) * (g(1) * (-t38 * rSges(2,1) - rSges(2,2) * t41) + g(2) * (rSges(2,1) * t41 - t38 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t41 + t30) + g(2) * (rSges(3,1) * t72 - t41 * t77 + t68) + (g(1) * (-pkin(1) - t52) + g(2) * rSges(3,3)) * t38) - m(4) * (g(1) * (-pkin(2) * t73 - t80 + t30 + (-t34 * t73 + t75) * rSges(4,1) + (t33 * t73 + t34 * t41) * rSges(4,2)) + g(2) * (pkin(2) * t72 + (t34 * t72 + t74) * rSges(4,1) + (-t33 * t72 + t38 * t34) * rSges(4,2) + t68) + (g(2) * t59 - t67 * t84) * t37) - m(5) * (g(1) * (-rSges(5,1) * t12 + rSges(5,2) * t11 + t64) + g(2) * (t14 * rSges(5,1) - t13 * rSges(5,2) + t37 * t60 + t58) + (-rSges(5,3) * t37 + t63) * t84) - m(6) * (g(1) * (-rSges(6,1) * t12 - rSges(6,3) * t11 + t44) + g(2) * (t14 * rSges(6,1) + t66 * t13 + t37 * t61 + t55) + (-rSges(6,2) * t37 + t63) * t84) - m(7) * (g(1) * (-t49 * rSges(7,1) - rSges(7,2) * t89 - t12 * pkin(5) - t38 * t19 + t44 - t80) + g(2) * (t3 * rSges(7,1) + t2 * rSges(7,2) + t14 * pkin(5) + t13 * qJ(5) + t55) + (g(2) * t56 + t78 * t84) * t37) -m(3) * (g(3) * t52 + t87 * (-rSges(3,1) * t37 - rSges(3,2) * t40)) - m(4) * ((g(1) * t59 + g(3) * t46 + t67 * t82) * t40 + (g(3) * t67 - t87 * t46) * t37) - m(5) * (g(3) * t19 + (g(1) * t60 + g(3) * t50 + t69 * t82) * t40 + (g(3) * t69 + t87 * (-t26 - t50)) * t37) - m(6) * (t57 + (g(3) * (rSges(6,1) * t28 + rSges(6,3) * t27) + g(1) * t61 + t70 * t82) * t40 + (g(3) * t70 + t87 * (-t66 * t27 + t79 * t28 - t26)) * t37) - m(7) * (t57 + (g(3) * (t47 * rSges(7,1) + t48 * rSges(7,2) + t28 * pkin(5)) + g(1) * t56 + t65 * t82) * t40 + (g(3) * t65 + t87 * (-t26 + (-rSges(7,1) * t36 - rSges(7,2) * t39 - qJ(5)) * t27 + (-rSges(7,1) * t39 + rSges(7,2) * t36 + t85) * t28)) * t37) (-m(4) - m(5) + t86) * (-g(3) * t40 + t87 * t37) -m(5) * (g(1) * (-rSges(5,1) * t13 - rSges(5,2) * t14) + g(2) * (-rSges(5,1) * t11 - rSges(5,2) * t12)) - m(6) * (g(1) * (-rSges(6,1) * t13 + t66 * t14 - t9) + g(2) * (-rSges(6,1) * t11 + t66 * t12 - t7) + g(3) * t17) - m(7) * (g(1) * (-t13 * pkin(5) + t14 * qJ(5) - t54 - t9) + g(2) * (-t11 * pkin(5) + t12 * qJ(5) - t7 - t90) + g(3) * (t17 - t53)) + ((m(5) * rSges(5,2) - m(6) * rSges(6,3)) * t28 + (m(5) * rSges(5,1) - m(6) * t79 - m(7) * t85) * t27) * t81, t86 * (g(1) * t13 + g(2) * t11 + t27 * t81) -m(7) * (g(1) * t54 + g(2) * t90 + g(3) * t53)];
taug  = t1(:);
