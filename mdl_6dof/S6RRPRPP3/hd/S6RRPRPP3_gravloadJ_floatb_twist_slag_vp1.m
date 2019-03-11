% Calculate Gravitation load on the joints for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPRPP3_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP3_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPRPP3_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:54:15
% EndTime: 2019-03-09 09:54:17
% DurationCPUTime: 0.86s
% Computational Cost: add. (452->151), mult. (613->199), div. (0->0), fcn. (610->8), ass. (0->58)
t31 = sin(qJ(2));
t71 = g(3) * t31;
t27 = pkin(9) + qJ(4);
t23 = cos(t27);
t78 = t23 * t71;
t34 = cos(qJ(1));
t32 = sin(qJ(1));
t73 = g(2) * t32;
t77 = g(1) * t34 + t73;
t57 = rSges(7,2) + qJ(5);
t54 = rSges(7,3) + qJ(6);
t76 = -m(6) - m(7);
t75 = g(1) * t32;
t72 = qJ(5) * t78;
t70 = t32 * pkin(1);
t69 = -rSges(7,1) - pkin(5);
t68 = rSges(6,2) - pkin(4);
t67 = rSges(3,2) * t31;
t28 = sin(pkin(9));
t65 = t28 * t34;
t64 = t32 * t28;
t33 = cos(qJ(2));
t63 = t32 * t33;
t29 = cos(pkin(9));
t21 = pkin(3) * t29 + pkin(2);
t14 = t33 * t21;
t62 = t33 * t34;
t22 = sin(t27);
t61 = t34 * t22;
t30 = -pkin(8) - qJ(3);
t60 = rSges(6,1) - t30;
t59 = rSges(5,3) - t30;
t58 = t34 * pkin(1) + t32 * pkin(7);
t56 = rSges(4,3) + qJ(3);
t55 = rSges(6,3) + qJ(5);
t53 = -t30 - t69;
t52 = -pkin(4) - t54;
t25 = t34 * pkin(7);
t51 = t32 * t31 * t30 + pkin(3) * t65 + t25;
t50 = -pkin(1) - t14;
t49 = t60 * t34;
t48 = t59 * t34;
t47 = t34 * t56;
t46 = pkin(3) * t64 + t21 * t62 + t58;
t45 = g(3) * (t14 + (pkin(4) * t23 + qJ(5) * t22) * t33);
t44 = t53 * t34;
t9 = t32 * t22 + t23 * t62;
t43 = t9 * pkin(4) + t46;
t42 = rSges(3,1) * t33 - t67;
t40 = rSges(5,1) * t23 - rSges(5,2) * t22;
t39 = rSges(4,1) * t29 - rSges(4,2) * t28 + pkin(2);
t6 = t22 * t63 + t23 * t34;
t7 = t23 * t63 - t61;
t37 = -t7 * pkin(4) - qJ(5) * t6 + t51;
t8 = -t32 * t23 + t33 * t61;
t4 = t8 * pkin(4);
t2 = t6 * pkin(4);
t1 = [-m(2) * (g(1) * (-t32 * rSges(2,1) - rSges(2,2) * t34) + g(2) * (rSges(2,1) * t34 - t32 * rSges(2,2))) - m(3) * (g(1) * (rSges(3,3) * t34 + t25) + g(2) * (rSges(3,1) * t62 - t34 * t67 + t58) + (g(1) * (-pkin(1) - t42) + g(2) * rSges(3,3)) * t32) - m(4) * (g(1) * (-pkin(2) * t63 - t70 + t25 + (-t29 * t63 + t65) * rSges(4,1) + (t28 * t63 + t29 * t34) * rSges(4,2)) + g(2) * (pkin(2) * t62 + (t29 * t62 + t64) * rSges(4,1) + (-t28 * t62 + t32 * t29) * rSges(4,2) + t58) + (g(2) * t47 - t56 * t75) * t31) - m(5) * (g(1) * (-rSges(5,1) * t7 + rSges(5,2) * t6 + t51) + g(2) * (t9 * rSges(5,1) - t8 * rSges(5,2) + t31 * t48 + t46) + (-rSges(5,3) * t31 + t50) * t75) - m(6) * (g(1) * (rSges(6,2) * t7 - rSges(6,3) * t6 + t37) + g(2) * (-t9 * rSges(6,2) + t31 * t49 + t55 * t8 + t43) + (-rSges(6,1) * t31 + t50) * t75) - m(7) * (g(1) * (-rSges(7,2) * t6 - t32 * t14 - t54 * t7 + t37 - t70) + g(2) * (t54 * t9 + t57 * t8 + t43) + (g(2) * t44 + t69 * t75) * t31) -m(3) * (g(3) * t42 + t77 * (-rSges(3,1) * t31 - rSges(3,2) * t33)) - m(4) * ((g(1) * t47 + g(3) * t39 + t56 * t73) * t33 + (g(3) * t56 - t77 * t39) * t31) - m(5) * (g(3) * t14 + (g(1) * t48 + g(3) * t40 + t59 * t73) * t33 + (g(3) * t59 + t77 * (-t21 - t40)) * t31) - m(6) * (t45 + (g(3) * (-rSges(6,2) * t23 + rSges(6,3) * t22) + g(1) * t49 + t60 * t73) * t33 + (g(3) * t60 + t77 * (-t55 * t22 + t68 * t23 - t21)) * t31) - m(7) * (t45 + (g(3) * (rSges(7,2) * t22 + t54 * t23) + g(1) * t44 + t53 * t73) * t33 + (g(3) * t53 + t77 * (-t57 * t22 + t52 * t23 - t21)) * t31) (-m(4) - m(5) + t76) * (-g(3) * t33 + t77 * t31) -m(5) * (g(1) * (-rSges(5,1) * t8 - rSges(5,2) * t9) + g(2) * (-rSges(5,1) * t6 - rSges(5,2) * t7)) - m(6) * (g(1) * (rSges(6,2) * t8 + t55 * t9 - t4) + g(2) * (rSges(6,2) * t6 + t55 * t7 - t2) + t72) - m(7) * (g(1) * (-t54 * t8 + t57 * t9 - t4) + g(2) * (-t54 * t6 + t57 * t7 - t2) + t72) + ((m(5) * rSges(5,2) - m(6) * rSges(6,3) - m(7) * rSges(7,2)) * t23 + (m(5) * rSges(5,1) - m(6) * t68 - m(7) * t52) * t22) * t71, t76 * (g(1) * t8 + g(2) * t6 + t22 * t71) -m(7) * (g(1) * t9 + g(2) * t7 + t78)];
taug  = t1(:);
