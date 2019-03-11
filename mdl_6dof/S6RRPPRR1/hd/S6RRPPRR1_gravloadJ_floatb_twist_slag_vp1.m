% Calculate Gravitation load on the joints for
% S6RRPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
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
% Datum: 2019-03-09 08:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6RRPPRR1_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR1_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR1_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR1_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPRR1_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6RRPPRR1_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:45:53
% EndTime: 2019-03-09 08:45:54
% DurationCPUTime: 0.82s
% Computational Cost: add. (471->125), mult. (578->164), div. (0->0), fcn. (591->10), ass. (0->64)
t31 = sin(qJ(6));
t34 = cos(qJ(6));
t33 = sin(qJ(1));
t29 = qJ(2) + pkin(10);
t26 = sin(t29);
t27 = cos(t29);
t69 = sin(qJ(5));
t70 = cos(qJ(5));
t9 = t26 * t70 - t27 * t69;
t4 = t9 * t33;
t36 = cos(qJ(1));
t57 = t36 * t69;
t58 = t36 * t70;
t7 = -t26 * t58 + t27 * t57;
t8 = t26 * t69 + t27 * t70;
t86 = (g(1) * t7 - g(2) * t4 + g(3) * t8) * (rSges(7,1) * t34 - rSges(7,2) * t31 + pkin(5));
t85 = g(3) * t9;
t84 = rSges(7,3) + pkin(9);
t30 = -qJ(3) - pkin(7);
t71 = -pkin(8) - t30;
t22 = t26 * qJ(4);
t83 = -t27 * pkin(3) - t22;
t75 = g(1) * t36;
t82 = g(2) * t33 + t75;
t81 = t82 * t26;
t32 = sin(qJ(2));
t76 = pkin(2) * t32;
t23 = t27 * pkin(4);
t73 = rSges(3,3) + pkin(7);
t68 = t27 * t36;
t67 = rSges(5,2) - t30;
t66 = rSges(4,3) - t30;
t65 = qJ(4) * t27;
t64 = -m(5) - m(6) - m(7);
t63 = -rSges(6,3) + t71;
t35 = cos(qJ(2));
t28 = t35 * pkin(2);
t25 = t28 + pkin(1);
t20 = t36 * t25;
t62 = pkin(3) * t68 + t36 * t22 + t20;
t61 = t28 - t83;
t56 = pkin(4) * t68 + t62;
t55 = t23 + t61;
t14 = t33 * t65;
t54 = -t33 * t76 + t14;
t16 = t36 * t65;
t53 = -t36 * t76 + t16;
t5 = t8 * t33;
t52 = rSges(6,1) * t4 - rSges(6,2) * t5;
t6 = -t26 * t57 - t27 * t58;
t51 = -t7 * rSges(6,1) + t6 * rSges(6,2);
t50 = -rSges(6,1) * t8 - rSges(6,2) * t9;
t49 = t5 * t31 - t34 * t36;
t48 = -t31 * t36 - t5 * t34;
t47 = rSges(3,1) * t35 - rSges(3,2) * t32;
t45 = rSges(4,1) * t27 - rSges(4,2) * t26;
t44 = rSges(5,1) * t27 + rSges(5,3) * t26;
t43 = pkin(1) + t47;
t41 = -t25 + t83;
t38 = g(1) * (t41 - t23);
t37 = (-pkin(3) - pkin(4)) * t81;
t2 = -t31 * t33 - t34 * t6;
t1 = t31 * t6 - t33 * t34;
t3 = [-m(2) * (g(1) * (-t33 * rSges(2,1) - rSges(2,2) * t36) + g(2) * (rSges(2,1) * t36 - t33 * rSges(2,2))) - m(3) * ((g(1) * t73 + g(2) * t43) * t36 + (-g(1) * t43 + g(2) * t73) * t33) - m(4) * (g(2) * t20 + (g(1) * t66 + g(2) * t45) * t36 + (g(1) * (-t25 - t45) + g(2) * t66) * t33) - m(5) * (g(2) * t62 + (g(1) * t67 + g(2) * t44) * t36 + (g(1) * (t41 - t44) + g(2) * t67) * t33) - m(6) * (g(1) * (-t5 * rSges(6,1) - t4 * rSges(6,2)) + g(2) * (-rSges(6,1) * t6 - rSges(6,2) * t7 + t56) + t63 * t75 + (g(2) * t63 + t38) * t33) - m(7) * (g(1) * (t48 * rSges(7,1) + t49 * rSges(7,2) - t5 * pkin(5) + t71 * t36 + t84 * t4) + g(2) * (rSges(7,1) * t2 + rSges(7,2) * t1 - pkin(5) * t6 + t84 * t7 + t56) + (g(2) * t71 + t38) * t33) -m(3) * (g(3) * t47 + t82 * (-rSges(3,1) * t32 - rSges(3,2) * t35)) - m(4) * (g(3) * (t28 + t45) + t82 * (-rSges(4,1) * t26 - rSges(4,2) * t27 - t76)) - m(5) * (g(1) * t16 + g(2) * t14 + g(3) * (t44 + t61) + t82 * (rSges(5,3) * t27 - t76 + (-rSges(5,1) - pkin(3)) * t26)) - m(6) * (g(1) * (-t51 + t53) + g(2) * (-t52 + t54) + g(3) * (-t50 + t55) + t37) - m(7) * (g(1) * (t84 * t6 + t53) + g(2) * (-t84 * t5 + t54) + g(3) * (-t84 * t9 + t55) + t37 + t86) (-m(4) + t64) * (g(1) * t33 - g(2) * t36) t64 * (-g(3) * t27 + t81) -m(6) * (g(1) * t51 + g(2) * t52 + g(3) * t50) - m(7) * ((-g(1) * t6 + g(2) * t5 + t85) * t84 - t86) -m(7) * (g(1) * (rSges(7,1) * t1 - rSges(7,2) * t2) + g(2) * (-t49 * rSges(7,1) + t48 * rSges(7,2)) + (-t31 * rSges(7,1) - t34 * rSges(7,2)) * t85)];
taug  = t3(:);
