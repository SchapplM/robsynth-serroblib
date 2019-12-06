% Calculate Gravitation load on the joints for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% rSges [6x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (starting with base)
%   columns: x-, y-, z-coordinates
% 
% Output:
% taug [5x1]
%   joint torques required to compensate gravitation load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5PRRRP7_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP7_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5PRRRP7_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:54:07
% EndTime: 2019-12-05 16:54:09
% DurationCPUTime: 0.56s
% Computational Cost: add. (312->111), mult. (763->178), div. (0->0), fcn. (897->10), ass. (0->62)
t37 = sin(qJ(2));
t40 = cos(qJ(2));
t53 = cos(pkin(9));
t54 = cos(pkin(5));
t45 = t54 * t53;
t52 = sin(pkin(9));
t19 = t52 * t37 - t40 * t45;
t44 = t54 * t52;
t21 = t53 * t37 + t40 * t44;
t78 = g(1) * t21 + g(2) * t19;
t20 = t37 * t45 + t52 * t40;
t36 = sin(qJ(3));
t39 = cos(qJ(3));
t33 = sin(pkin(5));
t48 = t33 * t53;
t10 = t20 * t39 - t36 * t48;
t22 = -t37 * t44 + t53 * t40;
t47 = t33 * t52;
t12 = t22 * t39 + t36 * t47;
t77 = g(1) * t12 + g(2) * t10;
t38 = cos(qJ(4));
t32 = t38 * pkin(4) + pkin(3);
t56 = rSges(6,3) + qJ(5) + pkin(8);
t76 = t32 * t39 + t56 * t36;
t11 = t22 * t36 - t39 * t47;
t62 = t33 * t37;
t23 = t36 * t62 - t54 * t39;
t9 = t20 * t36 + t39 * t48;
t75 = g(1) * t11 + g(2) * t9 + g(3) * t23;
t74 = pkin(3) * t39;
t68 = g(3) * t33;
t67 = rSges(4,3) + pkin(7);
t66 = rSges(5,3) + pkin(8);
t35 = sin(qJ(4));
t65 = t20 * t35;
t64 = t22 * t35;
t61 = t33 * t40;
t60 = t35 * t37;
t59 = t35 * t39;
t58 = t38 * t39;
t57 = t39 * t40;
t55 = pkin(2) * t61 + pkin(7) * t62;
t51 = g(3) * t66;
t17 = t19 * pkin(2);
t50 = t20 * pkin(7) - t17;
t18 = t21 * pkin(2);
t49 = t22 * pkin(7) - t18;
t46 = rSges(4,1) * t39 - rSges(4,2) * t36;
t1 = -t10 * t35 + t19 * t38;
t3 = -t12 * t35 + t21 * t38;
t24 = t54 * t36 + t39 * t62;
t13 = -t24 * t35 - t38 * t61;
t16 = (t38 * t57 + t60) * t33;
t15 = (-t35 * t57 + t37 * t38) * t33;
t14 = -t24 * t38 + t35 * t61;
t8 = -t21 * t58 + t64;
t7 = t21 * t59 + t22 * t38;
t6 = -t19 * t58 + t65;
t5 = t19 * t59 + t20 * t38;
t4 = -t12 * t38 - t21 * t35;
t2 = -t10 * t38 - t19 * t35;
t25 = [(-m(2) - m(3) - m(4) - m(5) - m(6)) * g(3), -m(3) * (g(1) * (-t21 * rSges(3,1) - t22 * rSges(3,2)) + g(2) * (-t19 * rSges(3,1) - t20 * rSges(3,2)) + (rSges(3,1) * t40 - rSges(3,2) * t37) * t68) - m(4) * (g(1) * (-t46 * t21 + t67 * t22 - t18) + g(2) * (-t46 * t19 + t67 * t20 - t17) + g(3) * t55 + (rSges(4,3) * t37 + t46 * t40) * t68) - m(5) * (g(1) * (t8 * rSges(5,1) + t7 * rSges(5,2) - t21 * t74 + t49) + g(2) * (t6 * rSges(5,1) + t5 * rSges(5,2) - t19 * t74 + t50) + g(3) * (t33 * pkin(3) * t57 + t16 * rSges(5,1) + t15 * rSges(5,2) + t55) + (t51 * t61 - t78 * t66) * t36) - m(6) * (g(1) * (t8 * rSges(6,1) + t7 * rSges(6,2) + pkin(4) * t64 + t49) + g(2) * (t6 * rSges(6,1) + t5 * rSges(6,2) + pkin(4) * t65 + t50) + g(3) * (t16 * rSges(6,1) + t15 * rSges(6,2) + t55) + (pkin(4) * t60 + t76 * t40) * t68 - t78 * t76), -m(4) * (g(1) * (-t11 * rSges(4,1) - t12 * rSges(4,2)) + g(2) * (-t9 * rSges(4,1) - t10 * rSges(4,2)) + g(3) * (-t23 * rSges(4,1) - t24 * rSges(4,2))) - m(5) * (t24 * t51 + t77 * t66 + t75 * (-rSges(5,1) * t38 + rSges(5,2) * t35 - pkin(3))) - m(6) * ((g(3) * t24 + t77) * t56 + t75 * (-rSges(6,1) * t38 + rSges(6,2) * t35 - t32)), -m(5) * (g(1) * (t3 * rSges(5,1) + t4 * rSges(5,2)) + g(2) * (t1 * rSges(5,1) + t2 * rSges(5,2)) + g(3) * (t13 * rSges(5,1) + t14 * rSges(5,2))) + (-g(1) * (t3 * rSges(6,1) + t4 * rSges(6,2)) - g(2) * (t1 * rSges(6,1) + t2 * rSges(6,2)) - g(3) * (t13 * rSges(6,1) + t14 * rSges(6,2)) - (g(1) * t3 + g(2) * t1 + g(3) * t13) * pkin(4)) * m(6), -m(6) * t75];
taug = t25(:);
