% Calculate Gravitation load on the joints for
% S5RRPRR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRPRR16_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR16_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR16_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR16_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR16_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRPRR16_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:44:44
% EndTime: 2019-12-31 20:44:46
% DurationCPUTime: 0.68s
% Computational Cost: add. (333->127), mult. (783->184), div. (0->0), fcn. (910->10), ass. (0->51)
t71 = rSges(5,3) + pkin(8);
t37 = sin(qJ(2));
t38 = sin(qJ(1));
t41 = cos(qJ(2));
t42 = cos(qJ(1));
t60 = cos(pkin(5));
t56 = t42 * t60;
t20 = t37 * t56 + t38 * t41;
t57 = t38 * t60;
t22 = -t37 * t57 + t41 * t42;
t80 = g(1) * t22 + g(2) * t20;
t19 = t37 * t38 - t41 * t56;
t21 = t42 * t37 + t41 * t57;
t79 = g(1) * t21 + g(2) * t19;
t35 = sin(qJ(5));
t39 = cos(qJ(5));
t36 = sin(qJ(4));
t40 = cos(qJ(4));
t34 = sin(pkin(5));
t64 = t34 * t42;
t49 = -t19 * t36 + t40 * t64;
t78 = t20 * t39 + t35 * t49;
t77 = -t20 * t35 + t39 * t49;
t72 = g(3) * t34;
t70 = rSges(6,3) + pkin(9);
t67 = t34 * t37;
t66 = t34 * t38;
t65 = t34 * t41;
t63 = pkin(2) * t65 + qJ(3) * t67;
t62 = t42 * pkin(1) + pkin(7) * t66;
t61 = rSges(4,3) + qJ(3);
t59 = t22 * pkin(2) + t62;
t58 = -t38 * pkin(1) + pkin(7) * t64;
t55 = g(3) * (pkin(8) * t65 + t63);
t54 = -t20 * pkin(2) + t58;
t53 = rSges(5,1) * t36 + rSges(5,2) * t40;
t52 = t35 * rSges(6,1) + t39 * rSges(6,2);
t51 = t39 * rSges(6,1) - t35 * rSges(6,2) + pkin(4);
t48 = t19 * t40 + t36 * t64;
t46 = pkin(3) * t66 + qJ(3) * t21 + t59;
t45 = pkin(3) * t64 - t19 * qJ(3) + t54;
t44 = t51 * t36 - t70 * t40;
t18 = -t36 * t65 + t60 * t40;
t17 = -t60 * t36 - t40 * t65;
t15 = t21 * pkin(2);
t13 = t19 * pkin(2);
t5 = t21 * t36 + t40 * t66;
t4 = -t21 * t40 + t36 * t66;
t2 = t22 * t35 + t39 * t5;
t1 = t22 * t39 - t35 * t5;
t3 = [-m(2) * (g(1) * (-t38 * rSges(2,1) - rSges(2,2) * t42) + g(2) * (rSges(2,1) * t42 - t38 * rSges(2,2))) - m(3) * (g(1) * (-t20 * rSges(3,1) + t19 * rSges(3,2) + rSges(3,3) * t64 + t58) + g(2) * (rSges(3,1) * t22 - rSges(3,2) * t21 + rSges(3,3) * t66 + t62)) - m(4) * (g(1) * (rSges(4,1) * t64 + t20 * rSges(4,2) - t61 * t19 + t54) + g(2) * (rSges(4,1) * t66 - rSges(4,2) * t22 + t61 * t21 + t59)) - m(5) * (g(1) * (rSges(5,1) * t49 - rSges(5,2) * t48 - t71 * t20 + t45) + g(2) * (rSges(5,1) * t5 - rSges(5,2) * t4 + t71 * t22 + t46)) - m(6) * (g(1) * (rSges(6,1) * t77 - rSges(6,2) * t78 + t49 * pkin(4) - t20 * pkin(8) + t70 * t48 + t45) + g(2) * (rSges(6,1) * t2 + rSges(6,2) * t1 + pkin(4) * t5 + pkin(8) * t22 + t70 * t4 + t46)), -m(3) * (g(1) * (-rSges(3,1) * t21 - rSges(3,2) * t22) + g(2) * (-rSges(3,1) * t19 - rSges(3,2) * t20) + (rSges(3,1) * t41 - rSges(3,2) * t37) * t72) - m(4) * (g(1) * (rSges(4,2) * t21 + t61 * t22 - t15) + g(2) * (rSges(4,2) * t19 + t61 * t20 - t13) + g(3) * ((-rSges(4,2) * t41 + rSges(4,3) * t37) * t34 + t63)) - m(5) * (g(1) * (-t21 * t71 - t15) + g(2) * (-t19 * t71 - t13) + t55 + (rSges(5,3) * t41 + t53 * t37) * t72 + t80 * (qJ(3) + t53)) - m(6) * (-g(1) * t15 - g(2) * t13 + t55 + (t44 * t37 + t52 * t41) * t72 + t79 * (-pkin(8) - t52) + t80 * (qJ(3) + t44)), (-m(4) - m(5) - m(6)) * (-g(3) * t65 + t79), -m(5) * (g(1) * (-rSges(5,1) * t4 - rSges(5,2) * t5) + g(2) * (rSges(5,1) * t48 + rSges(5,2) * t49) + g(3) * (rSges(5,1) * t17 - rSges(5,2) * t18)) - m(6) * (g(2) * (t48 * t51 - t49 * t70) + (t51 * t17 + t70 * t18) * g(3) + (-t51 * t4 + t70 * t5) * g(1)), -m(6) * (g(1) * (rSges(6,1) * t1 - rSges(6,2) * t2) + g(2) * (rSges(6,1) * t78 + rSges(6,2) * t77) + g(3) * ((-t18 * t35 + t39 * t67) * rSges(6,1) + (-t18 * t39 - t35 * t67) * rSges(6,2)))];
taug = t3(:);
