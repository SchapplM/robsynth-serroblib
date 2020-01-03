% Calculate Gravitation load on the joints for
% S5RRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-31 21:42
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RRRPR12_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR12_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR12_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR12_gravloadJ_floatb_twist_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR12_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RRRPR12_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:37:05
% EndTime: 2019-12-31 21:37:08
% DurationCPUTime: 0.87s
% Computational Cost: add. (468->145), mult. (1016->217), div. (0->0), fcn. (1205->12), ass. (0->60)
t63 = pkin(9) + qJ(4) + rSges(6,3);
t62 = qJ(4) + rSges(5,3);
t38 = sin(qJ(3));
t39 = sin(qJ(2));
t61 = cos(pkin(5));
t35 = sin(pkin(5));
t41 = cos(qJ(3));
t66 = t35 * t41;
t15 = t38 * t61 + t39 * t66;
t40 = sin(qJ(1));
t53 = t39 * t61;
t69 = cos(qJ(2));
t70 = cos(qJ(1));
t17 = t40 * t69 + t53 * t70;
t58 = t35 * t70;
t5 = t17 * t41 - t38 * t58;
t19 = -t40 * t53 + t69 * t70;
t67 = t35 * t40;
t9 = t19 * t41 + t38 * t67;
t79 = g(1) * t9 + g(2) * t5 + g(3) * t15;
t68 = t35 * t39;
t14 = t38 * t68 - t41 * t61;
t4 = t17 * t38 + t41 * t58;
t8 = t19 * t38 - t40 * t66;
t78 = g(1) * t8 + g(2) * t4 + g(3) * t14;
t73 = g(3) * t35;
t34 = sin(pkin(10));
t72 = t34 * pkin(4);
t71 = rSges(4,3) + pkin(8);
t65 = t70 * pkin(1) + pkin(7) * t67;
t57 = t35 * t69;
t64 = pkin(2) * t57 + pkin(8) * t68;
t60 = t19 * pkin(2) + t65;
t59 = pkin(8) + t72;
t56 = t38 * t69;
t55 = t41 * t69;
t54 = -t40 * pkin(1) + pkin(7) * t58;
t52 = t35 * t56;
t51 = t35 * t55;
t50 = -t17 * pkin(2) + t54;
t49 = t61 * t69;
t48 = -rSges(4,1) * t41 + rSges(4,2) * t38;
t36 = cos(pkin(10));
t47 = -rSges(5,1) * t36 + rSges(5,2) * t34 - pkin(3);
t46 = rSges(5,1) * t34 + rSges(5,2) * t36 + pkin(8);
t29 = pkin(4) * t36 + pkin(3);
t33 = pkin(10) + qJ(5);
t30 = sin(t33);
t31 = cos(t33);
t45 = rSges(6,1) * t31 - rSges(6,2) * t30 + t29;
t44 = rSges(6,1) * t30 + rSges(6,2) * t31 + t59;
t43 = -t38 * t62 + t41 * t47;
t42 = -t38 * t63 - t41 * t45;
t18 = t39 * t70 + t40 * t49;
t16 = t39 * t40 - t49 * t70;
t12 = t18 * pkin(2);
t10 = t16 * pkin(2);
t3 = t18 * t30 + t31 * t9;
t2 = t18 * t31 - t30 * t9;
t1 = [-m(2) * (g(1) * (-t40 * rSges(2,1) - rSges(2,2) * t70) + g(2) * (rSges(2,1) * t70 - t40 * rSges(2,2))) - m(3) * (g(1) * (-t17 * rSges(3,1) + t16 * rSges(3,2) + rSges(3,3) * t58 + t54) + g(2) * (rSges(3,1) * t19 - rSges(3,2) * t18 + rSges(3,3) * t67 + t65)) - m(4) * (g(1) * (-rSges(4,1) * t5 + rSges(4,2) * t4 - t16 * t71 + t50) + g(2) * (rSges(4,1) * t9 - rSges(4,2) * t8 + t18 * t71 + t60)) - m(5) * (g(1) * (-t5 * pkin(3) - t16 * pkin(8) + (-t16 * t34 - t36 * t5) * rSges(5,1) + (-t16 * t36 + t34 * t5) * rSges(5,2) - t62 * t4 + t50) + g(2) * (t18 * pkin(8) + t9 * pkin(3) + (t18 * t34 + t36 * t9) * rSges(5,1) + (t18 * t36 - t34 * t9) * rSges(5,2) + t62 * t8 + t60)) - m(6) * (g(1) * (-t16 * t44 - t4 * t63 - t45 * t5 + t50) + g(2) * (rSges(6,1) * t3 + rSges(6,2) * t2 + t18 * t59 + t29 * t9 + t63 * t8 + t60)), -m(3) * (g(1) * (-rSges(3,1) * t18 - rSges(3,2) * t19) + g(2) * (-rSges(3,1) * t16 - rSges(3,2) * t17) + (rSges(3,1) * t69 - rSges(3,2) * t39) * t73) - m(4) * (g(1) * (t18 * t48 + t19 * t71 - t12) + g(2) * (t16 * t48 + t17 * t71 - t10) + g(3) * t64 + (rSges(4,1) * t55 - rSges(4,2) * t56 + rSges(4,3) * t39) * t73) - m(5) * (g(1) * (t18 * t43 + t19 * t46 - t12) + g(2) * (t16 * t43 + t17 * t46 - t10) + g(3) * (pkin(3) * t51 + t64 + t62 * t52 + ((t34 * t39 + t36 * t55) * rSges(5,1) + (-t34 * t55 + t36 * t39) * rSges(5,2)) * t35)) - m(6) * (g(1) * (t18 * t42 + t19 * t44 - t12) + g(2) * (t16 * t42 + t17 * t44 - t10) + g(3) * (t29 * t51 + t68 * t72 + t64 + t63 * t52 + ((t30 * t39 + t31 * t55) * rSges(6,1) + (-t30 * t55 + t31 * t39) * rSges(6,2)) * t35)), -m(4) * (g(1) * (-rSges(4,1) * t8 - rSges(4,2) * t9) + g(2) * (-rSges(4,1) * t4 - rSges(4,2) * t5) + g(3) * (-rSges(4,1) * t14 - rSges(4,2) * t15)) - m(5) * (t47 * t78 + t62 * t79) - m(6) * (-t45 * t78 + t63 * t79), (-m(5) - m(6)) * t78, -m(6) * (g(1) * (rSges(6,1) * t2 - rSges(6,2) * t3) + g(2) * ((t16 * t31 - t30 * t5) * rSges(6,1) + (-t16 * t30 - t31 * t5) * rSges(6,2)) + g(3) * ((-t15 * t30 - t31 * t57) * rSges(6,1) + (-t15 * t31 + t30 * t57) * rSges(6,2)))];
taug = t1(:);
