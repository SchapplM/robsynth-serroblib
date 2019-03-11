% Calculate Gravitation load on the joints for
% S6PRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S6PRPRPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [7x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [7,3]), ...
  'S6PRPRPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [7x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:43:00
% EndTime: 2019-03-08 19:43:01
% DurationCPUTime: 0.67s
% Computational Cost: add. (459->122), mult. (765->180), div. (0->0), fcn. (881->12), ass. (0->56)
t35 = pkin(11) + qJ(4);
t33 = sin(t35);
t34 = cos(t35);
t77 = pkin(4) * t34 + qJ(5) * t33;
t69 = rSges(7,3) + pkin(9);
t76 = t69 * t34;
t42 = sin(qJ(2));
t44 = cos(qJ(2));
t37 = sin(pkin(10));
t58 = cos(pkin(6));
t55 = t37 * t58;
t57 = cos(pkin(10));
t22 = t57 * t42 + t44 * t55;
t47 = t58 * t57;
t20 = t37 * t42 - t44 * t47;
t71 = g(2) * t20;
t75 = -g(1) * t22 - t71;
t74 = -m(6) - m(7);
t38 = sin(pkin(6));
t70 = g(3) * t38;
t41 = sin(qJ(6));
t68 = t33 * t41;
t43 = cos(qJ(6));
t67 = t33 * t43;
t66 = t37 * t38;
t65 = t38 * t42;
t64 = t38 * t44;
t21 = t37 * t44 + t42 * t47;
t39 = cos(pkin(11));
t32 = t39 * pkin(3) + pkin(2);
t40 = -pkin(8) - qJ(3);
t63 = -t20 * t32 - t21 * t40;
t23 = -t42 * t55 + t57 * t44;
t62 = -t22 * t32 - t23 * t40;
t60 = rSges(4,3) + qJ(3);
t59 = rSges(6,3) + qJ(5);
t56 = -m(4) - m(5) + t74;
t54 = t38 * t57;
t53 = -t77 * t20 + t63;
t52 = -t77 * t22 + t62;
t25 = t32 * t64;
t51 = g(3) * (t77 * t64 + t25);
t50 = rSges(5,1) * t34 - rSges(5,2) * t33;
t49 = t41 * rSges(7,1) + t43 * rSges(7,2);
t48 = rSges(6,2) * t34 - rSges(6,3) * t33;
t46 = rSges(4,1) * t39 - rSges(4,2) * sin(pkin(11)) + pkin(2);
t17 = t58 * t33 + t34 * t65;
t16 = t33 * t65 - t58 * t34;
t15 = t16 * pkin(4);
t8 = t23 * t34 + t33 * t66;
t7 = t23 * t33 - t34 * t66;
t6 = t21 * t34 - t33 * t54;
t5 = t21 * t33 + t34 * t54;
t4 = t7 * pkin(4);
t3 = t5 * pkin(4);
t1 = [(-m(2) - m(3) + t56) * g(3), -m(3) * (g(1) * (-t22 * rSges(3,1) - t23 * rSges(3,2)) + g(2) * (-t20 * rSges(3,1) - t21 * rSges(3,2)) + (rSges(3,1) * t44 - rSges(3,2) * t42) * t70) - m(4) * (g(1) * (-t46 * t22 + t60 * t23) + g(2) * t60 * t21 - t46 * t71 + (t60 * t42 + t46 * t44) * t70) - m(5) * (g(1) * (t23 * rSges(5,3) - t50 * t22 + t62) + g(2) * (t21 * rSges(5,3) - t50 * t20 + t63) + g(3) * t25 + (t50 * t44 + (rSges(5,3) - t40) * t42) * t70) - m(6) * (g(1) * (t23 * rSges(6,1) + t48 * t22 + t52) + g(2) * (t21 * rSges(6,1) + t48 * t20 + t53) + t51 + (-t48 * t44 + (rSges(6,1) - t40) * t42) * t70) - m(7) * (g(1) * (t23 * pkin(5) + (-t22 * t68 + t23 * t43) * rSges(7,1) + (-t22 * t67 - t23 * t41) * rSges(7,2) + t52) + g(2) * (t21 * pkin(5) + (-t20 * t68 + t21 * t43) * rSges(7,1) + (-t20 * t67 - t21 * t41) * rSges(7,2) + t53) + t51 + t75 * t76 + ((t43 * rSges(7,1) - t41 * rSges(7,2) + pkin(5) - t40) * t42 + (t49 * t33 + t76) * t44) * t70) t56 * (-g(3) * t64 - t75) -m(5) * (g(1) * (-t7 * rSges(5,1) - t8 * rSges(5,2)) + g(2) * (-t5 * rSges(5,1) - t6 * rSges(5,2)) + g(3) * (-t16 * rSges(5,1) - t17 * rSges(5,2))) - m(6) * (g(1) * (t7 * rSges(6,2) + t59 * t8 - t4) + g(2) * (t5 * rSges(6,2) + t59 * t6 - t3) + g(3) * (t16 * rSges(6,2) + t59 * t17 - t15)) + (-g(1) * (-t69 * t7 - t4) - g(2) * (-t69 * t5 - t3) - g(3) * (-t69 * t16 - t15) - (g(1) * t8 + g(2) * t6 + g(3) * t17) * (qJ(5) + t49)) * m(7), t74 * (g(1) * t7 + g(2) * t5 + g(3) * t16) -m(7) * (g(1) * ((-t22 * t41 + t7 * t43) * rSges(7,1) + (-t22 * t43 - t7 * t41) * rSges(7,2)) + g(2) * ((-t20 * t41 + t5 * t43) * rSges(7,1) + (-t20 * t43 - t5 * t41) * rSges(7,2)) + g(3) * ((t16 * t43 + t41 * t64) * rSges(7,1) + (-t16 * t41 + t43 * t64) * rSges(7,2)))];
taug  = t1(:);
