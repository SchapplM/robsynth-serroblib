% Calculate Gravitation load on the joints for
% S5RPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug = S5RPRPR5_gravloadJ_floatb_twist_slag_vp1(qJ, g, ...
  pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp1: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp1: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp1: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp1: m has to be [6x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [6,3]), ...
  'S5RPRPR5_gravloadJ_floatb_twist_slag_vp1: rSges has to be [6x3] (double)');

%% Symbolic Calculation
% From gravload_joint_floatb_twist_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:41:26
% EndTime: 2020-01-03 11:41:28
% DurationCPUTime: 0.46s
% Computational Cost: add. (234->105), mult. (294->143), div. (0->0), fcn. (285->10), ass. (0->55)
t39 = sin(qJ(1));
t41 = cos(qJ(1));
t69 = g(2) * t41 + g(3) * t39;
t68 = -m(5) - m(6);
t36 = cos(pkin(8));
t34 = qJ(3) + pkin(9);
t28 = qJ(5) + t34;
t24 = cos(t28);
t50 = t41 * t24;
t23 = sin(t28);
t57 = t39 * t23;
t5 = -t36 * t57 - t50;
t51 = t41 * t23;
t56 = t39 * t24;
t6 = t36 * t56 - t51;
t67 = t5 * rSges(6,1) - t6 * rSges(6,2);
t7 = t36 * t51 - t56;
t8 = t36 * t50 + t57;
t66 = t7 * rSges(6,1) + t8 * rSges(6,2);
t35 = sin(pkin(8));
t65 = g(1) * t35;
t38 = sin(qJ(3));
t62 = t38 * pkin(3);
t60 = rSges(3,2) * t35;
t59 = t36 * t39;
t58 = t36 * t41;
t26 = sin(t34);
t55 = t39 * t26;
t27 = cos(t34);
t54 = t39 * t27;
t53 = t39 * t38;
t40 = cos(qJ(3));
t52 = t39 * t40;
t49 = t41 * t26;
t48 = t41 * t27;
t47 = t41 * t38;
t46 = t41 * t40;
t37 = -qJ(4) - pkin(6);
t31 = t40 * pkin(3);
t20 = pkin(4) * t27 + t31;
t45 = t41 * pkin(1) + t39 * qJ(2);
t44 = -rSges(6,1) * t23 - rSges(6,2) * t24;
t16 = t36 * t47 - t52;
t14 = -t36 * t53 - t46;
t43 = (t31 + pkin(2)) * t36 + (rSges(5,3) - t37) * t35;
t42 = (pkin(2) + t20) * t36 + (rSges(6,3) + pkin(7) - t37) * t35;
t30 = t39 * pkin(1);
t19 = pkin(4) * t26 + t62;
t17 = t36 * t46 + t53;
t15 = t36 * t52 - t47;
t12 = t36 * t48 + t55;
t11 = t36 * t49 - t54;
t10 = t36 * t54 - t49;
t9 = -t36 * t55 - t48;
t1 = [-m(2) * (g(2) * (t41 * rSges(2,1) - t39 * rSges(2,2)) + g(3) * (t39 * rSges(2,1) + t41 * rSges(2,2))) - m(3) * (g(2) * (t39 * rSges(3,3) + t45) + g(3) * (rSges(3,1) * t59 - t39 * t60 + t30) + (g(2) * (rSges(3,1) * t36 - t60) + g(3) * (-rSges(3,3) - qJ(2))) * t41) - m(4) * (g(2) * (t17 * rSges(4,1) - t16 * rSges(4,2) + pkin(2) * t58 + t45) + g(3) * (t15 * rSges(4,1) + t14 * rSges(4,2) + pkin(2) * t59 - t41 * qJ(2) + t30) + t69 * t35 * (rSges(4,3) + pkin(6))) - m(5) * (g(2) * (t12 * rSges(5,1) - t11 * rSges(5,2) + t45) + g(3) * (t10 * rSges(5,1) + t9 * rSges(5,2) + t30) + (g(2) * t62 + g(3) * t43) * t39 + (g(2) * t43 + g(3) * (-qJ(2) - t62)) * t41) - m(6) * (g(2) * (t8 * rSges(6,1) - t7 * rSges(6,2) + t45) + g(3) * (t6 * rSges(6,1) + t5 * rSges(6,2) + t30) + (g(2) * t19 + g(3) * t42) * t39 + (g(2) * t42 + g(3) * (-qJ(2) - t19)) * t41), -(-m(3) - m(4) + t68) * t69, -m(4) * (g(2) * (t14 * rSges(4,1) - t15 * rSges(4,2)) + g(3) * (t16 * rSges(4,1) + t17 * rSges(4,2))) - m(5) * (g(2) * (t9 * rSges(5,1) - t10 * rSges(5,2) + t14 * pkin(3)) + g(3) * (t11 * rSges(5,1) + t12 * rSges(5,2) + t16 * pkin(3))) - m(6) * (g(2) * (-t19 * t59 - t41 * t20 + t67) + g(3) * (t19 * t58 - t39 * t20 + t66)) + (-m(4) * (-rSges(4,1) * t38 - rSges(4,2) * t40) - m(5) * (-rSges(5,1) * t26 - rSges(5,2) * t27 - t62) - m(6) * (-t19 + t44)) * t65, t68 * (-g(1) * t36 + (g(2) * t39 - g(3) * t41) * t35), -m(6) * (g(2) * t67 + g(3) * t66 + t44 * t65)];
taug = t1(:);
