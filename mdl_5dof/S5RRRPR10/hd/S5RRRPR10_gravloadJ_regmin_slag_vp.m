% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:29:51
% EndTime: 2019-12-31 21:29:52
% DurationCPUTime: 0.34s
% Computational Cost: add. (210->74), mult. (444->136), div. (0->0), fcn. (542->12), ass. (0->48)
t27 = sin(qJ(2));
t28 = sin(qJ(1));
t31 = cos(qJ(2));
t32 = cos(qJ(1));
t44 = cos(pkin(5));
t40 = t32 * t44;
t12 = t28 * t27 - t31 * t40;
t25 = sin(qJ(5));
t29 = cos(qJ(5));
t13 = t27 * t40 + t28 * t31;
t22 = qJ(3) + pkin(10);
t20 = sin(t22);
t21 = cos(t22);
t23 = sin(pkin(5));
t47 = t23 * t32;
t5 = -t13 * t21 + t20 * t47;
t62 = t12 * t29 + t5 * t25;
t61 = -t12 * t25 + t5 * t29;
t60 = g(1) * t32 + g(2) * t28;
t26 = sin(qJ(3));
t30 = cos(qJ(3));
t38 = t13 * t26 + t30 * t47;
t50 = t23 * t27;
t41 = t28 * t44;
t15 = -t27 * t41 + t32 * t31;
t48 = t23 * t30;
t8 = -t15 * t26 + t28 * t48;
t59 = -g(3) * (-t26 * t50 + t44 * t30) + g(2) * t38 - g(1) * t8;
t55 = g(3) * t23;
t52 = t21 * t25;
t51 = t21 * t29;
t49 = t23 * t28;
t46 = t25 * t31;
t45 = t29 * t31;
t42 = t13 * t30 - t26 * t47;
t14 = t32 * t27 + t31 * t41;
t39 = g(1) * t12 - g(2) * t14;
t37 = g(1) * (-t15 * t20 + t21 * t49) + g(2) * (-t13 * t20 - t21 * t47) + g(3) * (-t20 * t50 + t44 * t21);
t35 = -g(1) * t14 - g(2) * t12 + t31 * t55;
t34 = g(1) * t15 + g(2) * t13 + g(3) * t50;
t24 = -qJ(4) - pkin(8);
t19 = t30 * pkin(3) + pkin(2);
t11 = t44 * t20 + t21 * t50;
t9 = t15 * t30 + t26 * t49;
t7 = t15 * t21 + t20 * t49;
t2 = t14 * t25 + t7 * t29;
t1 = t14 * t29 - t7 * t25;
t3 = [0, g(1) * t28 - g(2) * t32, t60, 0, 0, 0, 0, 0, g(1) * t13 - g(2) * t15, -t39, 0, 0, 0, 0, 0, g(1) * t42 - g(2) * t9, -g(1) * t38 - g(2) * t8, t39, -g(1) * (-t28 * pkin(1) + t12 * t24 - t13 * t19) - g(2) * (t32 * pkin(1) - t14 * t24 + t15 * t19) - t60 * t23 * (pkin(3) * t26 + pkin(7)), 0, 0, 0, 0, 0, -g(1) * t61 - g(2) * t2, g(1) * t62 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, -t35, t34, 0, 0, 0, 0, 0, -t35 * t30, t35 * t26, -t34, -g(1) * (-t14 * t19 - t15 * t24) - g(2) * (-t12 * t19 - t13 * t24) - (t19 * t31 - t24 * t27) * t55, 0, 0, 0, 0, 0, -g(1) * (-t14 * t51 + t15 * t25) - g(2) * (-t12 * t51 + t13 * t25) - (t21 * t45 + t25 * t27) * t55, -g(1) * (t14 * t52 + t15 * t29) - g(2) * (t12 * t52 + t13 * t29) - (-t21 * t46 + t27 * t29) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, g(1) * t9 + g(2) * t42 - g(3) * (-t44 * t26 - t27 * t48), 0, t59 * pkin(3), 0, 0, 0, 0, 0, -t37 * t29, t37 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t62 - g(3) * (-t11 * t25 - t23 * t45), g(1) * t2 - g(2) * t61 - g(3) * (-t11 * t29 + t23 * t46);];
taug_reg = t3;
