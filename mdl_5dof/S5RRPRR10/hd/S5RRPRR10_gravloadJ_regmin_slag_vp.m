% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:26:39
% EndTime: 2019-12-31 20:26:40
% DurationCPUTime: 0.33s
% Computational Cost: add. (232->65), mult. (616->126), div. (0->0), fcn. (795->12), ass. (0->52)
t34 = cos(pkin(5));
t37 = sin(qJ(2));
t42 = cos(qJ(1));
t55 = t42 * t37;
t38 = sin(qJ(1));
t41 = cos(qJ(2));
t58 = t38 * t41;
t23 = -t34 * t58 - t55;
t54 = t42 * t41;
t59 = t38 * t37;
t46 = t34 * t54 - t59;
t32 = sin(pkin(5));
t66 = g(3) * t32;
t71 = -g(1) * t23 - g(2) * t46 - t41 * t66;
t35 = sin(qJ(5));
t39 = cos(qJ(5));
t36 = sin(qJ(4));
t40 = cos(qJ(4));
t31 = sin(pkin(10));
t33 = cos(pkin(10));
t48 = t41 * t31 + t37 * t33;
t19 = t48 * t34;
t25 = t37 * t31 - t41 * t33;
t50 = t42 * t19 - t38 * t25;
t62 = t32 * t42;
t4 = -t36 * t62 + t40 * t50;
t45 = t25 * t34;
t9 = -t38 * t48 - t42 * t45;
t70 = t4 * t35 + t9 * t39;
t69 = -t9 * t35 + t4 * t39;
t63 = t32 * t38;
t61 = t35 * t40;
t57 = t39 * t40;
t52 = g(1) * t42 + g(2) * t38;
t51 = g(1) * t38 - g(2) * t42;
t49 = -t38 * t19 - t42 * t25;
t47 = t36 * t50 + t40 * t62;
t18 = t48 * t32;
t6 = -t36 * t49 + t40 * t63;
t44 = g(1) * t6 - g(2) * t47 + g(3) * (-t18 * t36 + t34 * t40);
t12 = t38 * t45 - t42 * t48;
t17 = t25 * t32;
t43 = g(1) * t12 + g(2) * t9 - g(3) * t17;
t30 = t41 * pkin(2) + pkin(1);
t24 = -t34 * t59 + t54;
t22 = -t34 * t55 - t58;
t20 = t34 * t37 * pkin(2) + (-pkin(7) - qJ(3)) * t32;
t15 = t18 * t40 + t34 * t36;
t7 = t36 * t63 + t40 * t49;
t2 = -t12 * t35 + t7 * t39;
t1 = -t12 * t39 - t7 * t35;
t3 = [0, t51, t52, 0, 0, 0, 0, 0, -g(1) * t22 - g(2) * t24, g(1) * t46 - g(2) * t23, -t52 * t32, -g(1) * (-t42 * t20 - t38 * t30) - g(2) * (-t38 * t20 + t42 * t30), 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t7, -g(1) * t47 - g(2) * t6, 0, 0, 0, 0, 0, g(1) * t69 - g(2) * t2, -g(1) * t70 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, t71, g(1) * t24 - g(2) * t22 + t37 * t66, 0, t71 * pkin(2), 0, 0, 0, 0, 0, -t43 * t40, t43 * t36, 0, 0, 0, 0, 0, -g(1) * (t12 * t57 + t35 * t49) - g(2) * (t35 * t50 + t9 * t57) - g(3) * (-t17 * t57 + t18 * t35), -g(1) * (-t12 * t61 + t39 * t49) - g(2) * (t39 * t50 - t9 * t61) - g(3) * (t17 * t61 + t18 * t39); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t34 - t51 * t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, g(1) * t7 + g(2) * t4 + g(3) * t15, 0, 0, 0, 0, 0, -t44 * t39, t44 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t70 - g(3) * (-t15 * t35 + t17 * t39), g(1) * t2 + g(2) * t69 - g(3) * (-t15 * t39 - t17 * t35);];
taug_reg = t3;
