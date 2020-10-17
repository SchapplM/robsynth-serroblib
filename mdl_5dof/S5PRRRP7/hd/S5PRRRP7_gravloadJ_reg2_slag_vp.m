% Calculate inertial parameters regressor of gravitation load for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:26
% EndTime: 2019-12-05 16:56:27
% DurationCPUTime: 0.35s
% Computational Cost: add. (286->86), mult. (739->131), div. (0->0), fcn. (897->10), ass. (0->54)
t29 = sin(pkin(5));
t61 = g(3) * t29;
t33 = sin(qJ(2));
t60 = t29 * t33;
t36 = cos(qJ(2));
t59 = t29 * t36;
t31 = sin(qJ(4));
t58 = t31 * t33;
t35 = cos(qJ(3));
t57 = t31 * t35;
t34 = cos(qJ(4));
t56 = t34 * t35;
t55 = t35 * t36;
t54 = pkin(2) * t59 + pkin(7) * t60;
t53 = cos(pkin(5));
t52 = cos(pkin(9));
t51 = sin(pkin(9));
t50 = pkin(4) * t31 + pkin(7);
t49 = g(3) * t54;
t42 = t53 * t52;
t15 = t51 * t33 - t36 * t42;
t13 = t15 * pkin(2);
t16 = t33 * t42 + t51 * t36;
t48 = t16 * pkin(7) - t13;
t41 = t53 * t51;
t17 = t52 * t33 + t36 * t41;
t14 = t17 * pkin(2);
t18 = -t33 * t41 + t52 * t36;
t47 = t18 * pkin(7) - t14;
t46 = t29 * t52;
t45 = t29 * t51;
t32 = sin(qJ(3));
t44 = pkin(3) * t35 + pkin(8) * t32;
t28 = pkin(4) * t34 + pkin(3);
t30 = -qJ(5) - pkin(8);
t43 = t28 * t35 - t30 * t32;
t11 = t18 * t32 - t35 * t45;
t19 = t32 * t60 - t53 * t35;
t9 = t16 * t32 + t35 * t46;
t40 = g(1) * t11 + g(2) * t9 + g(3) * t19;
t10 = t16 * t35 - t32 * t46;
t12 = t18 * t35 + t32 * t45;
t20 = t53 * t32 + t35 * t60;
t39 = g(1) * t12 + g(2) * t10 + g(3) * t20;
t38 = -g(1) * t17 - g(2) * t15 + g(3) * t59;
t37 = g(1) * t18 + g(2) * t16 + g(3) * t60;
t1 = -g(1) * (-t12 * t31 + t17 * t34) - g(2) * (-t10 * t31 + t15 * t34) - g(3) * (-t20 * t31 - t34 * t59);
t8 = t38 * t32;
t6 = t40 * t34;
t5 = t40 * t31;
t4 = -g(1) * (-t17 * t56 + t18 * t31) - g(2) * (-t15 * t56 + t16 * t31) - (t34 * t55 + t58) * t61;
t3 = -g(1) * (t17 * t57 + t18 * t34) - g(2) * (t15 * t57 + t16 * t34) - (-t31 * t55 + t33 * t34) * t61;
t2 = -g(1) * (-t12 * t34 - t17 * t31) - g(2) * (-t10 * t34 - t15 * t31) - g(3) * (-t20 * t34 + t31 * t59);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t37, 0, 0, 0, 0, 0, 0, 0, 0, -t38 * t35, t8, -t37, -g(1) * t47 - g(2) * t48 - t49, 0, 0, 0, 0, 0, 0, t4, t3, -t8, -g(1) * (-t44 * t17 + t47) - g(2) * (-t44 * t15 + t48) - g(3) * (t44 * t59 + t54), 0, 0, 0, 0, 0, 0, t4, t3, -t8, -g(1) * (-t43 * t17 + t50 * t18 - t14) - g(2) * (-t43 * t15 + t50 * t16 - t13) - t49 - (pkin(4) * t58 + t43 * t36) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t39, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t39, -g(1) * (-pkin(3) * t11 + pkin(8) * t12) - g(2) * (-pkin(3) * t9 + pkin(8) * t10) - g(3) * (-pkin(3) * t19 + pkin(8) * t20), 0, 0, 0, 0, 0, 0, t6, -t5, -t39, -g(1) * (-t11 * t28 - t12 * t30) - g(2) * (-t10 * t30 - t28 * t9) - g(3) * (-t19 * t28 - t20 * t30); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40;];
taug_reg = t7;
