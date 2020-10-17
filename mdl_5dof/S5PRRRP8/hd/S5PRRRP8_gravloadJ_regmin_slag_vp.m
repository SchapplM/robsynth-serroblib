% Calculate minimal parameter regressor of gravitation load for
% S5PRRRP8
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
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:35
% EndTime: 2019-12-05 17:00:36
% DurationCPUTime: 0.26s
% Computational Cost: add. (262->71), mult. (703->113), div. (0->0), fcn. (879->10), ass. (0->52)
t37 = sin(pkin(9));
t41 = sin(qJ(2));
t44 = cos(qJ(2));
t56 = cos(pkin(9));
t57 = cos(pkin(5));
t52 = t57 * t56;
t26 = t37 * t41 - t44 * t52;
t54 = t37 * t57;
t28 = t56 * t41 + t44 * t54;
t69 = g(1) * t28 + g(2) * t26;
t27 = t37 * t44 + t41 * t52;
t29 = -t41 * t54 + t56 * t44;
t40 = sin(qJ(3));
t43 = cos(qJ(3));
t38 = sin(pkin(5));
t53 = t38 * t56;
t62 = t38 * t43;
t63 = t38 * t41;
t47 = g(3) * (-t40 * t63 + t57 * t43) + g(2) * (-t27 * t40 - t43 * t53) + g(1) * (-t29 * t40 + t37 * t62);
t61 = t38 * t44;
t39 = sin(qJ(4));
t60 = t39 * t43;
t42 = cos(qJ(4));
t59 = t42 * t43;
t58 = t42 * t44;
t55 = t39 * t61;
t51 = pkin(3) * t43 + pkin(8) * t40 + pkin(2);
t31 = t57 * t40 + t41 * t62;
t18 = t31 * t39 + t38 * t58;
t15 = t27 * t43 - t40 * t53;
t5 = t15 * t39 - t26 * t42;
t17 = t37 * t38 * t40 + t29 * t43;
t7 = t17 * t39 - t28 * t42;
t1 = g(1) * t7 + g(2) * t5 + g(3) * t18;
t19 = t31 * t42 - t55;
t6 = t15 * t42 + t26 * t39;
t8 = t17 * t42 + t28 * t39;
t49 = g(1) * t8 + g(2) * t6 + g(3) * t19;
t10 = -t26 * t60 - t27 * t42;
t12 = -t28 * t60 - t29 * t42;
t20 = -t42 * t63 + t43 * t55;
t48 = g(1) * t12 + g(2) * t10 + g(3) * t20;
t46 = g(1) * t17 + g(2) * t15 + g(3) * t31;
t45 = g(3) * t61 - t69;
t21 = (t39 * t41 + t43 * t58) * t38;
t13 = -t28 * t59 + t29 * t39;
t11 = -t26 * t59 + t27 * t39;
t9 = t45 * t40;
t4 = t47 * t42;
t3 = t47 * t39;
t2 = -g(1) * t13 - g(2) * t11 - g(3) * t21;
t14 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -t45, g(1) * t29 + g(2) * t27 + g(3) * t63, 0, 0, 0, 0, 0, -t45 * t43, t9, 0, 0, 0, 0, 0, t2, t48, t2, -t9, -t48, -g(1) * (t13 * pkin(4) + t29 * pkin(7) + t12 * qJ(5)) - g(2) * (t11 * pkin(4) + t27 * pkin(7) + t10 * qJ(5)) + t69 * t51 + (-t21 * pkin(4) - t20 * qJ(5) - (pkin(7) * t41 + t51 * t44) * t38) * g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t46, 0, 0, 0, 0, 0, -t4, t3, -t4, -t46, -t3, -t46 * pkin(8) - t47 * (pkin(4) * t42 + qJ(5) * t39 + pkin(3)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t49, t1, 0, -t49, -g(1) * (-t7 * pkin(4) + t8 * qJ(5)) - g(2) * (-t5 * pkin(4) + t6 * qJ(5)) - g(3) * (-t18 * pkin(4) + t19 * qJ(5)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t14;
