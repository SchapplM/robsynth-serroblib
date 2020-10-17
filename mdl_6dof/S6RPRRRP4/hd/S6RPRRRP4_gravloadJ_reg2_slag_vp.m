% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP4_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:25:29
% EndTime: 2019-05-06 01:25:30
% DurationCPUTime: 0.33s
% Computational Cost: add. (438->81), mult. (400->100), div. (0->0), fcn. (388->10), ass. (0->52)
t34 = pkin(10) + qJ(3);
t30 = qJ(4) + t34;
t24 = sin(t30);
t25 = cos(t30);
t41 = cos(qJ(5));
t27 = t41 * pkin(5) + pkin(4);
t37 = -qJ(6) - pkin(9);
t69 = -t24 * t37 + t25 * t27;
t68 = t25 * pkin(4) + t24 * pkin(9);
t40 = sin(qJ(1));
t42 = cos(qJ(1));
t20 = g(1) * t42 + g(2) * t40;
t53 = t42 * t41;
t39 = sin(qJ(5));
t56 = t40 * t39;
t10 = t25 * t56 + t53;
t54 = t42 * t39;
t55 = t40 * t41;
t12 = -t25 * t54 + t55;
t59 = g(3) * t24;
t1 = -g(1) * t12 + g(2) * t10 + t39 * t59;
t7 = -g(3) * t25 + t20 * t24;
t28 = sin(t34);
t67 = pkin(3) * t28;
t66 = pkin(4) * t24;
t65 = pkin(9) * t25;
t29 = cos(t34);
t23 = pkin(3) * t29;
t36 = cos(pkin(10));
t26 = t36 * pkin(2) + pkin(1);
t15 = t23 + t26;
t14 = t42 * t15;
t61 = g(2) * t14;
t38 = -pkin(7) - qJ(2);
t33 = -pkin(8) + t38;
t50 = pkin(5) * t39 - t33;
t48 = -t66 - t67;
t19 = g(1) * t40 - g(2) * t42;
t45 = t24 * t27 + t25 * t37;
t43 = -g(3) * t29 + t20 * t28;
t18 = t42 * t65;
t17 = t40 * t65;
t13 = t25 * t53 + t56;
t11 = -t25 * t55 + t54;
t9 = t19 * t24;
t8 = t20 * t25 + t59;
t6 = t7 * t41;
t5 = t7 * t39;
t4 = -g(1) * t11 - g(2) * t13;
t3 = -g(1) * t10 - g(2) * t12;
t2 = g(1) * t13 - g(2) * t11 + t41 * t59;
t16 = [0, 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, 0, 0, 0, t19 * t36, -t19 * sin(pkin(10)) -t20, -g(1) * (-t40 * pkin(1) + t42 * qJ(2)) - g(2) * (t42 * pkin(1) + t40 * qJ(2)) 0, 0, 0, 0, 0, 0, t19 * t29, -t19 * t28, -t20, -g(1) * (-t40 * t26 - t42 * t38) - g(2) * (t42 * t26 - t40 * t38) 0, 0, 0, 0, 0, 0, t19 * t25, -t9, -t20, -g(1) * (-t40 * t15 - t42 * t33) - g(2) * (-t40 * t33 + t14) 0, 0, 0, 0, 0, 0, t4, t3, t9, -t61 + (g(1) * t33 - g(2) * t68) * t42 + (-g(1) * (-t15 - t68) + g(2) * t33) * t40, 0, 0, 0, 0, 0, 0, t4, t3, t9, -t61 + (-g(1) * t50 - g(2) * t69) * t42 + (-g(1) * (-t15 - t69) - g(2) * t50) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, g(3) * t28 + t20 * t29, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, t43 * pkin(3), 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (t48 * t42 + t18) - g(2) * (t48 * t40 + t17) - g(3) * (t23 + t68) 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * (t23 + t69) + t20 * (t45 + t67); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (-t42 * t66 + t18) - g(2) * (-t40 * t66 + t17) - g(3) * t68, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(3) * t69 + t20 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7;];
taug_reg  = t16;
