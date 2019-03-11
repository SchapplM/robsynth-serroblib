% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRP7
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
% Datum: 2019-03-09 06:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t45 = cos(qJ(4));
t32 = t45 * pkin(4) + pkin(3);
t38 = pkin(10) + qJ(3);
t34 = cos(t38);
t21 = t34 * t32;
t33 = sin(t38);
t47 = -pkin(9) - pkin(8);
t77 = -t33 * t47 + t21;
t44 = sin(qJ(1));
t46 = cos(qJ(1));
t25 = g(1) * t46 + g(2) * t44;
t49 = -g(3) * t34 + t25 * t33;
t74 = g(3) * t33;
t39 = qJ(4) + qJ(5);
t35 = sin(t39);
t72 = t33 * t35;
t36 = cos(t39);
t71 = t33 * t36;
t69 = t44 * t35;
t68 = t44 * t36;
t43 = sin(qJ(4));
t67 = t44 * t43;
t66 = t44 * t45;
t65 = t46 * t35;
t64 = t46 * t36;
t42 = -pkin(7) - qJ(2);
t63 = t46 * t42;
t62 = t46 * t43;
t61 = t46 * t45;
t60 = t34 * t62;
t10 = t34 * t69 + t64;
t11 = t34 * t68 - t65;
t59 = -t10 * pkin(5) + t11 * qJ(6);
t12 = t34 * t65 - t68;
t13 = t34 * t64 + t69;
t58 = -t12 * pkin(5) + t13 * qJ(6);
t41 = cos(pkin(10));
t31 = t41 * pkin(2) + pkin(1);
t23 = t46 * t31;
t57 = -t44 * t42 + t23;
t56 = t34 * pkin(3) + t33 * pkin(8);
t54 = g(1) * t10 - g(2) * t12;
t24 = g(1) * t44 - g(2) * t46;
t53 = pkin(5) * t36 + qJ(6) * t35;
t15 = t34 * t67 + t61;
t1 = g(1) * t12 + g(2) * t10 + g(3) * t72;
t3 = g(1) * t13 + g(2) * t11 + g(3) * t71;
t50 = pkin(4) * t67 + t77 * t46 + t57;
t7 = t25 * t34 + t74;
t48 = pkin(4) * t62 - t63 + (-t31 - t77) * t44;
t29 = pkin(4) * t66;
t20 = qJ(6) * t71;
t18 = t34 * t61 + t67;
t17 = -t60 + t66;
t16 = -t34 * t66 + t62;
t14 = t24 * t33;
t6 = t49 * t36;
t5 = t49 * t35;
t4 = g(1) * t11 - g(2) * t13;
t2 = [0, 0, 0, 0, 0, 0, t24, t25, 0, 0, 0, 0, 0, 0, 0, 0, t24 * t41, -t24 * sin(pkin(10)) -t25, -g(1) * (-t44 * pkin(1) + t46 * qJ(2)) - g(2) * (t46 * pkin(1) + t44 * qJ(2)) 0, 0, 0, 0, 0, 0, t24 * t34, -t14, -t25, -g(1) * (-t44 * t31 - t63) - g(2) * t57, 0, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, -g(1) * t15 - g(2) * t17, t14, -g(2) * t23 + (g(1) * t42 - g(2) * t56) * t46 + (-g(1) * (-t31 - t56) + g(2) * t42) * t44, 0, 0, 0, 0, 0, 0, t4, -t54, t14, -g(1) * t48 - g(2) * t50, 0, 0, 0, 0, 0, 0, t4, t14, t54, -g(1) * (-t11 * pkin(5) - t10 * qJ(6) + t48) - g(2) * (t13 * pkin(5) + t12 * qJ(6) + t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t7, 0, 0, 0, 0, 0, 0, 0, 0, t49 * t45, -t49 * t43, -t7, -g(3) * t56 + t25 * (pkin(3) * t33 - pkin(8) * t34) 0, 0, 0, 0, 0, 0, t6, -t5, -t7, -g(3) * t77 + t25 * (t32 * t33 + t34 * t47) 0, 0, 0, 0, 0, 0, t6, -t7, t5, -g(3) * t21 + (-g(3) * t53 + t25 * t47) * t34 + (g(3) * t47 + t25 * (t32 + t53)) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t17 + g(2) * t15 + t43 * t74, g(1) * t18 - g(2) * t16 + t45 * t74, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, -g(1) * t29 + (g(2) * t61 + t43 * t7) * pkin(4), 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * (-pkin(4) * t60 + t29 + t58) - g(2) * (-pkin(4) * t15 + t59) - g(3) * (t20 + (-pkin(4) * t43 - pkin(5) * t35) * t33); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * t58 - g(2) * t59 - g(3) * (-pkin(5) * t72 + t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
