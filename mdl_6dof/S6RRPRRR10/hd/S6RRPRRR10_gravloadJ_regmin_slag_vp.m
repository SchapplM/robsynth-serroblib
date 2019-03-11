% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 14:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR10_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRRR10_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t31 = sin(qJ(2));
t32 = sin(qJ(1));
t34 = cos(qJ(2));
t43 = cos(pkin(6));
t58 = cos(qJ(1));
t39 = t43 * t58;
t14 = t32 * t31 - t34 * t39;
t26 = qJ(5) + qJ(6);
t23 = sin(t26);
t24 = cos(t26);
t15 = t31 * t39 + t32 * t34;
t25 = pkin(12) + qJ(4);
t21 = sin(t25);
t22 = cos(t25);
t28 = sin(pkin(6));
t42 = t28 * t58;
t8 = t15 * t22 - t21 * t42;
t63 = -t14 * t24 + t8 * t23;
t62 = t14 * t23 + t8 * t24;
t30 = sin(qJ(5));
t33 = cos(qJ(5));
t61 = -t14 * t33 + t8 * t30;
t60 = t14 * t30 + t8 * t33;
t59 = g(3) * t28;
t53 = t22 * t23;
t52 = t22 * t24;
t51 = t22 * t30;
t50 = t22 * t33;
t49 = t22 * t34;
t48 = t28 * t31;
t47 = t28 * t32;
t46 = t28 * t34;
t45 = t30 * t34;
t44 = t33 * t34;
t41 = t32 * t43;
t16 = t58 * t31 + t34 * t41;
t40 = g(1) * t14 - g(2) * t16;
t17 = -t31 * t41 + t58 * t34;
t10 = -t17 * t21 + t22 * t47;
t37 = t15 * t21 + t22 * t42;
t38 = g(1) * t10 - g(2) * t37 + g(3) * (-t21 * t48 + t43 * t22);
t36 = -g(1) * t16 - g(2) * t14 + g(3) * t46;
t35 = g(1) * t17 + g(2) * t15 + g(3) * t48;
t29 = cos(pkin(12));
t27 = sin(pkin(12));
t13 = t43 * t21 + t22 * t48;
t11 = t17 * t22 + t21 * t47;
t6 = t11 * t33 + t16 * t30;
t5 = -t11 * t30 + t16 * t33;
t4 = t11 * t24 + t16 * t23;
t3 = -t11 * t23 + t16 * t24;
t2 = g(1) * t4 + g(2) * t62 - g(3) * (-t13 * t24 + t23 * t46);
t1 = -g(1) * t3 + g(2) * t63 - g(3) * (-t13 * t23 - t24 * t46);
t7 = [0, g(1) * t32 - g(2) * t58, g(1) * t58 + g(2) * t32, 0, 0, 0, 0, 0, g(1) * t15 - g(2) * t17, -t40, -g(1) * (-t15 * t29 + t27 * t42) - g(2) * (t17 * t29 + t27 * t47) -g(1) * (t15 * t27 + t29 * t42) - g(2) * (-t17 * t27 + t29 * t47) t40, -g(1) * (-t32 * pkin(1) - t15 * pkin(2) + pkin(8) * t42 - t14 * qJ(3)) - g(2) * (t58 * pkin(1) + t17 * pkin(2) + pkin(8) * t47 + t16 * qJ(3)) 0, 0, 0, 0, 0, g(1) * t8 - g(2) * t11, -g(1) * t37 - g(2) * t10, 0, 0, 0, 0, 0, g(1) * t60 - g(2) * t6, -g(1) * t61 - g(2) * t5, 0, 0, 0, 0, 0, g(1) * t62 - g(2) * t4, -g(1) * t63 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, -t36, t35, -t36 * t29, t36 * t27, -t35, -g(1) * (-t16 * pkin(2) + t17 * qJ(3)) - g(2) * (-t14 * pkin(2) + t15 * qJ(3)) - (pkin(2) * t34 + qJ(3) * t31) * t59, 0, 0, 0, 0, 0, -t36 * t22, t36 * t21, 0, 0, 0, 0, 0, -g(1) * (-t16 * t50 + t17 * t30) - g(2) * (-t14 * t50 + t15 * t30) - (t22 * t44 + t30 * t31) * t59, -g(1) * (t16 * t51 + t17 * t33) - g(2) * (t14 * t51 + t15 * t33) - (-t22 * t45 + t31 * t33) * t59, 0, 0, 0, 0, 0, -g(1) * (-t16 * t52 + t17 * t23) - g(2) * (-t14 * t52 + t15 * t23) - (t23 * t31 + t24 * t49) * t59, -g(1) * (t16 * t53 + t17 * t24) - g(2) * (t14 * t53 + t15 * t24) - (-t23 * t49 + t24 * t31) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, g(1) * t11 + g(2) * t8 + g(3) * t13, 0, 0, 0, 0, 0, -t38 * t33, t38 * t30, 0, 0, 0, 0, 0, -t38 * t24, t38 * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t61 - g(3) * (-t13 * t30 - t28 * t44) g(1) * t6 + g(2) * t60 - g(3) * (-t13 * t33 + t28 * t45) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t7;
