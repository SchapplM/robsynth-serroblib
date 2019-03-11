% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR1_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t37 = cos(qJ(2));
t27 = t37 * pkin(2);
t31 = qJ(2) + pkin(10);
t26 = qJ(4) + t31;
t21 = sin(t26);
t22 = cos(t26);
t45 = t22 * pkin(4) + t21 * qJ(5);
t57 = pkin(3) * cos(t31) + t27 + t45;
t36 = sin(qJ(1));
t38 = cos(qJ(1));
t17 = g(1) * t38 + g(2) * t36;
t5 = -g(3) * t22 + t17 * t21;
t56 = pkin(4) * t21;
t55 = g(3) * t21;
t30 = pkin(11) + qJ(6);
t24 = sin(t30);
t53 = t36 * t24;
t25 = cos(t30);
t52 = t36 * t25;
t32 = sin(pkin(11));
t51 = t36 * t32;
t33 = cos(pkin(11));
t50 = t36 * t33;
t49 = t38 * t24;
t48 = t38 * t25;
t47 = t38 * t32;
t46 = t38 * t33;
t34 = -qJ(3) - pkin(7);
t43 = qJ(5) * t22;
t35 = sin(qJ(2));
t42 = -pkin(3) * sin(t31) - t35 * pkin(2) - t56;
t16 = g(1) * t36 - g(2) * t38;
t41 = pkin(1) + t57;
t39 = -g(3) * t37 + t17 * t35;
t29 = -pkin(8) + t34;
t23 = t27 + pkin(1);
t15 = t38 * t43;
t14 = t36 * t43;
t11 = t16 * t21;
t10 = t22 * t48 + t53;
t9 = -t22 * t49 + t52;
t8 = -t22 * t52 + t49;
t7 = t22 * t53 + t48;
t6 = t17 * t22 + t55;
t4 = t5 * t33;
t3 = t5 * t32;
t2 = t5 * t25;
t1 = t5 * t24;
t12 = [0, t16, t17, 0, 0, 0, 0, 0, t16 * t37, -t16 * t35, -t17, -g(1) * (-t36 * t23 - t38 * t34) - g(2) * (t38 * t23 - t36 * t34) 0, 0, 0, 0, 0, t16 * t22, -t11, -g(1) * (-t22 * t50 + t47) - g(2) * (t22 * t46 + t51) -g(1) * (t22 * t51 + t46) - g(2) * (-t22 * t47 + t50) t11 (g(1) * t29 - g(2) * t41) * t38 + (g(1) * t41 + g(2) * t29) * t36, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9; 0, 0, 0, 0, 0, 0, 0, 0, t39, g(3) * t35 + t17 * t37, 0, t39 * pkin(2), 0, 0, 0, 0, 0, t5, t6, t4, -t3, -t6, -g(1) * (t42 * t38 + t15) - g(2) * (t42 * t36 + t14) - g(3) * t57, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, t4, -t3, -t6, -g(1) * (-t38 * t56 + t15) - g(2) * (-t36 * t56 + t14) - g(3) * t45, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t24 * t55, g(1) * t10 - g(2) * t8 + t25 * t55;];
taug_reg  = t12;
