% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:07
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR4_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t37 = pkin(11) + qJ(3);
t32 = qJ(4) + t37;
t29 = qJ(5) + t32;
t23 = sin(t29);
t24 = cos(t29);
t60 = pkin(5) * t24 + pkin(10) * t23;
t42 = sin(qJ(1));
t44 = cos(qJ(1));
t19 = g(1) * t44 + g(2) * t42;
t3 = -g(3) * t24 + t19 * t23;
t26 = sin(t32);
t59 = pkin(4) * t26;
t58 = pkin(5) * t23;
t57 = pkin(10) * t24;
t56 = g(3) * t23;
t39 = cos(pkin(11));
t28 = pkin(2) * t39 + pkin(1);
t41 = sin(qJ(6));
t54 = t42 * t41;
t43 = cos(qJ(6));
t53 = t42 * t43;
t52 = t44 * t41;
t51 = t44 * t43;
t40 = -pkin(7) - qJ(2);
t36 = -pkin(8) + t40;
t31 = cos(t37);
t25 = pkin(3) * t31;
t15 = t25 + t28;
t27 = cos(t32);
t22 = pkin(4) * t27;
t49 = t22 + t60;
t30 = sin(t37);
t14 = -pkin(3) * t30 - t59;
t48 = t14 - t58;
t47 = -t58 - t59;
t18 = g(1) * t42 - g(2) * t44;
t5 = -g(3) * t27 + t19 * t26;
t45 = -g(3) * t31 + t19 * t30;
t33 = -pkin(9) + t36;
t17 = t44 * t57;
t16 = t42 * t57;
t13 = t24 * t51 + t54;
t12 = -t24 * t52 + t53;
t11 = -t24 * t53 + t52;
t10 = t24 * t54 + t51;
t9 = t22 + t15;
t8 = t44 * t9;
t7 = t18 * t23;
t6 = g(3) * t26 + t19 * t27;
t4 = t19 * t24 + t56;
t2 = t3 * t43;
t1 = t3 * t41;
t20 = [0, 0, 0, 0, 0, 0, t18, t19, 0, 0, 0, 0, 0, 0, 0, 0, t18 * t39, -t18 * sin(pkin(11)) -t19, -g(1) * (-pkin(1) * t42 + qJ(2) * t44) - g(2) * (pkin(1) * t44 + qJ(2) * t42) 0, 0, 0, 0, 0, 0, t18 * t31, -t18 * t30, -t19, -g(1) * (-t28 * t42 - t40 * t44) - g(2) * (t28 * t44 - t40 * t42) 0, 0, 0, 0, 0, 0, t18 * t27, -t18 * t26, -t19, -g(1) * (-t15 * t42 - t36 * t44) - g(2) * (t15 * t44 - t36 * t42) 0, 0, 0, 0, 0, 0, t18 * t24, -t7, -t19, -g(1) * (-t44 * t33 - t42 * t9) - g(2) * (-t42 * t33 + t8) 0, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t13, -g(1) * t10 - g(2) * t12, t7, -g(2) * t8 + (g(1) * t33 - g(2) * t60) * t44 + (-g(1) * (-t60 - t9) + g(2) * t33) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, g(3) * t30 + t19 * t31, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t45 * pkin(3), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * (t22 + t25) - t19 * t14, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t44 * t48 + t17) - g(2) * (t42 * t48 + t16) - g(3) * (t25 + t49); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(4), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t44 * t47 + t17) - g(2) * (t42 * t47 + t16) - g(3) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (-t44 * t58 + t17) - g(2) * (-t42 * t58 + t16) - g(3) * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t12 + g(2) * t10 + t41 * t56, g(1) * t13 - g(2) * t11 + t43 * t56, 0, 0;];
taug_reg  = t20;
