% Calculate inertial parameters regressor of gravitation load for
% S6RPRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t34 = cos(pkin(10));
t23 = t34 * pkin(4) + pkin(3);
t31 = pkin(9) + qJ(3);
t28 = cos(t31);
t13 = t28 * t23;
t26 = sin(t31);
t36 = -pkin(8) - qJ(4);
t64 = -t26 * t36 + t13;
t38 = sin(qJ(1));
t39 = cos(qJ(1));
t17 = g(1) * t39 + g(2) * t38;
t5 = -g(3) * t28 + t17 * t26;
t61 = g(3) * t26;
t30 = pkin(10) + qJ(5);
t25 = sin(t30);
t58 = t38 * t25;
t27 = cos(t30);
t57 = t38 * t27;
t32 = sin(pkin(10));
t56 = t38 * t32;
t55 = t38 * t34;
t54 = t39 * t25;
t53 = t39 * t27;
t52 = t39 * t32;
t51 = t39 * t34;
t37 = -pkin(7) - qJ(2);
t50 = t39 * t37;
t35 = cos(pkin(9));
t24 = t35 * pkin(2) + pkin(1);
t15 = t39 * t24;
t49 = -t38 * t37 + t15;
t7 = t28 * t58 + t53;
t9 = t28 * t54 - t57;
t48 = g(1) * t7 - g(2) * t9;
t16 = g(1) * t38 - g(2) * t39;
t47 = t28 * pkin(3) + t26 * qJ(4);
t45 = pkin(5) * t27 + qJ(6) * t25;
t1 = g(1) * t9 + g(2) * t7 + t25 * t61;
t10 = t28 * t53 + t58;
t8 = t28 * t57 - t54;
t42 = g(1) * t10 + g(2) * t8 + t27 * t61;
t41 = pkin(4) * t56 + t64 * t39 + t49;
t40 = pkin(4) * t52 - t50 + (-t24 - t64) * t38;
t11 = t16 * t26;
t6 = t17 * t28 + t61;
t4 = t5 * t27;
t3 = t5 * t25;
t2 = g(1) * t8 - g(2) * t10;
t12 = [0, 0, 0, 0, 0, 0, t16, t17, 0, 0, 0, 0, 0, 0, 0, 0, t16 * t35, -t16 * sin(pkin(9)) -t17, -g(1) * (-t38 * pkin(1) + t39 * qJ(2)) - g(2) * (t39 * pkin(1) + t38 * qJ(2)) 0, 0, 0, 0, 0, 0, t16 * t28, -t11, -t17, -g(1) * (-t38 * t24 - t50) - g(2) * t49, 0, 0, 0, 0, 0, 0, -g(1) * (-t28 * t55 + t52) - g(2) * (t28 * t51 + t56) -g(1) * (t28 * t56 + t51) - g(2) * (-t28 * t52 + t55) t11, -g(2) * t15 + (g(1) * t37 - g(2) * t47) * t39 + (-g(1) * (-t24 - t47) + g(2) * t37) * t38, 0, 0, 0, 0, 0, 0, t2, -t48, t11, -g(1) * t40 - g(2) * t41, 0, 0, 0, 0, 0, 0, t2, t11, t48, -g(1) * (-t8 * pkin(5) - t7 * qJ(6) + t40) - g(2) * (t10 * pkin(5) + t9 * qJ(6) + t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t34, -t5 * t32, -t6, -g(3) * t47 + t17 * (pkin(3) * t26 - qJ(4) * t28) 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(3) * t64 + t17 * (t23 * t26 + t28 * t36) 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(3) * t13 + (-g(3) * t45 + t17 * t36) * t28 + (g(3) * t36 + t17 * (t23 + t45)) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t42, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t42, -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (-t7 * pkin(5) + t8 * qJ(6)) - (-pkin(5) * t25 + qJ(6) * t27) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t12;
