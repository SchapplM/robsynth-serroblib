% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 07:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR7_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t33 = qJ(3) + qJ(4);
t29 = qJ(5) + t33;
t23 = sin(t29);
t24 = cos(t29);
t47 = -t23 * pkin(5) + t24 * pkin(10);
t36 = sin(qJ(1));
t39 = cos(qJ(1));
t60 = g(2) * t39;
t19 = g(1) * t36 - t60;
t64 = -g(3) * t23 + t19 * t24;
t40 = -pkin(8) - pkin(7);
t25 = sin(t33);
t63 = pkin(4) * t25;
t26 = cos(t33);
t62 = pkin(4) * t26;
t61 = pkin(5) * t24;
t58 = g(3) * t24;
t35 = sin(qJ(3));
t56 = t35 * pkin(3);
t55 = t23 * t36;
t34 = sin(qJ(6));
t54 = t36 * t34;
t37 = cos(qJ(6));
t53 = t36 * t37;
t52 = t39 * t34;
t51 = t39 * t37;
t50 = pkin(10) * t55 + t36 * t61;
t49 = t39 * pkin(1) + t36 * qJ(2);
t28 = t39 * qJ(2);
t48 = -t36 * pkin(1) + t28;
t46 = -pkin(10) * t23 - t61;
t20 = g(1) * t39 + g(2) * t36;
t14 = t56 + t63;
t32 = -pkin(9) + t40;
t44 = t39 * t14 + t36 * t32 + t48;
t43 = t36 * t14 - t39 * t32 + t49;
t42 = t47 - t63;
t6 = g(3) * t25 - t19 * t26;
t38 = cos(qJ(3));
t41 = g(3) * t35 - t19 * t38;
t15 = t38 * pkin(3) + t62;
t11 = t23 * t51 - t54;
t10 = t23 * t52 + t53;
t9 = t23 * t53 + t52;
t8 = -t23 * t54 + t51;
t7 = t20 * t24;
t5 = g(3) * t26 + t19 * t25;
t3 = g(1) * t55 - t23 * t60 + t58;
t2 = t64 * t37;
t1 = t64 * t34;
t4 = [0, 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t20, -g(1) * t48 - g(2) * t49, 0, 0, 0, 0, 0, 0, -t20 * t35, -t20 * t38, t19, -g(1) * (t28 + (-pkin(1) - pkin(7)) * t36) - g(2) * (t39 * pkin(7) + t49) 0, 0, 0, 0, 0, 0, -t20 * t25, -t20 * t26, t19, -g(1) * (t39 * t56 + t28 + (-pkin(1) + t40) * t36) - g(2) * (t36 * t56 - t39 * t40 + t49) 0, 0, 0, 0, 0, 0, -t20 * t23, -t7, t19, -g(1) * t44 - g(2) * t43, 0, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t9, g(1) * t10 - g(2) * t8, t7, -g(1) * (-t39 * t47 + t44) - g(2) * (-t36 * t47 + t43); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, g(3) * t38 + t19 * t35, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, t41 * pkin(3), 0, 0, 0, 0, 0, 0, -t64, t3, 0, g(3) * t14 - t19 * t15, 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(1) * (t36 * t15 + t50) - g(3) * (t42 - t56) - (-t15 + t46) * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t64, t3, 0, t6 * pkin(4), 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(1) * (t36 * t62 + t50) - g(3) * t42 - (t46 - t62) * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, t3, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(1) * t50 - g(3) * t47 - t46 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10 + t34 * t58, g(1) * t9 - g(2) * t11 + t37 * t58, 0, 0;];
taug_reg  = t4;
