% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR1_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:56:49
% EndTime: 2019-05-05 21:56:50
% DurationCPUTime: 0.31s
% Computational Cost: add. (410->81), mult. (294->98), div. (0->0), fcn. (277->12), ass. (0->58)
t35 = qJ(3) + qJ(4);
t28 = pkin(11) + t35;
t22 = sin(t28);
t19 = t22 * pkin(9);
t23 = cos(t28);
t61 = -t23 * pkin(5) - t19;
t34 = qJ(1) + pkin(10);
t26 = sin(t34);
t27 = cos(t34);
t14 = g(1) * t27 + g(2) * t26;
t3 = -g(3) * t23 + t14 * t22;
t42 = -pkin(8) - pkin(7);
t29 = sin(t35);
t60 = pkin(4) * t29;
t59 = pkin(5) * t22;
t58 = g(3) * t22;
t38 = sin(qJ(1));
t56 = t38 * pkin(1);
t55 = t23 * t27;
t36 = sin(qJ(6));
t54 = t26 * t36;
t39 = cos(qJ(6));
t53 = t26 * t39;
t52 = t27 * t36;
t51 = t27 * t39;
t30 = cos(t35);
t24 = pkin(4) * t30;
t40 = cos(qJ(3));
t31 = t40 * pkin(3);
t49 = t24 + t31;
t17 = pkin(2) + t49;
t41 = cos(qJ(1));
t32 = t41 * pkin(1);
t50 = t27 * t17 + t32;
t48 = t24 - t61;
t37 = sin(qJ(3));
t18 = -t37 * pkin(3) - t60;
t47 = t18 - t59;
t46 = -t59 - t60;
t13 = g(1) * t26 - g(2) * t27;
t45 = g(1) * t38 - g(2) * t41;
t33 = -qJ(5) + t42;
t44 = -t27 * t33 - t56;
t5 = -g(3) * t30 + t14 * t29;
t43 = -g(3) * t40 + t14 * t37;
t25 = t31 + pkin(2);
t16 = pkin(9) * t55;
t15 = t26 * t23 * pkin(9);
t11 = t23 * t51 + t54;
t10 = -t23 * t52 + t53;
t9 = -t23 * t53 + t52;
t8 = t23 * t54 + t51;
t7 = t13 * t22;
t6 = g(3) * t29 + t14 * t30;
t4 = t14 * t23 + t58;
t2 = t3 * t39;
t1 = t3 * t36;
t12 = [0, 0, 0, 0, 0, 0, t45, g(1) * t41 + g(2) * t38, 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, t45 * pkin(1), 0, 0, 0, 0, 0, 0, t13 * t40, -t13 * t37, -t14, -g(1) * (-t26 * pkin(2) + t27 * pkin(7) - t56) - g(2) * (t27 * pkin(2) + t26 * pkin(7) + t32) 0, 0, 0, 0, 0, 0, t13 * t30, -t13 * t29, -t14, -g(1) * (-t26 * t25 - t27 * t42 - t56) - g(2) * (t27 * t25 - t26 * t42 + t32) 0, 0, 0, 0, 0, 0, t13 * t23, -t7, -t14, -g(1) * (-t26 * t17 + t44) - g(2) * (-t26 * t33 + t50) 0, 0, 0, 0, 0, 0, -g(1) * t9 - g(2) * t11, -g(1) * t8 - g(2) * t10, t7, -g(1) * t44 - g(2) * (pkin(5) * t55 + t27 * t19 + t50) + (-g(1) * (-t17 + t61) + g(2) * t33) * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, g(3) * t37 + t14 * t40, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t43 * pkin(3), 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(3) * t49 - t14 * t18, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t47 * t27 + t16) - g(2) * (t47 * t26 + t15) - g(3) * (t31 + t48); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(4), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t46 * t27 + t16) - g(2) * (t46 * t26 + t15) - g(3) * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t10 + g(2) * t8 + t36 * t58, g(1) * t11 - g(2) * t9 + t39 * t58, 0, 0;];
taug_reg  = t12;
