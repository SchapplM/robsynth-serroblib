% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t43 = cos(qJ(4));
t29 = t43 * pkin(4) + pkin(3);
t35 = pkin(9) + qJ(3);
t32 = cos(t35);
t17 = t32 * t29;
t30 = sin(t35);
t39 = -qJ(5) - pkin(8);
t70 = -t30 * t39 + t17;
t42 = sin(qJ(1));
t44 = cos(qJ(1));
t21 = g(1) * t44 + g(2) * t42;
t5 = -g(3) * t32 + t21 * t30;
t67 = g(3) * t30;
t36 = qJ(4) + pkin(10);
t31 = sin(t36);
t64 = t42 * t31;
t33 = cos(t36);
t63 = t42 * t33;
t41 = sin(qJ(4));
t62 = t42 * t41;
t61 = t42 * t43;
t60 = t44 * t31;
t59 = t44 * t33;
t40 = -pkin(7) - qJ(2);
t58 = t44 * t40;
t57 = t44 * t41;
t56 = t44 * t43;
t55 = t32 * t57;
t38 = cos(pkin(9));
t28 = t38 * pkin(2) + pkin(1);
t19 = t44 * t28;
t54 = -t42 * t40 + t19;
t7 = t32 * t64 + t59;
t9 = t32 * t60 - t63;
t53 = g(1) * t7 - g(2) * t9;
t52 = t32 * pkin(3) + t30 * pkin(8);
t20 = g(1) * t42 - g(2) * t44;
t50 = pkin(5) * t33 + qJ(6) * t31;
t12 = t32 * t62 + t56;
t1 = g(1) * t9 + g(2) * t7 + t31 * t67;
t10 = t32 * t59 + t64;
t8 = t32 * t63 - t60;
t47 = g(1) * t10 + g(2) * t8 + t33 * t67;
t46 = pkin(4) * t62 + t70 * t44 + t54;
t6 = t21 * t32 + t67;
t45 = pkin(4) * t57 - t58 + (-t28 - t70) * t42;
t25 = pkin(4) * t61;
t15 = t32 * t56 + t62;
t14 = -t55 + t61;
t13 = -t32 * t61 + t57;
t11 = t20 * t30;
t4 = t5 * t33;
t3 = t5 * t31;
t2 = g(1) * t8 - g(2) * t10;
t16 = [0, 0, 0, 0, 0, 0, t20, t21, 0, 0, 0, 0, 0, 0, 0, 0, t20 * t38, -t20 * sin(pkin(9)) -t21, -g(1) * (-t42 * pkin(1) + t44 * qJ(2)) - g(2) * (t44 * pkin(1) + t42 * qJ(2)) 0, 0, 0, 0, 0, 0, t20 * t32, -t11, -t21, -g(1) * (-t42 * t28 - t58) - g(2) * t54, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t15, -g(1) * t12 - g(2) * t14, t11, -g(2) * t19 + (g(1) * t40 - g(2) * t52) * t44 + (-g(1) * (-t28 - t52) + g(2) * t40) * t42, 0, 0, 0, 0, 0, 0, t2, -t53, t11, -g(1) * t45 - g(2) * t46, 0, 0, 0, 0, 0, 0, t2, t11, t53, -g(1) * (-t8 * pkin(5) - t7 * qJ(6) + t45) - g(2) * (t10 * pkin(5) + t9 * qJ(6) + t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t43, -t5 * t41, -t6, -g(3) * t52 + t21 * (pkin(3) * t30 - pkin(8) * t32) 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(3) * t70 + t21 * (t29 * t30 + t32 * t39) 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(3) * t17 + (-g(3) * t50 + t21 * t39) * t32 + (g(3) * t39 + t21 * (t29 + t50)) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t14 + g(2) * t12 + t41 * t67, g(1) * t15 - g(2) * t13 + t43 * t67, 0, 0, 0, 0, 0, 0, 0, 0, t1, t47, 0, -g(1) * t25 + (g(2) * t56 + t6 * t41) * pkin(4), 0, 0, 0, 0, 0, 0, t1, 0, -t47, -g(1) * (-pkin(4) * t55 - t9 * pkin(5) + t10 * qJ(6) + t25) - g(2) * (-t12 * pkin(4) - t7 * pkin(5) + t8 * qJ(6)) - (-pkin(4) * t41 - pkin(5) * t31 + qJ(6) * t33) * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t16;
