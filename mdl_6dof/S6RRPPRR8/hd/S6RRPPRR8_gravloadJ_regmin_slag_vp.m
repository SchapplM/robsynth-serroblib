% Calculate minimal parameter regressor of gravitation load for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t35 = sin(qJ(2));
t38 = cos(qJ(2));
t58 = t38 * pkin(2) + t35 * qJ(3);
t71 = -pkin(1) - t58;
t36 = sin(qJ(1));
t39 = cos(qJ(1));
t54 = g(1) * t39 + g(2) * t36;
t69 = t54 * t35;
t70 = -g(3) * t38 + t69;
t67 = pkin(2) * t35;
t66 = g(1) * t36;
t63 = g(3) * t35;
t61 = t36 * t38;
t32 = sin(pkin(10));
t60 = t39 * t32;
t33 = cos(pkin(10));
t59 = t39 * t33;
t57 = qJ(3) * t38;
t56 = t36 * pkin(7) - t71 * t39;
t12 = t32 * t61 + t59;
t14 = -t36 * t33 + t38 * t60;
t55 = g(1) * t12 - g(2) * t14;
t53 = -g(2) * t39 + t66;
t52 = pkin(3) * t33 + qJ(4) * t32;
t13 = t33 * t61 - t60;
t31 = qJ(5) + qJ(6);
t23 = sin(t31);
t24 = cos(t31);
t51 = t12 * t24 - t13 * t23;
t50 = t12 * t23 + t13 * t24;
t34 = sin(qJ(5));
t37 = cos(qJ(5));
t49 = t12 * t37 - t13 * t34;
t48 = t12 * t34 + t13 * t37;
t47 = t23 * t33 - t24 * t32;
t46 = t23 * t32 + t24 * t33;
t45 = t32 * t37 - t33 * t34;
t44 = t32 * t34 + t33 * t37;
t42 = g(3) * t47;
t41 = g(3) * t45;
t40 = t71 * t66;
t28 = t39 * pkin(7);
t19 = t39 * t57;
t17 = t36 * t57;
t16 = t53 * t35;
t15 = t36 * t32 + t38 * t59;
t11 = t54 * t38 + t63;
t9 = t70 * t33;
t8 = t70 * t32;
t7 = g(1) * t13 - g(2) * t15;
t6 = t14 * t34 + t15 * t37;
t5 = t14 * t37 - t15 * t34;
t4 = t14 * t23 + t15 * t24;
t3 = t14 * t24 - t15 * t23;
t2 = g(1) * t4 + g(2) * t50 + t46 * t63;
t1 = -g(1) * t3 - g(2) * t51 + t35 * t42;
t10 = [0, t53, t54, 0, 0, 0, 0, 0, t53 * t38, -t16, t7, -t55, t16, -g(1) * t28 - g(2) * t56 - t40, t7, t16, t55, -g(1) * (-t13 * pkin(3) - t12 * qJ(4) + t28) - g(2) * (t15 * pkin(3) + t14 * qJ(4) + t56) - t40, 0, 0, 0, 0, 0, g(1) * t48 - g(2) * t6, g(1) * t49 - g(2) * t5, 0, 0, 0, 0, 0, g(1) * t50 - g(2) * t4, g(1) * t51 - g(2) * t3; 0, 0, 0, 0, 0, 0, 0, 0, t70, t11, t9, -t8, -t11, -g(1) * (-t39 * t67 + t19) - g(2) * (-t36 * t67 + t17) - g(3) * t58, t9, -t11, t8, -g(1) * t19 - g(2) * t17 - g(3) * (t52 * t38 + t58) + (pkin(2) + t52) * t69, 0, 0, 0, 0, 0, t70 * t44, -t38 * t41 + t45 * t69, 0, 0, 0, 0, 0, t70 * t46, t38 * t42 - t47 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, 0, 0, 0, -t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t12 - t32 * t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 - g(2) * t49 - t35 * t41, g(1) * t6 + g(2) * t48 + t44 * t63, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t10;
