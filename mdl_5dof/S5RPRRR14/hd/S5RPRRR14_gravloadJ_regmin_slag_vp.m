% Calculate minimal parameter regressor of gravitation load for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR14_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR14_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t35 = cos(pkin(5));
t39 = sin(qJ(1));
t56 = cos(pkin(11));
t52 = t39 * t56;
t32 = sin(pkin(11));
t42 = cos(qJ(1));
t59 = t42 * t32;
t25 = t35 * t59 + t52;
t38 = sin(qJ(3));
t33 = sin(pkin(6));
t34 = sin(pkin(5));
t66 = cos(qJ(3));
t55 = t34 * t66;
t50 = t33 * t55;
t51 = t42 * t56;
t61 = t39 * t32;
t24 = -t35 * t51 + t61;
t57 = cos(pkin(6));
t54 = t24 * t57;
t10 = t25 * t38 + t42 * t50 + t54 * t66;
t36 = sin(qJ(5));
t63 = t34 * t38;
t11 = -t42 * t33 * t63 + t25 * t66 - t38 * t54;
t53 = t34 * t57;
t19 = t24 * t33 - t42 * t53;
t37 = sin(qJ(4));
t41 = cos(qJ(4));
t4 = t11 * t41 + t19 * t37;
t40 = cos(qJ(5));
t71 = -t10 * t40 + t36 * t4;
t70 = t10 * t36 + t4 * t40;
t67 = t11 * t37 - t19 * t41;
t65 = t33 * t35;
t64 = t34 * t33;
t62 = t36 * t41;
t60 = t40 * t41;
t58 = qJ(2) * t34;
t49 = g(1) * t42 + g(2) * t39;
t48 = g(1) * t39 - g(2) * t42;
t47 = t57 * t56;
t17 = t38 * t65 + (t32 * t66 + t38 * t47) * t34;
t23 = t35 * t57 - t56 * t64;
t26 = -t35 * t61 + t51;
t44 = t35 * t52 + t59;
t43 = t44 * t57;
t15 = t26 * t66 + (t39 * t64 - t43) * t38;
t20 = t33 * t44 + t39 * t53;
t6 = -t15 * t37 + t20 * t41;
t46 = g(1) * t6 - g(2) * t67 + g(3) * (-t17 * t37 + t23 * t41);
t14 = t26 * t38 - t39 * t50 + t43 * t66;
t16 = t32 * t63 - t47 * t55 - t65 * t66;
t45 = g(1) * t14 + g(2) * t10 + g(3) * t16;
t9 = t17 * t41 + t23 * t37;
t7 = t15 * t41 + t20 * t37;
t2 = t14 * t36 + t40 * t7;
t1 = t14 * t40 - t36 * t7;
t3 = [0, t48, t49, g(1) * t25 - g(2) * t26, -g(1) * t24 + g(2) * t44, -t49 * t34, -g(1) * (-pkin(1) * t39 + t42 * t58) - g(2) * (pkin(1) * t42 + t39 * t58), 0, 0, 0, 0, 0, g(1) * t11 - g(2) * t15, -g(1) * t10 + g(2) * t14, 0, 0, 0, 0, 0, g(1) * t4 - g(2) * t7, -g(1) * t67 - g(2) * t6, 0, 0, 0, 0, 0, g(1) * t70 - g(2) * t2, -g(1) * t71 - g(2) * t1; 0, 0, 0, 0, 0, 0, -g(3) * t35 - t34 * t48, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, g(1) * t15 + g(2) * t11 + g(3) * t17, 0, 0, 0, 0, 0, t45 * t41, -t45 * t37, 0, 0, 0, 0, 0, -g(1) * (-t14 * t60 + t15 * t36) - g(2) * (-t10 * t60 + t11 * t36) - g(3) * (-t16 * t60 + t17 * t36), -g(1) * (t14 * t62 + t15 * t40) - g(2) * (t10 * t62 + t11 * t40) - g(3) * (t16 * t62 + t17 * t40); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, g(1) * t7 + g(2) * t4 + g(3) * t9, 0, 0, 0, 0, 0, -t46 * t40, t46 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t71 - g(3) * (t16 * t40 - t36 * t9), g(1) * t2 + g(2) * t70 - g(3) * (-t16 * t36 - t40 * t9);];
taug_reg = t3;
