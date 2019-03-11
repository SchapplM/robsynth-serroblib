% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRR3
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
% Datum: 2019-03-09 07:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR3_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t43 = cos(qJ(4));
t35 = t43 * pkin(4);
t29 = t35 + pkin(3);
t41 = sin(qJ(3));
t44 = cos(qJ(3));
t46 = -pkin(9) - pkin(8);
t50 = t44 * t29 - t41 * t46;
t39 = qJ(4) + qJ(5);
t33 = cos(t39);
t22 = pkin(5) * t33 + t35;
t20 = pkin(3) + t22;
t38 = -pkin(10) + t46;
t52 = t44 * t20 - t41 * t38;
t37 = qJ(1) + pkin(11);
t30 = sin(t37);
t31 = cos(t37);
t19 = g(1) * t31 + g(2) * t30;
t40 = sin(qJ(4));
t66 = t40 * t44;
t14 = t30 * t66 + t31 * t43;
t16 = t30 * t43 - t31 * t66;
t74 = g(3) * t41;
t82 = -g(1) * t16 + g(2) * t14 + t40 * t74;
t32 = sin(t39);
t68 = t32 * t44;
t11 = t30 * t33 - t31 * t68;
t9 = t30 * t68 + t31 * t33;
t3 = -g(1) * t11 + g(2) * t9 + t32 * t74;
t48 = -g(3) * t44 + t19 * t41;
t78 = g(1) * t30;
t72 = t40 * pkin(4);
t71 = t30 * t44;
t70 = t31 * t40;
t69 = t31 * t44;
t67 = t33 * t44;
t63 = t43 * t44;
t45 = cos(qJ(1));
t58 = t45 * pkin(1) + t31 * pkin(2) + t30 * pkin(7);
t42 = sin(qJ(1));
t57 = -t42 * pkin(1) + t31 * pkin(7);
t56 = t44 * pkin(3) + t41 * pkin(8);
t54 = -g(2) * t31 + t78;
t53 = g(1) * t42 - g(2) * t45;
t34 = qJ(6) + t39;
t28 = cos(t34);
t27 = sin(t34);
t21 = pkin(5) * t32 + t72;
t18 = t54 * t41;
t17 = t30 * t40 + t31 * t63;
t15 = -t30 * t63 + t70;
t13 = t19 * t44 + t74;
t12 = t30 * t32 + t31 * t67;
t10 = -t30 * t67 + t31 * t32;
t8 = t30 * t27 + t28 * t69;
t7 = -t27 * t69 + t30 * t28;
t6 = t31 * t27 - t28 * t71;
t5 = t27 * t71 + t31 * t28;
t4 = g(1) * t12 - g(2) * t10 + t33 * t74;
t2 = g(1) * t8 - g(2) * t6 + t28 * t74;
t1 = -g(1) * t7 + g(2) * t5 + t27 * t74;
t23 = [0, 0, 0, 0, 0, 0, t53, g(1) * t45 + g(2) * t42, 0, 0, 0, 0, 0, 0, 0, 0, t54, t19, 0, t53 * pkin(1), 0, 0, 0, 0, 0, 0, t54 * t44, -t18, -t19, -g(1) * (-t30 * pkin(2) + t57) - g(2) * t58, 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t18, -g(1) * t57 - g(2) * (t56 * t31 + t58) - (-pkin(2) - t56) * t78, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t18, -g(1) * (pkin(4) * t70 + t57) - g(2) * (t50 * t31 + t58) + (-g(1) * (-pkin(2) - t50) - g(2) * t72) * t30, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7, t18, -g(1) * (t31 * t21 + t57) - g(2) * (t52 * t31 + t58) + (-g(1) * (-pkin(2) - t52) - g(2) * t21) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t13, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t43, -t48 * t40, -t13, -g(3) * t56 + t19 * (pkin(3) * t41 - pkin(8) * t44) 0, 0, 0, 0, 0, 0, t48 * t33, -t48 * t32, -t13, -g(3) * t50 + t19 * (t29 * t41 + t44 * t46) 0, 0, 0, 0, 0, 0, t48 * t28, -t48 * t27, -t13, -g(3) * t52 + t19 * (t20 * t41 + t38 * t44); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, g(1) * t17 - g(2) * t15 + t43 * t74, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t82 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t21 * t69 + t30 * t22) - g(2) * (-t21 * t71 - t31 * t22) + t21 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t23;
