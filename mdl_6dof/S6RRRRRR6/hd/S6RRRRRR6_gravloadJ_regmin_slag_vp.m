% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5,d6]';
% 
% Output:
% taug_reg [6x38]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 04:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR6_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR6_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR6_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR6_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t38 = sin(qJ(2));
t39 = sin(qJ(1));
t42 = cos(qJ(2));
t51 = cos(pkin(6));
t67 = cos(qJ(1));
t47 = t51 * t67;
t23 = t38 * t47 + t39 * t42;
t34 = qJ(3) + qJ(4);
t30 = sin(t34);
t32 = cos(t34);
t35 = sin(pkin(6));
t50 = t35 * t67;
t14 = t23 * t32 - t30 * t50;
t22 = t39 * t38 - t42 * t47;
t33 = qJ(5) + qJ(6);
t29 = sin(t33);
t31 = cos(t33);
t72 = t14 * t29 - t22 * t31;
t71 = t14 * t31 + t22 * t29;
t36 = sin(qJ(5));
t40 = cos(qJ(5));
t70 = t14 * t36 - t22 * t40;
t69 = t14 * t40 + t22 * t36;
t68 = g(3) * t35;
t62 = t29 * t32;
t61 = t31 * t32;
t60 = t32 * t36;
t59 = t32 * t40;
t58 = t32 * t42;
t57 = t35 * t38;
t56 = t35 * t39;
t41 = cos(qJ(3));
t55 = t35 * t41;
t54 = t35 * t42;
t53 = t36 * t42;
t52 = t40 * t42;
t37 = sin(qJ(3));
t49 = t23 * t41 - t37 * t50;
t48 = t39 * t51;
t25 = -t38 * t48 + t42 * t67;
t16 = -t25 * t30 + t32 * t56;
t45 = t23 * t30 + t32 * t50;
t46 = g(1) * t16 - g(2) * t45 + g(3) * (-t30 * t57 + t32 * t51);
t44 = t23 * t37 + t41 * t50;
t24 = t38 * t67 + t42 * t48;
t43 = -g(1) * t24 - g(2) * t22 + g(3) * t54;
t21 = t30 * t51 + t32 * t57;
t19 = t25 * t41 + t37 * t56;
t18 = -t25 * t37 + t39 * t55;
t17 = t25 * t32 + t30 * t56;
t12 = t17 * t40 + t24 * t36;
t11 = -t17 * t36 + t24 * t40;
t10 = t17 * t31 + t24 * t29;
t9 = -t17 * t29 + t24 * t31;
t8 = g(1) * t17 + g(2) * t14 + g(3) * t21;
t6 = t46 * t40;
t5 = t46 * t36;
t4 = t46 * t31;
t3 = t46 * t29;
t2 = g(1) * t10 + g(2) * t71 - g(3) * (-t21 * t31 + t29 * t54);
t1 = -g(1) * t9 + g(2) * t72 - g(3) * (-t21 * t29 - t31 * t54);
t7 = [0, g(1) * t39 - g(2) * t67, g(1) * t67 + g(2) * t39, 0, 0, 0, 0, 0, g(1) * t23 - g(2) * t25, -g(1) * t22 + g(2) * t24, 0, 0, 0, 0, 0, g(1) * t49 - g(2) * t19, -g(1) * t44 - g(2) * t18, 0, 0, 0, 0, 0, g(1) * t14 - g(2) * t17, -g(1) * t45 - g(2) * t16, 0, 0, 0, 0, 0, g(1) * t69 - g(2) * t12, -g(1) * t70 - g(2) * t11, 0, 0, 0, 0, 0, g(1) * t71 - g(2) * t10, -g(1) * t72 - g(2) * t9; 0, 0, 0, 0, 0, 0, 0, 0, -t43, g(1) * t25 + g(2) * t23 + g(3) * t57, 0, 0, 0, 0, 0, -t43 * t41, t43 * t37, 0, 0, 0, 0, 0, -t43 * t32, t43 * t30, 0, 0, 0, 0, 0, -g(1) * (-t24 * t59 + t25 * t36) - g(2) * (-t22 * t59 + t23 * t36) - (t32 * t52 + t36 * t38) * t68, -g(1) * (t24 * t60 + t25 * t40) - g(2) * (t22 * t60 + t23 * t40) - (-t32 * t53 + t38 * t40) * t68, 0, 0, 0, 0, 0, -g(1) * (-t24 * t61 + t25 * t29) - g(2) * (-t22 * t61 + t23 * t29) - (t29 * t38 + t31 * t58) * t68, -g(1) * (t24 * t62 + t25 * t31) - g(2) * (t22 * t62 + t23 * t31) - (-t29 * t58 + t31 * t38) * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t18 + g(2) * t44 - g(3) * (-t37 * t57 + t41 * t51) g(1) * t19 + g(2) * t49 - g(3) * (-t37 * t51 - t38 * t55) 0, 0, 0, 0, 0, -t46, t8, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, t8, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 + g(2) * t70 - g(3) * (-t21 * t36 - t35 * t52) g(1) * t12 + g(2) * t69 - g(3) * (-t21 * t40 + t35 * t53) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t7;
