% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRR5
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
% Datum: 2019-03-10 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRRR5_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t37 = sin(qJ(2));
t38 = sin(qJ(1));
t41 = cos(qJ(2));
t42 = cos(qJ(1));
t52 = cos(pkin(6));
t48 = t42 * t52;
t21 = t37 * t48 + t38 * t41;
t33 = qJ(3) + qJ(4);
t32 = qJ(5) + t33;
t28 = sin(t32);
t29 = cos(t32);
t34 = sin(pkin(6));
t55 = t34 * t42;
t10 = t21 * t29 - t28 * t55;
t20 = t38 * t37 - t41 * t48;
t35 = sin(qJ(6));
t39 = cos(qJ(6));
t65 = t10 * t35 - t20 * t39;
t64 = t10 * t39 + t20 * t35;
t63 = g(3) * t34;
t60 = t29 * t35;
t59 = t29 * t39;
t58 = t34 * t37;
t57 = t34 * t38;
t40 = cos(qJ(3));
t56 = t34 * t40;
t54 = t35 * t41;
t53 = t39 * t41;
t30 = sin(t33);
t31 = cos(t33);
t51 = t21 * t31 - t30 * t55;
t36 = sin(qJ(3));
t50 = t21 * t40 - t36 * t55;
t49 = t38 * t52;
t47 = t21 * t28 + t29 * t55;
t46 = t21 * t30 + t31 * t55;
t45 = t21 * t36 + t40 * t55;
t23 = -t37 * t49 + t42 * t41;
t12 = -t23 * t28 + t29 * t57;
t44 = g(1) * t12 - g(2) * t47 + g(3) * (-t28 * t58 + t52 * t29);
t22 = t42 * t37 + t41 * t49;
t43 = -g(1) * t22 - g(2) * t20 + t41 * t63;
t19 = t52 * t28 + t29 * t58;
t17 = t23 * t40 + t36 * t57;
t16 = -t23 * t36 + t38 * t56;
t15 = t23 * t31 + t30 * t57;
t14 = -t23 * t30 + t31 * t57;
t13 = t23 * t29 + t28 * t57;
t8 = t13 * t39 + t22 * t35;
t7 = -t13 * t35 + t22 * t39;
t6 = g(1) * t15 + g(2) * t51 - g(3) * (-t52 * t30 - t31 * t58);
t5 = -g(1) * t14 + g(2) * t46 - g(3) * (-t30 * t58 + t52 * t31);
t4 = g(1) * t13 + g(2) * t10 + g(3) * t19;
t2 = t44 * t39;
t1 = t44 * t35;
t3 = [0, g(1) * t38 - g(2) * t42, g(1) * t42 + g(2) * t38, 0, 0, 0, 0, 0, g(1) * t21 - g(2) * t23, -g(1) * t20 + g(2) * t22, 0, 0, 0, 0, 0, g(1) * t50 - g(2) * t17, -g(1) * t45 - g(2) * t16, 0, 0, 0, 0, 0, g(1) * t51 - g(2) * t15, -g(1) * t46 - g(2) * t14, 0, 0, 0, 0, 0, g(1) * t10 - g(2) * t13, -g(1) * t47 - g(2) * t12, 0, 0, 0, 0, 0, g(1) * t64 - g(2) * t8, -g(1) * t65 - g(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, -t43, g(1) * t23 + g(2) * t21 + g(3) * t58, 0, 0, 0, 0, 0, -t43 * t40, t43 * t36, 0, 0, 0, 0, 0, -t43 * t31, t43 * t30, 0, 0, 0, 0, 0, -t43 * t29, t43 * t28, 0, 0, 0, 0, 0, -g(1) * (-t22 * t59 + t23 * t35) - g(2) * (-t20 * t59 + t21 * t35) - (t29 * t53 + t35 * t37) * t63, -g(1) * (t22 * t60 + t23 * t39) - g(2) * (t20 * t60 + t21 * t39) - (-t29 * t54 + t37 * t39) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t16 + g(2) * t45 - g(3) * (-t36 * t58 + t52 * t40) g(1) * t17 + g(2) * t50 - g(3) * (-t52 * t36 - t37 * t56) 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, -t44, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, -t44, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t4, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 + g(2) * t65 - g(3) * (-t19 * t35 - t34 * t53) g(1) * t8 + g(2) * t64 - g(3) * (-t19 * t39 + t34 * t54);];
taug_reg  = t3;
