% Calculate minimal parameter regressor of gravitation load for
% S6RRRRRR10V2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
% 
% Output:
% taug_reg [6x38]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRRRR10V2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRRRR10V2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_gravloadJ_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t29 = qJ(2) + qJ(3);
t28 = cos(t29);
t32 = sin(qJ(4));
t39 = cos(qJ(1));
t49 = t39 * t32;
t34 = sin(qJ(1));
t37 = cos(qJ(4));
t52 = t34 * t37;
t22 = t28 * t52 - t49;
t36 = cos(qJ(5));
t27 = sin(t29);
t31 = sin(qJ(5));
t60 = t27 * t31;
t10 = t22 * t36 + t34 * t60;
t48 = t39 * t37;
t53 = t34 * t32;
t21 = t28 * t53 + t48;
t30 = sin(qJ(6));
t35 = cos(qJ(6));
t69 = t10 * t30 - t21 * t35;
t68 = t10 * t35 + t21 * t30;
t67 = -g(1) * t39 - g(2) * t34;
t14 = -g(3) * t28 - t27 * t67;
t64 = g(3) * t27;
t59 = t27 * t36;
t58 = t27 * t39;
t57 = t30 * t32;
t56 = t30 * t36;
t55 = t31 * t37;
t54 = t32 * t35;
t51 = t35 * t36;
t50 = t36 * t37;
t47 = t27 * t57;
t46 = t27 * t54;
t45 = t27 * t49;
t43 = g(1) * t34 - g(2) * t39;
t9 = -t22 * t31 + t34 * t59;
t19 = t27 * t50 - t28 * t31;
t42 = t27 * t55 + t28 * t36;
t24 = t28 * t48 + t53;
t12 = -t24 * t31 + t36 * t58;
t41 = g(1) * t12 + g(2) * t9 - g(3) * t42;
t23 = t28 * t49 - t52;
t40 = g(1) * t23 + g(2) * t21 + t32 * t64;
t38 = cos(qJ(2));
t33 = sin(qJ(2));
t20 = t28 * t50 + t60;
t17 = t19 * t39;
t16 = t19 * t34;
t15 = -t28 * t67 + t64;
t13 = t24 * t36 + t31 * t58;
t8 = t14 * t37;
t7 = t14 * t32;
t6 = t13 * t35 + t23 * t30;
t5 = -t13 * t30 + t23 * t35;
t4 = g(1) * t17 + g(2) * t16 - g(3) * t20;
t3 = -g(3) * (-t28 * t55 + t59) + t67 * t42;
t2 = -g(1) * (-t17 * t35 - t30 * t45) - g(2) * (-t16 * t35 - t34 * t47) - g(3) * (t20 * t35 + t28 * t57);
t1 = -g(1) * (t17 * t30 - t35 * t45) - g(2) * (t16 * t30 - t34 * t46) - g(3) * (-t20 * t30 + t28 * t54);
t11 = [0, t43, -t67, 0, 0, 0, 0, 0, t43 * t38, -t43 * t33, 0, 0, 0, 0, 0, t43 * t28, -t43 * t27, 0, 0, 0, 0, 0, g(1) * t22 - g(2) * t24, -g(1) * t21 + g(2) * t23, 0, 0, 0, 0, 0, g(1) * t10 - g(2) * t13, g(1) * t9 - g(2) * t12, 0, 0, 0, 0, 0, g(1) * t68 - g(2) * t6, -g(1) * t69 - g(2) * t5; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t38 - t33 * t67, g(3) * t33 - t38 * t67, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, t8, -t7, 0, 0, 0, 0, 0, t4, t3, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, t8, -t7, 0, 0, 0, 0, 0, t4, t3, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, g(1) * t24 + g(2) * t22 + t37 * t64, 0, 0, 0, 0, 0, t40 * t36, -t40 * t31, 0, 0, 0, 0, 0, -g(1) * (-t23 * t51 + t24 * t30) - g(2) * (-t21 * t51 + t22 * t30) - (t30 * t37 - t32 * t51) * t64, -g(1) * (t23 * t56 + t24 * t35) - g(2) * (t21 * t56 + t22 * t35) - (t32 * t56 + t35 * t37) * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, g(1) * t13 + g(2) * t10 + g(3) * t19, 0, 0, 0, 0, 0, -t41 * t35, t41 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t69 - g(3) * (-t19 * t30 + t46) g(1) * t6 + g(2) * t68 - g(3) * (-t19 * t35 - t47);];
taug_reg  = t11;
