% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:19
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR4_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t39 = cos(qJ(2));
t34 = qJ(2) + qJ(3);
t30 = sin(t34);
t31 = cos(t34);
t48 = pkin(3) * t31 + qJ(4) * t30;
t64 = pkin(2) * t39 + t48;
t38 = sin(qJ(1));
t40 = cos(qJ(1));
t45 = g(1) * t40 + g(2) * t38;
t9 = -g(3) * t31 + t30 * t45;
t63 = pkin(3) * t30;
t62 = g(3) * t30;
t33 = pkin(11) + qJ(5);
t29 = qJ(6) + t33;
t23 = sin(t29);
t60 = t38 * t23;
t24 = cos(t29);
t59 = t38 * t24;
t27 = sin(t33);
t58 = t38 * t27;
t28 = cos(t33);
t57 = t38 * t28;
t35 = sin(pkin(11));
t56 = t38 * t35;
t36 = cos(pkin(11));
t55 = t38 * t36;
t54 = t40 * t23;
t53 = t40 * t24;
t52 = t40 * t27;
t51 = t40 * t28;
t50 = t40 * t35;
t49 = t40 * t36;
t47 = qJ(4) * t31;
t37 = sin(qJ(2));
t46 = -pkin(2) * t37 - t63;
t44 = g(1) * t38 - g(2) * t40;
t43 = pkin(1) + t64;
t41 = -pkin(8) - pkin(7);
t21 = t40 * t47;
t20 = t38 * t47;
t19 = t44 * t30;
t18 = t31 * t51 + t58;
t17 = -t31 * t52 + t57;
t16 = -t31 * t57 + t52;
t15 = t31 * t58 + t51;
t14 = t31 * t53 + t60;
t13 = -t31 * t54 + t59;
t12 = -t31 * t59 + t54;
t11 = t31 * t60 + t53;
t10 = t31 * t45 + t62;
t8 = t9 * t36;
t7 = t9 * t35;
t6 = t9 * t28;
t5 = t9 * t27;
t4 = t9 * t24;
t3 = t9 * t23;
t2 = g(1) * t14 - g(2) * t12 + t24 * t62;
t1 = -g(1) * t13 + g(2) * t11 + t23 * t62;
t22 = [0, t44, t45, 0, 0, 0, 0, 0, t44 * t39, -t44 * t37, 0, 0, 0, 0, 0, t44 * t31, -t19, -g(1) * (-t31 * t55 + t50) - g(2) * (t31 * t49 + t56) -g(1) * (t31 * t56 + t49) - g(2) * (-t31 * t50 + t55) t19 (g(1) * t41 - g(2) * t43) * t40 + (g(1) * t43 + g(2) * t41) * t38, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, -g(1) * t15 - g(2) * t17, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t14, -g(1) * t11 - g(2) * t13; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t39 + t37 * t45, g(3) * t37 + t39 * t45, 0, 0, 0, 0, 0, t9, t10, t8, -t7, -t10, -g(1) * (t40 * t46 + t21) - g(2) * (t38 * t46 + t20) - g(3) * t64, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, t8, -t7, -t10, -g(1) * (-t40 * t63 + t21) - g(2) * (-t38 * t63 + t20) - g(3) * t48, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t17 + g(2) * t15 + t27 * t62, g(1) * t18 - g(2) * t16 + t28 * t62, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t22;
