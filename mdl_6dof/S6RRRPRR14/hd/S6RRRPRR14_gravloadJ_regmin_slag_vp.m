% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6]';
% 
% Output:
% taug_reg [6x35]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRR14_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR14_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRR14_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR14_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t33 = sin(qJ(2));
t34 = sin(qJ(1));
t37 = cos(qJ(2));
t50 = cos(pkin(6));
t65 = cos(qJ(1));
t43 = t50 * t65;
t19 = t33 * t43 + t34 * t37;
t32 = sin(qJ(3));
t36 = cos(qJ(3));
t30 = sin(pkin(6));
t49 = t30 * t65;
t10 = t19 * t32 + t36 * t49;
t18 = t34 * t33 - t37 * t43;
t29 = qJ(5) + qJ(6);
t27 = sin(t29);
t28 = cos(t29);
t73 = t10 * t27 + t18 * t28;
t72 = t10 * t28 - t18 * t27;
t31 = sin(qJ(5));
t35 = cos(qJ(5));
t71 = t10 * t31 + t18 * t35;
t70 = t10 * t35 - t18 * t31;
t48 = t34 * t50;
t20 = t33 * t65 + t37 * t48;
t69 = -g(1) * t20 - g(2) * t18;
t66 = g(3) * t30;
t60 = t27 * t32;
t59 = t28 * t32;
t58 = t30 * t33;
t57 = t30 * t34;
t56 = t30 * t36;
t55 = t30 * t37;
t54 = t31 * t32;
t53 = t32 * t35;
t52 = t32 * t37;
t51 = t35 * t37;
t11 = t19 * t36 - t32 * t49;
t21 = -t33 * t48 + t37 * t65;
t14 = t21 * t32 - t34 * t56;
t47 = -g(1) * t10 + g(2) * t14;
t15 = t21 * t36 + t32 * t57;
t46 = -g(1) * t11 + g(2) * t15;
t45 = g(1) * t18 - g(2) * t20;
t44 = -g(1) * t21 - g(2) * t19;
t16 = t32 * t58 - t36 * t50;
t41 = g(1) * t14 + g(2) * t10 + g(3) * t16;
t17 = t32 * t50 + t33 * t56;
t40 = g(1) * t15 + g(2) * t11 + g(3) * t17;
t39 = g(3) * t55 + t69;
t38 = g(3) * t58 - t44;
t9 = t39 * t36;
t8 = t39 * t32;
t7 = t14 * t31 + t20 * t35;
t6 = t14 * t35 - t20 * t31;
t5 = t14 * t27 + t20 * t28;
t4 = t14 * t28 - t20 * t27;
t2 = g(1) * t5 + g(2) * t73 - g(3) * (-t16 * t27 + t28 * t55);
t1 = -g(1) * t4 - g(2) * t72 - g(3) * (t16 * t28 + t27 * t55);
t3 = [0, g(1) * t34 - g(2) * t65, g(1) * t65 + g(2) * t34, 0, 0, 0, 0, 0, g(1) * t19 - g(2) * t21, -t45, 0, 0, 0, 0, 0, -t46, t47, t45, t46, -t47, -g(1) * (-t34 * pkin(1) - t19 * pkin(2) - pkin(3) * t11 + pkin(8) * t49 - t18 * pkin(9) - qJ(4) * t10) - g(2) * (pkin(1) * t65 + t21 * pkin(2) + t15 * pkin(3) + pkin(8) * t57 + t20 * pkin(9) + t14 * qJ(4)) 0, 0, 0, 0, 0, g(1) * t71 - g(2) * t7, g(1) * t70 - g(2) * t6, 0, 0, 0, 0, 0, g(1) * t73 - g(2) * t5, g(1) * t72 - g(2) * t4; 0, 0, 0, 0, 0, 0, 0, 0, -t39, t38, 0, 0, 0, 0, 0, -t9, t8, -t38, t9, -t8 (-t33 * t66 + t44) * pkin(9) + (-t37 * t66 - t69) * (pkin(3) * t36 + qJ(4) * t32 + pkin(2)) 0, 0, 0, 0, 0, -g(1) * (-t20 * t54 + t21 * t35) - g(2) * (-t18 * t54 + t19 * t35) - (t31 * t52 + t33 * t35) * t66, -g(1) * (-t20 * t53 - t21 * t31) - g(2) * (-t18 * t53 - t19 * t31) - (-t31 * t33 + t32 * t51) * t66, 0, 0, 0, 0, 0, -g(1) * (-t20 * t60 + t21 * t28) - g(2) * (-t18 * t60 + t19 * t28) - (t27 * t52 + t28 * t33) * t66, -g(1) * (-t20 * t59 - t21 * t27) - g(2) * (-t18 * t59 - t19 * t27) - (-t27 * t33 + t28 * t52) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t40, 0, -t41, -t40, -g(1) * (-t14 * pkin(3) + t15 * qJ(4)) - g(2) * (-t10 * pkin(3) + t11 * qJ(4)) - g(3) * (-t16 * pkin(3) + t17 * qJ(4)) 0, 0, 0, 0, 0, -t40 * t31, -t40 * t35, 0, 0, 0, 0, 0, -t40 * t27, -t40 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t70 - g(3) * (t16 * t35 + t31 * t55) g(1) * t7 + g(2) * t71 - g(3) * (-t16 * t31 + t30 * t51) 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t3;
