% Calculate minimal parameter regressor of gravitation load for
% S6PRPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
% 
% Output:
% taug_reg [6x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:33
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRPR2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_gravloadJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t33 = sin(pkin(11));
t39 = sin(qJ(2));
t41 = cos(qJ(2));
t58 = cos(pkin(11));
t20 = -t41 * t33 - t39 * t58;
t34 = sin(pkin(10));
t36 = cos(pkin(10));
t37 = cos(pkin(6));
t61 = t37 * t41;
t75 = -t34 * t61 - t36 * t39;
t47 = -t39 * t33 + t41 * t58;
t43 = t47 * t37;
t10 = t36 * t20 - t34 * t43;
t57 = sin(pkin(6));
t49 = t58 * t57;
t52 = t39 * t57;
t17 = t33 * t52 - t41 * t49;
t7 = t34 * t20 + t36 * t43;
t74 = -g(1) * t10 - g(2) * t7 + g(3) * t17;
t59 = t20 * t37;
t11 = t34 * t59 + t36 * t47;
t6 = -t34 * t47 + t36 * t59;
t31 = pkin(12) + qJ(6);
t29 = sin(t31);
t40 = cos(qJ(4));
t70 = t29 * t40;
t30 = cos(t31);
t69 = t30 * t40;
t32 = sin(pkin(12));
t68 = t32 * t40;
t66 = t34 * t39;
t35 = cos(pkin(12));
t65 = t35 * t40;
t62 = t37 * t39;
t55 = t36 * t61;
t38 = sin(qJ(4));
t54 = t38 * t57;
t51 = t40 * t57;
t50 = t41 * t57;
t18 = t33 * t50 + t39 * t49;
t12 = t18 * t38 - t37 * t40;
t2 = t36 * t51 - t38 * t6;
t4 = t11 * t38 - t34 * t51;
t46 = g(1) * t4 + g(2) * t2 + g(3) * t12;
t13 = t18 * t40 + t37 * t38;
t3 = -t36 * t54 - t40 * t6;
t5 = t11 * t40 + t34 * t54;
t45 = g(1) * t5 + g(2) * t3 + g(3) * t13;
t42 = -g(1) * t75 - g(3) * t50;
t21 = pkin(2) * t55;
t16 = -g(3) * t37 + (-g(1) * t34 + g(2) * t36) * t57;
t1 = t74 * t38;
t8 = [-g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, -g(2) * (t55 - t66) + t42, -g(1) * (t34 * t62 - t36 * t41) - g(2) * (-t34 * t41 - t36 * t62) + g(3) * t52, -g(2) * t21 + (g(2) * t66 + t42) * pkin(2), 0, 0, 0, 0, 0, t74 * t40, -t1, -g(1) * (t10 * t65 + t11 * t32) - g(2) * (-t6 * t32 + t7 * t65) - g(3) * (-t17 * t65 + t18 * t32) -g(1) * (-t10 * t68 + t11 * t35) - g(2) * (-t6 * t35 - t7 * t68) - g(3) * (t17 * t68 + t18 * t35) t1, -g(1) * (t75 * pkin(2) + pkin(8) * t11) - g(2) * (-pkin(2) * t66 - t6 * pkin(8) + t21) - g(3) * (pkin(2) * t50 + t18 * pkin(8)) + t74 * (pkin(4) * t40 + qJ(5) * t38 + pkin(3)) 0, 0, 0, 0, 0, -g(1) * (t10 * t69 + t11 * t29) - g(2) * (-t6 * t29 + t7 * t69) - g(3) * (-t17 * t69 + t18 * t29) -g(1) * (-t10 * t70 + t11 * t30) - g(2) * (-t6 * t30 - t7 * t70) - g(3) * (t17 * t70 + t18 * t30); 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t45, t46 * t35, -t46 * t32, -t45, -g(1) * (-t4 * pkin(4) + t5 * qJ(5)) - g(2) * (-t2 * pkin(4) + t3 * qJ(5)) - g(3) * (-t12 * pkin(4) + t13 * qJ(5)) 0, 0, 0, 0, 0, t46 * t30, -t46 * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t10 * t30 - t5 * t29) - g(2) * (-t3 * t29 - t7 * t30) - g(3) * (-t13 * t29 + t17 * t30) -g(1) * (t10 * t29 - t5 * t30) - g(2) * (t7 * t29 - t3 * t30) - g(3) * (-t13 * t30 - t17 * t29);];
taug_reg  = t8;
