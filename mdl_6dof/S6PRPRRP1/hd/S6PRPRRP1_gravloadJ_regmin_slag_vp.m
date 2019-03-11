% Calculate minimal parameter regressor of gravitation load for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [6x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRP1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t29 = sin(pkin(11));
t37 = sin(qJ(2));
t40 = cos(qJ(2));
t52 = cos(pkin(11));
t47 = -t37 * t29 + t40 * t52;
t30 = sin(pkin(10));
t32 = cos(pkin(10));
t33 = cos(pkin(6));
t56 = t33 * t40;
t70 = -t30 * t56 - t32 * t37;
t21 = -t40 * t29 - t37 * t52;
t43 = t47 * t33;
t10 = t32 * t21 - t30 * t43;
t31 = sin(pkin(6));
t17 = t47 * t31;
t7 = t30 * t21 + t32 * t43;
t69 = -g(1) * t10 - g(2) * t7 - g(3) * t17;
t35 = sin(qJ(5));
t19 = t21 * t33;
t8 = -t32 * t19 + t30 * t47;
t65 = t8 * t35;
t9 = -t30 * t19 - t32 * t47;
t64 = t9 * t35;
t18 = t21 * t31;
t63 = t18 * t35;
t62 = t30 * t37;
t36 = sin(qJ(4));
t61 = t31 * t36;
t39 = cos(qJ(4));
t60 = t31 * t39;
t59 = t31 * t40;
t57 = t33 * t37;
t55 = t35 * t39;
t38 = cos(qJ(5));
t53 = t38 * t39;
t50 = t32 * t56;
t12 = -t18 * t36 - t33 * t39;
t2 = t32 * t60 + t8 * t36;
t4 = -t30 * t60 - t36 * t9;
t46 = g(1) * t4 + g(2) * t2 + g(3) * t12;
t13 = -t18 * t39 + t33 * t36;
t3 = -t32 * t61 + t8 * t39;
t5 = t30 * t61 - t39 * t9;
t45 = g(1) * t5 + g(2) * t3 + g(3) * t13;
t42 = -g(1) * t70 - g(3) * t59;
t41 = -g(1) * (-t10 * t38 - t5 * t35) - g(2) * (-t3 * t35 - t7 * t38) - g(3) * (-t13 * t35 - t17 * t38);
t34 = -qJ(6) - pkin(9);
t28 = t38 * pkin(5) + pkin(4);
t22 = pkin(2) * t50;
t16 = -g(3) * t33 + (-g(1) * t30 + g(2) * t32) * t31;
t1 = t69 * t36;
t6 = [-g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -g(2) * (t50 - t62) + t42, -g(1) * (t30 * t57 - t32 * t40) - g(2) * (-t30 * t40 - t32 * t57) + g(3) * t31 * t37, -g(2) * t22 + (g(2) * t62 + t42) * pkin(2), 0, 0, 0, 0, 0, t69 * t39, -t1, 0, 0, 0, 0, 0, -g(1) * (t10 * t53 - t64) - g(2) * (t7 * t53 + t65) - g(3) * (t17 * t53 - t63) -g(1) * (-t10 * t55 - t9 * t38) - g(2) * (t38 * t8 - t7 * t55) - g(3) * (-t17 * t55 - t18 * t38) t1, -g(1) * (t70 * pkin(2) - pkin(5) * t64 - t9 * pkin(8)) - g(2) * (-pkin(2) * t62 + pkin(5) * t65 + pkin(8) * t8 + t22) - g(3) * (pkin(2) * t59 - pkin(5) * t63 - t18 * pkin(8)) + t69 * (t28 * t39 - t34 * t36 + pkin(3)); 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, t45, 0, 0, 0, 0, 0, t46 * t38, -t46 * t35, -t45, -g(1) * (-t4 * t28 - t5 * t34) - g(2) * (-t2 * t28 - t3 * t34) - g(3) * (-t12 * t28 - t13 * t34); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -g(1) * (t10 * t35 - t5 * t38) - g(2) * (-t3 * t38 + t7 * t35) - g(3) * (-t13 * t38 + t17 * t35) 0, t41 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46;];
taug_reg  = t6;
