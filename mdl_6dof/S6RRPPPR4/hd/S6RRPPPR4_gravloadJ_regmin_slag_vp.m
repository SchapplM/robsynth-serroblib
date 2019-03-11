% Calculate minimal parameter regressor of gravitation load for
% S6RRPPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t34 = sin(qJ(2));
t35 = sin(qJ(1));
t38 = cos(qJ(1));
t14 = g(1) * t38 + g(2) * t35;
t37 = cos(qJ(2));
t69 = t14 * t37;
t7 = g(3) * t34 + t69;
t70 = t14 * t34;
t22 = t34 * qJ(3);
t26 = t37 * pkin(2);
t55 = t26 + t22;
t53 = qJ(3) * t37;
t17 = t35 * t53;
t20 = t38 * t53;
t23 = t37 * qJ(4);
t56 = pkin(2) + qJ(4);
t67 = t56 * t70 - g(1) * t20 - g(2) * t17 - g(3) * (t23 + t55);
t66 = pkin(2) * t34;
t65 = g(1) * t35;
t61 = g(3) * t37;
t31 = sin(pkin(9));
t60 = t35 * t31;
t32 = cos(pkin(9));
t59 = t35 * t32;
t58 = t38 * t31;
t57 = t38 * t32;
t27 = t38 * pkin(7);
t54 = t38 * pkin(3) + t27;
t52 = -pkin(1) - t22;
t51 = t35 * pkin(7) + (pkin(1) + t55) * t38;
t10 = t34 * t59 + t58;
t8 = -t34 * t57 + t60;
t50 = g(1) * t10 + g(2) * t8;
t49 = -g(2) * t38 + t65;
t11 = -t34 * t60 + t57;
t33 = sin(qJ(6));
t36 = cos(qJ(6));
t47 = t10 * t36 - t11 * t33;
t46 = t10 * t33 + t11 * t36;
t45 = t31 * t36 - t32 * t33;
t44 = t31 * t33 + t32 * t36;
t43 = t35 * pkin(3) + t38 * t23 + t51;
t42 = g(3) * t45;
t39 = (-t56 * t37 + t52) * t65;
t13 = t49 * t37;
t12 = t49 * t34;
t9 = t34 * t58 + t59;
t6 = -t61 + t70;
t5 = t7 * t32;
t4 = t7 * t31;
t3 = -g(1) * t11 - g(2) * t9;
t2 = t8 * t33 + t9 * t36;
t1 = -t9 * t33 + t8 * t36;
t15 = [0, t49, t14, 0, 0, 0, 0, 0, t13, -t12, -t14, -t13, t12, -g(1) * t27 - g(2) * t51 - (t52 - t26) * t65, t3, t50, t13, -g(1) * t54 - g(2) * t43 - t39, t3, t13, -t50, -g(1) * (t11 * pkin(4) + t10 * qJ(5) + t54) - g(2) * (t9 * pkin(4) + t8 * qJ(5) + t43) - t39, 0, 0, 0, 0, 0, -g(1) * t46 - g(2) * t2, -g(1) * t47 - g(2) * t1; 0, 0, 0, 0, 0, 0, 0, 0, t6, t7, 0, -t6, -t7, -g(1) * (-t38 * t66 + t20) - g(2) * (-t35 * t66 + t17) - g(3) * t55, -t4, -t5, t6, t67, -t4, t6, t5, -t7 * (pkin(4) * t31 - qJ(5) * t32) + t67, 0, 0, 0, 0, 0, -t34 * t42 - t45 * t69, t7 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, -t6, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t8 + g(2) * t10 - t32 * t61, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t47 - t44 * t61, g(1) * t2 - g(2) * t46 - t37 * t42;];
taug_reg  = t15;
