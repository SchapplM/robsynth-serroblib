% Calculate minimal parameter regressor of gravitation load for
% S6RRPRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [6x29]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRPP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t43 = cos(qJ(1));
t41 = sin(qJ(1));
t74 = g(2) * t41;
t53 = g(1) * t43 + t74;
t40 = sin(qJ(2));
t42 = cos(qJ(2));
t16 = -g(3) * t42 + t53 * t40;
t36 = pkin(9) + qJ(4);
t31 = cos(t36);
t77 = pkin(4) * t31;
t76 = g(1) * t41;
t39 = -pkin(8) - qJ(3);
t72 = pkin(5) - t39;
t30 = sin(t36);
t71 = t30 * t40;
t70 = t31 * t40;
t69 = t40 * t39;
t37 = sin(pkin(9));
t68 = t41 * t37;
t67 = t41 * t42;
t38 = cos(pkin(9));
t29 = t38 * pkin(3) + pkin(2);
t22 = t42 * t29;
t66 = t43 * t30;
t65 = t43 * t31;
t64 = t43 * t37;
t63 = t43 * t38;
t62 = -pkin(4) - qJ(6);
t61 = t43 * pkin(1) + t41 * pkin(7);
t60 = qJ(5) * t30;
t59 = t72 * t43;
t58 = -pkin(1) - t22;
t12 = t30 * t67 + t65;
t13 = t31 * t67 - t66;
t57 = -t12 * pkin(4) + t13 * qJ(5);
t14 = -t41 * t31 + t42 * t66;
t15 = t41 * t30 + t42 * t65;
t56 = -t14 * pkin(4) + t15 * qJ(5);
t55 = -t29 - t60;
t54 = g(3) * (t22 + (t60 + t77) * t42);
t4 = g(1) * t12 - g(2) * t14;
t5 = g(1) * t13 - g(2) * t15;
t52 = -g(2) * t43 + t76;
t51 = t42 * pkin(2) + t40 * qJ(3);
t48 = t53 * t42;
t33 = t43 * pkin(7);
t47 = pkin(3) * t64 - t13 * pkin(4) - t12 * qJ(5) + t41 * t69 + t33;
t46 = pkin(3) * t68 + t15 * pkin(4) + t14 * qJ(5) + t43 * t22 + t61;
t2 = g(1) * t14 + g(2) * t12 + g(3) * t71;
t44 = g(1) * t15 + g(2) * t13 + g(3) * t70;
t20 = qJ(5) * t70;
t18 = t52 * t40;
t17 = g(3) * t40 + t48;
t7 = t16 * t31;
t6 = t16 * t30;
t1 = [0, t52, t53, 0, 0, 0, 0, 0, t52 * t42, -t18, -g(1) * (-t38 * t67 + t64) - g(2) * (t42 * t63 + t68) -g(1) * (t37 * t67 + t63) - g(2) * (t41 * t38 - t42 * t64) t18, -g(1) * t33 - g(2) * (t51 * t43 + t61) - (-pkin(1) - t51) * t76, 0, 0, 0, 0, 0, t5, -t4, t18, -t5, t4, -g(1) * (t58 * t41 + t47) - g(2) * (-t43 * t69 + t46) t18, t4, t5, -g(1) * (-t13 * qJ(6) + t47) - g(2) * (t15 * qJ(6) + t40 * t59 + t46) - (-t40 * pkin(5) + t58) * t76; 0, 0, 0, 0, 0, 0, 0, 0, t16, t17, t16 * t38, -t16 * t37, -t17, -g(3) * t51 + t53 * (pkin(2) * t40 - qJ(3) * t42) 0, 0, 0, 0, 0, t7, -t6, -t17, -t7, t6, -t54 + t39 * t48 + (g(3) * t39 + t53 * (-t55 + t77)) * t40, -t17, t6, t7, -t54 + (-g(3) * qJ(6) * t31 - g(1) * t59 - t72 * t74) * t42 + (-g(3) * t72 + t53 * (-t62 * t31 - t55)) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t44, 0, -t2, -t44, -g(1) * t56 - g(2) * t57 - g(3) * (-pkin(4) * t71 + t20) 0, -t44, t2, -g(1) * (-t14 * qJ(6) + t56) - g(2) * (-t12 * qJ(6) + t57) - g(3) * (t62 * t71 + t20); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44;];
taug_reg  = t1;
