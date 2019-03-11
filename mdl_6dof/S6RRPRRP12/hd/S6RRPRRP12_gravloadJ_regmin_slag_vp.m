% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5]';
% 
% Output:
% taug_reg [6x32]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP12_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRRP12_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t38 = sin(qJ(2));
t41 = cos(qJ(2));
t76 = t41 * pkin(2) + t38 * qJ(3);
t47 = -pkin(1) - t76;
t39 = sin(qJ(1));
t42 = cos(qJ(1));
t22 = g(1) * t42 + g(2) * t39;
t40 = cos(qJ(4));
t55 = t42 * t40;
t37 = sin(qJ(4));
t61 = t39 * t37;
t15 = t38 * t55 - t61;
t56 = t42 * t37;
t60 = t39 * t40;
t17 = t38 * t60 + t56;
t75 = -g(1) * t15 - g(2) * t17;
t14 = g(3) * t38 + t22 * t41;
t74 = pkin(2) * t38;
t73 = pkin(4) * t37;
t71 = g(1) * t39;
t66 = g(3) * t41;
t65 = t40 * pkin(4);
t36 = qJ(4) + qJ(5);
t28 = sin(t36);
t63 = t39 * t28;
t29 = cos(t36);
t62 = t39 * t29;
t43 = -pkin(9) - pkin(8);
t59 = t41 * t43;
t58 = t42 * t28;
t57 = t42 * t29;
t54 = qJ(3) * t41;
t53 = t38 * t56;
t52 = g(3) * t76;
t51 = t39 * pkin(7) - t42 * t47;
t11 = t38 * t62 + t58;
t9 = -t38 * t57 + t63;
t50 = g(1) * t11 + g(2) * t9;
t49 = -g(2) * t42 + t71;
t48 = -pkin(5) * t29 - qJ(6) * t28;
t46 = pkin(5) * t28 - qJ(6) * t29 + t73;
t1 = g(1) * t9 - g(2) * t11 + t29 * t66;
t10 = t38 * t58 + t62;
t12 = -t38 * t63 + t57;
t2 = -g(1) * t10 + g(2) * t12 + t28 * t66;
t44 = -g(1) * (-t9 * pkin(5) + t10 * qJ(6)) - g(2) * (t11 * pkin(5) - t12 * qJ(6));
t33 = t42 * pkin(7);
t27 = pkin(3) + t65;
t25 = t42 * t54;
t23 = t39 * t54;
t20 = t49 * t41;
t19 = t49 * t38;
t18 = -t38 * t61 + t55;
t16 = t53 + t60;
t13 = t22 * t38 - t66;
t6 = t14 * t29;
t5 = t14 * t28;
t4 = -g(1) * t12 - g(2) * t10;
t3 = [0, t49, t22, 0, 0, 0, 0, 0, t20, -t19, -t22, -t20, t19, -g(1) * t33 - g(2) * t51 - t47 * t71, 0, 0, 0, 0, 0, -g(1) * t18 - g(2) * t16, g(1) * t17 - g(2) * t15, 0, 0, 0, 0, 0, t4, t50, t4, t20, -t50, -g(1) * (t12 * pkin(5) + t11 * qJ(6) + t42 * t27 + t33) - g(2) * (pkin(4) * t53 + t10 * pkin(5) + t9 * qJ(6) - t42 * t59 + t51) + (-g(1) * (-t38 * t73 + t47 + t59) - g(2) * t27) * t39; 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, -t13, -t14, -g(1) * (-t42 * t74 + t25) - g(2) * (-t39 * t74 + t23) - t52, 0, 0, 0, 0, 0, -t14 * t37, -t14 * t40, 0, 0, 0, 0, 0, -t5, -t6, -t5, t13, t6, -g(1) * t25 - g(2) * t23 - t52 + (g(3) * t43 - t22 * t46) * t41 + (-g(3) * t46 + t22 * (pkin(2) - t43)) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40 * t66 + t75, g(1) * t16 - g(2) * t18 - t37 * t66, 0, 0, 0, 0, 0, t1, -t2, t1, 0, t2, t75 * pkin(4) - (t48 - t65) * t66 + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, t1, 0, t2, -t48 * t66 + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t3;
