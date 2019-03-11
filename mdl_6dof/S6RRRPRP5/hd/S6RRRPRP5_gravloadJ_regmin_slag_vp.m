% Calculate minimal parameter regressor of gravitation load for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t41 = sin(qJ(1));
t44 = cos(qJ(1));
t49 = g(1) * t44 + g(2) * t41;
t39 = sin(qJ(3));
t42 = cos(qJ(3));
t55 = t44 * t42;
t43 = cos(qJ(2));
t61 = t41 * t43;
t15 = t39 * t61 + t55;
t56 = t44 * t39;
t17 = t41 * t42 - t43 * t56;
t40 = sin(qJ(2));
t67 = g(3) * t40;
t72 = -g(1) * t17 + g(2) * t15 + t39 * t67;
t13 = -g(3) * t43 + t49 * t40;
t65 = t39 * pkin(3);
t37 = qJ(3) + pkin(10);
t30 = qJ(5) + t37;
t27 = sin(t30);
t64 = t27 * t40;
t28 = cos(t30);
t63 = t28 * t40;
t62 = t40 * t44;
t60 = t43 * t44;
t21 = pkin(4) * sin(t37) + t65;
t59 = t44 * t21;
t58 = t44 * t27;
t57 = t44 * t28;
t38 = -qJ(4) - pkin(8);
t32 = t42 * pkin(3);
t22 = pkin(4) * cos(t37) + t32;
t54 = t44 * pkin(1) + t41 * pkin(7);
t10 = t28 * t61 - t58;
t9 = t27 * t61 + t57;
t52 = -t9 * pkin(5) + t10 * qJ(6);
t11 = -t41 * t28 + t43 * t58;
t12 = t41 * t27 + t43 * t57;
t51 = -t11 * pkin(5) + t12 * qJ(6);
t50 = g(1) * t9 - g(2) * t11;
t48 = g(1) * t41 - g(2) * t44;
t29 = t32 + pkin(2);
t47 = t43 * t29 - t40 * t38;
t20 = pkin(2) + t22;
t45 = pkin(5) * t28 + qJ(6) * t27 + t20;
t1 = g(1) * t11 + g(2) * t9 + g(3) * t64;
t3 = g(1) * t12 + g(2) * t10 + g(3) * t63;
t36 = -pkin(9) + t38;
t33 = t44 * pkin(7);
t23 = qJ(6) * t63;
t19 = t48 * t40;
t18 = t41 * t39 + t43 * t55;
t16 = -t42 * t61 + t56;
t14 = t49 * t43 + t67;
t6 = t13 * t28;
t5 = t13 * t27;
t4 = g(1) * t10 - g(2) * t12;
t2 = [0, t48, t49, 0, 0, 0, 0, 0, t48 * t43, -t19, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t18, -g(1) * t15 - g(2) * t17, t19, -g(1) * (pkin(3) * t56 + t33) - g(2) * (t29 * t60 - t38 * t62 + t54) + (-g(1) * (-pkin(1) - t47) - g(2) * t65) * t41, 0, 0, 0, 0, 0, t4, -t50, t4, t19, t50, -g(1) * (-t10 * pkin(5) - t9 * qJ(6) + t33 + t59) - g(2) * (t12 * pkin(5) + t11 * qJ(6) + t20 * t60 - t36 * t62 + t54) + (-g(1) * (-t43 * t20 + t40 * t36 - pkin(1)) - g(2) * t21) * t41; 0, 0, 0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, t13 * t42, -t13 * t39, -t14, -g(3) * t47 + t49 * (t29 * t40 + t38 * t43) 0, 0, 0, 0, 0, t6, -t5, t6, -t14, t5 (-g(3) * t45 + t49 * t36) * t43 + (g(3) * t36 + t49 * t45) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, g(1) * t18 - g(2) * t16 + t42 * t67, 0, t72 * pkin(3), 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, -g(1) * (t41 * t22 - t43 * t59 + t51) - g(2) * (-t21 * t61 - t44 * t22 + t52) - g(3) * (t23 + (-pkin(5) * t27 - t21) * t40); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, t1, 0, -t3, -g(1) * t51 - g(2) * t52 - g(3) * (-pkin(5) * t64 + t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
