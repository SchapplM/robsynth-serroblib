% Calculate inertial parameters regressor of gravitation load for
% S5RRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP9_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t39 = cos(qJ(3));
t29 = t39 * pkin(3) + pkin(2);
t40 = cos(qJ(2));
t22 = t40 * t29;
t37 = sin(qJ(2));
t42 = -pkin(8) - pkin(7);
t72 = -t37 * t42 + t22;
t38 = sin(qJ(1));
t41 = cos(qJ(1));
t20 = g(1) * t41 + g(2) * t38;
t43 = -g(3) * t40 + t20 * t37;
t71 = g(1) * t38;
t68 = g(3) * t37;
t35 = qJ(3) + qJ(4);
t30 = sin(t35);
t66 = t30 * t37;
t31 = cos(t35);
t65 = t31 * t37;
t36 = sin(qJ(3));
t63 = t38 * t36;
t62 = t38 * t39;
t61 = t38 * t40;
t60 = t41 * t30;
t59 = t41 * t31;
t58 = t41 * t36;
t57 = t41 * t39;
t56 = t41 * pkin(1) + t38 * pkin(6);
t55 = t40 * t58;
t10 = t31 * t61 - t60;
t9 = t30 * t61 + t59;
t54 = -t9 * pkin(4) + t10 * qJ(5);
t11 = -t38 * t31 + t40 * t60;
t12 = t38 * t30 + t40 * t59;
t53 = -t11 * pkin(4) + t12 * qJ(5);
t52 = g(1) * t9 - g(2) * t11;
t51 = t40 * pkin(2) + t37 * pkin(7);
t49 = -g(2) * t41 + t71;
t48 = pkin(4) * t31 + qJ(5) * t30;
t14 = t36 * t61 + t57;
t45 = pkin(3) * t63 + t72 * t41 + t56;
t1 = g(1) * t11 + g(2) * t9 + g(3) * t66;
t3 = g(1) * t12 + g(2) * t10 + g(3) * t65;
t33 = t41 * pkin(6);
t44 = pkin(3) * t58 + t33 + (-pkin(1) - t72) * t38;
t13 = t20 * t40 + t68;
t27 = pkin(3) * t62;
t21 = qJ(5) * t65;
t18 = t49 * t37;
t17 = t40 * t57 + t63;
t16 = -t55 + t62;
t15 = -t39 * t61 + t58;
t6 = t43 * t31;
t5 = t43 * t30;
t4 = g(1) * t10 - g(2) * t12;
t2 = [0, 0, 0, 0, 0, 0, t49, t20, 0, 0, 0, 0, 0, 0, 0, 0, t49 * t40, -t18, -t20, -g(1) * (-t38 * pkin(1) + t33) - g(2) * t56, 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t18, -g(1) * t33 - g(2) * (t51 * t41 + t56) - (-pkin(1) - t51) * t71, 0, 0, 0, 0, 0, 0, t4, -t52, t18, -g(1) * t44 - g(2) * t45, 0, 0, 0, 0, 0, 0, t4, t18, t52, -g(1) * (-t10 * pkin(4) - t9 * qJ(5) + t44) - g(2) * (t12 * pkin(4) + t11 * qJ(5) + t45); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t13, 0, 0, 0, 0, 0, 0, 0, 0, t43 * t39, -t43 * t36, -t13, -g(3) * t51 + t20 * (pkin(2) * t37 - pkin(7) * t40), 0, 0, 0, 0, 0, 0, t6, -t5, -t13, -g(3) * t72 + t20 * (t29 * t37 + t40 * t42), 0, 0, 0, 0, 0, 0, t6, -t13, t5, -g(3) * t22 + (-g(3) * t48 + t20 * t42) * t40 + (g(3) * t42 + t20 * (t29 + t48)) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t16 + g(2) * t14 + t36 * t68, g(1) * t17 - g(2) * t15 + t39 * t68, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, -g(1) * t27 + (g(2) * t57 + t13 * t36) * pkin(3), 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * (-pkin(3) * t55 + t27 + t53) - g(2) * (-t14 * pkin(3) + t54) - g(3) * (t21 + (-pkin(3) * t36 - pkin(4) * t30) * t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * t53 - g(2) * t54 - g(3) * (-pkin(4) * t66 + t21); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t2;
