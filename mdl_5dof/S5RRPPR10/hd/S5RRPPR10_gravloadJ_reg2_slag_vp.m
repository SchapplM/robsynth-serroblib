% Calculate inertial parameters regressor of gravitation load for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPPR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t38 = cos(qJ(2));
t36 = sin(qJ(1));
t39 = cos(qJ(1));
t16 = g(1) * t39 + g(2) * t36;
t35 = sin(qJ(2));
t72 = t16 * t35;
t7 = -g(3) * t38 + t72;
t71 = g(1) * t36;
t68 = g(3) * t35;
t28 = t38 * pkin(2);
t33 = cos(pkin(8));
t66 = t33 * t38;
t65 = t35 * t36;
t64 = t35 * t39;
t63 = t36 * t38;
t62 = t38 * t39;
t32 = sin(pkin(8));
t61 = t39 * t32;
t60 = t39 * t33;
t26 = t35 * qJ(3);
t59 = t28 + t26;
t58 = t39 * pkin(1) + t36 * pkin(6);
t57 = qJ(3) * t38;
t56 = qJ(4) * t32;
t55 = -pkin(1) - t28;
t54 = -pkin(2) - t56;
t53 = pkin(3) * t66 + t38 * t56 + t59;
t52 = pkin(2) * t62 + t39 * t26 + t58;
t11 = t32 * t63 + t60;
t12 = t33 * t63 - t61;
t29 = t39 * pkin(6);
t51 = -t12 * pkin(3) - t11 * qJ(4) + t29;
t13 = -t36 * t33 + t38 * t61;
t50 = g(1) * t11 - g(2) * t13;
t49 = -g(2) * t39 + t71;
t34 = sin(qJ(5));
t37 = cos(qJ(5));
t48 = t11 * t37 - t12 * t34;
t47 = t11 * t34 + t12 * t37;
t46 = t32 * t37 - t33 * t34;
t45 = t32 * t34 + t33 * t37;
t43 = g(3) * t46;
t14 = t36 * t32 + t38 * t60;
t41 = t14 * pkin(3) + t13 * qJ(4) + t52;
t40 = (t55 - t26) * t71;
t21 = t39 * t57;
t18 = t36 * t57;
t15 = g(1) * t65 - g(2) * t64;
t8 = t16 * t38 + t68;
t6 = t7 * t33;
t5 = t7 * t32;
t4 = g(1) * t12 - g(2) * t14;
t3 = t13 * t34 + t14 * t37;
t2 = t13 * t37 - t14 * t34;
t1 = -g(1) * t13 - g(2) * t11 - t32 * t68;
t9 = [0, 0, 0, 0, 0, 0, t49, t16, 0, 0, 0, 0, 0, 0, 0, 0, t49 * t38, -t15, -t16, -g(1) * (-t36 * pkin(1) + t29) - g(2) * t58, 0, 0, 0, 0, 0, 0, t4, -t50, t15, -g(1) * t29 - g(2) * t52 - t40, 0, 0, 0, 0, 0, 0, t4, t15, t50, -g(1) * t51 - g(2) * t41 - t40, 0, 0, 0, 0, 0, 0, g(1) * t47 - g(2) * t3, g(1) * t48 - g(2) * t2, -t15, -g(1) * (-t12 * pkin(4) + t51) - g(2) * (t14 * pkin(4) - pkin(7) * t64 + t41) - ((pkin(7) - qJ(3)) * t35 + t55) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t8, -g(1) * (-pkin(2) * t64 + t21) - g(2) * (-pkin(2) * t65 + t18) - g(3) * t59, 0, 0, 0, 0, 0, 0, t6, -t8, t5, -g(1) * t21 - g(2) * t18 - g(3) * t53 + (pkin(3) * t33 - t54) * t72, 0, 0, 0, 0, 0, 0, t7 * t45, -t38 * t43 + t46 * t72, t8, -g(1) * (-pkin(7) * t62 + t21) - g(2) * (-pkin(7) * t63 + t18) - g(3) * (pkin(4) * t66 + t53) + (g(3) * pkin(7) + t16 * (-(-pkin(3) - pkin(4)) * t33 - t54)) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 - g(2) * t48 - t35 * t43, g(1) * t3 + g(2) * t47 + t45 * t68, 0, 0;];
taug_reg = t9;
