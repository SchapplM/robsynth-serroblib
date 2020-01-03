% Calculate inertial parameters regressor of gravitation load for
% S5RRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:03
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t37 = cos(qJ(3));
t25 = t37 * pkin(3) + pkin(2);
t38 = cos(qJ(2));
t18 = t38 * t25;
t33 = -qJ(4) - pkin(7);
t35 = sin(qJ(2));
t65 = -t35 * t33 + t18;
t36 = sin(qJ(1));
t39 = cos(qJ(1));
t17 = g(1) * t39 + g(2) * t36;
t9 = -g(3) * t38 + t17 * t35;
t64 = g(1) * t36;
t61 = g(3) * t35;
t34 = sin(qJ(3));
t58 = t36 * t34;
t57 = t36 * t37;
t56 = t36 * t38;
t32 = qJ(3) + pkin(8);
t26 = sin(t32);
t55 = t39 * t26;
t27 = cos(t32);
t54 = t39 * t27;
t53 = t39 * t34;
t52 = t39 * t37;
t51 = t39 * pkin(1) + t36 * pkin(6);
t50 = t38 * t53;
t5 = t26 * t56 + t54;
t7 = -t36 * t27 + t38 * t55;
t49 = g(1) * t5 - g(2) * t7;
t48 = t38 * pkin(2) + t35 * pkin(7);
t46 = -g(2) * t39 + t64;
t45 = pkin(4) * t27 + qJ(5) * t26;
t11 = t34 * t56 + t52;
t42 = pkin(3) * t58 + t65 * t39 + t51;
t1 = g(1) * t7 + g(2) * t5 + t26 * t61;
t6 = t27 * t56 - t55;
t8 = t36 * t26 + t38 * t54;
t41 = g(1) * t8 + g(2) * t6 + t27 * t61;
t29 = t39 * pkin(6);
t40 = pkin(3) * t53 + t29 + (-pkin(1) - t65) * t36;
t10 = t17 * t38 + t61;
t23 = pkin(3) * t57;
t15 = t46 * t35;
t14 = t38 * t52 + t58;
t13 = -t50 + t57;
t12 = -t37 * t56 + t53;
t4 = t9 * t27;
t3 = t9 * t26;
t2 = g(1) * t6 - g(2) * t8;
t16 = [0, 0, 0, 0, 0, 0, t46, t17, 0, 0, 0, 0, 0, 0, 0, 0, t46 * t38, -t15, -t17, -g(1) * (-t36 * pkin(1) + t29) - g(2) * t51, 0, 0, 0, 0, 0, 0, -g(1) * t12 - g(2) * t14, -g(1) * t11 - g(2) * t13, t15, -g(1) * t29 - g(2) * (t48 * t39 + t51) - (-pkin(1) - t48) * t64, 0, 0, 0, 0, 0, 0, t2, -t49, t15, -g(1) * t40 - g(2) * t42, 0, 0, 0, 0, 0, 0, t2, t15, t49, -g(1) * (-t6 * pkin(4) - t5 * qJ(5) + t40) - g(2) * (t8 * pkin(4) + t7 * qJ(5) + t42); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t10, 0, 0, 0, 0, 0, 0, 0, 0, t9 * t37, -t9 * t34, -t10, -g(3) * t48 + t17 * (pkin(2) * t35 - pkin(7) * t38), 0, 0, 0, 0, 0, 0, t4, -t3, -t10, -g(3) * t65 + t17 * (t25 * t35 + t33 * t38), 0, 0, 0, 0, 0, 0, t4, -t10, t3, -g(3) * t18 + (-g(3) * t45 + t17 * t33) * t38 + (g(3) * t33 + t17 * (t25 + t45)) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 + g(2) * t11 + t34 * t61, g(1) * t14 - g(2) * t12 + t37 * t61, 0, 0, 0, 0, 0, 0, 0, 0, t1, t41, 0, -g(1) * t23 + (g(2) * t52 + t10 * t34) * pkin(3), 0, 0, 0, 0, 0, 0, t1, 0, -t41, -g(1) * (-pkin(3) * t50 - t7 * pkin(4) + t8 * qJ(5) + t23) - g(2) * (-t11 * pkin(3) - t5 * pkin(4) + t6 * qJ(5)) - (-pkin(3) * t34 - pkin(4) * t26 + qJ(5) * t27) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t16;
