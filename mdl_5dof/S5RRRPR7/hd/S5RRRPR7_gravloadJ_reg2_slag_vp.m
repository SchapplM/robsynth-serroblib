% Calculate inertial parameters regressor of gravitation load for
% S5RRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR7_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR7_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t30 = cos(pkin(9));
t18 = pkin(4) * t30 + pkin(3);
t28 = qJ(2) + qJ(3);
t24 = sin(t28);
t25 = cos(t28);
t31 = -pkin(8) - qJ(4);
t64 = t25 * t18 - t24 * t31;
t63 = t25 * pkin(3) + t24 * qJ(4);
t33 = sin(qJ(1));
t35 = cos(qJ(1));
t13 = g(1) * t35 + g(2) * t33;
t5 = -g(3) * t25 + t13 * t24;
t32 = sin(qJ(2));
t62 = pkin(2) * t32;
t61 = pkin(3) * t24;
t34 = cos(qJ(2));
t26 = t34 * pkin(2);
t21 = t26 + pkin(1);
t16 = t35 * t21;
t59 = g(2) * t16;
t57 = g(3) * t24;
t54 = t25 * t33;
t53 = t25 * t35;
t29 = sin(pkin(9));
t52 = t29 * t33;
t51 = t29 * t35;
t50 = t30 * t33;
t49 = t30 * t35;
t47 = qJ(4) * t25;
t36 = -pkin(7) - pkin(6);
t46 = pkin(4) * t29 - t36;
t44 = -t61 - t62;
t43 = g(1) * t33 - g(2) * t35;
t40 = t18 * t24 + t25 * t31;
t37 = -g(3) * t34 + t13 * t32;
t27 = pkin(9) + qJ(5);
t23 = cos(t27);
t22 = sin(t27);
t15 = t35 * t47;
t14 = t33 * t47;
t11 = t43 * t24;
t10 = t22 * t33 + t23 * t53;
t9 = -t22 * t53 + t23 * t33;
t8 = t22 * t35 - t23 * t54;
t7 = t22 * t54 + t23 * t35;
t6 = t13 * t25 + t57;
t4 = t5 * t30;
t3 = t5 * t29;
t2 = t5 * t23;
t1 = t5 * t22;
t12 = [0, 0, 0, 0, 0, 0, t43, t13, 0, 0, 0, 0, 0, 0, 0, 0, t43 * t34, -t43 * t32, -t13, -g(1) * (-pkin(1) * t33 + pkin(6) * t35) - g(2) * (pkin(1) * t35 + pkin(6) * t33), 0, 0, 0, 0, 0, 0, t43 * t25, -t11, -t13, -g(1) * (-t33 * t21 - t35 * t36) - g(2) * (-t33 * t36 + t16), 0, 0, 0, 0, 0, 0, -g(1) * (-t25 * t50 + t51) - g(2) * (t25 * t49 + t52), -g(1) * (t25 * t52 + t49) - g(2) * (-t25 * t51 + t50), t11, -t59 + (g(1) * t36 - g(2) * t63) * t35 + (-g(1) * (-t21 - t63) + g(2) * t36) * t33, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, t11, -t59 + (-g(1) * t46 - g(2) * t64) * t35 + (-g(1) * (-t21 - t64) - g(2) * t46) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, g(3) * t32 + t13 * t34, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t37 * pkin(2), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t44 * t35 + t15) - g(2) * (t44 * t33 + t14) - g(3) * (t26 + t63), 0, 0, 0, 0, 0, 0, t2, -t1, -t6, -g(3) * (t26 + t64) + t13 * (t40 + t62); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-t35 * t61 + t15) - g(2) * (-t33 * t61 + t14) - g(3) * t63, 0, 0, 0, 0, 0, 0, t2, -t1, -t6, -g(3) * t64 + t13 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 + t22 * t57, g(1) * t10 - g(2) * t8 + t23 * t57, 0, 0;];
taug_reg = t12;
