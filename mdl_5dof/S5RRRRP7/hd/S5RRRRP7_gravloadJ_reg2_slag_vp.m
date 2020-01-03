% Calculate inertial parameters regressor of gravitation load for
% S5RRRRP7
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
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t30 = sin(qJ(4));
t33 = cos(qJ(4));
t64 = pkin(4) * t33 + qJ(5) * t30;
t32 = sin(qJ(1));
t35 = cos(qJ(1));
t12 = g(1) * t35 + g(2) * t32;
t29 = qJ(2) + qJ(3);
t26 = sin(t29);
t41 = t12 * t26;
t27 = cos(t29);
t49 = t27 * pkin(3) + t26 * pkin(8);
t31 = sin(qJ(2));
t63 = pkin(2) * t31;
t62 = pkin(3) * t26;
t58 = g(3) * t26;
t57 = g(3) * t30;
t56 = t26 * t35;
t55 = t27 * t35;
t54 = t30 * t32;
t53 = t32 * t33;
t52 = t33 * t35;
t51 = t35 * t30;
t36 = -pkin(7) - pkin(6);
t50 = t35 * t36;
t34 = cos(qJ(2));
t28 = t34 * pkin(2);
t25 = t28 + pkin(1);
t14 = t35 * t25;
t47 = pkin(3) * t55 + pkin(8) * t56 + t14;
t46 = t64 * t27 + t49;
t10 = t27 * t51 - t53;
t8 = t27 * t54 + t52;
t45 = g(1) * t8 - g(2) * t10;
t44 = -t62 - t63;
t43 = g(1) * t32 - g(2) * t35;
t1 = g(1) * t10 + g(2) * t8 + t26 * t57;
t11 = t27 * t52 + t54;
t9 = t27 * t53 - t51;
t40 = g(1) * t11 + g(2) * t9 + t33 * t58;
t5 = -g(3) * t27 + t41;
t39 = -g(3) * t34 + t12 * t31;
t38 = (-g(1) * (-t25 - t49) + g(2) * t36) * t32;
t37 = (pkin(3) + t64) * t41;
t18 = pkin(8) * t55;
t15 = t32 * t27 * pkin(8);
t7 = t43 * t26;
t6 = t12 * t27 + t58;
t4 = t5 * t33;
t3 = -t27 * t57 + t30 * t41;
t2 = g(1) * t9 - g(2) * t11;
t13 = [0, 0, 0, 0, 0, 0, t43, t12, 0, 0, 0, 0, 0, 0, 0, 0, t43 * t34, -t43 * t31, -t12, -g(1) * (-pkin(1) * t32 + pkin(6) * t35) - g(2) * (pkin(1) * t35 + pkin(6) * t32), 0, 0, 0, 0, 0, 0, t43 * t27, -t7, -t12, -g(1) * (-t32 * t25 - t50) - g(2) * (-t32 * t36 + t14), 0, 0, 0, 0, 0, 0, t2, -t45, t7, g(1) * t50 - g(2) * t47 + t38, 0, 0, 0, 0, 0, 0, t2, t7, t45, -g(1) * (-t9 * pkin(4) - t8 * qJ(5) - t50) - g(2) * (t11 * pkin(4) + t10 * qJ(5) + t47) + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, g(3) * t31 + t12 * t34, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, t39 * pkin(2), 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (t44 * t35 + t18) - g(2) * (t44 * t32 + t15) - g(3) * (t28 + t49), 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * (-t35 * t63 + t18) - g(2) * (-t32 * t63 + t15) - g(3) * (t28 + t46) + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, -t6, -g(1) * (-pkin(3) * t56 + t18) - g(2) * (-t32 * t62 + t15) - g(3) * t49, 0, 0, 0, 0, 0, 0, t4, -t6, t3, -g(1) * t18 - g(2) * t15 - g(3) * t46 + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t40, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t40, -g(1) * (-pkin(4) * t10 + qJ(5) * t11) - g(2) * (-pkin(4) * t8 + qJ(5) * t9) - (-pkin(4) * t30 + qJ(5) * t33) * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t13;
