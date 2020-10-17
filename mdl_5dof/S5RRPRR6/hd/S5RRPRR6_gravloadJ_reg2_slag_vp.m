% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:05:54
% EndTime: 2020-01-03 12:05:55
% DurationCPUTime: 0.25s
% Computational Cost: add. (298->63), mult. (283->85), div. (0->0), fcn. (294->10), ass. (0->51)
t38 = qJ(1) + qJ(2);
t32 = sin(t38);
t34 = cos(t38);
t43 = cos(qJ(4));
t40 = cos(pkin(9));
t41 = sin(qJ(4));
t55 = t40 * t41;
t11 = -t32 * t55 - t34 * t43;
t13 = -t32 * t43 + t34 * t55;
t39 = sin(pkin(9));
t62 = g(1) * t39;
t64 = -g(2) * t11 - g(3) * t13 + t41 * t62;
t63 = pkin(7) * t39;
t59 = t32 * t40;
t58 = t32 * t41;
t57 = t34 * t40;
t56 = t39 * (-pkin(8) - pkin(7));
t54 = t40 * t43;
t53 = t34 * pkin(2) + t32 * qJ(3);
t27 = t32 * pkin(2);
t51 = -t34 * qJ(3) + t27;
t50 = pkin(3) * t57 + t34 * t63 + t53;
t20 = g(2) * t34 + g(3) * t32;
t42 = sin(qJ(1));
t44 = cos(qJ(1));
t49 = -g(2) * t44 - g(3) * t42;
t48 = pkin(3) * t59 + t32 * t63 + t51;
t30 = t43 * pkin(4) + pkin(3);
t47 = pkin(4) * t58 + t30 * t57 - t34 * t56 + t53;
t46 = -t32 * t56 + t30 * t59 + t27 + (-pkin(4) * t41 - qJ(3)) * t34;
t37 = qJ(4) + qJ(5);
t36 = t44 * pkin(1);
t35 = t42 * pkin(1);
t33 = cos(t37);
t31 = sin(t37);
t19 = g(2) * t32 - g(3) * t34;
t16 = t20 * t40;
t15 = t20 * t39;
t14 = t34 * t54 + t58;
t12 = t32 * t54 - t34 * t41;
t10 = t32 * t31 + t33 * t57;
t9 = t31 * t57 - t32 * t33;
t8 = -t34 * t31 + t33 * t59;
t7 = -t31 * t59 - t34 * t33;
t6 = -g(2) * t14 - g(3) * t12;
t5 = g(2) * t13 - g(3) * t11;
t4 = -g(2) * t10 - g(3) * t8;
t3 = g(2) * t9 - g(3) * t7;
t2 = g(2) * t8 - g(3) * t10 + t33 * t62;
t1 = -g(2) * t7 - g(3) * t9 + t31 * t62;
t17 = [0, 0, 0, 0, 0, 0, t49, g(2) * t42 - g(3) * t44, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t19, 0, t49 * pkin(1), 0, 0, 0, 0, 0, 0, -t16, t15, -t19, -g(2) * (t36 + t53) - g(3) * (t35 + t51), 0, 0, 0, 0, 0, 0, t6, t5, -t15, -g(2) * (t36 + t50) - g(3) * (t35 + t48), 0, 0, 0, 0, 0, 0, t4, t3, -t15, -g(2) * (t36 + t47) - g(3) * (t35 + t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t19, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t15, -t19, -g(2) * t53 - g(3) * t51, 0, 0, 0, 0, 0, 0, t6, t5, -t15, -g(2) * t50 - g(3) * t48, 0, 0, 0, 0, 0, 0, t4, t3, -t15, -g(2) * t47 - g(3) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, g(2) * t12 - g(3) * t14 + t43 * t62, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t64 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t17;
