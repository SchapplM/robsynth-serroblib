% Calculate inertial parameters regressor of gravitation load for
% S5PRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:27:54
% EndTime: 2019-12-05 16:27:56
% DurationCPUTime: 0.52s
% Computational Cost: add. (294->92), mult. (556->148), div. (0->0), fcn. (662->12), ass. (0->54)
t32 = sin(qJ(3));
t35 = cos(qJ(3));
t52 = cos(pkin(5));
t29 = sin(pkin(5));
t33 = sin(qJ(2));
t58 = t29 * t33;
t67 = -t32 * t58 + t52 * t35;
t36 = cos(qJ(2));
t28 = sin(pkin(9));
t48 = t28 * t52;
t51 = cos(pkin(9));
t17 = -t33 * t48 + t51 * t36;
t57 = t29 * t35;
t66 = -t17 * t32 + t28 * t57;
t65 = g(3) * t29;
t42 = t52 * t51;
t14 = t28 * t33 - t36 * t42;
t15 = t28 * t36 + t33 * t42;
t24 = pkin(3) * t35 + pkin(2);
t30 = -qJ(4) - pkin(7);
t64 = -t14 * t24 - t15 * t30;
t16 = t51 * t33 + t36 * t48;
t63 = -t16 * t24 - t17 * t30;
t27 = qJ(3) + pkin(10);
t26 = cos(t27);
t31 = sin(qJ(5));
t61 = t26 * t31;
t34 = cos(qJ(5));
t60 = t26 * t34;
t59 = t28 * t29;
t56 = t29 * t36;
t55 = t30 * t33;
t54 = t31 * t36;
t53 = t34 * t36;
t47 = t29 * t51;
t45 = t66 * pkin(3);
t25 = sin(t27);
t44 = pkin(4) * t26 + pkin(8) * t25;
t43 = t67 * pkin(3);
t10 = -t25 * t58 + t52 * t26;
t4 = -t15 * t25 - t26 * t47;
t6 = -t17 * t25 + t26 * t59;
t41 = g(1) * t6 + g(2) * t4 + g(3) * t10;
t11 = t52 * t25 + t26 * t58;
t5 = t15 * t26 - t25 * t47;
t7 = t17 * t26 + t25 * t59;
t40 = g(1) * t7 + g(2) * t5 + g(3) * t11;
t39 = -t15 * t32 - t35 * t47;
t2 = -g(1) * t16 - g(2) * t14 + g(3) * t56;
t38 = g(1) * t17 + g(2) * t15 + g(3) * t58;
t37 = t39 * pkin(3);
t18 = t24 * t56;
t1 = t2 * t25;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t38, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * t35, t2 * t32, -t38, -g(1) * (-pkin(2) * t16 + pkin(7) * t17) - g(2) * (-pkin(2) * t14 + pkin(7) * t15) - (pkin(2) * t36 + pkin(7) * t33) * t65, 0, 0, 0, 0, 0, 0, -t2 * t26, t1, -t38, -g(1) * t63 - g(2) * t64 - g(3) * (-t29 * t55 + t18), 0, 0, 0, 0, 0, 0, -g(1) * (-t16 * t60 + t17 * t31) - g(2) * (-t14 * t60 + t15 * t31) - (t26 * t53 + t31 * t33) * t65, -g(1) * (t16 * t61 + t17 * t34) - g(2) * (t14 * t61 + t15 * t34) - (-t26 * t54 + t33 * t34) * t65, -t1, -g(1) * (-t44 * t16 + t63) - g(2) * (-t44 * t14 + t64) - g(3) * t18 - (t44 * t36 - t55) * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t66 - g(2) * t39 - g(3) * t67, -g(1) * (-t17 * t35 - t32 * t59) - g(2) * (-t15 * t35 + t32 * t47) - g(3) * (-t52 * t32 - t33 * t57), 0, 0, 0, 0, 0, 0, 0, 0, -t41, t40, 0, -g(1) * t45 - g(2) * t37 - g(3) * t43, 0, 0, 0, 0, 0, 0, -t41 * t34, t41 * t31, -t40, -g(1) * (pkin(4) * t6 + pkin(8) * t7 + t45) - g(2) * (t4 * pkin(4) + t5 * pkin(8) + t37) - g(3) * (pkin(4) * t10 + pkin(8) * t11 + t43); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t16 * t34 - t31 * t7) - g(2) * (t14 * t34 - t31 * t5) - g(3) * (-t11 * t31 - t29 * t53), -g(1) * (-t16 * t31 - t34 * t7) - g(2) * (-t14 * t31 - t34 * t5) - g(3) * (-t11 * t34 + t29 * t54), 0, 0;];
taug_reg = t3;
