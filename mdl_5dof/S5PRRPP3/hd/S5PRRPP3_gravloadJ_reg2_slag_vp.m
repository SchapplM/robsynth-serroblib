% Calculate inertial parameters regressor of gravitation load for
% S5PRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:13:32
% EndTime: 2019-12-05 16:13:34
% DurationCPUTime: 0.28s
% Computational Cost: add. (153->69), mult. (416->98), div. (0->0), fcn. (456->8), ass. (0->51)
t35 = sin(qJ(2));
t31 = sin(pkin(7));
t33 = cos(pkin(7));
t47 = g(1) * t33 + g(2) * t31;
t64 = t47 * t35;
t63 = pkin(2) * t35;
t37 = cos(qJ(2));
t62 = pkin(6) * t37;
t34 = sin(qJ(3));
t59 = t34 * t35;
t58 = t34 * t37;
t32 = cos(pkin(8));
t57 = t35 * t32;
t36 = cos(qJ(3));
t56 = t35 * t36;
t55 = t36 * t37;
t54 = t37 * t32;
t53 = t37 * pkin(2) + t35 * pkin(6);
t52 = qJ(4) * t34;
t51 = g(3) * t59;
t14 = t31 * t58 + t33 * t36;
t15 = t31 * t55 - t33 * t34;
t50 = -t14 * pkin(3) + t15 * qJ(4);
t16 = -t31 * t36 + t33 * t58;
t17 = t31 * t34 + t33 * t55;
t49 = -t16 * pkin(3) + t17 * qJ(4);
t48 = pkin(3) * t55 + t37 * t52 + t53;
t30 = sin(pkin(8));
t46 = -pkin(4) * t32 - qJ(5) * t30;
t44 = t30 * t37 - t32 * t56;
t43 = t30 * t56 + t54;
t18 = t30 * t55 - t57;
t7 = t43 * t31;
t9 = t43 * t33;
t42 = g(1) * t9 + g(2) * t7 - g(3) * t18;
t41 = g(1) * t16 + g(2) * t14 + t51;
t40 = g(1) * t17 + g(2) * t15 + g(3) * t56;
t39 = -g(3) * t37 + t64;
t38 = (pkin(3) * t36 + pkin(2) + t52) * t64;
t22 = t33 * t62;
t21 = qJ(4) * t56;
t20 = t31 * t62;
t19 = t35 * t30 + t36 * t54;
t11 = g(3) * t35 + t47 * t37;
t10 = t44 * t33;
t8 = t44 * t31;
t6 = t39 * t34;
t3 = t41 * t32;
t2 = t41 * t30;
t1 = -g(1) * t10 - g(2) * t8 - g(3) * t19;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t11, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t36, -t6, -t11, -g(1) * (-t33 * t63 + t22) - g(2) * (-t31 * t63 + t20) - g(3) * t53, 0, 0, 0, 0, 0, 0, t1, -t42, t6, -g(1) * t22 - g(2) * t20 - g(3) * t48 + t38, 0, 0, 0, 0, 0, 0, t1, t6, t42, -g(1) * (t10 * pkin(4) - t9 * qJ(5) + t22) - g(2) * (t8 * pkin(4) - t7 * qJ(5) + t20) - g(3) * (t19 * pkin(4) + t18 * qJ(5) + t48) + t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t40, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t2, -t40, -g(1) * t49 - g(2) * t50 - g(3) * (-pkin(3) * t59 + t21), 0, 0, 0, 0, 0, 0, t3, -t40, t2, -g(1) * (t46 * t16 + t49) - g(2) * (t46 * t14 + t50) - g(3) * t21 - (-pkin(3) + t46) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t17 * t30 - t33 * t57) - g(2) * (t15 * t30 - t31 * t57) - g(3) * t43;];
taug_reg = t4;
