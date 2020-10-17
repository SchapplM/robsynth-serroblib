% Calculate inertial parameters regressor of gravitation load for
% S5RPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPPR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR2_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:22:55
% EndTime: 2020-01-03 11:22:58
% DurationCPUTime: 0.26s
% Computational Cost: add. (142->53), mult. (345->87), div. (0->0), fcn. (411->10), ass. (0->45)
t32 = sin(pkin(7));
t34 = cos(pkin(7));
t61 = pkin(2) * t34 + qJ(3) * t32;
t36 = sin(qJ(1));
t48 = cos(pkin(8));
t47 = t36 * t48;
t31 = sin(pkin(8));
t38 = cos(qJ(1));
t51 = t38 * t31;
t15 = t34 * t51 - t47;
t35 = sin(qJ(5));
t37 = cos(qJ(5));
t46 = t38 * t48;
t52 = t36 * t31;
t16 = t34 * t46 + t52;
t30 = sin(pkin(9));
t33 = cos(pkin(9));
t53 = t32 * t38;
t7 = t16 * t33 + t30 * t53;
t60 = -t15 * t37 + t7 * t35;
t59 = t15 * t35 + t7 * t37;
t55 = t31 * t32;
t54 = t32 * t36;
t50 = t38 * pkin(1) + t36 * qJ(2);
t45 = t36 * pkin(1) - t38 * qJ(2);
t44 = t61 * t38 + t50;
t14 = t34 * t47 - t51;
t4 = t14 * t30 - t33 * t54;
t6 = t16 * t30 - t33 * t53;
t43 = g(2) * t6 + g(3) * t4;
t13 = t34 * t52 + t46;
t42 = g(2) * t15 + g(3) * t13;
t21 = g(2) * t38 + g(3) * t36;
t20 = g(2) * t36 - g(3) * t38;
t41 = t61 * t36 + t45;
t40 = t16 * pkin(3) + t15 * qJ(4) + t44;
t39 = t14 * pkin(3) + t13 * qJ(4) + t41;
t17 = t21 * t32;
t12 = t32 * t48 * t33 - t34 * t30;
t9 = g(1) * t34 - t20 * t32;
t5 = t14 * t33 + t30 * t54;
t3 = -g(1) * t55 - g(2) * t13 + g(3) * t15;
t2 = t13 * t35 + t5 * t37;
t1 = t13 * t37 - t5 * t35;
t8 = [0, 0, 0, 0, 0, 0, -t21, t20, 0, 0, 0, 0, 0, 0, 0, 0, -t21 * t34, t17, -t20, -g(2) * t50 - g(3) * t45, 0, 0, 0, 0, 0, 0, -g(2) * t16 - g(3) * t14, t42, -t17, -g(2) * t44 - g(3) * t41, 0, 0, 0, 0, 0, 0, -g(2) * t7 - g(3) * t5, t43, -t42, -g(2) * t40 - g(3) * t39, 0, 0, 0, 0, 0, 0, -g(2) * t59 - g(3) * t2, g(2) * t60 - g(3) * t1, -t43, -g(2) * (t7 * pkin(4) + t6 * pkin(6) + t40) - g(3) * (t5 * pkin(4) + t4 * pkin(6) + t39); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t12 * t35 + t37 * t55) - g(2) * t1 - g(3) * t60, -g(1) * (-t12 * t37 - t35 * t55) + g(2) * t2 - g(3) * t59, 0, 0;];
taug_reg = t8;
