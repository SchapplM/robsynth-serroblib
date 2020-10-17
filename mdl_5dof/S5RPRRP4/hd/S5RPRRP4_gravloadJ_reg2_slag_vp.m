% Calculate inertial parameters regressor of gravitation load for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:50:24
% EndTime: 2020-01-03 11:50:25
% DurationCPUTime: 0.27s
% Computational Cost: add. (201->63), mult. (301->92), div. (0->0), fcn. (316->8), ass. (0->50)
t30 = cos(pkin(8));
t31 = sin(qJ(3));
t34 = cos(qJ(1));
t42 = t34 * t31;
t32 = sin(qJ(1));
t33 = cos(qJ(3));
t45 = t32 * t33;
t11 = t30 * t42 - t45;
t29 = sin(pkin(8));
t55 = g(1) * t29;
t41 = t34 * t33;
t46 = t32 * t31;
t9 = -t30 * t46 - t41;
t59 = -g(2) * t9 - g(3) * t11 + t31 * t55;
t28 = qJ(3) + qJ(4);
t21 = sin(t28);
t22 = cos(t28);
t43 = t34 * t22;
t48 = t32 * t21;
t5 = -t30 * t48 - t43;
t44 = t34 * t21;
t47 = t32 * t22;
t7 = t30 * t44 - t47;
t1 = -g(2) * t5 - g(3) * t7 + t21 * t55;
t35 = -pkin(7) - pkin(6);
t53 = t31 * pkin(3);
t52 = (-qJ(5) + t35) * t29;
t51 = t29 * t35;
t50 = t30 * t32;
t15 = pkin(4) * t21 + t53;
t49 = t32 * t15;
t25 = t33 * pkin(3);
t16 = pkin(4) * t22 + t25;
t40 = t34 * pkin(1) + t32 * qJ(2);
t24 = t32 * pkin(1);
t37 = -t34 * qJ(2) + t24;
t36 = pkin(2) * t30 + pkin(6) * t29;
t18 = g(2) * t34 + g(3) * t32;
t17 = g(2) * t32 - g(3) * t34;
t20 = t25 + pkin(2);
t14 = pkin(2) + t16;
t13 = t18 * t29;
t12 = t30 * t41 + t46;
t10 = t30 * t45 - t42;
t8 = t30 * t43 + t48;
t6 = t30 * t47 - t44;
t4 = -g(2) * t8 - g(3) * t6;
t3 = g(2) * t7 - g(3) * t5;
t2 = g(2) * t6 - g(3) * t8 + t22 * t55;
t19 = [0, 0, 0, 0, 0, 0, -t18, t17, 0, 0, 0, 0, 0, 0, 0, 0, -t18 * t30, t13, -t17, -g(2) * t40 - g(3) * t37, 0, 0, 0, 0, 0, 0, -g(2) * t12 - g(3) * t10, g(2) * t11 - g(3) * t9, -t13, -g(2) * (t36 * t34 + t40) - g(3) * (t36 * t32 + t37), 0, 0, 0, 0, 0, 0, t4, t3, -t13, -g(2) * (pkin(3) * t46 + t40) - g(3) * (t20 * t50 - t32 * t51 + t24) + (-g(2) * (t20 * t30 - t51) - g(3) * (-qJ(2) - t53)) * t34, 0, 0, 0, 0, 0, 0, t4, t3, -t13, -g(2) * (t40 + t49) - g(3) * (t14 * t50 - t32 * t52 + t24) + (-g(2) * (t14 * t30 - t52) - g(3) * (-qJ(2) - t15)) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, g(2) * t10 - g(3) * t12 + t33 * t55, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t59 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, t15 * t55 - g(2) * (-t34 * t16 - t30 * t49) - g(3) * (t34 * t30 * t15 - t32 * t16); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t30 - t17 * t29;];
taug_reg = t19;
