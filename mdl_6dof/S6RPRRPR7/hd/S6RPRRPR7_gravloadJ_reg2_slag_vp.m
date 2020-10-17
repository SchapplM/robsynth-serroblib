% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR7_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:01:48
% EndTime: 2019-05-05 23:01:49
% DurationCPUTime: 0.31s
% Computational Cost: add. (298->81), mult. (312->92), div. (0->0), fcn. (295->10), ass. (0->51)
t33 = qJ(3) + qJ(4);
t25 = pkin(10) + t33;
t23 = sin(t25);
t24 = cos(t25);
t45 = t23 * pkin(5) - t24 * pkin(9);
t36 = sin(qJ(1));
t39 = cos(qJ(1));
t59 = g(2) * t39;
t19 = g(1) * t36 - t59;
t63 = -g(3) * t23 + t19 * t24;
t40 = -pkin(8) - pkin(7);
t26 = sin(t33);
t62 = pkin(4) * t26;
t27 = cos(t33);
t61 = pkin(4) * t27;
t60 = pkin(5) * t24;
t57 = g(3) * t24;
t35 = sin(qJ(3));
t55 = t35 * pkin(3);
t54 = t23 * t36;
t34 = sin(qJ(6));
t53 = t36 * t34;
t37 = cos(qJ(6));
t52 = t36 * t37;
t51 = t39 * t34;
t50 = t39 * t37;
t49 = pkin(9) * t54 + t36 * t60;
t48 = t39 * pkin(1) + t36 * qJ(2);
t29 = t39 * qJ(2);
t47 = -t36 * pkin(1) + t29;
t46 = -pkin(9) * t23 - t60;
t20 = g(1) * t39 + g(2) * t36;
t14 = t55 + t62;
t32 = -qJ(5) + t40;
t44 = t39 * t14 + t36 * t32 + t47;
t43 = t36 * t14 - t39 * t32 + t48;
t42 = -t62 - t45;
t6 = g(3) * t26 - t19 * t27;
t38 = cos(qJ(3));
t41 = g(3) * t35 - t19 * t38;
t15 = t38 * pkin(3) + t61;
t11 = t23 * t50 - t53;
t10 = t23 * t51 + t52;
t9 = t23 * t52 + t51;
t8 = -t23 * t53 + t50;
t7 = t20 * t24;
t5 = g(3) * t27 + t19 * t26;
t3 = g(1) * t54 - t23 * t59 + t57;
t2 = t63 * t37;
t1 = t63 * t34;
t4 = [0, 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t20, -g(1) * t47 - g(2) * t48, 0, 0, 0, 0, 0, 0, -t20 * t35, -t20 * t38, t19, -g(1) * (t29 + (-pkin(1) - pkin(7)) * t36) - g(2) * (t39 * pkin(7) + t48) 0, 0, 0, 0, 0, 0, -t20 * t26, -t20 * t27, t19, -g(1) * (t39 * t55 + t29 + (-pkin(1) + t40) * t36) - g(2) * (t36 * t55 - t39 * t40 + t48) 0, 0, 0, 0, 0, 0, -t20 * t23, -t7, t19, -g(1) * t44 - g(2) * t43, 0, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t9, g(1) * t10 - g(2) * t8, t7, -g(1) * (t45 * t39 + t44) - g(2) * (t45 * t36 + t43); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, g(3) * t38 + t19 * t35, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, t41 * pkin(3), 0, 0, 0, 0, 0, 0, -t63, t3, 0, g(3) * t14 - t19 * t15, 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(1) * (t36 * t15 + t49) - g(3) * (t42 - t55) - (-t15 + t46) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t3, 0, t6 * pkin(4), 0, 0, 0, 0, 0, 0, -t2, t1, -t3, -g(1) * (t36 * t61 + t49) - g(3) * t42 - (t46 - t61) * t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10 + t34 * t57, g(1) * t9 - g(2) * t11 + t37 * t57, 0, 0;];
taug_reg  = t4;
