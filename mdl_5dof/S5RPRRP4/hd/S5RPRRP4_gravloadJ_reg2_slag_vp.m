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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:32:52
% EndTime: 2022-01-23 09:32:53
% DurationCPUTime: 0.26s
% Computational Cost: add. (201->60), mult. (293->88), div. (0->0), fcn. (308->8), ass. (0->49)
t32 = qJ(3) + qJ(4);
t24 = cos(t32);
t37 = cos(qJ(3));
t28 = t37 * pkin(3);
t17 = pkin(4) * t24 + t28;
t33 = sin(pkin(8));
t34 = cos(pkin(8));
t58 = pkin(7) + pkin(6);
t62 = -(-qJ(5) - t58) * t33 + (pkin(2) + t17) * t34;
t38 = cos(qJ(1));
t42 = t38 * t37;
t35 = sin(qJ(3));
t36 = sin(qJ(1));
t48 = t36 * t35;
t10 = t34 * t48 + t42;
t43 = t38 * t35;
t47 = t36 * t37;
t12 = -t34 * t43 + t47;
t55 = g(3) * t33;
t61 = -g(1) * t12 + g(2) * t10 + t35 * t55;
t23 = sin(t32);
t44 = t38 * t24;
t50 = t36 * t23;
t6 = t34 * t50 + t44;
t45 = t38 * t23;
t49 = t36 * t24;
t8 = -t34 * t45 + t49;
t1 = -g(1) * t8 + g(2) * t6 + t23 * t55;
t54 = t35 * pkin(3);
t53 = pkin(2) * t34 + pkin(1);
t16 = pkin(4) * t23 + t54;
t46 = t38 * t16;
t25 = t36 * qJ(2);
t41 = t38 * pkin(1) + t25;
t20 = g(1) * t38 + g(2) * t36;
t19 = g(1) * t36 - g(2) * t38;
t26 = t38 * qJ(2);
t21 = qJ(2) + t54;
t18 = t33 * pkin(6) + t53;
t14 = t19 * t33;
t13 = t34 * t42 + t48;
t11 = -t34 * t47 + t43;
t9 = t34 * t44 + t50;
t7 = -t34 * t49 + t45;
t5 = t34 * t28 + t58 * t33 + t53;
t4 = -g(1) * t7 - g(2) * t9;
t3 = -g(1) * t6 - g(2) * t8;
t2 = g(1) * t9 - g(2) * t7 + t24 * t55;
t15 = [0, 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, 0, 0, 0, t19 * t34, -t14, -t20, -g(1) * (-t36 * pkin(1) + t26) - g(2) * t41, 0, 0, 0, 0, 0, 0, -g(1) * t11 - g(2) * t13, -g(1) * t10 - g(2) * t12, t14, -g(1) * (-t18 * t36 + t26) - g(2) * (t18 * t38 + t25), 0, 0, 0, 0, 0, 0, t4, t3, t14, -g(1) * (t21 * t38 - t5 * t36) - g(2) * (t21 * t36 + t5 * t38), 0, 0, 0, 0, 0, 0, t4, t3, t14, -g(1) * (t26 + t46) - g(2) * (t62 * t38 + t41) + (-g(1) * (-pkin(1) - t62) - g(2) * t16) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, g(1) * t13 - g(2) * t11 + t37 * t55, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t61 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (t36 * t17 - t34 * t46) - g(2) * (-t36 * t34 * t16 - t38 * t17) + t16 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(3) * t34 - t20 * t33;];
taug_reg = t15;
