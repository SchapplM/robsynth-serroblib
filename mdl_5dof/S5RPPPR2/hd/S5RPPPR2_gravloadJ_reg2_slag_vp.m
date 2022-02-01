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
% Datum: 2022-01-23 09:00
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 08:59:46
% EndTime: 2022-01-23 08:59:47
% DurationCPUTime: 0.32s
% Computational Cost: add. (142->74), mult. (303->120), div. (0->0), fcn. (353->10), ass. (0->50)
t33 = cos(pkin(9));
t34 = cos(pkin(8));
t36 = sin(qJ(5));
t31 = sin(pkin(8));
t38 = cos(qJ(5));
t53 = t31 * t38;
t13 = t33 * t53 - t36 * t34;
t37 = sin(qJ(1));
t39 = cos(qJ(1));
t30 = sin(pkin(9));
t32 = sin(pkin(7));
t35 = cos(pkin(7));
t50 = t33 * t34;
t11 = t32 * t30 + t35 * t50;
t54 = t31 * t36;
t41 = t11 * t38 + t35 * t54;
t59 = -t39 * t13 + t41 * t37;
t12 = t33 * t54 + t38 * t34;
t3 = -t11 * t36 + t35 * t53;
t58 = t37 * t12 - t3 * t39;
t56 = t32 * qJ(3) + pkin(1);
t55 = t31 * qJ(4) + pkin(2);
t52 = t32 * t37;
t51 = t32 * t39;
t48 = t37 * t31;
t47 = t37 * t34;
t46 = t39 * t31;
t45 = t39 * t34;
t44 = qJ(4) * t34 - qJ(2);
t15 = -t35 * t47 + t46;
t17 = t35 * t45 + t48;
t43 = g(1) * (t15 * t30 + t33 * t52) + g(2) * (t17 * t30 - t33 * t51);
t14 = t35 * t48 + t45;
t16 = t35 * t46 - t47;
t42 = g(1) * t14 - g(2) * t16;
t24 = g(1) * t39 + g(2) * t37;
t23 = g(1) * t37 - g(2) * t39;
t21 = t33 * pkin(4) + t30 * pkin(6) + pkin(3);
t40 = -t21 * t31 + t44;
t28 = t39 * qJ(2);
t27 = t37 * qJ(2);
t20 = pkin(2) * t35 + t56;
t19 = -t31 * pkin(3) + t44;
t18 = t23 * t32;
t10 = -t35 * t30 + t32 * t50;
t8 = g(3) * t35 - t24 * t32;
t7 = (t34 * pkin(3) + t55) * t35 + t56;
t2 = (t21 * t34 + t55) * t35 + pkin(1) + (t30 * pkin(4) - t33 * pkin(6) + qJ(3)) * t32;
t1 = -g(3) * t32 * t31 - g(1) * t16 - g(2) * t14;
t4 = [0, 0, 0, 0, 0, 0, t23, t24, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t35, -t18, -t24, -g(1) * (-t37 * pkin(1) + t28) - g(2) * (t39 * pkin(1) + t27), 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -t42, t18, -g(1) * (-t20 * t37 + t28) - g(2) * (t20 * t39 + t27), 0, 0, 0, 0, 0, 0, -g(1) * (t15 * t33 - t30 * t52) - g(2) * (t17 * t33 + t30 * t51), t43, t42, -g(1) * (-t19 * t39 - t7 * t37) - g(2) * (-t19 * t37 + t7 * t39), 0, 0, 0, 0, 0, 0, g(1) * t59 - g(2) * ((t11 * t39 + t33 * t48) * t38 + t16 * t36), -g(1) * (-t39 * t12 - t3 * t37) + g(2) * t58, -t43, -g(1) * (-t2 * t37 - t40 * t39) - g(2) * (t2 * t39 - t40 * t37); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t58 - g(2) * (-(t11 * t37 - t33 * t46) * t36 + t14 * t38) - g(3) * (-t36 * t10 + t32 * t53), -g(1) * (-t37 * t13 - t39 * t41) + g(2) * t59 - g(3) * (-t10 * t38 - t32 * t54), 0, 0;];
taug_reg = t4;
