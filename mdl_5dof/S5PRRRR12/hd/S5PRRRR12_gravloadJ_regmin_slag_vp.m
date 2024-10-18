% Calculate minimal parameter regressor of gravitation load for
% S5PRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-28 18:09
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR12_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR12_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-28 18:08:20
% EndTime: 2024-09-28 18:08:20
% DurationCPUTime: 0.26s
% Computational Cost: add. (441->73), mult. (390->126), div. (0->0), fcn. (446->22), ass. (0->59)
t39 = sin(pkin(5));
t58 = g(3) * t39;
t38 = sin(pkin(6));
t57 = t38 * t39;
t41 = cos(pkin(6));
t43 = sin(qJ(5));
t56 = t41 * t43;
t45 = cos(qJ(5));
t55 = t41 * t45;
t42 = cos(pkin(5));
t54 = t42 * t43;
t44 = sin(qJ(2));
t53 = t42 * t44;
t52 = t42 * t45;
t46 = cos(qJ(2));
t51 = t42 * t46;
t36 = qJ(2) + qJ(3);
t50 = t43 * t57;
t49 = t45 * t57;
t48 = t41 * t54;
t47 = t41 * t52;
t32 = pkin(5) - t36;
t31 = pkin(5) + t36;
t40 = cos(pkin(11));
t37 = sin(pkin(11));
t35 = qJ(4) + t36;
t34 = cos(t36);
t33 = sin(t36);
t30 = cos(t35);
t29 = sin(t35);
t28 = -qJ(4) + t32;
t27 = qJ(4) + t31;
t26 = cos(t32);
t25 = sin(t31);
t24 = cos(t28);
t23 = sin(t27);
t22 = cos(t31) / 0.2e1;
t21 = sin(t32) / 0.2e1;
t20 = cos(t27) / 0.2e1;
t19 = sin(t28) / 0.2e1;
t18 = t26 / 0.2e1 + t22;
t17 = t21 - t25 / 0.2e1;
t16 = t24 / 0.2e1 + t20;
t15 = t19 - t23 / 0.2e1;
t14 = t37 * t54 - t40 * t55;
t13 = t37 * t52 + t40 * t56;
t12 = -t37 * t55 - t40 * t54;
t11 = t37 * t56 - t40 * t52;
t10 = -t37 * t43 + t40 * t47;
t9 = t37 * t45 + t40 * t48;
t8 = -t37 * t47 - t40 * t43;
t7 = t37 * t48 - t40 * t45;
t6 = -g(1) * (-t37 * t17 - t40 * t34) - g(2) * (t40 * t17 - t37 * t34) - g(3) * (t22 - t26 / 0.2e1);
t5 = -g(1) * (-t37 * t18 - t40 * t33) - g(2) * (t40 * t18 - t37 * t33) - g(3) * (t25 / 0.2e1 + t21);
t4 = -g(1) * (-t37 * t15 - t40 * t30) - g(2) * (t40 * t15 - t37 * t30) - g(3) * (t20 - t24 / 0.2e1);
t3 = -g(1) * (-t37 * t16 - t40 * t29) - g(2) * (t40 * t16 - t37 * t29) - g(3) * (t23 / 0.2e1 + t19);
t2 = -g(1) * (t14 * t30 - t8 * t29) - g(2) * (-t10 * t29 + t12 * t30) - (-t29 * t55 - t30 * t43) * t58;
t1 = -g(1) * (-t13 * t30 + t7 * t29) - g(2) * (-t11 * t30 - t9 * t29) - (-t29 * t56 + t30 * t45) * t58;
t59 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -g(1) * (-t37 * t51 - t40 * t44) - g(2) * (-t37 * t44 + t40 * t51) - t46 * t58, -g(1) * (t37 * t53 - t40 * t46) - g(2) * (-t37 * t46 - t40 * t53) + t44 * t58, 0, t5, t6, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, t5, t6, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t14 * t29 + t8 * t30 + t37 * t49) - g(2) * (t10 * t30 + t12 * t29 - t40 * t49) - g(3) * (t38 * t52 + (-t29 * t43 + t30 * t55) * t39), -g(1) * (t13 * t29 + t7 * t30 - t37 * t50) - g(2) * (t11 * t29 - t9 * t30 + t40 * t50) - g(3) * (-t38 * t54 + (-t29 * t45 - t30 * t56) * t39);];
taug_reg = t59;
