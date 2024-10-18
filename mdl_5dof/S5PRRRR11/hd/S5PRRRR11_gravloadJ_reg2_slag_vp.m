% Calculate inertial parameters regressor of gravitation load for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:57
% EndTime: 2024-09-27 21:45:57
% DurationCPUTime: 0.50s
% Computational Cost: add. (518->81), mult. (339->122), div. (0->0), fcn. (322->20), ass. (0->66)
t44 = pkin(10) + qJ(2);
t37 = sin(t44);
t38 = cos(t44);
t49 = sin(qJ(3));
t47 = cos(pkin(5));
t51 = cos(qJ(3));
t69 = t47 * t51;
t14 = -t37 * t69 - t38 * t49;
t53 = -t37 * t49 + t38 * t69;
t46 = sin(pkin(5));
t71 = g(3) * t46;
t76 = -g(1) * t14 - g(2) * t53 - t51 * t71;
t74 = pkin(7) + pkin(8);
t73 = pkin(7) * t46;
t36 = t51 * pkin(3) + pkin(2);
t70 = t47 * t49;
t48 = sin(qJ(4));
t68 = t48 * t49;
t45 = qJ(3) + qJ(4);
t39 = pkin(5) + t45;
t32 = qJ(5) + t39;
t66 = sin(t32) / 0.2e1;
t65 = sin(t39) / 0.2e1;
t64 = pkin(5) - t45;
t63 = -qJ(5) + t64;
t62 = g(1) * t38 + g(2) * t37;
t61 = g(1) * t37 - g(2) * t38;
t60 = sin(t64);
t52 = sin(t63);
t17 = t66 - t52 / 0.2e1;
t42 = qJ(5) + t45;
t34 = cos(t42);
t59 = t38 * t17 + t37 * t34;
t58 = t37 * t17 - t38 * t34;
t19 = t65 - t60 / 0.2e1;
t41 = cos(t45);
t57 = t38 * t19 + t37 * t41;
t56 = t37 * t19 - t38 * t41;
t50 = cos(qJ(4));
t35 = t50 * pkin(4) + pkin(3);
t54 = -pkin(4) * t68 + t35 * t51;
t40 = sin(t45);
t33 = sin(t42);
t31 = cos(t64);
t30 = cos(t63);
t29 = cos(t39) / 0.2e1;
t26 = cos(t32) / 0.2e1;
t23 = -t49 * pkin(3) - pkin(4) * t40;
t22 = pkin(4) * t41 + t36;
t21 = pkin(3) * t70 - t46 * t74;
t20 = t31 / 0.2e1 + t29;
t18 = t30 / 0.2e1 + t26;
t16 = t62 * t46;
t15 = -t37 * t70 + t38 * t51;
t13 = -t37 * t51 - t38 * t70;
t10 = t54 * t47;
t9 = -t46 * (pkin(9) + t74) + (t51 * t48 * pkin(4) + t49 * t35) * t47;
t8 = -t37 * t20 - t38 * t40;
t7 = -t38 * t20 + t37 * t40;
t6 = -t37 * t18 - t38 * t33;
t5 = -t38 * t18 + t37 * t33;
t4 = -g(1) * t56 + g(2) * t57 - g(3) * (t29 - t31 / 0.2e1);
t3 = -g(1) * t8 + g(2) * t7 - g(3) * (t65 + t60 / 0.2e1);
t2 = -g(1) * t58 + g(2) * t59 - g(3) * (t26 - t30 / 0.2e1);
t1 = -g(1) * t6 + g(2) * t5 - g(3) * (t66 + t52 / 0.2e1);
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t62, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t13 - g(2) * t15, g(1) * t53 - g(2) * t14, -t16, -g(1) * (-t37 * pkin(2) + t38 * t73) - g(2) * (t38 * pkin(2) + t37 * t73), 0, 0, 0, 0, 0, 0, g(1) * t57 + g(2) * t56, -g(1) * t7 - g(2) * t8, -t16, -g(1) * (-t38 * t21 - t37 * t36) - g(2) * (-t37 * t21 + t38 * t36), 0, 0, 0, 0, 0, 0, g(1) * t59 + g(2) * t58, -g(1) * t5 - g(2) * t6, -t16, -g(1) * (-t37 * t22 - t38 * t9) - g(2) * (t38 * t22 - t37 * t9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, g(1) * t15 - g(2) * t13 + t49 * t71, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t76 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t37 * t10 + t38 * t23) - g(2) * (t38 * t10 + t37 * t23) - t54 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, (t62 * t40 + (t47 * t61 - t71) * (t50 * t51 - t68)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t11;
