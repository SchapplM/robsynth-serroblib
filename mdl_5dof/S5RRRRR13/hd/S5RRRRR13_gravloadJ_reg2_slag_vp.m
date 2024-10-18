% Calculate inertial parameters regressor of gravitation load for
% S5RRRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha4,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 17:33
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR13_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR13_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR13_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRRR13_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 17:32:16
% EndTime: 2024-09-27 17:32:16
% DurationCPUTime: 0.34s
% Computational Cost: add. (543->65), mult. (335->91), div. (0->0), fcn. (326->16), ass. (0->63)
t41 = qJ(1) + qJ(2);
t38 = qJ(3) + t41;
t30 = sin(t38);
t31 = cos(t38);
t44 = sin(qJ(4));
t43 = cos(pkin(5));
t46 = cos(qJ(4));
t64 = t43 * t46;
t11 = -t30 * t64 - t31 * t44;
t49 = -t30 * t44 + t31 * t64;
t42 = sin(pkin(5));
t67 = g(3) * t42;
t72 = -g(1) * t11 - g(2) * t49 - t46 * t67;
t35 = sin(t41);
t70 = pkin(2) * t35;
t69 = pkin(9) * t42;
t45 = sin(qJ(1));
t66 = t45 * pkin(1);
t65 = t43 * t44;
t63 = t31 * pkin(3) + t30 * t69;
t40 = qJ(4) + qJ(5);
t33 = pkin(5) + t40;
t61 = sin(t33) / 0.2e1;
t37 = cos(t41);
t29 = pkin(2) * t37;
t60 = t29 + t63;
t59 = -t30 * pkin(3) + t31 * t69;
t18 = pkin(4) * t65 - t42 * (pkin(9) + pkin(10));
t32 = t46 * pkin(4) + pkin(3);
t58 = -t30 * t18 + t31 * t32;
t57 = pkin(5) - t40;
t56 = t29 + t58;
t15 = g(1) * t31 + g(2) * t30;
t19 = g(1) * t35 - g(2) * t37;
t47 = cos(qJ(1));
t55 = g(1) * t45 - g(2) * t47;
t54 = sin(t57);
t16 = t61 - t54 / 0.2e1;
t36 = cos(t40);
t53 = t31 * t16 + t30 * t36;
t52 = t30 * t16 - t31 * t36;
t51 = -t31 * t18 - t30 * t32;
t50 = t59 - t70;
t48 = t51 - t70;
t39 = t47 * pkin(1);
t34 = sin(t40);
t28 = cos(t57);
t26 = cos(t33) / 0.2e1;
t20 = g(1) * t37 + g(2) * t35;
t17 = t28 / 0.2e1 + t26;
t14 = g(1) * t30 - g(2) * t31;
t13 = t15 * t42;
t12 = -t30 * t65 + t31 * t46;
t10 = -t30 * t46 - t31 * t65;
t8 = -t30 * t17 - t31 * t34;
t7 = -t31 * t17 + t30 * t34;
t6 = -g(1) * t10 - g(2) * t12;
t5 = g(1) * t49 - g(2) * t11;
t4 = g(1) * t53 + g(2) * t52;
t3 = -g(1) * t7 - g(2) * t8;
t2 = -g(1) * t52 + g(2) * t53 - g(3) * (t26 - t28 / 0.2e1);
t1 = -g(1) * t8 + g(2) * t7 - g(3) * (t61 + t54 / 0.2e1);
t9 = [0, 0, 0, 0, 0, 0, t55, g(1) * t47 + g(2) * t45, 0, 0, 0, 0, 0, 0, 0, 0, t19, t20, 0, t55 * pkin(1), 0, 0, 0, 0, 0, 0, t14, t15, 0, -g(1) * (-t66 - t70) - g(2) * (t29 + t39), 0, 0, 0, 0, 0, 0, t6, t5, -t13, -g(1) * (t50 - t66) - g(2) * (t39 + t60), 0, 0, 0, 0, 0, 0, t4, t3, -t13, -g(1) * (t48 - t66) - g(2) * (t39 + t56); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t20, 0, 0, 0, 0, 0, 0, 0, 0, t14, t15, 0, t19 * pkin(2), 0, 0, 0, 0, 0, 0, t6, t5, -t13, -g(1) * t50 - g(2) * t60, 0, 0, 0, 0, 0, 0, t4, t3, -t13, -g(1) * t48 - g(2) * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, 0, 0, 0, t6, t5, -t13, -g(1) * t59 - g(2) * t63, 0, 0, 0, 0, 0, 0, t4, t3, -t13, -g(1) * t51 - g(2) * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, g(1) * t12 - g(2) * t10 + t44 * t67, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t72 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg = t9;
