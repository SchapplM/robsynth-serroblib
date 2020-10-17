% Calculate inertial parameters regressor of gravitation load for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:57:01
% EndTime: 2019-05-04 23:57:03
% DurationCPUTime: 0.44s
% Computational Cost: add. (335->93), mult. (852->135), div. (0->0), fcn. (1029->10), ass. (0->62)
t32 = sin(pkin(10));
t34 = cos(pkin(10));
t41 = cos(qJ(2));
t38 = sin(qJ(2));
t57 = cos(pkin(6));
t53 = t38 * t57;
t20 = t32 * t41 + t34 * t53;
t22 = -t32 * t53 + t34 * t41;
t75 = -g(1) * t22 - g(2) * t20;
t33 = sin(pkin(6));
t72 = g(3) * t33;
t52 = t41 * t57;
t19 = t32 * t38 - t34 * t52;
t71 = t19 * pkin(8);
t21 = t32 * t52 + t34 * t38;
t70 = t21 * pkin(8);
t36 = sin(qJ(5));
t69 = t19 * t36;
t68 = t21 * t36;
t37 = sin(qJ(4));
t67 = t33 * t37;
t66 = t33 * t38;
t40 = cos(qJ(4));
t65 = t33 * t40;
t64 = t33 * t41;
t63 = t36 * t37;
t62 = t36 * t38;
t61 = t36 * t41;
t39 = cos(qJ(5));
t60 = t37 * t39;
t59 = t38 * t39;
t58 = pkin(2) * t64 + qJ(3) * t66;
t56 = pkin(8) * t64 + t58;
t17 = t19 * pkin(2);
t55 = -t17 - t71;
t18 = t21 * pkin(2);
t54 = -t18 - t70;
t51 = t20 * qJ(3) - t17;
t50 = t22 * qJ(3) - t18;
t49 = g(3) * t56;
t48 = pkin(4) * t37 - pkin(9) * t40;
t31 = pkin(5) * t39 + pkin(4);
t35 = -qJ(6) - pkin(9);
t47 = t31 * t37 + t35 * t40;
t11 = -t21 * t40 + t32 * t67;
t13 = t19 * t40 + t34 * t67;
t23 = t37 * t57 + t40 * t64;
t44 = g(1) * t11 - g(2) * t13 + g(3) * t23;
t12 = t21 * t37 + t32 * t65;
t14 = -t19 * t37 + t34 * t65;
t24 = -t37 * t64 + t40 * t57;
t43 = g(1) * t12 - g(2) * t14 + g(3) * t24;
t9 = -g(1) * t21 - g(2) * t19 + g(3) * t64;
t42 = g(3) * t66 - t75;
t1 = -g(1) * (-t12 * t36 + t22 * t39) - g(2) * (t14 * t36 + t20 * t39) - g(3) * (-t24 * t36 + t33 * t59);
t8 = t42 * t40;
t6 = t44 * t39;
t5 = t44 * t36;
t4 = -g(1) * (t22 * t60 - t68) - g(2) * (t20 * t60 - t69) - (t37 * t59 + t61) * t72;
t3 = -g(1) * (-t21 * t39 - t22 * t63) - g(2) * (-t19 * t39 - t20 * t63) - (-t37 * t62 + t39 * t41) * t72;
t2 = -g(1) * (-t12 * t39 - t22 * t36) - g(2) * (t14 * t39 - t20 * t36) - g(3) * (-t24 * t39 - t33 * t62);
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t42, -g(1) * t50 - g(2) * t51 - g(3) * t58, 0, 0, 0, 0, 0, 0, -t42 * t37, -t8, -t9, -g(1) * (t50 - t70) - g(2) * (t51 - t71) - t49, 0, 0, 0, 0, 0, 0, t4, t3, t8, -g(1) * t54 - g(2) * t55 - g(3) * (t48 * t66 + t56) + t75 * (qJ(3) + t48) 0, 0, 0, 0, 0, 0, t4, t3, t8, -g(1) * (-pkin(5) * t68 + t54) - g(2) * (-pkin(5) * t69 + t55) - t49 - (pkin(5) * t61 + t38 * t47) * t72 + t75 * (qJ(3) + t47); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t43, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t43, -g(1) * (-pkin(4) * t11 + pkin(9) * t12) - g(2) * (pkin(4) * t13 - pkin(9) * t14) - g(3) * (-pkin(4) * t23 + pkin(9) * t24) 0, 0, 0, 0, 0, 0, t6, -t5, -t43, -g(1) * (-t11 * t31 - t12 * t35) - g(2) * (t13 * t31 + t14 * t35) - g(3) * (-t23 * t31 - t24 * t35); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44;];
taug_reg  = t7;
