% Calculate minimal parameter regressor of gravitation load for
% S6PRPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% taug_reg [6x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRPRRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:38:15
% EndTime: 2019-05-04 23:38:16
% DurationCPUTime: 0.38s
% Computational Cost: add. (496->85), mult. (1312->138), div. (0->0), fcn. (1705->12), ass. (0->61)
t47 = sin(pkin(11));
t54 = sin(qJ(2));
t57 = cos(qJ(2));
t70 = cos(pkin(11));
t64 = -t54 * t47 + t57 * t70;
t48 = sin(pkin(10));
t50 = cos(pkin(10));
t51 = cos(pkin(6));
t74 = t51 * t57;
t88 = -t48 * t74 - t50 * t54;
t39 = t64 * t51;
t42 = -t57 * t47 - t54 * t70;
t27 = t50 * t39 + t48 * t42;
t30 = -t48 * t39 + t50 * t42;
t49 = sin(pkin(6));
t37 = t64 * t49;
t87 = -g(1) * t30 - g(2) * t27 - g(3) * t37;
t40 = t42 * t51;
t28 = -t50 * t40 + t48 * t64;
t29 = -t48 * t40 - t50 * t64;
t38 = t42 * t49;
t53 = sin(qJ(4));
t56 = cos(qJ(4));
t78 = t49 * t56;
t61 = g(3) * (t38 * t53 + t51 * t56) + g(2) * (-t28 * t53 - t50 * t78) + g(1) * (t29 * t53 + t48 * t78);
t80 = t48 * t54;
t79 = t49 * t53;
t77 = t49 * t57;
t75 = t51 * t54;
t52 = sin(qJ(5));
t73 = t52 * t56;
t55 = cos(qJ(5));
t71 = t55 * t56;
t68 = t50 * t74;
t33 = -t38 * t56 + t51 * t53;
t14 = t33 * t52 + t37 * t55;
t17 = t28 * t56 - t50 * t79;
t6 = t17 * t52 + t27 * t55;
t19 = -t29 * t56 + t48 * t79;
t8 = t19 * t52 + t30 * t55;
t1 = g(1) * t8 + g(2) * t6 + g(3) * t14;
t15 = t33 * t55 - t37 * t52;
t7 = t17 * t55 - t27 * t52;
t9 = t19 * t55 - t30 * t52;
t63 = g(1) * t9 + g(2) * t7 + g(3) * t15;
t10 = t27 * t73 - t28 * t55;
t12 = t29 * t55 + t30 * t73;
t20 = t37 * t73 + t38 * t55;
t62 = g(1) * t12 + g(2) * t10 + g(3) * t20;
t60 = g(1) * t19 + g(2) * t17 + g(3) * t33;
t58 = -g(1) * t88 - g(3) * t77;
t43 = pkin(2) * t68;
t36 = -g(3) * t51 + (-g(1) * t48 + g(2) * t50) * t49;
t21 = t37 * t71 - t38 * t52;
t13 = -t29 * t52 + t30 * t71;
t11 = t27 * t71 + t28 * t52;
t5 = t87 * t53;
t4 = t61 * t55;
t3 = t61 * t52;
t2 = -g(1) * t13 - g(2) * t11 - g(3) * t21;
t16 = [-g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -g(2) * (t68 - t80) + t58, -g(1) * (t48 * t75 - t50 * t57) - g(2) * (-t48 * t57 - t50 * t75) + g(3) * t49 * t54, -g(2) * t43 + (g(2) * t80 + t58) * pkin(2), 0, 0, 0, 0, 0, t87 * t56, -t5, 0, 0, 0, 0, 0, t2, t62, t2, t5, -t62, -g(1) * (t88 * pkin(2) + t13 * pkin(5) - t29 * pkin(8) + t12 * qJ(6)) - g(2) * (-pkin(2) * t80 + t11 * pkin(5) + pkin(8) * t28 + t10 * qJ(6) + t43) - g(3) * (pkin(2) * t77 + t21 * pkin(5) - t38 * pkin(8) + t20 * qJ(6)) + t87 * (pkin(4) * t56 + pkin(9) * t53 + pkin(3)); 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t60, 0, 0, 0, 0, 0, -t4, t3, -t4, -t60, -t3, -t60 * pkin(9) - t61 * (pkin(5) * t55 + qJ(6) * t52 + pkin(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t63, t1, 0, -t63, -g(1) * (-t8 * pkin(5) + t9 * qJ(6)) - g(2) * (-t6 * pkin(5) + t7 * qJ(6)) - g(3) * (-t14 * pkin(5) + t15 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t16;
