% Calculate inertial parameters regressor of gravitation load for
% S6RRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR2_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR2_gravloadJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:47:03
% EndTime: 2019-05-06 09:47:04
% DurationCPUTime: 0.45s
% Computational Cost: add. (417->96), mult. (414->127), div. (0->0), fcn. (413->12), ass. (0->61)
t35 = sin(pkin(11));
t37 = -qJ(3) - pkin(7);
t51 = t35 * pkin(4) - t37;
t40 = sin(qJ(1));
t42 = cos(qJ(1));
t17 = g(1) * t42 + g(2) * t40;
t34 = qJ(2) + pkin(10);
t27 = cos(t34);
t33 = pkin(11) + qJ(5);
t24 = sin(t33);
t57 = t42 * t24;
t26 = cos(t33);
t62 = t40 * t26;
t11 = -t27 * t57 + t62;
t25 = sin(t34);
t68 = g(3) * t25;
t56 = t42 * t26;
t63 = t40 * t24;
t9 = t27 * t63 + t56;
t75 = -g(1) * t11 + g(2) * t9 + t24 * t68;
t3 = -g(3) * t27 + t17 * t25;
t39 = sin(qJ(2));
t73 = pkin(2) * t39;
t41 = cos(qJ(2));
t30 = t41 * pkin(2);
t23 = t30 + pkin(1);
t18 = t42 * t23;
t70 = g(2) * t18;
t36 = cos(pkin(11));
t22 = t36 * pkin(4) + pkin(3);
t28 = qJ(6) + t33;
t20 = sin(t28);
t65 = t40 * t20;
t21 = cos(t28);
t64 = t40 * t21;
t61 = t40 * t35;
t60 = t40 * t36;
t59 = t42 * t20;
t58 = t42 * t21;
t55 = t42 * t35;
t54 = t42 * t36;
t38 = -pkin(8) - qJ(4);
t53 = pkin(5) * t24 + t51;
t16 = g(1) * t40 - g(2) * t42;
t50 = t27 * pkin(3) + t25 * qJ(4);
t14 = pkin(5) * t26 + t22;
t32 = -pkin(9) + t38;
t49 = t27 * t14 - t25 * t32;
t48 = t27 * t22 - t25 * t38;
t43 = -g(3) * t41 + t17 * t39;
t13 = t16 * t25;
t12 = t27 * t56 + t63;
t10 = -t27 * t62 + t57;
t8 = t27 * t58 + t65;
t7 = -t27 * t59 + t64;
t6 = -t27 * t64 + t59;
t5 = t27 * t65 + t58;
t4 = t17 * t27 + t68;
t2 = g(1) * t8 - g(2) * t6 + t21 * t68;
t1 = -g(1) * t7 + g(2) * t5 + t20 * t68;
t15 = [0, 0, 0, 0, 0, 0, t16, t17, 0, 0, 0, 0, 0, 0, 0, 0, t16 * t41, -t16 * t39, -t17, -g(1) * (-t40 * pkin(1) + t42 * pkin(7)) - g(2) * (t42 * pkin(1) + t40 * pkin(7)) 0, 0, 0, 0, 0, 0, t16 * t27, -t13, -t17, -g(1) * (-t40 * t23 - t42 * t37) - g(2) * (-t40 * t37 + t18) 0, 0, 0, 0, 0, 0, -g(1) * (-t27 * t60 + t55) - g(2) * (t27 * t54 + t61) -g(1) * (t27 * t61 + t54) - g(2) * (-t27 * t55 + t60) t13, -t70 + (g(1) * t37 - g(2) * t50) * t42 + (-g(1) * (-t23 - t50) + g(2) * t37) * t40, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, t13, -t70 + (-g(1) * t51 - g(2) * t48) * t42 + (-g(1) * (-t23 - t48) - g(2) * t51) * t40, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t8, -g(1) * t5 - g(2) * t7, t13, -t70 + (-g(1) * t53 - g(2) * t49) * t42 + (-g(1) * (-t23 - t49) - g(2) * t53) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, g(3) * t39 + t17 * t41, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t43 * pkin(2), 0, 0, 0, 0, 0, 0, t3 * t36, -t3 * t35, -t4, -g(3) * (t30 + t50) + t17 * (pkin(3) * t25 - qJ(4) * t27 + t73) 0, 0, 0, 0, 0, 0, t3 * t26, -t3 * t24, -t4, -g(3) * (t30 + t48) + t17 * (t22 * t25 + t27 * t38 + t73) 0, 0, 0, 0, 0, 0, t3 * t21, -t3 * t20, -t4, -g(3) * (t30 + t49) + t17 * (t14 * t25 + t27 * t32 + t73); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, g(1) * t12 - g(2) * t10 + t26 * t68, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t75 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t15;
