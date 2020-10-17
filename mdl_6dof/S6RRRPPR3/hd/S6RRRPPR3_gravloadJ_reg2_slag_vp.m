% Calculate inertial parameters regressor of gravitation load for
% S6RRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPPR3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPPR3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 04:32:43
% EndTime: 2019-05-07 04:32:43
% DurationCPUTime: 0.37s
% Computational Cost: add. (364->98), mult. (438->111), div. (0->0), fcn. (415->8), ass. (0->67)
t37 = cos(qJ(1));
t34 = sin(qJ(1));
t67 = g(2) * t34;
t12 = g(1) * t37 + t67;
t31 = qJ(2) + qJ(3);
t28 = sin(t31);
t73 = t12 * t28;
t22 = t28 * pkin(5);
t29 = cos(t31);
t72 = t29 * pkin(9) + t22;
t4 = g(3) * t28 + t12 * t29;
t71 = pkin(3) + pkin(4);
t33 = sin(qJ(2));
t70 = pkin(2) * t33;
t69 = pkin(3) * t28;
t65 = g(3) * t29;
t25 = t29 * pkin(3);
t24 = t29 * pkin(4);
t64 = t29 * t37;
t32 = sin(qJ(6));
t63 = t34 * t32;
t35 = cos(qJ(6));
t62 = t34 * t35;
t61 = t37 * t32;
t60 = t37 * t35;
t38 = -pkin(8) - pkin(7);
t59 = t37 * t38;
t21 = t28 * qJ(4);
t58 = t25 + t21;
t57 = qJ(4) * t29;
t56 = qJ(5) + t38;
t55 = pkin(9) + t71;
t36 = cos(qJ(2));
t30 = t36 * pkin(2);
t27 = t30 + pkin(1);
t16 = t37 * t27;
t54 = pkin(3) * t64 + t37 * t21 + t16;
t53 = t24 + t58;
t52 = t30 + t58;
t51 = g(1) * t56;
t50 = g(2) * t56;
t49 = -t27 - t21;
t48 = g(1) * t55;
t13 = t34 * t57;
t47 = -t34 * t70 + t13;
t15 = t37 * t57;
t46 = -t37 * t70 + t15;
t45 = t53 + t72;
t44 = -t69 - t70;
t11 = g(1) * t34 - g(2) * t37;
t43 = g(2) * (pkin(4) * t64 + t54);
t42 = t49 - t25;
t41 = -g(3) * t36 + t12 * t33;
t40 = t71 * t73;
t39 = (t37 * t48 + t55 * t67) * t28;
t18 = pkin(5) * t64;
t17 = t34 * t29 * pkin(5);
t10 = t28 * t60 - t63;
t9 = -t28 * t61 - t62;
t8 = -t28 * t62 - t61;
t7 = t28 * t63 - t60;
t6 = t11 * t29;
t5 = t11 * t28;
t3 = -t65 + t73;
t2 = t4 * t35;
t1 = t4 * t32;
t14 = [0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t36, -t11 * t33, -t12, -g(1) * (-t34 * pkin(1) + t37 * pkin(7)) - g(2) * (t37 * pkin(1) + t34 * pkin(7)) 0, 0, 0, 0, 0, 0, t6, -t5, -t12, -g(1) * (-t34 * t27 - t59) - g(2) * (-t34 * t38 + t16) 0, 0, 0, 0, 0, 0, t6, -t12, t5, g(1) * t59 - g(2) * t54 + (-g(1) * t42 + g(2) * t38) * t34, 0, 0, 0, 0, 0, 0, t5, -t6, t12, -t43 + t37 * t51 + (-g(1) * (t42 - t24) + t50) * t34, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t10, -g(1) * t7 - g(2) * t9, t6, -t43 + (-g(2) * t72 + t51) * t37 + (-g(1) * (t49 - t22) + t50 + t29 * t48) * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, g(3) * t33 + t12 * t36, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t41 * pkin(2), 0, 0, 0, 0, 0, 0, t3, 0, -t4, -g(1) * (t44 * t37 + t15) - g(2) * (t44 * t34 + t13) - g(3) * t52, 0, 0, 0, 0, 0, 0, -t4, -t3, 0, -g(1) * t46 - g(2) * t47 - g(3) * (t24 + t52) + t40, 0, 0, 0, 0, 0, 0, -t2, t1, t3, -g(1) * (t18 + t46) - g(2) * (t17 + t47) - g(3) * (t30 + t45) + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, -t4, -g(1) * (-t37 * t69 + t15) - g(2) * (-t34 * t69 + t13) - g(3) * t58, 0, 0, 0, 0, 0, 0, -t4, -t3, 0, -g(1) * t15 - g(2) * t13 - g(3) * t53 + t40, 0, 0, 0, 0, 0, 0, -t2, t1, t3, -g(1) * (t15 + t18) - g(2) * (t13 + t17) - g(3) * t45 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t9 + g(2) * t7 - t32 * t65, g(1) * t10 - g(2) * t8 - t35 * t65, 0, 0;];
taug_reg  = t14;
