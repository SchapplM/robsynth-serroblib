% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPR8
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
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPR8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 23:10:28
% EndTime: 2019-05-05 23:10:30
% DurationCPUTime: 0.44s
% Computational Cost: add. (294->94), mult. (412->136), div. (0->0), fcn. (415->10), ass. (0->66)
t42 = sin(qJ(1));
t45 = cos(qJ(1));
t82 = -g(1) * t42 + g(2) * t45;
t41 = sin(qJ(3));
t43 = cos(qJ(4));
t57 = t45 * t43;
t40 = sin(qJ(4));
t66 = t42 * t40;
t13 = -t41 * t66 + t57;
t58 = t45 * t40;
t65 = t42 * t43;
t15 = t41 * t58 + t65;
t44 = cos(qJ(3));
t74 = g(3) * t44;
t81 = -g(1) * t13 - g(2) * t15 + t40 * t74;
t12 = -g(3) * t41 - t44 * t82;
t80 = -pkin(1) - pkin(7);
t73 = t40 * pkin(4);
t72 = t41 * t42;
t71 = t41 * t45;
t38 = qJ(4) + pkin(10);
t30 = qJ(6) + t38;
t25 = sin(t30);
t70 = t42 * t25;
t26 = cos(t30);
t69 = t42 * t26;
t28 = sin(t38);
t68 = t42 * t28;
t29 = cos(t38);
t67 = t42 * t29;
t64 = t44 * t45;
t19 = pkin(5) * t28 + t73;
t63 = t45 * t19;
t62 = t45 * t25;
t61 = t45 * t26;
t60 = t45 * t28;
t59 = t45 * t29;
t39 = -qJ(5) - pkin(8);
t33 = t43 * pkin(4);
t20 = pkin(5) * t29 + t33;
t56 = t45 * pkin(1) + t42 * qJ(2);
t54 = t45 * pkin(7) + t56;
t53 = g(2) * t54;
t51 = t41 * pkin(3) - t44 * pkin(8);
t22 = g(1) * t45 + g(2) * t42;
t18 = pkin(3) + t20;
t37 = -pkin(9) + t39;
t49 = t41 * t18 + t44 * t37;
t27 = t33 + pkin(3);
t47 = t41 * t27 + t44 * t39;
t32 = t45 * qJ(2);
t17 = t22 * t44;
t16 = t41 * t57 - t66;
t14 = t41 * t65 + t58;
t11 = g(1) * t72 - g(2) * t71 + t74;
t10 = t41 * t59 - t68;
t9 = t41 * t60 + t67;
t8 = t41 * t67 + t60;
t7 = -t41 * t68 + t59;
t6 = t41 * t61 - t70;
t5 = t41 * t62 + t69;
t4 = t41 * t69 + t62;
t3 = -t41 * t70 + t61;
t2 = g(1) * t4 - g(2) * t6 + t26 * t74;
t1 = -g(1) * t3 - g(2) * t5 + t25 * t74;
t21 = [0, 0, 0, 0, 0, 0, -t82, t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, -t22, -g(1) * (-t42 * pkin(1) + t32) - g(2) * t56, 0, 0, 0, 0, 0, 0, -t22 * t41, -t17, -t82, -g(1) * (t80 * t42 + t32) - t53, 0, 0, 0, 0, 0, 0, -g(1) * t16 - g(2) * t14, g(1) * t15 - g(2) * t13, t17, -g(1) * (pkin(3) * t71 - pkin(8) * t64 + t32) - t53 + (-g(1) * t80 - g(2) * t51) * t42, 0, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7, t17, -g(1) * (t27 * t71 + t39 * t64 + t32) - g(2) * (pkin(4) * t58 + t54) + (-g(1) * (-t73 + t80) - g(2) * t47) * t42, 0, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3, t17, -g(1) * (t18 * t71 + t37 * t64 + t32) - g(2) * (t54 + t63) + (-g(1) * (-t19 + t80) - g(2) * t49) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, t11, 0, 0, 0, 0, 0, 0, 0, 0, -t12 * t43, t12 * t40, -t11, g(3) * t51 + t82 * (pkin(3) * t44 + pkin(8) * t41) 0, 0, 0, 0, 0, 0, -t12 * t29, t12 * t28, -t11, g(3) * t47 + t82 * (t27 * t44 - t39 * t41) 0, 0, 0, 0, 0, 0, -t12 * t26, t12 * t25, -t11, g(3) * t49 + t82 * (t18 * t44 - t37 * t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, g(1) * t14 - g(2) * t16 + t43 * t74, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t28 * t74, g(1) * t8 - g(2) * t10 + t29 * t74, 0, t81 * pkin(4), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(1) * (-t19 * t72 + t45 * t20) - g(2) * (t42 * t20 + t41 * t63) + t19 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0;];
taug_reg  = t21;
