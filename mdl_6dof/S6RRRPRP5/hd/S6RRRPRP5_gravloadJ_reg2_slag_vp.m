% Calculate inertial parameters regressor of gravitation load for
% S6RRRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRRPRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRRPRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP5_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 07:55:18
% EndTime: 2019-05-07 07:55:20
% DurationCPUTime: 0.56s
% Computational Cost: add. (547->114), mult. (618->159), div. (0->0), fcn. (638->10), ass. (0->73)
t52 = qJ(3) + pkin(10);
t44 = cos(t52);
t57 = cos(qJ(3));
t47 = t57 * pkin(3);
t33 = pkin(4) * t44 + t47;
t31 = pkin(2) + t33;
t58 = cos(qJ(2));
t22 = t58 * t31;
t53 = -qJ(4) - pkin(8);
t51 = -pkin(9) + t53;
t55 = sin(qJ(2));
t99 = -t55 * t51 + t22;
t56 = sin(qJ(1));
t59 = cos(qJ(1));
t38 = g(1) * t59 + g(2) * t56;
t54 = sin(qJ(3));
t78 = t59 * t57;
t85 = t56 * t58;
t23 = t54 * t85 + t78;
t79 = t59 * t54;
t25 = t56 * t57 - t58 * t79;
t92 = g(3) * t55;
t98 = -g(1) * t25 + g(2) * t23 + t54 * t92;
t17 = -g(3) * t58 + t38 * t55;
t96 = g(1) * t56;
t90 = t54 * pkin(3);
t45 = qJ(5) + t52;
t40 = sin(t45);
t89 = t40 * t55;
t41 = cos(t45);
t88 = t41 * t55;
t86 = t55 * t59;
t84 = t58 * t59;
t43 = sin(t52);
t32 = pkin(4) * t43 + t90;
t29 = t59 * t32;
t83 = t59 * t40;
t82 = t59 * t41;
t81 = t59 * t43;
t80 = t59 * t44;
t77 = -t58 * t29 + t56 * t33;
t76 = t59 * pkin(1) + t56 * pkin(7);
t10 = t41 * t85 - t83;
t9 = t40 * t85 + t82;
t74 = -t9 * pkin(5) + t10 * qJ(6);
t11 = -t56 * t41 + t58 * t83;
t12 = t56 * t40 + t58 * t82;
t73 = -t11 * pkin(5) + t12 * qJ(6);
t72 = -t32 * t85 - t59 * t33;
t71 = g(1) * t9 - g(2) * t11;
t70 = t58 * pkin(2) + t55 * pkin(8);
t68 = -g(2) * t59 + t96;
t67 = pkin(5) * t41 + qJ(6) * t40;
t42 = t47 + pkin(2);
t65 = t58 * t42 - t55 * t53;
t62 = t31 * t84 + t56 * t32 - t51 * t86 + t76;
t1 = g(1) * t11 + g(2) * t9 + g(3) * t89;
t3 = g(1) * t12 + g(2) * t10 + g(3) * t88;
t48 = t59 * pkin(7);
t61 = t29 + t48 + (-pkin(1) - t99) * t56;
t34 = qJ(6) * t88;
t30 = t68 * t55;
t26 = t56 * t54 + t58 * t78;
t24 = -t57 * t85 + t79;
t18 = t38 * t58 + t92;
t16 = t56 * t43 + t58 * t80;
t15 = t56 * t44 - t58 * t81;
t14 = -t44 * t85 + t81;
t13 = t43 * t85 + t80;
t6 = t17 * t41;
t5 = t17 * t40;
t4 = g(1) * t10 - g(2) * t12;
t2 = [0, 0, 0, 0, 0, 0, t68, t38, 0, 0, 0, 0, 0, 0, 0, 0, t68 * t58, -t30, -t38, -g(1) * (-t56 * pkin(1) + t48) - g(2) * t76, 0, 0, 0, 0, 0, 0, -g(1) * t24 - g(2) * t26, -g(1) * t23 - g(2) * t25, t30, -g(1) * t48 - g(2) * (t70 * t59 + t76) - (-pkin(1) - t70) * t96, 0, 0, 0, 0, 0, 0, -g(1) * t14 - g(2) * t16, -g(1) * t13 - g(2) * t15, t30, -g(1) * (pkin(3) * t79 + t48) - g(2) * (t42 * t84 - t53 * t86 + t76) + (-g(1) * (-pkin(1) - t65) - g(2) * t90) * t56, 0, 0, 0, 0, 0, 0, t4, -t71, t30, -g(1) * t61 - g(2) * t62, 0, 0, 0, 0, 0, 0, t4, t30, t71, -g(1) * (-t10 * pkin(5) - t9 * qJ(6) + t61) - g(2) * (t12 * pkin(5) + t11 * qJ(6) + t62); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, t18, 0, 0, 0, 0, 0, 0, 0, 0, t17 * t57, -t17 * t54, -t18, -g(3) * t70 + t38 * (pkin(2) * t55 - pkin(8) * t58) 0, 0, 0, 0, 0, 0, t17 * t44, -t17 * t43, -t18, -g(3) * t65 + t38 * (t42 * t55 + t53 * t58) 0, 0, 0, 0, 0, 0, t6, -t5, -t18, -g(3) * t99 + t38 * (t31 * t55 + t51 * t58) 0, 0, 0, 0, 0, 0, t6, -t18, t5, -g(3) * t22 + (-g(3) * t67 + t38 * t51) * t58 + (g(3) * t51 + t38 * (t31 + t67)) * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, g(1) * t26 - g(2) * t24 + t57 * t92, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 + g(2) * t13 + t43 * t92, g(1) * t16 - g(2) * t14 + t44 * t92, 0, t98 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t3, 0, -g(1) * t77 - g(2) * t72 + t32 * t92, 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * (t73 + t77) - g(2) * (t72 + t74) - g(3) * (t34 + (-pkin(5) * t40 - t32) * t55); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * t73 - g(2) * t74 - g(3) * (-pkin(5) * t89 + t34); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
