% Calculate inertial parameters regressor of gravitation load for
% S5RRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR10_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR10_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR10_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRPRR10_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:26:39
% EndTime: 2019-12-31 20:26:41
% DurationCPUTime: 0.60s
% Computational Cost: add. (481->110), mult. (1233->184), div. (0->0), fcn. (1563->12), ass. (0->70)
t46 = sin(pkin(10));
t51 = sin(qJ(2));
t55 = cos(qJ(2));
t79 = cos(pkin(10));
t36 = -t55 * t46 - t51 * t79;
t52 = sin(qJ(1));
t56 = cos(qJ(1));
t48 = cos(pkin(5));
t66 = -t51 * t46 + t55 * t79;
t59 = t66 * t48;
t15 = t52 * t36 + t56 * t59;
t49 = sin(qJ(5));
t28 = t36 * t48;
t16 = -t56 * t28 + t52 * t66;
t50 = sin(qJ(4));
t54 = cos(qJ(4));
t47 = sin(pkin(5));
t87 = t47 * t56;
t5 = t16 * t54 - t50 * t87;
t53 = cos(qJ(5));
t95 = t15 * t53 + t5 * t49;
t94 = -t15 * t49 + t5 * t53;
t81 = t56 * t51;
t83 = t52 * t55;
t32 = -t48 * t83 - t81;
t88 = t47 * t55;
t93 = -g(1) * t32 - g(3) * t88;
t89 = t47 * t52;
t86 = t49 * t54;
t84 = t52 * t51;
t82 = t53 * t54;
t80 = t56 * t55;
t77 = t48 * t80;
t76 = -t16 * t50 - t54 * t87;
t29 = t48 * t51 * pkin(2) + (-pkin(7) - qJ(3)) * t47;
t45 = t55 * pkin(2) + pkin(1);
t75 = -t52 * t29 + t56 * t45;
t26 = t66 * t47;
t27 = t36 * t47;
t73 = pkin(2) * t88 + t26 * pkin(3) - t27 * pkin(8);
t17 = -t52 * t28 - t56 * t66;
t8 = -t17 * t50 - t54 * t89;
t72 = g(1) * t76 + g(2) * t8;
t71 = pkin(4) * t54 + pkin(9) * t50;
t18 = t56 * t36 - t52 * t59;
t70 = g(1) * t15 - g(2) * t18;
t69 = g(1) * t56 + g(2) * t52;
t68 = g(1) * t52 - g(2) * t56;
t67 = -t56 * t29 - t52 * t45;
t20 = t27 * t50 + t48 * t54;
t65 = g(1) * t8 - g(2) * t76 - g(3) * t20;
t21 = -t27 * t54 + t48 * t50;
t9 = -t17 * t54 + t50 * t89;
t64 = g(1) * t9 + g(2) * t5 + g(3) * t21;
t63 = -pkin(3) * t17 - t18 * pkin(8) + t75;
t62 = g(1) * t17 - g(2) * t16 + g(3) * t27;
t61 = g(1) * t18 + g(2) * t15 + g(3) * t26;
t37 = pkin(2) * t77;
t60 = -pkin(2) * t84 + t15 * pkin(3) + pkin(8) * t16 + t37;
t58 = -t16 * pkin(3) + t15 * pkin(8) + t67;
t57 = pkin(2) * t32 + t18 * pkin(3) - t17 * pkin(8);
t34 = t69 * t47;
t33 = -t48 * t84 + t80;
t31 = -t48 * t81 - t83;
t30 = -t77 + t84;
t25 = -g(3) * t48 - t47 * t68;
t3 = -t18 * t49 + t9 * t53;
t2 = -t18 * t53 - t9 * t49;
t1 = t61 * t50;
t4 = [0, 0, 0, 0, 0, 0, t68, t69, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t31 - g(2) * t33, -g(1) * t30 - g(2) * t32, -t34, -g(1) * (-t52 * pkin(1) + pkin(7) * t87) - g(2) * (t56 * pkin(1) + pkin(7) * t89), 0, 0, 0, 0, 0, 0, g(1) * t16 + g(2) * t17, t70, -t34, -g(1) * t67 - g(2) * t75, 0, 0, 0, 0, 0, 0, g(1) * t5 - g(2) * t9, t72, -t70, -g(1) * t58 - g(2) * t63, 0, 0, 0, 0, 0, 0, g(1) * t94 - g(2) * t3, -g(1) * t95 - g(2) * t2, -t72, -g(1) * (-pkin(4) * t5 + pkin(9) * t76 + t58) - g(2) * (t9 * pkin(4) + t8 * pkin(9) + t63); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t30 + t93, g(3) * t47 * t51 + g(1) * t33 - g(2) * t31, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t62, 0, -g(2) * t37 + (g(2) * t84 + t93) * pkin(2), 0, 0, 0, 0, 0, 0, -t61 * t54, t1, t62, -g(1) * t57 - g(2) * t60 - g(3) * t73, 0, 0, 0, 0, 0, 0, -g(1) * (-t17 * t49 + t18 * t82) - g(2) * (t15 * t82 + t16 * t49) - g(3) * (t26 * t82 - t27 * t49), -g(1) * (-t17 * t53 - t18 * t86) - g(2) * (-t15 * t86 + t16 * t53) - g(3) * (-t26 * t86 - t27 * t53), -t1, -g(1) * (t18 * t71 + t57) - g(2) * (t15 * t71 + t60) - g(3) * (t26 * t71 + t73); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t64, 0, 0, 0, 0, 0, 0, 0, 0, t65 * t53, -t65 * t49, -t64, -g(1) * (-t8 * pkin(4) + t9 * pkin(9)) - g(2) * (pkin(4) * t76 + t5 * pkin(9)) - g(3) * (t20 * pkin(4) + t21 * pkin(9)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t2 + g(2) * t95 - g(3) * (-t21 * t49 - t26 * t53), g(1) * t3 + g(2) * t94 - g(3) * (-t21 * t53 + t26 * t49), 0, 0;];
taug_reg = t4;
