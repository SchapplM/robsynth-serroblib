% Calculate inertial parameters regressor of gravitation load for
% S6RRPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPPRR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR5_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR5_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 10:42:34
% EndTime: 2019-05-06 10:42:36
% DurationCPUTime: 0.63s
% Computational Cost: add. (387->119), mult. (963->176), div. (0->0), fcn. (1149->10), ass. (0->72)
t51 = sin(qJ(2));
t52 = sin(qJ(1));
t55 = cos(qJ(2));
t56 = cos(qJ(1));
t79 = cos(pkin(6));
t74 = t56 * t79;
t33 = t51 * t74 + t52 * t55;
t50 = sin(qJ(5));
t54 = cos(qJ(5));
t48 = sin(pkin(6));
t89 = t48 * t56;
t12 = t33 * t54 + t50 * t89;
t32 = t52 * t51 - t55 * t74;
t49 = sin(qJ(6));
t53 = cos(qJ(6));
t98 = t12 * t49 + t32 * t53;
t97 = t12 * t53 - t32 * t49;
t96 = g(3) * t48;
t93 = t48 * t51;
t92 = t48 * t52;
t91 = t48 * t54;
t90 = t48 * t55;
t88 = t49 * t54;
t87 = t53 * t54;
t86 = t54 * t55;
t85 = qJ(3) - pkin(9);
t84 = pkin(2) * t90 + qJ(3) * t93;
t83 = t56 * pkin(1) + pkin(8) * t92;
t82 = qJ(4) * t48;
t81 = t32 * qJ(3);
t75 = t52 * t79;
t34 = t56 * t51 + t55 * t75;
t80 = t34 * qJ(3);
t35 = -t51 * t75 + t56 * t55;
t78 = t35 * pkin(2) + t83;
t77 = pkin(3) * t90 + t84;
t76 = -t52 * pkin(1) + pkin(8) * t89;
t20 = t32 * pkin(2);
t73 = t33 * qJ(3) - t20;
t26 = t34 * pkin(2);
t72 = t35 * qJ(3) - t26;
t71 = pkin(4) * t90 + t77;
t70 = -t33 * pkin(2) + t76;
t69 = pkin(5) * t54 + pkin(10) * t50;
t15 = t35 * t50 + t52 * t91;
t65 = -t33 * t50 + t54 * t89;
t68 = g(1) * t65 + g(2) * t15;
t9 = g(1) * t32 - g(2) * t34;
t67 = g(1) * t56 + g(2) * t52;
t66 = g(1) * t52 - g(2) * t56;
t19 = t32 * pkin(3);
t64 = -t32 * pkin(4) + t85 * t33 - t19 - t20;
t25 = t34 * pkin(3);
t63 = -t34 * pkin(4) + t85 * t35 - t25 - t26;
t30 = -t50 * t93 - t79 * t54;
t62 = g(1) * t15 - g(2) * t65 - g(3) * t30;
t16 = t35 * t54 - t50 * t92;
t31 = -t79 * t50 + t51 * t91;
t61 = g(1) * t16 + g(2) * t12 + g(3) * t31;
t60 = t35 * pkin(3) - t52 * t82 + t78;
t4 = -g(1) * t34 - g(2) * t32 + g(3) * t90;
t7 = g(1) * t35 + g(2) * t33 + g(3) * t93;
t59 = -t33 * pkin(3) - t56 * t82 + t70;
t58 = t35 * pkin(4) + t85 * t34 + t60;
t57 = -t33 * pkin(4) - t85 * t32 + t59;
t36 = t67 * t48;
t17 = g(3) * t79 + t66 * t48;
t10 = g(1) * t33 - g(2) * t35;
t3 = t4 * t50;
t2 = t16 * t53 - t34 * t49;
t1 = -t16 * t49 - t34 * t53;
t5 = [0, 0, 0, 0, 0, 0, t66, t67, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, -t36, -g(1) * t76 - g(2) * t83, 0, 0, 0, 0, 0, 0, t10, -t36, t9, -g(1) * (t70 - t81) - g(2) * (t78 + t80) 0, 0, 0, 0, 0, 0, t10, t9, t36, -g(1) * (t59 - t81) - g(2) * (t60 + t80) 0, 0, 0, 0, 0, 0, g(1) * t12 - g(2) * t16, t68, -t9, -g(1) * t57 - g(2) * t58, 0, 0, 0, 0, 0, 0, g(1) * t97 - g(2) * t2, -g(1) * t98 - g(2) * t1, -t68, -g(1) * (-pkin(5) * t12 + pkin(10) * t65 + t57) - g(2) * (t16 * pkin(5) + t15 * pkin(10) + t58); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t4, 0, -t7, -g(1) * t72 - g(2) * t73 - g(3) * t84, 0, 0, 0, 0, 0, 0, -t4, -t7, 0, -g(1) * (-t25 + t72) - g(2) * (-t19 + t73) - g(3) * t77, 0, 0, 0, 0, 0, 0, -t4 * t54, t3, t7, -g(1) * t63 - g(2) * t64 - g(3) * (-pkin(9) * t93 + t71) 0, 0, 0, 0, 0, 0, -g(1) * (-t34 * t87 - t35 * t49) - g(2) * (-t32 * t87 - t33 * t49) - (-t49 * t51 + t53 * t86) * t96, -g(1) * (t34 * t88 - t35 * t53) - g(2) * (t32 * t88 - t33 * t53) - (-t49 * t86 - t51 * t53) * t96, -t3, -g(1) * (-t69 * t34 + t63) - g(2) * (-t69 * t32 + t64) - g(3) * t71 - (-pkin(9) * t51 + t69 * t55) * t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t61, 0, 0, 0, 0, 0, 0, 0, 0, t62 * t53, -t62 * t49, -t61, -g(1) * (-t15 * pkin(5) + t16 * pkin(10)) - g(2) * (pkin(5) * t65 + t12 * pkin(10)) - g(3) * (t30 * pkin(5) + t31 * pkin(10)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 + g(2) * t98 - g(3) * (-t31 * t49 + t53 * t90) g(1) * t2 + g(2) * t97 - g(3) * (-t31 * t53 - t49 * t90) 0, 0;];
taug_reg  = t5;
