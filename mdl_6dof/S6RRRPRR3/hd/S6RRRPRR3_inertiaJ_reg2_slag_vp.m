% Calculate inertial parameters regressor of joint inertia matrix for
% S6RRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 18:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRR3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR3_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRR3_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t69 = sin(qJ(6));
t63 = t69 ^ 2;
t73 = cos(qJ(6));
t66 = t73 ^ 2;
t99 = t63 + t66;
t93 = t99 * pkin(10);
t71 = sin(qJ(3));
t60 = t71 * pkin(2);
t51 = t60 + qJ(4);
t70 = sin(qJ(5));
t74 = cos(qJ(5));
t75 = cos(qJ(3));
t110 = t75 * pkin(2);
t54 = pkin(3) + t110;
t88 = -pkin(4) - t54;
t30 = t51 * t74 + t70 * t88;
t27 = -pkin(10) + t30;
t124 = t99 * t27;
t114 = -pkin(3) - pkin(4);
t45 = qJ(4) * t74 + t114 * t70;
t42 = -pkin(10) + t45;
t123 = t99 * t42;
t113 = -pkin(8) - pkin(7);
t72 = sin(qJ(2));
t94 = t113 * t72;
t76 = cos(qJ(2));
t95 = t113 * t76;
t24 = t71 * t94 - t75 * t95;
t35 = t71 * t72 - t75 * t76;
t13 = t35 * pkin(9) + t24;
t22 = -t71 * t95 - t75 * t94;
t38 = t71 * t76 + t72 * t75;
t82 = -t38 * pkin(9) + t22;
t7 = t13 * t70 - t74 * t82;
t122 = t7 ^ 2;
t18 = -t35 * t74 + t38 * t70;
t121 = t18 ^ 2;
t120 = t35 ^ 2;
t119 = 0.2e1 * t18;
t20 = t35 * t70 + t38 * t74;
t118 = 0.2e1 * t20;
t56 = -pkin(2) * t76 - pkin(1);
t117 = 0.2e1 * t56;
t116 = 0.2e1 * t73;
t115 = 0.2e1 * t76;
t112 = t7 * t69;
t111 = t7 * t74;
t47 = t74 * t88;
t28 = t51 * t70 - t47;
t26 = pkin(5) + t28;
t109 = t26 * t69;
t108 = t38 * t35;
t58 = t74 * t114;
t43 = qJ(4) * t70 - t58;
t41 = pkin(5) + t43;
t107 = t41 * t69;
t106 = t69 * t18;
t105 = t69 * t20;
t104 = t69 * t73;
t103 = t73 * t18;
t102 = t73 * t20;
t101 = t74 * t73;
t100 = t26 + t41;
t65 = t72 ^ 2;
t68 = t76 ^ 2;
t98 = t65 + t68;
t97 = -0.2e1 * t20 * t18;
t50 = -0.2e1 * t104;
t96 = t22 ^ 2 + t24 ^ 2;
t37 = t99 * t70;
t16 = pkin(3) * t35 - qJ(4) * t38 + t56;
t87 = -pkin(5) * t20 - pkin(10) * t18;
t11 = -pkin(4) * t35 - t16;
t5 = pkin(5) * t18 - pkin(10) * t20 + t11;
t9 = t74 * t13 + t70 * t82;
t3 = t5 * t73 - t69 * t9;
t4 = t5 * t69 + t73 * t9;
t1 = -t3 * t69 + t4 * t73;
t86 = -t18 * t27 + t20 * t26;
t85 = -t18 * t42 + t20 * t41;
t84 = -t18 * t70 - t20 * t74;
t83 = 0.2e1 * t22 * t38 - 0.2e1 * t24 * t35;
t78 = 0.2e1 * pkin(3);
t77 = 0.2e1 * qJ(4);
t67 = t74 ^ 2;
t64 = t70 ^ 2;
t62 = pkin(5) * t69;
t52 = t74 * t69;
t49 = 0.2e1 * t104;
t34 = t38 ^ 2;
t17 = t20 ^ 2;
t15 = t69 * t102;
t10 = (t63 - t66) * t20;
t6 = t7 * t73;
t2 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t65, t72 * t115, 0, t68, 0, 0, pkin(1) * t115, -0.2e1 * pkin(1) * t72, 0.2e1 * t98 * pkin(7), pkin(7) ^ 2 * t98 + pkin(1) ^ 2, t34, -0.2e1 * t108, 0, t120, 0, 0, t35 * t117, t38 * t117, t83, t56 ^ 2 + t96, t34, 0, 0.2e1 * t108, 0, 0, t120, 0.2e1 * t16 * t35, t83, -0.2e1 * t16 * t38, t16 ^ 2 + t96, t17, t97, 0, t121, 0, 0, t11 * t119, t11 * t118, -0.2e1 * t18 * t9 + 0.2e1 * t20 * t7, t11 ^ 2 + t9 ^ 2 + t122, t66 * t17, t17 * t50, t102 * t119, t63 * t17, t69 * t97, t121, 0.2e1 * t105 * t7 + 0.2e1 * t18 * t3, 0.2e1 * t102 * t7 - 0.2e1 * t18 * t4 (-t3 * t73 - t4 * t69) * t118, t3 ^ 2 + t4 ^ 2 + t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, 0, t76, 0, -t72 * pkin(7), -t76 * pkin(7), 0, 0, 0, 0, t38, 0, -t35, 0, -t22, -t24 (-t35 * t71 - t38 * t75) * pkin(2) (-t22 * t75 + t24 * t71) * pkin(2), 0, t38, 0, 0, t35, 0, -t22, -t35 * t51 - t38 * t54, t24, -t22 * t54 + t24 * t51, 0, 0, -t20, 0, t18, 0, t7, t9, -t18 * t30 + t20 * t28, t28 * t7 + t30 * t9, -t15, t10, -t106, t15, -t103, 0, t69 * t86 + t6, t73 * t86 - t112, -t1, t1 * t27 + t7 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t110, -0.2e1 * t60, 0 (t71 ^ 2 + t75 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t54, 0, 0.2e1 * t51, t51 ^ 2 + t54 ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t28, 0.2e1 * t30, 0, t28 ^ 2 + t30 ^ 2, t63, t49, 0, t66, 0, 0, t26 * t116, -0.2e1 * t109, -0.2e1 * t124, t27 ^ 2 * t99 + t26 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, -t35, 0, -t22, -t24, 0, 0, 0, t38, 0, 0, t35, 0, -t22, -pkin(3) * t38 - qJ(4) * t35, t24, -pkin(3) * t22 + qJ(4) * t24, 0, 0, -t20, 0, t18, 0, t7, t9, -t18 * t45 + t20 * t43, t43 * t7 + t45 * t9, -t15, t10, -t106, t15, -t103, 0, t69 * t85 + t6, t73 * t85 - t112, -t1, t1 * t42 + t7 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t110, -t60, 0, 0, 0, 0, 0, 1, 0, 0, t78 + t110, 0, t77 + t60, pkin(3) * t54 + qJ(4) * t51, 0, 0, 0, 0, 0, 1, -t47 - t58 + (qJ(4) + t51) * t70, t45 + t30, 0, t28 * t43 + t30 * t45, t63, t49, 0, t66, 0, 0, t100 * t73, -t100 * t69, -t123 - t124, t123 * t27 + t26 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t78, 0, t77, pkin(3) ^ 2 + qJ(4) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t43, 0.2e1 * t45, 0, t43 ^ 2 + t45 ^ 2, t63, t49, 0, t66, 0, 0, t41 * t116, -0.2e1 * t107, -0.2e1 * t123, t42 ^ 2 * t99 + t41 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t22, 0, 0, 0, 0, 0, 0, 0, 0, t84, t70 * t9 - t111, 0, 0, 0, 0, 0, 0, t84 * t69, t84 * t73, 0, t1 * t70 - t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t54, 0, 0, 0, 0, 0, 0, -t74, t70, 0, -t28 * t74 + t30 * t70, 0, 0, 0, 0, 0, 0, -t101, t52, -t37, -t26 * t74 + t27 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t74, t70, 0, -t43 * t74 + t45 * t70, 0, 0, 0, 0, 0, 0, -t101, t52, -t37, t37 * t42 - t41 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64 + t67, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64 * t99 + t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, -t18, 0, -t7, -t9, 0, 0, t15, -t10, t106, -t15, t103, 0, t69 * t87 - t6, t73 * t87 + t112, t1, -t7 * pkin(5) + pkin(10) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t28, -t30, 0, 0, -t63, t50, 0, -t66, 0, 0 (-pkin(5) - t26) * t73, t62 + t109, t124 - t93, -t26 * pkin(5) + pkin(10) * t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t43, -t45, 0, 0, -t63, t50, 0, -t66, 0, 0 (-pkin(5) - t41) * t73, t62 + t107, t123 - t93, -t41 * pkin(5) + pkin(10) * t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t70, 0, 0, 0, 0, 0, 0, 0, 0, t101, -t52, t37, t74 * pkin(5) + pkin(10) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t63, t49, 0, t66, 0, 0, pkin(5) * t116, -0.2e1 * t62, 0.2e1 * t93, pkin(10) ^ 2 * t99 + pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, 0, -t105, t18, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, 0, -t73, 0, -t69 * t27, -t73 * t27, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69, 0, -t73, 0, -t69 * t42, -t73 * t42, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69 * t70, -t73 * t70, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, t73, 0, -t69 * pkin(10), -t73 * pkin(10), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t2;
