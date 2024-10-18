% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha5,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x31]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 22:28
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR15_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR15_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR15_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 22:26:25
% EndTime: 2024-09-27 22:26:26
% DurationCPUTime: 1.00s
% Computational Cost: add. (1663->122), mult. (4495->220), div. (0->0), fcn. (5137->12), ass. (0->101)
t76 = sin(pkin(5));
t86 = cos(qJ(2));
t101 = t76 * t86;
t82 = sin(qJ(2));
t102 = t76 * t82;
t81 = sin(qJ(3));
t85 = cos(qJ(3));
t42 = -t85 * t101 + t81 * t102;
t43 = (t81 * t86 + t82 * t85) * t76;
t80 = sin(qJ(4));
t84 = cos(qJ(4));
t30 = -t80 * t42 + t84 * t43;
t79 = sin(qJ(5));
t83 = cos(qJ(5));
t29 = t84 * t42 + t80 * t43;
t75 = sin(pkin(6));
t77 = cos(pkin(6));
t78 = cos(pkin(5));
t87 = -t29 * t77 + t75 * t78;
t16 = t83 * t30 + t79 * t87;
t116 = -0.2e1 * t16;
t115 = -0.2e1 * t78;
t114 = 0.2e1 * t78;
t113 = pkin(8) + pkin(9);
t112 = pkin(1) * t82;
t51 = (-pkin(2) * t86 - pkin(1)) * t76;
t31 = t42 * pkin(3) + t51;
t68 = t75 * pkin(11);
t13 = t29 * pkin(4) - t30 * t68 + t31;
t109 = t78 * pkin(3);
t110 = t78 * pkin(2);
t61 = t78 * t86 * pkin(1);
t34 = -t113 * t102 + t110 + t61;
t91 = t78 * t112;
t38 = t113 * t101 + t91;
t21 = t85 * t34 - t81 * t38;
t18 = -t43 * pkin(10) + t109 + t21;
t97 = t85 * t38;
t22 = t81 * t34 + t97;
t19 = -t42 * pkin(10) + t22;
t11 = t84 * t18 - t80 * t19;
t8 = -t30 * t77 * pkin(11) + t78 * pkin(4) + t11;
t5 = t77 * t13 - t75 * t8;
t111 = t5 * t83;
t108 = t80 * pkin(3);
t107 = t81 * pkin(2);
t69 = t84 * pkin(3);
t70 = t85 * pkin(2);
t106 = t16 * t75;
t71 = t75 ^ 2;
t105 = t71 * t79;
t104 = t71 * t83;
t72 = t76 ^ 2;
t103 = t72 * t86;
t64 = t75 * t79;
t65 = t75 * t83;
t100 = t77 * t79;
t99 = t77 * t83;
t98 = t84 * t19;
t67 = t70 + pkin(3);
t90 = t84 * t107;
t46 = t80 * t67 + t90;
t40 = t46 + t68;
t45 = -t80 * t107 + t84 * t67;
t44 = pkin(4) + t45;
t27 = -t79 * t40 + t44 * t99;
t96 = t44 * t104 + t27 * t77;
t56 = t68 + t108;
t66 = t69 + pkin(4);
t35 = -t79 * t56 + t66 * t99;
t95 = t66 * t104 + t35 * t77;
t47 = pkin(4) * t99 - pkin(11) * t64;
t94 = pkin(4) * t104 + t47 * t77;
t93 = 0.2e1 * t75 * t77;
t92 = t76 * t114;
t12 = t80 * t18 + t98;
t7 = pkin(11) * t87 + t12;
t88 = t13 * t75 + t77 * t8;
t3 = t83 * t7 + t79 * t88;
t89 = -t3 * t77 + t5 * t64;
t74 = t78 ^ 2;
t73 = t77 ^ 2;
t63 = t71 * t79 ^ 2;
t55 = 0.2e1 * t79 * t104;
t54 = t83 * t93;
t53 = t79 * t93;
t50 = pkin(8) * t101 + t91;
t49 = pkin(4) * t100 + pkin(11) * t65;
t48 = -pkin(8) * t102 + t61;
t36 = t66 * t100 + t83 * t56;
t28 = t44 * t100 + t83 * t40;
t23 = -t75 * t29 - t77 * t78;
t20 = t23 * t77;
t15 = t29 * t99 + t79 * t30 - t78 * t65;
t14 = t16 * t64;
t10 = -t15 * t77 - t23 * t65;
t9 = t16 * t77 - t23 * t64;
t6 = (-t15 * t79 + t16 * t83) * t75;
t2 = -t79 * t7 + t83 * t88;
t1 = t2 * t77;
t4 = [1, 0, 0, t72 * t82 ^ 2, 0.2e1 * t82 * t103, t82 * t92, t86 * t92, t74, 0.2e1 * pkin(1) * t103 + 0.2e1 * t48 * t78, -0.2e1 * t72 * t112 - 0.2e1 * t50 * t78, t43 ^ 2, -0.2e1 * t43 * t42, t43 * t114, t42 * t115, t74, 0.2e1 * t21 * t78 + 0.2e1 * t51 * t42, -0.2e1 * t22 * t78 + 0.2e1 * t51 * t43, t30 ^ 2, -0.2e1 * t30 * t29, t30 * t114, t29 * t115, t74, 0.2e1 * t11 * t78 + 0.2e1 * t31 * t29, -0.2e1 * t12 * t78 + 0.2e1 * t31 * t30, t16 ^ 2, t15 * t116, t23 * t116, 0.2e1 * t15 * t23, t23 ^ 2, 0.2e1 * t5 * t15 - 0.2e1 * t2 * t23, 0.2e1 * t5 * t16 + 0.2e1 * t3 * t23; 0, 0, 0, 0, 0, t102, t101, t78, t48, -t50, 0, 0, t43, -t42, t78, t78 * t70 + t21, -t97 + (-t34 - t110) * t81, 0, 0, t30, -t29, t78, t45 * t78 + t11, -t46 * t78 - t12, t14, t6, t9, t10, -t20, -t27 * t23 + t1 + (-t15 * t44 - t111) * t75, -t44 * t106 + t28 * t23 + t89; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t70, -0.2e1 * t107, 0, 0, 0, 0, 1, 0.2e1 * t45, -0.2e1 * t46, t63, t55, t53, t54, t73, 0.2e1 * t96, -0.2e1 * t44 * t105 - 0.2e1 * t28 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t42, t78, t21, -t22, 0, 0, t30, -t29, t78, t78 * t69 + t11, -t98 + (-t18 - t109) * t80, t14, t6, t9, t10, -t20, -t35 * t23 + t1 + (-t15 * t66 - t111) * t75, -t66 * t106 + t36 * t23 + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t70, -t107, 0, 0, 0, 0, 1, t45 + t69, -t90 + (-pkin(3) - t67) * t80, t63, t55, t53, t54, t73, t95 + t96, (-t28 - t36) * t77 + (-t44 - t66) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t69, -0.2e1 * t108, t63, t55, t53, t54, t73, 0.2e1 * t95, -0.2e1 * t66 * t105 - 0.2e1 * t36 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, t78, t11, -t12, t14, t6, t9, t10, -t20, -t47 * t23 + t1 + (-pkin(4) * t15 - t111) * t75, -pkin(4) * t106 + t49 * t23 + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t45, -t46, t63, t55, t53, t54, t73, t94 + t96, (-t28 - t49) * t77 + (-pkin(4) - t44) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t69, -t108, t63, t55, t53, t54, t73, t94 + t95, (-t36 - t49) * t77 + (-pkin(4) - t66) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t63, t55, t53, t54, t73, 0.2e1 * t94, -0.2e1 * pkin(4) * t105 - 0.2e1 * t49 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, -t23, t2, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t65, t77, t27, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t65, t77, t35, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t65, t77, t47, -t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t4;
