% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:13
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPRP7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRP7_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 08:22:09
% EndTime: 2019-05-07 08:22:13
% DurationCPUTime: 1.26s
% Computational Cost: add. (2440->171), mult. (5435->336), div. (0->0), fcn. (6281->10), ass. (0->102)
t71 = sin(pkin(6));
t79 = cos(qJ(2));
t108 = t71 * t79;
t76 = sin(qJ(2));
t109 = t71 * t76;
t73 = cos(pkin(6));
t75 = sin(qJ(3));
t78 = cos(qJ(3));
t49 = t78 * t109 + t73 * t75;
t70 = sin(pkin(11));
t72 = cos(pkin(11));
t58 = t75 * t109;
t89 = t73 * t78 - t58;
t32 = t72 * t49 + t70 * t89;
t74 = sin(qJ(5));
t77 = cos(qJ(5));
t25 = t77 * t108 + t32 * t74;
t121 = -0.2e1 * t25;
t120 = -0.2e1 * t77;
t119 = 0.2e1 * t78;
t118 = pkin(1) * t76;
t117 = pkin(1) * t79;
t31 = t49 * t70 - t72 * t89;
t116 = pkin(5) * t31;
t54 = t70 * t75 - t72 * t78;
t115 = pkin(5) * t54;
t60 = pkin(8) * t109;
t66 = -pkin(3) * t78 - pkin(2);
t35 = t58 * pkin(3) + t60 + (t66 - t117) * t73;
t13 = t31 * pkin(4) - t32 * pkin(10) + t35;
t97 = pkin(8) * t108;
t43 = t97 + (pkin(9) + t118) * t73;
t44 = (-pkin(2) * t79 - pkin(9) * t76 - pkin(1)) * t71;
t27 = -t43 * t75 + t78 * t44;
t17 = -pkin(3) * t108 - qJ(4) * t49 + t27;
t28 = t78 * t43 + t75 * t44;
t22 = t89 * qJ(4) + t28;
t10 = t70 * t17 + t72 * t22;
t8 = -pkin(10) * t108 + t10;
t4 = t74 * t13 + t77 * t8;
t114 = t25 * t77;
t26 = -t74 * t108 + t32 * t77;
t113 = t26 * t74;
t112 = t26 * t77;
t64 = pkin(3) * t70 + pkin(10);
t111 = t54 * t64;
t67 = t71 ^ 2;
t110 = t67 * t79;
t29 = t74 * t31;
t46 = t74 * t54;
t55 = t70 * t78 + t72 * t75;
t107 = t74 * t55;
t106 = t74 * t64;
t105 = t74 * t77;
t30 = t77 * t31;
t47 = t77 * t54;
t48 = t77 * t55;
t104 = t77 * t64;
t103 = -qJ(4) - pkin(9);
t36 = pkin(4) * t54 - pkin(10) * t55 + t66;
t57 = t103 * t78;
t91 = t103 * t75;
t40 = -t72 * t57 + t70 * t91;
t19 = t74 * t36 + t77 * t40;
t68 = t74 ^ 2;
t69 = t77 ^ 2;
t102 = t68 + t69;
t101 = qJ(6) * t31;
t100 = qJ(6) * t54;
t99 = 0.2e1 * t71 * t73;
t98 = -0.2e1 * t108;
t96 = t31 * t106;
t95 = t31 * t104;
t94 = t75 * t108;
t93 = t78 * t108;
t65 = -pkin(3) * t72 - pkin(4);
t92 = -t77 * t13 + t74 * t8;
t9 = t17 * t72 - t70 * t22;
t90 = -t77 * t36 + t40 * t74;
t38 = -t57 * t70 - t72 * t91;
t7 = pkin(4) * t108 - t9;
t1 = t101 + t4;
t2 = t92 - t116;
t88 = t1 * t77 + t2 * t74;
t87 = t1 * t74 - t2 * t77;
t86 = pkin(5) * t77 + qJ(6) * t74;
t85 = pkin(5) * t74 - qJ(6) * t77;
t14 = t100 + t19;
t15 = t90 - t115;
t84 = t14 * t77 + t15 * t74;
t83 = t14 * t74 - t15 * t77;
t52 = t65 - t86;
t82 = t52 * t55 - t111;
t81 = t55 * t65 - t111;
t53 = t55 ^ 2;
t51 = t73 * t118 + t97;
t50 = t73 * t117 - t60;
t42 = t60 + (-pkin(2) - t117) * t73;
t24 = t74 * t25;
t23 = t85 * t55 + t38;
t5 = pkin(5) * t25 - qJ(6) * t26 + t7;
t3 = [1, 0, 0, t67 * t76 ^ 2, 0.2e1 * t76 * t110, t76 * t99, t79 * t99, t73 ^ 2, 0.2e1 * pkin(1) * t110 + 0.2e1 * t50 * t73, -0.2e1 * t67 * t118 - 0.2e1 * t51 * t73, t49 ^ 2, 0.2e1 * t49 * t89, t49 * t98, t89 * t98, t67 * t79 ^ 2, -0.2e1 * t27 * t108 - 0.2e1 * t42 * t89, 0.2e1 * t28 * t108 + 0.2e1 * t42 * t49, -0.2e1 * t10 * t31 - 0.2e1 * t32 * t9, t10 ^ 2 + t35 ^ 2 + t9 ^ 2, t26 ^ 2, t26 * t121, 0.2e1 * t26 * t31, t31 * t121, t31 ^ 2, 0.2e1 * t25 * t7 - 0.2e1 * t31 * t92, 0.2e1 * t26 * t7 - 0.2e1 * t31 * t4, -0.2e1 * t2 * t31 + 0.2e1 * t25 * t5, -0.2e1 * t1 * t25 + 0.2e1 * t2 * t26, 0.2e1 * t1 * t31 - 0.2e1 * t26 * t5, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t109, t108, t73, t50, -t51, t49 * t75, t49 * t78 + t75 * t89, -t94, -t93, 0, pkin(2) * t89 + pkin(9) * t94 - t42 * t78, -pkin(2) * t49 + pkin(9) * t93 + t42 * t75, -t10 * t54 - t31 * t40 + t32 * t38 - t55 * t9, t10 * t40 + t35 * t66 - t38 * t9, t26 * t48 (-t113 - t114) * t55, t26 * t54 + t31 * t48, -t31 * t107 - t25 * t54, t31 * t54, t7 * t107 + t25 * t38 - t31 * t90 - t54 * t92, -t19 * t31 + t26 * t38 - t4 * t54 + t7 * t48, t5 * t107 - t15 * t31 - t2 * t54 + t23 * t25, -t14 * t25 + t15 * t26 - t55 * t87, t1 * t54 + t14 * t31 - t23 * t26 - t48 * t5, t1 * t14 + t15 * t2 + t23 * t5; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t75 ^ 2, t75 * t119, 0, 0, 0, pkin(2) * t119, -0.2e1 * pkin(2) * t75, 0.2e1 * t38 * t55 - 0.2e1 * t40 * t54, t38 ^ 2 + t40 ^ 2 + t66 ^ 2, t69 * t53, -0.2e1 * t53 * t105, 0.2e1 * t54 * t48, -0.2e1 * t54 * t107, t54 ^ 2, 0.2e1 * t38 * t107 - 0.2e1 * t54 * t90, -0.2e1 * t19 * t54 + 0.2e1 * t38 * t48, 0.2e1 * t23 * t107 - 0.2e1 * t15 * t54, -0.2e1 * t83 * t55, 0.2e1 * t14 * t54 - 0.2e1 * t23 * t48, t14 ^ 2 + t15 ^ 2 + t23 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t89, -t108, t27, -t28 (-t31 * t70 - t32 * t72) * pkin(3) (t10 * t70 + t72 * t9) * pkin(3), t113, -t24 + t112, t29, t30, 0, t25 * t65 - t7 * t77 - t96, t26 * t65 + t7 * t74 - t95, t25 * t52 - t5 * t77 - t96 (t113 - t114) * t64 + t88, -t26 * t52 - t5 * t74 + t95, t5 * t52 + t64 * t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t78, 0, -t75 * pkin(9), -t78 * pkin(9) (-t54 * t70 - t55 * t72) * pkin(3) (-t38 * t72 + t40 * t70) * pkin(3), t74 * t48 (-t68 + t69) * t55, t46, t47, 0, -t38 * t77 + t81 * t74, t38 * t74 + t81 * t77, -t23 * t77 + t82 * t74, t84, -t23 * t74 - t77 * t82, t23 * t52 + t64 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t70 ^ 2 + t72 ^ 2) * pkin(3) ^ 2, t68, 0.2e1 * t105, 0, 0, 0, t65 * t120, 0.2e1 * t65 * t74, t52 * t120, 0.2e1 * t102 * t64, -0.2e1 * t52 * t74, t102 * t64 ^ 2 + t52 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, t30, -t29, t30, -t24 - t112, t29, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, 0, 0, 0, 0, 0, t47, -t46, t47, -t102 * t55, t46, t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, -t25, t31, -t92, -t4, -t92 + 0.2e1 * t116, -pkin(5) * t26 - qJ(6) * t25, 0.2e1 * t101 + t4, -pkin(5) * t2 + qJ(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t107, t54, -t90, -t19, -t90 + 0.2e1 * t115, -t86 * t55, 0.2e1 * t100 + t19, -pkin(5) * t15 + qJ(6) * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t77, 0, -t106, -t104, -t106, -t85, t104, -t85 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, -t74, t77, 0, t74, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t26, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t54, t48, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0, t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
