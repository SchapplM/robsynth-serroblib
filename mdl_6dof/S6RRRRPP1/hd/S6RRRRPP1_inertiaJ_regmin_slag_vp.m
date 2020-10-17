% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP1_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 17:56:55
% EndTime: 2019-05-07 17:56:58
% DurationCPUTime: 1.00s
% Computational Cost: add. (1769->136), mult. (3246->232), div. (0->0), fcn. (3728->8), ass. (0->93)
t85 = cos(qJ(4));
t117 = t85 * pkin(9);
t77 = t85 * qJ(5);
t65 = t77 + t117;
t80 = sin(pkin(10));
t81 = cos(pkin(10));
t82 = sin(qJ(4));
t94 = (-qJ(5) - pkin(9)) * t82;
t44 = t80 * t65 - t81 * t94;
t46 = t81 * t65 + t80 * t94;
t57 = t80 * t82 - t81 * t85;
t58 = t80 * t85 + t81 * t82;
t125 = t44 * t58 - t46 * t57;
t129 = 0.2e1 * t125;
t83 = sin(qJ(3));
t119 = t83 * pkin(2);
t95 = pkin(9) + t119;
t90 = t85 * t95;
t56 = t90 + t77;
t88 = (-qJ(5) - t95) * t82;
t35 = t80 * t56 - t81 * t88;
t37 = t81 * t56 + t80 * t88;
t126 = t35 * t58 - t37 * t57;
t128 = 0.2e1 * t126;
t127 = t125 + t126;
t124 = 0.2e1 * t57;
t123 = -0.2e1 * t58;
t115 = cos(qJ(2));
t76 = -t115 * pkin(2) - pkin(1);
t122 = 0.2e1 * t76;
t114 = cos(qJ(3));
t84 = sin(qJ(2));
t61 = -t114 * t115 + t83 * t84;
t62 = t114 * t84 + t83 * t115;
t40 = t61 * pkin(3) - t62 * pkin(9) + t76;
t66 = (-pkin(8) - pkin(7)) * t84;
t97 = t115 * pkin(7);
t67 = t115 * pkin(8) + t97;
t48 = t114 * t67 + t83 * t66;
t18 = t85 * t40 - t82 * t48;
t14 = t61 * pkin(4) - t62 * t77 + t18;
t108 = t85 * t48;
t16 = t108 + (-qJ(5) * t62 + t40) * t82;
t8 = t80 * t14 + t81 * t16;
t3 = t61 * qJ(6) + t8;
t7 = t81 * t14 - t80 * t16;
t4 = -t61 * pkin(5) - t7;
t121 = -t3 * t57 + t4 * t58;
t120 = -t8 * t57 - t7 * t58;
t118 = t85 * pkin(4);
t96 = t114 * pkin(2);
t74 = -t96 - pkin(3);
t116 = pkin(3) - t74;
t47 = -t114 * t66 + t83 * t67;
t111 = t47 * t85;
t110 = t82 * t62;
t109 = t82 * t85;
t107 = t85 * t62;
t75 = -pkin(3) - t118;
t39 = t57 * pkin(5) - t58 * qJ(6) + t75;
t30 = -t96 + t39;
t102 = t30 + t39;
t101 = 0.2e1 * t115;
t100 = -0.2e1 * t62 * t61;
t99 = t35 ^ 2 + t37 ^ 2;
t98 = t44 ^ 2 + t46 ^ 2;
t31 = t58 * t62;
t32 = t81 * t107 - t80 * t110;
t93 = -t37 * t31 + t35 * t32;
t92 = -t46 * t31 + t44 * t32;
t91 = t35 * t44 + t37 * t46;
t25 = pkin(4) * t110 + t47;
t89 = -pkin(3) * t62 - pkin(9) * t61;
t87 = -t95 * t61 + t74 * t62;
t79 = t85 ^ 2;
t78 = t82 ^ 2;
t71 = t81 * pkin(4) + pkin(5);
t69 = t80 * pkin(4) + qJ(6);
t68 = 0.2e1 * t109;
t64 = t74 - t118;
t59 = t62 ^ 2;
t54 = t85 * t61;
t53 = t82 * t61;
t50 = t82 * t107;
t42 = t47 * t82;
t41 = (-t78 + t79) * t62;
t38 = (-t57 * t80 - t58 * t81) * pkin(4);
t26 = -t69 * t57 - t71 * t58;
t19 = t82 * t40 + t108;
t11 = t31 * pkin(5) - t32 * qJ(6) + t25;
t10 = t11 * t58;
t9 = t11 * t57;
t1 = [1, 0, 0, t84 ^ 2, t84 * t101, 0, 0, 0, pkin(1) * t101, -0.2e1 * pkin(1) * t84, t59, t100, 0, 0, 0, t61 * t122, t62 * t122, t79 * t59, -0.2e1 * t59 * t109, 0.2e1 * t61 * t107, t82 * t100, t61 ^ 2, 0.2e1 * t47 * t110 + 0.2e1 * t18 * t61, 0.2e1 * t47 * t107 - 0.2e1 * t19 * t61, -0.2e1 * t8 * t31 - 0.2e1 * t7 * t32, t25 ^ 2 + t7 ^ 2 + t8 ^ 2, 0.2e1 * t11 * t31 - 0.2e1 * t4 * t61, -0.2e1 * t3 * t31 + 0.2e1 * t4 * t32, -0.2e1 * t11 * t32 + 0.2e1 * t3 * t61, t11 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, t84, t115, 0, -t84 * pkin(7), -t97, 0, 0, t62, -t61, 0, -t47, -t48, t50, t41, t53, t54, 0, t87 * t82 - t111, t85 * t87 + t42, t93 + t120, t25 * t64 - t7 * t35 + t8 * t37, t30 * t31 - t35 * t61 + t9, t93 + t121, -t30 * t32 + t37 * t61 - t10, t11 * t30 + t3 * t37 + t4 * t35; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t96, -0.2e1 * t119, t78, t68, 0, 0, 0, -0.2e1 * t74 * t85, 0.2e1 * t74 * t82, t128, t64 ^ 2 + t99, t30 * t124, t128, t30 * t123, t30 ^ 2 + t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t61, 0, -t47, -t48, t50, t41, t53, t54, 0, t89 * t82 - t111, t89 * t85 + t42, t92 + t120, t25 * t75 - t7 * t44 + t8 * t46, t39 * t31 - t44 * t61 + t9, t92 + t121, -t39 * t32 + t46 * t61 - t10, t11 * t39 + t3 * t46 + t4 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t96, -t119, t78, t68, 0, 0, 0, t116 * t85, -t116 * t82, t127, t64 * t75 + t91, t102 * t57, t127, -t102 * t58, t30 * t39 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t78, t68, 0, 0, 0, 0.2e1 * pkin(3) * t85, -0.2e1 * pkin(3) * t82, t129, t75 ^ 2 + t98, t39 * t124, t129, t39 * t123, t39 ^ 2 + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t107, -t110, t61, t18, -t19 (-t31 * t80 - t32 * t81) * pkin(4) (t7 * t81 + t8 * t80) * pkin(4) (pkin(5) + t71) * t61 + t7, -t69 * t31 - t71 * t32, t69 * t61 + t3, t3 * t69 - t4 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t85, 0, -t82 * t95, -t90, t38 (-t35 * t81 + t37 * t80) * pkin(4), -t35, t26, t37, -t35 * t71 + t37 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t85, 0, -t82 * pkin(9), -t117, t38 (-t44 * t81 + t46 * t80) * pkin(4), -t44, t26, t46, -t44 * t71 + t46 * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t80 ^ 2 + t81 ^ 2) * pkin(4) ^ 2, 0.2e1 * t71, 0, 0.2e1 * t69, t69 ^ 2 + t71 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t31, 0, -t32, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t57, 0, -t58, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t57, 0, -t58, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t32, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t1;
