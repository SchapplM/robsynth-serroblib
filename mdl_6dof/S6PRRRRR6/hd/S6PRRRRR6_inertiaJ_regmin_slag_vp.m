% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d2,d3,d4,d5,d6,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRRR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRR6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PRRRRR6_inertiaJ_regmin_slag_vp: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t75 = sin(pkin(7));
t88 = cos(qJ(3));
t116 = t75 * t88;
t83 = sin(qJ(3));
t131 = pkin(2) * t83;
t78 = cos(pkin(7));
t61 = pkin(10) * t116 + t78 * t131;
t74 = sin(pkin(8));
t77 = cos(pkin(8));
t112 = t77 * t88;
t93 = t75 * t112;
t36 = (t74 * t78 + t93) * pkin(11) + t61;
t82 = sin(qJ(4));
t87 = cos(qJ(4));
t117 = t75 * t83;
t66 = t78 * t88 * pkin(2);
t41 = t78 * pkin(3) + t66 + (-pkin(11) * t77 - pkin(10)) * t117;
t46 = (-pkin(11) * t74 * t83 - pkin(3) * t88 - pkin(2)) * t75;
t90 = t41 * t77 + t46 * t74;
t17 = -t82 * t36 + t90 * t87;
t119 = t74 * t82;
t40 = t78 * t119 + (t82 * t112 + t83 * t87) * t75;
t53 = t74 * t116 - t77 * t78;
t81 = sin(qJ(5));
t86 = cos(qJ(5));
t26 = t81 * t40 + t86 * t53;
t138 = -0.2e1 * t26;
t118 = t74 * t87;
t37 = t82 * t117 - t78 * t118 - t87 * t93;
t137 = 0.2e1 * t37;
t136 = -0.2e1 * t40;
t57 = t86 * t119 + t81 * t77;
t80 = sin(qJ(6));
t85 = cos(qJ(6));
t42 = t85 * t118 + t80 * t57;
t135 = -0.2e1 * t42;
t134 = -0.2e1 * t57;
t133 = -0.2e1 * t81;
t132 = 0.2e1 * t86;
t130 = pkin(3) * t82;
t129 = pkin(5) * t85;
t128 = pkin(12) * t80;
t23 = -t74 * t41 + t77 * t46;
t12 = t37 * pkin(4) - t40 * pkin(12) + t23;
t18 = t87 * t36 + t90 * t82;
t14 = -t53 * pkin(12) + t18;
t7 = t86 * t12 - t81 * t14;
t3 = -t37 * pkin(5) - t7;
t127 = t3 * t80;
t126 = t3 * t85;
t27 = t86 * t40 - t81 * t53;
t22 = t85 * t27 + t80 * t37;
t125 = t22 * t80;
t94 = pkin(11) * t118;
t51 = t94 + (pkin(12) + t130) * t77;
t52 = (-pkin(4) * t87 - pkin(12) * t82 - pkin(3)) * t74;
t31 = -t81 * t51 + t86 * t52;
t29 = pkin(5) * t118 - t31;
t124 = t29 * t80;
t123 = t29 * t85;
t43 = -t80 * t118 + t85 * t57;
t122 = t43 * t80;
t69 = t74 ^ 2;
t121 = t69 * t87;
t70 = t75 ^ 2;
t120 = t70 * t88;
t76 = sin(pkin(6));
t89 = cos(qJ(2));
t115 = t76 * t89;
t114 = t77 * t82;
t113 = t77 * t87;
t111 = t78 * t89;
t110 = t80 * t26;
t56 = t81 * t119 - t86 * t77;
t109 = t80 * t56;
t108 = t80 * t81;
t107 = t80 * t85;
t106 = t80 * t86;
t105 = t81 * t26;
t104 = t81 * t37;
t103 = t81 * t56;
t102 = t85 * t26;
t101 = t85 * t56;
t100 = t85 * t81;
t99 = t85 * t86;
t98 = t86 * t37;
t97 = 0.2e1 * t118;
t96 = 0.2e1 * t75 * t78;
t95 = t81 * t132;
t92 = t81 * t118;
t91 = t86 * t118;
t8 = t81 * t12 + t86 * t14;
t32 = t86 * t51 + t81 * t52;
t65 = pkin(11) * t119;
t50 = t65 + (-pkin(3) * t87 - pkin(4)) * t77;
t13 = t53 * pkin(4) - t17;
t84 = sin(qJ(2));
t79 = cos(pkin(6));
t73 = t85 ^ 2;
t72 = t81 ^ 2;
t71 = t80 ^ 2;
t63 = -t86 * pkin(5) - t81 * pkin(13) - pkin(4);
t60 = pkin(3) * t114 + t94;
t59 = -pkin(10) * t117 + t66;
t58 = pkin(3) * t113 - t65;
t55 = -t75 * t115 + t79 * t78;
t48 = pkin(12) * t99 + t80 * t63;
t47 = -pkin(12) * t106 + t85 * t63;
t39 = t79 * t117 + (t83 * t111 + t84 * t88) * t76;
t38 = t79 * t116 + (t88 * t111 - t83 * t84) * t76;
t30 = -pkin(13) * t118 + t32;
t28 = t56 * pkin(5) - t57 * pkin(13) + t50;
t25 = -t38 * t74 + t55 * t77;
t21 = t80 * t27 - t85 * t37;
t20 = t39 * t87 + (t38 * t77 + t55 * t74) * t82;
t19 = -t38 * t113 - t55 * t118 + t39 * t82;
t16 = t80 * t28 + t85 * t30;
t15 = t85 * t28 - t80 * t30;
t11 = t20 * t86 + t25 * t81;
t10 = t20 * t81 - t25 * t86;
t9 = t26 * pkin(5) - t27 * pkin(13) + t13;
t6 = t11 * t85 + t19 * t80;
t5 = -t11 * t80 + t19 * t85;
t4 = t37 * pkin(13) + t8;
t2 = t85 * t4 + t80 * t9;
t1 = -t80 * t4 + t85 * t9;
t24 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t115, -t76 * t84, 0, 0, 0, 0, 0, -t55 * t116 + t38 * t78, t55 * t117 - t39 * t78, 0, 0, 0, 0, 0, t19 * t53 + t25 * t37, t20 * t53 + t25 * t40, 0, 0, 0, 0, 0, -t10 * t37 + t19 * t26, -t11 * t37 + t19 * t27, 0, 0, 0, 0, 0, t10 * t21 + t5 * t26, t10 * t22 - t6 * t26; 0, 1, 0, 0, t70 * t83 ^ 2, 0.2e1 * t83 * t120, t83 * t96, t88 * t96, t78 ^ 2, 0.2e1 * pkin(2) * t120 + 0.2e1 * t59 * t78, -0.2e1 * t70 * t131 - 0.2e1 * t61 * t78, t40 ^ 2, t37 * t136, t53 * t136, t53 * t137, t53 ^ 2, -0.2e1 * t17 * t53 + 0.2e1 * t23 * t37, 0.2e1 * t18 * t53 + 0.2e1 * t23 * t40, t27 ^ 2, t27 * t138, t27 * t137, t37 * t138, t37 ^ 2, 0.2e1 * t13 * t26 + 0.2e1 * t7 * t37, 0.2e1 * t13 * t27 - 0.2e1 * t8 * t37, t22 ^ 2, -0.2e1 * t22 * t21, 0.2e1 * t22 * t26, t21 * t138, t26 ^ 2, 0.2e1 * t1 * t26 + 0.2e1 * t3 * t21, -0.2e1 * t2 * t26 + 0.2e1 * t3 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t39, 0, 0, 0, 0, 0, -t25 * t118 - t19 * t77, t25 * t119 - t20 * t77, 0, 0, 0, 0, 0, t10 * t118 + t19 * t56, t11 * t118 + t19 * t57, 0, 0, 0, 0, 0, t10 * t42 + t5 * t56, t10 * t43 - t6 * t56; 0, 0, 0, 0, 0, 0, t117, t116, t78, t59, -t61, t40 * t119 (-t37 * t82 + t40 * t87) * t74, -t53 * t119 + t77 * t40, -t53 * t118 - t77 * t37, -t53 * t77, t17 * t77 - t58 * t53 + (-pkin(3) * t37 - t23 * t87) * t74, -t18 * t77 + t60 * t53 + (-pkin(3) * t40 + t23 * t82) * t74, t27 * t57, -t57 * t26 - t27 * t56, -t27 * t118 + t57 * t37, t26 * t118 - t56 * t37, -t37 * t118, -t7 * t118 + t13 * t56 + t50 * t26 + t31 * t37, t8 * t118 + t13 * t57 + t50 * t27 - t32 * t37, t22 * t43, -t43 * t21 - t22 * t42, t22 * t56 + t43 * t26, -t21 * t56 - t42 * t26, t26 * t56, t1 * t56 + t15 * t26 + t29 * t21 + t3 * t42, -t16 * t26 - t2 * t56 + t29 * t22 + t3 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t69 * t82 ^ 2, 0.2e1 * t82 * t121, 0.2e1 * t74 * t114, t77 * t97, t77 ^ 2, 0.2e1 * pkin(3) * t121 + 0.2e1 * t58 * t77, -0.2e1 * t69 * t130 - 0.2e1 * t60 * t77, t57 ^ 2, t56 * t134, t118 * t134, t56 * t97, t69 * t87 ^ 2, -0.2e1 * t31 * t118 + 0.2e1 * t50 * t56, 0.2e1 * t32 * t118 + 0.2e1 * t50 * t57, t43 ^ 2, t43 * t135, 0.2e1 * t43 * t56, t56 * t135, t56 ^ 2, 0.2e1 * t15 * t56 + 0.2e1 * t29 * t42, -0.2e1 * t16 * t56 + 0.2e1 * t29 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t20, 0, 0, 0, 0, 0, -t19 * t86, t19 * t81, 0, 0, 0, 0, 0, t10 * t108 - t5 * t86, t10 * t100 + t6 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t37, -t53, t17, -t18, t27 * t81, t27 * t86 - t105, t104, t98, 0, -pkin(4) * t26 - pkin(12) * t104 - t13 * t86, -pkin(4) * t27 - pkin(12) * t98 + t13 * t81, t22 * t100 (-t21 * t85 - t125) * t81, t100 * t26 - t22 * t86, -t105 * t80 + t21 * t86, -t26 * t86, -t1 * t86 + t47 * t26 + (pkin(12) * t21 + t127) * t81, t2 * t86 - t48 * t26 + (pkin(12) * t22 + t126) * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t119, t118, t77, t58, -t60, t57 * t81, t57 * t86 - t103, -t92, -t91, 0, -pkin(4) * t56 + pkin(12) * t92 - t50 * t86, -pkin(4) * t57 + pkin(12) * t91 + t50 * t81, t43 * t100 (-t42 * t85 - t122) * t81, t100 * t56 - t43 * t86, -t103 * t80 + t42 * t86, -t56 * t86, -t15 * t86 + t47 * t56 + (pkin(12) * t42 + t124) * t81, t16 * t86 - t48 * t56 + (pkin(12) * t43 + t123) * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t72, t95, 0, 0, 0, pkin(4) * t132, pkin(4) * t133, t73 * t72, -0.2e1 * t72 * t107, t99 * t133, t80 * t95, t86 ^ 2, 0.2e1 * t128 * t72 - 0.2e1 * t47 * t86, 0.2e1 * t72 * pkin(12) * t85 + 0.2e1 * t48 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t10 * t85, t10 * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, t37, t7, -t8, t125, -t80 * t21 + t22 * t85, t110, t102, 0, -pkin(5) * t21 - pkin(13) * t110 - t126, -pkin(5) * t22 - pkin(13) * t102 + t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t56, -t118, t31, -t32, t122, -t80 * t42 + t43 * t85, t109, t101, 0, -pkin(5) * t42 - pkin(13) * t109 - t123, -pkin(5) * t43 - pkin(13) * t101 + t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t86, 0, -t81 * pkin(12), -t86 * pkin(12), t80 * t100 (-t71 + t73) * t81, -t106, -t99, 0, -pkin(12) * t100 + (-pkin(5) * t81 + pkin(13) * t86) * t80, pkin(13) * t99 + (t128 - t129) * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t71, 0.2e1 * t107, 0, 0, 0, 0.2e1 * t129, -0.2e1 * pkin(5) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, t26, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, -t42, t56, t15, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, -t108, -t86, t47, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t85, 0, -t80 * pkin(13), -t85 * pkin(13); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t24;
