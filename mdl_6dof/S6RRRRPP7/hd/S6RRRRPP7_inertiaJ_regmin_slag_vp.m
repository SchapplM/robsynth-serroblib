% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 21:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPP7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPP7_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t79 = sin(pkin(6));
t87 = cos(qJ(2));
t111 = t79 * t87;
t84 = sin(qJ(2));
t112 = t79 * t84;
t81 = cos(pkin(6));
t83 = sin(qJ(3));
t86 = cos(qJ(3));
t53 = t86 * t112 + t81 * t83;
t82 = sin(qJ(4));
t85 = cos(qJ(4));
t36 = t85 * t111 + t53 * t82;
t125 = -0.2e1 * t36;
t124 = -0.2e1 * t53;
t123 = -0.2e1 * t83;
t122 = 0.2e1 * t86;
t120 = pkin(1) * t87;
t64 = pkin(8) * t112;
t46 = t64 + (-pkin(2) - t120) * t81;
t52 = t83 * t112 - t81 * t86;
t22 = t52 * pkin(3) - t53 * pkin(10) + t46;
t121 = pkin(1) * t84;
t95 = pkin(8) * t111;
t47 = t95 + (pkin(9) + t121) * t81;
t48 = (-pkin(2) * t87 - pkin(9) * t84 - pkin(1)) * t79;
t27 = t86 * t47 + t83 * t48;
t24 = -pkin(10) * t111 + t27;
t12 = t82 * t22 + t85 * t24;
t10 = -t36 * qJ(5) + t12;
t78 = sin(pkin(11));
t11 = t85 * t22 - t82 * t24;
t37 = -t82 * t111 + t53 * t85;
t8 = t52 * pkin(4) - t37 * qJ(5) + t11;
t80 = cos(pkin(11));
t4 = t80 * t10 + t78 * t8;
t119 = pkin(3) * t85;
t118 = pkin(9) * t82;
t70 = t80 * pkin(4) + pkin(5);
t117 = pkin(5) + t70;
t26 = -t83 * t47 + t86 * t48;
t23 = pkin(3) * t111 - t26;
t116 = t23 * t82;
t115 = t23 * t85;
t114 = t37 * t82;
t74 = t79 ^ 2;
t113 = t74 * t87;
t110 = t81 * t84;
t109 = t82 * t52;
t108 = t82 * t83;
t107 = t82 * t85;
t106 = t82 * t86;
t105 = t83 * t52;
t104 = t85 * t52;
t103 = t85 * t83;
t102 = t85 * t86;
t101 = -qJ(5) - pkin(10);
t100 = qJ(5) * t83;
t61 = -t86 * pkin(3) - t83 * pkin(10) - pkin(2);
t58 = t85 * t61;
t34 = -t85 * t100 + t58 + (-pkin(4) - t118) * t86;
t96 = pkin(9) * t102;
t38 = t96 + (t61 - t100) * t82;
t18 = t78 * t34 + t80 * t38;
t73 = t83 * pkin(9);
t60 = pkin(4) * t108 + t73;
t99 = 0.2e1 * t111;
t98 = t83 * t122;
t1 = t52 * qJ(6) + t4;
t62 = t101 * t85;
t92 = t101 * t82;
t40 = -t62 * t78 - t80 * t92;
t42 = -t80 * t62 + t78 * t92;
t97 = t40 ^ 2 + t42 ^ 2;
t94 = t83 * t111;
t93 = t86 * t111;
t72 = -t85 * pkin(4) - pkin(3);
t3 = -t10 * t78 + t80 * t8;
t19 = t80 * t36 + t37 * t78;
t20 = -t78 * t36 + t80 * t37;
t91 = -t42 * t19 + t20 * t40;
t57 = t78 * t85 + t80 * t82;
t49 = t57 * t83;
t50 = t80 * t103 - t78 * t108;
t90 = t40 * t50 - t42 * t49;
t17 = t80 * t34 - t78 * t38;
t56 = t78 * t82 - t80 * t85;
t89 = 0.2e1 * t40 * t57 - 0.2e1 * t42 * t56;
t14 = t36 * pkin(4) + t23;
t77 = t85 ^ 2;
t76 = t83 ^ 2;
t75 = t82 ^ 2;
t67 = pkin(4) * t78 + qJ(6);
t55 = pkin(1) * t110 + t95;
t54 = t81 * t120 - t64;
t45 = t82 * t61 + t96;
t44 = -pkin(9) * t106 + t58;
t30 = t56 * pkin(5) - t57 * qJ(6) + t72;
t25 = pkin(5) * t49 - qJ(6) * t50 + t60;
t16 = t86 * pkin(5) - t17;
t15 = -t86 * qJ(6) + t18;
t5 = t19 * pkin(5) - t20 * qJ(6) + t14;
t2 = -t52 * pkin(5) - t3;
t6 = [1, 0, 0, t74 * t84 ^ 2, 0.2e1 * t84 * t113, 0.2e1 * t79 * t110, t81 * t99, t81 ^ 2, 0.2e1 * pkin(1) * t113 + 0.2e1 * t54 * t81, -0.2e1 * t74 * t121 - 0.2e1 * t55 * t81, t53 ^ 2, t52 * t124, t111 * t124, t52 * t99, t74 * t87 ^ 2, -0.2e1 * t26 * t111 + 0.2e1 * t46 * t52, 0.2e1 * t27 * t111 + 0.2e1 * t46 * t53, t37 ^ 2, t37 * t125, 0.2e1 * t37 * t52, t52 * t125, t52 ^ 2, 0.2e1 * t11 * t52 + 0.2e1 * t23 * t36, -0.2e1 * t12 * t52 + 0.2e1 * t23 * t37, -0.2e1 * t19 * t4 - 0.2e1 * t20 * t3, t14 ^ 2 + t3 ^ 2 + t4 ^ 2, 0.2e1 * t19 * t5 - 0.2e1 * t2 * t52, -0.2e1 * t1 * t19 + 0.2e1 * t2 * t20, 0.2e1 * t1 * t52 - 0.2e1 * t20 * t5, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, t112, t111, t81, t54, -t55, t53 * t83, t53 * t86 - t105, -t94, -t93, 0, -pkin(2) * t52 + pkin(9) * t94 - t46 * t86, -pkin(2) * t53 + pkin(9) * t93 + t46 * t83, t37 * t103 (-t36 * t85 - t114) * t83, t52 * t103 - t37 * t86, -t82 * t105 + t36 * t86, -t52 * t86, -t11 * t86 + t44 * t52 + (pkin(9) * t36 + t116) * t83, t12 * t86 - t45 * t52 + (pkin(9) * t37 + t115) * t83, -t17 * t20 - t18 * t19 - t3 * t50 - t4 * t49, t14 * t60 + t17 * t3 + t18 * t4, -t16 * t52 + t25 * t19 + t2 * t86 + t5 * t49, -t1 * t49 - t15 * t19 + t16 * t20 + t2 * t50, -t1 * t86 + t15 * t52 - t25 * t20 - t5 * t50, t1 * t15 + t16 * t2 + t25 * t5; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t76, t98, 0, 0, 0, pkin(2) * t122, pkin(2) * t123, t77 * t76, -0.2e1 * t76 * t107, t102 * t123, t82 * t98, t86 ^ 2, 0.2e1 * t76 * t118 - 0.2e1 * t44 * t86, 0.2e1 * t76 * pkin(9) * t85 + 0.2e1 * t45 * t86, -0.2e1 * t17 * t50 - 0.2e1 * t18 * t49, t17 ^ 2 + t18 ^ 2 + t60 ^ 2, 0.2e1 * t16 * t86 + 0.2e1 * t25 * t49, -0.2e1 * t15 * t49 + 0.2e1 * t16 * t50, -0.2e1 * t15 * t86 - 0.2e1 * t25 * t50, t15 ^ 2 + t16 ^ 2 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t52, -t111, t26, -t27, t114, -t82 * t36 + t37 * t85, t109, t104, 0, -pkin(3) * t36 - pkin(10) * t109 - t115, -pkin(3) * t37 - pkin(10) * t104 + t116, -t3 * t57 - t4 * t56 + t91, t14 * t72 - t3 * t40 + t4 * t42, t19 * t30 - t40 * t52 + t5 * t56, -t1 * t56 + t2 * t57 + t91, -t20 * t30 + t42 * t52 - t5 * t57, t1 * t42 + t2 * t40 + t30 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t86, 0, -t73, -t86 * pkin(9), t82 * t103 (-t75 + t77) * t83, -t106, -t102, 0, -pkin(9) * t103 + (-pkin(3) * t83 + pkin(10) * t86) * t82, pkin(10) * t102 + (t118 - t119) * t83, -t17 * t57 - t18 * t56 + t90, -t17 * t40 + t18 * t42 + t60 * t72, t25 * t56 + t30 * t49 + t40 * t86, -t15 * t56 + t16 * t57 + t90, -t25 * t57 - t30 * t50 - t42 * t86, t15 * t42 + t16 * t40 + t25 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t75, 0.2e1 * t107, 0, 0, 0, 0.2e1 * t119, -0.2e1 * pkin(3) * t82, t89, t72 ^ 2 + t97, 0.2e1 * t30 * t56, t89, -0.2e1 * t30 * t57, t30 ^ 2 + t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t36, t52, t11, -t12 (-t19 * t78 - t20 * t80) * pkin(4) (t3 * t80 + t4 * t78) * pkin(4), t117 * t52 + t3, -t19 * t67 - t20 * t70, t52 * t67 + t1, t1 * t67 - t2 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103, -t108, -t86, t44, -t45 (-t49 * t78 - t50 * t80) * pkin(4) (t17 * t80 + t18 * t78) * pkin(4), -t117 * t86 + t17, -t49 * t67 - t50 * t70 (-qJ(6) - t67) * t86 + t18, t15 * t67 - t16 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82, t85, 0, -t82 * pkin(10), -t85 * pkin(10) (-t56 * t78 - t57 * t80) * pkin(4) (-t40 * t80 + t42 * t78) * pkin(4), -t40, -t56 * t67 - t57 * t70, t42, -t40 * t70 + t42 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t78 ^ 2 + t80 ^ 2) * pkin(4) ^ 2, 0.2e1 * t70, 0, 0.2e1 * t67, t67 ^ 2 + t70 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t19, 0, -t20, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, t49, 0, -t50, t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t56, 0, -t57, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, t20, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t50, 0, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
