% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d2,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x31]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR12_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR12_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RRRRR12_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t68 = sin(qJ(2));
t111 = pkin(1) * t68;
t64 = cos(pkin(5));
t62 = sin(pkin(5));
t72 = cos(qJ(2));
t95 = t62 * t72;
t48 = pkin(8) * t95 + t64 * t111;
t61 = sin(pkin(6));
t63 = cos(pkin(6));
t93 = t63 * t72;
t76 = t62 * t93;
t27 = (t61 * t64 + t76) * pkin(9) + t48;
t67 = sin(qJ(3));
t71 = cos(qJ(3));
t53 = t64 * t72 * pkin(1);
t96 = t62 * t68;
t30 = t64 * pkin(2) + t53 + (-pkin(9) * t63 - pkin(8)) * t96;
t34 = (-pkin(9) * t61 * t68 - pkin(2) * t72 - pkin(1)) * t62;
t73 = t30 * t63 + t34 * t61;
t13 = -t67 * t27 + t73 * t71;
t98 = t61 * t67;
t29 = t64 * t98 + (t67 * t93 + t68 * t71) * t62;
t41 = t61 * t95 - t64 * t63;
t66 = sin(qJ(4));
t70 = cos(qJ(4));
t18 = t29 * t66 + t41 * t70;
t118 = -0.2e1 * t18;
t97 = t61 * t71;
t28 = -t64 * t97 + t67 * t96 - t71 * t76;
t117 = 0.2e1 * t28;
t116 = -0.2e1 * t29;
t44 = t66 * t63 + t70 * t98;
t65 = sin(qJ(5));
t69 = cos(qJ(5));
t31 = t65 * t44 + t69 * t97;
t115 = -0.2e1 * t31;
t114 = -0.2e1 * t44;
t113 = -0.2e1 * t66;
t112 = 0.2e1 * t70;
t110 = pkin(2) * t67;
t109 = pkin(2) * t71;
t108 = pkin(4) * t69;
t107 = pkin(10) * t65;
t14 = t71 * t27 + t73 * t67;
t10 = -t41 * pkin(10) + t14;
t17 = -t61 * t30 + t63 * t34;
t8 = t28 * pkin(3) - t29 * pkin(10) + t17;
t5 = -t66 * t10 + t70 * t8;
t3 = -t28 * pkin(4) - t5;
t106 = t3 * t65;
t105 = t3 * t69;
t19 = t29 * t70 - t41 * t66;
t16 = t19 * t69 + t28 * t65;
t104 = t16 * t65;
t77 = pkin(9) * t97;
t39 = t77 + (pkin(10) + t110) * t63;
t40 = (-pkin(3) * t71 - pkin(10) * t67 - pkin(2)) * t61;
t23 = -t66 * t39 + t70 * t40;
t21 = pkin(4) * t97 - t23;
t103 = t21 * t65;
t102 = t21 * t69;
t32 = t69 * t44 - t65 * t97;
t101 = t32 * t65;
t56 = t61 ^ 2;
t100 = t56 * t71;
t57 = t62 ^ 2;
t99 = t57 * t72;
t94 = t63 * t67;
t92 = t65 * t18;
t43 = -t70 * t63 + t66 * t98;
t91 = t65 * t43;
t90 = t65 * t69;
t89 = t65 * t70;
t88 = t66 * t18;
t87 = t66 * t28;
t86 = t66 * t43;
t85 = t69 * t18;
t84 = t69 * t43;
t83 = t69 * t66;
t82 = t69 * t70;
t81 = t70 * t28;
t80 = 0.2e1 * t97;
t79 = 0.2e1 * t62 * t64;
t78 = t66 * t112;
t75 = t66 * t97;
t74 = t70 * t97;
t6 = t70 * t10 + t66 * t8;
t24 = t70 * t39 + t66 * t40;
t52 = pkin(9) * t98;
t38 = t52 + (-pkin(3) - t109) * t63;
t9 = t41 * pkin(3) - t13;
t60 = t69 ^ 2;
t59 = t66 ^ 2;
t58 = t65 ^ 2;
t50 = -t70 * pkin(4) - t66 * pkin(11) - pkin(3);
t47 = pkin(2) * t94 + t77;
t46 = -pkin(8) * t96 + t53;
t45 = t63 * t109 - t52;
t36 = pkin(10) * t82 + t65 * t50;
t35 = -pkin(10) * t89 + t69 * t50;
t22 = -pkin(11) * t97 + t24;
t20 = t43 * pkin(4) - t44 * pkin(11) + t38;
t15 = t19 * t65 - t28 * t69;
t12 = t65 * t20 + t69 * t22;
t11 = t69 * t20 - t65 * t22;
t7 = t18 * pkin(4) - t19 * pkin(11) + t9;
t4 = t28 * pkin(11) + t6;
t2 = t69 * t4 + t65 * t7;
t1 = -t65 * t4 + t69 * t7;
t25 = [1, 0, 0, t57 * t68 ^ 2, 0.2e1 * t68 * t99, t68 * t79, t72 * t79, t64 ^ 2, 0.2e1 * pkin(1) * t99 + 0.2e1 * t46 * t64, -0.2e1 * t57 * t111 - 0.2e1 * t48 * t64, t29 ^ 2, t28 * t116, t41 * t116, t41 * t117, t41 ^ 2, -0.2e1 * t13 * t41 + 0.2e1 * t17 * t28, 0.2e1 * t14 * t41 + 0.2e1 * t17 * t29, t19 ^ 2, t19 * t118, t19 * t117, t28 * t118, t28 ^ 2, 0.2e1 * t9 * t18 + 0.2e1 * t5 * t28, 0.2e1 * t9 * t19 - 0.2e1 * t6 * t28, t16 ^ 2, -0.2e1 * t16 * t15, 0.2e1 * t16 * t18, t15 * t118, t18 ^ 2, 0.2e1 * t1 * t18 + 0.2e1 * t3 * t15, 0.2e1 * t3 * t16 - 0.2e1 * t2 * t18; 0, 0, 0, 0, 0, t96, t95, t64, t46, -t48, t29 * t98, (-t28 * t67 + t29 * t71) * t61, t29 * t63 - t41 * t98, -t28 * t63 - t41 * t97, -t41 * t63, t13 * t63 - t45 * t41 + (-pkin(2) * t28 - t17 * t71) * t61, -t14 * t63 + t47 * t41 + (-pkin(2) * t29 + t17 * t67) * t61, t19 * t44, -t44 * t18 - t19 * t43, -t19 * t97 + t44 * t28, t18 * t97 - t43 * t28, -t28 * t97, t38 * t18 + t23 * t28 + t9 * t43 - t5 * t97, t38 * t19 - t24 * t28 + t9 * t44 + t6 * t97, t16 * t32, -t32 * t15 - t16 * t31, t16 * t43 + t32 * t18, -t15 * t43 - t31 * t18, t18 * t43, t1 * t43 + t11 * t18 + t21 * t15 + t3 * t31, -t12 * t18 + t21 * t16 - t2 * t43 + t3 * t32; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t56 * t67 ^ 2, 0.2e1 * t67 * t100, 0.2e1 * t61 * t94, t63 * t80, t63 ^ 2, 0.2e1 * pkin(2) * t100 + 0.2e1 * t45 * t63, -0.2e1 * t56 * t110 - 0.2e1 * t47 * t63, t44 ^ 2, t43 * t114, t97 * t114, t43 * t80, t56 * t71 ^ 2, -0.2e1 * t23 * t97 + 0.2e1 * t38 * t43, 0.2e1 * t24 * t97 + 0.2e1 * t38 * t44, t32 ^ 2, t32 * t115, 0.2e1 * t32 * t43, t43 * t115, t43 ^ 2, 0.2e1 * t11 * t43 + 0.2e1 * t21 * t31, -0.2e1 * t12 * t43 + 0.2e1 * t21 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t28, -t41, t13, -t14, t19 * t66, t19 * t70 - t88, t87, t81, 0, -pkin(3) * t18 - pkin(10) * t87 - t9 * t70, -pkin(3) * t19 - pkin(10) * t81 + t9 * t66, t16 * t83, (-t15 * t69 - t104) * t66, -t16 * t70 + t18 * t83, t15 * t70 - t65 * t88, -t18 * t70, -t1 * t70 + t35 * t18 + (pkin(10) * t15 + t106) * t66, -t36 * t18 + t2 * t70 + (pkin(10) * t16 + t105) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t97, t63, t45, -t47, t44 * t66, t44 * t70 - t86, -t75, -t74, 0, -pkin(3) * t43 + pkin(10) * t75 - t38 * t70, -pkin(3) * t44 + pkin(10) * t74 + t38 * t66, t32 * t83, (-t31 * t69 - t101) * t66, -t32 * t70 + t43 * t83, t31 * t70 - t65 * t86, -t43 * t70, -t11 * t70 + t35 * t43 + (pkin(10) * t31 + t103) * t66, t12 * t70 - t36 * t43 + (pkin(10) * t32 + t102) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t59, t78, 0, 0, 0, pkin(3) * t112, pkin(3) * t113, t60 * t59, -0.2e1 * t59 * t90, t82 * t113, t65 * t78, t70 ^ 2, 0.2e1 * t107 * t59 - 0.2e1 * t35 * t70, 0.2e1 * t59 * pkin(10) * t69 + 0.2e1 * t36 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, t28, t5, -t6, t104, -t65 * t15 + t16 * t69, t92, t85, 0, -pkin(4) * t15 - pkin(11) * t92 - t105, -pkin(4) * t16 - pkin(11) * t85 + t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t43, -t97, t23, -t24, t101, -t65 * t31 + t32 * t69, t91, t84, 0, -pkin(4) * t31 - pkin(11) * t91 - t102, -pkin(4) * t32 - pkin(11) * t84 + t103; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t70, 0, -t66 * pkin(10), -t70 * pkin(10), t65 * t83, (-t58 + t60) * t66, -t89, -t82, 0, -pkin(10) * t83 + (-pkin(4) * t66 + pkin(11) * t70) * t65, pkin(11) * t82 + (t107 - t108) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t58, 0.2e1 * t90, 0, 0, 0, 0.2e1 * t108, -0.2e1 * pkin(4) * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, t18, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, t43, t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t65 * t66, -t70, t35, -t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t69, 0, -t65 * pkin(11), -t69 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t25;
