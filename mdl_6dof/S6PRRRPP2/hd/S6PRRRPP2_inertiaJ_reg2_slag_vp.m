% Calculate inertial parameters regressor of joint inertia matrix for
% S6PRRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRPP2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP2_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRRPP2_inertiaJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t62 = sin(qJ(4));
t56 = t62 ^ 2;
t65 = cos(qJ(4));
t58 = t65 ^ 2;
t118 = t56 + t58;
t61 = cos(pkin(6));
t63 = sin(qJ(3));
t66 = cos(qJ(3));
t60 = sin(pkin(6));
t64 = sin(qJ(2));
t99 = t60 * t64;
t25 = t61 * t63 + t66 * t99;
t67 = cos(qJ(2));
t98 = t60 * t67;
t11 = t25 * t62 + t65 * t98;
t13 = t25 * t65 - t62 * t98;
t2 = t11 * t62 + t13 * t65;
t23 = -t61 * t66 + t63 * t99;
t114 = t23 ^ 2;
t117 = t11 ^ 2 + t13 ^ 2 + t114;
t1 = (t11 * t65 - t13 * t62) * t63;
t68 = pkin(4) + pkin(5);
t88 = t65 * qJ(5);
t116 = t68 * t62 - t88;
t89 = t62 * qJ(5);
t115 = t65 * t68 + t89;
t26 = pkin(3) + t115;
t113 = 0.2e1 * t26;
t80 = -t65 * pkin(4) - t89;
t31 = -pkin(3) + t80;
t112 = -0.2e1 * t31;
t111 = 0.2e1 * t63;
t110 = pkin(3) * t62;
t109 = pkin(3) * t65;
t57 = t63 ^ 2;
t108 = t57 * pkin(8);
t107 = t62 * pkin(9);
t106 = t63 * pkin(8);
t97 = t62 * t63;
t105 = t11 * t66 + t23 * t97;
t102 = t23 * t62;
t101 = t23 * t63;
t100 = t23 * t65;
t96 = t62 * t65;
t95 = t62 * t66;
t94 = t63 * t66;
t48 = t65 * t63;
t93 = t65 * t66;
t32 = -t66 * pkin(3) - t63 * pkin(9) - pkin(2);
t91 = pkin(8) * t95 - t65 * t32;
t19 = pkin(8) * t93 + t62 * t32;
t90 = t118 * pkin(9) ^ 2;
t87 = t65 * qJ(6);
t86 = t66 * qJ(5);
t85 = t62 * t94;
t84 = t57 * t96;
t83 = t63 * t93;
t54 = t66 * pkin(4);
t16 = t54 + t91;
t82 = -0.2e1 * t86 + t19;
t81 = t2 * pkin(9);
t15 = -t86 + t19;
t79 = -pkin(4) * t62 + t88;
t77 = t15 * t65 + t16 * t62;
t76 = t19 * t65 + t62 * t91;
t75 = t25 * t66 + t101;
t3 = t13 * t66 + t23 * t48;
t74 = t63 * t87 - t16;
t73 = pkin(8) ^ 2;
t71 = qJ(5) ^ 2;
t70 = 0.2e1 * qJ(5);
t59 = t66 ^ 2;
t55 = t60 ^ 2;
t53 = t65 * pkin(9);
t51 = t57 * t73;
t47 = t58 * t57;
t46 = t56 * t57;
t44 = t55 * t67 ^ 2;
t43 = -0.2e1 * t96;
t41 = pkin(9) * t95;
t39 = qJ(6) * t97;
t38 = t62 * t48;
t37 = -0.2e1 * t83;
t36 = 0.2e1 * t84;
t35 = 0.2e1 * t85;
t34 = t53 - t87;
t33 = (pkin(9) - qJ(6)) * t62;
t30 = 0.2e1 * t118 * pkin(9);
t29 = (t56 - t58) * t63;
t20 = (pkin(8) - t79) * t63;
t14 = (-pkin(8) - t116) * t63;
t8 = t13 * qJ(5);
t6 = t15 + t39;
t4 = t66 * pkin(5) - t74;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 * t64 ^ 2 + t61 ^ 2 + t44, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25 ^ 2 + t114 + t44, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, -t99, 0, 0, 0, 0, 0, 0, 0, 0, t66 * t98, -t63 * t98, t75, pkin(2) * t98 + pkin(8) * t75, 0, 0, 0, 0, 0, 0, t105, t3, t1, pkin(8) * t101 + t11 * t91 + t13 * t19, 0, 0, 0, 0, 0, 0, t105, t1, -t3, t11 * t16 + t13 * t15 + t23 * t20, 0, 0, 0, 0, 0, 0, t105, -t3, -t1, t11 * t4 + t13 * t6 - t23 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t57, 0.2e1 * t94, 0, t59, 0, 0, 0.2e1 * pkin(2) * t66, -0.2e1 * pkin(2) * t63, 0.2e1 * (t57 + t59) * pkin(8), pkin(2) ^ 2 + t59 * t73 + t51, t47, -0.2e1 * t84, t37, t46, t35, t59, 0.2e1 * t62 * t108 + 0.2e1 * t66 * t91, 0.2e1 * t65 * t108 + 0.2e1 * t19 * t66 (-t19 * t62 + t65 * t91) * t111, t19 ^ 2 + t91 ^ 2 + t51, t47, t37, t36, t59, -0.2e1 * t85, t46, 0.2e1 * t16 * t66 + 0.2e1 * t20 * t97 (-t15 * t62 + t16 * t65) * t111, -0.2e1 * t15 * t66 - 0.2e1 * t20 * t48, t15 ^ 2 + t16 ^ 2 + t20 ^ 2, t47, t36, 0.2e1 * t83, t46, t35, t59, -0.2e1 * t14 * t97 + 0.2e1 * t4 * t66, 0.2e1 * t14 * t48 - 0.2e1 * t6 * t66 (-t4 * t65 + t6 * t62) * t111, t14 ^ 2 + t4 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t25, 0, 0, 0, 0, 0, 0, 0, 0, -t100, t102, t2, -t23 * pkin(3) + t81, 0, 0, 0, 0, 0, 0, -t100, t2, -t102, t23 * t31 + t81, 0, 0, 0, 0, 0, 0, -t100, -t102, -t2, t11 * t33 + t13 * t34 - t23 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, t66, 0, -t106, -t66 * pkin(8), 0, 0, t38, -t29, -t95, -t38, -t93, 0, t41 + (-pkin(8) * t65 - t110) * t63, pkin(9) * t93 + (pkin(8) * t62 - t109) * t63, t76, -pkin(3) * t106 + pkin(9) * t76, t38, -t95, t29, 0, t93, -t38, -t20 * t65 + t31 * t97 + t41, t77, -t20 * t62 + (-pkin(9) * t66 - t31 * t63) * t65, pkin(9) * t77 + t20 * t31, t38, t29, t95, -t38, -t93, 0, t14 * t65 - t26 * t97 + t33 * t66, t14 * t62 + t26 * t48 - t34 * t66 (-t33 * t63 - t6) * t65 + (t34 * t63 - t4) * t62, t14 * t26 + t4 * t33 + t6 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t56, 0.2e1 * t96, 0, t58, 0, 0, 0.2e1 * t109, -0.2e1 * t110, t30, pkin(3) ^ 2 + t90, t56, 0, t43, 0, 0, t58, t65 * t112, t30, t62 * t112, t31 ^ 2 + t90, t56, t43, 0, t58, 0, 0, t65 * t113, t62 * t113, -0.2e1 * t33 * t62 - 0.2e1 * t34 * t65, t26 ^ 2 + t33 ^ 2 + t34 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t13, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, t13, -t11 * pkin(4) + t8, 0, 0, 0, 0, 0, 0, -t11, t13, 0, -t11 * t68 + t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, 0, -t97, -t66, -t91, -t19, 0, 0, 0, t48, 0, -t66, t97, 0, -0.2e1 * t54 - t91, t80 * t63, t82, -t16 * pkin(4) + t15 * qJ(5), 0, 0, -t48, 0, -t97, -t66 (-pkin(5) - t68) * t66 + t74, t39 + t82, t115 * t63, t6 * qJ(5) - t4 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, t65, 0, -t107, -t53, 0, 0, 0, t62, 0, 0, -t65, 0, -t107, t79, t53, t79 * pkin(9), 0, 0, -t62, 0, t65, 0, -t33, t34, t116, t34 * qJ(5) - t33 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, t70, pkin(4) ^ 2 + t71, 0, 0, 0, 0, 0, 1, 0.2e1 * t68, t70, 0, t68 ^ 2 + t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t48, 0, t16, 0, 0, 0, 0, 0, 0, t66, 0, -t48, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, t107, 0, 0, 0, 0, 0, 0, 0, 0, -t62, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), 0, 0, 0, 0, 0, 0, -1, 0, 0, -t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, t48, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t62, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t5;
