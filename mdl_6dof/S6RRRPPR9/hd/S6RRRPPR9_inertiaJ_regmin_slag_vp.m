% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPPR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 06:38:53
% EndTime: 2019-05-07 06:38:58
% DurationCPUTime: 1.53s
% Computational Cost: add. (1668->198), mult. (3870->368), div. (0->0), fcn. (4322->10), ass. (0->107)
t78 = sin(pkin(11));
t80 = cos(pkin(11));
t127 = t78 ^ 2 + t80 ^ 2;
t82 = sin(qJ(6));
t85 = cos(qJ(6));
t53 = t85 * t78 - t82 * t80;
t79 = sin(pkin(6));
t87 = cos(qJ(2));
t107 = t79 * t87;
t84 = sin(qJ(2));
t108 = t79 * t84;
t81 = cos(pkin(6));
t83 = sin(qJ(3));
t86 = cos(qJ(3));
t47 = t86 * t108 + t81 * t83;
t27 = t80 * t107 + t47 * t78;
t28 = -t78 * t107 + t47 * t80;
t14 = t27 * t82 + t28 * t85;
t126 = -0.2e1 * t14;
t42 = t53 * t83;
t125 = 0.2e1 * t42;
t124 = -0.2e1 * t47;
t119 = pkin(4) + pkin(5);
t92 = t78 * qJ(5) + pkin(3);
t48 = t119 * t80 + t92;
t123 = 0.2e1 * t48;
t55 = -t80 * pkin(4) - t92;
t122 = -0.2e1 * t55;
t121 = 0.2e1 * t83;
t120 = 0.2e1 * t86;
t118 = pkin(1) * t84;
t117 = pkin(1) * t87;
t116 = pkin(3) * t78;
t115 = pkin(9) * t80;
t45 = t83 * t108 - t81 * t86;
t63 = pkin(8) * t108;
t38 = t63 + (-pkin(2) - t117) * t81;
t18 = t45 * pkin(3) - t47 * qJ(4) + t38;
t97 = pkin(8) * t107;
t39 = t97 + (pkin(9) + t118) * t81;
t41 = (-pkin(2) * t87 - pkin(9) * t84 - pkin(1)) * t79;
t22 = t86 * t39 + t83 * t41;
t19 = -qJ(4) * t107 + t22;
t9 = t78 * t18 + t80 * t19;
t6 = t45 * qJ(5) + t9;
t114 = t6 * t80;
t72 = t83 * pkin(9);
t113 = t86 * pkin(9);
t112 = t9 * t80;
t21 = -t83 * t39 + t86 * t41;
t20 = pkin(3) * t107 - t21;
t111 = t20 * t78;
t110 = t20 * t80;
t75 = t79 ^ 2;
t109 = t75 * t87;
t66 = t78 * t83;
t67 = t80 * t83;
t106 = t81 * t84;
t56 = -t86 * pkin(3) - t83 * qJ(4) - pkin(2);
t37 = t80 * t113 + t78 * t56;
t103 = t127 * qJ(4) ^ 2;
t102 = qJ(4) * t28;
t101 = qJ(4) * t86;
t100 = qJ(5) * t80;
t70 = t78 * qJ(4);
t99 = t80 * qJ(4);
t98 = 0.2e1 * t107;
t96 = t83 * t107;
t95 = t86 * t107;
t94 = t45 * t70;
t93 = t45 * t99;
t8 = t80 * t18 - t78 * t19;
t62 = t78 * t113;
t36 = t80 * t56 - t62;
t32 = -t86 * qJ(5) + t37;
t73 = t86 * pkin(4);
t33 = -t36 + t73;
t91 = t32 * t80 + t33 * t78;
t90 = -t36 * t78 + t37 * t80;
t52 = t82 * t78 + t85 * t80;
t89 = t28 * qJ(5) - t20;
t77 = t83 ^ 2;
t61 = t86 * t70;
t58 = (-pkin(10) + qJ(4)) * t80;
t57 = -t78 * pkin(10) + t70;
t54 = 0.2e1 * t127 * qJ(4);
t50 = pkin(1) * t106 + t97;
t49 = t81 * t117 - t63;
t43 = t52 * t83;
t40 = t72 + (pkin(4) * t78 - t100) * t83;
t31 = -t72 + (-t119 * t78 + t100) * t83;
t30 = t82 * t57 + t85 * t58;
t29 = t85 * t57 - t82 * t58;
t25 = pkin(10) * t66 + t32;
t24 = t27 * t99;
t23 = t86 * pkin(5) + t62 + t73 + (-pkin(10) * t83 - t56) * t80;
t13 = -t27 * t85 + t28 * t82;
t12 = t82 * t23 + t85 * t25;
t11 = t85 * t23 - t82 * t25;
t10 = t27 * pkin(4) - t89;
t7 = -t45 * pkin(4) - t8;
t5 = -t119 * t27 + t89;
t4 = t27 * pkin(10) + t6;
t3 = -t28 * pkin(10) - t119 * t45 - t8;
t2 = t82 * t3 + t85 * t4;
t1 = t85 * t3 - t82 * t4;
t15 = [1, 0, 0, t75 * t84 ^ 2, 0.2e1 * t84 * t109, 0.2e1 * t79 * t106, t81 * t98, t81 ^ 2, 0.2e1 * pkin(1) * t109 + 0.2e1 * t49 * t81, -0.2e1 * t75 * t118 - 0.2e1 * t50 * t81, t47 ^ 2, t45 * t124, t107 * t124, t45 * t98, t75 * t87 ^ 2, -0.2e1 * t21 * t107 + 0.2e1 * t38 * t45, 0.2e1 * t22 * t107 + 0.2e1 * t38 * t47, 0.2e1 * t20 * t27 + 0.2e1 * t8 * t45, 0.2e1 * t20 * t28 - 0.2e1 * t9 * t45, -0.2e1 * t9 * t27 - 0.2e1 * t8 * t28, t20 ^ 2 + t8 ^ 2 + t9 ^ 2, 0.2e1 * t10 * t27 - 0.2e1 * t7 * t45, -0.2e1 * t6 * t27 + 0.2e1 * t7 * t28, -0.2e1 * t10 * t28 + 0.2e1 * t6 * t45, t10 ^ 2 + t6 ^ 2 + t7 ^ 2, t14 ^ 2, t13 * t126, t45 * t126, 0.2e1 * t13 * t45, t45 ^ 2, -0.2e1 * t1 * t45 + 0.2e1 * t5 * t13, 0.2e1 * t5 * t14 + 0.2e1 * t2 * t45; 0, 0, 0, 0, 0, t108, t107, t81, t49, -t50, t47 * t83, -t83 * t45 + t47 * t86, -t96, -t95, 0, -pkin(2) * t45 + pkin(9) * t96 - t38 * t86, -pkin(2) * t47 + pkin(9) * t95 + t38 * t83, t36 * t45 - t8 * t86 + (pkin(9) * t27 + t111) * t83, -t37 * t45 + t9 * t86 + (pkin(9) * t28 + t110) * t83, -t37 * t27 - t36 * t28 + (-t78 * t9 - t8 * t80) * t83, t20 * t72 + t8 * t36 + t9 * t37, t10 * t66 + t40 * t27 - t33 * t45 + t7 * t86, -t32 * t27 + t33 * t28 + (-t6 * t78 + t7 * t80) * t83, -t10 * t67 - t40 * t28 + t32 * t45 - t6 * t86, t10 * t40 + t6 * t32 + t7 * t33, t14 * t43, -t43 * t13 + t14 * t42, t14 * t86 - t43 * t45, -t13 * t86 - t42 * t45, -t45 * t86, t1 * t86 - t11 * t45 + t31 * t13 - t5 * t42, t12 * t45 + t31 * t14 - t2 * t86 + t5 * t43; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t77, t83 * t120, 0, 0, 0, pkin(2) * t120, -0.2e1 * pkin(2) * t83, 0.2e1 * t77 * pkin(9) * t78 - 0.2e1 * t36 * t86, 0.2e1 * t77 * t115 + 0.2e1 * t37 * t86 (-t36 * t80 - t37 * t78) * t121, t77 * pkin(9) ^ 2 + t36 ^ 2 + t37 ^ 2, 0.2e1 * t33 * t86 + 0.2e1 * t40 * t66 (-t32 * t78 + t33 * t80) * t121, -0.2e1 * t32 * t86 - 0.2e1 * t40 * t67, t32 ^ 2 + t33 ^ 2 + t40 ^ 2, t43 ^ 2, t43 * t125, t43 * t120, t86 * t125, t86 ^ 2, 0.2e1 * t11 * t86 - 0.2e1 * t31 * t42, -0.2e1 * t12 * t86 + 0.2e1 * t31 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, -t45, -t107, t21, -t22, -pkin(3) * t27 - t110 - t94, -pkin(3) * t28 + t111 - t93, t112 - t24 + (-t8 + t102) * t78, -t20 * pkin(3) + (-t8 * t78 + t112) * qJ(4), -t10 * t80 + t55 * t27 - t94, t114 - t24 + (t7 + t102) * t78, -t10 * t78 - t55 * t28 + t93, t10 * t55 + (t7 * t78 + t114) * qJ(4), t14 * t53, -t53 * t13 - t14 * t52, -t53 * t45, t52 * t45, 0, t48 * t13 - t29 * t45 + t5 * t52, t48 * t14 + t30 * t45 + t5 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, t86, 0, -t72, -t113, t61 + (-t115 - t116) * t83, pkin(9) * t66 + (-pkin(3) * t83 + t101) * t80, t90, -pkin(3) * t72 + t90 * qJ(4), -t40 * t80 + t55 * t66 + t61, t91, -t40 * t78 + (-t55 * t83 - t101) * t80, qJ(4) * t91 + t40 * t55, t43 * t53, t53 * t42 - t43 * t52, t53 * t86, -t52 * t86, 0, t29 * t86 + t31 * t52 - t48 * t42, -t30 * t86 + t31 * t53 + t48 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t80, -0.2e1 * t116, t54, pkin(3) ^ 2 + t103, t80 * t122, t54, t78 * t122, t55 ^ 2 + t103, t53 ^ 2, -0.2e1 * t53 * t52, 0, 0, 0, t52 * t123, t53 * t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t28, 0, t20, t27, 0, -t28, t10, 0, 0, 0, 0, 0, -t13, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t66, t67, 0, t72, t66, 0, -t67, t40, 0, 0, 0, 0, 0, t42, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, t78, 0, -pkin(3), -t80, 0, -t78, t55, 0, 0, 0, 0, 0, -t52, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45, t28, 0, t7, 0, 0, 0, 0, 0, -t85 * t45, t82 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t67, 0, t33, 0, 0, 0, 0, 0, t85 * t86, -t82 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, 0, t70, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, -t45, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t42, t86, t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t52, 0, t29, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, -t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t15;
