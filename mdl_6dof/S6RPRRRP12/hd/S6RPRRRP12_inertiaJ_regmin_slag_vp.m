% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRP12_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP12_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 02:22:04
% EndTime: 2019-05-06 02:22:08
% DurationCPUTime: 1.49s
% Computational Cost: add. (3283->189), mult. (8815->363), div. (0->0), fcn. (10247->12), ass. (0->112)
t72 = cos(pkin(6));
t125 = pkin(1) * t72;
t67 = sin(pkin(12));
t70 = cos(pkin(12));
t69 = sin(pkin(6));
t96 = qJ(2) * t69;
t46 = t67 * t125 + t70 * t96;
t68 = sin(pkin(7));
t71 = cos(pkin(7));
t107 = t70 * t71;
t90 = t69 * t107;
t29 = (t68 * t72 + t90) * pkin(9) + t46;
t75 = sin(qJ(3));
t78 = cos(qJ(3));
t111 = t67 * t69;
t58 = t70 * t125;
t32 = t72 * pkin(2) + t58 + (-pkin(9) * t71 - qJ(2)) * t111;
t38 = (-pkin(9) * t67 * t68 - pkin(2) * t70 - pkin(1)) * t69;
t82 = t32 * t71 + t38 * t68;
t19 = -t75 * t29 + t78 * t82;
t110 = t68 * t75;
t31 = t72 * t110 + (t75 * t107 + t67 * t78) * t69;
t108 = t69 * t70;
t44 = -t68 * t108 + t72 * t71;
t74 = sin(qJ(4));
t77 = cos(qJ(4));
t24 = t31 * t74 - t44 * t77;
t131 = -0.2e1 * t24;
t109 = t68 * t78;
t30 = -t72 * t109 + t75 * t111 - t78 * t90;
t130 = -0.2e1 * t30;
t73 = sin(qJ(5));
t76 = cos(qJ(5));
t85 = -t76 * pkin(5) - t73 * qJ(6);
t52 = -pkin(4) + t85;
t129 = -0.2e1 * t52;
t128 = -0.2e1 * t74;
t127 = 0.2e1 * t77;
t63 = t69 ^ 2;
t126 = pkin(1) * t63;
t124 = pkin(4) * t73;
t123 = pkin(4) * t76;
t122 = pkin(10) * t73;
t121 = pkin(10) * t76;
t120 = t24 * pkin(5);
t23 = -t68 * t32 + t71 * t38;
t14 = t30 * pkin(3) - t31 * pkin(10) + t23;
t20 = t78 * t29 + t75 * t82;
t17 = t44 * pkin(10) + t20;
t9 = t77 * t14 - t74 * t17;
t7 = -t30 * pkin(4) - t9;
t119 = t7 * t73;
t118 = t7 * t76;
t117 = t73 * pkin(11);
t116 = t76 * pkin(11);
t16 = -t44 * pkin(3) - t19;
t25 = t31 * t77 + t44 * t74;
t13 = t24 * pkin(4) - t25 * pkin(11) + t16;
t10 = t74 * t14 + t77 * t17;
t8 = t30 * pkin(11) + t10;
t4 = t73 * t13 + t76 * t8;
t21 = t25 * t73 - t30 * t76;
t115 = t21 * t76;
t22 = t25 * t76 + t30 * t73;
t114 = t22 * t73;
t47 = t74 * t110 - t77 * t71;
t113 = t47 * t73;
t112 = t47 * t76;
t106 = t73 * t24;
t105 = t73 * t74;
t104 = t73 * t76;
t103 = t73 * t77;
t102 = t74 * t24;
t101 = t74 * t30;
t100 = t76 * t24;
t62 = t76 * t74;
t99 = t76 * t77;
t98 = t77 * t30;
t53 = -t77 * pkin(4) - t74 * pkin(11) - pkin(3);
t41 = pkin(10) * t99 + t73 * t53;
t64 = t73 ^ 2;
t66 = t76 ^ 2;
t97 = t64 + t66;
t95 = t24 * qJ(6);
t94 = t77 * qJ(6);
t93 = t74 * t127;
t92 = pkin(11) * t106;
t91 = pkin(11) * t100;
t89 = -t76 * t13 + t73 * t8;
t48 = t77 * t110 + t74 * t71;
t33 = t76 * t109 + t73 * t48;
t88 = t47 * t21 - t33 * t24;
t87 = t47 * t105 + t33 * t77;
t1 = t95 + t4;
t2 = t89 - t120;
t86 = t1 * t76 + t2 * t73;
t84 = -pkin(5) * t73 + t76 * qJ(6);
t34 = -t73 * t109 + t76 * t48;
t83 = t47 * t22 - t34 * t24;
t81 = t33 * t73 + t34 * t76;
t36 = -t94 + t41;
t50 = t76 * t53;
t37 = -t50 + (pkin(5) + t122) * t77;
t80 = t36 * t76 + t37 * t73;
t65 = t74 ^ 2;
t59 = pkin(11) * t103;
t45 = -t67 * t96 + t58;
t43 = (pkin(10) - t84) * t74;
t40 = -pkin(10) * t103 + t50;
t26 = t34 * t77 + t47 * t62;
t5 = t21 * pkin(5) - t22 * qJ(6) + t7;
t3 = [1, 0, 0, 0.2e1 * t70 * t126 + 0.2e1 * t45 * t72, -0.2e1 * t67 * t126 - 0.2e1 * t46 * t72, 0.2e1 * (-t45 * t67 + t46 * t70) * t69, t63 * pkin(1) ^ 2 + t45 ^ 2 + t46 ^ 2, t31 ^ 2, t31 * t130, 0.2e1 * t31 * t44, t44 * t130, t44 ^ 2, 0.2e1 * t19 * t44 + 0.2e1 * t23 * t30, -0.2e1 * t20 * t44 + 0.2e1 * t23 * t31, t25 ^ 2, t25 * t131, 0.2e1 * t25 * t30, t24 * t130, t30 ^ 2, 0.2e1 * t16 * t24 + 0.2e1 * t9 * t30, -0.2e1 * t10 * t30 + 0.2e1 * t16 * t25, t22 ^ 2, -0.2e1 * t22 * t21, 0.2e1 * t22 * t24, t21 * t131, t24 ^ 2, 0.2e1 * t7 * t21 - 0.2e1 * t24 * t89, 0.2e1 * t7 * t22 - 0.2e1 * t4 * t24, -0.2e1 * t2 * t24 + 0.2e1 * t5 * t21, -0.2e1 * t1 * t21 + 0.2e1 * t2 * t22, 0.2e1 * t1 * t24 - 0.2e1 * t5 * t22, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, -t108, t111, 0, -t69 * pkin(1), 0, 0, 0, 0, 0, t44 * t109 + t71 * t30, -t44 * t110 + t71 * t31, 0, 0, 0, 0, 0, -t24 * t109 - t47 * t30, -t25 * t109 - t48 * t30, 0, 0, 0, 0, 0, t88, t83, t88, -t34 * t21 + t33 * t22, -t83, t1 * t34 + t2 * t33 + t5 * t47; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 ^ 2 + t34 ^ 2 + t47 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, t44, t19, -t20, t25 * t74, t25 * t77 - t102, t101, t98, 0, -pkin(3) * t24 - pkin(10) * t101 - t16 * t77, -pkin(3) * t25 - pkin(10) * t98 + t16 * t74, t22 * t62 (-t114 - t115) * t74, -t22 * t77 + t24 * t62, -t102 * t73 + t21 * t77, -t24 * t77, t40 * t24 + t89 * t77 + (pkin(10) * t21 + t119) * t74, -t41 * t24 + t4 * t77 + (pkin(10) * t22 + t118) * t74, t105 * t5 + t2 * t77 + t43 * t21 - t37 * t24, -t36 * t21 + t37 * t22 + (-t1 * t73 + t2 * t76) * t74, -t1 * t77 - t43 * t22 + t36 * t24 - t5 * t62, t1 * t36 + t2 * t37 + t5 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, -t110, 0, 0, 0, 0, 0, t77 * t109, -t74 * t109, 0, 0, 0, 0, 0, t87, t26, t87 (t33 * t76 - t34 * t73) * t74, -t26, t33 * t37 + t34 * t36 + t47 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t65, t93, 0, 0, 0, pkin(3) * t127, pkin(3) * t128, t66 * t65, -0.2e1 * t65 * t104, t99 * t128, t73 * t93, t77 ^ 2, 0.2e1 * t65 * t122 - 0.2e1 * t40 * t77, 0.2e1 * t65 * t121 + 0.2e1 * t41 * t77, 0.2e1 * t105 * t43 + 0.2e1 * t37 * t77, 0.2e1 * (-t36 * t73 + t37 * t76) * t74, -0.2e1 * t36 * t77 - 0.2e1 * t43 * t62, t36 ^ 2 + t37 ^ 2 + t43 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, t30, t9, -t10, t114, -t73 * t21 + t22 * t76, t106, t100, 0, -pkin(4) * t21 - t118 - t92, -pkin(4) * t22 + t119 - t91, t52 * t21 - t5 * t76 - t92 (t114 - t115) * pkin(11) + t86, -t52 * t22 - t5 * t73 + t91, pkin(11) * t86 + t5 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t48, 0, 0, 0, 0, 0, -t112, t113, -t112, t81, -t113, pkin(11) * t81 + t47 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, t77, 0, -t74 * pkin(10), -t77 * pkin(10), t73 * t62 (-t64 + t66) * t74, -t103, -t99, 0, t59 + (-t121 - t124) * t74, pkin(11) * t99 + (t122 - t123) * t74, t105 * t52 - t43 * t76 + t59, t80, -t43 * t73 + (-pkin(11) * t77 - t52 * t74) * t76, pkin(11) * t80 + t43 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t64, 0.2e1 * t104, 0, 0, 0, 0.2e1 * t123, -0.2e1 * t124, t76 * t129, 0.2e1 * t97 * pkin(11), t73 * t129, pkin(11) ^ 2 * t97 + t52 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, t24, -t89, -t4, -t89 + 0.2e1 * t120, -t22 * pkin(5) - t21 * qJ(6), 0.2e1 * t95 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t34, -t33, 0, t34, -t33 * pkin(5) + t34 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t105, -t77, t40, -t41, t50 + (-0.2e1 * pkin(5) - t122) * t77, t85 * t74, -0.2e1 * t94 + t41, -t37 * pkin(5) + t36 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t76, 0, -t117, -t116, -t117, t84, t116, t84 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t22, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t62, 0, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
