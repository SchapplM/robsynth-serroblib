% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
% 
% Output:
% MM_reg [((6+1)*6/2)x35]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 22:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR5_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 20:30:53
% EndTime: 2019-05-07 20:30:58
% DurationCPUTime: 1.21s
% Computational Cost: add. (1006->137), mult. (1919->224), div. (0->0), fcn. (2223->8), ass. (0->110)
t84 = sin(qJ(4));
t81 = t84 ^ 2;
t88 = cos(qJ(4));
t82 = t88 ^ 2;
t102 = t81 + t82;
t85 = sin(qJ(3));
t119 = t85 * pkin(2);
t69 = pkin(9) + t119;
t104 = t102 * t69;
t132 = -t88 * pkin(4) - t84 * qJ(5);
t115 = cos(qJ(3));
t116 = cos(qJ(2));
t86 = sin(qJ(2));
t48 = -t115 * t116 + t85 * t86;
t131 = t48 ^ 2;
t83 = sin(qJ(6));
t87 = cos(qJ(6));
t47 = t83 * t84 + t87 * t88;
t130 = 0.2e1 * t47;
t129 = -0.2e1 * t48;
t128 = 0.2e1 * t48;
t108 = t87 * t84;
t50 = -t83 * t88 + t108;
t127 = 0.2e1 * t50;
t71 = -t116 * pkin(2) - pkin(1);
t126 = 0.2e1 * t71;
t125 = -0.2e1 * t84;
t124 = -0.2e1 * t88;
t123 = -pkin(4) - pkin(5);
t122 = pkin(9) * t48;
t121 = t48 * pkin(4);
t120 = t84 * pkin(10);
t118 = t88 * pkin(10);
t80 = t115 * pkin(2);
t70 = -t80 - pkin(3);
t117 = pkin(3) - t70;
t58 = (-pkin(8) - pkin(7)) * t86;
t97 = t116 * pkin(7);
t60 = t116 * pkin(8) + t97;
t31 = -t115 * t58 + t85 * t60;
t51 = t115 * t86 + t85 * t116;
t101 = t88 * qJ(5);
t92 = pkin(4) * t84 - t101;
t14 = t92 * t51 + t31;
t114 = t14 * t84;
t113 = t14 * t88;
t112 = t31 * t88;
t111 = t48 * t69;
t110 = t84 * t51;
t109 = t84 * t88;
t40 = t88 * t51;
t22 = t48 * pkin(3) - t51 * pkin(9) + t71;
t33 = t115 * t60 + t85 * t58;
t107 = -t88 * t22 + t84 * t33;
t13 = t84 * t22 + t88 * t33;
t98 = pkin(3) - t132;
t41 = -t80 - t98;
t78 = t88 * pkin(5);
t34 = t78 - t41;
t42 = t78 + t98;
t106 = t34 + t42;
t105 = -t41 + t98;
t103 = t102 * pkin(9);
t100 = 0.2e1 * t116;
t99 = t51 * t129;
t43 = t48 * qJ(5);
t9 = t43 + t13;
t62 = t84 * t69;
t96 = t62 - t120;
t75 = t84 * pkin(9);
t95 = t75 - t120;
t94 = -pkin(3) * t51 - t122;
t4 = -pkin(10) * t40 + t123 * t48 + t107;
t6 = pkin(10) * t110 + t9;
t1 = t87 * t4 - t83 * t6;
t2 = t83 * t4 + t87 * t6;
t93 = t51 * t98 + t122;
t10 = t107 - t121;
t3 = t10 * t84 + t9 * t88;
t91 = t41 * t51 - t111;
t90 = t51 * t70 - t111;
t77 = t88 * pkin(9);
t65 = 0.2e1 * t109;
t64 = t88 * t69;
t59 = t77 - t118;
t55 = t87 * qJ(5) + t83 * t123;
t54 = t83 * qJ(5) - t87 * t123;
t46 = t51 ^ 2;
t45 = t50 ^ 2;
t44 = t64 - t118;
t39 = t88 * t48;
t38 = t84 * t48;
t36 = t84 * t40;
t32 = t87 * t59 + t83 * t95;
t30 = t83 * t59 - t87 * t95;
t29 = t48 * t50;
t28 = t47 * t48;
t27 = -0.2e1 * t47 * t50;
t26 = t31 * t84;
t23 = (-t81 + t82) * t51;
t21 = t87 * t44 + t83 * t96;
t20 = t83 * t44 - t87 * t96;
t17 = t47 * t51;
t16 = -t51 * t108 + t83 * t40;
t15 = t17 * t50;
t11 = (t123 * t84 + t101) * t51 - t31;
t8 = t11 * t50;
t7 = t11 * t47;
t5 = -t16 * t50 - t47 * t17;
t12 = [1, 0, 0, t86 ^ 2, t86 * t100, 0, 0, 0, pkin(1) * t100, -0.2e1 * pkin(1) * t86, t46, t99, 0, 0, 0, t48 * t126, t51 * t126, t82 * t46, -0.2e1 * t46 * t109, t40 * t128, t84 * t99, t131, -0.2e1 * t107 * t48 + 0.2e1 * t31 * t110, -0.2e1 * t13 * t48 + 0.2e1 * t31 * t40, -0.2e1 * t10 * t48 + 0.2e1 * t14 * t110, 0.2e1 * (t10 * t88 - t84 * t9) * t51, -0.2e1 * t14 * t40 + 0.2e1 * t9 * t48, t10 ^ 2 + t14 ^ 2 + t9 ^ 2, t17 ^ 2, -0.2e1 * t16 * t17, t17 * t129, t16 * t128, t131, -0.2e1 * t1 * t48 + 0.2e1 * t11 * t16, 0.2e1 * t11 * t17 + 0.2e1 * t2 * t48; 0, 0, 0, 0, 0, t86, t116, 0, -t86 * pkin(7), -t97, 0, 0, t51, -t48, 0, -t31, -t33, t36, t23, t38, t39, 0, t90 * t84 - t112, t90 * t88 + t26, t91 * t84 - t113, t3, -t91 * t88 - t114, t14 * t41 + t3 * t69, t15, t5, -t29, t28, 0, t34 * t16 + t20 * t48 + t7, t34 * t17 + t21 * t48 + t8; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t80, -0.2e1 * t119, t81, t65, 0, 0, 0, t70 * t124, 0.2e1 * t70 * t84, t41 * t124, 0.2e1 * t104, t41 * t125, t102 * t69 ^ 2 + t41 ^ 2, t45, t27, 0, 0, 0, t34 * t130, t34 * t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t48, 0, -t31, -t33, t36, t23, t38, t39, 0, t94 * t84 - t112, t94 * t88 + t26, -t93 * t84 - t113, t3, t93 * t88 - t114, t3 * pkin(9) - t14 * t98, t15, t5, -t29, t28, 0, t42 * t16 + t30 * t48 + t7, t42 * t17 + t32 * t48 + t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t80, -t119, t81, t65, 0, 0, 0, t117 * t88, -t117 * t84, t105 * t88, t103 + t104, t105 * t84, pkin(9) * t104 - t41 * t98, t45, t27, 0, 0, 0, t106 * t47, t106 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t81, t65, 0, 0, 0, 0.2e1 * pkin(3) * t88, pkin(3) * t125, -t98 * t124, 0.2e1 * t103, -t98 * t125, t102 * pkin(9) ^ 2 + t98 ^ 2, t45, t27, 0, 0, 0, t42 * t130, t42 * t127; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t110, t48, -t107, -t13, -t107 + 0.2e1 * t121, t132 * t51, 0.2e1 * t43 + t13, -t10 * pkin(4) + t9 * qJ(5), 0, 0, -t17, t16, t48, t54 * t48 - t1, t55 * t48 + t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t88, 0, -t62, -t64, -t62, -t92, t64, -t92 * t69, 0, 0, -t50, t47, 0, t20, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t88, 0, -t75, -t77, -t75, -t92, t77, -t92 * pkin(9), 0, 0, -t50, t47, 0, t30, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t54, 0.2e1 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, t40, 0, t10, 0, 0, 0, 0, 0, -t87 * t48, t83 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, t62, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, 0, t75, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4), 0, 0, 0, 0, 0, -t87, t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, -t48, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t47, 0, -t20, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50, -t47, 0, -t30, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t54, -t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, -t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t12;
