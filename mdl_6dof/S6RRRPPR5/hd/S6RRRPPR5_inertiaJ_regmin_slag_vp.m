% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 15:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRPPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPPR5_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 05:06:18
% EndTime: 2019-05-07 05:06:22
% DurationCPUTime: 1.20s
% Computational Cost: add. (2587->168), mult. (5841->350), div. (0->0), fcn. (6840->12), ass. (0->103)
t107 = -qJ(4) - pkin(9);
t88 = sin(qJ(3));
t100 = t107 * t88;
t90 = cos(qJ(3));
t67 = t107 * t90;
t82 = sin(pkin(11));
t85 = cos(pkin(11));
t48 = -t85 * t100 - t67 * t82;
t125 = t48 ^ 2;
t117 = cos(qJ(6));
t83 = sin(pkin(6));
t91 = cos(qJ(2));
t110 = t83 * t91;
t89 = sin(qJ(2));
t111 = t83 * t89;
t86 = cos(pkin(6));
t56 = t90 * t111 + t86 * t88;
t69 = t88 * t111;
t99 = t86 * t90 - t69;
t39 = t85 * t56 + t82 * t99;
t81 = sin(pkin(12));
t84 = cos(pkin(12));
t31 = t84 * t110 + t81 * t39;
t32 = -t81 * t110 + t39 * t84;
t87 = sin(qJ(6));
t18 = t117 * t31 + t87 * t32;
t124 = -0.2e1 * t18;
t63 = t82 * t90 + t85 * t88;
t64 = t117 * t81 + t87 * t84;
t36 = t64 * t63;
t123 = -0.2e1 * t36;
t76 = -pkin(3) * t85 - pkin(4);
t66 = -t84 * pkin(5) + t76;
t122 = 0.2e1 * t66;
t121 = 0.2e1 * t90;
t120 = pkin(1) * t89;
t119 = t91 * pkin(1);
t73 = pkin(3) * t82 + qJ(5);
t118 = pkin(10) + t73;
t38 = t56 * t82 - t85 * t99;
t116 = t64 * t38;
t61 = t82 * t88 - t85 * t90;
t115 = t64 * t61;
t79 = t83 ^ 2;
t114 = t79 * t91;
t113 = t81 * t38;
t112 = t81 * t63;
t109 = t84 * t38;
t108 = t84 * t63;
t103 = pkin(8) * t110;
t53 = t103 + (pkin(9) + t120) * t86;
t54 = (-pkin(2) * t91 - pkin(9) * t89 - pkin(1)) * t83;
t33 = -t88 * t53 + t90 * t54;
t23 = -pkin(3) * t110 - t56 * qJ(4) + t33;
t34 = t90 * t53 + t88 * t54;
t28 = t99 * qJ(4) + t34;
t14 = t82 * t23 + t85 * t28;
t11 = -qJ(5) * t110 + t14;
t70 = pkin(8) * t111;
t77 = -t90 * pkin(3) - pkin(2);
t44 = t69 * pkin(3) + t70 + (t77 - t119) * t86;
t17 = t38 * pkin(4) - t39 * qJ(5) + t44;
t6 = t84 * t11 + t81 * t17;
t45 = t61 * pkin(4) - t63 * qJ(5) + t77;
t50 = t82 * t100 - t85 * t67;
t25 = t81 * t45 + t84 * t50;
t106 = t81 ^ 2 + t84 ^ 2;
t105 = 0.2e1 * t83 * t86;
t104 = -0.2e1 * t110;
t102 = t88 * t110;
t101 = t90 * t110;
t5 = -t11 * t81 + t84 * t17;
t13 = t85 * t23 - t82 * t28;
t24 = t84 * t45 - t50 * t81;
t12 = pkin(4) * t110 - t13;
t98 = t5 * t84 + t6 * t81;
t97 = -t5 * t81 + t6 * t84;
t96 = t24 * t84 + t25 * t81;
t95 = -t24 * t81 + t25 * t84;
t94 = -t61 * t73 + t63 * t76;
t93 = t117 * t84 - t87 * t81;
t60 = t86 * t120 + t103;
t59 = t86 * t119 - t70;
t58 = t118 * t84;
t57 = t118 * t81;
t52 = t70 + (-pkin(2) - t119) * t86;
t47 = t93 * t61;
t43 = t117 * t58 - t87 * t57;
t42 = -t117 * t57 - t87 * t58;
t37 = t93 * t63;
t35 = pkin(5) * t112 + t48;
t30 = t93 * t38;
t21 = -pkin(10) * t112 + t25;
t20 = pkin(5) * t61 - pkin(10) * t108 + t24;
t19 = t117 * t32 - t87 * t31;
t9 = t117 * t21 + t87 * t20;
t8 = t117 * t20 - t87 * t21;
t7 = t31 * pkin(5) + t12;
t4 = -pkin(10) * t31 + t6;
t3 = pkin(5) * t38 - pkin(10) * t32 + t5;
t2 = t117 * t4 + t87 * t3;
t1 = t117 * t3 - t87 * t4;
t10 = [1, 0, 0, t79 * t89 ^ 2, 0.2e1 * t89 * t114, t89 * t105, t91 * t105, t86 ^ 2, 0.2e1 * pkin(1) * t114 + 0.2e1 * t59 * t86, -0.2e1 * t79 * t120 - 0.2e1 * t60 * t86, t56 ^ 2, 0.2e1 * t56 * t99, t56 * t104, t99 * t104, t79 * t91 ^ 2, -0.2e1 * t33 * t110 - 0.2e1 * t52 * t99, 0.2e1 * t34 * t110 + 0.2e1 * t52 * t56, -0.2e1 * t13 * t39 - 0.2e1 * t14 * t38, t13 ^ 2 + t14 ^ 2 + t44 ^ 2, 0.2e1 * t12 * t31 + 0.2e1 * t38 * t5, 0.2e1 * t12 * t32 - 0.2e1 * t38 * t6, -0.2e1 * t31 * t6 - 0.2e1 * t32 * t5, t12 ^ 2 + t5 ^ 2 + t6 ^ 2, t19 ^ 2, t19 * t124, 0.2e1 * t19 * t38, t38 * t124, t38 ^ 2, 0.2e1 * t1 * t38 + 0.2e1 * t18 * t7, 0.2e1 * t19 * t7 - 0.2e1 * t2 * t38; 0, 0, 0, 0, 0, t111, t110, t86, t59, -t60, t56 * t88, t56 * t90 + t88 * t99, -t102, -t101, 0, pkin(2) * t99 + pkin(9) * t102 - t52 * t90, -pkin(2) * t56 + pkin(9) * t101 + t52 * t88, -t13 * t63 - t14 * t61 - t38 * t50 + t39 * t48, -t13 * t48 + t14 * t50 + t44 * t77, t12 * t112 + t24 * t38 + t31 * t48 + t5 * t61, t12 * t108 - t25 * t38 + t32 * t48 - t6 * t61, -t24 * t32 - t25 * t31 - t98 * t63, t12 * t48 + t24 * t5 + t25 * t6, t19 * t37, -t18 * t37 - t19 * t36, t19 * t61 + t37 * t38, -t18 * t61 - t36 * t38, t38 * t61, t1 * t61 + t18 * t35 + t36 * t7 + t38 * t8, t19 * t35 - t2 * t61 + t37 * t7 - t38 * t9; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t88 ^ 2, t88 * t121, 0, 0, 0, pkin(2) * t121, -0.2e1 * pkin(2) * t88, 0.2e1 * t48 * t63 - 0.2e1 * t50 * t61, t50 ^ 2 + t77 ^ 2 + t125, 0.2e1 * t48 * t112 + 0.2e1 * t24 * t61, 0.2e1 * t48 * t108 - 0.2e1 * t25 * t61, -0.2e1 * t96 * t63, t24 ^ 2 + t25 ^ 2 + t125, t37 ^ 2, t37 * t123, 0.2e1 * t37 * t61, t61 * t123, t61 ^ 2, 0.2e1 * t35 * t36 + 0.2e1 * t61 * t8, 0.2e1 * t35 * t37 - 0.2e1 * t9 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, t99, -t110, t33, -t34 (-t38 * t82 - t39 * t85) * pkin(3) (t13 * t85 + t14 * t82) * pkin(3), -t73 * t113 - t12 * t84 + t31 * t76, -t73 * t109 + t12 * t81 + t32 * t76 (-t31 * t84 + t32 * t81) * t73 + t97, t12 * t76 + t97 * t73, t19 * t64, -t18 * t64 + t19 * t93, t116, t30, 0, t18 * t66 + t38 * t42 - t7 * t93, t19 * t66 - t38 * t43 + t64 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t90, 0, -t88 * pkin(9), -t90 * pkin(9) (-t61 * t82 - t63 * t85) * pkin(3) (-t48 * t85 + t50 * t82) * pkin(3), -t48 * t84 + t94 * t81, t48 * t81 + t94 * t84, t95, t48 * t76 + t95 * t73, t37 * t64, -t36 * t64 + t37 * t93, t115, t47, 0, -t35 * t93 + t36 * t66 + t42 * t61, t35 * t64 + t37 * t66 - t43 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t82 ^ 2 + t85 ^ 2) * pkin(3) ^ 2, -0.2e1 * t76 * t84, 0.2e1 * t76 * t81, 0.2e1 * t106 * t73, t106 * t73 ^ 2 + t76 ^ 2, t64 ^ 2, 0.2e1 * t64 * t93, 0, 0, 0, -t93 * t122, t64 * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, t109, -t113, -t31 * t81 - t84 * t32, t98, 0, 0, 0, 0, 0, t30, -t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, t84 * t61, -t81 * t61, -t106 * t63, t96, 0, 0, 0, 0, 0, t47, -t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t106, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t32, 0, t12, 0, 0, 0, 0, 0, t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, t108, 0, t48, 0, 0, 0, 0, 0, t36, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, t81, 0, t76, 0, 0, 0, 0, 0, -t93, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, t38, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t36, t61, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, t93, 0, t42, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, -t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t10;
