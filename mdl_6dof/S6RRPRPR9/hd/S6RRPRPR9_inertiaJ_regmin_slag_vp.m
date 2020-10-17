% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 11:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR9_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR9_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 15:18:28
% EndTime: 2019-05-06 15:18:33
% DurationCPUTime: 1.36s
% Computational Cost: add. (2510->170), mult. (5740->352), div. (0->0), fcn. (6770->12), ass. (0->101)
t116 = cos(qJ(6));
t82 = sin(pkin(6));
t89 = cos(qJ(2));
t110 = t82 * t89;
t117 = cos(qJ(4));
t88 = sin(qJ(2));
t111 = t82 * t88;
t81 = sin(pkin(11));
t84 = cos(pkin(11));
t85 = cos(pkin(6));
t55 = t84 * t111 + t81 * t85;
t87 = sin(qJ(4));
t67 = t81 * t111;
t99 = t84 * t85 - t67;
t35 = t117 * t55 + t87 * t99;
t80 = sin(pkin(12));
t83 = cos(pkin(12));
t29 = t83 * t110 + t35 * t80;
t30 = -t80 * t110 + t35 * t83;
t86 = sin(qJ(6));
t18 = t116 * t29 + t30 * t86;
t124 = -0.2e1 * t18;
t123 = -0.2e1 * t35;
t60 = -t117 * t84 + t81 * t87;
t122 = -0.2e1 * t60;
t73 = -pkin(5) * t83 - pkin(4);
t121 = 0.2e1 * t73;
t74 = -t84 * pkin(3) - pkin(2);
t120 = 0.2e1 * t74;
t119 = pkin(1) * t88;
t118 = t89 * pkin(1);
t34 = -t117 * t99 + t55 * t87;
t61 = t116 * t80 + t86 * t83;
t115 = t61 * t34;
t114 = t61 * t60;
t77 = t82 ^ 2;
t113 = t77 * t89;
t62 = t117 * t81 + t87 * t84;
t112 = t80 * t62;
t109 = t83 * t62;
t108 = t85 * t88;
t107 = pkin(9) + qJ(3);
t106 = pkin(10) + qJ(5);
t101 = pkin(8) * t110;
t51 = t101 + (qJ(3) + t119) * t85;
t52 = (-pkin(2) * t89 - qJ(3) * t88 - pkin(1)) * t82;
t31 = -t51 * t81 + t84 * t52;
t22 = -pkin(3) * t110 - pkin(9) * t55 + t31;
t32 = t84 * t51 + t81 * t52;
t26 = t99 * pkin(9) + t32;
t14 = t117 * t26 + t87 * t22;
t11 = -qJ(5) * t110 + t14;
t69 = pkin(8) * t111;
t40 = t67 * pkin(3) + t69 + (t74 - t118) * t85;
t17 = t34 * pkin(4) - t35 * qJ(5) + t40;
t6 = t83 * t11 + t80 * t17;
t41 = pkin(4) * t60 - qJ(5) * t62 + t74;
t64 = t107 * t81;
t66 = t107 * t84;
t46 = t117 * t66 - t87 * t64;
t24 = t80 * t41 + t83 * t46;
t105 = t80 ^ 2 + t83 ^ 2;
t104 = t81 ^ 2 + t84 ^ 2;
t103 = qJ(5) * t34;
t102 = 0.2e1 * t110;
t100 = qJ(3) * t110;
t5 = -t11 * t80 + t83 * t17;
t23 = t83 * t41 - t46 * t80;
t44 = t117 * t64 + t66 * t87;
t98 = t5 * t83 + t6 * t80;
t97 = -t5 * t80 + t6 * t83;
t13 = t117 * t22 - t87 * t26;
t96 = -pkin(4) * t62 - qJ(5) * t60;
t95 = t23 * t83 + t24 * t80;
t94 = -t23 * t80 + t24 * t83;
t93 = -t31 * t81 + t32 * t84;
t12 = pkin(4) * t110 - t13;
t92 = t116 * t83 - t80 * t86;
t65 = t106 * t83;
t63 = t106 * t80;
t57 = pkin(1) * t108 + t101;
t56 = t85 * t118 - t69;
t54 = t69 + (-pkin(2) - t118) * t85;
t47 = t92 * t60;
t45 = t116 * t65 - t86 * t63;
t43 = -t116 * t63 - t86 * t65;
t37 = t92 * t62;
t36 = t61 * t62;
t33 = pkin(5) * t112 + t44;
t28 = t92 * t34;
t21 = -pkin(10) * t112 + t24;
t20 = pkin(5) * t60 - pkin(10) * t109 + t23;
t19 = t116 * t30 - t86 * t29;
t9 = t116 * t21 + t86 * t20;
t8 = t116 * t20 - t86 * t21;
t7 = t29 * pkin(5) + t12;
t4 = -pkin(10) * t29 + t6;
t3 = pkin(5) * t34 - pkin(10) * t30 + t5;
t2 = t116 * t4 + t86 * t3;
t1 = t116 * t3 - t86 * t4;
t10 = [1, 0, 0, t77 * t88 ^ 2, 0.2e1 * t88 * t113, 0.2e1 * t82 * t108, t85 * t102, t85 ^ 2, 0.2e1 * pkin(1) * t113 + 0.2e1 * t56 * t85, -0.2e1 * t77 * t119 - 0.2e1 * t57 * t85, -0.2e1 * t31 * t110 - 0.2e1 * t54 * t99, 0.2e1 * t32 * t110 + 0.2e1 * t54 * t55, -0.2e1 * t31 * t55 + 0.2e1 * t32 * t99, t31 ^ 2 + t32 ^ 2 + t54 ^ 2, t35 ^ 2, t34 * t123, t110 * t123, t34 * t102, t77 * t89 ^ 2, -0.2e1 * t13 * t110 + 0.2e1 * t34 * t40, 0.2e1 * t14 * t110 + 0.2e1 * t35 * t40, 0.2e1 * t12 * t29 + 0.2e1 * t34 * t5, 0.2e1 * t12 * t30 - 0.2e1 * t34 * t6, -0.2e1 * t29 * t6 - 0.2e1 * t30 * t5, t12 ^ 2 + t5 ^ 2 + t6 ^ 2, t19 ^ 2, t19 * t124, 0.2e1 * t19 * t34, t34 * t124, t34 ^ 2, 0.2e1 * t1 * t34 + 0.2e1 * t18 * t7, 0.2e1 * t19 * t7 - 0.2e1 * t2 * t34; 0, 0, 0, 0, 0, t111, t110, t85, t56, -t57, pkin(2) * t99 + t81 * t100 - t54 * t84, -pkin(2) * t55 + t84 * t100 + t54 * t81 (t81 * t55 + t84 * t99) * qJ(3) + t93, -pkin(2) * t54 + t93 * qJ(3), t35 * t62, -t34 * t62 - t35 * t60, -t62 * t110, t60 * t110, 0, t44 * t110 + t34 * t74 + t40 * t60, t46 * t110 + t35 * t74 + t40 * t62, t12 * t112 + t23 * t34 + t29 * t44 + t5 * t60, t12 * t109 - t24 * t34 + t30 * t44 - t6 * t60, -t23 * t30 - t24 * t29 - t98 * t62, t12 * t44 + t23 * t5 + t24 * t6, t19 * t37, -t18 * t37 - t19 * t36, t19 * t60 + t34 * t37, -t18 * t60 - t34 * t36, t34 * t60, t1 * t60 + t18 * t33 + t34 * t8 + t36 * t7, t19 * t33 - t2 * t60 - t34 * t9 + t37 * t7; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t84, -0.2e1 * pkin(2) * t81, 0.2e1 * t104 * qJ(3), t104 * qJ(3) ^ 2 + pkin(2) ^ 2, t62 ^ 2, t62 * t122, 0, 0, 0, t60 * t120, t62 * t120, 0.2e1 * t44 * t112 + 0.2e1 * t23 * t60, 0.2e1 * t44 * t109 - 0.2e1 * t24 * t60, -0.2e1 * t95 * t62, t23 ^ 2 + t24 ^ 2 + t44 ^ 2, t37 ^ 2, -0.2e1 * t37 * t36, 0.2e1 * t37 * t60, t36 * t122, t60 ^ 2, 0.2e1 * t33 * t36 + 0.2e1 * t60 * t8, 0.2e1 * t33 * t37 - 0.2e1 * t60 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, t55, 0, t54, 0, 0, 0, 0, 0, t34, t35, t83 * t34, -t80 * t34, -t29 * t80 - t30 * t83, t98, 0, 0, 0, 0, 0, t28, -t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, t81, 0, -pkin(2), 0, 0, 0, 0, 0, t60, t62, t83 * t60, -t80 * t60, -t105 * t62, t95, 0, 0, 0, 0, 0, t47, -t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t34, -t110, t13, -t14, -pkin(4) * t29 - t80 * t103 - t12 * t83, -pkin(4) * t30 - t83 * t103 + t12 * t80 (-t29 * t83 + t30 * t80) * qJ(5) + t97, -pkin(4) * t12 + qJ(5) * t97, t19 * t61, -t18 * t61 + t19 * t92, t115, t28, 0, t18 * t73 + t34 * t43 - t7 * t92, t19 * t73 - t34 * t45 + t61 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t60, 0, -t44, -t46, -t44 * t83 + t96 * t80, t44 * t80 + t96 * t83, t94, -pkin(4) * t44 + qJ(5) * t94, t37 * t61, -t36 * t61 + t37 * t92, t114, t47, 0, -t33 * t92 + t36 * t73 + t43 * t60, t33 * t61 + t37 * t73 - t45 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4) * t83, -0.2e1 * pkin(4) * t80, 0.2e1 * t105 * qJ(5), t105 * qJ(5) ^ 2 + pkin(4) ^ 2, t61 ^ 2, 0.2e1 * t61 * t92, 0, 0, 0, -t92 * t121, t61 * t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t30, 0, t12, 0, 0, 0, 0, 0, t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, t109, 0, t44, 0, 0, 0, 0, 0, t36, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, t80, 0, -pkin(4), 0, 0, 0, 0, 0, -t92, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t18, t34, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t36, t60, t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t92, -t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, t92, 0, t43, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t10;
