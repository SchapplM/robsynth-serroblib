% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPR12_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR12_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRPR12_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t70 = sin(pkin(6));
t100 = cos(pkin(6));
t69 = sin(pkin(7));
t89 = t100 * t69;
t71 = cos(pkin(7));
t99 = cos(pkin(12));
t90 = t71 * t99;
t126 = t70 * t90 + t89;
t68 = sin(pkin(12));
t118 = pkin(1) * t68;
t91 = t70 * t99;
t45 = qJ(2) * t91 + t100 * t118;
t29 = pkin(9) * t126 + t45;
t74 = sin(qJ(3));
t77 = cos(qJ(3));
t114 = t68 * t70;
t93 = pkin(1) * t99;
t57 = t100 * t93;
t33 = t100 * pkin(2) + t57 + (-pkin(9) * t71 - qJ(2)) * t114;
t39 = (-pkin(9) * t68 * t69 - pkin(2) * t99 - pkin(1)) * t70;
t85 = t33 * t71 + t39 * t69;
t17 = -t29 * t74 + t77 * t85;
t32 = t74 * t89 + (t68 * t77 + t74 * t90) * t70;
t42 = -t100 * t71 + t69 * t91;
t73 = sin(qJ(4));
t76 = cos(qJ(4));
t27 = t32 * t76 - t42 * t73;
t125 = -0.2e1 * t27;
t31 = t74 * t114 - t126 * t77;
t124 = 0.2e1 * t31;
t123 = -0.2e1 * t32;
t122 = -0.2e1 * t73;
t121 = 0.2e1 * t76;
t120 = 2 * qJ(5);
t119 = pkin(4) + pkin(11);
t117 = t31 * pkin(4);
t26 = t32 * t73 + t42 * t76;
t72 = sin(qJ(6));
t75 = cos(qJ(6));
t20 = t26 * t72 + t31 * t75;
t116 = t20 * t75;
t115 = t27 * t72;
t23 = t27 * t73;
t113 = t69 * t74;
t112 = t69 * t77;
t111 = t72 * t73;
t110 = t72 * t76;
t109 = t72 * t119;
t108 = t73 * t31;
t107 = t73 * t76;
t106 = t75 * t72;
t105 = t75 * t76;
t104 = t75 * t119;
t103 = t76 * t31;
t21 = -t33 * t69 + t39 * t71;
t12 = pkin(3) * t31 - pkin(10) * t32 + t21;
t18 = t77 * t29 + t74 * t85;
t16 = -t42 * pkin(10) + t18;
t9 = t12 * t73 + t16 * t76;
t65 = t73 ^ 2;
t67 = t76 ^ 2;
t102 = t65 + t67;
t101 = qJ(5) * t76;
t30 = t31 * qJ(5);
t98 = -0.2e1 * t107;
t97 = pkin(10) * t108;
t96 = pkin(10) * t103;
t95 = t73 * t112;
t94 = t76 * t112;
t6 = -t30 - t9;
t92 = -t73 * qJ(5) - pkin(3);
t8 = t76 * t12 - t16 * t73;
t7 = -t8 - t117;
t87 = -t6 * t76 + t7 * t73;
t86 = -pkin(4) * t73 + t101;
t46 = t113 * t73 - t71 * t76;
t47 = t113 * t76 + t71 * t73;
t84 = t46 * t73 + t47 * t76;
t83 = -t119 * t73 + t101;
t82 = t112 * t26 + t31 * t46;
t81 = t112 * t27 + t31 * t47;
t15 = t42 * pkin(3) - t17;
t80 = -t27 * qJ(5) + t15;
t66 = t75 ^ 2;
t64 = t72 ^ 2;
t63 = t70 ^ 2;
t62 = t76 * pkin(10);
t61 = t73 * pkin(10);
t60 = t75 * t73;
t54 = pkin(5) * t76 + t62;
t53 = pkin(5) * t73 + t61;
t51 = -pkin(4) * t76 + t92;
t49 = -t119 * t76 + t92;
t44 = -qJ(2) * t114 + t57;
t37 = t112 * t75 - t46 * t72;
t36 = t112 * t72 + t46 * t75;
t35 = t49 * t75 + t53 * t72;
t34 = -t49 * t72 + t53 * t75;
t25 = t27 ^ 2;
t24 = t27 * t75;
t19 = -t26 * t75 + t31 * t72;
t10 = t26 * pkin(4) + t80;
t5 = t119 * t26 + t80;
t4 = -pkin(5) * t26 - t6;
t3 = t27 * pkin(5) - t119 * t31 - t8;
t2 = t3 * t72 + t5 * t75;
t1 = t3 * t75 - t5 * t72;
t11 = [1, 0, 0, 0.2e1 * t100 * t44 + 0.2e1 * t63 * t93, -0.2e1 * t100 * t45 - 0.2e1 * t118 * t63, 0.2e1 * (-t44 * t68 + t45 * t99) * t70, pkin(1) ^ 2 * t63 + t44 ^ 2 + t45 ^ 2, t32 ^ 2, t31 * t123, t42 * t123, t42 * t124, t42 ^ 2, -0.2e1 * t17 * t42 + 0.2e1 * t21 * t31, 0.2e1 * t18 * t42 + 0.2e1 * t21 * t32, t25, t26 * t125, t27 * t124, -0.2e1 * t26 * t31, t31 ^ 2, 0.2e1 * t15 * t26 + 0.2e1 * t31 * t8, 0.2e1 * t15 * t27 - 0.2e1 * t31 * t9, 0.2e1 * t26 * t6 + 0.2e1 * t27 * t7, -0.2e1 * t10 * t26 + 0.2e1 * t31 * t7, -0.2e1 * t10 * t27 - 0.2e1 * t31 * t6, t10 ^ 2 + t6 ^ 2 + t7 ^ 2, t20 ^ 2, -0.2e1 * t20 * t19, 0.2e1 * t20 * t27, t19 * t125, t25, 0.2e1 * t1 * t27 + 0.2e1 * t19 * t4, -0.2e1 * t2 * t27 + 0.2e1 * t20 * t4; 0, 0, 0, -t91, t114, 0, -t70 * pkin(1), 0, 0, 0, 0, 0, -t112 * t42 + t31 * t71, t113 * t42 + t32 * t71, 0, 0, 0, 0, 0, -t82, -t81, -t26 * t47 + t27 * t46, t82, t81, -t10 * t112 + t46 * t7 - t47 * t6, 0, 0, 0, 0, 0, t19 * t47 + t27 * t36, t20 * t47 + t27 * t37; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69 ^ 2 * t77 ^ 2 + t46 ^ 2 + t47 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, -t42, t17, -t18, t23, -t26 * t73 + t27 * t76, t108, t103, 0, -pkin(3) * t26 - t15 * t76 - t97, -pkin(3) * t27 + t15 * t73 - t96 (-t26 * t76 + t23) * pkin(10) + t87, t10 * t76 - t26 * t51 + t97, -t10 * t73 - t27 * t51 + t96, pkin(10) * t87 + t10 * t51, -t20 * t110 (t19 * t72 - t116) * t76, -t110 * t27 + t20 * t73, -t105 * t27 - t19 * t73, t23, t1 * t73 + t105 * t4 + t19 * t54 + t27 * t34, -t110 * t4 - t2 * t73 + t20 * t54 - t27 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, -t113, 0, 0, 0, 0, 0, t94, -t95, t84, -t94, t95, pkin(10) * t84 - t112 * t51, 0, 0, 0, 0, 0, t105 * t47 + t36 * t73, -t110 * t47 + t37 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t65, 0.2e1 * t107, 0, 0, 0, pkin(3) * t121, pkin(3) * t122, 0.2e1 * t102 * pkin(10), t51 * t121, t51 * t122, pkin(10) ^ 2 * t102 + t51 ^ 2, t64 * t67, 0.2e1 * t67 * t106, t72 * t98, t75 * t98, t65, 0.2e1 * t105 * t54 + 0.2e1 * t34 * t73, -0.2e1 * t110 * t54 - 0.2e1 * t35 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, t31, t8, -t9, -pkin(4) * t27 - qJ(5) * t26, -t8 - 0.2e1 * t117, -t6 + t30, -pkin(4) * t7 - qJ(5) * t6, t116, -t19 * t75 - t20 * t72, t24, -t115, 0, qJ(5) * t19 - t104 * t27 + t4 * t72, qJ(5) * t20 + t109 * t27 + t4 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t47, 0, t46, t47, -pkin(4) * t46 + qJ(5) * t47, 0, 0, 0, 0, 0, t47 * t72, t47 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, t76, 0, -t61, -t62, t86, t61, t62, t86 * pkin(10), -t72 * t105 (t64 - t66) * t76, t60, -t111, 0, t54 * t72 + t75 * t83, t54 * t75 - t72 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t120, pkin(4) ^ 2 + (qJ(5) ^ 2) t66, -0.2e1 * t106, 0, 0, 0, t72 * t120, t75 * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t31, 0, t7, 0, 0, 0, 0, 0, t24, -t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, 0, 0, t61, 0, 0, 0, 0, 0, t60, -t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, t27, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t110, -t105, t73, t34, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t72, 0, -t104, t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, -t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t11;
