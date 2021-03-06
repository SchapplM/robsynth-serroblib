% Calculate inertial parameters regressor of joint inertia matrix for
% S6PRRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 23:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRRRPR2_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_inertiaJ_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 07:16:24
% EndTime: 2019-05-05 07:16:32
% DurationCPUTime: 2.16s
% Computational Cost: add. (1636->175), mult. (3459->313), div. (0->0), fcn. (4188->12), ass. (0->115)
t95 = sin(pkin(12));
t90 = t95 ^ 2;
t97 = cos(pkin(12));
t92 = t97 ^ 2;
t119 = t90 + t92;
t120 = t119 * qJ(5);
t103 = cos(qJ(6));
t99 = sin(qJ(6));
t144 = t103 * t97 - t95 * t99;
t100 = sin(qJ(4));
t104 = cos(qJ(4));
t101 = sin(qJ(3));
t105 = cos(qJ(3));
t102 = sin(qJ(2));
t96 = sin(pkin(6));
t116 = t96 * t102;
t98 = cos(pkin(6));
t57 = -t101 * t116 + t105 * t98;
t58 = t101 * t98 + t105 * t116;
t30 = t100 * t58 - t104 * t57;
t29 = t30 ^ 2;
t137 = -pkin(9) - pkin(8);
t113 = t137 * t101;
t77 = t137 * t105;
t49 = -t100 * t77 - t104 * t113;
t143 = t49 ^ 2;
t68 = t100 * t101 - t104 * t105;
t66 = t68 ^ 2;
t142 = 0.2e1 * t68;
t135 = pkin(5) * t97;
t132 = t104 * pkin(3);
t85 = -pkin(4) - t132;
t72 = t85 - t135;
t141 = 0.2e1 * t72;
t83 = -pkin(4) - t135;
t140 = 0.2e1 * t83;
t86 = -pkin(3) * t105 - pkin(2);
t139 = 0.2e1 * t86;
t138 = 0.2e1 * t105;
t70 = t100 * t105 + t101 * t104;
t125 = t97 * t70;
t40 = pkin(4) * t68 - qJ(5) * t70 + t86;
t51 = t100 * t113 - t104 * t77;
t13 = t97 * t40 - t51 * t95;
t11 = pkin(5) * t68 - pkin(10) * t125 + t13;
t128 = t95 * t70;
t14 = t95 * t40 + t97 * t51;
t12 = -pkin(10) * t128 + t14;
t4 = t103 * t11 - t12 * t99;
t5 = t103 * t12 + t11 * t99;
t65 = t103 * t95 + t97 * t99;
t136 = t144 * t5 - t4 * t65;
t134 = pkin(4) - t85;
t133 = t100 * pkin(3);
t131 = t30 * t49;
t130 = t30 * t97;
t129 = t49 * t97;
t127 = t95 * t97;
t82 = qJ(5) + t133;
t59 = (-pkin(10) - t82) * t95;
t89 = t97 * pkin(10);
t60 = t82 * t97 + t89;
t38 = t103 * t59 - t60 * t99;
t39 = t103 * t60 + t59 * t99;
t124 = t144 * t39 - t38 * t65;
t73 = (-pkin(10) - qJ(5)) * t95;
t74 = qJ(5) * t97 + t89;
t45 = t103 * t73 - t74 * t99;
t46 = t103 * t74 + t73 * t99;
t123 = t144 * t46 - t45 * t65;
t122 = t72 + t83;
t121 = t119 * t82;
t93 = t101 ^ 2;
t94 = t105 ^ 2;
t118 = t93 + t94;
t106 = cos(qJ(2));
t115 = t96 * t106;
t114 = -0.2e1 * t70 * t68;
t112 = -pkin(4) * t70 - qJ(5) * t68;
t6 = -t13 * t95 + t14 * t97;
t32 = t100 * t57 + t104 * t58;
t19 = -t97 * t115 - t32 * t95;
t20 = -t95 * t115 + t32 * t97;
t7 = -t19 * t95 + t20 * t97;
t111 = -t68 * t82 + t70 * t85;
t110 = -t101 * t57 + t105 * t58;
t91 = t96 ^ 2;
t80 = t91 * t106 ^ 2;
t78 = 0.2e1 * t127;
t67 = t70 ^ 2;
t62 = t65 ^ 2;
t61 = t144 ^ 2;
t56 = t97 * t68;
t55 = t95 * t68;
t52 = t95 * t125;
t48 = t65 * t68;
t47 = t144 * t68;
t44 = 0.2e1 * t65 * t144;
t43 = t49 * t95;
t41 = (-t90 + t92) * t70;
t35 = t144 * t70;
t33 = t65 * t70;
t26 = pkin(5) * t128 + t49;
t25 = t30 * t95;
t22 = t35 * t65;
t21 = t33 * t144;
t18 = t30 * t65;
t17 = t30 * t144;
t16 = t26 * t65;
t15 = t26 * t144;
t10 = t144 * t35 - t33 * t65;
t9 = t103 * t20 + t19 * t99;
t8 = t103 * t19 - t20 * t99;
t1 = t144 * t9 - t65 * t8;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102 ^ 2 * t91 + t98 ^ 2 + t80, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 ^ 2 + t58 ^ 2 + t80, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 ^ 2 + t29 + t80, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19 ^ 2 + t20 ^ 2 + t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 ^ 2 + t9 ^ 2 + t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t115, -t116, 0, 0, 0, 0, 0, 0, 0, 0, t105 * t115, -t101 * t115, t110, pkin(2) * t115 + t110 * pkin(8), 0, 0, 0, 0, 0, 0, -t68 * t115, -t70 * t115, t30 * t70 - t32 * t68, -t86 * t115 + t32 * t51 + t131, 0, 0, 0, 0, 0, 0, t30 * t128 + t19 * t68, t30 * t125 - t20 * t68 (-t19 * t97 - t20 * t95) * t70, t13 * t19 + t14 * t20 + t131, 0, 0, 0, 0, 0, 0, t30 * t33 + t68 * t8, t30 * t35 - t68 * t9, -t33 * t9 - t35 * t8, t26 * t30 + t4 * t8 + t5 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t93, t101 * t138, 0, t94, 0, 0, pkin(2) * t138, -0.2e1 * pkin(2) * t101, 0.2e1 * t118 * pkin(8), t118 * pkin(8) ^ 2 + pkin(2) ^ 2, t67, t114, 0, t66, 0, 0, t68 * t139, t70 * t139, 0.2e1 * t49 * t70 - 0.2e1 * t51 * t68, t51 ^ 2 + t86 ^ 2 + t143, t92 * t67, -0.2e1 * t67 * t127, t125 * t142, t90 * t67, t95 * t114, t66, 0.2e1 * t49 * t128 + 0.2e1 * t13 * t68, 0.2e1 * t49 * t125 - 0.2e1 * t14 * t68, 0.2e1 * (-t13 * t97 - t14 * t95) * t70, t13 ^ 2 + t14 ^ 2 + t143, t35 ^ 2, -0.2e1 * t35 * t33, t35 * t142, t33 ^ 2, -t33 * t142, t66, 0.2e1 * t26 * t33 + 0.2e1 * t4 * t68, 0.2e1 * t26 * t35 - 0.2e1 * t5 * t68, -0.2e1 * t33 * t5 - 0.2e1 * t35 * t4, t26 ^ 2 + t4 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, -t58, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t32, 0 (t100 * t32 - t104 * t30) * pkin(3), 0, 0, 0, 0, 0, 0, -t130, t25, t7, t30 * t85 + t7 * t82, 0, 0, 0, 0, 0, 0, -t17, t18, t1, t30 * t72 + t38 * t8 + t39 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, 0, t105, 0, -t101 * pkin(8), -t105 * pkin(8), 0, 0, 0, 0, t70, 0, -t68, 0, -t49, -t51 (-t100 * t68 - t104 * t70) * pkin(3) (t100 * t51 - t104 * t49) * pkin(3), t52, t41, t55, -t52, t56, 0, t111 * t95 - t129, t111 * t97 + t43, t6, t49 * t85 + t6 * t82, t22, t10, t48, -t21, t47, 0, t33 * t72 + t38 * t68 - t15, t35 * t72 - t39 * t68 + t16, -t33 * t39 - t35 * t38 + t136, t26 * t72 + t38 * t4 + t39 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t132, -0.2e1 * t133, 0 (t100 ^ 2 + t104 ^ 2) * pkin(3) ^ 2, t90, t78, 0, t92, 0, 0, -0.2e1 * t85 * t97, 0.2e1 * t85 * t95, 0.2e1 * t121, t119 * t82 ^ 2 + t85 ^ 2, t62, t44, 0, t61, 0, 0, -t144 * t141, t65 * t141, 0.2e1 * t124, t38 ^ 2 + t39 ^ 2 + t72 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t32, 0, 0, 0, 0, 0, 0, 0, 0, -t130, t25, t7, -t30 * pkin(4) + t7 * qJ(5), 0, 0, 0, 0, 0, 0, -t17, t18, t1, t30 * t83 + t45 * t8 + t46 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, -t68, 0, -t49, -t51, 0, 0, t52, t41, t55, -t52, t56, 0, t112 * t95 - t129, t112 * t97 + t43, t6, -t49 * pkin(4) + t6 * qJ(5), t22, t10, t48, -t21, t47, 0, t33 * t83 + t45 * t68 - t15, t35 * t83 - t46 * t68 + t16, -t33 * t46 - t35 * t45 + t136, t26 * t83 + t4 * t45 + t46 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t132, -t133, 0, 0, t90, t78, 0, t92, 0, 0, t134 * t97, -t134 * t95, t120 + t121, -t85 * pkin(4) + t120 * t82, t62, t44, 0, t61, 0, 0, -t122 * t144, t122 * t65, t123 + t124, t38 * t45 + t39 * t46 + t72 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t90, t78, 0, t92, 0, 0, 0.2e1 * pkin(4) * t97, -0.2e1 * pkin(4) * t95, 0.2e1 * t120, t119 * qJ(5) ^ 2 + pkin(4) ^ 2, t62, t44, 0, t61, 0, 0, -t144 * t140, t65 * t140, 0.2e1 * t123, t45 ^ 2 + t46 ^ 2 + t83 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t125, 0, t49, 0, 0, 0, 0, 0, 0, t33, t35, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, t95, 0, t85, 0, 0, 0, 0, 0, 0, -t144, t65, 0, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, t95, 0, -pkin(4), 0, 0, 0, 0, 0, 0, -t144, t65, 0, t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, -t33, t68, t4, -t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, t144, 0, t38, -t39, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, t144, 0, t45, -t46, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t2;
