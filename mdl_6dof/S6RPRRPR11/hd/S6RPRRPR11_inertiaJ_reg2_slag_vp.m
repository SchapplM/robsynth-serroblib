% Calculate inertial parameters regressor of joint inertia matrix for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRPR11_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_inertiaJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_inertiaJ_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 00:13:36
% EndTime: 2019-05-06 00:13:47
% DurationCPUTime: 3.43s
% Computational Cost: add. (6501->259), mult. (17406->538), div. (0->0), fcn. (20320->14), ass. (0->132)
t105 = sin(qJ(3));
t107 = cos(qJ(3));
t101 = cos(pkin(7));
t102 = cos(pkin(6));
t96 = sin(pkin(12));
t98 = sin(pkin(6));
t140 = t96 * t98;
t100 = cos(pkin(12));
t146 = pkin(1) * t102;
t81 = t100 * t146;
t45 = t102 * pkin(2) + t81 + (-pkin(9) * t101 - qJ(2)) * t140;
t97 = sin(pkin(7));
t52 = (-pkin(9) * t96 * t97 - pkin(2) * t100 - pkin(1)) * t98;
t113 = t101 * t45 + t52 * t97;
t126 = t98 * t100;
t118 = t101 * t126;
t110 = t102 * t97 + t118;
t62 = qJ(2) * t126 + t96 * t146;
t40 = t110 * pkin(9) + t62;
t23 = -t105 * t40 + t113 * t107;
t104 = sin(qJ(4));
t106 = cos(qJ(4));
t44 = t110 * t105 + t107 * t140;
t59 = t102 * t101 - t97 * t126;
t35 = t44 * t104 - t59 * t106;
t34 = t35 ^ 2;
t127 = t97 * t107;
t42 = -t102 * t127 + t105 * t140 - t107 * t118;
t158 = t42 ^ 2;
t128 = t97 * t105;
t63 = -t106 * t101 + t104 * t128;
t60 = t63 ^ 2;
t157 = -0.2e1 * t35;
t156 = 0.2e1 * t35;
t155 = -0.2e1 * t42;
t99 = cos(pkin(13));
t86 = -t99 * pkin(5) - pkin(4);
t154 = 0.2e1 * t86;
t153 = 0.2e1 * t98;
t152 = -0.2e1 * t106;
t95 = sin(pkin(13));
t151 = pkin(4) * t95;
t150 = pkin(4) * t99;
t149 = pkin(10) * t95;
t93 = t104 ^ 2;
t148 = t93 * pkin(10);
t21 = -t59 * pkin(3) - t23;
t37 = t59 * t104 + t44 * t106;
t15 = t35 * pkin(4) - t37 * qJ(5) + t21;
t33 = t101 * t52 - t97 * t45;
t19 = t42 * pkin(3) - t44 * pkin(10) + t33;
t24 = t113 * t105 + t107 * t40;
t22 = t59 * pkin(10) + t24;
t12 = t104 * t19 + t106 * t22;
t9 = t42 * qJ(5) + t12;
t6 = t95 * t15 + t99 * t9;
t147 = cos(qJ(6));
t11 = -t104 * t22 + t106 * t19;
t10 = -t42 * pkin(4) - t11;
t145 = t10 * t95;
t144 = t10 * t99;
t88 = t104 * pkin(10);
t26 = t37 * t95 - t42 * t99;
t143 = t26 * t99;
t28 = t37 * t99 + t42 * t95;
t142 = t28 * t95;
t141 = t95 * t99;
t139 = pkin(11) + qJ(5);
t124 = t99 * t106;
t73 = -t106 * pkin(4) - t104 * qJ(5) - pkin(3);
t54 = pkin(10) * t124 + t95 * t73;
t89 = t95 ^ 2;
t92 = t99 ^ 2;
t138 = t89 + t92;
t137 = qJ(5) * t35;
t91 = t98 ^ 2;
t136 = t100 * t91;
t135 = t104 * t42;
t134 = t106 * t42;
t133 = t35 * t106;
t132 = t37 * t104;
t131 = t63 * t104;
t130 = t95 * t104;
t129 = t95 * t106;
t125 = t99 * t104;
t123 = qJ(5) * t106;
t122 = t104 * t106;
t121 = t102 * t153;
t120 = 0.2e1 * t122;
t119 = t95 * t125;
t117 = t147 * t99;
t5 = t99 * t15 - t95 * t9;
t116 = -t5 * t95 + t6 * t99;
t65 = t104 * t101 + t106 * t128;
t46 = -t99 * t127 - t95 * t65;
t47 = -t95 * t127 + t99 * t65;
t115 = -t46 * t95 + t47 * t99;
t67 = t99 * t73;
t53 = -pkin(10) * t129 + t67;
t114 = -t53 * t95 + t54 * t99;
t112 = -t11 * t104 + t12 * t106;
t111 = t65 * t106 + t131;
t103 = sin(qJ(6));
t70 = t103 * t99 + t147 * t95;
t109 = pkin(10) ^ 2;
t94 = t106 ^ 2;
t90 = t97 ^ 2;
t87 = t93 * t109;
t83 = t90 * t107 ^ 2;
t75 = t139 * t99;
t74 = t139 * t95;
t72 = pkin(5) * t130 + t88;
t68 = t103 * t95 - t117;
t61 = -qJ(2) * t140 + t81;
t58 = -t103 * t130 + t104 * t117;
t56 = t70 * t104;
t50 = -t103 * t74 + t147 * t75;
t49 = -t103 * t75 - t147 * t74;
t48 = -pkin(11) * t130 + t54;
t41 = -pkin(11) * t125 + t67 + (-pkin(5) - t149) * t106;
t32 = t103 * t46 + t147 * t47;
t31 = -t103 * t47 + t147 * t46;
t30 = t103 * t41 + t147 * t48;
t29 = -t103 * t48 + t147 * t41;
t18 = -t103 * t26 + t147 * t28;
t16 = t103 * t28 + t147 * t26;
t7 = t26 * pkin(5) + t10;
t4 = -t26 * pkin(11) + t6;
t3 = t35 * pkin(5) - t28 * pkin(11) + t5;
t2 = t103 * t3 + t147 * t4;
t1 = -t103 * t4 + t147 * t3;
t8 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t91 * t96 ^ 2, 0.2e1 * t96 * t136, t96 * t121, t91 * t100 ^ 2, t100 * t121, t102 ^ 2, 0.2e1 * pkin(1) * t136 + 0.2e1 * t61 * t102, -0.2e1 * t91 * pkin(1) * t96 - 0.2e1 * t62 * t102 (t100 * t62 - t61 * t96) * t153, t91 * pkin(1) ^ 2 + t61 ^ 2 + t62 ^ 2, t44 ^ 2, t44 * t155, 0.2e1 * t44 * t59, t158, t59 * t155, t59 ^ 2, 0.2e1 * t23 * t59 + 0.2e1 * t33 * t42, -0.2e1 * t24 * t59 + 0.2e1 * t33 * t44, -0.2e1 * t23 * t44 - 0.2e1 * t24 * t42, t23 ^ 2 + t24 ^ 2 + t33 ^ 2, t37 ^ 2, t37 * t157, 0.2e1 * t37 * t42, t34, t35 * t155, t158, 0.2e1 * t11 * t42 + 0.2e1 * t21 * t35, -0.2e1 * t12 * t42 + 0.2e1 * t21 * t37, -0.2e1 * t11 * t37 - 0.2e1 * t12 * t35, t11 ^ 2 + t12 ^ 2 + t21 ^ 2, t28 ^ 2, -0.2e1 * t28 * t26, t28 * t156, t26 ^ 2, t26 * t157, t34, 0.2e1 * t10 * t26 + 0.2e1 * t5 * t35, 0.2e1 * t10 * t28 - 0.2e1 * t6 * t35, -0.2e1 * t6 * t26 - 0.2e1 * t5 * t28, t10 ^ 2 + t5 ^ 2 + t6 ^ 2, t18 ^ 2, -0.2e1 * t18 * t16, t18 * t156, t16 ^ 2, t16 * t157, t34, 0.2e1 * t1 * t35 + 0.2e1 * t7 * t16, 0.2e1 * t7 * t18 - 0.2e1 * t2 * t35, -0.2e1 * t1 * t18 - 0.2e1 * t2 * t16, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126, t140, 0, -t98 * pkin(1), 0, 0, 0, 0, 0, 0, t101 * t42 + t59 * t127, t101 * t44 - t59 * t128 (-t105 * t42 - t107 * t44) * t97, t33 * t101 + (t105 * t24 + t107 * t23) * t97, 0, 0, 0, 0, 0, 0, -t35 * t127 - t63 * t42, -t37 * t127 - t65 * t42, -t65 * t35 + t63 * t37, -t11 * t63 + t12 * t65 - t21 * t127, 0, 0, 0, 0, 0, 0, t63 * t26 + t46 * t35, t63 * t28 - t47 * t35, -t47 * t26 - t46 * t28, t10 * t63 + t5 * t46 + t6 * t47, 0, 0, 0, 0, 0, 0, t63 * t16 + t31 * t35, t63 * t18 - t32 * t35, -t32 * t16 - t31 * t18, t1 * t31 + t2 * t32 + t7 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90 * t105 ^ 2 + t101 ^ 2 + t83, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65 ^ 2 + t60 + t83, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46 ^ 2 + t47 ^ 2 + t60, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 ^ 2 + t32 ^ 2 + t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, -t42, t59, t23, -t24, 0, 0, t132, -t104 * t35 + t37 * t106, t135, -t133, t134, 0, -pkin(3) * t35 - pkin(10) * t135 - t21 * t106, -pkin(3) * t37 - pkin(10) * t134 + t21 * t104 (t132 - t133) * pkin(10) + t112, -t21 * pkin(3) + t112 * pkin(10), t28 * t125 (-t142 - t143) * t104, -t28 * t106 + t125 * t35, t26 * t130, t26 * t106 - t130 * t35, -t133, -t5 * t106 + t53 * t35 + (pkin(10) * t26 + t145) * t104, t6 * t106 - t54 * t35 + (pkin(10) * t28 + t144) * t104, -t54 * t26 - t53 * t28 + (-t5 * t99 - t6 * t95) * t104, t10 * t88 + t5 * t53 + t6 * t54, t18 * t58, -t58 * t16 - t18 * t56, -t18 * t106 + t58 * t35, t16 * t56, t16 * t106 - t56 * t35, -t133, -t1 * t106 + t72 * t16 + t29 * t35 + t7 * t56, t2 * t106 + t72 * t18 - t30 * t35 + t7 * t58, -t1 * t58 - t30 * t16 - t29 * t18 - t2 * t56, t1 * t29 + t2 * t30 + t7 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, -t128, 0, 0, 0, 0, 0, 0, 0, 0, t106 * t127, -t104 * t127, t111, pkin(3) * t127 + t111 * pkin(10), 0, 0, 0, 0, 0, 0, -t46 * t106 + t130 * t63, t47 * t106 + t125 * t63 (-t46 * t99 - t47 * t95) * t104, pkin(10) * t131 + t46 * t53 + t47 * t54, 0, 0, 0, 0, 0, 0, -t31 * t106 + t63 * t56, t32 * t106 + t63 * t58, -t31 * t58 - t32 * t56, t31 * t29 + t32 * t30 + t63 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t93, t120, 0, t94, 0, 0, 0.2e1 * pkin(3) * t106, -0.2e1 * pkin(3) * t104, 0.2e1 * (t93 + t94) * pkin(10), pkin(3) ^ 2 + t94 * t109 + t87, t92 * t93, -0.2e1 * t93 * t141, -0.2e1 * t99 * t122, t89 * t93, t95 * t120, t94, -0.2e1 * t53 * t106 + 0.2e1 * t95 * t148, 0.2e1 * t54 * t106 + 0.2e1 * t99 * t148, 0.2e1 * (-t53 * t99 - t54 * t95) * t104, t53 ^ 2 + t54 ^ 2 + t87, t58 ^ 2, -0.2e1 * t58 * t56, t58 * t152, t56 ^ 2, -t56 * t152, t94, -0.2e1 * t29 * t106 + 0.2e1 * t72 * t56, 0.2e1 * t30 * t106 + 0.2e1 * t72 * t58, -0.2e1 * t29 * t58 - 0.2e1 * t30 * t56, t29 ^ 2 + t30 ^ 2 + t72 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, 0, -t35, t42, t11, -t12, 0, 0, t142, -t95 * t26 + t28 * t99, t95 * t35, -t143, t99 * t35, 0, -pkin(4) * t26 - t95 * t137 - t144, -pkin(4) * t28 - t99 * t137 + t145 (t142 - t143) * qJ(5) + t116, -t10 * pkin(4) + qJ(5) * t116, t18 * t70, -t70 * t16 - t18 * t68, t70 * t35, t16 * t68, -t68 * t35, 0, t86 * t16 + t49 * t35 + t7 * t68, t86 * t18 - t50 * t35 + t7 * t70, -t1 * t70 - t50 * t16 - t49 * t18 - t2 * t68, t1 * t49 + t2 * t50 + t7 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t65, 0, 0, 0, 0, 0, 0, 0, 0, -t63 * t99, t63 * t95, t115, -t63 * pkin(4) + qJ(5) * t115, 0, 0, 0, 0, 0, 0, t63 * t68, t63 * t70, -t31 * t70 - t32 * t68, t31 * t49 + t32 * t50 + t63 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, 0, t106, 0, -t88, -t106 * pkin(10), 0, 0, t119 (-t89 + t92) * t104, -t129, -t119, -t124, 0, t95 * t123 + (-pkin(10) * t99 - t151) * t104, t99 * t123 + (t149 - t150) * t104, t114, -pkin(4) * t88 + qJ(5) * t114, t58 * t70, -t70 * t56 - t58 * t68, -t70 * t106, t56 * t68, t68 * t106, 0, -t49 * t106 + t86 * t56 + t72 * t68, t50 * t106 + t86 * t58 + t72 * t70, -t29 * t70 - t30 * t68 - t49 * t58 - t50 * t56, t29 * t49 + t30 * t50 + t72 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t89, 0.2e1 * t141, 0, t92, 0, 0, 0.2e1 * t150, -0.2e1 * t151, 0.2e1 * t138 * qJ(5), t138 * qJ(5) ^ 2 + pkin(4) ^ 2, t70 ^ 2, -0.2e1 * t70 * t68, 0, t68 ^ 2, 0, 0, t68 * t154, t70 * t154, -0.2e1 * t49 * t70 - 0.2e1 * t50 * t68, t49 ^ 2 + t50 ^ 2 + t86 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t28, 0, t10, 0, 0, 0, 0, 0, 0, t16, t18, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t125, 0, t88, 0, 0, 0, 0, 0, 0, t56, t58, 0, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, t95, 0, -pkin(4), 0, 0, 0, 0, 0, 0, t68, t70, 0, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t16, t35, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t32, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, -t56, -t106, t29, -t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, -t68, 0, t49, -t50, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg  = t8;
