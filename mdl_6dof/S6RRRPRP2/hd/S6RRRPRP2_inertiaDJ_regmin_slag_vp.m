% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPRP2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:37:18
% EndTime: 2019-03-09 16:37:23
% DurationCPUTime: 1.78s
% Computational Cost: add. (5122->205), mult. (11190->365), div. (0->0), fcn. (10917->8), ass. (0->140)
t86 = sin(qJ(3));
t87 = sin(qJ(2));
t89 = cos(qJ(3));
t90 = cos(qJ(2));
t104 = t86 * t87 - t89 * t90;
t140 = cos(pkin(10));
t66 = t86 * t90 + t89 * t87;
t84 = sin(pkin(10));
t46 = t140 * t104 + t84 * t66;
t47 = -t84 * t104 + t140 * t66;
t78 = -t90 * pkin(2) - pkin(1);
t93 = t104 * pkin(3) + t78;
t30 = t46 * pkin(4) - t47 * pkin(9) + t93;
t162 = pkin(7) + pkin(8);
t70 = t162 * t87;
t71 = t162 * t90;
t105 = t86 * t70 - t89 * t71;
t41 = -t104 * qJ(4) - t105;
t106 = -t89 * t70 - t86 * t71;
t98 = -t66 * qJ(4) + t106;
t32 = t140 * t41 + t84 * t98;
t85 = sin(qJ(5));
t88 = cos(qJ(5));
t166 = t85 * t30 + t88 * t32;
t82 = t85 ^ 2;
t83 = t88 ^ 2;
t145 = t82 - t83;
t119 = t145 * qJD(5);
t100 = t66 * qJD(2);
t124 = qJD(2) * t162;
t67 = t87 * t124;
t68 = t90 * t124;
t107 = t89 * t67 + t86 * t68;
t24 = -qJ(4) * t100 + t98 * qJD(3) - t104 * qJD(4) - t107;
t36 = t105 * qJD(3) + t86 * t67 - t89 * t68;
t49 = (-qJD(2) - qJD(3)) * t104;
t91 = -t49 * qJ(4) - t66 * qJD(4) + t36;
t12 = t140 * t24 + t84 * t91;
t92 = -t66 * qJD(3) - t100;
t33 = -t140 * t92 + t84 * t49;
t34 = t140 * t49 + t84 * t92;
t139 = qJD(2) * t87;
t80 = pkin(2) * t139;
t42 = -t92 * pkin(3) + t80;
t16 = t33 * pkin(4) - t34 * pkin(9) + t42;
t5 = -qJD(5) * t166 - t85 * t12 + t88 * t16;
t115 = t88 * pkin(5) + t85 * qJ(6);
t165 = t115 * qJD(5) - t88 * qJD(6);
t113 = t88 * t30 - t85 * t32;
t10 = -t46 * pkin(5) - t113;
t9 = t46 * qJ(6) + t166;
t116 = t10 * t85 + t88 * t9;
t135 = t46 * qJD(6);
t142 = t33 * qJ(6);
t136 = qJD(5) * t85;
t81 = qJD(5) * t88;
t4 = -t88 * t12 + t32 * t136 - t85 * t16 - t30 * t81;
t2 = t135 - t4 + t142;
t161 = t33 * pkin(5);
t3 = -t161 - t5;
t164 = qJD(5) * t116 + t2 * t85 - t3 * t88;
t163 = 0.2e1 * qJD(6);
t11 = -t140 * t91 + t84 * t24;
t31 = -t140 * t98 + t84 * t41;
t160 = t11 * t85 + t31 * t81;
t120 = t140 * t86;
t77 = t89 * pkin(2) + pkin(3);
t60 = pkin(2) * t120 + t84 * t77;
t56 = pkin(9) + t60;
t159 = t33 * t56;
t75 = t84 * pkin(3) + pkin(9);
t158 = t33 * t75;
t157 = t46 * t56;
t156 = t46 * t75;
t155 = t47 * t85;
t154 = t47 * t88;
t143 = pkin(2) * qJD(3);
t152 = t84 * t86;
t58 = (t140 * t89 - t152) * t143;
t153 = t58 * t46;
t151 = t85 * t33;
t150 = t88 * t33;
t149 = t88 * t34;
t57 = (t84 * t89 + t120) * t143;
t134 = t85 * qJD(6);
t62 = -pkin(5) * t136 + qJ(6) * t81 + t134;
t43 = t57 - t62;
t147 = -t43 + t62;
t59 = -pkin(2) * t152 + t140 * t77;
t55 = -pkin(4) - t59;
t146 = t55 * t81 + t57 * t85;
t144 = t82 + t83;
t138 = qJD(2) * t90;
t114 = pkin(5) * t85 - qJ(6) * t88;
t18 = t114 * t47 + t31;
t137 = qJD(5) * t18;
t132 = -0.2e1 * pkin(1) * qJD(2);
t131 = t86 * t143;
t130 = t89 * t143;
t129 = t75 * t136;
t128 = t75 * t81;
t76 = -t140 * pkin(3) - pkin(4);
t127 = t76 * t136;
t126 = t76 * t81;
t125 = t85 * t81;
t40 = t144 * t58;
t123 = -0.4e1 * t85 * t154;
t121 = t55 * t136 - t57 * t88;
t117 = t10 * t88 - t85 * t9;
t111 = t34 * t76 - t158;
t110 = -t47 * t55 + t157;
t109 = t57 * t47 - t153;
t108 = -t47 * t76 + t156;
t22 = t46 * t81 + t151;
t20 = -t46 * t136 + t150;
t103 = t85 * t34 + t47 * t81;
t102 = t47 * t136 - t149;
t101 = t78 * t66;
t64 = -t115 + t76;
t99 = t34 * t64 - t47 * t62 - t158;
t50 = -t115 + t55;
t6 = t114 * t34 + t165 * t47 + t11;
t97 = -t6 + (t47 * t50 - t157) * qJD(5);
t96 = -t6 + (t47 * t64 - t156) * qJD(5);
t95 = t34 * t50 + t43 * t47 - t153 - t159;
t94 = t34 * t55 + t109 - t159;
t1 = t117 * qJD(5) + t2 * t88 + t3 * t85;
t72 = 0.2e1 * t125;
t65 = -0.2e1 * t119;
t54 = t64 * t136;
t48 = t50 * t136;
t44 = t47 ^ 2;
t38 = t56 * t81 + t85 * t58;
t37 = t56 * t136 - t88 * t58;
t35 = -t106 * qJD(3) + t107;
t26 = t31 * t136;
t19 = -t47 * t119 + t85 * t149;
t17 = t18 * t136;
t13 = qJD(5) * t123 - t145 * t34;
t7 = [0, 0, 0, 0.2e1 * t87 * t138, 0.2e1 * (-t87 ^ 2 + t90 ^ 2) * qJD(2), 0, 0, 0, t87 * t132, t90 * t132, 0.2e1 * t66 * t49, -0.2e1 * t49 * t104 + 0.2e1 * t66 * t92, 0, 0, 0, 0.2e1 * qJD(3) * t101 + 0.2e1 * (t87 * pkin(2) * t104 + t101) * qJD(2), 0.2e1 * t78 * t49 + 0.2e1 * t66 * t80, 0.2e1 * t11 * t47 - 0.2e1 * t12 * t46 + 0.2e1 * t31 * t34 - 0.2e1 * t32 * t33, 0.2e1 * t31 * t11 + 0.2e1 * t32 * t12 + 0.2e1 * t42 * t93, 0.2e1 * t83 * t47 * t34 - 0.2e1 * t125 * t44, 0.2e1 * t44 * t119 + t123 * t34, -0.2e1 * t102 * t46 + 0.2e1 * t150 * t47, -0.2e1 * t103 * t46 - 0.2e1 * t151 * t47, 0.2e1 * t46 * t33, 0.2e1 * t103 * t31 + 0.2e1 * t11 * t155 + 0.2e1 * t113 * t33 + 0.2e1 * t5 * t46, -0.2e1 * t102 * t31 + 0.2e1 * t11 * t154 - 0.2e1 * t166 * t33 + 0.2e1 * t4 * t46, -0.2e1 * t10 * t33 + 0.2e1 * t103 * t18 + 0.2e1 * t6 * t155 - 0.2e1 * t3 * t46, 0.2e1 * t117 * t34 - 0.2e1 * t164 * t47, 0.2e1 * t102 * t18 - 0.2e1 * t6 * t154 + 0.2e1 * t2 * t46 + 0.2e1 * t9 * t33, 0.2e1 * t10 * t3 + 0.2e1 * t18 * t6 + 0.2e1 * t9 * t2; 0, 0, 0, 0, 0, t138, -t139, 0, -pkin(7) * t138, pkin(7) * t139, 0, 0, t49, t92, 0, t36, t35, -t60 * t33 - t59 * t34 + t109, -t11 * t59 + t12 * t60 + t31 * t57 + t32 * t58, t19, t13, t22, t20, 0, t26 + (-qJD(5) * t110 - t11) * t88 + t94 * t85, t110 * t136 + t88 * t94 + t160, t85 * t95 + t88 * t97 + t17, t1, t97 * t85 + (-t95 - t137) * t88, t1 * t56 + t116 * t58 + t18 * t43 + t6 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t131, -0.2e1 * t130, 0, -0.2e1 * t59 * t57 + 0.2e1 * t60 * t58, t72, t65, 0, 0, 0, 0.2e1 * t121, 0.2e1 * t146, -0.2e1 * t43 * t88 + 0.2e1 * t48, 0.2e1 * t40, -0.2e1 * t43 * t85 - 0.2e1 * t50 * t81, 0.2e1 * t40 * t56 + 0.2e1 * t50 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t92, 0, t36, t35 (-t140 * t34 - t33 * t84) * pkin(3) (-t11 * t140 + t12 * t84) * pkin(3), t19, t13, t22, t20, 0, t26 + t111 * t85 + (-qJD(5) * t108 - t11) * t88, t108 * t136 + t111 * t88 + t160, t85 * t99 + t88 * t96 + t17, t1, t96 * t85 + (-t99 - t137) * t88, t1 * t75 - t18 * t62 + t6 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, -t130, 0 (-t140 * t57 + t58 * t84) * pkin(3), t72, t65, 0, 0, 0, t121 + t127, t126 + t146, t147 * t88 + t48 + t54, t40, t147 * t85 + (-t50 - t64) * t81, t40 * t75 + t43 * t64 - t50 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72, t65, 0, 0, 0, 0.2e1 * t127, 0.2e1 * t126, 0.2e1 * t62 * t88 + 0.2e1 * t54, 0, 0.2e1 * t62 * t85 - 0.2e1 * t64 * t81, -0.2e1 * t64 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, t20, -t22, t20, -t144 * t34, t22, t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t102, -t103, t33, t5, t4, t5 + 0.2e1 * t161, -t115 * t34 + (qJD(5) * t114 - t134) * t47, 0.2e1 * t135 - t4 + 0.2e1 * t142, -t3 * pkin(5) + t2 * qJ(6) + t9 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t136, 0, -t38, t37, -t38, -t165, -t37, -t114 * t58 - t165 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, -t136, 0, -t128, t129, -t128, -t165, -t129, -t165 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, -t81, -t136, 0, t81, t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, qJ(6) * t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t102, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, 0, t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t7;
