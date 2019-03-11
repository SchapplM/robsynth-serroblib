% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRPRP1
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
% MMD_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 16:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPRP1_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_inertiaDJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:33:01
% EndTime: 2019-03-09 16:33:06
% DurationCPUTime: 1.43s
% Computational Cost: add. (4043->216), mult. (8840->372), div. (0->0), fcn. (8568->8), ass. (0->133)
t151 = pkin(7) + pkin(8);
t93 = sin(qJ(2));
t73 = t151 * t93;
t96 = cos(qJ(2));
t74 = t151 * t96;
t92 = sin(qJ(3));
t95 = cos(qJ(3));
t106 = -t95 * t73 - t92 * t74;
t68 = t92 * t96 + t95 * t93;
t100 = -t68 * qJ(4) + t106;
t104 = t92 * t93 - t95 * t96;
t105 = t92 * t73 - t95 * t74;
t40 = -t104 * qJ(4) - t105;
t89 = sin(pkin(10));
t90 = cos(pkin(10));
t31 = t89 * t100 + t90 * t40;
t94 = cos(qJ(5));
t28 = t94 * t31;
t44 = t90 * t104 + t89 * t68;
t45 = -t89 * t104 + t90 * t68;
t81 = -t96 * pkin(2) - pkin(1);
t98 = t104 * pkin(3) + t81;
t29 = t44 * pkin(4) - t45 * pkin(9) + t98;
t91 = sin(qJ(5));
t140 = t91 * t29 + t28;
t87 = t91 ^ 2;
t88 = t94 ^ 2;
t138 = t87 - t88;
t116 = t138 * qJD(5);
t153 = qJD(2) + qJD(3);
t121 = qJD(2) * t151;
t69 = t93 * t121;
t70 = t96 * t121;
t34 = -t106 * qJD(3) + t95 * t69 + t92 * t70;
t46 = t153 * t104;
t47 = t153 * t68;
t33 = -t90 * t46 - t89 * t47;
t103 = -qJ(6) * t33 - qJD(6) * t45;
t22 = -t47 * qJ(4) - t104 * qJD(4) - t34;
t35 = t105 * qJD(3) + t92 * t69 - t95 * t70;
t97 = t46 * qJ(4) - t68 * qJD(4) + t35;
t12 = t90 * t22 + t89 * t97;
t32 = -t89 * t46 + t90 * t47;
t133 = qJD(2) * t93;
t83 = pkin(2) * t133;
t41 = t47 * pkin(3) + t83;
t16 = t32 * pkin(4) - t33 * pkin(9) + t41;
t119 = -t91 * t12 + t94 * t16;
t136 = qJ(6) * t45;
t1 = t32 * pkin(5) + t103 * t94 + (-t28 + (-t29 + t136) * t91) * qJD(5) + t119;
t85 = qJD(5) * t94;
t125 = t45 * t85;
t129 = t94 * t12 + t91 * t16 + t29 * t85;
t3 = -qJ(6) * t125 + (-qJD(5) * t31 + t103) * t91 + t129;
t118 = t94 * t29 - t91 * t31;
t86 = t94 * qJ(6);
t7 = t44 * pkin(5) - t45 * t86 + t118;
t8 = -t91 * t136 + t140;
t152 = -t1 * t94 - t3 * t91 + (t7 * t91 - t8 * t94) * qJD(5);
t150 = t94 * pkin(5);
t149 = t45 * t91;
t148 = t45 * t94;
t147 = t89 * t92;
t146 = t90 * t92;
t145 = t91 * t32;
t144 = t94 * t32;
t143 = t94 * t33;
t137 = pkin(2) * qJD(3);
t61 = (t90 * t95 - t147) * t137;
t142 = t94 * t61;
t11 = t89 * t22 - t90 * t97;
t30 = -t90 * t100 + t89 * t40;
t141 = t11 * t91 + t30 * t85;
t80 = t95 * pkin(2) + pkin(3);
t62 = -pkin(2) * t147 + t90 * t80;
t58 = -pkin(4) - t62;
t60 = (t89 * t95 + t146) * t137;
t139 = t58 * t85 + t60 * t91;
t63 = pkin(2) * t146 + t89 * t80;
t59 = pkin(9) + t63;
t135 = -qJ(6) - t59;
t78 = t89 * pkin(3) + pkin(9);
t134 = -qJ(6) - t78;
t132 = qJD(2) * t96;
t131 = qJD(5) * t91;
t130 = -0.2e1 * pkin(1) * qJD(2);
t128 = t92 * t137;
t127 = t95 * t137;
t82 = pkin(5) * t131;
t126 = pkin(5) * t85;
t79 = -t90 * pkin(3) - pkin(4);
t124 = t79 * t131;
t123 = t79 * t85;
t122 = t91 * t85;
t120 = -0.4e1 * t91 * t148;
t117 = t58 * t131 - t60 * t94;
t115 = qJD(5) * t135;
t114 = qJD(5) * t134;
t111 = -t32 * t78 + t33 * t79;
t110 = t44 * t59 - t45 * t58;
t109 = -t61 * t44 + t60 * t45;
t108 = t44 * t78 - t45 * t79;
t20 = t44 * t85 + t145;
t102 = t91 * t33 + t125;
t101 = t45 * t131 - t143;
t99 = -t32 * t59 + t33 * t58 + t109;
t84 = t94 * qJD(6);
t75 = 0.2e1 * t122;
t71 = t79 - t150;
t67 = -0.2e1 * t116;
t66 = t94 * t78 + t86;
t65 = t134 * t91;
t56 = -t91 * qJD(6) + t94 * t114;
t55 = t91 * t114 + t84;
t54 = t58 - t150;
t51 = t82 + t60;
t50 = t55 * t94;
t49 = t94 * t59 + t86;
t48 = t135 * t91;
t42 = t45 ^ 2;
t38 = (-qJD(6) - t61) * t91 + t94 * t115;
t37 = t91 * t115 + t142 + t84;
t36 = t37 * t94;
t24 = t30 * t131;
t19 = -t44 * t131 + t144;
t18 = pkin(5) * t149 + t30;
t17 = -t45 * t116 + t91 * t143;
t13 = qJD(5) * t120 - t138 * t33;
t6 = t102 * pkin(5) + t11;
t5 = -t140 * qJD(5) + t119;
t4 = t31 * t131 - t129;
t2 = t3 * t94;
t9 = [0, 0, 0, 0.2e1 * t93 * t132, 0.2e1 * (-t93 ^ 2 + t96 ^ 2) * qJD(2), 0, 0, 0, t93 * t130, t96 * t130, -0.2e1 * t68 * t46, 0.2e1 * t46 * t104 - 0.2e1 * t68 * t47, 0, 0, 0, 0.2e1 * t104 * t83 + 0.2e1 * t81 * t47, -0.2e1 * t81 * t46 + 0.2e1 * t68 * t83, 0.2e1 * t11 * t45 - 0.2e1 * t12 * t44 + 0.2e1 * t30 * t33 - 0.2e1 * t31 * t32, 0.2e1 * t30 * t11 + 0.2e1 * t31 * t12 + 0.2e1 * t41 * t98, 0.2e1 * t88 * t45 * t33 - 0.2e1 * t122 * t42, 0.2e1 * t42 * t116 + t33 * t120, -0.2e1 * t101 * t44 + 0.2e1 * t45 * t144, -0.2e1 * t102 * t44 - 0.2e1 * t45 * t145, 0.2e1 * t44 * t32, 0.2e1 * t102 * t30 + 0.2e1 * t11 * t149 + 0.2e1 * t118 * t32 + 0.2e1 * t5 * t44, -0.2e1 * t101 * t30 + 0.2e1 * t11 * t148 - 0.2e1 * t140 * t32 + 0.2e1 * t4 * t44, 0.2e1 * (-t7 * t94 - t8 * t91) * t33 + 0.2e1 * t152 * t45, 0.2e1 * t7 * t1 + 0.2e1 * t18 * t6 + 0.2e1 * t8 * t3; 0, 0, 0, 0, 0, t132, -t133, 0, -pkin(7) * t132, pkin(7) * t133, 0, 0, -t46, -t47, 0, t35, t34, -t63 * t32 - t62 * t33 + t109, -t11 * t62 + t12 * t63 + t30 * t60 + t31 * t61, t17, t13, t20, t19, 0, t24 + (-qJD(5) * t110 - t11) * t94 + t99 * t91, t110 * t131 + t94 * t99 + t141, t2 + (-t33 * t48 - t38 * t45 + (-t45 * t49 - t7) * qJD(5)) * t94 + (-t33 * t49 - t37 * t45 - t1 + (t45 * t48 - t8) * qJD(5)) * t91, t1 * t48 + t18 * t51 + t3 * t49 + t8 * t37 + t7 * t38 + t6 * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t128, -0.2e1 * t127, 0, -0.2e1 * t62 * t60 + 0.2e1 * t63 * t61, t75, t67, 0, 0, 0, 0.2e1 * t117, 0.2e1 * t139, -0.2e1 * t38 * t91 + 0.2e1 * t36 + 0.2e1 * (-t48 * t94 - t49 * t91) * qJD(5), 0.2e1 * t49 * t37 + 0.2e1 * t48 * t38 + 0.2e1 * t54 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t46, -t47, 0, t35, t34 (-t32 * t89 - t33 * t90) * pkin(3) (-t11 * t90 + t12 * t89) * pkin(3), t17, t13, t20, t19, 0, t24 + t111 * t91 + (-qJD(5) * t108 - t11) * t94, t108 * t131 + t111 * t94 + t141, t2 + (-t33 * t65 - t45 * t56 + (-t45 * t66 - t7) * qJD(5)) * t94 + (-t33 * t66 - t45 * t55 - t1 + (t45 * t65 - t8) * qJD(5)) * t91, t1 * t65 + t18 * t82 + t3 * t66 + t8 * t55 + t7 * t56 + t6 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, -t127, 0 (-t60 * t90 + t61 * t89) * pkin(3), t75, t67, 0, 0, 0, t117 + t124, t123 + t139, t36 + t50 + (-t38 - t56) * t91 + ((-t48 - t65) * t94 + (-t49 - t66) * t91) * qJD(5), t37 * t66 + t38 * t65 + t48 * t56 + t49 * t55 + t51 * t71 + t54 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t67, 0, 0, 0, 0.2e1 * t124, 0.2e1 * t123, -0.2e1 * t56 * t91 + 0.2e1 * t50 + 0.2e1 * (-t65 * t94 - t66 * t91) * qJD(5), 0.2e1 * t66 * t55 + 0.2e1 * t65 * t56 + 0.2e1 * t71 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, t19, -t20 (-t87 - t88) * t33, -t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37 * t91 + t38 * t94 + (-t48 * t91 + t49 * t94) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 * t91 + t56 * t94 + (-t65 * t91 + t66 * t94) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t101, -t102, t32, t5, t4, t101 * pkin(5), t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, -t131, 0, -t59 * t85 - t91 * t61, t131 * t59 - t142, -t126, t38 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, -t131, 0, -t78 * t85, t78 * t131, -t126, t56 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131, -t85, 0, -t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
