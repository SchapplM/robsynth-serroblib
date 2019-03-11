% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPRR8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR8_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR8_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR8_inertiaDJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:40:37
% EndTime: 2019-03-08 22:40:44
% DurationCPUTime: 2.08s
% Computational Cost: add. (1665->278), mult. (5052->516), div. (0->0), fcn. (4884->12), ass. (0->159)
t99 = cos(qJ(3));
t131 = -pkin(2) * t99 - pkin(3);
t89 = sin(pkin(7));
t95 = sin(qJ(3));
t176 = t89 * t95;
t78 = pkin(9) * t176;
t91 = cos(pkin(7));
t34 = pkin(4) * t176 + t78 + (-pkin(10) + t131) * t91;
t101 = -pkin(3) - pkin(10);
t168 = qJ(4) * t95;
t44 = (t101 * t99 - pkin(2) - t168) * t89;
t94 = sin(qJ(5));
t98 = cos(qJ(5));
t185 = t94 * t34 + t98 * t44;
t175 = t89 * t99;
t179 = pkin(2) * t91;
t81 = t95 * t179;
t184 = pkin(9) * t175 + t81;
t121 = pkin(5) * t98 + pkin(11) * t94;
t97 = cos(qJ(6));
t183 = t97 * t121;
t87 = t97 ^ 2;
t93 = sin(qJ(6));
t171 = t93 ^ 2 - t87;
t127 = t171 * qJD(6);
t160 = qJD(4) * t95;
t162 = qJD(3) * t95;
t135 = t89 * t162;
t75 = pkin(3) * t135;
t32 = t75 + (-t160 + (pkin(10) * t95 - qJ(4) * t99) * qJD(3)) * t89;
t180 = pkin(4) + pkin(9);
t45 = (t180 * t175 + t81) * qJD(3);
t10 = -qJD(5) * t185 - t94 * t32 + t98 * t45;
t182 = 0.2e1 * t89;
t181 = 0.2e1 * qJD(4);
t146 = t94 * t175;
t60 = t91 * t98 - t146;
t113 = t97 * t176 - t60 * t93;
t59 = t98 * t175 + t91 * t94;
t39 = -qJD(5) * t59 + t94 * t135;
t161 = qJD(3) * t99;
t77 = t89 * t161;
t16 = t113 * qJD(6) + t97 * t39 + t93 * t77;
t178 = t16 * t93;
t56 = t184 * qJD(3);
t177 = t56 * t91;
t96 = sin(qJ(2));
t174 = t96 * t99;
t173 = t97 * t98;
t86 = t94 ^ 2;
t88 = t98 ^ 2;
t170 = t86 - t88;
t169 = t86 + t88;
t100 = cos(qJ(2));
t167 = t100 * t91;
t166 = t100 * t95;
t165 = t101 * t94;
t164 = t101 * t98;
t90 = sin(pkin(6));
t163 = qJD(2) * t90;
t159 = qJD(5) * t113;
t42 = t93 * t176 + t60 * t97;
t158 = qJD(5) * t42;
t157 = qJD(5) * t93;
t156 = qJD(5) * t94;
t155 = qJD(5) * t97;
t154 = qJD(5) * t98;
t153 = qJD(6) * t93;
t152 = qJD(6) * t97;
t151 = qJD(6) * t98;
t150 = qJ(4) * qJD(5);
t149 = qJD(5) * t101;
t148 = qJD(6) * t101;
t147 = -0.2e1 * pkin(5) * qJD(6);
t52 = -t91 * qJ(4) - t184;
t145 = t101 * t176;
t144 = t93 * t165;
t143 = t97 * t165;
t142 = t93 * t164;
t84 = t89 ^ 2;
t141 = t84 * t161;
t140 = t59 * t157;
t139 = t59 * t155;
t138 = t94 * t155;
t137 = t93 * t151;
t136 = t97 * t151;
t134 = t96 * t163;
t133 = t93 * t152;
t132 = t94 * t154;
t130 = t98 * t149;
t129 = t93 * t148;
t128 = t100 * t163;
t126 = t170 * qJD(5);
t43 = pkin(4) * t175 - t52;
t125 = t89 * t134;
t124 = t95 * t134;
t123 = t93 * t138;
t122 = t101 * t77;
t120 = pkin(5) * t94 - pkin(11) * t98;
t119 = -pkin(3) * t99 - t168;
t19 = pkin(11) * t176 + t185;
t21 = pkin(5) * t59 - pkin(11) * t60 + t43;
t8 = t19 * t97 + t21 * t93;
t92 = cos(pkin(6));
t35 = -t92 * t175 + (-t99 * t167 + t95 * t96) * t90;
t57 = -t100 * t89 * t90 + t91 * t92;
t26 = t35 * t94 + t57 * t98;
t110 = t91 * t166 + t174;
t36 = t110 * t90 + t92 * t176;
t15 = t26 * t97 + t36 * t93;
t14 = -t26 * t93 + t36 * t97;
t117 = t34 * t98 - t44 * t94;
t115 = t35 * t98 - t57 * t94;
t114 = -t113 * t97 + t42 * t93;
t76 = t161 * t179;
t55 = pkin(9) * t135 - t76;
t18 = -pkin(5) * t176 - t117;
t6 = -pkin(5) * t77 - t10;
t112 = t18 * t152 + t6 * t93;
t111 = t18 * t153 - t6 * t97;
t68 = qJ(4) + t120;
t49 = t93 * t68 + t143;
t23 = t92 * t135 + (t110 * qJD(3) + (t91 * t174 + t166) * qJD(2)) * t90;
t11 = t26 * qJD(5) + t94 * t125 - t23 * t98;
t109 = t11 * t93 - t115 * t152;
t108 = -t11 * t97 - t115 * t153;
t40 = -qJD(5) * t146 - t98 * t135 + t91 * t154;
t107 = t59 * t152 + t40 * t93;
t106 = t59 * t153 - t40 * t97;
t9 = -t34 * t154 + t44 * t156 - t98 * t32 - t94 * t45;
t62 = -t137 - t138;
t82 = t91 * qJD(4);
t33 = -t180 * t135 + t76 + t82;
t105 = pkin(11) * t77 - t9;
t104 = pkin(5) * t40 - pkin(11) * t39 + t33;
t103 = t84 * t99 * t134 - t57 * t135 + t23 * t91;
t24 = -t91 * t124 - t90 * t96 * t162 + (t128 + (t90 * t167 + t89 * t92) * qJD(3)) * t99;
t102 = t84 * t124 - t24 * t91 + t57 * t77;
t67 = 0.2e1 * t95 * t141;
t64 = t93 * t156 - t136;
t63 = t94 * t152 + t93 * t154;
t61 = t94 * t153 - t97 * t154;
t54 = t131 * t91 + t78;
t53 = (-pkin(2) + t119) * t89;
t51 = (-t95 * t154 - t94 * t161) * t89;
t50 = (-t95 * t156 + t98 * t161) * t89;
t48 = t97 * t68 - t144;
t47 = t55 - t82;
t46 = t75 + (-qJ(4) * t161 - t160) * t89;
t28 = t97 * qJD(4) - t49 * qJD(6) + (-t142 + t183) * qJD(5);
t27 = t94 * t129 - t93 * (t121 * qJD(5) + qJD(4)) - t68 * t152 - t97 * t130;
t17 = t42 * qJD(6) + t93 * t39 - t97 * t77;
t12 = t115 * qJD(5) + t98 * t125 + t23 * t94;
t7 = -t19 * t93 + t21 * t97;
t4 = t14 * qJD(6) + t12 * t97 + t24 * t93;
t3 = -t15 * qJD(6) - t12 * t93 + t24 * t97;
t2 = -t8 * qJD(6) + t97 * t104 - t93 * t105;
t1 = -t93 * t104 - t97 * t105 - t21 * t152 + t19 * t153;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t57 * t125 + 0.2e1 * t23 * t35 + 0.2e1 * t24 * t36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t134, -t128, 0, 0, 0, 0, 0, -t103, t102 (t23 * t95 + t24 * t99 + (t35 * t99 - t36 * t95) * qJD(3)) * t89, t103, -t102, t53 * t125 + t23 * t54 - t24 * t52 + t35 * t56 - t36 * t47 + t46 * t57, 0, 0, 0, 0, 0, t24 * t59 + t36 * t40 + (-t11 * t95 + t115 * t161) * t89, t24 * t60 + t36 * t39 + (-t12 * t95 - t26 * t161) * t89, 0, 0, 0, 0, 0, -t11 * t113 - t115 * t17 + t14 * t40 + t3 * t59, t11 * t42 - t115 * t16 - t15 * t40 - t4 * t59; 0, 0, 0, 0, t67, 0.2e1 * (-t95 ^ 2 + t99 ^ 2) * t84 * qJD(3), 0.2e1 * t91 * t77, -0.2e1 * t91 * t135, 0, -0.2e1 * pkin(2) * t84 * t162 - 0.2e1 * t177, -0.2e1 * pkin(2) * t141 + 0.2e1 * t55 * t91 (-t47 * t99 + t56 * t95 + (t52 * t95 + t54 * t99) * qJD(3)) * t182, 0.2e1 * t177 + 0.2e1 * (-t53 * t162 + t46 * t99) * t89, -0.2e1 * t47 * t91 + 0.2e1 * (-t53 * t161 - t46 * t95) * t89, 0.2e1 * t46 * t53 + 0.2e1 * t47 * t52 + 0.2e1 * t54 * t56, 0.2e1 * t60 * t39, -0.2e1 * t39 * t59 - 0.2e1 * t40 * t60 (t60 * t161 + t39 * t95) * t182 (-t59 * t161 - t40 * t95) * t182, t67, 0.2e1 * t33 * t59 + 0.2e1 * t43 * t40 + 0.2e1 * (t10 * t95 + t117 * t161) * t89, 0.2e1 * t33 * t60 + 0.2e1 * t43 * t39 + 0.2e1 * (-t161 * t185 + t9 * t95) * t89, 0.2e1 * t42 * t16, 0.2e1 * t113 * t16 - 0.2e1 * t17 * t42, 0.2e1 * t16 * t59 + 0.2e1 * t40 * t42, 0.2e1 * t113 * t40 - 0.2e1 * t17 * t59, 0.2e1 * t59 * t40, -0.2e1 * t113 * t6 + 0.2e1 * t17 * t18 + 0.2e1 * t2 * t59 + 0.2e1 * t40 * t7, 0.2e1 * t1 * t59 + 0.2e1 * t16 * t18 - 0.2e1 * t40 * t8 + 0.2e1 * t42 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t24, 0, t23, t24, -pkin(3) * t23 + qJ(4) * t24 + qJD(4) * t36, 0, 0, 0, 0, 0, t36 * t154 + t24 * t94, -t36 * t156 + t24 * t98, 0, 0, 0, 0, 0 (t115 * t157 + t3) * t94 + (qJD(5) * t14 + t109) * t98 (t115 * t155 - t4) * t94 + (-qJD(5) * t15 - t108) * t98; 0, 0, 0, 0, 0, 0, t77, -t135, 0, -t56, t55 (t119 * qJD(3) + qJD(4) * t99) * t89, t56, -t55 + 0.2e1 * t82, -pkin(3) * t56 - qJ(4) * t47 - qJD(4) * t52, -t60 * t156 + t39 * t98, -t39 * t94 - t98 * t40 + (t59 * t94 - t60 * t98) * qJD(5), t50, t51, 0, t98 * t122 + qJ(4) * t40 + qJD(4) * t59 + t33 * t94 + (-t94 * t145 + t43 * t98) * qJD(5), -t94 * t122 + qJ(4) * t39 + qJD(4) * t60 + t33 * t98 + (-t98 * t145 - t43 * t94) * qJD(5), t16 * t173 + t42 * t62, t114 * t156 + (-t178 - t17 * t97 + (-t113 * t93 - t42 * t97) * qJD(6)) * t98 (t16 - t139) * t94 + (-t106 + t158) * t98 (-t17 + t140) * t94 + (-t107 + t159) * t98, t59 * t154 + t40 * t94, t28 * t59 + t48 * t40 + (t2 + (-t101 * t113 - t18 * t93) * qJD(5)) * t94 + (qJD(5) * t7 - t101 * t17 + t112) * t98, t27 * t59 - t49 * t40 + (t1 + (t101 * t42 - t18 * t97) * qJD(5)) * t94 + (-qJD(5) * t8 - t101 * t16 - t111) * t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181, qJ(4) * t181, -0.2e1 * t132, 0.2e1 * t126, 0, 0, 0, 0.2e1 * qJD(4) * t94 + 0.2e1 * t98 * t150, 0.2e1 * qJD(4) * t98 - 0.2e1 * t94 * t150, -0.2e1 * t87 * t132 - 0.2e1 * t88 * t133, 0.4e1 * t98 * t123 + 0.2e1 * t88 * t127, -0.2e1 * t94 * t137 - 0.2e1 * t170 * t155, 0.2e1 * t93 * t126 - 0.2e1 * t94 * t136, 0.2e1 * t132, -0.2e1 * t88 * t97 * t148 + 0.2e1 * t28 * t94 + 0.2e1 * (t48 + 0.2e1 * t144) * t154, 0.2e1 * t88 * t129 + 0.2e1 * t27 * t94 + 0.2e1 * (-t49 + 0.2e1 * t143) * t154; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t77, 0, 0, t56, 0, 0, 0, 0, 0, t50, t51, 0, 0, 0, 0, 0 (-t17 - t140) * t98 + (-t107 - t159) * t94 (-t16 - t139) * t98 + (t106 + t158) * t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169 * t152, t169 * t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, 0, 0, 0, 0, 0, t108, t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, -t40, t77, t10, t9, t42 * t152 + t178, -qJD(6) * t114 + t16 * t97 - t93 * t17, t107, -t106, 0, -pkin(5) * t17 - pkin(11) * t107 + t111, -pkin(5) * t16 + pkin(11) * t106 + t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, -t154, 0, -t94 * t149, -t130, -t98 * t127 - t123, -0.4e1 * t98 * t133 + t171 * t156, t63, -t61, 0 (-t142 - t183) * qJD(6) + (t120 * t93 - t143) * qJD(5) (t121 * t93 - t164 * t97) * qJD(6) + (-pkin(11) * t173 + (pkin(5) * t97 + t101 * t93) * t94) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t156, -t154, 0, 0, 0, 0, 0, t62, t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t133, -0.2e1 * t127, 0, 0, 0, t93 * t147, t97 * t147; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t17, t40, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, t64, t154, t28, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, -t153, 0, -pkin(11) * t152, pkin(11) * t153; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
