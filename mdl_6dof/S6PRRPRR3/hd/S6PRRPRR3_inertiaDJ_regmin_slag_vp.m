% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6PRRPRR3_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR3_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR3_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_inertiaDJ_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 22:07:43
% EndTime: 2019-03-08 22:07:50
% DurationCPUTime: 2.08s
% Computational Cost: add. (2771->263), mult. (8616->509), div. (0->0), fcn. (8868->14), ass. (0->158)
t107 = sin(qJ(5));
t193 = -0.4e1 * t107;
t110 = cos(qJ(6));
t104 = cos(pkin(7));
t111 = cos(qJ(5));
t100 = sin(pkin(13));
t101 = sin(pkin(7));
t103 = cos(pkin(13));
t108 = sin(qJ(3));
t112 = cos(qJ(3));
t74 = (t100 * t112 + t103 * t108) * t101;
t132 = t111 * t104 - t107 * t74;
t106 = sin(qJ(6));
t165 = qJD(6) * t106;
t63 = t107 * t104 + t111 * t74;
t170 = qJD(3) * t112;
t150 = t101 * t170;
t171 = qJD(3) * t108;
t151 = t101 * t171;
t71 = -t100 * t151 + t103 * t150;
t46 = t63 * qJD(5) + t107 * t71;
t125 = -t110 * t46 - t132 * t165;
t192 = t125 * t107;
t102 = sin(pkin(6));
t105 = cos(pkin(6));
t109 = sin(qJ(2));
t175 = t109 * t112;
t113 = cos(qJ(2));
t176 = t108 * t113;
t121 = t104 * t176 + t175;
t119 = t121 * qJD(3);
t122 = -t104 * t175 - t176;
t114 = (t122 * qJD(2) - t119) * t102 - t105 * t151;
t180 = t101 * t108;
t190 = pkin(2) * t104;
t160 = t108 * t190;
t179 = t101 * t112;
t189 = pkin(9) + qJ(4);
t68 = t189 * t179 + t160;
t191 = -t68 * qJD(3) - qJD(4) * t180;
t144 = t189 * t108;
t64 = (pkin(2) * t112 + pkin(3)) * t104 - t101 * t144;
t48 = t100 * t64 + t103 * t68;
t37 = t104 * pkin(10) + t48;
t73 = t100 * t180 - t103 * t179;
t87 = (-pkin(3) * t112 - pkin(2)) * t101;
t53 = t73 * pkin(4) - t74 * pkin(10) + t87;
t134 = t107 * t53 + t111 * t37;
t91 = t170 * t190;
t56 = t91 + (-qJD(3) * t144 + qJD(4) * t112) * t101;
t31 = t191 * t100 + t103 * t56;
t70 = qJD(3) * t74;
t90 = pkin(3) * t151;
t50 = t70 * pkin(4) - t71 * pkin(10) + t90;
t10 = -t134 * qJD(5) - t107 * t31 + t111 * t50;
t98 = t110 ^ 2;
t188 = t106 ^ 2 - t98;
t97 = t107 ^ 2;
t187 = -t111 ^ 2 + t97;
t93 = t100 * pkin(3) + pkin(10);
t186 = t107 * t93;
t185 = t111 * t93;
t45 = t132 * qJD(5) + t111 * t71;
t51 = t106 * t63 - t110 * t73;
t17 = -t51 * qJD(6) + t106 * t70 + t110 * t45;
t184 = t17 * t106;
t183 = t17 * t110;
t182 = qJD(5) * t51;
t181 = qJD(6) * t97;
t178 = t105 * t104;
t177 = t108 * t109;
t174 = t110 * t111;
t173 = t112 * t113;
t172 = qJD(2) * t102;
t169 = qJD(5) * t106;
t168 = qJD(5) * t107;
t167 = qJD(5) * t110;
t166 = qJD(5) * t111;
t164 = qJD(6) * t110;
t163 = qJD(6) * t111;
t162 = -0.2e1 * pkin(5) * qJD(6);
t94 = -t103 * pkin(3) - pkin(4);
t161 = 0.2e1 * qJD(5) * t94;
t159 = t93 * t181;
t158 = t106 * t185;
t157 = t93 * t174;
t95 = t101 ^ 2;
t156 = t95 * t170;
t155 = t132 * t169;
t154 = t132 * t167;
t153 = t106 * t163;
t152 = t110 * t163;
t148 = t109 * t172;
t147 = t106 * t164;
t146 = t107 * t166;
t145 = t110 * t166;
t30 = t100 * t56 - t103 * t191;
t47 = -t100 * t68 + t103 * t64;
t143 = t188 * qJD(6);
t142 = t187 * qJD(5);
t141 = t101 * t148;
t139 = t106 * t145;
t138 = -t111 * pkin(5) - t107 * pkin(11);
t137 = pkin(5) * t107 - pkin(11) * t111;
t20 = t73 * pkin(11) + t134;
t36 = -t104 * pkin(4) - t47;
t26 = -pkin(5) * t132 - t63 * pkin(11) + t36;
t8 = t106 * t26 + t110 * t20;
t123 = t104 * t173 - t177;
t116 = t123 * t102 + t105 * t179;
t58 = t121 * t102 + t105 * t180;
t34 = t100 * t116 + t103 * t58;
t78 = -t102 * t113 * t101 + t178;
t28 = t78 * t107 + t111 * t34;
t33 = t100 * t58 - t103 * t116;
t16 = t106 * t33 + t110 * t28;
t52 = t106 * t73 + t110 * t63;
t136 = -t106 * t52 - t110 * t51;
t133 = -t107 * t37 + t111 * t53;
t19 = -t73 * pkin(5) - t133;
t6 = -t70 * pkin(5) - t10;
t131 = t6 * t106 + t19 * t164;
t130 = -t6 * t110 + t19 * t165;
t84 = t138 + t94;
t61 = t106 * t84 + t157;
t129 = t107 * t70 + t73 * t166;
t44 = t105 * t150 + (t123 * qJD(3) + (-t104 * t177 + t173) * qJD(2)) * t102;
t25 = t114 * t100 + t103 * t44;
t13 = t28 * qJD(5) + t107 * t25 - t111 * t141;
t27 = t107 * t34 - t78 * t111;
t128 = t13 * t106 + t27 * t164;
t127 = -t13 * t110 + t27 * t165;
t126 = t106 * t46 - t132 * t164;
t124 = t137 * t106;
t9 = -t107 * t50 - t111 * t31 - t53 * t166 + t37 * t168;
t80 = t107 * t167 + t153;
t120 = t46 * pkin(5) - t45 * pkin(11) + t30;
t118 = -t70 * pkin(11) + t9;
t82 = t106 * t168 - t152;
t81 = -t106 * t166 - t107 * t164;
t79 = t107 * t165 - t145;
t77 = (-pkin(9) * t179 - t160) * qJD(3);
t76 = pkin(9) * t151 - t91;
t60 = t110 * t84 - t158;
t54 = t111 * t70 - t73 * t168;
t41 = t52 * t168;
t40 = -t61 * qJD(6) + (t106 * t186 + t110 * t137) * qJD(5);
t39 = -qJD(5) * t124 - t84 * t164 + t80 * t93;
t24 = t100 * t44 - t103 * t114;
t18 = qJD(6) * t52 + t106 * t45 - t110 * t70;
t15 = -t106 * t28 + t110 * t33;
t12 = -t107 * t141 - t111 * t25 - t78 * t166 + t34 * t168;
t7 = -t106 * t20 + t110 * t26;
t4 = -t16 * qJD(6) + t106 * t12 + t110 * t24;
t3 = -t106 * t24 + t110 * t12 - t33 * t164 + t28 * t165;
t2 = -t8 * qJD(6) + t106 * t118 + t110 * t120;
t1 = -t106 * t120 + t110 * t118 - t26 * t164 + t20 * t165;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t141 * t78 + 0.2e1 * t33 * t24 + 0.2e1 * t34 * t25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, -t148, -t113 * t172, 0, 0, 0, 0, 0 (t78 - t178) * t151 + (-t104 * t119 + (t104 * t122 - t175 * t95) * qJD(2)) * t102, t95 * t108 * t148 - t44 * t104 + t150 * t78, t24 * t74 - t25 * t73 + t33 * t71 - t34 * t70, -t24 * t47 + t25 * t48 + t33 * t30 + t34 * t31 + (pkin(3) * t171 * t78 + t148 * t87) * t101, 0, 0, 0, 0, 0, -t13 * t73 - t132 * t24 - t27 * t70 + t33 * t46, t12 * t73 + t24 * t63 - t28 * t70 + t33 * t45, 0, 0, 0, 0, 0, t13 * t51 - t132 * t4 + t15 * t46 + t27 * t18, t13 * t52 - t132 * t3 - t16 * t46 + t27 * t17; 0, 0, 0, 0, 0.2e1 * t108 * t156, 0.2e1 * (-t108 ^ 2 + t112 ^ 2) * t95 * qJD(3), 0.2e1 * t104 * t150, -0.2e1 * t104 * t151, 0, -0.2e1 * pkin(2) * t171 * t95 + 0.2e1 * t77 * t104, -0.2e1 * pkin(2) * t156 + 0.2e1 * t76 * t104, 0.2e1 * t30 * t74 - 0.2e1 * t31 * t73 - 0.2e1 * t47 * t71 - 0.2e1 * t48 * t70, -0.2e1 * t47 * t30 + 0.2e1 * t48 * t31 + 0.2e1 * t87 * t90, 0.2e1 * t63 * t45, 0.2e1 * t132 * t45 - 0.2e1 * t63 * t46, 0.2e1 * t45 * t73 + 0.2e1 * t63 * t70, 0.2e1 * t132 * t70 - 0.2e1 * t46 * t73, 0.2e1 * t73 * t70, 0.2e1 * t10 * t73 - 0.2e1 * t132 * t30 + 0.2e1 * t133 * t70 + 0.2e1 * t36 * t46, -0.2e1 * t134 * t70 + 0.2e1 * t30 * t63 + 0.2e1 * t36 * t45 + 0.2e1 * t9 * t73, 0.2e1 * t52 * t17, -0.2e1 * t17 * t51 - 0.2e1 * t52 * t18, -0.2e1 * t132 * t17 + 0.2e1 * t52 * t46, 0.2e1 * t132 * t18 - 0.2e1 * t51 * t46, -0.2e1 * t132 * t46, -0.2e1 * t132 * t2 + 0.2e1 * t19 * t18 + 0.2e1 * t7 * t46 + 0.2e1 * t6 * t51, -0.2e1 * t1 * t132 + 0.2e1 * t19 * t17 - 0.2e1 * t8 * t46 + 0.2e1 * t6 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, -t44, 0 (t100 * t25 - t103 * t24) * pkin(3), 0, 0, 0, 0, 0, -t24 * t111 + t168 * t33, t24 * t107 + t166 * t33, 0, 0, 0, 0, 0 (t169 * t27 - t4) * t111 + (qJD(5) * t15 + t128) * t107 (t167 * t27 - t3) * t111 + (-qJD(5) * t16 - t127) * t107; 0, 0, 0, 0, 0, 0, t150, -t151, 0, t77, t76 (-t100 * t70 - t103 * t71) * pkin(3) (t100 * t31 - t103 * t30) * pkin(3), t45 * t107 + t166 * t63, -t107 * t46 + t45 * t111 + (-t107 * t63 + t111 * t132) * qJD(5), t129, t54, 0, -t70 * t186 - t30 * t111 + t94 * t46 + (t107 * t36 - t73 * t185) * qJD(5), -t70 * t185 + t30 * t107 + t94 * t45 + (t111 * t36 + t73 * t186) * qJD(5), t52 * t145 + (-t165 * t52 + t183) * t107, t136 * t166 + (-t184 - t110 * t18 + (t106 * t51 - t110 * t52) * qJD(6)) * t107, t41 + (-t17 - t154) * t111 - t192 (t18 + t155) * t111 + (-t126 - t182) * t107, -t46 * t111 - t132 * t168, -t40 * t132 + t60 * t46 + (-t2 + (t106 * t19 + t51 * t93) * qJD(5)) * t111 + (qJD(5) * t7 + t18 * t93 + t131) * t107, -t39 * t132 - t61 * t46 + (-t1 + (t110 * t19 + t52 * t93) * qJD(5)) * t111 + (-qJD(5) * t8 + t17 * t93 - t130) * t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t146, -0.2e1 * t142, 0, 0, 0, t107 * t161, t111 * t161, 0.2e1 * t146 * t98 - 0.2e1 * t147 * t97, t139 * t193 + 0.2e1 * t188 * t181, 0.2e1 * t107 * t153 + 0.2e1 * t187 * t167, -0.2e1 * t106 * t142 + 0.2e1 * t107 * t152, -0.2e1 * t146, 0.2e1 * t110 * t159 - 0.2e1 * t40 * t111 + 0.2e1 * (t60 + 0.2e1 * t158) * t168, -0.2e1 * t106 * t159 - 0.2e1 * t39 * t111 + 0.2e1 * (-t61 + 0.2e1 * t157) * t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, 0, 0, 0, 0, 0, t54, -t129, 0, 0, 0, 0, 0 (-t18 + t155) * t111 + (-t126 + t182) * t107, t41 + (-t17 + t154) * t111 + t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, t12, 0, 0, 0, 0, 0, t127, t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t46, t70, t10, t9, t164 * t52 + t184, qJD(6) * t136 - t106 * t18 + t183, t126, -t125, 0, -pkin(5) * t18 - pkin(11) * t126 + t130, -pkin(5) * t17 + pkin(11) * t125 + t131; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166, -t168, 0, -t93 * t166, t93 * t168, -t107 * t143 + t139, t147 * t193 - t188 * t166, t82, t80, 0 (pkin(11) * t174 + (-t110 * pkin(5) + t106 * t93) * t107) * qJD(6) + (t106 * t138 - t157) * qJD(5) (t110 * t186 + t124) * qJD(6) + (t110 * t138 + t158) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t168, -t166, 0, 0, 0, 0, 0, -t80, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t147, -0.2e1 * t143, 0, 0, 0, t106 * t162, t110 * t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t18, t46, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, t81, t168, t40, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t164, -t165, 0, -pkin(11) * t164, pkin(11) * t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t5;
