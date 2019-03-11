% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RRRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRRPRP8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:19:40
% EndTime: 2019-03-09 17:19:46
% DurationCPUTime: 2.19s
% Computational Cost: add. (2151->266), mult. (5062->464), div. (0->0), fcn. (4117->6), ass. (0->136)
t109 = sin(qJ(3));
t111 = cos(qJ(3));
t112 = cos(qJ(2));
t166 = qJD(3) * t112;
t150 = t111 * t166;
t110 = sin(qJ(2));
t97 = t110 * qJD(2);
t189 = t109 * t97 - t150;
t188 = -0.4e1 * t110;
t100 = t109 * qJ(4);
t187 = -t111 * pkin(3) - t100;
t106 = t111 ^ 2;
t171 = t109 ^ 2 - t106;
t133 = t171 * qJD(3);
t182 = cos(qJ(5));
t141 = qJD(5) * t182;
t142 = qJD(3) * t182;
t108 = sin(qJ(5));
t99 = qJD(3) * t111;
t147 = t108 * t99;
t165 = qJD(5) * t108;
t36 = -t111 * t165 + t147 + (t141 - t142) * t109;
t128 = t111 * t142;
t168 = qJD(3) * t109;
t59 = t108 * t109 + t182 * t111;
t35 = t59 * qJD(5) - t108 * t168 - t128;
t151 = t109 * t166;
t119 = t111 * t97 + t151;
t181 = pkin(8) * t112;
t125 = pkin(2) * t110 - t181;
t66 = t125 * qJD(2);
t179 = t110 * pkin(8);
t126 = -pkin(2) * t112 - t179;
t69 = -pkin(1) + t126;
t178 = t109 * t66 + t69 * t99;
t24 = t119 * pkin(7) - t178;
t173 = t110 * t111;
t88 = qJ(4) * t173;
t42 = -t88 + (pkin(3) * t109 + pkin(7)) * t110;
t175 = qJ(4) * t99 + t109 * qJD(4);
t47 = pkin(3) * t168 - t175;
t68 = -pkin(2) + t187;
t186 = (-t112 * t68 + t179) * qJD(2) - qJD(3) * t42 - t110 * t47;
t183 = pkin(8) - pkin(9);
t75 = t183 * t109;
t76 = t183 * t111;
t120 = t108 * t75 + t182 * t76;
t65 = t183 * t168;
t19 = t120 * qJD(5) - t108 * t65 - t183 * t128;
t103 = t112 * pkin(3);
t180 = pkin(9) * t110;
t174 = t109 * t112;
t90 = pkin(7) * t174;
t32 = t112 * pkin(4) + t103 + t90 + (-t69 - t180) * t111;
t172 = t111 * t112;
t91 = pkin(7) * t172;
t176 = t109 * t69 + t91;
t39 = -qJ(4) * t112 + t176;
t34 = t109 * t180 + t39;
t121 = t108 * t32 + t182 * t34;
t135 = -t189 * pkin(7) - t111 * t66 + t69 * t168;
t167 = qJD(3) * t110;
t152 = t109 * t167;
t184 = -pkin(3) - pkin(4);
t13 = pkin(9) * t152 + (-pkin(9) * t172 + t184 * t110) * qJD(2) + t135;
t93 = qJ(4) * t97;
t14 = t93 + (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t173 + (-qJD(4) + (-pkin(7) * qJD(3) + pkin(9) * qJD(2)) * t109) * t112 + t178;
t4 = -t121 * qJD(5) - t108 * t14 + t182 * t13;
t185 = 0.2e1 * qJD(4);
t163 = t112 * qJD(2);
t153 = t111 * t163;
t164 = t111 * qJD(4);
t177 = qJ(4) * t153 + t110 * t164;
t105 = t110 ^ 2;
t170 = -t112 ^ 2 + t105;
t169 = qJD(2) * t111;
t162 = t112 * qJD(4);
t161 = -0.2e1 * pkin(1) * qJD(2);
t160 = -0.2e1 * pkin(2) * qJD(3);
t159 = pkin(3) * t97;
t158 = pkin(8) * t168;
t157 = pkin(8) * t99;
t156 = pkin(7) * t163;
t155 = t184 * t109;
t154 = t182 * t109;
t145 = t109 * t99;
t144 = t110 * t163;
t143 = qJD(2) * t182;
t139 = -t108 * t34 + t182 * t32;
t137 = -t108 * t76 + t182 * t75;
t136 = t111 * t69 - t90;
t95 = t182 * t184;
t134 = -qJ(4) * t108 + t95;
t132 = t170 * qJD(2);
t131 = 0.2e1 * t144;
t56 = t111 * pkin(4) - t68;
t130 = t109 * t153;
t127 = -pkin(7) + t155;
t40 = t103 - t136;
t123 = -t109 * t39 + t111 * t40;
t3 = -t108 * t13 - t182 * t14 - t32 * t141 + t34 * t165;
t18 = -t75 * t141 - t183 * t147 + t76 * t165 + t182 * t65;
t41 = qJD(3) * t155 + t175;
t38 = t127 * t110 + t88;
t67 = t182 * qJ(4) + t108 * t184;
t118 = t109 * t163 + t110 * t99;
t115 = t187 * qJD(3) + t164;
t20 = -t24 + t93 - t162;
t22 = t135 - t159;
t113 = t123 * qJD(3) + t22 * t109 + t20 * t111;
t15 = (t184 * t111 - t100) * t167 + t127 * t163 + t177;
t82 = -0.2e1 * t144;
t81 = pkin(8) * t150;
t64 = -pkin(5) + t134;
t60 = -t108 * t111 + t154;
t48 = -t152 + t153;
t46 = t59 * t110;
t45 = t108 * t173 - t110 * t154;
t44 = t108 * qJD(4) + t67 * qJD(5);
t43 = qJ(4) * t165 - t182 * qJD(4) - qJD(5) * t95;
t37 = pkin(5) * t59 + t56;
t29 = -qJ(6) * t59 + t120;
t28 = -qJ(6) * t60 + t137;
t26 = t45 * pkin(5) + t38;
t23 = t118 * pkin(3) + qJ(4) * t152 + t156 - t177;
t21 = t36 * pkin(5) + t41;
t17 = t36 * t110 + t59 * t163;
t16 = t108 * t153 + t35 * t110 - t143 * t174;
t9 = -qJ(6) * t45 + t121;
t8 = pkin(5) * t112 - t46 * qJ(6) + t139;
t7 = t35 * qJ(6) - t60 * qJD(6) - t19;
t6 = -qJ(6) * t36 - qJD(6) * t59 - t18;
t5 = t16 * pkin(5) + t15;
t2 = -qJ(6) * t16 - qJD(6) * t45 - t3;
t1 = -pkin(5) * t97 - t17 * qJ(6) - t46 * qJD(6) + t4;
t10 = [0, 0, 0, t131, -0.2e1 * t132, 0, 0, 0, t110 * t161, t112 * t161, -0.2e1 * t105 * t145 + 0.2e1 * t106 * t144, 0.2e1 * t105 * t133 + t130 * t188, 0.2e1 * t110 * t151 + 0.2e1 * t170 * t169, -0.2e1 * t109 * t132 + 0.2e1 * t110 * t150, t82, 0.2e1 * t135 * t112 + 0.2e1 * t136 * t97 + 0.2e1 * (t105 * t99 + t109 * t131) * pkin(7), -0.2e1 * t24 * t112 - 0.2e1 * t176 * t97 + 0.2e1 * (-t105 * t168 + t111 * t131) * pkin(7), 0.2e1 * (qJD(2) * t109 * t42 + t22) * t112 + 0.2e1 * (-qJD(2) * t40 + t23 * t109 + t42 * t99) * t110, 0.2e1 * t123 * t163 + 0.2e1 * (-t109 * t20 + t111 * t22 + (-t109 * t40 - t111 * t39) * qJD(3)) * t110, 0.2e1 * (-t42 * t169 - t20) * t112 + 0.2e1 * (qJD(2) * t39 - t23 * t111 + t42 * t168) * t110, 0.2e1 * t20 * t39 + 0.2e1 * t22 * t40 + 0.2e1 * t23 * t42, 0.2e1 * t46 * t17, -0.2e1 * t16 * t46 - 0.2e1 * t17 * t45, 0.2e1 * t112 * t17 - 0.2e1 * t46 * t97, -0.2e1 * t112 * t16 + 0.2e1 * t45 * t97, t82, 0.2e1 * t4 * t112 - 0.2e1 * t139 * t97 + 0.2e1 * t15 * t45 + 0.2e1 * t38 * t16, 0.2e1 * t3 * t112 + 0.2e1 * t121 * t97 + 0.2e1 * t15 * t46 + 0.2e1 * t38 * t17, -0.2e1 * t1 * t46 - 0.2e1 * t16 * t9 - 0.2e1 * t17 * t8 - 0.2e1 * t2 * t45, 0.2e1 * t1 * t8 + 0.2e1 * t2 * t9 + 0.2e1 * t26 * t5; 0, 0, 0, 0, 0, t163, -t97, 0, -t156, pkin(7) * t97, -t110 * t133 + t130, t145 * t188 - t171 * t163, t189, t119, 0, t81 + (-pkin(2) * t111 + pkin(7) * t109) * t167 + (t126 * t109 - t91) * qJD(2) (pkin(7) * t173 + t125 * t109) * qJD(3) + (t126 * t111 + t90) * qJD(2), t81 + (t68 * t167 - t23) * t111 - t186 * t109, t113 (-t23 + (t110 * t68 + t181) * qJD(3)) * t109 + t186 * t111, t113 * pkin(8) + t23 * t68 + t42 * t47, t17 * t60 - t35 * t46, -t16 * t60 - t17 * t59 + t35 * t45 - t36 * t46, -t112 * t35 - t60 * t97, -t112 * t36 + t59 * t97, 0, -t19 * t112 - t137 * t97 + t15 * t59 + t56 * t16 + t38 * t36 + t41 * t45, t18 * t112 + t120 * t97 + t15 * t60 + t56 * t17 - t38 * t35 + t41 * t46, -t1 * t60 - t16 * t29 - t17 * t28 - t2 * t59 + t35 * t8 - t36 * t9 - t45 * t6 - t46 * t7, t1 * t28 + t2 * t29 + t21 * t26 + t37 * t5 + t6 * t9 + t7 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t145, -0.2e1 * t133, 0, 0, 0, t109 * t160, t111 * t160, -0.2e1 * t111 * t47 + 0.2e1 * t68 * t168, 0, -0.2e1 * t109 * t47 - 0.2e1 * t68 * t99, 0.2e1 * t68 * t47, -0.2e1 * t60 * t35, 0.2e1 * t35 * t59 - 0.2e1 * t36 * t60, 0, 0, 0, 0.2e1 * t36 * t56 + 0.2e1 * t41 * t59, -0.2e1 * t35 * t56 + 0.2e1 * t41 * t60, 0.2e1 * t28 * t35 - 0.2e1 * t29 * t36 - 0.2e1 * t59 * t6 - 0.2e1 * t60 * t7, 0.2e1 * t21 * t37 + 0.2e1 * t28 * t7 + 0.2e1 * t29 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t118, t97, -t135, t24, -t135 + 0.2e1 * t159 (-pkin(3) * t163 - qJ(4) * t167) * t111 + (-qJ(4) * t163 + (pkin(3) * qJD(3) - qJD(4)) * t110) * t109, -t24 + 0.2e1 * t93 - 0.2e1 * t162, -pkin(3) * t22 + qJ(4) * t20 + qJD(4) * t39, 0, 0, -t17, t16, t97, -t44 * t112 - t134 * t97 - t4, t112 * t43 + t67 * t97 - t3, -t16 * t67 - t17 * t64 + t43 * t45 + t44 * t46, t1 * t64 + t2 * t67 - t43 * t9 - t44 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, -t168, 0, -t157, t158, -t157, t115, -t158, t115 * pkin(8), 0, 0, t35, t36, 0, t19, -t18, t35 * t64 - t36 * t67 + t43 * t59 + t44 * t60, -t28 * t44 - t29 * t43 + t6 * t67 + t64 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, qJ(4) * t185, 0, 0, 0, 0, 0, 0.2e1 * t44, -0.2e1 * t43, 0, -0.2e1 * t43 * t67 - 0.2e1 * t44 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t97, t48, 0, t22, 0, 0, 0, 0, 0, -t110 * t143 - t112 * t165, t108 * t97 - t112 * t141, -t182 * t17 - t108 * t16 + (t108 * t46 - t182 * t45) * qJD(5), t1 * t182 + t2 * t108 + (-t108 * t8 + t182 * t9) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, 0, t157, 0, 0, 0, 0, 0, 0, 0, t182 * t35 - t108 * t36 + (t108 * t60 - t182 * t59) * qJD(5), t7 * t182 + t6 * t108 + (-t108 * t28 + t182 * t29) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, t141, 0, -t44 * t182 - t43 * t108 + (-t108 * t64 + t182 * t67) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, -t97, t4, t3, -pkin(5) * t17, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35, -t36, 0, -t19, t18, pkin(5) * t35, t7 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, t43, 0, -t44 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t165, -t141, 0, -pkin(5) * t165; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t10;
