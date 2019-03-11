% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:31:29
% EndTime: 2019-03-09 02:31:35
% DurationCPUTime: 2.06s
% Computational Cost: add. (1785->297), mult. (3798->426), div. (0->0), fcn. (2474->6), ass. (0->159)
t108 = sin(qJ(6));
t111 = cos(qJ(6));
t109 = sin(qJ(5));
t112 = cos(qJ(5));
t163 = t112 * qJD(4);
t113 = cos(qJ(4));
t173 = qJD(1) * t113;
t67 = t109 * t173 - t163;
t150 = t112 * t173;
t69 = qJD(4) * t109 + t150;
t128 = t108 * t67 - t111 * t69;
t23 = t108 * t69 + t111 * t67;
t213 = t128 * t23;
t160 = qJD(5) + qJD(6);
t71 = t108 * t112 + t109 * t111;
t116 = t160 * t71;
t106 = -pkin(7) + qJ(2);
t110 = sin(qJ(4));
t170 = qJD(4) * t113;
t212 = qJD(2) * t110 + t106 * t170;
t211 = t128 ^ 2 - t23 ^ 2;
t164 = qJD(6) * t111;
t165 = qJD(6) * t108;
t166 = qJD(5) * t113;
t145 = t109 * t166;
t205 = -t110 * t163 - t145;
t35 = qJD(1) * t205 + qJD(5) * t163;
t161 = qJD(1) * qJD(4);
t143 = t110 * t161;
t36 = qJD(5) * t69 - t109 * t143;
t6 = -t108 * t36 + t111 * t35 - t67 * t164 - t69 * t165;
t174 = qJD(1) * t110;
t95 = qJD(5) + t174;
t91 = qJD(6) + t95;
t210 = t23 * t91 + t6;
t107 = pkin(1) + qJ(3);
t76 = pkin(4) * t110 - pkin(8) * t113 + t107;
t47 = qJD(1) * t76 - qJD(2);
t100 = qJD(1) * qJ(2) + qJD(3);
t89 = -qJD(1) * pkin(7) + t100;
t75 = t110 * t89;
t57 = qJD(4) * pkin(8) + t75;
t19 = t109 * t47 + t112 * t57;
t15 = -pkin(9) * t67 + t19;
t11 = t15 * t165;
t58 = -qJD(4) * pkin(4) - t113 * t89;
t30 = pkin(5) * t67 + t58;
t209 = t23 * t30 + t11;
t117 = qJD(6) * t128 - t108 * t35 - t111 * t36;
t208 = -t128 * t91 + t117;
t131 = pkin(4) * t113 + pkin(8) * t110;
t65 = qJD(4) * t131 + qJD(3);
t48 = t65 * qJD(1);
t41 = t112 * t48;
t162 = qJD(1) * qJD(2);
t46 = t110 * t162 + t170 * t89;
t118 = -qJD(5) * t19 - t109 * t46 + t41;
t96 = t113 * t161;
t4 = pkin(5) * t96 - pkin(9) * t35 + t118;
t167 = qJD(5) * t112;
t159 = -t109 * t48 - t112 * t46 - t47 * t167;
t169 = qJD(5) * t109;
t123 = -t57 * t169 - t159;
t5 = -pkin(9) * t36 + t123;
t154 = -t108 * t5 + t111 * t4;
t18 = -t109 * t57 + t112 * t47;
t14 = -pkin(9) * t69 + t18;
t10 = pkin(5) * t95 + t14;
t195 = t111 * t15;
t2 = t10 * t108 + t195;
t207 = -t2 * qJD(6) + t30 * t128 + t154;
t121 = t71 * qJD(1);
t199 = t110 * t121 + t116;
t206 = t199 * t91;
t204 = qJD(1) * t107;
t203 = pkin(8) + pkin(9);
t202 = t67 * t95;
t201 = t69 * t95;
t151 = t109 * t174;
t179 = t111 * t112;
t183 = t108 * t109;
t200 = t108 * t151 - t111 * t167 - t112 * t164 + t160 * t183 - t174 * t179;
t178 = t112 * t113;
t72 = t131 * qJD(1);
t198 = t109 * t72 + t89 * t178;
t180 = t110 * t112;
t197 = t106 * t180 + t109 * t76;
t196 = t109 * t95;
t194 = t112 * t58;
t193 = t112 * t95;
t192 = t113 * t35;
t191 = t113 * t69;
t190 = t35 * t109;
t171 = qJD(4) * t110;
t45 = -t113 * t162 + t89 * t171;
t189 = t45 * t109;
t188 = t45 * t112;
t70 = -t179 + t183;
t187 = qJD(4) * t70;
t186 = qJD(4) * t71;
t185 = qJD(4) * t91;
t184 = t106 * t109;
t182 = t109 * t110;
t181 = t109 * t113;
t90 = -qJD(2) + t204;
t177 = qJD(2) - t90;
t105 = t113 ^ 2;
t176 = t110 ^ 2 - t105;
t114 = qJD(4) ^ 2;
t115 = qJD(1) ^ 2;
t175 = -t114 - t115;
t168 = qJD(5) * t110;
t158 = t95 * t182;
t157 = t95 * t180;
t102 = 0.2e1 * t162;
t156 = 0.2e1 * qJD(3) * qJD(1);
t155 = t95 * t169;
t153 = qJD(5) * t203;
t152 = t200 * t91;
t148 = t109 * t171;
t146 = t106 * t168;
t144 = t112 * t166;
t142 = qJD(6) * t10 + t5;
t141 = t106 * t95 + t57;
t139 = t109 * t65 + t212 * t112 + t76 * t167;
t138 = t177 * qJD(1);
t137 = qJD(1) + t168;
t136 = -t75 + (t151 + t169) * pkin(5);
t135 = t112 * t72 - t89 * t181;
t125 = pkin(5) * t113 + pkin(9) * t180;
t84 = t203 * t112;
t134 = qJD(1) * t125 + qJD(6) * t84 + t112 * t153 + t135;
t83 = t203 * t109;
t133 = pkin(9) * t151 + qJD(6) * t83 + t109 * t153 + t198;
t132 = t109 * t96 + t95 * t167;
t130 = qJD(2) + t90 + t204;
t62 = t112 * t76;
t21 = -pkin(9) * t178 + t62 + (pkin(5) - t184) * t110;
t27 = -pkin(9) * t181 + t197;
t129 = t108 * t21 + t111 * t27;
t127 = qJD(1) * t105 - t110 * t95;
t126 = -t106 * t114 + t156;
t124 = -pkin(8) * t170 + t110 * t58;
t122 = qJD(1) * t70;
t120 = -t144 + t148;
t99 = -pkin(5) * t112 - pkin(4);
t87 = t110 * t96;
t64 = (pkin(5) * t109 - t106) * t113;
t53 = t112 * t65;
t50 = t70 * t113;
t49 = t71 * t113;
t31 = -pkin(5) * t120 - qJD(2) * t113 + t106 * t171;
t16 = pkin(5) * t36 + t45;
t13 = -t165 * t181 + (t160 * t178 - t148) * t111 + t205 * t108;
t12 = -t113 * t116 + t171 * t70;
t9 = pkin(9) * t120 - t109 * t146 + t139;
t8 = -t112 * t146 + t53 + t125 * qJD(4) + ((pkin(9) * t113 - t76) * qJD(5) - t212) * t109;
t1 = t10 * t111 - t108 * t15;
t3 = [0, 0, 0, 0, t102, qJ(2) * t102, t102, t156, qJD(2) * t100 + qJD(3) * t90 + (qJ(2) * qJD(2) + qJD(3) * t107) * qJD(1), -0.2e1 * t87, 0.2e1 * t176 * t161, -t114 * t110, -t114 * t113, 0, t110 * t126 + t130 * t170, t113 * t126 - t130 * t171, -t69 * t145 + (-t171 * t69 + t192) * t112 (t109 * t69 + t112 * t67) * t171 + (-t190 - t112 * t36 + (t109 * t67 - t112 * t69) * qJD(5)) * t113, -t95 * t145 + t110 * t35 + (t112 * t127 + t191) * qJD(4), -t95 * t144 - t110 * t36 + (-t109 * t127 - t113 * t67) * qJD(4), t170 * t95 + t87 (-t169 * t76 + t53) * t95 + (qJD(4) * t106 * t67 + t41 - t141 * t167 + (-qJD(2) * t95 - qJD(4) * t58 - qJD(5) * t47 - t46) * t109) * t110 + (t58 * t167 - qJD(2) * t67 - t106 * t36 + t189 + (-t95 * t184 + (-t106 * t182 + t62) * qJD(1) + t18) * qJD(4)) * t113, -t139 * t95 + (t141 * t169 + (t106 * t69 - t194) * qJD(4) + t159) * t110 + (-t58 * t169 - qJD(2) * t69 - t106 * t35 + t188 + (-t197 * qJD(1) - t19) * qJD(4)) * t113, -t12 * t128 - t50 * t6, -t117 * t50 - t12 * t23 + t128 * t13 - t49 * t6, t110 * t6 + t12 * t91 + (-qJD(1) * t50 - t128) * t170, t110 * t117 - t13 * t91 + (-qJD(1) * t49 - t23) * t170, t170 * t91 + t87 (-t108 * t9 + t111 * t8) * t91 + t154 * t110 + t31 * t23 - t64 * t117 + t16 * t49 + t30 * t13 + (-t110 * t2 - t129 * t91) * qJD(6) + ((-t108 * t27 + t111 * t21) * qJD(1) + t1) * t170, t11 * t110 + t30 * t12 - t16 * t50 - t31 * t128 + t64 * t6 + (-(-qJD(6) * t27 + t8) * t91 - t4 * t110) * t108 + (-(qJD(6) * t21 + t9) * t91 - t142 * t110) * t111 + (-qJD(1) * t129 - t2) * t170; 0, 0, 0, 0, -t115, -t115 * qJ(2), -t115, 0 (-qJD(3) - t100) * qJD(1), 0, 0, 0, 0, 0, -0.2e1 * t96, 0.2e1 * t143, 0, 0, 0, 0, 0, t155 + (t158 + (t67 - t163) * t113) * qJD(1) (t157 + t191) * qJD(1) + t132, 0, 0, 0, 0, 0, t206 + (t23 + t187) * t173, -t152 + (-t128 + t186) * t173; 0, 0, 0, 0, 0, 0, 0, -t115, t138, 0, 0, 0, 0, 0, t175 * t110, t175 * t113, 0, 0, 0, 0, 0, -t113 * t36 - t137 * t193 + (t110 * t67 + (-t95 - t174) * t181) * qJD(4), -t192 + t137 * t196 + (-t95 * t178 + (t69 - t150) * t110) * qJD(4), 0, 0, 0, 0, 0, t91 * t122 + (-t185 * t71 + t117) * t113 + ((-t173 * t71 + t23) * qJD(4) + t160 * t91 * t70) * t110, t91 * t121 + (t185 * t70 - t6) * t113 + (t116 * t91 + (t113 * t122 - t128) * qJD(4)) * t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, t113 * t115 * t110, -t176 * t115, 0, 0, 0, t113 * t138, -t177 * t174, t69 * t193 + t190 (t35 - t202) * t112 + (-t36 - t201) * t109 (t157 - t191) * qJD(1) + t132, -t155 + (-t158 + (t67 + t163) * t113) * qJD(1), -t95 * t173, -pkin(4) * t36 - t188 - t135 * t95 - t67 * t75 + (-pkin(8) * t193 + t109 * t58) * qJD(5) + (t109 * t124 - t18 * t113) * qJD(1), -pkin(4) * t35 + t189 + t198 * t95 - t69 * t75 + (pkin(8) * t196 + t194) * qJD(5) + (t112 * t124 + t19 * t113) * qJD(1), t128 * t200 + t6 * t71, t117 * t71 + t128 * t199 + t200 * t23 - t6 * t70, -t152 + (t128 + t186) * t173, -t206 + (t23 - t187) * t173, -t91 * t173, t16 * t70 - t99 * t117 + (t108 * t133 - t111 * t134) * t91 + t199 * t30 + t136 * t23 + ((-t108 * t84 - t111 * t83) * qJD(4) - t1) * t173, t16 * t71 + t99 * t6 + (t108 * t134 + t111 * t133) * t91 - t200 * t30 - t136 * t128 + (-(-t108 * t83 + t111 * t84) * qJD(4) + t2) * t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69 * t67, -t67 ^ 2 + t69 ^ 2, t35 + t202, t201 - t36, t96, t19 * t95 - t58 * t69 + t118, t18 * t95 + t58 * t67 - t123, -t213, t211, t210, t208, t96 -(-t108 * t14 - t195) * t91 + (t111 * t96 - t165 * t91 - t23 * t69) * pkin(5) + t207 (-t15 * t91 - t4) * t108 + (t14 * t91 - t142) * t111 + (-t108 * t96 + t128 * t69 - t164 * t91) * pkin(5) + t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t213, t211, t210, t208, t96, t2 * t91 + t207, t1 * t91 - t108 * t4 - t111 * t142 + t209;];
tauc_reg  = t3;
