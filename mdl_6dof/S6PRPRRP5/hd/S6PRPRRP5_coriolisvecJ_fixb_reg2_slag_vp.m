% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6PRPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRP5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:16:43
% EndTime: 2019-03-08 20:16:50
% DurationCPUTime: 3.07s
% Computational Cost: add. (3201->358), mult. (7227->491), div. (0->0), fcn. (4873->8), ass. (0->205)
t116 = sin(qJ(4));
t181 = t116 * qJD(2);
t106 = qJD(5) + t181;
t115 = sin(qJ(5));
t118 = cos(qJ(5));
t180 = t118 * qJD(4);
t119 = cos(qJ(4));
t190 = qJD(2) * t119;
t90 = t115 * t190 - t180;
t210 = t90 * t106;
t183 = qJD(5) * t119;
t166 = t115 * t183;
t168 = t116 * t180;
t62 = qJD(2) * (t166 + t168) - qJD(5) * t180;
t247 = -t62 - t210;
t184 = qJD(5) * t118;
t185 = qJD(5) * t115;
t114 = cos(pkin(6));
t113 = sin(pkin(6));
t117 = sin(qJ(2));
t171 = qJD(2) * t113 * t117;
t187 = qJD(4) * t119;
t192 = qJD(1) * t116;
t121 = -pkin(2) - pkin(8);
t120 = cos(qJ(2));
t193 = qJD(1) * t113;
t162 = t120 * t193;
t142 = qJD(3) - t162;
t79 = t121 * qJD(2) + t142;
t41 = t79 * t187 + (-qJD(4) * t114 + t171) * t192;
t203 = t114 * t119;
t159 = qJD(1) * t203;
t59 = t116 * t79 + t159;
t51 = qJD(4) * pkin(9) + t59;
t144 = pkin(4) * t119 + pkin(9) * t116;
t86 = t144 * qJD(4) + qJD(3);
t60 = (t86 + t162) * qJD(2);
t163 = t117 * t193;
t97 = t116 * pkin(4) - t119 * pkin(9) + qJ(3);
t69 = qJD(2) * t97 + t163;
t155 = -t115 * t60 - t118 * t41 - t69 * t184 + t51 * t185;
t28 = -t115 * t51 + t118 * t69;
t246 = -t28 * t106 - t155;
t29 = t115 * t69 + t118 * t51;
t9 = -qJD(5) * t29 - t115 * t41 + t118 * t60;
t245 = -t29 * t106 - t9;
t198 = t116 * t121;
t71 = t115 * t97 + t118 * t198;
t201 = t115 * t117;
t244 = -t71 * qJD(5) + t118 * t86 - (-t116 * t201 + t118 * t120) * t193;
t186 = qJD(4) * t121;
t167 = t119 * t186;
t197 = t117 * t118;
t243 = t115 * t86 + t118 * t167 + t97 * t184 - (t115 * t120 + t116 * t197) * t193;
t160 = t118 * t190;
t182 = t115 * qJD(4);
t92 = t160 + t182;
t209 = t92 * t106;
t169 = t116 * t182;
t208 = qJD(5) * t92;
t63 = -qJD(2) * t169 + t208;
t242 = t63 + t209;
t58 = -t114 * t192 + t119 * t79;
t146 = qJD(2) * t163;
t96 = t119 * t146;
t42 = t59 * qJD(4) - t96;
t125 = -(t116 * t58 - t119 * t59) * qJD(4) + t41 * t116 - t42 * t119;
t241 = t92 ^ 2;
t240 = t90 * pkin(5);
t239 = pkin(5) * t115;
t205 = t113 * t120;
t77 = t114 * t116 + t119 * t205;
t238 = t42 * t77;
t237 = t92 * t90;
t236 = -qJ(6) - pkin(9);
t18 = -t92 * qJ(6) + t28;
t13 = t106 * pkin(5) + t18;
t235 = t13 - t18;
t157 = -t115 * t121 + pkin(5);
t179 = t118 * qJD(6);
t234 = qJ(6) * t168 + (qJ(6) * t185 + t157 * qJD(4) - t179) * t119 + t244;
t165 = t118 * t183;
t233 = -qJ(6) * t165 + (-qJD(6) * t119 + (qJ(6) * qJD(4) - qJD(5) * t121) * t116) * t115 + t243;
t156 = qJD(5) * t236;
t199 = t116 * t118;
t94 = t144 * qJD(2);
t35 = -t115 * t58 + t118 * t94;
t232 = (pkin(5) * t119 + qJ(6) * t199) * qJD(2) + t35 + t115 * qJD(6) - t118 * t156;
t36 = t115 * t94 + t118 * t58;
t231 = -t179 + t36 + (qJ(6) * t181 - t156) * t115;
t173 = t115 * t198;
t230 = -qJD(5) * t173 + t243;
t229 = -t115 * t167 + t244;
t228 = qJD(2) * pkin(2);
t50 = -qJD(4) * pkin(4) - t58;
t227 = t115 * t50;
t226 = t115 * t92;
t87 = (qJD(3) + t162) * qJD(2);
t225 = t117 * t87;
t224 = t118 * t50;
t223 = t119 * t62;
t222 = t119 * t63;
t191 = qJD(2) * qJ(3);
t95 = t163 + t191;
t220 = t120 * t95;
t19 = -t90 * qJ(6) + t29;
t219 = t19 * t106;
t24 = t63 * pkin(5) + t42;
t218 = t24 * t115;
t217 = t24 * t118;
t214 = t42 * t115;
t213 = t42 * t118;
t212 = t62 * t115;
t211 = t63 * t118;
t207 = t106 * t115;
t206 = t106 * t118;
t123 = qJD(2) ^ 2;
t204 = t113 * t123;
t202 = t115 * t116;
t200 = t115 * t119;
t196 = t118 * t119;
t111 = t116 ^ 2;
t112 = t119 ^ 2;
t195 = t111 - t112;
t122 = qJD(4) ^ 2;
t194 = -t122 - t123;
t189 = qJD(2) * t120;
t188 = qJD(4) * t116;
t178 = qJD(2) * qJD(4);
t176 = t116 * t205;
t175 = t117 * t204;
t174 = t120 * t204;
t172 = t119 * t123 * t116;
t170 = t113 * t189;
t164 = t106 * t190;
t161 = t116 * t186;
t107 = t119 * t178;
t158 = -qJD(6) - t240;
t154 = qJD(5) * t90 - t62;
t153 = -t63 + t208;
t152 = t106 + t181;
t151 = qJD(5) * t116 + qJD(2);
t150 = pkin(5) * t107;
t149 = t119 * t163;
t148 = t119 * t171;
t147 = t116 * t171;
t145 = t116 * t107;
t143 = -t95 + t163;
t141 = t115 * t19 + t118 * t13;
t140 = t115 * t13 - t118 * t19;
t139 = t115 * t29 + t118 * t28;
t138 = t115 * t28 - t118 * t29;
t135 = t87 * qJ(3) + t95 * qJD(3);
t134 = qJD(2) * t112 - t106 * t116;
t133 = -pkin(9) * t187 + t116 * t50;
t132 = t63 * qJ(6) + t155;
t78 = -t176 + t203;
t55 = t113 * t197 - t78 * t115;
t56 = t113 * t201 + t78 * t118;
t131 = t143 - t191;
t128 = t142 * qJD(2) - t121 * t122 + t87;
t126 = -t139 * qJD(5) - t9 * t115 - t118 * t155;
t124 = t62 * qJ(6) + t9;
t109 = -t118 * pkin(5) - pkin(4);
t102 = t236 * t118;
t101 = t236 * t115;
t89 = t142 - t228;
t88 = (-t121 + t239) * t119;
t85 = t90 ^ 2;
t84 = t118 * t97;
t72 = t152 * t187;
t70 = t84 - t173;
t68 = t92 * t149;
t67 = t90 * t149;
t64 = t161 + (t165 - t169) * pkin(5);
t54 = -qJD(4) * t176 + t114 * t187 - t148;
t53 = -qJD(4) * t77 + t147;
t52 = -qJ(6) * t200 + t71;
t46 = -qJ(6) * t196 + t157 * t116 + t84;
t44 = t159 + (-qJD(2) * t239 + t79) * t116;
t43 = -t85 + t241;
t40 = t209 - t63;
t39 = -t62 + t210;
t37 = -t158 + t50;
t32 = t106 * t184 + (t106 * t199 + (-t92 + t182) * t119) * qJD(2);
t31 = -t106 * t185 + (-t106 * t202 + (t90 + t180) * t119) * qJD(2);
t26 = t90 * t207 - t211;
t25 = t92 * t206 - t212;
t23 = t90 * t165 + (-t90 * t188 + t222) * t115;
t22 = -t92 * t166 + (-t92 * t188 - t223) * t118;
t21 = qJD(5) * t55 + t115 * t170 + t53 * t118;
t20 = -qJD(5) * t56 - t53 * t115 + t118 * t170;
t15 = -t106 * t165 - t63 * t116 + (-t115 * t134 - t119 * t90) * qJD(4);
t14 = -t106 * t166 - t62 * t116 + (t118 * t134 + t119 * t92) * qJD(4);
t12 = -t222 - t151 * t206 + (t116 * t90 - t152 * t200) * qJD(4);
t11 = t223 + t151 * t207 + (-t106 * t196 + (t92 - t160) * t116) * qJD(4);
t10 = -t242 * t115 + t247 * t118;
t7 = (t118 * t90 + t226) * t188 + (t212 - t211 + (t115 * t90 - t118 * t92) * qJD(5)) * t119;
t6 = (qJD(2) * t92 + t153 * t116 - t90 * t187) * t118 + (qJD(2) * t90 + t154 * t116 + t92 * t187) * t115;
t5 = -t90 * qJD(6) - t132;
t4 = -t92 * qJD(6) + t124 + t150;
t3 = -t21 * t106 - t56 * t107 + t54 * t92 - t77 * t62;
t2 = t20 * t106 + t55 * t107 + t54 * t90 + t77 * t63;
t1 = -t20 * t92 - t21 * t90 + t55 * t62 - t56 * t63;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t175, -t174, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, t174 (t225 + (t220 + (t89 - t162) * t117) * qJD(2)) * t113, 0, 0, 0, 0, 0, 0, t116 * t174 + (-t54 + t148) * qJD(4), t119 * t174 + (-t53 - t147) * qJD(4) (-t116 * t53 + t119 * t54 + (-t116 * t77 - t119 * t78) * qJD(4)) * qJD(2), t41 * t78 + t238 + t59 * t53 - t58 * t54 + (t189 * t95 + t225) * t113, 0, 0, 0, 0, 0, 0, t2, t3, t1, -t155 * t56 + t28 * t20 + t29 * t21 + t50 * t54 + t9 * t55 + t238, 0, 0, 0, 0, 0, 0, t2, t3, t1, t13 * t20 + t19 * t21 + t24 * t77 + t37 * t54 + t4 * t55 + t5 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * qJD(2) * qJD(3) (-t220 + (-t89 - t228) * t117) * t193 + t135, -0.2e1 * t145, 0.2e1 * t195 * t178, -t122 * t116, 0.2e1 * t145, -t122 * t119, 0, t116 * t128 - t131 * t187, t119 * t128 + t131 * t188 (t111 + t112) * t146 - t125 (-t220 + (-t116 * t59 - t119 * t58) * t117) * t193 + t125 * t121 + t135, t22, t7, t14, t23, t15, t72, t67 + t229 * t106 + (t9 + (t121 * t90 - t227) * qJD(4)) * t116 + (t50 * t184 + t214 - t121 * t63 + (qJD(2) * t70 + t28) * qJD(4)) * t119, t68 - t230 * t106 + (t155 + (t121 * t92 - t224) * qJD(4)) * t116 + (-t50 * t185 + t213 + t121 * t62 + (-qJD(2) * t71 - t29) * qJD(4)) * t119, t70 * t62 - t71 * t63 - t229 * t92 - t230 * t90 + t139 * t188 + (qJD(5) * t138 + t115 * t155 - t118 * t9) * t119, t50 * t161 + t9 * t70 - t155 * t71 + t230 * t29 + t229 * t28 + (-t121 * t42 + t163 * t50) * t119, t22, t7, t14, t23, t15, t72, t88 * t63 + t64 * t90 + t67 + (-t182 * t37 + t4) * t116 + t234 * t106 + (t37 * t184 + t218 + (qJD(2) * t46 + t13) * qJD(4)) * t119, -t88 * t62 + t64 * t92 + t68 + (-t180 * t37 - t5) * t116 - t233 * t106 + (-t37 * t185 + t217 + (-qJD(2) * t52 - t19) * qJD(4)) * t119, t46 * t62 - t52 * t63 - t234 * t92 - t233 * t90 + t141 * t188 + (qJD(5) * t140 - t115 * t5 - t118 * t4) * t119, t24 * t88 + t4 * t46 + t5 * t52 + (t64 + t149) * t37 + t233 * t19 + t234 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t123, t143 * qJD(2), 0, 0, 0, 0, 0, 0, t194 * t116, t194 * t119, 0, -t95 * qJD(2) + t125, 0, 0, 0, 0, 0, 0, t12, t11, t6, -t139 * qJD(2) + (-qJD(4) * t138 - t42) * t119 + (qJD(4) * t50 + t126) * t116, 0, 0, 0, 0, 0, 0, t12, t11, t6, -t141 * qJD(2) + (-qJD(4) * t140 - t24) * t119 + (qJD(4) * t37 - qJD(5) * t141 - t4 * t115 + t5 * t118) * t116; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, -t195 * t123, 0, -t172, 0, 0, -t190 * t95 + t96, -t143 * t181, 0, 0, t25, t10, t32, t26, t31, -t164, -pkin(4) * t63 - t35 * t106 - t213 - t59 * t90 + (-pkin(9) * t206 + t227) * qJD(5) + (t115 * t133 - t119 * t28) * qJD(2), pkin(4) * t62 + t36 * t106 + t214 - t59 * t92 + (pkin(9) * t207 + t224) * qJD(5) + (t118 * t133 + t119 * t29) * qJD(2), t35 * t92 + t36 * t90 + (pkin(9) * t153 + t246) * t118 + (pkin(9) * t154 + t245) * t115, -t42 * pkin(4) + pkin(9) * t126 - t28 * t35 - t29 * t36 - t50 * t59, t25, t10, t32, t26, t31, -t164, t109 * t63 - t217 - t44 * t90 - t232 * t106 + (t37 + t240) * t185 + (t37 * t202 + (qJD(4) * t101 - t13) * t119) * qJD(2), -t109 * t62 + t218 - t44 * t92 + t231 * t106 + (pkin(5) * t226 + t118 * t37) * qJD(5) + (t37 * t199 + (qJD(4) * t102 + t19) * t119) * qJD(2), t101 * t62 + t102 * t63 + t232 * t92 + t231 * t90 + (-t106 * t13 + t5) * t118 + (-t4 - t219) * t115, t4 * t101 - t5 * t102 + t24 * t109 + (pkin(5) * t185 - t44) * t37 - t231 * t19 - t232 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237, t43, t39, -t237, t40, t107, -t50 * t92 - t245, t50 * t90 - t246, 0, 0, t237, t43, t39, -t237, t40, t107, 0.2e1 * t150 + t219 + (t158 - t37) * t92 + t124, -t241 * pkin(5) + t18 * t106 + (qJD(6) + t37) * t90 + t132, t62 * pkin(5) - t235 * t90, t235 * t19 + (-t37 * t92 + t4) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t242, t247, -t85 - t241, t13 * t92 + t19 * t90 + t24;];
tauc_reg  = t8;
