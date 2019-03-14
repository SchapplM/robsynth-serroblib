% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S6RRPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x(6*10)]
%   inertial parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RRPPRR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR4_inertiaDJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR4_inertiaDJ_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPPRR4_inertiaDJ_reg2_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:05:03
% EndTime: 2019-03-09 09:05:16
% DurationCPUTime: 4.56s
% Computational Cost: add. (6706->340), mult. (18239->637), div. (0->0), fcn. (17992->10), ass. (0->199)
t232 = pkin(3) + pkin(9);
t85 = sin(qJ(6));
t80 = t85 ^ 2;
t88 = cos(qJ(6));
t82 = t88 ^ 2;
t200 = t80 + t82;
t194 = sin(pkin(11));
t87 = sin(qJ(2));
t150 = t87 * t194;
t90 = cos(qJ(2));
t157 = -t90 * pkin(2) - pkin(1);
t195 = cos(pkin(11));
t118 = t194 * t90 + t195 * t87;
t84 = sin(pkin(6));
t60 = t118 * t84;
t197 = t60 * qJ(4);
t209 = t84 * t90;
t73 = t195 * t209;
t100 = t197 + t232 * t73 - (t150 * t232 + t157) * t84;
t196 = cos(pkin(6));
t156 = pkin(1) * t196;
t141 = t90 * t156;
t203 = pkin(8) + qJ(3);
t155 = t84 * t203;
t114 = -t155 * t87 + t141;
t190 = qJD(3) * t84;
t101 = qJD(2) * t114 + t190 * t90;
t142 = t87 * t156;
t102 = -t87 * t190 + (-t155 * t90 - t142) * qJD(2);
t25 = t101 * t194 - t102 * t195;
t128 = t150 * t84 - t73;
t58 = t128 * qJD(2);
t229 = -t58 * pkin(4) + qJD(5) * t100 + t25;
t112 = qJD(2) * t60;
t193 = qJD(2) * t84;
t113 = t118 * t193;
t192 = qJD(2) * t87;
t162 = t84 * t192;
t146 = pkin(2) * t162;
t24 = pkin(3) * t113 + t58 * qJ(4) - t60 * qJD(4) + t146;
t104 = pkin(2) * t196 + t114;
t115 = t203 * t209 + t142;
t35 = t104 * t195 - t115 * t194;
t30 = -pkin(3) * t196 - t35;
t96 = t60 * pkin(4) - pkin(9) * t196 + t30;
t230 = -pkin(9) * t112 - qJD(5) * t96 - t24;
t86 = sin(qJ(5));
t89 = cos(qJ(5));
t5 = -t229 * t86 + t230 * t89;
t36 = t104 * t194 + t115 * t195;
t29 = -qJ(4) * t196 - t36;
t21 = -pkin(4) * t128 - t29;
t48 = -t128 * t89 + t196 * t86;
t49 = t128 * t86 + t196 * t89;
t94 = t48 * pkin(5) - t49 * pkin(10) + t21;
t231 = t58 * pkin(10) - qJD(6) * t94 + t5;
t228 = t128 * t58;
t223 = pkin(5) * t89;
t139 = pkin(10) * t86 + t223;
t227 = t139 * t88;
t183 = qJD(6) * t85;
t185 = qJD(5) * t89;
t64 = t183 * t86 - t185 * t88;
t201 = t80 - t82;
t148 = qJD(6) * t201;
t34 = qJD(5) * t49 - t112 * t89;
t31 = t89 * t34;
t186 = qJD(5) * t88;
t160 = t86 * t186;
t181 = qJD(6) * t89;
t164 = t85 * t181;
t65 = t160 + t164;
t226 = t31 * t88 - t48 * t65;
t13 = -t100 * t89 + t86 * t96;
t11 = pkin(10) * t60 + t13;
t7 = -t85 * t11 + t88 * t94;
t8 = t88 * t11 + t85 * t94;
t135 = t7 * t85 - t8 * t88;
t6 = t229 * t89 + t230 * t86;
t4 = t58 * pkin(5) - t6;
t225 = qJD(5) * t135 + t4;
t224 = 0.2e1 * qJD(4);
t222 = t86 * pkin(5);
t33 = qJD(5) * t48 - t86 * t112;
t38 = t49 * t88 + t60 * t85;
t17 = qJD(6) * t38 - t33 * t85 + t58 * t88;
t221 = t17 * t86;
t220 = t17 * t88;
t18 = -t85 * t58 - t49 * t183 + (qJD(6) * t60 - t33) * t88;
t219 = t18 * t85;
t218 = t25 * t60;
t217 = t33 * t86;
t216 = t33 * t89;
t215 = t34 * t86;
t37 = t49 * t85 - t60 * t88;
t214 = t37 * t88;
t213 = t38 * t85;
t212 = t48 * t86;
t211 = t49 * t89;
t77 = -pkin(2) * t195 - pkin(3);
t75 = -pkin(9) + t77;
t210 = t75 * t86;
t208 = t85 * t75;
t207 = t85 * t89;
t206 = t86 * t58;
t205 = t88 * t89;
t204 = t89 * t58;
t202 = t18 * t86 + t185 * t38;
t81 = t86 ^ 2;
t83 = t89 ^ 2;
t199 = t81 - t83;
t198 = t81 + t83;
t191 = qJD(2) * t90;
t189 = qJD(5) * t37;
t188 = qJD(5) * t85;
t187 = qJD(5) * t86;
t184 = qJD(6) * t37;
t182 = qJD(6) * t88;
t180 = 0.2e1 * t48 * t34;
t42 = -0.2e1 * t60 * t58;
t179 = -0.2e1 * pkin(5) * qJD(6);
t177 = t86 * t208;
t176 = t75 * t207;
t175 = t88 * t210;
t76 = pkin(2) * t194 + qJ(4);
t174 = t76 * t224;
t79 = t84 ^ 2;
t173 = t79 * t191;
t172 = t37 * t187;
t171 = t38 * t187;
t170 = t85 * t187;
t169 = t75 * t187;
t168 = t75 * t185;
t167 = t38 * t181;
t166 = qJD(6) * t75 * t83;
t163 = t88 * t181;
t161 = t85 * t182;
t158 = t86 * t185;
t154 = t200 * t86;
t153 = t200 * t89;
t152 = t25 * t196;
t151 = t58 * t196;
t147 = t196 * qJD(4);
t145 = t85 * t160;
t144 = t83 * t161;
t143 = t87 * t173;
t138 = -pkin(10) * t89 + t222;
t137 = -t5 * t86 + t6 * t89;
t136 = t7 * t88 + t8 * t85;
t134 = t196 * t193;
t12 = t100 * t86 + t89 * t96;
t133 = -t12 * t86 + t13 * t89;
t132 = -t215 + t216;
t116 = t138 + t76;
t111 = t88 * t116;
t46 = t111 - t177;
t47 = t116 * t85 + t175;
t131 = t46 * t88 + t47 * t85;
t130 = t46 * t85 - t47 * t88;
t129 = -t48 * t89 + t49 * t86;
t10 = -t60 * pkin(5) - t12;
t127 = t10 * t182 + t4 * t85;
t126 = t10 * t183 - t4 * t88;
t125 = t185 * t48 + t215;
t40 = -t185 * t60 + t206;
t124 = t187 * t60 + t204;
t123 = t182 * t48 + t34 * t85;
t122 = t183 * t48 - t34 * t88;
t67 = -t163 + t170;
t120 = -pkin(8) * t209 - t142;
t119 = pkin(8) * t84 * t87 - t141;
t117 = -t123 + t189;
t110 = -0.2e1 * t113;
t109 = 0.2e1 * t112;
t26 = t101 * t195 + t102 * t194;
t22 = -t147 - t26;
t19 = -pkin(4) * t112 - t22;
t92 = t34 * pkin(5) + t33 * pkin(10) + t19;
t1 = t11 * t183 + t231 * t88 - t85 * t92;
t2 = -t11 * t182 + t231 * t85 + t88 * t92;
t108 = -qJD(6) * t136 - t1 * t88 - t2 * t85;
t107 = qJD(5) * t133 + t137;
t106 = -t220 + t219 + (t37 * t85 + t38 * t88) * qJD(6);
t27 = -qJD(6) * t111 - t85 * (qJD(5) * t139 + qJD(4)) + t64 * t75;
t28 = t88 * qJD(4) - t47 * qJD(6) + (-t176 + t227) * qJD(5);
t105 = -qJD(6) * t131 - t27 * t88 - t28 * t85;
t103 = qJD(5) * t10 + t108;
t78 = qJD(5) * t81;
t74 = 0.2e1 * t158;
t71 = t157 * t84;
t66 = t182 * t86 + t185 * t85;
t62 = t120 * qJD(2);
t61 = t119 * qJD(2);
t51 = t148 * t89 + t145;
t50 = (-0.1e1 + t200) * t158;
t44 = t48 * t170;
t39 = -t73 * pkin(3) - t197 + (pkin(3) * t150 + t157) * t84;
t14 = t17 * t205;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t143, 0.2e1 * (-t87 ^ 2 + t90 ^ 2) * t79 * qJD(2), 0.2e1 * t90 * t134, -0.2e1 * t143, -0.2e1 * t87 * t134, 0, -0.2e1 * pkin(1) * t192 * t79 + 0.2e1 * t196 * t62, -0.2e1 * pkin(1) * t173 + 0.2e1 * t196 * t61, 0.2e1 * (-t61 * t90 - t62 * t87 + (t119 * t90 + t120 * t87) * qJD(2)) * t84, -0.2e1 * t119 * t62 + 0.2e1 * t120 * t61, t42, -0.2e1 * t113 * t60 + 0.2e1 * t228, -0.2e1 * t151, -t128 * t110, t196 * t110, 0, 0.2e1 * t112 * t71 + 0.2e1 * t128 * t146 - 0.2e1 * t152, 0.2e1 * t146 * t60 - 0.2e1 * t196 * t26 - 0.2e1 * t58 * t71, -0.2e1 * t113 * t36 - 0.2e1 * t128 * t26 + 0.2e1 * t35 * t58 + 0.2e1 * t218, 0.2e1 * t146 * t71 - 0.2e1 * t25 * t35 + 0.2e1 * t26 * t36, 0, 0.2e1 * t151, t196 * t109, t42, -0.2e1 * t112 * t60 + 0.2e1 * t228, t128 * t109, 0.2e1 * t112 * t29 + 0.2e1 * t128 * t22 - 0.2e1 * t30 * t58 + 0.2e1 * t218, -0.2e1 * t113 * t39 - 0.2e1 * t128 * t24 + 0.2e1 * t152, -0.2e1 * t196 * t22 - 0.2e1 * t24 * t60 + 0.2e1 * t39 * t58, 0.2e1 * t22 * t29 + 0.2e1 * t24 * t39 + 0.2e1 * t25 * t30, -0.2e1 * t49 * t33, 0.2e1 * t33 * t48 - 0.2e1 * t34 * t49, -0.2e1 * t33 * t60 - 0.2e1 * t49 * t58, t180, -0.2e1 * t34 * t60 + 0.2e1 * t48 * t58, t42, -0.2e1 * t12 * t58 + 0.2e1 * t19 * t48 + 0.2e1 * t21 * t34 + 0.2e1 * t6 * t60, 0.2e1 * t13 * t58 + 0.2e1 * t19 * t49 - 0.2e1 * t21 * t33 + 0.2e1 * t5 * t60, 0.2e1 * t12 * t33 - 0.2e1 * t13 * t34 + 0.2e1 * t48 * t5 - 0.2e1 * t49 * t6, 0.2e1 * t12 * t6 - 0.2e1 * t13 * t5 + 0.2e1 * t19 * t21, 0.2e1 * t38 * t18, -0.2e1 * t17 * t38 - 0.2e1 * t18 * t37, 0.2e1 * t18 * t48 + 0.2e1 * t34 * t38, 0.2e1 * t37 * t17, -0.2e1 * t17 * t48 - 0.2e1 * t34 * t37, t180, 0.2e1 * t10 * t17 + 0.2e1 * t2 * t48 + 0.2e1 * t34 * t7 + 0.2e1 * t37 * t4, 0.2e1 * t1 * t48 + 0.2e1 * t10 * t18 - 0.2e1 * t34 * t8 + 0.2e1 * t38 * t4, 0.2e1 * t1 * t37 - 0.2e1 * t17 * t8 - 0.2e1 * t18 * t7 - 0.2e1 * t2 * t38, -0.2e1 * t1 * t8 + 0.2e1 * t10 * t4 + 0.2e1 * t2 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84 * t191, 0, -t162, 0, t62, t61, 0, 0, 0, 0, -t58, 0, -t113, 0, -t25, -t26 (-t113 * t194 + t195 * t58) * pkin(2) (t194 * t26 - t195 * t25) * pkin(2), 0, t58, t112, 0, 0, 0, -qJD(4) * t128 - t112 * t76 - t77 * t58, t25, 0.2e1 * t147 + t26, -qJD(4) * t29 - t22 * t76 + t25 * t77, -t187 * t49 - t216, t217 - t31 + (-t211 + t212) * qJD(5), -t124, t125, t40, 0, -t75 * t204 + qJD(4) * t48 + t19 * t86 + t34 * t76 + (t21 * t89 - t210 * t60) * qJD(5), t75 * t206 + qJD(4) * t49 + t19 * t89 - t33 * t76 + (-t60 * t75 * t89 - t21 * t86) * qJD(5), t132 * t75 + (t129 * t75 - t133) * qJD(5) - t137, qJD(4) * t21 + t107 * t75 + t19 * t76, t18 * t205 - t38 * t65, -t14 + (-t167 + t172) * t88 + (t171 + (-t18 + t184) * t89) * t85, t202 + t226, t17 * t207 - t37 * t67, -t221 + t44 + (-t123 - t189) * t89, t125, t28 * t48 + t34 * t46 + (t2 + (-t10 * t85 + t37 * t75) * qJD(5)) * t86 + (qJD(5) * t7 - t17 * t75 + t127) * t89, t27 * t48 - t34 * t47 + (t1 + (-t10 * t88 + t38 * t75) * qJD(5)) * t86 + (-qJD(5) * t8 - t18 * t75 - t126) * t89, -t17 * t47 - t18 * t46 + t27 * t37 - t28 * t38 + t136 * t187 + (qJD(6) * t135 + t1 * t85 - t2 * t88) * t89, -t1 * t47 + t2 * t46 - t27 * t8 + t28 * t7 + (t10 * t187 - t4 * t89) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, t174, -0.2e1 * t158, -0.2e1 * qJD(5) * t83 + 0.2e1 * t78, 0, t74, 0, 0, 0.2e1 * qJD(4) * t86 + 0.2e1 * t185 * t76, 0.2e1 * qJD(4) * t89 - 0.2e1 * t187 * t76, 0, t174, -0.2e1 * t158 * t82 - 0.2e1 * t144, 0.4e1 * t145 * t89 + 0.2e1 * t148 * t83, -0.2e1 * t164 * t86 - 0.2e1 * t186 * t199, -0.2e1 * t158 * t80 + 0.2e1 * t144, -0.2e1 * t163 * t86 + 0.2e1 * t188 * t199, t74, -0.2e1 * t88 * t166 + 0.2e1 * t28 * t86 + 0.2e1 * (t46 + 0.2e1 * t177) * t185, 0.2e1 * t85 * t166 + 0.2e1 * t27 * t86 + 0.2e1 * (-t47 + 0.2e1 * t175) * t185, 0.2e1 * t131 * t187 + 0.2e1 * (qJD(6) * t130 + t27 * t85 - t28 * t88) * t89, -0.2e1 * t158 * t75 ^ 2 - 0.2e1 * t27 * t47 + 0.2e1 * t28 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, -t58, 0, t146, 0, 0, 0, 0, 0, 0, 0, -t113, t58, t24, 0, 0, 0, 0, 0, 0, t40, t124, -t217 - t31 + (t211 + t212) * qJD(5), -t5 * t89 - t6 * t86 + (-t12 * t89 - t13 * t86) * qJD(5), 0, 0, 0, 0, 0, 0, t117 * t89 + t221 + t44, t202 - t226, -t14 + (t167 + t172) * t88 + (-t171 + (t18 + t184) * t89) * t85, t103 * t89 + t225 * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105 * t89 + (t130 * t86 + t199 * t75) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, 0, 0, t25, 0, 0, 0, 0, 0, 0, -t124, t40, qJD(5) * t129 + t132, t107, 0, 0, 0, 0, 0, 0 (-t188 * t48 - t17) * t89 + t117 * t86 (-t186 * t48 - t18) * t89 + (qJD(5) * t38 + t122) * t86 (t213 - t214) * t185 + t106 * t86, t103 * t86 - t225 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198 * t182, t198 * t183, 0, -t130 * t185 + (t105 - 0.2e1 * t168) * t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78 + (-t199 * t200 - t83) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, 0, -t34, -t58, t6, t5, 0, 0, t182 * t38 + t219, -t17 * t85 + t18 * t88 + (-t213 - t214) * qJD(6), t123, t183 * t37 - t220, -t122, 0, -pkin(5) * t17 - pkin(10) * t123 + t126, -pkin(5) * t18 + pkin(10) * t122 + t127, pkin(10) * t106 + t108, -pkin(5) * t4 + pkin(10) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t187, 0, -t185, 0, -t169, -t168, 0, 0, -t51, -0.4e1 * t161 * t89 + t187 * t201, t66, t51, -t64, 0 (-t176 - t227) * qJD(6) + (t138 * t85 - t175) * qJD(5) (t139 * t85 - t205 * t75) * qJD(6) + (-pkin(10) * t205 + (pkin(5) * t88 + t208) * t86) * qJD(5), t105, -pkin(5) * t169 + pkin(10) * t105; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t185, t187, 0, 0, 0, 0, 0, 0, 0, 0, t64, t66, -qJD(5) * t154 (-pkin(10) * t154 - t223) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t187, -t185, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t67, qJD(5) * t153 (pkin(10) * t153 - t222) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t161, -0.2e1 * t148, 0, -0.2e1 * t161, 0, 0, t85 * t179, t88 * t179, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t17, t34, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, 0, t67, t185, t28, t27, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67, t65, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t66, t64, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t182, 0, -t183, 0, -pkin(10) * t182, pkin(10) * t183, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;