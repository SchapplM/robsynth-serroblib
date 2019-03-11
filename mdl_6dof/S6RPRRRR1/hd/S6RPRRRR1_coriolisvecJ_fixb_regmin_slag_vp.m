% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:55:16
% EndTime: 2019-03-09 06:55:23
% DurationCPUTime: 2.83s
% Computational Cost: add. (5772->301), mult. (13768->405), div. (0->0), fcn. (10330->10), ass. (0->189)
t156 = cos(qJ(6));
t213 = qJD(6) * t156;
t154 = sin(qJ(4));
t155 = sin(qJ(3));
t219 = qJD(1) * t155;
t207 = t154 * t219;
t158 = cos(qJ(4));
t159 = cos(qJ(3));
t218 = qJD(1) * t159;
t208 = t158 * t218;
t119 = -t207 + t208;
t120 = -t154 * t218 - t158 * t219;
t153 = sin(qJ(5));
t157 = cos(qJ(5));
t80 = t157 * t119 + t120 * t153;
t267 = t156 * t80;
t270 = t213 - t267;
t138 = -cos(pkin(11)) * pkin(1) - pkin(2);
t126 = -pkin(3) * t159 + t138;
t121 = t126 * qJD(1);
t90 = -t119 * pkin(4) + t121;
t247 = t80 * t90;
t216 = qJD(5) * t153;
t137 = sin(pkin(11)) * pkin(1) + pkin(7);
t245 = pkin(8) + t137;
t194 = t245 * qJD(1);
t104 = t155 * qJD(2) + t194 * t159;
t100 = t158 * t104;
t103 = t159 * qJD(2) - t194 * t155;
t238 = qJD(3) * pkin(3);
t101 = t103 + t238;
t178 = -t101 * t154 - t100;
t96 = t103 * qJD(3);
t97 = t104 * qJD(3);
t197 = -t154 * t96 - t158 * t97;
t166 = t178 * qJD(4) + t197;
t147 = qJD(3) + qJD(4);
t211 = qJD(1) * qJD(3);
t204 = t159 * t211;
t85 = qJD(4) * t208 - t147 * t207 + t158 * t204;
t31 = -t85 * pkin(9) + t166;
t254 = pkin(9) * t119;
t58 = -t178 + t254;
t199 = t153 * t31 - t58 * t216;
t125 = t154 * t159 + t155 * t158;
t261 = qJD(1) * t125;
t162 = t147 * t261;
t217 = qJD(4) * t154;
t196 = -t104 * t217 - t154 * t97;
t262 = (qJD(4) * t101 + t96) * t158;
t30 = -pkin(9) * t162 + t196 + t262;
t117 = t120 * pkin(9);
t98 = t154 * t104;
t192 = t158 * t101 - t98;
t57 = t117 + t192;
t54 = pkin(4) * t147 + t57;
t3 = (qJD(5) * t54 + t30) * t157 + t199;
t269 = -t247 - t3;
t221 = -qJD(6) + t80;
t268 = qJD(6) + t221;
t152 = sin(qJ(6));
t145 = qJD(5) + t147;
t177 = t119 * t153 - t157 * t120;
t214 = qJD(6) * t152;
t215 = qJD(5) * t157;
t41 = t119 * t215 + t120 * t216 - t153 * t162 + t157 * t85;
t19 = t145 * t213 + t156 * t41 - t177 * t214;
t17 = t19 * t152;
t67 = t145 * t152 + t156 * t177;
t9 = t270 * t67 + t17;
t42 = t177 * qJD(5) + t153 * t85 + t157 * t162;
t72 = t221 * t213;
t241 = t152 * t42 - t72;
t8 = -t177 * t67 + t221 * t267 + t241;
t266 = t152 * t221;
t40 = t156 * t42;
t65 = -t156 * t145 + t152 * t177;
t7 = t177 * t65 - t221 * t266 + t40;
t20 = qJD(6) * t67 + t152 * t41;
t1 = -t152 * t20 + t19 * t156 + t266 * t67 - t270 * t65;
t235 = t153 * t58;
t26 = t157 * t54 - t235;
t23 = -pkin(5) * t145 - t26;
t253 = t23 * t80;
t246 = t177 * t80;
t248 = t177 * t90;
t200 = t153 * t30 - t157 * t31;
t232 = t157 * t58;
t27 = t153 * t54 + t232;
t4 = t27 * qJD(5) + t200;
t265 = -t4 - t248;
t34 = t177 ^ 2 - t80 ^ 2;
t55 = pkin(5) * t177 - pkin(10) * t80;
t28 = -t145 * t80 + t41;
t249 = t221 * t177;
t24 = pkin(10) * t145 + t27;
t43 = -pkin(5) * t80 - pkin(10) * t177 + t90;
t180 = t152 * t24 - t156 * t43;
t205 = t177 * t180 + t23 * t214;
t12 = t152 * t43 + t156 * t24;
t184 = t12 * t177 + t4 * t152 + t23 * t213;
t29 = t145 * t177 - t42;
t124 = t154 * t155 - t158 * t159;
t89 = -t124 * t153 + t125 * t157;
t251 = t42 * t89;
t88 = t157 * t124 + t125 * t153;
t93 = t147 * t124;
t171 = t125 * qJD(4);
t94 = t125 * qJD(3) + t171;
t50 = -t88 * qJD(5) - t153 * t94 - t157 * t93;
t181 = -t221 * t50 + t251;
t209 = t89 * t214;
t258 = -t181 * t156 - t209 * t221;
t123 = t245 * t159;
t122 = t245 * t155;
t226 = t122 * t158;
t68 = -pkin(9) * t125 - t123 * t154 - t226;
t176 = t122 * t154 - t123 * t158;
t69 = -pkin(9) * t124 - t176;
t47 = t153 * t68 + t157 * t69;
t195 = qJD(3) * t245;
t115 = t155 * t195;
t116 = t159 * t195;
t172 = -qJD(4) * t226 - t158 * t115 - t154 * t116 - t123 * t217;
t44 = -pkin(9) * t94 + t172;
t164 = t176 * qJD(4) + t154 * t115 - t158 * t116;
t45 = t93 * pkin(9) + t164;
t46 = t153 * t69 - t157 * t68;
t5 = -t46 * qJD(5) + t153 * t45 + t157 * t44;
t102 = pkin(4) * t124 + t126;
t53 = pkin(5) * t88 - pkin(10) * t89 + t102;
t257 = (qJD(6) * t53 + t5) * t221 - (qJD(6) * t43 + t3) * t88 + t23 * t50 + t4 * t89 - t47 * t42;
t255 = pkin(4) * t120;
t252 = t23 * t89;
t250 = t53 * t42;
t51 = t89 * qJD(5) - t153 * t93 + t157 * t94;
t244 = t19 * t88 + t67 * t51;
t142 = pkin(3) * t158 + pkin(4);
t225 = t153 * t154;
t186 = -t103 * t154 - t100;
t59 = t186 - t254;
t239 = t158 * t103 - t98;
t60 = t117 + t239;
t242 = t153 * t59 + t157 * t60 - t142 * t215 - (-t154 * t216 + (t157 * t158 - t225) * qJD(4)) * pkin(3);
t224 = t154 * t157;
t240 = -t153 * t60 + t157 * t59 + t142 * t216 + (t154 * t215 + (t153 * t158 + t224) * qJD(4)) * pkin(3);
t231 = t50 * t145;
t230 = t93 * t147;
t228 = t120 * t119;
t227 = t121 * t120;
t160 = qJD(3) ^ 2;
t223 = t160 * t155;
t222 = t160 * t159;
t220 = t155 ^ 2 - t159 ^ 2;
t128 = qJD(1) * t138;
t144 = t155 * t238;
t143 = pkin(3) * t219;
t82 = pkin(4) * t94 + t144;
t206 = -pkin(4) * t145 - t54;
t203 = -pkin(3) * t147 - t101;
t114 = pkin(3) * t224 + t142 * t153 + pkin(10);
t52 = -t255 + t55;
t188 = qJD(6) * t114 + t143 + t52;
t140 = pkin(4) * t153 + pkin(10);
t187 = qJD(6) * t140 + t52;
t32 = t153 * t57 + t232;
t183 = pkin(4) * t216 - t32;
t182 = -t20 * t88 - t51 * t65;
t175 = 0.2e1 * qJD(3) * t128;
t174 = -t121 * t119 - t196;
t170 = -t114 * t42 - t221 * t242 - t253;
t33 = t157 * t57 - t235;
t165 = -t140 * t42 - t253 - (-pkin(4) * t215 + t33) * t221;
t163 = -t181 * t152 + t89 * t72;
t71 = pkin(4) * t162 + qJD(3) * t143;
t161 = qJD(1) ^ 2;
t141 = -pkin(4) * t157 - pkin(5);
t113 = pkin(3) * t225 - t142 * t157 - pkin(5);
t105 = t143 - t255;
t86 = t94 * t147;
t70 = -t119 ^ 2 + t120 ^ 2;
t62 = -t120 * t147 - t162;
t61 = -t119 * t147 + t85;
t48 = t51 * t145;
t14 = pkin(5) * t51 - pkin(10) * t50 + t82;
t13 = t42 * pkin(5) - t41 * pkin(10) + t71;
t10 = t156 * t13;
t6 = t47 * qJD(5) + t153 * t44 - t157 * t45;
t2 = [0, 0, 0, 0, 0.2e1 * t155 * t204, -0.2e1 * t220 * t211, t222, -t223, 0, -t137 * t222 + t155 * t175, t137 * t223 + t159 * t175, t120 * t93 + t125 * t85, -t93 * t119 + t120 * t94 - t85 * t124 - t125 * t162, -t230, -t86, 0, -t119 * t144 + t121 * t94 + t164 * t147 + (t126 * t171 + (t155 * pkin(3) * t124 + t125 * t126) * qJD(3)) * qJD(1), t126 * t85 - t121 * t93 - t172 * t147 + (-t120 + t261) * t144, t177 * t50 + t41 * t89, -t177 * t51 - t41 * t88 + t50 * t80 - t251, t231, -t48, 0, t102 * t42 - t6 * t145 + t90 * t51 + t71 * t88 - t80 * t82, t102 * t41 - t145 * t5 + t177 * t82 + t50 * t90 + t71 * t89, -t67 * t209 + (t19 * t89 + t50 * t67) * t156 (-t152 * t67 - t156 * t65) * t50 + (-t17 - t156 * t20 + (t152 * t65 - t156 * t67) * qJD(6)) * t89, t244 - t258, t163 + t182, -t221 * t51 + t42 * t88, t10 * t88 - t180 * t51 + t46 * t20 + t6 * t65 + (-t14 * t221 + t250 + (t221 * t47 - t24 * t88 + t252) * qJD(6)) * t156 + t257 * t152, -t12 * t51 + t46 * t19 + t6 * t67 + ((-qJD(6) * t47 + t14) * t221 - t250 - (-qJD(6) * t24 + t13) * t88 - qJD(6) * t252) * t152 + t257 * t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t223, -t222, 0, 0, 0, 0, 0, -t86, t230, 0, 0, 0, 0, 0, -t48, -t231, 0, 0, 0, 0, 0, t163 - t182, t244 + t258; 0, 0, 0, 0, -t155 * t161 * t159, t220 * t161, 0, 0, 0, -t128 * t219, -t128 * t218, t228, t70, t61, t62, 0, t119 * t143 + t227 - t186 * t147 + (t154 * t203 - t100) * qJD(4) + t197, t120 * t143 + t239 * t147 + (qJD(4) * t203 - t96) * t158 + t174, -t246, t34, t28, t29, 0, t105 * t80 - t240 * t145 + t265, -t105 * t177 + t242 * t145 + t269, t9, t1, t8, t7, t249, t113 * t20 + t240 * t65 + (t188 * t221 - t4) * t156 + t170 * t152 + t205, t113 * t19 + t170 * t156 - t188 * t266 + t240 * t67 + t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t228, t70, t61, t62, 0, -t147 * t178 + t166 + t227, t147 * t192 + t174 - t262, -t246, t34, t28, t29, 0, -t80 * t255 + t145 * t32 - t248 + (t153 * t206 - t232) * qJD(5) - t200, t177 * t255 + t145 * t33 - t247 + (qJD(5) * t206 - t30) * t157 - t199, t9, t1, t8, t7, t249, t141 * t20 + t183 * t65 + (t187 * t221 - t4) * t156 + t165 * t152 + t205, t141 * t19 + t165 * t156 + t183 * t67 - t187 * t266 + t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t246, t34, t28, t29, 0, t145 * t27 + t265, t145 * t26 + t269, t9, t1, t8, t7, t249, -pkin(5) * t20 - t4 * t156 + (-t152 * t26 + t156 * t55) * t221 - t27 * t65 - t152 * t253 - t241 * pkin(10) + t205, -pkin(5) * t19 - (t152 * t55 + t156 * t26) * t221 - t27 * t67 - t23 * t267 + (-t214 * t221 - t40) * pkin(10) + t184; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67 * t65, -t65 ^ 2 + t67 ^ 2, -t221 * t65 + t19, -t221 * t67 - t20, t42, -t268 * t12 - t152 * t3 - t23 * t67 + t10, -t152 * t13 - t156 * t3 + t268 * t180 + t23 * t65;];
tauc_reg  = t2;
