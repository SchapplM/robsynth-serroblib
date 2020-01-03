% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR15_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR15_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:11
% EndTime: 2019-12-31 20:43:19
% DurationCPUTime: 2.83s
% Computational Cost: add. (2282->364), mult. (4831->491), div. (0->0), fcn. (3167->10), ass. (0->210)
t144 = sin(qJ(2));
t211 = qJD(1) * qJD(2);
t197 = t144 * t211;
t148 = cos(qJ(2));
t209 = t148 * qJDD(1);
t274 = t197 - t209;
t261 = pkin(3) + pkin(6);
t142 = sin(qJ(5));
t143 = sin(qJ(4));
t146 = cos(qJ(5));
t147 = cos(qJ(4));
t92 = t142 * t147 + t143 * t146;
t168 = t92 * t144;
t264 = qJD(4) + qJD(5);
t254 = -qJD(1) * t168 - t264 * t92;
t222 = qJD(2) * t143;
t223 = qJD(1) * t148;
t89 = t147 * t223 + t222;
t199 = t143 * t223;
t220 = qJD(2) * t147;
t91 = -t199 + t220;
t181 = t142 * t89 - t146 * t91;
t37 = t142 * t91 + t146 * t89;
t273 = t181 * t37;
t234 = t144 * qJ(3);
t184 = pkin(2) * t148 + t234;
t175 = pkin(1) + t184;
t74 = t175 * qJD(1);
t272 = qJDD(1) * t175;
t271 = t181 ^ 2 - t37 ^ 2;
t224 = qJD(1) * t144;
t120 = qJD(4) + t224;
t113 = qJD(5) + t120;
t213 = qJD(5) * t146;
t214 = qJD(5) * t142;
t31 = -qJD(4) * t89 + t147 * qJDD(2) + t274 * t143;
t32 = -qJD(4) * t199 + t143 * qJDD(2) + (qJD(2) * qJD(4) - t274) * t147;
t6 = -t142 * t32 + t146 * t31 - t89 * t213 - t91 * t214;
t270 = t113 * t37 + t6;
t150 = -pkin(2) - pkin(7);
t84 = t150 * t148 - pkin(1) - t234;
t51 = t84 * qJD(1);
t125 = pkin(6) * t224;
t265 = qJD(3) + t125;
t226 = pkin(3) * t224 + t265;
t57 = t150 * qJD(2) + t226;
t21 = t143 * t57 + t147 * t51;
t14 = -pkin(8) * t89 + t21;
t12 = t14 * t214;
t141 = qJ(4) + qJ(5);
t131 = sin(t141);
t135 = g(3) * t148;
t138 = qJD(2) * qJ(3);
t126 = pkin(6) * t223;
t99 = pkin(3) * t223 + t126;
t73 = t138 + t99;
t43 = pkin(4) * t89 + t73;
t132 = cos(t141);
t145 = sin(qJ(1));
t149 = cos(qJ(1));
t231 = t144 * t149;
t59 = t131 * t231 + t132 * t145;
t233 = t144 * t145;
t61 = -t131 * t233 + t132 * t149;
t269 = g(1) * t59 - g(2) * t61 - t131 * t135 + t37 * t43 + t12;
t119 = pkin(2) * t197;
t241 = qJ(3) * t148;
t182 = pkin(7) * t144 - t241;
t212 = t144 * qJD(3);
t160 = qJD(2) * t182 - t212;
t19 = qJD(1) * t160 + qJDD(1) * t84 + t119;
t196 = t148 * t211;
t210 = t144 * qJDD(1);
t167 = t196 + t210;
t118 = pkin(6) * t196;
t122 = pkin(6) * t210;
t208 = qJDD(3) + t122;
t195 = t118 + t208;
t34 = t167 * pkin(3) + t150 * qJDD(2) + t195;
t192 = -t143 * t19 + t147 * t34;
t158 = -qJD(4) * t21 + t192;
t88 = qJDD(4) + t167;
t2 = pkin(4) * t88 - pkin(8) * t31 + t158;
t217 = qJD(4) * t147;
t207 = -t143 * t34 - t147 * t19 - t57 * t217;
t218 = qJD(4) * t143;
t3 = -pkin(8) * t32 - t51 * t218 - t207;
t203 = -t142 * t3 + t146 * t2;
t20 = -t143 * t51 + t147 * t57;
t13 = -pkin(8) * t91 + t20;
t11 = pkin(4) * t120 + t13;
t248 = t14 * t146;
t5 = t11 * t142 + t248;
t58 = -t131 * t145 + t132 * t231;
t60 = t131 * t149 + t132 * t233;
t268 = -g(1) * t58 - g(2) * t60 - qJD(5) * t5 + t132 * t135 + t43 * t181 + t203;
t157 = qJD(5) * t181 - t142 * t31 - t146 * t32;
t267 = -t113 * t181 + t157;
t69 = t147 * t88;
t266 = -t120 * t218 + t69;
t185 = g(1) * t145 - g(2) * t149;
t186 = g(1) * t149 + g(2) * t145;
t263 = t120 * t73 + t150 * t88;
t262 = t264 * t148;
t256 = g(3) * t144;
t255 = pkin(8) - t150;
t202 = t147 * t224;
t229 = t146 * t147;
t236 = t142 * t143;
t253 = -t142 * t218 - t143 * t214 + t146 * t202 - t224 * t236 + t264 * t229;
t129 = pkin(2) * t224;
t63 = qJD(1) * t182 + t129;
t252 = t143 * t99 + t147 * t63;
t107 = t261 * t144;
t94 = t143 * t107;
t251 = t147 * t84 + t94;
t250 = t120 * t89;
t249 = t120 * t91;
t247 = t143 * t88;
t246 = t147 * t91;
t244 = t31 * t147;
t204 = -pkin(4) * t147 - pkin(3);
t243 = pkin(4) * t217 - t204 * t224 + t265;
t242 = pkin(6) * qJDD(2);
t240 = qJD(2) * t89;
t239 = qJD(2) * t91;
t238 = qJD(4) * t51;
t237 = qJDD(2) * pkin(2);
t235 = t143 * t144;
t232 = t144 * t147;
t230 = t145 * t147;
t228 = t147 * t148;
t227 = t147 * t149;
t108 = t261 * t148;
t139 = t144 ^ 2;
t140 = t148 ^ 2;
t225 = t139 - t140;
t221 = qJD(2) * t144;
t219 = qJD(2) * t148;
t216 = qJD(4) * t148;
t215 = qJD(4) * t150;
t152 = qJD(1) ^ 2;
t206 = t144 * t152 * t148;
t123 = pkin(6) * t209;
t136 = qJDD(2) * qJ(3);
t137 = qJD(2) * qJD(3);
t205 = t123 + t136 + t137;
t200 = t143 * t216;
t198 = pkin(8) * t148 - t84;
t103 = t255 * t147;
t194 = qJD(5) * t11 + t3;
t100 = t261 * t219;
t128 = pkin(2) * t221;
t48 = t128 + t160;
t191 = t147 * t100 - t143 * t48;
t190 = -qJD(2) * pkin(2) + qJD(3);
t98 = t261 * t221;
t178 = -t229 + t236;
t83 = qJDD(5) + t88;
t189 = t254 * t113 - t178 * t83;
t102 = t255 * t143;
t176 = pkin(4) * t148 - pkin(8) * t235;
t81 = t147 * t99;
t188 = qJD(1) * t176 - qJD(5) * t102 - t143 * t63 - t255 * t218 + t81;
t187 = pkin(8) * t202 + t264 * t103 + t252;
t183 = pkin(2) * t144 - t241;
t101 = t125 + t190;
t106 = -t126 - t138;
t179 = t101 * t148 + t106 * t144;
t177 = t120 * t143;
t173 = -t253 * t113 - t92 * t83;
t172 = -0.2e1 * pkin(1) * t211 - t242;
t171 = -t120 * t217 - t247;
t170 = t143 * t100 + t107 * t217 + t147 * t48 - t84 * t218;
t169 = -qJ(3) * t219 - t212;
t166 = pkin(1) * t152 + t186;
t151 = qJD(2) ^ 2;
t165 = pkin(6) * t151 - t185;
t164 = 0.2e1 * t74 * qJD(2) + t242;
t162 = -t148 * t186 - t256;
t161 = 0.2e1 * qJDD(1) * pkin(1) - t165;
t35 = pkin(3) * t209 - qJD(1) * t98 + t205;
t159 = t35 + t162;
t156 = -t144 * t186 - t74 * t224 + t135 + t208;
t33 = qJD(1) * t169 + t119 - t272;
t68 = t128 + t169;
t155 = qJD(1) * t68 + t165 - t272 + t33;
t54 = pkin(6) * t197 - t205;
t62 = t195 - t237;
t153 = qJD(2) * t179 + t62 * t144 - t54 * t148 - t186;
t121 = pkin(4) * t143 + qJ(3);
t96 = -qJ(3) * t223 + t129;
t95 = t147 * t107;
t78 = -t143 * t233 + t227;
t77 = t143 * t149 + t144 * t230;
t76 = t143 * t231 + t230;
t75 = -t143 * t145 + t144 * t227;
t72 = pkin(4) * t228 + t108;
t65 = t92 * t148;
t64 = t178 * t148;
t44 = -pkin(4) * t200 + (-pkin(6) + t204) * t221;
t28 = -pkin(8) * t228 + t251;
t23 = t144 * pkin(4) + t198 * t143 + t95;
t16 = -t178 * t221 + t92 * t262;
t15 = qJD(2) * t168 + t178 * t262;
t10 = pkin(4) * t32 + t35;
t9 = (t144 * t220 + t200) * pkin(8) + t170;
t8 = t176 * qJD(2) + (t198 * t147 - t94) * qJD(4) + t191;
t4 = t11 * t146 - t14 * t142;
t1 = [qJDD(1), t185, t186, qJDD(1) * t139 + 0.2e1 * t144 * t196, 0.2e1 * t144 * t209 - 0.2e1 * t225 * t211, qJDD(2) * t144 + t148 * t151, qJDD(2) * t148 - t144 * t151, 0, t144 * t172 + t148 * t161, -t144 * t161 + t148 * t172, (t139 + t140) * qJDD(1) * pkin(6) + t153, t144 * t164 + t148 * t155, -t144 * t155 + t148 * t164, pkin(6) * t153 - t74 * t68 + (t185 - t33) * t175, -t216 * t246 + (-t148 * t31 + t91 * t221) * t143, (-t143 * t89 + t246) * t221 + (t143 * t32 - t244 + (t143 * t91 + t147 * t89) * qJD(4)) * t148, (t120 * t222 + t31) * t144 + (t171 + t239) * t148, (t120 * t220 - t32) * t144 + (-t240 - t266) * t148, t120 * t219 + t144 * t88, t191 * t120 + (-t143 * t84 + t95) * t88 + t192 * t144 - t98 * t89 + t108 * t32 + t35 * t228 - g(1) * t78 - g(2) * t76 + (t148 * t20 - t232 * t73) * qJD(2) + (-t73 * t143 * t148 - t251 * t120 - t21 * t144) * qJD(4), -t170 * t120 - t251 * t88 - t98 * t91 + t108 * t31 + g(1) * t77 - g(2) * t75 + ((qJD(2) * t73 + t238) * t143 + t207) * t144 + (-qJD(2) * t21 - t35 * t143 - t217 * t73) * t148, -t15 * t181 - t6 * t65, -t15 * t37 - t157 * t65 - t16 * t181 + t6 * t64, t113 * t15 + t144 * t6 - t181 * t219 - t65 * t83, t113 * t16 + t144 * t157 - t219 * t37 + t64 * t83, t113 * t219 + t144 * t83, (-t142 * t9 + t146 * t8) * t113 + (-t142 * t28 + t146 * t23) * t83 + t203 * t144 + t4 * t219 + t44 * t37 - t72 * t157 - t10 * t64 - t43 * t16 - g(1) * t61 - g(2) * t59 + ((-t142 * t23 - t146 * t28) * t113 - t5 * t144) * qJD(5), -t5 * t219 + g(1) * t60 - g(2) * t58 - t10 * t65 + t12 * t144 + t43 * t15 - t44 * t181 + t72 * t6 + (-(-qJD(5) * t28 + t8) * t113 - t23 * t83 - t2 * t144) * t142 + (-(qJD(5) * t23 + t9) * t113 - t28 * t83 - t194 * t144) * t146; 0, 0, 0, -t206, t225 * t152, t210, t209, qJDD(2), t144 * t166 - t122 - t135, t148 * t166 - t123 + t256, -t183 * qJDD(1) + ((-t106 - t138) * t144 + (-t101 + t190) * t148) * qJD(1), -t96 * t223 + t156 - 0.2e1 * t237, t123 + 0.2e1 * t136 + 0.2e1 * t137 + (qJD(1) * t96 - g(3)) * t144 + (-qJD(1) * t74 - t186) * t148, -pkin(6) * qJD(1) * t179 - t62 * pkin(2) - g(3) * t184 - t54 * qJ(3) - t106 * qJD(3) + t186 * t183 + t74 * t96, -t177 * t91 + t244, (-t32 - t249) * t147 + (-t31 + t250) * t143, (-t120 * t235 - t148 * t91) * qJD(1) + t266, (-t120 * t232 + t148 * t89) * qJD(1) + t171, -t120 * t223, -t20 * t223 + qJ(3) * t32 - t81 * t120 + t226 * t89 + t263 * t147 + ((t63 - t215) * t120 + t159) * t143, qJ(3) * t31 + t252 * t120 + t21 * t223 + t226 * t91 - t263 * t143 + (-t120 * t215 + t159) * t147, -t178 * t6 - t181 * t254, -t157 * t178 + t181 * t253 - t254 * t37 - t6 * t92, t181 * t223 + t189, t223 * t37 + t173, -t113 * t223, (t102 * t142 - t103 * t146) * t83 - t121 * t157 + t10 * t92 - t4 * t223 + t253 * t43 + t243 * t37 + (t142 * t187 - t146 * t188) * t113 + t162 * t131, -(-t102 * t146 - t103 * t142) * t83 + t121 * t6 - t10 * t178 + t5 * t223 + t254 * t43 - t243 * t181 + (t142 * t188 + t146 * t187) * t113 + t162 * t132; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t210, qJDD(2) + t206, -t139 * t152 - t151, qJD(2) * t106 + t118 + t156 - t237, 0, 0, 0, 0, 0, -t120 * t177 - t240 + t69, -t120 ^ 2 * t147 - t239 - t247, 0, 0, 0, 0, 0, -qJD(2) * t37 + t189, qJD(2) * t181 + t173; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91 * t89, -t89 ^ 2 + t91 ^ 2, t31 + t250, t249 - t32, t88, -g(1) * t75 - g(2) * t77 + g(3) * t228 + t120 * t21 - t73 * t91 + t158, g(1) * t76 - g(2) * t78 + t120 * t20 + t73 * t89 + (t238 - t135) * t143 + t207, -t273, t271, t270, t267, t83, -(-t13 * t142 - t248) * t113 + (-t113 * t214 + t146 * t83 - t37 * t91) * pkin(4) + t268, (-t113 * t14 - t2) * t142 + (t113 * t13 - t194) * t146 + (-t113 * t213 - t142 * t83 + t181 * t91) * pkin(4) + t269; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t273, t271, t270, t267, t83, t113 * t5 + t268, t113 * t4 - t142 * t2 - t146 * t194 + t269;];
tau_reg = t1;
