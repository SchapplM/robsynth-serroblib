% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP8_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP8_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:09:27
% EndTime: 2019-12-31 21:09:35
% DurationCPUTime: 3.08s
% Computational Cost: add. (2778->449), mult. (6029->537), div. (0->0), fcn. (3767->6), ass. (0->210)
t263 = pkin(3) + qJ(5);
t142 = cos(qJ(2));
t126 = t142 * qJDD(1);
t139 = sin(qJ(2));
t215 = qJD(1) * qJD(2);
t282 = -t139 * t215 + t126;
t71 = qJDD(3) - t282;
t205 = t263 * t71;
t143 = cos(qJ(1));
t236 = t139 * t143;
t140 = sin(qJ(1));
t238 = t139 * t140;
t287 = g(1) * t236 + g(2) * t238;
t138 = sin(qJ(3));
t141 = cos(qJ(3));
t128 = t139 * pkin(7);
t132 = t142 * pkin(2);
t208 = -pkin(1) - t132;
t173 = t208 - t128;
t65 = t173 * qJD(1);
t223 = qJD(1) * t142;
t121 = pkin(6) * t223;
t89 = qJD(2) * pkin(7) + t121;
t32 = t138 * t89 - t141 * t65;
t222 = qJD(2) * t138;
t224 = qJD(1) * t139;
t75 = t141 * t224 + t222;
t172 = t75 * pkin(4) + t32;
t229 = qJD(4) + t172;
t218 = qJD(3) * t139;
t286 = qJD(1) * t218 - qJDD(2);
t214 = t139 * qJDD(1);
t26 = ((qJD(3) + t223) * qJD(2) + t214) * t138 + t286 * t141;
t63 = t71 * qJ(4);
t104 = -qJD(3) + t223;
t92 = qJD(4) * t104;
t285 = t92 - t63;
t100 = t104 ^ 2;
t70 = t75 ^ 2;
t284 = -t70 - t100;
t239 = t138 * t142;
t213 = pkin(3) * t239;
t219 = qJD(3) * t138;
t283 = pkin(3) * t219 - qJD(1) * t213 - t138 * qJD(4) - t121;
t179 = g(1) * t143 + g(2) * t140;
t281 = -pkin(4) * t26 + qJDD(5);
t279 = -0.2e1 * pkin(1);
t216 = t141 * qJD(2);
t73 = t138 * t224 - t216;
t278 = t73 ^ 2;
t277 = 0.2e1 * t63;
t276 = pkin(4) + pkin(7);
t275 = pkin(3) * t71;
t273 = pkin(4) * t73;
t272 = pkin(7) * t71;
t199 = t142 * t215;
t25 = -qJD(3) * t216 + (-t199 - t214) * t141 + t286 * t138;
t271 = t25 * pkin(4);
t270 = g(1) * t140;
t267 = g(2) * t143;
t266 = g(3) * t139;
t265 = g(3) * t142;
t264 = t75 * t73;
t242 = qJ(4) * t141;
t174 = qJ(5) * t138 - t242;
t162 = t174 * t142;
t262 = -qJD(1) * t162 + qJD(3) * t174 - t141 * qJD(5) + t283;
t217 = qJD(3) * t141;
t261 = -qJ(4) * t217 + t223 * t242 + t283;
t180 = pkin(2) * t139 - pkin(7) * t142;
t81 = t180 * qJD(2);
t226 = -t128 - t132;
t83 = -pkin(1) + t226;
t260 = t138 * t81 + t83 * t217;
t33 = t138 * t65 + t141 * t89;
t196 = pkin(6) * t141 - qJ(4);
t155 = -pkin(4) * t239 - t196 * t139;
t78 = t180 * qJD(1);
t61 = t138 * t78;
t259 = -qJD(1) * t155 - t276 * t219 - t61;
t207 = -pkin(6) * t138 - pkin(3);
t234 = t141 * t142;
t153 = pkin(4) * t234 + (-qJ(5) + t207) * t139;
t248 = t141 * t78;
t91 = t276 * t141;
t258 = -qJD(1) * t153 + qJD(3) * t91 + t248;
t257 = t287 * t138;
t256 = t287 * t141;
t255 = pkin(7) * qJD(3);
t254 = qJ(4) * t26;
t22 = qJ(4) * t104 - t33;
t253 = t104 * t22;
t252 = t104 * t33;
t251 = t104 * t73;
t250 = t104 * t75;
t249 = t141 * t75;
t247 = t142 * t75;
t246 = t25 * t138;
t245 = t73 * qJ(4);
t109 = pkin(6) * t234;
t244 = t138 * t83 + t109;
t241 = qJD(2) * t75;
t240 = t138 * t139;
t237 = t139 * t141;
t235 = t140 * t142;
t233 = t142 * qJ(5);
t232 = t142 * t143;
t231 = t143 * t138;
t230 = -qJD(4) - t32;
t18 = t33 - t273;
t228 = -qJD(5) - t18;
t227 = pkin(3) * t240 - qJ(4) * t237;
t135 = t139 ^ 2;
t225 = -t142 ^ 2 + t135;
t221 = qJD(2) * t139;
t220 = qJD(2) * t142;
t107 = pkin(6) * t239;
t88 = -qJD(2) * pkin(2) + pkin(6) * t224;
t167 = -t75 * qJ(4) + t88;
t14 = t263 * t73 + t167;
t210 = t14 * t219;
t209 = t14 * t217;
t119 = pkin(6) * t214;
t56 = -qJDD(2) * pkin(2) + pkin(6) * t199 + t119;
t152 = t25 * qJ(4) - t75 * qJD(4) + t56;
t2 = t73 * qJD(5) + t263 * t26 + t152;
t206 = -t2 - t265;
t204 = t104 * t222;
t203 = t104 * t216;
t202 = t104 * t219;
t201 = t138 * t218;
t197 = -t138 * qJ(4) - pkin(2);
t57 = t138 * t235 + t141 * t143;
t58 = t140 * t234 - t231;
t195 = -t57 * pkin(3) + qJ(4) * t58;
t59 = -t140 * t141 + t142 * t231;
t60 = t138 * t140 + t141 * t232;
t194 = -t59 * pkin(3) + qJ(4) * t60;
t37 = qJD(1) * t81 + qJDD(1) * t173;
t55 = t282 * pkin(6) + qJDD(2) * pkin(7);
t193 = t138 * t55 - t141 * t37 + t89 * t217 + t65 * t219;
t192 = -t138 * t37 - t141 * t55 - t65 * t217 + t89 * t219;
t191 = t141 * t83 - t107;
t190 = t139 * pkin(3) * t217 + pkin(6) * t220 + qJ(4) * t201 + qJD(2) * t213;
t187 = pkin(3) * t234 + qJ(4) * t239 - t226;
t186 = -qJD(3) * t109 + t141 * t81 - t83 * t219;
t185 = g(1) * t57 - g(2) * t59;
t184 = g(1) * t58 - g(2) * t60;
t183 = -t58 * pkin(3) + t143 * pkin(6) - t57 * qJ(4);
t40 = qJ(4) * t142 - t244;
t182 = -t142 * qJD(4) + t260;
t181 = t139 * t207;
t178 = -qJDD(4) - t193;
t177 = qJD(3) * t88 - t272;
t10 = t263 * t104 + t229;
t13 = qJD(5) - t22 - t273;
t176 = t10 * t141 - t13 * t138;
t21 = pkin(3) * t104 - t230;
t175 = t138 * t22 + t141 * t21;
t4 = t192 + t285;
t171 = pkin(3) * t141 - t197;
t170 = -qJ(5) * t240 - t227;
t169 = t104 * t255 - t265;
t166 = -pkin(6) * qJDD(2) + t215 * t279;
t165 = -t104 * t217 + t138 * t71;
t164 = t141 * t71 + t202;
t5 = pkin(3) * t26 + t152;
t163 = t169 - t5;
t145 = qJD(1) ^ 2;
t161 = pkin(1) * t145 + t179;
t160 = t263 * t141 - t197;
t144 = qJD(2) ^ 2;
t159 = pkin(6) * t144 + qJDD(1) * t279 + t267;
t158 = t143 * pkin(1) + pkin(2) * t232 + t60 * pkin(3) + t140 * pkin(6) + pkin(7) * t236 + qJ(4) * t59;
t24 = pkin(3) * t73 + t167;
t157 = t104 * t24 + t272;
t156 = t71 - t264;
t154 = g(1) * t59 + g(2) * t57 + g(3) * t240 - t193;
t151 = -qJDD(4) + t154;
t12 = -t25 - t251;
t150 = g(1) * t60 + g(2) * t58 + g(3) * t237 + t192;
t149 = t24 * t75 - t151;
t148 = t14 * t75 - t151 - t271;
t147 = -t14 * t73 - t150 + t281;
t146 = -t250 - t26;
t131 = t142 * pkin(3);
t129 = t139 * pkin(6);
t115 = g(1) * t238;
t112 = pkin(7) * t232;
t108 = pkin(7) * t235;
t90 = t276 * t138;
t47 = t129 + t227;
t41 = t131 - t191;
t39 = t129 - t170;
t36 = pkin(3) * t75 + t245;
t35 = qJD(1) * t181 - t248;
t34 = t196 * t224 - t61;
t30 = -pkin(4) * t240 - t40;
t27 = t233 + t107 + t131 + (pkin(4) * t139 - t83) * t141;
t19 = t263 * t75 + t245;
t16 = (-qJ(4) * t220 - qJD(4) * t139) * t141 + t190;
t15 = qJD(2) * t181 - t186;
t11 = -qJ(4) * t221 + (t139 * t216 + t142 * t219) * pkin(6) - t182;
t9 = qJD(2) * t162 + (qJD(5) * t138 + (qJ(5) * qJD(3) - qJD(4)) * t141) * t139 + t190;
t8 = (-pkin(4) * t237 - t107) * qJD(3) + t155 * qJD(2) + t182;
t7 = -pkin(4) * t201 + qJD(2) * t153 + t142 * qJD(5) - t186;
t6 = -t178 - t275;
t3 = -t4 + t281;
t1 = qJD(5) * t104 - t178 - t205 - t271;
t17 = [qJDD(1), -t267 + t270, t179, qJDD(1) * t135 + 0.2e1 * t139 * t199, 0.2e1 * t139 * t126 - 0.2e1 * t225 * t215, qJDD(2) * t139 + t142 * t144, qJDD(2) * t142 - t139 * t144, 0, t166 * t139 + (-t159 + t270) * t142, t139 * t159 + t142 * t166 - t115, t216 * t247 + (-t25 * t141 - t75 * t219) * t139, (-t138 * t75 - t141 * t73) * t220 + (t246 - t141 * t26 + (t138 * t73 - t249) * qJD(3)) * t139, (t25 - t203) * t142 + (t164 + t241) * t139, (t26 + t204) * t142 + (-qJD(2) * t73 - t165) * t139, -t104 * t221 - t142 * t71, -t186 * t104 + t191 * t71 + ((pkin(6) * t73 + t138 * t88) * qJD(2) + t193) * t142 + (t88 * t217 - t32 * qJD(2) + t56 * t138 + (t26 - t204) * pkin(6)) * t139 + t184, t260 * t104 - t244 * t71 + (t88 * t216 + (-t202 + t241) * pkin(6) - t192) * t142 + (-t88 * t219 - t33 * qJD(2) + t56 * t141 + (-t25 - t203) * pkin(6)) * t139 - t185, t11 * t73 + t15 * t75 - t41 * t25 + t40 * t26 + t115 + t175 * t220 + (-t267 + t138 * t4 + t141 * t6 + (-t138 * t21 + t141 * t22) * qJD(3)) * t139, -t15 * t104 - t16 * t73 - t47 * t26 + t41 * t71 + (-t24 * t222 - t6) * t142 + (qJD(2) * t21 - t5 * t138 - t24 * t217) * t139 - t184, t11 * t104 - t16 * t75 + t47 * t25 - t40 * t71 + (-t24 * t216 + t4) * t142 + (-qJD(2) * t22 - t5 * t141 + t24 * t219) * t139 + t185, -g(1) * t183 - g(2) * t158 + t22 * t11 + t21 * t15 + t24 * t16 - t173 * t270 + t4 * t40 + t6 * t41 + t5 * t47, -t27 * t25 - t30 * t26 + t7 * t75 - t8 * t73 + t115 + t176 * t220 + (-t267 + t1 * t141 - t138 * t3 + (-t10 * t138 - t13 * t141) * qJD(3)) * t139, -t8 * t104 + t39 * t25 + t30 * t71 - t9 * t75 + (-t14 * t216 - t3) * t142 + (qJD(2) * t13 - t2 * t141 + t210) * t139 + t185, t7 * t104 + t39 * t26 - t27 * t71 + t9 * t73 + (t14 * t222 + t1) * t142 + (-qJD(2) * t10 + t2 * t138 + t209) * t139 + t184, t2 * t39 + t14 * t9 + t1 * t27 + t10 * t7 + t3 * t30 + t13 * t8 - g(1) * (-t58 * qJ(5) + t183) - g(2) * (pkin(4) * t236 + qJ(5) * t60 + t158) - (-t276 * t139 + t208) * t270; 0, 0, 0, -t139 * t145 * t142, t225 * t145, t214, t126, qJDD(2), t139 * t161 - t119 - t265, t266 + (-pkin(6) * qJDD(1) + t161) * t142, -t104 * t249 - t246, (-t25 + t251) * t141 + (-t26 + t250) * t138, (t104 * t234 - t139 * t75) * qJD(1) + t165, (-t104 * t239 + t139 * t73) * qJD(1) + t164, t104 * t224, -pkin(2) * t26 + t177 * t138 + (-t265 - t56 + (t78 + t255) * t104) * t141 + (-t88 * t239 + t32 * t139 + (t104 * t240 - t142 * t73) * pkin(6)) * qJD(1) + t256, pkin(2) * t25 - t61 * t104 + t177 * t141 + (-t169 + t56) * t138 + (-t88 * t234 + t33 * t139 + (t104 * t237 - t247) * pkin(6)) * qJD(1) - t257, -t266 - t34 * t73 - t35 * t75 - t179 * t142 + (-t4 - t104 * t21 + (qJD(3) * t75 - t26) * pkin(7)) * t141 + (t6 - t253 + (qJD(3) * t73 - t25) * pkin(7)) * t138, t35 * t104 + t157 * t138 - t163 * t141 + t171 * t26 - t21 * t224 - t261 * t73 - t256, -t34 * t104 + t163 * t138 + t157 * t141 - t171 * t25 + t22 * t224 - t261 * t75 + t257, -t22 * t34 - t21 * t35 - g(1) * t112 - g(2) * t108 - g(3) * t187 + t261 * t24 + (qJD(3) * t175 + t6 * t138 - t4 * t141) * pkin(7) + (t179 * t139 - t5) * t171, -t266 + t1 * t138 + t141 * t3 - t25 * t90 - t26 * t91 + t258 * t75 - t259 * t73 + t176 * qJD(3) + (-qJD(1) * t176 - t179) * t142, -t209 - t160 * t25 + t91 * t71 - t262 * t75 + t206 * t138 - t259 * t104 + (-t13 * t139 + t14 * t234) * qJD(1) + t257, t210 - t160 * t26 - t90 * t71 + t262 * t73 + t206 * t141 + t258 * t104 + (t10 * t139 - t14 * t239) * qJD(1) + t256, -t2 * t160 + t1 * t90 + t3 * t91 - g(1) * (pkin(4) * t232 + t112) - g(2) * (pkin(4) * t235 + t108) - g(3) * (t141 * t233 + t187) + t262 * t14 + t259 * t13 + t258 * t10 + (-g(3) * pkin(4) + t179 * t160) * t139; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t264, t70 - t278, t12, t146, t71, -t75 * t88 + t154 - t252, t104 * t32 + t73 * t88 + t150, pkin(3) * t25 - t254 + (-t22 - t33) * t75 + (t21 + t230) * t73, t36 * t73 + t149 + t252 - 0.2e1 * t275, t230 * t104 - t24 * t73 + t36 * t75 - t150 + t277 - t92, -t6 * pkin(3) - g(1) * t194 - g(2) * t195 + g(3) * t227 - t4 * qJ(4) - t21 * t33 + t230 * t22 - t24 * t36, -t254 + t263 * t25 + (t13 + t228) * t75 + (t10 - t229) * t73, -t104 * t172 + t19 * t75 + t147 + t277 - 0.2e1 * t92, -t19 * t73 + (-0.2e1 * qJD(5) - t18) * t104 + 0.2e1 * t205 - t148, -t1 * t263 + t3 * qJ(4) - t14 * t19 - g(1) * (-qJ(5) * t59 + t194) - g(2) * (-qJ(5) * t57 + t195) - g(3) * t170 + t229 * t13 + t228 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t156, t284, t149 - t253 - t275, t12, t284, -t156, -t205 + (qJD(5) + t13) * t104 + t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, t71 + t264, -t100 - t278, -t10 * t104 + t147 - t285;];
tau_reg = t17;
