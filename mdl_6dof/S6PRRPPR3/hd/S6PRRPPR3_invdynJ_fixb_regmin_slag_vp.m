% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 21:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRRPPR3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPPR3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPPR3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPPR3_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:11:05
% EndTime: 2019-03-08 21:11:10
% DurationCPUTime: 3.25s
% Computational Cost: add. (1861->399), mult. (4224->512), div. (0->0), fcn. (3127->10), ass. (0->218)
t145 = cos(qJ(3));
t138 = cos(pkin(6));
t241 = qJD(1) * t138;
t100 = t145 * t241;
t142 = sin(qJ(3));
t143 = sin(qJ(2));
t136 = sin(pkin(6));
t242 = qJD(1) * t136;
t77 = qJD(2) * pkin(8) + t143 * t242;
t265 = t142 * t77 - t100;
t290 = qJD(4) + t265;
t240 = qJD(2) * t142;
t107 = qJD(6) + t240;
t289 = t107 - qJD(6);
t288 = t142 * qJ(4) + pkin(2);
t287 = qJDD(3) * qJ(4) + qJD(3) * qJD(4);
t37 = -qJD(3) * pkin(3) + t290;
t230 = qJD(2) * qJD(3);
t210 = t145 * t230;
t226 = t142 * qJDD(2);
t286 = t210 + t226;
t211 = t142 * t230;
t225 = t145 * qJDD(2);
t285 = -t211 + t225;
t147 = -pkin(3) - pkin(4);
t218 = qJD(3) * t147;
t33 = -qJ(5) * t240 + t265;
t247 = -qJD(4) - t33;
t19 = t218 - t247;
t238 = qJD(2) * t145;
t146 = cos(qJ(2));
t237 = qJD(2) * t146;
t216 = t136 * t237;
t255 = t136 * t143;
t220 = t142 * t255;
t31 = -qJD(3) * t220 + (qJD(3) * t138 + t216) * t145;
t254 = t136 * t145;
t66 = t138 * t142 + t143 * t254;
t32 = qJD(3) * t66 + t142 * t216;
t65 = -t138 * t145 + t220;
t284 = qJD(2) * (-t142 * t32 - t145 * t31 + (t142 * t66 - t145 * t65) * qJD(3)) - (t142 * t65 + t145 * t66) * qJDD(2);
t283 = t107 ^ 2;
t277 = pkin(8) - qJ(5);
t132 = qJD(3) * qJ(4);
t42 = t142 * t241 + t145 * t77;
t34 = -qJ(5) * t238 + t42;
t22 = -t132 - t34;
t38 = t132 + t42;
t282 = 0.2e1 * t287;
t214 = qJDD(3) * t147;
t188 = pkin(5) * t142 + pkin(9) * t145;
t139 = qJD(2) * pkin(2);
t99 = t146 * t242;
t78 = -t99 - t139;
t44 = -pkin(3) * t238 - qJ(4) * t240 + t78;
t25 = pkin(4) * t238 + qJD(5) - t44;
t15 = qJD(2) * t188 + t25;
t264 = sin(pkin(10));
t202 = t264 * t146;
t137 = cos(pkin(10));
t252 = t137 * t143;
t60 = t138 * t252 + t202;
t203 = t264 * t143;
t251 = t137 * t146;
t62 = -t138 * t203 + t251;
t189 = g(1) * t62 + g(2) * t60;
t20 = qJD(3) * pkin(5) - t22;
t131 = -pkin(9) + t147;
t235 = qJD(3) * t142;
t212 = qJD(1) * t235;
t228 = qJDD(1) * t138;
t234 = qJD(3) * t145;
t213 = qJD(1) * t237;
t227 = qJDD(1) * t143;
t258 = qJDD(2) * pkin(8);
t47 = t258 + (t213 + t227) * t136;
t201 = t138 * t212 + t142 * t47 - t145 * t228 + t77 * t234;
t187 = -qJDD(4) - t201;
t229 = qJD(2) * qJD(5);
t150 = -qJ(5) * t286 - t142 * t229 - t187;
t3 = qJDD(3) * t131 + t150;
t79 = -t145 * pkin(3) - t288;
t69 = t145 * pkin(4) - t79;
t43 = t188 + t69;
t57 = -t142 * qJD(5) + t277 * t234;
t102 = qJ(5) * t211;
t200 = -qJD(3) * t100 - t142 * t228 - t145 * t47 + t77 * t235;
t12 = -t200 + t287;
t8 = t145 * (qJ(5) * qJDD(2) + t229) - t102 - t12;
t6 = qJDD(3) * pkin(5) - t8;
t70 = qJDD(6) + t286;
t84 = t277 * t142;
t280 = -(qJD(6) * t43 + t57) * t107 + (qJD(3) * t20 - qJD(6) * t15 - t3) * t142 - t6 * t145 - t84 * t70 + t189;
t113 = qJ(4) * t238;
t172 = pkin(5) * t145 + t131 * t142;
t27 = -t137 * t136 * t142 + t145 * t60;
t204 = t136 * t264;
t29 = t142 * t204 + t62 * t145;
t173 = g(1) * t29 + g(2) * t27 + g(3) * t66;
t279 = t107 * (qJD(2) * t172 + qJD(6) * t131 + t113) + t173 - t6;
t253 = t136 * t146;
t239 = qJD(2) * t143;
t217 = t136 * t239;
t89 = qJD(1) * t217;
t276 = t145 * t89 + t212 * t253;
t141 = sin(qJ(6));
t275 = t141 * t70;
t144 = cos(qJ(6));
t274 = t144 * t70;
t273 = t145 * t22;
t59 = t138 * t251 - t203;
t272 = t145 * t59;
t61 = -t138 * t202 - t252;
t271 = t145 * t61;
t231 = t144 * qJD(3);
t71 = t141 * t238 - t231;
t270 = t145 * t71;
t236 = qJD(3) * t141;
t72 = t144 * t238 + t236;
t269 = t145 * t72;
t232 = qJD(6) * t145;
t215 = t141 * t232;
t17 = qJD(2) * t215 - t141 * qJDD(3) + (-qJD(3) * qJD(6) - t285) * t144;
t268 = t17 * t141;
t267 = t71 * t107;
t266 = t72 * t107;
t263 = pkin(8) * qJDD(3);
t262 = qJ(4) * t145;
t261 = qJD(3) * t22;
t260 = qJD(3) * t71;
t259 = qJD(3) * t72;
t135 = qJDD(2) * pkin(2);
t257 = qJDD(3) * pkin(3);
t256 = t107 * t144;
t250 = t142 * t146;
t249 = t144 * t146;
t248 = t145 * t146;
t246 = -qJD(5) - t25;
t121 = t142 * qJD(4);
t245 = qJ(4) * t234 + t121;
t133 = t142 ^ 2;
t134 = t145 ^ 2;
t244 = t133 - t134;
t243 = t133 + t134;
t233 = qJD(6) * t107;
t224 = pkin(3) * t272 + t288 * t59;
t223 = pkin(3) * t271 + t288 * t61;
t96 = qJDD(1) * t253;
t46 = t89 - t96 - t135;
t222 = t107 * t141 * t142;
t221 = t142 * t256;
t149 = qJD(2) ^ 2;
t219 = t142 * t149 * t145;
t209 = t146 * t230;
t208 = t78 - t139;
t26 = t137 * t254 + t142 * t60;
t207 = -t26 * pkin(3) + qJ(4) * t27;
t28 = t142 * t62 - t145 * t204;
t206 = -t28 * pkin(3) + qJ(4) * t29;
t205 = -t65 * pkin(3) + qJ(4) * t66;
t199 = qJD(2) * t69 + t25;
t198 = qJD(2) * t79 + t44;
t196 = t246 * t142;
t193 = 0.2e1 * t210;
t191 = t142 * t218;
t190 = g(1) * t61 + g(2) * t59;
t16 = qJD(3) * t131 - t247;
t186 = t141 * t16 - t144 * t15;
t5 = t141 * t15 + t144 * t16;
t182 = g(3) * (pkin(2) * t253 + pkin(8) * t255 + (pkin(3) * t248 + qJ(4) * t250) * t136);
t106 = g(3) * t255;
t181 = t106 + t189;
t180 = qJDD(2) * t146 - t143 * t149;
t35 = t136 * t249 - t141 * t65;
t36 = t141 * t253 + t144 * t65;
t178 = -t144 * t233 - t275;
t177 = t141 * t233 - t274;
t176 = -t141 * t143 + t142 * t249;
t175 = t141 * t250 + t143 * t144;
t171 = pkin(3) * t225 + t286 * qJ(4) + qJD(2) * t121 - t46;
t169 = -qJD(6) * t16 - t190;
t167 = g(3) * t253 + t190;
t166 = t172 * qJD(3);
t165 = pkin(4) * t225 + qJDD(5) + t171;
t148 = qJD(3) ^ 2;
t164 = pkin(8) * t148 + t167;
t163 = g(1) * t28 + g(2) * t26 + g(3) * t65 - t201;
t162 = -t173 - t200;
t161 = -t131 * t70 + (-t20 + t34) * t107;
t160 = -qJDD(4) + t163;
t159 = -(-qJD(6) * t84 + t166 + t245) * t107 - t43 * t70 + t20 * t232;
t158 = qJD(3) * t42 + t163;
t157 = qJD(3) * t265 + t162;
t156 = t164 + t46 - t135;
t11 = qJD(2) * t191 + t165;
t45 = t191 + t245;
t155 = -qJD(2) * t45 - qJDD(2) * t69 - t11 + t167;
t18 = qJD(6) * t72 - t144 * qJDD(3) + t285 * t141;
t154 = -qJ(5) * t226 - t160;
t14 = pkin(3) * t211 - t171;
t58 = pkin(3) * t235 - t245;
t153 = -qJD(2) * t58 - qJDD(2) * t79 - t14 - t164;
t13 = -t187 - t257;
t152 = t12 * t145 + t13 * t142 + (-t142 * t38 + t145 * t37) * qJD(3) - t189;
t10 = t31 * qJD(3) + t66 * qJDD(3) + (t142 * t180 + t145 * t209) * t136;
t9 = -t32 * qJD(3) - t65 * qJDD(3) + (-t142 * t209 + t145 * t180) * t136;
t140 = qJ(4) + pkin(5);
t101 = -t133 * t149 - t148;
t85 = t277 * t145;
t83 = qJDD(3) + t219;
t75 = t142 * t89;
t73 = pkin(3) * t240 - t113;
t56 = qJD(5) * t145 + t235 * t277;
t52 = t147 * t240 + t113;
t7 = t214 + t150;
t2 = qJD(2) * t166 + qJDD(2) * t188 + t165;
t1 = t144 * t2;
t4 = [qJDD(1) - g(3), 0, t180 * t136 (-qJDD(2) * t143 - t146 * t149) * t136, 0, 0, 0, 0, 0, t9, -t10, t9, -t284, t10, t12 * t66 + t13 * t65 + t38 * t31 + t37 * t32 - g(3) + (-t14 * t146 + t239 * t44) * t136, t10, -t9, t284, t19 * t32 - t22 * t31 + t7 * t65 - t8 * t66 - g(3) + (t11 * t146 - t239 * t25) * t136, 0, 0, 0, 0, 0 (-qJD(6) * t36 - t32 * t141 - t144 * t217) * t107 + t35 * t70 - t31 * t71 - t66 * t18 -(qJD(6) * t35 - t141 * t217 + t32 * t144) * t107 - t36 * t70 - t31 * t72 + t66 * t17; 0, qJDD(2), -t167 + t96, -t136 * t227 + t181, qJDD(2) * t133 + t142 * t193, 0.2e1 * t142 * t225 - 0.2e1 * t230 * t244, qJDD(3) * t142 + t145 * t148, qJDD(3) * t145 - t142 * t148, 0 (qJD(3) * t208 - t263) * t142 - t156 * t145 + t276, -t75 + (-t263 + (t208 + t99) * qJD(3)) * t145 + t156 * t142 (qJD(3) * t198 - t263) * t142 + t153 * t145 + t276, -t106 + t152 + (-t136 * t213 + t258) * t243, t75 + (t263 + (-t198 - t99) * qJD(3)) * t145 + t153 * t142, t14 * t79 + t44 * t58 - g(1) * t223 - g(2) * t224 - t182 + (-t143 * t44 + (-t142 * t37 - t145 * t38) * t146) * t242 + t152 * pkin(8), qJDD(3) * t85 + t75 + (-t56 + (t199 - t99) * t145) * qJD(3) - t155 * t142, qJDD(3) * t84 + (t142 * t199 + t57) * qJD(3) + t155 * t145 - t276 (-qJD(3) * t19 - qJDD(2) * t85 + t8) * t145 + (-qJDD(2) * t84 - t261 - t7) * t142 + (-t142 * t57 + t145 * t56 + (t142 * t85 - t145 * t84) * qJD(3) + t243 * t99) * qJD(2) + t181, t7 * t84 + t19 * t57 - t8 * t85 + t22 * t56 + t11 * t69 + t25 * t45 - g(1) * (pkin(4) * t271 + t277 * t62 + t223) - g(2) * (pkin(4) * t272 + t277 * t60 + t224) - t182 + (-g(3) * (pkin(4) * t248 - qJ(5) * t143) + (t25 * t143 + (-t142 * t19 + t273) * t146) * qJD(1)) * t136, -t72 * t215 + (-t145 * t17 - t235 * t72) * t144 (t141 * t72 + t144 * t71) * t235 + (t268 - t144 * t18 + (t141 * t71 - t144 * t72) * qJD(6)) * t145 (t107 * t231 + t17) * t142 + (t177 - t259) * t145 (-t107 * t236 + t18) * t142 + (-t178 + t260) * t145, t107 * t234 + t142 * t70, -t186 * t234 + t1 * t142 - t85 * t18 + t56 * t71 + (t142 * t169 - t159) * t144 + t280 * t141 + (-g(3) * t176 + (t107 * t175 + t248 * t71) * qJD(1)) * t136, -t5 * t234 + t85 * t17 + t56 * t72 + ((-t169 - t2) * t142 + t159) * t141 + t280 * t144 + (g(3) * t175 + (t107 * t176 + t248 * t72) * qJD(1)) * t136; 0, 0, 0, 0, -t219, t244 * t149, t226, t225, qJDD(3), -t240 * t78 + t158, -t238 * t78 - t157, 0.2e1 * t257 - qJDD(4) + (-t142 * t44 + t145 * t73) * qJD(2) + t158 (-pkin(3) * t142 + t262) * qJDD(2) (t142 * t73 + t145 * t44) * qJD(2) + t157 + t282, -t13 * pkin(3) - g(1) * t206 - g(2) * t207 - g(3) * t205 + t12 * qJ(4) + t290 * t38 - t37 * t42 - t44 * t73, -qJ(5) * t225 + qJD(3) * t33 + t102 + (-t142 * t52 + t145 * t246) * qJD(2) + t162 + t282, -qJD(3) * t34 + 0.2e1 * t214 + ((-qJ(5) * qJD(3) + t52) * t145 + t196) * qJD(2) + t154 (-t142 * t147 - t262) * qJDD(2), t7 * t147 - t8 * qJ(4) - t19 * t34 - t25 * t52 - g(1) * (-pkin(4) * t28 + t206) - g(2) * (-pkin(4) * t26 + t207) - g(3) * (-pkin(4) * t65 + t205) + t247 * t22, t256 * t72 - t268 (-t17 - t267) * t144 + (-t18 - t266) * t141 (-t221 + t269) * qJD(2) + t178 (t222 - t270) * qJD(2) + t177, -t107 * t238, -t140 * t18 + t161 * t141 - t144 * t279 + t186 * t238 + t247 * t71, t140 * t17 + t141 * t279 + t161 * t144 + t5 * t238 + t247 * t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, t226, t101, -qJD(3) * t38 + t240 * t44 - t160 - t257, t101, t83, -t226, t261 + t214 + (-qJ(5) * t234 + t196) * qJD(2) + t154, 0, 0, 0, 0, 0, -t144 * t283 + t260 - t275, t141 * t283 + t259 - t274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193 + t226, 0.2e1 * t211 - t225, -t243 * t149 (-t273 + (t19 + t218) * t142) * qJD(2) + t165 - t167, 0, 0, 0, 0, 0 (-t222 - t270) * qJD(2) - t177 (-t221 - t269) * qJD(2) + t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 * t71, -t71 ^ 2 + t72 ^ 2, t17 - t267, t18 - t266, t70, -t141 * t3 + t1 + t20 * t72 - g(1) * (-t141 * t28 + t144 * t61) - g(2) * (-t141 * t26 + t144 * t59) - g(3) * t35 + t289 * t5, -t144 * t3 - t141 * t2 - t20 * t71 - g(1) * (-t141 * t61 - t144 * t28) - g(2) * (-t141 * t59 - t144 * t26) + g(3) * t36 - t289 * t186;];
tau_reg  = t4;
