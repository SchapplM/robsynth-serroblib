% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 18:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPRR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:44:54
% EndTime: 2019-05-05 18:45:04
% DurationCPUTime: 3.25s
% Computational Cost: add. (14497->349), mult. (29159->464), div. (0->0), fcn. (17669->10), ass. (0->207)
t206 = sin(qJ(6));
t212 = cos(qJ(3));
t253 = qJD(1) * qJD(3);
t188 = t212 * t253;
t208 = sin(qJ(3));
t190 = t208 * qJDD(1);
t165 = t190 + t188;
t151 = qJDD(5) + t165;
t145 = qJDD(6) + t151;
t207 = sin(qJ(5));
t211 = cos(qJ(5));
t256 = qJD(1) * t212;
t156 = qJD(3) * t207 + t211 * t256;
t158 = qJD(3) * t211 - t207 * t256;
t210 = cos(qJ(6));
t121 = t156 * t210 + t158 * t206;
t123 = -t156 * t206 + t158 * t210;
t93 = t123 * t121;
t288 = -t93 + t145;
t293 = t206 * t288;
t128 = t158 * t156;
t287 = -t128 + t151;
t292 = t207 * t287;
t291 = t210 * t288;
t290 = t211 * t287;
t289 = t165 + t188;
t215 = qJD(1) ^ 2;
t244 = t208 * t253;
t251 = t212 * qJDD(1);
t166 = -t244 + t251;
t117 = -qJD(5) * t156 + qJDD(3) * t211 - t166 * t207;
t255 = t208 * qJD(1);
t183 = qJD(5) + t255;
t138 = t183 * t156;
t286 = -t138 + t117;
t100 = t138 + t117;
t214 = qJD(3) ^ 2;
t209 = sin(qJ(1));
t213 = cos(qJ(1));
t243 = g(1) * t209 - g(2) * t213;
t159 = qJDD(1) * pkin(1) + t243;
t235 = g(1) * t213 + g(2) * t209;
t161 = -pkin(1) * t215 - t235;
t202 = sin(pkin(10));
t203 = cos(pkin(10));
t257 = t159 * t202 + t161 * t203;
t114 = -pkin(2) * t215 + qJDD(1) * pkin(7) + t257;
t265 = t208 * qJ(4);
t278 = t212 * pkin(3);
t239 = t215 * (-t265 - t278) + t114;
t285 = -pkin(3) * t214 + t212 * t239;
t198 = t212 ^ 2;
t193 = t198 * t215;
t258 = t212 * t215;
t247 = t208 * t258;
t172 = qJDD(3) - t247;
t259 = t212 * t172;
t284 = t259 + t208 * (t193 - t214);
t197 = t208 ^ 2;
t192 = t197 * t215;
t177 = -t192 - t214;
t283 = t177 * t208 + t259;
t119 = t121 ^ 2;
t120 = t123 ^ 2;
t149 = t156 ^ 2;
t150 = t158 ^ 2;
t175 = qJD(6) + t183;
t174 = t175 ^ 2;
t180 = t183 ^ 2;
t282 = 2 * qJD(4);
t281 = -pkin(3) - pkin(8);
t173 = pkin(4) * t255 - qJD(3) * pkin(8);
t240 = t159 * t203 - t202 * t161;
t113 = -qJDD(1) * pkin(2) - pkin(7) * t215 - t240;
t220 = -t166 * pkin(3) - qJ(4) * t289 + t113;
t242 = pkin(3) * qJD(3) - (2 * qJD(4));
t71 = -pkin(4) * t193 - t166 * pkin(8) + (-t173 + t242) * t255 + t220;
t199 = -g(3) + qJDD(2);
t187 = t212 * t199;
t236 = -qJDD(3) * pkin(3) - qJ(4) * t214 + qJDD(4) - t187;
t76 = -qJDD(3) * pkin(8) + (t165 - t188) * pkin(4) + (-pkin(8) * t258 + t239) * t208 + t236;
t40 = t207 * t71 - t211 * t76;
t29 = pkin(5) * t287 - pkin(9) * t100 - t40;
t238 = t207 * qJDD(3) + t166 * t211;
t116 = -qJD(5) * t158 - t238;
t135 = pkin(5) * t183 - pkin(9) * t158;
t41 = t207 * t76 + t211 * t71;
t36 = -pkin(5) * t149 + pkin(9) * t116 - t135 * t183 + t41;
t13 = t206 * t36 - t210 * t29;
t14 = t206 * t29 + t210 * t36;
t9 = -t13 * t210 + t14 * t206;
t280 = t207 * t9;
t279 = t211 * t9;
t252 = qJDD(3) * qJ(4);
t261 = t208 * t199;
t75 = t252 + t261 + t166 * pkin(4) - pkin(8) * t193 + (t282 + t173) * qJD(3) + t285;
t42 = -pkin(5) * t116 - pkin(9) * t149 + t158 * t135 + t75;
t277 = t206 * t42;
t83 = t93 + t145;
t276 = t206 * t83;
t275 = t207 * t75;
t221 = (-qJD(5) + t183) * t158 - t238;
t69 = -t100 * t211 + t207 * t221;
t274 = t208 * t69;
t273 = t210 * t42;
t272 = t210 * t83;
t271 = t211 * t75;
t270 = t175 * t206;
t269 = t175 * t210;
t268 = t183 * t207;
t267 = t183 * t211;
t111 = t128 + t151;
t266 = t207 * t111;
t171 = qJDD(3) + t247;
t264 = t208 * t171;
t260 = t211 * t111;
t254 = qJD(6) + t175;
t250 = -t150 - t180;
t249 = t208 * t93;
t248 = t208 * t128;
t168 = (t197 + t198) * qJDD(1);
t169 = t192 + t193;
t246 = pkin(1) * (t168 * t202 + t169 * t203) + pkin(7) * t168 + pkin(2) * t169;
t245 = pkin(1) * t202 + pkin(7);
t10 = t13 * t206 + t14 * t210;
t103 = t208 * t114 - t187;
t104 = t212 * t114 + t261;
t77 = t103 * t208 + t104 * t212;
t241 = -t116 * t210 + t206 * t117;
t87 = -t174 - t119;
t46 = t206 * t87 + t291;
t237 = pkin(5) * t46 - t13;
t21 = t207 * t41 - t211 * t40;
t234 = t207 * t40 + t211 * t41;
t232 = t206 * t116 + t210 * t117;
t231 = t212 * (-t192 + t214) + t264;
t179 = -t193 - t214;
t230 = t171 * t212 + t179 * t208;
t229 = -t179 * t212 + t264;
t228 = t172 * t208 - t177 * t212;
t226 = -pkin(1) * t203 - pkin(2) - t265;
t102 = -t120 - t174;
t62 = t102 * t210 - t276;
t225 = pkin(5) * t62 - t14;
t223 = (-qJD(6) + t175) * t123 - t241;
t73 = -qJD(6) * t121 + t232;
t222 = -t226 + t278;
t219 = t212 * t281 + t226;
t218 = qJD(3) * t282 + t285;
t91 = t208 * t239 + t236;
t217 = t218 + t252;
t216 = t242 * t255 + t220;
t170 = t192 - t193;
t164 = t190 + 0.2e1 * t188;
t163 = -0.2e1 * t244 + t251;
t137 = -t150 + t180;
t136 = t149 - t180;
t130 = t289 * t208;
t129 = (t166 - t244) * t212;
t127 = t150 - t149;
t125 = t163 * t208 + t164 * t212;
t118 = -t180 - t149;
t109 = -t149 - t150;
t108 = t175 * t121;
t106 = -t120 + t174;
t105 = t119 - t174;
t95 = (qJD(5) + t183) * t158 + t238;
t92 = t120 - t119;
t89 = t211 * t250 - t266;
t88 = t217 + t261;
t85 = t118 * t207 + t290;
t80 = (-t121 * t210 + t123 * t206) * t175;
t79 = (-t121 * t206 - t123 * t210) * t175;
t78 = -t119 - t120;
t72 = -qJD(6) * t123 - t241;
t68 = t105 * t210 - t276;
t67 = -t106 * t206 + t291;
t66 = t105 * t206 + t272;
t65 = t106 * t210 + t293;
t63 = -t102 * t206 - t272;
t58 = -t121 * t254 + t232;
t57 = t108 + t73;
t56 = -t108 + t73;
t53 = t123 * t254 + t241;
t51 = -t123 * t270 + t210 * t73;
t50 = t123 * t269 + t206 * t73;
t49 = t121 * t269 - t206 * t72;
t48 = t121 * t270 + t210 * t72;
t47 = t210 * t87 - t293;
t37 = t207 * t63 + t211 * t62;
t35 = t206 * t57 + t210 * t223;
t34 = -t206 * t56 - t210 * t53;
t33 = t206 * t223 - t210 * t57;
t32 = -t206 * t53 + t210 * t56;
t31 = pkin(5) * t33;
t26 = t207 * t47 + t211 * t46;
t25 = -pkin(9) * t62 + t273;
t24 = -pkin(9) * t46 + t277;
t20 = -pkin(5) * t58 + pkin(9) * t63 + t277;
t18 = -pkin(5) * t53 + pkin(9) * t47 - t273;
t15 = t207 * t35 + t211 * t33;
t8 = pkin(5) * t9;
t6 = -pkin(5) * t42 + pkin(9) * t10;
t5 = -pkin(9) * t33 - t9;
t4 = -pkin(5) * t78 + pkin(9) * t35 + t10;
t2 = t10 * t207 + t279;
t1 = [0, 0, 0, 0, 0, qJDD(1), t243, t235, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t203 - t202 * t215) + t240, pkin(1) * (-qJDD(1) * t202 - t203 * t215) - t257, 0, pkin(1) * (t202 * t257 + t203 * t240), t130, t125, t231, t129, t284, 0, -t212 * t113 + pkin(2) * t163 - pkin(7) * t229 + pkin(1) * (t163 * t203 - t202 * t229), t208 * t113 - pkin(2) * t164 - pkin(7) * t283 + pkin(1) * (-t164 * t203 - t202 * t283), t77 + t246, -pkin(2) * t113 + pkin(7) * t77 + pkin(1) * (-t113 * t203 + t202 * t77), 0, -t231, -t284, t130, t125, t129, (pkin(3) * t169 + t217) * t212 + (qJ(4) * t169 + t187 + t91) * t208 + t246, -t163 * t222 + t212 * t216 + t229 * t245, t208 * (-pkin(3) * t244 + t255 * t282 - t220) + t245 * t283 + t222 * t164, t245 * (t208 * t91 + t212 * t88) - t222 * t216, t248 + t212 * (-t117 * t207 - t158 * t267), t208 * t127 + t212 * (t207 * t95 - t211 * t286), t208 * t100 + t212 * (-t137 * t211 - t292), -t248 + t212 * (-t116 * t211 - t156 * t268), t208 * t221 + t212 * (-t136 * t207 - t260), t208 * t151 + t212 * (t156 * t207 + t158 * t211) * t183, t208 * (pkin(4) * t85 - t40) + t212 * (pkin(4) * t95 + t271) + t245 * (t208 * t85 + t212 * t95) + t219 * (t118 * t211 - t292), t208 * (pkin(4) * t89 - t41) + t212 * (pkin(4) * t286 - t275) + t245 * (t208 * t89 + t212 * t286) + t219 * (-t207 * t250 - t260), pkin(4) * t274 + t212 * (pkin(4) * t109 - t234) + t245 * (t109 * t212 + t274) + t219 * (t100 * t207 + t211 * t221), t219 * t234 + (pkin(4) + t245) * (t208 * t21 + t212 * t75), t249 + t212 * (-t207 * t51 - t211 * t50), t208 * t92 + t212 * (-t207 * t34 - t211 * t32), t208 * t57 + t212 * (-t207 * t67 - t211 * t65), -t249 + t212 * (-t207 * t49 - t211 * t48), t208 * t223 + t212 * (-t207 * t68 - t211 * t66), t208 * t145 + t212 * (-t207 * t80 - t211 * t79), t208 * (pkin(4) * t26 + t237) + t212 * (pkin(4) * t53 - t211 * t18 - t207 * t24) + t245 * (t208 * t26 + t212 * t53) + t219 * (-t207 * t46 + t211 * t47), t208 * (pkin(4) * t37 + t225) + t212 * (pkin(4) * t58 - t211 * t20 - t207 * t25) + t245 * (t208 * t37 + t212 * t58) + t219 * (-t207 * t62 + t211 * t63), t208 * (pkin(4) * t15 + t31) + t212 * (pkin(4) * t78 - t207 * t5 - t211 * t4) + t245 * (t15 * t208 + t212 * t78) + t219 * (-t207 * t33 + t211 * t35), t208 * (pkin(4) * t2 + t8) + t212 * (pkin(4) * t42 + pkin(9) * t280 - t211 * t6) + t245 * (t2 * t208 + t212 * t42) + t219 * (t10 * t211 - t280); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t199, 0, 0, 0, 0, 0, 0, t230, -t228, 0, -t103 * t212 + t104 * t208, 0, 0, 0, 0, 0, 0, 0, -t230, t228, t208 * t88 - t212 * t91, 0, 0, 0, 0, 0, 0, t208 * t95 - t212 * t85, t208 * t286 - t212 * t89, t109 * t208 - t212 * t69, t208 * t75 - t21 * t212, 0, 0, 0, 0, 0, 0, t208 * t53 - t212 * t26, t208 * t58 - t212 * t37, -t15 * t212 + t208 * t78, -t2 * t212 + t208 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t247, t170, t190, t247, t251, qJDD(3), -t103, -t104, 0, 0, qJDD(3), -t190, -t251, -t247, t170, t247, (-pkin(3) * t208 + qJ(4) * t212) * qJDD(1), -pkin(3) * t171 - qJ(4) * t179 + t91, -pkin(3) * t177 + t261 + (qJDD(3) + t172) * qJ(4) + t218, -pkin(3) * t91 + qJ(4) * t88, t117 * t211 - t158 * t268, -t207 * t286 - t211 * t95, -t137 * t207 + t290, -t116 * t207 + t156 * t267, t136 * t211 - t266, (-t156 * t211 + t158 * t207) * t183, qJ(4) * t95 + t281 * t85 + t275, qJ(4) * t286 + t281 * t89 + t271, qJ(4) * t109 + t281 * t69 - t21, qJ(4) * t75 + t21 * t281, -t207 * t50 + t211 * t51, -t207 * t32 + t211 * t34, -t207 * t65 + t211 * t67, -t207 * t48 + t211 * t49, -t207 * t66 + t211 * t68, -t207 * t79 + t211 * t80, qJ(4) * t53 - t207 * t18 + t211 * t24 + t26 * t281, qJ(4) * t58 - t207 * t20 + t211 * t25 + t281 * t37, qJ(4) * t78 + t15 * t281 - t207 * t4 + t211 * t5, -pkin(9) * t279 + qJ(4) * t42 + t2 * t281 - t207 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t190, t171, t177, t91, 0, 0, 0, 0, 0, 0, t85, t89, t69, t21, 0, 0, 0, 0, 0, 0, t26, t37, t15, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128, t127, t100, -t128, t221, t151, -t40, -t41, 0, 0, t93, t92, t57, -t93, t223, t145, t237, t225, t31, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, t92, t57, -t93, t223, t145, -t13, -t14, 0, 0;];
tauJ_reg  = t1;
