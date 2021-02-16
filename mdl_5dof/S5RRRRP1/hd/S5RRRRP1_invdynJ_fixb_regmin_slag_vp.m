% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:52:44
% EndTime: 2021-01-15 23:52:53
% DurationCPUTime: 2.71s
% Computational Cost: add. (4732->337), mult. (11275->420), div. (0->0), fcn. (8236->12), ass. (0->209)
t169 = sin(qJ(3));
t172 = cos(qJ(2));
t282 = cos(qJ(3));
t232 = qJD(1) * t282;
t170 = sin(qJ(2));
t244 = qJD(1) * t170;
t101 = -t169 * t244 + t172 * t232;
t251 = t169 * t172;
t102 = -qJD(1) * t251 - t170 * t232;
t168 = sin(qJ(4));
t281 = cos(qJ(4));
t205 = -t168 * t101 + t281 * t102;
t286 = t205 ^ 2;
t69 = t281 * t101 + t168 * t102;
t65 = t69 ^ 2;
t13 = -t65 + t286;
t164 = qJD(2) + qJD(3);
t297 = t101 * t164;
t235 = t282 * t170;
t112 = t235 + t251;
t295 = t164 * t112;
t178 = t295 * qJD(1);
t207 = t169 * t170 - t282 * t172;
t176 = -t207 * qJDD(1) - t178;
t239 = t172 * qJDD(1);
t52 = qJDD(1) * t235 + t169 * t239 + t297;
t191 = qJD(4) * t205 - t168 * t52 + t281 * t176;
t157 = qJD(4) + t164;
t261 = t205 * t157;
t7 = t191 - t261;
t230 = qJD(4) * t281;
t242 = qJD(4) * t168;
t204 = t101 * t230 + t102 * t242 + t168 * t176 + t281 * t52;
t263 = t69 * t157;
t6 = t204 - t263;
t162 = t172 * pkin(2);
t294 = -pkin(1) - t162;
t227 = -pkin(4) * t69 + qJD(5);
t128 = t294 * qJD(1);
t84 = -t101 * pkin(3) + t128;
t34 = t227 + t84;
t296 = t205 * t34;
t273 = t205 * t69;
t262 = t69 * qJ(5);
t62 = t205 * qJ(5);
t171 = sin(qJ(1));
t173 = cos(qJ(1));
t212 = g(1) * t173 + g(2) * t171;
t167 = qJ(2) + qJ(3);
t161 = qJ(4) + t167;
t149 = sin(t161);
t163 = qJDD(2) + qJDD(3);
t285 = pkin(6) + pkin(7);
t130 = t285 * t172;
t118 = qJD(1) * t130;
t107 = t282 * t118;
t129 = t285 * t170;
t116 = qJD(1) * t129;
t267 = qJD(2) * pkin(2);
t109 = -t116 + t267;
t206 = -t169 * t109 - t107;
t241 = qJD(1) * qJD(2);
t228 = t172 * t241;
t240 = t170 * qJDD(1);
t82 = qJDD(2) * pkin(2) + t285 * (-t228 - t240);
t229 = t170 * t241;
t83 = t285 * (-t229 + t239);
t190 = qJD(3) * t206 - t169 * t83 + t282 * t82;
t10 = t163 * pkin(3) - t52 * pkin(8) + t190;
t231 = qJD(3) * t282;
t243 = qJD(3) * t169;
t192 = t109 * t231 - t118 * t243 + t169 * t82 + t282 * t83;
t14 = t176 * pkin(8) + t192;
t103 = t169 * t118;
t221 = t282 * t109 - t103;
t96 = t102 * pkin(8);
t50 = t221 + t96;
t38 = t164 * pkin(3) + t50;
t276 = t101 * pkin(8);
t51 = -t206 + t276;
t193 = t168 * t10 + t281 * t14 + t38 * t230 - t51 * t242;
t150 = cos(t161);
t255 = t150 * t173;
t256 = t150 * t171;
t188 = g(1) * t255 + g(2) * t256 + g(3) * t149 - t193;
t186 = -t84 * t69 + t188;
t42 = t281 * t51;
t210 = -t168 * t38 - t42;
t292 = t281 * t10 - t168 * t14;
t195 = qJD(4) * t210 + t292;
t257 = t149 * t173;
t258 = t149 * t171;
t277 = g(3) * t150;
t288 = g(1) * t257 + g(2) * t258 - t277;
t187 = t195 + t288;
t181 = t205 * t84 + t187;
t156 = qJDD(4) + t163;
t152 = t282 * pkin(2) + pkin(3);
t252 = t168 * t169;
t73 = t152 * t230 + (-t169 * t242 + (t281 * t282 - t252) * qJD(3)) * pkin(2);
t234 = t281 * t169;
t98 = pkin(2) * t234 + t168 * t152;
t293 = -t98 * t156 - t73 * t157;
t247 = -t169 * t129 + t282 * t130;
t280 = pkin(3) * t157;
t291 = -t168 * pkin(3) * t156 - t230 * t280;
t146 = t156 * pkin(4);
t290 = qJD(5) * t205 + t146;
t88 = pkin(3) * t207 + t294;
t289 = -pkin(4) * t191 + qJDD(5);
t287 = t290 + t296;
t40 = t168 * t51;
t226 = t281 * t38 - t40;
t11 = t226 + t62;
t8 = t157 * pkin(4) + t11;
t283 = t11 - t8;
t275 = t102 * pkin(3);
t272 = t281 * t50 - t40;
t220 = t169 * t116 - t107;
t54 = t220 - t276;
t248 = -t282 * t116 - t103;
t55 = t96 + t248;
t271 = t168 * t54 + t281 * t55;
t219 = -t282 * t129 - t169 * t130;
t63 = -t112 * pkin(8) + t219;
t64 = -pkin(8) * t207 + t247;
t270 = t168 * t63 + t281 * t64;
t20 = t62 + t271;
t269 = t73 - t20;
t223 = -t168 * t55 + t281 * t54;
t19 = t223 - t262;
t74 = -t152 * t242 + (-t169 * t230 + (-t282 * t168 - t234) * qJD(3)) * pkin(2);
t268 = t74 - t19;
t265 = t204 * qJ(5);
t264 = t191 * qJ(5);
t260 = qJDD(1) * pkin(1);
t259 = t102 * t101;
t159 = cos(t167);
t254 = t159 * t171;
t253 = t159 * t173;
t250 = t69 * qJD(5);
t246 = pkin(3) * t159 + pkin(4) * t150;
t165 = t170 ^ 2;
t245 = -t172 ^ 2 + t165;
t155 = t170 * t267;
t238 = t162 + t246;
t236 = qJD(2) * t285;
t225 = -t168 * t50 - t42;
t222 = -t168 * t64 + t281 * t63;
t217 = -pkin(2) * t252 + t281 * t152;
t216 = -g(1) * t258 + g(2) * t257;
t215 = g(1) * t256 - g(2) * t255;
t37 = -pkin(4) * t205 - t275;
t12 = -t210 + t262;
t214 = -t12 * t205 + t8 * t69;
t158 = sin(t167);
t213 = -pkin(3) * t158 - pkin(4) * t149;
t211 = g(1) * t171 - g(2) * t173;
t209 = -0.2e1 * pkin(1) * t241 - pkin(6) * qJDD(2);
t117 = t170 * t236;
t119 = t172 * t236;
t201 = -t282 * t117 - t169 * t119 - t129 * t231 - t130 * t243;
t29 = -pkin(8) * t295 + t201;
t189 = -t247 * qJD(3) + t169 * t117 - t282 * t119;
t81 = t164 * t207;
t30 = t81 * pkin(8) + t189;
t208 = t168 * t30 + t63 * t230 - t64 * t242 + t281 * t29;
t198 = t281 * t207;
t174 = qJD(2) ^ 2;
t197 = -pkin(6) * t174 + t211 + 0.2e1 * t260;
t175 = qJD(1) ^ 2;
t196 = pkin(1) * t175 - pkin(6) * qJDD(1) + t212;
t194 = -t270 * qJD(4) - t168 * t29 + t281 * t30;
t80 = t281 * t112 - t168 * t207;
t185 = t188 - t264;
t182 = g(1) * t253 + g(2) * t254 + g(3) * t158 - t128 * t101 - t192;
t180 = t187 - t265;
t179 = -g(3) * t159 + t128 * t102 + t212 * t158 + t190;
t72 = pkin(3) * t295 + t155;
t177 = -t34 * t69 + t185 - t250;
t145 = pkin(2) * t229;
t33 = pkin(3) * t178 + t88 * qJDD(1) + t145;
t160 = -qJ(5) - pkin(8) - t285;
t154 = pkin(2) * t244;
t151 = t281 * pkin(3) + pkin(4);
t97 = qJDD(1) * t294 + t145;
t94 = pkin(4) + t217;
t89 = pkin(1) + t238;
t85 = t154 - t275;
t79 = t168 * t112 + t198;
t61 = t74 * t157;
t53 = -t101 ^ 2 + t102 ^ 2;
t43 = t79 * pkin(4) + t88;
t35 = t154 + t37;
t32 = -t102 * t164 + t176;
t31 = t52 - t297;
t25 = t80 * qJD(4) - t168 * t81 + t281 * t295;
t24 = qJD(4) * t198 + t112 * t242 + t168 * t295 + t281 * t81;
t23 = -t79 * qJ(5) + t270;
t22 = -t80 * qJ(5) + t222;
t21 = t25 * pkin(4) + t72;
t16 = t62 + t272;
t15 = t225 - t262;
t5 = t33 + t289;
t4 = t24 * qJ(5) - t80 * qJD(5) + t194;
t3 = -t25 * qJ(5) - t79 * qJD(5) + t208;
t2 = t193 + t250 + t264;
t1 = t195 - t265 + t290;
t9 = [qJDD(1), t211, t212, t165 * qJDD(1) + 0.2e1 * t170 * t228, 0.2e1 * t170 * t239 - 0.2e1 * t245 * t241, qJDD(2) * t170 + t174 * t172, qJDD(2) * t172 - t174 * t170, 0, t170 * t209 + t172 * t197, -t170 * t197 + t172 * t209, t102 * t81 + t52 * t112, -t81 * t101 + t102 * t295 + t176 * t112 - t207 * t52, t112 * t163 - t81 * t164, -t163 * t207 - t164 * t295, 0, g(1) * t254 - g(2) * t253 - t101 * t155 + t128 * t295 + t219 * t163 + t189 * t164 - t294 * t176 + t97 * t207, -t102 * t155 + t97 * t112 - t128 * t81 - t211 * t158 - t247 * t163 - t201 * t164 + t294 * t52, t204 * t80 + t205 * t24, t191 * t80 - t204 * t79 + t205 * t25 - t24 * t69, t80 * t156 - t24 * t157, -t79 * t156 - t25 * t157, 0, t222 * t156 + t194 * t157 - t191 * t88 + t84 * t25 + t33 * t79 - t69 * t72 + t215, -t270 * t156 - t208 * t157 + t204 * t88 - t205 * t72 - t84 * t24 + t33 * t80 + t216, t22 * t156 + t4 * t157 - t191 * t43 - t21 * t69 + t34 * t25 + t5 * t79 + t215, -t23 * t156 - t3 * t157 + t204 * t43 - t205 * t21 - t34 * t24 + t5 * t80 + t216, -t1 * t80 - t12 * t25 + t191 * t23 - t2 * t79 - t204 * t22 + t205 * t4 + t8 * t24 + t3 * t69 - t212, t2 * t23 + t12 * t3 + t1 * t22 + t8 * t4 + t5 * t43 + t34 * t21 - g(1) * (-t173 * t160 - t171 * t89) - g(2) * (-t171 * t160 + t173 * t89); 0, 0, 0, -t170 * t175 * t172, t245 * t175, t240, t239, qJDD(2), -g(3) * t172 + t170 * t196, g(3) * t170 + t172 * t196, t259, t53, t31, t32, t163, -t220 * t164 + (t101 * t244 + t282 * t163 - t164 * t243) * pkin(2) + t179, t248 * t164 + (t102 * t244 - t169 * t163 - t164 * t231) * pkin(2) + t182, t273, t13, t6, t7, t156, t217 * t156 - t223 * t157 + t69 * t85 + t181 + t61, t271 * t157 + t205 * t85 + t186 + t293, t94 * t156 - t19 * t157 + t35 * t69 + t180 + t287 + t61, t20 * t157 + t205 * t35 + t177 + t293, t191 * t98 - t204 * t94 + t205 * t268 + t269 * t69 + t214, t2 * t98 + t1 * t94 - t34 * t35 - g(3) * t238 - t212 * (-t170 * pkin(2) + t213) + t268 * t8 + t269 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t259, t53, t31, t32, t163, -t164 * t206 + t179, t164 * t221 + t182, t273, t13, t6, t7, t156, -t225 * t157 + (-t102 * t69 + t281 * t156 - t157 * t242) * pkin(3) + t181, t272 * t157 - t205 * t275 + t186 + t291, -t265 - t15 * t157 + t151 * t156 + t37 * t69 + (-t42 + (-t38 - t280) * t168) * qJD(4) + t287 + t288 + t292, t16 * t157 + t205 * t37 + t177 + t291, -t15 * t205 - t151 * t204 - t16 * t69 + (t168 * t191 + (-t168 * t205 + t281 * t69) * qJD(4)) * pkin(3) + t214, t1 * t151 - t12 * t16 - t8 * t15 - t34 * t37 - g(3) * t246 - t212 * t213 + (t2 * t168 + (t281 * t12 - t168 * t8) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t273, t13, t6, t7, t156, -t157 * t210 + t181, t226 * t157 + t186, t12 * t157 + 0.2e1 * t146 - (-t227 - t34) * t205 + t180, -t286 * pkin(4) + t11 * t157 - (qJD(5) + t34) * t69 + t185, -pkin(4) * t204 - t283 * t69, -t283 * t12 + (t149 * t212 + t1 - t277 + t296) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191 - t261, t204 + t263, -t65 - t286, -pkin(2) * t239 - pkin(3) * t176 - t12 * t69 - t205 * t8 + t145 - t211 - t260 + t289;];
tau_reg = t9;
