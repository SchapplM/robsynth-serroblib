% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRP9
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
% Datum: 2019-12-31 22:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP9_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:06:08
% EndTime: 2019-12-31 22:06:18
% DurationCPUTime: 4.24s
% Computational Cost: add. (5064->461), mult. (11090->601), div. (0->0), fcn. (7546->10), ass. (0->226)
t178 = sin(qJ(2));
t181 = cos(qJ(2));
t220 = pkin(2) * t178 - pkin(7) * t181;
t117 = t220 * qJD(1);
t177 = sin(qJ(3));
t102 = t177 * t117;
t306 = pkin(8) + pkin(7);
t242 = qJD(3) * t306;
t180 = cos(qJ(3));
t277 = t178 * t180;
t278 = t177 * t181;
t329 = -t177 * t242 - t102 - (-pkin(6) * t277 - pkin(8) * t278) * qJD(1);
t272 = t180 * t181;
t215 = pkin(3) * t178 - pkin(8) * t272;
t261 = qJD(1) * t178;
t237 = t177 * t261;
t264 = pkin(6) * t237 + t117 * t180;
t328 = qJD(1) * t215 + t180 * t242 + t264;
t252 = t180 * qJD(2);
t113 = -t237 + t252;
t253 = t177 * qJD(2);
t114 = t180 * t261 + t253;
t176 = sin(qJ(4));
t304 = cos(qJ(4));
t206 = t113 * t176 + t114 * t304;
t58 = -t113 * t304 + t114 * t176;
t129 = -qJD(2) * pkin(2) + pkin(6) * t261;
t76 = -pkin(3) * t113 + t129;
t22 = pkin(4) * t58 - qJ(5) * t206 + t76;
t327 = t22 * t58;
t326 = t58 * t76;
t280 = t176 * t177;
t204 = t180 * t304 - t280;
t251 = t181 * qJD(1);
t311 = qJD(3) + qJD(4);
t235 = t304 * qJD(4);
t314 = qJD(3) * t304 + t235;
t292 = -t180 * t314 + t204 * t251 + t280 * t311;
t116 = t176 * t180 + t177 * t304;
t71 = t311 * t116;
t291 = -t116 * t251 + t71;
t297 = t206 * t58;
t256 = qJD(3) * t180;
t239 = t178 * t256;
t241 = t181 * t253;
t325 = t239 + t241;
t169 = t181 * qJDD(1);
t250 = qJD(1) * qJD(2);
t232 = t178 * t250;
t315 = -t232 + t169;
t110 = qJDD(3) - t315;
t106 = qJDD(4) + t110;
t175 = qJ(3) + qJ(4);
t170 = sin(t175);
t179 = sin(qJ(1));
t182 = cos(qJ(1));
t219 = g(1) * t182 + g(2) * t179;
t298 = g(3) * t181;
t194 = t178 * t219 - t298;
t131 = t306 * t177;
t132 = t306 * t180;
t75 = -t131 * t176 + t132 * t304;
t324 = t75 * t106 + t194 * t170;
t307 = t206 ^ 2;
t323 = -t58 ^ 2 + t307;
t234 = t181 * t250;
t247 = t178 * qJDD(1);
t322 = qJD(2) * qJD(3) + t234 + t247;
t155 = -qJD(3) + t251;
t139 = -qJD(4) + t155;
t227 = t177 * qJDD(2) + t180 * t322;
t249 = qJD(1) * qJD(3);
t233 = t178 * t249;
t195 = -t177 * t233 + t227;
t228 = t177 * t322 + t180 * t233;
t200 = t180 * qJDD(2) - t228;
t255 = qJD(4) * t176;
t16 = -t113 * t235 + t114 * t255 - t176 * t200 - t195 * t304;
t7 = -t139 * t58 - t16;
t30 = pkin(4) * t206 + qJ(5) * t58;
t164 = pkin(6) * t247;
t286 = qJDD(2) * pkin(2);
t94 = pkin(6) * t234 + t164 - t286;
t320 = qJD(3) * pkin(7) * t155 - t94;
t205 = -t131 * t304 - t132 * t176;
t319 = t205 * qJD(4) - t176 * t328 + t304 * t329;
t318 = t75 * qJD(4) + t176 * t329 + t304 * t328;
t124 = t139 * qJD(5);
t92 = t106 * qJ(5);
t317 = t92 - t124;
t166 = pkin(6) * t251;
t257 = qJD(3) * t177;
t303 = pkin(3) * t177;
t221 = pkin(3) * t257 - t251 * t303 - t166;
t269 = t182 * t177;
t100 = t179 * t180 - t181 * t269;
t268 = t182 * t180;
t275 = t179 * t181;
t98 = t177 * t275 + t268;
t313 = -g(1) * t100 + g(2) * t98;
t95 = t106 * pkin(4);
t312 = t95 - qJDD(5);
t125 = -pkin(2) * t181 - pkin(7) * t178 - pkin(1);
t107 = t125 * qJD(1);
t130 = qJD(2) * pkin(7) + t166;
t273 = t180 * t130;
t120 = t220 * qJD(2);
t72 = qJD(1) * t120 + qJDD(1) * t125;
t64 = t180 * t72;
t93 = pkin(6) * t315 + qJDD(2) * pkin(7);
t10 = -t177 * t93 + t64 - t227 * pkin(8) + t110 * pkin(3) + (-t273 + (pkin(8) * t261 - t107) * t177) * qJD(3);
t202 = t107 * t256 - t130 * t257 + t177 * t72 + t180 * t93;
t13 = pkin(8) * t200 + t202;
t67 = t107 * t180 - t130 * t177;
t46 = -pkin(8) * t114 + t67;
t42 = -pkin(3) * t155 + t46;
t68 = t107 * t177 + t273;
t47 = pkin(8) * t113 + t68;
t230 = -t10 * t304 + t13 * t176 + t235 * t47 + t255 * t42;
t282 = t170 * t178;
t171 = cos(t175);
t270 = t182 * t171;
t84 = t170 * t275 + t270;
t271 = t182 * t170;
t276 = t179 * t171;
t86 = t181 * t271 - t276;
t197 = g(1) * t86 + g(2) * t84 + g(3) * t282 - t230;
t187 = t206 * t22 - t197 - t312;
t310 = -t206 * t76 + t197;
t17 = qJD(4) * t206 + t176 * t195 - t200 * t304;
t309 = -t139 * t206 - t17;
t112 = t180 * t125;
t66 = -pkin(8) * t277 + t112 + (-pkin(6) * t177 - pkin(3)) * t181;
t157 = pkin(6) * t272;
t263 = t125 * t177 + t157;
t279 = t177 * t178;
t73 = -pkin(8) * t279 + t263;
t209 = t176 * t66 + t304 * t73;
t172 = t178 * pkin(6);
t265 = t120 * t180 + t172 * t253;
t31 = t215 * qJD(2) + (-t157 + (pkin(8) * t178 - t125) * t177) * qJD(3) + t265;
t199 = -t178 * t252 - t181 * t257;
t266 = t120 * t177 + t125 * t256;
t33 = t199 * pkin(6) - pkin(8) * t325 + t266;
t308 = -qJD(4) * t209 - t176 * t33 + t304 * t31;
t299 = g(3) * t178;
t296 = pkin(4) * t291 + qJ(5) * t292 - qJD(5) * t116 + t221;
t295 = -qJ(5) * t261 + t319;
t294 = pkin(4) * t261 + t318;
t244 = t304 * t47;
t19 = t176 * t42 + t244;
t290 = t139 * t19;
t289 = t176 * t47;
t21 = t304 * t46 - t289;
t287 = pkin(3) * t235 + qJD(5) - t21;
t285 = t113 * t155;
t284 = t114 * t180;
t283 = t114 * t181;
t281 = t171 * t178;
t274 = t180 * t110;
t18 = t304 * t42 - t289;
t267 = qJD(5) - t18;
t156 = pkin(3) * t279;
t121 = t172 + t156;
t173 = t178 ^ 2;
t262 = -t181 ^ 2 + t173;
t260 = qJD(2) * t113;
t259 = qJD(2) * t178;
t258 = qJD(2) * t181;
t254 = t114 * qJD(2);
t77 = pkin(3) * t325 + pkin(6) * t258;
t163 = pkin(3) * t180 + pkin(2);
t243 = pkin(6) + t303;
t240 = t178 * t257;
t238 = t155 * t256;
t231 = t10 * t176 + t13 * t304 + t235 * t42 - t255 * t47;
t229 = -qJD(3) * t107 - t93;
t226 = t304 * t258;
t20 = t176 * t46 + t244;
t225 = pkin(3) * t255 - t20;
t224 = -pkin(4) * t282 + qJ(5) * t281;
t223 = -g(1) * t84 + g(2) * t86;
t85 = t171 * t275 - t271;
t87 = t170 * t179 + t181 * t270;
t222 = g(1) * t85 - g(2) * t87;
t218 = g(1) * t179 - g(2) * t182;
t217 = t130 * t256 - t64;
t216 = -pkin(7) * t110 + qJD(3) * t129;
t214 = t163 * t181 + t178 * t306 + pkin(1);
t213 = pkin(4) * t171 + qJ(5) * t170 + t163;
t211 = -t176 * t73 + t304 * t66;
t208 = -0.2e1 * pkin(1) * t250 - pkin(6) * qJDD(2);
t207 = t176 * t31 + t235 * t66 - t255 * t73 + t304 * t33;
t203 = t177 * t110 - t238;
t185 = qJD(1) ^ 2;
t201 = pkin(1) * t185 + t219;
t198 = t205 * t106 - t171 * t298 + (g(1) * t270 + g(2) * t276) * t178;
t196 = g(1) * t87 + g(2) * t85 + g(3) * t281 - t231;
t184 = qJD(2) ^ 2;
t193 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t184 + t218;
t192 = -g(1) * (-pkin(4) * t86 + qJ(5) * t87) - g(2) * (-pkin(4) * t84 + t85 * qJ(5));
t191 = t195 * t178;
t189 = -t139 * t18 + t196;
t41 = -pkin(3) * t200 + t94;
t162 = -pkin(3) * t304 - pkin(4);
t158 = pkin(3) * t176 + qJ(5);
t101 = t177 * t179 + t181 * t268;
t99 = -t179 * t272 + t269;
t90 = t204 * t178;
t89 = t116 * t178;
t56 = -pkin(4) * t204 - qJ(5) * t116 - t163;
t43 = pkin(4) * t89 - qJ(5) * t90 + t121;
t35 = t177 * t226 - t176 * t240 - t255 * t279 + (t176 * t258 + t178 * t314) * t180;
t34 = t176 * t241 + t178 * t71 - t180 * t226;
t29 = pkin(4) * t181 - t211;
t28 = -qJ(5) * t181 + t209;
t24 = pkin(3) * t114 + t30;
t15 = -qJ(5) * t139 + t19;
t14 = pkin(4) * t139 + t267;
t6 = pkin(4) * t35 + qJ(5) * t34 - qJD(5) * t90 + t77;
t5 = -pkin(4) * t259 - t308;
t4 = qJ(5) * t259 - qJD(5) * t181 + t207;
t3 = t17 * pkin(4) + t16 * qJ(5) - qJD(5) * t206 + t41;
t2 = t230 - t312;
t1 = t231 + t317;
t8 = [qJDD(1), t218, t219, qJDD(1) * t173 + 0.2e1 * t181 * t232, 0.2e1 * t169 * t178 - 0.2e1 * t250 * t262, qJDD(2) * t178 + t181 * t184, qJDD(2) * t181 - t178 * t184, 0, t178 * t208 + t181 * t193, -t178 * t193 + t181 * t208, t180 * t191 + (t181 * t252 - t240) * t114, (t113 * t180 - t114 * t177) * t258 + (t180 * t200 - t227 * t177 + (-t284 + (-t113 + t237) * t177) * qJD(3)) * t178, (-t155 * t252 - t227) * t181 + (t254 + t274 + (t155 + t251) * t257) * t178, (t155 * t253 - t200) * t181 + (-t203 + t260) * t178, -t110 * t181 - t155 * t259, -(-t125 * t257 + t265) * t155 + t112 * t110 - g(1) * t99 - g(2) * t101 + ((t238 - t260) * pkin(6) + (-pkin(6) * t110 + qJD(2) * t129 - t229) * t177 + t217) * t181 + (-pkin(6) * t200 + t67 * qJD(2) + t129 * t256 + t94 * t177) * t178, t266 * t155 - t263 * t110 - g(1) * t98 - g(2) * t100 + (t129 * t252 + t202) * t181 + (-qJD(2) * t68 - t129 * t257 + t94 * t180) * t178 + (t155 * t199 + t181 * t254 + t191) * pkin(6), -t16 * t90 - t206 * t34, t16 * t89 - t17 * t90 - t206 * t35 + t34 * t58, t106 * t90 + t139 * t34 + t16 * t181 + t206 * t259, -t106 * t89 + t139 * t35 + t17 * t181 - t259 * t58, -t106 * t181 - t139 * t259, t106 * t211 + t121 * t17 - t139 * t308 + t18 * t259 + t181 * t230 + t76 * t35 + t41 * t89 + t77 * t58 + t222, -t106 * t209 - t121 * t16 + t139 * t207 + t181 * t231 - t19 * t259 + t206 * t77 - t34 * t76 + t41 * t90 + t223, -t106 * t29 + t139 * t5 - t14 * t259 + t17 * t43 + t181 * t2 + t22 * t35 + t3 * t89 + t58 * t6 + t222, -t1 * t89 - t14 * t34 - t15 * t35 - t29 * t16 - t28 * t17 + t178 * t218 + t2 * t90 + t206 * t5 - t4 * t58, -t1 * t181 + t106 * t28 - t139 * t4 + t15 * t259 + t16 * t43 - t206 * t6 + t22 * t34 - t3 * t90 - t223, t1 * t28 + t15 * t4 + t3 * t43 + t22 * t6 + t2 * t29 + t14 * t5 - g(1) * (-t85 * pkin(4) - t84 * qJ(5)) - g(2) * (t87 * pkin(4) + t86 * qJ(5)) + (-g(1) * t243 - g(2) * t214) * t182 + (g(1) * t214 - g(2) * t243) * t179; 0, 0, 0, -t178 * t185 * t181, t262 * t185, t247, t169, qJDD(2), t178 * t201 - t164 - t298, t299 + (-pkin(6) * qJDD(1) + t201) * t181, -t155 * t284 + t177 * t195, (t227 - t285) * t180 + (-t114 * qJD(3) + (-t239 + t283) * qJD(1) + t200) * t177, (-t114 * t178 + t155 * t272) * qJD(1) + t203, t155 * t257 + t274 + (-t113 * t178 - t155 * t278) * qJD(1), t155 * t261, -pkin(2) * t228 + t264 * t155 + t216 * t177 + (-t67 * t178 + (pkin(6) * t113 - t129 * t177) * t181) * qJD(1) + (t194 + t286 + t320) * t180, -pkin(2) * t227 - t102 * t155 + t216 * t180 + (-t129 * t272 + t68 * t178 + (t155 * t277 - t283) * pkin(6)) * qJD(1) + (t298 + (pkin(2) * t249 - t219) * t178 - t320) * t177, -t16 * t116 - t206 * t292, -t116 * t17 - t16 * t204 - t206 * t291 + t292 * t58, t116 * t106 + t139 * t292 - t206 * t261, t106 * t204 + t139 * t291 + t261 * t58, t139 * t261, t139 * t318 - t163 * t17 - t18 * t261 - t204 * t41 + t221 * t58 + t291 * t76 + t198, t41 * t116 + t139 * t319 + t163 * t16 + t19 * t261 + t206 * t221 - t292 * t76 - t324, t139 * t294 + t14 * t261 + t56 * t17 - t204 * t3 + t22 * t291 + t296 * t58 + t198, t1 * t204 + t2 * t116 - t14 * t292 - t15 * t291 + t16 * t205 - t75 * t17 - t181 * t219 + t206 * t294 - t295 * t58 - t299, -t3 * t116 - t139 * t295 - t15 * t261 + t56 * t16 - t206 * t296 + t22 * t292 + t324, t1 * t75 - t2 * t205 + t3 * t56 + t296 * t22 + t295 * t15 + t294 * t14 + (-g(3) * t213 - t219 * t306) * t181 + (-g(3) * t306 + t213 * t219) * t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114 * t113, -t113 ^ 2 + t114 ^ 2, t195 + t285, -t114 * t155 + t200, t110, -t129 * t114 - t68 * t155 + (t229 + t299) * t177 - t217 + t313, g(1) * t101 - g(2) * t99 + g(3) * t277 - t113 * t129 - t155 * t67 - t202, t297, t323, t7, t309, t106, -t20 * t139 + (t106 * t304 - t114 * t58 + t139 * t255) * pkin(3) + t310, -t21 * t139 + t326 + (-t106 * t176 - t114 * t206 + t139 * t235) * pkin(3) + t196, -t162 * t106 + t139 * t225 - t24 * t58 - t187, -t158 * t17 - t162 * t16 + (t15 + t225) * t206 + (t14 - t287) * t58, t158 * t106 - t139 * t287 + t206 * t24 - t196 + t317 - t327, t1 * t158 + t2 * t162 - t22 * t24 - t14 * t20 - g(3) * (-t156 + t224) + t287 * t15 + (t14 * t255 + t313) * pkin(3) + t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t297, t323, t7, t309, t106, -t290 + t310, t189 + t326, -t30 * t58 - t187 - t290 + t95, pkin(4) * t16 - t17 * qJ(5) + (t15 - t19) * t206 + (t14 - t267) * t58, t206 * t30 - 0.2e1 * t124 - t189 - t327 + 0.2e1 * t92, -t2 * pkin(4) - g(3) * t224 + t1 * qJ(5) - t14 * t19 + t15 * t267 - t22 * t30 + t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106 + t297, t7, -t139 ^ 2 - t307, t139 * t15 + t187;];
tau_reg = t8;
