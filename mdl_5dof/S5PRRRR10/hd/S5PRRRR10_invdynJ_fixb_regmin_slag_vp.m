% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRRR10
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:26:09
% EndTime: 2019-12-05 17:26:24
% DurationCPUTime: 5.27s
% Computational Cost: add. (3710->437), mult. (9694->667), div. (0->0), fcn. (8480->14), ass. (0->228)
t169 = sin(pkin(5));
t180 = cos(qJ(2));
t271 = qJD(1) * t180;
t141 = qJD(2) * pkin(2) + t169 * t271;
t171 = cos(pkin(6));
t168 = sin(pkin(6));
t172 = cos(pkin(5));
t272 = qJD(1) * t172;
t249 = t168 * t272;
t322 = t141 * t171 + t249;
t179 = cos(qJ(3));
t276 = t179 * t180;
t175 = sin(qJ(3));
t176 = sin(qJ(2));
t281 = t175 * t176;
t202 = -t171 * t281 + t276;
t104 = t202 * t169;
t290 = t168 * t175;
t161 = pkin(8) * t290;
t282 = t171 * t179;
t315 = pkin(2) * t282 - t161;
t321 = -qJD(1) * t104 + t315 * qJD(3);
t205 = (pkin(3) * t175 - pkin(9) * t179) * t168;
t287 = t169 * t176;
t248 = qJD(1) * t287;
t320 = qJD(3) * t205 - t168 * t248;
t270 = qJD(2) * t168;
t246 = t179 * t270;
t319 = qJD(4) - t246;
t178 = cos(qJ(4));
t174 = sin(qJ(4));
t247 = t175 * t270;
t221 = t174 * t247;
t269 = qJD(2) * t171;
t227 = qJD(3) + t269;
t100 = -t178 * t227 + t221;
t93 = qJD(5) + t100;
t318 = -pkin(8) * qJDD(2) * t168 - (qJD(2) * t271 + qJDD(1) * t176) * t169 - t322 * qJD(3);
t258 = qJDD(2) * t171;
t160 = qJDD(3) + t258;
t207 = qJD(4) * t227;
t267 = qJD(3) * t179;
t241 = t174 * t267;
t256 = t175 * qJDD(2);
t263 = qJD(4) * t178;
t45 = -t178 * t160 + (qJD(2) * (t175 * t263 + t241) + t174 * t256) * t168 + t174 * t207;
t279 = t176 * t179;
t280 = t175 * t180;
t204 = t171 * t279 + t280;
t103 = t204 * t169;
t283 = t171 * t175;
t289 = t168 * t179;
t274 = pkin(2) * t283 + pkin(8) * t289;
t301 = -qJD(1) * t103 + qJD(3) * t274;
t111 = pkin(9) * t171 + t274;
t214 = -pkin(3) * t179 - pkin(9) * t175;
t112 = (-pkin(2) + t214) * t168;
t265 = qJD(4) * t174;
t317 = -t111 * t265 + t112 * t263 + t174 * t320 + t178 * t321;
t307 = t111 * t178 + t112 * t174;
t316 = t307 * qJD(4) + t174 * t321 - t178 * t320;
t314 = t171 * t276 - t281;
t102 = t174 * t227 + t178 * t247;
t286 = t169 * t180;
t124 = -t168 * t286 + t171 * t172;
t203 = t171 * t280 + t279;
t75 = t169 * t203 + t172 * t290;
t209 = t124 * t178 - t174 * t75;
t298 = sin(pkin(11));
t230 = t298 * t176;
t170 = cos(pkin(11));
t284 = t170 * t180;
t125 = t172 * t284 - t230;
t229 = t298 * t180;
t285 = t170 * t176;
t126 = t172 * t285 + t229;
t288 = t169 * t170;
t252 = t168 * t288;
t49 = t126 * t179 + (t125 * t171 - t252) * t175;
t127 = -t172 * t229 - t285;
t128 = -t172 * t230 + t284;
t231 = t169 * t298;
t215 = t168 * t231;
t51 = t128 * t179 + (t127 * t171 + t215) * t175;
t76 = -t125 * t168 - t171 * t288;
t77 = -t127 * t168 + t171 * t231;
t199 = g(1) * (-t174 * t51 + t178 * t77) + g(2) * (-t174 * t49 + t178 * t76) + g(3) * t209;
t257 = qJDD(2) * t179;
t158 = t168 * t257;
t260 = qJD(2) * qJD(3);
t239 = t175 * t260;
t113 = t168 * t239 + qJDD(4) - t158;
t159 = qJDD(1) * t286;
t245 = qJD(2) * t287;
t219 = qJD(1) * t245;
t108 = qJDD(2) * pkin(2) + t159 - t219;
t132 = pkin(8) * t270 + t248;
t259 = qJDD(1) * t172;
t237 = t168 * t259;
t268 = qJD(3) * t175;
t188 = -t108 * t283 + t132 * t268 - t175 * t237 + t179 * t318;
t12 = pkin(9) * t160 - t188;
t60 = t132 * t179 + t141 * t283 + t175 * t249;
t47 = pkin(9) * t227 + t60;
t157 = t171 * t272;
t68 = t157 + (qJD(2) * t214 - t141) * t168;
t19 = t174 * t68 + t178 * t47;
t156 = t171 * t259;
t193 = t239 - t257;
t238 = t179 * t260;
t194 = t238 + t256;
t36 = t156 + (pkin(3) * t193 - pkin(9) * t194 - t108) * t168;
t184 = -qJD(4) * t19 - t12 * t174 + t178 * t36;
t4 = -pkin(4) * t113 - t184;
t313 = (pkin(4) * t102 + pkin(10) * t93) * t93 + t199 + t4;
t152 = -pkin(4) * t178 - pkin(10) * t174 - pkin(3);
t41 = qJDD(5) + t45;
t312 = (-t60 + t319 * (pkin(4) * t174 - pkin(10) * t178)) * t93 + t152 * t41;
t173 = sin(qJ(5));
t177 = cos(qJ(5));
t242 = t168 * t267;
t220 = t178 * t242;
t236 = t168 * t256;
t44 = qJD(2) * t220 - qJD(4) * t221 + t174 * t160 + (t207 + t236) * t178;
t71 = t102 * t177 + t173 * t319;
t15 = qJD(5) * t71 - t113 * t177 + t173 * t44;
t59 = -t175 * t132 + t179 * t322;
t17 = pkin(10) * t319 + t19;
t46 = -pkin(3) * t227 - t59;
t20 = t100 * pkin(4) - t102 * pkin(10) + t46;
t210 = t17 * t173 - t177 * t20;
t201 = t12 * t178 + t174 * t36 + t263 * t68 - t265 * t47;
t3 = pkin(10) * t113 + t201;
t206 = t108 * t282 - t132 * t267 + t175 * t318 + t179 * t237;
t13 = -pkin(3) * t160 - t206;
t6 = pkin(4) * t45 - pkin(10) * t44 + t13;
t1 = -t210 * qJD(5) + t173 * t6 + t177 * t3;
t181 = qJD(2) ^ 2;
t69 = t102 * t173 - t177 * t319;
t310 = t69 * t93;
t309 = t71 * t93;
t243 = t168 * t268;
t308 = -pkin(4) * t243 + t316;
t306 = pkin(9) * qJD(4);
t261 = qJD(5) * t177;
t262 = qJD(5) * t173;
t14 = -t102 * t262 + t113 * t173 + t177 * t44 + t261 * t319;
t305 = t14 * t173;
t303 = t173 * t41;
t302 = t177 * t41;
t117 = qJD(2) * t205;
t299 = t117 * t174 + t178 * t59;
t297 = t100 * t319;
t296 = t102 * t319;
t295 = t124 * t168;
t293 = t319 * t174;
t165 = t168 ^ 2;
t292 = t165 * t181;
t291 = t168 * t174;
t278 = t176 * t181;
t277 = t178 * t179;
t275 = qJDD(1) - g(3);
t166 = t175 ^ 2;
t273 = -t179 ^ 2 + t166;
t266 = qJD(4) * t173;
t264 = qJD(4) * t177;
t254 = t93 * t266;
t253 = t93 * t264;
t251 = t173 * t289;
t240 = t168 * t171 * t181;
t228 = t177 * t93;
t226 = qJD(3) + 0.2e1 * t269;
t225 = t160 + t258;
t224 = t165 * t169 * t278;
t222 = t168 * t245;
t110 = t161 + (-pkin(2) * t179 - pkin(3)) * t171;
t129 = -t171 * t178 + t174 * t290;
t130 = t171 * t174 + t178 * t290;
t54 = pkin(4) * t129 - pkin(10) * t130 + t110;
t217 = -pkin(10) * t243 - qJD(5) * t54 - t317;
t56 = -pkin(10) * t289 + t307;
t78 = -qJD(4) * t129 + t220;
t79 = qJD(4) * t130 + t168 * t241;
t216 = -pkin(4) * t79 + pkin(10) * t78 + qJD(5) * t56 - t301;
t212 = g(1) * t128 + g(2) * t126;
t8 = t17 * t177 + t173 * t20;
t53 = t124 * t174 + t178 * t75;
t74 = -t169 * t314 - t172 * t289;
t24 = t173 * t74 + t177 * t53;
t23 = -t173 * t53 + t177 * t74;
t18 = -t174 * t47 + t178 * t68;
t208 = -t111 * t174 + t112 * t178;
t80 = t130 * t173 + t177 * t289;
t48 = -t125 * t282 + t126 * t175 + t179 * t252;
t50 = -t127 * t282 + t128 * t175 - t179 * t215;
t198 = g(1) * t50 + g(2) * t48 + g(3) * t74;
t197 = -g(1) * t51 - g(2) * t49 - g(3) * t75;
t63 = t125 * t179 - t126 * t283;
t65 = t127 * t179 - t128 * t283;
t196 = g(1) * t65 + g(2) * t63 + g(3) * t104;
t191 = -t13 + t198;
t16 = -pkin(4) * t319 - t18;
t187 = -pkin(10) * t41 + (t16 + t18) * t93;
t186 = -pkin(9) * t113 + t319 * t46;
t185 = pkin(9) * qJD(5) * t93 - t198;
t2 = -qJD(5) * t8 - t173 * t3 + t177 * t6;
t183 = (pkin(10) * t247 - qJD(5) * t152 + t299) * t93 + t197;
t92 = -t141 * t168 + t157;
t90 = (t173 * t175 + t177 * t277) * t270;
t89 = t173 * t178 * t246 - t177 * t247;
t81 = t130 * t177 - t251;
t73 = -t108 * t168 + t156;
t72 = t104 * t178 + t287 * t291;
t64 = t127 * t175 + t128 * t282;
t62 = t125 * t175 + t126 * t282;
t55 = pkin(4) * t289 - t208;
t43 = t172 * t242 + (qJD(2) * t202 + qJD(3) * t314) * t169;
t42 = t172 * t243 + (qJD(2) * t204 + qJD(3) * t203) * t169;
t38 = t128 * t291 + t178 * t65;
t37 = t126 * t291 + t178 * t63;
t33 = -qJD(5) * t251 + t130 * t261 + t173 * t78 - t177 * t243;
t32 = -qJD(5) * t80 + t173 * t243 + t177 * t78;
t29 = -pkin(4) * t247 - t117 * t178 + t174 * t59;
t28 = t174 * t77 + t178 * t51;
t26 = t174 * t76 + t178 * t49;
t10 = qJD(4) * t209 + t174 * t222 + t178 * t43;
t9 = qJD(4) * t53 + t174 * t43 - t178 * t222;
t5 = [t275, 0, (qJDD(2) * t180 - t278) * t169, (-qJDD(2) * t176 - t180 * t181) * t169, 0, 0, 0, 0, 0, -t74 * t160 - t179 * t224 + t193 * t295 - t227 * t42, -t160 * t75 + t175 * t224 + t194 * t295 - t227 * t43, 0, 0, 0, 0, 0, t100 * t42 + t113 * t209 - t319 * t9 + t45 * t74, -t10 * t319 + t102 * t42 - t113 * t53 + t44 * t74, 0, 0, 0, 0, 0, (-qJD(5) * t24 - t10 * t173 + t177 * t42) * t93 + t23 * t41 + t9 * t69 - t209 * t15, -(qJD(5) * t23 + t10 * t177 + t173 * t42) * t93 - t24 * t41 + t9 * t71 - t209 * t14; 0, qJDD(2), -g(1) * t127 - g(2) * t125 - g(3) * t286 + t159, -t275 * t287 + t212, (qJDD(2) * t166 + 0.2e1 * t175 * t238) * t165, 0.2e1 * (t179 * t256 - t260 * t273) * t165, (t175 * t225 + t226 * t267) * t168, (t179 * t225 - t226 * t268) * t168, t160 * t171, t315 * t160 + t206 * t171 + (-t179 * t73 + t268 * t92) * t168 + (-pkin(2) * t193 + t179 * t219) * t165 - t196 - t301 * t227, -t274 * t160 + t188 * t171 + g(1) * t64 + g(2) * t62 + g(3) * t103 + (t175 * t73 + t267 * t92) * t168 + (-pkin(2) * t194 - t175 * t219) * t165 - t321 * t227, t102 * t78 + t130 * t44, -t100 * t78 - t102 * t79 - t129 * t44 - t130 * t45, t113 * t130 + t319 * t78 + (t102 * t268 - t179 * t44) * t168, -t113 * t129 - t319 * t79 + (-t100 * t268 + t179 * t45) * t168, (-t113 * t179 + t268 * t319) * t168, t208 * t113 + t110 * t45 + t13 * t129 + t46 * t79 - g(1) * t38 - g(2) * t37 - g(3) * t72 + (-t179 * t184 + t18 * t268) * t168 - t316 * t319 + t301 * t100, -t307 * t113 + t110 * t44 + t13 * t130 + t46 * t78 + t196 * t174 + (t201 * t179 - t19 * t268 + (-g(3) * t287 - t212) * t178) * t168 - t317 * t319 + t301 * t102, t14 * t81 + t32 * t71, -t14 * t80 - t15 * t81 - t32 * t69 - t33 * t71, t129 * t14 + t32 * t93 + t41 * t81 + t71 * t79, -t129 * t15 - t33 * t93 - t41 * t80 - t69 * t79, t129 * t41 + t79 * t93, (-t173 * t56 + t177 * t54) * t41 + t2 * t129 - t210 * t79 + t55 * t15 + t4 * t80 + t16 * t33 - g(1) * (t173 * t64 + t177 * t38) - g(2) * (t173 * t62 + t177 * t37) - g(3) * (t103 * t173 + t177 * t72) + (t173 * t217 - t177 * t216) * t93 + t308 * t69, -(t173 * t54 + t177 * t56) * t41 - t1 * t129 - t8 * t79 + t55 * t14 + t4 * t81 + t16 * t32 - g(1) * (-t173 * t38 + t177 * t64) - g(2) * (-t173 * t37 + t177 * t62) - g(3) * (t103 * t177 - t173 * t72) + (t173 * t216 + t177 * t217) * t93 + t308 * t71; 0, 0, 0, 0, -t175 * t179 * t292, t273 * t292, -t179 * t240 + t236, t175 * t240 + t158, t160, t227 * t60 - t247 * t92 + t198 + t206, t227 * t59 - t246 * t92 + t188 - t197, t174 * t44 + t178 * t296, (t44 - t297) * t178 + (-t296 - t45) * t174, t319 * t263 + t174 * t113 + (-t102 * t175 - t277 * t319) * t270, -t319 * t265 + t178 * t113 + (t100 * t175 + t179 * t293) * t270, -t319 * t247, -t18 * t247 - pkin(3) * t45 - t60 * t100 + (t319 * t59 + t186) * t174 + (-(t117 + t306) * t319 + t191) * t178, -pkin(3) * t44 + t299 * t319 + t19 * t247 - t60 * t102 + t186 * t178 + (t306 * t319 - t191) * t174, t14 * t174 * t177 + (-t174 * t262 + t177 * t263 - t90) * t71, t69 * t90 + t71 * t89 + (-t173 * t71 - t177 * t69) * t263 + (-t305 - t15 * t177 + (t173 * t69 - t177 * t71) * qJD(5)) * t174, -t90 * t93 + (-t14 + t253) * t178 + (-t262 * t93 + t319 * t71 + t302) * t174, t89 * t93 + (t15 - t254) * t178 + (-t261 * t93 - t319 * t69 - t303) * t174, -t178 * t41 + t293 * t93, -t16 * t89 - t29 * t69 + t312 * t177 + t183 * t173 + (t16 * t266 - t2 + (qJD(4) * t69 - t303) * pkin(9) - t185 * t177) * t178 + (t16 * t261 + t4 * t173 - t319 * t210 + (t15 + t254) * pkin(9)) * t174, -t16 * t90 - t29 * t71 - t312 * t173 + t183 * t177 + (t16 * t264 + t1 + (qJD(4) * t71 - t302) * pkin(9) + t185 * t173) * t178 + (-t16 * t262 + t4 * t177 - t319 * t8 + (t14 + t253) * pkin(9)) * t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102 * t100, -t100 ^ 2 + t102 ^ 2, t44 + t297, -t45 + t296, t113, -t102 * t46 + t19 * t319 + t184 - t199, g(1) * t28 + g(2) * t26 + g(3) * t53 + t100 * t46 + t18 * t319 - t201, t228 * t71 + t305, (t14 - t310) * t177 + (-t15 - t309) * t173, -t102 * t71 + t228 * t93 + t303, -t173 * t93 ^ 2 + t102 * t69 + t302, -t93 * t102, -pkin(4) * t15 + t102 * t210 + t173 * t187 - t177 * t313 - t19 * t69, -pkin(4) * t14 + t8 * t102 + t173 * t313 + t177 * t187 - t19 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t69, -t69 ^ 2 + t71 ^ 2, t14 + t310, -t15 + t309, t41, t8 * t93 - t16 * t71 - g(1) * (-t173 * t28 + t177 * t50) - g(2) * (-t173 * t26 + t177 * t48) - g(3) * t23 + t2, -t210 * t93 + t16 * t69 - g(1) * (-t173 * t50 - t177 * t28) - g(2) * (-t173 * t48 - t177 * t26) + g(3) * t24 - t1;];
tau_reg = t5;
