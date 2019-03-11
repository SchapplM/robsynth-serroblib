% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRRR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6,theta2]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 06:55:25
% EndTime: 2019-03-09 06:55:33
% DurationCPUTime: 3.87s
% Computational Cost: add. (7093->368), mult. (15539->485), div. (0->0), fcn. (11763->18), ass. (0->226)
t203 = cos(qJ(6));
t280 = qJD(6) * t203;
t200 = sin(qJ(4));
t201 = sin(qJ(3));
t287 = qJD(1) * t201;
t269 = t200 * t287;
t205 = cos(qJ(4));
t206 = cos(qJ(3));
t286 = qJD(1) * t206;
t270 = t205 * t286;
t129 = -t269 + t270;
t130 = -t200 * t286 - t205 * t287;
t199 = sin(qJ(5));
t204 = cos(qJ(5));
t87 = t204 * t129 + t130 * t199;
t338 = t203 * t87;
t343 = t280 - t338;
t198 = sin(qJ(6));
t190 = qJDD(3) + qJDD(4);
t181 = qJDD(5) + t190;
t191 = qJD(3) + qJD(4);
t184 = qJD(5) + t191;
t233 = t129 * t199 - t204 * t130;
t281 = qJD(6) * t198;
t282 = qJD(5) * t204;
t283 = qJD(5) * t199;
t278 = qJD(1) * qJD(3);
t268 = t206 * t278;
t276 = t206 * qJDD(1);
t277 = t201 * qJDD(1);
t74 = qJD(4) * t270 - t191 * t269 + t200 * t276 + (t268 + t277) * t205;
t239 = t200 * t277 - t205 * t276;
t137 = t200 * t206 + t201 * t205;
t98 = t191 * t137;
t75 = qJD(1) * t98 + t239;
t42 = t129 * t282 + t130 * t283 - t199 * t75 + t204 * t74;
t22 = t198 * t181 + t184 * t280 + t203 * t42 - t233 * t281;
t78 = t184 * t198 + t203 * t233;
t23 = qJD(6) * t78 - t203 * t181 + t198 * t42;
t337 = qJD(6) - t87;
t342 = t198 * t337;
t76 = -t203 * t184 + t198 * t233;
t1 = -t198 * t23 + t22 * t203 - t342 * t78 - t343 * t76;
t43 = qJD(5) * t233 + t199 * t74 + t204 * t75;
t41 = qJDD(6) + t43;
t7 = t203 * t41 + t233 * t76 - t337 * t342;
t195 = qJ(3) + qJ(4);
t189 = qJ(5) + t195;
t174 = sin(t189);
t192 = qJ(1) + pkin(11);
t182 = sin(t192);
t183 = cos(t192);
t335 = g(1) * t183 + g(2) * t182;
t341 = t335 * t174;
t20 = t22 * t198;
t9 = t343 * t78 + t20;
t82 = t337 * t280;
t8 = t198 * t41 - t233 * t78 - t337 * t338 + t82;
t196 = sin(pkin(11));
t170 = pkin(1) * t196 + pkin(7);
t312 = pkin(8) + t170;
t259 = t312 * qJD(1);
t106 = t201 * qJD(2) + t206 * t259;
t102 = t205 * t106;
t105 = t206 * qJD(2) - t259 * t201;
t306 = qJD(3) * pkin(3);
t103 = t105 + t306;
t234 = -t103 * t200 - t102;
t321 = pkin(9) * t129;
t59 = -t234 + t321;
t305 = t199 * t59;
t127 = t130 * pkin(9);
t100 = t200 * t106;
t256 = t205 * t103 - t100;
t58 = t127 + t256;
t55 = pkin(4) * t191 + t58;
t30 = t204 * t55 - t305;
t28 = -pkin(5) * t184 - t30;
t317 = t28 * t87;
t175 = cos(t189);
t318 = g(3) * t175;
t185 = t206 * qJDD(2);
t149 = t170 * qJDD(1);
t252 = pkin(8) * qJDD(1) + t149;
t66 = qJDD(3) * pkin(3) - qJD(3) * t106 - t252 * t201 + t185;
t72 = qJD(3) * t105 + t201 * qJDD(2) + t252 * t206;
t222 = qJD(4) * t234 - t200 * t72 + t205 * t66;
t17 = t190 * pkin(4) - t74 * pkin(9) + t222;
t285 = qJD(4) * t200;
t324 = (qJD(4) * t103 + t72) * t205 - t106 * t285 + t200 * t66;
t18 = -t75 * pkin(9) + t324;
t302 = t204 * t59;
t31 = t199 * t55 + t302;
t327 = qJD(5) * t31 - t204 * t17 + t199 * t18;
t4 = -t181 * pkin(5) + t327;
t339 = t4 + t318;
t313 = t233 * t87;
t34 = t233 ^ 2 - t87 ^ 2;
t24 = -t184 * t87 + t42;
t167 = g(3) * t174;
t325 = (qJD(5) * t55 + t18) * t204 + t199 * t17 - t59 * t283;
t197 = cos(pkin(11));
t171 = -pkin(1) * t197 - pkin(2);
t146 = -pkin(3) * t206 + t171;
t131 = t146 * qJD(1);
t95 = -t129 * pkin(4) + t131;
t215 = t335 * t175 - t95 * t87 + t167 - t325;
t333 = pkin(5) * t233;
t176 = pkin(4) * t199 + pkin(10);
t322 = pkin(4) * t130;
t53 = -pkin(10) * t87 - t322 + t333;
t332 = (qJD(6) * t176 + t53) * t337;
t178 = pkin(3) * t205 + pkin(4);
t294 = t200 * t204;
t289 = pkin(3) * t294 + t199 * t178;
t124 = pkin(10) + t289;
t179 = pkin(3) * t287;
t331 = (qJD(6) * t124 + t179 + t53) * t337;
t330 = (t337 * pkin(10) + t333) * t337;
t314 = t337 * t233;
t134 = t312 * t201;
t135 = t312 * t206;
t290 = -t200 * t134 + t205 * t135;
t29 = pkin(10) * t184 + t31;
t45 = -pkin(5) * t87 - pkin(10) * t233 + t95;
t12 = -t198 * t29 + t203 * t45;
t232 = -t12 * t233 + t203 * t341 + t28 * t281;
t13 = t198 * t45 + t203 * t29;
t237 = t13 * t233 + t339 * t198 + t28 * t280;
t212 = -t233 * t95 - t318 - t327 + t341;
t25 = t184 * t233 - t43;
t136 = t200 * t201 - t205 * t206;
t93 = t204 * t136 + t137 * t199;
t97 = t191 * t136;
t51 = -qJD(5) * t93 - t199 * t98 - t204 * t97;
t94 = -t136 * t199 + t137 * t204;
t240 = t337 * t51 + t41 * t94;
t272 = t94 * t281;
t326 = -t203 * t240 + t272 * t337;
t260 = qJD(3) * t312;
t125 = t201 * t260;
t126 = t206 * t260;
t284 = qJD(4) * t205;
t226 = -t205 * t125 - t200 * t126 - t134 * t284 - t135 * t285;
t46 = -pkin(9) * t98 + t226;
t219 = -qJD(4) * t290 + t200 * t125 - t205 * t126;
t47 = t97 * pkin(9) + t219;
t247 = -t205 * t134 - t135 * t200;
t79 = -pkin(9) * t137 + t247;
t80 = -pkin(9) * t136 + t290;
t48 = t199 * t80 - t204 * t79;
t10 = -qJD(5) * t48 + t199 * t47 + t204 * t46;
t265 = t181 * pkin(10) + qJD(6) * t45 + t325;
t49 = t199 * t79 + t204 * t80;
t104 = pkin(4) * t136 + t146;
t54 = pkin(5) * t93 - pkin(10) * t94 + t104;
t323 = -(qJD(6) * t54 + t10) * t337 - t265 * t93 + t28 * t51 + t4 * t94 - t49 * t41;
t316 = t28 * t94;
t315 = t54 * t41;
t52 = qJD(5) * t94 - t199 * t97 + t204 * t98;
t311 = t22 * t93 + t78 * t52;
t295 = t199 * t200;
t248 = -t105 * t200 - t102;
t61 = t248 - t321;
t291 = t205 * t105 - t100;
t62 = t127 + t291;
t308 = t199 * t61 + t204 * t62 - t178 * t282 - (-t200 * t283 + (t204 * t205 - t295) * qJD(4)) * pkin(3);
t307 = -t199 * t62 + t204 * t61 + t178 * t283 + (t200 * t282 + (t199 * t205 + t294) * qJD(4)) * pkin(3);
t300 = t130 * t129;
t299 = t182 * t198;
t298 = t182 * t203;
t297 = t183 * t198;
t296 = t183 * t203;
t292 = qJDD(2) - g(3);
t193 = t201 ^ 2;
t288 = -t206 ^ 2 + t193;
t153 = qJD(1) * t171;
t180 = t201 * t306;
t89 = pkin(4) * t98 + t180;
t112 = qJD(3) * t179 + qJDD(1) * t146;
t60 = t75 * pkin(4) + t112;
t6 = pkin(5) * t43 - pkin(10) * t42 + t60;
t266 = qJD(6) * t29 - t6;
t32 = t199 * t58 + t302;
t245 = pkin(4) * t283 - t32;
t243 = g(1) * t182 - g(2) * t183;
t202 = sin(qJ(1));
t207 = cos(qJ(1));
t242 = g(1) * t202 - g(2) * t207;
t241 = -t93 * t23 - t52 * t76;
t238 = t181 * t94 + t184 * t51;
t235 = t137 * t190 - t191 * t97;
t231 = -t265 + t167;
t230 = -pkin(3) * t295 + t178 * t204;
t228 = -pkin(10) * t41 + t30 * t337 - t317;
t225 = -qJD(1) * t153 - t149 + t335;
t224 = -t124 * t41 + t308 * t337 - t317;
t223 = 0.2e1 * qJD(3) * t153 - qJDD(3) * t170;
t33 = t204 * t58 - t305;
t220 = -t176 * t41 - t317 + (-pkin(4) * t282 + t33) * t337;
t208 = qJD(3) ^ 2;
t218 = -0.2e1 * qJDD(1) * t171 - t170 * t208 + t243;
t217 = -t198 * t240 - t82 * t94;
t187 = sin(t195);
t188 = cos(t195);
t214 = g(3) * t187 - t131 * t129 + t335 * t188 - t324;
t211 = -g(3) * t188 + t131 * t130 + t335 * t187 + t222;
t209 = qJD(1) ^ 2;
t177 = -pkin(4) * t204 - pkin(5);
t148 = qJDD(3) * t206 - t201 * t208;
t147 = qJDD(3) * t201 + t206 * t208;
t123 = -pkin(5) - t230;
t111 = t175 * t296 + t299;
t110 = -t175 * t297 + t298;
t109 = -t175 * t298 + t297;
t108 = t175 * t299 + t296;
t107 = t179 - t322;
t81 = -t129 ^ 2 + t130 ^ 2;
t73 = -t136 * t190 - t191 * t98;
t64 = -t239 + (-qJD(1) * t137 - t130) * t191;
t63 = -t129 * t191 + t74;
t44 = -t181 * t93 - t184 * t52;
t15 = pkin(5) * t52 - pkin(10) * t51 + t89;
t11 = qJD(5) * t49 + t199 * t46 - t204 * t47;
t5 = t203 * t6;
t2 = [qJDD(1), t242, g(1) * t207 + g(2) * t202 (t242 + (t196 ^ 2 + t197 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t193 + 0.2e1 * t201 * t268, 0.2e1 * t201 * t276 - 0.2e1 * t278 * t288, t147, t148, 0, t201 * t223 + t206 * t218, -t201 * t218 + t206 * t223, t130 * t97 + t137 * t74, -t129 * t97 + t130 * t98 - t136 * t74 - t137 * t75, t235, t73, 0, t112 * t136 - t129 * t180 + t131 * t98 + t146 * t75 + t188 * t243 + t190 * t247 + t191 * t219, t112 * t137 - t130 * t180 - t131 * t97 + t146 * t74 - t187 * t243 - t190 * t290 - t191 * t226, t233 * t51 + t42 * t94, -t233 * t52 - t42 * t93 - t43 * t94 + t51 * t87, t238, t44, 0, t104 * t43 - t11 * t184 + t175 * t243 - t181 * t48 + t52 * t95 + t60 * t93 - t87 * t89, -t10 * t184 + t104 * t42 - t174 * t243 - t181 * t49 + t233 * t89 + t51 * t95 + t60 * t94, -t78 * t272 + (t22 * t94 + t51 * t78) * t203 (-t198 * t78 - t203 * t76) * t51 + (-t20 - t203 * t23 + (t198 * t76 - t203 * t78) * qJD(6)) * t94, t311 - t326, t217 + t241, t337 * t52 + t41 * t93, -g(1) * t109 - g(2) * t111 + t11 * t76 + t12 * t52 + t48 * t23 + t5 * t93 + (t15 * t337 + t315 + (-t29 * t93 - t337 * t49 + t316) * qJD(6)) * t203 + t323 * t198, -g(1) * t108 - g(2) * t110 + t11 * t78 - t13 * t52 + t48 * t22 + (-(-qJD(6) * t49 + t15) * t337 - t315 + t266 * t93 - qJD(6) * t316) * t198 + t323 * t203; 0, 0, 0, t292, 0, 0, 0, 0, 0, t148, -t147, 0, 0, 0, 0, 0, t73, -t235, 0, 0, 0, 0, 0, t44, -t238, 0, 0, 0, 0, 0, t217 - t241, t311 + t326; 0, 0, 0, 0, -t201 * t209 * t206, t288 * t209, t277, t276, qJDD(3), -g(3) * t206 + t201 * t225 + t185, -t201 * t292 + t225 * t206, t300, t81, t63, t64, t190, -t248 * t191 + (t129 * t287 + t190 * t205 - t191 * t285) * pkin(3) + t211, t291 * t191 + (t130 * t287 - t190 * t200 - t191 * t284) * pkin(3) + t214, -t313, t34, t24, t25, t181, t107 * t87 + t230 * t181 - t307 * t184 + t212, -t107 * t233 - t289 * t181 + t308 * t184 + t215, t9, t1, t8, t7, -t314, t123 * t23 + t307 * t76 + (-t339 - t331) * t203 + t224 * t198 + t232, t123 * t22 + t307 * t78 + t224 * t203 + (-t341 + t331) * t198 + t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t300, t81, t63, t64, t190, -t191 * t234 + t211, t191 * t256 + t214, -t313, t34, t24, t25, t181, t184 * t32 + (-t130 * t87 + t181 * t204 - t184 * t283) * pkin(4) + t212, t184 * t33 + (t130 * t233 - t181 * t199 - t184 * t282) * pkin(4) + t215, t9, t1, t8, t7, -t314, t177 * t23 + t245 * t76 + (-t339 - t332) * t203 + t220 * t198 + t232, t177 * t22 + t245 * t78 + t220 * t203 + (-t341 + t332) * t198 + t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t313, t34, t24, t25, t181, t184 * t31 + t212, t184 * t30 + t215, t9, t1, t8, t7, -t314, -pkin(5) * t23 - t31 * t76 + t228 * t198 + (-t339 - t330) * t203 + t232, -pkin(5) * t22 - t31 * t78 + t228 * t203 + (-t341 + t330) * t198 + t237; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78 * t76, -t76 ^ 2 + t78 ^ 2, t337 * t76 + t22, t337 * t78 - t23, t41, -g(1) * t110 + g(2) * t108 + t13 * t337 + t198 * t231 - t28 * t78 - t280 * t29 + t5, g(1) * t111 - g(2) * t109 + t12 * t337 + t198 * t266 + t203 * t231 + t28 * t76;];
tau_reg  = t2;
