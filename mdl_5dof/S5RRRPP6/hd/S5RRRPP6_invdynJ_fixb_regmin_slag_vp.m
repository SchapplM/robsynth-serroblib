% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRPP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,theta4]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 22:38
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRPP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:37:07
% EndTime: 2021-01-15 22:37:25
% DurationCPUTime: 4.48s
% Computational Cost: add. (4680->463), mult. (10501->589), div. (0->0), fcn. (6963->10), ass. (0->225)
t193 = sin(qJ(2));
t248 = t193 * qJDD(1);
t173 = pkin(6) * t248;
t196 = cos(qJ(2));
t249 = qJD(1) * qJD(2);
t235 = t196 * t249;
t101 = -qJDD(2) * pkin(2) + pkin(6) * t235 + t173;
t250 = t196 * qJD(1);
t160 = -qJD(3) + t250;
t185 = g(3) * t196;
t194 = sin(qJ(1));
t197 = cos(qJ(1));
t223 = g(1) * t197 + g(2) * t194;
t205 = t223 * t193 - t185;
t320 = qJD(3) * pkin(7) * t160 - t101 + t205;
t181 = t196 * qJDD(1);
t233 = t193 * t249;
t313 = -t233 + t181;
t121 = qJDD(3) - t313;
t186 = qJ(3) + pkin(8);
t177 = sin(t186);
t191 = qJ(4) + pkin(7);
t195 = cos(qJ(3));
t143 = t191 * t195;
t189 = sin(pkin(8));
t190 = cos(pkin(8));
t192 = sin(qJ(3));
t236 = t191 * t192;
t83 = t190 * t143 - t189 * t236;
t319 = -t83 * t121 - t177 * t205;
t262 = qJD(1) * t193;
t238 = t192 * t262;
t252 = t195 * qJD(2);
t127 = t238 - t252;
t253 = t192 * qJD(2);
t129 = t195 * t262 + t253;
t66 = t190 * t127 + t189 * t129;
t318 = t160 * t66;
t123 = t189 * t195 + t190 * t192;
t104 = t123 * qJD(3);
t288 = -t123 * t250 + t104;
t255 = qJD(3) * t195;
t257 = qJD(3) * t192;
t105 = -t189 * t257 + t190 * t255;
t245 = t192 * t250;
t281 = t190 * t195;
t287 = t189 * t245 - t250 * t281 + t105;
t256 = qJD(3) * t193;
t317 = -qJD(1) * t256 + qJDD(2);
t61 = ((qJD(3) + t250) * qJD(2) + t248) * t192 - t317 * t195;
t219 = -t189 * t127 + t190 * t129;
t316 = t219 ^ 2;
t135 = pkin(4) * t190 + qJ(5) * t189 + pkin(3);
t141 = -t189 * pkin(4) + qJ(5) * t190;
t80 = -t135 * t192 + t141 * t195;
t315 = -pkin(6) + t80;
t232 = qJD(3) * t191;
t207 = -t192 * qJD(4) - t195 * t232;
t275 = t195 * t196;
t305 = pkin(3) * t193;
t217 = -qJ(4) * t275 + t305;
t224 = pkin(2) * t193 - pkin(7) * t196;
t130 = t224 * qJD(1);
t265 = pkin(6) * t238 + t195 * t130;
t59 = t217 * qJD(1) + t265;
t111 = t192 * t130;
t278 = t193 * t195;
t279 = t192 * t196;
t70 = t111 + (-pkin(6) * t278 - qJ(4) * t279) * qJD(1);
t251 = t195 * qJD(4);
t98 = -t192 * t232 + t251;
t296 = (t207 - t59) * t190 + (t70 - t98) * t189;
t144 = -qJD(2) * pkin(2) + pkin(6) * t262;
t84 = t127 * pkin(3) + qJD(4) + t144;
t24 = t66 * pkin(4) - qJ(5) * t219 + t84;
t314 = t24 * t219;
t175 = pkin(6) * t250;
t225 = -t175 + (-t245 + t257) * pkin(3);
t311 = -t135 * t195 - t141 * t192;
t270 = t197 * t195;
t276 = t194 * t196;
t107 = t192 * t276 + t270;
t271 = t197 * t192;
t109 = t194 * t195 - t196 * t271;
t310 = -g(1) * t109 + g(2) * t107;
t60 = qJD(3) * t252 + (t235 + t248) * t195 + t317 * t192;
t28 = t189 * t60 + t190 * t61;
t29 = -t189 * t61 + t190 * t60;
t309 = t28 * pkin(4) - t29 * qJ(5) - t219 * qJD(5);
t307 = -2 * pkin(1);
t100 = t313 * pkin(6) + qJDD(2) * pkin(7);
t131 = t224 * qJD(2);
t142 = -t196 * pkin(2) - t193 * pkin(7) - pkin(1);
t78 = qJD(1) * t131 + t142 * qJDD(1);
t72 = t195 * t78;
t116 = t142 * qJD(1);
t145 = qJD(2) * pkin(7) + t175;
t77 = t192 * t116 + t195 * t145;
t12 = t121 * pkin(3) - t60 * qJ(4) - t77 * qJD(3) - t129 * qJD(4) - t192 * t100 + t72;
t211 = t195 * t100 + t116 * t255 - t145 * t257 + t192 * t78;
t15 = -t61 * qJ(4) - t127 * qJD(4) + t211;
t3 = t190 * t12 - t189 * t15;
t4 = t189 * t12 + t190 * t15;
t303 = g(1) * t194;
t301 = g(2) * t197;
t300 = g(3) * t193;
t184 = t193 * pkin(6);
t50 = -t127 * qJ(4) + t77;
t47 = t190 * t50;
t76 = t195 * t116 - t192 * t145;
t49 = -t129 * qJ(4) + t76;
t22 = t189 * t49 + t47;
t299 = t22 * t219;
t169 = t191 * t193;
t298 = pkin(1) + t169;
t163 = pkin(6) * t275;
t266 = t195 * t131 + t253 * t184;
t33 = -t193 * t251 + t217 * qJD(2) + (-t163 + (qJ(4) * t193 - t142) * t192) * qJD(3) + t266;
t267 = t192 * t131 + t142 * t255;
t38 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t278 + (-qJD(4) * t193 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t196) * t192 + t267;
t11 = t189 * t33 + t190 * t38;
t297 = t288 * pkin(4) - t287 * qJ(5) - t123 * qJD(5) + t225;
t45 = -t160 * pkin(3) + t49;
t21 = t189 * t45 + t47;
t295 = pkin(4) * t262 - t296;
t32 = t189 * t59 + t190 * t70;
t26 = qJ(5) * t262 + t32;
t54 = t189 * t207 + t190 * t98;
t294 = t54 - t26;
t293 = t54 - t32;
t125 = t195 * t142;
t74 = -qJ(4) * t278 + t125 + (-pkin(6) * t192 - pkin(3)) * t196;
t264 = t192 * t142 + t163;
t280 = t192 * t193;
t81 = -qJ(4) * t280 + t264;
t40 = t189 * t74 + t190 * t81;
t292 = t189 * t50;
t290 = t60 * t192;
t286 = t127 * t160;
t285 = t129 * t160;
t284 = t129 * t195;
t178 = cos(t186);
t277 = t194 * t178;
t274 = t196 * t197;
t273 = t197 * t177;
t272 = t197 * t178;
t23 = t190 * t49 - t292;
t268 = qJD(5) - t23;
t132 = pkin(3) * t280 + t184;
t187 = t193 ^ 2;
t263 = -t196 ^ 2 + t187;
t261 = qJD(2) * t127;
t260 = qJD(2) * t129;
t259 = qJD(2) * t193;
t258 = qJD(2) * t196;
t254 = t144 * qJD(3);
t243 = t196 * t253;
t85 = pkin(3) * t243 + pkin(6) * t258 + t255 * t305;
t172 = t195 * pkin(3) + pkin(2);
t244 = t160 * t252;
t242 = t160 * t257;
t241 = t160 * t255;
t237 = t196 * t252;
t230 = -qJD(3) * t116 - t100;
t90 = t177 * t276 + t272;
t92 = t196 * t273 - t277;
t228 = g(1) * t90 - g(2) * t92;
t91 = t178 * t276 - t273;
t93 = t194 * t177 + t196 * t272;
t227 = g(1) * t91 - g(2) * t93;
t164 = t193 * t303;
t226 = -t193 * t301 + t164;
t222 = t145 * t255 - t72;
t221 = -t66 ^ 2 - t316;
t10 = -t189 * t38 + t190 * t33;
t20 = t190 * t45 - t292;
t39 = -t189 * t81 + t190 * t74;
t220 = -pkin(7) * t121 + t254;
t2 = -t121 * pkin(4) + qJDD(5) - t3;
t218 = -t160 * t219 + t28;
t214 = -pkin(6) * qJDD(2) + t249 * t307;
t213 = t192 * t121 - t241;
t212 = t195 * t121 + t242;
t199 = qJD(1) ^ 2;
t210 = pkin(1) * t199 + t223;
t198 = qJD(2) ^ 2;
t209 = pkin(6) * t198 + qJDD(1) * t307 + t301;
t208 = -t29 - t318;
t82 = t189 * t143 + t190 * t236;
t206 = -t82 * t121 - t178 * t185 + (g(1) * t272 + g(2) * t277) * t193;
t44 = t61 * pkin(3) + qJDD(4) + t101;
t203 = g(1) * t93 + g(2) * t91 + t178 * t300 - t4;
t202 = g(1) * t92 + g(2) * t90 - t22 * t160 + t177 * t300 + t3;
t201 = -t223 * t196 - t83 * t28 + t82 * t29 - t54 * t66 - t300;
t200 = t44 - t205;
t171 = t192 * pkin(3) + pkin(6);
t170 = t191 * t196;
t167 = -t190 * pkin(3) - pkin(4);
t165 = t189 * pkin(3) + qJ(5);
t147 = t172 * t196;
t122 = t189 * t192 - t281;
t115 = t121 * qJ(5);
t110 = t194 * t192 + t196 * t270;
t108 = -t194 * t275 + t271;
t99 = t147 + t298;
t96 = -t189 * t280 + t190 * t278;
t95 = t123 * t193;
t79 = pkin(2) - t311;
t73 = t79 * t196;
t63 = t122 * pkin(4) - t123 * qJ(5) - t172;
t56 = t73 + t298;
t52 = t193 * t104 + t189 * t243 - t190 * t237;
t51 = -t105 * t193 - t189 * t237 - t190 * t243;
t46 = t95 * pkin(4) - t96 * qJ(5) + t132;
t37 = t196 * pkin(4) - t39;
t36 = -t196 * qJ(5) + t40;
t25 = t129 * pkin(3) + pkin(4) * t219 + qJ(5) * t66;
t18 = -t160 * qJ(5) + t21;
t17 = t160 * pkin(4) + qJD(5) - t20;
t16 = -t51 * pkin(4) + t52 * qJ(5) - t96 * qJD(5) + t85;
t7 = -pkin(4) * t259 - t10;
t6 = qJ(5) * t259 - t196 * qJD(5) + t11;
t5 = t44 + t309;
t1 = -t160 * qJD(5) + t115 + t4;
t8 = [qJDD(1), -t301 + t303, t223, t187 * qJDD(1) + 0.2e1 * t196 * t233, 0.2e1 * t193 * t181 - 0.2e1 * t263 * t249, qJDD(2) * t193 + t198 * t196, qJDD(2) * t196 - t198 * t193, 0, t214 * t193 + (-t209 + t303) * t196, t209 * t193 + t214 * t196 - t164, t60 * t278 + (-t192 * t256 + t237) * t129, (-t127 * t195 - t129 * t192) * t258 + (-t290 - t195 * t61 + (t127 * t192 - t284) * qJD(3)) * t193, (-t60 - t244) * t196 + (t212 + t260) * t193, (t160 * t253 + t61) * t196 + (-t213 - t261) * t193, -t121 * t196 - t160 * t259, -(-t142 * t257 + t266) * t160 + t125 * t121 - g(1) * t108 - g(2) * t110 + ((t241 + t261) * pkin(6) + (-pkin(6) * t121 + qJD(2) * t144 - t230) * t192 + t222) * t196 + (pkin(6) * t61 + t76 * qJD(2) + t101 * t192 + t195 * t254) * t193, t267 * t160 - t264 * t121 - g(1) * t107 - g(2) * t109 + (t144 * t252 + (-t242 + t260) * pkin(6) + t211) * t196 + (-t192 * t254 - t77 * qJD(2) + t101 * t195 + (t60 - t244) * pkin(6)) * t193, -t10 * t160 + t39 * t121 + t132 * t28 - t3 * t196 + t20 * t259 + t44 * t95 - t84 * t51 + t85 * t66 + t227, t11 * t160 - t40 * t121 + t132 * t29 + t4 * t196 - t21 * t259 + t219 * t85 + t44 * t96 - t84 * t52 - t228, -t10 * t219 - t11 * t66 + t20 * t52 + t21 * t51 - t40 * t28 - t39 * t29 - t3 * t96 - t4 * t95 + t226, t4 * t40 + t21 * t11 + t3 * t39 + t20 * t10 + t44 * t132 + t84 * t85 - g(1) * (t171 * t197 - t99 * t194) - g(2) * (t171 * t194 + t99 * t197), -t37 * t121 + t16 * t66 + t7 * t160 - t17 * t259 + t2 * t196 - t24 * t51 + t46 * t28 + t5 * t95 + t227, -t1 * t95 - t17 * t52 + t18 * t51 + t2 * t96 + t219 * t7 - t36 * t28 + t37 * t29 - t6 * t66 + t226, -t1 * t196 + t36 * t121 - t16 * t219 - t6 * t160 + t18 * t259 + t24 * t52 - t46 * t29 - t5 * t96 + t228, t1 * t36 + t18 * t6 + t5 * t46 + t24 * t16 + t2 * t37 + t17 * t7 - g(1) * (-t56 * t194 - t315 * t197) - g(2) * (-t315 * t194 + t56 * t197); 0, 0, 0, -t193 * t199 * t196, t263 * t199, t248, t181, qJDD(2), t210 * t193 - t173 - t185, t300 + (-pkin(6) * qJDD(1) + t210) * t196, -t160 * t284 + t290, (t60 + t286) * t195 + (-t61 + t285) * t192, (-t129 * t193 + t160 * t275) * qJD(1) + t213, (t127 * t193 - t160 * t279) * qJD(1) + t212, t160 * t262, -pkin(2) * t61 + t265 * t160 + t220 * t192 + (-t76 * t193 + (-pkin(6) * t127 - t144 * t192) * t196) * qJD(1) + t320 * t195, -pkin(2) * t60 - t111 * t160 + t220 * t195 + (-t144 * t275 + t77 * t193 + (-t129 * t196 + t160 * t278) * pkin(6)) * qJD(1) - t320 * t192, t44 * t122 - t296 * t160 - t172 * t28 - t20 * t262 + t225 * t66 + t288 * t84 + t206, t44 * t123 + t293 * t160 - t172 * t29 + t21 * t262 + t219 * t225 + t287 * t84 + t319, -t4 * t122 - t3 * t123 - t287 * t20 - t288 * t21 - t219 * t296 + t32 * t66 + t201, t4 * t83 - t3 * t82 - t44 * t172 - g(3) * (t147 + t169) + t225 * t84 + t293 * t21 + t296 * t20 - t223 * (-t193 * t172 + t170), t5 * t122 + t295 * t160 + t17 * t262 + t288 * t24 + t63 * t28 + t297 * t66 + t206, -t1 * t122 + t2 * t123 + t287 * t17 - t288 * t18 + t219 * t295 + t26 * t66 + t201, -t5 * t123 - t160 * t294 - t18 * t262 - t219 * t297 - t24 * t287 - t63 * t29 - t319, t1 * t83 + t5 * t63 + t2 * t82 - g(3) * (t73 + t169) - t223 * (-t79 * t193 + t170) + t297 * t24 + t294 * t18 + t295 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129 * t127, -t127 ^ 2 + t129 ^ 2, t60 - t286, -t285 - t61, t121, -t144 * t129 - t77 * t160 + (t230 + t300) * t192 - t222 + t310, g(1) * t110 - g(2) * t108 + g(3) * t278 + t144 * t127 - t76 * t160 - t211, -t84 * t219 + (t121 * t190 - t129 * t66) * pkin(3) + t202, -t23 * t160 + t84 * t66 + (-t121 * t189 - t129 * t219) * pkin(3) + t203, t21 * t219 - t299 + (-t189 * t28 - t190 * t29) * pkin(3) + (-t20 + t23) * t66, t20 * t22 - t21 * t23 + (g(3) * t280 - t84 * t129 + t4 * t189 + t3 * t190 + t310) * pkin(3), -t314 - t25 * t66 - qJDD(5) + (pkin(4) - t167) * t121 + t202, -t165 * t28 + t167 * t29 + t18 * t219 - t299 + (t17 - t268) * t66, t165 * t121 - t24 * t66 + t25 * t219 + t115 + (-0.2e1 * qJD(5) + t23) * t160 - t203, t1 * t165 + t2 * t167 - t24 * t25 - t17 * t22 - g(1) * (-t194 * t311 + t274 * t80) - g(2) * (t311 * t197 + t80 * t276) - t80 * t300 + t268 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t218, -t208, t221, t20 * t219 + t21 * t66 + t200, t218, t221, t208, -t17 * t219 + t18 * t66 + t200 + t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t219 * t66 - t121, t29 - t318, -t160 ^ 2 - t316, t18 * t160 + t314 - g(1) * (t122 * t194 + t123 * t274) - g(2) * (-t122 * t197 + t123 * t276) - t123 * t300 + t2;];
tau_reg = t8;
