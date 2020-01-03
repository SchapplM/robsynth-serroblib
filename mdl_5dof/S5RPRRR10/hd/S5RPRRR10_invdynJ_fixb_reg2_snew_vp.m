% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRR10
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRR10_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR10_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR10_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:10:47
% EndTime: 2019-12-31 19:10:59
% DurationCPUTime: 5.14s
% Computational Cost: add. (21757->415), mult. (51813->599), div. (0->0), fcn. (38986->10), ass. (0->252)
t211 = sin(pkin(9));
t208 = t211 ^ 2;
t212 = cos(pkin(9));
t209 = t212 ^ 2;
t258 = t208 + t209;
t284 = sin(qJ(1));
t285 = cos(qJ(1));
t227 = t285 * g(1) + t284 * g(2);
t305 = (2 * qJD(2) * qJD(1)) - t227;
t214 = sin(qJ(5));
t219 = cos(qJ(3));
t216 = sin(qJ(3));
t261 = t212 * t216;
t232 = t211 * t219 + t261;
t198 = t232 * qJD(1);
t215 = sin(qJ(4));
t218 = cos(qJ(4));
t180 = -t218 * qJD(3) + t198 * t215;
t182 = qJD(3) * t215 + t198 * t218;
t217 = cos(qJ(5));
t149 = t217 * t180 + t182 * t214;
t151 = -t180 * t214 + t182 * t217;
t116 = t151 * t149;
t250 = t212 * qJDD(1);
t251 = t211 * qJDD(1);
t234 = t216 * t251 - t219 * t250;
t256 = t198 * qJD(3);
t171 = -t234 - t256;
t163 = qJDD(4) - t171;
t162 = qJDD(5) + t163;
t295 = -t116 + t162;
t304 = t214 * t295;
t155 = t182 * t180;
t294 = -t155 + t163;
t303 = t215 * t294;
t262 = t211 * t216;
t196 = (-t212 * t219 + t262) * qJD(1);
t175 = t198 * t196;
t293 = qJDD(3) - t175;
t302 = t216 * t293;
t301 = t217 * t295;
t300 = t218 * t294;
t299 = t219 * t293;
t220 = qJD(1) ^ 2;
t298 = -(t220 * pkin(1)) + qJDD(1) * qJ(2) + t305;
t235 = t284 * g(1) - t285 * g(2);
t229 = -qJDD(2) + t235;
t242 = pkin(2) * t212 + pkin(1);
t166 = t242 * qJDD(1) + (t258 * pkin(6) + qJ(2)) * t220 + t229;
t297 = qJ(2) + pkin(6);
t191 = qJD(4) + t196;
t188 = qJD(5) + t191;
t138 = t188 * t149;
t195 = t232 * qJDD(1);
t257 = qJD(3) * t196;
t173 = t195 - t257;
t231 = -qJDD(3) * t215 - t173 * t218;
t142 = -qJD(4) * t180 - t231;
t238 = -t218 * qJDD(3) + t173 * t215;
t228 = qJD(4) * t182 + t238;
t91 = -t149 * qJD(5) + t217 * t142 - t214 * t228;
t296 = -t138 + t91;
t161 = t191 * t180;
t121 = t142 + t161;
t259 = t220 * qJ(2);
t277 = qJDD(1) * pkin(1);
t192 = t229 + t259 + t277;
t292 = t258 * t259 - t192 - t277;
t239 = t142 * t214 + t217 * t228;
t72 = (qJD(5) - t188) * t151 + t239;
t117 = (qJD(4) - t191) * t182 + t238;
t147 = t149 ^ 2;
t148 = t151 ^ 2;
t289 = t180 ^ 2;
t179 = t182 ^ 2;
t186 = t188 ^ 2;
t190 = t191 ^ 2;
t193 = t196 ^ 2;
t194 = t198 ^ 2;
t288 = qJD(3) ^ 2;
t164 = pkin(3) * t196 - pkin(7) * t198;
t224 = (-t297 * qJDD(1) + t242 * t220 - t305) * t211;
t283 = t212 * g(3);
t223 = t224 - t283;
t236 = -g(3) * t211 + t298 * t212;
t160 = -pkin(2) * t209 * t220 + pkin(6) * t250 + t236;
t260 = t219 * t160;
t103 = -t288 * pkin(3) + qJDD(3) * pkin(7) - t196 * t164 + t216 * t223 + t260;
t109 = (-t173 + t257) * pkin(7) + (-t171 + t256) * pkin(3) - t166;
t63 = t103 * t215 - t218 * t109;
t50 = pkin(4) * t294 - t121 * pkin(8) - t63;
t157 = pkin(4) * t191 - pkin(8) * t182;
t64 = t218 * t103 + t215 * t109;
t54 = -t289 * pkin(4) - t228 * pkin(8) - t191 * t157 + t64;
t24 = t214 * t54 - t217 * t50;
t25 = t214 * t50 + t217 * t54;
t13 = t214 * t25 - t217 * t24;
t287 = pkin(4) * t13;
t75 = t138 + t91;
t44 = -t214 * t72 - t217 * t75;
t286 = pkin(4) * t44;
t282 = t13 * t215;
t281 = t13 * t218;
t125 = t160 * t216 - t219 * t223;
t126 = -g(3) * t261 + t216 * t224 + t260;
t92 = -t125 * t219 + t126 * t216;
t280 = t211 * t92;
t102 = -qJDD(3) * pkin(3) - t288 * pkin(7) + t164 * t198 + t125;
t57 = t228 * pkin(4) - t289 * pkin(8) + t157 * t182 + t102;
t279 = t214 * t57;
t278 = t217 * t57;
t276 = t102 * t215;
t275 = t102 * t218;
t105 = t116 + t162;
t274 = t105 * t214;
t273 = t105 * t217;
t129 = t155 + t163;
t272 = t129 * t215;
t271 = t129 * t218;
t270 = t166 * t216;
t269 = t166 * t219;
t168 = qJDD(3) + t175;
t268 = t168 * t216;
t267 = t168 * t219;
t266 = t188 * t214;
t265 = t188 * t217;
t264 = t191 * t215;
t263 = t191 * t218;
t254 = qJD(4) + t191;
t246 = t216 * t116;
t245 = t219 * t116;
t244 = t216 * t155;
t243 = t219 * t155;
t241 = -pkin(3) * t219 - pkin(2);
t14 = t214 * t24 + t217 * t25;
t37 = t215 * t63 + t218 * t64;
t93 = t125 * t216 + t219 * t126;
t237 = t211 * (t298 * t211 + t283) + t212 * t236;
t36 = t215 * t64 - t218 * t63;
t112 = -t186 - t147;
t65 = t112 * t214 + t301;
t230 = pkin(4) * t65 - t24;
t124 = -t148 - t186;
t81 = t124 * t217 - t274;
t226 = pkin(4) * t81 - t25;
t204 = t209 * qJDD(1);
t203 = t208 * qJDD(1);
t199 = t258 * t220;
t185 = -t194 - t288;
t184 = -t194 + t288;
t183 = t193 - t288;
t172 = t195 - 0.2e1 * t257;
t170 = t234 + 0.2e1 * t256;
t165 = -t288 - t193;
t159 = -t179 + t190;
t158 = -t190 + t289;
t154 = -t193 - t194;
t153 = t179 - t289;
t146 = -t179 - t190;
t145 = -t185 * t216 - t267;
t144 = t185 * t219 - t268;
t143 = -t190 - t289;
t137 = t179 + t289;
t136 = t195 * t216 - t219 * t234;
t135 = -t195 * t219 - t216 * t234;
t134 = -t148 + t186;
t133 = t147 - t186;
t132 = t165 * t219 - t302;
t131 = t165 * t216 + t299;
t127 = (-t180 * t218 + t182 * t215) * t191;
t122 = t254 * t180 + t231;
t120 = t142 - t161;
t118 = -t254 * t182 - t238;
t115 = t148 - t147;
t114 = t142 * t218 - t182 * t264;
t113 = t180 * t263 + t215 * t228;
t111 = t158 * t218 - t272;
t110 = -t159 * t215 + t300;
t101 = -t146 * t215 - t271;
t100 = t146 * t218 - t272;
t98 = (-t149 * t217 + t151 * t214) * t188;
t97 = (-t149 * t214 - t151 * t217) * t188;
t96 = t143 * t218 - t303;
t95 = t143 * t215 + t300;
t94 = -t147 - t148;
t90 = -qJD(5) * t151 - t239;
t89 = -t117 * t218 + t121 * t215;
t88 = t118 * t218 - t120 * t215;
t87 = -t117 * t215 - t121 * t218;
t86 = t133 * t217 - t274;
t85 = -t134 * t214 + t301;
t84 = t133 * t214 + t273;
t83 = t134 * t217 + t304;
t82 = -t124 * t214 - t273;
t80 = t101 * t219 - t122 * t216;
t79 = t101 * t216 + t122 * t219;
t78 = -t118 * t216 + t219 * t96;
t77 = t118 * t219 + t216 * t96;
t71 = (qJD(5) + t188) * t151 + t239;
t70 = -t151 * t266 + t217 * t91;
t69 = t151 * t265 + t214 * t91;
t68 = t149 * t265 - t214 * t90;
t67 = t149 * t266 + t217 * t90;
t66 = t112 * t217 - t304;
t61 = -pkin(7) * t100 + t275;
t60 = -t137 * t216 + t219 * t89;
t59 = t137 * t219 + t216 * t89;
t58 = -pkin(7) * t95 + t276;
t56 = -t215 * t97 + t218 * t98;
t55 = -pkin(3) * t100 + t64;
t53 = -pkin(3) * t95 + t63;
t52 = -t215 * t84 + t218 * t86;
t51 = -t215 * t83 + t218 * t85;
t48 = -t215 * t81 + t218 * t82;
t47 = t215 * t82 + t218 * t81;
t46 = t214 * t75 - t217 * t72;
t45 = -t214 * t296 - t217 * t71;
t43 = -t214 * t71 + t217 * t296;
t42 = -pkin(8) * t81 + t278;
t41 = -t215 * t69 + t218 * t70;
t40 = -t215 * t67 + t218 * t68;
t39 = -t215 * t65 + t218 * t66;
t38 = t215 * t66 + t218 * t65;
t35 = -pkin(8) * t65 + t279;
t32 = t216 * t296 + t219 * t48;
t31 = t216 * t48 - t219 * t296;
t30 = -pkin(7) * t87 - t36;
t29 = -pkin(4) * t296 + pkin(8) * t82 + t279;
t28 = t216 * t71 + t219 * t39;
t27 = t216 * t39 - t219 * t71;
t26 = -pkin(4) * t71 + pkin(8) * t66 - t278;
t22 = -t215 * t44 + t218 * t46;
t21 = -t215 * t43 + t218 * t45;
t20 = t215 * t46 + t218 * t44;
t19 = t216 * t94 + t219 * t22;
t18 = t216 * t22 - t219 * t94;
t17 = -pkin(3) * t20 - t286;
t16 = -pkin(3) * t47 - t226;
t15 = -pkin(3) * t38 - t230;
t12 = -pkin(7) * t47 - t215 * t29 + t218 * t42;
t11 = -pkin(7) * t38 - t215 * t26 + t218 * t35;
t10 = -pkin(4) * t57 + pkin(8) * t14;
t9 = -pkin(8) * t44 - t13;
t8 = -pkin(4) * t94 + pkin(8) * t46 + t14;
t7 = t14 * t218 - t282;
t6 = t14 * t215 + t281;
t5 = t216 * t57 + t219 * t7;
t4 = t216 * t7 - t219 * t57;
t3 = -pkin(3) * t6 - t287;
t2 = -pkin(7) * t20 - t215 * t8 + t218 * t9;
t1 = -pkin(7) * t6 - pkin(8) * t281 - t10 * t215;
t23 = [0, 0, 0, 0, 0, qJDD(1), t235, t227, 0, 0, t203, 0.2e1 * t211 * t250, 0, t204, 0, 0, -t292 * t212, t292 * t211, pkin(1) * t199 + qJ(2) * (t204 + t203) + t237, pkin(1) * t192 + qJ(2) * t237, t211 * (t173 * t219 - t216 * t256) + t212 * (t173 * t216 + t219 * t256), t211 * (-t170 * t219 - t172 * t216) + t212 * (-t170 * t216 + t172 * t219), t211 * (-t184 * t216 + t299) + t212 * (t184 * t219 + t302), t211 * (-t171 * t216 + t219 * t257) + t212 * (t171 * t219 + t216 * t257), t211 * (t183 * t219 - t268) + t212 * (t183 * t216 + t267), (t211 * (-t196 * t219 + t198 * t216) + t212 * (-t196 * t216 - t198 * t219)) * qJD(3), t211 * (-pkin(6) * t131 - t270) + t212 * (-pkin(2) * t170 + pkin(6) * t132 + t269) - pkin(1) * t170 + qJ(2) * (-t131 * t211 + t132 * t212), t211 * (-pkin(6) * t144 - t269) + t212 * (-pkin(2) * t172 + pkin(6) * t145 - t270) - pkin(1) * t172 + qJ(2) * (-t144 * t211 + t145 * t212), t211 * (-pkin(6) * t135 - t92) + t212 * (-pkin(2) * t154 + pkin(6) * t136 + t93) - pkin(1) * t154 + qJ(2) * (-t135 * t211 + t136 * t212), -pkin(6) * t280 + t212 * (pkin(2) * t166 + pkin(6) * t93) + pkin(1) * t166 + qJ(2) * (t212 * t93 - t280), t211 * (t114 * t219 + t244) + t212 * (t114 * t216 - t243), t211 * (t153 * t216 + t219 * t88) + t212 * (-t153 * t219 + t216 * t88), t211 * (t110 * t219 + t121 * t216) + t212 * (t110 * t216 - t121 * t219), t211 * (t113 * t219 - t244) + t212 * (t113 * t216 + t243), t211 * (t111 * t219 - t117 * t216) + t212 * (t111 * t216 + t117 * t219), t211 * (t127 * t219 + t163 * t216) + t212 * (t127 * t216 - t163 * t219), t211 * (-pkin(6) * t77 - t216 * t53 + t219 * t58) + t212 * (-pkin(2) * t95 + pkin(6) * t78 + t216 * t58 + t219 * t53) - pkin(1) * t95 + qJ(2) * (-t211 * t77 + t212 * t78), t211 * (-pkin(6) * t79 - t216 * t55 + t219 * t61) + t212 * (-pkin(2) * t100 + pkin(6) * t80 + t216 * t61 + t219 * t55) - pkin(1) * t100 + qJ(2) * (-t211 * t79 + t212 * t80), t211 * (-pkin(6) * t59 + t219 * t30) + t212 * (pkin(6) * t60 + t216 * t30) + qJ(2) * (-t211 * t59 + t212 * t60) + (pkin(3) * t262 + t212 * t241 - pkin(1)) * t87, (t211 * (pkin(3) * t216 - pkin(7) * t219) + t212 * (-pkin(7) * t216 + t241) - pkin(1)) * t36 + t297 * (-t211 * (-t102 * t219 + t216 * t37) + t212 * (t102 * t216 + t219 * t37)), t211 * (t219 * t41 + t246) + t212 * (t216 * t41 - t245), t211 * (t115 * t216 + t21 * t219) + t212 * (-t115 * t219 + t21 * t216), t211 * (t216 * t75 + t219 * t51) + t212 * (t216 * t51 - t219 * t75), t211 * (t219 * t40 - t246) + t212 * (t216 * t40 + t245), t211 * (-t216 * t72 + t219 * t52) + t212 * (t216 * t52 + t219 * t72), t211 * (t162 * t216 + t219 * t56) + t212 * (-t162 * t219 + t216 * t56), t211 * (-pkin(6) * t27 + t11 * t219 - t15 * t216) + t212 * (-pkin(2) * t38 + pkin(6) * t28 + t11 * t216 + t15 * t219) - pkin(1) * t38 + qJ(2) * (-t211 * t27 + t212 * t28), t211 * (-pkin(6) * t31 + t12 * t219 - t16 * t216) + t212 * (-pkin(2) * t47 + pkin(6) * t32 + t12 * t216 + t16 * t219) - pkin(1) * t47 + qJ(2) * (-t211 * t31 + t212 * t32), t211 * (-pkin(6) * t18 - t17 * t216 + t2 * t219) + t212 * (-pkin(2) * t20 + pkin(6) * t19 + t17 * t219 + t2 * t216) - pkin(1) * t20 + qJ(2) * (-t18 * t211 + t19 * t212), t211 * (-pkin(6) * t4 + t1 * t219 - t216 * t3) + t212 * (-pkin(2) * t6 + pkin(6) * t5 + t1 * t216 + t219 * t3) - pkin(1) * t6 + qJ(2) * (-t211 * t4 + t212 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t250, t251, -t199, -t192, 0, 0, 0, 0, 0, 0, t170, t172, t154, -t166, 0, 0, 0, 0, 0, 0, t95, t100, t87, t36, 0, 0, 0, 0, 0, 0, t38, t47, t20, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, t194 - t193, t195, -t175, -t234, qJDD(3), -t125, -t126, 0, 0, t142 * t215 + t182 * t263, t118 * t215 + t120 * t218, t159 * t218 + t303, t180 * t264 - t218 * t228, t158 * t215 + t271, (-t180 * t215 - t182 * t218) * t191, pkin(3) * t118 + pkin(7) * t96 - t275, pkin(3) * t122 + pkin(7) * t101 + t276, pkin(3) * t137 + pkin(7) * t89 + t37, -pkin(3) * t102 + pkin(7) * t37, t215 * t70 + t218 * t69, t215 * t45 + t218 * t43, t215 * t85 + t218 * t83, t215 * t68 + t218 * t67, t215 * t86 + t218 * t84, t215 * t98 + t218 * t97, -pkin(3) * t71 + pkin(7) * t39 + t215 * t35 + t218 * t26, -pkin(3) * t296 + pkin(7) * t48 + t215 * t42 + t218 * t29, -pkin(3) * t94 + pkin(7) * t22 + t215 * t9 + t218 * t8, -pkin(3) * t57 + pkin(7) * t7 - pkin(8) * t282 + t10 * t218; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t155, t153, t121, -t155, -t117, t163, -t63, -t64, 0, 0, t116, t115, t75, -t116, -t72, t162, t230, t226, t286, t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t116, t115, t75, -t116, -t72, t162, -t24, -t25, 0, 0;];
tauJ_reg = t23;
