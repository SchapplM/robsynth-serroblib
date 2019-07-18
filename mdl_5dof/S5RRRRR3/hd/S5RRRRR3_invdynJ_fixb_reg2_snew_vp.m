% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRR3
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:19
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRR3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_invdynJ_fixb_reg2_snew_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:18:44
% EndTime: 2019-07-18 17:18:58
% DurationCPUTime: 3.56s
% Computational Cost: add. (17209->344), mult. (35206->504), div. (0->0), fcn. (27294->10), ass. (0->239)
t221 = sin(qJ(5));
t226 = cos(qJ(5));
t223 = sin(qJ(3));
t224 = sin(qJ(2));
t228 = cos(qJ(3));
t229 = cos(qJ(2));
t201 = (t223 * t229 + t224 * t228) * qJD(1);
t218 = qJD(2) + qJD(3);
t222 = sin(qJ(4));
t227 = cos(qJ(4));
t183 = t201 * t222 - t218 * t227;
t185 = t201 * t227 + t218 * t222;
t165 = t185 * t183;
t215 = t224 * qJDD(1);
t259 = qJD(1) * qJD(2);
t249 = t229 * t259;
t205 = t215 + t249;
t216 = t229 * qJDD(1);
t250 = t224 * t259;
t206 = t216 - t250;
t243 = t223 * t205 - t206 * t228;
t168 = -qJD(3) * t201 - t243;
t167 = qJDD(4) - t168;
t287 = -t165 + t167;
t268 = t223 * t224;
t199 = (-t228 * t229 + t268) * qJD(1);
t169 = -qJD(3) * t199 + t205 * t228 + t206 * t223;
t189 = t218 * t199;
t154 = -t189 + t169;
t225 = sin(qJ(1));
t283 = cos(qJ(1));
t207 = t225 * g(1) - g(2) * t283;
t291 = t206 - t250;
t180 = pkin(1) * t291 + t207;
t273 = t218 * t201;
t108 = -t154 * pkin(5) + (-t168 + t273) * pkin(2) - t180;
t237 = g(1) * t283 + g(2) * t225;
t194 = t229 * g(3) - t224 * t237;
t231 = qJD(1) ^ 2;
t209 = t229 * t231 * t224;
t241 = qJDD(2) + t209;
t232 = pkin(1) * t241 - t194;
t195 = -t224 * g(3) - t229 * t237;
t230 = qJD(2) ^ 2;
t220 = t229 ^ 2;
t272 = t220 * t231;
t186 = (-t230 - t272) * pkin(1) + t195;
t264 = t228 * t186;
t157 = t223 * t232 + t264;
t175 = pkin(2) * t199 - pkin(5) * t201;
t258 = qJDD(2) + qJDD(3);
t286 = t218 ^ 2;
t133 = -pkin(2) * t286 + pkin(5) * t258 - t199 * t175 + t157;
t77 = -t108 * t227 + t222 * t133;
t55 = pkin(3) * t287 - t77;
t181 = t183 ^ 2;
t196 = qJD(4) + t199;
t193 = t196 ^ 2;
t151 = -t181 - t193;
t78 = t108 * t222 + t133 * t227;
t56 = pkin(3) * t151 + t78;
t25 = t221 * t56 - t226 * t55;
t26 = t221 * t55 + t226 * t56;
t13 = t221 * t26 - t226 * t25;
t297 = t222 * t13;
t296 = t227 * t13;
t161 = t183 * t226 + t185 * t221;
t163 = -t183 * t221 + t185 * t226;
t135 = t163 * t161;
t166 = qJDD(5) + t167;
t288 = -t135 + t166;
t295 = t221 * t288;
t294 = t222 * t287;
t293 = t226 * t288;
t292 = t227 * t287;
t40 = t222 * t77 + t227 * t78;
t191 = qJD(5) + t196;
t149 = t191 * t161;
t233 = -t227 * t169 - t222 * t258;
t144 = -t183 * qJD(4) - t233;
t244 = t222 * t169 - t227 * t258;
t235 = qJD(4) * t185 + t244;
t94 = -qJD(5) * t161 + t144 * t226 - t221 * t235;
t290 = -t149 + t94;
t178 = t201 * t199;
t173 = -t178 + t258;
t289 = t228 * t173;
t245 = t221 * t144 + t226 * t235;
t66 = (qJD(5) - t191) * t163 + t245;
t122 = (qJD(4) - t196) * t185 + t244;
t159 = t161 ^ 2;
t160 = t163 ^ 2;
t182 = t185 ^ 2;
t190 = t191 ^ 2;
t197 = t199 ^ 2;
t198 = t201 ^ 2;
t285 = pkin(3) * t13;
t69 = t149 + t94;
t32 = -t221 * t66 - t226 * t69;
t284 = pkin(3) * t32;
t282 = pkin(2) * t223;
t156 = t223 * t186 - t228 * t232;
t132 = -pkin(2) * t258 - pkin(5) * t286 + t175 * t201 + t156;
t92 = t132 + (t185 * t196 + t235) * pkin(3);
t281 = t222 * t92;
t280 = t226 * t92;
t279 = t227 * t92;
t278 = -pkin(2) * t132 + pkin(5) * t40;
t277 = t191 * t221;
t276 = t191 * t226;
t275 = t196 * t222;
t274 = t196 * t227;
t106 = t135 + t166;
t271 = t221 * t106;
t128 = t222 * t132;
t137 = t165 + t167;
t270 = t222 * t137;
t269 = t223 * t199;
t267 = t226 * t106;
t129 = t227 * t132;
t266 = t227 * t137;
t265 = t228 * t180;
t263 = t228 * t201;
t261 = qJD(4) + t196;
t123 = -t185 * t261 - t244;
t98 = t151 * t227 - t294;
t257 = pkin(2) * t123 + pkin(5) * t98 - t129;
t256 = t223 * t135;
t255 = t223 * t165;
t254 = t228 * t135;
t253 = t228 * t165;
t158 = -t182 - t193;
t104 = -t158 * t222 - t266;
t127 = t183 * t261 + t233;
t252 = pkin(2) * t127 + pkin(5) * t104 + t128;
t251 = -pkin(2) * t228 - pkin(1);
t15 = t221 * t25 + t226 * t26;
t99 = -t159 - t160;
t10 = -pkin(3) * t99 + t15;
t34 = t221 * t69 - t226 * t66;
t22 = -t222 * t32 + t227 * t34;
t248 = -pkin(2) * t99 + pkin(5) * t22 + t10 * t227 - t297;
t43 = -pkin(3) * t290 + t221 * t92;
t139 = -t160 - t190;
t79 = t139 * t226 - t271;
t80 = -t139 * t221 - t267;
t47 = -t222 * t79 + t227 * t80;
t247 = -pkin(2) * t290 + pkin(5) * t47 + t222 * t280 + t227 * t43;
t146 = t181 + t182;
t172 = t196 * t183;
t126 = t172 + t144;
t85 = -t122 * t227 + t126 * t222;
t246 = pkin(2) * t146 + pkin(5) * t85 + t40;
t6 = t15 * t227 - t297;
t242 = -pkin(2) * t92 - pkin(3) * t279 + pkin(5) * t6;
t240 = t222 * t78 - t227 * t77;
t239 = t156 * t228 - t157 * t223;
t115 = -t190 - t159;
t71 = t115 * t221 + t293;
t238 = pkin(3) * t71 - t25;
t72 = t115 * t226 - t295;
t37 = -t222 * t71 + t227 * t72;
t65 = (qJD(5) + t191) * t163 + t245;
t44 = -pkin(3) * t65 - t280;
t236 = -pkin(2) * t65 + pkin(5) * t37 + t221 * t281 + t227 * t44;
t234 = pkin(3) * t79 - t26;
t153 = (-qJD(3) + t218) * t201 - t243;
t219 = t224 ^ 2;
t188 = -t198 + t286;
t187 = t197 - t286;
t176 = t198 - t197;
t174 = t178 + t258;
t171 = -t182 + t193;
t170 = t181 - t193;
t164 = t182 - t181;
t155 = t189 + t169;
t152 = (qJD(3) + t218) * t201 + t243;
t148 = -t160 + t190;
t147 = t159 - t190;
t141 = (-t183 * t227 + t185 * t222) * t196;
t140 = (-t183 * t222 - t185 * t227) * t196;
t134 = t160 - t159;
t125 = -t172 + t144;
t119 = t144 * t227 - t185 * t275;
t118 = t144 * t222 + t185 * t274;
t117 = t183 * t274 + t222 * t235;
t116 = t183 * t275 - t227 * t235;
t114 = t170 * t227 - t270;
t113 = -t171 * t222 + t292;
t112 = t170 * t222 + t266;
t111 = t171 * t227 + t294;
t110 = (-t161 * t226 + t163 * t221) * t191;
t109 = (-t161 * t221 - t163 * t226) * t191;
t103 = t158 * t227 - t270;
t97 = t151 * t222 + t292;
t93 = -qJD(5) * t163 - t245;
t90 = t147 * t226 - t271;
t89 = -t148 * t221 + t293;
t88 = t147 * t221 + t267;
t87 = t148 * t226 + t295;
t84 = t123 * t227 - t125 * t222;
t83 = -t122 * t222 - t126 * t227;
t82 = t123 * t222 + t125 * t227;
t76 = -pkin(5) * t103 + t129;
t73 = -pkin(5) * t97 + t128;
t62 = -t163 * t277 + t226 * t94;
t61 = t163 * t276 + t221 * t94;
t60 = t161 * t276 - t221 * t93;
t59 = t161 * t277 + t226 * t93;
t58 = -t109 * t222 + t110 * t227;
t57 = t109 * t227 + t110 * t222;
t53 = -pkin(2) * t103 + t78;
t52 = -pkin(2) * t97 + t77;
t51 = -t222 * t88 + t227 * t90;
t50 = -t222 * t87 + t227 * t89;
t49 = t222 * t90 + t227 * t88;
t48 = t222 * t89 + t227 * t87;
t46 = t222 * t80 + t227 * t79;
t36 = t222 * t72 + t227 * t71;
t33 = -t221 * t290 - t226 * t65;
t31 = -t221 * t65 + t226 * t290;
t30 = -t222 * t61 + t227 * t62;
t29 = -t222 * t59 + t227 * t60;
t28 = t222 * t62 + t227 * t61;
t27 = t222 * t60 + t227 * t59;
t23 = -pkin(5) * t83 - t240;
t21 = -t222 * t31 + t227 * t33;
t20 = t222 * t34 + t227 * t32;
t19 = t222 * t33 + t227 * t31;
t17 = -pkin(5) * t46 - t222 * t43 + t226 * t279;
t16 = -pkin(5) * t36 + t221 * t279 - t222 * t44;
t11 = -pkin(2) * t46 - t234;
t8 = -pkin(2) * t36 - t238;
t7 = -pkin(2) * t20 - t284;
t5 = t15 * t222 + t296;
t3 = pkin(3) * t281 - pkin(5) * t5;
t2 = -pkin(5) * t20 - t10 * t222 - t296;
t1 = -pkin(2) * t5 - t285;
t4 = [0, 0, 0, 0, 0, qJDD(1), t207, t237, 0, 0, (t205 + t249) * t224, t224 * (t216 - 0.2e1 * t250) + t229 * (t215 + 0.2e1 * t249), t224 * t241 + t229 * (-t219 * t231 + t230), t291 * t229, t224 * (-t230 + t272) + t229 * (qJDD(2) - t209), 0, t229 * t207, -t224 * t207, t194 * t224 + t195 * t229, 0, t224 * (t169 * t228 - t223 * t273) + t229 * (t169 * t223 + t218 * t263), t224 * (-t152 * t228 - t154 * t223) + t229 * (-t152 * t223 + t154 * t228), t224 * (-t188 * t223 + t289) + t229 * (t173 * t223 + t188 * t228), t224 * (-t168 * t223 + t189 * t228) + t229 * (t168 * t228 + t218 * t269), t224 * (-t223 * t174 + t187 * t228) + t229 * (t174 * t228 + t187 * t223), (t224 * (-t199 * t228 + t201 * t223) + t229 * (-t263 - t269)) * t218, -t180 * t268 + t229 * (-pkin(1) * t152 + t265), -t224 * t265 + t229 * (-pkin(1) * t154 - t223 * t180), t224 * t239 + t229 * (t223 * t156 + t228 * t157 - pkin(1) * (-t197 - t198)), t229 * pkin(1) * t180, t224 * (t119 * t228 + t255) + t229 * (t119 * t223 - t253), t224 * (t164 * t223 + t228 * t84) + t229 * (-t164 * t228 + t223 * t84), t224 * (t113 * t228 + t126 * t223) + t229 * (t113 * t223 - t126 * t228), t224 * (t117 * t228 - t255) + t229 * (t117 * t223 + t253), t224 * (t114 * t228 - t122 * t223) + t229 * (t114 * t223 + t122 * t228), t224 * (t141 * t228 + t167 * t223) + t229 * (t141 * t223 - t167 * t228), t224 * (-t223 * t52 + t228 * t73) + t229 * (-pkin(1) * t97 + t223 * t73 + t228 * t52), t224 * (-t223 * t53 + t228 * t76) + t229 * (-pkin(1) * t103 + t223 * t76 + t228 * t53), t224 * (t228 * t23 + t282 * t83) + t229 * (t223 * t23 + t251 * t83), (t224 * (-pkin(5) * t228 + t282) + t229 * (-pkin(5) * t223 + t251)) * t240, t224 * (t228 * t30 + t256) + t229 * (t223 * t30 - t254), t224 * (t134 * t223 + t21 * t228) + t229 * (-t134 * t228 + t21 * t223), t224 * (t223 * t69 + t228 * t50) + t229 * (t223 * t50 - t228 * t69), t224 * (t228 * t29 - t256) + t229 * (t223 * t29 + t254), t224 * (-t223 * t66 + t228 * t51) + t229 * (t223 * t51 + t228 * t66), t224 * (t166 * t223 + t228 * t58) + t229 * (-t166 * t228 + t223 * t58), t224 * (t16 * t228 - t223 * t8) + t229 * (-pkin(1) * t36 + t16 * t223 + t228 * t8), t224 * (-t11 * t223 + t17 * t228) + t229 * (-pkin(1) * t46 + t11 * t228 + t17 * t223), t224 * (t2 * t228 - t223 * t7) + t229 * (-pkin(1) * t20 + t2 * t223 + t228 * t7), t224 * (-t1 * t223 + t228 * t3) + t229 * (-pkin(1) * t5 + t1 * t228 + t223 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t209, (t219 - t220) * t231, t215, t209, t216, qJDD(2), -t194, -t195, 0, 0, t178, t176, t155, -t178, t153, t258, pkin(1) * (t223 * (-t286 - t197) + t289) - t156, -t264 + t223 * t194 + (t228 * (-t198 - t286) + (-t174 - t241) * t223) * pkin(1), pkin(1) * (t153 * t223 - t228 * t155), -pkin(1) * t239, t118, t82, t111, t116, t112, t140, pkin(1) * (t123 * t228 + t223 * t98) + t257, pkin(1) * (t104 * t223 + t127 * t228) + t252, pkin(1) * (t146 * t228 + t223 * t85) + t246, pkin(1) * (-t132 * t228 + t223 * t40) + t278, t28, t19, t48, t27, t49, t57, pkin(1) * (t223 * t37 - t228 * t65) + t236, pkin(1) * (t223 * t47 - t228 * t290) + t247, pkin(1) * (t22 * t223 - t228 * t99) + t248, pkin(1) * (t223 * t6 - t228 * t92) + t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178, t176, t155, -t178, t153, t258, -t156, -t157, 0, 0, t118, t82, t111, t116, t112, t140, t257, t252, t246, t278, t28, t19, t48, t27, t49, t57, t236, t247, t248, t242; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, t164, t126, -t165, -t122, t167, -t77, -t78, 0, 0, t135, t134, t69, -t135, -t66, t166, t238, t234, t284, t285; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, t134, t69, -t135, -t66, t166, -t25, -t26, 0, 0;];
tauJ_reg  = t4;
