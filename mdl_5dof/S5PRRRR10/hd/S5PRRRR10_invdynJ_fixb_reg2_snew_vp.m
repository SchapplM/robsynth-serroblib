% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
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
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRR10_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR10_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR10_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PRRRR10_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:26:11
% EndTime: 2019-12-05 17:26:23
% DurationCPUTime: 3.90s
% Computational Cost: add. (17469->374), mult. (36401->583), div. (0->0), fcn. (29370->14), ass. (0->256)
t194 = sin(qJ(5));
t193 = cos(pkin(6));
t187 = qJD(2) * t193 + qJD(3);
t195 = sin(qJ(4));
t199 = cos(qJ(4));
t196 = sin(qJ(3));
t191 = sin(pkin(6));
t250 = qJD(2) * t191;
t236 = t196 * t250;
t165 = t187 * t195 + t199 * t236;
t200 = cos(qJ(3));
t249 = qJD(2) * t200;
t235 = t191 * t249;
t182 = -qJD(4) + t235;
t198 = cos(qJ(5));
t146 = t165 * t194 + t182 * t198;
t148 = t165 * t198 - t182 * t194;
t120 = t148 * t146;
t243 = qJDD(2) * t196;
t171 = (qJD(3) * t249 + t243) * t191;
t186 = qJDD(2) * t193 + qJDD(3);
t233 = t195 * t171 - t186 * t199;
t132 = -qJD(4) * t165 - t233;
t131 = qJDD(5) - t132;
t280 = -t120 + t131;
t287 = t194 * t280;
t286 = t198 * t280;
t197 = sin(qJ(2));
t201 = cos(qJ(2));
t192 = sin(pkin(5));
t252 = -g(3) + qJDD(1);
t264 = sin(pkin(11));
t265 = cos(pkin(11));
t211 = g(1) * t264 - g(2) * t265;
t266 = cos(pkin(5));
t281 = t266 * t211;
t207 = t192 * t252 + t281;
t212 = -g(1) * t265 - g(2) * t264;
t142 = -t197 * t212 + t201 * t207;
t202 = qJD(2) ^ 2;
t274 = pkin(8) * t191;
t205 = qJDD(2) * pkin(2) + t202 * t274 + t142;
t208 = -t192 * t211 + t252 * t266;
t285 = t191 * t208 + t193 * t205;
t188 = t191 ^ 2;
t284 = t188 * (qJD(2) * t187 - t193 * t202);
t244 = qJDD(2) * t191;
t229 = -qJD(3) * t236 + t200 * t244;
t167 = -qJDD(4) + t229;
t163 = -t187 * t199 + t195 * t236;
t261 = t165 * t163;
t209 = -t167 - t261;
t283 = t195 * t209;
t282 = t199 * t209;
t216 = -t199 * t171 - t195 * t186;
t133 = -qJD(4) * t163 - t216;
t155 = t163 * t182;
t114 = t133 + t155;
t176 = t187 * t235;
t279 = t171 + t176;
t175 = t187 * t236;
t278 = t175 - t229;
t159 = qJD(5) + t163;
t234 = t194 * t133 + t167 * t198;
t76 = (qJD(5) - t159) * t148 + t234;
t111 = (qJD(4) + t182) * t165 + t233;
t143 = t197 * t207 + t201 * t212;
t137 = -t202 * pkin(2) + pkin(8) * t244 + t143;
t94 = t200 * t137 + t196 * t285;
t144 = t146 ^ 2;
t145 = t148 ^ 2;
t157 = t159 ^ 2;
t161 = t163 ^ 2;
t162 = t165 ^ 2;
t277 = t182 ^ 2;
t276 = t187 ^ 2;
t275 = pkin(4) * t195;
t273 = t200 * pkin(3);
t121 = t191 * t205 - t193 * t208;
t203 = pkin(3) * t278 - pkin(9) * t279 - t121;
t231 = -pkin(9) * t196 - t273;
t170 = t231 * t250;
t87 = -t276 * pkin(3) + t186 * pkin(9) + t170 * t235 + t94;
t55 = t195 * t203 + t199 * t87;
t138 = pkin(4) * t163 - pkin(10) * t165;
t54 = t195 * t87 - t199 * t203;
t37 = pkin(4) * t167 - pkin(10) * t277 + t138 * t165 + t54;
t272 = t194 * t37;
t91 = t120 + t131;
t271 = t194 * t91;
t251 = t285 * t200;
t86 = -t186 * pkin(3) - t276 * pkin(9) + (t170 * t250 + t137) * t196 - t251;
t270 = t195 * t86;
t269 = t198 * t37;
t268 = t198 * t91;
t267 = t199 * t86;
t263 = t159 * t194;
t262 = t159 * t198;
t260 = t182 * t195;
t259 = t182 * t199;
t258 = t188 * t202;
t127 = t167 - t261;
t257 = t195 * t127;
t181 = t196 * t200 * t258;
t168 = t181 + t186;
t256 = t196 * t168;
t255 = t199 * t127;
t169 = -t181 + t186;
t253 = t200 * t169;
t248 = qJD(4) - t182;
t245 = qJD(5) + t159;
t189 = t196 ^ 2;
t242 = t189 * t258;
t190 = t200 ^ 2;
t241 = t190 * t258;
t240 = t195 * t120;
t239 = t199 * t120;
t238 = t200 * t261;
t237 = -pkin(4) * t199 - pkin(3);
t38 = -pkin(4) * t277 - pkin(10) * t167 - t138 * t163 + t55;
t51 = -t114 * pkin(10) + (-t165 * t182 - t132) * pkin(4) + t86;
t23 = t194 * t38 - t198 * t51;
t24 = t194 * t51 + t198 * t38;
t11 = t194 * t23 + t198 * t24;
t28 = t195 * t54 + t199 * t55;
t232 = -pkin(4) * t37 + pkin(10) * t11;
t10 = t194 * t24 - t198 * t23;
t8 = t11 * t199 + t195 * t37;
t230 = -t10 * t200 + t196 * t8;
t27 = t195 * t55 - t199 * t54;
t228 = t196 * t28 - t200 * t86;
t101 = t144 + t145;
t220 = -t198 * t133 + t194 * t167;
t100 = -qJD(5) * t146 - t220;
t130 = t159 * t146;
t80 = t100 + t130;
t47 = t194 * t80 - t198 * t76;
t32 = -t101 * t195 + t199 * t47;
t45 = -t194 * t76 - t198 * t80;
t227 = t196 * t32 - t200 * t45;
t106 = -t157 - t144;
t62 = t106 * t198 - t287;
t78 = -t148 * t245 - t234;
t36 = -t195 * t78 + t199 * t62;
t61 = t106 * t194 + t286;
t226 = t196 * t36 - t200 * t61;
t117 = -t145 - t157;
t64 = -t117 * t194 - t268;
t81 = t146 * t245 + t220;
t40 = -t195 * t81 + t199 * t64;
t63 = t117 * t198 - t271;
t225 = t196 * t40 - t200 * t63;
t93 = t196 * t137 - t251;
t224 = t196 * t94 - t200 * t93;
t59 = t196 * t93 + t200 * t94;
t112 = -t165 * t248 - t233;
t134 = -t277 - t161;
t98 = t134 * t199 - t283;
t223 = t112 * t200 + t196 * t98;
t124 = t161 + t162;
t115 = t133 - t155;
t84 = -t111 * t199 + t115 * t195;
t222 = t124 * t200 + t196 * t84;
t140 = -t162 - t277;
t103 = -t140 * t195 + t255;
t116 = t163 * t248 + t216;
t221 = t103 * t196 + t116 * t200;
t152 = -t176 + t171;
t153 = t175 + t229;
t219 = -t152 * t200 + t153 * t196;
t158 = -t242 - t276;
t218 = t158 * t200 - t169 * t196;
t172 = -t241 - t276;
t217 = t168 * t200 + t172 * t196;
t214 = pkin(4) * t81 + pkin(10) * t64 + t272;
t213 = pkin(4) * t78 + pkin(10) * t62 - t269;
t210 = pkin(4) * t101 + pkin(10) * t47 + t11;
t174 = (-t189 - t190) * t258;
t173 = (t189 - t190) * t258;
t151 = (t243 + (qJD(3) + t187) * t249) * t191;
t150 = -t162 + t277;
t149 = t161 - t277;
t141 = t172 * t200 - t256;
t139 = t162 - t161;
t136 = -t158 * t196 - t253;
t126 = -t145 + t157;
t125 = t144 - t157;
t123 = (t163 * t195 + t165 * t199) * t182;
t122 = t152 * t196 + t153 * t200;
t119 = t145 - t144;
t118 = -t191 * t278 + t193 * t217;
t110 = -t191 * t151 + t193 * t218;
t109 = t133 * t195 - t165 * t259;
t108 = t132 * t199 - t163 * t260;
t107 = -t191 * t174 + t193 * t219;
t105 = t149 * t195 - t255;
t104 = t150 * t199 + t283;
t102 = t140 * t199 + t257;
t99 = -qJD(5) * t148 - t234;
t97 = t134 * t195 + t282;
t96 = (-t146 * t198 + t148 * t194) * t159;
t95 = (-t146 * t194 - t148 * t198) * t159;
t83 = -t111 * t195 - t115 * t199;
t82 = t112 * t195 + t114 * t199;
t79 = t100 - t130;
t75 = t100 * t198 - t148 * t263;
t74 = t100 * t194 + t148 * t262;
t73 = t146 * t262 - t194 * t99;
t72 = -t146 * t263 - t198 * t99;
t71 = -t131 * t199 + t195 * t96;
t70 = t125 * t198 - t271;
t69 = -t126 * t194 + t286;
t68 = t125 * t194 + t268;
t67 = t126 * t198 + t287;
t66 = t103 * t200 - t116 * t196;
t65 = -t112 * t196 + t200 * t98;
t60 = -t124 * t196 + t200 * t84;
t58 = t195 * t75 - t239;
t57 = t195 * t73 + t239;
t56 = -t191 * t102 + t193 * t221;
t52 = -t191 * t97 + t193 * t223;
t49 = t191 * t121 + t193 * t224;
t48 = -t194 * t79 + t198 * t78;
t46 = t194 * t78 + t198 * t79;
t44 = pkin(3) * t116 + pkin(9) * t103 + t270;
t43 = pkin(3) * t112 + pkin(9) * t98 - t267;
t42 = t195 * t70 + t199 * t76;
t41 = t195 * t69 - t199 * t80;
t39 = t195 * t64 + t199 * t81;
t35 = t195 * t62 + t199 * t78;
t34 = -t191 * t83 + t193 * t222;
t33 = -t119 * t199 + t195 * t48;
t31 = t101 * t199 + t195 * t47;
t30 = t196 * t63 + t200 * t40;
t29 = -pkin(10) * t63 + t269;
t26 = -pkin(10) * t61 + t272;
t25 = t196 * t61 + t200 * t36;
t21 = -pkin(3) * t86 + pkin(9) * t28;
t20 = t196 * t86 + t200 * t28;
t19 = t196 * t45 + t200 * t32;
t18 = pkin(3) * t124 + pkin(9) * t84 + t28;
t17 = -pkin(4) * t63 + t24;
t16 = -pkin(4) * t61 + t23;
t15 = -t191 * t39 + t193 * t225;
t14 = -t191 * t35 + t193 * t226;
t13 = -t191 * t31 + t193 * t227;
t12 = -t191 * t27 + t193 * t228;
t9 = -pkin(10) * t45 - t10;
t7 = t11 * t195 - t199 * t37;
t6 = -pkin(3) * t63 + pkin(9) * t40 + t17 * t199 + t195 * t29;
t5 = -pkin(3) * t61 + pkin(9) * t36 + t16 * t199 + t195 * t26;
t4 = pkin(9) * t32 + t195 * t9 + t237 * t45;
t3 = t10 * t196 + t200 * t8;
t2 = -t191 * t7 + t193 * t230;
t1 = pkin(9) * t8 + (-pkin(10) * t195 + t237) * t10;
t22 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t252, 0, 0, 0, 0, 0, 0, (qJDD(2) * t201 - t197 * t202) * t192, (-qJDD(2) * t197 - t201 * t202) * t192, 0, t266 ^ 2 * t252 + (t201 * t142 + t197 * t143 - t281) * t192, 0, 0, 0, 0, 0, 0, t266 * (t191 * t217 + t193 * t278) + (t118 * t201 + t141 * t197) * t192, t266 * (t193 * t151 + t191 * t218) + (t110 * t201 + t136 * t197) * t192, t266 * (t193 * t174 + t191 * t219) + (t107 * t201 + t122 * t197) * t192, t266 * (-t193 * t121 + t191 * t224) + (t197 * t59 + t201 * t49) * t192, 0, 0, 0, 0, 0, 0, t266 * (t191 * t223 + t193 * t97) + (t197 * t65 + t201 * t52) * t192, t266 * (t193 * t102 + t191 * t221) + (t197 * t66 + t201 * t56) * t192, t266 * (t191 * t222 + t193 * t83) + (t197 * t60 + t201 * t34) * t192, t266 * (t191 * t228 + t193 * t27) + (t12 * t201 + t197 * t20) * t192, 0, 0, 0, 0, 0, 0, t266 * (t191 * t226 + t193 * t35) + (t14 * t201 + t197 * t25) * t192, t266 * (t191 * t225 + t193 * t39) + (t15 * t201 + t197 * t30) * t192, t266 * (t191 * t227 + t193 * t31) + (t13 * t201 + t19 * t197) * t192, t266 * (t191 * t230 + t193 * t7) + (t197 * t3 + t2 * t201) * t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t142, -t143, 0, 0, (t171 * t191 + t200 * t284) * t196, t193 * t173 + (-t196 * t278 + t200 * t279) * t191, t193 * t152 + (t256 + t200 * (-t242 + t276)) * t191, (t191 * t229 - t196 * t284) * t200, t193 * t153 + (t196 * (t241 - t276) + t253) * t191, t193 * t186, pkin(2) * t118 - t193 * t93 + (pkin(8) * t141 + t121 * t200) * t191, pkin(2) * t110 - t193 * t94 + (pkin(8) * t136 - t121 * t196) * t191, pkin(2) * t107 + (pkin(8) * t122 + t59) * t191, pkin(2) * t49 + t274 * t59, t193 * t109 + (t196 * (t133 * t199 + t165 * t260) - t238) * t191, t193 * t82 + (t196 * (t112 * t199 - t114 * t195) - t200 * t139) * t191, t193 * t104 + (t196 * (-t150 * t195 + t282) - t200 * t115) * t191, t193 * t108 + (t196 * (-t132 * t195 - t163 * t259) + t238) * t191, t193 * t105 + (t196 * (t149 * t199 + t257) + t200 * t111) * t191, t193 * t123 + (t200 * t167 + t196 * (t163 * t199 - t165 * t195) * t182) * t191, pkin(2) * t52 + t193 * t43 + (t196 * (-pkin(9) * t97 + t270) + t200 * (-pkin(3) * t97 + t54) + pkin(8) * t65) * t191, pkin(2) * t56 + t193 * t44 + (t196 * (-pkin(9) * t102 + t267) + t200 * (-pkin(3) * t102 + t55) + pkin(8) * t66) * t191, pkin(2) * t34 + t193 * t18 + (t196 * (-pkin(9) * t83 - t27) - t83 * t273 + pkin(8) * t60) * t191, pkin(2) * t12 + t193 * t21 + (pkin(8) * t20 + t231 * t27) * t191, t193 * t58 + (t196 * (t199 * t75 + t240) - t200 * t74) * t191, t193 * t33 + (t196 * (t119 * t195 + t199 * t48) - t200 * t46) * t191, t193 * t41 + (t196 * (t195 * t80 + t199 * t69) - t200 * t67) * t191, t193 * t57 + (t196 * (t199 * t73 - t240) + t200 * t72) * t191, t193 * t42 + (t196 * (-t195 * t76 + t199 * t70) - t200 * t68) * t191, t193 * t71 + (t196 * (t131 * t195 + t199 * t96) - t200 * t95) * t191, pkin(2) * t14 + t193 * t5 + (t196 * (-pkin(9) * t35 - t16 * t195 + t199 * t26) + t200 * (-pkin(3) * t35 - t213) + pkin(8) * t25) * t191, pkin(2) * t15 + t193 * t6 + (t196 * (-pkin(9) * t39 - t17 * t195 + t199 * t29) + t200 * (-pkin(3) * t39 - t214) + pkin(8) * t30) * t191, pkin(2) * t13 + t193 * t4 + (t196 * (-pkin(9) * t31 + t199 * t9 + t275 * t45) + t200 * (-pkin(3) * t31 - t210) + pkin(8) * t19) * t191, pkin(2) * t2 + t193 * t1 + (t196 * (-pkin(9) * t7 + (-pkin(10) * t199 + t275) * t10) + t200 * (-pkin(3) * t7 - t232) + pkin(8) * t3) * t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181, t173, t152, t181, t153, t186, -t93, -t94, 0, 0, t109, t82, t104, t108, t105, t123, t43, t44, t18, t21, t58, t33, t41, t57, t42, t71, t5, t6, t4, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t261, t139, t115, -t261, -t111, -t167, -t54, -t55, 0, 0, t74, t46, t67, -t72, t68, t95, t213, t214, t210, t232; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t120, t119, t80, -t120, -t76, t131, -t23, -t24, 0, 0;];
tauJ_reg = t22;
