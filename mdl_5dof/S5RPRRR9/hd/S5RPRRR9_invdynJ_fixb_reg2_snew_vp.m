% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRR9
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRR9_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR9_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR9_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:08:02
% EndTime: 2019-12-31 19:08:14
% DurationCPUTime: 5.43s
% Computational Cost: add. (22120->423), mult. (53831->617), div. (0->0), fcn. (41796->10), ass. (0->268)
t229 = sin(pkin(9));
t226 = t229 ^ 2;
t230 = cos(pkin(9));
t227 = t230 ^ 2;
t278 = t226 + t227;
t231 = sin(qJ(5));
t233 = sin(qJ(3));
t237 = cos(qJ(3));
t210 = (-t229 * t233 + t230 * t237) * qJD(1);
t249 = t229 * t237 + t230 * t233;
t211 = t249 * qJD(1);
t232 = sin(qJ(4));
t236 = cos(qJ(4));
t185 = t232 * t210 + t236 * t211;
t270 = t230 * qJDD(1);
t271 = t229 * qJDD(1);
t168 = -t233 * t271 + t237 * t270;
t276 = t211 * qJD(3);
t194 = t168 - t276;
t209 = t249 * qJDD(1);
t277 = t210 * qJD(3);
t196 = t209 + t277;
t255 = -t236 * t194 + t232 * t196;
t143 = -t185 * qJD(4) - t255;
t142 = qJDD(5) - t143;
t228 = qJD(3) + qJD(4);
t235 = cos(qJ(5));
t171 = t231 * t185 - t235 * t228;
t173 = t235 * t185 + t231 * t228;
t146 = t173 * t171;
t316 = t142 - t146;
t323 = t231 * t316;
t183 = -t236 * t210 + t232 * t211;
t161 = t185 * t183;
t225 = qJDD(3) + qJDD(4);
t315 = -t161 + t225;
t322 = t232 * t315;
t198 = t210 * t211;
t314 = qJDD(3) + t198;
t321 = t233 * t314;
t320 = t235 * t316;
t319 = t236 * t315;
t318 = t237 * t314;
t239 = qJD(1) ^ 2;
t234 = sin(qJ(1));
t307 = cos(qJ(1));
t260 = t234 * g(1) - t307 * g(2);
t253 = -qJDD(2) + t260;
t302 = t230 * pkin(2);
t189 = (t278 * pkin(6) + qJ(2)) * t239 + (pkin(1) + t302) * qJDD(1) + t253;
t207 = t210 ^ 2;
t252 = qJD(3) * pkin(3) - t211 * pkin(7);
t140 = t194 * pkin(3) + t207 * pkin(7) - t211 * t252 + t189;
t248 = t307 * g(1) + t234 * g(2);
t212 = -t239 * pkin(1) + qJDD(1) * qJ(2) - t248;
t272 = qJD(1) * qJD(2);
t317 = t212 + 0.2e1 * t272;
t159 = t183 * pkin(4) - t185 * pkin(8);
t308 = t228 ^ 2;
t279 = t317 * t230;
t245 = -t227 * t239 * pkin(2) + pkin(6) * t270 + t279;
t247 = -pkin(6) * qJDD(1) + t239 * t302 - t212;
t301 = t230 * g(3);
t259 = -0.2e1 * t229 * t272 - t301;
t148 = -(t233 * g(3) + t237 * t247) * t229 + t233 * t245 - t237 * t259;
t311 = (-t196 + t277) * pkin(7);
t242 = pkin(3) * t314 - t148 + t311;
t149 = t237 * t245 + t233 * t259 + (-t237 * g(3) + t233 * t247) * t229;
t119 = -t207 * pkin(3) + t194 * pkin(7) - qJD(3) * t252 + t149;
t285 = t236 * t119;
t75 = t232 * t242 + t285;
t61 = -t308 * pkin(4) + t225 * pkin(8) - t183 * t159 + t75;
t250 = t232 * t194 + t236 * t196;
t144 = -t183 * qJD(4) + t250;
t178 = t228 * t183;
t128 = t144 - t178;
t64 = -t128 * pkin(8) + (t228 * t185 - t143) * pkin(4) - t140;
t31 = t231 * t61 - t235 * t64;
t32 = t231 * t64 + t235 * t61;
t16 = t231 * t31 + t235 * t32;
t280 = t239 * qJ(2);
t297 = qJDD(1) * pkin(1);
t206 = t253 + t280 + t297;
t313 = t278 * t280 - t206 - t297;
t180 = qJD(5) + t183;
t256 = t231 * t144 - t235 * t225;
t100 = (qJD(5) - t180) * t173 + t256;
t169 = t171 ^ 2;
t170 = t173 ^ 2;
t179 = t180 ^ 2;
t181 = t183 ^ 2;
t182 = t185 ^ 2;
t208 = t211 ^ 2;
t306 = pkin(4) * t232;
t303 = t229 * g(3);
t74 = t232 * t119 - t236 * t242;
t60 = -t225 * pkin(4) - t308 * pkin(8) + t185 * t159 + t74;
t300 = -pkin(4) * t60 + pkin(8) * t16;
t57 = t231 * t60;
t36 = t232 * t75 - t236 * t74;
t299 = t233 * t36;
t58 = t235 * t60;
t298 = t237 * t36;
t296 = t180 * t231;
t295 = t180 * t235;
t294 = t228 * t232;
t293 = t228 * t236;
t110 = -t237 * t148 + t233 * t149;
t292 = t229 * t110;
t108 = t142 + t146;
t291 = t231 * t108;
t290 = t232 * t140;
t157 = t161 + t225;
t289 = t232 * t157;
t288 = t233 * t189;
t191 = qJDD(3) - t198;
t287 = t233 * t191;
t286 = t235 * t108;
t284 = t236 * t140;
t283 = t236 * t157;
t282 = t237 * t189;
t281 = t237 * t191;
t275 = qJD(4) + t228;
t273 = qJD(5) + t180;
t251 = -t235 * t144 - t231 * t225;
t105 = t273 * t171 + t251;
t139 = -t170 - t179;
t78 = -t231 * t139 - t286;
t269 = pkin(4) * t105 + pkin(8) * t78 + t57;
t101 = -t273 * t173 - t256;
t132 = -t179 - t169;
t73 = t235 * t132 - t323;
t268 = pkin(4) * t101 + pkin(8) * t73 - t58;
t264 = t232 * t146;
t263 = t236 * t146;
t261 = -pkin(4) * t236 - pkin(3);
t37 = t232 * t74 + t236 * t75;
t131 = t169 + t170;
t116 = -t171 * qJD(5) - t251;
t154 = t180 * t171;
t104 = t116 + t154;
t56 = -t100 * t235 + t231 * t104;
t258 = pkin(4) * t131 + pkin(8) * t56 + t16;
t111 = t233 * t148 + t237 * t149;
t254 = t229 * (t317 * t229 + t301) + t230 * (t279 - t303);
t15 = t231 * t32 - t235 * t31;
t246 = (-qJD(4) + t228) * t185 - t255;
t238 = qJD(3) ^ 2;
t222 = t227 * qJDD(1);
t221 = t226 * qJDD(1);
t213 = t278 * t239;
t202 = -t208 - t238;
t201 = -t208 + t238;
t200 = t207 - t238;
t195 = t209 + 0.2e1 * t277;
t193 = -t168 + 0.2e1 * t276;
t188 = -t238 - t207;
t177 = -t182 + t308;
t176 = t181 - t308;
t175 = -t182 - t308;
t174 = -t207 - t208;
t167 = -t233 * t202 - t281;
t166 = t237 * t202 - t287;
t165 = t237 * t168 + t233 * t209;
t164 = t233 * t168 - t237 * t209;
t163 = t237 * t188 - t321;
t162 = t233 * t188 + t318;
t160 = t182 - t181;
t155 = -t308 - t181;
t153 = -t170 + t179;
t152 = t169 - t179;
t151 = (-t183 * t236 + t185 * t232) * t228;
t150 = (-t183 * t232 - t185 * t236) * t228;
t145 = t170 - t169;
t141 = -t181 - t182;
t138 = t236 * t176 - t289;
t137 = -t232 * t177 + t319;
t136 = t232 * t176 + t283;
t135 = t236 * t177 + t322;
t134 = -t232 * t175 - t283;
t133 = t236 * t175 - t289;
t129 = t144 + t178;
t127 = -t275 * t183 + t250;
t124 = t275 * t185 + t255;
t123 = t236 * t144 - t185 * t294;
t122 = t232 * t144 + t185 * t293;
t121 = -t232 * t143 + t183 * t293;
t120 = t236 * t143 + t183 * t294;
t118 = t236 * t155 - t322;
t117 = t232 * t155 + t319;
t115 = -t173 * qJD(5) - t256;
t114 = (-t171 * t235 + t173 * t231) * t180;
t113 = (-t171 * t231 - t173 * t235) * t180;
t106 = -pkin(7) * t133 - t284;
t103 = t116 - t154;
t97 = t235 * t116 - t173 * t296;
t96 = t231 * t116 + t173 * t295;
t95 = -t231 * t115 + t171 * t295;
t94 = t235 * t115 + t171 * t296;
t93 = -t233 * t133 + t237 * t134;
t92 = t237 * t133 + t233 * t134;
t91 = -pkin(7) * t117 - t290;
t90 = t236 * t114 + t232 * t142;
t89 = t232 * t114 - t236 * t142;
t88 = t235 * t152 - t291;
t87 = -t231 * t153 + t320;
t86 = t231 * t152 + t286;
t85 = t235 * t153 + t323;
t84 = t232 * t129 + t236 * t246;
t83 = -t236 * t124 - t232 * t128;
t82 = -t236 * t129 + t232 * t246;
t81 = -t232 * t124 + t236 * t128;
t80 = -t233 * t117 + t237 * t118;
t79 = t237 * t117 + t233 * t118;
t77 = t235 * t139 - t291;
t72 = t231 * t132 + t320;
t69 = t236 * t97 + t264;
t68 = t236 * t95 - t264;
t67 = t232 * t97 - t263;
t66 = t232 * t95 + t263;
t65 = -pkin(3) * t127 + pkin(7) * t134 - t290;
t62 = -pkin(3) * t124 + pkin(7) * t118 + t284;
t55 = t235 * t101 - t231 * t103;
t54 = -t100 * t231 - t235 * t104;
t53 = t231 * t101 + t235 * t103;
t51 = -t232 * t100 + t236 * t88;
t50 = t232 * t104 + t236 * t87;
t49 = t236 * t100 + t232 * t88;
t48 = -t236 * t104 + t232 * t87;
t47 = -t232 * t105 + t236 * t78;
t46 = t236 * t105 + t232 * t78;
t45 = -t232 * t101 + t236 * t73;
t44 = t236 * t101 + t232 * t73;
t43 = t232 * t145 + t236 * t55;
t42 = -t236 * t145 + t232 * t55;
t41 = -t233 * t82 + t237 * t84;
t40 = t233 * t84 + t237 * t82;
t39 = -t232 * t131 + t236 * t56;
t38 = t236 * t131 + t232 * t56;
t35 = -pkin(8) * t77 + t58;
t34 = pkin(3) * t140 + pkin(7) * t37;
t33 = -pkin(8) * t72 + t57;
t28 = -pkin(7) * t82 - t36;
t27 = -pkin(3) * t141 + pkin(7) * t84 + t37;
t26 = -t233 * t46 + t237 * t47;
t25 = t233 * t47 + t237 * t46;
t24 = -t233 * t44 + t237 * t45;
t23 = t233 * t45 + t237 * t44;
t22 = -pkin(4) * t77 + t32;
t21 = -pkin(4) * t72 + t31;
t20 = -t233 * t38 + t237 * t39;
t19 = t233 * t39 + t237 * t38;
t18 = t237 * t37 - t299;
t17 = t233 * t37 + t298;
t13 = t236 * t16 + t232 * t60;
t12 = t232 * t16 - t236 * t60;
t11 = -pkin(8) * t54 - t15;
t10 = -pkin(7) * t46 - t232 * t22 + t236 * t35;
t9 = -pkin(7) * t44 - t232 * t21 + t236 * t33;
t8 = -pkin(3) * t77 + pkin(7) * t47 + t236 * t22 + t232 * t35;
t7 = -pkin(3) * t72 + pkin(7) * t45 + t236 * t21 + t232 * t33;
t6 = -pkin(7) * t38 + t236 * t11 + t54 * t306;
t5 = pkin(7) * t39 + t232 * t11 + t261 * t54;
t4 = -t233 * t12 + t237 * t13;
t3 = t237 * t12 + t233 * t13;
t2 = -pkin(7) * t12 + (-pkin(8) * t236 + t306) * t15;
t1 = pkin(7) * t13 + (-pkin(8) * t232 + t261) * t15;
t14 = [0, 0, 0, 0, 0, qJDD(1), t260, t248, 0, 0, t221, 0.2e1 * t229 * t270, 0, t222, 0, 0, -t313 * t230, t313 * t229, pkin(1) * t213 + qJ(2) * (t222 + t221) + t254, pkin(1) * t206 + qJ(2) * t254, t229 * (t237 * t196 - t233 * t276) + t230 * (t233 * t196 + t237 * t276), t229 * (-t237 * t193 - t233 * t195) + t230 * (-t233 * t193 + t237 * t195), t229 * (-t233 * t201 + t318) + t230 * (t237 * t201 + t321), t229 * (-t233 * t194 - t237 * t277) + t230 * (t237 * t194 - t233 * t277), t229 * (t237 * t200 - t287) + t230 * (t233 * t200 + t281), (t229 * (t210 * t237 + t211 * t233) + t230 * (t210 * t233 - t211 * t237)) * qJD(3), t229 * (-pkin(6) * t162 - t288) + t230 * (-pkin(2) * t193 + pkin(6) * t163 + t282) - pkin(1) * t193 + qJ(2) * (-t229 * t162 + t230 * t163), t229 * (-pkin(6) * t166 - t282) + t230 * (-pkin(2) * t195 + pkin(6) * t167 - t288) - pkin(1) * t195 + qJ(2) * (-t229 * t166 + t230 * t167), t229 * (-pkin(6) * t164 - t110) + t230 * (-pkin(2) * t174 + pkin(6) * t165 + t111) - pkin(1) * t174 + qJ(2) * (-t229 * t164 + t230 * t165), -pkin(6) * t292 + t230 * (pkin(2) * t189 + pkin(6) * t111) + pkin(1) * t189 + qJ(2) * (t230 * t111 - t292), t229 * (-t233 * t122 + t237 * t123) + t230 * (t237 * t122 + t233 * t123), t229 * (-t233 * t81 + t237 * t83) + t230 * (t233 * t83 + t237 * t81), t229 * (-t233 * t135 + t237 * t137) + t230 * (t237 * t135 + t233 * t137), t229 * (-t233 * t120 + t237 * t121) + t230 * (t237 * t120 + t233 * t121), t229 * (-t233 * t136 + t237 * t138) + t230 * (t237 * t136 + t233 * t138), t229 * (-t233 * t150 + t237 * t151) + t230 * (t237 * t150 + t233 * t151), t229 * (-pkin(6) * t79 - t233 * t62 + t237 * t91) + t230 * (-pkin(2) * t124 + pkin(6) * t80 + t233 * t91 + t237 * t62) - pkin(1) * t124 + qJ(2) * (-t229 * t79 + t230 * t80), t229 * (-pkin(6) * t92 + t237 * t106 - t233 * t65) + t230 * (-pkin(2) * t127 + pkin(6) * t93 + t233 * t106 + t237 * t65) - pkin(1) * t127 + qJ(2) * (-t229 * t92 + t230 * t93), t229 * (-pkin(6) * t40 - t233 * t27 + t237 * t28) + t230 * (-pkin(2) * t141 + pkin(6) * t41 + t233 * t28 + t237 * t27) - pkin(1) * t141 + qJ(2) * (-t229 * t40 + t230 * t41), t229 * (-pkin(6) * t17 - pkin(7) * t298 - t233 * t34) + t230 * (pkin(2) * t140 + pkin(6) * t18 - pkin(7) * t299 + t237 * t34) + pkin(1) * t140 + qJ(2) * (-t229 * t17 + t230 * t18), t229 * (-t233 * t67 + t237 * t69) + t230 * (t233 * t69 + t237 * t67), t229 * (-t233 * t42 + t237 * t43) + t230 * (t233 * t43 + t237 * t42), t229 * (-t233 * t48 + t237 * t50) + t230 * (t233 * t50 + t237 * t48), t229 * (-t233 * t66 + t237 * t68) + t230 * (t233 * t68 + t237 * t66), t229 * (-t233 * t49 + t237 * t51) + t230 * (t233 * t51 + t237 * t49), t229 * (-t233 * t89 + t237 * t90) + t230 * (t233 * t90 + t237 * t89), t229 * (-pkin(6) * t23 - t233 * t7 + t237 * t9) + t230 * (-pkin(2) * t72 + pkin(6) * t24 + t233 * t9 + t237 * t7) - pkin(1) * t72 + qJ(2) * (-t229 * t23 + t230 * t24), t229 * (-pkin(6) * t25 + t237 * t10 - t233 * t8) + t230 * (-pkin(2) * t77 + pkin(6) * t26 + t233 * t10 + t237 * t8) - pkin(1) * t77 + qJ(2) * (-t229 * t25 + t230 * t26), t229 * (-pkin(6) * t19 - t233 * t5 + t237 * t6) + t230 * (-pkin(2) * t54 + pkin(6) * t20 + t233 * t6 + t237 * t5) - pkin(1) * t54 + qJ(2) * (-t229 * t19 + t230 * t20), t229 * (-pkin(6) * t3 - t233 * t1 + t237 * t2) + t230 * (-pkin(2) * t15 + pkin(6) * t4 + t237 * t1 + t233 * t2) - pkin(1) * t15 + qJ(2) * (-t229 * t3 + t230 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t270, t271, -t213, -t206, 0, 0, 0, 0, 0, 0, t193, t195, t174, -t189, 0, 0, 0, 0, 0, 0, t124, t127, t141, -t140, 0, 0, 0, 0, 0, 0, t72, t77, t54, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t198, t208 - t207, t209, t198, t168, qJDD(3), -t148, -t149, 0, 0, t161, t160, t129, -t161, t246, t225, pkin(3) * t117 - t74, -t285 - t232 * (-t233 * (t245 - t303) + t237 * (t229 * t247 + t259) + t311) + (-t232 * t314 + t133) * pkin(3), pkin(3) * t82, pkin(3) * t36, t96, t53, t85, t94, t86, t113, pkin(3) * t44 + t268, pkin(3) * t46 + t269, pkin(3) * t38 + t258, pkin(3) * t12 + t300; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t161, t160, t129, -t161, t246, t225, -t74, -t75, 0, 0, t96, t53, t85, t94, t86, t113, t268, t269, t258, t300; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t146, t145, t104, -t146, -t100, t142, -t31, -t32, 0, 0;];
tauJ_reg = t14;
