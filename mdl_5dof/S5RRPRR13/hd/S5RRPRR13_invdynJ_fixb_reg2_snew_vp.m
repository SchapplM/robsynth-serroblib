% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRR13
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRR13_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR13_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR13_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR13_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:33:55
% EndTime: 2019-12-31 20:34:08
% DurationCPUTime: 5.61s
% Computational Cost: add. (32049->441), mult. (70284->628), div. (0->0), fcn. (50247->10), ass. (0->269)
t246 = sin(qJ(5));
t244 = sin(pkin(9));
t245 = cos(pkin(9));
t248 = sin(qJ(2));
t279 = qJD(1) * t248;
t216 = -qJD(2) * t245 + t244 * t279;
t217 = qJD(2) * t244 + t245 * t279;
t247 = sin(qJ(4));
t251 = cos(qJ(4));
t191 = t216 * t251 + t217 * t247;
t193 = -t216 * t247 + t217 * t251;
t250 = cos(qJ(5));
t161 = t191 * t250 + t193 * t246;
t163 = -t191 * t246 + t193 * t250;
t118 = t163 * t161;
t275 = qJD(1) * qJD(2);
t233 = t248 * t275;
t252 = cos(qJ(2));
t274 = t252 * qJDD(1);
t222 = -t233 + t274;
t218 = -qJDD(4) + t222;
t214 = -qJDD(5) + t218;
t312 = -t118 - t214;
t319 = t246 * t312;
t318 = t250 * t312;
t254 = qJD(1) ^ 2;
t296 = t217 * t216;
t257 = -t222 - t296;
t317 = t244 * t257;
t316 = t245 * t257;
t297 = t193 * t191;
t255 = -t218 - t297;
t315 = t247 * t255;
t314 = t251 * t255;
t236 = t248 * qJDD(1);
t270 = t252 * t275;
t221 = t236 + t270;
t200 = qJDD(2) * t244 + t221 * t245;
t265 = -qJDD(2) * t245 + t221 * t244;
t266 = t247 * t200 + t251 * t265;
t147 = -qJD(4) * t193 - t266;
t148 = -qJD(4) * t191 + t200 * t251 - t247 * t265;
t102 = -qJD(5) * t161 + t147 * t246 + t148 * t250;
t278 = t252 * qJD(1);
t231 = -qJD(4) + t278;
t227 = -qJD(5) + t231;
t151 = t161 * t227;
t313 = t102 + t151;
t179 = t191 * t231;
t135 = -t179 + t148;
t311 = t179 + t148;
t205 = t216 * t278;
t183 = -t200 + t205;
t206 = t217 * t278;
t181 = -t265 - t206;
t268 = -t147 * t250 + t246 * t148;
t75 = (qJD(5) + t227) * t163 + t268;
t131 = (qJD(4) + t231) * t193 + t266;
t159 = t161 ^ 2;
t160 = t163 ^ 2;
t189 = t191 ^ 2;
t190 = t193 ^ 2;
t310 = t216 ^ 2;
t215 = t217 ^ 2;
t226 = t227 ^ 2;
t229 = t231 ^ 2;
t309 = qJD(2) ^ 2;
t249 = sin(qJ(1));
t253 = cos(qJ(1));
t269 = t249 * g(1) - t253 * g(2);
t210 = qJDD(1) * pkin(1) + t254 * pkin(6) + t269;
t261 = t221 + t270;
t169 = -t261 * qJ(3) + (-t222 + t233) * pkin(2) - t210;
t263 = g(1) * t253 + g(2) * t249;
t298 = qJDD(1) * pkin(6);
t211 = -pkin(1) * t254 - t263 + t298;
t262 = -pkin(2) * t252 - qJ(3) * t248;
t264 = t254 * t262 + t211;
t306 = t248 * g(3);
t173 = -pkin(2) * t309 + qJDD(2) * qJ(3) + t252 * t264 - t306;
t136 = 0.2e1 * qJD(3) * t217 - t169 * t245 + t244 * t173;
t106 = t257 * pkin(3) + pkin(7) * t183 - t136;
t137 = -0.2e1 * qJD(3) * t216 + t169 * t244 + t173 * t245;
t201 = -pkin(3) * t278 - pkin(7) * t217;
t108 = -pkin(3) * t310 - pkin(7) * t265 + t201 * t278 + t137;
t65 = -t106 * t251 + t247 * t108;
t52 = t255 * pkin(4) - pkin(8) * t135 - t65;
t174 = -pkin(4) * t231 - pkin(8) * t193;
t66 = t106 * t247 + t108 * t251;
t54 = -pkin(4) * t189 + pkin(8) * t147 + t174 * t231 + t66;
t29 = t246 * t54 - t250 * t52;
t30 = t246 * t52 + t250 * t54;
t16 = t246 * t30 - t250 * t29;
t308 = pkin(4) * t16;
t78 = t102 - t151;
t45 = -t246 * t75 - t250 * t78;
t307 = pkin(4) * t45;
t305 = t252 * g(3);
t38 = t247 * t66 - t251 * t65;
t304 = t244 * t38;
t303 = t245 * t38;
t172 = -qJDD(2) * pkin(2) - t309 * qJ(3) + t248 * t264 + qJDD(3) + t305;
t138 = pkin(3) * t265 - pkin(7) * t310 + t217 * t201 + t172;
t84 = -pkin(4) * t147 - pkin(8) * t189 + t193 * t174 + t138;
t302 = t246 * t84;
t301 = t247 * t16;
t300 = t250 * t84;
t299 = t251 * t16;
t295 = t227 * t246;
t294 = t227 * t250;
t293 = t231 * t247;
t292 = t231 * t251;
t291 = t244 * t172;
t184 = t222 - t296;
t290 = t244 * t184;
t289 = t245 * t172;
t288 = t245 * t184;
t111 = -t118 + t214;
t287 = t246 * t111;
t286 = t247 * t138;
t153 = t218 - t297;
t285 = t247 * t153;
t230 = t252 * t254 * t248;
t284 = t248 * (qJDD(2) + t230);
t283 = t250 * t111;
t282 = t251 * t138;
t281 = t251 * t153;
t280 = t252 * (-t230 + qJDD(2));
t273 = t252 * t118;
t272 = t252 * t297;
t271 = t252 * t296;
t17 = t246 * t29 + t250 * t30;
t39 = t247 * t65 + t251 * t66;
t96 = t136 * t244 + t137 * t245;
t197 = t211 * t248 + t305;
t198 = t211 * t252 - t306;
t267 = t248 * t197 + t198 * t252;
t260 = t136 * t245 - t137 * t244;
t116 = -t226 - t159;
t82 = t116 * t246 + t318;
t259 = pkin(4) * t82 - t29;
t258 = -pkin(1) + t262;
t140 = -t160 - t226;
t88 = t140 * t250 + t287;
t256 = pkin(4) * t88 - t30;
t242 = t252 ^ 2;
t241 = t248 ^ 2;
t238 = t242 * t254;
t237 = t241 * t254;
t223 = -0.2e1 * t233 + t274;
t220 = t236 + 0.2e1 * t270;
t212 = t252 * t222;
t204 = -t215 - t238;
t203 = -t215 + t238;
t202 = -t238 + t310;
t194 = -t238 - t310;
t182 = t200 + t205;
t180 = -t206 + t265;
t177 = -t215 - t310;
t176 = -t190 + t229;
t175 = t189 - t229;
t171 = -t190 - t229;
t166 = -t204 * t244 + t288;
t165 = t204 * t245 + t290;
t164 = t190 - t189;
t158 = -t229 - t189;
t157 = t194 * t245 - t317;
t156 = t194 * t244 + t316;
t150 = -t160 + t226;
t149 = t159 - t226;
t146 = t181 * t245 - t183 * t244;
t143 = (t191 * t251 - t193 * t247) * t231;
t142 = (t191 * t247 + t193 * t251) * t231;
t139 = -t189 - t190;
t130 = (qJD(4) - t231) * t193 + t266;
t128 = t175 * t251 + t285;
t127 = -t176 * t247 + t314;
t126 = t175 * t247 - t281;
t125 = t176 * t251 + t315;
t124 = t148 * t251 + t193 * t293;
t123 = t148 * t247 - t193 * t292;
t122 = -t147 * t247 - t191 * t292;
t121 = t147 * t251 - t191 * t293;
t120 = -t171 * t247 + t281;
t119 = t171 * t251 + t285;
t117 = t160 - t159;
t115 = t158 * t251 - t315;
t114 = t158 * t247 + t314;
t110 = (t161 * t250 - t163 * t246) * t227;
t109 = (t161 * t246 + t163 * t250) * t227;
t103 = -t159 - t160;
t101 = -qJD(5) * t163 - t268;
t100 = t149 * t250 + t287;
t99 = -t150 * t246 + t318;
t98 = t149 * t246 - t283;
t97 = t150 * t250 + t319;
t94 = -t131 * t251 + t135 * t247;
t93 = -t130 * t251 - t247 * t311;
t92 = -t131 * t247 - t135 * t251;
t91 = -t130 * t247 + t251 * t311;
t90 = -pkin(7) * t119 + t282;
t89 = -t140 * t246 + t283;
t87 = -t119 * t244 + t120 * t245;
t86 = t119 * t245 + t120 * t244;
t85 = -pkin(7) * t114 + t286;
t83 = t116 * t250 - t319;
t81 = -t114 * t244 + t115 * t245;
t80 = t114 * t245 + t115 * t244;
t74 = (qJD(5) - t227) * t163 + t268;
t73 = t102 * t250 + t163 * t295;
t72 = t102 * t246 - t163 * t294;
t71 = -t101 * t246 - t161 * t294;
t70 = t101 * t250 - t161 * t295;
t69 = -t109 * t247 + t110 * t251;
t68 = t109 * t251 + t110 * t247;
t67 = -pkin(3) * t311 + pkin(7) * t120 + t286;
t63 = -pkin(3) * t130 + pkin(7) * t115 - t282;
t62 = t100 * t251 - t247 * t98;
t61 = -t247 * t97 + t251 * t99;
t60 = t100 * t247 + t251 * t98;
t59 = t247 * t99 + t251 * t97;
t58 = -t244 * t92 + t245 * t94;
t57 = t244 * t94 + t245 * t92;
t56 = -t247 * t88 + t251 * t89;
t55 = t247 * t89 + t251 * t88;
t53 = -pkin(8) * t88 + t300;
t50 = -pkin(8) * t82 + t302;
t49 = -t247 * t82 + t251 * t83;
t48 = t247 * t83 + t251 * t82;
t47 = t246 * t78 - t250 * t75;
t46 = -t246 * t313 - t250 * t74;
t44 = -t246 * t74 + t250 * t313;
t43 = -t247 * t72 + t251 * t73;
t42 = -t247 * t70 + t251 * t71;
t41 = t247 * t73 + t251 * t72;
t40 = t247 * t71 + t251 * t70;
t37 = -pkin(4) * t313 + pkin(8) * t89 + t302;
t36 = -pkin(3) * t138 + pkin(7) * t39;
t35 = -pkin(4) * t74 + pkin(8) * t83 - t300;
t34 = -pkin(7) * t92 - t38;
t33 = -t244 * t55 + t245 * t56;
t32 = t244 * t56 + t245 * t55;
t31 = -pkin(3) * t139 + pkin(7) * t94 + t39;
t27 = -t244 * t48 + t245 * t49;
t26 = t244 * t49 + t245 * t48;
t25 = -t247 * t45 + t251 * t47;
t24 = -t247 * t44 + t251 * t46;
t23 = t247 * t47 + t251 * t45;
t22 = t247 * t46 + t251 * t44;
t21 = t245 * t39 - t304;
t20 = t244 * t39 + t303;
t19 = -pkin(7) * t55 - t247 * t37 + t251 * t53;
t18 = -pkin(7) * t48 - t247 * t35 + t251 * t50;
t15 = -pkin(3) * t313 + pkin(7) * t56 + t247 * t53 + t251 * t37;
t14 = -pkin(3) * t74 + pkin(7) * t49 + t247 * t50 + t251 * t35;
t13 = -pkin(4) * t84 + pkin(8) * t17;
t12 = -t23 * t244 + t245 * t25;
t11 = t23 * t245 + t244 * t25;
t10 = -pkin(8) * t45 - t16;
t9 = -pkin(4) * t103 + pkin(8) * t47 + t17;
t8 = t17 * t251 - t301;
t7 = t17 * t247 + t299;
t6 = -pkin(7) * t23 + t10 * t251 - t247 * t9;
t5 = -pkin(3) * t103 + pkin(7) * t25 + t10 * t247 + t251 * t9;
t4 = -t244 * t7 + t245 * t8;
t3 = t244 * t8 + t245 * t7;
t2 = -pkin(7) * t7 - pkin(8) * t299 - t13 * t247;
t1 = -pkin(3) * t84 + pkin(7) * t8 - pkin(8) * t301 + t13 * t251;
t28 = [0, 0, 0, 0, 0, qJDD(1), t269, t263, 0, 0, t261 * t248, t220 * t252 + t223 * t248, t284 + t252 * (-t237 + t309), -t248 * t270 + t212, t248 * (t238 - t309) + t280, 0, t252 * t210 + pkin(1) * t223 + pkin(6) * (t252 * (-t238 - t309) - t284), -t248 * t210 - pkin(1) * t220 + pkin(6) * (-t280 - t248 * (-t237 - t309)), pkin(1) * (t237 + t238) + (t241 + t242) * t298 + t267, pkin(1) * t210 + pkin(6) * t267, t248 * (t200 * t245 + t206 * t244) - t271, t248 * (-t180 * t245 - t182 * t244) + t252 * (-t215 + t310), t248 * (-t203 * t244 + t316) + t252 * t183, t248 * (-t205 * t245 + t244 * t265) + t271, t248 * (t202 * t245 + t290) - t252 * t181, t212 + t248 * (t216 * t245 - t217 * t244) * t278, t248 * (-qJ(3) * t156 + t291) + t252 * (-pkin(2) * t156 + t136) - pkin(1) * t156 + pkin(6) * (t157 * t252 + t180 * t248), t248 * (-qJ(3) * t165 + t289) + t252 * (-pkin(2) * t165 + t137) - pkin(1) * t165 + pkin(6) * (t166 * t252 + t182 * t248), t248 * t260 + pkin(6) * (t146 * t252 + t177 * t248) + t258 * (t181 * t244 + t183 * t245), pkin(6) * (t172 * t248 + t252 * t96) - t258 * t260, t248 * (-t123 * t244 + t124 * t245) - t272, t248 * (-t244 * t91 + t245 * t93) - t252 * t164, t248 * (-t125 * t244 + t127 * t245) - t252 * t135, t248 * (-t121 * t244 + t122 * t245) + t272, t248 * (-t126 * t244 + t128 * t245) + t252 * t131, t248 * (-t142 * t244 + t143 * t245) + t252 * t218, t248 * (-qJ(3) * t80 - t244 * t63 + t245 * t85) + t252 * (-pkin(2) * t80 - pkin(3) * t114 + t65) - pkin(1) * t80 + pkin(6) * (t130 * t248 + t252 * t81), t248 * (-qJ(3) * t86 - t244 * t67 + t245 * t90) + t252 * (-pkin(2) * t86 - pkin(3) * t119 + t66) - pkin(1) * t86 + pkin(6) * (t248 * t311 + t252 * t87), t248 * (-qJ(3) * t57 - t244 * t31 + t245 * t34) + t252 * (-pkin(2) * t57 - pkin(3) * t92) - pkin(1) * t57 + pkin(6) * (t139 * t248 + t252 * t58), t248 * (-pkin(7) * t303 - qJ(3) * t20 - t244 * t36) + t252 * (-pkin(2) * t20 - pkin(3) * t38) - pkin(1) * t20 + pkin(6) * (t138 * t248 + t21 * t252), t248 * (-t244 * t41 + t245 * t43) - t273, t248 * (-t22 * t244 + t24 * t245) - t252 * t117, t248 * (-t244 * t59 + t245 * t61) - t252 * t78, t248 * (-t244 * t40 + t245 * t42) + t273, t248 * (-t244 * t60 + t245 * t62) + t252 * t75, t248 * (-t244 * t68 + t245 * t69) + t252 * t214, t248 * (-qJ(3) * t26 - t14 * t244 + t18 * t245) + t252 * (-pkin(2) * t26 - pkin(3) * t48 - t259) - pkin(1) * t26 + pkin(6) * (t248 * t74 + t252 * t27), t248 * (-qJ(3) * t32 - t15 * t244 + t19 * t245) + t252 * (-pkin(2) * t32 - pkin(3) * t55 - t256) - pkin(1) * t32 + pkin(6) * (t248 * t313 + t252 * t33), t248 * (-qJ(3) * t11 - t244 * t5 + t245 * t6) + t252 * (-pkin(2) * t11 - pkin(3) * t23 - t307) - pkin(1) * t11 + pkin(6) * (t103 * t248 + t12 * t252), t248 * (-qJ(3) * t3 - t1 * t244 + t2 * t245) + t252 * (-pkin(2) * t3 - pkin(3) * t7 - t308) - pkin(1) * t3 + pkin(6) * (t248 * t84 + t252 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t230, t237 - t238, t236, t230, t274, qJDD(2), -t197, -t198, 0, 0, t200 * t244 - t206 * t245, -t180 * t244 + t182 * t245, t203 * t245 + t317, -t205 * t244 - t245 * t265, t202 * t244 - t288, (t216 * t244 + t217 * t245) * t278, -pkin(2) * t180 + qJ(3) * t157 - t289, -pkin(2) * t182 + qJ(3) * t166 + t291, -pkin(2) * t177 + qJ(3) * t146 + t96, -pkin(2) * t172 + qJ(3) * t96, t123 * t245 + t124 * t244, t244 * t93 + t245 * t91, t125 * t245 + t127 * t244, t121 * t245 + t122 * t244, t126 * t245 + t128 * t244, t142 * t245 + t143 * t244, -pkin(2) * t130 + qJ(3) * t81 + t244 * t85 + t245 * t63, -pkin(2) * t311 + qJ(3) * t87 + t244 * t90 + t245 * t67, -pkin(2) * t139 + qJ(3) * t58 + t244 * t34 + t245 * t31, -pkin(2) * t138 - pkin(7) * t304 + qJ(3) * t21 + t245 * t36, t244 * t43 + t245 * t41, t22 * t245 + t24 * t244, t244 * t61 + t245 * t59, t244 * t42 + t245 * t40, t244 * t62 + t245 * t60, t244 * t69 + t245 * t68, -pkin(2) * t74 + qJ(3) * t27 + t14 * t245 + t18 * t244, -pkin(2) * t313 + qJ(3) * t33 + t15 * t245 + t19 * t244, -pkin(2) * t103 + qJ(3) * t12 + t244 * t6 + t245 * t5, -pkin(2) * t84 + qJ(3) * t4 + t1 * t245 + t2 * t244; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180, t182, t177, t172, 0, 0, 0, 0, 0, 0, t130, t311, t139, t138, 0, 0, 0, 0, 0, 0, t74, t313, t103, t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t297, t164, t135, -t297, -t131, -t218, -t65, -t66, 0, 0, t118, t117, t78, -t118, -t75, -t214, t259, t256, t307, t308; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t117, t78, -t118, -t75, -t214, -t29, -t30, 0, 0;];
tauJ_reg = t28;
