% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:22:30
% EndTime: 2019-12-31 22:22:43
% DurationCPUTime: 6.10s
% Computational Cost: add. (31527->452), mult. (68970->643), div. (0->0), fcn. (51050->10), ass. (0->287)
t260 = cos(qJ(2));
t244 = t260 * qJDD(1);
t255 = sin(qJ(2));
t297 = qJD(1) * qJD(2);
t287 = t255 * t297;
t231 = t244 - t287;
t251 = t260 ^ 2;
t262 = qJD(1) ^ 2;
t256 = sin(qJ(1));
t329 = cos(qJ(1));
t286 = g(1) * t256 - t329 * g(2);
t269 = qJDD(1) * pkin(1) + t286;
t302 = qJD(1) * t255;
t270 = qJD(2) * pkin(2) - pkin(7) * t302;
t195 = pkin(2) * t231 - t270 * t302 + t269 + (pkin(7) * t251 + pkin(6)) * t262;
t252 = sin(qJ(5));
t254 = sin(qJ(3));
t259 = cos(qJ(3));
t223 = qJD(1) * t259 * t260 - t254 * t302;
t224 = (t260 * t254 + t255 * t259) * qJD(1);
t253 = sin(qJ(4));
t258 = cos(qJ(4));
t206 = t223 * t253 + t224 * t258;
t243 = t255 * qJDD(1);
t288 = t260 * t297;
t230 = t243 + t288;
t279 = t230 * t254 - t259 * t231;
t192 = -qJD(3) * t224 - t279;
t272 = t230 * t259 + t231 * t254;
t193 = qJD(3) * t223 + t272;
t281 = -t258 * t192 + t193 * t253;
t141 = -qJD(4) * t206 - t281;
t140 = qJDD(5) - t141;
t249 = qJD(2) + qJD(3);
t242 = qJD(4) + t249;
t257 = cos(qJ(5));
t188 = t206 * t252 - t257 * t242;
t190 = t206 * t257 + t242 * t252;
t159 = t190 * t188;
t333 = t140 - t159;
t339 = t252 * t333;
t204 = -t258 * t223 + t224 * t253;
t173 = t206 * t204;
t296 = qJDD(2) + qJDD(3);
t241 = qJDD(4) + t296;
t332 = -t173 + t241;
t338 = t253 * t332;
t212 = t223 * t224;
t331 = t212 + t296;
t337 = t254 * t331;
t336 = t257 * t333;
t335 = t258 * t332;
t334 = t259 * t331;
t221 = t223 ^ 2;
t276 = pkin(3) * t249 - pkin(8) * t224;
t138 = pkin(3) * t192 + pkin(8) * t221 - t224 * t276 + t195;
t169 = pkin(4) * t204 - pkin(9) * t206;
t330 = t242 ^ 2;
t304 = t255 * t262;
t271 = t329 * g(1) + t256 * g(2);
t324 = qJDD(1) * pkin(6);
t226 = -t262 * pkin(1) - t271 + t324;
t310 = t226 * t255;
t183 = qJDD(2) * pkin(2) - pkin(7) * t230 - t310 + (pkin(2) * t304 + pkin(7) * t297 - g(3)) * t260;
t216 = -t255 * g(3) + t260 * t226;
t246 = t251 * t262;
t184 = -pkin(2) * t246 + t231 * pkin(7) - qJD(2) * t270 + t216;
t157 = t254 * t183 + t259 * t184;
t119 = -t221 * pkin(3) + t192 * pkin(8) - t249 * t276 + t157;
t156 = -t259 * t183 + t254 * t184;
t219 = t249 * t223;
t179 = t193 - t219;
t264 = pkin(3) * t331 - t179 * pkin(8) - t156;
t73 = t258 * t119 + t253 * t264;
t66 = -t330 * pkin(4) + t241 * pkin(9) - t204 * t169 + t73;
t273 = t192 * t253 + t193 * t258;
t142 = -qJD(4) * t204 + t273;
t198 = t242 * t204;
t129 = t142 - t198;
t68 = -t129 * pkin(9) + (t206 * t242 - t141) * pkin(4) - t138;
t32 = t252 * t66 - t257 * t68;
t33 = t252 * t68 + t257 * t66;
t17 = t252 * t32 + t257 * t33;
t200 = qJD(5) + t204;
t282 = t142 * t252 - t257 * t241;
t104 = (qJD(5) - t200) * t190 + t282;
t186 = t188 ^ 2;
t187 = t190 ^ 2;
t199 = t200 ^ 2;
t202 = t204 ^ 2;
t203 = t206 ^ 2;
t222 = t224 ^ 2;
t248 = t249 ^ 2;
t328 = pkin(4) * t253;
t72 = t119 * t253 - t258 * t264;
t65 = -t241 * pkin(4) - t330 * pkin(9) + t169 * t206 + t72;
t327 = -pkin(4) * t65 + pkin(9) * t17;
t62 = t252 * t65;
t38 = t253 * t73 - t258 * t72;
t326 = t254 * t38;
t63 = t257 * t65;
t325 = t259 * t38;
t114 = t140 + t159;
t323 = t114 * t252;
t322 = t114 * t257;
t116 = -t156 * t259 + t157 * t254;
t321 = t116 * t255;
t320 = t138 * t253;
t319 = t138 * t258;
t167 = t173 + t241;
t318 = t167 * t253;
t317 = t167 * t258;
t316 = t195 * t254;
t315 = t195 * t259;
t314 = t200 * t252;
t313 = t200 * t257;
t209 = -t212 + t296;
t312 = t209 * t254;
t311 = t209 * t259;
t309 = t242 * t253;
t308 = t242 * t258;
t307 = t249 * t254;
t306 = t249 * t259;
t238 = t260 * t304;
t305 = t255 * (qJDD(2) + t238);
t303 = t260 * (qJDD(2) - t238);
t301 = qJD(3) + t249;
t300 = qJD(4) + t242;
t298 = qJD(5) + t200;
t13 = t17 * t253 - t258 * t65;
t295 = pkin(3) * t13 + t327;
t274 = -t142 * t257 - t241 * t252;
t109 = t298 * t188 + t274;
t153 = -t187 - t199;
t83 = -t153 * t252 - t322;
t294 = pkin(4) * t109 + pkin(9) * t83 + t62;
t105 = -t298 * t190 - t282;
t145 = -t199 - t186;
t80 = t145 * t257 - t339;
t293 = pkin(4) * t105 + pkin(9) * t80 - t63;
t292 = t253 * t159;
t291 = t258 * t159;
t290 = -pkin(4) * t258 - pkin(3);
t39 = t253 * t72 + t258 * t73;
t139 = t186 + t187;
t121 = -qJD(5) * t188 - t274;
t164 = t200 * t188;
t108 = t121 + t164;
t61 = -t104 * t257 + t108 * t252;
t285 = pkin(4) * t139 + pkin(9) * t61 + t17;
t51 = t109 * t258 + t253 * t83;
t284 = pkin(3) * t51 + t294;
t48 = t105 * t258 + t253 * t80;
t283 = pkin(3) * t48 + t293;
t117 = t156 * t254 + t259 * t157;
t215 = g(3) * t260 + t310;
t280 = t255 * t215 + t260 * t216;
t43 = t139 * t258 + t253 * t61;
t278 = pkin(3) * t43 + t285;
t165 = -t330 - t202;
t135 = t165 * t253 + t335;
t275 = pkin(3) * t135 - t72;
t16 = t252 * t33 - t257 * t32;
t268 = (-qJD(4) + t242) * t206 - t281;
t267 = (-qJD(3) + t249) * t224 - t279;
t191 = -t203 - t330;
t147 = t191 * t258 - t318;
t263 = pkin(3) * t147 - t73;
t261 = qJD(2) ^ 2;
t250 = t255 ^ 2;
t245 = t250 * t262;
t232 = t244 - 0.2e1 * t287;
t229 = t243 + 0.2e1 * t288;
t225 = pkin(6) * t262 + t269;
t218 = -t222 + t248;
t217 = t221 - t248;
t214 = -t222 - t248;
t211 = t222 - t221;
t207 = -t248 - t221;
t197 = -t203 + t330;
t196 = t202 - t330;
t194 = -t221 - t222;
t181 = -t214 * t254 - t311;
t180 = t214 * t259 - t312;
t178 = t193 + t219;
t177 = t301 * t223 + t272;
t174 = t301 * t224 + t279;
t172 = t203 - t202;
t171 = t207 * t259 - t337;
t170 = t207 * t254 + t334;
t163 = -t187 + t199;
t162 = t186 - t199;
t161 = (-t204 * t258 + t206 * t253) * t242;
t160 = (-t204 * t253 - t206 * t258) * t242;
t158 = t187 - t186;
t155 = -t202 - t203;
t152 = t196 * t258 - t318;
t151 = -t197 * t253 + t335;
t150 = t196 * t253 + t317;
t149 = t197 * t258 + t338;
t148 = -t191 * t253 - t317;
t144 = t179 * t254 + t259 * t267;
t143 = -t179 * t259 + t254 * t267;
t136 = t165 * t258 - t338;
t133 = (-t188 * t257 + t190 * t252) * t200;
t132 = (-t188 * t252 - t190 * t257) * t200;
t131 = -t300 * t204 + t273;
t130 = t142 + t198;
t126 = t300 * t206 + t281;
t125 = t142 * t258 - t206 * t309;
t124 = t142 * t253 + t206 * t308;
t123 = -t141 * t253 + t204 * t308;
t122 = t141 * t258 + t204 * t309;
t120 = -qJD(5) * t190 - t282;
t112 = -t147 * t254 + t148 * t259;
t111 = t147 * t259 + t148 * t254;
t110 = -pkin(8) * t147 - t319;
t107 = t121 - t164;
t101 = t121 * t257 - t190 * t314;
t100 = t121 * t252 + t190 * t313;
t99 = -t120 * t252 + t188 * t313;
t98 = t120 * t257 + t188 * t314;
t97 = -pkin(8) * t135 - t320;
t96 = t133 * t258 + t140 * t253;
t95 = t133 * t253 - t140 * t258;
t94 = t162 * t257 - t323;
t93 = -t163 * t252 + t336;
t92 = t162 * t252 + t322;
t91 = t163 * t257 + t339;
t90 = -t135 * t254 + t136 * t259;
t89 = t135 * t259 + t136 * t254;
t88 = t130 * t253 + t258 * t268;
t87 = -t126 * t258 - t129 * t253;
t86 = -t130 * t258 + t253 * t268;
t85 = -t126 * t253 + t129 * t258;
t84 = pkin(3) * t86;
t82 = t153 * t257 - t323;
t79 = t145 * t252 + t336;
t77 = t101 * t258 + t292;
t76 = t258 * t99 - t292;
t75 = t101 * t253 - t291;
t74 = t253 * t99 + t291;
t70 = -pkin(3) * t131 + pkin(8) * t148 - t320;
t69 = -pkin(3) * t126 + pkin(8) * t136 + t319;
t60 = t105 * t257 - t107 * t252;
t59 = -t104 * t252 - t108 * t257;
t58 = t105 * t252 + t107 * t257;
t56 = -t104 * t253 + t258 * t94;
t55 = t108 * t253 + t258 * t93;
t54 = t104 * t258 + t253 * t94;
t53 = -t108 * t258 + t253 * t93;
t52 = -t109 * t253 + t258 * t83;
t49 = -t105 * t253 + t258 * t80;
t46 = t158 * t253 + t258 * t60;
t45 = -t158 * t258 + t253 * t60;
t44 = -t139 * t253 + t258 * t61;
t41 = -t254 * t86 + t259 * t88;
t40 = t254 * t88 + t259 * t86;
t37 = pkin(3) * t38;
t36 = -pkin(9) * t82 + t63;
t35 = -pkin(9) * t79 + t62;
t34 = pkin(3) * t138 + pkin(8) * t39;
t29 = -pkin(8) * t86 - t38;
t28 = -pkin(3) * t155 + pkin(8) * t88 + t39;
t27 = -t254 * t51 + t259 * t52;
t26 = t254 * t52 + t259 * t51;
t25 = -t254 * t48 + t259 * t49;
t24 = t254 * t49 + t259 * t48;
t23 = -pkin(4) * t82 + t33;
t22 = -pkin(4) * t79 + t32;
t21 = -t254 * t43 + t259 * t44;
t20 = t254 * t44 + t259 * t43;
t19 = t259 * t39 - t326;
t18 = t254 * t39 + t325;
t14 = t17 * t258 + t253 * t65;
t11 = -pkin(9) * t59 - t16;
t10 = -pkin(8) * t51 - t23 * t253 + t258 * t36;
t9 = -pkin(8) * t48 - t22 * t253 + t258 * t35;
t8 = -pkin(3) * t82 + pkin(8) * t52 + t23 * t258 + t253 * t36;
t7 = -pkin(3) * t79 + pkin(8) * t49 + t22 * t258 + t253 * t35;
t6 = -pkin(8) * t43 + t11 * t258 + t59 * t328;
t5 = pkin(8) * t44 + t11 * t253 + t290 * t59;
t4 = -t13 * t254 + t14 * t259;
t3 = t13 * t259 + t14 * t254;
t2 = -pkin(8) * t13 + (-pkin(9) * t258 + t328) * t16;
t1 = pkin(8) * t14 + (-pkin(9) * t253 + t290) * t16;
t12 = [0, 0, 0, 0, 0, qJDD(1), t286, t271, 0, 0, (t230 + t288) * t255, t229 * t260 + t232 * t255, t305 + t260 * (-t245 + t261), (t231 - t287) * t260, t255 * (t246 - t261) + t303, 0, t260 * t225 + pkin(1) * t232 + pkin(6) * (t260 * (-t246 - t261) - t305), -t255 * t225 - pkin(1) * t229 + pkin(6) * (-t303 - t255 * (-t245 - t261)), pkin(1) * (t245 + t246) + (t250 + t251) * t324 + t280, pkin(1) * t225 + pkin(6) * t280, t255 * (t193 * t259 - t224 * t307) + t260 * (t193 * t254 + t224 * t306), t255 * (-t174 * t259 - t178 * t254) + t260 * (-t174 * t254 + t178 * t259), t255 * (-t218 * t254 + t334) + t260 * (t218 * t259 + t337), t255 * (-t192 * t254 - t223 * t306) + t260 * (t192 * t259 - t223 * t307), t255 * (t217 * t259 - t312) + t260 * (t217 * t254 + t311), (t255 * (t223 * t259 + t224 * t254) + t260 * (t223 * t254 - t224 * t259)) * t249, t255 * (-pkin(7) * t170 - t316) + t260 * (-pkin(2) * t174 + pkin(7) * t171 + t315) - pkin(1) * t174 + pkin(6) * (-t170 * t255 + t171 * t260), t255 * (-pkin(7) * t180 - t315) + t260 * (-pkin(2) * t177 + pkin(7) * t181 - t316) - pkin(1) * t177 + pkin(6) * (-t180 * t255 + t181 * t260), t255 * (-pkin(7) * t143 - t116) + t260 * (-pkin(2) * t194 + pkin(7) * t144 + t117) - pkin(1) * t194 + pkin(6) * (-t143 * t255 + t144 * t260), -pkin(7) * t321 + t260 * (pkin(2) * t195 + pkin(7) * t117) + pkin(1) * t195 + pkin(6) * (t117 * t260 - t321), t255 * (-t124 * t254 + t125 * t259) + t260 * (t124 * t259 + t125 * t254), t255 * (-t254 * t85 + t259 * t87) + t260 * (t254 * t87 + t259 * t85), t255 * (-t149 * t254 + t151 * t259) + t260 * (t149 * t259 + t151 * t254), t255 * (-t122 * t254 + t123 * t259) + t260 * (t122 * t259 + t123 * t254), t255 * (-t150 * t254 + t152 * t259) + t260 * (t150 * t259 + t152 * t254), t255 * (-t160 * t254 + t161 * t259) + t260 * (t160 * t259 + t161 * t254), t255 * (-pkin(7) * t89 - t254 * t69 + t259 * t97) + t260 * (-pkin(2) * t126 + pkin(7) * t90 + t254 * t97 + t259 * t69) - pkin(1) * t126 + pkin(6) * (-t255 * t89 + t260 * t90), t255 * (-pkin(7) * t111 + t110 * t259 - t254 * t70) + t260 * (-pkin(2) * t131 + pkin(7) * t112 + t110 * t254 + t259 * t70) - pkin(1) * t131 + pkin(6) * (-t111 * t255 + t112 * t260), t255 * (-pkin(7) * t40 - t254 * t28 + t259 * t29) + t260 * (-pkin(2) * t155 + pkin(7) * t41 + t254 * t29 + t259 * t28) - pkin(1) * t155 + pkin(6) * (-t255 * t40 + t260 * t41), t255 * (-pkin(7) * t18 - pkin(8) * t325 - t254 * t34) + t260 * (pkin(2) * t138 + pkin(7) * t19 - pkin(8) * t326 + t259 * t34) + pkin(1) * t138 + pkin(6) * (-t18 * t255 + t19 * t260), t255 * (-t254 * t75 + t259 * t77) + t260 * (t254 * t77 + t259 * t75), t255 * (-t254 * t45 + t259 * t46) + t260 * (t254 * t46 + t259 * t45), t255 * (-t254 * t53 + t259 * t55) + t260 * (t254 * t55 + t259 * t53), t255 * (-t254 * t74 + t259 * t76) + t260 * (t254 * t76 + t259 * t74), t255 * (-t254 * t54 + t259 * t56) + t260 * (t254 * t56 + t259 * t54), t255 * (-t254 * t95 + t259 * t96) + t260 * (t254 * t96 + t259 * t95), t255 * (-pkin(7) * t24 - t254 * t7 + t259 * t9) + t260 * (-pkin(2) * t79 + pkin(7) * t25 + t254 * t9 + t259 * t7) - pkin(1) * t79 + pkin(6) * (-t24 * t255 + t25 * t260), t255 * (-pkin(7) * t26 + t10 * t259 - t254 * t8) + t260 * (-pkin(2) * t82 + pkin(7) * t27 + t10 * t254 + t259 * t8) - pkin(1) * t82 + pkin(6) * (-t255 * t26 + t260 * t27), t255 * (-pkin(7) * t20 - t254 * t5 + t259 * t6) + t260 * (-pkin(2) * t59 + pkin(7) * t21 + t254 * t6 + t259 * t5) - pkin(1) * t59 + pkin(6) * (-t20 * t255 + t21 * t260), t255 * (-pkin(7) * t3 - t1 * t254 + t2 * t259) + t260 * (-pkin(2) * t16 + pkin(7) * t4 + t1 * t259 + t2 * t254) - pkin(1) * t16 + pkin(6) * (-t255 * t3 + t260 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t238, t245 - t246, t243, t238, t244, qJDD(2), -t215, -t216, 0, 0, -t212, t211, t179, t212, t267, t296, pkin(2) * t170 - t156, pkin(2) * t180 - t157, pkin(2) * t143, pkin(2) * t116, t173, t172, t130, -t173, t268, t241, pkin(2) * t89 + t275, pkin(2) * t111 + t263, pkin(2) * t40 + t84, pkin(2) * t18 + t37, t100, t58, t91, t98, t92, t132, pkin(2) * t24 + t283, pkin(2) * t26 + t284, pkin(2) * t20 + t278, pkin(2) * t3 + t295; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t212, t211, t179, t212, t267, t296, -t156, -t157, 0, 0, t173, t172, t130, -t173, t268, t241, t275, t263, t84, t37, t100, t58, t91, t98, t92, t132, t283, t284, t278, t295; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t173, t172, t130, -t173, t268, t241, -t72, -t73, 0, 0, t100, t58, t91, t98, t92, t132, t293, t294, t285, t327; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t159, t158, t108, -t159, -t104, t140, -t32, -t33, 0, 0;];
tauJ_reg = t12;
