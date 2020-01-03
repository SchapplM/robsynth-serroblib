% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPRR8
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
% Datum: 2019-12-31 20:19
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPRR8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR8_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:18:16
% EndTime: 2019-12-31 20:18:28
% DurationCPUTime: 5.71s
% Computational Cost: add. (25749->448), mult. (60332->643), div. (0->0), fcn. (44170->10), ass. (0->275)
t249 = cos(qJ(2));
t233 = t249 * qJDD(1);
t245 = sin(qJ(2));
t276 = qJD(1) * qJD(2);
t269 = t245 * t276;
t222 = t233 - t269;
t240 = t249 ^ 2;
t251 = qJD(1) ^ 2;
t284 = qJD(1) * t245;
t256 = qJD(2) * pkin(2) - qJ(3) * t284;
t246 = sin(qJ(1));
t311 = cos(qJ(1));
t267 = g(1) * t246 - t311 * g(2);
t257 = qJDD(1) * pkin(1) + t267;
t174 = pkin(2) * t222 - t256 * t284 - qJDD(3) + t257 + (qJ(3) * t240 + pkin(6)) * t251;
t323 = 2 * qJD(3);
t241 = sin(pkin(9));
t242 = cos(pkin(9));
t214 = qJD(1) * t242 * t249 - t241 * t284;
t215 = (t249 * t241 + t245 * t242) * qJD(1);
t200 = t214 * t215;
t314 = qJDD(2) + t200;
t322 = t241 * t314;
t321 = t242 * t314;
t243 = sin(qJ(5));
t244 = sin(qJ(4));
t248 = cos(qJ(4));
t193 = t214 * t244 + t215 * t248;
t232 = t245 * qJDD(1);
t270 = t249 * t276;
t221 = t232 + t270;
t201 = -t221 * t241 + t242 * t222;
t202 = t221 * t242 + t222 * t241;
t264 = -t248 * t201 + t202 * t244;
t146 = -qJD(4) * t193 - t264;
t145 = qJDD(5) - t146;
t238 = qJD(2) + qJD(4);
t247 = cos(qJ(5));
t175 = t193 * t243 - t247 * t238;
t177 = t193 * t247 + t238 * t243;
t151 = t177 * t175;
t316 = t145 - t151;
t320 = t243 * t316;
t191 = -t248 * t214 + t215 * t244;
t163 = t193 * t191;
t237 = qJDD(2) + qJDD(4);
t315 = -t163 + t237;
t319 = t244 * t315;
t318 = t247 * t316;
t317 = t248 * t315;
t212 = t214 ^ 2;
t262 = qJD(2) * pkin(3) - pkin(7) * t215;
t136 = pkin(3) * t201 + pkin(7) * t212 - t215 * t262 + t174;
t161 = pkin(4) * t191 - pkin(8) * t193;
t312 = t238 ^ 2;
t288 = t245 * t251;
t258 = t311 * g(1) + t246 * g(2);
t305 = qJDD(1) * pkin(6);
t218 = -t251 * pkin(1) - t258 + t305;
t292 = t218 * t245;
t170 = qJDD(2) * pkin(2) - qJ(3) * t221 - t292 + (pkin(2) * t288 + qJ(3) * t276 - g(3)) * t249;
t205 = -t245 * g(3) + t249 * t218;
t235 = t240 * t251;
t171 = -pkin(2) * t235 + t222 * qJ(3) - qJD(2) * t256 + t205;
t137 = -t170 * t242 + t171 * t241 + t215 * t323;
t211 = qJD(2) * t214;
t182 = t202 - t211;
t313 = -t182 * pkin(7) - t137;
t252 = pkin(3) * t314 + t313;
t138 = t241 * t170 + t242 * t171 + t214 * t323;
t113 = -t212 * pkin(3) + t201 * pkin(7) - qJD(2) * t262 + t138;
t287 = t248 * t113;
t67 = t244 * t252 + t287;
t56 = -t312 * pkin(4) + t237 * pkin(8) - t191 * t161 + t67;
t259 = t201 * t244 + t202 * t248;
t147 = -qJD(4) * t191 + t259;
t186 = t238 * t191;
t128 = t147 - t186;
t64 = -t128 * pkin(8) + (t193 * t238 - t146) * pkin(4) - t136;
t31 = t243 * t56 - t247 * t64;
t32 = t243 * t64 + t247 * t56;
t16 = t243 * t31 + t247 * t32;
t188 = qJD(5) + t191;
t265 = t147 * t243 - t247 * t237;
t103 = (qJD(5) - t188) * t177 + t265;
t172 = t175 ^ 2;
t173 = t177 ^ 2;
t187 = t188 ^ 2;
t189 = t191 ^ 2;
t190 = t193 ^ 2;
t213 = t215 ^ 2;
t310 = pkin(4) * t244;
t66 = t113 * t244 - t248 * t252;
t55 = -t237 * pkin(4) - t312 * pkin(8) + t161 * t193 + t66;
t309 = -pkin(4) * t55 + pkin(8) * t16;
t36 = t244 * t67 - t248 * t66;
t308 = t241 * t36;
t307 = t242 * t36;
t52 = t243 * t55;
t95 = -t137 * t242 + t138 * t241;
t306 = t245 * t95;
t53 = t247 * t55;
t111 = t145 + t151;
t304 = t111 * t243;
t303 = t111 * t247;
t302 = t136 * t244;
t301 = t136 * t248;
t159 = t163 + t237;
t300 = t159 * t244;
t299 = t159 * t248;
t298 = t174 * t241;
t297 = t174 * t242;
t296 = t188 * t243;
t295 = t188 * t247;
t197 = qJDD(2) - t200;
t294 = t197 * t241;
t293 = t197 * t242;
t291 = t238 * t244;
t290 = t238 * t248;
t227 = t249 * t288;
t289 = t245 * (qJDD(2) + t227);
t286 = t249 * (qJDD(2) - t227);
t283 = qJD(2) * t215;
t282 = qJD(2) * t241;
t281 = qJD(2) * t242;
t279 = qJD(4) + t238;
t277 = qJD(5) + t188;
t260 = -t147 * t247 - t237 * t243;
t108 = t277 * t175 + t260;
t143 = -t173 - t187;
t80 = -t143 * t243 - t303;
t275 = pkin(4) * t108 + pkin(8) * t80 + t52;
t104 = -t277 * t177 - t265;
t132 = -t187 - t172;
t75 = t132 * t247 - t320;
t274 = pkin(4) * t104 + pkin(8) * t75 - t53;
t273 = t244 * t151;
t272 = t248 * t151;
t271 = -pkin(4) * t248 - pkin(3);
t37 = t244 * t66 + t248 * t67;
t131 = t172 + t173;
t117 = -qJD(5) * t175 - t260;
t156 = t188 * t175;
t107 = t117 + t156;
t61 = -t103 * t247 + t107 * t243;
t266 = pkin(4) * t131 + pkin(8) * t61 + t16;
t96 = t137 * t241 + t242 * t138;
t204 = g(3) * t249 + t292;
t263 = t245 * t204 + t249 * t205;
t15 = t243 * t32 - t247 * t31;
t180 = t201 + t283;
t255 = (-qJD(4) + t238) * t193 - t264;
t250 = qJD(2) ^ 2;
t239 = t245 ^ 2;
t234 = t239 * t251;
t223 = t233 - 0.2e1 * t269;
t220 = t232 + 0.2e1 * t270;
t217 = pkin(6) * t251 + t257;
t208 = -t213 - t250;
t207 = -t213 + t250;
t206 = t212 - t250;
t195 = -t250 - t212;
t185 = -t190 + t312;
t184 = t189 - t312;
t183 = -t190 - t312;
t181 = t202 + t211;
t179 = -t201 + t283;
t178 = -t212 - t213;
t167 = -t208 * t241 - t293;
t166 = t208 * t242 - t294;
t165 = t195 * t242 - t322;
t164 = t195 * t241 + t321;
t162 = t190 - t189;
t157 = -t312 - t189;
t155 = -t173 + t187;
t154 = t172 - t187;
t153 = (-t191 * t248 + t193 * t244) * t238;
t152 = (-t191 * t244 - t193 * t248) * t238;
t150 = t173 - t172;
t149 = t180 * t242 + t182 * t241;
t148 = t180 * t241 - t182 * t242;
t144 = -t189 - t190;
t142 = t184 * t248 - t300;
t141 = -t185 * t244 + t317;
t140 = t184 * t244 + t299;
t139 = t185 * t248 + t319;
t135 = -t183 * t244 - t299;
t134 = t183 * t248 - t300;
t129 = t147 + t186;
t127 = -t279 * t191 + t259;
t124 = t279 * t193 + t264;
t123 = t147 * t248 - t193 * t291;
t122 = t147 * t244 + t193 * t290;
t121 = -t146 * t244 + t191 * t290;
t120 = t146 * t248 + t191 * t291;
t119 = t157 * t248 - t319;
t118 = t157 * t244 + t317;
t116 = -qJD(5) * t177 - t265;
t115 = (-t175 * t247 + t177 * t243) * t188;
t114 = (-t175 * t243 - t177 * t247) * t188;
t106 = t117 - t156;
t100 = t117 * t247 - t177 * t296;
t99 = t117 * t243 + t177 * t295;
t98 = -t116 * t243 + t175 * t295;
t97 = t116 * t247 + t175 * t296;
t94 = -pkin(7) * t134 - t301;
t93 = -t134 * t241 + t135 * t242;
t92 = t134 * t242 + t135 * t241;
t91 = t115 * t248 + t145 * t244;
t90 = t115 * t244 - t145 * t248;
t89 = t154 * t247 - t304;
t88 = -t155 * t243 + t318;
t87 = t154 * t243 + t303;
t86 = t155 * t247 + t320;
t85 = t129 * t244 + t248 * t255;
t84 = -t124 * t248 - t128 * t244;
t83 = -t129 * t248 + t244 * t255;
t82 = -t124 * t244 + t128 * t248;
t81 = -pkin(7) * t118 - t302;
t79 = t143 * t247 - t304;
t77 = -t118 * t241 + t119 * t242;
t76 = t118 * t242 + t119 * t241;
t74 = t132 * t243 + t318;
t72 = t100 * t248 + t273;
t71 = t248 * t98 - t273;
t70 = t100 * t244 - t272;
t69 = t244 * t98 + t272;
t68 = -pkin(3) * t127 + pkin(7) * t135 - t302;
t62 = -pkin(3) * t124 + pkin(7) * t119 + t301;
t60 = t104 * t247 - t106 * t243;
t59 = -t103 * t243 - t107 * t247;
t58 = t104 * t243 + t106 * t247;
t51 = -t103 * t244 + t248 * t89;
t50 = t107 * t244 + t248 * t88;
t49 = t103 * t248 + t244 * t89;
t48 = -t107 * t248 + t244 * t88;
t47 = -t108 * t244 + t248 * t80;
t46 = t108 * t248 + t244 * t80;
t45 = -t104 * t244 + t248 * t75;
t44 = t104 * t248 + t244 * t75;
t43 = -t241 * t83 + t242 * t85;
t42 = t241 * t85 + t242 * t83;
t41 = t150 * t244 + t248 * t60;
t40 = -t150 * t248 + t244 * t60;
t39 = -t131 * t244 + t248 * t61;
t38 = t131 * t248 + t244 * t61;
t35 = -pkin(8) * t79 + t53;
t34 = -pkin(8) * t74 + t52;
t33 = pkin(3) * t136 + pkin(7) * t37;
t28 = -pkin(7) * t83 - t36;
t27 = -t241 * t46 + t242 * t47;
t26 = t241 * t47 + t242 * t46;
t25 = -t241 * t44 + t242 * t45;
t24 = t241 * t45 + t242 * t44;
t23 = -pkin(3) * t144 + pkin(7) * t85 + t37;
t22 = -pkin(4) * t79 + t32;
t21 = -pkin(4) * t74 + t31;
t20 = -t241 * t38 + t242 * t39;
t19 = t241 * t39 + t242 * t38;
t18 = t242 * t37 - t308;
t17 = t241 * t37 + t307;
t13 = -pkin(8) * t59 - t15;
t12 = t16 * t248 + t244 * t55;
t11 = t16 * t244 - t248 * t55;
t10 = -pkin(7) * t46 - t22 * t244 + t248 * t35;
t9 = -pkin(7) * t44 - t21 * t244 + t248 * t34;
t8 = -pkin(3) * t79 + pkin(7) * t47 + t22 * t248 + t244 * t35;
t7 = -pkin(3) * t74 + pkin(7) * t45 + t21 * t248 + t244 * t34;
t6 = -pkin(7) * t38 + t13 * t248 + t59 * t310;
t5 = pkin(7) * t39 + t13 * t244 + t271 * t59;
t4 = -t11 * t241 + t12 * t242;
t3 = t11 * t242 + t12 * t241;
t2 = -pkin(7) * t11 + (-pkin(8) * t248 + t310) * t15;
t1 = pkin(7) * t12 + (-pkin(8) * t244 + t271) * t15;
t14 = [0, 0, 0, 0, 0, qJDD(1), t267, t258, 0, 0, (t221 + t270) * t245, t220 * t249 + t223 * t245, t289 + t249 * (-t234 + t250), (t222 - t269) * t249, t245 * (t235 - t250) + t286, 0, t249 * t217 + pkin(1) * t223 + pkin(6) * (t249 * (-t235 - t250) - t289), -t245 * t217 - pkin(1) * t220 + pkin(6) * (-t286 - t245 * (-t234 - t250)), pkin(1) * (t234 + t235) + (t239 + t240) * t305 + t263, pkin(1) * t217 + pkin(6) * t263, t245 * (t202 * t242 - t215 * t282) + t249 * (t202 * t241 + t215 * t281), t245 * (-t179 * t242 - t181 * t241) + t249 * (-t179 * t241 + t181 * t242), t245 * (-t207 * t241 + t321) + t249 * (t207 * t242 + t322), t245 * (-t201 * t241 - t214 * t281) + t249 * (t201 * t242 - t214 * t282), t245 * (t206 * t242 - t294) + t249 * (t206 * t241 + t293), (t245 * (t214 * t242 + t215 * t241) + t249 * (t214 * t241 - t215 * t242)) * qJD(2), t245 * (-qJ(3) * t164 - t298) + t249 * (-pkin(2) * t179 + qJ(3) * t165 + t297) - pkin(1) * t179 + pkin(6) * (-t164 * t245 + t165 * t249), t245 * (-qJ(3) * t166 - t297) + t249 * (-pkin(2) * t181 + qJ(3) * t167 - t298) - pkin(1) * t181 + pkin(6) * (-t166 * t245 + t167 * t249), t245 * (-qJ(3) * t148 - t95) + t249 * (-pkin(2) * t178 + qJ(3) * t149 + t96) - pkin(1) * t178 + pkin(6) * (-t148 * t245 + t149 * t249), -qJ(3) * t306 + t249 * (pkin(2) * t174 + qJ(3) * t96) + pkin(1) * t174 + pkin(6) * (t249 * t96 - t306), t245 * (-t122 * t241 + t123 * t242) + t249 * (t122 * t242 + t123 * t241), t245 * (-t241 * t82 + t242 * t84) + t249 * (t241 * t84 + t242 * t82), t245 * (-t139 * t241 + t141 * t242) + t249 * (t139 * t242 + t141 * t241), t245 * (-t120 * t241 + t121 * t242) + t249 * (t120 * t242 + t121 * t241), t245 * (-t140 * t241 + t142 * t242) + t249 * (t140 * t242 + t142 * t241), t245 * (-t152 * t241 + t153 * t242) + t249 * (t152 * t242 + t153 * t241), t245 * (-qJ(3) * t76 - t241 * t62 + t242 * t81) + t249 * (-pkin(2) * t124 + qJ(3) * t77 + t241 * t81 + t242 * t62) - pkin(1) * t124 + pkin(6) * (-t245 * t76 + t249 * t77), t245 * (-qJ(3) * t92 - t241 * t68 + t242 * t94) + t249 * (-pkin(2) * t127 + qJ(3) * t93 + t241 * t94 + t242 * t68) - pkin(1) * t127 + pkin(6) * (-t245 * t92 + t249 * t93), t245 * (-qJ(3) * t42 - t23 * t241 + t242 * t28) + t249 * (-pkin(2) * t144 + qJ(3) * t43 + t23 * t242 + t241 * t28) - pkin(1) * t144 + pkin(6) * (-t245 * t42 + t249 * t43), t245 * (-pkin(7) * t307 - qJ(3) * t17 - t241 * t33) + t249 * (pkin(2) * t136 - pkin(7) * t308 + qJ(3) * t18 + t242 * t33) + pkin(1) * t136 + pkin(6) * (-t17 * t245 + t18 * t249), t245 * (-t241 * t70 + t242 * t72) + t249 * (t241 * t72 + t242 * t70), t245 * (-t241 * t40 + t242 * t41) + t249 * (t241 * t41 + t242 * t40), t245 * (-t241 * t48 + t242 * t50) + t249 * (t241 * t50 + t242 * t48), t245 * (-t241 * t69 + t242 * t71) + t249 * (t241 * t71 + t242 * t69), t245 * (-t241 * t49 + t242 * t51) + t249 * (t241 * t51 + t242 * t49), t245 * (-t241 * t90 + t242 * t91) + t249 * (t241 * t91 + t242 * t90), t245 * (-qJ(3) * t24 - t241 * t7 + t242 * t9) + t249 * (-pkin(2) * t74 + qJ(3) * t25 + t241 * t9 + t242 * t7) - pkin(1) * t74 + pkin(6) * (-t24 * t245 + t249 * t25), t245 * (-qJ(3) * t26 + t10 * t242 - t241 * t8) + t249 * (-pkin(2) * t79 + qJ(3) * t27 + t10 * t241 + t242 * t8) - pkin(1) * t79 + pkin(6) * (-t245 * t26 + t249 * t27), t245 * (-qJ(3) * t19 - t241 * t5 + t242 * t6) + t249 * (-pkin(2) * t59 + qJ(3) * t20 + t241 * t6 + t242 * t5) - pkin(1) * t59 + pkin(6) * (-t19 * t245 + t20 * t249), t245 * (-qJ(3) * t3 - t1 * t241 + t2 * t242) + t249 * (-pkin(2) * t15 + qJ(3) * t4 + t1 * t242 + t2 * t241) - pkin(1) * t15 + pkin(6) * (-t245 * t3 + t249 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t227, t234 - t235, t232, t227, t233, qJDD(2), -t204, -t205, 0, 0, -t200, t213 - t212, t182, t200, t180, qJDD(2), pkin(2) * t164 - t137, pkin(2) * t166 - t138, pkin(2) * t148, pkin(2) * t95, t163, t162, t129, -t163, t255, t237, pkin(2) * t76 + pkin(3) * t118 - t66, pkin(2) * t92 - t287 - t244 * t313 + (-t244 * t314 + t134) * pkin(3), pkin(2) * t42 + pkin(3) * t83, pkin(2) * t17 + pkin(3) * t36, t99, t58, t86, t97, t87, t114, pkin(2) * t24 + pkin(3) * t44 + t274, pkin(2) * t26 + pkin(3) * t46 + t275, pkin(2) * t19 + pkin(3) * t38 + t266, pkin(2) * t3 + pkin(3) * t11 + t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, t181, t178, -t174, 0, 0, 0, 0, 0, 0, t124, t127, t144, -t136, 0, 0, 0, 0, 0, 0, t74, t79, t59, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t163, t162, t129, -t163, t255, t237, -t66, -t67, 0, 0, t99, t58, t86, t97, t87, t114, t274, t275, t266, t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t151, t150, t107, -t151, -t103, t145, -t31, -t32, 0, 0;];
tauJ_reg = t14;
