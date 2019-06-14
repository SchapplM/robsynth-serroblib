% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPPRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPPRPR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:32:12
% EndTime: 2019-05-05 14:32:27
% DurationCPUTime: 6.29s
% Computational Cost: add. (25505->422), mult. (60730->608), div. (0->0), fcn. (44156->10), ass. (0->265)
t237 = sin(qJ(4));
t233 = sin(pkin(9));
t240 = cos(qJ(4));
t235 = cos(pkin(9));
t274 = t235 * t237;
t254 = t233 * t240 + t274;
t213 = t254 * qJD(1);
t275 = t233 * t237;
t215 = (t235 * t240 - t275) * qJD(1);
t276 = t215 * t213;
t319 = qJDD(4) - t276;
t322 = t237 * t319;
t321 = t240 * t319;
t266 = t235 * qJDD(1);
t267 = t233 * qJDD(1);
t212 = -t237 * t267 + t240 * t266;
t271 = qJD(4) * t213;
t189 = t212 - t271;
t232 = sin(pkin(10));
t234 = cos(pkin(10));
t171 = qJDD(4) * t232 + t189 * t234;
t197 = -t234 * qJD(4) + t215 * t232;
t279 = t213 * t197;
t139 = -t171 - t279;
t320 = t171 - t279;
t199 = qJD(4) * t232 + t215 * t234;
t165 = t199 * t197;
t270 = t215 * qJD(4);
t307 = t254 * qJDD(1);
t187 = t307 + t270;
t308 = -t165 + t187;
t318 = t232 * t308;
t317 = t234 * t308;
t236 = sin(qJ(6));
t239 = cos(qJ(6));
t159 = t239 * t197 + t199 * t236;
t161 = -t197 * t236 + t199 * t239;
t126 = t161 * t159;
t179 = qJDD(6) + t187;
t310 = -t126 + t179;
t316 = t236 * t310;
t315 = t239 * t310;
t242 = qJD(1) ^ 2;
t238 = sin(qJ(1));
t241 = cos(qJ(1));
t255 = t238 * g(1) - t241 * g(2);
t252 = qJDD(2) - t255;
t249 = -t242 * qJ(2) + t252;
t265 = -0.2e1 * qJD(3) * qJD(1);
t300 = pkin(1) + qJ(3);
t314 = -qJDD(1) * t300 + t249 + t265;
t227 = t233 ^ 2;
t228 = t235 ^ 2;
t272 = t227 + t228;
t313 = pkin(3) * t267 - (pkin(7) * t272 + t300) * t242;
t312 = -pkin(7) - t300;
t311 = t272 * t242;
t258 = -t234 * qJDD(4) + t189 * t232;
t114 = -t159 * qJD(6) + t239 * t171 - t236 * t258;
t207 = qJD(6) + t213;
t151 = t207 * t159;
t309 = -t151 + t114;
t259 = t171 * t236 + t239 * t258;
t91 = (qJD(6) - t207) * t161 + t259;
t157 = t159 ^ 2;
t158 = t161 ^ 2;
t304 = t197 ^ 2;
t196 = t199 ^ 2;
t206 = t207 ^ 2;
t303 = t213 ^ 2;
t211 = t215 ^ 2;
t302 = qJD(4) ^ 2;
t301 = t233 * g(3);
t178 = pkin(4) * t213 - qJ(5) * t215;
t246 = (t265 + (-t233 * pkin(3) - qJ(2)) * t242 + t312 * qJDD(1) + t252) * t235;
t245 = t246 + t301;
t181 = -g(3) * t235 + t233 * t314;
t167 = -pkin(3) * t227 * t242 - pkin(7) * t267 + t181;
t273 = t240 * t167;
t104 = -pkin(4) * t302 + qJDD(4) * qJ(5) - t213 * t178 + t237 * t245 + t273;
t268 = qJD(2) * qJD(1);
t226 = 0.2e1 * t268;
t229 = qJDD(1) * qJ(2);
t256 = g(1) * t241 + g(2) * t238;
t253 = -t229 + t256;
t251 = -qJDD(3) + t253;
t120 = t226 + (-t189 + t271) * qJ(5) + (t187 + t270) * pkin(4) - t251 + t313;
t66 = 0.2e1 * qJD(5) * t199 + t104 * t232 - t234 * t120;
t48 = t308 * pkin(5) + pkin(8) * t139 - t66;
t170 = pkin(5) * t213 - pkin(8) * t199;
t67 = -0.2e1 * qJD(5) * t197 + t234 * t104 + t232 * t120;
t55 = -pkin(5) * t304 - pkin(8) * t258 - t213 * t170 + t67;
t26 = t236 * t55 - t239 * t48;
t27 = t236 * t48 + t239 * t55;
t14 = t236 * t27 - t239 * t26;
t299 = t14 * t232;
t298 = t14 * t234;
t128 = t167 * t237 - t240 * t245;
t129 = g(3) * t275 + t237 * t246 + t273;
t96 = -t128 * t240 + t129 * t237;
t297 = t235 * t96;
t103 = -qJDD(4) * pkin(4) - t302 * qJ(5) + t178 * t215 + qJDD(5) + t128;
t71 = pkin(5) * t258 - pkin(8) * t304 + t170 * t199 + t103;
t296 = t236 * t71;
t295 = t239 * t71;
t294 = qJDD(1) * pkin(1);
t293 = t103 * t232;
t292 = t103 * t234;
t118 = t126 + t179;
t291 = t118 * t236;
t290 = t118 * t239;
t141 = t165 + t187;
t289 = t141 * t232;
t288 = t141 * t234;
t248 = t251 - 0.2e1 * t268;
t176 = t248 - t313;
t287 = t176 * t237;
t286 = t176 * t240;
t184 = qJDD(4) + t276;
t285 = t184 * t237;
t284 = t184 * t240;
t283 = t187 * t237;
t282 = t199 * t213;
t281 = t207 * t236;
t280 = t207 * t239;
t278 = t213 * t232;
t277 = t213 * t234;
t264 = t237 * t126;
t263 = t240 * t126;
t262 = t237 * t165;
t261 = t240 * t165;
t260 = -pkin(4) * t240 - pkin(3);
t41 = t232 * t66 + t234 * t67;
t15 = t236 * t26 + t239 * t27;
t97 = t128 * t237 + t240 * t129;
t201 = t242 * t300 + t248;
t257 = -t201 + t229;
t40 = t232 * t67 - t234 * t66;
t19 = t233 * (t103 * t237 + t240 * t41) + t235 * (-t103 * t240 + t237 * t41);
t143 = (t235 * t314 + t301) * t235 + t181 * t233;
t135 = t258 - t282;
t218 = t272 * qJDD(1);
t217 = t233 * t311;
t216 = t235 * t311;
t208 = -t249 + t294;
t204 = -t211 - t302;
t203 = -t211 + t302;
t202 = t303 - t302;
t188 = t212 - 0.2e1 * t271;
t186 = t307 + 0.2e1 * t270;
t182 = -t303 - t302;
t177 = t240 * t187;
t173 = -t196 + t303;
t172 = -t303 + t304;
t166 = -t303 - t211;
t163 = -t196 + t304;
t156 = -t196 - t303;
t155 = -t204 * t237 - t284;
t154 = t204 * t240 - t285;
t153 = -t303 - t304;
t150 = t212 * t237 - t240 * t307;
t149 = -t212 * t240 - t237 * t307;
t148 = -t158 + t206;
t147 = t157 - t206;
t146 = t182 * t240 - t322;
t145 = t182 * t237 + t321;
t144 = -t196 - t304;
t134 = t258 + t282;
t133 = (-t197 * t234 + t199 * t232) * t213;
t132 = t171 * t234 - t199 * t278;
t131 = t197 * t277 + t232 * t258;
t130 = -t158 - t206;
t125 = t158 - t157;
t124 = t154 * t235 + t155 * t233;
t123 = -t206 - t157;
t122 = t172 * t234 - t289;
t121 = -t173 * t232 + t317;
t116 = -t156 * t232 - t288;
t115 = t156 * t234 - t289;
t113 = -qJD(6) * t161 - t259;
t110 = (-t159 * t239 + t161 * t236) * t207;
t109 = (-t159 * t236 - t161 * t239) * t207;
t108 = t149 * t235 + t150 * t233;
t107 = t153 * t234 - t318;
t106 = t153 * t232 + t317;
t105 = t145 * t235 + t146 * t233;
t101 = -t157 - t158;
t100 = -t135 * t234 - t139 * t232;
t99 = -t134 * t234 - t232 * t320;
t98 = -t135 * t232 + t139 * t234;
t94 = t151 + t114;
t90 = (qJD(6) + t207) * t161 + t259;
t89 = t147 * t239 - t291;
t88 = -t148 * t236 + t315;
t87 = t147 * t236 + t290;
t86 = t148 * t239 + t316;
t85 = t114 * t239 - t161 * t281;
t84 = t114 * t236 + t161 * t280;
t83 = -t113 * t236 + t159 * t280;
t82 = t113 * t239 + t159 * t281;
t81 = t116 * t240 + t237 * t320;
t80 = t116 * t237 - t240 * t320;
t79 = -t130 * t236 - t290;
t78 = t130 * t239 - t291;
t77 = t107 * t240 + t134 * t237;
t76 = t107 * t237 - t134 * t240;
t75 = t100 * t240 + t144 * t237;
t74 = t100 * t237 - t144 * t240;
t73 = t123 * t239 - t316;
t72 = t123 * t236 + t315;
t70 = -t109 * t232 + t110 * t234;
t69 = -qJ(5) * t115 + t292;
t68 = -qJ(5) * t106 + t293;
t64 = t233 * t97 + t297;
t63 = t236 * t94 - t239 * t91;
t62 = -t236 * t309 - t239 * t90;
t61 = -t236 * t91 - t239 * t94;
t60 = -t236 * t90 + t239 * t309;
t59 = -t232 * t87 + t234 * t89;
t58 = -t232 * t86 + t234 * t88;
t57 = -t232 * t84 + t234 * t85;
t56 = -t232 * t82 + t234 * t83;
t54 = t233 * t81 + t235 * t80;
t53 = -pkin(4) * t115 + t67;
t52 = -pkin(4) * t106 + t66;
t51 = -t232 * t78 + t234 * t79;
t50 = t232 * t79 + t234 * t78;
t49 = t233 * t77 + t235 * t76;
t46 = t233 * t75 + t235 * t74;
t45 = -pkin(8) * t78 + t295;
t44 = -t232 * t72 + t234 * t73;
t43 = t232 * t73 + t234 * t72;
t42 = -pkin(8) * t72 + t296;
t39 = t237 * t309 + t240 * t51;
t38 = t237 * t51 - t240 * t309;
t37 = -pkin(5) * t309 + pkin(8) * t79 + t296;
t36 = t237 * t90 + t240 * t44;
t35 = t237 * t44 - t240 * t90;
t34 = -pkin(5) * t90 + pkin(8) * t73 - t295;
t31 = -qJ(5) * t98 - t40;
t30 = -t232 * t61 + t234 * t63;
t29 = -t232 * t60 + t234 * t62;
t28 = t232 * t63 + t234 * t61;
t24 = t101 * t237 + t240 * t30;
t23 = -t101 * t240 + t237 * t30;
t22 = -pkin(4) * t28 - pkin(5) * t61;
t21 = t233 * t39 + t235 * t38;
t20 = t233 * t36 + t235 * t35;
t18 = -pkin(4) * t50 - pkin(5) * t78 + t27;
t17 = -qJ(5) * t50 - t232 * t37 + t234 * t45;
t16 = -pkin(4) * t43 - pkin(5) * t72 + t26;
t13 = -qJ(5) * t43 - t232 * t34 + t234 * t42;
t12 = t23 * t235 + t233 * t24;
t11 = -pkin(5) * t71 + pkin(8) * t15;
t10 = -pkin(8) * t61 - t14;
t9 = -pkin(5) * t101 + pkin(8) * t63 + t15;
t8 = t15 * t234 - t299;
t7 = t15 * t232 + t298;
t6 = t237 * t71 + t240 * t8;
t5 = t237 * t8 - t240 * t71;
t4 = -pkin(4) * t7 - pkin(5) * t14;
t3 = -qJ(5) * t28 + t10 * t234 - t232 * t9;
t2 = -pkin(8) * t298 - qJ(5) * t7 - t11 * t232;
t1 = t233 * t6 + t235 * t5;
t25 = [0, 0, 0, 0, 0, qJDD(1), t255, t256, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t252 - 0.2e1 * t294, t226 + 0.2e1 * t229 - t256, pkin(1) * t208 + qJ(2) * (-pkin(1) * t242 + t226 - t253), t228 * qJDD(1), -0.2e1 * t233 * t266, 0, t227 * qJDD(1), 0, 0, t217 * t300 + t233 * t257, t216 * t300 + t235 * t257, -qJ(2) * t311 + t218 * t300 - t143, -qJ(2) * t201 - t143 * t300, t235 * (t189 * t240 - t237 * t270) - t233 * (t189 * t237 + t240 * t270), t235 * (-t186 * t240 - t188 * t237) - t233 * (-t186 * t237 + t188 * t240), t235 * (-t203 * t237 + t321) - t233 * (t203 * t240 + t322), t235 * (t240 * t271 + t283) - t233 * (t237 * t271 - t177), t235 * (t202 * t240 - t285) - t233 * (t202 * t237 + t284), (t235 * (-t213 * t240 + t215 * t237) - t233 * (-t213 * t237 - t215 * t240)) * qJD(4), t235 * (-pkin(7) * t145 - t287) - t233 * (-pkin(3) * t186 + pkin(7) * t146 + t286) + qJ(2) * t186 - t300 * t105, t235 * (-pkin(7) * t154 - t286) - t233 * (-pkin(3) * t188 + pkin(7) * t155 - t287) + qJ(2) * t188 - t300 * t124, t235 * (-pkin(7) * t149 - t96) - t233 * (-pkin(3) * t166 + pkin(7) * t150 + t97) + qJ(2) * t166 - t300 * t108, -pkin(7) * t297 - t233 * (pkin(3) * t176 + pkin(7) * t97) - qJ(2) * t176 - t300 * t64, t235 * (t132 * t240 + t262) - t233 * (t132 * t237 - t261), t235 * (-t163 * t237 + t240 * t99) - t233 * (t163 * t240 + t237 * t99), t235 * (t121 * t240 - t139 * t237) - t233 * (t121 * t237 + t139 * t240), t235 * (t131 * t240 - t262) - t233 * (t131 * t237 + t261), t235 * (t122 * t240 - t135 * t237) - t233 * (t122 * t237 + t135 * t240), t235 * (t133 * t240 + t283) - t233 * (t133 * t237 - t177), t235 * (-pkin(7) * t76 - t237 * t52 + t240 * t68) - t233 * (-pkin(3) * t106 + pkin(7) * t77 + t237 * t68 + t240 * t52) + qJ(2) * t106 - t300 * t49, t235 * (-pkin(7) * t80 - t237 * t53 + t240 * t69) - t233 * (-pkin(3) * t115 + pkin(7) * t81 + t237 * t69 + t240 * t53) + qJ(2) * t115 - t300 * t54, t235 * (-pkin(7) * t74 + t240 * t31) - t233 * (pkin(7) * t75 + t237 * t31) + (pkin(4) * t274 - t233 * t260 + qJ(2)) * t98 - t300 * t46, (t235 * (pkin(4) * t237 - qJ(5) * t240) - t233 * (-qJ(5) * t237 + t260) + qJ(2)) * t40 + t312 * t19, t235 * (t240 * t57 + t264) - t233 * (t237 * t57 - t263), t235 * (t125 * t237 + t240 * t29) - t233 * (-t125 * t240 + t237 * t29), t235 * (t237 * t94 + t240 * t58) - t233 * (t237 * t58 - t240 * t94), t235 * (t240 * t56 - t264) - t233 * (t237 * t56 + t263), t235 * (-t237 * t91 + t240 * t59) - t233 * (t237 * t59 + t240 * t91), t235 * (t179 * t237 + t240 * t70) - t233 * (-t179 * t240 + t237 * t70), t235 * (-pkin(7) * t35 + t13 * t240 - t16 * t237) - t233 * (-pkin(3) * t43 + pkin(7) * t36 + t13 * t237 + t16 * t240) + qJ(2) * t43 - t300 * t20, t235 * (-pkin(7) * t38 + t17 * t240 - t18 * t237) - t233 * (-pkin(3) * t50 + pkin(7) * t39 + t17 * t237 + t18 * t240) + qJ(2) * t50 - t300 * t21, t235 * (-pkin(7) * t23 - t22 * t237 + t240 * t3) - t233 * (-pkin(3) * t28 + pkin(7) * t24 + t22 * t240 + t237 * t3) + qJ(2) * t28 - t300 * t12, t235 * (-pkin(7) * t5 + t2 * t240 - t237 * t4) - t233 * (-pkin(3) * t7 + pkin(7) * t6 + t2 * t237 + t240 * t4) + qJ(2) * t7 - t300 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t242, -t208, 0, 0, 0, 0, 0, 0, -t217, -t216, -t218, t143, 0, 0, 0, 0, 0, 0, t105, t124, t108, t64, 0, 0, 0, 0, 0, 0, t49, t54, t46, t19, 0, 0, 0, 0, 0, 0, t20, t21, t12, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t267, t266, -t311, -t201, 0, 0, 0, 0, 0, 0, t186, t188, t166, -t176, 0, 0, 0, 0, 0, 0, t106, t115, t98, t40, 0, 0, 0, 0, 0, 0, t43, t50, t28, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t276, t211 - t303, t212, -t276, -t307, qJDD(4), -t128, -t129, 0, 0, t171 * t232 + t199 * t277, -t134 * t232 + t234 * t320, t173 * t234 + t318, t197 * t278 - t234 * t258, t172 * t232 + t288, (-t197 * t232 - t199 * t234) * t213, -pkin(4) * t134 + qJ(5) * t107 - t292, -pkin(4) * t320 + qJ(5) * t116 + t293, -pkin(4) * t144 + qJ(5) * t100 + t41, -pkin(4) * t103 + qJ(5) * t41, t232 * t85 + t234 * t84, t232 * t62 + t234 * t60, t232 * t88 + t234 * t86, t232 * t83 + t234 * t82, t232 * t89 + t234 * t87, t109 * t234 + t110 * t232, -pkin(4) * t90 + qJ(5) * t44 + t232 * t42 + t234 * t34, -pkin(4) * t309 + qJ(5) * t51 + t232 * t45 + t234 * t37, -pkin(4) * t101 + qJ(5) * t30 + t10 * t232 + t234 * t9, -pkin(4) * t71 - pkin(8) * t299 + qJ(5) * t8 + t11 * t234; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, t320, t144, t103, 0, 0, 0, 0, 0, 0, t90, t309, t101, t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t126, t125, t94, -t126, -t91, t179, -t26, -t27, 0, 0;];
tauJ_reg  = t25;
