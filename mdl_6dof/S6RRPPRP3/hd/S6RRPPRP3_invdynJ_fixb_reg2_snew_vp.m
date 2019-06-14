% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RRPPRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-06 09:17
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RRPPRP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPRP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPRP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 09:16:36
% EndTime: 2019-05-06 09:16:43
% DurationCPUTime: 2.45s
% Computational Cost: add. (6554->311), mult. (13879->352), div. (0->0), fcn. (7330->6), ass. (0->199)
t299 = pkin(2) + pkin(8);
t272 = pkin(4) + qJ(3);
t181 = cos(qJ(2));
t234 = qJD(1) * qJD(2);
t224 = t181 * t234;
t178 = sin(qJ(2));
t231 = t178 * qJDD(1);
t135 = t224 + t231;
t123 = qJDD(5) + t135;
t177 = sin(qJ(5));
t180 = cos(qJ(5));
t240 = qJD(1) * t181;
t129 = -t180 * qJD(2) + t177 * t240;
t130 = qJD(2) * t177 + t180 * t240;
t94 = t129 * t130;
t289 = t123 - t94;
t298 = pkin(5) * t289;
t235 = pkin(3) + t299;
t276 = pkin(2) + pkin(3);
t297 = pkin(8) + t276;
t296 = t177 * t289;
t295 = t180 * t289;
t184 = qJD(1) ^ 2;
t287 = t181 * t184;
t228 = t178 * t287;
t144 = qJDD(2) + t228;
t134 = 0.2e1 * t224 + t231;
t183 = qJD(2) ^ 2;
t171 = t178 ^ 2;
t250 = t171 * t184;
t148 = t183 + t250;
t145 = -qJDD(2) + t228;
t245 = t181 * t145;
t294 = pkin(1) * t134 - pkin(7) * (t148 * t178 + t245);
t172 = t181 ^ 2;
t249 = t172 * t184;
t149 = t183 + t249;
t246 = t178 * t144;
t293 = pkin(7) * (t149 * t181 + t246);
t292 = pkin(7) - qJ(4);
t241 = qJD(1) * t178;
t153 = qJD(5) + t241;
t108 = t153 * t129;
t159 = t178 * t234;
t230 = t181 * qJDD(1);
t136 = -t159 + t230;
t84 = qJD(5) * t129 - t177 * qJDD(2) - t136 * t180;
t291 = t108 + t84;
t238 = qJD(6) * t129;
t122 = t130 ^ 2;
t151 = t153 ^ 2;
t90 = -t122 - t151;
t290 = pkin(5) * t90 - 0.2e1 * t238;
t288 = qJDD(1) * pkin(7);
t286 = pkin(2) * t148 - qJ(3) * t145;
t233 = qJD(1) * qJD(4);
t225 = t181 * t233;
t232 = qJD(3) * qJD(2);
t285 = 0.2e1 * t225 - 0.2e1 * t232;
t164 = 0.2e1 * t232;
t284 = t164 - 0.2e1 * t225;
t179 = sin(qJ(1));
t182 = cos(qJ(1));
t243 = t179 * g(1) - t182 * g(2);
t115 = qJDD(1) * pkin(1) + t184 * pkin(7) + t243;
t196 = -pkin(2) * t159 + t115;
t191 = t136 * pkin(3) - qJ(4) * t249 + qJDD(4) + t196;
t143 = -qJD(2) * pkin(3) - qJ(4) * t241;
t221 = (0.2e1 * qJD(3) + t143) * t178;
t24 = t299 * t136 + t272 * t135 + (t221 + (-pkin(8) * t178 + t272 * t181) * qJD(2)) * qJD(1) + t191;
t215 = g(1) * t182 + g(2) * t179;
t116 = -pkin(1) * t184 - t215 + t288;
t273 = t181 * g(3);
t102 = t178 * t116 + t273;
t223 = qJDD(2) * pkin(2) + t183 * qJ(3) - qJDD(3);
t204 = pkin(3) * t144 + t135 * qJ(4) + t223;
t190 = t102 - t204;
t216 = pkin(4) * t178 + pkin(8) * t181;
t239 = qJD(2) * t181;
t227 = qJ(4) * t239;
t248 = t178 * qJ(3);
t213 = -pkin(2) * t181 - t248;
t132 = t213 * qJD(1);
t236 = -0.2e1 * qJD(4) + t132;
t36 = -t183 * pkin(4) - qJDD(2) * pkin(8) + (t227 + (-t216 * qJD(1) + t236) * t178) * qJD(1) + t190;
t14 = t177 * t36 - t180 * t24;
t283 = -qJ(6) * t108 - 0.2e1 * qJD(6) * t130 + t14 - t298;
t218 = -t181 * t276 - pkin(1);
t282 = t218 - t248;
t281 = t245 - (-t183 + t249) * t178;
t101 = pkin(5) * t153 + qJ(6) * t130;
t121 = t129 ^ 2;
t208 = t180 * qJDD(2) - t177 * t136;
t83 = qJD(5) * t130 - t208;
t279 = -t121 * pkin(5) + t83 * qJ(6) - t153 * t101;
t85 = -t151 - t121;
t47 = t177 * t85 + t295;
t48 = t180 * t85 - t296;
t278 = -qJ(4) * t48 + t272 * t47;
t277 = qJDD(2) * pkin(4) - t183 * pkin(8) - t216 * t287;
t220 = qJD(1) * t132 + t116;
t72 = t220 * t178 - t223 + t273;
t242 = t171 + t172;
t133 = t242 * t288;
t274 = t136 * pkin(2);
t194 = (qJD(5) - t153) * t130 - t208;
t70 = -t108 + t84;
t31 = t177 * t194 - t180 * t70;
t33 = t177 * t70 + t180 * t194;
t77 = -t121 - t122;
t271 = pkin(7) * (t178 * t33 + t181 * t77) + pkin(1) * t31;
t65 = (-qJD(5) - t153) * t130 + t208;
t270 = pkin(7) * (t178 * t48 + t181 * t65) + pkin(1) * t47;
t79 = t123 + t94;
t265 = t177 * t79;
t53 = t180 * t90 - t265;
t263 = t180 * t79;
t54 = -t177 * t90 - t263;
t269 = pkin(7) * (t178 * t54 + t181 * t291) + pkin(1) * t53;
t15 = t177 * t24 + t180 * t36;
t103 = -t178 * g(3) + t181 * t116;
t210 = -t183 * pkin(2) + qJDD(2) * qJ(3) + t132 * t240 + t103;
t189 = pkin(3) * t249 + t136 * qJ(4) - qJD(2) * t143 - t210;
t41 = -t189 + t284;
t35 = t41 + t277;
t266 = t177 * t35;
t264 = t180 * t35;
t137 = -0.2e1 * t159 + t230;
t262 = pkin(1) * t137 - t293;
t141 = t242 * t184;
t261 = qJ(3) * t141;
t260 = qJ(3) * t149;
t259 = qJ(3) * t181;
t258 = qJ(6) * t177;
t257 = qJ(6) * t180;
t255 = t134 * t181;
t252 = t153 * t177;
t251 = t153 * t180;
t247 = t178 * t137;
t244 = pkin(1) * t141 + t133;
t229 = t178 * t94;
t226 = qJD(3) * t241;
t222 = t102 * t178 + t181 * t103;
t217 = t15 + t279;
t212 = t135 + t224;
t4 = -t180 * t14 + t177 * t15;
t5 = t177 * t14 + t180 * t15;
t88 = t247 + t255;
t11 = 0.2e1 * t238 + t217;
t206 = -qJ(4) * t33 + t272 * t31;
t71 = t164 + t210;
t203 = -qJ(4) * t77 + t297 * t31;
t202 = -qJ(4) * t65 + t297 * t47;
t201 = -qJ(4) * t291 + t297 * t53;
t200 = -t235 * t48 + t272 * t65;
t199 = -t235 * t54 + t272 * t291;
t198 = -t235 * t33 + t272 * t77;
t197 = -qJ(4) * t54 + t272 * t53 - t15;
t195 = t178 * t272 + t181 * t235 + pkin(1);
t10 = -qJ(6) * t84 - t283;
t193 = t10 + t298;
t192 = (t134 + t212) * t248 + t294;
t188 = t191 + t274;
t187 = t196 + 0.2e1 * t226 + t274;
t186 = t135 * qJ(3) + t188;
t42 = (t236 * t178 + t227) * qJD(1) + t190;
t185 = -t83 * pkin(5) - t121 * qJ(6) - t101 * t130 + qJDD(6) - t189 + t277;
t17 = t185 + t284;
t142 = (t171 - t172) * t184;
t105 = -t122 + t151;
t104 = t121 - t151;
t99 = t212 * t178;
t98 = t246 + t181 * (t183 - t250);
t97 = (t136 - t159) * t181;
t91 = t122 - t121;
t73 = (-t129 * t177 - t130 * t180) * t153;
t64 = pkin(5) * t70;
t60 = t130 * t251 - t177 * t84;
t59 = t129 * t252 - t180 * t83;
t57 = t178 * t123 + t181 * (-t129 * t180 + t130 * t177) * t153;
t56 = -t104 * t177 - t263;
t55 = -t105 * t180 - t296;
t40 = t229 + t181 * (-t130 * t252 - t180 * t84);
t39 = -t229 + t181 * (t129 * t251 + t177 * t83);
t38 = (qJ(3) * t239 + t221) * qJD(1) + t186;
t37 = -pkin(5) * t291 - qJ(6) * t79;
t32 = t177 * t65 - t180 * t291;
t26 = t178 * t70 + t181 * (t105 * t177 - t295);
t25 = t178 * t194 + t181 * (-t104 * t180 + t265);
t19 = t178 * t91 + t181 * (t177 * t291 + t180 * t65);
t16 = -qJ(6) * t90 + t17;
t12 = -pkin(5) * t65 + qJ(6) * t85 - t185 + t285;
t9 = pkin(5) * t10;
t7 = (t70 + t84) * qJ(6) + t283;
t6 = -pkin(5) * t77 + qJ(6) * t194 + t11;
t3 = -pkin(5) * t17 + qJ(6) * t11;
t2 = -t177 * t10 + t11 * t180;
t1 = t10 * t180 + t11 * t177;
t8 = [0, 0, 0, 0, 0, qJDD(1), t243, t215, 0, 0, t99, t88, t98, t97, -t281, 0, t115 * t181 + t262, -t178 * t115 - t294, t222 + t244, pkin(1) * t115 + pkin(7) * t222, t99, t98, -t88, 0, t281, t97, t181 * (pkin(2) * t137 + t187) + (t181 * t212 + t247) * qJ(3) + t262, t181 * (pkin(2) * t141 + t71) + (t72 + t261) * t178 + t244, pkin(2) * t255 + t178 * t187 + t192, pkin(7) * (t178 * t72 + t181 * t71) + (pkin(1) - t213) * (qJ(3) * t212 + t187), t97, t88, -t281, t99, t98, 0, t181 * (qJ(4) * t145 + t276 * t134) + (qJ(4) * t148 + qJD(1) * t221 + t188) * t178 + t192, -qJ(4) * t246 + t181 * (-qJ(3) * t224 - qJ(4) * t149 - t143 * t241 - t186 - 0.2e1 * t226) + t293 + t282 * t137, t181 * (qJ(4) * t230 + t189 + t285) - t133 + t218 * t141 + (-t261 + (-qJ(4) * t234 - g(3)) * t181 + t204 + (qJ(4) * qJDD(1) - t220 + 0.2e1 * t233) * t178) * t178, -t282 * t38 + t292 * (t178 * t42 + t181 * t41), t40, t19, t26, t39, t25, t57, t178 * (-t14 + t278) + t181 * (t202 - t266) + t270, t178 * t197 + t181 * (t201 - t264) + t269, t178 * t206 + t181 * (t203 + t4) + t271, t195 * t4 + t292 * (t178 * t5 + t181 * t35), t40, t19, t26, t39, t25, t57, t178 * (t193 + t278) + t181 * (t177 * t12 + t257 * t289 + t202) + t270, t178 * (t197 - t279 + t290) + t181 * (-t180 * t16 + t177 * t37 + t201) + t269, t178 * (t206 - t64) + t181 * (t177 * t6 - t180 * t7 + t203) + t271, t178 * (-qJ(4) * t2 + t9) + t181 * (-qJ(4) * t17 + t10 * t257 + t177 * t3) + pkin(7) * (t17 * t181 + t178 * t2) + t195 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t228, t142, t231, t228, t230, qJDD(2), -t102, -t103, 0, 0, -t228, t231, -t142, qJDD(2), -t230, t228, pkin(2) * t144 - t260 - t72, (-pkin(2) * t178 + t259) * qJDD(1), t71 + t286, -pkin(2) * t72 + qJ(3) * t71, t228, t142, t230, -t228, t231, qJDD(2), pkin(3) * t148 + t286 + t41, -t276 * t144 + t260 + t42, (t276 * t178 - t259) * qJDD(1), qJ(3) * t41 - t276 * t42, t60, t32, t55, t59, t56, t73, t200 + t264, t199 - t266, t198 - t5, -t235 * t5 + t272 * t35, t60, t32, t55, t59, t56, t73, -t180 * t12 + t258 * t289 + t200, -t177 * t16 - t180 * t37 + t199, -t177 * t7 - t180 * t6 + t198, t10 * t258 + t272 * t17 - t180 * t3 - t235 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t144, t231, -t148, t72, 0, 0, 0, 0, 0, 0, -t148, t144, -t231, t42, 0, 0, 0, 0, 0, 0, t48, t54, t33, t5, 0, 0, 0, 0, 0, 0, t48, t54, t33, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, -t137, -t141, t38, 0, 0, 0, 0, 0, 0, t47, t53, t31, t4, 0, 0, 0, 0, 0, 0, t47, t53, t31, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t91, t70, -t94, t194, t123, -t14, -t15, 0, 0, t94, t91, t70, -t94, t194, t123, t193, -t217 + t290, -t64, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t291, t77, t17;];
tauJ_reg  = t8;
