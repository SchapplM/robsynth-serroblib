% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRRR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 20:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:29:58
% EndTime: 2019-03-08 20:30:08
% DurationCPUTime: 4.12s
% Computational Cost: add. (3348->420), mult. (7755->612), div. (0->0), fcn. (6538->16), ass. (0->218)
t184 = sin(qJ(4));
t188 = cos(qJ(4));
t226 = pkin(4) * t184 - pkin(9) * t188;
t137 = t226 * qJD(4);
t179 = cos(pkin(12));
t185 = sin(qJ(2));
t178 = sin(pkin(6));
t267 = qJD(1) * t178;
t245 = t185 * t267;
t143 = t179 * t245;
t176 = sin(pkin(12));
t189 = cos(qJ(2));
t244 = t189 * t267;
t95 = t176 * t244 + t143;
t318 = t137 - t95;
t183 = sin(qJ(5));
t187 = cos(qJ(5));
t253 = t187 * qJD(4);
t265 = qJD(2) * t184;
t126 = t183 * t265 - t253;
t262 = qJD(4) * t183;
t128 = t187 * t265 + t262;
t182 = sin(qJ(6));
t186 = cos(qJ(6));
t216 = t126 * t182 - t186 * t128;
t78 = t186 * t126 + t128 * t182;
t317 = t216 * t78;
t266 = qJD(2) * t178;
t238 = qJD(1) * t266;
t251 = qJDD(1) * t178;
t316 = t185 * t251 + t189 * t238;
t257 = qJD(5) * t184;
t315 = -qJD(2) * t257 + qJDD(4);
t314 = t216 ^ 2 - t78 ^ 2;
t254 = qJD(6) * t186;
t255 = qJD(6) * t182;
t252 = qJD(2) * qJD(4);
t237 = t188 * t252;
t250 = t184 * qJDD(2);
t73 = qJD(5) * t253 + (t237 + t250) * t187 + t315 * t183;
t264 = qJD(2) * t188;
t74 = ((qJD(5) + t264) * qJD(4) + t250) * t183 - t315 * t187;
t11 = -t126 * t254 - t128 * t255 - t182 * t74 + t186 * t73;
t158 = -qJD(5) + t264;
t153 = -qJD(6) + t158;
t313 = -t153 * t78 + t11;
t273 = t189 * t179;
t281 = t178 * t185;
t102 = t176 * t281 - t178 * t273;
t181 = cos(pkin(6));
t157 = qJD(1) * t181 + qJD(3);
t139 = qJD(2) * pkin(2) + t244;
t92 = t176 * t139 + t143;
t89 = qJD(2) * pkin(8) + t92;
t58 = t184 * t157 + t188 * t89;
t55 = qJD(4) * pkin(9) + t58;
t214 = -pkin(4) * t188 - pkin(9) * t184 - pkin(3);
t142 = t176 * t245;
t91 = t139 * t179 - t142;
t61 = qJD(2) * t214 - t91;
t23 = t183 * t61 + t187 * t55;
t17 = -pkin(10) * t126 + t23;
t15 = t17 * t255;
t175 = qJ(5) + qJ(6);
t171 = sin(t175);
t172 = cos(t175);
t307 = t157 * t188 - t184 * t89;
t54 = -qJD(4) * pkin(4) - t307;
t32 = pkin(5) * t126 + t54;
t180 = cos(pkin(11));
t215 = t176 * t189 + t179 * t185;
t104 = t215 * t181;
t124 = t176 * t185 - t273;
t177 = sin(pkin(11));
t218 = t104 * t180 - t124 * t177;
t282 = t178 * t184;
t45 = -t180 * t282 + t188 * t218;
t217 = -t104 * t177 - t124 * t180;
t47 = t177 * t282 + t188 * t217;
t204 = t124 * t181;
t63 = -t177 * t215 - t180 * t204;
t66 = t177 * t204 - t180 * t215;
t103 = t215 * t178;
t87 = t103 * t188 + t181 * t184;
t312 = t32 * t78 - g(1) * (t171 * t66 - t172 * t47) - g(2) * (t171 * t63 - t172 * t45) - g(3) * (-t102 * t171 - t172 * t87) + t15;
t169 = t188 * qJDD(2);
t123 = t184 * t252 + qJDD(5) - t169;
t155 = t181 * qJDD(1) + qJDD(3);
t288 = t155 * t184;
t154 = t189 * t251;
t99 = qJDD(2) * pkin(2) - t185 * t238 + t154;
t52 = t176 * t99 + t316 * t179;
t50 = qJDD(2) * pkin(8) + t52;
t13 = qJDD(4) * pkin(9) + qJD(4) * t307 + t188 * t50 + t288;
t51 = -t316 * t176 + t179 * t99;
t29 = qJD(2) * t137 + qJDD(2) * t214 - t51;
t26 = t187 * t29;
t197 = -qJD(5) * t23 - t183 * t13 + t26;
t2 = pkin(5) * t123 - pkin(10) * t73 + t197;
t256 = qJD(5) * t187;
t249 = t187 * t13 + t183 * t29 + t61 * t256;
t258 = qJD(5) * t183;
t209 = t55 * t258 - t249;
t3 = -pkin(10) * t74 - t209;
t247 = -t182 * t3 + t186 * t2;
t22 = -t183 * t55 + t187 * t61;
t16 = -pkin(10) * t128 + t22;
t10 = -pkin(5) * t158 + t16;
t295 = t17 * t186;
t5 = t182 * t10 + t295;
t311 = t32 * t216 - g(1) * (-t171 * t47 - t172 * t66) - g(2) * (-t171 * t45 - t172 * t63) - g(3) * (t102 * t172 - t171 * t87) - t5 * qJD(6) + t247;
t196 = qJD(6) * t216 - t182 * t73 - t186 * t74;
t310 = t153 * t216 + t196;
t261 = qJD(4) * t184;
t276 = t183 * t188;
t161 = pkin(2) * t176 + pkin(8);
t285 = t161 * t183;
t98 = t179 * t244 - t142;
t309 = -t318 * t187 - t261 * t285 - t98 * t276;
t300 = pkin(2) * t179;
t121 = t214 - t300;
t274 = t187 * t188;
t308 = t121 * t256 + t318 * t183 - t98 * t274;
t130 = t182 * t187 + t183 * t186;
t106 = t130 * t184;
t201 = -t183 * t257 + t188 * t253;
t306 = qJD(5) + qJD(6);
t275 = t184 * t187;
t305 = t123 * t275 - t158 * t201;
t234 = qJD(6) * t10 + t3;
t304 = t182 * t2 + t186 * t234;
t206 = g(1) * t66 + g(2) * t63 - g(3) * t102;
t302 = (t158 * t161 + t55) * qJD(5) - t206;
t301 = pkin(9) + pkin(10);
t133 = t161 * t274;
t212 = pkin(5) * t184 - pkin(10) * t274;
t298 = -t212 * qJD(4) - (-t133 + (pkin(10) * t184 - t121) * t183) * qJD(5) + t309;
t260 = qJD(4) * t188;
t240 = t183 * t260;
t200 = t184 * t256 + t240;
t297 = (-t184 * t253 - t188 * t258) * t161 - t200 * pkin(10) + t308;
t120 = qJDD(6) + t123;
t277 = t183 * t184;
t34 = -t255 * t277 + (t306 * t275 + t240) * t186 + t201 * t182;
t296 = -t106 * t120 + t34 * t153;
t294 = t183 * t73;
t134 = t226 * qJD(2);
t293 = t183 * t134 + t187 * t307;
t129 = t182 * t183 - t186 * t187;
t203 = t129 * t188;
t292 = qJD(2) * t203 - t306 * t129;
t291 = (-t264 + t306) * t130;
t290 = t126 * t158;
t289 = t128 * t158;
t286 = t158 * t187;
t284 = t171 * t188;
t283 = t172 * t188;
t280 = t178 * t188;
t279 = t181 * t185;
t278 = t181 * t189;
t272 = qJDD(1) - g(3);
t269 = t183 * t121 + t133;
t173 = t184 ^ 2;
t268 = -t188 ^ 2 + t173;
t263 = qJD(4) * t126;
t259 = qJD(5) * t158;
t246 = qJD(5) * t301;
t243 = t183 * t264;
t242 = t158 * t262;
t233 = -t11 * t188 - t216 * t261;
t231 = t128 * t261 - t188 * t73;
t227 = -t58 + (-t243 + t258) * pkin(5);
t150 = t301 * t183;
t225 = pkin(10) * t243 - qJD(6) * t150 - t183 * t246 - t293;
t117 = t187 * t134;
t151 = t301 * t187;
t224 = qJD(2) * t212 + qJD(6) * t151 - t183 * t307 + t187 * t246 + t117;
t39 = t102 * t187 - t183 * t87;
t40 = t102 * t183 + t187 * t87;
t223 = -t182 * t40 + t186 * t39;
t222 = t182 * t39 + t186 * t40;
t109 = t187 * t121;
t60 = -pkin(10) * t275 + t109 + (-pkin(5) - t285) * t188;
t72 = -pkin(10) * t277 + t269;
t221 = t182 * t60 + t186 * t72;
t107 = t129 * t184;
t33 = -qJD(4) * t203 - t306 * t106;
t220 = t107 * t120 + t153 * t33;
t86 = t103 * t184 - t181 * t188;
t213 = -t155 * t188 + t157 * t261 + t184 * t50 + t89 * t260;
t210 = -t188 * t196 - t78 * t261;
t208 = -t123 * t183 + t158 * t256;
t207 = g(1) * (t177 * t280 - t184 * t217) + g(2) * (-t180 * t280 - t184 * t218) - g(3) * t86;
t205 = -g(1) * t217 - g(2) * t218 - g(3) * t103;
t14 = -qJDD(4) * pkin(4) + t213;
t199 = -pkin(9) * t123 - t158 * t54;
t162 = -pkin(3) - t300;
t88 = -qJD(2) * pkin(3) - t91;
t198 = -qJDD(4) * t161 + (qJD(2) * t162 + t88 + t98) * qJD(4);
t195 = pkin(9) * t259 - t14 - t207;
t190 = qJD(4) ^ 2;
t194 = -qJD(2) * t95 + t161 * t190 + t206 - t51 + (-pkin(3) + t162) * qJDD(2);
t192 = -g(1) * (-t177 * t278 - t180 * t185) - g(2) * (-t177 * t185 + t180 * t278) - g(3) * t178 * t189;
t191 = qJD(2) ^ 2;
t165 = -pkin(5) * t187 - pkin(4);
t146 = qJDD(4) * t188 - t184 * t190;
t145 = qJDD(4) * t184 + t188 * t190;
t112 = (pkin(5) * t183 + t161) * t184;
t97 = t124 * t266;
t96 = qJD(2) * t103;
t90 = pkin(5) * t200 + t161 * t260;
t38 = -qJD(4) * t86 - t188 * t97;
t37 = qJD(4) * t87 - t184 * t97;
t8 = qJD(5) * t39 + t183 * t96 + t187 * t38;
t7 = -qJD(5) * t40 - t183 * t38 + t187 * t96;
t6 = pkin(5) * t74 + t14;
t4 = t186 * t10 - t17 * t182;
t1 = [t272, 0 (qJDD(2) * t189 - t185 * t191) * t178 (-qJDD(2) * t185 - t189 * t191) * t178, -t102 * t51 + t103 * t52 + t155 * t181 - t91 * t96 - t92 * t97 - g(3), 0, 0, 0, 0, 0, -t102 * t169 - qJD(4) * t37 - qJDD(4) * t86 + (t102 * t261 - t188 * t96) * qJD(2), t102 * t250 - qJD(4) * t38 - qJDD(4) * t87 + (t102 * t260 + t184 * t96) * qJD(2), 0, 0, 0, 0, 0, t123 * t39 + t126 * t37 - t158 * t7 + t74 * t86, -t123 * t40 + t128 * t37 + t158 * t8 + t73 * t86, 0, 0, 0, 0, 0 -(-qJD(6) * t222 - t182 * t8 + t186 * t7) * t153 + t223 * t120 + t37 * t78 - t86 * t196 (qJD(6) * t223 + t182 * t7 + t186 * t8) * t153 - t222 * t120 - t37 * t216 + t86 * t11; 0, qJDD(2), t154 + t192, -g(1) * (t177 * t279 - t180 * t189) - g(2) * (-t177 * t189 - t180 * t279) - t272 * t281, t91 * t95 - t92 * t98 + (t52 * t176 + t51 * t179 + t192) * pkin(2), qJDD(2) * t173 + 0.2e1 * t184 * t237, 0.2e1 * t169 * t184 - 0.2e1 * t252 * t268, t145, t146, 0, t184 * t198 - t188 * t194, t184 * t194 + t188 * t198, t128 * t201 + t275 * t73 (-t126 * t187 - t128 * t183) * t260 + (-t294 - t187 * t74 + (t126 * t183 - t128 * t187) * qJD(5)) * t184, t231 + t305 (t74 + t242) * t188 + (t208 - t263) * t184, -t123 * t188 - t158 * t261, t109 * t123 + t309 * t158 + (t121 * t259 + t205) * t183 + (t161 * t263 - t26 + (qJD(4) * t54 + qJD(5) * t61 - t123 * t161 + t13) * t183 + t302 * t187) * t188 + (t22 * qJD(4) - t98 * t126 + t14 * t183 + t161 * t74 + t256 * t54) * t184, -t269 * t123 + t308 * t158 + t205 * t187 + ((t128 * t161 + t187 * t54) * qJD(4) - t302 * t183 + t249) * t188 + (-t54 * t258 - t98 * t128 + t14 * t187 + t161 * t73 + (-t161 * t286 - t23) * qJD(4)) * t184, -t107 * t11 - t216 * t33, -t106 * t11 - t107 * t196 + t216 * t34 - t33 * t78, -t220 + t233, t210 + t296, -t120 * t188 - t153 * t261 (-t182 * t72 + t186 * t60) * t120 - t247 * t188 + t90 * t78 - t112 * t196 + t6 * t106 + t32 * t34 - g(1) * (t171 * t217 + t283 * t66) - g(2) * (t171 * t218 + t283 * t63) - g(3) * (-t102 * t283 + t103 * t171) + (t4 * qJD(4) - t98 * t78) * t184 + (t297 * t182 + t298 * t186) * t153 + (t153 * t221 + t188 * t5) * qJD(6), -t221 * t120 + (-t15 + t304) * t188 - t90 * t216 + t112 * t11 - t6 * t107 + t32 * t33 - g(1) * (t172 * t217 - t284 * t66) - g(2) * (t172 * t218 - t284 * t63) - g(3) * (t102 * t284 + t103 * t172) + (-t5 * qJD(4) + t216 * t98) * t184 + ((qJD(6) * t60 + t297) * t186 + (-qJD(6) * t72 - t298) * t182) * t153; 0, 0, 0, 0, -g(3) * t181 + (-g(1) * t177 + g(2) * t180) * t178 + t155, 0, 0, 0, 0, 0, t146, -t145, 0, 0, 0, 0, 0 (-t74 + t242) * t188 + (t208 + t263) * t184, t231 - t305, 0, 0, 0, 0, 0, -t210 + t296, t220 + t233; 0, 0, 0, 0, 0, -t184 * t191 * t188, t268 * t191, t250, t169, qJDD(4), qJD(4) * t58 - t265 * t88 - t207 - t213, g(1) * t47 + g(2) * t45 + g(3) * t87 - t288 + (-qJD(2) * t88 - t50) * t188, -t128 * t286 + t294 (t73 + t290) * t187 + (-t74 + t289) * t183 (-t128 * t184 + t158 * t274) * qJD(2) - t208, t158 * t258 + t123 * t187 + (t126 * t184 - t158 * t276) * qJD(2), t158 * t265, -t22 * t265 - pkin(4) * t74 + t117 * t158 - t58 * t126 + (-t158 * t307 + t199) * t183 + t195 * t187, -pkin(4) * t73 - t58 * t128 - t293 * t158 - t195 * t183 + t199 * t187 + t23 * t265, t11 * t130 - t216 * t292, -t11 * t129 + t130 * t196 + t216 * t291 - t292 * t78, t120 * t130 - t292 * t153 + t216 * t265, -t120 * t129 + t291 * t153 + t78 * t265, t153 * t265 (-t150 * t186 - t151 * t182) * t120 - t165 * t196 + t6 * t129 - t4 * t265 + t227 * t78 + t291 * t32 + (t182 * t225 + t186 * t224) * t153 - t207 * t172 -(-t150 * t182 + t151 * t186) * t120 + t165 * t11 + t6 * t130 + t5 * t265 - t227 * t216 + t292 * t32 + (-t182 * t224 + t186 * t225) * t153 + t207 * t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t128 * t126, -t126 ^ 2 + t128 ^ 2, t73 - t290, -t289 - t74, t123, -t23 * t158 - t54 * t128 - g(1) * (-t183 * t47 - t187 * t66) - g(2) * (-t183 * t45 - t187 * t63) - g(3) * t39 + t197, -t22 * t158 + t54 * t126 - g(1) * (t183 * t66 - t187 * t47) - g(2) * (t183 * t63 - t187 * t45) + g(3) * t40 + t209, -t317, t314, t313, t310, t120 (-t16 * t182 - t295) * t153 + (t186 * t120 - t128 * t78 + t153 * t255) * pkin(5) + t311 (t17 * t153 - t2) * t182 + (-t16 * t153 - t234) * t186 + (-t182 * t120 + t128 * t216 + t153 * t254) * pkin(5) + t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t317, t314, t313, t310, t120, -t153 * t5 + t311, -t153 * t4 - t304 + t312;];
tau_reg  = t1;
