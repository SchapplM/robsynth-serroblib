% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tau_reg [6x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PPRRRP2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_invdynJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:58:12
% EndTime: 2019-03-08 18:58:21
% DurationCPUTime: 3.92s
% Computational Cost: add. (4692->426), mult. (11859->594), div. (0->0), fcn. (11040->14), ass. (0->213)
t154 = sin(qJ(4));
t157 = cos(qJ(4));
t209 = pkin(4) * t154 - pkin(10) * t157;
t126 = t209 * qJD(4);
t271 = cos(pkin(6));
t136 = t271 * qJD(1) + qJD(2);
t151 = sin(pkin(6));
t149 = sin(pkin(12));
t155 = sin(qJ(3));
t158 = cos(qJ(3));
t152 = cos(pkin(12));
t270 = cos(pkin(7));
t217 = t152 * t270;
t184 = t149 * t158 + t155 * t217;
t176 = t184 * t151;
t150 = sin(pkin(7));
t258 = t150 * t155;
t79 = qJD(1) * t176 + t136 * t258;
t304 = t126 - t79;
t133 = t271 * qJDD(1) + qJDD(2);
t213 = t151 * t217;
t199 = qJD(1) * t213;
t247 = qJD(3) * t155;
t229 = t150 * t247;
t249 = qJD(1) * t151;
t230 = t149 * t249;
t259 = t149 * t155;
t232 = t151 * t259;
t245 = qJD(3) * t158;
t161 = -t158 * (qJDD(1) * t213 + t133 * t150) + qJDD(1) * t232 + t136 * t229 + t199 * t247 + t230 * t245;
t269 = cos(pkin(11));
t208 = t271 * t269;
t268 = sin(pkin(11));
t106 = t149 * t208 + t268 * t152;
t169 = t268 * t149 - t152 * t208;
t219 = t151 * t269;
t297 = t150 * t219 + t169 * t270;
t62 = t106 * t155 + t297 * t158;
t207 = t271 * t268;
t107 = -t149 * t207 + t269 * t152;
t170 = t269 * t149 + t152 * t207;
t218 = t151 * t268;
t296 = -t150 * t218 + t170 * t270;
t64 = t107 * t155 + t296 * t158;
t211 = t158 * t217;
t220 = t150 * t271;
t212 = t158 * t220;
t86 = -t151 * t211 - t212 + t232;
t187 = g(1) * t64 + g(2) * t62 + g(3) * t86;
t303 = t79 * qJD(3) - t161 + t187;
t246 = qJD(3) * t157;
t302 = qJD(5) - t246;
t238 = qJD(5) * t154;
t301 = qJD(3) * t238 - qJDD(4);
t153 = sin(qJ(5));
t156 = cos(qJ(5));
t233 = t154 * qJDD(3);
t81 = ((qJD(5) + t246) * qJD(4) + t233) * t153 + t301 * t156;
t201 = pkin(4) * t157 + pkin(10) * t154 + pkin(3);
t239 = qJD(5) * t153;
t256 = t153 * t157;
t78 = -t155 * t230 + t158 * (t136 * t150 + t199);
t300 = -t304 * t156 - t201 * t239 - t78 * t256;
t237 = qJD(5) * t156;
t254 = t156 * t157;
t299 = t304 * t153 - t201 * t237 - t78 * t254;
t231 = t151 * t152 * t150;
t102 = -qJD(1) * t231 + t270 * t136;
t77 = qJD(3) * pkin(9) + t79;
t298 = t102 * t157 - t154 * t77;
t295 = qJDD(1) * t184;
t63 = t106 * t158 - t297 * t155;
t65 = t107 * t158 - t296 * t155;
t178 = t271 * t270 - t231;
t87 = t155 * t220 + t176;
t66 = t87 * t154 - t157 * t178;
t88 = t169 * t150 - t270 * t219;
t89 = t170 * t150 + t270 * t218;
t294 = g(3) * t66 - g(2) * (-t154 * t63 + t157 * t88) - g(1) * (-t154 * t65 + t157 * t89);
t235 = t156 * qJD(4);
t248 = qJD(3) * t154;
t117 = t153 * t248 - t235;
t242 = qJD(4) * t153;
t119 = t156 * t248 + t242;
t48 = -qJD(4) * pkin(4) - t298;
t26 = pkin(5) * t117 - qJ(6) * t119 + t48;
t144 = t157 * qJDD(3);
t234 = qJD(3) * qJD(4);
t115 = t154 * t234 + qJDD(5) - t144;
t283 = pkin(10) * t115;
t293 = -t26 * t302 + t283;
t291 = t119 ^ 2;
t284 = pkin(5) * t115;
t236 = qJD(5) * t157;
t241 = qJD(4) * t154;
t282 = qJ(6) * t241 - qJD(6) * t157 + (-t153 * t236 - t154 * t235) * pkin(9) + t299;
t281 = -pkin(5) * t241 + (-t153 * t241 + t156 * t236) * pkin(9) + t300;
t53 = t154 * t102 + t157 * t77;
t280 = pkin(10) * qJD(5);
t279 = qJD(3) * pkin(3);
t49 = qJD(4) * pkin(10) + t53;
t69 = -qJD(3) * t201 - t78;
t17 = t153 * t69 + t156 * t49;
t15 = qJ(6) * t302 + t17;
t278 = t302 * t15;
t277 = t302 * t17;
t224 = t157 * t234;
t80 = -qJD(5) * t235 + (-t224 - t233) * t156 + t301 * t153;
t276 = t153 * t80;
t275 = t78 * t117;
t274 = t78 * t119;
t205 = pkin(5) * t153 - qJ(6) * t156;
t273 = -qJD(6) * t153 + t205 * t302 - t53;
t125 = t209 * qJD(3);
t272 = t153 * t125 + t156 * t298;
t267 = qJ(6) * t115;
t101 = -qJDD(1) * t231 + t270 * t133;
t265 = t101 * t154;
t263 = t117 * t302;
t262 = t119 * t117;
t261 = t119 * t302;
t260 = t119 * t156;
t257 = t150 * t158;
t255 = t156 * t115;
t16 = -t153 * t49 + t156 * t69;
t253 = qJD(6) - t16;
t251 = pkin(9) * t254 - t153 * t201;
t147 = t154 ^ 2;
t250 = -t157 ^ 2 + t147;
t244 = qJD(4) * t117;
t243 = qJD(4) * t119;
t240 = qJD(4) * t157;
t228 = t150 * t245;
t227 = t302 * t242;
t226 = t302 * t235;
t225 = t302 * t239;
t223 = t158 * t234;
t25 = qJD(3) * t126 - qJDD(3) * t201 + t161;
t183 = t211 - t259;
t42 = qJDD(3) * pkin(9) + (t133 * t155 + t136 * t245) * t150 + (qJD(1) * qJD(3) * t183 + t295) * t151;
t8 = qJDD(4) * pkin(10) + qJD(4) * t298 + t157 * t42 + t265;
t221 = t153 * t8 - t156 * t25 + t49 * t237 + t69 * t239;
t206 = pkin(5) * t156 + qJ(6) * t153;
t14 = -pkin(5) * t302 + t253;
t204 = t14 * t156 - t15 * t153;
t67 = t154 * t178 + t87 * t157;
t36 = t153 * t86 + t156 * t67;
t35 = t153 * t67 - t86 * t156;
t202 = -t101 * t157 + t102 * t241 + t154 * t42 + t77 * t240;
t160 = qJD(3) ^ 2;
t200 = qJDD(3) * t158 - t155 * t160;
t198 = pkin(4) + t206;
t197 = pkin(9) + t205;
t195 = t153 * t25 + t156 * t8 + t69 * t237 - t49 * t239;
t109 = t154 * t270 + t157 * t258;
t93 = t109 * t153 + t156 * t257;
t94 = t109 * t156 - t153 * t257;
t194 = -t153 * t115 - t237 * t302;
t193 = -t225 + t255;
t18 = -t63 * t156 - t62 * t256;
t20 = -t65 * t156 - t64 * t256;
t50 = -t87 * t156 - t86 * t256;
t191 = g(1) * t20 + g(2) * t18 + g(3) * t50;
t19 = t153 * t63 - t62 * t254;
t21 = t153 * t65 - t64 * t254;
t51 = t153 * t87 - t86 * t254;
t190 = -g(1) * t21 - g(2) * t19 - g(3) * t51;
t32 = t154 * t88 + t157 * t63;
t34 = t154 * t89 + t157 * t65;
t188 = g(1) * t34 + g(2) * t32 + g(3) * t67;
t82 = (t183 * t151 + t212) * qJD(3);
t27 = qJD(4) * t67 + t82 * t154;
t28 = -qJD(4) * t66 + t82 * t157;
t83 = t87 * qJD(3);
t4 = qJD(5) * t36 + t153 * t28 - t156 * t83;
t186 = -t115 * t35 + t27 * t117 - t302 * t4 + t66 * t81;
t108 = t154 * t258 - t157 * t270;
t91 = -qJD(4) * t108 + t157 * t228;
t47 = qJD(5) * t94 + t153 * t91 - t156 * t229;
t92 = qJD(4) * t109 + t154 * t228;
t185 = t108 * t81 - t115 * t93 + t92 * t117 - t302 * t47;
t9 = -qJDD(4) * pkin(4) + t202;
t182 = t302 * t48 - t283;
t76 = -t78 - t279;
t179 = -pkin(9) * qJDD(4) + (t76 + t78 - t279) * qJD(4);
t175 = -g(1) * t218 + g(2) * t219 - g(3) * t271;
t5 = -qJD(5) * t35 + t153 * t83 + t156 * t28;
t174 = t115 * t36 - t119 * t27 + t302 * t5 + t66 * t80;
t46 = -qJD(5) * t93 + t153 * t229 + t156 * t91;
t173 = t108 * t80 + t115 * t94 - t119 * t92 + t302 * t46;
t10 = t153 * t32 - t62 * t156;
t12 = t153 * t34 - t64 * t156;
t172 = g(1) * t12 + g(2) * t10 + g(3) * t35 - t221;
t171 = -t280 * t302 + t294;
t3 = pkin(5) * t81 + qJ(6) * t80 - qJD(6) * t119 + t9;
t168 = t171 - t3;
t11 = t153 * t62 + t156 * t32;
t13 = t153 * t64 + t156 * t34;
t166 = -g(1) * t13 - g(2) * t11 - g(3) * t36 + t195;
t163 = t119 * t26 + qJDD(6) - t172;
t159 = qJD(4) ^ 2;
t162 = 0.2e1 * qJDD(3) * pkin(3) - pkin(9) * t159 + t303;
t104 = t197 * t154;
t98 = t201 * t156 + (pkin(9) * t153 + pkin(5)) * t157;
t97 = -qJ(6) * t157 + t251;
t90 = pkin(5) * t119 + qJ(6) * t117;
t71 = (t206 * qJD(5) - qJD(6) * t156) * t154 + t197 * t240;
t61 = -t80 + t263;
t30 = -pkin(5) * t248 - t125 * t156 + t153 * t298;
t29 = qJ(6) * t248 + t272;
t2 = qJDD(6) + t221 - t284;
t1 = qJD(6) * t302 + t195 + t267;
t6 = [qJDD(1) - g(3), t133 * t271 - g(3) + (t149 ^ 2 + t152 ^ 2) * t151 ^ 2 * qJDD(1), 0, -qJD(3) * t83 - qJDD(3) * t86, -qJD(3) * t82 - qJDD(3) * t87, 0, 0, 0, 0, 0, -t86 * t144 - qJD(4) * t27 - qJDD(4) * t66 + (-t157 * t83 + t86 * t241) * qJD(3), t86 * t233 - qJD(4) * t28 - qJDD(4) * t67 + (t154 * t83 + t240 * t86) * qJD(3), 0, 0, 0, 0, 0, t186, -t174, t186, -t117 * t5 + t119 * t4 - t35 * t80 - t36 * t81, t174, t1 * t36 + t14 * t4 + t15 * t5 + t2 * t35 + t26 * t27 + t3 * t66 - g(3); 0, t175 + t133, 0, t200 * t150 (-qJDD(3) * t155 - t158 * t160) * t150, 0, 0, 0, 0, 0, -qJD(4) * t92 - qJDD(4) * t108 + (-t154 * t223 + t157 * t200) * t150, -qJD(4) * t91 - qJDD(4) * t109 + (-t154 * t200 - t157 * t223) * t150, 0, 0, 0, 0, 0, t185, -t173, t185, -t117 * t46 + t119 * t47 - t80 * t93 - t81 * t94, t173, t1 * t94 + t3 * t108 + t14 * t47 + t15 * t46 + t2 * t93 + t26 * t92 + t175; 0, 0, qJDD(3), t303, -t133 * t258 + g(1) * t65 + g(2) * t63 + g(3) * t87 - t151 * t295 + (-t136 * t257 - t183 * t249 + t78) * qJD(3), qJDD(3) * t147 + 0.2e1 * t154 * t224, 0.2e1 * t154 * t144 - 0.2e1 * t250 * t234, qJDD(4) * t154 + t157 * t159, qJDD(4) * t157 - t154 * t159, 0, t154 * t179 + t157 * t162, -t154 * t162 + t157 * t179, -t154 * t156 * t80 + (-t153 * t238 + t157 * t235) * t119 (-t117 * t156 - t119 * t153) * t240 + (t276 - t156 * t81 + (t117 * t153 - t260) * qJD(5)) * t154 (t80 + t226) * t157 + (t193 + t243) * t154 (t81 - t227) * t157 + (t194 - t244) * t154, -t115 * t157 + t241 * t302, -t201 * t255 - t300 * t302 + (t48 * t242 + (t194 + t244) * pkin(9) + t221) * t157 + (t48 * t237 + t16 * qJD(4) - t275 + t9 * t153 + (t81 + t227) * pkin(9)) * t154 + t190, -t251 * t115 - t299 * t302 + (t48 * t235 + (t225 + t243) * pkin(9) + t195) * t157 + (-t48 * t239 - t17 * qJD(4) - t274 + t9 * t156 + (-t80 + t226) * pkin(9)) * t154 + t191, t104 * t81 - t115 * t98 + t117 * t71 + (t242 * t26 + t2) * t157 - t281 * t302 + (-qJD(4) * t14 + t153 * t3 + t237 * t26 - t275) * t154 + t190, -t80 * t98 - t81 * t97 + t281 * t119 - t282 * t117 + t204 * t240 + (-t1 * t153 + t156 * t2 + (-t14 * t153 - t15 * t156) * qJD(5) + t187) * t154, t104 * t80 + t115 * t97 - t119 * t71 + (-t235 * t26 - t1) * t157 + t282 * t302 + (qJD(4) * t15 - t156 * t3 + t239 * t26 + t274) * t154 - t191, t1 * t97 + t3 * t104 + t2 * t98 - g(1) * (pkin(5) * t21 + pkin(9) * t65 + qJ(6) * t20) - g(2) * (pkin(5) * t19 + pkin(9) * t63 + qJ(6) * t18) - g(3) * (pkin(5) * t51 + pkin(9) * t87 + qJ(6) * t50) + (-t154 * t78 + t71) * t26 + t282 * t15 + t281 * t14 + t187 * t201; 0, 0, 0, 0, 0, -t154 * t160 * t157, t250 * t160, t233, t144, qJDD(4), qJD(4) * t53 - t248 * t76 - t202 + t294, -t265 + (-qJD(3) * t76 - t42) * t157 + t188, t260 * t302 - t276 (-t80 - t263) * t156 + (-t261 - t81) * t153 (-t119 * t154 - t254 * t302) * qJD(3) - t194 (t117 * t154 + t256 * t302) * qJD(3) + t193, -t302 * t248, -t16 * t248 - pkin(4) * t81 - t53 * t117 + (t298 * t302 + t182) * t153 + (-t9 - (t125 + t280) * t302 + t294) * t156, pkin(4) * t80 + t272 * t302 + t17 * t248 - t53 * t119 + t182 * t156 + (-t171 + t9) * t153, t273 * t117 + t14 * t248 - t293 * t153 + t168 * t156 - t198 * t81 + t30 * t302, t117 * t29 - t119 * t30 + (t1 + t302 * t14 + (qJD(5) * t119 - t81) * pkin(10)) * t156 + (t2 - t278 + (qJD(5) * t117 - t80) * pkin(10)) * t153 - t188, -t273 * t119 - t15 * t248 + t168 * t153 + t293 * t156 - t198 * t80 - t29 * t302, -t14 * t30 - t15 * t29 + t273 * t26 + (qJD(5) * t204 + t1 * t156 + t2 * t153 - t188) * pkin(10) + (-t3 + t294) * t198; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t262, -t117 ^ 2 + t291, t61, -t81 + t261, t115, -t119 * t48 + t172 + t277, t117 * t48 + t16 * t302 - t166, -t117 * t90 - t163 + t277 + 0.2e1 * t284, pkin(5) * t80 - qJ(6) * t81 + (t15 - t17) * t119 + (t14 - t253) * t117, 0.2e1 * t267 - t117 * t26 + t119 * t90 - (-0.2e1 * qJD(6) + t16) * t302 + t166, t1 * qJ(6) - t2 * pkin(5) - t26 * t90 - t14 * t17 - g(1) * (-pkin(5) * t12 + qJ(6) * t13) - g(2) * (-pkin(5) * t10 + qJ(6) * t11) - g(3) * (-pkin(5) * t35 + qJ(6) * t36) + t253 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115 + t262, t61, -t302 ^ 2 - t291, t163 - t278 - t284;];
tau_reg  = t6;
