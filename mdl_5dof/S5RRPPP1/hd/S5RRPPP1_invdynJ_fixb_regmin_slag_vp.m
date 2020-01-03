% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPPP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d1,d2,theta3]';
% 
% Output:
% tau_reg [5x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:24:33
% EndTime: 2019-12-31 19:24:42
% DurationCPUTime: 3.72s
% Computational Cost: add. (3861->453), mult. (10838->584), div. (0->0), fcn. (7982->8), ass. (0->221)
t179 = sin(qJ(2));
t181 = cos(qJ(2));
t251 = t181 * qJDD(1);
t254 = qJD(1) * qJD(2);
t313 = t179 * t254 - t251;
t176 = sin(pkin(5));
t178 = cos(pkin(5));
t180 = sin(qJ(1));
t272 = t180 * t178;
t182 = cos(qJ(1));
t273 = t179 * t182;
t123 = t176 * t273 + t272;
t281 = t176 * t179;
t265 = t181 * pkin(2) + qJ(3) * t281;
t312 = -pkin(1) - t265;
t175 = sin(pkin(8));
t271 = t181 * t175;
t177 = cos(pkin(8));
t275 = t179 * t177;
t116 = t178 * t275 + t271;
t103 = t116 * qJD(2);
t270 = t181 * t177;
t245 = t178 * t270;
t252 = t179 * qJDD(1);
t253 = qJDD(2) * t176;
t210 = -qJDD(1) * t245 + t175 * t252 - t177 * t253;
t49 = qJD(1) * t103 + t210;
t262 = qJD(1) * t181;
t237 = t178 * t262;
t263 = qJD(1) * t179;
t238 = t175 * t263;
t260 = qJD(2) * t177;
t72 = -t176 * t260 - t177 * t237 + t238;
t311 = t49 * qJ(5) + t72 * qJD(5);
t257 = qJD(3) * t177;
t279 = t176 * t181;
t248 = qJ(3) * t279;
t300 = pkin(2) * t179;
t211 = -t248 + t300;
t125 = t211 * qJD(1);
t232 = qJ(3) * t178 + pkin(7);
t219 = qJD(1) * t232;
t126 = t179 * t219;
t127 = t181 * t219;
t283 = t175 * t178;
t284 = t175 * t176;
t38 = t125 * t284 - t177 * t126 - t127 * t283;
t288 = -t178 * qJD(4) + t38 + (qJ(4) * t263 - t257) * t176;
t131 = -t178 * qJD(2) + t176 * t262;
t115 = t178 * t271 + t275;
t231 = t178 * t238;
t234 = t181 * t254;
t50 = -qJD(2) * t231 + qJDD(1) * t115 + t175 * t253 + t177 * t234;
t261 = qJD(2) * t176;
t133 = t237 + t261;
t306 = t176 ^ 2 + t178 ^ 2;
t56 = t177 * t133 - t306 * t238;
t310 = -t56 * t131 + t50;
t128 = t131 ^ 2;
t73 = qJD(1) * t115 + t175 * t261;
t309 = -t73 ^ 2 - t128;
t135 = t232 * t179;
t308 = (-t135 * t178 + t176 * t312) * t177;
t294 = pkin(3) + qJ(5);
t85 = t178 * qJDD(2) + t313 * t176;
t307 = t50 * pkin(4) - t294 * t85;
t268 = t182 * t178;
t274 = t179 * t180;
t122 = t176 * t274 - t268;
t305 = -g(1) * t123 - g(2) * t122 + g(3) * t279;
t259 = qJD(2) * t179;
t256 = qJD(3) * t179;
t86 = t211 * qJD(2) - t176 * t256;
t218 = qJD(2) * t232;
t87 = t181 * t178 * qJD(3) - t179 * t218;
t88 = -t178 * t256 - t181 * t218;
t28 = t177 * t87 + t88 * t283 + t86 * t284;
t17 = -t176 * (qJ(4) * t259 - qJD(4) * t181) - t28;
t109 = qJD(2) * pkin(2) - t126;
t110 = t312 * qJD(1);
t92 = pkin(7) * t262 + qJ(3) * t133;
t31 = -t175 * t92 + (t109 * t178 + t110 * t176) * t177;
t44 = t133 * qJD(3) - t313 * pkin(7) + (-t178 * t313 + t253) * qJ(3);
t53 = qJDD(2) * pkin(2) + t88 * qJD(1) - qJDD(1) * t135;
t54 = qJD(1) * t86 + qJDD(1) * t312;
t7 = -t175 * t44 + (t176 * t54 + t178 * t53) * t177;
t301 = t85 * pkin(3);
t298 = g(1) * t180;
t296 = g(3) * t179;
t295 = t72 * t73;
t293 = t177 * t86;
t55 = t306 * t177 * t263 + t175 * t133;
t292 = t55 * t131;
t101 = t116 * qJD(1);
t102 = t177 * t262 - t231;
t57 = t178 * t125 + t176 * t127;
t208 = -t102 * qJ(4) + t57;
t255 = qJD(4) * t175;
t290 = -t294 * t101 + (-qJD(5) * t177 - t255) * t176 - t208;
t106 = t175 * t126;
t277 = t177 * t178;
t222 = t127 * t277 - t106;
t236 = t294 * t179;
t258 = qJD(3) * t175;
t239 = t176 * t258;
t285 = t125 * t177;
t289 = -t178 * qJD(5) + t239 - t102 * pkin(4) - (-qJD(1) * t236 - t285) * t176 - t222;
t287 = pkin(4) * t101 - t288;
t286 = qJD(4) * t73;
t282 = t176 * t177;
t280 = t176 * t180;
t278 = t176 * t182;
t276 = t179 * t175;
t269 = t181 * t182;
t136 = t232 * t181;
t120 = t175 * t136;
t267 = pkin(3) * t279 + t120;
t119 = pkin(2) * t283 + qJ(3) * t282;
t266 = t182 * pkin(7) + qJ(3) * t268;
t173 = t179 ^ 2;
t264 = -t181 ^ 2 + t173;
t8 = t177 * t44 + t53 * t283 + t54 * t284;
t32 = t109 * t283 + t110 * t284 + t177 * t92;
t250 = pkin(2) * t274;
t249 = pkin(2) * t273;
t247 = t178 * t276;
t246 = t175 * t274;
t244 = t175 * t273;
t242 = t177 * t269;
t43 = -t135 * t283 + t177 * t136 + t284 * t312;
t241 = -pkin(2) * t177 - pkin(3);
t233 = -qJ(4) * t175 - pkin(2);
t45 = -t176 * t88 + t178 * t86;
t58 = t176 * t135 + t178 * t312;
t62 = t180 * t271 + (t179 * t272 + t278) * t177;
t64 = t116 * t182 - t177 * t280;
t230 = g(1) * t62 - g(2) * t64;
t63 = -t175 * t278 - t178 * t246 + t180 * t270;
t65 = t175 * t280 - t178 * t244 + t242;
t229 = g(1) * t63 - g(2) * t65;
t75 = t175 * t87;
t228 = -t88 * t277 + t75;
t148 = t180 * t248;
t93 = -t180 * t245 + t246;
t94 = t115 * t180;
t227 = -t94 * pkin(3) - t93 * qJ(4) + t148;
t150 = t182 * t248;
t95 = -t178 * t242 + t244;
t96 = t115 * t182;
t226 = -t96 * pkin(3) - t95 * qJ(4) + t150;
t5 = -t85 * qJ(4) + t131 * qJD(4) - t8;
t225 = g(1) * t122 - g(2) * t123;
t224 = g(1) * t182 + g(2) * t180;
t223 = -g(2) * t182 + t298;
t221 = -t55 * t73 + t56 * t72;
t97 = -t178 * qJ(4) - t119;
t52 = -t176 * t109 + t178 * t110 + qJD(3);
t20 = -t176 * t53 + t178 * t54 + qJDD(3);
t217 = t182 * pkin(1) + pkin(2) * t269 + t180 * pkin(7) + t123 * qJ(3);
t24 = qJ(4) * t131 - t32;
t214 = -t63 * pkin(3) - t62 * qJ(4) + t266;
t213 = pkin(4) * t279 - t300;
t117 = -t247 + t270;
t212 = t117 * pkin(3) + qJ(4) * t116 + t265;
t209 = -0.2e1 * pkin(1) * t254 - pkin(7) * qJDD(2);
t207 = -t115 * qJ(4) + t58;
t204 = g(1) * t95 + g(2) * t93 - g(3) * t116;
t203 = g(1) * t96 + g(2) * t94 - g(3) * t117;
t36 = qJ(4) * t279 - t43;
t201 = -t73 * qJ(4) + t52;
t3 = -pkin(4) * t49 + qJDD(5) - t5;
t200 = t65 * pkin(3) + qJ(4) * t64 + t217;
t104 = -qJD(2) * t247 + t181 * t260;
t199 = -t104 * qJ(4) - t115 * qJD(4) + t45;
t198 = t312 * t298;
t197 = -t224 * t181 - t296;
t196 = qJDD(4) - t7;
t195 = qJD(4) - t31;
t183 = qJD(2) ^ 2;
t194 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t183 + t223;
t184 = qJD(1) ^ 2;
t193 = pkin(1) * t184 - pkin(7) * qJDD(1) + t224;
t191 = t85 - t295;
t4 = pkin(3) * t49 - qJ(4) * t50 + t20 - t286;
t190 = -t116 * t254 - t210;
t114 = -t245 + t276;
t189 = -g(1) * t64 - g(2) * t62 - g(3) * t114 + t196;
t188 = t4 + t305;
t187 = t49 - t292;
t185 = -t131 * t72 + t50;
t157 = qJ(3) * t284;
t118 = pkin(2) * t277 - t157;
t99 = (-pkin(3) * t177 + t233) * t176;
t98 = t241 * t178 + t157;
t70 = (-t294 * t177 + t233) * t176;
t69 = pkin(4) * t282 - t97;
t59 = pkin(4) * t284 + t157 + (-qJ(5) + t241) * t178;
t42 = -t120 + t308;
t39 = t267 - t308;
t37 = t106 + (t125 * t176 - t127 * t178) * t177;
t35 = pkin(3) * t114 + t207;
t34 = (-pkin(3) * t263 - t285) * t176 + t222;
t30 = pkin(3) * t101 + t208;
t29 = -pkin(4) * t114 - t36;
t27 = -t75 + (t176 * t86 + t178 * t88) * t177;
t26 = t294 * t114 + t207;
t25 = t135 * t277 + t115 * pkin(4) + (qJ(5) * t181 - t177 * t312) * t176 + t267;
t23 = t131 * pkin(3) + t195;
t21 = (-pkin(3) * t259 - t293) * t176 + t228;
t16 = pkin(3) * t72 + t201;
t15 = pkin(3) * t103 + t199;
t14 = -pkin(4) * t72 + qJD(5) - t24;
t13 = -t103 * pkin(4) - t17;
t12 = t104 * pkin(4) + (-qJD(2) * t236 + qJD(5) * t181 - t293) * t176 + t228;
t11 = t294 * t72 + t201;
t10 = t73 * pkin(4) + t294 * t131 + t195;
t9 = t114 * qJD(5) + t294 * t103 + t199;
t6 = t196 - t301;
t2 = t131 * qJD(5) + t196 + t307;
t1 = t4 + t311;
t18 = [qJDD(1), t223, t224, qJDD(1) * t173 + 0.2e1 * t179 * t234, 0.2e1 * t179 * t251 - 0.2e1 * t264 * t254, qJDD(2) * t179 + t181 * t183, qJDD(2) * t181 - t179 * t183, 0, t209 * t179 + t194 * t181, -t194 * t179 + t209 * t181, t52 * t103 + t20 * t114 - t27 * t131 + t42 * t85 + t45 * t72 + t58 * t49 + (-t181 * t7 + t31 * t259) * t176 + t229, t52 * t104 + t20 * t115 + t28 * t131 - t43 * t85 + t45 * t73 + t58 * t50 + (t181 * t8 - t32 * t259) * t176 - t230, -t103 * t32 - t104 * t31 - t114 * t8 - t115 * t7 - t27 * t73 - t28 * t72 - t42 * t50 - t43 * t49 + t225, -g(1) * t266 - g(2) * t217 + t20 * t58 + t31 * t27 + t32 * t28 + t7 * t42 + t8 * t43 + t52 * t45 - t198, t103 * t24 + t104 * t23 + t114 * t5 + t115 * t6 + t17 * t72 + t21 * t73 + t36 * t49 + t39 * t50 + t225, -t16 * t103 - t4 * t114 - t21 * t131 - t15 * t72 - t35 * t49 + t39 * t85 + (-t181 * t6 + t23 * t259) * t176 - t229, -t16 * t104 - t4 * t115 + t17 * t131 - t15 * t73 - t35 * t50 - t36 * t85 + (t181 * t5 - t24 * t259) * t176 + t230, -g(1) * t214 - g(2) * t200 + t16 * t15 + t24 * t17 + t23 * t21 + t4 * t35 + t5 * t36 + t6 * t39 - t198, t10 * t104 - t103 * t14 - t114 * t3 + t115 * t2 + t12 * t73 - t13 * t72 + t25 * t50 - t29 * t49 + t225, -t1 * t115 - t11 * t104 - t13 * t131 - t26 * t50 + t29 * t85 - t9 * t73 + (t14 * t259 - t181 * t3) * t176 + t230, t1 * t114 + t11 * t103 + t12 * t131 - t25 * t85 + t26 * t49 + t9 * t72 + (-t10 * t259 + t181 * t2) * t176 + t229, t1 * t26 + t11 * t9 + t2 * t25 + t10 * t12 + t3 * t29 + t14 * t13 - g(1) * (-t122 * pkin(4) - t63 * qJ(5) + t214) - g(2) * (pkin(4) * t123 + qJ(5) * t65 + t200) - t198; 0, 0, 0, -t179 * t184 * t181, t264 * t184, t252, t251, qJDD(2), -g(3) * t181 + t179 * t193, t181 * t193 + t296, -t52 * t101 + t118 * t85 + t37 * t131 + t7 * t178 - t57 * t72 + (-pkin(2) * t49 + t131 * t258 - t177 * t20 - t31 * t263) * t176 + t203, -t52 * t102 - t119 * t85 - t38 * t131 - t8 * t178 - t57 * t73 + (-pkin(2) * t50 + t131 * t257 + t175 * t20 + t32 * t263) * t176 - t204, t32 * t101 + t31 * t102 - t118 * t50 - t119 * t49 + t37 * t73 + t38 * t72 + (-t175 * t7 + t177 * t8 + (t175 * t73 - t177 * t72) * qJD(3) + t197) * t176, t8 * t119 + t7 * t118 - t32 * t38 - t31 * t37 - t52 * t57 - g(1) * (t150 - t249) - g(2) * (t148 - t250) - g(3) * t265 + (-t20 * pkin(2) + (-t175 * t31 + t177 * t32) * qJD(3)) * t176, -t24 * t101 - t23 * t102 - t34 * t73 + t97 * t49 + t98 * t50 + t288 * t72 + (-t177 * t5 + (qJD(3) * t73 + t6) * t175 + t197) * t176, t16 * t101 + t34 * t131 + t6 * t178 + t30 * t72 - t99 * t49 + t98 * t85 + (-t23 * t263 + t177 * t4 + (-qJD(3) * t131 + qJD(4) * t72) * t175) * t176 - t203, t16 * t102 - t5 * t178 + t30 * t73 - t99 * t50 - t97 * t85 + t288 * t131 + (t24 * t263 + (-t4 + t286) * t175) * t176 + t204, t4 * t99 + t5 * t97 + t6 * t98 - g(1) * (t226 - t249) - g(2) * (t227 - t250) - g(3) * t212 + t288 * t24 + (-t34 + t239) * t23 + (-t176 * t255 - t30) * t16, -t10 * t102 + t14 * t101 - t69 * t49 + t59 * t50 + t289 * t73 - t287 * t72 + (t175 * t2 + t177 * t3 + t197) * t176, t11 * t102 + t3 * t178 - t70 * t50 + t69 * t85 - t290 * t73 + (-t1 * t175 - t14 * t263) * t176 - t287 * t131 + t204, -t11 * t101 - t2 * t178 + t70 * t49 - t59 * t85 + t290 * t72 + (-t1 * t177 + t10 * t263) * t176 + t289 * t131 + t203, t1 * t70 + t2 * t59 + t3 * t69 - g(1) * (-t96 * qJ(5) + t182 * t213 + t226) - g(2) * (-t94 * qJ(5) + t180 * t213 + t227) - g(3) * (pkin(4) * t281 + qJ(5) * t117 + t212) + t287 * t14 + t290 * t11 + t289 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t187, t310, t221, t31 * t55 - t32 * t56 + t20 + t305, t221, t190 + t292, -t310, -t23 * t55 + t24 * t56 + t188, t221, -t310, t187, -t10 * t55 - t14 * t56 + t188 + t311; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, t191, t309, -t24 * t131 + t16 * t73 + t189 - t301, t185, t309, -t191, t11 * t73 + (qJD(5) + t14) * t131 + t189 + t307; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t131 * t73 + t190, t85 + t295, -t72 ^ 2 - t128, -g(1) * t65 - g(2) * t63 - g(3) * t115 - t10 * t131 - t11 * t72 + t3;];
tau_reg = t18;
