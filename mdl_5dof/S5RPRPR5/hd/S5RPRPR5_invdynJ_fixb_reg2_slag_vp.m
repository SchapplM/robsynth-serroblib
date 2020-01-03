% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRPR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:43:13
% EndTime: 2020-01-03 11:43:22
% DurationCPUTime: 4.66s
% Computational Cost: add. (5356->435), mult. (13301->576), div. (0->0), fcn. (9801->14), ass. (0->247)
t201 = sin(pkin(8));
t200 = sin(pkin(9));
t202 = cos(pkin(9));
t206 = sin(qJ(3));
t209 = cos(qJ(3));
t142 = t200 * t209 + t202 * t206;
t224 = qJD(1) * t142;
t106 = t201 * t224;
t282 = qJD(1) * t201;
t260 = t206 * t282;
t237 = t200 * t260;
t280 = qJD(1) * t209;
t259 = t201 * t280;
t110 = t202 * t259 - t237;
t205 = sin(qJ(5));
t208 = cos(qJ(5));
t274 = qJD(5) * t208;
t275 = qJD(5) * t205;
t267 = t201 * qJDD(1);
t249 = t209 * t267;
t250 = t206 * t267;
t271 = qJD(1) * qJD(3);
t252 = t209 * t271;
t340 = t201 * t252 + t250;
t65 = qJD(3) * t237 - t200 * t249 - t340 * t202;
t270 = qJDD(1) * t206;
t332 = t142 * qJD(3);
t66 = (qJD(1) * t332 + t200 * t270) * t201 - t202 * t249;
t16 = t106 * t274 + t110 * t275 - t205 * t65 + t208 * t66;
t203 = cos(pkin(8));
t281 = qJD(1) * t203;
t169 = -qJD(3) + t281;
t162 = -qJD(5) + t169;
t51 = t208 * t106 + t110 * t205;
t314 = t162 * t51;
t345 = -t16 - t314;
t319 = t51 ^ 2;
t230 = -t106 * t205 + t208 * t110;
t320 = t230 ^ 2;
t344 = -t319 + t320;
t318 = t51 * t230;
t343 = t203 * t224 - t332;
t229 = t200 * t206 - t202 * t209;
t132 = t229 * qJD(3);
t342 = t229 * t281 - t132;
t17 = qJD(5) * t230 - t205 * t66 - t208 * t65;
t315 = t162 * t230;
t341 = -t17 - t315;
t306 = qJDD(1) * pkin(1);
t181 = qJDD(2) - t306;
t207 = sin(qJ(1));
t210 = cos(qJ(1));
t335 = -g(2) * t210 - g(3) * t207;
t339 = t335 - t181;
t234 = pkin(2) * t203 + pkin(6) * t201;
t151 = -pkin(1) - t234;
t307 = qJ(2) * t209;
t168 = t203 * t307;
t105 = t206 * t151 + t168;
t338 = qJD(3) * t105;
t325 = pkin(7) * t110;
t129 = qJD(1) * t151 + qJD(2);
t117 = t209 * t129;
t298 = t201 * t209;
t262 = qJ(4) * t298;
t308 = qJ(2) * t206;
t264 = t203 * t308;
t222 = -t262 - t264;
t74 = qJD(1) * t222 + t117;
t64 = -pkin(3) * t169 + t74;
t261 = qJ(2) * t281;
t89 = t129 * t206 + t209 * t261;
t75 = -qJ(4) * t260 + t89;
t69 = t200 * t75;
t35 = t202 * t64 - t69;
t20 = -pkin(4) * t169 - t325 + t35;
t326 = pkin(7) * t106;
t311 = t202 * t75;
t36 = t200 * t64 + t311;
t22 = t36 - t326;
t266 = t203 * qJDD(1);
t167 = -qJDD(3) + t266;
t279 = qJD(2) * t203;
t256 = t206 * t279;
t276 = qJD(4) * t201;
t221 = -t209 * t276 - t256;
t126 = qJDD(1) * t151 + qJDD(2);
t116 = t209 * t126;
t278 = qJD(3) * t129;
t231 = -t206 * t278 + t116;
t300 = t201 * t206;
t263 = qJ(4) * t300;
t25 = -pkin(3) * t167 + t222 * qJDD(1) + ((-t168 + t263) * qJD(3) + t221) * qJD(1) + t231;
t216 = qJD(3) * t222 - t206 * t276;
t272 = qJD(1) * qJD(2);
t253 = t209 * t272;
t277 = qJD(3) * t209;
t240 = qJDD(1) * t168 + t206 * t126 + t129 * t277 + t203 * t253;
t31 = -qJ(4) * t250 + qJD(1) * t216 + t240;
t12 = -t200 * t31 + t202 * t25;
t6 = -pkin(4) * t167 + pkin(7) * t66 + t12;
t13 = t200 * t25 + t202 * t31;
t7 = pkin(7) * t65 + t13;
t1 = (qJD(5) * t20 + t7) * t208 + t205 * t6 - t22 * t275;
t197 = qJ(3) + pkin(9);
t186 = qJ(5) + t197;
t178 = cos(t186);
t177 = sin(t186);
t291 = t210 * t177;
t297 = t203 * t207;
t100 = t178 * t297 - t291;
t296 = t203 * t210;
t102 = t177 * t207 + t178 * t296;
t324 = g(1) * t201;
t130 = pkin(3) * t260 + qJ(2) * t282 + qJD(4);
t73 = pkin(4) * t106 + t130;
t337 = g(2) * t100 - g(3) * t102 + t178 * t324 + t51 * t73 - t1;
t265 = 0.2e1 * qJ(2);
t334 = qJ(2) * qJDD(1) + t272;
t292 = t209 * t210;
t295 = t206 * t207;
t133 = -t203 * t295 - t292;
t293 = t207 * t209;
t294 = t206 * t210;
t135 = t203 * t294 - t293;
t333 = -g(2) * t133 - g(3) * t135;
t101 = -t178 * t207 + t203 * t291;
t9 = t20 * t205 + t208 * t22;
t2 = -qJD(5) * t9 - t205 * t7 + t208 * t6;
t99 = -t177 * t297 - t178 * t210;
t331 = -g(2) * t99 - g(3) * t101 + t177 * t324 - t230 * t73 + t2;
t330 = t110 ^ 2;
t329 = pkin(4) * t65;
t328 = pkin(3) * t200;
t327 = pkin(3) * t206;
t321 = g(3) * t210;
t204 = -qJ(4) - pkin(6);
t81 = -t142 * t205 - t208 * t229;
t317 = qJD(5) * t81 + t343 * t205 + t342 * t208;
t82 = t142 * t208 - t205 * t229;
t316 = -qJD(5) * t82 - t342 * t205 + t343 * t208;
t288 = t151 * t277 + t209 * t279;
t62 = t216 + t288;
t63 = (-t168 + (qJ(4) * t201 - t151) * t206) * qJD(3) + t221;
t29 = t200 * t63 + t202 * t62;
t39 = t202 * t74 - t69;
t140 = t209 * t151;
t79 = -t262 + t140 + (-pkin(3) - t308) * t203;
t90 = -t263 + t105;
t44 = t200 * t79 + t202 * t90;
t88 = -t206 * t261 + t117;
t313 = t169 * t88;
t312 = t169 * t89;
t179 = pkin(3) * t202 + pkin(4);
t127 = t179 * t208 - t205 * t328;
t38 = -t200 * t74 - t311;
t26 = t38 + t326;
t27 = t39 - t325;
t310 = t127 * qJD(5) - t205 * t26 - t208 * t27;
t128 = t179 * t205 + t208 * t328;
t309 = -t128 * qJD(5) + t205 * t27 - t208 * t26;
t305 = t106 * t169;
t304 = t110 * t106;
t303 = t110 * t169;
t302 = t167 * t203;
t195 = t201 ^ 2;
t211 = qJD(1) ^ 2;
t301 = t195 * t211;
t299 = t201 * t207;
t212 = qJ(2) ^ 2;
t243 = t272 * t265;
t268 = t195 * qJDD(1);
t287 = t195 * t243 + t212 * t268;
t184 = cos(t197);
t190 = t209 * pkin(3);
t147 = pkin(4) * t184 + t190;
t138 = (pkin(3) * t277 + qJD(2)) * t201;
t143 = pkin(3) * t300 + t201 * qJ(2);
t286 = t210 * pkin(1) + t207 * qJ(2);
t196 = t203 ^ 2;
t284 = t195 + t196;
t198 = t206 ^ 2;
t199 = t209 ^ 2;
t283 = t198 - t199;
t269 = qJDD(1) * t209;
t258 = t106 * t282;
t257 = t110 * t282;
t254 = t206 * t272;
t28 = -t200 * t62 + t202 * t63;
t43 = -t200 * t90 + t202 * t79;
t246 = t284 * t211;
t189 = t207 * pkin(1);
t245 = -t210 * qJ(2) + t189;
t244 = qJD(1) * (-qJD(3) - t169);
t242 = t167 + t266;
t241 = t209 * t206 * t301;
t238 = qJD(3) * t264;
t235 = t206 * t252;
t232 = g(2) * t207 - t321;
t125 = t229 * t201;
t34 = -pkin(4) * t203 + pkin(7) * t125 + t43;
t124 = t142 * t201;
t37 = -pkin(7) * t124 + t44;
t14 = -t205 * t37 + t208 * t34;
t15 = t205 * t34 + t208 * t37;
t68 = -t124 * t205 - t125 * t208;
t80 = t340 * pkin(3) + qJ(2) * t267 + t201 * t272 + qJDD(4);
t228 = qJD(3) * (t169 + t281);
t227 = t335 * t201;
t226 = t203 * t271 + t301;
t225 = t252 + t270;
t220 = t306 + t339;
t219 = -t169 ^ 2 - t301;
t218 = g(1) * t203 - g(2) * t299 + t201 * t321 + t80;
t194 = -pkin(7) + t204;
t183 = sin(t197);
t180 = t190 + pkin(2);
t157 = -qJDD(5) + t167;
t146 = pkin(4) * t183 + t327;
t144 = pkin(2) + t147;
t136 = t203 * t292 + t295;
t134 = t203 * t293 - t294;
t123 = t183 * t207 + t184 * t296;
t122 = t183 * t296 - t184 * t207;
t121 = -t183 * t210 + t184 * t297;
t120 = -t183 * t297 - t184 * t210;
t113 = t201 * t332;
t109 = t201 * t132;
t104 = t140 - t264;
t98 = t106 ^ 2;
t87 = -t256 - t338;
t86 = -t238 + t288;
t85 = pkin(3) * t259 + pkin(4) * t110;
t83 = pkin(4) * t124 + t143;
t76 = -pkin(4) * t109 + t138;
t67 = t208 * t124 - t125 * t205;
t46 = (-qJ(2) * t225 - t254) * t203 + t231;
t45 = -qJD(1) * t238 + t240;
t40 = t80 - t329;
t33 = qJD(5) * t68 - t208 * t109 - t205 * t113;
t32 = -t205 * t109 + t208 * t113 + t124 * t274 - t125 * t275;
t19 = pkin(7) * t109 + t29;
t18 = pkin(7) * t113 + t28;
t8 = t20 * t208 - t205 * t22;
t4 = -qJD(5) * t15 + t208 * t18 - t205 * t19;
t3 = qJD(5) * t14 + t205 * t18 + t208 * t19;
t5 = [0, 0, 0, 0, 0, qJDD(1), t335, t232, 0, 0, t268, 0.2e1 * t201 * t266, 0, t196 * qJDD(1), 0, 0, t220 * t203, -t220 * t201, 0.2e1 * t334 * t284 - t232, -t181 * pkin(1) - g(2) * t286 - g(3) * t245 + (qJDD(1) * t212 + t243) * t196 + t287, (qJDD(1) * t199 - 0.2e1 * t235) * t195, 0.2e1 * (-t206 * t269 + t271 * t283) * t195, (t206 * t228 - t209 * t242) * t201, (qJDD(1) * t198 + 0.2e1 * t235) * t195, (t206 * t242 + t209 * t228) * t201, t302, -g(2) * t136 - g(3) * t134 - t104 * t167 - t169 * t87 - t203 * t46 + (t225 * t265 + 0.2e1 * t254) * t195, g(2) * t135 - g(3) * t133 + t105 * t167 + t169 * t86 + t203 * t45 + (0.2e1 * t253 + (-t206 * t271 + t269) * t265) * t195, ((-qJD(3) * t89 - qJDD(1) * t104 - t46 + (-t87 - t338) * qJD(1)) * t209 + (qJD(3) * t88 - qJDD(1) * t105 - t45 + (qJD(3) * t104 - t86) * qJD(1)) * t206 + t335) * t201, t45 * t105 + t89 * t86 + t46 * t104 + t88 * t87 - g(2) * (t210 * t234 + t286) - g(3) * (t207 * t234 + t245) + t287, -t110 * t113 + t125 * t66, t106 * t113 + t109 * t110 + t124 * t66 - t125 * t65, t113 * t169 + t125 * t167 + t203 * t66, -t106 * t109 - t124 * t65, -t109 * t169 + t124 * t167 - t203 * t65, t302, -g(2) * t123 - g(3) * t121 + t106 * t138 - t109 * t130 - t12 * t203 + t124 * t80 - t143 * t65 - t167 * t43 - t169 * t28, g(2) * t122 - g(3) * t120 + t110 * t138 - t113 * t130 - t125 * t80 + t13 * t203 - t143 * t66 + t167 * t44 + t169 * t29, -t106 * t29 + t109 * t36 - t110 * t28 + t113 * t35 + t12 * t125 - t124 * t13 + t43 * t66 + t44 * t65 + t227, t13 * t44 + t36 * t29 + t12 * t43 + t35 * t28 + t80 * t143 + t130 * t138 - g(2) * (pkin(3) * t295 + t286) - g(3) * (t180 * t297 - t204 * t299 + t189) + (-g(2) * (t180 * t203 - t201 * t204) - g(3) * (-qJ(2) - t327)) * t210, -t16 * t68 - t230 * t32, t16 * t67 - t17 * t68 - t230 * t33 + t32 * t51, -t157 * t68 + t16 * t203 + t162 * t32, t17 * t67 + t33 * t51, t157 * t67 + t162 * t33 + t17 * t203, t157 * t203, -g(2) * t102 - g(3) * t100 - t14 * t157 - t162 * t4 + t17 * t83 - t2 * t203 + t33 * t73 + t40 * t67 + t51 * t76, g(2) * t101 - g(3) * t99 + t1 * t203 + t15 * t157 - t16 * t83 + t162 * t3 + t230 * t76 - t32 * t73 + t40 * t68, -t1 * t67 + t14 * t16 - t15 * t17 - t2 * t68 - t230 * t4 - t3 * t51 + t32 * t8 - t33 * t9 + t227, t1 * t15 + t9 * t3 + t2 * t14 + t8 * t4 + t40 * t83 + t73 * t76 - g(2) * (t146 * t207 + t286) - g(3) * (t144 * t297 - t194 * t299 + t189) + (-g(2) * (t144 * t203 - t194 * t201) - g(3) * (-qJ(2) - t146)) * t210; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t266, t267, -t246, -qJ(2) * t246 - t339, 0, 0, 0, 0, 0, 0, -t167 * t209 + t206 * t219, t167 * t206 + t209 * t219, (-t198 - t199) * t267, -qJ(2) * t301 + (t46 - t312) * t209 + (t45 + t313) * t206 - t335, 0, 0, 0, 0, 0, 0, t167 * t229 - t169 * t343 - t258, t142 * t167 + t169 * t342 - t257, -t106 * t342 - t110 * t343 + t142 * t65 - t229 * t66, -t12 * t229 + t13 * t142 - t130 * t282 + t342 * t36 + t343 * t35 - t335, 0, 0, 0, 0, 0, 0, -t157 * t81 - t162 * t316 - t51 * t282, t157 * t82 + t317 * t162 - t230 * t282, t16 * t81 - t17 * t82 - t230 * t316 - t317 * t51, t1 * t82 + t2 * t81 - t73 * t282 + t316 * t8 + t317 * t9 - t335; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, -t283 * t301, (t206 * t244 + t269) * t201, -t241, (t209 * t244 - t270) * t201, -t167, -t312 + t116 - t226 * t307 + (-t334 * t203 - t278 + t324) * t206 + t333, g(1) * t298 + g(2) * t134 - g(3) * t136 + t226 * t308 - t240 - t313, 0, 0, t304, -t98 + t330, -t66 - t305, -t304, t65 - t303, -t167, t183 * t324 - g(2) * t120 - g(3) * t122 - t110 * t130 + t169 * t38 + (-t167 * t202 - t209 * t258) * pkin(3) + t12, t184 * t324 + g(2) * t121 - g(3) * t123 + t106 * t130 - t169 * t39 + (t167 * t200 - t209 * t257) * pkin(3) - t13, (t36 + t38) * t110 + (-t35 + t39) * t106 + (t200 * t65 + t202 * t66) * pkin(3), -t35 * t38 - t36 * t39 + (t13 * t200 + t12 * t202 + (g(1) * t206 - t130 * t280) * t201 + t333) * pkin(3), t318, t344, t345, -t318, t341, -t157, -t127 * t157 - t162 * t309 - t51 * t85 + t331, t128 * t157 + t162 * t310 - t230 * t85 + t337, t127 * t16 - t128 * t17 + (-t309 + t9) * t230 + (-t310 - t8) * t51, t1 * t128 + t2 * t127 - t73 * t85 + t146 * t324 - g(2) * (-t146 * t297 - t147 * t210) - g(3) * (t146 * t296 - t147 * t207) + t310 * t9 + t309 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65 - t303, -t66 + t305, -t98 - t330, t106 * t36 + t110 * t35 + t218, 0, 0, 0, 0, 0, 0, t17 - t315, -t16 + t314, -t319 - t320, t230 * t8 + t51 * t9 + t218 - t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t318, t344, t345, -t318, t341, -t157, -t162 * t9 + t331, -t162 * t8 + t337, 0, 0;];
tau_reg = t5;
