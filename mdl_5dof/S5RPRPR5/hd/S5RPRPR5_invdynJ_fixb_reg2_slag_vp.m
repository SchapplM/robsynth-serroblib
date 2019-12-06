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
% Datum: 2019-12-05 17:58
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:57:20
% EndTime: 2019-12-05 17:57:31
% DurationCPUTime: 4.47s
% Computational Cost: add. (5356->436), mult. (13301->575), div. (0->0), fcn. (9801->14), ass. (0->242)
t199 = sin(pkin(8));
t198 = sin(pkin(9));
t200 = cos(pkin(9));
t204 = sin(qJ(3));
t207 = cos(qJ(3));
t142 = t198 * t207 + t200 * t204;
t222 = qJD(1) * t142;
t106 = t199 * t222;
t278 = qJD(1) * t199;
t256 = t204 * t278;
t236 = t198 * t256;
t276 = qJD(1) * t207;
t255 = t199 * t276;
t110 = t200 * t255 - t236;
t203 = sin(qJ(5));
t206 = cos(qJ(5));
t270 = qJD(5) * t206;
t271 = qJD(5) * t203;
t263 = t199 * qJDD(1);
t247 = t207 * t263;
t248 = t204 * t263;
t267 = qJD(1) * qJD(3);
t250 = t207 * t267;
t335 = t199 * t250 + t248;
t65 = qJD(3) * t236 - t198 * t247 - t335 * t200;
t266 = qJDD(1) * t204;
t329 = t142 * qJD(3);
t66 = (qJD(1) * t329 + t198 * t266) * t199 - t200 * t247;
t16 = t106 * t270 + t110 * t271 - t203 * t65 + t206 * t66;
t201 = cos(pkin(8));
t277 = qJD(1) * t201;
t169 = -qJD(3) + t277;
t162 = -qJD(5) + t169;
t51 = t206 * t106 + t110 * t203;
t310 = t162 * t51;
t340 = -t16 - t310;
t315 = t51 ^ 2;
t230 = -t106 * t203 + t206 * t110;
t316 = t230 ^ 2;
t339 = -t315 + t316;
t314 = t51 * t230;
t338 = t201 * t222 - t329;
t229 = t198 * t204 - t200 * t207;
t132 = t229 * qJD(3);
t337 = t229 * t277 - t132;
t17 = qJD(5) * t230 - t203 * t66 - t206 * t65;
t311 = t162 * t230;
t336 = -t17 - t311;
t303 = qJ(2) * t207;
t168 = t201 * t303;
t228 = pkin(2) * t201 + pkin(6) * t199 + pkin(1);
t105 = -t204 * t228 + t168;
t334 = qJD(3) * t105;
t322 = pkin(7) * t110;
t129 = -qJD(1) * t228 + qJD(2);
t117 = t207 * t129;
t294 = t199 * t207;
t258 = qJ(4) * t294;
t304 = qJ(2) * t204;
t260 = t201 * t304;
t220 = -t258 - t260;
t74 = qJD(1) * t220 + t117;
t64 = -pkin(3) * t169 + t74;
t257 = qJ(2) * t277;
t89 = t129 * t204 + t207 * t257;
t75 = -qJ(4) * t256 + t89;
t69 = t198 * t75;
t35 = t200 * t64 - t69;
t20 = -pkin(4) * t169 - t322 + t35;
t323 = pkin(7) * t106;
t307 = t200 * t75;
t36 = t198 * t64 + t307;
t22 = t36 - t323;
t262 = t201 * qJDD(1);
t167 = -qJDD(3) + t262;
t275 = qJD(2) * t201;
t254 = t204 * t275;
t272 = qJD(4) * t199;
t218 = -t207 * t272 - t254;
t126 = -qJDD(1) * t228 + qJDD(2);
t116 = t207 * t126;
t274 = qJD(3) * t129;
t231 = -t204 * t274 + t116;
t296 = t199 * t204;
t259 = qJ(4) * t296;
t25 = -pkin(3) * t167 + t220 * qJDD(1) + ((-t168 + t259) * qJD(3) + t218) * qJD(1) + t231;
t214 = qJD(3) * t220 - t204 * t272;
t268 = qJD(1) * qJD(2);
t251 = t207 * t268;
t273 = qJD(3) * t207;
t239 = qJDD(1) * t168 + t204 * t126 + t129 * t273 + t201 * t251;
t31 = -qJ(4) * t248 + qJD(1) * t214 + t239;
t12 = -t198 * t31 + t200 * t25;
t6 = -pkin(4) * t167 + pkin(7) * t66 + t12;
t13 = t198 * t25 + t200 * t31;
t7 = pkin(7) * t65 + t13;
t1 = (qJD(5) * t20 + t7) * t206 + t203 * t6 - t22 * t271;
t195 = qJ(3) + pkin(9);
t188 = qJ(5) + t195;
t179 = sin(t188);
t208 = cos(qJ(1));
t180 = cos(t188);
t205 = sin(qJ(1));
t288 = t205 * t180;
t100 = -t179 * t208 + t201 * t288;
t291 = t201 * t208;
t102 = -t179 * t205 - t180 * t291;
t321 = g(1) * t199;
t130 = pkin(3) * t256 + qJ(2) * t278 + qJD(4);
t73 = pkin(4) * t106 + t130;
t333 = -g(2) * t100 - g(3) * t102 + t180 * t321 + t51 * t73 - t1;
t261 = 0.2e1 * qJ(2);
t331 = qJ(2) * qJDD(1) + t268;
t286 = t207 * t208;
t290 = t204 * t205;
t133 = t201 * t290 + t286;
t287 = t205 * t207;
t289 = t204 * t208;
t135 = t201 * t289 - t287;
t330 = -g(2) * t133 + g(3) * t135;
t101 = t179 * t291 - t288;
t9 = t20 * t203 + t206 * t22;
t2 = -qJD(5) * t9 - t203 * t7 + t206 * t6;
t292 = t201 * t205;
t99 = t179 * t292 + t180 * t208;
t328 = -g(2) * t99 + g(3) * t101 + t179 * t321 - t230 * t73 + t2;
t327 = t110 ^ 2;
t326 = pkin(4) * t65;
t325 = pkin(3) * t198;
t324 = pkin(3) * t204;
t319 = g(2) * t208;
t189 = t208 * qJ(2);
t317 = g(3) * t189;
t202 = -qJ(4) - pkin(6);
t81 = -t142 * t203 - t206 * t229;
t313 = qJD(5) * t81 + t338 * t203 + t337 * t206;
t82 = t142 * t206 - t203 * t229;
t312 = -qJD(5) * t82 - t337 * t203 + t338 * t206;
t283 = t207 * t275 - t228 * t273;
t62 = t214 + t283;
t63 = (-t168 + (qJ(4) * t199 + t228) * t204) * qJD(3) + t218;
t29 = t198 * t63 + t200 * t62;
t39 = t200 * t74 - t69;
t140 = t207 * t228;
t79 = -t258 - t140 + (-pkin(3) - t304) * t201;
t90 = -t259 + t105;
t44 = t198 * t79 + t200 * t90;
t88 = -t204 * t257 + t117;
t309 = t169 * t88;
t308 = t169 * t89;
t181 = pkin(3) * t200 + pkin(4);
t127 = t181 * t206 - t203 * t325;
t38 = -t198 * t74 - t307;
t26 = t38 + t323;
t27 = t39 - t322;
t306 = t127 * qJD(5) - t203 * t26 - t206 * t27;
t128 = t181 * t203 + t206 * t325;
t305 = -t128 * qJD(5) + t203 * t27 - t206 * t26;
t302 = qJDD(1) * pkin(1);
t301 = t106 * t169;
t300 = t110 * t106;
t299 = t110 * t169;
t298 = t167 * t201;
t193 = t199 ^ 2;
t209 = qJD(1) ^ 2;
t297 = t193 * t209;
t295 = t199 * t205;
t293 = t199 * t208;
t210 = qJ(2) ^ 2;
t242 = t268 * t261;
t264 = t193 * qJDD(1);
t282 = t193 * t242 + t210 * t264;
t281 = g(2) * t293 + g(3) * t295;
t186 = cos(t195);
t191 = t207 * pkin(3);
t147 = pkin(4) * t186 + t191;
t138 = (pkin(3) * t273 + qJD(2)) * t199;
t143 = pkin(3) * t296 + t199 * qJ(2);
t194 = t201 ^ 2;
t280 = t193 + t194;
t196 = t204 ^ 2;
t197 = t207 ^ 2;
t279 = t196 - t197;
t265 = qJDD(1) * t207;
t252 = t204 * t268;
t28 = -t198 * t62 + t200 * t63;
t43 = -t198 * t90 + t200 * t79;
t244 = t280 * t209;
t243 = qJD(1) * (-qJD(3) - t169);
t241 = t167 + t262;
t240 = t207 * t204 * t297;
t237 = qJD(3) * t260;
t234 = t204 * t250;
t233 = g(3) * t205 + t319;
t232 = -g(2) * t205 + g(3) * t208;
t125 = t229 * t199;
t34 = -pkin(4) * t201 + pkin(7) * t125 + t43;
t124 = t142 * t199;
t37 = -pkin(7) * t124 + t44;
t14 = -t203 * t37 + t206 * t34;
t15 = t203 * t34 + t206 * t37;
t68 = -t124 * t203 - t125 * t206;
t80 = t335 * pkin(3) + qJ(2) * t263 + t199 * t268 + qJDD(4);
t227 = qJD(3) * (t169 + t277);
t226 = (pkin(2) + t147) * t201 - (-pkin(7) + t202) * t199 + pkin(1);
t225 = (t191 + pkin(2)) * t201 - t199 * t202 + pkin(1);
t224 = t201 * t267 + t297;
t223 = t250 + t266;
t219 = -t233 - t302;
t217 = -t169 ^ 2 - t297;
t216 = g(1) * t201 + g(2) * t295 - g(3) * t293 + t80;
t185 = sin(t195);
t183 = qJDD(2) - t302;
t157 = -qJDD(5) + t167;
t146 = pkin(4) * t185 + t324;
t136 = -t201 * t286 - t290;
t134 = t201 * t287 - t289;
t123 = -t185 * t205 - t186 * t291;
t122 = t185 * t291 - t186 * t205;
t121 = -t185 * t208 + t186 * t292;
t120 = t185 * t292 + t186 * t208;
t113 = t199 * t329;
t109 = t199 * t132;
t104 = -t140 - t260;
t98 = t106 ^ 2;
t87 = -t254 - t334;
t86 = -t237 + t283;
t85 = pkin(3) * t255 + pkin(4) * t110;
t83 = pkin(4) * t124 + t143;
t76 = -pkin(4) * t109 + t138;
t67 = t206 * t124 - t125 * t203;
t46 = (-qJ(2) * t223 - t252) * t201 + t231;
t45 = -qJD(1) * t237 + t239;
t40 = t80 - t326;
t33 = qJD(5) * t68 - t206 * t109 - t113 * t203;
t32 = -t203 * t109 + t206 * t113 + t124 * t270 - t125 * t271;
t19 = pkin(7) * t109 + t29;
t18 = pkin(7) * t113 + t28;
t8 = t20 * t206 - t203 * t22;
t4 = -qJD(5) * t15 + t18 * t206 - t19 * t203;
t3 = qJD(5) * t14 + t18 * t203 + t19 * t206;
t5 = [0, 0, 0, 0, 0, qJDD(1), t233, t232, 0, 0, t264, 0.2e1 * t199 * t262, 0, t194 * qJDD(1), 0, 0, (-t183 - t219) * t201, (t183 - t302) * t199 - t281, 0.2e1 * t331 * t280 - t232, -t183 * pkin(1) - g(2) * (-pkin(1) * t208 - qJ(2) * t205) - g(3) * (-pkin(1) * t205 + t189) + (qJDD(1) * t210 + t242) * t194 + t282, (qJDD(1) * t197 - 0.2e1 * t234) * t193, 0.2e1 * (-t204 * t265 + t267 * t279) * t193, (t204 * t227 - t207 * t241) * t199, (qJDD(1) * t196 + 0.2e1 * t234) * t193, (t204 * t241 + t207 * t227) * t199, t298, -g(2) * t136 + g(3) * t134 - t104 * t167 - t169 * t87 - t201 * t46 + (t223 * t261 + 0.2e1 * t252) * t193, -g(2) * t135 - g(3) * t133 + t105 * t167 + t169 * t86 + t201 * t45 + (0.2e1 * t251 + (-t204 * t267 + t265) * t261) * t193, ((-qJD(3) * t89 - qJDD(1) * t104 - t46 + (-t87 - t334) * qJD(1)) * t207 + (qJD(3) * t88 - qJDD(1) * t105 - t45 + (qJD(3) * t104 - t86) * qJD(1)) * t204) * t199 + t281, -t317 + t46 * t104 + t45 * t105 + t89 * t86 + t88 * t87 + t228 * t319 + (g(2) * qJ(2) + g(3) * t228) * t205 + t282, -t110 * t113 + t125 * t66, t106 * t113 + t109 * t110 + t124 * t66 - t125 * t65, t113 * t169 + t125 * t167 + t201 * t66, -t106 * t109 - t124 * t65, -t109 * t169 + t124 * t167 - t201 * t65, t298, -g(2) * t123 + g(3) * t121 + t106 * t138 - t109 * t130 - t12 * t201 + t124 * t80 - t143 * t65 - t167 * t43 - t169 * t28, -g(2) * t122 - g(3) * t120 + t110 * t138 - t113 * t130 - t125 * t80 + t13 * t201 - t143 * t66 + t167 * t44 + t169 * t29, -t106 * t29 + t109 * t36 - t110 * t28 + t113 * t35 + t12 * t125 - t124 * t13 + t43 * t66 + t44 * t65 + t281, -t317 + t12 * t43 + t13 * t44 + t130 * t138 + t80 * t143 + t35 * t28 + t36 * t29 + (g(2) * t225 - g(3) * t324) * t208 + (-g(2) * (-qJ(2) - t324) + g(3) * t225) * t205, -t16 * t68 - t230 * t32, t16 * t67 - t17 * t68 - t230 * t33 + t32 * t51, -t157 * t68 + t16 * t201 + t162 * t32, t17 * t67 + t33 * t51, t157 * t67 + t162 * t33 + t17 * t201, t157 * t201, -g(2) * t102 + g(3) * t100 - t14 * t157 - t162 * t4 + t17 * t83 - t2 * t201 + t33 * t73 + t40 * t67 + t51 * t76, -g(2) * t101 - g(3) * t99 + t1 * t201 + t15 * t157 - t16 * t83 + t162 * t3 + t230 * t76 - t32 * t73 + t40 * t68, -t1 * t67 + t14 * t16 - t15 * t17 - t2 * t68 - t230 * t4 - t3 * t51 + t32 * t8 - t33 * t9 + t281, -t317 + t1 * t15 + t2 * t14 + t9 * t3 + t8 * t4 + t40 * t83 + t73 * t76 + (g(2) * t226 - g(3) * t146) * t208 + (-g(2) * (-qJ(2) - t146) + g(3) * t226) * t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t262, t263, -t244, -qJ(2) * t244 + qJDD(2) + t219, 0, 0, 0, 0, 0, 0, -t167 * t207 + t204 * t217, t167 * t204 + t207 * t217, (-t196 - t197) * t263, -qJ(2) * t297 + (t46 - t308) * t207 + (t45 + t309) * t204 - t233, 0, 0, 0, 0, 0, 0, -t106 * t278 + t167 * t229 - t169 * t338, -t110 * t278 + t142 * t167 + t169 * t337, -t106 * t337 - t110 * t338 + t142 * t65 - t229 * t66, -t12 * t229 + t13 * t142 - t130 * t278 + t337 * t36 + t338 * t35 - t233, 0, 0, 0, 0, 0, 0, -t157 * t81 - t312 * t162 - t51 * t278, t157 * t82 + t313 * t162 - t230 * t278, t16 * t81 - t17 * t82 - t230 * t312 - t313 * t51, t1 * t82 + t2 * t81 - t278 * t73 + t312 * t8 + t313 * t9 - t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t240, -t279 * t297, (t204 * t243 + t265) * t199, -t240, (t207 * t243 - t266) * t199, -t167, -t308 + t116 - t224 * t303 + (-t331 * t201 - t274 + t321) * t204 + t330, g(1) * t294 - g(2) * t134 - g(3) * t136 + t224 * t304 - t239 - t309, 0, 0, t300, -t98 + t327, -t66 - t301, -t300, t65 - t299, -t167, t185 * t321 - g(2) * t120 + g(3) * t122 - t110 * t130 + t169 * t38 + (-t106 * t255 - t167 * t200) * pkin(3) + t12, t186 * t321 - g(2) * t121 - g(3) * t123 + t106 * t130 - t169 * t39 + (-t110 * t255 + t167 * t198) * pkin(3) - t13, (t36 + t38) * t110 + (-t35 + t39) * t106 + (t198 * t65 + t200 * t66) * pkin(3), -t35 * t38 - t36 * t39 + (t13 * t198 + t12 * t200 + (g(1) * t204 - t130 * t276) * t199 + t330) * pkin(3), t314, t339, t340, -t314, t336, -t157, -t127 * t157 - t305 * t162 - t51 * t85 + t328, t128 * t157 + t306 * t162 - t230 * t85 + t333, t127 * t16 - t128 * t17 + (-t305 + t9) * t230 + (-t306 - t8) * t51, t1 * t128 + t2 * t127 - t73 * t85 + t146 * t321 - g(2) * (t146 * t292 + t147 * t208) - g(3) * (-t146 * t291 + t147 * t205) + t306 * t9 + t305 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65 - t299, -t66 + t301, -t98 - t327, t106 * t36 + t110 * t35 + t216, 0, 0, 0, 0, 0, 0, t17 - t311, -t16 + t310, -t315 - t316, t230 * t8 + t51 * t9 + t216 - t326; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t314, t339, t340, -t314, t336, -t157, -t162 * t9 + t328, -t162 * t8 + t333, 0, 0;];
tau_reg = t5;
