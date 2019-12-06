% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRRRP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:48
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP2_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP2_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:48:19
% EndTime: 2019-12-05 18:48:26
% DurationCPUTime: 3.01s
% Computational Cost: add. (4719->373), mult. (7259->430), div. (0->0), fcn. (4699->12), ass. (0->222)
t185 = qJ(3) + qJ(4);
t172 = sin(t185);
t174 = cos(t185);
t186 = qJ(1) + qJ(2);
t173 = sin(t186);
t159 = g(2) * t173;
t175 = cos(t186);
t311 = -g(3) * t175 + t159;
t309 = -g(1) * t174 - t172 * t311;
t187 = sin(qJ(4));
t191 = cos(qJ(3));
t188 = sin(qJ(3));
t307 = cos(qJ(4));
t249 = t307 * t188;
t120 = t187 * t191 + t249;
t182 = qJD(1) + qJD(2);
t98 = t120 * t182;
t285 = t98 * qJ(5);
t189 = sin(qJ(2));
t290 = pkin(1) * qJD(1);
t255 = t189 * t290;
t128 = t182 * pkin(7) + t255;
t242 = pkin(8) * t182 + t128;
t87 = t242 * t191;
t76 = t187 * t87;
t86 = t242 * t188;
t79 = qJD(3) * pkin(3) - t86;
t44 = t307 * t79 - t76;
t28 = t44 - t285;
t192 = cos(qJ(2));
t183 = t188 ^ 2;
t184 = t191 ^ 2;
t268 = t183 + t184;
t313 = t268 * t192;
t316 = t182 * t313;
t179 = qJDD(1) + qJDD(2);
t261 = qJDD(1) * t189;
t265 = qJD(2) * t192;
t102 = t179 * pkin(7) + (qJD(1) * t265 + t261) * pkin(1);
t238 = t268 * t102;
t222 = g(2) * t175 + g(3) * t173;
t244 = qJD(4) * t307;
t248 = t307 * t191;
t315 = -qJD(3) * t248 - t191 * t244;
t267 = qJD(1) * t192;
t254 = pkin(1) * t267;
t302 = t182 * pkin(2);
t129 = -t254 - t302;
t314 = t128 * t313 + t129 * t189;
t194 = -pkin(8) - pkin(7);
t250 = qJD(3) * t194;
t125 = t188 * t250;
t126 = t191 * t250;
t144 = t194 * t188;
t176 = t191 * pkin(8);
t145 = t191 * pkin(7) + t176;
t276 = t187 * t188;
t211 = t248 - t276;
t262 = qJD(4) * t187;
t292 = t307 * t125 + t187 * t126 + t144 * t244 - t145 * t262 - t211 * t254;
t81 = t187 * t144 + t307 * t145;
t291 = -qJD(4) * t81 + t120 * t254 - t187 * t125 + t307 * t126;
t164 = t189 * pkin(1) + pkin(7);
t297 = -pkin(8) - t164;
t116 = t297 * t188;
t117 = t191 * t164 + t176;
t64 = t187 * t116 + t307 * t117;
t178 = qJDD(3) + qJDD(4);
t181 = qJD(3) + qJD(4);
t229 = pkin(3) * t244;
t306 = pkin(3) * t187;
t312 = -t178 * t306 - t181 * t229;
t171 = t178 * pkin(4);
t218 = t181 * t276;
t274 = t191 * t179;
t228 = -t179 * t249 + t315 * t182 - t187 * t274;
t35 = t182 * t218 + t228;
t289 = t35 * qJ(5);
t310 = t171 + t289;
t308 = t98 ^ 2;
t304 = g(1) * t191;
t303 = t179 * pkin(2);
t190 = sin(qJ(1));
t301 = t190 * pkin(1);
t300 = t192 * pkin(1);
t193 = cos(qJ(1));
t299 = t193 * pkin(1);
t253 = t182 * t276;
t96 = -t182 * t248 + t253;
t298 = t98 * t96;
t68 = t181 * t120;
t282 = -t68 * qJ(5) + qJD(5) * t211;
t296 = t282 + t292;
t67 = t218 + t315;
t215 = t67 * qJ(5) - t120 * qJD(5);
t295 = t215 + t291;
t27 = t181 * pkin(4) + t28;
t294 = t27 - t28;
t275 = t188 * t179;
t219 = -t179 * t248 + t187 * t275;
t36 = t182 * t68 + t219;
t293 = -t96 * t229 - t36 * t306;
t48 = -t307 * t86 - t76;
t288 = t36 * qJ(5);
t287 = t96 * qJ(5);
t286 = t96 * t181;
t270 = -qJD(2) * t255 + qJDD(1) * t300;
t101 = -t270 - t303;
t263 = qJD(3) * t191;
t283 = t101 * t188 + t129 * t263;
t281 = t120 * qJ(5);
t279 = t173 * t174;
t278 = t174 * t175;
t277 = t182 * t188;
t166 = t191 * pkin(3) + pkin(2);
t100 = -t166 * t182 - t254;
t243 = t96 * pkin(4) + qJD(5);
t58 = t100 + t243;
t273 = qJD(5) + t58;
t272 = g(2) * t278 + g(3) * t279;
t269 = t183 - t184;
t266 = qJD(2) * t189;
t264 = qJD(3) * t188;
t259 = pkin(3) * t277;
t258 = t307 * pkin(3);
t257 = pkin(1) * t265;
t256 = pkin(3) * t262;
t169 = pkin(3) * t264;
t78 = t307 * t87;
t177 = t182 ^ 2;
t252 = t188 * t177 * t191;
t251 = t129 * t264 + t222 * t191;
t247 = t182 * t266;
t246 = t182 * t264;
t245 = t182 * t263;
t61 = t68 * pkin(4) + t169;
t41 = -t128 * t263 + qJDD(3) * pkin(3) - t188 * t102 + (-t245 - t275) * pkin(8);
t46 = -t128 * t264 + t191 * t102 + (-t246 + t274) * pkin(8);
t241 = -t187 * t46 + t307 * t41;
t47 = t187 * t86 - t78;
t240 = qJD(3) * t297;
t8 = t187 * t41 + t79 * t244 - t87 * t262 + t307 * t46;
t237 = t268 * t179;
t63 = t307 * t116 - t187 * t117;
t80 = t307 * t144 - t187 * t145;
t162 = pkin(4) * t174;
t124 = t162 + t166;
t180 = -qJ(5) + t194;
t235 = -t175 * t124 + t173 * t180;
t234 = -t175 * t166 + t173 * t194;
t62 = pkin(3) * t246 - t166 * t179 - t270;
t20 = t36 * pkin(4) + qJDD(5) + t62;
t233 = -t20 * t211 + t58 * t68 + t272;
t232 = t100 * t68 - t211 * t62 + t272;
t231 = t311 + t238;
t230 = t182 * t255;
t227 = -t270 - t222;
t226 = t188 * t245;
t223 = t61 - t255;
t220 = g(2) * t193 + g(3) * t190;
t45 = t187 * t79 + t78;
t29 = t45 - t287;
t9 = -t45 * qJD(4) + t241;
t3 = -t98 * qJD(5) + t310 + t9;
t4 = -t96 * qJD(5) - t288 + t8;
t217 = -t3 * t120 + t211 * t4 + t27 * t67 - t29 * t68 + t311;
t216 = -t9 * t120 + t211 * t8 + t44 * t67 - t45 * t68 + t311;
t214 = -t173 * t124 - t175 * t180;
t213 = -t173 * t166 - t175 * t194;
t95 = -pkin(4) * t211 - t166;
t212 = t222 * t172;
t83 = t188 * t240 + t191 * t257;
t84 = -t188 * t257 + t191 * t240;
t21 = t116 * t244 - t117 * t262 + t187 * t84 + t307 * t83;
t210 = -t255 + t169;
t209 = -t129 * t182 - t102 - t311;
t208 = t20 * t120 - t58 * t67 - t212;
t207 = -t100 * t67 + t62 * t120 - t212;
t206 = g(1) * t172 - g(2) * t279 + g(3) * t278 - t8;
t195 = qJD(3) ^ 2;
t205 = -pkin(7) * t195 + t230 + t303;
t167 = -pkin(2) - t300;
t204 = -pkin(1) * t247 - t164 * t195 - t167 * t179;
t203 = -t181 * t253 - t228;
t202 = -pkin(7) * qJDD(3) + (t254 - t302) * qJD(3);
t201 = -qJDD(3) * t164 + (t167 * t182 - t257) * qJD(3);
t200 = t100 * t96 + t206;
t22 = -qJD(4) * t64 - t187 * t83 + t307 * t84;
t199 = t273 * t96 + t206 + t288;
t198 = t9 + t309;
t197 = -t100 * t98 + t198;
t170 = pkin(1) * t266;
t165 = t258 + pkin(4);
t156 = t175 * pkin(7);
t139 = -t166 - t300;
t138 = t188 * qJDD(3) + t195 * t191;
t137 = qJDD(3) * t191 - t195 * t188;
t127 = t170 + t169;
t115 = t211 * qJ(5);
t104 = t184 * t179 - 0.2e1 * t226;
t103 = t183 * t179 + 0.2e1 * t226;
t94 = t96 ^ 2;
t85 = t95 - t300;
t70 = -0.2e1 * t269 * t182 * qJD(3) + 0.2e1 * t188 * t274;
t69 = t98 * pkin(4) + t259;
t60 = t115 + t81;
t59 = t80 - t281;
t57 = t170 + t61;
t52 = t115 + t64;
t51 = t63 - t281;
t50 = t178 * t211 - t68 * t181;
t49 = t120 * t178 - t67 * t181;
t39 = -t94 + t308;
t31 = -t285 + t48;
t30 = t47 + t287;
t23 = t203 + t286;
t13 = -t211 * t36 + t96 * t68;
t12 = -t35 * t120 - t98 * t67;
t11 = t215 + t22;
t10 = t21 + t282;
t5 = -t120 * t36 - t211 * t35 + t67 * t96 - t98 * t68;
t1 = [0, 0, 0, 0, 0, qJDD(1), t220, -g(2) * t190 + g(3) * t193, 0, 0, 0, 0, 0, 0, 0, t179, (t179 * t192 - t247) * pkin(1) - t227, ((-qJDD(1) - t179) * t189 + (-qJD(1) - t182) * t265) * pkin(1) - t311, 0, (t220 + (t189 ^ 2 + t192 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t103, t70, t138, t104, t137, 0, t201 * t188 + (-t101 + t204) * t191 + t251, t201 * t191 + (-t204 - t222) * t188 + t283, pkin(1) * qJD(2) * t316 + t164 * t237 + t231, t101 * t167 - g(2) * (-t175 * pkin(2) - t173 * pkin(7)) - g(3) * (-t173 * pkin(2) + t156) + t164 * t238 + (t314 * qJD(2) + t220) * pkin(1), t12, t5, t49, t13, t50, 0, t127 * t96 + t139 * t36 + t63 * t178 + t22 * t181 + t232, t127 * t98 - t139 * t35 - t64 * t178 - t21 * t181 + t207, -t21 * t96 - t22 * t98 + t63 * t35 - t64 * t36 + t216, t8 * t64 + t45 * t21 + t9 * t63 + t44 * t22 + t62 * t139 + t100 * t127 - g(2) * (t234 - t299) - g(3) * (t213 - t301), t12, t5, t49, t13, t50, 0, t11 * t181 + t51 * t178 + t85 * t36 + t57 * t96 + t233, -t10 * t181 - t52 * t178 - t85 * t35 + t57 * t98 + t208, -t10 * t96 - t11 * t98 + t51 * t35 - t52 * t36 + t217, t4 * t52 + t29 * t10 + t3 * t51 + t27 * t11 + t20 * t85 + t58 * t57 - g(2) * (t235 - t299) - g(3) * (t214 - t301); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, -t227 + t230, (-t261 + (-qJD(2) + t182) * t267) * pkin(1) - t311, 0, 0, t103, t70, t138, t104, t137, 0, t202 * t188 + (-t101 + t205) * t191 + t251, t202 * t191 + (-t205 - t222) * t188 + t283, pkin(7) * t237 - t290 * t316 + t231, -g(3) * t156 + (-t101 + t222) * pkin(2) + (t238 + t159) * pkin(7) - t314 * t290, t12, t5, t49, t13, t50, 0, -t166 * t36 + t80 * t178 + t181 * t291 + t210 * t96 + t232, t166 * t35 - t81 * t178 - t181 * t292 + t210 * t98 + t207, -t291 * t98 - t292 * t96 + t80 * t35 - t81 * t36 + t216, -g(2) * t234 - g(3) * t213 + t100 * t210 - t62 * t166 + t291 * t44 + t292 * t45 + t8 * t81 + t9 * t80, t12, t5, t49, t13, t50, 0, t59 * t178 + t181 * t295 + t223 * t96 + t95 * t36 + t233, -t60 * t178 - t181 * t296 + t223 * t98 - t95 * t35 + t208, -t295 * t98 - t296 * t96 + t59 * t35 - t60 * t36 + t217, -g(2) * t235 - g(3) * t214 + t20 * t95 + t223 * t58 + t27 * t295 + t29 * t296 + t3 * t59 + t4 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t252, t269 * t177, t275, t252, t274, qJDD(3), t188 * t209 - t304, g(1) * t188 + t191 * t209, 0, 0, t298, t39, t23, -t298, -t219, t178, -t47 * t181 + (t178 * t307 - t181 * t262 - t277 * t96) * pkin(3) + t197, t48 * t181 - t259 * t98 + t200 + t312, t35 * t258 + (-t44 + t48) * t96 + (t45 + t47 + t256) * t98 + t293, -t44 * t47 - t45 * t48 + (t307 * t9 - t304 + t187 * t8 + (-t187 * t44 + t307 * t45) * qJD(4) + (-t100 * t182 - t311) * t188) * pkin(3), t298, t39, t23, -t298, -t219, t178, t165 * t178 - t30 * t181 - t69 * t96 - t273 * t98 + (-t78 + (-pkin(3) * t181 - t79) * t187) * qJD(4) + t241 + t309 + t310, t31 * t181 - t69 * t98 + t199 + t312, t165 * t35 + (-t27 + t31) * t96 + (t29 + t30 + t256) * t98 + t293, -g(1) * t162 + t3 * t165 - t27 * t30 - t29 * t31 - t58 * t69 + t311 * (-t188 * pkin(3) - pkin(4) * t172) + (-t304 + t4 * t187 + (-t187 * t27 + t29 * t307) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t298, t39, t23, -t298, -t219, t178, t45 * t181 + t197, t44 * t181 + t200, 0, 0, t298, t39, t23, -t298, -t219, t178, t289 + t29 * t181 + 0.2e1 * t171 + (-t243 - t58) * t98 + t198, -t308 * pkin(4) + t28 * t181 + t199, t35 * pkin(4) - t294 * t96, t294 * t29 + (-t58 * t98 + t3 + t309) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98 * t181 + t36, t203 - t286, -t94 - t308, t27 * t98 + t29 * t96 + t20 - t222;];
tau_reg = t1;
