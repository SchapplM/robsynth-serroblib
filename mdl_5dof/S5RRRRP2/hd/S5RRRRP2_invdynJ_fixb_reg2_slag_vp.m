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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:11:51
% EndTime: 2020-01-03 12:11:55
% DurationCPUTime: 2.49s
% Computational Cost: add. (4719->375), mult. (7259->430), div. (0->0), fcn. (4699->12), ass. (0->226)
t193 = sin(qJ(4));
t197 = cos(qJ(3));
t194 = sin(qJ(3));
t314 = cos(qJ(4));
t253 = t314 * t194;
t122 = t193 * t197 + t253;
t188 = qJD(1) + qJD(2);
t98 = t122 * t188;
t290 = t98 * qJ(5);
t195 = sin(qJ(2));
t296 = pkin(1) * qJD(1);
t258 = t195 * t296;
t130 = pkin(7) * t188 + t258;
t246 = pkin(8) * t188 + t130;
t87 = t246 * t197;
t76 = t193 * t87;
t86 = t246 * t194;
t79 = qJD(3) * pkin(3) - t86;
t44 = t314 * t79 - t76;
t28 = t44 - t290;
t198 = cos(qJ(2));
t305 = t198 * pkin(1);
t273 = -qJD(2) * t258 + qJDD(1) * t305;
t185 = qJDD(1) + qJDD(2);
t307 = t185 * pkin(2);
t101 = -t273 - t307;
t192 = qJ(1) + qJ(2);
t177 = sin(t192);
t179 = cos(t192);
t226 = g(2) * t179 + g(3) * t177;
t324 = t101 + t226;
t189 = t194 ^ 2;
t190 = t197 ^ 2;
t271 = t189 + t190;
t320 = t271 * t198;
t323 = t188 * t320;
t264 = qJDD(1) * t195;
t268 = qJD(2) * t198;
t102 = t185 * pkin(7) + (qJD(1) * t268 + t264) * pkin(1);
t242 = t271 * t102;
t163 = g(3) * t179;
t274 = -g(2) * t177 + t163;
t248 = qJD(4) * t314;
t252 = t314 * t197;
t322 = -qJD(3) * t252 - t197 * t248;
t270 = qJD(1) * t198;
t257 = pkin(1) * t270;
t306 = t188 * pkin(2);
t131 = -t257 - t306;
t321 = t130 * t320 + t131 * t195;
t200 = -pkin(8) - pkin(7);
t254 = qJD(3) * t200;
t127 = t194 * t254;
t128 = t197 * t254;
t148 = t200 * t194;
t181 = t197 * pkin(8);
t149 = pkin(7) * t197 + t181;
t282 = t193 * t194;
t217 = t252 - t282;
t265 = qJD(4) * t193;
t298 = t127 * t314 + t128 * t193 + t148 * t248 - t149 * t265 - t217 * t257;
t81 = t148 * t193 + t149 * t314;
t297 = -qJD(4) * t81 + t122 * t257 - t193 * t127 + t128 * t314;
t168 = pkin(1) * t195 + pkin(7);
t303 = -pkin(8) - t168;
t118 = t303 * t194;
t119 = t168 * t197 + t181;
t64 = t118 * t193 + t119 * t314;
t184 = qJDD(3) + qJDD(4);
t187 = qJD(3) + qJD(4);
t232 = pkin(3) * t248;
t313 = pkin(3) * t193;
t319 = -t184 * t313 - t187 * t232;
t175 = t184 * pkin(4);
t222 = t187 * t282;
t280 = t197 * t185;
t231 = -t185 * t253 + t188 * t322 - t193 * t280;
t35 = t188 * t222 + t231;
t294 = t35 * qJ(5);
t318 = t175 + t294;
t201 = qJD(3) ^ 2;
t317 = pkin(7) * t201 - t307;
t191 = qJ(3) + qJ(4);
t176 = sin(t191);
t284 = t176 * t179;
t285 = t176 * t177;
t178 = cos(t191);
t311 = g(1) * t178;
t316 = g(2) * t285 - g(3) * t284 - t311;
t315 = t98 ^ 2;
t310 = g(1) * t197;
t256 = t188 * t282;
t96 = -t188 * t252 + t256;
t304 = t98 * t96;
t68 = t187 * t122;
t288 = -qJ(5) * t68 + qJD(5) * t217;
t302 = t288 + t298;
t67 = t222 + t322;
t219 = t67 * qJ(5) - t122 * qJD(5);
t301 = t219 + t297;
t27 = pkin(4) * t187 + t28;
t300 = t27 - t28;
t281 = t194 * t185;
t223 = -t185 * t252 + t193 * t281;
t36 = t188 * t68 + t223;
t299 = -t232 * t96 - t313 * t36;
t48 = -t314 * t86 - t76;
t295 = pkin(1) * qJD(2);
t293 = t36 * qJ(5);
t292 = t96 * qJ(5);
t291 = t96 * t187;
t287 = t122 * qJ(5);
t283 = t188 * t194;
t170 = pkin(3) * t197 + pkin(2);
t100 = -t170 * t188 - t257;
t247 = t96 * pkin(4) + qJD(5);
t58 = t100 + t247;
t279 = qJD(5) + t58;
t166 = pkin(4) * t178;
t126 = t166 + t170;
t186 = -qJ(5) + t200;
t278 = t126 * t177 + t179 * t186;
t277 = t170 * t177 + t179 * t200;
t276 = g(2) * t284 + g(3) * t285;
t275 = pkin(2) * t179 + pkin(7) * t177;
t272 = t189 - t190;
t269 = qJD(2) * t195;
t267 = qJD(3) * t194;
t266 = qJD(3) * t197;
t262 = pkin(3) * t283;
t261 = t314 * pkin(3);
t260 = pkin(1) * t268;
t259 = pkin(3) * t265;
t173 = pkin(3) * t267;
t78 = t314 * t87;
t183 = t188 ^ 2;
t255 = t194 * t183 * t197;
t251 = t188 * t269;
t250 = t188 * t267;
t249 = t188 * t266;
t61 = pkin(4) * t68 + t173;
t41 = -t130 * t266 + qJDD(3) * pkin(3) - t194 * t102 + (-t249 - t281) * pkin(8);
t46 = -t130 * t267 + t197 * t102 + (-t250 + t280) * pkin(8);
t245 = -t193 * t46 + t314 * t41;
t47 = t193 * t86 - t78;
t244 = qJD(3) * t303;
t8 = t193 * t41 + t248 * t79 - t265 * t87 + t314 * t46;
t241 = t271 * t185;
t63 = t118 * t314 - t119 * t193;
t239 = t126 * t179 - t177 * t186;
t80 = t148 * t314 - t149 * t193;
t238 = t170 * t179 - t177 * t200;
t62 = pkin(3) * t250 - t170 * t185 - t273;
t20 = t36 * pkin(4) + qJDD(5) + t62;
t237 = t122 * t20 - t58 * t67 + t276;
t236 = -t100 * t67 + t122 * t62 + t276;
t235 = t274 + t242;
t234 = t188 * t258;
t233 = t131 * t266 + t194 * t324;
t230 = t194 * t249;
t227 = t61 - t258;
t196 = sin(qJ(1));
t199 = cos(qJ(1));
t224 = -g(2) * t199 - g(3) * t196;
t45 = t193 * t79 + t78;
t29 = t45 - t292;
t9 = -qJD(4) * t45 + t245;
t3 = -t98 * qJD(5) + t318 + t9;
t4 = -qJD(5) * t96 - t293 + t8;
t221 = -t122 * t3 + t217 * t4 + t27 * t67 - t29 * t68 + t274;
t220 = -t122 * t9 + t217 * t8 + t44 * t67 - t45 * t68 + t274;
t95 = -pkin(4) * t217 - t170;
t218 = t226 * t178;
t83 = t194 * t244 + t197 * t260;
t84 = -t194 * t260 + t197 * t244;
t21 = t118 * t248 - t119 * t265 + t193 * t84 + t314 * t83;
t216 = -t258 + t173;
t215 = -t131 * t188 - t102 - t274;
t214 = -t20 * t217 + t58 * t68 - t218;
t213 = t100 * t68 - t217 * t62 - t218;
t212 = g(1) * t176 - t178 * t274 - t8;
t211 = -t226 + t234;
t171 = -pkin(2) - t305;
t210 = pkin(1) * t251 + t168 * t201 + t171 * t185;
t209 = -t187 * t256 - t231;
t208 = -pkin(7) * qJDD(3) + (t257 - t306) * qJD(3);
t207 = -qJDD(3) * t168 + (t171 * t188 - t260) * qJD(3);
t206 = t100 * t96 + t212;
t22 = -qJD(4) * t64 - t193 * t83 + t314 * t84;
t205 = t279 * t96 + t212 + t293;
t204 = t9 + t316;
t203 = -t100 * t98 + t204;
t182 = t199 * pkin(1);
t180 = t196 * pkin(1);
t174 = pkin(1) * t269;
t169 = t261 + pkin(4);
t161 = t177 * pkin(2);
t143 = -t170 - t305;
t142 = qJDD(3) * t197 - t194 * t201;
t141 = qJDD(3) * t194 + t197 * t201;
t129 = t174 + t173;
t117 = t217 * qJ(5);
t111 = t131 * t267;
t104 = t185 * t190 - 0.2e1 * t230;
t103 = t185 * t189 + 0.2e1 * t230;
t94 = t96 ^ 2;
t85 = t95 - t305;
t70 = -0.2e1 * qJD(3) * t188 * t272 + 0.2e1 * t194 * t280;
t69 = pkin(4) * t98 + t262;
t60 = t117 + t81;
t59 = t80 - t287;
t57 = t174 + t61;
t52 = t117 + t64;
t51 = t63 - t287;
t50 = t184 * t217 - t187 * t68;
t49 = t122 * t184 - t187 * t67;
t39 = -t94 + t315;
t31 = -t290 + t48;
t30 = t47 + t292;
t23 = t209 + t291;
t13 = -t217 * t36 + t68 * t96;
t12 = -t122 * t35 - t67 * t98;
t11 = t219 + t22;
t10 = t21 + t288;
t5 = -t122 * t36 - t217 * t35 + t67 * t96 - t68 * t98;
t1 = [0, 0, 0, 0, 0, qJDD(1), t224, g(2) * t196 - g(3) * t199, 0, 0, 0, 0, 0, 0, 0, t185, (t185 * t198 - t251) * pkin(1) - t226 + t273, ((-qJDD(1) - t185) * t195 + (-qJD(1) - t188) * t268) * pkin(1) - t274, 0, (t224 + (t195 ^ 2 + t198 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t103, t70, t141, t104, t142, 0, t111 + t207 * t194 + (-t210 - t324) * t197, t194 * t210 + t197 * t207 + t233, t168 * t241 + t295 * t323 + t235, t101 * t171 - g(2) * (t182 + t275) - g(3) * (-pkin(7) * t179 + t161 + t180) + t168 * t242 + t321 * t295, t12, t5, t49, t13, t50, 0, t129 * t96 + t143 * t36 + t63 * t184 + t22 * t187 + t213, t129 * t98 - t143 * t35 - t184 * t64 - t187 * t21 + t236, -t21 * t96 - t22 * t98 + t35 * t63 - t36 * t64 + t220, t8 * t64 + t45 * t21 + t9 * t63 + t44 * t22 + t62 * t143 + t100 * t129 - g(2) * (t182 + t238) - g(3) * (t180 + t277), t12, t5, t49, t13, t50, 0, t11 * t187 + t51 * t184 + t85 * t36 + t57 * t96 + t214, -t10 * t187 - t184 * t52 - t35 * t85 + t57 * t98 + t237, -t10 * t96 - t11 * t98 + t35 * t51 - t36 * t52 + t221, t4 * t52 + t29 * t10 + t3 * t51 + t27 * t11 + t20 * t85 + t58 * t57 - g(2) * (t182 + t239) - g(3) * (t180 + t278); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t185, t211 + t273, (-t264 + (-qJD(2) + t188) * t270) * pkin(1) - t274, 0, 0, t103, t70, t141, t104, t142, 0, t111 + t208 * t194 + (-t101 + t211 - t317) * t197, t208 * t197 + (-t234 + t317) * t194 + t233, pkin(7) * t241 - t296 * t323 + t235, -t101 * pkin(2) - g(2) * t275 - g(3) * t161 + (t242 + t163) * pkin(7) - t321 * t296, t12, t5, t49, t13, t50, 0, -t170 * t36 + t80 * t184 + t187 * t297 + t216 * t96 + t213, t170 * t35 - t81 * t184 - t187 * t298 + t216 * t98 + t236, -t297 * t98 - t298 * t96 + t80 * t35 - t81 * t36 + t220, -g(2) * t238 - g(3) * t277 + t100 * t216 - t62 * t170 + t297 * t44 + t298 * t45 + t8 * t81 + t9 * t80, t12, t5, t49, t13, t50, 0, t59 * t184 + t187 * t301 + t227 * t96 + t95 * t36 + t214, -t60 * t184 - t187 * t302 + t227 * t98 - t95 * t35 + t237, -t301 * t98 - t302 * t96 + t59 * t35 - t60 * t36 + t221, -g(2) * t239 - g(3) * t278 + t20 * t95 + t227 * t58 + t27 * t301 + t29 * t302 + t3 * t59 + t4 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t255, t272 * t183, t281, t255, t280, qJDD(3), t194 * t215 - t310, g(1) * t194 + t197 * t215, 0, 0, t304, t39, t23, -t304, -t223, t184, -t47 * t187 + (t184 * t314 - t187 * t265 - t283 * t96) * pkin(3) + t203, t187 * t48 - t262 * t98 + t206 + t319, t35 * t261 + (-t44 + t48) * t96 + (t45 + t47 + t259) * t98 + t299, -t44 * t47 - t45 * t48 + (t314 * t9 - t310 + t193 * t8 + (-t193 * t44 + t314 * t45) * qJD(4) + (-t100 * t188 - t274) * t194) * pkin(3), t304, t39, t23, -t304, -t223, t184, t169 * t184 - t30 * t187 - t69 * t96 - t279 * t98 + (-t78 + (-pkin(3) * t187 - t79) * t193) * qJD(4) + t245 + t316 + t318, t31 * t187 - t69 * t98 + t205 + t319, t169 * t35 + (-t27 + t31) * t96 + (t29 + t30 + t259) * t98 + t299, -g(1) * t166 + t3 * t169 - t27 * t30 - t29 * t31 - t58 * t69 + t274 * (-pkin(3) * t194 - pkin(4) * t176) + (-t310 + t4 * t193 + (-t193 * t27 + t29 * t314) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t304, t39, t23, -t304, -t223, t184, t45 * t187 + t203, t187 * t44 + t206, 0, 0, t304, t39, t23, -t304, -t223, t184, t294 + t29 * t187 + 0.2e1 * t175 + (-t247 - t58) * t98 + t204, -pkin(4) * t315 + t28 * t187 + t205, t35 * pkin(4) - t300 * t96, t300 * t29 + (-t176 * t274 - t58 * t98 + t3 - t311) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98 * t187 + t36, t209 - t291, -t94 - t315, t27 * t98 + t29 * t96 + t20 + t226;];
tau_reg = t1;
