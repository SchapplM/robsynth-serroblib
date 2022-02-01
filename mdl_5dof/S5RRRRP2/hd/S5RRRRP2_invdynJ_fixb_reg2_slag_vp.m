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
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 11:49:28
% EndTime: 2022-01-20 11:49:35
% DurationCPUTime: 2.77s
% Computational Cost: add. (4719->377), mult. (7259->435), div. (0->0), fcn. (4699->12), ass. (0->229)
t192 = sin(qJ(4));
t196 = cos(qJ(3));
t193 = sin(qJ(3));
t316 = cos(qJ(4));
t256 = t316 * t193;
t121 = t192 * t196 + t256;
t187 = qJD(1) + qJD(2);
t98 = t121 * t187;
t293 = t98 * qJ(5);
t194 = sin(qJ(2));
t298 = pkin(1) * qJD(1);
t262 = t194 * t298;
t129 = t187 * pkin(7) + t262;
t249 = pkin(8) * t187 + t129;
t87 = t249 * t196;
t76 = t192 * t87;
t86 = t249 * t193;
t79 = qJD(3) * pkin(3) - t86;
t44 = t316 * t79 - t76;
t28 = t44 - t293;
t197 = cos(qJ(2));
t188 = t193 ^ 2;
t189 = t196 ^ 2;
t274 = t188 + t189;
t323 = t274 * t197;
t327 = t187 * t323;
t184 = qJDD(1) + qJDD(2);
t267 = qJDD(1) * t194;
t271 = qJD(2) * t197;
t102 = t184 * pkin(7) + (qJD(1) * t271 + t267) * pkin(1);
t244 = t274 * t102;
t307 = t197 * pkin(1);
t276 = -qJD(2) * t262 + qJDD(1) * t307;
t310 = t184 * pkin(2);
t101 = -t276 - t310;
t191 = qJ(1) + qJ(2);
t179 = cos(t191);
t313 = g(2) * t179;
t326 = t101 + t313;
t251 = t316 * qJD(4);
t255 = t316 * t196;
t325 = -qJD(3) * t255 - t196 * t251;
t273 = qJD(1) * t197;
t261 = pkin(1) * t273;
t309 = t187 * pkin(2);
t130 = -t261 - t309;
t324 = t129 * t323 + t130 * t194;
t199 = -pkin(8) - pkin(7);
t257 = qJD(3) * t199;
t126 = t193 * t257;
t127 = t196 * t257;
t148 = t199 * t193;
t180 = t196 * pkin(8);
t149 = t196 * pkin(7) + t180;
t283 = t192 * t193;
t214 = t255 - t283;
t268 = qJD(4) * t192;
t300 = t316 * t126 + t192 * t127 + t148 * t251 - t149 * t268 - t214 * t261;
t81 = t192 * t148 + t316 * t149;
t299 = -t81 * qJD(4) + t121 * t261 - t192 * t126 + t316 * t127;
t168 = t194 * pkin(1) + pkin(7);
t305 = -pkin(8) - t168;
t117 = t305 * t193;
t118 = t196 * t168 + t180;
t64 = t192 * t117 + t316 * t118;
t183 = qJDD(3) + qJDD(4);
t186 = qJD(3) + qJD(4);
t236 = pkin(3) * t251;
t315 = pkin(3) * t192;
t322 = -t183 * t315 - t186 * t236;
t177 = sin(t191);
t321 = g(1) * t179 + g(2) * t177;
t175 = t183 * pkin(4);
t224 = t186 * t283;
t281 = t196 * t184;
t235 = -t184 * t256 + t325 * t187 - t192 * t281;
t35 = t187 * t224 + t235;
t297 = t35 * qJ(5);
t320 = t175 + t297;
t164 = g(1) * t177;
t319 = t313 - t164;
t190 = qJ(3) + qJ(4);
t176 = sin(t190);
t287 = t176 * t179;
t288 = t176 * t177;
t178 = cos(t190);
t312 = g(3) * t178;
t318 = g(1) * t287 + g(2) * t288 - t312;
t317 = t98 ^ 2;
t195 = sin(qJ(1));
t314 = g(1) * t195;
t311 = g(3) * t196;
t308 = t195 * pkin(1);
t259 = t187 * t283;
t96 = -t187 * t255 + t259;
t306 = t98 * t96;
t68 = t186 * t121;
t291 = -t68 * qJ(5) + qJD(5) * t214;
t304 = t291 + t300;
t67 = t224 + t325;
t221 = t67 * qJ(5) - t121 * qJD(5);
t303 = t221 + t299;
t27 = t186 * pkin(4) + t28;
t302 = t27 - t28;
t282 = t193 * t184;
t225 = -t184 * t255 + t192 * t282;
t36 = t187 * t68 + t225;
t301 = -t96 * t236 - t36 * t315;
t48 = -t316 * t86 - t76;
t296 = t36 * qJ(5);
t295 = t96 * qJ(5);
t294 = t96 * t186;
t290 = t121 * qJ(5);
t286 = t177 * t178;
t285 = t178 * t179;
t284 = t187 * t193;
t170 = t196 * pkin(3) + pkin(2);
t100 = -t170 * t187 - t261;
t250 = t96 * pkin(4) + qJD(5);
t58 = t100 + t250;
t280 = qJD(5) + t58;
t270 = qJD(3) * t193;
t279 = t130 * t270 + t196 * t164;
t278 = t179 * pkin(2) + t177 * pkin(7);
t275 = t188 - t189;
t272 = qJD(2) * t194;
t269 = qJD(3) * t196;
t266 = pkin(3) * t284;
t265 = t316 * pkin(3);
t264 = pkin(1) * t271;
t263 = pkin(3) * t268;
t173 = pkin(3) * t270;
t260 = t130 * t269 + t326 * t193;
t78 = t316 * t87;
t182 = t187 ^ 2;
t258 = t193 * t182 * t196;
t254 = t187 * t272;
t253 = t187 * t270;
t252 = t187 * t269;
t61 = t68 * pkin(4) + t173;
t41 = -t129 * t269 + qJDD(3) * pkin(3) - t193 * t102 + (-t252 - t282) * pkin(8);
t46 = -t129 * t270 + t196 * t102 + (-t253 + t281) * pkin(8);
t247 = -t192 * t46 + t316 * t41;
t47 = t192 * t86 - t78;
t246 = qJD(3) * t305;
t8 = t192 * t41 + t79 * t251 - t87 * t268 + t316 * t46;
t243 = t274 * t184;
t63 = t316 * t117 - t192 * t118;
t166 = pkin(4) * t178;
t125 = t166 + t170;
t185 = qJ(5) - t199;
t241 = t179 * t125 + t177 * t185;
t80 = t316 * t148 - t192 * t149;
t240 = -t177 * t125 + t185 * t179;
t239 = t179 * t170 - t177 * t199;
t238 = -t321 + t244;
t237 = t187 * t262;
t234 = t193 * t252;
t231 = t61 - t262;
t230 = -g(1) * t288 + g(2) * t287;
t229 = g(1) * t286 - g(2) * t285;
t228 = g(1) * (-t177 * pkin(2) + t179 * pkin(7));
t198 = cos(qJ(1));
t226 = -g(2) * t198 + t314;
t45 = t192 * t79 + t78;
t29 = t45 - t295;
t9 = -qJD(4) * t45 + t247;
t3 = -t98 * qJD(5) + t320 + t9;
t4 = -t96 * qJD(5) - t296 + t8;
t223 = -t3 * t121 + t214 * t4 + t27 * t67 - t29 * t68 - t321;
t222 = -t9 * t121 + t214 * t8 + t44 * t67 - t45 * t68 - t321;
t220 = -t177 * t170 - t179 * t199;
t95 = -pkin(4) * t214 - t170;
t219 = -t276 + t319;
t62 = pkin(3) * t253 - t170 * t184 - t276;
t20 = t36 * pkin(4) + qJDD(5) + t62;
t218 = t20 * t121 - t58 * t67 + t230;
t217 = -t100 * t67 + t62 * t121 + t230;
t216 = -t20 * t214 + t58 * t68 + t229;
t215 = t100 * t68 - t214 * t62 + t229;
t213 = g(1) * t285 + g(2) * t286 + g(3) * t176 - t8;
t83 = t193 * t246 + t196 * t264;
t84 = -t193 * t264 + t196 * t246;
t21 = t117 * t251 - t118 * t268 + t192 * t84 + t316 * t83;
t212 = -t262 + t173;
t211 = -t130 * t187 - t102 + t321;
t210 = t100 * t96 + t213;
t200 = qJD(3) ^ 2;
t209 = pkin(7) * t200 - t237 - t310;
t171 = -pkin(2) - t307;
t208 = pkin(1) * t254 + t168 * t200 + t171 * t184;
t207 = -t186 * t259 - t235;
t206 = -pkin(7) * qJDD(3) + (t261 - t309) * qJD(3);
t205 = -qJDD(3) * t168 + (t171 * t187 - t264) * qJD(3);
t22 = -t64 * qJD(4) - t192 * t83 + t316 * t84;
t204 = t280 * t96 + t213 + t296;
t203 = t9 + t318;
t202 = -t100 * t98 + t203;
t181 = t198 * pkin(1);
t174 = pkin(1) * t272;
t169 = t265 + pkin(4);
t143 = -t170 - t307;
t142 = qJDD(3) * t196 - t200 * t193;
t141 = qJDD(3) * t193 + t200 * t196;
t128 = t174 + t173;
t116 = t214 * qJ(5);
t104 = t189 * t184 - 0.2e1 * t234;
t103 = t188 * t184 + 0.2e1 * t234;
t94 = t96 ^ 2;
t85 = t95 - t307;
t70 = -0.2e1 * t275 * t187 * qJD(3) + 0.2e1 * t193 * t281;
t69 = t98 * pkin(4) + t266;
t60 = t116 + t81;
t59 = t80 - t290;
t57 = t174 + t61;
t52 = t116 + t64;
t51 = t63 - t290;
t50 = t183 * t214 - t68 * t186;
t49 = t121 * t183 - t67 * t186;
t39 = -t94 + t317;
t31 = -t293 + t48;
t30 = t47 + t295;
t23 = t207 + t294;
t13 = -t214 * t36 + t96 * t68;
t12 = -t35 * t121 - t98 * t67;
t11 = t22 + t221;
t10 = t21 + t291;
t5 = -t121 * t36 - t214 * t35 + t67 * t96 - t98 * t68;
t1 = [0, 0, 0, 0, 0, qJDD(1), t226, g(1) * t198 + g(2) * t195, 0, 0, 0, 0, 0, 0, 0, t184, (t184 * t197 - t254) * pkin(1) - t219, ((-qJDD(1) - t184) * t194 + (-qJD(1) - t187) * t271) * pkin(1) + t321, 0, (t226 + (t194 ^ 2 + t197 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t103, t70, t141, t104, t142, 0, t205 * t193 + (-t208 - t326) * t196 + t279, t205 * t196 + (t208 - t164) * t193 + t260, pkin(1) * qJD(2) * t327 + t168 * t243 + t238, t101 * t171 - t228 - g(2) * (t181 + t278) + t168 * t244 + (t324 * qJD(2) + t314) * pkin(1), t12, t5, t49, t13, t50, 0, t128 * t96 + t143 * t36 + t63 * t183 + t22 * t186 + t215, t128 * t98 - t143 * t35 - t64 * t183 - t21 * t186 + t217, -t21 * t96 - t22 * t98 + t63 * t35 - t64 * t36 + t222, t8 * t64 + t45 * t21 + t9 * t63 + t44 * t22 + t62 * t143 + t100 * t128 - g(1) * (t220 - t308) - g(2) * (t181 + t239), t12, t5, t49, t13, t50, 0, t11 * t186 + t51 * t183 + t85 * t36 + t57 * t96 + t216, -t10 * t186 - t52 * t183 - t85 * t35 + t57 * t98 + t218, -t10 * t96 - t11 * t98 + t51 * t35 - t52 * t36 + t223, t4 * t52 + t29 * t10 + t3 * t51 + t27 * t11 + t20 * t85 + t58 * t57 - g(1) * (t240 - t308) - g(2) * (t181 + t241); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t184, -t219 + t237, (-t267 + (-qJD(2) + t187) * t273) * pkin(1) + t321, 0, 0, t103, t70, t141, t104, t142, 0, t206 * t193 + (-t209 - t326) * t196 + t279, t206 * t196 + (t209 - t164) * t193 + t260, pkin(7) * t243 - t298 * t327 + t238, -t101 * pkin(2) + pkin(7) * t244 - g(2) * t278 - t324 * t298 - t228, t12, t5, t49, t13, t50, 0, -t170 * t36 + t80 * t183 + t299 * t186 + t212 * t96 + t215, t170 * t35 - t81 * t183 - t300 * t186 + t212 * t98 + t217, -t299 * t98 - t300 * t96 + t80 * t35 - t81 * t36 + t222, -g(1) * t220 - g(2) * t239 + t212 * t100 - t62 * t170 + t299 * t44 + t300 * t45 + t8 * t81 + t9 * t80, t12, t5, t49, t13, t50, 0, t59 * t183 + t303 * t186 + t231 * t96 + t95 * t36 + t216, -t60 * t183 - t304 * t186 + t231 * t98 - t95 * t35 + t218, -t303 * t98 - t304 * t96 + t59 * t35 - t60 * t36 + t223, -g(1) * t240 - g(2) * t241 + t20 * t95 + t231 * t58 + t303 * t27 + t304 * t29 + t3 * t59 + t4 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t258, t275 * t182, t282, t258, t281, qJDD(3), t193 * t211 - t311, g(3) * t193 + t196 * t211, 0, 0, t306, t39, t23, -t306, -t225, t183, -t47 * t186 + (t316 * t183 - t186 * t268 - t96 * t284) * pkin(3) + t202, t48 * t186 - t98 * t266 + t210 + t322, t35 * t265 + (-t44 + t48) * t96 + (t45 + t47 + t263) * t98 + t301, -t44 * t47 - t45 * t48 + (t316 * t9 - t311 + t192 * t8 + (-t192 * t44 + t316 * t45) * qJD(4) + (-t100 * t187 + t321) * t193) * pkin(3), t306, t39, t23, -t306, -t225, t183, t169 * t183 - t30 * t186 - t69 * t96 - t280 * t98 + (-t78 + (-pkin(3) * t186 - t79) * t192) * qJD(4) + t247 + t318 + t320, t31 * t186 - t69 * t98 + t204 + t322, t169 * t35 + (-t27 + t31) * t96 + (t29 + t30 + t263) * t98 + t301, -g(3) * t166 + t3 * t169 - t27 * t30 - t29 * t31 - t58 * t69 - t321 * (-t193 * pkin(3) - pkin(4) * t176) + (-t311 + t4 * t192 + (-t192 * t27 + t316 * t29) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t306, t39, t23, -t306, -t225, t183, t45 * t186 + t202, t44 * t186 + t210, 0, 0, t306, t39, t23, -t306, -t225, t183, t297 + t29 * t186 + 0.2e1 * t175 + (-t250 - t58) * t98 + t203, -t317 * pkin(4) + t28 * t186 + t204, t35 * pkin(4) - t302 * t96, t302 * t29 + (t176 * t321 - t58 * t98 + t3 - t312) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98 * t186 + t36, t207 - t294, -t94 - t317, t27 * t98 + t29 * t96 + t20 + t319;];
tau_reg = t1;
