% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRRR9
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,d5,theta1]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRRR9_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR9_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR9_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR9_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:21:18
% EndTime: 2019-12-05 17:21:32
% DurationCPUTime: 6.26s
% Computational Cost: add. (5633->520), mult. (13255->748), div. (0->0), fcn. (10123->14), ass. (0->236)
t187 = sin(qJ(3));
t191 = cos(qJ(3));
t226 = pkin(3) * t187 - pkin(8) * t191;
t137 = t226 * qJD(3);
t227 = t191 * pkin(3) + t187 * pkin(8);
t141 = -pkin(2) - t227;
t186 = sin(qJ(4));
t190 = cos(qJ(4));
t264 = t190 * qJD(3);
t268 = qJD(4) * t190;
t270 = qJD(4) * t186;
t183 = sin(pkin(5));
t278 = qJD(1) * t183;
t192 = cos(qJ(2));
t283 = t191 * t192;
t188 = sin(qJ(2));
t287 = t186 * t188;
t315 = t186 * t137 + t141 * t268 + (-t187 * t264 - t191 * t270) * pkin(7) - (t190 * t283 + t287) * t278;
t265 = t186 * qJD(3);
t342 = t187 * pkin(7) * t265 + t190 * t137 - (-t186 * t283 + t188 * t190) * t278;
t284 = t190 * t191;
t167 = pkin(7) * t284;
t218 = pkin(4) * t187 - pkin(9) * t284;
t341 = t218 * qJD(3) + (-t167 + (pkin(9) * t187 - t141) * t186) * qJD(4) + t342;
t247 = t191 * t265;
t209 = t187 * t268 + t247;
t340 = -pkin(9) * t209 + t315;
t274 = qJD(2) * t187;
t129 = t186 * t274 - t264;
t131 = t190 * t274 + t265;
t185 = sin(qJ(5));
t189 = cos(qJ(5));
t221 = t185 * t129 - t189 * t131;
t65 = t189 * t129 + t185 * t131;
t318 = t65 * t221;
t325 = pkin(9) + pkin(8);
t253 = qJD(4) * t325;
t263 = t191 * qJD(2);
t134 = t226 * qJD(2);
t139 = qJD(2) * pkin(7) + t188 * t278;
t125 = t187 * t139;
t184 = cos(pkin(5));
t277 = qJD(1) * t184;
t95 = t191 * t277 - t125;
t58 = t186 * t134 + t190 * t95;
t339 = t58 + (-pkin(9) * t263 + t253) * t186;
t57 = t190 * t134 - t186 * t95;
t338 = -qJD(2) * t218 - t190 * t253 - t57;
t262 = qJD(1) * qJD(2);
t241 = t192 * t262;
t303 = qJDD(2) * pkin(7);
t337 = t303 + (qJDD(1) * t188 + t241) * t183 + qJD(3) * t277;
t166 = -qJD(4) + t263;
t259 = qJDD(1) * t184;
t256 = -t187 * t259 - t337 * t191;
t272 = qJD(3) * t187;
t40 = -t139 * t272 - t256;
t34 = qJDD(3) * pkin(8) + t40;
t242 = t188 * t262;
t290 = t183 * t192;
t223 = -qJDD(1) * t290 + t183 * t242;
t56 = qJD(2) * t137 + qJDD(2) * t141 + t223;
t276 = qJD(1) * t187;
t246 = t184 * t276;
t96 = t191 * t139 + t246;
t87 = qJD(3) * pkin(8) + t96;
t275 = qJD(1) * t192;
t248 = t183 * t275;
t98 = qJD(2) * t141 - t248;
t215 = -t186 * t56 - t190 * t34 - t98 * t268 + t87 * t270;
t47 = -t186 * t87 + t190 * t98;
t336 = t47 * t166 - t215;
t261 = qJD(2) * qJD(3);
t239 = t191 * t261;
t258 = t187 * qJDD(2);
t335 = qJD(3) * qJD(4) + t239 + t258;
t334 = t221 ^ 2 - t65 ^ 2;
t156 = -qJD(5) + t166;
t269 = qJD(4) * t187;
t238 = qJD(2) * t269;
t230 = t335 * t186 + t190 * t238;
t210 = t190 * qJDD(3) - t230;
t266 = qJD(5) * t189;
t267 = qJD(5) * t185;
t63 = (-qJDD(3) + t238) * t186 - t335 * t190;
t16 = t129 * t266 + t131 * t267 - t185 * t210 + t189 * t63;
t333 = -t65 * t156 - t16;
t25 = -t131 * pkin(9) + t47;
t21 = -t166 * pkin(4) + t25;
t48 = t186 * t98 + t190 * t87;
t26 = -t129 * pkin(9) + t48;
t11 = -t48 * qJD(4) - t186 * t34 + t190 * t56;
t174 = t191 * qJDD(2);
t237 = t187 * t261;
t126 = qJDD(4) - t174 + t237;
t4 = t126 * pkin(4) + t63 * pkin(9) + t11;
t7 = pkin(9) * t210 - t215;
t1 = (qJD(5) * t21 + t7) * t189 + t185 * t4 - t26 * t267;
t305 = cos(pkin(10));
t232 = t305 * t192;
t182 = sin(pkin(10));
t294 = t182 * t188;
t110 = -t184 * t232 + t294;
t233 = t305 * t188;
t293 = t182 * t192;
t112 = t184 * t293 + t233;
t291 = t183 * t191;
t118 = t184 * t187 + t188 * t291;
t181 = qJ(4) + qJ(5);
t176 = sin(t181);
t177 = cos(t181);
t86 = -qJD(3) * pkin(3) - t95;
t59 = t129 * pkin(4) + t86;
t111 = t184 * t233 + t293;
t234 = t183 * t305;
t76 = t111 * t191 - t187 * t234;
t113 = -t184 * t294 + t232;
t78 = t182 * t183 * t187 + t113 * t191;
t332 = t59 * t65 - g(1) * (-t112 * t176 - t78 * t177) - g(2) * (-t110 * t176 - t76 * t177) - g(3) * (-t118 * t177 + t176 * t290) - t1;
t311 = t189 * t26;
t13 = t185 * t21 + t311;
t2 = -t13 * qJD(5) - t185 * t7 + t189 * t4;
t331 = t59 * t221 - g(1) * (t112 * t177 - t78 * t176) - g(2) * (t110 * t177 - t76 * t176) - g(3) * (-t118 * t176 - t177 * t290) + t2;
t202 = t221 * qJD(5) + t185 * t63 + t189 * t210;
t330 = t156 * t221 + t202;
t329 = t48 * t166 - t11;
t102 = t186 * t141 + t167;
t81 = -t118 * t186 - t190 * t290;
t328 = -g(1) * (t112 * t190 - t78 * t186) - g(2) * (t110 * t190 - t76 * t186) - g(3) * t81;
t245 = t191 * t264;
t327 = -t186 * t269 + t245;
t257 = qJD(4) + qJD(5);
t194 = qJD(3) ^ 2;
t225 = g(1) * t112 + g(2) * t110;
t304 = qJDD(2) * pkin(2);
t99 = t223 - t304;
t326 = -pkin(7) * t194 + t183 * (-g(3) * t192 + t242) + t225 + t304 - t99;
t128 = t190 * t141;
t285 = t187 * t190;
t71 = -pkin(9) * t285 + t128 + (-pkin(7) * t186 - pkin(4)) * t191;
t288 = t186 * t187;
t83 = -pkin(9) * t288 + t102;
t28 = -t185 * t83 + t189 * t71;
t323 = t28 * qJD(5) + t341 * t185 + t340 * t189;
t29 = t185 * t71 + t189 * t83;
t322 = -t29 * qJD(5) - t340 * t185 + t341 * t189;
t321 = pkin(4) * t186;
t148 = t325 * t186;
t149 = t325 * t190;
t89 = -t189 * t148 - t185 * t149;
t317 = t89 * qJD(5) + t338 * t185 - t339 * t189;
t90 = -t185 * t148 + t189 * t149;
t316 = -t90 * qJD(5) + t339 * t185 + t338 * t189;
t314 = -qJD(4) * t102 + t342;
t313 = qJD(2) * pkin(2);
t312 = t185 * t26;
t308 = t63 * t186;
t133 = t185 * t190 + t189 * t186;
t74 = t257 * t133;
t307 = t133 * t263 - t74;
t289 = t185 * t186;
t132 = -t189 * t190 + t289;
t306 = -t132 * t263 - t189 * t268 - t190 * t266 + t257 * t289;
t302 = qJDD(3) * pkin(3);
t301 = t129 * t166;
t300 = t129 * t186;
t299 = t131 * t129;
t298 = t131 * t166;
t297 = t131 * t190;
t296 = t176 * t191;
t295 = t177 * t191;
t292 = t183 * t188;
t286 = t186 * t191;
t282 = qJDD(1) - g(3);
t179 = t187 ^ 2;
t180 = t191 ^ 2;
t280 = t179 - t180;
t279 = t179 + t180;
t273 = qJD(2) * t188;
t271 = qJD(3) * t191;
t195 = qJD(2) ^ 2;
t255 = t187 * t195 * t191;
t254 = pkin(7) + t321;
t252 = t187 * t275;
t251 = t183 * t273;
t250 = qJD(2) * t290;
t243 = g(3) * (pkin(2) * t290 + pkin(7) * t292);
t231 = t139 * t271 + t337 * t187 - t191 * t259;
t229 = t191 * t237;
t228 = pkin(4) * t270 - t246 - (qJD(2) * t321 + t139) * t191;
t224 = g(1) * t113 + g(2) * t111;
t216 = -t118 * t190 + t186 * t290;
t30 = t185 * t216 + t189 * t81;
t31 = t185 * t81 - t189 * t216;
t222 = -t186 * t48 - t190 * t47;
t172 = t190 * pkin(4) + pkin(3);
t220 = t172 * t191 + t187 * t325;
t219 = qJDD(2) * t192 - t188 * t195;
t117 = -t184 * t191 + t187 * t292;
t214 = t186 * t126 - t166 * t268;
t213 = t190 * t126 + t166 * t270;
t35 = t231 - t302;
t75 = -t111 * t187 - t191 * t234;
t77 = -t113 * t187 + t182 * t291;
t212 = g(1) * t77 + g(2) * t75 - g(3) * t117;
t211 = g(1) * t78 + g(2) * t76 + g(3) * t118;
t208 = -t212 - t35;
t207 = g(3) * t290 - t225;
t206 = -g(3) * t292 - t224;
t205 = -pkin(8) * t126 - t166 * t86;
t204 = t210 * t190;
t201 = -qJD(4) * pkin(8) * t166 - t208;
t140 = -t248 - t313;
t200 = -pkin(7) * qJDD(3) + (t140 + t248 - t313) * qJD(3);
t197 = t231 * t187 + t40 * t191 + (-t187 * t96 - t191 * t95) * qJD(3) - t224;
t138 = t254 * t187;
t123 = qJDD(5) + t126;
t108 = t112 * pkin(2);
t107 = t110 * pkin(2);
t106 = t132 * t187;
t105 = t133 * t187;
t101 = -pkin(7) * t286 + t128;
t97 = pkin(4) * t209 + pkin(7) * t271;
t80 = qJD(3) * t118 + t187 * t250;
t79 = -qJD(3) * t117 + t191 * t250;
t37 = -t267 * t288 + (t257 * t285 + t247) * t189 + t327 * t185;
t36 = t185 * t247 + t187 * t74 - t189 * t245;
t24 = qJD(4) * t81 + t186 * t251 + t79 * t190;
t23 = qJD(4) * t216 - t79 * t186 + t190 * t251;
t20 = -pkin(4) * t210 + t35;
t15 = t189 * t25 - t312;
t14 = -t185 * t25 - t311;
t12 = t189 * t21 - t312;
t6 = -t31 * qJD(5) - t185 * t24 + t189 * t23;
t5 = t30 * qJD(5) + t185 * t23 + t189 * t24;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t282, 0, 0, 0, 0, 0, 0, t219 * t183, (-qJDD(2) * t188 - t192 * t195) * t183, 0, -g(3) + (t184 ^ 2 + (t188 ^ 2 + t192 ^ 2) * t183 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -t80 * qJD(3) - t117 * qJDD(3) + (t191 * t219 - t192 * t237) * t183, -t79 * qJD(3) - t118 * qJDD(3) + (-t187 * t219 - t192 * t239) * t183, (t117 * t187 + t118 * t191) * qJDD(2) + (t187 * t80 + t191 * t79 + (t117 * t191 - t118 * t187) * qJD(3)) * qJD(2), t231 * t117 + t40 * t118 + t96 * t79 - t95 * t80 - g(3) + (t140 * t273 - t192 * t99) * t183, 0, 0, 0, 0, 0, 0, -t117 * t210 + t81 * t126 + t80 * t129 - t23 * t166, -t117 * t63 + t126 * t216 + t80 * t131 + t24 * t166, -t24 * t129 - t23 * t131 - t210 * t216 + t81 * t63, t11 * t81 + t117 * t35 + t215 * t216 + t23 * t47 + t24 * t48 + t80 * t86 - g(3), 0, 0, 0, 0, 0, 0, -t117 * t202 + t30 * t123 - t6 * t156 + t80 * t65, -t117 * t16 - t31 * t123 + t5 * t156 - t221 * t80, t16 * t30 + t202 * t31 + t221 * t6 - t5 * t65, t1 * t31 + t117 * t20 + t12 * t6 + t13 * t5 + t2 * t30 + t59 * t80 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t282 * t290 + t225, -t282 * t292 + t224, 0, 0, t179 * qJDD(2) + 0.2e1 * t229, 0.2e1 * t187 * t174 - 0.2e1 * t280 * t261, t187 * qJDD(3) + t194 * t191, t180 * qJDD(2) - 0.2e1 * t229, qJDD(3) * t191 - t194 * t187, 0, t200 * t187 + t191 * t326, -t187 * t326 + t200 * t191, t279 * t303 + (-g(3) * t188 - t241 * t279) * t183 + t197, -t99 * pkin(2) + g(1) * t108 + g(2) * t107 - t243 + (-t140 * t188 + (t187 * t95 - t191 * t96) * t192) * t278 + t197 * pkin(7), t131 * t327 - t63 * t285, (-t190 * t129 - t131 * t186) * t271 + (t204 + t308 + (-t297 + t300) * qJD(4)) * t187, (-t166 * t264 + t63) * t191 + (qJD(3) * t131 + t213) * t187, t129 * t209 - t210 * t288, (t166 * t265 - t210) * t191 + (-t129 * qJD(3) - t214) * t187, -t126 * t191 - t166 * t272, t101 * t126 - t314 * t166 + t206 * t186 + (-t11 + (pkin(7) * t129 + t86 * t186) * qJD(3) - t207 * t190) * t191 + (-pkin(7) * t210 + t47 * qJD(3) - t129 * t248 + t35 * t186 + t268 * t86) * t187, -t102 * t126 + t315 * t166 + t206 * t190 + (-t215 + (pkin(7) * t131 + t190 * t86) * qJD(3) + t207 * t186) * t191 + (-pkin(7) * t63 - t48 * qJD(3) - t131 * t248 + t35 * t190 - t270 * t86) * t187, t102 * t210 + t101 * t63 - t314 * t131 - t315 * t129 + t222 * t271 + (t215 * t186 - t11 * t190 + (t186 * t47 - t190 * t48) * qJD(4) - t207) * t187, -t215 * t102 + t11 * t101 - g(1) * (-t112 * t227 - t108) - g(2) * (-t110 * t227 - t107) - t243 + t315 * t48 + t314 * t47 + (-g(3) * t227 - t276 * t86) * t290 + (t35 * t187 + t271 * t86 - t224) * pkin(7), t106 * t16 + t221 * t36, t105 * t16 - t106 * t202 + t221 * t37 + t36 * t65, -t106 * t123 + t36 * t156 + t16 * t191 - t221 * t272, -t105 * t202 + t37 * t65, -t105 * t123 + t37 * t156 - t191 * t202 - t272 * t65, -t123 * t191 - t156 * t272, t97 * t65 - t138 * t202 + t20 * t105 + t59 * t37 + t28 * t123 - t2 * t191 + t12 * t272 - g(1) * (-t112 * t295 + t113 * t176) - g(2) * (-t110 * t295 + t111 * t176) - t322 * t156 + (-t65 * t252 - g(3) * (t176 * t188 + t177 * t283)) * t183, -t97 * t221 - t138 * t16 - t20 * t106 - t59 * t36 - t29 * t123 + t1 * t191 - t13 * t272 - g(1) * (t112 * t296 + t113 * t177) - g(2) * (t110 * t296 + t111 * t177) + t323 * t156 + (t221 * t252 - g(3) * (-t176 * t283 + t177 * t188)) * t183, -t1 * t105 + t2 * t106 + t12 * t36 - t13 * t37 + t28 * t16 - t207 * t187 + t202 * t29 + t221 * t322 - t323 * t65, t1 * t29 + t2 * t28 + t20 * t138 + t59 * t97 - g(1) * (-t112 * t220 + t113 * t254 - t108) - g(2) * (-t110 * t220 + t111 * t254 - t107) - t243 + t323 * t13 + t322 * t12 + (-g(3) * pkin(4) * t287 + (-g(3) * t220 - t276 * t59) * t192) * t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t255, t280 * t195, t258, t255, t174, qJDD(3), t96 * qJD(3) - t140 * t274 - t212 - t231, -t140 * t263 + (t95 + t125) * qJD(3) + t211 + t256, 0, 0, -t166 * t297 - t308, (-t63 + t301) * t190 + (t210 + t298) * t186, (-t131 * t187 + t166 * t284) * qJD(2) + t214, -t166 * t300 + t204, (t129 * t187 - t166 * t286) * qJD(2) + t213, t166 * t274, -pkin(3) * t230 + t57 * t166 - t47 * t274 - t96 * t129 + t205 * t186 + (-t201 + t302) * t190, pkin(3) * t63 - t96 * t131 - t58 * t166 + t186 * t201 + t190 * t205 + t274 * t48, t58 * t129 + t57 * t131 + ((t131 * qJD(4) + t210) * pkin(8) + t336) * t190 + ((t129 * qJD(4) - t63) * pkin(8) + t329) * t186 - t211, -t47 * t57 - t48 * t58 - t86 * t96 + t208 * pkin(3) + (qJD(4) * t222 - t11 * t186 - t190 * t215 - t211) * pkin(8), -t133 * t16 + t221 * t306, t132 * t16 + t133 * t202 - t221 * t307 + t306 * t65, t133 * t123 + t306 * t156 + t221 * t274, -t132 * t202 - t307 * t65, -t132 * t123 - t307 * t156 + t65 * t274, t156 * t274, -t12 * t274 + t89 * t123 + t20 * t132 - t316 * t156 + t172 * t202 - t212 * t177 + t228 * t65 - t307 * t59, -t90 * t123 + t13 * t274 + t20 * t133 + t317 * t156 + t172 * t16 + t212 * t176 - t221 * t228 - t306 * t59, -t1 * t132 + t306 * t12 + t307 * t13 - t133 * t2 + t16 * t89 + t202 * t90 + t221 * t316 - t317 * t65 - t211, t1 * t90 + t2 * t89 - t20 * t172 - g(1) * (t77 * t172 + t325 * t78) - g(2) * (t75 * t172 + t325 * t76) - g(3) * (-t117 * t172 + t118 * t325) + t228 * t59 + t317 * t13 + t316 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t299, -t129 ^ 2 + t131 ^ 2, -t63 - t301, -t299, t210 - t298, t126, -t86 * t131 + t328 - t329, t86 * t129 - g(1) * (-t112 * t186 - t78 * t190) - g(2) * (-t110 * t186 - t76 * t190) - g(3) * t216 - t336, 0, 0, -t318, t334, t333, t318, t330, t123, t14 * t156 + (t123 * t189 - t131 * t65 + t156 * t267) * pkin(4) + t331, -t15 * t156 + (-t123 * t185 + t131 * t221 + t156 * t266) * pkin(4) + t332, -t12 * t65 - t13 * t221 - t14 * t221 + t15 * t65 + (t16 * t189 + t202 * t185 + (-t185 * t221 - t189 * t65) * qJD(5)) * pkin(4), -t12 * t14 - t13 * t15 + (t1 * t185 + t2 * t189 - t59 * t131 + (-t12 * t185 + t13 * t189) * qJD(5) + t328) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t318, t334, t333, t318, t330, t123, -t13 * t156 + t331, -t12 * t156 + t332, 0, 0;];
tau_reg = t3;
