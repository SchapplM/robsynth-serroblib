% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRP9
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP9_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:07:34
% EndTime: 2019-12-31 20:07:43
% DurationCPUTime: 5.21s
% Computational Cost: add. (5129->480), mult. (11785->611), div. (0->0), fcn. (8256->10), ass. (0->234)
t188 = sin(qJ(2));
t267 = t188 * qJDD(1);
t169 = pkin(6) * t267;
t190 = cos(qJ(2));
t269 = qJD(1) * qJD(2);
t255 = t190 * t269;
t296 = qJDD(2) * pkin(2);
t108 = pkin(6) * t255 + qJDD(3) + t169 - t296;
t180 = g(3) * t190;
t189 = sin(qJ(1));
t191 = cos(qJ(1));
t241 = g(1) * t191 + g(2) * t189;
t207 = t241 * t188 - t180;
t198 = t108 - t207;
t175 = t190 * qJDD(1);
t254 = t188 * t269;
t215 = t254 - t175;
t127 = qJDD(4) + t215;
t181 = pkin(8) + qJ(4);
t173 = sin(t181);
t184 = sin(pkin(8));
t312 = pkin(7) + qJ(3);
t141 = t312 * t184;
t185 = cos(pkin(8));
t142 = t312 * t185;
t187 = sin(qJ(4));
t323 = cos(qJ(4));
t84 = -t187 * t141 + t323 * t142;
t334 = -t127 * t84 - t207 * t173;
t277 = qJD(1) * t188;
t257 = t184 * t277;
t274 = qJD(2) * t185;
t123 = -t257 + t274;
t261 = t185 * t277;
t275 = qJD(2) * t184;
t124 = t261 + t275;
t69 = -t323 * t123 + t124 * t187;
t333 = t69 ^ 2;
t276 = qJD(1) * t190;
t164 = -qJD(4) + t276;
t332 = t164 * t69;
t264 = t323 * t185;
t247 = t190 * t264;
t262 = t184 * t276;
t256 = qJD(4) * t323;
t271 = qJD(4) * t187;
t330 = -t184 * t271 + t185 * t256;
t300 = qJD(1) * t247 - t187 * t262 - t330;
t129 = t323 * t184 + t187 * t185;
t114 = t129 * qJD(4);
t210 = t190 * t129;
t299 = -qJD(1) * t210 + t114;
t220 = t187 * t123 + t323 * t124;
t325 = t220 ^ 2;
t291 = t185 * t190;
t229 = pkin(3) * t188 - pkin(7) * t291;
t236 = pkin(2) * t188 - qJ(3) * t190;
t131 = t236 * qJD(1);
t89 = pkin(6) * t257 + t185 * t131;
t61 = qJD(1) * t229 + t89;
t115 = t184 * t131;
t292 = t185 * t188;
t293 = t184 * t190;
t216 = -pkin(6) * t292 - pkin(7) * t293;
t74 = qJD(1) * t216 + t115;
t308 = -qJD(3) * t129 - qJD(4) * t84 + t187 * t74 - t323 * t61;
t155 = t184 * t267;
t281 = t184 * t255 + t155;
t233 = qJDD(2) * t185 - t281;
t268 = qJDD(2) * t184;
t92 = t268 + (t255 + t267) * t185;
t28 = qJD(4) * t220 + t187 * t92 - t323 * t233;
t302 = t220 * t164;
t331 = t28 + t302;
t218 = -t187 * t184 + t264;
t182 = t188 ^ 2;
t183 = t190 ^ 2;
t279 = t182 - t183;
t248 = qJD(1) * t279;
t27 = -t123 * t256 + t124 * t271 - t187 * t233 - t323 * t92;
t329 = t28 * pkin(4) + t27 * qJ(5) - t220 * qJD(5);
t328 = t127 * t218 + t299 * t164 + t69 * t277;
t237 = pkin(2) * t190 + qJ(3) * t188;
t137 = -pkin(1) - t237;
t122 = t185 * t137;
t75 = -pkin(7) * t292 + t122 + (-pkin(6) * t184 - pkin(3)) * t190;
t294 = t184 * t188;
t94 = pkin(6) * t291 + t184 * t137;
t81 = -pkin(7) * t294 + t94;
t306 = t187 * t75 + t323 * t81;
t110 = qJD(2) * t236 - qJD(3) * t188;
t273 = qJD(2) * t188;
t266 = pkin(6) * t273;
t79 = t185 * t110 + t184 * t266;
t53 = qJD(2) * t229 + t79;
t97 = t184 * t110;
t62 = qJD(2) * t216 + t97;
t12 = -qJD(4) * t306 - t187 * t62 + t323 * t53;
t327 = t129 * t28 + t218 * t27 + t220 * t299 - t300 * t69;
t326 = -0.2e1 * pkin(1);
t322 = pkin(3) * t184;
t321 = pkin(4) * t127;
t320 = pkin(6) * t123;
t319 = pkin(6) * t124;
t318 = g(1) * t189;
t315 = g(2) * t191;
t314 = g(3) * t188;
t313 = t220 * t69;
t171 = pkin(6) * t276;
t118 = pkin(3) * t262 + t171;
t311 = t299 * pkin(4) + t300 * qJ(5) - qJD(5) * t129 - t118;
t34 = t187 * t61 + t323 * t74;
t30 = qJ(5) * t277 + t34;
t219 = -t323 * t141 - t187 * t142;
t51 = qJD(3) * t218 + qJD(4) * t219;
t310 = t51 - t30;
t309 = t51 - t34;
t307 = pkin(4) * t277 - t308;
t66 = qJD(1) * t110 + qJDD(1) * t137;
t99 = -pkin(6) * t215 + qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t41 = t184 * t66 + t185 * t99;
t117 = t137 * qJD(1);
t143 = qJD(2) * qJ(3) + t171;
t76 = t185 * t117 - t143 * t184;
t45 = -pkin(3) * t276 - pkin(7) * t124 + t76;
t77 = t184 * t117 + t185 * t143;
t49 = pkin(7) * t123 + t77;
t20 = t187 * t45 + t323 * t49;
t304 = t164 * t20;
t301 = t92 * t184;
t298 = pkin(6) * qJDD(1);
t297 = qJ(5) * t127;
t193 = qJD(1) ^ 2;
t295 = t183 * t193;
t289 = t188 * t189;
t288 = t188 * t191;
t174 = cos(t181);
t287 = t189 * t174;
t286 = t189 * t190;
t168 = t185 * pkin(3) + pkin(2);
t144 = t190 * t168;
t285 = t190 * t191;
t284 = t191 * t173;
t19 = -t187 * t49 + t323 * t45;
t282 = qJD(5) - t19;
t272 = qJD(2) * t190;
t260 = t184 * t272;
t119 = pkin(3) * t260 + pkin(6) * t272;
t132 = pkin(3) * t294 + t188 * pkin(6);
t280 = t191 * pkin(1) + t189 * pkin(6);
t278 = t182 + t183;
t136 = -qJD(2) * pkin(2) + pkin(6) * t277 + qJD(3);
t270 = qJD(3) - t136;
t265 = -t325 + t333;
t263 = t124 * t276;
t258 = t164 * t277;
t253 = t185 * t267;
t252 = t185 * t175;
t251 = t184 * t175;
t40 = -t184 * t99 + t185 * t66;
t23 = pkin(3) * t215 - pkin(7) * t92 + t40;
t29 = pkin(7) * t233 + t41;
t249 = t187 * t29 - t323 * t23 + t49 * t256 + t45 * t271;
t245 = t190 * t254;
t166 = g(1) * t289;
t244 = -g(2) * t288 + t166;
t100 = t173 * t286 + t174 * t191;
t102 = t190 * t284 - t287;
t243 = g(1) * t100 - g(2) * t102;
t101 = t174 * t286 - t284;
t103 = t173 * t189 + t174 * t285;
t242 = g(1) * t101 - g(2) * t103;
t240 = -t315 + t318;
t239 = -t325 - t333;
t106 = t129 * t188;
t57 = qJD(2) * t210 + t330 * t188;
t238 = t106 * t28 + t57 * t69;
t235 = pkin(4) * t174 + qJ(5) * t173;
t232 = t185 * t123 - t184 * t124;
t230 = qJD(1) * (-t123 + t274);
t38 = -t187 * t81 + t323 * t75;
t224 = t233 * pkin(3);
t223 = -pkin(6) * qJDD(2) + t269 * t326;
t3 = t187 * t23 + t45 * t256 - t49 * t271 + t323 * t29;
t11 = t187 * t53 + t75 * t256 - t81 * t271 + t323 * t62;
t221 = -t218 * t28 + t299 * t69;
t217 = t233 * t185;
t214 = pkin(1) * t193 + t241;
t88 = -pkin(3) * t123 + t136;
t192 = qJD(2) ^ 2;
t213 = pkin(6) * t192 + qJDD(1) * t326 + t315;
t212 = t168 * t285 + t189 * t322 + t288 * t312 + t280;
t211 = g(2) * t188 * t287 + t127 * t219 + (g(1) * t288 - t180) * t174;
t178 = t191 * pkin(6);
t209 = -t312 * t289 + t191 * t322 + t178 + (-pkin(1) - t144) * t189;
t208 = t27 - t332;
t205 = -t190 * t241 - t314;
t204 = g(1) * t102 + g(2) * t100 + t173 * t314 - t249;
t107 = t218 * t188;
t56 = -qJD(2) * t247 + t114 * t188 + t187 * t260;
t203 = t106 * t27 - t107 * t28 - t220 * t57 + t56 * t69;
t201 = t106 * t127 - t164 * t57 - t190 * t28 + t69 * t273;
t54 = t108 - t224;
t24 = pkin(4) * t69 - qJ(5) * t220 + t88;
t200 = t220 * t24 + qJDD(5) - t204;
t199 = t219 * t27 - t84 * t28 - t51 * t69 + t205;
t197 = -g(1) * t103 - g(2) * t101 - t174 * t314 + t3;
t196 = t28 - t302;
t195 = -t224 + t198;
t158 = t188 * t193 * t190;
t120 = qJDD(1) * t183 - 0.2e1 * t245;
t93 = -pkin(6) * t293 + t122;
t90 = -pkin(6) * t261 + t115;
t82 = -t127 * t190 - t164 * t273;
t80 = -t185 * t266 + t97;
t67 = -pkin(4) * t218 - qJ(5) * t129 - t168;
t46 = pkin(4) * t106 - qJ(5) * t107 + t132;
t37 = pkin(4) * t220 + qJ(5) * t69;
t36 = t190 * pkin(4) - t38;
t35 = -qJ(5) * t190 + t306;
t18 = -t164 * qJ(5) + t20;
t16 = t164 * pkin(4) + t282;
t15 = t127 * t129 + t300 * t164 - t220 * t277;
t14 = pkin(4) * t57 + qJ(5) * t56 - qJD(5) * t107 + t119;
t13 = -t27 - t332;
t10 = -pkin(4) * t273 - t12;
t9 = qJ(5) * t273 - qJD(5) * t190 + t11;
t8 = -t107 * t27 - t220 * t56;
t7 = -t129 * t27 - t220 * t300;
t6 = t107 * t127 + t164 * t56 + t190 * t27 + t220 * t273;
t5 = t54 + t329;
t2 = qJDD(5) + t249 - t321;
t1 = -qJD(5) * t164 + t297 + t3;
t4 = [0, 0, 0, 0, 0, qJDD(1), t240, t241, 0, 0, qJDD(1) * t182 + 0.2e1 * t245, -0.2e1 * qJD(2) * t248 + 0.2e1 * t188 * t175, qJDD(2) * t188 + t190 * t192, t120, qJDD(2) * t190 - t188 * t192, 0, t223 * t188 + (-t213 + t318) * t190, t188 * t213 + t190 * t223 - t166, 0.2e1 * t278 * t298 - t241, -g(1) * (-pkin(1) * t189 + t178) - g(2) * t280 + (pkin(6) ^ 2 * t278 + pkin(1) ^ 2) * qJDD(1), (t124 * t272 + t188 * t92) * t185, (t217 - t301) * t188 + t232 * t272, (-t92 - t253) * t190 + (t124 * t188 + t185 * t248) * qJD(2), (-t123 * t272 - t188 * t233) * t184, (-t233 + t155) * t190 + (t188 * t123 - t184 * t248) * qJD(2), t120, -t241 * t184 + (-pkin(6) * t233 + t108 * t184 + (qJD(1) * t93 + t76) * qJD(2)) * t188 + (-t79 * qJD(1) - t93 * qJDD(1) - t40 + t240 * t185 + (t136 * t184 - t320) * qJD(2)) * t190, -t241 * t185 + (pkin(6) * t92 + t108 * t185 + (-qJD(1) * t94 - t77) * qJD(2)) * t188 + (t80 * qJD(1) + t94 * qJDD(1) + t41 - t240 * t184 + (t136 * t185 + t319) * qJD(2)) * t190, t80 * t123 + t94 * t233 - t79 * t124 - t93 * t92 + t166 + (-t184 * t77 - t185 * t76) * t272 + (-t184 * t41 - t185 * t40 - t315) * t188, t41 * t94 + t77 * t80 + t40 * t93 + t76 * t79 - g(1) * t178 - g(2) * (t191 * t237 + t280) - t137 * t318 + (t108 * t188 + t136 * t272) * pkin(6), t8, t203, t6, t238, -t201, t82, t106 * t54 + t119 * t69 - t12 * t164 + t127 * t38 + t132 * t28 + t19 * t273 + t190 * t249 + t57 * t88 + t242, t107 * t54 + t11 * t164 + t119 * t220 - t127 * t306 - t132 * t27 + t190 * t3 - t20 * t273 - t56 * t88 - t243, -t106 * t3 + t107 * t249 - t11 * t69 - t12 * t220 + t19 * t56 - t20 * t57 + t27 * t38 - t28 * t306 + t244, -g(1) * t209 - g(2) * t212 + t20 * t11 + t88 * t119 + t19 * t12 + t54 * t132 - t249 * t38 + t3 * t306, t8, t6, -t203, t82, t201, t238, t10 * t164 + t106 * t5 - t127 * t36 + t14 * t69 - t16 * t273 + t190 * t2 + t24 * t57 + t28 * t46 + t242, -t1 * t106 + t10 * t220 + t107 * t2 - t16 * t56 - t18 * t57 - t27 * t36 - t28 * t35 - t69 * t9 + t244, -t1 * t190 - t107 * t5 + t127 * t35 - t14 * t220 - t164 * t9 + t18 * t273 + t24 * t56 + t27 * t46 + t243, t1 * t35 + t18 * t9 + t5 * t46 + t24 * t14 + t2 * t36 + t16 * t10 - g(1) * (-pkin(4) * t101 - qJ(5) * t100 + t209) - g(2) * (pkin(4) * t103 + qJ(5) * t102 + t212); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t158, t279 * t193, t267, t158, t175, qJDD(2), t188 * t214 - t169 - t180, t314 + (t214 - t298) * t190, 0, 0, -t185 * t263 + t301, t184 * t233 + t92 * t185 - t232 * t276, -t251 + t185 * t295 + (-t124 + t275) * t277, t123 * t262 + t217, -t184 * t295 + t188 * t230 - t252, t158, qJ(3) * t251 - pkin(2) * t281 + (-t198 + t296) * t185 + ((-qJ(3) * t275 - t76) * t188 + (t184 * t270 + t320 + t89) * t190) * qJD(1), qJ(3) * t252 - pkin(2) * t92 + t198 * t184 + ((-qJ(3) * t274 + t77) * t188 + (t185 * t270 - t319 - t90) * t190) * qJD(1), -t90 * t123 + t89 * t124 + (qJ(3) * t233 + qJD(3) * t123 + t276 * t76 + t41) * t185 + (qJ(3) * t92 + qJD(3) * t124 + t276 * t77 - t40) * t184 + t205, -t136 * t171 - t76 * t89 - t77 * t90 + (-t184 * t76 + t185 * t77) * qJD(3) - t198 * pkin(2) + (-t40 * t184 + t41 * t185 + t205) * qJ(3), t7, -t327, t15, t221, t328, t258, -t118 * t69 - t308 * t164 - t168 * t28 - t19 * t277 - t218 * t54 + t299 * t88 + t211, -t118 * t220 + t129 * t54 + t309 * t164 + t168 * t27 + t20 * t277 - t300 * t88 + t334, t129 * t249 + t300 * t19 - t299 * t20 + t218 * t3 - t220 * t308 + t34 * t69 + t199, t3 * t84 - t249 * t219 - t54 * t168 - t88 * t118 - g(3) * (t188 * t312 + t144) + t309 * t20 + t308 * t19 + t241 * (t168 * t188 - t190 * t312), t7, t15, t327, t258, -t328, t221, t16 * t277 + t307 * t164 - t218 * t5 + t299 * t24 + t28 * t67 + t311 * t69 + t211, t1 * t218 + t129 * t2 - t300 * t16 - t299 * t18 + t220 * t307 + t30 * t69 + t199, -t129 * t5 - t310 * t164 - t18 * t277 - t220 * t311 + t300 * t24 + t27 * t67 - t334, -g(3) * t144 + t1 * t84 - t2 * t219 + t5 * t67 + t311 * t24 + t310 * t18 + t307 * t16 + (-g(3) * t235 - t241 * t312) * t190 + (-g(3) * t312 + t241 * (t168 + t235)) * t188; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t233 - t263, t190 * t230 + t253 + t268, -t123 ^ 2 - t124 ^ 2, -t123 * t77 + t124 * t76 + t198, 0, 0, 0, 0, 0, 0, t196, -t208, t239, t19 * t220 + t20 * t69 + t195, 0, 0, 0, 0, 0, 0, t196, t239, t208, -t16 * t220 + t18 * t69 + t195 + t329; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t313, -t265, t13, -t313, -t331, t127, -t220 * t88 + t204 - t304, -t164 * t19 + t69 * t88 - t197, 0, 0, t313, t13, t265, t127, t331, -t313, -t37 * t69 - t200 - t304 + 0.2e1 * t321, pkin(4) * t27 - qJ(5) * t28 + (t18 - t20) * t220 + (t16 - t282) * t69, 0.2e1 * t297 - t24 * t69 + t37 * t220 + (-0.2e1 * qJD(5) + t19) * t164 + t197, t1 * qJ(5) - t2 * pkin(4) - t24 * t37 - t16 * t20 - g(1) * (-pkin(4) * t102 + qJ(5) * t103) - g(2) * (-pkin(4) * t100 + qJ(5) * t101) - (-pkin(4) * t173 + qJ(5) * t174) * t314 + t282 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t127 + t313, t13, -t164 ^ 2 - t325, t164 * t18 + t200 - t321;];
tau_reg = t4;
