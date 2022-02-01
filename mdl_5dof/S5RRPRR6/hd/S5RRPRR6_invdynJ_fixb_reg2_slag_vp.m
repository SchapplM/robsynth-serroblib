% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:18
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:17:25
% EndTime: 2022-01-20 11:17:34
% DurationCPUTime: 4.06s
% Computational Cost: add. (5893->406), mult. (8918->543), div. (0->0), fcn. (5766->14), ass. (0->250)
t181 = qJDD(1) + qJDD(2);
t184 = qJD(1) + qJD(2);
t189 = sin(pkin(9));
t195 = cos(qJ(5));
t196 = cos(qJ(4));
t291 = t189 * t196;
t248 = t181 * t291;
t191 = sin(qJ(5));
t192 = sin(qJ(4));
t284 = t191 * t192;
t283 = t191 * t196;
t130 = t195 * t192 + t283;
t263 = qJD(4) + qJD(5);
t335 = t263 * t130;
t31 = (t181 * t284 + t184 * t335) * t189 - t195 * t248;
t193 = sin(qJ(2));
t190 = cos(pkin(9));
t269 = qJD(3) * t190;
t197 = cos(qJ(2));
t285 = t190 * t197;
t317 = pkin(1) * qJD(1);
t137 = -t190 * pkin(3) - t189 * pkin(7) - pkin(2);
t286 = t190 * t196;
t100 = qJ(3) * t286 + t192 * t137;
t332 = qJD(4) * t100;
t310 = -t192 * t269 - (-t192 * t285 + t193 * t196) * t317 - t332;
t267 = qJD(4) * t196;
t339 = (t192 * t193 + t196 * t285) * t317 - t137 * t267 - t196 * t269;
t268 = qJD(4) * t192;
t243 = t189 * t268;
t153 = pkin(8) * t243;
t337 = t153 + t310;
t287 = t190 * t192;
t250 = qJ(3) * t287;
t260 = pkin(8) * t291;
t336 = -(-t250 - t260) * qJD(4) + t339;
t251 = t189 * t284;
t334 = -qJD(5) * t251 - t191 * t243;
t258 = t193 * t317;
t133 = t184 * qJ(3) + t258;
t264 = qJDD(1) * t193;
t270 = qJD(2) * t197;
t90 = t181 * qJ(3) + t184 * qJD(3) + (qJD(1) * t270 + t264) * pkin(1);
t333 = t90 * qJ(3) + t133 * qJD(3);
t272 = qJD(1) * t197;
t257 = pkin(1) * t272;
t222 = qJD(3) - t257;
t182 = t189 ^ 2;
t331 = t133 * (qJD(4) * t190 + t182 * t184);
t324 = t197 * pkin(1);
t118 = t137 - t324;
t166 = t193 * pkin(1) + qJ(3);
t73 = t192 * t118 + t166 * t286;
t188 = qJ(1) + qJ(2);
t176 = sin(t188);
t178 = cos(t188);
t106 = t176 * t287 + t178 * t196;
t108 = t176 * t196 - t178 * t287;
t330 = -g(1) * t108 + g(2) * t106;
t289 = t190 * t181;
t142 = -qJDD(4) + t289;
t275 = -qJD(2) * t258 + qJDD(1) * t324;
t241 = qJDD(3) - t275;
t65 = t137 * t181 + t241;
t62 = t196 * t65;
t223 = -t90 * t287 + t62;
t247 = t133 * t286;
t296 = t184 * t189;
t81 = t137 * t184 + t222;
t18 = -pkin(8) * t248 - t142 * pkin(4) + (-t247 + (pkin(8) * t296 - t81) * t192) * qJD(4) + t223;
t244 = t184 * t267;
t299 = t181 * t192;
t215 = t244 + t299;
t242 = t190 * t268;
t262 = t192 * t65 + t81 * t267 + t90 * t286;
t25 = -t133 * t242 + t262;
t19 = -pkin(8) * t189 * t215 + t25;
t266 = qJD(5) * t191;
t288 = t190 * t184;
t143 = -qJD(4) + t288;
t252 = t184 * t291;
t52 = -t133 * t287 + t196 * t81;
t44 = -pkin(8) * t252 + t52;
t34 = -t143 * pkin(4) + t44;
t292 = t189 * t192;
t261 = pkin(8) * t292;
t53 = t192 * t81 + t247;
t45 = -t184 * t261 + t53;
t4 = (qJD(5) * t34 + t19) * t195 + t191 * t18 - t45 * t266;
t170 = g(1) * t176;
t169 = g(2) * t178;
t327 = g(3) * t189;
t326 = t181 * pkin(2);
t194 = sin(qJ(1));
t325 = t194 * pkin(1);
t213 = t184 * t130;
t85 = t189 * t213;
t88 = -t184 * t251 + t195 * t252;
t323 = t88 * t85;
t123 = t196 * t137;
t67 = -t260 + t123 + (-qJ(3) * t192 - pkin(4)) * t190;
t80 = -t261 + t100;
t37 = -t191 * t80 + t195 * t67;
t322 = qJD(5) * t37 + t337 * t191 - t336 * t195;
t38 = t191 * t67 + t195 * t80;
t321 = -qJD(5) * t38 + t336 * t191 + t337 * t195;
t157 = pkin(1) * t270 + qJD(3);
t308 = t133 * t182;
t313 = t90 * t166;
t320 = t157 * t308 + t182 * t313;
t129 = t195 * t196 - t284;
t319 = (t263 - t288) * t129;
t318 = t190 * t213 - t335;
t82 = t182 * t90;
t316 = t191 * t45;
t315 = t195 * t45;
t312 = t333 * t182;
t311 = -qJ(3) * t242 - t339;
t309 = qJD(4) * t53;
t134 = -qJDD(5) + t142;
t307 = t134 * t190;
t306 = t142 * t190;
t305 = t157 * t184;
t304 = t176 * t190;
t303 = t176 * t192;
t302 = t178 * t189;
t301 = t178 * t190;
t300 = t178 * t192;
t298 = t181 * t196;
t180 = t184 ^ 2;
t297 = t182 * t180;
t295 = t184 * t192;
t294 = t184 * t196;
t293 = t189 * t181;
t290 = t189 * (-pkin(8) - pkin(7));
t282 = t192 * t196;
t280 = t334 * t184;
t279 = -g(2) * t302 + t189 * t170;
t278 = t178 * pkin(2) + t176 * qJ(3);
t277 = g(1) * t178 + g(2) * t176;
t276 = t170 - t169;
t183 = t190 ^ 2;
t274 = t182 + t183;
t185 = t192 ^ 2;
t186 = t196 ^ 2;
t273 = t185 - t186;
t271 = qJD(2) * t193;
t259 = pkin(1) * t271;
t256 = t52 * t243 + t279;
t102 = t241 - t326;
t255 = t102 * t189 - t279;
t253 = t181 * t283;
t249 = t166 * t287;
t246 = t118 * t267 + t157 * t286 + t192 * t259;
t245 = t184 * t271;
t165 = t178 * qJ(3);
t240 = -t176 * pkin(2) + t165;
t239 = -t102 - t169;
t236 = t274 * t197;
t235 = t274 * t181;
t234 = t142 + t289;
t233 = t184 * (-qJD(4) - t143);
t232 = t183 * t90 - t277 + t82;
t231 = t263 * t196;
t230 = t184 * t258;
t229 = t282 * t297;
t228 = pkin(3) * t301 + pkin(7) * t302 + t278;
t227 = -t275 - t276;
t226 = t192 * t244;
t198 = cos(qJ(1));
t221 = g(1) * t194 - g(2) * t198;
t21 = t191 * t34 + t315;
t114 = t196 * t118;
t56 = -t260 + t114 + (-t166 * t192 - pkin(4)) * t190;
t63 = -t261 + t73;
t29 = -t191 * t63 + t195 * t56;
t30 = t191 * t56 + t195 * t63;
t220 = t192 * t52 - t196 * t53;
t111 = t130 * t189;
t112 = t129 * t189;
t20 = t195 * t34 - t316;
t5 = -t21 * qJD(5) + t195 * t18 - t191 * t19;
t57 = t335 * t189;
t58 = t195 * t189 * t231 + t334;
t219 = -t4 * t111 - t5 * t112 + t20 * t57 - t21 * t58 + t279;
t218 = qJD(4) * (t143 + t288);
t154 = t189 * pkin(4) * t267;
t216 = t222 * t189 + t154;
t214 = -g(1) * t106 - g(2) * t108 + t25 * t190 + t196 * t82;
t212 = -t230 - t326;
t51 = (pkin(4) * t215 + t90) * t189;
t84 = (pkin(4) * t295 + t133) * t189;
t187 = qJ(4) + qJ(5);
t175 = sin(t187);
t177 = cos(t187);
t94 = t175 * t304 + t178 * t177;
t96 = -t175 * t301 + t176 * t177;
t211 = -g(1) * t94 - g(2) * t96 + t51 * t112 + t4 * t190 - t84 * t57;
t173 = -pkin(2) - t324;
t210 = pkin(1) * t245 + t173 * t181;
t209 = t137 * t170;
t172 = t196 * pkin(4) + pkin(3);
t208 = pkin(4) * t303 + t172 * t301 - t178 * t290 + t278;
t95 = t178 * t175 - t177 * t304;
t97 = t176 * t175 + t177 * t301;
t207 = -g(1) * t95 - g(2) * t97 + t51 * t111 - t5 * t190 + t84 * t58;
t206 = -t143 ^ 2 - t297;
t205 = pkin(4) * t300 + t165 + (-t172 * t190 - pkin(2) + t290) * t176;
t107 = -t176 * t286 + t300;
t109 = t178 * t286 + t303;
t26 = t223 - t309;
t204 = -g(1) * t107 - g(2) * t109 - t26 * t190 + t192 * t82 + t267 * t308;
t50 = -t73 * qJD(4) - t157 * t287 + t196 * t259;
t202 = g(1) * t97 - g(2) * t95 + t177 * t327 + t84 * t85 - t4;
t201 = -g(1) * t96 + g(2) * t94 + t175 * t327 - t84 * t88 + t5;
t179 = t198 * pkin(1);
t162 = pkin(4) * t292;
t161 = t183 * t181;
t160 = t182 * t181;
t150 = g(1) * t304;
t136 = -qJD(5) + t143;
t135 = 0.2e1 * t189 * t289;
t127 = t189 * qJ(3) + t162;
t126 = -t184 * pkin(2) + t222;
t115 = t189 * t166 + t162;
t101 = t189 * t157 + t154;
t99 = t123 - t250;
t92 = (t181 * t186 - 0.2e1 * t226) * t182;
t91 = (t181 * t185 + 0.2e1 * t226) * t182;
t72 = t114 - t249;
t64 = 0.2e1 * (t273 * t184 * qJD(4) - t181 * t282) * t182;
t49 = -t166 * t242 + t246;
t47 = (t234 * t192 + t196 * t218) * t189;
t46 = (t192 * t218 - t234 * t196) * t189;
t43 = t153 + t50;
t42 = (-t249 - t260) * qJD(4) + t246;
t33 = -t85 ^ 2 + t88 ^ 2;
t32 = (t253 + (t184 * t231 + t299) * t195) * t189 + t280;
t28 = -t88 * t136 + (-t253 + (-t263 * t294 - t299) * t195) * t189 - t280;
t27 = -t85 * t136 - t31;
t23 = t195 * t44 - t316;
t22 = -t191 * t44 - t315;
t14 = t32 * t111 + t85 * t58;
t13 = -t31 * t112 - t88 * t57;
t12 = t111 * t134 + t58 * t136 + t32 * t190;
t11 = -t112 * t134 + t57 * t136 + t31 * t190;
t8 = -t30 * qJD(5) - t191 * t42 + t195 * t43;
t7 = t29 * qJD(5) + t191 * t43 + t195 * t42;
t6 = t31 * t111 - t112 * t32 + t57 * t85 - t88 * t58;
t1 = [0, 0, 0, 0, 0, qJDD(1), t221, g(1) * t198 + g(2) * t194, 0, 0, 0, 0, 0, 0, 0, t181, (t181 * t197 - t245) * pkin(1) - t227, ((-qJDD(1) - t181) * t193 + (-qJD(1) - t184) * t270) * pkin(1) + t277, 0, (t221 + (t193 ^ 2 + t197 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t160, t135, 0, t161, 0, 0, t150 + (-t210 + t239) * t190, t189 * t210 + t255, t166 * t235 + t274 * t305 + t232, t102 * t173 + t126 * t259 - g(1) * (t240 - t325) - g(2) * (t179 + t278) + (t133 * t157 + t313) * t183 + t320, t92, t64, t46, t91, t47, t306, -t72 * t142 - t50 * t143 + (t157 * t295 + t166 * t215) * t182 + t204, t73 * t142 + t49 * t143 + ((t166 * t181 + t305) * t196 + (-t166 * t184 - t133) * t268) * t182 + t214, ((-t181 * t73 - t25 + (qJD(4) * t72 - t49) * t184) * t192 + (-t181 * t72 - t184 * t50 - t26 + (-t184 * t73 - t53) * qJD(4)) * t196) * t189 + t256, t25 * t73 + t53 * t49 + t26 * t72 + t52 * t50 - g(1) * (t165 - t325) - g(2) * (t179 + t228) - t209 + t320, t13, t6, t11, t14, t12, t307, t101 * t85 + t115 * t32 - t29 * t134 - t8 * t136 + t207, t101 * t88 - t115 * t31 + t30 * t134 + t7 * t136 + t211, t29 * t31 - t30 * t32 - t7 * t85 - t8 * t88 + t219, t4 * t30 + t21 * t7 + t5 * t29 + t20 * t8 + t51 * t115 + t84 * t101 - g(1) * (t205 - t325) - g(2) * (t179 + t208); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181, -t227 + t230, (-t264 + (-qJD(2) + t184) * t272) * pkin(1) + t277, 0, 0, t160, t135, 0, t161, 0, 0, t150 + (-t212 + t239) * t190, t189 * t212 + t255, qJ(3) * t235 + (qJD(3) * t274 - t236 * t317) * t184 + t232, -t102 * pkin(2) - g(1) * t240 - g(2) * t278 + t333 * t183 + (-t126 * t193 - t133 * t236) * t317 + t312, t92, t64, t46, t91, t47, t306, -t99 * t142 - t310 * t143 + (qJ(3) * t299 + (qJ(3) * t267 + t222 * t192) * t184) * t182 + t204, t100 * t142 + t311 * t143 + (qJ(3) * t298 - t133 * t268 + (-qJ(3) * t268 + t196 * t222) * t184) * t182 + t214, ((-t100 * t181 - t25) * t192 + (-t181 * t99 - t26 - t309) * t196 + ((-t310 - t332) * t196 + (qJD(4) * t99 - t311) * t192) * t184) * t189 + t256, -g(1) * t165 - g(2) * t228 + t25 * t100 - t257 * t308 + t26 * t99 + t310 * t52 + t311 * t53 - t209 + t312, t13, t6, t11, t14, t12, t307, t127 * t32 - t37 * t134 - t321 * t136 + t216 * t85 + t207, -t127 * t31 + t38 * t134 + t322 * t136 + t216 * t88 + t211, t37 * t31 - t38 * t32 - t321 * t88 - t322 * t85 + t219, -g(1) * t205 - g(2) * t208 + t51 * t127 + t321 * t20 + t322 * t21 + t216 * t84 + t5 * t37 + t4 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t289, t293, -t274 * t180, -t133 * t184 * t274 + t102 - t276, 0, 0, 0, 0, 0, 0, -t196 * t142 + t192 * t206, t192 * t142 + t196 * t206, (-t185 - t186) * t293, t25 * t192 + t26 * t196 - t220 * qJD(4) + (t190 * t220 - t308) * t184 - t276, 0, 0, 0, 0, 0, 0, -t129 * t134 - t318 * t136 - t85 * t296, t130 * t134 + t319 * t136 - t88 * t296, t129 * t31 - t130 * t32 - t318 * t88 - t319 * t85, t5 * t129 + t4 * t130 + t318 * t20 + t319 * t21 - t84 * t296 - t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, -t273 * t297, (t192 * t233 + t298) * t189, -t229, (t196 * t233 - t299) * t189, -t142, -t53 * t143 + t62 - t196 * t331 + (-qJD(4) * t81 - t190 * t90 + t327) * t192 + t330, g(1) * t109 - g(2) * t107 + g(3) * t291 - t52 * t143 + t192 * t331 - t262, 0, 0, t323, t33, t27, -t323, t28, -t134, t22 * t136 + (-t134 * t195 + t136 * t266 - t252 * t85) * pkin(4) + t201, -t23 * t136 + (qJD(5) * t136 * t195 + t134 * t191 - t252 * t88) * pkin(4) + t202, (t21 + t22) * t88 + (-t20 + t23) * t85 + (-t191 * t32 + t195 * t31 + (t191 * t88 - t195 * t85) * qJD(5)) * pkin(4), -t20 * t22 - t21 * t23 + (t4 * t191 + t5 * t195 + (g(3) * t192 - t294 * t84) * t189 + (-t20 * t191 + t21 * t195) * qJD(5) + t330) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t323, t33, t27, -t323, t28, -t134, -t21 * t136 + t201, -t20 * t136 + t202, 0, 0;];
tau_reg = t1;
