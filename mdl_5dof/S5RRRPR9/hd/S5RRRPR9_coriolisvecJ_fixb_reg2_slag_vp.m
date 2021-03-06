% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR9_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR9_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:24:39
% EndTime: 2019-12-31 21:24:56
% DurationCPUTime: 7.03s
% Computational Cost: add. (8146->458), mult. (20406->643), div. (0->0), fcn. (14322->8), ass. (0->220)
t220 = sin(qJ(5));
t221 = sin(qJ(3));
t222 = sin(qJ(2));
t272 = qJD(1) * t222;
t256 = t221 * t272;
t223 = cos(qJ(3));
t263 = t223 * qJD(2);
t182 = t256 - t263;
t255 = t223 * t272;
t270 = qJD(2) * t221;
t184 = t255 + t270;
t218 = sin(pkin(9));
t219 = cos(pkin(9));
t242 = t182 * t219 + t184 * t218;
t319 = cos(qJ(5));
t230 = t319 * t242;
t233 = -t182 * t218 + t184 * t219;
t264 = qJD(5) * t220;
t224 = cos(qJ(2));
t261 = qJD(1) * qJD(2);
t247 = t224 * t261;
t266 = qJD(3) * t222;
t252 = t221 * t266;
t260 = qJD(2) * qJD(3);
t129 = qJD(1) * t252 + (-t247 - t260) * t223;
t268 = qJD(2) * t224;
t254 = t221 * t268;
t265 = qJD(3) * t223;
t229 = t222 * t265 + t254;
t130 = qJD(1) * t229 + t221 * t260;
t76 = -t129 * t218 + t130 * t219;
t77 = -t129 * t219 - t130 * t218;
t20 = qJD(5) * t230 + t220 * t76 + t233 * t264 - t319 * t77;
t271 = qJD(1) * t224;
t203 = -qJD(3) + t271;
t196 = -qJD(5) + t203;
t58 = t220 * t233 + t230;
t305 = t196 * t58;
t339 = -t20 - t305;
t315 = t58 ^ 2;
t330 = -t220 * t242 + t233 * t319;
t316 = t330 ^ 2;
t338 = -t315 + t316;
t314 = t58 * t330;
t248 = qJD(5) * t319;
t325 = pkin(8) * t233;
t189 = -pkin(2) * t224 - pkin(7) * t222 - pkin(1);
t171 = t189 * qJD(1);
t212 = pkin(6) * t271;
t194 = qJD(2) * pkin(7) + t212;
t121 = t171 * t223 - t194 * t221;
t88 = -qJ(4) * t184 + t121;
t78 = -pkin(3) * t203 + t88;
t122 = t171 * t221 + t194 * t223;
t89 = -qJ(4) * t182 + t122;
t82 = t218 * t89;
t36 = t219 * t78 - t82;
t26 = -pkin(4) * t203 - t325 + t36;
t332 = pkin(8) * t242;
t304 = t219 * t89;
t37 = t218 * t78 + t304;
t28 = t37 - t332;
t208 = t222 * t261;
t236 = pkin(2) * t222 - pkin(7) * t224;
t186 = t236 * qJD(2);
t172 = qJD(1) * t186;
t239 = pkin(6) * t208;
t67 = -qJD(3) * t122 + t172 * t223 + t221 * t239;
t35 = pkin(3) * t208 + t129 * qJ(4) - t184 * qJD(4) + t67;
t267 = qJD(3) * t221;
t66 = t171 * t265 + t172 * t221 - t194 * t267 - t223 * t239;
t40 = -qJ(4) * t130 - qJD(4) * t182 + t66;
t12 = -t218 * t40 + t219 * t35;
t8 = pkin(4) * t208 - pkin(8) * t77 + t12;
t13 = t218 * t35 + t219 * t40;
t9 = -pkin(8) * t76 + t13;
t228 = -t220 * t8 - t248 * t26 + t264 * t28 - t319 * t9;
t211 = pkin(6) * t272;
t306 = qJD(2) * pkin(2);
t193 = t211 - t306;
t127 = t182 * pkin(3) + qJD(4) + t193;
t68 = pkin(4) * t242 + t127;
t337 = t58 * t68 + t228;
t185 = t236 * qJD(1);
t131 = pkin(6) * t256 + t185 * t223;
t282 = t223 * t224;
t231 = pkin(3) * t222 - qJ(4) * t282;
t313 = -qJ(4) - pkin(7);
t244 = qJD(3) * t313;
t336 = qJD(1) * t231 + t221 * qJD(4) - t223 * t244 + t131;
t167 = t221 * t185;
t262 = t223 * qJD(4);
t283 = t222 * t223;
t284 = t221 * t224;
t335 = t167 + (-pkin(6) * t283 - qJ(4) * t284) * qJD(1) - t221 * t244 - t262;
t21 = qJD(5) * t330 + t220 * t77 + t319 * t76;
t303 = t330 * t196;
t334 = -t21 - t303;
t176 = t218 * t223 + t219 * t221;
t165 = t176 * qJD(3);
t277 = -t176 * t271 + t165;
t175 = t218 * t221 - t219 * t223;
t327 = t203 * t175;
t6 = t220 * t26 + t28 * t319;
t2 = -qJD(5) * t6 - t220 * t9 + t319 * t8;
t333 = -t330 * t68 + t2;
t307 = -t218 * t335 + t219 * t336;
t324 = t218 * t336 + t219 * t335;
t331 = t242 * t233;
t329 = -pkin(4) * t272 - pkin(8) * t327 - t307;
t328 = pkin(8) * t277 + t324;
t326 = -0.2e1 * t261;
t323 = t121 * t203 + t66;
t322 = t122 * t203 - t67;
t237 = -t212 + (-t221 * t271 + t267) * pkin(3);
t205 = pkin(6) * t282;
t141 = t189 * t221 + t205;
t318 = pkin(3) * t218;
t317 = pkin(6) * t221;
t191 = t313 * t221;
t192 = t313 * t223;
t125 = t191 * t219 + t192 * t218;
t101 = -pkin(8) * t176 + t125;
t126 = t191 * t218 - t192 * t219;
t102 = -pkin(8) * t175 + t126;
t46 = t101 * t319 - t102 * t220;
t312 = t46 * qJD(5) + t220 * t329 - t319 * t328;
t47 = t101 * t220 + t102 * t319;
t311 = -t47 * qJD(5) + t220 * t328 + t319 * t329;
t309 = t175 * t248 + t176 * t264 + t220 * t277 - t319 * t327;
t112 = -t175 * t220 + t176 * t319;
t308 = t112 * qJD(5) + t220 * t327 + t277 * t319;
t269 = qJD(2) * t222;
t274 = t186 * t223 + t269 * t317;
t53 = -t222 * t262 + t231 * qJD(2) + (-t205 + (qJ(4) * t222 - t189) * t221) * qJD(3) + t274;
t275 = t186 * t221 + t189 * t265;
t62 = (-pkin(6) * qJD(2) - qJ(4) * qJD(3)) * t283 + (-qJD(4) * t222 + (-pkin(6) * qJD(3) - qJ(4) * qJD(2)) * t224) * t221 + t275;
t25 = t218 * t53 + t219 * t62;
t43 = t219 * t88 - t82;
t207 = pkin(3) * t219 + pkin(4);
t159 = t207 * t220 + t318 * t319;
t42 = -t218 * t88 - t304;
t31 = t42 + t332;
t32 = t43 - t325;
t302 = qJD(5) * t159 - t220 * t32 + t31 * t319;
t158 = t207 * t319 - t220 * t318;
t301 = -qJD(5) * t158 + t220 * t31 + t319 * t32;
t300 = t233 * t203;
t299 = t242 * t203;
t298 = t233 ^ 2;
t295 = t129 * t221;
t294 = t130 * t223;
t293 = t182 * t203;
t291 = t184 * t182;
t290 = t184 * t203;
t289 = t193 * t221;
t288 = t193 * t223;
t287 = t203 * t221;
t286 = t203 * t223;
t285 = t221 * t222;
t226 = qJD(1) ^ 2;
t281 = t224 * t226;
t225 = qJD(2) ^ 2;
t280 = t225 * t222;
t279 = t225 * t224;
t278 = pkin(4) * t277 + t237;
t178 = t223 * t189;
t118 = -qJ(4) * t283 + t178 + (-pkin(3) - t317) * t224;
t123 = -qJ(4) * t285 + t141;
t64 = t118 * t218 + t123 * t219;
t187 = pkin(3) * t285 + pkin(6) * t222;
t216 = t222 ^ 2;
t273 = -t224 ^ 2 + t216;
t259 = pkin(6) * t284;
t213 = pkin(6) * t268;
t257 = t222 * t281;
t135 = pkin(3) * t229 + t213;
t210 = -pkin(3) * t223 - pkin(2);
t253 = t224 * t263;
t250 = t203 * t265;
t249 = t203 * t272;
t113 = pkin(3) * t130 + pkin(6) * t247;
t24 = -t218 * t62 + t219 * t53;
t243 = pkin(1) * t326;
t63 = t118 * t219 - t123 * t218;
t241 = t182 + t263;
t240 = -t184 + t270;
t238 = t224 * t208;
t235 = t242 ^ 2;
t234 = -t121 * t223 - t122 * t221;
t232 = qJD(1) * t216 - t203 * t224;
t45 = pkin(4) * t76 + t113;
t153 = t175 * t222;
t44 = -pkin(4) * t224 + pkin(8) * t153 + t63;
t152 = t176 * t222;
t48 = -pkin(8) * t152 + t64;
t22 = -t220 * t48 + t319 * t44;
t23 = t220 * t44 + t319 * t48;
t94 = -t152 * t220 - t153 * t319;
t142 = pkin(4) * t175 + t210;
t140 = t178 - t259;
t139 = (-t203 - t271) * t269;
t132 = -pkin(6) * t255 + t167;
t119 = pkin(4) * t152 + t187;
t111 = t175 * t319 + t176 * t220;
t96 = t165 * t222 + t218 * t254 - t219 * t253;
t95 = t175 * t266 - t176 * t268;
t93 = t152 * t319 - t153 * t220;
t87 = -qJD(3) * t141 + t274;
t86 = (-t222 * t263 - t224 * t267) * pkin(6) + t275;
t79 = pkin(3) * t184 + pkin(4) * t233;
t65 = -pkin(4) * t95 + t135;
t30 = qJD(5) * t94 - t220 * t96 - t319 * t95;
t29 = t152 * t248 - t153 * t264 - t220 * t95 + t319 * t96;
t19 = pkin(8) * t95 + t25;
t18 = pkin(4) * t269 + pkin(8) * t96 + t24;
t5 = -t220 * t28 + t26 * t319;
t4 = -qJD(5) * t23 + t18 * t319 - t19 * t220;
t3 = qJD(5) * t22 + t220 * t18 + t19 * t319;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t238, t273 * t326, t279, -0.2e1 * t238, -t280, 0, -pkin(6) * t279 + t222 * t243, pkin(6) * t280 + t224 * t243, 0, 0, -t129 * t283 + (-t252 + t253) * t184, (-t182 * t223 - t184 * t221) * t268 + (t295 - t294 + (t182 * t221 - t184 * t223) * qJD(3)) * t222, t203 * t252 + t129 * t224 + (t184 * t222 + t223 * t232) * qJD(2), t130 * t285 + t182 * t229, t222 * t250 + t130 * t224 + (-t182 * t222 - t221 * t232) * qJD(2), t139, -t87 * t203 - t67 * t224 + (pkin(6) * t130 + t193 * t265) * t222 + ((pkin(6) * t182 + t289) * t224 + (t121 + (t140 + t259) * qJD(1)) * t222) * qJD(2), t86 * t203 + t66 * t224 + (-pkin(6) * t129 - t193 * t267) * t222 + ((pkin(6) * t184 + t288) * t224 + (-t122 + (-t141 + t205) * qJD(1)) * t222) * qJD(2), t140 * t129 - t141 * t130 - t86 * t182 - t87 * t184 + t234 * t268 + (-t221 * t66 - t223 * t67 + (t121 * t221 - t122 * t223) * qJD(3)) * t222, t121 * t87 + t122 * t86 + t67 * t140 + t66 * t141 + (t193 + t211) * t213, -t153 * t77 - t233 * t96, -t152 * t77 + t153 * t76 + t233 * t95 + t242 * t96, t203 * t96 - t224 * t77 + (-qJD(1) * t153 + t233) * t269, t152 * t76 - t242 * t95, -t203 * t95 + t224 * t76 + (-qJD(1) * t152 - t242) * t269, t139, -t24 * t203 - t12 * t224 + t135 * t242 + t187 * t76 + t113 * t152 - t127 * t95 + (qJD(1) * t63 + t36) * t269, -t113 * t153 + t135 * t233 - t127 * t96 + t13 * t224 + t187 * t77 + t25 * t203 + (-qJD(1) * t64 - t37) * t269, t12 * t153 - t13 * t152 - t233 * t24 - t242 * t25 + t36 * t96 + t37 * t95 - t63 * t77 - t64 * t76, t113 * t187 + t12 * t63 + t127 * t135 + t13 * t64 + t24 * t36 + t25 * t37, -t20 * t94 - t29 * t330, t20 * t93 - t21 * t94 + t29 * t58 - t30 * t330, t29 * t196 + t20 * t224 + (qJD(1) * t94 + t330) * t269, t21 * t93 + t30 * t58, t30 * t196 + t21 * t224 + (-qJD(1) * t93 - t58) * t269, (-t196 - t271) * t269, t119 * t21 - t4 * t196 - t2 * t224 + t68 * t30 + t45 * t93 + t65 * t58 + (qJD(1) * t22 + t5) * t269, -t228 * t224 - t119 * t20 + t3 * t196 - t68 * t29 + t45 * t94 + t65 * t330 + (-qJD(1) * t23 - t6) * t269, -t2 * t94 + t20 * t22 - t21 * t23 + t228 * t93 + t29 * t5 - t3 * t58 - t30 * t6 - t330 * t4, t119 * t45 + t2 * t22 - t228 * t23 + t3 * t6 + t4 * t5 + t65 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t257, t273 * t226, 0, t257, 0, 0, t226 * pkin(1) * t222, pkin(1) * t281, 0, 0, -t184 * t286 - t295, (-t129 + t293) * t223 + (-t130 + t290) * t221, -t250 + (t203 * t282 + t222 * t240) * qJD(1), -t182 * t287 - t294, t203 * t267 + (-t203 * t284 + t222 * t241) * qJD(1), t249, -pkin(2) * t130 + t131 * t203 + (pkin(7) * t286 + t289) * qJD(3) + ((-pkin(7) * t270 - t121) * t222 + (-pkin(6) * t241 - t289) * t224) * qJD(1), pkin(2) * t129 - t132 * t203 + (-pkin(7) * t287 + t288) * qJD(3) + ((-pkin(7) * t263 + t122) * t222 + (pkin(6) * t240 - t288) * t224) * qJD(1), t131 * t184 + t132 * t182 + ((qJD(3) * t184 - t130) * pkin(7) + t323) * t223 + ((qJD(3) * t182 - t129) * pkin(7) + t322) * t221, -t121 * t131 - t122 * t132 + (-t193 - t306) * t212 + (qJD(3) * t234 - t67 * t221 + t66 * t223) * pkin(7), t77 * t176 + t233 * t327, -t77 * t175 - t176 * t76 - t233 * t277 - t242 * t327, -t327 * t203 + (qJD(2) * t176 - t233) * t272, t76 * t175 + t242 * t277, t277 * t203 + (-qJD(2) * t175 + t242) * t272, t249, t210 * t76 + t113 * t175 + t307 * t203 + t277 * t127 + (qJD(2) * t125 - t36) * t272 + t237 * t242, t113 * t176 + t210 * t77 - t324 * t203 + t327 * t127 + t237 * t233 + (-qJD(2) * t126 + t37) * t272, -t12 * t176 - t125 * t77 - t126 * t76 - t13 * t175 + t233 * t307 + t242 * t324 - t277 * t37 - t327 * t36, t113 * t210 + t12 * t125 + t13 * t126 + t127 * t237 - t307 * t36 - t324 * t37, -t20 * t112 - t309 * t330, t20 * t111 - t112 * t21 - t308 * t330 + t309 * t58, t309 * t196 + (qJD(2) * t112 - t330) * t272, t21 * t111 + t308 * t58, t308 * t196 + (-qJD(2) * t111 + t58) * t272, t196 * t272, t45 * t111 + t142 * t21 + t308 * t68 + t278 * t58 - t311 * t196 + (qJD(2) * t46 - t5) * t272, t45 * t112 - t142 * t20 - t309 * t68 + t278 * t330 + t312 * t196 + (-qJD(2) * t47 + t6) * t272, t111 * t228 - t2 * t112 + t46 * t20 - t47 * t21 - t308 * t6 + t309 * t5 - t311 * t330 - t312 * t58, t45 * t142 + t2 * t46 - t228 * t47 + t278 * t68 + t311 * t5 + t312 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t291, -t182 ^ 2 + t184 ^ 2, -t129 - t293, -t291, -t130 - t290, t208, -t193 * t184 - t322, t182 * t193 - t323, 0, 0, t331, -t235 + t298, t77 - t299, -t331, -t76 - t300, t208, -t127 * t233 + t42 * t203 + (-t184 * t242 + t208 * t219) * pkin(3) + t12, t127 * t242 - t43 * t203 + (-t184 * t233 - t208 * t218) * pkin(3) - t13, (-t218 * t76 - t219 * t77) * pkin(3) + (t37 + t42) * t233 + (t43 - t36) * t242, -t36 * t42 - t37 * t43 + (t12 * t219 - t127 * t184 + t13 * t218) * pkin(3), t314, t338, t339, -t314, t334, t208, t158 * t208 + t196 * t302 - t79 * t58 + t333, -t159 * t208 - t196 * t301 - t330 * t79 + t337, t158 * t20 - t159 * t21 + (t302 + t6) * t330 + (t301 - t5) * t58, t2 * t158 - t159 * t228 - t301 * t6 - t302 * t5 - t68 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76 - t300, t77 + t299, -t235 - t298, t233 * t36 + t242 * t37 + t113, 0, 0, 0, 0, 0, 0, t21 - t303, -t20 + t305, -t315 - t316, t330 * t5 + t58 * t6 + t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t314, t338, t339, -t314, t334, t208, -t6 * t196 + t333, -t196 * t5 + t337, 0, 0;];
tauc_reg = t1;
