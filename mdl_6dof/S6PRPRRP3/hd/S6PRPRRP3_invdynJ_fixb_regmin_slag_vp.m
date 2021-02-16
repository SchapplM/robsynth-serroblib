% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tau_reg [6x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:40
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:38:30
% EndTime: 2021-01-16 01:38:48
% DurationCPUTime: 5.11s
% Computational Cost: add. (4996->451), mult. (11724->601), div. (0->0), fcn. (9645->14), ass. (0->235)
t336 = 2 * qJD(4);
t169 = sin(pkin(11));
t172 = cos(pkin(11));
t177 = sin(qJ(4));
t308 = cos(qJ(4));
t200 = -t177 * t169 + t308 * t172;
t171 = sin(pkin(6));
t180 = cos(qJ(2));
t262 = t171 * t180;
t188 = t200 * t262;
t105 = qJD(1) * t188;
t175 = qJ(3) + pkin(8);
t142 = t175 * t169;
t143 = t175 * t172;
t201 = -t308 * t142 - t177 * t143;
t52 = qJD(3) * t200 + qJD(4) * t201;
t335 = t105 - t52;
t133 = t200 * qJD(4);
t137 = t308 * t169 + t177 * t172;
t134 = t137 * qJD(4);
t178 = sin(qJ(2));
t254 = qJD(1) * t171;
t236 = t178 * t254;
t334 = pkin(4) * t134 - pkin(9) * t133 - t236;
t139 = qJD(2) * qJ(3) + t236;
t173 = cos(pkin(6));
t253 = qJD(1) * t173;
t152 = t172 * t253;
t299 = pkin(8) * qJD(2);
t97 = t152 + (-t139 - t299) * t169;
t114 = t172 * t139 + t169 * t253;
t98 = t172 * t299 + t114;
t40 = t177 * t97 + t308 * t98;
t333 = t40 * qJD(4);
t332 = t200 * qJD(2);
t176 = sin(qJ(5));
t179 = cos(qJ(5));
t331 = t105 * t176 + t334 * t179;
t249 = qJD(5) * t179;
t158 = t172 * pkin(3) + pkin(2);
t82 = -pkin(4) * t200 - pkin(9) * t137 - t158;
t330 = t334 * t176 - t335 * t179 + t82 * t249;
t103 = -t177 * t142 + t308 * t143;
t189 = t137 * t262;
t288 = -qJD(1) * t189 + qJD(3) * t137 + qJD(4) * t103;
t132 = t137 * qJD(2);
t251 = qJD(4) * t176;
t109 = t132 * t179 + t251;
t280 = t109 * t176;
t190 = t137 * qJDD(2);
t183 = qJD(2) * t133 + t190;
t329 = -qJD(4) * qJD(5) - t183;
t168 = pkin(11) + qJ(4);
t161 = sin(t168);
t162 = cos(t168);
t263 = t171 * t178;
t118 = t161 * t263 - t173 * t162;
t287 = cos(pkin(10));
t224 = t287 * t178;
t170 = sin(pkin(10));
t264 = t170 * t180;
t128 = t173 * t224 + t264;
t228 = t171 * t287;
t87 = t128 * t161 + t162 * t228;
t223 = t287 * t180;
t265 = t170 * t178;
t126 = t173 * t265 - t223;
t268 = t170 * t171;
t89 = t126 * t161 + t162 * t268;
t195 = -g(1) * t89 + g(2) * t87 + g(3) * t118;
t211 = -qJ(6) * t133 - qJD(6) * t137;
t96 = t179 * t103;
t310 = pkin(5) * t134 - t176 * t52 + t211 * t179 + (-t96 + (qJ(6) * t137 - t82) * t176) * qJD(5) + t331;
t233 = t137 * t249;
t328 = -qJ(6) * t233 + (-qJD(5) * t103 + t211) * t176 + t330;
t327 = g(3) * t171;
t235 = t180 * t254;
t245 = qJDD(2) * qJ(3);
t247 = qJDD(1) * t171;
t111 = t178 * t247 + t245 + (qJD(3) + t235) * qJD(2);
t246 = qJDD(1) * t173;
t150 = t172 * t246;
t59 = t150 + (-pkin(8) * qJDD(2) - t111) * t169;
t244 = t172 * qJDD(2);
t77 = t172 * t111 + t169 * t246;
t60 = pkin(8) * t244 + t77;
t208 = t177 * t59 + t308 * t60;
t324 = -t177 * t98 + t308 * t97;
t10 = qJDD(4) * pkin(9) + qJD(4) * t324 + t208;
t36 = qJD(4) * pkin(9) + t40;
t215 = qJD(3) - t235;
t120 = -qJD(2) * t158 + t215;
t47 = -pkin(4) * t332 - pkin(9) * t132 + t120;
t20 = t176 * t47 + t179 * t36;
t231 = t180 * t247;
t252 = qJD(2) * t178;
t232 = qJD(1) * t252;
t321 = t171 * t232 + qJDD(3);
t204 = -t231 + t321;
t286 = qJDD(2) * pkin(2);
t117 = t204 - t286;
t214 = t200 * qJDD(2);
t85 = qJD(2) * t134 - t214;
t25 = -pkin(3) * t244 + t85 * pkin(4) - pkin(9) * t183 + t117;
t23 = t179 * t25;
t185 = -qJD(5) * t20 - t176 * t10 + t23;
t250 = qJD(5) * t176;
t33 = -t176 * qJDD(4) + t132 * t250 + t329 * t179;
t298 = qJ(6) * t33;
t79 = qJDD(5) + t85;
t317 = pkin(5) * t79;
t1 = -qJD(6) * t109 + t185 + t298 + t317;
t121 = qJD(5) - t332;
t107 = -t179 * qJD(4) + t132 * t176;
t13 = -qJ(6) * t107 + t20;
t296 = t121 * t13;
t326 = t1 + t296;
t205 = t177 * t60 - t308 * t59 + t333;
t11 = -qJDD(4) * pkin(4) + t205;
t216 = -t179 * qJDD(4) + t132 * t249;
t34 = t176 * t190 + (qJD(5) + t332) * t251 + t216;
t5 = t34 * pkin(5) + qJDD(6) + t11;
t325 = t195 - t5;
t261 = t176 * t133;
t199 = t233 + t261;
t289 = pkin(5) * t199 + t288;
t323 = t121 * t280;
t127 = -t173 * t223 + t265;
t129 = t173 * t264 + t224;
t218 = g(1) * t129 + g(2) * t127;
t192 = -g(3) * t262 + t218;
t322 = t192 * t161;
t119 = t161 * t173 + t162 * t263;
t258 = t179 * t180;
t239 = t171 * t258;
t86 = t126 * t162 - t161 * t268;
t88 = t128 * t162 - t161 * t228;
t320 = -g(3) * (-t119 * t176 - t239) - g(2) * (t127 * t179 - t176 * t88) - g(1) * (t129 * t179 + t176 * t86);
t213 = (-t139 * t169 + t152) * t169 - t114 * t172;
t319 = t180 * t213 - (-qJD(2) * pkin(2) + t215) * t178;
t318 = t109 ^ 2;
t19 = -t176 * t36 + t179 * t47;
t12 = -qJ(6) * t109 + t19;
t8 = pkin(5) * t121 + t12;
t311 = t12 - t8;
t304 = t107 * pkin(5);
t303 = t176 * t5;
t174 = -qJ(6) - pkin(9);
t302 = -t107 * t249 - t176 * t34;
t80 = pkin(4) * t132 - pkin(9) * t332;
t301 = t176 * t80 + t179 * t324;
t300 = t176 * t82 + t96;
t297 = qJ(6) * t34;
t295 = t176 * t33;
t294 = t176 * t79;
t76 = -t111 * t169 + t150;
t293 = t76 * t169;
t292 = t77 * t172;
t229 = qJD(5) * t174;
t276 = t332 * t176;
t291 = qJ(6) * t276 + qJD(6) * t179 + t176 * t229 - t301;
t275 = t332 * t179;
t69 = t179 * t80;
t290 = -pkin(5) * t132 + qJ(6) * t275 + t179 * t229 - t69 + (-qJD(6) + t324) * t176;
t285 = t107 * t121;
t284 = t107 * t332;
t283 = t107 * t132;
t282 = t109 * t121;
t281 = t109 * t132;
t278 = t126 * t176;
t277 = t128 * t176;
t274 = t133 * t179;
t273 = t137 * t176;
t272 = t137 * t179;
t270 = t162 * t176;
t269 = t162 * t179;
t267 = t170 * t173;
t266 = t170 * t175;
t260 = t176 * t180;
t181 = qJD(2) ^ 2;
t257 = t180 * t181;
t256 = qJDD(1) - g(3);
t255 = t169 ^ 2 + t172 ^ 2;
t242 = t195 * t176;
t240 = t171 * t260;
t237 = t287 * pkin(2);
t234 = t171 * t252;
t230 = t179 * t10 + t176 * t25 + t47 * t249 - t36 * t250;
t227 = t287 * qJ(3);
t226 = t287 * t158;
t225 = t287 * t175;
t222 = qJD(6) + t304;
t221 = t121 * t179;
t219 = -g(1) * t126 + g(2) * t128;
t2 = -qJD(6) * t107 + t230 - t297;
t217 = -t121 * t8 + t2;
t160 = pkin(5) * t179 + pkin(4);
t212 = -t160 * t162 + t161 * t174;
t210 = t179 * t79 + (-t250 + t276) * t121;
t124 = -t169 * t263 + t172 * t173;
t125 = t169 * t173 + t172 * t263;
t64 = t177 * t124 + t308 * t125;
t48 = -t176 * t64 - t239;
t206 = -t179 * t64 + t240;
t35 = -qJD(4) * pkin(4) - t324;
t203 = -pkin(9) * t79 + t121 * t35;
t202 = t308 * t124 - t177 * t125;
t198 = -t137 * t250 + t274;
t197 = -g(1) * (-t126 * t179 + t129 * t270) - g(2) * (t127 * t270 + t128 * t179) - (-t162 * t260 + t178 * t179) * t327;
t196 = -g(1) * (-t129 * t269 - t278) - g(2) * (-t127 * t269 + t277) - (t162 * t258 + t176 * t178) * t327;
t194 = -g(1) * t86 + g(2) * t88 + g(3) * t119;
t193 = (qJDD(2) * t180 - t178 * t181) * t171;
t187 = t192 + t231;
t186 = -g(1) * (-t129 * t176 + t179 * t86) - g(2) * (-t127 * t176 - t179 * t88) - g(3) * (-t119 * t179 + t240) - t230;
t184 = t185 + t320;
t182 = t329 * t176 - t216;
t145 = t174 * t179;
t144 = t174 * t176;
t106 = t107 ^ 2;
t100 = -qJDD(2) * t158 + t204;
t72 = t179 * t82;
t54 = pkin(5) * t273 - t201;
t38 = qJD(2) * t189 + qJD(4) * t64;
t37 = qJD(2) * t188 + qJD(4) * t202;
t28 = pkin(5) * t276 + t40;
t27 = -qJ(6) * t273 + t300;
t26 = t222 + t35;
t24 = -pkin(5) * t200 - qJ(6) * t272 - t103 * t176 + t72;
t18 = qJD(5) * t206 - t176 * t37 + t179 * t234;
t17 = qJD(5) * t48 + t176 * t234 + t179 * t37;
t16 = -t121 ^ 2 * t179 - t281 - t294;
t15 = t210 - t283;
t4 = t107 * t38 + t121 * t18 - t202 * t34 + t48 * t79;
t3 = t109 * t38 - t121 * t17 + t202 * t33 + t206 * t79;
t6 = [t256, 0, t193, (-qJDD(2) * t178 - t257) * t171, t172 * t193, t255 * t171 * t257 + (-t124 * t169 + t125 * t172) * qJDD(2), t124 * t76 + t125 * t77 - g(3) + (-t319 * qJD(2) - t117 * t180) * t171, 0, 0, 0, 0, 0, -qJD(4) * t38 + qJDD(4) * t202 + (-t180 * t85 - t252 * t332) * t171, -t37 * qJD(4) - t64 * qJDD(4) + (-t180 * t190 + (t178 * t132 - t133 * t180) * qJD(2)) * t171, 0, 0, 0, 0, 0, t4, t3, t4, t3, -t107 * t17 - t109 * t18 + t206 * t34 + t33 * t48, t1 * t48 + t13 * t17 + t18 * t8 - t2 * t206 - t202 * t5 + t26 * t38 - g(3); 0, qJDD(2), t187, -t256 * t263 + t219, (t286 - t117 + (-g(3) * t180 + t232) * t171 + t218) * t172, -g(3) * t263 - t219 + t292 - t293 + (t215 * qJD(2) + t245) * t255, qJ(3) * t292 - qJ(3) * t293 - t117 * pkin(2) - g(1) * (-(qJ(3) * t267 + t237) * t178 + (-pkin(2) * t267 + t227) * t180) - g(2) * (-(t170 * pkin(2) - t173 * t227) * t178 + (t170 * qJ(3) + t173 * t237) * t180) - t213 * qJD(3) + (-g(3) * (pkin(2) * t180 + qJ(3) * t178) + t319 * qJD(1)) * t171, t132 * t133 + t137 * t183, -t132 * t134 + t133 * t332 - t137 * t85 + t183 * t200, qJD(4) * t133 + qJDD(4) * t137, -qJD(4) * t134 + qJDD(4) * t200, 0, -qJD(4) * t288 + qJDD(4) * t201 - t100 * t200 + t120 * t134 - t158 * t85 + t162 * t192 + t236 * t332, qJD(4) * t335 - t103 * qJDD(4) + t100 * t137 + t120 * t133 - t132 * t236 - t158 * t183 - t322, t109 * t198 - t272 * t33, (-t107 * t179 - t280) * t133 + (t295 - t179 * t34 + (t107 * t176 - t109 * t179) * qJD(5)) * t137, t109 * t134 + t121 * t198 + t200 * t33 + t272 * t79, -t107 * t134 - t121 * t199 + t200 * t34 - t273 * t79, t121 * t134 - t200 * t79, t72 * t79 - (-t249 * t36 + t23) * t200 + t19 * t134 - t201 * t34 + t35 * t233 + (-t103 * t249 + t331) * t121 + t288 * t107 + ((-qJD(5) * t82 - t52) * t121 - t103 * t79 - (-qJD(5) * t47 - t10) * t200 + t11 * t137 + t35 * t133) * t176 + t196, -t300 * t79 + t230 * t200 - t20 * t134 + t201 * t33 + t35 * t274 + (t11 * t179 - t250 * t35) * t137 + (t103 * t250 - t330) * t121 + t288 * t109 + t197, t26 * t261 - t1 * t200 + t134 * t8 + t24 * t79 + t34 * t54 + (t249 * t26 + t303) * t137 + t310 * t121 + t289 * t107 + t196, t26 * t274 - t13 * t134 + t200 * t2 - t27 * t79 - t33 * t54 + (t179 * t5 - t250 * t26) * t137 - t328 * t121 + t289 * t109 + t197, t24 * t33 - t27 * t34 + (-t13 * t176 - t179 * t8) * t133 - t310 * t109 - t328 * t107 + t322 + (-t1 * t179 - t176 * t2 + (-t13 * t179 + t176 * t8) * qJD(5)) * t137, t2 * t27 + t1 * t24 + t5 * t54 - g(1) * (-pkin(5) * t278 - (t173 * t266 + t226) * t178 + (-t158 * t267 + t225) * t180 + t212 * t129) - g(2) * (pkin(5) * t277 - (t170 * t158 - t173 * t225) * t178 + (t173 * t226 + t266) * t180 + t212 * t127) - ((pkin(5) * t176 + t175) * t178 + (t158 - t212) * t180) * t327 + t310 * t8 + t289 * t26 + t328 * t13; 0, 0, 0, 0, -t244, -t255 * t181, qJD(2) * t213 - t187 - t286 + t321, 0, 0, 0, 0, 0, t132 * t336 - t214, t332 * t336 + t190, 0, 0, 0, 0, 0, t15, t16, t15, t16, (t33 + t284) * t179 + t323 + t302, -t132 * t26 + t217 * t176 + t326 * t179 - t192; 0, 0, 0, 0, 0, 0, 0, -t132 * t332, t132 ^ 2 - t332 ^ 2, t190, t214, qJDD(4), -t120 * t132 + t195 - t205 + t333, -t120 * t332 + t194 - t208, t109 * t221 - t295, (-t33 + t284) * t179 - t323 + t302, t121 * t221 - t281 + t294, t210 + t283, -t121 * t132, -pkin(4) * t34 - t40 * t107 - t69 * t121 - t19 * t132 + (t121 * t324 + t203) * t176 + (-pkin(9) * qJD(5) * t121 - t11 + t195) * t179, pkin(4) * t33 - t40 * t109 + t11 * t176 + t20 * t132 + (pkin(9) * t250 + t301) * t121 + t203 * t179 - t242, -t107 * t28 - t132 * t8 + t144 * t79 - t160 * t34 + t290 * t121 + (-t332 * t26 + (t26 + t304) * qJD(5)) * t176 + t325 * t179, -t26 * t275 - t109 * t28 + t13 * t132 + t145 * t79 + t160 * t33 + t303 - t291 * t121 + (pkin(5) * t280 + t179 * t26) * qJD(5) - t242, -t291 * t107 - t290 * t109 + t144 * t33 + t145 * t34 - t326 * t176 + t217 * t179 - t194, -t2 * t145 + t1 * t144 - t5 * t160 - g(1) * (t160 * t89 + t174 * t86) - g(2) * (-t160 * t87 - t174 * t88) - g(3) * (-t118 * t160 - t119 * t174) + t290 * t8 + (pkin(5) * t250 - t28) * t26 + t291 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109 * t107, -t106 + t318, -t33 + t285, t182 + t282, t79, -t109 * t35 + t121 * t20 + t184, t107 * t35 + t121 * t19 + t186, 0.2e1 * t317 + t298 + t296 + (-t222 - t26) * t109 + t184, -pkin(5) * t318 + t297 + t12 * t121 + (qJD(6) + t26) * t107 + t186, pkin(5) * t33 + t311 * t107, -t311 * t13 + (-t26 * t109 + t1 + t320) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182 + t282, -t33 - t285, -t106 - t318, t13 * t107 + t8 * t109 - t325;];
tau_reg = t6;
