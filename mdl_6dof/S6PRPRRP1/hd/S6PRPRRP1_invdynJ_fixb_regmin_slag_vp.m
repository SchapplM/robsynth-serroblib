% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPRRP1
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
% tau_reg [6x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPRRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:24:18
% EndTime: 2021-01-16 01:24:34
% DurationCPUTime: 4.64s
% Computational Cost: add. (3724->432), mult. (7918->585), div. (0->0), fcn. (6245->18), ass. (0->239)
t183 = sin(qJ(4));
t186 = cos(qJ(4));
t214 = pkin(4) * t183 - pkin(9) * t186;
t129 = t214 * qJD(4);
t178 = cos(pkin(11));
t184 = sin(qJ(2));
t177 = sin(pkin(6));
t255 = qJD(1) * t177;
t234 = t184 * t255;
t136 = t178 * t234;
t175 = sin(pkin(11));
t187 = cos(qJ(2));
t228 = t187 * t255;
t87 = t175 * t228 + t136;
t326 = t129 - t87;
t159 = t175 * pkin(2) + pkin(8);
t185 = cos(qJ(5));
t182 = sin(qJ(5));
t244 = t182 * qJD(4);
t265 = t182 * t186;
t135 = t175 * t234;
t89 = t178 * t228 - t135;
t325 = t183 * t159 * t244 + t326 * t185 + t89 * t265;
t208 = -t186 * pkin(4) - t183 * pkin(9) - pkin(3);
t298 = t178 * pkin(2);
t115 = t208 - t298;
t245 = qJD(5) * t185;
t261 = t185 * t186;
t324 = t115 * t245 + t326 * t182 - t89 * t261;
t254 = qJD(2) * t177;
t226 = qJD(1) * t254;
t239 = qJDD(1) * t177;
t323 = t184 * t239 + t187 * t226;
t246 = qJD(5) * t183;
t322 = qJD(2) * t246 - qJDD(4);
t180 = cos(pkin(6));
t158 = t180 * t186;
t172 = qJ(2) + pkin(11);
t166 = sin(t172);
t272 = t177 * t183;
t100 = t166 * t272 - t158;
t179 = cos(pkin(10));
t271 = t177 * t186;
t167 = cos(t172);
t176 = sin(pkin(10));
t146 = t176 * t167;
t269 = t179 * t180;
t95 = t166 * t269 + t146;
t65 = t179 * t271 + t95 * t183;
t147 = t179 * t167;
t273 = t176 * t180;
t92 = t166 * t273 - t147;
t67 = t176 * t271 + t92 * t183;
t201 = -g(1) * t67 + g(2) * t65 + g(3) * t100;
t238 = t183 * qJDD(2);
t241 = t186 * qJD(2);
t63 = ((qJD(5) + t241) * qJD(4) + t238) * t182 + t322 * t185;
t127 = t159 * t261;
t206 = pkin(5) * t183 - qJ(6) * t261;
t242 = t185 * qJD(6);
t279 = qJ(6) * t183;
t321 = -t183 * t242 + t206 * qJD(4) + (-t127 + (-t115 + t279) * t182) * qJD(5) + t325;
t251 = qJD(4) * t159;
t263 = t183 * t185;
t294 = (-qJ(6) * qJD(5) - t251) * t263 + (-qJD(6) * t183 + (-qJ(6) * qJD(4) - qJD(5) * t159) * t186) * t182 + t324;
t149 = t180 * qJDD(1) + qJDD(3);
t152 = t180 * qJD(1) + qJD(3);
t249 = qJD(4) * t186;
t250 = qJD(4) * t183;
t148 = t187 * t239;
t91 = qJDD(2) * pkin(2) - t184 * t226 + t148;
t46 = t175 * t91 + t323 * t178;
t44 = qJDD(2) * pkin(8) + t46;
t131 = qJD(2) * pkin(2) + t228;
t84 = t175 * t131 + t136;
t77 = qJD(2) * pkin(8) + t84;
t207 = -t186 * t149 + t152 * t250 + t183 * t44 + t77 * t249;
t11 = -qJDD(4) * pkin(4) + t207;
t5 = t63 * pkin(5) + qJDD(6) + t11;
t320 = t201 - t5;
t243 = t185 * qJD(4);
t253 = qJD(2) * t183;
t121 = t182 * t253 - t243;
t153 = -qJD(5) + t241;
t277 = t121 * t153;
t227 = t186 * t243;
t62 = -qJD(2) * t227 - qJD(5) * t243 + t322 * t182 - t185 * t238;
t319 = -t62 + t277;
t123 = t185 * t253 + t244;
t276 = t123 * t153;
t318 = t63 - t276;
t317 = t186 * t152 - t183 * t77;
t98 = (t175 * t187 + t178 * t184) * t177;
t316 = -t183 * t98 + t158;
t268 = t180 * t183;
t101 = t166 * t271 + t268;
t164 = pkin(6) + t172;
t154 = sin(t164);
t165 = pkin(6) - t172;
t155 = sin(t165);
t112 = -t155 / 0.2e1 - t154 / 0.2e1;
t64 = -t176 * t272 + t92 * t186;
t66 = -t179 * t272 + t95 * t186;
t156 = cos(t164);
t157 = cos(t165);
t117 = t156 + t157;
t274 = t176 * t166;
t305 = -t179 / 0.2e1;
t79 = t117 * t305 + t274;
t270 = t179 * t166;
t306 = t176 / 0.2e1;
t81 = t117 * t306 + t270;
t315 = -g(3) * (-t101 * t182 + t112 * t185) - g(2) * (-t66 * t182 + t79 * t185) - g(1) * (t182 * t64 + t81 * t185);
t198 = -t182 * t246 + t227;
t170 = t186 * qJDD(2);
t240 = qJD(2) * qJD(4);
t223 = t183 * t240;
t120 = qJDD(5) - t170 + t223;
t262 = t185 * t120;
t314 = -t153 * t198 + t183 * t262;
t312 = t123 ^ 2;
t52 = t183 * t152 + t186 * t77;
t49 = qJD(4) * pkin(9) + t52;
t83 = t178 * t131 - t135;
t56 = qJD(2) * t208 - t83;
t16 = -t182 * t49 + t185 * t56;
t12 = -t123 * qJ(6) + t16;
t9 = -t153 * pkin(5) + t12;
t304 = t12 - t9;
t303 = pkin(5) * t182;
t301 = g(3) * t177;
t300 = t120 * pkin(5);
t299 = t121 * pkin(5);
t297 = t5 * t182;
t296 = qJ(6) + pkin(9);
t293 = -t121 * t227 - t63 * t263;
t292 = t11 * t182;
t17 = t182 * t56 + t185 * t49;
t13 = -t121 * qJ(6) + t17;
t291 = t13 * t153;
t289 = t62 * qJ(6);
t288 = t63 * qJ(6);
t116 = -t154 + t155;
t80 = t116 * t305 + t146;
t287 = t80 * t182;
t82 = t116 * t306 + t147;
t286 = t82 * t182;
t285 = t89 * t121;
t284 = t89 * t123;
t221 = qJD(5) * t296;
t233 = t182 * t241;
t128 = t214 * qJD(2);
t281 = -t182 * t128 - t185 * t317;
t283 = qJ(6) * t233 - t182 * t221 + t242 + t281;
t110 = t185 * t128;
t282 = -qJD(2) * t206 - t185 * t221 - t110 + (-qJD(6) + t317) * t182;
t113 = t157 / 0.2e1 - t156 / 0.2e1;
t278 = t113 * t182;
t275 = t153 * t185;
t267 = t180 * t184;
t266 = t180 * t187;
t264 = t183 * t149;
t259 = qJDD(1) - g(3);
t257 = t182 * t115 + t127;
t173 = t183 ^ 2;
t256 = -t186 ^ 2 + t173;
t252 = qJD(4) * t121;
t248 = qJD(5) * t121;
t247 = qJD(5) * t182;
t237 = t201 * t182;
t236 = t167 * t271;
t232 = t153 * t244;
t231 = t123 * t249;
t230 = t153 * t247;
t229 = t183 * t245;
t10 = qJDD(4) * pkin(9) + qJD(4) * t317 + t186 * t44 + t264;
t45 = -t323 * t175 + t178 * t91;
t23 = qJD(2) * t129 + qJDD(2) * t208 - t45;
t222 = t185 * t10 + t182 * t23 + t56 * t245 - t49 * t247;
t220 = -qJD(6) - t299;
t219 = t123 * t250 + t62 * t186;
t217 = t123 * t229;
t213 = t13 * t185 - t182 * t9;
t212 = -t13 * t182 - t185 * t9;
t75 = t186 * t98 + t268;
t209 = t175 * t184 - t178 * t187;
t97 = t209 * t177;
t37 = t182 * t97 + t185 * t75;
t36 = -t182 * t75 + t185 * t97;
t162 = t185 * pkin(5) + pkin(4);
t210 = t186 * t162 + t183 * t296;
t48 = -qJD(4) * pkin(4) - t317;
t204 = -t182 * t120 + t153 * t245;
t93 = t167 * t273 + t270;
t94 = t167 * t269 - t274;
t203 = -g(1) * (t82 * t185 + t93 * t265) - g(2) * (t80 * t185 - t94 * t265) - g(3) * (t113 * t185 - t182 * t236);
t202 = -g(1) * (-t93 * t261 + t286) - g(2) * (t94 * t261 + t287) - g(3) * (t185 * t236 + t278);
t200 = -g(1) * t64 + g(2) * t66 + g(3) * t101;
t199 = -g(3) * t180 + (-g(1) * t176 + g(2) * t179) * t177;
t197 = g(1) * t93 - g(2) * t94 - t167 * t301;
t196 = -pkin(9) * t120 - t153 * t48;
t160 = -pkin(3) - t298;
t76 = -qJD(2) * pkin(3) - t83;
t195 = -qJDD(4) * t159 + (qJD(2) * t160 + t76 + t89) * qJD(4);
t194 = -g(1) * (-t81 * t182 + t185 * t64) - g(2) * (-t79 * t182 - t66 * t185) - g(3) * (-t101 * t185 - t112 * t182) - t222;
t21 = t185 * t23;
t193 = -t17 * qJD(5) - t182 * t10 + t21;
t188 = qJD(4) ^ 2;
t192 = -qJD(2) * t87 + t159 * t188 - t197 - t45 + (-pkin(3) + t160) * qJDD(2);
t191 = t193 + t315;
t190 = -g(1) * (-t176 * t266 - t179 * t184) - g(2) * (-t176 * t184 + t179 * t266) - t187 * t301;
t189 = qJD(2) ^ 2;
t142 = t296 * t185;
t141 = t296 * t182;
t140 = -t175 * pkin(3) + t178 * pkin(8);
t139 = qJDD(4) * t186 - t188 * t183;
t138 = qJDD(4) * t183 + t188 * t186;
t134 = t178 * pkin(3) + t175 * pkin(8) + pkin(2);
t119 = t121 ^ 2;
t107 = (t159 + t303) * t183;
t103 = t185 * t115;
t90 = t209 * t254;
t88 = qJD(2) * t98;
t78 = t159 * t249 + (t186 * t244 + t229) * pkin(5);
t61 = -t182 * t279 + t257;
t55 = -qJ(6) * t263 + t103 + (-t159 * t182 - pkin(5)) * t186;
t41 = pkin(5) * t233 + t52;
t35 = t75 * qJD(4) - t183 * t90;
t34 = t316 * qJD(4) - t186 * t90;
t29 = -t220 + t48;
t15 = t219 - t314;
t14 = (-t63 + t232) * t186 + (t204 + t252) * t183;
t7 = -t37 * qJD(5) - t182 * t34 + t185 * t88;
t6 = t36 * qJD(5) + t182 * t88 + t185 * t34;
t4 = -t121 * qJD(6) + t222 - t288;
t3 = -t37 * t120 + t35 * t123 + t6 * t153 + t316 * t62;
t2 = t36 * t120 + t35 * t121 - t7 * t153 - t316 * t63;
t1 = -t123 * qJD(6) + t193 + t289 + t300;
t8 = [t259, 0, (qJDD(2) * t187 - t184 * t189) * t177, (-qJDD(2) * t184 - t187 * t189) * t177, t149 * t180 - t45 * t97 + t46 * t98 - t83 * t88 - t84 * t90 - g(3), 0, 0, 0, 0, 0, -t97 * t170 - t35 * qJD(4) + t316 * qJDD(4) + (-t186 * t88 + t97 * t250) * qJD(2), t97 * t238 - t34 * qJD(4) - t75 * qJDD(4) + (t183 * t88 + t97 * t249) * qJD(2), 0, 0, 0, 0, 0, t2, t3, t2, t3, -t6 * t121 - t7 * t123 + t36 * t62 - t37 * t63, t1 * t36 + t13 * t6 + t29 * t35 - t316 * t5 + t37 * t4 + t7 * t9 - g(3); 0, qJDD(2), t148 + t190, -g(1) * (t176 * t267 - t179 * t187) - g(2) * (-t176 * t187 - t179 * t267) - t259 * t184 * t177, t83 * t87 - t84 * t89 + (t46 * t175 + t45 * t178 + t190) * pkin(2), t173 * qJDD(2) + 0.2e1 * t186 * t223, 0.2e1 * t183 * t170 - 0.2e1 * t256 * t240, t138, t139, 0, t183 * t195 - t186 * t192, t183 * t192 + t186 * t195, t123 * t198 - t62 * t263, -t217 + (-t231 + (t62 + t248) * t183) * t182 + t293, t219 + t314, (t63 + t232) * t186 + (t204 - t252) * t183, -t120 * t186 - t153 * t250, t103 * t120 + (t115 * t247 - t325) * t153 + (t121 * t251 - t21 + (t153 * t159 + t49) * t245 + (qJD(4) * t48 + qJD(5) * t56 - t120 * t159 + t10) * t182) * t186 + (t16 * qJD(4) + t159 * t63 + t48 * t245 - t285 + t292) * t183 + t202, -t257 * t120 + t324 * t153 + (-t159 * t230 + (t123 * t159 + t48 * t185) * qJD(4) + t222) * t186 + (-t48 * t247 + t11 * t185 - t284 - t159 * t62 + (-t159 * t275 - t17) * qJD(4)) * t183 + t203, t107 * t63 + t55 * t120 + t78 * t121 + (t29 * t244 - t1) * t186 - t321 * t153 + (qJD(4) * t9 + t29 * t245 - t285 + t297) * t183 + t202, -t107 * t62 - t61 * t120 + t78 * t123 + (t29 * t243 + t4) * t186 + t294 * t153 + (-qJD(4) * t13 + t5 * t185 - t29 * t247 - t284) * t183 + t203, t55 * t62 - t61 * t63 - t321 * t123 - t294 * t121 + t212 * t249 + (-t213 * qJD(5) - t1 * t185 - t182 * t4 + t197) * t183, t4 * t61 + t1 * t55 + t5 * t107 - g(1) * (pkin(5) * t286 - (t179 * t134 + t140 * t273) * t184 + (-t134 * t273 + t179 * t140) * t187 - t210 * t93) - g(2) * (pkin(5) * t287 - (t176 * t134 - t140 * t269) * t184 + (t134 * t269 + t176 * t140) * t187 + t210 * t94) - g(3) * (pkin(5) * t278 + (t134 * t187 + t140 * t184 + t167 * t210) * t177) + t321 * t9 + (-t183 * t89 + t78) * t29 + t294 * t13; 0, 0, 0, 0, t199 + t149, 0, 0, 0, 0, 0, t139, -t138, 0, 0, 0, 0, 0, t14, t15, t14, t15, t217 + (t231 + (-t62 + t248) * t183) * t182 + t293, (qJD(4) * t213 - t5) * t186 + (qJD(4) * t29 + t212 * qJD(5) - t1 * t182 + t4 * t185) * t183 + t199; 0, 0, 0, 0, 0, -t183 * t189 * t186, t256 * t189, t238, t170, qJDD(4), t52 * qJD(4) - t76 * t253 + t201 - t207, -t264 + (-qJD(2) * t76 - t44) * t186 + t200, -t123 * t275 - t62 * t182, -t318 * t182 + t319 * t185, (-t123 * t183 + t153 * t261) * qJD(2) - t204, t230 + t262 + (t121 * t183 - t153 * t265) * qJD(2), t153 * t253, -t16 * t253 - pkin(4) * t63 + t110 * t153 - t52 * t121 + (-t153 * t317 + t196) * t182 + (qJD(5) * pkin(9) * t153 - t11 + t201) * t185, t17 * t253 + pkin(4) * t62 + t292 - t52 * t123 + (-pkin(9) * t247 + t281) * t153 + t196 * t185 - t237, -t9 * t253 - t141 * t120 - t41 * t121 - t162 * t63 - t282 * t153 + (-t29 * t241 + (t29 + t299) * qJD(5)) * t182 + t320 * t185, -t142 * t120 - t41 * t123 + t162 * t62 + t297 + t283 * t153 + (t123 * t303 + t185 * t29) * qJD(5) + (t13 * t183 - t29 * t261) * qJD(2) - t237, -t141 * t62 - t142 * t63 - t282 * t123 - t283 * t121 + (t153 * t9 + t4) * t185 + (-t1 + t291) * t182 - t200, t4 * t142 - t1 * t141 - t5 * t162 - g(1) * (t67 * t162 - t296 * t64) - g(2) * (-t65 * t162 + t296 * t66) - g(3) * (-t100 * t162 + t101 * t296) + t282 * t9 + (pkin(5) * t247 - t41) * t29 + t283 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123 * t121, -t119 + t312, -t62 - t277, -t276 - t63, t120, -t48 * t123 - t17 * t153 + t191, t48 * t121 - t16 * t153 + t194, 0.2e1 * t300 + t289 - t291 + (t220 - t29) * t123 + t191, -t312 * pkin(5) + t288 - t12 * t153 + (qJD(6) + t29) * t121 + t194, t62 * pkin(5) + t304 * t121, -t304 * t13 + (-t29 * t123 + t1 + t315) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t318, t319, -t119 - t312, t13 * t121 + t9 * t123 - t320;];
tau_reg = t8;
