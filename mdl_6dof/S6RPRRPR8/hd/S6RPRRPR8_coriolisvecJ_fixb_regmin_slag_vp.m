% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta5]';
% 
% Output:
% tauc_reg [6x29]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:24:55
% EndTime: 2019-03-09 05:25:04
% DurationCPUTime: 4.10s
% Computational Cost: add. (4589->381), mult. (10333->555), div. (0->0), fcn. (7152->8), ass. (0->205)
t186 = sin(qJ(4));
t189 = cos(qJ(4));
t240 = t189 * qJD(3);
t190 = cos(qJ(3));
t248 = qJD(1) * t190;
t156 = t186 * t248 - t240;
t229 = t189 * t248;
t247 = qJD(3) * t186;
t158 = t229 + t247;
t183 = sin(pkin(10));
t184 = cos(pkin(10));
t100 = t156 * t183 - t158 * t184;
t185 = sin(qJ(6));
t188 = cos(qJ(6));
t241 = qJD(6) * t185;
t206 = -t156 * t184 - t158 * t183;
t260 = t188 * t206;
t242 = qJD(4) * t190;
t225 = t186 * t242;
t187 = sin(qJ(3));
t227 = t187 * t240;
t196 = -t225 - t227;
t111 = qJD(1) * t196 + qJD(4) * t240;
t246 = qJD(3) * t187;
t228 = t186 * t246;
t112 = -qJD(1) * t228 + qJD(4) * t158;
t61 = -t111 * t183 - t112 * t184;
t62 = t111 * t184 - t112 * t183;
t10 = qJD(6) * t260 + t100 * t241 + t185 * t61 + t188 * t62;
t249 = qJD(1) * t187;
t175 = qJD(4) + t249;
t171 = qJD(6) + t175;
t47 = t100 * t185 + t260;
t280 = t171 * t47;
t302 = t10 - t280;
t290 = -t188 * t100 + t185 * t206;
t301 = t290 * t47;
t300 = t290 ^ 2 - t47 ^ 2;
t162 = pkin(3) * t187 - pkin(8) * t190 + qJ(2);
t141 = t162 * qJD(1);
t191 = -pkin(1) - pkin(7);
t173 = t191 * qJD(1) + qJD(2);
t161 = t187 * t173;
t144 = qJD(3) * pkin(8) + t161;
t94 = t141 * t186 + t144 * t189;
t68 = -qJ(5) * t156 + t94;
t279 = t184 * t68;
t93 = t189 * t141 - t144 * t186;
t67 = -qJ(5) * t158 + t93;
t55 = pkin(4) * t175 + t67;
t26 = t183 * t55 + t279;
t289 = pkin(9) * t206;
t16 = t26 + t289;
t15 = t16 * t241;
t238 = qJD(1) * qJD(3);
t178 = t190 * t238;
t212 = pkin(3) * t190 + pkin(8) * t187;
t150 = qJD(3) * t212 + qJD(2);
t130 = t150 * qJD(1);
t116 = t189 * t130;
t195 = -qJD(4) * t94 + t116;
t245 = qJD(3) * t190;
t23 = -t111 * qJ(5) - t158 * qJD(5) + (pkin(4) * qJD(1) - t173 * t186) * t245 + t195;
t226 = t190 * t240;
t243 = qJD(4) * t189;
t244 = qJD(4) * t186;
t198 = t186 * t130 + t141 * t243 - t144 * t244 + t173 * t226;
t28 = -qJ(5) * t112 - qJD(5) * t156 + t198;
t6 = -t183 * t28 + t184 * t23;
t4 = pkin(5) * t178 - pkin(9) * t62 + t6;
t269 = t173 * t190;
t145 = -qJD(3) * pkin(3) - t269;
t106 = pkin(4) * t156 + qJD(5) + t145;
t52 = -pkin(5) * t206 + t106;
t299 = -t185 * t4 - t52 * t47 + t15;
t11 = qJD(6) * t290 + t185 * t62 - t188 * t61;
t277 = t290 * t171;
t297 = -t11 + t277;
t285 = -qJ(5) - pkin(8);
t219 = qJD(4) * t285;
t230 = t186 * t249;
t239 = t189 * qJD(5);
t160 = t212 * qJD(1);
t259 = t189 * t190;
t253 = t186 * t160 + t173 * t259;
t296 = -qJ(5) * t230 + t186 * t219 + t239 - t253;
t143 = t189 * t160;
t264 = t186 * t190;
t234 = t173 * t264;
t262 = t187 * t189;
t295 = -t186 * qJD(5) + t189 * t219 + t234 - t143 - (pkin(4) * t190 + qJ(5) * t262) * qJD(1);
t7 = t183 * t23 + t184 * t28;
t5 = pkin(9) * t61 + t7;
t231 = -t185 * t5 + t188 * t4;
t294 = -t52 * t290 + t231;
t293 = pkin(9) * t100;
t152 = t183 * t189 + t184 * t186;
t137 = t152 * qJD(4);
t138 = t152 * qJD(1);
t292 = t187 * t138 + t137;
t265 = t184 * t189;
t205 = t183 * t186 - t265;
t287 = qJD(4) * t205;
t291 = -t183 * t230 + t249 * t265 - t287;
t237 = 0.2e1 * qJD(1);
t282 = -t296 * t183 + t295 * t184;
t281 = t295 * t183 + t296 * t184;
t261 = t187 * t191;
t252 = t186 * t162 + t189 * t261;
t286 = pkin(4) * t183;
t135 = t189 * t150;
t263 = t186 * t191;
t222 = pkin(4) - t263;
t39 = qJ(5) * t227 + t135 - t252 * qJD(4) + (qJ(5) * t244 + t222 * qJD(3) - t239) * t190;
t224 = t189 * t242;
t233 = t186 * t150 + t162 * t243 + t191 * t226;
t41 = -qJ(5) * t224 + (-qJD(5) * t190 + (qJ(5) * qJD(3) - qJD(4) * t191) * t187) * t186 + t233;
t13 = t183 * t39 + t184 * t41;
t207 = -t152 * t185 - t188 * t205;
t284 = qJD(6) * t207 - t292 * t185 + t291 * t188;
t99 = t152 * t188 - t185 * t205;
t283 = qJD(6) * t99 + t291 * t185 + t292 * t188;
t63 = t183 * t68;
t31 = t184 * t67 - t63;
t25 = t184 * t55 - t63;
t14 = pkin(5) * t175 + t25 + t293;
t278 = t188 * t14;
t107 = -qJ(5) * t264 + t252;
t149 = t189 * t162;
t97 = -qJ(5) * t259 + t222 * t187 + t149;
t50 = t184 * t107 + t183 * t97;
t276 = -t161 + t292 * pkin(5) + (t230 + t244) * pkin(4);
t126 = t152 * t190;
t275 = -t205 * qJD(1) + qJD(3) * t126 - t187 * t287;
t128 = t205 * t190;
t274 = -qJD(3) * t128 - t137 * t187 - t138;
t273 = t111 * t186;
t272 = t145 * t186;
t271 = t156 * t175;
t270 = t158 * t175;
t268 = t175 * t186;
t267 = t175 * t187;
t266 = t175 * t189;
t258 = t190 * t111;
t257 = t190 * t191;
t192 = qJD(3) ^ 2;
t256 = t192 * t187;
t255 = t192 * t190;
t193 = qJD(1) ^ 2;
t254 = t193 * qJ(2);
t166 = t285 * t186;
t167 = t285 * t189;
t109 = t183 * t166 - t184 * t167;
t182 = t190 ^ 2;
t251 = t187 ^ 2 - t182;
t250 = -t192 - t193;
t236 = qJD(2) * t237;
t235 = t186 * t261;
t232 = -pkin(4) * t189 - pkin(3);
t88 = pkin(4) * t112 + t173 * t246;
t221 = qJD(6) * t14 + t5;
t12 = -t183 * t41 + t184 * t39;
t30 = -t183 * t67 - t279;
t49 = -t107 * t183 + t184 * t97;
t108 = t184 * t166 + t167 * t183;
t218 = pkin(4) * t264 - t257;
t217 = t156 + t240;
t216 = -t158 + t247;
t215 = qJD(4) * t187 + qJD(1);
t87 = -pkin(9) * t205 + t109;
t214 = pkin(5) * t248 + t291 * pkin(9) + qJD(6) * t87 - t282;
t86 = -pkin(9) * t152 + t108;
t213 = t292 * pkin(9) - qJD(6) * t86 - t281;
t125 = t152 * t187;
t211 = qJD(6) * t125 - t274;
t127 = t205 * t187;
t210 = -qJD(6) * t127 + t275;
t2 = t185 * t14 + t188 * t16;
t33 = pkin(5) * t187 + pkin(9) * t128 + t49;
t34 = -pkin(9) * t126 + t50;
t209 = t185 * t33 + t188 * t34;
t208 = -t188 * t126 + t128 * t185;
t76 = -t126 * t185 - t128 * t188;
t204 = qJD(1) * t182 - t267;
t177 = pkin(4) * t184 + pkin(5);
t203 = t177 * t185 + t188 * t286;
t202 = t177 * t188 - t185 * t286;
t201 = -pkin(8) * t245 + t145 * t187;
t199 = t191 * t246 + (t224 - t228) * pkin(4);
t169 = t187 * t178;
t119 = pkin(5) * t205 + t232;
t105 = pkin(5) * t126 + t218;
t83 = t152 * t242 - t183 * t228 + t184 * t227;
t81 = t152 * t246 + t190 * t287;
t69 = pkin(4) * t158 - pkin(5) * t100;
t51 = -pkin(5) * t81 + t199;
t32 = -pkin(5) * t61 + t88;
t20 = qJD(6) * t76 - t185 * t83 - t188 * t81;
t19 = qJD(6) * t208 + t185 * t81 - t188 * t83;
t18 = t31 + t293;
t17 = t30 - t289;
t9 = pkin(9) * t81 + t13;
t8 = pkin(5) * t245 + pkin(9) * t83 + t12;
t1 = -t16 * t185 + t278;
t3 = [0, 0, 0, 0, t236, qJ(2) * t236, -0.2e1 * t169, 0.2e1 * t251 * t238, -t256, -t255, 0, -t191 * t256 + (qJ(2) * t245 + qJD(2) * t187) * t237, -t191 * t255 + (-qJ(2) * t246 + qJD(2) * t190) * t237, t158 * t196 + t189 * t258 (t156 * t189 + t158 * t186) * t246 + (-t273 - t112 * t189 + (t156 * t186 - t158 * t189) * qJD(4)) * t190, -t175 * t225 + t111 * t187 + (t158 * t190 + t189 * t204) * qJD(3), -t175 * t224 - t112 * t187 + (-t156 * t190 - t186 * t204) * qJD(3), t175 * t245 + t169, -t112 * t257 + t116 * t187 + t135 * t175 + (t145 * t259 - t175 * t252 - t187 * t94) * qJD(4) + ((t156 * t191 - t272) * t187 + (-t175 * t263 + (t149 - t235) * qJD(1) + t93) * t190) * qJD(3) -(-qJD(4) * t235 + t233) * t175 - t198 * t187 + (-t191 * t111 - t145 * t244) * t190 + ((-t252 * qJD(1) - t94) * t190 + (t191 * t158 + (-t145 + t269) * t189) * t187) * qJD(3), t100 * t12 - t126 * t7 + t128 * t6 + t13 * t206 + t25 * t83 + t26 * t81 - t49 * t62 + t50 * t61, t106 * t199 + t25 * t12 + t26 * t13 + t218 * t88 + t6 * t49 + t7 * t50, t10 * t76 + t19 * t290, t10 * t208 - t11 * t76 + t19 * t47 - t20 * t290, t10 * t187 + t19 * t171 + (qJD(1) * t76 + t290) * t245, -t11 * t187 - t20 * t171 + (qJD(1) * t208 + t47) * t245, t171 * t245 + t169 (-t185 * t9 + t188 * t8) * t171 + t231 * t187 - t51 * t47 + t105 * t11 - t32 * t208 + t52 * t20 + (-t171 * t209 - t187 * t2) * qJD(6) + ((-t185 * t34 + t188 * t33) * qJD(1) + t1) * t245, t105 * t10 + t15 * t187 + t52 * t19 + t32 * t76 + t51 * t290 + (-(-qJD(6) * t34 + t8) * t171 - t4 * t187) * t185 + (-(qJD(6) * t33 + t9) * t171 - t221 * t187) * t188 + (-qJD(1) * t209 - t2) * t245; 0, 0, 0, 0, -t193, -t254, 0, 0, 0, 0, 0, t250 * t187, t250 * t190, 0, 0, 0, 0, 0, -t190 * t112 - t215 * t266 + (t187 * t156 + (-t175 - t249) * t264) * qJD(3), -t258 + t215 * t268 + (-t175 * t259 + (t158 - t229) * t187) * qJD(3), -t100 * t275 + t125 * t62 - t127 * t61 + t206 * t274, t106 * t246 - t6 * t125 - t7 * t127 - t88 * t190 - t275 * t25 + t274 * t26, 0, 0, 0, 0, 0, -t190 * t11 + (t185 * t211 - t188 * t210) * t171 + ((-t125 * t188 + t127 * t185) * t248 - t187 * t47) * qJD(3), -t190 * t10 + (t185 * t210 + t188 * t211) * t171 + (-(-t125 * t185 - t127 * t188) * t248 + t187 * t290) * qJD(3); 0, 0, 0, 0, 0, 0, t190 * t193 * t187, -t251 * t193, 0, 0, 0, -t190 * t254, t187 * t254, t158 * t266 + t273 (t111 - t271) * t189 + (-t112 - t270) * t186, t175 * t243 + (t175 * t262 + t190 * t216) * qJD(1), -t175 * t244 + (-t186 * t267 + t190 * t217) * qJD(1), -t175 * t248, -pkin(3) * t112 - t143 * t175 + (t175 * t264 - t187 * t217) * t173 + (-pkin(8) * t266 + t272) * qJD(4) + (t186 * t201 - t93 * t190) * qJD(1), -pkin(3) * t111 + t253 * t175 + t216 * t161 + (pkin(8) * t268 + t145 * t189) * qJD(4) + (t189 * t201 + t94 * t190) * qJD(1), t282 * t100 - t108 * t62 + t109 * t61 - t6 * t152 - t7 * t205 + t281 * t206 - t291 * t25 - t292 * t26, t7 * t109 + t6 * t108 + t88 * t232 + t281 * t26 + t282 * t25 + (pkin(4) * t268 - t161) * t106, t10 * t99 + t284 * t290, t10 * t207 - t99 * t11 - t283 * t290 + t284 * t47, t284 * t171 + (qJD(3) * t99 - t290) * t248, -t283 * t171 + (qJD(3) * t207 - t47) * t248, -t171 * t248, t119 * t11 - t32 * t207 + t283 * t52 - t276 * t47 + (t185 * t213 - t188 * t214) * t171 + ((-t185 * t87 + t188 * t86) * qJD(3) - t1) * t248, t119 * t10 + t32 * t99 + t284 * t52 + t276 * t290 + (t185 * t214 + t188 * t213) * t171 + (-(t185 * t86 + t188 * t87) * qJD(3) + t2) * t248; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158 * t156, -t156 ^ 2 + t158 ^ 2, t111 + t271, -t112 + t270, t178, -qJD(3) * t234 - t145 * t158 + t94 * t175 + t195, t145 * t156 + t175 * t93 - t198 (t183 * t61 - t184 * t62) * pkin(4) + (-t31 + t25) * t206 + (-t26 - t30) * t100, -t25 * t30 - t26 * t31 + (-t106 * t158 + t183 * t7 + t184 * t6) * pkin(4), -t301, t300, t302, t297, t178, t202 * t178 - (t17 * t188 - t18 * t185) * t171 + t69 * t47 + (-t171 * t203 - t2) * qJD(6) + t294, -t203 * t178 - t188 * t5 + (t17 * t185 + t18 * t188) * t171 - t69 * t290 + (-t171 * t202 - t278) * qJD(6) + t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100 ^ 2 - t206 ^ 2, -t100 * t25 - t206 * t26 + t88, 0, 0, 0, 0, 0, t11 + t277, t10 + t280; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t301, t300, t302, t297, t178 (-qJD(6) + t171) * t2 + t294, t1 * t171 - t188 * t221 + t299;];
tauc_reg  = t3;
