% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRP7
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
% tau_reg [5x28]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRP7_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP7_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP7_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:57:47
% EndTime: 2019-12-31 21:57:57
% DurationCPUTime: 4.19s
% Computational Cost: add. (5172->409), mult. (11537->520), div. (0->0), fcn. (8081->10), ass. (0->216)
t158 = sin(qJ(3));
t162 = cos(qJ(2));
t289 = cos(qJ(3));
t222 = t289 * t162;
t159 = sin(qJ(2));
t234 = t159 * qJDD(1);
t247 = t158 * t162;
t110 = t289 * t159 + t247;
t232 = qJD(2) + qJD(3);
t77 = t232 * t110;
t52 = t77 * qJD(1) - qJDD(1) * t222 + t158 * t234;
t51 = qJDD(4) + t52;
t160 = sin(qJ(1));
t163 = cos(qJ(1));
t206 = g(1) * t163 + g(2) * t160;
t161 = cos(qJ(4));
t236 = qJD(4) * t161;
t220 = qJD(1) * t289;
t240 = qJD(1) * t159;
t102 = t158 * t240 - t162 * t220;
t255 = t102 * t161;
t307 = t236 + t255;
t157 = sin(qJ(4));
t237 = qJD(4) * t157;
t306 = t102 * t157 + t237;
t156 = qJ(2) + qJ(3);
t152 = sin(t156);
t252 = t152 * t163;
t253 = t152 * t160;
t303 = g(1) * t252 + g(2) * t253;
t204 = t161 * pkin(4) + t157 * qJ(5);
t153 = cos(t156);
t302 = t206 * t153;
t292 = pkin(7) + pkin(6);
t235 = qJD(1) * qJD(2);
t218 = t159 * t235;
t140 = pkin(2) * t218;
t190 = -t158 * t159 + t222;
t76 = t232 * t190;
t171 = t76 * qJD(1);
t150 = t162 * pkin(2) + pkin(1);
t295 = -pkin(8) * t110 - t150;
t17 = t52 * pkin(3) - pkin(8) * t171 + t295 * qJDD(1) + t140;
t121 = t292 * t159;
t111 = qJD(1) * t121;
t270 = qJD(2) * pkin(2);
t107 = -t111 + t270;
t122 = t292 * t162;
t113 = qJD(1) * t122;
t219 = qJD(3) * t289;
t239 = qJD(3) * t158;
t217 = t162 * t235;
t79 = qJDD(2) * pkin(2) + t292 * (-t217 - t234);
t233 = t162 * qJDD(1);
t81 = t292 * (-t218 + t233);
t177 = t107 * t219 - t113 * t239 + t158 * t79 + t289 * t81;
t231 = qJDD(2) + qJDD(3);
t25 = t231 * pkin(8) + t177;
t104 = -qJD(1) * t247 - t159 * t220;
t120 = t150 * qJD(1);
t60 = t102 * pkin(3) + t104 * pkin(8) - t120;
t106 = t289 * t113;
t72 = t158 * t107 + t106;
t63 = t232 * pkin(8) + t72;
t193 = t157 * t17 + t161 * t25 + t60 * t236 - t63 * t237;
t259 = t51 * qJ(5);
t98 = qJD(4) + t102;
t2 = qJD(5) * t98 + t193 + t259;
t216 = t157 * t25 - t161 * t17 + t63 * t236 + t60 * t237;
t291 = pkin(4) * t51;
t4 = qJDD(5) + t216 - t291;
t301 = t4 * t157 + t2 * t161;
t263 = t161 * t51;
t195 = t98 * t237 - t263;
t298 = t195 * pkin(8);
t74 = -t158 * t111 + t106;
t211 = pkin(2) * t239 - t74;
t297 = t306 * pkin(4) - t307 * qJ(5) - qJD(5) * t157;
t296 = t153 * pkin(3) + t152 * pkin(8);
t169 = t110 * qJDD(1) + t171;
t186 = t161 * t104 - t157 * t232;
t34 = -t186 * qJD(4) + t169 * t157 - t161 * t231;
t294 = t186 ^ 2;
t293 = t98 ^ 2;
t290 = pkin(8) * t51;
t288 = pkin(4) * t104;
t144 = g(3) * t152;
t284 = g(3) * t153;
t283 = g(3) * t157;
t282 = g(3) * t162;
t36 = t157 * t60 + t161 * t63;
t280 = t36 * t98;
t213 = t161 * t232;
t82 = -t104 * t157 - t213;
t279 = t82 * t98;
t278 = t186 * t82;
t277 = t186 * t98;
t276 = -t297 - t211;
t68 = -pkin(3) * t104 + pkin(8) * t102;
t105 = t158 * t113;
t71 = t289 * t107 - t105;
t275 = t157 * t68 + t161 * t71;
t61 = pkin(2) * t240 + t68;
t75 = -t289 * t111 - t105;
t274 = t157 * t61 + t161 * t75;
t70 = -pkin(3) * t190 + t295;
t87 = -t158 * t121 + t289 * t122;
t273 = t157 * t70 + t161 * t87;
t272 = -t72 + t297;
t271 = pkin(8) * qJD(4);
t62 = -t232 * pkin(3) - t71;
t30 = t82 * pkin(4) + qJ(5) * t186 + t62;
t269 = t102 * t30;
t148 = pkin(2) * t158 + pkin(8);
t268 = t148 * t51;
t267 = t157 * t51;
t266 = t157 * t76;
t265 = t157 * t82;
t264 = t161 * t34;
t262 = t161 * t76;
t261 = t161 * t186;
t214 = t161 * t98;
t33 = -qJD(4) * t213 - t104 * t237 - t157 * t231 - t161 * t169;
t260 = t33 * t157;
t258 = t62 * t102;
t257 = t98 * t104;
t254 = t104 * t102;
t248 = t157 * t160;
t246 = t160 * t161;
t245 = t161 * t163;
t244 = t163 * t157;
t35 = -t157 * t63 + t161 * t60;
t243 = qJD(5) - t35;
t242 = t303 * t161;
t154 = t159 ^ 2;
t241 = -t162 ^ 2 + t154;
t238 = qJD(4) * t148;
t230 = t289 * pkin(2);
t229 = t159 * t270;
t27 = t30 * t237;
t56 = t62 * t237;
t227 = t153 * t283 - t303 * t157;
t226 = t144 + t302;
t225 = qJD(2) * t292;
t224 = t157 * t289;
t223 = t161 * t289;
t215 = t107 * t239 + t113 * t219 + t158 * t81 - t289 * t79;
t26 = -t231 * pkin(3) + t215;
t221 = -t26 - t284;
t212 = pkin(2) * t219;
t93 = t153 * t248 + t245;
t95 = t153 * t244 - t246;
t210 = -g(1) * t93 + g(2) * t95;
t94 = t153 * t246 - t244;
t96 = t153 * t245 + t248;
t209 = g(1) * t94 - g(2) * t96;
t6 = t34 * pkin(4) + t33 * qJ(5) + qJD(5) * t186 + t26;
t208 = -t98 * t238 - t6;
t205 = g(1) * t160 - g(2) * t163;
t203 = pkin(4) * t157 - qJ(5) * t161;
t20 = -pkin(4) * t98 + t243;
t21 = qJ(5) * t98 + t36;
t201 = -t157 * t21 + t161 * t20;
t200 = -t20 * t104 + t242 + t27;
t199 = t35 * t104 + t242 + t56;
t198 = pkin(3) + t204;
t197 = t150 + t296;
t196 = t6 * t157 + t30 * t236;
t194 = -0.2e1 * pkin(1) * t235 - pkin(6) * qJDD(2);
t44 = pkin(3) * t77 - pkin(8) * t76 + t229;
t112 = t159 * t225;
t114 = t162 * t225;
t191 = -t289 * t121 - t158 * t122;
t47 = t191 * qJD(3) - t289 * t112 - t158 * t114;
t192 = t157 * t44 + t161 * t47 + t70 * t236 - t87 * t237;
t189 = t21 * t104 - t30 * t255 - t227;
t185 = -t36 * t104 + t26 * t157 + t62 * t236 + t227;
t184 = -t98 * t212 - t268;
t182 = g(1) * t95 + g(2) * t93 + t152 * t283 - t216;
t165 = qJD(2) ^ 2;
t181 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t165 + t205;
t166 = qJD(1) ^ 2;
t180 = pkin(1) * t166 - pkin(6) * qJDD(1) + t206;
t179 = -t120 * t104 - t215 - t284 + t303;
t178 = t201 * qJD(4) + t301;
t176 = -t260 - t264 + (-t261 + t265) * qJD(4);
t175 = -t186 * t30 + qJDD(5) - t182;
t174 = t307 * t20 - t306 * t21 - t226 + t301;
t173 = -g(1) * t96 - g(2) * t94 - t161 * t144 + t193;
t48 = t87 * qJD(3) - t158 * t112 + t289 * t114;
t172 = -t120 * t102 - t177 + t226;
t168 = -g(3) * (t204 * t153 + t296) + t206 * t152 * t198 - t302 * pkin(8);
t149 = -t230 - pkin(3);
t108 = -t230 - t198;
t99 = -t150 * qJDD(1) + t140;
t92 = t104 * qJ(5);
t53 = -t102 ^ 2 + t104 ^ 2;
t49 = -pkin(4) * t186 + qJ(5) * t82;
t46 = t203 * t110 - t191;
t43 = -t104 * t232 - t52;
t42 = t102 * t232 + t169;
t38 = pkin(4) * t190 + t157 * t87 - t161 * t70;
t37 = -qJ(5) * t190 + t273;
t32 = t157 * t71 - t161 * t68 + t288;
t31 = -t92 + t275;
t29 = t157 * t75 - t161 * t61 + t288;
t28 = -t92 + t274;
t13 = -t33 + t279;
t12 = -t104 * t186 + t98 * t214 + t267;
t11 = -t82 * t104 - t293 * t157 + t263;
t10 = -t186 * t214 - t260;
t9 = t203 * t76 + (t204 * qJD(4) - qJD(5) * t161) * t110 + t48;
t8 = -t77 * pkin(4) + t273 * qJD(4) + t157 * t47 - t161 * t44;
t7 = qJ(5) * t77 - qJD(5) * t190 + t192;
t5 = (-t33 - t279) * t161 + (-t34 + t277) * t157;
t1 = [qJDD(1), t205, t206, qJDD(1) * t154 + 0.2e1 * t159 * t217, 0.2e1 * t159 * t233 - 0.2e1 * t241 * t235, qJDD(2) * t159 + t162 * t165, qJDD(2) * t162 - t159 * t165, 0, t194 * t159 + t181 * t162, -t181 * t159 + t162 * t194, -t104 * t76 + t169 * t110, -t76 * t102 + t104 * t77 - t110 * t52 + t169 * t190, t110 * t231 + t232 * t76, t190 * t231 - t232 * t77, 0, t102 * t229 - t120 * t77 - t150 * t52 + t153 * t205 - t190 * t99 + t191 * t231 - t232 * t48, -g(1) * t253 + g(2) * t252 - t104 * t229 + t99 * t110 - t120 * t76 - t150 * t169 - t231 * t87 - t232 * t47, -t76 * t261 + (-t33 * t161 + t186 * t237) * t110, (t157 * t186 - t161 * t82) * t76 + (t260 - t264 + (t261 + t265) * qJD(4)) * t110, -t110 * t195 - t186 * t77 + t190 * t33 + t214 * t76, -t98 * t266 + t34 * t190 - t82 * t77 + (-t236 * t98 - t267) * t110, -t190 * t51 + t77 * t98, t216 * t190 + t35 * t77 + t48 * t82 - t191 * t34 + ((-qJD(4) * t87 + t44) * t98 + t70 * t51 + t62 * qJD(4) * t110) * t161 + ((-qJD(4) * t70 - t47) * t98 - t87 * t51 + t26 * t110 + t62 * t76) * t157 + t209, -t192 * t98 - t273 * t51 + t193 * t190 - t36 * t77 - t48 * t186 + t191 * t33 + t62 * t262 + (t26 * t161 - t56) * t110 + t210, t110 * t196 + t190 * t4 - t20 * t77 + t30 * t266 + t46 * t34 - t38 * t51 - t8 * t98 + t9 * t82 + t209, -t38 * t33 - t37 * t34 - t7 * t82 - t8 * t186 + t201 * t76 + t205 * t152 + (-t157 * t2 + t161 * t4 + (-t157 * t20 - t161 * t21) * qJD(4)) * t110, -t30 * t262 - t2 * t190 + t21 * t77 + t46 * t33 + t37 * t51 + t7 * t98 + t9 * t186 + (-t6 * t161 + t27) * t110 - t210, t2 * t37 + t21 * t7 + t6 * t46 + t30 * t9 + t4 * t38 + t20 * t8 - g(1) * (-t94 * pkin(4) - t93 * qJ(5)) - g(2) * (t96 * pkin(4) + t95 * qJ(5)) + (-g(1) * t292 - g(2) * t197) * t163 + (g(1) * t197 - g(2) * t292) * t160; 0, 0, 0, -t159 * t166 * t162, t241 * t166, t234, t233, qJDD(2), t180 * t159 - t282, g(3) * t159 + t162 * t180, -t254, t53, t42, t43, t231, t74 * t232 + (-t102 * t240 + t289 * t231 - t232 * t239) * pkin(2) + t179, t75 * t232 + (t104 * t240 - t158 * t231 - t219 * t232) * pkin(2) + t172, t10, t5, t12, t11, t257, t149 * t34 + t211 * t82 + ((-t61 - t238) * t98 + t221) * t161 + (t258 - t268 + (-t212 + t75) * t98) * t157 + t199, -t149 * t33 + (t148 * t237 + t274) * t98 - t211 * t186 + (t184 + t258) * t161 + t185, t108 * t34 + t29 * t98 - t276 * t82 + (t208 - t284) * t161 + (t184 + t269) * t157 + t200, t28 * t82 + t29 * t186 + (-t186 * t224 - t223 * t82) * qJD(3) * pkin(2) + t176 * t148 + t174, t108 * t33 - t28 * t98 - t276 * t186 + t208 * t157 + (-qJD(4) * t30 - t184) * t161 + t189, t6 * t108 - t21 * t28 - t20 * t29 - t276 * t30 + (-t282 + t206 * t159 + (t20 * t224 + t21 * t223) * qJD(3)) * pkin(2) + t178 * t148 + t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t254, t53, t42, t43, t231, t232 * t72 + t179, t232 * t71 + t172, t10, t5, t12, t11, t257, -pkin(3) * t34 - t72 * t82 + (t71 * t98 + t258 - t290) * t157 + ((-t68 - t271) * t98 + t221) * t161 + t199, pkin(3) * t33 + t186 * t72 + t62 * t255 + t275 * t98 + t185 + t298, -t198 * t34 + t32 * t98 + t272 * t82 + (t269 - t290) * t157 + (-t98 * t271 - t284 - t6) * t161 + t200, pkin(8) * t176 + t186 * t32 + t31 * t82 + t174, t186 * t272 - t198 * t33 - t31 * t98 + t189 - t196 - t298, t178 * pkin(8) - t198 * t6 - t20 * t32 - t21 * t31 + t272 * t30 + t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t278, -t82 ^ 2 + t294, t13, -t34 - t277, t51, t186 * t62 + t182 + t280, t35 * t98 + t62 * t82 - t173, -t49 * t82 - t175 + t280 + 0.2e1 * t291, pkin(4) * t33 - t34 * qJ(5) - (t21 - t36) * t186 + (t20 - t243) * t82, 0.2e1 * t259 - t30 * t82 - t49 * t186 + (0.2e1 * qJD(5) - t35) * t98 + t173, t2 * qJ(5) - t4 * pkin(4) - t30 * t49 - t20 * t36 - g(1) * (-pkin(4) * t95 + qJ(5) * t96) - g(2) * (-pkin(4) * t93 + qJ(5) * t94) + t243 * t21 + t203 * t144; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t278 - t51, t13, -t293 - t294, -t21 * t98 + t175 - t291;];
tau_reg = t1;
