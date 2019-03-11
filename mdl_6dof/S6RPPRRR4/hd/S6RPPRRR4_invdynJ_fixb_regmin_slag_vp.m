% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPRRR4
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% tau_reg [6x30]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:26:50
% EndTime: 2019-03-09 02:26:59
% DurationCPUTime: 3.70s
% Computational Cost: add. (3187->421), mult. (6005->570), div. (0->0), fcn. (4110->12), ass. (0->222)
t165 = cos(qJ(4));
t250 = qJD(1) * t165;
t124 = qJD(5) + t250;
t162 = sin(qJ(4));
t244 = qJD(4) * t165;
t158 = cos(pkin(10));
t237 = qJD(1) * qJD(2);
t128 = t158 * t237;
t166 = -pkin(1) - pkin(2);
t120 = t166 * qJDD(1) + qJDD(2);
t157 = sin(pkin(10));
t232 = qJDD(1) * t158;
t257 = qJ(2) * t232 + t157 * t120;
t65 = t128 + t257;
t60 = -qJDD(1) * pkin(7) + t65;
t301 = -qJD(3) * qJD(4) - t60;
t123 = t166 * qJD(1) + qJD(2);
t252 = qJ(2) * qJD(1);
t83 = t157 * t123 + t158 * t252;
t73 = -qJD(1) * pkin(7) + t83;
t226 = t162 * t301 - t73 * t244;
t17 = -qJDD(4) * pkin(4) - t165 * qJDD(3) - t226;
t283 = sin(qJ(1));
t284 = cos(qJ(1));
t98 = t284 * t157 - t283 * t158;
t285 = g(2) * t98;
t97 = -t283 * t157 - t284 * t158;
t288 = g(1) * t97;
t199 = t285 + t288;
t174 = -g(3) * t165 + t199 * t162;
t303 = qJD(5) * pkin(8) * t124 + t17 + t174;
t164 = cos(qJ(5));
t161 = sin(qJ(5));
t239 = t161 * qJD(4);
t251 = qJD(1) * t162;
t100 = t164 * t251 - t239;
t160 = sin(qJ(6));
t163 = cos(qJ(6));
t245 = qJD(4) * t164;
t99 = t161 * t251 + t245;
t188 = t163 * t100 - t160 * t99;
t52 = t100 * t160 + t163 * t99;
t302 = t188 * t52;
t55 = qJD(3) * t165 - t162 * t73;
t300 = t55 * qJD(4);
t299 = t188 ^ 2 - t52 ^ 2;
t122 = qJD(6) + t124;
t240 = qJD(6) * t163;
t241 = qJD(6) * t160;
t236 = qJD(1) * qJD(4);
t209 = t165 * t236;
t229 = t162 * qJDD(1);
t180 = t209 + t229;
t243 = qJD(5) * t161;
t212 = t162 * t243;
t234 = qJD(4) * qJD(5);
t40 = qJD(1) * t212 + t161 * qJDD(4) + (-t180 + t234) * t164;
t211 = t165 * t239;
t242 = qJD(5) * t164;
t177 = t162 * t242 + t211;
t41 = t177 * qJD(1) + t164 * qJDD(4) + (t229 - t234) * t161;
t7 = t100 * t241 + t160 * t41 + t163 * t40 + t99 * t240;
t298 = -t122 * t52 + t7;
t56 = t162 * qJD(3) + t165 * t73;
t47 = qJD(4) * pkin(8) + t56;
t195 = pkin(4) * t165 + pkin(8) * t162;
t82 = t158 * t123 - t157 * t252;
t72 = qJD(1) * pkin(3) - t82;
t48 = t195 * qJD(1) + t72;
t19 = t161 * t48 + t164 * t47;
t14 = pkin(9) * t99 + t19;
t12 = t14 * t241;
t155 = qJ(5) + qJ(6);
t147 = cos(t155);
t282 = g(3) * t162;
t46 = -qJD(4) * pkin(4) - t55;
t29 = -pkin(5) * t99 + t46;
t146 = sin(t155);
t264 = t147 * t165;
t32 = t146 * t97 + t98 * t264;
t34 = t146 * t98 - t97 * t264;
t297 = g(1) * t34 - g(2) * t32 - t147 * t282 - t29 * t52 + t12;
t16 = qJDD(4) * pkin(8) + t162 * qJDD(3) + t165 * t60 + t300;
t143 = t165 * qJDD(1);
t210 = t162 * t236;
t293 = -t210 + t143;
t126 = t157 * t237;
t233 = qJDD(1) * t157;
t204 = -qJ(2) * t233 + t158 * t120;
t64 = -t126 + t204;
t59 = qJDD(1) * pkin(3) - t64;
t26 = pkin(4) * t293 + t180 * pkin(8) + t59;
t25 = t164 * t26;
t96 = -qJDD(5) - t293;
t2 = -t96 * pkin(5) - t40 * pkin(9) - t19 * qJD(5) - t161 * t16 + t25;
t228 = -t164 * t16 - t161 * t26 - t48 * t242;
t182 = t47 * t243 + t228;
t3 = pkin(9) * t41 - t182;
t222 = -t160 * t3 + t163 * t2;
t265 = t146 * t165;
t31 = -t147 * t97 + t98 * t265;
t33 = t147 * t98 + t97 * t265;
t18 = -t161 * t47 + t164 * t48;
t13 = pkin(9) * t100 + t18;
t11 = pkin(5) * t124 + t13;
t273 = t14 * t163;
t5 = t11 * t160 + t273;
t296 = -g(1) * t33 - g(2) * t31 - t5 * qJD(6) - t146 * t282 + t29 * t188 + t222;
t172 = t188 * qJD(6) - t160 * t40 + t163 * t41;
t295 = -t122 * t188 + t172;
t292 = t164 * t244 - t212;
t291 = qJD(5) + qJD(6);
t110 = t158 * qJ(2) + t157 * t166;
t102 = -pkin(7) + t110;
t287 = g(1) * t98;
t290 = qJD(5) * (t102 * t124 + t47) + t287;
t289 = pkin(8) + pkin(9);
t286 = g(2) * t97;
t194 = -pkin(4) * t162 + pkin(8) * t165;
t106 = t194 * qJD(1);
t280 = t161 * t106 + t164 * t55;
t260 = t163 * t164;
t103 = t160 * t161 - t260;
t181 = t103 * t165;
t279 = -qJD(1) * t181 - t103 * t291;
t104 = t160 * t164 + t161 * t163;
t58 = t291 * t104;
t278 = t104 * t250 + t58;
t214 = t162 * t245;
t259 = t164 * t165;
t262 = t161 * t165;
t91 = -t157 * t262 - t158 * t164;
t277 = -t91 * qJD(5) + t157 * t214 + (t157 * t161 + t158 * t259) * qJD(1);
t215 = t162 * t239;
t92 = t157 * t259 - t158 * t161;
t276 = -t92 * qJD(5) + t157 * t215 - (t157 * t164 - t158 * t262) * qJD(1);
t88 = t157 * qJD(2) + t194 * qJD(4);
t275 = t102 * t215 + t164 * t88;
t109 = -t157 * qJ(2) + t158 * t166;
t101 = pkin(3) - t109;
t74 = t101 + t195;
t84 = t102 * t259;
t274 = t161 * t74 + t84;
t272 = t161 * t97;
t271 = t40 * t161;
t270 = t99 * t124;
t269 = pkin(1) * qJDD(1);
t268 = qJD(4) * t99;
t267 = t100 * t124;
t266 = t124 * t164;
t263 = t161 * t162;
t261 = t162 * t164;
t258 = qJDD(3) + g(3);
t256 = t284 * pkin(1) + t283 * qJ(2);
t255 = g(1) * t283 - g(2) * t284;
t153 = t162 ^ 2;
t254 = -t165 ^ 2 + t153;
t167 = qJD(4) ^ 2;
t168 = qJD(1) ^ 2;
t253 = t167 + t168;
t249 = qJD(2) * t158;
t247 = qJD(4) * t100;
t246 = qJD(4) * t162;
t238 = qJ(2) * qJDD(1);
t231 = qJDD(4) * t162;
t230 = qJDD(4) * t165;
t218 = t165 * t249;
t227 = t161 * t88 + t164 * t218 + t74 * t242;
t224 = 0.2e1 * t237;
t223 = t160 * t263;
t221 = qJD(5) * t289;
t220 = t158 * t251;
t219 = t161 * t250;
t217 = t124 * t239;
t216 = t124 * t245;
t208 = qJD(6) * t11 + t3;
t206 = -qJD(5) * t48 - t16;
t203 = 0.2e1 * t209;
t202 = qJDD(2) - t269;
t201 = -t56 + (t219 + t243) * pkin(5);
t200 = t286 - t287;
t198 = -qJD(6) * t91 + t277;
t197 = qJD(6) * t92 - t276;
t196 = -t283 * pkin(1) + t284 * qJ(2);
t117 = t289 * t161;
t193 = pkin(9) * t219 + qJD(6) * t117 + t161 * t221 + t280;
t118 = t289 * t164;
t187 = -pkin(5) * t162 + pkin(9) * t259;
t94 = t164 * t106;
t192 = t187 * qJD(1) + qJD(6) * t118 - t161 * t55 + t164 * t221 + t94;
t21 = qJD(4) * t181 + t58 * t162;
t86 = t162 * t260 - t223;
t95 = -qJDD(6) + t96;
t191 = t122 * t21 + t86 * t95;
t22 = -qJD(6) * t223 + (t261 * t291 + t211) * t163 + t292 * t160;
t85 = t104 * t162;
t190 = -t122 * t22 + t85 * t95;
t189 = t157 * t82 - t158 * t83;
t186 = t165 * t172 - t246 * t52;
t185 = t165 * t7 + t188 * t246;
t184 = t124 * t242 - t161 * t96;
t183 = t124 * t243 + t164 * t96;
t179 = g(1) * t284 + g(2) * t283;
t178 = qJD(1) * t72 - t199;
t176 = pkin(8) * t96 + t124 * t46;
t171 = -qJDD(4) * t102 + (-qJD(1) * t101 - t249 - t72) * qJD(4);
t170 = qJDD(1) * t101 - t102 * t167 + t126 + t200 + t59;
t139 = -pkin(5) * t164 - pkin(4);
t113 = -t162 * t167 + t230;
t112 = -t165 * t167 - t231;
t71 = (-pkin(5) * t161 + t102) * t162;
t70 = t164 * t74;
t43 = t161 * t98 - t97 * t259;
t42 = t164 * t98 + t97 * t262;
t30 = -t177 * pkin(5) + t102 * t244 + t162 * t249;
t28 = pkin(9) * t263 + t274;
t27 = pkin(9) * t261 + t70 + (-t102 * t161 + pkin(5)) * t165;
t10 = -pkin(5) * t41 + t17;
t9 = (-t165 * t243 - t214) * t102 + t177 * pkin(9) + t227;
t6 = -t161 * t218 + t187 * qJD(4) + (-t84 + (-pkin(9) * t162 - t74) * t161) * qJD(5) + t275;
t4 = t11 * t163 - t14 * t160;
t1 = [qJDD(1), t255, t179, -qJDD(2) + t255 + 0.2e1 * t269, -t179 + t224 + 0.2e1 * t238, -t202 * pkin(1) - g(1) * t196 - g(2) * t256 + (t224 + t238) * qJ(2), -qJDD(1) * t109 + 0.2e1 * t126 + t200 - t204, qJDD(1) * t110 + 0.2e1 * t128 + t199 + t257, t65 * t110 + t64 * t109 - g(1) * (-t283 * pkin(2) + t196) - g(2) * (t284 * pkin(2) + t256) - t189 * qJD(2), qJDD(1) * t153 + t162 * t203, 0.2e1 * t162 * t143 - 0.2e1 * t254 * t236, t112, -t113, 0, t171 * t162 + t170 * t165, -t170 * t162 + t171 * t165, t100 * t292 - t40 * t261 (-t100 * t161 - t164 * t99) * t244 + (t271 - t164 * t41 + (-t100 * t164 + t161 * t99) * qJD(5)) * t162 (t40 - t216) * t165 + (t183 + t247) * t162 (t41 + t217) * t165 + (t184 - t268) * t162, -t124 * t246 - t165 * t96 (-t243 * t74 + t275) * t124 - t70 * t96 - g(1) * t272 - g(2) * t43 + (-t102 * t268 + t25 - t290 * t164 + (-qJD(4) * t46 + t102 * t96 - t124 * t249 + t206) * t161) * t165 + (-qJD(4) * t18 - t102 * t41 - t17 * t161 - t242 * t46 - t249 * t99) * t162, -t227 * t124 + t274 * t96 - t164 * t288 - g(2) * t42 + ((-t100 * t102 - t164 * t46) * qJD(4) + t290 * t161 + t228) * t165 + (-t100 * t249 + t46 * t243 + t102 * t40 - t17 * t164 + (t102 * t266 + t19) * qJD(4)) * t162, -t188 * t21 - t7 * t86, -t172 * t86 - t188 * t22 + t21 * t52 + t7 * t85, t185 + t191, t186 - t190, -t122 * t246 - t165 * t95 (-t160 * t9 + t163 * t6) * t122 - (-t160 * t28 + t163 * t27) * t95 + t222 * t165 - t4 * t246 - t30 * t52 - t71 * t172 - t10 * t85 - t29 * t22 - g(1) * t32 - g(2) * t34 + ((-t160 * t27 - t163 * t28) * t122 - t5 * t165) * qJD(6), t12 * t165 + t5 * t246 - t30 * t188 + t71 * t7 - t10 * t86 + t29 * t21 + g(1) * t31 - g(2) * t33 + (-(-qJD(6) * t28 + t6) * t122 + t27 * t95 - t2 * t165) * t160 + (-(qJD(6) * t27 + t9) * t122 + t28 * t95 - t208 * t165) * t163; 0, 0, 0, -qJDD(1), -t168, -qJ(2) * t168 + t202 - t255, -t157 * t168 - t232, -t158 * t168 + t233, t189 * qJD(1) + t65 * t157 + t64 * t158 - t255, 0, 0, 0, 0, 0 (0.2e1 * t210 - t143) * t158 + (-t253 * t165 - t231) * t157 (t203 + t229) * t158 + (t253 * t162 - t230) * t157, 0, 0, 0, 0, 0, t99 * t220 - t91 * t96 + (-t162 * t41 - t244 * t99) * t157 + t276 * t124, t162 * t157 * t40 + t92 * t96 + t277 * t124 + (-t157 * t244 + t220) * t100, 0, 0, 0, 0, 0 -(-t160 * t92 + t163 * t91) * t95 + t52 * t220 + (-t162 * t172 - t244 * t52) * t157 + (t160 * t198 - t163 * t197) * t122 (t160 * t91 + t163 * t92) * t95 + t188 * t220 + (t162 * t7 - t188 * t244) * t157 + (t160 * t197 + t163 * t198) * t122; 0, 0, 0, 0, 0, 0, 0, 0, t258, 0, 0, 0, 0, 0, t113, t112, 0, 0, 0, 0, 0 (t41 - t217) * t165 + (-t184 - t268) * t162 (-t40 - t216) * t165 + (t183 - t247) * t162, 0, 0, 0, 0, 0, t186 + t190, -t185 + t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162 * t168 * t165, t254 * t168, -t229, -t143, qJDD(4), t56 * qJD(4) + t178 * t162 + t258 * t165 + t226, t300 + (qJD(4) * t73 - t258) * t162 + (t178 + t301) * t165, -t100 * t266 + t271 (t40 + t270) * t164 + (t41 + t267) * t161 (-t100 * t162 + t124 * t259) * qJD(1) + t184 (-t124 * t262 + t162 * t99) * qJD(1) - t183, t124 * t251, t18 * t251 + pkin(4) * t41 - t94 * t124 + t56 * t99 + (t55 * t124 + t176) * t161 - t303 * t164, -pkin(4) * t40 + t56 * t100 + t280 * t124 + t303 * t161 + t176 * t164 - t19 * t251, t7 * t104 - t188 * t279, -t7 * t103 + t104 * t172 + t188 * t278 + t279 * t52, -t104 * t95 + t122 * t279 - t188 * t251, t103 * t95 - t122 * t278 + t251 * t52, t122 * t251 -(-t117 * t163 - t118 * t160) * t95 - t139 * t172 + t10 * t103 + t4 * t251 - t201 * t52 + t278 * t29 + (t160 * t193 - t163 * t192) * t122 - t174 * t147 (-t117 * t160 + t118 * t163) * t95 + t139 * t7 + t10 * t104 - t5 * t251 - t201 * t188 + t279 * t29 + (t160 * t192 + t163 * t193) * t122 + t174 * t146; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100 * t99, t100 ^ 2 - t99 ^ 2, t40 - t270, t41 - t267, -t96, -g(1) * t42 + t46 * t100 + t19 * t124 + t25 + (-qJD(5) * t47 + t286) * t164 + (-t165 * t285 + t206 - t282) * t161, t18 * t124 - t46 * t99 + g(1) * t43 - g(2) * (t259 * t98 + t272) - g(3) * t261 + t182, t302, t299, t298, t295, -t95 -(-t13 * t160 - t273) * t122 + (-t100 * t52 - t122 * t241 - t163 * t95) * pkin(5) + t296 (-t122 * t14 - t2) * t160 + (t122 * t13 - t208) * t163 + (-t100 * t188 - t122 * t240 + t160 * t95) * pkin(5) + t297; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t302, t299, t298, t295, -t95, t5 * t122 + t296, t4 * t122 - t160 * t2 - t163 * t208 + t297;];
tau_reg  = t1;
