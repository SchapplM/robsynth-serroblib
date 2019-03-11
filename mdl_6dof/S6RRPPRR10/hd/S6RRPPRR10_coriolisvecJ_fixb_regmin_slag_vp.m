% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta4]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRR10_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR10_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR10_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR10_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:37:01
% EndTime: 2019-03-09 09:37:12
% DurationCPUTime: 3.89s
% Computational Cost: add. (4351->394), mult. (10289->554), div. (0->0), fcn. (7289->8), ass. (0->213)
t202 = sin(qJ(2));
t272 = qJD(1) * t202;
t182 = qJD(5) + t272;
t175 = qJD(6) + t182;
t200 = sin(qJ(6));
t203 = cos(qJ(6));
t198 = cos(pkin(10));
t205 = cos(qJ(2));
t271 = qJD(1) * t205;
t254 = t198 * t271;
t197 = sin(pkin(10));
t270 = qJD(2) * t197;
t149 = t254 + t270;
t252 = t197 * t271;
t269 = qJD(2) * t198;
t151 = -t252 + t269;
t201 = sin(qJ(5));
t204 = cos(qJ(5));
t88 = t149 * t201 - t151 * t204;
t89 = t204 * t149 + t151 * t201;
t305 = t200 * t88 - t203 * t89;
t316 = t305 * t175;
t153 = t197 * t204 + t198 * t201;
t138 = t153 * qJD(5);
t216 = t153 * t202;
t276 = -qJD(1) * t216 - t138;
t298 = pkin(3) + pkin(7);
t315 = t182 * t89;
t228 = t200 * t89 + t203 * t88;
t314 = t228 * t305;
t313 = t88 * t182;
t312 = t175 * t228;
t255 = t198 * t272;
t265 = qJD(5) * t204;
t266 = qJD(5) * t201;
t284 = t197 * t201;
t275 = -t197 * t266 + t198 * t265 + t204 * t255 - t272 * t284;
t223 = -t198 * t204 + t284;
t185 = pkin(7) * t272;
t311 = qJD(3) + t185;
t310 = t228 ^ 2 - t305 ^ 2;
t199 = -pkin(2) - qJ(4);
t250 = -t202 * qJ(3) - pkin(1);
t148 = t199 * t205 + t250;
t114 = t148 * qJD(1);
t260 = pkin(3) * t272 + t311;
t120 = t199 * qJD(2) + t260;
t61 = -t114 * t197 + t198 * t120;
t37 = pkin(4) * t272 - pkin(8) * t151 + t61;
t62 = t198 * t114 + t197 * t120;
t44 = -pkin(8) * t149 + t62;
t15 = t201 * t37 + t204 * t44;
t13 = -pkin(9) * t89 + t15;
t264 = qJD(6) * t200;
t11 = t13 * t264;
t186 = pkin(7) * t271;
t187 = pkin(3) * t271;
t160 = t186 + t187;
t194 = qJD(2) * qJ(3);
t301 = qJD(4) + t194;
t133 = t160 + t301;
t97 = pkin(4) * t149 + t133;
t38 = pkin(5) * t89 + t97;
t309 = -t305 * t38 + t11;
t263 = qJD(6) * t203;
t259 = qJD(1) * qJD(2);
t251 = t202 * t259;
t237 = t204 * t251;
t238 = t201 * t251;
t51 = -t149 * t265 - t151 * t266 + t197 * t237 + t198 * t238;
t52 = -qJD(5) * t88 + t197 * t238 - t198 * t237;
t8 = -t200 * t52 + t203 * t51 - t89 * t263 + t264 * t88;
t308 = t8 - t316;
t184 = t205 * t259;
t283 = t197 * t202;
t221 = pkin(4) * t205 - pkin(8) * t283;
t214 = t221 * qJD(2);
t180 = pkin(7) * t184;
t128 = t180 + (-qJD(4) + t187) * qJD(2);
t181 = pkin(2) * t251;
t225 = -qJ(3) * t205 + qJ(4) * t202;
t262 = t202 * qJD(3);
t209 = qJD(2) * t225 - t205 * qJD(4) - t262;
t87 = qJD(1) * t209 + t181;
t49 = t198 * t128 - t197 * t87;
t25 = qJD(1) * t214 + t49;
t268 = qJD(2) * t202;
t239 = pkin(8) * t198 * t268;
t50 = t197 * t128 + t198 * t87;
t36 = qJD(1) * t239 + t50;
t247 = -t201 * t36 + t204 * t25;
t212 = -qJD(5) * t15 + t247;
t4 = pkin(5) * t184 - t51 * pkin(9) + t212;
t219 = t201 * t25 + t204 * t36 + t37 * t265 - t44 * t266;
t5 = -pkin(9) * t52 + t219;
t256 = -t200 * t5 + t203 * t4;
t14 = -t201 * t44 + t204 * t37;
t12 = pkin(9) * t88 + t14;
t10 = pkin(5) * t182 + t12;
t291 = t13 * t203;
t3 = t10 * t200 + t291;
t307 = -qJD(6) * t3 + t38 * t228 + t256;
t211 = qJD(6) * t228 - t200 * t51 - t203 * t52;
t306 = t211 - t312;
t304 = -0.2e1 * t259;
t303 = t276 * t203;
t94 = -t153 * t200 - t203 * t223;
t297 = -pkin(8) + t199;
t162 = t297 * t197;
t163 = t297 * t198;
t274 = t204 * t162 + t201 * t163;
t189 = pkin(2) * t272;
t129 = qJD(1) * t225 + t189;
t82 = -t197 * t129 + t198 * t160;
t57 = qJD(1) * t221 + t82;
t83 = t198 * t129 + t197 * t160;
t66 = pkin(8) * t255 + t83;
t302 = t223 * qJD(4) - t274 * qJD(5) + t201 * t66 - t204 * t57;
t257 = -pkin(4) * t198 - pkin(3);
t261 = -t257 * t272 + t311;
t267 = qJD(2) * t205;
t300 = t202 * (-t133 + t301) - t199 * t267;
t299 = qJD(4) * t153 + t162 * t266 - t163 * t265 + t201 * t57 + t204 * t66;
t93 = t203 * t153 - t200 * t223;
t296 = -qJD(6) * t93 - t275 * t200 + t303;
t295 = t94 * qJD(6) + t276 * t200 + t275 * t203;
t172 = t298 * t202;
t156 = t198 * t172;
t75 = t202 * pkin(4) + t156 + (pkin(8) * t205 - t148) * t197;
t280 = t198 * t205;
t96 = t198 * t148 + t197 * t172;
t79 = -pkin(8) * t280 + t96;
t293 = t201 * t75 + t204 * t79;
t292 = qJD(2) * pkin(2);
t290 = t205 * t61;
t289 = t205 * t62;
t288 = t275 * pkin(5) + t261;
t188 = pkin(2) * t268;
t100 = t188 + t209;
t161 = t298 * t267;
t65 = t198 * t100 + t197 * t161;
t159 = t298 * t268;
t193 = qJD(2) * qJD(3);
t127 = -qJD(1) * t159 + t193;
t287 = t127 * t197;
t195 = t202 ^ 2;
t207 = qJD(1) ^ 2;
t285 = t195 * t207;
t282 = t198 * t202;
t279 = t205 * t207;
t206 = qJD(2) ^ 2;
t278 = t206 * t202;
t277 = t206 * t205;
t183 = t197 * pkin(4) + qJ(3);
t173 = t298 * t205;
t273 = -t205 ^ 2 + t195;
t169 = -pkin(2) * t205 + t250;
t144 = qJD(1) * t169;
t258 = t202 * t279;
t137 = pkin(4) * t280 + t173;
t249 = qJD(6) * t10 + t5;
t64 = -t197 * t100 + t198 * t161;
t48 = t214 + t64;
t54 = t239 + t65;
t246 = -t201 * t54 + t204 * t48;
t244 = -t201 * t79 + t204 * t75;
t243 = t275 * t182;
t242 = pkin(1) * t304;
t241 = qJD(3) - t292;
t240 = -t162 * t201 + t204 * t163;
t70 = -pkin(9) * t153 + t274;
t236 = pkin(5) * t271 + t276 * pkin(9) + qJD(6) * t70 - t302;
t69 = pkin(9) * t223 + t240;
t235 = t275 * pkin(9) - qJD(6) * t69 + t299;
t234 = t276 * t182 - t223 * t184;
t232 = qJD(6) * t223 - t275;
t231 = t50 * t197 + t49 * t198;
t230 = -t197 * t61 + t198 * t62;
t125 = t153 * t205;
t20 = pkin(5) * t202 + pkin(9) * t125 + t244;
t124 = t223 * t205;
t21 = pkin(9) * t124 + t293;
t229 = t20 * t200 + t203 * t21;
t224 = t203 * t124 + t125 * t200;
t73 = t124 * t200 - t125 * t203;
t222 = -0.2e1 * qJD(2) * t144;
t217 = -qJ(3) * t267 - t262;
t115 = qJD(1) * t217 + t181;
t131 = t188 + t217;
t220 = pkin(7) * t206 + qJD(1) * t131 + t115;
t218 = t201 * t48 + t204 * t54 + t75 * t265 - t79 * t266;
t117 = (-pkin(7) + t257) * t268;
t99 = qJD(1) * t117 + t193;
t167 = pkin(7) * t251 - t193;
t168 = t185 + t241;
t171 = -t186 - t194;
t208 = -t167 * t205 + (t168 * t205 + (t171 + t186) * t202) * qJD(2);
t174 = t202 * t184;
t157 = -qJ(3) * t271 + t189;
t118 = t144 * t272;
t111 = pkin(5) * t153 + t183;
t95 = -t148 * t197 + t156;
t86 = -pkin(5) * t124 + t137;
t78 = t205 * t138 - t223 * t268;
t77 = qJD(2) * t216 + qJD(5) * t124;
t43 = -t78 * pkin(5) + t117;
t22 = t52 * pkin(5) + t99;
t17 = qJD(6) * t73 + t200 * t77 - t203 * t78;
t16 = qJD(6) * t224 + t200 * t78 + t203 * t77;
t7 = pkin(9) * t78 + t218;
t6 = pkin(5) * t267 - t77 * pkin(9) - qJD(5) * t293 + t246;
t2 = t10 * t203 - t13 * t200;
t1 = [0, 0, 0, 0.2e1 * t174, t273 * t304, t277, -t278, 0, -pkin(7) * t277 + t202 * t242, pkin(7) * t278 + t205 * t242, t208, t202 * t222 + t205 * t220, -t202 * t220 + t205 * t222, pkin(7) * t208 + t115 * t169 + t144 * t131, t127 * t280 - t159 * t149 + (qJD(1) * t64 + t49) * t202 + (-t133 * t282 + t290 + (-t173 * t282 + t205 * t95) * qJD(1)) * qJD(2), -t205 * t287 - t159 * t151 + (-qJD(1) * t65 - t50) * t202 + (t133 * t283 - t289 + (t173 * t283 - t205 * t96) * qJD(1)) * qJD(2), -t65 * t149 - t64 * t151 + (t197 * t49 - t198 * t50) * t205 + ((-t197 * t95 + t198 * t96) * qJD(1) + t230) * t268, t127 * t173 - t133 * t159 + t49 * t95 + t50 * t96 + t61 * t64 + t62 * t65, -t125 * t51 - t77 * t88, t124 * t51 + t125 * t52 - t77 * t89 - t78 * t88, t77 * t182 + t51 * t202 + (-qJD(1) * t125 - t88) * t267, t78 * t182 - t52 * t202 + (qJD(1) * t124 - t89) * t267, t182 * t267 + t174, t246 * t182 + t247 * t202 + t117 * t89 + t137 * t52 - t99 * t124 - t97 * t78 + (-t15 * t202 - t182 * t293) * qJD(5) + (qJD(1) * t244 + t14) * t267, -t218 * t182 - t219 * t202 - t117 * t88 + t137 * t51 - t99 * t125 + t97 * t77 + (-t293 * qJD(1) - t15) * t267, -t16 * t228 + t73 * t8, t16 * t305 + t17 * t228 + t211 * t73 + t224 * t8, t16 * t175 + t8 * t202 + (qJD(1) * t73 - t228) * t267, -t17 * t175 + t211 * t202 + (qJD(1) * t224 + t305) * t267, t175 * t267 + t174 (-t200 * t7 + t203 * t6) * t175 + t256 * t202 - t43 * t305 - t86 * t211 - t22 * t224 + t38 * t17 + (-t175 * t229 - t202 * t3) * qJD(6) + ((t20 * t203 - t200 * t21) * qJD(1) + t2) * t267, t11 * t202 + t38 * t16 + t22 * t73 - t43 * t228 + t86 * t8 + (-(-qJD(6) * t21 + t6) * t175 - t4 * t202) * t200 + (-(qJD(6) * t20 + t7) * t175 - t249 * t202) * t203 + (-qJD(1) * t229 - t3) * t267; 0, 0, 0, -t258, t273 * t207, 0, 0, 0, t207 * pkin(1) * t202, pkin(1) * t279 ((-t171 - t194) * t202 + (-t168 + t241) * t205) * qJD(1), -t157 * t271 + t118, 0.2e1 * t193 + (t144 * t205 + t157 * t202) * qJD(1), -t167 * qJ(3) - t171 * qJD(3) - t144 * t157 + (-t171 * t202 + (-t168 - t292) * t205) * qJD(1) * pkin(7), t287 + t260 * t149 + (-t300 * t198 - t202 * t82 - t290) * qJD(1), t127 * t198 + t260 * t151 + (t300 * t197 + t202 * t83 + t289) * qJD(1), t83 * t149 + t82 * t151 + (qJD(4) * t151 - t272 * t62 - t49) * t198 + (qJD(4) * t149 + t272 * t61 - t50) * t197, t127 * qJ(3) - t61 * t82 - t62 * t83 + t231 * t199 + t260 * t133 + (-t197 * t62 - t198 * t61) * qJD(4), -t223 * t51 - t276 * t88, -t51 * t153 + t223 * t52 + t275 * t88 - t276 * t89, t271 * t88 + t234, -t243 + (-qJD(2) * t153 + t89) * t271, -t182 * t271, t99 * t153 + t183 * t52 + t275 * t97 + t261 * t89 + t302 * t182 + (qJD(2) * t240 - t14) * t271, -t99 * t223 + t183 * t51 + t276 * t97 - t261 * t88 + t299 * t182 + (-qJD(2) * t274 + t15) * t271, -t228 * t296 + t8 * t94, t211 * t94 + t228 * t295 + t296 * t305 - t8 * t93, t296 * t175 + (qJD(2) * t94 + t228) * t271, -t295 * t175 + (-qJD(2) * t93 - t305) * t271, -t175 * t271, -t111 * t211 + t22 * t93 + t295 * t38 - t288 * t305 + (t200 * t235 - t203 * t236) * t175 + ((-t200 * t70 + t203 * t69) * qJD(2) - t2) * t271, t111 * t8 + t22 * t94 + t296 * t38 - t288 * t228 + (t200 * t236 + t203 * t235) * t175 + (-(t200 * t69 + t203 * t70) * qJD(2) + t3) * t271; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t258, -t206 - t285, qJD(2) * t171 + t118 + t180, -t197 * t285 + (-t149 + t254) * qJD(2), -t198 * t285 + (-t151 - t252) * qJD(2) (-t149 * t198 + t151 * t197) * t272, -t133 * qJD(2) + t230 * t272 + t231, 0, 0, 0, 0, 0, -qJD(2) * t89 + t234, -t243 + (-t153 * t271 + t88) * qJD(2), 0, 0, 0, 0, 0 (-t153 * t263 + t200 * t232 + t303) * t175 + (t271 * t94 + t305) * qJD(2) (t232 * t203 + (qJD(6) * t153 - t276) * t200) * t175 + (-t271 * t93 + t228) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t151 - t269) * t272 (-t149 + t270) * t272, -t149 ^ 2 - t151 ^ 2, t62 * t149 + t61 * t151 + t127, 0, 0, 0, 0, 0, t52 - t313, t51 - t315, 0, 0, 0, 0, 0, -t211 - t312, t8 + t316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88 * t89, t88 ^ 2 - t89 ^ 2, t51 + t315, -t52 - t313, t184, t15 * t182 + t88 * t97 + t212, t14 * t182 + t89 * t97 - t219, t314, t310, t308, t306, t184 -(-t12 * t200 - t291) * t175 + (-t175 * t264 + t184 * t203 - t305 * t88) * pkin(5) + t307 (-t13 * t175 - t4) * t200 + (t12 * t175 - t249) * t203 + (-t175 * t263 - t184 * t200 - t228 * t88) * pkin(5) + t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t314, t310, t308, t306, t184, t3 * t175 + t307, t2 * t175 - t200 * t4 - t203 * t249 + t309;];
tauc_reg  = t1;
