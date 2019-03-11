% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6PRPPRR1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:16
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPPRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPPRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:16:21
% EndTime: 2019-03-08 19:16:27
% DurationCPUTime: 4.62s
% Computational Cost: add. (6013->446), mult. (13993->598), div. (0->0), fcn. (11880->16), ass. (0->235)
t179 = sin(pkin(11));
t189 = sin(qJ(2));
t181 = sin(pkin(6));
t267 = qJD(1) * t181;
t246 = t189 * t267;
t146 = t179 * t246;
t183 = cos(pkin(11));
t191 = cos(qJ(2));
t245 = t191 * t267;
t115 = t183 * t245 - t146;
t258 = qJD(4) - t115;
t182 = cos(pkin(12));
t315 = cos(qJ(5));
t247 = t315 * t182;
t178 = sin(pkin(12));
t188 = sin(qJ(5));
t272 = t188 * t178;
t211 = t247 - t272;
t162 = pkin(2) * t179 + qJ(4);
t312 = pkin(8) + t162;
t128 = t312 * t178;
t129 = t312 * t182;
t212 = -t315 * t128 - t188 * t129;
t308 = qJD(5) * t212 + t211 * t258;
t219 = t179 * t191 + t183 * t189;
t118 = t219 * t181;
t112 = qJD(1) * t118;
t242 = qJD(5) * t315;
t263 = qJD(5) * t188;
t243 = t178 * t263;
t126 = -t182 * t242 + t243;
t133 = t315 * t178 + t188 * t182;
t127 = t133 * qJD(5);
t322 = pkin(5) * t127 + pkin(9) * t126 - t112;
t264 = qJD(2) * t181;
t241 = qJD(1) * t264;
t257 = qJDD(1) * t181;
t321 = t189 * t257 + t191 * t241;
t165 = pkin(4) * t182 + pkin(3);
t225 = t165 * qJDD(2);
t153 = t191 * t257;
t297 = qJDD(2) * pkin(2);
t116 = -t189 * t241 + t153 + t297;
t249 = -t183 * t116 + t179 * t321;
t233 = qJDD(4) + t249;
t44 = -t225 + t233;
t156 = qJD(2) * t247;
t238 = qJDD(2) * t315;
t255 = t182 * qJDD(2);
t248 = qJD(5) * t156 + t178 * t238 + t188 * t255;
t86 = qJD(2) * t243 - t248;
t256 = t178 * qJDD(2);
t224 = -t182 * t238 + t188 * t256;
t87 = qJD(2) * t127 + t224;
t16 = pkin(5) * t87 + pkin(9) * t86 + t44;
t187 = sin(qJ(6));
t190 = cos(qJ(6));
t185 = cos(pkin(6));
t160 = qJD(1) * t185 + qJD(3);
t144 = t182 * t160;
t303 = pkin(8) * qJD(2);
t138 = qJD(2) * pkin(2) + t245;
t101 = t179 * t138 + t183 * t246;
t98 = qJD(2) * qJ(4) + t101;
t57 = t144 + (-t98 - t303) * t178;
t65 = t178 * t160 + t182 * t98;
t58 = t182 * t303 + t65;
t22 = t188 * t57 + t315 * t58;
t20 = qJD(5) * pkin(9) + t22;
t244 = qJD(2) * t272;
t122 = -t156 + t244;
t124 = t133 * qJD(2);
t100 = t138 * t183 - t146;
t226 = qJD(4) - t100;
t82 = -t165 * qJD(2) + t226;
t37 = pkin(5) * t122 - pkin(9) * t124 + t82;
t222 = t187 * t20 - t190 * t37;
t158 = t185 * qJDD(1) + qJDD(3);
t140 = t182 * t158;
t61 = t179 * t116 + t183 * t321;
t55 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t61;
t34 = t140 + (-pkin(8) * qJDD(2) - t55) * t178;
t39 = t178 * t158 + t182 * t55;
t35 = pkin(8) * t255 + t39;
t254 = -t188 * t34 - t57 * t242 - t315 * t35;
t5 = -t58 * t263 - t254;
t3 = qJDD(5) * pkin(9) + t5;
t1 = -t222 * qJD(6) + t187 * t16 + t190 * t3;
t120 = qJD(6) + t122;
t230 = t120 * t222 + t1;
t10 = t187 * t37 + t190 * t20;
t2 = -qJD(6) * t10 + t190 * t16 - t187 * t3;
t320 = t10 * t120 + t2;
t105 = qJD(5) * t187 + t124 * t190;
t235 = t120 * t187;
t319 = t105 * t235;
t177 = pkin(12) + qJ(5);
t169 = sin(t177);
t271 = t191 * t183;
t277 = t181 * t189;
t117 = t179 * t277 - t181 * t271;
t180 = sin(pkin(10));
t184 = cos(pkin(10));
t132 = t179 * t189 - t271;
t203 = t132 * t185;
t71 = -t180 * t219 - t184 * t203;
t74 = t180 * t203 - t184 * t219;
t206 = g(1) * t74 + g(2) * t71 - g(3) * t117;
t318 = t206 * t169;
t111 = t117 * t190;
t96 = -t118 * t178 + t182 * t185;
t97 = t118 * t182 + t178 * t185;
t43 = t188 * t96 + t315 * t97;
t23 = -t187 * t43 + t111;
t262 = qJD(6) * t105;
t41 = -t190 * qJDD(5) - t187 * t86 + t262;
t239 = t188 * t35 - t315 * t34;
t6 = -qJD(5) * t22 - t239;
t261 = qJD(6) * t187;
t286 = t126 * t190;
t209 = t133 * t261 + t286;
t284 = t133 * t190;
t83 = qJDD(6) + t87;
t317 = -t120 * t209 + t83 * t284;
t273 = t185 * t191;
t274 = t185 * t189;
t269 = -t179 * t273 - t183 * t274;
t75 = -t184 * t132 + t180 * t269;
t70 = t180 * t132 + t184 * t269;
t316 = t124 ^ 2;
t314 = pkin(2) * t183;
t145 = -t165 - t314;
t76 = -pkin(5) * t211 - pkin(9) * t133 + t145;
t79 = -t188 * t128 + t315 * t129;
t28 = -t187 * t79 + t190 * t76;
t311 = qJD(6) * t28 + t322 * t187 + t308 * t190;
t29 = t187 * t76 + t190 * t79;
t310 = -qJD(6) * t29 - t308 * t187 + t322 * t190;
t259 = t190 * qJD(5);
t103 = t124 * t187 - t259;
t260 = qJD(6) * t190;
t309 = -t103 * t260 - t187 * t41;
t307 = qJD(5) * t79 + t133 * t258;
t306 = t103 * t286 - t41 * t284;
t40 = -qJD(6) * t259 - t187 * qJDD(5) + t124 * t261 + t190 * t86;
t305 = t105 * t127 + t211 * t40;
t304 = t126 * t122 - t133 * t87;
t301 = t187 * t40;
t299 = t187 * t83;
t298 = t188 * t58;
t296 = qJDD(2) * pkin(3);
t295 = t103 * t122;
t294 = t103 * t124;
t293 = t103 * t187;
t292 = t105 * t103;
t291 = t105 * t124;
t290 = t105 * t190;
t289 = t117 * t187;
t288 = t124 * t122;
t287 = t126 * t187;
t285 = t133 * t187;
t170 = cos(t177);
t283 = t170 * t187;
t282 = t170 * t190;
t280 = t180 * t181;
t279 = t180 * t189;
t278 = t181 * t184;
t276 = t181 * t191;
t270 = qJDD(1) - g(3);
t174 = t178 ^ 2;
t176 = t182 ^ 2;
t268 = t174 + t176;
t266 = qJD(2) * t112;
t114 = t132 * t264;
t265 = qJD(2) * t114;
t252 = t105 * t287;
t251 = t184 * t273;
t161 = pkin(2) * t276;
t186 = -pkin(8) - qJ(4);
t250 = -t117 * t165 - t118 * t186 + t161;
t236 = t120 * t190;
t147 = pkin(2) * t251;
t232 = -pkin(2) * t279 + t147;
t231 = pkin(5) * t170 + pkin(9) * t169;
t229 = -t10 * t190 - t187 * t222;
t228 = t10 * t187 - t190 * t222;
t64 = -t178 * t98 + t144;
t223 = t178 * t64 - t182 * t65;
t221 = -t103 * t127 + t211 * t41;
t24 = t190 * t43 + t289;
t220 = t124 * t127 + t211 * t86;
t218 = -t120 * t261 - t122 * t235 + t190 * t83;
t113 = qJD(2) * t118;
t217 = qJD(2) * t113 + qJDD(2) * t117;
t21 = t315 * t57 - t298;
t216 = -t188 * t97 + t315 * t96;
t215 = t71 * t165 + t70 * t186 + t232;
t214 = -t180 * t273 - t184 * t189;
t19 = -qJD(5) * pkin(5) - t21;
t213 = -pkin(9) * t83 + t120 * t19;
t210 = t133 * t260 - t287;
t208 = -g(1) * (-t169 * t75 + t170 * t280) - g(2) * (t169 * t70 - t170 * t278) - g(3) * (-t118 * t169 + t170 * t185);
t47 = -t169 * t278 - t170 * t70;
t49 = t169 * t280 + t170 * t75;
t91 = t118 * t170 + t169 * t185;
t207 = -g(1) * t49 - g(2) * t47 - g(3) * t91;
t205 = -g(1) * t75 + g(2) * t70 - g(3) * t118;
t204 = -g(1) * t280 + g(2) * t278 - g(3) * t185;
t4 = -qJDD(5) * pkin(5) - t6;
t202 = t208 - t4;
t201 = t214 * pkin(2);
t200 = t206 + t249;
t199 = t74 * t165 - t186 * t75 + t201;
t198 = pkin(9) * qJD(6) * t120 - t202;
t197 = qJDD(4) + t200;
t196 = -g(1) * t214 - g(3) * t276;
t195 = -t228 * qJD(6) + t1 * t190 - t2 * t187;
t194 = -t120 * t210 - t83 * t285;
t166 = -pkin(3) - t314;
t56 = t233 - t296;
t193 = qJDD(2) * t166 + t206 - t266 + t56;
t192 = qJD(2) ^ 2;
t121 = t122 ^ 2;
t95 = -qJD(2) * pkin(3) + t226;
t89 = -qJD(5) * t127 + qJDD(5) * t211;
t88 = -qJD(5) * t126 + qJDD(5) * t133;
t84 = pkin(5) * t124 + pkin(9) * t122;
t38 = -t178 * t55 + t140;
t18 = qJD(5) * t43 - t114 * t133;
t17 = qJD(5) * t216 - t114 * t211;
t14 = t187 * t84 + t190 * t21;
t13 = -t187 * t21 + t190 * t84;
t8 = -qJD(6) * t24 + t113 * t190 - t17 * t187;
t7 = t23 * qJD(6) + t113 * t187 + t17 * t190;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t270, 0, 0, 0, 0, 0, 0 (qJDD(2) * t191 - t189 * t192) * t181 (-qJDD(2) * t189 - t191 * t192) * t181, 0, -g(3) + (t185 ^ 2 + (t189 ^ 2 + t191 ^ 2) * t181 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -t217, -qJDD(2) * t118 + t265, 0, -t100 * t113 - t101 * t114 + t117 * t249 + t118 * t61 + t158 * t185 - g(3), 0, 0, 0, 0, 0, 0, -t217 * t182, t217 * t178 (-t178 * t96 + t182 * t97) * qJDD(2) - t268 * t265, t113 * t95 + t114 * t223 + t117 * t56 + t38 * t96 + t39 * t97 - g(3), 0, 0, 0, 0, 0, 0, -qJD(5) * t18 + qJDD(5) * t216 + t113 * t122 + t117 * t87, -qJD(5) * t17 - qJDD(5) * t43 + t113 * t124 - t117 * t86, -t122 * t17 + t124 * t18 + t216 * t86 - t43 * t87, t113 * t82 + t117 * t44 + t17 * t22 - t18 * t21 + t216 * t6 + t43 * t5 - g(3), 0, 0, 0, 0, 0, 0, t103 * t18 + t120 * t8 - t216 * t41 + t23 * t83, t105 * t18 - t120 * t7 + t216 * t40 - t24 * t83, -t103 * t7 - t105 * t8 + t23 * t40 - t24 * t41, t1 * t24 + t10 * t7 + t18 * t19 + t2 * t23 - t216 * t4 - t222 * t8 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t153 - g(2) * (t251 - t279) + t196, -g(1) * (t180 * t274 - t184 * t191) - g(2) * (-t180 * t191 - t184 * t274) - t270 * t277, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t183 * t297 - t200 + t266, qJD(2) * t115 - t179 * t297 - t205 - t61, 0, -g(2) * t147 + t100 * t112 - t101 * t115 + (g(2) * t279 + t61 * t179 - t183 * t249 + t196) * pkin(2), t174 * qJDD(2), 0.2e1 * t178 * t255, 0, t176 * qJDD(2), 0, 0, -t193 * t182, t193 * t178, -t178 * t38 + t182 * t39 + t205 + (t258 * qJD(2) + t162 * qJDD(2)) * t268, t56 * t166 - t95 * t112 - g(1) * (pkin(3) * t74 + qJ(4) * t75 + t201) - g(2) * (pkin(3) * t71 - qJ(4) * t70 + t232) - g(3) * (-pkin(3) * t117 + qJ(4) * t118 + t161) + (t39 * t162 + t258 * t65) * t182 + (-t38 * t162 - t258 * t64) * t178, -t124 * t126 - t133 * t86, -t220 + t304, t88, t122 * t127 - t211 * t87, t89, 0, -t307 * qJD(5) + qJDD(5) * t212 - t112 * t122 + t127 * t82 + t145 * t87 - t206 * t170 - t211 * t44, -t308 * qJD(5) - qJDD(5) * t79 - t112 * t124 - t126 * t82 + t133 * t44 - t145 * t86 + t318, -t308 * t122 + t307 * t124 + t126 * t21 - t127 * t22 - t133 * t6 + t211 * t5 + t212 * t86 - t79 * t87 + t205, -g(1) * t199 - g(2) * t215 - g(3) * t250 - t82 * t112 + t44 * t145 - t307 * t21 + t212 * t6 + t308 * t22 + t5 * t79, -t105 * t209 - t284 * t40, t252 + (t301 + (-t290 + t293) * qJD(6)) * t133 + t306, t305 + t317, t103 * t210 + t285 * t41, t194 + t221, t120 * t127 - t211 * t83, t28 * t83 - t2 * t211 - t222 * t127 - t212 * t41 - t19 * t287 - g(1) * (t187 * t75 + t282 * t74) - g(2) * (-t187 * t70 + t282 * t71) - g(3) * (-t117 * t282 + t118 * t187) + (t187 * t4 + t19 * t260) * t133 + t310 * t120 + t307 * t103, -t29 * t83 + t1 * t211 - t10 * t127 + t212 * t40 - t19 * t286 - g(1) * (t190 * t75 - t283 * t74) - g(2) * (-t190 * t70 - t283 * t71) - g(3) * (t117 * t283 + t118 * t190) + (-t19 * t261 + t190 * t4) * t133 - t311 * t120 + t307 * t105, t28 * t40 - t29 * t41 + t228 * t126 - t310 * t105 - t311 * t103 - t318 + (qJD(6) * t229 - t1 * t187 - t190 * t2) * t133, t1 * t29 + t2 * t28 - t4 * t212 - g(1) * (t231 * t74 + t199) - g(2) * (t231 * t71 + t215) - g(3) * (-t117 * t231 + t250) - t310 * t222 + t307 * t19 + t311 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204 + t158, 0, 0, 0, 0, 0, 0, 0, 0, 0, t178 * t39 + t182 * t38 + t204, 0, 0, 0, 0, 0, 0, t89, -t88, t220 + t304, -t126 * t22 - t127 * t21 + t133 * t5 + t211 * t6 + t204, 0, 0, 0, 0, 0, 0, t194 - t221, t305 - t317, -t252 + (-t301 + (t290 + t293) * qJD(6)) * t133 + t306, t126 * t229 + t127 * t19 + t133 * t195 - t211 * t4 + t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t255, t256, -t268 * t192, qJD(2) * t223 + t197 - t296, 0, 0, 0, 0, 0, 0, 0.2e1 * t124 * qJD(5) + t224 (-t122 - t244) * qJD(5) + t248, -t121 - t316, t122 * t22 + t124 * t21 + t197 - t225, 0, 0, 0, 0, 0, 0, t218 - t294, -t120 ^ 2 * t190 - t291 - t299 (t40 - t295) * t190 + t319 + t309, -t124 * t19 + t230 * t187 + t320 * t190 + t206; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t288, -t121 + t316 (t122 - t244) * qJD(5) + t248, -t288, -t224, qJDD(5), -t82 * t124 + t208 - t239, t122 * t82 + (t21 + t298) * qJD(5) - t207 + t254, 0, 0, t105 * t236 - t301 (-t40 - t295) * t190 - t319 + t309, t120 * t236 - t291 + t299, t103 * t235 - t190 * t41, t218 + t294, -t120 * t124, -pkin(5) * t41 - t103 * t22 - t120 * t13 + t124 * t222 + t187 * t213 - t190 * t198, pkin(5) * t40 + t10 * t124 - t105 * t22 + t120 * t14 + t187 * t198 + t190 * t213, t103 * t14 + t105 * t13 + ((-t41 + t262) * pkin(9) + t230) * t190 + ((qJD(6) * t103 - t40) * pkin(9) - t320) * t187 + t207, -t10 * t14 + t222 * t13 - t19 * t22 + t202 * pkin(5) + (t195 + t207) * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t292, -t103 ^ 2 + t105 ^ 2, t103 * t120 - t40, -t292, t105 * t120 - t41, t83, -t19 * t105 - g(1) * (-t187 * t49 - t190 * t74) - g(2) * (-t187 * t47 - t190 * t71) - g(3) * (-t187 * t91 + t111) + t320, t19 * t103 - g(1) * (t187 * t74 - t190 * t49) - g(2) * (t187 * t71 - t190 * t47) - g(3) * (-t190 * t91 - t289) - t230, 0, 0;];
tau_reg  = t9;
