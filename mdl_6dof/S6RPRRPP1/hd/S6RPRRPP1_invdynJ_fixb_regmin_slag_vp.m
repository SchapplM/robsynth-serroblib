% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% 
% Output:
% tau_reg [6x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:29:54
% EndTime: 2019-03-09 04:30:05
% DurationCPUTime: 4.61s
% Computational Cost: add. (5822->462), mult. (12151->595), div. (0->0), fcn. (7987->14), ass. (0->227)
t204 = cos(qJ(3));
t265 = t204 * qJD(1);
t166 = -qJD(4) + t265;
t190 = g(3) * t204;
t201 = sin(qJ(3));
t192 = qJ(1) + pkin(9);
t181 = sin(t192);
t183 = cos(t192);
t234 = g(1) * t183 + g(2) * t181;
t328 = -t234 * t201 + t190;
t196 = sin(pkin(9));
t174 = t196 * pkin(1) + pkin(7);
t150 = t174 * qJD(1);
t272 = qJD(3) * t204;
t148 = t174 * qJDD(1);
t334 = -qJD(2) * qJD(3) - t148;
t226 = t204 * qJDD(2) - t150 * t272 + t334 * t201;
t62 = -qJDD(3) * pkin(3) - t226;
t336 = -qJD(4) * pkin(8) * t166 + t328 + t62;
t200 = sin(qJ(4));
t203 = cos(qJ(4));
t267 = t203 * qJD(3);
t276 = qJD(1) * t201;
t131 = t200 * t276 - t267;
t268 = t200 * qJD(3);
t133 = t203 * t276 + t268;
t195 = sin(pkin(10));
t197 = cos(pkin(10));
t79 = t197 * t131 + t195 * t133;
t335 = t166 * t79;
t126 = t195 * t203 + t197 * t200;
t113 = t126 * qJD(4);
t304 = -t126 * t265 + t113;
t269 = qJD(4) * t203;
t271 = qJD(4) * t200;
t114 = -t195 * t271 + t197 * t269;
t251 = t200 * t265;
t292 = t197 * t203;
t303 = t195 * t251 - t265 * t292 + t114;
t179 = t203 * pkin(4) + pkin(3);
t312 = qJ(5) + pkin(8);
t238 = t204 * t179 + t201 * t312;
t105 = t204 * qJD(2) - t201 * t150;
t333 = t105 * qJD(3);
t270 = qJD(4) * t201;
t332 = -qJD(1) * t270 + qJDD(3);
t331 = -pkin(2) - t238;
t261 = t201 * qJDD(1);
t75 = (qJD(3) * (qJD(4) + t265) + t261) * t200 - t332 * t203;
t230 = -t195 * t131 + t197 * t133;
t330 = t230 ^ 2;
t242 = qJD(4) * t312;
t266 = t203 * qJD(5);
t110 = -t200 * t242 + t266;
t219 = -t200 * qJD(5) - t203 * t242;
t235 = pkin(3) * t201 - pkin(8) * t204;
t135 = t235 * qJD(1);
t117 = t203 * t135;
t285 = t203 * t204;
t320 = pkin(4) * t201;
t228 = -qJ(5) * t285 + t320;
t54 = qJD(1) * t228 - t200 * t105 + t117;
t302 = t203 * t105 + t200 * t135;
t58 = -qJ(5) * t251 + t302;
t310 = (t219 - t54) * t197 + (-t110 + t58) * t195;
t106 = t201 * qJD(2) + t204 * t150;
t327 = -t106 + (-t251 + t271) * pkin(4);
t246 = t204 * t267;
t74 = qJD(1) * t246 + qJD(4) * t267 + t332 * t200 + t203 * t261;
t40 = t195 * t74 + t197 * t75;
t41 = -t195 * t75 + t197 * t74;
t326 = t40 * pkin(5) - t41 * qJ(6) - t230 * qJD(6);
t191 = qJ(4) + pkin(10);
t180 = sin(t191);
t61 = qJDD(3) * pkin(8) + t201 * qJDD(2) + t204 * t148 + t333;
t198 = cos(pkin(9));
t176 = -t198 * pkin(1) - pkin(2);
t120 = -t204 * pkin(3) - t201 * pkin(8) + t176;
t136 = t235 * qJD(3);
t77 = qJD(1) * t136 + qJDD(1) * t120;
t96 = t120 * qJD(1);
t259 = t200 * t77 + t203 * t61 + t96 * t269;
t93 = qJD(3) * pkin(8) + t106;
t224 = -t271 * t93 + t259;
t13 = -t75 * qJ(5) - t131 * qJD(5) + t224;
t186 = t204 * qJDD(1);
t263 = qJD(1) * qJD(3);
t243 = t201 * t263;
t124 = qJDD(4) - t186 + t243;
t57 = t200 * t96 + t203 * t93;
t71 = t203 * t77;
t8 = t124 * pkin(4) - t74 * qJ(5) - qJD(4) * t57 - t133 * qJD(5) - t200 * t61 + t71;
t3 = -t195 * t13 + t197 * t8;
t256 = -qJDD(6) + t3;
t92 = -qJD(3) * pkin(3) - t105;
t73 = t131 * pkin(4) + qJD(5) + t92;
t26 = t79 * pkin(5) - qJ(6) * t230 + t73;
t315 = g(3) * t201;
t182 = cos(t191);
t296 = t181 * t204;
t86 = t180 * t296 + t183 * t182;
t293 = t183 * t204;
t88 = t180 * t293 - t181 * t182;
t325 = g(1) * t88 + g(2) * t86 + t180 * t315 - t26 * t230 + t256;
t220 = -t200 * t270 + t246;
t286 = t203 * t124;
t324 = -t166 * t220 + t201 * t286;
t323 = t204 * t234 + t315;
t4 = t197 * t13 + t195 * t8;
t319 = g(1) * t181;
t316 = g(2) * t183;
t314 = t124 * pkin(5);
t50 = -t131 * qJ(5) + t57;
t47 = t197 * t50;
t56 = -t200 * t93 + t203 * t96;
t49 = -t133 * qJ(5) + t56;
t20 = t195 * t49 + t47;
t313 = t20 * t230;
t134 = t174 * t285;
t153 = t201 * t174;
t280 = t203 * t136 + t268 * t153;
t32 = -t201 * t266 + t228 * qJD(3) + (-t134 + (qJ(5) * t201 - t120) * t200) * qJD(4) + t280;
t274 = qJD(3) * t174;
t281 = t120 * t269 + t200 * t136;
t287 = t201 * t203;
t37 = (-qJ(5) * qJD(4) - t274) * t287 + (-qJD(5) * t201 + (-qJ(5) * qJD(3) - qJD(4) * t174) * t204) * t200 + t281;
t15 = t195 * t32 + t197 * t37;
t46 = -t166 * pkin(4) + t49;
t19 = t195 * t46 + t47;
t311 = t304 * pkin(5) - t303 * qJ(6) - t126 * qJD(6) + t327;
t28 = t195 * t54 + t197 * t58;
t309 = pkin(5) * t276 - t310;
t23 = qJ(6) * t276 + t28;
t67 = t197 * t110 + t195 * t219;
t308 = t67 - t23;
t108 = t203 * t120;
t65 = -qJ(5) * t287 + t108 + (-t174 * t200 - pkin(4)) * t204;
t279 = t200 * t120 + t134;
t290 = t200 * t201;
t72 = -qJ(5) * t290 + t279;
t36 = t195 * t65 + t197 * t72;
t307 = t195 * t50;
t305 = t74 * t200;
t301 = t131 * t166;
t300 = t133 * t166;
t299 = t166 * t203;
t298 = t181 * t200;
t297 = t181 * t203;
t295 = t183 * t200;
t294 = t183 * t203;
t291 = t312 * t204;
t289 = t200 * t204;
t21 = t197 * t49 - t307;
t283 = qJD(6) - t21;
t282 = qJDD(2) - g(3);
t168 = pkin(4) * t290;
t278 = t153 + t168;
t193 = t201 ^ 2;
t277 = -t204 ^ 2 + t193;
t151 = qJD(1) * t176;
t275 = qJD(3) * t131;
t273 = qJD(3) * t201;
t260 = t124 * qJ(6) + t4;
t255 = t183 * t289;
t249 = t204 * t268;
t252 = pkin(4) * t249 + t174 * t272 + t269 * t320;
t250 = t166 * t268;
t245 = t312 * t200;
t241 = -qJD(4) * t96 - t61;
t240 = t133 * t273 - t74 * t204;
t239 = t166 * t174 + t93;
t157 = t201 * t319;
t236 = -t201 * t316 + t157;
t202 = sin(qJ(1));
t205 = cos(qJ(1));
t233 = g(1) * t202 - g(2) * t205;
t232 = -t79 ^ 2 - t330;
t231 = pkin(5) * t182 + qJ(6) * t180;
t14 = -t195 * t37 + t197 * t32;
t18 = t197 * t46 - t307;
t35 = -t195 * t72 + t197 * t65;
t98 = t181 * t289 + t294;
t223 = -t200 * t124 + t166 * t269;
t103 = t126 * t201;
t104 = -t195 * t290 + t197 * t287;
t63 = -t114 * t201 - t195 * t246 - t197 * t249;
t64 = t113 * t201 + t195 * t249 - t197 * t246;
t222 = t103 * t41 - t104 * t40 - t230 * t63 + t64 * t79;
t221 = -qJD(1) * t151 + t234;
t218 = -pkin(8) * t124 - t166 * t92;
t217 = t205 * pkin(1) + pkin(4) * t298 + t181 * pkin(7) - t331 * t183;
t206 = qJD(3) ^ 2;
t216 = 0.2e1 * qJDD(1) * t176 + t174 * t206 + t316;
t215 = 0.2e1 * qJD(3) * t151 - qJDD(3) * t174;
t211 = -t202 * pkin(1) + pkin(4) * t295 + t183 * pkin(7) + t331 * t181;
t39 = t75 * pkin(4) + qJDD(5) + t62;
t152 = t312 * t203;
t83 = t195 * t152 + t197 * t245;
t84 = t197 * t152 - t195 * t245;
t210 = -t84 * t40 + t83 * t41 - t67 * t79 - t323;
t209 = t39 + t328;
t207 = qJD(1) ^ 2;
t175 = -t197 * pkin(4) - pkin(5);
t169 = t195 * pkin(4) + qJ(6);
t155 = pkin(4) * t297;
t142 = qJDD(3) * t204 - t206 * t201;
t141 = qJDD(3) * t201 + t206 * t204;
t125 = t195 * t200 - t292;
t101 = t183 * t285 + t298;
t100 = -t255 + t297;
t99 = -t181 * t285 + t295;
t89 = t181 * t180 + t182 * t293;
t87 = -t183 * t180 + t182 * t296;
t76 = t125 * pkin(5) - t126 * qJ(6) - t179;
t51 = t103 * pkin(5) - t104 * qJ(6) + t278;
t38 = t133 * pkin(4) + pkin(5) * t230 + qJ(6) * t79;
t33 = t204 * pkin(5) - t35;
t31 = -t204 * qJ(6) + t36;
t22 = -t63 * pkin(5) + t64 * qJ(6) - t104 * qJD(6) + t252;
t17 = -t166 * qJ(6) + t19;
t16 = t166 * pkin(5) + qJD(6) - t18;
t12 = -pkin(5) * t273 - t14;
t9 = qJ(6) * t273 - t204 * qJD(6) + t15;
t5 = t39 + t326;
t2 = -t256 - t314;
t1 = -t166 * qJD(6) + t260;
t6 = [qJDD(1), t233, g(1) * t205 + g(2) * t202 (t233 + (t196 ^ 2 + t198 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t193 * qJDD(1) + 0.2e1 * t204 * t243, 0.2e1 * t186 * t201 - 0.2e1 * t263 * t277, t141, t142, 0, t215 * t201 + (-t216 + t319) * t204, t201 * t216 + t204 * t215 - t157, t133 * t220 + t287 * t74 (-t131 * t203 - t133 * t200) * t272 + (-t305 - t203 * t75 + (t131 * t200 - t133 * t203) * qJD(4)) * t201, t240 + t324 (t75 + t250) * t204 + (t223 - t275) * t201, -t124 * t204 - t166 * t273 -(-t120 * t271 + t280) * t166 + t108 * t124 - g(1) * t99 - g(2) * t101 + (t131 * t274 - t71 + t239 * t269 + (qJD(3) * t92 - t124 * t174 - t241) * t200) * t204 + (t56 * qJD(3) + t174 * t75 + t62 * t200 + t269 * t92) * t201, t281 * t166 - t279 * t124 - g(1) * t98 - g(2) * t100 + (-t239 * t271 + (t133 * t174 + t203 * t92) * qJD(3) + t259) * t204 + (-t92 * t271 + t174 * t74 + t62 * t203 + (-t174 * t299 - t57) * qJD(3)) * t201, -t4 * t103 - t3 * t104 - t14 * t230 - t15 * t79 + t18 * t64 + t19 * t63 - t35 * t41 - t36 * t40 + t236, -g(1) * t211 - g(2) * t217 + t18 * t14 + t19 * t15 + t252 * t73 + t278 * t39 + t3 * t35 + t4 * t36, g(1) * t87 - g(2) * t89 + t5 * t103 + t12 * t166 - t33 * t124 - t16 * t273 + t2 * t204 + t22 * t79 - t26 * t63 + t51 * t40, -t1 * t103 + t2 * t104 + t12 * t230 - t16 * t64 + t17 * t63 - t31 * t40 + t33 * t41 - t9 * t79 + t236, g(1) * t86 - g(2) * t88 - t1 * t204 - t5 * t104 + t31 * t124 - t9 * t166 + t17 * t273 - t22 * t230 + t26 * t64 - t51 * t41, t1 * t31 + t17 * t9 + t5 * t51 + t26 * t22 + t2 * t33 + t16 * t12 - g(1) * (-t87 * pkin(5) - t86 * qJ(6) + t211) - g(2) * (t89 * pkin(5) + t88 * qJ(6) + t217); 0, 0, 0, t282, 0, 0, 0, 0, 0, t142, -t141, 0, 0, 0, 0, 0 (-t75 + t250) * t204 + (t223 + t275) * t201, t240 - t324, t222, -t3 * t103 + t4 * t104 + t18 * t63 - t19 * t64 - t39 * t204 + t273 * t73 - g(3), -t103 * t124 - t63 * t166 - t204 * t40 + t273 * t79, t222, t104 * t124 + t64 * t166 + t204 * t41 - t230 * t273, t1 * t104 + t2 * t103 - t16 * t63 - t17 * t64 - t5 * t204 + t26 * t273 - g(3); 0, 0, 0, 0, -t204 * t207 * t201, t277 * t207, t261, t186, qJDD(3), t106 * qJD(3) + t201 * t221 - t190 + t226, t333 + (qJD(3) * t150 - t282) * t201 + (t221 + t334) * t204, -t133 * t299 + t305 (t74 + t301) * t203 + (-t75 + t300) * t200 (-t133 * t201 + t166 * t285) * qJD(1) - t223, t166 * t271 + t286 + (t131 * t201 - t166 * t289) * qJD(1), t166 * t276, -t56 * t276 - pkin(3) * t75 - t106 * t131 + t117 * t166 + (-t105 * t166 + t218) * t200 - t336 * t203, -pkin(3) * t74 - t106 * t133 - t302 * t166 + t200 * t336 + t218 * t203 + t57 * t276, -t4 * t125 - t3 * t126 - t303 * t18 - t304 * t19 - t230 * t310 + t28 * t79 + t210, t4 * t84 - t3 * t83 - t39 * t179 - g(3) * t238 + t327 * t73 + (t67 - t28) * t19 + t310 * t18 + t234 * (t179 * t201 - t291) -t83 * t124 + t5 * t125 + t16 * t276 + t309 * t166 - t182 * t328 + t304 * t26 + t311 * t79 + t76 * t40, -t1 * t125 + t2 * t126 + t303 * t16 - t304 * t17 + t23 * t79 + t230 * t309 + t210, t84 * t124 - t5 * t126 - t308 * t166 - t17 * t276 - t180 * t328 - t230 * t311 - t303 * t26 - t76 * t41, t1 * t84 + t5 * t76 + t2 * t83 - g(3) * (t204 * t231 + t238) + t311 * t26 + t308 * t17 + t309 * t16 + t234 * (-t291 - (-t179 - t231) * t201); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133 * t131, -t131 ^ 2 + t133 ^ 2, t74 - t301, -t300 - t75, t124, -t93 * t269 - g(1) * t100 + g(2) * t98 - t92 * t133 - t57 * t166 + t71 + (t241 + t315) * t200, g(1) * t101 - g(2) * t99 + g(3) * t287 + t92 * t131 - t56 * t166 - t224, t19 * t230 - t313 + (-t195 * t40 - t197 * t41) * pkin(4) + (-t18 + t21) * t79, -g(1) * t155 + t18 * t20 - t19 * t21 + (g(2) * t294 - t73 * t133 + t4 * t195 + t3 * t197 + t323 * t200) * pkin(4), -t20 * t166 - t38 * t79 + (pkin(5) - t175) * t124 + t325, -t169 * t40 + t17 * t230 + t175 * t41 - t313 + (t16 - t283) * t79, -t182 * t315 - g(1) * t89 - g(2) * t87 + t169 * t124 - t26 * t79 + t38 * t230 + (-0.2e1 * qJD(6) + t21) * t166 + t260, t1 * t169 + t2 * t175 - t26 * t38 - t16 * t20 - g(1) * (-pkin(4) * t255 - t88 * pkin(5) + t89 * qJ(6) + t155) - g(2) * (-pkin(4) * t98 - t86 * pkin(5) + t87 * qJ(6)) - g(3) * (-t168 + (-pkin(5) * t180 + qJ(6) * t182) * t201) + t283 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t232, t18 * t230 + t19 * t79 + t209, -t166 * t230 + t40, t232, -t41 - t335, -t16 * t230 + t17 * t79 + t209 + t326; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230 * t79 - t124, t41 - t335, -t166 ^ 2 - t330, t17 * t166 - t314 - t325;];
tau_reg  = t6;
