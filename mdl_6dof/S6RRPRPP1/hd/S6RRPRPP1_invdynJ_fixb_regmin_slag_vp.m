% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
% 
% Output:
% tau_reg [6x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:48
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPRPP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPRPP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRPP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:47:55
% EndTime: 2019-03-09 09:48:07
% DurationCPUTime: 5.01s
% Computational Cost: add. (9779->502), mult. (22730->637), div. (0->0), fcn. (16735->14), ass. (0->246)
t218 = sin(pkin(9));
t224 = sin(qJ(2));
t289 = qJD(1) * t224;
t220 = cos(pkin(9));
t227 = cos(qJ(2));
t303 = t220 * t227;
t158 = qJD(1) * t303 - t218 * t289;
t149 = qJD(4) - t158;
t198 = pkin(2) * t218 + pkin(8);
t214 = qJ(2) + pkin(9);
t205 = sin(t214);
t207 = cos(t214);
t225 = sin(qJ(1));
t228 = cos(qJ(1));
t262 = g(1) * t228 + g(2) * t225;
t346 = g(3) * t207 - t262 * t205;
t330 = qJ(3) + pkin(7);
t270 = qJD(2) * t330;
t156 = -t224 * qJD(3) - t227 * t270;
t180 = t330 * t224;
t115 = qJDD(2) * pkin(2) + t156 * qJD(1) - qJDD(1) * t180;
t155 = t227 * qJD(3) - t224 * t270;
t181 = t330 * t227;
t123 = t155 * qJD(1) + qJDD(1) * t181;
t69 = t220 * t115 - t218 * t123;
t67 = -qJDD(2) * pkin(3) - t69;
t352 = qJD(4) * t198 * t149 + t346 + t67;
t173 = t218 * t227 + t220 * t224;
t217 = sin(pkin(10));
t219 = cos(pkin(10));
t223 = sin(qJ(4));
t226 = cos(qJ(4));
t344 = -t217 * t223 + t219 * t226;
t101 = t344 * t173;
t288 = qJD(4) * t223;
t314 = t158 * t223;
t351 = t288 - t314;
t161 = t173 * qJD(1);
t286 = t226 * qJD(2);
t132 = t161 * t223 - t286;
t134 = qJD(2) * t223 + t161 * t226;
t77 = t219 * t132 + t134 * t217;
t350 = t149 * t77;
t172 = t217 * t226 + t219 * t223;
t157 = t172 * qJD(4);
t320 = -t172 * t158 + t157;
t287 = qJD(4) * t226;
t319 = -t344 * t158 - t217 * t288 + t219 * t287;
t171 = t218 * t224 - t303;
t163 = t171 * qJD(2);
t275 = t173 * t287;
t349 = -t163 * t223 + t275;
t253 = -t132 * t217 + t219 * t134;
t348 = t253 ^ 2;
t211 = t227 * pkin(2);
t203 = t211 + pkin(1);
t293 = qJ(5) + t198;
t265 = qJD(4) * t293;
t137 = t226 * qJD(5) - t223 * t265;
t241 = -t223 * qJD(5) - t226 * t265;
t176 = qJD(1) * t181;
t164 = t218 * t176;
t175 = qJD(1) * t180;
t125 = -t175 * t220 - t164;
t96 = pkin(2) * t289 + pkin(3) * t161 - pkin(8) * t158;
t92 = t226 * t96;
t49 = -qJ(5) * t158 * t226 + pkin(4) * t161 - t125 * t223 + t92;
t321 = t226 * t125 + t223 * t96;
t56 = -qJ(5) * t314 + t321;
t328 = (t241 - t49) * t219 + (-t137 + t56) * t217;
t271 = g(1) * t225 - g(2) * t228;
t304 = t220 * t176;
t124 = -t175 * t218 + t304;
t345 = t351 * pkin(4) - t124;
t285 = qJD(1) * qJD(2);
t272 = t227 * t285;
t273 = t224 * t285;
t122 = qJDD(1) * t173 - t218 * t273 + t220 * t272;
t73 = qJD(4) * t286 + t223 * qJDD(2) + t226 * t122 - t161 * t288;
t74 = qJD(4) * t134 - t226 * qJDD(2) + t223 * t122;
t39 = t217 * t73 + t219 * t74;
t40 = -t217 * t74 + t219 * t73;
t343 = t39 * pkin(5) - t40 * qJ(6) - t253 * qJD(6);
t213 = qJ(4) + pkin(10);
t206 = cos(t213);
t204 = sin(t213);
t310 = t204 * t225;
t138 = t206 * t228 + t207 * t310;
t296 = t228 * t204;
t299 = t225 * t206;
t140 = t207 * t296 - t299;
t159 = t173 * qJD(2);
t283 = t227 * qJDD(1);
t284 = t224 * qJDD(1);
t254 = -t218 * t284 + t220 * t283;
t116 = qJD(1) * t159 + qJDD(4) - t254;
t121 = -qJD(2) * t161 + t254;
t242 = pkin(2) * t273 - qJDD(1) * t203 + qJDD(3);
t62 = -t121 * pkin(3) - t122 * pkin(8) + t242;
t58 = t226 * t62;
t325 = qJD(2) * pkin(2);
t167 = -t175 + t325;
t120 = t218 * t167 + t304;
t107 = qJD(2) * pkin(8) + t120;
t179 = -qJD(1) * t203 + qJD(3);
t86 = -t158 * pkin(3) - t161 * pkin(8) + t179;
t60 = t107 * t226 + t223 * t86;
t70 = t218 * t115 + t220 * t123;
t68 = qJDD(2) * pkin(8) + t70;
t12 = t116 * pkin(4) - t73 * qJ(5) - qJD(4) * t60 - t134 * qJD(5) - t223 * t68 + t58;
t245 = -t107 * t288 + t223 * t62 + t226 * t68 + t86 * t287;
t15 = -qJ(5) * t74 - qJD(5) * t132 + t245;
t3 = t219 * t12 - t217 * t15;
t277 = -qJDD(6) + t3;
t333 = g(3) * t205;
t119 = t167 * t220 - t164;
t106 = -qJD(2) * pkin(3) - t119;
t75 = pkin(4) * t132 + qJD(5) + t106;
t36 = pkin(5) * t77 - qJ(6) * t253 + t75;
t342 = g(1) * t140 + g(2) * t138 + t204 * t333 - t36 * t253 + t277;
t341 = t262 * t207 + t333;
t340 = t149 ^ 2;
t338 = pkin(2) * t220;
t337 = pkin(5) * t116;
t332 = g(3) * t227;
t51 = -qJ(5) * t132 + t60;
t47 = t219 * t51;
t59 = -t107 * t223 + t226 * t86;
t50 = -qJ(5) * t134 + t59;
t23 = t217 * t50 + t47;
t331 = t23 * t253;
t4 = t217 * t12 + t219 * t15;
t118 = pkin(3) * t171 - pkin(8) * t173 - t203;
t130 = -t180 * t218 + t181 * t220;
t126 = t226 * t130;
t252 = qJ(5) * t163 - qJD(5) * t173;
t280 = t224 * t325;
t97 = pkin(3) * t159 + pkin(8) * t163 + t280;
t93 = t226 * t97;
t95 = t155 * t220 + t156 * t218;
t26 = t159 * pkin(4) - t223 * t95 + t93 + t252 * t226 + (-t126 + (qJ(5) * t173 - t118) * t223) * qJD(4);
t282 = t118 * t287 + t223 * t97 + t226 * t95;
t30 = -qJ(5) * t275 + (-qJD(4) * t130 + t252) * t223 + t282;
t9 = t217 * t26 + t219 * t30;
t44 = pkin(4) * t149 + t50;
t20 = t217 * t44 + t47;
t28 = t217 * t49 + t219 * t56;
t105 = t226 * t118;
t311 = t173 * t226;
t55 = pkin(4) * t171 - qJ(5) * t311 - t130 * t223 + t105;
t291 = t223 * t118 + t126;
t312 = t173 * t223;
t63 = -qJ(5) * t312 + t291;
t34 = t217 * t55 + t219 * t63;
t329 = t320 * pkin(5) - t319 * qJ(6) - qJD(6) * t172 + t345;
t327 = pkin(5) * t161 - t328;
t21 = qJ(6) * t161 + t28;
t83 = t219 * t137 + t217 * t241;
t326 = t83 - t21;
t324 = t217 * t51;
t322 = t73 * t223;
t318 = t132 * t149;
t317 = t132 * t161;
t316 = t134 * t149;
t315 = t134 * t161;
t309 = t205 * t225;
t308 = t205 * t228;
t202 = pkin(4) * t226 + pkin(3);
t177 = t207 * t202;
t307 = t207 * t228;
t302 = t223 * t116;
t301 = t223 * t225;
t300 = t223 * t228;
t298 = t225 * t226;
t103 = t226 * t116;
t297 = t226 * t228;
t295 = t228 * t330;
t24 = t219 * t50 - t324;
t292 = qJD(6) - t24;
t215 = t224 ^ 2;
t290 = -t227 ^ 2 + t215;
t281 = t116 * qJ(6) + t4;
t278 = t207 * t300;
t269 = -qJD(4) * t86 - t68;
t268 = t293 * t223;
t94 = t155 * t218 - t220 * t156;
t129 = t220 * t180 + t181 * t218;
t266 = t228 * t203 + t225 * t330;
t264 = t149 * t226;
t263 = g(1) * t309 - g(2) * t308;
t261 = -t107 * t287 + t58;
t260 = -t77 ^ 2 - t348;
t258 = pkin(4) * t312 + t129;
t221 = -qJ(5) - pkin(8);
t257 = -t205 * t221 + t177 + t211;
t256 = -pkin(2) * t224 - t207 * t221;
t255 = pkin(5) * t206 + qJ(6) * t204;
t8 = -t217 * t30 + t219 * t26;
t19 = t219 * t44 - t324;
t33 = -t217 * t63 + t219 * t55;
t251 = -t202 - t338;
t250 = t349 * pkin(4) + t94;
t249 = -t351 * t149 + t103;
t247 = -0.2e1 * pkin(1) * t285 - pkin(7) * qJDD(2);
t150 = t207 * t301 + t297;
t246 = -t163 * t226 - t173 * t288;
t244 = t149 * t106 - t198 * t116;
t239 = pkin(4) * t301 + t202 * t307 - t221 * t308 + t266;
t229 = qJD(2) ^ 2;
t238 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t229 + t271;
t230 = qJD(1) ^ 2;
t237 = pkin(1) * t230 - pkin(7) * qJDD(1) + t262;
t38 = t74 * pkin(4) + qJDD(5) + t67;
t236 = t295 + t221 * t309 + pkin(4) * t300 + (-t203 - t177) * t225;
t235 = -t172 * t39 + t253 * t320 - t319 * t77 - t344 * t40;
t168 = t293 * t226;
t110 = t217 * t168 + t219 * t268;
t111 = t219 * t168 - t217 * t268;
t232 = t110 * t40 - t111 * t39 - t83 * t77 - t341;
t231 = t38 + t346;
t200 = -pkin(3) - t338;
t199 = -pkin(4) * t219 - pkin(5);
t194 = pkin(4) * t217 + qJ(6);
t192 = pkin(4) * t298;
t153 = t207 * t297 + t301;
t152 = -t278 + t298;
t151 = -t207 * t298 + t300;
t141 = t206 * t307 + t310;
t139 = t207 * t299 - t296;
t100 = t172 * t173;
t98 = -pkin(5) * t344 - qJ(6) * t172 + t251;
t65 = t157 * t173 + t344 * t163;
t64 = -qJD(4) * t101 + t163 * t172;
t45 = pkin(5) * t100 - qJ(6) * t101 + t258;
t42 = pkin(4) * t134 + pkin(5) * t253 + qJ(6) * t77;
t32 = -pkin(5) * t171 - t33;
t31 = qJ(6) * t171 + t34;
t18 = qJ(6) * t149 + t20;
t17 = -pkin(5) * t149 + qJD(6) - t19;
t16 = -pkin(5) * t64 + qJ(6) * t65 - qJD(6) * t101 + t250;
t7 = -pkin(5) * t159 - t8;
t6 = qJ(6) * t159 + qJD(6) * t171 + t9;
t5 = t38 + t343;
t2 = -t277 - t337;
t1 = qJD(6) * t149 + t281;
t10 = [qJDD(1), t271, t262, qJDD(1) * t215 + 0.2e1 * t224 * t272, 0.2e1 * t224 * t283 - 0.2e1 * t290 * t285, qJDD(2) * t224 + t227 * t229, qJDD(2) * t227 - t224 * t229, 0, t224 * t247 + t227 * t238, -t224 * t238 + t227 * t247, t119 * t163 - t120 * t159 + t121 * t130 + t122 * t129 + t158 * t95 + t161 * t94 - t171 * t70 - t173 * t69 - t262, t70 * t130 + t120 * t95 - t69 * t129 - t119 * t94 - t242 * t203 + t179 * t280 - g(1) * (-t203 * t225 + t295) - g(2) * t266, t134 * t246 + t311 * t73 -(-t132 * t226 - t134 * t223) * t163 + (-t322 - t226 * t74 + (t132 * t223 - t134 * t226) * qJD(4)) * t173, t103 * t173 + t134 * t159 + t149 * t246 + t73 * t171, -t132 * t159 - t349 * t149 - t74 * t171 - t173 * t302, t116 * t171 + t149 * t159 (-t130 * t287 + t93) * t149 + t105 * t116 + t261 * t171 + t59 * t159 + t94 * t132 + t129 * t74 + t106 * t275 - g(1) * t151 - g(2) * t153 + ((-qJD(4) * t118 - t95) * t149 - t130 * t116 + t269 * t171 + t67 * t173 - t106 * t163) * t223 -(-t130 * t288 + t282) * t149 - t291 * t116 - t245 * t171 - t60 * t159 + t94 * t134 + t129 * t73 + t67 * t311 - g(1) * t150 - g(2) * t152 + t246 * t106, -t100 * t4 - t101 * t3 + t19 * t65 + t20 * t64 - t253 * t8 - t33 * t40 - t34 * t39 - t77 * t9 + t263, -g(1) * t236 - g(2) * t239 + t19 * t8 + t20 * t9 + t250 * t75 + t258 * t38 + t3 * t33 + t4 * t34, g(1) * t139 - g(2) * t141 + t100 * t5 - t116 * t32 - t149 * t7 - t159 * t17 + t16 * t77 - t171 * t2 - t36 * t64 + t39 * t45, -t1 * t100 + t101 * t2 - t17 * t65 + t18 * t64 + t253 * t7 - t31 * t39 + t32 * t40 - t6 * t77 + t263, g(1) * t138 - g(2) * t140 + t1 * t171 - t101 * t5 + t116 * t31 + t149 * t6 + t159 * t18 - t16 * t253 + t36 * t65 - t40 * t45, t1 * t31 + t18 * t6 + t5 * t45 + t36 * t16 + t2 * t32 + t17 * t7 - g(1) * (-t139 * pkin(5) - t138 * qJ(6) + t236) - g(2) * (pkin(5) * t141 + qJ(6) * t140 + t239); 0, 0, 0, -t224 * t230 * t227, t290 * t230, t284, t283, qJDD(2), t224 * t237 - t332, g(3) * t224 + t227 * t237 (t120 - t124) * t161 + (t119 - t125) * t158 + (t121 * t218 - t122 * t220) * pkin(2), t119 * t124 - t120 * t125 + (-t332 + t218 * t70 + t220 * t69 + (-qJD(1) * t179 + t262) * t224) * pkin(2), t134 * t264 + t322 (t73 - t318) * t226 + (-t74 - t316) * t223, t149 * t264 + t302 - t315, t249 + t317, -t149 * t161, -t124 * t132 - t92 * t149 - t59 * t161 + t200 * t74 + (t125 * t149 + t244) * t223 - t352 * t226, -t124 * t134 + t321 * t149 + t60 * t161 + t200 * t73 + t352 * t223 + t244 * t226, -t3 * t172 - t19 * t319 - t20 * t320 - t253 * t328 + t28 * t77 + t344 * t4 + t232, t4 * t111 - t3 * t110 + t38 * t251 - g(3) * t257 + t345 * t75 + (t83 - t28) * t20 + t328 * t19 + t262 * (t202 * t205 - t256) -t110 * t116 - t149 * t327 + t17 * t161 - t206 * t346 + t320 * t36 + t329 * t77 - t344 * t5 + t98 * t39, t1 * t344 + t17 * t319 + t2 * t172 - t18 * t320 + t21 * t77 + t253 * t327 + t232, t111 * t116 + t149 * t326 - t18 * t161 - t5 * t172 - t204 * t346 - t253 * t329 - t319 * t36 - t98 * t40, t1 * t111 + t5 * t98 + t2 * t110 - g(3) * (t207 * t255 + t257) + t329 * t36 + t326 * t18 + t327 * t17 + t262 * (-(-t202 - t255) * t205 - t256); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t158 ^ 2 - t161 ^ 2, t119 * t161 - t120 * t158 + t242 - t271, 0, 0, 0, 0, 0, t249 - t317, -t340 * t226 - t302 - t315, t235, -t75 * t161 + t4 * t172 - t19 * t320 + t20 * t319 + t3 * t344 - t271, t116 * t344 - t149 * t320 - t161 * t77, t235, t172 * t116 + t149 * t319 + t161 * t253, t1 * t172 - t36 * t161 + t17 * t320 + t18 * t319 - t2 * t344 - t271; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134 * t132, -t132 ^ 2 + t134 ^ 2, t73 + t318, t316 - t74, t116, -g(1) * t152 + g(2) * t150 - t106 * t134 + t60 * t149 + (t269 + t333) * t223 + t261, g(1) * t153 - g(2) * t151 + t106 * t132 + t149 * t59 + t226 * t333 - t245, t20 * t253 - t331 + (-t217 * t39 - t219 * t40) * pkin(4) + (-t19 + t24) * t77, -g(1) * t192 + t19 * t23 - t20 * t24 + (g(2) * t297 - t75 * t134 + t4 * t217 + t3 * t219 + t341 * t223) * pkin(4), t23 * t149 - t42 * t77 + (pkin(5) - t199) * t116 + t342, t18 * t253 - t194 * t39 + t199 * t40 - t331 + (t17 - t292) * t77, -t206 * t333 - g(1) * t141 - g(2) * t139 + t194 * t116 - t36 * t77 + t42 * t253 + (0.2e1 * qJD(6) - t24) * t149 + t281, t1 * t194 + t2 * t199 - t36 * t42 - t17 * t23 - g(1) * (-pkin(4) * t278 - pkin(5) * t140 + qJ(6) * t141 + t192) - g(2) * (-pkin(4) * t150 - t138 * pkin(5) + t139 * qJ(6)) + t292 * t18 - (-pkin(4) * t223 - pkin(5) * t204 + qJ(6) * t206) * t333; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t260, t19 * t253 + t20 * t77 + t231, t149 * t253 + t39, t260, -t40 + t350, -t17 * t253 + t18 * t77 + t231 + t343; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253 * t77 - qJDD(4) + t121, t40 + t350, -t340 - t348, -t149 * t18 - t337 - t342;];
tau_reg  = t10;
