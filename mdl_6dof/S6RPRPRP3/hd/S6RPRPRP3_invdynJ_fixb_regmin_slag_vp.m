% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:09:40
% EndTime: 2019-03-09 03:09:50
% DurationCPUTime: 5.04s
% Computational Cost: add. (5232->452), mult. (10914->581), div. (0->0), fcn. (7484->14), ass. (0->220)
t284 = qJDD(3) * pkin(3);
t189 = cos(qJ(3));
t176 = g(3) * t189;
t187 = sin(qJ(3));
t178 = qJ(1) + pkin(9);
t172 = cos(t178);
t170 = sin(t178);
t301 = g(2) * t170;
t235 = g(1) * t172 + t301;
t203 = t235 * t187 - t176;
t182 = sin(pkin(9));
t165 = pkin(1) * t182 + pkin(7);
t147 = t165 * qJD(1);
t262 = qJD(3) * t189;
t145 = t165 * qJDD(1);
t316 = -qJD(2) * qJD(3) - t145;
t223 = t189 * qJDD(2) - t147 * t262 + t316 * t187;
t58 = qJDD(4) - t223 - t284;
t318 = t58 - t203;
t322 = -t318 + t284;
t173 = t189 * qJDD(1);
t260 = qJD(1) * qJD(3);
t310 = t187 * t260 - t173;
t128 = qJDD(5) + t310;
t177 = pkin(10) + qJ(5);
t169 = sin(t177);
t181 = sin(pkin(10));
t296 = pkin(8) + qJ(4);
t143 = t296 * t181;
t183 = cos(pkin(10));
t144 = t296 * t183;
t186 = sin(qJ(5));
t307 = cos(qJ(5));
t85 = -t186 * t143 + t307 * t144;
t321 = -t85 * t128 - t203 * t169;
t264 = qJD(3) * t183;
t267 = qJD(1) * t187;
t125 = -t181 * t267 + t264;
t253 = t183 * t267;
t126 = qJD(3) * t181 + t253;
t75 = -t307 * t125 + t126 * t186;
t320 = t75 ^ 2;
t215 = -t186 * t125 - t307 * t126;
t308 = t215 ^ 2;
t266 = qJD(1) * t189;
t159 = -qJD(5) + t266;
t319 = t159 * t75;
t317 = t159 * t215;
t213 = -t186 * t181 + t307 * t183;
t248 = qJD(5) * t307;
t261 = qJD(5) * t186;
t312 = -t181 * t261 + t183 * t248;
t271 = -t213 * t266 + t312;
t130 = t307 * t181 + t186 * t183;
t117 = t130 * qJD(5);
t208 = t189 * t130;
t270 = -qJD(1) * t208 + t117;
t246 = t189 * t260;
t257 = t187 * qJDD(1);
t212 = t246 + t257;
t258 = t181 * qJDD(3);
t198 = t183 * t212 + t258;
t155 = t181 * t257;
t269 = t181 * t246 + t155;
t227 = qJDD(3) * t183 - t269;
t26 = -t125 * t248 + t126 * t261 - t186 * t227 - t307 * t198;
t263 = qJD(3) * t187;
t244 = t189 * t26 - t215 * t263;
t107 = t213 * t187;
t250 = t183 * t262;
t251 = t181 * t262;
t63 = t187 * t117 + t186 * t251 - t307 * t250;
t292 = t107 * t128 + t63 * t159;
t315 = t244 - t292;
t214 = -t307 * t143 - t186 * t144;
t274 = t183 * t189;
t226 = pkin(4) * t187 - pkin(8) * t274;
t110 = qJD(2) * t189 - t187 * t147;
t230 = pkin(3) * t187 - qJ(4) * t189;
t135 = t230 * qJD(1);
t68 = -t181 * t110 + t183 * t135;
t39 = qJD(1) * t226 + t68;
t254 = t181 * t266;
t69 = t183 * t110 + t181 * t135;
t50 = -pkin(8) * t254 + t69;
t314 = qJD(4) * t213 + qJD(5) * t214 - t186 * t39 - t307 * t50;
t313 = qJD(4) * t130 + qJD(5) * t85 - t186 * t50 + t307 * t39;
t27 = -qJD(5) * t215 + t186 * t198 - t307 * t227;
t243 = -t189 * t27 + t75 * t263;
t106 = t130 * t187;
t64 = qJD(3) * t208 + t312 * t187;
t291 = -t106 * t128 + t64 * t159;
t311 = t243 + t291;
t231 = pkin(3) * t189 + qJ(4) * t187;
t225 = -pkin(2) - t231;
t184 = cos(pkin(9));
t306 = pkin(1) * t184;
t120 = t225 - t306;
t109 = t183 * t120;
t275 = t183 * t187;
t60 = -pkin(8) * t275 + t109 + (-t165 * t181 - pkin(4)) * t189;
t277 = t181 * t187;
t83 = t181 * t120 + t165 * t274;
t67 = -pkin(8) * t277 + t83;
t219 = t186 * t60 + t307 * t67;
t114 = qJD(3) * t230 - t187 * qJD(4);
t252 = t165 * t263;
t72 = t183 * t114 + t181 * t252;
t48 = qJD(3) * t226 + t72;
t103 = t181 * t114;
t149 = t187 * t165;
t276 = t181 * t189;
t59 = t103 + (-pkin(8) * t276 - t183 * t149) * qJD(3);
t309 = -qJD(5) * t219 - t186 * t59 + t307 * t48;
t180 = t189 ^ 2;
t305 = pkin(4) * t181;
t304 = pkin(5) * t128;
t303 = g(1) * t170;
t300 = g(2) * t172;
t299 = g(3) * t187;
t102 = t120 * qJD(1);
t111 = t187 * qJD(2) + t189 * t147;
t99 = qJD(3) * qJ(4) + t111;
t42 = t183 * t102 - t181 * t99;
t33 = -pkin(4) * t266 - pkin(8) * t126 + t42;
t43 = t181 * t102 + t183 * t99;
t36 = pkin(8) * t125 + t43;
t9 = t186 * t33 + t307 * t36;
t298 = t159 * t9;
t297 = t215 * t75;
t295 = -t107 * t27 + t63 * t75;
t86 = pkin(4) * t254 + t111;
t293 = t270 * pkin(5) - t271 * qJ(6) - qJD(6) * t130 - t86;
t52 = qJDD(3) * qJ(4) + t187 * qJDD(2) + t189 * t145 + (qJD(4) + t110) * qJD(3);
t61 = qJD(1) * t114 + qJDD(1) * t120;
t23 = t181 * t61 + t183 * t52;
t290 = -qJ(6) * t267 + t314;
t289 = pkin(5) * t267 + t313;
t8 = -t186 * t36 + t307 * t33;
t285 = qJD(6) - t8;
t283 = t128 * qJ(6);
t282 = t169 * t170;
t171 = cos(t177);
t281 = t171 * t187;
t280 = t171 * t189;
t279 = t172 * t187;
t278 = t172 * t189;
t273 = t187 * t296;
t272 = qJDD(2) - g(3);
t105 = pkin(4) * t251 + t165 * t262;
t113 = pkin(4) * t277 + t149;
t179 = t187 ^ 2;
t268 = t179 - t180;
t167 = -pkin(2) - t306;
t148 = qJD(1) * t167;
t190 = cos(qJ(1));
t256 = t190 * pkin(1) + t172 * pkin(2) + t170 * pkin(7);
t166 = pkin(4) * t183 + pkin(3);
t255 = qJ(4) * t263;
t188 = sin(qJ(1));
t245 = -t188 * pkin(1) + t172 * pkin(7);
t22 = -t181 * t52 + t183 * t61;
t12 = t310 * pkin(4) - t198 * pkin(8) + t22;
t17 = pkin(8) * t227 + t23;
t242 = -t307 * t12 + t186 * t17 + t36 * t248 + t33 * t261;
t241 = -qJD(3) * pkin(3) + qJD(4);
t240 = qJD(1) * t268;
t93 = t171 * t172 + t189 * t282;
t95 = t169 * t278 - t170 * t171;
t238 = -g(1) * t93 + g(2) * t95;
t94 = -t172 * t169 + t170 * t280;
t96 = t171 * t278 + t282;
t237 = g(1) * t94 - g(2) * t96;
t152 = t187 * t303;
t236 = -g(2) * t279 + t152;
t234 = -t300 + t303;
t233 = g(1) * t188 - g(2) * t190;
t232 = -t106 * t26 - t215 * t64;
t229 = -t22 * t181 + t23 * t183;
t228 = -t181 * t42 + t183 * t43;
t224 = pkin(5) * t171 + qJ(6) * t169 + t166;
t221 = -t186 * t67 + t307 * t60;
t217 = t186 * t12 + t307 * t17 + t33 * t248 - t261 * t36;
t216 = t186 * t48 + t60 * t248 - t261 * t67 + t307 * t59;
t211 = t183 * t257 + t258;
t210 = -qJD(1) * t148 + t235;
t97 = -t110 + t241;
t209 = g(1) * t171 * t279 - g(3) * t280 + t128 * t214 + t281 * t301;
t191 = qJD(3) ^ 2;
t207 = 0.2e1 * qJDD(1) * t167 + t165 * t191 + t300;
t206 = g(1) * t95 + g(2) * t93 + t169 * t299 - t242;
t205 = 0.2e1 * qJD(3) * t148 - qJDD(3) * t165;
t204 = t26 - t319;
t201 = -t189 * t235 - t299;
t70 = -pkin(4) * t125 + t97;
t20 = pkin(5) * t75 + qJ(6) * t215 + t70;
t197 = -t20 * t215 + qJDD(6) - t206;
t196 = -g(1) * t96 - g(2) * t94 - g(3) * t281 + t217;
t35 = -pkin(4) * t227 + t58;
t3 = t27 * pkin(5) + t26 * qJ(6) + qJD(6) * t215 + t35;
t194 = t27 + t317;
t192 = qJD(1) ^ 2;
t140 = qJDD(3) * t189 - t187 * t191;
t139 = qJDD(3) * t187 + t189 * t191;
t82 = -t165 * t276 + t109;
t73 = -t183 * t252 + t103;
t71 = -pkin(5) * t213 - qJ(6) * t130 - t166;
t37 = pkin(5) * t106 - qJ(6) * t107 + t113;
t32 = -pkin(5) * t215 + qJ(6) * t75;
t25 = t189 * pkin(5) - t221;
t24 = -qJ(6) * t189 + t219;
t14 = pkin(5) * t64 + qJ(6) * t63 - qJD(6) * t107 + t105;
t13 = -t26 - t319;
t7 = -t159 * qJ(6) + t9;
t6 = t159 * pkin(5) + t285;
t5 = -pkin(5) * t263 - t309;
t4 = qJ(6) * t263 - qJD(6) * t189 + t216;
t2 = qJDD(6) + t242 - t304;
t1 = -qJD(6) * t159 + t217 + t283;
t10 = [qJDD(1), t233, g(1) * t190 + g(2) * t188 (t233 + (t182 ^ 2 + t184 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t179 + 0.2e1 * t187 * t246, -0.2e1 * qJD(3) * t240 + 0.2e1 * t173 * t187, t139, t140, 0, t205 * t187 + (-t207 + t303) * t189, t187 * t207 + t189 * t205 - t152, -t235 * t181 + (-t165 * t227 + t58 * t181 + (qJD(1) * t82 + t42) * qJD(3)) * t187 + (-t72 * qJD(1) - t82 * qJDD(1) - t22 + t234 * t183 + (-t125 * t165 + t181 * t97) * qJD(3)) * t189, -t235 * t183 + (t58 * t183 + t211 * t165 + (-qJD(1) * t83 - t43) * qJD(3)) * t187 + (t73 * qJD(1) + t83 * qJDD(1) + t23 - t234 * t181 + (t183 * t97 + (t126 + t253) * t165) * qJD(3)) * t189, t73 * t125 - t83 * t269 - t72 * t126 + (-qJDD(3) * t82 - t187 * t23 - t262 * t43) * t181 + (t83 * qJDD(3) - t22 * t187 - t212 * t82 - t262 * t42) * t183 + t236, t23 * t83 + t43 * t73 + t22 * t82 + t42 * t72 - g(1) * t245 - g(2) * (t172 * t231 + t256) - t225 * t303 + (t58 * t187 + t262 * t97) * t165, -t107 * t26 + t215 * t63, -t232 + t295, t244 + t292, t291 - t243, -t128 * t189 - t159 * t263, t105 * t75 + t35 * t106 + t113 * t27 + t221 * t128 - t309 * t159 + t242 * t189 + t8 * t263 + t70 * t64 + t237, -t105 * t215 + t35 * t107 - t113 * t26 - t219 * t128 + t216 * t159 + t217 * t189 - t9 * t263 - t70 * t63 + t238, t106 * t3 - t128 * t25 + t14 * t75 + t159 * t5 + t189 * t2 + t20 * t64 - t263 * t6 + t27 * t37 + t237, -t1 * t106 + t107 * t2 - t215 * t5 - t24 * t27 - t25 * t26 - t4 * t75 - t6 * t63 - t64 * t7 + t236, -t1 * t189 - t107 * t3 + t128 * t24 + t14 * t215 - t159 * t4 + t20 * t63 + t26 * t37 + t263 * t7 - t238, t1 * t24 + t7 * t4 + t3 * t37 + t20 * t14 + t2 * t25 + t6 * t5 - g(1) * (-t94 * pkin(5) - t93 * qJ(6) + t172 * t305 + t245) - g(2) * (t96 * pkin(5) + t95 * qJ(6) + t166 * t278 + t172 * t273 + t256) + (-g(1) * (-t166 * t189 - pkin(2) - t273) - g(2) * t305) * t170; 0, 0, 0, t272, 0, 0, 0, 0, 0, t140, -t139 (t227 + t155) * t189 + (-t187 * t125 - t181 * t240) * qJD(3), -t189 * t258 + (t187 * t126 + (-t268 - t180) * t183 * qJD(1)) * qJD(3), t125 * t250 + t126 * t251 + t198 * t277 + t227 * t275, -t58 * t189 - g(3) + t229 * t187 + (t187 * t97 + t189 * t228) * qJD(3), 0, 0, 0, 0, 0, t311, t315, t311, t232 + t295, -t315, t1 * t107 + t106 * t2 - t189 * t3 + t20 * t263 + t6 * t64 - t63 * t7 - g(3); 0, 0, 0, 0, -t187 * t192 * t189, t268 * t192, t257, t173, qJDD(3), t111 * qJD(3) + t187 * t210 - t176 + t223, t110 * qJD(3) + (qJD(3) * t147 - t272) * t187 + (t210 + t316) * t189, t181 * qJ(4) * t173 - pkin(3) * t269 + t111 * t125 + t322 * t183 + (-t42 * t187 + t68 * t189 + (-t255 + (qJD(4) - t97) * t189) * t181) * qJD(1), -t111 * t126 - t230 * t183 * qJDD(1) - t322 * t181 + (t43 * t187 - t69 * t189 + (-t255 + (t241 - t97) * t189) * t183) * qJD(1), -t69 * t125 + t68 * t126 + (qJ(4) * t227 + qJD(4) * t125 + t266 * t42 + t23) * t183 + (qJ(4) * t198 + qJD(4) * t126 + t266 * t43 - t22) * t181 + t201, -t97 * t111 - t42 * t68 - t43 * t69 + t228 * qJD(4) - t318 * pkin(3) + (t201 + t229) * qJ(4), -t26 * t130 - t215 * t271, -t130 * t27 - t213 * t26 + t215 * t270 - t271 * t75, t130 * t128 - t271 * t159 + t215 * t267, t128 * t213 + t270 * t159 + t75 * t267, t159 * t267, t313 * t159 - t166 * t27 - t213 * t35 - t8 * t267 + t270 * t70 - t86 * t75 + t209, t35 * t130 + t314 * t159 + t166 * t26 + t215 * t86 + t9 * t267 + t271 * t70 + t321, t289 * t159 + t270 * t20 - t213 * t3 + t6 * t267 + t71 * t27 + t293 * t75 + t209, t1 * t213 + t2 * t130 + t214 * t26 - t215 * t289 - t85 * t27 - t270 * t7 + t271 * t6 - t290 * t75 + t201, -t3 * t130 - t290 * t159 - t271 * t20 + t215 * t293 + t71 * t26 - t7 * t267 - t321, -g(3) * t273 + t1 * t85 - t224 * t176 - t2 * t214 + t293 * t20 + t289 * t6 + t290 * t7 + t3 * t71 + t235 * (t187 * t224 - t189 * t296); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126 * t266 - t227 (-t125 + t264) * t266 + t211, -t125 ^ 2 - t126 ^ 2, -t43 * t125 + t42 * t126 + t318, 0, 0, 0, 0, 0, t194, -t204, t194, -t308 - t320, t204, t215 * t6 + t7 * t75 - t203 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t297, t308 - t320, t13, -t27 + t317, t128, t215 * t70 + t206 - t298, -t159 * t8 + t70 * t75 - t196, -t32 * t75 - t197 - t298 + 0.2e1 * t304, pkin(5) * t26 - t27 * qJ(6) - (t7 - t9) * t215 + (t6 - t285) * t75, 0.2e1 * t283 - t20 * t75 - t32 * t215 + (-0.2e1 * qJD(6) + t8) * t159 + t196, t1 * qJ(6) - t2 * pkin(5) - t20 * t32 - t6 * t9 - g(1) * (-pkin(5) * t95 + qJ(6) * t96) - g(2) * (-pkin(5) * t93 + qJ(6) * t94) + t285 * t7 - (-pkin(5) * t169 + qJ(6) * t171) * t299; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128 - t297, t13, -t159 ^ 2 - t308, t159 * t7 + t197 - t304;];
tau_reg  = t10;
