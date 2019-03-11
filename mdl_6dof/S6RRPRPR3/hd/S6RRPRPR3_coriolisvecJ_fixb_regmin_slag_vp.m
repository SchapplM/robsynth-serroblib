% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPRPR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:19:07
% EndTime: 2019-03-09 10:19:21
% DurationCPUTime: 4.92s
% Computational Cost: add. (8041->399), mult. (20653->557), div. (0->0), fcn. (15929->10), ass. (0->220)
t220 = sin(pkin(10));
t225 = sin(qJ(2));
t274 = qJD(1) * t225;
t222 = cos(pkin(10));
t228 = cos(qJ(2));
t286 = t222 * t228;
t186 = qJD(1) * t286 - t220 * t274;
t180 = qJD(4) - t186;
t223 = sin(qJ(6));
t226 = cos(qJ(6));
t270 = qJD(6) * t223;
t201 = t220 * t228 + t222 * t225;
t189 = t201 * qJD(1);
t224 = sin(qJ(4));
t227 = cos(qJ(4));
t269 = t227 * qJD(2);
t156 = t189 * t224 - t269;
t158 = qJD(2) * t224 + t189 * t227;
t219 = sin(pkin(11));
t221 = cos(pkin(11));
t244 = -t156 * t221 - t158 * t219;
t305 = t226 * t244;
t268 = qJD(1) * qJD(2);
t259 = t228 * t268;
t260 = t225 * t268;
t174 = -t220 * t260 + t222 * t259;
t272 = qJD(4) * t224;
t112 = qJD(4) * t269 + t227 * t174 - t189 * t272;
t113 = qJD(4) * t158 + t224 * t174;
t66 = -t112 * t219 - t113 * t221;
t67 = t112 * t221 - t113 * t219;
t97 = t156 * t219 - t158 * t221;
t12 = qJD(6) * t305 + t223 * t66 + t226 * t67 + t270 * t97;
t172 = qJD(6) + t180;
t54 = t223 * t97 + t305;
t309 = t172 * t54;
t334 = t12 - t309;
t320 = t223 * t244 - t226 * t97;
t333 = t320 * t54;
t200 = t219 * t227 + t221 * t224;
t323 = t180 * t200;
t242 = t219 * t224 - t221 * t227;
t322 = t180 * t242;
t332 = t320 ^ 2 - t54 ^ 2;
t264 = -pkin(2) * t228 - pkin(1);
t247 = t264 * qJD(1);
t205 = qJD(3) + t247;
t117 = -t186 * pkin(3) - t189 * pkin(8) + t205;
t313 = -qJ(3) - pkin(7);
t206 = t313 * t225;
t203 = qJD(1) * t206;
t310 = qJD(2) * pkin(2);
t195 = t203 + t310;
t207 = t313 * t228;
t204 = qJD(1) * t207;
t287 = t222 * t204;
t143 = t220 * t195 - t287;
t136 = qJD(2) * pkin(8) + t143;
t79 = t117 * t224 + t136 * t227;
t62 = -qJ(5) * t156 + t79;
t307 = t221 * t62;
t78 = t227 * t117 - t136 * t224;
t61 = -qJ(5) * t158 + t78;
t47 = pkin(4) * t180 + t61;
t29 = t219 * t47 + t307;
t318 = pkin(9) * t244;
t16 = t29 + t318;
t15 = t16 * t270;
t187 = t201 * qJD(2);
t173 = qJD(1) * t187;
t210 = pkin(2) * t260;
t116 = pkin(3) * t173 - pkin(8) * t174 + t210;
t108 = t227 * t116;
t256 = qJD(2) * t313;
t182 = t228 * qJD(3) + t225 * t256;
t166 = t182 * qJD(1);
t183 = -t225 * qJD(3) + t228 * t256;
t167 = t183 * qJD(1);
t111 = t166 * t222 + t167 * t220;
t231 = -qJD(4) * t79 - t224 * t111 + t108;
t21 = t173 * pkin(4) - t112 * qJ(5) - t158 * qJD(5) + t231;
t271 = qJD(4) * t227;
t234 = t227 * t111 + t224 * t116 + t117 * t271 - t136 * t272;
t24 = -qJ(5) * t113 - qJD(5) * t156 + t234;
t6 = t221 * t21 - t219 * t24;
t2 = pkin(5) * t173 - pkin(9) * t67 + t6;
t192 = t220 * t204;
t142 = t195 * t222 + t192;
t135 = -qJD(2) * pkin(3) - t142;
t93 = pkin(4) * t156 + qJD(5) + t135;
t48 = -pkin(5) * t244 + t93;
t331 = -t223 * t2 - t48 * t54 + t15;
t13 = qJD(6) * t320 + t223 * t67 - t226 * t66;
t304 = t320 * t172;
t329 = -t13 + t304;
t212 = pkin(2) * t220 + pkin(8);
t280 = qJ(5) + t212;
t254 = qJD(4) * t280;
t128 = pkin(2) * t274 + pkin(3) * t189 - pkin(8) * t186;
t148 = t203 * t222 + t192;
t277 = t224 * t128 + t227 * t148;
t293 = t186 * t224;
t328 = qJ(5) * t293 + t227 * qJD(5) - t224 * t254 - t277;
t123 = t227 * t128;
t327 = -pkin(4) * t189 - t123 + (qJ(5) * t186 - t254) * t227 + (-qJD(5) + t148) * t224;
t326 = t272 - t293;
t7 = t219 * t21 + t221 * t24;
t3 = pkin(9) * t66 + t7;
t263 = t226 * t2 - t223 * t3;
t325 = -t48 * t320 + t263;
t324 = pkin(9) * t97;
t145 = t200 * t226 - t223 * t242;
t311 = qJD(6) * t145 - t223 * t322 + t226 * t323;
t199 = t220 * t225 - t286;
t191 = t199 * qJD(2);
t262 = t201 * t271;
t321 = -t191 * t224 + t262;
t319 = -0.2e1 * t268;
t302 = -t328 * t219 + t327 * t221;
t301 = t327 * t219 + t328 * t221;
t147 = t203 * t220 - t287;
t316 = t326 * pkin(4) - t147;
t243 = -t200 * t223 - t226 * t242;
t312 = qJD(6) * t243 - t223 * t323 - t226 * t322;
t315 = -t145 * t173 - t312 * t172;
t314 = pkin(4) * t219;
t267 = t225 * t310;
t129 = pkin(3) * t187 + pkin(8) * t191 + t267;
t124 = t227 * t129;
t127 = t182 * t222 + t183 * t220;
t141 = pkin(3) * t199 - pkin(8) * t201 + t264;
t151 = t206 * t220 - t207 * t222;
t149 = t227 * t151;
t241 = qJ(5) * t191 - qJD(5) * t201;
t36 = t187 * pkin(4) - t224 * t127 + t124 + t241 * t227 + (-t149 + (qJ(5) * t201 - t141) * t224) * qJD(4);
t265 = t227 * t127 + t224 * t129 + t141 * t271;
t40 = -qJ(5) * t262 + (-qJD(4) * t151 + t241) * t224 + t265;
t11 = t219 * t36 + t221 * t40;
t56 = t219 * t62;
t32 = t221 * t61 - t56;
t134 = t227 * t141;
t290 = t201 * t227;
t71 = pkin(4) * t199 - qJ(5) * t290 - t151 * t224 + t134;
t276 = t224 * t141 + t149;
t291 = t201 * t224;
t82 = -qJ(5) * t291 + t276;
t43 = t219 * t71 + t221 * t82;
t308 = t189 * t54;
t28 = t221 * t47 - t56;
t14 = pkin(5) * t180 + t28 + t324;
t306 = t226 * t14;
t303 = t320 * t189;
t300 = pkin(5) * t323 + t316;
t299 = t112 * t224;
t297 = t156 * t180;
t296 = t156 * t189;
t295 = t158 * t180;
t294 = t158 * t189;
t285 = t224 * t173;
t164 = t227 * t173;
t230 = qJD(1) ^ 2;
t283 = t228 * t230;
t229 = qJD(2) ^ 2;
t282 = t229 * t225;
t281 = t229 * t228;
t196 = t280 * t224;
t197 = t280 * t227;
t138 = -t219 * t196 + t221 * t197;
t275 = t225 ^ 2 - t228 ^ 2;
t273 = qJD(4) * t201;
t214 = -pkin(2) * t222 - pkin(3);
t258 = qJD(6) * t14 + t3;
t10 = -t219 * t40 + t221 * t36;
t31 = -t219 * t61 - t307;
t42 = -t219 * t82 + t221 * t71;
t255 = pkin(1) * t319;
t110 = t166 * t220 - t222 * t167;
t126 = t182 * t220 - t222 * t183;
t137 = -t221 * t196 - t197 * t219;
t150 = -t222 * t206 - t207 * t220;
t253 = t180 * t227;
t252 = -t172 * t311 + t243 * t173;
t105 = -pkin(9) * t200 + t137;
t251 = pkin(9) * t323 - qJD(6) * t105 - t301;
t106 = -pkin(9) * t242 + t138;
t250 = pkin(5) * t189 - pkin(9) * t322 + qJD(6) * t106 - t302;
t248 = pkin(4) * t291 + t150;
t5 = t223 * t14 + t226 * t16;
t246 = t110 * t201 - t151 * t173;
t131 = t200 * t201;
t132 = t242 * t201;
t245 = -t226 * t131 + t132 * t223;
t86 = -t131 * t223 - t132 * t226;
t240 = -pkin(4) * t227 + t214;
t239 = pkin(4) * t321 + t126;
t238 = -t326 * t180 + t164;
t213 = pkin(4) * t221 + pkin(5);
t237 = t213 * t223 + t226 * t314;
t236 = t213 * t226 - t223 * t314;
t73 = pkin(4) * t113 + t110;
t235 = -t191 * t227 - t201 * t272;
t233 = t135 * t180 - t212 * t173;
t155 = pkin(5) * t242 + t240;
t146 = t173 * t199;
t88 = pkin(5) * t131 + t248;
t87 = pkin(4) * t158 - pkin(5) * t97;
t84 = -t242 * t191 + t200 * t273;
t83 = t191 * t200 + t242 * t273;
t44 = -pkin(5) * t83 + t239;
t41 = -pkin(5) * t66 + t73;
t35 = -pkin(9) * t131 + t43;
t30 = pkin(5) * t199 + pkin(9) * t132 + t42;
t27 = qJD(6) * t86 - t223 * t84 - t226 * t83;
t26 = qJD(6) * t245 + t223 * t83 - t226 * t84;
t18 = t32 + t324;
t17 = t31 - t318;
t9 = pkin(9) * t83 + t11;
t8 = pkin(5) * t187 + pkin(9) * t84 + t10;
t4 = -t16 * t223 + t306;
t1 = [0, 0, 0, 0.2e1 * t225 * t259, t275 * t319, t281, -t282, 0, -pkin(7) * t281 + t225 * t255, pkin(7) * t282 + t228 * t255, -t111 * t199 + t126 * t189 + t127 * t186 + t142 * t191 - t143 * t187 + t150 * t174 + t246, t110 * t150 + t111 * t151 - t142 * t126 + t143 * t127 + (t205 + t247) * t267, t112 * t290 + t158 * t235 -(-t156 * t227 - t158 * t224) * t191 + (-t299 - t113 * t227 + (t156 * t224 - t158 * t227) * qJD(4)) * t201, t112 * t199 + t158 * t187 + t201 * t164 + t180 * t235, -t113 * t199 - t156 * t187 - t180 * t321 - t201 * t285, t180 * t187 + t146 (-t151 * t271 + t124) * t180 + t134 * t173 + (-t136 * t271 + t108) * t199 + t78 * t187 + t126 * t156 + t150 * t113 + t135 * t262 + ((-qJD(4) * t141 - t127) * t180 + (-qJD(4) * t117 - t111) * t199 - t135 * t191 + t246) * t224 -(-t151 * t272 + t265) * t180 - t276 * t173 - t234 * t199 - t79 * t187 + t126 * t158 + t150 * t112 + t110 * t290 + t235 * t135, t10 * t97 + t11 * t244 - t131 * t7 + t132 * t6 + t28 * t84 + t29 * t83 - t42 * t67 + t43 * t66, t28 * t10 + t29 * t11 + t239 * t93 + t248 * t73 + t6 * t42 + t7 * t43, t12 * t86 + t26 * t320, t12 * t245 - t13 * t86 + t26 * t54 - t27 * t320, t12 * t199 + t172 * t26 + t173 * t86 + t187 * t320, -t13 * t199 - t172 * t27 + t173 * t245 + t187 * t54, t172 * t187 + t146 (-t223 * t9 + t226 * t8) * t172 + (-t223 * t35 + t226 * t30) * t173 + t263 * t199 + t4 * t187 - t44 * t54 + t88 * t13 - t41 * t245 + t48 * t27 + ((-t223 * t30 - t226 * t35) * t172 - t5 * t199) * qJD(6), t88 * t12 + t15 * t199 - t5 * t187 + t48 * t26 + t41 * t86 + t44 * t320 + (-(-qJD(6) * t35 + t8) * t172 - t30 * t173 - t2 * t199) * t223 + (-(qJD(6) * t30 + t9) * t172 - t35 * t173 - t258 * t199) * t226; 0, 0, 0, -t225 * t283, t275 * t230, 0, 0, 0, t230 * pkin(1) * t225, pkin(1) * t283 (t143 - t147) * t189 + (t142 - t148) * t186 + (-t173 * t220 - t174 * t222) * pkin(2), t142 * t147 - t143 * t148 + (-t110 * t222 + t111 * t220 - t205 * t274) * pkin(2), t158 * t253 + t299 (t112 - t297) * t227 + (-t113 - t295) * t224, t180 * t253 + t285 - t294, t238 + t296, -t180 * t189, -t110 * t227 + t214 * t113 - t147 * t156 - t78 * t189 + (-t212 * t271 - t123) * t180 + (t148 * t180 + t233) * t224, t110 * t224 + t214 * t112 - t147 * t158 + t79 * t189 + (t212 * t272 + t277) * t180 + t233 * t227, -t137 * t67 + t138 * t66 - t6 * t200 - t242 * t7 + t244 * t301 + t28 * t322 - t29 * t323 + t302 * t97, t6 * t137 + t7 * t138 + t73 * t240 + t302 * t28 + t301 * t29 + t316 * t93, t12 * t145 + t312 * t320, t12 * t243 - t145 * t13 - t311 * t320 + t312 * t54, -t303 - t315, t252 - t308, -t172 * t189 (t105 * t226 - t106 * t223) * t173 + t155 * t13 - t41 * t243 - t4 * t189 - t300 * t54 + t311 * t48 + (t223 * t251 - t226 * t250) * t172 -(t105 * t223 + t106 * t226) * t173 + t155 * t12 + t41 * t145 + t5 * t189 + t300 * t320 + t312 * t48 + (t223 * t250 + t226 * t251) * t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t186 ^ 2 - t189 ^ 2, t142 * t189 - t143 * t186 + t210, 0, 0, 0, 0, 0, t238 - t296, -t180 ^ 2 * t227 - t285 - t294, t200 * t66 + t242 * t67 - t244 * t322 - t323 * t97, -t93 * t189 + t7 * t200 - t242 * t6 - t28 * t323 - t29 * t322, 0, 0, 0, 0, 0, t252 + t308, -t303 + t315; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158 * t156, -t156 ^ 2 + t158 ^ 2, t112 + t297, -t113 + t295, t173, -t135 * t158 + t79 * t180 + t231, t135 * t156 + t180 * t78 - t234 (t219 * t66 - t221 * t67) * pkin(4) + (t28 - t32) * t244 + (-t31 - t29) * t97, -t28 * t31 - t29 * t32 + (-t158 * t93 + t219 * t7 + t221 * t6) * pkin(4), -t333, t332, t334, t329, t173, t236 * t173 - (t17 * t226 - t18 * t223) * t172 + t87 * t54 + (-t172 * t237 - t5) * qJD(6) + t325, -t237 * t173 - t226 * t3 + (t17 * t223 + t18 * t226) * t172 - t87 * t320 + (-t172 * t236 - t306) * qJD(6) + t331; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t244 ^ 2 - t97 ^ 2, -t244 * t29 - t28 * t97 + t73, 0, 0, 0, 0, 0, t13 + t304, t12 + t309; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t333, t332, t334, t329, t173 (-qJD(6) + t172) * t5 + t325, t4 * t172 - t226 * t258 + t331;];
tauc_reg  = t1;
