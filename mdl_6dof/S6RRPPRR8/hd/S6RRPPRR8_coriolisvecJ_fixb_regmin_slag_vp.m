% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
% 
% Output:
% tauc_reg [6x32]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRR8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:26:08
% EndTime: 2019-03-09 09:26:22
% DurationCPUTime: 4.89s
% Computational Cost: add. (4143->454), mult. (10374->627), div. (0->0), fcn. (7535->8), ass. (0->224)
t212 = cos(qJ(2));
t279 = qJD(1) * t212;
t191 = qJD(5) + t279;
t180 = qJD(6) + t191;
t207 = sin(qJ(6));
t210 = cos(qJ(6));
t204 = sin(pkin(10));
t209 = sin(qJ(2));
t280 = qJD(1) * t209;
t257 = t204 * t280;
t205 = cos(pkin(10));
t266 = t205 * qJD(2);
t149 = t257 - t266;
t255 = t205 * t280;
t278 = qJD(2) * t204;
t151 = t255 + t278;
t208 = sin(qJ(5));
t211 = cos(qJ(5));
t229 = t149 * t208 + t151 * t211;
t230 = -t211 * t149 + t151 * t208;
t267 = qJD(6) * t210;
t268 = qJD(6) * t207;
t264 = qJD(1) * qJD(2);
t252 = t212 * t264;
t176 = t204 * t252;
t240 = t205 * t252;
t269 = qJD(5) * t211;
t270 = qJD(5) * t208;
t45 = t149 * t269 - t151 * t270 + t208 * t176 + t211 * t240;
t46 = qJD(5) * t229 - t211 * t176 + t208 * t240;
t224 = t207 * t46 - t210 * t45 + t229 * t268 + t230 * t267;
t26 = t207 * t229 + t210 * t230;
t318 = t180 * t26;
t326 = -t224 + t318;
t232 = t207 * t230 - t210 * t229;
t325 = t232 * t26;
t324 = t180 * t232;
t156 = t204 * t208 + t205 * t211;
t222 = t156 * t212;
t286 = -qJD(1) * t222 - t156 * qJD(5);
t254 = t205 * t279;
t256 = t204 * t279;
t285 = t204 * t269 - t205 * t270 - t208 * t254 + t211 * t256;
t323 = t232 ^ 2 - t26 ^ 2;
t299 = qJ(3) * t209;
t167 = -pkin(2) * t212 - pkin(1) - t299;
t141 = t167 * qJD(1);
t195 = pkin(7) * t279;
t174 = qJD(2) * qJ(3) + t195;
t98 = t141 * t205 - t204 * t174;
t81 = pkin(3) * t279 + qJD(4) - t98;
t44 = pkin(4) * t279 - pkin(8) * t151 + t81;
t99 = t204 * t141 + t205 * t174;
t84 = -qJ(4) * t279 + t99;
t50 = pkin(8) * t149 + t84;
t17 = t208 * t44 + t211 * t50;
t13 = -pkin(9) * t230 + t17;
t11 = t13 * t268;
t194 = pkin(7) * t280;
t206 = qJD(2) * pkin(2);
t166 = qJD(3) + t194 - t206;
t74 = t149 * pkin(3) - t151 * qJ(4) + t166;
t47 = -pkin(4) * t149 - t74;
t23 = pkin(5) * t230 + t47;
t322 = t23 * t26 + t11;
t290 = t205 * t212;
t263 = pkin(8) * t290;
t309 = -pkin(3) - pkin(4);
t234 = pkin(2) * t209 - qJ(3) * t212;
t132 = qJD(2) * t234 - t209 * qJD(3);
t122 = t132 * qJD(1);
t164 = (qJD(3) - t194) * qJD(2);
t78 = t205 * t122 - t204 * t164;
t36 = (t309 * t209 - t263) * t264 - t78;
t251 = t209 * t264;
t79 = t204 * t122 + t205 * t164;
t261 = qJ(4) * t251 + t79;
t37 = (pkin(8) * t278 - qJD(4)) * t279 + t261;
t248 = -t208 * t37 + t211 * t36;
t217 = -qJD(5) * t17 + t248;
t4 = -pkin(5) * t251 - t45 * pkin(9) + t217;
t226 = t208 * t36 + t211 * t37 + t44 * t269 - t50 * t270;
t5 = -pkin(9) * t46 + t226;
t258 = -t207 * t5 + t210 * t4;
t16 = -t208 * t50 + t211 * t44;
t12 = -pkin(9) * t229 + t16;
t10 = pkin(5) * t191 + t12;
t303 = t13 * t210;
t3 = t10 * t207 + t303;
t321 = -qJD(6) * t3 + t23 * t232 + t258;
t216 = qJD(6) * t232 - t207 * t45 - t210 * t46;
t320 = t216 - t324;
t319 = -0.2e1 * t264;
t272 = qJD(3) * t211;
t260 = -pkin(7) * t204 - pkin(3);
t219 = -t263 + (-pkin(4) + t260) * t209;
t160 = t234 * qJD(1);
t292 = t205 * t160;
t64 = qJD(1) * t219 - t292;
t317 = t204 * t272 - t211 * t64;
t316 = t191 * t229;
t315 = t191 * t230;
t271 = qJD(4) * t204;
t283 = qJ(4) * t254 - t195;
t242 = -t309 * t256 + t271 - t283;
t147 = t151 ^ 2;
t314 = -t149 ^ 2 - t147;
t308 = -pkin(8) + qJ(3);
t171 = t308 * t204;
t172 = t308 * t205;
t284 = t208 * t171 + t211 * t172;
t265 = t212 * qJD(4);
t277 = qJD(2) * t209;
t313 = qJ(4) * t277 - t265;
t273 = qJD(3) * t208;
t137 = t204 * t160;
t192 = qJ(4) * t280;
t291 = t205 * t209;
t294 = t204 * t212;
t223 = -pkin(7) * t291 + pkin(8) * t294;
t82 = qJD(1) * t223 + t137 + t192;
t312 = -t171 * t269 + t172 * t270 - t204 * t273 - t205 * t272 + t208 * t64 + t211 * t82;
t311 = qJD(1) * t277;
t310 = t180 ^ 2;
t295 = t204 * t211;
t157 = -t205 * t208 + t295;
t95 = t210 * t156 + t157 * t207;
t307 = -qJD(6) * t95 - t285 * t207 + t286 * t210;
t96 = -t156 * t207 + t157 * t210;
t306 = qJD(6) * t96 + t286 * t207 + t285 * t210;
t188 = pkin(7) * t294;
t201 = t212 * pkin(3);
t87 = t212 * pkin(4) + t188 + t201 + (-pkin(8) * t209 - t167) * t205;
t189 = pkin(7) * t290;
t118 = t204 * t167 + t189;
t109 = -qJ(4) * t212 + t118;
t296 = t204 * t209;
t97 = pkin(8) * t296 + t109;
t304 = t208 * t87 + t211 * t97;
t302 = t209 * t81;
t301 = t209 * t84;
t300 = t285 * pkin(5) + t242;
t297 = t191 * t212;
t293 = t205 * t132;
t214 = qJD(1) ^ 2;
t289 = t212 * t214;
t213 = qJD(2) ^ 2;
t288 = t213 * t209;
t287 = t213 * t212;
t253 = t212 * t266;
t282 = -qJ(4) * t253 - qJD(4) * t291;
t203 = t212 ^ 2;
t281 = t209 ^ 2 - t203;
t276 = qJD(2) * t212;
t275 = qJD(3) * t151;
t274 = qJD(3) * t205;
t165 = -t205 * pkin(3) - t204 * qJ(4) - pkin(2);
t262 = pkin(7) * t277;
t259 = pkin(3) * t204 + pkin(7);
t250 = qJD(6) * t10 + t5;
t57 = qJD(2) * t219 - t293;
t123 = t204 * t132;
t58 = qJD(2) * t223 + t123 + t313;
t247 = -t208 * t58 + t211 * t57;
t246 = -t208 * t97 + t211 * t87;
t245 = -t166 - t206;
t244 = pkin(1) * t319;
t243 = t211 * t171 - t172 * t208;
t117 = t167 * t205 - t188;
t138 = t205 * pkin(4) - t165;
t115 = pkin(3) * t256 - t283;
t241 = t115 + t271;
t190 = pkin(7) * t252;
t67 = pkin(3) * t176 - qJ(4) * t240 - t151 * qJD(4) + t190;
t239 = t309 * t204 - pkin(7);
t70 = -pkin(9) * t156 + t284;
t238 = -pkin(5) * t280 + t286 * pkin(9) + t284 * qJD(5) + qJD(6) * t70 + t205 * t273 - t208 * t82 - t317;
t69 = -pkin(9) * t157 + t243;
t237 = t285 * pkin(9) - qJD(6) * t69 + t312;
t236 = t260 * t209;
t235 = -t151 * t279 + t176;
t128 = t156 * t209;
t21 = pkin(5) * t212 - pkin(9) * t128 + t246;
t127 = t208 * t291 - t209 * t295;
t22 = -pkin(9) * t127 + t304;
t233 = t207 * t21 + t210 * t22;
t65 = t210 * t127 + t128 * t207;
t66 = -t127 * t207 + t128 * t210;
t228 = t207 * t211 + t208 * t210;
t227 = t207 * t208 - t210 * t211;
t112 = -pkin(7) * t255 + t137;
t106 = -t205 * t262 + t123;
t225 = t208 * t57 + t211 * t58 + t87 * t269 - t97 * t270;
t186 = qJ(4) * t291;
t108 = t209 * t239 + t186;
t80 = t239 * t276 - t282;
t51 = -pkin(4) * t176 - t67;
t175 = qJD(3) * t256;
t130 = t149 * t279;
t129 = t149 * t274;
t124 = t259 * t209 - t186;
t111 = pkin(7) * t257 + t292;
t110 = -t117 + t201;
t107 = t130 + t240;
t105 = t204 * t262 + t293;
t104 = pkin(5) * t156 + t138;
t103 = qJD(1) * t236 - t292;
t102 = t112 + t192;
t101 = t259 * t276 + t282;
t89 = qJD(2) * t236 - t293;
t77 = t106 + t313;
t73 = qJD(5) * t157 * t209 + qJD(2) * t222;
t72 = qJD(5) * t128 + t208 * t253 - t276 * t295;
t62 = t127 * pkin(5) + t108;
t61 = -pkin(3) * t251 - t78;
t49 = -qJD(1) * t265 + t261;
t24 = t72 * pkin(5) + t80;
t19 = pkin(5) * t46 + t51;
t15 = qJD(6) * t66 + t207 * t73 + t210 * t72;
t14 = -qJD(6) * t65 - t207 * t72 + t210 * t73;
t7 = -pkin(9) * t72 + t225;
t6 = -pkin(5) * t277 - t73 * pkin(9) - qJD(5) * t304 + t247;
t2 = t10 * t210 - t13 * t207;
t1 = [0, 0, 0, 0.2e1 * t212 * t251, t281 * t319, t287, -t288, 0, -pkin(7) * t287 + t209 * t244, pkin(7) * t288 + t212 * t244 (-qJD(1) * t105 - t78) * t212 + ((pkin(7) * t149 + t166 * t204) * t212 + (t98 + (t117 + 0.2e1 * t188) * qJD(1)) * t209) * qJD(2) (qJD(1) * t106 + t79) * t212 + ((pkin(7) * t151 + t166 * t205) * t212 + (-t99 + (-t118 + 0.2e1 * t189) * qJD(1)) * t209) * qJD(2), -t105 * t151 - t106 * t149 + (-t204 * t79 - t205 * t78) * t209 + (-t204 * t99 - t205 * t98 + (-t117 * t205 - t118 * t204) * qJD(1)) * t276, t98 * t105 + t99 * t106 + t78 * t117 + t79 * t118 + (t166 + t194) * pkin(7) * t276, t67 * t296 + t101 * t149 + (qJD(1) * t89 + t61) * t212 + (t74 * t294 - t302 + (-t110 * t209 + t124 * t294) * qJD(1)) * qJD(2), -t77 * t149 + t89 * t151 + (-t204 * t49 + t205 * t61) * t209 + (-t204 * t84 + t205 * t81 + (-t109 * t204 + t110 * t205) * qJD(1)) * t276, -t67 * t291 - t101 * t151 + (-qJD(1) * t77 - t49) * t212 + (-t74 * t290 + t301 + (t109 * t209 - t124 * t290) * qJD(1)) * qJD(2), t101 * t74 + t109 * t49 + t110 * t61 + t124 * t67 + t77 * t84 + t81 * t89, t128 * t45 + t229 * t73, -t127 * t45 - t128 * t46 - t229 * t72 - t230 * t73, t73 * t191 + t45 * t212 + (-qJD(1) * t128 - t229) * t277, -t72 * t191 - t46 * t212 + (qJD(1) * t127 + t230) * t277 (-t191 - t279) * t277, t247 * t191 + t248 * t212 + t80 * t230 + t108 * t46 + t51 * t127 + t47 * t72 + (-t17 * t212 - t191 * t304) * qJD(5) + (-qJD(1) * t246 - t16) * t277, -t225 * t191 - t226 * t212 + t80 * t229 + t108 * t45 + t51 * t128 + t47 * t73 + (t304 * qJD(1) + t17) * t277, -t14 * t232 - t224 * t66, -t14 * t26 + t15 * t232 + t216 * t66 + t224 * t65, t14 * t180 - t224 * t212 + (-qJD(1) * t66 + t232) * t277, -t15 * t180 + t216 * t212 + (qJD(1) * t65 + t26) * t277 (-t180 - t279) * t277 (-t207 * t7 + t210 * t6) * t180 + t258 * t212 + t24 * t26 - t62 * t216 + t19 * t65 + t23 * t15 + (-t180 * t233 - t212 * t3) * qJD(6) + (-(-t207 * t22 + t21 * t210) * qJD(1) - t2) * t277, t11 * t212 + t23 * t14 + t19 * t66 - t24 * t232 - t62 * t224 + (-(-qJD(6) * t22 + t6) * t180 - t4 * t212) * t207 + (-(qJD(6) * t21 + t7) * t180 - t250 * t212) * t210 + (qJD(1) * t233 + t3) * t277; 0, 0, 0, -t209 * t289, t281 * t214, 0, 0, 0, t214 * pkin(1) * t209, pkin(1) * t289, t175 + ((-qJ(3) * t278 - t98) * t209 + (t111 + t245 * t204 + (-t149 - t266) * pkin(7)) * t212) * qJD(1) ((-qJ(3) * t266 + t99) * t209 + (-t112 + (-t151 + t278) * pkin(7) + (qJD(3) + t245) * t205) * t212) * qJD(1), t111 * t151 + t112 * t149 - t129 + (t98 * t279 + t79) * t205 + (t99 * t279 + t275 - t78) * t204, -t98 * t111 - t99 * t112 + (-t204 * t98 + t205 * t99) * qJD(3) + (-t78 * t204 + t79 * t205) * qJ(3) + t245 * t195, -t205 * t67 + t175 - t241 * t149 + (-t103 * t212 + t302 + (-t212 * t74 + (t165 * t212 - t299) * qJD(2)) * t204) * qJD(1), t102 * t149 - t103 * t151 - t129 + (-t279 * t81 + t49) * t205 + (t279 * t84 + t275 + t61) * t204, -t204 * t67 + t241 * t151 + (t102 * t212 - t301 + (qJ(3) * t277 + (-qJD(2) * t165 - qJD(3) + t74) * t212) * t205) * qJD(1), -t102 * t84 - t103 * t81 - t115 * t74 + t165 * t67 + (qJ(3) * t49 + qJD(3) * t84) * t205 + (qJ(3) * t61 + qJD(3) * t81 - qJD(4) * t74) * t204, t45 * t157 + t229 * t286, -t45 * t156 - t157 * t46 - t229 * t285 - t230 * t286, t286 * t191 + (-qJD(2) * t157 + t229) * t280, -t285 * t191 + (qJD(2) * t156 - t230) * t280, t191 * t280, t138 * t46 + t51 * t156 + t242 * t230 + t285 * t47 + (-t172 * t269 + (-qJD(5) * t171 - t274 + t82) * t208 + t317) * t191 + (-qJD(2) * t243 + t16) * t280, t138 * t45 + t51 * t157 + t242 * t229 + t286 * t47 + t312 * t191 + (qJD(2) * t284 - t17) * t280, -t224 * t96 - t232 * t307, t216 * t96 + t224 * t95 + t232 * t306 - t307 * t26, t307 * t180 + (-qJD(2) * t96 - t232) * t280, -t306 * t180 + (qJD(2) * t95 - t26) * t280, t180 * t280, -t104 * t216 + t19 * t95 + t300 * t26 + t306 * t23 + (t207 * t237 - t210 * t238) * t180 + (-(-t207 * t70 + t210 * t69) * qJD(2) + t2) * t280, -t104 * t224 + t19 * t96 - t300 * t232 + t307 * t23 + (t207 * t238 + t210 * t237) * t180 + ((t207 * t69 + t210 * t70) * qJD(2) - t3) * t280; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t235, t107, t314, t149 * t99 + t151 * t98 + t190, t235, t314, -t107, t149 * t84 - t151 * t81 + t67, 0, 0, 0, 0, 0, -t46 - t316, -t45 + t315, 0, 0, 0, 0, 0, t216 + t324, t224 + t318; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t149 * t151 - t251, -t130 + t240, -t203 * t214 - t147, t74 * t151 + (-pkin(3) * t277 + t212 * t84) * qJD(1) - t78, 0, 0, 0, 0, 0, -t191 * t270 - t151 * t230 + (-t208 * t297 - t211 * t277) * qJD(1), -t191 * t269 - t151 * t229 + (t208 * t277 - t211 * t297) * qJD(1), 0, 0, 0, 0, 0, -t151 * t26 + t227 * t311 - t228 * t310, t151 * t232 + t227 * t310 + t228 * t311; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229 * t230, t229 ^ 2 - t230 ^ 2, t45 + t315, -t46 + t316, -t251, t17 * t191 - t229 * t47 + t217, t16 * t191 + t230 * t47 - t226, -t325, t323, t326, t320, -t251 -(-t12 * t207 - t303) * t180 + (-t180 * t268 - t210 * t251 - t229 * t26) * pkin(5) + t321 (-t13 * t180 - t4) * t207 + (t12 * t180 - t250) * t210 + (-t180 * t267 + t207 * t251 + t229 * t232) * pkin(5) + t322; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t325, t323, t326, t320, -t251, t3 * t180 + t321, t2 * t180 - t207 * t4 - t210 * t250 + t322;];
tauc_reg  = t1;
