% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 00:27
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRRRP5_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP5_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRRP5_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 00:25:53
% EndTime: 2019-03-09 00:26:12
% DurationCPUTime: 7.11s
% Computational Cost: add. (6350->454), mult. (17148->663), div. (0->0), fcn. (13950->12), ass. (0->230)
t177 = sin(pkin(7));
t183 = sin(qJ(3));
t289 = t177 * t183;
t170 = pkin(9) * t289;
t179 = cos(pkin(7));
t187 = cos(qJ(3));
t188 = cos(qJ(2));
t278 = t187 * t188;
t184 = sin(qJ(2));
t283 = t183 * t184;
t203 = -t179 * t283 + t278;
t178 = sin(pkin(6));
t270 = qJD(1) * t178;
t285 = t179 * t187;
t277 = t203 * t270 - (pkin(2) * t285 - t170) * qJD(3);
t216 = pkin(3) * t183 - pkin(10) * t187;
t201 = t216 * qJD(3);
t248 = t184 * t270;
t331 = (t201 - t248) * t177;
t264 = qJD(2) * t187;
t169 = t177 * t264;
t213 = t169 - qJD(4);
t288 = t177 * t187;
t254 = pkin(9) * t288;
t132 = t254 + (pkin(2) * t183 + pkin(10)) * t179;
t217 = -pkin(3) * t187 - pkin(10) * t183;
t133 = (-pkin(2) + t217) * t177;
t182 = sin(qJ(4));
t186 = cos(qJ(4));
t259 = qJD(4) * t186;
t261 = qJD(4) * t182;
t330 = -t132 * t261 + t133 * t259 + t331 * t182 - t186 * t277;
t281 = t184 * t187;
t282 = t183 * t188;
t205 = t179 * t281 + t282;
t286 = t179 * t183;
t275 = -t205 * t270 + (pkin(2) * t286 + t254) * qJD(3);
t263 = qJD(3) * t183;
t246 = t177 * t263;
t329 = -pkin(11) * t246 - t330;
t144 = -t186 * t179 + t182 * t289;
t262 = qJD(3) * t187;
t245 = t177 * t262;
t107 = -t144 * qJD(4) + t186 * t245;
t145 = t179 * t182 + t186 * t289;
t108 = t145 * qJD(4) + t182 * t245;
t328 = -pkin(4) * t108 + pkin(11) * t107 - t275;
t180 = cos(pkin(6));
t269 = qJD(1) * t180;
t249 = t177 * t269;
t158 = t183 * t249;
t268 = qJD(2) * t177;
t150 = pkin(9) * t268 + t248;
t304 = qJD(2) * pkin(2);
t157 = t188 * t270 + t304;
t319 = t187 * t150 + t157 * t286;
t86 = t158 + t319;
t327 = -t86 - t213 * (pkin(4) * t182 - pkin(11) * t186);
t266 = qJD(2) * t179;
t229 = qJD(3) + t266;
t265 = qJD(2) * t183;
t247 = t177 * t265;
t192 = -t182 * t247 + t186 * t229;
t255 = qJD(2) * qJD(3);
t240 = t177 * t255;
t218 = t187 * t240;
t210 = t186 * t218;
t190 = qJD(4) * t192 + t210;
t326 = -qJD(5) * t213 + t190;
t312 = -qJ(6) - pkin(11);
t325 = qJ(6) * t192 + qJD(5) * t312;
t181 = sin(qJ(5));
t185 = cos(qJ(5));
t257 = qJD(5) * t185;
t258 = qJD(5) * t181;
t131 = t170 + (-pkin(2) * t187 - pkin(3)) * t179;
t77 = pkin(4) * t144 - pkin(11) * t145 + t131;
t276 = t186 * t132 + t182 * t133;
t79 = -pkin(11) * t288 + t276;
t324 = t328 * t181 + t185 * t329 - t77 * t257 + t79 * t258;
t323 = t132 * t259 + t133 * t261 - t277 * t182 - t186 * t331;
t305 = t181 * t77 + t185 * t79;
t322 = -t305 * qJD(5) + t181 * t329 - t328 * t185;
t253 = pkin(10) * t261;
t321 = t181 * t253 + t185 * t327;
t320 = t183 * t187;
t223 = t186 * t169;
t318 = t223 - t259;
t164 = -pkin(4) * t186 - pkin(11) * t182 - pkin(3);
t134 = t216 * t268;
t139 = t183 * t150;
t291 = t157 * t179;
t85 = (t249 + t291) * t187 - t139;
t294 = t182 * t134 + t186 * t85;
t51 = pkin(11) * t247 + t294;
t317 = -t164 * t257 - t181 * t327 + t185 * t51;
t126 = t182 * t229 + t186 * t247;
t97 = t185 * t126 - t181 * t213;
t316 = t97 ^ 2;
t189 = qJD(2) ^ 2;
t252 = t181 * t288;
t110 = t145 * t185 - t252;
t279 = t185 * t187;
t109 = t145 * t181 + t177 * t279;
t56 = t109 * qJD(5) - t185 * t107 - t181 * t246;
t315 = pkin(5) * t108 + qJ(6) * t56 - qJD(6) * t110 + t322;
t57 = -qJD(5) * t252 + t107 * t181 + t145 * t257 - t185 * t246;
t314 = -qJ(6) * t57 - qJD(6) * t109 - t324;
t313 = pkin(5) * t181;
t74 = t229 * pkin(10) + t86;
t168 = t179 * t269;
t92 = t168 + (t217 * qJD(2) - t157) * t177;
t36 = t182 * t92 + t186 * t74;
t31 = -t213 * pkin(11) + t36;
t73 = -t229 * pkin(3) - t85;
t38 = -pkin(4) * t192 - t126 * pkin(11) + t73;
t16 = -t181 * t31 + t185 * t38;
t11 = -qJ(6) * t97 + t16;
t120 = qJD(5) - t192;
t10 = pkin(5) * t120 + t11;
t311 = t10 - t11;
t35 = -t182 * t74 + t186 * t92;
t84 = pkin(4) * t126 - pkin(11) * t192;
t310 = t181 * t84 + t185 * t35;
t309 = -pkin(4) * t246 + t323;
t118 = (t181 * t183 + t186 * t279) * t268;
t280 = t185 * t186;
t172 = pkin(10) * t280;
t224 = t182 * t169;
t256 = qJD(6) * t185;
t292 = qJ(6) * t182;
t307 = -pkin(5) * t224 + qJ(6) * t118 + t181 * t51 - t182 * t256 + (pkin(5) * t182 - qJ(6) * t280) * qJD(4) + (-t172 + (-t164 + t292) * t181) * qJD(5) + t321;
t117 = t181 * t223 - t185 * t247;
t284 = t182 * t185;
t306 = qJ(6) * t117 + (-pkin(10) * qJD(4) - qJ(6) * qJD(5)) * t284 + (-qJD(6) * t182 + (-pkin(10) * qJD(5) - qJ(6) * qJD(4)) * t186) * t181 - t317;
t95 = t126 * t181 + t185 * t213;
t303 = t120 * t95;
t219 = t183 * t240;
t102 = (t201 + t248) * t268;
t195 = t203 * qJD(2);
t221 = t180 * t245;
t53 = (t157 * t285 - t139) * qJD(3) + (t178 * t195 + t221) * qJD(1);
t237 = t186 * t102 - t182 * t53 - t74 * t259 - t92 * t261;
t15 = -pkin(4) * t219 - t237;
t302 = t15 * t181;
t301 = t15 * t185;
t43 = t126 * t258 - t181 * t219 - t185 * t326;
t300 = t181 * t43;
t211 = t182 * t218;
t99 = t126 * qJD(4) + t211;
t299 = t181 * t99;
t298 = t185 * t99;
t297 = t97 * t120;
t296 = t181 * t325 + t256 - t310;
t81 = t185 * t84;
t295 = -pkin(5) * t126 - t81 + t325 * t185 + (-qJD(6) + t35) * t181;
t174 = t177 ^ 2;
t290 = t174 * t189;
t287 = t178 * t189;
t272 = t181 * t164 + t172;
t271 = t183 ^ 2 - t187 ^ 2;
t267 = qJD(2) * t178;
t260 = qJD(4) * t185;
t251 = t184 * t287;
t250 = pkin(10) + t313;
t244 = t120 * t258;
t243 = t177 * t179 * t189;
t242 = qJD(3) * t291;
t241 = qJD(1) * t267;
t239 = -t181 * t79 + t185 * t77;
t236 = -t182 * t132 + t133 * t186;
t235 = t187 * t213;
t234 = t120 * t185;
t233 = t213 * t177;
t232 = qJD(4) * t213;
t230 = 0.2e1 * t174 * t255;
t228 = qJD(3) + 0.2e1 * t266;
t227 = t174 * t251;
t225 = t177 * t184 * t267;
t222 = t180 * t246;
t78 = pkin(4) * t288 - t236;
t214 = t185 * t259 - t118;
t17 = t181 * t38 + t185 * t31;
t206 = t179 * t278 - t283;
t105 = -t206 * t178 - t180 * t288;
t204 = t179 * t282 + t281;
t106 = t178 * t204 + t180 * t289;
t143 = -t177 * t178 * t188 + t179 * t180;
t76 = t106 * t186 + t143 * t182;
t41 = t105 * t185 - t181 * t76;
t42 = t105 * t181 + t185 * t76;
t75 = t106 * t182 - t143 * t186;
t212 = t179 * t184 * t241;
t119 = -t157 * t177 + t168;
t208 = t119 * t177 - t174 * t304;
t207 = -t120 * t257 - t299;
t198 = -t182 * t102 - t186 * t53 - t92 * t259 + t74 * t261;
t14 = pkin(11) * t219 - t198;
t26 = t99 * pkin(4) - pkin(11) * t190 + qJD(3) * t158 + t150 * t262 + t183 * t242 + t187 * t212 + t241 * t282;
t3 = t185 * t14 + t181 * t26 + t38 * t257 - t31 * t258;
t30 = t213 * pkin(4) - t35;
t200 = -pkin(11) * t99 + t120 * t30;
t196 = t205 * qJD(2);
t82 = t182 * t85;
t50 = -pkin(4) * t247 - t134 * t186 + t82;
t194 = qJD(3) * t150 + t212;
t4 = -t17 * qJD(5) - t14 * t181 + t185 * t26;
t44 = t126 * t257 + t181 * t326 - t185 * t219;
t7 = pkin(5) * t44 + t15;
t191 = -t119 * t268 - t242 + (-qJD(3) * t177 * t180 - t188 * t267) * qJD(1);
t166 = t312 * t185;
t165 = t312 * t181;
t152 = t185 * t164;
t111 = -t181 * t292 + t272;
t104 = -qJ(6) * t284 + t152 + (-pkin(10) * t181 - pkin(5)) * t186;
t94 = t95 ^ 2;
t69 = t221 + (t206 * qJD(3) + t195) * t178;
t68 = t222 + (qJD(3) * t204 + t196) * t178;
t54 = t319 * qJD(3) + (t178 * t196 + t222) * qJD(1);
t29 = -t75 * qJD(4) + t182 * t225 + t186 * t69;
t28 = t76 * qJD(4) + t182 * t69 - t186 * t225;
t27 = -qJ(6) * t109 + t305;
t23 = t95 * pkin(5) + qJD(6) + t30;
t22 = pkin(5) * t144 - qJ(6) * t110 + t239;
t12 = -qJ(6) * t95 + t17;
t9 = t41 * qJD(5) + t181 * t68 + t185 * t29;
t8 = -t42 * qJD(5) - t181 * t29 + t185 * t68;
t2 = -qJ(6) * t44 - qJD(6) * t95 + t3;
t1 = pkin(5) * t99 + qJ(6) * t43 - qJD(6) * t97 + t4;
t5 = [0, 0, -t251, -t188 * t287, 0, 0, 0, 0, 0, t143 * t219 - t187 * t227 - t68 * t229, t143 * t218 + t183 * t227 - t69 * t229, 0, 0, 0, 0, 0, t105 * t99 - t192 * t68 + t28 * t213 - t75 * t219, t105 * t190 + t68 * t126 + t29 * t213 - t76 * t219, 0, 0, 0, 0, 0, t120 * t8 + t28 * t95 + t41 * t99 + t44 * t75, -t120 * t9 + t28 * t97 - t42 * t99 - t43 * t75, t41 * t43 - t42 * t44 - t8 * t97 - t9 * t95, t1 * t41 + t10 * t8 + t12 * t9 + t2 * t42 + t23 * t28 + t7 * t75; 0, 0, 0, 0, t230 * t320, -t271 * t230, t228 * t245, -t228 * t246, 0 (-t275 * qJD(2) - t54) * t179 + (t208 * t183 - t275) * qJD(3) (t277 * qJD(2) - t53) * t179 + (t208 * t187 + t277) * qJD(3), t126 * t107 + t145 * t190, t107 * t192 - t126 * t108 - t144 * t190 - t145 * t99, -t107 * t213 + t126 * t246 + t145 * t219 - t190 * t288, t108 * t213 + (t99 * t187 + (-qJD(2) * t144 + t192) * t263) * t177 (-t174 * t264 - t233) * t263, t131 * t99 + t54 * t144 + t73 * t108 - t275 * t192 + (-t237 * t187 + (t236 * qJD(2) + t35) * t263) * t177 + t323 * t213, t73 * t107 + t275 * t126 + t131 * t190 + t54 * t145 - t198 * t288 + t213 * t330 - t276 * t219 - t36 * t246, -t110 * t43 - t56 * t97, t109 * t43 - t110 * t44 + t56 * t95 - t57 * t97, t108 * t97 + t110 * t99 - t120 * t56 - t144 * t43, -t108 * t95 - t109 * t99 - t120 * t57 - t144 * t44, t108 * t120 + t144 * t99, t16 * t108 + t15 * t109 + t322 * t120 + t4 * t144 + t239 * t99 + t30 * t57 + t309 * t95 + t78 * t44, -t17 * t108 + t15 * t110 + t324 * t120 - t3 * t144 - t30 * t56 - t305 * t99 + t309 * t97 - t78 * t43, -t1 * t110 + t10 * t56 - t109 * t2 - t12 * t57 + t22 * t43 - t27 * t44 - t314 * t95 - t315 * t97, t2 * t27 + t1 * t22 + t7 * (pkin(5) * t109 + t78) + (pkin(5) * t57 + t309) * t23 + t314 * t12 + t315 * t10; 0, 0, 0, 0, -t290 * t320, t271 * t290, -t187 * t243, t183 * t243, 0, t191 * t183 - t194 * t187 + t86 * t229, t194 * t183 + t191 * t187 + t85 * t229, -qJD(4) * t182 ^ 2 * t247 + ((qJD(4) * t229 + t218) * t182 - t213 * t126) * t186, -t182 * t99 + t186 * t190 + (t224 - t261) * t126 - t318 * t192, -t186 * t232 + (t186 * t235 + (qJD(3) * t182 - t126) * t183) * t268, t182 * t232 + (-t182 * t235 + (t186 * qJD(3) - t192) * t183) * t268, t233 * t265, -pkin(3) * t99 + t73 * t261 - t82 * t213 + t86 * t192 + (pkin(10) * t232 + t134 * t213 - t54) * t186 + (-t183 * t35 + (-pkin(10) * t263 - t187 * t73) * t182) * t268, -t186 * pkin(10) * t219 - pkin(3) * t190 - t86 * t126 + t54 * t182 + t36 * t247 - t318 * t73 + (-t253 - t294) * t213, -t43 * t284 + (-t182 * t258 + t214) * t97, t117 * t97 + t118 * t95 + (-t181 * t97 - t185 * t95) * t259 + (t300 - t185 * t44 + (t181 * t95 - t185 * t97) * qJD(5)) * t182, t186 * t43 + t214 * t120 + (-t213 * t97 - t244 + t298) * t182, t186 * t44 + (-t181 * t259 + t117) * t120 + (t213 * t95 + t207) * t182, -t120 * t182 * t213 - t186 * t99, -t30 * t117 + t152 * t99 - t50 * t95 + ((-qJD(5) * t164 + t51) * t181 + t321) * t120 + (t30 * t181 * qJD(4) - t4 + (qJD(4) * t95 + t207) * pkin(10)) * t186 + (pkin(10) * t44 - t16 * t213 + t257 * t30 + t302) * t182, -t272 * t99 - t50 * t97 - t30 * t118 + t317 * t120 + (t30 * t260 + t3 + (qJD(4) * t97 + t244) * pkin(10)) * t186 + (-t30 * t258 + t301 + t213 * t17 + (t120 * t260 - t43) * pkin(10)) * t182, t10 * t118 + t104 * t43 - t111 * t44 + t117 * t12 - t307 * t97 - t306 * t95 + (-t10 * t185 - t12 * t181) * t259 + (-t1 * t185 - t181 * t2 + (t10 * t181 - t12 * t185) * qJD(5)) * t182, t250 * t7 * t182 + t1 * t104 + t307 * t10 + t2 * t111 + t306 * t12 + (t250 * t259 - t50 + (t257 * t182 - t117) * pkin(5)) * t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126 * t192, t126 ^ 2 - t192 ^ 2, t192 * t169 + t210, -t126 * t169 - t211, t219, -t73 * t126 - t36 * t213 + t237, -t192 * t73 - t35 * t213 + t198, t97 * t234 - t300 (-t43 - t303) * t185 + (-t44 - t297) * t181, t120 * t234 - t126 * t97 + t299, -t120 ^ 2 * t181 + t126 * t95 + t298, -t120 * t126, -pkin(4) * t44 - t16 * t126 - t301 - t36 * t95 + (-pkin(11) * t257 - t81) * t120 + (t35 * t120 + t200) * t181, pkin(4) * t43 + t17 * t126 + t302 - t36 * t97 + (pkin(11) * t258 + t310) * t120 + t200 * t185, t165 * t43 + t166 * t44 - t295 * t97 - t296 * t95 + (-t10 * t120 + t2) * t185 + (-t12 * t120 - t1) * t181, -t2 * t166 + t1 * t165 + t7 * (-pkin(5) * t185 - pkin(4)) + (t120 * t313 - t36) * t23 + t296 * t12 + t295 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97 * t95, -t94 + t316, -t43 + t303, -t44 + t297, t99, t120 * t17 - t30 * t97 + t4, t120 * t16 + t30 * t95 - t3, pkin(5) * t43 - t311 * t95, t311 * t12 + (-t23 * t97 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94 - t316, t10 * t97 + t12 * t95 + t7;];
tauc_reg  = t5;
