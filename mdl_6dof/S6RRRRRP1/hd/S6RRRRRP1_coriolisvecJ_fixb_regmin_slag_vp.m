% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
% 
% Output:
% tauc_reg [6x33]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRRRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:58:19
% EndTime: 2019-03-10 00:58:31
% DurationCPUTime: 5.18s
% Computational Cost: add. (10454->421), mult. (26344->550), div. (0->0), fcn. (19685->8), ass. (0->246)
t207 = cos(qJ(5));
t279 = qJD(5) * t207;
t205 = sin(qJ(3));
t206 = sin(qJ(2));
t287 = qJD(1) * t206;
t267 = t205 * t287;
t209 = cos(qJ(3));
t210 = cos(qJ(2));
t286 = qJD(1) * t210;
t269 = t209 * t286;
t156 = -t267 + t269;
t157 = -t205 * t286 - t209 * t287;
t204 = sin(qJ(4));
t208 = cos(qJ(4));
t349 = t208 * t156 + t204 * t157;
t358 = t349 * t207;
t364 = t279 - t358;
t203 = sin(qJ(5));
t236 = t204 * t156 - t208 * t157;
t200 = qJD(2) + qJD(3);
t260 = qJD(4) + t200;
t108 = t203 * t236 - t207 * t260;
t110 = t203 * t260 + t207 * t236;
t277 = qJD(5) - t349;
t363 = t277 * t203;
t280 = qJD(5) * t203;
t298 = t205 * t210;
t171 = t206 * t209 + t298;
t344 = qJD(1) * t171;
t216 = t200 * t344;
t215 = t204 * t216;
t276 = qJD(1) * qJD(2);
t266 = t210 * t276;
t132 = qJD(3) * t269 - t200 * t267 + t209 * t266;
t296 = t208 * t132;
t214 = -t215 + t296;
t213 = t349 * qJD(4) + t214;
t357 = -qJD(5) * t260 - t213;
t42 = t357 * t207 + t236 * t280;
t43 = t110 * qJD(5) + t213 * t203;
t6 = -t364 * t108 - t110 * t363 - t203 * t43 - t42 * t207;
t257 = t204 * t132 + t208 * t216;
t55 = qJD(4) * t236 + t257;
t52 = t207 * t55;
t13 = t108 * t236 - t277 * t363 + t52;
t340 = pkin(7) + pkin(8);
t270 = qJD(2) * t340;
t245 = qJD(1) * t270;
t166 = t210 * t245;
t181 = t340 * t210;
t174 = qJD(1) * t181;
t284 = qJD(3) * t205;
t254 = -t205 * t166 - t174 * t284;
t180 = t340 * t206;
t172 = qJD(1) * t180;
t164 = qJD(2) * pkin(2) - t172;
t165 = t206 * t245;
t350 = t209 * (qJD(3) * t164 - t165);
t63 = -pkin(9) * t216 + t254 + t350;
t153 = t157 * pkin(9);
t158 = t205 * t174;
t256 = t209 * t164 - t158;
t105 = t153 + t256;
t96 = pkin(3) * t200 + t105;
t259 = qJD(4) * t96 + t63;
t162 = t209 * t174;
t235 = -t164 * t205 - t162;
t337 = pkin(9) * t156;
t106 = -t235 + t337;
t282 = qJD(4) * t204;
t255 = t205 * t165 - t209 * t166;
t221 = qJD(3) * t235 + t255;
t64 = -pkin(9) * t132 + t221;
t262 = -t106 * t282 + t204 * t64;
t11 = t208 * t259 + t262;
t194 = -pkin(2) * t210 - pkin(1);
t179 = t194 * qJD(1);
t140 = -pkin(3) * t156 + t179;
t304 = t140 * t349;
t362 = -t11 - t304;
t253 = t172 * t205 - t162;
t113 = t253 - t337;
t290 = -t209 * t172 - t158;
t114 = t153 + t290;
t193 = pkin(2) * t209 + pkin(3);
t281 = qJD(4) * t208;
t300 = t204 * t205;
t351 = t113 * t204 + t114 * t208 - t193 * t281 - (-t205 * t282 + (t208 * t209 - t300) * qJD(3)) * pkin(2);
t103 = t204 * t106;
t66 = t105 * t208 - t103;
t360 = -pkin(3) * t281 + t66;
t359 = t349 * t203;
t348 = qJ(6) * t359 + t207 * qJD(6);
t40 = t42 * t203;
t17 = t364 * t110 - t40;
t320 = t203 * t55 + t277 * t279;
t14 = -t110 * t236 - t277 * t358 + t320;
t307 = t236 * t349;
t299 = t205 * t208;
t356 = (t205 * t281 + (t204 * t209 + t299) * qJD(3)) * pkin(2) + t113 * t208;
t47 = t236 ^ 2 - t349 ^ 2;
t91 = pkin(4) * t236 - pkin(10) * t349;
t199 = t207 * qJ(6);
t243 = pkin(5) * t236 - t199 * t349;
t45 = -t349 * t200 + t214;
t354 = -0.2e1 * t276;
t319 = -t114 * t204 + t193 * t282 + t356;
t12 = t106 * t281 + t204 * t63 - t208 * t64 + t96 * t282;
t57 = t208 * t96 - t103;
t53 = -pkin(4) * t260 - t57;
t264 = -t12 * t207 + t53 * t280;
t312 = t277 * t236;
t119 = -pkin(9) * t171 - t180 * t209 - t181 * t205;
t170 = t205 * t206 - t209 * t210;
t234 = t180 * t205 - t181 * t209;
t120 = -pkin(9) * t170 - t234;
t87 = -t119 * t208 + t204 * t120;
t223 = t156 * t281 + t157 * t282 + t296;
t347 = -t215 + t223;
t195 = pkin(2) * t287;
t339 = pkin(3) * t157;
t79 = -t339 + t91;
t74 = t195 + t79;
t346 = t203 * t74 + t207 * t351;
t345 = t203 * t79 + t207 * t360;
t231 = -t140 * t236 - t12;
t104 = t208 * t106;
t58 = t204 * t96 + t104;
t54 = pkin(10) * t260 + t58;
t71 = -pkin(4) * t349 - pkin(10) * t236 + t140;
t36 = -t203 * t54 + t207 * t71;
t233 = -t36 * t236 + t264;
t37 = t203 * t71 + t207 * t54;
t246 = t12 * t203 + t37 * t236 + t53 * t279;
t46 = t236 * t200 - t257;
t341 = t110 ^ 2;
t338 = pkin(5) * t207;
t336 = t208 * pkin(3);
t335 = -qJ(6) - pkin(10);
t22 = -qJ(6) * t110 + t36;
t16 = pkin(5) * t277 + t22;
t333 = t16 - t22;
t332 = t203 * t91 + t207 * t57;
t88 = t119 * t204 + t120 * t208;
t85 = t207 * t88;
t134 = t208 * t170 + t171 * t204;
t135 = -t170 * t204 + t171 * t208;
t144 = pkin(3) * t170 + t194;
t86 = pkin(4) * t134 - pkin(10) * t135 + t144;
t329 = t203 * t86 + t85;
t150 = pkin(2) * t299 + t193 * t204 + pkin(10);
t292 = -qJ(6) - t150;
t252 = qJD(5) * t292;
t328 = t203 * t252 - t346 + t348;
t70 = t207 * t74;
t327 = t207 * t252 - t243 - t70 + (-qJD(6) + t351) * t203;
t326 = pkin(3) * qJD(4);
t324 = t16 * t207;
t138 = t200 * t170;
t225 = t171 * qJD(3);
t139 = qJD(2) * t171 + t225;
t75 = -qJD(4) * t134 - t138 * t208 - t139 * t204;
t323 = t207 * t75;
t322 = t53 * t349;
t191 = t204 * pkin(3) + pkin(10);
t291 = -qJ(6) - t191;
t251 = qJD(5) * t291;
t318 = t203 * t251 - t345 + t348;
t78 = t207 * t79;
t317 = t207 * t251 - t243 - t78 + (-qJD(6) + t360) * t203;
t261 = qJD(5) * t335;
t316 = t203 * t261 - t332 + t348;
t263 = -t203 * t57 + t207 * t91;
t315 = -qJD(6) * t203 + t207 * t261 - t243 - t263;
t306 = t135 * t203;
t305 = t135 * t207;
t303 = t157 * t156;
t302 = t179 * t157;
t212 = qJD(1) ^ 2;
t295 = t210 * t212;
t211 = qJD(2) ^ 2;
t294 = t211 * t206;
t293 = t211 * t210;
t289 = (t280 - t359) * pkin(5);
t288 = t206 ^ 2 - t210 ^ 2;
t285 = qJD(2) * t206;
t283 = qJD(3) * t209;
t173 = t206 * t270;
t175 = t210 * t270;
t227 = -t209 * t173 - t205 * t175 - t180 * t283 - t181 * t284;
t83 = -pkin(9) * t139 + t227;
t220 = qJD(3) * t234 + t173 * t205 - t209 * t175;
t84 = pkin(9) * t138 + t220;
t27 = -t87 * qJD(4) + t204 * t84 + t208 * t83;
t197 = pkin(2) * t285;
t129 = pkin(3) * t139 + t197;
t76 = qJD(4) * t135 - t138 * t204 + t208 * t139;
t34 = pkin(4) * t76 - pkin(10) * t75 + t129;
t275 = t203 * t34 + t207 * t27 + t86 * t279;
t271 = -pkin(4) - t338;
t268 = t135 * t279;
t265 = -pkin(2) * t200 - t164;
t258 = pkin(1) * t354;
t65 = t105 * t204 + t104;
t244 = pkin(3) * t282 - t65;
t149 = pkin(2) * t300 - t193 * t208 - pkin(4);
t242 = -t150 * t55 - t322;
t241 = -t191 * t55 - t322;
t23 = -qJ(6) * t108 + t37;
t240 = -t203 * t23 - t324;
t238 = -qJ(6) * t75 - qJD(6) * t135;
t7 = pkin(5) * t43 + t12;
t188 = qJD(2) * t195;
t24 = t55 * pkin(4) - t223 * pkin(10) + t188 + (pkin(10) * t204 + pkin(3)) * qJD(1) * (qJD(2) * t298 + t206 * t283 + t209 * t285 + t210 * t284);
t232 = t207 * t11 + t203 * t24 + t71 * t279 - t280 * t54;
t230 = t203 * t75 + t268;
t229 = -t135 * t280 + t323;
t228 = -t179 * t156 - t254;
t21 = t207 * t24;
t222 = -qJD(5) * t37 - t11 * t203 + t21;
t1 = pkin(5) * t55 + qJ(6) * t42 - qJD(6) * t110 + t222;
t3 = -qJ(6) * t43 - qJD(6) * t108 + t232;
t2 = t3 * t207;
t219 = qJD(5) * t240 - t1 * t203 + t16 * t358 + t23 * t359 + t2;
t28 = qJD(4) * t88 + t204 * t83 - t208 * t84;
t192 = -pkin(4) - t336;
t178 = pkin(10) * t207 + t199;
t177 = t335 * t203;
t168 = t191 * t207 + t199;
t167 = t291 * t203;
t143 = t195 - t339;
t142 = t150 * t207 + t199;
t141 = t292 * t203;
t115 = pkin(3) * t216 + t188;
t112 = -t156 ^ 2 + t157 ^ 2;
t107 = t108 ^ 2;
t100 = -t157 * t200 - t216;
t99 = -t156 * t200 + t132;
t82 = t207 * t86;
t44 = t108 * pkin(5) + qJD(6) + t53;
t38 = -qJ(6) * t306 + t329;
t33 = t207 * t34;
t31 = pkin(5) * t134 - t135 * t199 - t203 * t88 + t82;
t5 = -qJ(6) * t268 + (-qJD(5) * t88 + t238) * t203 + t275;
t4 = pkin(5) * t76 - t203 * t27 + t33 + t238 * t207 + (-t85 + (qJ(6) * t135 - t86) * t203) * qJD(5);
t8 = [0, 0, 0, 0.2e1 * t206 * t266, t288 * t354, t293, -t294, 0, -pkin(7) * t293 + t206 * t258, pkin(7) * t294 + t210 * t258, t132 * t171 + t138 * t157, -t132 * t170 - t138 * t156 + t157 * t139 - t171 * t216, -t138 * t200, -t139 * t200, 0, -t156 * t197 + t179 * t139 + t220 * t200 + (t194 * t225 + (t206 * pkin(2) * t170 + t171 * t194) * qJD(2)) * qJD(1), t194 * t132 - t179 * t138 - t227 * t200 + (-t157 + t344) * t197, t347 * t135 + t236 * t75, -t347 * t134 - t135 * t55 - t236 * t76 + t349 * t75, t75 * t260, -t76 * t260, 0, t115 * t134 - t129 * t349 + t140 * t76 + t144 * t55 - t260 * t28, t115 * t135 + t129 * t236 + t140 * t75 + t144 * t213 - t260 * t27, t110 * t229 - t305 * t42 (-t108 * t207 - t110 * t203) * t75 + (t40 - t207 * t43 + (t108 * t203 - t110 * t207) * qJD(5)) * t135, t110 * t76 - t134 * t42 + t229 * t277 + t305 * t55, -t108 * t76 - t134 * t43 - t230 * t277 - t306 * t55, t134 * t55 + t277 * t76 (-t279 * t88 + t33) * t277 + t82 * t55 + (-t279 * t54 + t21) * t134 + t36 * t76 + t28 * t108 + t87 * t43 + t53 * t268 + ((-qJD(5) * t86 - t27) * t277 - t88 * t55 + (-qJD(5) * t71 - t11) * t134 + t12 * t135 + t53 * t75) * t203 -(-t280 * t88 + t275) * t277 - t329 * t55 - t232 * t134 - t37 * t76 + t28 * t110 - t87 * t42 + t53 * t323 - t264 * t135, -t108 * t5 - t110 * t4 + t31 * t42 - t38 * t43 + t240 * t75 + (-t1 * t207 - t203 * t3 + (t16 * t203 - t207 * t23) * qJD(5)) * t135, t3 * t38 + t23 * t5 + t1 * t31 + t16 * t4 + t7 * (pkin(5) * t306 + t87) + t44 * (pkin(5) * t230 + t28); 0, 0, 0, -t206 * t295, t288 * t212, 0, 0, 0, t212 * pkin(1) * t206, pkin(1) * t295, t303, t112, t99, t100, 0, t156 * t195 + t302 - t253 * t200 + (t205 * t265 - t162) * qJD(3) + t255, t157 * t195 + t290 * t200 + (qJD(3) * t265 + t165) * t209 + t228, -t307, t47, t45, t46, 0, t143 * t349 - t319 * t260 + t231, -t143 * t236 + t351 * t260 + t362, t17, t6, t14, t13, -t312, t149 * t43 + t242 * t203 + t319 * t108 + (-t150 * t279 + t351 * t203 - t70) * t277 + t233, -t149 * t42 + t242 * t207 + t319 * t110 + (t150 * t280 + t346) * t277 + t246, -t108 * t328 - t110 * t327 + t141 * t42 - t142 * t43 + t219, t3 * t142 + t1 * t141 + t7 * (t149 - t338) + t328 * t23 + t327 * t16 + ((qJD(4) * t193 - t114) * t204 + t289 + t356) * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t303, t112, t99, t100, 0, -t200 * t235 + t221 + t302, t200 * t256 + t228 - t350, -t307, t47, t45, t46, 0, t65 * t260 + (-t157 * t349 - t260 * t282) * pkin(3) + t231, t236 * t339 - t304 + t66 * t260 + (-t260 * t326 - t259) * t208 - t262, t17, t6, t14, t13, -t312, t192 * t43 + t241 * t203 + t244 * t108 + (-t191 * t279 + t360 * t203 - t78) * t277 + t233, -t192 * t42 + t241 * t207 + t244 * t110 + (t191 * t280 + t345) * t277 + t246, -t108 * t318 - t110 * t317 + t167 * t42 - t168 * t43 + t219, t3 * t168 + t1 * t167 + t7 * (t271 - t336) + (-t104 + (-t105 + t326) * t204 + t289) * t44 + t318 * t23 + t317 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t307, t47, t45, t46, 0, t260 * t58 + t231, t260 * t57 + t362, t17, t6, t14, t13, -t312, -pkin(4) * t43 - pkin(10) * t320 - t58 * t108 - t263 * t277 - t359 * t53 + t233, pkin(4) * t42 + t332 * t277 - t58 * t110 - t53 * t358 + (t277 * t280 - t52) * pkin(10) + t246, t177 * t42 - t178 * t43 + t2 - t277 * t324 - t315 * t110 - t316 * t108 + (-t23 * t277 - t1) * t203, t3 * t178 + t1 * t177 + t7 * t271 + (pkin(5) * t363 - t58) * t44 + t316 * t23 + t315 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t110 * t108, -t107 + t341, t108 * t277 - t42, t110 * t277 + t357 * t203 - t236 * t279, t55, -t110 * t53 + t277 * t37 + t222, t108 * t53 + t277 * t36 - t232, pkin(5) * t42 - t333 * t108, t333 * t23 + (-t110 * t44 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107 - t341, t108 * t23 + t110 * t16 + t7;];
tauc_reg  = t8;
