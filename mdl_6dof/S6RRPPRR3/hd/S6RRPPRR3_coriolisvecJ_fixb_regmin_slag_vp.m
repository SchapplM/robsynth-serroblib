% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d5,d6,theta3,theta4]';
% 
% Output:
% tauc_reg [6x30]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 09:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPRR3_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR3_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPPRR3_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:58:38
% EndTime: 2019-03-09 08:58:57
% DurationCPUTime: 7.14s
% Computational Cost: add. (9511->436), mult. (29068->620), div. (0->0), fcn. (24316->12), ass. (0->230)
t215 = sin(qJ(6));
t218 = cos(qJ(6));
t219 = cos(qJ(5));
t214 = cos(pkin(6));
t274 = t214 * qJD(1);
t202 = qJD(2) + t274;
t213 = cos(pkin(12));
t211 = sin(pkin(11));
t217 = sin(qJ(2));
t220 = cos(qJ(2));
t296 = cos(pkin(11));
t231 = t211 * t220 + t217 * t296;
t212 = sin(pkin(6));
t280 = qJD(1) * t212;
t183 = t231 * t280;
t210 = sin(pkin(12));
t289 = t210 * t183;
t323 = t213 * t202 - t289;
t131 = t219 * t323;
t143 = -t213 * t183 - t210 * t202;
t216 = sin(qJ(5));
t94 = t216 * t143 + t131;
t87 = qJD(6) - t94;
t328 = t87 ^ 2;
t259 = t296 * t220;
t250 = t212 * t259;
t196 = qJD(1) * t250;
t279 = qJD(2) * t212;
t266 = t217 * t279;
t251 = qJD(1) * t266;
t174 = qJD(2) * t196 - t211 * t251;
t193 = t219 * t210 + t216 * t213;
t319 = t174 * t193;
t324 = -t219 * t143 + t216 * t323;
t59 = qJD(5) * t324 + t319;
t329 = -t215 * t328 + t218 * t59;
t267 = t217 * t280;
t179 = t211 * t267 - t196;
t175 = qJD(5) + t179;
t327 = t94 * t175;
t69 = -t218 * t175 + t215 * t324;
t326 = t324 * t69;
t71 = t215 * t175 + t218 * t324;
t325 = t324 * t71;
t284 = t219 * t213;
t192 = t216 * t210 - t284;
t122 = t192 * t179;
t188 = t192 * qJD(5);
t283 = -t188 - t122;
t282 = t175 * t193;
t320 = pkin(2) * t266;
t204 = t211 * pkin(2) + qJ(4);
t304 = pkin(9) + t204;
t190 = t304 * t210;
t191 = t304 * t213;
t236 = -t219 * t190 - t216 * t191;
t308 = pkin(9) * t213;
t310 = pkin(1) * t217;
t271 = t214 * t310;
t287 = t212 * t220;
t305 = pkin(8) + qJ(3);
t177 = t305 * t287 + t271;
t164 = t177 * qJD(1);
t153 = t211 * t164;
t309 = pkin(1) * t220;
t270 = t214 * t309;
t200 = qJD(1) * t270;
t265 = t305 * t217;
t253 = t212 * t265;
t163 = -qJD(1) * t253 + t200;
t114 = t163 * t296 - t153;
t123 = pkin(2) * t267 + t183 * pkin(3) + t179 * qJ(4);
t63 = -t210 * t114 + t213 * t123;
t46 = t183 * pkin(4) + t179 * t308 + t63;
t295 = t179 * t210;
t64 = t213 * t114 + t210 * t123;
t53 = pkin(9) * t295 + t64;
t318 = qJD(4) * t192 - qJD(5) * t236 + t216 * t46 + t219 * t53;
t139 = -t216 * t190 + t219 * t191;
t317 = -qJD(4) * t193 - qJD(5) * t139 + t216 * t53 - t219 * t46;
t316 = -qJD(5) + t175;
t249 = (-pkin(2) * t220 - pkin(1)) * t212;
t235 = qJD(1) * t249;
t187 = qJD(3) + t235;
t110 = t179 * pkin(3) - t183 * qJ(4) + t187;
t147 = t202 * pkin(2) + t163;
t260 = t296 * t164;
t103 = t211 * t147 + t260;
t98 = t202 * qJ(4) + t103;
t54 = t213 * t110 - t210 * t98;
t35 = t179 * pkin(4) + pkin(9) * t143 + t54;
t55 = t210 * t110 + t213 * t98;
t41 = pkin(9) * t323 + t55;
t14 = t216 * t35 + t219 * t41;
t186 = t231 * t212;
t181 = qJD(2) * t186;
t173 = qJD(1) * t181;
t286 = t213 * t174;
t198 = qJD(2) * t200;
t225 = (-qJD(2) * t265 + qJD(3) * t220) * t212;
t134 = qJD(1) * t225 + t198;
t288 = t212 * t217;
t149 = -t177 * qJD(2) - qJD(3) * t288;
t135 = t149 * qJD(1);
t82 = t296 * t134 + t211 * t135;
t74 = t202 * qJD(4) + t82;
t197 = pkin(2) * t251;
t83 = t173 * pkin(3) - t174 * qJ(4) - t183 * qJD(4) + t197;
t39 = -t210 * t74 + t213 * t83;
t28 = t173 * pkin(4) - pkin(9) * t286 + t39;
t290 = t210 * t174;
t40 = t210 * t83 + t213 * t74;
t31 = -pkin(9) * t290 + t40;
t262 = t216 * t31 - t219 * t28;
t4 = -t173 * pkin(5) + qJD(5) * t14 + t262;
t315 = (pkin(5) * t324 + t87 * pkin(10)) * t87 + t4;
t276 = qJD(6) * t215;
t85 = t218 * t122 + t215 * t183;
t228 = t218 * t188 + t193 * t276 + t85;
t292 = t193 * t218;
t314 = -t228 * t87 + t59 * t292;
t206 = -pkin(2) * t296 - pkin(3);
t195 = -t213 * pkin(4) + t206;
t132 = t192 * pkin(5) - t193 * pkin(10) + t195;
t113 = t211 * t163 + t260;
t80 = -pkin(4) * t295 + t113;
t13 = -t216 * t41 + t219 * t35;
t9 = -t175 * pkin(5) - t13;
t313 = -t9 * qJD(6) * t193 - t132 * t59 + (-t282 * pkin(5) + t283 * pkin(10) + qJD(6) * t139 + t80) * t87;
t58 = qJD(5) * t131 + t174 * t284 + (qJD(5) * t143 - t290) * t216;
t22 = qJD(6) * t71 - t218 * t173 + t215 * t58;
t312 = -t192 * t22 - t282 * t69;
t311 = -t193 * t173 - t175 * t283;
t81 = t211 * t134 - t296 * t135;
t68 = pkin(4) * t290 + t81;
t12 = t59 * pkin(5) - t58 * pkin(10) + t68;
t10 = t175 * pkin(10) + t14;
t102 = t147 * t296 - t153;
t97 = -t202 * pkin(3) + qJD(4) - t102;
t67 = -pkin(4) * t323 + t97;
t23 = -pkin(5) * t94 - pkin(10) * t324 + t67;
t245 = t215 * t10 - t218 * t23;
t277 = qJD(5) * t219;
t278 = qJD(5) * t216;
t233 = t216 * t28 + t219 * t31 + t35 * t277 - t278 * t41;
t3 = t173 * pkin(10) + t233;
t1 = -qJD(6) * t245 + t215 * t12 + t218 * t3;
t176 = t179 ^ 2;
t307 = t69 * t87;
t306 = t71 * t87;
t157 = t186 * t213 + t214 * t210;
t185 = t211 * t288 - t250;
t162 = (pkin(2) + t309) * t214 - t253;
t120 = t211 * t162 + t296 * t177;
t111 = t214 * qJ(4) + t120;
t126 = t185 * pkin(3) - t186 * qJ(4) + t249;
t65 = -t210 * t111 + t213 * t126;
t45 = t185 * pkin(4) - t157 * pkin(9) + t65;
t156 = t186 * t210 - t214 * t213;
t66 = t213 * t111 + t210 * t126;
t51 = -t156 * pkin(9) + t66;
t239 = t216 * t45 + t219 * t51;
t201 = qJD(2) * t270;
t148 = t201 + t225;
t101 = t296 * t148 + t211 * t149;
t89 = t214 * qJD(4) + t101;
t182 = (-t211 * t217 + t259) * t279;
t99 = t181 * pkin(3) - t182 * qJ(4) - t186 * qJD(4) + t320;
t48 = t210 * t99 + t213 * t89;
t302 = t183 * t94;
t275 = qJD(6) * t218;
t21 = t215 * t173 + t175 * t275 + t218 * t58 - t276 * t324;
t300 = t21 * t215;
t299 = t215 * t59;
t298 = t324 * t183;
t297 = t183 * pkin(5) - t317;
t294 = t182 * t210;
t207 = t212 ^ 2;
t221 = qJD(1) ^ 2;
t291 = t207 * t221;
t281 = t217 ^ 2 - t220 ^ 2;
t273 = qJD(2) - t202;
t272 = t207 * t310;
t268 = t220 * t291;
t264 = qJD(1) * qJD(2) * t207;
t47 = -t210 * t89 + t213 * t99;
t257 = t218 * t87;
t100 = t211 * t148 - t296 * t149;
t255 = t202 + t274;
t254 = t21 * t192 + t282 * t71;
t252 = t220 * t264;
t247 = -t192 * t173 - t282 * t175;
t72 = pkin(4) * t294 + t100;
t8 = t218 * t10 + t215 * t23;
t16 = t185 * pkin(10) + t239;
t104 = t219 * t156 + t216 * t157;
t105 = -t216 * t156 + t219 * t157;
t119 = t162 * t296 - t211 * t177;
t112 = -t214 * pkin(3) - t119;
t77 = t156 * pkin(4) + t112;
t29 = t104 * pkin(5) - t105 * pkin(10) + t77;
t244 = t218 * t16 + t215 * t29;
t243 = -t215 * t16 + t218 * t29;
t242 = -t210 * t54 + t213 * t55;
t34 = t181 * pkin(4) - t182 * t308 + t47;
t38 = -pkin(9) * t294 + t48;
t241 = -t216 * t38 + t219 * t34;
t240 = -t216 * t51 + t219 * t45;
t237 = t112 * t174 + t97 * t182;
t76 = t218 * t105 + t185 * t215;
t75 = t215 * t105 - t185 * t218;
t232 = t216 * t34 + t219 * t38 + t45 * t277 - t278 * t51;
t230 = -pkin(8) * t287 - t271;
t229 = -pkin(8) * t251 + t198;
t227 = t230 * t202;
t226 = -pkin(10) * t59 + (t13 + t9) * t87;
t224 = -t204 * t173 + t206 * t174 + (-qJD(4) + t97) * t179;
t2 = -qJD(6) * t8 + t218 * t12 - t215 * t3;
t84 = t215 * t122 - t218 * t183;
t223 = (t215 * t188 - t193 * t275 + t84) * t87 - t193 * t299;
t222 = -t139 * t59 - t9 * t188 + t4 * t193 + (t183 * pkin(10) - qJD(6) * t132 + t318) * t87;
t62 = qJD(5) * t105 + t182 * t193;
t61 = -qJD(5) * t104 - t182 * t192;
t25 = qJD(6) * t76 - t181 * t218 + t215 * t61;
t24 = -qJD(6) * t75 + t181 * t215 + t218 * t61;
t17 = t62 * pkin(5) - t61 * pkin(10) + t72;
t15 = -t185 * pkin(5) - t240;
t6 = -t181 * pkin(5) + qJD(5) * t239 - t241;
t5 = t181 * pkin(10) + t232;
t7 = [0, 0, 0, 0.2e1 * t217 * t252, -0.2e1 * t281 * t264, t255 * t220 * t279, -t255 * t266, 0 (t227 + (t214 * t230 - 0.2e1 * t272) * qJD(1)) * qJD(2), -0.2e1 * pkin(1) * t252 - (-pkin(8) * t266 + t201) * t202 - t229 * t214, t100 * t183 - t101 * t179 - t102 * t182 - t103 * t181 - t119 * t174 - t120 * t173 - t82 * t185 + t81 * t186, -t102 * t100 + t103 * t101 - t81 * t119 + t82 * t120 + (t187 + t235) * t320, -t100 * t323 + t81 * t156 + t65 * t173 + t47 * t179 + t54 * t181 + t39 * t185 + t210 * t237, -t100 * t143 + t81 * t157 - t66 * t173 - t48 * t179 - t55 * t181 - t40 * t185 + t213 * t237, t48 * t323 - t40 * t156 + t47 * t143 - t39 * t157 + (-t55 * t210 - t54 * t213) * t182 + (-t66 * t210 - t65 * t213) * t174, t97 * t100 + t81 * t112 + t39 * t65 + t40 * t66 + t54 * t47 + t55 * t48, t58 * t105 + t324 * t61, -t58 * t104 - t105 * t59 - t324 * t62 + t61 * t94, t105 * t173 + t61 * t175 + t181 * t324 + t58 * t185, -t104 * t173 - t62 * t175 + t181 * t94 - t59 * t185, t173 * t185 + t175 * t181, t241 * t175 + t240 * t173 - t262 * t185 + t13 * t181 - t72 * t94 + t77 * t59 + t68 * t104 + t67 * t62 + (-t14 * t185 - t175 * t239) * qJD(5), t68 * t105 - t14 * t181 - t239 * t173 - t232 * t175 - t233 * t185 + t324 * t72 + t77 * t58 + t67 * t61, t21 * t76 + t71 * t24, -t21 * t75 - t22 * t76 - t24 * t69 - t25 * t71, t21 * t104 + t24 * t87 + t76 * t59 + t71 * t62, -t22 * t104 - t25 * t87 - t75 * t59 - t69 * t62, t59 * t104 + t87 * t62 (-qJD(6) * t244 + t218 * t17 - t215 * t5) * t87 + t243 * t59 + t2 * t104 - t245 * t62 + t6 * t69 + t15 * t22 + t4 * t75 + t9 * t25 -(qJD(6) * t243 + t215 * t17 + t218 * t5) * t87 - t244 * t59 - t1 * t104 - t8 * t62 + t6 * t71 + t15 * t21 + t4 * t76 + t9 * t24; 0, 0, 0, -t217 * t268, t281 * t291, t273 * t220 * t280, -t273 * t267, 0, t221 * t272 + (qJD(2) * t230 - t227) * qJD(1), pkin(1) * t268 + (-pkin(8) * t267 + t200) * t202 - t229 (t103 - t113) * t183 + (-t102 + t114) * t179 + (-t173 * t211 - t174 * t296) * pkin(2), t102 * t113 - t103 * t114 + (-t187 * t267 + t211 * t82 - t296 * t81) * pkin(2), t113 * t323 - t63 * t179 - t54 * t183 + t210 * t224 - t81 * t213, t113 * t143 + t64 * t179 + t55 * t183 + t81 * t210 + t213 * t224, -t63 * t143 + t64 * t289 + (qJD(4) * t323 - t54 * t179 - t64 * t202 + t40) * t213 + (-qJD(4) * t143 - t55 * t179 - t39) * t210, -t97 * t113 + t81 * t206 - t54 * t63 - t55 * t64 + (-t39 * t210 + t40 * t213) * t204 + t242 * qJD(4), t58 * t193 + t283 * t324, -t58 * t192 - t193 * t59 - t282 * t324 + t283 * t94, -t298 - t311, t247 - t302, -t175 * t183, -t13 * t183 + t173 * t236 + t317 * t175 + t68 * t192 + t195 * t59 + t282 * t67 + t80 * t94, -t139 * t173 + t14 * t183 + t318 * t175 + t68 * t193 + t195 * t58 + t283 * t67 - t324 * t80, t21 * t292 - t228 * t71, t85 * t69 + t71 * t84 - (-t215 * t71 - t218 * t69) * t188 + (-t300 - t218 * t22 + (t215 * t69 - t218 * t71) * qJD(6)) * t193, t254 + t314, t223 + t312, t59 * t192 + t282 * t87, t2 * t192 + t222 * t215 - t313 * t218 - t22 * t236 - t245 * t282 + t297 * t69 - t9 * t84, -t1 * t192 - t21 * t236 + t313 * t215 + t222 * t218 - t282 * t8 + t297 * t71 - t9 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t183 ^ 2 - t176, t102 * t183 + t103 * t179 + t197, t213 * t173 - t210 * t176 + t183 * t323, t143 * t183 - t210 * t173 - t213 * t176 (-t143 * t210 + t213 * t323) * t179 + (-t210 ^ 2 - t213 ^ 2) * t174, t179 * t242 - t97 * t183 + t40 * t210 + t39 * t213, 0, 0, 0, 0, 0, t247 + t302, -t298 + t311, 0, 0, 0, 0, 0, t223 - t312, t254 - t314; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t143 * t179 + t290, t179 * t323 + t286, -t143 ^ 2 - t323 ^ 2, -t54 * t143 - t323 * t55 + t81, 0, 0, 0, 0, 0, t175 * t324 + t59, t58 + t327, 0, 0, 0, 0, 0, -t326 + t329, -t218 * t328 - t299 - t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t324 * t94, t324 ^ 2 - t94 ^ 2, t58 - t327, t316 * t324 - t319, t173, t316 * t14 - t324 * t67 - t262, t13 * t175 - t67 * t94 - t233, t257 * t71 + t300 (t21 - t307) * t218 + (-t22 - t306) * t215, t257 * t87 + t299 - t325, t326 + t329, -t87 * t324, -pkin(5) * t22 - t14 * t69 + t226 * t215 - t315 * t218 + t245 * t324, -pkin(5) * t21 - t14 * t71 + t315 * t215 + t226 * t218 + t324 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t69, -t69 ^ 2 + t71 ^ 2, t21 + t307, -t22 + t306, t59, -t9 * t71 + t8 * t87 + t2, -t245 * t87 + t9 * t69 - t1;];
tauc_reg  = t7;
