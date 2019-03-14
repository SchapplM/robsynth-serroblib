% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPRRPP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRRPP6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:48:28
% EndTime: 2019-03-09 04:48:40
% DurationCPUTime: 4.32s
% Computational Cost: add. (5910->457), mult. (13063->601), div. (0->0), fcn. (8354->6), ass. (0->226)
t172 = sin(qJ(3));
t254 = qJD(1) * t172;
t158 = qJD(4) + t254;
t171 = sin(qJ(4));
t173 = cos(qJ(4));
t243 = t173 * qJD(3);
t174 = cos(qJ(3));
t253 = qJD(1) * t174;
t134 = t171 * t253 - t243;
t231 = t173 * t253;
t252 = qJD(3) * t171;
t136 = t231 + t252;
t170 = sin(pkin(9));
t275 = cos(pkin(9));
t79 = t275 * t134 + t136 * t170;
t284 = t79 * t158;
t244 = qJD(4) * t174;
t226 = t171 * t244;
t229 = t172 * t243;
t186 = t226 + t229;
t91 = qJD(1) * t186 - qJD(4) * t243;
t251 = qJD(3) * t172;
t230 = t171 * t251;
t248 = qJD(4) * t136;
t92 = -qJD(1) * t230 + t248;
t53 = -t170 * t92 - t275 * t91;
t28 = t53 + t284;
t188 = -t170 * t134 + t275 * t136;
t310 = t79 * t188;
t220 = t275 * t171;
t131 = t170 * t173 + t220;
t116 = t131 * qJD(4);
t117 = t131 * qJD(1);
t278 = t172 * t117 + t116;
t219 = t275 * t173;
t265 = t170 * t171;
t187 = t219 - t265;
t246 = qJD(4) * t172;
t250 = qJD(3) * t174;
t277 = t131 * t246 - t187 * t250 + t117;
t247 = qJD(4) * t171;
t118 = qJD(4) * t219 - t170 * t247;
t204 = qJD(1) * t219;
t232 = t171 * t254;
t276 = -t170 * t232 + t172 * t204 + t118;
t225 = t173 * t244;
t309 = t225 - t230;
t300 = t188 ^ 2;
t240 = 0.2e1 * qJD(1);
t296 = -qJ(5) - pkin(8);
t221 = qJD(4) * t296;
t242 = t173 * qJD(5);
t110 = t171 * t221 + t242;
t185 = -t171 * qJD(5) + t173 * t221;
t263 = t172 * t173;
t203 = pkin(3) * t174 + pkin(8) * t172;
t138 = t203 * qJD(1);
t175 = -pkin(1) - pkin(7);
t153 = t175 * qJD(1) + qJD(2);
t264 = t171 * t174;
t235 = t153 * t264;
t86 = t173 * t138 - t235;
t62 = (pkin(4) * t174 + qJ(5) * t263) * qJD(1) + t86;
t261 = t173 * t174;
t87 = t171 * t138 + t153 * t261;
t71 = qJ(5) * t232 + t87;
t293 = (-t185 + t62) * t275 + (t110 - t71) * t170;
t269 = t153 * t174;
t290 = qJD(3) * pkin(3);
t124 = -t269 - t290;
t84 = pkin(4) * t134 + qJD(5) + t124;
t27 = pkin(5) * t79 - qJ(6) * t188 + t84;
t308 = t27 * t188;
t129 = qJD(3) * t203 + qJD(2);
t109 = t129 * qJD(1);
t140 = pkin(3) * t172 - pkin(8) * t174 + qJ(2);
t120 = t140 * qJD(1);
t139 = t172 * t153;
t123 = qJD(3) * pkin(8) + t139;
t228 = t174 * t243;
t245 = qJD(4) * t173;
t38 = t171 * t109 + t120 * t245 - t123 * t247 + t153 * t228;
t74 = t173 * t120 - t123 * t171;
t307 = -t158 * t74 + t38;
t75 = t120 * t171 + t123 * t173;
t184 = -qJD(4) * t75 + t173 * t109;
t39 = -qJD(3) * t235 + t184;
t306 = -t158 * t75 - t39;
t207 = -t139 + (t232 + t247) * pkin(4);
t262 = t172 * t175;
t95 = t171 * t140 + t173 * t262;
t105 = t131 * t174;
t52 = -t170 * t91 + t275 * t92;
t104 = t131 * t172;
t65 = qJD(3) * t104 - t118 * t174;
t305 = -t65 * t158 + t172 * t52 + (qJD(1) * t105 + t79) * t250;
t304 = t278 * t158 + (-qJD(3) * t187 - t79) * t253;
t303 = t131 * t52 - t187 * t53 + t188 * t278 + t276 * t79;
t106 = t187 * t172;
t302 = (t106 * t253 - t172 * t188) * qJD(3) - t277 * t158 + t174 * t53;
t58 = -qJ(5) * t134 + t75;
t54 = t275 * t58;
t57 = -qJ(5) * t136 + t74;
t24 = t170 * t57 + t54;
t301 = t24 * t188;
t298 = t79 ^ 2;
t18 = qJ(5) * t91 - qJD(5) * t136 + (pkin(4) * qJD(1) - t153 * t171) * t250 + t184;
t23 = -qJ(5) * t92 - qJD(5) * t134 + t38;
t3 = -t170 * t23 + t275 * t18;
t4 = t170 * t18 + t275 * t23;
t32 = t170 * t62 + t275 * t71;
t29 = qJ(6) * t253 + t32;
t69 = t275 * t110 + t170 * t185;
t295 = -t29 + t69;
t294 = -pkin(5) * t253 - t293;
t292 = t32 - t69;
t182 = -t95 * qJD(4) + t173 * t129;
t222 = -t171 * t175 + pkin(4);
t36 = qJ(5) * t229 + (qJ(5) * t247 + t222 * qJD(3) - t242) * t174 + t182;
t249 = qJD(3) * t175;
t227 = t174 * t249;
t233 = t171 * t129 + t140 * t245 + t173 * t227;
t41 = -qJ(5) * t225 + (-qJD(5) * t174 + (qJ(5) * qJD(3) - qJD(4) * t175) * t172) * t171 + t233;
t11 = t170 * t36 + t275 * t41;
t291 = -pkin(5) * t278 + qJ(6) * t276 + qJD(6) * t131 - t207;
t50 = pkin(4) * t158 + t57;
t22 = t170 * t50 + t54;
t128 = t173 * t140;
t77 = -qJ(5) * t261 + t222 * t172 + t128;
t85 = -qJ(5) * t264 + t95;
t47 = t170 * t77 + t275 * t85;
t287 = t158 * t188;
t286 = t170 * t58;
t283 = t91 * t171;
t282 = t91 * t172;
t281 = t92 * t172;
t280 = t92 * t173;
t218 = qJD(3) * t275;
t206 = t174 * t218;
t279 = qJD(1) * t265 - t118 * t172 - t170 * t228 - t171 * t206 - t204;
t148 = t296 * t173;
t89 = -t148 * t170 - t296 * t220;
t274 = qJD(3) * t89;
t90 = -t275 * t148 + t296 * t265;
t273 = qJD(3) * t90;
t272 = t134 * t158;
t271 = t136 * t134;
t270 = t136 * t158;
t268 = t158 * t171;
t267 = t158 * t172;
t266 = t158 * t173;
t176 = qJD(3) ^ 2;
t260 = t176 * t172;
t259 = t176 * t174;
t177 = qJD(1) ^ 2;
t258 = t177 * qJ(2);
t25 = t275 * t57 - t286;
t257 = qJD(6) - t25;
t169 = t174 ^ 2;
t256 = t172 ^ 2 - t169;
t255 = -t176 - t177;
t241 = qJD(1) * qJD(3);
t163 = t174 * t241;
t239 = qJ(6) * t163 + t4;
t237 = qJD(2) * t240;
t236 = t171 * t262;
t234 = t174 * t177 * t172;
t165 = -pkin(4) * t173 - pkin(3);
t224 = t158 * t253;
t70 = pkin(4) * t92 + t153 * t251;
t215 = -t124 + t269;
t132 = pkin(4) * t264 - t174 * t175;
t214 = t158 + t254;
t213 = t134 + t243;
t212 = -t136 + t252;
t211 = qJD(1) + t246;
t210 = -t90 * t52 + t53 * t89 - t69 * t79;
t209 = t24 * t158 + t3;
t208 = t172 * t163;
t202 = -t298 - t300;
t201 = -t298 + t300;
t200 = t105 * t52 - t65 * t79;
t197 = t171 * t75 + t173 * t74;
t196 = t171 * t74 - t173 * t75;
t195 = t215 * qJD(3);
t194 = qJD(1) * t169 - t267;
t193 = t52 + t287;
t192 = -t52 + t287;
t191 = -pkin(8) * t250 + t124 * t172;
t10 = -t170 * t41 + t275 * t36;
t21 = t275 * t50 - t286;
t46 = -t170 * t85 + t275 * t77;
t189 = -t187 * t52 + t278 * t79;
t2 = -pkin(5) * t163 - t3;
t93 = pkin(4) * t309 + t172 * t249;
t183 = -t53 + t284;
t107 = -t170 * t264 + t174 * t219;
t67 = t174 * t116 - t170 * t230 + t218 * t263;
t181 = t105 * t53 + t107 * t52 - t188 * t65 - t67 * t79;
t5 = pkin(5) * t52 - qJ(6) * t53 - qJD(6) * t188 + t70;
t180 = t104 * t53 - t106 * t52 - t188 * t279 + t277 * t79;
t179 = t79 * t251 + (-t104 * t241 - t52) * t174 + t279 * t158;
t178 = -qJD(4) * t197 - t171 * t39 + t173 * t38;
t166 = qJ(2) * t237;
t162 = -t275 * pkin(4) - pkin(5);
t160 = pkin(4) * t170 + qJ(6);
t100 = t214 * t250;
t94 = t128 - t236;
t76 = -pkin(5) * t187 - qJ(6) * t131 + t165;
t60 = -t171 * t227 + t182;
t59 = -qJD(4) * t236 + t233;
t56 = pkin(5) * t105 - qJ(6) * t107 + t132;
t42 = -t172 * pkin(5) - t46;
t40 = qJ(6) * t172 + t47;
t33 = pkin(4) * t136 + pkin(5) * t188 + qJ(6) * t79;
t26 = t276 * t158 + (qJD(3) * t131 - t188) * t253;
t15 = qJ(6) * t158 + t22;
t14 = -t158 * pkin(5) + qJD(6) - t21;
t13 = -pkin(5) * t65 + qJ(6) * t67 - qJD(6) * t107 + t93;
t12 = t107 * t53 - t188 * t67;
t9 = -pkin(5) * t250 - t10;
t8 = -t67 * t158 + t53 * t172 + (qJD(1) * t107 + t188) * t250;
t7 = qJ(6) * t250 + qJD(6) * t172 + t11;
t6 = t53 * t131 + t188 * t276;
t1 = qJD(6) * t158 + t239;
t16 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t237, t166, -0.2e1 * t208, 0.2e1 * t256 * t241, -t260, 0.2e1 * t208, -t259, 0, -t175 * t260 + (qJ(2) * t250 + qJD(2) * t172) * t240, -t175 * t259 + (-qJ(2) * t251 + qJD(2) * t174) * t240, 0, t166, -t136 * t186 - t261 * t91 (t134 * t173 + t136 * t171) * t251 + (t283 - t280 + (t134 * t171 - t136 * t173) * qJD(4)) * t174, -t158 * t226 - t282 + (t136 * t174 + t173 * t194) * qJD(3), t134 * t309 + t92 * t264, -t158 * t225 - t281 + (-t134 * t174 - t171 * t194) * qJD(3), t100, t158 * t60 + t172 * t39 + (t124 * t245 - t175 * t92) * t174 + ((qJD(1) * t94 + t74) * t174 + (t134 * t175 + t215 * t171) * t172) * qJD(3), -t158 * t59 - t172 * t38 + (-t124 * t247 + t175 * t91) * t174 + ((-qJD(1) * t95 - t75) * t174 + (t136 * t175 + t173 * t215) * t172) * qJD(3), -t134 * t59 - t136 * t60 + t91 * t94 - t92 * t95 + t197 * t251 + (qJD(4) * t196 - t171 * t38 - t173 * t39) * t174, -t195 * t262 + t38 * t95 + t39 * t94 + t59 * t75 + t60 * t74, t12, -t181, t8, t200, -t305, t100, t10 * t158 + t105 * t70 + t132 * t52 + t172 * t3 - t65 * t84 + t79 * t93 + (qJD(1) * t46 + t21) * t250, t107 * t70 - t11 * t158 + t132 * t53 - t172 * t4 - t67 * t84 + t188 * t93 + (-qJD(1) * t47 - t22) * t250, -t10 * t188 - t105 * t4 - t107 * t3 - t11 * t79 + t21 * t67 + t22 * t65 - t46 * t53 - t47 * t52, t10 * t21 + t11 * t22 + t132 * t70 + t3 * t46 + t4 * t47 + t84 * t93, t12, t8, t181, t100, t305, t200, t105 * t5 + t13 * t79 - t158 * t9 - t172 * t2 - t27 * t65 + t52 * t56 + (-qJD(1) * t42 - t14) * t250, -t1 * t105 + t107 * t2 - t14 * t67 + t15 * t65 + t188 * t9 - t40 * t52 + t42 * t53 - t7 * t79, t1 * t172 - t107 * t5 - t13 * t188 + t158 * t7 + t27 * t67 - t53 * t56 + (qJD(1) * t40 + t15) * t250, t1 * t40 + t13 * t27 + t14 * t9 + t15 * t7 + t2 * t42 + t5 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t177, -t258, 0, 0, 0, 0, 0, 0, t255 * t172, t255 * t174, 0, -t258, 0, 0, 0, 0, 0, 0, -t174 * t92 - t211 * t266 + (t134 * t172 - t214 * t264) * qJD(3), t174 * t91 + t211 * t268 + (-t158 * t261 + (t136 - t231) * t172) * qJD(3) (-t134 * t250 + t136 * t211 - t281) * t173 + (t134 * t211 + t136 * t250 - t282) * t171, -t196 * t250 - t197 * qJD(1) + (-t195 + t178) * t172, 0, 0, 0, 0, 0, 0, t179, -t302, t180, -t104 * t3 + t106 * t4 - t174 * t70 + t279 * t21 - t277 * t22 + t84 * t251, 0, 0, 0, 0, 0, 0, t179, t180, t302, t1 * t106 + t104 * t2 - t14 * t279 - t15 * t277 - t174 * t5 + t251 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t234, -t256 * t177, 0, -t234, 0, 0, -t174 * t258, t172 * t258, 0, 0, t136 * t266 - t283 (-t91 - t272) * t173 + (-t92 - t270) * t171, t158 * t245 + (t158 * t263 + t174 * t212) * qJD(1), t134 * t268 - t280, -t158 * t247 + (-t171 * t267 + t174 * t213) * qJD(1), -t224, -pkin(3) * t92 - t158 * t86 - t213 * t139 + (-pkin(8) * t266 + t124 * t171) * qJD(4) + (t171 * t191 - t174 * t74) * qJD(1), pkin(3) * t91 + t158 * t87 + t212 * t139 + (pkin(8) * t268 + t124 * t173) * qJD(4) + (t173 * t191 + t174 * t75) * qJD(1), t134 * t87 + t136 * t86 + ((-t92 + t248) * pkin(8) + t307) * t173 + ((qJD(4) * t134 - t91) * pkin(8) + t306) * t171, -t74 * t86 - t75 * t87 + (-t124 - t290) * t139 + t178 * pkin(8), t6, -t303, t26, t189, -t304, -t224, -t187 * t70 + t165 * t52 + t278 * t84 + t207 * t79 - t293 * t158 + (-t21 - t274) * t253, t131 * t70 + t165 * t53 + t276 * t84 + t207 * t188 + t292 * t158 + (t22 - t273) * t253, -t131 * t3 + t187 * t4 + t188 * t293 - t276 * t21 - t278 * t22 + t32 * t79 + t210, t165 * t70 + t207 * t84 - t293 * t21 - t292 * t22 - t3 * t89 + t4 * t90, t6, t26, t303, -t224, t304, t189, -t187 * t5 + t52 * t76 - t291 * t79 + t278 * t27 + t294 * t158 + (t14 - t274) * t253, t1 * t187 + t131 * t2 + t14 * t276 - t15 * t278 - t188 * t294 + t29 * t79 + t210, -t131 * t5 - t53 * t76 + t291 * t188 - t276 * t27 + t295 * t158 + (-t15 + t273) * t253, t1 * t90 - t14 * t294 + t15 * t295 + t2 * t89 - t27 * t291 + t5 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t271, -t134 ^ 2 + t136 ^ 2, t272 - t91, -t271, t270 - t92, t163, -t124 * t136 - t306, t124 * t134 - t307, 0, 0, t310, t201, t28, -t310, t192, t163, -t84 * t188 + (qJD(1) * t206 - t136 * t79) * pkin(4) + t209, t158 * t25 + t79 * t84 + (-t136 * t188 - t163 * t170) * pkin(4) - t4, t22 * t188 - t301 + (-t170 * t52 - t275 * t53) * pkin(4) + (-t21 + t25) * t79, t21 * t24 - t22 * t25 + (-t136 * t84 + t170 * t4 + t275 * t3) * pkin(4), t310, t28, -t201, t163, -t192, -t310, -t308 - t33 * t79 + (pkin(5) - t162) * t163 + t209, t15 * t188 - t160 * t52 + t162 * t53 - t301 + (t14 - t257) * t79, t160 * t163 - t27 * t79 + t33 * t188 + (0.2e1 * qJD(6) - t25) * t158 + t239, t1 * t160 - t14 * t24 + t15 * t257 + t162 * t2 - t27 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, -t183, t202, t188 * t21 + t22 * t79 + t70, 0, 0, 0, 0, 0, 0, t193, t202, t183, -t14 * t188 + t15 * t79 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t163 + t310, t28, -t158 ^ 2 - t300, -t15 * t158 + t2 + t308;];
tauc_reg  = t16;