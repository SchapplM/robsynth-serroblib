% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6PPRPRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PPRPRR1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_invdynJ_fixb_reg2_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:43:43
% EndTime: 2019-03-08 18:43:52
% DurationCPUTime: 4.99s
% Computational Cost: add. (6641->462), mult. (17447->674), div. (0->0), fcn. (16682->16), ass. (0->230)
t172 = sin(qJ(5));
t175 = cos(qJ(5));
t214 = pkin(5) * t172 - pkin(10) * t175;
t136 = t214 * qJD(5);
t215 = pkin(5) * t175 + pkin(10) * t172;
t205 = -pkin(4) - t215;
t161 = sin(pkin(13));
t166 = cos(pkin(13));
t285 = qJDD(3) * pkin(3);
t173 = sin(qJ(3));
t162 = sin(pkin(12));
t165 = sin(pkin(6));
t279 = t162 * t165;
t225 = qJDD(1) * t279;
t170 = cos(pkin(6));
t143 = t170 * qJDD(1) + qJDD(2);
t167 = cos(pkin(12));
t169 = cos(pkin(7));
t176 = cos(qJ(3));
t268 = t169 * t176;
t235 = t167 * t268;
t219 = t165 * t235;
t164 = sin(pkin(7));
t274 = t164 * t176;
t262 = qJDD(1) * t219 + t143 * t274;
t146 = qJD(1) * t170 + qJD(2);
t275 = t164 * t173;
t270 = t167 * t169;
t202 = t162 * t176 + t173 * t270;
t299 = t202 * t165;
t84 = qJD(1) * t299 + t146 * t275;
t61 = -t84 * qJD(3) - t173 * t225 + t262;
t59 = t61 + t285;
t278 = t162 * t173;
t203 = t235 - t278;
t247 = qJD(1) * qJD(3);
t60 = (qJD(3) * t176 * t146 + t173 * t143) * t164 + (qJDD(1) * t202 + t203 * t247) * t165;
t293 = t161 * t60 - t166 * t59;
t12 = qJD(3) * t136 + qJDD(3) * t205 + t293;
t171 = sin(qJ(6));
t174 = cos(qJ(6));
t273 = t165 * t167;
t233 = qJD(1) * t273;
t106 = t146 * t169 - t164 * t233 + qJD(4);
t79 = t166 * t84;
t218 = t169 * t233;
t241 = t165 * t278;
t83 = -qJD(1) * t241 + t146 * t274 + t176 * t218;
t81 = qJD(3) * pkin(3) + t83;
t52 = t161 * t81 + t79;
t50 = qJD(3) * pkin(9) + t52;
t35 = t106 * t172 + t175 * t50;
t33 = qJD(5) * pkin(10) + t35;
t78 = t161 * t84;
t51 = t166 * t81 - t78;
t36 = qJD(3) * t205 - t51;
t209 = t171 * t33 - t174 * t36;
t240 = t164 * t273;
t107 = -qJDD(1) * t240 + t169 * t143;
t103 = qJDD(4) + t107;
t23 = t161 * t59 + t166 * t60;
t19 = qJDD(3) * pkin(9) + t23;
t253 = qJD(5) * t175;
t244 = -t172 * t103 - t106 * t253 - t175 * t19;
t254 = qJD(5) * t172;
t7 = -t50 * t254 - t244;
t5 = qJDD(5) * pkin(10) + t7;
t1 = -t209 * qJD(6) + t171 * t12 + t174 * t5;
t257 = qJD(3) * t175;
t147 = -qJD(6) + t257;
t303 = -t209 * t147 + t1;
t55 = t161 * t83 + t79;
t302 = t136 - t55;
t91 = t170 * t275 + t299;
t250 = qJD(6) * t172;
t301 = qJD(3) * t250 - qJDD(5);
t245 = t172 * qJDD(3);
t86 = t171 * ((qJD(6) + t257) * qJD(5) + t245) + t301 * t174;
t10 = t171 * t36 + t174 * t33;
t2 = -qJD(6) * t10 + t174 * t12 - t171 * t5;
t300 = -t10 * t147 + t2;
t117 = t169 * t170 - t240;
t110 = t117 * t175;
t237 = t170 * t274;
t90 = t203 * t165 + t237;
t65 = t161 * t90 + t166 * t91;
t298 = -t172 * t65 + t110;
t100 = t175 * t103;
t8 = -qJD(5) * t35 - t172 * t19 + t100;
t156 = t175 * qJDD(3);
t246 = qJD(3) * qJD(5);
t127 = t172 * t246 + qJDD(6) - t156;
t248 = t174 * qJD(5);
t228 = t175 * t248;
t193 = -t171 * t250 + t228;
t266 = t172 * t174;
t297 = t127 * t266 - t147 * t193;
t206 = t161 * t176 + t166 * t173;
t113 = t206 * t164;
t115 = t206 * t169;
t264 = t176 * t166;
t128 = t161 * t173 - t264;
t69 = t113 * t170 + (t115 * t167 - t128 * t162) * t165;
t295 = pkin(3) * t166;
t150 = pkin(3) * t161 + pkin(9);
t231 = t150 * t254;
t265 = t174 * t175;
t56 = t166 * t83 - t78;
t122 = t205 - t295;
t267 = t171 * t175;
t94 = t122 * t174 - t150 * t267;
t292 = t94 * qJD(6) + t171 * t302 - t174 * t231 - t56 * t265;
t95 = t122 * t171 + t150 * t265;
t291 = -t95 * qJD(6) + t171 * t231 + t174 * t302 + t56 * t267;
t289 = t172 * t50;
t258 = qJD(3) * t172;
t132 = t171 * t258 - t248;
t287 = -t132 * t228 - t86 * t266;
t286 = qJD(3) * t56;
t283 = t117 * t172;
t282 = t132 * t147;
t255 = qJD(5) * t171;
t134 = t174 * t258 + t255;
t281 = t134 * t132;
t280 = t134 * t147;
t163 = sin(pkin(11));
t277 = t163 * t165;
t276 = t163 * t170;
t168 = cos(pkin(11));
t272 = t165 * t168;
t271 = t165 * t169;
t269 = t168 * t170;
t120 = -t162 * t168 - t167 * t276;
t239 = t165 * t274;
t243 = pkin(3) * t268;
t263 = t163 * pkin(3) * t239 + t120 * t243;
t261 = (t219 + t237) * pkin(3);
t159 = t172 ^ 2;
t160 = t175 ^ 2;
t260 = t159 - t160;
t256 = qJD(5) * t132;
t252 = qJD(6) * t132;
t251 = qJD(6) * t171;
t249 = qJD(6) * t174;
t178 = qJD(3) ^ 2;
t234 = t172 * t178 * t175;
t232 = t147 * t255;
t230 = t134 * t253;
t229 = t172 * t249;
t227 = t175 * t246;
t85 = -qJD(6) * t248 + (-t227 - t245) * t174 + t301 * t171;
t223 = t134 * t254 + t175 * t85;
t222 = -t85 + t252;
t220 = t168 * t239;
t217 = t134 * t229;
t216 = t172 * t227;
t119 = t162 * t269 + t163 * t167;
t121 = -t162 * t276 + t167 * t168;
t213 = g(1) * t121 + g(2) * t119;
t212 = t10 * t174 + t171 * t209;
t211 = -t10 * t171 + t174 * t209;
t40 = t175 * t65 + t283;
t64 = t161 * t91 - t166 * t90;
t21 = t171 * t64 + t174 * t40;
t20 = -t171 * t40 + t174 * t64;
t34 = t106 * t175 - t289;
t208 = t34 * t172 - t35 * t175;
t112 = t161 * t275 - t164 * t264;
t97 = t113 * t175 + t169 * t172;
t74 = t112 * t174 - t171 * t97;
t75 = t112 * t171 + t174 * t97;
t96 = t113 * t172 - t175 * t169;
t201 = -t127 * t171 + t147 * t249;
t118 = -t162 * t163 + t167 * t269;
t41 = t113 * t272 - t115 * t118 + t119 * t128;
t46 = t113 * t277 + t115 * t120 - t121 * t128;
t92 = -t118 * t164 - t168 * t271;
t93 = -t120 * t164 + t163 * t271;
t200 = -g(1) * (-t172 * t46 + t175 * t93) - g(2) * (t172 * t41 + t175 * t92) - g(3) * (-t172 * t69 + t110);
t29 = t172 * t92 - t175 * t41;
t31 = t172 * t93 + t175 * t46;
t48 = t175 * t69 + t283;
t199 = -g(1) * t31 - g(2) * t29 - g(3) * t48;
t198 = -g(1) * t46 + g(2) * t41 - g(3) * t69;
t114 = t128 * t169;
t42 = t112 * t272 - t114 * t118 - t119 * t206;
t45 = -t112 * t277 - t114 * t120 - t121 * t206;
t68 = -t112 * t170 + (-t114 * t167 - t162 * t206) * t165;
t197 = g(1) * t45 + g(2) * t42 + g(3) * t68;
t196 = -g(1) * t93 - g(2) * t92 - g(3) * t117;
t195 = -g(1) * t277 + g(2) * t272 - g(3) * t170;
t6 = -qJDD(5) * pkin(5) - t8;
t194 = t200 - t6;
t32 = -qJD(5) * pkin(5) - t34;
t190 = -pkin(10) * t127 - t147 * t32;
t189 = -pkin(3) * t121 * t173 + t45 * pkin(4) + pkin(9) * t46 + t263;
t188 = t213 - t225;
t187 = qJD(3) * t55 - t197;
t186 = -pkin(3) * t241 + t68 * pkin(4) + pkin(9) * t69 + t261;
t151 = -pkin(4) - t295;
t49 = -qJD(3) * pkin(4) - t51;
t185 = -qJDD(5) * t150 + (qJD(3) * t151 + t49 + t56) * qJD(5);
t184 = pkin(10) * qJD(6) * t147 + t194;
t183 = t211 * qJD(6) + t1 * t174 - t2 * t171;
t182 = -t8 * t172 + t7 * t175 + (-t172 * t35 - t175 * t34) * qJD(5);
t104 = t118 * t243;
t180 = -pkin(9) * t41 + t104 + t42 * pkin(4) + (-t119 * t173 - t220) * pkin(3);
t177 = qJD(5) ^ 2;
t18 = -qJDD(3) * pkin(4) + t293;
t179 = qJDD(3) * t151 + t150 * t177 + t18 - t187;
t139 = qJDD(5) * t175 - t172 * t177;
t138 = qJDD(5) * t172 + t175 * t177;
t135 = t214 * qJD(3);
t109 = t128 * t164 * qJD(3);
t108 = qJD(3) * t113;
t88 = t91 * qJD(3);
t87 = t90 * qJD(3);
t73 = t97 * qJD(5) - t109 * t172;
t72 = -t96 * qJD(5) - t109 * t175;
t63 = -t161 * t88 + t166 * t87;
t62 = t161 * t87 + t166 * t88;
t27 = t135 * t171 + t174 * t34;
t26 = t135 * t174 - t171 * t34;
t25 = -t75 * qJD(6) + t108 * t174 - t171 * t72;
t24 = t74 * qJD(6) + t108 * t171 + t174 * t72;
t14 = qJD(5) * t40 + t172 * t63;
t13 = qJD(5) * t298 + t175 * t63;
t4 = -t21 * qJD(6) - t13 * t171 + t174 * t62;
t3 = t20 * qJD(6) + t13 * t174 + t171 * t62;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t143 * t170 - g(3) + (t162 ^ 2 + t167 ^ 2) * t165 ^ 2 * qJDD(1), 0, 0, 0, 0, 0, 0, -qJD(3) * t88 + qJDD(3) * t90, -qJD(3) * t87 - qJDD(3) * t91, 0, t107 * t117 + t60 * t91 + t61 * t90 - t83 * t88 + t84 * t87 - g(3), 0, 0, 0, 0, 0, 0, -qJD(3) * t62 - qJDD(3) * t64, -qJD(3) * t63 - qJDD(3) * t65, 0, t103 * t117 + t23 * t65 + t293 * t64 - t51 * t62 + t52 * t63 - g(3), 0, 0, 0, 0, 0, 0, -t64 * t156 - qJD(5) * t14 + qJDD(5) * t298 + (-t175 * t62 + t64 * t254) * qJD(3), t64 * t245 - qJD(5) * t13 - qJDD(5) * t40 + (t172 * t62 + t253 * t64) * qJD(3) (-t172 * t298 + t175 * t40) * qJDD(3) + (t13 * t175 + t14 * t172 + (-t172 * t40 - t175 * t298) * qJD(5)) * qJD(3), t13 * t35 - t14 * t34 + t18 * t64 + t298 * t8 + t40 * t7 + t49 * t62 - g(3), 0, 0, 0, 0, 0, 0, t127 * t20 + t132 * t14 - t147 * t4 - t298 * t86, -t127 * t21 + t134 * t14 + t147 * t3 + t298 * t85, -t132 * t3 - t134 * t4 + t20 * t85 - t21 * t86, t1 * t21 + t10 * t3 + t14 * t32 + t2 * t20 - t209 * t4 - t298 * t6 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t195 + t143, 0, 0, 0, 0, 0, 0 (qJDD(3) * t176 - t173 * t178) * t164 (-qJDD(3) * t173 - t176 * t178) * t164, 0, t107 * t169 + (t173 * t60 + t176 * t61 + (-t173 * t83 + t176 * t84) * qJD(3)) * t164 + t195, 0, 0, 0, 0, 0, 0, -qJD(3) * t108 - qJDD(3) * t112, qJD(3) * t109 - qJDD(3) * t113, 0, t103 * t169 - t108 * t51 - t109 * t52 + t112 * t293 + t113 * t23 + t195, 0, 0, 0, 0, 0, 0, -t112 * t156 - qJD(5) * t73 - qJDD(5) * t96 + (-t108 * t175 + t112 * t254) * qJD(3), t112 * t245 - qJD(5) * t72 - qJDD(5) * t97 + (t108 * t172 + t112 * t253) * qJD(3) (t172 * t96 + t175 * t97) * qJDD(3) + (t172 * t73 + t175 * t72 + (-t172 * t97 + t175 * t96) * qJD(5)) * qJD(3), t108 * t49 + t112 * t18 - t34 * t73 + t35 * t72 + t7 * t97 - t8 * t96 + t195, 0, 0, 0, 0, 0, 0, t127 * t74 + t132 * t73 - t147 * t25 + t86 * t96, -t127 * t75 + t134 * t73 + t147 * t24 - t85 * t96, -t132 * t24 - t134 * t25 + t74 * t85 - t75 * t86, t1 * t75 + t10 * t24 + t2 * t74 - t209 * t25 + t32 * t73 + t6 * t96 + t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -g(3) * t90 + t188 * t173 + (-g(1) * (t120 * t169 + t164 * t277) - g(2) * (t118 * t169 - t164 * t272)) * t176 + t262, g(3) * t91 + t83 * qJD(3) + ((-t146 * t164 - t218) * qJD(3) + t188) * t176 + (-t164 * t143 + (g(1) * t120 + g(2) * t118) * t169 + (t162 * t247 - qJDD(1) * t270 + (g(1) * t163 - g(2) * t168) * t164) * t165) * t173, 0, 0, 0, 0, 0, 0, 0, qJDD(3), t166 * t285 + t187 - t293, -t161 * t285 - t198 - t23 + t286, 0, -t52 * t56 + t51 * t55 - g(1) * t263 - g(2) * t104 - g(3) * t261 + (g(2) * t220 + t23 * t161 - t293 * t166 + (g(3) * t279 + t213) * t173) * pkin(3), qJDD(3) * t159 + 0.2e1 * t216, 0.2e1 * t172 * t156 - 0.2e1 * t260 * t246, t138, qJDD(3) * t160 - 0.2e1 * t216, t139, 0, t172 * t185 - t175 * t179, t172 * t179 + t175 * t185, t182 + t198 + (qJDD(3) * t150 - t286) * (t159 + t160) -g(1) * t189 - g(2) * t180 - g(3) * t186 + t150 * t182 + t18 * t151 + t208 * t56 - t49 * t55, t134 * t193 - t266 * t85, -t217 + (-t230 + (t85 + t252) * t172) * t171 + t287, t223 + t297, t171 * t172 * t86 + (t171 * t253 + t229) * t132 (t86 + t232) * t175 + (t201 - t256) * t172, -t127 * t175 - t147 * t254, t94 * t127 - t291 * t147 + t198 * t171 + (-t2 + (t132 * t150 + t171 * t32) * qJD(5) - t197 * t174) * t175 + (-qJD(5) * t209 - t56 * t132 + t150 * t86 + t6 * t171 + t249 * t32) * t172, -t95 * t127 + t292 * t147 + t198 * t174 + (t1 + (t134 * t150 + t174 * t32) * qJD(5) + t197 * t171) * t175 + (-t10 * qJD(5) - t56 * t134 - t150 * t85 + t6 * t174 - t251 * t32) * t172, t85 * t94 - t86 * t95 - t291 * t134 - t292 * t132 + t211 * t253 + (-qJD(6) * t212 - t1 * t171 - t174 * t2 - t197) * t172, t1 * t95 + t2 * t94 - t32 * t172 * t56 - g(1) * (t215 * t45 + t189) - g(2) * (t215 * t42 + t180) - g(3) * (t215 * t68 + t186) - t291 * t209 + (t6 * t172 + t253 * t32) * t150 + t292 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t103 + t196, 0, 0, 0, 0, 0, 0, t139, -t138, 0, -qJD(5) * t208 + t172 * t7 + t175 * t8 + t196, 0, 0, 0, 0, 0, 0 (-t86 + t232) * t175 + (t201 + t256) * t172, t223 - t297, t217 + (t172 * t222 + t230) * t171 + t287 (qJD(5) * t212 - t6) * t175 + (qJD(5) * t32 + t183) * t172 + t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, t260 * t178, t245, t234, t156, qJDD(5), t100 + (-qJD(3) * t49 - t19) * t172 + t200, -t49 * t257 + (t34 + t289) * qJD(5) - t199 + t244, 0, 0, -t171 * t85 - t174 * t280 (-t85 + t282) * t174 + (-t86 + t280) * t171 (-t134 * t172 + t147 * t265) * qJD(3) - t201, -t171 * t282 - t174 * t86, t147 * t251 + t127 * t174 + (t132 * t172 - t147 * t267) * qJD(3), t147 * t258, -pkin(5) * t86 - t132 * t35 + t147 * t26 + t171 * t190 + t174 * t184 + t209 * t258, pkin(5) * t85 + t10 * t258 - t134 * t35 - t147 * t27 - t171 * t184 + t174 * t190, t132 * t27 + t134 * t26 + ((qJD(6) * t134 - t86) * pkin(10) + t303) * t174 + (pkin(10) * t222 - t300) * t171 + t199, -t10 * t27 + t209 * t26 - t32 * t35 + t194 * pkin(5) + (t183 + t199) * pkin(10); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t281, -t132 ^ 2 + t134 ^ 2, -t85 - t282, -t281, -t280 - t86, t127, -t32 * t134 - g(1) * (-t171 * t31 - t174 * t45) - g(2) * (-t171 * t29 - t174 * t42) - g(3) * (-t171 * t48 - t174 * t68) + t300, t32 * t132 - g(1) * (t171 * t45 - t174 * t31) - g(2) * (t171 * t42 - t174 * t29) - g(3) * (t171 * t68 - t174 * t48) - t303, 0, 0;];
tau_reg  = t9;
