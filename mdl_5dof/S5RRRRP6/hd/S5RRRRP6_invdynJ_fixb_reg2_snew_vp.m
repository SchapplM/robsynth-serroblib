% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRP6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP6_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:54:35
% EndTime: 2019-12-31 21:54:44
% DurationCPUTime: 3.14s
% Computational Cost: add. (14016->357), mult. (28663->461), div. (0->0), fcn. (20103->8), ass. (0->235)
t206 = sin(qJ(3));
t207 = sin(qJ(2));
t210 = cos(qJ(3));
t211 = cos(qJ(2));
t181 = (t211 * t206 + t207 * t210) * qJD(1);
t202 = qJD(2) + qJD(3);
t205 = sin(qJ(4));
t209 = cos(qJ(4));
t166 = t205 * t181 - t209 * t202;
t168 = t209 * t181 + t205 * t202;
t139 = t168 * t166;
t197 = t207 * qJDD(1);
t250 = qJD(1) * qJD(2);
t239 = t211 * t250;
t186 = t197 + t239;
t198 = t211 * qJDD(1);
t240 = t207 * t250;
t187 = t198 - t240;
t232 = t206 * t186 - t210 * t187;
t142 = -t181 * qJD(3) - t232;
t141 = qJDD(4) - t142;
t301 = -t139 + t141;
t307 = pkin(4) * t301;
t255 = qJD(1) * t207;
t179 = -t210 * t211 * qJD(1) + t206 * t255;
t229 = t210 * t186 + t206 * t187;
t143 = -t179 * qJD(3) + t229;
t249 = qJDD(2) + qJDD(3);
t120 = -t166 * qJD(4) + t209 * t143 + t205 * t249;
t176 = qJD(4) + t179;
t150 = t176 * t166;
t106 = t150 + t120;
t306 = qJ(5) * t106;
t267 = t205 * t301;
t158 = t181 * t179;
t299 = -t158 + t249;
t305 = t206 * t299;
t261 = t209 * t301;
t304 = t210 * t299;
t164 = t168 ^ 2;
t175 = t176 ^ 2;
t134 = -t164 - t175;
t163 = t166 ^ 2;
t234 = t205 * t143 - t209 * t249;
t119 = -t168 * qJD(4) - t234;
t146 = t176 * pkin(4) - t168 * qJ(5);
t174 = t202 * t179;
t130 = t143 - t174;
t204 = t211 ^ 2;
t213 = qJD(1) ^ 2;
t208 = sin(qJ(1));
t287 = cos(qJ(1));
t238 = t208 * g(1) - t287 * g(2);
t225 = qJDD(1) * pkin(1) + t238;
t226 = qJD(2) * pkin(2) - pkin(7) * t255;
t145 = t187 * pkin(2) + (pkin(7) * t204 + pkin(6)) * t213 - t226 * t255 + t225;
t75 = -t130 * pkin(8) + (t202 * t181 - t142) * pkin(3) - t145;
t227 = t287 * g(1) + t208 * g(2);
t273 = qJDD(1) * pkin(6);
t216 = -t213 * pkin(1) - t227 + t273;
t171 = -t207 * g(3) + t211 * t216;
t200 = t204 * t213;
t138 = -pkin(2) * t200 + t187 * pkin(7) - qJD(2) * t226 + t171;
t260 = t210 * t138;
t215 = t207 * t216;
t263 = t207 * t213;
t284 = t186 * pkin(7);
t298 = qJDD(2) * pkin(2) - t215 + (pkin(2) * t263 + pkin(7) * t250 - g(3)) * t211 - t284;
t110 = t298 * t206 + t260;
t156 = t179 * pkin(3) - t181 * pkin(8);
t297 = t202 ^ 2;
t85 = -t297 * pkin(3) + t249 * pkin(8) - t179 * t156 + t110;
t45 = t205 * t75 + t209 * t85;
t222 = t119 * qJ(5) - 0.2e1 * qJD(5) * t166 - t176 * t146 + t45;
t303 = -t222 + (t134 + t163) * pkin(4);
t44 = t205 * t85 - t209 * t75;
t24 = t205 * t44 + t209 * t45;
t300 = -t150 + t120;
t103 = (qJD(4) - t176) * t168 + t234;
t177 = t179 ^ 2;
t178 = t181 ^ 2;
t125 = -t175 - t163;
t78 = t205 * t125 + t261;
t296 = pkin(3) * t78;
t112 = t139 + t141;
t268 = t205 * t112;
t88 = t209 * t134 - t268;
t295 = pkin(3) * t88;
t253 = qJD(5) * t168;
t160 = -0.2e1 * t253;
t221 = -t306 - t44 + t307;
t28 = t160 + t221;
t294 = pkin(4) * t28;
t122 = -t163 - t164;
t68 = -t103 * t209 + t205 * t106;
t51 = -t210 * t122 + t206 * t68;
t293 = pkin(7) * t51;
t102 = (qJD(4) + t176) * t168 + t234;
t79 = t209 * t125 - t267;
t56 = -t210 * t102 + t206 * t79;
t292 = pkin(7) * t56;
t262 = t209 * t112;
t89 = -t205 * t134 - t262;
t60 = t206 * t89 - t210 * t300;
t291 = pkin(7) * t60;
t66 = -t103 * t205 - t209 * t106;
t290 = pkin(8) * t66;
t289 = pkin(8) * t78;
t288 = pkin(8) * t88;
t286 = pkin(3) * t206;
t285 = pkin(4) * t106;
t52 = t206 * t122 + t210 * t68;
t283 = pkin(6) * (-t207 * t51 + t211 * t52) - pkin(1) * t66;
t57 = t206 * t102 + t210 * t79;
t282 = pkin(6) * (-t207 * t56 + t211 * t57) - pkin(1) * t78;
t61 = t206 * t300 + t210 * t89;
t281 = pkin(6) * (-t207 * t60 + t211 * t61) - pkin(1) * t88;
t109 = t206 * t138 - t210 * t298;
t84 = -t249 * pkin(3) - t297 * pkin(8) + t181 * t156 + t109;
t280 = -pkin(3) * t84 + pkin(8) * t24;
t279 = t205 * t28;
t80 = t205 * t84;
t70 = -t210 * t109 + t206 * t110;
t278 = t207 * t70;
t277 = t209 * t28;
t81 = t209 * t84;
t276 = -pkin(3) * t122 + pkin(8) * t68;
t275 = -pkin(3) * t102 + pkin(8) * t79;
t274 = -pkin(3) * t300 + pkin(8) * t89;
t272 = t176 * t205;
t271 = t176 * t209;
t270 = t202 * t206;
t269 = t202 * t210;
t266 = t206 * t145;
t154 = t158 + t249;
t265 = t206 * t154;
t192 = t211 * t263;
t189 = qJDD(2) + t192;
t264 = t207 * t189;
t259 = t210 * t145;
t258 = t210 * t154;
t257 = t211 * (qJDD(2) - t192);
t252 = qJD(3) + t202;
t248 = t80 + t274;
t247 = -t81 + t275;
t246 = t206 * t139;
t245 = t210 * t139;
t244 = -pkin(3) * t210 - pkin(2);
t243 = -pkin(2) * t78 + pkin(7) * t57;
t242 = -pkin(2) * t88 + pkin(7) * t61;
t19 = -qJ(5) * t103 + (-t122 - t163) * pkin(4) + t222;
t161 = 0.2e1 * t253;
t21 = t161 - t221 + t306;
t237 = t209 * t19 + t205 * t21 + t276;
t47 = -t119 * pkin(4) - t163 * qJ(5) + t168 * t146 + qJDD(5) + t84;
t38 = -qJ(5) * t134 + t47;
t72 = -pkin(4) * t300 - qJ(5) * t112;
t236 = t205 * t38 + t209 * t72 + t274;
t235 = t276 + t24;
t71 = t206 * t109 + t210 * t110;
t170 = t211 * g(3) + t215;
t233 = t207 * t170 + t211 * t171;
t230 = t205 * t45 - t209 * t44;
t32 = -t163 * pkin(4) + t222;
t11 = -pkin(4) * t47 + qJ(5) * t32;
t8 = t209 * t32 - t279;
t224 = -pkin(3) * t47 + pkin(8) * t8 - qJ(5) * t279 + t209 * t11;
t34 = -pkin(4) * t102 + qJ(5) * t125 - t47;
t223 = -qJ(5) * t267 + t209 * t34 + t275;
t220 = (-qJD(3) + t202) * t181 - t232;
t217 = t221 + t307;
t212 = qJD(2) ^ 2;
t203 = t207 ^ 2;
t199 = t203 * t213;
t188 = t198 - 0.2e1 * t240;
t185 = t197 + 0.2e1 * t239;
t182 = t213 * pkin(6) + t225;
t173 = -t178 + t297;
t172 = t177 - t297;
t169 = -t178 - t297;
t157 = t178 - t177;
t152 = -t297 - t177;
t148 = -t164 + t175;
t147 = t163 - t175;
t144 = -t177 - t178;
t136 = t164 - t163;
t133 = -t206 * t169 - t258;
t132 = t210 * t169 - t265;
t131 = t143 + t174;
t129 = -t252 * t179 + t229;
t126 = t252 * t181 + t232;
t124 = t210 * t152 - t305;
t123 = t206 * t152 + t304;
t115 = (-t166 * t209 + t168 * t205) * t176;
t114 = (-t166 * t205 - t168 * t209) * t176;
t99 = t209 * t120 - t168 * t272;
t98 = t205 * t120 + t168 * t271;
t97 = -t205 * t119 + t166 * t271;
t96 = t209 * t119 + t166 * t272;
t95 = t206 * t131 + t210 * t220;
t94 = -t210 * t131 + t206 * t220;
t93 = t209 * t147 - t268;
t92 = -t205 * t148 + t261;
t91 = t205 * t147 + t262;
t90 = t209 * t148 + t267;
t67 = -t209 * t102 - t205 * t300;
t65 = -t205 * t102 + t209 * t300;
t62 = t207 * (t210 * t115 + t206 * t141) + t211 * (t206 * t115 - t210 * t141);
t59 = pkin(2) * t60;
t55 = pkin(2) * t56;
t53 = t81 - t288;
t50 = pkin(2) * t51;
t49 = pkin(7) * t52;
t48 = t80 - t289;
t41 = t207 * (t210 * t99 + t246) + t211 * (t206 * t99 - t245);
t40 = t207 * (t210 * t97 - t246) + t211 * (t206 * t97 + t245);
t39 = -pkin(3) * t66 + t285;
t36 = t45 - t295;
t35 = t44 - t296;
t31 = t207 * (-t206 * t103 + t210 * t93) + t211 * (t210 * t103 + t206 * t93);
t30 = t207 * (t206 * t106 + t210 * t92) + t211 * (-t210 * t106 + t206 * t92);
t26 = t207 * (t206 * t136 + t210 * t67) + t211 * (-t210 * t136 + t206 * t67);
t17 = -t295 - t303;
t16 = -t205 * t72 + t209 * t38 - t288;
t15 = -qJ(5) * t261 - t205 * t34 - t289;
t14 = t161 - t217 - t296;
t12 = t206 * t24 - t210 * t84;
t9 = -t230 - t290;
t7 = t205 * t32 + t277;
t5 = t206 * t47 + t210 * t8;
t4 = t206 * t8 - t210 * t47;
t3 = -t205 * t19 + t209 * t21 - t290;
t2 = -pkin(3) * t7 - t294;
t1 = -pkin(8) * t7 - qJ(5) * t277 - t205 * t11;
t6 = [0, 0, 0, 0, 0, qJDD(1), t238, t227, 0, 0, (t186 + t239) * t207, t211 * t185 + t207 * t188, t264 + t211 * (-t199 + t212), (t187 - t240) * t211, t207 * (t200 - t212) + t257, 0, t211 * t182 + pkin(1) * t188 + pkin(6) * (t211 * (-t200 - t212) - t264), -t207 * t182 - pkin(1) * t185 + pkin(6) * (-t257 - t207 * (-t199 - t212)), pkin(1) * (t199 + t200) + (t203 + t204) * t273 + t233, pkin(1) * t182 + pkin(6) * t233, t207 * (t210 * t143 - t181 * t270) + t211 * (t206 * t143 + t181 * t269), t207 * (-t210 * t126 - t206 * t130) + t211 * (-t206 * t126 + t210 * t130), t207 * (-t206 * t173 + t304) + t211 * (t210 * t173 + t305), t207 * (-t206 * t142 + t179 * t269) + t211 * (t210 * t142 + t179 * t270), t207 * (t210 * t172 - t265) + t211 * (t206 * t172 + t258), (t207 * (-t179 * t210 + t181 * t206) + t211 * (-t179 * t206 - t181 * t210)) * t202, t207 * (-pkin(7) * t123 - t266) + t211 * (-pkin(2) * t126 + pkin(7) * t124 + t259) - pkin(1) * t126 + pkin(6) * (-t207 * t123 + t211 * t124), t207 * (-pkin(7) * t132 - t259) + t211 * (-pkin(2) * t129 + pkin(7) * t133 - t266) - pkin(1) * t129 + pkin(6) * (-t207 * t132 + t211 * t133), t207 * (-pkin(7) * t94 - t70) + t211 * (-pkin(2) * t144 + pkin(7) * t95 + t71) - pkin(1) * t144 + pkin(6) * (-t207 * t94 + t211 * t95), -pkin(7) * t278 + t211 * (pkin(2) * t145 + pkin(7) * t71) + pkin(1) * t145 + pkin(6) * (t211 * t71 - t278), t41, t26, t30, t40, t31, t62, t207 * (-t206 * t35 + t210 * t48 - t292) + t211 * (t206 * t48 + t210 * t35 + t243) + t282, t207 * (-t206 * t36 + t210 * t53 - t291) + t211 * (t206 * t53 + t210 * t36 + t242) + t281, t207 * (t210 * t9 + t66 * t286 - t293) + t211 * (t206 * t9 + t244 * t66 + t49) + t283, (t207 * (-pkin(8) * t210 + t286) + t211 * (-pkin(8) * t206 + t244) - pkin(1)) * t230 + (pkin(6) + pkin(7)) * (-t207 * t12 + t211 * (t206 * t84 + t210 * t24)), t41, t26, t30, t40, t31, t62, t207 * (-t206 * t14 + t210 * t15 - t292) + t211 * (t210 * t14 + t206 * t15 + t243) + t282, t207 * (t210 * t16 - t206 * t17 - t291) + t211 * (t206 * t16 + t210 * t17 + t242) + t281, t207 * (-t206 * t39 + t210 * t3 - t293) + t211 * (-pkin(2) * t66 + t206 * t3 + t210 * t39 + t49) + t283, t207 * (-pkin(7) * t4 + t210 * t1 - t206 * t2) + t211 * (-pkin(2) * t7 + pkin(7) * t5 + t206 * t1 + t210 * t2) - pkin(1) * t7 + pkin(6) * (-t207 * t4 + t211 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t192, t199 - t200, t197, t192, t198, qJDD(2), -t170, -t171, 0, 0, t158, t157, t131, -t158, t220, t249, pkin(2) * t123 - t109, -t260 - t206 * (pkin(7) * t239 - t170 - t284) + (-t189 * t206 + t132) * pkin(2), pkin(2) * t94, pkin(2) * t70, t98, t65, t90, t96, t91, t114, t55 + t247, t59 + t248, t50 + t235, pkin(2) * t12 + t280, t98, t65, t90, t96, t91, t114, t223 + t55, t59 + t236, t50 + t237, pkin(2) * t4 + t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t157, t131, -t158, t220, t249, -t109, -t110, 0, 0, t98, t65, t90, t96, t91, t114, t247, t248, t235, t280, t98, t65, t90, t96, t91, t114, t223, t236, t237, t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t136, t106, -t139, -t103, t141, -t44, -t45, 0, 0, t139, t136, t106, -t139, -t103, t141, t160 + t217, t303, -t285, t294; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t102, t300, t122, t47;];
tauJ_reg = t6;
