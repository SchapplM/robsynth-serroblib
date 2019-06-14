% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PPRRRP2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 20:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PPRRRP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:34:03
% EndTime: 2019-05-04 20:34:13
% DurationCPUTime: 4.03s
% Computational Cost: add. (9973->310), mult. (18290->437), div. (0->0), fcn. (14653->14), ass. (0->222)
t161 = sin(pkin(12));
t162 = sin(pkin(7));
t163 = sin(pkin(6));
t164 = cos(pkin(12));
t165 = cos(pkin(7));
t166 = cos(pkin(6));
t171 = sin(qJ(3));
t174 = cos(qJ(3));
t169 = sin(qJ(5));
t172 = cos(qJ(5));
t170 = sin(qJ(4));
t228 = qJD(3) * qJD(4);
t152 = t170 * t228;
t173 = cos(qJ(4));
t226 = t173 * qJDD(3);
t137 = -t152 + t226;
t130 = -qJDD(5) + t137;
t231 = qJD(3) * t170;
t131 = -t172 * qJD(4) + t169 * t231;
t133 = qJD(4) * t169 + t172 * t231;
t240 = t133 * t131;
t88 = t130 + t240;
t249 = t172 * t88;
t230 = qJD(3) * t173;
t149 = -qJD(5) + t230;
t263 = t149 ^ 2;
t264 = t131 ^ 2;
t267 = -t263 - t264;
t277 = t169 * t267 - t249;
t118 = t133 * t149;
t222 = t173 * t228;
t227 = t170 * qJDD(3);
t136 = t222 + t227;
t218 = t172 * qJDD(4) - t169 * t136;
t191 = qJD(5) * t133 - t218;
t274 = -t118 + t191;
t255 = t169 * t88;
t276 = t172 * t267 + t255;
t290 = t170 * t274 + t173 * t276;
t291 = t170 * t276 - t173 * t274;
t304 = t171 * t290 - t174 * t277;
t311 = t162 * t304 + t165 * t291;
t320 = t166 * t311 + (t161 * (t171 * t277 + t174 * t290) + t164 * (-t162 * t291 + t165 * t304)) * t163;
t186 = t130 - t240;
t256 = t169 * t186;
t129 = t133 ^ 2;
t273 = -t129 - t263;
t60 = -t172 * t273 - t256;
t319 = pkin(3) * t60;
t318 = pkin(4) * t60;
t317 = pkin(10) * t60;
t250 = t172 * t186;
t68 = -t169 * t273 + t250;
t316 = pkin(10) * t68;
t315 = t170 * t68;
t314 = t171 * t60;
t313 = t173 * t68;
t312 = t174 * t60;
t307 = -pkin(3) * t277 + pkin(9) * t290;
t112 = t264 - t263;
t75 = t118 + t191;
t306 = t170 * (t112 * t172 + t256) + t173 * t75;
t193 = -t169 * qJDD(4) - t172 * t136;
t185 = -qJD(5) * t131 - t193;
t241 = t131 * t149;
t269 = t241 + t185;
t258 = t169 * t269;
t272 = t129 - t264;
t305 = t170 * (t172 * t274 + t258) + t173 * t272;
t301 = pkin(4) * t277;
t300 = pkin(10) * t276;
t299 = pkin(10) * t277;
t298 = qJ(6) * t269;
t113 = -t129 + t263;
t295 = t172 * t113 - t255;
t293 = t112 * t169 - t250;
t270 = -t241 + t185;
t289 = t170 * (-t113 * t169 - t249) - t173 * t270;
t271 = t129 + t264;
t288 = pkin(4) * t271;
t213 = -pkin(4) * t173 - pkin(10) * t170;
t134 = t213 * qJD(3);
t262 = qJD(4) ^ 2;
t175 = qJD(3) ^ 2;
t244 = sin(pkin(11));
t245 = cos(pkin(11));
t188 = t244 * g(1) - t245 * g(2);
t232 = -g(3) + qJDD(1);
t182 = t163 * t232 + t166 * t188;
t189 = -t245 * g(1) - t244 * g(2);
t178 = -t161 * t189 + t182 * t164;
t180 = -t163 * t188 + t166 * t232 + qJDD(2);
t278 = t162 * t180 + t165 * t178;
t84 = t182 * t161 + t164 * t189;
t55 = t171 * t278 + t174 * t84;
t53 = -t175 * pkin(3) + qJDD(3) * pkin(9) + t55;
t176 = -t162 * t178 + t165 * t180;
t70 = t173 * t176;
t32 = (qJD(3) * t134 + t53) * t170 - qJDD(4) * pkin(4) - t262 * pkin(10) - t70;
t287 = pkin(5) * t191 - t298 + t32;
t285 = t170 * t271;
t281 = t173 * t271;
t275 = -t169 * t274 + t172 * t269;
t217 = t171 * t84 - t174 * t278;
t268 = t171 * t55 - t174 * t217;
t36 = t170 * t176 + t173 * t53;
t33 = -t262 * pkin(4) + qJDD(4) * pkin(10) + t134 * t230 + t36;
t207 = -t137 + t152;
t208 = t136 + t222;
t52 = -qJDD(3) * pkin(3) - t175 * pkin(9) + t217;
t39 = t207 * pkin(4) - t208 * pkin(10) + t52;
t24 = t169 * t39 + t172 * t33;
t99 = pkin(5) * t131 - qJ(6) * t133;
t220 = -t130 * qJ(6) - t131 * t99 + t24;
t265 = -(t263 + t273) * pkin(5) - qJ(6) * t186 + t220;
t261 = pkin(5) * t172;
t260 = t169 * t32;
t257 = t169 * t270;
t253 = t172 * t32;
t251 = t172 * t270;
t243 = qJ(6) * t172;
t239 = t149 * t169;
t238 = t149 * t172;
t237 = t161 * t163;
t236 = t163 * t164;
t235 = t164 * t165;
t148 = t170 * t175 * t173;
t141 = qJDD(4) + t148;
t234 = t170 * t141;
t142 = qJDD(4) - t148;
t233 = t173 * t142;
t229 = qJD(6) * t149;
t224 = t131 * t238;
t223 = t173 * t240;
t221 = -qJ(6) * t169 - pkin(4);
t23 = t169 * t33 - t172 * t39;
t11 = t169 * t23 + t172 * t24;
t35 = t170 * t53 - t70;
t21 = t170 * t35 + t173 * t36;
t111 = t133 * t239;
t216 = t170 * (t172 * t185 + t111) - t223;
t215 = -t131 * t239 - t172 * t191;
t143 = -0.2e1 * t229;
t214 = t143 + t220;
t14 = -pkin(5) * t263 + t214;
t15 = t130 * pkin(5) - qJ(6) * t263 + t133 * t99 + qJDD(6) + t23;
t212 = -pkin(5) * t15 + qJ(6) * t14;
t211 = -pkin(5) * t270 - qJ(6) * t75;
t29 = (-pkin(5) * t149 - 0.2e1 * qJD(6)) * t133 + t287;
t9 = t14 * t172 + t15 * t169;
t4 = t170 * t29 + t173 * t9;
t8 = t14 * t169 - t15 * t172;
t210 = t171 * t4 - t174 * t8;
t10 = t169 * t24 - t172 * t23;
t6 = t11 * t173 + t170 * t32;
t209 = -t10 * t174 + t171 * t6;
t206 = t171 * t21 - t174 * t52;
t58 = -t172 * t75 + t257;
t42 = t173 * t58 - t285;
t56 = -t169 * t75 - t251;
t205 = t171 * t42 - t174 * t56;
t76 = (-qJD(5) - t149) * t133 + t218;
t59 = t172 * t76 + t257;
t43 = t173 * t59 - t285;
t57 = t169 * t76 - t251;
t204 = t171 * t43 - t174 * t57;
t46 = -t170 * t269 - t313;
t203 = t171 * t46 - t312;
t81 = (qJD(5) - t149) * t131 + t193;
t50 = -t170 * t81 + t313;
t201 = t171 * t50 + t312;
t158 = t173 ^ 2;
t156 = t158 * t175;
t146 = -t156 - t262;
t109 = t146 * t173 - t234;
t138 = -0.2e1 * t152 + t226;
t199 = t109 * t171 + t138 * t174;
t157 = t170 ^ 2;
t154 = t157 * t175;
t145 = -t154 - t262;
t110 = -t145 * t170 - t233;
t135 = 0.2e1 * t222 + t227;
t198 = t110 * t171 - t135 * t174;
t139 = (t157 + t158) * qJDD(3);
t140 = t154 + t156;
t197 = t139 * t171 + t140 * t174;
t196 = -pkin(3) + t213;
t195 = qJDD(3) * t174 - t171 * t175;
t194 = -qJDD(3) * t171 - t174 * t175;
t190 = (t131 * t169 + t133 * t172) * t149;
t187 = t170 * (-t111 + t224) + t173 * t130;
t184 = t170 * (t169 * t191 - t224) + t223;
t183 = 0.2e1 * qJD(6) * t133 - t287;
t181 = -pkin(5) * t88 + qJ(6) * t267 - t15;
t122 = t195 * t162;
t121 = t194 * t162;
t108 = -t142 * t170 + t145 * t173;
t107 = t141 * t173 + t146 * t170;
t98 = t197 * t162;
t72 = -t133 * t238 + t169 * t185;
t65 = t165 * t108 + t198 * t162;
t64 = t165 * t107 + t199 * t162;
t48 = t173 * t81 + t315;
t44 = t173 * t269 - t315;
t41 = t170 * t59 + t281;
t40 = t170 * t58 + t281;
t30 = t268 * t162 + t165 * t176;
t27 = t201 * t162 + t165 * t48;
t25 = t203 * t162 + t165 * t44;
t20 = t170 * t36 - t173 * t35;
t19 = (-t274 + t118) * pkin(5) + t183;
t18 = pkin(5) * t118 + t183 + t298;
t17 = t204 * t162 + t165 * t41;
t16 = t205 * t162 + t165 * t40;
t13 = qJ(6) * t271 + t15;
t12 = (-t263 + t271) * pkin(5) + t214;
t7 = t206 * t162 + t165 * t20;
t5 = t11 * t170 - t173 * t32;
t3 = t170 * t9 - t173 * t29;
t2 = t209 * t162 + t165 * t5;
t1 = t210 * t162 + t165 * t3;
t22 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t232, 0, 0, 0, 0, 0, 0, 0, 0, 0, t166 * t180 + t178 * t236 + t84 * t237, 0, 0, 0, 0, 0, 0, t166 * t122 + (t161 * t194 + t195 * t235) * t163, t166 * t121 + (-t161 * t195 + t194 * t235) * t163, 0, (t171 * t217 + t174 * t55) * t237 + (-t162 * t176 + t268 * t165) * t236 + t166 * t30, 0, 0, 0, 0, 0, 0, t166 * t64 + (t161 * (t109 * t174 - t138 * t171) + t164 * (-t162 * t107 + t199 * t165)) * t163, t166 * t65 + (t161 * (t110 * t174 + t135 * t171) + t164 * (-t162 * t108 + t198 * t165)) * t163, t166 * t98 + (t161 * (t139 * t174 - t140 * t171) + t197 * t235) * t163, t166 * t7 + (t161 * (t171 * t52 + t174 * t21) + t164 * (-t162 * t20 + t206 * t165)) * t163, 0, 0, 0, 0, 0, 0, t320, t166 * t27 + (t161 * (t174 * t50 - t314) + t164 * (-t162 * t48 + t201 * t165)) * t163, t166 * t17 + (t161 * (t171 * t57 + t174 * t43) + t164 * (-t162 * t41 + t204 * t165)) * t163, t166 * t2 + (t161 * (t10 * t171 + t174 * t6) + t164 * (-t162 * t5 + t209 * t165)) * t163, 0, 0, 0, 0, 0, 0, t320, t166 * t16 + (t161 * (t171 * t56 + t174 * t42) + t164 * (-t162 * t40 + t205 * t165)) * t163, t166 * t25 + (t161 * (t174 * t46 + t314) + t164 * (-t162 * t44 + t203 * t165)) * t163, t166 * t1 + (t161 * (t171 * t8 + t174 * t4) + t164 * (-t162 * t3 + t210 * t165)) * t163; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180, 0, 0, 0, 0, 0, 0, t122, t121, 0, t30, 0, 0, 0, 0, 0, 0, t64, t65, t98, t7, 0, 0, 0, 0, 0, 0, t311, t27, t17, t2, 0, 0, 0, 0, 0, 0, t311, t16, t25, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t217, -t55, 0, 0, t208 * t170, t135 * t173 + t138 * t170, t234 + t173 * (-t154 + t262), -t207 * t173, t170 * (t156 - t262) + t233, 0, pkin(3) * t138 + pkin(9) * t109 - t173 * t52, -pkin(3) * t135 + pkin(9) * t110 + t170 * t52, pkin(3) * t140 + pkin(9) * t139 + t21, -pkin(3) * t52 + pkin(9) * t21, t216, -t305, t289, t184, t306, t187, t170 * (t260 - t299) + t173 * (t23 - t301) + t307, t170 * (t253 + t317) + t173 * (t24 + t318) + t319 + pkin(9) * t50, pkin(9) * t43 - t170 * t10 + t196 * t57, pkin(9) * t6 + t196 * t10, t216, t289, t305, t187, -t306, t184, t170 * (-t169 * t19 - t243 * t274 - t299) + t173 * (-t181 - t301) + t307, t170 * (-pkin(10) * t56 - t12 * t169 + t13 * t172) + t173 * (-pkin(4) * t56 - t211) - pkin(3) * t56 + pkin(9) * t42, t170 * (-pkin(5) * t258 + t172 * t18 - t317) + t173 * (0.2e1 * t229 - t265 - t318) - t319 + pkin(9) * t46, t170 * (-pkin(10) * t8 + (pkin(5) * t169 - t243) * t29) + t173 * (-pkin(4) * t8 - t212) - pkin(3) * t8 + pkin(9) * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t148, t154 - t156, t227, t148, t226, qJDD(4), -t35, -t36, 0, 0, t72, t275, t295, t215, t293, t190, -pkin(4) * t274 - t253 + t300, pkin(4) * t81 + t260 + t316, pkin(10) * t59 + t11 + t288, -pkin(4) * t32 + pkin(10) * t11, t72, t295, -t275, t190, -t293, t215, t172 * t19 + t221 * t274 + t300, pkin(10) * t58 + t12 * t172 + t13 * t169 + t288, -t316 + t169 * t18 + (pkin(4) + t261) * t269, pkin(10) * t9 + (t221 - t261) * t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t240, t272, t270, -t240, -t75, -t130, -t23, -t24, 0, 0, t240, t270, -t272, -t130, t75, -t240, t181, t211, t143 + t265, t212; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88, t270, t273, t15;];
tauJ_reg  = t22;
