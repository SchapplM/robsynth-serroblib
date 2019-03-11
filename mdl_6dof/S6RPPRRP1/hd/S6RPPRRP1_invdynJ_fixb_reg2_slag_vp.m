% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRRP1_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP1_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP1_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP1_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:58:56
% EndTime: 2019-03-09 01:58:59
% DurationCPUTime: 3.51s
% Computational Cost: add. (6872->464), mult. (14880->547), div. (0->0), fcn. (10875->14), ass. (0->247)
t184 = cos(pkin(9));
t302 = t184 * pkin(1);
t162 = -pkin(2) - t302;
t252 = qJDD(1) * t162;
t145 = qJDD(3) + t252;
t180 = qJ(1) + pkin(9);
t172 = cos(t180);
t170 = sin(t180);
t310 = g(1) * t170;
t238 = g(2) * t172 - t310;
t327 = -t145 - t238;
t183 = cos(pkin(10));
t314 = cos(qJ(4));
t244 = t314 * t183;
t156 = qJD(1) * t244;
t181 = sin(pkin(10));
t188 = sin(qJ(4));
t143 = t314 * t181 + t188 * t183;
t258 = qJD(4) * t188;
t243 = t181 * t258;
t199 = -qJD(1) * t243 + qJDD(1) * t143;
t320 = qJD(4) * t156 + t199;
t326 = qJD(4) * qJD(5) + t320;
t132 = t143 * qJD(1);
t135 = t143 * qJD(4);
t262 = t188 * t181;
t207 = t244 - t262;
t325 = t132 * t135 - t207 * t320;
t187 = sin(qJ(5));
t190 = cos(qJ(5));
t112 = qJD(4) * t187 + t132 * t190;
t256 = qJD(5) * t187;
t204 = t187 * qJDD(4) - t132 * t256 + t190 * t326;
t295 = qJ(6) * t204;
t251 = t181 * qJDD(1);
t222 = -qJDD(1) * t244 + t188 * t251;
t92 = qJD(1) * t135 + t222;
t89 = qJDD(5) + t92;
t315 = pkin(5) * t89;
t182 = sin(pkin(9));
t158 = pkin(1) * t182 + qJ(3);
t148 = t158 * qJD(1);
t121 = t181 * qJD(2) + t183 * t148;
t296 = pkin(7) * qJD(1);
t109 = t183 * t296 + t121;
t140 = qJD(1) * qJD(3) + qJDD(1) * t158;
t166 = t183 * qJDD(2);
t103 = t166 + (-pkin(7) * qJDD(1) - t140) * t181;
t115 = t181 * qJDD(2) + t183 * t140;
t250 = t183 * qJDD(1);
t104 = pkin(7) * t250 + t115;
t168 = t183 * qJD(2);
t108 = t168 + (-t148 - t296) * t181;
t240 = qJD(4) * t314;
t248 = -t188 * t103 - t314 * t104 - t108 * t240;
t33 = -t109 * t258 - t248;
t31 = qJDD(4) * pkin(8) + t33;
t65 = t188 * t108 + t314 * t109;
t62 = qJD(4) * pkin(8) + t65;
t161 = t183 * pkin(3) + pkin(2);
t147 = -t161 - t302;
t126 = qJD(1) * t147 + qJD(3);
t130 = qJD(1) * t262 - t156;
t71 = pkin(4) * t130 - pkin(8) * t132 + t126;
t37 = t187 * t71 + t190 * t62;
t285 = pkin(1) * qJDD(1);
t245 = t184 * t285;
t47 = -qJDD(1) * pkin(2) - pkin(3) * t250 + t92 * pkin(4) - pkin(8) * t320 + qJDD(3) - t245;
t4 = -qJD(5) * t37 - t187 * t31 + t190 * t47;
t1 = -qJD(6) * t112 - t295 + t315 + t4;
t124 = qJD(5) + t130;
t110 = -t190 * qJD(4) + t132 * t187;
t24 = -qJ(6) * t110 + t37;
t293 = t124 * t24;
t324 = t1 + t293;
t291 = t124 * t37;
t323 = t4 + t291;
t233 = t124 * t187;
t322 = t112 * t233;
t179 = pkin(10) + qJ(4);
t169 = sin(t179);
t228 = g(1) * t172 + g(2) * t170;
t209 = t228 * t169;
t301 = pkin(7) + t158;
t136 = t301 * t181;
t137 = t301 * t183;
t321 = -t314 * t136 - t188 * t137;
t264 = t172 * t190;
t171 = cos(t179);
t266 = t171 * t187;
t116 = t170 * t266 + t264;
t265 = t172 * t187;
t267 = t170 * t190;
t118 = -t171 * t265 + t267;
t305 = g(3) * t169;
t319 = -g(1) * t118 + g(2) * t116 + t187 * t305;
t134 = -t183 * t240 + t243;
t271 = t134 * t190;
t205 = t143 * t256 + t271;
t269 = t143 * t190;
t318 = -t124 * t205 + t89 * t269;
t202 = -t228 * t171 - t305;
t317 = t112 ^ 2;
t316 = t132 ^ 2;
t189 = sin(qJ(1));
t313 = pkin(1) * t189;
t312 = pkin(4) * t171;
t308 = g(1) * t189;
t304 = g(3) * t171;
t303 = t110 * pkin(5);
t185 = -qJ(6) - pkin(8);
t36 = -t187 * t62 + t190 * t71;
t23 = -qJ(6) * t112 + t36;
t20 = pkin(5) * t124 + t23;
t300 = -t23 + t20;
t106 = t188 * t109;
t64 = t314 * t108 - t106;
t90 = pkin(4) * t132 + pkin(8) * t130;
t43 = t187 * t90 + t190 * t64;
t255 = qJD(5) * t190;
t58 = -t190 * qJDD(4) + t132 * t255 + t187 * t326;
t299 = t110 * t271 - t58 * t269;
t298 = t112 * t135 - t204 * t207;
t86 = -t188 * t136 + t314 * t137;
t78 = t190 * t86;
t79 = -pkin(4) * t207 - pkin(8) * t143 + t147;
t51 = t187 * t79 + t78;
t297 = t134 * t130 - t143 * t92;
t294 = qJ(6) * t58;
t292 = t124 * t36;
t290 = t187 * t204;
t289 = t187 * t89;
t288 = -t110 * t255 - t187 * t58;
t236 = qJD(5) * t185;
t275 = t130 * t187;
t287 = -qJ(6) * t275 + qJD(6) * t190 + t187 * t236 - t43;
t42 = -t187 * t64 + t190 * t90;
t286 = -pkin(5) * t132 - qJD(6) * t187 - t42 + (-qJ(6) * t130 + t236) * t190;
t284 = t110 * t124;
t283 = t110 * t130;
t282 = t110 * t132;
t281 = t110 * t187;
t280 = t112 * t110;
t279 = t112 * t124;
t278 = t112 * t132;
t277 = t112 * t190;
t276 = t124 * t132;
t274 = t132 * t130;
t272 = t134 * t187;
t270 = t143 * t187;
t268 = t169 * t172;
t261 = (g(1) * t264 + g(2) * t267) * t169;
t191 = cos(qJ(1));
t176 = t191 * pkin(1);
t260 = t172 * t161 + t176;
t177 = t181 ^ 2;
t178 = t183 ^ 2;
t259 = t177 + t178;
t257 = qJD(5) * t112;
t67 = t207 * qJD(3) + qJD(4) * t321;
t91 = pkin(4) * t135 + pkin(8) * t134;
t249 = t187 * t91 + t190 * t67 + t79 * t255;
t246 = t112 * t272;
t242 = t143 * t255;
t34 = t314 * t103 - t188 * t104 - t108 * t258 - t109 * t240;
t32 = -qJDD(4) * pkin(4) - t34;
t15 = pkin(5) * t58 + qJDD(6) + t32;
t241 = -t15 - t304;
t186 = -pkin(7) - qJ(3);
t239 = pkin(5) * t187 - t186;
t237 = -t187 * t67 + t190 * t91;
t50 = -t187 * t86 + t190 * t79;
t3 = t187 * t47 + t190 * t31 + t71 * t255 - t62 * t256;
t235 = qJD(6) + t303;
t234 = t124 * t190;
t232 = pkin(8) * qJD(5) * t124 + t32;
t231 = g(2) * t268 - t169 * t310;
t230 = -g(1) * t116 - g(2) * t118;
t117 = -t171 * t267 + t265;
t119 = t170 * t187 + t171 * t264;
t229 = -g(1) * t117 - g(2) * t119;
t226 = -g(2) * t191 + t308;
t2 = -qJD(6) * t110 - t294 + t3;
t225 = -t124 * t20 + t2;
t224 = t3 - t292;
t223 = -t172 * t186 - t313;
t221 = t187 * t24 + t190 * t20;
t220 = t187 * t20 - t190 * t24;
t219 = t187 * t37 + t190 * t36;
t218 = t187 * t36 - t190 * t37;
t217 = -t110 * t135 + t207 * t58;
t114 = -t140 * t181 + t166;
t216 = -t114 * t181 + t115 * t183;
t215 = (-t148 * t181 + t168) * t181 - t121 * t183;
t164 = pkin(5) * t190 + pkin(4);
t214 = t164 * t171 - t169 * t185;
t212 = qJ(6) * t134 - qJD(6) * t143;
t210 = t190 * t89 + (-t256 - t275) * t124;
t61 = -qJD(4) * pkin(4) - t64;
t208 = -pkin(8) * t89 + t124 * t61;
t206 = t242 - t272;
t203 = -t252 + t327;
t123 = qJDD(1) * t147 + qJDD(3);
t201 = -t304 + t209;
t200 = g(1) * t119 - g(2) * t117 + t190 * t305 - t3;
t198 = -qJD(5) * t219 - t4 * t187 + t3 * t190;
t197 = -t124 * t206 - t89 * t270;
t194 = t4 + t319;
t68 = qJD(3) * t143 + qJD(4) * t86;
t154 = g(3) * t266;
t150 = t185 * t190;
t149 = t185 * t187;
t129 = t130 ^ 2;
t107 = t110 ^ 2;
t94 = -qJD(4) * t135 + qJDD(4) * t207;
t93 = -qJD(4) * t134 + qJDD(4) * t143;
t69 = pkin(5) * t270 - t321;
t55 = -t107 + t317;
t54 = -pkin(5) * t275 + t65;
t52 = t124 * t135 - t207 * t89;
t48 = t235 + t61;
t44 = pkin(5) * t206 + t68;
t41 = -qJ(6) * t270 + t51;
t40 = -t58 + t279;
t39 = t204 + t284;
t38 = -pkin(5) * t207 - qJ(6) * t269 + t50;
t28 = -t124 ^ 2 * t190 - t278 - t289;
t27 = t124 * t234 - t278 + t289;
t26 = t210 + t282;
t25 = t210 - t282;
t22 = t110 * t233 - t190 * t58;
t21 = t112 * t234 + t290;
t19 = t110 * t206 + t58 * t270;
t18 = -t112 * t205 + t204 * t269;
t17 = -t51 * qJD(5) + t237;
t16 = -t256 * t86 + t249;
t14 = -qJ(6) * t242 + (-qJD(5) * t86 + t212) * t187 + t249;
t13 = t197 - t217;
t12 = t197 + t217;
t11 = t298 - t318;
t10 = t298 + t318;
t9 = pkin(5) * t135 + t212 * t190 + (-t78 + (qJ(6) * t143 - t79) * t187) * qJD(5) + t237;
t8 = (t204 - t283) * t190 - t322 + t288;
t7 = (-t204 - t283) * t190 + t322 + t288;
t6 = t246 + (-t290 + (-t277 + t281) * qJD(5)) * t143 + t299;
t5 = -t246 + (t290 + (t277 + t281) * qJD(5)) * t143 + t299;
t29 = [0, 0, 0, 0, 0, qJDD(1), t226, g(1) * t191 + g(2) * t189, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t238 + 0.2e1 * t245, -0.2e1 * t182 * t285 + t228, 0 (t226 + (t182 ^ 2 + t184 ^ 2) * t285) * pkin(1), t177 * qJDD(1), 0.2e1 * t181 * t250, 0, t178 * qJDD(1), 0, 0, t203 * t183, -t203 * t181, t140 * t259 + t216 - t228, t145 * t162 - g(1) * (-pkin(2) * t170 + qJ(3) * t172 - t313) - g(2) * (pkin(2) * t172 + qJ(3) * t170 + t176) + t216 * t158 - t215 * qJD(3), -t132 * t134 + t143 * t320, t297 - t325, t93, t130 * t135 - t207 * t92, t94, 0, -qJD(4) * t68 + qJDD(4) * t321 - t123 * t207 + t126 * t135 + t147 * t92 - t171 * t238, -t67 * qJD(4) - t86 * qJDD(4) + t123 * t143 - t126 * t134 + t147 * t320 + t231, -t67 * t130 + t68 * t132 + t64 * t134 - t65 * t135 - t34 * t143 + t207 * t33 - t320 * t321 - t86 * t92 - t228, t33 * t86 + t65 * t67 + t34 * t321 - t64 * t68 + t123 * t147 - g(1) * (-t161 * t170 + t223) - g(2) * (-t170 * t186 + t260) t18, t6, t10, t19, t12, t52, -t61 * t272 + t110 * t68 + t124 * t17 + t135 * t36 - t207 * t4 + t50 * t89 - t58 * t321 + (t187 * t32 + t255 * t61) * t143 + t229, -t61 * t271 + t112 * t68 - t124 * t16 - t135 * t37 + t207 * t3 - t51 * t89 - t204 * t321 + (t190 * t32 - t256 * t61) * t143 + t230, -t110 * t16 - t112 * t17 - t50 * t204 - t51 * t58 + t219 * t134 + (qJD(5) * t218 - t187 * t3 - t190 * t4) * t143 - t231, t3 * t51 + t37 * t16 + t4 * t50 + t36 * t17 - t32 * t321 + t61 * t68 - g(1) * t223 - g(2) * (pkin(8) * t268 + t172 * t312 + t260) + (-g(1) * (-pkin(8) * t169 - t161 - t312) + g(2) * t186) * t170, t18, t6, t10, t19, t12, t52, -t48 * t272 - t1 * t207 + t110 * t44 + t124 * t9 + t135 * t20 + t38 * t89 + t58 * t69 + (t15 * t187 + t255 * t48) * t143 + t229, -t48 * t271 + t112 * t44 - t124 * t14 - t135 * t24 + t207 * t2 - t41 * t89 + t204 * t69 + (t15 * t190 - t256 * t48) * t143 + t230, -t110 * t14 - t112 * t9 - t38 * t204 - t41 * t58 + t221 * t134 + (qJD(5) * t220 - t1 * t190 - t187 * t2) * t143 - t231, t2 * t41 + t24 * t14 + t1 * t38 + t20 * t9 + t15 * t69 + t48 * t44 + pkin(1) * t308 - g(2) * t260 + (-g(1) * t239 - g(2) * t214) * t172 + (-g(1) * (-t161 - t214) - g(2) * t239) * t170; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2) - g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, t114 * t183 + t115 * t181 - g(3), 0, 0, 0, 0, 0, 0, t94, -t93, t297 + t325, -t134 * t65 - t135 * t64 + t143 * t33 + t207 * t34 - g(3), 0, 0, 0, 0, 0, 0, t13, t11, t5, t134 * t218 + t135 * t61 + t143 * t198 - t207 * t32 - g(3), 0, 0, 0, 0, 0, 0, t13, t11, t5, t135 * t48 - t207 * t15 - g(3) + t220 * t134 + (-qJD(5) * t221 - t1 * t187 + t190 * t2) * t143; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t250, t251, -t259 * qJD(1) ^ 2, qJD(1) * t215 - t327, 0, 0, 0, 0, 0, 0, 0.2e1 * t132 * qJD(4) + t222 (t156 - t130) * qJD(4) + t199, -t129 - t316, t130 * t65 + t132 * t64 + t123 + t238, 0, 0, 0, 0, 0, 0, t25, t28, t7, -t132 * t61 + t224 * t187 + t190 * t323 + t238, 0, 0, 0, 0, 0, 0, t25, t28, t7, -t132 * t48 + t225 * t187 + t190 * t324 + t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t274, -t129 + t316 (t156 + t130) * qJD(4) + t199, -t274, -t222, qJDD(4), qJD(4) * t65 - t126 * t132 + t201 + t34, t126 * t130 + (t64 + t106) * qJD(4) + t248 - t202, 0, 0, t21, t8, t27, t22, t26, -t276, -pkin(4) * t58 - t110 * t65 - t124 * t42 - t132 * t36 + (-t232 - t304) * t190 + t208 * t187 + t261, -pkin(4) * t204 - t112 * t65 + t124 * t43 + t132 * t37 + t154 + t208 * t190 + (-t209 + t232) * t187, t110 * t43 + t112 * t42 + ((-t58 + t257) * pkin(8) + t224) * t190 + ((qJD(5) * t110 + t204) * pkin(8) - t323) * t187 + t202, -t36 * t42 - t37 * t43 - t61 * t65 + (-t32 + t201) * pkin(4) + (t198 + t202) * pkin(8), t21, t8, t27, t22, t26, -t276, -t110 * t54 - t132 * t20 + t149 * t89 - t164 * t58 + t241 * t190 + t286 * t124 + (t130 * t48 + (t48 + t303) * qJD(5)) * t187 + t261, -t112 * t54 + t132 * t24 + t150 * t89 - t164 * t204 + t154 + t48 * t234 - t287 * t124 + (pkin(5) * t257 + t15 - t209) * t187, -t287 * t110 - t286 * t112 - t149 * t204 + t150 * t58 - t187 * t324 + t225 * t190 + t202, -t2 * t150 + t1 * t149 - t15 * t164 - g(3) * t214 + (pkin(5) * t256 - t54) * t48 + t287 * t24 + t286 * t20 + t228 * (t164 * t169 + t171 * t185); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t280, t55, t39, -t280, t40, t89, -t112 * t61 + t194 + t291, t110 * t61 + t200 + t292, 0, 0, t280, t55, t39, -t280, t40, t89, 0.2e1 * t315 - t295 + t293 + (-t235 - t48) * t112 + t194, -pkin(5) * t317 + t294 + t124 * t23 + (qJD(6) + t48) * t110 + t200, -pkin(5) * t204 - t110 * t300, t300 * t24 + (-t48 * t112 + t1 + t319) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 + t279, t204 - t284, -t107 - t317, t110 * t24 + t112 * t20 - t209 - t241;];
tau_reg  = t29;
