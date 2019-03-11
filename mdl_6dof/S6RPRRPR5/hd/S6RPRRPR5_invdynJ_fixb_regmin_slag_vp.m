% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
% 
% Output:
% tau_reg [6x32]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:13:38
% EndTime: 2019-03-09 05:13:48
% DurationCPUTime: 4.24s
% Computational Cost: add. (6423->417), mult. (15579->506), div. (0->0), fcn. (12502->14), ass. (0->231)
t164 = qJD(3) + qJD(4);
t168 = sin(qJ(6));
t171 = cos(qJ(6));
t166 = sin(pkin(10));
t288 = sin(qJ(3));
t243 = t288 * t166;
t167 = cos(pkin(10));
t290 = cos(qJ(3));
t245 = t290 * t167;
t201 = t243 - t245;
t115 = t201 * qJD(1);
t121 = t290 * t166 + t288 * t167;
t116 = t121 * qJD(1);
t169 = sin(qJ(4));
t289 = cos(qJ(4));
t87 = t289 * t115 + t116 * t169;
t76 = t164 * t168 - t171 * t87;
t283 = t76 * t87;
t208 = -t169 * t115 + t289 * t116;
t300 = qJD(6) + t208;
t319 = t171 * t300;
t190 = t121 * qJDD(1);
t301 = t115 * qJD(3);
t179 = t190 - t301;
t239 = qJD(3) * t288;
t240 = qJD(3) * t290;
t246 = -qJDD(1) * t243 + (-t166 * t240 - t167 * t239) * qJD(1);
t199 = -qJDD(1) * t245 - t246;
t238 = qJD(4) * t289;
t253 = qJD(4) * t169;
t46 = t115 * t238 + t116 * t253 + t169 * t199 - t289 * t179;
t43 = -qJDD(6) + t46;
t325 = t168 * t43 - t300 * t319;
t331 = t325 - t283;
t78 = t164 * t171 + t168 * t87;
t282 = t78 * t87;
t320 = t168 * t300;
t327 = -t171 * t43 - t300 * t320;
t330 = t327 + t282;
t274 = t164 * t87;
t196 = -t46 + t274;
t269 = qJDD(1) * pkin(1);
t170 = sin(qJ(1));
t172 = cos(qJ(1));
t302 = -g(1) * t170 + g(2) * t172;
t211 = -qJDD(2) + t269 - t302;
t280 = t87 ^ 2;
t311 = t208 ^ 2;
t329 = -t280 + t311;
t159 = qJDD(3) + qJDD(4);
t251 = qJD(6) * t171;
t252 = qJD(6) * t168;
t178 = t169 * t179 + t289 * t199;
t47 = qJD(4) * t208 + t178;
t20 = t171 * t159 - t164 * t252 + t168 * t47 + t87 * t251;
t19 = t20 * t171;
t328 = -t320 * t78 + t19;
t232 = t159 * t168 - t171 * t47;
t21 = qJD(6) * t78 + t232;
t323 = t300 * t76;
t326 = -t319 * t78 + (-t20 + t323) * t168 - t171 * t21;
t292 = pkin(5) * t87;
t277 = qJ(5) * t87;
t279 = t208 * t87;
t278 = pkin(7) + qJ(2);
t133 = t278 * t166;
t122 = qJD(1) * t133;
t134 = t278 * t167;
t123 = qJD(1) * t134;
t75 = -t115 * pkin(8) - t288 * t122 + t290 * t123;
t72 = t289 * t75;
t304 = -t290 * t122 - t288 * t123;
t74 = -t116 * pkin(8) + t304;
t73 = qJD(3) * pkin(3) + t74;
t50 = t169 * t73 + t72;
t48 = -qJ(5) * t164 - t50;
t24 = -t48 - t292;
t324 = t24 * t300;
t322 = t300 * t87;
t71 = t169 * t75;
t49 = -t289 * t73 + t71;
t259 = qJD(5) + t49;
t163 = pkin(10) + qJ(3);
t156 = cos(t163);
t157 = qJ(4) + t163;
t147 = sin(t157);
t148 = cos(t157);
t255 = t148 * pkin(4) + t147 * qJ(5);
t318 = pkin(3) * t156 + t255;
t177 = t115 * t253 - t116 * t238 - t178;
t272 = t208 * t164;
t317 = t177 + t272;
t309 = pkin(5) * t208;
t260 = t309 + t259;
t293 = pkin(4) + pkin(9);
t22 = -t293 * t164 + t260;
t150 = -t167 * pkin(2) - pkin(1);
t126 = t150 * qJD(1) + qJD(2);
t97 = pkin(3) * t115 + t126;
t186 = -qJ(5) * t208 + t97;
t27 = t293 * t87 + t186;
t11 = t168 * t22 + t171 * t27;
t249 = qJD(1) * qJD(2);
t296 = t278 * qJDD(1) + t249;
t101 = t296 * t166;
t102 = t296 * t167;
t222 = -t290 * t101 - t288 * t102;
t40 = qJDD(3) * pkin(3) - t179 * pkin(8) + t122 * t239 - t123 * t240 + t222;
t202 = -t288 * t101 + t290 * t102;
t45 = -t199 * pkin(8) + qJD(3) * t304 + t202;
t233 = -t169 * t40 - t73 * t238 + t75 * t253 - t289 * t45;
t152 = t159 * qJ(5);
t303 = -t164 * qJD(5) - t152;
t7 = t233 + t303;
t5 = -pkin(5) * t47 - t7;
t316 = -t11 * t87 + t5 * t171;
t10 = -t168 * t27 + t171 * t22;
t315 = t10 * t87 + t5 * t168 + t24 * t251;
t142 = g(3) * t147;
t265 = t148 * t172;
t266 = t148 * t170;
t206 = -g(1) * t265 - g(2) * t266 - t142 - t233;
t51 = pkin(4) * t87 + t186;
t314 = -t51 * t87 + t206 - t303;
t313 = t87 * t97 - t206;
t310 = pkin(4) * t208;
t308 = t208 * t24;
t53 = t289 * t74 - t71;
t270 = pkin(3) * t238 + qJD(5) - t53;
t52 = t169 * t74 + t72;
t223 = pkin(3) * t253 - t52;
t306 = t208 * t293;
t305 = -qJD(6) + t300;
t299 = qJ(2) * qJDD(1);
t143 = g(3) * t148;
t234 = t169 * t45 + t75 * t238 + t73 * t253 - t289 * t40;
t267 = t147 * t172;
t268 = t147 * t170;
t207 = -g(1) * t267 - g(2) * t268 + t143 + t234;
t188 = t51 * t208 + qJDD(5) + t207;
t298 = -t97 * t208 - t207;
t151 = -t289 * pkin(3) - pkin(4);
t145 = -pkin(9) + t151;
t297 = (t223 + t292) * t300 - t145 * t43;
t294 = qJD(3) ^ 2;
t287 = pkin(4) * t159;
t187 = t289 * t201;
t95 = t121 * t169 + t187;
t285 = t24 * t95;
t194 = t169 * t201;
t96 = t289 * t121 - t194;
t99 = pkin(3) * t201 + t150;
t181 = -t96 * qJ(5) + t99;
t31 = t293 * t95 + t181;
t284 = t31 * t43;
t275 = t164 * t50;
t271 = t270 + t309;
t264 = t168 * t170;
t263 = t168 * t172;
t262 = t170 * t171;
t261 = t171 * t172;
t256 = -t288 * t133 + t290 * t134;
t254 = t166 ^ 2 + t167 ^ 2;
t250 = t116 * qJD(3);
t247 = t95 * t251;
t242 = qJD(2) * t290;
t241 = qJD(2) * t288;
t79 = qJDD(2) - t246 * pkin(3) + (-pkin(1) + (-t290 * pkin(3) - pkin(2)) * t167) * qJDD(1);
t180 = t46 * qJ(5) - qJD(5) * t208 + t79;
t6 = t293 * t47 + t180;
t236 = qJD(6) * t22 + t6;
t216 = qJDD(5) + t234;
t4 = -pkin(5) * t46 - t293 * t159 + t216;
t235 = -qJD(6) * t27 + t4;
t227 = t254 * qJD(1) ^ 2;
t225 = 0.2e1 * t254;
t155 = sin(t163);
t221 = -pkin(3) * t155 - pkin(4) * t147;
t220 = g(1) * t172 + g(2) * t170;
t191 = qJD(3) * t201;
t192 = qJD(3) * t121;
t60 = -qJD(4) * t194 + t121 * t238 - t169 * t191 + t289 * t192;
t218 = t300 * t60 - t43 * t95;
t217 = -t290 * t133 - t288 * t134;
t215 = pkin(3) * t116 + t277;
t214 = -t293 * t43 + (t50 - t292) * t300;
t184 = -t133 * t240 - t134 * t239 - t166 * t241 + t167 * t242;
t63 = -pkin(8) * t192 + t184;
t64 = pkin(8) * t191 + t133 * t239 - t134 * t240 - t166 * t242 - t167 * t241;
t83 = -t121 * pkin(8) + t217;
t84 = -pkin(8) * t201 + t256;
t15 = -t169 * t64 - t83 * t238 + t84 * t253 - t289 * t63;
t57 = t169 * t83 + t289 * t84;
t213 = -t15 * t164 + t57 * t159;
t16 = t57 * qJD(4) + t169 * t63 - t289 * t64;
t56 = t169 * t84 - t289 * t83;
t212 = t56 * t159 + t16 * t164;
t209 = -t150 + t318;
t32 = t96 * pkin(5) + t56;
t205 = t24 * t60 + t32 * t43 + t5 * t95;
t198 = t211 + t269;
t195 = t46 + t274;
t193 = -t148 * t220 - t142;
t189 = pkin(3) * t192;
t185 = t225 * t249 - t220;
t183 = (-qJD(6) * t145 + t215 + t306) * t300 + t193;
t182 = (qJD(6) * t293 + t277 + t306) * t300 + t193;
t59 = t121 * t253 + t164 * t187 + t169 * t192;
t17 = t60 * pkin(4) + t59 * qJ(5) - t96 * qJD(5) + t189;
t9 = t47 * pkin(4) + t180;
t176 = -t177 + t272;
t160 = -pkin(8) - t278;
t149 = pkin(3) * t169 + qJ(5);
t128 = qJ(5) * t265;
t127 = qJ(5) * t266;
t125 = t150 * qJDD(1) + qJDD(2);
t111 = -t147 * t264 + t261;
t110 = t147 * t262 + t263;
t109 = t147 * t263 + t262;
t108 = t147 * t261 - t264;
t58 = t277 + t310;
t55 = t95 * pkin(4) + t181;
t54 = t215 + t310;
t44 = -pkin(4) * t164 + t259;
t33 = -t95 * pkin(5) + t57;
t14 = t60 * pkin(9) + t17;
t13 = -t59 * pkin(5) + t16;
t12 = -pkin(5) * t60 - t15;
t8 = t216 - t287;
t1 = t171 * t4;
t2 = [qJDD(1), -t302, t220, t198 * t167, -t198 * t166, t225 * t299 + t185, t211 * pkin(1) + (t254 * t299 + t185) * qJ(2), t121 * t190 + (-t115 * t121 - t116 * t201) * qJD(3), t115 * t191 - t116 * t192 - t121 * t199 - t179 * t201, t121 * qJDD(3) - t294 * t201, -qJDD(3) * t201 - t294 * t121, 0, t217 * qJDD(3) + t125 * t201 + t150 * t199 - t302 * t156 + (-t256 * qJD(3) + (-qJD(2) + t126) * t121) * qJD(3), -t184 * qJD(3) - t256 * qJDD(3) + t125 * t121 - t126 * t191 + t150 * t179 + t155 * t302, -t208 * t59 - t46 * t96, -t208 * t60 + t46 * t95 - t47 * t96 + t59 * t87, t159 * t96 - t164 * t59, -t159 * t95 - t164 * t60, 0, g(1) * t266 - g(2) * t265 + t189 * t87 + t99 * t47 + t97 * t60 + t79 * t95 - t212, -g(1) * t268 + g(2) * t267 + t189 * t208 - t99 * t46 - t97 * t59 + t79 * t96 - t213, t15 * t87 + t16 * t208 - t44 * t59 - t46 * t56 - t47 * t57 + t48 * t60 + t7 * t95 + t8 * t96 - t220, t148 * t302 - t17 * t87 - t47 * t55 - t51 * t60 - t9 * t95 + t212, -t147 * t302 - t17 * t208 + t46 * t55 + t51 * t59 - t9 * t96 + t213, t48 * t15 + t44 * t16 + t51 * t17 + t9 * t55 + t8 * t56 - t7 * t57 + (g(1) * t160 - g(2) * t209) * t172 + (g(1) * t209 + g(2) * t160) * t170, t78 * t247 + (t20 * t95 + t60 * t78) * t168 (-t168 * t76 + t171 * t78) * t60 + (-t168 * t21 + t19 + (-t168 * t78 - t171 * t76) * qJD(6)) * t95, t168 * t218 + t20 * t96 + t247 * t300 - t59 * t78, -t252 * t300 * t95 + t171 * t218 - t21 * t96 + t59 * t76, -t300 * t59 - t43 * t96, -g(1) * t111 - g(2) * t109 + t1 * t96 - t10 * t59 + t12 * t76 + t33 * t21 + (-t14 * t300 - t6 * t96 + t284) * t168 + (t13 * t300 - t205) * t171 + ((-t168 * t32 - t171 * t31) * t300 - t11 * t96 + t168 * t285) * qJD(6), g(1) * t110 - g(2) * t108 + t11 * t59 + t12 * t78 + t33 * t20 + (-(qJD(6) * t32 + t14) * t300 + t284 - t236 * t96 + qJD(6) * t285) * t171 + (-(-qJD(6) * t31 + t13) * t300 - t235 * t96 + t205) * t168; 0, 0, 0, -t167 * qJDD(1), t166 * qJDD(1), -t227, -qJ(2) * t227 - t211, 0, 0, 0, 0, 0, t199 + t250, t190 - 0.2e1 * t301, 0, 0, 0, 0, 0, t176, -t195, -t280 - t311, -t176, t195, -t208 * t44 - t48 * t87 + t302 + t9, 0, 0, 0, 0, 0, t325 + t283, t282 - t327; 0, 0, 0, 0, 0, 0, 0, t116 * t115, -t115 ^ 2 + t116 ^ 2, t190, -t199 + t250, qJDD(3), -g(3) * t156 - t126 * t116 + t220 * t155 + t222, g(3) * t155 + t126 * t115 + t220 * t156 - t202, t279, t329, t196, t317, t159, t52 * t164 + (-t116 * t87 + t289 * t159 - t164 * t253) * pkin(3) + t298, t53 * t164 + (-t116 * t208 - t159 * t169 - t164 * t238) * pkin(3) + t313, -t149 * t47 - t151 * t46 + (t223 - t48) * t208 + (-t270 + t44) * t87, t54 * t87 + t223 * t164 + (-pkin(4) + t151) * t159 + t188, t149 * t159 + t270 * t164 + t208 * t54 + t314, -t7 * t149 + t8 * t151 - t51 * t54 - g(1) * (t172 * t221 + t128) - g(2) * (t170 * t221 + t127) - g(3) * t318 - t270 * t48 + t223 * t44, t328, t326, t330, t331, t322, t149 * t21 + t271 * t76 + (t297 + t308) * t171 + t183 * t168 + t315, t149 * t20 + t271 * t78 + t183 * t171 + (-t297 - t324) * t168 + t316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t279, t329, t196, t317, t159, t275 + t298, -t164 * t49 + t313, pkin(4) * t46 - qJ(5) * t47 + (-t48 - t50) * t208 + (t44 - t259) * t87, t58 * t87 + t188 - t275 - 0.2e1 * t287, t259 * t164 + t208 * t58 + t152 + t314, -t7 * qJ(5) - t8 * pkin(4) - t51 * t58 - t44 * t50 - g(1) * (-pkin(4) * t267 + t128) - g(2) * (-pkin(4) * t268 + t127) - g(3) * t255 - t259 * t48, t328, t326, t330, t331, t322, qJ(5) * t21 + t260 * t76 + (-t214 + t308) * t171 + t182 * t168 + t315, qJ(5) * t20 + t260 * t78 + (t214 - t324) * t168 + t182 * t171 + t316; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, t159 - t279, -t164 ^ 2 - t311, t164 * t48 + t188 - t287, 0, 0, 0, 0, 0, -t164 * t76 + t327, -t164 * t78 + t325; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78 * t76, -t76 ^ 2 + t78 ^ 2, t20 + t323, t305 * t78 - t232, -t43, -g(1) * t108 - g(2) * t110 + t11 * t305 + t171 * t143 - t168 * t6 - t24 * t78 + t1, g(1) * t109 - g(2) * t111 + t10 * t300 + t24 * t76 - t236 * t171 + (-t235 - t143) * t168;];
tau_reg  = t2;
