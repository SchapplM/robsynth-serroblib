% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRPRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 04:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRPRP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 04:10:56
% EndTime: 2019-05-05 04:11:04
% DurationCPUTime: 3.00s
% Computational Cost: add. (6218->314), mult. (12476->373), div. (0->0), fcn. (7884->10), ass. (0->210)
t185 = cos(qJ(5));
t182 = sin(qJ(5));
t186 = cos(qJ(3));
t227 = qJD(2) * qJD(3);
t165 = t186 * t227;
t183 = sin(qJ(3));
t224 = t183 * qJDD(2);
t139 = t165 + t224;
t128 = qJDD(5) + t139;
t229 = qJD(2) * t186;
t133 = qJD(3) * t182 + t185 * t229;
t135 = qJD(3) * t185 - t182 * t229;
t246 = t135 * t133;
t274 = t128 + t246;
t260 = t182 * t274;
t127 = t135 ^ 2;
t230 = qJD(2) * t183;
t161 = qJD(5) + t230;
t269 = t161 ^ 2;
t276 = -t127 - t269;
t41 = -t185 * t276 + t260;
t325 = pkin(4) * t41;
t324 = t183 * t41;
t323 = t186 * t41;
t187 = cos(qJ(2));
t253 = t185 * t274;
t43 = t182 * t276 + t253;
t322 = t187 * t43;
t236 = t183 * qJ(4);
t266 = -pkin(3) - pkin(9);
t199 = t186 * t266 - pkin(2) - t236;
t321 = t199 * t43;
t320 = t266 * t41;
t247 = t133 * t161;
t218 = t183 * t227;
t223 = t186 * qJDD(2);
t140 = -t218 + t223;
t203 = t185 * qJDD(3) - t182 * t140;
t86 = -qJD(5) * t133 + t203;
t306 = t247 - t86;
t319 = t306 * qJ(6);
t270 = t133 ^ 2;
t109 = t270 - t269;
t116 = t135 * t161;
t214 = t182 * qJDD(3) + t185 * t140;
t200 = qJD(5) * t135 + t214;
t58 = -t116 + t200;
t318 = t183 * t58 + t186 * (t109 * t182 + t253);
t316 = -t109 * t185 + t260;
t171 = t183 ^ 2;
t189 = qJD(2) ^ 2;
t167 = t171 * t189;
t188 = qJD(3) ^ 2;
t155 = -t167 - t188;
t233 = t186 * t189;
t219 = t183 * t233;
t149 = qJDD(3) - t219;
t234 = t186 * t149;
t104 = t155 * t183 + t234;
t138 = 0.2e1 * t165 + t224;
t177 = sin(pkin(6));
t179 = cos(pkin(6));
t184 = sin(qJ(2));
t315 = t179 * (t149 * t183 - t155 * t186) + (t104 * t184 + t138 * t187) * t177;
t255 = t185 * t306;
t57 = t116 + t200;
t95 = t127 - t270;
t314 = t186 * (-t182 * t57 - t255) - t183 * t95;
t268 = 2 * qJD(4);
t72 = -t270 - t127;
t313 = pkin(4) * t72;
t275 = t128 - t246;
t252 = t185 * t275;
t272 = -t269 - t270;
t283 = t182 * t272 + t252;
t312 = pkin(4) * t283;
t311 = pkin(8) * t104;
t310 = qJ(4) * t72;
t309 = t183 * t72;
t307 = t186 * t72;
t303 = t183 * t283;
t302 = t186 * t283;
t69 = t182 * t275;
t285 = t185 * t272 - t69;
t301 = t187 * t285;
t300 = t199 * t285;
t299 = t266 * t283;
t111 = -t127 + t269;
t298 = -t111 * t182 + t252;
t172 = t186 ^ 2;
t168 = t172 * t189;
t157 = -t168 - t188;
t148 = qJDD(3) + t219;
t242 = t148 * t183;
t103 = -t157 * t186 + t242;
t141 = -0.2e1 * t218 + t223;
t295 = (t103 * t184 - t141 * t187) * t177 - t179 * (t148 * t186 + t157 * t183);
t279 = t247 + t86;
t294 = t183 * t279 + t186 * (-t111 * t185 - t69);
t293 = pkin(8) * t103;
t288 = t139 + t165;
t208 = -pkin(3) * t186 - t236;
t176 = sin(pkin(10));
t178 = cos(pkin(10));
t147 = -g(1) * t178 - g(2) * t176;
t146 = g(1) * t176 - g(2) * t178;
t173 = -g(3) + qJDD(1);
t286 = t146 * t179 + t173 * t177;
t80 = t187 * t147 + t286 * t184;
t67 = -t189 * pkin(2) + qJDD(2) * pkin(8) + t80;
t215 = t189 * t208 + t67;
t207 = t215 * t186;
t287 = qJD(3) * t268 + t207;
t282 = -t182 * t306 + t185 * t57;
t281 = qJ(6) * t182 + pkin(4);
t273 = t234 + (t168 - t188) * t183;
t271 = pkin(5) * t200 + t319;
t267 = 2 * qJD(6);
t113 = -t146 * t177 + t173 * t179;
t106 = t186 * t113;
t213 = -qJDD(3) * pkin(3) - t188 * qJ(4) + qJDD(4) - t106;
t31 = -qJDD(3) * pkin(9) + (t139 - t165) * pkin(4) + (-pkin(9) * t233 + t215) * t183 + t213;
t152 = pkin(4) * t230 - qJD(3) * pkin(9);
t209 = t184 * t147 - t286 * t187;
t66 = -qJDD(2) * pkin(2) - t189 * pkin(8) + t209;
t194 = -t140 * pkin(3) - t288 * qJ(4) + t66;
t217 = pkin(3) * qJD(3) - (2 * qJD(4));
t33 = -pkin(4) * t168 - t140 * pkin(9) + (-t152 + t217) * t230 + t194;
t18 = t182 * t31 + t185 * t33;
t93 = pkin(5) * t133 - qJ(6) * t135;
t211 = t128 * qJ(6) - t133 * t93 + t161 * t267 + t18;
t10 = -pkin(5) * t269 + t211;
t17 = t182 * t33 - t185 * t31;
t13 = -t128 * pkin(5) - qJ(6) * t269 + t135 * t93 + qJDD(6) + t17;
t265 = -pkin(5) * t13 + qJ(6) * t10;
t264 = pkin(5) * t182;
t263 = -pkin(5) * t279 - qJ(6) * t58;
t180 = t188 * pkin(3);
t225 = qJDD(3) * qJ(4);
t235 = t183 * t113;
t197 = t140 * pkin(4) - pkin(9) * t168 - t180 + t225 + t235;
t30 = t207 + (t268 + t152) * qJD(3) + t197;
t262 = t182 * t30;
t254 = t185 * t279;
t60 = (-qJD(5) + t161) * t135 - t214;
t34 = t182 * t60 - t254;
t259 = t183 * t34;
t257 = t185 * t30;
t19 = (pkin(5) * t161 - (2 * qJD(6))) * t135 + t30 + t271;
t251 = t186 * t19;
t239 = t161 * t182;
t238 = t161 * t185;
t143 = (t171 + t172) * qJDD(2);
t144 = t167 + t168;
t231 = pkin(2) * t144 + pkin(8) * t143;
t228 = qJD(5) + t161;
t222 = t133 * t239;
t221 = t133 * t238;
t220 = t183 * t246;
t50 = t183 * t67 - t106;
t51 = t186 * t67 + t235;
t23 = t183 * t50 + t186 * t51;
t216 = -qJ(6) * t185 + qJ(4);
t108 = t135 * t238;
t212 = t186 * (-t182 * t86 - t108) + t220;
t107 = t135 * t239;
t210 = t107 - t221;
t5 = -t185 * t17 + t182 * t18;
t4 = t183 * t5 + t186 * t30;
t6 = t182 * t17 + t185 * t18;
t206 = (-t167 + t188) * t186 + t242;
t202 = pkin(2) - t208;
t201 = t182 * t200 + t221;
t198 = t183 * t128 + t186 * (t108 + t222);
t196 = -pkin(5) * t276 + qJ(6) * t274 + t10;
t195 = -t180 + t287;
t39 = t215 * t183 + t213;
t193 = pkin(5) * t275 + qJ(6) * t272 - t13;
t192 = -t220 + t186 * (t185 * t200 - t222);
t191 = t195 + t225;
t40 = t217 * t230 + t194;
t190 = -qJD(3) * t152 + t135 * t267 - t197 - t271 - t287;
t145 = t167 - t168;
t101 = t288 * t183;
t100 = (t140 - t218) * t186;
t91 = t138 * t186 + t141 * t183;
t89 = (t143 * t184 + t144 * t187) * t177;
t63 = -t228 * t133 + t203;
t59 = t228 * t135 + t214;
t54 = t182 * t279;
t53 = t185 * t86 - t107;
t38 = t191 + t235;
t37 = -t185 * t58 + t54;
t36 = t185 * t60 + t54;
t35 = -t182 * t58 - t254;
t27 = t186 * t63 - t324;
t26 = t186 * t59 + t303;
t25 = t186 * t306 + t324;
t24 = t186 * t57 + t303;
t22 = t183 * t35 + t307;
t21 = t259 + t307;
t20 = t183 * t39 + t186 * t38;
t15 = t190 + (-t59 - t116) * pkin(5);
t14 = -pkin(5) * t116 + t190 - t319;
t8 = -qJ(6) * t72 + t13;
t7 = (-t269 - t72) * pkin(5) + t211;
t3 = t10 * t185 + t182 * t13;
t2 = t10 * t182 - t13 * t185;
t1 = t183 * t2 + t251;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t173, 0, 0, 0, 0, 0, 0, (qJDD(2) * t187 - t184 * t189) * t177, (-qJDD(2) * t184 - t187 * t189) * t177, 0, t179 * t113 + (t184 * t80 - t187 * t209) * t177, 0, 0, 0, 0, 0, 0, -t295, -t315, t89, t179 * (t183 * t51 - t186 * t50) + (t184 * t23 - t187 * t66) * t177, 0, 0, 0, 0, 0, 0, t89, t295, t315, t179 * (t183 * t38 - t186 * t39) + (t184 * t20 - t187 * t40) * t177, 0, 0, 0, 0, 0, 0, t179 * (t183 * t57 - t302) + (t184 * t24 - t301) * t177, t179 * (t183 * t63 + t323) + (t184 * t27 + t322) * t177, t179 * (-t186 * t34 + t309) + (t184 * t21 - t187 * t36) * t177, t179 * (t183 * t30 - t186 * t5) + (t184 * t4 - t187 * t6) * t177, 0, 0, 0, 0, 0, 0, t179 * (t183 * t59 - t302) + (t184 * t26 - t301) * t177, t179 * (-t186 * t35 + t309) + (t184 * t22 - t187 * t37) * t177, t179 * (t183 * t306 - t323) + (t184 * t25 - t322) * t177, t179 * (t183 * t19 - t186 * t2) + (t1 * t184 - t187 * t3) * t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t209, -t80, 0, 0, t101, t91, t206, t100, t273, 0, pkin(2) * t141 - t186 * t66 - t293, -pkin(2) * t138 + t183 * t66 - t311, t23 + t231, -pkin(2) * t66 + pkin(8) * t23, 0, -t206, -t273, t101, t91, t100, (pkin(3) * t144 + t191) * t186 + (qJ(4) * t144 + t106 + t39) * t183 + t231, -t141 * t202 + t186 * t40 + t293, t183 * (-pkin(3) * t218 + t230 * t268 - t194) + t311 + t202 * t138, pkin(8) * t20 - t202 * t40, t212, -t314, t294, t192, -t318, t198, t183 * (-t17 + t312) + t186 * (pkin(4) * t57 + t257) + pkin(8) * t24 + t300, t183 * (-t18 - t325) + t186 * (pkin(4) * t63 - t262) + pkin(8) * t27 - t321, pkin(4) * t259 + t186 * (-t6 + t313) + pkin(8) * t21 + t199 * t36, t199 * t6 + (pkin(4) + pkin(8)) * t4, t212, t294, t314, t198, t318, t192, t183 * (t193 + t312) + t186 * (-t185 * t15 + t281 * t59) + pkin(8) * t26 + t300, t183 * (pkin(4) * t35 + t263) + t186 * (-t182 * t8 - t185 * t7 + t313) + pkin(8) * t22 + t199 * t37, t183 * (t196 + t325) + t186 * (pkin(4) * t306 + pkin(5) * t255 - t182 * t14) + pkin(8) * t25 + t321, t183 * (pkin(4) * t2 + t265) + pkin(8) * t1 + t199 * t3 + (pkin(5) * t185 + t281) * t251; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t219, t145, t224, t219, t223, qJDD(3), -t50, -t51, 0, 0, qJDD(3), -t224, -t223, -t219, t145, t219, (-pkin(3) * t183 + qJ(4) * t186) * qJDD(2), -pkin(3) * t148 - qJ(4) * t157 + t39, -pkin(3) * t155 + t235 + (qJDD(3) + t149) * qJ(4) + t195, -pkin(3) * t39 + qJ(4) * t38, t53, -t282, t298, t201, -t316, t210, qJ(4) * t57 + t262 + t299, qJ(4) * t63 + t257 - t320, t266 * t34 + t310 - t5, qJ(4) * t30 + t266 * t5, t53, t298, t282, t210, t316, t201, -t182 * t15 + t216 * t59 + t299, -t182 * t7 + t185 * t8 + t266 * t35 + t310, t185 * t14 - (-qJ(4) - t264) * t306 + t320, t266 * t2 + (t216 + t264) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t224, t148, t155, t39, 0, 0, 0, 0, 0, 0, t283, -t41, t34, t5, 0, 0, 0, 0, 0, 0, t283, t35, t41, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t246, t95, t279, -t246, -t58, t128, -t17, -t18, 0, 0, t246, t279, -t95, t128, t58, -t246, t193, t263, t196, t265; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t275, t279, t276, t13;];
tauJ_reg  = t9;
