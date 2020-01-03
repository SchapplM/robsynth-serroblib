% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPP7
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPP7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP7_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:05:54
% EndTime: 2019-12-31 21:06:04
% DurationCPUTime: 3.31s
% Computational Cost: add. (4933->309), mult. (9872->327), div. (0->0), fcn. (6214->6), ass. (0->188)
t165 = sin(qJ(3));
t168 = cos(qJ(3));
t166 = sin(qJ(2));
t210 = qJD(1) * t166;
t135 = -t168 * qJD(2) + t165 * t210;
t155 = t166 * qJDD(1);
t169 = cos(qJ(2));
t206 = qJD(1) * qJD(2);
t202 = t169 * t206;
t140 = t155 + t202;
t186 = -t165 * qJDD(2) - t168 * t140;
t179 = t135 * qJD(3) + t186;
t152 = t169 * qJD(1) - qJD(3);
t216 = t152 * t135;
t263 = t216 - t179;
t137 = t165 * qJD(2) + t168 * t210;
t101 = t137 * t135;
t154 = t166 * t206;
t205 = t169 * qJDD(1);
t141 = -t154 + t205;
t134 = -qJDD(3) + t141;
t268 = -t101 + t134;
t275 = t168 * t268;
t133 = t137 ^ 2;
t252 = t152 ^ 2;
t90 = -t252 - t133;
t43 = -t165 * t90 + t275;
t313 = pkin(7) * t43;
t319 = pkin(2) * t263 - t313;
t312 = t169 * t43;
t279 = t165 * t268;
t40 = t168 * t90 + t279;
t316 = pkin(1) * t40;
t318 = t316 + pkin(6) * (-t166 * t263 - t312);
t315 = pkin(2) * t40;
t314 = pkin(7) * t40;
t264 = -t216 - t179;
t239 = t165 * t264;
t199 = -t168 * qJDD(2) + t165 * t140;
t55 = (qJD(3) + t152) * t137 + t199;
t27 = t168 * t55 - t239;
t253 = t135 ^ 2;
t71 = -t253 - t133;
t297 = t166 * t71;
t317 = pkin(6) * (t169 * t27 - t297);
t240 = t165 * t263;
t261 = -t253 + t133;
t217 = t137 * t152;
t87 = t137 * qJD(3) + t199;
t270 = -t217 + t87;
t302 = t169 * t261 + t166 * (t168 * t270 + t240);
t301 = pkin(2) * t71;
t311 = -pkin(7) * t27 - t301;
t260 = t253 - t252;
t269 = t217 + t87;
t287 = t166 * (t168 * t260 + t279) + t169 * t269;
t109 = t133 - t252;
t267 = t101 + t134;
t276 = t168 * t267;
t309 = t166 * (-t165 * t109 + t276) + t169 * t264;
t262 = -t252 - t253;
t291 = t165 * t262 - t276;
t308 = pkin(1) * t291;
t300 = pkin(2) * t291;
t299 = pkin(7) * t291;
t229 = t168 * t263;
t294 = -t165 * t270 + t229;
t306 = -t165 * t260 + t275;
t280 = t165 * t267;
t303 = t168 * t109 + t280;
t290 = t168 * t262 + t280;
t289 = -pkin(2) * t270 + pkin(7) * t290;
t298 = qJ(4) * t71;
t296 = qJ(4) * t263;
t228 = t168 * t264;
t24 = t165 * t55 + t228;
t232 = t166 * t270;
t286 = pkin(6) * (t169 * t290 + t232) - t308;
t171 = qJD(1) ^ 2;
t167 = sin(qJ(1));
t170 = cos(qJ(1));
t191 = t170 * g(1) + t167 * g(2);
t218 = qJDD(1) * pkin(6);
t125 = -t171 * pkin(1) - t191 + t218;
t194 = -t169 * pkin(2) - t166 * pkin(7);
t198 = t171 * t194 + t125;
t246 = t169 * g(3);
t250 = qJD(2) ^ 2;
t67 = -qJDD(2) * pkin(2) - t250 * pkin(7) + t198 * t166 + t246;
t285 = t87 * pkin(3) - t296 + t67;
t284 = qJ(4) * t268;
t283 = qJ(5) * t264;
t271 = -pkin(3) * t267 + qJ(4) * t262;
t214 = t152 * t168;
t204 = t135 * t214;
t94 = t169 * t101;
t178 = t166 * (t165 * t87 - t204) + t94;
t209 = qJD(4) * t152;
t145 = -0.2e1 * t209;
t208 = qJD(5) * t135;
t266 = 0.2e1 * t208 + t145;
t146 = 0.2e1 * t209;
t265 = -0.2e1 * t208 + t146;
t103 = t152 * pkin(4) - t137 * qJ(5);
t259 = t137 * t103 + qJDD(5);
t257 = -t87 * pkin(4) + t259;
t201 = t167 * g(1) - t170 * g(2);
t124 = qJDD(1) * pkin(1) + t171 * pkin(6) + t201;
t189 = -t141 + t154;
t190 = t140 + t202;
t52 = t189 * pkin(2) - t190 * pkin(7) - t124;
t247 = t166 * g(3);
t68 = -t250 * pkin(2) + qJDD(2) * pkin(7) + t198 * t169 - t247;
t31 = t165 * t68 - t168 * t52;
t184 = t134 * pkin(3) - qJ(4) * t252 + qJDD(4) + t31;
t176 = t134 * pkin(4) + t184 - t283;
t93 = t135 * pkin(3) - t137 * qJ(4);
t203 = -pkin(4) * t135 - t93;
t185 = (-0.2e1 * qJD(5) - t203) * t137;
t11 = t185 + t176;
t32 = t165 * t52 + t168 * t68;
t195 = -pkin(3) * t252 - t134 * qJ(4) - t135 * t93 + t32;
t182 = -pkin(4) * t253 + t87 * qJ(5) - t152 * t103 + t195;
t12 = t182 + t266;
t249 = pkin(3) + pkin(4);
t256 = qJ(4) * t12 - t249 * t11;
t215 = t152 * t165;
t106 = t137 * t215;
t245 = t166 * (-t168 * t179 + t106) - t94;
t255 = -t249 * t90 + t182 - t284;
t254 = qJ(4) * t55 + t249 * t264;
t251 = 0.2e1 * t137;
t242 = t137 * t93;
t238 = t165 * t67;
t227 = t168 * t67;
t220 = qJ(4) * t165;
t219 = qJ(4) * t168;
t151 = t169 * t171 * t166;
t212 = t166 * (qJDD(2) + t151);
t211 = t169 * (-t151 + qJDD(2));
t18 = t165 * t31 + t168 * t32;
t104 = t166 * t125 + t246;
t105 = t169 * t125 - t247;
t200 = t166 * t104 + t169 * t105;
t196 = -t135 * t215 - t168 * t87;
t47 = -t137 * t214 - t165 * t179;
t20 = t145 + t195;
t22 = t184 + t242;
t193 = -pkin(3) * t22 + qJ(4) * t20;
t192 = -pkin(3) * t264 - qJ(4) * t269;
t188 = t165 * t32 - t168 * t31;
t187 = -pkin(1) + t194;
t183 = (t135 * t165 + t137 * t168) * t152;
t119 = t169 * t134;
t180 = t166 * (-t106 + t204) + t119;
t177 = -pkin(3) * t90 + t195 - t284;
t174 = -t22 + t271;
t173 = qJD(4) * t251 - t285;
t172 = -t176 + t271;
t21 = (-pkin(3) * t152 - 0.2e1 * qJD(4)) * t137 + t285;
t15 = (-t270 + t217) * pkin(3) + t173;
t14 = pkin(3) * t217 + t173 + t296;
t162 = t169 ^ 2;
t161 = t166 ^ 2;
t159 = t162 * t171;
t157 = t161 * t171;
t142 = -0.2e1 * t154 + t205;
t139 = t155 + 0.2e1 * t202;
t122 = qJD(5) * t251;
t64 = (qJD(3) - t152) * t135 + t186;
t33 = -qJ(4) * t270 - qJ(5) * t267;
t28 = -t168 * t269 + t239;
t25 = -t165 * t269 - t228;
t23 = qJ(5) * t268 + t249 * t263;
t19 = t22 - t298;
t16 = -pkin(3) * t71 + t20;
t13 = qJ(5) * t253 + t21 - t257;
t10 = (-t253 - t90) * qJ(5) + t14 + t257;
t9 = t165 * t22 + t168 * t20;
t8 = t165 * t20 - t168 * t22;
t7 = t203 * t137 + t122 - t176 + t283 + t298;
t6 = (-t253 - t262) * qJ(5) + (-t270 - t87) * pkin(4) + t15 + t259;
t5 = -qJ(5) * t55 + t249 * t71 - t182 + t265;
t4 = -qJ(4) * t13 - qJ(5) * t11;
t3 = t165 * t11 + t168 * t12;
t2 = -t168 * t11 + t165 * t12;
t1 = -qJ(5) * t12 - t249 * t13;
t17 = [0, 0, 0, 0, 0, qJDD(1), t201, t191, 0, 0, t190 * t166, t169 * t139 + t166 * t142, t212 + t169 * (-t157 + t250), -t189 * t169, t166 * (t159 - t250) + t211, 0, t169 * t124 + pkin(1) * t142 + pkin(6) * (t169 * (-t159 - t250) - t212), -t166 * t124 - pkin(1) * t139 + pkin(6) * (-t211 - t166 * (-t157 - t250)), pkin(1) * (t157 + t159) + (t161 + t162) * t218 + t200, pkin(1) * t124 + pkin(6) * t200, t245, -t302, -t309, t178, t287, t180, t166 * (t238 - t299) + t169 * (t31 - t300) + t286, t166 * (t227 - t314) + t169 * (t32 - t315) - t316 + pkin(6) * (-t166 * t64 + t312), -t166 * t188 - t187 * t24 - t317, pkin(6) * (t166 * t67 + t169 * t18) + t187 * t188, t245, -t309, t302, t180, -t287, t178, t166 * (-t165 * t15 - t219 * t270 - t299) + t169 * (-t174 - t300) + t286, t166 * (-pkin(7) * t25 - t165 * t16 + t168 * t19) + t169 * (-pkin(2) * t25 - t192) - pkin(1) * t25 + pkin(6) * (t169 * t28 + t297), t166 * (-pkin(3) * t240 + t168 * t14 + t314) + t169 * (t146 - t177 + t315) + t318, t166 * (-pkin(7) * t8 + (pkin(3) * t165 - t219) * t21) + t169 * (-pkin(2) * t8 - t193) - pkin(1) * t8 + pkin(6) * (t166 * t21 + t169 * t9), t245, t302, t309, t178, t287, t119 + t166 * (t135 * t168 - t137 * t165) * t152, t166 * (-t165 * t6 + t168 * t33 - t299) - t308 + pkin(6) * t232 + (pkin(4) * t267 + pkin(6) * t290 - t172 + t185 - t300) * t169, t166 * (t168 * t10 - t165 * t23 + t314) + t169 * (-t255 + t265 + t315) + t318, t166 * (-pkin(7) * t24 - t165 * t5 + t168 * t7) + t169 * (-pkin(2) * t24 - t254) - pkin(1) * t24 + t317, t166 * (-pkin(7) * t2 - t165 * t1 + t168 * t4) + t169 * (-pkin(2) * t2 - t256) - pkin(1) * t2 + pkin(6) * (t166 * t13 + t169 * t3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t151, t157 - t159, t155, t151, t205, qJDD(2), -t104, -t105, 0, 0, t47, t294, -t303, t196, -t306, t183, -t227 + t289, pkin(2) * t64 + t238 + t313, t18 + t311, -pkin(2) * t67 + pkin(7) * t18, t47, -t303, -t294, t183, t306, t196, t168 * t15 - t220 * t270 + t289, pkin(7) * t28 + t168 * t16 + t165 * t19 - t301, pkin(3) * t229 + t165 * t14 + t319, pkin(7) * t9 + (-pkin(3) * t168 - pkin(2) - t220) * t21, t47, -t294, t303, t196, -t306, t183, t165 * t33 + t168 * t6 + t289, t165 * t10 + t168 * t23 + t319, t165 * t7 + t168 * t5 - t311, -pkin(2) * t13 + pkin(7) * t3 + t168 * t1 + t165 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t261, t264, -t101, -t269, -t134, -t31, -t32, 0, 0, t101, t264, -t261, -t134, t269, -t101, t174, t192, t145 + t177, t193, t101, -t261, -t264, -t101, -t269, -t134, -t242 + t122 + (-t267 - t101) * pkin(4) + t172, t255 + t266, t254, t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t267, t264, t90, t22, 0, 0, 0, 0, 0, 0, t267, t90, -t264, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t270, t263, t71, -t13;];
tauJ_reg = t17;
