% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPP8
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPP8_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP8_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP8_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:09:24
% EndTime: 2019-12-31 21:09:33
% DurationCPUTime: 3.01s
% Computational Cost: add. (4934->303), mult. (9835->327), div. (0->0), fcn. (6172->6), ass. (0->185)
t170 = cos(qJ(2));
t153 = t170 * qJD(1) - qJD(3);
t212 = qJD(4) * t153;
t169 = cos(qJ(3));
t151 = t153 ^ 2;
t166 = sin(qJ(3));
t167 = sin(qJ(2));
t214 = qJD(1) * t167;
t137 = t166 * qJD(2) + t169 * t214;
t253 = t137 ^ 2;
t259 = -t253 - t151;
t209 = qJD(1) * qJD(2);
t155 = t167 * t209;
t208 = t170 * qJDD(1);
t142 = -t155 + t208;
t134 = -qJDD(3) + t142;
t135 = -t169 * qJD(2) + t166 * t214;
t221 = t137 * t135;
t183 = t134 - t221;
t264 = t166 * t183;
t40 = t169 * t259 + t264;
t305 = pkin(2) * t40;
t312 = 0.2e1 * t212 + t305;
t112 = -t253 + t151;
t182 = t134 + t221;
t229 = t169 * t182;
t280 = t167 * (t166 * t112 + t229);
t156 = t167 * qJDD(1);
t204 = t170 * t209;
t141 = t156 + t204;
t187 = t166 * qJDD(2) + t169 * t141;
t61 = (-qJD(3) - t153) * t135 + t187;
t311 = t170 * t61 + t280;
t220 = t137 * t153;
t89 = t137 * qJD(3) - t169 * qJDD(2) + t166 * t141;
t266 = t89 - t220;
t236 = t166 * t182;
t254 = t135 ^ 2;
t260 = -t151 - t254;
t284 = t169 * t260 + t236;
t36 = -t166 * t260 + t229;
t310 = pkin(6) * (t167 * t266 + t170 * t284) + pkin(1) * t36;
t262 = t169 * t183;
t277 = t166 * t259 - t262;
t65 = (qJD(3) - t153) * t135 - t187;
t309 = pkin(6) * (t167 * t65 + t170 * t277) + pkin(1) * t40;
t306 = pkin(2) * t36;
t304 = pkin(7) * t36;
t303 = pkin(7) * t40;
t302 = pkin(7) * t277;
t301 = pkin(7) * t284;
t285 = t169 * t112 - t236;
t299 = -pkin(2) * t266 + t301;
t298 = pkin(2) * t65 - t302;
t111 = t254 - t151;
t267 = t220 + t89;
t297 = t167 * (t169 * t111 + t264) + t170 * t267;
t234 = t169 * t266;
t98 = t254 - t253;
t286 = t170 * t98;
t296 = t167 * (-t166 * t65 + t234) - t286;
t74 = t253 + t254;
t290 = pkin(2) * t74;
t289 = (t266 - t220) * pkin(3);
t288 = qJ(4) * t74;
t287 = t167 * t74;
t117 = t135 * t153;
t90 = -t135 * qJD(3) + t187;
t63 = t117 + t90;
t242 = t166 * t267;
t64 = -t117 + t90;
t25 = -t169 * t64 - t242;
t243 = t166 * t266;
t283 = t169 * t65 + t243;
t279 = t166 * t111 - t262;
t52 = t169 * t267;
t273 = t166 * t64 - t52;
t276 = pkin(7) * t273 + t290;
t275 = pkin(6) * (t170 * t273 - t287);
t272 = qJ(4) * t260;
t172 = qJD(2) ^ 2;
t173 = qJD(1) ^ 2;
t168 = sin(qJ(1));
t171 = cos(qJ(1));
t192 = t171 * g(1) + t168 * g(2);
t222 = qJDD(1) * pkin(6);
t124 = -t173 * pkin(1) - t192 + t222;
t194 = -t170 * pkin(2) - t167 * pkin(7);
t199 = t173 * t194 + t124;
t247 = t170 * g(3);
t67 = -qJDD(2) * pkin(2) - t172 * pkin(7) + t199 * t167 + t247;
t178 = t89 * pkin(3) - t63 * qJ(4) + t67;
t176 = 0.2e1 * qJD(4) * t137 - t178;
t265 = qJ(4) * t65;
t261 = t254 * pkin(4) - 0.2e1 * qJD(5) * t135;
t246 = pkin(3) + qJ(5);
t53 = qJ(4) * t267;
t258 = -t246 * t64 - t53;
t71 = qJ(4) * t183;
t257 = -t246 * t259 - t71;
t202 = t168 * g(1) - t171 * g(2);
t123 = qJDD(1) * pkin(1) + t173 * pkin(6) + t202;
t190 = -t142 + t155;
t191 = t141 + t204;
t54 = t190 * pkin(2) - t191 * pkin(7) - t123;
t248 = t167 * g(3);
t68 = -t172 * pkin(2) + qJDD(2) * pkin(7) + t199 * t170 - t248;
t31 = t166 * t68 - t169 * t54;
t94 = t135 * pkin(3) - t137 * qJ(4);
t22 = t134 * pkin(3) - t151 * qJ(4) + t137 * t94 + qJDD(4) + t31;
t179 = t90 * pkin(4) + t182 * qJ(5) + t22;
t252 = 0.2e1 * qJD(5);
t11 = (-pkin(4) * t135 + t252) * t153 + t179;
t148 = -0.2e1 * t212;
t107 = t137 * pkin(4) + t153 * qJ(5);
t32 = t166 * t54 + t169 * t68;
t195 = -t151 * pkin(3) - t134 * qJ(4) - t135 * t94 + t32;
t181 = -t89 * pkin(4) - qJ(5) * t254 - t153 * t107 + qJDD(5) + t195;
t13 = t148 + t181;
t256 = qJ(4) * t13 - t11 * t246;
t219 = t153 * t166;
t110 = t137 * t219;
t206 = t170 * t221;
t196 = t167 * (t169 * t90 + t110) - t206;
t255 = t182 * t246 - t272;
t251 = pkin(3) * t153;
t250 = pkin(3) * t166;
t249 = pkin(3) * t169;
t239 = t166 * t67;
t231 = t169 * t67;
t225 = t89 * qJ(5);
t223 = qJ(4) * t169;
t218 = t153 * t169;
t152 = t170 * t173 * t167;
t216 = t167 * (qJDD(2) + t152);
t215 = t170 * (-t152 + qJDD(2));
t207 = t135 * t218;
t205 = -pkin(3) * t61 - t53;
t203 = qJ(4) * t166 + pkin(2);
t18 = t166 * t31 + t169 * t32;
t201 = -0.2e1 * qJD(4) - t251;
t108 = t167 * t124 + t247;
t109 = t170 * t124 - t248;
t200 = t167 * t108 + t170 * t109;
t197 = t167 * (t166 * t89 - t207) + t206;
t48 = -t135 * t219 - t169 * t89;
t49 = -t137 * t218 + t166 * t90;
t20 = t148 + t195;
t193 = -pkin(3) * t22 + qJ(4) * t20;
t189 = t166 * t32 - t169 * t31;
t188 = -pkin(1) + t194;
t186 = (t135 * t166 + t137 * t169) * t153;
t185 = -pkin(3) * t259 + t195 - t71;
t120 = t170 * t134;
t184 = t167 * (-t110 + t207) + t120;
t177 = pkin(3) * t182 + t22 - t272;
t175 = t153 * t252 + t179;
t174 = t176 + t261;
t163 = t170 ^ 2;
t162 = t167 ^ 2;
t160 = t163 * t173;
t158 = t162 * t173;
t143 = -0.2e1 * t155 + t208;
t140 = t156 + 0.2e1 * t204;
t33 = pkin(4) * t182 - qJ(4) * t266;
t27 = t166 * t61 - t52;
t24 = -t169 * t61 - t242;
t23 = -pkin(4) * t183 - t246 * t65;
t21 = t201 * t137 + t178;
t19 = t22 + t288;
t16 = pkin(3) * t74 + t20;
t15 = -t176 + t289;
t14 = pkin(3) * t220 + t176 - t265;
t12 = t225 + (-t107 + t201) * t137 + t178 - t261;
t10 = t166 * t22 + t169 * t20;
t9 = t166 * t20 - t169 * t22;
t8 = (t107 + t251) * t137 - t265 + t174 - t225 + pkin(4) * t259;
t7 = t288 + (t64 - t117) * pkin(4) + t175;
t6 = -pkin(4) * t267 + t246 * t74 + t13;
t5 = (-t266 - t89) * qJ(5) - t289 + t174 + pkin(4) * t260 + t137 * t107;
t4 = t166 * t11 + t169 * t13;
t3 = -t169 * t11 + t166 * t13;
t2 = pkin(4) * t11 - qJ(4) * t12;
t1 = pkin(4) * t13 - t246 * t12;
t17 = [0, 0, 0, 0, 0, qJDD(1), t202, t192, 0, 0, t191 * t167, t170 * t140 + t167 * t143, t216 + t170 * (-t158 + t172), -t190 * t170, t167 * (t160 - t172) + t215, 0, t170 * t123 + pkin(1) * t143 + pkin(6) * (t170 * (-t160 - t172) - t216), -t167 * t123 - pkin(1) * t140 + pkin(6) * (-t215 - t167 * (-t158 - t172)), pkin(1) * (t158 + t160) + (t162 + t163) * t222 + t200, pkin(1) * t123 + pkin(6) * t200, t196, t167 * (-t166 * t63 - t234) + t286, -t170 * t64 - t280, t197, t297, t184, t167 * (t239 + t304) + t170 * (t31 + t306) + t310, t167 * (t231 - t303) + t170 * (t32 - t305) - t309, -t167 * t189 + t188 * t25 + t275, pkin(6) * (t167 * t67 + t170 * t18) + t188 * t189, t120 + t167 * (t135 * t169 - t137 * t166) * t153, t311, -t297, t196, -t296, t197, t167 * (-pkin(7) * t24 - t166 * t16 + t169 * t19) + t170 * (-pkin(2) * t24 - t205) - pkin(1) * t24 + pkin(6) * (t170 * t27 - t287), t167 * (-t166 * t15 + t223 * t266 - t304) + t170 * (-t177 - t306) - t310, t167 * (t169 * t14 + t250 * t65 + t303) + t170 * (-t185 + t312) + t309, t167 * (-pkin(7) * t9 + (-t223 + t250) * t21) + t170 * (-pkin(2) * t9 - t193) - pkin(1) * t9 + pkin(6) * (t170 * t10 + t167 * t21), t184, -t297, -t311, t197, t296, t196, t167 * (-pkin(7) * t25 - t166 * t6 + t169 * t7) + t170 * (-pkin(2) * t25 - t258) - pkin(1) * t25 + t275, t167 * (-t166 * t23 + t169 * t8 + t303) + t170 * (-t181 - t257 + t312) + t309, t167 * (-t166 * t5 + t169 * t33 + t304) + t170 * (t11 + t255 + t306) + t310, t167 * (-pkin(7) * t3 - t166 * t1 + t169 * t2) + t170 * (-pkin(2) * t3 - t256) - pkin(1) * t3 + pkin(6) * (t167 * t12 + t170 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t152, t158 - t160, t156, t152, t208, qJDD(2), -t108, -t109, 0, 0, t49, t169 * t63 - t243, t285, t48, t279, t186, -t231 + t299, t239 + t298, t18 + t276, -pkin(2) * t67 + pkin(7) * t18, t186, -t285, -t279, t49, -t283, t48, pkin(7) * t27 + t169 * t16 + t166 * t19 + t290, t169 * t15 + t203 * t266 - t301, t302 + t166 * t14 - (pkin(2) + t249) * t65, pkin(7) * t10 + (-t203 - t249) * t21, t186, -t279, t285, t48, t283, t49, t166 * t7 + t169 * t6 + t276, t166 * t8 + t169 * t23 - t298, t166 * t33 + t169 * t5 + t299, -pkin(2) * t12 + pkin(7) * t4 + t169 * t1 + t166 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t221, -t98, t64, -t221, -t267, -t134, -t31, -t32, 0, 0, -t134, -t61, t267, t221, -t98, -t221, t205, t177, t148 + t185, t193, -t134, t267, t61, -t221, t98, t221, t258, t13 + t257, pkin(4) * t117 - t175 - t255, t256; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t182, t259, t22, 0, 0, 0, 0, 0, 0, t64, t259, t182, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t267, -t183, t260, t13;];
tauJ_reg = t17;
