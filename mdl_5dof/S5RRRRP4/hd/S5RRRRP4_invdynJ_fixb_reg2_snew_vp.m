% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRRP4
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
% Datum: 2019-12-31 21:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:51:14
% EndTime: 2019-12-31 21:51:22
% DurationCPUTime: 2.91s
% Computational Cost: add. (8917->243), mult. (11444->296), div. (0->0), fcn. (7275->8), ass. (0->167)
t203 = sin(qJ(2));
t206 = cos(qJ(2));
t198 = qJD(1) + qJD(2);
t201 = sin(qJ(4));
t204 = cos(qJ(4));
t205 = cos(qJ(3));
t202 = sin(qJ(3));
t247 = t198 * t202;
t161 = -t204 * t205 * t198 + t201 * t247;
t236 = qJD(3) * t198;
t231 = t205 * t236;
t196 = qJDD(1) + qJDD(2);
t241 = t202 * t196;
t170 = t231 + t241;
t232 = t202 * t236;
t237 = t205 * t196;
t219 = -t232 + t237;
t103 = -t161 * qJD(4) + t204 * t170 + t201 * t219;
t197 = qJD(3) + qJD(4);
t250 = t197 * t161;
t275 = t103 - t250;
t195 = qJDD(3) + qJDD(4);
t240 = t202 * t204;
t163 = (t205 * t201 + t240) * t198;
t252 = t163 * t161;
t122 = -t252 - t195;
t246 = t201 * t122;
t160 = t163 ^ 2;
t268 = t197 ^ 2;
t273 = -t160 - t268;
t92 = t204 * t273 + t246;
t239 = t204 * t122;
t94 = -t201 * t273 + t239;
t40 = t202 * t92 - t205 * t94;
t307 = pkin(1) * (t203 * t40 + t206 * t275);
t306 = -pkin(2) * t275 - pkin(7) * t40;
t304 = pkin(3) * t92;
t303 = pkin(8) * t92;
t302 = pkin(8) * t94;
t269 = t161 ^ 2;
t148 = t269 - t268;
t61 = t202 * (-t204 * t148 - t246) + t205 * (-t201 * t148 + t239);
t156 = t197 * t163;
t223 = -t201 * t170 + t204 * t219;
t218 = t163 * qJD(4) - t223;
t274 = t156 + t218;
t271 = -t252 + t195;
t245 = t201 * t271;
t270 = -t268 - t269;
t277 = t204 * t270 - t245;
t110 = t204 * t271;
t278 = t201 * t270 + t110;
t287 = -t202 * t278 + t205 * t277;
t299 = -pkin(2) * t274 + pkin(7) * t287;
t259 = t201 * t275;
t24 = t202 * (t204 * t274 + t259) - t205 * (-t201 * t274 + t204 * t275);
t298 = pkin(1) * (t203 * t287 - t206 * t274);
t104 = -t269 - t160;
t297 = pkin(2) * t104;
t296 = pkin(3) * t104;
t295 = pkin(3) * t278;
t294 = pkin(8) * t277;
t293 = pkin(8) * t278;
t289 = t206 * t104;
t149 = -t160 + t268;
t286 = t202 * (-t201 * t149 + t110) + t205 * (t204 * t149 + t245);
t283 = qJ(5) * t275;
t265 = sin(qJ(1));
t266 = cos(qJ(1));
t217 = t266 * g(1) + t265 * g(2);
t176 = -qJD(1) ^ 2 * pkin(1) - t217;
t216 = t265 * g(1) - t266 * g(2);
t214 = qJDD(1) * pkin(1) + t216;
t134 = t206 * t176 + t203 * t214;
t194 = t198 ^ 2;
t126 = -t194 * pkin(2) + t196 * pkin(7) + t134;
t243 = t202 * t126;
t108 = t205 * g(3) + t243;
t109 = -t202 * g(3) + t205 * t126;
t65 = t202 * t108 + t205 * t109;
t127 = t160 - t269;
t272 = t250 + t103;
t267 = 2 * qJD(5);
t264 = pkin(4) * t204;
t124 = t161 * pkin(4) - t163 * qJ(5);
t251 = t194 * t202;
t208 = qJDD(3) * pkin(3) - t170 * pkin(8) - t243 + (pkin(3) * t251 + pkin(8) * t236 - g(3)) * t205;
t179 = qJD(3) * pkin(3) - pkin(8) * t247;
t200 = t205 ^ 2;
t186 = t200 * t194;
t73 = -pkin(3) * t186 + t219 * pkin(8) - qJD(3) * t179 + t109;
t48 = t201 * t208 + t204 * t73;
t221 = t195 * qJ(5) - t161 * t124 + t197 * t267 + t48;
t34 = -pkin(4) * t268 + t221;
t47 = t201 * t73 - t204 * t208;
t36 = -t195 * pkin(4) - qJ(5) * t268 + t163 * t124 + qJDD(5) + t47;
t263 = -pkin(4) * t36 + qJ(5) * t34;
t85 = -t156 + t218;
t262 = -pkin(4) * t272 - qJ(5) * t85;
t133 = -t203 * t176 + t206 * t214;
t125 = -t196 * pkin(2) - t194 * pkin(7) - t133;
t82 = -t219 * pkin(3) - pkin(8) * t186 + t179 * t247 + t125;
t261 = t201 * t82;
t258 = t201 * t272;
t18 = t201 * t48 - t204 * t47;
t257 = t202 * t18;
t256 = t204 * t82;
t254 = -pkin(2) * t125 + pkin(7) * t65;
t253 = qJ(5) * t204;
t249 = t197 * t201;
t248 = t197 * t204;
t181 = t205 * t251;
t242 = t202 * (qJDD(3) + t181);
t238 = t205 * (qJDD(3) - t181);
t199 = t202 ^ 2;
t185 = t199 * t194;
t207 = qJD(3) ^ 2;
t143 = -t238 - t202 * (-t185 - t207);
t169 = 0.2e1 * t231 + t241;
t235 = -pkin(2) * t169 + pkin(7) * t143 + t202 * t125;
t142 = t205 * (-t186 - t207) - t242;
t171 = -0.2e1 * t232 + t237;
t234 = pkin(2) * t171 + pkin(7) * t142 - t205 * t125;
t15 = t201 * t34 - t204 * t36;
t16 = t201 * t36 + t204 * t34;
t228 = -qJ(5) * t201 - pkin(3);
t212 = t218 * pkin(4) - t283 + t82;
t31 = (pkin(4) * t197 - (2 * qJD(5))) * t163 + t212;
t4 = -t202 * t15 + t205 * t16;
t233 = t202 * (-pkin(8) * t15 + (pkin(4) * t201 - t253) * t31) + t205 * (pkin(8) * t16 + (t228 - t264) * t31) - pkin(2) * t31 + pkin(7) * t4;
t210 = t163 * t267 - t212;
t20 = -pkin(4) * t156 + t210 + t283;
t230 = t202 * (-pkin(4) * t259 + t204 * t20 + t303) + t205 * (-t302 + t201 * t20 + (pkin(3) + t264) * t275) - t306;
t74 = t204 * t272;
t49 = -t201 * t85 - t74;
t51 = -t204 * t85 + t258;
t26 = -t202 * t49 + t205 * t51;
t28 = (-t104 - t268) * pkin(4) + t221;
t29 = -qJ(5) * t104 + t36;
t229 = t202 * (-pkin(8) * t49 - t201 * t28 + t204 * t29) + t205 * (pkin(8) * t51 + t201 * t29 + t204 * t28 - t296) - t297 + pkin(7) * t26;
t19 = t201 * t47 + t204 * t48;
t21 = (-t274 - t156) * pkin(4) + t210;
t227 = t202 * (-t201 * t21 - t253 * t274 - t293) + t205 * (t204 * t21 + t228 * t274 + t294) + t299;
t226 = t202 * (t261 - t293) + t205 * (-pkin(3) * t274 - t256 + t294) + t299;
t225 = t202 * (t256 - t303) + t205 * (-pkin(3) * t275 + t261 + t302) + t306;
t86 = (-qJD(4) + t197) * t163 + t223;
t50 = t201 * t86 - t74;
t52 = t204 * t86 + t258;
t27 = -t202 * t50 + t205 * t52;
t224 = t202 * (-pkin(8) * t50 - t18) + t205 * (pkin(8) * t52 + t19 - t296) - t297 + pkin(7) * t27;
t174 = (t199 + t200) * t196;
t175 = t185 + t186;
t222 = pkin(2) * t175 + pkin(7) * t174 + t65;
t8 = t205 * t19 - t257;
t220 = pkin(7) * t8 - pkin(8) * t257 + t205 * (-pkin(3) * t82 + pkin(8) * t19) - pkin(2) * t82;
t215 = -pkin(4) * t273 - qJ(5) * t122 + t34;
t213 = pkin(4) * t271 + qJ(5) * t270 - t36;
t211 = t202 * (t161 * t248 + t201 * t218) + t205 * (t161 * t249 - t204 * t218);
t147 = t163 * t249;
t209 = t202 * t147 + (-t161 * t240 + t205 * (-t161 * t201 - t163 * t204)) * t197;
t141 = t242 + t205 * (-t185 + t207);
t140 = t202 * (t186 - t207) + t238;
t136 = (t170 + t231) * t202;
t135 = t171 * t205;
t132 = t205 * t169 + t202 * t171;
t44 = t202 * (t204 * t103 - t147) + t205 * (t201 * t103 + t163 * t248);
t1 = [0, 0, 0, 0, 0, qJDD(1), t216, t217, 0, 0, 0, 0, 0, 0, 0, t196, pkin(1) * (-t203 * t194 + t206 * t196) + t133, pkin(1) * (-t206 * t194 - t203 * t196) - t134, 0, pkin(1) * (t206 * t133 + t203 * t134), t136, t132, t141, t135, t140, 0, pkin(1) * (t203 * t142 + t206 * t171) + t234, pkin(1) * (t203 * t143 - t206 * t169) + t235, pkin(1) * (t203 * t174 + t206 * t175) + t222, pkin(1) * (-t206 * t125 + t203 * t65) + t254, t44, -t24, t286, t211, -t61, t209, t298 + t226, t225 - t307, pkin(1) * (t203 * t27 - t289) + t224, pkin(1) * (t203 * t8 - t206 * t82) + t220, t44, t286, t24, t209, t61, t211, t298 + t227, pkin(1) * (t203 * t26 - t289) + t229, t230 + t307, pkin(1) * (t203 * t4 - t206 * t31) + t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196, t133, -t134, 0, 0, t136, t132, t141, t135, t140, 0, t234, t235, t222, t254, t44, -t24, t286, t211, -t61, t209, t226, t225, t224, t220, t44, t286, t24, t209, t61, t211, t227, t229, t230, t233; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t181, t185 - t186, t241, t181, t237, qJDD(3), -t108, -t109, 0, 0, t252, t127, t272, -t252, -t85, t195, -t47 + t295, -t48 + t304, pkin(3) * t50, pkin(3) * t18, t252, t272, -t127, t195, t85, -t252, t213 + t295, pkin(3) * t49 + t262, t215 - t304, pkin(3) * t15 + t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t252, t127, t272, -t252, -t85, t195, -t47, -t48, 0, 0, t252, t272, -t127, t195, t85, -t252, t213, t262, t215, t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t271, t272, t273, t36;];
tauJ_reg = t1;
