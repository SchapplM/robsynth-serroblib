% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPRP4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 17:45
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPRP4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:44:27
% EndTime: 2019-05-05 17:44:34
% DurationCPUTime: 2.67s
% Computational Cost: add. (5992->298), mult. (11986->337), div. (0->0), fcn. (6784->8), ass. (0->191)
t182 = cos(qJ(5));
t179 = sin(qJ(5));
t183 = cos(qJ(3));
t232 = qJD(1) * qJD(3);
t163 = t183 * t232;
t180 = sin(qJ(3));
t229 = t180 * qJDD(1);
t139 = t163 + t229;
t126 = qJDD(5) + t139;
t235 = qJD(1) * t183;
t131 = t179 * qJD(3) + t182 * t235;
t133 = t182 * qJD(3) - t179 * t235;
t249 = t133 * t131;
t277 = t126 + t249;
t261 = t179 * t277;
t125 = t133 ^ 2;
t234 = t180 * qJD(1);
t158 = qJD(5) + t234;
t271 = t158 ^ 2;
t279 = -t125 - t271;
t39 = -t182 * t279 + t261;
t317 = pkin(4) * t39;
t316 = t180 * t39;
t315 = t183 * t39;
t268 = -pkin(3) - pkin(8);
t314 = t268 * t39;
t250 = t131 * t158;
t221 = t180 * t232;
t228 = t183 * qJDD(1);
t140 = -t221 + t228;
t202 = t182 * qJDD(3) - t179 * t140;
t84 = -t131 * qJD(5) + t202;
t301 = t250 - t84;
t313 = qJ(6) * t301;
t272 = t131 ^ 2;
t106 = t272 - t271;
t254 = t182 * t277;
t112 = t158 * t133;
t215 = t179 * qJDD(3) + t182 * t140;
t199 = t133 * qJD(5) + t215;
t55 = -t112 + t199;
t312 = t180 * t55 + t183 * (t179 * t106 + t254);
t176 = cos(pkin(9));
t245 = t180 * qJ(4);
t201 = -pkin(1) * t176 - pkin(2) - t245;
t193 = t183 * t268 + t201;
t311 = t193 * (t179 * t279 + t254);
t309 = -t182 * t106 + t261;
t256 = t182 * t301;
t54 = t112 + t199;
t93 = t125 - t272;
t308 = t183 * (-t179 * t54 - t256) - t180 * t93;
t270 = 2 * qJD(4);
t70 = -t272 - t125;
t307 = pkin(4) * t70;
t278 = t126 - t249;
t253 = t182 * t278;
t274 = -t271 - t272;
t285 = t179 * t274 + t253;
t306 = pkin(4) * t285;
t305 = qJ(4) * t70;
t304 = t180 * t70;
t302 = t183 * t70;
t299 = t180 * t285;
t297 = t183 * t285;
t296 = t268 * t285;
t108 = -t125 + t271;
t295 = -t179 * t108 + t253;
t280 = t250 + t84;
t67 = t179 * t278;
t292 = t180 * t280 + t183 * (-t182 * t108 - t67);
t287 = t139 + t163;
t186 = qJD(1) ^ 2;
t265 = t183 * pkin(3);
t181 = sin(qJ(1));
t184 = cos(qJ(1));
t220 = t181 * g(1) - t184 * g(2);
t134 = qJDD(1) * pkin(1) + t220;
t211 = t184 * g(1) + t181 * g(2);
t136 = -t186 * pkin(1) - t211;
t175 = sin(pkin(9));
t237 = t175 * t134 + t176 * t136;
t78 = -t186 * pkin(2) + qJDD(1) * pkin(7) + t237;
t217 = t186 * (-t245 - t265) + t78;
t207 = t217 * t183;
t286 = qJD(3) * t270 + t207;
t284 = -t179 * t301 + t182 * t54;
t283 = t193 * (t182 * t274 - t67);
t282 = qJ(6) * t179 + pkin(4);
t171 = t183 ^ 2;
t166 = t171 * t186;
t185 = qJD(3) ^ 2;
t238 = t183 * t186;
t223 = t180 * t238;
t147 = qJDD(3) - t223;
t239 = t183 * t147;
t276 = t239 + t180 * (t166 - t185);
t170 = t180 ^ 2;
t165 = t170 * t186;
t152 = -t165 - t185;
t275 = t180 * t152 + t239;
t273 = pkin(5) * t199 + t313;
t269 = 2 * qJD(6);
t149 = pkin(4) * t234 - qJD(3) * pkin(8);
t216 = t176 * t134 - t175 * t136;
t77 = -qJDD(1) * pkin(2) - t186 * pkin(7) - t216;
t195 = -t140 * pkin(3) - t287 * qJ(4) + t77;
t219 = pkin(3) * qJD(3) - (2 * qJD(4));
t32 = -pkin(4) * t166 - t140 * pkin(8) + (-t149 + t219) * t234 + t195;
t172 = -g(3) + qJDD(2);
t162 = t183 * t172;
t214 = -qJDD(3) * pkin(3) - t185 * qJ(4) + qJDD(4) - t162;
t36 = -qJDD(3) * pkin(8) + (t139 - t163) * pkin(4) + (-pkin(8) * t238 + t217) * t180 + t214;
t19 = t179 * t36 + t182 * t32;
t91 = t131 * pkin(5) - t133 * qJ(6);
t212 = t126 * qJ(6) - t131 * t91 + t158 * t269 + t19;
t10 = -pkin(5) * t271 + t212;
t18 = t179 * t32 - t182 * t36;
t13 = -t126 * pkin(5) - qJ(6) * t271 + t133 * t91 + qJDD(6) + t18;
t267 = -pkin(5) * t13 + qJ(6) * t10;
t266 = pkin(5) * t179;
t264 = -pkin(5) * t280 - qJ(6) * t55;
t177 = t185 * pkin(3);
t230 = qJDD(3) * qJ(4);
t241 = t180 * t172;
t196 = t140 * pkin(4) - pkin(8) * t166 - t177 + t230 + t241;
t35 = t207 + (t270 + t149) * qJD(3) + t196;
t263 = t179 * t35;
t255 = t182 * t280;
t57 = (-qJD(5) + t158) * t133 - t215;
t29 = t179 * t57 - t255;
t260 = t180 * t29;
t258 = t182 * t35;
t16 = (pkin(5) * t158 - (2 * qJD(6))) * t133 + t35 + t273;
t252 = t183 * t16;
t248 = t158 * t179;
t247 = t158 * t182;
t146 = qJDD(3) + t223;
t244 = t180 * t146;
t233 = qJD(5) + t158;
t143 = (t170 + t171) * qJDD(1);
t144 = t165 + t166;
t227 = pkin(7) * t143 + pkin(2) * t144 + pkin(1) * (t175 * t143 + t176 * t144);
t226 = t131 * t248;
t225 = t180 * t249;
t224 = t131 * t247;
t222 = pkin(1) * t175 + pkin(7);
t64 = t180 * t78 - t162;
t65 = t183 * t78 + t241;
t37 = t180 * t64 + t183 * t65;
t218 = -qJ(6) * t182 + qJ(4);
t105 = t133 * t247;
t213 = t183 * (-t179 * t84 - t105) + t225;
t104 = t133 * t248;
t210 = t104 - t224;
t5 = t179 * t19 - t182 * t18;
t208 = t179 * t18 + t182 * t19;
t206 = t183 * (-t165 + t185) + t244;
t154 = -t166 - t185;
t205 = t183 * t146 + t180 * t154;
t204 = -t183 * t154 + t244;
t203 = t180 * t147 - t183 * t152;
t200 = t179 * t199 + t224;
t198 = t180 * t126 + t183 * (t105 + t226);
t197 = -t201 + t265;
t194 = -pkin(5) * t279 + qJ(6) * t277 + t10;
t192 = -t177 + t286;
t48 = t217 * t180 + t214;
t191 = pkin(5) * t278 + qJ(6) * t274 - t13;
t190 = -t225 + t183 * (t182 * t199 - t226);
t189 = t192 + t230;
t188 = t219 * t234 + t195;
t187 = -qJD(3) * t149 + t133 * t269 - t196 - t273 - t286;
t145 = t165 - t166;
t141 = -0.2e1 * t221 + t228;
t138 = 0.2e1 * t163 + t229;
t99 = t287 * t180;
t98 = (t140 - t221) * t183;
t88 = t183 * t138 + t180 * t141;
t60 = -t233 * t131 + t202;
t56 = t233 * t133 + t215;
t51 = t179 * t280;
t50 = t182 * t84 - t104;
t43 = t189 + t241;
t28 = -t179 * t55 - t255;
t15 = t187 + (-t56 - t112) * pkin(5);
t14 = -pkin(5) * t112 + t187 - t313;
t8 = -qJ(6) * t70 + t13;
t7 = (-t271 - t70) * pkin(5) + t212;
t2 = t179 * t10 - t182 * t13;
t1 = [0, 0, 0, 0, 0, qJDD(1), t220, t211, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t176 * qJDD(1) - t175 * t186) + t216, pkin(1) * (-t175 * qJDD(1) - t176 * t186) - t237, 0, pkin(1) * (t175 * t237 + t176 * t216), t99, t88, t206, t98, t276, 0, -t183 * t77 + pkin(2) * t141 - pkin(7) * t204 + pkin(1) * (t176 * t141 - t175 * t204), t180 * t77 - pkin(2) * t138 - pkin(7) * t275 + pkin(1) * (-t176 * t138 - t175 * t275), t37 + t227, -pkin(2) * t77 + pkin(7) * t37 + pkin(1) * (t175 * t37 - t176 * t77), 0, -t206, -t276, t99, t88, t98, (pkin(3) * t144 + t189) * t183 + (qJ(4) * t144 + t162 + t48) * t180 + t227, -t141 * t197 + t183 * t188 + t204 * t222, t180 * (-pkin(3) * t221 + t234 * t270 - t195) + t222 * t275 + t197 * t138, t222 * (t180 * t48 + t183 * t43) - t197 * t188, t213, -t308, t292, t190, -t312, t198, t180 * (-t18 + t306) + t183 * (pkin(4) * t54 + t258) + t222 * (t183 * t54 + t299) + t283, t180 * (-t19 - t317) + t183 * (pkin(4) * t60 - t263) + t222 * (t183 * t60 - t316) - t311, pkin(4) * t260 + t183 * (-t208 + t307) + t222 * (t260 + t302) + t193 * (t182 * t57 + t51), t193 * t208 + (pkin(4) + t222) * (t180 * t5 + t183 * t35), t213, t292, t308, t198, t312, t190, t180 * (t191 + t306) + t183 * (-t182 * t15 + t282 * t56) + t222 * (t183 * t56 + t299) + t283, t180 * (pkin(4) * t28 + t264) + t183 * (-t179 * t8 - t182 * t7 + t307) + t222 * (t180 * t28 + t302) + t193 * (-t182 * t55 + t51), t180 * (t194 + t317) + t183 * (pkin(4) * t301 + pkin(5) * t256 - t179 * t14) + t222 * (t183 * t301 + t316) + t311, t180 * (pkin(4) * t2 + t267) + t222 * (t180 * t2 + t252) + t193 * (t182 * t10 + t179 * t13) + (pkin(5) * t182 + t282) * t252; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t172, 0, 0, 0, 0, 0, 0, t205, -t203, 0, t180 * t65 - t183 * t64, 0, 0, 0, 0, 0, 0, 0, -t205, t203, t180 * t43 - t183 * t48, 0, 0, 0, 0, 0, 0, t180 * t54 - t297, t180 * t60 + t315, -t183 * t29 + t304, t180 * t35 - t183 * t5, 0, 0, 0, 0, 0, 0, t180 * t56 - t297, -t183 * t28 + t304, t180 * t301 - t315, t180 * t16 - t183 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t223, t145, t229, t223, t228, qJDD(3), -t64, -t65, 0, 0, qJDD(3), -t229, -t228, -t223, t145, t223, (-pkin(3) * t180 + qJ(4) * t183) * qJDD(1), -pkin(3) * t146 - qJ(4) * t154 + t48, -pkin(3) * t152 + t241 + (qJDD(3) + t147) * qJ(4) + t192, -pkin(3) * t48 + qJ(4) * t43, t50, -t284, t295, t200, -t309, t210, qJ(4) * t54 + t263 + t296, qJ(4) * t60 + t258 - t314, t268 * t29 + t305 - t5, qJ(4) * t35 + t268 * t5, t50, t295, t284, t210, t309, t200, -t179 * t15 + t218 * t56 + t296, -t179 * t7 + t182 * t8 + t268 * t28 + t305, t182 * t14 - (-qJ(4) - t266) * t301 + t314, t268 * t2 + (t218 + t266) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t229, t146, t152, t48, 0, 0, 0, 0, 0, 0, t285, -t39, t29, t5, 0, 0, 0, 0, 0, 0, t285, t28, t39, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t249, t93, t280, -t249, -t55, t126, -t18, -t19, 0, 0, t249, t280, -t93, t126, t55, -t249, t191, t264, t194, t267; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t278, t280, t279, t13;];
tauJ_reg  = t1;
