% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6RPRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 17:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6RPRPRP7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:58:30
% EndTime: 2019-05-05 17:58:41
% DurationCPUTime: 4.14s
% Computational Cost: add. (13878->361), mult. (30741->470), div. (0->0), fcn. (20721->8), ass. (0->243)
t202 = sin(pkin(9));
t203 = cos(pkin(9));
t206 = sin(qJ(3));
t209 = cos(qJ(3));
t181 = (t202 * t209 + t203 * t206) * qJD(1);
t262 = qJD(1) * t209;
t183 = -qJD(1) * t202 * t206 + t203 * t262;
t272 = t183 * t181;
t308 = qJDD(3) - t272;
t310 = t202 * t308;
t309 = t203 * t308;
t252 = qJD(1) * qJD(3);
t238 = t209 * t252;
t250 = t206 * qJDD(1);
t187 = -t238 - t250;
t195 = t209 * qJDD(1);
t239 = t206 * t252;
t188 = t195 - t239;
t154 = t187 * t202 + t188 * t203;
t261 = qJD(3) * t181;
t137 = t154 - t261;
t300 = pkin(7) + pkin(1);
t205 = sin(qJ(5));
t208 = cos(qJ(5));
t163 = -qJD(3) * t208 + t183 * t205;
t165 = qJD(3) * t205 + t183 * t208;
t134 = t165 * t163;
t235 = -t187 * t203 + t188 * t202;
t151 = qJDD(5) + t235;
t302 = -t134 + t151;
t307 = pkin(5) * t302;
t121 = -qJD(5) * t163 + qJDD(3) * t205 + t154 * t208;
t177 = qJD(5) + t181;
t143 = t177 * t163;
t103 = t121 + t143;
t306 = qJ(6) * t103;
t280 = t302 * t205;
t279 = t302 * t208;
t162 = t165 ^ 2;
t176 = t177 ^ 2;
t128 = -t162 - t176;
t161 = t163 ^ 2;
t236 = -qJDD(3) * t208 + t205 * t154;
t120 = -qJD(5) * t165 - t236;
t139 = pkin(5) * t177 - qJ(6) * t165;
t145 = pkin(4) * t181 - pkin(8) * t183;
t211 = qJD(3) ^ 2;
t207 = sin(qJ(1));
t210 = cos(qJ(1));
t228 = t207 * g(1) - t210 * g(2);
t223 = qJDD(2) - t228;
t212 = qJD(1) ^ 2;
t266 = t212 * qJ(2);
t219 = t223 - t266;
t215 = -qJDD(1) * t300 + t219;
t156 = t209 * g(3) - t206 * t215;
t222 = qJD(3) * pkin(3) - qJ(4) * t262;
t199 = t206 ^ 2;
t271 = t199 * t212;
t129 = -pkin(3) * t271 + t187 * qJ(4) - qJD(3) * t222 - t156;
t214 = t209 * t215;
t267 = t209 * t212;
t213 = t214 - t188 * qJ(4) + qJDD(3) * pkin(3) + (-pkin(3) * t267 - qJ(4) * t252 + g(3)) * t206;
t87 = -0.2e1 * qJD(4) * t181 + t129 * t203 + t202 * t213;
t71 = -pkin(4) * t211 + qJDD(3) * pkin(8) - t145 * t181 + t87;
t260 = qJD(3) * t183;
t135 = t235 + t260;
t251 = qJD(2) * qJD(1);
t197 = 0.2e1 * t251;
t198 = qJDD(1) * qJ(2);
t229 = t210 * g(1) + t207 * g(2);
t225 = -t198 + t229;
t301 = -t187 * pkin(3) - (qJ(4) * t199 + t300) * t212 + t222 * t262 + qJDD(4) - t225;
t78 = pkin(4) * t135 - pkin(8) * t137 + t197 + t301;
t41 = t205 * t78 + t208 * t71;
t221 = qJ(6) * t120 - 0.2e1 * qJD(6) * t163 - t177 * t139 + t41;
t305 = -t221 + (t128 + t161) * pkin(5);
t303 = t121 - t143;
t100 = (qJD(5) - t177) * t165 + t236;
t179 = t181 ^ 2;
t180 = t183 ^ 2;
t116 = -t176 - t161;
t79 = t116 * t205 + t279;
t299 = pkin(4) * t79;
t110 = t134 + t151;
t282 = t110 * t205;
t83 = t128 * t208 - t282;
t298 = pkin(4) * t83;
t254 = qJD(6) * t165;
t158 = -0.2e1 * t254;
t40 = t205 * t71 - t208 * t78;
t220 = -t306 - t40 + t307;
t25 = t158 + t220;
t297 = pkin(5) * t25;
t65 = -t100 * t205 - t103 * t208;
t296 = pkin(8) * t65;
t295 = pkin(8) * t79;
t294 = pkin(8) * t83;
t293 = pkin(4) * t202;
t292 = pkin(5) * t103;
t115 = -t161 - t162;
t67 = -t100 * t208 + t103 * t205;
t49 = -t115 * t203 + t202 * t67;
t291 = qJ(4) * t49;
t80 = t116 * t208 - t280;
t99 = (qJD(5) + t177) * t165 + t236;
t55 = t202 * t80 - t203 * t99;
t290 = qJ(4) * t55;
t281 = t110 * t208;
t84 = -t128 * t205 - t281;
t60 = t202 * t84 - t203 * t303;
t289 = qJ(4) * t60;
t288 = t205 * t25;
t237 = t202 * t129 - t203 * t213;
t224 = -qJDD(3) * pkin(4) - pkin(8) * t211 + t237;
t234 = (0.2e1 * qJD(4) + t145) * t183;
t70 = t234 + t224;
t287 = t205 * t70;
t286 = t208 * t25;
t285 = t208 * t70;
t256 = qJD(4) * t183;
t86 = t237 + 0.2e1 * t256;
t51 = t202 * t87 - t203 * t86;
t284 = t209 * t51;
t283 = qJDD(1) * pkin(1);
t246 = -0.2e1 * t251;
t130 = t246 - t301;
t278 = t130 * t202;
t277 = t130 * t203;
t148 = qJDD(3) + t272;
t276 = t148 * t202;
t275 = t148 * t203;
t274 = t177 * t205;
t273 = t177 * t208;
t200 = t209 ^ 2;
t270 = t200 * t212;
t243 = t206 * t267;
t269 = t206 * (qJDD(3) + t243);
t268 = t209 * (qJDD(3) - t243);
t263 = t199 + t200;
t259 = qJD(3) * t202;
t258 = qJD(3) * t203;
t249 = pkin(3) * t55 - pkin(4) * t99 + pkin(8) * t80;
t248 = pkin(3) * t60 - pkin(4) * t303 + pkin(8) * t84;
t247 = pkin(3) * t49 - pkin(4) * t115 + pkin(8) * t67;
t245 = t202 * t134;
t244 = t203 * t134;
t242 = -pkin(4) * t203 - pkin(3);
t56 = t202 * t99 + t203 * t80;
t241 = -pkin(3) * t79 + qJ(4) * t56;
t61 = t202 * t303 + t203 * t84;
t240 = -pkin(3) * t83 + qJ(4) * t61;
t52 = t202 * t86 + t203 * t87;
t19 = t205 * t40 + t208 * t41;
t50 = t115 * t202 + t203 * t67;
t23 = t206 * t50 + t209 * t49;
t233 = qJ(2) * t65 - t23 * t300;
t28 = t206 * t56 + t209 * t55;
t232 = qJ(2) * t79 - t28 * t300;
t32 = t206 * t61 + t209 * t60;
t231 = qJ(2) * t83 - t300 * t32;
t11 = t19 * t202 - t203 * t70;
t3 = t11 * t209 + (t19 * t203 + t202 * t70) * t206;
t18 = t205 * t41 - t208 * t40;
t155 = t206 * g(3) + t214;
t124 = t209 * t155 - t206 * t156;
t136 = -t235 + t260;
t218 = t220 + t307;
t216 = -pkin(5) * t120 - qJ(6) * t161 + t165 * t139 + qJDD(6) + t224;
t39 = t234 + t216;
t190 = t263 * qJDD(1);
t189 = t195 - 0.2e1 * t239;
t186 = 0.2e1 * t238 + t250;
t178 = -t219 + t283;
t174 = -0.2e1 * t256;
t172 = -t180 - t211;
t171 = -t180 + t211;
t170 = t179 - t211;
t169 = t212 * t300 + t225 + t246;
t167 = -t269 + t209 * (-t211 - t270);
t166 = t206 * (-t211 - t271) + t268;
t159 = 0.2e1 * t254;
t146 = -t211 - t179;
t141 = -t162 + t176;
t140 = t161 - t176;
t138 = t154 + t261;
t133 = -t179 - t180;
t131 = t162 - t161;
t127 = -t172 * t202 - t275;
t126 = t172 * t203 - t276;
t113 = t146 * t203 - t310;
t112 = t146 * t202 + t309;
t108 = (-t163 * t208 + t165 * t205) * t177;
t107 = (-t163 * t205 - t165 * t208) * t177;
t106 = t136 * t203 + t138 * t202;
t105 = t136 * t202 - t138 * t203;
t96 = t121 * t208 - t165 * t274;
t95 = t121 * t205 + t165 * t273;
t94 = -t120 * t205 + t163 * t273;
t93 = t120 * t208 + t163 * t274;
t92 = t126 * t209 + t127 * t206;
t91 = t140 * t208 - t282;
t90 = -t141 * t205 + t279;
t89 = t140 * t205 + t281;
t88 = t141 * t208 + t280;
t73 = t112 * t209 + t113 * t206;
t72 = -pkin(5) * t303 - qJ(6) * t110;
t68 = t105 * t209 + t106 * t206;
t66 = -t205 * t303 - t208 * t99;
t64 = -t205 * t99 + t208 * t303;
t57 = t209 * (t108 * t203 + t151 * t202) - t206 * (t108 * t202 - t151 * t203);
t47 = qJ(4) * t50;
t46 = -pkin(4) * t65 + t292;
t45 = t209 * (t203 * t96 + t245) - t206 * (t202 * t96 - t244);
t44 = t209 * (t203 * t94 - t245) - t206 * (t202 * t94 + t244);
t43 = t285 - t294;
t42 = t287 - t295;
t37 = -qJ(6) * t128 + t39;
t36 = t209 * (-t100 * t202 + t203 * t91) - t206 * (t100 * t203 + t202 * t91);
t35 = t209 * (t103 * t202 + t203 * t90) - t206 * (-t103 * t203 + t202 * t90);
t34 = t41 - t298;
t33 = t40 - t299;
t30 = -pkin(5) * t99 + qJ(6) * t116 - t145 * t183 + t174 - t216;
t29 = -pkin(5) * t161 + t221;
t26 = t209 * (t131 * t202 + t203 * t66) - t206 * (-t131 * t203 + t202 * t66);
t24 = t206 * t52 + t284;
t21 = t159 - t220 + t306;
t20 = -qJ(6) * t100 + (-t115 - t161) * pkin(5) + t221;
t17 = -t298 - t305;
t16 = -t205 * t72 + t208 * t37 - t294;
t15 = -qJ(6) * t279 - t205 * t30 - t295;
t14 = t159 - t218 - t299;
t13 = -pkin(5) * t39 + qJ(6) * t29;
t10 = -t18 - t296;
t9 = t208 * t29 - t288;
t8 = t205 * t29 + t286;
t7 = t202 * t39 + t203 * t9;
t6 = t202 * t9 - t203 * t39;
t5 = -t20 * t205 + t208 * t21 - t296;
t4 = -pkin(4) * t8 - t297;
t2 = -pkin(8) * t8 - qJ(6) * t286 - t13 * t205;
t1 = t206 * t7 + t209 * t6;
t12 = [0, 0, 0, 0, 0, qJDD(1), t228, t229, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t223 - 0.2e1 * t283, t197 + 0.2e1 * t198 - t229, pkin(1) * t178 + qJ(2) * (-t212 * pkin(1) + t197 - t225), (t188 - t239) * t209, -t186 * t209 - t189 * t206, t268 - t206 * (t211 - t270), (-t187 + t238) * t206, t209 * (-t211 + t271) - t269, 0, qJ(2) * t186 - t166 * t300 - t206 * t169, qJ(2) * t189 - t167 * t300 - t209 * t169, t190 * t300 - t263 * t266 - t124, -qJ(2) * t169 - t124 * t300, t209 * (t154 * t203 - t183 * t259) - t206 * (t154 * t202 + t183 * t258), t209 * (-t135 * t203 - t137 * t202) - t206 * (-t135 * t202 + t137 * t203), t209 * (-t171 * t202 + t309) - t206 * (t171 * t203 + t310), t209 * (t181 * t258 + t202 * t235) - t206 * (t181 * t259 - t203 * t235), t209 * (t170 * t203 - t276) - t206 * (t170 * t202 + t275), (t209 * (-t181 * t203 + t183 * t202) - t206 * (-t181 * t202 - t183 * t203)) * qJD(3), t209 * (-qJ(4) * t112 - t278) - t206 * (-pkin(3) * t135 + qJ(4) * t113 + t277) + qJ(2) * t135 - t300 * t73, t209 * (-qJ(4) * t126 - t277) - t206 * (-pkin(3) * t137 + qJ(4) * t127 - t278) + qJ(2) * t137 - t300 * t92, t209 * (-qJ(4) * t105 - t51) - t206 * (-pkin(3) * t133 + qJ(4) * t106 + t52) + qJ(2) * t133 - t300 * t68, -qJ(4) * t284 - t206 * (pkin(3) * t130 + qJ(4) * t52) - qJ(2) * t130 - t300 * t24, t45, t26, t35, t44, t36, t57, t209 * (-t202 * t33 + t203 * t42 - t290) - t206 * (t202 * t42 + t203 * t33 + t241) + t232, t209 * (-t202 * t34 + t203 * t43 - t289) - t206 * (t202 * t43 + t203 * t34 + t240) + t231, t209 * (t10 * t203 + t293 * t65 - t291) - t206 * (t202 * t10 + t242 * t65 + t47) + t233, (t209 * (-pkin(8) * t203 + t293) - t206 * (-pkin(8) * t202 + t242) + qJ(2)) * t18 + (-t300 - qJ(4)) * t3, t45, t26, t35, t44, t36, t57, t209 * (-t14 * t202 + t15 * t203 - t290) - t206 * (t14 * t203 + t15 * t202 + t241) + t232, t209 * (t16 * t203 - t17 * t202 - t289) - t206 * (t16 * t202 + t17 * t203 + t240) + t231, t209 * (-t202 * t46 + t203 * t5 - t291) - t206 * (-pkin(3) * t65 + t202 * t5 + t203 * t46 + t47) + t233, t209 * (-qJ(4) * t6 + t2 * t203 - t202 * t4) - t206 * (-pkin(3) * t8 + qJ(4) * t7 + t2 * t202 + t203 * t4) + qJ(2) * t8 - t300 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t212, -t178, 0, 0, 0, 0, 0, 0, t166, t167, -t190, t124, 0, 0, 0, 0, 0, 0, t73, t92, t68, t24, 0, 0, 0, 0, 0, 0, t28, t32, t23, t3, 0, 0, 0, 0, 0, 0, t28, t32, t23, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t243, (-t199 + t200) * t212, t195, -t243, -t250, qJDD(3), t155, t156, 0, 0, t272, t180 - t179, t138, -t272, t136, qJDD(3), pkin(3) * t112 + t174 - t237, pkin(3) * t126 - t87, pkin(3) * t105, pkin(3) * t51, t95, t64, t88, t93, t89, t107, t249 - t285, t248 + t287, t19 + t247, pkin(3) * t11 - pkin(4) * t70 + pkin(8) * t19, t95, t64, t88, t93, t89, t107, -qJ(6) * t280 + t208 * t30 + t249, t205 * t37 + t208 * t72 + t248, t20 * t208 + t205 * t21 + t247, pkin(3) * t6 - pkin(4) * t39 + pkin(8) * t9 - qJ(6) * t288 + t13 * t208; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, t137, t133, -t130, 0, 0, 0, 0, 0, 0, t79, t83, t65, t18, 0, 0, 0, 0, 0, 0, t79, t83, t65, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, t131, t103, -t134, -t100, t151, -t40, -t41, 0, 0, t134, t131, t103, -t134, -t100, t151, t158 + t218, t305, -t292, t297; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t99, t303, t115, t39;];
tauJ_reg  = t12;
