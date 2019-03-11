% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRPRP6
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
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPRP6_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRP6_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP6_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:19:49
% EndTime: 2019-03-09 03:19:55
% DurationCPUTime: 3.40s
% Computational Cost: add. (3966->388), mult. (9118->466), div. (0->0), fcn. (6745->10), ass. (0->207)
t146 = sin(pkin(9));
t147 = cos(pkin(9));
t151 = sin(qJ(3));
t261 = cos(qJ(3));
t105 = t261 * t146 + t151 * t147;
t277 = t105 * qJD(1);
t284 = qJD(5) + t277;
t153 = cos(qJ(5));
t150 = sin(qJ(5));
t216 = t150 * qJD(3);
t205 = t261 * t147;
t189 = qJD(1) * t205;
t233 = t151 * t146;
t204 = qJD(1) * t233;
t97 = -t189 + t204;
t70 = -t153 * t97 + t216;
t198 = t284 * t70;
t102 = t105 * qJD(3);
t199 = qJDD(1) * t261;
t213 = t146 * qJDD(1);
t183 = -t147 * t199 + t151 * t213;
t158 = qJD(1) * t102 + t183;
t218 = qJD(5) * t153;
t27 = qJD(5) * t216 - t153 * qJDD(3) - t150 * t158 - t97 * t218;
t287 = t27 - t198;
t142 = pkin(9) + qJ(3);
t137 = cos(t142);
t131 = g(3) * t137;
t136 = sin(t142);
t154 = cos(qJ(1));
t237 = t136 * t154;
t152 = sin(qJ(1));
t238 = t136 * t152;
t206 = -g(1) * t237 - g(2) * t238 + t131;
t242 = qJDD(1) * pkin(1);
t135 = qJDD(2) - t242;
t281 = t277 * qJD(3);
t159 = -t183 - t281;
t212 = t147 * qJDD(1);
t207 = qJD(3) * t189 + t146 * t199 + t151 * t212;
t64 = qJD(3) * t204 - t207;
t21 = -pkin(2) * t212 - pkin(3) * t159 + t64 * qJ(4) - qJD(4) * t277 + t135;
t11 = pkin(8) * t158 + t21;
t251 = pkin(7) + qJ(2);
t112 = t251 * t146;
t106 = qJD(1) * t112;
t113 = t251 * t147;
t107 = qJD(1) * t113;
t203 = qJD(3) * t261;
t220 = qJD(3) * t151;
t214 = qJD(1) * qJD(2);
t270 = t251 * qJDD(1) + t214;
t78 = t270 * t146;
t79 = t270 * t147;
t175 = -t106 * t220 + t107 * t203 + t151 * t79 + t261 * t78;
t170 = qJDD(4) + t175;
t266 = pkin(3) + pkin(8);
t18 = -t64 * pkin(4) - t266 * qJDD(3) + t170;
t65 = t261 * t106 + t151 * t107;
t280 = qJD(4) + t65;
t182 = pkin(4) * t277 + t280;
t37 = -t266 * qJD(3) + t182;
t211 = -t153 * t11 - t150 * t18 - t37 * t218;
t219 = qJD(5) * t150;
t133 = t147 * pkin(2) + pkin(1);
t109 = -t133 * qJD(1) + qJD(2);
t167 = -qJ(4) * t277 + t109;
t32 = t266 * t97 + t167;
t174 = -t32 * t219 - t211;
t191 = t150 * qJDD(3) - t153 * t158;
t72 = t153 * qJD(3) + t150 * t97;
t28 = qJD(5) * t72 + t191;
t2 = -t28 * qJ(6) - t70 * qJD(6) + t174;
t14 = -t150 * t32 + t153 * t37;
t7 = -t72 * qJ(6) + t14;
t6 = pkin(5) * t284 + t7;
t273 = t284 * t6 - t2;
t15 = t150 * t37 + t153 * t32;
t200 = -t150 * t11 + t153 * t18;
t162 = -qJD(5) * t15 + t200;
t62 = -qJDD(5) + t64;
t1 = -t62 * pkin(5) + t27 * qJ(6) - t72 * qJD(6) + t162;
t8 = -t70 * qJ(6) + t15;
t274 = t284 * t8 + t1;
t289 = -t273 * t150 + t274 * t153 + t206;
t282 = t153 * t284;
t288 = t72 * t282;
t195 = t150 * t284;
t56 = t153 * t62;
t176 = -t195 * t284 - t56;
t279 = g(1) * t152 - g(2) * t154;
t177 = -t135 + t279;
t286 = -t183 - 0.2e1 * t281;
t258 = g(1) * t154;
t185 = g(2) * t152 + t258;
t256 = g(3) * t136;
t166 = -t137 * t185 - t256;
t197 = t106 * t203 + t107 * t220 + t151 * t78 - t261 * t79;
t285 = -t65 * qJD(3) - t166 + t197;
t223 = t137 * pkin(3) + t136 * qJ(4);
t276 = qJ(2) * qJDD(1);
t229 = t154 * t153;
t232 = t152 * t150;
t89 = t136 * t229 - t232;
t230 = t154 * t150;
t231 = t152 * t153;
t91 = t136 * t231 + t230;
t275 = -g(1) * t89 - g(2) * t91 + t153 * t131;
t145 = qJD(3) * qJ(4);
t66 = -t151 * t106 + t261 * t107;
t47 = -t97 * pkin(4) + t66;
t40 = t145 + t47;
t272 = t266 * t62 + t284 * t40;
t50 = (qJD(2) * t146 + qJD(3) * t113) * t151 - qJD(2) * t205 + t112 * t203;
t68 = -t151 * t112 + t261 * t113;
t271 = t50 * qJD(3) - t68 * qJDD(3) - t136 * t279;
t269 = t72 ^ 2;
t268 = t97 ^ 2;
t267 = t277 ^ 2;
t265 = -t7 + t6;
t215 = t153 * qJD(6);
t226 = qJ(6) + t266;
t42 = t153 * t47;
t244 = t97 * qJ(4);
t43 = t266 * t277 + t244;
t262 = t226 * t219 - t215 + t97 * pkin(5) - t42 - (-qJ(6) * t277 - t43) * t150;
t260 = pkin(5) * t150;
t255 = t72 * t97;
t254 = t97 * t70;
t253 = t277 * t97;
t148 = -qJ(6) - pkin(8);
t252 = pkin(3) - t148;
t250 = t150 * t47 + t153 * t43;
t104 = -t205 + t233;
t178 = -t105 * qJ(4) - t133;
t45 = t266 * t104 + t178;
t67 = t261 * t112 + t151 * t113;
t54 = t105 * pkin(4) + t67;
t249 = t150 * t54 + t153 * t45;
t111 = t226 * t153;
t243 = qJ(6) * t153;
t248 = -qJD(5) * t111 - t150 * qJD(6) - t243 * t277 - t250;
t247 = t150 * t62;
t246 = t153 * t27;
t241 = qJDD(3) * pkin(3);
t240 = t102 * t150;
t239 = t102 * t153;
t236 = t137 * t148;
t235 = t137 * t152;
t234 = t137 * t154;
t227 = t66 * qJD(3);
t134 = t153 * pkin(5) + pkin(4);
t222 = t134 + t251;
t221 = t146 ^ 2 + t147 ^ 2;
t217 = qJD(5) * t266;
t210 = t137 * t260;
t208 = pkin(3) * t234 + qJ(4) * t237 + t154 * t133;
t202 = qJ(4) + t260;
t196 = -qJ(6) * t104 - t45;
t192 = t221 * qJD(1) ^ 2;
t190 = 0.2e1 * t221;
t188 = -g(1) * t235 + g(2) * t234;
t101 = t146 * t220 - t147 * t203;
t181 = t101 * qJ(4) - t105 * qJD(4);
t143 = qJDD(3) * qJ(4);
t144 = qJD(3) * qJD(4);
t22 = -t143 - t144 + t197;
t29 = t266 * t102 + t181;
t51 = qJD(2) * t105 + qJD(3) * t68;
t34 = -t101 * pkin(4) + t51;
t173 = t150 * t34 + t153 * t29 + t54 * t218 - t45 * t219;
t172 = t104 * t219 - t239;
t171 = -t282 * t284 + t247;
t168 = t177 + t242;
t165 = -t175 - t206;
t164 = -t51 * qJD(3) - t67 * qJDD(3) - t188;
t19 = -pkin(4) * t158 - t22;
t163 = t19 + t166;
t161 = t190 * t214 - t185;
t33 = -t102 * pkin(4) - t50;
t49 = t97 * pkin(3) + t167;
t160 = t277 * t49 + qJDD(4) - t165;
t5 = t28 * pkin(5) + qJDD(6) + t19;
t116 = qJ(4) * t234;
t114 = qJ(4) * t235;
t110 = t226 * t150;
t108 = -t133 * qJDD(1) + qJDD(2);
t92 = -t136 * t232 + t229;
t90 = t136 * t230 + t231;
t83 = qJD(3) * t97;
t69 = t70 ^ 2;
t63 = t104 * pkin(3) + t178;
t61 = pkin(3) * t277 + t244;
t60 = -t145 - t66;
t59 = -qJD(3) * pkin(3) + t280;
t55 = -t104 * pkin(4) + t68;
t53 = t153 * t54;
t44 = t102 * pkin(3) + t181;
t31 = t153 * t34;
t25 = t153 * t28;
t24 = t70 * pkin(5) + qJD(6) + t40;
t23 = t170 - t241;
t20 = t104 * t243 + t249;
t12 = t105 * pkin(5) + t196 * t150 + t53;
t4 = -qJ(6) * t172 + t104 * t215 + t173;
t3 = -t101 * pkin(5) + t31 + t196 * t218 + (-qJ(6) * t102 - qJD(5) * t54 - qJD(6) * t104 - t29) * t150;
t9 = [qJDD(1), t279, t185, t168 * t147, -t168 * t146, t190 * t276 + t161, pkin(1) * t177 + (t221 * t276 + t161) * qJ(2), -t101 * t277 - t64 * t105, t101 * t97 - t102 * t277 + t64 * t104 + t105 * t159, -t101 * qJD(3) + t105 * qJDD(3), -t102 * qJD(3) - t104 * qJDD(3), 0, t109 * t102 + t108 * t104 - t133 * t158 + t164, -t109 * t101 + t108 * t105 + t133 * t64 + t271, -t59 * t101 + t60 * t102 + t22 * t104 + t23 * t105 - t158 * t68 + t277 * t51 + t50 * t97 - t67 * t64 - t185, -t49 * t102 - t21 * t104 + t159 * t63 - t44 * t97 - t164, t49 * t101 - t21 * t105 - t277 * t44 + t63 * t64 - t271, t21 * t63 + t49 * t44 - t22 * t68 + t60 * t50 + t23 * t67 + t59 * t51 - t251 * t258 - g(2) * t208 + (-g(1) * (-t133 - t223) - g(2) * t251) * t152, t72 * t240 + (-t150 * t27 + t218 * t72) * t104 (-t150 * t70 + t153 * t72) * t102 + (-t150 * t28 - t246 + (-t150 * t72 - t153 * t70) * qJD(5)) * t104, t284 * t240 - t72 * t101 - t27 * t105 + (t218 * t284 - t247) * t104, t284 * t239 + t70 * t101 - t28 * t105 + (-t219 * t284 - t56) * t104, -t101 * t284 - t105 * t62 (-t150 * t29 + t31) * t284 - (-t150 * t45 + t53) * t62 + t200 * t105 - t14 * t101 + t33 * t70 + t55 * t28 - g(1) * t92 - g(2) * t90 + (-t102 * t40 - t104 * t19) * t153 + (t40 * t150 * t104 - t105 * t15 - t249 * t284) * qJD(5), -t173 * t284 + t249 * t62 - t174 * t105 + t15 * t101 + t33 * t72 - t55 * t27 + t40 * t240 + g(1) * t91 - g(2) * t89 + (t19 * t150 + t218 * t40) * t104, t12 * t27 - t20 * t28 - t3 * t72 - t4 * t70 + (-t150 * t6 + t153 * t8) * t102 + (-t1 * t150 + t2 * t153 + (-t150 * t8 - t153 * t6) * qJD(5)) * t104 - t188, t2 * t20 + t8 * t4 + t1 * t12 + t6 * t3 + t5 * (-t104 * t134 + t68) + t24 * (pkin(5) * t172 + t33) - g(1) * (t222 * t154 + (-t136 * t202 - t137 * t252 - t133) * t152) - g(2) * ((t136 * t260 - t236) * t154 + t222 * t152 + t208); 0, 0, 0, -t212, t213, -t192, -qJ(2) * t192 - t177, 0, 0, 0, 0, 0, -t286 (-t97 - t204) * qJD(3) + t207, -t267 - t268, t286, t64 + t83, -t277 * t59 - t60 * t97 + t21 - t279, 0, 0, 0, 0, 0, t171 + t254, t255 - t176, -t150 * t287 - t25 + t288, -t150 * t274 - t153 * t273 + t24 * t97 - t279; 0, 0, 0, 0, 0, 0, 0, t253, t267 - t268 (t97 - t204) * qJD(3) + t207, -t183, qJDD(3), -t109 * t277 + t165 + t227, t109 * t97 + t285, pkin(3) * t64 - qJ(4) * t158 - (t60 + t66) * t277 + (t59 - t280) * t97, t61 * t97 + t160 - t227 - 0.2e1 * t241, t277 * t61 - t49 * t97 + 0.2e1 * t143 + 0.2e1 * t144 - t285, -t22 * qJ(4) - t23 * pkin(3) - t49 * t61 - t59 * t66 - g(1) * (-pkin(3) * t237 + t116) - g(2) * (-pkin(3) * t238 + t114) - g(3) * t223 - t280 * t60, -t195 * t72 - t246, -t25 - t288 + (t27 + t198) * t150, t176 + t255, t171 - t254, t284 * t97, qJ(4) * t28 + t14 * t97 - t42 * t284 + t182 * t70 + t272 * t153 + ((t43 + t217) * t284 + t163) * t150, -qJ(4) * t27 + t250 * t284 - t15 * t97 + t182 * t72 - t272 * t150 + (t217 * t284 + t163) * t153, t110 * t28 - t111 * t27 - t248 * t70 - t262 * t72 - t289, -t2 * t110 - t1 * t111 + t5 * t202 - g(1) * (t154 * t210 + t116) - g(2) * (t152 * t210 + t114) - g(3) * (t223 - t236) + t248 * t8 + t262 * t6 + (pkin(5) * t282 + t182) * t24 + (-g(3) * t260 + t185 * t252) * t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64 + t83, qJDD(3) - t253, -qJD(3) ^ 2 - t267, t60 * qJD(3) + t160 - t241, 0, 0, 0, 0, 0, -qJD(3) * t70 + t176, -qJD(3) * t72 + t171, t287 * t153 + (t284 * t72 - t28) * t150, -t24 * qJD(3) + t289; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 * t70, -t69 + t269, -t287, -t191 + (-qJD(5) + t284) * t72, -t62, t15 * t284 - t40 * t72 + t162 + t275, g(1) * t90 - g(2) * t92 + t14 * t284 + t40 * t70 + (qJD(5) * t32 - t131) * t150 + t211, pkin(5) * t27 - t265 * t70, t265 * t8 + (-t24 * t72 + t1 + t275) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t69 - t269, -g(1) * t234 - g(2) * t235 + t6 * t72 + t8 * t70 - t256 + t5;];
tau_reg  = t9;
