% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPRRP3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 23:46
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPRRP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:44:49
% EndTime: 2019-05-04 23:44:59
% DurationCPUTime: 3.66s
% Computational Cost: add. (14278->362), mult. (32119->518), div. (0->0), fcn. (24714->12), ass. (0->231)
t204 = sin(pkin(11));
t206 = cos(pkin(11));
t210 = sin(qJ(4));
t213 = cos(qJ(4));
t229 = t204 * t213 + t206 * t210;
t181 = t229 * qJD(2);
t209 = sin(qJ(5));
t212 = cos(qJ(5));
t165 = -t212 * qJD(4) + t181 * t209;
t167 = qJD(4) * t209 + t181 * t212;
t136 = t167 * t165;
t240 = t206 * qJDD(2);
t241 = t204 * qJDD(2);
t230 = t210 * t241 - t213 * t240;
t244 = t181 * qJD(4);
t155 = -t230 - t244;
t148 = qJDD(5) - t155;
t292 = -t136 + t148;
t302 = pkin(5) * t292;
t178 = t229 * qJDD(2);
t248 = qJD(2) * t206;
t179 = qJD(2) * t204 * t210 - t213 * t248;
t245 = t179 * qJD(4);
t157 = t178 - t245;
t126 = -qJD(5) * t165 + qJDD(4) * t209 + t157 * t212;
t175 = qJD(5) + t179;
t143 = t175 * t165;
t104 = t126 + t143;
t301 = qJ(6) * t104;
t158 = t181 * t179;
t290 = qJDD(4) - t158;
t300 = t210 * t290;
t299 = t213 * t290;
t199 = t204 ^ 2;
t200 = t206 ^ 2;
t289 = qJD(2) ^ 2;
t186 = (t199 + t200) * t289;
t260 = t292 * t209;
t259 = t292 * t212;
t205 = sin(pkin(6));
t207 = cos(pkin(6));
t263 = sin(pkin(10));
t264 = cos(pkin(10));
t225 = t263 * g(1) - t264 * g(2);
t224 = t207 * t225;
t250 = -g(3) + qJDD(1);
t298 = t205 * t250 + t224;
t164 = t167 ^ 2;
t174 = t175 ^ 2;
t131 = -t164 - t174;
t163 = t165 ^ 2;
t232 = -t212 * qJDD(4) + t157 * t209;
t125 = -qJD(5) * t167 - t232;
t139 = pkin(5) * t175 - qJ(6) * t167;
t149 = pkin(4) * t179 - pkin(9) * t181;
t288 = qJD(4) ^ 2;
t187 = -t264 * g(1) - t263 * g(2);
t211 = sin(qJ(2));
t214 = cos(qJ(2));
t146 = t214 * t187 + t298 * t211;
t219 = qJDD(2) * qJ(3) + t146;
t222 = -t205 * t225 + t207 * t250;
t287 = 2 * qJD(3);
t115 = t206 * (-t289 * pkin(2) + t219) + t204 * t222 + t248 * t287;
t243 = t200 * t289;
t107 = -pkin(3) * t243 + pkin(8) * t240 + t115;
t220 = t206 * t222;
t295 = qJ(3) + pkin(8);
t217 = t220 + (-t295 * qJDD(2) + (-(2 * qJD(3)) + (pkin(3) * t206 + pkin(2)) * qJD(2)) * qJD(2) - t146) * t204;
t74 = t213 * t107 + t210 * t217;
t64 = -t288 * pkin(4) + qJDD(4) * pkin(9) - t179 * t149 + t74;
t203 = qJDD(2) * pkin(2);
t231 = t187 * t211 - t298 * t214;
t138 = -t289 * qJ(3) + qJDD(3) - t203 + t231;
t128 = -pkin(3) * t240 + t138 + (-t199 * t289 - t243) * pkin(8);
t77 = (-t157 + t245) * pkin(9) + (-t155 + t244) * pkin(4) + t128;
t43 = t209 * t77 + t212 * t64;
t227 = t125 * qJ(6) - 0.2e1 * qJD(6) * t165 - t139 * t175 + t43;
t297 = -t227 + (t131 + t163) * pkin(5);
t293 = t126 - t143;
t101 = (qJD(5) - t175) * t167 + t232;
t176 = t179 ^ 2;
t177 = t181 ^ 2;
t127 = -t174 - t163;
t82 = t127 * t209 + t259;
t286 = pkin(4) * t82;
t111 = t136 + t148;
t262 = t111 * t209;
t87 = t131 * t212 - t262;
t285 = pkin(4) * t87;
t246 = qJD(6) * t167;
t160 = -0.2e1 * t246;
t42 = t209 * t64 - t212 * t77;
t226 = -t301 - t42 + t302;
t26 = t160 + t226;
t284 = pkin(5) * t26;
t121 = -t163 - t164;
t69 = -t101 * t212 + t104 * t209;
t53 = -t121 * t213 + t210 * t69;
t283 = pkin(8) * t53;
t100 = (qJD(5) + t175) * t167 + t232;
t83 = t127 * t212 - t260;
t56 = -t100 * t213 + t210 * t83;
t282 = pkin(8) * t56;
t261 = t111 * t212;
t88 = -t131 * t209 - t261;
t60 = t210 * t88 - t213 * t293;
t281 = pkin(8) * t60;
t67 = -t101 * t209 - t104 * t212;
t280 = pkin(9) * t67;
t279 = pkin(9) * t82;
t278 = pkin(9) * t87;
t277 = pkin(4) * t210;
t276 = pkin(5) * t104;
t54 = t121 * t210 + t213 * t69;
t28 = -t204 * t53 + t206 * t54;
t275 = -pkin(2) * t67 + qJ(3) * t28;
t57 = t100 * t210 + t213 * t83;
t35 = -t204 * t56 + t206 * t57;
t274 = -pkin(2) * t82 + qJ(3) * t35;
t61 = t210 * t293 + t213 * t88;
t37 = -t204 * t60 + t206 * t61;
t273 = -pkin(2) * t87 + qJ(3) * t37;
t272 = -pkin(4) * t100 + pkin(9) * t83;
t271 = -pkin(4) * t293 + pkin(9) * t88;
t73 = t107 * t210 - t213 * t217;
t47 = t210 * t74 - t213 * t73;
t270 = t204 * t47;
t269 = t209 * t26;
t63 = -qJDD(4) * pkin(4) - t288 * pkin(9) + t149 * t181 + t73;
t268 = t209 * t63;
t267 = t212 * t26;
t266 = t212 * t63;
t265 = -pkin(4) * t121 + pkin(9) * t69;
t258 = t128 * t210;
t257 = t128 * t213;
t152 = qJDD(4) + t158;
t256 = t152 * t210;
t255 = t152 * t213;
t254 = t175 * t209;
t253 = t175 * t212;
t239 = t210 * t136;
t238 = t213 * t136;
t237 = -pkin(4) * t213 - pkin(3);
t236 = -pkin(3) * t82 + pkin(8) * t57;
t235 = -pkin(3) * t87 + pkin(8) * t61;
t19 = t209 * t42 + t212 * t43;
t48 = t210 * t73 + t213 * t74;
t233 = -t138 + t203;
t114 = -t220 + ((-pkin(2) * qJD(2) + t287) * qJD(2) + t219) * t204;
t78 = t114 * t204 + t206 * t115;
t13 = t19 * t210 - t213 * t63;
t14 = t19 * t213 + t210 * t63;
t3 = -t13 * t204 + t14 * t206;
t18 = t209 * t43 - t212 * t42;
t223 = t226 + t302;
t44 = -t125 * pkin(5) - t163 * qJ(6) + t139 * t167 + qJDD(6) + t63;
t196 = t200 * qJDD(2);
t195 = t199 * qJDD(2);
t185 = t196 + t195;
t184 = t206 * t186;
t183 = t204 * t186;
t171 = -t177 - t288;
t170 = -t177 + t288;
t169 = t176 - t288;
t161 = 0.2e1 * t246;
t156 = t178 - 0.2e1 * t245;
t154 = t230 + 0.2e1 * t244;
t150 = -t288 - t176;
t141 = -t164 + t174;
t140 = t163 - t174;
t135 = -t176 - t177;
t133 = t164 - t163;
t130 = -t171 * t210 - t255;
t129 = t171 * t213 - t256;
t120 = t178 * t210 - t213 * t230;
t119 = -t178 * t213 - t210 * t230;
t117 = t150 * t213 - t300;
t116 = t150 * t210 + t299;
t109 = (-t165 * t212 + t167 * t209) * t175;
t108 = (-t165 * t209 - t167 * t212) * t175;
t97 = t126 * t212 - t167 * t254;
t96 = t126 * t209 + t167 * t253;
t95 = -t125 * t209 + t165 * t253;
t94 = t125 * t212 + t165 * t254;
t93 = -t129 * t204 + t130 * t206;
t92 = t140 * t212 - t262;
t91 = -t141 * t209 + t259;
t90 = t140 * t209 + t261;
t89 = t141 * t212 + t260;
t84 = -t119 * t204 + t120 * t206;
t81 = -t116 * t204 + t117 * t206;
t71 = -pkin(5) * t293 - qJ(6) * t111;
t70 = -t100 * t212 - t209 * t293;
t68 = -t100 * t209 + t212 * t293;
t58 = t204 * (t109 * t213 + t148 * t210) + t206 * (t109 * t210 - t148 * t213);
t52 = pkin(8) * t54;
t51 = -pkin(4) * t67 + t276;
t50 = t204 * (t213 * t97 + t239) + t206 * (t210 * t97 - t238);
t49 = t204 * (t213 * t95 - t239) + t206 * (t210 * t95 + t238);
t46 = t266 - t278;
t45 = t268 - t279;
t40 = -qJ(6) * t131 + t44;
t39 = t204 * (-t101 * t210 + t213 * t92) + t206 * (t101 * t213 + t210 * t92);
t38 = t204 * (t104 * t210 + t213 * t91) + t206 * (-t104 * t213 + t210 * t91);
t33 = t204 * (t133 * t210 + t213 * t70) + t206 * (-t133 * t213 + t210 * t70);
t32 = t43 - t285;
t31 = t42 - t286;
t30 = -pkin(5) * t100 + qJ(6) * t127 - t44;
t29 = -pkin(5) * t163 + t227;
t25 = t161 - t226 + t301;
t24 = t206 * t48 - t270;
t23 = -qJ(6) * t101 + (-t121 - t163) * pkin(5) + t227;
t22 = -t285 - t297;
t21 = -t209 * t71 + t212 * t40 - t278;
t20 = -qJ(6) * t259 - t209 * t30 - t279;
t17 = t161 - t223 - t286;
t16 = -pkin(5) * t44 + qJ(6) * t29;
t15 = -t18 - t280;
t12 = t207 * (t204 * t61 + t206 * t60) + (t211 * t37 - t214 * t87) * t205;
t11 = t212 * t29 - t269;
t10 = t209 * t29 + t267;
t9 = t207 * (t204 * t57 + t206 * t56) + (t211 * t35 - t214 * t82) * t205;
t8 = t207 * (t204 * t54 + t206 * t53) + (t211 * t28 - t214 * t67) * t205;
t7 = t11 * t213 + t210 * t44;
t6 = t11 * t210 - t213 * t44;
t5 = -t209 * t23 + t212 * t25 - t280;
t4 = -pkin(4) * t10 - t284;
t2 = -pkin(9) * t10 - qJ(6) * t267 - t16 * t209;
t1 = -t204 * t6 + t206 * t7;
t27 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t250, 0, 0, 0, 0, 0, 0, (qJDD(2) * t214 - t289 * t211) * t205, (-qJDD(2) * t211 - t289 * t214) * t205, 0, t207 ^ 2 * t250 + (t211 * t146 - t214 * t231 - t224) * t205, 0, 0, 0, 0, 0, 0, (-t184 * t211 + t214 * t240) * t205, (t183 * t211 - t214 * t241) * t205, (t185 * t211 + t186 * t214) * t205, t207 * (-t114 * t206 + t115 * t204) + (-t138 * t214 + t211 * t78) * t205, 0, 0, 0, 0, 0, 0, t207 * (t116 * t206 + t117 * t204) + (-t154 * t214 + t211 * t81) * t205, t207 * (t129 * t206 + t130 * t204) + (-t156 * t214 + t211 * t93) * t205, t207 * (t119 * t206 + t120 * t204) + (-t135 * t214 + t211 * t84) * t205, t207 * (t204 * t48 + t206 * t47) + (-t128 * t214 + t211 * t24) * t205, 0, 0, 0, 0, 0, 0, t9, t12, t8, t207 * (t13 * t206 + t14 * t204) + (-t18 * t214 + t211 * t3) * t205, 0, 0, 0, 0, 0, 0, t9, t12, t8, t207 * (t204 * t7 + t206 * t6) + (t1 * t211 - t10 * t214) * t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t231, -t146, 0, 0, t195, 0.2e1 * t204 * t240, 0, t196, 0, 0, -qJ(3) * t184 + t233 * t206, qJ(3) * t183 - t233 * t204, pkin(2) * t186 + qJ(3) * t185 + t78, -pkin(2) * t138 + qJ(3) * t78, t204 * (t157 * t213 - t210 * t244) + t206 * (t157 * t210 + t213 * t244), t204 * (-t154 * t213 - t156 * t210) + t206 * (-t154 * t210 + t156 * t213), t204 * (-t170 * t210 + t299) + t206 * (t170 * t213 + t300), t204 * (-t155 * t210 + t213 * t245) + t206 * (t155 * t213 + t210 * t245), t204 * (t169 * t213 - t256) + t206 * (t169 * t210 + t255), (t204 * (-t179 * t213 + t181 * t210) + t206 * (-t179 * t210 - t181 * t213)) * qJD(4), t204 * (-pkin(8) * t116 + t258) + t206 * (-pkin(3) * t154 + pkin(8) * t117 - t257) - pkin(2) * t154 + qJ(3) * t81, t204 * (-pkin(8) * t129 + t257) + t206 * (-pkin(3) * t156 + pkin(8) * t130 + t258) - pkin(2) * t156 + qJ(3) * t93, t204 * (-pkin(8) * t119 - t47) + t206 * (-pkin(3) * t135 + pkin(8) * t120 + t48) - pkin(2) * t135 + qJ(3) * t84, -pkin(8) * t270 + t206 * (-pkin(3) * t128 + pkin(8) * t48) - pkin(2) * t128 + qJ(3) * t24, t50, t33, t38, t49, t39, t58, t204 * (-t210 * t31 + t213 * t45 - t282) + t206 * (t210 * t45 + t213 * t31 + t236) + t274, t204 * (-t210 * t32 + t213 * t46 - t281) + t206 * (t210 * t46 + t213 * t32 + t235) + t273, t204 * (t15 * t213 + t67 * t277 - t283) + t206 * (t15 * t210 + t237 * t67 + t52) + t275, (t204 * (-pkin(9) * t213 + t277) + t206 * (-pkin(9) * t210 + t237) - pkin(2)) * t18 + t295 * t3, t50, t33, t38, t49, t39, t58, t204 * (-t17 * t210 + t20 * t213 - t282) + t206 * (t17 * t213 + t20 * t210 + t236) + t274, t204 * (t21 * t213 - t210 * t22 - t281) + t206 * (t21 * t210 + t213 * t22 + t235) + t273, t204 * (-t210 * t51 + t213 * t5 - t283) + t206 * (-pkin(3) * t67 + t210 * t5 + t213 * t51 + t52) + t275, t204 * (-pkin(8) * t6 + t2 * t213 - t210 * t4) + t206 * (-pkin(3) * t10 + pkin(8) * t7 + t2 * t210 + t213 * t4) - pkin(2) * t10 + qJ(3) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t240, t241, -t186, t138, 0, 0, 0, 0, 0, 0, t154, t156, t135, t128, 0, 0, 0, 0, 0, 0, t82, t87, t67, t18, 0, 0, 0, 0, 0, 0, t82, t87, t67, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t158, t177 - t176, t178, -t158, -t230, qJDD(4), -t73, -t74, 0, 0, t96, t68, t89, t94, t90, t108, -t266 + t272, t268 + t271, t19 + t265, -pkin(4) * t63 + pkin(9) * t19, t96, t68, t89, t94, t90, t108, -qJ(6) * t260 + t212 * t30 + t272, t209 * t40 + t212 * t71 + t271, t209 * t25 + t212 * t23 + t265, -pkin(4) * t44 + pkin(9) * t11 - qJ(6) * t269 + t16 * t212; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, t133, t104, -t136, -t101, t148, -t42, -t43, 0, 0, t136, t133, t104, -t136, -t101, t148, t160 + t223, t297, -t276, t284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t293, t121, t44;];
tauJ_reg  = t27;
