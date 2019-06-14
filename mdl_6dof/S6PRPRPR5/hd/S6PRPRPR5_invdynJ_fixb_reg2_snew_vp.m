% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 23:01
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRPRPR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:59:09
% EndTime: 2019-05-04 22:59:22
% DurationCPUTime: 6.24s
% Computational Cost: add. (10797->371), mult. (24728->525), div. (0->0), fcn. (18547->12), ass. (0->212)
t199 = sin(pkin(11));
t201 = cos(pkin(11));
t210 = qJD(4) ^ 2;
t205 = sin(qJ(4));
t208 = cos(qJ(4));
t250 = t199 * t208;
t225 = t201 * t205 + t250;
t173 = t225 * qJD(2);
t268 = t173 ^ 2;
t153 = t268 + t210;
t251 = t199 * t205;
t171 = (-t201 * t208 + t251) * qJD(2);
t253 = t173 * t171;
t285 = qJDD(4) + t253;
t299 = t205 * t285;
t95 = t208 * t153 + t299;
t297 = t208 * t285;
t97 = -t205 * t153 + t297;
t61 = t199 * t95 - t201 * t97;
t314 = qJ(3) * t61;
t206 = sin(qJ(2));
t313 = t206 * t61;
t202 = cos(pkin(6));
t312 = t202 * (t199 * t97 + t201 * t95);
t311 = pkin(8) * t95;
t310 = pkin(8) * t97;
t154 = t268 - t210;
t286 = qJDD(4) - t253;
t296 = t208 * t286;
t298 = t205 * t286;
t308 = t199 * (t205 * t154 + t296) - t201 * (t208 * t154 - t298);
t269 = t171 ^ 2;
t149 = t269 - t210;
t305 = t199 * (-t208 * t149 + t299) - t201 * (t205 * t149 + t297);
t125 = -t210 - t269;
t81 = t205 * t125 + t296;
t84 = -t208 * t125 + t298;
t51 = t199 * t81 + t201 * t84;
t304 = qJ(3) * t51;
t303 = t206 * t51;
t300 = t202 * (t199 * t84 - t201 * t81);
t160 = qJD(4) * t171;
t170 = t225 * qJDD(2);
t135 = t170 - t160;
t277 = t160 - t135;
t295 = t277 * qJ(5);
t294 = pkin(8) * t81;
t293 = pkin(8) * t84;
t256 = sin(pkin(10));
t257 = cos(pkin(10));
t180 = -t257 * g(1) - t256 * g(2);
t209 = cos(qJ(2));
t200 = sin(pkin(6));
t224 = t256 * g(1) - t257 * g(2);
t243 = -g(3) + qJDD(1);
t276 = t200 * t243 + t202 * t224;
t114 = t209 * t180 + t276 * t206;
t211 = qJD(2) ^ 2;
t278 = -t211 * pkin(2) + qJDD(2) * qJ(3) + 0.2e1 * qJD(2) * qJD(3) + t114;
t270 = -t269 - t268;
t283 = pkin(2) * t270;
t282 = pkin(3) * t270;
t204 = sin(qJ(6));
t207 = cos(qJ(6));
t142 = t204 * qJD(4) - t207 * t171;
t144 = t207 * qJD(4) + t204 * t171;
t106 = t144 * t142;
t122 = qJDD(6) + t135;
t273 = -t106 + t122;
t281 = t204 * t273;
t280 = t207 * t273;
t279 = t209 * t270;
t164 = qJD(6) + t173;
t112 = t164 * t142;
t238 = t201 * qJDD(2);
t239 = t199 * qJDD(2);
t226 = t205 * t239 - t208 * t238;
t242 = t173 * qJD(4);
t133 = t226 + t242;
t91 = -t142 * qJD(6) + t207 * qJDD(4) + t204 * t133;
t275 = -t112 + t91;
t274 = t201 * t211;
t101 = t160 + t135;
t193 = t199 ^ 2;
t194 = t201 ^ 2;
t272 = t193 + t194;
t271 = t268 - t269;
t179 = t272 * t211;
t140 = t142 ^ 2;
t141 = t144 ^ 2;
t162 = t164 ^ 2;
t267 = 2 * qJD(5);
t266 = -pkin(4) - pkin(9);
t265 = pkin(4) * t208;
t148 = -t200 * t224 + t202 * t243;
t146 = t201 * t148;
t72 = t146 + (pkin(3) * t274 - pkin(8) * qJDD(2) - t278) * t199;
t252 = t194 * t211;
t80 = t199 * t148 + t278 * t201;
t73 = -pkin(3) * t252 + pkin(8) * t238 + t80;
t46 = t205 * t73 - t208 * t72;
t47 = t205 * t72 + t208 * t73;
t24 = t205 * t47 - t208 * t46;
t264 = t199 * t24;
t147 = t173 * pkin(5) - qJD(4) * pkin(9);
t121 = t171 * pkin(4) - t173 * qJ(5);
t222 = -t210 * pkin(4) - t171 * t121 + t47;
t240 = qJDD(4) * qJ(5);
t28 = t240 - t133 * pkin(5) - t269 * pkin(9) + (t267 + t147) * qJD(4) + t222;
t263 = t204 * t28;
t77 = t106 + t122;
t262 = t204 * t77;
t198 = qJDD(2) * pkin(2);
t227 = t206 * t180 - t276 * t209;
t108 = -t211 * qJ(3) + qJDD(3) - t198 + t227;
t93 = -pkin(3) * t238 + t108 + (-t193 * t211 - t252) * pkin(8);
t261 = t205 * t93;
t260 = t207 * t28;
t259 = t207 * t77;
t258 = t208 * t93;
t255 = t164 * t204;
t254 = t164 * t207;
t237 = -t141 - t162;
t235 = t205 * t106;
t234 = t208 * t106;
t233 = qJ(5) * t205 + pkin(3);
t79 = t278 * t199 - t146;
t49 = t199 * t79 + t201 * t80;
t39 = -qJDD(4) * pkin(4) - t210 * qJ(5) + t173 * t121 + qJDD(5) + t46;
t27 = pkin(5) * t101 - pkin(9) * t286 + t39;
t218 = t133 * pkin(4) + t295 + t93;
t230 = pkin(4) * qJD(4) - (2 * qJD(5));
t33 = (-t147 + t230) * t173 + t218 - t269 * pkin(5) + t133 * pkin(9);
t16 = t204 * t33 - t207 * t27;
t25 = t205 * t46 + t208 * t47;
t231 = -t108 + t198;
t229 = t204 * qJDD(4) - t207 * t133;
t18 = t204 * t27 + t207 * t33;
t7 = -t207 * t16 + t204 * t18;
t8 = t204 * t16 + t207 * t18;
t220 = qJD(4) * t267 + t222;
t38 = t220 + t240;
t21 = t205 * t38 - t208 * t39;
t22 = t205 * t39 + t208 * t38;
t9 = -t199 * t21 + t201 * t22;
t223 = (-qJD(6) + t164) * t144 - t229;
t219 = t199 * (t208 * t135 - t205 * t242) + t201 * (t205 * t135 + t208 * t242);
t216 = t199 * (t205 * t133 + t208 * t160) + t201 * (-t208 * t133 + t205 * t160);
t215 = (t199 * (-t171 * t208 + t173 * t205) + t201 * (-t171 * t205 - t173 * t208)) * qJD(4);
t214 = -t173 * t267 + t218;
t190 = t194 * qJDD(2);
t189 = t193 * qJDD(2);
t178 = t190 + t189;
t177 = t272 * t274;
t176 = t199 * t179;
t134 = t170 - 0.2e1 * t160;
t132 = t226 + 0.2e1 * t242;
t111 = -t141 + t162;
t110 = t140 - t162;
t103 = t141 - t140;
t100 = t133 + t242;
t99 = t133 - t242;
t92 = -t162 - t140;
t90 = -t144 * qJD(6) - t229;
t89 = -t140 - t141;
t88 = t205 * t101 - t208 * t226;
t87 = t205 * t170 - t208 * t99;
t86 = -t208 * t101 - t205 * t226;
t85 = -t208 * t170 - t205 * t99;
t74 = (t142 * t204 + t144 * t207) * t164;
t69 = t112 + t91;
t65 = (qJD(6) + t164) * t144 + t229;
t63 = -t144 * t254 - t204 * t91;
t62 = -t142 * t255 - t207 * t90;
t59 = -t204 * t110 - t259;
t58 = -t207 * t111 - t281;
t57 = -t204 * t237 - t259;
t56 = t207 * t237 - t262;
t55 = -t199 * t86 + t201 * t88;
t54 = -t199 * t85 + t201 * t87;
t53 = t207 * t92 - t281;
t52 = t204 * t92 + t280;
t48 = t230 * t173 + t218;
t44 = t204 * t69 + t207 * t223;
t43 = t204 * t65 - t207 * t275;
t42 = t204 * t223 - t207 * t69;
t41 = (t100 + t242) * pkin(4) + t214;
t40 = -pkin(4) * t242 - t214 - t295;
t37 = t205 * t56 + t208 * t275;
t36 = t205 * t275 - t208 * t56;
t35 = t205 * t52 + t208 * t65;
t34 = t205 * t65 - t208 * t52;
t32 = -qJ(5) * t270 + t39;
t31 = -pkin(4) * t270 + t38;
t30 = t205 * t42 + t208 * t89;
t29 = t205 * t89 - t208 * t42;
t23 = pkin(5) * t42 - qJ(5) * t44;
t20 = -t199 * t36 + t201 * t37;
t19 = -t199 * t34 + t201 * t35;
t17 = -t199 * t29 + t201 * t30;
t14 = pkin(5) * t275 + t266 * t57 - t263;
t13 = pkin(5) * t65 + t266 * t53 + t260;
t12 = t201 * t25 - t264;
t11 = pkin(5) * t56 - qJ(5) * t57 - t18;
t10 = pkin(5) * t52 - qJ(5) * t53 - t16;
t6 = t205 * t7 + t208 * t28;
t5 = t205 * t28 - t208 * t7;
t4 = pkin(5) * t89 + t266 * t44 - t8;
t3 = pkin(5) * t7 - qJ(5) * t8;
t2 = pkin(5) * t28 + t266 * t8;
t1 = -t199 * t5 + t201 * t6;
t15 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t243, 0, 0, 0, 0, 0, 0, (qJDD(2) * t209 - t206 * t211) * t200, (-qJDD(2) * t206 - t209 * t211) * t200, 0, t202 * t148 + (t114 * t206 - t209 * t227) * t200, 0, 0, 0, 0, 0, 0, (-t177 * t206 + t209 * t238) * t200, (t176 * t206 - t209 * t239) * t200, (t178 * t206 + t179 * t209) * t200, t202 * (t199 * t80 - t201 * t79) + (-t209 * t108 + t206 * t49) * t200, 0, 0, 0, 0, 0, 0, -t300 + (-t209 * t132 - t303) * t200, -t312 + (-t209 * t134 + t313) * t200, t202 * (t199 * t87 + t201 * t85) + (t206 * t54 - t279) * t200, t202 * (t199 * t25 + t201 * t24) + (t206 * t12 - t209 * t93) * t200, 0, 0, 0, 0, 0, 0, t202 * (t199 * t88 + t201 * t86) + (t206 * t55 - t279) * t200, t300 + (t209 * t100 + t303) * t200, t312 + (-t209 * t277 - t313) * t200, t202 * (t199 * t22 + t201 * t21) + (t206 * t9 - t209 * t48) * t200, 0, 0, 0, 0, 0, 0, t202 * (t199 * t35 + t201 * t34) + (t206 * t19 - t209 * t53) * t200, t202 * (t199 * t37 + t201 * t36) + (t206 * t20 - t209 * t57) * t200, t202 * (t199 * t30 + t201 * t29) + (t206 * t17 - t209 * t44) * t200, t202 * (t199 * t6 + t201 * t5) + (t206 * t1 - t209 * t8) * t200; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t227, -t114, 0, 0, t189, 0.2e1 * t199 * t238, 0, t190, 0, 0, -qJ(3) * t177 + t201 * t231, qJ(3) * t176 - t199 * t231, pkin(2) * t179 + qJ(3) * t178 + t49, -pkin(2) * t108 + qJ(3) * t49, t219, t199 * (-t208 * t132 - t205 * t134) + t201 * (-t205 * t132 + t208 * t134), t308, t216, -t305, t215, t199 * (t261 - t294) + t201 * (-pkin(3) * t132 - t258 - t293) - pkin(2) * t132 - t304, t199 * (t258 + t311) + t201 * (-pkin(3) * t134 + t261 - t310) - pkin(2) * t134 + t314, t199 * (-pkin(8) * t85 - t24) + t201 * (pkin(8) * t87 + t25 - t282) - t283 + qJ(3) * t54, -pkin(8) * t264 + t201 * (-pkin(3) * t93 + pkin(8) * t25) - pkin(2) * t93 + qJ(3) * t12, t215, -t308, t305, t219, t199 * (-t208 * t100 + t205 * t277) + t201 * (-t205 * t100 - t208 * t277), t216, t199 * (-pkin(8) * t86 - t205 * t31 + t208 * t32) + t201 * (pkin(8) * t88 + t205 * t32 + t208 * t31 - t282) - t283 + qJ(3) * t55, t199 * (-t205 * t41 + t294) + t201 * (t208 * t41 + t293) + t304 + (qJ(5) * t250 + t201 * t233 + pkin(2)) * t100, t199 * (t208 * t40 - t311) + t201 * (t205 * t40 + t310) - t314 - (-pkin(4) * t251 + t201 * (pkin(3) + t265) + pkin(2)) * t277, (t199 * (pkin(4) * t205 - qJ(5) * t208) + t201 * (-t233 - t265) - pkin(2)) * t48 + (qJ(3) + pkin(8)) * t9, t199 * (-t205 * t63 + t234) + t201 * (t208 * t63 + t235), t199 * (t208 * t103 - t205 * t43) + t201 * (t205 * t103 + t208 * t43), t199 * (-t205 * t58 + t208 * t69) + t201 * (t205 * t69 + t208 * t58), t199 * (-t205 * t62 - t234) + t201 * (t208 * t62 - t235), t199 * (-t205 * t59 + t208 * t223) + t201 * (t205 * t223 + t208 * t59), t199 * (t208 * t122 - t205 * t74) + t201 * (t205 * t122 + t208 * t74), t199 * (-pkin(8) * t34 + t208 * t10 - t205 * t13) + t201 * (-pkin(3) * t53 + pkin(8) * t35 + t205 * t10 + t208 * t13) - pkin(2) * t53 + qJ(3) * t19, t199 * (-pkin(8) * t36 + t208 * t11 - t205 * t14) + t201 * (-pkin(3) * t57 + pkin(8) * t37 + t205 * t11 + t208 * t14) - pkin(2) * t57 + qJ(3) * t20, t199 * (-pkin(8) * t29 - t205 * t4 + t208 * t23) + t201 * (-pkin(3) * t44 + pkin(8) * t30 + t205 * t23 + t208 * t4) - pkin(2) * t44 + qJ(3) * t17, t199 * (-pkin(8) * t5 - t205 * t2 + t208 * t3) + t201 * (-pkin(3) * t8 + pkin(8) * t6 + t208 * t2 + t205 * t3) - pkin(2) * t8 + qJ(3) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t238, t239, -t179, t108, 0, 0, 0, 0, 0, 0, t132, t134, t270, t93, 0, 0, 0, 0, 0, 0, t270, -t100, t277, t48, 0, 0, 0, 0, 0, 0, t53, t57, t44, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t253, t271, t170, -t253, -t226, qJDD(4), -t46, -t47, 0, 0, qJDD(4), -t101, t99, t253, t271, -t253, -pkin(4) * t101 - qJ(5) * t226, -pkin(4) * t286 - qJ(5) * t125 + t39, pkin(4) * t153 + (qJDD(4) + t285) * qJ(5) + t220, -pkin(4) * t39 + qJ(5) * t38, -t144 * t255 + t207 * t91, -t204 * t275 - t207 * t65, -t204 * t111 + t280, t142 * t254 - t204 * t90, t207 * t110 - t262, (-t142 * t207 + t144 * t204) * t164, qJ(5) * t65 + t266 * t52 + t263, qJ(5) * t275 + t266 * t56 + t260, qJ(5) * t89 + t266 * t42 - t7, qJ(5) * t28 + t266 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t101, t286, -t153, t39, 0, 0, 0, 0, 0, 0, t52, t56, t42, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, t103, t69, -t106, t223, t122, -t16, -t18, 0, 0;];
tauJ_reg  = t15;
