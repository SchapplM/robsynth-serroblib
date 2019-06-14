% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRPRR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 06:06
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRPRR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 06:04:49
% EndTime: 2019-05-05 06:05:02
% DurationCPUTime: 4.30s
% Computational Cost: add. (14949->365), mult. (30101->502), div. (0->0), fcn. (19893->12), ass. (0->226)
t211 = sin(qJ(3));
t198 = t211 ^ 2;
t218 = qJD(2) ^ 2;
t194 = t198 * t218;
t217 = qJD(3) ^ 2;
t180 = -t194 - t217;
t215 = cos(qJ(3));
t253 = t215 * t218;
t242 = t211 * t253;
t174 = qJDD(3) - t242;
t254 = t215 * t174;
t135 = t180 * t211 + t254;
t248 = qJD(2) * qJD(3);
t190 = t215 * t248;
t192 = t211 * qJDD(2);
t163 = t192 + 0.2e1 * t190;
t204 = sin(pkin(6));
t206 = cos(pkin(6));
t212 = sin(qJ(2));
t216 = cos(qJ(2));
t299 = (t135 * t212 + t163 * t216) * t204 + t206 * (t174 * t211 - t180 * t215);
t298 = pkin(8) * t135;
t199 = t215 ^ 2;
t195 = t199 * t218;
t182 = -t195 - t217;
t173 = qJDD(3) + t242;
t266 = t173 * t211;
t134 = -t182 * t215 + t266;
t241 = t211 * t248;
t246 = t215 * qJDD(2);
t166 = -0.2e1 * t241 + t246;
t296 = (t134 * t212 - t166 * t216) * t204 - t206 * (t173 * t215 + t182 * t211);
t295 = pkin(8) * t134;
t209 = sin(qJ(6));
t164 = t192 + t190;
t153 = qJDD(5) + t164;
t147 = qJDD(6) + t153;
t210 = sin(qJ(5));
t214 = cos(qJ(5));
t250 = qJD(2) * t215;
t158 = qJD(3) * t210 + t214 * t250;
t160 = qJD(3) * t214 - t210 * t250;
t213 = cos(qJ(6));
t124 = t213 * t158 + t160 * t209;
t126 = -t158 * t209 + t160 * t213;
t96 = t126 * t124;
t287 = -t96 + t147;
t294 = t209 * t287;
t130 = t160 * t158;
t284 = -t130 + t153;
t293 = t210 * t284;
t292 = t213 * t287;
t291 = t214 * t284;
t290 = t164 + t190;
t203 = sin(pkin(11));
t205 = cos(pkin(11));
t171 = g(1) * t203 - g(2) * t205;
t200 = -g(3) + qJDD(1);
t289 = t171 * t206 + t200 * t204;
t165 = -t241 + t246;
t119 = -t158 * qJD(5) + t214 * qJDD(3) - t210 * t165;
t251 = qJD(2) * t211;
t186 = qJD(5) + t251;
t142 = t186 * t158;
t283 = -t142 + t119;
t103 = t142 + t119;
t172 = -g(1) * t205 - g(2) * t203;
t115 = t216 * t172 + t212 * t289;
t106 = -t218 * pkin(2) + qJDD(2) * pkin(8) + t115;
t257 = t211 * qJ(4);
t233 = -pkin(3) * t215 - t257;
t238 = t218 * t233 + t106;
t282 = -t217 * pkin(3) + t238 * t215;
t281 = t254 + (t195 - t217) * t211;
t122 = t124 ^ 2;
t123 = t126 ^ 2;
t151 = t158 ^ 2;
t152 = t160 ^ 2;
t178 = qJD(6) + t186;
t177 = t178 ^ 2;
t183 = t186 ^ 2;
t280 = 2 * qJD(4);
t279 = -pkin(3) - pkin(9);
t141 = -t171 * t204 + t200 * t206;
t137 = t215 * t141;
t235 = -qJDD(3) * pkin(3) - t217 * qJ(4) + qJDD(4) - t137;
t72 = -qJDD(3) * pkin(9) + (t164 - t190) * pkin(4) + (-pkin(9) * t253 + t238) * t211 + t235;
t176 = pkin(4) * t251 - qJD(3) * pkin(9);
t234 = t212 * t172 - t216 * t289;
t105 = -qJDD(2) * pkin(2) - t218 * pkin(8) + t234;
t221 = -t165 * pkin(3) - qJ(4) * t290 + t105;
t240 = pkin(3) * qJD(3) - (2 * qJD(4));
t73 = -pkin(4) * t195 - t165 * pkin(9) + (-t176 + t240) * t251 + t221;
t40 = t210 * t73 - t214 * t72;
t27 = pkin(5) * t284 - pkin(10) * t103 - t40;
t237 = t210 * qJDD(3) + t214 * t165;
t118 = -qJD(5) * t160 - t237;
t138 = pkin(5) * t186 - pkin(10) * t160;
t41 = t210 * t72 + t214 * t73;
t28 = -pkin(5) * t151 + pkin(10) * t118 - t138 * t186 + t41;
t13 = t209 * t28 - t213 * t27;
t14 = t209 * t27 + t213 * t28;
t9 = -t13 * t213 + t14 * t209;
t278 = t210 * t9;
t277 = t214 * t9;
t247 = qJDD(3) * qJ(4);
t256 = t211 * t141;
t71 = t247 + t256 + t165 * pkin(4) - pkin(9) * t195 + (t280 + t176) * qJD(3) + t282;
t42 = -t118 * pkin(5) - t151 * pkin(10) + t160 * t138 + t71;
t276 = t209 * t42;
t85 = t96 + t147;
t275 = t209 * t85;
t274 = t210 * t71;
t222 = (-qJD(5) + t186) * t160 - t237;
t74 = -t103 * t214 + t210 * t222;
t273 = t211 * t74;
t272 = t213 * t42;
t271 = t213 * t85;
t270 = t214 * t71;
t265 = t178 * t209;
t264 = t178 * t213;
t261 = t186 * t210;
t260 = t186 * t214;
t113 = t130 + t153;
t258 = t210 * t113;
t255 = t214 * t113;
t168 = (t198 + t199) * qJDD(2);
t169 = t194 + t195;
t252 = pkin(2) * t169 + pkin(8) * t168;
t249 = qJD(6) + t178;
t245 = -t152 - t183;
t244 = t211 * t96;
t243 = t211 * t130;
t10 = t13 * t209 + t213 * t14;
t93 = t211 * t106 - t137;
t94 = t215 * t106 + t256;
t59 = t211 * t93 + t215 * t94;
t239 = -t213 * t118 + t209 * t119;
t90 = -t177 - t122;
t47 = t209 * t90 + t292;
t236 = pkin(5) * t47 - t13;
t19 = t210 * t41 - t214 * t40;
t15 = t19 * t211 + t215 * t71;
t20 = t210 * t40 + t214 * t41;
t232 = t209 * t118 + t213 * t119;
t231 = (-t194 + t217) * t215 + t266;
t227 = pkin(2) - t233;
t104 = -t123 - t177;
t63 = t104 * t213 - t275;
t226 = pkin(5) * t63 - t14;
t224 = t215 * t279 - pkin(2) - t257;
t223 = (-qJD(6) + t178) * t126 - t239;
t77 = -qJD(6) * t124 + t232;
t220 = qJD(3) * t280 + t282;
t80 = t238 * t211 + t235;
t219 = t220 + t247;
t81 = t240 * t251 + t221;
t170 = t194 - t195;
t140 = -t152 + t183;
t139 = t151 - t183;
t132 = t290 * t211;
t131 = (t165 - t241) * t215;
t129 = t152 - t151;
t128 = t163 * t215 + t166 * t211;
t121 = (t168 * t212 + t169 * t216) * t204;
t120 = -t183 - t151;
t111 = -t151 - t152;
t110 = t178 * t124;
t108 = -t123 + t177;
t107 = t122 - t177;
t98 = (qJD(5) + t186) * t160 + t237;
t95 = -t122 + t123;
t92 = -t210 * t245 - t255;
t91 = t214 * t245 - t258;
t88 = t214 * t120 - t293;
t87 = t210 * t120 + t291;
t83 = (-t124 * t213 + t126 * t209) * t178;
t82 = (-t124 * t209 - t126 * t213) * t178;
t79 = -t122 - t123;
t78 = t219 + t256;
t76 = -qJD(6) * t126 - t239;
t75 = t210 * t103 + t214 * t222;
t70 = t107 * t213 - t275;
t69 = -t108 * t209 + t292;
t68 = t107 * t209 + t271;
t67 = t108 * t213 + t294;
t65 = t211 * t91 + t215 * t283;
t64 = -t104 * t209 - t271;
t60 = t211 * t87 + t215 * t98;
t58 = -t249 * t124 + t232;
t57 = t110 + t77;
t56 = -t110 + t77;
t53 = t249 * t126 + t239;
t52 = -t126 * t265 + t213 * t77;
t51 = t126 * t264 + t209 * t77;
t50 = t124 * t264 - t209 * t76;
t49 = t124 * t265 + t213 * t76;
t48 = t213 * t90 - t294;
t45 = t111 * t215 + t273;
t43 = t211 * t80 + t215 * t78;
t38 = -t210 * t63 + t214 * t64;
t37 = t210 * t64 + t214 * t63;
t36 = -t209 * t56 - t213 * t53;
t35 = t209 * t57 + t213 * t223;
t34 = -t209 * t53 + t213 * t56;
t33 = t209 * t223 - t213 * t57;
t32 = pkin(5) * t33;
t30 = -t210 * t47 + t214 * t48;
t29 = t210 * t48 + t214 * t47;
t26 = -pkin(10) * t63 + t272;
t24 = -pkin(10) * t47 + t276;
t23 = t211 * t37 + t215 * t58;
t22 = t211 * t29 + t215 * t53;
t21 = -pkin(5) * t58 + pkin(10) * t64 + t276;
t18 = -pkin(5) * t53 + pkin(10) * t48 - t272;
t17 = -t210 * t33 + t214 * t35;
t16 = t210 * t35 + t214 * t33;
t11 = t16 * t211 + t215 * t79;
t8 = pkin(5) * t9;
t6 = -pkin(5) * t42 + pkin(10) * t10;
t5 = -pkin(10) * t33 - t9;
t4 = -pkin(5) * t79 + pkin(10) * t35 + t10;
t3 = t10 * t214 - t278;
t2 = t10 * t210 + t277;
t1 = t2 * t211 + t215 * t42;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t200, 0, 0, 0, 0, 0, 0, (qJDD(2) * t216 - t212 * t218) * t204, (-qJDD(2) * t212 - t216 * t218) * t204, 0, t206 * t141 + (t115 * t212 - t216 * t234) * t204, 0, 0, 0, 0, 0, 0, -t296, -t299, t121, t206 * (t211 * t94 - t215 * t93) + (-t105 * t216 + t212 * t59) * t204, 0, 0, 0, 0, 0, 0, t121, t296, t299, t206 * (t211 * t78 - t215 * t80) + (t212 * t43 - t216 * t81) * t204, 0, 0, 0, 0, 0, 0, t206 * (t211 * t98 - t215 * t87) + (t212 * t60 - t216 * t88) * t204, t206 * (t211 * t283 - t215 * t91) + (t212 * t65 - t216 * t92) * t204, t206 * (t111 * t211 - t215 * t74) + (t212 * t45 - t216 * t75) * t204, t206 * (-t19 * t215 + t211 * t71) + (t15 * t212 - t20 * t216) * t204, 0, 0, 0, 0, 0, 0, t206 * (t211 * t53 - t215 * t29) + (t212 * t22 - t216 * t30) * t204, t206 * (t211 * t58 - t215 * t37) + (t212 * t23 - t216 * t38) * t204, t206 * (-t16 * t215 + t211 * t79) + (t11 * t212 - t17 * t216) * t204, t206 * (-t2 * t215 + t211 * t42) + (t1 * t212 - t216 * t3) * t204; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t234, -t115, 0, 0, t132, t128, t231, t131, t281, 0, pkin(2) * t166 - t105 * t215 - t295, -pkin(2) * t163 + t105 * t211 - t298, t59 + t252, -pkin(2) * t105 + pkin(8) * t59, 0, -t231, -t281, t132, t128, t131, (pkin(3) * t169 + t219) * t215 + (qJ(4) * t169 + t137 + t80) * t211 + t252, -t227 * t166 + t215 * t81 + t295, t211 * (-pkin(3) * t241 + t251 * t280 - t221) + t298 + t227 * t163, pkin(8) * t43 - t227 * t81, t243 + t215 * (-t119 * t210 - t160 * t260), t211 * t129 + t215 * (t210 * t98 - t214 * t283), t211 * t103 + t215 * (-t140 * t214 - t293), -t243 + t215 * (-t118 * t214 - t158 * t261), t211 * t222 + t215 * (-t139 * t210 - t255), t211 * t153 + t215 * (t158 * t210 + t160 * t214) * t186, t211 * (pkin(4) * t87 - t40) + t215 * (pkin(4) * t98 + t270) + pkin(8) * t60 + t224 * t88, t211 * (pkin(4) * t91 - t41) + t215 * (pkin(4) * t283 - t274) + pkin(8) * t65 + t224 * t92, pkin(4) * t273 + t215 * (pkin(4) * t111 - t20) + pkin(8) * t45 + t224 * t75, t224 * t20 + (pkin(4) + pkin(8)) * t15, t244 + t215 * (-t210 * t52 - t214 * t51), t211 * t95 + t215 * (-t210 * t36 - t214 * t34), t211 * t57 + t215 * (-t210 * t69 - t214 * t67), -t244 + t215 * (-t210 * t50 - t214 * t49), t211 * t223 + t215 * (-t210 * t70 - t214 * t68), t211 * t147 + t215 * (-t210 * t83 - t214 * t82), t211 * (pkin(4) * t29 + t236) + t215 * (pkin(4) * t53 - t214 * t18 - t210 * t24) + pkin(8) * t22 + t224 * t30, t211 * (pkin(4) * t37 + t226) + t215 * (pkin(4) * t58 - t214 * t21 - t210 * t26) + pkin(8) * t23 + t224 * t38, t211 * (pkin(4) * t16 + t32) + t215 * (pkin(4) * t79 - t210 * t5 - t214 * t4) + pkin(8) * t11 + t224 * t17, t211 * (pkin(4) * t2 + t8) + t215 * (pkin(4) * t42 + pkin(10) * t278 - t214 * t6) + pkin(8) * t1 + t224 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t242, t170, t192, t242, t246, qJDD(3), -t93, -t94, 0, 0, qJDD(3), -t192, -t246, -t242, t170, t242, (-pkin(3) * t211 + qJ(4) * t215) * qJDD(2), -pkin(3) * t173 - qJ(4) * t182 + t80, -pkin(3) * t180 + t256 + (qJDD(3) + t174) * qJ(4) + t220, -pkin(3) * t80 + qJ(4) * t78, t119 * t214 - t160 * t261, -t210 * t283 - t214 * t98, -t140 * t210 + t291, -t118 * t210 + t158 * t260, t139 * t214 - t258, (-t158 * t214 + t160 * t210) * t186, qJ(4) * t98 + t279 * t87 + t274, qJ(4) * t283 + t279 * t91 + t270, qJ(4) * t111 + t279 * t74 - t19, qJ(4) * t71 + t279 * t19, -t210 * t51 + t214 * t52, -t210 * t34 + t214 * t36, -t210 * t67 + t214 * t69, -t210 * t49 + t214 * t50, -t210 * t68 + t214 * t70, -t210 * t82 + t214 * t83, qJ(4) * t53 - t210 * t18 + t214 * t24 + t279 * t29, qJ(4) * t58 - t210 * t21 + t214 * t26 + t279 * t37, qJ(4) * t79 + t279 * t16 - t210 * t4 + t214 * t5, -pkin(10) * t277 + qJ(4) * t42 + t279 * t2 - t210 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t192, t173, t180, t80, 0, 0, 0, 0, 0, 0, t87, t91, t74, t19, 0, 0, 0, 0, 0, 0, t29, t37, t16, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t129, t103, -t130, t222, t153, -t40, -t41, 0, 0, t96, t95, t57, -t96, t223, t147, t236, t226, t32, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t95, t57, -t96, t223, t147, -t13, -t14, 0, 0;];
tauJ_reg  = t7;
