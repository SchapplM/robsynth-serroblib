% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PPRRRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
% 
% Output:
% tauJ_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 20:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PPRRRR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:40:38
% EndTime: 2019-05-04 20:40:46
% DurationCPUTime: 4.09s
% Computational Cost: add. (21237->378), mult. (39756->582), div. (0->0), fcn. (32265->16), ass. (0->246)
t200 = sin(qJ(6));
t202 = sin(qJ(4));
t205 = cos(qJ(5));
t201 = sin(qJ(5));
t206 = cos(qJ(4));
t263 = t201 * t206;
t165 = (t202 * t205 + t263) * qJD(3);
t189 = qJD(4) + qJD(5);
t204 = cos(qJ(6));
t145 = t165 * t200 - t204 * t189;
t147 = t165 * t204 + t189 * t200;
t124 = t147 * t145;
t254 = qJD(3) * qJD(4);
t245 = t206 * t254;
t253 = t202 * qJDD(3);
t169 = t245 + t253;
t246 = t202 * t254;
t252 = t206 * qJDD(3);
t223 = -t246 + t252;
t241 = t201 * t169 - t205 * t223;
t126 = -qJD(5) * t165 - t241;
t125 = qJDD(6) - t126;
t283 = -t124 + t125;
t288 = t200 * t283;
t257 = qJD(3) * t202;
t163 = -t205 * t206 * qJD(3) + t201 * t257;
t141 = t165 * t163;
t188 = qJDD(4) + qJDD(5);
t282 = -t141 + t188;
t287 = t201 * t282;
t286 = t204 * t283;
t285 = t205 * t282;
t194 = sin(pkin(7));
t197 = cos(pkin(7));
t193 = sin(pkin(13));
t196 = cos(pkin(13));
t195 = sin(pkin(6));
t198 = cos(pkin(6));
t273 = sin(pkin(12));
t274 = cos(pkin(12));
t221 = g(1) * t273 - g(2) * t274;
t258 = -g(3) + qJDD(1);
t219 = t195 * t258 + t198 * t221;
t222 = -g(1) * t274 - g(2) * t273;
t216 = -t193 * t222 + t196 * t219;
t218 = -t195 * t221 + t198 * t258 + qJDD(2);
t284 = t194 * t218 + t197 * t216;
t138 = pkin(5) * t163 - pkin(11) * t165;
t209 = qJD(3) ^ 2;
t261 = t202 * t209;
t178 = t206 * t261;
t133 = t193 * t219 + t196 * t222;
t203 = sin(qJ(3));
t207 = cos(qJ(3));
t98 = t207 * t133 + t284 * t203;
t214 = -t209 * pkin(3) + qJDD(3) * pkin(9) + t98;
t213 = t202 * t214;
t212 = qJDD(4) * pkin(4) - t169 * pkin(10) - t213;
t113 = -t194 * t216 + t197 * t218;
t260 = t206 * t113;
t211 = pkin(4) * t178 + pkin(10) * t245 + t212 + t260;
t175 = qJD(4) * pkin(4) - pkin(10) * t257;
t191 = t206 ^ 2;
t187 = t191 * t209;
t66 = t202 * t113 + t206 * t214;
t63 = -pkin(4) * t187 + pkin(10) * t223 - qJD(4) * t175 + t66;
t276 = t205 * t63;
t210 = -t201 * t211 - t276;
t281 = t189 ^ 2;
t38 = -pkin(5) * t281 + t188 * pkin(11) - t163 * t138 - t210;
t127 = -t163 * qJD(5) + t205 * t169 + t201 * t223;
t156 = t189 * t163;
t118 = -t156 + t127;
t240 = t203 * t133 - t284 * t207;
t96 = -qJDD(3) * pkin(3) - t209 * pkin(9) + t240;
t75 = -t223 * pkin(4) - pkin(10) * t187 + t175 * t257 + t96;
t54 = -t118 * pkin(11) + (t165 * t189 - t126) * pkin(5) + t75;
t23 = t200 * t38 - t204 * t54;
t24 = t200 * t54 + t204 * t38;
t13 = t200 * t23 + t204 * t24;
t158 = qJD(6) + t163;
t242 = t200 * t127 - t204 * t188;
t90 = (qJD(6) - t158) * t147 + t242;
t143 = t145 ^ 2;
t144 = t147 ^ 2;
t157 = t158 ^ 2;
t161 = t163 ^ 2;
t162 = t165 ^ 2;
t280 = pkin(5) * t201;
t40 = t201 * t63 - t205 * t211;
t37 = -t188 * pkin(5) - pkin(11) * t281 + t138 * t165 + t40;
t279 = -pkin(5) * t37 + pkin(11) * t13;
t41 = t276 + t201 * t212 + (pkin(4) * t261 + pkin(10) * t254 + t113) * t263;
t19 = t201 * t41 - t205 * t40;
t278 = t19 * t202;
t34 = t200 * t37;
t277 = t201 * t75;
t35 = t204 * t37;
t275 = t205 * t75;
t100 = t124 + t125;
t272 = t100 * t200;
t271 = t100 * t204;
t136 = t141 + t188;
t270 = t136 * t201;
t269 = t136 * t205;
t268 = t158 * t200;
t267 = t158 * t204;
t266 = t189 * t201;
t265 = t189 * t205;
t264 = t196 * t197;
t173 = qJDD(4) + t178;
t262 = t202 * t173;
t174 = qJDD(4) - t178;
t259 = t206 * t174;
t255 = qJD(6) + t158;
t122 = -t144 - t157;
t72 = -t122 * t200 - t271;
t229 = -t204 * t127 - t200 * t188;
t95 = t145 * t255 + t229;
t251 = pkin(5) * t95 + pkin(11) * t72 + t34;
t112 = -t157 - t143;
t69 = t112 * t204 - t288;
t92 = -t147 * t255 - t242;
t250 = pkin(5) * t92 + pkin(11) * t69 - t35;
t249 = t201 * t124;
t248 = t205 * t124;
t247 = -pkin(5) * t205 - pkin(4);
t20 = t201 * t40 + t205 * t41;
t65 = t213 - t260;
t43 = t202 * t65 + t206 * t66;
t107 = t143 + t144;
t105 = -qJD(6) * t145 - t229;
t132 = t158 * t145;
t94 = t105 + t132;
t60 = t200 * t94 - t204 * t90;
t244 = pkin(5) * t107 + pkin(11) * t60 + t13;
t5 = t13 * t201 - t205 * t37;
t6 = t13 * t205 + t201 * t37;
t3 = -t202 * t5 + t206 * t6;
t12 = t200 * t24 - t204 * t23;
t239 = -t12 * t207 + t203 * t3;
t9 = t20 * t206 - t278;
t238 = t203 * t9 - t207 * t75;
t46 = t107 * t205 + t201 * t60;
t47 = -t107 * t201 + t205 * t60;
t27 = -t202 * t46 + t206 * t47;
t58 = -t200 * t90 - t204 * t94;
t237 = t203 * t27 - t207 * t58;
t49 = t201 * t69 + t205 * t92;
t50 = -t201 * t92 + t205 * t69;
t31 = -t202 * t49 + t206 * t50;
t68 = t112 * t200 + t286;
t236 = t203 * t31 - t207 * t68;
t51 = t201 * t72 + t205 * t95;
t52 = -t201 * t95 + t205 * t72;
t33 = -t202 * t51 + t206 * t52;
t71 = t122 * t204 - t272;
t235 = t203 * t33 - t207 * t71;
t234 = t203 * t43 - t207 * t96;
t233 = t203 * t98 - t207 * t240;
t114 = (qJD(5) + t189) * t165 + t241;
t134 = -t281 - t161;
t108 = t134 * t201 + t285;
t109 = t134 * t205 - t287;
t74 = -t108 * t202 + t109 * t206;
t232 = -t114 * t207 + t203 * t74;
t152 = -t162 - t281;
t120 = t152 * t205 - t270;
t121 = -t152 * t201 - t269;
t83 = -t120 * t202 + t121 * t206;
t231 = -t118 * t207 + t203 * t83;
t128 = -t161 - t162;
t119 = t156 + t127;
t220 = (-qJD(5) + t189) * t165 - t241;
t80 = -t119 * t205 + t201 * t220;
t81 = t119 * t201 + t205 * t220;
t56 = -t202 * t80 + t206 * t81;
t230 = -t128 * t207 + t203 * t56;
t208 = qJD(4) ^ 2;
t177 = -t187 - t208;
t150 = t177 * t206 - t262;
t170 = -0.2e1 * t246 + t252;
t228 = t150 * t203 + t170 * t207;
t190 = t202 ^ 2;
t186 = t190 * t209;
t176 = -t186 - t208;
t151 = -t176 * t202 - t259;
t168 = 0.2e1 * t245 + t253;
t227 = t151 * t203 - t168 * t207;
t171 = (t190 + t191) * qJDD(3);
t172 = t186 + t187;
t226 = t171 * t203 + t172 * t207;
t225 = qJDD(3) * t207 - t203 * t209;
t224 = -qJDD(3) * t203 - t207 * t209;
t160 = t225 * t194;
t159 = t224 * t194;
t154 = -t162 + t281;
t153 = t161 - t281;
t149 = -t174 * t202 + t176 * t206;
t148 = t173 * t206 + t177 * t202;
t140 = t162 - t161;
t139 = t226 * t194;
t131 = -t144 + t157;
t130 = t143 - t157;
t123 = t144 - t143;
t111 = t197 * t149 + t194 * t227;
t110 = t197 * t148 + t194 * t228;
t104 = -qJD(6) * t147 - t242;
t103 = (-t145 * t204 + t147 * t200) * t158;
t102 = (-t145 * t200 - t147 * t204) * t158;
t93 = t105 - t132;
t87 = t105 * t204 - t147 * t268;
t86 = t105 * t200 + t147 * t267;
t85 = -t104 * t200 + t145 * t267;
t84 = t104 * t204 + t145 * t268;
t82 = t120 * t206 + t121 * t202;
t79 = t130 * t204 - t272;
t78 = -t131 * t200 + t286;
t77 = t130 * t200 + t271;
t76 = t131 * t204 + t288;
t73 = t108 * t206 + t109 * t202;
t61 = -t200 * t93 + t204 * t92;
t59 = t200 * t92 + t204 * t93;
t55 = t202 * t81 + t206 * t80;
t48 = t197 * t113 + t194 * t233;
t45 = t194 * t231 + t197 * t82;
t44 = t194 * t232 + t197 * t73;
t42 = t202 * t66 - t206 * t65;
t32 = t202 * t52 + t206 * t51;
t30 = t202 * t50 + t206 * t49;
t29 = t194 * t230 + t197 * t55;
t28 = -pkin(11) * t71 + t35;
t26 = t202 * t47 + t206 * t46;
t25 = -pkin(11) * t68 + t34;
t18 = t194 * t234 + t197 * t42;
t17 = -pkin(5) * t71 + t24;
t16 = -pkin(5) * t68 + t23;
t15 = t194 * t235 + t197 * t32;
t14 = t194 * t236 + t197 * t30;
t10 = t194 * t237 + t197 * t26;
t8 = t19 * t206 + t20 * t202;
t7 = -pkin(11) * t58 - t12;
t4 = t194 * t238 + t197 * t8;
t2 = t202 * t6 + t206 * t5;
t1 = t194 * t239 + t197 * t2;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t258, 0, 0, 0, 0, 0, 0, 0, 0, 0, t198 * t218 + (t193 * t133 + t196 * t216) * t195, 0, 0, 0, 0, 0, 0, t198 * t160 + (t193 * t224 + t225 * t264) * t195, t198 * t159 + (-t193 * t225 + t224 * t264) * t195, 0, t198 * t48 + (t193 * (t203 * t240 + t207 * t98) + t196 * (-t194 * t113 + t197 * t233)) * t195, 0, 0, 0, 0, 0, 0, t198 * t110 + (t193 * (t150 * t207 - t170 * t203) + t196 * (-t194 * t148 + t197 * t228)) * t195, t198 * t111 + (t193 * (t151 * t207 + t168 * t203) + t196 * (-t194 * t149 + t197 * t227)) * t195, t198 * t139 + (t193 * (t171 * t207 - t172 * t203) + t226 * t264) * t195, t198 * t18 + (t193 * (t203 * t96 + t207 * t43) + t196 * (-t194 * t42 + t197 * t234)) * t195, 0, 0, 0, 0, 0, 0, t198 * t44 + (t193 * (t114 * t203 + t207 * t74) + t196 * (-t194 * t73 + t197 * t232)) * t195, t198 * t45 + (t193 * (t118 * t203 + t207 * t83) + t196 * (-t194 * t82 + t197 * t231)) * t195, t198 * t29 + (t193 * (t128 * t203 + t207 * t56) + t196 * (-t194 * t55 + t197 * t230)) * t195, t198 * t4 + (t193 * (t203 * t75 + t207 * t9) + t196 * (-t194 * t8 + t197 * t238)) * t195, 0, 0, 0, 0, 0, 0, t198 * t14 + (t193 * (t203 * t68 + t207 * t31) + t196 * (-t194 * t30 + t197 * t236)) * t195, t198 * t15 + (t193 * (t203 * t71 + t207 * t33) + t196 * (-t194 * t32 + t197 * t235)) * t195, t198 * t10 + (t193 * (t203 * t58 + t207 * t27) + t196 * (-t194 * t26 + t197 * t237)) * t195, t198 * t1 + (t193 * (t12 * t203 + t207 * t3) + t196 * (-t194 * t2 + t197 * t239)) * t195; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t218, 0, 0, 0, 0, 0, 0, t160, t159, 0, t48, 0, 0, 0, 0, 0, 0, t110, t111, t139, t18, 0, 0, 0, 0, 0, 0, t44, t45, t29, t4, 0, 0, 0, 0, 0, 0, t14, t15, t10, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3), -t240, -t98, 0, 0, (t169 + t245) * t202, t168 * t206 + t170 * t202, t262 + t206 * (-t186 + t208), t170 * t206, t202 * (t187 - t208) + t259, 0, pkin(3) * t170 + pkin(9) * t150 - t206 * t96, -pkin(3) * t168 + pkin(9) * t151 + t202 * t96, pkin(3) * t172 + pkin(9) * t171 + t43, -pkin(3) * t96 + pkin(9) * t43, t202 * (t127 * t205 - t165 * t266) + t206 * (t127 * t201 + t165 * t265), t202 * (-t114 * t205 - t118 * t201) + t206 * (-t114 * t201 + t118 * t205), t202 * (-t154 * t201 + t285) + t206 * (t154 * t205 + t287), t202 * (-t126 * t201 + t163 * t265) + t206 * (t126 * t205 + t163 * t266), t202 * (t153 * t205 - t270) + t206 * (t153 * t201 + t269), (t202 * (-t163 * t205 + t165 * t201) + t206 * (-t163 * t201 - t165 * t205)) * t189, t202 * (-pkin(10) * t108 + t277) + t206 * (-pkin(4) * t114 + pkin(10) * t109 - t275) - pkin(3) * t114 + pkin(9) * t74, t202 * (-pkin(10) * t120 + t275) + t206 * (-pkin(4) * t118 + pkin(10) * t121 + t277) - pkin(3) * t118 + pkin(9) * t83, t202 * (-pkin(10) * t80 - t19) + t206 * (-pkin(4) * t128 + pkin(10) * t81 + t20) - pkin(3) * t128 + pkin(9) * t56, -pkin(10) * t278 + t206 * (-pkin(4) * t75 + pkin(10) * t20) - pkin(3) * t75 + pkin(9) * t9, t202 * (t205 * t87 + t249) + t206 * (t201 * t87 - t248), t202 * (t123 * t201 + t205 * t61) + t206 * (-t123 * t205 + t201 * t61), t202 * (t201 * t94 + t205 * t78) + t206 * (t201 * t78 - t205 * t94), t202 * (t205 * t85 - t249) + t206 * (t201 * t85 + t248), t202 * (-t201 * t90 + t205 * t79) + t206 * (t201 * t79 + t205 * t90), t202 * (t103 * t205 + t125 * t201) + t206 * (t103 * t201 - t125 * t205), t202 * (-pkin(10) * t49 - t16 * t201 + t205 * t25) + t206 * (-pkin(4) * t68 + pkin(10) * t50 + t16 * t205 + t201 * t25) - pkin(3) * t68 + pkin(9) * t31, t202 * (-pkin(10) * t51 - t17 * t201 + t205 * t28) + t206 * (-pkin(4) * t71 + pkin(10) * t52 + t17 * t205 + t201 * t28) - pkin(3) * t71 + pkin(9) * t33, t202 * (-pkin(10) * t46 + t205 * t7) + t206 * (pkin(10) * t47 + t201 * t7) + pkin(9) * t27 + (t202 * t280 + t206 * t247 - pkin(3)) * t58, (t202 * (-pkin(11) * t205 + t280) + t206 * (-pkin(11) * t201 + t247) - pkin(3)) * t12 + (pkin(9) + pkin(10)) * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t178, t186 - t187, t253, t178, t252, qJDD(4), -t65, -t66, 0, 0, t141, t140, t119, -t141, t220, t188, pkin(4) * t108 - t40, pkin(4) * t120 + t210, pkin(4) * t80, pkin(4) * t19, t86, t59, t76, t84, t77, t102, pkin(4) * t49 + t250, pkin(4) * t51 + t251, pkin(4) * t46 + t244, pkin(4) * t5 + t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, t140, t119, -t141, t220, t188, -t40, -t41, 0, 0, t86, t59, t76, t84, t77, t102, t250, t251, t244, t279; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t124, t123, t94, -t124, -t90, t125, -t23, -t24, 0, 0;];
tauJ_reg  = t11;
