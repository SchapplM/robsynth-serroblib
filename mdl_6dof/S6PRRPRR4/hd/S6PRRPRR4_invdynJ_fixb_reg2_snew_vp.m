% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S6PRRPRR4
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
% Datum: 2019-05-05 05:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S6PRRPRR4_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR4_invdynJ_fixb_reg2_snew_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR4_invdynJ_fixb_reg2_snew_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRPRR4_invdynJ_fixb_reg2_snew_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRR4_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR4_invdynJ_fixb_reg2_snew_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 05:10:52
% EndTime: 2019-05-05 05:11:04
% DurationCPUTime: 4.15s
% Computational Cost: add. (12755->373), mult. (25671->499), div. (0->0), fcn. (17314->12), ass. (0->225)
t189 = sin(qJ(5));
t178 = qJDD(3) - qJDD(5);
t190 = sin(qJ(3));
t193 = cos(qJ(5));
t194 = cos(qJ(3));
t140 = (t190 * t189 + t194 * t193) * qJD(2);
t244 = t190 * qJD(2);
t142 = -t189 * t194 * qJD(2) + t193 * t244;
t259 = t142 * t140;
t285 = -t178 - t259;
t289 = t189 * t285;
t288 = t193 * t285;
t196 = qJD(3) ^ 2;
t181 = t190 ^ 2;
t197 = qJD(2) ^ 2;
t258 = t181 * t197;
t163 = t196 + t258;
t168 = t190 * t197 * t194;
t161 = qJDD(3) - t168;
t250 = t194 * t161;
t125 = -t190 * t163 + t250;
t241 = qJD(2) * qJD(3);
t233 = t194 * t241;
t240 = t190 * qJDD(2);
t151 = 0.2e1 * t233 + t240;
t185 = sin(pkin(6));
t186 = cos(pkin(6));
t191 = sin(qJ(2));
t195 = cos(qJ(2));
t287 = t186 * (t190 * t161 + t194 * t163) + (t191 * t125 + t195 * t151) * t185;
t286 = pkin(8) * t125;
t238 = qJD(3) - qJD(5);
t225 = t140 * t238;
t152 = t233 + t240;
t174 = t190 * t241;
t239 = t194 * qJDD(2);
t153 = -t174 + t239;
t98 = -t140 * qJD(5) + t193 * t152 - t189 * t153;
t88 = t98 + t225;
t188 = sin(qJ(6));
t192 = cos(qJ(6));
t117 = t188 * t142 + t192 * t238;
t119 = t192 * t142 - t188 * t238;
t94 = t119 * t117;
t228 = t189 * t152 + t193 * t153;
t97 = -t142 * qJD(5) - t228;
t95 = qJDD(6) - t97;
t278 = -t94 + t95;
t283 = t188 * t278;
t282 = t192 * t278;
t281 = t152 + t233;
t263 = sin(pkin(11));
t264 = cos(pkin(11));
t205 = t263 * g(1) - t264 * g(2);
t203 = t186 * t205;
t248 = -g(3) + qJDD(1);
t280 = t185 * t248 + t203;
t236 = t238 ^ 2;
t279 = pkin(8) - pkin(9);
t182 = t194 ^ 2;
t257 = t182 * t197;
t275 = t250 + t190 * (-t196 + t257);
t211 = -qJD(3) * pkin(4) - pkin(9) * t244;
t274 = -t153 * pkin(4) + pkin(9) * t257 - t211 * t244;
t136 = qJD(6) + t140;
t230 = t192 * t178 + t188 * t98;
t58 = (qJD(6) - t136) * t119 + t230;
t262 = qJ(4) * t190;
t272 = pkin(3) + pkin(4);
t207 = t194 * t272 + pkin(2) + t262;
t115 = t117 ^ 2;
t116 = t119 ^ 2;
t134 = t136 ^ 2;
t138 = t140 ^ 2;
t139 = t142 ^ 2;
t273 = 2 * qJD(4);
t160 = qJDD(3) + t168;
t159 = -t264 * g(1) - t263 * g(2);
t106 = t195 * t159 + t280 * t191;
t101 = -t197 * pkin(2) + qJDD(2) * pkin(8) + t106;
t200 = -t185 * t205 + t186 * t248;
t81 = t190 * t101 - t194 * t200;
t206 = -qJDD(3) * pkin(3) - t196 * qJ(4) + qJDD(4) + t81;
t217 = -pkin(3) * t194 - t262;
t198 = -t152 * pkin(9) - t160 * pkin(4) + (qJD(3) * t194 * pkin(9) + t217 * t244) * qJD(2) + t206;
t210 = t197 * t217;
t82 = t194 * t101 + t190 * t200;
t218 = qJDD(3) * qJ(4) + qJD(3) * t273 + t194 * t210 + t82;
t69 = -t196 * pkin(3) + t218;
t64 = -pkin(4) * t257 - t153 * pkin(9) + qJD(3) * t211 + t69;
t35 = t189 * t198 + t193 * t64;
t111 = t140 * pkin(5) - t142 * pkin(10);
t34 = t189 * t64 - t193 * t198;
t22 = t178 * pkin(5) - t236 * pkin(10) + t142 * t111 + t34;
t270 = t188 * t22;
t67 = t94 + t95;
t269 = t188 * t67;
t223 = t191 * t159 - t280 * t195;
t100 = -qJDD(2) * pkin(2) - t197 * pkin(8) + t223;
t202 = -t153 * pkin(3) - t281 * qJ(4) + t100;
t73 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t244 + t202;
t65 = t73 + t274;
t268 = t189 * t65;
t267 = t192 * t22;
t266 = t192 * t67;
t265 = t193 * t65;
t261 = t136 * t188;
t260 = t136 * t192;
t108 = -t259 + t178;
t256 = t189 * t108;
t154 = -0.2e1 * t174 + t239;
t255 = t190 * t154;
t254 = t190 * t160;
t251 = t193 * t108;
t165 = -t196 - t257;
t124 = t194 * t165 - t254;
t247 = pkin(2) * t154 + pkin(8) * t124;
t245 = t181 + t182;
t156 = t245 * qJDD(2);
t157 = t245 * t197;
t246 = pkin(2) * t157 + pkin(8) * t156;
t242 = qJD(6) + t136;
t235 = t189 * t94;
t234 = t193 * t94;
t231 = pkin(5) * t189 + qJ(4);
t23 = -t236 * pkin(5) - t178 * pkin(10) - t140 * t111 + t35;
t201 = t244 * t273 - t202;
t199 = -pkin(3) * t174 + t201;
t226 = t142 * t238;
t36 = t199 - t88 * pkin(10) + (-t97 - t226) * pkin(5) - t274;
t14 = t188 * t23 - t192 * t36;
t15 = t188 * t36 + t192 * t23;
t6 = t188 * t14 + t192 * t15;
t47 = t190 * t81 + t194 * t82;
t227 = pkin(5) * t193 + t272;
t224 = -pkin(5) * t22 + pkin(10) * t6;
t2 = t189 * t6 - t193 * t22;
t3 = t189 * t22 + t193 * t6;
t1 = t190 * t2 + t194 * t3;
t222 = t189 * t226;
t221 = t189 * t225;
t220 = t193 * t226;
t219 = t193 * t225;
t5 = -t192 * t14 + t188 * t15;
t16 = t189 * t35 - t193 * t34;
t17 = t189 * t34 + t193 * t35;
t7 = t190 * t16 + t194 * t17;
t216 = t188 * t178 - t192 * t98;
t215 = t194 * t151 + t255;
t212 = pkin(2) - t217;
t85 = t142 * qJD(3) + t228;
t92 = -t116 - t134;
t41 = -t188 * t92 - t266;
t63 = t242 * t117 + t216;
t209 = pkin(5) * t63 + pkin(10) * t41 + t270;
t83 = -t134 - t115;
t39 = t192 * t83 - t283;
t60 = -t242 * t119 - t230;
t208 = pkin(5) * t60 + pkin(10) * t39 - t267;
t104 = t136 * t117;
t75 = -t117 * qJD(6) - t216;
t62 = t104 + t75;
t31 = t188 * t62 - t192 * t58;
t76 = t115 + t116;
t204 = pkin(5) * t76 + pkin(10) * t31 + t6;
t70 = t190 * t210 + t206;
t158 = (-t181 + t182) * t197;
t131 = -t139 + t236;
t130 = t138 - t236;
t127 = -t139 - t236;
t123 = t254 + t194 * (t196 - t258);
t122 = t281 * t190;
t121 = (t153 - t174) * t194;
t113 = t139 - t138;
t112 = (t156 * t191 + t157 * t195) * t185;
t107 = -t236 - t138;
t103 = -t116 + t134;
t102 = t115 - t134;
t99 = -t138 - t139;
t93 = t116 - t115;
t91 = -t189 * t127 + t251;
t90 = t193 * t127 + t256;
t89 = -t225 + t98;
t84 = (0.2e1 * qJD(5) - qJD(3)) * t142 + t228;
t79 = t186 * (t194 * t160 + t190 * t165) + (t191 * t124 + t195 * t154) * t185;
t78 = t193 * t107 - t289;
t77 = t189 * t107 + t288;
t74 = -t119 * qJD(6) - t230;
t72 = (-t117 * t192 + t119 * t188) * t136;
t71 = (-t117 * t188 - t119 * t192) * t136;
t61 = -t104 + t75;
t54 = -t119 * t261 + t192 * t75;
t53 = t119 * t260 + t188 * t75;
t52 = t117 * t260 - t188 * t74;
t51 = -t117 * t261 - t192 * t74;
t50 = t190 * t90 + t194 * t91;
t49 = t189 * t89 - t193 * t85;
t48 = -t189 * t85 - t193 * t89;
t46 = t192 * t102 - t269;
t45 = -t188 * t103 + t282;
t44 = t188 * t102 + t266;
t43 = t192 * t103 + t283;
t42 = t190 * t77 + t194 * t78;
t40 = t192 * t92 - t269;
t38 = t188 * t83 + t282;
t37 = t190 * t70 + t194 * t69;
t32 = -t188 * t61 + t192 * t60;
t30 = t188 * t60 + t192 * t61;
t29 = -t188 * t58 - t192 * t62;
t28 = t190 * t48 + t194 * t49;
t27 = -t189 * t63 + t193 * t41;
t26 = t189 * t41 + t193 * t63;
t25 = -t189 * t60 + t193 * t39;
t24 = t189 * t39 + t193 * t60;
t21 = -t189 * t76 + t193 * t31;
t20 = t189 * t31 + t193 * t76;
t19 = -pkin(10) * t40 + t267;
t18 = -pkin(10) * t38 + t270;
t12 = t190 * t26 + t194 * t27;
t11 = t190 * t24 + t194 * t25;
t10 = t190 * t20 + t194 * t21;
t9 = -pkin(5) * t40 + t15;
t8 = -pkin(5) * t38 + t14;
t4 = -pkin(10) * t29 - t5;
t13 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t248, 0, 0, 0, 0, 0, 0, (qJDD(2) * t195 - t191 * t197) * t185, (-qJDD(2) * t191 - t195 * t197) * t185, 0, t186 ^ 2 * t248 + (t191 * t106 - t195 * t223 - t203) * t185, 0, 0, 0, 0, 0, 0, t79, -t287, t112, t186 * (t190 * t82 - t194 * t81) + (-t195 * t100 + t191 * t47) * t185, 0, 0, 0, 0, 0, 0, t79, t112, t287, t186 * (t190 * t69 - t194 * t70) + (t191 * t37 - t195 * t73) * t185, 0, 0, 0, 0, 0, 0, t186 * (t190 * t78 - t194 * t77) + (t191 * t42 + t195 * t84) * t185, t186 * (t190 * t91 - t194 * t90) + (t191 * t50 + t195 * t88) * t185, t186 * (t190 * t49 - t194 * t48) + (t191 * t28 + t195 * t99) * t185, t186 * (-t194 * t16 + t190 * t17) + (t191 * t7 - t195 * t65) * t185, 0, 0, 0, 0, 0, 0, t186 * (t190 * t25 - t194 * t24) + (t191 * t11 + t195 * t38) * t185, t186 * (t190 * t27 - t194 * t26) + (t191 * t12 + t195 * t40) * t185, t186 * (t190 * t21 - t194 * t20) + (t191 * t10 + t195 * t29) * t185, t186 * (t190 * t3 - t194 * t2) + (t191 * t1 + t195 * t5) * t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t223, -t106, 0, 0, t122, t215, t123, t121, t275, 0, -t194 * t100 + t247, -pkin(2) * t151 + t190 * t100 - t286, t47 + t246, -pkin(2) * t100 + pkin(8) * t47, t122, t123, -t215, 0, -t275, t121, qJ(4) * t255 + t194 * ((t154 - t174) * pkin(3) + t201) + t247, t190 * (qJ(4) * t157 + t70) + t194 * ((t157 - t196) * pkin(3) + t218) + t246, t151 * t212 + t190 * t199 + t286, pkin(8) * t37 - t212 * t73, t190 * (t193 * t98 + t222) + t194 * (-t189 * t98 + t220), t190 * (-t189 * t88 - t193 * t84) + t194 * (t189 * t84 - t193 * t88), t190 * (-t189 * t131 + t288) + t194 * (-t193 * t131 - t289), t190 * (-t189 * t97 - t219) + t194 * (-t193 * t97 + t221), t190 * (t193 * t130 + t256) + t194 * (-t189 * t130 + t251), t190 * (t219 - t222) + t194 * (-t221 - t220), t190 * (-pkin(9) * t77 - t268) + t194 * (-pkin(9) * t78 - t265) + pkin(8) * t42 + t207 * t84, t190 * (-pkin(9) * t90 - t265) + t194 * (-pkin(9) * t91 + t268) + pkin(8) * t50 + t207 * t88, t190 * (-pkin(9) * t48 - t16) + t194 * (-pkin(9) * t49 - t17) + pkin(8) * t28 + t207 * t99, -t207 * t65 + t279 * t7, t190 * (t193 * t54 + t235) + t194 * (-t189 * t54 + t234), t190 * (t189 * t93 + t193 * t32) + t194 * (-t189 * t32 + t193 * t93), t190 * (t189 * t62 + t193 * t45) + t194 * (-t189 * t45 + t193 * t62), t190 * (t193 * t52 - t235) + t194 * (-t189 * t52 - t234), t190 * (-t189 * t58 + t193 * t46) + t194 * (-t189 * t46 - t193 * t58), t190 * (t189 * t95 + t193 * t72) + t194 * (-t189 * t72 + t193 * t95), t190 * (-pkin(9) * t24 + t193 * t18 - t189 * t8) + t194 * (-pkin(9) * t25 - t189 * t18 - t193 * t8) + pkin(8) * t11 + t207 * t38, t190 * (-pkin(9) * t26 - t189 * t9 + t193 * t19) + t194 * (-pkin(9) * t27 - t189 * t19 - t193 * t9) + pkin(8) * t12 + t207 * t40, t190 * (-pkin(9) * t20 + t193 * t4) + t194 * (-pkin(9) * t21 - t189 * t4) + pkin(8) * t10 + (t190 * t231 + t194 * t227 + pkin(2)) * t29, (t190 * (-pkin(10) * t193 + t231) + t194 * (pkin(10) * t189 + t227) + pkin(2)) * t5 + t279 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t168, -t158, t240, t168, t239, qJDD(3), -t81, -t82, 0, 0, -t168, t240, t158, qJDD(3), -t239, t168, pkin(3) * t160 + qJ(4) * t165 - t70, (-t190 * pkin(3) + qJ(4) * t194) * qJDD(2), qJ(4) * t161 + (t163 - t196) * pkin(3) + t218, -pkin(3) * t70 + qJ(4) * t69, -t259, -t113, -t89, t259, t85, t178, qJ(4) * t78 - t272 * t77 + t34, qJ(4) * t91 - t272 * t90 + t35, qJ(4) * t49 - t272 * t48, qJ(4) * t17 - t272 * t16, -t53, -t30, -t43, t51, -t44, -t71, qJ(4) * t25 - t272 * t24 - t208, qJ(4) * t27 - t272 * t26 - t209, qJ(4) * t21 - t272 * t20 - t204, qJ(4) * t3 - t272 * t2 - t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t160, t240, -t163, t70, 0, 0, 0, 0, 0, 0, t77, t90, t48, t16, 0, 0, 0, 0, 0, 0, t24, t26, t20, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t259, t113, t89, -t259, -t85, -t178, -t34, -t35, 0, 0, t53, t30, t43, -t51, t44, t71, t208, t209, t204, t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t93, t62, -t94, -t58, t95, -t14, -t15, 0, 0;];
tauJ_reg  = t13;
