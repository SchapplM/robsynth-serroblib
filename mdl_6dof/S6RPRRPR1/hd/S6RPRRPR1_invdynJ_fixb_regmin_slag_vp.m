% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPRRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRRPR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRPR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:58:48
% EndTime: 2019-03-09 04:58:56
% DurationCPUTime: 2.88s
% Computational Cost: add. (6007->354), mult. (13421->470), div. (0->0), fcn. (9889->18), ass. (0->217)
t182 = qJ(3) + qJ(4);
t173 = sin(t182);
t174 = cos(t182);
t179 = qJ(1) + pkin(10);
t168 = sin(t179);
t169 = cos(t179);
t230 = g(1) * t169 + g(2) * t168;
t311 = -g(3) * t174 + t230 * t173;
t184 = sin(pkin(10));
t159 = pkin(1) * t184 + pkin(7);
t287 = pkin(8) + t159;
t171 = pkin(11) + t182;
t157 = cos(t171);
t176 = qJDD(3) + qJDD(4);
t187 = sin(qJ(4));
t191 = cos(qJ(4));
t192 = cos(qJ(3));
t263 = qJD(1) * t192;
t188 = sin(qJ(3));
t264 = qJD(1) * t188;
t119 = -t187 * t263 - t191 * t264;
t242 = t287 * qJD(1);
t92 = t188 * qJD(2) + t192 * t242;
t89 = t191 * t92;
t283 = qJD(3) * pkin(3);
t91 = t192 * qJD(2) - t242 * t188;
t90 = t91 + t283;
t223 = -t187 * t90 - t89;
t170 = t192 * qJDD(2);
t137 = t159 * qJDD(1);
t236 = pkin(8) * qJDD(1) + t137;
t62 = qJDD(3) * pkin(3) - qJD(3) * t92 - t236 * t188 + t170;
t65 = qJD(3) * t91 + t188 * qJDD(2) + t236 * t192;
t204 = qJD(4) * t223 - t187 * t65 + t191 * t62;
t178 = qJD(3) + qJD(4);
t257 = qJD(1) * qJD(3);
t249 = t192 * t257;
t251 = t191 * t263;
t252 = t187 * t264;
t255 = t192 * qJDD(1);
t256 = t188 * qJDD(1);
t67 = qJD(4) * t251 - t178 * t252 + t187 * t255 + (t249 + t256) * t191;
t17 = t176 * pkin(4) - t67 * qJ(5) + t119 * qJD(5) + t204;
t183 = sin(pkin(11));
t117 = -t251 + t252;
t262 = qJD(4) * t187;
t300 = t191 * (qJD(4) * t90 + t65) + t187 * t62 - t92 * t262;
t224 = t187 * t256 - t191 * t255;
t126 = t187 * t192 + t188 * t191;
t85 = t178 * t126;
t68 = qJD(1) * t85 + t224;
t19 = -t68 * qJ(5) - t117 * qJD(5) + t300;
t277 = cos(pkin(11));
t5 = t277 * t17 - t183 * t19;
t3 = -t176 * pkin(5) - t5;
t310 = g(3) * t157 + t3;
t233 = -t277 * t117 + t119 * t183;
t302 = qJD(6) - t233;
t186 = sin(qJ(6));
t190 = cos(qJ(6));
t210 = -t183 * t117 - t277 * t119;
t70 = -t190 * t178 + t186 * t210;
t309 = t302 * t70;
t237 = t190 * t302;
t41 = -t183 * t67 - t277 * t68;
t40 = qJDD(6) - t41;
t307 = -t186 * t40 - t237 * t302;
t158 = pkin(4) * t183 + pkin(9);
t298 = pkin(4) * t119;
t46 = pkin(5) * t210 - pkin(9) * t233 - t298;
t306 = (qJD(6) * t158 + t46) * t302;
t165 = pkin(3) * t191 + pkin(4);
t243 = t277 * t187;
t112 = pkin(3) * t243 + t183 * t165;
t105 = pkin(9) + t112;
t166 = pkin(3) * t264;
t305 = (qJD(6) * t105 + t166 + t46) * t302;
t123 = t287 * t188;
t124 = t287 * t192;
t267 = -t187 * t123 + t191 * t124;
t113 = t119 * qJ(5);
t87 = t187 * t92;
t246 = t191 * t90 - t87;
t54 = t113 + t246;
t185 = cos(pkin(10));
t161 = -pkin(1) * t185 - pkin(2);
t175 = t192 * pkin(3);
t303 = t161 - t175;
t125 = t187 * t188 - t191 * t192;
t84 = t178 * t125;
t60 = -t183 * t85 - t277 * t84;
t82 = -t183 * t125 + t277 * t126;
t225 = t302 * t60 + t40 * t82;
t260 = qJD(6) * t186;
t254 = t82 * t260;
t301 = -t190 * t225 + t254 * t302;
t244 = qJD(3) * t287;
t114 = t188 * t244;
t115 = t192 * t244;
t203 = -qJD(4) * t267 + t187 * t114 - t191 * t115;
t198 = t84 * qJ(5) - t126 * qJD(5) + t203;
t261 = qJD(4) * t191;
t211 = -t191 * t114 - t187 * t115 - t123 * t261 - t124 * t262;
t35 = -qJ(5) * t85 - qJD(5) * t125 + t211;
t15 = t183 * t198 + t277 * t35;
t120 = t303 * qJD(1);
t80 = t117 * pkin(4) + qJD(5) + t120;
t37 = -pkin(5) * t233 - pkin(9) * t210 + t80;
t6 = t183 * t17 + t277 * t19;
t248 = pkin(9) * t176 + qJD(6) * t37 + t6;
t276 = qJ(5) * t117;
t55 = -t223 - t276;
t282 = t183 * t55;
t50 = pkin(4) * t178 + t54;
t28 = t277 * t50 - t282;
t26 = -t178 * pkin(5) - t28;
t232 = -t191 * t123 - t124 * t187;
t213 = -qJ(5) * t126 + t232;
t69 = -qJ(5) * t125 + t267;
t44 = t183 * t213 + t277 * t69;
t212 = pkin(4) * t125 + t303;
t81 = t277 * t125 + t126 * t183;
t47 = pkin(5) * t81 - pkin(9) * t82 + t212;
t299 = -(qJD(6) * t47 + t15) * t302 - t248 * t81 + t26 * t60 + t3 * t82 - t44 * t40;
t293 = t26 * t233;
t292 = t26 * t82;
t291 = t47 * t40;
t290 = t70 * t210;
t72 = t178 * t186 + t190 * t210;
t289 = t72 * t210;
t288 = t302 * t210;
t259 = qJD(6) * t190;
t42 = -t183 * t68 + t277 * t67;
t24 = t186 * t176 + t178 * t259 + t190 * t42 - t210 * t260;
t59 = -t183 * t84 + t277 * t85;
t286 = t24 * t81 + t72 * t59;
t51 = t277 * t55;
t29 = t183 * t50 + t51;
t285 = t191 * t91 - t87;
t284 = pkin(3) * qJD(4);
t280 = t24 * t186;
t245 = -t187 * t91 - t89;
t215 = t245 + t276;
t56 = t113 + t285;
t279 = -t183 * t56 + t277 * t215 + (t183 * t191 + t243) * t284;
t270 = t183 * t187;
t278 = -t183 * t215 - t277 * t56 + (t277 * t191 - t270) * t284;
t275 = t119 * t117;
t274 = t168 * t186;
t273 = t168 * t190;
t272 = t169 * t186;
t271 = t169 * t190;
t268 = qJDD(2) - g(3);
t266 = pkin(4) * t174 + t175;
t180 = t188 ^ 2;
t265 = -t192 ^ 2 + t180;
t140 = qJD(1) * t161;
t167 = t188 * t283;
t250 = pkin(4) * t85 + t167;
t97 = qJD(3) * t166 + qJDD(1) * t303;
t200 = t68 * pkin(4) + qJDD(5) + t97;
t11 = -t41 * pkin(5) - t42 * pkin(9) + t200;
t27 = pkin(9) * t178 + t29;
t240 = qJD(6) * t27 - t11;
t238 = -t190 * t176 + t186 * t42;
t229 = g(1) * t168 - g(2) * t169;
t189 = sin(qJ(1));
t193 = cos(qJ(1));
t228 = g(1) * t189 - g(2) * t193;
t25 = qJD(6) * t72 + t238;
t227 = -t81 * t25 - t59 * t70;
t226 = t210 * t29 + t233 * t28;
t13 = t186 * t37 + t190 * t27;
t222 = t13 * t210 + t186 * t310 + t26 * t259;
t220 = t126 * t176 - t178 * t84;
t12 = -t186 * t27 + t190 * t37;
t156 = sin(t171);
t219 = -t12 * t210 + t26 * t260 + (g(1) * t271 + g(2) * t273) * t156;
t217 = t190 * t40 + (t186 * t233 - t260) * t302;
t216 = g(3) * t156 - t248;
t214 = t230 * t156;
t31 = t277 * t54 - t282;
t208 = -t158 * t40 + t302 * t31 - t293;
t111 = -pkin(3) * t270 + t277 * t165;
t207 = -qJD(1) * t140 - t137 + t230;
t206 = -t105 * t40 - t278 * t302 - t293;
t205 = 0.2e1 * t140 * qJD(3) - qJDD(3) * t159;
t194 = qJD(3) ^ 2;
t202 = -0.2e1 * qJDD(1) * t161 - t159 * t194 + t229;
t201 = -t259 * t302 * t82 - t186 * t225;
t199 = g(3) * t173 + t120 * t117 + t174 * t230 - t300;
t197 = t120 * t119 + t204 + t311;
t195 = qJD(1) ^ 2;
t177 = -qJ(5) - pkin(8) - pkin(7);
t160 = -t277 * pkin(4) - pkin(5);
t135 = qJDD(3) * t192 - t188 * t194;
t134 = qJDD(3) * t188 + t192 * t194;
t127 = pkin(2) + t266;
t104 = -pkin(5) - t111;
t96 = t157 * t271 + t274;
t95 = -t157 * t272 + t273;
t94 = -t157 * t273 + t272;
t93 = t157 * t274 + t271;
t73 = -t117 ^ 2 + t119 ^ 2;
t66 = -t125 * t176 - t178 * t85;
t58 = -t224 + (-qJD(1) * t126 - t119) * t178;
t57 = t117 * t178 + t67;
t43 = t183 * t69 - t277 * t213;
t30 = t183 * t54 + t51;
t21 = pkin(5) * t59 - pkin(9) * t60 + t250;
t14 = t183 * t35 - t277 * t198;
t10 = t190 * t11;
t9 = t237 * t72 + t280;
t8 = -t289 - t307;
t7 = t217 + t290;
t1 = (t24 - t309) * t190 + (-t302 * t72 - t25) * t186;
t2 = [qJDD(1), t228, g(1) * t193 + g(2) * t189 (t228 + (t184 ^ 2 + t185 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), qJDD(1) * t180 + 0.2e1 * t188 * t249, 0.2e1 * t188 * t255 - 0.2e1 * t265 * t257, t134, t135, 0, t188 * t205 + t192 * t202, -t188 * t202 + t192 * t205, t119 * t84 + t126 * t67, t117 * t84 + t119 * t85 - t125 * t67 - t126 * t68, t220, t66, 0, t117 * t167 + t120 * t85 + t97 * t125 + t174 * t229 + t176 * t232 + t178 * t203 + t303 * t68, -t119 * t167 - t120 * t84 + t97 * t126 - t229 * t173 - t267 * t176 - t211 * t178 + t303 * t67, t14 * t210 + t15 * t233 - t28 * t60 - t29 * t59 + t41 * t44 + t42 * t43 - t5 * t82 - t6 * t81 - t230, t6 * t44 + t29 * t15 - t5 * t43 - t28 * t14 + t200 * t212 + t80 * t250 - g(1) * (-pkin(1) * t189 - t127 * t168 - t169 * t177) - g(2) * (pkin(1) * t193 + t127 * t169 - t168 * t177) -t72 * t254 + (t24 * t82 + t60 * t72) * t190 (-t186 * t72 - t190 * t70) * t60 + (-t280 - t190 * t25 + (t186 * t70 - t190 * t72) * qJD(6)) * t82, t286 - t301, t201 + t227, t302 * t59 + t40 * t81, -g(1) * t94 - g(2) * t96 + t10 * t81 + t12 * t59 + t14 * t70 + t43 * t25 + (t21 * t302 + t291 + (-t27 * t81 - t302 * t44 + t292) * qJD(6)) * t190 + t299 * t186, -g(1) * t93 - g(2) * t95 - t13 * t59 + t14 * t72 + t43 * t24 + (-(-qJD(6) * t44 + t21) * t302 - t291 + t240 * t81 - qJD(6) * t292) * t186 + t299 * t190; 0, 0, 0, t268, 0, 0, 0, 0, 0, t135, -t134, 0, 0, 0, 0, 0, t66, -t220, t210 * t59 + t233 * t60 + t41 * t82 + t42 * t81, -t28 * t59 + t29 * t60 - t5 * t81 + t6 * t82 - g(3), 0, 0, 0, 0, 0, t201 - t227, t286 + t301; 0, 0, 0, 0, -t188 * t195 * t192, t265 * t195, t256, t255, qJDD(3), -g(3) * t192 + t188 * t207 + t170, -t268 * t188 + t207 * t192, -t275, t73, t57, t58, t176, -t245 * t178 + (-t117 * t264 + t176 * t191 - t178 * t262) * pkin(3) + t197, t285 * t178 + (t119 * t264 - t176 * t187 - t178 * t261) * pkin(3) + t199, -t111 * t42 + t112 * t41 + t210 * t279 + t233 * t278 + t226, t6 * t112 + t5 * t111 - t80 * (t166 - t298) - g(3) * t266 + t278 * t29 - t279 * t28 - t230 * (-pkin(3) * t188 - pkin(4) * t173) t9, t1, t8, t7, -t288, t104 * t25 + t279 * t70 + (-t310 - t305) * t190 + t206 * t186 + t219, t104 * t24 + t279 * t72 + t206 * t190 + (-t214 + t305) * t186 + t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t275, t73, t57, t58, t176, -t178 * t223 + t197, t178 * t246 + t199, -t30 * t210 - t31 * t233 + (t183 * t41 - t277 * t42) * pkin(4) + t226, t28 * t30 - t29 * t31 + (t119 * t80 + t183 * t6 + t277 * t5 + t311) * pkin(4), t9, t1, t8, t7, -t288, t160 * t25 - t30 * t70 + t208 * t186 + (-t310 - t306) * t190 + t219, t160 * t24 - t30 * t72 + t208 * t190 + (-t214 + t306) * t186 + t222; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t210 ^ 2 - t233 ^ 2, t210 * t28 - t233 * t29 + t200 - t229, 0, 0, 0, 0, 0, t217 - t290, -t289 + t307; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 * t70, -t70 ^ 2 + t72 ^ 2, t24 + t309, -t238 + (-qJD(6) + t302) * t72, t40, -g(1) * t95 + g(2) * t93 + t13 * t302 + t186 * t216 - t27 * t259 - t26 * t72 + t10, g(1) * t96 - g(2) * t94 + t12 * t302 + t186 * t240 + t190 * t216 + t26 * t70;];
tau_reg  = t2;
