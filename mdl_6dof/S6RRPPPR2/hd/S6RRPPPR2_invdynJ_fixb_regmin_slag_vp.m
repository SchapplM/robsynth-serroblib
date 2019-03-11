% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RRPPPR2
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% 
% Output:
% tau_reg [6x27]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RRPPPR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RRPPPR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPPPR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:12:10
% EndTime: 2019-03-09 08:12:21
% DurationCPUTime: 5.15s
% Computational Cost: add. (4726->446), mult. (10923->562), div. (0->0), fcn. (8073->14), ass. (0->228)
t189 = sin(pkin(9));
t191 = cos(pkin(9));
t194 = sin(qJ(2));
t197 = cos(qJ(2));
t145 = t189 * t197 + t191 * t194;
t134 = t145 * qJD(1);
t124 = qJD(6) + t134;
t188 = sin(pkin(10));
t190 = cos(pkin(10));
t193 = sin(qJ(6));
t196 = cos(qJ(6));
t144 = t188 * t196 + t190 * t193;
t277 = t124 * t144;
t306 = -t188 * t193 + t190 * t196;
t254 = qJD(1) * qJD(2);
t245 = t197 * t254;
t246 = t194 * t254;
t99 = qJDD(1) * t145 - t189 * t246 + t191 * t245;
t94 = qJDD(6) + t99;
t312 = -t124 * t277 + t306 * t94;
t265 = t191 * t197;
t143 = t189 * t194 - t265;
t86 = t306 * t143;
t256 = qJD(6) * t196;
t257 = qJD(6) * t193;
t276 = t306 * t134 - t188 * t257 + t190 * t256;
t238 = -t124 * t276 - t144 * t94;
t259 = qJD(1) * t194;
t131 = -qJD(1) * t265 + t189 * t259;
t109 = qJD(2) * t188 - t190 * t131;
t111 = qJD(2) * t190 + t131 * t188;
t223 = t109 * t193 - t111 * t196;
t311 = t124 * t223;
t179 = t197 * pkin(2);
t172 = t179 + pkin(1);
t182 = qJ(2) + pkin(9);
t178 = cos(t182);
t195 = sin(qJ(1));
t198 = cos(qJ(1));
t307 = g(1) * t195 - g(2) * t198;
t310 = t178 * t307;
t176 = sin(t182);
t236 = g(1) * t198 + g(2) * t195;
t309 = t236 * t176;
t308 = -t178 * pkin(3) - t176 * qJ(4);
t299 = t134 ^ 2;
t305 = -t131 ^ 2 - t299;
t285 = qJ(3) + pkin(7);
t153 = t285 * t197;
t149 = qJD(1) * t153;
t139 = t189 * t149;
t152 = t285 * t194;
t148 = qJD(1) * t152;
t104 = -t191 * t148 - t139;
t304 = qJD(4) - t104;
t303 = -qJD(6) + t124;
t133 = t145 * qJD(2);
t252 = t197 * qJDD(1);
t253 = t194 * qJDD(1);
t229 = -t189 * t253 + t191 * t252;
t98 = qJD(1) * t133 - t229;
t302 = -t98 * pkin(4) + qJDD(5);
t169 = g(3) * t178;
t213 = t169 - t309;
t282 = t188 * t99;
t301 = -t190 * t299 - t282;
t78 = qJDD(2) * t188 - t190 * t98;
t79 = qJDD(2) * t190 + t188 * t98;
t19 = -qJD(6) * t223 + t193 * t79 + t196 * t78;
t298 = t98 * pkin(3);
t296 = pkin(2) * t194;
t295 = pkin(4) * t131;
t294 = pkin(8) * t190;
t290 = g(3) * t197;
t289 = t134 * pkin(4);
t288 = pkin(3) + qJ(5);
t287 = pkin(4) + t285;
t171 = -pkin(2) * t191 - pkin(3);
t164 = -qJ(5) + t171;
t286 = -pkin(8) + t164;
t212 = pkin(2) * t246 - qJDD(1) * t172 + qJDD(3);
t279 = t99 * qJ(4);
t203 = -t134 * qJD(4) + t212 - t279;
t13 = t131 * qJD(5) + t288 * t98 + t203;
t242 = qJD(2) * t285;
t125 = t197 * qJD(3) - t194 * t242;
t102 = qJD(1) * t125 + qJDD(1) * t153;
t126 = -t194 * qJD(3) - t197 * t242;
t93 = qJDD(2) * pkin(2) + qJD(1) * t126 - qJDD(1) * t152;
t42 = -t189 * t102 + t191 * t93;
t231 = qJDD(4) - t42;
t26 = t99 * pkin(4) - qJD(2) * qJD(5) - t288 * qJDD(2) + t231;
t7 = t190 * t13 + t188 * t26;
t258 = qJD(2) * t194;
t136 = qJD(2) * t265 - t189 * t258;
t174 = pkin(2) * t258;
t217 = -t136 * qJ(4) - t145 * qJD(4) + t174;
t35 = t143 * qJD(5) + t288 * t133 + t217;
t81 = t125 * t189 - t191 * t126;
t55 = pkin(4) * t136 + t81;
t16 = t188 * t55 + t190 * t35;
t151 = -qJD(1) * t172 + qJD(3);
t210 = -t134 * qJ(4) + t151;
t48 = t288 * t131 + t210;
t142 = qJD(2) * pkin(2) - t148;
t96 = t191 * t142 - t139;
t230 = qJD(4) - t96;
t52 = -t288 * qJD(2) + t230 + t289;
t23 = t188 * t52 + t190 * t48;
t241 = pkin(2) * t259 + t131 * qJ(4);
t53 = t288 * t134 + t241;
t266 = t191 * t149;
t103 = -t148 * t189 + t266;
t66 = t103 - t295;
t30 = t188 * t66 + t190 * t53;
t239 = -t145 * qJ(4) - t172;
t65 = t288 * t143 + t239;
t107 = t191 * t152 + t153 * t189;
t83 = pkin(4) * t145 + t107;
t33 = t188 * t83 + t190 * t65;
t43 = t191 * t102 + t189 * t93;
t57 = t196 * t109 + t111 * t193;
t284 = t124 * t57;
t283 = t131 * t57;
t281 = t190 * t99;
t280 = t223 * t131;
t275 = qJ(5) * t178;
t274 = qJDD(2) * pkin(3);
t181 = pkin(10) + qJ(6);
t175 = sin(t181);
t273 = t175 * t195;
t272 = t176 * t198;
t177 = cos(t181);
t271 = t177 * t195;
t270 = t177 * t198;
t269 = t178 * t198;
t264 = t198 * t285;
t248 = -pkin(5) * t190 - pkin(4);
t263 = -t134 * t248 + t304;
t262 = t289 + t304;
t97 = t189 * t142 + t266;
t89 = -qJD(2) * qJ(4) - t97;
t60 = qJD(5) - t89 - t295;
t261 = qJD(5) - t60;
t186 = t194 ^ 2;
t260 = -t197 ^ 2 + t186;
t251 = qJDD(2) * qJ(4) + t43;
t250 = t179 - t308;
t6 = -t13 * t188 + t190 * t26;
t2 = pkin(5) * t99 - pkin(8) * t79 + t6;
t5 = -pkin(8) * t78 + t7;
t247 = -t193 * t5 + t196 * t2;
t167 = pkin(2) * t189 + qJ(4);
t22 = -t188 * t48 + t190 * t52;
t184 = qJD(2) * qJD(4);
t39 = -t184 - t251;
t159 = t198 * t172;
t240 = g(2) * (pkin(3) * t269 + qJ(4) * t272 + t159);
t237 = -pkin(3) * t176 - t296;
t234 = t7 * t188 + t6 * t190;
t233 = -t6 * t188 + t7 * t190;
t232 = t193 * t2 + t196 * t5;
t10 = pkin(5) * t134 - pkin(8) * t111 + t22;
t14 = -pkin(8) * t109 + t23;
t3 = t10 * t196 - t14 * t193;
t4 = t10 * t193 + t14 * t196;
t228 = -t188 * t23 - t190 * t22;
t227 = -t188 * t22 + t190 * t23;
t77 = t190 * t83;
t21 = t145 * pkin(5) + t77 + (-pkin(8) * t143 - t65) * t188;
t28 = t143 * t294 + t33;
t226 = -t193 * t28 + t196 * t21;
t225 = t193 * t21 + t196 * t28;
t224 = -t299 * t188 + t281;
t82 = t191 * t125 + t189 * t126;
t108 = -t189 * t152 + t191 * t153;
t220 = -t172 + t308;
t219 = t307 * t176;
t218 = -0.2e1 * pkin(1) * t254 - pkin(7) * qJDD(2);
t129 = t286 * t188;
t64 = t190 * t66;
t216 = qJD(5) * t190 + qJD(6) * t129 - t131 * pkin(5) + t64 + (-pkin(8) * t134 - t53) * t188;
t130 = t286 * t190;
t215 = qJD(5) * t188 - qJD(6) * t130 + t134 * t294 + t30;
t18 = -t109 * t256 - t111 * t257 - t193 * t78 + t196 * t79;
t87 = t144 * t143;
t27 = -t39 + t302;
t211 = -g(3) * t176 - t178 * t236;
t199 = qJD(2) ^ 2;
t209 = 0.2e1 * qJDD(1) * pkin(1) - pkin(7) * t199 + t307;
t200 = qJD(1) ^ 2;
t208 = pkin(1) * t200 - pkin(7) * qJDD(1) + t236;
t207 = t133 * t60 + t143 * t27 + t236;
t206 = t27 + t211;
t205 = t212 - t307;
t204 = t211 + t251;
t202 = t107 * t99 - t108 * t98 - t131 * t82 + t134 * t81 - t236;
t73 = t131 * pkin(3) + t210;
t201 = t73 * t134 + t213 + t231;
t156 = qJ(4) * t269;
t154 = t195 * t178 * qJ(4);
t150 = pkin(5) * t188 + t167;
t123 = qJD(2) * t131;
t118 = -t176 * t273 + t270;
t117 = t175 * t198 + t176 * t271;
t116 = t175 * t272 + t271;
t115 = t176 * t270 - t273;
t95 = pkin(3) * t143 + t239;
t85 = -qJD(2) * pkin(3) + t230;
t84 = -pkin(4) * t143 + t108;
t80 = pkin(3) * t134 + t241;
t61 = pkin(3) * t133 + t217;
t56 = -pkin(4) * t133 + t82;
t54 = t143 * t248 + t108;
t51 = t190 * t55;
t41 = t133 * t248 + t82;
t40 = t231 - t274;
t38 = pkin(5) * t109 + t60;
t37 = qJD(6) * t87 - t306 * t133;
t36 = qJD(6) * t86 + t133 * t144;
t32 = -t188 * t65 + t77;
t31 = t203 + t298;
t29 = -t188 * t53 + t64;
t15 = -t188 * t35 + t51;
t11 = pkin(5) * t78 + t27;
t9 = t133 * t294 + t16;
t8 = t136 * pkin(5) + t51 + (-pkin(8) * t133 - t35) * t188;
t1 = [qJDD(1), t307, t236, qJDD(1) * t186 + 0.2e1 * t194 * t245, 0.2e1 * t194 * t252 - 0.2e1 * t254 * t260, qJDD(2) * t194 + t197 * t199, qJDD(2) * t197 - t194 * t199, 0, t194 * t218 + t197 * t209, -t194 * t209 + t197 * t218, -t133 * t97 - t136 * t96 - t143 * t43 - t145 * t42 + t202, t43 * t108 + t97 * t82 - t42 * t107 - t96 * t81 - t212 * t172 + t151 * t174 - g(1) * (-t195 * t172 + t264) - g(2) * (t195 * t285 + t159) t133 * t89 + t136 * t85 + t143 * t39 + t145 * t40 + t202, t81 * qJD(2) + t107 * qJDD(2) - t61 * t131 - t73 * t133 - t31 * t143 - t95 * t98 - t310, t82 * qJD(2) + t108 * qJDD(2) - t61 * t134 - t73 * t136 - t31 * t145 - t95 * t99 + t219, t31 * t95 + t73 * t61 - t39 * t108 - t89 * t82 + t40 * t107 + t85 * t81 - g(1) * t264 - t240 + (-g(1) * t220 - g(2) * t285) * t195, t56 * t109 + t15 * t134 + t22 * t136 + t6 * t145 + t188 * t219 - t190 * t207 + t32 * t99 + t84 * t78, t56 * t111 - t16 * t134 - t23 * t136 - t7 * t145 + t188 * t207 + t190 * t219 - t33 * t99 + t84 * t79, -t16 * t109 - t15 * t111 + t133 * t227 + t143 * t233 - t32 * t79 - t33 * t78 + t310, t7 * t33 + t23 * t16 + t6 * t32 + t22 * t15 + t27 * t84 + t60 * t56 - t240 + (-g(1) * t287 - g(2) * t275) * t198 + (-g(1) * (t220 - t275) - g(2) * t287) * t195, t18 * t87 - t223 * t36, t18 * t86 - t19 * t87 + t223 * t37 - t36 * t57, t124 * t36 - t136 * t223 + t145 * t18 + t87 * t94, -t124 * t37 - t136 * t57 - t145 * t19 + t86 * t94, t124 * t136 + t145 * t94 (-t193 * t9 + t196 * t8) * t124 + t226 * t94 + t247 * t145 + t3 * t136 + t41 * t57 + t54 * t19 - t11 * t86 + t38 * t37 - g(1) * t118 - g(2) * t116 + (-t124 * t225 - t145 * t4) * qJD(6) -(t193 * t8 + t196 * t9) * t124 - t225 * t94 - t232 * t145 - t4 * t136 - t41 * t223 + t54 * t18 + t11 * t87 + t38 * t36 + g(1) * t117 - g(2) * t115 + (-t124 * t226 - t145 * t3) * qJD(6); 0, 0, 0, -t194 * t200 * t197, t260 * t200, t253, t252, qJDD(2), t194 * t208 - t290, g(3) * t194 + t197 * t208 (-t103 + t97) * t134 + (t104 - t96) * t131 + (-t189 * t98 - t191 * t99) * pkin(2), t96 * t103 - t97 * t104 + (-t290 + t189 * t43 + t191 * t42 + (-qJD(1) * t151 + t236) * t194) * pkin(2), -t167 * t98 + t171 * t99 + (-t103 - t89) * t134 + (t85 - t304) * t131, -t103 * qJD(2) + t80 * t131 + (-pkin(3) + t171) * qJDD(2) + t201, -t104 * qJD(2) + t167 * qJDD(2) - t73 * t131 + t80 * t134 + 0.2e1 * t184 + t204, -t39 * t167 + t40 * t171 - t73 * t80 - t85 * t103 - g(1) * (t198 * t237 + t156) - g(2) * (t195 * t237 + t154) - g(3) * t250 - t304 * t89, t164 * t281 + t22 * t131 + t167 * t78 + t262 * t109 + (-t190 * t261 - t29) * t134 + t206 * t188, -t164 * t282 - t23 * t131 + t167 * t79 + t262 * t111 + (t188 * t261 + t30) * t134 + t206 * t190, t30 * t109 + t29 * t111 + (qJD(5) * t111 - t134 * t23 - t164 * t79 - t6) * t190 + (qJD(5) * t109 + t134 * t22 - t164 * t78 - t7) * t188 - t213, t27 * t167 - t23 * t30 - t22 * t29 - g(1) * (-t198 * t296 + t156) - g(2) * (-t195 * t296 + t154) - g(3) * (t250 + t275) + t262 * t60 + t234 * t164 + t228 * qJD(5) + t288 * t309, t18 * t306 + t223 * t277, -t18 * t144 - t19 * t306 + t223 * t276 + t277 * t57, -t280 + t312, t238 - t283, t124 * t131 (-t129 * t193 + t130 * t196) * t94 + t150 * t19 + t11 * t144 + t3 * t131 + t263 * t57 + t276 * t38 + (t193 * t215 - t196 * t216) * t124 + t211 * t175 -(t129 * t196 + t130 * t193) * t94 + t150 * t18 + t11 * t306 - t4 * t131 - t263 * t223 - t277 * t38 + (t193 * t216 + t196 * t215) * t124 + t211 * t177; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t305, t97 * t131 + t96 * t134 + t205, t305, -0.2e1 * t134 * qJD(2) + t229, t123 - t99, t298 - t279 - t89 * t131 + (-qJD(4) - t85) * t134 + t205, t109 * t131 + t301, t111 * t131 - t224, t188 * t79 - t190 * t78 + (t109 * t188 + t111 * t190) * t134, t60 * t131 + t134 * t228 + t233 - t307, 0, 0, 0, 0, 0, t238 + t283, -t280 - t312; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t123 + t99, -t131 * t134 + qJDD(2), -t299 - t199, t89 * qJD(2) + t201 - t274, -qJD(2) * t109 + t224, -qJD(2) * t111 + t301, -t188 * t78 - t190 * t79 + (-t109 * t190 + t111 * t188) * t134, -t60 * qJD(2) + t134 * t227 + t213 + t234, 0, 0, 0, 0, 0, -qJD(2) * t57 + t312, qJD(2) * t223 + t238; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111 * t134 + t78, -t109 * t134 + t79, -t109 ^ 2 - t111 ^ 2, t23 * t109 + t22 * t111 + t184 + t204 + t302, 0, 0, 0, 0, 0, t19 - t311, t18 - t284; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t223 * t57, t223 ^ 2 - t57 ^ 2, t18 + t284, -t19 - t311, t94, -g(1) * t115 - g(2) * t117 + t177 * t169 + t223 * t38 + t303 * t4 + t247, g(1) * t116 - g(2) * t118 - t175 * t169 + t303 * t3 + t38 * t57 - t232;];
tau_reg  = t1;
