% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRPR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:47:16
% EndTime: 2019-03-09 01:47:22
% DurationCPUTime: 3.58s
% Computational Cost: add. (6127->471), mult. (11414->583), div. (0->0), fcn. (7465->12), ass. (0->231)
t168 = cos(qJ(4));
t267 = cos(pkin(10));
t223 = t267 * t168;
t118 = qJD(1) * t223;
t160 = sin(pkin(10));
t166 = sin(qJ(4));
t251 = qJD(1) * t166;
t95 = t160 * t251 - t118;
t260 = qJD(6) - t95;
t165 = sin(qJ(6));
t167 = cos(qJ(6));
t246 = t166 * qJD(3);
t169 = -pkin(1) - pkin(2);
t115 = t169 * qJD(1) + qJD(2);
t161 = sin(pkin(9));
t162 = cos(pkin(9));
t252 = qJD(1) * t162;
t89 = qJ(2) * t252 + t161 * t115;
t82 = -qJD(1) * pkin(7) + t89;
t303 = (qJ(5) * qJD(1) - t82) * t168;
t56 = t246 - t303;
t46 = t267 * t56;
t274 = t166 * t82;
t63 = t168 * qJD(3) - t274;
t55 = qJ(5) * t251 + t63;
t49 = qJD(4) * pkin(4) + t55;
t25 = t160 * t49 + t46;
t23 = qJD(4) * pkin(8) + t25;
t151 = t168 * pkin(4);
t253 = qJD(1) * t161;
t88 = -qJ(2) * t253 + t115 * t162;
t81 = qJD(1) * pkin(3) - t88;
t70 = qJD(1) * t151 + qJD(5) + t81;
t101 = t160 * t168 + t267 * t166;
t97 = t101 * qJD(1);
t35 = -pkin(5) * t95 + pkin(8) * t97 + t70;
t9 = -t165 * t23 + t167 * t35;
t309 = t260 * t9;
t186 = -t160 * t166 + t223;
t96 = t101 * qJD(4);
t280 = t161 * t96 + t186 * t252;
t10 = t165 * t35 + t167 * t23;
t308 = t10 * t260;
t243 = qJD(1) * qJD(2);
t122 = t162 * t243;
t114 = t169 * qJDD(1) + qJDD(2);
t238 = qJDD(1) * t162;
t259 = qJ(2) * t238 + t161 * t114;
t72 = t122 + t259;
t69 = -qJDD(1) * pkin(7) + t72;
t215 = -qJD(3) * qJD(4) - t69;
t221 = t165 * t260;
t242 = qJD(1) * qJD(4);
t225 = t166 * t242;
t234 = t168 * qJDD(1);
t307 = t225 - t234;
t220 = t167 * t260;
t235 = t166 * qJDD(1);
t197 = -qJDD(1) * t223 + t160 * t235;
t58 = qJD(1) * t96 + t197;
t57 = -qJDD(6) + t58;
t278 = t165 * t57;
t306 = t220 * t260 - t278;
t248 = qJD(6) * t165;
t304 = t167 * t57 + t248 * t260;
t156 = qJ(4) + pkin(10);
t139 = sin(t156);
t296 = sin(qJ(1));
t297 = cos(qJ(1));
t102 = t297 * t161 - t296 * t162;
t100 = -t296 * t161 - t297 * t162;
t294 = g(2) * t100;
t207 = g(1) * t102 - t294;
t302 = t207 * t139;
t249 = qJD(4) * t166;
t98 = qJD(4) * t223 - t160 * t249;
t293 = g(2) * t102;
t206 = g(1) * t100 + t293;
t129 = pkin(4) * t160 + pkin(8);
t140 = cos(t156);
t180 = g(3) * t140 - t139 * t206;
t146 = t168 * qJDD(3);
t241 = qJD(1) * qJD(5);
t20 = qJDD(4) * pkin(4) + t146 + qJD(4) * t303 + (qJ(5) * qJDD(1) + t215 + t241) * t166;
t231 = -t166 * qJDD(3) + t215 * t168;
t33 = -t249 * t82 - t231;
t21 = t307 * qJ(5) - t168 * t241 + t33;
t7 = -t160 * t21 + t267 * t20;
t5 = -qJDD(4) * pkin(5) - t7;
t301 = -qJD(6) * t129 * t260 + t180 - t5;
t247 = qJD(6) * t167;
t275 = t165 * t98;
t300 = -t101 * (-t247 * t260 + t278) + t260 * t275;
t292 = g(3) * t139;
t299 = -t206 * t140 - t292;
t270 = t168 * t82;
t64 = t246 + t270;
t34 = -t64 * qJD(4) - t166 * t69 + t146;
t173 = -(t166 * t64 + t168 * t63) * qJD(4) - t34 * t166 + t33 * t168;
t178 = -qJDD(1) * t101 + t160 * t225;
t298 = t97 ^ 2;
t290 = g(3) * t168;
t245 = t167 * qJD(4);
t73 = -t165 * t97 - t245;
t289 = t73 * t95;
t288 = t73 * t97;
t75 = qJD(4) * t165 - t167 * t97;
t287 = t75 * t73;
t286 = t75 * t97;
t285 = t95 * t97;
t8 = t160 * t20 + t267 * t21;
t59 = qJD(4) * t118 - t178;
t222 = -t167 * qJDD(4) - t165 * t59;
t32 = qJD(6) * t75 + t222;
t284 = -t165 * t32 - t73 * t247;
t87 = t186 * t161;
t65 = -t162 * t167 - t165 * t87;
t283 = qJD(6) * t65 - t165 * t253 - t280 * t167;
t192 = t162 * t165 - t167 * t87;
t282 = qJD(6) * t192 + t280 * t165 - t167 * t253;
t281 = -t98 * t161 + t162 * t97;
t279 = t160 * t56;
t277 = t165 * t73;
t276 = t165 * t75;
t273 = t167 * t73;
t272 = t167 * t75;
t271 = t167 * t98;
t31 = -qJD(6) * t245 - t165 * qJDD(4) + t167 * t59 - t248 * t97;
t269 = t31 * t165;
t268 = t32 * t167;
t266 = pkin(1) * qJDD(1);
t265 = qJD(1) * t81;
t264 = t140 * t165;
t263 = t140 * t167;
t171 = qJD(1) ^ 2;
t262 = t162 * t171;
t107 = t162 * qJ(2) + t161 * t169;
t104 = -pkin(7) + t107;
t261 = qJ(5) - t104;
t258 = t297 * pkin(1) + t296 * qJ(2);
t257 = g(1) * t296 - g(2) * t297;
t157 = t166 ^ 2;
t158 = t168 ^ 2;
t256 = t157 - t158;
t255 = t157 + t158;
t170 = qJD(4) ^ 2;
t254 = t170 + t171;
t250 = qJD(2) * t162;
t244 = qJ(2) * qJDD(1);
t239 = qJDD(1) * t161;
t237 = qJDD(4) * t166;
t236 = qJDD(4) * t168;
t232 = 0.2e1 * t243;
t228 = t166 * t171 * t168;
t227 = t297 * pkin(2) + t258;
t224 = t168 * t242;
t120 = t161 * t243;
t219 = t261 * t166;
t217 = -qJ(2) * t239 + t114 * t162;
t106 = -t161 * qJ(2) + t162 * t169;
t216 = qJD(4) * t261;
t214 = qJDD(1) * t255;
t213 = qJDD(2) - t266;
t212 = -qJD(5) + t250;
t211 = t166 * t224;
t103 = pkin(3) - t106;
t210 = -t296 * pkin(1) + t297 * qJ(2);
t109 = -pkin(4) * t249 + t161 * qJD(2);
t208 = pkin(5) * t140 + pkin(8) * t139;
t205 = t186 * t31 + t96 * t75;
t204 = -t186 * t32 + t96 * t73;
t203 = t186 * t59 - t96 * t97;
t71 = -t120 + t217;
t135 = t151 + pkin(3);
t164 = -qJ(5) - pkin(7);
t202 = -t100 * t135 - t102 * t164 + t227;
t201 = t10 * t167 - t165 * t9;
t200 = t10 * t165 + t167 * t9;
t199 = -t101 * t58 - t95 * t98;
t198 = -qJD(6) * t23 + t294;
t196 = t161 * t88 - t162 * t89;
t83 = t261 * t168;
t39 = t160 * t219 - t267 * t83;
t91 = t103 + t151;
t41 = pkin(5) * t186 + pkin(8) * t101 + t91;
t16 = t165 * t41 + t167 * t39;
t15 = -t165 * t39 + t167 * t41;
t193 = t166 * t63 - t168 * t64;
t191 = qJD(4) * t96 - qJDD(4) * t186;
t190 = t95 * t221 - t304;
t189 = -qJD(4) * t98 - qJDD(4) * t101;
t68 = qJDD(1) * pkin(3) - t71;
t24 = t267 * t49 - t279;
t22 = -qJD(4) * pkin(5) - t24;
t187 = t129 * t57 + t260 * t22;
t185 = g(1) * t297 + g(2) * t296;
t183 = -t206 + t265;
t182 = -t296 * pkin(2) + t210;
t6 = qJDD(4) * pkin(8) + t8;
t181 = -qJD(6) * t35 - t140 * t293 - t292 - t6;
t179 = -t100 * t164 + t102 * t135 + t182;
t177 = -t166 * t212 + t168 * t216;
t176 = -qJDD(4) * t104 + (-qJD(1) * t103 - t250 - t81) * qJD(4);
t175 = t304 * t101 - t260 * t271;
t42 = -t307 * pkin(4) + qJDD(5) + t68;
t14 = -pkin(5) * t58 + pkin(8) * t59 + t42;
t1 = qJD(6) * t9 + t165 * t14 + t167 * t6;
t13 = t167 * t14;
t2 = -qJD(6) * t10 - t165 * t6 + t13;
t174 = -qJD(6) * t200 + t1 * t167 - t2 * t165;
t172 = qJDD(1) * t103 - t104 * t170 + t120 - t207 + t68;
t130 = -t267 * pkin(4) - pkin(5);
t112 = -t166 * t170 + t236;
t111 = -t168 * t170 - t237;
t94 = t95 ^ 2;
t86 = t101 * t161;
t54 = t166 * t216 + t168 * t212;
t51 = -t100 * t263 + t102 * t165;
t50 = t100 * t264 + t102 * t167;
t44 = -pkin(4) * t251 - pkin(5) * t97 - pkin(8) * t95;
t40 = -pkin(5) * t96 + pkin(8) * t98 + t109;
t38 = -t160 * t83 - t267 * t219;
t29 = t267 * t55 - t279;
t28 = t160 * t177 + t267 * t54;
t27 = t160 * t55 + t46;
t26 = t160 * t54 - t267 * t177;
t12 = t165 * t44 + t167 * t29;
t11 = -t165 * t29 + t167 * t44;
t4 = -qJD(6) * t16 - t165 * t28 + t167 * t40;
t3 = qJD(6) * t15 + t165 * t40 + t167 * t28;
t17 = [0, 0, 0, 0, 0, qJDD(1), t257, t185, 0, 0, 0, 0, 0, qJDD(1), 0, 0, -qJDD(2) + t257 + 0.2e1 * t266, 0, -t185 + t232 + 0.2e1 * t244, -t213 * pkin(1) - g(1) * t210 - g(2) * t258 + (t232 + t244) * qJ(2), 0, 0, 0, 0, 0, qJDD(1), -qJDD(1) * t106 + 0.2e1 * t120 - t207 - t217, qJDD(1) * t107 + 0.2e1 * t122 + t206 + t259, 0, -g(1) * t182 - g(2) * t227 - qJD(2) * t196 + t71 * t106 + t72 * t107, qJDD(1) * t157 + 0.2e1 * t211, 0.2e1 * t166 * t234 - 0.2e1 * t242 * t256, t111, qJDD(1) * t158 - 0.2e1 * t211, -t112, 0, t166 * t176 + t168 * t172, -t166 * t172 + t168 * t176, -t104 * t214 - t122 * t255 - t173 - t206, t68 * t103 - g(1) * (t102 * pkin(3) + t100 * pkin(7) + t182) - g(2) * (-pkin(3) * t100 + pkin(7) * t102 + t227) + (t81 * t161 - t162 * t193) * qJD(2) + t173 * t104, t101 * t59 + t97 * t98, t199 + t203, t189, -t186 * t58 + t95 * t96, t191, 0, -t26 * qJD(4) - t38 * qJDD(4) - t109 * t95 - t140 * t207 + t186 * t42 - t91 * t58 - t70 * t96, -t28 * qJD(4) - t39 * qJDD(4) - t42 * t101 - t109 * t97 - t91 * t59 - t70 * t98 + t302, t101 * t7 - t186 * t8 + t24 * t98 + t25 * t96 - t26 * t97 + t28 * t95 - t38 * t59 + t39 * t58 - t206, -g(1) * t179 - g(2) * t202 + t70 * t109 - t24 * t26 + t25 * t28 - t7 * t38 + t8 * t39 + t42 * t91, -t75 * t271 + (t167 * t31 + t248 * t75) * t101 (t273 + t276) * t98 + (-t269 + t268 + (t272 - t277) * qJD(6)) * t101, t175 - t205, t284 * t101 - t73 * t275, t204 + t300, -t186 * t57 - t260 * t96, t4 * t260 - t15 * t57 + t2 * t186 - t9 * t96 + t26 * t73 + t38 * t32 - t22 * t275 - g(1) * (t100 * t165 + t102 * t263) - g(2) * t51 + (-t5 * t165 - t22 * t247) * t101, -t3 * t260 + t16 * t57 - t1 * t186 + t10 * t96 + t26 * t75 - t38 * t31 - t22 * t271 - g(1) * (t100 * t167 - t102 * t264) - g(2) * t50 + (-t5 * t167 + t22 * t248) * t101, t15 * t31 - t16 * t32 - t3 * t73 - t4 * t75 + t200 * t98 - t302 + (qJD(6) * t201 + t1 * t165 + t2 * t167) * t101, t1 * t16 + t10 * t3 + t2 * t15 + t9 * t4 + t5 * t38 + t22 * t26 - g(1) * (t102 * t208 + t179) - g(2) * (-t100 * t208 + t202); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), 0, -t171, -qJ(2) * t171 + t213 - t257, 0, 0, 0, 0, 0, 0, -t161 * t171 - t238, t239 - t262, 0, qJD(1) * t196 + t72 * t161 + t71 * t162 - t257, 0, 0, 0, 0, 0, 0 (0.2e1 * t225 - t234) * t162 + (-t168 * t254 - t237) * t161 (0.2e1 * t224 + t235) * t162 + (t166 * t254 - t236) * t161, -t161 * t214 + t255 * t262 (qJD(1) * t193 - t68) * t162 + (t173 - t265) * t161 - t257, 0, 0, 0, 0, 0, 0, t281 * qJD(4) - t86 * qJDD(4) + t162 * t58 + t95 * t253, t280 * qJD(4) - t87 * qJDD(4) + t162 * t59 + t97 * t253, -t280 * t95 + t281 * t97 + t87 * t58 - t86 * t59, -t42 * t162 + t281 * t24 - t280 * t25 - t70 * t253 - t7 * t86 + t8 * t87 - t257, 0, 0, 0, 0, 0, 0, t260 * t282 - t281 * t73 + t86 * t32 - t65 * t57, -t192 * t57 - t260 * t283 - t281 * t75 - t86 * t31, t192 * t32 - t282 * t75 - t283 * t73 + t65 * t31, -t1 * t192 + t10 * t283 + t2 * t65 - t22 * t281 + t282 * t9 + t5 * t86 - t257; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) + g(3), 0, 0, 0, 0, 0, 0, t112, t111, 0, -qJD(4) * t193 + t33 * t166 + t34 * t168 + g(3), 0, 0, 0, 0, 0, 0, -t191, t189, -t199 + t203, t101 * t8 + t186 * t7 - t24 * t96 + t25 * t98 + g(3), 0, 0, 0, 0, 0, 0, t204 - t300, t175 + t205 (-t273 + t276) * t98 + (-t269 - t268 + (t272 + t277) * qJD(6)) * t101, t101 * t174 - t186 * t5 + t201 * t98 + t22 * t96 + g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t228, t256 * t171, -t235, t228, -t234, qJDD(4), t290 + t146 + (t64 - t270) * qJD(4) + (t183 + t215) * t166, -g(3) * t166 + (t63 + t274) * qJD(4) + t183 * t168 + t231, 0, 0, t285, -t94 + t298 (-t118 - t95) * qJD(4) + t178, -t285, t197, qJDD(4), t27 * qJD(4) + t70 * t97 + (qJDD(4) * t267 - t251 * t95) * pkin(4) + t180 + t7, t29 * qJD(4) - t70 * t95 + (-qJDD(4) * t160 - t251 * t97) * pkin(4) - t8 + t299 (-t25 + t27) * t97 + (t24 - t29) * t95 + (t160 * t58 + t267 * t59) * pkin(4), t24 * t27 - t25 * t29 + (t267 * t7 + t290 + t160 * t8 + (qJD(1) * t70 - t206) * t166) * pkin(4), t220 * t75 - t269 (-t31 + t289) * t167 - t260 * t276 + t284, t286 + t306, t221 * t73 - t268, t190 - t288, t260 * t97, -t11 * t260 + t130 * t32 + t187 * t165 + t167 * t301 - t27 * t73 + t9 * t97, -t10 * t97 + t12 * t260 - t130 * t31 - t165 * t301 + t187 * t167 - t27 * t75, t11 * t75 + t12 * t73 + (-t129 * t32 + t9 * t95 + t1 + (t129 * t75 - t9) * qJD(6)) * t167 + (t10 * t95 - t129 * t31 - t2 + (t129 * t73 - t10) * qJD(6)) * t165 - t299, t5 * t130 - t10 * t12 - t9 * t11 - t22 * t27 - g(3) * (-t151 - t208) + t174 * t129 - t206 * (pkin(4) * t166 + pkin(5) * t139 - pkin(8) * t140); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t97 * qJD(4) - t197 (-t118 + t95) * qJD(4) + t178, -t94 - t298, -t24 * t97 - t25 * t95 - t207 + t42, 0, 0, 0, 0, 0, 0, t190 + t288, t286 - t306 (t31 + t289) * t167 + t75 * t221 + t284, t22 * t97 + (t2 + t308) * t167 + (t1 - t309) * t165 - t207; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t287, -t73 ^ 2 + t75 ^ 2, t260 * t73 - t31, -t287, -t222 + (-qJD(6) + t260) * t75, -t57, -g(1) * t50 + t165 * t181 + t167 * t198 - t22 * t75 + t13 + t308, g(1) * t51 + t22 * t73 + t309 + (-t14 - t198) * t165 + t181 * t167, 0, 0;];
tau_reg  = t17;
