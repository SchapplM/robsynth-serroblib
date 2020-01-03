% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP7_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:01:37
% EndTime: 2019-12-31 20:01:44
% DurationCPUTime: 3.66s
% Computational Cost: add. (5082->469), mult. (11887->567), div. (0->0), fcn. (8439->10), ass. (0->222)
t154 = cos(qJ(4));
t149 = sin(pkin(8));
t155 = cos(qJ(2));
t152 = sin(qJ(2));
t248 = cos(pkin(8));
t211 = t248 * t152;
t112 = t149 * t155 + t211;
t103 = t112 * qJD(1);
t151 = sin(qJ(4));
t210 = t248 * t155;
t129 = qJD(1) * t210;
t228 = qJD(1) * qJD(2);
t215 = t152 * t228;
t169 = qJDD(1) * t112 - t149 * t215;
t162 = qJD(2) * t129 + t169;
t227 = qJD(2) * qJD(4);
t230 = qJD(4) * t151;
t40 = -t151 * qJDD(2) + t103 * t230 + (-t162 - t227) * t154;
t86 = qJD(2) * t151 + t103 * t154;
t80 = t86 * t230;
t295 = -t154 * t40 - t80;
t231 = qJD(1) * t152;
t100 = t149 * t231 - t129;
t92 = qJD(4) + t100;
t180 = -t149 * t152 + t210;
t105 = t180 * qJD(2);
t84 = -t154 * qJD(2) + t103 * t151;
t251 = t154 * t84;
t254 = t151 * t86;
t190 = t251 + t254;
t161 = -t154 * qJDD(2) + t151 * t162;
t41 = qJD(4) * t86 + t161;
t252 = t154 * t41;
t255 = t151 * t40;
t297 = t112 * ((t151 * t84 - t154 * t86) * qJD(4) - t252 + t255) - t190 * t105;
t229 = qJD(4) * t154;
t102 = t112 * qJD(2);
t226 = t152 * qJDD(1);
t196 = -qJDD(1) * t210 + t149 * t226;
t75 = qJD(1) * t102 + t196;
t70 = qJDD(4) + t75;
t60 = t151 * t70;
t183 = t92 * t229 + t60;
t61 = t154 * t70;
t296 = -t92 * t230 + t61;
t266 = -t151 * t41 - t84 * t229;
t294 = (-t251 + t254) * t100 + t266 - t295;
t146 = qJ(2) + pkin(8);
t140 = sin(t146);
t153 = sin(qJ(1));
t156 = cos(qJ(1));
t200 = g(1) * t156 + g(2) * t153;
t292 = t200 * t140;
t141 = cos(t146);
t291 = t141 * pkin(3) + t140 * pkin(7);
t278 = g(1) * t153;
t290 = -g(2) * t156 + t278;
t289 = t103 * qJD(2);
t267 = qJ(3) + pkin(6);
t119 = t267 * t155;
t115 = qJD(1) * t119;
t106 = t149 * t115;
t216 = t267 * t152;
t114 = qJD(1) * t216;
t262 = qJD(2) * pkin(2);
t109 = -t114 + t262;
t73 = t248 * t109 - t106;
t63 = -qJD(2) * pkin(3) - t73;
t24 = t84 * pkin(4) - t86 * qJ(5) + t63;
t280 = pkin(2) * t149;
t134 = pkin(7) + t280;
t258 = t134 * t70;
t288 = t24 * t92 - t258;
t244 = t105 * t151;
t287 = t102 * t84 + t112 * t183 - t180 * t41 + t92 * t244;
t275 = g(3) * t140;
t286 = t200 * t141 + t275;
t272 = t155 * pkin(2);
t139 = pkin(1) + t272;
t72 = -pkin(3) * t180 - pkin(7) * t112 - t139;
t83 = t248 * t119 - t149 * t216;
t265 = t151 * t72 + t154 * t83;
t213 = qJD(2) * t267;
t177 = -qJD(3) * t152 - t155 * t213;
t97 = qJD(3) * t155 - t152 * t213;
t54 = t149 * t177 + t248 * t97;
t223 = t152 * t262;
t56 = pkin(3) * t102 - pkin(7) * t105 + t223;
t11 = -qJD(4) * t265 - t151 * t54 + t154 * t56;
t284 = t86 ^ 2;
t283 = t92 ^ 2;
t282 = t103 ^ 2;
t281 = pkin(4) * t70;
t279 = pkin(2) * t152;
t274 = g(3) * t141;
t273 = g(3) * t155;
t118 = -t139 * qJD(1) + qJD(3);
t48 = pkin(3) * t100 - pkin(7) * t103 + t118;
t212 = t248 * t115;
t74 = t149 * t109 + t212;
t64 = qJD(2) * pkin(7) + t74;
t23 = t151 * t48 + t154 * t64;
t17 = qJ(5) * t92 + t23;
t271 = t17 * t92;
t22 = -t151 * t64 + t154 * t48;
t270 = t22 * t92;
t269 = t23 * t92;
t207 = t84 * t92;
t268 = t86 * t84;
t69 = qJDD(2) * pkin(2) + t177 * qJD(1) - qJDD(1) * t216;
t76 = t97 * qJD(1) + qJDD(1) * t119;
t36 = -t149 * t76 + t248 * t69;
t37 = t149 * t69 + t248 * t76;
t55 = pkin(2) * t231 + pkin(3) * t103 + pkin(7) * t100;
t78 = -t248 * t114 - t106;
t32 = t151 * t55 + t154 * t78;
t197 = pkin(4) * t151 - qJ(5) * t154;
t77 = -t114 * t149 + t212;
t264 = -qJD(5) * t151 + t92 * t197 - t77;
t263 = qJ(5) * t70;
t261 = t103 * t84;
t260 = t103 * t86;
t259 = t134 * t40;
t257 = t134 * t84;
t256 = t134 * t86;
t250 = t154 * t92;
t249 = t92 * t103;
t247 = pkin(6) * qJDD(1);
t246 = qJDD(1) * pkin(1);
t245 = t103 * t100;
t243 = t105 * t154;
t242 = t140 * t156;
t241 = t141 * t156;
t240 = t267 * t156;
t239 = t151 * t153;
t238 = t153 * t154;
t237 = t154 * t156;
t236 = t156 * t151;
t235 = qJD(5) - t22;
t234 = (g(1) * t237 + g(2) * t238) * t140;
t147 = t152 ^ 2;
t148 = t155 ^ 2;
t233 = t147 - t148;
t232 = t147 + t148;
t225 = t155 * qJDD(1);
t224 = pkin(2) * t215 + qJDD(3);
t222 = t84 ^ 2 - t284;
t221 = qJD(4) * t134 * t92;
t158 = qJD(1) ^ 2;
t219 = t152 * t158 * t155;
t127 = t156 * t139;
t218 = pkin(3) * t241 + pkin(7) * t242 + t127;
t217 = t248 * pkin(2);
t53 = t149 * t97 - t248 * t177;
t25 = -pkin(2) * t225 + t75 * pkin(3) - pkin(7) * t162 + t224 - t246;
t34 = qJDD(2) * pkin(7) + t37;
t209 = t151 * t34 - t154 * t25 + t64 * t229 + t48 * t230;
t82 = t119 * t149 + t267 * t211;
t206 = t155 * t215;
t93 = t141 * t239 + t237;
t95 = t141 * t236 - t238;
t205 = g(1) * t93 - g(2) * t95;
t94 = t141 * t238 - t236;
t96 = t141 * t237 + t239;
t204 = g(1) * t94 - g(2) * t96;
t203 = -g(2) * t242 + t140 * t278;
t202 = t272 + t291;
t33 = -qJDD(2) * pkin(3) - t36;
t201 = -pkin(3) * t140 - t279;
t198 = t154 * pkin(4) + t151 * qJ(5);
t16 = -pkin(4) * t92 + t235;
t195 = -t151 * t17 + t154 * t16;
t194 = -t151 * t23 - t154 * t22;
t31 = -t151 * t78 + t154 * t55;
t42 = -t151 * t83 + t154 * t72;
t188 = t100 * t250 + t183;
t187 = -t151 * t100 * t92 + t296;
t186 = pkin(3) + t198;
t184 = t221 + t274;
t182 = -0.2e1 * pkin(1) * t228 - pkin(6) * qJDD(2);
t3 = t151 * t25 + t154 * t34 + t48 * t229 - t64 * t230;
t10 = t151 * t56 + t154 * t54 + t72 * t229 - t83 * t230;
t181 = t92 * t63 - t258;
t179 = t184 + t33;
t178 = t151 * t207 - t252;
t91 = -t139 * qJDD(1) + t224;
t176 = -t187 - t261;
t174 = g(1) * t95 + g(2) * t93 + t151 * t275 - t209;
t173 = -t274 + t292;
t157 = qJD(2) ^ 2;
t172 = -pkin(6) * t157 + 0.2e1 * t246 + t290;
t171 = pkin(1) * t158 + t200 - t247;
t170 = -t134 * t252 - t286;
t167 = (-g(1) * (-t139 - t291) - g(2) * t267) * t153;
t166 = -t266 * t112 + t84 * t244;
t165 = t24 * t86 + qJDD(5) - t174;
t164 = -g(1) * t96 - g(2) * t94 - t154 * t275 + t3;
t160 = -t103 * t229 - t151 * t227 + t86 * t92 - t161;
t135 = -t217 - pkin(3);
t123 = pkin(7) * t241;
t121 = t153 * t141 * pkin(7);
t110 = -t217 - t186;
t99 = t100 ^ 2;
t45 = pkin(4) * t86 + qJ(5) * t84;
t44 = t112 * t197 + t82;
t28 = pkin(4) * t180 - t42;
t27 = -qJ(5) * t180 + t265;
t26 = t102 * t92 - t180 * t70;
t19 = -pkin(4) * t103 - t31;
t18 = qJ(5) * t103 + t32;
t15 = -t40 + t207;
t14 = t197 * t105 + (qJD(4) * t198 - qJD(5) * t154) * t112 + t53;
t13 = t188 - t260;
t12 = t86 * t250 - t255;
t9 = t112 * t295 + t86 * t243;
t8 = -pkin(4) * t102 - t11;
t7 = qJ(5) * t102 - qJD(5) * t180 + t10;
t6 = pkin(4) * t41 + qJ(5) * t40 - qJD(5) * t86 + t33;
t5 = t102 * t86 + t112 * t296 + t180 * t40 + t92 * t243;
t2 = qJDD(5) + t209 - t281;
t1 = qJD(5) * t92 + t263 + t3;
t4 = [0, 0, 0, 0, 0, qJDD(1), t290, t200, 0, 0, qJDD(1) * t147 + 0.2e1 * t206, 0.2e1 * t152 * t225 - 0.2e1 * t233 * t228, qJDD(2) * t152 + t155 * t157, qJDD(1) * t148 - 0.2e1 * t206, qJDD(2) * t155 - t152 * t157, 0, t152 * t182 + t155 * t172, -t152 * t172 + t155 * t182, 0.2e1 * t232 * t247 - t200, -g(1) * (-pkin(1) * t153 + pkin(6) * t156) - g(2) * (pkin(1) * t156 + pkin(6) * t153) + (t232 * pkin(6) ^ 2 + pkin(1) ^ 2) * qJDD(1), t103 * t105 + t112 * t162, -t105 * t100 - t103 * t102 - t112 * t75 + t162 * t180, qJD(2) * t105 + qJDD(2) * t112, t100 * t102 - t180 * t75, -qJD(2) * t102 + qJDD(2) * t180, 0, -qJDD(2) * t82 + t102 * t118 - t180 * t91 - t139 * t75 + t290 * t141 + (t100 * t279 - t53) * qJD(2), -t54 * qJD(2) - t83 * qJDD(2) + t103 * t223 + t118 * t105 + t91 * t112 - t139 * t162 - t203, -t54 * t100 - t74 * t102 + t53 * t103 - t73 * t105 - t36 * t112 + t162 * t82 + t180 * t37 - t83 * t75 - t200, t37 * t83 + t74 * t54 - t36 * t82 - t73 * t53 - t91 * t139 + t118 * t223 - g(1) * (-t139 * t153 + t240) - g(2) * (t153 * t267 + t127), t9, t297, t5, t166, -t287, t26, t63 * t244 + t102 * t22 + t11 * t92 + t180 * t209 + t41 * t82 + t42 * t70 + t53 * t84 + (t151 * t33 + t229 * t63) * t112 + t204, t63 * t243 - t10 * t92 - t102 * t23 + t180 * t3 - t40 * t82 - t265 * t70 + t53 * t86 + (t154 * t33 - t230 * t63) * t112 - t205, -t10 * t84 - t11 * t86 + t40 * t42 - t41 * t265 + t194 * t105 + (-t151 * t3 + t154 * t209 + (t151 * t22 - t154 * t23) * qJD(4)) * t112 + t203, -g(1) * t240 - g(2) * t218 + t23 * t10 + t22 * t11 - t209 * t42 + t265 * t3 + t33 * t82 + t63 * t53 + t167, t9, t5, -t297, t26, t287, t166, t24 * t244 - t102 * t16 + t180 * t2 + t14 * t84 - t28 * t70 + t41 * t44 - t8 * t92 + (t151 * t6 + t229 * t24) * t112 + t204, -t27 * t41 - t28 * t40 - t7 * t84 + t8 * t86 + t195 * t105 + (-t1 * t151 + t154 * t2 + (-t151 * t16 - t154 * t17) * qJD(4)) * t112 + t203, -t24 * t243 - t1 * t180 + t102 * t17 - t14 * t86 + t27 * t70 + t40 * t44 + t7 * t92 + (-t154 * t6 + t230 * t24) * t112 + t205, t1 * t27 + t17 * t7 + t6 * t44 + t24 * t14 + t2 * t28 + t16 * t8 - g(1) * (-pkin(4) * t94 - qJ(5) * t93 + t240) - g(2) * (pkin(4) * t96 + qJ(5) * t95 + t218) + t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t219, t233 * t158, t226, t219, t225, qJDD(2), t152 * t171 - t273, g(3) * t152 + t155 * t171, 0, 0, t245, -t99 + t282, (t129 + t100) * qJD(2) + t169, -t245, -t196, qJDD(2), t77 * qJD(2) - t118 * t103 + (t248 * qJDD(2) - t100 * t231) * pkin(2) + t173 + t36, qJD(2) * t78 + t100 * t118 + (-qJDD(2) * t149 - t103 * t231) * pkin(2) - t37 + t286, -t162 * t217 - t75 * t280 - (-t74 + t77) * t103 + (-t73 + t78) * t100, t73 * t77 - t74 * t78 + (t248 * t36 - t273 + t149 * t37 + (-qJD(1) * t118 + t200) * t152) * pkin(2), t12, -t100 * t190 + t266 + t295, t13, t178, -t176, -t249, -t103 * t22 + t135 * t41 + t151 * t181 - t154 * t179 - t31 * t92 - t77 * t84 + t234, t103 * t23 - t135 * t40 + t32 * t92 - t77 * t86 + t181 * t154 + (t179 - t292) * t151, t31 * t86 + t32 * t84 + (-t100 * t22 + t3 + (-t22 + t256) * qJD(4)) * t154 + (-t100 * t23 - t259 + t209 + (-t23 + t257) * qJD(4)) * t151 + t170, t33 * t135 - t23 * t32 - t22 * t31 - t63 * t77 - g(1) * (t156 * t201 + t123) - g(2) * (t153 * t201 + t121) - g(3) * t202 + (qJD(4) * t194 + t151 * t209 + t3 * t154) * t134, t12, t13, t80 + (t100 * t86 + t41) * t151 + (t40 + t207) * t154, -t249, t176, t178, t103 * t16 + t110 * t41 + t19 * t92 + t264 * t84 + (-t184 - t6) * t154 + t288 * t151 + t234, t18 * t84 - t19 * t86 + (t100 * t16 + t1 + (t16 + t256) * qJD(4)) * t154 + (-t100 * t17 - t259 + t2 + (-t17 + t257) * qJD(4)) * t151 + t170, -t103 * t17 + t110 * t40 - t18 * t92 - t264 * t86 - t288 * t154 + (t173 - t6 - t221) * t151, t6 * t110 - t17 * t18 - t16 * t19 - g(1) * (-t156 * t279 + t123) - g(2) * (-t153 * t279 + t121) - g(3) * (t141 * t198 + t202) + t264 * t24 + (qJD(4) * t195 + t1 * t154 + t2 * t151) * t134 + t186 * t292; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t196 + 0.2e1 * t289, (t129 - t100) * qJD(2) + t169, -t99 - t282, t100 * t74 + t103 * t73 - t290 + t91, 0, 0, 0, 0, 0, 0, t187 - t261, -t154 * t283 - t260 - t60, t294, -t103 * t63 + (-t209 + t269) * t154 + (t3 - t270) * t151 - t290, 0, 0, 0, 0, 0, 0, -t151 * t283 - t261 + t61, t294, t188 + t260, -t103 * t24 + (-t2 + t271) * t154 + (t16 * t92 + t1) * t151 - t290; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t268, -t222, t15, -t268, t160, t70, -t63 * t86 + t174 + t269, t63 * t84 - t164 + t270, 0, 0, t268, t15, t222, t70, -t160, -t268, -t45 * t84 - t165 + t269 + 0.2e1 * t281, pkin(4) * t40 - qJ(5) * t41 + (t17 - t23) * t86 + (t16 - t235) * t84, 0.2e1 * t263 - t24 * t84 + t45 * t86 + (0.2e1 * qJD(5) - t22) * t92 + t164, t1 * qJ(5) - t2 * pkin(4) - t24 * t45 - t16 * t23 - g(1) * (-pkin(4) * t95 + qJ(5) * t96) - g(2) * (-pkin(4) * t93 + qJ(5) * t94) + t235 * t17 + t197 * t275; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(4) - t196 + t268 - t289, t15, -t283 - t284, t165 - t271 - t281;];
tau_reg = t4;
