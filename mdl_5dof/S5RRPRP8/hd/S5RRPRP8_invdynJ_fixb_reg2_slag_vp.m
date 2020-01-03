% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRP8
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_invdynJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:04:32
% EndTime: 2019-12-31 20:04:37
% DurationCPUTime: 3.04s
% Computational Cost: add. (2498->376), mult. (5478->420), div. (0->0), fcn. (3462->6), ass. (0->201)
t169 = qJD(2) - qJD(4);
t181 = cos(qJ(4));
t182 = cos(qJ(2));
t178 = sin(qJ(4));
t179 = sin(qJ(2));
t253 = t178 * t179;
t97 = t181 * t182 + t253;
t81 = t97 * qJD(1);
t259 = t169 * t81;
t231 = qJD(1) * qJD(2);
t221 = t179 * t231;
t230 = t182 * qJDD(1);
t281 = t221 - t230;
t156 = t179 * qJDD(1);
t220 = t182 * t231;
t282 = -t220 - t156;
t26 = qJD(4) * t81 - t178 * t281 + t181 * t282;
t12 = t26 + t259;
t232 = qJD(4) * t182;
t233 = qJD(4) * t181;
t235 = qJD(2) * t182;
t283 = t179 * t233 + (-t232 + t235) * t178;
t243 = t182 * pkin(2) + t179 * qJ(3);
t280 = -pkin(1) - t243;
t273 = pkin(2) + pkin(3);
t109 = qJ(3) * t181 - t178 * t273;
t168 = qJDD(2) - qJDD(4);
t108 = -qJ(3) * t178 - t181 * t273;
t69 = qJD(3) * t181 + t108 * qJD(4);
t279 = t109 * t168 + t69 * t169;
t27 = t283 * qJD(1) + t97 * qJDD(1) - t181 * t221;
t278 = -t27 * qJ(5) - t81 * qJD(5);
t237 = qJD(1) * t182;
t238 = qJD(1) * t179;
t83 = -t178 * t237 + t181 * t238;
t258 = t169 * t83;
t14 = t27 + t258;
t274 = t83 ^ 2;
t79 = t81 ^ 2;
t32 = -t79 + t274;
t180 = sin(qJ(1));
t183 = cos(qJ(1));
t276 = g(1) * t183 + g(2) * t180;
t257 = pkin(6) * qJDD(2);
t86 = -qJD(1) * pkin(1) - pkin(2) * t237 - qJ(3) * t238;
t275 = (qJD(1) * t280 + t86) * qJD(2) - t257;
t272 = pkin(6) - pkin(7);
t271 = g(3) * t97;
t270 = pkin(4) * t168;
t166 = g(1) * t180;
t269 = g(2) * t183;
t160 = t182 * pkin(3);
t268 = t83 * t81;
t263 = qJ(5) * t83;
t152 = pkin(6) * t238;
t104 = pkin(7) * t238 - t152;
t223 = t273 * qJD(2);
t68 = qJD(3) - t223 - t104;
t153 = pkin(6) * t237;
t106 = -pkin(7) * t237 + t153;
t172 = qJD(2) * qJ(3);
t85 = t106 + t172;
t35 = -t178 * t85 + t181 * t68;
t18 = t35 - t263;
t16 = -pkin(4) * t169 + t18;
t266 = -t18 + t16;
t51 = t181 * t104 + t178 * t106;
t265 = qJ(5) * t26;
t264 = qJ(5) * t81;
t36 = t178 * t68 + t181 * t85;
t19 = t36 - t264;
t262 = t169 * t19;
t261 = t169 * t35;
t260 = t169 * t36;
t116 = t272 * t179;
t117 = t272 * t182;
t56 = t178 * t116 + t181 * t117;
t175 = qJDD(1) * pkin(1);
t255 = qJDD(2) * pkin(2);
t146 = pkin(4) * t181 + pkin(3);
t254 = t146 * t182;
t252 = t179 * t180;
t251 = t179 * t181;
t250 = t179 * t183;
t186 = qJD(1) ^ 2;
t249 = t179 * t186;
t248 = t180 * t182;
t247 = t182 * t178;
t246 = t182 * t183;
t157 = t179 * qJD(3);
t244 = qJ(3) * t235 + t157;
t242 = t183 * pkin(1) + t180 * pkin(6);
t173 = t179 ^ 2;
t174 = t182 ^ 2;
t240 = -t173 + t174;
t239 = t173 + t174;
t236 = qJD(2) * t179;
t234 = qJD(4) * t178;
t70 = -qJD(3) * t178 - t109 * qJD(4);
t229 = -t109 * t27 - t69 * t81 - t70 * t83;
t228 = pkin(4) * t253;
t227 = t180 * t247;
t226 = t178 * t246;
t225 = t160 + t243;
t224 = -g(1) * t250 - g(2) * t252 + g(3) * t182;
t147 = pkin(6) * t156;
t219 = pkin(6) * t220 + qJDD(3) + t147;
t218 = -pkin(4) * t81 - qJD(5);
t217 = t166 - t269;
t48 = pkin(7) * t282 - t273 * qJDD(2) + t219;
t148 = pkin(6) * t230;
t170 = qJDD(2) * qJ(3);
t171 = qJD(2) * qJD(3);
t64 = -pkin(6) * t221 + t148 + t170 + t171;
t49 = pkin(7) * t281 + t64;
t5 = t178 * t48 + t181 * t49 + t68 * t233 - t85 * t234;
t6 = -t178 * t49 + t181 * t48 - t85 * t233 - t68 * t234;
t50 = -t104 * t178 + t181 * t106;
t216 = -qJD(2) * pkin(2) + qJD(3);
t55 = t181 * t116 - t117 * t178;
t93 = pkin(1) + t225;
t214 = t169 ^ 2;
t213 = pkin(2) * t246 + qJ(3) * t250 + t242;
t212 = -t147 - t224;
t63 = pkin(3) * t237 - t86;
t211 = t179 * t223;
t210 = t179 * t220;
t73 = -t180 * t251 + t227;
t75 = -t181 * t250 + t226;
t209 = -g(1) * t73 + g(2) * t75;
t74 = t97 * t180;
t76 = t97 * t183;
t208 = g(1) * t74 - g(2) * t76;
t207 = t239 * qJDD(1) * pkin(6);
t185 = qJD(2) ^ 2;
t206 = pkin(6) * t185 + t269;
t204 = pkin(2) * t230 - qJ(3) * t282 + qJD(1) * t157 + t175;
t110 = t152 + t216;
t115 = t153 + t172;
t203 = t110 * t182 - t115 * t179;
t202 = t247 - t251;
t143 = qJ(3) * t237;
t72 = -t273 * t238 + t143;
t71 = t219 - t255;
t200 = -0.2e1 * pkin(1) * t231 - t257;
t105 = t272 * t236;
t107 = qJD(2) * t117;
t20 = -t181 * t105 + t178 * t107 + t116 * t233 - t117 * t234;
t61 = -t211 + t244;
t198 = -t206 + 0.2e1 * t175;
t196 = g(1) * t76 + g(2) * t74 - g(3) * t202 - t5;
t195 = g(1) * t75 + g(2) * t73 + t271 + t6;
t21 = -t56 * qJD(4) + t105 * t178 + t181 * t107;
t44 = pkin(2) * t221 - t204;
t78 = pkin(2) * t236 - t244;
t194 = -qJD(1) * t78 - qJDD(1) * t280 - t206 - t44;
t28 = pkin(3) * t230 - qJD(1) * t211 + t204;
t193 = t63 * t81 + t196;
t192 = -t63 * t83 + t195;
t191 = t203 * qJD(2) + t71 * t179 + t64 * t182;
t190 = t195 + t265;
t39 = -t218 + t63;
t189 = t39 * t81 + t196 - t278;
t11 = pkin(4) * t27 + qJDD(5) + t28;
t177 = -qJ(5) - pkin(7);
t162 = t183 * pkin(6);
t138 = g(1) * t248;
t134 = qJ(3) * t246;
t132 = qJ(3) * t248;
t128 = t182 * t249;
t114 = t240 * t186;
t113 = qJDD(2) * t182 - t179 * t185;
t112 = qJDD(2) * t179 + t182 * t185;
t103 = -pkin(4) + t108;
t102 = pkin(2) * t238 - t143;
t95 = qJDD(1) * t174 - 0.2e1 * t210;
t94 = qJDD(1) * t173 + 0.2e1 * t210;
t65 = t179 * t230 + t240 * t231;
t60 = t70 * t169;
t54 = pkin(4) * t97 + t93;
t53 = t97 * qJD(2) - t179 * t234 - t181 * t232;
t52 = -t181 * t236 + t283;
t45 = -pkin(4) * t83 + t72;
t38 = -qJ(5) * t97 + t56;
t37 = qJ(5) * t202 + t55;
t34 = t168 * t178 - t181 * t214 - t83 * t238;
t33 = -t168 * t181 - t178 * t214 - t81 * t238;
t30 = t263 + t51;
t29 = t50 - t264;
t25 = t168 * t202 - t169 * t53;
t24 = t168 * t97 + t169 * t52;
t22 = pkin(4) * t52 + t61;
t10 = -qJ(5) * t53 + qJD(5) * t202 + t21;
t9 = -qJ(5) * t52 - qJD(5) * t97 + t20;
t8 = t27 * t97 + t52 * t81;
t7 = t202 * t26 + t53 * t83;
t4 = t5 + t278;
t3 = t12 * t181 - t14 * t178;
t2 = -qJD(5) * t83 + t265 - t270 + t6;
t1 = t202 * t27 + t26 * t97 - t52 * t83 - t53 * t81;
t13 = [0, 0, 0, 0, 0, qJDD(1), t217, t276, 0, 0, t94, 0.2e1 * t65, t112, t95, t113, 0, t179 * t200 + t182 * t198 + t138, t200 * t182 + (-t198 - t166) * t179, 0.2e1 * t207 - t276, -g(1) * (-pkin(1) * t180 + t162) - g(2) * t242 + (t239 * pkin(6) ^ 2 + pkin(1) ^ 2) * qJDD(1), t94, t112, -0.2e1 * t65, 0, -t113, t95, t275 * t179 + t194 * t182 + t138, t207 + t191 - t276, -t275 * t182 + (t194 + t166) * t179, pkin(6) * t191 - g(1) * t162 - g(2) * t213 + t86 * t78 + (-t166 + t44) * t280, t7, t1, t25, t8, t24, 0, -t168 * t55 - t169 * t21 + t27 * t93 + t28 * t97 + t52 * t63 + t61 * t81 + t208, t168 * t56 + t169 * t20 - t202 * t28 - t26 * t93 + t53 * t63 + t61 * t83 + t209, -t20 * t81 + t202 * t6 - t21 * t83 + t26 * t55 - t27 * t56 - t35 * t53 - t36 * t52 - t5 * t97 + t276, t5 * t56 + t36 * t20 + t6 * t55 + t35 * t21 + t28 * t93 + t63 * t61 - g(1) * (-pkin(7) * t183 + t162) - g(2) * (pkin(3) * t246 + t213) + (-g(1) * (t280 - t160) + g(2) * pkin(7)) * t180, t7, t1, t25, t8, t24, 0, -t10 * t169 + t11 * t97 - t168 * t37 + t22 * t81 + t27 * t54 + t39 * t52 + t208, -t11 * t202 + t168 * t38 + t169 * t9 + t22 * t83 - t26 * t54 + t39 * t53 + t209, -t10 * t83 - t16 * t53 - t19 * t52 + t2 * t202 + t26 * t37 - t27 * t38 - t4 * t97 - t81 * t9 + t276, t4 * t38 + t19 * t9 + t2 * t37 + t16 * t10 + t11 * t54 + t39 * t22 - g(1) * (t177 * t183 + t162) - g(2) * (t146 * t246 + t183 * t228 + t213) + (-g(1) * (t280 - t228 - t254) - g(2) * t177) * t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t128, -t114, t156, t128, t230, qJDD(2), pkin(1) * t249 + t212, g(3) * t179 - t148 + (pkin(1) * t186 + t276) * t182, 0, 0, -t128, t156, t114, qJDD(2), -t230, t128, 0.2e1 * t255 - qJDD(3) + (t102 * t182 - t179 * t86) * qJD(1) + t212, (-pkin(2) * t179 + qJ(3) * t182) * qJDD(1) + ((t115 - t172) * t179 + (-t110 + t216) * t182) * qJD(1), t148 + 0.2e1 * t170 + 0.2e1 * t171 + (qJD(1) * t102 - g(3)) * t179 + (qJD(1) * t86 - t276) * t182, t64 * qJ(3) + t115 * qJD(3) - t71 * pkin(2) - t86 * t102 - g(1) * (-pkin(2) * t250 + t134) - g(2) * (-pkin(2) * t252 + t132) - g(3) * t243 - t203 * qJD(1) * pkin(6), -t268, -t32, t12, t268, t14, t168, -t108 * t168 + t169 * t50 - t72 * t81 - t192 - t60, -t169 * t51 - t72 * t83 - t193 + t279, t108 * t26 + (-t36 + t50) * t83 + (t35 + t51) * t81 + t229, t5 * t109 + t6 * t108 - t63 * t72 - g(1) * t134 - g(2) * t132 - g(3) * t225 + (t69 - t51) * t36 + (t70 - t50) * t35 + t276 * t179 * t273, -t268, -t32, t12, t268, t14, t168, t169 * t29 - t45 * t81 - t60 + (qJD(5) + t39) * t83 + (pkin(4) - t103) * t168 - t190, -t169 * t30 - t45 * t83 - t189 + t279, t103 * t26 + (-t19 + t29) * t83 + (t16 + t30) * t81 + t229, t4 * t109 + t2 * t103 - t39 * t45 - g(1) * (pkin(4) * t226 + t134) - g(2) * (pkin(4) * t227 + t132) - g(3) * (t243 + t254) + (t69 - t30) * t19 + (-g(3) * pkin(4) * t178 + t276 * (pkin(2) + t146)) * t179 + (t70 - t29) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2) - t128, t156, -t173 * t186 - t185, -qJD(2) * t115 + t238 * t86 + t224 + t71, 0, 0, 0, 0, 0, 0, t33, t34, t3, -t63 * t238 + (t6 - t260) * t181 + (t5 + t261) * t178 + t224, 0, 0, 0, 0, 0, 0, t33, t34, t3, -t39 * t238 + (t2 - t262) * t181 + (t16 * t169 + t4) * t178 + t224; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t268, t32, -t12, -t268, -t14, -t168, t192 - t260, t193 - t261, 0, 0, t268, t32, -t12, -t268, -t14, -t168, -0.2e1 * t270 - t262 + (t218 - t39) * t83 + t190, -pkin(4) * t274 - t169 * t18 + t189, pkin(4) * t26 - t266 * t81, t266 * t19 + (t276 * t202 - t39 * t83 + t2 + t271) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 - t258, -t26 + t259, -t79 - t274, t16 * t83 + t19 * t81 + t11 + t217;];
tau_reg = t13;
