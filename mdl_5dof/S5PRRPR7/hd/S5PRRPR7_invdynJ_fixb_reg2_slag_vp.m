% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR7_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR7_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR7_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:37:41
% EndTime: 2019-12-05 16:37:55
% DurationCPUTime: 5.86s
% Computational Cost: add. (4018->515), mult. (9739->748), div. (0->0), fcn. (7647->12), ass. (0->231)
t173 = sin(qJ(3));
t176 = cos(qJ(3));
t201 = pkin(3) * t173 - qJ(4) * t176;
t111 = qJD(3) * t201 - t173 * qJD(4);
t168 = sin(pkin(10));
t169 = sin(pkin(5));
t170 = cos(pkin(10));
t174 = sin(qJ(2));
t177 = cos(qJ(2));
t266 = t176 * t177;
t93 = (t168 * t174 + t170 * t266) * t169;
t316 = -qJD(1) * t93 + t168 * t111;
t245 = t176 * qJDD(2);
t249 = qJD(2) * qJD(3);
t186 = t173 * t249 - t245;
t255 = qJD(3) * t173;
t315 = (-pkin(7) * t170 + pkin(8)) * t255 + t316;
t205 = pkin(4) * t168 - pkin(8) * t170;
t197 = pkin(7) + t205;
t261 = qJD(1) * t169;
t232 = t177 * t261;
t208 = t173 * t232;
t254 = qJD(3) * t176;
t314 = t197 * t254 - t208;
t250 = qJD(1) * qJD(2);
t228 = t177 * t250;
t171 = cos(pkin(5));
t260 = qJD(1) * t171;
t280 = qJDD(2) * pkin(7);
t313 = qJD(3) * t260 + t280 + (qJDD(1) * t174 + t228) * t169;
t259 = qJD(2) * t173;
t130 = t201 * qJD(2);
t236 = t174 * t261;
t131 = qJD(2) * pkin(7) + t236;
t123 = t173 * t131;
t85 = t176 * t260 - t123;
t42 = t130 * t168 + t170 * t85;
t312 = -pkin(8) * t259 + qJD(4) * t170 - t42;
t311 = -pkin(3) * t176 - qJ(4) * t173;
t248 = qJDD(1) * t171;
t211 = t131 * t254 + t173 * t313 - t176 * t248;
t305 = -qJDD(3) * pkin(3) + qJDD(4);
t23 = t211 + t305;
t163 = t170 * qJDD(3);
t225 = t176 * t249;
t246 = t173 * qJDD(2);
t187 = t225 + t246;
t88 = t168 * t187 - t163;
t247 = t168 * qJDD(3);
t89 = t170 * t187 + t247;
t14 = pkin(4) * t88 - pkin(8) * t89 + t23;
t172 = sin(qJ(5));
t175 = cos(qJ(5));
t257 = qJD(2) * t176;
t230 = t173 * t260;
t86 = t131 * t176 + t230;
t72 = qJD(3) * qJ(4) + t86;
t136 = -pkin(2) + t311;
t87 = qJD(2) * t136 - t232;
t28 = t168 * t87 + t170 * t72;
t22 = -pkin(8) * t257 + t28;
t164 = t170 * qJD(3);
t126 = t168 * t259 - t164;
t256 = qJD(3) * t168;
t128 = t170 * t259 + t256;
t69 = -qJD(3) * pkin(3) + qJD(4) - t85;
t26 = pkin(4) * t126 - pkin(8) * t128 + t69;
t199 = t172 * t22 - t175 * t26;
t241 = t173 * t248 + t176 * t313;
t20 = qJDD(3) * qJ(4) + (qJD(4) - t123) * qJD(3) + t241;
t229 = t174 * t250;
t275 = t169 * t177;
t200 = -qJDD(1) * t275 + t169 * t229;
t32 = qJD(2) * t111 + qJDD(2) * t136 + t200;
t10 = t168 * t32 + t170 * t20;
t4 = pkin(8) * t186 + t10;
t1 = -t199 * qJD(5) + t172 * t14 + t175 * t4;
t122 = qJD(5) + t126;
t310 = t122 * t199 + t1;
t8 = t172 * t26 + t175 * t22;
t2 = -qJD(5) * t8 + t14 * t175 - t172 * t4;
t309 = t8 * t122 + t2;
t76 = t128 * t172 + t175 * t257;
t17 = qJD(5) * t76 - t172 * t186 - t175 * t89;
t308 = -t122 * t76 + t17;
t231 = t172 * t257;
t18 = -qJD(5) * t231 + t172 * t89 + (qJD(5) * t128 - t186) * t175;
t78 = t128 * t175 - t231;
t307 = t122 * t78 - t18;
t239 = t169 * t266;
t210 = t168 * t239;
t306 = qJD(1) * t210 + (t111 - t236) * t170;
t166 = t173 ^ 2;
t167 = t176 ^ 2;
t263 = t166 - t167;
t213 = qJD(2) * t263;
t178 = qJD(3) ^ 2;
t283 = sin(pkin(9));
t215 = t283 * t174;
t284 = cos(pkin(9));
t216 = t284 * t177;
t112 = -t171 * t216 + t215;
t214 = t283 * t177;
t217 = t284 * t174;
t114 = t171 * t214 + t217;
t204 = g(1) * t114 + g(2) * t112;
t281 = qJDD(2) * pkin(2);
t96 = t200 - t281;
t304 = -pkin(7) * t178 + t169 * (-g(3) * t177 + t229) + t204 + t281 - t96;
t302 = pkin(4) * t170;
t300 = t78 * t76;
t105 = t197 * t173;
t270 = t170 * t176;
t95 = pkin(7) * t270 + t136 * t168;
t75 = -pkin(8) * t176 + t95;
t35 = t105 * t175 - t172 * t75;
t298 = qJD(5) * t35 + t172 * t314 + t175 * t315;
t36 = t105 * t172 + t175 * t75;
t297 = -qJD(5) * t36 - t172 * t315 + t175 * t314;
t237 = pkin(7) * t168 + pkin(4);
t296 = -t237 * t255 - t306;
t49 = t230 + (qJD(2) * t205 + t131) * t176;
t133 = -pkin(8) * t168 - pkin(3) - t302;
t272 = t170 * t172;
t90 = -qJ(4) * t272 + t133 * t175;
t295 = qJD(5) * t90 - t172 * t49 + t175 * t312;
t271 = t170 * t175;
t91 = qJ(4) * t271 + t133 * t172;
t294 = -qJD(5) * t91 - t172 * t312 - t175 * t49;
t244 = pkin(7) * t255;
t293 = t168 * t244 + t306;
t292 = -t170 * t244 + t316;
t291 = qJD(2) * pkin(2);
t288 = t168 * t89;
t287 = t170 * t88;
t81 = qJDD(5) + t88;
t286 = t172 * t81;
t285 = t175 * t81;
t179 = qJD(2) ^ 2;
t278 = t167 * t179;
t277 = t168 * t176;
t276 = t169 * t174;
t273 = t170 * t136;
t269 = t172 * t173;
t268 = t173 * t175;
t267 = t175 * t176;
t265 = qJDD(1) - g(3);
t264 = pkin(2) * t275 + pkin(7) * t276;
t262 = t166 + t167;
t258 = qJD(2) * t174;
t252 = qJD(5) * t172;
t251 = qJD(5) * t175;
t108 = t112 * pkin(2);
t243 = t112 * t311 - t108;
t109 = t114 * pkin(2);
t242 = t114 * t311 - t109;
t240 = t173 * t275;
t238 = t170 * t267;
t235 = t168 * t257;
t234 = t169 * t258;
t233 = qJD(2) * t275;
t224 = t177 * t249;
t223 = t168 * t246;
t222 = t168 * t245;
t221 = t170 * t246;
t220 = t170 * t245;
t219 = t169 * t284;
t218 = t169 * t283;
t212 = t122 ^ 2;
t209 = pkin(3) * t239 + qJ(4) * t240 + t264;
t207 = t173 * t225;
t113 = t171 * t217 + t214;
t115 = -t171 * t215 + t216;
t203 = g(1) * t115 + g(2) * t113;
t9 = -t168 * t20 + t170 * t32;
t27 = -t168 * t72 + t170 * t87;
t117 = -t171 * t176 + t173 * t276;
t118 = t171 * t173 + t176 * t276;
t62 = t118 * t170 - t168 * t275;
t33 = t117 * t175 - t172 * t62;
t34 = t117 * t172 + t175 * t62;
t41 = t130 * t170 - t168 * t85;
t198 = t126 * t170 + t128 * t168;
t196 = qJD(2) * (t126 + t164);
t195 = qJD(2) * (-t128 + t256);
t194 = qJDD(2) * t177 - t174 * t179;
t120 = t170 * t269 + t267;
t43 = -t112 * t277 - t113 * t170;
t45 = -t114 * t277 - t115 * t170;
t92 = -t170 * t276 + t210;
t190 = g(1) * t45 + g(2) * t43 + g(3) * t92;
t63 = t113 * t173 + t176 * t219;
t65 = t115 * t173 - t176 * t218;
t189 = g(1) * t65 + g(2) * t63 + g(3) * t117;
t64 = t113 * t176 - t173 * t219;
t66 = t115 * t176 + t173 * t218;
t188 = g(1) * t66 + g(2) * t64 + g(3) * t118;
t185 = t189 - t23;
t184 = -qJ(4) * t255 + (qJD(4) - t69) * t176;
t182 = t189 - t211;
t132 = -t232 - t291;
t181 = -pkin(7) * qJDD(3) + (t132 + t232 - t291) * qJD(3);
t24 = -t131 * t255 + t241;
t180 = t211 * t173 + t24 * t176 + (-t173 * t86 - t176 * t85) * qJD(3) - t203;
t155 = t173 * t179 * t176;
t124 = qJDD(2) * t167 - 0.2e1 * t207;
t121 = t170 * t268 - t172 * t176;
t110 = t117 * pkin(3);
t107 = (t238 + t269) * qJD(2);
t106 = t170 * t231 - t175 * t259;
t94 = -pkin(7) * t277 + t273;
t74 = t176 * t237 - t273;
t68 = qJD(3) * t118 + t173 * t233;
t67 = -qJD(3) * t117 + t176 * t233;
t61 = t118 * t168 + t170 * t275;
t57 = t65 * pkin(3);
t56 = t63 * pkin(3);
t51 = -t175 * t255 - t176 * t252 + (t172 * t254 + t173 * t251) * t170;
t50 = -qJD(3) * t238 + qJD(5) * t120 - t172 * t255;
t46 = -t114 * t270 + t115 * t168;
t44 = -t112 * t270 + t113 * t168;
t40 = t168 * t234 + t170 * t67;
t39 = t168 * t67 - t170 * t234;
t37 = -pkin(4) * t259 - t41;
t31 = t114 * t168 + t170 * t66;
t30 = t112 * t168 + t170 * t64;
t21 = pkin(4) * t257 - t27;
t6 = qJD(5) * t33 + t68 * t172 + t40 * t175;
t5 = -qJD(5) * t34 - t40 * t172 + t68 * t175;
t3 = -pkin(4) * t186 - t9;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t265, 0, 0, 0, 0, 0, 0, t194 * t169, (-qJDD(2) * t174 - t177 * t179) * t169, 0, -g(3) + (t171 ^ 2 + (t174 ^ 2 + t177 ^ 2) * t169 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, -t68 * qJD(3) - t117 * qJDD(3) + (-t173 * t224 + t176 * t194) * t169, -qJD(3) * t67 - qJDD(3) * t118 + (-t173 * t194 - t176 * t224) * t169, (t117 * t173 + t118 * t176) * qJDD(2) + (t173 * t68 + t176 * t67 + (t117 * t176 - t118 * t173) * qJD(3)) * qJD(2), t117 * t211 + t118 * t24 + t67 * t86 - t68 * t85 - g(3) + (t132 * t258 - t177 * t96) * t169, 0, 0, 0, 0, 0, 0, t61 * t245 + t117 * t88 + t68 * t126 + (t176 * t39 - t255 * t61) * qJD(2), t62 * t245 + t117 * t89 + t68 * t128 + (t176 * t40 - t255 * t62) * qJD(2), -t126 * t40 + t128 * t39 + t61 * t89 - t62 * t88, t10 * t62 + t117 * t23 - t27 * t39 + t28 * t40 - t61 * t9 + t68 * t69 - g(3), 0, 0, 0, 0, 0, 0, t122 * t5 + t18 * t61 + t33 * t81 + t39 * t76, -t122 * t6 - t17 * t61 - t34 * t81 + t39 * t78, t17 * t33 - t18 * t34 - t5 * t78 - t6 * t76, t1 * t34 - t199 * t5 + t2 * t33 + t21 * t39 + t3 * t61 + t6 * t8 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t265 * t275 + t204, -t265 * t276 + t203, 0, 0, qJDD(2) * t166 + 0.2e1 * t207, -0.2e1 * qJD(3) * t213 + 0.2e1 * t173 * t245, qJDD(3) * t173 + t176 * t178, t124, qJDD(3) * t176 - t173 * t178, 0, t173 * t181 + t176 * t304, -t173 * t304 + t176 * t181, t262 * t280 + (-g(3) * t174 - t228 * t262) * t169 + t180, -t96 * pkin(2) + g(1) * t109 + g(2) * t108 - g(3) * t264 + (-t132 * t174 + (t173 * t85 - t176 * t86) * t177) * t261 + t180 * pkin(7), (t128 * t254 + t173 * t89) * t170, (-t287 - t288) * t173 - t198 * t254, (-t89 - t221) * t176 + (t128 * t173 + t170 * t213) * qJD(3), (t126 * t254 + t173 * t88) * t168, (t88 + t223) * t176 + (-t126 * t173 - t168 * t213) * qJD(3), t124, -g(1) * t46 - g(2) * t44 - g(3) * t93 + (-t126 * t232 + pkin(7) * t88 + t168 * t23 + (qJD(2) * t94 + t27) * qJD(3)) * t173 + (-qJDD(2) * t94 - t9 + (pkin(7) * t126 + t168 * t69) * qJD(3) - t293 * qJD(2)) * t176, (-t128 * t232 + pkin(7) * t89 + t170 * t23 + (-qJD(2) * t95 - t28) * qJD(3)) * t173 + (qJDD(2) * t95 + t10 + (pkin(7) * t128 + t170 * t69) * qJD(3) + t292 * qJD(2)) * t176 + t190, -t88 * t95 - t89 * t94 - t293 * t128 - t292 * t126 + (-t168 * t28 - t170 * t27) * t254 + (-g(3) * t275 - t10 * t168 - t170 * t9 + t204) * t173, t10 * t95 + t9 * t94 - t69 * t208 - g(1) * t242 - g(2) * t243 - g(3) * t209 + t292 * t28 + t293 * t27 + (t173 * t23 + t254 * t69 - t203) * pkin(7), -t121 * t17 - t50 * t78, t120 * t17 - t121 * t18 + t50 * t76 - t51 * t78, t121 * t81 - t50 * t122 + (-t17 * t173 + t254 * t78) * t168, t120 * t18 + t51 * t76, -t120 * t81 - t51 * t122 + (-t173 * t18 - t254 * t76) * t168, (t122 * t254 + t173 * t81) * t168, t35 * t81 + t74 * t18 + t3 * t120 + t21 * t51 - g(1) * (-t114 * t269 + t175 * t46) - g(2) * (-t112 * t269 + t175 * t44) - g(3) * (t172 * t240 + t175 * t93) + t296 * t76 + (t173 * t2 - t199 * t254) * t168 + t297 * t122, -t36 * t81 - t74 * t17 + t3 * t121 - t21 * t50 - g(1) * (-t114 * t268 - t172 * t46) - g(2) * (-t112 * t268 - t172 * t44) - g(3) * (-t172 * t93 + t175 * t240) + t296 * t78 + (-t1 * t173 - t254 * t8) * t168 - t298 * t122, -t1 * t120 - t121 * t2 + t17 * t35 - t18 * t36 - t199 * t50 - t297 * t78 - t298 * t76 - t51 * t8 - t190, t1 * t36 + t2 * t35 + t3 * t74 - g(1) * (pkin(4) * t46 + pkin(7) * t115 + pkin(8) * t45 + t242) - g(2) * (pkin(4) * t44 + pkin(7) * t113 + pkin(8) * t43 + t243) - g(3) * (pkin(4) * t93 + pkin(8) * t92 + t209) + t298 * t8 - t297 * t199 + t296 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t155, t263 * t179, t246, t155, t245, qJDD(3), qJD(3) * t86 - t132 * t259 + t182, -t132 * t257 + (t85 + t123) * qJD(3) + t188 - t241, 0, 0, -t128 * t170 * t257 + t288, -t168 * t88 + t89 * t170 + t198 * t257, t170 * t278 + t173 * t195 - t222, -t126 * t235 - t287, -t168 * t278 + t173 * t196 - t220, t155, qJ(4) * t222 - pkin(3) * t88 - t126 * t86 + t185 * t170 + (t168 * t184 - t173 * t27 + t176 * t41) * qJD(2), qJ(4) * t220 - pkin(3) * t89 - t128 * t86 - t185 * t168 + (t170 * t184 + t173 * t28 - t176 * t42) * qJD(2), t126 * t42 + t128 * t41 + (-qJ(4) * t88 - qJD(4) * t126 + t257 * t27 + t10) * t170 + (qJ(4) * t89 + qJD(4) * t128 + t257 * t28 - t9) * t168 - t188, -t23 * pkin(3) + g(1) * t57 + g(2) * t56 + g(3) * t110 - t27 * t41 - t28 * t42 - t69 * t86 + (-t168 * t27 + t170 * t28) * qJD(4) + (t10 * t170 - t168 * t9 - t188) * qJ(4), -t78 * t107 + (-t17 * t175 - t252 * t78) * t168, t106 * t78 + t107 * t76 + (t17 * t172 - t175 * t18 + (t172 * t76 - t175 * t78) * qJD(5)) * t168, -t107 * t122 + t17 * t170 + (-t122 * t252 - t257 * t78 + t285) * t168, -t76 * t106 + (t172 * t18 + t251 * t76) * t168, t106 * t122 + t170 * t18 + (-t122 * t251 + t257 * t76 - t286) * t168, -t122 * t235 - t170 * t81, t90 * t81 - t2 * t170 - t37 * t76 - t21 * t106 - g(1) * (t172 * t66 - t271 * t65) - g(2) * (t172 * t64 - t271 * t63) - g(3) * (-t117 * t271 + t118 * t172) + t294 * t122 + (qJ(4) * t18 + qJD(4) * t76 + t172 * t3 + t199 * t257 + t21 * t251) * t168, -t91 * t81 + t1 * t170 - t37 * t78 - t21 * t107 - g(1) * (t175 * t66 + t272 * t65) - g(2) * (t175 * t64 + t272 * t63) - g(3) * (t117 * t272 + t118 * t175) - t295 * t122 + (-qJ(4) * t17 + qJD(4) * t78 + t175 * t3 - t21 * t252 + t257 * t8) * t168, t106 * t8 - t107 * t199 + t17 * t90 - t18 * t91 - t294 * t78 - t295 * t76 + (-t1 * t172 - t175 * t2 + (-t172 * t199 - t175 * t8) * qJD(5) + t189) * t168, t1 * t91 + t2 * t90 - t21 * t37 - g(1) * (qJ(4) * t66 - t302 * t65 - t57) - g(2) * (t64 * qJ(4) - t302 * t63 - t56) - g(3) * (t118 * qJ(4) - t117 * t302 - t110) + t295 * t8 - t294 * t199 + (pkin(8) * t189 + t3 * qJ(4) + t21 * qJD(4)) * t168; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176 * t195 - t163 + t223, t176 * t196 + t221 + t247, -t126 ^ 2 - t128 ^ 2, t126 * t28 + t128 * t27 - t182 + t305, 0, 0, 0, 0, 0, 0, -t128 * t76 - t172 * t212 + t285, -t128 * t78 - t175 * t212 - t286, t172 * t307 + t175 * t308, -t128 * t21 + t172 * t310 + t175 * t309 - t189; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t300, -t76 ^ 2 + t78 ^ 2, -t308, -t300, t307, t81, -t21 * t78 - g(1) * (-t172 * t31 + t175 * t65) - g(2) * (-t172 * t30 + t175 * t63) - g(3) * t33 + t309, t21 * t76 - g(1) * (-t172 * t65 - t175 * t31) - g(2) * (-t172 * t63 - t175 * t30) + g(3) * t34 - t310, 0, 0;];
tau_reg = t7;
