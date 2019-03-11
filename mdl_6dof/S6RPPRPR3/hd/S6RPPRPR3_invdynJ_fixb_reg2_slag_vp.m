% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPRPR3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:45:05
% EndTime: 2019-03-09 01:45:08
% DurationCPUTime: 2.74s
% Computational Cost: add. (5087->411), mult. (9423->499), div. (0->0), fcn. (6151->14), ass. (0->216)
t282 = 2 * qJD(4);
t142 = sin(pkin(10));
t147 = sin(qJ(4));
t150 = cos(qJ(4));
t250 = cos(pkin(10));
t167 = -t142 * t150 - t147 * t250;
t82 = t167 * qJD(1);
t275 = qJD(6) - t82;
t146 = sin(qJ(6));
t149 = cos(qJ(6));
t227 = t150 * qJD(2);
t144 = cos(pkin(9));
t119 = -t144 * pkin(1) - pkin(2);
t108 = -pkin(7) + t119;
t90 = qJD(1) * t108 + qJD(3);
t270 = (qJ(5) * qJD(1) - t90) * t147 - t227;
t56 = t250 * t270;
t233 = qJD(1) * t150;
t68 = -t147 * qJD(2) + t150 * t90;
t60 = -qJ(5) * t233 + t68;
t58 = (qJD(4) * pkin(4)) + t60;
t29 = t142 * t58 - t56;
t24 = (qJD(4) * pkin(8)) + t29;
t234 = qJD(1) * t147;
t143 = sin(pkin(9));
t112 = t143 * pkin(1) + qJ(3);
t97 = qJD(1) * t112;
t80 = pkin(4) * t234 + qJD(5) + t97;
t209 = t250 * t150;
t192 = qJD(1) * t209;
t212 = t142 * t234;
t85 = t192 - t212;
t35 = -t82 * pkin(5) - t85 * pkin(8) + t80;
t9 = -t146 * t24 + t149 * t35;
t281 = t9 * t275;
t10 = t146 * t35 + t149 * t24;
t280 = t10 * t275;
t136 = qJ(1) + pkin(9);
t126 = sin(t136);
t128 = cos(t136);
t274 = g(1) * t126 - g(2) * t128;
t180 = -t97 * qJD(1) - t274;
t206 = t149 * t275;
t202 = qJDD(1) * t250;
t220 = t150 * qJDD(1);
t218 = qJD(4) * t192 + t142 * t220 + t147 * t202;
t51 = qJD(4) * t212 - t218;
t49 = -qJDD(6) + t51;
t252 = t146 * t49;
t279 = t206 * t275 - t252;
t228 = t149 * qJD(4);
t230 = qJD(6) * t146;
t221 = t147 * qJDD(1);
t177 = -t142 * t221 + t150 * t202;
t87 = t167 * qJD(4);
t52 = -qJD(1) * t87 - t177;
t26 = -qJD(6) * t228 - t146 * qJDD(4) + t149 * t52 + t230 * t85;
t66 = t146 * qJD(4) + t149 * t85;
t232 = qJD(4) * t147;
t84 = -qJD(4) * t209 + t142 * t232;
t257 = t167 * t26 - t66 * t84;
t207 = t146 * t275;
t277 = t66 * t207;
t249 = pkin(1) * qJDD(1);
t135 = qJ(4) + pkin(10);
t127 = cos(t135);
t190 = g(1) * t128 + g(2) * t126;
t276 = t190 * t127;
t273 = qJDD(1) * t119;
t117 = t142 * pkin(4) + pkin(8);
t125 = sin(t135);
t161 = g(3) * t125 - t127 * t274;
t170 = -qJ(5) * qJDD(1) - qJD(1) * qJD(5);
t89 = qJDD(1) * t108 + qJDD(3);
t77 = t150 * t89;
t200 = -t147 * qJDD(2) + t77;
t25 = qJDD(4) * pkin(4) + t270 * qJD(4) + t170 * t150 + t200;
t226 = qJD(1) * qJD(4);
t211 = t150 * t226;
t231 = qJD(4) * t150;
t219 = t150 * qJDD(2) + t147 * t89 + t90 * t231;
t225 = qJD(2) * qJD(4);
t30 = -qJ(5) * t211 + (t170 - t225) * t147 + t219;
t7 = -t142 * t30 + t250 * t25;
t5 = -qJDD(4) * pkin(5) - t7;
t272 = -qJD(6) * t117 * t275 + t161 - t5;
t265 = g(3) * t127;
t271 = -t125 * t274 - t265;
t93 = -t142 * t147 + t209;
t184 = -t275 * t87 + t49 * t93;
t216 = t93 * t230;
t269 = -t149 * t184 - t216 * t275;
t37 = -t147 * t225 + t219;
t69 = t147 * t90 + t227;
t243 = t69 * qJD(4);
t38 = t200 - t243;
t156 = -(t147 * t68 - t150 * t69) * qJD(4) + t37 * t147 + t38 * t150;
t268 = qJDD(4) * t108 + t97 * t282;
t267 = t85 ^ 2;
t133 = t147 * pkin(4);
t64 = t146 * t85 - t228;
t264 = t64 * t82;
t263 = t66 * t64;
t262 = t66 * t85;
t261 = t85 * t64;
t260 = t85 * t82;
t259 = t87 * t64;
t258 = t87 * t66;
t229 = qJD(6) * t149;
t208 = -t149 * qJDD(4) - t146 * t52;
t27 = qJD(6) * t66 + t208;
t256 = -t146 * t27 - t64 * t229;
t8 = t142 * t25 + t250 * t30;
t255 = t93 * t51 + t87 * t82;
t251 = t27 * t149;
t254 = -t149 * t259 - t93 * t251;
t253 = t142 * t270;
t44 = t149 * t49;
t248 = qJD(6) * t64;
t247 = t126 * t146;
t246 = t126 * t149;
t245 = t128 * t146;
t244 = t128 * t149;
t242 = t80 * qJD(1);
t240 = qJ(5) - t108;
t141 = qJDD(2) - g(3);
t120 = t143 * t249;
t237 = qJDD(1) * qJ(3) + t120;
t139 = t147 ^ 2;
t140 = t150 ^ 2;
t236 = t139 - t140;
t152 = qJD(4) ^ 2;
t153 = qJD(1) ^ 2;
t235 = -t152 - t153;
t105 = pkin(4) * t231 + qJD(3);
t138 = qJD(3) * qJD(1);
t223 = qJDD(4) * t147;
t222 = qJDD(4) * t150;
t215 = t93 * t229;
t214 = t150 * t153 * t147;
t151 = cos(qJ(1));
t213 = t151 * pkin(1) + t128 * pkin(2) + t126 * qJ(3);
t91 = t138 + t237;
t148 = sin(qJ(1));
t210 = -t148 * pkin(1) + t128 * qJ(3);
t205 = t240 * t150;
t201 = -qJD(6) * t167 + qJD(1);
t199 = (-t139 - t140) * qJDD(1);
t197 = t66 * t215;
t196 = t147 * t211;
t62 = qJDD(5) + t91 + (t211 + t221) * pkin(4);
t14 = -t51 * pkin(5) + t52 * pkin(8) + t62;
t13 = t149 * t14;
t6 = qJDD(4) * pkin(8) + t8;
t2 = -qJD(6) * t10 - t146 * t6 + t13;
t195 = t167 * t2 + t9 * t84;
t94 = t112 + t133;
t1 = qJD(6) * t9 + t146 * t14 + t149 * t6;
t194 = -t1 * t167 - t10 * t84;
t28 = t250 * t58 + t253;
t23 = -qJD(4) * pkin(5) - t28;
t193 = t23 * t87 + t5 * t93;
t191 = t125 * pkin(5) - t127 * pkin(8);
t188 = g(1) * t148 - g(2) * t151;
t187 = -t93 * t26 + t258;
t186 = t167 * t27 + t84 * t64;
t185 = t93 * t27 + t259;
t183 = t167 * t51 + t82 * t84;
t182 = t167 * t52 - t85 * t84;
t181 = -t93 * t52 + t85 * t87;
t179 = -t10 * t149 + t146 * t9;
t178 = t10 * t146 + t149 * t9;
t43 = -pkin(5) * t167 - t93 * pkin(8) + t94;
t88 = t240 * t147;
t46 = -t142 * t205 - t250 * t88;
t16 = -t146 * t46 + t149 * t43;
t17 = t146 * t43 + t149 * t46;
t174 = t97 * qJD(3) + t91 * t112;
t53 = t84 * qJD(4) + qJDD(4) * t167;
t54 = t87 * qJD(4) + t93 * qJDD(4);
t173 = -t44 + (t146 * t82 - t230) * t275;
t172 = -t126 * pkin(2) + t210;
t171 = -qJD(6) * t35 + t265 - t6;
t169 = t117 * t49 + t23 * t275;
t145 = -qJ(5) - pkin(7);
t168 = t126 * t133 - t128 * t145 + t213;
t166 = qJDD(3) + t273;
t164 = t112 * qJDD(1) - t190;
t163 = t126 * t145 + t128 * t133 + t172;
t162 = -t150 * qJD(5) + t232 * t240;
t160 = t225 - t180;
t159 = -t167 * t8 + t28 * t87 - t29 * t84 + t7 * t93 - t274;
t158 = t146 * t184 - t215 * t275;
t157 = -qJD(6) * t178 + t1 * t149 - t2 * t146;
t155 = -t108 * t152 + t138 + t164 + t91;
t118 = -pkin(4) * t250 - pkin(5);
t96 = -t152 * t147 + t222;
t95 = -t152 * t150 - t223;
t81 = t82 ^ 2;
t73 = t125 * t244 - t247;
t72 = t125 * t245 + t246;
t71 = t125 * t246 + t245;
t70 = -t125 * t247 + t244;
t67 = -qJD(4) * t205 - t147 * qJD(5);
t45 = -t142 * t88 + t205 * t250;
t41 = pkin(4) * t233 + t85 * pkin(5) - t82 * pkin(8);
t40 = -t84 * pkin(5) - t87 * pkin(8) + t105;
t34 = t142 * t162 + t250 * t67;
t33 = t142 * t67 - t162 * t250;
t32 = t250 * t60 + t253;
t31 = t142 * t60 - t56;
t12 = t146 * t41 + t149 * t32;
t11 = -t146 * t32 + t149 * t41;
t4 = -qJD(6) * t17 - t146 * t34 + t149 * t40;
t3 = qJD(6) * t16 + t146 * t40 + t149 * t34;
t15 = [0, 0, 0, 0, 0, qJDD(1), t188, g(1) * t151 + g(2) * t148, 0, 0, 0, 0, 0, 0, 0, qJDD(1), 0.2e1 * t144 * t249 + t274, -0.2e1 * t120 + t190, 0 (t188 + (t143 ^ 2 + t144 ^ 2) * t249) * pkin(1), qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(3) - t274 + 0.2e1 * t273, 0.2e1 * t138 + t164 + t237, -g(1) * t172 - g(2) * t213 + t119 * t166 + t174, t140 * qJDD(1) - 0.2e1 * t196, -0.2e1 * t147 * t220 + 0.2e1 * t226 * t236, t96, t139 * qJDD(1) + 0.2e1 * t196, t95, 0, t155 * t147 + t268 * t150, -t268 * t147 + t155 * t150, t108 * t199 - t156 + t274, -g(1) * ((-pkin(2) - pkin(7)) * t126 + t210) - g(2) * (t128 * pkin(7) + t213) + t156 * t108 + t174, t181, -t182 + t255, t54, t183, t53, 0, -t33 * qJD(4) - t45 * qJDD(4) - t105 * t82 - t125 * t190 - t167 * t62 - t94 * t51 - t80 * t84, -t34 * qJD(4) - t46 * qJDD(4) + t105 * t85 - t94 * t52 + t62 * t93 + t80 * t87 - t276, t33 * t85 + t34 * t82 - t45 * t52 + t46 * t51 - t159, -g(1) * t163 - g(2) * t168 + t80 * t105 - t28 * t33 + t29 * t34 - t7 * t45 + t8 * t46 + t62 * t94, t149 * t187 - t216 * t66, -t197 + (-t258 + (t26 + t248) * t93) * t146 + t254, t257 + t269, t146 * t185 + t215 * t64, t158 + t186, t167 * t49 - t275 * t84, -g(1) * t73 - g(2) * t71 + t146 * t193 - t16 * t49 + t215 * t23 + t45 * t27 + t275 * t4 + t33 * t64 - t195, g(1) * t72 - g(2) * t70 + t149 * t193 + t17 * t49 - t216 * t23 - t45 * t26 - t275 * t3 + t33 * t66 - t194, t16 * t26 - t17 * t27 - t3 * t64 - t4 * t66 - t178 * t87 + t276 + (qJD(6) * t179 - t1 * t146 - t2 * t149) * t93, t1 * t17 + t10 * t3 + t2 * t16 + t9 * t4 + t5 * t45 + t23 * t33 - g(1) * (t128 * t191 + t163) - g(2) * (t126 * t191 + t168); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, 0, 0, 0, 0, 0, 0, 0, 0, 0, t141, 0, 0, 0, 0, 0, 0, t95, -t96, 0, -t38 * t147 + t37 * t150 - g(3) + (-t147 * t69 - t150 * t68) * qJD(4), 0, 0, 0, 0, 0, 0, t53, -t54, t182 + t255, t167 * t7 + t28 * t84 + t29 * t87 + t8 * t93 - g(3), 0, 0, 0, 0, 0, 0, t158 - t186, t257 - t269, t197 + (t258 + (-t26 + t248) * t93) * t146 + t254, t157 * t93 - t167 * t5 - t179 * t87 - t23 * t84 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t153, t166 + t180, 0, 0, 0, 0, 0, 0, t147 * t235 + t222, t150 * t235 - t223, t199, t156 + t180, 0, 0, 0, 0, 0, 0, qJD(1) * t82 + t54, -qJD(1) * t85 + t53, -t181 - t183, t159 - t242, 0, 0, 0, 0, 0, 0, -t167 * t252 + (t146 * t84 - t149 * t201) * t275 - t185, -t167 * t44 + (t146 * t201 + t149 * t84) * t275 - t187 (t201 * t66 + t186) * t149 + (t201 * t64 + t257) * t146 (-t201 * t9 + t194) * t149 + (-t10 * t201 + t195) * t146 - t193 - t274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t214, -t236 * t153, t220, -t214, -t221, qJDD(4), t243 + t77 + (-qJD(4) * t90 - t141) * t147 - t160 * t150, g(3) * t150 + t68 * qJD(4) + t147 * t160 - t219, 0, 0, -t260, -t81 + t267, t177, t260 (t85 + t212) * qJD(4) - t218, qJDD(4), t31 * qJD(4) - t80 * t85 + (qJDD(4) * t250 + t233 * t82) * pkin(4) + t161 + t7, t32 * qJD(4) - t80 * t82 + (-qJDD(4) * t142 - t233 * t85) * pkin(4) - t8 - t271 (t29 - t31) * t85 - (-t28 + t32) * t82 + (t142 * t51 + t250 * t52) * pkin(4), t28 * t31 - t29 * t32 + (t250 * t7 + g(3) * t147 + t142 * t8 + (-t274 - t242) * t150) * pkin(4), -t26 * t146 + t206 * t66 (-t26 + t264) * t149 - t277 + t256, -t262 + t279, t207 * t64 - t251, t173 + t261, -t275 * t85, -t11 * t275 + t118 * t27 + t169 * t146 + t272 * t149 - t31 * t64 - t9 * t85, t10 * t85 - t118 * t26 + t12 * t275 - t272 * t146 + t169 * t149 - t31 * t66, t11 * t66 + t12 * t64 + (-t117 * t27 + t82 * t9 + t1 + (t117 * t66 - t9) * qJD(6)) * t149 + (t10 * t82 - t117 * t26 - t2 + (t117 * t64 - t10) * qJD(6)) * t146 + t271, t5 * t118 - t10 * t12 - t9 * t11 - t23 * t31 - g(3) * (-t133 - t191) + t157 * t117 - t274 * (pkin(4) * t150 + pkin(5) * t127 + pkin(8) * t125); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t85 - t212) * qJD(4) + t218, t82 * t282 + t177, -t81 - t267, t28 * t85 - t29 * t82 - t190 + t62, 0, 0, 0, 0, 0, 0, t173 - t261, -t262 - t279 (t26 + t264) * t149 + t277 + t256, -t23 * t85 + (t2 + t280) * t149 + (t1 - t281) * t146 - t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t263, -t64 ^ 2 + t66 ^ 2, t275 * t64 - t26, -t263, -t208 + (-qJD(6) + t275) * t66, -t49, -g(1) * t70 - g(2) * t72 + t146 * t171 - t229 * t24 - t23 * t66 + t13 + t280, g(1) * t71 - g(2) * t73 + t23 * t64 + t281 + (qJD(6) * t24 - t14) * t146 + t171 * t149, 0, 0;];
tau_reg  = t15;
