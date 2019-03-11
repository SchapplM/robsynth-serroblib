% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S6RPRPPR8
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% tau_reg [6x(6*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPRPPR8_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPPR8_invdynJ_fixb_reg2_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:00:06
% EndTime: 2019-03-09 03:00:12
% DurationCPUTime: 3.68s
% Computational Cost: add. (3057->444), mult. (5326->497), div. (0->0), fcn. (2698->6), ass. (0->245)
t135 = cos(qJ(3));
t234 = qJD(1) * t135;
t132 = sin(qJ(3));
t103 = t132 * qJDD(1);
t225 = qJD(1) * qJD(3);
t205 = t135 * t225;
t294 = -t103 - t205;
t137 = -pkin(3) - pkin(4);
t138 = -pkin(1) - pkin(7);
t74 = t138 * qJD(1) + qJD(2);
t40 = (qJ(5) * qJD(1) + t74) * t135;
t292 = qJD(4) - t40;
t28 = qJD(3) * t137 + t292;
t131 = sin(qJ(6));
t134 = cos(qJ(6));
t232 = qJD(3) * t134;
t236 = qJD(1) * t132;
t55 = t131 * t236 + t232;
t82 = qJD(6) + t234;
t173 = t55 * t82;
t256 = qJD(6) * t55;
t19 = t131 * qJDD(3) + t294 * t134 + t256;
t293 = t19 - t173;
t126 = -pkin(8) + t137;
t130 = qJ(4) + pkin(5);
t152 = t126 * t135 - t130 * t132;
t143 = qJD(3) * t152 - qJD(2);
t115 = t135 * pkin(5);
t171 = t126 * t132 - qJ(2);
t154 = t171 + t115;
t105 = t135 * qJDD(1);
t108 = t135 * qJD(4);
t272 = -qJ(4) * t105 - qJD(1) * t108;
t210 = qJDD(5) - t272;
t6 = qJD(1) * t143 + qJDD(1) * t154 + t210;
t100 = qJ(4) * t234;
t227 = qJD(5) + t100;
t22 = qJD(1) * t154 + t227;
t27 = qJD(3) * t126 + t292;
t7 = -t131 * t27 + t134 * t22;
t233 = qJD(3) * t132;
t58 = t74 * t233;
t284 = t138 * qJDD(1);
t73 = qJDD(2) + t284;
t62 = t135 * t73;
t211 = qJDD(4) + t58 - t62;
t224 = qJD(1) * qJD(5);
t226 = qJ(5) * qJDD(1);
t206 = t132 * t225;
t77 = qJ(5) * t206;
t142 = t77 + (-t224 - t226) * t135 + t211;
t9 = qJDD(3) * t126 + t142;
t1 = qJD(6) * t7 + t131 * t6 + t134 * t9;
t8 = t131 * t22 + t134 * t27;
t182 = t131 * t7 - t134 * t8;
t5 = t134 * t6;
t2 = -qJD(6) * t8 - t131 * t9 + t5;
t133 = sin(qJ(1));
t136 = cos(qJ(1));
t286 = g(1) * t136 + g(2) * t133;
t291 = t182 * qJD(6) - t1 * t131 - t134 * t2 - t286;
t290 = qJDD(1) * qJ(2);
t128 = t132 ^ 2;
t129 = t135 ^ 2;
t238 = t128 + t129;
t202 = t238 * t73;
t252 = t132 * t136;
t289 = pkin(4) * t252 + t133 * qJ(5);
t228 = t131 * qJD(3);
t235 = qJD(1) * t134;
t57 = t132 * t235 - t228;
t172 = t57 * t82;
t229 = qJD(6) * t134;
t208 = t132 * t229;
t20 = qJD(1) * (t135 * t228 + t208) - qJD(6) * t228 + t134 * qJDD(3) + t103 * t131;
t161 = -t20 + t172;
t111 = t135 * qJ(4);
t242 = t111 + t115;
t120 = g(2) * t136;
t121 = g(1) * t133;
t287 = t121 - t120;
t207 = qJDD(3) * t137;
t283 = t294 * qJ(5) - t132 * t224;
t123 = qJDD(3) * qJ(4);
t125 = qJD(3) * qJD(4);
t231 = qJD(3) * t135;
t59 = t74 * t231;
t61 = t132 * t73;
t23 = t123 + t59 + t61 + t125;
t15 = -t23 + t283;
t10 = qJDD(3) * pkin(5) - t15;
t276 = g(3) * t135;
t282 = qJD(6) * t126 * t82 + t132 * t287 - t10 + t276;
t261 = qJ(4) * t132;
t164 = t135 * t137 - t261;
t255 = qJDD(3) * pkin(3);
t24 = t211 - t255;
t197 = -t135 * t74 + qJD(4);
t271 = qJD(3) * pkin(3);
t42 = t197 - t271;
t127 = qJD(3) * qJ(4);
t63 = t132 * t74;
t44 = t63 + t127;
t141 = (t132 * t42 + t135 * t44) * qJD(3) + t23 * t132 - t24 * t135;
t183 = t131 * t8 + t134 * t7;
t280 = -qJD(6) * t183 + t1 * t134 - t2 * t131;
t279 = t7 * t82;
t278 = t8 * t82;
t277 = pkin(3) * t132;
t118 = g(3) * t132;
t275 = t57 * t55;
t251 = t133 * t135;
t253 = t132 * t133;
t274 = pkin(3) * t251 + qJ(4) * t253;
t249 = t135 * t136;
t273 = g(1) * t249 + g(2) * t251;
t270 = t131 * t19;
t54 = -qJDD(6) - t105 + t206;
t269 = t131 * t54;
t268 = t131 * t82;
t267 = t132 * t55;
t266 = t132 * t57;
t265 = t134 * t54;
t264 = t134 * t57;
t263 = t20 * t134;
t262 = pkin(1) * qJDD(1);
t99 = qJ(5) * t236;
t34 = -t99 - t44;
t29 = qJD(3) * pkin(5) - t34;
t260 = qJD(3) * t29;
t259 = qJD(3) * t34;
t258 = qJD(3) * t55;
t257 = qJD(3) * t57;
t140 = qJD(1) ^ 2;
t254 = t128 * t140;
t250 = t134 * t135;
t139 = qJD(3) ^ 2;
t248 = t138 * t139;
t247 = t140 * qJ(2);
t246 = qJ(5) + t138;
t214 = t137 * t132;
t190 = -qJ(2) + t214;
t33 = qJD(1) * t190 + t227;
t244 = -qJD(5) - t33;
t216 = 0.2e1 * qJD(1) * qJD(2);
t243 = (t216 + t290) * qJ(2);
t241 = t136 * pkin(1) + t133 * qJ(2);
t237 = t139 + t140;
t230 = qJD(6) * t131;
t223 = qJDD(3) * t132;
t222 = qJDD(3) * t135;
t221 = qJDD(3) * t138;
t220 = pkin(4) * t251 + t274;
t219 = g(1) * t251 - g(2) * t249 - t118;
t218 = t135 * t268;
t217 = t82 * t250;
t215 = t28 * t233 - t287;
t213 = t82 * t228;
t212 = t82 * t232;
t209 = t136 * pkin(7) + t241;
t204 = qJ(2) + t277;
t203 = t111 - t277;
t189 = t111 + t214;
t51 = -qJ(2) + t189;
t201 = qJD(1) * t51 + t33;
t14 = t207 + t142;
t200 = -t14 - t259;
t199 = -t62 + t219;
t196 = qJD(3) * t246;
t195 = 0.2e1 * t205;
t194 = pkin(3) * t253 + t209;
t193 = qJDD(2) - t262;
t192 = -g(2) * t252 + t276 - t61;
t191 = t132 * t205;
t112 = t136 * qJ(2);
t188 = t138 * t133 + t112;
t181 = -t287 - t247;
t179 = -qJDD(4) - t199;
t178 = pkin(3) * t135 + t261;
t31 = t171 + t242;
t65 = t246 * t135;
t16 = t131 * t65 + t134 * t31;
t17 = t131 * t31 - t134 * t65;
t43 = qJD(1) * t204 - t100;
t66 = qJ(2) - t203;
t175 = (qJD(1) * t66 + t43) * qJD(3);
t170 = t134 * t82;
t169 = -qJD(6) * t22 + t118 - t9;
t168 = t58 - t179;
t167 = t216 + 0.2e1 * t290;
t166 = -t229 * t82 + t269;
t165 = t230 * t82 + t265;
t144 = qJDD(1) * t190 + t210;
t145 = qJD(3) * t164 - qJD(2);
t11 = qJD(1) * t145 + t144;
t30 = t108 + t145;
t163 = qJD(1) * t30 + qJDD(1) * t51 + t11;
t162 = t77 + t168;
t158 = -pkin(8) * t132 + t242;
t157 = -t238 * t284 + t287;
t155 = -t133 * t111 + t194;
t151 = qJD(3) * t178 + qJD(2);
t150 = t126 * t54 - t29 * t82;
t149 = t167 - t248;
t18 = qJD(1) * t151 + qJDD(1) * t204 + t272;
t38 = -t108 + t151;
t148 = -qJD(1) * t38 - qJDD(1) * t66 - t18 + t248;
t88 = pkin(3) * t252;
t147 = -qJ(4) * t249 + t188 + t88;
t146 = -g(1) * t253 + 0.2e1 * t123 + 0.2e1 * t125 - t192;
t107 = t129 * t140;
t104 = t129 * qJDD(1);
t102 = t128 * qJDD(1);
t83 = pkin(4) * t253;
t80 = t135 * t140 * t132;
t79 = t135 * t221;
t75 = -t107 - t139;
t71 = qJDD(3) - t80;
t70 = -t107 + t254;
t69 = -t132 * t139 + t222;
t68 = t135 * t139 + t223;
t67 = -t102 - t104;
t64 = t246 * t132;
t60 = t178 * qJD(1);
t53 = t104 - 0.2e1 * t191;
t52 = t102 + 0.2e1 * t191;
t50 = -t132 * t237 + t222;
t49 = t135 * t237 + t223;
t48 = t131 * t133 - t134 * t249;
t47 = t131 * t249 + t133 * t134;
t46 = t131 * t136 + t133 * t250;
t45 = t131 * t251 - t134 * t136;
t41 = t164 * qJD(1);
t39 = t63 + t99;
t37 = t132 * qJD(5) + t135 * t196;
t36 = -t135 * qJD(5) + t132 * t196;
t35 = -t132 * t105 + (t128 - t129) * t225;
t32 = 0.2e1 * t35;
t26 = t152 * qJD(1);
t21 = t108 + t143;
t13 = t131 * t26 + t134 * t39;
t12 = -t131 * t39 + t134 * t26;
t4 = -qJD(6) * t17 - t131 * t36 + t134 * t21;
t3 = qJD(6) * t16 + t131 * t21 + t134 * t36;
t25 = [0, 0, 0, 0, 0, qJDD(1), t287, t286, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, qJDD(2) - t287 - 0.2e1 * t262, t167 - t286, -t193 * pkin(1) - g(1) * (-pkin(1) * t133 + t112) - g(2) * t241 + t243, t53, t32, t69, t52, -t68, 0, qJ(2) * t195 + t79 + (t149 - t286) * t132 (-0.2e1 * qJ(2) * t225 - t221) * t132 + t149 * t135 - t273, t157 - t202, -g(1) * t188 - g(2) * t209 + t138 * t202 + t243, t53, t69, -0.2e1 * t35, 0, t68, t52, t79 + t135 * t175 + (-t148 - t286) * t132, -t141 + t157 (t175 + t221) * t132 + t148 * t135 + t273, -g(1) * t147 - g(2) * t155 + t138 * t141 + t18 * t66 + t43 * t38, t52, t32, -t68, t53, t69, 0, qJDD(3) * t64 + t163 * t135 + (-t132 * t201 + t37) * qJD(3) + t273, -qJDD(3) * t65 + (t135 * t201 + t36) * qJD(3) + (t163 + t286) * t132 (qJDD(1) * t64 - t15 + (-qJD(3) * t65 + t37) * qJD(1)) * t132 + (qJDD(1) * t65 + (qJD(3) * t64 - t36) * qJD(1) + t200) * t135 + t215, -t14 * t65 + t28 * t36 - t15 * t64 - t34 * t37 + t11 * t51 + t33 * t30 - g(1) * (t147 + t289) - g(2) * (-qJ(5) * t136 + t155 + t83) t231 * t264 + (-t134 * t19 - t230 * t57) * t132 (-t131 * t57 - t134 * t55) * t231 + (t270 - t263 + (t131 * t55 - t264) * qJD(6)) * t132 (-t19 + t212) * t135 + (-t165 - t257) * t132, t55 * t208 + (t132 * t20 + t231 * t55) * t131 (-t20 - t213) * t135 + (t166 + t258) * t132, -t135 * t54 - t233 * t82, -g(1) * t48 + g(2) * t46 - t16 * t54 + t20 * t64 + t37 * t55 + t4 * t82 + (t228 * t29 + t2) * t135 + (-qJD(3) * t7 + t10 * t131 + t229 * t29) * t132, -g(1) * t47 - g(2) * t45 + t17 * t54 - t19 * t64 - t3 * t82 + t37 * t57 + (t232 * t29 - t1) * t135 + (qJD(3) * t8 + t10 * t134 - t230 * t29) * t132, t132 * t291 + t16 * t19 - t17 * t20 - t183 * t231 - t3 * t55 - t4 * t57, t1 * t17 + t8 * t3 + t2 * t16 + t7 * t4 + t10 * t64 + t29 * t37 - g(1) * (t112 + t88 + t289) - g(2) * (t83 + t194) + (g(1) * t158 + g(2) * qJ(5)) * t136 + (-g(1) * t138 + g(2) * t158) * t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t140, t181 + t193, 0, 0, 0, 0, 0, 0, t50, -t49, t67, t202 + t181, 0, 0, 0, 0, 0, 0, t50, t67, t49, -qJD(1) * t43 + t141 - t287, 0, 0, 0, 0, 0, 0, t49, -t50, -t67, qJD(1) * t33 - t132 * t15 + t135 * t200 + t215, 0, 0, 0, 0, 0, 0, t82 * t235 + (t20 - t213) * t132 + (-t166 + t258) * t135, -qJD(1) * t268 + (-t19 - t212) * t132 + (-t165 + t257) * t135 (-t55 * t233 - qJD(1) * t57 + (-qJD(6) * t57 + t20) * t135) * t134 + (t57 * t233 - qJD(1) * t55 + (t19 - t256) * t135) * t131, t183 * qJD(1) + (-qJD(3) * t182 + t10) * t132 + (t260 - t280) * t135 - t287; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, -t70, t105, -t80, -t103, qJDD(3), -t135 * t247 - t199 (t247 + t121) * t132 + t192, 0, 0, t80, t105, t70, qJDD(3), t103, -t80, 0.2e1 * t255 + (-t132 * t60 - t135 * t43) * qJD(1) + t179, -t178 * qJDD(1) + ((t44 - t127) * t135 + (-qJD(4) + t42 + t271) * t132) * qJD(1) (-t132 * t43 + t135 * t60) * qJD(1) + t146, -t24 * pkin(3) - g(1) * t274 - g(3) * t203 + t23 * qJ(4) + t178 * t120 + t197 * t44 - t42 * t63 - t43 * t60, -t80, -t70, -t103, t80, t105, qJDD(3), -qJD(3) * t40 + t59 + (t132 * t33 - t135 * t41) * qJD(1) + t146 - t283, -qJ(5) * t105 - qJD(3) * t39 + 0.2e1 * t207 + (-t132 * t41 + t135 * t244) * qJD(1) + t162, -t164 * qJDD(1) + (t34 + t39 + t127) * t234, -g(1) * t220 - g(3) * t189 - t15 * qJ(4) - t164 * t120 + t14 * t137 - t28 * t39 - t292 * t34 - t33 * t41, -t170 * t57 + t270 (t19 + t173) * t134 + (t20 + t172) * t131 (-t217 + t266) * qJD(1) + t166, -t131 * t173 + t263 (t218 - t267) * qJD(1) + t165, t82 * t236, -t12 * t82 + t130 * t20 + t150 * t131 - t282 * t134 + t7 * t236 + t292 * t55, t13 * t82 - t130 * t19 + t282 * t131 + t150 * t134 - t8 * t236 + t292 * t57, t12 * t57 + t13 * t55 + (t7 * t234 - t126 * t20 - t1 + (t126 * t57 + t7) * qJD(6)) * t134 + (t8 * t234 - t126 * t19 + t2 + (t126 * t55 + t8) * qJD(6)) * t131 - t219, t10 * t130 - t8 * t13 - t7 * t12 - g(1) * (pkin(8) * t251 + t220) - g(3) * t242 + t292 * t29 + (-pkin(5) * t121 - g(3) * t126) * t132 + t280 * t126 - t152 * t120; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, t105, t75, -qJD(3) * t44 + t234 * t43 + t168 - t255, 0, 0, 0, 0, 0, 0, t75, t71, -t105, t259 + t207 + (qJD(1) * t244 - t226) * t135 + t162, 0, 0, 0, 0, 0, 0, -t170 * t82 - t258 + t269, t131 * t82 ^ 2 - t257 + t265, -t131 * t293 + t161 * t134, -t260 + (t1 - t279) * t134 + (-t2 - t278) * t131 + t219; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105 - 0.2e1 * t206, t103 + t195, -t107 - t254 (t132 * t34 + t135 * t28 + t145) * qJD(1) + t144 + t286, 0, 0, 0, 0, 0, 0 (-t218 - t267) * qJD(1) - t165 (-t217 - t266) * qJD(1) + t166, t161 * t131 + t134 * t293 (-t132 * t29 - t135 * t182) * qJD(1) - t291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t275, -t55 ^ 2 + t57 ^ 2, -t293, -t275, t161, -t54, -g(1) * t45 + g(2) * t47 + t131 * t169 - t229 * t27 - t29 * t57 + t278 + t5, -g(1) * t46 - g(2) * t48 + t29 * t55 + t279 + (qJD(6) * t27 - t6) * t131 + t169 * t134, 0, 0;];
tau_reg  = t25;
