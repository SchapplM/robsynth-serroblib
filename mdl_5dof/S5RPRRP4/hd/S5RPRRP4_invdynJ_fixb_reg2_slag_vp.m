% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPRRP4
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:32:52
% EndTime: 2022-01-23 09:33:00
% DurationCPUTime: 4.02s
% Computational Cost: add. (4081->404), mult. (9943->498), div. (0->0), fcn. (7047->10), ass. (0->240)
t168 = sin(pkin(8));
t170 = sin(qJ(4));
t173 = cos(qJ(4));
t174 = cos(qJ(3));
t258 = t173 * t174;
t222 = t168 * t258;
t171 = sin(qJ(3));
t234 = qJDD(1) * t171;
t114 = t170 * t174 + t173 * t171;
t229 = qJD(3) + qJD(4);
t303 = t229 * t114;
t41 = -qJDD(1) * t222 + (qJD(1) * t303 + t170 * t234) * t168;
t228 = 2 * qJ(2);
t169 = cos(pkin(8));
t238 = t169 * qJD(1);
t142 = -qJD(3) + t238;
t272 = qJ(2) * t171;
t223 = t169 * t272;
t265 = t168 * t174;
t227 = pkin(7) * t265;
t186 = -t223 - t227;
t120 = pkin(2) * t169 + t168 * pkin(6) + pkin(1);
t100 = -t120 * qJD(1) + qJD(2);
t91 = t174 * t100;
t61 = qJD(1) * t186 + t91;
t47 = -t142 * pkin(3) + t61;
t245 = qJD(1) * t168;
t220 = t171 * t245;
t221 = qJ(2) * t238;
t75 = t171 * t100 + t174 * t221;
t62 = -pkin(7) * t220 + t75;
t54 = t170 * t62;
t30 = t173 * t47 - t54;
t198 = t170 * t220;
t244 = qJD(1) * t174;
t219 = t168 * t244;
t86 = t173 * t219 - t198;
t78 = t86 * qJ(5);
t18 = t30 - t78;
t271 = qJ(2) * t174;
t141 = t169 * t271;
t80 = -t171 * t120 + t141;
t302 = qJD(3) * t80;
t230 = t169 * qJDD(1);
t140 = -qJDD(3) + t230;
t127 = -qJDD(4) + t140;
t134 = -qJD(4) + t142;
t239 = qJD(4) * t173;
t226 = pkin(3) * t239;
t292 = pkin(3) * t170;
t301 = t127 * t292 + t134 * t226;
t236 = qJD(1) * qJD(2);
t300 = qJ(2) * qJDD(1) + t236;
t121 = t127 * pkin(4);
t279 = t41 * qJ(5);
t299 = t279 - t121;
t175 = cos(qJ(1));
t253 = t175 * t174;
t172 = sin(qJ(1));
t260 = t172 * t171;
t102 = t169 * t260 + t253;
t254 = t175 * t171;
t259 = t172 * t174;
t104 = -t169 * t254 + t259;
t298 = -g(1) * t104 + g(2) * t102;
t167 = qJ(3) + qJ(4);
t155 = sin(t167);
t287 = g(3) * t168;
t156 = cos(t167);
t255 = t175 * t156;
t262 = t172 * t155;
t93 = t169 * t262 + t255;
t256 = t175 * t155;
t261 = t172 * t156;
t95 = -t169 * t256 + t261;
t297 = -g(1) * t95 + g(2) * t93 + t155 * t287;
t296 = t86 ^ 2;
t295 = pkin(7) + pkin(6);
t291 = pkin(7) * t168;
t289 = g(1) * t172;
t161 = g(2) * t175;
t286 = t171 * pkin(3);
t285 = t173 * pkin(3);
t159 = t174 * pkin(3);
t187 = qJD(1) * t114;
t83 = t168 * t187;
t284 = t86 * t83;
t15 = -t134 * pkin(4) + t18;
t283 = -t18 + t15;
t204 = t229 * t174;
t233 = qJDD(1) * t174;
t211 = t170 * t233;
t251 = t229 * t198;
t42 = (t211 + (qJD(1) * t204 + t234) * t173) * t168 - t251;
t282 = -t83 * t226 - t42 * t292;
t35 = t173 * t61 - t54;
t110 = t174 * t120;
t65 = -t227 - t110 + (-pkin(3) - t272) * t169;
t266 = t168 * t171;
t76 = -pkin(7) * t266 + t80;
t40 = t170 * t65 + t173 * t76;
t263 = t170 * t171;
t113 = t258 - t263;
t281 = (t229 - t238) * t113;
t280 = t169 * t187 - t303;
t56 = t173 * t62;
t278 = t42 * qJ(5);
t74 = -t171 * t221 + t91;
t277 = t74 * t142;
t276 = t75 * t142;
t275 = t83 * qJ(5);
t274 = t83 * t134;
t273 = t86 * t134;
t270 = qJDD(1) * pkin(1);
t119 = pkin(4) * t156 + t159;
t269 = (pkin(2) + t119) * t169;
t268 = t127 * t169;
t163 = t168 ^ 2;
t176 = qJD(1) ^ 2;
t267 = t163 * t176;
t264 = t168 * t175;
t118 = pkin(4) * t155 + t286;
t257 = t175 * t118;
t101 = pkin(3) * t220 + qJ(2) * t245;
t210 = -t83 * pkin(4) - qJD(5);
t59 = -t210 + t101;
t252 = qJD(5) + t59;
t177 = qJ(2) ^ 2;
t202 = t236 * t228;
t232 = t163 * qJDD(1);
t250 = t163 * t202 + t177 * t232;
t241 = qJD(3) * t174;
t243 = qJD(2) * t169;
t249 = -t120 * t241 + t174 * t243;
t107 = (pkin(3) * t241 + qJD(2)) * t168;
t143 = pkin(3) * t266;
t112 = t168 * qJ(2) + t143;
t157 = t172 * qJ(2);
t248 = t175 * pkin(1) + t157;
t164 = t169 ^ 2;
t247 = t163 + t164;
t165 = t171 ^ 2;
t166 = t174 ^ 2;
t246 = t165 - t166;
t242 = qJD(3) * t100;
t240 = qJD(4) * t170;
t235 = qJD(1) * qJD(3);
t231 = t168 * qJDD(1);
t225 = pkin(3) * t240;
t224 = t168 * t263;
t218 = t171 * t243;
t217 = qJ(2) * t230;
t216 = t171 * t236;
t215 = t174 * t236;
t214 = t169 * t235;
t213 = t174 * t235;
t209 = -t161 + t289;
t197 = qJ(2) * t214;
t99 = -t120 * qJDD(1) + qJDD(2);
t90 = t174 * t99;
t25 = -t140 * pkin(3) + t90 + (-pkin(7) * t231 - t197) * t174 + (-t217 - t242 + (qJD(3) * t291 - t243) * qJD(1)) * t171;
t188 = t213 + t234;
t206 = t100 * t241 + t169 * t215 + t171 * t99 + t174 * t217;
t37 = -t171 * t197 + t206;
t33 = -t188 * t291 + t37;
t208 = -t170 * t33 + t173 * t25;
t34 = -t170 * t61 - t56;
t39 = -t170 * t76 + t173 * t65;
t5 = t170 * t25 + t173 * t33 + t47 * t239 - t62 * t240;
t207 = t247 * t176;
t205 = qJD(1) * (-qJD(3) - t142);
t203 = pkin(3) * t219;
t152 = qJDD(2) - t270;
t201 = t140 + t230;
t200 = t174 * t171 * t267;
t70 = qJ(2) * t231 + qJD(3) * t203 + qJDD(1) * t143 + t168 * t236;
t196 = t171 * t213;
t195 = -g(1) * t93 - g(2) * t95;
t94 = -t169 * t261 + t256;
t96 = t169 * t255 + t262;
t194 = -g(1) * t94 - g(2) * t96;
t144 = t168 * t289;
t193 = -g(2) * t264 + t144;
t192 = g(1) * t175 + g(2) * t172;
t31 = t170 * t47 + t56;
t191 = qJD(3) * (t142 + t238);
t190 = t152 - t270 + t161;
t57 = t186 * qJD(3) + t249;
t58 = -t218 + (-t141 + (t120 + t291) * t171) * qJD(3);
t9 = t170 * t58 + t173 * t57 + t65 * t239 - t76 * t240;
t189 = t214 + t267;
t22 = t42 * pkin(4) + qJDD(5) + t70;
t185 = -t142 ^ 2 - t267;
t184 = g(1) * t96 - g(2) * t94 + t156 * t287 - t5;
t6 = -t31 * qJD(4) + t208;
t10 = -t40 * qJD(4) - t170 * t57 + t173 * t58;
t183 = t101 * t83 + t184;
t181 = t252 * t83 + t184 + t278;
t180 = t6 + t297;
t179 = -t101 * t86 + t180;
t162 = -qJ(5) - t295;
t158 = t175 * qJ(2);
t151 = pkin(4) + t285;
t149 = qJ(2) + t286;
t105 = t169 * t253 + t260;
t103 = -t169 * t259 + t254;
t98 = t222 - t224;
t97 = t114 * t168;
t92 = t295 * t168 + pkin(1) + (pkin(2) + t159) * t169;
t81 = t83 ^ 2;
t79 = -t110 - t223;
t71 = t86 * pkin(4) + t203;
t68 = -t218 - t302;
t67 = -qJD(3) * t223 + t249;
t66 = t97 * pkin(4) + t112;
t53 = t173 * t168 * t204 - t229 * t224;
t52 = t303 * t168;
t44 = t53 * pkin(4) + t107;
t43 = -t81 + t296;
t38 = -t171 * t242 + t90 + (-qJ(2) * t188 - t216) * t169;
t32 = -t97 * qJ(5) + t40;
t28 = -t169 * pkin(4) - t98 * qJ(5) + t39;
t27 = -t273 + (-t211 + (-t229 * t244 - t234) * t173) * t168 + t251;
t26 = -t41 - t274;
t21 = -t78 + t35;
t20 = t34 + t275;
t19 = t31 - t275;
t17 = -t113 * t127 - t280 * t134 - t83 * t245;
t16 = t114 * t127 + t281 * t134 - t86 * t245;
t14 = t42 * t97 + t83 * t53;
t13 = -t41 * t98 - t86 * t52;
t12 = t97 * t127 + t53 * t134 + t42 * t169;
t11 = -t98 * t127 + t52 * t134 + t41 * t169;
t8 = t52 * qJ(5) - t98 * qJD(5) + t10;
t7 = -t53 * qJ(5) - t97 * qJD(5) + t9;
t4 = t41 * t97 - t98 * t42 + t52 * t83 - t86 * t53;
t3 = t113 * t41 - t114 * t42 - t280 * t86 - t281 * t83;
t2 = -t83 * qJD(5) - t278 + t5;
t1 = -t86 * qJD(5) + t299 + t6;
t23 = [0, 0, 0, 0, 0, qJDD(1), t209, t192, 0, 0, t232, 0.2e1 * t168 * t230, 0, t164 * qJDD(1), 0, 0, (-t190 + t289) * t169, t168 * t190 - t144, 0.2e1 * t300 * t247 - t192, -t152 * pkin(1) - g(1) * (-t172 * pkin(1) + t158) - g(2) * t248 + (qJDD(1) * t177 + t202) * t164 + t250, (qJDD(1) * t166 - 0.2e1 * t196) * t163, 0.2e1 * (-t171 * t233 + t246 * t235) * t163, (t171 * t191 - t201 * t174) * t168, (qJDD(1) * t165 + 0.2e1 * t196) * t163, (t201 * t171 + t174 * t191) * t168, t140 * t169, -g(1) * t103 - g(2) * t105 - t79 * t140 - t68 * t142 - t38 * t169 + (t188 * t228 + 0.2e1 * t216) * t163, -g(1) * t102 - g(2) * t104 + t80 * t140 + t67 * t142 + t37 * t169 + (0.2e1 * t215 + (-t171 * t235 + t233) * t228) * t163, t144 + (-t161 + (-qJD(3) * t75 - qJDD(1) * t79 - t38 + (-t68 - t302) * qJD(1)) * t174 + (qJD(3) * t74 - qJDD(1) * t80 - t37 + (qJD(3) * t79 - t67) * qJD(1)) * t171) * t168, t37 * t80 + t75 * t67 + t38 * t79 + t74 * t68 - g(1) * (-t120 * t172 + t158) - g(2) * (t120 * t175 + t157) + t250, t13, t4, t11, t14, t12, t268, -t10 * t134 + t101 * t53 + t107 * t83 + t112 * t42 - t39 * t127 - t6 * t169 + t70 * t97 + t194, -t101 * t52 + t107 * t86 - t112 * t41 + t40 * t127 + t9 * t134 + t5 * t169 + t70 * t98 + t195, -t10 * t86 + t30 * t52 - t31 * t53 + t39 * t41 - t40 * t42 - t5 * t97 - t6 * t98 - t9 * t83 + t193, t5 * t40 + t31 * t9 + t6 * t39 + t30 * t10 + t70 * t112 + t101 * t107 - g(1) * (t149 * t175 - t92 * t172) - g(2) * (t149 * t172 + t92 * t175), t13, t4, t11, t14, t12, t268, -t1 * t169 - t28 * t127 - t8 * t134 + t22 * t97 + t66 * t42 + t44 * t83 + t59 * t53 + t194, t32 * t127 + t7 * t134 + t2 * t169 + t22 * t98 - t66 * t41 + t44 * t86 - t59 * t52 + t195, -t1 * t98 + t15 * t52 - t19 * t53 - t2 * t97 + t28 * t41 - t32 * t42 - t7 * t83 - t8 * t86 + t193, t2 * t32 + t19 * t7 + t1 * t28 + t15 * t8 + t22 * t66 + t59 * t44 - g(1) * (t158 + t257) - g(2) * (-t162 * t264 + t175 * t269 + t248) + (-g(1) * (t162 * t168 - pkin(1) - t269) - g(2) * t118) * t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t230, t231, -t207, -qJ(2) * t207 + t152 - t209, 0, 0, 0, 0, 0, 0, -t174 * t140 + t171 * t185, t171 * t140 + t174 * t185, (-t165 - t166) * t231, -qJ(2) * t267 + (t38 - t276) * t174 + (t37 + t277) * t171 - t209, 0, 0, 0, 0, 0, 0, t17, t16, t3, -t101 * t245 + t6 * t113 + t5 * t114 + t280 * t30 + t281 * t31 - t209, 0, 0, 0, 0, 0, 0, t17, t16, t3, t1 * t113 + t2 * t114 + t280 * t15 + t281 * t19 - t59 * t245 - t209; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t200, -t246 * t267, (t171 * t205 + t233) * t168, -t200, (t174 * t205 - t234) * t168, -t140, -t276 + t90 - t189 * t271 + (-t300 * t169 - t242 + t287) * t171 + t298, g(1) * t105 - g(2) * t103 + g(3) * t265 + t189 * t272 - t206 - t277, 0, 0, t284, t43, t26, -t284, t27, -t127, t34 * t134 + (-t127 * t173 + t134 * t240 - t219 * t83) * pkin(3) + t179, -t35 * t134 - t203 * t86 + t183 + t301, t41 * t285 + (-t30 + t35) * t83 + (t31 + t34 + t225) * t86 + t282, -t30 * t34 - t31 * t35 + (t5 * t170 + t6 * t173 + (g(3) * t171 - t101 * t244) * t168 + (-t30 * t170 + t31 * t173) * qJD(4) + t298) * pkin(3), t284, t43, t26, -t284, t27, -t127, -t151 * t127 + t20 * t134 - t71 * t83 - t252 * t86 + (-t56 + (pkin(3) * t134 - t47) * t170) * qJD(4) + t208 + t297 + t299, -t21 * t134 - t71 * t86 + t181 + t301, t151 * t41 + (-t15 + t21) * t83 + (t19 + t20 + t225) * t86 + t282, t1 * t151 - t19 * t21 - t15 * t20 - t59 * t71 - g(1) * (t172 * t119 - t169 * t257) - g(2) * (-t172 * t169 * t118 - t175 * t119) + t118 * t287 + (t2 * t170 + (-t15 * t170 + t173 * t19) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t284, t43, t26, -t284, t27, -t127, -t31 * t134 + t179, -t30 * t134 + t183, 0, 0, t284, t43, t26, -t284, t27, -t127, t279 - t19 * t134 - 0.2e1 * t121 + (t210 - t59) * t86 + t180, -t296 * pkin(4) - t18 * t134 + t181, t41 * pkin(4) - t283 * t83, t283 * t19 + (-t59 * t86 + t1 + t297) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 - t273, -t41 + t274, -t81 - t296, g(3) * t169 + t15 * t86 - t168 * t192 + t19 * t83 + t22;];
tau_reg = t23;
