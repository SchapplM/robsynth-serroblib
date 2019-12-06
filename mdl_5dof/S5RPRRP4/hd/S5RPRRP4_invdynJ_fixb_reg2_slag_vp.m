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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:07:00
% EndTime: 2019-12-05 18:07:10
% DurationCPUTime: 3.93s
% Computational Cost: add. (4081->404), mult. (9951->495), div. (0->0), fcn. (7055->10), ass. (0->234)
t167 = sin(pkin(8));
t169 = sin(qJ(4));
t172 = cos(qJ(4));
t173 = cos(qJ(3));
t258 = t172 * t173;
t224 = t167 * t258;
t170 = sin(qJ(3));
t236 = qJDD(1) * t170;
t113 = t169 * t173 + t172 * t170;
t231 = qJD(3) + qJD(4);
t304 = t231 * t113;
t41 = -qJDD(1) * t224 + t167 * (qJD(1) * t304 + t169 * t236);
t171 = sin(qJ(1));
t174 = cos(qJ(1));
t290 = g(2) * t174;
t196 = g(3) * t171 + t290;
t230 = 2 * qJ(2);
t168 = cos(pkin(8));
t240 = t168 * qJD(1);
t143 = -qJD(3) + t240;
t272 = qJ(2) * t170;
t225 = t168 * t272;
t264 = t167 * t173;
t229 = pkin(7) * t264;
t187 = -t225 - t229;
t194 = t168 * pkin(2) + t167 * pkin(6) + pkin(1);
t99 = -qJD(1) * t194 + qJD(2);
t91 = t173 * t99;
t61 = qJD(1) * t187 + t91;
t47 = -t143 * pkin(3) + t61;
t246 = qJD(1) * t167;
t222 = t170 * t246;
t223 = qJ(2) * t240;
t75 = t170 * t99 + t173 * t223;
t62 = -pkin(7) * t222 + t75;
t54 = t169 * t62;
t30 = t172 * t47 - t54;
t201 = t169 * t222;
t245 = qJD(1) * t173;
t221 = t167 * t245;
t86 = t172 * t221 - t201;
t78 = t86 * qJ(5);
t18 = t30 - t78;
t271 = qJ(2) * t173;
t142 = t168 * t271;
t80 = -t170 * t194 + t142;
t303 = qJD(3) * t80;
t232 = t168 * qJDD(1);
t141 = -qJDD(3) + t232;
t126 = -qJDD(4) + t141;
t133 = -qJD(4) + t143;
t241 = qJD(4) * t172;
t227 = pkin(3) * t241;
t294 = pkin(3) * t169;
t302 = t126 * t294 + t133 * t227;
t238 = qJD(1) * qJD(2);
t301 = qJ(2) * qJDD(1) + t238;
t120 = t126 * pkin(4);
t279 = t41 * qJ(5);
t300 = t279 - t120;
t254 = t174 * t173;
t260 = t171 * t170;
t101 = t168 * t260 + t254;
t255 = t174 * t170;
t259 = t171 * t173;
t103 = t168 * t255 - t259;
t299 = -g(2) * t101 + g(3) * t103;
t166 = qJ(3) + qJ(4);
t157 = sin(t166);
t292 = g(1) * t167;
t158 = cos(t166);
t256 = t174 * t158;
t262 = t171 * t157;
t92 = t168 * t262 + t256;
t257 = t174 * t157;
t261 = t171 * t158;
t94 = t168 * t257 - t261;
t298 = -g(2) * t92 + g(3) * t94 + t157 * t292;
t297 = t86 ^ 2;
t175 = -pkin(7) - pkin(6);
t293 = pkin(7) * t167;
t159 = t174 * qJ(2);
t288 = g(3) * t159;
t286 = t170 * pkin(3);
t285 = t172 * pkin(3);
t188 = qJD(1) * t113;
t83 = t167 * t188;
t284 = t86 * t83;
t15 = -t133 * pkin(4) + t18;
t283 = -t18 + t15;
t207 = t231 * t173;
t235 = qJDD(1) * t173;
t213 = t169 * t235;
t251 = t231 * t201;
t42 = (t213 + (qJD(1) * t207 + t236) * t172) * t167 - t251;
t282 = -t83 * t227 - t42 * t294;
t35 = t172 * t61 - t54;
t109 = t173 * t194;
t65 = -t229 - t109 + (-pkin(3) - t272) * t168;
t265 = t167 * t170;
t76 = -pkin(7) * t265 + t80;
t40 = t169 * t65 + t172 * t76;
t263 = t169 * t170;
t112 = t258 - t263;
t281 = (t231 - t240) * t112;
t280 = t168 * t188 - t304;
t56 = t172 * t62;
t278 = t42 * qJ(5);
t74 = -t170 * t223 + t91;
t277 = t74 * t143;
t276 = t75 * t143;
t275 = t83 * qJ(5);
t274 = t83 * t133;
t273 = t86 * t133;
t270 = qJD(3) * t99;
t269 = qJDD(1) * pkin(1);
t117 = pkin(4) * t157 + t286;
t268 = t117 * t168;
t267 = t126 * t168;
t162 = t167 ^ 2;
t176 = qJD(1) ^ 2;
t266 = t162 * t176;
t100 = pkin(3) * t222 + qJ(2) * t246;
t212 = -t83 * pkin(4) - qJD(5);
t59 = -t212 + t100;
t253 = qJD(5) + t59;
t243 = qJD(3) * t173;
t244 = qJD(2) * t168;
t252 = t173 * t244 - t194 * t243;
t177 = qJ(2) ^ 2;
t205 = t238 * t230;
t234 = t162 * qJDD(1);
t250 = t162 * t205 + t177 * t234;
t249 = t196 * t167;
t160 = t173 * pkin(3);
t118 = pkin(4) * t158 + t160;
t106 = (pkin(3) * t243 + qJD(2)) * t167;
t144 = pkin(3) * t265;
t111 = t167 * qJ(2) + t144;
t163 = t168 ^ 2;
t248 = t162 + t163;
t164 = t170 ^ 2;
t165 = t173 ^ 2;
t247 = t164 - t165;
t242 = qJD(4) * t169;
t237 = qJD(1) * qJD(3);
t233 = t167 * qJDD(1);
t228 = pkin(3) * t242;
t226 = t167 * t263;
t220 = t170 * t244;
t219 = qJ(2) * t232;
t218 = t170 * t238;
t217 = t173 * t238;
t216 = t168 * t237;
t215 = t173 * t237;
t200 = qJ(2) * t216;
t98 = -qJDD(1) * t194 + qJDD(2);
t90 = t173 * t98;
t25 = -t141 * pkin(3) + t90 + (-pkin(7) * t233 - t200) * t173 + (-t219 - t270 + (qJD(3) * t293 - t244) * qJD(1)) * t170;
t189 = t215 + t236;
t209 = t168 * t217 + t170 * t98 + t173 * t219 + t99 * t243;
t37 = -t170 * t200 + t209;
t33 = -t189 * t293 + t37;
t211 = -t169 * t33 + t172 * t25;
t34 = -t169 * t61 - t56;
t39 = -t169 * t76 + t172 * t65;
t5 = t169 * t25 + t172 * t33 + t47 * t241 - t62 * t242;
t210 = t248 * t176;
t208 = qJD(1) * (-qJD(3) - t143);
t206 = pkin(3) * t221;
t204 = t141 + t232;
t203 = t173 * t170 * t266;
t70 = qJ(2) * t233 + qJD(3) * t206 + qJDD(1) * t144 + t167 * t238;
t199 = t170 * t215;
t198 = -g(2) * t94 - g(3) * t92;
t93 = t168 * t261 - t257;
t95 = -t168 * t256 - t262;
t197 = -g(2) * t95 + g(3) * t93;
t195 = g(2) * t171 - g(3) * t174;
t31 = t169 * t47 + t56;
t193 = qJD(3) * (t143 + t240);
t192 = (pkin(2) + t118) * t168 - (-qJ(5) + t175) * t167 + pkin(1);
t191 = (t160 + pkin(2)) * t168 - t167 * t175 + pkin(1);
t57 = qJD(3) * t187 + t252;
t58 = -t220 + (-t142 + (t194 + t293) * t170) * qJD(3);
t9 = t169 * t58 + t172 * t57 + t65 * t241 - t76 * t242;
t190 = t216 + t266;
t186 = -t196 - t269;
t22 = t42 * pkin(4) + qJDD(5) + t70;
t185 = -t143 ^ 2 - t266;
t184 = -g(2) * t93 - g(3) * t95 + t158 * t292 - t5;
t6 = -qJD(4) * t31 + t211;
t10 = -t40 * qJD(4) - t169 * t57 + t172 * t58;
t183 = t100 * t83 + t184;
t181 = t253 * t83 + t184 + t278;
t180 = t6 + t298;
t179 = -t100 * t86 + t180;
t154 = qJDD(2) - t269;
t152 = pkin(4) + t285;
t104 = -t168 * t254 - t260;
t102 = t168 * t259 - t255;
t97 = t224 - t226;
t96 = t113 * t167;
t81 = t83 ^ 2;
t79 = -t109 - t225;
t71 = t86 * pkin(4) + t206;
t68 = -t220 - t303;
t67 = -qJD(3) * t225 + t252;
t66 = t96 * pkin(4) + t111;
t53 = t167 * t172 * t207 - t231 * t226;
t52 = t304 * t167;
t44 = t53 * pkin(4) + t106;
t43 = -t81 + t297;
t38 = -t170 * t270 + t90 + (-qJ(2) * t189 - t218) * t168;
t32 = -t96 * qJ(5) + t40;
t28 = -t168 * pkin(4) - t97 * qJ(5) + t39;
t27 = -t273 + (-t213 + (-t231 * t245 - t236) * t172) * t167 + t251;
t26 = -t41 - t274;
t21 = -t78 + t35;
t20 = t34 + t275;
t19 = t31 - t275;
t17 = -t112 * t126 - t280 * t133 - t83 * t246;
t16 = t113 * t126 + t281 * t133 - t86 * t246;
t14 = t42 * t96 + t83 * t53;
t13 = -t41 * t97 - t86 * t52;
t12 = t96 * t126 + t53 * t133 + t42 * t168;
t11 = -t97 * t126 + t52 * t133 + t41 * t168;
t8 = t52 * qJ(5) - t97 * qJD(5) + t10;
t7 = -t53 * qJ(5) - t96 * qJD(5) + t9;
t4 = t41 * t96 - t97 * t42 + t52 * t83 - t86 * t53;
t3 = t112 * t41 - t113 * t42 - t280 * t86 - t281 * t83;
t2 = -t83 * qJD(5) - t278 + t5;
t1 = -t86 * qJD(5) + t300 + t6;
t23 = [0, 0, 0, 0, 0, qJDD(1), t196, -t195, 0, 0, t234, 0.2e1 * t167 * t232, 0, t163 * qJDD(1), 0, 0, (-t154 - t186) * t168, (t154 - t269) * t167 - t249, 0.2e1 * t301 * t248 + t195, -t154 * pkin(1) - g(2) * (-t174 * pkin(1) - t171 * qJ(2)) - g(3) * (-t171 * pkin(1) + t159) + (qJDD(1) * t177 + t205) * t163 + t250, (qJDD(1) * t165 - 0.2e1 * t199) * t162, 0.2e1 * (-t170 * t235 + t247 * t237) * t162, (t170 * t193 - t173 * t204) * t167, (qJDD(1) * t164 + 0.2e1 * t199) * t162, (t170 * t204 + t173 * t193) * t167, t141 * t168, -g(2) * t104 + g(3) * t102 - t79 * t141 - t68 * t143 - t38 * t168 + (t189 * t230 + 0.2e1 * t218) * t162, -g(2) * t103 - g(3) * t101 + t80 * t141 + t67 * t143 + t37 * t168 + (0.2e1 * t217 + (-t170 * t237 + t235) * t230) * t162, ((-qJD(3) * t75 - qJDD(1) * t79 - t38 + (-t68 - t303) * qJD(1)) * t173 + (qJD(3) * t74 - qJDD(1) * t80 - t37 + (qJD(3) * t79 - t67) * qJD(1)) * t170) * t167 + t249, -t288 + t37 * t80 + t38 * t79 + t75 * t67 + t74 * t68 + t194 * t290 + (g(2) * qJ(2) + g(3) * t194) * t171 + t250, t13, t4, t11, t14, t12, t267, -t10 * t133 + t100 * t53 + t106 * t83 + t111 * t42 - t39 * t126 - t6 * t168 + t70 * t96 + t197, -t100 * t52 + t106 * t86 - t111 * t41 + t40 * t126 + t9 * t133 + t5 * t168 + t70 * t97 + t198, -t10 * t86 + t30 * t52 - t31 * t53 + t39 * t41 - t40 * t42 - t5 * t96 - t6 * t97 - t9 * t83 + t249, -t288 + t30 * t10 + t100 * t106 + t70 * t111 + t31 * t9 + t6 * t39 + t5 * t40 + (g(2) * t191 - g(3) * t286) * t174 + (-g(2) * (-qJ(2) - t286) + g(3) * t191) * t171, t13, t4, t11, t14, t12, t267, -t1 * t168 - t28 * t126 - t8 * t133 + t22 * t96 + t66 * t42 + t44 * t83 + t59 * t53 + t197, t32 * t126 + t7 * t133 + t2 * t168 + t22 * t97 - t66 * t41 + t44 * t86 - t59 * t52 + t198, -t1 * t97 + t15 * t52 - t19 * t53 - t2 * t96 + t28 * t41 - t32 * t42 - t7 * t83 - t8 * t86 + t249, -t288 + t1 * t28 + t15 * t8 + t19 * t7 + t2 * t32 + t22 * t66 + t59 * t44 + (g(2) * t192 - g(3) * t117) * t174 + (-g(2) * (-qJ(2) - t117) + g(3) * t192) * t171; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t232, t233, -t210, -qJ(2) * t210 + qJDD(2) + t186, 0, 0, 0, 0, 0, 0, -t173 * t141 + t170 * t185, t170 * t141 + t173 * t185, (-t164 - t165) * t233, -qJ(2) * t266 + (t38 - t276) * t173 + (t37 + t277) * t170 - t196, 0, 0, 0, 0, 0, 0, t17, t16, t3, -t100 * t246 + t6 * t112 + t5 * t113 + t280 * t30 + t281 * t31 - t196, 0, 0, 0, 0, 0, 0, t17, t16, t3, t1 * t112 + t2 * t113 + t15 * t280 + t19 * t281 - t246 * t59 - t196; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t203, -t247 * t266, (t170 * t208 + t235) * t167, -t203, (t173 * t208 - t236) * t167, -t141, -t276 + t90 - t190 * t271 + (-t301 * t168 - t270 + t292) * t170 + t299, g(1) * t264 - g(2) * t102 - g(3) * t104 + t190 * t272 - t209 - t277, 0, 0, t284, t43, t26, -t284, t27, -t126, t34 * t133 + (-t126 * t172 + t133 * t242 - t221 * t83) * pkin(3) + t179, -t35 * t133 - t206 * t86 + t183 + t302, t41 * t285 + (-t30 + t35) * t83 + (t31 + t34 + t228) * t86 + t282, -t30 * t34 - t31 * t35 + (t5 * t169 + t6 * t172 + (g(1) * t170 - t100 * t245) * t167 + (-t169 * t30 + t172 * t31) * qJD(4) + t299) * pkin(3), t284, t43, t26, -t284, t27, -t126, -t152 * t126 + t20 * t133 - t71 * t83 - t253 * t86 + (-t56 + (pkin(3) * t133 - t47) * t169) * qJD(4) + t211 + t298 + t300, -t21 * t133 - t71 * t86 + t181 + t302, t152 * t41 + (-t15 + t21) * t83 + (t19 + t20 + t228) * t86 + t282, t1 * t152 - t19 * t21 - t15 * t20 - t59 * t71 + t117 * t292 - g(2) * (t118 * t174 + t171 * t268) - g(3) * (t118 * t171 - t174 * t268) + (t2 * t169 + (-t15 * t169 + t172 * t19) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t284, t43, t26, -t284, t27, -t126, -t31 * t133 + t179, -t133 * t30 + t183, 0, 0, t284, t43, t26, -t284, t27, -t126, t279 - t19 * t133 - 0.2e1 * t120 + (t212 - t59) * t86 + t180, -t297 * pkin(4) - t18 * t133 + t181, t41 * pkin(4) - t283 * t83, t283 * t19 + (-t59 * t86 + t1 + t298) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 - t273, -t41 + t274, -t81 - t297, g(1) * t168 + t15 * t86 + t167 * t195 + t19 * t83 + t22;];
tau_reg = t23;
