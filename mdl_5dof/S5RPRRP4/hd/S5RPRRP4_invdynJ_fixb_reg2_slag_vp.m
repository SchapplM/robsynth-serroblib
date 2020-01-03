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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:50:31
% EndTime: 2020-01-03 11:50:39
% DurationCPUTime: 4.00s
% Computational Cost: add. (4081->403), mult. (9951->499), div. (0->0), fcn. (7055->10), ass. (0->239)
t169 = sin(pkin(8));
t171 = sin(qJ(4));
t174 = cos(qJ(4));
t175 = cos(qJ(3));
t261 = t174 * t175;
t226 = t169 * t261;
t172 = sin(qJ(3));
t238 = qJDD(1) * t172;
t113 = t171 * t175 + t174 * t172;
t233 = qJD(3) + qJD(4);
t309 = t233 * t113;
t41 = -qJDD(1) * t226 + (qJD(1) * t309 + t171 * t238) * t169;
t275 = qJDD(1) * pkin(1);
t152 = qJDD(2) - t275;
t173 = sin(qJ(1));
t176 = cos(qJ(1));
t305 = -g(2) * t176 - g(3) * t173;
t308 = t305 - t152;
t232 = 2 * qJ(2);
t170 = cos(pkin(8));
t242 = t170 * qJD(1);
t143 = -qJD(3) + t242;
t278 = qJ(2) * t172;
t227 = t170 * t278;
t270 = t169 * t175;
t231 = pkin(7) * t270;
t189 = -t227 - t231;
t197 = t170 * pkin(2) + t169 * pkin(6);
t119 = -pkin(1) - t197;
t99 = t119 * qJD(1) + qJD(2);
t91 = t175 * t99;
t61 = t189 * qJD(1) + t91;
t47 = -t143 * pkin(3) + t61;
t248 = qJD(1) * t169;
t224 = t172 * t248;
t225 = qJ(2) * t242;
t75 = t172 * t99 + t175 * t225;
t62 = -pkin(7) * t224 + t75;
t54 = t171 * t62;
t30 = t174 * t47 - t54;
t202 = t171 * t224;
t247 = qJD(1) * t175;
t223 = t169 * t247;
t86 = t174 * t223 - t202;
t78 = t86 * qJ(5);
t18 = t30 - t78;
t277 = qJ(2) * t175;
t142 = t170 * t277;
t80 = t172 * t119 + t142;
t307 = qJD(3) * t80;
t234 = t170 * qJDD(1);
t141 = -qJDD(3) + t234;
t126 = -qJDD(4) + t141;
t133 = -qJD(4) + t143;
t243 = qJD(4) * t174;
t229 = pkin(3) * t243;
t297 = pkin(3) * t171;
t306 = t126 * t297 + t133 * t229;
t240 = qJD(1) * qJD(2);
t304 = qJ(2) * qJDD(1) + t240;
t120 = t126 * pkin(4);
t285 = t41 * qJ(5);
t303 = t285 - t120;
t257 = t176 * t175;
t263 = t173 * t172;
t101 = -t170 * t263 - t257;
t258 = t176 * t172;
t262 = t173 * t175;
t103 = t170 * t258 - t262;
t302 = -g(2) * t101 - g(3) * t103;
t168 = qJ(3) + qJ(4);
t155 = sin(t168);
t295 = g(1) * t169;
t156 = cos(t168);
t259 = t176 * t156;
t265 = t173 * t155;
t92 = -t170 * t265 - t259;
t260 = t176 * t155;
t264 = t173 * t156;
t94 = t170 * t260 - t264;
t301 = -g(2) * t92 - g(3) * t94 + t155 * t295;
t300 = t86 ^ 2;
t177 = -pkin(7) - pkin(6);
t296 = pkin(7) * t169;
t292 = t172 * pkin(3);
t291 = t174 * pkin(3);
t190 = qJD(1) * t113;
t83 = t169 * t190;
t290 = t86 * t83;
t15 = -t133 * pkin(4) + t18;
t289 = -t18 + t15;
t208 = t175 * t233;
t237 = qJDD(1) * t175;
t215 = t171 * t237;
t254 = t233 * t202;
t42 = (t215 + (qJD(1) * t208 + t238) * t174) * t169 - t254;
t288 = -t83 * t229 - t42 * t297;
t35 = t174 * t61 - t54;
t109 = t175 * t119;
t65 = -t231 + t109 + (-pkin(3) - t278) * t170;
t271 = t169 * t172;
t76 = -pkin(7) * t271 + t80;
t40 = t171 * t65 + t174 * t76;
t267 = t171 * t172;
t112 = t261 - t267;
t287 = (t233 - t242) * t112;
t286 = t170 * t190 - t309;
t56 = t174 * t62;
t284 = t42 * qJ(5);
t74 = -t172 * t225 + t91;
t283 = t74 * t143;
t282 = t75 * t143;
t281 = t83 * qJ(5);
t280 = t83 * t133;
t279 = t86 * t133;
t276 = qJD(3) * t99;
t274 = t126 * t170;
t273 = (-qJ(5) + t177) * t169;
t164 = t169 ^ 2;
t178 = qJD(1) ^ 2;
t272 = t164 * t178;
t269 = t169 * t177;
t268 = t170 * t173;
t117 = pkin(4) * t155 + t292;
t266 = t173 * t117;
t100 = pkin(3) * t224 + qJ(2) * t248;
t214 = -t83 * pkin(4) - qJD(5);
t59 = -t214 + t100;
t256 = qJD(5) + t59;
t245 = qJD(3) * t175;
t246 = qJD(2) * t170;
t255 = t119 * t245 + t175 * t246;
t179 = qJ(2) ^ 2;
t206 = t240 * t232;
t236 = t164 * qJDD(1);
t253 = t164 * t206 + t179 * t236;
t159 = t175 * pkin(3);
t118 = pkin(4) * t156 + t159;
t106 = (pkin(3) * t245 + qJD(2)) * t169;
t144 = pkin(3) * t271;
t111 = t169 * qJ(2) + t144;
t252 = t176 * pkin(1) + t173 * qJ(2);
t165 = t170 ^ 2;
t250 = t164 + t165;
t166 = t172 ^ 2;
t167 = t175 ^ 2;
t249 = t166 - t167;
t244 = qJD(4) * t171;
t239 = qJD(1) * qJD(3);
t235 = t169 * qJDD(1);
t230 = pkin(3) * t244;
t228 = t169 * t267;
t222 = t172 * t246;
t221 = qJ(2) * t234;
t220 = t172 * t240;
t219 = t175 * t240;
t218 = t170 * t239;
t217 = t175 * t239;
t201 = qJ(2) * t218;
t98 = t119 * qJDD(1) + qJDD(2);
t90 = t175 * t98;
t25 = -t141 * pkin(3) + t90 + (-pkin(7) * t235 - t201) * t175 + (-t221 - t276 + (qJD(3) * t296 - t246) * qJD(1)) * t172;
t191 = t217 + t238;
t210 = t170 * t219 + t172 * t98 + t175 * t221 + t99 * t245;
t37 = -t172 * t201 + t210;
t33 = -t191 * t296 + t37;
t213 = -t171 * t33 + t174 * t25;
t34 = -t171 * t61 - t56;
t39 = -t171 * t76 + t174 * t65;
t5 = t171 * t25 + t174 * t33 + t47 * t243 - t62 * t244;
t212 = t250 * t178;
t158 = t173 * pkin(1);
t211 = -t176 * qJ(2) + t158;
t209 = qJD(1) * (-qJD(3) - t143);
t207 = pkin(3) * t223;
t205 = t141 + t234;
t204 = t175 * t172 * t272;
t70 = qJ(2) * t235 + qJD(3) * t207 + qJDD(1) * t144 + t169 * t240;
t200 = t172 * t217;
t199 = g(2) * t94 - g(3) * t92;
t93 = t170 * t264 - t260;
t95 = t170 * t259 + t265;
t198 = -g(2) * t95 - g(3) * t93;
t195 = g(2) * t173 - g(3) * t176;
t31 = t171 * t47 + t56;
t194 = qJD(3) * (t143 + t242);
t193 = t305 * t169;
t57 = t189 * qJD(3) + t255;
t58 = -t222 + (-t142 + (-t119 + t296) * t172) * qJD(3);
t9 = t171 * t58 + t174 * t57 + t65 * t243 - t76 * t244;
t192 = t218 + t272;
t22 = t42 * pkin(4) + qJDD(5) + t70;
t188 = t275 + t308;
t187 = -t143 ^ 2 - t272;
t186 = g(2) * t93 - g(3) * t95 + t156 * t295 - t5;
t6 = -t31 * qJD(4) + t213;
t10 = -t40 * qJD(4) - t171 * t57 + t174 * t58;
t185 = t100 * t83 + t186;
t183 = t256 * t83 + t186 + t284;
t182 = t6 + t301;
t181 = -t100 * t86 + t182;
t151 = t159 + pkin(2);
t150 = pkin(4) + t291;
t114 = pkin(2) + t118;
t104 = t170 * t257 + t263;
t102 = t170 * t262 - t258;
t97 = t226 - t228;
t96 = t113 * t169;
t81 = t83 ^ 2;
t79 = t109 - t227;
t71 = t86 * pkin(4) + t207;
t68 = -t222 - t307;
t67 = -qJD(3) * t227 + t255;
t66 = t96 * pkin(4) + t111;
t53 = t174 * t169 * t208 - t233 * t228;
t52 = t309 * t169;
t44 = t53 * pkin(4) + t106;
t43 = -t81 + t300;
t38 = -t172 * t276 + t90 + (-t191 * qJ(2) - t220) * t170;
t32 = -t96 * qJ(5) + t40;
t28 = -t170 * pkin(4) - t97 * qJ(5) + t39;
t27 = -t279 + (-t215 + (-t233 * t247 - t238) * t174) * t169 + t254;
t26 = -t41 - t280;
t21 = -t78 + t35;
t20 = t34 + t281;
t19 = t31 - t281;
t17 = -t112 * t126 - t286 * t133 - t83 * t248;
t16 = t113 * t126 + t287 * t133 - t86 * t248;
t14 = t42 * t96 + t83 * t53;
t13 = -t41 * t97 - t86 * t52;
t12 = t96 * t126 + t53 * t133 + t42 * t170;
t11 = -t97 * t126 + t52 * t133 + t41 * t170;
t8 = t52 * qJ(5) - t97 * qJD(5) + t10;
t7 = -qJ(5) * t53 - qJD(5) * t96 + t9;
t4 = t41 * t96 - t42 * t97 + t52 * t83 - t53 * t86;
t3 = t112 * t41 - t113 * t42 - t286 * t86 - t287 * t83;
t2 = -qJD(5) * t83 - t284 + t5;
t1 = -t86 * qJD(5) + t303 + t6;
t23 = [0, 0, 0, 0, 0, qJDD(1), t305, t195, 0, 0, t236, 0.2e1 * t169 * t234, 0, t165 * qJDD(1), 0, 0, t188 * t170, -t188 * t169, 0.2e1 * t304 * t250 - t195, -t152 * pkin(1) - g(2) * t252 - g(3) * t211 + (qJDD(1) * t179 + t206) * t165 + t253, (qJDD(1) * t167 - 0.2e1 * t200) * t164, 0.2e1 * (-t172 * t237 + t249 * t239) * t164, (t172 * t194 - t205 * t175) * t169, (qJDD(1) * t166 + 0.2e1 * t200) * t164, (t205 * t172 + t175 * t194) * t169, t141 * t170, -g(2) * t104 - g(3) * t102 - t79 * t141 - t68 * t143 - t38 * t170 + (t191 * t232 + 0.2e1 * t220) * t164, g(2) * t103 - g(3) * t101 + t80 * t141 + t67 * t143 + t37 * t170 + (0.2e1 * t219 + (-t172 * t239 + t237) * t232) * t164, ((-qJD(3) * t75 - qJDD(1) * t79 - t38 + (-t68 - t307) * qJD(1)) * t175 + (qJD(3) * t74 - qJDD(1) * t80 - t37 + (qJD(3) * t79 - t67) * qJD(1)) * t172 + t305) * t169, t37 * t80 + t75 * t67 + t38 * t79 + t74 * t68 - g(2) * (t197 * t176 + t252) - g(3) * (t197 * t173 + t211) + t253, t13, t4, t11, t14, t12, t274, -t10 * t133 + t100 * t53 + t106 * t83 + t111 * t42 - t126 * t39 - t170 * t6 + t70 * t96 + t198, -t100 * t52 + t106 * t86 - t111 * t41 + t126 * t40 + t133 * t9 + t170 * t5 + t70 * t97 + t199, -t10 * t86 + t30 * t52 - t31 * t53 + t39 * t41 - t40 * t42 - t5 * t96 - t6 * t97 - t9 * t83 + t193, t5 * t40 + t31 * t9 + t6 * t39 + t30 * t10 + t70 * t111 + t100 * t106 - g(2) * (pkin(3) * t263 + t252) - g(3) * (t151 * t268 - t173 * t269 + t158) + (-g(2) * (t151 * t170 - t269) - g(3) * (-qJ(2) - t292)) * t176, t13, t4, t11, t14, t12, t274, -t1 * t170 - t126 * t28 - t133 * t8 + t22 * t96 + t42 * t66 + t44 * t83 + t53 * t59 + t198, t126 * t32 + t133 * t7 + t170 * t2 + t22 * t97 - t41 * t66 + t44 * t86 - t52 * t59 + t199, -t1 * t97 + t15 * t52 - t19 * t53 - t2 * t96 + t28 * t41 - t32 * t42 - t7 * t83 - t8 * t86 + t193, t2 * t32 + t19 * t7 + t1 * t28 + t15 * t8 + t22 * t66 + t59 * t44 - g(2) * (t252 + t266) - g(3) * (t114 * t268 - t173 * t273 + t158) + (-g(2) * (t114 * t170 - t273) - g(3) * (-qJ(2) - t117)) * t176; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t234, t235, -t212, -qJ(2) * t212 - t308, 0, 0, 0, 0, 0, 0, -t175 * t141 + t187 * t172, t172 * t141 + t187 * t175, (-t166 - t167) * t235, -qJ(2) * t272 + (t38 - t282) * t175 + (t37 + t283) * t172 - t305, 0, 0, 0, 0, 0, 0, t17, t16, t3, -t100 * t248 + t6 * t112 + t5 * t113 + t286 * t30 + t287 * t31 - t305, 0, 0, 0, 0, 0, 0, t17, t16, t3, t1 * t112 + t2 * t113 + t286 * t15 + t287 * t19 - t59 * t248 - t305; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t204, -t249 * t272, (t172 * t209 + t237) * t169, -t204, (t175 * t209 - t238) * t169, -t141, -t282 + t90 - t192 * t277 + (-t304 * t170 - t276 + t295) * t172 + t302, g(1) * t270 + g(2) * t102 - g(3) * t104 + t192 * t278 - t210 - t283, 0, 0, t290, t43, t26, -t290, t27, -t126, t34 * t133 + (-t126 * t174 + t133 * t244 - t83 * t223) * pkin(3) + t181, -t133 * t35 - t86 * t207 + t185 + t306, t41 * t291 + (-t30 + t35) * t83 + (t31 + t34 + t230) * t86 + t288, -t30 * t34 - t31 * t35 + (t5 * t171 + t6 * t174 + (g(1) * t172 - t100 * t247) * t169 + (-t171 * t30 + t174 * t31) * qJD(4) + t302) * pkin(3), t290, t43, t26, -t290, t27, -t126, -t150 * t126 + t20 * t133 - t71 * t83 - t256 * t86 + (-t56 + (pkin(3) * t133 - t47) * t171) * qJD(4) + t213 + t301 + t303, -t21 * t133 - t71 * t86 + t183 + t306, t150 * t41 + (-t15 + t21) * t83 + (t19 + t20 + t230) * t86 + t288, t1 * t150 - t19 * t21 - t15 * t20 - t59 * t71 + t117 * t295 - g(2) * (-t118 * t176 - t170 * t266) - g(3) * (t117 * t170 * t176 - t118 * t173) + (t2 * t171 + (-t15 * t171 + t174 * t19) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t290, t43, t26, -t290, t27, -t126, -t31 * t133 + t181, -t133 * t30 + t185, 0, 0, t290, t43, t26, -t290, t27, -t126, t285 - t19 * t133 - 0.2e1 * t120 + (t214 - t59) * t86 + t182, -t300 * pkin(4) - t18 * t133 + t183, t41 * pkin(4) - t289 * t83, t289 * t19 + (-t59 * t86 + t1 + t301) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42 - t279, -t41 + t280, -t81 - t300, g(1) * t170 + t15 * t86 - t169 * t195 + t19 * t83 + t22;];
tau_reg = t23;
