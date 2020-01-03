% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPPR7
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPPR7_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR7_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR7_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:36:23
% EndTime: 2019-12-31 19:36:33
% DurationCPUTime: 4.83s
% Computational Cost: add. (8053->351), mult. (18945->447), div. (0->0), fcn. (12643->8), ass. (0->201)
t177 = sin(qJ(2));
t180 = cos(qJ(2));
t181 = qJD(2) ^ 2;
t173 = sin(pkin(8));
t174 = cos(pkin(8));
t226 = t177 * t174;
t149 = (t180 * t173 + t226) * qJD(1);
t252 = t149 ^ 2;
t129 = t252 + t181;
t222 = qJD(1) * t177;
t147 = -t174 * t180 * qJD(1) + t173 * t222;
t231 = t149 * t147;
t273 = qJDD(2) + t231;
t289 = t173 * t273;
t71 = t174 * t129 + t289;
t287 = t174 * t273;
t73 = -t173 * t129 + t287;
t307 = pkin(6) * (t177 * t71 - t180 * t73);
t306 = pkin(2) * t71;
t305 = qJ(3) * t71;
t304 = qJ(3) * t73;
t130 = t252 - t181;
t274 = qJDD(2) - t231;
t286 = t174 * t274;
t288 = t173 * t274;
t302 = t177 * (t173 * t130 + t286) - t180 * (t174 * t130 - t288);
t163 = t177 * qJDD(1);
t216 = qJD(1) * qJD(2);
t211 = t180 * t216;
t156 = t163 + t211;
t164 = t180 * qJDD(1);
t212 = t177 * t216;
t157 = t164 - t212;
t115 = t174 * t156 + t173 * t157;
t218 = t147 * qJD(2);
t257 = t115 + t218;
t114 = t173 * t156 - t174 * t157;
t221 = qJD(2) * t149;
t87 = t114 - t221;
t292 = -t173 * t87 - t174 * t257;
t301 = pkin(2) * t292;
t300 = qJ(3) * t292;
t253 = t147 ^ 2;
t254 = -t253 - t252;
t291 = t173 * t257 - t174 * t87;
t297 = -pkin(2) * t254 + qJ(3) * t291;
t125 = t253 - t181;
t296 = t177 * (-t174 * t125 + t289) - t180 * (t173 * t125 + t287);
t295 = pkin(6) * (-t177 * t292 + t180 * t291) - pkin(1) * t254;
t98 = -t181 - t253;
t63 = t173 * t98 + t286;
t66 = -t174 * t98 + t288;
t290 = pkin(6) * (t177 * t63 + t180 * t66);
t258 = t115 - t218;
t275 = t114 + t221;
t285 = t177 * (-t173 * t258 - t174 * t275) + t180 * (-t173 * t275 + t174 * t258);
t284 = pkin(2) * t63;
t283 = qJ(3) * t63;
t282 = qJ(3) * t66;
t259 = t258 * qJ(4);
t251 = 2 * qJD(4);
t171 = t180 ^ 2;
t182 = qJD(1) ^ 2;
t195 = qJD(2) * pkin(2) - qJ(3) * t222;
t178 = sin(qJ(1));
t249 = cos(qJ(1));
t208 = t178 * g(1) - t249 * g(2);
t196 = qJDD(1) * pkin(1) + t208;
t82 = t157 * pkin(2) + (qJ(3) * t171 + pkin(6)) * t182 - t195 * t222 - qJDD(3) + t196;
t269 = pkin(3) * t221 - t149 * t251 - t82;
t176 = sin(qJ(5));
t109 = qJDD(5) + t115;
t179 = cos(qJ(5));
t119 = t176 * qJD(2) - t179 * t147;
t121 = t179 * qJD(2) + t176 * t147;
t85 = t121 * t119;
t260 = t109 - t85;
t264 = t176 * t260;
t263 = t179 * t260;
t225 = t177 * t182;
t198 = t249 * g(1) + t178 * g(2);
t234 = qJDD(1) * pkin(6);
t153 = -t182 * pkin(1) - t198 + t234;
t228 = t177 * t153;
t79 = qJDD(2) * pkin(2) - t156 * qJ(3) - t228 + (pkin(2) * t225 + qJ(3) * t216 - g(3)) * t180;
t123 = -t177 * g(3) + t180 * t153;
t166 = t171 * t182;
t80 = -pkin(2) * t166 + t157 * qJ(3) - qJD(2) * t195 + t123;
t207 = t173 * t80 - t174 * t79;
t193 = -qJDD(2) * pkin(3) - t181 * qJ(4) + qJDD(4) + t207;
t190 = -qJDD(2) * pkin(7) + t193;
t262 = t257 * pkin(4) + t190;
t261 = 2 * qJD(3);
t44 = -0.2e1 * qJD(3) * t147 + t173 * t79 + t174 * t80;
t255 = t252 - t253;
t117 = t119 ^ 2;
t118 = t121 ^ 2;
t142 = qJD(5) + t149;
t140 = t142 ^ 2;
t250 = pkin(3) + pkin(7);
t248 = pkin(3) * t173;
t247 = pkin(3) * t174;
t245 = t173 * t82;
t242 = t174 * t82;
t124 = t149 * pkin(4) - qJD(2) * pkin(7);
t97 = t147 * pkin(3) - t149 * qJ(4);
t194 = -t181 * pkin(3) - t147 * t97 + t44;
t215 = qJDD(2) * qJ(4);
t24 = t215 - t114 * pkin(4) - t253 * pkin(7) + (t251 + t124) * qJD(2) + t194;
t239 = t176 * t24;
t61 = t109 + t85;
t238 = t176 * t61;
t43 = t149 * t261 + t207;
t20 = t173 * t44 - t174 * t43;
t237 = t177 * t20;
t236 = t179 * t24;
t235 = t179 * t61;
t233 = t142 * t176;
t232 = t142 * t179;
t162 = t180 * t225;
t227 = t177 * (qJDD(2) + t162);
t224 = t180 * (qJDD(2) - t162);
t223 = t261 + t97;
t217 = qJD(5) + t142;
t214 = t173 * t85;
t213 = t174 * t85;
t210 = qJ(4) * t173 + pkin(2);
t21 = t173 * t43 + t174 * t44;
t197 = (pkin(7) * t147 + t223) * t149;
t183 = -t259 + t269;
t27 = -pkin(4) * t253 + t250 * t114 - t149 * t124 + t183;
t12 = t176 * t27 - t179 * (t197 + t262);
t122 = t180 * g(3) + t228;
t206 = t177 * t122 + t180 * t123;
t203 = t176 * qJDD(2) - t179 * t114;
t189 = t176 * t197 + t179 * t27;
t13 = t176 * t262 + t189;
t6 = -t179 * t12 + t176 * t13;
t7 = t176 * t12 + t179 * t13;
t199 = t179 * qJDD(2) + t176 * t114;
t192 = (-qJD(5) + t142) * t121 - t203;
t70 = -t119 * qJD(5) + t199;
t191 = qJD(2) * t251 + t194;
t36 = t223 * t149 + t193;
t35 = t191 + t215;
t187 = t177 * (t174 * t115 - t173 * t221) + t180 * (t173 * t115 + t174 * t221);
t186 = t177 * (t173 * t114 + t174 * t218) + t180 * (-t174 * t114 + t173 * t218);
t185 = (t177 * (-t147 * t174 + t149 * t173) + t180 * (-t147 * t173 - t149 * t174)) * qJD(2);
t184 = -t114 * pkin(3) - t269;
t170 = t177 ^ 2;
t165 = t170 * t182;
t158 = t164 - 0.2e1 * t212;
t155 = t163 + 0.2e1 * t211;
t152 = t182 * pkin(6) + t196;
t96 = t142 * t119;
t95 = -t118 + t140;
t94 = t117 - t140;
t81 = t118 - t117;
t75 = -t118 - t140;
t69 = -t121 * qJD(5) - t203;
t68 = -t140 - t117;
t67 = -t117 - t118;
t59 = (t119 * t176 + t121 * t179) * t142;
t54 = -t217 * t119 + t199;
t53 = t70 + t96;
t52 = t70 - t96;
t49 = t217 * t121 + t203;
t48 = -t121 * t232 - t176 * t70;
t47 = -t119 * t233 - t179 * t69;
t46 = -t176 * t94 - t235;
t45 = -t179 * t95 - t264;
t41 = -t176 * t75 - t235;
t40 = t179 * t75 - t238;
t39 = t179 * t68 - t264;
t38 = t176 * t68 + t263;
t37 = t184 + t259;
t34 = (t114 + t275) * pkin(3) + t183;
t33 = t184 + 0.2e1 * t259;
t32 = t176 * t53 + t179 * t192;
t31 = t176 * t192 - t179 * t53;
t30 = t176 * t49 - t179 * t52;
t29 = -qJ(4) * t254 + t36;
t28 = -pkin(3) * t254 + t35;
t26 = t173 * t40 + t174 * t54;
t25 = t173 * t54 - t174 * t40;
t23 = t173 * t38 + t174 * t49;
t22 = t173 * t49 - t174 * t38;
t18 = t173 * t31 + t174 * t67;
t17 = t173 * t67 - t174 * t31;
t15 = t173 * t35 - t174 * t36;
t14 = pkin(4) * t31 - qJ(4) * t32;
t11 = pkin(4) * t54 - t250 * t41 - t239;
t10 = pkin(4) * t49 - t250 * t39 + t236;
t9 = -t176 * t190 - qJ(4) * t41 + (-t176 * t257 + t40) * pkin(4) - t189;
t8 = pkin(4) * t38 - qJ(4) * t39 - t12;
t5 = t173 * t6 + t174 * t24;
t4 = t173 * t24 - t174 * t6;
t3 = pkin(4) * t67 - t250 * t32 - t7;
t2 = pkin(4) * t6 - qJ(4) * t7;
t1 = pkin(4) * t24 - t250 * t7;
t16 = [0, 0, 0, 0, 0, qJDD(1), t208, t198, 0, 0, (t156 + t211) * t177, t180 * t155 + t177 * t158, t227 + t180 * (-t165 + t181), (t157 - t212) * t180, t177 * (t166 - t181) + t224, 0, t180 * t152 + pkin(1) * t158 + pkin(6) * (t180 * (-t166 - t181) - t227), -t177 * t152 - pkin(1) * t155 + pkin(6) * (-t224 - t177 * (-t165 - t181)), pkin(1) * (t165 + t166) + (t170 + t171) * t234 + t206, pkin(1) * t152 + pkin(6) * t206, t187, t285, t302, t186, -t296, t185, t177 * (-t245 - t283) + t180 * (-pkin(2) * t275 + t242 - t282) - pkin(1) * t275 - t290, t177 * (-t242 + t305) + t180 * (-pkin(2) * t258 - t245 - t304) - pkin(1) * t258 + t307, t177 * (-t20 - t300) + t180 * (t21 + t297) + t295, -qJ(3) * t237 + t180 * (pkin(2) * t82 + qJ(3) * t21) + pkin(1) * t82 + pkin(6) * (t180 * t21 - t237), t185, -t302, t296, t187, t285, t186, t177 * (-t173 * t28 + t174 * t29 - t300) + t180 * (t173 * t29 + t174 * t28 + t297) + t295, t177 * (-t173 * t34 + t283) + t180 * (t174 * t34 + t282) + t290 + (qJ(4) * t226 + t180 * t210 + pkin(1)) * t275, t177 * (t174 * t33 - t305) + t180 * (t173 * t33 + t304) - t307 + (-t177 * t248 + t180 * (pkin(2) + t247) + pkin(1)) * t258, (t177 * (qJ(4) * t174 - t248) + t180 * (t210 + t247) + pkin(1)) * t37 + (pkin(6) + qJ(3)) * (-t177 * t15 + t180 * (t173 * t36 + t174 * t35)), t177 * (-t173 * t48 + t213) + t180 * (t174 * t48 + t214), t177 * (-t173 * t30 + t174 * t81) + t180 * (t173 * t81 + t174 * t30), t177 * (-t173 * t45 + t174 * t53) + t180 * (t173 * t53 + t174 * t45), t177 * (-t173 * t47 - t213) + t180 * (t174 * t47 - t214), t177 * (-t173 * t46 + t174 * t192) + t180 * (t173 * t192 + t174 * t46), t177 * (t174 * t109 - t173 * t59) + t180 * (t173 * t109 + t174 * t59), t177 * (-qJ(3) * t22 - t173 * t10 + t174 * t8) + t180 * (-pkin(2) * t39 + qJ(3) * t23 + t174 * t10 + t173 * t8) - pkin(1) * t39 + pkin(6) * (-t177 * t22 + t180 * t23), t177 * (-qJ(3) * t25 - t173 * t11 + t174 * t9) + t180 * (-pkin(2) * t41 + qJ(3) * t26 + t174 * t11 + t173 * t9) - pkin(1) * t41 + pkin(6) * (-t177 * t25 + t180 * t26), t177 * (-qJ(3) * t17 + t174 * t14 - t173 * t3) + t180 * (-pkin(2) * t32 + qJ(3) * t18 + t173 * t14 + t174 * t3) - pkin(1) * t32 + pkin(6) * (-t177 * t17 + t180 * t18), t177 * (-qJ(3) * t4 - t173 * t1 + t174 * t2) + t180 * (-pkin(2) * t7 + qJ(3) * t5 + t174 * t1 + t173 * t2) - pkin(1) * t7 + pkin(6) * (-t177 * t4 + t180 * t5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t162, t165 - t166, t163, t162, t164, qJDD(2), -t122, -t123, 0, 0, t231, t255, t257, -t231, -t87, qJDD(2), -t43 + t284, -t44 - t306, t301, pkin(2) * t20, qJDD(2), -t257, t87, t231, t255, -t231, -pkin(3) * t257 - qJ(4) * t87 + t301, -pkin(3) * t274 - qJ(4) * t98 - t284 + t36, t306 + pkin(3) * t129 + (qJDD(2) + t273) * qJ(4) + t191, pkin(2) * t15 - pkin(3) * t36 + qJ(4) * t35, -t121 * t233 + t179 * t70, -t176 * t52 - t179 * t49, -t176 * t95 + t263, t119 * t232 - t176 * t69, t179 * t94 - t238, (-t119 * t179 + t121 * t176) * t142, pkin(2) * t22 + qJ(4) * t49 - t250 * t38 + t239, pkin(2) * t25 + qJ(4) * t54 - t250 * t40 + t236, pkin(2) * t17 + qJ(4) * t67 - t250 * t31 - t6, pkin(2) * t4 + qJ(4) * t24 - t250 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t275, t258, t254, -t82, 0, 0, 0, 0, 0, 0, t254, -t275, -t258, -t37, 0, 0, 0, 0, 0, 0, t39, t41, t32, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t257, t274, -t129, t36, 0, 0, 0, 0, 0, 0, t38, t40, t31, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t81, t53, -t85, t192, t109, -t12, -t13, 0, 0;];
tauJ_reg = t16;
