% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RPPRR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR4_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:31:33
% EndTime: 2020-01-03 11:31:45
% DurationCPUTime: 4.16s
% Computational Cost: add. (5009->363), mult. (12548->501), div. (0->0), fcn. (9488->14), ass. (0->199)
t176 = cos(pkin(8));
t235 = qJD(1) * t176;
t282 = qJD(4) - t235;
t178 = sin(qJ(5));
t181 = cos(qJ(5));
t228 = qJD(5) * t181;
t229 = qJD(5) * t178;
t174 = sin(pkin(8));
t179 = sin(qJ(4));
t175 = cos(pkin(9));
t182 = cos(qJ(4));
t245 = t175 * t182;
t220 = t174 * t245;
t173 = sin(pkin(9));
t225 = qJDD(1) * t173;
t123 = t173 * t182 + t175 * t179;
t275 = t123 * qJD(4);
t56 = -qJDD(1) * t220 + (qJD(1) * t275 + t179 * t225) * t174;
t236 = qJD(1) * t174;
t218 = t173 * t236;
t206 = t179 * t218;
t230 = qJD(4) * t182;
t217 = t175 * t230;
t57 = -qJD(4) * t206 + (qJD(1) * t217 + t123 * qJDD(1)) * t174;
t196 = qJD(1) * t123;
t94 = t174 * t196;
t97 = qJD(1) * t220 - t206;
t16 = t178 * t57 + t181 * t56 + t94 * t228 + t97 * t229;
t141 = -qJD(5) - t282;
t43 = t178 * t97 + t181 * t94;
t257 = t141 * t43;
t281 = -t16 - t257;
t249 = qJDD(1) * pkin(1);
t154 = qJDD(2) - t249;
t180 = sin(qJ(1));
t183 = cos(qJ(1));
t239 = -g(2) * t183 - g(3) * t180;
t280 = t239 - t154;
t262 = t43 ^ 2;
t202 = -t178 * t94 + t181 * t97;
t263 = t202 ^ 2;
t279 = -t262 + t263;
t261 = t43 * t202;
t122 = -t173 * t179 + t245;
t251 = t282 * t122;
t250 = t176 * t196 - t275;
t17 = t202 * qJD(5) - t178 * t56 + t181 * t57;
t254 = t202 * t141;
t278 = -t17 - t254;
t128 = -pkin(2) * t176 - qJ(3) * t174 - pkin(1);
t113 = t128 * qJD(1) + qJD(2);
t103 = t175 * t113;
t193 = -t175 * t174 * pkin(6) + (-qJ(2) * t173 - pkin(3)) * t176;
t58 = t193 * qJD(1) + t103;
t219 = qJ(2) * t235;
t72 = t173 * t113 + t175 * t219;
t63 = -pkin(6) * t218 + t72;
t30 = -t179 * t63 + t182 * t58;
t22 = -pkin(7) * t97 + t30;
t20 = pkin(4) * t282 + t22;
t31 = t179 * t58 + t182 * t63;
t23 = -pkin(7) * t94 + t31;
t226 = qJD(1) * qJD(2);
t215 = t176 * t226;
t232 = qJD(3) * t174;
t83 = -qJD(1) * t232 + t128 * qJDD(1) + qJDD(2);
t76 = t175 * t83;
t39 = t193 * qJDD(1) - t173 * t215 + t76;
t223 = t174 * qJDD(1);
t214 = t173 * t223;
t222 = t176 * qJDD(1);
t212 = t175 * t222;
t55 = qJ(2) * t212 + t173 * t83 + t175 * t215;
t41 = -pkin(6) * t214 + t55;
t13 = -t31 * qJD(4) - t179 * t41 + t182 * t39;
t144 = -qJDD(4) + t222;
t6 = -t144 * pkin(4) + t56 * pkin(7) + t13;
t231 = qJD(4) * t179;
t12 = t179 * t39 + t182 * t41 + t58 * t230 - t63 * t231;
t7 = -pkin(7) * t57 + t12;
t1 = (qJD(5) * t20 + t7) * t181 + t178 * t6 - t23 * t229;
t172 = pkin(9) + qJ(4);
t159 = qJ(5) + t172;
t152 = cos(t159);
t268 = g(1) * t174;
t142 = qJ(2) * t236 + qJD(3);
t114 = pkin(3) * t218 + t142;
t61 = pkin(4) * t94 + t114;
t151 = sin(t159);
t241 = t183 * t151;
t244 = t176 * t180;
t90 = t152 * t244 - t241;
t243 = t176 * t183;
t92 = t151 * t180 + t152 * t243;
t277 = g(2) * t90 - g(3) * t92 + t152 * t268 + t43 * t61 - t1;
t227 = qJ(2) * qJDD(1);
t156 = sin(t172);
t157 = cos(t172);
t104 = -t156 * t244 - t157 * t183;
t106 = t156 * t243 - t157 * t180;
t274 = -g(2) * t104 - g(3) * t106 + t156 * t268;
t255 = t181 * t23;
t9 = t178 * t20 + t255;
t2 = -t9 * qJD(5) - t178 * t7 + t181 * t6;
t89 = -t151 * t244 - t152 * t183;
t91 = -t152 * t180 + t176 * t241;
t273 = -g(2) * t89 - g(3) * t91 + t151 * t268 - t202 * t61 + t2;
t271 = t97 ^ 2;
t270 = pkin(4) * t57;
t269 = pkin(3) * t173;
t266 = g(2) * t180;
t264 = g(3) * t183;
t260 = t97 * t94;
t177 = -pkin(6) - qJ(3);
t153 = t175 * pkin(3) + pkin(2);
t69 = t122 * t181 - t123 * t178;
t259 = t69 * qJD(5) + t250 * t178 + t251 * t181;
t70 = t122 * t178 + t123 * t181;
t258 = -t70 * qJD(5) - t251 * t178 + t250 * t181;
t121 = t175 * t128;
t66 = t121 + t193;
t248 = t173 * t174;
t87 = t175 * t176 * qJ(2) + t173 * t128;
t73 = -pkin(6) * t248 + t87;
t35 = t179 * t66 + t182 * t73;
t256 = t178 * t23;
t253 = t94 * t282;
t252 = t97 * t282;
t247 = t173 * t176;
t246 = t174 * t180;
t184 = qJD(1) ^ 2;
t242 = t176 * t184;
t124 = pkin(3) * t248 + t174 * qJ(2);
t240 = t183 * pkin(1) + t180 * qJ(2);
t168 = t173 ^ 2;
t170 = t175 ^ 2;
t238 = -t168 - t170;
t169 = t174 ^ 2;
t171 = t176 ^ 2;
t237 = t169 + t171;
t234 = qJD(2) * t174;
t233 = qJD(2) * t176;
t224 = t169 * qJDD(1);
t119 = qJ(2) * t223 + t174 * t226 + qJDD(3);
t34 = -t179 * t73 + t182 * t66;
t209 = t237 * t184;
t208 = 0.2e1 * t174 * t222;
t207 = 0.2e1 * t237;
t139 = pkin(3) * t214;
t84 = t139 + t119;
t203 = -t264 + t266;
t112 = t122 * t174;
t28 = -pkin(4) * t176 - pkin(7) * t112 + t34;
t111 = t123 * t174;
t29 = -pkin(7) * t111 + t35;
t14 = -t178 * t29 + t181 * t28;
t15 = t178 * t28 + t181 * t29;
t60 = -t111 * t178 + t112 * t181;
t201 = g(1) * t176 + t174 * t264 + t119;
t200 = t226 + t227;
t199 = t239 * t174;
t115 = -t173 * t233 - t175 * t232;
t54 = -t200 * t247 + t76;
t86 = -qJ(2) * t247 + t121;
t198 = -t115 * qJD(1) - t86 * qJDD(1) - t54;
t116 = -t173 * t232 + t175 * t233;
t197 = t116 * qJD(1) + t87 * qJDD(1) + t55;
t24 = t179 * t115 + t182 * t116 + t66 * t230 - t73 * t231;
t194 = t249 + t280;
t192 = t207 * t226 + t264;
t191 = -g(2) * t246 + t139 + t201;
t25 = -t35 * qJD(4) + t182 * t115 - t179 * t116;
t188 = t119 * t174 + t200 * t169 - t203;
t167 = -pkin(7) + t177;
t163 = t180 * pkin(1);
t155 = t171 * qJDD(1);
t138 = -qJDD(5) + t144;
t127 = pkin(4) * t156 + t269;
t125 = pkin(4) * t157 + t153;
t107 = t156 * t180 + t157 * t243;
t105 = -t156 * t183 + t157 * t244;
t101 = t174 * t217 - t231 * t248;
t100 = t174 * t275;
t88 = t94 ^ 2;
t74 = pkin(4) * t101 + t234;
t71 = -t173 * t219 + t103;
t68 = pkin(4) * t111 + t124;
t59 = t181 * t111 + t112 * t178;
t36 = t84 + t270;
t27 = t60 * qJD(5) - t178 * t100 + t181 * t101;
t26 = t181 * t100 + t178 * t101 + t111 * t228 + t112 * t229;
t19 = t100 * pkin(7) + t25;
t18 = -pkin(7) * t101 + t24;
t11 = t181 * t22 - t256;
t10 = -t178 * t22 - t255;
t8 = t181 * t20 - t256;
t4 = -t15 * qJD(5) - t178 * t18 + t181 * t19;
t3 = t14 * qJD(5) + t178 * t19 + t181 * t18;
t5 = [0, 0, 0, 0, 0, qJDD(1), t239, t203, 0, 0, t224, t208, 0, t155, 0, 0, t194 * t176, -t194 * t174, t207 * t227 + t192 - t266, -t154 * pkin(1) - g(2) * t240 - g(3) * t163 + (t237 * t227 + t192) * qJ(2), t170 * t224, -0.2e1 * t175 * t173 * t224, -0.2e1 * t174 * t212, t168 * t224, t173 * t208, t155, (t175 * t239 + t198) * t176 + t188 * t173, (-t173 * t239 + t197) * t176 + t188 * t175, (-t173 * t197 + t175 * t198 + t239) * t174, t55 * t87 + t72 * t116 + t54 * t86 + t71 * t115 - g(2) * (pkin(2) * t243 + t240) - g(3) * (pkin(2) * t244 - t183 * qJ(2) + t163) + (t119 * qJ(2) + qJ(3) * t239 + t142 * qJD(2)) * t174, -t100 * t97 - t112 * t56, t100 * t94 - t101 * t97 + t111 * t56 - t112 * t57, -t100 * t282 - t112 * t144 + t176 * t56, t101 * t94 + t111 * t57, -t101 * t282 + t111 * t144 + t176 * t57, t144 * t176, -g(2) * t107 - g(3) * t105 + t101 * t114 + t111 * t84 + t124 * t57 - t13 * t176 - t144 * t34 + t94 * t234 + t25 * t282, g(2) * t106 - g(3) * t104 - t100 * t114 + t112 * t84 + t12 * t176 - t124 * t56 + t144 * t35 + t97 * t234 - t24 * t282, t100 * t30 - t101 * t31 - t111 * t12 - t112 * t13 - t24 * t94 - t25 * t97 + t34 * t56 - t35 * t57 + t199, t12 * t35 + t31 * t24 + t13 * t34 + t30 * t25 + t84 * t124 + t114 * t234 - g(2) * (t180 * t269 + t240) - g(3) * (t153 * t244 - t177 * t246 + t163) + (-g(2) * (t153 * t176 - t174 * t177) - g(3) * (-qJ(2) - t269)) * t183, -t16 * t60 - t202 * t26, t16 * t59 - t17 * t60 - t202 * t27 + t26 * t43, -t138 * t60 + t141 * t26 + t16 * t176, t17 * t59 + t27 * t43, t138 * t59 + t141 * t27 + t17 * t176, t138 * t176, -g(2) * t92 - g(3) * t90 - t138 * t14 - t141 * t4 + t17 * t68 - t176 * t2 + t27 * t61 + t36 * t59 + t43 * t74, g(2) * t91 - g(3) * t89 + t1 * t176 + t138 * t15 + t141 * t3 - t16 * t68 + t202 * t74 - t26 * t61 + t36 * t60, -t1 * t59 + t14 * t16 - t15 * t17 - t2 * t60 - t202 * t4 + t26 * t8 - t27 * t9 - t3 * t43 + t199, t1 * t15 + t9 * t3 + t2 * t14 + t8 * t4 + t36 * t68 + t61 * t74 - g(2) * (t180 * t127 + t240) - g(3) * (t125 * t244 - t167 * t246 + t163) + (-g(2) * (t125 * t176 - t167 * t174) - g(3) * (-qJ(2) - t127)) * t183; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t222, t223, -t209, -qJ(2) * t209 - t280, 0, 0, 0, 0, 0, 0, -t173 * t209 - t212, t173 * t222 - t175 * t209, t238 * t223, t55 * t173 + t54 * t175 + (-t142 * t174 + (t173 * t71 - t175 * t72) * t176) * qJD(1) - t239, 0, 0, 0, 0, 0, 0, -t122 * t144 - t94 * t236 + t250 * t282, t123 * t144 - t97 * t236 - t251 * t282, t122 * t56 - t123 * t57 - t250 * t97 - t251 * t94, -t114 * t236 + t12 * t123 + t122 * t13 + t250 * t30 + t251 * t31 - t239, 0, 0, 0, 0, 0, 0, -t69 * t138 - t258 * t141 - t43 * t236, t70 * t138 + t259 * t141 - t202 * t236, t69 * t16 - t70 * t17 - t202 * t258 - t259 * t43, t1 * t70 + t2 * t69 - t61 * t236 + t258 * t8 + t259 * t9 - t239; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t175 * t242 + t225) * t174, (qJDD(1) * t175 + t173 * t242) * t174, t238 * t184 * t169, (-t266 + (t173 * t72 + t175 * t71) * qJD(1)) * t174 + t201, 0, 0, 0, 0, 0, 0, t252 + t57, -t56 - t253, -t88 - t271, t30 * t97 + t31 * t94 + t191, 0, 0, 0, 0, 0, 0, t17 - t254, -t16 + t257, -t262 - t263, t202 * t8 + t43 * t9 + t191 + t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t260, -t88 + t271, -t56 + t253, -t260, -t57 + t252, -t144, -t114 * t97 + t282 * t31 + t13 + t274, g(2) * t105 - g(3) * t107 + t114 * t94 + t157 * t268 + t282 * t30 - t12, 0, 0, t261, t279, t281, -t261, t278, -t138, t10 * t141 + (-t138 * t181 + t141 * t229 - t43 * t97) * pkin(4) + t273, -t11 * t141 + (t138 * t178 + t141 * t228 - t202 * t97) * pkin(4) + t277, t10 * t202 + t11 * t43 + t9 * t202 - t8 * t43 + (t16 * t181 - t17 * t178 + (t178 * t202 - t181 * t43) * qJD(5)) * pkin(4), -t8 * t10 - t9 * t11 + (t1 * t178 + t2 * t181 - t61 * t97 + (-t178 * t8 + t181 * t9) * qJD(5) + t274) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t261, t279, t281, -t261, t278, -t138, -t141 * t9 + t273, -t141 * t8 + t277, 0, 0;];
tau_reg = t5;
