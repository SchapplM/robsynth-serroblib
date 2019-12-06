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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 17:45:05
% EndTime: 2019-12-05 17:45:15
% DurationCPUTime: 3.84s
% Computational Cost: add. (5009->365), mult. (12548->498), div. (0->0), fcn. (9488->14), ass. (0->199)
t174 = cos(pkin(8));
t234 = qJD(1) * t174;
t281 = qJD(4) - t234;
t176 = sin(qJ(5));
t179 = cos(qJ(5));
t227 = qJD(5) * t179;
t228 = qJD(5) * t176;
t172 = sin(pkin(8));
t177 = sin(qJ(4));
t173 = cos(pkin(9));
t180 = cos(qJ(4));
t243 = t173 * t180;
t219 = t172 * t243;
t171 = sin(pkin(9));
t224 = qJDD(1) * t171;
t123 = t171 * t180 + t173 * t177;
t275 = t123 * qJD(4);
t56 = -qJDD(1) * t219 + (qJD(1) * t275 + t177 * t224) * t172;
t235 = qJD(1) * t172;
t217 = t171 * t235;
t205 = t177 * t217;
t229 = qJD(4) * t180;
t216 = t173 * t229;
t57 = -qJD(4) * t205 + (qJD(1) * t216 + t123 * qJDD(1)) * t172;
t194 = qJD(1) * t123;
t94 = t172 * t194;
t97 = qJD(1) * t219 - t205;
t16 = t176 * t57 + t179 * t56 + t94 * t227 + t97 * t228;
t141 = -qJD(5) - t281;
t43 = t176 * t97 + t179 * t94;
t256 = t141 * t43;
t280 = -t16 - t256;
t261 = t43 ^ 2;
t201 = -t176 * t94 + t179 * t97;
t262 = t201 ^ 2;
t279 = -t261 + t262;
t260 = t43 * t201;
t122 = -t171 * t177 + t243;
t250 = t281 * t122;
t249 = t174 * t194 - t275;
t17 = t201 * qJD(5) - t176 * t56 + t179 * t57;
t253 = t201 * t141;
t278 = -t17 - t253;
t128 = -pkin(2) * t174 - qJ(3) * t172 - pkin(1);
t113 = t128 * qJD(1) + qJD(2);
t103 = t173 * t113;
t191 = -t173 * t172 * pkin(6) + (-qJ(2) * t171 - pkin(3)) * t174;
t58 = t191 * qJD(1) + t103;
t218 = qJ(2) * t234;
t72 = t171 * t113 + t173 * t218;
t63 = -pkin(6) * t217 + t72;
t30 = -t177 * t63 + t180 * t58;
t22 = -pkin(7) * t97 + t30;
t20 = pkin(4) * t281 + t22;
t31 = t177 * t58 + t180 * t63;
t23 = -pkin(7) * t94 + t31;
t225 = qJD(1) * qJD(2);
t214 = t174 * t225;
t231 = qJD(3) * t172;
t83 = -qJD(1) * t231 + t128 * qJDD(1) + qJDD(2);
t76 = t173 * t83;
t39 = t191 * qJDD(1) - t171 * t214 + t76;
t222 = t172 * qJDD(1);
t213 = t171 * t222;
t221 = t174 * qJDD(1);
t211 = t173 * t221;
t55 = qJ(2) * t211 + t171 * t83 + t173 * t214;
t41 = -pkin(6) * t213 + t55;
t13 = -t31 * qJD(4) - t177 * t41 + t180 * t39;
t144 = -qJDD(4) + t221;
t6 = -t144 * pkin(4) + t56 * pkin(7) + t13;
t230 = qJD(4) * t177;
t12 = t177 * t39 + t180 * t41 + t58 * t229 - t63 * t230;
t7 = -pkin(7) * t57 + t12;
t1 = (qJD(5) * t20 + t7) * t179 + t176 * t6 - t23 * t228;
t170 = pkin(9) + qJ(4);
t161 = qJ(5) + t170;
t154 = cos(t161);
t268 = g(1) * t172;
t142 = qJ(2) * t235 + qJD(3);
t114 = pkin(3) * t217 + t142;
t61 = pkin(4) * t94 + t114;
t153 = sin(t161);
t181 = cos(qJ(1));
t178 = sin(qJ(1));
t239 = t178 * t154;
t90 = -t153 * t181 + t174 * t239;
t241 = t174 * t181;
t92 = -t153 * t178 - t154 * t241;
t277 = -g(2) * t90 - g(3) * t92 + t154 * t268 + t43 * t61 - t1;
t226 = qJ(2) * qJDD(1);
t158 = sin(t170);
t159 = cos(t170);
t242 = t174 * t178;
t104 = t158 * t242 + t159 * t181;
t106 = t158 * t241 - t159 * t178;
t274 = -g(2) * t104 + g(3) * t106 + t158 * t268;
t254 = t179 * t23;
t9 = t176 * t20 + t254;
t2 = -t9 * qJD(5) - t176 * t7 + t179 * t6;
t89 = t153 * t242 + t154 * t181;
t91 = t153 * t241 - t239;
t273 = -g(2) * t89 + g(3) * t91 + t153 * t268 - t201 * t61 + t2;
t271 = t97 ^ 2;
t270 = pkin(4) * t57;
t269 = pkin(3) * t171;
t266 = g(2) * t178;
t163 = t181 * qJ(2);
t264 = g(3) * t163;
t263 = g(3) * t181;
t259 = t97 * t94;
t175 = -pkin(6) - qJ(3);
t155 = t173 * pkin(3) + pkin(2);
t69 = t122 * t179 - t123 * t176;
t258 = t69 * qJD(5) + t249 * t176 + t250 * t179;
t70 = t122 * t176 + t123 * t179;
t257 = -t70 * qJD(5) - t250 * t176 + t249 * t179;
t121 = t173 * t128;
t66 = t121 + t191;
t247 = t171 * t172;
t87 = t173 * t174 * qJ(2) + t171 * t128;
t73 = -pkin(6) * t247 + t87;
t35 = t177 * t66 + t180 * t73;
t255 = t176 * t23;
t252 = t94 * t281;
t251 = t97 * t281;
t248 = qJDD(1) * pkin(1);
t246 = t171 * t174;
t245 = t172 * t178;
t244 = t172 * t181;
t182 = qJD(1) ^ 2;
t240 = t174 * t182;
t238 = g(2) * t244 + g(3) * t245;
t124 = pkin(3) * t247 + t172 * qJ(2);
t166 = t171 ^ 2;
t168 = t173 ^ 2;
t237 = -t166 - t168;
t167 = t172 ^ 2;
t169 = t174 ^ 2;
t236 = t167 + t169;
t233 = qJD(2) * t172;
t232 = qJD(2) * t174;
t223 = t167 * qJDD(1);
t119 = qJ(2) * t222 + t172 * t225 + qJDD(3);
t34 = -t177 * t73 + t180 * t66;
t208 = t236 * t182;
t207 = 0.2e1 * t172 * t221;
t206 = 0.2e1 * t236;
t139 = pkin(3) * t213;
t84 = t139 + t119;
t203 = g(2) * t181 + g(3) * t178;
t202 = -t263 + t266;
t112 = t122 * t172;
t28 = -pkin(4) * t174 - pkin(7) * t112 + t34;
t111 = t123 * t172;
t29 = -pkin(7) * t111 + t35;
t14 = -t176 * t29 + t179 * t28;
t15 = t176 * t28 + t179 * t29;
t60 = -t111 * t176 + t112 * t179;
t200 = g(1) * t174 + g(2) * t245 + t119;
t199 = t225 + t226;
t198 = (pkin(4) * t159 + t155) * t174 - (-pkin(7) + t175) * t172 + pkin(1);
t197 = t155 * t174 - t172 * t175 + pkin(1);
t115 = -t171 * t232 - t173 * t231;
t54 = -t199 * t246 + t76;
t86 = -qJ(2) * t246 + t121;
t196 = -t115 * qJD(1) - t86 * qJDD(1) - t54;
t116 = -t171 * t231 + t173 * t232;
t195 = t116 * qJD(1) + t87 * qJDD(1) + t55;
t24 = t177 * t115 + t180 * t116 + t66 * t229 - t73 * t230;
t192 = -t203 - t248;
t190 = t206 * t225 + t266;
t189 = -g(3) * t244 + t139 + t200;
t25 = -t35 * qJD(4) + t180 * t115 - t177 * t116;
t186 = t119 * t172 + t199 * t167 + t202;
t157 = t169 * qJDD(1);
t156 = qJDD(2) - t248;
t138 = -qJDD(5) + t144;
t127 = pkin(4) * t158 + t269;
t107 = -t158 * t178 - t159 * t241;
t105 = -t158 * t181 + t159 * t242;
t101 = t172 * t216 - t230 * t247;
t100 = t172 * t275;
t88 = t94 ^ 2;
t74 = pkin(4) * t101 + t233;
t71 = -t171 * t218 + t103;
t68 = pkin(4) * t111 + t124;
t59 = t179 * t111 + t112 * t176;
t36 = t84 + t270;
t27 = t60 * qJD(5) - t176 * t100 + t179 * t101;
t26 = t179 * t100 + t176 * t101 + t111 * t227 + t112 * t228;
t19 = t100 * pkin(7) + t25;
t18 = -pkin(7) * t101 + t24;
t11 = t179 * t22 - t255;
t10 = -t176 * t22 - t254;
t8 = t179 * t20 - t255;
t4 = -t15 * qJD(5) - t176 * t18 + t179 * t19;
t3 = t14 * qJD(5) + t176 * t19 + t179 * t18;
t5 = [0, 0, 0, 0, 0, qJDD(1), t203, -t202, 0, 0, t223, t207, 0, t157, 0, 0, (-t156 - t192) * t174, (t156 - t248) * t172 - t238, t206 * t226 + t190 - t263, -t264 + (-t156 + t203) * pkin(1) + (t236 * t226 + t190) * qJ(2), t168 * t223, -0.2e1 * t173 * t171 * t223, -0.2e1 * t172 * t211, t166 * t223, t171 * t207, t157, (t203 * t173 + t196) * t174 + t186 * t171, (-t203 * t171 + t195) * t174 + t186 * t173, (-t195 * t171 + t196 * t173) * t172 + t238, t55 * t87 + t72 * t116 + t54 * t86 + t71 * t115 - g(2) * (-t181 * pkin(1) - pkin(2) * t241 - t178 * qJ(2)) - g(3) * (-pkin(1) * t178 - pkin(2) * t242 + t163) + (t119 * qJ(2) + t203 * qJ(3) + t142 * qJD(2)) * t172, -t100 * t97 - t112 * t56, t100 * t94 - t101 * t97 + t111 * t56 - t112 * t57, -t100 * t281 - t112 * t144 + t174 * t56, t101 * t94 + t111 * t57, -t101 * t281 + t111 * t144 + t174 * t57, t144 * t174, -g(2) * t107 + g(3) * t105 + t101 * t114 + t111 * t84 + t124 * t57 - t13 * t174 - t144 * t34 + t94 * t233 + t25 * t281, -g(2) * t106 - g(3) * t104 - t100 * t114 + t112 * t84 + t12 * t174 - t124 * t56 + t144 * t35 + t97 * t233 - t24 * t281, t100 * t30 - t101 * t31 - t111 * t12 - t112 * t13 - t24 * t94 - t25 * t97 + t34 * t56 - t35 * t57 + t238, t114 * t233 - t264 + t12 * t35 + t84 * t124 + t13 * t34 + t31 * t24 + t30 * t25 + (g(2) * t197 - g(3) * t269) * t181 + (-g(2) * (-qJ(2) - t269) + g(3) * t197) * t178, -t16 * t60 - t201 * t26, t16 * t59 - t17 * t60 - t201 * t27 + t26 * t43, -t138 * t60 + t141 * t26 + t16 * t174, t17 * t59 + t27 * t43, t138 * t59 + t141 * t27 + t17 * t174, t138 * t174, -g(2) * t92 + g(3) * t90 - t138 * t14 - t141 * t4 + t17 * t68 - t174 * t2 + t27 * t61 + t36 * t59 + t43 * t74, -g(2) * t91 - g(3) * t89 + t1 * t174 + t138 * t15 + t141 * t3 - t16 * t68 + t201 * t74 - t26 * t61 + t36 * t60, -t1 * t59 + t14 * t16 - t15 * t17 - t2 * t60 - t201 * t4 + t26 * t8 - t27 * t9 - t3 * t43 + t238, -t264 + t1 * t15 + t2 * t14 + t9 * t3 + t36 * t68 + t8 * t4 + t61 * t74 + (g(2) * t198 - g(3) * t127) * t181 + (-g(2) * (-qJ(2) - t127) + g(3) * t198) * t178; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t221, t222, -t208, -qJ(2) * t208 + qJDD(2) + t192, 0, 0, 0, 0, 0, 0, -t171 * t208 - t211, t171 * t221 - t173 * t208, t237 * t222, t55 * t171 + t54 * t173 + (-t142 * t172 + (t171 * t71 - t173 * t72) * t174) * qJD(1) - t203, 0, 0, 0, 0, 0, 0, -t122 * t144 - t94 * t235 + t249 * t281, t123 * t144 - t97 * t235 - t250 * t281, t122 * t56 - t123 * t57 - t249 * t97 - t250 * t94, -t114 * t235 + t12 * t123 + t122 * t13 + t249 * t30 + t250 * t31 - t203, 0, 0, 0, 0, 0, 0, -t69 * t138 - t257 * t141 - t43 * t235, t70 * t138 + t258 * t141 - t201 * t235, t69 * t16 - t70 * t17 - t201 * t257 - t258 * t43, t1 * t70 + t2 * t69 - t235 * t61 + t257 * t8 + t258 * t9 - t203; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t173 * t240 + t224) * t172, (qJDD(1) * t173 + t171 * t240) * t172, t237 * t182 * t167, (-t263 + (t171 * t72 + t173 * t71) * qJD(1)) * t172 + t200, 0, 0, 0, 0, 0, 0, t251 + t57, -t56 - t252, -t88 - t271, t30 * t97 + t31 * t94 + t189, 0, 0, 0, 0, 0, 0, t17 - t253, -t16 + t256, -t261 - t262, t201 * t8 + t43 * t9 + t189 + t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t259, -t88 + t271, -t56 + t252, -t259, -t57 + t251, -t144, -t114 * t97 + t281 * t31 + t13 + t274, -g(2) * t105 - g(3) * t107 + t114 * t94 + t159 * t268 + t281 * t30 - t12, 0, 0, t260, t279, t280, -t260, t278, -t138, t10 * t141 + (-t138 * t179 + t141 * t228 - t43 * t97) * pkin(4) + t273, -t11 * t141 + (t138 * t176 + t141 * t227 - t201 * t97) * pkin(4) + t277, t10 * t201 + t11 * t43 + t201 * t9 - t43 * t8 + (t16 * t179 - t17 * t176 + (t176 * t201 - t179 * t43) * qJD(5)) * pkin(4), -t8 * t10 - t9 * t11 + (t1 * t176 + t2 * t179 - t61 * t97 + (-t176 * t8 + t179 * t9) * qJD(5) + t274) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t260, t279, t280, -t260, t278, -t138, -t141 * t9 + t273, -t141 * t8 + t277, 0, 0;];
tau_reg = t5;
