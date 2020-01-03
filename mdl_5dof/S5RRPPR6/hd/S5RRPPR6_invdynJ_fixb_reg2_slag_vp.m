% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPPR6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_invdynJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:33:12
% EndTime: 2019-12-31 19:33:26
% DurationCPUTime: 5.95s
% Computational Cost: add. (6353->454), mult. (15099->588), div. (0->0), fcn. (11229->14), ass. (0->218)
t175 = sin(pkin(8));
t182 = cos(qJ(2));
t180 = sin(qJ(2));
t261 = cos(pkin(8));
t225 = t261 * t180;
t143 = t175 * t182 + t225;
t130 = t143 * qJD(1);
t174 = sin(pkin(9));
t176 = cos(pkin(9));
t112 = t174 * qJD(2) + t176 * t130;
t179 = sin(qJ(5));
t288 = cos(qJ(5));
t251 = t174 * t130;
t303 = t176 * qJD(2) - t251;
t200 = t288 * t303;
t55 = -t179 * t112 + t200;
t306 = t55 ^ 2;
t224 = t261 * t182;
t155 = qJD(1) * t224;
t242 = qJD(1) * t180;
t127 = t175 * t242 - t155;
t124 = qJD(5) + t127;
t305 = t55 * t124;
t54 = t288 * t112 + t179 * t303;
t304 = t54 ^ 2;
t234 = t288 * t176;
t201 = -t179 * t174 + t234;
t241 = qJD(5) * t179;
t294 = qJD(5) * t234 - t174 * t241;
t263 = -t201 * t127 - t294;
t144 = t288 * t174 + t179 * t176;
t134 = t144 * qJD(5);
t262 = t144 * t127 + t134;
t181 = sin(qJ(1));
t183 = cos(qJ(1));
t300 = -g(1) * t181 + g(2) * t183;
t171 = qJ(2) + pkin(8);
t166 = sin(t171);
t168 = cos(t171);
t217 = g(1) * t183 + g(2) * t181;
t299 = -g(3) * t168 + t217 * t166;
t275 = qJ(3) + pkin(6);
t227 = qJD(2) * t275;
t125 = t182 * qJD(3) - t180 * t227;
t150 = t275 * t182;
t103 = t125 * qJD(1) + qJDD(1) * t150;
t196 = -t180 * qJD(3) - t182 * t227;
t232 = t275 * t180;
t93 = qJDD(2) * pkin(2) + t196 * qJD(1) - qJDD(1) * t232;
t42 = -t175 * t103 + t261 * t93;
t41 = -qJDD(2) * pkin(3) + qJDD(4) - t42;
t187 = t41 - t299;
t203 = t300 * t168;
t297 = t303 * t130;
t146 = qJD(1) * t232;
t147 = qJD(1) * t150;
t226 = t261 * t147;
t104 = -t175 * t146 + t226;
t295 = t104 * qJD(2) + t299;
t239 = qJD(1) * qJD(2);
t231 = t180 * t239;
t188 = t143 * qJDD(1) - t175 * t231;
t101 = qJD(2) * t155 + t188;
t74 = -t176 * qJDD(2) + t174 * t101;
t75 = t174 * qJDD(2) + t176 * t101;
t17 = -qJD(5) * t200 + t112 * t241 + t179 * t74 - t288 * t75;
t293 = -t17 * t201 - t262 * t54;
t129 = t143 * qJD(2);
t238 = t180 * qJDD(1);
t214 = -qJDD(1) * t224 + t175 * t238;
t100 = qJD(1) * t129 + t214;
t94 = qJDD(5) + t100;
t292 = -t263 * t124 + t144 * t94;
t126 = t127 ^ 2;
t252 = t174 * t100;
t291 = -t176 * t126 - t252;
t282 = g(3) * t166;
t193 = -t217 * t168 - t282;
t278 = t182 * pkin(2);
t164 = pkin(1) + t278;
t149 = -t164 * qJD(1) + qJD(3);
t67 = t127 * pkin(3) - t130 * qJ(4) + t149;
t271 = qJD(2) * pkin(2);
t140 = -t146 + t271;
t99 = t175 * t140 + t226;
t87 = qJD(2) * qJ(4) + t99;
t32 = -t174 * t87 + t176 * t67;
t20 = t127 * pkin(4) - t112 * pkin(7) + t32;
t33 = t174 * t67 + t176 * t87;
t24 = pkin(7) * t303 + t33;
t206 = t179 * t24 - t288 * t20;
t123 = pkin(2) * t231 - t164 * qJDD(1) + qJDD(3);
t27 = t100 * pkin(3) - t101 * qJ(4) - t130 * qJD(4) + t123;
t43 = t261 * t103 + t175 * t93;
t38 = qJDD(2) * qJ(4) + qJD(2) * qJD(4) + t43;
t12 = -t174 * t38 + t176 * t27;
t6 = t100 * pkin(4) - t75 * pkin(7) + t12;
t13 = t174 * t27 + t176 * t38;
t9 = -pkin(7) * t74 + t13;
t1 = -t206 * qJD(5) + t179 * t6 + t288 * t9;
t290 = t130 ^ 2;
t289 = t74 * pkin(4);
t287 = pkin(2) * t180;
t153 = t183 * t164;
t284 = g(2) * t153;
t281 = g(3) * t182;
t280 = t176 * pkin(4);
t279 = t176 * pkin(7);
t277 = t54 * t55;
t158 = t175 * pkin(2) + qJ(4);
t276 = pkin(7) + t158;
t255 = t127 * t176;
t135 = t175 * t147;
t105 = -t261 * t146 - t135;
t77 = pkin(2) * t242 + t130 * pkin(3) + t127 * qJ(4);
t39 = -t174 * t105 + t176 * t77;
t23 = t130 * pkin(4) + pkin(7) * t255 + t39;
t256 = t127 * t174;
t40 = t176 * t105 + t174 * t77;
t31 = pkin(7) * t256 + t40;
t138 = t276 * t174;
t139 = t276 * t176;
t89 = -t288 * t138 - t179 * t139;
t274 = t201 * qJD(4) + t89 * qJD(5) - t179 * t23 - t288 * t31;
t90 = -t179 * t138 + t288 * t139;
t273 = -t144 * qJD(4) - t90 * qJD(5) + t179 * t31 - t288 * t23;
t199 = -t175 * t180 + t224;
t132 = t199 * qJD(2);
t236 = t180 * t271;
t56 = t129 * pkin(3) - t132 * qJ(4) - t143 * qJD(4) + t236;
t79 = t261 * t125 + t175 * t196;
t29 = t174 * t56 + t176 * t79;
t88 = t176 * t100;
t272 = -t174 * t126 + t88;
t270 = t130 * t55;
t267 = t54 * t130;
t266 = t74 * t176;
t265 = t75 * t174;
t264 = t75 * t176;
t110 = t261 * t150 - t175 * t232;
t96 = -pkin(3) * t199 - t143 * qJ(4) - t164;
t46 = t176 * t110 + t174 * t96;
t260 = pkin(6) * qJDD(1);
t259 = t112 * t130;
t258 = t112 * t174;
t257 = t127 * t130;
t254 = t132 * t174;
t253 = t143 * t174;
t170 = pkin(9) + qJ(5);
t165 = sin(t170);
t249 = t181 * t165;
t167 = cos(t170);
t248 = t181 * t167;
t247 = t183 * t165;
t246 = t183 * t167;
t98 = t261 * t140 - t135;
t80 = -qJD(2) * pkin(3) + qJD(4) - t98;
t245 = -qJD(4) + t80;
t172 = t180 ^ 2;
t173 = t182 ^ 2;
t244 = t172 - t173;
t243 = t172 + t173;
t237 = t182 * qJDD(1);
t185 = qJD(1) ^ 2;
t235 = t180 * t185 * t182;
t230 = pkin(4) * t174 + t275;
t28 = -t174 * t79 + t176 * t56;
t228 = t179 * t75 + t288 * t74;
t45 = -t174 * t110 + t176 * t96;
t78 = t175 * t125 - t261 * t196;
t109 = t175 * t150 + t275 * t225;
t18 = qJD(5) * t54 + t228;
t222 = -t144 * t18 - t263 * t55;
t220 = t182 * t231;
t219 = t300 * t166;
t218 = -t262 * t124 + t201 * t94;
t163 = -t261 * pkin(2) - pkin(3);
t215 = t168 * pkin(3) + t166 * qJ(4);
t213 = -t12 * t176 - t13 * t174;
t212 = -t32 * t174 + t33 * t176;
t211 = -t100 * t199 + t127 * t129;
t210 = t143 * t100 + t132 * t127;
t162 = pkin(3) + t280;
t178 = -pkin(7) - qJ(4);
t209 = t168 * t162 - t166 * t178;
t208 = t176 * t303;
t30 = -pkin(4) * t199 - t143 * t279 + t45;
t34 = -pkin(7) * t253 + t46;
t14 = -t179 * t34 + t288 * t30;
t8 = t179 * t20 + t288 * t24;
t15 = t179 * t30 + t288 * t34;
t202 = -0.2e1 * pkin(1) * t239 - pkin(6) * qJDD(2);
t195 = t208 - t258;
t184 = qJD(2) ^ 2;
t191 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t184 - t300;
t190 = pkin(1) * t185 + t217 - t260;
t189 = t80 * t132 + t41 * t143 - t217;
t2 = -t8 * qJD(5) - t179 * t9 + t288 * t6;
t148 = t163 - t280;
t118 = t168 * t246 + t249;
t117 = -t168 * t247 + t248;
t116 = -t168 * t248 + t247;
t115 = t168 * t249 + t246;
t82 = t201 * t143;
t81 = t144 * t143;
t73 = pkin(4) * t253 + t109;
t63 = t174 * t74;
t59 = -pkin(4) * t256 + t104;
t48 = pkin(4) * t254 + t78;
t47 = -pkin(4) * t303 + t80;
t37 = t144 * t132 + t143 * t294;
t36 = -t132 * t201 + t143 * t134;
t22 = t41 + t289;
t21 = -pkin(7) * t254 + t29;
t19 = t129 * pkin(4) - t132 * t279 + t28;
t4 = -t15 * qJD(5) - t179 * t21 + t288 * t19;
t3 = t14 * qJD(5) + t179 * t19 + t288 * t21;
t5 = [0, 0, 0, 0, 0, qJDD(1), -t300, t217, 0, 0, t172 * qJDD(1) + 0.2e1 * t220, 0.2e1 * t180 * t237 - 0.2e1 * t244 * t239, qJDD(2) * t180 + t184 * t182, t173 * qJDD(1) - 0.2e1 * t220, qJDD(2) * t182 - t184 * t180, 0, t202 * t180 + t191 * t182, -t191 * t180 + t202 * t182, 0.2e1 * t243 * t260 - t217, -g(1) * (-t181 * pkin(1) + t183 * pkin(6)) - g(2) * (t183 * pkin(1) + t181 * pkin(6)) + (t243 * pkin(6) ^ 2 + pkin(1) ^ 2) * qJDD(1), t101 * t143 + t130 * t132, t101 * t199 - t130 * t129 - t210, t132 * qJD(2) + t143 * qJDD(2), t211, -t129 * qJD(2) + qJDD(2) * t199, 0, -t109 * qJDD(2) - t164 * t100 - t123 * t199 + t149 * t129 - t203 + (t127 * t287 - t78) * qJD(2), -t110 * qJDD(2) - t164 * t101 + t123 * t143 + t149 * t132 + (t130 * t287 - t79) * qJD(2) + t219, -t110 * t100 + t109 * t101 - t79 * t127 - t99 * t129 + t78 * t130 - t98 * t132 - t42 * t143 + t199 * t43 - t217, t43 * t110 + t99 * t79 - t42 * t109 - t98 * t78 - t123 * t164 + t149 * t236 - g(1) * (-t181 * t164 + t183 * t275) - g(2) * (t181 * t275 + t153), (t112 * t132 + t75 * t143) * t176, (-t265 - t266) * t143 + t195 * t132, t112 * t129 + t210 * t176 - t199 * t75, (-t132 * t303 + t74 * t143) * t174, t129 * t303 - t210 * t174 + t199 * t74, t211, t45 * t100 + t109 * t74 - t12 * t199 + t28 * t127 + t32 * t129 + t189 * t174 - t176 * t203 - t303 * t78, -t46 * t100 + t109 * t75 + t78 * t112 - t29 * t127 - t33 * t129 + t13 * t199 + t174 * t203 + t189 * t176, t29 * t303 - t46 * t74 - t28 * t112 - t45 * t75 + t213 * t143 + (-t33 * t174 - t32 * t176) * t132 - t219, -t284 + t41 * t109 + t12 * t45 + t13 * t46 + t32 * t28 + t33 * t29 + t80 * t78 + (-g(1) * t275 - g(2) * t215) * t183 + (-g(1) * (-t164 - t215) - g(2) * t275) * t181, -t17 * t82 - t36 * t54, t17 * t81 - t18 * t82 - t36 * t55 - t37 * t54, -t36 * t124 + t54 * t129 + t17 * t199 + t82 * t94, t18 * t81 - t37 * t55, -t37 * t124 + t129 * t55 + t18 * t199 - t81 * t94, t124 * t129 - t199 * t94, -g(1) * t116 - g(2) * t118 + t4 * t124 - t129 * t206 + t14 * t94 + t73 * t18 - t199 * t2 + t22 * t81 + t47 * t37 - t48 * t55, -g(1) * t115 - g(2) * t117 + t1 * t199 - t3 * t124 - t8 * t129 - t15 * t94 - t73 * t17 + t22 * t82 - t47 * t36 + t48 * t54, -t1 * t81 + t14 * t17 - t15 * t18 - t2 * t82 - t206 * t36 + t3 * t55 - t8 * t37 - t4 * t54 - t219, -t284 + t1 * t15 + t2 * t14 + t22 * t73 + t8 * t3 - t206 * t4 + t47 * t48 + (-g(1) * t230 - g(2) * t209) * t183 + (-g(1) * (-t164 - t209) - g(2) * t230) * t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t235, t244 * t185, t238, t235, t237, qJDD(2), t180 * t190 - t281, g(3) * t180 + t182 * t190, 0, 0, t257, -t126 + t290, (t155 + t127) * qJD(2) + t188, -t257, -t214, qJDD(2), -t149 * t130 + (t261 * qJDD(2) - t127 * t242) * pkin(2) + t42 + t295, t105 * qJD(2) + t149 * t127 + (-qJDD(2) * t175 - t130 * t242) * pkin(2) - t43 - t193, (-t104 + t99) * t130 + (t105 - t98) * t127 + (-t100 * t175 - t261 * t101) * pkin(2), t98 * t104 - t99 * t105 + (t261 * t42 - t281 + t175 * t43 + (-qJD(1) * t149 + t217) * t180) * pkin(2), t112 * t255 + t265, t127 * t195 + t264 - t63, -t259 - t291, -t256 * t303 - t266, t272 - t297, -t257, -t158 * t252 - t104 * t251 - t32 * t130 + t163 * t74 + (t245 * t174 - t39) * t127 + (-t41 + t295) * t176, -t158 * t88 - t104 * t112 + t33 * t130 + t163 * t75 + (t245 * t176 + t40) * t127 + t187 * t174, t39 * t112 + t40 * t251 + (-qJD(4) * t251 - t32 * t127 - t158 * t74 + t13 + (qJD(4) * t176 - t40) * qJD(2)) * t176 + (qJD(4) * t112 - t33 * t127 + t158 * t75 - t12) * t174 + t193, t41 * t163 - t33 * t40 - t32 * t39 - t80 * t104 - g(3) * (t215 + t278) + (-t12 * t174 + t13 * t176) * t158 + t212 * qJD(4) + t217 * (pkin(3) * t166 - qJ(4) * t168 + t287), -t17 * t144 - t263 * t54, t222 + t293, -t267 + t292, -t18 * t201 - t262 * t55, t218 - t270, -t124 * t130, t273 * t124 + t206 * t130 + t148 * t18 + t299 * t167 - t22 * t201 + t262 * t47 + t55 * t59 + t89 * t94, -t274 * t124 + t8 * t130 + t22 * t144 - t148 * t17 - t165 * t299 - t263 * t47 - t59 * t54 - t90 * t94, t1 * t201 - t2 * t144 + t89 * t17 - t90 * t18 - t206 * t263 - t262 * t8 - t273 * t54 + t274 * t55 + t193, t1 * t90 + t2 * t89 + t22 * t148 - t47 * t59 - g(3) * (t209 + t278) + t274 * t8 - t273 * t206 + t217 * (t162 * t166 + t168 * t178 + t287); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t130 * qJD(2) + t214, (t155 - t127) * qJD(2) + t188, -t126 - t290, t99 * t127 + t98 * t130 + t123 + t300, 0, 0, 0, 0, 0, 0, t272 + t297, -t259 + t291, -t264 - t63 + (t208 + t258) * t127, t212 * t127 - t80 * t130 - t213 + t300, 0, 0, 0, 0, 0, 0, t218 + t270, -t267 - t292, t222 - t293, t1 * t144 - t47 * t130 + t2 * t201 + t206 * t262 - t263 * t8 + t300; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112 * t127 + t74, t127 * t303 + t75, -t112 ^ 2 - t303 ^ 2, t112 * t32 - t303 * t33 + t187, 0, 0, 0, 0, 0, 0, t54 * t124 + t18, -t17 + t305, -t304 - t306, -t206 * t54 - t8 * t55 + t187 + t289; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t277, t304 - t306, -t17 - t305, t277, -t228 + (-qJD(5) + t124) * t54, t94, -g(1) * t117 + g(2) * t115 + t8 * t124 + t165 * t282 - t47 * t54 + t2, g(1) * t118 - g(2) * t116 - t124 * t206 + t167 * t282 - t47 * t55 - t1, 0, 0;];
tau_reg = t5;
