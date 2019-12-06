% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRPR5
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
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:28:04
% EndTime: 2019-12-05 16:28:15
% DurationCPUTime: 4.73s
% Computational Cost: add. (4003->452), mult. (9663->637), div. (0->0), fcn. (7603->14), ass. (0->220)
t149 = sin(qJ(2));
t144 = sin(pkin(5));
t227 = qJD(1) * t144;
t203 = t149 * t227;
t148 = sin(qJ(3));
t222 = qJD(3) * t148;
t285 = pkin(3) * t222 - t203;
t142 = sin(pkin(10));
t151 = cos(qJ(3));
t264 = qJ(4) + pkin(7);
t192 = qJD(3) * t264;
t161 = -t148 * qJD(4) - t151 * t192;
t249 = cos(pkin(10));
t187 = t249 * t151;
t241 = t142 * t148;
t163 = t187 - t241;
t152 = cos(qJ(2));
t200 = t152 * t227;
t92 = t151 * qJD(4) - t148 * t192;
t261 = t142 * t161 - t163 * t200 + t249 * t92;
t145 = cos(pkin(5));
t226 = qJD(1) * t145;
t218 = qJD(1) * qJD(2);
t197 = t152 * t218;
t247 = qJDD(2) * pkin(7);
t83 = t247 + (qJDD(1) * t149 + t197) * t144;
t284 = qJD(3) * t226 + t83;
t188 = t249 * t148;
t107 = t142 * t151 + t188;
t100 = t107 * qJD(3);
t103 = t163 * qJD(3);
t283 = pkin(4) * t100 - pkin(8) * t103 + t285;
t237 = t144 * t149;
t104 = t145 * t151 - t148 * t237;
t147 = sin(qJ(5));
t150 = cos(qJ(5));
t132 = pkin(3) * t151 + pkin(2);
t235 = t144 * t152;
t193 = qJDD(1) * t235;
t198 = t149 * t218;
t176 = t144 * t198 - t193;
t217 = qJD(2) * qJD(3);
t195 = t148 * t217;
t54 = pkin(3) * t195 - t132 * qJDD(2) + qJDD(4) + t176;
t214 = t148 * qJDD(2);
t175 = -qJDD(2) * t187 + t142 * t214;
t55 = qJD(2) * t100 + t175;
t123 = qJD(2) * t187;
t159 = qJDD(2) * t107 - t142 * t195;
t56 = qJD(3) * t123 + t159;
t16 = pkin(4) * t55 - pkin(8) * t56 + t54;
t199 = t148 * t226;
t111 = qJD(2) * pkin(7) + t203;
t183 = qJ(4) * qJD(2) + t111;
t277 = t183 * t151;
t69 = t199 + t277;
t62 = t249 * t69;
t124 = t151 * t226;
t68 = -t183 * t148 + t124;
t64 = qJD(3) * pkin(3) + t68;
t30 = t142 * t64 + t62;
t24 = qJD(3) * pkin(8) + t30;
t101 = t107 * qJD(2);
t90 = -t132 * qJD(2) + qJD(4) - t200;
t224 = qJD(2) * t148;
t98 = t142 * t224 - t123;
t37 = pkin(4) * t98 - pkin(8) * t101 + t90;
t174 = t147 * t24 - t150 * t37;
t215 = qJDD(1) * t145;
t122 = t151 * t215;
t216 = qJD(2) * qJD(4);
t20 = qJDD(3) * pkin(3) + t122 - qJD(3) * t277 + (-qJ(4) * qJDD(2) - t216 - t284) * t148;
t213 = t151 * qJDD(2);
t210 = t148 * t215 + t284 * t151;
t35 = -t111 * t222 + t210;
t21 = t151 * t216 + (-t195 + t213) * qJ(4) + t35;
t6 = t142 * t20 + t249 * t21;
t4 = qJDD(3) * pkin(8) + t6;
t1 = -t174 * qJD(5) + t147 * t16 + t150 * t4;
t91 = qJD(5) + t98;
t282 = t174 * t91 + t1;
t12 = t147 * t37 + t150 * t24;
t2 = -qJD(5) * t12 - t147 * t4 + t150 * t16;
t281 = t12 * t91 + t2;
t185 = t147 * t91;
t76 = qJD(3) * t147 + t101 * t150;
t280 = t76 * t185;
t221 = qJD(5) * t147;
t51 = qJDD(5) + t55;
t279 = t150 * t51 - t91 * t221;
t139 = qJ(3) + pkin(10);
t133 = sin(t139);
t250 = cos(pkin(9));
t189 = t250 * t152;
t143 = sin(pkin(9));
t239 = t143 * t149;
t94 = -t145 * t189 + t239;
t190 = t250 * t149;
t238 = t143 * t152;
t96 = t145 * t238 + t190;
t180 = g(1) * t96 + g(2) * t94;
t276 = g(3) * t235 - t180;
t278 = t276 * t133;
t153 = qJD(3) ^ 2;
t248 = qJDD(2) * pkin(2);
t82 = t176 - t248;
t173 = t180 - t82;
t275 = -pkin(7) * t153 + t144 * (-g(3) * t152 + t198) + t173 + t248;
t274 = t101 ^ 2;
t53 = -pkin(4) * t163 - pkin(8) * t107 - t132;
t114 = t264 * t151;
t73 = t249 * t114 - t264 * t241;
t27 = -t147 * t73 + t150 * t53;
t273 = qJD(5) * t27 + t283 * t147 + t261 * t150;
t28 = t147 * t53 + t150 * t73;
t272 = -qJD(5) * t28 - t261 * t147 + t283 * t150;
t271 = pkin(3) * t142;
t270 = pkin(3) * t148;
t269 = g(3) * t144;
t219 = t150 * qJD(3);
t74 = t101 * t147 - t219;
t266 = t74 * t98;
t265 = t76 * t74;
t220 = qJD(5) * t150;
t186 = -t150 * qJDD(3) + t147 * t56;
t26 = qJD(5) * t76 + t186;
t263 = -t147 * t26 - t74 * t220;
t262 = -t107 * t200 + t142 * t92 - t249 * t161;
t95 = t145 * t190 + t238;
t260 = -t94 * t132 + t264 * t95;
t97 = -t145 * t239 + t189;
t259 = -t96 * t132 + t264 * t97;
t258 = qJD(2) * pkin(2);
t257 = t101 * t74;
t256 = t101 * t98;
t255 = t142 * t69;
t254 = t147 * t51;
t25 = -qJD(5) * t219 - t147 * qJDD(3) + t101 * t221 - t150 * t56;
t253 = t25 * t147;
t252 = t26 * t150;
t251 = t76 * t101;
t246 = t103 * t147;
t245 = t103 * t150;
t244 = t111 * t148;
t134 = cos(t139);
t243 = t134 * t147;
t242 = t134 * t150;
t240 = t143 * t144;
t236 = t144 * t151;
t233 = t264 * t149;
t232 = t147 * t152;
t231 = t150 * t152;
t230 = qJDD(1) - g(3);
t140 = t148 ^ 2;
t141 = t151 ^ 2;
t229 = t140 - t141;
t228 = t140 + t141;
t112 = -t200 - t258;
t225 = qJD(2) * t112;
t223 = qJD(2) * t149;
t212 = pkin(3) * t224;
t208 = t144 * t232;
t206 = t144 * t231;
t154 = qJD(2) ^ 2;
t205 = t148 * t154 * t151;
t204 = t249 * pkin(3);
t202 = t144 * t223;
t201 = qJD(2) * t235;
t194 = t152 * t217;
t5 = -t142 * t21 + t249 * t20;
t191 = t144 * t250;
t184 = t150 * t91;
t182 = t143 * pkin(3) * t236 - t97 * t270;
t181 = t151 * t195;
t179 = g(1) * t97 + g(2) * t95;
t178 = pkin(4) * t134 + pkin(8) * t133;
t172 = t104 * pkin(3);
t170 = -t98 * t185 + t279;
t169 = qJDD(2) * t152 - t149 * t154;
t105 = t145 * t148 + t149 * t236;
t47 = t142 * t104 + t249 * t105;
t38 = -t147 * t47 - t206;
t168 = -t150 * t47 + t208;
t29 = t249 * t64 - t255;
t79 = t111 * t151 + t199;
t129 = pkin(8) + t271;
t23 = -qJD(3) * pkin(4) - t29;
t166 = -t129 * t51 + t91 * t23;
t57 = -t95 * t133 - t134 * t191;
t59 = -t133 * t97 + t134 * t240;
t85 = -t133 * t237 + t134 * t145;
t165 = -g(1) * t59 - g(2) * t57 - g(3) * t85;
t58 = -t133 * t191 + t134 * t95;
t60 = t133 * t240 + t134 * t97;
t86 = t133 * t145 + t134 * t237;
t164 = -g(1) * t60 - g(2) * t58 - g(3) * t86;
t160 = (-t95 * t148 - t151 * t191) * pkin(3);
t3 = -qJDD(3) * pkin(4) - t5;
t157 = qJD(5) * t129 * t91 - t165 + t3;
t156 = -pkin(7) * qJDD(3) + (t112 + t200 - t258) * qJD(3);
t36 = -t79 * qJD(3) - t148 * t83 + t122;
t78 = t124 - t244;
t155 = -t36 * t148 + t35 * t151 + (-t148 * t79 - t151 * t78) * qJD(3) - t179;
t130 = -t204 - pkin(4);
t109 = t132 * t235;
t93 = t98 ^ 2;
t72 = t114 * t142 + t264 * t188;
t67 = -qJD(3) * t105 - t148 * t201;
t66 = qJD(3) * t104 + t151 * t201;
t46 = -t249 * t104 + t105 * t142;
t42 = pkin(4) * t101 + pkin(8) * t98 + t212;
t34 = t249 * t68 - t255;
t33 = t142 * t67 + t249 * t66;
t32 = t142 * t68 + t62;
t31 = t142 * t66 - t249 * t67;
t14 = t147 * t42 + t150 * t34;
t13 = -t147 * t34 + t150 * t42;
t10 = qJD(5) * t168 - t147 * t33 + t150 * t202;
t9 = qJD(5) * t38 + t147 * t202 + t150 * t33;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t230, 0, 0, 0, 0, 0, 0, t169 * t144, (-qJDD(2) * t149 - t152 * t154) * t144, 0, -g(3) + (t145 ^ 2 + (t149 ^ 2 + t152 ^ 2) * t144 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t67 * qJD(3) + t104 * qJDD(3) + (-t148 * t194 + t151 * t169) * t144, -t66 * qJD(3) - t105 * qJDD(3) + (-t148 * t169 - t151 * t194) * t144, (-t104 * t148 + t105 * t151) * qJDD(2) + (-t148 * t67 + t151 * t66 + (-t104 * t151 - t105 * t148) * qJD(3)) * qJD(2), t104 * t36 + t105 * t35 + t66 * t79 + t67 * t78 - g(3) + (t112 * t223 - t152 * t82) * t144, 0, 0, 0, 0, 0, 0, -t31 * qJD(3) - t46 * qJDD(3) + (-t152 * t55 + t98 * t223) * t144, -t33 * qJD(3) - t47 * qJDD(3) + (t101 * t223 - t152 * t56) * t144, t101 * t31 - t33 * t98 + t46 * t56 - t47 * t55, -t29 * t31 + t30 * t33 - t46 * t5 + t47 * t6 - g(3) + (-t152 * t54 + t90 * t223) * t144, 0, 0, 0, 0, 0, 0, t10 * t91 + t26 * t46 + t31 * t74 + t38 * t51, t168 * t51 - t25 * t46 + t31 * t76 - t9 * t91, -t10 * t76 + t168 * t26 + t25 * t38 - t74 * t9, -t1 * t168 - t10 * t174 + t12 * t9 + t2 * t38 + t23 * t31 + t3 * t46 - g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), -t276 + t193, -t230 * t237 + t179, 0, 0, qJDD(2) * t140 + 0.2e1 * t181, 0.2e1 * t148 * t213 - 0.2e1 * t229 * t217, qJDD(3) * t148 + t151 * t153, qJDD(2) * t141 - 0.2e1 * t181, qJDD(3) * t151 - t148 * t153, 0, t156 * t148 + t275 * t151, -t275 * t148 + t156 * t151, t228 * t247 + (-g(3) * t149 - t228 * t197) * t144 + t155, t173 * pkin(2) + t155 * pkin(7) + (-g(3) * (pkin(2) * t152 + pkin(7) * t149) + (-t112 * t149 + (t148 * t78 - t151 * t79) * t152) * qJD(1)) * t144, t101 * t103 + t107 * t56, -t100 * t101 - t103 * t98 - t107 * t55 + t163 * t56, qJD(3) * t103 + qJDD(3) * t107, t100 * t98 - t163 * t55, -qJD(3) * t100 + qJDD(3) * t163, 0, -t98 * t203 - qJDD(3) * t72 + t100 * t90 - t163 * t54 - t132 * t55 - t276 * t134 + (t98 * t270 - t262) * qJD(3), -t101 * t203 - qJDD(3) * t73 + t103 * t90 + t107 * t54 - t132 * t56 + t278 + (t101 * t270 - t261) * qJD(3), -g(3) * t237 - t100 * t30 + t262 * t101 - t103 * t29 - t107 * t5 + t163 * t6 - t261 * t98 - t55 * t73 + t56 * t72 - t179, t6 * t73 - t5 * t72 - t54 * t132 - g(1) * t259 - g(2) * t260 - g(3) * (t144 * t233 + t109) + t285 * t90 + t261 * t30 - t262 * t29, t76 * t245 + (-t150 * t25 - t76 * t221) * t107, (-t147 * t76 - t150 * t74) * t103 + (t253 - t252 + (t147 * t74 - t150 * t76) * qJD(5)) * t107, t76 * t100 + t279 * t107 + t163 * t25 + t91 * t245, -t263 * t107 + t74 * t246, -t91 * t246 - t74 * t100 + t26 * t163 + (-t91 * t220 - t254) * t107, t100 * t91 - t163 * t51, t27 * t51 - t2 * t163 - t174 * t100 + t72 * t26 + t23 * t246 - g(1) * (t147 * t97 - t96 * t242) - g(2) * (t147 * t95 - t94 * t242) + t272 * t91 + t262 * t74 - (t134 * t231 + t147 * t149) * t269 + (t147 * t3 + t23 * t220) * t107, -t28 * t51 + t1 * t163 - t12 * t100 - t72 * t25 + t23 * t245 - g(1) * (t150 * t97 + t96 * t243) - g(2) * (t150 * t95 + t94 * t243) - t273 * t91 + t262 * t76 - (-t134 * t232 + t149 * t150) * t269 + (t150 * t3 - t23 * t221) * t107, t25 * t27 - t26 * t28 - t272 * t76 - t273 * t74 + (-t12 * t147 + t150 * t174) * t103 - t278 + (-t1 * t147 - t150 * t2 + (-t12 * t150 - t147 * t174) * qJD(5)) * t107, t1 * t28 + t2 * t27 + t3 * t72 - g(1) * (-t178 * t96 + t259) - g(2) * (-t178 * t94 + t260) - g(3) * t109 + t262 * t23 - (t152 * t178 + t233) * t269 + t273 * t12 - t272 * t174; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t205, t229 * t154, t214, t205, t213, qJDD(3), -g(3) * t104 + t122 + (-g(1) * t143 + g(2) * t250) * t236 + (t179 - t83 - t225) * t148, -t151 * t225 - g(1) * (-t148 * t240 - t151 * t97) - g(2) * (t148 * t191 - t95 * t151) + g(3) * t105 + (t78 + t244) * qJD(3) - t210, 0, 0, t256, -t93 + t274, (t123 + t98) * qJD(3) + t159, -t256, -t175, qJDD(3), t32 * qJD(3) - t90 * t101 + (t249 * qJDD(3) - t98 * t224) * pkin(3) + t165 + t5, qJD(3) * t34 + t90 * t98 + (-qJDD(3) * t142 - t101 * t224) * pkin(3) - t164 - t6, (-t29 + t34) * t98 + (t30 - t32) * t101 + (-t142 * t55 - t249 * t56) * pkin(3), -g(1) * t182 - g(2) * t160 - g(3) * t172 + t204 * t5 - t90 * t212 + t6 * t271 + t29 * t32 - t30 * t34, t184 * t76 - t253, (-t25 - t266) * t150 - t280 + t263, t184 * t91 - t251 + t254, t185 * t74 - t252, t170 + t257, -t91 * t101, t101 * t174 - t13 * t91 + t130 * t26 + t147 * t166 - t150 * t157 - t32 * t74, t101 * t12 - t130 * t25 + t14 * t91 + t147 * t157 + t150 * t166 - t32 * t76, t13 * t76 + t14 * t74 + (t174 * t98 - t129 * t26 + t1 + (t129 * t76 + t174) * qJD(5)) * t150 + (-t12 * t98 - t129 * t25 - t2 + (t129 * t74 - t12) * qJD(5)) * t147 + t164, t3 * t130 - t12 * t14 + t174 * t13 - t23 * t32 - g(1) * (pkin(4) * t59 + pkin(8) * t60 + t182) - g(2) * (t57 * pkin(4) + t58 * pkin(8) + t160) - g(3) * (pkin(4) * t85 + pkin(8) * t86 + t172) + (t1 * t150 - t12 * t221 - t2 * t147 + t174 * t220) * t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t101 * qJD(3) + t175, (t123 - t98) * qJD(3) + t159, -t93 - t274, t101 * t29 + t30 * t98 + t276 + t54, 0, 0, 0, 0, 0, 0, t170 - t257, -t91 ^ 2 * t150 - t251 - t254, (t25 - t266) * t150 + t280 + t263, -t101 * t23 + t282 * t147 + t281 * t150 + t276; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t265, -t74 ^ 2 + t76 ^ 2, t74 * t91 - t25, -t265, -t186 + (-qJD(5) + t91) * t76, t51, -t23 * t76 - g(1) * (-t147 * t60 + t150 * t96) - g(2) * (-t147 * t58 + t150 * t94) - g(3) * (-t147 * t86 - t206) + t281, t23 * t74 - g(1) * (-t147 * t96 - t150 * t60) - g(2) * (-t147 * t94 - t150 * t58) - g(3) * (-t150 * t86 + t208) - t282, 0, 0;];
tau_reg = t7;
