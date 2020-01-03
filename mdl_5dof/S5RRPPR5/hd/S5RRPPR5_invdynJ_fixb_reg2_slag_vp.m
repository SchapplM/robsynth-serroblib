% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPPR5
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
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPPR5_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR5_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR5_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:30:03
% EndTime: 2019-12-31 19:30:07
% DurationCPUTime: 2.75s
% Computational Cost: add. (3372->372), mult. (7812->454), div. (0->0), fcn. (5527->10), ass. (0->201)
t162 = sin(qJ(5));
t165 = cos(qJ(5));
t159 = sin(pkin(8));
t160 = cos(pkin(8));
t166 = cos(qJ(2));
t234 = t160 * t166;
t213 = qJD(1) * t234;
t163 = sin(qJ(2));
t227 = qJD(1) * t163;
t103 = t159 * t227 - t213;
t116 = t159 * t166 + t160 * t163;
t277 = t116 * qJD(1);
t56 = t162 * t103 + t165 * t277;
t105 = t116 * qJD(2);
t219 = t166 * qJDD(1);
t220 = t163 * qJDD(1);
t196 = t159 * t220 - t160 * t219;
t67 = qJD(1) * t105 + t196;
t221 = qJD(1) * qJD(2);
t210 = t163 * t221;
t180 = qJDD(1) * t116 - t159 * t210;
t209 = t166 * t221;
t68 = t160 * t209 + t180;
t14 = qJD(5) * t56 + t162 * t68 - t165 * t67;
t152 = qJD(2) - qJD(5);
t290 = t56 * t152;
t292 = t14 + t290;
t252 = qJ(3) + pkin(6);
t211 = t252 * t163;
t120 = qJD(1) * t211;
t124 = t252 * t166;
t121 = qJD(1) * t124;
t236 = t159 * t121;
t74 = -t160 * t120 - t236;
t231 = qJD(4) - t74;
t257 = t56 ^ 2;
t194 = -t165 * t103 + t162 * t277;
t258 = t194 ^ 2;
t291 = -t258 + t257;
t256 = t56 * t194;
t244 = t68 * qJ(4);
t259 = t166 * pkin(2);
t145 = pkin(1) + t259;
t94 = pkin(2) * t210 - t145 * qJDD(1) + qJDD(3);
t289 = t244 - t94;
t261 = t277 * pkin(7);
t288 = -t261 + t231;
t123 = -t145 * qJD(1) + qJD(3);
t287 = -t277 * qJ(4) + t123;
t224 = qJD(5) * t165;
t225 = qJD(5) * t162;
t187 = t103 * t224 + t162 * t67 + t165 * t68 - t225 * t277;
t281 = t194 * t152;
t286 = t187 - t281;
t153 = qJ(2) + pkin(8);
t146 = sin(t153);
t147 = cos(t153);
t193 = t146 * t162 + t147 * t165;
t270 = pkin(3) + pkin(4);
t29 = -t270 * t103 - t287;
t164 = sin(qJ(1));
t97 = t146 * t165 - t147 * t162;
t79 = t97 * t164;
t167 = cos(qJ(1));
t237 = t147 * t167;
t239 = t146 * t167;
t81 = t162 * t237 - t165 * t239;
t285 = g(1) * t81 - g(2) * t79 + g(3) * t193 - t29 * t56;
t206 = qJD(2) * t252;
t183 = -t163 * qJD(3) - t166 * t206;
t63 = qJDD(2) * pkin(2) + qJD(1) * t183 - qJDD(1) * t211;
t98 = t166 * qJD(3) - t163 * t206;
t72 = qJD(1) * t98 + qJDD(1) * t124;
t251 = t159 * t72 - t160 * t63;
t214 = -qJDD(4) - t251;
t11 = -t68 * pkin(7) - t270 * qJDD(2) - t214;
t154 = qJDD(2) * qJ(4);
t155 = qJD(2) * qJD(4);
t24 = t159 * t63 + t160 * t72;
t20 = t154 + t155 + t24;
t12 = t67 * pkin(7) + t20;
t114 = qJD(2) * pkin(2) - t120;
t65 = t160 * t114 - t236;
t197 = qJD(4) - t65;
t31 = -t270 * qJD(2) + t197 - t261;
t262 = t103 * pkin(7);
t235 = t160 * t121;
t66 = t159 * t114 + t235;
t57 = qJD(2) * qJ(4) + t66;
t37 = t57 + t262;
t9 = t162 * t31 + t165 * t37;
t2 = -t9 * qJD(5) + t165 * t11 - t162 * t12;
t284 = -t9 * t152 + t2;
t271 = t277 ^ 2;
t99 = t103 ^ 2;
t283 = -t99 - t271;
t282 = -t99 + t271;
t280 = -t147 * pkin(3) - t146 * qJ(4);
t279 = g(1) * t164 - g(2) * t167;
t148 = g(2) * t164;
t278 = g(1) * t167 + t148;
t1 = t162 * t11 + t165 * t12 + t31 * t224 - t225 * t37;
t80 = t193 * t164;
t82 = t193 * t167;
t275 = g(1) * t82 + g(2) * t80 + g(3) * t97 + t194 * t29 - t1;
t223 = t277 * qJD(4);
t274 = -t270 * t67 + t223;
t78 = t160 * t124 - t159 * t211;
t273 = -t78 * qJDD(2) - t146 * t279;
t272 = g(3) * t146 + t74 * qJD(2) + t147 * t278 - t24;
t266 = t67 * pkin(3);
t265 = pkin(2) * t163;
t263 = g(3) * t166;
t260 = t147 * pkin(4);
t8 = -t162 * t37 + t165 * t31;
t255 = t8 * t152;
t253 = pkin(7) - t252;
t73 = -t159 * t120 + t235;
t41 = t73 + t262;
t144 = -t160 * pkin(2) - pkin(3);
t137 = -pkin(4) + t144;
t140 = t159 * pkin(2) + qJ(4);
t83 = t165 * t137 - t162 * t140;
t250 = qJD(5) * t83 - t162 * t41 + t288 * t165;
t84 = t162 * t137 + t165 * t140;
t249 = -qJD(5) * t84 - t288 * t162 - t165 * t41;
t46 = t159 * t183 + t160 * t98;
t243 = pkin(6) * qJDD(1);
t242 = qJDD(2) * pkin(3);
t240 = t277 * t103;
t238 = t147 * t164;
t233 = t167 * t252;
t232 = t73 * qJD(2);
t157 = t163 ^ 2;
t158 = t166 ^ 2;
t229 = t157 - t158;
t228 = t157 + t158;
t226 = qJD(2) * t163;
t217 = pkin(2) * t226;
t169 = qJD(1) ^ 2;
t216 = t163 * t169 * t166;
t215 = -g(1) * t239 + g(3) * t147 - t146 * t148;
t45 = t159 * t98 - t160 * t183;
t77 = t159 * t124 + t160 * t211;
t205 = t152 ^ 2;
t133 = t167 * t145;
t204 = g(2) * (pkin(3) * t237 + qJ(4) * t239 + t133);
t203 = t163 * t209;
t202 = t259 - t280;
t201 = -pkin(3) * t146 - t265;
t198 = -t215 - t251;
t47 = -t116 * pkin(7) + t77;
t115 = t159 * t163 - t234;
t48 = t115 * pkin(7) + t78;
t18 = -t162 * t48 + t165 * t47;
t19 = t162 * t47 + t165 * t48;
t195 = t103 * t105 + t67 * t115;
t71 = t162 * t115 + t165 * t116;
t192 = qJD(2) * t105 + qJDD(2) * t115;
t191 = t116 * qJ(4) + t145;
t190 = -t145 + t280;
t189 = -pkin(2) * t227 - t103 * qJ(4);
t188 = -0.2e1 * pkin(1) * t221 - pkin(6) * qJDD(2);
t184 = g(1) * t238 - g(2) * t237 - t77 * qJDD(2);
t182 = -t46 * t103 + t277 * t45 - t78 * t67 + t77 * t68 - t278;
t43 = t103 * pkin(3) + t287;
t181 = -t277 * t43 - qJDD(4) + t198;
t108 = qJD(2) * t234 - t159 * t226;
t179 = t108 * qJ(4) + t116 * qJD(4) - t217;
t168 = qJD(2) ^ 2;
t178 = 0.2e1 * qJDD(1) * pkin(1) - pkin(6) * t168 + t279;
t177 = pkin(1) * t169 - t243 + t278;
t176 = t108 * t103 + t105 * t277 + t68 * t115 + t116 * t67;
t174 = -t279 + t94;
t172 = t174 - t244;
t171 = 0.2e1 * t277 * qJD(2) + t196;
t151 = qJDD(2) - qJDD(5);
t127 = qJ(4) * t237;
t125 = qJ(4) * t238;
t70 = -t165 * t115 + t162 * t116;
t69 = t108 * qJD(2) + t116 * qJDD(2);
t64 = t115 * pkin(3) - t191;
t49 = -qJD(2) * pkin(3) + t197;
t44 = pkin(3) * t277 - t189;
t40 = -t270 * t115 + t191;
t39 = (t103 + t213) * qJD(2) + t180;
t38 = (-t103 + t213) * qJD(2) + t180;
t36 = t105 * pkin(3) - t179;
t34 = t105 * pkin(7) + t46;
t33 = -t108 * pkin(7) + t45;
t32 = -t270 * t277 + t189;
t27 = -t270 * t105 + t179;
t26 = qJD(5) * t71 - t165 * t105 + t162 * t108;
t25 = -t162 * t105 - t165 * t108 - t115 * t224 + t116 * t225;
t22 = -t214 - t242;
t21 = t108 * t277 + t68 * t116;
t15 = -t223 + t266 - t289;
t5 = t274 + t289;
t4 = -qJD(5) * t19 - t162 * t34 + t165 * t33;
t3 = qJD(5) * t18 + t162 * t33 + t165 * t34;
t6 = [0, 0, 0, 0, 0, qJDD(1), t279, t278, 0, 0, t157 * qJDD(1) + 0.2e1 * t203, 0.2e1 * t163 * t219 - 0.2e1 * t221 * t229, qJDD(2) * t163 + t168 * t166, t158 * qJDD(1) - 0.2e1 * t203, qJDD(2) * t166 - t168 * t163, 0, t163 * t188 + t166 * t178, -t163 * t178 + t166 * t188, 0.2e1 * t228 * t243 - t278, -g(1) * (-t164 * pkin(1) + t167 * pkin(6)) - g(2) * (t167 * pkin(1) + t164 * pkin(6)) + (pkin(6) ^ 2 * t228 + pkin(1) ^ 2) * qJDD(1), t21, -t176, t69, t195, -t192, 0, t123 * t105 + t94 * t115 - t145 * t67 + (t103 * t265 - t45) * qJD(2) + t184, t123 * t108 + t94 * t116 - t145 * t68 + (t265 * t277 - t46) * qJD(2) + t273, -t66 * t105 - t65 * t108 - t24 * t115 + t116 * t251 + t182, t24 * t78 + t66 * t46 + t251 * t77 - t65 * t45 - t94 * t145 + t123 * t217 - g(1) * (-t164 * t145 + t233) - g(2) * (t164 * t252 + t133), t21, t69, t176, 0, t192, t195, -t45 * qJD(2) + t36 * t103 + t43 * t105 + t15 * t115 + t64 * t67 + t184, -t57 * t105 + t49 * t108 - t20 * t115 + t22 * t116 + t182, t46 * qJD(2) - t43 * t108 - t15 * t116 - t277 * t36 - t64 * t68 - t273, t20 * t78 + t57 * t46 + t15 * t64 + t43 * t36 + t22 * t77 + t49 * t45 - g(1) * t233 - t204 + (-g(1) * t190 - g(2) * t252) * t164, t187 * t71 - t56 * t25, -t71 * t14 - t187 * t70 + t194 * t25 - t56 * t26, -t71 * t151 + t25 * t152, t14 * t70 + t194 * t26, t70 * t151 + t26 * t152, 0, g(1) * t80 - g(2) * t82 + t40 * t14 - t18 * t151 - t4 * t152 + t194 * t27 + t29 * t26 + t5 * t70, g(1) * t79 + g(2) * t81 + t19 * t151 + t3 * t152 + t187 * t40 - t29 * t25 + t27 * t56 + t5 * t71, -t1 * t70 - t19 * t14 - t18 * t187 - t194 * t3 - t2 * t71 + t8 * t25 - t9 * t26 - t4 * t56 + t278, t1 * t19 + t9 * t3 + t2 * t18 + t8 * t4 + t5 * t40 + t29 * t27 - t204 + (g(1) * t253 - g(2) * t260) * t167 + (-g(1) * (t190 - t260) + g(2) * t253) * t164; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t216, t229 * t169, t220, t216, t219, qJDD(2), t163 * t177 - t263, g(3) * t163 + t166 * t177, 0, 0, t240, t282, t39, -t240, -t196, qJDD(2), t232 - t123 * t277 + (qJDD(2) * t160 - t103 * t227) * pkin(2) + t198, t123 * t103 + (-qJDD(2) * t159 - t227 * t277) * pkin(2) + t272, (t66 - t73) * t277 + (-t65 + t74) * t103 + (-t159 * t67 - t160 * t68) * pkin(2), t65 * t73 - t66 * t74 + (-t263 + t159 * t24 - t160 * t251 + (-qJD(1) * t123 + t278) * t163) * pkin(2), t240, t39, -t282, qJDD(2), t196, -t240, t232 - t44 * t103 + (pkin(3) - t144) * qJDD(2) + t181, -t140 * t67 + t144 * t68 + (t57 - t73) * t277 + (t49 - t231) * t103, t140 * qJDD(2) - t43 * t103 + t277 * t44 + t154 + 0.2e1 * t155 - t272, t20 * t140 + t22 * t144 - t43 * t44 - t49 * t73 - g(1) * (t167 * t201 + t127) - g(2) * (t164 * t201 + t125) - g(3) * t202 + t231 * t57, -t256, -t291, -t286, t256, t292, t151, -t83 * t151 - t152 * t249 - t194 * t32 - t2 - t285, t84 * t151 + t152 * t250 - t32 * t56 - t275, -t83 * t187 - t84 * t14 + (-t250 + t8) * t194 + (-t249 - t9) * t56, t1 * t84 + t2 * t83 - t29 * t32 - g(1) * (-t167 * t265 + t127) - g(2) * (-t164 * t265 + t125) - g(3) * (t202 + t260) + t250 * t9 + t249 * t8 + t278 * t146 * t270; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t171, t38, t283, t66 * t103 + t277 * t65 + t174, 0, 0, 0, 0, 0, 0, t171, t283, -t38, t266 + t57 * t103 + (-qJD(4) - t49) * t277 + t172, 0, 0, 0, 0, 0, 0, -t14 + t290, -t187 - t281, t257 + t258, -t194 * t9 - t56 * t8 + t172 - t274; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(2) + t240, t39, -t271 - t168, -t57 * qJD(2) - t181 - t242, 0, 0, 0, 0, 0, 0, -t165 * t151 - t162 * t205 - t194 * t277, t162 * t151 - t165 * t205 - t277 * t56, -t292 * t162 - t286 * t165, -t29 * t277 + t284 * t165 + (t1 + t255) * t162 + t215; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t256, t291, t286, -t256, -t292, -t151, t284 + t285, -t255 + t275, 0, 0;];
tau_reg = t6;
