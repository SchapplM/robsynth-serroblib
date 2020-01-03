% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5RRPRP6
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRP6_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP6_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP6_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP6_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:58:33
% EndTime: 2019-12-31 19:58:39
% DurationCPUTime: 3.35s
% Computational Cost: add. (5032->452), mult. (11864->556), div. (0->0), fcn. (8455->10), ass. (0->227)
t162 = sin(qJ(1));
t165 = cos(qJ(1));
t291 = g(1) * t162 - g(2) * t165;
t157 = sin(pkin(8));
t164 = cos(qJ(2));
t161 = sin(qJ(2));
t249 = cos(pkin(8));
t208 = t249 * t161;
t124 = t157 * t164 + t208;
t114 = t124 * qJD(1);
t160 = sin(qJ(4));
t163 = cos(qJ(4));
t207 = t249 * t164;
t140 = qJD(1) * t207;
t227 = qJD(1) * qJD(2);
t215 = t161 * t227;
t174 = qJDD(1) * t124 - t157 * t215;
t171 = qJD(2) * t140 + t174;
t226 = qJD(2) * qJD(4);
t229 = qJD(4) * t160;
t182 = t160 * qJDD(2) - t114 * t229 + (t171 + t226) * t163;
t266 = qJ(5) * t182;
t113 = t124 * qJD(2);
t225 = t161 * qJDD(1);
t195 = -qJDD(1) * t207 + t157 * t225;
t86 = qJD(1) * t113 + t195;
t81 = qJDD(4) + t86;
t287 = pkin(4) * t81;
t230 = qJD(1) * t161;
t111 = t157 * t230 - t140;
t272 = t164 * pkin(2);
t147 = pkin(1) + t272;
t130 = -t147 * qJD(1) + qJD(3);
t56 = pkin(3) * t111 - pkin(7) * t114 + t130;
t270 = qJ(3) + pkin(6);
t216 = t270 * t161;
t126 = qJD(1) * t216;
t264 = qJD(2) * pkin(2);
t120 = -t126 + t264;
t132 = t270 * t164;
t127 = qJD(1) * t132;
t209 = t249 * t127;
t85 = t157 * t120 + t209;
t75 = qJD(2) * pkin(7) + t85;
t34 = t160 * t56 + t163 * t75;
t223 = pkin(2) * t215 + qJDD(3);
t224 = t164 * qJDD(1);
t246 = qJDD(1) * pkin(1);
t35 = -pkin(2) * t224 + t86 * pkin(3) - pkin(7) * t171 + t223 - t246;
t210 = qJD(2) * t270;
t179 = -qJD(3) * t161 - t164 * t210;
t80 = qJDD(2) * pkin(2) + t179 * qJD(1) - qJDD(1) * t216;
t109 = qJD(3) * t164 - t161 * t210;
t87 = t109 * qJD(1) + qJDD(1) * t132;
t45 = t157 * t80 + t249 * t87;
t43 = qJDD(2) * pkin(7) + t45;
t7 = -qJD(4) * t34 - t160 * t43 + t163 * t35;
t97 = qJD(2) * t160 + t114 * t163;
t1 = -qJD(5) * t97 - t266 + t287 + t7;
t104 = qJD(4) + t111;
t95 = -t163 * qJD(2) + t114 * t160;
t26 = -qJ(5) * t95 + t34;
t263 = t104 * t26;
t292 = t1 + t263;
t154 = qJ(2) + pkin(8);
t148 = sin(t154);
t199 = g(1) * t165 + g(2) * t162;
t187 = t199 * t148;
t149 = cos(t154);
t235 = t163 * t165;
t238 = t160 * t162;
t105 = t149 * t238 + t235;
t236 = t162 * t163;
t237 = t160 * t165;
t107 = -t149 * t237 + t236;
t275 = g(3) * t160;
t290 = -g(1) * t107 + g(2) * t105 + t148 * t275;
t277 = g(3) * t148;
t178 = -t199 * t149 - t277;
t289 = t97 ^ 2;
t288 = t114 ^ 2;
t286 = t95 * pkin(4);
t285 = pkin(2) * t157;
t284 = pkin(2) * t161;
t138 = t165 * t147;
t279 = g(2) * t138;
t276 = g(3) * t149;
t274 = g(3) * t164;
t273 = t163 * pkin(4);
t271 = t97 * t95;
t33 = -t160 * t75 + t163 * t56;
t25 = -qJ(5) * t97 + t33;
t23 = pkin(4) * t104 + t25;
t269 = -t25 + t23;
t228 = qJD(4) * t163;
t170 = -t163 * qJDD(2) + t160 * t171;
t247 = qJD(4) * t97;
t49 = t170 + t247;
t268 = -t160 * t49 - t95 * t228;
t44 = -t157 * t87 + t249 * t80;
t66 = pkin(2) * t230 + pkin(3) * t114 + pkin(7) * t111;
t117 = t157 * t127;
t89 = -t249 * t126 - t117;
t40 = t160 * t66 + t163 * t89;
t183 = -t157 * t161 + t207;
t83 = -pkin(3) * t183 - pkin(7) * t124 - t147;
t93 = t249 * t132 - t157 * t216;
t90 = t163 * t93;
t51 = t160 * t83 + t90;
t143 = pkin(7) + t285;
t234 = qJ(5) + t143;
t206 = qJD(4) * t234;
t244 = t111 * t160;
t267 = -qJ(5) * t244 + qJD(5) * t163 - t160 * t206 - t40;
t265 = qJ(5) * t49;
t262 = t104 * t33;
t261 = t104 * t34;
t260 = t104 * t95;
t259 = t114 * t95;
t258 = t114 * t97;
t257 = t160 * t182;
t256 = t160 * t81;
t255 = t160 * t95;
t254 = t160 * t97;
t253 = t163 * t49;
t252 = t163 * t95;
t251 = t97 * t104;
t39 = -t160 * t89 + t163 * t66;
t250 = -pkin(4) * t114 - qJD(5) * t160 - t39 + (-qJ(5) * t111 - t206) * t163;
t248 = pkin(6) * qJDD(1);
t245 = t104 * t114;
t243 = t114 * t111;
t116 = t183 * qJD(2);
t242 = t116 * t160;
t241 = t116 * t163;
t240 = t124 * t160;
t239 = t124 * t163;
t233 = (g(1) * t235 + g(2) * t236) * t148;
t155 = t161 ^ 2;
t156 = t164 ^ 2;
t232 = t155 - t156;
t231 = t155 + t156;
t65 = t249 * t109 + t157 * t179;
t221 = t161 * t264;
t67 = pkin(3) * t113 - pkin(7) * t116 + t221;
t222 = t160 * t67 + t163 * t65 + t83 * t228;
t167 = qJD(1) ^ 2;
t220 = t161 * t167 * t164;
t219 = t249 * pkin(2);
t218 = t124 * t228;
t42 = -qJDD(2) * pkin(3) - t44;
t16 = pkin(4) * t49 + qJDD(5) + t42;
t217 = -t16 - t276;
t214 = qJD(5) + t286;
t213 = pkin(4) * t160 + t270;
t211 = -t160 * t65 + t163 * t67;
t50 = -t160 * t93 + t163 * t83;
t6 = t160 * t35 + t163 * t43 + t56 * t228 - t75 * t229;
t64 = t109 * t157 - t249 * t179;
t88 = -t126 * t157 + t209;
t92 = t132 * t157 + t270 * t208;
t205 = t104 * t163;
t204 = t164 * t215;
t203 = t291 * t148;
t144 = -t219 - pkin(3);
t202 = pkin(3) * t149 + pkin(7) * t148;
t201 = -g(1) * t105 - g(2) * t107;
t106 = -t149 * t236 + t237;
t108 = t149 * t235 + t238;
t200 = -g(1) * t106 - g(2) * t108;
t197 = qJD(4) * t104 * t143 + t42;
t2 = -qJD(5) * t95 - t265 + t6;
t196 = -t104 * t23 + t2;
t194 = -t160 * t34 - t163 * t33;
t193 = -t252 - t254;
t84 = t249 * t120 - t117;
t146 = pkin(3) + t273;
t158 = -qJ(5) - pkin(7);
t192 = t146 * t149 - t148 * t158;
t191 = -qJ(5) * t116 - qJD(5) * t124;
t190 = t163 * t81 + (-t229 - t244) * t104;
t188 = t163 * t182 - t97 * t229;
t186 = -0.2e1 * pkin(1) * t227 - pkin(6) * qJDD(2);
t185 = t218 + t242;
t74 = -qJD(2) * pkin(3) - t84;
t184 = t104 * t74 - t143 * t81;
t103 = -t147 * qJDD(1) + t223;
t177 = g(1) * t108 - g(2) * t106 + t163 * t277 - t6;
t166 = qJD(2) ^ 2;
t176 = -pkin(6) * t166 + 0.2e1 * t246 + t291;
t175 = pkin(1) * t167 + t199 - t248;
t173 = t7 + t290;
t169 = t114 * t228 + t160 * t226 + t170;
t136 = t149 * t275;
t131 = t144 - t273;
t122 = t234 * t163;
t121 = t234 * t160;
t110 = t111 ^ 2;
t94 = t95 ^ 2;
t62 = pkin(4) * t240 + t92;
t53 = -pkin(4) * t244 + t88;
t52 = t214 + t74;
t46 = -t94 + t289;
t41 = pkin(4) * t185 + t64;
t37 = t104 * t113 - t183 * t81;
t36 = -qJ(5) * t240 + t51;
t29 = -pkin(4) * t183 - qJ(5) * t239 + t50;
t28 = -t169 + t251;
t27 = t182 + t260;
t22 = -t104 ^ 2 * t163 - t256 - t258;
t21 = t104 * t205 + t256 - t258;
t20 = t190 + t259;
t19 = t190 - t259;
t18 = t104 * t255 - t253;
t17 = t97 * t205 + t257;
t15 = -t51 * qJD(4) + t211;
t14 = -t93 * t229 + t222;
t13 = -t268 * t124 + t95 * t242;
t12 = t124 * t188 + t97 * t241;
t11 = -qJ(5) * t218 + (-qJD(4) * t93 + t191) * t160 + t222;
t10 = pkin(4) * t113 + t191 * t163 + (-t90 + (qJ(5) * t124 - t83) * t160) * qJD(4) + t211;
t9 = -t104 * t185 - t113 * t95 + t183 * t49 - t81 * t240;
t8 = t81 * t239 + t113 * t97 - t183 * t182 + (-t124 * t229 + t241) * t104;
t5 = t193 * t111 + t188 + t268;
t4 = (-t252 + t254) * t111 - t188 + t268;
t3 = t193 * t116 + (-t257 - t253 + (-t163 * t97 + t255) * qJD(4)) * t124;
t24 = [0, 0, 0, 0, 0, qJDD(1), t291, t199, 0, 0, qJDD(1) * t155 + 0.2e1 * t204, 0.2e1 * t161 * t224 - 0.2e1 * t232 * t227, qJDD(2) * t161 + t164 * t166, qJDD(1) * t156 - 0.2e1 * t204, qJDD(2) * t164 - t161 * t166, 0, t161 * t186 + t164 * t176, -t161 * t176 + t164 * t186, 0.2e1 * t231 * t248 - t199, -g(1) * (-pkin(1) * t162 + pkin(6) * t165) - g(2) * (pkin(1) * t165 + pkin(6) * t162) + (t231 * pkin(6) ^ 2 + pkin(1) ^ 2) * qJDD(1), t114 * t116 + t124 * t171, -t116 * t111 - t114 * t113 - t124 * t86 + t171 * t183, qJD(2) * t116 + qJDD(2) * t124, t111 * t113 - t183 * t86, -qJD(2) * t113 + qJDD(2) * t183, 0, -qJDD(2) * t92 - t103 * t183 + t113 * t130 - t147 * t86 + t291 * t149 + (t111 * t284 - t64) * qJD(2), -t65 * qJD(2) - t93 * qJDD(2) + t103 * t124 + t114 * t221 + t130 * t116 - t147 * t171 - t203, -t65 * t111 - t85 * t113 + t64 * t114 - t84 * t116 - t44 * t124 + t171 * t92 + t183 * t45 - t93 * t86 - t199, t45 * t93 + t85 * t65 - t44 * t92 - t84 * t64 - t103 * t147 + t130 * t221 - g(1) * (-t147 * t162 + t165 * t270) - g(2) * (t162 * t270 + t138), t12, t3, t8, t13, t9, t37, t74 * t242 + t104 * t15 + t113 * t33 - t183 * t7 + t49 * t92 + t50 * t81 + t64 * t95 + (t160 * t42 + t228 * t74) * t124 + t200, t74 * t241 - t104 * t14 - t113 * t34 + t183 * t6 + t182 * t92 - t51 * t81 + t64 * t97 + (t163 * t42 - t229 * t74) * t124 + t201, -t14 * t95 - t15 * t97 - t182 * t50 - t49 * t51 + t194 * t116 + (-t160 * t6 - t163 * t7 + (t160 * t33 - t163 * t34) * qJD(4)) * t124 + t203, -t279 + t34 * t14 + t33 * t15 + t42 * t92 + t7 * t50 + t6 * t51 + t74 * t64 + (-g(1) * t270 - g(2) * t202) * t165 + (-g(1) * (-t147 - t202) - g(2) * t270) * t162, t12, t3, t8, t13, t9, t37, t52 * t242 - t1 * t183 + t10 * t104 + t113 * t23 + t29 * t81 + t41 * t95 + t49 * t62 + (t16 * t160 + t228 * t52) * t124 + t200, t52 * t241 - t104 * t11 - t113 * t26 + t183 * t2 - t36 * t81 + t41 * t97 + t182 * t62 + (t16 * t163 - t229 * t52) * t124 + t201, -t10 * t97 - t11 * t95 - t29 * t182 - t36 * t49 + (-t160 * t26 - t163 * t23) * t116 + (-t1 * t163 - t160 * t2 + (t160 * t23 - t163 * t26) * qJD(4)) * t124 + t203, -t279 + t1 * t29 + t23 * t10 + t26 * t11 + t16 * t62 + t2 * t36 + t52 * t41 + (-g(1) * t213 - g(2) * t192) * t165 + (-g(1) * (-t147 - t192) - g(2) * t213) * t162; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t220, t232 * t167, t225, t220, t224, qJDD(2), t161 * t175 - t274, g(3) * t161 + t164 * t175, 0, 0, t243, -t110 + t288, (t140 + t111) * qJD(2) + t174, -t243, -t195, qJDD(2), -t276 + t88 * qJD(2) - t130 * t114 + t187 + (t249 * qJDD(2) - t111 * t230) * pkin(2) + t44, qJD(2) * t89 + t111 * t130 + (-qJDD(2) * t157 - t114 * t230) * pkin(2) - t45 - t178, -t171 * t219 - t86 * t285 - (-t85 + t88) * t114 + (t89 - t84) * t111, t84 * t88 - t85 * t89 + (t249 * t44 - t274 + t157 * t45 + (-qJD(1) * t130 + t199) * t161) * pkin(2), t17, t5, t21, t18, t20, -t245, -t104 * t39 - t114 * t33 + t144 * t49 - t88 * t95 + (-t197 - t276) * t163 + t184 * t160 + t233, t104 * t40 + t114 * t34 + t144 * t182 - t88 * t97 + t136 + t184 * t163 + (-t187 + t197) * t160, t39 * t97 + t40 * t95 + (-t111 * t33 - t143 * t49 + t6 + (t143 * t97 - t33) * qJD(4)) * t163 + (-t111 * t34 + t143 * t182 - t7 + (t143 * t95 - t34) * qJD(4)) * t160 + t178, t42 * t144 - t34 * t40 - t33 * t39 - t74 * t88 - g(3) * (t202 + t272) + (qJD(4) * t194 - t7 * t160 + t6 * t163) * t143 + t199 * (pkin(3) * t148 - pkin(7) * t149 + t284), t17, t5, t21, t18, t20, -t245, -t114 * t23 - t121 * t81 + t131 * t49 - t53 * t95 + t217 * t163 + t250 * t104 + (t111 * t52 + (t52 + t286) * qJD(4)) * t160 + t233, t114 * t26 - t122 * t81 + t131 * t182 - t53 * t97 + t136 + t52 * t205 - t267 * t104 + (pkin(4) * t247 + t16 - t187) * t160, t121 * t182 - t122 * t49 - t292 * t160 + t196 * t163 - t250 * t97 - t267 * t95 + t178, t2 * t122 - t1 * t121 + t16 * t131 - g(3) * (t192 + t272) + (pkin(4) * t229 - t53) * t52 + t267 * t26 + t250 * t23 + t199 * (t146 * t148 + t149 * t158 + t284); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t114 * qJD(2) + t195, (t140 - t111) * qJD(2) + t174, -t110 - t288, t111 * t85 + t114 * t84 + t103 - t291, 0, 0, 0, 0, 0, 0, t19, t22, t4, -t114 * t74 + (t7 + t261) * t163 + (t6 - t262) * t160 - t291, 0, 0, 0, 0, 0, 0, t19, t22, t4, -t114 * t52 + t196 * t160 + t292 * t163 - t291; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t271, t46, t27, -t271, t28, t81, -t74 * t97 + t173 + t261, t74 * t95 + t177 + t262, 0, 0, t271, t46, t27, -t271, t28, t81, 0.2e1 * t287 - t266 + t263 + (-t214 - t52) * t97 + t173, -pkin(4) * t289 + t265 + t104 * t25 + (qJD(5) + t52) * t95 + t177, -pkin(4) * t182 - t269 * t95, t269 * t26 + (-t52 * t97 + t1 + t290) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t169 + t251, t182 - t260, -t94 - t289, t23 * t97 + t26 * t95 - t187 - t217;];
tau_reg = t24;
