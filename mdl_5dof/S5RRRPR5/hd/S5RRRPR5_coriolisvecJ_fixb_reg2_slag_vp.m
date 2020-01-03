% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:14:17
% EndTime: 2019-12-31 21:14:26
% DurationCPUTime: 3.33s
% Computational Cost: add. (7895->373), mult. (20651->505), div. (0->0), fcn. (15067->8), ass. (0->212)
t167 = sin(qJ(5));
t170 = cos(qJ(5));
t163 = qJD(2) + qJD(3);
t166 = sin(pkin(9));
t168 = sin(qJ(3));
t169 = sin(qJ(2));
t171 = cos(qJ(3));
t172 = cos(qJ(2));
t135 = t168 * t172 + t171 * t169;
t125 = qJD(1) * t135;
t120 = t125 * qJ(4);
t265 = -pkin(7) - pkin(6);
t146 = t265 * t172;
t141 = qJD(1) * t146;
t126 = t168 * t141;
t145 = t265 * t169;
t139 = qJD(1) * t145;
t247 = qJD(2) * pkin(2);
t130 = t139 + t247;
t93 = t171 * t130 + t126;
t74 = -t120 + t93;
t67 = t163 * pkin(3) + t74;
t232 = cos(pkin(9));
t225 = t171 * t172;
t205 = qJD(1) * t225;
t218 = qJD(1) * t169;
t206 = t168 * t218;
t123 = -t205 + t206;
t231 = t123 * qJ(4);
t226 = t171 * t141;
t94 = t168 * t130 - t226;
t75 = t94 - t231;
t70 = t232 * t75;
t42 = t166 * t67 + t70;
t40 = t163 * pkin(8) + t42;
t160 = -t172 * pkin(2) - pkin(1);
t144 = qJD(1) * t160;
t102 = t123 * pkin(3) + qJD(4) + t144;
t182 = -t166 * t123 + t232 * t125;
t87 = t232 * t123 + t166 * t125;
t46 = t87 * pkin(4) - pkin(8) * t182 + t102;
t17 = t167 * t46 + t170 * t40;
t271 = qJD(5) + t87;
t207 = qJD(2) * t265;
t195 = qJD(1) * t207;
t131 = t169 * t195;
t132 = t172 * t195;
t198 = -t168 * t131 + t171 * t132;
t61 = -qJD(3) * t94 + t198;
t227 = t168 * t169;
t190 = t163 * t227;
t213 = qJD(1) * qJD(2);
t203 = t172 * t213;
t220 = -qJD(3) * t205 - t171 * t203;
t91 = qJD(1) * t190 + t220;
t177 = t91 * qJ(4) - t125 * qJD(4) + t61;
t216 = qJD(3) * t171;
t217 = qJD(3) * t168;
t197 = -t130 * t216 - t171 * t131 - t168 * t132 - t141 * t217;
t101 = t163 * t135;
t92 = t101 * qJD(1);
t34 = -t92 * qJ(4) - t123 * qJD(4) - t197;
t15 = t166 * t177 + t232 * t34;
t58 = -t166 * t91 + t232 * t92;
t59 = -t166 * t92 - t232 * t91;
t161 = pkin(2) * t218;
t81 = t92 * pkin(3) + qJD(2) * t161;
t25 = t58 * pkin(4) - t59 * pkin(8) + t81;
t3 = -qJD(5) * t17 - t167 * t15 + t170 * t25;
t274 = t17 * t271 + t3;
t186 = t167 * t40 - t170 * t46;
t2 = -t186 * qJD(5) + t170 * t15 + t167 * t25;
t273 = t186 * t271 + t2;
t272 = t87 * t182;
t200 = t170 * t271;
t245 = t167 * t58;
t270 = -t200 * t271 - t245;
t162 = t169 * t247;
t269 = 0.2e1 * t162;
t268 = -0.2e1 * t213;
t215 = qJD(5) * t167;
t266 = -t170 * t58 + t215 * t271;
t104 = t168 * t145 - t171 * t146;
t264 = t125 * pkin(3);
t14 = t166 * t34 - t232 * t177;
t103 = t171 * t145 + t168 * t146;
t183 = -t135 * qJ(4) + t103;
t134 = -t225 + t227;
t84 = -t134 * qJ(4) + t104;
t53 = t166 * t84 - t232 * t183;
t263 = t14 * t53;
t97 = -t166 * t134 + t232 * t135;
t262 = t14 * t97;
t1 = t2 * t170;
t261 = t3 * t167;
t96 = t232 * t134 + t166 * t135;
t260 = t58 * t96;
t76 = -t170 * t163 + t167 * t182;
t259 = t76 * t87;
t78 = t167 * t163 + t170 * t182;
t258 = t78 * t76;
t257 = t78 * t182;
t256 = t271 * t182;
t255 = t182 * t76;
t253 = t87 ^ 2;
t252 = t182 ^ 2;
t250 = t97 * t58;
t214 = qJD(5) * t170;
t244 = t167 * t59;
t32 = qJD(5) * t78 + t244;
t249 = -t167 * t32 - t76 * t214;
t248 = pkin(2) * qJD(3);
t246 = t166 * t75;
t243 = t167 * t76;
t242 = t167 * t78;
t241 = t167 * t87;
t240 = t170 * t78;
t239 = t170 * t87;
t31 = -t163 * t214 - t170 * t59 + t182 * t215;
t238 = t31 * t167;
t237 = t32 * t170;
t236 = t182 * t163;
t235 = t87 * t163;
t98 = -t168 * t139 + t226;
t179 = t98 + t231;
t202 = t232 * t168;
t99 = t171 * t139 + t126;
t79 = -t120 + t99;
t234 = -t166 * t79 + t232 * t179 + (t166 * t171 + t202) * t248;
t228 = t166 * t168;
t117 = (t232 * t171 - t228) * t248;
t49 = t166 * t179 + t232 * t79;
t233 = t117 - t49;
t230 = t125 * t123;
t229 = t144 * t125;
t174 = qJD(1) ^ 2;
t224 = t172 * t174;
t173 = qJD(2) ^ 2;
t223 = t173 * t169;
t222 = t173 * t172;
t159 = t171 * pkin(2) + pkin(3);
t119 = pkin(2) * t202 + t166 * t159;
t219 = t169 ^ 2 - t172 ^ 2;
t212 = -t17 * t241 + t186 * t239 + t1;
t210 = t97 * t215;
t209 = t97 * t214;
t208 = t169 * t224;
t41 = t232 * t67 - t246;
t39 = -t163 * pkin(4) - t41;
t35 = t39 * t215;
t36 = t39 * t214;
t204 = t182 * t186 + t35;
t90 = t101 * pkin(3) + t162;
t201 = t167 * t271;
t199 = pkin(1) * t268;
t196 = t14 * t167 + t17 * t182 + t36;
t194 = t169 * t203;
t100 = -qJD(2) * t225 - t172 * t216 + t190;
t63 = -t232 * t100 - t166 * t101;
t193 = t39 * t63 + t262;
t192 = t182 * t42 - t41 * t87;
t191 = t271 * t63 + t250;
t188 = -t167 * t17 + t170 * t186;
t187 = -t167 * t186 - t17 * t170;
t106 = t134 * pkin(3) + t160;
t52 = t96 * pkin(4) - t97 * pkin(8) + t106;
t54 = t166 * t183 + t232 * t84;
t26 = -t167 * t54 + t170 * t52;
t27 = t167 * t52 + t170 * t54;
t185 = -t241 * t271 - t266;
t184 = -t102 * t182 - t14;
t51 = pkin(4) * t182 + pkin(8) * t87 + t264;
t181 = t144 * t123 + t197;
t140 = t169 * t207;
t142 = t172 * t207;
t65 = t171 * t140 + t168 * t142 + t145 * t216 + t146 * t217;
t111 = pkin(8) + t119;
t180 = -t111 * t58 - t117 * t271 + t39 * t87;
t118 = -pkin(2) * t228 + t232 * t159;
t178 = qJD(5) * t188 + t1 - t261;
t66 = -t104 * qJD(3) - t168 * t140 + t171 * t142;
t176 = t100 * qJ(4) - t135 * qJD(4) + t66;
t175 = t102 * t87 - t15;
t157 = -t232 * pkin(3) - pkin(4);
t156 = t166 * pkin(3) + pkin(8);
t110 = -pkin(4) - t118;
t105 = t161 + t264;
t80 = -t123 ^ 2 + t125 ^ 2;
t68 = -t220 + (t123 - t206) * t163;
t62 = -t166 * t100 + t232 * t101;
t50 = t161 + t51;
t47 = -t101 * qJ(4) - t134 * qJD(4) + t65;
t44 = t232 * t74 - t246;
t43 = t166 * t74 + t70;
t38 = t59 + t235;
t37 = -t58 + t236;
t29 = t252 - t253;
t28 = t62 * pkin(4) - t63 * pkin(8) + t90;
t23 = t167 * t50 + t170 * t49;
t22 = -t167 * t49 + t170 * t50;
t21 = t167 * t51 + t170 * t44;
t20 = -t167 * t44 + t170 * t51;
t19 = t166 * t176 + t232 * t47;
t18 = t166 * t47 - t232 * t176;
t10 = t201 * t76 - t237;
t9 = t200 * t78 - t238;
t8 = -t257 - t270;
t7 = t185 + t255;
t6 = (-t31 - t259) * t170 - t271 * t242 + t249;
t5 = -qJD(5) * t27 - t167 * t19 + t170 * t28;
t4 = qJD(5) * t26 + t167 * t28 + t170 * t19;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t194, t219 * t268, t222, -0.2e1 * t194, -t223, 0, -pkin(6) * t222 + t169 * t199, pkin(6) * t223 + t172 * t199, 0, 0, -t125 * t100 - t91 * t135, t100 * t123 - t125 * t101 + t91 * t134 - t135 * t92, -t100 * t163, t123 * t101 + t92 * t134, -t101 * t163, 0, t144 * t101 + t160 * t92 + t66 * t163 + (qJD(1) * t134 + t123) * t162, -t144 * t100 + t125 * t269 - t160 * t91 - t65 * t163, t93 * t100 - t94 * t101 + t103 * t91 - t104 * t92 - t65 * t123 - t66 * t125 + t134 * t197 - t61 * t135, t61 * t103 - t104 * t197 + t144 * t269 + t94 * t65 + t93 * t66, t182 * t63 + t59 * t97, -t182 * t62 - t59 * t96 - t63 * t87 - t250, t63 * t163, t87 * t62 + t260, -t62 * t163, 0, t102 * t62 + t106 * t58 - t18 * t163 + t81 * t96 + t90 * t87, t102 * t63 + t106 * t59 - t19 * t163 + t182 * t90 + t81 * t97, -t15 * t96 + t18 * t182 - t19 * t87 - t41 * t63 - t42 * t62 + t53 * t59 - t54 * t58 + t262, t102 * t90 + t81 * t106 + t15 * t54 - t41 * t18 + t42 * t19 + t263, -t78 * t210 + (-t31 * t97 + t63 * t78) * t170, (-t170 * t76 - t242) * t63 + (t238 - t237 + (-t240 + t243) * qJD(5)) * t97, t170 * t191 - t210 * t271 - t31 * t96 + t78 * t62, t76 * t209 + (t32 * t97 + t63 * t76) * t167, -t167 * t191 - t209 * t271 - t32 * t96 - t76 * t62, t271 * t62 + t260, t167 * t193 + t18 * t76 - t186 * t62 + t26 * t58 + t271 * t5 + t3 * t96 + t53 * t32 + t36 * t97, -t17 * t62 + t170 * t193 + t18 * t78 - t2 * t96 - t27 * t58 - t271 * t4 - t53 * t31 - t35 * t97, t26 * t31 - t27 * t32 - t4 * t76 - t5 * t78 + t188 * t63 + (qJD(5) * t187 - t2 * t167 - t3 * t170) * t97, t17 * t4 + t39 * t18 - t186 * t5 + t2 * t27 + t3 * t26 + t263; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t208, t219 * t174, 0, t208, 0, 0, t174 * pkin(1) * t169, pkin(1) * t224, 0, 0, t230, t80, t68, -t230, 0, 0, -t123 * t161 - t229 - t98 * t163 + (t226 + (-pkin(2) * t163 - t130) * t168) * qJD(3) + t198, t99 * t163 + (-t125 * t218 - t163 * t216) * pkin(2) + t181, (t94 + t98) * t125 + (-t93 + t99) * t123 + (-t168 * t92 + t171 * t91 + (-t123 * t171 + t125 * t168) * qJD(3)) * pkin(2), -t93 * t98 - t94 * t99 + (-t144 * t218 - t168 * t197 + t171 * t61 + (-t168 * t93 + t171 * t94) * qJD(3)) * pkin(2), t272, t29, t38, -t272, t37, 0, -t105 * t87 - t234 * t163 + t184, -t105 * t182 - t233 * t163 + t175, -t118 * t59 - t119 * t58 + t182 * t234 - t233 * t87 + t192, -t102 * t105 - t14 * t118 + t15 * t119 + t233 * t42 - t234 * t41, t9, t6, t8, t10, t7, -t256, t110 * t32 - t22 * t271 + t234 * t76 + (-qJD(5) * t111 * t271 - t14) * t170 + t180 * t167 + t204, -t110 * t31 + (t111 * t215 + t23) * t271 + t234 * t78 + t180 * t170 + t196, t22 * t78 + t23 * t76 + (-t111 * t32 - t117 * t76 + (t111 * t78 + t186) * qJD(5)) * t170 + (-t111 * t31 + t117 * t78 - t3 + (t111 * t76 - t17) * qJD(5)) * t167 + t212, t14 * t110 + t111 * t178 - t117 * t187 - t17 * t23 + t186 * t22 + t234 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t230, t80, t68, -t230, 0, 0, t94 * t163 - t229 + t61, t93 * t163 + t181, 0, 0, t272, t29, t38, -t272, t37, 0, t43 * t163 - t87 * t264 + t184, t44 * t163 - t182 * t264 + t175, -t43 * t182 + t44 * t87 + (-t166 * t58 - t232 * t59) * pkin(3) + t192, t41 * t43 - t42 * t44 + (-t102 * t125 - t14 * t232 + t15 * t166) * pkin(3), t9, t6, t8, t10, t7, -t256, t39 * t241 - t14 * t170 + t157 * t32 - t20 * t271 - t43 * t76 + (-t214 * t271 - t245) * t156 + t204, t266 * t156 - t157 * t31 + t21 * t271 + t39 * t239 - t43 * t78 + t196, -t261 + t20 * t78 + t21 * t76 + (-t237 - t238) * t156 + ((t240 + t243) * t156 + t188) * qJD(5) + t212, t14 * t157 + t156 * t178 - t17 * t21 + t186 * t20 - t39 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 + t236, t59 - t235, -t252 - t253, t182 * t41 + t42 * t87 + t81, 0, 0, 0, 0, 0, 0, t185 - t255, -t257 + t270, (t31 - t259) * t170 + t78 * t201 + t249, t273 * t167 + t274 * t170 - t39 * t182; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t258, -t76 ^ 2 + t78 ^ 2, t271 * t76 - t31, -t258, -t244 + (-qJD(5) + t271) * t78, t58, -t39 * t78 + t274, t39 * t76 - t273, 0, 0;];
tauc_reg = t11;
