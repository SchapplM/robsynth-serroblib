% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRR15_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:43:09
% EndTime: 2019-12-31 20:43:18
% DurationCPUTime: 3.08s
% Computational Cost: add. (4428->382), mult. (9911->527), div. (0->0), fcn. (5972->6), ass. (0->203)
t246 = pkin(3) + pkin(6);
t153 = cos(qJ(2));
t149 = sin(qJ(4));
t150 = sin(qJ(2));
t217 = t149 * t150;
t167 = pkin(4) * t153 - pkin(8) * t217;
t205 = qJD(4) * t149;
t154 = -pkin(2) - pkin(7);
t244 = pkin(8) - t154;
t208 = qJD(1) * t153;
t138 = pkin(6) * t208;
t108 = pkin(3) * t208 + t138;
t152 = cos(qJ(4));
t199 = t150 * qJD(1);
t141 = pkin(2) * t199;
t173 = pkin(7) * t150 - qJ(3) * t153;
t77 = qJD(1) * t173 + t141;
t40 = t152 * t108 - t149 * t77;
t263 = qJD(1) * t167 - t244 * t205 + t40;
t113 = t244 * t152;
t191 = t152 * t199;
t41 = t149 * t108 + t152 * t77;
t262 = pkin(8) * t191 + qJD(4) * t113 + t41;
t148 = sin(qJ(5));
t151 = cos(qJ(5));
t101 = t148 * t152 + t151 * t149;
t163 = t101 * t150;
t252 = qJD(4) + qJD(5);
t50 = t252 * t101;
t241 = qJD(1) * t163 + t50;
t187 = t149 * t208;
t197 = t152 * qJD(2);
t100 = -t187 + t197;
t200 = t149 * qJD(2);
t98 = t152 * t208 + t200;
t172 = t151 * t100 - t148 * t98;
t43 = t148 * t100 + t151 * t98;
t245 = t43 * t172;
t261 = t172 ^ 2 - t43 ^ 2;
t133 = qJD(4) + t199;
t125 = qJD(5) + t133;
t136 = qJD(4) * t197;
t203 = qJD(4) * t153;
t190 = t149 * t203;
t162 = t150 * t197 + t190;
t158 = qJD(1) * t162 - t136;
t201 = qJD(5) * t151;
t202 = qJD(5) * t148;
t195 = qJD(1) * qJD(2);
t183 = t150 * t195;
t122 = t149 * t183;
t61 = qJD(4) * t98 - t122;
t18 = t100 * t202 - t148 * t158 + t151 * t61 + t98 * t201;
t260 = t43 * t125 - t18;
t204 = qJD(4) * t152;
t132 = pkin(2) * t183;
t198 = t150 * qJD(3);
t160 = qJD(2) * t173 - t198;
t56 = qJD(1) * t160 + t132;
t216 = t150 * qJ(3);
t182 = -pkin(1) - t216;
t254 = t153 * t154;
t94 = t182 + t254;
t67 = t94 * qJD(1);
t137 = pkin(6) * t199;
t253 = qJD(3) + t137;
t196 = pkin(3) * t199 + t253;
t72 = t154 * qJD(2) + t196;
t135 = t153 * t195;
t131 = pkin(6) * t135;
t93 = pkin(3) * t135 + t131;
t14 = t149 * t93 + t152 * t56 + t72 * t204 - t205 * t67;
t11 = pkin(8) * t158 + t14;
t33 = -t149 * t67 + t152 * t72;
t28 = -t100 * pkin(8) + t33;
t24 = t133 * pkin(4) + t28;
t34 = t149 * t72 + t152 * t67;
t29 = -t98 * pkin(8) + t34;
t15 = -qJD(4) * t34 - t149 * t56 + t152 * t93;
t8 = pkin(4) * t135 + t61 * pkin(8) + t15;
t1 = (qJD(5) * t24 + t11) * t151 + t148 * t8 - t29 * t202;
t145 = qJD(2) * qJ(3);
t87 = t145 + t108;
t52 = t98 * pkin(4) + t87;
t259 = t52 * t43 - t1;
t257 = -0.2e1 * t195;
t256 = -t33 * t133 + t14;
t255 = t34 * t133 + t15;
t120 = t246 * t150;
t103 = t149 * t120;
t49 = t152 * t94 + t103;
t214 = t151 * t152;
t218 = t148 * t149;
t170 = -t214 + t218;
t78 = t170 * t153;
t146 = t150 ^ 2;
t147 = t153 ^ 2;
t209 = t146 - t147;
t232 = t151 * t29;
t6 = t148 * t24 + t232;
t2 = -qJD(5) * t6 - t148 * t11 + t151 * t8;
t251 = -t52 * t172 + t2;
t19 = qJD(5) * t172 - t148 * t61 - t151 * t158;
t240 = -t148 * t205 - t149 * t202 + t151 * t191 - t199 * t218 + t252 * t214;
t250 = t101 * t19 + t240 * t43;
t249 = t125 * t172 - t19;
t248 = t170 * t18 - t172 * t241;
t235 = t148 * t29;
t5 = t151 * t24 - t235;
t247 = t1 * t101 - t170 * t2 + t240 * t6 - t241 * t5;
t112 = t244 * t149;
t54 = -t151 * t112 - t148 * t113;
t243 = qJD(5) * t54 - t262 * t148 + t263 * t151;
t53 = t148 * t112 - t151 * t113;
t242 = -qJD(5) * t53 + t263 * t148 + t262 * t151;
t239 = qJD(2) * pkin(2);
t238 = t100 * t98;
t234 = t149 * t33;
t233 = t149 * t98;
t231 = t152 * t61;
t228 = t61 * t149;
t207 = qJD(2) * t150;
t107 = t246 * t207;
t144 = qJD(2) * qJD(3);
t74 = -qJD(1) * t107 + t144;
t227 = t74 * t149;
t226 = t74 * t152;
t225 = t87 * t150;
t224 = t98 * t152;
t192 = -pkin(4) * t152 - pkin(3);
t223 = pkin(4) * t204 - t192 * t199 + t253;
t222 = t100 * t153;
t221 = t133 * t150;
t220 = t133 * t154;
t219 = t136 * t149;
t215 = t150 * t152;
t213 = t152 * t153;
t156 = qJD(1) ^ 2;
t212 = t153 * t156;
t155 = qJD(2) ^ 2;
t211 = t155 * t150;
t210 = t155 * t153;
t121 = t246 * t153;
t114 = -t153 * pkin(2) + t182;
t88 = qJD(1) * t114;
t206 = qJD(2) * t153;
t194 = t98 * t215;
t193 = t133 * t215;
t189 = t133 * t204;
t188 = t152 * t203;
t186 = pkin(8) * t153 - t94;
t184 = t240 * t125;
t109 = t246 * t206;
t140 = pkin(2) * t207;
t63 = t140 + t160;
t180 = t152 * t109 - t149 * t63;
t179 = qJD(1) * t49 + t34;
t177 = pkin(1) * t257;
t176 = qJD(3) - t239;
t175 = -t241 * t125 - t170 * t135;
t104 = t152 * t120;
t37 = t150 * pkin(4) + t149 * t186 + t104;
t39 = -pkin(8) * t213 + t49;
t20 = -t148 * t39 + t151 * t37;
t21 = t148 * t37 + t151 * t39;
t171 = t100 * t152 - t233;
t169 = -0.2e1 * qJD(2) * t88;
t168 = t133 * t149;
t164 = -qJ(3) * t206 - t198;
t65 = qJD(1) * t164 + t132;
t82 = t140 + t164;
t166 = pkin(6) * t155 + qJD(1) * t82 + t65;
t165 = t152 * t183 - t136;
t22 = t149 * t109 + t120 * t204 + t152 * t63 - t205 * t94;
t55 = -pkin(4) * t190 + (-pkin(6) + t192) * t207;
t110 = pkin(6) * t183 - t144;
t111 = t137 + t176;
t117 = -t138 - t145;
t157 = -t110 * t153 + (t111 * t153 + (t117 + t138) * t150) * qJD(2);
t134 = t149 * pkin(4) + qJ(3);
t129 = t150 * t212;
t124 = t152 * t135;
t123 = t150 * t135;
t119 = -0.2e1 * t123;
t118 = 0.2e1 * t123;
t115 = t209 * t156;
t105 = -qJ(3) * t208 + t141;
t86 = t209 * t257;
t85 = pkin(4) * t213 + t121;
t79 = t101 * t153;
t71 = t88 * t199;
t48 = -t149 * t94 + t104;
t36 = t136 * pkin(4) + qJD(1) * t55 + t144;
t31 = t153 * t50 - t170 * t207;
t30 = qJD(2) * t163 + t252 * t78;
t23 = -t49 * qJD(4) + t180;
t17 = pkin(8) * t162 + t22;
t16 = t167 * qJD(2) + (t152 * t186 - t103) * qJD(4) + t180;
t10 = t151 * t28 - t235;
t9 = -t148 * t28 - t232;
t4 = -qJD(5) * t21 - t148 * t17 + t151 * t16;
t3 = qJD(5) * t20 + t148 * t16 + t151 * t17;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t118, t86, t210, t119, -t211, 0, -pkin(6) * t210 + t150 * t177, pkin(6) * t211 + t153 * t177, 0, 0, 0, -t210, t211, t118, t86, t119, t157, t150 * t169 + t153 * t166, -t150 * t166 + t153 * t169, pkin(6) * t157 + t65 * t114 + t88 * t82, t153 * t228 + (t150 * t200 - t188) * t100, t171 * t207 + (-t165 * t149 + t231 + (t224 + (t100 - t187) * t149) * qJD(4)) * t153, -t133 * t188 - t61 * t150 + (t222 + (-qJD(1) * t147 + t221) * t149) * qJD(2), -t158 * t213 - t162 * t98, t133 * t190 + (qJD(4) * t187 - t136) * t150 + (-t98 * t153 + (t209 * qJD(1) + t221) * t152) * qJD(2), t133 * t206 + t123, -t107 * t98 + t121 * t136 + t23 * t133 + (t15 + (-qJD(1) * t121 - t87) * t197) * t150 + (-t87 * t205 + t33 * qJD(2) + t226 + (qJD(2) * t48 - t121 * t205) * qJD(1)) * t153, -t107 * t100 - t121 * t61 - t22 * t133 + (t200 * t87 - t14) * t150 + (-qJD(2) * t179 - t204 * t87 - t227) * t153, -t23 * t100 - t49 * t136 - t22 * t98 + t48 * t61 + (t152 * t179 - t234) * t207 + (-t14 * t152 + t15 * t149 + (t149 * t179 + t152 * t33) * qJD(4)) * t153, -t87 * t107 + t74 * t121 + t14 * t49 + t15 * t48 + t34 * t22 + t33 * t23, t172 * t30 + t18 * t79, t172 * t31 - t18 * t78 + t79 * t19 - t30 * t43, t30 * t125 - t18 * t150 + (-qJD(1) * t79 + t172) * t206, -t19 * t78 - t43 * t31, t31 * t125 - t19 * t150 + (qJD(1) * t78 - t43) * t206, t125 * t206 + t123, t4 * t125 + t2 * t150 + t85 * t19 - t52 * t31 - t36 * t78 + t55 * t43 + (qJD(1) * t20 + t5) * t206, -t1 * t150 - t3 * t125 - t85 * t18 + t52 * t30 - t36 * t79 + t55 * t172 + (-qJD(1) * t21 - t6) * t206, t1 * t78 - t172 * t4 + t20 * t18 - t21 * t19 + t2 * t79 - t3 * t43 - t5 * t30 + t6 * t31, t1 * t21 + t2 * t20 + t6 * t3 + t36 * t85 + t5 * t4 + t52 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t129, t115, 0, t129, 0, 0, t156 * pkin(1) * t150, pkin(1) * t212, 0, 0, 0, 0, 0, -t129, t115, t129, ((-t117 - t145) * t150 + (-t111 + t176) * t153) * qJD(1), -t105 * t208 + t71, 0.2e1 * t144 + (t105 * t150 + t153 * t88) * qJD(1), -t110 * qJ(3) - t117 * qJD(3) - t88 * t105 + (-t117 * t150 + (-t111 - t239) * t153) * qJD(1) * pkin(6), -t100 * t168 - t231, -t136 * t152 + t228 - t171 * qJD(4) + (t149 * t188 + (t233 + (-t100 + t197) * t152) * t150) * qJD(1), -t133 * t205 + t124 + (-t133 * t217 - t222) * qJD(1), t98 * t204 + t219 + (-t149 * t162 + t194) * qJD(1), -t189 + (-t193 + (t98 - t200) * t153) * qJD(1), -t133 * t208, qJ(3) * t136 - t40 * t133 + t227 + t196 * t98 + (-t149 * t220 + t152 * t87) * qJD(4) + ((-qJ(3) * t205 - t33) * t153 + (t225 + (-t216 + t254) * qJD(2)) * t152) * qJD(1), -qJ(3) * t61 + t41 * t133 + t226 + t196 * t100 + (-t149 * t87 - t152 * t220) * qJD(4) + (t153 * t34 + (-t154 * t206 - t225) * t149) * qJD(1), t40 * t100 + t41 * t98 + (-t34 * t199 + t154 * t61 - t15 + (-t154 * t98 - t34) * qJD(4)) * t152 + (t154 * t165 - t14 + t33 * t199 + (t33 + (t100 + t187) * t154) * qJD(4)) * t149, t74 * qJ(3) - t33 * t40 - t34 * t41 + t196 * t87 + (t14 * t149 + t15 * t152 + (t152 * t34 - t234) * qJD(4)) * t154, t248, t18 * t101 + t170 * t19 - t172 * t240 + t241 * t43, -t172 * t208 + t175, t250, -t184 + (-qJD(2) * t101 + t43) * t208, -t125 * t208, t36 * t101 + t134 * t19 + t240 * t52 + t223 * t43 - t243 * t125 + (qJD(2) * t53 - t5) * t208, -t36 * t170 - t134 * t18 - t241 * t52 + t223 * t172 + t242 * t125 + (-qJD(2) * t54 + t6) * t208, t172 * t243 + t53 * t18 - t54 * t19 + t242 * t43 - t247, t1 * t54 + t36 * t134 + t2 * t53 + t223 * t52 - t242 * t6 - t243 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t129, -t146 * t156 - t155, t117 * qJD(2) + t131 + t71, 0, 0, 0, 0, 0, 0, -qJD(2) * t98 - t133 * t168 + t124, -t189 - qJD(2) * t100 + (-t153 * t200 - t193) * qJD(1), -t219 + t231 + (t100 * t149 - t224) * qJD(4) + (-t194 + (t190 + (t100 + t197) * t150) * t149) * qJD(1), -t87 * qJD(2) + t256 * t149 + t255 * t152, 0, 0, 0, 0, 0, 0, -qJD(2) * t43 + t175, -t184 + (-t101 * t208 - t172) * qJD(2), -t248 - t250, -t52 * qJD(2) + t247; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t238, t100 ^ 2 - t98 ^ 2, t122 + (-qJD(4) + t133) * t98, -t238, t100 * t133 + t158, t135, -t87 * t100 + t255, t87 * t98 - t256, 0, 0, t245, t261, t260, -t245, t249, t135, -t9 * t125 + (-t100 * t43 - t125 * t202 + t135 * t151) * pkin(4) + t251, t10 * t125 + (-t100 * t172 - t125 * t201 - t135 * t148) * pkin(4) + t259, t10 * t43 + t6 * t172 - t5 * t43 + t9 * t172 + (-t148 * t19 + t151 * t18 + (t148 * t172 - t151 * t43) * qJD(5)) * pkin(4), -t6 * t10 - t5 * t9 + (t1 * t148 - t100 * t52 + t151 * t2 + (-t148 * t5 + t151 * t6) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t245, t261, t260, -t245, t249, t135, t6 * t125 + t251, t5 * t125 + t259, 0, 0;];
tauc_reg = t7;
