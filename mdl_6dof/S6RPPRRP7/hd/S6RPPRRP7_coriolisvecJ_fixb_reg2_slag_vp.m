% Calculate inertial parameters regressor of coriolis joint torque vector for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% 
% Output:
% tauc_reg [6x(6*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRRP7_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:14:04
% EndTime: 2019-03-09 02:14:07
% DurationCPUTime: 2.88s
% Computational Cost: add. (5313->347), mult. (11709->442), div. (0->0), fcn. (8216->6), ass. (0->188)
t133 = sin(pkin(9));
t226 = sin(qJ(4));
t176 = t226 * t133;
t134 = cos(pkin(9));
t138 = cos(qJ(4));
t189 = t134 * t138;
t107 = -t176 + t189;
t136 = sin(qJ(5));
t137 = cos(qJ(5));
t135 = -pkin(1) - qJ(3);
t235 = t135 * qJD(1);
t117 = qJD(2) + t235;
t166 = -pkin(7) * qJD(1) + t117;
t92 = t166 * t133;
t93 = t166 * t134;
t67 = t138 * t92 + t226 * t93;
t63 = qJD(4) * pkin(8) + t67;
t175 = qJD(1) * t189;
t101 = -qJD(1) * t176 + t175;
t132 = qJD(1) * qJ(2);
t126 = qJD(3) + t132;
t128 = t133 * pkin(3);
t110 = qJD(1) * t128 + t126;
t106 = t138 * t133 + t226 * t134;
t143 = qJD(1) * t106;
t64 = pkin(4) * t143 - pkin(8) * t101 + t110;
t33 = t136 * t64 + t137 * t63;
t181 = t137 * qJD(4);
t84 = t101 * t136 - t181;
t21 = -qJ(6) * t84 + t33;
t239 = qJD(5) + t143;
t242 = t21 * t239;
t144 = t106 * qJD(3);
t236 = t138 * t93 - t226 * t92;
t46 = -qJD(1) * t144 + qJD(4) * t236;
t131 = qJD(1) * qJD(2);
t172 = qJD(4) * t226;
t159 = qJD(1) * t172;
t185 = qJD(4) * t138;
t171 = qJD(1) * t185;
t187 = t133 * t171 + t134 * t159;
t114 = t133 * t159;
t91 = t134 * t171 - t114;
t65 = t91 * pkin(4) + t187 * pkin(8) + t131;
t12 = -qJD(5) * t33 - t136 * t46 + t137 * t65;
t232 = t239 * t33 + t12;
t183 = qJD(5) * t137;
t184 = qJD(5) * t136;
t167 = -t136 * t65 - t137 * t46 - t64 * t183 + t63 * t184;
t32 = -t136 * t63 + t137 * t64;
t157 = -t239 * t32 - t167;
t241 = t107 * qJD(3);
t163 = t137 * t239;
t210 = t136 * t91;
t240 = -t163 * t239 - t210;
t237 = t137 * t91 - t184 * t239;
t86 = qJD(4) * t136 + t101 * t137;
t194 = qJD(5) * t86;
t55 = -t136 * t187 + t194;
t220 = -pkin(7) + t135;
t108 = t220 * t133;
t109 = t220 * t134;
t80 = t226 * t108 - t138 * t109;
t186 = t133 ^ 2 + t134 ^ 2;
t234 = t186 * qJD(3);
t54 = -qJD(5) * t181 + t101 * t184 + t137 * t187;
t140 = t54 * qJ(6) + t12;
t227 = t91 * pkin(5);
t2 = -t86 * qJD(6) + t140 + t227;
t233 = t2 + t242;
t103 = t106 * qJD(4);
t231 = -t101 * t103 - t107 * t187;
t230 = t86 ^ 2;
t229 = t101 ^ 2;
t179 = 0.2e1 * t131;
t228 = pkin(5) * t84;
t148 = t241 * qJD(1);
t47 = t67 * qJD(4) + t148;
t225 = t47 * t80;
t224 = t84 * t239;
t223 = t84 * t143;
t222 = t86 * t84;
t221 = t86 * t239;
t219 = -qJ(6) - pkin(8);
t20 = -qJ(6) * t86 + t32;
t19 = pkin(5) * t239 + t20;
t218 = t19 - t20;
t168 = qJD(5) * t219;
t196 = qJ(6) * t137;
t78 = pkin(4) * t101 + pkin(8) * t143;
t38 = -t136 * t236 + t137 * t78;
t217 = pkin(5) * t101 + t136 * qJD(6) - t137 * t168 + t143 * t196 + t38;
t209 = t136 * t143;
t39 = t136 * t78 + t137 * t236;
t216 = qJ(6) * t209 - t137 * qJD(6) - t136 * t168 + t39;
t215 = -t136 * t55 - t84 * t183;
t121 = qJ(2) + t128;
t77 = pkin(4) * t106 - pkin(8) * t107 + t121;
t81 = t138 * t108 + t226 * t109;
t79 = t137 * t81;
t43 = t136 * t77 + t79;
t214 = t101 * t84;
t213 = t101 * t143;
t212 = t136 * t84;
t211 = t136 * t86;
t170 = -qJD(6) - t228;
t62 = -qJD(4) * pkin(4) - t236;
t41 = -t170 + t62;
t208 = t137 * t41;
t207 = t137 * t84;
t206 = t137 * t86;
t29 = t55 * pkin(5) + t47;
t204 = t29 * t136;
t203 = t29 * t137;
t202 = t47 * t107;
t201 = t54 * t136;
t200 = t55 * t137;
t199 = t86 * t101;
t198 = t91 * t106;
t197 = t239 * t101;
t195 = qJD(4) * t143;
t192 = t103 * t136;
t191 = t103 * t137;
t190 = t107 * t136;
t104 = -t133 * t172 + t134 * t185;
t182 = t104 * qJD(4);
t56 = -t80 * qJD(4) - t144;
t75 = pkin(4) * t104 + pkin(8) * t103 + qJD(2);
t180 = t136 * t75 + t137 * t56 + t77 * t183;
t174 = t107 * t183;
t169 = -t136 * t56 + t137 * t75;
t42 = -t136 * t81 + t137 * t77;
t165 = t136 * t239;
t162 = qJD(1) * t186;
t161 = qJD(5) * t106 + qJD(1);
t160 = pkin(8) * qJD(5) * t239 + t47;
t147 = t55 * qJ(6) + t167;
t4 = -qJD(6) * t84 - t147;
t158 = -t19 * t239 + t4;
t156 = t104 * t143 + t198;
t155 = t136 * t21 + t137 * t19;
t154 = t136 * t19 - t137 * t21;
t153 = t136 * t33 + t137 * t32;
t152 = t136 * t32 - t137 * t33;
t151 = t206 + t212;
t150 = qJ(6) * t103 - qJD(6) * t107;
t149 = -t209 * t239 + t237;
t145 = -pkin(8) * t91 + t239 * t62;
t142 = t103 * t236 - t67 * t104 - t46 * t106 + t202;
t141 = -t153 * qJD(5) - t12 * t136 - t137 * t167;
t57 = t81 * qJD(4) + t241;
t139 = qJD(1) ^ 2;
t125 = -pkin(5) * t137 - pkin(4);
t112 = t219 * t137;
t111 = t219 * t136;
t98 = t143 ^ 2;
t94 = t103 * qJD(4);
t83 = t84 ^ 2;
t60 = pkin(5) * t190 + t80;
t48 = -pkin(5) * t209 + t67;
t45 = t104 * t239 + t198;
t40 = -t83 + t230;
t37 = (t174 - t192) * pkin(5) + t57;
t36 = -qJ(6) * t190 + t43;
t35 = t221 - t55;
t34 = -t54 + t224;
t31 = pkin(5) * t106 - t107 * t196 + t42;
t28 = -t199 + t240;
t27 = -t199 - t240;
t26 = t149 + t214;
t25 = t149 - t214;
t23 = t84 * t165 - t200;
t22 = t86 * t163 - t201;
t18 = -t215 * t107 - t84 * t192;
t17 = -t86 * t191 + (-t137 * t54 - t86 * t184) * t107;
t16 = -t43 * qJD(5) + t169;
t15 = -t81 * t184 + t180;
t14 = t239 * t192 - t84 * t104 - t55 * t106 + (-t183 * t239 - t210) * t107;
t13 = t86 * t104 - t54 * t106 + t237 * t107 - t191 * t239;
t10 = -qJ(6) * t174 + (-qJD(5) * t81 + t150) * t136 + t180;
t9 = -t136 * t198 + t103 * t84 - t107 * t55 + (-t104 * t136 - t161 * t137) * t239;
t8 = -t137 * t198 + t103 * t86 + t107 * t54 + (-t104 * t137 + t161 * t136) * t239;
t7 = t104 * pkin(5) + t150 * t137 + (-t79 + (qJ(6) * t107 - t77) * t136) * qJD(5) + t169;
t6 = (-t54 - t223) * t137 - t239 * t211 + t215;
t5 = (t54 - t223) * t137 + t86 * t165 + t215;
t3 = (t207 + t211) * t103 + (t201 - t200 + (-t206 + t212) * qJD(5)) * t107;
t1 = (-t207 + t211) * t104 + t151 * qJD(1) + (t151 * qJD(5) - t200 - t201) * t106;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t179, qJ(2) * t179, 0, 0, 0, 0, 0, 0, t133 * t179, t134 * t179, 0.2e1 * qJD(3) * t162 (t126 + t132) * qJD(2) + (-t117 - t235) * t234, t231, -t101 * t104 + t103 * t143 + t187 * t106 - t107 * t91, -t94, t156, -t182, 0, 0.2e1 * t143 * qJD(2) - t57 * qJD(4) + t110 * t104 + t121 * t91, -t121 * t187 - t110 * t103 - t56 * qJD(4) + (qJD(1) * t107 + t101) * qJD(2), t57 * t101 - t143 * t56 - t80 * t187 - t81 * t91 + t142, t46 * t81 + t225 + t67 * t56 - t236 * t57 + (qJD(1) * t121 + t110) * qJD(2), t17, t3, t13, t18, t14, t45, -t62 * t192 + t32 * t104 + t12 * t106 + t16 * t239 + t42 * t91 + t80 * t55 + t57 * t84 + (t47 * t136 + t183 * t62) * t107, -t62 * t191 - t33 * t104 + t167 * t106 - t15 * t239 - t43 * t91 - t80 * t54 + t57 * t86 + (t47 * t137 - t184 * t62) * t107, -t15 * t84 - t16 * t86 + t42 * t54 - t43 * t55 + t153 * t103 + (qJD(5) * t152 - t12 * t137 + t136 * t167) * t107, t12 * t42 + t15 * t33 + t16 * t32 - t167 * t43 + t57 * t62 + t225, t17, t3, t13, t18, t14, t45, -t41 * t192 + t19 * t104 + t2 * t106 + t31 * t91 + t37 * t84 + t60 * t55 + t7 * t239 + (t183 * t41 + t204) * t107, -t41 * t191 - t10 * t239 - t21 * t104 - t4 * t106 - t36 * t91 + t37 * t86 - t60 * t54 + (-t184 * t41 + t203) * t107, -t10 * t84 + t31 * t54 - t36 * t55 - t7 * t86 + t155 * t103 + (qJD(5) * t154 - t4 * t136 - t2 * t137) * t107, t10 * t21 + t19 * t7 + t2 * t31 + t29 * t60 + t36 * t4 + t37 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t139, -t139 * qJ(2), 0, 0, 0, 0, 0, 0, -t139 * t133, -t139 * t134, 0 (-t126 - t234) * qJD(1), 0, 0, 0, 0, 0, 0, -qJD(1) * t143 - t94, -qJD(1) * t101 - t182, -t156 - t231, -qJD(1) * t110 - t142, 0, 0, 0, 0, 0, 0, t9, t8, t1, -qJD(1) * t153 + t62 * t103 - t104 * t152 + t106 * t141 - t202, 0, 0, 0, 0, 0, 0, t9, t8, t1, t41 * t103 - t29 * t107 - t154 * t104 - t155 * qJD(1) + (-qJD(5) * t155 - t2 * t136 + t4 * t137) * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t186 * t139, t117 * t162 + t131, 0, 0, 0, 0, 0, 0, -t114 + (t101 + t175) * qJD(4), -t187 - t195, -t98 - t229, t101 * t236 + t143 * t67 + t131, 0, 0, 0, 0, 0, 0, t25, t28, t5, -t62 * t101 + t157 * t136 + t232 * t137, 0, 0, 0, 0, 0, 0, t25, t28, t5, -t41 * t101 + t158 * t136 + t233 * t137; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t213, -t98 + t229, -t187 + t195, -t213, t114 + (t101 - t175) * qJD(4), 0, -t110 * t101 - t148 (qJD(3) + t110) * t143, 0, 0, t22, t6, t27, t23, t26, -t197, -pkin(4) * t55 - t32 * t101 + t136 * t145 - t137 * t160 - t239 * t38 - t67 * t84, pkin(4) * t54 + t33 * t101 + t136 * t160 + t137 * t145 + t239 * t39 - t67 * t86, t38 * t86 + t39 * t84 + ((-t55 + t194) * pkin(8) + t157) * t137 + ((qJD(5) * t84 - t54) * pkin(8) - t232) * t136, -t47 * pkin(4) + pkin(8) * t141 - t32 * t38 - t33 * t39 - t62 * t67, t22, t6, t27, t23, t26, -t197, -t19 * t101 + t111 * t91 + t125 * t55 - t203 - t48 * t84 - t217 * t239 + (t41 * t143 + (t41 + t228) * qJD(5)) * t136, t143 * t208 + t21 * t101 + t112 * t91 - t125 * t54 + t204 - t48 * t86 + t216 * t239 + (pkin(5) * t211 + t208) * qJD(5), t111 * t54 + t112 * t55 - t233 * t136 + t158 * t137 + t216 * t84 + t217 * t86, t2 * t111 - t4 * t112 + t29 * t125 + (pkin(5) * t184 - t48) * t41 - t216 * t21 - t217 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t222, t40, t34, -t222, t35, t91, -t62 * t86 + t232, t62 * t84 - t157, 0, 0, t222, t40, t34, -t222, t35, t91, 0.2e1 * t227 + t242 + (t170 - t41) * t86 + t140, -t230 * pkin(5) + t20 * t239 + (qJD(6) + t41) * t84 + t147, t54 * pkin(5) - t218 * t84, t218 * t21 + (-t41 * t86 + t2) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55 + t221, -t54 - t224, -t83 - t230, t19 * t86 + t21 * t84 + t29;];
tauc_reg  = t11;
