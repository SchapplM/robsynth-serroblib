% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRRPP3
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
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:14
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPP3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP3_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:13:41
% EndTime: 2019-12-05 16:13:49
% DurationCPUTime: 3.27s
% Computational Cost: add. (1851->390), mult. (4258->513), div. (0->0), fcn. (2956->8), ass. (0->197)
t119 = sin(pkin(8));
t125 = cos(qJ(3));
t123 = sin(qJ(3));
t191 = t123 * qJDD(2);
t175 = t119 * t191;
t121 = cos(pkin(8));
t197 = t121 * qJD(3);
t204 = qJD(2) * t123;
t80 = t119 * t204 - t197;
t224 = t123 * t80;
t117 = t123 ^ 2;
t118 = t125 ^ 2;
t209 = t117 - t118;
t249 = qJD(2) * t209;
t112 = t121 * qJDD(3);
t194 = qJD(2) * qJD(3);
t176 = t125 * t194;
t140 = t176 + t191;
t41 = t119 * t140 - t112;
t265 = (t119 * t249 + t224) * qJD(3) - (t41 + t175) * t125;
t124 = sin(qJ(2));
t126 = cos(qJ(2));
t208 = t117 + t118;
t171 = t126 * t208;
t206 = qJD(1) * t124;
t96 = qJD(2) * pkin(6) + t206;
t205 = qJD(1) * t126;
t230 = qJD(2) * pkin(2);
t97 = -t205 - t230;
t264 = t97 * t124 + t171 * t96;
t120 = sin(pkin(7));
t122 = cos(pkin(7));
t160 = g(1) * t122 + g(2) * t120;
t201 = qJD(3) * t119;
t82 = t121 * t204 + t201;
t254 = (t119 * t82 + t121 * t80) * t125;
t42 = qJDD(3) * t119 + t121 * t140;
t262 = qJD(2) * t254 - t119 * t41 + t121 * t42;
t114 = t125 * qJDD(2);
t199 = qJD(3) * t125;
t146 = t123 * t42 + t82 * t199;
t200 = qJD(3) * t123;
t180 = t124 * t200;
t212 = t125 * t126;
t71 = t119 * t124 + t121 * t212;
t35 = qJD(2) * t71 - t121 * t180;
t213 = t124 * t125;
t67 = -t119 * t126 + t121 * t213;
t261 = (-t125 * t35 + (qJD(3) * t67 - t126 * t82) * t123) * qJD(2) - t67 * t114 - t146 * t124;
t228 = t121 * t41;
t229 = t119 * t42;
t260 = qJD(3) * t254 + t123 * (t228 + t229);
t158 = pkin(3) * t123 - qJ(4) * t125;
t62 = qJD(3) * t158 - qJD(4) * t123;
t259 = -t71 * qJD(1) + t119 * t62;
t78 = t82 ^ 2;
t258 = -t80 ^ 2 - t78;
t215 = t123 * t126;
t64 = t120 * t215 + t122 * t125;
t68 = -t120 * t125 + t122 * t215;
t162 = g(1) * t68 + g(2) * t64;
t216 = t123 * t124;
t136 = g(3) * t216 + t162;
t195 = qJD(1) * qJD(2);
t75 = qJDD(2) * pkin(6) + qJDD(1) * t124 + t126 * t195;
t59 = t123 * t75;
t26 = -qJDD(3) * pkin(3) + t96 * t199 + qJDD(4) + t59;
t131 = -t136 + t26;
t186 = t119 * t212;
t256 = qJD(1) * t186 + (-t206 + t62) * t121;
t252 = t160 * t124;
t177 = t123 * t194;
t251 = t177 - t114;
t250 = pkin(4) * t41 - qJ(5) * t42;
t173 = t121 * t114;
t128 = qJD(2) ^ 2;
t218 = t118 * t128;
t246 = -t119 * t218 + (t80 + t197) * t204 - t173;
t149 = pkin(3) * t125 + qJ(4) * t123 + pkin(2);
t110 = t124 * t195;
t166 = -qJDD(1) * t126 + t110;
t12 = qJD(2) * t62 - qJDD(2) * t149 + t166;
t86 = t123 * t96;
t23 = qJDD(3) * qJ(4) + t125 * t75 + (qJD(4) - t86) * qJD(3);
t5 = t119 * t12 + t121 * t23;
t244 = pkin(2) * t124;
t243 = pkin(4) * t121;
t242 = pkin(6) * t125;
t241 = pkin(6) * t126;
t238 = g(3) * t124;
t237 = g(3) * t126;
t198 = qJD(5) * t125;
t236 = -t198 + (-pkin(6) * t121 + qJ(5)) * t200 + t259;
t184 = pkin(6) * t119 + pkin(4);
t235 = -t184 * t200 - t256;
t189 = pkin(6) * t200;
t234 = t119 * t189 + t256;
t233 = -t121 * t189 + t259;
t52 = -qJD(2) * t149 - t205;
t87 = t125 * t96;
t73 = qJD(3) * qJ(4) + t87;
t15 = t119 * t52 + t121 * t73;
t232 = t160 * t216;
t85 = t158 * qJD(2);
t226 = t121 * t85;
t225 = t121 * t149;
t48 = -t119 * t149 + t121 * t242;
t221 = qJ(4) * t121;
t220 = qJD(5) * t82;
t219 = qJDD(2) * pkin(2);
t217 = t121 * t123;
t214 = t124 * t121;
t211 = qJDD(1) - g(3);
t210 = t126 * pkin(2) + t124 * pkin(6);
t127 = qJD(3) ^ 2;
t207 = t127 + t128;
t203 = qJD(2) * t125;
t202 = qJD(2) * t126;
t196 = t121 * qJD(4);
t193 = qJDD(2) * t126;
t192 = qJDD(3) * t123;
t174 = t119 * t114;
t181 = t119 * t203;
t185 = t123 * t214;
t190 = g(3) * t185 + qJ(4) * t174 + qJD(4) * t181;
t188 = t80 * t205;
t187 = t82 * t205;
t183 = qJ(4) * t200;
t182 = t123 * t205;
t179 = t124 * t199;
t4 = -t119 * t23 + t121 * t12;
t172 = t208 * t75;
t167 = t208 * qJDD(2);
t165 = pkin(3) * t212 + qJ(4) * t215 + t210;
t164 = t123 * t176;
t163 = -pkin(3) * t216 + qJ(4) * t213;
t65 = t120 * t212 - t122 * t123;
t69 = t120 * t123 + t122 * t212;
t161 = g(1) * t69 + g(2) * t65;
t63 = -qJD(3) * pkin(3) + qJD(4) + t86;
t157 = pkin(4) * t119 - qJ(5) * t121;
t156 = -qJD(2) * t97 + t238;
t14 = -t119 * t73 + t121 * t52;
t152 = qJ(4) * t42 + qJD(4) * t82;
t150 = qJD(2) * (-t82 + t201);
t148 = pkin(6) + t157;
t147 = pkin(4) * t114 + qJDD(5) - t4;
t66 = t119 * t213 + t121 * t126;
t34 = -qJD(2) * t214 - t119 * t180 + t126 * t181;
t145 = t34 * t82 - t35 * t80 - t67 * t41 + t42 * t66;
t43 = t66 * t120;
t45 = t66 * t122;
t70 = t186 - t214;
t142 = -g(1) * t45 - g(2) * t43 + g(3) * t70;
t44 = t67 * t120;
t46 = t67 * t122;
t141 = g(1) * t46 + g(2) * t44 - g(3) * t71;
t139 = -t80 * t181 - t228;
t137 = (t123 * t41 + t80 * t199) * t119;
t74 = t166 - t219;
t135 = -pkin(6) * t127 + t219 - t237 - t74;
t134 = -g(3) * t213 - t80 * t196 - t41 * t221 - t161;
t133 = -pkin(6) * qJDD(3) + (t205 + t97 - t230) * qJD(3);
t132 = t125 * t150 - t112 + t175;
t130 = t66 * t114 + t41 * t216 + t202 * t224 + t80 * t179 + (t125 * t34 - t66 * t200) * qJD(2);
t129 = t149 * t252;
t104 = t122 * t241;
t101 = t120 * t241;
t100 = t123 * t128 * t125;
t88 = -qJ(5) * t119 - pkin(3) - t243;
t76 = qJDD(2) * t118 - 0.2e1 * t164;
t72 = t119 * t85;
t61 = t68 * pkin(3);
t60 = t64 * pkin(3);
t58 = t80 * t203;
t53 = t148 * t123;
t47 = -t119 * t242 - t225;
t39 = t184 * t125 + t225;
t38 = -qJ(5) * t125 + t48;
t31 = t157 * t203 + t87;
t30 = -t96 * t217 + t72;
t29 = t119 * t86 + t226;
t28 = -qJD(5) * t217 + t148 * t199;
t25 = -t226 + (-pkin(4) * qJD(2) - t119 * t96) * t123;
t24 = t72 + (qJ(5) * qJD(2) - t121 * t96) * t123;
t20 = t58 + t42;
t16 = t121 * t218 + t123 * t150 - t174;
t13 = pkin(4) * t80 - qJ(5) * t82 + t63;
t11 = -t121 * t82 * t203 + t229;
t10 = -qJ(5) * t203 + t15;
t9 = pkin(4) * t203 + qJD(5) - t14;
t7 = t146 * t121;
t6 = (-t121 * t191 - t42) * t125 + (t121 * t249 + t123 * t82) * qJD(3);
t3 = -t220 + t26 + t250;
t2 = -pkin(4) * t177 + t147;
t1 = t251 * qJ(5) - qJD(2) * t198 + t5;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t211, 0, 0, 0, 0, 0, 0, -t124 * t128 + t193, -qJDD(2) * t124 - t126 * t128, 0, -g(3) + (t124 ^ 2 + t126 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, (-0.2e1 * t177 + t114) * t126 + (-t125 * t207 - t192) * t124, (-qJDD(3) * t124 - 0.2e1 * t126 * t194) * t125 + (t124 * t207 - t193) * t123, t124 * t167 + t128 * t171, qJD(2) * t264 + t124 * t172 - t126 * t74 - g(3), 0, 0, 0, 0, 0, 0, t130, -t261, t145, t63 * t179 - t14 * t34 + t15 * t35 - t4 * t66 + t5 * t67 - g(3) + (t124 * t26 + t202 * t63) * t123, 0, 0, 0, 0, 0, 0, t130, t145, t261, t13 * t179 + t1 * t67 + t10 * t35 + t2 * t66 + t34 * t9 - g(3) + (t124 * t3 + t13 * t202) * t123; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t211 * t126 + t252, -t124 * t211 + t126 * t160, 0, 0, qJDD(2) * t117 + 0.2e1 * t164, -0.2e1 * qJD(3) * t249 + 0.2e1 * t114 * t123, t125 * t127 + t192, t76, qJDD(3) * t125 - t123 * t127, 0, t133 * t123 + ((t160 + t195) * t124 + t135) * t125, t133 * t125 + (-t135 - t110) * t123 - t232, -t238 + t172 + pkin(6) * t167 + (-t195 * t208 - t160) * t126, -t74 * pkin(2) - g(1) * (-t122 * t244 + t104) - g(2) * (-t120 * t244 + t101) - g(3) * t210 + pkin(6) * t172 - t264 * qJD(1), t7, -t260, t6, t137, -t265, t76, (-t188 + pkin(6) * t41 + t119 * t26 + (qJD(2) * t47 + t14) * qJD(3)) * t123 + (-qJDD(2) * t47 - t4 + (pkin(6) * t80 + t119 * t63) * qJD(3) - t234 * qJD(2)) * t125 + t141, (-t187 + pkin(6) * t42 + t121 * t26 + (-qJD(2) * t48 - t15) * qJD(3)) * t123 + (qJDD(2) * t48 + t5 + (pkin(6) * t82 + t121 * t63) * qJD(3) + t233 * qJD(2)) * t125 + t142, -t41 * t48 - t42 * t47 - t234 * t82 - t233 * t80 + (-t119 * t15 - t121 * t14) * t199 + (-t119 * t5 - t121 * t4 - t237) * t123 + t232, t5 * t48 + t4 * t47 - t63 * t182 - g(1) * t104 - g(2) * t101 - g(3) * t165 + t233 * t15 + t234 * t14 + (t26 * t123 + t199 * t63) * pkin(6) + t129, t7, t6, t260, t76, t265, t137, t28 * t80 + t41 * t53 + (-t188 + t119 * t3 + (-qJD(2) * t39 - t9) * qJD(3)) * t123 + (qJD(2) * t235 + qJDD(2) * t39 + t13 * t201 + t2) * t125 + t141, -t38 * t41 + t39 * t42 + t235 * t82 - t236 * t80 + (-t10 * t119 + t121 * t9) * t199 + (-t1 * t119 + t121 * t2 - t237) * t123 + t232, -t28 * t82 - t42 * t53 + (t187 - t121 * t3 + (qJD(2) * t38 + t10) * qJD(3)) * t123 + (-qJD(2) * t236 - qJDD(2) * t38 - t13 * t197 - t1) * t125 - t142, t1 * t38 + t3 * t53 + t2 * t39 - g(1) * (-pkin(4) * t46 - qJ(5) * t45 + t104) - g(2) * (-pkin(4) * t44 - qJ(5) * t43 + t101) - g(3) * (pkin(4) * t71 + t70 * qJ(5) + t165) + t235 * t9 + (t28 - t182) * t13 + t236 * t10 + t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, t209 * t128, t191, t100, t114, qJDD(3), t123 * t156 + t162 - t59, (t156 - t75) * t125 + t161, 0, 0, t11, t262, t16, t139, t246, t100, -t80 * t87 - pkin(3) * t41 + (t162 - t26) * t121 + (-t123 * t14 + t125 * t29 + (-t125 * t63 - t183) * t119) * qJD(2) + t190, -pkin(3) * t42 + (qJDD(2) * t221 - t82 * t96) * t125 + t131 * t119 + (t123 * t15 - t125 * t30 + (-t183 + (qJD(4) - t63) * t125) * t121) * qJD(2), t29 * t82 + t30 * t80 + (t14 * t203 + t5) * t121 + (t15 * t203 + t152 - t4) * t119 + t134, -t26 * pkin(3) - t15 * t30 - t14 * t29 - t63 * t87 + g(1) * t61 + g(2) * t60 - g(3) * t163 + (-t119 * t14 + t121 * t15) * qJD(4) + (-t119 * t4 + t121 * t5 - t161) * qJ(4), t11, t16, -t262, t100, -t246, t139, t41 * t88 + (-t119 * qJD(5) - t31) * t80 + (t162 - t3) * t121 + (t123 * t9 - t125 * t25 + (-t125 * t13 - t183) * t119) * qJD(2) + t190, t24 * t80 - t25 * t82 + (-t203 * t9 + t1) * t121 + (t10 * t203 + t152 + t2) * t119 + t134, -qJ(4) * t173 + t31 * t82 - t42 * t88 + (t136 - t3 + t220) * t119 + (-t10 * t123 + t125 * t24 + (t183 + (-qJD(4) + t13) * t125) * t121) * qJD(2), t1 * t221 + t3 * t88 - t13 * t31 - t9 * t25 - g(1) * (qJ(4) * t69 - t68 * t243 - t61) - g(2) * (qJ(4) * t65 - t64 * t243 - t60) - g(3) * (-pkin(4) * t185 + t163) + (-t24 + t196) * t10 + (t2 * qJ(4) + qJ(5) * t136 + t9 * qJD(4) - t13 * qJD(5)) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t20, t258, t14 * t82 + t15 * t80 + t131, 0, 0, 0, 0, 0, 0, t132, t258, -t20, t10 * t80 + (-qJD(5) - t9) * t82 + t131 + t250; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80 * t82 - t251, -t58 + t42, -t78 - t218, t13 * t82 - g(1) * (t119 * t69 - t122 * t214) - g(2) * (t119 * t65 - t120 * t214) - g(3) * t66 + (-pkin(4) * t200 + t10 * t125) * qJD(2) + t147;];
tau_reg = t8;
