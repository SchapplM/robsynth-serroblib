% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
% 
% Output:
% tauc_reg [6x26]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 03:14
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRRPRP4_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 03:12:47
% EndTime: 2021-01-16 03:12:58
% DurationCPUTime: 2.66s
% Computational Cost: add. (2417->328), mult. (5757->445), div. (0->0), fcn. (3797->8), ass. (0->188)
t122 = sin(qJ(3));
t199 = qJD(2) * t122;
t108 = qJD(5) + t199;
t124 = cos(qJ(5));
t121 = sin(qJ(5));
t195 = qJD(3) * t121;
t125 = cos(qJ(3));
t197 = qJD(2) * t125;
t77 = t124 * t197 + t195;
t219 = t77 * t108;
t186 = qJD(2) * qJD(3);
t168 = t122 * t186;
t51 = t77 * qJD(5) - t121 * t168;
t250 = t51 - t219;
t249 = t51 + t219;
t234 = pkin(4) + pkin(8);
t120 = cos(pkin(6));
t201 = qJD(1) * t120;
t123 = sin(qJ(2));
t119 = sin(pkin(6));
t202 = qJD(1) * t119;
t179 = t123 * t202;
t88 = qJD(2) * pkin(8) + t179;
t55 = t122 * t88 - t125 * t201;
t248 = qJD(4) + t55;
t126 = cos(qJ(2));
t200 = qJD(2) * t119;
t170 = qJD(1) * t200;
t153 = t126 * t170;
t196 = qJD(3) * t120;
t247 = qJD(1) * t196 + t153;
t127 = -pkin(3) - pkin(9);
t204 = pkin(4) * t199 + t248;
t33 = t127 * qJD(3) + t204;
t173 = t126 * t202;
t166 = -qJ(4) * t122 - pkin(2);
t74 = t127 * t125 + t166;
t45 = t74 * qJD(2) - t173;
t14 = t121 * t33 + t124 * t45;
t167 = t125 * t186;
t192 = qJD(3) * t125;
t30 = t247 * t122 + t88 * t192;
t23 = pkin(4) * t167 + t30;
t151 = pkin(9) * t122 - qJ(4) * t125;
t191 = qJD(4) * t122;
t133 = t151 * qJD(3) - t191;
t217 = pkin(3) * t168 + t123 * t170;
t31 = t133 * qJD(2) + t217;
t164 = -t121 * t31 + t124 * t23;
t132 = -t14 * qJD(5) + t164;
t130 = qJ(6) * t51 + t132;
t159 = pkin(5) * t167;
t172 = t121 * t197;
t193 = qJD(3) * t124;
t79 = -t172 + t193;
t1 = -qJD(6) * t79 + t130 + t159;
t189 = qJD(5) * t124;
t190 = qJD(5) * t121;
t162 = -t121 * t23 - t124 * t31 - t33 * t189 + t45 * t190;
t52 = qJD(3) * t189 - qJD(5) * t172 - t124 * t168;
t143 = qJ(6) * t52 + t162;
t2 = -qJD(6) * t77 - t143;
t9 = -qJ(6) * t77 + t14;
t230 = t108 * t9;
t13 = -t121 * t45 + t124 * t33;
t8 = -qJ(6) * t79 + t13;
t7 = pkin(5) * t108 + t8;
t246 = -(t108 * t7 - t2) * t121 + (t1 + t230) * t124;
t208 = t124 * t126;
t85 = t234 * t192;
t244 = -(-t121 * t123 + t122 * t208) * t202 + t124 * t85;
t116 = qJD(3) * qJ(4);
t56 = t122 * t201 + t125 * t88;
t49 = -t116 - t56;
t211 = t121 * t126;
t194 = qJD(3) * t122;
t111 = pkin(3) * t194;
t58 = t111 + t133;
t94 = t234 * t122;
t243 = -t121 * t85 - t124 * t58 - t94 * t189 + t74 * t190 + (t122 * t211 + t123 * t124) * t202;
t48 = -qJD(3) * pkin(3) + t248;
t224 = t108 * t79;
t242 = -t52 + t224;
t241 = t52 + t224;
t128 = qJD(3) ^ 2;
t137 = -qJ(4) * t192 - t191;
t38 = t137 * qJD(2) + t217;
t64 = t111 + t137;
t238 = (-t64 + t179) * qJD(2) - pkin(8) * t128 - t38;
t177 = t126 * t200;
t129 = qJD(2) ^ 2;
t213 = t119 * t129;
t182 = t123 * t213;
t214 = t119 * t123;
t183 = t122 * t214;
t39 = -qJD(3) * t183 + (t177 + t196) * t125;
t237 = (t125 * t177 + t39) * qJD(3) - t122 * t182;
t154 = t122 * t177;
t69 = t120 * t122 + t125 * t214;
t40 = t69 * qJD(3) + t154;
t236 = (t40 + t154) * qJD(3) + t125 * t182;
t235 = t79 ^ 2;
t233 = t7 - t8;
t161 = qJ(6) * t125 - t74;
t232 = pkin(5) * t192 + t161 * t189 + (-qJ(6) * t194 - qJD(5) * t94 + qJD(6) * t125 - t58) * t121 + t244;
t188 = qJD(5) * t125;
t176 = t121 * t188;
t187 = t124 * qJD(6);
t231 = t125 * t187 - (t122 * t193 + t176) * qJ(6) + t243;
t44 = pkin(4) * t197 + t56;
t112 = pkin(3) * t199;
t62 = t151 * qJD(2) + t112;
t163 = -t121 * t62 + t124 * t44;
t205 = qJ(6) - t127;
t212 = t121 * t122;
t229 = (pkin(5) * t125 - qJ(6) * t212) * qJD(2) + t163 - t205 * t190 + t187;
t227 = t121 * t44 + t124 * t62;
t91 = t205 * t124;
t228 = qJ(6) * t124 * t199 + qJD(5) * t91 + t121 * qJD(6) + t227;
t226 = t121 * t94 + t124 * t74;
t225 = qJD(2) * pkin(2);
t223 = t124 * t51;
t222 = t125 * t79;
t115 = qJD(3) * qJD(4);
t185 = t247 * t125 - t88 * t194;
t25 = -t115 - t185;
t19 = -pkin(4) * t168 - t25;
t221 = t19 * t121;
t220 = t19 * t124;
t180 = -pkin(5) * t124 - pkin(4);
t218 = pkin(5) * t189 - t180 * t199 + t248;
t92 = -pkin(3) * t125 + t166;
t216 = qJD(2) * t92;
t215 = t108 * t127;
t210 = t122 * t124;
t209 = t124 * t125;
t207 = t128 * t122;
t206 = t128 * t125;
t95 = t234 * t125;
t117 = t122 ^ 2;
t118 = t125 ^ 2;
t203 = t117 - t118;
t198 = qJD(2) * t123;
t184 = t108 * t210;
t181 = t122 * t129 * t125;
t178 = t119 * t198;
t175 = t108 * t189;
t174 = t124 * t188;
t165 = -pkin(5) * t77 - qJD(6);
t35 = t116 + t44;
t158 = t77 * t173;
t157 = t79 * t173;
t150 = qJD(3) * t55 + t185;
t149 = qJD(3) * t56 - t30;
t148 = -qJD(2) * t118 + t108 * t122;
t146 = t108 * t121;
t12 = pkin(5) * t52 + t19;
t22 = -t165 + t35;
t142 = -t12 * t121 - t22 * t189;
t141 = t12 * t124 - t22 * t190;
t68 = -t120 * t125 + t183;
t140 = t119 * t208 - t121 * t68;
t41 = t119 * t211 + t124 * t68;
t139 = t122 * t35 + t127 * t192;
t89 = -t173 - t225;
t136 = qJD(3) * (t173 + t89 - t225);
t57 = -t173 + t216;
t135 = qJD(3) * (-t173 - t57 - t216);
t131 = t122 * t30 - t125 * t25 + (t122 * t49 + t125 * t48) * qJD(3);
t109 = pkin(5) * t121 + qJ(4);
t100 = t124 * t167;
t90 = t205 * t121;
t84 = t234 * t194;
t83 = -qJ(4) * t197 + t112;
t82 = t124 * t94;
t76 = t77 ^ 2;
t67 = pkin(5) * t209 + t95;
t47 = t57 * t199;
t46 = -pkin(5) * t176 + (-pkin(8) + t180) * t194;
t29 = -qJ(6) * t209 + t226;
t24 = t122 * pkin(5) + t161 * t121 + t82;
t18 = -t175 - qJD(3) * t79 + (-t121 * t192 - t184) * qJD(2);
t17 = -qJD(3) * t77 - t108 * t146 + t100;
t11 = t41 * qJD(5) + t40 * t121 + t124 * t178;
t10 = t140 * qJD(5) - t121 * t178 + t40 * t124;
t4 = -t108 * t11 + t140 * t167 + t39 * t79 - t51 * t69;
t3 = t10 * t108 + t41 * t167 + t39 * t77 + t52 * t69;
t5 = [0, 0, -t182, -t126 * t213, 0, 0, 0, 0, 0, -t236, -t237, (t122 * t40 + t125 * t39 + (-t122 * t69 + t125 * t68) * qJD(3)) * qJD(2), t236, t237, -t25 * t69 + t30 * t68 - t39 * t49 + t40 * t48 + (-t126 * t38 + t57 * t198) * t119, 0, 0, 0, 0, 0, t3, t4, t3, t4, -t10 * t79 - t11 * t77 + t140 * t52 + t41 * t51, t1 * t41 + t10 * t7 + t11 * t9 + t12 * t69 - t140 * t2 + t22 * t39; 0, 0, 0, 0, 0.2e1 * t122 * t167, -0.2e1 * t203 * t186, t206, -t207, 0, -pkin(8) * t206 + t122 * t136, pkin(8) * t207 + t125 * t136, (-t117 - t118) * t153 + t131, t122 * t135 - t125 * t238, t122 * t238 + t125 * t135, t38 * t92 + t57 * t64 + (-t123 * t57 + (-t122 * t48 + t125 * t49) * t126) * t202 + t131 * pkin(8), -t79 * t174 + (t125 * t51 + t79 * t194) * t121, (-t121 * t77 + t124 * t79) * t194 + (t121 * t52 + t223 + (t121 * t79 + t124 * t77) * qJD(5)) * t125, -t108 * t174 - t51 * t122 + (t148 * t121 + t222) * qJD(3), t108 * t176 - t52 * t122 + (t148 * t124 - t125 * t77) * qJD(3), (t108 + t199) * t192, t95 * t52 - t84 * t77 + (-t35 * t193 + t164) * t122 + (-t121 * t58 + t244) * t108 + (-t226 * t108 - t14 * t122) * qJD(5) + (-t158 - t35 * t190 + t220 + ((-t121 * t74 + t82) * qJD(2) + t13) * qJD(3)) * t125, -t95 * t51 - t84 * t79 + (t35 * t195 + t162) * t122 + t243 * t108 + (-t157 - t35 * t189 - t221 + (-t226 * qJD(2) - t14) * qJD(3)) * t125, t46 * t77 + t52 * t67 + (-t22 * t193 + t1) * t122 + t232 * t108 + (-t158 + (qJD(2) * t24 + t7) * qJD(3) + t141) * t125, t46 * t79 - t51 * t67 + (t22 * t195 - t2) * t122 + t231 * t108 + (-t157 + (-qJD(2) * t29 - t9) * qJD(3) + t142) * t125, t24 * t51 - t29 * t52 - t232 * t79 + t231 * t77 + (-t121 * t7 + t124 * t9) * t194 + (t1 * t121 - t124 * t2 + (t121 * t9 + t124 * t7) * qJD(5)) * t125, t1 * t24 + t12 * t67 + t2 * t29 - t231 * t9 + t232 * t7 + (-t125 * t173 + t46) * t22; 0, 0, 0, 0, -t181, t203 * t129, 0, 0, 0, -t89 * t199 + t149, -t89 * t197 - t150, 0, -t83 * t197 - t149 + t47, 0.2e1 * t115 + (t122 * t83 + t125 * t57) * qJD(2) + t150, -pkin(3) * t30 - qJ(4) * t25 - t248 * t49 - t48 * t56 - t57 * t83, -t79 * t146 - t223, t249 * t121 - t241 * t124, -t108 * t190 + t100 + (-t108 * t212 - t222) * qJD(2), -t175 + (-t184 + (t77 - t195) * t125) * qJD(2), -t108 * t197, qJ(4) * t52 + t221 - t163 * t108 + t204 * t77 + (-t121 * t215 + t124 * t35) * qJD(5) + (t124 * t139 - t13 * t125) * qJD(2), -qJ(4) * t51 + t220 + t227 * t108 + t204 * t79 + (-t121 * t35 - t124 * t215) * qJD(5) + (-t121 * t139 + t14 * t125) * qJD(2), t109 * t52 + t218 * t77 - t229 * t108 + (t22 * t210 + (-qJD(3) * t91 - t7) * t125) * qJD(2) - t142, -t109 * t51 + t218 * t79 + t228 * t108 + (-t22 * t212 + (qJD(3) * t90 + t9) * t125) * qJD(2) + t141, t228 * t77 + t229 * t79 - t51 * t91 + t52 * t90 - t246, -t1 * t91 + t109 * t12 - t2 * t90 + t218 * t22 - t228 * t9 - t229 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t181, -t117 * t129 - t128, qJD(3) * t49 + t30 + t47, 0, 0, 0, 0, 0, t17, t18, t17, t18, t242 * t121 + t250 * t124, -qJD(3) * t22 + t246; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79 * t77, -t76 + t235, -t250, t242, t167, t108 * t14 - t35 * t79 + t132, t108 * t13 + t35 * t77 + t162, 0.2e1 * t159 + t230 + (t165 - t22) * t79 + t130, -pkin(5) * t235 + t108 * t8 + (qJD(6) + t22) * t77 + t143, t51 * pkin(5) - t233 * t77, t233 * t9 + (-t22 * t79 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t241, -t249, -t76 - t235, t7 * t79 + t77 * t9 + t12;];
tauc_reg = t5;
