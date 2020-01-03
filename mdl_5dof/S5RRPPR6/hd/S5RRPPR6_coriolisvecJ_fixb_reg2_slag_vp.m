% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR6_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR6_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:33:05
% EndTime: 2019-12-31 19:33:18
% DurationCPUTime: 3.74s
% Computational Cost: add. (4967->330), mult. (12948->464), div. (0->0), fcn. (9562->8), ass. (0->179)
t157 = sin(qJ(5));
t159 = cos(qJ(5));
t156 = cos(pkin(9));
t154 = sin(pkin(9));
t155 = sin(pkin(8));
t158 = sin(qJ(2));
t160 = cos(qJ(2));
t214 = cos(pkin(8));
t132 = t155 * t160 + t214 * t158;
t195 = qJD(1) * t132;
t203 = t154 * t195;
t247 = t156 * qJD(2) - t203;
t170 = t159 * t247;
t96 = t154 * qJD(2) + t156 * t195;
t50 = -t157 * t96 + t170;
t251 = t50 ^ 2;
t182 = t214 * t160;
t142 = qJD(1) * t182;
t193 = qJD(1) * t158;
t116 = t155 * t193 - t142;
t113 = qJD(5) + t116;
t250 = t50 * t113;
t49 = t157 * t247 + t159 * t96;
t249 = t49 ^ 2;
t191 = qJD(1) * qJD(2);
t185 = t158 * t191;
t141 = t155 * t185;
t111 = qJD(2) * t142 - t141;
t202 = t156 * t111;
t248 = t116 * t247 + t202;
t200 = t159 * t156;
t246 = -t157 * t154 + t200;
t238 = t246 * qJD(5);
t216 = -t116 * t246 - t238;
t133 = t159 * t154 + t157 * t156;
t123 = t133 * qJD(5);
t215 = t133 * t116 + t123;
t204 = t154 * t111;
t243 = t116 * t96 + t204;
t242 = -0.2e1 * t191;
t118 = t132 * qJD(2);
t110 = qJD(1) * t118;
t169 = -t155 * t158 + t182;
t86 = t110 * t169;
t241 = t116 * t118 - t86;
t240 = t111 * t133;
t239 = t195 * t247;
t237 = -qJD(5) + t113;
t22 = t157 * (qJD(5) * t96 + t204) - qJD(5) * t170 - t111 * t200;
t236 = -t215 * t49 - t22 * t246;
t235 = t133 * t110 - t216 * t113;
t115 = t116 ^ 2;
t234 = -t154 * t110 - t156 * t115;
t233 = t195 ^ 2;
t229 = pkin(7) * t156;
t189 = pkin(2) * t193;
t69 = pkin(3) * t195 + t116 * qJ(4) + t189;
t225 = -qJ(3) - pkin(6);
t140 = t225 * t160;
t136 = qJD(1) * t140;
t124 = t155 * t136;
t139 = t225 * t158;
t135 = qJD(1) * t139;
t89 = t214 * t135 + t124;
t33 = -t154 * t89 + t156 * t69;
t20 = pkin(4) * t195 + t116 * t229 + t33;
t210 = t116 * t154;
t34 = t154 * t69 + t156 * t89;
t27 = pkin(7) * t210 + t34;
t145 = t155 * pkin(2) + qJ(4);
t226 = pkin(7) + t145;
t127 = t226 * t154;
t128 = t226 * t156;
t78 = -t159 * t127 - t157 * t128;
t232 = qJD(4) * t246 + qJD(5) * t78 - t157 * t20 - t159 * t27;
t79 = -t157 * t127 + t159 * t128;
t231 = -qJD(4) * t133 - qJD(5) * t79 + t157 * t27 - t159 * t20;
t230 = pkin(2) * t158;
t228 = t49 * t50;
t184 = qJD(2) * t225;
t114 = t160 * qJD(3) + t158 * t184;
t105 = t114 * qJD(1);
t167 = -t158 * qJD(3) + t160 * t184;
t106 = t167 * qJD(1);
t57 = t155 * t105 - t214 * t106;
t93 = -t214 * t139 - t155 * t140;
t227 = t57 * t93;
t144 = pkin(2) * t185;
t41 = t110 * pkin(3) - t111 * qJ(4) - qJD(4) * t195 + t144;
t58 = t214 * t105 + t155 * t106;
t55 = qJD(2) * qJD(4) + t58;
t17 = t154 * t41 + t156 * t55;
t121 = t169 * qJD(2);
t223 = qJD(2) * pkin(2);
t190 = t158 * t223;
t52 = t118 * pkin(3) - t121 * qJ(4) - t132 * qJD(4) + t190;
t71 = t214 * t114 + t155 * t167;
t25 = t154 * t52 + t156 * t71;
t149 = -t160 * pkin(2) - pkin(1);
t194 = qJD(1) * t149;
t138 = qJD(3) + t194;
t61 = t116 * pkin(3) - qJ(4) * t195 + t138;
t129 = t135 + t223;
t183 = t214 * t136;
t84 = t155 * t129 - t183;
t77 = qJD(2) * qJ(4) + t84;
t29 = t154 * t61 + t156 * t77;
t81 = -pkin(3) * t169 - t132 * qJ(4) + t149;
t94 = t155 * t139 - t214 * t140;
t38 = t154 * t81 + t156 * t94;
t222 = t195 * t50;
t220 = t49 * t195;
t219 = t96 * t195;
t218 = t96 * t154;
t217 = t156 * t110 - t154 * t115;
t213 = t111 * t132;
t211 = t116 * t195;
t209 = t121 * t154;
t208 = t121 * t156;
t207 = t132 * t154;
t162 = qJD(1) ^ 2;
t199 = t160 * t162;
t161 = qJD(2) ^ 2;
t198 = t161 * t158;
t197 = t161 * t160;
t196 = t158 ^ 2 - t160 ^ 2;
t188 = t158 * t199;
t16 = -t154 * t55 + t156 * t41;
t24 = -t154 * t71 + t156 * t52;
t28 = -t154 * t77 + t156 * t61;
t37 = -t154 * t94 + t156 * t81;
t181 = pkin(1) * t242;
t70 = t155 * t114 - t214 * t167;
t88 = t155 * t135 - t183;
t180 = 0.2e1 * t195;
t23 = qJD(5) * t49 + t240;
t178 = -t133 * t23 - t216 * t50;
t177 = t110 * t246 - t215 * t113;
t176 = t160 * t185;
t148 = -t214 * pkin(2) - pkin(3);
t36 = pkin(4) * t204 + t57;
t175 = t93 * t111 + t57 * t132;
t12 = t110 * pkin(4) - pkin(7) * t202 + t16;
t13 = -pkin(7) * t204 + t17;
t174 = t157 * t12 + t159 * t13;
t15 = t116 * pkin(4) - t96 * pkin(7) + t28;
t21 = pkin(7) * t247 + t29;
t5 = t159 * t15 - t157 * t21;
t6 = t157 * t15 + t159 * t21;
t173 = -t154 * t28 + t156 * t29;
t26 = -pkin(4) * t169 - t132 * t229 + t37;
t30 = -pkin(7) * t207 + t38;
t9 = -t157 * t30 + t159 * t26;
t10 = t157 * t26 + t159 * t30;
t83 = t214 * t129 + t124;
t171 = t156 * t247;
t72 = -qJD(2) * pkin(3) + qJD(4) - t83;
t168 = t72 * t121 + t175;
t166 = t132 * t110 - t111 * t169 + t121 * t116;
t165 = t171 - t218;
t163 = -t145 * t110 + t148 * t111 + (-qJD(4) + t72) * t116;
t2 = -qJD(5) * t6 + t159 * t12 - t157 * t13;
t151 = t156 ^ 2;
t150 = t154 ^ 2;
t137 = -t156 * pkin(4) + t148;
t74 = t246 * t132;
t73 = t133 * t132;
t67 = pkin(4) * t207 + t93;
t56 = -pkin(4) * t210 + t88;
t43 = pkin(4) * t209 + t70;
t42 = -pkin(4) * t247 + t72;
t32 = t121 * t133 + t238 * t132;
t31 = -t121 * t246 + t123 * t132;
t18 = -pkin(7) * t209 + t25;
t14 = t118 * pkin(4) - pkin(7) * t208 + t24;
t4 = -qJD(5) * t10 + t159 * t14 - t157 * t18;
t3 = qJD(5) * t9 + t157 * t14 + t159 * t18;
t1 = qJD(5) * t5 + t174;
t7 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t176, t196 * t242, t197, -0.2e1 * t176, -t198, 0, -pkin(6) * t197 + t158 * t181, pkin(6) * t198 + t160 * t181, 0, 0, t121 * t195 + t213, -t118 * t195 - t166, t121 * qJD(2), t241, -t118 * qJD(2), 0, t149 * t110 + t138 * t118 + (-t70 + (-qJD(1) * t169 + t116) * t230) * qJD(2), t149 * t111 + t138 * t121 + (t180 * t230 - t71) * qJD(2), -t94 * t110 - t71 * t116 - t84 * t118 - t83 * t121 + t169 * t58 + t195 * t70 + t175, t227 + t58 * t94 - t83 * t70 + t84 * t71 + (t138 + t194) * t190, t151 * t213 + t96 * t208, t121 * t165 - 0.2e1 * t202 * t207, t96 * t118 + t156 * t166, t150 * t213 - t209 * t247, t118 * t247 - t154 * t166, t241, t37 * t110 + t24 * t116 + t28 * t118 + t154 * t168 - t16 * t169 - t247 * t70, -t38 * t110 - t25 * t116 - t29 * t118 + t156 * t168 + t169 * t17 + t70 * t96, -t25 * t203 - t24 * t96 + (t25 * qJD(2) - t37 * t111 - t28 * t121 - t16 * t132) * t156 + (-t38 * t111 - t29 * t121 - t17 * t132) * t154, t16 * t37 + t17 * t38 + t28 * t24 + t29 * t25 + t72 * t70 + t227, -t22 * t74 - t49 * t31, t22 * t73 - t74 * t23 - t31 * t50 - t49 * t32, t74 * t110 - t31 * t113 + t49 * t118 + t169 * t22, t23 * t73 - t32 * t50, -t73 * t110 - t32 * t113 + t118 * t50 + t169 * t23, t113 * t118 - t86, t9 * t110 + t4 * t113 + t5 * t118 - t169 * t2 + t67 * t23 + t42 * t32 + t36 * t73 - t43 * t50, t1 * t169 - t10 * t110 - t3 * t113 - t6 * t118 - t67 * t22 - t42 * t31 + t36 * t74 + t43 * t49, -t1 * t73 - t10 * t23 - t2 * t74 + t9 * t22 + t3 * t50 + t5 * t31 - t6 * t32 - t4 * t49, t1 * t10 + t2 * t9 + t6 * t3 + t36 * t67 + t5 * t4 + t42 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t188, t196 * t162, 0, t188, 0, 0, t162 * pkin(1) * t158, pkin(1) * t199, 0, 0, t211, -t115 + t233, -t141 + (t142 + t116) * qJD(2), -t211, 0, 0, t88 * qJD(2) - t116 * t189 - t138 * t195 - t57, t89 * qJD(2) + t138 * t116 - t189 * t195 - t58, (t84 - t88) * t195 + (-t83 + t89) * t116 + (-t110 * t155 - t214 * t111) * pkin(2), t83 * t88 - t84 * t89 + (-t138 * t193 + t155 * t58 - t214 * t57) * pkin(2), t243 * t156, t165 * t116 + (-t150 + t151) * t111, -t219 - t234, -t248 * t154, t217 - t239, -t211, -t33 * t116 + t154 * t163 - t57 * t156 - t195 * t28 + t247 * t88, t34 * t116 + t57 * t154 + t156 * t163 + t195 * t29 - t88 * t96, t34 * t203 + t33 * t96 + (-qJD(4) * t203 - t28 * t116 + t17 + (t156 * qJD(4) - t34) * qJD(2)) * t156 + (qJD(4) * t96 - t29 * t116 - t16) * t154, t57 * t148 - t28 * t33 - t29 * t34 - t72 * t88 + (-t16 * t154 + t17 * t156) * t145 + t173 * qJD(4), -t22 * t133 - t216 * t49, t178 + t236, -t220 + t235, -t215 * t50 - t23 * t246, t177 - t222, -t113 * t195, t78 * t110 + t231 * t113 + t137 * t23 - t195 * t5 + t215 * t42 - t246 * t36 + t50 * t56, -t79 * t110 - t232 * t113 + t36 * t133 - t137 * t22 + t195 * t6 - t216 * t42 - t56 * t49, t1 * t246 - t2 * t133 - t215 * t6 + t216 * t5 + t78 * t22 - t79 * t23 - t231 * t49 + t232 * t50, t1 * t79 + t36 * t137 + t2 * t78 + t231 * t5 + t232 * t6 - t42 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t180 * qJD(2), -t141 + (t142 - t116) * qJD(2), -t115 - t233, t84 * t116 + t195 * t83 + t144, 0, 0, 0, 0, 0, 0, t217 + t239, -t219 + t234, (t171 + t218) * t116 + (-t150 - t151) * t111, t116 * t173 + t17 * t154 + t16 * t156 - t195 * t72, 0, 0, 0, 0, 0, 0, t177 + t222, -t220 - t235, t178 - t236, t1 * t133 - t195 * t42 + t2 * t246 - t215 * t5 - t216 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t243, t248, -t247 ^ 2 - t96 ^ 2, -t247 * t29 + t28 * t96 + t57, 0, 0, 0, 0, 0, 0, t49 * t113 + t23, -t22 + t250, -t249 - t251, t49 * t5 - t6 * t50 + t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t228, t249 - t251, -t22 - t250, t228, t237 * t49 - t240, t110, t6 * t113 - t42 * t49 + t2, t237 * t5 - t42 * t50 - t174, 0, 0;];
tauc_reg = t7;
