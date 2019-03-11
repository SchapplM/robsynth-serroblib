% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
% 
% Output:
% tauc_reg [6x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RRPPPR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:12:06
% EndTime: 2019-03-09 08:12:16
% DurationCPUTime: 3.36s
% Computational Cost: add. (3544->332), mult. (9139->452), div. (0->0), fcn. (6674->8), ass. (0->177)
t161 = sin(pkin(9));
t163 = cos(pkin(9));
t165 = sin(qJ(2));
t167 = cos(qJ(2));
t136 = t161 * t167 + t163 * t165;
t125 = t136 * qJD(1);
t115 = qJD(6) + t125;
t203 = qJD(1) * qJD(2);
t197 = t165 * t203;
t220 = t163 * t167;
t198 = qJD(1) * t220;
t112 = qJD(2) * t198 - t161 * t197;
t160 = sin(pkin(10));
t162 = cos(pkin(10));
t164 = sin(qJ(6));
t166 = cos(qJ(6));
t135 = t166 * t160 + t164 * t162;
t226 = t115 * t135;
t245 = -t164 * t160 + t166 * t162;
t249 = t245 * t112 - t226 * t115;
t134 = t161 * t165 - t220;
t78 = t245 * t134;
t205 = qJD(6) * t166;
t206 = qJD(6) * t164;
t225 = t245 * t125 - t160 * t206 + t162 * t205;
t192 = -t135 * t112 - t225 * t115;
t204 = t165 * qJD(1);
t122 = t161 * t204 - t198;
t94 = t160 * qJD(2) - t162 * t122;
t96 = t162 * qJD(2) + t160 * t122;
t183 = t164 * t94 - t166 * t96;
t248 = t115 * t183;
t247 = -0.2e1 * t203;
t231 = -qJ(3) - pkin(7);
t144 = t231 * t167;
t140 = qJD(1) * t144;
t130 = t161 * t140;
t143 = t231 * t165;
t139 = qJD(1) * t143;
t89 = t163 * t139 + t130;
t246 = qJD(4) - t89;
t239 = t125 ^ 2;
t244 = -t122 ^ 2 - t239;
t243 = -qJD(6) + t115;
t153 = -t163 * pkin(2) - pkin(3);
t149 = -qJ(5) + t153;
t124 = t136 * qJD(2);
t111 = qJD(1) * t124;
t151 = t161 * pkin(2) + qJ(4);
t224 = t151 * t111;
t237 = t122 * pkin(4);
t133 = qJD(2) * pkin(2) + t139;
t221 = t163 * t140;
t83 = t161 * t133 - t221;
t80 = -qJD(2) * qJ(4) - t83;
t56 = qJD(5) - t80 - t237;
t242 = -t112 * t149 + (qJD(5) - t56) * t125 + t224;
t241 = -t160 * t112 - t162 * t239;
t23 = -t183 * qJD(6) - t245 * t111;
t238 = t111 * pkin(3);
t236 = t125 * pkin(4);
t235 = t162 * pkin(8);
t196 = qJD(2) * t231;
t116 = t167 * qJD(3) + t165 * t196;
t106 = t116 * qJD(1);
t117 = -t165 * qJD(3) + t167 * t196;
t107 = t117 * qJD(1);
t63 = t161 * t106 - t163 * t107;
t92 = -t163 * t143 - t161 * t144;
t234 = t63 * t92;
t233 = pkin(3) + qJ(5);
t232 = -pkin(8) + t149;
t150 = pkin(2) * t197;
t194 = -t112 * qJ(4) + t150;
t177 = -t125 * qJD(4) + t194;
t25 = t122 * qJD(5) + t233 * t111 + t177;
t37 = t112 * pkin(4) - qJD(2) * qJD(5) + t63;
t8 = t160 * t37 + t162 * t25;
t207 = qJD(2) * t165;
t127 = qJD(2) * t220 - t161 * t207;
t155 = pkin(2) * t207;
t176 = -t127 * qJ(4) - t136 * qJD(4) + t155;
t30 = t134 * qJD(5) + t233 * t124 + t176;
t73 = t161 * t116 - t163 * t117;
t51 = t127 * pkin(4) + t73;
t13 = t160 * t51 + t162 * t30;
t201 = -t167 * pkin(2) - pkin(1);
t188 = t201 * qJD(1);
t142 = qJD(3) + t188;
t171 = -t125 * qJ(4) + t142;
t44 = t233 * t122 + t171;
t82 = t163 * t133 + t130;
t187 = qJD(4) - t82;
t48 = -t233 * qJD(2) + t187 + t236;
t18 = t160 * t48 + t162 * t44;
t193 = pkin(2) * t204 + t122 * qJ(4);
t49 = t233 * t125 + t193;
t88 = t161 * t139 - t221;
t65 = t88 - t237;
t21 = t160 * t65 + t162 * t49;
t179 = -t136 * qJ(4) + t201;
t62 = t233 * t134 + t179;
t75 = t136 * pkin(4) + t92;
t27 = t160 * t75 + t162 * t62;
t53 = t164 * t96 + t166 * t94;
t229 = t122 * t53;
t228 = t53 * t115;
t227 = t183 * t122;
t222 = t162 * t111;
t169 = qJD(1) ^ 2;
t215 = t167 * t169;
t168 = qJD(2) ^ 2;
t214 = t168 * t165;
t213 = t168 * t167;
t200 = -pkin(5) * t162 - pkin(4);
t211 = -t200 * t125 + t246;
t210 = t236 + t246;
t64 = t163 * t106 + t161 * t107;
t208 = t165 ^ 2 - t167 ^ 2;
t156 = qJD(2) * qJD(4);
t202 = t156 + t64;
t36 = t162 * t37;
t4 = t112 * pkin(5) + t36 + (-pkin(8) * t111 - t25) * t160;
t5 = pkin(8) * t222 + t8;
t199 = -t164 * t5 + t166 * t4;
t17 = -t160 * t44 + t162 * t48;
t195 = pkin(1) * t247;
t7 = -t160 * t25 + t36;
t190 = t8 * t160 + t7 * t162;
t189 = t164 * t4 + t166 * t5;
t10 = t125 * pkin(5) - t96 * pkin(8) + t17;
t11 = -t94 * pkin(8) + t18;
t1 = t166 * t10 - t164 * t11;
t2 = t164 * t10 + t166 * t11;
t71 = t162 * t75;
t16 = t136 * pkin(5) + t71 + (-pkin(8) * t134 - t62) * t160;
t19 = t134 * t235 + t27;
t186 = t166 * t16 - t164 * t19;
t185 = t164 * t16 + t166 * t19;
t184 = -t160 * t18 - t162 * t17;
t181 = t162 * t112 - t160 * t239;
t74 = t163 * t116 + t161 * t117;
t93 = t161 * t143 - t163 * t144;
t38 = -t111 * pkin(4) + t202;
t67 = t122 * pkin(3) + t171;
t178 = t67 * t125 + t63;
t22 = t135 * t111 - t94 * t205 - t96 * t206;
t120 = t232 * t160;
t61 = t162 * t65;
t175 = qJD(5) * t162 + qJD(6) * t120 - t122 * pkin(5) + t61 + (-pkin(8) * t125 - t49) * t160;
t121 = t232 * t162;
t174 = qJD(5) * t160 - qJD(6) * t121 + t125 * t235 + t21;
t79 = t135 * t134;
t76 = -t134 * pkin(4) + t93;
t172 = t111 * t76 + t124 * t56 + t134 * t38;
t170 = -t93 * t111 + t92 * t112 - t74 * t122 + t73 * t125 + t63 * t136;
t141 = t160 * pkin(5) + t151;
t114 = qJD(2) * t122;
t81 = t134 * pkin(3) + t179;
t77 = -qJD(2) * pkin(3) + t187;
t72 = t125 * pkin(3) + t193;
t57 = t124 * pkin(3) + t176;
t52 = -t124 * pkin(4) + t74;
t50 = t200 * t134 + t93;
t47 = t162 * t51;
t41 = t177 + t238;
t34 = t200 * t124 + t74;
t33 = t94 * pkin(5) + t56;
t32 = qJD(6) * t79 - t245 * t124;
t31 = qJD(6) * t78 + t135 * t124;
t28 = t200 * t111 + t202;
t26 = -t160 * t62 + t71;
t20 = -t160 * t49 + t61;
t12 = -t160 * t30 + t47;
t9 = t124 * t235 + t13;
t6 = t127 * pkin(5) + t47 + (-pkin(8) * t124 - t30) * t160;
t3 = [0, 0, 0, 0.2e1 * t167 * t197, t208 * t247, t213, -t214, 0, -pkin(7) * t213 + t165 * t195, pkin(7) * t214 + t167 * t195, -t83 * t124 - t82 * t127 - t64 * t134 + t170, t234 + t64 * t93 - t82 * t73 + t83 * t74 + (t142 + t188) * t155, t80 * t124 + t77 * t127 - t134 * t202 + t170, t73 * qJD(2) - t81 * t111 - t57 * t122 - t67 * t124 - t41 * t134, t74 * qJD(2) - t81 * t112 - t57 * t125 - t67 * t127 - t41 * t136, t202 * t93 + t41 * t81 + t67 * t57 + t77 * t73 - t80 * t74 + t234, t26 * t112 + t12 * t125 + t17 * t127 + t7 * t136 - t162 * t172 + t52 * t94, -t27 * t112 - t13 * t125 - t18 * t127 - t8 * t136 + t160 * t172 + t52 * t96, -t12 * t96 - t13 * t94 + (t111 * t27 + t124 * t18 + t134 * t8) * t162 + (-t111 * t26 - t124 * t17 - t134 * t7) * t160, t17 * t12 + t18 * t13 + t7 * t26 + t8 * t27 + t38 * t76 + t56 * t52, -t183 * t31 + t22 * t79, t183 * t32 + t22 * t78 - t79 * t23 - t31 * t53, t79 * t112 + t31 * t115 - t127 * t183 + t22 * t136, t78 * t112 - t32 * t115 - t53 * t127 - t23 * t136, t112 * t136 + t115 * t127 (-t164 * t9 + t166 * t6) * t115 + t186 * t112 + t199 * t136 + t1 * t127 + t34 * t53 + t50 * t23 - t28 * t78 + t33 * t32 + (-t115 * t185 - t136 * t2) * qJD(6) -(t164 * t6 + t166 * t9) * t115 - t185 * t112 - t189 * t136 - t2 * t127 - t34 * t183 + t50 * t22 + t28 * t79 + t33 * t31 + (-t1 * t136 - t115 * t186) * qJD(6); 0, 0, 0, -t165 * t215, t208 * t169, 0, 0, 0, t169 * pkin(1) * t165, pkin(1) * t215 (t83 - t88) * t125 + (-t82 + t89) * t122 + (-t111 * t161 - t112 * t163) * pkin(2), t82 * t88 - t83 * t89 + (-t142 * t204 + t161 * t64 - t163 * t63) * pkin(2), -t224 + t153 * t112 + (-t80 - t88) * t125 + (t77 - t246) * t122, -t88 * qJD(2) + t72 * t122 + t178, -t89 * qJD(2) - t67 * t122 + t72 * t125 + 0.2e1 * t156 + t64, t151 * t202 + t63 * t153 - t246 * t80 - t67 * t72 - t77 * t88, t17 * t122 - t20 * t125 + t38 * t160 - t242 * t162 + t210 * t94, -t18 * t122 + t21 * t125 + t242 * t160 + t38 * t162 + t210 * t96, t20 * t96 + t21 * t94 + (qJD(5) * t96 - t125 * t18 - t7) * t162 + (qJD(5) * t94 + t125 * t17 - t8) * t160, qJD(5) * t184 + t149 * t190 + t38 * t151 - t17 * t20 - t18 * t21 + t210 * t56, t183 * t226 + t22 * t245, -t22 * t135 + t183 * t225 + t226 * t53 - t23 * t245, -t227 + t249, t192 - t229, t115 * t122 (-t164 * t120 + t166 * t121) * t112 + t141 * t23 + t28 * t135 + t1 * t122 + t211 * t53 + t225 * t33 + (t164 * t174 - t166 * t175) * t115 -(t166 * t120 + t164 * t121) * t112 + t141 * t22 + t28 * t245 - t2 * t122 - t211 * t183 - t226 * t33 + (t164 * t175 + t166 * t174) * t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t244, t83 * t122 + t82 * t125 + t150, t244, -0.2e1 * t125 * qJD(2), -t112 + t114, t238 - t80 * t122 + (-qJD(4) - t77) * t125 + t194, t122 * t94 + t241, t122 * t96 - t181 (t160 * t94 + t162 * t96) * t125 + (t160 ^ 2 + t162 ^ 2) * t111, t56 * t122 + t125 * t184 - t7 * t160 + t8 * t162, 0, 0, 0, 0, 0, t192 + t229, -t227 - t249; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112 + t114, -t125 * t122, -t239 - t168, t80 * qJD(2) + t178, -qJD(2) * t94 + t181, -qJD(2) * t96 + t241 (t160 * t96 - t162 * t94) * t125, -t56 * qJD(2) + (-t160 * t17 + t162 * t18) * t125 + t190, 0, 0, 0, 0, 0, -qJD(2) * t53 + t249, qJD(2) * t183 + t192; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96 * t125 - t222, t160 * t111 - t94 * t125, -t94 ^ 2 - t96 ^ 2, t17 * t96 + t18 * t94 + t38, 0, 0, 0, 0, 0, t23 - t248, t22 - t228; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t183 * t53, t183 ^ 2 - t53 ^ 2, t22 + t228, -t23 - t248, t112, t183 * t33 + t243 * t2 + t199, t243 * t1 + t33 * t53 - t189;];
tauc_reg  = t3;
