% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% tauc_reg [6x23]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 01:25
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPRRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 01:24:18
% EndTime: 2021-01-16 01:24:27
% DurationCPUTime: 2.15s
% Computational Cost: add. (2541->297), mult. (6369->411), div. (0->0), fcn. (4769->10), ass. (0->166)
t106 = sin(pkin(11));
t115 = cos(qJ(2));
t107 = sin(pkin(6));
t165 = qJD(1) * t107;
t150 = t115 * t165;
t108 = cos(pkin(11));
t112 = sin(qJ(2));
t151 = t112 * t165;
t92 = t108 * t151;
t60 = t106 * t150 + t92;
t111 = sin(qJ(4));
t114 = cos(qJ(4));
t133 = pkin(4) * t111 - pkin(9) * t114;
t88 = t133 * qJD(4);
t213 = -t60 + t88;
t113 = cos(qJ(5));
t161 = qJD(4) * t111;
t110 = sin(qJ(5));
t171 = t110 * t114;
t99 = pkin(2) * t106 + pkin(8);
t182 = t110 * t99;
t91 = t106 * t151;
t62 = t108 * t150 - t91;
t212 = t213 * t113 + t161 * t182 + t62 * t171;
t158 = qJD(5) * t113;
t169 = t113 * t114;
t126 = -pkin(4) * t114 - pkin(9) * t111 - pkin(3);
t202 = pkin(2) * t108;
t79 = t126 - t202;
t211 = t213 * t110 + t79 * t158 - t62 * t169;
t90 = qJD(2) * pkin(2) + t150;
t55 = t106 * t90 + t92;
t51 = qJD(2) * pkin(8) + t55;
t109 = cos(pkin(6));
t97 = qJD(1) * t109 + qJD(3);
t210 = -t111 * t51 + t114 * t97;
t157 = t113 * qJD(4);
t164 = qJD(2) * t111;
t82 = t110 * t164 - t157;
t163 = qJD(2) * t114;
t98 = -qJD(5) + t163;
t198 = t82 * t98;
t159 = qJD(5) * t110;
t147 = t111 * t159;
t148 = t114 * t157;
t154 = qJD(4) * qJD(5);
t58 = -t113 * t154 + (t147 - t148) * qJD(2);
t209 = -t58 + t198;
t162 = qJD(4) * t110;
t84 = t113 * t164 + t162;
t197 = t84 * t98;
t146 = t111 * t158;
t160 = qJD(4) * t114;
t120 = t110 * t160 + t146;
t59 = qJD(2) * t120 + t110 * t154;
t208 = t59 - t197;
t65 = (t106 * t112 - t108 * t115) * t107;
t104 = t111 ^ 2;
t128 = qJD(2) * t104 - t114 * t98;
t207 = -t128 * t157 - t98 * t147;
t206 = t84 ^ 2;
t36 = t111 * t97 + t114 * t51;
t33 = qJD(4) * pkin(9) + t36;
t54 = t108 * t90 - t91;
t41 = qJD(2) * t126 - t54;
t11 = -t110 * t33 + t113 * t41;
t8 = -qJ(6) * t84 + t11;
t7 = -pkin(5) * t98 + t8;
t205 = t7 - t8;
t204 = pkin(5) * t82;
t183 = t110 * t41;
t12 = t113 * t33 + t183;
t9 = -qJ(6) * t82 + t12;
t203 = t9 * t98;
t201 = pkin(5) * t110;
t200 = t62 * t82;
t199 = t62 * t84;
t196 = -qJ(6) - pkin(9);
t125 = pkin(5) * t111 - qJ(6) * t169;
t156 = t113 * qJD(6);
t174 = qJ(6) * t111;
t86 = t99 * t169;
t195 = -t111 * t156 + t125 * qJD(4) + (-t86 + (-t79 + t174) * t110) * qJD(5) + t212;
t170 = t111 * t113;
t194 = (-qJ(6) * qJD(5) - qJD(4) * t99) * t170 + (-qJD(6) * t111 + (-qJ(6) * qJD(4) - qJD(5) * t99) * t114) * t110 + t211;
t139 = qJD(5) * t196;
t87 = t133 * qJD(2);
t141 = -t110 * t210 + t113 * t87;
t193 = qJD(2) * t125 + t110 * qJD(6) - t113 * t139 + t141;
t149 = t110 * t163;
t191 = t110 * t87 + t113 * t210;
t192 = -qJ(6) * t149 - t110 * t139 - t156 + t191;
t190 = -t82 * t148 - t59 * t170;
t187 = t110 * t79 + t86;
t63 = qJD(2) * t65;
t57 = qJD(1) * t63;
t20 = -t111 * t57 + t51 * t160 + t97 * t161;
t10 = pkin(5) * t59 + t20;
t186 = t10 * t110;
t185 = t10 * t113;
t32 = -qJD(4) * pkin(4) - t210;
t184 = t110 * t32;
t181 = t111 * t82;
t180 = t113 * t32;
t179 = t113 * t98;
t178 = t114 * t59;
t176 = t20 * t110;
t175 = t20 * t113;
t173 = qJD(5) * t82;
t117 = qJD(2) ^ 2;
t172 = t107 * t117;
t116 = qJD(4) ^ 2;
t168 = t116 * t111;
t167 = t116 * t114;
t166 = -t114 ^ 2 + t104;
t155 = qJD(2) * qJD(4);
t153 = t84 * t160;
t152 = t98 * t159;
t144 = t111 * t155;
t143 = -qJD(6) - t204;
t19 = qJD(4) * t210 - t114 * t57;
t66 = (t106 * t115 + t108 * t112) * t107;
t42 = (qJD(1) * t66 + t88) * qJD(2);
t142 = t110 * t19 - t113 * t42;
t140 = t58 * t114 + t84 * t161;
t138 = -t110 * t42 - t113 * t19 - t41 * t158 + t33 * t159;
t137 = pkin(5) * t144;
t135 = t84 * t146;
t134 = t98 * t146;
t132 = -t110 * t9 - t113 * t7;
t131 = t110 * t7 - t113 * t9;
t49 = t109 * t111 + t114 * t66;
t26 = t110 * t65 + t113 * t49;
t25 = -t110 * t49 + t113 * t65;
t129 = t109 * t114 - t111 * t66;
t61 = qJD(2) * t66;
t56 = qJD(1) * t61;
t124 = qJD(2) * t60 - t116 * t99 - t56;
t123 = qJ(6) * t59 + t138;
t50 = -qJD(2) * pkin(3) - t54;
t122 = qJD(4) * (qJD(2) * (-pkin(3) - t202) + t50 + t62);
t121 = t128 * t110;
t119 = -qJD(5) * t12 - t142;
t118 = qJ(6) * t58 + t119;
t102 = -pkin(5) * t113 - pkin(4);
t95 = t196 * t113;
t94 = t196 * t110;
t81 = t82 ^ 2;
t73 = (t99 + t201) * t111;
t69 = t113 * t79;
t53 = pkin(5) * t120 + t160 * t99;
t43 = -t110 * t174 + t187;
t40 = -qJ(6) * t170 + t69 + (-pkin(5) - t182) * t114;
t29 = pkin(5) * t149 + t36;
t24 = qJD(4) * t49 - t111 * t63;
t23 = qJD(4) * t129 - t114 * t63;
t22 = -t143 + t32;
t16 = t140 + t207;
t15 = t134 - t178 + (-t121 + t181) * qJD(4);
t6 = -qJD(5) * t26 - t110 * t23 + t113 * t61;
t5 = qJD(5) * t25 + t110 * t61 + t113 * t23;
t4 = -qJD(6) * t82 - t123;
t3 = -qJD(6) * t84 + t118 + t137;
t2 = t129 * t58 - t144 * t26 + t24 * t84 + t5 * t98;
t1 = -t129 * t59 + t144 * t25 + t24 * t82 - t6 * t98;
t13 = [0, 0, -t112 * t172, -t115 * t172, -t54 * t61 - t55 * t63 + t56 * t65 - t57 * t66, 0, 0, 0, 0, 0, -t24 * qJD(4) + (-t114 * t61 + t161 * t65) * qJD(2), -t23 * qJD(4) + (t111 * t61 + t160 * t65) * qJD(2), 0, 0, 0, 0, 0, t1, t2, t1, t2, t25 * t58 - t26 * t59 - t5 * t82 - t6 * t84, -t10 * t129 + t22 * t24 + t25 * t3 + t26 * t4 + t5 * t9 + t6 * t7; 0, 0, 0, 0, t54 * t60 - t55 * t62 + (-t106 * t57 - t108 * t56) * pkin(2), 0.2e1 * t114 * t144, -0.2e1 * t166 * t155, t167, -t168, 0, t111 * t122 + t114 * t124, -t111 * t124 + t114 * t122, t84 * t148 + (-t113 * t58 - t159 * t84) * t111, -t135 + (-t153 + (t58 + t173) * t111) * t110 + t190, t140 - t207, t134 + t178 + (-t121 - t181) * qJD(4), (-t98 - t163) * t161, (t159 * t79 - t212) * t98 + ((t82 * t99 + t184) * qJD(4) + (t183 + (t98 * t99 + t33) * t113) * qJD(5) + t142) * t114 + (t32 * t158 + t176 + t99 * t59 - t200 + ((-t99 * t171 + t69) * qJD(2) + t11) * qJD(4)) * t111, t211 * t98 + (-t99 * t152 + (t84 * t99 + t180) * qJD(4) - t138) * t114 + (-t32 * t159 + t175 - t99 * t58 - t199 + (-t187 * qJD(2) - t99 * t179 - t12) * qJD(4)) * t111, t53 * t82 + t59 * t73 - t195 * t98 + (t162 * t22 - t3) * t114 + (t22 * t158 + t186 - t200 + (qJD(2) * t40 + t7) * qJD(4)) * t111, t53 * t84 - t58 * t73 + t194 * t98 + (t157 * t22 + t4) * t114 + (-t22 * t159 + t185 - t199 + (-qJD(2) * t43 - t9) * qJD(4)) * t111, t40 * t58 - t43 * t59 - t195 * t84 - t194 * t82 + t132 * t160 + (qJD(5) * t131 - t110 * t4 - t113 * t3) * t111, t10 * t73 + t3 * t40 + t4 * t43 + t194 * t9 + t195 * t7 + (-t111 * t62 + t53) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t168, -t167, 0, 0, 0, 0, 0, t15, t16, t15, t16, t135 + (t153 + (-t58 + t173) * t111) * t110 + t190, (-qJD(4) * t131 - t10) * t114 + (qJD(4) * t22 + qJD(5) * t132 - t110 * t3 + t113 * t4) * t111; 0, 0, 0, 0, 0, -t111 * t117 * t114, t166 * t117, 0, 0, 0, qJD(4) * t36 - t164 * t50 - t20, (-qJD(2) * t50 + t57) * t114, -t58 * t110 - t84 * t179, -t208 * t110 + t209 * t113, -t98 * t158 + (t98 * t169 + (-t84 + t162) * t111) * qJD(2), t152 + (-t98 * t171 + (t82 + t157) * t111) * qJD(2), t98 * t164, -pkin(4) * t59 - t175 + t141 * t98 - t36 * t82 + (pkin(9) * t179 + t184) * qJD(5) + (-t11 * t111 + (-pkin(9) * t161 - t114 * t32) * t110) * qJD(2), pkin(4) * t58 + t176 - t36 * t84 - t191 * t98 + (-pkin(9) * t110 * t98 + t180) * qJD(5) + (-t32 * t169 + (-pkin(9) * t157 + t12) * t111) * qJD(2), -t185 + t102 * t59 - t29 * t82 + t193 * t98 + (t22 + t204) * t159 + (-t22 * t171 + (qJD(4) * t94 - t7) * t111) * qJD(2), t186 - t102 * t58 - t29 * t84 - t192 * t98 + (t113 * t22 + t84 * t201) * qJD(5) + (-t22 * t169 + (qJD(4) * t95 + t9) * t111) * qJD(2), t58 * t94 + t59 * t95 + t193 * t84 + t192 * t82 + (t7 * t98 + t4) * t113 + (-t3 + t203) * t110, t10 * t102 + t3 * t94 - t4 * t95 - t192 * t9 - t193 * t7 + (pkin(5) * t159 - t29) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84 * t82, -t81 + t206, -t58 - t198, -t197 - t59, t144, -t12 * t98 - t32 * t84 + t119, -t11 * t98 + t32 * t82 + t138, 0.2e1 * t137 - t203 + (t143 - t22) * t84 + t118, -pkin(5) * t206 - t8 * t98 + (qJD(6) + t22) * t82 + t123, t58 * pkin(5) - t205 * t82, t205 * t9 + (-t22 * t84 + t3) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, t209, -t81 - t206, t7 * t84 + t82 * t9 + t10;];
tauc_reg = t13;
