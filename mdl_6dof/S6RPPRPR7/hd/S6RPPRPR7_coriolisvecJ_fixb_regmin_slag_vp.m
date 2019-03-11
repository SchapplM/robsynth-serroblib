% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
% 
% Output:
% tauc_reg [6x28]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPPRPR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:53:44
% EndTime: 2019-03-09 01:53:49
% DurationCPUTime: 1.94s
% Computational Cost: add. (2994->286), mult. (6904->403), div. (0->0), fcn. (5164->8), ass. (0->151)
t129 = sin(pkin(9));
t131 = cos(pkin(9));
t134 = sin(qJ(4));
t136 = cos(qJ(4));
t102 = t129 * t136 + t131 * t134;
t142 = qJD(1) * t102;
t90 = qJD(6) + t142;
t208 = qJD(6) - t90;
t130 = cos(pkin(10));
t133 = sin(qJ(6));
t128 = sin(pkin(10));
t135 = cos(qJ(6));
t176 = t128 * t135;
t103 = t130 * t133 + t176;
t96 = t103 * qJD(6);
t189 = t103 * t142 + t96;
t171 = qJD(1) * t129;
t161 = t134 * t171;
t174 = t136 * t131;
t163 = qJD(1) * t174;
t203 = t161 - t163;
t74 = -qJD(4) * t130 - t128 * t203;
t207 = t135 * t74;
t206 = t142 * t74;
t76 = qJD(4) * t128 - t130 * t203;
t205 = t142 * t76;
t150 = t133 * t74 - t135 * t76;
t204 = t150 * t90;
t144 = t103 * t90;
t132 = -pkin(1) - qJ(3);
t199 = t132 * qJD(1);
t111 = qJD(2) + t199;
t160 = -pkin(7) * qJD(1) + t111;
t87 = t160 * t129;
t88 = t160 * t131;
t202 = -t134 * t87 + t136 * t88;
t193 = -pkin(7) + t132;
t104 = t193 * t129;
t105 = t193 * t131;
t70 = t104 * t134 - t105 * t136;
t201 = qJD(4) * t142;
t175 = t135 * t130;
t177 = t128 * t133;
t100 = -t175 + t177;
t200 = t100 * qJD(6);
t172 = t129 ^ 2 + t131 ^ 2;
t198 = t172 * qJD(3);
t110 = qJD(4) * t161;
t169 = qJD(4) * t136;
t162 = t131 * t169;
t86 = qJD(1) * t162 - t110;
t186 = t103 * t86;
t190 = -t100 * t142 - t200;
t197 = -t190 * t90 - t186;
t14 = -qJD(6) * t150 - t103 * t201;
t91 = t142 ^ 2;
t126 = qJD(1) * qJD(2);
t167 = 0.2e1 * t126;
t196 = pkin(8) * t130;
t38 = t133 * t76 + t207;
t195 = t38 * t203;
t194 = t150 * t203;
t192 = pkin(8) + qJ(5);
t143 = t102 * qJD(3);
t32 = -qJD(1) * t143 + (qJD(5) + t202) * qJD(4);
t37 = pkin(4) * t86 + qJ(5) * t201 + qJD(5) * t203 + t126;
t9 = t128 * t37 + t130 * t32;
t101 = t129 * t134 - t174;
t170 = qJD(4) * t134;
t97 = -t129 * t169 - t131 * t170;
t98 = -t129 * t170 + t162;
t43 = pkin(4) * t98 - qJ(5) * t97 + qJD(5) * t101 + qJD(2);
t49 = -qJD(4) * t70 - t143;
t16 = t128 * t43 + t130 * t49;
t58 = t134 * t88 + t136 * t87;
t53 = qJD(4) * qJ(5) + t58;
t127 = qJD(1) * qJ(2);
t121 = qJD(3) + t127;
t106 = pkin(3) * t171 + t121;
t54 = pkin(4) * t142 + qJ(5) * t203 + t106;
t20 = t128 * t54 + t130 * t53;
t66 = -pkin(4) * t203 + qJ(5) * t142;
t24 = t128 * t66 + t130 * t202;
t116 = pkin(3) * t129 + qJ(2);
t65 = pkin(4) * t102 + qJ(5) * t101 + t116;
t71 = t104 * t136 + t105 * t134;
t29 = t128 * t65 + t130 * t71;
t168 = qJD(6) * t135;
t191 = -t168 * t74 - t175 * t201;
t69 = t100 * t86;
t188 = t101 * t201;
t187 = t102 * t86;
t185 = t128 * t201;
t184 = t128 * t142;
t183 = t128 * t97;
t182 = t130 * t201;
t180 = qJD(1) * t142;
t179 = t101 * t128;
t173 = t98 * qJD(4);
t8 = -t128 * t32 + t130 * t37;
t4 = pkin(5) * t86 + pkin(8) * t182 + t8;
t5 = pkin(8) * t185 + t9;
t166 = -t133 * t5 + t135 * t4;
t19 = -t128 * t53 + t130 * t54;
t165 = -t142 * t19 + t9;
t164 = t142 * t20 + t8;
t15 = -t128 * t49 + t130 * t43;
t23 = -t128 * t202 + t130 * t66;
t28 = -t128 * t71 + t130 * t65;
t159 = qJD(1) * t172;
t158 = -t189 * t90 - t69;
t157 = t102 * t8 + t19 * t98;
t156 = -t102 * t9 - t20 * t98;
t155 = t133 * t4 + t135 * t5;
t139 = qJD(3) * t203 - t87 * t169 - t88 * t170;
t48 = -qJD(4) * pkin(4) + qJD(5) - t202;
t154 = -t101 * t139 - t48 * t97;
t11 = -pkin(8) * t74 + t20;
t7 = pkin(5) * t142 - pkin(8) * t76 + t19;
t2 = t11 * t135 + t133 * t7;
t153 = t11 * t133 - t135 * t7;
t18 = pkin(5) * t102 + t101 * t196 + t28;
t21 = pkin(8) * t179 + t29;
t152 = -t133 * t21 + t135 * t18;
t151 = t133 * t18 + t135 * t21;
t148 = -qJD(6) * t76 + t185;
t107 = t192 * t128;
t147 = pkin(8) * t184 - qJD(5) * t130 + qJD(6) * t107 + t24;
t108 = t192 * t130;
t146 = -pkin(5) * t203 + qJD(5) * t128 + qJD(6) * t108 + t142 * t196 + t23;
t145 = t90 * t100;
t141 = -t201 * t70 - t154;
t140 = -t142 * t98 - t187 - t188;
t138 = pkin(4) * t201 - qJ(5) * t86 + (-qJD(5) + t48) * t142;
t13 = t133 * t148 + t191;
t50 = -qJD(3) * t101 + qJD(4) * t71;
t137 = qJD(1) ^ 2;
t119 = -pkin(5) * t130 - pkin(4);
t89 = t97 * qJD(4);
t64 = t100 * t101;
t63 = t103 * t101;
t51 = -pkin(5) * t179 + t70;
t36 = -pkin(5) * t184 + t58;
t31 = pkin(5) * t183 + t50;
t27 = pkin(5) * t74 + t48;
t26 = t97 * t176 + qJD(6) * t101 * t177 + (-t101 * t168 + t133 * t97) * t130;
t25 = -t100 * t97 + t101 * t96;
t22 = -pkin(5) * t185 - t139;
t10 = -pkin(8) * t183 + t16;
t6 = pkin(5) * t98 - t196 * t97 + t15;
t1 = [0, 0, 0, 0, t167, qJ(2) * t167, t129 * t167, t131 * t167, 0.2e1 * qJD(3) * t159 (t121 + t127) * qJD(2) + (-t111 - t199) * t198, -t203 * t97 + t188, t101 * t86 + t102 * t201 - t142 * t97 + t203 * t98, t89, -t173, 0, 0.2e1 * qJD(2) * t142 - qJD(4) * t50 + t106 * t98 + t116 * t86, -qJD(4) * t49 + t106 * t97 - t116 * t201 + (-qJD(1) * t101 - t203) * qJD(2), t128 * t141 + t142 * t15 + t28 * t86 + t50 * t74 + t157, t130 * t141 - t142 * t16 - t29 * t86 + t50 * t76 + t156, -t15 * t76 - t16 * t74 + (t101 * t8 - t19 * t97 + t201 * t28) * t130 + (t101 * t9 - t20 * t97 + t201 * t29) * t128, -t139 * t70 + t15 * t19 + t16 * t20 + t28 * t8 + t29 * t9 + t48 * t50, t13 * t64 - t150 * t25, t13 * t63 - t14 * t64 + t150 * t26 - t25 * t38, t102 * t13 - t150 * t98 + t25 * t90 + t64 * t86, -t102 * t14 - t26 * t90 - t38 * t98 + t63 * t86, t90 * t98 + t187 (-t10 * t133 + t135 * t6) * t90 + t152 * t86 + t166 * t102 - t153 * t98 + t31 * t38 + t51 * t14 - t22 * t63 + t27 * t26 + (-t102 * t2 - t151 * t90) * qJD(6) -(t10 * t135 + t133 * t6) * t90 - t151 * t86 - t155 * t102 - t2 * t98 - t31 * t150 + t51 * t13 + t22 * t64 + t27 * t25 + (t102 * t153 - t152 * t90) * qJD(6); 0, 0, 0, 0, -t137, -t137 * qJ(2), -t137 * t129, -t137 * t131, 0 (-t121 - t198) * qJD(1), 0, 0, 0, 0, 0, t89 - t180, qJD(1) * t203 - t173, t128 * t140 - t130 * t180 - t74 * t97, t128 * t180 + t130 * t140 - t76 * t97 (t128 * t76 - t130 * t74) * t98 + (t128 * t74 + t130 * t76) * qJD(1) (-qJD(1) * t19 - t156) * t130 + (-qJD(1) * t20 - t157) * t128 + t154, 0, 0, 0, 0, 0, t101 * t14 - t97 * t38 - t98 * t144 + qJD(1) * t145 + (t200 * t90 - t186) * t102, t101 * t13 + t97 * t150 + t98 * t145 + qJD(1) * t144 + (qJD(6) * t144 + t69) * t102; 0, 0, 0, 0, 0, 0, 0, 0, -t172 * t137, t111 * t159 + t126, 0, 0, 0, 0, 0, -t110 + (-t203 + t163) * qJD(4), -0.2e1 * t201, -t128 * t91 + t130 * t86 + t203 * t74, -t128 * t86 - t130 * t91 + t203 * t76 (t182 - t206) * t130 + (t185 + t205) * t128, t128 * t165 + t130 * t164 + t203 * t48, 0, 0, 0, 0, 0, t158 + t195, -t194 + t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t203 * t142, t203 ^ 2 - t91, 0, t110 + (-t203 - t163) * qJD(4), 0, qJD(4) * t58 + t106 * t203 + t139 (qJD(3) + t106) * t142, t128 * t138 + t130 * t139 - t142 * t23 + t19 * t203 - t58 * t74, -t128 * t139 + t130 * t138 + t142 * t24 - t20 * t203 - t58 * t76, t23 * t76 + t24 * t74 + (-qJD(5) * t74 + t165) * t130 + (qJD(5) * t76 - t164) * t128, pkin(4) * t139 - t19 * t23 - t20 * t24 - t48 * t58 + (-t128 * t19 + t130 * t20) * qJD(5) + (-t128 * t8 + t130 * t9) * qJ(5), t103 * t13 - t150 * t190, -t100 * t13 - t103 * t14 + t150 * t189 - t190 * t38, -t194 - t197, t158 - t195, t90 * t203 (-t107 * t135 - t108 * t133) * t86 + t119 * t14 + t22 * t100 - t153 * t203 - t36 * t38 + (t133 * t147 - t135 * t146) * t90 + t189 * t27 -(-t107 * t133 + t108 * t135) * t86 + t119 * t13 + t22 * t103 - t2 * t203 + t36 * t150 + (t133 * t146 + t135 * t147) * t90 + t190 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t185 + t205, -t182 - t206, -t74 ^ 2 - t76 ^ 2, t19 * t76 + t20 * t74 - t139, 0, 0, 0, 0, 0, t14 - t204, -t90 * t207 + (-t76 * t90 + t148) * t133 + t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t150 * t38, t150 ^ 2 - t38 ^ 2, t38 * t90 + t13, -t14 - t204, t86, t150 * t27 - t2 * t208 + t166, t153 * t208 + t27 * t38 - t155;];
tauc_reg  = t1;
