% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:55
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP5_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP5_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:54:57
% EndTime: 2019-12-31 19:55:03
% DurationCPUTime: 1.99s
% Computational Cost: add. (3624->260), mult. (9636->325), div. (0->0), fcn. (7049->6), ass. (0->152)
t133 = sin(qJ(4));
t131 = sin(pkin(8));
t132 = cos(pkin(8));
t135 = cos(qJ(2));
t181 = t132 * t135;
t167 = qJD(1) * t181;
t134 = sin(qJ(2));
t173 = t134 * qJD(1);
t101 = -t131 * t173 + t167;
t108 = t131 * t135 + t132 * t134;
t150 = qJD(1) * t108;
t204 = cos(qJ(4));
t153 = -t101 * t133 - t150 * t204;
t88 = qJD(2) * t150;
t172 = qJD(2) * qJD(1);
t165 = t134 * t172;
t119 = t131 * t165;
t89 = qJD(2) * t167 - t119;
t142 = qJD(4) * t153 - t133 * t89 - t204 * t88;
t128 = qJD(2) + qJD(4);
t210 = t128 * t153;
t217 = t142 - t210;
t166 = qJD(4) * t204;
t174 = qJD(4) * t133;
t152 = t101 * t166 - t133 * t88 - t150 * t174 + t204 * t89;
t53 = t101 * t204 - t133 * t150;
t184 = t53 * t128;
t216 = t152 - t184;
t197 = t53 ^ 2;
t215 = t153 ^ 2;
t156 = t215 - t197;
t124 = -t135 * pkin(2) - pkin(1);
t175 = qJD(1) * t124;
t115 = qJD(3) + t175;
t62 = -t101 * pkin(3) + t115;
t17 = -pkin(4) * t53 + qJ(5) * t153 + t62;
t214 = t17 * t53;
t213 = t62 * t53;
t212 = t153 * t62;
t199 = t17 * t153;
t211 = t53 * t153;
t26 = -pkin(4) * t153 - qJ(5) * t53;
t209 = -0.2e1 * t172;
t123 = pkin(2) * t132 + pkin(3);
t201 = t101 * pkin(7);
t191 = -qJ(3) - pkin(6);
t117 = t191 * t135;
t113 = qJD(1) * t117;
t106 = t132 * t113;
t116 = t191 * t134;
t112 = qJD(1) * t116;
t60 = -t112 * t131 + t106;
t148 = t60 - t201;
t203 = pkin(2) * t131;
t171 = t133 * t203;
t200 = t150 * pkin(7);
t104 = t131 * t113;
t61 = t112 * t132 + t104;
t44 = t61 - t200;
t188 = qJD(4) * t171 - t123 * t166 + t133 * t148 + t204 * t44;
t208 = 0.2e1 * t150;
t96 = t123 * t133 + t203 * t204;
t187 = qJD(2) * pkin(2);
t107 = t112 + t187;
t56 = t107 * t132 + t104;
t39 = qJD(2) * pkin(3) - t200 + t56;
t57 = t107 * t131 - t106;
t41 = t57 + t201;
t16 = t133 * t39 + t204 * t41;
t162 = qJD(2) * t191;
t97 = t135 * qJD(3) + t134 * t162;
t75 = t97 * qJD(1);
t98 = -t134 * qJD(3) + t135 * t162;
t76 = t98 * qJD(1);
t42 = -t131 * t75 + t132 * t76;
t33 = -pkin(7) * t89 + t42;
t43 = t131 * t76 + t132 * t75;
t34 = -pkin(7) * t88 + t43;
t3 = qJD(4) * t16 + t133 * t34 - t204 * t33;
t206 = t150 ^ 2;
t136 = qJD(2) ^ 2;
t64 = t116 * t132 + t117 * t131;
t48 = -pkin(7) * t108 + t64;
t155 = t131 * t134 - t181;
t65 = t116 * t131 - t117 * t132;
t49 = -pkin(7) * t155 + t65;
t22 = t133 * t49 - t204 * t48;
t205 = t3 * t22;
t202 = pkin(2) * t134;
t149 = t155 * qJD(2);
t46 = -t131 * t97 + t132 * t98;
t141 = pkin(7) * t149 + t46;
t151 = t108 * qJD(2);
t47 = t131 * t98 + t132 * t97;
t37 = -pkin(7) * t151 + t47;
t7 = t133 * t141 + t166 * t48 - t174 * t49 + t204 * t37;
t193 = t7 * t128;
t23 = t133 * t48 + t204 * t49;
t8 = qJD(4) * t23 + t133 * t37 - t141 * t204;
t192 = t8 * t128;
t190 = -qJD(4) * t96 + t133 * t44 - t148 * t204;
t189 = -qJD(5) + t188;
t59 = t108 * t204 - t133 * t155;
t29 = qJD(4) * t59 - t133 * t149 + t151 * t204;
t186 = t128 * t29;
t182 = t150 * t101;
t137 = qJD(1) ^ 2;
t180 = t135 * t137;
t179 = t136 * t134;
t178 = t136 * t135;
t15 = -t133 * t41 + t204 * t39;
t177 = qJD(5) - t15;
t176 = t134 ^ 2 - t135 ^ 2;
t126 = t134 * t187;
t125 = pkin(2) * t173;
t169 = t134 * t180;
t168 = t190 * t153;
t122 = pkin(2) * t165;
t63 = pkin(3) * t88 + t122;
t69 = pkin(3) * t150 + t125;
t161 = -t133 * t33 - t166 * t39 + t174 * t41 - t204 * t34;
t160 = pkin(1) * t209;
t159 = t135 * t165;
t127 = t128 * qJD(5);
t1 = t127 - t161;
t147 = t204 * t155;
t58 = t108 * t133 + t147;
t158 = -t142 * t58 - t29 * t53;
t157 = -t215 - t197;
t154 = t128 * t15 + t161;
t95 = t123 * t204 - t171;
t146 = t142 * t23 + t152 * t22 - t153 * t8 + t3 * t59 + t53 * t7;
t145 = t152 + t184;
t70 = pkin(3) * t151 + t126;
t28 = t108 * t174 + t128 * t147 + t133 * t151;
t144 = t142 * t59 - t152 * t58 + t153 * t29 - t28 * t53;
t78 = pkin(3) * t155 + t124;
t5 = -pkin(4) * t142 - qJ(5) * t152 + qJD(5) * t153 + t63;
t140 = t16 * t128 - t3;
t139 = -t142 - t210;
t138 = t128 * t190 - t3;
t99 = t101 ^ 2;
t93 = -pkin(4) - t95;
t92 = qJ(5) + t96;
t27 = t28 * t128;
t21 = t58 * pkin(4) - t59 * qJ(5) + t78;
t20 = t26 + t69;
t14 = qJ(5) * t128 + t16;
t13 = -pkin(4) * t128 + t177;
t9 = t29 * pkin(4) + t28 * qJ(5) - t59 * qJD(5) + t70;
t6 = t152 * t59 + t153 * t28;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t159, t176 * t209, t178, -0.2e1 * t159, -t179, 0, -pkin(6) * t178 + t134 * t160, pkin(6) * t179 + t135 * t160, 0, 0, t89 * t108 - t149 * t150, -t108 * t88 - t89 * t155 + (-t101 * t155 - t108 * t150) * qJD(2), -t155 * t136, -t101 * t151 + t155 * t88, -t108 * t136, 0, t124 * t88 + (t115 * t108 + t46 + (qJD(1) * t155 - t101) * t202) * qJD(2), t124 * t89 + (-t115 * t155 + t202 * t208 - t47) * qJD(2), t47 * t101 - t65 * t88 - t43 * t155 - t46 * t150 - t64 * t89 - t42 * t108 + (-t108 * t57 + t155 * t56) * qJD(2), t42 * t64 + t43 * t65 + t56 * t46 + t57 * t47 + (t115 + t175) * t126, t6, t144, -t27, t158, -t186, 0, -t142 * t78 + t29 * t62 - t53 * t70 + t58 * t63 - t192, t152 * t78 - t153 * t70 - t28 * t62 + t59 * t63 - t193, t15 * t28 - t16 * t29 + t161 * t58 + t146, -t15 * t8 + t16 * t7 - t161 * t23 + t62 * t70 + t63 * t78 + t205, t6, -t27, -t144, 0, t186, t158, -t142 * t21 + t17 * t29 + t5 * t58 - t53 * t9 - t192, -t1 * t58 - t13 * t28 - t14 * t29 + t146, -t152 * t21 + t153 * t9 + t17 * t28 - t5 * t59 + t193, t1 * t23 + t13 * t8 + t14 * t7 + t17 * t9 + t21 * t5 + t205; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t169, t176 * t137, 0, t169, 0, 0, t137 * pkin(1) * t134, pkin(1) * t180, 0, 0, -t182, -t99 + t206, -t119 + (-t101 + t167) * qJD(2), t182, 0, 0, -qJD(2) * t60 + t101 * t125 - t115 * t150 + t42, qJD(2) * t61 - t101 * t115 - t125 * t150 - t43, (t57 + t60) * t150 + (t56 - t61) * t101 + (-t131 * t88 - t132 * t89) * pkin(2), -t56 * t60 - t57 * t61 + (-t115 * t173 + t131 * t43 + t132 * t42) * pkin(2), t211, t156, t216, -t211, t217, 0, t53 * t69 + t138 + t212, t128 * t188 + t153 * t69 + t161 - t213, -t152 * t95 - t153 * t16 + t96 * t142 + t168 + (t15 - t188) * t53, t15 * t190 - t16 * t188 - t161 * t96 - t3 * t95 - t62 * t69, t211, t216, -t156, 0, -t217, -t211, t20 * t53 + t138 + t199, -t14 * t153 + t152 * t93 + t92 * t142 + t168 + (-t13 - t189) * t53, -t128 * t189 - t153 * t20 + t1 + t214, t1 * t92 - t13 * t190 - t14 * t189 - t17 * t20 + t3 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208 * qJD(2), -t119 + (t101 + t167) * qJD(2), -t99 - t206, -t101 * t57 + t150 * t56 + t122, 0, 0, 0, 0, 0, 0, t139, t145, t157, -t15 * t153 - t16 * t53 + t63, 0, 0, 0, 0, 0, 0, t139, t157, -t145, t13 * t153 - t14 * t53 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, t156, t216, -t211, t217, 0, t140 + t212, t154 - t213, 0, 0, t211, t216, -t156, 0, -t217, -t211, t26 * t53 + t140 + t199, -pkin(4) * t152 + t142 * qJ(5) - (t14 - t16) * t153 - (t13 - t177) * t53, -t153 * t26 + 0.2e1 * t127 - t154 + t214, -t3 * pkin(4) + t1 * qJ(5) - t13 * t16 + t14 * t177 - t17 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t211, t216, -t128 ^ 2 - t215, -t14 * t128 - t199 + t3;];
tauc_reg = t2;
