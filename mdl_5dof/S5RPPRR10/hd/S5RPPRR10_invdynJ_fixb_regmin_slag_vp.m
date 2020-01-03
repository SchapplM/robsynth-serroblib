% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRR10
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRR10_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR10_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR10_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR10_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:04:20
% EndTime: 2019-12-31 18:04:23
% DurationCPUTime: 1.26s
% Computational Cost: add. (1107->224), mult. (2571->308), div. (0->0), fcn. (2029->10), ass. (0->138)
t118 = sin(qJ(5));
t121 = cos(qJ(5));
t164 = qJD(5) * t121;
t165 = qJD(5) * t118;
t119 = sin(qJ(4));
t122 = cos(qJ(4));
t116 = cos(pkin(8));
t159 = t116 * qJDD(1);
t115 = sin(pkin(8));
t99 = t115 * qJDD(1);
t148 = -t119 * t159 + t122 * t99;
t67 = t115 * t119 + t116 * t122;
t55 = t67 * qJD(4);
t26 = -qJD(1) * t55 + t148;
t171 = qJD(1) * t116;
t155 = t119 * t171;
t131 = qJD(4) * t155 - t67 * qJDD(1);
t166 = qJD(4) * t122;
t156 = t115 * t166;
t27 = qJD(1) * t156 - t131;
t57 = t67 * qJD(1);
t172 = qJD(1) * t115;
t157 = t122 * t172;
t59 = -t155 + t157;
t134 = t118 * t27 - t121 * t26 + t57 * t164 + t59 * t165;
t112 = qJD(4) + qJD(5);
t19 = t118 * t59 + t121 * t57;
t189 = t112 * t19;
t195 = -t134 + t189;
t194 = t112 ^ 2;
t144 = -t118 * t57 + t121 * t59;
t193 = t144 * t19;
t192 = t144 ^ 2 - t19 ^ 2;
t76 = qJ(2) * t172 + qJD(3);
t63 = -pkin(6) * t172 + t76;
t182 = -pkin(6) + qJ(2);
t74 = t182 * t116;
t69 = qJD(1) * t74;
t142 = -t119 * t63 - t122 * t69;
t14 = -pkin(7) * t57 - t142;
t51 = -qJD(1) * pkin(1) - pkin(2) * t171 - qJ(3) * t172 + qJD(2);
t35 = pkin(3) * t171 - t51;
t17 = pkin(4) * t57 + t35;
t162 = qJD(1) * qJD(2);
t61 = qJ(2) * t99 + t115 * t162 + qJDD(3);
t40 = -pkin(6) * t99 + t61;
t42 = (t182 * qJDD(1) + t162) * t116;
t151 = -t119 * t42 + t122 * t40;
t2 = qJDD(4) * pkin(4) - pkin(7) * t26 + qJD(4) * t142 + t151;
t120 = sin(qJ(1));
t113 = qJ(4) + qJ(5);
t101 = sin(t113);
t102 = cos(t113);
t141 = t101 * t115 + t102 * t116;
t37 = t141 * t120;
t123 = cos(qJ(1));
t39 = t141 * t123;
t175 = t102 * t115;
t53 = t101 * t116 - t175;
t191 = t17 * t19 + t14 * t165 + g(1) * t39 + (-t14 * t112 - t2) * t118 + g(2) * t37 - g(3) * t53;
t188 = -t119 * t69 + t122 * t63;
t187 = t112 * t144;
t108 = g(2) * t123;
t183 = g(1) * t120;
t153 = -t108 + t183;
t161 = qJD(1) * qJD(3);
t89 = t115 * t161;
t186 = -pkin(2) * t159 - qJ(3) * t99 - t89;
t13 = -pkin(7) * t59 + t188;
t12 = qJD(4) * pkin(4) + t13;
t177 = t121 * t14;
t147 = -t118 * t12 - t177;
t143 = t119 * t40 + t122 * t42;
t3 = -pkin(7) * t27 + t188 * qJD(4) + t143;
t174 = t116 * t120;
t36 = t101 * t174 - t120 * t175;
t38 = t53 * t123;
t185 = g(1) * t38 + g(2) * t36 + g(3) * t141 + t147 * qJD(5) - t118 * t3 + t121 * t2 - t17 * t144;
t5 = t144 * qJD(5) + t118 * t26 + t121 * t27;
t184 = -t5 + t187;
t105 = t116 * pkin(2);
t73 = t182 * t115;
t181 = t119 * t73 + t122 * t74;
t111 = t116 ^ 2;
t160 = qJDD(1) * qJ(2) ^ 2;
t158 = 0.2e1 * t162;
t85 = t111 * t158;
t180 = qJ(2) * t85 + t111 * t160;
t114 = qJDD(1) * pkin(1);
t173 = t123 * pkin(1) + t120 * qJ(2);
t170 = qJD(2) * t119;
t169 = qJD(2) * t122;
t168 = qJD(3) * t115;
t167 = qJD(4) * t119;
t163 = qJ(2) * qJDD(1);
t98 = qJDD(2) - t114;
t70 = -t115 * qJ(3) - pkin(1) - t105;
t154 = -pkin(1) * t120 + t123 * qJ(2);
t150 = -t119 * t74 + t122 * t73;
t110 = t115 ^ 2;
t125 = qJD(1) ^ 2;
t71 = (-t110 - t111) * t125;
t60 = t116 * pkin(3) - t70;
t149 = g(1) * t123 + g(2) * t120;
t33 = t98 + t186;
t68 = t115 * t122 - t116 * t119;
t15 = -pkin(7) * t68 + t150;
t16 = -pkin(7) * t67 + t181;
t146 = -t118 * t16 + t121 * t15;
t145 = t118 * t15 + t121 * t16;
t28 = t118 * t68 + t121 * t67;
t29 = -t118 * t67 + t121 * t68;
t140 = t118 * t122 + t119 * t121;
t139 = -t118 * t119 + t121 * t122;
t138 = -t153 + t98;
t30 = pkin(3) * t159 - t33;
t137 = t114 - t98 - t108;
t136 = -qJDD(1) * t70 - t108 - t33;
t135 = t115 * t170 + t116 * t169 + t73 * t166 - t74 * t167;
t133 = 0.2e1 * t111 * t163 - t149 + t85;
t132 = (t162 + t163) * t110;
t129 = -t181 * qJD(4) + t115 * t169 - t116 * t170;
t124 = qJD(4) ^ 2;
t109 = qJDD(4) + qJDD(5);
t86 = g(1) * t174;
t56 = -t116 * t167 + t156;
t49 = t67 * t123;
t48 = t68 * t123;
t47 = t67 * t120;
t46 = t68 * t120;
t34 = pkin(4) * t56 + t168;
t31 = pkin(4) * t67 + t60;
t10 = pkin(7) * t55 + t129;
t9 = -pkin(7) * t56 + t135;
t8 = pkin(4) * t27 + t30;
t7 = qJD(5) * t29 - t118 * t55 + t121 * t56;
t6 = -qJD(5) * t28 - t118 * t56 - t121 * t55;
t1 = [qJDD(1), t153, t149, t116 * t137 + t86, (-t137 - t183) * t115, 0.2e1 * t132 + t133, -t98 * pkin(1) - g(1) * t154 - g(2) * t173 + (qJ(2) * t158 + t160) * t110 + t180, t86 + (t136 + t89) * t116, t115 * t61 + t132 + t133, t110 * t161 + (t136 + t183) * t115, t33 * t70 - g(1) * (-pkin(2) * t174 + t154) - g(2) * (t123 * t105 + t173) + (t61 * qJ(2) + t153 * qJ(3) + t76 * qJD(2) - t51 * qJD(3)) * t115 + t180, t26 * t68 - t55 * t59, -t26 * t67 - t27 * t68 + t55 * t57 - t56 * t59, -qJD(4) * t55 + qJDD(4) * t68, -qJD(4) * t56 - qJDD(4) * t67, 0, g(1) * t47 - g(2) * t49 + t129 * qJD(4) + t150 * qJDD(4) + t57 * t168 + t60 * t27 + t30 * t67 + t35 * t56, g(1) * t46 - g(2) * t48 - t135 * qJD(4) - t181 * qJDD(4) + t59 * t168 + t60 * t26 + t30 * t68 - t35 * t55, -t134 * t29 + t144 * t6, t134 * t28 - t144 * t7 - t19 * t6 - t29 * t5, t109 * t29 + t112 * t6, -t109 * t28 - t112 * t7, 0, t34 * t19 + t31 * t5 + t8 * t28 + t17 * t7 + (-qJD(5) * t145 + t10 * t121 - t118 * t9) * t112 + t146 * t109 + g(1) * t37 - g(2) * t39, t34 * t144 - t31 * t134 + t8 * t29 + t17 * t6 - (qJD(5) * t146 + t10 * t118 + t121 * t9) * t112 - t145 * t109 - g(1) * t36 + g(2) * t38; 0, 0, 0, -t159, t99, t71, qJ(2) * t71 + t138, -t159, t71, -t99, -qJ(2) * t111 * t125 - t76 * t172 + t138 + t186, 0, 0, 0, 0, 0, (-t59 - t157) * qJD(4) + t131, 0.2e1 * qJD(4) * t57 - t148, 0, 0, 0, 0, 0, -t5 - t187, t134 + t189; 0, 0, 0, 0, 0, 0, 0, -t115 * t125 * t116, t99, -t110 * t125, g(3) * t116 + (qJD(1) * t51 - t149) * t115 + t61, 0, 0, 0, 0, 0, qJDD(4) * t122 - t119 * t124 - t57 * t172, -t119 * qJDD(4) - t122 * t124 - t59 * t172, 0, 0, 0, 0, 0, t139 * t109 - t194 * t140 - t19 * t172, -t140 * t109 - t194 * t139 - t144 * t172; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59 * t57, -t57 ^ 2 + t59 ^ 2, t148, (t59 - t157) * qJD(4) + t131, qJDD(4), -g(1) * t48 - g(2) * t46 + g(3) * t67 - t35 * t59 + t151, g(1) * t49 + g(2) * t47 + g(3) * t68 + t35 * t57 - t143, t193, t192, t195, t184, t109, -(-t118 * t13 - t177) * t112 + (t121 * t109 - t112 * t165 - t59 * t19) * pkin(4) + t185, (-qJD(5) * t12 + t13 * t112 - t3) * t121 + (-t118 * t109 - t112 * t164 - t144 * t59) * pkin(4) + t191; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t193, t192, t195, t184, t109, -t112 * t147 + t185, (-t3 + (-qJD(5) + t112) * t12) * t121 + t191;];
tau_reg = t1;
