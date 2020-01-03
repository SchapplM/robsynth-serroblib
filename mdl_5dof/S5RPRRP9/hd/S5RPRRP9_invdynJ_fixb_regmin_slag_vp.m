% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP9
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tau_reg [5x25]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP9_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP9_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP9_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP9_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:49:14
% EndTime: 2019-12-31 18:49:18
% DurationCPUTime: 1.66s
% Computational Cost: add. (3006->269), mult. (7185->328), div. (0->0), fcn. (5517->12), ass. (0->156)
t124 = sin(qJ(4));
t127 = cos(qJ(3));
t122 = sin(pkin(8));
t123 = cos(pkin(8));
t125 = sin(qJ(3));
t175 = qJD(1) * qJD(3);
t165 = t127 * t175;
t166 = t125 * t175;
t174 = t122 * qJDD(1);
t172 = -t122 * t165 - t123 * t166 - t125 * t174;
t173 = t123 * qJDD(1);
t144 = t127 * t173 + t172;
t204 = cos(qJ(4));
t167 = qJD(4) * t204;
t178 = qJD(4) * t124;
t171 = t123 * t165 + t125 * t173 + t127 * t174;
t51 = -t122 * t166 + t171;
t182 = t127 * t123;
t183 = t125 * t122;
t75 = -t182 + t183;
t68 = t75 * qJD(1);
t76 = t127 * t122 + t125 * t123;
t69 = t76 * qJD(1);
t15 = -t124 * t144 + t68 * t167 + t69 * t178 - t204 * t51;
t121 = qJD(3) + qJD(4);
t43 = t124 * t69 + t204 * t68;
t190 = t43 * t121;
t7 = -t15 + t190;
t201 = t43 ^ 2;
t146 = -t124 * t68 + t204 * t69;
t216 = t146 ^ 2;
t221 = -t201 + t216;
t169 = t123 * pkin(2) + pkin(1);
t81 = -t169 * qJD(1) + qJD(2);
t54 = t68 * pkin(3) + t81;
t20 = t43 * pkin(4) - qJ(5) * t146 + t54;
t220 = t20 * t43;
t219 = t54 * t43;
t200 = t146 * t43;
t120 = pkin(8) + qJ(3);
t113 = cos(t120);
t114 = qJ(4) + t120;
t103 = sin(t114);
t104 = cos(t114);
t189 = t104 * pkin(4) + t103 * qJ(5);
t218 = pkin(3) * t113 + t189;
t188 = qJDD(1) * pkin(1);
t126 = sin(qJ(1));
t128 = cos(qJ(1));
t212 = g(1) * t126 - g(2) * t128;
t149 = -qJDD(2) + t188 + t212;
t16 = t146 * qJD(4) + t124 * t51 - t204 * t144;
t191 = t146 * t121;
t217 = -t16 + t191;
t27 = pkin(4) * t146 + t43 * qJ(5);
t199 = pkin(6) + qJ(2);
t89 = t199 * t122;
t77 = qJD(1) * t89;
t90 = t199 * t123;
t78 = qJD(1) * t90;
t214 = -t125 * t78 - t127 * t77;
t116 = qJDD(3) + qJDD(4);
t108 = t116 * qJ(5);
t109 = t121 * qJD(5);
t213 = t108 + t109;
t211 = qJ(2) * qJDD(1);
t111 = t116 * pkin(4);
t210 = t111 - qJDD(5);
t150 = t125 * t77 - t127 * t78;
t176 = qJD(1) * qJD(2);
t207 = t199 * qJDD(1) + t176;
t58 = t207 * t122;
t59 = t207 * t123;
t160 = -t125 * t59 - t127 * t58;
t12 = qJDD(3) * pkin(3) - t51 * pkin(7) + t150 * qJD(3) + t160;
t151 = -t125 * t58 + t127 * t59;
t14 = t144 * pkin(7) + t214 * qJD(3) + t151;
t37 = -t69 * pkin(7) + t214;
t36 = qJD(3) * pkin(3) + t37;
t38 = -t68 * pkin(7) - t150;
t164 = -t204 * t12 + t124 * t14 + t38 * t167 + t36 * t178;
t186 = t103 * t128;
t187 = t103 * t126;
t142 = g(1) * t186 + g(2) * t187 - g(3) * t104 - t164;
t133 = t146 * t20 - t142 - t210;
t209 = -t54 * t146 + t142;
t193 = t125 * t90;
t73 = t127 * t89;
t159 = -t73 - t193;
t205 = t76 * pkin(7);
t40 = t159 - t205;
t197 = -t125 * t89 + t127 * t90;
t41 = -t75 * pkin(7) + t197;
t26 = t124 * t40 + t204 * t41;
t147 = -t124 * t41 + t204 * t40;
t177 = t122 * qJD(2);
t198 = qJD(2) * t182 - qJD(3) * t73;
t30 = -t125 * t177 + (-t193 - t205) * qJD(3) + t198;
t131 = -t76 * qJD(2) - t197 * qJD(3);
t70 = t75 * qJD(3);
t31 = t70 * pkin(7) + t131;
t4 = t147 * qJD(4) + t124 * t31 + t204 * t30;
t208 = -t103 * t212 - t26 * t116 - t4 * t121;
t195 = t124 * t38;
t22 = t204 * t37 - t195;
t196 = pkin(3) * t167 + qJD(5) - t22;
t170 = t204 * t38;
t19 = t124 * t36 + t170;
t192 = t19 * t121;
t185 = t104 * t126;
t184 = t104 * t128;
t181 = t69 * qJD(3);
t18 = t204 * t36 - t195;
t180 = qJD(5) - t18;
t179 = t122 ^ 2 + t123 ^ 2;
t168 = qJD(1) * t183;
t163 = t124 * t12 + t204 * t14 + t36 * t167 - t38 * t178;
t158 = t179 * qJD(1) ^ 2;
t157 = 0.2e1 * t179;
t21 = t124 * t37 + t170;
t156 = pkin(3) * t178 - t21;
t112 = sin(t120);
t155 = -pkin(3) * t112 - pkin(4) * t103;
t154 = g(1) * t128 + g(2) * t126;
t56 = t75 * pkin(3) - t169;
t148 = t169 + t218;
t53 = -t124 * t75 + t204 * t76;
t145 = -g(1) * t184 - g(2) * t185 - g(3) * t103 + t163;
t143 = t76 * qJD(3);
t141 = pkin(3) * t143;
t140 = t149 + t188;
t139 = t18 * t121 - t145;
t138 = t15 + t190;
t5 = t26 * qJD(4) + t124 * t30 - t204 * t31;
t137 = g(1) * t185 - g(2) * t184 + t116 * t147 - t5 * t121;
t134 = t157 * t176 - t154;
t39 = qJDD(2) - t172 * pkin(3) + (-pkin(1) + (-pkin(3) * t127 - pkin(2)) * t123) * qJDD(1);
t132 = t16 + t191;
t3 = t16 * pkin(4) + t15 * qJ(5) - qJD(5) * t146 + t39;
t117 = -pkin(7) - t199;
t107 = -t204 * pkin(3) - pkin(4);
t105 = t124 * pkin(3) + qJ(5);
t83 = qJ(5) * t184;
t82 = qJ(5) * t185;
t80 = -t169 * qJDD(1) + qJDD(2);
t52 = t124 * t76 + t204 * t75;
t29 = t53 * qJD(4) - t124 * t70 + t204 * t143;
t28 = t124 * t143 + t75 * t167 + t76 * t178 + t204 * t70;
t24 = t52 * pkin(4) - t53 * qJ(5) + t56;
t23 = t69 * pkin(3) + t27;
t17 = t121 * qJ(5) + t19;
t13 = -t121 * pkin(4) + t180;
t6 = t29 * pkin(4) + t28 * qJ(5) - t53 * qJD(5) + t141;
t2 = t164 - t210;
t1 = t163 + t213;
t8 = [qJDD(1), t212, t154, t140 * t123, -t140 * t122, t157 * t211 + t134, t149 * pkin(1) + (t179 * t211 + t134) * qJ(2), t51 * t76 - t69 * t70, -t69 * t143 + t76 * t144 - t51 * t75 + t70 * t68, -t70 * qJD(3) + t76 * qJDD(3), -t76 * qJD(3) ^ 2 - t75 * qJDD(3), 0, t169 * t144 + t80 * t75 + t159 * qJDD(3) + t212 * t113 + (t81 * t76 + t131) * qJD(3), -t169 * t51 + t80 * t76 - t81 * t70 - ((-qJD(3) * t90 - t177) * t125 + t198) * qJD(3) - t197 * qJDD(3) - t212 * t112, -t146 * t28 - t15 * t53, -t146 * t29 + t15 * t52 - t53 * t16 + t28 * t43, t53 * t116 - t28 * t121, -t52 * t116 - t29 * t121, 0, t141 * t43 + t56 * t16 + t54 * t29 + t39 * t52 + t137, t141 * t146 - t56 * t15 - t54 * t28 + t39 * t53 + t208, t24 * t16 + t20 * t29 + t3 * t52 + t6 * t43 + t137, -t1 * t52 - t13 * t28 + t146 * t5 + t147 * t15 - t26 * t16 - t17 * t29 + t2 * t53 - t4 * t43 - t154, -t146 * t6 + t24 * t15 + t20 * t28 - t3 * t53 - t208, t1 * t26 + t13 * t5 + t17 * t4 - t2 * t147 + t20 * t6 + t3 * t24 + (g(1) * t117 - g(2) * t148) * t128 + (g(1) * t148 + g(2) * t117) * t126; 0, 0, 0, -t173, t174, -t158, -qJ(2) * t158 - t149, 0, 0, 0, 0, 0, -t144 + t181, (-t68 - t168) * qJD(3) + t171, 0, 0, 0, 0, 0, t132, -t138, t132, -t201 - t216, t138, -t13 * t146 + t17 * t43 - t212 + t3; 0, 0, 0, 0, 0, 0, 0, t69 * t68, -t68 ^ 2 + t69 ^ 2, (t68 - t168) * qJD(3) + t171, t144 + t181, qJDD(3), -g(3) * t113 + t154 * t112 - t81 * t69 + t160, g(3) * t112 + t154 * t113 + t81 * t68 - t151, t200, t221, t7, t217, t116, t21 * t121 + (t204 * t116 - t121 * t178 - t43 * t69) * pkin(3) + t209, t22 * t121 + t219 + (-t116 * t124 - t121 * t167 - t146 * t69) * pkin(3) - t145, -t107 * t116 - t121 * t156 - t23 * t43 - t133, -t105 * t16 - t107 * t15 + (t156 + t17) * t146 + (t13 - t196) * t43, t105 * t116 + t121 * t196 + t146 * t23 + t145 + t213 - t220, t1 * t105 + t2 * t107 - t20 * t23 - g(1) * (t128 * t155 + t83) - g(2) * (t126 * t155 + t82) - g(3) * t218 + t196 * t17 + t156 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t200, t221, t7, t217, t116, t192 + t209, t139 + t219, -t27 * t43 + t111 - t133 + t192, pkin(4) * t15 - t16 * qJ(5) + (t17 - t19) * t146 + (t13 - t180) * t43, t146 * t27 + 0.2e1 * t108 + 0.2e1 * t109 - t139 - t220, t1 * qJ(5) - t2 * pkin(4) - t20 * t27 - t13 * t19 - g(1) * (-pkin(4) * t186 + t83) - g(2) * (-pkin(4) * t187 + t82) - g(3) * t189 + t180 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116 + t200, t7, -t121 ^ 2 - t216, -t17 * t121 + t133;];
tau_reg = t8;
