% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRRRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% tau_reg [5x23]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRRRR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:25
% EndTime: 2020-01-03 12:13:29
% DurationCPUTime: 1.45s
% Computational Cost: add. (2163->238), mult. (2960->318), div. (0->0), fcn. (1907->16), ass. (0->173)
t122 = sin(qJ(3));
t127 = cos(qJ(3));
t128 = cos(qJ(2));
t187 = qJD(1) * t128;
t168 = qJD(2) * t187;
t123 = sin(qJ(2));
t181 = qJDD(1) * t123;
t139 = (t168 + t181) * pkin(1);
t213 = pkin(1) * t128;
t101 = qJDD(1) * t213;
t113 = qJDD(1) + qJDD(2);
t205 = pkin(1) * qJD(1);
t180 = t123 * t205;
t50 = pkin(2) * t113 - qJD(2) * t180 + t101;
t115 = qJD(1) + qJD(2);
t179 = pkin(1) * t187;
t69 = pkin(2) * t115 + t179;
t224 = -t122 * t50 - (qJD(3) * t69 + t139) * t127;
t105 = qJD(3) + t115;
t120 = sin(qJ(5));
t126 = cos(qJ(4));
t121 = sin(qJ(4));
t125 = cos(qJ(5));
t193 = t121 * t125;
t59 = t120 * t126 + t193;
t46 = t59 * t105;
t161 = qJD(3) * t180;
t186 = qJD(3) * t122;
t192 = t122 * t123;
t90 = pkin(1) * t192;
t160 = t122 * pkin(1) * t168 + qJDD(1) * t90 + t69 * t186 + (t161 - t50) * t127;
t104 = qJDD(3) + t113;
t211 = pkin(3) * t104;
t13 = t160 - t211;
t119 = qJ(1) + qJ(2);
t110 = qJ(3) + t119;
t95 = sin(t110);
t96 = cos(t110);
t157 = -g(2) * t96 - g(3) * t95;
t150 = -t13 + t157;
t185 = qJD(3) * t127;
t88 = t122 * t180;
t54 = t127 * t179 - t88;
t223 = pkin(2) * t185 - t54;
t221 = g(2) * t95 - g(3) * t96;
t191 = t123 * t127;
t151 = t122 * t128 + t191;
t53 = t151 * t205;
t159 = pkin(2) * t186 - t53;
t114 = qJD(4) + qJD(5);
t130 = qJD(4) ^ 2;
t97 = pkin(2) * t122 + pkin(8);
t212 = pkin(2) * t127;
t98 = -pkin(3) - t212;
t220 = t104 * t98 + t105 * t159 + t130 * t97;
t219 = -pkin(8) - pkin(9);
t99 = pkin(2) + t213;
t206 = pkin(1) * t191 + t122 * t99;
t52 = pkin(8) + t206;
t215 = -pkin(9) - t52;
t214 = -pkin(9) - t97;
t210 = pkin(3) * t105;
t209 = pkin(4) * t126;
t190 = t125 * t126;
t175 = t105 * t190;
t196 = t120 * t121;
t176 = t105 * t196;
t44 = -t175 + t176;
t208 = t46 * t44;
t184 = qJD(4) * t121;
t102 = pkin(4) * t184;
t207 = t102 + t159;
t43 = t122 * t69 + t127 * t180;
t204 = t105 * t43;
t118 = qJ(4) + qJ(5);
t106 = sin(t118);
t203 = t106 * t95;
t202 = t106 * t96;
t37 = pkin(8) * t105 + t43;
t172 = pkin(9) * t105 + t37;
t26 = t172 * t126;
t200 = t125 * t26;
t30 = t99 * t186 + (qJD(2) * t151 + t123 * t185) * pkin(1);
t199 = t30 * t105;
t197 = t105 * t121;
t194 = t121 * t104;
t189 = t126 * t104;
t116 = t121 ^ 2;
t188 = -t126 ^ 2 + t116;
t183 = qJD(4) * t126;
t182 = qJD(5) * t120;
t100 = -pkin(3) - t209;
t174 = qJD(4) * t219;
t173 = t105 * t183;
t107 = sin(t119);
t109 = cos(t119);
t171 = g(2) * t107 - g(3) * t109;
t170 = qJD(4) * t215;
t169 = qJD(4) * t214;
t144 = t105 * t184 - t189;
t10 = pkin(4) * t144 + t13;
t42 = t127 * t69 - t88;
t28 = t100 * t105 - t42;
t58 = -t190 + t196;
t34 = t114 * t58;
t167 = g(2) * t202 + g(3) * t203 + t10 * t59 - t28 * t34;
t166 = t127 * t99 - t90;
t36 = -t42 - t210;
t165 = -t150 * t121 + t36 * t183;
t164 = qJD(1) * (-qJD(2) + t115);
t163 = qJD(2) * (-qJD(1) - t115);
t82 = t122 * t161;
t162 = t82 + t221;
t51 = -pkin(3) - t166;
t158 = -t43 + t102;
t156 = t120 * t194 - t125 * t189;
t25 = t172 * t121;
t24 = qJD(4) * pkin(4) - t25;
t154 = -t120 * t24 - t200;
t38 = t215 * t121;
t111 = t126 * pkin(9);
t39 = t126 * t52 + t111;
t153 = -t120 * t39 + t125 * t38;
t152 = t120 * t38 + t125 * t39;
t149 = -g(2) * t109 - g(3) * t107 + t101;
t57 = t126 * t97 + t111;
t148 = qJD(5) * t57 + t223 * t121 - t126 * t169;
t85 = pkin(8) * t126 + t111;
t147 = qJD(5) * t85 - t121 * t42 - t126 * t174;
t56 = t214 * t121;
t146 = -qJD(5) * t56 - t121 * t169 - t223 * t126;
t84 = t219 * t121;
t145 = -qJD(5) * t84 - t121 * t174 + t126 * t42;
t142 = pkin(8) * t130 - t204 - t211;
t141 = t104 * t51 + t130 * t52 + t199;
t12 = pkin(8) * t104 - t224 - t82;
t140 = -t105 * t36 - t12 + t221;
t138 = -pkin(8) * qJDD(4) + (t42 - t210) * qJD(4);
t108 = cos(t118);
t35 = t114 * t59;
t137 = t10 * t58 + t108 * t157 + t28 * t35;
t29 = t99 * t185 + (-t123 * t186 + (t127 * t128 - t192) * qJD(2)) * pkin(1);
t136 = -qJDD(4) * t52 + (t105 * t51 - t29) * qJD(4);
t135 = t157 - t160;
t14 = qJD(5) * t175 + t104 * t193 - t114 * t176 + t120 * t189 + t125 * t173;
t134 = -qJDD(4) * t97 + (t105 * t98 - t223) * qJD(4);
t4 = -t37 * t183 + qJDD(4) * pkin(4) - t12 * t121 + (-t173 - t194) * pkin(9);
t133 = t28 * t44 + t26 * t182 + g(1) * t106 + (-t26 * t114 - t4) * t120 + t221 * t108;
t5 = -pkin(9) * t144 + t12 * t126 - t37 * t184;
t132 = -g(1) * t108 + g(2) * t203 - g(3) * t202 + qJD(5) * t154 - t120 * t5 + t125 * t4 - t28 * t46;
t131 = t162 + t224;
t129 = cos(qJ(1));
t124 = sin(qJ(1));
t112 = qJDD(4) + qJDD(5);
t103 = t105 ^ 2;
t78 = t100 - t212;
t77 = qJDD(4) * t126 - t121 * t130;
t76 = qJDD(4) * t121 + t126 * t130;
t49 = t104 * t116 + 0.2e1 * t121 * t173;
t48 = t51 - t209;
t33 = -0.2e1 * t188 * t105 * qJD(4) + 0.2e1 * t121 * t189;
t31 = t36 * t184;
t27 = t102 + t30;
t22 = -t112 * t58 - t114 * t35;
t21 = t112 * t59 - t114 * t34;
t20 = -t121 * t29 + t126 * t170;
t19 = t121 * t170 + t126 * t29;
t18 = -t44 ^ 2 + t46 ^ 2;
t15 = t105 * t35 + t156;
t8 = t114 * t44 + t14;
t2 = t14 * t59 - t34 * t46;
t1 = -t14 * t58 - t15 * t59 + t34 * t44 - t35 * t46;
t3 = [qJDD(1), -g(2) * t129 - g(3) * t124, g(2) * t124 - g(3) * t129, t113, (t113 * t128 + t123 * t163) * pkin(1) + t149, ((-qJDD(1) - t113) * t123 + t128 * t163) * pkin(1) + t171, t104, t104 * t166 + t135 - t199, -t206 * t104 - t29 * t105 + t131, t49, t33, t76, t77, 0, t31 + t136 * t121 + (-t141 + t150) * t126, t121 * t141 + t126 * t136 + t165, t2, t1, t21, t22, 0, t27 * t44 + t48 * t15 + (-qJD(5) * t152 - t120 * t19 + t125 * t20) * t114 + t153 * t112 + t137, t27 * t46 + t48 * t14 - (qJD(5) * t153 + t120 * t20 + t125 * t19) * t114 - t152 * t112 + t167; 0, 0, 0, t113, pkin(1) * t123 * t164 + t149, (t128 * t164 - t181) * pkin(1) + t171, t104, t105 * t53 + (t104 * t127 - t105 * t186) * pkin(2) + t135, t105 * t54 + (-pkin(2) * t104 - t50) * t122 + ((-pkin(2) * t105 - t69) * qJD(3) - t139) * t127 + t162, t49, t33, t76, t77, 0, t31 + t134 * t121 + (t150 - t220) * t126, t220 * t121 + t134 * t126 + t165, t2, t1, t21, t22, 0, t78 * t15 + (-t120 * t57 + t125 * t56) * t112 + t207 * t44 + (t120 * t146 - t125 * t148) * t114 + t137, t78 * t14 - (t120 * t56 + t125 * t57) * t112 + t207 * t46 + (t120 * t148 + t125 * t146) * t114 + t167; 0, 0, 0, 0, 0, 0, t104, t135 + t204, t105 * t42 + t131, t49, t33, t76, t77, 0, t31 + t138 * t121 + (-t142 + t150) * t126, t121 * t142 + t126 * t138 + t165, t2, t1, t21, t22, 0, t100 * t15 + (-t120 * t85 + t125 * t84) * t112 + t158 * t44 + (t120 * t145 - t125 * t147) * t114 + t137, t100 * t14 - (t120 * t84 + t125 * t85) * t112 + t158 * t46 + (t120 * t147 + t125 * t145) * t114 + t167; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121 * t103 * t126, t188 * t103, t194, t189, qJDD(4), -g(1) * t126 + t121 * t140, g(1) * t121 + t126 * t140, t208, t18, t8, -t156, t112, -(t120 * t25 - t200) * t114 + (t125 * t112 - t114 * t182 - t44 * t197) * pkin(4) + t132, (-qJD(5) * t24 - t25 * t114 - t5) * t125 + (-qJD(5) * t125 * t114 - t120 * t112 - t46 * t197) * pkin(4) + t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t208, t18, t8, -t156, t112, -t114 * t154 + t132, (-t5 + (-qJD(5) + t114) * t24) * t125 + t133;];
tau_reg = t3;
