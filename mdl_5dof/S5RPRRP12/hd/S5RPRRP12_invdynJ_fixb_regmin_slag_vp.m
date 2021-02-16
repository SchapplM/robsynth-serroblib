% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP12
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tau_reg [5x24]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP12_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_invdynJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:25:53
% EndTime: 2021-01-15 19:26:03
% DurationCPUTime: 2.03s
% Computational Cost: add. (1923->308), mult. (3692->398), div. (0->0), fcn. (2194->6), ass. (0->165)
t89 = cos(qJ(4));
t149 = t89 * qJD(3);
t87 = sin(qJ(3));
t132 = t87 * t149;
t90 = cos(qJ(3));
t153 = qJD(4) * t90;
t86 = sin(qJ(4));
t107 = t86 * t153 + t132;
t145 = t90 * qJDD(1);
t15 = t107 * qJD(1) - qJD(4) * t149 - t86 * qJDD(3) - t89 * t145;
t161 = qJD(1) * t90;
t52 = t86 * t161 - t149;
t150 = t87 * qJD(1);
t70 = qJD(4) + t150;
t186 = t52 * t70;
t207 = -t15 - t186;
t88 = sin(qJ(1));
t194 = g(1) * t88;
t91 = cos(qJ(1));
t82 = g(2) * t91;
t200 = -t82 + t194;
t206 = -qJDD(2) + t200;
t77 = t89 * pkin(4) + pkin(3);
t85 = -qJ(5) - pkin(7);
t205 = -t77 * t87 - t85 * t90;
t181 = t87 * t89;
t59 = t87 * pkin(3) - t90 * pkin(7) + qJ(2);
t92 = -pkin(1) - pkin(6);
t168 = t92 * t181 + t86 * t59;
t35 = t59 * qJD(1);
t68 = t92 * qJD(1) + qJD(2);
t58 = t87 * t68;
t42 = qJD(3) * pkin(7) + t58;
t14 = t86 * t35 + t89 * t42;
t10 = -t52 * qJ(5) + t14;
t13 = t89 * t35 - t86 * t42;
t151 = t86 * qJD(3);
t54 = t89 * t161 + t151;
t9 = -t54 * qJ(5) + t13;
t8 = t70 * pkin(4) + t9;
t116 = t10 * t86 + t8 * t89;
t204 = -t116 * qJD(1) - t200;
t152 = qJDD(1) * pkin(1);
t203 = t152 + t206;
t202 = -qJ(2) + t205;
t133 = t87 * t151;
t199 = t54 * qJD(4);
t16 = -qJD(1) * t133 - t89 * qJDD(3) + t86 * t145 + t199;
t185 = t54 * t70;
t201 = t16 + t185;
t74 = t86 * pkin(4) - t92;
t191 = g(3) * t90;
t173 = t91 * t89;
t180 = t88 * t86;
t36 = -t87 * t180 + t173;
t174 = t91 * t86;
t179 = t88 * t89;
t38 = t87 * t174 + t179;
t198 = -g(1) * t36 - g(2) * t38 + t86 * t191;
t197 = t54 ^ 2;
t196 = -t9 + t8;
t192 = g(3) * t87;
t144 = qJD(1) * qJD(3);
t130 = t90 * t144;
t146 = t87 * qJDD(1);
t49 = qJDD(4) + t130 + t146;
t190 = t49 * pkin(4);
t189 = t52 * pkin(4);
t187 = t15 * t86;
t184 = t54 * t89;
t182 = t86 * t49;
t178 = t89 * t49;
t177 = t90 * t15;
t66 = t92 * qJDD(1) + qJDD(2);
t176 = t90 * t66;
t175 = t90 * t68;
t127 = qJD(4) * t85;
t138 = t86 * t150;
t148 = t89 * qJD(5);
t121 = pkin(3) * t90 + pkin(7) * t87;
t57 = t121 * qJD(1);
t170 = t89 * t175 + t86 * t57;
t172 = -qJ(5) * t138 + t86 * t127 + t148 - t170;
t41 = t89 * t57;
t171 = t89 * t127 - t41 - (pkin(4) * t90 + qJ(5) * t181) * qJD(1) + (-qJD(5) + t175) * t86;
t141 = t90 * t82;
t169 = g(3) * t181 + t89 * t141;
t84 = t90 ^ 2;
t167 = t87 ^ 2 - t84;
t93 = qJD(3) ^ 2;
t94 = qJD(1) ^ 2;
t166 = -t93 - t94;
t165 = qJ(5) * t90;
t164 = t15 * qJ(5);
t163 = t16 * qJ(5);
t162 = t94 * qJ(2);
t160 = qJD(3) * t52;
t159 = qJD(3) * t54;
t158 = qJD(3) * t87;
t157 = qJD(3) * t90;
t156 = qJD(3) * t92;
t155 = qJD(4) * t86;
t154 = qJD(4) * t89;
t147 = qJDD(3) * t87;
t143 = qJDD(1) * qJ(2);
t142 = t90 * t194;
t137 = t90 * t156;
t50 = t121 * qJD(3) + qJD(2);
t140 = t89 * t137 + t59 * t154 + t86 * t50;
t136 = t70 * t155;
t135 = t89 * t153;
t129 = -qJD(5) - t189;
t43 = -qJD(3) * pkin(3) - t175;
t22 = -t129 + t43;
t134 = t22 * t154;
t131 = -t86 * t92 + pkin(4);
t20 = t50 * qJD(1) + t59 * qJDD(1);
t25 = qJDD(3) * pkin(7) + t68 * t157 + t87 * t66;
t126 = -t35 * t154 + t42 * t155 - t86 * t20 - t89 * t25;
t125 = -qJDD(3) * pkin(3) + t68 * t158;
t123 = qJD(4) * t87 + qJD(1);
t24 = t125 - t176;
t122 = -qJD(4) * pkin(7) * t70 - t24;
t120 = g(1) * t38 - g(2) * t36;
t37 = t87 * t179 + t174;
t39 = t87 * t173 - t180;
t119 = -g(1) * t39 - g(2) * t37;
t118 = g(1) * t91 + g(2) * t88;
t117 = t10 * t89 - t8 * t86;
t114 = -t141 - t192;
t113 = -t66 + t200;
t111 = t70 * t154 + t182;
t110 = -t136 + t178;
t108 = 0.2e1 * qJ(2) * t144 + qJDD(3) * t92;
t106 = t16 * pkin(4) + qJDD(5) + t125;
t105 = 0.2e1 * qJD(1) * qJD(2) - t118;
t104 = t113 + t162;
t103 = -pkin(7) * t49 + t70 * t43;
t102 = g(1) * t37 - g(2) * t39 + t89 * t191 + t126;
t18 = t89 * t20;
t100 = -t14 * qJD(4) - t86 * t25 + t18;
t99 = t105 + 0.2e1 * t143;
t1 = -t54 * qJD(5) + t100 + t164 + t190;
t2 = -t52 * qJD(5) - t126 - t163;
t98 = -t116 * qJD(4) - t1 * t86 + t2 * t89;
t96 = -t92 * t93 + t99;
t95 = t100 + t198;
t80 = qJDD(3) * t90;
t64 = t86 * t142;
t62 = t85 * t89;
t61 = t85 * t86;
t51 = t74 * t90;
t48 = t52 ^ 2;
t47 = t89 * t59;
t34 = t89 * t50;
t27 = -pkin(4) * t138 + t58;
t26 = t87 * t156 + (-t133 + t135) * pkin(4);
t23 = -t86 * t165 + t168;
t19 = t131 * t87 - t89 * t165 + t47;
t7 = t106 - t176;
t6 = -qJ(5) * t135 + (-qJD(5) * t90 + (qJ(5) * qJD(3) - qJD(4) * t92) * t87) * t86 + t140;
t5 = qJ(5) * t132 + t34 - t168 * qJD(4) + (qJ(5) * t155 + t131 * qJD(3) - t148) * t90;
t4 = -t90 * t16 + (t160 - t182) * t87 + (-t123 * t89 - t90 * t151) * t70;
t3 = t177 + (t159 - t178) * t87 + (t123 * t86 - t90 * t149) * t70;
t11 = [qJDD(1), t200, t118, -0.2e1 * t152 - t206, t99, t203 * pkin(1) + (t105 + t143) * qJ(2), t84 * qJDD(1) - 0.2e1 * t87 * t130, 0.2e1 * t167 * t144 - 0.2e1 * t87 * t145, -t93 * t87 + t80, -t93 * t90 - t147, 0, t108 * t90 + t96 * t87, -t108 * t87 + t96 * t90, -t107 * t54 - t89 * t177, (t52 * t89 + t54 * t86) * t158 + (t187 - t16 * t89 + (t52 * t86 - t184) * qJD(4)) * t90, (-t70 * t149 - t15) * t87 + (t110 + t159) * t90, (t151 * t70 - t16) * t87 + (-t111 - t160) * t90, t157 * t70 + t49 * t87, t34 * t70 + t47 * t49 + (t52 * t156 + t18 + (-t70 * t92 - t42) * t154) * t87 + (t13 * qJD(3) + t154 * t43 - t92 * t16) * t90 + ((-qJD(4) * t59 - t137) * t70 + t24 * t90 + (-t43 * qJD(3) - qJD(4) * t35 - t92 * t49 - t25) * t87) * t86 + t119, -t140 * t70 - t168 * t49 + (t92 * t136 + (-t43 * t89 + t54 * t92) * qJD(3) + t126) * t87 + (-t14 * qJD(3) + t92 * t15 - t155 * t43 + t24 * t89) * t90 + t120, t51 * t16 + t19 * t49 + t26 * t52 + t5 * t70 + (-t151 * t22 + t1) * t87 + (qJD(3) * t8 + t7 * t86 + t134) * t90 + t119, -t51 * t15 - t23 * t49 + t26 * t54 - t6 * t70 + (-t149 * t22 - t2) * t87 + (-qJD(3) * t10 - t155 * t22 + t7 * t89) * t90 + t120, t19 * t15 - t23 * t16 - t5 * t54 - t6 * t52 + t116 * t158 + (-qJD(4) * t117 - t1 * t89 - t2 * t86 + t118) * t90, t2 * t23 + t10 * t6 + t1 * t19 + t8 * t5 + t7 * t51 + t22 * t26 - g(1) * (-t202 * t91 - t74 * t88) - g(2) * (-t202 * t88 + t74 * t91); 0, 0, 0, qJDD(1), -t94, -t162 - t203, 0, 0, 0, 0, 0, t166 * t87 + t80, t166 * t90 - t147, 0, 0, 0, 0, 0, t4, t3, t4, t3, (t123 * t54 - t157 * t52 - t16 * t87) * t89 + (t123 * t52 - t15 * t87 + t157 * t54) * t86, (t117 * qJD(3) - t7) * t90 + (qJD(3) * t22 + t98) * t87 + t204; 0, 0, 0, 0, 0, 0, t90 * t94 * t87, -t167 * t94, t145, -t146, qJDD(3), -t104 * t90 + t192, t104 * t87 + t191, t70 * t184 - t187, -t201 * t86 + t207 * t89, (t70 * t181 - t54 * t90) * qJD(1) + t111, (-t70 * t86 * t87 + t52 * t90) * qJD(1) + t110, -t70 * t161, -t13 * t161 - t52 * t58 - pkin(3) * t16 - t41 * t70 + (t122 - t142) * t89 + (t70 * t175 + t103) * t86 + t169, pkin(3) * t15 + t170 * t70 + t14 * t161 - t54 * t58 + t64 + t103 * t89 + (t114 - t122) * t86, -t8 * t161 - t77 * t16 - t27 * t52 + t61 * t49 + (-t7 - t142) * t89 + t171 * t70 + (t22 * t150 + (t22 + t189) * qJD(4)) * t86 + t169, t134 + t77 * t15 - t27 * t54 + t62 * t49 + t64 - t172 * t70 + (t10 * t90 + t22 * t181) * qJD(1) + (pkin(4) * t199 + t114 + t7) * t86, t61 * t15 + t62 * t16 - t171 * t54 - t172 * t52 + t204 * t87 - t191 + t98, -t2 * t62 + t1 * t61 - t7 * t77 - g(3) * t205 + t171 * t8 + (pkin(4) * t155 - t27) * t22 + t172 * t10 - t200 * (t90 * t77 - t85 * t87); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t52, -t48 + t197, -t15 + t186, -t16 + t185, t49, t14 * t70 - t43 * t54 + t95, t13 * t70 + t43 * t52 + t102, 0.2e1 * t190 + t164 + t10 * t70 + (t129 - t22) * t54 + t95, -t197 * pkin(4) + t163 + t9 * t70 + (qJD(5) + t22) * t52 + t102, t15 * pkin(4) - t196 * t52, t196 * t10 + (-t22 * t54 + t1 + t198) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t201, t207, -t48 - t197, t10 * t52 + t113 * t90 + t8 * t54 + t106 - t192;];
tau_reg = t11;
