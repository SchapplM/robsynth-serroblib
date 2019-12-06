% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:33
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRPR6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:32:47
% EndTime: 2019-12-05 16:32:56
% DurationCPUTime: 1.92s
% Computational Cost: add. (1391->257), mult. (3788->402), div. (0->0), fcn. (2866->10), ass. (0->147)
t112 = cos(qJ(3));
t154 = t112 * qJD(2);
t98 = -qJD(5) + t154;
t192 = qJD(5) + t98;
t108 = sin(qJ(5));
t111 = cos(qJ(5));
t104 = sin(pkin(10));
t106 = cos(pkin(10));
t155 = t106 * qJD(3);
t109 = sin(qJ(3));
t161 = qJD(2) * t109;
t78 = t104 * t161 - t155;
t156 = t104 * qJD(3);
t80 = t106 * t161 + t156;
t126 = t108 * t78 - t111 * t80;
t191 = t126 * t98;
t110 = sin(qJ(2));
t105 = sin(pkin(5));
t164 = qJD(1) * t105;
t113 = cos(qJ(2));
t168 = t112 * t113;
t131 = pkin(3) * t109 - qJ(4) * t112;
t64 = qJD(3) * t131 - t109 * qJD(4);
t160 = qJD(3) * t109;
t152 = pkin(7) * t160;
t95 = t104 * t152;
t184 = t106 * t64 - (-t104 * t168 + t106 * t110) * t164 + t95;
t190 = t104 * t64 - (t104 * t110 + t106 * t168) * t164;
t107 = cos(pkin(5));
t163 = qJD(1) * t107;
t145 = t110 * t164;
t88 = qJD(2) * pkin(7) + t145;
t189 = -t109 * t88 + t112 * t163;
t188 = qJD(5) * t126;
t187 = pkin(8) + qJ(4);
t144 = t113 * t164;
t136 = qJD(2) * t144;
t26 = t112 * t136 + (qJD(4) + t189) * qJD(3);
t36 = (t64 + t145) * qJD(2);
t7 = t104 * t36 + t106 * t26;
t171 = t106 * t112;
t125 = pkin(4) * t109 - pkin(8) * t171;
t119 = t125 * qJD(3);
t186 = t119 + t184;
t172 = t106 * t109;
t176 = t104 * t112;
t185 = (-pkin(7) * t172 - pkin(8) * t176) * qJD(3) + t190;
t96 = t109 * t163;
t53 = t112 * t88 + t96;
t49 = qJD(3) * qJ(4) + t53;
t90 = -t112 * pkin(3) - t109 * qJ(4) - pkin(2);
t54 = qJD(2) * t90 - t144;
t14 = t104 * t54 + t106 * t49;
t139 = t106 * t152;
t183 = -t139 + t190;
t84 = t131 * qJD(2);
t22 = t104 * t84 + t106 * t189;
t169 = t111 * t106;
t170 = t108 * t104;
t82 = -t169 + t170;
t120 = t82 * t112;
t182 = qJD(2) * t120 - t82 * qJD(5);
t83 = t111 * t104 + t108 * t106;
t121 = t112 * t83;
t181 = -qJD(2) * t121 + t83 * qJD(5);
t153 = qJD(2) * qJD(3);
t142 = t112 * t153;
t157 = qJD(5) * t111;
t180 = t142 * t169 - t78 * t157;
t56 = pkin(7) * t171 + t104 * t90;
t179 = qJD(2) * pkin(2);
t159 = qJD(3) * t112;
t28 = qJD(3) * t96 + t109 * t136 + t88 * t159;
t178 = t28 * t104;
t177 = t28 * t106;
t175 = t105 * t110;
t174 = t105 * t113;
t115 = qJD(2) ^ 2;
t173 = t105 * t115;
t114 = qJD(3) ^ 2;
t167 = t114 * t109;
t166 = t114 * t112;
t165 = t109 ^ 2 - t112 ^ 2;
t162 = qJD(2) * t105;
t158 = qJD(5) * t109;
t151 = t110 * t173;
t150 = pkin(4) * t104 + pkin(7);
t6 = -t104 * t26 + t106 * t36;
t4 = qJD(2) * t119 + t6;
t135 = t104 * t142;
t5 = -pkin(8) * t135 + t7;
t149 = -t108 * t5 + t111 * t4;
t148 = t104 * t154;
t147 = t110 * t162;
t146 = t113 * t162;
t141 = t109 * t153;
t13 = -t104 * t49 + t106 * t54;
t21 = -t104 * t189 + t106 * t84;
t140 = -qJD(3) * pkin(3) + qJD(4);
t138 = t109 * t146;
t137 = t112 * t146;
t89 = -t144 - t179;
t134 = -t89 - t144;
t133 = t108 * t4 + t111 * t5;
t8 = -pkin(4) * t154 - t80 * pkin(8) + t13;
t9 = -t78 * pkin(8) + t14;
t132 = t108 * t9 - t111 * t8;
t2 = t108 * t8 + t111 * t9;
t77 = t106 * t90;
t34 = -pkin(8) * t172 + t77 + (-pkin(7) * t104 - pkin(4)) * t112;
t44 = -t104 * t109 * pkin(8) + t56;
t130 = -t108 * t44 + t111 * t34;
t129 = t108 * t34 + t111 * t44;
t70 = t107 * t109 + t112 * t175;
t40 = -t70 * t104 - t106 * t174;
t41 = -t104 * t174 + t70 * t106;
t128 = -t108 * t41 + t111 * t40;
t127 = t108 * t40 + t111 * t41;
t69 = -t107 * t112 + t109 * t175;
t93 = t187 * t106;
t123 = qJD(2) * t125 + qJD(4) * t104 + qJD(5) * t93 + t21;
t92 = t187 * t104;
t122 = pkin(8) * t148 + qJD(4) * t106 - qJD(5) * t92 - t22;
t118 = qJD(3) * t121;
t45 = t140 - t189;
t117 = qJD(3) * (-t134 - t179);
t116 = -qJ(4) * t160 + (t140 - t45) * t112;
t10 = (-qJD(5) * t80 - t135) * t108 + t180;
t11 = qJD(2) * t118 - t188;
t100 = -t106 * pkin(4) - pkin(3);
t85 = t150 * t109;
t73 = t150 * t159;
t65 = t111 * t78;
t62 = t82 * t109;
t61 = t83 * t109;
t55 = -pkin(7) * t176 + t77;
t43 = qJD(3) * t70 + t138;
t42 = -qJD(3) * t69 + t137;
t35 = pkin(4) * t148 + t53;
t29 = t108 * t80 + t65;
t25 = t157 * t172 - t158 * t170 + t118;
t24 = -qJD(3) * t120 - t158 * t83;
t23 = t78 * pkin(4) + t45;
t20 = t104 * t147 + t42 * t106;
t19 = -t42 * t104 + t106 * t147;
t16 = pkin(4) * t135 + t28;
t1 = [0, 0, -t151, -t113 * t173, 0, 0, 0, 0, 0, -t112 * t151 + (-t43 - t138) * qJD(3), t109 * t151 + (-t42 - t137) * qJD(3), t43 * t78 + (-t112 * t19 + (t40 * t109 + t176 * t69) * qJD(3)) * qJD(2), t43 * t80 + (t112 * t20 + (-t41 * t109 + t171 * t69) * qJD(3)) * qJD(2), -t19 * t80 - t20 * t78 + (-t104 * t41 - t106 * t40) * t142, t13 * t19 + t14 * t20 + t28 * t69 + t6 * t40 + t7 * t41 + t45 * t43, 0, 0, 0, 0, 0, -(-qJD(5) * t127 - t108 * t20 + t111 * t19) * t98 + t128 * t141 + t43 * t29 + t69 * t11, -t43 * t126 + t69 * t10 + (qJD(5) * t128 + t108 * t19 + t111 * t20) * t98 - t127 * t141; 0, 0, 0, 0, 0.2e1 * t112 * t141, -0.2e1 * t165 * t153, t166, -t167, 0, -pkin(7) * t166 + t109 * t117, pkin(7) * t167 + t112 * t117, (-t78 * t144 + t178 + (qJD(2) * t55 + t13) * qJD(3)) * t109 + (-t6 + (pkin(7) * t78 + t104 * t45) * qJD(3) + (t95 - t184) * qJD(2)) * t112, (-t80 * t144 + t177 + (-qJD(2) * t56 - t14) * qJD(3)) * t109 + (t7 + (pkin(7) * t80 + t106 * t45) * qJD(3) + (t139 + t183) * qJD(2)) * t112, -t184 * t80 - t183 * t78 + (-t104 * t7 - t106 * t6) * t109 + (-t104 * t14 - t106 * t13 + (-t104 * t56 - t106 * t55) * qJD(2)) * t159, -t45 * t109 * t144 + t6 * t55 + t7 * t56 + t183 * t14 + t184 * t13 + (t109 * t28 + t159 * t45) * pkin(7), -t10 * t62 - t126 * t24, -t10 * t61 + t62 * t11 + t126 * t25 - t24 * t29, -t10 * t112 - t24 * t98 + (-qJD(2) * t62 - t126) * t160, t11 * t112 + t25 * t98 + (-qJD(2) * t61 - t29) * t160, (-t98 - t154) * t160, -t149 * t112 + t73 * t29 + t85 * t11 + t16 * t61 + t23 * t25 + (t185 * t108 - t186 * t111) * t98 + (t112 * t2 + t129 * t98) * qJD(5) + (-t29 * t144 + (qJD(2) * t130 - t132) * qJD(3)) * t109, -t73 * t126 + t85 * t10 - t16 * t62 + t23 * t24 + t133 * t112 + (t186 * t108 + t185 * t111) * t98 + (-t112 * t132 + t130 * t98) * qJD(5) + (t126 * t144 + (-qJD(2) * t129 - t2) * qJD(3)) * t109; 0, 0, 0, 0, -t109 * t115 * t112, t165 * t115, 0, 0, 0, t53 * qJD(3) - t161 * t89 - t28, t134 * t154, -t177 - t53 * t78 + (t104 * t116 - t109 * t13 + t112 * t21) * qJD(2), t178 - t53 * t80 + (t106 * t116 + t109 * t14 - t112 * t22) * qJD(2), t21 * t80 + t22 * t78 + (-qJD(4) * t78 + t13 * t154 + t7) * t106 + (qJD(4) * t80 + t14 * t154 - t6) * t104, -t28 * pkin(3) - t13 * t21 - t14 * t22 - t45 * t53 + (-t104 * t13 + t106 * t14) * qJD(4) + (-t6 * t104 + t7 * t106) * qJ(4), t10 * t83 - t126 * t182, -t10 * t82 - t83 * t11 + t126 * t181 - t182 * t29, -t182 * t98 + (qJD(3) * t83 + t126) * t161, t181 * t98 + (-qJD(3) * t82 + t29) * t161, t98 * t161, t100 * t11 + t16 * t82 - t35 * t29 + (t108 * t122 + t111 * t123) * t98 + t181 * t23 + ((-t108 * t93 - t111 * t92) * qJD(3) + t132) * t161, t100 * t10 + t16 * t83 + t35 * t126 + (-t108 * t123 + t111 * t122) * t98 + t182 * t23 + (-(-t108 * t92 + t111 * t93) * qJD(3) + t2) * t161; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t80 + t156) * t154, (t78 + t155) * t154, -t78 ^ 2 - t80 ^ 2, t13 * t80 + t14 * t78 + t28, 0, 0, 0, 0, 0, t11 + t191, t65 * t98 + (-t135 + (-qJD(5) + t98) * t80) * t108 + t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t126 * t29, t126 ^ 2 - t29 ^ 2, -t29 * t98 + t10, -t142 * t83 + t188 + t191, t141, t23 * t126 - t192 * t2 + t149, t132 * t192 + t23 * t29 - t133;];
tauc_reg = t1;
