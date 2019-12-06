% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PRRPR5
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% tau_reg [5x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:28
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRRPR5_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR5_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR5_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR5_invdynJ_fixb_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:27:58
% EndTime: 2019-12-05 16:28:08
% DurationCPUTime: 2.00s
% Computational Cost: add. (1783->284), mult. (4379->421), div. (0->0), fcn. (3520->14), ass. (0->159)
t106 = sin(pkin(10));
t114 = sin(qJ(3));
t174 = qJD(2) * t114;
t109 = cos(pkin(10));
t117 = cos(qJ(3));
t181 = t109 * t117;
t71 = qJD(2) * t181 - t106 * t174;
t65 = qJD(5) - t71;
t207 = t65 - qJD(5);
t196 = qJ(4) + pkin(7);
t111 = cos(pkin(5));
t176 = qJD(1) * t111;
t108 = sin(pkin(5));
t115 = sin(qJ(2));
t118 = cos(qJ(2));
t170 = qJD(1) * qJD(2);
t59 = qJDD(2) * pkin(7) + (qJDD(1) * t115 + t118 * t170) * t108;
t123 = qJ(4) * qJDD(2) + qJD(2) * qJD(4) + qJD(3) * t176 + t59;
t175 = qJD(1) * t115;
t161 = t108 * t175;
t145 = qJD(2) * t196 + t161;
t137 = t145 * qJD(3);
t168 = t111 * qJDD(1);
t89 = t117 * t168;
t11 = qJDD(3) * pkin(3) - t114 * t123 - t117 * t137 + t89;
t12 = (-t137 + t168) * t114 + t123 * t117;
t3 = -t106 * t12 + t109 * t11;
t1 = -qJDD(3) * pkin(4) - t3;
t50 = t114 * t176 + t117 * t145;
t191 = t106 * t50;
t192 = qJD(3) * pkin(3);
t49 = -t114 * t145 + t117 * t176;
t46 = t49 + t192;
t17 = t109 * t46 - t191;
t13 = -qJD(3) * pkin(4) - t17;
t107 = sin(pkin(9));
t110 = cos(pkin(9));
t180 = t111 * t115;
t68 = t107 * t118 + t110 * t180;
t70 = -t107 * t180 + t110 * t118;
t143 = g(1) * t70 + g(2) * t68;
t184 = t108 * t115;
t130 = -g(3) * t184 - t143;
t150 = qJD(3) * t196;
t128 = -t114 * qJD(4) - t117 * t150;
t182 = t108 * t118;
t158 = qJD(1) * t182;
t66 = t117 * qJD(4) - t114 * t150;
t78 = t106 * t114 - t181;
t194 = t106 * t128 + t109 * t66 + t158 * t78;
t4 = t106 * t11 + t109 * t12;
t2 = qJDD(3) * pkin(8) + t4;
t97 = pkin(3) * t117 + pkin(2);
t64 = -qJD(2) * t97 + qJD(4) - t158;
t79 = t106 * t117 + t109 * t114;
t73 = t79 * qJD(2);
t23 = -pkin(4) * t71 - pkin(8) * t73 + t64;
t166 = t117 * qJDD(2);
t167 = t114 * qJDD(2);
t139 = -t106 * t167 + t109 * t166;
t72 = t79 * qJD(3);
t33 = qJD(2) * t72 + qJDD(5) - t139;
t35 = pkin(4) * t78 - pkin(8) * t79 - t97;
t157 = t196 * t114;
t83 = t196 * t117;
t54 = -t106 * t157 + t109 * t83;
t75 = t78 * qJD(3);
t206 = -(qJD(5) * t23 + t2) * t78 + t1 * t79 - t13 * t75 + (-qJD(5) * t35 - t194) * t65 - t54 * t33 + t130;
t179 = t111 * t118;
t67 = t107 * t115 - t110 * t179;
t69 = t107 * t179 + t110 * t115;
t144 = g(1) * t69 + g(2) * t67;
t103 = qJ(3) + pkin(10);
t99 = cos(t103);
t205 = -(g(3) * t182 - t144) * t99 + t35 * t33;
t185 = t108 * t110;
t186 = t107 * t108;
t94 = pkin(3) * t106 + pkin(8);
t98 = sin(t103);
t204 = (pkin(3) * t174 + pkin(4) * t73 - pkin(8) * t71 + qJD(5) * t94) * t65 + g(1) * (t186 * t99 - t70 * t98) + g(2) * (-t185 * t99 - t68 * t98) + g(3) * (t111 * t99 - t184 * t98) + t1;
t113 = sin(qJ(5));
t116 = cos(qJ(5));
t169 = qJD(2) * qJD(3);
t154 = t117 * t169;
t155 = t114 * t169;
t38 = qJDD(2) * t79 - t106 * t155 + t109 * t154;
t57 = qJD(3) * t113 + t116 * t73;
t16 = qJD(5) * t57 - qJDD(3) * t116 + t113 * t38;
t119 = qJD(3) ^ 2;
t152 = qJDD(1) * t182;
t156 = t115 * t170;
t87 = t108 * t156;
t203 = 0.2e1 * qJDD(2) * pkin(2) - pkin(7) * t119 + t108 * (-g(3) * t118 + t156) + t144 + t152 - t87;
t76 = t111 * t117 - t114 * t184;
t202 = g(3) * t76;
t171 = t116 * qJD(3);
t55 = t113 * t73 - t171;
t200 = t55 * t65;
t199 = t55 * t73;
t198 = t57 * t65;
t197 = t57 * t73;
t195 = t106 * t66 - t109 * t128 - t158 * t79;
t44 = t109 * t50;
t18 = t106 * t46 + t44;
t193 = qJD(2) * pkin(2);
t190 = t113 * t33;
t146 = t116 * t65;
t172 = qJD(5) * t113;
t15 = qJD(5) * t171 + qJDD(3) * t113 + t116 * t38 - t172 * t73;
t189 = t15 * t113;
t188 = qJD(5) * t79;
t183 = t108 * t117;
t178 = qJDD(1) - g(3);
t104 = t114 ^ 2;
t177 = -t117 ^ 2 + t104;
t173 = qJD(2) * t115;
t165 = t114 * t192;
t164 = t79 * t172;
t163 = t113 * t182;
t162 = t116 * t182;
t160 = t108 * t173;
t159 = qJD(2) * t182;
t153 = t118 * t169;
t147 = t108 * t178;
t142 = g(1) * t107 - g(2) * t110;
t141 = pkin(4) * t72 + pkin(8) * t75 - t161 + t165;
t140 = t33 * t79 - t65 * t75;
t14 = qJD(3) * pkin(8) + t18;
t6 = t113 * t23 + t116 * t14;
t138 = t113 * t14 - t116 * t23;
t136 = t116 * t33 + (t113 * t71 - t172) * t65;
t120 = qJD(2) ^ 2;
t135 = qJDD(2) * t118 - t115 * t120;
t77 = t111 * t114 + t115 * t183;
t29 = t106 * t76 + t109 * t77;
t134 = -t113 * t29 - t162;
t133 = -t116 * t29 + t163;
t82 = -t158 - t193;
t127 = -qJD(2) * t82 + t143 - t59;
t126 = pkin(3) * t155 - qJDD(2) * t97 + qJDD(4) + t87;
t22 = t109 * t49 - t191;
t125 = -t94 * t33 + (t13 + t22) * t65;
t122 = -pkin(7) * qJDD(3) + (t158 + t82 - t193) * qJD(3);
t36 = t126 - t152;
t95 = -pkin(3) * t109 - pkin(4);
t62 = t111 * t98 + t184 * t99;
t53 = t106 * t83 + t109 * t157;
t48 = -qJD(3) * t77 - t114 * t159;
t47 = qJD(3) * t76 + t117 * t159;
t42 = t186 * t98 + t70 * t99;
t40 = -t185 * t98 + t68 * t99;
t37 = -qJD(3) * t73 + t139;
t28 = t106 * t77 - t109 * t76;
t21 = t106 * t48 + t109 * t47;
t20 = t106 * t49 + t44;
t19 = t106 * t47 - t109 * t48;
t8 = -pkin(4) * t37 - pkin(8) * t38 + t36;
t7 = t116 * t8;
t5 = [t178, 0, t135 * t108, (-qJDD(2) * t115 - t118 * t120) * t108, 0, 0, 0, 0, 0, t48 * qJD(3) + t76 * qJDD(3) + (-t114 * t153 + t117 * t135) * t108, -t47 * qJD(3) - t77 * qJDD(3) + (-t114 * t135 - t117 * t153) * t108, t19 * t73 + t21 * t71 + t28 * t38 + t29 * t37, -t17 * t19 + t18 * t21 - t28 * t3 + t29 * t4 - g(3) + (-t118 * t36 + t173 * t64) * t108, 0, 0, 0, 0, 0, (qJD(5) * t133 - t113 * t21 + t116 * t160) * t65 + t134 * t33 + t19 * t55 + t28 * t16, -(qJD(5) * t134 + t113 * t160 + t116 * t21) * t65 + t133 * t33 + t19 * t57 + t28 * t15; 0, qJDD(2), t178 * t182 + t144, -t115 * t147 + t143, qJDD(2) * t104 + 0.2e1 * t114 * t154, 0.2e1 * t114 * t166 - 0.2e1 * t169 * t177, qJDD(3) * t114 + t117 * t119, qJDD(3) * t117 - t114 * t119, 0, t122 * t114 + t117 * t203, -t114 * t203 + t122 * t117, t17 * t75 - t18 * t72 + t194 * t71 + t195 * t73 - t3 * t79 + t37 * t54 + t38 * t53 - t4 * t78 + t130, t4 * t54 - t3 * t53 - t36 * t97 + t64 * t165 - g(1) * (t196 * t70 - t69 * t97) - g(2) * (t196 * t68 - t67 * t97) + t194 * t18 - t195 * t17 + (-t64 * t175 - g(3) * (t115 * t196 + t118 * t97)) * t108, -t57 * t164 + (t15 * t79 - t57 * t75) * t116, -(-t113 * t57 - t116 * t55) * t75 + (-t189 - t116 * t16 + (t113 * t55 - t116 * t57) * qJD(5)) * t79, t116 * t140 + t15 * t78 - t164 * t65 + t57 * t72, -t113 * t140 - t146 * t188 - t16 * t78 - t55 * t72, t33 * t78 + t65 * t72, t53 * t16 - t138 * t72 + t7 * t78 + t195 * t55 + (t141 * t65 + (t13 * t79 - t14 * t78 - t54 * t65) * qJD(5) + t205) * t116 + t206 * t113, t53 * t15 - t6 * t72 + t195 * t57 + (-(-qJD(5) * t14 + t8) * t78 - t13 * t188 + (qJD(5) * t54 - t141) * t65 - t205) * t113 + t206 * t116; 0, 0, 0, 0, -t114 * t120 * t117, t177 * t120, t167, t166, qJDD(3), t114 * t127 - t142 * t183 - t202 + t89, g(3) * t77 + (t108 * t142 - t168) * t114 + t127 * t117, (t18 - t20) * t73 + (t17 - t22) * t71 + (t106 * t37 - t109 * t38) * pkin(3), t17 * t20 - t18 * t22 + (t4 * t106 + t3 * t109 - t64 * t174 - g(1) * (t107 * t183 - t114 * t70) - g(2) * (-t110 * t183 - t114 * t68) - t202) * pkin(3), t146 * t57 + t189, (t15 - t200) * t116 + (-t16 - t198) * t113, t146 * t65 + t190 - t197, t136 + t199, -t65 * t73, t125 * t113 - t116 * t204 + t138 * t73 + t95 * t16 - t20 * t55, t113 * t204 + t125 * t116 + t95 * t15 - t20 * t57 + t6 * t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71 ^ 2 - t73 ^ 2, -t118 * t147 + t17 * t73 - t18 * t71 + t126 - t144, 0, 0, 0, 0, 0, t136 - t199, -t116 * t65 ^ 2 - t190 - t197; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57 * t55, -t55 ^ 2 + t57 ^ 2, t15 + t200, -t16 + t198, t33, -t113 * t2 + t7 - t13 * t57 - g(1) * (-t113 * t42 + t116 * t69) - g(2) * (-t113 * t40 + t116 * t67) - g(3) * (-t113 * t62 - t162) + t207 * t6, -t116 * t2 - t113 * t8 + t13 * t55 - g(1) * (-t113 * t69 - t116 * t42) - g(2) * (-t113 * t67 - t116 * t40) - g(3) * (-t116 * t62 + t163) - t207 * t138;];
tau_reg = t5;
