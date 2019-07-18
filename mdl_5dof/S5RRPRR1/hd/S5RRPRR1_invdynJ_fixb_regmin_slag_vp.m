% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RRPRR1
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d4,d5]';
% 
% Output:
% tau_reg [5x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 17:22
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RRPRR1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(4,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S5RRPRR1_invdynJ_fixb_regmin_slag_vp: pkin has to be [4x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 17:22:29
% EndTime: 2019-07-18 17:22:35
% DurationCPUTime: 1.91s
% Computational Cost: add. (1720->251), mult. (3914->352), div. (0->0), fcn. (2736->10), ass. (0->153)
t104 = qJD(2) + qJD(4);
t111 = sin(qJ(4));
t115 = cos(qJ(2));
t198 = cos(qJ(4));
t154 = qJD(1) * t198;
t112 = sin(qJ(2));
t171 = qJD(1) * t112;
t206 = -t111 * t171 + t115 * t154;
t180 = t104 * t206;
t107 = qJ(2) + qJ(4);
t101 = sin(t107);
t113 = sin(qJ(1));
t116 = cos(qJ(1));
t138 = g(1) * t116 + g(2) * t113;
t207 = t138 * t101;
t146 = qJDD(1) * t198;
t162 = t112 * qJDD(1);
t63 = t111 * t115 + t198 * t112;
t38 = t104 * t63;
t22 = t38 * qJD(1) + t111 * t162 - t115 * t146;
t117 = pkin(1) + pkin(2);
t57 = qJD(5) - t206;
t145 = t57 ^ 2;
t205 = pkin(4) * t145;
t102 = cos(t107);
t195 = g(3) * t102;
t108 = qJDD(2) * pkin(1);
t183 = pkin(3) + qJ(3);
t147 = qJD(2) * t183;
t55 = -t112 * qJD(3) - t115 * t147;
t71 = t183 * t112;
t29 = qJDD(2) * pkin(2) + t55 * qJD(1) - qJDD(1) * t71 + t108;
t54 = t115 * qJD(3) - t112 * t147;
t72 = t183 * t115;
t33 = t54 * qJD(1) + qJDD(1) * t72;
t148 = t111 * t33 - t198 * t29;
t65 = qJD(1) * t72;
t159 = t198 * t65;
t109 = qJD(2) * pkin(1);
t64 = qJD(1) * t71;
t49 = qJD(2) * pkin(2) + t109 - t64;
t31 = t111 * t49 + t159;
t7 = t31 * qJD(4) + t148;
t204 = t7 + t195;
t67 = t117 * t171;
t86 = t111 * t117 + pkin(4);
t202 = (-pkin(4) * t206 + qJD(5) * t86 + t67) * t57;
t73 = t117 * t115;
t61 = -qJD(1) * t73 + qJD(3);
t95 = g(3) * t101;
t201 = t102 * t138 - t206 * t61 + t95;
t103 = qJDD(2) + qJDD(4);
t110 = sin(qJ(5));
t114 = cos(qJ(5));
t161 = t115 * qJDD(1);
t21 = t111 * t161 + t112 * t146 + t180;
t170 = qJD(1) * t115;
t60 = -t111 * t170 - t112 * t154;
t44 = t110 * t104 - t114 * t60;
t9 = t44 * qJD(5) - t114 * t103 + t110 * t21;
t127 = -t111 * t112 + t198 * t115;
t168 = qJD(4) * t111;
t139 = -t65 * t168 + t198 * t33;
t153 = t198 * qJD(4);
t122 = t111 * t29 + t49 * t153 + t139;
t36 = t60 * pkin(4) + t61;
t150 = t103 * pkin(4) + qJD(5) * t36 + t122;
t129 = -t111 * t72 - t198 * t71;
t16 = t129 * qJD(4) + t111 * t55 + t198 * t54;
t20 = qJDD(5) + t22;
t160 = t198 * t49;
t181 = t111 * t65;
t30 = -t160 + t181;
t37 = t104 * t127;
t45 = -t111 * t71 + t198 * t72;
t47 = -t63 * pkin(4) - t73;
t200 = -(qJD(5) * t47 + t16) * t57 + t150 * t127 - t45 * t20 + t30 * t37 + t7 * t63;
t106 = t115 ^ 2;
t199 = 0.2e1 * t106;
t194 = g(3) * t115;
t193 = t30 * t63;
t41 = -t114 * t104 - t110 * t60;
t192 = t41 * t57;
t191 = t44 * t57;
t190 = t44 * t60;
t189 = t47 * t20;
t188 = t57 * t60;
t187 = t60 * t41;
t186 = t60 * t206;
t166 = qJD(5) * t114;
t167 = qJD(5) * t110;
t8 = t110 * t103 + t104 * t166 + t114 * t21 + t60 * t167;
t184 = t8 * t110;
t66 = t117 * qJD(2) * t112;
t182 = t110 * t20;
t179 = t60 * t104;
t119 = qJD(1) ^ 2;
t178 = t112 * t119;
t177 = t113 * t110;
t176 = t113 * t114;
t175 = t116 * t110;
t174 = t116 * t114;
t105 = t112 ^ 2;
t173 = t105 - t106;
t172 = t105 + t199;
t164 = qJD(1) * qJD(2);
t152 = t112 * t164;
t165 = pkin(1) * t152 + qJDD(3);
t163 = qJD(1) * qJD(3);
t158 = t63 * t167;
t156 = t117 * t198;
t151 = t115 * t164;
t40 = pkin(2) * t152 - qJDD(1) * t73 + t165;
t13 = -t21 * pkin(4) + t40;
t25 = t104 * pkin(4) + t31;
t144 = qJD(5) * t25 - t13;
t141 = t114 * t57;
t140 = t112 * t151;
t137 = g(1) * t113 - g(2) * t116;
t136 = t20 * t63 + t37 * t57;
t12 = t110 * t36 + t114 * t25;
t134 = t204 * t110 - t12 * t60 + t30 * t166;
t11 = -t110 * t25 + t114 * t36;
t133 = t11 * t60 + t114 * t207 + t30 * t167;
t132 = t114 * t20 + (t110 * t206 - t167) * t57;
t131 = -t150 + t95;
t130 = -pkin(4) * t20 + (-t57 - t206) * t30;
t70 = -qJ(3) * t171 + t109;
t125 = -t70 * qJD(2) * t115 - t138;
t124 = t61 * t60 - t148 - t195 + t207;
t123 = pkin(1) * t161 + t137 - t165;
t35 = -t198 * t64 - t181;
t121 = -t86 * t20 - t30 * t206 + (-t117 * t153 + t35) * t57;
t120 = qJ(3) ^ 2;
t118 = qJD(2) ^ 2;
t82 = -pkin(1) * t170 + qJD(3);
t53 = t102 * t174 + t177;
t52 = -t102 * t175 + t176;
t51 = -t102 * t176 + t175;
t50 = t102 * t177 + t174;
t46 = -t112 * t163 + t108 + (-t151 - t162) * qJ(3);
t34 = -t111 * t64 + t159;
t24 = -t37 * pkin(4) + t66;
t23 = -t206 ^ 2 + t60 ^ 2;
t17 = t45 * qJD(4) + t111 * t54 - t198 * t55;
t15 = -t179 - t22;
t14 = t21 - t180;
t10 = t114 * t13;
t4 = t57 * t141 + t182 + t190;
t3 = t132 - t187;
t2 = t44 * t141 + t184;
t1 = (t8 - t192) * t114 + (-t9 - t191) * t110;
t5 = [qJDD(1), t137, t138, t105 * qJDD(1) + 0.2e1 * t140, 0.2e1 * t112 * t161 - 0.2e1 * t173 * t164, qJDD(2) * t112 + t118 * t115, qJDD(2) * t115 - t118 * t112, 0, t137 * t115, -t137 * t112, -t46 * t112 + t172 * t163 + (t172 * qJDD(1) - 0.2e1 * t140) * qJ(3) + t125, t106 * qJDD(1) * t120 + t123 * t115 * pkin(1) + (t163 * t199 + t125) * qJ(3) + (-t46 * qJ(3) - t70 * qJD(3) + (pkin(1) * t82 - 0.2e1 * t120 * t170) * qJD(2)) * t112, t21 * t63 - t60 * t37, t127 * t21 + t206 * t37 - t63 * t22 + t60 * t38, t63 * t103 + t37 * t104, t103 * t127 - t38 * t104, 0, t137 * t102 + t103 * t129 - t17 * t104 - t127 * t40 - t206 * t66 - t73 * t22 + t61 * t38, -t137 * t101 - t45 * t103 - t16 * t104 - t73 * t21 + t61 * t37 + t40 * t63 - t66 * t60, -t44 * t158 + (t37 * t44 + t63 * t8) * t114, (-t110 * t44 - t114 * t41) * t37 + (-t184 - t114 * t9 + (t110 * t41 - t114 * t44) * qJD(5)) * t63, t136 * t114 - t127 * t8 - t57 * t158 + t44 * t38, -t63 * t57 * t166 - t136 * t110 + t127 * t9 - t41 * t38, -t127 * t20 + t57 * t38, -g(1) * t51 - g(2) * t53 - t10 * t127 + t11 * t38 + t17 * t41 - t129 * t9 + (t189 + t24 * t57 + (t127 * t25 - t45 * t57 + t193) * qJD(5)) * t114 + t200 * t110, -g(1) * t50 - g(2) * t52 - t12 * t38 + t17 * t44 - t129 * t8 + (-(-qJD(5) * t45 + t24) * t57 - t189 - t144 * t127 - qJD(5) * t193) * t110 + t200 * t114; 0, 0, 0, -t115 * t178, t173 * t119, t162, t161, qJDD(2), t138 * t112 - t194, g(3) * t112 + t138 * t115, -pkin(1) * t162 + (qJ(3) * t178 + (t70 - t109) * qJD(1)) * t115, (qJ(3) * qJD(1) * t70 + t120 * t178) * t115 + (-t194 + t46 + (-qJD(1) * t82 + t138) * t112) * pkin(1), t186, t23, t14, t15, t103, t103 * t156 + t34 * t104 + t67 * t206 + (-t159 + (-t104 * t117 - t49) * t111) * qJD(4) + t124, t35 * t104 + t67 * t60 + (-t103 * t117 - t29) * t111 + (-t104 * t156 - t160) * qJD(4) - t139 + t201, t2, t1, t4, t3, t188, -t34 * t41 + (t41 * t168 - t198 * t9) * t117 + (-t204 - t202) * t114 + t121 * t110 + t133, -t34 * t44 + (t44 * t168 - t198 * t8) * t117 + t121 * t114 + (-t207 + t202) * t110 + t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t105 - t106) * t119, -t106 * t119 * qJ(3) + t70 * t171 - t123, 0, 0, 0, 0, 0, t22 - t179, t21 + t180, 0, 0, 0, 0, 0, t132 + t187, -t114 * t145 - t182 + t190; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t186, t23, t14, t15, t103, t124 + (-qJD(4) + t104) * t31, -t30 * t104 - t122 + t201, t2, t1, t4, t3, t188, -t31 * t41 + t130 * t110 + (-t204 - t205) * t114 + t133, -t31 * t44 + t130 * t114 + (-t207 + t205) * t110 + t134; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44 * t41, -t41 ^ 2 + t44 ^ 2, t8 + t192, t191 - t9, t20, -g(1) * t52 + g(2) * t50 + t110 * t131 + t12 * t57 - t25 * t166 - t30 * t44 + t10, g(1) * t53 - g(2) * t51 + t11 * t57 + t110 * t144 + t114 * t131 + t30 * t41;];
tau_reg  = t5;
