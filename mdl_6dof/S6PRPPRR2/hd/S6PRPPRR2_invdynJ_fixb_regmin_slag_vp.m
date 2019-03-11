% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% tau_reg [6x22]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6PRPPRR2_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR2_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR2_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:20:03
% EndTime: 2019-03-08 19:20:08
% DurationCPUTime: 2.02s
% Computational Cost: add. (1340->281), mult. (3058->393), div. (0->0), fcn. (2624->12), ass. (0->155)
t104 = sin(qJ(2));
t107 = cos(qJ(2));
t96 = sin(pkin(11));
t99 = cos(pkin(11));
t208 = -t104 * t96 + t107 * t99;
t103 = sin(qJ(5));
t85 = qJD(2) * t103 + qJD(6);
t207 = t85 - qJD(6);
t98 = sin(pkin(6));
t174 = qJD(1) * t98;
t145 = qJD(2) * t174;
t156 = t98 * qJDD(1);
t206 = t104 * t156 + t107 * t145;
t205 = t208 * t98;
t100 = cos(pkin(10));
t106 = cos(qJ(5));
t184 = t103 * t98;
t101 = cos(pkin(6));
t121 = t208 * t101;
t181 = t104 * t99;
t129 = t107 * t96 + t181;
t97 = sin(pkin(10));
t27 = t100 * t121 - t129 * t97;
t30 = -t100 * t129 - t97 * t121;
t38 = t101 * t103 + t106 * t205;
t120 = g(1) * (-t106 * t30 - t97 * t184) + g(2) * (t100 * t184 - t106 * t27) - g(3) * t38;
t132 = pkin(5) * t106 + pkin(9) * t103;
t80 = t107 * t156;
t51 = qJDD(2) * pkin(2) - t104 * t145 + t80;
t22 = -t206 * t96 + t51 * t99;
t128 = qJDD(4) - t22;
t199 = -pkin(3) - pkin(8);
t11 = t199 * qJDD(2) + t128;
t81 = t101 * qJDD(1) + qJDD(3);
t185 = t103 * t81;
t147 = t107 * t174;
t67 = qJD(2) * pkin(2) + t147;
t148 = t104 * t174;
t71 = t96 * t148;
t41 = t67 * t99 - t71;
t134 = qJD(4) - t41;
t35 = t199 * qJD(2) + t134;
t83 = qJD(1) * t101 + qJD(3);
t21 = t103 * t35 + t106 * t83;
t2 = -qJDD(5) * pkin(5) + qJD(5) * t21 - t106 * t11 + t185;
t204 = (pkin(9) * qJD(6) + t132 * qJD(2)) * t85 + t120 + t2;
t168 = t101 * t107;
t190 = -t101 * t181 - t96 * t168;
t31 = t100 * t208 + t97 * t190;
t50 = t99 * t147 - t71;
t127 = pkin(5) * t103 - pkin(9) * t106 + qJ(4);
t198 = pkin(2) * t96;
t55 = t127 + t198;
t155 = qJD(2) * qJD(5);
t140 = t106 * t155;
t152 = t103 * qJDD(2);
t57 = qJDD(6) + t140 + t152;
t58 = t132 * qJD(5) + qJD(4);
t203 = -t55 * t57 + (t50 - t58) * t85;
t102 = sin(qJ(6));
t105 = cos(qJ(6));
t162 = qJD(5) * t103;
t144 = t102 * t162;
t151 = t106 * qJDD(2);
t163 = qJD(5) * t102;
t164 = qJD(2) * t106;
t63 = t105 * t164 + t163;
t34 = -qJD(2) * t144 + t63 * qJD(6) - t105 * qJDD(5) + t102 * t151;
t26 = t100 * t190 - t208 * t97;
t160 = qJD(6) * t102;
t143 = t106 * t160;
t202 = -t105 * (-t106 * t57 + t85 * t162) - t85 * t143;
t42 = t99 * t148 + t67 * t96;
t37 = qJD(2) * qJ(4) + t42;
t53 = t129 * t98;
t47 = qJD(1) * t53;
t88 = -pkin(2) * t99 - pkin(3);
t84 = -pkin(8) + t88;
t86 = qJ(4) + t198;
t201 = qJDD(5) * t84 + (qJD(2) * t86 + t37 - t47) * qJD(5);
t119 = -g(1) * t31 + g(2) * t26 - g(3) * t53;
t19 = qJD(5) * pkin(9) + t21;
t200 = -(t84 * t85 + t19) * qJD(6) + t119;
t197 = t81 * t101 - g(3);
t157 = t105 * qJD(5);
t61 = t102 * t164 - t157;
t195 = t61 * t85;
t194 = t63 * t85;
t161 = qJD(5) * t106;
t33 = t105 * t151 + qJD(6) * t157 + t102 * qJDD(5) + (-t103 * t157 - t143) * qJD(2);
t192 = t33 * t103 + t63 * t161;
t95 = t106 ^ 2;
t189 = t103 ^ 2 - t95;
t187 = t102 * t85;
t186 = t103 * t34;
t182 = t104 * t97;
t180 = t105 * t85;
t179 = t106 * t81;
t178 = t106 * t98;
t177 = t107 * t98;
t175 = t33 * t102;
t173 = qJD(2) * t37;
t172 = qJD(2) * t85;
t171 = qJD(5) * t61;
t170 = qJD(5) * t84;
t169 = qJDD(2) * pkin(3);
t167 = qJD(4) - t50;
t166 = qJDD(1) - g(3);
t108 = qJD(5) ^ 2;
t109 = qJD(2) ^ 2;
t165 = -t108 - t109;
t159 = qJD(6) * t105;
t154 = qJDD(5) * t103;
t153 = qJDD(5) * t106;
t23 = t206 * t99 + t96 * t51;
t146 = t100 * t168;
t138 = -t11 + t173;
t133 = -g(1) * t97 + g(2) * t100;
t25 = t127 * qJD(2) + t42;
t4 = t102 * t25 + t105 * t19;
t131 = t102 * t19 - t105 * t25;
t39 = t101 * t106 - t103 * t205;
t10 = t102 * t53 + t105 * t39;
t9 = -t102 * t39 + t105 * t53;
t130 = t103 * t83 - t106 * t35;
t124 = t102 * t57 + t85 * t159;
t123 = t105 * t57 - t85 * t160;
t122 = -t100 * t104 - t97 * t168;
t118 = g(1) * t30 + g(2) * t27 + g(3) * t205;
t116 = -t124 + t171;
t115 = -qJD(6) * t55 * t85 - t118;
t18 = -qJD(5) * pkin(5) + t130;
t114 = -pkin(9) * t57 + (t18 - t130) * t85;
t1 = qJDD(5) * pkin(9) - t130 * qJD(5) + t103 * t11 + t179;
t113 = -qJD(5) * t18 - qJD(6) * t25 + t47 * t85 - t57 * t84 - t1;
t112 = -g(1) * t122 - g(3) * t177;
t111 = t118 + t128;
t12 = qJDD(2) * qJ(4) + qJD(4) * qJD(2) + t23;
t110 = t167 * qJD(2) + qJDD(2) * t86 - t108 * t84 + t119 + t12;
t75 = pkin(2) * t146;
t74 = -t103 * t108 + t153;
t73 = -t106 * t108 - t154;
t56 = t85 * t144;
t49 = t205 * qJD(2);
t48 = qJD(2) * t53;
t40 = -g(3) * t101 + t133 * t98 + t81;
t36 = -qJD(2) * pkin(3) + t134;
t17 = t128 - t169;
t16 = t100 * t178 + t27 * t103;
t14 = -t103 * t30 + t97 * t178;
t8 = t39 * qJD(5) - t48 * t106;
t7 = -t38 * qJD(5) + t48 * t103;
t6 = t58 * qJD(2) + t127 * qJDD(2) + t23;
t5 = t105 * t6;
t3 = [t166, 0 (qJDD(2) * t107 - t104 * t109) * t98 (-qJDD(2) * t104 - t107 * t109) * t98, t205 * t22 + t23 * t53 - t41 * t48 + t42 * t49 + t197, qJD(2) * t48 - qJDD(2) * t205, qJD(2) * t49 + qJDD(2) * t53, t12 * t53 - t17 * t205 + t36 * t48 + t37 * t49 + t197, 0, 0, 0, 0, 0, t53 * t152 - qJD(5) * t8 - qJDD(5) * t38 + (t103 * t49 + t53 * t161) * qJD(2), t53 * t151 - qJD(5) * t7 - qJDD(5) * t39 + (t106 * t49 - t53 * t162) * qJD(2), 0, 0, 0, 0, 0 (-qJD(6) * t10 - t7 * t102 + t49 * t105) * t85 + t9 * t57 + t8 * t61 + t38 * t34 -(qJD(6) * t9 + t49 * t102 + t7 * t105) * t85 - t10 * t57 + t8 * t63 + t38 * t33; 0, qJDD(2), t80 - g(2) * (t146 - t182) + t112 (g(1) * t100 + g(2) * t97) * t107 + (t133 * t101 - t166 * t98) * t104, -g(2) * t75 + t41 * t47 - t42 * t50 + (g(2) * t182 + t22 * t99 + t23 * t96 + t112) * pkin(2), -qJD(2) * t47 + (-pkin(3) + t88) * qJDD(2) + t111 (qJ(4) + t86) * qJDD(2) + (0.2e1 * qJD(4) - t50) * qJD(2) + t119 + t23, t12 * t86 + t17 * t88 - t36 * t47 - g(1) * (t122 * pkin(2) + pkin(3) * t30 + qJ(4) * t31) - g(2) * (-pkin(2) * t182 + pkin(3) * t27 - qJ(4) * t26 + t75) - g(3) * (pkin(2) * t177 + pkin(3) * t205 + qJ(4) * t53) + t167 * t37, qJDD(2) * t95 - 0.2e1 * t103 * t140, -0.2e1 * t103 * t151 + 0.2e1 * t189 * t155, t74, t73, 0, t110 * t103 + t201 * t106, -t201 * t103 + t110 * t106, -t63 * t143 + (t106 * t33 - t63 * t162) * t105 (t102 * t63 + t105 * t61) * t162 + (-t175 - t105 * t34 + (t102 * t61 - t105 * t63) * qJD(6)) * t106, t192 + t202, -t186 + t56 + (-t124 - t171) * t106, t103 * t57 + t85 * t161, -t203 * t105 + t115 * t102 + (t113 * t102 + t200 * t105 + t61 * t170 + t5) * t103 + (t18 * t159 + t2 * t102 - t84 * t34 + t47 * t61 + (-t84 * t187 - t131) * qJD(5)) * t106, t203 * t102 + t115 * t105 + (t63 * t170 + t113 * t105 + (-t200 - t6) * t102) * t103 + (-t18 * t160 + t2 * t105 - t84 * t33 + t47 * t63 + (-t84 * t180 - t4) * qJD(5)) * t106; 0, 0, 0, 0, t40, 0, 0, t40, 0, 0, 0, 0, 0, t73, -t74, 0, 0, 0, 0, 0, t106 * t116 + t186 + t56, t192 - t202; 0, 0, 0, 0, 0, qJDD(2), -t109, t111 - t169 - t173, 0, 0, 0, 0, 0, t165 * t103 + t153, t165 * t106 - t154, 0, 0, 0, 0, 0, -t105 * t172 + (-t85 * t163 - t34) * t106 + t116 * t103, t102 * t172 + (-t85 * t157 - t33) * t106 + (qJD(5) * t63 - t123) * t103; 0, 0, 0, 0, 0, 0, 0, 0, t106 * t109 * t103, -t189 * t109, t151, -t152, qJDD(5), -t106 * t138 - t120 - t185, g(1) * t14 - g(2) * t16 + g(3) * t39 + t138 * t103 - t179, t63 * t180 + t175 (t33 - t195) * t105 + (-t34 - t194) * t102 (t103 * t180 - t106 * t63) * qJD(2) + t124 (-t103 * t187 + t106 * t61) * qJD(2) + t123, -t85 * t164, -pkin(5) * t34 + t114 * t102 - t204 * t105 + t131 * t164 - t21 * t61, -pkin(5) * t33 + t204 * t102 + t114 * t105 + t4 * t164 - t21 * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63 * t61, -t61 ^ 2 + t63 ^ 2, t33 + t195, t194 - t34, t57, -t102 * t1 + t5 - t18 * t63 - g(1) * (-t102 * t14 + t105 * t31) - g(2) * (t102 * t16 - t105 * t26) - g(3) * t9 + t207 * t4, -t105 * t1 - t102 * t6 + t18 * t61 - g(1) * (-t102 * t31 - t105 * t14) - g(2) * (t102 * t26 + t105 * t16) + g(3) * t10 - t207 * t131;];
tau_reg  = t3;
