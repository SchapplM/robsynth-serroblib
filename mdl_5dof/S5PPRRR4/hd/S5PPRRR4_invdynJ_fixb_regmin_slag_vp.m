% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5PPRRR4
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PPRRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5PPRRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:19:46
% EndTime: 2019-12-05 15:19:53
% DurationCPUTime: 1.84s
% Computational Cost: add. (1333->248), mult. (3537->397), div. (0->0), fcn. (3395->14), ass. (0->149)
t91 = cos(qJ(4));
t153 = t91 * qJD(3);
t70 = -qJD(5) + t153;
t194 = t70 + qJD(5);
t168 = sin(pkin(10));
t80 = sin(pkin(11));
t133 = t168 * t80;
t84 = cos(pkin(10));
t86 = cos(pkin(5));
t174 = t84 * t86;
t83 = cos(pkin(11));
t106 = -t83 * t174 + t133;
t82 = sin(pkin(5));
t176 = t82 * t84;
t81 = sin(pkin(6));
t85 = cos(pkin(6));
t191 = t106 * t85 + t81 * t176;
t132 = t168 * t83;
t47 = t80 * t174 + t132;
t89 = sin(qJ(3));
t92 = cos(qJ(3));
t21 = t191 * t92 + t47 * t89;
t107 = t86 * t132 + t84 * t80;
t134 = t82 * t168;
t190 = t107 * t85 - t81 * t134;
t48 = -t86 * t133 + t84 * t83;
t23 = t190 * t92 + t48 * t89;
t177 = t81 * t92;
t146 = t86 * t177;
t175 = t83 * t85;
t148 = t92 * t175;
t179 = t80 * t89;
t149 = t82 * t179;
t36 = -t82 * t148 - t146 + t149;
t112 = g(1) * t23 + g(2) * t21 + g(3) * t36;
t117 = t89 * t175 + t80 * t92;
t108 = t117 * t82;
t178 = t81 * t89;
t69 = t86 * qJD(1) + qJD(2);
t31 = qJD(1) * t108 + t69 * t178;
t167 = qJD(1) * t82;
t144 = t83 * t167;
t128 = t85 * t144;
t152 = qJDD(1) * t82;
t139 = t83 * t152;
t165 = qJD(3) * t89;
t143 = t81 * t165;
t145 = t80 * t167;
t164 = qJD(3) * t92;
t67 = t86 * qJDD(1) + qJDD(2);
t96 = -(t85 * t139 + t67 * t81) * t92 + qJDD(1) * t149 + t128 * t165 + t69 * t143 + t145 * t164;
t193 = t31 * qJD(3) + t112 - t96;
t88 = sin(qJ(4));
t158 = qJD(5) * t88;
t192 = -qJD(3) * t158 + qJDD(4);
t151 = t88 * qJDD(3);
t87 = sin(qJ(5));
t90 = cos(qJ(5));
t33 = ((qJD(5) + t153) * qJD(4) + t151) * t87 - t192 * t90;
t30 = (t69 * t81 + t128) * t92 - t89 * t145;
t22 = -t191 * t89 + t47 * t92;
t24 = -t190 * t89 + t48 * t92;
t116 = -t82 * t83 * t81 + t86 * t85;
t37 = t86 * t178 + t108;
t25 = -t116 * t91 + t37 * t88;
t38 = t106 * t81 - t85 * t176;
t39 = t107 * t81 + t85 * t134;
t113 = g(1) * (-t24 * t88 + t39 * t91) + g(2) * (-t22 * t88 + t38 * t91) - g(3) * t25;
t125 = pkin(4) * t88 - pkin(9) * t91;
t118 = t148 - t179;
t15 = qJDD(3) * pkin(8) + (t69 * t164 + t67 * t89) * t81 + (qJD(1) * qJD(3) * t118 + qJDD(1) * t117) * t82;
t44 = -t81 * t139 + t85 * t67;
t172 = t91 * t44;
t29 = qJD(3) * pkin(8) + t31;
t45 = -t81 * t144 + t85 * t69;
t20 = t91 * t29 + t88 * t45;
t2 = -qJDD(4) * pkin(4) + t20 * qJD(4) + t88 * t15 - t172;
t189 = (pkin(9) * qJD(5) + t125 * qJD(3)) * t70 - t113 - t2;
t150 = qJD(3) * qJD(4);
t135 = t88 * t150;
t75 = t91 * qJDD(3);
t53 = qJDD(5) - t75 + t135;
t61 = t125 * qJD(4);
t64 = -t91 * pkin(4) - t88 * pkin(9) - pkin(3);
t188 = t64 * t53 + (t31 - t61) * t70;
t18 = qJD(4) * pkin(9) + t20;
t186 = (pkin(8) * t70 + t18) * qJD(5) + t112;
t137 = t91 * t150;
t154 = t90 * qJD(4);
t32 = qJD(5) * t154 + (t137 + t151) * t90 + t192 * t87;
t185 = t32 * t87;
t166 = qJD(3) * t88;
t54 = t87 * t166 - t154;
t184 = t54 * t70;
t155 = t87 * qJD(4);
t56 = t90 * t166 + t155;
t183 = t56 * t70;
t182 = t56 * t90;
t180 = t70 * t91;
t173 = t88 * t44;
t78 = t88 ^ 2;
t170 = -t91 ^ 2 + t78;
t169 = qJD(3) * pkin(3);
t163 = qJD(4) * t54;
t162 = qJD(4) * t56;
t161 = qJD(4) * t88;
t160 = qJD(4) * t91;
t159 = qJD(5) * t87;
t157 = qJD(5) * t90;
t142 = t81 * t164;
t141 = t70 * t155;
t140 = t70 * t154;
t28 = -t30 - t169;
t129 = -qJD(3) * t28 - t15;
t27 = qJD(3) * t64 - t30;
t4 = t90 * t18 + t87 * t27;
t124 = t87 * t18 - t90 * t27;
t26 = t116 * t88 + t37 * t91;
t14 = t26 * t90 + t36 * t87;
t13 = -t26 * t87 + t36 * t90;
t123 = t88 * t29 - t91 * t45;
t94 = qJD(3) ^ 2;
t122 = qJDD(3) * t92 - t89 * t94;
t50 = t91 * t178 + t88 * t85;
t120 = -t90 * t177 - t87 * t50;
t119 = t87 * t177 - t90 * t50;
t49 = t88 * t178 - t91 * t85;
t115 = -t70 * t157 + t87 * t53;
t114 = t70 * t159 + t90 * t53;
t111 = g(1) * t24 + g(2) * t22 + g(3) * t37;
t102 = -pkin(8) * qJDD(4) + (t28 + t30 - t169) * qJD(4);
t99 = qJD(5) * t64 * t70 - t111;
t17 = -qJD(4) * pkin(4) + t123;
t98 = -pkin(9) * t53 + (-t17 + t123) * t70;
t1 = qJDD(4) * pkin(9) - qJD(4) * t123 + t91 * t15 + t173;
t97 = -pkin(8) * t53 + qJD(4) * t17 + qJD(5) * t27 - t30 * t70 + t1;
t93 = qJD(4) ^ 2;
t95 = 0.2e1 * qJDD(3) * pkin(3) - pkin(8) * t93 + t193;
t41 = qJD(4) * t50 + t88 * t142;
t40 = -qJD(4) * t49 + t91 * t142;
t35 = t37 * qJD(3);
t34 = (t118 * t82 + t146) * qJD(3);
t12 = t24 * t91 + t39 * t88;
t10 = t22 * t91 + t38 * t88;
t8 = -qJD(4) * t25 + t34 * t91;
t7 = qJD(4) * t26 + t34 * t88;
t6 = qJD(3) * t61 + qJDD(3) * t64 + t96;
t5 = t90 * t6;
t3 = [qJDD(1) - g(3), t67 * t86 - g(3) + (t80 ^ 2 + t83 ^ 2) * t82 ^ 2 * qJDD(1), 0, -t35 * qJD(3) - t36 * qJDD(3), -t34 * qJD(3) - t37 * qJDD(3), 0, 0, 0, 0, 0, -t36 * t75 - t7 * qJD(4) - t25 * qJDD(4) + (t161 * t36 - t35 * t91) * qJD(3), t36 * t151 - t8 * qJD(4) - t26 * qJDD(4) + (t160 * t36 + t35 * t88) * qJD(3), 0, 0, 0, 0, 0, -(-qJD(5) * t14 + t35 * t90 - t8 * t87) * t70 + t13 * t53 + t7 * t54 + t25 * t33, (qJD(5) * t13 + t35 * t87 + t8 * t90) * t70 - t14 * t53 + t7 * t56 + t25 * t32; 0, -g(3) * t86 + (-t168 * g(1) + g(2) * t84) * t82 + t67, 0, t122 * t81, (-qJDD(3) * t89 - t92 * t94) * t81, 0, 0, 0, 0, 0, -t41 * qJD(4) - t49 * qJDD(4) + (t122 * t91 - t135 * t92) * t81, -t40 * qJD(4) - t50 * qJDD(4) + (-t122 * t88 - t137 * t92) * t81, 0, 0, 0, 0, 0, -(qJD(5) * t119 + t143 * t90 - t87 * t40) * t70 + t120 * t53 + t41 * t54 + t49 * t33, (qJD(5) * t120 + t143 * t87 + t90 * t40) * t70 + t119 * t53 + t41 * t56 + t49 * t32; 0, 0, qJDD(3), t193, -t67 * t178 - t117 * t152 + (-t118 * t167 - t69 * t177 + t30) * qJD(3) + t111, t78 * qJDD(3) + 0.2e1 * t91 * t135, -0.2e1 * t170 * t150 + 0.2e1 * t88 * t75, qJDD(4) * t88 + t93 * t91, qJDD(4) * t91 - t93 * t88, 0, t102 * t88 + t91 * t95, t102 * t91 - t88 * t95, t32 * t90 * t88 + (t154 * t91 - t158 * t87) * t56, (-t54 * t90 - t56 * t87) * t160 + (-t185 - t33 * t90 + (t54 * t87 - t182) * qJD(5)) * t88, (-t32 - t140) * t91 + (t114 + t162) * t88, (t33 + t141) * t91 + (-t115 - t163) * t88, -t161 * t70 - t53 * t91, t188 * t90 + t99 * t87 + (pkin(8) * t163 + t186 * t90 + t97 * t87 - t5) * t91 + (t17 * t157 - t124 * qJD(4) + t2 * t87 - t30 * t54 + (t33 - t141) * pkin(8)) * t88, -t188 * t87 + t99 * t90 + (pkin(8) * t162 + t97 * t90 + (-t186 + t6) * t87) * t91 + (-t17 * t159 - t4 * qJD(4) + t2 * t90 - t30 * t56 + (t32 - t140) * pkin(8)) * t88; 0, 0, 0, 0, 0, -t88 * t94 * t91, t170 * t94, t151, t75, qJDD(4), t129 * t88 - t113 + t172, g(1) * t12 + g(2) * t10 + g(3) * t26 + t129 * t91 - t173, -t182 * t70 + t185, (t32 + t184) * t90 + (-t33 + t183) * t87, (t180 * t90 - t88 * t56) * qJD(3) + t115, (-t180 * t87 + t88 * t54) * qJD(3) + t114, t70 * t166, -pkin(4) * t33 + t124 * t166 + t189 * t90 - t20 * t54 + t98 * t87, -pkin(4) * t32 + t4 * t166 - t189 * t87 - t20 * t56 + t98 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56 * t54, -t54 ^ 2 + t56 ^ 2, t32 - t184, -t183 - t33, t53, -t87 * t1 + t5 - t17 * t56 - g(1) * (-t12 * t87 + t23 * t90) - g(2) * (-t10 * t87 + t21 * t90) - g(3) * t13 - t194 * t4, -t90 * t1 - t87 * t6 + t17 * t54 - g(1) * (-t12 * t90 - t23 * t87) - g(2) * (-t10 * t90 - t21 * t87) + g(3) * t14 + t194 * t124;];
tau_reg = t3;
