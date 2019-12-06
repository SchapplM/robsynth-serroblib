% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRRP3
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
% tau_reg [5x20]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRRP3_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP3_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP3_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:03:57
% EndTime: 2019-12-05 18:04:02
% DurationCPUTime: 1.08s
% Computational Cost: add. (1366->186), mult. (2905->249), div. (0->0), fcn. (1922->12), ass. (0->124)
t91 = qJ(1) + pkin(8);
t81 = sin(t91);
t82 = cos(t91);
t167 = -g(2) * t81 + g(3) * t82;
t94 = qJ(3) + qJ(4);
t85 = sin(t94);
t86 = cos(t94);
t169 = -g(1) * t86 + t167 * t85;
t159 = cos(qJ(4));
t100 = cos(qJ(3));
t95 = sin(pkin(8));
t75 = t95 * pkin(1) + pkin(6);
t160 = pkin(7) + t75;
t129 = t160 * qJD(1);
t98 = sin(qJ(3));
t33 = t98 * qJD(2) + t129 * t100;
t97 = sin(qJ(4));
t28 = t97 * t33;
t148 = qJD(3) * pkin(3);
t32 = t100 * qJD(2) - t129 * t98;
t31 = t32 + t148;
t134 = t159 * t31 - t28;
t54 = t97 * t100 + t159 * t98;
t48 = t54 * qJD(1);
t146 = t48 * qJ(5);
t168 = t146 - t134;
t90 = qJD(3) + qJD(4);
t51 = t160 * t98;
t52 = t160 * t100;
t151 = t159 * t52 - t97 * t51;
t96 = cos(pkin(8));
t76 = -t96 * pkin(1) - pkin(2);
t87 = t100 * pkin(3);
t166 = t76 - t87;
t165 = t48 ^ 2;
t5 = t90 * pkin(4) - t168;
t164 = t5 + t168;
t136 = t159 * t100;
t125 = qJD(1) * t136;
t145 = qJD(1) * t98;
t137 = t97 * t145;
t46 = -t125 + t137;
t50 = t166 * qJD(1);
t23 = t46 * pkin(4) + qJD(5) + t50;
t158 = t23 * t48;
t157 = t48 * t46;
t156 = t81 * t86;
t155 = t82 * t86;
t154 = t97 * t98;
t128 = qJDD(1) * t159;
t141 = t98 * qJDD(1);
t121 = -t100 * t128 + t97 * t141;
t26 = t90 * t54;
t18 = t26 * qJD(1) + t121;
t120 = t90 * t154;
t130 = t159 * qJD(4);
t25 = -qJD(3) * t136 - t100 * t130 + t120;
t153 = -t54 * t18 + t25 * t46;
t152 = t159 * t32 - t28;
t150 = pkin(4) * t86 + t87;
t92 = t98 ^ 2;
t149 = -t100 ^ 2 + t92;
t147 = t46 * qJ(5);
t65 = qJD(1) * t76;
t144 = qJD(4) * t97;
t142 = qJDD(2) - g(1);
t140 = qJD(1) * qJD(3);
t139 = t100 * qJDD(1);
t138 = t98 * t148;
t30 = t159 * t33;
t135 = t98 * t140;
t133 = -t97 * t32 - t30;
t132 = -t159 * t51 - t97 * t52;
t131 = qJD(3) * t160;
t127 = -t90 * t125 - t98 * t128 - t97 * t139;
t62 = t75 * qJDD(1);
t126 = pkin(7) * qJDD(1) + t62;
t124 = -g(2) * t82 - g(3) * t81;
t101 = cos(qJ(1));
t99 = sin(qJ(1));
t122 = g(2) * t101 + g(3) * t99;
t17 = qJD(1) * t120 + t127;
t53 = -t136 + t154;
t119 = -t53 * t17 + t48 * t26;
t88 = qJDD(3) + qJDD(4);
t118 = t25 * t90 - t54 * t88;
t115 = -t97 * t31 - t30;
t43 = t98 * t131;
t44 = t100 * t131;
t114 = -t51 * t130 - t52 * t144 - t159 * t43 - t97 * t44;
t113 = -t65 * qJD(1) + t167 - t62;
t34 = pkin(3) * t135 + qJDD(1) * t166;
t112 = 0.2e1 * t65 * qJD(3) - qJDD(3) * t75;
t102 = qJD(3) ^ 2;
t111 = -0.2e1 * qJDD(1) * t76 - t102 * t75 - t124;
t83 = t100 * qJDD(2);
t14 = qJDD(3) * pkin(3) - qJD(3) * t33 - t126 * t98 + t83;
t15 = qJD(3) * t32 + t98 * qJDD(2) + t126 * t100;
t110 = t115 * qJD(4) + t159 * t14 - t97 * t15;
t109 = -qJD(4) * t151 - t159 * t44 + t97 * t43;
t108 = t31 * t130 + t97 * t14 - t33 * t144 + t159 * t15;
t107 = t18 * pkin(4) + qJDD(5) + t34;
t106 = g(1) * t85 - g(2) * t156 + g(3) * t155 + t50 * t46 - t108;
t105 = -t50 * t48 + t110 + t169;
t103 = qJD(1) ^ 2;
t89 = -qJ(5) - pkin(7) - pkin(6);
t80 = t159 * pkin(3) + pkin(4);
t61 = qJDD(3) * t100 - t102 * t98;
t60 = qJDD(3) * t98 + t102 * t100;
t55 = pkin(2) + t150;
t45 = t46 ^ 2;
t21 = -t45 + t165;
t20 = -t53 * qJ(5) + t151;
t19 = -t54 * qJ(5) + t132;
t16 = -t26 * t90 - t53 * t88;
t10 = -t127 + (-t137 + t46) * t90;
t9 = -t146 + t152;
t8 = t133 + t147;
t7 = -t115 - t147;
t4 = t25 * qJ(5) - t54 * qJD(5) + t109;
t3 = -t26 * qJ(5) - t53 * qJD(5) + t114;
t2 = -t18 * qJ(5) - t46 * qJD(5) + t108;
t1 = t88 * pkin(4) + t17 * qJ(5) - t48 * qJD(5) + t110;
t6 = [qJDD(1), t122, -g(2) * t99 + g(3) * t101, (t122 + (t95 ^ 2 + t96 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t92 * qJDD(1) + 0.2e1 * t100 * t135, 0.2e1 * t98 * t139 - 0.2e1 * t149 * t140, t60, t61, 0, t111 * t100 + t112 * t98, t112 * t100 - t111 * t98, -t17 * t54 - t48 * t25, -t119 + t153, -t118, t16, 0, g(2) * t155 + g(3) * t156 + t109 * t90 + t132 * t88 + t46 * t138 + t166 * t18 + t50 * t26 + t34 * t53, -t114 * t90 + t124 * t85 + t48 * t138 - t151 * t88 - t166 * t17 - t50 * t25 + t34 * t54, -t1 * t54 + t19 * t17 - t20 * t18 - t2 * t53 + t5 * t25 - t7 * t26 - t3 * t46 - t4 * t48 - t167, t2 * t20 + t7 * t3 + t1 * t19 + t5 * t4 + t107 * (t53 * pkin(4) + t166) + t23 * (t26 * pkin(4) + t138) - g(2) * (-t101 * pkin(1) - t82 * t55 + t81 * t89) - g(3) * (-t99 * pkin(1) - t81 * t55 - t82 * t89); 0, 0, 0, t142, 0, 0, 0, 0, 0, t61, -t60, 0, 0, 0, 0, 0, t16, t118, t119 + t153, -t1 * t53 + t2 * t54 - t7 * t25 - t5 * t26 - g(1); 0, 0, 0, 0, -t98 * t103 * t100, t149 * t103, t141, t139, qJDD(3), -g(1) * t100 + t113 * t98 + t83, t113 * t100 - t142 * t98, t157, t21, t10, -t121, t88, -t133 * t90 + (-t90 * t144 - t46 * t145 + t159 * t88) * pkin(3) + t105, t152 * t90 + (-t90 * t130 - t48 * t145 - t97 * t88) * pkin(3) + t106, t80 * t17 + (t7 + t8) * t48 + (-t5 + t9) * t46 + (-t18 * t97 + (-t159 * t46 + t48 * t97) * qJD(4)) * pkin(3), t1 * t80 - t7 * t9 - t5 * t8 - pkin(4) * t158 - g(1) * t150 - t167 * (-t98 * pkin(3) - pkin(4) * t85) + (-t23 * t145 + t2 * t97 + (t159 * t7 - t5 * t97) * qJD(4)) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t157, t21, t10, -t121, t88, -t115 * t90 + t105, t134 * t90 + t106, pkin(4) * t17 - t164 * t46, t164 * t7 + (t1 - t158 + t169) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45 - t165, t7 * t46 + t5 * t48 + t107 + t124;];
tau_reg = t6;
