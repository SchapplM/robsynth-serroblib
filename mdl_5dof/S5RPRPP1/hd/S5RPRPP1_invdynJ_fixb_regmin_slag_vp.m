% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% tau_reg [5x19]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:14
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPRPP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:14:04
% EndTime: 2021-01-15 11:14:08
% DurationCPUTime: 0.96s
% Computational Cost: add. (1538->227), mult. (3149->281), div. (0->0), fcn. (2020->12), ass. (0->123)
t162 = 2 * qJD(3);
t93 = sin(pkin(7));
t75 = t93 * pkin(1) + pkin(6);
t145 = qJ(4) + t75;
t88 = qJ(1) + pkin(7);
t80 = sin(t88);
t82 = cos(t88);
t126 = g(1) * t82 + g(2) * t80;
t92 = sin(pkin(8));
t94 = cos(pkin(8));
t97 = sin(qJ(3));
t99 = cos(qJ(3));
t55 = t92 * t99 + t94 * t97;
t50 = t55 * qJD(1);
t45 = t50 ^ 2;
t149 = t94 * t99;
t133 = qJD(1) * t149;
t144 = qJD(1) * t97;
t47 = t92 * t144 - t133;
t161 = -t47 ^ 2 - t45;
t95 = cos(pkin(7));
t77 = -t95 * pkin(1) - pkin(2);
t85 = t99 * pkin(3);
t160 = t77 - t85;
t155 = g(1) * t80;
t159 = g(2) * t82 - t155;
t130 = t145 * t97;
t53 = t145 * t99;
t29 = -t92 * t130 + t94 * t53;
t87 = qJ(3) + pkin(8);
t79 = sin(t87);
t158 = -t29 * qJDD(3) + t159 * t79;
t63 = t75 * qJDD(1);
t108 = qJ(4) * qJDD(1) + qJD(1) * qJD(4) + (qJD(3) * qJD(2)) + t63;
t128 = t145 * qJD(1);
t119 = t128 * qJD(3);
t83 = t99 * qJDD(2);
t14 = qJDD(3) * pkin(3) - t108 * t97 - t99 * t119 + t83;
t17 = (qJDD(2) - t119) * t97 + t108 * t99;
t4 = t92 * t14 + t94 * t17;
t81 = cos(t87);
t157 = g(3) * t79 + t126 * t81 - t4;
t156 = pkin(3) * t97;
t152 = g(3) * t99;
t3 = t94 * t14 - t92 * t17;
t151 = t81 * t82;
t41 = t97 * qJD(2) + t128 * t99;
t150 = t92 * t41;
t35 = t94 * t41;
t40 = t99 * qJD(2) - t128 * t97;
t37 = qJD(3) * pkin(3) + t40;
t16 = t92 * t37 + t35;
t100 = cos(qJ(1));
t78 = t85 + pkin(2);
t148 = t100 * pkin(1) + t82 * t78;
t90 = t97 ^ 2;
t147 = -t99 ^ 2 + t90;
t146 = t79 * qJ(5);
t66 = qJD(1) * t77;
t143 = qJDD(3) * pkin(4);
t19 = t92 * t40 + t35;
t142 = t19 * qJD(3);
t141 = t97 * qJD(3);
t20 = t94 * t40 - t150;
t140 = qJD(5) - t20;
t139 = qJDD(2) - g(3);
t137 = t97 * qJDD(1);
t136 = t99 * qJDD(1);
t135 = qJD(1) * qJD(3);
t134 = pkin(3) * t141;
t132 = t97 * t135;
t131 = t99 * t135;
t129 = qJD(3) * t145;
t96 = -qJ(4) - pkin(6);
t98 = sin(qJ(1));
t124 = -t98 * pkin(1) - t82 * t96;
t123 = g(1) * t98 - g(2) * t100;
t122 = -t94 * t136 + t92 * t137;
t121 = -t81 * pkin(4) - t146;
t15 = t94 * t37 - t150;
t49 = t55 * qJD(3);
t54 = t92 * t97 - t149;
t118 = -t49 * qJD(3) - t54 * qJDD(3);
t117 = -g(3) * t81 + t126 * t79 + t3;
t30 = qJD(1) * t49 + t122;
t111 = t55 * qJDD(1) - t92 * t132;
t31 = t94 * t131 + t111;
t52 = qJD(3) * t149 - t92 * t141;
t116 = -t55 * t30 + t54 * t31 - t52 * t47 + t49 * t50;
t28 = t94 * t130 + t92 * t53;
t115 = -g(2) * t151 - t28 * qJDD(3) + t81 * t155;
t113 = -t97 * qJD(4) - t99 * t129;
t112 = -t66 * qJD(1) + t126 - t63;
t46 = qJD(1) * t160 + qJD(4);
t110 = -qJDD(3) * t75 + t66 * t162;
t24 = t47 * pkin(4) - t50 * qJ(5) + t46;
t109 = -t24 * t50 - qJDD(5) + t117;
t39 = pkin(3) * t132 + qJDD(1) * t160 + qJDD(4);
t101 = qJD(3) ^ 2;
t107 = -0.2e1 * qJDD(1) * t77 - t101 * t75 - t159;
t42 = t99 * qJD(4) - t97 * t129;
t22 = -t94 * t113 + t92 * t42;
t23 = t92 * t113 + t94 * t42;
t106 = t22 * t50 - t23 * t47 + t28 * t31 - t29 * t30 - t126;
t105 = t30 * pkin(4) - t31 * qJ(5) + t39;
t104 = t50 * t162 + t122;
t102 = qJD(1) ^ 2;
t89 = qJDD(3) * qJ(5);
t76 = -t94 * pkin(3) - pkin(4);
t72 = t92 * pkin(3) + qJ(5);
t62 = qJDD(3) * t99 - t101 * t97;
t61 = qJDD(3) * t97 + t101 * t99;
t32 = t52 * qJD(3) + t55 * qJDD(3);
t27 = t54 * pkin(4) - t55 * qJ(5) + t160;
t26 = pkin(3) * t144 + t50 * pkin(4) + t47 * qJ(5);
t21 = (-t47 + t133) * qJD(3) + t111;
t18 = t49 * pkin(4) - t52 * qJ(5) - t55 * qJD(5) + t134;
t11 = qJD(3) * qJ(5) + t16;
t10 = -qJD(3) * pkin(4) + qJD(5) - t15;
t5 = -t50 * qJD(5) + t105;
t2 = qJDD(5) - t143 - t3;
t1 = qJD(3) * qJD(5) + t4 + t89;
t6 = [qJDD(1), t123, g(1) * t100 + g(2) * t98, (t123 + (t93 ^ 2 + t95 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), t90 * qJDD(1) + 0.2e1 * t97 * t131, -0.2e1 * t147 * t135 + 0.2e1 * t97 * t136, t61, t62, 0, t107 * t99 + t110 * t97, -t107 * t97 + t110 * t99, t160 * t30 + t39 * t54 + t46 * t49 + (t47 * t156 - t22) * qJD(3) + t115, t160 * t31 + t39 * t55 + t46 * t52 + (t50 * t156 - t23) * qJD(3) + t158, -t15 * t52 - t16 * t49 - t3 * t55 - t4 * t54 + t106, t4 * t29 + t16 * t23 - t3 * t28 - t15 * t22 + t39 * t160 + t46 * t134 - g(1) * (-t80 * t78 + t124) - g(2) * (-t80 * t96 + t148), -t22 * qJD(3) + t18 * t47 + t24 * t49 + t27 * t30 + t5 * t54 + t115, -t1 * t54 + t10 * t52 - t11 * t49 + t2 * t55 + t106, t23 * qJD(3) - t18 * t50 - t24 * t52 - t27 * t31 - t5 * t55 - t158, t1 * t29 + t11 * t23 + t5 * t27 + t24 * t18 + t2 * t28 + t10 * t22 - g(1) * t124 - g(2) * (pkin(4) * t151 + t82 * t146 + t148) + (-g(1) * (t121 - t78) + g(2) * t96) * t80; 0, 0, 0, t139, 0, 0, 0, 0, 0, t62, -t61, t118, -t32, t116, -t15 * t49 + t16 * t52 - t3 * t54 + t4 * t55 - g(3), t118, t116, t32, t1 * t55 + t10 * t49 + t11 * t52 + t2 * t54 - g(3); 0, 0, 0, 0, -t97 * t102 * t99, t147 * t102, t137, t136, qJDD(3), t112 * t97 - t152 + t83, t112 * t99 - t139 * t97, t142 - t46 * t50 + (qJDD(3) * t94 - t47 * t144) * pkin(3) + t117, t20 * qJD(3) + t46 * t47 + (-qJDD(3) * t92 - t50 * t144) * pkin(3) + t157, (t16 - t19) * t50 + (-t15 + t20) * t47 + (-t30 * t92 - t31 * t94) * pkin(3), t15 * t19 - t16 * t20 + (-t152 + t3 * t94 + t4 * t92 + (-qJD(1) * t46 + t126) * t97) * pkin(3), t142 - t26 * t47 + (pkin(4) - t76) * qJDD(3) + t109, -t72 * t30 + t76 * t31 + (t11 - t19) * t50 + (t10 - t140) * t47, t72 * qJDD(3) - t24 * t47 + t26 * t50 + t89 + (0.2e1 * qJD(5) - t20) * qJD(3) - t157, t1 * t72 + t2 * t76 - t24 * t26 - t10 * t19 - g(3) * (-t121 + t85) + t140 * t11 + t126 * (pkin(4) * t79 - qJ(5) * t81 + t156); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, t21, t161, t15 * t50 + t16 * t47 + t159 + t39, t104, t161, -t21, t11 * t47 + (-qJD(5) - t10) * t50 + t105 + t159; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50 * t47 - qJDD(3), (t47 + t133) * qJD(3) + t111, -t45 - t101, -t11 * qJD(3) - t109 - t143;];
tau_reg = t6;
